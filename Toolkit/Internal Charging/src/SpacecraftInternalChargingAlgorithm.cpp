#include "SpacecraftInternalChargingAlgorithm.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{
namespace
{
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kMinMassDensityKgPerM3 = 1.0;
constexpr double kMinTimeScaleS = 1.0e-18;
}

bool SpacecraftInternalChargingAlgorithm::initialize(const InternalChargingConfiguration& config)
{
    config_ = config;
    material_ = material_database_.findByName(config.material_name);
    if (material_ == nullptr)
    {
        return false;
    }
    if (!transport_model_.initialize(config.thickness_m, config.layers))
    {
        return false;
    }

    status_ = InternalChargingStatus{};
    status_.volume_charge_density_c_per_m3.assign(config.layers, 0.0);
    status_.electric_field_v_per_m.assign(config.layers, 0.0);
    initialized_ = true;
    return true;
}

bool SpacecraftInternalChargingAlgorithm::advance(double dt)
{
    if (!initialized_ || material_ == nullptr || dt <= 0.0)
    {
        return false;
    }

    double incident_current_density_a_per_m2 = config_.incident_current_density_a_per_m2;
    double incident_energy_ev = config_.incident_energy_ev;
    double incident_charge_state_abs = config_.incident_charge_state_abs;
    if (config_.source_mode == InternalChargingSourceMode::Radiation)
    {
        if (!radiation_drive_.has_value())
        {
            return false;
        }
        incident_current_density_a_per_m2 = radiation_drive_->incident_current_density_a_per_m2;
        incident_energy_ev = radiation_drive_->incident_energy_ev;
        incident_charge_state_abs = radiation_drive_->incident_charge_state_abs;
    }

    if (incident_current_density_a_per_m2 < 0.0 || incident_energy_ev < 0.0 ||
        incident_charge_state_abs <= 0.0)
    {
        return false;
    }

    const auto& layers_before = transport_model_.getLayerStates();
    std::vector<double> previous_charge_density_c_per_m3;
    std::vector<double> previous_deposited_energy_j_per_m3;
    previous_charge_density_c_per_m3.reserve(layers_before.size());
    previous_deposited_energy_j_per_m3.reserve(layers_before.size());
    for (const auto& layer : layers_before)
    {
        previous_charge_density_c_per_m3.push_back(layer.charge_density_c_per_m3);
        previous_deposited_energy_j_per_m3.push_back(layer.deposited_energy_j_per_m3);
    }

    transport_model_.depositFlux(incident_current_density_a_per_m2, incident_energy_ev,
                                 incident_charge_state_abs, dt);
    updateRadiationInducedConductivity(dt, previous_charge_density_c_per_m3,
                                       previous_deposited_energy_j_per_m3,
                                       incident_charge_state_abs);
    updateFieldEstimate();

    const double max_field =
        status_.electric_field_v_per_m.empty()
            ? 0.0
            : *std::max_element(status_.electric_field_v_per_m.begin(), status_.electric_field_v_per_m.end());
    status_.max_electric_field_v_per_m = max_field;
    const auto assessment =
        breakdown_model_.assess(max_field, material_->getBreakdownFieldVPerM(), status_.time_s);
    status_.discharge_probability = assessment.probability;
    status_.breakdown_active = assessment.triggered;

    if (assessment.triggered)
    {
        auto relieved = discharge_model_.relieve(transport_model_.getLayerStates(), 0.35);
        status_.total_stored_energy_j =
            std::max(0.0, status_.total_stored_energy_j - relieved.released_energy_j_per_m2 * config_.area_m2);
        status_.discharge_events += 1;
        updateFieldEstimate();
    }

    status_.time_s += dt;
    return true;
}

void SpacecraftInternalChargingAlgorithm::reset()
{
    *this = SpacecraftInternalChargingAlgorithm{};
}

bool SpacecraftInternalChargingAlgorithm::exportResults(const std::filesystem::path& csv_path) const
{
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "depth_m";
    const auto& layers = transport_model_.getLayerStates();
    data_set.scalar_series["volume_charge_density_c_per_m3"] = {};
    data_set.scalar_series["electric_field_v_per_m"] = status_.electric_field_v_per_m;
    data_set.scalar_series["deposited_energy_j_per_m3"] = {};
    data_set.scalar_series["dose_gy"] = {};
    data_set.scalar_series["effective_conductivity_s_per_m"] = {};
    data_set.scalar_series["deposited_charge_rate_c_per_m3_s"] = {};
    data_set.scalar_series["deposited_particle_rate_per_m3_s"] = {};

    for (const auto& layer : layers)
    {
        data_set.axis_values.push_back(layer.depth_m);
        data_set.scalar_series["volume_charge_density_c_per_m3"].push_back(layer.charge_density_c_per_m3);
        data_set.scalar_series["deposited_energy_j_per_m3"].push_back(layer.deposited_energy_j_per_m3);
        data_set.scalar_series["dose_gy"].push_back(layer.dose_gy);
        data_set.scalar_series["effective_conductivity_s_per_m"].push_back(
            layer.effective_conductivity_s_per_m);
        data_set.scalar_series["deposited_charge_rate_c_per_m3_s"].push_back(
            layer.deposited_charge_rate_c_per_m3_s);
        data_set.scalar_series["deposited_particle_rate_per_m3_s"].push_back(
            layer.deposited_particle_rate_per_m3_s);
    }
    data_set.metadata["module"] = "Internal Charging";
    data_set.metadata["source_mode"] = toString(config_.source_mode);
    data_set.metadata["max_electric_field_v_per_m"] = std::to_string(status_.max_electric_field_v_per_m);
    data_set.metadata["breakdown_field_v_per_m"] =
        std::to_string(material_ != nullptr ? material_->getBreakdownFieldVPerM() : 0.0);
    data_set.metadata["average_dose_gy"] = std::to_string(status_.average_dose_gy);
    data_set.metadata["effective_conductivity_s_per_m"] =
        std::to_string(status_.effective_conductivity_s_per_m);
    data_set.metadata["deposited_charge_rate_c_per_m3_s"] =
        std::to_string(status_.deposited_charge_rate_c_per_m3_s);
    data_set.metadata["deposited_particle_rate_per_m3_s"] =
        std::to_string(status_.deposited_particle_rate_per_m3_s);
    data_set.metadata["effective_incident_charge_state_abs"] =
        std::to_string(status_.effective_incident_charge_state_abs);
    return static_cast<bool>(exporter_.exportDataSet(csv_path, data_set));
}

void SpacecraftInternalChargingAlgorithm::updateRadiationInducedConductivity(
    double dt, const std::vector<double>& previous_charge_density_c_per_m3,
    const std::vector<double>& previous_deposited_energy_j_per_m3,
    double incident_charge_state_abs)
{
    auto& layers = transport_model_.getLayerStates();
    if (layers.empty() || dt <= 0.0 || material_ == nullptr)
    {
        status_.average_dose_gy = 0.0;
        status_.effective_conductivity_s_per_m = 0.0;
        status_.deposited_charge_rate_c_per_m3_s = 0.0;
        status_.deposited_particle_rate_per_m3_s = 0.0;
        status_.effective_incident_charge_state_abs = 1.0;
        return;
    }

    const double effective_charge_state_abs = std::max(incident_charge_state_abs, 1.0e-12);

    const double eps_r = std::max(1.0, material_->getPermittivity());
    const double dielectric_permittivity_f_per_m = kEpsilon0 * eps_r;
    const double base_conductivity_s_per_m = std::max(0.0, material_->getConductivity());
    const double min_conductivity_s_per_m =
        std::max(1.0e-20, config_.min_effective_conductivity_s_per_m);
    const double max_conductivity_s_per_m =
        std::max(min_conductivity_s_per_m, config_.max_effective_conductivity_s_per_m);
    const double mass_density_kg_per_m3 =
        std::max(kMinMassDensityKgPerM3, material_->getMassDensityKgPerM3());

    double sum_dose_gy = 0.0;
    double sum_effective_conductivity_s_per_m = 0.0;
    double sum_charge_rate_c_per_m3_s = 0.0;
    double sum_particle_rate_per_m3_s = 0.0;

    for (std::size_t i = 0; i < layers.size(); ++i)
    {
        const double previous_charge =
            i < previous_charge_density_c_per_m3.size() ? previous_charge_density_c_per_m3[i] : 0.0;
        const double previous_energy =
            i < previous_deposited_energy_j_per_m3.size() ? previous_deposited_energy_j_per_m3[i] : 0.0;

        const double deposited_charge_delta_c_per_m3 =
            std::max(0.0, layers[i].charge_density_c_per_m3 - previous_charge);
        const double deposited_energy_delta_j_per_m3 =
            std::max(0.0, layers[i].deposited_energy_j_per_m3 - previous_energy);

        const double deposited_charge_rate_c_per_m3_s = deposited_charge_delta_c_per_m3 / dt;
        const double deposited_particle_rate_per_m3_s =
            deposited_charge_rate_c_per_m3_s / (kElementaryChargeC * effective_charge_state_abs);

        const double dose_gy = std::max(0.0, layers[i].deposited_energy_j_per_m3 / mass_density_kg_per_m3);
        const double dose_rate_gy_per_s =
            std::max(0.0, deposited_energy_delta_j_per_m3 / (mass_density_kg_per_m3 * dt));

        const double radiation_conductivity_s_per_m =
            config_.radiation_conductivity_charge_gain_s_per_m_per_c_per_m3 *
                std::abs(layers[i].charge_density_c_per_m3) +
            config_.radiation_conductivity_dose_gain_s_per_m_per_gy * dose_gy +
            config_.radiation_conductivity_dose_rate_gain_s2_per_m_per_gy * dose_rate_gy_per_s;
        const double effective_conductivity_s_per_m =
            std::clamp(base_conductivity_s_per_m + radiation_conductivity_s_per_m,
                       min_conductivity_s_per_m, max_conductivity_s_per_m);

        const double dielectric_relaxation_time_s =
            dielectric_permittivity_f_per_m /
            std::max(effective_conductivity_s_per_m, min_conductivity_s_per_m);
        const double relaxation_factor =
            std::exp(-dt / std::max(dielectric_relaxation_time_s, kMinTimeScaleS));

        // Radiation-induced conductivity drains trapped charge and moderates field growth.
        layers[i].charge_density_c_per_m3 *= relaxation_factor;

        layers[i].dose_gy = dose_gy;
        layers[i].effective_conductivity_s_per_m = effective_conductivity_s_per_m;
        layers[i].deposited_charge_rate_c_per_m3_s = deposited_charge_rate_c_per_m3_s;
        layers[i].deposited_particle_rate_per_m3_s = deposited_particle_rate_per_m3_s;

        sum_dose_gy += dose_gy;
        sum_effective_conductivity_s_per_m += effective_conductivity_s_per_m;
        sum_charge_rate_c_per_m3_s += deposited_charge_rate_c_per_m3_s;
        sum_particle_rate_per_m3_s += deposited_particle_rate_per_m3_s;
    }

    const double layer_count = static_cast<double>(layers.size());
    status_.average_dose_gy = sum_dose_gy / layer_count;
    status_.effective_conductivity_s_per_m = sum_effective_conductivity_s_per_m / layer_count;
    status_.deposited_charge_rate_c_per_m3_s = sum_charge_rate_c_per_m3_s / layer_count;
    status_.deposited_particle_rate_per_m3_s = sum_particle_rate_per_m3_s / layer_count;
    status_.effective_incident_charge_state_abs = effective_charge_state_abs;
}

void SpacecraftInternalChargingAlgorithm::updateFieldEstimate()
{
    const auto& layers = transport_model_.getLayerStates();
    status_.volume_charge_density_c_per_m3.clear();
    status_.electric_field_v_per_m.assign(layers.size(), 0.0);

    if (layers.empty())
    {
        status_.total_stored_energy_j = 0.0;
        return;
    }

    const double eps_r = std::max(1.0, material_->getPermittivity());
    const double layer_thickness_m = config_.thickness_m / static_cast<double>(layers.size());
    double cumulative_charge = 0.0;
    double stored_energy_j_per_m2 = 0.0;
    for (std::size_t i = 0; i < layers.size(); ++i)
    {
        cumulative_charge += layers[i].charge_density_c_per_m3 * layer_thickness_m;
        status_.electric_field_v_per_m[i] = cumulative_charge / (kEpsilon0 * eps_r);
        status_.volume_charge_density_c_per_m3.push_back(layers[i].charge_density_c_per_m3);
        const double energy_density_j_per_m3 =
            0.5 * kEpsilon0 * eps_r * status_.electric_field_v_per_m[i] * status_.electric_field_v_per_m[i];
        stored_energy_j_per_m2 += energy_density_j_per_m3 * layer_thickness_m;
    }
    status_.total_stored_energy_j = stored_energy_j_per_m2 * config_.area_m2;
}

std::string toString(InternalChargingSourceMode source_mode)
{
    switch (source_mode)
    {
    case InternalChargingSourceMode::Preset:
        return "preset";
    case InternalChargingSourceMode::Radiation:
        return "radiation";
    }
    return "unknown";
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
