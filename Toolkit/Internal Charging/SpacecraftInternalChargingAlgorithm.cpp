#include "SpacecraftInternalChargingAlgorithm.h"

#include <algorithm>
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
    if (!initialized_ || material_ == nullptr)
    {
        return false;
    }

    transport_model_.depositFlux(config_.incident_current_density_a_per_m2, config_.incident_energy_ev, dt);
    updateFieldEstimate();

    const double max_field =
        status_.electric_field_v_per_m.empty()
            ? 0.0
            : *std::max_element(status_.electric_field_v_per_m.begin(), status_.electric_field_v_per_m.end());
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

    for (const auto& layer : layers)
    {
        data_set.axis_values.push_back(layer.depth_m);
        data_set.scalar_series["volume_charge_density_c_per_m3"].push_back(layer.charge_density_c_per_m3);
    }
    data_set.metadata["module"] = "Internal Charging";
    return static_cast<bool>(exporter_.exportDataSet(csv_path, data_set));
}

void SpacecraftInternalChargingAlgorithm::updateFieldEstimate()
{
    const auto& layers = transport_model_.getLayerStates();
    status_.volume_charge_density_c_per_m3.clear();
    status_.electric_field_v_per_m.assign(layers.size(), 0.0);

    const double eps_r = std::max(1.0, material_->getPermittivity());
    double cumulative_charge = 0.0;
    double stored_energy_density = 0.0;
    for (std::size_t i = 0; i < layers.size(); ++i)
    {
        cumulative_charge += layers[i].charge_density_c_per_m3 * (config_.thickness_m / config_.layers);
        status_.electric_field_v_per_m[i] = cumulative_charge / (kEpsilon0 * eps_r);
        status_.volume_charge_density_c_per_m3.push_back(layers[i].charge_density_c_per_m3);
        stored_energy_density += 0.5 * kEpsilon0 * eps_r * status_.electric_field_v_per_m[i] *
                                 status_.electric_field_v_per_m[i];
    }
    status_.total_stored_energy_j = stored_energy_density * config_.area_m2 * config_.thickness_m;
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
