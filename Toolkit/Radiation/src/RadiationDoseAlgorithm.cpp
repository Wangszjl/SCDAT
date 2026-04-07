#include "RadiationDoseAlgorithm.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace Radiation
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kNumericalFloor = 1.0e-12;
constexpr double kElectronMassMeV = 0.51099895;
constexpr double kProtonMassMeV = 938.2720813;
constexpr double kAlphaMassMeV = 3727.379378;
constexpr double kBetheKMeVCm2PerG = 0.307075;
constexpr double kTrackTransportEnergyFloorEv = 20.0;

std::string normalizeToken(std::string text)
{
    std::string token;
    token.reserve(text.size());
    for (const char ch : text)
    {
        if (!std::isspace(static_cast<unsigned char>(ch)) && ch != '_' && ch != '-')
        {
            token.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
        }
    }
    return token;
}
}

bool RadiationDoseAlgorithm::initialize(const RadiationConfiguration& config)
{
    if (config.layers == 0 || config.thickness_m <= 0.0 || config.area_m2 <= 0.0 ||
        config.particle_flux_m2_s < 0.0 || config.mean_energy_ev < 0.0 ||
        config.monte_carlo_histories_per_step == 0 ||
        config.monte_carlo_max_steps_per_track == 0 ||
        config.monte_carlo_max_recorded_points == 0 ||
        config.monte_carlo_lateral_span_m <= 0.0 ||
        config.electron_process.elastic_scattering_scale < 0.0 ||
        config.electron_process.inelastic_energy_loss_scale <= 0.0 ||
        config.electron_process.secondary_energy_fraction < 0.0 ||
        config.electron_process.secondary_energy_fraction > 0.45 ||
        config.electron_process.secondary_depth_attenuation_fraction <= 0.0 ||
        config.proton_process.nuclear_reaction_energy_scale <= 0.0 ||
        config.proton_process.elastic_scattering_scale < 0.0 ||
        config.proton_process.secondary_energy_fraction < 0.0 ||
        config.proton_process.secondary_energy_fraction > 0.40 ||
        config.proton_process.secondary_depth_buildup_fraction <= 0.0 ||
        config.heavy_ion_process.effective_charge_scale <= 0.0 ||
        config.heavy_ion_process.stopping_power_scale <= 0.0 ||
        config.heavy_ion_process.elastic_scattering_scale < 0.0 ||
        config.heavy_ion_process.fragmentation_secondary_energy_fraction < 0.0 ||
        config.heavy_ion_process.fragmentation_secondary_energy_fraction > 0.55 ||
        config.heavy_ion_process.fragmentation_depth_buildup_fraction <= 0.0 ||
        config.heavy_ion_process.fragmentation_event_rate_scale < 0.0)
    {
        return false;
    }

    config_ = config;
    material_ = material_database_.findByAliasOrName(config.material_name);
    if (material_ == nullptr)
    {
        return false;
    }
    resolved_material_name_ = material_->getName();

    stopping_data_resolution_ =
        stopping_data_governance_.resolve(config_.particle_species, config_.physics_list,
                                          resolved_material_name_,
                                          config_.stopping_data_version,
                                          config_.stopping_data_allow_rollback);
    if (!stopping_data_resolution_.valid)
    {
        return false;
    }

    layer_thickness_m_ = config_.thickness_m / static_cast<double>(config_.layers);

    status_ = RadiationStatus{};
    status_.layers.assign(config_.layers, {});
    for (std::size_t i = 0; i < config_.layers; ++i)
    {
        status_.layers[i].depth_m = (static_cast<double>(i) + 0.5) * layer_thickness_m_;
    }

    next_track_id_ = 0;
    random_engine_.seed(config_.monte_carlo_seed);
    electron_process_inelastic_accumulator_ = 0.0;
    electron_process_scattering_accumulator_ = 0.0;
    electron_process_sample_count_ = 0.0;
    proton_process_nuclear_accumulator_ = 0.0;
    proton_process_scattering_accumulator_ = 0.0;
    proton_process_sample_count_ = 0.0;
    heavy_ion_process_effective_charge_accumulator_ = 0.0;
    heavy_ion_process_stopping_accumulator_ = 0.0;
    heavy_ion_process_scattering_accumulator_ = 0.0;
    heavy_ion_process_sample_count_ = 0.0;
    track_scattering_sigma_accumulator_ = 0.0;
    track_terminal_energy_accumulator_ev_ = 0.0;

    initialized_ = true;
    return true;
}

bool RadiationDoseAlgorithm::advance(double dt)
{
    if (!initialized_ || dt <= 0.0)
    {
        return false;
    }

    const double incident_step_energy_j_per_m2 =
        config_.particle_flux_m2_s * config_.mean_energy_ev * kElementaryCharge * dt;
    double deposited_step_energy_j_per_m2 = 0.0;

    if (config_.enable_monte_carlo_transport)
    {
        if (!advanceMonteCarlo(dt, incident_step_energy_j_per_m2, deposited_step_energy_j_per_m2))
        {
            return false;
        }
    }
    else
    {
        const bool electron_process_enabled =
            (config_.particle_species == ParticleSpecies::Electron) && config_.electron_process.enable;
        const bool proton_process_enabled =
            (config_.particle_species == ParticleSpecies::Proton) && config_.proton_process.enable;
        const bool heavy_ion_process_enabled =
            (config_.particle_species == ParticleSpecies::HeavyIon) && config_.heavy_ion_process.enable;
        const double electron_inelastic_multiplier =
            electron_process_enabled ? electronInelasticEnergyMultiplier(config_.mean_energy_ev) : 1.0;
        const double proton_nuclear_multiplier =
            proton_process_enabled ? protonNuclearEnergyMultiplier(config_.mean_energy_ev) : 1.0;
        const double heavy_ion_effective_charge_multiplier =
            heavy_ion_process_enabled ? heavyIonEffectiveChargeMultiplier(config_.mean_energy_ev) : 1.0;
        const double heavy_ion_stopping_multiplier =
            heavy_ion_process_enabled ? heavyIonStoppingPowerMultiplier(config_.mean_energy_ev) : 1.0;

        deposited_step_energy_j_per_m2 = incident_step_energy_j_per_m2 * depositionEfficiency();
        if (electron_process_enabled)
        {
            deposited_step_energy_j_per_m2 = std::min(
                incident_step_energy_j_per_m2,
                deposited_step_energy_j_per_m2 * electron_inelastic_multiplier);
        }
        if (proton_process_enabled)
        {
            deposited_step_energy_j_per_m2 = std::min(
                incident_step_energy_j_per_m2,
                deposited_step_energy_j_per_m2 * proton_nuclear_multiplier);
        }
        if (heavy_ion_process_enabled)
        {
            deposited_step_energy_j_per_m2 = std::min(
                incident_step_energy_j_per_m2,
                deposited_step_energy_j_per_m2 *
                    heavy_ion_effective_charge_multiplier * heavy_ion_stopping_multiplier);
        }

        const double base_attenuation_length = std::max(
            layer_thickness_m_, config_.thickness_m * std::max(0.05, config_.attenuation_length_fraction));
        // Use physics-list scattering scale to shape depth profile while preserving
        // deterministic energy-conservation behavior in the surrogate model.
        const double attenuation_length =
            std::max(layer_thickness_m_,
                     base_attenuation_length /
                         std::max(0.5, physicsListScatteringScale()));

        double weight_sum = 0.0;
        for (const auto& layer : status_.layers)
        {
            double weight = std::exp(-layer.depth_m / attenuation_length);
            if (electron_process_enabled)
            {
                const double depth_ratio = 1.0 -
                    std::clamp(layer.depth_m / std::max(config_.thickness_m, kNumericalFloor), 0.0, 1.0);
                const double secondary_shape =
                    1.0 + config_.electron_process.secondary_energy_fraction *
                              electronSecondaryDepthWeight(layer.depth_m);
                const double inelastic_shape =
                    1.0 + 0.35 * std::max(0.0, electron_inelastic_multiplier - 1.0) * depth_ratio;
                weight *= secondary_shape * inelastic_shape;
            }
            if (proton_process_enabled)
            {
                const double depth_ratio =
                    std::clamp(layer.depth_m / std::max(config_.thickness_m, kNumericalFloor), 0.0, 1.0);
                const double secondary_shape =
                    1.0 + config_.proton_process.secondary_energy_fraction *
                              protonSecondaryDepthWeight(layer.depth_m);
                const double nuclear_shape =
                    1.0 + 0.32 * std::max(0.0, proton_nuclear_multiplier - 1.0) * depth_ratio;
                weight *= secondary_shape * nuclear_shape;
            }
            if (heavy_ion_process_enabled)
            {
                const double depth_ratio =
                    std::clamp(layer.depth_m / std::max(config_.thickness_m, kNumericalFloor), 0.0, 1.0);
                const double secondary_shape =
                    1.0 + config_.heavy_ion_process.fragmentation_secondary_energy_fraction *
                              heavyIonFragmentationDepthWeight(layer.depth_m);
                const double stopping_shape =
                    1.0 + 0.34 * std::max(0.0, heavy_ion_stopping_multiplier - 1.0) *
                              (0.25 + 0.75 * depth_ratio);
                const double charge_shape =
                    1.0 + 0.24 * std::max(0.0, heavy_ion_effective_charge_multiplier - 1.0) * depth_ratio;
                weight *= secondary_shape * stopping_shape * charge_shape;
            }
            weight_sum += weight;
        }
        weight_sum = std::max(weight_sum, kNumericalFloor);

        const double density_kg_per_m3 = effectiveMassDensityKgPerM3();
        for (auto& layer : status_.layers)
        {
            double weight = std::exp(-layer.depth_m / attenuation_length);
            if (electron_process_enabled)
            {
                const double depth_ratio = 1.0 -
                    std::clamp(layer.depth_m / std::max(config_.thickness_m, kNumericalFloor), 0.0, 1.0);
                const double secondary_shape =
                    1.0 + config_.electron_process.secondary_energy_fraction *
                              electronSecondaryDepthWeight(layer.depth_m);
                const double inelastic_shape =
                    1.0 + 0.35 * std::max(0.0, electron_inelastic_multiplier - 1.0) * depth_ratio;
                weight *= secondary_shape * inelastic_shape;
            }
            if (proton_process_enabled)
            {
                const double depth_ratio =
                    std::clamp(layer.depth_m / std::max(config_.thickness_m, kNumericalFloor), 0.0, 1.0);
                const double secondary_shape =
                    1.0 + config_.proton_process.secondary_energy_fraction *
                              protonSecondaryDepthWeight(layer.depth_m);
                const double nuclear_shape =
                    1.0 + 0.32 * std::max(0.0, proton_nuclear_multiplier - 1.0) * depth_ratio;
                weight *= secondary_shape * nuclear_shape;
            }
            if (heavy_ion_process_enabled)
            {
                const double depth_ratio =
                    std::clamp(layer.depth_m / std::max(config_.thickness_m, kNumericalFloor), 0.0, 1.0);
                const double secondary_shape =
                    1.0 + config_.heavy_ion_process.fragmentation_secondary_energy_fraction *
                              heavyIonFragmentationDepthWeight(layer.depth_m);
                const double stopping_shape =
                    1.0 + 0.34 * std::max(0.0, heavy_ion_stopping_multiplier - 1.0) *
                              (0.25 + 0.75 * depth_ratio);
                const double charge_shape =
                    1.0 + 0.24 * std::max(0.0, heavy_ion_effective_charge_multiplier - 1.0) * depth_ratio;
                weight *= secondary_shape * stopping_shape * charge_shape;
            }
            const double deposited_layer_energy_j_per_m2 =
                deposited_step_energy_j_per_m2 * weight / weight_sum;
            const double deposited_layer_energy_j_per_m3 =
                deposited_layer_energy_j_per_m2 / std::max(layer_thickness_m_, kNumericalFloor);

            layer.deposited_energy_j_per_m3 += deposited_layer_energy_j_per_m3;
            layer.dose_gy = layer.deposited_energy_j_per_m3 / density_kg_per_m3;

            if (electron_process_enabled)
            {
                const double secondary_shape =
                    config_.electron_process.secondary_energy_fraction *
                    electronSecondaryDepthWeight(layer.depth_m);
                const double secondary_component_j_per_m2 =
                    deposited_layer_energy_j_per_m2 * secondary_shape /
                    std::max(1.0 + secondary_shape, kNumericalFloor);
                status_.electron_process_secondary_energy_j_per_m2 += secondary_component_j_per_m2;
            }
            if (proton_process_enabled)
            {
                const double secondary_shape =
                    config_.proton_process.secondary_energy_fraction *
                    protonSecondaryDepthWeight(layer.depth_m);
                const double secondary_component_j_per_m2 =
                    deposited_layer_energy_j_per_m2 * secondary_shape /
                    std::max(1.0 + secondary_shape, kNumericalFloor);
                status_.proton_process_secondary_energy_j_per_m2 += secondary_component_j_per_m2;
            }
            if (heavy_ion_process_enabled)
            {
                const double secondary_shape =
                    config_.heavy_ion_process.fragmentation_secondary_energy_fraction *
                    heavyIonFragmentationDepthWeight(layer.depth_m);
                const double secondary_component_j_per_m2 =
                    deposited_layer_energy_j_per_m2 * secondary_shape /
                    std::max(1.0 + secondary_shape, kNumericalFloor);
                status_.heavy_ion_process_secondary_energy_j_per_m2 += secondary_component_j_per_m2;
            }
        }

        if (electron_process_enabled)
        {
            status_.electron_process_secondary_event_count +=
                config_.particle_flux_m2_s * config_.area_m2 * dt *
                std::max(0.0, config_.electron_process.secondary_energy_fraction);
            updateElectronProcessRunningMeans(
                electron_inelastic_multiplier,
                electronElasticScatteringSigmaRad(config_.mean_energy_ev, layer_thickness_m_));
        }
        if (proton_process_enabled)
        {
            status_.proton_process_secondary_event_count +=
                config_.particle_flux_m2_s * config_.area_m2 * dt *
                std::max(0.0, config_.proton_process.secondary_energy_fraction);
            updateProtonProcessRunningMeans(
                proton_nuclear_multiplier,
                protonElasticScatteringSigmaRad(config_.mean_energy_ev, layer_thickness_m_));
        }
        if (heavy_ion_process_enabled)
        {
            status_.heavy_ion_process_fragmentation_event_count +=
                config_.particle_flux_m2_s * config_.area_m2 * dt *
                std::max(0.0, config_.heavy_ion_process.fragmentation_secondary_energy_fraction) *
                std::max(0.0, config_.heavy_ion_process.fragmentation_event_rate_scale);
            updateHeavyIonProcessRunningMeans(
                heavy_ion_effective_charge_multiplier,
                heavy_ion_stopping_multiplier,
                heavyIonElasticScatteringSigmaRad(config_.mean_energy_ev, layer_thickness_m_));
        }
    }

    updateConservationMetrics(incident_step_energy_j_per_m2, deposited_step_energy_j_per_m2);
    status_.time_s += dt;

    return true;
}

void RadiationDoseAlgorithm::reset()
{
    *this = RadiationDoseAlgorithm{};
}

bool RadiationDoseAlgorithm::exportResults(const std::filesystem::path& csv_path) const
{
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "depth_m";
    data_set.scalar_series["deposited_energy_j_per_m3"] = {};
    data_set.scalar_series["dose_gy"] = {};

    for (const auto& layer : status_.layers)
    {
        data_set.axis_values.push_back(layer.depth_m);
        data_set.scalar_series["deposited_energy_j_per_m3"].push_back(layer.deposited_energy_j_per_m3);
        data_set.scalar_series["dose_gy"].push_back(layer.dose_gy);
    }

    data_set.metadata["module"] = "Radiation";
    data_set.metadata["particle_species"] = toString(config_.particle_species);
    data_set.metadata["physics_list"] = toString(config_.physics_list);
    data_set.metadata["material"] = config_.material_name;
    data_set.metadata["material_resolved"] = resolved_material_name_;
    data_set.metadata["time_s"] = std::to_string(status_.time_s);
    data_set.metadata["incident_energy_j_per_m2"] = std::to_string(status_.incident_energy_j_per_m2);
    data_set.metadata["deposited_energy_j_per_m2"] = std::to_string(status_.deposited_energy_j_per_m2);
    data_set.metadata["escaped_energy_j_per_m2"] = std::to_string(status_.escaped_energy_j_per_m2);
    data_set.metadata["energy_conservation_error"] = std::to_string(status_.energy_conservation_error);
    data_set.metadata["stopping_model"] = "geant4-aligned-bethe-bragg-migration";
    data_set.metadata["material_z_over_a"] = std::to_string(effectiveMaterialZOverA());
    data_set.metadata["mean_excitation_ev"] = std::to_string(effectiveMeanExcitationEnergyEv());
    data_set.metadata["mass_stopping_power_mev_cm2_per_g"] =
        std::to_string(geant4AlignedMassStoppingPowerMeVcm2PerG(config_.mean_energy_ev));
    data_set.metadata["stopping_data_requested_version"] =
        stopping_data_resolution_.requested_version.empty()
            ? "latest"
            : stopping_data_resolution_.requested_version;
    data_set.metadata["stopping_data_resolved_version"] =
        stopping_data_resolution_.resolved_version;
    data_set.metadata["stopping_data_dataset_id"] = stopping_data_resolution_.dataset_id;
    data_set.metadata["stopping_data_dataset_source_uri"] =
        stopping_data_resolution_.dataset_source_uri;
    data_set.metadata["stopping_data_dataset_created_utc"] =
        stopping_data_resolution_.dataset_created_utc;
    data_set.metadata["stopping_data_material_source_hash"] =
        stopping_data_resolution_.material_source_hash;
    data_set.metadata["stopping_data_physics_source_hash"] =
        stopping_data_resolution_.physics_source_hash;
    data_set.metadata["stopping_data_cache_hit"] =
        stopping_data_resolution_.cache_hit ? "true" : "false";
    data_set.metadata["stopping_data_rollback_applied"] =
        stopping_data_resolution_.rollback_applied ? "true" : "false";
    data_set.metadata["stopping_data_rollback_reason"] =
        stopping_data_resolution_.rollback_reason;
    data_set.metadata["stopping_data_validation_error"] =
        stopping_data_resolution_.validation_error;
    data_set.metadata["enable_monte_carlo_transport"] =
        config_.enable_monte_carlo_transport ? "true" : "false";
    data_set.metadata["electron_process_enabled"] =
        config_.electron_process.enable ? "true" : "false";
    data_set.metadata["electron_process_elastic_scattering_scale"] =
        std::to_string(config_.electron_process.elastic_scattering_scale);
    data_set.metadata["electron_process_inelastic_energy_loss_scale"] =
        std::to_string(config_.electron_process.inelastic_energy_loss_scale);
    data_set.metadata["electron_process_secondary_energy_fraction"] =
        std::to_string(config_.electron_process.secondary_energy_fraction);
    data_set.metadata["electron_process_secondary_depth_attenuation_fraction"] =
        std::to_string(config_.electron_process.secondary_depth_attenuation_fraction);
    data_set.metadata["electron_process_secondary_energy_j_per_m2"] =
        std::to_string(status_.electron_process_secondary_energy_j_per_m2);
    data_set.metadata["electron_process_mean_inelastic_multiplier"] =
        std::to_string(status_.electron_process_mean_inelastic_multiplier);
    data_set.metadata["electron_process_mean_scattering_sigma_rad"] =
        std::to_string(status_.electron_process_mean_scattering_sigma_rad);
    data_set.metadata["electron_process_secondary_event_count"] =
        std::to_string(status_.electron_process_secondary_event_count);
    data_set.metadata["proton_process_enabled"] =
        config_.proton_process.enable ? "true" : "false";
    data_set.metadata["proton_process_nuclear_reaction_energy_scale"] =
        std::to_string(config_.proton_process.nuclear_reaction_energy_scale);
    data_set.metadata["proton_process_elastic_scattering_scale"] =
        std::to_string(config_.proton_process.elastic_scattering_scale);
    data_set.metadata["proton_process_secondary_energy_fraction"] =
        std::to_string(config_.proton_process.secondary_energy_fraction);
    data_set.metadata["proton_process_secondary_depth_buildup_fraction"] =
        std::to_string(config_.proton_process.secondary_depth_buildup_fraction);
    data_set.metadata["proton_process_secondary_energy_j_per_m2"] =
        std::to_string(status_.proton_process_secondary_energy_j_per_m2);
    data_set.metadata["proton_process_mean_nuclear_multiplier"] =
        std::to_string(status_.proton_process_mean_nuclear_multiplier);
    data_set.metadata["proton_process_mean_scattering_sigma_rad"] =
        std::to_string(status_.proton_process_mean_scattering_sigma_rad);
    data_set.metadata["proton_process_secondary_event_count"] =
        std::to_string(status_.proton_process_secondary_event_count);
    data_set.metadata["heavy_ion_process_enabled"] =
        config_.heavy_ion_process.enable ? "true" : "false";
    data_set.metadata["heavy_ion_process_effective_charge_scale"] =
        std::to_string(config_.heavy_ion_process.effective_charge_scale);
    data_set.metadata["heavy_ion_process_stopping_power_scale"] =
        std::to_string(config_.heavy_ion_process.stopping_power_scale);
    data_set.metadata["heavy_ion_process_elastic_scattering_scale"] =
        std::to_string(config_.heavy_ion_process.elastic_scattering_scale);
    data_set.metadata["heavy_ion_process_fragmentation_secondary_energy_fraction"] =
        std::to_string(config_.heavy_ion_process.fragmentation_secondary_energy_fraction);
    data_set.metadata["heavy_ion_process_fragmentation_depth_buildup_fraction"] =
        std::to_string(config_.heavy_ion_process.fragmentation_depth_buildup_fraction);
    data_set.metadata["heavy_ion_process_fragmentation_event_rate_scale"] =
        std::to_string(config_.heavy_ion_process.fragmentation_event_rate_scale);
    data_set.metadata["heavy_ion_process_secondary_energy_j_per_m2"] =
        std::to_string(status_.heavy_ion_process_secondary_energy_j_per_m2);
    data_set.metadata["heavy_ion_process_mean_effective_charge_multiplier"] =
        std::to_string(status_.heavy_ion_process_mean_effective_charge_multiplier);
    data_set.metadata["heavy_ion_process_mean_stopping_multiplier"] =
        std::to_string(status_.heavy_ion_process_mean_stopping_multiplier);
    data_set.metadata["heavy_ion_process_mean_scattering_sigma_rad"] =
        std::to_string(status_.heavy_ion_process_mean_scattering_sigma_rad);
    data_set.metadata["heavy_ion_process_fragmentation_event_count"] =
        std::to_string(status_.heavy_ion_process_fragmentation_event_count);
    data_set.metadata["track_count"] = std::to_string(status_.track_count);
    data_set.metadata["track_schema_version"] = "v2";
    data_set.metadata["track_step_count"] = std::to_string(status_.track_step_count);
    data_set.metadata["track_total_path_length_m"] = std::to_string(status_.track_total_path_length_m);
    data_set.metadata["track_mean_step_length_m"] = std::to_string(status_.track_mean_step_length_m);
    data_set.metadata["track_mean_scattering_sigma_rad"] =
        std::to_string(status_.track_mean_scattering_sigma_rad);
    data_set.metadata["track_secondary_event_count"] =
        std::to_string(status_.track_secondary_event_count);
    data_set.metadata["track_secondary_energy_j_per_m2"] =
        std::to_string(status_.track_secondary_energy_j_per_m2);
    data_set.metadata["track_secondary_yield_per_primary"] =
        std::to_string(status_.track_secondary_yield_per_primary);
    data_set.metadata["track_mean_terminal_energy_ev"] =
        std::to_string(status_.track_mean_terminal_energy_ev);
    data_set.metadata["track_max_depth_m"] = std::to_string(status_.track_max_depth_m);
    data_set.metadata["recorded_track_points"] = std::to_string(status_.track_points.size());
    data_set.metadata["dropped_track_points"] = std::to_string(status_.dropped_track_points);

    if (!static_cast<bool>(exporter_.exportDataSet(csv_path, data_set)))
    {
        return false;
    }

    if (config_.enable_monte_carlo_transport)
    {
        return exportTrackCsv(csv_path);
    }
    return true;
}

bool RadiationDoseAlgorithm::advanceMonteCarlo(double dt, double incident_step_energy_j_per_m2,
                                               double& deposited_step_energy_j_per_m2)
{
    deposited_step_energy_j_per_m2 = 0.0;
    if (status_.layers.empty())
    {
        return true;
    }

    const double total_particles = config_.particle_flux_m2_s * config_.area_m2 * dt;
    if (total_particles <= 0.0 || config_.mean_energy_ev <= 0.0)
    {
        return true;
    }

    const std::size_t history_count = std::max<std::size_t>(1, config_.monte_carlo_histories_per_step);
    const double particle_weight = total_particles / static_cast<double>(history_count);
    const bool electron_process_enabled =
        (config_.particle_species == ParticleSpecies::Electron) && config_.electron_process.enable;
    const bool proton_process_enabled =
        (config_.particle_species == ParticleSpecies::Proton) && config_.proton_process.enable;
    const bool heavy_ion_process_enabled =
        (config_.particle_species == ParticleSpecies::HeavyIon) && config_.heavy_ion_process.enable;

    const double density_kg_per_m3 = effectiveMassDensityKgPerM3();
    const double density_g_per_cm3 = density_kg_per_m3 / 1000.0;

    const double half_lateral_span_m = 0.5 * config_.monte_carlo_lateral_span_m;
    constexpr double kPi = 3.14159265358979323846;
    const double angular_sigma_rad =
        std::max(1.0e-6, config_.monte_carlo_angular_sigma_deg * kPi / 180.0);
    std::uniform_real_distribution<double> lateral_distribution(-half_lateral_span_m,
                                                                half_lateral_span_m);
    std::normal_distribution<double> angle_distribution(0.0, angular_sigma_rad);
    std::uniform_real_distribution<double> step_jitter_distribution(0.65, 1.35);

    for (std::size_t history = 0; history < history_count; ++history)
    {
        const std::size_t track_id = next_track_id_++;
        double x_m = lateral_distribution(random_engine_);
        double z_m = 0.0;
        double theta_rad = angle_distribution(random_engine_);
        double energy_ev = config_.mean_energy_ev;

        for (std::size_t step = 0;
             step < config_.monte_carlo_max_steps_per_track && z_m < config_.thickness_m &&
             energy_ev > kTrackTransportEnergyFloorEv;
             ++step)
        {
            const double remaining_depth_m = std::max(0.0, config_.thickness_m - z_m);
            if (remaining_depth_m <= 0.0)
            {
                break;
            }

            const double base_step_m = std::max(layer_thickness_m_ * 0.25, 1.0e-7);
            double step_m = std::min(remaining_depth_m, base_step_m * step_jitter_distribution(random_engine_));
            step_m = std::max(step_m, std::min(remaining_depth_m, 1.0e-7));

            const double stopping_power_mev_cm2_per_g =
                geant4AlignedMassStoppingPowerMeVcm2PerG(energy_ev) *
                physicsListStoppingScale();
            const double d_edx_mev_per_cm = stopping_power_mev_cm2_per_g * density_g_per_cm3;
            double deposited_energy_ev = std::max(0.0, d_edx_mev_per_cm * step_m * 100.0 * 1.0e6);

            const double straggling_sigma_ev =
                0.20 * physicsListScatteringScale() * deposited_energy_ev;
            std::normal_distribution<double> straggling_distribution(0.0, straggling_sigma_ev);
            deposited_energy_ev += straggling_distribution(random_engine_);
            deposited_energy_ev = std::clamp(deposited_energy_ev, 0.0, energy_ev);

            const double mid_depth_m = std::min(config_.thickness_m, z_m + 0.5 * step_m);
            double electron_scattering_sigma_rad = 0.0;
            double proton_scattering_sigma_rad = 0.0;
            double heavy_ion_scattering_sigma_rad = 0.0;
            double total_secondary_energy_ev = 0.0;
            if (electron_process_enabled)
            {
                const double inelastic_multiplier = electronInelasticEnergyMultiplier(energy_ev);
                deposited_energy_ev = std::min(energy_ev, deposited_energy_ev * inelastic_multiplier);

                const double secondary_energy_ev = std::min(
                    std::max(0.0, energy_ev - deposited_energy_ev),
                    electronSecondaryEnergyEv(deposited_energy_ev, mid_depth_m));
                deposited_energy_ev = std::min(energy_ev, deposited_energy_ev + secondary_energy_ev);

                if (secondary_energy_ev > 0.0)
                {
                    status_.electron_process_secondary_event_count += particle_weight;
                    status_.electron_process_secondary_energy_j_per_m2 +=
                        secondary_energy_ev * kElementaryCharge * particle_weight /
                        std::max(config_.area_m2, kNumericalFloor);
                    total_secondary_energy_ev += secondary_energy_ev;
                }

                electron_scattering_sigma_rad =
                    electronElasticScatteringSigmaRad(energy_ev, step_m);
                updateElectronProcessRunningMeans(inelastic_multiplier,
                                                  electron_scattering_sigma_rad);
            }
            if (proton_process_enabled)
            {
                const double nuclear_multiplier = protonNuclearEnergyMultiplier(energy_ev);
                deposited_energy_ev = std::min(energy_ev, deposited_energy_ev * nuclear_multiplier);

                const double secondary_energy_ev = std::min(
                    std::max(0.0, energy_ev - deposited_energy_ev),
                    protonSecondaryEnergyEv(deposited_energy_ev, mid_depth_m));
                deposited_energy_ev = std::min(energy_ev, deposited_energy_ev + secondary_energy_ev);

                if (secondary_energy_ev > 0.0)
                {
                    status_.proton_process_secondary_event_count += particle_weight;
                    status_.proton_process_secondary_energy_j_per_m2 +=
                        secondary_energy_ev * kElementaryCharge * particle_weight /
                        std::max(config_.area_m2, kNumericalFloor);
                    total_secondary_energy_ev += secondary_energy_ev;
                }

                proton_scattering_sigma_rad =
                    protonElasticScatteringSigmaRad(energy_ev, step_m);
                updateProtonProcessRunningMeans(nuclear_multiplier,
                                                proton_scattering_sigma_rad);
            }
            if (heavy_ion_process_enabled)
            {
                const double effective_charge_multiplier = heavyIonEffectiveChargeMultiplier(energy_ev);
                const double stopping_multiplier = heavyIonStoppingPowerMultiplier(energy_ev);
                deposited_energy_ev = std::min(
                    energy_ev,
                    deposited_energy_ev * effective_charge_multiplier * stopping_multiplier);

                const double secondary_energy_ev = std::min(
                    std::max(0.0, energy_ev - deposited_energy_ev),
                    heavyIonFragmentationSecondaryEnergyEv(deposited_energy_ev, mid_depth_m));
                deposited_energy_ev = std::min(energy_ev, deposited_energy_ev + secondary_energy_ev);

                if (secondary_energy_ev > 0.0)
                {
                    status_.heavy_ion_process_fragmentation_event_count +=
                        particle_weight * std::max(0.0, config_.heavy_ion_process.fragmentation_event_rate_scale);
                    status_.heavy_ion_process_secondary_energy_j_per_m2 +=
                        secondary_energy_ev * kElementaryCharge * particle_weight /
                        std::max(config_.area_m2, kNumericalFloor);
                    total_secondary_energy_ev += secondary_energy_ev;
                }

                heavy_ion_scattering_sigma_rad =
                    heavyIonElasticScatteringSigmaRad(energy_ev, step_m);
                updateHeavyIonProcessRunningMeans(effective_charge_multiplier,
                                                  stopping_multiplier,
                                                  heavy_ion_scattering_sigma_rad);
            }

            const std::size_t layer_index = layerIndexForDepth(mid_depth_m);

            const double deposited_energy_j = deposited_energy_ev * kElementaryCharge * particle_weight;
            const double deposited_energy_j_per_m2 =
                deposited_energy_j / std::max(config_.area_m2, kNumericalFloor);
            const double deposited_energy_j_per_m3 =
                deposited_energy_j_per_m2 / std::max(layer_thickness_m_, kNumericalFloor);
            status_.layers[layer_index].deposited_energy_j_per_m3 += deposited_energy_j_per_m3;
            deposited_step_energy_j_per_m2 += deposited_energy_j_per_m2;

            energy_ev -= deposited_energy_ev;

            double scattering_sigma_rad =
                std::clamp(0.08 * std::sqrt(step_m / std::max(config_.thickness_m, kNumericalFloor)) *
                               std::sqrt(1.0e5 / std::max(energy_ev, 1.0)) *
                               physicsListScatteringScale(),
                           1.0e-4, 0.25);

            if (electron_process_enabled)
            {
                scattering_sigma_rad =
                    std::clamp(std::hypot(scattering_sigma_rad, electron_scattering_sigma_rad),
                               1.0e-4, 0.35);
            }
            if (proton_process_enabled)
            {
                scattering_sigma_rad =
                    std::clamp(std::hypot(scattering_sigma_rad, proton_scattering_sigma_rad),
                               1.0e-4, 0.35);
            }
            if (heavy_ion_process_enabled)
            {
                scattering_sigma_rad =
                    std::clamp(std::hypot(scattering_sigma_rad, heavy_ion_scattering_sigma_rad),
                               1.0e-4, 0.35);
            }

            std::normal_distribution<double> scattering_distribution(0.0, scattering_sigma_rad);
            theta_rad = std::clamp(theta_rad + scattering_distribution(random_engine_), -1.3, 1.3);

            const double previous_x_m = x_m;
            const double previous_z_m = z_m;
            z_m = std::min(config_.thickness_m, z_m + step_m);
            x_m += std::tan(theta_rad) * step_m;

            const double step_length_m =
                std::hypot(x_m - previous_x_m, z_m - previous_z_m);
            status_.track_step_count += 1;
            status_.track_total_path_length_m += step_length_m;
            status_.track_max_depth_m = std::max(status_.track_max_depth_m, z_m);
            track_scattering_sigma_accumulator_ += scattering_sigma_rad;
            status_.track_mean_step_length_m =
                status_.track_total_path_length_m /
                std::max(1.0, static_cast<double>(status_.track_step_count));
            status_.track_mean_scattering_sigma_rad =
                track_scattering_sigma_accumulator_ /
                std::max(1.0, static_cast<double>(status_.track_step_count));

            if (total_secondary_energy_ev > 0.0)
            {
                status_.track_secondary_event_count += particle_weight;
                status_.track_secondary_energy_j_per_m2 +=
                    total_secondary_energy_ev * kElementaryCharge * particle_weight /
                    std::max(config_.area_m2, kNumericalFloor);
            }

            appendTrackPoint(RadiationStatus::TrackPoint{track_id,
                                                         step,
                                                         status_.time_s + dt,
                                                         x_m,
                                                         z_m,
                                                         energy_ev,
                                                         deposited_energy_ev,
                                                         total_secondary_energy_ev,
                                                         scattering_sigma_rad,
                                                         step_length_m,
                                                         particle_weight});
        }

        track_terminal_energy_accumulator_ev_ += std::max(0.0, energy_ev);
    }

    deposited_step_energy_j_per_m2 =
        std::clamp(deposited_step_energy_j_per_m2, 0.0, std::max(0.0, incident_step_energy_j_per_m2));

    for (auto& layer : status_.layers)
    {
        layer.dose_gy = layer.deposited_energy_j_per_m3 / std::max(density_kg_per_m3, kNumericalFloor);
    }

    status_.track_count = next_track_id_;
    if (status_.track_count > 0)
    {
        const double track_count = static_cast<double>(status_.track_count);
        status_.track_mean_terminal_energy_ev =
            track_terminal_energy_accumulator_ev_ / std::max(track_count, 1.0);
        status_.track_secondary_yield_per_primary =
            status_.track_secondary_event_count / std::max(track_count, 1.0);
    }
    return true;
}

std::size_t RadiationDoseAlgorithm::layerIndexForDepth(double depth_m) const
{
    if (status_.layers.empty())
    {
        return 0;
    }

    if (depth_m <= 0.0)
    {
        return 0;
    }

    const double capped_depth_m = std::min(depth_m, std::max(0.0, config_.thickness_m - 1.0e-15));
    const std::size_t index = static_cast<std::size_t>(
        capped_depth_m / std::max(layer_thickness_m_, kNumericalFloor));
    return std::min(index, status_.layers.size() - 1);
}

void RadiationDoseAlgorithm::appendTrackPoint(const RadiationStatus::TrackPoint& point)
{
    if (status_.track_points.size() < config_.monte_carlo_max_recorded_points)
    {
        status_.track_points.push_back(point);
    }
    else
    {
        status_.dropped_track_points += 1;
    }
}

bool RadiationDoseAlgorithm::exportTrackCsv(const std::filesystem::path& csv_path) const
{
    std::filesystem::path stem_path = csv_path;
    stem_path.replace_extension();
    const auto track_csv_path = stem_path.string() + ".tracks.csv";

    std::ofstream output(track_csv_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        return false;
    }

    output << "track_id,step_id,time_s,x_m,z_m,energy_ev,deposited_energy_ev,secondary_energy_ev,scattering_sigma_rad,step_length_m,particle_weight\n";
    for (const auto& point : status_.track_points)
    {
        output << point.track_id << ',' << point.step_id << ',' << point.time_s << ',' << point.x_m
               << ',' << point.z_m << ',' << point.energy_ev << ',' << point.deposited_energy_ev
               << ',' << point.secondary_energy_ev << ',' << point.scattering_sigma_rad
               << ',' << point.step_length_m << ',' << point.particle_weight << '\n';
    }

    return static_cast<bool>(output);
}

double RadiationDoseAlgorithm::depositionEfficiency() const
{
    const double physics_scale = physicsListStoppingScale();
    const double mass_stopping_power_mev_cm2_per_g =
        geant4AlignedMassStoppingPowerMeVcm2PerG(config_.mean_energy_ev) * physics_scale;
    const double rho_kg_per_m3 = effectiveMassDensityKgPerM3();
    const double areal_density_g_per_cm2 = rho_kg_per_m3 * config_.thickness_m / 10.0;
    const double mean_energy_mev = std::max(config_.mean_energy_ev * 1.0e-6, kNumericalFloor);
    const double energy_loss_mev =
        std::max(0.0, mass_stopping_power_mev_cm2_per_g * areal_density_g_per_cm2);
    const double geant4_aligned_efficiency = std::clamp(energy_loss_mev / mean_energy_mev, 0.0, 1.0);

    switch (config_.particle_species)
    {
    case ParticleSpecies::Electron:
        return std::clamp(std::max(geant4_aligned_efficiency, 0.05) * physics_scale, 0.0, 1.0);
    case ParticleSpecies::Proton:
        return std::clamp(std::max(geant4_aligned_efficiency, 0.08) * physics_scale, 0.0, 1.0);
    case ParticleSpecies::HeavyIon:
        return std::clamp(std::max(geant4_aligned_efficiency, 0.10) * physics_scale, 0.0, 1.0);
    }
    return std::clamp(std::max(geant4_aligned_efficiency, 0.05) * physics_scale, 0.0, 1.0);
}

double RadiationDoseAlgorithm::electronInelasticEnergyMultiplier(double energy_ev) const
{
    if (config_.particle_species != ParticleSpecies::Electron || !config_.electron_process.enable)
    {
        return 1.0;
    }

    const double normalized_energy =
        std::clamp(energy_ev / std::max(config_.mean_energy_ev, 1.0), 0.05, 4.0);
    const double low_energy_boost = 1.0 + 0.30 / std::sqrt(normalized_energy);
    return std::clamp(low_energy_boost * config_.electron_process.inelastic_energy_loss_scale,
                      0.5, 1.9);
}

double RadiationDoseAlgorithm::electronElasticScatteringSigmaRad(double energy_ev,
                                                                  double step_m) const
{
    if (config_.particle_species != ParticleSpecies::Electron || !config_.electron_process.enable)
    {
        return 0.0;
    }

    const double step_ratio =
        std::sqrt(std::max(step_m, 1.0e-9) / std::max(config_.thickness_m, kNumericalFloor));
    const double energy_term = std::sqrt(5.0e4 / std::max(energy_ev, 1.0));
    const double base_sigma = 0.012 * step_ratio * energy_term * physicsListScatteringScale();
    return std::clamp(base_sigma * config_.electron_process.elastic_scattering_scale,
                      0.0, 0.18);
}

double RadiationDoseAlgorithm::electronSecondaryDepthWeight(double depth_m) const
{
    const double attenuation_fraction =
        std::clamp(config_.electron_process.secondary_depth_attenuation_fraction, 0.05, 1.0);
    const double attenuation_length =
        std::max(layer_thickness_m_, config_.thickness_m * attenuation_fraction);
    return std::exp(-depth_m / attenuation_length);
}

double RadiationDoseAlgorithm::electronSecondaryEnergyEv(double deposited_energy_ev,
                                                         double depth_m) const
{
    if (config_.particle_species != ParticleSpecies::Electron || !config_.electron_process.enable)
    {
        return 0.0;
    }

    const double secondary_fraction =
        std::clamp(config_.electron_process.secondary_energy_fraction, 0.0, 0.45);
    return std::max(0.0,
                    deposited_energy_ev * secondary_fraction * electronSecondaryDepthWeight(depth_m));
}

void RadiationDoseAlgorithm::updateElectronProcessRunningMeans(double inelastic_multiplier,
                                                               double scattering_sigma_rad)
{
    electron_process_sample_count_ += 1.0;
    electron_process_inelastic_accumulator_ += inelastic_multiplier;
    electron_process_scattering_accumulator_ += scattering_sigma_rad;

    const double denominator = std::max(electron_process_sample_count_, 1.0);
    status_.electron_process_mean_inelastic_multiplier =
        electron_process_inelastic_accumulator_ / denominator;
    status_.electron_process_mean_scattering_sigma_rad =
        electron_process_scattering_accumulator_ / denominator;
}

double RadiationDoseAlgorithm::protonNuclearEnergyMultiplier(double energy_ev) const
{
    if (config_.particle_species != ParticleSpecies::Proton || !config_.proton_process.enable)
    {
        return 1.0;
    }

    const double normalized_energy =
        std::clamp(energy_ev / std::max(config_.mean_energy_ev, 1.0), 0.10, 4.0);
    const double high_energy_boost = 1.0 + 0.22 * std::sqrt(normalized_energy);
    return std::clamp(high_energy_boost * config_.proton_process.nuclear_reaction_energy_scale,
                      0.6, 2.0);
}

double RadiationDoseAlgorithm::protonElasticScatteringSigmaRad(double energy_ev,
                                                                double step_m) const
{
    if (config_.particle_species != ParticleSpecies::Proton || !config_.proton_process.enable)
    {
        return 0.0;
    }

    const double step_ratio =
        std::sqrt(std::max(step_m, 1.0e-9) / std::max(config_.thickness_m, kNumericalFloor));
    const double energy_term = std::sqrt(1.5e5 / std::max(energy_ev, 1.0));
    const double base_sigma = 0.006 * step_ratio * energy_term * physicsListScatteringScale();
    return std::clamp(base_sigma * config_.proton_process.elastic_scattering_scale,
                      0.0, 0.12);
}

double RadiationDoseAlgorithm::protonSecondaryDepthWeight(double depth_m) const
{
    const double buildup_fraction =
        std::clamp(config_.proton_process.secondary_depth_buildup_fraction, 0.05, 1.2);
    const double buildup_length =
        std::max(layer_thickness_m_, config_.thickness_m * buildup_fraction);
    const double remaining_depth_m = std::max(0.0, config_.thickness_m - depth_m);
    return std::exp(-remaining_depth_m / buildup_length);
}

double RadiationDoseAlgorithm::protonSecondaryEnergyEv(double deposited_energy_ev,
                                                       double depth_m) const
{
    if (config_.particle_species != ParticleSpecies::Proton || !config_.proton_process.enable)
    {
        return 0.0;
    }

    const double secondary_fraction =
        std::clamp(config_.proton_process.secondary_energy_fraction, 0.0, 0.40);
    return std::max(0.0,
                    deposited_energy_ev * secondary_fraction * protonSecondaryDepthWeight(depth_m));
}

void RadiationDoseAlgorithm::updateProtonProcessRunningMeans(double nuclear_multiplier,
                                                             double scattering_sigma_rad)
{
    proton_process_sample_count_ += 1.0;
    proton_process_nuclear_accumulator_ += nuclear_multiplier;
    proton_process_scattering_accumulator_ += scattering_sigma_rad;

    const double denominator = std::max(proton_process_sample_count_, 1.0);
    status_.proton_process_mean_nuclear_multiplier =
        proton_process_nuclear_accumulator_ / denominator;
    status_.proton_process_mean_scattering_sigma_rad =
        proton_process_scattering_accumulator_ / denominator;
}

double RadiationDoseAlgorithm::heavyIonEffectiveChargeMultiplier(double energy_ev) const
{
    if (config_.particle_species != ParticleSpecies::HeavyIon || !config_.heavy_ion_process.enable)
    {
        return 1.0;
    }

    const double normalized_energy =
        std::clamp(energy_ev / std::max(config_.mean_energy_ev, 1.0), 0.10, 5.0);
    const double charge_state_boost = 1.0 + 0.18 / std::sqrt(normalized_energy);
    return std::clamp(charge_state_boost * config_.heavy_ion_process.effective_charge_scale,
                      0.5, 2.4);
}

double RadiationDoseAlgorithm::heavyIonStoppingPowerMultiplier(double energy_ev) const
{
    if (config_.particle_species != ParticleSpecies::HeavyIon || !config_.heavy_ion_process.enable)
    {
        return 1.0;
    }

    const double normalized_energy =
        std::clamp(energy_ev / std::max(config_.mean_energy_ev, 1.0), 0.10, 6.0);
    const double stopping_boost = 1.0 + 0.26 * std::sqrt(normalized_energy);
    return std::clamp(stopping_boost * config_.heavy_ion_process.stopping_power_scale,
                      0.6, 3.0);
}

double RadiationDoseAlgorithm::heavyIonElasticScatteringSigmaRad(double energy_ev,
                                                                  double step_m) const
{
    if (config_.particle_species != ParticleSpecies::HeavyIon || !config_.heavy_ion_process.enable)
    {
        return 0.0;
    }

    const double step_ratio =
        std::sqrt(std::max(step_m, 1.0e-9) / std::max(config_.thickness_m, kNumericalFloor));
    const double energy_term = std::sqrt(3.0e5 / std::max(energy_ev, 1.0));
    const double base_sigma = 0.004 * step_ratio * energy_term * physicsListScatteringScale();
    return std::clamp(base_sigma * config_.heavy_ion_process.elastic_scattering_scale,
                      0.0, 0.10);
}

double RadiationDoseAlgorithm::heavyIonFragmentationDepthWeight(double depth_m) const
{
    const double buildup_fraction =
        std::clamp(config_.heavy_ion_process.fragmentation_depth_buildup_fraction, 0.05, 1.6);
    const double buildup_length =
        std::max(layer_thickness_m_, config_.thickness_m * buildup_fraction);
    const double remaining_depth_m = std::max(0.0, config_.thickness_m - depth_m);
    return std::exp(-remaining_depth_m / buildup_length);
}

double RadiationDoseAlgorithm::heavyIonFragmentationSecondaryEnergyEv(double deposited_energy_ev,
                                                                      double depth_m) const
{
    if (config_.particle_species != ParticleSpecies::HeavyIon || !config_.heavy_ion_process.enable)
    {
        return 0.0;
    }

    const double secondary_fraction =
        std::clamp(config_.heavy_ion_process.fragmentation_secondary_energy_fraction, 0.0, 0.55);
    return std::max(0.0,
                    deposited_energy_ev * secondary_fraction *
                        heavyIonFragmentationDepthWeight(depth_m));
}

void RadiationDoseAlgorithm::updateHeavyIonProcessRunningMeans(double effective_charge_multiplier,
                                                                double stopping_multiplier,
                                                                double scattering_sigma_rad)
{
    heavy_ion_process_sample_count_ += 1.0;
    heavy_ion_process_effective_charge_accumulator_ += effective_charge_multiplier;
    heavy_ion_process_stopping_accumulator_ += stopping_multiplier;
    heavy_ion_process_scattering_accumulator_ += scattering_sigma_rad;

    const double denominator = std::max(heavy_ion_process_sample_count_, 1.0);
    status_.heavy_ion_process_mean_effective_charge_multiplier =
        heavy_ion_process_effective_charge_accumulator_ / denominator;
    status_.heavy_ion_process_mean_stopping_multiplier =
        heavy_ion_process_stopping_accumulator_ / denominator;
    status_.heavy_ion_process_mean_scattering_sigma_rad =
        heavy_ion_process_scattering_accumulator_ / denominator;
}

double RadiationDoseAlgorithm::geant4AlignedMassStoppingPowerMeVcm2PerG(double mean_energy_ev) const
{
    const double z_over_a = effectiveMaterialZOverA();
    const double mean_excitation_ev = effectiveMeanExcitationEnergyEv();

    switch (config_.particle_species)
    {
    case ParticleSpecies::Electron:
        return geant4AlignedElectronMassStoppingPowerMeVcm2PerG(mean_energy_ev, z_over_a,
                                                                 mean_excitation_ev);
    case ParticleSpecies::Proton:
        return geant4AlignedIonMassStoppingPowerMeVcm2PerG(mean_energy_ev, z_over_a,
                                                            mean_excitation_ev, kProtonMassMeV,
                                                            1.0);
    case ParticleSpecies::HeavyIon:
        return geant4AlignedIonMassStoppingPowerMeVcm2PerG(mean_energy_ev, z_over_a,
                                                            mean_excitation_ev, kAlphaMassMeV,
                                                            2.0);
    }
    return 0.0;
}

double RadiationDoseAlgorithm::geant4AlignedElectronMassStoppingPowerMeVcm2PerG(
    double mean_energy_ev, double z_over_a, double mean_excitation_ev) const
{
    const double kinetic_mev = std::max(mean_energy_ev * 1.0e-6, 1.0e-9);
    const double gamma = 1.0 + kinetic_mev / kElectronMassMeV;
    const double beta2 = std::clamp(1.0 - 1.0 / (gamma * gamma), 1.0e-12, 1.0 - 1.0e-12);
    const double tau = kinetic_mev / kElectronMassMeV;
    const double excitation_mev = std::max(mean_excitation_ev * 1.0e-6, 1.0e-9);
    const double argument = 1.0 + tau * tau / std::max(excitation_mev * excitation_mev, 1.0e-18);
    const double log_term = std::max(0.0, std::log(argument));

    const double stopping = 0.1535 * z_over_a * (log_term - beta2) / beta2;
    return std::max(0.0, stopping);
}

double RadiationDoseAlgorithm::geant4AlignedIonMassStoppingPowerMeVcm2PerG(
    double mean_energy_ev, double z_over_a, double mean_excitation_ev, double particle_mass_mev,
    double projectile_charge_number) const
{
    const double kinetic_mev = std::max(mean_energy_ev * 1.0e-6, 1.0e-9);
    const double gamma = 1.0 + kinetic_mev / std::max(particle_mass_mev, kNumericalFloor);
    const double beta2 = std::clamp(1.0 - 1.0 / (gamma * gamma), 1.0e-12, 1.0 - 1.0e-12);
    const double beta = std::sqrt(beta2);
    const double gamma2 = gamma * gamma;

    const double z_pow = std::pow(std::max(projectile_charge_number, 1.0), -2.0 / 3.0);
    const double effective_charge = projectile_charge_number * (1.0 - std::exp(-125.0 * beta * z_pow));

    const double mass_ratio = kElectronMassMeV / std::max(particle_mass_mev, kNumericalFloor);
    const double wmax_mev = (2.0 * kElectronMassMeV * beta2 * gamma2) /
                            (1.0 + 2.0 * gamma * mass_ratio + mass_ratio * mass_ratio);
    const double excitation_mev = std::max(mean_excitation_ev * 1.0e-6, 1.0e-9);
    const double log_argument =
        std::max((2.0 * kElectronMassMeV * beta2 * gamma2 * std::max(wmax_mev, 1.0e-12)) /
                     std::max(excitation_mev * excitation_mev, 1.0e-24),
                 1.0 + 1.0e-12);
    const double stopping =
        kBetheKMeVCm2PerG * z_over_a * (effective_charge * effective_charge / beta2) *
        (0.5 * std::log(log_argument) - beta2);
    return std::max(0.0, stopping);
}

double RadiationDoseAlgorithm::effectiveMaterialZOverA() const
{
    if (stopping_data_resolution_.valid)
    {
        return stopping_data_resolution_.z_over_a;
    }

    if (material_ != nullptr)
    {
        const double configured = material_->getScalarProperty("z_over_a", 0.0);
        if (configured > 0.0)
        {
            return configured;
        }

        const auto name = normalizeToken(material_->getName());
        if (name.find("kapton") != std::string::npos)
        {
            return 0.52;
        }
        if (name.find("ptfe") != std::string::npos || name.find("teflon") != std::string::npos)
        {
            return 0.48;
        }
        if (name.find("aluminum") != std::string::npos)
        {
            return 0.481;
        }
    }
    return 0.50;
}

double RadiationDoseAlgorithm::effectiveMeanExcitationEnergyEv() const
{
    if (stopping_data_resolution_.valid)
    {
        return stopping_data_resolution_.mean_excitation_ev;
    }

    if (material_ != nullptr)
    {
        const double configured = material_->getScalarProperty("mean_excitation_ev", 0.0);
        if (configured > 0.0)
        {
            return configured;
        }

        const auto name = normalizeToken(material_->getName());
        if (name.find("kapton") != std::string::npos)
        {
            return 79.0;
        }
        if (name.find("ptfe") != std::string::npos || name.find("teflon") != std::string::npos)
        {
            return 99.0;
        }
        if (name.find("aluminum") != std::string::npos)
        {
            return 166.0;
        }
    }
    return 85.0;
}

double RadiationDoseAlgorithm::physicsListStoppingScale() const
{
    if (stopping_data_resolution_.valid)
    {
        return stopping_data_resolution_.stopping_scale;
    }

    return 1.00;
}

double RadiationDoseAlgorithm::physicsListScatteringScale() const
{
    if (stopping_data_resolution_.valid)
    {
        return stopping_data_resolution_.scattering_scale;
    }

    return 1.00;
}

double RadiationDoseAlgorithm::effectiveMassDensityKgPerM3() const
{
    if (material_ != nullptr)
    {
        const double rho = material_->getMassDensityKgPerM3();
        if (std::isfinite(rho) && rho > 0.0)
        {
            return rho;
        }
    }
    return 1200.0;
}

void RadiationDoseAlgorithm::updateConservationMetrics(double incident_step_energy_j_per_m2,
                                                       double deposited_step_energy_j_per_m2)
{
    status_.incident_energy_j_per_m2 += incident_step_energy_j_per_m2;
    status_.deposited_energy_j_per_m2 += deposited_step_energy_j_per_m2;
    status_.escaped_energy_j_per_m2 =
        std::max(0.0, status_.incident_energy_j_per_m2 - status_.deposited_energy_j_per_m2);

    const double denominator = std::max(std::abs(status_.incident_energy_j_per_m2), kNumericalFloor);
    const double closure = std::abs(status_.incident_energy_j_per_m2 -
                                    (status_.deposited_energy_j_per_m2 + status_.escaped_energy_j_per_m2));
    status_.energy_conservation_error = closure / denominator;
}

std::string toString(ParticleSpecies species)
{
    switch (species)
    {
    case ParticleSpecies::Electron:
        return "electron";
    case ParticleSpecies::Proton:
        return "proton";
    case ParticleSpecies::HeavyIon:
        return "heavy_ion";
    }
    return "unknown";
}

std::string toString(RadiationPhysicsList physics_list)
{
    switch (physics_list)
    {
    case RadiationPhysicsList::Geant4EmStandard:
        return "geant4_em_standard";
    case RadiationPhysicsList::Geant4EmLivermore:
        return "geant4_em_livermore";
    case RadiationPhysicsList::Geant4EmPenelope:
        return "geant4_em_penelope";
    case RadiationPhysicsList::Geant4SpaceShielding:
        return "geant4_space_shielding";
    }
    return "unknown";
}

} // namespace Radiation
} // namespace Toolkit
} // namespace SCDAT
