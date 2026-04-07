#pragma once

#include "RadiationStoppingDataGovernance.h"

#include "../../Tools/Material/include/MaterialDatabase.h"
#include "../../Tools/Output/include/ResultExporter.h"

#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace Radiation
{

enum class ParticleSpecies
{
    Electron,
    Proton,
    HeavyIon
};

enum class RadiationPhysicsList
{
    Geant4EmStandard,
    Geant4EmLivermore,
    Geant4EmPenelope,
    Geant4SpaceShielding
};

struct ElectronProcessEnhancementConfig
{
    bool enable = true;
    double elastic_scattering_scale = 1.0;
    double inelastic_energy_loss_scale = 1.0;
    double secondary_energy_fraction = 0.08;
    double secondary_depth_attenuation_fraction = 0.30;
};

struct ProtonProcessEnhancementConfig
{
    bool enable = true;
    double nuclear_reaction_energy_scale = 1.0;
    double elastic_scattering_scale = 1.0;
    double secondary_energy_fraction = 0.05;
    double secondary_depth_buildup_fraction = 0.45;
};

struct HeavyIonProcessEnhancementConfig
{
    bool enable = true;
    double effective_charge_scale = 1.0;
    double stopping_power_scale = 1.0;
    double elastic_scattering_scale = 1.0;
    double fragmentation_secondary_energy_fraction = 0.08;
    double fragmentation_depth_buildup_fraction = 0.60;
    double fragmentation_event_rate_scale = 1.0;
};

struct RadiationConfiguration
{
    std::size_t layers = 16;
    double thickness_m = 2.0e-3;
    double area_m2 = 1.0e-2;
    double particle_flux_m2_s = 5.0e10;
    double mean_energy_ev = 1.0e5;
    ParticleSpecies particle_species = ParticleSpecies::Proton;
    RadiationPhysicsList physics_list = RadiationPhysicsList::Geant4EmStandard;
    std::string material_name = "kapton";
    double attenuation_length_fraction = 0.25;

    // RD-003 governance controls.
    std::string stopping_data_version = "2026.04-rd003";
    bool stopping_data_allow_rollback = true;

    // Optional Geant4-style Monte Carlo transport. Disabled by default to preserve
    // historical deterministic behavior and existing regression baselines.
    bool enable_monte_carlo_transport = false;
    std::size_t monte_carlo_histories_per_step = 512;
    std::size_t monte_carlo_max_steps_per_track = 256;
    std::size_t monte_carlo_max_recorded_points = 200000;
    double monte_carlo_lateral_span_m = 2.0e-2;
    double monte_carlo_angular_sigma_deg = 3.0;
    unsigned int monte_carlo_seed = 20260403u;

    // RD-004 electron process enhancement controls.
    ElectronProcessEnhancementConfig electron_process;

    // RD-005 proton process enhancement controls.
    ProtonProcessEnhancementConfig proton_process;

    // RD-006 heavy-ion process enhancement controls.
    HeavyIonProcessEnhancementConfig heavy_ion_process;
};

struct RadiationLayerState
{
    double depth_m = 0.0;
    double deposited_energy_j_per_m3 = 0.0;
    double dose_gy = 0.0;
};

struct RadiationStatus
{
    double time_s = 0.0;
    double incident_energy_j_per_m2 = 0.0;
    double deposited_energy_j_per_m2 = 0.0;
    double escaped_energy_j_per_m2 = 0.0;
    double energy_conservation_error = 0.0;
    std::vector<RadiationLayerState> layers;

    struct TrackPoint
    {
        std::size_t track_id = 0;
        std::size_t step_id = 0;
        double time_s = 0.0;
        double x_m = 0.0;
        double z_m = 0.0;
        double energy_ev = 0.0;
        double deposited_energy_ev = 0.0;
        double secondary_energy_ev = 0.0;
        double scattering_sigma_rad = 0.0;
        double step_length_m = 0.0;
        double particle_weight = 0.0;
    };

    std::vector<TrackPoint> track_points;
    std::size_t track_count = 0;
    std::size_t dropped_track_points = 0;
    std::size_t track_step_count = 0;
    double track_total_path_length_m = 0.0;
    double track_mean_step_length_m = 0.0;
    double track_mean_scattering_sigma_rad = 0.0;
    double track_secondary_event_count = 0.0;
    double track_secondary_energy_j_per_m2 = 0.0;
    double track_secondary_yield_per_primary = 0.0;
    double track_mean_terminal_energy_ev = 0.0;
    double track_max_depth_m = 0.0;

    // RD-004 diagnostics.
    double electron_process_secondary_energy_j_per_m2 = 0.0;
    double electron_process_mean_inelastic_multiplier = 1.0;
    double electron_process_mean_scattering_sigma_rad = 0.0;
    double electron_process_secondary_event_count = 0.0;

    // RD-005 diagnostics.
    double proton_process_secondary_energy_j_per_m2 = 0.0;
    double proton_process_mean_nuclear_multiplier = 1.0;
    double proton_process_mean_scattering_sigma_rad = 0.0;
    double proton_process_secondary_event_count = 0.0;

    // RD-006 diagnostics.
    double heavy_ion_process_secondary_energy_j_per_m2 = 0.0;
    double heavy_ion_process_mean_effective_charge_multiplier = 1.0;
    double heavy_ion_process_mean_stopping_multiplier = 1.0;
    double heavy_ion_process_mean_scattering_sigma_rad = 0.0;
    double heavy_ion_process_fragmentation_event_count = 0.0;
};

class RadiationDoseAlgorithm
{
  public:
    bool initialize(const RadiationConfiguration& config);
    bool advance(double dt);
    const RadiationStatus& getStatus() const { return status_; }
        const RadiationStoppingResolution& getStoppingDataResolution() const
        {
                return stopping_data_resolution_;
        }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

  private:
    double depositionEfficiency() const;
        bool advanceMonteCarlo(double dt, double incident_step_energy_j_per_m2,
                                                     double& deposited_step_energy_j_per_m2);
        double electronInelasticEnergyMultiplier(double energy_ev) const;
        double electronElasticScatteringSigmaRad(double energy_ev, double step_m) const;
        double electronSecondaryEnergyEv(double deposited_energy_ev, double depth_m) const;
        double electronSecondaryDepthWeight(double depth_m) const;
        void updateElectronProcessRunningMeans(double inelastic_multiplier,
                                                                                     double scattering_sigma_rad);
        double protonNuclearEnergyMultiplier(double energy_ev) const;
        double protonElasticScatteringSigmaRad(double energy_ev, double step_m) const;
        double protonSecondaryEnergyEv(double deposited_energy_ev, double depth_m) const;
        double protonSecondaryDepthWeight(double depth_m) const;
        void updateProtonProcessRunningMeans(double nuclear_multiplier,
                                                                                 double scattering_sigma_rad);
        double heavyIonEffectiveChargeMultiplier(double energy_ev) const;
        double heavyIonStoppingPowerMultiplier(double energy_ev) const;
        double heavyIonElasticScatteringSigmaRad(double energy_ev, double step_m) const;
        double heavyIonFragmentationSecondaryEnergyEv(double deposited_energy_ev, double depth_m) const;
        double heavyIonFragmentationDepthWeight(double depth_m) const;
        void updateHeavyIonProcessRunningMeans(double effective_charge_multiplier,
                                                                                     double stopping_multiplier,
                                                                                     double scattering_sigma_rad);
        std::size_t layerIndexForDepth(double depth_m) const;
        void appendTrackPoint(const RadiationStatus::TrackPoint& point);
        bool exportTrackCsv(const std::filesystem::path& csv_path) const;
        double geant4AlignedMassStoppingPowerMeVcm2PerG(double mean_energy_ev) const;
        double geant4AlignedElectronMassStoppingPowerMeVcm2PerG(double mean_energy_ev,
                                                                                                                        double z_over_a,
                                                                                                                        double mean_excitation_ev) const;
        double geant4AlignedIonMassStoppingPowerMeVcm2PerG(double mean_energy_ev,
                                                                                                             double z_over_a,
                                                                                                             double mean_excitation_ev,
                                                                                                             double particle_mass_mev,
                                                                                                             double projectile_charge_number) const;
        double effectiveMaterialZOverA() const;
        double effectiveMeanExcitationEnergyEv() const;
        double physicsListStoppingScale() const;
        double physicsListScatteringScale() const;
    double effectiveMassDensityKgPerM3() const;
    void updateConservationMetrics(double incident_step_energy_j_per_m2,
                                   double deposited_step_energy_j_per_m2);

    RadiationConfiguration config_;
    RadiationStatus status_;
    Material::UnifiedMaterialDatabase material_database_;
    RadiationStoppingDataGovernance stopping_data_governance_;
    RadiationStoppingResolution stopping_data_resolution_;
    const Material::MaterialProperty* material_ = nullptr;
    std::string resolved_material_name_;
    Output::ResultExporter exporter_;
    bool initialized_ = false;
    double layer_thickness_m_ = 0.0;
    std::mt19937 random_engine_;
    std::size_t next_track_id_ = 0;
    double electron_process_inelastic_accumulator_ = 0.0;
    double electron_process_scattering_accumulator_ = 0.0;
    double electron_process_sample_count_ = 0.0;
    double proton_process_nuclear_accumulator_ = 0.0;
    double proton_process_scattering_accumulator_ = 0.0;
    double proton_process_sample_count_ = 0.0;
    double heavy_ion_process_effective_charge_accumulator_ = 0.0;
    double heavy_ion_process_stopping_accumulator_ = 0.0;
    double heavy_ion_process_scattering_accumulator_ = 0.0;
    double heavy_ion_process_sample_count_ = 0.0;
    double track_scattering_sigma_accumulator_ = 0.0;
    double track_terminal_energy_accumulator_ev_ = 0.0;
};

std::string toString(ParticleSpecies species);
std::string toString(RadiationPhysicsList physics_list);

} // namespace Radiation
} // namespace Toolkit
} // namespace SCDAT
