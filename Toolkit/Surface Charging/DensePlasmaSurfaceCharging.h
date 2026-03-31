#pragma once

#include "ChargeAccumulationModel.h"

#include "../Plasma Analysis/FluidAlgorithmConfig.h"

#include "../../Tools/Material/include/MaterialDatabase.h"
#include "../../Tools/Material/include/SurfaceInteraction.h"
#include "../../Tools/Output/include/ResultExporter.h"

#include <filesystem>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

enum class SurfaceChargingRegime
{
    Generic,
    LeoFlowingPlasma,
    GeoKineticPicLike,
    ThrusterPlume
};

struct SurfaceCurrents
{
    double electron_current_a_per_m2 = 0.0;
    double ion_current_a_per_m2 = 0.0;
    double secondary_emission_a_per_m2 = 0.0;
    double photo_emission_a_per_m2 = 0.0;
    double thermionic_emission_a_per_m2 = 0.0;
    double field_emission_a_per_m2 = 0.0;
    double total_current_a_per_m2 = 0.0;
};

struct EmissionModelParameters
{
    double surface_temperature_k = 300.0;
    double enhancement_factor = 1.0;
    double photon_flux_m2_s = 0.0;
};

struct SurfaceChargingConfig
{
    double surface_area_m2 = 1.0e-4;
    bool floating = true;
    bool derive_capacitance_from_material = true;
    bool derive_sheath_length_from_plasma = true;
    double capacitance_per_area_f_per_m2 = 8.0e-10;
    double dielectric_thickness_m = 1.25e-4;
    double sheath_length_m = 1.0e-3;
    double minimum_sheath_length_m = 5.0e-5;
    double maximum_sheath_length_m = 10.0;
    double radiation_conductivity_coefficient = 0.0;
    double radiation_conductivity_exponent = 1.0;
    SurfaceChargingRegime regime = SurfaceChargingRegime::Generic;
    double bulk_flow_velocity_m_per_s = 0.0;
    double flow_alignment_cosine = 1.0;
    double electron_flow_coupling = 0.0;
    double electron_collection_coefficient = 1.0;
    double ion_collection_coefficient = 1.0;
    double ion_directed_velocity_m_per_s = 0.0;
    bool enable_pic_calibration = false;
    std::size_t pic_calibration_samples = 4096;
    double max_abs_potential_v = 5.0e4;
    double max_abs_current_density_a_per_m2 = 5.0e3;
    std::size_t internal_substeps = 16;
    PlasmaAnalysis::PlasmaParameters plasma;
    Material::MaterialProperty material{2, Mesh::MaterialType::DIELECTRIC, "kapton"};
    EmissionModelParameters emission;
};

struct SurfaceChargingStatus
{
    double time_s = 0.0;
    SurfaceCurrents currents;
    ChargeAccumulationState state;
    double equilibrium_error = 0.0;
    bool equilibrium_reached = false;
    std::size_t steps_completed = 0;
};

class DensePlasmaSurfaceCharging
{
  public:
    bool initialize(const SurfaceChargingConfig& config);
    bool advance(double dt);
    const SurfaceChargingStatus& getStatus() const { return status_; }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

    SurfaceCurrents computeSurfaceCurrents(double surface_potential_v) const;
    double computeFloatingPotential() const;
    double recommendTimeStep(double remaining_time_s, double minimum_dt_s, double maximum_dt_s) const;

  private:
    double computeCapacitancePerArea() const;
    double computeEffectiveSheathLength() const;
    double computeEffectiveConductivity(double surface_potential_v) const;
    double computeLeakageCurrentDensity(double surface_potential_v) const;
    double computeNetCurrentDensity(double surface_potential_v) const;
    double computeRadiationInducedConductivity() const;
    void applyGeoPicCalibration();
    double estimateElectronFluxPicLike(double surface_potential_v, std::size_t samples) const;
    double estimateIonFluxPicLike(double surface_potential_v, std::size_t samples) const;
    double estimateCurrentDerivative(double surface_potential_v) const;
    double advancePotentialImplicit(double surface_potential_v, double dt) const;

    SurfaceChargingConfig config_;
    SurfaceChargingStatus status_;
    ChargeAccumulationModel accumulation_model_;
    Material::SurfaceInteraction surface_interaction_;
    Output::ResultExporter exporter_;
    std::vector<double> history_time_;
    std::vector<double> history_potential_;
    std::vector<double> history_charge_;
    std::vector<double> history_current_;
    std::vector<double> history_electron_current_;
    std::vector<double> history_ion_current_;
    std::vector<double> history_secondary_current_;
    std::vector<double> history_photo_current_;
    std::vector<double> history_thermionic_current_;
    std::vector<double> history_field_emission_current_;
    std::vector<double> history_leakage_current_;
    std::vector<double> history_net_current_;
    std::vector<double> history_capacitance_;
    std::vector<double> history_effective_conductivity_;
    std::vector<double> history_effective_sheath_length_;
    std::vector<double> history_adaptive_time_step_;
    std::vector<double> history_electron_calibration_factor_;
    std::vector<double> history_ion_calibration_factor_;
    std::vector<double> history_equilibrium_potential_;
    std::vector<double> history_equilibrium_error_;
    double electron_pic_calibration_factor_ = 1.0;
    double ion_pic_calibration_factor_ = 1.0;
    bool initialized_ = false;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
