#pragma once

#include "DielectricBreakdownModel.h"
#include "DielectricMaterialDatabase.h"
#include "DischargeChannelModel.h"
#include "ParticleTransportModel.h"

#include "../../Tools/Output/include/ResultExporter.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

enum class InternalChargingSourceMode
{
  Preset,
  Radiation
};

struct InternalChargingRadiationDrive
{
  double incident_current_density_a_per_m2 = 0.0;
  double incident_energy_ev = 0.0;
  double incident_charge_state_abs = 1.0;
};

struct InternalChargingConfiguration
{
    std::size_t layers = 16;
    double thickness_m = 2.0e-3;
    double area_m2 = 1.0e-2;
    double incident_current_density_a_per_m2 = 2.0e-8;
    double incident_energy_ev = 5.0e3;
  InternalChargingSourceMode source_mode = InternalChargingSourceMode::Preset;
    std::string material_name = "kapton";
    double incident_charge_state_abs = 1.0;

  // Radiation-induced conductivity model coefficients.
  double radiation_conductivity_charge_gain_s_per_m_per_c_per_m3 = 1.0e-6;
  double radiation_conductivity_dose_gain_s_per_m_per_gy = 1.0e-9;
  double radiation_conductivity_dose_rate_gain_s2_per_m_per_gy = 1.0e-12;
  double min_effective_conductivity_s_per_m = 1.0e-15;
  double max_effective_conductivity_s_per_m = 1.0e-6;
};

struct InternalChargingStatus
{
    double time_s = 0.0;
    double total_stored_energy_j = 0.0;
  double max_electric_field_v_per_m = 0.0;
    double discharge_probability = 0.0;
  double average_dose_gy = 0.0;
  double effective_conductivity_s_per_m = 0.0;
  double deposited_charge_rate_c_per_m3_s = 0.0;
  double deposited_particle_rate_per_m3_s = 0.0;
  double effective_incident_charge_state_abs = 1.0;
    std::vector<double> volume_charge_density_c_per_m3;
    std::vector<double> electric_field_v_per_m;
    std::size_t discharge_events = 0;
    bool breakdown_active = false;
};

class SpacecraftInternalChargingAlgorithm
{
  public:
    bool initialize(const InternalChargingConfiguration& config);
    bool advance(double dt);
    void setRadiationDrive(const InternalChargingRadiationDrive& drive) { radiation_drive_ = drive; }
    void clearRadiationDrive() { radiation_drive_.reset(); }
    const InternalChargingStatus& getStatus() const { return status_; }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

  private:
    void updateRadiationInducedConductivity(double dt,
                                            const std::vector<double>& previous_charge_density_c_per_m3,
                                            const std::vector<double>& previous_deposited_energy_j_per_m3,
                                            double incident_charge_state_abs);
    void updateFieldEstimate();

    InternalChargingConfiguration config_;
    InternalChargingStatus status_;
    ParticleTransportModel transport_model_;
    DielectricMaterialDatabase material_database_;
    DielectricBreakdownModel breakdown_model_;
    DischargeChannelModel discharge_model_;
    const Material::MaterialProperty* material_ = nullptr;
    std::optional<InternalChargingRadiationDrive> radiation_drive_;
    Output::ResultExporter exporter_;
    bool initialized_ = false;
};

  std::string toString(InternalChargingSourceMode source_mode);

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
