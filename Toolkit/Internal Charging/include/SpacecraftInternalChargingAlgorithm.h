#pragma once

#include "DielectricBreakdownModel.h"
#include "DielectricMaterialDatabase.h"
#include "DischargeChannelModel.h"
#include "ParticleTransportModel.h"

#include "../../Tools/Output/include/ResultExporter.h"

#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

struct InternalChargingConfiguration
{
    std::size_t layers = 16;
    double thickness_m = 2.0e-3;
    double area_m2 = 1.0e-2;
    double incident_current_density_a_per_m2 = 2.0e-8;
    double incident_energy_ev = 5.0e3;
    std::string material_name = "kapton";
};

struct InternalChargingStatus
{
    double time_s = 0.0;
    double total_stored_energy_j = 0.0;
    double discharge_probability = 0.0;
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
    const InternalChargingStatus& getStatus() const { return status_; }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

  private:
    void updateFieldEstimate();

    InternalChargingConfiguration config_;
    InternalChargingStatus status_;
    ParticleTransportModel transport_model_;
    DielectricMaterialDatabase material_database_;
    DielectricBreakdownModel breakdown_model_;
    DischargeChannelModel discharge_model_;
    const Material::MaterialProperty* material_ = nullptr;
    Output::ResultExporter exporter_;
    bool initialized_ = false;
};

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
