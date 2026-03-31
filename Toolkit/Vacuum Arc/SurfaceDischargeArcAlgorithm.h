#pragma once

#include "AnodeHeatingModel.h"
#include "ArcDischargeDetector.h"
#include "ArcPICEmissionModel.h"
#include "ArcPICIntegrator.h"
#include "CathodeSpotModel.h"
#include "CylindricalPICSolver.h"
#include "PlasmaChannelModel.h"
#include "SurfaceArcCoupling.h"
#include "TownsendAvalancheModel.h"

#include "../../Tools/Output/include/ResultExporter.h"

#include <filesystem>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

struct DischargeConfiguration
{
    double gap_distance_m = 5.0e-4;
    double applied_field_v_per_m = 2.5e7;
    double surface_potential_v = 0.0;
    double surface_charge_density_c_per_m2 = 0.0;
    double cathode_temperature_k = 450.0;
    double anode_temperature_k = 320.0;
    double channel_radius_m = 2.0e-5;
};

struct DischargeStatus
{
    bool discharge_active = false;
    double current_time_s = 0.0;
    double total_discharge_current_a = 0.0;
    double peak_current_density_a_per_m2 = 0.0;
    std::size_t active_arc_channels = 0;
    double cathode_temperature_k = 0.0;
    double anode_temperature_k = 0.0;
    double channel_conductivity_s_per_m = 0.0;
};

class SurfaceDischargeArcAlgorithm
{
  public:
    bool initialize(const DischargeConfiguration& config);
    bool advance(double dt);
    const DischargeStatus& getStatus() const { return status_; }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

  private:
    DischargeConfiguration config_;
    DischargeStatus status_;
    TownsendAvalancheModel avalanche_model_;
    CathodeSpotModel cathode_spot_model_;
    AnodeHeatingModel anode_heating_model_;
    PlasmaChannelModel plasma_channel_model_;
    ArcDischargeDetector detector_;
    ArcPICEmissionModel emission_model_;
    ArcPICIntegrator integrator_;
    CylindricalPICSolver cylindrical_solver_;
    SurfaceArcCoupling coupling_;
    Output::ResultExporter exporter_;
    ArcChannelState channel_state_;
    CylindricalProfile radial_profile_;
    std::vector<double> history_time_;
    std::vector<double> history_current_;
    std::vector<double> history_current_density_;
    std::vector<double> history_cathode_temperature_;
    std::vector<double> history_anode_temperature_;
    std::vector<double> history_conductivity_;
    bool initialized_ = false;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
