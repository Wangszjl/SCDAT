#include "SurfaceDischargeArcAlgorithm.h"

#include <algorithm>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

bool SurfaceDischargeArcAlgorithm::initialize(const DischargeConfiguration& config)
{
    config_ = config;
    status_ = DischargeStatus{};
    status_.cathode_temperature_k = config.cathode_temperature_k;
    status_.anode_temperature_k = config.anode_temperature_k;
    initialized_ = true;
    history_time_.clear();
    history_current_.clear();
    history_current_density_.clear();
    history_cathode_temperature_.clear();
    history_anode_temperature_.clear();
    history_conductivity_.clear();
    return true;
}

bool SurfaceDischargeArcAlgorithm::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }

    const double effective_field = coupling_.computeEffectiveField(
        config_.applied_field_v_per_m, config_.surface_potential_v,
        config_.surface_charge_density_c_per_m2);
    const double avalanche_gain = avalanche_model_.computeGain(effective_field, config_.gap_distance_m);

    if (!status_.discharge_active)
    {
        status_.discharge_active = detector_.shouldTrigger(effective_field, avalanche_gain);
        status_.active_arc_channels = status_.discharge_active ? 1u : 0u;
    }

    if (status_.discharge_active)
    {
        const double cathode_emission =
            cathode_spot_model_.emitCurrentDensity(effective_field, status_.cathode_temperature_k);
        const double emitted_current_density =
            emission_model_.computeEmissionCurrentDensity(effective_field, status_.cathode_temperature_k) +
            cathode_emission;
        channel_state_ =
            plasma_channel_model_.update(effective_field, emitted_current_density, config_.gap_distance_m);
        channel_state_ = integrator_.advance(channel_state_, dt);
        radial_profile_ =
            cylindrical_solver_.solve(config_.channel_radius_m, channel_state_.current_density_a_per_m2, 24);

        status_.peak_current_density_a_per_m2 =
            std::max(status_.peak_current_density_a_per_m2, channel_state_.current_density_a_per_m2);
        status_.total_discharge_current_a =
            channel_state_.current_density_a_per_m2 * 3.14159265358979323846 * config_.channel_radius_m *
            config_.channel_radius_m;
        status_.cathode_temperature_k += emitted_current_density * dt * 2.0e2;
        status_.anode_temperature_k = anode_heating_model_.advanceTemperature(
            status_.anode_temperature_k, channel_state_.current_density_a_per_m2, dt);
        status_.channel_conductivity_s_per_m = channel_state_.conductivity_s_per_m;
    }

    status_.current_time_s += dt;
    history_time_.push_back(status_.current_time_s);
    history_current_.push_back(status_.total_discharge_current_a);
    history_current_density_.push_back(channel_state_.current_density_a_per_m2);
    history_cathode_temperature_.push_back(status_.cathode_temperature_k);
    history_anode_temperature_.push_back(status_.anode_temperature_k);
    history_conductivity_.push_back(status_.channel_conductivity_s_per_m);
    return true;
}

void SurfaceDischargeArcAlgorithm::reset()
{
    *this = SurfaceDischargeArcAlgorithm{};
}

bool SurfaceDischargeArcAlgorithm::exportResults(const std::filesystem::path& csv_path) const
{
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "time_s";
    data_set.axis_values = history_time_;
    data_set.scalar_series["discharge_current_a"] = history_current_;
    data_set.scalar_series["current_density_a_per_m2"] = history_current_density_;
    data_set.scalar_series["cathode_temperature_k"] = history_cathode_temperature_;
    data_set.scalar_series["anode_temperature_k"] = history_anode_temperature_;
    data_set.scalar_series["channel_conductivity_s_per_m"] = history_conductivity_;
    data_set.metadata["module"] = "Vacuum Arc";
    return static_cast<bool>(exporter_.exportDataSet(csv_path, data_set));
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
