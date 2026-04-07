#include "PlasmaChannelModel.h"

#include <algorithm>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

ArcChannelState PlasmaChannelModel::update(double electric_field_v_per_m,
                                           double emission_current_density_a_per_m2,
                                           double gap_distance_m) const
{
    ArcChannelState state;
    const double conductivity = electric_field_v_per_m * parameters_.field_conductivity_gain +
                                emission_current_density_a_per_m2 *
                                    parameters_.emission_conductivity_gain;
    state.conductivity_s_per_m = std::clamp(
        conductivity, parameters_.min_conductivity_s_per_m, parameters_.max_conductivity_s_per_m);

    const double current_density =
        emission_current_density_a_per_m2 *
        (1.0 + electric_field_v_per_m * parameters_.field_current_gain);
    state.current_density_a_per_m2 =
        std::clamp(current_density, parameters_.min_current_density_a_per_m2,
                   parameters_.max_current_density_a_per_m2);

    const double resistance = gap_distance_m /
                              std::max(parameters_.min_conductivity_s_per_m,
                                       state.conductivity_s_per_m) *
                              parameters_.resistance_scale;
    state.resistance_ohm =
        std::clamp(resistance, parameters_.min_resistance_ohm, parameters_.max_resistance_ohm);
    return state;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
