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
    state.conductivity_s_per_m =
        std::max(1.0, electric_field_v_per_m * 1.0e-5 + emission_current_density_a_per_m2 * 10.0);
    state.current_density_a_per_m2 =
        emission_current_density_a_per_m2 * (1.0 + electric_field_v_per_m * 1.0e-8);
    state.resistance_ohm =
        gap_distance_m / std::max(1.0, state.conductivity_s_per_m) * 1.0e3;
    return state;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
