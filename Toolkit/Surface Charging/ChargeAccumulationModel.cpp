#include "ChargeAccumulationModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

ChargeAccumulationState ChargeAccumulationModel::advance(const ChargeAccumulationState& state,
                                                         double total_current_density_a_per_m2,
                                                         double dt) const
{
    ChargeAccumulationState updated = state;
    updated.surface_charge_density_c_per_m2 += total_current_density_a_per_m2 * dt;
    updated.surface_potential_v =
        updated.surface_charge_density_c_per_m2 /
        (updated.capacitance_per_area_f_per_m2 > 0.0 ? updated.capacitance_per_area_f_per_m2 : 1.0);
    updated.stored_energy_j_per_m2 =
        0.5 * updated.capacitance_per_area_f_per_m2 * updated.surface_potential_v *
        updated.surface_potential_v;
    return updated;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
