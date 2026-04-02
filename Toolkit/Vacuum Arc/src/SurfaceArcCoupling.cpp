#include "SurfaceArcCoupling.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

double SurfaceArcCoupling::computeEffectiveField(double applied_field_v_per_m,
                                                 double surface_potential_v,
                                                 double surface_charge_density_c_per_m2) const
{
    return applied_field_v_per_m + surface_potential_v * 5.0e4 +
           surface_charge_density_c_per_m2 * 1.0e9;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
