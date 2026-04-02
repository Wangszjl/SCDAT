#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class SurfaceArcCoupling
{
  public:
    double computeEffectiveField(double applied_field_v_per_m, double surface_potential_v,
                                 double surface_charge_density_c_per_m2) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
