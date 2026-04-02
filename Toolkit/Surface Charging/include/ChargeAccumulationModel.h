#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct ChargeAccumulationState
{
    double surface_charge_density_c_per_m2 = 0.0;
    double surface_potential_v = 0.0;
    double capacitance_per_area_f_per_m2 = 1.0e-9;
    double stored_energy_j_per_m2 = 0.0;
};

class ChargeAccumulationModel
{
  public:
    ChargeAccumulationState advance(const ChargeAccumulationState& state,
                                    double total_current_density_a_per_m2, double dt) const;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
