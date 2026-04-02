#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

struct ArcChannelState
{
    double conductivity_s_per_m = 0.0;
    double current_density_a_per_m2 = 0.0;
    double resistance_ohm = 0.0;
};

class PlasmaChannelModel
{
  public:
    ArcChannelState update(double electric_field_v_per_m, double emission_current_density_a_per_m2,
                           double gap_distance_m) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
