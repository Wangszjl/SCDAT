#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class TownsendAvalancheModel
{
  public:
    double computeGain(double electric_field_v_per_m, double gap_distance_m) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
