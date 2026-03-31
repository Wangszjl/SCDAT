#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class ArcPICEmissionModel
{
  public:
    double computeEmissionCurrentDensity(double electric_field_v_per_m,
                                         double cathode_temperature_k) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
