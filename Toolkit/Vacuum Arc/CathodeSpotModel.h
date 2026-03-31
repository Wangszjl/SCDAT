#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class CathodeSpotModel
{
  public:
    double emitCurrentDensity(double electric_field_v_per_m, double surface_temperature_k) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
