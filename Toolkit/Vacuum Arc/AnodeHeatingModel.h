#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class AnodeHeatingModel
{
  public:
    double advanceTemperature(double current_temperature_k, double current_density_a_per_m2,
                              double dt) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
