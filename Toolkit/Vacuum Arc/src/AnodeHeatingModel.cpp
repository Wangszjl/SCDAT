#include "AnodeHeatingModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

double AnodeHeatingModel::advanceTemperature(double current_temperature_k,
                                             double current_density_a_per_m2, double dt) const
{
    return current_temperature_k + current_density_a_per_m2 * dt * 5.0e2;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
