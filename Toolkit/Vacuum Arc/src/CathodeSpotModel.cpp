#include "CathodeSpotModel.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

double CathodeSpotModel::emitCurrentDensity(double electric_field_v_per_m,
                                            double surface_temperature_k) const
{
    const double field_term = std::max(0.0, electric_field_v_per_m - 1.0e7) * 5.0e-7;
    const double thermal_term = std::max(0.0, surface_temperature_k - 300.0) * 2.0e-3;
    return field_term + thermal_term;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
