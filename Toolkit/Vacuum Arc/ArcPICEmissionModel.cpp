#include "ArcPICEmissionModel.h"

#include <algorithm>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

double ArcPICEmissionModel::computeEmissionCurrentDensity(double electric_field_v_per_m,
                                                          double cathode_temperature_k) const
{
    return std::max(0.0, electric_field_v_per_m - 1.5e7) * 1.0e-6 +
           std::max(0.0, cathode_temperature_k - 500.0) * 1.0e-2;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
