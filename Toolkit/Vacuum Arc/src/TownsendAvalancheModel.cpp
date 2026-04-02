#include "TownsendAvalancheModel.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

double TownsendAvalancheModel::computeGain(double electric_field_v_per_m, double gap_distance_m) const
{
    const double exponent = std::clamp((std::max(0.0, electric_field_v_per_m - 2.0e7) * 5.0e-6) *
                                           std::max(0.0, gap_distance_m) * 1.0e3,
                                       0.0, 40.0);
    return std::exp(exponent);
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
