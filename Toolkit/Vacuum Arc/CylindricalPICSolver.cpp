#include "CylindricalPICSolver.h"

#include <algorithm>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

CylindricalProfile CylindricalPICSolver::solve(double channel_radius_m, double current_density_a_per_m2,
                                               std::size_t points) const
{
    CylindricalProfile profile;
    points = std::max<std::size_t>(points, 4);
    for (std::size_t i = 0; i < points; ++i)
    {
        const double radius =
            channel_radius_m * static_cast<double>(i) / static_cast<double>(points - 1);
        profile.radius_m.push_back(radius);
        profile.potential_v.push_back(current_density_a_per_m2 * 1.0e-3 *
                                      (1.0 - radius / std::max(1.0e-12, channel_radius_m)));
    }
    return profile;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
