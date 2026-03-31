#pragma once

#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

struct CylindricalProfile
{
    std::vector<double> radius_m;
    std::vector<double> potential_v;
};

class CylindricalPICSolver
{
  public:
    CylindricalProfile solve(double channel_radius_m, double current_density_a_per_m2,
                             std::size_t points) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
