#include "CylindricalPICSolver.h"

#include "../../Tools/Solver/include/SuperLUSolver.h"

#include <algorithm>
#include <cmath>

namespace
{

std::vector<double> buildLinearFallbackPotential(double channel_radius_m,
                                                 double current_density_a_per_m2,
                                                 std::size_t points)
{
    std::vector<double> potential(points, 0.0);
    for (std::size_t i = 0; i < points; ++i)
    {
        const double radius =
            channel_radius_m * static_cast<double>(i) / static_cast<double>(points - 1);
        potential[i] = current_density_a_per_m2 * 1.0e-3 *
                       (1.0 - radius / std::max(1.0e-12, channel_radius_m));
    }
    return potential;
}

} // namespace

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
    const double radius = std::max(1.0e-12, channel_radius_m);
    const double dr = radius / static_cast<double>(points - 1);

    Solver::CSRMatrix matrix(static_cast<int>(points), static_cast<int>(points),
                             static_cast<int>(3 * points));
    std::vector<double> rhs(points, 0.0);
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_pointers(points + 1, 0);

    values.reserve(3 * points);
    col_indices.reserve(3 * points);

    auto append = [&](std::size_t row, std::size_t col, double value)
    {
        values.push_back(value);
        col_indices.push_back(static_cast<int>(col));
        row_pointers[row + 1] = static_cast<int>(values.size());
    };

    // Axis boundary: dphi/dr = 0 at r = 0.
    append(0, 0, 1.0);
    append(0, 1, -1.0);
    rhs[0] = 0.0;

    const double source_term =
        4.0 * current_density_a_per_m2 * 1.0e-3 / std::max(1.0e-18, radius * radius);

    for (std::size_t i = 1; i + 1 < points; ++i)
    {
        const double ri = std::max(dr * static_cast<double>(i), 1.0e-12);
        const double inv_dr2 = 1.0 / std::max(1.0e-24, dr * dr);
        const double coeff_left = inv_dr2 - 1.0 / (2.0 * ri * dr);
        const double coeff_center = -2.0 * inv_dr2;
        const double coeff_right = inv_dr2 + 1.0 / (2.0 * ri * dr);

        append(i, i - 1, coeff_left);
        append(i, i, coeff_center);
        append(i, i + 1, coeff_right);
        rhs[i] = -source_term;
    }

    // Edge boundary: phi(R) = 0.
    append(points - 1, points - 1, 1.0);
    rhs[points - 1] = 0.0;

    matrix.setData(values, col_indices, row_pointers);

    std::vector<double> solved_potential;
    Solver::SuperLUSolver solver;
    bool solved_ok = matrix.validate() && solver.factorize(matrix) && solver.solve(rhs, solved_potential) &&
                     solved_potential.size() == points;

    if (!solved_ok)
    {
        solved_potential = buildLinearFallbackPotential(channel_radius_m, current_density_a_per_m2, points);
    }

    for (std::size_t i = 0; i < points; ++i)
    {
        const double ri = radius * static_cast<double>(i) / static_cast<double>(points - 1);
        profile.radius_m.push_back(ri);
        const double phi = solved_potential[i];
        profile.potential_v.push_back(std::isfinite(phi) ? phi : 0.0);
    }

    return profile;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
