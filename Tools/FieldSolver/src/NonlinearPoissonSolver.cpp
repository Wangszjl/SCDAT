#include "../include/NonlinearPoissonSolver.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace FieldSolver
{
namespace
{
constexpr double kEpsilon0 = 8.8541878128e-12;
}

bool NonlinearPoissonSolver::initializeGrid(const Geometry::Vector3D& grid_size,
                                            const Geometry::Vector3D& resolution)
{
    nx_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(resolution.x())));
    ny_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(resolution.y())));
    nz_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(resolution.z())));

    grid_size_ = grid_size;
    spacing_ = Geometry::Vector3D(
        grid_size.x() / static_cast<double>(std::max<std::size_t>(1, nx_ - 1)),
        grid_size.y() / static_cast<double>(std::max<std::size_t>(1, ny_ - 1)),
        grid_size.z() / static_cast<double>(std::max<std::size_t>(1, nz_ - 1)));

    potential_.assign(nx_ * ny_ * nz_, 0.0);
    electric_field_.assign(nx_ * ny_ * nz_, Geometry::Vector3D{});
    initialized_ = true;
    return true;
}

void NonlinearPoissonSolver::setChargeDensityFunction(ChargeDensityFunction charge_density)
{
    charge_density_ = std::move(charge_density);
}

void NonlinearPoissonSolver::setPermittivityFunction(PermittivityFunction permittivity)
{
    permittivity_ = std::move(permittivity);
}

void NonlinearPoissonSolver::setBoundaryPotential(BoundaryFunction boundary)
{
    boundary_ = std::move(boundary);
}

bool NonlinearPoissonSolver::solve()
{
    if (!initialized_)
    {
        return false;
    }

    const double dx2 = spacing_.x() * spacing_.x();
    const double dy2 = spacing_.y() * spacing_.y();
    const double dz2 = spacing_.z() * spacing_.z();
    const double denom = 2.0 / dx2 + 2.0 / dy2 + 2.0 / dz2;

    iterations_used_ = 0;
    last_max_delta_ = 0.0;
    last_residual_norm_ = 0.0;
    converged_ = false;

    const double derivative_step =
        std::max(1.0e-9, std::abs(parameters_.derivative_step_v));
    const double effective_diag_floor =
        std::max(1.0e-18, std::abs(parameters_.minimum_effective_diagonal));
    const double relaxation_floor =
        std::clamp(parameters_.adaptive_relaxation_floor, 0.01, 1.0);
    const bool use_quasi_newton = parameters_.solve_policy == NonlinearSolvePolicy::QuasiNewton;

    for (int iteration = 1; iteration <= parameters_.max_iterations; ++iteration)
    {
        for (std::size_t k = 0; k < nz_; ++k)
        {
            for (std::size_t j = 0; j < ny_; ++j)
            {
                for (std::size_t i = 0; i < nx_; ++i)
                {
                    if (isBoundary(i, j, k))
                    {
                        const auto point = gridPoint(i, j, k);
                        potential_[flattenIndex(i, j, k)] = boundary_ ? boundary_(point) : 0.0;
                    }
                }
            }
        }

        double max_delta = 0.0;
        double residual_sum = 0.0;
        std::size_t interior_count = 0;
        for (std::size_t k = 0; k < nz_; ++k)
        {
            for (std::size_t j = 0; j < ny_; ++j)
            {
                for (std::size_t i = 0; i < nx_; ++i)
                {
                    const std::size_t index = flattenIndex(i, j, k);
                    const auto point = gridPoint(i, j, k);
                    if (isBoundary(i, j, k))
                    {
                        potential_[index] = boundary_ ? boundary_(point) : 0.0;
                        continue;
                    }

                    const double current_phi = potential_[index];
                    const double rho = charge_density_ ? charge_density_(point, current_phi) : 0.0;
                    const double eps_r =
                        std::max(1.0, permittivity_ ? permittivity_(point, current_phi) : 1.0);

                    double effective_denom = denom;
                    if (use_quasi_newton && charge_density_)
                    {
                        const double rho_probe = charge_density_(point, current_phi + derivative_step);
                        const double drho_dphi = (rho_probe - rho) / derivative_step;
                        effective_denom -= drho_dphi / (kEpsilon0 * eps_r);
                        if (std::abs(effective_denom) < effective_diag_floor)
                        {
                            effective_denom = effective_denom < 0.0
                                                  ? -effective_diag_floor
                                                  : effective_diag_floor;
                        }
                    }

                    const double rhs =
                        (potential_[flattenIndex(i + 1, j, k)] +
                         potential_[flattenIndex(i - 1, j, k)]) /
                            dx2 +
                        (potential_[flattenIndex(i, j + 1, k)] +
                         potential_[flattenIndex(i, j - 1, k)]) /
                            dy2 +
                        (potential_[flattenIndex(i, j, k + 1)] +
                         potential_[flattenIndex(i, j, k - 1)]) /
                            dz2 +
                        rho / (kEpsilon0 * eps_r);

                    const double updated = rhs / effective_denom;
                    const double raw_delta = updated - current_phi;
                    double local_relaxation = parameters_.relaxation;
                    if (use_quasi_newton)
                    {
                        const double scale = 1.0 + std::abs(raw_delta) / std::max(1.0, std::abs(current_phi));
                        local_relaxation = std::clamp(parameters_.relaxation / scale,
                                                      relaxation_floor,
                                                      parameters_.relaxation);
                    }

                    const double relaxed = current_phi + local_relaxation * raw_delta;
                    const double step_delta = std::abs(relaxed - current_phi);
                    max_delta = std::max(max_delta, step_delta);
                    residual_sum += step_delta * step_delta;
                    ++interior_count;
                    potential_[index] = relaxed;
                }
            }
        }

        iterations_used_ = iteration;
        const double residual_norm =
            interior_count > 0
                ? std::sqrt(residual_sum / static_cast<double>(interior_count))
                : 0.0;
        last_max_delta_ = max_delta;
        last_residual_norm_ = residual_norm;

        if (max_delta <= parameters_.tolerance &&
            residual_norm <= std::max(parameters_.tolerance, 1.0e-12))
        {
            converged_ = true;
            updateElectricField();
            return true;
        }
    }

    converged_ = false;
    updateElectricField();
    return false;
}

double NonlinearPoissonSolver::getPotential(const Geometry::Point3D& position) const
{
    const auto clamp_index = [](double value, double spacing, std::size_t count) -> std::size_t
    {
        return static_cast<std::size_t>(std::clamp<long long>(
            static_cast<long long>(std::llround(value / std::max(1.0e-12, spacing))), 0,
            static_cast<long long>(count - 1)));
    };

    const std::size_t i = clamp_index(position.x(), spacing_.x(), nx_);
    const std::size_t j = clamp_index(position.y(), spacing_.y(), ny_);
    const std::size_t k = clamp_index(position.z(), spacing_.z(), nz_);
    return potential_[flattenIndex(i, j, k)];
}

std::size_t NonlinearPoissonSolver::flattenIndex(std::size_t i, std::size_t j, std::size_t k) const
{
    return i + nx_ * (j + ny_ * k);
}

Geometry::Point3D NonlinearPoissonSolver::gridPoint(std::size_t i, std::size_t j,
                                                    std::size_t k) const
{
    return Geometry::Point3D(static_cast<double>(i) * spacing_.x(),
                             static_cast<double>(j) * spacing_.y(),
                             static_cast<double>(k) * spacing_.z());
}

bool NonlinearPoissonSolver::isBoundary(std::size_t i, std::size_t j, std::size_t k) const
{
    return i == 0 || j == 0 || k == 0 || i + 1 == nx_ || j + 1 == ny_ || k + 1 == nz_;
}

void NonlinearPoissonSolver::updateElectricField()
{
    const auto safe_index = [&](std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) {
        return flattenIndex(static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(
                                i, 0, static_cast<std::ptrdiff_t>(nx_ - 1))),
                            static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(
                                j, 0, static_cast<std::ptrdiff_t>(ny_ - 1))),
                            static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(
                                k, 0, static_cast<std::ptrdiff_t>(nz_ - 1))));
    };

    for (std::size_t k = 0; k < nz_; ++k)
    {
        for (std::size_t j = 0; j < ny_; ++j)
        {
            for (std::size_t i = 0; i < nx_; ++i)
            {
                const std::size_t index = flattenIndex(i, j, k);
                const double dphidx =
                    (potential_[safe_index(static_cast<std::ptrdiff_t>(i) + 1, j, k)] -
                     potential_[safe_index(static_cast<std::ptrdiff_t>(i) - 1, j, k)]) /
                    (2.0 * std::max(1.0e-12, spacing_.x()));
                const double dphidy =
                    (potential_[safe_index(i, static_cast<std::ptrdiff_t>(j) + 1, k)] -
                     potential_[safe_index(i, static_cast<std::ptrdiff_t>(j) - 1, k)]) /
                    (2.0 * std::max(1.0e-12, spacing_.y()));
                const double dphidz =
                    (potential_[safe_index(i, j, static_cast<std::ptrdiff_t>(k) + 1)] -
                     potential_[safe_index(i, j, static_cast<std::ptrdiff_t>(k) - 1)]) /
                    (2.0 * std::max(1.0e-12, spacing_.z()));
                electric_field_[index] = Geometry::Vector3D(-dphidx, -dphidy, -dphidz);
            }
        }
    }
}

} // namespace FieldSolver
} // namespace SCDAT
