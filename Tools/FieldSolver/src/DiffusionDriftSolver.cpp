#include "../include/DiffusionDriftSolver.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace FieldSolver
{

VoidResult DiffusionDriftSolver::setParameters(const DiffusionDriftParameters& parameters)
{
    if (parameters.diffusion_coefficient < 0.0 || parameters.mobility < 0.0)
    {
        return VoidResult::failure("DiffusionDriftSolver parameters are invalid");
    }
    parameters_ = parameters;
    return VoidResult::success();
}

VoidResult DiffusionDriftSolver::initializeGrid(const Geometry::Vector3D& grid_size,
                                                const Geometry::Vector3D& resolution)
{
    nx_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(resolution.x())));
    ny_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(resolution.y())));
    nz_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(resolution.z())));

    state_.grid_size = grid_size;
    state_.resolution = Geometry::Vector3D(static_cast<double>(nx_), static_cast<double>(ny_),
                                           static_cast<double>(nz_));
    state_.spacing = Geometry::Vector3D(
        grid_size.x() / static_cast<double>(std::max<std::size_t>(1, nx_ - 1)),
        grid_size.y() / static_cast<double>(std::max<std::size_t>(1, ny_ - 1)),
        grid_size.z() / static_cast<double>(std::max<std::size_t>(1, nz_ - 1)));
    state_.density.assign(nx_ * ny_ * nz_, 0.0);
    state_.flux.assign(nx_ * ny_ * nz_, Geometry::Vector3D{});
    initialized_ = true;
    return VoidResult::success();
}

VoidResult DiffusionDriftSolver::setInitialDensity(const ScalarFieldFunction& initial_density)
{
    if (!initialized_)
    {
        return VoidResult::failure("DiffusionDriftSolver grid is not initialized");
    }

    for (std::size_t k = 0; k < nz_; ++k)
    {
        for (std::size_t j = 0; j < ny_; ++j)
        {
            for (std::size_t i = 0; i < nx_; ++i)
            {
                state_.density[flattenIndex(i, j, k)] = initial_density(gridPoint(i, j, k));
            }
        }
    }
    return VoidResult::success();
}

VoidResult DiffusionDriftSolver::setBoundaryDensity(const ScalarFieldFunction& boundary_density)
{
    boundary_density_ = boundary_density;
    return VoidResult::success();
}

VoidResult DiffusionDriftSolver::advance(const std::vector<Geometry::Vector3D>& electric_field,
                                         double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("DiffusionDriftSolver grid is not initialized");
    }
    if (electric_field.size() != state_.density.size())
    {
        return VoidResult::failure("DiffusionDriftSolver electric field size mismatch");
    }

    const auto old_density = state_.density;
    const double dx = std::max(1.0e-12, state_.spacing.x());
    const double dy = std::max(1.0e-12, state_.spacing.y());
    const double dz = std::max(1.0e-12, state_.spacing.z());

    for (std::size_t k = 0; k < nz_; ++k)
    {
        for (std::size_t j = 0; j < ny_; ++j)
        {
            for (std::size_t i = 0; i < nx_; ++i)
            {
                const std::size_t index = flattenIndex(i, j, k);
                if (isBoundary(i, j, k))
                {
                    if (boundary_density_)
                    {
                        state_.density[index] = boundary_density_(gridPoint(i, j, k));
                    }
                    continue;
                }

                const auto idx = [&](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return flattenIndex(ii, jj, kk);
                };

                const double center = old_density[index];
                const double laplacian =
                    (old_density[idx(i + 1, j, k)] - 2.0 * center + old_density[idx(i - 1, j, k)]) /
                        (dx * dx) +
                    (old_density[idx(i, j + 1, k)] - 2.0 * center + old_density[idx(i, j - 1, k)]) /
                        (dy * dy) +
                    (old_density[idx(i, j, k + 1)] - 2.0 * center + old_density[idx(i, j, k - 1)]) /
                        (dz * dz);

                const Geometry::Vector3D density_gradient(
                    (old_density[idx(i + 1, j, k)] - old_density[idx(i - 1, j, k)]) / (2.0 * dx),
                    (old_density[idx(i, j + 1, k)] - old_density[idx(i, j - 1, k)]) / (2.0 * dy),
                    (old_density[idx(i, j, k + 1)] - old_density[idx(i, j, k - 1)]) / (2.0 * dz));
                const double drift_term = electric_field[index].dot(density_gradient);

                const double updated =
                    center + dt * (parameters_.diffusion_coefficient * laplacian -
                                   parameters_.mobility * drift_term + parameters_.source_rate -
                                   parameters_.loss_rate * center);
                state_.density[index] = std::max(parameters_.floor_density, updated);
            }
        }
    }

    for (std::size_t k = 0; k < nz_; ++k)
    {
        for (std::size_t j = 0; j < ny_; ++j)
        {
            for (std::size_t i = 0; i < nx_; ++i)
            {
                const std::size_t index = flattenIndex(i, j, k);
                const auto safe = [&](std::ptrdiff_t ii, std::ptrdiff_t jj, std::ptrdiff_t kk) {
                    return state_.density[flattenIndex(static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(
                        ii, 0, static_cast<std::ptrdiff_t>(nx_ - 1))),
                                                       static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(
                        jj, 0, static_cast<std::ptrdiff_t>(ny_ - 1))),
                                                       static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(
                        kk, 0, static_cast<std::ptrdiff_t>(nz_ - 1))))];
                };

                const Geometry::Vector3D density_gradient(
                    (safe(static_cast<std::ptrdiff_t>(i) + 1, j, k) -
                     safe(static_cast<std::ptrdiff_t>(i) - 1, j, k)) /
                        (2.0 * dx),
                    (safe(i, static_cast<std::ptrdiff_t>(j) + 1, k) -
                     safe(i, static_cast<std::ptrdiff_t>(j) - 1, k)) /
                        (2.0 * dy),
                    (safe(i, j, static_cast<std::ptrdiff_t>(k) + 1) -
                     safe(i, j, static_cast<std::ptrdiff_t>(k) - 1)) /
                        (2.0 * dz));

                state_.flux[index] = electric_field[index] * (parameters_.mobility * state_.density[index]) -
                                     density_gradient * parameters_.diffusion_coefficient;
            }
        }
    }

    return VoidResult::success();
}

double DiffusionDriftSolver::getDensity(const Geometry::Point3D& position) const
{
    const auto clamp_index = [](double value, double spacing, std::size_t count) -> std::size_t
    {
        return static_cast<std::size_t>(std::clamp<long long>(
            static_cast<long long>(std::llround(value / std::max(1.0e-12, spacing))), 0,
            static_cast<long long>(count - 1)));
    };

    const std::size_t i = clamp_index(position.x(), state_.spacing.x(), nx_);
    const std::size_t j = clamp_index(position.y(), state_.spacing.y(), ny_);
    const std::size_t k = clamp_index(position.z(), state_.spacing.z(), nz_);
    return state_.density[flattenIndex(i, j, k)];
}

std::size_t DiffusionDriftSolver::flattenIndex(std::size_t i, std::size_t j, std::size_t k) const
{
    return i + nx_ * (j + ny_ * k);
}

Geometry::Point3D DiffusionDriftSolver::gridPoint(std::size_t i, std::size_t j, std::size_t k) const
{
    return Geometry::Point3D(static_cast<double>(i) * state_.spacing.x(),
                             static_cast<double>(j) * state_.spacing.y(),
                             static_cast<double>(k) * state_.spacing.z());
}

bool DiffusionDriftSolver::isBoundary(std::size_t i, std::size_t j, std::size_t k) const
{
    return i == 0 || j == 0 || k == 0 || i + 1 == nx_ || j + 1 == ny_ || k + 1 == nz_;
}

} // namespace FieldSolver
} // namespace SCDAT
