#include "DiffusionSolverAdapter.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

bool DiffusionSolverAdapter::initialize(const FluidAlgorithmConfig& config, double initial_density_m3)
{
    FieldSolver::DiffusionDriftParameters parameters;
    parameters.mobility = 0.08;
    parameters.diffusion_coefficient = 0.015;
    parameters.floor_density = 1.0e10;
    if (!solver_.setParameters(parameters))
    {
        return false;
    }
    if (!solver_.initializeGrid(config.domain_size, config.resolution))
    {
        return false;
    }
    if (!solver_.setInitialDensity(
            [initial_density_m3](const Geometry::Point3D&) { return initial_density_m3; }))
    {
        return false;
    }
    solver_.setBoundaryDensity(
        [initial_density_m3](const Geometry::Point3D&) { return 0.95 * initial_density_m3; });
    return true;
}

bool DiffusionSolverAdapter::advance(const std::vector<Geometry::Vector3D>& electric_field, double dt)
{
    return static_cast<bool>(solver_.advance(electric_field, dt));
}

const std::vector<double>& DiffusionSolverAdapter::getDensityDistribution() const
{
    return solver_.getDensityDistribution();
}

const std::vector<Geometry::Vector3D>& DiffusionSolverAdapter::getFluxDistribution() const
{
    return solver_.getFluxDistribution();
}

const FieldSolver::DiffusionDriftState& DiffusionSolverAdapter::getState() const
{
    return solver_.getState();
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
