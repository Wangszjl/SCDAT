#pragma once

#include "FluidAlgorithmConfig.h"

#include "../../Tools/FieldSolver/include/DiffusionDriftSolver.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

class DiffusionSolverAdapter
{
  public:
    bool initialize(const FluidAlgorithmConfig& config, double initial_density_m3);
    bool advance(const std::vector<Geometry::Vector3D>& electric_field, double dt);

    const std::vector<double>& getDensityDistribution() const;
    const std::vector<Geometry::Vector3D>& getFluxDistribution() const;
    const FieldSolver::DiffusionDriftState& getState() const;

  private:
    FieldSolver::DiffusionDriftSolver solver_;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
