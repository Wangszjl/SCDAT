#pragma once

#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"

#include <functional>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

enum class NonlinearSolvePolicy
{
  Picard,
  QuasiNewton
};

struct NonlinearPoissonParameters
{
    double tolerance = 1.0e-6;
    int max_iterations = 500;
    double relaxation = 0.8;
  double adaptive_relaxation_floor = 0.2;
  double derivative_step_v = 1.0e-3;
  double minimum_effective_diagonal = 1.0e-12;
  NonlinearSolvePolicy solve_policy = NonlinearSolvePolicy::Picard;
};

class NonlinearPoissonSolver
{
  public:
    using ChargeDensityFunction = std::function<double(const Geometry::Point3D&, double)>;
    using PermittivityFunction = std::function<double(const Geometry::Point3D&, double)>;
    using BoundaryFunction = std::function<double(const Geometry::Point3D&)>;

    void setParameters(const NonlinearPoissonParameters& parameters) { parameters_ = parameters; }
    bool initializeGrid(const Geometry::Vector3D& grid_size, const Geometry::Vector3D& resolution);
    void setChargeDensityFunction(ChargeDensityFunction charge_density);
    void setPermittivityFunction(PermittivityFunction permittivity);
    void setBoundaryPotential(BoundaryFunction boundary);

    bool solve();

    const std::vector<double>& getPotentialField() const { return potential_; }
    const std::vector<Geometry::Vector3D>& getElectricField() const { return electric_field_; }
    double getPotential(const Geometry::Point3D& position) const;
    int getIterationsUsed() const { return iterations_used_; }
    double getLastMaxDelta() const { return last_max_delta_; }
    double getLastResidualNorm() const { return last_residual_norm_; }
    bool hasConverged() const { return converged_; }

  private:
    std::size_t flattenIndex(std::size_t i, std::size_t j, std::size_t k) const;
    Geometry::Point3D gridPoint(std::size_t i, std::size_t j, std::size_t k) const;
    bool isBoundary(std::size_t i, std::size_t j, std::size_t k) const;
    void updateElectricField();

    NonlinearPoissonParameters parameters_;
    ChargeDensityFunction charge_density_;
    PermittivityFunction permittivity_;
    BoundaryFunction boundary_;
    Geometry::Vector3D grid_size_;
    Geometry::Vector3D spacing_;
    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
    std::size_t nz_ = 0;
    std::vector<double> potential_;
    std::vector<Geometry::Vector3D> electric_field_;
    int iterations_used_ = 0;
    double last_max_delta_ = 0.0;
    double last_residual_norm_ = 0.0;
    bool converged_ = false;
    bool initialized_ = false;
};

} // namespace FieldSolver
} // namespace SCDAT
