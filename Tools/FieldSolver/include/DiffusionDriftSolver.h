#pragma once

#include "../../Basic/include/VoidResult.h"
#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"

#include <functional>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

struct DiffusionDriftParameters
{
    double mobility = 1.0e-2;
    double diffusion_coefficient = 1.0e-3;
    double source_rate = 0.0;
    double loss_rate = 0.0;
    double floor_density = 0.0;
};

struct DiffusionDriftState
{
    Geometry::Vector3D grid_size;
    Geometry::Vector3D resolution;
    Geometry::Vector3D spacing;
    std::vector<double> density;
    std::vector<Geometry::Vector3D> flux;
};

class DiffusionDriftSolver
{
  public:
    using ScalarFieldFunction = std::function<double(const Geometry::Point3D&)>;

    VoidResult setParameters(const DiffusionDriftParameters& parameters);
    VoidResult initializeGrid(const Geometry::Vector3D& grid_size,
                              const Geometry::Vector3D& resolution);
    VoidResult setInitialDensity(const ScalarFieldFunction& initial_density);
    VoidResult setBoundaryDensity(const ScalarFieldFunction& boundary_density);
    VoidResult advance(const std::vector<Geometry::Vector3D>& electric_field, double dt);

    const DiffusionDriftState& getState() const { return state_; }
    const std::vector<double>& getDensityDistribution() const { return state_.density; }
    const std::vector<Geometry::Vector3D>& getFluxDistribution() const { return state_.flux; }
    double getDensity(const Geometry::Point3D& position) const;

  private:
    std::size_t flattenIndex(std::size_t i, std::size_t j, std::size_t k) const;
    Geometry::Point3D gridPoint(std::size_t i, std::size_t j, std::size_t k) const;
    bool isBoundary(std::size_t i, std::size_t j, std::size_t k) const;

    DiffusionDriftParameters parameters_;
    DiffusionDriftState state_;
    ScalarFieldFunction boundary_density_;
    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
    std::size_t nz_ = 0;
    bool initialized_ = false;
};

} // namespace FieldSolver
} // namespace SCDAT
