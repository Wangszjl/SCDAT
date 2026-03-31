#include "ElectronTransportSolver.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace SCDAT
{
namespace FieldSolver
{

namespace
{
constexpr double kBoltzmann = 1.380649e-23;
}

ElectronTransportSolver::ElectronTransportSolver()
    : initialized_(false), grid_size_(1.0, 1.0, 1.0), resolution_(1.0, 1.0, 1.0),
      spacing_(1.0, 1.0, 1.0), total_points_(1)
{
}

ElectronTransportSolver::~ElectronTransportSolver() = default;

VoidResult ElectronTransportSolver::setParameters(const ElectronTransportParameters& params)
{
    if (params.mobility <= 0.0 || params.diffusion_coefficient <= 0.0 ||
        params.thermal_velocity <= 0.0 || params.collision_frequency <= 0.0 || params.mass <= 0.0)
    {
        return VoidResult::failure("invalid electron transport parameters");
    }
    params_ = params;
    return VoidResult::success();
}

VoidResult ElectronTransportSolver::setConfiguration(const SolverConfiguration& config)
{
    if (config.time_step <= 0.0 || config.spatial_step <= 0.0 || config.tolerance <= 0.0 ||
        config.max_iterations <= 0)
    {
        return VoidResult::failure("invalid solver configuration");
    }
    config_ = config;
    return VoidResult::success();
}

VoidResult ElectronTransportSolver::initializeGrid(const Geometry::Vector3D& grid_size,
                                                   const Geometry::Vector3D& resolution)
{
    if (resolution.x() <= 0.0 || resolution.y() <= 0.0 || resolution.z() <= 0.0)
    {
        return VoidResult::failure("resolution must be positive");
    }

    grid_size_ = grid_size;
    resolution_ = resolution;
    spacing_ = Geometry::Vector3D(grid_size.x() / resolution.x(), grid_size.y() / resolution.y(),
                               grid_size.z() / resolution.z());

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution.z()));
    total_points_ = nx * ny * nz;

    electron_states_.assign(total_points_, ElectronState{});
    electron_states_old_ = electron_states_;

    initialized_ = true;
    updateStatistics();
    return VoidResult::success();
}

VoidResult ElectronTransportSolver::setInitialConditions(
    const std::function<ElectronState(const Geometry::Point3D&)>& initial_state)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    for (size_t k = 0; k < nz; ++k)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            for (size_t i = 0; i < nx; ++i)
            {
                electron_states_[getGridIndex(i, j, k)] = initial_state(getGridPosition(i, j, k));
            }
        }
    }

    electron_states_old_ = electron_states_;
    updateStatistics();
    return VoidResult::success();
}

VoidResult ElectronTransportSolver::setBoundaryConditions(
    const std::function<ElectronState(const Geometry::Point3D&, double)>& boundary_state)
{
    boundary_function_ = boundary_state;
    return VoidResult::success();
}

TransportCoefficients
ElectronTransportSolver::calculateTransportCoefficients(const ElectronState& state,
                                                        const Geometry::Vector3D& electric_field)
{
    TransportCoefficients c;
    const double emag = electric_field.magnitude();
    c.mobility = calculateMobility(state, emag);
    c.diffusion_coefficient = calculateDiffusionCoefficient(state);
    c.thermal_conductivity = calculateThermalConductivity(state);
    c.energy_mobility = 1.5 * c.mobility;
    c.energy_diffusion = 1.5 * c.diffusion_coefficient;
    return c;
}

VoidResult
ElectronTransportSolver::solveDriftDiffusion(const std::vector<Geometry::Vector3D>& electric_field,
                                             double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }
    if (dt <= 0.0)
    {
        return VoidResult::failure("dt must be positive");
    }

    if (config_.use_implicit_scheme)
    {
        return solveImplicitScheme(electric_field, dt);
    }
    return solveExplicitScheme(electric_field, dt);
}

VoidResult
ElectronTransportSolver::solveEnergyTransport(const std::vector<Geometry::Vector3D>& electric_field,
                                              double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }
    if (dt <= 0.0)
    {
        return VoidResult::failure("dt must be positive");
    }

    for (size_t i = 0; i < electron_states_.size(); ++i)
    {
        const Geometry::Vector3D ef =
            (i < electric_field.size()) ? electric_field[i] : Geometry::Vector3D(0.0, 0.0, 0.0);
        const double collision = calculateCollisionIntegral(electron_states_[i]);
        const double joule = std::abs(params_.charge) * electron_states_[i].density *
                             ef.magnitude() * params_.mobility;

        const double dEdt =
            joule - collision * kBoltzmann * std::max(1.0, electron_states_[i].temperature);
        electron_states_[i].energy_density =
            std::max(0.0, electron_states_[i].energy_density + dt * dEdt);

        const double denom = std::max(1.0, electron_states_[i].density) * 1.5 * kBoltzmann;
        electron_states_[i].temperature = std::max(1.0, electron_states_[i].energy_density / denom);
    }

    return VoidResult::success();
}

double ElectronTransportSolver::calculateCollisionIntegral(const ElectronState& state)
{
    if (!config_.enable_collision_terms)
    {
        return 0.0;
    }

    const double nu = calculateCollisionFrequency(state);
    const double sink = nu * std::max(0.0, state.density) * 1e-4;
    return sink;
}

VoidResult ElectronTransportSolver::timeStep(const std::vector<Geometry::Vector3D>& electric_field,
                                             double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }

    electron_states_old_ = electron_states_;

    auto drift = solveDriftDiffusion(electric_field, dt);
    if (!drift.isSuccess())
    {
        return drift;
    }

    if (config_.enable_energy_equation)
    {
        auto energy = solveEnergyTransport(electric_field, dt);
        if (!energy.isSuccess())
        {
            return energy;
        }
    }

    applyBoundaryConditions(0.0);
    updateStatistics();
    return VoidResult::success();
}

ElectronState ElectronTransportSolver::getElectronState(const Geometry::Point3D& position) const
{
    if (!initialized_ || electron_states_.empty())
    {
        return ElectronState{};
    }

    const size_t i = static_cast<size_t>(
        std::clamp(position.x() / std::max(1e-30, spacing_.x()), 0.0, resolution_.x() - 1.0));
    const size_t j = static_cast<size_t>(
        std::clamp(position.y() / std::max(1e-30, spacing_.y()), 0.0, resolution_.y() - 1.0));
    const size_t k = static_cast<size_t>(
        std::clamp(position.z() / std::max(1e-30, spacing_.z()), 0.0, resolution_.z() - 1.0));

    return electron_states_[getGridIndex(i, j, k)];
}

const std::vector<ElectronState>& ElectronTransportSolver::getAllElectronStates() const
{
    return electron_states_;
}

std::vector<double> ElectronTransportSolver::getDensityDistribution() const
{
    std::vector<double> out(electron_states_.size(), 0.0);
    for (size_t i = 0; i < electron_states_.size(); ++i)
    {
        out[i] = electron_states_[i].density;
    }
    return out;
}

std::vector<double> ElectronTransportSolver::getTemperatureDistribution() const
{
    std::vector<double> out(electron_states_.size(), 0.0);
    for (size_t i = 0; i < electron_states_.size(); ++i)
    {
        out[i] = electron_states_[i].temperature;
    }
    return out;
}

std::vector<Geometry::Vector3D> ElectronTransportSolver::getCurrentDensityDistribution(
    const std::vector<Geometry::Vector3D>& electric_field) const
{
    std::vector<Geometry::Vector3D> out(electron_states_.size(), Geometry::Vector3D(0.0, 0.0, 0.0));

    for (size_t i = 0; i < electron_states_.size(); ++i)
    {
        const Geometry::Vector3D ef =
            (i < electric_field.size()) ? electric_field[i] : Geometry::Vector3D(0.0, 0.0, 0.0);
        out[i] = ef * (params_.charge * params_.mobility * electron_states_[i].density);
    }

    return out;
}

std::string ElectronTransportSolver::getStatistics() const
{
    std::ostringstream oss;
    oss << "points=" << statistics_.grid_points << ", n_avg=" << statistics_.average_density
        << ", n_peak=" << statistics_.peak_density << ", E_avg=" << statistics_.average_energy;
    return oss.str();
}

void ElectronTransportSolver::reset()
{
    initialized_ = false;
    total_points_ = 1;
    electron_states_.clear();
    electron_states_old_.clear();
    statistics_ = ElectronTransportStatistics{};
}

size_t ElectronTransportSolver::getGridIndex(size_t i, size_t j, size_t k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    return k * nx * ny + j * nx + i;
}

Geometry::Point3D ElectronTransportSolver::getGridPosition(size_t i, size_t j, size_t k) const
{
    return Geometry::Point3D((static_cast<double>(i) + 0.5) * spacing_.x(),
                          (static_cast<double>(j) + 0.5) * spacing_.y(),
                          (static_cast<double>(k) + 0.5) * spacing_.z());
}

void ElectronTransportSolver::getGridIndices(size_t index, size_t& i, size_t& j, size_t& k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));

    k = index / (nx * ny);
    const size_t rem = index % (nx * ny);
    j = rem / nx;
    i = rem % nx;
}

VoidResult
ElectronTransportSolver::solveExplicitScheme(const std::vector<Geometry::Vector3D>& electric_field,
                                             double dt)
{
    std::vector<double> density(electron_states_.size(), 0.0);
    for (size_t i = 0; i < electron_states_.size(); ++i)
    {
        density[i] = electron_states_[i].density;
    }

    for (size_t i = 0; i < electron_states_.size(); ++i)
    {
        const auto grad = calculateGradient(i, density);
        const Geometry::Vector3D ef =
            (i < electric_field.size()) ? electric_field[i] : Geometry::Vector3D(0.0, 0.0, 0.0);
        const double source = -calculateCollisionIntegral(electron_states_[i]);

        const Geometry::Vector3D flux = ef * (params_.mobility * electron_states_[i].density) -
                                     grad * params_.diffusion_coefficient;
        electron_states_[i].flux = flux;

        const double div =
            calculateDivergence(i, std::vector<Geometry::Vector3D>(electron_states_.size(), flux));
        electron_states_[i].density =
            std::max(0.0, electron_states_[i].density + dt * (source - div));
        electron_states_[i].velocity = (std::abs(electron_states_[i].density) > 1e-30)
                                           ? (flux / electron_states_[i].density)
                                           : Geometry::Vector3D(0.0, 0.0, 0.0);
    }

    return VoidResult::success();
}

VoidResult
ElectronTransportSolver::solveImplicitScheme(const std::vector<Geometry::Vector3D>& electric_field,
                                             double dt)
{
    // Minimal robust fallback: perform two explicit sweeps for damping.
    auto r1 = solveExplicitScheme(electric_field, dt * 0.5);
    if (!r1.isSuccess())
    {
        return r1;
    }
    return solveExplicitScheme(electric_field, dt * 0.5);
}

double ElectronTransportSolver::calculateLaplacian(size_t index,
                                                   const std::vector<double>& field) const
{
    if (field.empty() || index >= field.size())
    {
        return 0.0;
    }

    size_t i = 0;
    size_t j = 0;
    size_t k = 0;
    getGridIndices(index, i, j, k);

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    auto at = [&](size_t ii, size_t jj, size_t kk) -> double
    {
        ii = std::min(ii, nx - 1);
        jj = std::min(jj, ny - 1);
        kk = std::min(kk, nz - 1);
        return field[getGridIndex(ii, jj, kk)];
    };

    const double c = at(i, j, k);
    const double dx = std::max(1e-30, spacing_.x());
    const double dy = std::max(1e-30, spacing_.y());
    const double dz = std::max(1e-30, spacing_.z());

    const double d2x =
        (at(std::min(i + 1, nx - 1), j, k) - 2.0 * c + at((i == 0) ? 0 : (i - 1), j, k)) /
        (dx * dx);
    const double d2y =
        (at(i, std::min(j + 1, ny - 1), k) - 2.0 * c + at(i, (j == 0) ? 0 : (j - 1), k)) /
        (dy * dy);
    const double d2z =
        (at(i, j, std::min(k + 1, nz - 1)) - 2.0 * c + at(i, j, (k == 0) ? 0 : (k - 1))) /
        (dz * dz);
    return d2x + d2y + d2z;
}

Geometry::Vector3D ElectronTransportSolver::calculateGradient(size_t index,
                                                           const std::vector<double>& field) const
{
    if (field.empty() || index >= field.size())
    {
        return Geometry::Vector3D(0.0, 0.0, 0.0);
    }

    size_t i = 0;
    size_t j = 0;
    size_t k = 0;
    getGridIndices(index, i, j, k);

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    auto at = [&](size_t ii, size_t jj, size_t kk) -> double
    {
        ii = std::min(ii, nx - 1);
        jj = std::min(jj, ny - 1);
        kk = std::min(kk, nz - 1);
        return field[getGridIndex(ii, jj, kk)];
    };

    const double dx = std::max(1e-30, spacing_.x());
    const double dy = std::max(1e-30, spacing_.y());
    const double dz = std::max(1e-30, spacing_.z());

    const double gx =
        (at(std::min(i + 1, nx - 1), j, k) - at((i == 0) ? 0 : (i - 1), j, k)) / (2.0 * dx);
    const double gy =
        (at(i, std::min(j + 1, ny - 1), k) - at(i, (j == 0) ? 0 : (j - 1), k)) / (2.0 * dy);
    const double gz =
        (at(i, j, std::min(k + 1, nz - 1)) - at(i, j, (k == 0) ? 0 : (k - 1))) / (2.0 * dz);

    return Geometry::Vector3D(gx, gy, gz);
}

double ElectronTransportSolver::calculateDivergence(size_t index,
                                                    const std::vector<Geometry::Vector3D>& field) const
{
    if (field.empty() || index >= field.size())
    {
        return 0.0;
    }

    size_t i = 0;
    size_t j = 0;
    size_t k = 0;
    getGridIndices(index, i, j, k);

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    auto at = [&](size_t ii, size_t jj, size_t kk) -> Geometry::Vector3D
    {
        ii = std::min(ii, nx - 1);
        jj = std::min(jj, ny - 1);
        kk = std::min(kk, nz - 1);
        return field[getGridIndex(ii, jj, kk)];
    };

    const double dx = std::max(1e-30, spacing_.x());
    const double dy = std::max(1e-30, spacing_.y());
    const double dz = std::max(1e-30, spacing_.z());

    const auto px = at(std::min(i + 1, nx - 1), j, k);
    const auto mx = at((i == 0) ? 0 : (i - 1), j, k);
    const auto py = at(i, std::min(j + 1, ny - 1), k);
    const auto my = at(i, (j == 0) ? 0 : (j - 1), k);
    const auto pz = at(i, j, std::min(k + 1, nz - 1));
    const auto mz = at(i, j, (k == 0) ? 0 : (k - 1));

    return (px.x() - mx.x()) / (2.0 * dx) + (py.y() - my.y()) / (2.0 * dy) +
           (pz.z() - mz.z()) / (2.0 * dz);
}

double ElectronTransportSolver::calculateMobility(const ElectronState& /*state*/,
                                                  double electric_field_magnitude)
{
    return params_.mobility / (1.0 + electric_field_magnitude * 1e-6);
}

double ElectronTransportSolver::calculateDiffusionCoefficient(const ElectronState& state)
{
    return std::max(params_.diffusion_coefficient, calculateMobility(state, 0.0) * kBoltzmann *
                                                       state.temperature /
                                                       std::abs(params_.charge));
}

double ElectronTransportSolver::calculateThermalConductivity(const ElectronState& state)
{
    return 1.5 * kBoltzmann * std::max(0.0, state.density) * calculateDiffusionCoefficient(state);
}

double ElectronTransportSolver::calculateCollisionFrequency(const ElectronState& state)
{
    return params_.collision_frequency * (1.0 + std::max(0.0, state.density) * 1e-20);
}

void ElectronTransportSolver::applyBoundaryConditions(double time)
{
    if (!boundary_function_)
    {
        return;
    }

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    for (size_t k = 0; k < nz; ++k)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            for (size_t i = 0; i < nx; ++i)
            {
                if (!isBoundaryPoint(i, j, k))
                {
                    continue;
                }
                const size_t idx = getGridIndex(i, j, k);
                electron_states_[idx] = boundary_function_(getGridPosition(i, j, k), time);
            }
        }
    }
}

bool ElectronTransportSolver::isBoundaryPoint(size_t i, size_t j, size_t k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));
    return i == 0 || j == 0 || k == 0 || i + 1 == nx || j + 1 == ny || k + 1 == nz;
}

void ElectronTransportSolver::updateStatistics()
{
    statistics_ = ElectronTransportStatistics{};
    statistics_.grid_points = electron_states_.size();

    if (electron_states_.empty())
    {
        return;
    }

    for (const auto& st : electron_states_)
    {
        statistics_.total_electrons += st.density;
        statistics_.average_density += st.density;
        statistics_.peak_density = std::max(statistics_.peak_density, st.density);
        statistics_.average_energy += st.energy_density;
        statistics_.total_current += st.flux.magnitude() * std::abs(params_.charge);
    }

    const double n = static_cast<double>(electron_states_.size());
    statistics_.average_density /= n;
    statistics_.average_energy /= n;
}

double ElectronTransportSolver::interpolateField(const Geometry::Point3D& position,
                                                 const std::vector<double>& field) const
{
    if (field.empty())
    {
        return 0.0;
    }

    const size_t i = static_cast<size_t>(
        std::clamp(position.x() / std::max(1e-30, spacing_.x()), 0.0, resolution_.x() - 1.0));
    const size_t j = static_cast<size_t>(
        std::clamp(position.y() / std::max(1e-30, spacing_.y()), 0.0, resolution_.y() - 1.0));
    const size_t k = static_cast<size_t>(
        std::clamp(position.z() / std::max(1e-30, spacing_.z()), 0.0, resolution_.z() - 1.0));
    const size_t idx = getGridIndex(i, j, k);
    return (idx < field.size()) ? field[idx] : field.front();
}

Geometry::Vector3D
ElectronTransportSolver::interpolateVectorField(const Geometry::Point3D& position,
                                                const std::vector<Geometry::Vector3D>& field) const
{
    if (field.empty())
    {
        return Geometry::Vector3D(0.0, 0.0, 0.0);
    }

    const size_t i = static_cast<size_t>(
        std::clamp(position.x() / std::max(1e-30, spacing_.x()), 0.0, resolution_.x() - 1.0));
    const size_t j = static_cast<size_t>(
        std::clamp(position.y() / std::max(1e-30, spacing_.y()), 0.0, resolution_.y() - 1.0));
    const size_t k = static_cast<size_t>(
        std::clamp(position.z() / std::max(1e-30, spacing_.z()), 0.0, resolution_.z() - 1.0));
    const size_t idx = getGridIndex(i, j, k);
    return (idx < field.size()) ? field[idx] : field.front();
}

} // namespace FieldSolver
} // namespace SCDAT
