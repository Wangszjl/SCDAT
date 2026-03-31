#include "IonTransportSolver.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace SCDAT
{
namespace FieldSolver
{

IonTransportSolver::IonTransportSolver()
    : initialized_(false), grid_size_(1.0, 1.0, 1.0), resolution_(1.0, 1.0, 1.0),
      spacing_(1.0, 1.0, 1.0), total_points_(1)
{
}

IonTransportSolver::~IonTransportSolver() = default;

VoidResult IonTransportSolver::setParameters(const IonTransportParameters& params)
{
    if (params.mobility <= 0.0 || params.diffusion_coefficient <= 0.0 ||
        params.thermal_velocity <= 0.0 || params.collision_frequency <= 0.0 || params.mass <= 0.0)
    {
        return VoidResult::failure("invalid ion transport parameters");
    }

    params_ = params;
    return VoidResult::success();
}

VoidResult IonTransportSolver::initializeGrid(const Utils::Vector3D& grid_size,
                                              const Utils::Vector3D& resolution)
{
    if (resolution.x() <= 0.0 || resolution.y() <= 0.0 || resolution.z() <= 0.0)
    {
        return VoidResult::failure("resolution must be positive");
    }

    grid_size_ = grid_size;
    resolution_ = resolution;
    spacing_ = Utils::Vector3D(grid_size.x() / resolution.x(), grid_size.y() / resolution.y(),
                               grid_size.z() / resolution.z());

    const size_t nx = static_cast<size_t>(std::max(1.0, resolution.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution.z()));
    total_points_ = nx * ny * nz;

    ion_states_.assign(total_points_, IonState{});
    ion_states_old_ = ion_states_;

    initialized_ = true;
    updateStatistics();
    return VoidResult::success();
}

VoidResult IonTransportSolver::setInitialConditions(
    const std::function<IonState(const Utils::Point3D&)>& initial_state)
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
                ion_states_[getGridIndex(i, j, k)] = initial_state(getGridPosition(i, j, k));
            }
        }
    }

    ion_states_old_ = ion_states_;
    updateStatistics();
    return VoidResult::success();
}

IonTransportCoefficients
IonTransportSolver::calculateTransportCoefficients(const IonState& state,
                                                   const Utils::Vector3D& electric_field)
{
    IonTransportCoefficients c;
    const double emag = electric_field.magnitude();
    c.mobility = calculateIonMobility(state, emag);
    c.diffusion_coefficient = calculateIonDiffusionCoefficient(state);
    c.thermal_conductivity = 1.5 * c.diffusion_coefficient * std::max(state.density, 0.0) * 1e-23;
    c.viscosity = 0.1 * params_.mass * std::max(state.density, 0.0) * c.mobility;
    return c;
}

VoidResult
IonTransportSolver::solveDriftDiffusion(const std::vector<Utils::Vector3D>& electric_field,
                                        const std::vector<double>& electron_density, double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }
    if (dt <= 0.0)
    {
        return VoidResult::failure("dt must be positive");
    }

    std::vector<double> density(ion_states_.size(), 0.0);
    for (size_t i = 0; i < ion_states_.size(); ++i)
    {
        density[i] = ion_states_[i].density;
    }

    for (size_t i = 0; i < ion_states_.size(); ++i)
    {
        const Utils::Vector3D ef =
            (i < electric_field.size()) ? electric_field[i] : Utils::Vector3D(0.0, 0.0, 0.0);
        const double ne = (i < electron_density.size()) ? std::max(0.0, electron_density[i]) : 0.0;

        const auto grad_n = calculateGradient(i, density);
        const auto coeff = calculateTransportCoefficients(ion_states_[i], ef);

        const Utils::Vector3D drift = ef * (coeff.mobility * ion_states_[i].density);
        const Utils::Vector3D diff = grad_n * coeff.diffusion_coefficient;
        const Utils::Vector3D flux = drift - diff;
        ion_states_[i].flux = flux;

        const double recomb = calculateRecombinationRate(ion_states_[i], ne);
        const double h = std::max({1e-30, spacing_.x(), spacing_.y(), spacing_.z()});
        const double div = flux.magnitude() / h;

        ion_states_[i].density = std::max(0.0, ion_states_[i].density + dt * (-div - recomb));
        ion_states_[i].velocity = (ion_states_[i].density > 1e-30) ? (flux / ion_states_[i].density)
                                                                   : Utils::Vector3D(0.0, 0.0, 0.0);
    }

    return VoidResult::success();
}

double IonTransportSolver::calculateRecombinationRate(const IonState& ion_state,
                                                      double electron_density)
{
    const double radiative =
        1e-19 * std::max(0.0, ion_state.density) * std::max(0.0, electron_density);
    const double three_body =
        1e-39 * std::max(0.0, ion_state.density) * std::pow(std::max(0.0, electron_density), 2.0);
    return radiative + three_body;
}

double IonTransportSolver::calculateIonNeutralCollision(const IonState& ion_state,
                                                        double neutral_density)
{
    return calculateIonCollisionFrequency(ion_state, neutral_density);
}

VoidResult IonTransportSolver::timeStep(const std::vector<Utils::Vector3D>& electric_field,
                                        const std::vector<double>& electron_density,
                                        const std::vector<double>& neutral_density, double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }

    ion_states_old_ = ion_states_;

    auto res = solveDriftDiffusion(electric_field, electron_density, dt);
    if (!res.isSuccess())
    {
        return res;
    }

    for (size_t i = 0; i < ion_states_.size(); ++i)
    {
        const double nn = (i < neutral_density.size()) ? std::max(0.0, neutral_density[i]) : 0.0;
        const double nu = calculateIonNeutralCollision(ion_states_[i], nn);
        ion_states_[i].temperature = std::max(1.0, ion_states_[i].temperature - dt * 1e-3 * nu);
    }

    applyBoundaryConditions();
    updateStatistics();
    return VoidResult::success();
}

IonState IonTransportSolver::getIonState(const Utils::Point3D& position) const
{
    if (!initialized_ || ion_states_.empty())
    {
        return IonState{};
    }

    const size_t i = static_cast<size_t>(
        std::clamp(position.x() / std::max(1e-30, spacing_.x()), 0.0, resolution_.x() - 1.0));
    const size_t j = static_cast<size_t>(
        std::clamp(position.y() / std::max(1e-30, spacing_.y()), 0.0, resolution_.y() - 1.0));
    const size_t k = static_cast<size_t>(
        std::clamp(position.z() / std::max(1e-30, spacing_.z()), 0.0, resolution_.z() - 1.0));

    return ion_states_[getGridIndex(i, j, k)];
}

const std::vector<IonState>& IonTransportSolver::getAllIonStates() const
{
    return ion_states_;
}

std::vector<double> IonTransportSolver::getDensityDistribution() const
{
    std::vector<double> out(ion_states_.size(), 0.0);
    for (size_t i = 0; i < ion_states_.size(); ++i)
    {
        out[i] = ion_states_[i].density;
    }
    return out;
}

std::vector<Utils::Vector3D> IonTransportSolver::getVelocityDistribution() const
{
    std::vector<Utils::Vector3D> out(ion_states_.size(), Utils::Vector3D(0.0, 0.0, 0.0));
    for (size_t i = 0; i < ion_states_.size(); ++i)
    {
        out[i] = ion_states_[i].velocity;
    }
    return out;
}

std::vector<Utils::Vector3D> IonTransportSolver::getCurrentDensityDistribution() const
{
    std::vector<Utils::Vector3D> out(ion_states_.size(), Utils::Vector3D(0.0, 0.0, 0.0));
    for (size_t i = 0; i < ion_states_.size(); ++i)
    {
        out[i] = ion_states_[i].flux * (params_.charge * ion_states_[i].charge_state);
    }
    return out;
}

std::string IonTransportSolver::getStatistics() const
{
    std::ostringstream oss;
    oss << "points=" << statistics_.grid_points << ", n_avg=" << statistics_.average_density
        << ", n_peak=" << statistics_.peak_density << ", v_avg=" << statistics_.average_velocity;
    return oss.str();
}

void IonTransportSolver::reset()
{
    initialized_ = false;
    total_points_ = 1;
    ion_states_.clear();
    ion_states_old_.clear();
    statistics_ = IonTransportStatistics{};
}

size_t IonTransportSolver::getGridIndex(size_t i, size_t j, size_t k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    return k * nx * ny + j * nx + i;
}

Utils::Point3D IonTransportSolver::getGridPosition(size_t i, size_t j, size_t k) const
{
    return Utils::Point3D((static_cast<double>(i) + 0.5) * spacing_.x(),
                          (static_cast<double>(j) + 0.5) * spacing_.y(),
                          (static_cast<double>(k) + 0.5) * spacing_.z());
}

void IonTransportSolver::getGridIndices(size_t index, size_t& i, size_t& j, size_t& k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));

    k = index / (nx * ny);
    const size_t rem = index % (nx * ny);
    j = rem / nx;
    i = rem % nx;
}

double IonTransportSolver::calculateLaplacian(size_t index, const std::vector<double>& field) const
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

Utils::Vector3D IonTransportSolver::calculateGradient(size_t index,
                                                      const std::vector<double>& field) const
{
    if (field.empty() || index >= field.size())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
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

    return Utils::Vector3D(gx, gy, gz);
}

double IonTransportSolver::calculateDivergence(size_t index,
                                               const std::vector<Utils::Vector3D>& field) const
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

    auto at = [&](size_t ii, size_t jj, size_t kk) -> Utils::Vector3D
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

double IonTransportSolver::calculateIonMobility(const IonState& /*state*/,
                                                double electric_field_magnitude)
{
    return params_.mobility / (1.0 + electric_field_magnitude * 1e-6);
}

double IonTransportSolver::calculateIonDiffusionCoefficient(const IonState& state)
{
    return std::max(params_.diffusion_coefficient, 1e-6 * std::max(1.0, state.temperature / 300.0));
}

double IonTransportSolver::calculateIonCollisionFrequency(const IonState& state,
                                                          double neutral_density)
{
    return params_.collision_frequency * (1.0 + neutral_density * 1e-20 + state.density * 1e-21);
}

double IonTransportSolver::calculateChargeExchangeRate(const IonState& state,
                                                       double neutral_density)
{
    return 1e-15 * neutral_density * std::max(state.velocity.magnitude(), 0.0);
}

void IonTransportSolver::applyBoundaryConditions()
{
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
                IonState& st = ion_states_[getGridIndex(i, j, k)];
                st.density = std::max(0.0, st.density);
                st.temperature = std::max(1.0, st.temperature);
            }
        }
    }
}

bool IonTransportSolver::isBoundaryPoint(size_t i, size_t j, size_t k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    return i == 0 || j == 0 || k == 0 || i + 1 == nx || j + 1 == ny || k + 1 == nz;
}

void IonTransportSolver::updateStatistics()
{
    statistics_ = IonTransportStatistics{};
    statistics_.grid_points = ion_states_.size();

    if (ion_states_.empty())
    {
        return;
    }

    for (const auto& st : ion_states_)
    {
        statistics_.total_ions += st.density;
        statistics_.average_density += st.density;
        statistics_.peak_density = std::max(statistics_.peak_density, st.density);
        statistics_.average_velocity += st.velocity.magnitude();
        statistics_.total_current += st.flux.magnitude() * params_.charge;
    }

    const double n = static_cast<double>(ion_states_.size());
    statistics_.average_density /= n;
    statistics_.average_velocity /= n;
}

double IonTransportSolver::interpolateField(const Utils::Point3D& position,
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

} // namespace FieldSolver
} // namespace SCDAT
