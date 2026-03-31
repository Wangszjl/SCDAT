#include "EnergyBalanceSolver.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace SCDAT
{
namespace FieldSolver
{

namespace
{
constexpr double kBoltzmann = 1.380649e-23;
}

EnergyBalanceSolver::EnergyBalanceSolver()
    : initialized_(false), grid_size_(1.0, 1.0, 1.0), resolution_(1.0, 1.0, 1.0),
      spacing_(1.0, 1.0, 1.0), total_points_(1)
{
}

EnergyBalanceSolver::~EnergyBalanceSolver() = default;

VoidResult EnergyBalanceSolver::setParameters(const EnergyBalanceParameters& params)
{
    if (params.electron_heat_capacity <= 0.0 || params.ion_heat_capacity <= 0.0)
    {
        return VoidResult::failure("invalid heat capacity");
    }
    if (params.electron_thermal_conductivity <= 0.0 || params.ion_thermal_conductivity <= 0.0)
    {
        return VoidResult::failure("invalid thermal conductivity");
    }

    params_ = params;
    return VoidResult::success();
}

VoidResult EnergyBalanceSolver::initializeGrid(const Utils::Vector3D& grid_size,
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

    const auto nx = static_cast<size_t>(std::max(1.0, resolution.x()));
    const auto ny = static_cast<size_t>(std::max(1.0, resolution.y()));
    const auto nz = static_cast<size_t>(std::max(1.0, resolution.z()));
    total_points_ = nx * ny * nz;

    energy_states_.assign(total_points_, EnergyState{});
    energy_states_old_ = energy_states_;

    initialized_ = true;
    updateStatistics();
    return VoidResult::success();
}

VoidResult EnergyBalanceSolver::setInitialConditions(
    const std::function<EnergyState(const Utils::Point3D&)>& initial_state)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }

    const auto nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const auto ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const auto nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    for (size_t k = 0; k < nz; ++k)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            for (size_t i = 0; i < nx; ++i)
            {
                const size_t idx = getGridIndex(i, j, k);
                energy_states_[idx] = initial_state(getGridPosition(i, j, k));
            }
        }
    }

    energy_states_old_ = energy_states_;
    updateStatistics();
    return VoidResult::success();
}

EnergySourceTerms EnergyBalanceSolver::calculateEnergySourceTerms(
    const Utils::Point3D& /*position*/, double electron_density, double ion_density,
    const Utils::Vector3D& electric_field, const Utils::Vector3D& current_density)
{
    EnergySourceTerms src;

    const double jmag = current_density.magnitude();
    const double emag = electric_field.magnitude();
    const double n = std::max(0.0, 0.5 * (electron_density + ion_density));

    src.electron_joule_heating =
        calculateJouleHeating(jmag, std::max(params_.electron_thermal_conductivity, 1e-20));
    src.ion_joule_heating = 0.1 * src.electron_joule_heating;
    src.elastic_energy_exchange = calculateEnergyExchange(1000.0, 900.0, n);
    src.inelastic_energy_loss =
        params_.inelastic_collision_frequency * n * kBoltzmann * std::max(300.0, emag * 1e-3);
    src.radiation_loss = calculateRadiationLoss(1000.0, 900.0, n);
    src.thermal_conduction_loss =
        (params_.electron_thermal_conductivity + params_.ion_thermal_conductivity) * emag * 1e-3;

    return src;
}

VoidResult EnergyBalanceSolver::solveElectronEnergyBalance(
    const std::vector<double>& electron_density, const std::vector<Utils::Vector3D>& electric_field,
    const std::vector<Utils::Vector3D>& current_density, double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }
    if (dt <= 0.0)
    {
        return VoidResult::failure("dt must be positive");
    }

    const auto conduction = calculateThermalConduction(getElectronTemperatureDistribution(),
                                                       params_.electron_thermal_conductivity);

    for (size_t i = 0; i < energy_states_.size(); ++i)
    {
        const double ne = (i < electron_density.size()) ? std::max(0.0, electron_density[i]) : 0.0;
        const Utils::Vector3D ef =
            (i < electric_field.size()) ? electric_field[i] : Utils::Vector3D(0.0, 0.0, 0.0);
        const Utils::Vector3D jf =
            (i < current_density.size()) ? current_density[i] : Utils::Vector3D(0.0, 0.0, 0.0);

        EnergySourceTerms src = calculateEnergySourceTerms(Utils::Point3D(), ne, ne, ef, jf);
        const double cap =
            std::max(1e-30, params_.electron_heat_capacity * kBoltzmann * std::max(ne, 1.0));
        const double dEdt = src.electron_joule_heating - src.inelastic_energy_loss -
                            src.radiation_loss - conduction[i];

        energy_states_[i].electron_energy_density =
            std::max(0.0, energy_states_[i].electron_energy_density + dt * dEdt);
        energy_states_[i].electron_temperature =
            std::max(1.0, energy_states_[i].electron_energy_density / cap);
        energy_states_[i].joule_heating_rate = src.electron_joule_heating;
        energy_states_[i].inelastic_cooling_rate = src.inelastic_energy_loss;
        energy_states_[i].radiation_loss_rate = src.radiation_loss;
    }

    return VoidResult::success();
}

VoidResult EnergyBalanceSolver::solveIonEnergyBalance(
    const std::vector<double>& ion_density, const std::vector<Utils::Vector3D>& electric_field,
    const std::vector<Utils::Vector3D>& current_density, double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }
    if (dt <= 0.0)
    {
        return VoidResult::failure("dt must be positive");
    }

    const auto conduction = calculateThermalConduction(getIonTemperatureDistribution(),
                                                       params_.ion_thermal_conductivity);

    for (size_t i = 0; i < energy_states_.size(); ++i)
    {
        const double ni = (i < ion_density.size()) ? std::max(0.0, ion_density[i]) : 0.0;
        const Utils::Vector3D ef =
            (i < electric_field.size()) ? electric_field[i] : Utils::Vector3D(0.0, 0.0, 0.0);
        const Utils::Vector3D jf =
            (i < current_density.size()) ? current_density[i] : Utils::Vector3D(0.0, 0.0, 0.0);

        EnergySourceTerms src = calculateEnergySourceTerms(Utils::Point3D(), ni, ni, ef, jf);
        const double cap =
            std::max(1e-30, params_.ion_heat_capacity * kBoltzmann * std::max(ni, 1.0));
        const double dEdt = src.ion_joule_heating + src.elastic_energy_exchange -
                            src.radiation_loss - conduction[i];

        energy_states_[i].ion_energy_density =
            std::max(0.0, energy_states_[i].ion_energy_density + dt * dEdt);
        energy_states_[i].ion_temperature =
            std::max(1.0, energy_states_[i].ion_energy_density / cap);
        energy_states_[i].elastic_cooling_rate = src.elastic_energy_exchange;
    }

    return VoidResult::success();
}

std::vector<double>
EnergyBalanceSolver::calculateThermalConduction(const std::vector<double>& temperature_field,
                                                double thermal_conductivity)
{
    std::vector<double> out(total_points_, 0.0);
    if (temperature_field.size() != total_points_)
    {
        return out;
    }

    for (size_t idx = 0; idx < total_points_; ++idx)
    {
        const double lap = calculateLaplacian(idx, temperature_field);
        out[idx] = -thermal_conductivity * lap;
    }

    return out;
}

VoidResult EnergyBalanceSolver::timeStep(const std::vector<double>& electron_density,
                                         const std::vector<double>& ion_density,
                                         const std::vector<Utils::Vector3D>& electric_field,
                                         const std::vector<Utils::Vector3D>& current_density,
                                         double dt)
{
    if (!initialized_)
    {
        return VoidResult::failure("solver not initialized");
    }

    energy_states_old_ = energy_states_;

    auto er = solveElectronEnergyBalance(electron_density, electric_field, current_density, dt);
    if (!er.isSuccess())
    {
        return er;
    }

    auto ir = solveIonEnergyBalance(ion_density, electric_field, current_density, dt);
    if (!ir.isSuccess())
    {
        return ir;
    }

    applyBoundaryConditions();
    updateStatistics();
    checkEnergyConservation();
    return VoidResult::success();
}

EnergyState EnergyBalanceSolver::getEnergyState(const Utils::Point3D& position) const
{
    if (!initialized_ || energy_states_.empty())
    {
        return EnergyState{};
    }

    const size_t i = static_cast<size_t>(
        std::clamp(position.x() / std::max(1e-30, spacing_.x()), 0.0, resolution_.x() - 1.0));
    const size_t j = static_cast<size_t>(
        std::clamp(position.y() / std::max(1e-30, spacing_.y()), 0.0, resolution_.y() - 1.0));
    const size_t k = static_cast<size_t>(
        std::clamp(position.z() / std::max(1e-30, spacing_.z()), 0.0, resolution_.z() - 1.0));

    return energy_states_[getGridIndex(i, j, k)];
}

const std::vector<EnergyState>& EnergyBalanceSolver::getAllEnergyStates() const
{
    return energy_states_;
}

std::vector<double> EnergyBalanceSolver::getElectronTemperatureDistribution() const
{
    std::vector<double> out(energy_states_.size(), 0.0);
    for (size_t i = 0; i < energy_states_.size(); ++i)
    {
        out[i] = energy_states_[i].electron_temperature;
    }
    return out;
}

std::vector<double> EnergyBalanceSolver::getIonTemperatureDistribution() const
{
    std::vector<double> out(energy_states_.size(), 0.0);
    for (size_t i = 0; i < energy_states_.size(); ++i)
    {
        out[i] = energy_states_[i].ion_temperature;
    }
    return out;
}

std::vector<double> EnergyBalanceSolver::getEnergyDensityDistribution() const
{
    std::vector<double> out(energy_states_.size(), 0.0);
    for (size_t i = 0; i < energy_states_.size(); ++i)
    {
        out[i] = energy_states_[i].electron_energy_density + energy_states_[i].ion_energy_density;
    }
    return out;
}

std::string EnergyBalanceSolver::getStatistics() const
{
    std::ostringstream oss;
    oss << "points=" << statistics_.grid_points << ", Te_avg=" << statistics_.average_electron_temp
        << ", Ti_avg=" << statistics_.average_ion_temp
        << ", Eerr=" << statistics_.energy_balance_error;
    return oss.str();
}

void EnergyBalanceSolver::reset()
{
    initialized_ = false;
    total_points_ = 1;
    energy_states_.clear();
    energy_states_old_.clear();
    statistics_ = EnergyBalanceStatistics{};
}

size_t EnergyBalanceSolver::getGridIndex(size_t i, size_t j, size_t k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    return k * nx * ny + j * nx + i;
}

Utils::Point3D EnergyBalanceSolver::getGridPosition(size_t i, size_t j, size_t k) const
{
    return Utils::Point3D((static_cast<double>(i) + 0.5) * spacing_.x(),
                          (static_cast<double>(j) + 0.5) * spacing_.y(),
                          (static_cast<double>(k) + 0.5) * spacing_.z());
}

void EnergyBalanceSolver::getGridIndices(size_t index, size_t& i, size_t& j, size_t& k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));

    k = index / (nx * ny);
    const size_t rem = index % (nx * ny);
    j = rem / nx;
    i = rem % nx;
}

double EnergyBalanceSolver::calculateLaplacian(size_t index, const std::vector<double>& field) const
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

Utils::Vector3D EnergyBalanceSolver::calculateGradient(size_t index,
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

double EnergyBalanceSolver::calculateDivergence(size_t index,
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

double EnergyBalanceSolver::calculateJouleHeating(double current_density_magnitude,
                                                  double conductivity)
{
    return (current_density_magnitude * current_density_magnitude) / std::max(conductivity, 1e-30);
}

double EnergyBalanceSolver::calculateElasticCooling(double electron_temp, double ion_temp,
                                                    double density)
{
    return params_.elastic_collision_frequency * density * kBoltzmann *
           std::max(0.0, electron_temp - ion_temp);
}

double EnergyBalanceSolver::calculateInelasticCooling(double electron_temp, double density)
{
    return params_.inelastic_collision_frequency * density * kBoltzmann *
           std::max(0.0, electron_temp - 300.0);
}

double EnergyBalanceSolver::calculateRadiationLoss(double electron_temp, double ion_temp,
                                                   double density)
{
    const double t = 0.5 * (electron_temp + ion_temp);
    return params_.radiation_coefficient * density * density * std::pow(std::max(1.0, t), 0.5);
}

double EnergyBalanceSolver::calculateEnergyExchange(double electron_temp, double ion_temp,
                                                    double density)
{
    return params_.energy_exchange_coefficient * density * (electron_temp - ion_temp);
}

void EnergyBalanceSolver::applyBoundaryConditions()
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
                EnergyState& st = energy_states_[getGridIndex(i, j, k)];
                st.electron_temperature = std::max(1.0, st.electron_temperature);
                st.ion_temperature = std::max(1.0, st.ion_temperature);
            }
        }
    }
}

bool EnergyBalanceSolver::isBoundaryPoint(size_t i, size_t j, size_t k) const
{
    const size_t nx = static_cast<size_t>(std::max(1.0, resolution_.x()));
    const size_t ny = static_cast<size_t>(std::max(1.0, resolution_.y()));
    const size_t nz = static_cast<size_t>(std::max(1.0, resolution_.z()));

    return i == 0 || j == 0 || k == 0 || i + 1 == nx || j + 1 == ny || k + 1 == nz;
}

void EnergyBalanceSolver::updateStatistics()
{
    statistics_ = EnergyBalanceStatistics{};
    statistics_.grid_points = energy_states_.size();

    if (energy_states_.empty())
    {
        return;
    }

    for (const auto& st : energy_states_)
    {
        statistics_.total_electron_energy += st.electron_energy_density;
        statistics_.total_ion_energy += st.ion_energy_density;
        statistics_.average_electron_temp += st.electron_temperature;
        statistics_.average_ion_temp += st.ion_temperature;
        statistics_.peak_electron_temp =
            std::max(statistics_.peak_electron_temp, st.electron_temperature);
        statistics_.peak_ion_temp = std::max(statistics_.peak_ion_temp, st.ion_temperature);
        statistics_.total_joule_heating += st.joule_heating_rate;
        statistics_.total_radiation_loss += st.radiation_loss_rate;
    }

    const double n = static_cast<double>(energy_states_.size());
    statistics_.average_electron_temp /= n;
    statistics_.average_ion_temp /= n;
}

double EnergyBalanceSolver::interpolateField(const Utils::Point3D& position,
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

void EnergyBalanceSolver::checkEnergyConservation()
{
    if (energy_states_old_.size() != energy_states_.size() || energy_states_.empty())
    {
        statistics_.energy_balance_error = 0.0;
        return;
    }

    double old_total = 0.0;
    double new_total = 0.0;
    for (size_t i = 0; i < energy_states_.size(); ++i)
    {
        old_total += energy_states_old_[i].electron_energy_density +
                     energy_states_old_[i].ion_energy_density;
        new_total +=
            energy_states_[i].electron_energy_density + energy_states_[i].ion_energy_density;
    }

    if (std::abs(old_total) < 1e-30)
    {
        statistics_.energy_balance_error = 0.0;
    }
    else
    {
        statistics_.energy_balance_error =
            std::abs(new_total - old_total) / std::abs(old_total) * 100.0;
    }
}

} // namespace FieldSolver
} // namespace SCDAT
