#include "../include/BoltzmannSolver.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace SCDAT
{
namespace FieldSolver
{

VoidResult BoltzmannSolver::setParameters(const BoltzmannParameters& parameters)
{
    if (parameters.energy_bins < 4 || parameters.max_energy_ev <= 0.0)
    {
        return VoidResult::failure("BoltzmannSolver parameters are invalid");
    }
    parameters_ = parameters;
    return VoidResult::success();
}

const BoltzmannState& BoltzmannSolver::solve(double electron_temperature_ev)
{
    const double effective_temperature = std::max(1.0e-3, electron_temperature_ev);
    const double dE = parameters_.max_energy_ev /
                      static_cast<double>(std::max<std::size_t>(1, parameters_.energy_bins - 1));

    state_.energy_grid_ev.assign(parameters_.energy_bins, 0.0);
    state_.distribution.assign(parameters_.energy_bins, 0.0);
    for (std::size_t i = 0; i < parameters_.energy_bins; ++i)
    {
        const double energy = static_cast<double>(i) * dE;
        state_.energy_grid_ev[i] = energy;
        state_.distribution[i] =
            std::sqrt(energy + dE) * std::exp(-energy / effective_temperature);
    }

    const double norm = std::accumulate(state_.distribution.begin(), state_.distribution.end(), 0.0);
    if (norm > 0.0)
    {
        for (double& value : state_.distribution)
        {
            value /= norm;
        }
    }

    state_.mean_energy_ev = 0.0;
    for (std::size_t i = 0; i < parameters_.energy_bins; ++i)
    {
        state_.mean_energy_ev += state_.energy_grid_ev[i] * state_.distribution[i];
    }

    const double collision_scale = 1.0 + parameters_.collision_frequency_hz * 1.0e-9;
    state_.mobility = parameters_.reduced_field_td / (collision_scale * 100.0);
    state_.ionization_rate = computeRateCoefficient(15.76);
    state_.excitation_rate = computeRateCoefficient(11.5);
    return state_;
}

double BoltzmannSolver::computeRateCoefficient(double threshold_ev) const
{
    double coefficient = 0.0;
    for (std::size_t i = 0; i < state_.energy_grid_ev.size(); ++i)
    {
        if (state_.energy_grid_ev[i] >= threshold_ev)
        {
            coefficient += std::sqrt(state_.energy_grid_ev[i]) * state_.distribution[i];
        }
    }

    return coefficient * parameters_.neutral_density_m3 * 1.0e-15;
}

} // namespace FieldSolver
} // namespace SCDAT
