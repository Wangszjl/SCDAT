#pragma once

#include "../../Basic/include/VoidResult.h"

#include <cstddef>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

struct BoltzmannParameters
{
    std::size_t energy_bins = 128;
    double max_energy_ev = 100.0;
    double reduced_field_td = 100.0;
    double neutral_density_m3 = 2.5e20;
    double collision_frequency_hz = 1.0e8;
};

struct BoltzmannState
{
    std::vector<double> energy_grid_ev;
    std::vector<double> distribution;
    double mean_energy_ev = 0.0;
    double mobility = 0.0;
    double ionization_rate = 0.0;
    double excitation_rate = 0.0;
};

class BoltzmannSolver
{
  public:
    VoidResult setParameters(const BoltzmannParameters& parameters);
    const BoltzmannParameters& getParameters() const { return parameters_; }

    const BoltzmannState& solve(double electron_temperature_ev);
    const BoltzmannState& getState() const { return state_; }
    double computeRateCoefficient(double threshold_ev) const;

  private:
    BoltzmannParameters parameters_;
    BoltzmannState state_;
};

} // namespace FieldSolver
} // namespace SCDAT
