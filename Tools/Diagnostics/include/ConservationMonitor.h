#pragma once

#include "Monitor.h"

#include <functional>

namespace SCDAT
{
namespace Diagnostics
{

class ConservationMonitor : public Monitor
{
  public:
    using ScalarProvider = std::function<double()>;

    ConservationMonitor(std::string name, ScalarProvider charge_provider,
                        ScalarProvider energy_provider);

    VoidResult sample(double time) override;

  private:
    ScalarProvider charge_provider_;
    ScalarProvider energy_provider_;
    bool baseline_initialized_ = false;
    double initial_charge_ = 0.0;
    double initial_energy_ = 0.0;
};

} // namespace Diagnostics
} // namespace SCDAT
