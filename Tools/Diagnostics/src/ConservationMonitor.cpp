#include "../include/ConservationMonitor.h"

#include <algorithm>

namespace SCDAT
{
namespace Diagnostics
{

ConservationMonitor::ConservationMonitor(std::string name, ScalarProvider charge_provider,
                                         ScalarProvider energy_provider)
    : Monitor(std::move(name)),
      charge_provider_(std::move(charge_provider)),
      energy_provider_(std::move(energy_provider))
{
}

VoidResult ConservationMonitor::sample(double time)
{
    if (!charge_provider_ || !energy_provider_)
    {
        return VoidResult::failure("ConservationMonitor requires charge and energy providers");
    }

    const double charge = charge_provider_();
    const double energy = energy_provider_();
    if (!baseline_initialized_)
    {
        initial_charge_ = charge;
        initial_energy_ = energy;
        baseline_initialized_ = true;
    }

    const double charge_scale = std::max(1.0, std::abs(initial_charge_));
    const double energy_scale = std::max(1.0, std::abs(initial_energy_));

    recordSample(time,
                 {
                     {"charge_c", charge},
                     {"energy_j", energy},
                     {"charge_drift_c", charge - initial_charge_},
                     {"energy_drift_j", energy - initial_energy_},
                     {"charge_drift_rel", (charge - initial_charge_) / charge_scale},
                     {"energy_drift_rel", (energy - initial_energy_) / energy_scale},
                 });
    return VoidResult::success();
}

} // namespace Diagnostics
} // namespace SCDAT
