#include "../include/ParticleMonitor.h"

namespace SCDAT
{
namespace Diagnostics
{

ParticleMonitor::ParticleMonitor(std::string name, StatisticsProvider statistics_provider)
    : Monitor(std::move(name)), statistics_provider_(std::move(statistics_provider))
{
}

VoidResult ParticleMonitor::sample(double time)
{
    if (!statistics_provider_)
    {
        return VoidResult::failure("ParticleMonitor has no statistics provider");
    }

    const auto stats = statistics_provider_();
    recordSample(time,
                 {
                     {"particles_total", static_cast<double>(stats.total_particles)},
                     {"particles_active", static_cast<double>(stats.active_particles)},
                     {"kinetic_energy_j", stats.total_kinetic_energy},
                     {"total_charge_c", stats.total_charge},
                     {"average_age_s", stats.average_age},
                 });
    return VoidResult::success();
}

} // namespace Diagnostics
} // namespace SCDAT
