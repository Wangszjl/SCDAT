#pragma once

#include "Monitor.h"

#include "../../Particle/include/ParticleManager.h"

#include <functional>

namespace SCDAT
{
namespace Diagnostics
{

class ParticleMonitor : public Monitor
{
  public:
    using StatisticsProvider = std::function<Particle::ParticleManager::Statistics()>;

    ParticleMonitor(std::string name, StatisticsProvider statistics_provider);

    VoidResult sample(double time) override;

  private:
    StatisticsProvider statistics_provider_;
};

} // namespace Diagnostics
} // namespace SCDAT
