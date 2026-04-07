#pragma once

#include "PlasmaChannelModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

enum class ArcIntegratorMode
{
  ExplicitGrowth,
  BoundedRelaxation,
};

struct ArcIntegratorParameters
{
  ArcIntegratorMode mode = ArcIntegratorMode::ExplicitGrowth;
  double conductivity_growth_per_ns = 0.1;
  double current_density_growth_per_ns = 0.05;
  double resistance_relaxation_per_ns = 0.05;
  double min_conductivity_s_per_m = 1.0;
  double max_conductivity_s_per_m = 1.0e7;
  double min_current_density_a_per_m2 = 0.0;
  double max_current_density_a_per_m2 = 1.0e8;
  double min_resistance_ohm = 1.0e-9;
  double max_resistance_ohm = 1.0e12;
};

class ArcPICIntegrator
{
  public:
  void setParameters(const ArcIntegratorParameters& parameters)
  {
    parameters_ = parameters;
  }

  const ArcIntegratorParameters& getParameters() const
  {
    return parameters_;
  }

    ArcChannelState advance(const ArcChannelState& channel, double dt) const;

  private:
  ArcIntegratorParameters parameters_{};
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
