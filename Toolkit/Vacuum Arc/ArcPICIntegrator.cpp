#include "ArcPICIntegrator.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

ArcChannelState ArcPICIntegrator::advance(const ArcChannelState& channel, double dt) const
{
    ArcChannelState updated = channel;
    updated.conductivity_s_per_m *= 1.0 + 0.1 * dt * 1.0e9;
    updated.current_density_a_per_m2 *= 1.0 + 0.05 * dt * 1.0e9;
    updated.resistance_ohm = updated.resistance_ohm / (1.0 + 0.05 * dt * 1.0e9);
    return updated;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
