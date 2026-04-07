#include "ArcPICIntegrator.h"

#include <algorithm>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

ArcChannelState ArcPICIntegrator::advance(const ArcChannelState& channel, double dt) const
{
    ArcChannelState updated = channel;

    const double dt_ns = std::max(0.0, dt * 1.0e9);
    if (parameters_.mode == ArcIntegratorMode::BoundedRelaxation)
    {
        const double conductivity_span =
            std::max(1.0e-18, parameters_.max_conductivity_s_per_m - parameters_.min_conductivity_s_per_m);
        const double current_span =
            std::max(1.0e-18, parameters_.max_current_density_a_per_m2 - parameters_.min_current_density_a_per_m2);
        const double resistance_span =
            std::max(1.0e-18, parameters_.max_resistance_ohm - parameters_.min_resistance_ohm);

        const double conductivity_ratio = std::clamp(
            (updated.conductivity_s_per_m - parameters_.min_conductivity_s_per_m) / conductivity_span,
            0.0, 1.0);
        const double current_ratio = std::clamp(
            (updated.current_density_a_per_m2 - parameters_.min_current_density_a_per_m2) / current_span,
            0.0, 1.0);
        const double resistance_ratio = std::clamp(
            (updated.resistance_ohm - parameters_.min_resistance_ohm) / resistance_span,
            0.0, 1.0);

        updated.conductivity_s_per_m += parameters_.conductivity_growth_per_ns * dt_ns *
                                        conductivity_span * (1.0 - conductivity_ratio);
        updated.current_density_a_per_m2 += parameters_.current_density_growth_per_ns * dt_ns *
                                            current_span * (1.0 - current_ratio);
        updated.resistance_ohm -= parameters_.resistance_relaxation_per_ns * dt_ns * resistance_span *
                                  resistance_ratio;
    }
    else
    {
        updated.conductivity_s_per_m *= 1.0 + parameters_.conductivity_growth_per_ns * dt_ns;
        updated.current_density_a_per_m2 *= 1.0 + parameters_.current_density_growth_per_ns * dt_ns;
        updated.resistance_ohm /= 1.0 + parameters_.resistance_relaxation_per_ns * dt_ns;
    }

    updated.conductivity_s_per_m =
        std::clamp(updated.conductivity_s_per_m, parameters_.min_conductivity_s_per_m,
                   parameters_.max_conductivity_s_per_m);
    updated.current_density_a_per_m2 =
        std::clamp(updated.current_density_a_per_m2, parameters_.min_current_density_a_per_m2,
                   parameters_.max_current_density_a_per_m2);
    updated.resistance_ohm = std::clamp(updated.resistance_ohm, parameters_.min_resistance_ohm,
                                        parameters_.max_resistance_ohm);

    return updated;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
