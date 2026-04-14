#include "SurfaceSimulationRunner.h"

#include "DensePlasmaSurfaceCharging.h"
#include "SurfaceRuntimePlan.h"

#include <algorithm>
#include <stdexcept>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

SurfaceSimulationRunResult makeFailure(std::string message)
{
    SurfaceSimulationRunResult result;
    result.success = false;
    result.error_message = std::move(message);
    return result;
}

SurfaceSimulationRunResult runAdaptiveSurfaceSimulation(DensePlasmaSurfaceCharging& charging,
                                                        const SurfaceChargingScenarioPreset& preset)
{
    double elapsed_time_s = 0.0;
    std::size_t sample_count = 0;
    const double minimum_dt_s = std::max(1.0e-12, preset.minimum_time_step_s);
    const double maximum_dt_s = std::max(minimum_dt_s, preset.maximum_time_step_s);
    while (elapsed_time_s + 1.0e-12 < preset.total_duration_s)
    {
        const double remaining_time_s = preset.total_duration_s - elapsed_time_s;
        double dt = 0.0;
        if (preset.steps > 0 && sample_count >= preset.steps)
        {
            // Keep tail steps bounded once adaptive sample budget is exhausted.
            dt = std::min(remaining_time_s, maximum_dt_s);
        }
        else
        {
            dt = charging.recommendTimeStep(remaining_time_s, minimum_dt_s, maximum_dt_s);
        }
        if (!charging.advance(dt))
        {
            return makeFailure("Surface charging adaptive advance failed");
        }
        elapsed_time_s += dt;
        sample_count += 1;
    }
    return SurfaceSimulationRunResult{true, {}};
}

} // namespace

SurfaceSimulationRunResult SurfaceSimulationRunner::run(
    const SurfaceChargingScenarioPreset& preset,
    const std::filesystem::path& output_path) const
{
    DensePlasmaSurfaceCharging charging;
    const auto runtime_plan = compileSurfaceRuntimePlan(preset);
    if (!charging.initialize(runtime_plan))
    {
        std::string message = "Failed to initialize surface charging toolkit";
        if (!charging.lastErrorMessage().empty())
        {
            message += ": " + charging.lastErrorMessage();
        }
        return makeFailure(message);
    }

    if (preset.adaptive_time_stepping && preset.total_duration_s > 0.0)
    {
        const auto adaptive_result = runAdaptiveSurfaceSimulation(charging, preset);
        if (!adaptive_result.success)
        {
            return adaptive_result;
        }
    }
    else
    {
        for (std::size_t i = 0; i < preset.steps; ++i)
        {
            if (!charging.advance(preset.time_step_s))
            {
                return makeFailure("Surface charging advance failed");
            }
        }
    }

    if (!charging.exportResults(output_path))
    {
        return makeFailure("Failed to export surface charging results");
    }

    return SurfaceSimulationRunResult{true, {}};
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
