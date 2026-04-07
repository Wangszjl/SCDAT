#include "MultiScaleTimeSynchronizer.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{

constexpr double kEpsilon = 1.0e-18;

} // namespace

bool MultiScaleTimeSynchronizer::initialize(const MultiScaleSynchronizationConfig& config)
{
    if (config.macro_time_step_s <= 0.0 || config.transition_substeps == 0 ||
        config.pic_substeps_per_transition == 0 || config.thresholds.max_iterations == 0)
    {
        return false;
    }

    config_ = config;
    status_ = MultiScaleSynchronizationStatus{};
    status_.convergence_reason.clear();
    initialized_ = true;
    return true;
}

bool MultiScaleTimeSynchronizer::advanceOneMacroStep(const AdvanceFunction& advance_macro,
                                                     const AdvanceFunction& advance_transition,
                                                     const AdvanceFunction& advance_pic,
                                                     const ExchangeSampler& sample_exchange,
                                                     const CorrectionFunction& apply_correction)
{
    if (!initialized_ || !advance_macro || !advance_transition || !advance_pic ||
        !sample_exchange)
    {
        return false;
    }

    const double transition_dt =
        config_.macro_time_step_s / static_cast<double>(config_.transition_substeps);
    const double pic_dt =
        transition_dt / static_cast<double>(config_.pic_substeps_per_transition);

    const MultiScaleExchangeState start_state = sample_exchange();
    MultiScaleExchangeState previous_state = start_state;

    bool converged = false;
    std::size_t iterations_used = 0;

    for (std::size_t iteration = 1; iteration <= config_.thresholds.max_iterations; ++iteration)
    {
        if (!advance_macro(config_.macro_time_step_s))
        {
            status_.convergence_reason = "macro-step-failed";
            return false;
        }

        for (std::size_t transition_step = 0; transition_step < config_.transition_substeps;
             ++transition_step)
        {
            if (!advance_transition(transition_dt))
            {
                status_.convergence_reason = "transition-step-failed";
                return false;
            }
            ++status_.transition_substeps_executed;

            for (std::size_t pic_step = 0;
                 pic_step < config_.pic_substeps_per_transition; ++pic_step)
            {
                if (!advance_pic(pic_dt))
                {
                    status_.convergence_reason = "pic-step-failed";
                    return false;
                }
                ++status_.pic_substeps_executed;
            }
        }

        const MultiScaleExchangeState current_state = sample_exchange();
        MultiScaleSynchronizationResidual residual;
        residual.relative_charge_error =
            relativeError(current_state.total_charge_c, start_state.total_charge_c);
        residual.relative_energy_error =
            relativeError(current_state.total_energy_j, start_state.total_energy_j);
        residual.relative_field_error =
            relativeError(current_state.field_l2_norm, previous_state.field_l2_norm);

        const bool meets_charge =
            residual.relative_charge_error <= config_.thresholds.relative_charge_tolerance;
        const bool meets_energy =
            residual.relative_energy_error <= config_.thresholds.relative_energy_tolerance;
        const bool meets_field =
            residual.relative_field_error <= config_.thresholds.relative_field_tolerance;

        residual.converged =
            (iteration >= config_.thresholds.min_iterations) && meets_charge && meets_energy &&
            meets_field;

        status_.last_residual = residual;
        previous_state = current_state;
        iterations_used = iteration;

        if (residual.converged)
        {
            converged = true;
            break;
        }

        if (apply_correction)
        {
            apply_correction(residual);
            ++status_.correction_iterations;
        }
    }

    status_.macro_steps_completed += 1;
    status_.simulated_time_s += config_.macro_time_step_s;

    if (converged)
    {
        status_.convergence_reason = "criteria-met";
        return true;
    }

    status_.convergence_reason =
        (iterations_used >= config_.thresholds.max_iterations) ? "max-iterations" : "not-converged";
    return false;
}

void MultiScaleTimeSynchronizer::reset()
{
    *this = MultiScaleTimeSynchronizer{};
}

double MultiScaleTimeSynchronizer::relativeError(double candidate, double reference)
{
    return std::abs(candidate - reference) / std::max(std::abs(reference), kEpsilon);
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
