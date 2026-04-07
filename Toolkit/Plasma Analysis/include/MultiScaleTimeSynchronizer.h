#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct MultiScaleSynchronizationThresholds
{
    double relative_charge_tolerance = 1.0e-3;
    double relative_energy_tolerance = 1.0e-3;
    double relative_field_tolerance = 1.0e-3;
    std::size_t max_iterations = 6;
    std::size_t min_iterations = 1;
};

struct MultiScaleSynchronizationConfig
{
    double macro_time_step_s = 1.0e-8;
    std::size_t transition_substeps = 2;
    std::size_t pic_substeps_per_transition = 2;
    MultiScaleSynchronizationThresholds thresholds;
};

struct MultiScaleExchangeState
{
    double total_charge_c = 0.0;
    double total_energy_j = 0.0;
    double field_l2_norm = 0.0;
};

struct MultiScaleSynchronizationResidual
{
    double relative_charge_error = 0.0;
    double relative_energy_error = 0.0;
    double relative_field_error = 0.0;
    bool converged = false;
};

struct MultiScaleSynchronizationStatus
{
    double simulated_time_s = 0.0;
    std::size_t macro_steps_completed = 0;
    std::size_t transition_substeps_executed = 0;
    std::size_t pic_substeps_executed = 0;
    std::size_t correction_iterations = 0;
    MultiScaleSynchronizationResidual last_residual{};
    std::string convergence_reason;
};

class MultiScaleTimeSynchronizer
{
  public:
    using AdvanceFunction = std::function<bool(double)>;
    using ExchangeSampler = std::function<MultiScaleExchangeState()>;
    using CorrectionFunction = std::function<void(const MultiScaleSynchronizationResidual&)>;

    bool initialize(const MultiScaleSynchronizationConfig& config);
    bool advanceOneMacroStep(const AdvanceFunction& advance_macro,
                             const AdvanceFunction& advance_transition,
                             const AdvanceFunction& advance_pic,
                             const ExchangeSampler& sample_exchange,
                             const CorrectionFunction& apply_correction = {});

    const MultiScaleSynchronizationConfig& getConfig() const { return config_; }
    const MultiScaleSynchronizationStatus& getStatus() const { return status_; }
    void reset();

  private:
    static double relativeError(double candidate, double reference);

    MultiScaleSynchronizationConfig config_{};
    MultiScaleSynchronizationStatus status_{};
    bool initialized_ = false;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
