#include "../include/SurfaceFieldVolumeBridge.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace FieldSolver
{
namespace
{

double finiteOr(double value, double fallback)
{
    return std::isfinite(value) ? value : fallback;
}

} // namespace

double blendFieldVolumeScalar(double base_value, double target_value, double blend_factor)
{
    const double base = finiteOr(base_value, 0.0);
    const double target = finiteOr(target_value, base);
    const double blend = std::isfinite(blend_factor)
                             ? std::clamp(blend_factor, 0.0, 1.0)
                             : 0.0;
    return (1.0 - blend) * base + blend * target;
}

double updateFieldVolumeMismatchMetric(double current_metric, double candidate_value,
                                       double baseline_value, double tolerance)
{
    double metric = finiteOr(current_metric, 0.0);
    if (!std::isfinite(candidate_value) || !std::isfinite(baseline_value))
    {
        return metric;
    }

    const double safe_tolerance = std::max(1.0e-12, std::abs(finiteOr(tolerance, 0.0)));
    const double normalized = std::abs(candidate_value - baseline_value) / safe_tolerance;
    if (!std::isfinite(normalized))
    {
        return metric;
    }
    return std::max(metric, normalized);
}

double computeFieldVolumeExternalBlendFactor(double mismatch_metric,
                                             double linear_residual_norm,
                                             double max_delta_v,
                                             bool linear_converged,
                                             double potential_tolerance_v,
                                             double volume_linear_relaxation)
{
    const double safe_tolerance = std::max(1.0e-9, std::abs(finiteOr(potential_tolerance_v, 0.0)));
    const double safe_mismatch = std::max(0.0, finiteOr(mismatch_metric, 0.0));
    const double safe_residual = std::abs(finiteOr(linear_residual_norm, 0.0));
    const double safe_delta = std::abs(finiteOr(max_delta_v, 0.0));

    const double residual_gate = 1.0 / (1.0 + safe_residual / safe_tolerance);
    const double delta_gate = 1.0 / (1.0 + safe_delta / safe_tolerance);
    const double convergence_gate = linear_converged ? 1.0 : 0.4;
    const double internal_confidence = convergence_gate * 0.5 * (residual_gate + delta_gate);
    const double stability_gate = 1.0 / (1.0 + safe_mismatch);
    const double external_trust =
        std::clamp(stability_gate + 0.5 * (1.0 - internal_confidence), 0.05, 1.0);
    const double relaxed_weight =
        std::clamp(finiteOr(volume_linear_relaxation, 0.2), 0.2, 1.0) * external_trust;
    return std::clamp(relaxed_weight, 0.05, 1.0);
}

} // namespace FieldSolver
} // namespace SCDAT
