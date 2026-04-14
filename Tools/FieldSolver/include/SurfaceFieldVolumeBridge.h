#pragma once

namespace SCDAT
{
namespace FieldSolver
{

double blendFieldVolumeScalar(double base_value, double target_value, double blend_factor);

double updateFieldVolumeMismatchMetric(double current_metric, double candidate_value,
                                       double baseline_value, double tolerance);

double computeFieldVolumeExternalBlendFactor(double mismatch_metric,
                                             double linear_residual_norm,
                                             double max_delta_v,
                                             bool linear_converged,
                                             double potential_tolerance_v,
                                             double volume_linear_relaxation);

} // namespace FieldSolver
} // namespace SCDAT
