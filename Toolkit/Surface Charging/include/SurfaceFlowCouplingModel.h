#pragma once

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceFlowProjection
{
    double body_alignment_cosine = 1.0;
    double patch_alignment_cosine = 1.0;
    double body_projected_speed_m_per_s = 0.0;
    double patch_projected_speed_m_per_s = 0.0;
};

inline double normalizeFlowAlignmentCosine(double flow_alignment_cosine)
{
    if (!std::isfinite(flow_alignment_cosine))
    {
        return 0.0;
    }
    return std::clamp(flow_alignment_cosine, -1.0, 1.0);
}

inline double flowAlignmentCosineFromAngleDeg(double flow_angle_deg)
{
    if (!std::isfinite(flow_angle_deg))
    {
        return 1.0;
    }

    constexpr double kPi = 3.14159265358979323846;
    const double angle_rad = flow_angle_deg * kPi / 180.0;
    return normalizeFlowAlignmentCosine(std::cos(angle_rad));
}

inline bool hasExplicitPatchFlowAngle(double patch_flow_angle_deg)
{
    if (!std::isfinite(patch_flow_angle_deg))
    {
        return false;
    }

    const double wrapped_deg = std::remainder(patch_flow_angle_deg, 360.0);
    return std::abs(wrapped_deg) > 1.0e-6;
}

inline SurfaceFlowProjection
resolveSurfaceFlowProjection(double bulk_flow_velocity_m_per_s, double flow_alignment_cosine,
                             double patch_flow_angle_deg)
{
    const double bulk_speed = std::max(0.0, bulk_flow_velocity_m_per_s);

    SurfaceFlowProjection projection;
    projection.body_alignment_cosine = normalizeFlowAlignmentCosine(flow_alignment_cosine);
    projection.patch_alignment_cosine = hasExplicitPatchFlowAngle(patch_flow_angle_deg)
                                          ? flowAlignmentCosineFromAngleDeg(patch_flow_angle_deg)
                                          : projection.body_alignment_cosine;

    projection.body_projected_speed_m_per_s =
        bulk_speed * projection.body_alignment_cosine;
    projection.patch_projected_speed_m_per_s =
        bulk_speed * projection.patch_alignment_cosine;
    return projection;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
