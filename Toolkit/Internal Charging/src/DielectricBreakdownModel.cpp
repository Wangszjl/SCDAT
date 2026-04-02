#include "DielectricBreakdownModel.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

BreakdownAssessment DielectricBreakdownModel::assess(double max_field_v_per_m,
                                                     double breakdown_field_v_per_m,
                                                     double time_s) const
{
    BreakdownAssessment assessment;
    const double field_ratio =
        max_field_v_per_m / std::max(1.0e3, breakdown_field_v_per_m);
    assessment.probability =
        std::clamp(1.0 - std::exp(-std::max(0.0, field_ratio - 0.8) * (1.0 + time_s * 1.0e6)),
                   0.0, 1.0);
    assessment.triggered = assessment.probability >= 0.5;
    return assessment;
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
