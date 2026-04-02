#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

struct BreakdownAssessment
{
    double probability = 0.0;
    bool triggered = false;
};

class DielectricBreakdownModel
{
  public:
    BreakdownAssessment assess(double max_field_v_per_m, double breakdown_field_v_per_m,
                               double time_s) const;
};

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
