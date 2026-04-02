#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class ArcDischargeDetector
{
  public:
    bool shouldTrigger(double effective_field_v_per_m, double avalanche_gain) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
