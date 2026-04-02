#include "ArcDischargeDetector.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

bool ArcDischargeDetector::shouldTrigger(double effective_field_v_per_m, double avalanche_gain) const
{
    return effective_field_v_per_m >= 3.0e7 && avalanche_gain >= 10.0;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
