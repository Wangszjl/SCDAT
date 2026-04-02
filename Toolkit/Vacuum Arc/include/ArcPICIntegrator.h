#pragma once

#include "PlasmaChannelModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class ArcPICIntegrator
{
  public:
    ArcChannelState advance(const ArcChannelState& channel, double dt) const;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
