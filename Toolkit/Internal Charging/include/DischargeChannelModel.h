#pragma once

#include "ParticleTransportModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

struct DischargeChannelResult
{
    double released_energy_j_per_m2 = 0.0;
    double released_charge_c_per_m2 = 0.0;
};

class DischargeChannelModel
{
  public:
    DischargeChannelResult relieve(std::vector<VolumeLayerState>& layers, double release_fraction) const;
};

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
