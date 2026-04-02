#include "DischargeChannelModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

DischargeChannelResult
DischargeChannelModel::relieve(std::vector<VolumeLayerState>& layers, double release_fraction) const
{
    DischargeChannelResult result;
    const double clamped_fraction = release_fraction < 0.0 ? 0.0 : release_fraction;
    for (auto& layer : layers)
    {
        const double released_charge = layer.charge_density_c_per_m3 * clamped_fraction;
        const double released_energy = layer.deposited_energy_j_per_m3 * clamped_fraction;
        layer.charge_density_c_per_m3 -= released_charge;
        layer.deposited_energy_j_per_m3 -= released_energy;
        result.released_charge_c_per_m2 += released_charge;
        result.released_energy_j_per_m2 += released_energy;
    }
    return result;
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
