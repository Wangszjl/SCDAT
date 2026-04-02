#include "ParticleTransportModel.h"

#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

bool ParticleTransportModel::initialize(double thickness_m, std::size_t layers)
{
    if (thickness_m <= 0.0 || layers == 0)
    {
        return false;
    }

    thickness_m_ = thickness_m;
    layer_thickness_m_ = thickness_m / static_cast<double>(layers);
    layers_.assign(layers, {});
    for (std::size_t i = 0; i < layers; ++i)
    {
        layers_[i].depth_m = (static_cast<double>(i) + 0.5) * layer_thickness_m_;
    }
    return true;
}

void ParticleTransportModel::depositFlux(double current_density_a_per_m2, double mean_energy_ev,
                                         double dt)
{
    if (layers_.empty())
    {
        return;
    }

    const double attenuation_length = std::max(layer_thickness_m_, thickness_m_ * 0.25);
    for (auto& layer : layers_)
    {
        const double weight = std::exp(-layer.depth_m / attenuation_length);
        const double volumetric_charge =
            current_density_a_per_m2 * dt * weight / std::max(layer_thickness_m_, 1.0e-12);
        layer.charge_density_c_per_m3 += volumetric_charge;
        layer.deposited_energy_j_per_m3 +=
            volumetric_charge * mean_energy_ev * 1.602176634e-19 / 1.602176634e-19;
    }
}

void ParticleTransportModel::reset()
{
    for (auto& layer : layers_)
    {
        layer.charge_density_c_per_m3 = 0.0;
        layer.deposited_energy_j_per_m3 = 0.0;
    }
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
