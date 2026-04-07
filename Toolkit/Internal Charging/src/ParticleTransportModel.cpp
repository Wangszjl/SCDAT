#include "ParticleTransportModel.h"

#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

namespace
{
constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kEvToJ = 1.602176634e-19;
constexpr double kMinLayerThicknessM = 1.0e-12;
}

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
                                         double incident_charge_state_abs, double dt)
{
    if (layers_.empty() || dt <= 0.0 || current_density_a_per_m2 < 0.0 || mean_energy_ev < 0.0 ||
        incident_charge_state_abs <= 0.0 ||
        !std::isfinite(current_density_a_per_m2) || !std::isfinite(mean_energy_ev) ||
        !std::isfinite(incident_charge_state_abs) || !std::isfinite(dt))
    {
        return;
    }

    const double layer_thickness_m = std::max(layer_thickness_m_, kMinLayerThicknessM);
    const double attenuation_length = std::max(layer_thickness_m, thickness_m_ * 0.25);
    double weight_sum = 0.0;
    for (const auto& layer : layers_)
    {
        weight_sum += std::exp(-layer.depth_m / attenuation_length);
    }
    if (weight_sum <= 0.0)
    {
        return;
    }

    const double incident_particle_charge_c =
        std::max(incident_charge_state_abs * kElementaryChargeC, kElementaryChargeC);
    const double mean_particle_energy_j = mean_energy_ev * kEvToJ;

    for (auto& layer : layers_)
    {
        const double weight = std::exp(-layer.depth_m / attenuation_length);
        const double volumetric_charge =
            current_density_a_per_m2 * dt * weight / (weight_sum * layer_thickness_m);
        const double particles_per_m3 = volumetric_charge / incident_particle_charge_c;
        const double deposited_energy_j_per_m3 = particles_per_m3 * mean_particle_energy_j;

        layer.charge_density_c_per_m3 += volumetric_charge;
        layer.deposited_energy_j_per_m3 += deposited_energy_j_per_m3;
    }
}

void ParticleTransportModel::reset()
{
    for (auto& layer : layers_)
    {
        layer.charge_density_c_per_m3 = 0.0;
        layer.deposited_energy_j_per_m3 = 0.0;
        layer.dose_gy = 0.0;
        layer.effective_conductivity_s_per_m = 0.0;
        layer.deposited_charge_rate_c_per_m3_s = 0.0;
        layer.deposited_particle_rate_per_m3_s = 0.0;
    }
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
