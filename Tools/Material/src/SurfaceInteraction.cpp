#include "../include/SurfaceInteraction.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Material
{

SurfaceInteractionResult SurfaceInteraction::evaluate(const MaterialProperty& material,
                                                      const SurfaceImpactState& impact) const
{
    const double clamped_angle = std::clamp(impact.incident_angle_rad, 0.0, 1.5707963267948966);
    const double cosine_factor = std::max(0.05, std::cos(clamped_angle));
    const double reflection =
        std::clamp(reflection_coefficient_ * cosine_factor, 0.0, 0.95);
    const double absorption = std::clamp(absorption_probability_ * (1.0 - 0.5 * reflection), 0.0,
                                         1.0);

    const double peak_energy_ev =
        std::max(1.0, material.getScalarProperty("secondary_yield_peak_energy_ev", 300.0));
    const double normalized_energy = std::max(1.0e-6, impact.incident_energy_ev / peak_energy_ev);
    const double emitted_secondaries =
        std::max(0.0, material.getSecondaryElectronYield() * emission_scaling_ * normalized_energy *
                          std::exp(1.0 - normalized_energy) * cosine_factor);

    SurfaceInteractionResult result;
    result.absorbed_charge_coulomb = impact.particle_charge_coulomb * absorption;
    result.reflected_energy_ev = impact.incident_energy_ev * reflection;
    result.secondary_emitted_electrons = emitted_secondaries;
    result.deposited_heat_j_per_m2 = impact.incident_energy_ev * 1.602176634e-19 *
                                     (absorption + 0.25 * emitted_secondaries);
    return result;
}

} // namespace Material
} // namespace SCDAT
