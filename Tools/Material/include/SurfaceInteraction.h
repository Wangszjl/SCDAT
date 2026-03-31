#pragma once

#include "MaterialProperty.h"

namespace SCDAT
{
namespace Material
{

struct SurfaceImpactState
{
    double incident_energy_ev = 0.0;
    double incident_angle_rad = 0.0;
    double particle_charge_coulomb = 0.0;
    double surface_temperature_k = 300.0;
};

struct SurfaceInteractionResult
{
    double absorbed_charge_coulomb = 0.0;
    double reflected_energy_ev = 0.0;
    double secondary_emitted_electrons = 0.0;
    double deposited_heat_j_per_m2 = 0.0;
};

class SurfaceInteraction
{
  public:
    void setAbsorptionProbability(double value) { absorption_probability_ = value; }
    void setReflectionCoefficient(double value) { reflection_coefficient_ = value; }
    void setEmissionScaling(double value) { emission_scaling_ = value; }

    SurfaceInteractionResult evaluate(const MaterialProperty& material,
                                      const SurfaceImpactState& impact) const;

  private:
    double absorption_probability_ = 0.9;
    double reflection_coefficient_ = 0.05;
    double emission_scaling_ = 1.0;
};

} // namespace Material
} // namespace SCDAT
