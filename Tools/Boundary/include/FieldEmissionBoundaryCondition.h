#pragma once

#include "BoundaryConditions.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace SCDAT
{
namespace Particle
{

class FieldEmissionBoundaryCondition : public BoundaryConditionBase
{
  public:
    struct Parameters
    {
        double emission_current_density_a_per_m2 = 0.0;
        double emission_area_m2 = 1.0e-8;
        double emitted_electron_speed_m_per_s = 2.0e6;
        double emitted_electron_weight = 1.0;
        std::size_t max_emitted_particles_per_event = 64;

      // Ion-driven secondary electron emission approximation.
      bool enable_secondary_electron_emission = false;
      double secondary_electron_yield_per_ion = 0.0;
      double secondary_electron_speed_m_per_s = 8.0e5;
      double secondary_electron_weight = 1.0;
      std::size_t max_secondary_particles_per_event = 64;
    };

    FieldEmissionBoundaryCondition();
    explicit FieldEmissionBoundaryCondition(const Parameters& parameters);

    void setEmissionCurrentDensity(double emission_current_density_a_per_m2);
    double getEmissionCurrentDensity() const
    {
        return parameters_.emission_current_density_a_per_m2;
    }

    void setEmissionArea(double emission_area_m2);

    // Clears the fractional carry used to avoid low-current over-emission bias.
    void resetEmissionAccumulator();

    std::vector<ParticleTypeDef>
    processParticle(ParticleTypeDef& particle,
                    const SCDAT::Geometry::Point3D& intersection_point,
                    const SCDAT::Geometry::Vector3D& normal, double dt) override;

  private:
    static bool isIonDrivenEmissionSource(const ParticleTypeDef& particle);

    Parameters parameters_{};
    std::uint64_t next_particle_id_ = 1000000000ULL;
    double pending_emission_particles_ = 0.0;
    double pending_secondary_particles_ = 0.0;
};

} // namespace Particle
} // namespace SCDAT
