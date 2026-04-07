#pragma once

#include "CollisionAlgorithm.h"

#include "../../../Particle/include/ParticleManager.h"

#include <cstddef>
#include <string>

namespace SCDAT
{
namespace Collision
{

struct CollisionBackgroundDensities
{
    double neutral_density_m3 = 1.0e20;
    double ion_density_m3 = 1.0e16;
    double electron_density_m3 = 1.0e16;
};

struct CollisionStepResult
{
    std::size_t active_particles = 0;
    std::size_t collision_events = 0;
    std::size_t generated_particles = 0;
    std::size_t ionization_events = 0;
    std::size_t excitation_events = 0;
    std::size_t charge_exchange_events = 0;
    double effective_neutral_density_m3 = 0.0;
    double effective_ion_density_m3 = 0.0;
    double effective_electron_density_m3 = 0.0;
};

CollisionBackgroundDensities
estimateCollisionBackgroundDensities(const Particle::ParticleContainer& container,
                                     double volume_m3,
                                     const CollisionBackgroundDensities& floor_densities);

CollisionStepResult executeBackgroundMccStep(Particle::ParticleContainer& container,
                                             MonteCarloCollisionHandler& handler,
                                             double dt,
                                             double volume_m3,
                                             const CollisionBackgroundDensities& floor_densities,
                                             bool reinitialize_configuration);

} // namespace Collision
} // namespace SCDAT
