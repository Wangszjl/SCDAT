#pragma once

#include "CollisionAlgorithm.h"

#include "../../../Particle/include/ParticleManager.h"

#include <cstddef>
#include <optional>
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

enum class CollisionCrossSectionSetId
{
    BackgroundMccV1,
    SurfacePicV1,
    VacuumArcPicV1,
};

struct BackgroundMccRuntimeConfig
{
    CollisionBackgroundDensities floor_densities;
    CollisionCrossSectionSetId cross_section_set_id = CollisionCrossSectionSetId::BackgroundMccV1;
    bool reinitialize_configuration = false;
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

std::optional<CollisionCrossSectionSetId> parseCollisionCrossSectionSetId(const std::string& text);
std::string toString(CollisionCrossSectionSetId set_id);

void initializeConfiguredCollisionHandler(MonteCarloCollisionHandler& handler,
                                          CollisionCrossSectionSetId set_id,
                                          double neutral_density_m3,
                                          double ion_density_m3,
                                          double electron_density_m3);

CollisionStepResult executeBackgroundMccStep(Particle::ParticleContainer& container,
                                             MonteCarloCollisionHandler& handler,
                                             double dt,
                                             double volume_m3,
                                             const CollisionBackgroundDensities& floor_densities,
                                             bool reinitialize_configuration);

CollisionStepResult executeBackgroundMccStep(Particle::ParticleContainer& container,
                                             MonteCarloCollisionHandler& handler,
                                             double dt,
                                             double volume_m3,
                                             const BackgroundMccRuntimeConfig& runtime_config);

} // namespace Collision
} // namespace SCDAT


