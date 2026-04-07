#include "../include/CollisionPicAdapter.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Collision
{

namespace
{

bool isElectronLike(const ParticleObject& particle)
{
    switch (particle.getType())
    {
    case ParticleType::ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::BACKSCATTERED_ELECTRON:
    case ParticleType::AUGER_ELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
    case ParticleType::BETA_PARTICLE:
        return true;
    default:
        return particle.getCharge() < 0.0;
    }
}

bool isIonLike(const ParticleObject& particle)
{
    switch (particle.getType())
    {
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
    case ParticleType::PROTON:
    case ParticleType::ALPHA:
    case ParticleType::HEAVY_ION:
        return true;
    default:
        return particle.getCharge() > 0.0;
    }
}

double sanitizeDensity(double density)
{
    if (!std::isfinite(density))
    {
        return 0.0;
    }
    return std::max(0.0, density);
}

std::string toLower(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return text;
}

void accumulateReactionCount(const std::string& process_name,
                             std::size_t delta,
                             CollisionStepResult& result)
{
    const std::string normalized = toLower(process_name);
    if (normalized.find("ionization") != std::string::npos)
    {
        result.ionization_events += delta;
    }
    if (normalized.find("excitation") != std::string::npos)
    {
        result.excitation_events += delta;
    }
    if (normalized.find("charge_exchange") != std::string::npos ||
        normalized.find("charge-exchange") != std::string::npos ||
        normalized.find("chargeexchange") != std::string::npos)
    {
        result.charge_exchange_events += delta;
    }
}

} // namespace

CollisionBackgroundDensities
estimateCollisionBackgroundDensities(const Particle::ParticleContainer& container,
                                     double volume_m3,
                                     const CollisionBackgroundDensities& floor_densities)
{
    CollisionBackgroundDensities effective{};
    effective.neutral_density_m3 = sanitizeDensity(floor_densities.neutral_density_m3);
    effective.ion_density_m3 = sanitizeDensity(floor_densities.ion_density_m3);
    effective.electron_density_m3 = sanitizeDensity(floor_densities.electron_density_m3);

    if (!(volume_m3 > 0.0) || !std::isfinite(volume_m3))
    {
        return effective;
    }

    double electron_weight_count = 0.0;
    double ion_weight_count = 0.0;
    double neutral_weight_count = 0.0;

    for (const auto& particle : container)
    {
        if (!particle.isActive())
        {
            continue;
        }

        const double weight = std::max(0.0, particle.getWeight());
        if (isElectronLike(particle))
        {
            electron_weight_count += weight;
        }
        else if (isIonLike(particle))
        {
            ion_weight_count += weight;
        }
        else
        {
            neutral_weight_count += weight;
        }
    }

    const double inv_volume = 1.0 / volume_m3;
    effective.electron_density_m3 =
        std::max(effective.electron_density_m3, electron_weight_count * inv_volume);
    effective.ion_density_m3 = std::max(effective.ion_density_m3, ion_weight_count * inv_volume);
    effective.neutral_density_m3 =
        std::max(effective.neutral_density_m3, neutral_weight_count * inv_volume);

    return effective;
}

CollisionStepResult executeBackgroundMccStep(Particle::ParticleContainer& container,
                                             MonteCarloCollisionHandler& handler,
                                             double dt,
                                             double volume_m3,
                                             const CollisionBackgroundDensities& floor_densities,
                                             bool reinitialize_configuration)
{
    CollisionStepResult result{};

    if (!(dt > 0.0) || !std::isfinite(dt) || !(volume_m3 > 0.0) || !std::isfinite(volume_m3))
    {
        return result;
    }

    const auto effective_densities =
        estimateCollisionBackgroundDensities(container, volume_m3, floor_densities);
    result.effective_neutral_density_m3 = effective_densities.neutral_density_m3;
    result.effective_ion_density_m3 = effective_densities.ion_density_m3;
    result.effective_electron_density_m3 = effective_densities.electron_density_m3;

    if (reinitialize_configuration)
    {
        handler.initializeDefaultConfiguration(effective_densities.neutral_density_m3,
                                               effective_densities.ion_density_m3,
                                               effective_densities.electron_density_m3);
    }

    std::vector<Particle::ParticleId> active_particle_ids;
    std::vector<ParticleObject> active_particles;
    active_particle_ids.reserve(container.size());
    active_particles.reserve(container.size());

    for (const auto& particle : container)
    {
        if (!particle.isActive())
        {
            continue;
        }
        active_particle_ids.push_back(particle.getId());
        active_particles.push_back(particle);
    }

    result.active_particles = active_particles.size();
    if (active_particles.empty())
    {
        return result;
    }

    const auto stats_before = handler.getStatistics();
    const std::size_t total_before = static_cast<std::size_t>(stats_before.total_collisions);
    const std::size_t created_before = static_cast<std::size_t>(stats_before.particles_created);

    auto generated_particles = handler.processCollisions(active_particles, dt);

    for (std::size_t i = 0; i < active_particle_ids.size(); ++i)
    {
        if (auto* destination = container.getParticle(active_particle_ids[i]); destination)
        {
            *destination = active_particles[i];
        }
    }

    if (!generated_particles.empty())
    {
        container.addParticles(generated_particles);
    }

    container.removeInactiveParticles();

    const auto stats_after = handler.getStatistics();
    const std::size_t total_after = static_cast<std::size_t>(stats_after.total_collisions);
    const std::size_t created_after = static_cast<std::size_t>(stats_after.particles_created);

    for (const auto& entry : stats_after.collision_counts)
    {
        const std::string& process_name = entry.first;
        const std::size_t after_count = static_cast<std::size_t>(std::max(0, entry.second));

        std::size_t before_count = 0;
        const auto before_it = stats_before.collision_counts.find(process_name);
        if (before_it != stats_before.collision_counts.end())
        {
            before_count = static_cast<std::size_t>(std::max(0, before_it->second));
        }

        const std::size_t delta = after_count >= before_count ? after_count - before_count : after_count;
        if (delta > 0)
        {
            accumulateReactionCount(process_name, delta, result);
        }
    }

    result.collision_events = total_after >= total_before ? total_after - total_before : total_after;
    result.generated_particles =
        created_after >= created_before ? created_after - created_before : created_after;
    return result;
}

} // namespace Collision
} // namespace SCDAT
