#include "../include/CollisionPicAdapter.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Collision
{

namespace
{

using CollisionParameters = SCDAT::Collision::CollisionParameters;
using CollisionProcessFactory = SCDAT::Collision::CollisionProcessFactory;
using ParticleSpecies = SCDAT::Collision::ParticleSpecies;

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

CollisionParameters makeElasticParams(double cross_section_max)
{
    CollisionParameters params;
    params.cross_section_max = cross_section_max;
    return params;
}

CollisionParameters makeExcitationParams(double threshold_ev,
                                         double cross_section_max,
                                         double energy_loss_ev)
{
    CollisionParameters params;
    params.energy_threshold = threshold_ev;
    params.cross_section_max = cross_section_max;
    params.energy_loss = energy_loss_ev;
    return params;
}

CollisionParameters makeIonizationParams(double threshold_ev,
                                         double cross_section_max,
                                         double energy_loss_ev,
                                         double ion_mass_amu,
                                         double ion_charge_state)
{
    CollisionParameters params;
    params.energy_threshold = threshold_ev;
    params.cross_section_max = cross_section_max;
    params.energy_loss = energy_loss_ev;
    params.additional_params["ion_mass_number"] = ion_mass_amu;
    params.additional_params["ion_charge_number"] = ion_charge_state;
    return params;
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

void initializeSurfacePicCollisionHandler(MonteCarloCollisionHandler& handler,
                                              double neutral_density_m3)
{
    handler.clearCollisionProcesses();
    handler.addCollisionProcess(
        "elastic_e_n",
        CollisionProcessFactory::createElasticProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            makeElasticParams(2.0e-19)));
    handler.addCollisionProcess(
        "excitation_e_n",
        CollisionProcessFactory::createExcitationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            makeExcitationParams(11.5, 2.5e-20, 11.5)));
    handler.addCollisionProcess(
        "ionization_e_n",
        CollisionProcessFactory::createIonizationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            makeIonizationParams(15.76, 3.0e-20, 15.76, 40.0, 1.0)));
    handler.resetStatistics();
}

void initializeVacuumArcCollisionHandler(MonteCarloCollisionHandler& handler,
                                         double neutral_density_m3,
                                         double ion_density_m3,
                                         double electron_density_m3)
{
    handler.clearCollisionProcesses();
    handler.initializeDefaultConfiguration(neutral_density_m3, ion_density_m3, electron_density_m3);

    auto excitation_params = handler.getDefaultCollisionParameters(CollisionType::EXCITATION);
    excitation_params.energy_threshold = 10.0;
    excitation_params.cross_section_max = 3.2e-20;
    excitation_params.energy_loss = 10.0;
    handler.setDefaultCollisionParameters(CollisionType::EXCITATION, excitation_params);

    auto charge_exchange_params =
        handler.getDefaultCollisionParameters(CollisionType::CHARGE_EXCHANGE);
    charge_exchange_params.cross_section_max = 1.8e-19;
    charge_exchange_params.energy_loss = 0.35;
    handler.setDefaultCollisionParameters(CollisionType::CHARGE_EXCHANGE, charge_exchange_params);

    handler.clearCollisionProcesses();
    handler.addCollisionProcess(
        "elastic_e_n",
        CollisionProcessFactory::createElasticProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::ELASTIC)));
    handler.addCollisionProcess(
        "inelastic_e_n",
        CollisionProcessFactory::createInelasticProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::INELASTIC)));
    handler.addCollisionProcess(
        "ionization_e_n",
        CollisionProcessFactory::createIonizationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::IONIZATION)));
    handler.addCollisionProcess(
        "excitation_e_n",
        CollisionProcessFactory::createExcitationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::EXCITATION)));
    handler.addCollisionProcess(
        "charge_exchange_i_n",
        CollisionProcessFactory::createChargeExchangeProcess(
            ParticleSpecies::ION, ParticleSpecies::NEUTRAL, neutral_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::CHARGE_EXCHANGE)));
    handler.addCollisionProcess(
        "recombination_e_i",
        CollisionProcessFactory::createRecombinationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::ION, ion_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::RECOMBINATION)));
    handler.addCollisionProcess(
        "coulomb_e_i",
        CollisionProcessFactory::createCoulombProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::ION, electron_density_m3,
            handler.getDefaultCollisionParameters(CollisionType::COULOMB)));
    handler.resetStatistics();
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

std::optional<CollisionCrossSectionSetId> parseCollisionCrossSectionSetId(const std::string& text)
{
    const std::string normalized = toLower(text);
    if (normalized == "background_mcc_v1" || normalized == "backgroundmccv1" ||
        normalized == "default" || normalized == "background")
    {
        return CollisionCrossSectionSetId::BackgroundMccV1;
    }
    if (normalized == "surface_pic_v1" || normalized == "surfacepicv1" ||
        normalized == "surface_pic" || normalized == "surfacepic")
    {
        return CollisionCrossSectionSetId::SurfacePicV1;
    }
    if (normalized == "vacuum_arc_pic_v1" || normalized == "vacuumarcpicv1" ||
        normalized == "arcpic_background_v1" || normalized == "arcpic")
    {
        return CollisionCrossSectionSetId::VacuumArcPicV1;
    }
    return std::nullopt;
}

std::string toString(CollisionCrossSectionSetId set_id)
{
    switch (set_id)
    {
    case CollisionCrossSectionSetId::SurfacePicV1:
        return "surface_pic_v1";
    case CollisionCrossSectionSetId::VacuumArcPicV1:
        return "vacuum_arc_pic_v1";
    case CollisionCrossSectionSetId::BackgroundMccV1:
    default:
        return "background_mcc_v1";
    }
}

void initializeConfiguredCollisionHandler(MonteCarloCollisionHandler& handler,
                                          CollisionCrossSectionSetId set_id,
                                          double neutral_density_m3,
                                          double ion_density_m3,
                                          double electron_density_m3)
{
    const double safe_neutral_density = sanitizeDensity(neutral_density_m3);
    const double safe_ion_density = sanitizeDensity(ion_density_m3);
    const double safe_electron_density = sanitizeDensity(electron_density_m3);

    switch (set_id)
    {
    case CollisionCrossSectionSetId::SurfacePicV1:
        initializeSurfacePicCollisionHandler(handler, safe_neutral_density);
        break;
    case CollisionCrossSectionSetId::VacuumArcPicV1:
        initializeVacuumArcCollisionHandler(handler, safe_neutral_density, safe_ion_density,
                                            safe_electron_density);
        break;
    case CollisionCrossSectionSetId::BackgroundMccV1:
    default:
        handler.initializeDefaultConfiguration(safe_neutral_density, safe_ion_density,
                                               safe_electron_density);
        break;
    }
}

CollisionStepResult executeBackgroundMccStep(Particle::ParticleContainer& container,
                                             MonteCarloCollisionHandler& handler,
                                             double dt,
                                             double volume_m3,
                                             const CollisionBackgroundDensities& floor_densities,
                                             bool reinitialize_configuration)
{
    BackgroundMccRuntimeConfig runtime_config;
    runtime_config.floor_densities = floor_densities;
    runtime_config.reinitialize_configuration = reinitialize_configuration;
    return executeBackgroundMccStep(container, handler, dt, volume_m3, runtime_config);
}

CollisionStepResult executeBackgroundMccStep(Particle::ParticleContainer& container,
                                             MonteCarloCollisionHandler& handler,
                                             double dt,
                                             double volume_m3,
                                             const BackgroundMccRuntimeConfig& runtime_config)
{
    CollisionStepResult result{};

    if (!(dt > 0.0) || !std::isfinite(dt) || !(volume_m3 > 0.0) || !std::isfinite(volume_m3))
    {
        return result;
    }

    const auto effective_densities =
        estimateCollisionBackgroundDensities(container, volume_m3, runtime_config.floor_densities);
    result.effective_neutral_density_m3 = effective_densities.neutral_density_m3;
    result.effective_ion_density_m3 = effective_densities.ion_density_m3;
    result.effective_electron_density_m3 = effective_densities.electron_density_m3;

    if (runtime_config.reinitialize_configuration)
    {
        initializeConfiguredCollisionHandler(handler, runtime_config.cross_section_set_id,
                                             effective_densities.neutral_density_m3,
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


