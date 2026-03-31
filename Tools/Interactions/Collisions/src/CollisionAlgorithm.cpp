/**
 * @file CollisionAlgorithm.cpp
 * @brief Monte Carlo 与 DSMC 碰撞算法实现
 */

#include "../include/CollisionAlgorithm.h"
#include "../../../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace SCDAT {
namespace Collision {

namespace {
constexpr double kEps = 1e-30;
constexpr double kPi = 3.14159265358979323846;
} // namespace

size_t processEnabledCollisionModels(std::vector<std::shared_ptr<CollisionModel>>& models,
                                     std::vector<ParticleObject>& particles,
                                     std::vector<ParticleObject>& new_particles, double dt,
                                     double volume)
{
    size_t total = 0;
    for (auto& model : models)
    {
        if (model && model->isEnabled())
        {
            total += model->processCollisions(particles, new_particles, dt, volume);
        }
    }
    return total;
}

size_t processAllCollisions(CollisionModelManager& manager,
                            std::vector<ParticleObject>& particles,
                            std::vector<ParticleObject>& new_particles, double dt,
                            double volume)
{
    size_t total = 0;
    for (const auto& model : manager.getModels())
    {
        if (model && model->isEnabled())
        {
            total += model->processCollisions(particles, new_particles, dt, volume);
        }
    }
    return total;
}

CollisionProcess::CollisionProcess(CollisionCrossSectionPtr cross_section,
                                   ParticleSpecies incident_species,
                                   ParticleSpecies target_species, double target_density)
    : cross_section_(std::move(cross_section)), incident_species_(incident_species),
      target_species_(target_species), target_density_(target_density)
{
    if (!cross_section_)
    {
        throw std::invalid_argument("CollisionProcess: cross_section is null");
    }
    if (target_density_ < 0.0)
    {
        throw std::invalid_argument("CollisionProcess: target_density must be non-negative");
    }
}

double CollisionProcess::calculateCollisionFrequency(const ParticleObject& particle) const
{
    const double energy = calculateKineticEnergyEv(particle);
    const double sigma = cross_section_->calculateCrossSection(energy);
    const double speed = particle.getVelocity().magnitude();
    return target_density_ * sigma * speed;
}

bool CollisionProcess::isApplicableTo(const ParticleObject& particle) const
{
    return classifyParticleSpecies(particle) == incident_species_;
}

ParticleSpecies CollisionProcess::getIncidentSpecies() const
{
    return incident_species_;
}

CollisionCrossSectionPtr CollisionProcess::getCrossSection() const
{
    return cross_section_;
}

ParticleSpecies CollisionProcess::getTargetSpecies() const
{
    return target_species_;
}

void CollisionProcess::setTargetDensity(double density)
{
    target_density_ = std::max(0.0, density);
}

double CollisionProcess::getTargetDensity() const
{
    return target_density_;
}

double CollisionProcess::calculateKineticEnergyEv(const ParticleObject& particle) const
{
    const double speed = particle.getVelocity().magnitude();
    const double energy_j = 0.5 * particle.getMass() * speed * speed;
    return energy_j / SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge;
}

double CollisionProcess::generateIsotropicScatteringAngle(std::mt19937& rng) const
{
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    return std::acos(u(rng));
}

double CollisionProcess::generateRandomAzimuthalAngle(std::mt19937& rng) const
{
    std::uniform_real_distribution<double> u(0.0, 2.0 * kPi);
    return u(rng);
}

ParticleSpecies CollisionProcess::classifyParticleSpecies(const ParticleObject& particle)
{
    using SCDAT::Particle::ParticleType;

    switch (particle.getType())
    {
    case ParticleType::ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::BACKSCATTERED_ELECTRON:
    case ParticleType::AUGER_ELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
    case ParticleType::POSITRON:
        return ParticleSpecies::ELECTRON;

    case ParticleType::PROTON:
    case ParticleType::ALPHA:
    case ParticleType::DEUTERON:
    case ParticleType::TRITON:
    case ParticleType::HELION:
    case ParticleType::HEAVY_ION:
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
        return ParticleSpecies::ION;

    case ParticleType::ATOM:
    case ParticleType::MOLECULE:
    case ParticleType::RADICAL:
        return ParticleSpecies::NEUTRAL;

    case ParticleType::PHOTON:
    case ParticleType::GAMMA:
    case ParticleType::XRAY:
    case ParticleType::OPTICAL_PHOTON:
    case ParticleType::VIRTUAL_PHOTON:
        return ParticleSpecies::PHOTON;

    default:
        if (particle.getCharge() < 0.0)
        {
            return ParticleSpecies::ELECTRON;
        }
        if (particle.getCharge() > 0.0)
        {
            return ParticleSpecies::ION;
        }
        return ParticleSpecies::NEUTRAL;
    }
}

ElasticCollisionProcess::ElasticCollisionProcess(CollisionCrossSectionPtr cross_section,
                                                 ParticleSpecies incident_species,
                                                 ParticleSpecies target_species,
                                                 double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> ElasticCollisionProcess::executeCollision(ParticleObject& particle,
                                                                      double dt,
                                                                      std::mt19937& rng)
{
    std::vector<ParticleObject> no_new_particles;

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return no_new_particles;
    }

    const double speed = particle.getVelocity().magnitude();
    const double theta = generateIsotropicScatteringAngle(rng);
    const double phi = generateRandomAzimuthalAngle(rng);

    const double sin_t = std::sin(theta);
    particle.setVelocity(Vector3D(speed * sin_t * std::cos(phi), speed * sin_t * std::sin(phi),
                                  speed * std::cos(theta)));
    return no_new_particles;
}

InelasticCollisionProcess::InelasticCollisionProcess(CollisionCrossSectionPtr cross_section,
                                                     ParticleSpecies incident_species,
                                                     ParticleSpecies target_species,
                                                     double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> InelasticCollisionProcess::executeCollision(ParticleObject& particle,
                                                                        double dt,
                                                                        std::mt19937& rng)
{
    std::vector<ParticleObject> no_new_particles;

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return no_new_particles;
    }

    const double initial_ev = calculateKineticEnergyEv(particle);
    const double loss_ev = std::max(0.0, cross_section_->getParameters().energy_loss);
    const double remaining_ev = std::max(1e-6, initial_ev - loss_ev);

    const double speed = std::sqrt(
        2.0 * remaining_ev * SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge /
        std::max(kEps, particle.getMass()));
    const double theta = generateIsotropicScatteringAngle(rng);
    const double phi = generateRandomAzimuthalAngle(rng);
    const double sin_t = std::sin(theta);
    particle.setVelocity(Vector3D(speed * sin_t * std::cos(phi), speed * sin_t * std::sin(phi),
                                  speed * std::cos(theta)));

    return no_new_particles;
}

IonizationCollisionProcess::IonizationCollisionProcess(CollisionCrossSectionPtr cross_section,
                                                       ParticleSpecies incident_species,
                                                       ParticleSpecies target_species,
                                                       double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> IonizationCollisionProcess::executeCollision(ParticleObject& particle,
                                                                         double dt,
                                                                         std::mt19937& rng)
{
    std::vector<ParticleObject> new_particles;

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return new_particles;
    }

    const double threshold = cross_section_->getParameters().energy_threshold;
    const double initial_ev = calculateKineticEnergyEv(particle);
    if (initial_ev <= threshold)
    {
        return new_particles;
    }

    const double available_ev = initial_ev - threshold;
    const double secondary_ev = u(rng) * available_ev;
    const double primary_ev = available_ev - secondary_ev;

    const double m = std::max(kEps, particle.getMass());
    const double primary_speed = std::sqrt(
        2.0 * primary_ev * SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge / m);
    const double secondary_speed =
        std::sqrt(2.0 * secondary_ev *
                  SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge /
                  SCDAT::Basic::Constants::PhysicsConstants::ElectronMass);

    const double theta1 = generateIsotropicScatteringAngle(rng);
    const double phi1 = generateRandomAzimuthalAngle(rng);
    const double theta2 = generateIsotropicScatteringAngle(rng);
    const double phi2 = generateRandomAzimuthalAngle(rng);

    particle.setVelocity(Vector3D(primary_speed * std::sin(theta1) * std::cos(phi1),
                                  primary_speed * std::sin(theta1) * std::sin(phi1),
                                  primary_speed * std::cos(theta1)));

    static SCDAT::Particle::ParticleId next_id = 2000000;
    new_particles.push_back(SCDAT::Particle::ParticleFactory::createElectron(
        next_id++, particle.getPosition(),
        Vector3D(secondary_speed * std::sin(theta2) * std::cos(phi2),
                 secondary_speed * std::sin(theta2) * std::sin(phi2),
                 secondary_speed * std::cos(theta2)),
        particle.getWeight()));

    const auto& params = cross_section_->getParameters();
    const int ion_mass_number =
        params.additional_params.contains("ion_mass_number")
            ? static_cast<int>(params.additional_params.at("ion_mass_number"))
            : 40;
    const int ion_charge_number =
        params.additional_params.contains("ion_charge_number")
            ? static_cast<int>(params.additional_params.at("ion_charge_number"))
            : 1;

    new_particles.push_back(SCDAT::Particle::ParticleFactory::createIon(
        next_id++, particle.getPosition(), Vector3D(0.0, 0.0, 0.0), ion_mass_number,
        ion_charge_number, particle.getWeight()));
    return new_particles;
}

ExcitationCollisionProcess::ExcitationCollisionProcess(CollisionCrossSectionPtr cross_section,
                                                       ParticleSpecies incident_species,
                                                       ParticleSpecies target_species,
                                                       double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> ExcitationCollisionProcess::executeCollision(ParticleObject& particle,
                                                                         double dt,
                                                                         std::mt19937& rng)
{
    std::vector<ParticleObject> no_new_particles;
    if (particle.getStatus() != SCDAT::Particle::ParticleStatus::ACTIVE)
    {
        return no_new_particles;
    }

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return no_new_particles;
    }

    const double initial_ev = calculateKineticEnergyEv(particle);
    const double threshold_ev = cross_section_->getParameters().energy_threshold;
    if (initial_ev <= threshold_ev)
    {
        return no_new_particles;
    }

    const double configured_loss_ev = cross_section_->getParameters().energy_loss;
    const double loss_ev = configured_loss_ev > 0.0 ? configured_loss_ev : threshold_ev;
    const double remaining_ev = std::max(1e-6, initial_ev - loss_ev);

    const double speed = std::sqrt(
        2.0 * remaining_ev * SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge /
        std::max(kEps, particle.getMass()));
    const double theta = generateIsotropicScatteringAngle(rng);
    const double phi = generateRandomAzimuthalAngle(rng);
    const double sin_t = std::sin(theta);

    particle.setVelocity(Vector3D(speed * sin_t * std::cos(phi), speed * sin_t * std::sin(phi),
                                  speed * std::cos(theta)));
    return no_new_particles;
}

RecombinationCollisionProcess::RecombinationCollisionProcess(
    CollisionCrossSectionPtr cross_section, ParticleSpecies incident_species,
    ParticleSpecies target_species, double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> RecombinationCollisionProcess::executeCollision(
    ParticleObject& particle, double dt, std::mt19937& rng)
{
    std::vector<ParticleObject> no_new_particles;
    if (particle.getStatus() != SCDAT::Particle::ParticleStatus::ACTIVE)
    {
        return no_new_particles;
    }

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return no_new_particles;
    }

    particle.setCharge(0.0);
    particle.setVelocity(Vector3D(0.0, 0.0, 0.0));
    particle.setStatus(SCDAT::Particle::ParticleStatus::RECOMBINED);
    return no_new_particles;
}

ChargeExchangeCollisionProcess::ChargeExchangeCollisionProcess(
    CollisionCrossSectionPtr cross_section, ParticleSpecies incident_species,
    ParticleSpecies target_species, double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> ChargeExchangeCollisionProcess::executeCollision(
    ParticleObject& particle, double dt, std::mt19937& rng)
{
    std::vector<ParticleObject> no_new_particles;
    if (particle.getStatus() != SCDAT::Particle::ParticleStatus::ACTIVE)
    {
        return no_new_particles;
    }

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return no_new_particles;
    }

    if (particle.getCharge() > 0.0)
    {
        particle.setCharge(0.0);
    }
    else if (particle.getCharge() < 0.0)
    {
        particle.setCharge(-SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge);
    }

    return no_new_particles;
}

CoulombCollisionProcess::CoulombCollisionProcess(CollisionCrossSectionPtr cross_section,
                                                 ParticleSpecies incident_species,
                                                 ParticleSpecies target_species,
                                                 double target_density)
    : CollisionProcess(std::move(cross_section), incident_species, target_species, target_density)
{
}

std::vector<ParticleObject> CoulombCollisionProcess::executeCollision(ParticleObject& particle,
                                                                      double dt,
                                                                      std::mt19937& rng)
{
    std::vector<ParticleObject> no_new_particles;
    if (particle.getCharge() == 0.0)
    {
        return no_new_particles;
    }

    const double p = std::clamp(calculateCollisionFrequency(particle) * dt, 0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) >= p)
    {
        return no_new_particles;
    }

    const double theta = generateIsotropicScatteringAngle(rng);
    const double phi = generateRandomAzimuthalAngle(rng);
    const double speed = particle.getVelocity().magnitude();
    const double sin_t = std::sin(theta);

    particle.setVelocity(Vector3D(speed * sin_t * std::cos(phi), speed * sin_t * std::sin(phi),
                                  speed * std::cos(theta)));
    return no_new_particles;
}

MonteCarloCollisionHandler::MonteCarloCollisionHandler(unsigned int random_seed)
    : rng_(random_seed), collision_enabled_(true)
{
    initializeDefaultConfiguration();
}

void MonteCarloCollisionHandler::addCollisionProcess(const std::string& process_name,
                                                     std::shared_ptr<CollisionProcess> process)
{
    if (!process)
    {
        throw std::invalid_argument("MonteCarloCollisionHandler: process must not be null");
    }
    collision_processes_[process_name] = std::move(process);
    statistics_.collision_counts[process_name] = 0;
    statistics_.total_rates[process_name] = 0.0;
}

void MonteCarloCollisionHandler::removeCollisionProcess(const std::string& process_name)
{
    collision_processes_.erase(process_name);
    statistics_.collision_counts.erase(process_name);
    statistics_.total_rates.erase(process_name);
}

void MonteCarloCollisionHandler::clearCollisionProcesses()
{
    collision_processes_.clear();
    resetStatistics();
}

std::vector<ParticleObject> MonteCarloCollisionHandler::processCollisions(
    std::vector<ParticleObject>& particles, double dt)
{
    std::vector<ParticleObject> new_particles;
    if (!collision_enabled_ || collision_processes_.empty())
    {
        return new_particles;
    }

    for (auto& p : particles)
    {
        if (p.getStatus() != SCDAT::Particle::ParticleStatus::ACTIVE)
        {
            continue;
        }

        for (const auto& kv : collision_processes_)
        {
            const std::string& name = kv.first;
            const auto& process = kv.second;
            if (!process || !process->isApplicableTo(p))
            {
                continue;
            }

            const double freq = process->calculateCollisionFrequency(p);
            const double prob = std::clamp(freq * dt, 0.0, 1.0);
            if (shouldCollisionOccur(prob))
            {
                auto generated = process->executeCollision(p, dt, rng_);
                updateStatistics(name, static_cast<int>(generated.size()));
                statistics_.total_rates[name] += freq;
                new_particles.insert(new_particles.end(), generated.begin(), generated.end());
            }
        }
    }

    return new_particles;
}

const MonteCarloCollisionHandler::CollisionStatistics&
MonteCarloCollisionHandler::getStatistics() const
{
    return statistics_;
}

void MonteCarloCollisionHandler::resetStatistics()
{
    statistics_.collision_counts.clear();
    statistics_.total_rates.clear();
    statistics_.total_collisions = 0;
    statistics_.particles_created = 0;

    for (const auto& kv : collision_processes_)
    {
        statistics_.collision_counts[kv.first] = 0;
        statistics_.total_rates[kv.first] = 0.0;
    }
}

void MonteCarloCollisionHandler::setCollisionEnabled(bool enabled)
{
    collision_enabled_ = enabled;
}

bool MonteCarloCollisionHandler::isCollisionEnabled() const
{
    return collision_enabled_;
}

const CollisionParameters&
MonteCarloCollisionHandler::getDefaultCollisionParameters(CollisionType collision_type) const
{
    return default_collision_parameters_.at(collision_type);
}

void MonteCarloCollisionHandler::setDefaultCollisionParameters(
    CollisionType collision_type, const CollisionParameters& params)
{
    default_collision_parameters_[collision_type] = params;
}

void MonteCarloCollisionHandler::initializeDefaultConfiguration(double neutral_density,
                                                               double ion_density,
                                                               double electron_density)
{
    initializeDefaultCollisionParameters();
    registerDefaultCollisionProcesses(neutral_density, ion_density, electron_density);
}

void MonteCarloCollisionHandler::initializeDefaultCollisionParameters()
{
    CollisionParameters elastic;
    elastic.energy_threshold = 0.0;
    elastic.cross_section_max = 2.0e-19;
    elastic.energy_loss = 0.0;
    elastic.reaction_equation = "A + B -> A + B";
    elastic.additional_params["reference_energy"] = 10.0;
    elastic.additional_params["power_index"] = 0.3;
    default_collision_parameters_[CollisionType::ELASTIC] = elastic;

    CollisionParameters inelastic;
    inelastic.energy_threshold = 1.0;
    inelastic.cross_section_max = 8.0e-20;
    inelastic.energy_loss = 1.0;
    inelastic.reaction_equation = "A + B -> A* + B";
    default_collision_parameters_[CollisionType::INELASTIC] = inelastic;

    CollisionParameters ionization;
    ionization.energy_threshold = 15.76;
    ionization.cross_section_max = 3.0e-20;
    ionization.energy_loss = 15.76;
    ionization.reaction_equation = "e + A -> 2e + A+";
    default_collision_parameters_[CollisionType::IONIZATION] = ionization;

    CollisionParameters excitation;
    excitation.energy_threshold = 11.5;
    excitation.cross_section_max = 2.5e-20;
    excitation.energy_loss = 11.5;
    excitation.reaction_equation = "e + A -> e + A*";
    excitation.additional_params["alpha"] = 1.2;
    excitation.additional_params["beta"] = 0.35;
    default_collision_parameters_[CollisionType::EXCITATION] = excitation;

    CollisionParameters charge_exchange;
    charge_exchange.energy_threshold = 0.0;
    charge_exchange.cross_section_max = 1.2e-19;
    charge_exchange.energy_loss = 0.2;
    charge_exchange.reaction_equation = "A+ + B -> A + B+";
    default_collision_parameters_[CollisionType::CHARGE_EXCHANGE] = charge_exchange;

    CollisionParameters recombination;
    recombination.energy_threshold = 0.0;
    recombination.cross_section_max = 5.0e-19;
    recombination.energy_loss = 0.0;
    recombination.reaction_equation = "e + A+ -> A";
    recombination.additional_params["reference_energy"] = 1.0;
    recombination.additional_params["power_index"] = 0.7;
    default_collision_parameters_[CollisionType::RECOMBINATION] = recombination;

    CollisionParameters coulomb;
    coulomb.energy_threshold = 0.0;
    coulomb.cross_section_max = 1.0e-18;
    coulomb.reaction_equation = "q1 + q2 -> q1 + q2";
    coulomb.additional_params["electron_density"] = 1.0e16;
    coulomb.additional_params["electron_temperature"] = 1000.0;
    default_collision_parameters_[CollisionType::COULOMB] = coulomb;
}

void MonteCarloCollisionHandler::registerDefaultCollisionProcesses(double neutral_density,
                                                                   double ion_density,
                                                                   double electron_density)
{
    collision_processes_.clear();

    addCollisionProcess(
        "elastic_e_n",
        CollisionProcessFactory::createElasticProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density,
            default_collision_parameters_.at(CollisionType::ELASTIC)));

    addCollisionProcess(
        "inelastic_e_n",
        CollisionProcessFactory::createInelasticProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density,
            default_collision_parameters_.at(CollisionType::INELASTIC)));

    addCollisionProcess(
        "ionization_e_n",
        CollisionProcessFactory::createIonizationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density,
            default_collision_parameters_.at(CollisionType::IONIZATION)));

    addCollisionProcess(
        "excitation_e_n",
        CollisionProcessFactory::createExcitationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::NEUTRAL, neutral_density,
            default_collision_parameters_.at(CollisionType::EXCITATION)));

    addCollisionProcess(
        "charge_exchange_i_n",
        CollisionProcessFactory::createChargeExchangeProcess(
            ParticleSpecies::ION, ParticleSpecies::NEUTRAL, neutral_density,
            default_collision_parameters_.at(CollisionType::CHARGE_EXCHANGE)));

    addCollisionProcess(
        "recombination_e_i",
        CollisionProcessFactory::createRecombinationProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::ION, ion_density,
            default_collision_parameters_.at(CollisionType::RECOMBINATION)));

    addCollisionProcess(
        "coulomb_e_i",
        CollisionProcessFactory::createCoulombProcess(
            ParticleSpecies::ELECTRON, ParticleSpecies::ION, electron_density,
            default_collision_parameters_.at(CollisionType::COULOMB)));

    resetStatistics();
}

bool MonteCarloCollisionHandler::shouldCollisionOccur(double collision_probability)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    return u(rng_) < collision_probability;
}

void MonteCarloCollisionHandler::updateStatistics(const std::string& process_name,
                                                  int new_particles_count)
{
    ++statistics_.collision_counts[process_name];
    ++statistics_.total_collisions;
    statistics_.particles_created += new_particles_count;
}

DSMCCollisionHandler::DSMCCollisionHandler(unsigned int random_seed) : rng_(random_seed)
{
}

void DSMCCollisionHandler::setModel(std::shared_ptr<CollisionModel> model)
{
    model_ = std::move(model);
}

size_t DSMCCollisionHandler::processCollisions(std::vector<ParticleObject>& particles,
                                               std::vector<ParticleObject>& new_particles,
                                               double dt, double cell_volume)
{
    if (!model_ || particles.size() < 2 || cell_volume <= 0.0)
    {
        return 0;
    }

    // DSMC 核心思想：每个时间步在单元内随机配对代表粒子，再按碰撞概率处理。
    std::shuffle(particles.begin(), particles.end(), rng_);
    return model_->processCollisions(particles, new_particles, dt, cell_volume);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createElasticProcess(ParticleSpecies incident_species,
                                              ParticleSpecies target_species,
                                              double target_density,
                                              const CollisionParameters& params)
{
    auto cs = std::make_shared<ElasticCollisionCrossSection>(params);
    return std::make_shared<ElasticCollisionProcess>(cs, incident_species, target_species,
                                                     target_density);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createIonizationProcess(ParticleSpecies incident_species,
                                                 ParticleSpecies target_species,
                                                 double target_density,
                                                 const CollisionParameters& params)
{
    auto cs = std::make_shared<IonizationCollisionCrossSection>(params);
    return std::make_shared<IonizationCollisionProcess>(cs, incident_species, target_species,
                                                        target_density);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createInelasticProcess(ParticleSpecies incident_species,
                                                ParticleSpecies target_species,
                                                double target_density,
                                                const CollisionParameters& params)
{
    auto cs = std::make_shared<ElasticCollisionCrossSection>(params);
    return std::make_shared<InelasticCollisionProcess>(cs, incident_species, target_species,
                                                       target_density);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createCoulombProcess(ParticleSpecies incident_species,
                                              ParticleSpecies target_species,
                                              double target_density,
                                              const CollisionParameters& params)
{
    auto cs = std::make_shared<CoulombCollisionCrossSection>(params);
    return std::make_shared<CoulombCollisionProcess>(cs, incident_species, target_species,
                                                     target_density);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createChargeExchangeProcess(ParticleSpecies incident_species,
                                                     ParticleSpecies target_species,
                                                     double target_density,
                                                     const CollisionParameters& params)
{
    auto cs = std::make_shared<ElasticCollisionCrossSection>(params);
    return std::make_shared<ChargeExchangeCollisionProcess>(cs, incident_species, target_species,
                                                            target_density);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createExcitationProcess(ParticleSpecies incident_species,
                                                 ParticleSpecies target_species,
                                                 double target_density,
                                                 const CollisionParameters& params)
{
    auto cs = std::make_shared<ExcitationCollisionCrossSection>(params);
    return std::make_shared<ExcitationCollisionProcess>(cs, incident_species, target_species,
                                                        target_density);
}

std::shared_ptr<CollisionProcess>
CollisionProcessFactory::createRecombinationProcess(ParticleSpecies incident_species,
                                                    ParticleSpecies target_species,
                                                    double target_density,
                                                    const CollisionParameters& params)
{
    auto cs = std::make_shared<RecombinationCollisionCrossSection>(params);
    return std::make_shared<RecombinationCollisionProcess>(cs, incident_species, target_species,
                                                           target_density);
}

} // namespace Collision
} // namespace SCDAT
