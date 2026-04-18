#include "VolInteract.h"

#include "../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <utility>

namespace SCDAT
{
namespace Field
{

namespace Constants = SCDAT::Basic::Constants;

namespace
{

bool isElectronLike(ParticleType type)
{
    switch (type)
    {
    case ParticleType::ELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::BACKSCATTERED_ELECTRON:
    case ParticleType::AUGER_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
        return true;
    default:
        return false;
    }
}

bool isIonLike(ParticleType type)
{
    switch (type)
    {
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
    case ParticleType::ALPHA:
    case ParticleType::BETA_PARTICLE:
        return true;
    default:
        return false;
    }
}

bool isNeutralLike(ParticleType type)
{
    switch (type)
    {
    case ParticleType::ATOM:
    case ParticleType::MOLECULE:
    case ParticleType::RADICAL:
        return true;
    default:
        return false;
    }
}

std::vector<ParticlePtr> collectActiveParticles(const std::vector<ParticlePtr>& particles)
{
    std::vector<ParticlePtr> active;
    active.reserve(particles.size());
    for (auto* particle : particles)
    {
        if (particle != nullptr && particle->isActive())
        {
            active.push_back(particle);
        }
    }
    return active;
}

std::vector<ParticlePtr> collectActiveParticlesByType(const std::vector<ParticlePtr>& particles,
                                                      ParticleType requested_type)
{
    std::vector<ParticlePtr> filtered;
    filtered.reserve(particles.size());
    for (auto* particle : particles)
    {
        if (particle == nullptr || !particle->isActive())
        {
            continue;
        }

        const auto type = particle->getType();
        const bool matched = (requested_type == ParticleType::ELECTRON && isElectronLike(type)) ||
                             (requested_type == ParticleType::ION && isIonLike(type)) ||
                             (requested_type == ParticleType::ATOM && isNeutralLike(type)) ||
                             type == requested_type;
        if (matched)
        {
            filtered.push_back(particle);
        }
    }
    return filtered;
}

} // namespace

VolInteractModel::VolInteractModel(std::string name) : name_(std::move(name))
{
}

VolInteractModel::~VolInteractModel() = default;

const std::string& VolInteractModel::getName() const
{
    return name_;
}

VolInteractor::VolInteractor(const std::string& name,
                             FieldSolver::NativeVolumeInteractionFamily family)
    : name_(name), enabled_(true), family_(family)
{
}

VolInteractor::~VolInteractor() = default;

void VolInteractor::setEnabled(bool enabled)
{
    enabled_ = enabled;
}

bool VolInteractor::isEnabled() const
{
    return enabled_;
}

const std::string& VolInteractor::getName() const
{
    return name_;
}

FieldSolver::NativeVolumeInteractionFamily VolInteractor::getFamily() const
{
    return family_;
}

const VolInteractModel* VolInteractor::getModel() const
{
    return model_.get();
}

void VolInteractor::setModel(std::shared_ptr<VolInteractModel> model)
{
    model_ = std::move(model);
}

ElasticCollisionInteractor::ElasticCollisionInteractor(
    const std::string& name, FieldSolver::NativeVolumeInteractionFamily family,
    const ParticleType& type1,
    const ParticleType& type2, double crossSection)
    : VolInteractor(name, family), type1_(type1), type2_(type2), crossSection_(crossSection),
      reducedMass_(0.0)
{
    const double m1 =
        (type1_ == ParticleType::ELECTRON) ? Constants::PhysicsConstants::ElectronMass
                                            : Constants::PhysicsConstants::ProtonMass;
    const double m2 =
        (type2_ == ParticleType::ELECTRON) ? Constants::PhysicsConstants::ElectronMass
                                            : Constants::PhysicsConstants::ProtonMass;
    reducedMass_ = (m1 + m2 > 0.0) ? (m1 * m2 / (m1 + m2)) : 0.0;
}

std::size_t ElasticCollisionInteractor::execute(std::vector<ParticlePtr>& particles, double dt)
{
    auto particles1 = collectActiveParticlesByType(particles, type1_);
    auto particles2 = collectActiveParticlesByType(particles, type2_);
    return computeInteraction(particles1, particles2, dt);
}

std::size_t ElasticCollisionInteractor::computeInteraction(
    std::vector<ParticlePtr>& particles1, std::vector<ParticlePtr>& particles2, double dt)
{
    if (!enabled_ || particles1.empty() || particles2.empty() || dt <= 0.0)
    {
        return 0;
    }

    const std::size_t interaction_budget =
        std::max<std::size_t>(1U, std::min(particles1.size(), particles2.size()));
    const double coupling = std::clamp(crossSection_ * std::max(0.0, dt) * 1.0e15, 0.0, 1.0);
    if (coupling <= 0.0)
    {
        return 0;
    }
    (void)interaction_budget;
    return 1;
}

void ElasticCollisionInteractor::performElasticCollision(ParticlePtr& p1, ParticlePtr& p2)
{
    if (!p1 || !p2)
    {
        return;
    }

    const auto v1 = p1->getVelocity();
    const auto v2 = p2->getVelocity();
    const double m1 = std::max(1e-30, p1->getMass());
    const double m2 = std::max(1e-30, p2->getMass());

    const auto v1n = ((m1 - m2) * v1 + (2.0 * m2) * v2) / (m1 + m2);
    const auto v2n = ((m2 - m1) * v2 + (2.0 * m1) * v1) / (m1 + m2);

    p1->setVelocity(v1n);
    p2->setVelocity(v2n);
}

ElasticCollisionsInteractor::ElasticCollisionsInteractor(
    const std::string& name, FieldSolver::NativeVolumeInteractionFamily family,
    const ParticleType& type1, const ParticleType& type2, double crossSection)
    : ElasticCollisionInteractor(name, family, type1, type2, crossSection)
{
    setModel(std::make_shared<VolInteractModel>("ElasticCollisions"));
}

MCCInteractor::MCCInteractor()
    : ElasticCollisionsInteractor("MCCInteractor",
                                  FieldSolver::NativeVolumeInteractionFamily::MCCInteractor,
                                  ParticleType::ELECTRON, ParticleType::ION, 1.0e-16)
{
}

CEXInteractor::CEXInteractor()
    : ElasticCollisionsInteractor("CEXInteractor",
                                  FieldSolver::NativeVolumeInteractionFamily::CEXInteractor,
                                  ParticleType::ION, ParticleType::ION, 5.0e-17)
{
}

ConstantIonizationInteractor::ConstantIonizationInteractor(
    const std::string& name, FieldSolver::NativeVolumeInteractionFamily family,
    double ionization_yield_scale)
    : VolInteractor(name, family),
      ionization_yield_scale_(ionization_yield_scale)
{
    setModel(std::make_shared<VolInteractModel>("ConstantIonizationInteractor"));
}

std::size_t ConstantIonizationInteractor::execute(std::vector<ParticlePtr>& particles, double dt)
{
    const auto active_particles = collectActiveParticles(particles);
    if (!enabled_ || active_particles.empty() || dt <= 0.0 || ionization_yield_scale_ <= 0.0)
    {
        return 0;
    }
    return 1;
}

double ConstantIonizationInteractor::ionizationYieldScale() const
{
    return ionization_yield_scale_;
}

PhotoIonizationInteractor::PhotoIonizationInteractor(double ionization_yield_scale)
    : ConstantIonizationInteractor("PhotoIonization",
                                   FieldSolver::NativeVolumeInteractionFamily::PhotoIonization,
                                   ionization_yield_scale)
{
}

VolInteractionAlongTrajectoryInteractor::VolInteractionAlongTrajectoryInteractor(
    const std::string& name, FieldSolver::NativeVolumeInteractionFamily family)
    : VolInteractor(name, family)
{
}

double VolInteractionAlongTrajectoryInteractor::integrationTimeSeconds() const
{
    return integration_time_s_;
}

void VolInteractionAlongTrajectoryInteractor::resetIntegrationTime()
{
    integration_time_s_ = 0.0;
}

void VolInteractionAlongTrajectoryInteractor::recordIntegrationTime(double dt)
{
    if (std::isfinite(dt) && dt > 0.0)
    {
        integration_time_s_ += dt;
    }
}

TrajectoryInteractionFromFieldInteractor::TrajectoryInteractionFromFieldInteractor(
    double field_coupling_scale)
    : VolInteractionAlongTrajectoryInteractor(
          "TrajectoryInteractionFromField",
          FieldSolver::NativeVolumeInteractionFamily::TrajectoryInteractionFromField),
      field_coupling_scale_(field_coupling_scale)
{
}

std::size_t TrajectoryInteractionFromFieldInteractor::execute(
    std::vector<ParticlePtr>& particles, double dt)
{
    const auto active_particles = collectActiveParticles(particles);
    if (!enabled_ || active_particles.empty() || dt <= 0.0 || field_coupling_scale_ <= 0.0)
    {
        return 0;
    }

    recordIntegrationTime(dt);
    for (auto* particle : active_particles)
    {
        auto velocity = particle->getVelocity();
        velocity.setZ(velocity.z() + field_coupling_scale_ * dt);
        particle->setVelocity(velocity);
    }
    return 1;
}

SpinningSpacecraftTrajectoryInteractor::SpinningSpacecraftTrajectoryInteractor(
    double angular_velocity_rad_per_s)
    : VolInteractionAlongTrajectoryInteractor(
          "SpinningSpacecraftTrajectory",
          FieldSolver::NativeVolumeInteractionFamily::SpinningSpacecraftTrajectory),
      angular_velocity_rad_per_s_(angular_velocity_rad_per_s)
{
}

std::size_t SpinningSpacecraftTrajectoryInteractor::execute(
    std::vector<ParticlePtr>& particles, double dt)
{
    const auto active_particles = collectActiveParticles(particles);
    if (!enabled_ || active_particles.empty() || dt <= 0.0 ||
        angular_velocity_rad_per_s_ <= 0.0)
    {
        return 0;
    }

    recordIntegrationTime(dt);
    const double angle = angular_velocity_rad_per_s_ * dt;
    const double cos_angle = std::cos(angle);
    const double sin_angle = std::sin(angle);
    for (auto* particle : active_particles)
    {
        const auto velocity = particle->getVelocity();
        const double rotated_x = cos_angle * velocity.x() - sin_angle * velocity.y();
        const double rotated_y = sin_angle * velocity.x() + cos_angle * velocity.y();
        particle->setVelocity(
            SCDAT::Geometry::Vector3D(rotated_x, rotated_y, velocity.z()));
    }
    return 1;
}

VolInteractionManager::VolInteractionManager() = default;

VolInteractionManager::~VolInteractionManager() = default;

void VolInteractionManager::addInteractor(std::shared_ptr<VolInteractor> interactor)
{
    if (interactor)
    {
        interactors_.push_back(std::move(interactor));
    }
}

void VolInteractionManager::removeInteractor(const std::string& name)
{
    interactors_.erase(std::remove_if(interactors_.begin(), interactors_.end(),
                                      [&name](const std::shared_ptr<VolInteractor>& p)
                                      { return p && p->getName() == name; }),
                       interactors_.end());
}

std::shared_ptr<VolInteractor> VolInteractionManager::getInteractor(const std::string& name)
{
    auto it = std::find_if(interactors_.begin(), interactors_.end(),
                           [&name](const std::shared_ptr<VolInteractor>& p)
                           { return p && p->getName() == name; });
    return (it == interactors_.end()) ? nullptr : *it;
}

std::size_t VolInteractionManager::executeAllInteractions(std::vector<ParticlePtr>& particles,
                                                          double dt)
{
    if (particles.empty() || dt <= 0.0)
    {
        return 0;
    }

    std::size_t executed = 0;
    for (auto& interactor : interactors_)
    {
        if (!interactor || !interactor->isEnabled())
        {
            continue;
        }
        executed += interactor->execute(particles, dt);
    }
    return executed;
}

void VolInteractionManager::setInteractorEnabled(const std::string& name, bool enabled)
{
    auto p = getInteractor(name);
    if (p)
    {
        p->setEnabled(enabled);
    }
}

size_t VolInteractionManager::getInteractorCount() const
{
    return interactors_.size();
}

std::vector<FieldSolver::NativeVolumeInteractionFamily>
VolInteractionManager::activeFamilies() const
{
    std::vector<FieldSolver::NativeVolumeInteractionFamily> families;
    std::unordered_set<int> seen;
    families.reserve(interactors_.size());
    for (const auto& interactor : interactors_)
    {
        if (!interactor || !interactor->isEnabled())
        {
            continue;
        }

        const auto family = interactor->getFamily();
        const auto key = static_cast<int>(family);
        if (seen.insert(key).second)
        {
            families.push_back(family);
        }
    }
    return families;
}

void VolInteractionManager::clear()
{
    interactors_.clear();
}

std::shared_ptr<VolInteractor> makeNativeVolumeInteractor(
    FieldSolver::NativeVolumeInteractionFamily family)
{
    switch (family)
    {
    case FieldSolver::NativeVolumeInteractionFamily::MCCInteractor:
        return std::make_shared<MCCInteractor>();
    case FieldSolver::NativeVolumeInteractionFamily::CEXInteractor:
        return std::make_shared<CEXInteractor>();
    case FieldSolver::NativeVolumeInteractionFamily::PhotoIonization:
        return std::make_shared<PhotoIonizationInteractor>();
    case FieldSolver::NativeVolumeInteractionFamily::TrajectoryInteractionFromField:
        return std::make_shared<TrajectoryInteractionFromFieldInteractor>();
    case FieldSolver::NativeVolumeInteractionFamily::SpinningSpacecraftTrajectory:
        return std::make_shared<SpinningSpacecraftTrajectoryInteractor>();
    }
    return nullptr;
}

VolInteractionManager buildNativeVolumeInteractionManager(
    const std::vector<FieldSolver::NativeVolumeInteractionFamily>& families)
{
    VolInteractionManager manager;
    for (const auto family : families)
    {
        manager.addInteractor(makeNativeVolumeInteractor(family));
    }
    return manager;
}

} // namespace Field
} // namespace SCDAT
