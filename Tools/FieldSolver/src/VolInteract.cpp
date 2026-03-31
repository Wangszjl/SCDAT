#include "VolInteract.h"

#include "../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace SCDAT
{
namespace Field
{

namespace Constants = SCDAT::Basic::Constants;

VolInteractor::VolInteractor(const std::string& name) : name_(name), enabled_(true) {}

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

ElasticCollisionInteractor::ElasticCollisionInteractor(const ParticleType& type1,
                                                       const ParticleType& type2,
                                                       double crossSection)
    : VolInteractor("ElasticCollision"), type1_(type1), type2_(type2), crossSection_(crossSection),
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

void ElasticCollisionInteractor::computeInteraction(std::vector<ParticlePtr>& particles1,
                                                    std::vector<ParticlePtr>& particles2, double dt)
{
    if (!enabled_ || particles1.empty() || particles2.empty() || dt <= 0.0)
    {
        return;
    }

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis01(0.0, 1.0);

    const size_t maxCollisions = std::min(particles1.size(), particles2.size()) / 10 + 1;
    const double pCollide = std::clamp(crossSection_ * dt * 1e15, 0.0, 1.0);

    for (size_t c = 0; c < maxCollisions; ++c)
    {
        const size_t i1 = static_cast<size_t>(dis01(gen) * static_cast<double>(particles1.size())) %
                          particles1.size();
        const size_t i2 = static_cast<size_t>(dis01(gen) * static_cast<double>(particles2.size())) %
                          particles2.size();

        auto* p1 = particles1[i1];
        auto* p2 = particles2[i2];
        if (!p1 || !p2 || !p1->isActive() || !p2->isActive())
        {
            continue;
        }

        if (dis01(gen) < pCollide)
        {
            ParticlePtr pp1 = p1;
            ParticlePtr pp2 = p2;
            performElasticCollision(pp1, pp2);
        }
    }
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

void VolInteractionManager::executeAllInteractions(std::vector<ParticlePtr>& particles, double dt)
{
    if (particles.empty() || dt <= 0.0)
    {
        return;
    }

    std::vector<ParticlePtr> electrons;
    std::vector<ParticlePtr> ions;
    std::vector<ParticlePtr> neutrals;

    for (auto* p : particles)
    {
        if (!p || !p->isActive())
        {
            continue;
        }

        switch (p->getType())
        {
        case ParticleType::ELECTRON:
        case ParticleType::SECONDARY_ELECTRON:
        case ParticleType::PHOTOELECTRON:
        case ParticleType::THERMAL_ELECTRON:
        case ParticleType::BACKSCATTERED_ELECTRON:
        case ParticleType::AUGER_ELECTRON:
        case ParticleType::FIELD_EMISSION_ELECTRON:
            electrons.push_back(p);
            break;
        case ParticleType::ION:
        case ParticleType::POSITIVE_ION:
        case ParticleType::NEGATIVE_ION:
        case ParticleType::MOLECULAR_ION:
        case ParticleType::CLUSTER_ION:
        case ParticleType::ALPHA:
        case ParticleType::BETA_PARTICLE:
            ions.push_back(p);
            break;
        case ParticleType::ATOM:
        case ParticleType::MOLECULE:
        case ParticleType::RADICAL:
            neutrals.push_back(p);
            break;
        default:
            break;
        }
    }

    for (auto& interactor : interactors_)
    {
        if (!interactor || !interactor->isEnabled())
        {
            continue;
        }

        auto elastic = std::dynamic_pointer_cast<ElasticCollisionInteractor>(interactor);
        if (!elastic)
        {
            continue;
        }

        if (!electrons.empty() && !ions.empty())
        {
            elastic->computeInteraction(electrons, ions, dt);
        }
        if (!electrons.empty() && !neutrals.empty())
        {
            elastic->computeInteraction(electrons, neutrals, dt);
        }
    }
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

void VolInteractionManager::clear()
{
    interactors_.clear();
}

} // namespace Field
} // namespace SCDAT
