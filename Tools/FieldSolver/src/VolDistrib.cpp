#include "VolDistrib.h"

#include "../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SCDAT
{
namespace Field
{

namespace Constants = SCDAT::Basic::Constants;

VolDistrib::VolDistrib(size_t nodeCount, ParticleType particleType)
    : nodeCount_(nodeCount), particleType_(particleType)
{
    if (nodeCount_ == 0)
    {
        throw std::invalid_argument("VolDistrib node count must be positive");
    }
}

VolDistrib::~VolDistrib() = default;

size_t VolDistrib::getNodeCount() const
{
    return nodeCount_;
}

ParticleType VolDistrib::getParticleType() const
{
    return particleType_;
}

ParticleDensityDistrib::ParticleDensityDistrib(size_t nodeCount, ParticleType particleType)
    : VolDistrib(nodeCount, particleType), density_(nodeCount)
{
}

void ParticleDensityDistrib::setUniformDensity(double density)
{
    for (size_t i = 0; i < nodeCount_; ++i)
    {
        density_.setValue(i, density);
    }
}

void ParticleDensityDistrib::setDensityProfile(std::function<double(size_t)> profile)
{
    for (size_t i = 0; i < nodeCount_; ++i)
    {
        density_.setValue(i, profile(i));
    }
}

void ParticleDensityDistrib::addGaussianDistribution(size_t centerIndex, double amplitude,
                                                     double sigma)
{
    if (sigma <= 0.0)
    {
        return;
    }

    for (size_t i = 0; i < nodeCount_; ++i)
    {
        const double d =
            static_cast<double>((i > centerIndex) ? (i - centerIndex) : (centerIndex - i));
        const double g = amplitude * std::exp(-0.5 * d * d / (sigma * sigma));
        density_.setValue(i, density_.getValue(i) + g);
    }
}

double ParticleDensityDistrib::getDensity(size_t nodeIndex) const
{
    return density_.getValue(nodeIndex);
}

const ScalarField& ParticleDensityDistrib::getDensityField() const
{
    return density_;
}

void ParticleDensityDistrib::computeMoment(ScalarField& result, int order) const
{
    if (result.size() != nodeCount_)
    {
        throw std::invalid_argument("moment result size mismatch");
    }

    for (size_t i = 0; i < nodeCount_; ++i)
    {
        result.setValue(i,
                        std::pow(std::max(0.0, density_.getValue(i)), static_cast<double>(order)));
    }
}

void ParticleDensityDistrib::updateFromParticles(const std::vector<ParticlePtr>& particles)
{
    density_.reset();

    if (particles.empty() || nodeCount_ == 0)
    {
        return;
    }

    const double uniform = static_cast<double>(particles.size()) / static_cast<double>(nodeCount_);
    for (size_t i = 0; i < nodeCount_; ++i)
    {
        density_.setValue(i, uniform);
    }
}

ChargeDensityDistrib::ChargeDensityDistrib(size_t nodeCount, ParticleType particleType)
    : VolDistrib(nodeCount, particleType), chargeDensity_(nodeCount)
{
}

void ChargeDensityDistrib::computeFromParticleDensity(const ParticleDensityDistrib& particleDensity)
{
    if (particleDensity.getNodeCount() != nodeCount_)
    {
        throw std::invalid_argument("particle density node count mismatch");
    }

    double q = 0.0;
    switch (particleType_)
    {
    case ParticleType::ELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::BACKSCATTERED_ELECTRON:
    case ParticleType::AUGER_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
        q = -Constants::PhysicsConstants::ElementaryCharge;
        break;
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
    case ParticleType::ALPHA:
        q = Constants::PhysicsConstants::ElementaryCharge;
        break;
    case ParticleType::NEGATIVE_ION:
    case ParticleType::BETA_PARTICLE:
        q = -Constants::PhysicsConstants::ElementaryCharge;
        break;
    default:
        q = 0.0;
        break;
    }

    for (size_t i = 0; i < nodeCount_; ++i)
    {
        chargeDensity_.setValue(i, particleDensity.getDensity(i) * q);
    }
}

void ChargeDensityDistrib::addContribution(const ParticleDensityDistrib& particleDensity,
                                           double chargeMultiplier)
{
    if (particleDensity.getNodeCount() != nodeCount_)
    {
        throw std::invalid_argument("particle density node count mismatch");
    }

    for (size_t i = 0; i < nodeCount_; ++i)
    {
        const double v =
            chargeDensity_.getValue(i) + particleDensity.getDensity(i) * chargeMultiplier;
        chargeDensity_.setValue(i, v);
    }
}

double ChargeDensityDistrib::getChargeDensity(size_t nodeIndex) const
{
    return chargeDensity_.getValue(nodeIndex);
}

const ScalarField& ChargeDensityDistrib::getChargeDensityField() const
{
    return chargeDensity_;
}

EnergyDensityDistrib::EnergyDensityDistrib(size_t nodeCount, ParticleType particleType)
    : VolDistrib(nodeCount, particleType), energyDensity_(nodeCount)
{
}

void EnergyDensityDistrib::computeKineticEnergy(const ParticleDensityDistrib& particleDensity,
                                                const VectorField& velocityField)
{
    if (particleDensity.getNodeCount() != nodeCount_ || velocityField.size() != nodeCount_)
    {
        throw std::invalid_argument("input field size mismatch");
    }

    double mass = Constants::PhysicsConstants::ElectronMass;
    switch (particleType_)
    {
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
    case ParticleType::ALPHA:
        mass = Constants::PhysicsConstants::ProtonMass;
        break;
    default:
        mass = Constants::PhysicsConstants::ElectronMass;
        break;
    }

    for (size_t i = 0; i < nodeCount_; ++i)
    {
        const double n = particleDensity.getDensity(i);
        const Vector3D v = velocityField.getValue(i);
        const double ek = 0.5 * mass * v.magnitudeSquared();
        energyDensity_.setValue(i, n * ek);
    }
}

double EnergyDensityDistrib::getEnergyDensity(size_t nodeIndex) const
{
    return energyDensity_.getValue(nodeIndex);
}

const ScalarField& EnergyDensityDistrib::getEnergyDensityField() const
{
    return energyDensity_;
}

void EnergyDensityDistrib::updateFromParticles(const std::vector<ParticlePtr>& particles)
{
    energyDensity_.reset();

    if (particles.empty() || nodeCount_ == 0)
    {
        return;
    }

    double total = 0.0;
    for (const auto* p : particles)
    {
        if (p && p->isActive())
        {
            total += p->getKineticEnergy();
        }
    }

    const double perNode = total / static_cast<double>(nodeCount_);
    for (size_t i = 0; i < nodeCount_; ++i)
    {
        energyDensity_.setValue(i, perNode);
    }
}

} // namespace Field
} // namespace SCDAT
