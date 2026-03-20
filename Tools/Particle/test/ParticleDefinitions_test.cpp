#include "../include/ParticleDefinitions.h"
#include "../include/ParticlePusher.h"
#include <gtest/gtest.h>

#include <cmath>
#include <string>

namespace SCDAT
{
namespace Particle
{

using SCDAT::Basic::Constants::PhysicsConstants;

TEST(ParticleDefinitionsTest, FactoryCreateElectronUsesGlobalConstants)
{
    Particle p = ParticleFactory::createElectron(1, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                                                 SCDAT::Geometry::Vector3D(1.0, 0.0, 0.0), 2.0);

    EXPECT_EQ(p.getType(), ParticleType::ELECTRON);
    EXPECT_DOUBLE_EQ(p.getMass(), PhysicsConstants::ElectronMass);
    EXPECT_DOUBLE_EQ(p.getCharge(), PhysicsConstants::ElectronCharge);
    EXPECT_DOUBLE_EQ(p.getWeight(), 2.0);
}

TEST(ParticleDefinitionsTest, FactoryCreateIonAssignsSpecies)
{
    Particle p = ParticleFactory::createIon(2, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                                            SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0), 4, 2, 1.5);

    EXPECT_EQ(p.getType(), ParticleType::ION);
    EXPECT_EQ(p.getMassNumber(), 4);
    EXPECT_EQ(p.getChargeNumber(), 2);
    EXPECT_EQ(p.getSpeciesName(), "He++");
    EXPECT_DOUBLE_EQ(p.getMass(), 4.0 * PhysicsConstants::ProtonMass);
    EXPECT_DOUBLE_EQ(p.getCharge(), 2.0 * PhysicsConstants::ElementaryCharge);
}

TEST(ParticleDefinitionsTest, LifecycleAndAgeManagement)
{
    Particle p(3, ParticleType::ELECTRON, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
               SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0), 1.0, -1.0, 1.0);

    p.updateAge(0.25);
    EXPECT_DOUBLE_EQ(p.getAge(), 0.25);
    EXPECT_TRUE(p.isActive());

    p.markAsRecombined();
    EXPECT_EQ(p.getStatus(), ParticleStatus::RECOMBINED);

    p.markAsAbsorbed();
    EXPECT_EQ(p.getStatus(), ParticleStatus::ABSORBED);

    p.markAsEscaped();
    EXPECT_EQ(p.getStatus(), ParticleStatus::ESCAPED);
}

TEST(ParticleDefinitionsTest, BorisPushPureElectricField)
{
    Particle p(4, ParticleType::ELECTRON, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
               SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);

    const SCDAT::Geometry::Vector3D e(2.0, 0.0, 0.0);
    const SCDAT::Geometry::Vector3D b(0.0, 0.0, 0.0);
    const double dt = 0.5;

    BorisAlgorithm pusher(dt);
    const FieldFunction ef = [&e](const SCDAT::Geometry::Point3D&) { return e; };
    const FieldFunction bf = [&b](const SCDAT::Geometry::Point3D&) { return b; };
    pusher.pushParticle(&p, ef, bf);

    EXPECT_NEAR(p.getVelocity().x(), 1.0, 1e-12);
    EXPECT_NEAR(p.getVelocity().y(), 0.0, 1e-12);
    EXPECT_NEAR(p.getPosition().x(), 0.5, 1e-12);
    EXPECT_NEAR(p.getPosition().y(), 0.0, 1e-12);
}

TEST(ParticleDefinitionsTest, RelativisticPathWithZeroFieldKeepsVelocity)
{
    const double c = PhysicsConstants::SpeedOfLight;
    const double vx = 0.2 * c;
    const double dt = 1e-9;

    Particle p(5, ParticleType::ELECTRON, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
               SCDAT::Geometry::Vector3D(vx, 0.0, 0.0), 1.0, 1.0, 1.0);

    const SCDAT::Geometry::Vector3D zero_field(0.0, 0.0, 0.0);
    BorisAlgorithm pusher(dt);
    pusher.setRelativistic(true);
    const FieldFunction ef = [&zero_field](const SCDAT::Geometry::Point3D&) { return zero_field; };
    const FieldFunction bf = [&zero_field](const SCDAT::Geometry::Point3D&) { return zero_field; };
    pusher.pushParticle(&p, ef, bf);

    EXPECT_GT(p.getVelocity().x(), 0.0);
    EXPECT_LT(std::abs(p.getVelocity().x()), c);
    EXPECT_NEAR(p.getVelocity().y(), 0.0, 1e-12);
    EXPECT_NEAR(p.getVelocity().z(), 0.0, 1e-12);
    EXPECT_NEAR(p.getPosition().x(), p.getVelocity().x() * dt, 1e-12);
    EXPECT_NEAR(p.getAge(), dt, 1e-15);
}

TEST(ParticleDefinitionsTest, ToStringContainsExtendedPhotoelectronFields)
{
    Particle p = ParticleFactory::createPhotoelectron(6, SCDAT::Geometry::Point3D(1.0, 2.0, 3.0),
                                                      SCDAT::Geometry::Vector3D(4.0, 5.0, 6.0), 7.5,
                                                      4.5, 1.0);

    const std::string text = p.toString();

    EXPECT_NE(text.find("Photoelectron"), std::string::npos);
    EXPECT_NE(text.find("photon_energy"), std::string::npos);
    EXPECT_NE(text.find("max_KE"), std::string::npos);
}

} // namespace Particle
} // namespace SCDAT
