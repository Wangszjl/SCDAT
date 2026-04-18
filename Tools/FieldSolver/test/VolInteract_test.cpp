#include "../include/SurfaceFieldVolumeBridge.h"
#include "../include/VolInteract.h"

#include <gtest/gtest.h>

#include <type_traits>
#include <vector>

namespace
{

std::vector<SCDAT::Particle::Particle> makeInteractionTestParticles()
{
    std::vector<SCDAT::Particle::Particle> particles;
    particles.emplace_back(1, SCDAT::Particle::ParticleType::ELECTRON,
                           SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                           SCDAT::Geometry::Vector3D(1.0e5, 0.0, 0.0), 9.1093837015e-31,
                           -1.602176634e-19, 1.0);
    particles.emplace_back(2, SCDAT::Particle::ParticleType::ION,
                           SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                           SCDAT::Geometry::Vector3D(0.0, 2.0e3, 0.0), 1.67262192369e-27,
                           1.602176634e-19, 1.0);
    particles.emplace_back(3, SCDAT::Particle::ParticleType::ELECTRON,
                           SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                           SCDAT::Geometry::Vector3D(0.0, 0.0, 3.0e4), 9.1093837015e-31,
                           -1.602176634e-19, 1.0);
    return particles;
}

std::vector<SCDAT::Particle::ParticlePtr> makeParticlePointers(
    std::vector<SCDAT::Particle::Particle>& particles)
{
    std::vector<SCDAT::Particle::ParticlePtr> pointers;
    pointers.reserve(particles.size());
    for (auto& particle : particles)
    {
        pointers.push_back(&particle);
    }
    return pointers;
}

} // namespace

TEST(VolInteractTest, CompatibilityInteractorsExposeExpectedInheritanceLayers)
{
    EXPECT_TRUE((std::is_base_of_v<SCDAT::Field::ElasticCollisionsInteractor,
                                   SCDAT::Field::MCCInteractor>));
    EXPECT_TRUE((std::is_base_of_v<SCDAT::Field::ElasticCollisionsInteractor,
                                   SCDAT::Field::CEXInteractor>));
    EXPECT_TRUE((std::is_base_of_v<SCDAT::Field::ConstantIonizationInteractor,
                                   SCDAT::Field::PhotoIonizationInteractor>));
    EXPECT_TRUE((std::is_base_of_v<SCDAT::Field::VolInteractionAlongTrajectoryInteractor,
                                   SCDAT::Field::TrajectoryInteractionFromFieldInteractor>));
    EXPECT_TRUE((std::is_base_of_v<SCDAT::Field::VolInteractionAlongTrajectoryInteractor,
                                   SCDAT::Field::SpinningSpacecraftTrajectoryInteractor>));

    const SCDAT::Field::MCCInteractor mcc;
    ASSERT_NE(mcc.getModel(), nullptr);
    EXPECT_EQ(mcc.getModel()->getName(), "ElasticCollisions");

    const SCDAT::Field::PhotoIonizationInteractor ionization;
    ASSERT_NE(ionization.getModel(), nullptr);
    EXPECT_EQ(ionization.getModel()->getName(), "ConstantIonizationInteractor");
}

TEST(VolInteractTest, AlongTrajectoryInteractorsAccumulateIntegrationTime)
{
    auto particles = makeInteractionTestParticles();
    auto pointers = makeParticlePointers(particles);

    SCDAT::Field::TrajectoryInteractionFromFieldInteractor field_interactor(2.5);
    ASSERT_EQ(field_interactor.execute(pointers, 2.0e-6), 1U);
    EXPECT_NEAR(field_interactor.integrationTimeSeconds(), 2.0e-6, 1.0e-15);

    SCDAT::Field::SpinningSpacecraftTrajectoryInteractor spinning_interactor(3.0);
    ASSERT_EQ(spinning_interactor.execute(pointers, 4.0e-6), 1U);
    EXPECT_NEAR(spinning_interactor.integrationTimeSeconds(), 4.0e-6, 1.0e-15);
    spinning_interactor.resetIntegrationTime();
    EXPECT_NEAR(spinning_interactor.integrationTimeSeconds(), 0.0, 1.0e-15);
}

TEST(VolInteractTest, ExecuteNativeVolumeInteractionsUsesConfiguredFamilySet)
{
    auto particles = makeInteractionTestParticles();
    auto pointers = makeParticlePointers(particles);

    const auto summary = SCDAT::FieldSolver::executeNativeVolumeInteractions(
        pointers, 1.0e-6,
        {SCDAT::FieldSolver::NativeVolumeInteractionFamily::MCCInteractor,
         SCDAT::FieldSolver::NativeVolumeInteractionFamily::PhotoIonization,
         SCDAT::FieldSolver::NativeVolumeInteractionFamily::SpinningSpacecraftTrajectory});

    EXPECT_EQ(summary.family_count, 3U);
    EXPECT_EQ(summary.interactor_count, 3U);
    EXPECT_EQ(summary.executed_interaction_count, 3U);
    EXPECT_EQ(summary.family_signature,
              "MCCInteractor+PhotoIonization+SpinningSpacecraftTrajectory");
    EXPECT_TRUE(summary.executed);
}
