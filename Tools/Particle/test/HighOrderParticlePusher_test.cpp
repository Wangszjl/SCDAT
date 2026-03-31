#include "../include/HighOrderParticlePusher.h"
#include "../include/PerThreadChargeBuffer.h"

#include <gtest/gtest.h>

namespace SCDAT
{
namespace Particle
{

TEST(HighOrderParticlePusherTest, RungeKuttaPushesParticleForward)
{
    Particle particle(1, ParticleType::ELECTRON, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                      SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);

    HighOrderParticlePusher pusher(0.5);
    pusher.pushParticle(
        &particle, [](const SCDAT::Geometry::Point3D&) { return SCDAT::Geometry::Vector3D(2.0, 0.0, 0.0); },
        [](const SCDAT::Geometry::Point3D&) { return SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0); });

    EXPECT_NEAR(particle.getVelocity().x(), 1.0, 1e-12);
    EXPECT_NEAR(particle.getPosition().x(), 0.25, 1e-12);
    EXPECT_EQ(pusher.getStatistics().pushed_particles, 1u);
}

TEST(HighOrderParticlePusherTest, PerThreadChargeBufferReducesAcrossThreads)
{
    PerThreadChargeBuffer buffer(3, 2);

    buffer.addCharge(0, 0, 1.5);
    buffer.addCharge(0, 2, 0.5);
    buffer.addCharge(1, 0, 2.0);
    buffer.addCharge(1, 1, -1.0);

    const auto reduced = buffer.reduce();

    ASSERT_EQ(reduced.size(), 3u);
    EXPECT_DOUBLE_EQ(reduced[0], 3.5);
    EXPECT_DOUBLE_EQ(reduced[1], -1.0);
    EXPECT_DOUBLE_EQ(reduced[2], 0.5);
}

} // namespace Particle
} // namespace SCDAT
