#include "../../Interactions/Collisions/include/CollisionPicAdapter.h"
#include "../../Particle/include/ParticleManager.h"

#include <gtest/gtest.h>

namespace
{

TEST(PICcoreCollisionRuntimeTest, SharedCollisionRuntimeSupportsToolkitCrossSectionSets)
{
    SCDAT::Particle::ParticleManager particle_manager;
    SCDAT::Collision::MonteCarloCollisionHandler handler(12345);

    SCDAT::Collision::BackgroundMccRuntimeConfig runtime_config;
    runtime_config.floor_densities.neutral_density_m3 = 2.4e20;
    runtime_config.floor_densities.ion_density_m3 = 1.4e14;
    runtime_config.floor_densities.electron_density_m3 = 1.4e14;
    runtime_config.cross_section_set_id =
        SCDAT::Collision::CollisionCrossSectionSetId::SurfacePicV1;
    runtime_config.reinitialize_configuration = true;

    const auto result = SCDAT::Collision::executeBackgroundMccStep(
        particle_manager.getContainer(), handler, 1.0e-10, 1.0e-6, runtime_config);

    EXPECT_EQ(result.active_particles, 0u);
    EXPECT_DOUBLE_EQ(result.effective_neutral_density_m3,
                     runtime_config.floor_densities.neutral_density_m3);
    EXPECT_DOUBLE_EQ(result.effective_ion_density_m3,
                     runtime_config.floor_densities.ion_density_m3);
    EXPECT_DOUBLE_EQ(result.effective_electron_density_m3,
                     runtime_config.floor_densities.electron_density_m3);

    const auto parsed_surface =
        SCDAT::Collision::parseCollisionCrossSectionSetId("surface_pic_v1");
    ASSERT_TRUE(parsed_surface.has_value());
    EXPECT_EQ(*parsed_surface, SCDAT::Collision::CollisionCrossSectionSetId::SurfacePicV1);

    const auto parsed_arc =
        SCDAT::Collision::parseCollisionCrossSectionSetId("vacuum_arc_pic_v1");
    ASSERT_TRUE(parsed_arc.has_value());
    EXPECT_EQ(*parsed_arc, SCDAT::Collision::CollisionCrossSectionSetId::VacuumArcPicV1);
}

} // namespace


