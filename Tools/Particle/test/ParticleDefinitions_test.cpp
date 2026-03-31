/**
 * @file ParticleDefinitions_test.cpp
 * @brief ParticleDefinitions / BorisAlgorithm 单元测试
 * @details
 * 覆盖粒子定义模块的关键路径：
 * - 粒子工厂参数映射（电子/离子）；
 * - 生命周期状态迁移；
 * - Boris 推进器在纯电场和相对论开关下的行为；
 * - toString 扩展字段输出。
 */

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
    // 电子工厂应使用全局物理常量中的电子质量/电荷。
    Particle p = ParticleFactory::createElectron(1, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
                                                 SCDAT::Geometry::Vector3D(1.0, 0.0, 0.0), 2.0);

    EXPECT_EQ(p.getType(), ParticleType::ELECTRON);
    EXPECT_DOUBLE_EQ(p.getMass(), PhysicsConstants::ElectronMass);
    EXPECT_DOUBLE_EQ(p.getCharge(), PhysicsConstants::ElectronCharge);
    EXPECT_DOUBLE_EQ(p.getWeight(), 2.0);
}

TEST(ParticleDefinitionsTest, FactoryCreateIonAssignsSpecies)
{
    // 检查离子种类信息、质量数/电荷数与派生质量、电荷值一致性。
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
    // 验证年龄累加与多个终态标记 API 行为。
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
    // 纯电场 + 无磁场场景：速度应沿电场方向线性增长。
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
    // 相对论路径在零场中应保持方向，速度不超过光速。
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
    // 字符串输出中应包含光电子扩展字段，便于诊断与日志追踪。
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
