#include "../include/Constants.h"
#include "../include/ErrorHandling.h"
#include "../include/MathUtils.h"
#include "../include/ModelRegistry.h"
#include "../include/StringTokenUtils.h"

#include <gtest/gtest.h>

namespace
{

enum class DemoModel
{
    First,
    Second
};

} // namespace

TEST(BasicTest, MathUtilsSafeDivideAndClampHandleFiniteRanges)
{
    EXPECT_DOUBLE_EQ(SCDAT::Basic::MathUtils::safeDivide(6.0, 2.0), 3.0);
    EXPECT_DOUBLE_EQ(SCDAT::Basic::MathUtils::safeDivide(6.0, 0.0, -1.0), -1.0);
    EXPECT_DOUBLE_EQ(SCDAT::Basic::MathUtils::clamp(5.0, 0.0, 3.0), 3.0);
    EXPECT_NEAR(SCDAT::Basic::MathUtils::smoothstep(0.5), 0.5, 1.0e-12);
}

TEST(BasicTest, StringTokenUtilsNormalizeAndTrimAsciiTokens)
{
    EXPECT_EQ(SCDAT::Basic::trimAscii("  surface util  "), "surface util");
    EXPECT_EQ(SCDAT::Basic::toLowerAscii("SurfacePic"), "surfacepic");
    EXPECT_EQ(SCDAT::Basic::normalizeAlnumToken(" Surface-Util/Matrix "), "surfaceutilmatrix");
}

TEST(BasicTest, StringModelRegistryParsesAliasesAfterNormalization)
{
    const SCDAT::Basic::StringModelRegistry<DemoModel> registry{
        {"surface_pic", DemoModel::First},
        {"surface-hybrid", DemoModel::Second},
    };

    ASSERT_TRUE(registry.tryParse("SurfacePic").has_value());
    EXPECT_EQ(*registry.tryParse("SurfacePic"), DemoModel::First);
    ASSERT_TRUE(registry.tryParse(" surface_hybrid ").has_value());
    EXPECT_EQ(*registry.tryParse(" surface_hybrid "), DemoModel::Second);
    EXPECT_FALSE(registry.tryParse("missing").has_value());
}

TEST(BasicTest, ErrorHandlerBuildsValueAndVoidResults)
{
    const auto success_value = SCDAT::Core::ErrorHandler::makeSuccess<int>(42);
    ASSERT_TRUE(success_value);
    EXPECT_EQ(success_value.value(), 42);

    const auto success_void = SCDAT::Core::ErrorHandler::makeSuccess();
    EXPECT_TRUE(success_void);

    const auto error_value =
        SCDAT::Core::ErrorHandler::makeError<int>(SCDAT::Core::ErrorCode::INVALID_ARGUMENT,
                                                  "invalid-input");
    EXPECT_FALSE(error_value);
    EXPECT_EQ(error_value.error(), SCDAT::Core::ErrorCode::INVALID_ARGUMENT);
    EXPECT_EQ(error_value.message(), "invalid-input");
}

TEST(BasicTest, PhysicsConstantsExposeStableConversions)
{
    using SCDAT::Basic::Constants::MathConstants;
    using SCDAT::Basic::Constants::PhysicsConstants;

    EXPECT_GT(MathConstants::Pi, 3.14);
    EXPECT_NEAR(PhysicsConstants::JouleToEV * PhysicsConstants::EVToJoule, 1.0, 1.0e-12);
    EXPECT_GT(PhysicsConstants::SpeedOfLight, 2.9e8);
    EXPECT_LT(PhysicsConstants::ZeroTolerance, 1.0e-9);
}
