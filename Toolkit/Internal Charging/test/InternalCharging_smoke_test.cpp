#include "InternalChargingCases.h"
#include "SpacecraftInternalChargingAlgorithm.h"

#include <gtest/gtest.h>

#include <filesystem>

using SCDAT::Toolkit::InternalCharging::InternalChargingScenarioPreset;
using SCDAT::Toolkit::InternalCharging::SpacecraftInternalChargingAlgorithm;

TEST(InternalChargingSmokeTest, ToolkitInitializesAdvancesAndExports)
{
    SpacecraftInternalChargingAlgorithm algorithm;
    InternalChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::InternalCharging::tryGetInternalChargingScenarioPreset(
        "geo_electron_belt", preset));

    ASSERT_TRUE(algorithm.initialize(preset.config));
    for (int i = 0; i < 4; ++i)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "internal_charging_smoke.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));
    EXPECT_EQ(algorithm.getStatus().volume_charge_density_c_per_m3.size(), preset.config.layers);
}
