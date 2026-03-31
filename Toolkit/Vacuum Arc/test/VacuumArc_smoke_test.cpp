#include "SurfaceDischargeArcAlgorithm.h"
#include "VacuumArcCases.h"

#include <gtest/gtest.h>

#include <filesystem>

using SCDAT::Toolkit::VacuumArc::SurfaceDischargeArcAlgorithm;
using SCDAT::Toolkit::VacuumArc::VacuumArcScenarioPreset;

TEST(VacuumArcSmokeTest, ToolkitInitializesAdvancesAndExports)
{
    SurfaceDischargeArcAlgorithm algorithm;
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("triple_junction_flashover", preset));

    ASSERT_TRUE(algorithm.initialize(preset.config));
    for (int i = 0; i < 6; ++i)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "vacuum_arc_smoke.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));
    EXPECT_TRUE(algorithm.getStatus().discharge_active);
}
