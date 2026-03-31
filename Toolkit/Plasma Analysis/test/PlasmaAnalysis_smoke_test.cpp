#include "PICFluidIntegration.h"
#include "PlasmaAnalysisCases.h"

#include <gtest/gtest.h>

#include <filesystem>

using SCDAT::Toolkit::PlasmaAnalysis::PICFluidIntegration;
using SCDAT::Toolkit::PlasmaAnalysis::PlasmaScenarioPreset;

TEST(PlasmaAnalysisSmokeTest, ToolkitInitializesAdvancesAndExports)
{
    PICFluidIntegration integration;
    PlasmaScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::PlasmaAnalysis::tryGetPlasmaScenarioPreset("hall_thruster_plume", preset));

    ASSERT_TRUE(integration.initialize(preset.config));
    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_TRUE(integration.advance(preset.config.time_step_s));
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "plasma_analysis_smoke.csv";
    ASSERT_TRUE(integration.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));
    EXPECT_GE(integration.getStatus().steps_completed, 5u);
}
