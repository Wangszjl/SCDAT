#include "DensePlasmaSurfaceCharging.h"
#include "SurfaceChargingCases.h"

#include <gtest/gtest.h>

#include <fstream>
#include <cmath>
#include <filesystem>
#include <sstream>

using SCDAT::Toolkit::SurfaceCharging::DensePlasmaSurfaceCharging;
using SCDAT::Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset;

TEST(SurfaceChargingSmokeTest, ToolkitInitializesAdvancesAndExports)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_eclipse_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    for (int i = 0; i < 8; ++i)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "surface_charging_smoke.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    EXPECT_NE(header.find("floating_equilibrium_potential_v"), std::string::npos);
    EXPECT_NE(header.find("net_current_density_a_per_m2"), std::string::npos);
    EXPECT_NE(header.find("effective_conductivity_s_per_m"), std::string::npos);
    EXPECT_NE(header.find("adaptive_time_step_s"), std::string::npos);
    EXPECT_NE(header.find("electron_pic_calibration_factor"), std::string::npos);

    EXPECT_GE(charging.getStatus().steps_completed, 8u);
    EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
    EXPECT_TRUE(std::isfinite(charging.getStatus().currents.total_current_a_per_m2));
}

TEST(SurfaceChargingSmokeTest, ThrusterPlumePresetRemainsFinite)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "thruster_plume_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    for (int i = 0; i < 12; ++i)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
        EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
        EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_charge_density_c_per_m2));
        EXPECT_TRUE(std::isfinite(charging.getStatus().currents.total_current_a_per_m2));
    }
}

TEST(SurfaceChargingSmokeTest, GeoEclipseFloatingPotentialIsFiniteAndNegative)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_eclipse_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const double floating_potential = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(floating_potential));
    EXPECT_LT(floating_potential, 0.0);
    EXPECT_GT(floating_potential, -5.0e3);
}

TEST(SurfaceChargingSmokeTest, LeoDaylightFloatingPotentialIsModeratelyPositive)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_daylight_kapton", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const double floating_potential = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(floating_potential));
    EXPECT_GT(floating_potential, 0.0);
    EXPECT_LT(floating_potential, 20.0);
}

TEST(SurfaceChargingSmokeTest, LeoWakePotentialExceedsRamFacingPotential)
{
    DensePlasmaSurfaceCharging ram_facing;
    DensePlasmaSurfaceCharging wake_facing;
    SurfaceChargingScenarioPreset ram_preset;
    SurfaceChargingScenarioPreset wake_preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_daylight_kapton", ram_preset));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_daylight_wake_kapton", wake_preset));

    ASSERT_TRUE(ram_facing.initialize(ram_preset.config));
    ASSERT_TRUE(wake_facing.initialize(wake_preset.config));

    const double ram_potential = ram_facing.computeFloatingPotential();
    const double wake_potential = wake_facing.computeFloatingPotential();

    EXPECT_TRUE(std::isfinite(ram_potential));
    EXPECT_TRUE(std::isfinite(wake_potential));
    EXPECT_GT(wake_potential, ram_potential);
}

TEST(SurfaceChargingSmokeTest, PositiveBiasSuppressesPhotoelectronEscape)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_daylight_kapton", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const auto near_zero = charging.computeSurfaceCurrents(0.0);
    const auto positive_bias = charging.computeSurfaceCurrents(2.0);

    EXPECT_LT(positive_bias.photo_emission_a_per_m2, near_zero.photo_emission_a_per_m2);
}

TEST(SurfaceChargingSmokeTest, GeoPicCalibrationProducesFiniteTrajectory)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_eclipse_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path = std::filesystem::temp_directory_path() / "surface_charging_geo_pic.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    EXPECT_NE(header.find("ion_pic_calibration_factor"), std::string::npos);
    EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
}

TEST(SurfaceChargingSmokeTest, ThrusterPlumeFloatingPotentialIsModeratelyNegative)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "thruster_plume_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const double floating_potential = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(floating_potential));
    EXPECT_LT(floating_potential, -1.0);
    EXPECT_GT(floating_potential, -25.0);
}
