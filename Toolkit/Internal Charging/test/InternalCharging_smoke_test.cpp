#include "InternalChargingCases.h"
#include "ParticleTransportModel.h"
#include "SpacecraftInternalChargingAlgorithm.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>

using SCDAT::Toolkit::InternalCharging::InternalChargingScenarioPreset;
using SCDAT::Toolkit::InternalCharging::ParticleTransportModel;
using SCDAT::Toolkit::InternalCharging::SpacecraftInternalChargingAlgorithm;

namespace
{
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kKaptonPermittivity = 3.4;
constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kEvToJ = 1.602176634e-19;
}

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

    const auto& status = algorithm.getStatus();
    EXPECT_EQ(status.volume_charge_density_c_per_m3.size(), preset.config.layers);
    EXPECT_GT(status.max_electric_field_v_per_m, 0.0);
    EXPECT_GT(status.average_dose_gy, 0.0);
    EXPECT_GT(status.effective_conductivity_s_per_m, 0.0);
    EXPECT_GT(status.deposited_charge_rate_c_per_m3_s, 0.0);
    EXPECT_GT(status.deposited_particle_rate_per_m3_s, 0.0);
    EXPECT_NEAR(status.effective_incident_charge_state_abs, 1.0, 1.0e-12);
}

TEST(InternalChargingSmokeTest, TransportDepositsConservedChargeAndEnergyPerArea)
{
    ParticleTransportModel model;
    constexpr double thickness_m = 4.0e-3;
    constexpr std::size_t layers = 16;
    constexpr double current_density_a_per_m2 = 1.0;
    constexpr double mean_energy_ev = 5.0e3;
    constexpr double dt = 2.5e-1;

    ASSERT_TRUE(model.initialize(thickness_m, layers));
    model.depositFlux(current_density_a_per_m2, mean_energy_ev, 1.0, dt);

    const double layer_thickness_m = thickness_m / static_cast<double>(layers);
    double total_charge_c_per_m2 = 0.0;
    double total_energy_j_per_m2 = 0.0;
    for (const auto& layer : model.getLayerStates())
    {
        total_charge_c_per_m2 += layer.charge_density_c_per_m3 * layer_thickness_m;
        total_energy_j_per_m2 += layer.deposited_energy_j_per_m3 * layer_thickness_m;
    }

    const double expected_charge_c_per_m2 = current_density_a_per_m2 * dt;
    const double expected_energy_j_per_m2 =
        (expected_charge_c_per_m2 / kElementaryChargeC) * (mean_energy_ev * kEvToJ);

    EXPECT_NEAR(total_charge_c_per_m2, expected_charge_c_per_m2, 1.0e-12);
    EXPECT_NEAR(total_energy_j_per_m2, expected_energy_j_per_m2, 1.0e-9);
}

TEST(InternalChargingSmokeTest, TransportChargeStateChangesParticleRateButKeepsChargeRate)
{
    ParticleTransportModel model_z1;
    ParticleTransportModel model_z2;
    constexpr double thickness_m = 2.0e-3;
    constexpr std::size_t layers = 10;
    constexpr double current_density_a_per_m2 = 2.0e-2;
    constexpr double mean_energy_ev = 2.0e4;
    constexpr double dt = 1.0e-3;

    ASSERT_TRUE(model_z1.initialize(thickness_m, layers));
    ASSERT_TRUE(model_z2.initialize(thickness_m, layers));

    model_z1.depositFlux(current_density_a_per_m2, mean_energy_ev, 1.0, dt);
    model_z2.depositFlux(current_density_a_per_m2, mean_energy_ev, 2.0, dt);

    const double layer_thickness_m = thickness_m / static_cast<double>(layers);
    double charge_rate_z1 = 0.0;
    double charge_rate_z2 = 0.0;
    double particle_rate_z1 = 0.0;
    double particle_rate_z2 = 0.0;
    for (std::size_t i = 0; i < layers; ++i)
    {
        const auto& l1 = model_z1.getLayerStates()[i];
        const auto& l2 = model_z2.getLayerStates()[i];
        charge_rate_z1 += l1.charge_density_c_per_m3 * layer_thickness_m / dt;
        charge_rate_z2 += l2.charge_density_c_per_m3 * layer_thickness_m / dt;
        particle_rate_z1 +=
            (l1.deposited_energy_j_per_m3 * layer_thickness_m / dt) / (mean_energy_ev * kEvToJ);
        particle_rate_z2 +=
            (l2.deposited_energy_j_per_m3 * layer_thickness_m / dt) / (mean_energy_ev * kEvToJ);
    }

    EXPECT_NEAR(charge_rate_z1, charge_rate_z2, std::abs(charge_rate_z1) * 1.0e-12 + 1.0e-18);
    EXPECT_NEAR(particle_rate_z2 / particle_rate_z1, 0.5, 1.0e-6);
}

TEST(InternalChargingSmokeTest, StoredEnergyUsesLayerVolumeIntegration)
{
    SpacecraftInternalChargingAlgorithm algorithm;
    SCDAT::Toolkit::InternalCharging::InternalChargingConfiguration config;
    config.layers = 12;
    config.thickness_m = 6.0e-3;
    config.area_m2 = 2.0e-2;
    config.incident_current_density_a_per_m2 = 3.0e-6;
    config.incident_energy_ev = 3.0e4;
    config.material_name = "kapton";

    ASSERT_TRUE(algorithm.initialize(config));
    ASSERT_TRUE(algorithm.advance(1.0e-4));

    const auto& status = algorithm.getStatus();
    ASSERT_EQ(status.electric_field_v_per_m.size(), config.layers);
    EXPECT_GT(status.average_dose_gy, 0.0);
    EXPECT_GT(status.effective_conductivity_s_per_m, 0.0);
    EXPECT_GT(status.deposited_particle_rate_per_m3_s, 0.0);
    EXPECT_NEAR(status.effective_incident_charge_state_abs, 1.0, 1.0e-12);

    const double layer_thickness_m = config.thickness_m / static_cast<double>(config.layers);
    double expected_total_stored_energy_j = 0.0;
    for (double electric_field_v_per_m : status.electric_field_v_per_m)
    {
        const double energy_density_j_per_m3 =
            0.5 * kEpsilon0 * kKaptonPermittivity * electric_field_v_per_m * electric_field_v_per_m;
        expected_total_stored_energy_j += energy_density_j_per_m3 * config.area_m2 * layer_thickness_m;
    }

    EXPECT_GT(status.total_stored_energy_j, 0.0);
    EXPECT_NEAR(status.total_stored_energy_j, expected_total_stored_energy_j,
                std::max(1.0e-18, std::abs(expected_total_stored_energy_j) * 1.0e-12));
}
