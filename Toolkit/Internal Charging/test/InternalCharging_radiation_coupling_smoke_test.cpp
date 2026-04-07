#include "InternalChargingCases.h"
#include "CouplingBase.h"
#include "MultiPhysicsManager.h"
#include "RadiationCases.h"
#include "RadiationDoseAlgorithm.h"
#include "SpacecraftInternalChargingAlgorithm.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

using SCDAT::Toolkit::InternalCharging::InternalChargingRadiationDrive;
using SCDAT::Toolkit::InternalCharging::InternalChargingScenarioPreset;
using SCDAT::Toolkit::InternalCharging::InternalChargingSourceMode;
using SCDAT::Toolkit::InternalCharging::SpacecraftInternalChargingAlgorithm;
using SCDAT::Toolkit::Radiation::RadiationDoseAlgorithm;
using SCDAT::Coupling::ConvergenceCriteria;
using SCDAT::Coupling::CouplingType;
using SCDAT::Coupling::FunctionalCoupling;
using SCDAT::Coupling::MultiPhysicsManager;
using SCDAT::Toolkit::Radiation::RadiationScenarioPreset;

namespace
{

double mapDepositedEnergyToCurrentDensity(double deposited_energy_j_per_m2, double dt,
                                          double mean_energy_ev,
                                          double incident_charge_state_abs)
{
    if (deposited_energy_j_per_m2 <= 0.0 || dt <= 0.0 || mean_energy_ev <= 0.0 ||
        incident_charge_state_abs <= 0.0)
    {
        return 0.0;
    }

    constexpr double kEvToJ = 1.602176634e-19;
    constexpr double kElementaryChargeC = 1.602176634e-19;
    const double particle_energy_j = mean_energy_ev * kEvToJ;
    const double particle_flux_m2_s = deposited_energy_j_per_m2 / (dt * particle_energy_j);
    const double particle_charge_c = incident_charge_state_abs * kElementaryChargeC;
    return particle_flux_m2_s * particle_charge_c;
}

double computeInternalToRadiationFeedbackScale(
    const SCDAT::Toolkit::InternalCharging::InternalChargingStatus& internal_status)
{
    const double field_shielding_term =
        1.0 / (1.0 + std::abs(internal_status.max_electric_field_v_per_m) / 1.0e8);
    const double conductivity_relief_term =
        1.0 + std::clamp(internal_status.effective_conductivity_s_per_m / 5.0e-8, 0.0, 0.15);
    return std::clamp(field_shielding_term * conductivity_relief_term, 0.85, 1.15);
}

} // namespace

TEST(InternalChargingRadiationCouplingSmokeTest, RadiationModeRequiresDriveBeforeAdvance)
{
    InternalChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::InternalCharging::tryGetInternalChargingScenarioPreset(
        "geo_electron_belt", preset));
    preset.config.source_mode = InternalChargingSourceMode::Radiation;

    SpacecraftInternalChargingAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));
    EXPECT_FALSE(algorithm.advance(preset.time_step_s));
}

TEST(InternalChargingRadiationCouplingSmokeTest, OnlineOneWayCouplingDrivesInternalState)
{
    InternalChargingScenarioPreset internal_preset;
    RadiationScenarioPreset radiation_preset;
    ASSERT_TRUE(SCDAT::Toolkit::InternalCharging::tryGetInternalChargingScenarioPreset(
        "geo_electron_belt", internal_preset));
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", radiation_preset));

    internal_preset.config.source_mode = InternalChargingSourceMode::Radiation;

    SpacecraftInternalChargingAlgorithm internal_algorithm;
    RadiationDoseAlgorithm radiation_algorithm;
    ASSERT_TRUE(internal_algorithm.initialize(internal_preset.config));
    ASSERT_TRUE(radiation_algorithm.initialize(radiation_preset.config));

    const double dt = internal_preset.time_step_s;
    ASSERT_GT(dt, 0.0);

    double deposited_energy_previous_j_per_m2 = 0.0;
    const std::size_t coupling_steps = std::min<std::size_t>(internal_preset.steps, 8);
    for (std::size_t i = 0; i < coupling_steps; ++i)
    {
        ASSERT_TRUE(radiation_algorithm.advance(dt));
        const auto& radiation_status = radiation_algorithm.getStatus();
        const double deposited_energy_delta_j_per_m2 =
            std::max(0.0, radiation_status.deposited_energy_j_per_m2 - deposited_energy_previous_j_per_m2);
        deposited_energy_previous_j_per_m2 = radiation_status.deposited_energy_j_per_m2;

        InternalChargingRadiationDrive drive;
        drive.incident_energy_ev = std::max(1.0, radiation_preset.config.mean_energy_ev);
        drive.incident_charge_state_abs = 1.0;
        drive.incident_current_density_a_per_m2 =
            mapDepositedEnergyToCurrentDensity(deposited_energy_delta_j_per_m2, dt,
                                               drive.incident_energy_ev,
                                               drive.incident_charge_state_abs);
        internal_algorithm.setRadiationDrive(drive);

        ASSERT_TRUE(internal_algorithm.advance(dt));
    }

    const auto& status = internal_algorithm.getStatus();
    EXPECT_GT(status.time_s, 0.0);
    EXPECT_GT(status.total_stored_energy_j, 0.0);
    EXPECT_GT(status.max_electric_field_v_per_m, 0.0);
    EXPECT_GT(status.average_dose_gy, 0.0);
    EXPECT_GT(status.effective_conductivity_s_per_m, 0.0);
    EXPECT_GT(status.deposited_charge_rate_c_per_m3_s, 0.0);
    EXPECT_GT(status.deposited_particle_rate_per_m3_s, 0.0);
    EXPECT_NEAR(status.effective_incident_charge_state_abs, 1.0, 1.0e-12);

    ASSERT_FALSE(status.volume_charge_density_c_per_m3.empty());
    ASSERT_FALSE(status.electric_field_v_per_m.empty());
    EXPECT_TRUE(std::all_of(status.volume_charge_density_c_per_m3.begin(),
                            status.volume_charge_density_c_per_m3.end(),
                            [](double value) { return std::isfinite(value); }));
    EXPECT_TRUE(std::all_of(status.electric_field_v_per_m.begin(),
                            status.electric_field_v_per_m.end(),
                            [](double value) { return std::isfinite(value); }));
}

TEST(InternalChargingRadiationCouplingSmokeTest, IterativeSubstepCouplingProvidesConvergenceHistory)
{
    constexpr std::size_t kCouplingSubstepsPerStep = 3;
    constexpr int kCouplingIterationsPerSubstep = 4;
    constexpr double kCouplingRelativeTolerance = 2.0;
    constexpr double kCouplingAbsoluteTolerance = 1.0e-3;

    InternalChargingScenarioPreset internal_preset;
    RadiationScenarioPreset radiation_preset;
    ASSERT_TRUE(SCDAT::Toolkit::InternalCharging::tryGetInternalChargingScenarioPreset(
        "geo_electron_belt", internal_preset));
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", radiation_preset));

    internal_preset.config.source_mode = InternalChargingSourceMode::Radiation;

    SpacecraftInternalChargingAlgorithm internal_algorithm;
    RadiationDoseAlgorithm radiation_algorithm;
    ASSERT_TRUE(internal_algorithm.initialize(internal_preset.config));
    ASSERT_TRUE(radiation_algorithm.initialize(radiation_preset.config));

    const double dt = internal_preset.time_step_s;
    ASSERT_GT(dt, 0.0);
    const double substep_dt = dt / static_cast<double>(kCouplingSubstepsPerStep);
    const double iteration_dt = substep_dt / static_cast<double>(kCouplingIterationsPerSubstep);
    ASSERT_GT(iteration_dt, 0.0);

    struct CouplingHistoryRow
    {
        std::size_t macro_step = 0;
        std::size_t substep = 0;
        int iteration = 0;
        double relative_error = 0.0;
        double absolute_error = 0.0;
        bool converged = false;
    };

    double deposited_energy_previous_j_per_m2 = 0.0;
    double previous_drive_current_density_a_per_m2 = 0.0;
    double current_relative_error = 0.0;
    double current_absolute_error = 0.0;
    std::size_t active_macro_step = 0;
    std::size_t active_substep = 0;
    int active_iteration = 0;
    std::vector<CouplingHistoryRow> history;

    auto coupling = std::make_shared<FunctionalCoupling>(
        "iterative_internal_radiation_test", CouplingType::ITERATIVE,
        [&](double step_dt) {
            ++active_iteration;
            if (!radiation_algorithm.advance(step_dt))
            {
                return SCDAT::VoidResult::failure("Radiation advance failed in iterative test coupling");
            }

            const auto& radiation_status = radiation_algorithm.getStatus();
            const double deposited_energy_delta_j_per_m2 =
                std::max(0.0, radiation_status.deposited_energy_j_per_m2 - deposited_energy_previous_j_per_m2);
            deposited_energy_previous_j_per_m2 = radiation_status.deposited_energy_j_per_m2;

            const auto& internal_status_before = internal_algorithm.getStatus();
            const double feedback_scale =
                computeInternalToRadiationFeedbackScale(internal_status_before);

            InternalChargingRadiationDrive drive;
            drive.incident_energy_ev = std::max(1.0, radiation_preset.config.mean_energy_ev);
            drive.incident_charge_state_abs = 1.0;
            const double base_drive_current_density_a_per_m2 =
                mapDepositedEnergyToCurrentDensity(deposited_energy_delta_j_per_m2, step_dt,
                                                   drive.incident_energy_ev,
                                                   drive.incident_charge_state_abs);
            const double drive_current_density_a_per_m2 =
                base_drive_current_density_a_per_m2 * feedback_scale;
            drive.incident_current_density_a_per_m2 = drive_current_density_a_per_m2;
            internal_algorithm.setRadiationDrive(drive);

            if (!internal_algorithm.advance(step_dt))
            {
                return SCDAT::VoidResult::failure(
                    "Internal charging advance failed in iterative test coupling");
            }

            current_absolute_error =
                std::abs(drive_current_density_a_per_m2 - previous_drive_current_density_a_per_m2);
            const double relative_denominator = std::max(
                std::max(std::abs(previous_drive_current_density_a_per_m2),
                         std::abs(drive_current_density_a_per_m2)),
                1.0e-12);
            current_relative_error = current_absolute_error / relative_denominator;
            previous_drive_current_density_a_per_m2 = drive_current_density_a_per_m2;

            history.push_back(CouplingHistoryRow{
                active_macro_step,
                active_substep,
                active_iteration,
                current_relative_error,
                current_absolute_error,
                current_relative_error <= kCouplingRelativeTolerance &&
                    current_absolute_error <= kCouplingAbsoluteTolerance,
            });

            return SCDAT::VoidResult::success();
        },
        [&]() { return std::pair<double, double>{current_relative_error, current_absolute_error}; });

    ConvergenceCriteria criteria;
    criteria.relative_tolerance = kCouplingRelativeTolerance;
    criteria.absolute_tolerance = kCouplingAbsoluteTolerance;
    criteria.max_iterations = kCouplingIterationsPerSubstep;
    criteria.min_iterations = kCouplingIterationsPerSubstep;
    coupling->setConvergenceCriteria(criteria);

    MultiPhysicsManager manager;
    manager.addCoupling(coupling);

    const std::size_t coupling_steps = std::min<std::size_t>(internal_preset.steps, 4);
    for (std::size_t macro_step = 0; macro_step < coupling_steps; ++macro_step)
    {
        for (std::size_t substep = 0; substep < kCouplingSubstepsPerStep; ++substep)
        {
            active_macro_step = macro_step;
            active_substep = substep;
            active_iteration = 0;

            ASSERT_TRUE(manager.executeIterativeCouplings(iteration_dt));
            const auto stats = manager.getAllStatistics();
            const auto stats_it = stats.find("iterative_internal_radiation_test");
            ASSERT_TRUE(stats_it != stats.end());
            EXPECT_TRUE(stats_it->second.converged);
            EXPECT_EQ(stats_it->second.total_iterations, kCouplingIterationsPerSubstep);
        }
    }

    ASSERT_EQ(history.size(), coupling_steps * kCouplingSubstepsPerStep *
                                  static_cast<std::size_t>(kCouplingIterationsPerSubstep));
    EXPECT_TRUE(std::all_of(history.begin(), history.end(),
                            [](const CouplingHistoryRow& row) {
                                return std::isfinite(row.relative_error) &&
                                       std::isfinite(row.absolute_error);
                            }));
    EXPECT_TRUE(std::any_of(history.begin(), history.end(),
                            [](const CouplingHistoryRow& row) { return row.converged; }));

    const auto& status = internal_algorithm.getStatus();
    EXPECT_GT(status.time_s, 0.0);
    EXPECT_GT(status.total_stored_energy_j, 0.0);
    EXPECT_GT(status.max_electric_field_v_per_m, 0.0);
    EXPECT_GT(status.average_dose_gy, 0.0);
    EXPECT_GT(status.effective_conductivity_s_per_m, 0.0);
    EXPECT_GT(status.deposited_charge_rate_c_per_m3_s, 0.0);
    EXPECT_GT(status.deposited_particle_rate_per_m3_s, 0.0);
}
