#include "SurfaceScenarioCatalog.h"
#include "SurfaceSimulationRunner.h"
#include "SurfaceRuntimePlan.h"
#include "SurfaceTransitionEngine.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

using SCDAT::Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset;
using SCDAT::Toolkit::SurfaceCharging::SurfaceScenarioCatalog;
using SCDAT::Toolkit::SurfaceCharging::SurfaceSimulationRunner;
using SCDAT::Toolkit::SurfaceCharging::SurfaceAdvanceTransitionInput;
using SCDAT::Toolkit::SurfaceCharging::evaluateSurfaceAdvanceTransition;
using SCDAT::Toolkit::SurfaceCharging::compileSurfaceRuntimePlan;

namespace
{

TEST(SurfaceScenarioCatalogTest, ListsPresetsAndResolvesAliases)
{
    const SurfaceScenarioCatalog catalog;
    const auto mainline_names = catalog.listMainlinePresetNames();
    const auto replay_names = catalog.listReplayPresetNames();

    EXPECT_FALSE(mainline_names.empty());
    EXPECT_FALSE(replay_names.empty());
    EXPECT_NE(std::find(mainline_names.begin(), mainline_names.end(), "geo_ecss_kapton_ref"),
              mainline_names.end());
    EXPECT_NE(std::find(replay_names.begin(), replay_names.end(),
                        "geo_ecss_kapton_ref_legacy_compatible"),
              replay_names.end());

    SurfaceChargingScenarioPreset alias_preset;
    EXPECT_TRUE(catalog.tryGetMainlinePreset("leo_daylight_kapton", alias_preset));
    EXPECT_EQ(alias_preset.name, "leo_pic_circuit_ram_facing");

    SurfaceChargingScenarioPreset missing_preset;
    EXPECT_FALSE(catalog.tryGetPreset("surface_preset_does_not_exist", missing_preset));
}

TEST(SurfaceSimulationRunnerTest, RunsDefaultPresetWithZeroSteps)
{
    const SurfaceScenarioCatalog catalog;
    auto preset = catalog.makeDefaultPreset();
    preset.steps = 0;
    preset.adaptive_time_stepping = false;

    const std::filesystem::path output_csv =
        std::filesystem::temp_directory_path() / "scdat_surface_object_layer_runner_test.csv";
    std::error_code ec;
    std::filesystem::remove(output_csv, ec);

    const SurfaceSimulationRunner runner;
    const auto result = runner.run(preset, output_csv);

    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_TRUE(std::filesystem::exists(output_csv));

    std::filesystem::remove(output_csv, ec);
}

TEST(SurfaceRuntimePlanCompilerTest, CompilesPresetMetadataAndSolverPolicy)
{
    const SurfaceScenarioCatalog catalog;
    auto preset = catalog.makeDefaultPreset();
    preset.name = "runtime_plan_test_case";
    preset.description = "object-layer runtime plan coverage";
    preset.config.solver_config.coupling_mode = "field_particle_implicit";
    preset.config.solver_config.convergence_policy = "residual_norm_guarded";
    preset.config.solver_config.max_iterations = 17;
    preset.config.solver_config.residual_tolerance = 5.0e-4;

    const auto plan = compileSurfaceRuntimePlan(preset);

    EXPECT_EQ(plan.source_name, preset.name);
    EXPECT_EQ(plan.source_description, preset.description);
    EXPECT_TRUE(plan.solver_policy_flags.implicit_coupling_requested);
    EXPECT_TRUE(plan.solver_policy_flags.residual_guarded_requested);
    EXPECT_GE(plan.compiled_config.internal_substeps, 24u);
    EXPECT_DOUBLE_EQ(
        plan.compiled_config.material.getScalarProperty("shared_surface_global_coupled_iterations",
                                                        0.0),
        17.0);
    EXPECT_DOUBLE_EQ(
        plan.compiled_config.material.getScalarProperty("shared_surface_global_coupled_tolerance_v",
                                                        0.0),
        5.0e-4);
}

TEST(SurfaceTransitionEngineTest, ResolvesLegacyAndSharedRuntimeDecisions)
{
    SurfaceAdvanceTransitionInput input;
    input.has_legacy_benchmark_replay = true;
    input.has_circuit_model = true;
    input.patch_count = 3;
    input.config.runtime_route = SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark;
    input.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ReplayFromReference;
    input.config.surface_pic_runtime_kind =
        SCDAT::Toolkit::SurfaceCharging::SurfacePicRuntimeKind::GraphCoupledSharedSurface;
    input.config.enable_live_pic_window = true;
    input.config.enable_pic_calibration = true;
    input.config.regime = SCDAT::Toolkit::SurfaceCharging::SurfaceChargingRegime::GeoKineticPicLike;
    input.config.enable_body_patch_circuit = true;
    input.config.use_reference_current_balance = true;
    input.config.pic_recalibration_interval_steps = 4;
    input.config.solver_config.coupling_mode = "field_particle_implicit";
    input.config.solver_config.convergence_policy = "fixed_iteration";
    input.config.solver_config.max_iterations = 5;
    input.config.material.setScalarProperty("shared_surface_global_coupled_iterations", 5.0);
    input.config.material.setScalarProperty("shared_surface_live_pic_coupled_refresh", 1.0);
    input.config.material.setScalarProperty("shared_surface_particle_transport_weight", 0.25);
    input.status.steps_completed = 4;

    const auto transition = evaluateSurfaceAdvanceTransition(input);

    EXPECT_TRUE(transition.use_legacy_benchmark_replay);
    EXPECT_TRUE(transition.use_reference_circuit);
    EXPECT_TRUE(transition.shared_surface_pic_runtime);
    EXPECT_TRUE(transition.shared_global_coupled_solve);
    EXPECT_TRUE(transition.shared_live_pic_coupled_refresh);
    EXPECT_TRUE(transition.shared_particle_transport_coupling);
    EXPECT_TRUE(transition.pic_recalibration_requested);
    EXPECT_TRUE(transition.fixed_iteration_policy);
    EXPECT_FALSE(transition.residual_guarded_policy);
    EXPECT_GE(transition.shared_global_coupled_iteration_limit, 5u);
}

TEST(SurfaceTransitionEngineTest, LegacyReplayRequiresLegacyRouteAndReplayMode)
{
    SurfaceAdvanceTransitionInput input;
    input.has_legacy_benchmark_replay = true;
    input.config.runtime_route = SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark;
    input.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;

    auto transition = evaluateSurfaceAdvanceTransition(input);
    EXPECT_FALSE(transition.use_legacy_benchmark_replay);

    input.config.runtime_route = SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::SCDATUnified;
    input.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ReplayFromReference;
    transition = evaluateSurfaceAdvanceTransition(input);
    EXPECT_FALSE(transition.use_legacy_benchmark_replay);

    input.config.runtime_route = SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark;
    input.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ReplayFromReference;
    transition = evaluateSurfaceAdvanceTransition(input);
    EXPECT_TRUE(transition.use_legacy_benchmark_replay);
}

TEST(SurfaceTransitionEngineTest, AppliesLocalTimeSpinningAndSunFluxEventChain)
{
    SurfaceAdvanceTransitionInput input;
    input.config.enable_pic_calibration = true;
    input.config.pic_recalibration_interval_steps = 4;
    input.status.steps_completed = 4;

    input.config.material.setScalarProperty("transition_local_time_enabled", 1.0);
    input.config.material.setScalarProperty("transition_local_time_base_hour", 0.0);
    input.config.material.setScalarProperty("transition_local_time_daylight_start_hour", 6.0);
    input.config.material.setScalarProperty("transition_local_time_daylight_end_hour", 18.0);

    input.config.material.setScalarProperty("transition_spinning_enabled", 1.0);
    input.config.material.setScalarProperty("transition_spin_period_s", 24.0 * 3600.0);

    input.config.material.setScalarProperty("transition_sun_flux_enabled", 1.0);
    input.config.material.setScalarProperty("transition_sun_flux_base_scale", 1.2);
    input.config.material.setScalarProperty("transition_sun_flux_night_scale", 0.1);
    input.config.material.setScalarProperty("transition_sun_flux_spin_amplitude", 0.2);
    input.config.material.setScalarProperty("transition_sun_flux_pic_recalibration_min_scale",
                                            0.5);

    input.status.time_s = 6.0 * 3600.0;
    const auto daytime_transition = evaluateSurfaceAdvanceTransition(input);
    EXPECT_TRUE(daytime_transition.local_time_transition_active);
    EXPECT_TRUE(daytime_transition.spinning_spacecraft_active);
    EXPECT_TRUE(daytime_transition.sun_flux_updater_active);
    EXPECT_NEAR(daytime_transition.local_time_hour, 6.0, 1.0e-9);
    EXPECT_GT(daytime_transition.sun_flux_scale, 1.0);
    EXPECT_TRUE(daytime_transition.pic_recalibration_requested);

    input.status.time_s = 23.0 * 3600.0;
    const auto nighttime_transition = evaluateSurfaceAdvanceTransition(input);
    EXPECT_FALSE(nighttime_transition.local_time_transition_active);
    EXPECT_TRUE(nighttime_transition.spinning_spacecraft_active);
    EXPECT_TRUE(nighttime_transition.sun_flux_updater_active);
    EXPECT_LT(nighttime_transition.sun_flux_scale, 0.5);
    EXPECT_LT(nighttime_transition.sun_flux_scale, daytime_transition.sun_flux_scale);
    EXPECT_FALSE(nighttime_transition.pic_recalibration_requested);
}

TEST(SurfaceTransitionEngineTest, EvaluatesExtendedEventBundleWhenSunFluxIsDisabled)
{
    SurfaceAdvanceTransitionInput input;
    input.config.material.setScalarProperty("transition_sun_flux_enabled", 0.0);
    input.config.material.setScalarProperty("transition_conductivity_evolution_enabled", 1.0);
    input.config.material.setScalarProperty("transition_source_flux_updater_enabled", 1.0);
    input.config.material.setScalarProperty("transition_simulation_param_updater_enabled", 1.0);
    input.config.material.setScalarProperty(
        "transition_sheath_or_presheath_poisson_bc_updater_enabled", 1.0);
    input.config.material.setScalarProperty("transition_rccabs_sc_updater_enabled", 1.0);
    input.config.material.setScalarProperty("transition_vcross_bfield_updater_enabled", 1.0);
    input.config.material.setScalarProperty("transition_basic_eclipse_exit_enabled", 1.0);
    input.config.material.setScalarProperty("transition_transient_artificial_sources_enabled",
                                            1.0);
    input.config.material.setScalarProperty("transition_langmuir_probe_transition_enabled", 1.0);

    const auto transition = evaluateSurfaceAdvanceTransition(input);

    EXPECT_FALSE(transition.sun_flux_updater_active);
    EXPECT_TRUE(transition.extended_events.conductivity_evolution_active);
    EXPECT_TRUE(transition.extended_events.source_flux_updater_active);
    EXPECT_TRUE(transition.extended_events.simulation_param_updater_active);
    EXPECT_TRUE(transition.extended_events.sheath_or_presheath_poisson_bc_updater_active);
    EXPECT_TRUE(transition.extended_events.rccabs_sc_updater_active);
    EXPECT_TRUE(transition.extended_events.vcross_bfield_updater_active);
    EXPECT_TRUE(transition.extended_events.basic_eclipse_exit_active);
    EXPECT_TRUE(transition.extended_events.transient_artificial_sources_active);
    EXPECT_TRUE(transition.extended_events.langmuir_probe_transition_active);
}

} // namespace
