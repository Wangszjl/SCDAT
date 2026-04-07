#include "PICFluidIntegration.h"
#include "FluidAlgorithmAdapter.h"
#include "DensePlasmaBoundaryLayer.h"
#include "PlasmaAdvancedClosureModel.h"
#include "PlasmaReactionCollisionLibrary.h"
#include "FluidPicHybridCouplingInterface.h"
#include "MultiScaleSpatialDecomposer.h"
#include "MultiScaleTimeSynchronizer.h"
#include "PlasmaAnalysisCases.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using SCDAT::Output::ColumnarDataSet;
using SCDAT::Toolkit::PlasmaAnalysis::FluidAlgorithmAdapter;
using SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaBoundaryLayer;
using SCDAT::Toolkit::PlasmaAnalysis::AdvancedClosureConfig;
using SCDAT::Toolkit::PlasmaAnalysis::FluidPicCouplingConfig;
using SCDAT::Toolkit::PlasmaAnalysis::FluidPicHybridCouplingInterface;
using SCDAT::Toolkit::PlasmaAnalysis::FluidPicInterfaceState;
using SCDAT::Toolkit::PlasmaAnalysis::IReactionCollisionProcess;
using SCDAT::Toolkit::PlasmaAnalysis::MultiScaleRegionType;
using SCDAT::Toolkit::PlasmaAnalysis::MultiScaleSpatialDecomposer;
using SCDAT::Toolkit::PlasmaAnalysis::MultiScaleSynchronizationConfig;
using SCDAT::Toolkit::PlasmaAnalysis::MultiScaleTimeSynchronizer;
using SCDAT::Toolkit::PlasmaAnalysis::NonEquilibriumClosureModel;
using SCDAT::Toolkit::PlasmaAnalysis::PICFluidIntegration;
using SCDAT::Toolkit::PlasmaAnalysis::PlasmaAdvancedClosureModel;
using SCDAT::Toolkit::PlasmaAnalysis::PlasmaReactionCollisionLibrary;
using SCDAT::Toolkit::PlasmaAnalysis::PlasmaScenarioPreset;
using SCDAT::Toolkit::PlasmaAnalysis::ReactionCollisionConfig;
using SCDAT::Toolkit::PlasmaAnalysis::TurbulenceClosureModel;

namespace
{

bool containsRegion(const std::vector<MultiScaleRegionType>& regions,
                    MultiScaleRegionType target)
{
    return std::find(regions.begin(), regions.end(), target) != regions.end();
}

ColumnarDataSet makeSyntheticMultiscaleProfile()
{
    ColumnarDataSet data_set;
    data_set.axis_name = "z_m";

    for (int i = 0; i < 12; ++i)
    {
        data_set.axis_values.push_back(static_cast<double>(i) * 1.0e-6);
    }

    data_set.scalar_series["electron_density_m3"] = {
        1.0e17, 1.0e17, 1.0e17, 1.0e17,
        1.2e20, 1.2e20, 1.2e20, 1.2e20,
        8.0e20, 8.0e20, 8.0e20, 8.0e20,
    };
    data_set.scalar_series["electric_field_z_v_per_m"] = {
        60.0, 62.0, 64.0, 66.0,
        120.0, 140.0, 160.0, 180.0,
        400.0, 460.0, 520.0, 580.0,
    };
    data_set.scalar_series["electron_temperature_ev"] = {
        0.8, 0.8, 0.8, 0.8,
        1.4, 1.4, 1.4, 1.4,
        3.0, 3.0, 3.0, 3.0,
    };
    data_set.scalar_series["ion_temperature_ev"] = {
        0.7, 0.7, 0.7, 0.7,
        1.0, 1.0, 1.0, 1.0,
        0.8, 0.8, 0.8, 0.8,
    };

    return data_set;
}

ColumnarDataSet makePerturbedProfile(const ColumnarDataSet& input)
{
    ColumnarDataSet output = input;
    auto& density = output.scalar_series["electron_density_m3"];
    auto& field = output.scalar_series["electric_field_z_v_per_m"];

    for (std::size_t i = 0; i < density.size(); ++i)
    {
        const double phase = static_cast<double>(i + 1);
        density[i] *= 1.0 + 1.0e-3 * std::sin(phase);
        field[i] *= 1.0 + 1.0e-3 * std::cos(phase);
    }
    return output;
}

ColumnarDataSet makeLayeredBenchmarkProfile(std::size_t cell_count,
                                            double dz_m,
                                            double density_bias)
{
    ColumnarDataSet data_set;
    data_set.axis_name = "z_m";

    data_set.axis_values.reserve(cell_count);
    auto& density = data_set.scalar_series["electron_density_m3"];
    auto& field = data_set.scalar_series["electric_field_z_v_per_m"];
    auto& te = data_set.scalar_series["electron_temperature_ev"];
    auto& ti = data_set.scalar_series["ion_temperature_ev"];

    density.reserve(cell_count);
    field.reserve(cell_count);
    te.reserve(cell_count);
    ti.reserve(cell_count);

    for (std::size_t i = 0; i < cell_count; ++i)
    {
        data_set.axis_values.push_back(static_cast<double>(i) * dz_m);

        if (i < cell_count / 3)
        {
            density.push_back((5.0e16 + density_bias) * (1.0 + 0.01 * std::sin(static_cast<double>(i))));
            field.push_back(45.0 + 4.0 * std::cos(static_cast<double>(i)));
            te.push_back(1.0);
            ti.push_back(1.0);
        }
        else if (i < 2 * cell_count / 3)
        {
            density.push_back((3.0e20 + density_bias) *
                              (1.0 + 0.015 * std::sin(static_cast<double>(i))));
            field.push_back(180.0 + 12.0 * std::cos(static_cast<double>(i)));
            te.push_back(2.2);
            ti.push_back(1.0);
        }
        else
        {
            density.push_back((3.0e21 + density_bias) *
                              (1.0 + 0.02 * std::sin(static_cast<double>(i))));
            field.push_back(760.0 + 24.0 * std::cos(static_cast<double>(i)));
            te.push_back(5.0);
            ti.push_back(1.0);
        }
    }

    return data_set;
}

SCDAT::Toolkit::PlasmaAnalysis::MultiScaleDecompositionSummary runBenchmarkSummary(
    const std::string& case_id,
    const ColumnarDataSet& profile,
    const SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters& plasma)
{
    MultiScaleSpatialDecomposer decomposer;
    const auto result = decomposer.decomposeFromProfile(profile, plasma);

    std::cout << "PLASMA_BENCHMARK_SUMMARY case_id=" << case_id
              << " fluid=" << result.summary.fluid_macro_cells
              << " hybrid=" << result.summary.hybrid_transition_cells
              << " pic=" << result.summary.local_pic_cells
              << " signature=" << result.summary.reproducibility_signature << std::endl;

    EXPECT_EQ(result.regions.size(), profile.axis_values.size());
    EXPECT_EQ(result.summary.fluid_macro_cells + result.summary.hybrid_transition_cells +
                  result.summary.local_pic_cells,
              profile.axis_values.size());
    EXPECT_GT(result.summary.fluid_macro_cells, 0u);
    EXPECT_GT(result.summary.hybrid_transition_cells, 0u);
    EXPECT_GT(result.summary.local_pic_cells, 0u);

    const auto rerun = decomposer.decomposeFromProfile(profile, plasma);
    EXPECT_EQ(result.summary.reproducibility_signature,
              rerun.summary.reproducibility_signature);
    EXPECT_EQ(result.summary.fluid_macro_cells, rerun.summary.fluid_macro_cells);
    EXPECT_EQ(result.summary.hybrid_transition_cells, rerun.summary.hybrid_transition_cells);
    EXPECT_EQ(result.summary.local_pic_cells, rerun.summary.local_pic_cells);

    return result.summary;
}

class SyntheticBoostProcess final : public IReactionCollisionProcess
{
    public:
        void evaluate(const SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters&,
                                    const SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment&,
                                    SCDAT::Toolkit::PlasmaAnalysis::ReactionCollisionState& state) const override
        {
                state.ionization_source_m3_per_s += 3.0e17;
                state.effective_collision_frequency_hz += 2.0e5;
                state.momentum_transfer_ratio += 0.2;
        }
};

} // namespace

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

TEST(PlasmaAnalysisSmokeTest, MultiScaleSpatialDecomposerProducesDeterministicSignature)
{
    PlasmaScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::PlasmaAnalysis::tryGetPlasmaScenarioPreset("hall_thruster_plume", preset));

    FluidAlgorithmAdapter adapter;
    ASSERT_TRUE(adapter.initialize(preset.config));
    for (std::size_t i = 0; i < 6; ++i)
    {
        ASSERT_TRUE(adapter.advance(preset.config.time_step_s));
    }

    const auto profile = adapter.buildProfileDataSet();
    MultiScaleSpatialDecomposer decomposer;
    const auto run_a = decomposer.decomposeFromProfile(profile, preset.config.initial_plasma);
    const auto run_b = decomposer.decomposeFromProfile(profile, preset.config.initial_plasma);

    ASSERT_EQ(run_a.regions.size(), profile.axis_values.size());
    ASSERT_EQ(run_a.regions, run_b.regions);
    EXPECT_EQ(run_a.summary.reproducibility_signature,
              run_b.summary.reproducibility_signature);

    const std::size_t total_cells =
        run_a.summary.fluid_macro_cells + run_a.summary.hybrid_transition_cells +
        run_a.summary.local_pic_cells;
    EXPECT_EQ(total_cells, run_a.regions.size());
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleSpatialDecomposerSeparatesThreeRegionTypes)
{
    const auto profile = makeSyntheticMultiscaleProfile();

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters plasma;
    plasma.neutral_density_m3 = 1.0e19;
    plasma.electron_temperature_ev = 1.2;
    plasma.ion_temperature_ev = 0.8;

    SCDAT::Toolkit::PlasmaAnalysis::MultiScaleSpatialThresholds thresholds;
    thresholds.smoothing_passes = 0;
    MultiScaleSpatialDecomposer decomposer(thresholds);
    const auto result = decomposer.decomposeFromProfile(profile, plasma);

    ASSERT_EQ(result.regions.size(), profile.axis_values.size());
    EXPECT_TRUE(containsRegion(result.regions, MultiScaleRegionType::FluidMacro));
    EXPECT_TRUE(containsRegion(result.regions, MultiScaleRegionType::HybridTransition));
    EXPECT_TRUE(containsRegion(result.regions, MultiScaleRegionType::LocalPIC));
    EXPECT_FALSE(result.summary.reproducibility_signature.empty());
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleSpatialDecomposerPartitionIsStableUnderSmallPerturbation)
{
    const auto baseline = makeSyntheticMultiscaleProfile();
    const auto perturbed = makePerturbedProfile(baseline);

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters plasma;
    plasma.neutral_density_m3 = 1.0e19;
    plasma.electron_temperature_ev = 1.2;
    plasma.ion_temperature_ev = 0.8;

    MultiScaleSpatialDecomposer decomposer;
    const auto baseline_result = decomposer.decomposeFromProfile(baseline, plasma);
    const auto perturbed_result = decomposer.decomposeFromProfile(perturbed, plasma);

    ASSERT_EQ(baseline_result.regions.size(), perturbed_result.regions.size());

    std::size_t changed = 0;
    for (std::size_t i = 0; i < baseline_result.regions.size(); ++i)
    {
        if (baseline_result.regions[i] != perturbed_result.regions[i])
        {
            ++changed;
        }
    }

    const double changed_ratio = static_cast<double>(changed) /
                                 static_cast<double>(baseline_result.regions.size());
    EXPECT_LE(changed_ratio, 0.2);
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleTimeSynchronizerConvergesWithCorrections)
{
    MultiScaleSynchronizationConfig config;
    config.macro_time_step_s = 2.0e-8;
    config.transition_substeps = 2;
    config.pic_substeps_per_transition = 3;
    config.thresholds.relative_charge_tolerance = 5.0e-3;
    config.thresholds.relative_energy_tolerance = 5.0e-3;
    config.thresholds.relative_field_tolerance = 5.0e-3;
    config.thresholds.min_iterations = 2;
    config.thresholds.max_iterations = 6;

    MultiScaleTimeSynchronizer synchronizer;
    ASSERT_TRUE(synchronizer.initialize(config));

    double charge_c = 10.0;
    double energy_j = 40.0;
    double field_norm = 15.0;
    const double base_charge_c = charge_c;
    const double base_energy_j = energy_j;
    const double base_field_norm = field_norm;

    auto advance_macro = [&](double dt) {
        charge_c += 30.0 * dt;
        energy_j += 60.0 * dt;
        field_norm += 15.0 * dt;
        return true;
    };
    auto advance_transition = [&](double dt) {
        charge_c += 12.0 * dt;
        energy_j += 24.0 * dt;
        field_norm += 8.0 * dt;
        return true;
    };
    auto advance_pic = [&](double dt) {
        charge_c += 8.0 * dt;
        energy_j += 16.0 * dt;
        field_norm += 5.0 * dt;
        return true;
    };
    auto sampler = [&]() {
        return SCDAT::Toolkit::PlasmaAnalysis::MultiScaleExchangeState{
            charge_c,
            energy_j,
            field_norm,
        };
    };

    auto correction = [&](const SCDAT::Toolkit::PlasmaAnalysis::MultiScaleSynchronizationResidual&) {
        charge_c = base_charge_c + 0.35 * (charge_c - base_charge_c);
        energy_j = base_energy_j + 0.35 * (energy_j - base_energy_j);
        field_norm = base_field_norm + 0.35 * (field_norm - base_field_norm);
    };

    ASSERT_TRUE(synchronizer.advanceOneMacroStep(advance_macro, advance_transition, advance_pic,
                                                 sampler, correction));

    const auto& status = synchronizer.getStatus();
    EXPECT_EQ(status.macro_steps_completed, 1u);
    EXPECT_GE(status.correction_iterations, 1u);
    EXPECT_EQ(status.convergence_reason, "criteria-met");
    EXPECT_TRUE(status.last_residual.converged);
    EXPECT_EQ(status.transition_substeps_executed % config.transition_substeps, 0u);
    EXPECT_EQ(status.pic_substeps_executed %
                  (config.transition_substeps * config.pic_substeps_per_transition),
              0u);
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleTimeSynchronizerReportsFailureOnNonConvergence)
{
    MultiScaleSynchronizationConfig config;
    config.macro_time_step_s = 1.0e-8;
    config.transition_substeps = 2;
    config.pic_substeps_per_transition = 2;
    config.thresholds.relative_charge_tolerance = 1.0e-8;
    config.thresholds.relative_energy_tolerance = 1.0e-8;
    config.thresholds.relative_field_tolerance = 1.0e-8;
    config.thresholds.min_iterations = 1;
    config.thresholds.max_iterations = 3;

    MultiScaleTimeSynchronizer synchronizer;
    ASSERT_TRUE(synchronizer.initialize(config));

    double charge_c = 1.0;
    double energy_j = 2.0;
    double field_norm = 3.0;

    auto advance_macro = [&](double dt) {
        charge_c += 3.0 * dt;
        energy_j += 6.0 * dt;
        field_norm += 9.0 * dt;
        return true;
    };
    auto advance_transition = [&](double dt) {
        charge_c += 1.5 * dt;
        energy_j += 3.0 * dt;
        field_norm += 4.5 * dt;
        return true;
    };
    auto advance_pic = [&](double dt) {
        charge_c += 1.0 * dt;
        energy_j += 2.0 * dt;
        field_norm += 3.0 * dt;
        return true;
    };
    auto sampler = [&]() {
        return SCDAT::Toolkit::PlasmaAnalysis::MultiScaleExchangeState{
            charge_c,
            energy_j,
            field_norm,
        };
    };

    EXPECT_FALSE(synchronizer.advanceOneMacroStep(advance_macro, advance_transition, advance_pic,
                                                  sampler));
    const auto& status = synchronizer.getStatus();
    EXPECT_EQ(status.convergence_reason, "max-iterations");
    EXPECT_FALSE(status.last_residual.converged);
}

TEST(PlasmaAnalysisSmokeTest, FluidPicHybridCouplingInterfaceClosesConservationGap)
{
    FluidPicCouplingConfig config;
    config.potential_jump_tolerance_v = 1.0e-9;
    config.charge_imbalance_tolerance_c = 1.0e-12;
    config.energy_imbalance_tolerance_j = 1.0e-12;
    config.correction_relaxation = 1.0;
    config.max_correction_iterations = 4;

    FluidPicHybridCouplingInterface coupling;
    ASSERT_TRUE(coupling.initialize(config));

    std::vector<FluidPicInterfaceState> interfaces(2);

    interfaces[0].interface_id = "if-0";
    interfaces[0].fluid_region_id = "fluid-macro-0";
    interfaces[0].pic_region_id = "local-pic-0";
    interfaces[0].transition_weight = 0.25;
    interfaces[0].fluid_potential_v = 120.0;
    interfaces[0].pic_potential_v = 95.0;
    interfaces[0].fluid_charge_current_a = 0.08;
    interfaces[0].pic_charge_current_a = 0.02;
    interfaces[0].fluid_energy_power_w = 15.0;
    interfaces[0].pic_energy_power_w = 9.0;
    interfaces[0].fluid_effective_conductivity_s_per_m = 1200.0;
    interfaces[0].pic_effective_conductivity_s_per_m = 2000.0;

    interfaces[1].interface_id = "if-1";
    interfaces[1].fluid_region_id = "fluid-macro-1";
    interfaces[1].pic_region_id = "local-pic-1";
    interfaces[1].transition_weight = 0.75;
    interfaces[1].fluid_potential_v = 85.0;
    interfaces[1].pic_potential_v = 112.0;
    interfaces[1].fluid_charge_current_a = 0.03;
    interfaces[1].pic_charge_current_a = 0.11;
    interfaces[1].fluid_energy_power_w = 7.0;
    interfaces[1].pic_energy_power_w = 18.0;
    interfaces[1].fluid_effective_conductivity_s_per_m = 900.0;
    interfaces[1].pic_effective_conductivity_s_per_m = 2500.0;

    const auto records = coupling.couple(2.5e-8, 2.5e-8, interfaces);
    ASSERT_EQ(records.size(), interfaces.size());

    const auto& metrics = coupling.getMetrics();
    EXPECT_GT(metrics.max_potential_jump_before_v, 0.0);
    EXPECT_LE(metrics.max_potential_jump_after_v, config.potential_jump_tolerance_v);
    EXPECT_TRUE(metrics.continuity_within_tolerance);
    EXPECT_TRUE(metrics.conservation_within_tolerance);
    EXPECT_LE(metrics.total_charge_imbalance_c, config.charge_imbalance_tolerance_c);
    EXPECT_LE(metrics.total_energy_imbalance_j, config.energy_imbalance_tolerance_j);
    EXPECT_GE(metrics.correction_iterations, 1u);

    for (const auto& state : interfaces)
    {
        EXPECT_NEAR(state.fluid_potential_v, state.pic_potential_v, 1.0e-12);
    }
    for (const auto& record : records)
    {
        EXPECT_FALSE(record.interface_id.empty());
        EXPECT_FALSE(record.fluid_region_id.empty());
        EXPECT_FALSE(record.pic_region_id.empty());
        EXPECT_EQ(record.units_tag, "SI[A,W,V,S/m,s,C,J]");
        EXPECT_TRUE(std::isfinite(record.boundary_potential_v));
    }
}

TEST(PlasmaAnalysisSmokeTest, FluidPicHybridCouplingInterfaceClampsTransitionWeight)
{
    FluidPicCouplingConfig config;
    config.potential_jump_tolerance_v = 1.0e-12;
    config.max_correction_iterations = 1;

    FluidPicHybridCouplingInterface coupling;
    ASSERT_TRUE(coupling.initialize(config));

    std::vector<FluidPicInterfaceState> interfaces(2);
    interfaces[0].interface_id = "if-fluid";
    interfaces[0].fluid_region_id = "fluid";
    interfaces[0].pic_region_id = "pic";
    interfaces[0].transition_weight = -0.5;
    interfaces[0].fluid_potential_v = 50.0;
    interfaces[0].pic_potential_v = 100.0;
    interfaces[0].fluid_effective_conductivity_s_per_m = 400.0;
    interfaces[0].pic_effective_conductivity_s_per_m = 800.0;

    interfaces[1].interface_id = "if-pic";
    interfaces[1].fluid_region_id = "fluid";
    interfaces[1].pic_region_id = "pic";
    interfaces[1].transition_weight = 2.5;
    interfaces[1].fluid_potential_v = 50.0;
    interfaces[1].pic_potential_v = 100.0;
    interfaces[1].fluid_effective_conductivity_s_per_m = 400.0;
    interfaces[1].pic_effective_conductivity_s_per_m = 800.0;

    const auto records = coupling.couple(1.0e-8, 1.0e-8, interfaces);
    ASSERT_EQ(records.size(), 2u);

    // weight < 0 clamps to fluid side.
    EXPECT_NEAR(records[0].boundary_potential_v, 50.0, 1.0e-12);
    EXPECT_NEAR(records[0].effective_conductivity_s_per_m, 400.0, 1.0e-12);

    // weight > 1 clamps to PIC side.
    EXPECT_NEAR(records[1].boundary_potential_v, 100.0, 1.0e-12);
    EXPECT_NEAR(records[1].effective_conductivity_s_per_m, 800.0, 1.0e-12);
}

TEST(PlasmaAnalysisSmokeTest, PlasmaReactionCollisionLibrarySupportsPluggableProcess)
{
    ReactionCollisionConfig config;
    config.enable_electron_impact_ionization = false;
    config.enable_radiative_recombination = false;
    config.enable_charge_exchange = false;
    config.enable_elastic_momentum_transfer = false;

    PlasmaReactionCollisionLibrary library;
    library.configure(config);

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters parameters;
    parameters.electron_density_m3 = 4.0e15;
    parameters.ion_density_m3 = 4.0e15;
    parameters.neutral_density_m3 = 8.0e18;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment assessment;
    assessment.collisionality = 0.4;

    const auto baseline = library.evaluate(parameters, assessment);
    EXPECT_EQ(baseline.active_processes, 0u);
    EXPECT_DOUBLE_EQ(baseline.ionization_source_m3_per_s, 0.0);

    library.registerProcess(std::make_shared<SyntheticBoostProcess>());
    const auto with_plugin = library.evaluate(parameters, assessment);

    EXPECT_EQ(with_plugin.active_processes, 1u);
    EXPECT_GT(with_plugin.ionization_source_m3_per_s, 0.0);
    EXPECT_GT(with_plugin.effective_collision_frequency_hz, 0.0);
    EXPECT_GT(with_plugin.momentum_transfer_ratio, 0.0);
}

TEST(PlasmaAnalysisSmokeTest, PlasmaReactionCollisionLibraryRespondsToThermalAndNeutralConditions)
{
    PlasmaReactionCollisionLibrary library;

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters low_state_parameters;
    low_state_parameters.electron_density_m3 = 3.0e15;
    low_state_parameters.ion_density_m3 = 3.0e15;
    low_state_parameters.electron_temperature_ev = 1.0;
    low_state_parameters.ion_temperature_ev = 0.2;
    low_state_parameters.neutral_density_m3 = 6.0e18;
    low_state_parameters.ion_mass_amu = 40.0;

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters high_state_parameters = low_state_parameters;
    high_state_parameters.electron_density_m3 = 1.4e16;
    high_state_parameters.ion_density_m3 = 1.4e16;
    high_state_parameters.electron_temperature_ev = 8.0;
    high_state_parameters.ion_temperature_ev = 1.0;
    high_state_parameters.neutral_density_m3 = 2.8e19;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment low_assessment;
    low_assessment.plasma_frequency_hz = 1.0e9;
    low_assessment.collisionality = 0.2;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment high_assessment = low_assessment;
    high_assessment.plasma_frequency_hz = 3.5e9;
    high_assessment.collisionality = 1.6;

    const auto low_state = library.evaluate(low_state_parameters, low_assessment);
    const auto high_state = library.evaluate(high_state_parameters, high_assessment);

    EXPECT_GT(high_state.ionization_source_m3_per_s, low_state.ionization_source_m3_per_s);
    EXPECT_GT(high_state.recombination_sink_m3_per_s, low_state.recombination_sink_m3_per_s);
    EXPECT_GT(high_state.effective_collision_frequency_hz, low_state.effective_collision_frequency_hz);
    EXPECT_GT(high_state.charge_exchange_frequency_hz, low_state.charge_exchange_frequency_hz);
    EXPECT_GE(high_state.active_processes, 4u);
}

TEST(PlasmaAnalysisSmokeTest, FluidAlgorithmAdapterExposesReactionCollisionDiagnostics)
{
    PlasmaScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::PlasmaAnalysis::tryGetPlasmaScenarioPreset("hall_thruster_plume", preset));

    FluidAlgorithmAdapter adapter;
    ASSERT_TRUE(adapter.initialize(preset.config));
    ASSERT_TRUE(adapter.advance(preset.config.time_step_s));

    const auto& status = adapter.getStatus();
    EXPECT_GT(status.reaction_active_processes, 0u);
    EXPECT_GE(status.ionization_source_m3_per_s, 0.0);
    EXPECT_GE(status.recombination_sink_m3_per_s, 0.0);
    EXPECT_GE(status.effective_collision_frequency_hz, 0.0);
    EXPECT_GE(status.charge_exchange_frequency_hz, 0.0);
    EXPECT_GE(status.reaction_momentum_transfer_ratio, 0.0);
    EXPECT_LE(status.reaction_momentum_transfer_ratio, 1.0);

    const auto data_set = adapter.buildProfileDataSet();
    EXPECT_TRUE(data_set.metadata.find("ionization_source_m3_per_s") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("recombination_sink_m3_per_s") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("effective_collision_frequency_hz") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("charge_exchange_frequency_hz") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("reaction_energy_loss_ev_per_s") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("reaction_momentum_transfer_ratio") != data_set.metadata.end());

    EXPECT_TRUE(std::isfinite(std::stod(data_set.metadata.at("effective_collision_frequency_hz"))));
}

TEST(PlasmaAnalysisSmokeTest, PlasmaAdvancedClosureModelRelaxesNonEquilibriumTemperature)
{
    AdvancedClosureConfig config;
    config.enable_non_equilibrium_closure = true;
    config.enable_turbulence_closure = false;
    config.non_equilibrium_model = NonEquilibriumClosureModel::TwoTemperatureRelaxation;
    config.non_equilibrium_relaxation_gain = 0.45;

    PlasmaAdvancedClosureModel closure_model;
    closure_model.configure(config);

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters parameters;
    parameters.electron_density_m3 = 6.0e15;
    parameters.ion_density_m3 = 6.0e15;
    parameters.electron_temperature_ev = 8.0;
    parameters.ion_temperature_ev = 0.7;
    parameters.electric_field_v_per_m = 600.0;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment assessment;
    assessment.collisionality = 1.5;

    SCDAT::Toolkit::PlasmaAnalysis::ReactionCollisionState reaction_state;
    reaction_state.effective_collision_frequency_hz = 3.0e6;

    const auto state =
        closure_model.evaluate(parameters, assessment, reaction_state, 2.0e-5, 2.0e-8);

    EXPECT_GT(state.non_equilibrium_ratio, 1.0);
    EXPECT_GT(state.non_equilibrium_relaxation_rate_hz, 0.0);
    EXPECT_LT(state.electron_temperature_delta_ev, 0.0);
    EXPECT_EQ(state.active_terms, 1u);
}

TEST(PlasmaAnalysisSmokeTest, PlasmaAdvancedClosureModelBuildsTurbulenceResponse)
{
    AdvancedClosureConfig config;
    config.enable_non_equilibrium_closure = false;
    config.enable_turbulence_closure = true;
    config.turbulence_model = TurbulenceClosureModel::CollisionalDamping;
    config.turbulence_gain = 0.9;

    PlasmaAdvancedClosureModel closure_model;
    closure_model.configure(config);

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters parameters;
    parameters.electron_density_m3 = 2.0e17;
    parameters.electron_temperature_ev = 3.5;
    parameters.ion_temperature_ev = 0.8;
    parameters.electric_field_v_per_m = 2400.0;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment assessment;
    assessment.collisionality = 2.0;

    SCDAT::Toolkit::PlasmaAnalysis::ReactionCollisionState reaction_state;
    const auto state =
        closure_model.evaluate(parameters, assessment, reaction_state, 5.0e-5, 2.5e-8);

    EXPECT_GT(state.turbulence_intensity, 0.0);
    EXPECT_GT(state.turbulence_eddy_diffusivity_m2_per_s, 0.0);
    EXPECT_GT(state.turbulence_dissipation_rate_w_per_m3, 0.0);
    EXPECT_LT(state.density_correction_m3, 0.0);
    EXPECT_EQ(state.active_terms, 1u);
}

TEST(PlasmaAnalysisSmokeTest, FluidAlgorithmAdapterExposesAdvancedClosureDiagnostics)
{
    PlasmaScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::PlasmaAnalysis::tryGetPlasmaScenarioPreset("hall_thruster_plume", preset));

    preset.config.advanced_closure.enable_non_equilibrium_closure = true;
    preset.config.advanced_closure.enable_turbulence_closure = true;
    preset.config.advanced_closure.non_equilibrium_model =
        NonEquilibriumClosureModel::CollisionalThermalization;
    preset.config.advanced_closure.turbulence_model =
        TurbulenceClosureModel::MixingLengthEddyDiffusivity;

    FluidAlgorithmAdapter adapter;
    ASSERT_TRUE(adapter.initialize(preset.config));
    ASSERT_TRUE(adapter.advance(preset.config.time_step_s));

    const auto& status = adapter.getStatus();
    EXPECT_TRUE(status.advanced_closure_enabled);
    EXPECT_GT(status.non_equilibrium_ratio, 1.0);
    EXPECT_GE(status.non_equilibrium_relaxation_rate_hz, 0.0);
    EXPECT_GE(status.turbulence_intensity, 0.0);
    EXPECT_GE(status.turbulence_eddy_diffusivity_m2_per_s, 0.0);
    EXPECT_GE(status.turbulence_dissipation_rate_w_per_m3, 0.0);
    EXPECT_NE(status.closure_active_terms, 0u);

    const auto data_set = adapter.buildProfileDataSet();
    EXPECT_TRUE(data_set.metadata.find("advanced_closure_enabled") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("non_equilibrium_ratio") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("non_equilibrium_relaxation_rate_hz") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("turbulence_eddy_diffusivity_m2_per_s") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("turbulence_dissipation_rate_w_per_m3") != data_set.metadata.end());
    EXPECT_TRUE(data_set.metadata.find("closure_active_terms") != data_set.metadata.end());

    EXPECT_TRUE(std::isfinite(std::stod(data_set.metadata.at("non_equilibrium_ratio"))));
    EXPECT_TRUE(
        std::isfinite(std::stod(data_set.metadata.at("turbulence_eddy_diffusivity_m2_per_s"))));
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleRegressionBenchmark1DLayeredBaseline)
{
    auto profile = makeLayeredBenchmarkProfile(36, 1.0e-6, 0.0);

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters plasma;
    plasma.neutral_density_m3 = 1.0e22;
    plasma.electron_temperature_ev = 1.2;
    plasma.ion_temperature_ev = 1.0;

    const auto summary = runBenchmarkSummary("plasma_benchmark_1d", profile, plasma);
    EXPECT_GE(summary.fluid_macro_cells, 8u);
    EXPECT_GE(summary.hybrid_transition_cells, 8u);
    EXPECT_GE(summary.local_pic_cells, 8u);
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleRegressionBenchmark2DProjectedBaseline)
{
    auto profile = makeLayeredBenchmarkProfile(42, 8.0e-7, 2.0e19);

    for (std::size_t i = 0; i < profile.axis_values.size(); ++i)
    {
        profile.scalar_series["electron_density_m3"][i] *=
            1.0 + 0.03 * std::sin(0.3 * static_cast<double>(i));
        profile.scalar_series["electric_field_z_v_per_m"][i] *=
            1.0 + 0.025 * std::cos(0.25 * static_cast<double>(i));
    }

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters plasma;
    plasma.neutral_density_m3 = 8.0e21;
    plasma.electron_temperature_ev = 1.5;
    plasma.ion_temperature_ev = 0.9;

    const auto summary = runBenchmarkSummary("plasma_benchmark_2d_projected", profile, plasma);
    EXPECT_GE(summary.fluid_macro_cells, 9u);
    EXPECT_GE(summary.hybrid_transition_cells, 9u);
    EXPECT_GE(summary.local_pic_cells, 9u);
}

TEST(PlasmaAnalysisSmokeTest, MultiScaleRegressionBenchmark3DProjectedBaseline)
{
    auto profile = makeLayeredBenchmarkProfile(48, 1.0e-6, 0.0);

    for (std::size_t i = 0; i < profile.axis_values.size(); ++i)
    {
        const double k = static_cast<double>(i);
        profile.scalar_series["electron_density_m3"][i] *=
            1.0 + 0.04 * std::sin(0.17 * k) * std::cos(0.11 * k);
        profile.scalar_series["electric_field_z_v_per_m"][i] *=
            1.0 + 0.03 * std::sin(0.21 * k);
        profile.scalar_series["electron_temperature_ev"][i] *=
            1.0 + 0.12 * std::cos(0.19 * k);

        if (i >= (2 * profile.axis_values.size()) / 3)
        {
            profile.scalar_series["electron_temperature_ev"][i] *= 1.15;
        }
        else if (i >= profile.axis_values.size() / 3)
        {
            profile.scalar_series["electron_temperature_ev"][i] *= 1.05;
        }
    }

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters plasma;
    plasma.neutral_density_m3 = 1.0e22;
    plasma.electron_temperature_ev = 1.8;
    plasma.ion_temperature_ev = 0.85;

    const auto summary = runBenchmarkSummary("plasma_benchmark_3d_projected", profile, plasma);
    EXPECT_GE(summary.fluid_macro_cells, 10u);
    EXPECT_GE(summary.hybrid_transition_cells, 10u);
    EXPECT_GE(summary.local_pic_cells, 10u);
}

TEST(PlasmaAnalysisSmokeTest, DensePlasmaBoundaryLayerBuildsSheathPresheathScales)
{
    DensePlasmaBoundaryLayer boundary_layer;

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters parameters;
    parameters.electron_temperature_ev = 3.0;
    parameters.ion_temperature_ev = 0.8;
    parameters.ion_mass_amu = 16.0;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment assessment;
    assessment.is_dense = true;
    assessment.debye_length_m = 2.0e-6;
    assessment.collisionality = 0.8;

    const double surface_potential_v = -48.0;
    const auto state = boundary_layer.analyze(parameters, assessment, surface_potential_v);

    EXPECT_GT(state.sheath_thickness_m, assessment.debye_length_m);
    EXPECT_GT(state.presheath_thickness_m, state.sheath_thickness_m);
    EXPECT_LT(state.debye_to_sheath_ratio, 1.0);
    EXPECT_LT(state.sheath_potential_drop_v, 0.0);
    EXPECT_LT(state.presheath_potential_drop_v, 0.0);
    EXPECT_NEAR(state.sheath_potential_drop_v + state.presheath_potential_drop_v,
                surface_potential_v,
                1.0e-9);
    EXPECT_GE(state.ion_mach_at_sheath_edge, 1.0);
    EXPECT_NEAR(state.wall_field_v_per_m * state.sheath_thickness_m,
                state.sheath_potential_drop_v,
                1.0e-9);
}

TEST(PlasmaAnalysisSmokeTest, DensePlasmaBoundaryLayerRespondsToCollisionalityAndIonMass)
{
    DensePlasmaBoundaryLayer boundary_layer;

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters light_ion_parameters;
    light_ion_parameters.electron_temperature_ev = 2.5;
    light_ion_parameters.ion_mass_amu = 1.0;

    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters heavy_ion_parameters = light_ion_parameters;
    heavy_ion_parameters.ion_mass_amu = 131.0;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment low_collisionality;
    low_collisionality.is_dense = false;
    low_collisionality.debye_length_m = 4.0e-6;
    low_collisionality.collisionality = 0.15;

    SCDAT::Toolkit::PlasmaAnalysis::DensePlasmaAssessment high_collisionality = low_collisionality;
    high_collisionality.collisionality = 2.2;

    const auto light_ion_state =
        boundary_layer.analyze(light_ion_parameters, low_collisionality, -30.0);
    const auto heavy_ion_state =
        boundary_layer.analyze(heavy_ion_parameters, low_collisionality, -30.0);
    const auto low_collisionality_state =
        boundary_layer.analyze(heavy_ion_parameters, low_collisionality, -30.0);
    const auto high_collisionality_state =
        boundary_layer.analyze(heavy_ion_parameters, high_collisionality, -30.0);

    EXPECT_GT(light_ion_state.bohm_velocity_m_per_s, heavy_ion_state.bohm_velocity_m_per_s);
    EXPECT_GT(high_collisionality_state.sheath_thickness_m,
              low_collisionality_state.sheath_thickness_m);
    EXPECT_GT(high_collisionality_state.presheath_thickness_m,
              low_collisionality_state.presheath_thickness_m);
    EXPECT_LT(std::abs(high_collisionality_state.wall_field_v_per_m),
              std::abs(low_collisionality_state.wall_field_v_per_m));
}
