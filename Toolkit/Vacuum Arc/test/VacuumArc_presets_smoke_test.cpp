#include "ArcPICEmissionModel.h"
#include "ArcFieldEmissionBoundaryCondition.h"
#include "SurfaceDischargeArcAlgorithm.h"
#include "VacuumArcCases.h"
#include "DataAnalyzer.h"

#include "../../../Tools/Basic/include/Constants.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using SCDAT::Toolkit::VacuumArc::ArcEmissionStrategyType;
using SCDAT::Toolkit::VacuumArc::ArcFieldEmissionBoundaryCondition;
using SCDAT::Toolkit::VacuumArc::ArcIntegratorMode;
using SCDAT::Toolkit::VacuumArc::ArcPicAlignmentMode;
using SCDAT::Toolkit::VacuumArc::ArcPICEmissionModel;
using SCDAT::Toolkit::VacuumArc::SurfaceCircuitLoadModel;
using SCDAT::Toolkit::VacuumArc::SurfaceDischargeArcAlgorithm;
using SCDAT::Toolkit::VacuumArc::VacuumArcScenarioPreset;

namespace
{
double relativeTolerance(double expected, double fraction)
{
    return std::max(1.0e-12, std::abs(expected) * fraction);
}

double lastValue(const SCDAT::Output::ColumnarDataSet& data_set, const std::string& key)
{
    const auto it = data_set.scalar_series.find(key);
    if (it == data_set.scalar_series.end() || it->second.empty())
    {
        return 0.0;
    }
    return it->second.back();
}

double maxAbsValue(const SCDAT::Output::ColumnarDataSet& data_set, const std::string& key)
{
    const auto it = data_set.scalar_series.find(key);
    if (it == data_set.scalar_series.end() || it->second.empty())
    {
        return 0.0;
    }

    double max_abs = 0.0;
    for (const double value : it->second)
    {
        max_abs = std::max(max_abs, std::abs(value));
    }
    return max_abs;
}

double sumValue(const SCDAT::Output::ColumnarDataSet& data_set, const std::string& key)
{
    const auto it = data_set.scalar_series.find(key);
    if (it == data_set.scalar_series.end() || it->second.empty())
    {
        return 0.0;
    }

    double sum = 0.0;
    for (const double value : it->second)
    {
        sum += value;
    }
    return sum;
}

bool hasColumn(const SCDAT::Output::ColumnarDataSet& data_set, const std::string& key)
{
    const auto it = data_set.scalar_series.find(key);
    return it != data_set.scalar_series.end() && !it->second.empty();
}

std::string readTextFile(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

double percentileAbsValue(const SCDAT::Output::ColumnarDataSet& data_set, const std::string& key,
                         double quantile)
{
    const auto it = data_set.scalar_series.find(key);
    if (it == data_set.scalar_series.end() || it->second.empty())
    {
        return 0.0;
    }

    std::vector<double> values;
    values.reserve(it->second.size());
    for (double value : it->second)
    {
        values.push_back(std::abs(value));
    }

    const double clamped_quantile = std::clamp(quantile, 0.0, 1.0);
    const std::size_t index = static_cast<std::size_t>(
        std::floor(clamped_quantile * static_cast<double>(values.size() - 1)));
    std::nth_element(values.begin(), values.begin() + index, values.end());
    return values[index];
}

std::size_t countValuesAbove(const SCDAT::Output::ColumnarDataSet& data_set,
                             const std::string& key, double threshold)
{
    const auto it = data_set.scalar_series.find(key);
    if (it == data_set.scalar_series.end() || it->second.empty())
    {
        return 0;
    }

    std::size_t count = 0;
    for (double value : it->second)
    {
        if (value > threshold)
        {
            ++count;
        }
    }
    return count;
}
} // namespace

TEST(VacuumArcPresetSmokeTest, AllPresetsAdvanceAndExport)
{
    const auto preset_names = SCDAT::Toolkit::VacuumArc::listVacuumArcScenarioPresetNames();
    ASSERT_FALSE(preset_names.empty());

    for (const auto& preset_name : preset_names)
    {
        VacuumArcScenarioPreset preset;
        ASSERT_TRUE(SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset(preset_name, preset));

        SurfaceDischargeArcAlgorithm algorithm;
        ASSERT_TRUE(algorithm.initialize(preset.config));

        for (std::size_t step = 0; step < preset.steps; ++step)
        {
            ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        }

        const auto csv_path =
            std::filesystem::temp_directory_path() / (std::string("vacuum_arc_") + preset.name + ".csv");
        ASSERT_TRUE(algorithm.exportResults(csv_path));
        EXPECT_TRUE(std::filesystem::exists(csv_path));
    }
}

TEST(VacuumArcPresetSmokeTest, EmissionStrategySwitchProducesFiniteDensity)
{
    ArcPICEmissionModel model;
    model.setStrategyType(ArcEmissionStrategyType::LinearThreshold);
    const double linear_value = model.computeEmissionCurrentDensity(5.5e7, 700.0);

    model.setStrategyType(ArcEmissionStrategyType::FowlerNordheimLike);
    const double fn_value = model.computeEmissionCurrentDensity(5.5e7, 700.0);

    EXPECT_GE(linear_value, 0.0);
    EXPECT_GE(fn_value, 0.0);
    EXPECT_NE(linear_value, fn_value);
}

TEST(VacuumArcPresetSmokeTest, EmissionStrategyExtendedModulationChangesResponse)
{
    ArcPICEmissionModel baseline_model;
    baseline_model.setStrategyType(ArcEmissionStrategyType::FowlerNordheimLike);

    SCDAT::Toolkit::VacuumArc::ArcEmissionStrategyParameters baseline_parameters;
    baseline_model.configureStrategy(baseline_parameters);
    const double baseline = baseline_model.computeEmissionCurrentDensity(4.5e7, 650.0);

    ArcPICEmissionModel tuned_model;
    tuned_model.setStrategyType(ArcEmissionStrategyType::FowlerNordheimLike);

    SCDAT::Toolkit::VacuumArc::ArcEmissionStrategyParameters tuned_parameters;
    tuned_parameters.field_enhancement_factor = 1.35;
    tuned_parameters.regional_gain = 1.2;
    tuned_parameters.thermal_activation_exponent = 1.1;
    tuned_model.configureStrategy(tuned_parameters);

    const double tuned = tuned_model.computeEmissionCurrentDensity(4.5e7, 650.0);
    EXPECT_TRUE(std::isfinite(baseline));
    EXPECT_TRUE(std::isfinite(tuned));
    EXPECT_GE(tuned, baseline);
}

TEST(VacuumArcPresetSmokeTest, FieldEmissionBoundaryAccumulatesSubParticleBudget)
{
    ArcFieldEmissionBoundaryCondition::Parameters parameters;
    parameters.emission_area_m2 = 1.0e-8;
    parameters.emitted_electron_speed_m_per_s = 1.0e6;
    parameters.max_emitted_particles_per_event = 16;

    ArcFieldEmissionBoundaryCondition boundary(parameters);

    const double dt = 1.0e-9;
    const double expected_emitted_per_step = 0.25;
    const double required_current_density =
        expected_emitted_per_step * SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge /
        (parameters.emission_area_m2 * dt);
    boundary.setEmissionCurrentDensity(required_current_density);

    std::size_t total_emitted = 0;
    for (std::size_t step = 0; step < 3; ++step)
    {
        SCDAT::Particle::ParticleTypeDef incident_ion(
            1000 + step, SCDAT::Particle::ParticleType::POSITIVE_ION,
            SCDAT::Geometry::Point3D(0.0, 0.0, 0.0), SCDAT::Geometry::Vector3D(0.0, 0.0, -1.0e5),
            SCDAT::Basic::Constants::PhysicsConstants::ProtonMass,
            SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge, 1.0);

        const auto emitted = boundary.processParticle(
            incident_ion, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
            SCDAT::Geometry::Vector3D(0.0, 0.0, 1.0), dt);
        EXPECT_FALSE(incident_ion.isActive());
        EXPECT_TRUE(emitted.empty());
        total_emitted += emitted.size();
    }

    SCDAT::Particle::ParticleTypeDef incident_ion(
        2000, SCDAT::Particle::ParticleType::POSITIVE_ION,
        SCDAT::Geometry::Point3D(0.0, 0.0, 0.0), SCDAT::Geometry::Vector3D(0.0, 0.0, -1.0e5),
        SCDAT::Basic::Constants::PhysicsConstants::ProtonMass,
        SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge, 1.0);
    const auto emitted = boundary.processParticle(
        incident_ion, SCDAT::Geometry::Point3D(0.0, 0.0, 0.0),
        SCDAT::Geometry::Vector3D(0.0, 0.0, 1.0), dt);
    total_emitted += emitted.size();

    ASSERT_EQ(emitted.size(), 1u);
    EXPECT_EQ(total_emitted, 1u);
    EXPECT_EQ(emitted.front().getType(), SCDAT::Particle::ParticleType::FIELD_EMISSION_ELECTRON);
    EXPECT_LT(emitted.front().getCharge(), 0.0);
}

TEST(VacuumArcPresetSmokeTest, PiccoreBoundaryCouplingProducesFiniteFeedback)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("triple_junction_flashover", preset));

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto& status = algorithm.getStatus();
    EXPECT_TRUE(std::isfinite(status.pic_boundary_net_current_a));
    EXPECT_TRUE(std::isfinite(status.pic_surface_charge_delta_c));
}

TEST(VacuumArcPresetSmokeTest, PiccoreNetCurrentDrivesSurfaceCircuitPotential)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 8);
    preset.config.surface_capacitance_f = 5.0e-13;
    preset.config.surface_leakage_conductance_s = 0.0;
    preset.config.max_surface_potential_step_v = 5.0;

    const double initial_potential_v = preset.config.surface_potential_v;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    bool potential_changed = false;
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        const auto& status = algorithm.getStatus();
        if (std::abs(status.surface_potential_v - initial_potential_v) > 1.0e-12)
        {
            potential_changed = true;
        }
    }

    const auto& status = algorithm.getStatus();
    EXPECT_TRUE(std::isfinite(status.surface_potential_v));
    EXPECT_TRUE(std::isfinite(status.surface_circuit_charge_c));
    EXPECT_TRUE(std::isfinite(status.surface_circuit_drive_current_a));
    EXPECT_TRUE(std::isfinite(status.surface_circuit_leak_current_a));
    EXPECT_NEAR(status.surface_circuit_drive_current_a, status.pic_boundary_net_current_a, 1.0e-18);
    EXPECT_TRUE(potential_changed);
}

TEST(VacuumArcPresetSmokeTest, LegacyAlignmentModeKeepsPresetRunnable)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("triple_junction_flashover", preset));

    preset.config.alignment_mode = ArcPicAlignmentMode::LegacyBaseline;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto& status = algorithm.getStatus();
    EXPECT_TRUE(std::isfinite(status.total_discharge_current_a));
    EXPECT_TRUE(std::isfinite(status.surface_potential_v));
    EXPECT_DOUBLE_EQ(status.pic_boundary_net_current_a, 0.0);
}

TEST(VacuumArcPresetSmokeTest, FrozenBaselineMetricsStayWithinTolerance)
{
    struct BaselineExpectation
    {
        const char* preset_name;
        double peak_discharge_current_a;
        double peak_current_density_a_per_m2;
        double peak_cathode_temperature_k;
        double final_discharge_current_a;
        double final_surface_potential_v;
        double final_surface_charge_density_c_per_m2;
        double max_abs_pic_boundary_net_current_a;
    };

    const std::vector<BaselineExpectation> expectations{
        {"microgap_flashover", 4.20220278e-08, 2.99976437e+01, 5.20000046e+02, 3.81767922e-08,
         1.49999850e+02, 2.99963470e-05, 4.32587691e-09},
        {"triple_junction_flashover", 8.19962561e-08, 1.55752087e+02, 6.50000124e+02,
         7.17423125e-08, 2.19999923e+02, 7.99942322e-05, 1.15356718e-08},
        {"restrike_recovery", 4.66585220e-08, 2.22942228e+01, 5.60000031e+02, 4.40950322e-08,
         7.99999400e+01, 2.49980133e-05, 2.88391794e-09},
    };

    SCDAT::Output::DataAnalyzer analyzer;
    for (const auto& expected : expectations)
    {
        VacuumArcScenarioPreset preset;
        ASSERT_TRUE(SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset(expected.preset_name,
                                                                              preset));

        SurfaceDischargeArcAlgorithm algorithm;
        ASSERT_TRUE(algorithm.initialize(preset.config));
        for (std::size_t step = 0; step < preset.steps; ++step)
        {
            ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        }

        const auto csv_path = std::filesystem::temp_directory_path() /
                              (std::string("vacuum_arc_baseline_anchor_") + expected.preset_name +
                               ".csv");
        ASSERT_TRUE(algorithm.exportResults(csv_path));

        const auto summary = analyzer.summarizeCsv(csv_path);
        const auto data_set = analyzer.loadCsv(csv_path);

        ASSERT_TRUE(data_set.isConsistent());
        ASSERT_TRUE(data_set.hasOnlyFiniteValues());

        const auto peak_current_it = summary.column_summaries.find("discharge_current_a");
        ASSERT_TRUE(peak_current_it != summary.column_summaries.end());
        EXPECT_NEAR(peak_current_it->second.maximum, expected.peak_discharge_current_a,
                    relativeTolerance(expected.peak_discharge_current_a, 0.05));

        const auto peak_density_it = summary.column_summaries.find("current_density_a_per_m2");
        ASSERT_TRUE(peak_density_it != summary.column_summaries.end());
        EXPECT_NEAR(peak_density_it->second.maximum, expected.peak_current_density_a_per_m2,
                    relativeTolerance(expected.peak_current_density_a_per_m2, 0.05));

        const auto peak_temp_it = summary.column_summaries.find("cathode_temperature_k");
        ASSERT_TRUE(peak_temp_it != summary.column_summaries.end());
        EXPECT_NEAR(peak_temp_it->second.maximum, expected.peak_cathode_temperature_k,
                    relativeTolerance(expected.peak_cathode_temperature_k, 0.01));

        EXPECT_NEAR(lastValue(data_set, "discharge_current_a"), expected.final_discharge_current_a,
                    relativeTolerance(expected.final_discharge_current_a, 0.05));
        EXPECT_NEAR(lastValue(data_set, "surface_potential_v"), expected.final_surface_potential_v,
                    relativeTolerance(expected.final_surface_potential_v, 0.01));
        EXPECT_NEAR(lastValue(data_set, "surface_charge_density_c_per_m2"),
                    expected.final_surface_charge_density_c_per_m2,
                    relativeTolerance(expected.final_surface_charge_density_c_per_m2, 0.05));
        EXPECT_NEAR(maxAbsValue(data_set, "pic_boundary_net_current_a"),
                    expected.max_abs_pic_boundary_net_current_a,
                    relativeTolerance(expected.max_abs_pic_boundary_net_current_a, 0.10));
    }
}

TEST(VacuumArcPresetSmokeTest, ConfigurableIntegratorAndChannelConstraintsTakeEffect)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("triple_junction_flashover", preset));

    preset.config.channel_parameters.max_current_density_a_per_m2 = 40.0;
    preset.config.channel_parameters.max_conductivity_s_per_m = 2.0e3;
    preset.config.integrator_parameters.mode = ArcIntegratorMode::BoundedRelaxation;
    preset.config.integrator_parameters.max_current_density_a_per_m2 = 35.0;
    preset.config.integrator_parameters.max_conductivity_s_per_m = 1.8e3;
    preset.config.integrator_parameters.current_density_growth_per_ns = 0.15;
    preset.config.integrator_parameters.conductivity_growth_per_ns = 0.12;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto& status = algorithm.getStatus();
    EXPECT_TRUE(std::isfinite(status.peak_current_density_a_per_m2));
    EXPECT_LE(status.peak_current_density_a_per_m2,
              preset.config.integrator_parameters.max_current_density_a_per_m2 + 1.0e-6);
    EXPECT_TRUE(std::isfinite(status.channel_conductivity_s_per_m));
    EXPECT_LE(status.channel_conductivity_s_per_m,
              preset.config.integrator_parameters.max_conductivity_s_per_m + 1.0e-6);
}

TEST(VacuumArcPresetSmokeTest, CollisionStageProducesCollisionMetrics)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 24);
    preset.config.enable_pic_mcc_collisions = true;
    preset.config.enable_collision_reaction_breakdown = true;
    preset.config.collision_reconfigure_interval_steps = 32;
    preset.config.collision_neutral_density_floor_m3 = 8.0e20;
    preset.config.collision_ion_density_floor_m3 = 1.0e16;
    preset.config.collision_electron_density_floor_m3 = 1.0e16;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    double total_collision_events = 0.0;
    double total_created_particles = 0.0;
    double total_ionization_events = 0.0;
    double total_excitation_events = 0.0;
    double total_charge_exchange_events = 0.0;
    double max_collision_event_rate = 0.0;
    bool fallback_triggered = false;
    double max_effective_neutral_density = 0.0;
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        const auto& status = algorithm.getStatus();
        total_collision_events += status.collision_events_step;
        total_created_particles += status.collision_particles_created_step;
        total_ionization_events += status.collision_ionization_events_step;
        total_excitation_events += status.collision_excitation_events_step;
        total_charge_exchange_events += status.collision_charge_exchange_events_step;
        max_collision_event_rate =
            std::max(max_collision_event_rate, status.collision_event_rate_per_s);
        fallback_triggered = fallback_triggered || status.collision_stage_fallback_triggered;
        max_effective_neutral_density =
            std::max(max_effective_neutral_density, status.collision_effective_neutral_density_m3);
    }

    EXPECT_GT(total_collision_events, 0.0);
    EXPECT_GE(total_created_particles, 0.0);
    EXPECT_GE(total_ionization_events, 0.0);
    EXPECT_GE(total_excitation_events, 0.0);
    EXPECT_GE(total_charge_exchange_events, 0.0);
    EXPECT_GE(total_collision_events,
              total_ionization_events + total_excitation_events + total_charge_exchange_events);
    EXPECT_TRUE(std::isfinite(max_collision_event_rate));
    EXPECT_FALSE(fallback_triggered);
    EXPECT_TRUE(std::isfinite(max_effective_neutral_density));
    EXPECT_GT(max_effective_neutral_density, 0.0);

    const auto csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_collision_metrics_smoke.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));

    SCDAT::Output::DataAnalyzer analyzer;
    const auto data_set = analyzer.loadCsv(csv_path);
    ASSERT_TRUE(data_set.isConsistent());
    ASSERT_TRUE(data_set.hasOnlyFiniteValues());

    const auto collision_events_it = data_set.scalar_series.find("collision_events_step");
    ASSERT_TRUE(collision_events_it != data_set.scalar_series.end());
    EXPECT_FALSE(collision_events_it->second.empty());
    EXPECT_TRUE(hasColumn(data_set, "collision_ionization_events_step"));
    EXPECT_TRUE(hasColumn(data_set, "collision_excitation_events_step"));
    EXPECT_TRUE(hasColumn(data_set, "collision_charge_exchange_events_step"));
    EXPECT_TRUE(hasColumn(data_set, "collision_event_rate_per_s"));
    EXPECT_TRUE(hasColumn(data_set, "collision_stage_fallback_triggered"));
}

TEST(VacuumArcPresetSmokeTest, CollisionStageStructuredReactionFeedbackIsFiniteAndExported)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 20);
    preset.config.enable_pic_mcc_collisions = true;
    preset.config.enable_collision_reaction_breakdown = true;
    preset.config.collision_emission_feedback_gain = 0.0;
    preset.config.collision_channel_feedback_gain = 0.0;
    preset.config.collision_ionization_emission_feedback_gain = 0.30;
    preset.config.collision_excitation_emission_feedback_gain = 0.20;
    preset.config.collision_charge_exchange_emission_feedback_gain = 0.10;
    preset.config.collision_ionization_channel_feedback_gain = 0.25;
    preset.config.collision_excitation_channel_feedback_gain = 0.15;
    preset.config.collision_charge_exchange_channel_feedback_gain = 0.05;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    bool observed_emission_feedback = false;
    bool observed_channel_feedback = false;
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        const auto& status = algorithm.getStatus();

        EXPECT_TRUE(std::isfinite(status.collision_ionization_fraction_step));
        EXPECT_TRUE(std::isfinite(status.collision_excitation_fraction_step));
        EXPECT_TRUE(std::isfinite(status.collision_charge_exchange_fraction_step));
        EXPECT_TRUE(std::isfinite(status.collision_reaction_weighted_emission_feedback_step));
        EXPECT_TRUE(std::isfinite(status.collision_reaction_weighted_channel_feedback_step));

        EXPECT_GE(status.collision_ionization_fraction_step, 0.0);
        EXPECT_GE(status.collision_excitation_fraction_step, 0.0);
        EXPECT_GE(status.collision_charge_exchange_fraction_step, 0.0);
        const double reaction_fraction_sum = status.collision_ionization_fraction_step +
                                             status.collision_excitation_fraction_step +
                                             status.collision_charge_exchange_fraction_step;
        EXPECT_LE(reaction_fraction_sum, 1.0 + 1.0e-9);

        observed_emission_feedback = observed_emission_feedback ||
                                     (status.collision_reaction_weighted_emission_feedback_step >
                                      1.0 + 1.0e-9);
        observed_channel_feedback = observed_channel_feedback ||
                                    (status.collision_reaction_weighted_channel_feedback_step >
                                     1.0 + 1.0e-9);
    }

    EXPECT_TRUE(observed_emission_feedback);
    EXPECT_TRUE(observed_channel_feedback);

    const auto csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_collision_feedback_stage_c.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));

    SCDAT::Output::DataAnalyzer analyzer;
    const auto data_set = analyzer.loadCsv(csv_path);
    ASSERT_TRUE(data_set.isConsistent());
    ASSERT_TRUE(data_set.hasOnlyFiniteValues());

    EXPECT_TRUE(hasColumn(data_set, "collision_ionization_fraction_step"));
    EXPECT_TRUE(hasColumn(data_set, "collision_excitation_fraction_step"));
      EXPECT_TRUE(hasColumn(data_set, "collision_charge_exchange_fraction_step"));
      EXPECT_TRUE(hasColumn(data_set, "collision_reaction_weighted_emission_feedback_step"));
      EXPECT_TRUE(hasColumn(data_set, "collision_reaction_weighted_channel_feedback_step"));

      const auto sidecar = readTextFile(csv_path.string() + ".metadata.json");
      EXPECT_NE(
          sidecar.find("\"collision_diagnostic_contract_id\": \"vacuum-arc-collision-emission-channel-v1\""),
          std::string::npos);
      EXPECT_NE(
          sidecar.find("\"benchmark_metrics_contract_id\": \"vacuum-arc-benchmark-metrics-v1\""),
          std::string::npos);
      EXPECT_NE(
          sidecar.find("\"simulation_artifact_contract_id\": \"simulation-artifact-v1\""),
          std::string::npos);

      auto benchmark_metrics_path = csv_path;
      benchmark_metrics_path.replace_extension(".benchmark_metrics.json");
      ASSERT_TRUE(std::filesystem::exists(benchmark_metrics_path));
      const auto benchmark_sidecar = readTextFile(benchmark_metrics_path);
      EXPECT_NE(
          benchmark_sidecar.find("\"schema_version\": \"scdat.vacuum_arc.benchmark_metrics.v1\""),
          std::string::npos);
      EXPECT_NE(
          benchmark_sidecar.find("\"contract_id\": \"vacuum-arc-benchmark-metrics-v1\""),
          std::string::npos);

      auto simulation_artifact_path = csv_path;
      simulation_artifact_path.replace_extension(".simulation_artifact.json");
      ASSERT_TRUE(std::filesystem::exists(simulation_artifact_path));
      const auto simulation_artifact = readTextFile(simulation_artifact_path);
      EXPECT_NE(
          simulation_artifact.find("\"schema_version\": \"scdat.simulation_artifact.v1\""),
          std::string::npos);
      EXPECT_NE(
          simulation_artifact.find("\"contract_id\": \"simulation-artifact-v1\""),
          std::string::npos);
  }

TEST(VacuumArcPresetSmokeTest, CollisionStageFallbackCanBeTriggeredByRateLimit)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 20);
    preset.config.enable_pic_mcc_collisions = true;
    preset.config.enable_collision_reaction_breakdown = true;
    preset.config.collision_reconfigure_interval_steps = 1;
    preset.config.collision_neutral_density_floor_m3 = 1.0e22;
    preset.config.collision_ion_density_floor_m3 = 5.0e16;
    preset.config.collision_electron_density_floor_m3 = 5.0e16;
    preset.config.collision_anomaly_event_rate_limit_per_s = 1.0;
    preset.config.collision_fallback_to_background_mcc_on_error = true;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    bool fallback_triggered = false;
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        const auto& status = algorithm.getStatus();
        fallback_triggered = fallback_triggered || status.collision_stage_fallback_triggered;
        EXPECT_TRUE(std::isfinite(status.collision_event_rate_per_s));
        EXPECT_TRUE(std::isfinite(status.collision_effective_neutral_density_m3));
        EXPECT_TRUE(std::isfinite(status.collision_effective_ion_density_m3));
        EXPECT_TRUE(std::isfinite(status.collision_effective_electron_density_m3));
    }

    EXPECT_TRUE(fallback_triggered);
}

TEST(VacuumArcPresetSmokeTest, ResidualMonitoringExportsDiagnosticsAndCanTriggerAlarm)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 24);
    preset.config.enable_residual_monitoring = true;
    preset.config.residual_field_tolerance_v_per_m = 1.0;
    preset.config.residual_particle_current_tolerance_a = 1.0e-18;
    preset.config.residual_circuit_charge_tolerance_c = 1.0e-24;
    preset.config.residual_alarm_consecutive_steps = 1;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    bool alarm_triggered = false;
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        const auto& status = algorithm.getStatus();
        EXPECT_TRUE(std::isfinite(status.field_residual_v_per_m));
        EXPECT_TRUE(std::isfinite(status.particle_residual_a));
        EXPECT_TRUE(std::isfinite(status.circuit_residual_c));
        alarm_triggered = alarm_triggered || status.residual_alarm_active;
    }

    EXPECT_TRUE(alarm_triggered);

    const auto csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_residual_monitoring_stage_a.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));

    SCDAT::Output::DataAnalyzer analyzer;
    const auto data_set = analyzer.loadCsv(csv_path);
    ASSERT_TRUE(data_set.isConsistent());
    ASSERT_TRUE(data_set.hasOnlyFiniteValues());

    EXPECT_TRUE(hasColumn(data_set, "field_residual_v_per_m"));
    EXPECT_TRUE(hasColumn(data_set, "particle_residual_a"));
    EXPECT_TRUE(hasColumn(data_set, "circuit_residual_c"));
    EXPECT_TRUE(hasColumn(data_set, "residual_alarm_active"));
}

TEST(VacuumArcPresetSmokeTest, ResidualMonitoringCalibratedThresholdSuppressesChatter)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 36);
    preset.config.enable_residual_monitoring = true;
    preset.config.enable_residual_ema_filter = true;
    preset.config.residual_ema_alpha = 0.25;
    preset.config.residual_alarm_consecutive_steps = 3;
    preset.config.residual_alarm_clear_consecutive_steps = 3;
    preset.config.residual_field_tolerance_v_per_m = 1.0e12;
    preset.config.residual_particle_current_tolerance_a = 1.0e3;
    preset.config.residual_circuit_charge_tolerance_c = 1.0;

    SurfaceDischargeArcAlgorithm calibration_algorithm;
    ASSERT_TRUE(calibration_algorithm.initialize(preset.config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(calibration_algorithm.advance(preset.time_step_s));
    }

    const auto calibration_csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_residual_monitoring_stage_b_calibration.csv";
    ASSERT_TRUE(calibration_algorithm.exportResults(calibration_csv_path));

    SCDAT::Output::DataAnalyzer analyzer;
    const auto calibration_data_set = analyzer.loadCsv(calibration_csv_path);
    ASSERT_TRUE(calibration_data_set.isConsistent());
    ASSERT_TRUE(calibration_data_set.hasOnlyFiniteValues());
    EXPECT_TRUE(hasColumn(calibration_data_set, "field_residual_filtered_v_per_m"));
    EXPECT_TRUE(hasColumn(calibration_data_set, "particle_residual_filtered_a"));
    EXPECT_TRUE(hasColumn(calibration_data_set, "circuit_residual_filtered_c"));

    const double calibrated_field_tolerance =
        std::max(1.0, 2.0 * percentileAbsValue(calibration_data_set,
                                               "field_residual_filtered_v_per_m", 0.95));
    const double calibrated_particle_tolerance =
        std::max(1.0e-12, 2.0 * percentileAbsValue(calibration_data_set,
                                                   "particle_residual_filtered_a", 0.95));
    const double calibrated_circuit_tolerance =
        std::max(1.0e-18, 2.0 * percentileAbsValue(calibration_data_set,
                                                   "circuit_residual_filtered_c", 0.95));

    auto calibrated_config = preset.config;
    calibrated_config.residual_field_tolerance_v_per_m = calibrated_field_tolerance;
    calibrated_config.residual_particle_current_tolerance_a = calibrated_particle_tolerance;
    calibrated_config.residual_circuit_charge_tolerance_c = calibrated_circuit_tolerance;

    SurfaceDischargeArcAlgorithm calibrated_algorithm;
    ASSERT_TRUE(calibrated_algorithm.initialize(calibrated_config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(calibrated_algorithm.advance(preset.time_step_s));
    }

    const auto calibrated_csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_residual_monitoring_stage_b_calibrated.csv";
    ASSERT_TRUE(calibrated_algorithm.exportResults(calibrated_csv_path));
    const auto calibrated_data_set = analyzer.loadCsv(calibrated_csv_path);
    ASSERT_TRUE(calibrated_data_set.isConsistent());
    ASSERT_TRUE(calibrated_data_set.hasOnlyFiniteValues());

    EXPECT_TRUE(hasColumn(calibrated_data_set, "field_residual_filtered_v_per_m"));
    EXPECT_TRUE(hasColumn(calibrated_data_set, "particle_residual_filtered_a"));
    EXPECT_TRUE(hasColumn(calibrated_data_set, "circuit_residual_filtered_c"));
    EXPECT_TRUE(hasColumn(calibrated_data_set, "residual_alarm_counter"));
    EXPECT_TRUE(hasColumn(calibrated_data_set, "residual_alarm_clear_counter"));

    const std::size_t calibrated_alarm_count =
        countValuesAbove(calibrated_data_set, "residual_alarm_active", 0.5);
    EXPECT_LE(calibrated_alarm_count, 2u);

    auto strict_config = calibrated_config;
    strict_config.residual_alarm_consecutive_steps = 1;
    strict_config.residual_alarm_clear_consecutive_steps = 1;
    strict_config.residual_field_tolerance_v_per_m =
        std::max(1.0e-12, calibrated_field_tolerance * 1.0e-6);
    strict_config.residual_particle_current_tolerance_a =
        std::max(1.0e-18, calibrated_particle_tolerance * 1.0e-6);
    strict_config.residual_circuit_charge_tolerance_c =
        std::max(1.0e-24, calibrated_circuit_tolerance * 1.0e-6);

    SurfaceDischargeArcAlgorithm strict_algorithm;
    ASSERT_TRUE(strict_algorithm.initialize(strict_config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(strict_algorithm.advance(preset.time_step_s));
    }

    const auto strict_csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_residual_monitoring_stage_b_strict.csv";
    ASSERT_TRUE(strict_algorithm.exportResults(strict_csv_path));
    const auto strict_data_set = analyzer.loadCsv(strict_csv_path);
    ASSERT_TRUE(strict_data_set.isConsistent());
    ASSERT_TRUE(strict_data_set.hasOnlyFiniteValues());

    const std::size_t strict_alarm_count =
        countValuesAbove(strict_data_set, "residual_alarm_active", 0.5);
    EXPECT_GT(strict_alarm_count, 0u);
}

TEST(VacuumArcPresetSmokeTest, OutgassingAndSurfaceLoadDiagnosticsAreExported)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("triple_junction_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 20);
    preset.config.surface_load_model =
        SCDAT::Toolkit::VacuumArc::SurfaceCircuitLoadModel::RlBranch;
    preset.config.surface_load_resistance_ohm = 1500.0;
    preset.config.surface_load_inductance_h = 2.0e-6;
    preset.config.surface_load_current_limit_a = 0.2;
    preset.config.enable_neutral_outgassing_feedback = true;
    preset.config.neutral_outgassing_gain_m3_per_a = 2.0e13;
    preset.config.neutral_outgassing_relaxation_per_s = 5.0e7;
    preset.config.neutral_outgassing_max_density_boost_m3 = 1.0e14;
    preset.config.enable_neutral_outgassing_reinjection = true;
    preset.config.neutral_outgassing_reinjection_gain = 1.0e8;
    preset.config.neutral_outgassing_reinjection_macro_weight = 1.0e5;
    preset.config.neutral_outgassing_reinjection_max_particles_per_step = 128;
    preset.config.neutral_outgassing_reinjection_speed_m_per_s = 2.0e3;
    preset.config.neutral_outgassing_reinjection_mass_amu = 28.0;
    preset.config.enable_secondary_electron_emission = true;
    preset.config.secondary_electron_yield_per_ion = 0.2;
    preset.config.secondary_electron_speed_m_per_s = 1.0e6;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_outgassing_load_diagnostics.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));

    SCDAT::Output::DataAnalyzer analyzer;
    const auto data_set = analyzer.loadCsv(csv_path);

    ASSERT_TRUE(data_set.isConsistent());
    ASSERT_TRUE(data_set.hasOnlyFiniteValues());
    EXPECT_TRUE(hasColumn(data_set, "surface_load_branch_current_a"));
    EXPECT_TRUE(hasColumn(data_set, "charge_conservation_error_c"));
    EXPECT_TRUE(hasColumn(data_set, "surface_relative_potential_v"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_capacitor_energy_j"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_drive_power_w"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_leak_power_w"));
    EXPECT_TRUE(hasColumn(data_set, "surface_load_branch_power_w"));
    EXPECT_TRUE(hasColumn(data_set, "surface_load_resistive_power_w"));
    EXPECT_TRUE(hasColumn(data_set, "surface_load_inductive_energy_j"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_power_balance_error_w"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_drive_energy_j"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_leak_dissipated_energy_j"));
    EXPECT_TRUE(hasColumn(data_set, "surface_load_dissipated_energy_j"));
    EXPECT_TRUE(hasColumn(data_set, "surface_circuit_energy_balance_error_j"));
    EXPECT_TRUE(hasColumn(data_set, "collision_effective_ion_density_m3"));
    EXPECT_TRUE(hasColumn(data_set, "collision_effective_electron_density_m3"));
    EXPECT_TRUE(hasColumn(data_set, "neutral_outgassing_density_boost_m3"));
    EXPECT_TRUE(hasColumn(data_set, "neutral_outgassing_reinjected_particles_step"));
    EXPECT_TRUE(hasColumn(data_set, "neutral_outgassing_reinjection_rate_per_s"));
    EXPECT_TRUE(hasColumn(data_set, "stability_substeps_used"));
    EXPECT_TRUE(hasColumn(data_set, "stability_rollbacks"));
    EXPECT_TRUE(hasColumn(data_set, "stability_effective_substep_s"));
    EXPECT_TRUE(hasColumn(data_set, "stability_adaptive_scale"));
    EXPECT_TRUE(hasColumn(data_set, "stability_adaptive_reductions"));
    EXPECT_TRUE(hasColumn(data_set, "stability_anomaly_isolation_triggered"));
    EXPECT_TRUE(hasColumn(data_set, "stability_anomaly_isolation_events"));

    EXPECT_GT(maxAbsValue(data_set, "surface_load_branch_current_a"), 0.0);
    EXPECT_GT(maxAbsValue(data_set, "neutral_outgassing_density_boost_m3"), 0.0);
    EXPECT_GT(maxAbsValue(data_set, "neutral_outgassing_reinjected_particles_step"), 0.0);
    EXPECT_GT(maxAbsValue(data_set, "surface_circuit_capacitor_energy_j"), 0.0);
    EXPECT_GE(lastValue(data_set, "surface_circuit_leak_dissipated_energy_j"), 0.0);
    EXPECT_GE(lastValue(data_set, "surface_load_dissipated_energy_j"), 0.0);
}

TEST(VacuumArcPresetSmokeTest, SurfaceCircuitLoadModelsExposePhaseAndEnergyDiagnostics)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("triple_junction_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 24);
    preset.config.max_surface_potential_step_v = 1.0e6;
    preset.config.surface_load_resistance_ohm = 1500.0;
    preset.config.surface_load_inductance_h = 2.0e-6;
    preset.config.surface_load_current_limit_a = 0.2;

    SCDAT::Output::DataAnalyzer analyzer;

    auto run_with_load_model = [&](SurfaceCircuitLoadModel model, const std::string& suffix) {
        auto model_preset = preset;
        model_preset.config.surface_load_model = model;

        SurfaceDischargeArcAlgorithm algorithm;
        EXPECT_TRUE(algorithm.initialize(model_preset.config));
        for (std::size_t step = 0; step < model_preset.steps; ++step)
        {
            EXPECT_TRUE(algorithm.advance(model_preset.time_step_s));
        }

        const auto csv_path = std::filesystem::temp_directory_path() /
                              ("vacuum_arc_surface_circuit_alignment_" + suffix + ".csv");
        EXPECT_TRUE(algorithm.exportResults(csv_path));
        return analyzer.loadCsv(csv_path);
    };

    const auto legacy_data =
        run_with_load_model(SurfaceCircuitLoadModel::LegacyLeakageCapacitor, "legacy");
    const auto resistive_data =
        run_with_load_model(SurfaceCircuitLoadModel::ResistiveShunt, "resistive");
    const auto rl_data =
        run_with_load_model(SurfaceCircuitLoadModel::RlBranch, "rl");

    auto verify_common_diagnostics = [&](const SCDAT::Output::ColumnarDataSet& data_set) {
        ASSERT_TRUE(data_set.isConsistent());
        ASSERT_TRUE(data_set.hasOnlyFiniteValues());
        EXPECT_TRUE(hasColumn(data_set, "surface_relative_potential_v"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_capacitor_energy_j"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_drive_power_w"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_leak_power_w"));
        EXPECT_TRUE(hasColumn(data_set, "surface_load_branch_power_w"));
        EXPECT_TRUE(hasColumn(data_set, "surface_load_resistive_power_w"));
        EXPECT_TRUE(hasColumn(data_set, "surface_load_inductive_energy_j"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_power_balance_error_w"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_drive_energy_j"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_leak_dissipated_energy_j"));
        EXPECT_TRUE(hasColumn(data_set, "surface_load_dissipated_energy_j"));
        EXPECT_TRUE(hasColumn(data_set, "surface_circuit_energy_balance_error_j"));

        EXPECT_GE(lastValue(data_set, "surface_circuit_capacitor_energy_j"), 0.0);
        EXPECT_GE(lastValue(data_set, "surface_circuit_leak_dissipated_energy_j"), 0.0);
        EXPECT_GE(lastValue(data_set, "surface_load_dissipated_energy_j"), 0.0);

        const auto& leak_dissipated_energy =
            data_set.scalar_series.at("surface_circuit_leak_dissipated_energy_j");
        const auto& load_dissipated_energy =
            data_set.scalar_series.at("surface_load_dissipated_energy_j");
        ASSERT_EQ(leak_dissipated_energy.size(), load_dissipated_energy.size());

        for (std::size_t i = 1; i < leak_dissipated_energy.size(); ++i)
        {
            EXPECT_GE(leak_dissipated_energy[i] + 1.0e-18, leak_dissipated_energy[i - 1]);
            EXPECT_GE(load_dissipated_energy[i] + 1.0e-18, load_dissipated_energy[i - 1]);
        }
    };

    verify_common_diagnostics(legacy_data);
    verify_common_diagnostics(resistive_data);
    verify_common_diagnostics(rl_data);

    EXPECT_LE(maxAbsValue(legacy_data, "surface_load_branch_current_a"), 1.0e-12);
    EXPECT_LE(maxAbsValue(legacy_data, "surface_load_branch_power_w"), 1.0e-12);
    EXPECT_LE(maxAbsValue(legacy_data, "surface_load_inductive_energy_j"), 1.0e-12);

    const auto& resistive_relative_potential =
        resistive_data.scalar_series.at("surface_relative_potential_v");
    const auto& resistive_load_current =
        resistive_data.scalar_series.at("surface_load_branch_current_a");
    ASSERT_EQ(resistive_relative_potential.size(), resistive_load_current.size());

    double resistive_ohmic_residual_v = 0.0;
    for (std::size_t i = 0; i < resistive_relative_potential.size(); ++i)
    {
        const double residual_v =
            resistive_relative_potential[i] -
            preset.config.surface_load_resistance_ohm * resistive_load_current[i];
        resistive_ohmic_residual_v = std::max(resistive_ohmic_residual_v, std::abs(residual_v));
    }
    const double resistive_voltage_scale =
        std::max(1.0e-12, maxAbsValue(resistive_data, "surface_relative_potential_v"));
    EXPECT_LE(resistive_ohmic_residual_v / resistive_voltage_scale, 5.0e-3);

    const auto& rl_relative_potential = rl_data.scalar_series.at("surface_relative_potential_v");
    const auto& rl_load_current = rl_data.scalar_series.at("surface_load_branch_current_a");
    ASSERT_EQ(rl_relative_potential.size(), rl_load_current.size());

    double rl_ohmic_residual_v = 0.0;
    for (std::size_t i = 0; i < rl_relative_potential.size(); ++i)
    {
        const double residual_v =
            rl_relative_potential[i] -
            preset.config.surface_load_resistance_ohm * rl_load_current[i];
        rl_ohmic_residual_v = std::max(rl_ohmic_residual_v, std::abs(residual_v));
    }
    const double rl_voltage_scale =
        std::max(1.0e-12, maxAbsValue(rl_data, "surface_relative_potential_v"));
    EXPECT_GT(maxAbsValue(rl_data, "surface_load_inductive_energy_j"), 0.0);
    EXPECT_GT(rl_ohmic_residual_v / rl_voltage_scale, 1.0e-4);
}

TEST(VacuumArcPresetSmokeTest, AdaptiveSubstepIsolationDiagnosticsRemainFinite)
{
    VacuumArcScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset("microgap_flashover", preset));

    preset.steps = std::max<std::size_t>(preset.steps, 20);
    preset.config.enable_residual_monitoring = true;
    preset.config.residual_alarm_consecutive_steps = 1;
    preset.config.residual_alarm_clear_consecutive_steps = 1;
    preset.config.residual_field_tolerance_v_per_m = 1.0;
    preset.config.residual_particle_current_tolerance_a = 1.0e-18;
    preset.config.residual_circuit_charge_tolerance_c = 1.0e-24;
    preset.config.enable_stability_rollback = true;
    preset.config.max_stability_rollbacks = 1;
    preset.config.enable_adaptive_internal_timestep = true;
    preset.config.adaptive_substep_shrink_factor = 0.5;
    preset.config.adaptive_substep_recovery_factor = 1.2;
    preset.config.adaptive_substep_min_scale = 0.125;
    preset.config.anomaly_current_density_limit_a_per_m2 = 1.0e-3;
    preset.config.anomaly_surface_potential_limit_v = 1.0e-3;
    preset.config.max_surface_potential_step_v = 1.0e6;

    SurfaceDischargeArcAlgorithm algorithm;
    ASSERT_TRUE(algorithm.initialize(preset.config));

    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "vacuum_arc_stability_stage_va008.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));

    SCDAT::Output::DataAnalyzer analyzer;
    const auto data_set = analyzer.loadCsv(csv_path);

    ASSERT_TRUE(data_set.isConsistent());
    ASSERT_TRUE(data_set.hasOnlyFiniteValues());
    EXPECT_TRUE(hasColumn(data_set, "stability_substeps_used"));
    EXPECT_TRUE(hasColumn(data_set, "stability_rollbacks"));
    EXPECT_TRUE(hasColumn(data_set, "stability_effective_substep_s"));
    EXPECT_TRUE(hasColumn(data_set, "stability_adaptive_scale"));
    EXPECT_TRUE(hasColumn(data_set, "stability_adaptive_reductions"));
    EXPECT_TRUE(hasColumn(data_set, "stability_anomaly_isolation_triggered"));
    EXPECT_TRUE(hasColumn(data_set, "stability_anomaly_isolation_events"));

    EXPECT_GT(maxAbsValue(data_set, "stability_effective_substep_s"), 0.0);
    EXPECT_GT(maxAbsValue(data_set, "stability_adaptive_reductions"), 0.0);
    EXPECT_GT(maxAbsValue(data_set, "stability_anomaly_isolation_events"), 0.0);
    EXPECT_GT(countValuesAbove(data_set, "stability_anomaly_isolation_triggered", 0.5), 0u);

    const auto& effective_substep_series =
        data_set.scalar_series.at("stability_effective_substep_s");
    ASSERT_FALSE(effective_substep_series.empty());
    double min_effective_substep_s = effective_substep_series.front();
    for (double sample : effective_substep_series)
    {
        min_effective_substep_s = std::min(min_effective_substep_s, sample);
    }
    EXPECT_LE(min_effective_substep_s, preset.time_step_s);
}

TEST(VacuumArcPresetSmokeTest, ArcPicBenchmarkGateTargetCasesCoverCoreMetrics)
{
    SCDAT::Output::DataAnalyzer analyzer;

    auto run_case_with_mode = [&](const std::string& preset_name,
                                  ArcPicAlignmentMode alignment_mode,
                                  const std::string& suffix) {
        VacuumArcScenarioPreset preset;
        EXPECT_TRUE(SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset(preset_name, preset));
        preset.config.alignment_mode = alignment_mode;
        preset.steps = std::max<std::size_t>(preset.steps, 16);
        preset.config.enable_pic_mcc_collisions = true;
        preset.config.collision_neutral_density_floor_m3 = 8.0e20;

        SurfaceDischargeArcAlgorithm algorithm;
        EXPECT_TRUE(algorithm.initialize(preset.config));
        for (std::size_t step = 0; step < preset.steps; ++step)
        {
            EXPECT_TRUE(algorithm.advance(preset.time_step_s));
        }

        const auto csv_path = std::filesystem::temp_directory_path() /
                              ("vacuum_arc_arcpic_gate_" + preset_name + "_" + suffix + ".csv");
        EXPECT_TRUE(algorithm.exportResults(csv_path));
        return analyzer.loadCsv(csv_path);
    };

    const std::vector<std::string> target_cases{"microgap_flashover", "triple_junction_flashover",
                                                 "restrike_recovery"};

    for (const auto& case_name : target_cases)
    {
          const auto aligned_data = run_case_with_mode(case_name, ArcPicAlignmentMode::ArcPicAligned,
                                                       "aligned");
          const auto legacy_data = run_case_with_mode(case_name, ArcPicAlignmentMode::LegacyBaseline,
                                                      "legacy");

        for (const auto* data_set : {&aligned_data, &legacy_data})
        {
            ASSERT_TRUE(data_set->isConsistent());
            ASSERT_TRUE(data_set->hasOnlyFiniteValues());

            const double peak_current_a = maxAbsValue(*data_set, "discharge_current_a");
            const double peak_current_density = maxAbsValue(*data_set, "current_density_a_per_m2");
            const double final_surface_potential_v = lastValue(*data_set, "surface_potential_v");
            const double collision_event_sum = sumValue(*data_set, "collision_events_step");

            EXPECT_GT(peak_current_a, 0.0);
            EXPECT_GT(peak_current_density, 0.0);
            EXPECT_TRUE(std::isfinite(final_surface_potential_v));
            EXPECT_GE(collision_event_sum, 0.0);
            EXPECT_TRUE(hasColumn(*data_set, "collision_effective_neutral_density_m3"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_effective_ion_density_m3"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_effective_electron_density_m3"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_ionization_fraction_step"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_excitation_fraction_step"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_charge_exchange_fraction_step"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_reaction_weighted_emission_feedback_step"));
            EXPECT_TRUE(hasColumn(*data_set, "collision_reaction_weighted_channel_feedback_step"));
            EXPECT_TRUE(hasColumn(*data_set, "field_residual_v_per_m"));
            EXPECT_TRUE(hasColumn(*data_set, "particle_residual_a"));
            EXPECT_TRUE(hasColumn(*data_set, "circuit_residual_c"));
            EXPECT_TRUE(hasColumn(*data_set, "neutral_outgassing_reinjected_particles_step"));
            EXPECT_TRUE(hasColumn(*data_set, "neutral_outgassing_reinjection_rate_per_s"));
            EXPECT_TRUE(hasColumn(*data_set, "field_residual_filtered_v_per_m"));
            EXPECT_TRUE(hasColumn(*data_set, "particle_residual_filtered_a"));
            EXPECT_TRUE(hasColumn(*data_set, "circuit_residual_filtered_c"));
            EXPECT_TRUE(hasColumn(*data_set, "residual_alarm_active"));
            EXPECT_TRUE(hasColumn(*data_set, "residual_alarm_counter"));
            EXPECT_TRUE(hasColumn(*data_set, "residual_alarm_clear_counter"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_relative_potential_v"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_capacitor_energy_j"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_drive_power_w"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_leak_power_w"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_load_branch_power_w"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_load_resistive_power_w"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_load_inductive_energy_j"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_power_balance_error_w"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_drive_energy_j"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_leak_dissipated_energy_j"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_load_dissipated_energy_j"));
            EXPECT_TRUE(hasColumn(*data_set, "surface_circuit_energy_balance_error_j"));
            EXPECT_TRUE(hasColumn(*data_set, "stability_effective_substep_s"));
            EXPECT_TRUE(hasColumn(*data_set, "stability_adaptive_scale"));
            EXPECT_TRUE(hasColumn(*data_set, "stability_adaptive_reductions"));
            EXPECT_TRUE(hasColumn(*data_set, "stability_anomaly_isolation_triggered"));
            EXPECT_TRUE(hasColumn(*data_set, "stability_anomaly_isolation_events"));
        }

        const double aligned_peak_current = maxAbsValue(aligned_data, "discharge_current_a");
        const double legacy_peak_current = maxAbsValue(legacy_data, "discharge_current_a");
        const double aligned_peak_density = maxAbsValue(aligned_data, "current_density_a_per_m2");
        const double legacy_peak_density = maxAbsValue(legacy_data, "current_density_a_per_m2");
          const double aligned_final_potential = lastValue(aligned_data, "surface_potential_v");
          const double legacy_final_potential = lastValue(legacy_data, "surface_potential_v");

        const double current_delta_ratio =
            std::abs(aligned_peak_current - legacy_peak_current) /
            std::max(1.0e-12, std::abs(aligned_peak_current));
          const double density_delta_ratio =
              std::abs(aligned_peak_density - legacy_peak_density) /
            std::max(1.0e-12, std::abs(aligned_peak_density));

        EXPECT_LE(current_delta_ratio, 0.75);
        EXPECT_LE(density_delta_ratio, 0.75);
        EXPECT_LE(std::abs(aligned_final_potential - legacy_final_potential), 15.0);
    }
}
