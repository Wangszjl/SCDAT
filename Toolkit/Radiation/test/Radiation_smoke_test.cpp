#include "RadiationCases.h"
#include "RadiationDoseAlgorithm.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

using SCDAT::Toolkit::Radiation::RadiationDoseAlgorithm;
using SCDAT::Toolkit::Radiation::RadiationPhysicsList;
using SCDAT::Toolkit::Radiation::RadiationScenarioPreset;
using SCDAT::Toolkit::Radiation::RadiationStatus;

namespace
{
std::string readTextFile(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::vector<double> normalizedDepthDoseProfile(const RadiationStatus& status,
                                               double layer_thickness_m)
{
    std::vector<double> profile;
    profile.reserve(status.layers.size());

    double total_energy_j_per_m2 = 0.0;
    for (const auto& layer : status.layers)
    {
        const double layer_energy_j_per_m2 =
            std::max(0.0, layer.deposited_energy_j_per_m3 * layer_thickness_m);
        profile.push_back(layer_energy_j_per_m2);
        total_energy_j_per_m2 += layer_energy_j_per_m2;
    }

    const double normalizer = std::max(total_energy_j_per_m2, 1.0e-30);
    for (auto& value : profile)
    {
        value /= normalizer;
    }
    return profile;
}

double normalizedProfileRmse(const RadiationStatus& candidate,
                             const RadiationStatus& reference,
                             double layer_thickness_m)
{
    const auto candidate_profile = normalizedDepthDoseProfile(candidate, layer_thickness_m);
    const auto reference_profile = normalizedDepthDoseProfile(reference, layer_thickness_m);

    if (candidate_profile.size() != reference_profile.size() || candidate_profile.empty())
    {
        return std::numeric_limits<double>::infinity();
    }

    double mse = 0.0;
    for (std::size_t i = 0; i < candidate_profile.size(); ++i)
    {
        const double diff = candidate_profile[i] - reference_profile[i];
        mse += diff * diff;
    }

    mse /= static_cast<double>(candidate_profile.size());
    return std::sqrt(mse);
}
}

TEST(RadiationSmokeTest, InitializesAdvancesAndExports)
{
    RadiationDoseAlgorithm algorithm;
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    ASSERT_TRUE(algorithm.initialize(preset.config));
    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto& status = algorithm.getStatus();
    ASSERT_EQ(status.layers.size(), preset.config.layers);
    EXPECT_GT(status.incident_energy_j_per_m2, 0.0);
    EXPECT_GT(status.deposited_energy_j_per_m2, 0.0);
    EXPECT_LE(status.energy_conservation_error, 1.0e-10);

    for (const auto& layer : status.layers)
    {
        EXPECT_TRUE(std::isfinite(layer.deposited_energy_j_per_m3));
        EXPECT_TRUE(std::isfinite(layer.dose_gy));
        EXPECT_GE(layer.deposited_energy_j_per_m3, 0.0);
        EXPECT_GE(layer.dose_gy, 0.0);
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "radiation_smoke.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));
}

TEST(RadiationSmokeTest, MonteCarloTransportExportsTrackCsv)
{
    RadiationDoseAlgorithm algorithm;
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    preset.config.enable_monte_carlo_transport = true;
    preset.config.monte_carlo_histories_per_step = 64;
    preset.config.monte_carlo_max_steps_per_track = 96;
    preset.config.monte_carlo_max_recorded_points = 20000;
    preset.config.monte_carlo_seed = 20260403u;
    preset.steps = 4;

    ASSERT_TRUE(algorithm.initialize(preset.config));
    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(algorithm.advance(preset.time_step_s));
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "radiation_mc_tracks_smoke.csv";
    ASSERT_TRUE(algorithm.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));

    auto track_path = csv_path;
    track_path.replace_extension();
    track_path += ".tracks.csv";
    EXPECT_TRUE(std::filesystem::exists(track_path));
    auto deposition_history_path = csv_path;
    deposition_history_path.replace_extension();
    deposition_history_path += ".deposition_history.json";
    auto process_history_path = csv_path;
    process_history_path.replace_extension();
    process_history_path += ".process_history.json";
    auto transport_benchmark_path = csv_path;
    transport_benchmark_path.replace_extension();
    transport_benchmark_path += ".radiation_transport_benchmark.json";
    auto simulation_artifact_path = csv_path;
    simulation_artifact_path.replace_extension();
    simulation_artifact_path += ".simulation_artifact.json";
    EXPECT_TRUE(std::filesystem::exists(deposition_history_path));
    EXPECT_TRUE(std::filesystem::exists(process_history_path));
    EXPECT_TRUE(std::filesystem::exists(transport_benchmark_path));
    EXPECT_TRUE(std::filesystem::exists(simulation_artifact_path));

    std::ifstream input(track_path);
    ASSERT_TRUE(input.is_open());
    std::string line;
    std::size_t line_count = 0;
    while (std::getline(input, line))
    {
        ++line_count;
    }
    EXPECT_GT(line_count, 1u);

    const auto sidecar = readTextFile(csv_path.string() + ".metadata.json");
    EXPECT_NE(
        sidecar.find("\"deposition_record_contract_id\": \"geant4-aligned-deposition-record-v1\""),
        std::string::npos);
    EXPECT_NE(
        sidecar.find("\"process_history_contract_id\": \"geant4-aligned-process-history-v1\""),
        std::string::npos);
    EXPECT_NE(
        sidecar.find("\"deposition_history_artifact_path\":"),
        std::string::npos);
    EXPECT_NE(
        sidecar.find("\"radiation_transport_benchmark_contract_id\": \"radiation-transport-benchmark-v1\""),
        std::string::npos);
    EXPECT_NE(
        sidecar.find("\"radiation_transport_benchmark_artifact_path\":"),
        std::string::npos);
    EXPECT_NE(
        sidecar.find("\"simulation_artifact_contract_id\": \"simulation-artifact-v1\""),
        std::string::npos);
    EXPECT_NE(
        sidecar.find("\"simulation_artifact_path\":"),
        std::string::npos);

    const auto transport_benchmark = readTextFile(transport_benchmark_path);
    EXPECT_NE(
        transport_benchmark.find("\"schema_version\": \"scdat.radiation_transport_benchmark.v1\""),
        std::string::npos);
    EXPECT_NE(
        transport_benchmark.find("\"contract_id\": \"radiation-transport-benchmark-v1\""),
        std::string::npos);
    const auto simulation_artifact = readTextFile(simulation_artifact_path);
    EXPECT_NE(
        simulation_artifact.find("\"schema_version\": \"scdat.simulation_artifact.v1\""),
        std::string::npos);
    EXPECT_NE(
        simulation_artifact.find("\"contract_id\": \"simulation-artifact-v1\""),
        std::string::npos);
}

TEST(RadiationSmokeTest, TrackSchemaV2StatisticsAreReproducible)
{
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    preset.steps = 5;
    preset.config.enable_monte_carlo_transport = true;
    preset.config.monte_carlo_histories_per_step = 96;
    preset.config.monte_carlo_max_steps_per_track = 96;
    preset.config.monte_carlo_max_recorded_points = 30000;
    preset.config.monte_carlo_seed = 20260403u;
    preset.config.electron_process.enable = true;

    RadiationDoseAlgorithm algorithm_a;
    RadiationDoseAlgorithm algorithm_b;
    ASSERT_TRUE(algorithm_a.initialize(preset.config));
    ASSERT_TRUE(algorithm_b.initialize(preset.config));

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(algorithm_a.advance(preset.time_step_s));
        ASSERT_TRUE(algorithm_b.advance(preset.time_step_s));
    }

    const auto& status_a = algorithm_a.getStatus();
    const auto& status_b = algorithm_b.getStatus();

    EXPECT_EQ(status_a.track_count, status_b.track_count);
    EXPECT_EQ(status_a.track_step_count, status_b.track_step_count);
    EXPECT_GT(status_a.track_count, 0u);
    EXPECT_GT(status_a.track_step_count, 0u);
    EXPECT_GT(status_a.track_total_path_length_m, 0.0);
    EXPECT_GT(status_a.track_mean_step_length_m, 0.0);
    EXPECT_GE(status_a.track_max_depth_m, 0.0);
    EXPECT_LE(status_a.track_max_depth_m, preset.config.thickness_m + 1.0e-12);

    EXPECT_NEAR(status_a.track_total_path_length_m, status_b.track_total_path_length_m, 1.0e-12);
    EXPECT_NEAR(status_a.track_mean_step_length_m, status_b.track_mean_step_length_m, 1.0e-12);
    EXPECT_NEAR(status_a.track_mean_scattering_sigma_rad,
                status_b.track_mean_scattering_sigma_rad, 1.0e-12);
    EXPECT_NEAR(status_a.track_secondary_event_count, status_b.track_secondary_event_count, 1.0e-12);
    EXPECT_NEAR(status_a.track_secondary_energy_j_per_m2,
                status_b.track_secondary_energy_j_per_m2, 1.0e-12);
    EXPECT_NEAR(status_a.track_secondary_yield_per_primary,
                status_b.track_secondary_yield_per_primary, 1.0e-12);
    EXPECT_NEAR(status_a.track_mean_terminal_energy_ev,
                status_b.track_mean_terminal_energy_ev, 1.0e-9);

    const auto csv_path = std::filesystem::temp_directory_path() / "radiation_track_schema_v2_smoke.csv";
    ASSERT_TRUE(algorithm_a.exportResults(csv_path));

    auto track_path = csv_path;
    track_path.replace_extension();
    track_path += ".tracks.csv";
    std::ifstream input(track_path);
    ASSERT_TRUE(input.is_open());

    std::string header;
    ASSERT_TRUE(std::getline(input, header));
    EXPECT_NE(header.find("secondary_energy_ev"), std::string::npos);
    EXPECT_NE(header.find("scattering_sigma_rad"), std::string::npos);
    EXPECT_NE(header.find("step_length_m"), std::string::npos);
}

TEST(RadiationSmokeTest, PhysicsListSwitchChangesDoseTrend)
{
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    preset.config.enable_monte_carlo_transport = false;
    preset.steps = 6;

    auto standard_preset = preset;
    standard_preset.config.physics_list = RadiationPhysicsList::Geant4EmStandard;

    auto livermore_preset = preset;
    livermore_preset.config.physics_list = RadiationPhysicsList::Geant4EmLivermore;

    RadiationDoseAlgorithm standard_algorithm;
    RadiationDoseAlgorithm livermore_algorithm;
    ASSERT_TRUE(standard_algorithm.initialize(standard_preset.config));
    ASSERT_TRUE(livermore_algorithm.initialize(livermore_preset.config));

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(standard_algorithm.advance(preset.time_step_s));
        ASSERT_TRUE(livermore_algorithm.advance(preset.time_step_s));
    }

    const auto& standard_status = standard_algorithm.getStatus();
    const auto& livermore_status = livermore_algorithm.getStatus();

    EXPECT_GT(standard_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(livermore_status.deposited_energy_j_per_m2, 0.0);

    ASSERT_FALSE(standard_status.layers.empty());
    ASSERT_EQ(standard_status.layers.size(), livermore_status.layers.size());

    const double layer_thickness_m =
        standard_preset.config.thickness_m / static_cast<double>(standard_preset.config.layers);
    const double standard_surface_fraction =
        standard_status.layers.front().deposited_energy_j_per_m3 * layer_thickness_m /
        std::max(standard_status.deposited_energy_j_per_m2, 1.0e-18);
    const double livermore_surface_fraction =
        livermore_status.layers.front().deposited_energy_j_per_m3 * layer_thickness_m /
        std::max(livermore_status.deposited_energy_j_per_m2, 1.0e-18);

    EXPECT_GT(
        std::abs(livermore_surface_fraction - standard_surface_fraction),
        std::abs(standard_surface_fraction) * 1.0e-3 + 1.0e-8);
}

TEST(RadiationSmokeTest, ElectronProcessEnhancementImprovesDepthDoseCurveAlignment)
{
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    preset.steps = 5;
    preset.config.enable_monte_carlo_transport = false;
    preset.config.electron_process.enable = true;
    preset.config.electron_process.inelastic_energy_loss_scale = 1.12;
    preset.config.electron_process.elastic_scattering_scale = 1.15;
    preset.config.electron_process.secondary_energy_fraction = 0.16;
    preset.config.electron_process.secondary_depth_attenuation_fraction = 0.22;

    auto reference_preset = preset;
    reference_preset.config.enable_monte_carlo_transport = true;
    reference_preset.config.monte_carlo_histories_per_step = 384;
    reference_preset.config.monte_carlo_max_steps_per_track = 128;
    reference_preset.config.monte_carlo_max_recorded_points = 40000;
    reference_preset.config.monte_carlo_seed = 20260403u;

    auto baseline_preset = preset;
    baseline_preset.config.electron_process.enable = false;

    auto enhanced_preset = preset;

    RadiationDoseAlgorithm reference_algorithm;
    RadiationDoseAlgorithm baseline_algorithm;
    RadiationDoseAlgorithm enhanced_algorithm;

    ASSERT_TRUE(reference_algorithm.initialize(reference_preset.config));
    ASSERT_TRUE(baseline_algorithm.initialize(baseline_preset.config));
    ASSERT_TRUE(enhanced_algorithm.initialize(enhanced_preset.config));

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(reference_algorithm.advance(reference_preset.time_step_s));
        ASSERT_TRUE(baseline_algorithm.advance(baseline_preset.time_step_s));
        ASSERT_TRUE(enhanced_algorithm.advance(enhanced_preset.time_step_s));
    }

    const auto& reference_status = reference_algorithm.getStatus();
    const auto& baseline_status = baseline_algorithm.getStatus();
    const auto& enhanced_status = enhanced_algorithm.getStatus();

    ASSERT_EQ(reference_status.layers.size(), baseline_status.layers.size());
    ASSERT_EQ(reference_status.layers.size(), enhanced_status.layers.size());

    const double layer_thickness_m =
        preset.config.thickness_m / static_cast<double>(preset.config.layers);
    const double baseline_rmse =
        normalizedProfileRmse(baseline_status, reference_status, layer_thickness_m);
    const double enhanced_rmse =
        normalizedProfileRmse(enhanced_status, reference_status, layer_thickness_m);

    EXPECT_GT(baseline_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.electron_process_secondary_event_count, 0.0);
    EXPECT_GT(enhanced_status.electron_process_secondary_energy_j_per_m2, 0.0);
    EXPECT_LT(enhanced_rmse, baseline_rmse);
    EXPECT_LT(enhanced_rmse, 0.25);
}

TEST(RadiationSmokeTest, ProtonProcessEnhancementMeetsDepthDoseCurveAlignmentTarget)
{
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("meo_proton_harness_dose", preset));

    preset.steps = 6;
    preset.config.enable_monte_carlo_transport = false;
    preset.config.proton_process.enable = true;
    preset.config.proton_process.nuclear_reaction_energy_scale = 1.06;
    preset.config.proton_process.elastic_scattering_scale = 1.10;
    preset.config.proton_process.secondary_energy_fraction = 0.06;
    preset.config.proton_process.secondary_depth_buildup_fraction = 0.55;

    auto reference_preset = preset;
    reference_preset.config.enable_monte_carlo_transport = true;
    reference_preset.config.monte_carlo_histories_per_step = 320;
    reference_preset.config.monte_carlo_max_steps_per_track = 128;
    reference_preset.config.monte_carlo_max_recorded_points = 40000;
    reference_preset.config.monte_carlo_seed = 20260403u;

    auto baseline_preset = preset;
    baseline_preset.config.proton_process.enable = false;

    auto enhanced_preset = preset;

    RadiationDoseAlgorithm reference_algorithm;
    RadiationDoseAlgorithm baseline_algorithm;
    RadiationDoseAlgorithm enhanced_algorithm;

    ASSERT_TRUE(reference_algorithm.initialize(reference_preset.config));
    ASSERT_TRUE(baseline_algorithm.initialize(baseline_preset.config));
    ASSERT_TRUE(enhanced_algorithm.initialize(enhanced_preset.config));

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(reference_algorithm.advance(reference_preset.time_step_s));
        ASSERT_TRUE(baseline_algorithm.advance(baseline_preset.time_step_s));
        ASSERT_TRUE(enhanced_algorithm.advance(enhanced_preset.time_step_s));
    }

    const auto& reference_status = reference_algorithm.getStatus();
    const auto& baseline_status = baseline_algorithm.getStatus();
    const auto& enhanced_status = enhanced_algorithm.getStatus();

    ASSERT_EQ(reference_status.layers.size(), baseline_status.layers.size());
    ASSERT_EQ(reference_status.layers.size(), enhanced_status.layers.size());

    const double layer_thickness_m =
        preset.config.thickness_m / static_cast<double>(preset.config.layers);
    const double baseline_rmse =
        normalizedProfileRmse(baseline_status, reference_status, layer_thickness_m);
    const double enhanced_rmse =
        normalizedProfileRmse(enhanced_status, reference_status, layer_thickness_m);

    EXPECT_GT(baseline_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.proton_process_secondary_event_count, 0.0);
    EXPECT_GT(enhanced_status.proton_process_secondary_energy_j_per_m2, 0.0);
    EXPECT_LE(enhanced_rmse, baseline_rmse + 0.01);
    EXPECT_LT(enhanced_rmse, 0.25);
}

TEST(RadiationSmokeTest, HeavyIonProcessEnhancementMeetsDepthDoseCurveAlignmentTarget)
{
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("heavy_ion_storm_dose", preset));

    preset.steps = 7;
    preset.config.enable_monte_carlo_transport = false;
    preset.config.heavy_ion_process.enable = true;
    preset.config.heavy_ion_process.effective_charge_scale = 1.08;
    preset.config.heavy_ion_process.stopping_power_scale = 1.05;
    preset.config.heavy_ion_process.elastic_scattering_scale = 1.12;
    preset.config.heavy_ion_process.fragmentation_secondary_energy_fraction = 0.10;
    preset.config.heavy_ion_process.fragmentation_depth_buildup_fraction = 0.72;
    preset.config.heavy_ion_process.fragmentation_event_rate_scale = 1.10;

    auto reference_preset = preset;
    reference_preset.config.enable_monte_carlo_transport = true;
    reference_preset.config.monte_carlo_histories_per_step = 320;
    reference_preset.config.monte_carlo_max_steps_per_track = 128;
    reference_preset.config.monte_carlo_max_recorded_points = 50000;
    reference_preset.config.monte_carlo_seed = 20260403u;

    auto baseline_preset = preset;
    baseline_preset.config.heavy_ion_process.enable = false;

    auto enhanced_preset = preset;

    RadiationDoseAlgorithm reference_algorithm;
    RadiationDoseAlgorithm baseline_algorithm;
    RadiationDoseAlgorithm enhanced_algorithm;

    ASSERT_TRUE(reference_algorithm.initialize(reference_preset.config));
    ASSERT_TRUE(baseline_algorithm.initialize(baseline_preset.config));
    ASSERT_TRUE(enhanced_algorithm.initialize(enhanced_preset.config));

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        ASSERT_TRUE(reference_algorithm.advance(reference_preset.time_step_s));
        ASSERT_TRUE(baseline_algorithm.advance(baseline_preset.time_step_s));
        ASSERT_TRUE(enhanced_algorithm.advance(enhanced_preset.time_step_s));
    }

    const auto& reference_status = reference_algorithm.getStatus();
    const auto& baseline_status = baseline_algorithm.getStatus();
    const auto& enhanced_status = enhanced_algorithm.getStatus();

    ASSERT_EQ(reference_status.layers.size(), baseline_status.layers.size());
    ASSERT_EQ(reference_status.layers.size(), enhanced_status.layers.size());

    const double layer_thickness_m =
        preset.config.thickness_m / static_cast<double>(preset.config.layers);
    const double baseline_rmse =
        normalizedProfileRmse(baseline_status, reference_status, layer_thickness_m);
    const double enhanced_rmse =
        normalizedProfileRmse(enhanced_status, reference_status, layer_thickness_m);

    EXPECT_GT(baseline_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.deposited_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.heavy_ion_process_fragmentation_event_count, 0.0);
    EXPECT_GT(enhanced_status.heavy_ion_process_secondary_energy_j_per_m2, 0.0);
    EXPECT_GT(enhanced_status.heavy_ion_process_mean_effective_charge_multiplier, 1.0);
    EXPECT_GT(enhanced_status.heavy_ion_process_mean_stopping_multiplier, 1.0);
    EXPECT_LE(enhanced_rmse, baseline_rmse + 0.02);
    EXPECT_LT(enhanced_rmse, 0.30);
}

TEST(RadiationSmokeTest, InitializeFailsForUnknownMaterialWithoutImplicitFallback)
{
    RadiationDoseAlgorithm algorithm;
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    preset.config.material_name = "unknown_material_for_rd003";
    EXPECT_FALSE(algorithm.initialize(preset.config));
}

TEST(RadiationSmokeTest, StoppingDataGovernanceSupportsRollbackAndCache)
{
    RadiationScenarioPreset preset;
    ASSERT_TRUE(
        SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset("geo_electron_belt_dose", preset));

    // Alias token should resolve to canonical PTFE material through unified alias table.
    preset.config.material_name = "teflon";
    preset.config.stopping_data_version = "2026.04-rd999";
    preset.config.stopping_data_allow_rollback = true;

    RadiationDoseAlgorithm rollback_algorithm;
    ASSERT_TRUE(rollback_algorithm.initialize(preset.config));
    const auto& rollback_resolution = rollback_algorithm.getStoppingDataResolution();
    EXPECT_TRUE(rollback_resolution.valid);
    EXPECT_TRUE(rollback_resolution.rollback_applied);
    EXPECT_EQ(rollback_resolution.resolved_version, "2026.04-rd003");
    EXPECT_EQ(rollback_resolution.material_name, "ptfe");
    EXPECT_FALSE(rollback_resolution.dataset_id.empty());

    ASSERT_TRUE(rollback_algorithm.initialize(preset.config));
    const auto& cached_resolution = rollback_algorithm.getStoppingDataResolution();
    EXPECT_TRUE(cached_resolution.valid);
    EXPECT_TRUE(cached_resolution.cache_hit);

    preset.config.stopping_data_allow_rollback = false;
    RadiationDoseAlgorithm strict_algorithm;
    EXPECT_FALSE(strict_algorithm.initialize(preset.config));
}
