#include "../include/BenchmarkContracts.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace
{

std::string readText(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

} // namespace

TEST(BenchmarkContractsTest, AppendSolverConfigMetadataWritesExpectedKeys)
{
    SCDAT::Coupling::Contracts::SolverConfig solver_config;
    solver_config.coupling_mode = "field_particle_implicit";
    solver_config.deposition_scheme = "pic_window_cic";
    solver_config.collision_set = "surface_pic_v1";
    solver_config.physics_process_set = "surface_process_v1";
    solver_config.convergence_policy = "residual_norm_guarded";
    solver_config.max_iterations = 96;
    solver_config.residual_tolerance = 1.0e-7;
    solver_config.relaxation_factor = 0.72;

    std::unordered_map<std::string, std::string> metadata;
    SCDAT::Coupling::Contracts::appendSolverConfigMetadata(metadata, solver_config, "solver_");

    EXPECT_EQ(metadata["solver_coupling_mode"], "field_particle_implicit");
    EXPECT_EQ(metadata["solver_deposition_scheme"], "pic_window_cic");
    EXPECT_EQ(metadata["solver_collision_set"], "surface_pic_v1");
    EXPECT_EQ(metadata["solver_physics_process_set"], "surface_process_v1");
    EXPECT_EQ(metadata["solver_convergence_policy"], "residual_norm_guarded");
    EXPECT_EQ(metadata["solver_max_iterations"], "96");
}

TEST(BenchmarkContractsTest, BenchmarkCaseCarriesReferenceAndTolerance)
{
    SCDAT::Coupling::Contracts::BenchmarkCase benchmark_case;
    benchmark_case.id = "surface_plate_sheath_surface";
    benchmark_case.module = "surface";
    benchmark_case.reference_family = "surface";
    benchmark_case.reference_datasets.push_back(
        SCDAT::Coupling::Contracts::ReferenceDataset{
            "surface_case_001",
            "surface",
            "surface benchmark dataset",
            "ref/surface/cases/sheath_case_001.csv",
            true,
        });
    benchmark_case.validation_metrics.push_back(
        SCDAT::Coupling::Contracts::ValidationMetric{
            "surface_potential_v",
            -12400.0,
            "V",
        });
    benchmark_case.tolerance_profile.relative_tolerance = 0.10;
    benchmark_case.tolerance_profile.absolute_tolerance = 1.0e-6;

    ASSERT_EQ(benchmark_case.reference_datasets.size(), 1u);
    ASSERT_EQ(benchmark_case.validation_metrics.size(), 1u);
    EXPECT_EQ(benchmark_case.reference_datasets.front().id, "surface_case_001");
    EXPECT_EQ(benchmark_case.validation_metrics.front().id, "surface_potential_v");
    EXPECT_NEAR(benchmark_case.tolerance_profile.relative_tolerance, 0.10, 1.0e-12);
}

TEST(BenchmarkContractsTest, SimulationArtifactJsonContainsAllCoreGroups)
{
    SCDAT::Coupling::Contracts::SimulationArtifact artifact;
    artifact.module = "Surface Charging";
    artifact.case_id = "geo_ecss_kapton_2000s";
    artifact.reference_family = "surface";
    artifact.seed = 20260408u;
    artifact.sampling_policy = "deterministic";
    artifact.field_metrics["normal_electric_field_v_per_m"] = 1.2e5;
    artifact.particle_metrics["live_pic_collision_event_total"] = 42.0;
    artifact.surface_metrics["final_surface_potential_v"] = -12450.0;
    artifact.radiation_metrics["deposited_energy_j_per_m2"] = 0.0;
    artifact.internal_metrics["max_electric_field_v_per_m"] = 0.0;
    artifact.metadata["pipeline_contract_id"] = "surface-pic-runtime-v1";

    const auto json_path =
        std::filesystem::temp_directory_path() / "benchmark_contracts_artifact_smoke.json";
    std::string error_message;
    ASSERT_TRUE(SCDAT::Coupling::Contracts::writeSimulationArtifactJson(
        json_path, artifact, &error_message));
    ASSERT_TRUE(std::filesystem::exists(json_path));

    const std::string payload = readText(json_path);
    EXPECT_NE(payload.find("\"schema_version\": \"scdat.simulation_artifact.v1\""),
              std::string::npos);
    EXPECT_NE(payload.find("\"contract_id\": \"simulation-artifact-v1\""), std::string::npos);
    EXPECT_NE(payload.find("\"field\": {"), std::string::npos);
    EXPECT_NE(payload.find("\"particle\": {"), std::string::npos);
    EXPECT_NE(payload.find("\"surface\": {"), std::string::npos);
    EXPECT_NE(payload.find("\"radiation\": {"), std::string::npos);
    EXPECT_NE(payload.find("\"internal\": {"), std::string::npos);
    EXPECT_NE(payload.find("\"report_fingerprint\": "), std::string::npos);
}

TEST(BenchmarkContractsTest, BenchmarkCaseJsonContainsDatasetsMetricsAndTolerance)
{
    SCDAT::Coupling::Contracts::BenchmarkCase benchmark_case;
    benchmark_case.id = "geo_ecss_kapton_2000s";
    benchmark_case.module = "surface";
    benchmark_case.reference_family = "surface_pic_reference_family";
    benchmark_case.reference_datasets.push_back(
        SCDAT::Coupling::Contracts::ReferenceDataset{
            "patch_curve",
            "surface",
            "public reference curve",
            "ref/benchmark/patch_curve.csv",
            true,
        });
    benchmark_case.validation_metrics.push_back(
        SCDAT::Coupling::Contracts::ValidationMetric{
            "benchmark_patch_rmse_v",
            4.2,
            "V",
        });
    benchmark_case.tolerance_profile.relative_tolerance = 0.10;
    benchmark_case.tolerance_profile.absolute_tolerance = 10.0;
    benchmark_case.tolerance_profile.rmse_tolerance = 10.0;
    benchmark_case.tolerance_profile.drift_tolerance = 10.0;

    const auto json_path =
        std::filesystem::temp_directory_path() / "benchmark_contracts_case_smoke.json";
    std::string error_message;
    ASSERT_TRUE(SCDAT::Coupling::Contracts::writeBenchmarkCaseJson(
        json_path, benchmark_case, &error_message));
    ASSERT_TRUE(std::filesystem::exists(json_path));

    const std::string payload = readText(json_path);
    EXPECT_NE(payload.find("\"schema_version\": \"scdat.benchmark_case.v1\""),
              std::string::npos);
    EXPECT_NE(payload.find("\"contract_id\": \"benchmark-case-v1\""), std::string::npos);
    EXPECT_NE(payload.find("\"reference_datasets\": ["), std::string::npos);
    EXPECT_NE(payload.find("\"validation_metrics\": ["), std::string::npos);
    EXPECT_NE(payload.find("\"benchmark_patch_rmse_v\""), std::string::npos);
    EXPECT_NE(payload.find("\"rmse_tolerance\": 10"), std::string::npos);
}


