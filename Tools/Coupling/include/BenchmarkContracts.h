#pragma once

#include <cstddef>
#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Coupling
{
namespace Contracts
{

struct ReferenceDataset
{
    std::string id;
    std::string source;
    std::string citation;
    std::string artifact_path;
    bool external_reference = true;
};

struct ValidationMetric
{
    std::string id;
    double value = 0.0;
    std::string unit;
};

struct ToleranceProfile
{
    double relative_tolerance = 0.0;
    double absolute_tolerance = 0.0;
    double rmse_tolerance = 0.0;
    double drift_tolerance = 0.0;
};

struct BenchmarkCase
{
    std::string id;
    std::string module;
    std::string reference_family;
    std::vector<ReferenceDataset> reference_datasets;
    std::vector<ValidationMetric> validation_metrics;
    ToleranceProfile tolerance_profile;
};

struct SolverConfig
{
    std::string coupling_mode;
    std::string deposition_scheme;
    std::string collision_set;
    std::string physics_process_set;
    std::string convergence_policy;
    std::size_t max_iterations = 64;
    double residual_tolerance = 1.0e-6;
    double relaxation_factor = 1.0;
};

struct SimulationArtifact
{
    std::string schema_version = "scdat.simulation_artifact.v1";
    std::string contract_id = "simulation-artifact-v1";
    std::string module;
    std::string case_id;
    std::string reference_family = "native";
    unsigned int seed = 0u;
    std::string sampling_policy = "deterministic";
    std::string report_fingerprint;
    std::map<std::string, double> field_metrics;
    std::map<std::string, double> particle_metrics;
    std::map<std::string, double> surface_metrics;
    std::map<std::string, double> radiation_metrics;
    std::map<std::string, double> internal_metrics;
    std::unordered_map<std::string, std::string> metadata;
};

std::string computeArtifactFingerprint(const SimulationArtifact& artifact);

bool writeSimulationArtifactJson(const std::filesystem::path& json_path,
                                 const SimulationArtifact& artifact,
                                 std::string* error_message = nullptr);

bool writeBenchmarkCaseJson(const std::filesystem::path& json_path,
                            const BenchmarkCase& benchmark_case,
                            std::string* error_message = nullptr);

void appendSolverConfigMetadata(std::unordered_map<std::string, std::string>& metadata,
                                const SolverConfig& solver_config,
                                const std::string& key_prefix = "solver_");

} // namespace Contracts
} // namespace Coupling
} // namespace SCDAT
