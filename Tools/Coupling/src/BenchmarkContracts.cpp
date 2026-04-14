#include "../include/BenchmarkContracts.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace SCDAT
{
namespace Coupling
{
namespace Contracts
{
namespace
{

std::string escapeJson(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size());
    for (const unsigned char ch : value)
    {
        switch (ch)
        {
        case '\"':
            escaped += "\\\"";
            break;
        case '\\':
            escaped += "\\\\";
            break;
        case '\b':
            escaped += "\\b";
            break;
        case '\f':
            escaped += "\\f";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            if (ch < 0x20)
            {
                std::ostringstream hex;
                hex << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<int>(ch);
                escaped += hex.str();
            }
            else
            {
                escaped.push_back(static_cast<char>(ch));
            }
            break;
        }
    }
    return escaped;
}

void appendMetricObject(std::ofstream& output, const std::map<std::string, double>& metrics)
{
    output << "{";
    if (!metrics.empty())
    {
        output << "\n";
    }
    std::size_t index = 0;
    for (const auto& [key, value] : metrics)
    {
        output << "    \"" << escapeJson(key) << "\": " << std::setprecision(15) << value;
        output << (++index < metrics.size() ? ",\n" : "\n");
    }
    output << "  }";
}

void appendStringMap(std::ofstream& output,
                     const std::unordered_map<std::string, std::string>& values)
{
    output << "{";
    if (!values.empty())
    {
        output << "\n";
    }

    std::vector<std::pair<std::string, std::string>> ordered(values.begin(), values.end());
    std::sort(ordered.begin(), ordered.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    for (std::size_t index = 0; index < ordered.size(); ++index)
    {
        output << "    \"" << escapeJson(ordered[index].first) << "\": \""
               << escapeJson(ordered[index].second) << "\"";
        output << (index + 1 < ordered.size() ? ",\n" : "\n");
    }
    output << "  }";
}

void appendReferenceDatasets(std::ofstream& output,
                             const std::vector<ReferenceDataset>& reference_datasets)
{
    output << "[";
    if (!reference_datasets.empty())
    {
        output << "\n";
    }
    for (std::size_t index = 0; index < reference_datasets.size(); ++index)
    {
        const auto& dataset = reference_datasets[index];
        output << "    {\n";
        output << "      \"id\": \"" << escapeJson(dataset.id) << "\",\n";
        output << "      \"source\": \"" << escapeJson(dataset.source) << "\",\n";
        output << "      \"citation\": \"" << escapeJson(dataset.citation) << "\",\n";
        output << "      \"artifact_path\": \"" << escapeJson(dataset.artifact_path) << "\",\n";
        output << "      \"external_reference\": "
               << (dataset.external_reference ? "true" : "false") << "\n";
        output << "    }";
        output << (index + 1 < reference_datasets.size() ? ",\n" : "\n");
    }
    output << "  ]";
}

void appendValidationMetrics(std::ofstream& output,
                             const std::vector<ValidationMetric>& validation_metrics)
{
    output << "[";
    if (!validation_metrics.empty())
    {
        output << "\n";
    }
    for (std::size_t index = 0; index < validation_metrics.size(); ++index)
    {
        const auto& metric = validation_metrics[index];
        output << "    {\n";
        output << "      \"id\": \"" << escapeJson(metric.id) << "\",\n";
        output << "      \"value\": " << std::setprecision(15) << metric.value << ",\n";
        output << "      \"unit\": \"" << escapeJson(metric.unit) << "\"\n";
        output << "    }";
        output << (index + 1 < validation_metrics.size() ? ",\n" : "\n");
    }
    output << "  ]";
}

std::uint64_t fnv1a64Update(std::uint64_t state, const std::string& text)
{
    constexpr std::uint64_t kOffsetBasis = 1469598103934665603ULL;
    constexpr std::uint64_t kPrime = 1099511628211ULL;
    state = state == 0 ? kOffsetBasis : state;
    for (const unsigned char ch : text)
    {
        state ^= static_cast<std::uint64_t>(ch);
        state *= kPrime;
    }
    return state;
}

} // namespace

std::string computeArtifactFingerprint(const SimulationArtifact& artifact)
{
    std::uint64_t hash = 0;
    hash = fnv1a64Update(hash, artifact.schema_version);
    hash = fnv1a64Update(hash, artifact.contract_id);
    hash = fnv1a64Update(hash, artifact.module);
    hash = fnv1a64Update(hash, artifact.case_id);
    hash = fnv1a64Update(hash, artifact.reference_family);
    hash = fnv1a64Update(hash, std::to_string(artifact.seed));
    hash = fnv1a64Update(hash, artifact.sampling_policy);
    for (const auto& [key, value] : artifact.field_metrics)
    {
        hash = fnv1a64Update(hash, key + "=" + std::to_string(value));
    }
    for (const auto& [key, value] : artifact.particle_metrics)
    {
        hash = fnv1a64Update(hash, key + "=" + std::to_string(value));
    }
    for (const auto& [key, value] : artifact.surface_metrics)
    {
        hash = fnv1a64Update(hash, key + "=" + std::to_string(value));
    }
    for (const auto& [key, value] : artifact.radiation_metrics)
    {
        hash = fnv1a64Update(hash, key + "=" + std::to_string(value));
    }
    for (const auto& [key, value] : artifact.internal_metrics)
    {
        hash = fnv1a64Update(hash, key + "=" + std::to_string(value));
    }

    std::ostringstream stream;
    stream << std::hex << std::setw(16) << std::setfill('0') << hash;
    return stream.str();
}

bool writeSimulationArtifactJson(const std::filesystem::path& json_path,
                                 const SimulationArtifact& artifact,
                                 std::string* error_message)
{
    std::ofstream output(json_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        if (error_message != nullptr)
        {
            *error_message = "Unable to open simulation artifact file: " + json_path.string();
        }
        return false;
    }

    const std::string fingerprint =
        artifact.report_fingerprint.empty() ? computeArtifactFingerprint(artifact)
                                            : artifact.report_fingerprint;

    output << "{\n";
    output << "  \"schema_version\": \"" << escapeJson(artifact.schema_version) << "\",\n";
    output << "  \"contract_id\": \"" << escapeJson(artifact.contract_id) << "\",\n";
    output << "  \"module\": \"" << escapeJson(artifact.module) << "\",\n";
    output << "  \"case_id\": \"" << escapeJson(artifact.case_id) << "\",\n";
    output << "  \"reference_family\": \"" << escapeJson(artifact.reference_family) << "\",\n";
    output << "  \"seed\": " << artifact.seed << ",\n";
    output << "  \"sampling_policy\": \"" << escapeJson(artifact.sampling_policy) << "\",\n";
    output << "  \"report_fingerprint\": \"" << escapeJson(fingerprint) << "\",\n";
    output << "  \"field\": ";
    appendMetricObject(output, artifact.field_metrics);
    output << ",\n";
    output << "  \"particle\": ";
    appendMetricObject(output, artifact.particle_metrics);
    output << ",\n";
    output << "  \"surface\": ";
    appendMetricObject(output, artifact.surface_metrics);
    output << ",\n";
    output << "  \"radiation\": ";
    appendMetricObject(output, artifact.radiation_metrics);
    output << ",\n";
    output << "  \"internal\": ";
    appendMetricObject(output, artifact.internal_metrics);
    output << ",\n";
    output << "  \"metadata\": ";
    appendStringMap(output, artifact.metadata);
    output << "\n";
    output << "}\n";

    if (!static_cast<bool>(output))
    {
        if (error_message != nullptr)
        {
            *error_message = "Failed to flush simulation artifact file: " + json_path.string();
        }
        return false;
    }
    return true;
}

bool writeBenchmarkCaseJson(const std::filesystem::path& json_path,
                            const BenchmarkCase& benchmark_case,
                            std::string* error_message)
{
    std::ofstream output(json_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        if (error_message != nullptr)
        {
            *error_message = "Unable to open benchmark case file: " + json_path.string();
        }
        return false;
    }

    output << "{\n";
    output << "  \"schema_version\": \"scdat.benchmark_case.v1\",\n";
    output << "  \"contract_id\": \"benchmark-case-v1\",\n";
    output << "  \"id\": \"" << escapeJson(benchmark_case.id) << "\",\n";
    output << "  \"module\": \"" << escapeJson(benchmark_case.module) << "\",\n";
    output << "  \"reference_family\": \"" << escapeJson(benchmark_case.reference_family)
           << "\",\n";
    output << "  \"reference_datasets\": ";
    appendReferenceDatasets(output, benchmark_case.reference_datasets);
    output << ",\n";
    output << "  \"validation_metrics\": ";
    appendValidationMetrics(output, benchmark_case.validation_metrics);
    output << ",\n";
    output << "  \"tolerance_profile\": {\n";
    output << "    \"relative_tolerance\": "
           << std::setprecision(15) << benchmark_case.tolerance_profile.relative_tolerance
           << ",\n";
    output << "    \"absolute_tolerance\": "
           << std::setprecision(15) << benchmark_case.tolerance_profile.absolute_tolerance
           << ",\n";
    output << "    \"rmse_tolerance\": "
           << std::setprecision(15) << benchmark_case.tolerance_profile.rmse_tolerance
           << ",\n";
    output << "    \"drift_tolerance\": "
           << std::setprecision(15) << benchmark_case.tolerance_profile.drift_tolerance
           << "\n";
    output << "  }\n";
    output << "}\n";

    if (!static_cast<bool>(output))
    {
        if (error_message != nullptr)
        {
            *error_message = "Failed to flush benchmark case file: " + json_path.string();
        }
        return false;
    }

    return true;
}

void appendSolverConfigMetadata(std::unordered_map<std::string, std::string>& metadata,
                                const SolverConfig& solver_config,
                                const std::string& key_prefix)
{
    metadata[key_prefix + "coupling_mode"] = solver_config.coupling_mode;
    metadata[key_prefix + "deposition_scheme"] = solver_config.deposition_scheme;
    metadata[key_prefix + "collision_set"] = solver_config.collision_set;
    metadata[key_prefix + "physics_process_set"] = solver_config.physics_process_set;
    metadata[key_prefix + "convergence_policy"] = solver_config.convergence_policy;
    metadata[key_prefix + "max_iterations"] = std::to_string(solver_config.max_iterations);
    metadata[key_prefix + "residual_tolerance"] = std::to_string(solver_config.residual_tolerance);
    metadata[key_prefix + "relaxation_factor"] = std::to_string(solver_config.relaxation_factor);
}

} // namespace Contracts
} // namespace Coupling
} // namespace SCDAT
