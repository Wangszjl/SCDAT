#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct LegacyBenchmarkPaths
{
    std::filesystem::path input_path;
    std::filesystem::path environment_path;
    std::filesystem::path structure_material_path;
    std::filesystem::path dielectric_patch_material_path;
    std::filesystem::path metal_patch_material_path;
    std::filesystem::path patch_reference_curve_path;
    std::filesystem::path body_reference_curve_path;
    std::filesystem::path matlab_generator_path;
};

struct LegacyBenchmarkInputConfig
{
    std::vector<double> raw_values;
};

struct LegacyEnvironmentRecord
{
    int case_id = 0;
    std::string name;
    std::vector<double> thermal_values;
    std::vector<double> electron_energy_ev;
    std::vector<double> electron_flux_per_m2_s_sr_ev;
    std::vector<double> ion_energy_ev;
    std::vector<double> ion_flux_per_m2_s_sr_ev;
};

struct LegacyMaterialRecord
{
    int material_id = 0;
    std::string name;
    std::vector<double> numeric_values;
};

struct LegacyBenchmarkCurveSample
{
    std::size_t cycle_index = 0;
    double time_s = 0.0;
    double potential_v = 0.0;
    double jnet_a_per_m2 = 0.0;
    double je_a_per_m2 = 0.0;
    double jse_a_per_m2 = 0.0;
    double jb_a_per_m2 = 0.0;
    double ji_a_per_m2 = 0.0;
    double jsi_a_per_m2 = 0.0;
    double jph_a_per_m2 = 0.0;
    double jcond_a_per_m2 = 0.0;
};

struct LegacyBenchmarkCaseDefinition
{
    LegacyBenchmarkPaths paths;
    LegacyBenchmarkInputConfig input;
    std::vector<LegacyEnvironmentRecord> environments;
    std::vector<LegacyMaterialRecord> structure_materials;
    std::vector<LegacyMaterialRecord> dielectric_patch_materials;
    std::vector<LegacyMaterialRecord> metal_patch_materials;
    std::vector<LegacyBenchmarkCurveSample> patch_reference_curve;
    std::vector<LegacyBenchmarkCurveSample> body_reference_curve;
};

struct LegacyBenchmarkMetrics
{
    bool valid = false;
    std::size_t compared_sample_count = 0;
    double rmse_v = 0.0;
    double terminal_potential_delta_v = 0.0;
    double terminal_time_delta_s = 0.0;
    double time_to_equilibrium_delta_s = 0.0;
};

struct LegacyBenchmarkBodyTailMetrics
{
    bool valid = false;
    std::size_t sample_count = 0;
    double threshold_v = 0.0;
    double je_rmse_a_per_m2 = 0.0;
    double jnet_rmse_a_per_m2 = 0.0;
};

struct LegacyBenchmarkConsistencyDiagnostics
{
    bool valid = false;
    bool input_reconstruction_supported = false;
    bool input_reference_consistent = true;
    std::string authority = "ReferenceCurve";
    std::string status = "NotApplicable";
    double reconstructed_body_initial_je_a_per_m2 = 0.0;
    double reference_body_initial_je_a_per_m2 = 0.0;
    double body_initial_je_delta_a_per_m2 = 0.0;
    double body_initial_je_ratio = 0.0;
};

struct LegacyBenchmarkAcceptanceCriteria
{
    std::string id = "none";
    std::string note;
    std::string focus_segment = "global";
    double patch_rmse_v_max = 0.0;
    double body_rmse_v_max = 0.0;
    double body_je_rmse_a_per_m2_max = 0.0;
    double body_jnet_rmse_a_per_m2_max = 0.0;
    double negative_tail_body_je_rmse_a_per_m2_max = 0.0;
    double negative_tail_body_jnet_rmse_a_per_m2_max = 0.0;
};

struct LegacyBenchmarkAcceptanceGateResult
{
    bool applicable = false;
    bool pass = true;
    std::size_t checks_total = 0;
    std::size_t checks_failed = 0;
    std::string status = "NOT_APPLICABLE";
    std::string failed_metric_keys = "none";
    std::string failure_details = "none";
};

double legacyBenchmarkExecutionModeId(LegacyBenchmarkExecutionMode mode);
std::string legacyBenchmarkExecutionModeName(LegacyBenchmarkExecutionMode mode);
std::string legacyBenchmarkBaselineFamilyName(SurfaceBenchmarkSource source);
std::string legacyBenchmarkBaselineOriginName(SurfaceBenchmarkSource source);
LegacyBenchmarkAcceptanceCriteria legacyBenchmarkAcceptanceCriteria(
    SurfaceBenchmarkSource source);

bool resolveLegacyBenchmarkPaths(SurfaceBenchmarkSource source, LegacyBenchmarkPaths& paths);
LegacyBenchmarkInputConfig parseLegacyBenchmarkInputFile(const std::filesystem::path& path);
std::vector<LegacyEnvironmentRecord>
parseLegacyEnvironmentFile(const std::filesystem::path& path);
std::vector<LegacyMaterialRecord>
parseLegacyMaterialTable(const std::filesystem::path& path, std::size_t trailing_numeric_count);
std::vector<LegacyBenchmarkCurveSample>
loadLegacyBenchmarkCurve(const std::filesystem::path& path, bool has_conduction_column);
double interpolateLegacyBenchmarkPotential(const std::vector<LegacyBenchmarkCurveSample>& samples,
                                          double time_s);
LegacyBenchmarkCaseDefinition loadLegacyBenchmarkCaseDefinition(SurfaceBenchmarkSource source);
SurfaceChargingConfig applyLegacyBenchmarkExecutionConfig(
    const SurfaceChargingConfig& fallback_config,
    LegacyBenchmarkCaseDefinition* resolved_definition = nullptr);
double firstLegacyEquilibriumTimeS(const std::vector<LegacyBenchmarkCurveSample>& curve,
                                   double equilibrium_tolerance_v = 1.0e-2);
LegacyBenchmarkMetrics computeLegacyBenchmarkMetrics(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    double equilibrium_tolerance_v = 1.0e-2);
LegacyBenchmarkBodyTailMetrics computeLegacyBenchmarkBodyTailMetrics(
    const std::vector<LegacyBenchmarkCurveSample>& actual_body,
    const std::vector<LegacyBenchmarkCurveSample>& reference_body);
LegacyBenchmarkAcceptanceGateResult evaluateLegacyBenchmarkAcceptanceGate(
    const LegacyBenchmarkAcceptanceCriteria& criteria,
    const LegacyBenchmarkMetrics& patch_metrics,
    const LegacyBenchmarkMetrics& body_metrics,
    double body_je_rmse_a_per_m2,
    double body_jnet_rmse_a_per_m2,
    const LegacyBenchmarkBodyTailMetrics& body_tail_metrics);
LegacyBenchmarkConsistencyDiagnostics analyzeLegacyBenchmarkConsistency(
    SurfaceBenchmarkSource source, const LegacyBenchmarkCaseDefinition& definition);
bool writeLegacyBenchmarkReport(const std::filesystem::path& report_path,
                                SurfaceRuntimeRoute runtime_route,
                                SurfaceBenchmarkSource source,
                                LegacyBenchmarkExecutionMode execution_mode,
                                const LegacyBenchmarkCaseDefinition& definition,
                                const std::vector<LegacyBenchmarkCurveSample>& actual_patch,
                                const std::vector<LegacyBenchmarkCurveSample>& actual_body,
                                const LegacyBenchmarkMetrics& patch_metrics,
                                const LegacyBenchmarkMetrics& body_metrics,
                                const LegacyBenchmarkConsistencyDiagnostics& consistency);
bool writeLegacyBenchmarkComparisonCsv(
    const std::filesystem::path& csv_path,
    const std::vector<LegacyBenchmarkCurveSample>& actual_patch,
    const std::vector<LegacyBenchmarkCurveSample>& reference_patch,
    const std::vector<LegacyBenchmarkCurveSample>& actual_body,
    const std::vector<LegacyBenchmarkCurveSample>& reference_body);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
