#include "SpacecraftInternalChargingAlgorithm.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <regex>
#include <sstream>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{
namespace
{
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kMinMassDensityKgPerM3 = 1.0;
constexpr double kMinTimeScaleS = 1.0e-18;
constexpr const char* kInternalDriveProvenanceSchema =
    "scdat.internal_radiation_drive_provenance.v1";
constexpr const char* kInternalDriveProvenanceContract =
    "internal-radiation-drive-provenance-v1";
constexpr const char* kInternalDriveLayerAlignmentContract =
    "internal-radiation-layer-alignment-v1";
constexpr const char* kInternalResponseCoupledSolveSchema =
    "scdat.internal_response_coupled_solve.v1";
constexpr const char* kInternalResponseCoupledSolveContract =
    "internal-response-coupled-solve-v1";

std::string normalizeToken(std::string text)
{
    std::string token;
    token.reserve(text.size());
    for (const unsigned char ch : text)
    {
        if (std::isalnum(ch))
        {
            token.push_back(static_cast<char>(std::tolower(ch)));
        }
    }
    return token;
}

struct RadiationDepositionLayerRecord
{
    std::size_t index = 0;
    double depth_m = 0.0;
    double deposited_energy_j_per_m3 = 0.0;
    double dose_gy = 0.0;
};

const char* materialStackModelName(InternalMaterialStackModelKind model)
{
    switch (model)
    {
    case InternalMaterialStackModelKind::SpisHarnessBundle:
        return "spis_harness_bundle";
    case InternalMaterialStackModelKind::SpisBacksheetStack:
        return "spis_backsheet_stack";
    case InternalMaterialStackModelKind::SpisLayeredStack:
    default:
        return "spis_layered_stack";
    }
}

const char* geometryModelName(InternalGeometryModelKind model)
{
    switch (model)
    {
    case InternalGeometryModelKind::ShieldedLayerStack1D:
        return "shielded_layer_stack_1d";
    case InternalGeometryModelKind::LayerStack1D:
    default:
        return "layer_stack_1d";
    }
}

const char* primarySourceModelName(InternalPrimarySourceModelKind model)
{
    switch (model)
    {
    case InternalPrimarySourceModelKind::RadiationDriveCoupled:
        return "radiation_drive_coupled";
    case InternalPrimarySourceModelKind::PresetMonoEnergeticFlux:
    default:
        return "preset_monoenergetic_flux";
    }
}

const char* physicsProcessListName(InternalPhysicsProcessListKind model)
{
    switch (model)
    {
    case InternalPhysicsProcessListKind::Geant4EmStandardLike:
        return "geant4_em_standard_like";
    case InternalPhysicsProcessListKind::Geant4ShieldingLike:
    default:
        return "geant4_shielding_like";
    }
}

const char* energyDepositionModelName(InternalEnergyDepositionModelKind model)
{
    switch (model)
    {
    case InternalEnergyDepositionModelKind::ContinuousSlabDeposition:
        return "continuous_slab_deposition";
    case InternalEnergyDepositionModelKind::Geant4StepRecorderLike:
    default:
        return "geant4_step_recorder_like";
    }
}

const char* chargeResponseModelName(InternalChargeResponseModelKind model)
{
    switch (model)
    {
    case InternalChargeResponseModelKind::SpisLayeredDielectric:
        return "spis_layered_dielectric";
    case InternalChargeResponseModelKind::RadiationInducedConductivityRelaxation:
    default:
        return "radiation_induced_conductivity_relaxation";
    }
}

InternalPrimarySourceModelKind effectivePrimarySourceModel(
    const InternalChargingConfiguration& config)
{
    if (config.source_mode == InternalChargingSourceMode::Radiation)
    {
        return InternalPrimarySourceModelKind::RadiationDriveCoupled;
    }
    return config.primary_source_model;
}

std::filesystem::path internalDriveProvenanceSummaryPath(
    const std::filesystem::path& csv_path)
{
    auto stem = csv_path;
    stem.replace_extension();
    return stem.string() + ".drive_provenance_summary.json";
}

std::filesystem::path internalDriveLayerAlignmentPath(
    const std::filesystem::path& csv_path)
{
    auto stem = csv_path;
    stem.replace_extension();
    return stem.string() + ".drive_layer_alignment.csv";
}

std::filesystem::path internalResponseCoupledSolvePath(
    const std::filesystem::path& csv_path)
{
    auto stem = csv_path;
    stem.replace_extension();
    return stem.string() + ".internal_response_coupled_solve.json";
}

std::optional<std::string> extractJsonStringField(const std::string& text,
                                                  const std::string& key)
{
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern) || match.size() < 2)
    {
        return std::nullopt;
    }
    return match[1].str();
}

std::string escapeJsonString(std::string text)
{
    std::string escaped;
    escaped.reserve(text.size() + 8);
    for (const char ch : text)
    {
        switch (ch)
        {
        case '\\':
            escaped += "\\\\";
            break;
        case '\"':
            escaped += "\\\"";
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
            escaped.push_back(ch);
            break;
        }
    }
    return escaped;
}

std::optional<double> extractJsonNumberField(const std::string& text,
                                             const std::string& key)
{
    const std::regex pattern("\"" + key + "\"\\s*:\\s*([-+0-9.eE]+)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern) || match.size() < 2)
    {
        return std::nullopt;
    }

    try
    {
        return std::stod(match[1].str());
    }
    catch (const std::exception&)
    {
        return std::nullopt;
    }
}

std::optional<bool> extractJsonBoolField(const std::string& text,
                                         const std::string& key)
{
    const std::regex pattern("\"" + key + "\"\\s*:\\s*(true|false)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern) || match.size() < 2)
    {
        return std::nullopt;
    }
    return match[1].str() == "true";
}

double sumJsonNumberFieldOccurrences(const std::string& text, const std::string& key)
{
    const std::regex pattern("\"" + key + "\"\\s*:\\s*([-+0-9.eE]+)");
    double total = 0.0;
    for (std::sregex_iterator it(text.begin(), text.end(), pattern), end; it != end; ++it)
    {
        try
        {
            total += std::stod((*it)[1].str());
        }
        catch (const std::exception&)
        {
        }
    }
    return total;
}

double maxJsonNumberFieldOccurrences(const std::string& text, const std::string& key)
{
    const std::regex pattern("\"" + key + "\"\\s*:\\s*([-+0-9.eE]+)");
    double maximum = 0.0;
    bool found = false;
    for (std::sregex_iterator it(text.begin(), text.end(), pattern), end; it != end; ++it)
    {
        try
        {
            maximum = found ? std::max(maximum, std::stod((*it)[1].str()))
                            : std::stod((*it)[1].str());
            found = true;
        }
        catch (const std::exception&)
        {
        }
    }
    return found ? maximum : 0.0;
}

std::size_t countJsonKeyOccurrences(const std::string& text, const std::string& key)
{
    const std::regex pattern("\"" + key + "\"\\s*:");
    return static_cast<std::size_t>(
        std::distance(std::sregex_iterator(text.begin(), text.end(), pattern),
                      std::sregex_iterator()));
}

std::vector<RadiationDepositionLayerRecord> parseRadiationDepositionLayers(
    const std::string& text)
{
    const std::regex pattern(
        "\\{\\s*\"index\"\\s*:\\s*(\\d+)\\s*,\\s*\"depth_m\"\\s*:\\s*([-+0-9.eE]+)"
        "\\s*,\\s*\"deposited_energy_j_per_m3\"\\s*:\\s*([-+0-9.eE]+)\\s*,\\s*\"dose_gy\""
        "\\s*:\\s*([-+0-9.eE]+)\\s*\\}");
    std::vector<RadiationDepositionLayerRecord> records;
    for (std::sregex_iterator it(text.begin(), text.end(), pattern), end; it != end; ++it)
    {
        try
        {
            RadiationDepositionLayerRecord record;
            record.index = static_cast<std::size_t>(std::stoull((*it)[1].str()));
            record.depth_m = std::stod((*it)[2].str());
            record.deposited_energy_j_per_m3 = std::stod((*it)[3].str());
            record.dose_gy = std::stod((*it)[4].str());
            records.push_back(record);
        }
        catch (const std::exception&)
        {
        }
    }
    return records;
}

std::size_t findNearestRadiationLayerByNormalizedDepth(
    const std::vector<RadiationDepositionLayerRecord>& radiation_layers,
    double normalized_internal_depth)
{
    if (radiation_layers.empty())
    {
        return 0;
    }

    const double max_radiation_depth =
        std::max(radiation_layers.back().depth_m, kMinTimeScaleS);
    std::size_t best_index = 0;
    double best_error = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < radiation_layers.size(); ++i)
    {
        const double normalized_radiation_depth =
            std::clamp(radiation_layers[i].depth_m / max_radiation_depth, 0.0, 1.0);
        const double error = std::abs(normalized_internal_depth - normalized_radiation_depth);
        if (error < best_error)
        {
            best_error = error;
            best_index = i;
        }
    }
    return best_index;
}

bool exportInternalDriveProvenanceSummary(
    const std::filesystem::path& csv_path, const InternalChargingStatus& status,
    const InternalChargingRadiationDrive& drive)
{
    if (drive.deposition_history_path.empty() || drive.process_history_path.empty())
    {
        return true;
    }

    const auto deposition_path = std::filesystem::path(drive.deposition_history_path);
    const auto process_path = std::filesystem::path(drive.process_history_path);
    if (!std::filesystem::is_regular_file(deposition_path) ||
        !std::filesystem::is_regular_file(process_path))
    {
        return false;
    }

    std::ifstream deposition_input(deposition_path);
    std::ifstream process_input(process_path);
    if (!deposition_input.is_open() || !process_input.is_open())
    {
        return false;
    }

    std::ostringstream deposition_buffer;
    std::ostringstream process_buffer;
    deposition_buffer << deposition_input.rdbuf();
    process_buffer << process_input.rdbuf();
    const std::string deposition_text = deposition_buffer.str();
    const std::string process_text = process_buffer.str();

    const auto summary_path = internalDriveProvenanceSummaryPath(csv_path);
    std::ofstream output(summary_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        return false;
    }

    const auto particle_species = extractJsonStringField(deposition_text, "particle_species");
    const auto physics_list = extractJsonStringField(deposition_text, "physics_list");
    const auto material = extractJsonStringField(deposition_text, "material");
    const auto deposition_time_s = extractJsonNumberField(deposition_text, "time_s");
    const auto dispatch_mode = extractJsonStringField(process_text, "dispatch_mode");
    const auto secondary_provenance_available =
        extractJsonBoolField(process_text, "secondary_provenance_available");
    const auto track_csv_path = extractJsonStringField(process_text, "track_csv_path");

    output << "{\n";
    output << "  \"schema\": \"" << kInternalDriveProvenanceSchema << "\",\n";
    output << "  \"contract_id\": \"" << kInternalDriveProvenanceContract << "\",\n";
    output << "  \"module\": \"Internal Charging\",\n";
    output << "  \"source_mode\": \"radiation\",\n";
    output << "  \"deposition_record_contract_id\": \"" << escapeJsonString(drive.deposition_record_contract_id)
           << "\",\n";
    output << "  \"process_history_contract_id\": \"" << escapeJsonString(drive.process_history_contract_id) << "\",\n";
    output << "  \"drive_provenance_source\": \"" << escapeJsonString(drive.provenance_source) << "\",\n";
    output << "  \"process_dispatch_mode\": \"" << escapeJsonString(drive.process_dispatch_mode) << "\",\n";
    output << "  \"deposition_history_path\": \"" << escapeJsonString(deposition_path.string()) << "\",\n";
    output << "  \"process_history_path\": \"" << escapeJsonString(process_path.string()) << "\",\n";
    output << "  \"deposition_summary\": {\n";
    output << "    \"particle_species\": \"" << escapeJsonString(particle_species.value_or("")) << "\",\n";
    output << "    \"physics_list\": \"" << escapeJsonString(physics_list.value_or("")) << "\",\n";
    output << "    \"material\": \"" << escapeJsonString(material.value_or("")) << "\",\n";
    output << "    \"time_s\": " << deposition_time_s.value_or(0.0) << ",\n";
    output << "    \"layer_count\": " << countJsonKeyOccurrences(deposition_text, "index") << ",\n";
    output << "    \"max_layer_dose_gy\": "
           << maxJsonNumberFieldOccurrences(deposition_text, "dose_gy") << ",\n";
    output << "    \"sum_layer_deposited_energy_j_per_m3\": "
           << sumJsonNumberFieldOccurrences(deposition_text, "deposited_energy_j_per_m3") << "\n";
    output << "  },\n";
    output << "  \"process_summary\": {\n";
    output << "    \"dispatch_mode\": \"" << escapeJsonString(dispatch_mode.value_or("")) << "\",\n";
    output << "    \"secondary_provenance_available\": "
           << (secondary_provenance_available.value_or(false) ? "true" : "false") << ",\n";
    output << "    \"track_csv_path\": \"" << escapeJsonString(track_csv_path.value_or("")) << "\",\n";
    output << "    \"electron_secondary_event_count\": "
           << extractJsonNumberField(process_text, "electron_secondary_event_count").value_or(0.0)
           << ",\n";
    output << "    \"proton_secondary_event_count\": "
           << extractJsonNumberField(process_text, "proton_secondary_event_count").value_or(0.0)
           << ",\n";
    output << "    \"heavy_ion_fragmentation_event_count\": "
           << extractJsonNumberField(process_text, "heavy_ion_fragmentation_event_count")
                  .value_or(0.0)
           << ",\n";
    output << "    \"track_secondary_event_count\": "
           << extractJsonNumberField(process_text, "track_secondary_event_count").value_or(0.0)
           << ",\n";
    output << "    \"track_secondary_energy_j_per_m2\": "
           << extractJsonNumberField(process_text, "track_secondary_energy_j_per_m2").value_or(0.0)
           << "\n";
    output << "  },\n";
    output << "  \"internal_response_snapshot\": {\n";
    output << "    \"average_dose_gy\": " << status.average_dose_gy << ",\n";
    output << "    \"effective_conductivity_s_per_m\": " << status.effective_conductivity_s_per_m
           << ",\n";
    output << "    \"deposited_charge_rate_c_per_m3_s\": "
           << status.deposited_charge_rate_c_per_m3_s << ",\n";
    output << "    \"deposited_particle_rate_per_m3_s\": "
           << status.deposited_particle_rate_per_m3_s << ",\n";
    output << "    \"max_electric_field_v_per_m\": " << status.max_electric_field_v_per_m << "\n";
    output << "  }\n";
    output << "}\n";
    return static_cast<bool>(output);
}

bool exportInternalDriveLayerAlignmentCsv(
    const std::filesystem::path& csv_path,
    const InternalChargingConfiguration& config,
    const InternalChargingStatus& status,
    const ParticleTransportModel& transport_model,
    const InternalChargingRadiationDrive& drive)
{
    if (drive.deposition_history_path.empty())
    {
        return true;
    }

    const auto deposition_path = std::filesystem::path(drive.deposition_history_path);
    if (!std::filesystem::is_regular_file(deposition_path))
    {
        return false;
    }

    std::ifstream deposition_input(deposition_path);
    if (!deposition_input.is_open())
    {
        return false;
    }

    std::ostringstream deposition_buffer;
    deposition_buffer << deposition_input.rdbuf();
    const auto radiation_layers = parseRadiationDepositionLayers(deposition_buffer.str());
    if (radiation_layers.empty())
    {
        return false;
    }

    const auto alignment_path = internalDriveLayerAlignmentPath(csv_path);
    std::ofstream output(alignment_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        return false;
    }

    const auto& internal_layers = transport_model.getLayerStates();
    const double internal_thickness_m = std::max(config.thickness_m, kMinTimeScaleS);
    const double max_radiation_depth =
        std::max(radiation_layers.back().depth_m, kMinTimeScaleS);

    output << "internal_layer_index,internal_depth_m,internal_dose_gy,internal_charge_density_c_per_m3,"
              "internal_field_v_per_m,radiation_layer_index,radiation_depth_m,radiation_dose_gy,"
              "radiation_deposited_energy_j_per_m3,normalized_depth_offset,dose_ratio_internal_to_radiation\n";

    for (std::size_t i = 0; i < internal_layers.size(); ++i)
    {
        const double normalized_internal_depth =
            std::clamp(internal_layers[i].depth_m / internal_thickness_m, 0.0, 1.0);
        const auto radiation_index =
            findNearestRadiationLayerByNormalizedDepth(radiation_layers, normalized_internal_depth);
        const auto& radiation = radiation_layers[radiation_index];
        const double normalized_radiation_depth =
            std::clamp(radiation.depth_m / max_radiation_depth, 0.0, 1.0);
        const double normalized_depth_offset =
            std::abs(normalized_internal_depth - normalized_radiation_depth);
        const double internal_dose =
            i < internal_layers.size() ? internal_layers[i].dose_gy : 0.0;
        const double dose_ratio =
            internal_dose / std::max(radiation.dose_gy, 1.0e-30);

        output << i << ','
               << internal_layers[i].depth_m << ','
               << internal_dose << ','
               << internal_layers[i].charge_density_c_per_m3 << ','
               << (i < status.electric_field_v_per_m.size() ? status.electric_field_v_per_m[i] : 0.0) << ','
               << radiation.index << ','
               << radiation.depth_m << ','
               << radiation.dose_gy << ','
               << radiation.deposited_energy_j_per_m3 << ','
               << normalized_depth_offset << ','
               << dose_ratio << '\n';
    }

    return static_cast<bool>(output);
}

bool exportInternalResponseCoupledSolveJson(
    const std::filesystem::path& csv_path,
    const InternalChargingConfiguration& config,
    const InternalChargingStatus& status,
    const ParticleTransportModel& transport_model,
    const InternalChargingRadiationDrive& drive)
{
    if (drive.deposition_history_path.empty() || drive.process_history_path.empty())
    {
        return true;
    }

    const auto deposition_path = std::filesystem::path(drive.deposition_history_path);
    const auto process_path = std::filesystem::path(drive.process_history_path);
    if (!std::filesystem::is_regular_file(deposition_path) ||
        !std::filesystem::is_regular_file(process_path))
    {
        return false;
    }

    std::ifstream deposition_input(deposition_path);
    std::ifstream process_input(process_path);
    if (!deposition_input.is_open() || !process_input.is_open())
    {
        return false;
    }

    std::ostringstream deposition_buffer;
    std::ostringstream process_buffer;
    deposition_buffer << deposition_input.rdbuf();
    process_buffer << process_input.rdbuf();
    const std::string deposition_text = deposition_buffer.str();
    const std::string process_text = process_buffer.str();

    const auto radiation_layers = parseRadiationDepositionLayers(deposition_text);
    const auto& internal_layers = transport_model.getLayerStates();
    if (radiation_layers.empty() || internal_layers.empty())
    {
        return false;
    }

    const double internal_thickness_m = std::max(config.thickness_m, kMinTimeScaleS);
    const double max_radiation_depth_m =
        std::max(radiation_layers.back().depth_m, kMinTimeScaleS);

    double max_normalized_depth_offset = 0.0;
    double sum_normalized_depth_offset = 0.0;
    double sum_dose_ratio = 0.0;
    std::size_t paired_layer_count = 0;
    for (std::size_t i = 0; i < internal_layers.size(); ++i)
    {
        const double normalized_internal_depth =
            std::clamp(internal_layers[i].depth_m / internal_thickness_m, 0.0, 1.0);
        const std::size_t radiation_index =
            findNearestRadiationLayerByNormalizedDepth(radiation_layers, normalized_internal_depth);
        const auto& radiation = radiation_layers[radiation_index];
        const double normalized_radiation_depth =
            std::clamp(radiation.depth_m / max_radiation_depth_m, 0.0, 1.0);
        const double depth_offset =
            std::abs(normalized_internal_depth - normalized_radiation_depth);
        max_normalized_depth_offset = std::max(max_normalized_depth_offset, depth_offset);
        sum_normalized_depth_offset += depth_offset;

        const double internal_dose = std::max(0.0, internal_layers[i].dose_gy);
        const double radiation_dose = std::max(1.0e-30, radiation.dose_gy);
        sum_dose_ratio += internal_dose / radiation_dose;
        ++paired_layer_count;
    }

    const double mean_normalized_depth_offset =
        paired_layer_count > 0
            ? (sum_normalized_depth_offset / static_cast<double>(paired_layer_count))
            : 0.0;
    const double mean_dose_ratio_internal_to_radiation =
        paired_layer_count > 0 ? (sum_dose_ratio / static_cast<double>(paired_layer_count)) : 0.0;
    const double charge_state_residual_abs = std::abs(
        status.effective_incident_charge_state_abs - std::max(0.0, drive.incident_charge_state_abs));
    const double dose_ratio_log10_abs =
        std::abs(std::log10(std::max(mean_dose_ratio_internal_to_radiation, 1.0e-30)));
    const double coupled_residual_metric =
        std::max({max_normalized_depth_offset, charge_state_residual_abs, dose_ratio_log10_abs});
    const bool coupled_solve_converged =
        std::isfinite(coupled_residual_metric) && coupled_residual_metric <= 2.0;

    const auto dispatch_mode = extractJsonStringField(process_text, "dispatch_mode");
    const auto secondary_provenance_available =
        extractJsonBoolField(process_text, "secondary_provenance_available");
    const auto track_secondary_event_count =
        extractJsonNumberField(process_text, "track_secondary_event_count");
    const auto track_secondary_energy_j_per_m2 =
        extractJsonNumberField(process_text, "track_secondary_energy_j_per_m2");

    const auto coupled_solve_path = internalResponseCoupledSolvePath(csv_path);
    std::ofstream output(coupled_solve_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        return false;
    }

    output << "{\n";
    output << "  \"schema_version\": \"" << kInternalResponseCoupledSolveSchema << "\",\n";
    output << "  \"contract_id\": \"" << kInternalResponseCoupledSolveContract << "\",\n";
    output << "  \"module\": \"Internal Charging\",\n";
    output << "  \"source_mode\": \"radiation\",\n";
    output << "  \"deposition_record_contract_id\": \""
           << escapeJsonString(drive.deposition_record_contract_id) << "\",\n";
    output << "  \"process_history_contract_id\": \""
           << escapeJsonString(drive.process_history_contract_id) << "\",\n";
    output << "  \"drive_provenance_source\": \""
           << escapeJsonString(drive.provenance_source) << "\",\n";
    output << "  \"process_dispatch_mode\": \""
           << escapeJsonString(dispatch_mode.value_or(drive.process_dispatch_mode)) << "\",\n";
    output << "  \"secondary_provenance_available\": "
           << (secondary_provenance_available.value_or(drive.secondary_provenance_available) ? "true"
                                                                                              : "false")
           << ",\n";
    output << "  \"coupled_solve_converged\": "
           << (coupled_solve_converged ? "true" : "false") << ",\n";
    output << "  \"coupled_residual_metric\": " << coupled_residual_metric << ",\n";
    output << "  \"charge_state_residual_abs\": " << charge_state_residual_abs << ",\n";
    output << "  \"internal_layer_count\": " << internal_layers.size() << ",\n";
    output << "  \"radiation_layer_count\": " << radiation_layers.size() << ",\n";
    output << "  \"alignment_max_normalized_depth_offset\": "
           << max_normalized_depth_offset << ",\n";
    output << "  \"alignment_mean_normalized_depth_offset\": "
           << mean_normalized_depth_offset << ",\n";
    output << "  \"alignment_mean_dose_ratio_internal_to_radiation\": "
           << mean_dose_ratio_internal_to_radiation << ",\n";
    output << "  \"track_secondary_event_count\": "
           << track_secondary_event_count.value_or(0.0) << ",\n";
    output << "  \"track_secondary_energy_j_per_m2\": "
           << track_secondary_energy_j_per_m2.value_or(0.0) << ",\n";
    output << "  \"internal_response\": {\n";
    output << "    \"time_s\": " << status.time_s << ",\n";
    output << "    \"average_dose_gy\": " << status.average_dose_gy << ",\n";
    output << "    \"effective_conductivity_s_per_m\": "
           << status.effective_conductivity_s_per_m << ",\n";
    output << "    \"deposited_charge_rate_c_per_m3_s\": "
           << status.deposited_charge_rate_c_per_m3_s << ",\n";
    output << "    \"deposited_particle_rate_per_m3_s\": "
           << status.deposited_particle_rate_per_m3_s << ",\n";
    output << "    \"max_electric_field_v_per_m\": "
           << status.max_electric_field_v_per_m << "\n";
    output << "  }\n";
    output << "}\n";
    return static_cast<bool>(output);
}
}

bool SpacecraftInternalChargingAlgorithm::initialize(const InternalChargingConfiguration& config)
{
    InternalChargingConfiguration resolved_config = config;
    const auto process_set = normalizeToken(resolved_config.solver_config.physics_process_set);
    if (process_set == "geant4emstandard" || process_set == "standard")
    {
        resolved_config.physics_process_list = InternalPhysicsProcessListKind::Geant4EmStandardLike;
    }
    else if (process_set == "geant4shielding" || process_set == "shielding")
    {
        resolved_config.physics_process_list = InternalPhysicsProcessListKind::Geant4ShieldingLike;
    }

    const auto deposition_scheme =
        normalizeToken(resolved_config.solver_config.deposition_scheme);
    if (deposition_scheme == "continuousslab" || deposition_scheme == "continuous")
    {
        resolved_config.energy_deposition_model =
            InternalEnergyDepositionModelKind::ContinuousSlabDeposition;
    }
    else if (deposition_scheme == "geant4step" || deposition_scheme == "steprecorder")
    {
        resolved_config.energy_deposition_model =
            InternalEnergyDepositionModelKind::Geant4StepRecorderLike;
    }

    config_ = resolved_config;
    material_ = material_database_.findByName(config_.material_name);
    if (material_ == nullptr)
    {
        return false;
    }
    if (!transport_model_.initialize(config_.thickness_m, config_.layers))
    {
        return false;
    }

    status_ = InternalChargingStatus{};
    status_.volume_charge_density_c_per_m3.assign(config_.layers, 0.0);
    status_.electric_field_v_per_m.assign(config_.layers, 0.0);
    initialized_ = true;
    return true;
}

bool SpacecraftInternalChargingAlgorithm::advance(double dt)
{
    if (!initialized_ || material_ == nullptr || dt <= 0.0)
    {
        return false;
    }

    double incident_current_density_a_per_m2 = config_.incident_current_density_a_per_m2;
    double incident_energy_ev = config_.incident_energy_ev;
    double incident_charge_state_abs = config_.incident_charge_state_abs;
    if (config_.source_mode == InternalChargingSourceMode::Radiation)
    {
        if (!radiation_drive_.has_value())
        {
            return false;
        }
        incident_current_density_a_per_m2 = radiation_drive_->incident_current_density_a_per_m2;
        incident_energy_ev = radiation_drive_->incident_energy_ev;
        incident_charge_state_abs = radiation_drive_->incident_charge_state_abs;
    }

    if (incident_current_density_a_per_m2 < 0.0 || incident_energy_ev < 0.0 ||
        incident_charge_state_abs <= 0.0)
    {
        return false;
    }

    const auto& layers_before = transport_model_.getLayerStates();
    std::vector<double> previous_charge_density_c_per_m3;
    std::vector<double> previous_deposited_energy_j_per_m3;
    previous_charge_density_c_per_m3.reserve(layers_before.size());
    previous_deposited_energy_j_per_m3.reserve(layers_before.size());
    for (const auto& layer : layers_before)
    {
        previous_charge_density_c_per_m3.push_back(layer.charge_density_c_per_m3);
        previous_deposited_energy_j_per_m3.push_back(layer.deposited_energy_j_per_m3);
    }

    transport_model_.depositFlux(incident_current_density_a_per_m2, incident_energy_ev,
                                 incident_charge_state_abs, dt);
    updateRadiationInducedConductivity(dt, previous_charge_density_c_per_m3,
                                       previous_deposited_energy_j_per_m3,
                                       incident_charge_state_abs);
    updateFieldEstimate();

    const double max_field =
        status_.electric_field_v_per_m.empty()
            ? 0.0
            : *std::max_element(status_.electric_field_v_per_m.begin(), status_.electric_field_v_per_m.end());
    status_.max_electric_field_v_per_m = max_field;
    const auto assessment =
        breakdown_model_.assess(max_field, material_->getBreakdownFieldVPerM(), status_.time_s);
    status_.discharge_probability = assessment.probability;
    status_.breakdown_active = assessment.triggered;

    if (assessment.triggered)
    {
        auto relieved = discharge_model_.relieve(transport_model_.getLayerStates(), 0.35);
        status_.total_stored_energy_j =
            std::max(0.0, status_.total_stored_energy_j - relieved.released_energy_j_per_m2 * config_.area_m2);
        status_.discharge_events += 1;
        updateFieldEstimate();
    }

    status_.time_s += dt;
    return true;
}

void SpacecraftInternalChargingAlgorithm::reset()
{
    *this = SpacecraftInternalChargingAlgorithm{};
}

bool SpacecraftInternalChargingAlgorithm::exportResults(const std::filesystem::path& csv_path) const
{
    const auto provenance_summary_path = internalDriveProvenanceSummaryPath(csv_path);
    const auto layer_alignment_path = internalDriveLayerAlignmentPath(csv_path);
    const auto coupled_solve_path = internalResponseCoupledSolvePath(csv_path);
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "depth_m";
    const auto& layers = transport_model_.getLayerStates();
    data_set.scalar_series["volume_charge_density_c_per_m3"] = {};
    data_set.scalar_series["electric_field_v_per_m"] = status_.electric_field_v_per_m;
    data_set.scalar_series["deposited_energy_j_per_m3"] = {};
    data_set.scalar_series["dose_gy"] = {};
    data_set.scalar_series["effective_conductivity_s_per_m"] = {};
    data_set.scalar_series["deposited_charge_rate_c_per_m3_s"] = {};
    data_set.scalar_series["deposited_particle_rate_per_m3_s"] = {};

    for (const auto& layer : layers)
    {
        data_set.axis_values.push_back(layer.depth_m);
        data_set.scalar_series["volume_charge_density_c_per_m3"].push_back(layer.charge_density_c_per_m3);
        data_set.scalar_series["deposited_energy_j_per_m3"].push_back(layer.deposited_energy_j_per_m3);
        data_set.scalar_series["dose_gy"].push_back(layer.dose_gy);
        data_set.scalar_series["effective_conductivity_s_per_m"].push_back(
            layer.effective_conductivity_s_per_m);
        data_set.scalar_series["deposited_charge_rate_c_per_m3_s"].push_back(
            layer.deposited_charge_rate_c_per_m3_s);
        data_set.scalar_series["deposited_particle_rate_per_m3_s"].push_back(
            layer.deposited_particle_rate_per_m3_s);
    }
    data_set.metadata["module"] = "Internal Charging";
    data_set.metadata["organization_family"] =
        config_.enable_spis_style_organization ? "spis_geant4_numeric_v1" : "native";
    data_set.metadata["source_mode"] = toString(config_.source_mode);
    data_set.metadata["material_stack_model"] =
        materialStackModelName(config_.material_stack_model);
    data_set.metadata["geometry_model"] = geometryModelName(config_.geometry_model);
    data_set.metadata["primary_source_model"] =
        primarySourceModelName(effectivePrimarySourceModel(config_));
    data_set.metadata["physics_process_list"] =
        physicsProcessListName(config_.physics_process_list);
    data_set.metadata["energy_deposition_model"] =
        energyDepositionModelName(config_.energy_deposition_model);
    data_set.metadata["charge_response_model"] =
        chargeResponseModelName(config_.charge_response_model);
    data_set.metadata["deposition_record_contract_id"] =
        radiation_drive_.has_value()
            ? radiation_drive_->deposition_record_contract_id
            : "aggregate-dose-drive-v1";
    data_set.metadata["process_history_contract_id"] =
        radiation_drive_.has_value()
            ? radiation_drive_->process_history_contract_id
            : "aggregate-process-history-v1";
    data_set.metadata["drive_provenance_source"] =
        radiation_drive_.has_value() ? radiation_drive_->provenance_source : "internal_preset";
    data_set.metadata["drive_process_dispatch_mode"] =
        radiation_drive_.has_value() ? radiation_drive_->process_dispatch_mode : "aggregate";
    data_set.metadata["drive_secondary_provenance_available"] =
        (radiation_drive_.has_value() && radiation_drive_->secondary_provenance_available)
            ? "true"
            : "false";
    data_set.metadata["deposition_history_path"] =
        radiation_drive_.has_value() ? radiation_drive_->deposition_history_path : "";
    data_set.metadata["process_history_path"] =
        radiation_drive_.has_value() ? radiation_drive_->process_history_path : "";
    data_set.metadata["internal_drive_provenance_contract_id"] =
        (radiation_drive_.has_value() &&
         !radiation_drive_->deposition_history_path.empty() &&
         !radiation_drive_->process_history_path.empty())
            ? kInternalDriveProvenanceContract
            : "";
    data_set.metadata["internal_drive_provenance_artifact_path"] =
        (radiation_drive_.has_value() &&
         !radiation_drive_->deposition_history_path.empty() &&
         !radiation_drive_->process_history_path.empty())
            ? provenance_summary_path.string()
            : "";
    data_set.metadata["internal_drive_layer_alignment_contract_id"] =
        (radiation_drive_.has_value() && !radiation_drive_->deposition_history_path.empty())
            ? kInternalDriveLayerAlignmentContract
            : "";
    data_set.metadata["internal_drive_layer_alignment_artifact_path"] =
        (radiation_drive_.has_value() && !radiation_drive_->deposition_history_path.empty())
            ? layer_alignment_path.string()
            : "";
    data_set.metadata["internal_response_coupled_solve_contract_id"] =
        (radiation_drive_.has_value() &&
         !radiation_drive_->deposition_history_path.empty() &&
         !radiation_drive_->process_history_path.empty())
            ? kInternalResponseCoupledSolveContract
            : "";
    data_set.metadata["internal_response_coupled_solve_artifact_path"] =
        (radiation_drive_.has_value() &&
         !radiation_drive_->deposition_history_path.empty() &&
         !radiation_drive_->process_history_path.empty())
            ? coupled_solve_path.string()
            : "";
    data_set.metadata["internal_consumes_deposition_history"] =
        (config_.source_mode == InternalChargingSourceMode::Radiation) ? "true" : "false";
    data_set.metadata["layer_count"] = std::to_string(config_.layers);
    data_set.metadata["thickness_m"] = std::to_string(config_.thickness_m);
    data_set.metadata["max_electric_field_v_per_m"] = std::to_string(status_.max_electric_field_v_per_m);
    data_set.metadata["breakdown_field_v_per_m"] =
        std::to_string(material_ != nullptr ? material_->getBreakdownFieldVPerM() : 0.0);
    data_set.metadata["average_dose_gy"] = std::to_string(status_.average_dose_gy);
    data_set.metadata["effective_conductivity_s_per_m"] =
        std::to_string(status_.effective_conductivity_s_per_m);
    data_set.metadata["deposited_charge_rate_c_per_m3_s"] =
        std::to_string(status_.deposited_charge_rate_c_per_m3_s);
    data_set.metadata["deposited_particle_rate_per_m3_s"] =
        std::to_string(status_.deposited_particle_rate_per_m3_s);
    data_set.metadata["effective_incident_charge_state_abs"] =
        std::to_string(status_.effective_incident_charge_state_abs);
    Coupling::Contracts::appendSolverConfigMetadata(data_set.metadata, config_.solver_config, "solver_");
    data_set.metadata["seed"] = std::to_string(config_.seed);
    data_set.metadata["sampling_policy"] = config_.sampling_policy;
    auto simulation_artifact_path = csv_path;
    simulation_artifact_path.replace_extension(".simulation_artifact.json");
    data_set.metadata["simulation_artifact_contract_id"] = "simulation-artifact-v1";
    data_set.metadata["simulation_artifact_path"] = simulation_artifact_path.filename().string();
    if (!static_cast<bool>(exporter_.exportDataSet(csv_path, data_set)))
    {
        return false;
    }

    if (radiation_drive_.has_value() &&
        !radiation_drive_->deposition_history_path.empty() &&
        !radiation_drive_->process_history_path.empty())
    {
        if (!exportInternalDriveProvenanceSummary(csv_path, status_, *radiation_drive_))
        {
            return false;
        }
    }

    if (radiation_drive_.has_value() && !radiation_drive_->deposition_history_path.empty())
    {
        if (!exportInternalDriveLayerAlignmentCsv(csv_path, config_, status_, transport_model_,
                                                  *radiation_drive_))
        {
            return false;
        }
    }

    if (radiation_drive_.has_value() &&
        !radiation_drive_->deposition_history_path.empty() &&
        !radiation_drive_->process_history_path.empty())
    {
        if (!exportInternalResponseCoupledSolveJson(csv_path, config_, status_, transport_model_,
                                                    *radiation_drive_))
        {
            return false;
        }
    }

    Coupling::Contracts::SimulationArtifact artifact;
    artifact.module = "Internal Charging";
    artifact.case_id = config_.material_name;
    artifact.reference_family = "geant4";
    artifact.seed = config_.seed;
    artifact.sampling_policy = config_.sampling_policy;
    artifact.internal_metrics["time_s"] = status_.time_s;
    artifact.internal_metrics["max_electric_field_v_per_m"] = status_.max_electric_field_v_per_m;
    artifact.internal_metrics["average_dose_gy"] = status_.average_dose_gy;
    artifact.internal_metrics["effective_conductivity_s_per_m"] =
        status_.effective_conductivity_s_per_m;
    artifact.internal_metrics["deposited_charge_rate_c_per_m3_s"] =
        status_.deposited_charge_rate_c_per_m3_s;
    artifact.internal_metrics["deposited_particle_rate_per_m3_s"] =
        status_.deposited_particle_rate_per_m3_s;
    artifact.internal_metrics["discharge_probability"] = status_.discharge_probability;
    artifact.metadata["source_mode"] = toString(config_.source_mode);
    artifact.metadata["physics_process_list"] = physicsProcessListName(config_.physics_process_list);
    artifact.metadata["energy_deposition_model"] =
        energyDepositionModelName(config_.energy_deposition_model);
    artifact.metadata["charge_response_model"] = chargeResponseModelName(config_.charge_response_model);
    std::string artifact_error;
    if (!Coupling::Contracts::writeSimulationArtifactJson(
            simulation_artifact_path, artifact, &artifact_error))
    {
        return false;
    }

    return true;
}

void SpacecraftInternalChargingAlgorithm::updateRadiationInducedConductivity(
    double dt, const std::vector<double>& previous_charge_density_c_per_m3,
    const std::vector<double>& previous_deposited_energy_j_per_m3,
    double incident_charge_state_abs)
{
    auto& layers = transport_model_.getLayerStates();
    if (layers.empty() || dt <= 0.0 || material_ == nullptr)
    {
        status_.average_dose_gy = 0.0;
        status_.effective_conductivity_s_per_m = 0.0;
        status_.deposited_charge_rate_c_per_m3_s = 0.0;
        status_.deposited_particle_rate_per_m3_s = 0.0;
        status_.effective_incident_charge_state_abs = 1.0;
        return;
    }

    const double effective_charge_state_abs = std::max(incident_charge_state_abs, 1.0e-12);

    const double eps_r = std::max(1.0, material_->getPermittivity());
    const double dielectric_permittivity_f_per_m = kEpsilon0 * eps_r;
    const double base_conductivity_s_per_m = std::max(0.0, material_->getConductivity());
    const double min_conductivity_s_per_m =
        std::max(1.0e-20, config_.min_effective_conductivity_s_per_m);
    const double max_conductivity_s_per_m =
        std::max(min_conductivity_s_per_m, config_.max_effective_conductivity_s_per_m);
    const double mass_density_kg_per_m3 =
        std::max(kMinMassDensityKgPerM3, material_->getMassDensityKgPerM3());

    double sum_dose_gy = 0.0;
    double sum_effective_conductivity_s_per_m = 0.0;
    double sum_charge_rate_c_per_m3_s = 0.0;
    double sum_particle_rate_per_m3_s = 0.0;

    for (std::size_t i = 0; i < layers.size(); ++i)
    {
        const double previous_charge =
            i < previous_charge_density_c_per_m3.size() ? previous_charge_density_c_per_m3[i] : 0.0;
        const double previous_energy =
            i < previous_deposited_energy_j_per_m3.size() ? previous_deposited_energy_j_per_m3[i] : 0.0;

        const double deposited_charge_delta_c_per_m3 =
            std::max(0.0, layers[i].charge_density_c_per_m3 - previous_charge);
        const double deposited_energy_delta_j_per_m3 =
            std::max(0.0, layers[i].deposited_energy_j_per_m3 - previous_energy);

        const double deposited_charge_rate_c_per_m3_s = deposited_charge_delta_c_per_m3 / dt;
        const double deposited_particle_rate_per_m3_s =
            deposited_charge_rate_c_per_m3_s / (kElementaryChargeC * effective_charge_state_abs);

        const double dose_gy = std::max(0.0, layers[i].deposited_energy_j_per_m3 / mass_density_kg_per_m3);
        const double dose_rate_gy_per_s =
            std::max(0.0, deposited_energy_delta_j_per_m3 / (mass_density_kg_per_m3 * dt));

        const double radiation_conductivity_s_per_m =
            config_.radiation_conductivity_charge_gain_s_per_m_per_c_per_m3 *
                std::abs(layers[i].charge_density_c_per_m3) +
            config_.radiation_conductivity_dose_gain_s_per_m_per_gy * dose_gy +
            config_.radiation_conductivity_dose_rate_gain_s2_per_m_per_gy * dose_rate_gy_per_s;
        const double effective_conductivity_s_per_m =
            std::clamp(base_conductivity_s_per_m + radiation_conductivity_s_per_m,
                       min_conductivity_s_per_m, max_conductivity_s_per_m);

        const double dielectric_relaxation_time_s =
            dielectric_permittivity_f_per_m /
            std::max(effective_conductivity_s_per_m, min_conductivity_s_per_m);
        const double relaxation_factor =
            std::exp(-dt / std::max(dielectric_relaxation_time_s, kMinTimeScaleS));

        // Radiation-induced conductivity drains trapped charge and moderates field growth.
        layers[i].charge_density_c_per_m3 *= relaxation_factor;

        layers[i].dose_gy = dose_gy;
        layers[i].effective_conductivity_s_per_m = effective_conductivity_s_per_m;
        layers[i].deposited_charge_rate_c_per_m3_s = deposited_charge_rate_c_per_m3_s;
        layers[i].deposited_particle_rate_per_m3_s = deposited_particle_rate_per_m3_s;

        sum_dose_gy += dose_gy;
        sum_effective_conductivity_s_per_m += effective_conductivity_s_per_m;
        sum_charge_rate_c_per_m3_s += deposited_charge_rate_c_per_m3_s;
        sum_particle_rate_per_m3_s += deposited_particle_rate_per_m3_s;
    }

    const double layer_count = static_cast<double>(layers.size());
    status_.average_dose_gy = sum_dose_gy / layer_count;
    status_.effective_conductivity_s_per_m = sum_effective_conductivity_s_per_m / layer_count;
    status_.deposited_charge_rate_c_per_m3_s = sum_charge_rate_c_per_m3_s / layer_count;
    status_.deposited_particle_rate_per_m3_s = sum_particle_rate_per_m3_s / layer_count;
    status_.effective_incident_charge_state_abs = effective_charge_state_abs;
}

void SpacecraftInternalChargingAlgorithm::updateFieldEstimate()
{
    const auto& layers = transport_model_.getLayerStates();
    status_.volume_charge_density_c_per_m3.clear();
    status_.electric_field_v_per_m.assign(layers.size(), 0.0);

    if (layers.empty())
    {
        status_.total_stored_energy_j = 0.0;
        return;
    }

    const double eps_r = std::max(1.0, material_->getPermittivity());
    const double layer_thickness_m = config_.thickness_m / static_cast<double>(layers.size());
    double cumulative_charge = 0.0;
    double stored_energy_j_per_m2 = 0.0;
    for (std::size_t i = 0; i < layers.size(); ++i)
    {
        cumulative_charge += layers[i].charge_density_c_per_m3 * layer_thickness_m;
        status_.electric_field_v_per_m[i] = cumulative_charge / (kEpsilon0 * eps_r);
        status_.volume_charge_density_c_per_m3.push_back(layers[i].charge_density_c_per_m3);
        const double energy_density_j_per_m3 =
            0.5 * kEpsilon0 * eps_r * status_.electric_field_v_per_m[i] * status_.electric_field_v_per_m[i];
        stored_energy_j_per_m2 += energy_density_j_per_m3 * layer_thickness_m;
    }
    status_.total_stored_energy_j = stored_energy_j_per_m2 * config_.area_m2;
}

std::string toString(InternalChargingSourceMode source_mode)
{
    switch (source_mode)
    {
    case InternalChargingSourceMode::Preset:
        return "preset";
    case InternalChargingSourceMode::Radiation:
        return "radiation";
    }
    return "unknown";
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
