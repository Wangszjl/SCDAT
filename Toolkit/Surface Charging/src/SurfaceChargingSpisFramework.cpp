#include "DensePlasmaSurfaceCharging.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kBoltzmannEvPerK = 8.617333262145e-5;
constexpr double kRichardsonConstant = 1.20173e6;
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kPi = 3.14159265358979323846;

using Coupling::SurfaceCircuitAdvanceResult;
using Coupling::SurfaceCircuitBranch;
using Coupling::SurfaceCircuitCoupling;
using Coupling::SurfaceCircuitLinearization;
using Coupling::SurfaceCircuitNode;

bool isReferenceRegime(const SurfaceChargingConfig& config)
{
    return config.use_reference_current_balance &&
           (config.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
            config.regime == SurfaceChargingRegime::GeoKineticPicLike);
}

bool isBodyNodeName(const std::string& name)
{
    return name.rfind("body:", 0) == 0 || name == "body";
}

bool isPatchNodeName(const std::string& name)
{
    return name.rfind("patch:", 0) == 0 || (!name.empty() && !isBodyNodeName(name));
}

std::string bodyNodeName(const std::string& body_id)
{
    return "body:" + body_id;
}

std::string patchNodeName(const std::string& patch_id)
{
    return "patch:" + patch_id;
}

std::string benchmarkSourceName(SurfaceBenchmarkSource source)
{
    switch (source)
    {
    case SurfaceBenchmarkSource::CGeo:
        return "C-GEO";
    case SurfaceBenchmarkSource::CLeoRam:
        return "C-LEO-RAM";
    case SurfaceBenchmarkSource::CLeoWake:
        return "C-LEO-WAKE";
    case SurfaceBenchmarkSource::MatlabGeo:
        return "MATLAB-GEO";
    case SurfaceBenchmarkSource::MatlabLeo:
        return "MATLAB-LEO";
    default:
        return "None";
    }
}

std::string runtimeRouteName(SurfaceRuntimeRoute route)
{
    return route == SurfaceRuntimeRoute::LegacyBenchmark ? "LegacyBenchmark"
                                                         : "SCDATUnified";
}

bool hasStructuredTopologyInput(const SurfaceChargingConfig& config)
{
    return !config.bodies.empty() || !config.patches.empty() || !config.interfaces.empty();
}

bool failValidation(const std::string& message, std::string* error_message)
{
    if (error_message != nullptr)
    {
        *error_message = message;
    }
    return false;
}

bool validateStructuredTopologyConfig(const SurfaceChargingConfig& config,
                                     std::string* error_message)
{
    if (!hasStructuredTopologyInput(config))
    {
        return true;
    }

    std::unordered_map<std::string, bool> body_ids;
    std::unordered_map<std::string, bool> patch_ids;
    std::unordered_map<std::string, bool> all_node_ids;
    std::unordered_map<std::string, bool> interface_ids;

    for (const auto& body : config.bodies)
    {
        if (body.id.empty())
        {
            return failValidation("Structured topology body id must not be empty.", error_message);
        }
        if (body_ids.contains(body.id) || all_node_ids.contains(body.id))
        {
            return failValidation("Duplicate structured topology body id: " + body.id,
                                  error_message);
        }
        if (body.area_m2 < 0.0)
        {
            return failValidation("Structured topology body area must be non-negative for body: " +
                                      body.id,
                                  error_message);
        }
        if (body.capacitance_f < 0.0)
        {
            return failValidation(
                "Structured topology body capacitance must be non-negative for body: " + body.id,
                error_message);
        }
        body_ids[body.id] = true;
        all_node_ids[body.id] = true;
    }

    for (const auto& patch : config.patches)
    {
        if (patch.id.empty())
        {
            return failValidation("Structured topology patch id must not be empty.", error_message);
        }
        if (patch_ids.contains(patch.id) || all_node_ids.contains(patch.id))
        {
            return failValidation("Duplicate structured topology patch id: " + patch.id,
                                  error_message);
        }
        if (patch.body_id.empty())
        {
            return failValidation("Structured topology patch must reference a body: " + patch.id,
                                  error_message);
        }
        if (!body_ids.contains(patch.body_id))
        {
            return failValidation("Structured topology patch references unknown body '" +
                                      patch.body_id + "': " + patch.id,
                                  error_message);
        }
        if (patch.area_m2 <= 0.0)
        {
            return failValidation("Structured topology patch area must be positive for patch: " +
                                      patch.id,
                                  error_message);
        }
        if (patch.capacitance_f < 0.0)
        {
            return failValidation(
                "Structured topology patch capacitance must be non-negative for patch: " +
                    patch.id,
                error_message);
        }
        patch_ids[patch.id] = true;
        all_node_ids[patch.id] = true;
    }

    if (config.bodies.empty() && !config.patches.empty())
    {
        return failValidation("Structured topology patches require at least one body.",
                              error_message);
    }

    std::unordered_set<std::string> boundary_group_ids;
    for (const auto& group : config.body_boundary_groups)
    {
        if (group.id.empty())
        {
            return failValidation("Structured topology body boundary group id must not be empty.",
                                  error_message);
        }
        if (!boundary_group_ids.insert(group.id).second)
        {
            return failValidation("Duplicate structured topology boundary group id: " + group.id,
                                  error_message);
        }
        if (group.body_id.empty() || !body_ids.contains(group.body_id))
        {
            return failValidation("Structured topology body boundary group references unknown body '" +
                                      group.body_id + "': " + group.id,
                                  error_message);
        }
    }

    for (const auto& group : config.patch_boundary_groups)
    {
        if (group.id.empty())
        {
            return failValidation("Structured topology patch boundary group id must not be empty.",
                                  error_message);
        }
        if (!boundary_group_ids.insert(group.id).second)
        {
            return failValidation("Duplicate structured topology boundary group id: " + group.id,
                                  error_message);
        }
        if (group.patch_id.empty() || !patch_ids.contains(group.patch_id))
        {
            return failValidation("Structured topology patch boundary group references unknown patch '" +
                                      group.patch_id + "': " + group.id,
                                  error_message);
        }
    }

    const auto mapping_node_is_known = [&](const std::string& node_id) {
        if (all_node_ids.contains(node_id))
        {
            return true;
        }
        if (node_id.rfind("body:", 0) == 0)
        {
            return all_node_ids.contains(node_id.substr(std::string("body:").size()));
        }
        if (node_id.rfind("patch:", 0) == 0)
        {
            return all_node_ids.contains(node_id.substr(std::string("patch:").size()));
        }
        return false;
    };

    for (const auto& mapping : config.boundary_mappings)
    {
        if (mapping.node_id.empty())
        {
            return failValidation("Structured topology boundary mapping node id must not be empty.",
                                  error_message);
        }
        if (!mapping_node_is_known(mapping.node_id))
        {
            return failValidation("Structured topology boundary mapping references unknown node '" +
                                      mapping.node_id + "'.",
                                  error_message);
        }
        if (mapping.boundary_group_id.empty() || !boundary_group_ids.contains(mapping.boundary_group_id))
        {
            return failValidation(
                "Structured topology boundary mapping references unknown boundary group '" +
                    mapping.boundary_group_id + "'.",
                error_message);
        }
    }

    for (const auto& interface_config : config.interfaces)
    {
        if (!interface_config.id.empty())
        {
            if (interface_ids.contains(interface_config.id))
            {
                return failValidation("Duplicate structured topology interface id: " +
                                          interface_config.id,
                                      error_message);
            }
            interface_ids[interface_config.id] = true;
        }
        if (interface_config.from_id.empty() || interface_config.to_id.empty())
        {
            return failValidation(
                "Structured topology interface endpoints must not be empty.",
                error_message);
        }
        if (!all_node_ids.contains(interface_config.from_id))
        {
            return failValidation("Structured topology interface references unknown from_id '" +
                                      interface_config.from_id + "'.",
                                  error_message);
        }
        if (!all_node_ids.contains(interface_config.to_id))
        {
            return failValidation("Structured topology interface references unknown to_id '" +
                                      interface_config.to_id + "'.",
                                  error_message);
        }
        if (interface_config.from_id == interface_config.to_id)
        {
            return failValidation(
                "Structured topology interface must connect distinct endpoints: " +
                    interface_config.from_id,
                error_message);
        }
        if (interface_config.conductance_s < 0.0 || interface_config.resistance_ohm < 0.0)
        {
            return failValidation(
                "Structured topology interface conductance/resistance must be non-negative.",
                error_message);
        }
    }

    return true;
}

std::string logicalNodeIdFromRuntimeName(const std::string& node_name)
{
    if (node_name.rfind("body:", 0) == 0)
    {
        return node_name.substr(std::string("body:").size());
    }
    if (node_name.rfind("patch:", 0) == 0)
    {
        return node_name.substr(std::string("patch:").size());
    }
    return node_name;
}

bool mappingTargetsNode(const SurfaceBoundaryMapping& mapping, const std::string& node_name)
{
    if (mapping.node_id == node_name)
    {
        return true;
    }
    return mapping.node_id == logicalNodeIdFromRuntimeName(node_name);
}

std::unordered_map<std::string, ExternalFieldSolveNodeResult>
loadExternalFieldSolveResultMap(const std::filesystem::path& result_path)
{
    std::unordered_map<std::string, ExternalFieldSolveNodeResult> results;
    if (result_path.empty() || !std::filesystem::exists(result_path))
    {
        return results;
    }

    std::ifstream input(result_path);
    if (!input.is_open())
    {
        return results;
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    const std::string content = buffer.str();

    const auto extract_string_field = [](const std::string& object_text,
                                         const std::string& key) -> std::string {
        const auto key_pos = object_text.find("\"" + key + "\"");
        if (key_pos == std::string::npos)
        {
            return {};
        }
        const auto colon_pos = object_text.find(':', key_pos);
        if (colon_pos == std::string::npos)
        {
            return {};
        }
        const auto first_quote = object_text.find('"', colon_pos + 1);
        if (first_quote == std::string::npos)
        {
            return {};
        }
        const auto second_quote = object_text.find('"', first_quote + 1);
        if (second_quote == std::string::npos)
        {
            return {};
        }
        return object_text.substr(first_quote + 1, second_quote - first_quote - 1);
    };
    const auto extract_number_field = [](const std::string& object_text, const std::string& key,
                                         double fallback_value) -> double {
        const auto key_pos = object_text.find("\"" + key + "\"");
        if (key_pos == std::string::npos)
        {
            return fallback_value;
        }
        const auto colon_pos = object_text.find(':', key_pos);
        if (colon_pos == std::string::npos)
        {
            return fallback_value;
        }
        const auto number_start = object_text.find_first_of("+-0123456789.", colon_pos + 1);
        if (number_start == std::string::npos)
        {
            return fallback_value;
        }
        const auto number_end =
            object_text.find_first_not_of("+-0123456789.eE", number_start);
        return std::stod(object_text.substr(number_start, number_end - number_start));
    };
    const std::string schema_version = extract_string_field(content, "schema_version");
    if (!schema_version.empty() && schema_version != "scdat.external_field_result.v1")
    {
        return results;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto node_key_pos = content.find("\"node_id\"", search_pos);
        const auto group_key_pos = content.find("\"boundary_group_id\"", search_pos);
        std::size_t key_pos = std::string::npos;
        if (node_key_pos != std::string::npos && group_key_pos != std::string::npos)
        {
            key_pos = std::min(node_key_pos, group_key_pos);
        }
        else if (node_key_pos != std::string::npos)
        {
            key_pos = node_key_pos;
        }
        else if (group_key_pos != std::string::npos)
        {
            key_pos = group_key_pos;
        }
        if (key_pos == std::string::npos)
        {
            break;
        }

        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        const bool has_reference_potential =
            object_text.find("\"reference_potential_v\"") != std::string::npos;
        const bool has_normal_field =
            object_text.find("\"normal_field_v_per_m\"") != std::string::npos;
        const bool has_local_charge =
            object_text.find("\"local_charge_density_c_per_m3\"") != std::string::npos;
        const bool has_capacitance_scale =
            object_text.find("\"capacitance_scale\"") != std::string::npos;
        if (!(has_reference_potential || has_normal_field || has_local_charge ||
              has_capacitance_scale))
        {
            search_pos = object_end + 1;
            continue;
        }

        ExternalFieldSolveNodeResult result;
        result.node_id = extract_string_field(object_text, "node_id");
        if (result.node_id.empty())
        {
            result.node_id = extract_string_field(object_text, "boundary_group_id");
        }
        result.reference_potential_v = extract_number_field(
            object_text, "reference_potential_v", std::numeric_limits<double>::quiet_NaN());
        result.normal_field_v_per_m = extract_number_field(
            object_text, "normal_field_v_per_m", std::numeric_limits<double>::quiet_NaN());
        result.local_charge_density_c_per_m3 = extract_number_field(
            object_text, "local_charge_density_c_per_m3", std::numeric_limits<double>::quiet_NaN());
        result.capacitance_scale = extract_number_field(
            object_text, "capacitance_scale", std::numeric_limits<double>::quiet_NaN());
        if (!result.node_id.empty())
        {
            results[result.node_id] = result;
        }
        search_pos = object_end + 1;
    }
    return results;
}

std::unordered_map<std::string, ExternalVolumeSolveCellResult>
loadExternalVolumeSolveResultMap(const std::filesystem::path& result_path)
{
    std::unordered_map<std::string, ExternalVolumeSolveCellResult> results;
    if (result_path.empty() || !std::filesystem::exists(result_path))
    {
        return results;
    }

    std::ifstream input(result_path);
    if (!input.is_open())
    {
        return results;
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    const std::string content = buffer.str();

    const auto extract_string_field = [](const std::string& object_text,
                                         const std::string& key) -> std::string {
        const auto key_pos = object_text.find("\"" + key + "\"");
        if (key_pos == std::string::npos)
        {
            return {};
        }
        const auto colon_pos = object_text.find(':', key_pos);
        if (colon_pos == std::string::npos)
        {
            return {};
        }
        const auto first_quote = object_text.find('"', colon_pos + 1);
        if (first_quote == std::string::npos)
        {
            return {};
        }
        const auto second_quote = object_text.find('"', first_quote + 1);
        if (second_quote == std::string::npos)
        {
            return {};
        }
        return object_text.substr(first_quote + 1, second_quote - first_quote - 1);
    };
    const auto extract_number_field = [](const std::string& object_text,
                                         const std::string& key,
                                         double fallback_value) -> double {
        const auto key_pos = object_text.find("\"" + key + "\"");
        if (key_pos == std::string::npos)
        {
            return fallback_value;
        }
        const auto colon_pos = object_text.find(':', key_pos);
        if (colon_pos == std::string::npos)
        {
            return fallback_value;
        }
        const auto number_start = object_text.find_first_of("+-0123456789.", colon_pos + 1);
        if (number_start == std::string::npos)
        {
            return fallback_value;
        }
        const auto number_end =
            object_text.find_first_not_of("+-0123456789.eE", number_start);
        return std::stod(object_text.substr(number_start, number_end - number_start));
    };

    const std::string schema_version = extract_string_field(content, "schema_version");
    if (!schema_version.empty() && schema_version != "scdat.external_volume_result.v1")
    {
        return results;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto cell_key_pos = content.find("\"cell_id\"", search_pos);
        const auto group_key_pos = content.find("\"boundary_group_id\"", search_pos);
        std::size_t key_pos = std::string::npos;
        if (cell_key_pos != std::string::npos && group_key_pos != std::string::npos)
        {
            key_pos = std::min(cell_key_pos, group_key_pos);
        }
        else if (cell_key_pos != std::string::npos)
        {
            key_pos = cell_key_pos;
        }
        else if (group_key_pos != std::string::npos)
        {
            key_pos = group_key_pos;
        }
        if (key_pos == std::string::npos)
        {
            break;
        }
        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        ExternalVolumeSolveCellResult result;
        result.cell_id = extract_string_field(object_text, "cell_id");
        result.boundary_group_id =
            extract_string_field(object_text, "boundary_group_id");
        result.potential_v = extract_number_field(object_text, "potential_v", 0.0);
        result.reference_potential_v =
            extract_number_field(object_text, "reference_potential_v", result.potential_v);
        result.normal_field_v_per_m =
            extract_number_field(object_text, "normal_field_v_per_m", 0.0);
        result.local_charge_density_c_per_m3 =
            extract_number_field(object_text, "local_charge_density_c_per_m3", 0.0);
        result.capacitance_scale =
            extract_number_field(object_text, "capacitance_scale", 1.0);
        result.coupling_gain = extract_number_field(object_text, "coupling_gain", 0.0);
        result.projection_weight_scale =
            extract_number_field(object_text, "projection_weight_scale", 1.0);
        result.sheath_length_scale =
            extract_number_field(object_text, "sheath_length_scale", 1.0);
        if (result.cell_id.empty())
        {
            result.cell_id = result.boundary_group_id;
        }
        if (!result.cell_id.empty())
        {
            results[result.cell_id] = result;
        }
        search_pos = object_end + 1;
    }
    return results;
}

std::string boundaryGroupIdForRuntimeNode(const SurfaceChargingConfig& config,
                                          const std::string& node_name)
{
    for (const auto& mapping : config.boundary_mappings)
    {
        if (mappingTargetsNode(mapping, node_name))
        {
            return mapping.boundary_group_id;
        }
    }
    return {};
}

std::size_t linkedBoundaryFaceCountForGroup(const SurfaceChargingConfig& config,
                                            const std::string& boundary_group_id)
{
    for (const auto& group : config.body_boundary_groups)
    {
        if (group.id == boundary_group_id)
        {
            return group.boundary_face_ids.size();
        }
    }
    for (const auto& group : config.patch_boundary_groups)
    {
        if (group.id == boundary_group_id)
        {
            return group.boundary_face_ids.size();
        }
    }
    return 0;
}

double projectionWeightSumForRuntimeNode(const SurfaceChargingConfig& config,
                                         const std::string& node_name)
{
    double projection_weight_sum = 1.0;
    std::size_t linked_boundary_faces = 0;
    for (const auto& mapping : config.boundary_mappings)
    {
        if (!mappingTargetsNode(mapping, node_name))
        {
            continue;
        }
        projection_weight_sum += 1.0;
        linked_boundary_faces +=
            linkedBoundaryFaceCountForGroup(config, mapping.boundary_group_id);
    }
    projection_weight_sum += 0.1 * static_cast<double>(linked_boundary_faces);
    return std::max(1.0, projection_weight_sum);
}

struct PseudoVolumePoissonCell
{
    std::size_t node_index = 0;
    std::string cell_id;
    std::string node_name;
    std::string boundary_group_id;
    double area_m2 = 0.0;
    double characteristic_length_m = 0.0;
    double sheath_length_m = 0.0;
    double pseudo_volume_m3 = 0.0;
    double projection_weight_sum = 1.0;
    double volume_to_surface_weight = 1.0;
    std::array<double, 3> center_m{0.0, 0.0, 0.0};
    std::array<double, 3> surface_normal{1.0, 0.0, 0.0};
    double surface_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double deposited_charge_c = 0.0;
    double deposited_charge_density_c_per_m3 = 0.0;
    double solved_potential_v = 0.0;
    double reconstructed_charge_density_c_per_m3 = 0.0;
    double normal_field_v_per_m = 0.0;
    double poisson_residual_v_m = 0.0;
    double coupling_gain = 0.0;
};

struct PseudoVolumePoissonEdge
{
    std::size_t from = 0;
    std::size_t to = 0;
    double coefficient_m = 0.0;
    double shared_face_area_m2 = 0.0;
    double distance_m = 0.0;
};

struct PseudoVolumePoissonSolution
{
    std::vector<PseudoVolumePoissonCell> cells;
    std::vector<PseudoVolumePoissonEdge> edges;
    std::string solver_mode = "dense";
    int self_consistent_iterations_completed = 0;
    int linear_iterations_completed = 0;
    bool converged = false;
    double linear_residual_norm = 0.0;
    double max_delta_v = 0.0;
    std::size_t matrix_nonzeros = 0;
};

struct LinearSolveResult
{
    std::vector<double> solution;
    int iterations = 0;
    bool converged = false;
    double residual_norm = 0.0;
};

struct SparseLinearEntry
{
    std::size_t column = 0;
    double value = 0.0;
};

struct SparseLinearRow
{
    double diagonal = 0.0;
    std::vector<SparseLinearEntry> off_diagonal;
};

struct ExternalVolumeMeshCellInput
{
    std::string cell_id;
    std::string node_id;
    std::string boundary_group_id;
    double node_area_m2 = 0.0;
    double characteristic_length_m = 0.0;
    double center_x_m = 0.0;
    double center_y_m = 0.0;
    double center_z_m = 0.0;
    double cell_volume_m3 = 0.0;
    double initial_potential_v = 0.0;
    double initial_charge_density_c_per_m3 = 0.0;
};

struct ExternalVolumeMeshNeighborInput
{
    std::string source_cell_id;
    std::string target_cell_id;
    double conductance_s = 0.0;
    double resistance_ohm = 0.0;
    double shared_face_area_m2 = 0.0;
    double face_distance_m = 0.0;
    double permittivity_scale = 1.0;
    double face_center_x_m = 0.0;
    double face_center_y_m = 0.0;
    double face_center_z_m = 0.0;
};

struct ExternalVolumeMeshFaceInput
{
    std::string face_id;
    std::string boundary_group_id;
    std::string node_id;
    double area_m2 = 0.0;
    double center_x_m = 0.0;
    double center_y_m = 0.0;
    double center_z_m = 0.0;
    double normal_x = 1.0;
    double normal_y = 0.0;
    double normal_z = 0.0;
};

struct ExternalVolumeMeshCellFaceInput
{
    std::string cell_id;
    std::string face_id;
    std::string role;
    double projection_weight = 1.0;
};

struct ExternalSurfaceProjectionInput
{
    std::string node_id;
    std::string cell_id;
    double surface_to_volume_weight = 1.0;
    double volume_to_surface_weight = 1.0;
};

double squaredNorm3(const std::array<double, 3>& value)
{
    return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

double norm3(const std::array<double, 3>& value)
{
    return std::sqrt(squaredNorm3(value));
}

std::array<double, 3> subtract3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::array<double, 3> add3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::array<double, 3> scale3(const std::array<double, 3>& value, double scale)
{
    return {value[0] * scale, value[1] * scale, value[2] * scale};
}

std::array<double, 3> normalize3(const std::array<double, 3>& value,
                                 const std::array<double, 3>& fallback)
{
    const double magnitude = norm3(value);
    if (magnitude <= 1.0e-15)
    {
        return fallback;
    }
    return scale3(value, 1.0 / magnitude);
}

double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> defaultPseudoCellCenter(const SurfaceCircuitModel* circuit_model,
                                              std::size_t node_index)
{
    const double radius_m = 0.15 + 0.03 * static_cast<double>(node_index);
    if (circuit_model == nullptr || node_index >= circuit_model->nodeCount())
    {
        return {radius_m, 0.0, 0.0};
    }
    if (!circuit_model->nodeIsPatch(node_index))
    {
        const double angle = 0.45 * static_cast<double>(node_index);
        return {radius_m * std::cos(angle), radius_m * std::sin(angle),
                0.02 * static_cast<double>(node_index)};
    }
    const double angle = 0.75 * static_cast<double>(node_index);
    return {radius_m * std::cos(angle), radius_m * std::sin(angle), 0.2 * radius_m};
}

std::array<double, 3> defaultPseudoSurfaceNormal(const SurfaceCircuitModel* circuit_model,
                                                 std::size_t node_index)
{
    if (circuit_model == nullptr || node_index >= circuit_model->nodeCount())
    {
        return {1.0, 0.0, 0.0};
    }
    if (!circuit_model->nodeIsPatch(node_index))
    {
        return {0.0, 0.0, 1.0};
    }
    const auto center = normalize3(defaultPseudoCellCenter(circuit_model, node_index), {1.0, 0.0, 0.0});
    return normalize3({center[0], center[1], std::max(0.1, center[2])}, {1.0, 0.0, 0.0});
}

std::string extractJsonStringField(const std::string& object_text, const std::string& key)
{
    const auto key_pos = object_text.find("\"" + key + "\"");
    if (key_pos == std::string::npos)
    {
        return {};
    }
    const auto colon_pos = object_text.find(':', key_pos);
    if (colon_pos == std::string::npos)
    {
        return {};
    }
    const auto first_quote = object_text.find('"', colon_pos + 1);
    if (first_quote == std::string::npos)
    {
        return {};
    }
    const auto second_quote = object_text.find('"', first_quote + 1);
    if (second_quote == std::string::npos)
    {
        return {};
    }
    return object_text.substr(first_quote + 1, second_quote - first_quote - 1);
}

double extractJsonNumberField(const std::string& object_text, const std::string& key,
                              double fallback_value)
{
    const auto key_pos = object_text.find("\"" + key + "\"");
    if (key_pos == std::string::npos)
    {
        return fallback_value;
    }
    const auto colon_pos = object_text.find(':', key_pos);
    if (colon_pos == std::string::npos)
    {
        return fallback_value;
    }
    const auto number_start = object_text.find_first_of("+-0123456789.", colon_pos + 1);
    if (number_start == std::string::npos)
    {
        return fallback_value;
    }
    const auto number_end = object_text.find_first_not_of("+-0123456789.eE", number_start);
    return std::stod(object_text.substr(number_start, number_end - number_start));
}

std::string loadFileToString(const std::filesystem::path& path)
{
    if (path.empty() || !std::filesystem::exists(path))
    {
        return {};
    }
    std::ifstream input(path);
    if (!input.is_open())
    {
        return {};
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::vector<ExternalVolumeMeshCellInput> loadExternalVolumeMeshCells(
    const std::filesystem::path& mesh_path)
{
    std::vector<ExternalVolumeMeshCellInput> cells;
    const std::string content = loadFileToString(mesh_path);
    if (content.empty())
    {
        return cells;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto key_pos = content.find("\"cell_id\"", search_pos);
        if (key_pos == std::string::npos)
        {
            break;
        }
        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        ExternalVolumeMeshCellInput cell;
        cell.cell_id = extractJsonStringField(object_text, "cell_id");
        cell.node_id = extractJsonStringField(object_text, "node_id");
        cell.boundary_group_id = extractJsonStringField(object_text, "boundary_group_id");
        cell.node_area_m2 = extractJsonNumberField(object_text, "node_area_m2", 0.0);
        cell.characteristic_length_m =
            extractJsonNumberField(object_text, "characteristic_length_m", 0.0);
        cell.center_x_m = extractJsonNumberField(object_text, "center_x_m", 0.0);
        cell.center_y_m = extractJsonNumberField(object_text, "center_y_m", 0.0);
        cell.center_z_m = extractJsonNumberField(object_text, "center_z_m", 0.0);
        cell.cell_volume_m3 = extractJsonNumberField(object_text, "cell_volume_m3", 0.0);
        cell.initial_potential_v = extractJsonNumberField(object_text, "initial_potential_v", 0.0);
        cell.initial_charge_density_c_per_m3 =
            extractJsonNumberField(object_text, "initial_charge_density_c_per_m3", 0.0);
        if (!cell.cell_id.empty())
        {
            cells.push_back(cell);
        }
        search_pos = object_end + 1;
    }
    return cells;
}

std::vector<ExternalVolumeMeshNeighborInput> loadExternalVolumeMeshNeighbors(
    const std::filesystem::path& mesh_path)
{
    std::vector<ExternalVolumeMeshNeighborInput> neighbors;
    const std::string content = loadFileToString(mesh_path);
    if (content.empty())
    {
        return neighbors;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto key_pos = content.find("\"source_cell_id\"", search_pos);
        if (key_pos == std::string::npos)
        {
            break;
        }
        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        ExternalVolumeMeshNeighborInput neighbor;
        neighbor.source_cell_id = extractJsonStringField(object_text, "source_cell_id");
        neighbor.target_cell_id = extractJsonStringField(object_text, "target_cell_id");
        neighbor.conductance_s = extractJsonNumberField(object_text, "conductance_s", 0.0);
        neighbor.resistance_ohm = extractJsonNumberField(object_text, "resistance_ohm", 0.0);
        neighbor.shared_face_area_m2 =
            extractJsonNumberField(object_text, "shared_face_area_m2", 0.0);
        neighbor.face_distance_m = extractJsonNumberField(object_text, "face_distance_m", 0.0);
        neighbor.permittivity_scale =
            extractJsonNumberField(object_text, "permittivity_scale", 1.0);
        neighbor.face_center_x_m = extractJsonNumberField(object_text, "face_center_x_m", 0.0);
        neighbor.face_center_y_m = extractJsonNumberField(object_text, "face_center_y_m", 0.0);
        neighbor.face_center_z_m = extractJsonNumberField(object_text, "face_center_z_m", 0.0);
        if (!neighbor.source_cell_id.empty() && !neighbor.target_cell_id.empty())
        {
            neighbors.push_back(neighbor);
        }
        search_pos = object_end + 1;
    }
    return neighbors;
}

std::vector<ExternalVolumeMeshFaceInput> loadExternalVolumeMeshFaces(
    const std::filesystem::path& mesh_path)
{
    std::vector<ExternalVolumeMeshFaceInput> faces;
    const std::string content = loadFileToString(mesh_path);
    if (content.empty())
    {
        return faces;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto key_pos = content.find("\"face_id\"", search_pos);
        if (key_pos == std::string::npos)
        {
            break;
        }
        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        ExternalVolumeMeshFaceInput face;
        face.face_id = extractJsonStringField(object_text, "face_id");
        face.boundary_group_id = extractJsonStringField(object_text, "boundary_group_id");
        face.node_id = extractJsonStringField(object_text, "node_id");
        face.area_m2 = extractJsonNumberField(object_text, "area_m2", 0.0);
        face.center_x_m = extractJsonNumberField(object_text, "center_x_m", 0.0);
        face.center_y_m = extractJsonNumberField(object_text, "center_y_m", 0.0);
        face.center_z_m = extractJsonNumberField(object_text, "center_z_m", 0.0);
        face.normal_x = extractJsonNumberField(object_text, "normal_x", 1.0);
        face.normal_y = extractJsonNumberField(object_text, "normal_y", 0.0);
        face.normal_z = extractJsonNumberField(object_text, "normal_z", 0.0);
        if (!face.face_id.empty())
        {
            faces.push_back(face);
        }
        search_pos = object_end + 1;
    }
    return faces;
}

std::vector<ExternalSurfaceProjectionInput> loadExternalSurfaceProjectionInputs(
    const std::filesystem::path& projection_path)
{
    std::vector<ExternalSurfaceProjectionInput> entries;
    const std::string content = loadFileToString(projection_path);
    if (content.empty())
    {
        return entries;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto key_pos = content.find("\"surface_to_volume_weight\"", search_pos);
        if (key_pos == std::string::npos)
        {
            break;
        }
        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        ExternalSurfaceProjectionInput entry;
        entry.node_id = extractJsonStringField(object_text, "node_id");
        entry.cell_id = extractJsonStringField(object_text, "cell_id");
        entry.surface_to_volume_weight =
            extractJsonNumberField(object_text, "surface_to_volume_weight", 1.0);
        entry.volume_to_surface_weight =
            extractJsonNumberField(object_text, "volume_to_surface_weight", 1.0);
        if (!entry.node_id.empty() && !entry.cell_id.empty())
        {
            entries.push_back(entry);
        }
        search_pos = object_end + 1;
    }
    return entries;
}

std::vector<ExternalVolumeMeshCellFaceInput> loadExternalVolumeMeshCellFaces(
    const std::filesystem::path& mesh_path)
{
    std::vector<ExternalVolumeMeshCellFaceInput> entries;
    const std::string content = loadFileToString(mesh_path);
    if (content.empty())
    {
        return entries;
    }

    std::size_t search_pos = 0;
    while (true)
    {
        const auto key_pos = content.find("\"role\"", search_pos);
        if (key_pos == std::string::npos)
        {
            break;
        }
        const auto object_start = content.rfind('{', key_pos);
        const auto object_end = content.find('}', key_pos);
        if (object_start == std::string::npos || object_end == std::string::npos ||
            object_end <= object_start)
        {
            break;
        }
        const std::string object_text =
            content.substr(object_start, object_end - object_start + 1);
        ExternalVolumeMeshCellFaceInput entry;
        entry.cell_id = extractJsonStringField(object_text, "cell_id");
        entry.face_id = extractJsonStringField(object_text, "face_id");
        entry.role = extractJsonStringField(object_text, "role");
        entry.projection_weight = extractJsonNumberField(object_text, "projection_weight", 1.0);
        if (!entry.cell_id.empty() && !entry.face_id.empty() && !entry.role.empty())
        {
            entries.push_back(entry);
        }
        search_pos = object_end + 1;
    }
    return entries;
}

std::vector<double> solveDenseLinearSystem(std::vector<std::vector<double>> matrix,
                                           std::vector<double> rhs)
{
    const std::size_t n = rhs.size();
    if (n == 0)
    {
        return {};
    }

    for (std::size_t pivot = 0; pivot < n; ++pivot)
    {
        std::size_t best_row = pivot;
        double best_value = std::abs(matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < n; ++row)
        {
            const double candidate = std::abs(matrix[row][pivot]);
            if (candidate > best_value)
            {
                best_value = candidate;
                best_row = row;
            }
        }
        if (best_row != pivot)
        {
            std::swap(matrix[pivot], matrix[best_row]);
            std::swap(rhs[pivot], rhs[best_row]);
        }

        const double diagonal = matrix[pivot][pivot];
        if (std::abs(diagonal) < 1.0e-20)
        {
            continue;
        }

        for (std::size_t row = pivot + 1; row < n; ++row)
        {
            const double factor = matrix[row][pivot] / diagonal;
            if (std::abs(factor) < 1.0e-20)
            {
                continue;
            }
            for (std::size_t column = pivot; column < n; ++column)
            {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    std::vector<double> solution(n, 0.0);
    for (std::size_t back = n; back-- > 0;)
    {
        double value = rhs[back];
        for (std::size_t column = back + 1; column < n; ++column)
        {
            value -= matrix[back][column] * solution[column];
        }
        const double diagonal = matrix[back][back];
        solution[back] = std::abs(diagonal) > 1.0e-20 ? value / diagonal : 0.0;
    }
    return solution;
}

double denseLinearResidualNorm(const std::vector<std::vector<double>>& matrix,
                               const std::vector<double>& rhs,
                               const std::vector<double>& solution)
{
    if (matrix.size() != rhs.size() || rhs.size() != solution.size())
    {
        return 0.0;
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i)
    {
        double residual = rhs[i];
        for (std::size_t j = 0; j < solution.size(); ++j)
        {
            residual -= matrix[i][j] * solution[j];
        }
        sum += residual * residual;
    }
    return std::sqrt(sum / std::max<std::size_t>(1, rhs.size()));
}

std::vector<SparseLinearRow> buildSparseLinearRows(const std::vector<std::vector<double>>& matrix,
                                                   std::size_t* nonzero_count)
{
    std::vector<SparseLinearRow> rows(matrix.size());
    std::size_t nnz = 0;
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        rows[i].diagonal = i < matrix[i].size() ? matrix[i][i] : 0.0;
        for (std::size_t j = 0; j < matrix[i].size(); ++j)
        {
            if (j == i)
            {
                if (std::abs(matrix[i][j]) > 1.0e-20)
                {
                    ++nnz;
                }
                continue;
            }
            if (std::abs(matrix[i][j]) <= 1.0e-20)
            {
                continue;
            }
            rows[i].off_diagonal.push_back({j, matrix[i][j]});
            ++nnz;
        }
    }
    if (nonzero_count != nullptr)
    {
        *nonzero_count = nnz;
    }
    return rows;
}

double sparseLinearResidualNorm(const std::vector<SparseLinearRow>& rows,
                                const std::vector<double>& rhs,
                                const std::vector<double>& solution)
{
    if (rows.size() != rhs.size() || rhs.size() != solution.size())
    {
        return 0.0;
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i)
    {
        double residual = rhs[i] - rows[i].diagonal * solution[i];
        for (const auto& entry : rows[i].off_diagonal)
        {
            residual -= entry.value * solution[entry.column];
        }
        sum += residual * residual;
    }
    return std::sqrt(sum / std::max<std::size_t>(1, rhs.size()));
}

LinearSolveResult solveIterativeLinearSystem(const std::vector<std::vector<double>>& matrix,
                                             const std::vector<double>& rhs,
                                             const std::vector<double>& initial_guess,
                                             int max_iterations, double tolerance,
                                             double relaxation)
{
    const std::size_t n = rhs.size();
    if (n == 0)
    {
        return {};
    }

    std::size_t sparse_nnz = 0;
    const auto sparse_rows = buildSparseLinearRows(matrix, &sparse_nnz);
    const double effective_tolerance = std::max(1.0e-12, tolerance);
    const double relaxation_factor = std::clamp(relaxation, 0.2, 1.0);
    LinearSolveResult result;
    result.solution = initial_guess.size() == n ? initial_guess : std::vector<double>(n, 0.0);
    for (int iteration = 0; iteration < std::max(1, max_iterations); ++iteration)
    {
        double max_delta = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const double diagonal = sparse_rows[i].diagonal;
            if (std::abs(diagonal) < 1.0e-20)
            {
                continue;
            }
            double sigma = rhs[i];
            for (const auto& entry : sparse_rows[i].off_diagonal)
            {
                sigma -= entry.value * result.solution[entry.column];
            }
            const double updated = sigma / diagonal;
            max_delta = std::max(max_delta, std::abs(updated - result.solution[i]));
            result.solution[i] =
                relaxation_factor * updated + (1.0 - relaxation_factor) * result.solution[i];
        }
        result.iterations = iteration + 1;
        if (max_delta <= effective_tolerance)
        {
            result.converged = true;
            break;
        }
    }
    result.residual_norm = sparseLinearResidualNorm(sparse_rows, rhs, result.solution);
    return result;
}

PseudoVolumePoissonSolution solvePseudoVolumePoisson(const SurfaceChargingConfig& config,
                                                     const SurfaceCircuitModel* circuit_model,
                                                     const SurfaceModelRuntimeState& state,
                                                     const std::unordered_map<std::string, double>* previous_potentials_v,
                                                     const std::unordered_map<std::string, double>* previous_charge_density_c_per_m3)
{
    PseudoVolumePoissonSolution solution;
    if (circuit_model == nullptr || circuit_model->nodeCount() == 0)
    {
        return solution;
    }

    const double global_sheath_length_m =
        std::clamp(std::max(1.0e-6, state.effective_sheath_length_m),
                   std::max(1.0e-6, config.minimum_sheath_length_m),
                   std::max(config.minimum_sheath_length_m, config.maximum_sheath_length_m));
    const double global_reference_potential_v =
        std::isfinite(state.field_solver_reference_potential_v)
            ? state.field_solver_reference_potential_v
            : state.reference_plasma_potential_v;
    const auto external_cells =
        loadExternalVolumeMeshCells(config.external_volume_mesh_path);
    const auto external_neighbors =
        loadExternalVolumeMeshNeighbors(config.external_volume_mesh_path);
    const auto external_faces =
        loadExternalVolumeMeshFaces(config.external_volume_mesh_path);
    const auto external_cell_faces =
        loadExternalVolumeMeshCellFaces(config.external_volume_mesh_path);
    const auto external_projections =
        loadExternalSurfaceProjectionInputs(config.external_surface_volume_projection_path);
    std::unordered_map<std::string, ExternalVolumeMeshCellInput> external_cell_by_node;
    std::unordered_map<std::string, ExternalVolumeMeshCellInput> external_cell_by_cell_id;
    for (const auto& cell : external_cells)
    {
        if (!cell.node_id.empty())
        {
            external_cell_by_node[cell.node_id] = cell;
            external_cell_by_node[logicalNodeIdFromRuntimeName(cell.node_id)] = cell;
        }
        external_cell_by_cell_id[cell.cell_id] = cell;
    }
    std::unordered_map<std::string, double> face_area_sum_by_group;
    std::unordered_map<std::string, std::array<double, 3>> face_center_sum_by_group;
    std::unordered_map<std::string, std::array<double, 3>> face_normal_sum_by_group;
    std::unordered_map<std::string, double> face_count_by_group;
    std::unordered_map<std::string, ExternalVolumeMeshFaceInput> external_face_by_id;
    for (const auto& face : external_faces)
    {
        external_face_by_id[face.face_id] = face;
        face_area_sum_by_group[face.boundary_group_id] += std::max(0.0, face.area_m2);
        face_center_sum_by_group[face.boundary_group_id] =
            add3(face_center_sum_by_group[face.boundary_group_id],
                 {face.center_x_m, face.center_y_m, face.center_z_m});
        face_normal_sum_by_group[face.boundary_group_id] =
            add3(face_normal_sum_by_group[face.boundary_group_id],
                 normalize3({face.normal_x, face.normal_y, face.normal_z}, {1.0, 0.0, 0.0}));
        face_count_by_group[face.boundary_group_id] += 1.0;
    }
    std::unordered_map<std::string, std::vector<ExternalVolumeMeshCellFaceInput>> cell_faces_by_cell_id;
    for (const auto& entry : external_cell_faces)
    {
        cell_faces_by_cell_id[entry.cell_id].push_back(entry);
    }
    std::unordered_map<std::string, ExternalSurfaceProjectionInput> projection_by_node;
    std::unordered_map<std::string, ExternalSurfaceProjectionInput> projection_by_cell_id;
    for (const auto& projection : external_projections)
    {
        projection_by_node[projection.node_id] = projection;
        projection_by_node[logicalNodeIdFromRuntimeName(projection.node_id)] = projection;
        projection_by_cell_id[projection.cell_id] = projection;
    }

    solution.cells.reserve(circuit_model->nodeCount());
    for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
    {
        const std::string node_name = circuit_model->nodeName(node_index);
        const double area_m2 = std::max(1.0e-12, circuit_model->nodeAreaM2(node_index));
        const auto external_cell_it = external_cell_by_node.find(node_name);
        const auto projection_it = projection_by_node.find(node_name);
        const double projection_weight_sum = projection_it != projection_by_node.end()
                                                 ? std::max(1.0, projection_it->second.surface_to_volume_weight)
                                                 : projectionWeightSumForRuntimeNode(config, node_name);
        const double volume_to_surface_weight =
            projection_it != projection_by_node.end()
                ? std::max(1.0, projection_it->second.volume_to_surface_weight)
                : projection_weight_sum;
        const std::string boundary_group_id =
            external_cell_it != external_cell_by_node.end() &&
                    !external_cell_it->second.boundary_group_id.empty()
                ? external_cell_it->second.boundary_group_id
                : boundaryGroupIdForRuntimeNode(config, node_name);
        const double face_derived_area_m2 =
            face_area_sum_by_group.contains(boundary_group_id)
                ? std::max(1.0e-12, face_area_sum_by_group[boundary_group_id])
                : 0.0;
        const double effective_area_m2 =
            external_cell_it != external_cell_by_node.end() &&
                    external_cell_it->second.node_area_m2 > 0.0
                ? external_cell_it->second.node_area_m2
                : (face_derived_area_m2 > 0.0 ? face_derived_area_m2 : area_m2);
        const std::array<double, 3> default_center_m =
            defaultPseudoCellCenter(circuit_model, node_index);
        const std::array<double, 3> default_normal =
            defaultPseudoSurfaceNormal(circuit_model, node_index);
        std::array<double, 3> center_m = default_center_m;
        const std::string external_cell_id =
            external_cell_it != external_cell_by_node.end() ? external_cell_it->second.cell_id
                                                            : "cell_" + std::to_string(node_index);
        auto cell_face_entries_it = cell_faces_by_cell_id.find(external_cell_id);
        std::array<double, 3> linked_face_center_sum{0.0, 0.0, 0.0};
        std::array<double, 3> linked_face_normal_sum{0.0, 0.0, 0.0};
        double linked_face_area_sum = 0.0;
        double linked_face_weight_sum = 0.0;
        if (cell_face_entries_it != cell_faces_by_cell_id.end())
        {
            for (const auto& entry : cell_face_entries_it->second)
            {
                const auto face_it = external_face_by_id.find(entry.face_id);
                if (face_it == external_face_by_id.end())
                {
                    continue;
                }
                const double face_weight =
                    std::max(0.0, face_it->second.area_m2) * std::max(1.0, entry.projection_weight);
                linked_face_center_sum =
                    add3(linked_face_center_sum,
                         scale3({face_it->second.center_x_m, face_it->second.center_y_m,
                                 face_it->second.center_z_m},
                                face_weight));
                linked_face_normal_sum =
                    add3(linked_face_normal_sum,
                         scale3(normalize3({face_it->second.normal_x, face_it->second.normal_y,
                                            face_it->second.normal_z},
                                           default_normal),
                                face_weight));
                linked_face_area_sum += std::max(0.0, face_it->second.area_m2);
                linked_face_weight_sum += face_weight;
            }
        }
        if (external_cell_it != external_cell_by_node.end() &&
            (std::abs(external_cell_it->second.center_x_m) > 0.0 ||
             std::abs(external_cell_it->second.center_y_m) > 0.0 ||
             std::abs(external_cell_it->second.center_z_m) > 0.0))
        {
            center_m = {external_cell_it->second.center_x_m, external_cell_it->second.center_y_m,
                        external_cell_it->second.center_z_m};
        }
        else if (linked_face_weight_sum > 0.0)
        {
            center_m = scale3(linked_face_center_sum, 1.0 / linked_face_weight_sum);
        }
        else if (face_count_by_group.contains(boundary_group_id) &&
                 face_count_by_group[boundary_group_id] > 0.0)
        {
            center_m = scale3(face_center_sum_by_group[boundary_group_id],
                              1.0 / face_count_by_group[boundary_group_id]);
        }
        std::array<double, 3> surface_normal = default_normal;
        if (linked_face_weight_sum > 0.0)
        {
            surface_normal =
                normalize3(scale3(linked_face_normal_sum, 1.0 / linked_face_weight_sum), default_normal);
        }
        if (face_count_by_group.contains(boundary_group_id) &&
            face_count_by_group[boundary_group_id] > 0.0)
        {
            surface_normal = normalize3(
                scale3(face_normal_sum_by_group[boundary_group_id],
                       1.0 / face_count_by_group[boundary_group_id]),
                default_normal);
        }
        double characteristic_length_m = std::max(1.0e-6, std::sqrt(effective_area_m2 / kPi));
        if (face_count_by_group.contains(boundary_group_id) &&
            face_count_by_group[boundary_group_id] > 0.0)
        {
            characteristic_length_m = std::max(
                1.0e-6,
                norm3(subtract3(
                    scale3(face_center_sum_by_group[boundary_group_id],
                           1.0 / face_count_by_group[boundary_group_id]),
                    center_m)));
        }
        if (linked_face_area_sum > 0.0)
        {
            characteristic_length_m = std::max(
                1.0e-6, std::sqrt(linked_face_area_sum / std::max(1.0, linked_face_weight_sum)));
        }
        if (external_cell_it != external_cell_by_node.end() &&
            external_cell_it->second.characteristic_length_m > 0.0)
        {
            characteristic_length_m =
                std::max(1.0e-6, external_cell_it->second.characteristic_length_m);
        }
        const double pseudo_volume_m3 =
            external_cell_it != external_cell_by_node.end() &&
                    external_cell_it->second.cell_volume_m3 > 0.0
                ? std::max(1.0e-18, external_cell_it->second.cell_volume_m3 * projection_weight_sum)
                : (external_cell_it != external_cell_by_node.end() && effective_area_m2 > 0.0
                       ? std::max(1.0e-18,
                                  effective_area_m2 * global_sheath_length_m * projection_weight_sum)
                       : std::max(1.0e-18, area_m2 * global_sheath_length_m * projection_weight_sum));
        const double surface_potential_v = circuit_model->nodePotential(node_index);
        const double reference_drop_v = surface_potential_v - global_reference_potential_v;
        const double surface_sigma_c_per_m2 =
            kEpsilon0 * reference_drop_v / std::max(1.0e-6, 0.5 * global_sheath_length_m);
        double deposited_charge_c = surface_sigma_c_per_m2 * area_m2 * projection_weight_sum;
        if (node_index == state.node_index)
        {
            deposited_charge_c += state.local_charge_density_c_per_m3 * pseudo_volume_m3;
            deposited_charge_c += kEpsilon0 * state.normal_electric_field_v_per_m * area_m2;
        }
        else if (external_cell_it != external_cell_by_node.end())
        {
            deposited_charge_c +=
                external_cell_it->second.initial_charge_density_c_per_m3 * pseudo_volume_m3;
        }

        PseudoVolumePoissonCell cell;
        cell.node_index = node_index;
        cell.cell_id = external_cell_id;
        cell.node_name = node_name;
        cell.boundary_group_id = boundary_group_id;
        cell.area_m2 = effective_area_m2;
        cell.characteristic_length_m = characteristic_length_m;
        cell.sheath_length_m = global_sheath_length_m;
        cell.pseudo_volume_m3 = pseudo_volume_m3;
        cell.projection_weight_sum = projection_weight_sum;
        cell.volume_to_surface_weight = volume_to_surface_weight;
        cell.center_m = center_m;
        cell.surface_normal = surface_normal;
        cell.surface_potential_v = surface_potential_v;
        cell.reference_potential_v = global_reference_potential_v;
        cell.deposited_charge_c = deposited_charge_c;
        cell.deposited_charge_density_c_per_m3 = deposited_charge_c / pseudo_volume_m3;
        cell.solved_potential_v =
            external_cell_it != external_cell_by_node.end() &&
                    std::isfinite(external_cell_it->second.initial_potential_v)
                ? external_cell_it->second.initial_potential_v
                : global_reference_potential_v;
        solution.cells.push_back(cell);
    }

    solution.edges.reserve(circuit_model->branchCount());
    if (!external_neighbors.empty())
    {
        std::unordered_map<std::string, std::size_t> cell_index_by_id;
        for (std::size_t i = 0; i < solution.cells.size(); ++i)
        {
            cell_index_by_id[solution.cells[i].cell_id] = i;
        }
        for (const auto& neighbor : external_neighbors)
        {
            const auto from_it = cell_index_by_id.find(neighbor.source_cell_id);
            const auto to_it = cell_index_by_id.find(neighbor.target_cell_id);
            if (from_it == cell_index_by_id.end() || to_it == cell_index_by_id.end())
            {
                continue;
            }
            const auto from_node = from_it->second;
            const auto to_node = to_it->second;
            const double conductance_s = std::max(
                0.0, neighbor.conductance_s > 0.0
                         ? neighbor.conductance_s
                         : (neighbor.resistance_ohm > 0.0 ? 1.0 / neighbor.resistance_ohm : 0.0));
            const auto center_delta =
                subtract3(solution.cells[to_node].center_m, solution.cells[from_node].center_m);
            const double center_distance_m = std::max(1.0e-6, norm3(center_delta));
            const auto center_direction =
                normalize3(center_delta, solution.cells[from_node].surface_normal);
            const double directional_coupling =
                std::max(0.1, std::abs(dot3(solution.cells[from_node].surface_normal, center_direction)));
            const double shared_face_area_m2 = neighbor.shared_face_area_m2 > 0.0
                                                   ? neighbor.shared_face_area_m2
                                                   : 0.5 * std::min(solution.cells[from_node].area_m2,
                                                                    solution.cells[to_node].area_m2) *
                                                         directional_coupling;
            const double distance_m = neighbor.face_distance_m > 0.0
                                          ? neighbor.face_distance_m
                                          : std::max(1.0e-6,
                                                     std::max(
                                                         center_distance_m,
                                                         0.5 * (solution.cells[from_node].characteristic_length_m +
                                                                solution.cells[to_node].characteristic_length_m) +
                                                             global_sheath_length_m));
            const double conductance_gain =
                1.0 + 0.25 * conductance_s / std::max(1.0e-12, conductance_s + 1.0e-9);
            const double permittivity_scale =
                std::clamp(neighbor.permittivity_scale, 0.1, 10.0);
            const double coefficient_m =
                std::max(1.0e-12, shared_face_area_m2 / distance_m) * conductance_gain *
                permittivity_scale *
                std::clamp(0.75 + 0.25 * directional_coupling, 0.25, 1.0);
            solution.edges.push_back(
                {from_node, to_node, coefficient_m, shared_face_area_m2, distance_m});
        }
    }
    else
    {
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            if (from_node >= solution.cells.size() || to_node >= solution.cells.size())
            {
                continue;
            }
            const double conductance_s = std::max(0.0, circuit_model->branchConductanceS(branch_index));
            const auto center_delta =
                subtract3(solution.cells[to_node].center_m, solution.cells[from_node].center_m);
            const double center_distance_m = std::max(1.0e-6, norm3(center_delta));
            const auto center_direction =
                normalize3(center_delta, solution.cells[from_node].surface_normal);
            const double directional_coupling =
                std::max(0.1, std::abs(dot3(solution.cells[from_node].surface_normal, center_direction)));
            const double shared_face_area_m2 =
                0.5 * std::min(solution.cells[from_node].area_m2, solution.cells[to_node].area_m2) *
                directional_coupling;
            const double distance_m =
                std::max(1.0e-6,
                         std::max(center_distance_m,
                                  0.5 * (solution.cells[from_node].characteristic_length_m +
                                         solution.cells[to_node].characteristic_length_m) +
                                      global_sheath_length_m));
            const double conductance_gain =
                1.0 + 0.25 * conductance_s / std::max(1.0e-12, conductance_s + 1.0e-9);
            const double coefficient_m =
                std::max(1.0e-12, shared_face_area_m2 / distance_m) * conductance_gain *
                std::clamp(0.75 + 0.25 * directional_coupling, 0.25, 1.0);
            solution.edges.push_back(
                {from_node, to_node, coefficient_m, shared_face_area_m2, distance_m});
        }
    }

    const std::size_t cell_count = solution.cells.size();
    std::vector<double> surface_coefficients(cell_count, 0.0);
    std::vector<double> reference_coefficients(cell_count, 0.0);
    std::vector<double> temporal_coefficients(cell_count, 0.0);
    std::vector<double> current_charge_density_c_per_m3(cell_count, 0.0);
    std::vector<double> previous_volume_potentials(cell_count, 0.0);

    for (std::size_t i = 0; i < cell_count; ++i)
    {
        auto& cell = solution.cells[i];
        const double surface_coefficient_m =
            cell.area_m2 / std::max(1.0e-6, 0.5 * cell.sheath_length_m);
        const double reference_coefficient_m =
            0.35 * cell.area_m2 /
            std::max(1.0e-6, cell.sheath_length_m + cell.characteristic_length_m);
        double temporal_coefficient_m = 0.0;
        double previous_volume_potential_v = cell.reference_potential_v;
        if (previous_potentials_v != nullptr)
        {
            if (const auto it = previous_potentials_v->find(cell.cell_id);
                it != previous_potentials_v->end() && std::isfinite(it->second))
            {
                previous_volume_potential_v = it->second;
                temporal_coefficient_m = 0.20 * surface_coefficient_m;
            }
        }
        if (previous_charge_density_c_per_m3 != nullptr)
        {
            if (const auto it = previous_charge_density_c_per_m3->find(cell.cell_id);
                it != previous_charge_density_c_per_m3->end() && std::isfinite(it->second))
            {
                cell.deposited_charge_density_c_per_m3 =
                    0.6 * cell.deposited_charge_density_c_per_m3 + 0.4 * it->second;
            }
        }
        surface_coefficients[i] = surface_coefficient_m;
        reference_coefficients[i] = reference_coefficient_m;
        temporal_coefficients[i] = temporal_coefficient_m;
        current_charge_density_c_per_m3[i] = cell.deposited_charge_density_c_per_m3;
        previous_volume_potentials[i] = previous_volume_potential_v;
    }

    const int configured_iterations =
        static_cast<int>(std::clamp<std::size_t>(config.volume_self_consistent_iterations, 1, 8));
    const double convergence_tolerance_v =
        std::clamp(config.volume_self_consistent_tolerance_v, 1.0e-9, 1.0);
    const double charge_relaxation =
        std::clamp(config.volume_charge_relaxation, 0.05, 1.0);
    const double potential_relaxation =
        std::clamp(config.volume_potential_relaxation, 0.05, 1.0);
    const int linear_max_iterations =
        static_cast<int>(std::clamp<std::size_t>(config.volume_linear_max_iterations, 8, 4096));
    const double linear_tolerance = convergence_tolerance_v *
                                    std::clamp(config.volume_linear_tolerance_scale, 0.01, 2.0);
    const double linear_relaxation =
        std::clamp(config.volume_linear_relaxation, 0.2, 1.0);
    const int geometry_bonus_iterations =
        std::clamp((external_cells.empty() ? 0 : 1) + (external_cell_faces.empty() ? 0 : 1), 0, 2);
    const int self_consistent_iterations =
        std::clamp(std::max(configured_iterations, 2 + geometry_bonus_iterations), 1, 8);
    std::vector<double> solved_potentials_v(cell_count, 0.0);
    std::vector<double> previous_iteration_potentials_v = previous_volume_potentials;
    int last_linear_iterations = 0;
    double last_linear_residual_norm = 0.0;
    bool last_linear_converged = false;
    double final_max_delta_v = 0.0;
    for (int iteration = 0; iteration < self_consistent_iterations; ++iteration)
    {
        std::vector<std::vector<double>> matrix(cell_count, std::vector<double>(cell_count, 0.0));
        std::vector<double> rhs(cell_count, 0.0);

        for (std::size_t i = 0; i < cell_count; ++i)
        {
            const auto& cell = solution.cells[i];
            matrix[i][i] +=
                surface_coefficients[i] + reference_coefficients[i] + temporal_coefficients[i];
            rhs[i] += surface_coefficients[i] * cell.surface_potential_v +
                      reference_coefficients[i] * cell.reference_potential_v +
                      temporal_coefficients[i] * previous_volume_potentials[i] -
                      current_charge_density_c_per_m3[i] * cell.pseudo_volume_m3 / kEpsilon0;
        }

        for (const auto& edge : solution.edges)
        {
            matrix[edge.from][edge.from] += edge.coefficient_m;
            matrix[edge.to][edge.to] += edge.coefficient_m;
            matrix[edge.from][edge.to] -= edge.coefficient_m;
            matrix[edge.to][edge.from] -= edge.coefficient_m;
        }
        std::size_t matrix_nonzeros = 0;
        for (const auto& row : matrix)
        {
            for (const double value : row)
            {
                if (std::abs(value) > 1.0e-20)
                {
                    ++matrix_nonzeros;
                }
            }
        }
        solution.matrix_nonzeros = matrix_nonzeros;

        const bool auto_prefers_iterative =
            cell_count >= 6 || !external_cells.empty() || !external_cell_faces.empty();
        const bool use_iterative_solver =
            config.volume_linear_solver_policy == VolumeLinearSolverPolicy::IterativeOnly ||
            (config.volume_linear_solver_policy == VolumeLinearSolverPolicy::Auto &&
             auto_prefers_iterative);
        solution.solver_mode = use_iterative_solver ? "iterative" : "dense";
        if (use_iterative_solver)
        {
            auto linear_result = solveIterativeLinearSystem(matrix, rhs, previous_iteration_potentials_v,
                                                            linear_max_iterations, linear_tolerance,
                                                            linear_relaxation);
            solved_potentials_v = std::move(linear_result.solution);
            last_linear_iterations = linear_result.iterations;
            last_linear_residual_norm = linear_result.residual_norm;
            last_linear_converged = linear_result.converged;
        }
        else
        {
            solved_potentials_v = solveDenseLinearSystem(matrix, rhs);
            last_linear_iterations = solved_potentials_v.empty() ? 0 : 1;
            last_linear_residual_norm = denseLinearResidualNorm(matrix, rhs, solved_potentials_v);
            last_linear_converged = !solved_potentials_v.empty();
        }
        if (solved_potentials_v.size() != cell_count)
        {
            break;
        }

        for (std::size_t i = 0; i < cell_count; ++i)
        {
            auto& cell = solution.cells[i];
            cell.solved_potential_v =
                std::isfinite(solved_potentials_v[i]) ? solved_potentials_v[i]
                                                      : cell.reference_potential_v;
        }

        double max_delta_v = 0.0;
        for (std::size_t i = 0; i < cell_count; ++i)
        {
            max_delta_v = std::max(max_delta_v,
                                   std::abs(solution.cells[i].solved_potential_v -
                                            previous_iteration_potentials_v[i]));
        }
        final_max_delta_v = max_delta_v;
        previous_iteration_potentials_v = solved_potentials_v;
        solution.self_consistent_iterations_completed = iteration + 1;

        if (iteration + 1 < self_consistent_iterations)
        {
            for (std::size_t i = 0; i < cell_count; ++i)
            {
                const auto& cell = solution.cells[i];
                double lhs = surface_coefficients[i] * (cell.solved_potential_v - cell.surface_potential_v) +
                             reference_coefficients[i] *
                                 (cell.solved_potential_v - cell.reference_potential_v) +
                             temporal_coefficients[i] *
                                 (cell.solved_potential_v - previous_volume_potentials[i]);
                for (const auto& edge : solution.edges)
                {
                    if (edge.from == i)
                    {
                        lhs += edge.coefficient_m *
                               (cell.solved_potential_v - solution.cells[edge.to].solved_potential_v);
                    }
                    else if (edge.to == i)
                    {
                        lhs += edge.coefficient_m *
                               (cell.solved_potential_v - solution.cells[edge.from].solved_potential_v);
                    }
                }
                const double reconstructed_charge_density =
                    (-kEpsilon0 * lhs / std::max(1.0e-18, cell.pseudo_volume_m3)) /
                    std::max(1.0, cell.volume_to_surface_weight);
                current_charge_density_c_per_m3[i] =
                    (1.0 - charge_relaxation) * current_charge_density_c_per_m3[i] +
                    charge_relaxation * reconstructed_charge_density;
                previous_volume_potentials[i] =
                    (1.0 - potential_relaxation) * previous_volume_potentials[i] +
                    potential_relaxation * cell.solved_potential_v;
            }
        }
        if (max_delta_v <= convergence_tolerance_v)
        {
            solution.converged = true;
            break;
        }
    }
    solution.linear_iterations_completed = last_linear_iterations;
    solution.linear_residual_norm = last_linear_residual_norm;
    solution.max_delta_v = final_max_delta_v;
    if (!solution.converged)
    {
        solution.converged = last_linear_converged && final_max_delta_v <= convergence_tolerance_v;
    }

    std::vector<double> edge_flux_accumulator(cell_count, 0.0);
    for (const auto& edge : solution.edges)
    {
        const double flux_term =
            edge.coefficient_m *
            (solution.cells[edge.from].solved_potential_v - solution.cells[edge.to].solved_potential_v);
        edge_flux_accumulator[edge.from] += std::abs(flux_term);
        edge_flux_accumulator[edge.to] += std::abs(flux_term);
    }

    for (std::size_t i = 0; i < cell_count; ++i)
    {
        auto& cell = solution.cells[i];
        const double center_offset_m =
            std::max(1.0e-6, std::max(0.5 * cell.sheath_length_m,
                                      std::abs(dot3(cell.center_m, cell.surface_normal))));
        const double solved_normal_field_v_per_m =
            (cell.surface_potential_v - cell.solved_potential_v) / center_offset_m;
        cell.normal_field_v_per_m = solved_normal_field_v_per_m;

        double lhs = surface_coefficients[i] * (cell.solved_potential_v - cell.surface_potential_v) +
                     reference_coefficients[i] *
                         (cell.solved_potential_v - cell.reference_potential_v) +
                     temporal_coefficients[i] *
                         (cell.solved_potential_v - previous_volume_potentials[i]);
        for (const auto& edge : solution.edges)
        {
            if (edge.from == i)
            {
                lhs += edge.coefficient_m *
                       (cell.solved_potential_v - solution.cells[edge.to].solved_potential_v);
            }
            else if (edge.to == i)
            {
                lhs += edge.coefficient_m *
                       (cell.solved_potential_v - solution.cells[edge.from].solved_potential_v);
            }
        }
        cell.reconstructed_charge_density_c_per_m3 =
            (-kEpsilon0 * lhs / std::max(1.0e-18, cell.pseudo_volume_m3)) /
            std::max(1.0, cell.volume_to_surface_weight);
        cell.poisson_residual_v_m = std::abs(lhs + cell.deposited_charge_density_c_per_m3 *
                                                       cell.pseudo_volume_m3 / kEpsilon0);

        const double reference_span_v =
            std::max(1.0, std::abs(cell.surface_potential_v - cell.reference_potential_v));
        const double cell_span_v =
            std::abs(cell.solved_potential_v - cell.reference_potential_v);
        const double lateral_flux_metric =
            edge_flux_accumulator[i] /
            std::max(1.0e-12, surface_coefficients[i] * reference_span_v);
        cell.coupling_gain = std::clamp(
            0.65 * (cell_span_v / reference_span_v) + 0.35 * std::min(1.0, lateral_flux_metric),
            0.0, 1.0);
    }

    return solution;
}

SurfaceChargingConfig normalizeTopologyConfig(const SurfaceChargingConfig& input)
{
    SurfaceChargingConfig normalized = input;
    if (!normalized.bodies.empty() || !normalized.patches.empty() || !normalized.interfaces.empty())
    {
        normalized.surface_nodes.clear();
        normalized.surface_branches.clear();
        normalized.patch_physics_overrides.clear();

        std::unordered_map<std::string, std::size_t> node_indices;
        double total_patch_area_m2 = 0.0;

        for (const auto& body : normalized.bodies)
        {
            SurfaceCircuitNodeConfig node;
            node.name = bodyNodeName(body.id);
            node.area_m2 = body.area_m2;
            node.is_patch = false;
            node.initial_potential_v = body.initial_potential_v;
            node.capacitance_f = body.capacitance_f;
            node.fixed_potential = body.fixed_potential || !body.floating;
            node.fixed_value_v = body.fixed_potential ? body.fixed_value_v : body.initial_potential_v;
            node_indices[body.id] = normalized.surface_nodes.size();
            normalized.surface_nodes.push_back(node);
        }

        for (const auto& patch : normalized.patches)
        {
            SurfaceCircuitNodeConfig node;
            node.name = patchNodeName(patch.id);
            node.area_m2 = patch.area_m2;
            node.is_patch = true;
            node.initial_potential_v = patch.initial_potential_v;
            node.capacitance_f = patch.capacitance_f;
            node.fixed_potential = false;
            node.fixed_value_v = 0.0;
            node_indices[patch.id] = normalized.surface_nodes.size();
            normalized.surface_nodes.push_back(node);
            total_patch_area_m2 += std::max(0.0, patch.area_m2);

            SurfacePatchPhysicsConfig physics;
            physics.node_name = node.name;
            if (patch.material.has_value())
            {
                physics.override_material = true;
                physics.material = *patch.material;
            }
            if (patch.reference_see_model.has_value())
            {
                physics.override_see_model = true;
                physics.reference_see_model = *patch.reference_see_model;
            }
            if (patch.patch_incidence_angle_deg.has_value())
            {
                physics.override_patch_incidence_angle = true;
                physics.patch_incidence_angle_deg = *patch.patch_incidence_angle_deg;
            }
            if (patch.patch_flow_angle_deg.has_value())
            {
                physics.override_patch_flow_angle = true;
                physics.patch_flow_angle_deg = *patch.patch_flow_angle_deg;
            }
            if (patch.photoelectron_temperature_ev.has_value())
            {
                physics.override_photoelectron_temperature = true;
                physics.photoelectron_temperature_ev = *patch.photoelectron_temperature_ev;
            }
            if (patch.patch_photo_current_density_a_per_m2.has_value())
            {
                physics.override_patch_photo_current_density = true;
                physics.patch_photo_current_density_a_per_m2 =
                    *patch.patch_photo_current_density_a_per_m2;
            }
            if (patch.electron_collection_coefficient.has_value())
            {
                physics.override_electron_collection_coefficient = true;
                physics.electron_collection_coefficient = *patch.electron_collection_coefficient;
            }
            if (patch.ion_collection_coefficient.has_value())
            {
                physics.override_ion_collection_coefficient = true;
                physics.ion_collection_coefficient = *patch.ion_collection_coefficient;
            }
            if (patch.plasma.has_value())
            {
                physics.override_plasma = true;
                physics.plasma = *patch.plasma;
            }
            if (patch.electron_spectrum.has_value())
            {
                physics.override_electron_spectrum = true;
                physics.electron_spectrum = *patch.electron_spectrum;
                physics.has_electron_spectrum = patch.has_electron_spectrum.value_or(true);
            }
            if (patch.ion_spectrum.has_value())
            {
                physics.override_ion_spectrum = true;
                physics.ion_spectrum = *patch.ion_spectrum;
                physics.has_ion_spectrum = patch.has_ion_spectrum.value_or(true);
            }
            if (patch.emission.has_value())
            {
                physics.override_emission = true;
                physics.emission = *patch.emission;
            }
            normalized.patch_physics_overrides.push_back(physics);

            if (!patch.body_id.empty() && node_indices.contains(patch.body_id))
            {
                SurfaceCircuitBranchConfig branch;
                branch.from_node = node_indices[patch.id];
                branch.to_node = node_indices[patch.body_id];
                branch.conductance_s = 0.0;
                branch.resistance_ohm = normalized.patch_body_resistance_ohm;
                branch.bias_v = 0.0;
                normalized.surface_branches.push_back(branch);
            }
        }

        for (const auto& interface_config : normalized.interfaces)
        {
            if (!node_indices.contains(interface_config.from_id) ||
                !node_indices.contains(interface_config.to_id))
            {
                continue;
            }
            SurfaceCircuitBranchConfig branch;
            branch.from_node = node_indices[interface_config.from_id];
            branch.to_node = node_indices[interface_config.to_id];
            branch.conductance_s = interface_config.conductance_s;
            branch.resistance_ohm = interface_config.resistance_ohm;
            branch.bias_v = interface_config.bias_v;
            normalized.surface_branches.push_back(branch);
        }

        if (!normalized.bodies.empty())
        {
            const auto& primary_body = normalized.bodies.front();
            normalized.body_initial_potential_v = primary_body.initial_potential_v;
            normalized.body_capacitance_f = primary_body.capacitance_f;
            normalized.body_leakage_conductance_s = primary_body.leakage_conductance_s;
            normalized.body_floating = primary_body.floating && !primary_body.fixed_potential;
        }
        if (!normalized.patches.empty())
        {
            const auto& primary_patch = normalized.patches.front();
            normalized.surface_area_m2 = total_patch_area_m2 > 0.0 ? total_patch_area_m2
                                                                   : primary_patch.area_m2;
            normalized.default_surface_physics.material = normalized.material;
            normalized.default_surface_physics.reference_see_model = normalized.reference_see_model;
            normalized.default_surface_physics.patch_incidence_angle_deg =
                normalized.patch_incidence_angle_deg;
            normalized.default_surface_physics.patch_flow_angle_deg =
                normalized.patch_flow_angle_deg;
            normalized.default_surface_physics.photoelectron_temperature_ev =
                normalized.photoelectron_temperature_ev;
            normalized.default_surface_physics.patch_photo_current_density_a_per_m2 =
                normalized.patch_photo_current_density_a_per_m2;
            normalized.default_surface_physics.electron_collection_coefficient =
                normalized.electron_collection_coefficient;
            normalized.default_surface_physics.ion_collection_coefficient =
                normalized.ion_collection_coefficient;
            normalized.default_surface_physics.plasma = normalized.plasma;
            normalized.default_surface_physics.electron_spectrum = normalized.electron_spectrum;
            normalized.default_surface_physics.ion_spectrum = normalized.ion_spectrum;
            normalized.default_surface_physics.has_electron_spectrum =
                normalized.has_electron_spectrum;
            normalized.default_surface_physics.has_ion_spectrum =
                normalized.has_ion_spectrum;
            normalized.default_surface_physics.emission = normalized.emission;
        }
    }
    else if (normalized.default_surface_physics.material.getName().empty())
    {
        normalized.default_surface_physics.material = normalized.material;
    }

    return normalized;
}

const SurfacePatchPhysicsConfig* findPatchPhysicsOverride(const SurfaceChargingConfig& config,
                                                          const SurfaceModelRuntimeState& state)
{
    for (const auto& override_config : config.patch_physics_overrides)
    {
        if (override_config.match_by_index && override_config.node_index == state.node_index)
        {
            return &override_config;
        }
        if (!override_config.node_name.empty() && override_config.node_name == state.node_name)
        {
            return &override_config;
        }
    }
    return nullptr;
}

SurfaceChargingConfig effectivePatchConfig(const SurfaceChargingConfig& base_config,
                                           const SurfaceModelRuntimeState& state)
{
    SurfaceChargingConfig effective = base_config;
    const auto* override_config = findPatchPhysicsOverride(base_config, state);
    if (override_config == nullptr)
    {
        return effective;
    }

    if (override_config->override_material)
    {
        effective.material = override_config->material;
    }
    if (override_config->override_see_model)
    {
        effective.reference_see_model = override_config->reference_see_model;
    }
    if (override_config->override_patch_incidence_angle)
    {
        effective.patch_incidence_angle_deg = override_config->patch_incidence_angle_deg;
    }
    if (override_config->override_patch_flow_angle)
    {
        effective.patch_flow_angle_deg = override_config->patch_flow_angle_deg;
    }
    if (override_config->override_photoelectron_temperature)
    {
        effective.photoelectron_temperature_ev = override_config->photoelectron_temperature_ev;
    }
    if (override_config->override_patch_photo_current_density)
    {
        effective.patch_photo_current_density_a_per_m2 =
            override_config->patch_photo_current_density_a_per_m2;
    }
    if (override_config->override_electron_collection_coefficient)
    {
        effective.electron_collection_coefficient =
            override_config->electron_collection_coefficient;
    }
    if (override_config->override_ion_collection_coefficient)
    {
        effective.ion_collection_coefficient = override_config->ion_collection_coefficient;
    }
    if (override_config->override_plasma)
    {
        effective.plasma = override_config->plasma;
    }
    if (override_config->override_electron_spectrum)
    {
        effective.electron_spectrum = override_config->electron_spectrum;
        effective.has_electron_spectrum = override_config->has_electron_spectrum;
    }
    if (override_config->override_ion_spectrum)
    {
        effective.ion_spectrum = override_config->ion_spectrum;
        effective.has_ion_spectrum = override_config->has_ion_spectrum;
    }
    if (override_config->override_emission)
    {
        effective.emission = override_config->emission;
    }
    effective.surface_area_m2 = state.surface_area_m2 > 0.0 ? state.surface_area_m2
                                                             : base_config.surface_area_m2;
    return effective;
}

double safeExp(double exponent)
{
    return std::exp(std::clamp(exponent, -200.0, 50.0));
}

double emittedElectronEscapeProbability(double surface_potential_v, double characteristic_energy_ev)
{
    if (surface_potential_v <= 0.0)
    {
        return 1.0;
    }

    const double escape_energy = std::max(1.0e-3, characteristic_energy_ev);
    return std::clamp(safeExp(-surface_potential_v / escape_energy), 0.0, 1.0);
}

double computePatchFieldEmission(const SurfaceChargingConfig& config,
                                 const SurfaceModelRuntimeState& state)
{
    const double work_function = std::max(0.1, config.material.getWorkFunctionEv());
    const double sheath_length = std::max(1.0e-6, state.effective_sheath_length_m);
    const double reference_field_v_per_m =
        std::abs(state.reference_plasma_potential_v - state.patch_potential_v) / sheath_length;
    const double local_field_v_per_m =
        std::max(std::abs(state.normal_electric_field_v_per_m), reference_field_v_per_m);
    if (local_field_v_per_m <= 1.0e7)
    {
        return 0.0;
    }

    const double fn_prefactor = 1.54e-6 * local_field_v_per_m * local_field_v_per_m / work_function;
    const double fn_exponent =
        -6.83e9 * std::pow(work_function, 1.5) / std::max(1.0e6, local_field_v_per_m);
    return fn_prefactor * safeExp(fn_exponent);
}

double computePatchThermionicEmission(const SurfaceChargingConfig& config,
                                      const SurfaceModelRuntimeState& state)
{
    const double work_function = std::max(0.1, config.material.getWorkFunctionEv());
    const double thermal_energy_ev =
        std::max(1.0e-6, kBoltzmannEvPerK * config.emission.surface_temperature_k);
    return kRichardsonConstant * config.emission.surface_temperature_k *
           config.emission.surface_temperature_k * safeExp(-work_function / thermal_energy_ev) *
           emittedElectronEscapeProbability(state.patch_potential_v,
                                            std::max(1.0e-3, thermal_energy_ev));
}

double legacyErfApprox(double x)
{
    return std::erf(x);
}

double legacyThermalCurrentDensity(double density_m3, double temperature_ev, double mass_kg)
{
    return kElementaryCharge * std::max(0.0, density_m3) *
           std::sqrt(std::max(1.0e-6, temperature_ev) * kElementaryCharge /
                     (2.0 * kPi * std::max(1.0e-30, mass_kg)));
}

double legacyWhippleYield(double peak_energy_ev, double peak_yield, double incident_energy_ev)
{
    if (incident_energy_ev <= 10.0)
    {
        return 0.0;
    }
    const double q = 2.28 * std::pow(incident_energy_ev / std::max(1.0, peak_energy_ev), 1.35);
    if (q <= 1.0e-12 || incident_energy_ev <= 0.0)
    {
        return 0.0;
    }
    return 2.228 * (q - 1.0 + std::exp(-q)) / q * peak_yield *
           std::pow(std::max(1.0, peak_energy_ev) / incident_energy_ev, 0.35);
}

double legacySimsYield(double peak_energy_ev, double peak_yield, double exponent_n,
                       double incident_energy_ev)
{
    if (incident_energy_ev <= 10.0)
    {
        return 0.0;
    }
    const double n = std::max(1.05, exponent_n);
    double xm = 2.5;
    for (int i = 1; i < 10001; ++i)
    {
        const double x = 0.5 + static_cast<double>(i) * (2.5 - 0.5) / 10000.0;
        const double y = (1.0 - 1.0 / n) * (std::exp(x) - 1.0);
        if (std::abs(x - y) <= 1.0e-3)
        {
            xm = x;
            break;
        }
    }
    const double x = xm * std::pow(incident_energy_ev / std::max(1.0, peak_energy_ev), n);
    if (x <= 1.0e-12 || std::abs(1.0 - std::exp(-xm)) < 1.0e-12)
    {
        return 0.0;
    }
    return peak_yield / (1.0 - std::exp(-xm)) *
           std::pow(std::max(1.0, peak_energy_ev) / incident_energy_ev, n - 1.0) * 2.0 *
           (x + std::exp(-x) - 1.0) / x;
}

double legacyKatzYield(const Material::MaterialProperty& material, double incident_energy_ev)
{
    if (incident_energy_ev <= 1.0e-9)
    {
        return 0.0;
    }
    const double peak_energy_ev =
        std::max(1.0, material.getScalarProperty("secondary_yield_peak_energy_ev", 400.0));
    const double peak_yield = std::max(0.0, material.getSecondaryElectronYield());
    const double r1 = material.getScalarProperty("katz_r1", 0.85);
    const double n1 = material.getScalarProperty("katz_n1", 0.35);
    const double r2 = material.getScalarProperty("katz_r2", 0.12);
    const double n2 = material.getScalarProperty("katz_n2", 1.05);
    const double energy_kev = incident_energy_ev * 1.0e-3;
    const double peak_energy_kev = peak_energy_ev * 1.0e-3;
    const double range = r1 * std::pow(energy_kev, n1) + r2 * std::pow(energy_kev, n2);
    const double peak_range =
        r1 * std::pow(peak_energy_kev, n1) + r2 * std::pow(peak_energy_kev, n2);
    if (range <= 1.0e-12 || peak_range <= 1.0e-12)
    {
        return 0.0;
    }
    const double attenuation = std::exp(-std::pow(range / peak_range - 1.0, 2.0));
    const double rolloff =
        1.0 / (1.0 + std::pow(incident_energy_ev / std::max(1.0, 3.0 * peak_energy_ev), 0.7));
    return peak_yield * attenuation * (0.55 + 0.45 * rolloff);
}

double legacySeeYield(const SurfaceChargingConfig& config, const Material::MaterialProperty& material,
                      double incident_energy_ev)
{
    const double peak_energy =
        std::max(1.0, material.getScalarProperty("secondary_yield_peak_energy_ev", 300.0));
    const double peak_yield = std::max(0.0, material.getSecondaryElectronYield());
    switch (config.reference_see_model)
    {
    case SecondaryElectronEmissionModel::Whipple:
        return legacyWhippleYield(peak_energy, peak_yield, incident_energy_ev);
    case SecondaryElectronEmissionModel::Sims:
        return legacySimsYield(peak_energy, peak_yield,
                               material.getScalarProperty("sims_exponent_n", 1.6),
                               incident_energy_ev);
    case SecondaryElectronEmissionModel::Katz:
        return legacyKatzYield(material, incident_energy_ev);
    }
    return 0.0;
}

double legacyIonSecondaryYield(const Material::MaterialProperty& material, double incident_energy_ev)
{
    const double peak_energy_kev =
        std::max(1.0e-6, material.getScalarProperty("ion_secondary_peak_energy_kev", 0.35));
    const double energy_kev = std::max(0.0, incident_energy_ev * 1.0e-3);
    const double peak_yield =
        std::max(0.0, material.getScalarProperty("ion_secondary_yield", 0.08));
    if (energy_kev <= 0.0)
    {
        return 0.0;
    }
    if (energy_kev < 0.476)
    {
        return 0.5 * std::pow(2.0, -2.0) * peak_yield * std::sqrt(energy_kev) /
               (1.0 + energy_kev / peak_energy_kev);
    }
    if (energy_kev <= 10.0)
    {
        const double exponent = -(1.0 / energy_kev - 0.1);
        return (2.0 - (1.0 / energy_kev - 0.1) / 2.0) * std::pow(2.0, exponent) * peak_yield *
               std::sqrt(energy_kev) / (1.0 + energy_kev / peak_energy_kev);
    }
    return 2.0 * peak_yield * std::sqrt(energy_kev) / (1.0 + energy_kev / peak_energy_kev);
}

double legacyBackscatterYield(const Material::MaterialProperty& material, double incident_energy_ev)
{
    const double z = material.getScalarProperty("atomic_number", 13.0);
    const double es = incident_energy_ev;
    double ybe = 0.0;
    if (es < 50.0)
    {
        ybe = 0.0;
    }
    else if (es < 1.0e3)
    {
        ybe = 0.3338 * std::log(es / 50.0) *
              (1.0 - std::pow(0.7358, 0.037 * z) + 0.1 * std::exp(-es / 5000.0));
    }
    else if (es < 1.0e4)
    {
        ybe = 1.0 - std::pow(0.7358, 0.037 * z) + 0.1 * std::exp(-es / 5000.0);
    }
    else if (es < 1.0e5)
    {
        ybe = 1.0 - std::pow(0.7358, 0.037 * z);
    }
    if (ybe <= 0.0 || std::abs(std::log(ybe)) < 1.0e-12)
    {
        return 0.0;
    }
    return 2.0 * (1.0 - ybe + ybe * std::log(ybe)) / (std::log(ybe) * std::log(ybe));
}

double legacyRamPatchCurrentDensity(double surface_potential_v, double ion_density_m3,
                                    double ion_temperature_ev, double ion_mass_amu,
                                    double flow_speed_m_per_s, double patch_flow_angle_deg)
{
    const double ni_cm3 = std::max(0.0, ion_density_m3 * 1.0e-6);
    const double ti_ev = std::max(1.0e-6, ion_temperature_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / std::max(1.0, ion_mass_amu));
    const double alpha = patch_flow_angle_deg * kPi / 180.0;
    const double qd = std::max(0.0, flow_speed_m_per_s * 1.0e-3 * std::cos(alpha) / std::max(1.0e-12, wt));
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    double current_na_per_m2 = 0.08011 * ni_cm3 * wt;
    if (surface_potential_v >= 0.0)
    {
        current_na_per_m2 *=
            0.5642 * std::exp(-(qd - qv) * (qd - qv)) + qd + qd * legacyErfApprox(qd - qv);
    }
    else
    {
        current_na_per_m2 *= 0.5642 * std::exp(-qd * qd) + qd + qd * legacyErfApprox(qd);
    }
    return std::max(0.0, current_na_per_m2) * 1.0e-9;
}

double legacyRamBodyCurrentDensity(double surface_potential_v, double ion_density_m3,
                                   double ion_temperature_ev, double ion_mass_amu,
                                   double flow_speed_m_per_s)
{
    const double ni_cm3 = std::max(0.0, ion_density_m3 * 1.0e-6);
    const double ti_ev = std::max(1.0e-6, ion_temperature_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / std::max(1.0, ion_mass_amu));
    const double qd = std::max(1.0e-12, flow_speed_m_per_s * 1.0e-3 / std::max(1.0e-12, wt));
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    double current_na_per_m2 = 0.0;
    if (surface_potential_v >= 0.0)
    {
        current_na_per_m2 =
            0.5642 / qd *
                ((qv + qd) * std::exp(-(qv - qd) * (qv - qd)) -
                 (qv - qd) * std::exp(-(qv + qd) * (qv + qd))) +
            (0.5 / qd + qd - qv * qv / qd) *
                (legacyErfApprox(qv + qd) - legacyErfApprox(qv - qd));
        current_na_per_m2 *= 0.02003 * ni_cm3 * wt;
    }
    else
    {
        current_na_per_m2 = 0.04005 * ni_cm3 * wt *
                            (0.5642 * std::exp(-qd * qd) +
                             (qd + 0.5 / qd) * legacyErfApprox(qd));
        current_na_per_m2 *= 4.0;
    }
    return std::max(0.0, current_na_per_m2) * 1.0e-9;
}

struct LegacyPopulation
{
    double density_m3 = 0.0;
    double temperature_ev = 0.0;
    double mass_amu = 1.0;
};

struct LegacyLeoFluxTables
{
    std::vector<double> energy_ev;
    std::vector<double> width_ev;
    std::vector<double> spectrum_energy_grid_ev;
    std::vector<double> smoothed_electron_flux;
    std::vector<double> smoothed_ion_flux;
    std::vector<double> total_electron_flux;
    std::vector<double> spectrum_electron_flux;
    std::vector<double> total_ion_flux;
    std::vector<double> spectrum_ion_flux;
    std::vector<std::vector<double>> electron_component_flux;
    std::vector<std::vector<double>> ion_component_flux;
    std::vector<LegacyPopulation> electron_populations;
    std::vector<LegacyPopulation> ion_populations;
    double average_spectrum_electron_temperature_ev = 0.0;
    double average_spectrum_ion_temperature_ev = 0.0;
};

std::vector<LegacyPopulation> collectLegacyElectronPopulations(const SurfaceChargingConfig& config)
{
    std::vector<LegacyPopulation> populations;
    if (config.has_electron_spectrum && !config.electron_spectrum.populations.empty())
    {
        for (const auto& population : config.electron_spectrum.populations)
        {
            populations.push_back(
                {population.density_m3, population.temperature_ev, population.mass_amu});
            if (populations.size() >= 3)
            {
                break;
            }
        }
    }
    if (populations.empty())
    {
        populations.push_back(
            {config.plasma.electron_density_m3, config.plasma.electron_temperature_ev,
             9.1093837015e-31 / 1.66053906660e-27});
    }
    return populations;
}

std::vector<LegacyPopulation> collectLegacyIonPopulations(const SurfaceChargingConfig& config)
{
    std::vector<LegacyPopulation> populations;
    if (config.has_ion_spectrum && !config.ion_spectrum.populations.empty())
    {
        for (const auto& population : config.ion_spectrum.populations)
        {
            populations.push_back(
                {population.density_m3, population.temperature_ev, population.mass_amu});
            if (populations.size() >= 3)
            {
                break;
            }
        }
    }
    if (populations.empty())
    {
        populations.push_back(
            {config.plasma.ion_density_m3, config.plasma.ion_temperature_ev, config.plasma.ion_mass_amu});
    }
    return populations;
}

double legacyInterpolateSpectrum(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 double target)
{
    if (x.empty() || y.empty())
    {
        return 0.0;
    }
    const std::size_t count = std::min(x.size(), y.size());
    if (count == 1)
    {
        return y.front();
    }
    if (count == 2)
    {
        const double x0 = x[0];
        const double x1 = x[1];
        if (x1 <= x0)
        {
            return y[0];
        }
        return (y[0] * (target - x1) - y[1] * (target - x0)) / (x0 - x1);
    }

    std::ptrdiff_t k = 0;
    std::ptrdiff_t m = 0;
    if (target <= x[1])
    {
        k = 0;
        m = 2;
    }
    else if (target >= x[count - 2])
    {
        k = static_cast<std::ptrdiff_t>(count) - 3;
        m = static_cast<std::ptrdiff_t>(count) - 1;
    }
    else
    {
        k = 1;
        m = static_cast<std::ptrdiff_t>(count);
        while (m - k != 1)
        {
            const auto i = (k + m) / 2;
            if (target < x[static_cast<std::size_t>(i - 1)])
            {
                m = i;
            }
            else
            {
                k = i;
            }
        }
        k -= 1;
        m -= 1;
        if (std::abs(target - x[static_cast<std::size_t>(k)]) <
            std::abs(target - x[static_cast<std::size_t>(m)]))
        {
            k -= 1;
        }
        else
        {
            m += 1;
        }
    }

    double value = 0.0;
    for (auto i = k; i <= m; ++i)
    {
        double weight = 1.0;
        for (auto j = k; j <= m; ++j)
        {
            if (j != i)
            {
                weight *=
                    (target - x[static_cast<std::size_t>(j)]) /
                    (x[static_cast<std::size_t>(i)] - x[static_cast<std::size_t>(j)]);
            }
        }
        value += weight * y[static_cast<std::size_t>(i)];
    }
    return value;
}

std::vector<double> legacySmoothLine3(const std::vector<double>& x, const std::vector<double>& y)
{
    const std::size_t count = std::min(x.size(), y.size());
    if (count <= 2)
    {
        return std::vector<double>(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(count));
    }

    std::vector<double> smoothed(count, 0.0);
    for (std::size_t i = 0; i < count; ++i)
    {
        std::size_t first = 0;
        if (i == 0)
        {
            first = 0;
        }
        else if (i >= count - 1)
        {
            first = count - 3;
        }
        else
        {
            first = i - 1;
        }
        const std::size_t last = std::min(count - 1, first + 2);
        const double x0 = x[first];
        const double x1 = x[first + 1];
        const double x2 = x[last];
        const double y0 = y[first];
        const double y1 = y[first + 1];
        const double y2 = y[last];
        const double t0 = 3.0;
        const double t1 = x0 + x1 + x2;
        const double t2 = x0 * x0 + x1 * x1 + x2 * x2;
        const double p1 = y0 + y1 + y2;
        const double p2 = y0 * x0 + y1 * x1 + y2 * x2;
        const double det = t0 * t2 - t1 * t1;
        if (std::abs(det) <= 1.0e-20)
        {
            smoothed[i] = y[i];
            continue;
        }
        const double a0 = (p1 * t2 - p2 * t1) / det;
        const double a1 = (t0 * p2 - t1 * p1) / det;
        smoothed[i] = a0 + a1 * x[i];
    }
    return smoothed;
}

double legacyElectronComponentFlux(double density_m3, double temperature_ev, double energy_ev)
{
    const double density_cm3 = std::max(0.0, density_m3) * 1.0e-6;
    const double temperature = std::max(1.0e-6, temperature_ev);
    return 5.325e10 * density_cm3 * energy_ev * std::pow(temperature, -1.5) *
           std::exp(-energy_ev / temperature);
}

double legacyIonComponentFlux(double density_m3, double temperature_ev, double mass_amu,
                              double energy_ev)
{
    const double density_cm3 = std::max(0.0, density_m3) * 1.0e-6;
    const double temperature = std::max(1.0e-6, temperature_ev);
    return 1.244e9 * density_cm3 * energy_ev * std::pow(temperature, -1.5) *
           std::exp(-energy_ev / temperature) / std::sqrt(std::max(1.0, mass_amu));
}

double legacyThermalIonCollectionNaPerM2(const std::vector<LegacyPopulation>& populations,
                                         double surface_potential_v,
                                         bool include_potential_barrier)
{
    double current_na_per_m2 = 0.0;
    for (const auto& population : populations)
    {
        if (population.density_m3 <= 0.0)
        {
            continue;
        }
        const double ni_cm3 = population.density_m3 * 1.0e-6;
        const double temperature = std::max(1.0e-6, population.temperature_ev);
        double term = 0.6255 * ni_cm3 * std::sqrt(temperature / std::max(1.0, population.mass_amu));
        if (include_potential_barrier && surface_potential_v > 0.0)
        {
            term *= std::exp(-surface_potential_v / temperature);
        }
        current_na_per_m2 += term;
    }
    return current_na_per_m2;
}

LegacyLeoFluxTables buildLegacyLeoFluxTables(const SurfaceChargingConfig& config)
{
    constexpr int kLegacyLeoSamples = 1000;
    LegacyLeoFluxTables tables;
    tables.electron_populations = collectLegacyElectronPopulations(config);
    tables.ion_populations = collectLegacyIonPopulations(config);

    double spectrum_min_ev = std::numeric_limits<double>::max();
    double spectrum_max_ev = 0.0;
    auto extend_spectrum_range = [&](const std::vector<double>& energies) {
        for (double energy_ev : energies)
        {
            if (energy_ev > 0.0)
            {
                spectrum_min_ev = std::min(spectrum_min_ev, energy_ev);
                spectrum_max_ev = std::max(spectrum_max_ev, energy_ev);
            }
        }
    };
    extend_spectrum_range(config.electron_spectrum.energy_grid_ev);
    extend_spectrum_range(config.ion_spectrum.energy_grid_ev);

    double min_thermal_ev = std::numeric_limits<double>::max();
    double max_thermal_ev = 0.0;
    auto extend_thermal_range = [&](const std::vector<LegacyPopulation>& populations) {
        for (const auto& population : populations)
        {
            if (population.temperature_ev > 0.0)
            {
                min_thermal_ev = std::min(min_thermal_ev, population.temperature_ev);
                max_thermal_ev = std::max(max_thermal_ev, population.temperature_ev);
            }
        }
    };
    extend_thermal_range(tables.electron_populations);
    extend_thermal_range(tables.ion_populations);

    double minimum_energy_ev = std::isfinite(spectrum_min_ev) ? spectrum_min_ev : 0.0;
    double maximum_energy_ev = spectrum_max_ev;
    if (std::isfinite(min_thermal_ev))
    {
        minimum_energy_ev = minimum_energy_ev > 0.0
                                ? std::min(minimum_energy_ev, 0.05 * min_thermal_ev)
                                : 0.05 * min_thermal_ev;
    }
    if (max_thermal_ev > 0.0)
    {
        maximum_energy_ev = std::max(maximum_energy_ev, 10.0 * max_thermal_ev);
    }
    minimum_energy_ev = std::max(1.0e-4, minimum_energy_ev > 0.0 ? minimum_energy_ev : 1.0e-2);
    maximum_energy_ev = std::max(10.0 * minimum_energy_ev, maximum_energy_ev);

    const double nmin = std::log10(minimum_energy_ev);
    const double nmax = std::log10(maximum_energy_ev);
    const double dn = (nmax - nmin) / static_cast<double>(kLegacyLeoSamples);

    tables.energy_ev.reserve(kLegacyLeoSamples);
    tables.width_ev.reserve(kLegacyLeoSamples);
    tables.total_electron_flux.reserve(kLegacyLeoSamples);
    tables.spectrum_electron_flux.reserve(kLegacyLeoSamples);
    tables.total_ion_flux.reserve(kLegacyLeoSamples);
    tables.spectrum_ion_flux.reserve(kLegacyLeoSamples);
    tables.electron_component_flux.assign(tables.electron_populations.size(),
                                          std::vector<double>{});
    tables.ion_component_flux.assign(tables.ion_populations.size(), std::vector<double>{});

    double spectrum_electron_energy_moment = 0.0;
    double spectrum_electron_flux_sum = 0.0;
    double spectrum_ion_energy_moment = 0.0;
    double spectrum_ion_flux_sum = 0.0;

    const bool has_electron_spectrum = config.has_electron_spectrum &&
                                       !config.electron_spectrum.energy_grid_ev.empty() &&
                                       !config.electron_spectrum.differential_number_flux.empty();
    const bool has_ion_spectrum = config.has_ion_spectrum &&
                                  !config.ion_spectrum.energy_grid_ev.empty() &&
                                  !config.ion_spectrum.differential_number_flux.empty();
    const double electron_spectrum_min_ev =
        has_electron_spectrum
            ? *std::min_element(config.electron_spectrum.energy_grid_ev.begin(),
                                config.electron_spectrum.energy_grid_ev.end())
            : 0.0;
    const double electron_spectrum_max_ev =
        has_electron_spectrum
            ? *std::max_element(config.electron_spectrum.energy_grid_ev.begin(),
                                config.electron_spectrum.energy_grid_ev.end())
            : 0.0;
    const double ion_spectrum_min_ev =
        has_ion_spectrum
            ? *std::min_element(config.ion_spectrum.energy_grid_ev.begin(),
                                config.ion_spectrum.energy_grid_ev.end())
            : 0.0;
    const double ion_spectrum_max_ev =
        has_ion_spectrum
            ? *std::max_element(config.ion_spectrum.energy_grid_ev.begin(),
                                config.ion_spectrum.energy_grid_ev.end())
            : 0.0;
    const std::vector<double> smoothed_electron_flux =
        has_electron_spectrum
            ? legacySmoothLine3(config.electron_spectrum.energy_grid_ev,
                                config.electron_spectrum.differential_number_flux)
            : std::vector<double>{};
    const std::vector<double> smoothed_ion_flux =
        has_ion_spectrum
            ? legacySmoothLine3(config.ion_spectrum.energy_grid_ev,
                                config.ion_spectrum.differential_number_flux)
            : std::vector<double>{};
    tables.spectrum_energy_grid_ev = config.electron_spectrum.energy_grid_ev;
    tables.smoothed_electron_flux = smoothed_electron_flux;
    tables.smoothed_ion_flux = smoothed_ion_flux;

    for (int i = 0; i < kLegacyLeoSamples; ++i)
    {
        const double energy_ev = std::pow(10.0, nmin + dn * static_cast<double>(i));
        const double next_energy_ev = std::pow(10.0, nmin + dn * static_cast<double>(i + 1));
        const double width_ev = next_energy_ev - energy_ev;
        tables.energy_ev.push_back(energy_ev);
        tables.width_ev.push_back(width_ev);

        double electron_component_sum = 0.0;
        for (std::size_t j = 0; j < tables.electron_populations.size(); ++j)
        {
            const auto& population = tables.electron_populations[j];
            const double value =
                legacyElectronComponentFlux(population.density_m3, population.temperature_ev, energy_ev);
            tables.electron_component_flux[j].push_back(value);
            electron_component_sum += value;
        }

        double ion_component_sum = 0.0;
        for (std::size_t j = 0; j < tables.ion_populations.size(); ++j)
        {
            const auto& population = tables.ion_populations[j];
            const double value = legacyIonComponentFlux(population.density_m3, population.temperature_ev,
                                                        population.mass_amu, energy_ev);
            tables.ion_component_flux[j].push_back(value);
            ion_component_sum += value;
        }

        double spectrum_electron_flux = 0.0;
        if (has_electron_spectrum && energy_ev >= electron_spectrum_min_ev &&
            energy_ev <= electron_spectrum_max_ev)
        {
            spectrum_electron_flux = std::max(
                0.0, legacyInterpolateSpectrum(config.electron_spectrum.energy_grid_ev,
                                               smoothed_electron_flux, energy_ev));
        }

        double spectrum_ion_flux = 0.0;
        if (has_ion_spectrum && energy_ev >= ion_spectrum_min_ev && energy_ev <= ion_spectrum_max_ev)
        {
            spectrum_ion_flux = std::max(
                0.0,
                legacyInterpolateSpectrum(config.ion_spectrum.energy_grid_ev, smoothed_ion_flux,
                                          energy_ev));
        }

        tables.total_electron_flux.push_back(spectrum_electron_flux > 0.0 ? spectrum_electron_flux
                                                                           : electron_component_sum);
        tables.spectrum_electron_flux.push_back(spectrum_electron_flux);
        tables.total_ion_flux.push_back(spectrum_ion_flux > 0.0 ? spectrum_ion_flux : ion_component_sum);
        tables.spectrum_ion_flux.push_back(spectrum_ion_flux);

        spectrum_electron_energy_moment += spectrum_electron_flux * energy_ev * width_ev;
        spectrum_electron_flux_sum += spectrum_electron_flux * width_ev;
        spectrum_ion_energy_moment += spectrum_ion_flux * energy_ev * width_ev;
        spectrum_ion_flux_sum += spectrum_ion_flux * width_ev;
    }

    if (spectrum_electron_flux_sum > 0.0)
    {
        tables.average_spectrum_electron_temperature_ev =
            2.0 * spectrum_electron_energy_moment / (3.0 * spectrum_electron_flux_sum);
    }
    if (spectrum_ion_flux_sum > 0.0)
    {
        tables.average_spectrum_ion_temperature_ev =
            2.0 * spectrum_ion_energy_moment / (3.0 * spectrum_ion_flux_sum);
    }
    return tables;
}

SurfaceCurrents evaluateLegacyLeoBodyCurrents(const SurfaceChargingConfig& effective_config,
                                              const Material::MaterialProperty& material,
                                              double potential_v)
{
    constexpr double kTseEv = 2.0;
    constexpr double kTsiEv = 5.0;
    constexpr double kTphEv = 3.0;
    constexpr double kLegacyLeoBodySecondaryEscapeFactor = 0.17;
    constexpr double kLegacyLeoBodyThermalIonCollectionFactor = 0.25;
    constexpr double kLegacyLeoBodySpectrumDecayScale = 8.0;
    const LegacyLeoFluxTables tables = buildLegacyLeoFluxTables(effective_config);
    SurfaceCurrents currents;

    double xe = 0.0;
    double x1 = 0.0;
    double xi = 0.0;
    double x2_na_per_m2 = 0.0;
    double xse = 0.0;
    double xsi = 0.0;
    double xb = 0.0;

    if (std::abs(potential_v) <= 1.0e-12)
    {
        for (std::size_t i = 0; i < tables.energy_ev.size(); ++i)
        {
            const double energy_ev = tables.energy_ev[i];
            const double width_ev = tables.width_ev[i];
            const double total_electron_flux = tables.total_electron_flux[i];
            const double total_ion_flux = tables.total_ion_flux[i];
            const double see_at_zero = legacySeeYield(effective_config, material, energy_ev);
            const double backscatter_at_zero = legacyBackscatterYield(material, energy_ev);
            const double ion_secondary_at_zero =
                legacyIonSecondaryYield(material, energy_ev);

            xe += total_electron_flux * width_ev;
            xi += total_ion_flux * width_ev;
            xse += see_at_zero * total_electron_flux * width_ev;
            xsi += ion_secondary_at_zero * total_ion_flux * width_ev;
            xb += backscatter_at_zero * total_electron_flux * width_ev;
        }

        currents.electron_current_a_per_m2 = -kPi * kElementaryCharge * xe;
        currents.ion_current_a_per_m2 = kPi * kElementaryCharge * xi;
        currents.secondary_emission_a_per_m2 = kPi * kElementaryCharge * xse;
        currents.ion_secondary_emission_a_per_m2 = kPi * kElementaryCharge * xsi;
        currents.backscatter_emission_a_per_m2 = kPi * kElementaryCharge * xb;
        currents.photo_emission_a_per_m2 =
            0.25 * std::max(0.0, effective_config.body_photo_current_density_a_per_m2);
        currents.ram_ion_current_a_per_m2 = 0.0;
        currents.conduction_current_a_per_m2 = 0.0;
        currents.current_derivative_a_per_m2_per_v = 0.0;
        currents.total_current_a_per_m2 =
            currents.electron_current_a_per_m2 + currents.ion_current_a_per_m2 +
            currents.secondary_emission_a_per_m2 + currents.ion_secondary_emission_a_per_m2 +
            currents.backscatter_emission_a_per_m2 + currents.photo_emission_a_per_m2;
        return currents;
    }

    for (std::size_t i = 0; i < tables.energy_ev.size(); ++i)
    {
        const double energy_ev = tables.energy_ev[i];
        const double width_ev = tables.width_ev[i];
        const double total_electron_flux = tables.total_electron_flux[i];
        const double total_ion_flux = tables.total_ion_flux[i];
        const double spectrum_electron_flux = tables.spectrum_electron_flux[i];
        const double spectrum_ion_flux = tables.spectrum_ion_flux[i];
        const double see_at_v = legacySeeYield(effective_config, material,
                                               potential_v > 0.0 ? energy_ev + potential_v : energy_ev);
        const double see_at_zero = legacySeeYield(effective_config, material, energy_ev);
        const double backscatter_at_v =
            legacyBackscatterYield(material, potential_v > 0.0 ? energy_ev + potential_v : energy_ev);
        const double backscatter_at_zero = legacyBackscatterYield(material, energy_ev);
        const double ion_secondary_at_v =
            legacyIonSecondaryYield(material, potential_v > 0.0 ? energy_ev : energy_ev - potential_v);

        if (potential_v > 0.0)
        {
            xe += spectrum_electron_flux * width_ev;
            for (std::size_t j = 0; j < tables.electron_component_flux.size(); ++j)
            {
                x1 += tables.electron_component_flux[j][i] * width_ev;
            }

            if (energy_ev >= std::max(tables.energy_ev.front(), potential_v))
            {
                xi += spectrum_ion_flux * width_ev * (1.0 - potential_v / energy_ev);
                xsi += ion_secondary_at_v * spectrum_ion_flux * width_ev *
                       (1.0 - potential_v / energy_ev);
            }

            for (std::size_t j = 0; j < tables.ion_component_flux.size(); ++j)
            {
                const auto& population = tables.ion_populations[j];
                const double barrier =
                    std::exp(-potential_v / std::max(1.0e-6, population.temperature_ev));
                xsi += ion_secondary_at_v * tables.ion_component_flux[j][i] * width_ev * barrier;
            }

            xse += see_at_v * spectrum_electron_flux * width_ev;
            xb += backscatter_at_v * spectrum_electron_flux * width_ev;
            for (std::size_t j = 0; j < tables.electron_component_flux.size(); ++j)
            {
                xse += see_at_zero * tables.electron_component_flux[j][i] * width_ev;
                xb += backscatter_at_zero * tables.electron_component_flux[j][i] * width_ev;
            }
        }
        else
        {
            if (tables.average_spectrum_electron_temperature_ev > 0.0)
            {
                xe += spectrum_electron_flux * width_ev *
                      std::exp(potential_v /
                               std::max(1.0e-6, kLegacyLeoBodySpectrumDecayScale *
                                                     tables.average_spectrum_electron_temperature_ev));
            }
            for (std::size_t j = 0; j < tables.electron_component_flux.size(); ++j)
            {
                const auto& population = tables.electron_populations[j];
                x1 += tables.electron_component_flux[j][i] * width_ev *
                      std::exp(potential_v / std::max(1.0e-6, population.temperature_ev));
            }

            const double shifted_energy_ev = energy_ev - potential_v;
            const bool shifted_within_spectrum =
                !tables.spectrum_energy_grid_ev.empty() &&
                shifted_energy_ev >= tables.spectrum_energy_grid_ev.front() &&
                shifted_energy_ev <= tables.spectrum_energy_grid_ev.back();
            const double shifted_spectrum_flux =
                shifted_within_spectrum
                    ? std::max(0.0, legacyInterpolateSpectrum(tables.spectrum_energy_grid_ev,
                                                              tables.smoothed_electron_flux,
                                                              shifted_energy_ev))
                    : 0.0;

            xi += spectrum_ion_flux * width_ev;
            xsi += ion_secondary_at_v * total_ion_flux * width_ev;

            const double shifted_weight =
                tables.average_spectrum_electron_temperature_ev > 0.0
                    ? std::exp(
                          potential_v /
                          std::max(1.0e-6, kLegacyLeoBodySpectrumDecayScale *
                                                tables.average_spectrum_electron_temperature_ev))
                    : 0.0;
            xse += see_at_zero * shifted_spectrum_flux * width_ev * shifted_weight;
            xb += backscatter_at_zero * shifted_spectrum_flux * width_ev * shifted_weight;
            for (std::size_t j = 0; j < tables.electron_component_flux.size(); ++j)
            {
                const auto& population = tables.electron_populations[j];
                const double barrier =
                    std::exp(potential_v / std::max(1.0e-6, population.temperature_ev));
                xse += see_at_zero * tables.electron_component_flux[j][i] * width_ev * barrier;
                xb += backscatter_at_zero * tables.electron_component_flux[j][i] * width_ev * barrier;
            }
        }
    }

    x2_na_per_m2 =
        kLegacyLeoBodyThermalIonCollectionFactor *
        legacyThermalIonCollectionNaPerM2(tables.ion_populations, potential_v, potential_v > 0.0);

    currents.electron_current_a_per_m2 = -kPi * kElementaryCharge * (x1 + xe);
    currents.ion_current_a_per_m2 = kPi * kElementaryCharge * xi + x2_na_per_m2 * 1.0e-9;
    currents.secondary_emission_a_per_m2 =
        kPi * kElementaryCharge * xse *
        (potential_v > 0.0 ? std::exp(-potential_v / (2.0 * kTseEv)) : 1.0) *
        kLegacyLeoBodySecondaryEscapeFactor;
    currents.ion_secondary_emission_a_per_m2 =
        kPi * kElementaryCharge * xsi *
        (potential_v > 0.0 ? std::exp(-potential_v / (2.0 * kTsiEv)) : 1.0);
    currents.backscatter_emission_a_per_m2 = kPi * kElementaryCharge * xb;
    currents.photo_emission_a_per_m2 =
        potential_v > 0.0
            ? 0.25 * std::max(0.0, effective_config.body_photo_current_density_a_per_m2) *
                  std::exp(-potential_v / (2.0 * kTphEv))
            : 0.25 * std::max(0.0, effective_config.body_photo_current_density_a_per_m2);
    currents.ram_ion_current_a_per_m2 = 0.0;
    currents.conduction_current_a_per_m2 = 0.0;
    currents.current_derivative_a_per_m2_per_v = 0.0;
    currents.total_current_a_per_m2 =
        currents.electron_current_a_per_m2 + currents.ion_current_a_per_m2 +
        currents.secondary_emission_a_per_m2 + currents.ion_secondary_emission_a_per_m2 +
        currents.backscatter_emission_a_per_m2 + currents.photo_emission_a_per_m2;
    return currents;
}

double legacySpectrumCharacteristicEnergyEv(const Particle::ResolvedSpectrum& spectrum,
                                            double fallback_energy_ev)
{
    if (spectrum.energy_grid_ev.size() >= 2 &&
        spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size())
    {
        double weighted_energy = 0.0;
        double total_flux = 0.0;
        for (std::size_t i = 1; i < spectrum.energy_grid_ev.size(); ++i)
        {
            const double e0 = std::max(0.0, spectrum.energy_grid_ev[i - 1]);
            const double e1 = std::max(e0, spectrum.energy_grid_ev[i]);
            const double f0 = std::max(0.0, spectrum.differential_number_flux[i - 1]);
            const double f1 = std::max(0.0, spectrum.differential_number_flux[i]);
            const double width = e1 - e0;
            if (width <= 0.0)
            {
                continue;
            }
            const double average_flux = 0.5 * (f0 + f1);
            const double average_energy = 0.5 * (e0 + e1);
            weighted_energy += average_energy * average_flux * width;
            total_flux += average_flux * width;
        }
        if (total_flux > 0.0)
        {
            return weighted_energy / total_flux;
        }
    }
    return std::max(1.0e-6, fallback_energy_ev);
}

std::vector<LegacyPopulation> benchmarkElectronPopulations(const SurfaceChargingConfig& config)
{
    std::vector<LegacyPopulation> populations;
    if (config.has_electron_spectrum && !config.electron_spectrum.populations.empty())
    {
        for (const auto& population : config.electron_spectrum.populations)
        {
            populations.push_back(
                {population.density_m3, population.temperature_ev, population.mass_amu});
            if (populations.size() >= 2)
            {
                break;
            }
        }
    }
    if (populations.empty())
    {
        populations.push_back(
            {config.plasma.electron_density_m3, config.plasma.electron_temperature_ev,
             9.1093837015e-31 / 1.66053906660e-27});
    }
    return populations;
}

std::vector<LegacyPopulation> benchmarkIonPopulations(const SurfaceChargingConfig& config)
{
    std::vector<LegacyPopulation> populations;
    if (config.has_ion_spectrum && !config.ion_spectrum.populations.empty())
    {
        for (const auto& population : config.ion_spectrum.populations)
        {
            populations.push_back(
                {population.density_m3, population.temperature_ev, population.mass_amu});
            if (populations.size() >= 2)
            {
                break;
            }
        }
    }
    if (populations.empty())
    {
        populations.push_back(
            {config.plasma.ion_density_m3, config.plasma.ion_temperature_ev, config.plasma.ion_mass_amu});
    }
    return populations;
}

double integrateLegacyMaxwellian(double temperature_ev,
                                 const std::function<double(double)>& kernel)
{
    const double scale = std::max(1.0e-6, temperature_ev);
    const int samples = 96;
    const double upper = 24.0 * scale;
    const double step = upper / static_cast<double>(samples);
    double sum = 0.0;
    for (int i = 0; i <= samples; ++i)
    {
        const double energy = i * step;
        const double weight =
            (i == 0 || i == samples) ? 1.0 : ((i % 2 == 0) ? 2.0 : 4.0);
        sum += weight * kernel(energy) * energy * std::exp(-energy / scale);
    }
    const double integral = step * sum / 3.0;
    return integral / (scale * scale);
}

ReferenceCurrentBalanceConfig makeReferenceConfig(const SurfaceChargingConfig& config,
                                                  const SurfaceModelRuntimeState& state)
{
    ReferenceCurrentBalanceConfig reference_config;
    reference_config.plasma = config.plasma;
    reference_config.electron_spectrum = config.electron_spectrum;
    reference_config.ion_spectrum = config.ion_spectrum;
    reference_config.has_electron_spectrum = config.has_electron_spectrum;
    reference_config.has_ion_spectrum = config.has_ion_spectrum;
    reference_config.patch_material = config.material;
    reference_config.see_model = config.reference_see_model;
    reference_config.patch_incidence_angle_deg = config.patch_incidence_angle_deg;
    reference_config.patch_flow_angle_deg = config.patch_flow_angle_deg;
    reference_config.patch_thickness_m = config.dielectric_thickness_m;
    reference_config.patch_conductivity_s_per_m = state.effective_conductivity_s_per_m;
    reference_config.electron_collection_coefficient = config.electron_collection_coefficient;
    reference_config.ion_collection_coefficient = config.ion_collection_coefficient;
    reference_config.bulk_flow_velocity_m_per_s = config.bulk_flow_velocity_m_per_s;
    reference_config.flow_alignment_cosine = config.flow_alignment_cosine;
    reference_config.ion_directed_velocity_m_per_s = config.ion_directed_velocity_m_per_s;
    reference_config.electron_calibration_factor = state.electron_calibration_factor;
    reference_config.ion_calibration_factor = state.ion_calibration_factor;
    reference_config.photoelectron_temperature_ev =
        std::max(1.0e-3, config.photoelectron_temperature_ev);
    reference_config.enable_ram_current = config.regime == SurfaceChargingRegime::LeoFlowingPlasma;

    const double derived_photo_current =
        kElementaryCharge * config.emission.photon_flux_m2_s *
        std::clamp(config.material.getScalarProperty("photoelectron_yield", 0.0) *
                       config.emission.enhancement_factor,
                   0.0, 1.0);
    reference_config.body_photo_current_density_a_per_m2 =
        config.body_photo_current_density_a_per_m2 > 0.0 ? config.body_photo_current_density_a_per_m2
                                                         : derived_photo_current;
    reference_config.patch_photo_current_density_a_per_m2 =
        config.patch_photo_current_density_a_per_m2 > 0.0 ? config.patch_photo_current_density_a_per_m2
                                                          : derived_photo_current;

    reference_config.body_material = config.material;
    reference_config.body_material.setType(Mesh::MaterialType::CONDUCTOR);
    reference_config.body_material.setName(config.material.getName() + "_body");
    reference_config.body_material.setConductivity(
        std::max(1.0e-6, config.material.getScalarProperty("body_conductivity_s_per_m", 1.0e4)));
    return reference_config;
}

class LegacyBenchmarkCurrentModelBase : public SurfaceCurrentModel
{
  public:
    LegacyBenchmarkCurrentModelBase(SurfaceChargingConfig config, std::string name)
        : config_(std::move(config)), name_(std::move(name))
    {
    }

    SurfaceCurrents evaluate(ReferenceSurfaceRole role,
                             const SurfaceModelRuntimeState& state) const override
    {
        return evaluateLegacy(role, state);
    }

    double solveEquilibriumPotential(ReferenceSurfaceRole role,
                                     const SurfaceModelRuntimeState& state,
                                     double minimum_potential_v,
                                     double maximum_potential_v) const override
    {
        const double initial_step =
            config_.regime == SurfaceChargingRegime::GeoKineticPicLike ? 1.0e4 : 1.0e3;
        const double minimum_step =
            config_.regime == SurfaceChargingRegime::GeoKineticPicLike ? 1.0e-1 : 1.0e-2;

        double potential = std::clamp(0.0, minimum_potential_v, maximum_potential_v);
        const double net0 = evaluateLegacy(role, withPotential(state, role, potential)).total_current_a_per_m2;
        double step = net0 >= 0.0 ? initial_step : -initial_step;
        double previous_potential = potential;
        double previous_net = net0;

        while (std::abs(step) >= minimum_step)
        {
            double candidate = potential;
            double net_candidate = previous_net;
            bool bracketed = false;
            for (int iteration = 0; iteration < 20000; ++iteration)
            {
                candidate = std::clamp(candidate + step, minimum_potential_v, maximum_potential_v);
                net_candidate =
                    evaluateLegacy(role, withPotential(state, role, candidate)).total_current_a_per_m2;
                if (!std::isfinite(net_candidate) ||
                    candidate == minimum_potential_v || candidate == maximum_potential_v)
                {
                    break;
                }
                if (net_candidate * net0 <= 0.0)
                {
                    bracketed = true;
                    break;
                }
            }

            if (!bracketed)
            {
                break;
            }

            previous_potential = candidate - step;
            previous_net = evaluateLegacy(role, withPotential(state, role, previous_potential))
                               .total_current_a_per_m2;
            potential = previous_potential;
            step /= 10.0;
        }

        return std::clamp(potential, minimum_potential_v, maximum_potential_v);
    }

    double computeCurrentDerivative(ReferenceSurfaceRole role,
                                    const SurfaceModelRuntimeState& state) const override
    {
        const double role_potential =
            role == ReferenceSurfaceRole::Body ? state.body_potential_v : state.patch_potential_v;
        const double step = std::max(1.0e-3, 1.0e-2 * std::max(1.0, std::abs(role_potential)));
        const double jp =
            evaluateLegacy(role, withPotential(state, role, role_potential + step)).total_current_a_per_m2;
        const double jm =
            evaluateLegacy(role, withPotential(state, role, role_potential - step)).total_current_a_per_m2;
        return (jp - jm) / (2.0 * step);
    }

    bool recalibrate(const SurfaceModelRuntimeState&, std::size_t) override
    {
        latest_sample_ = PicMccCurrentSample{};
        return false;
    }

    double electronCalibrationFactor() const override
    {
        return 1.0;
    }

    double electronCalibrationFactor(const SurfaceModelRuntimeState&) const override
    {
        return 1.0;
    }

    double ionCalibrationFactor() const override
    {
        return 1.0;
    }

    double ionCalibrationFactor(const SurfaceModelRuntimeState&) const override
    {
        return 1.0;
    }

    const PicMccCurrentSample& latestLivePicSample() const override
    {
        return latest_sample_;
    }

    const PicMccCurrentSample& latestLivePicSample(
        const SurfaceModelRuntimeState&) const override
    {
        return latest_sample_;
    }

    std::string algorithmName() const override
    {
        return name_;
    }

    SurfaceCurrentAlgorithmMode algorithmMode() const override
    {
        return SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    }

  protected:
    virtual SurfaceCurrents evaluateLegacy(ReferenceSurfaceRole role,
                                           const SurfaceModelRuntimeState& state) const = 0;

    const SurfaceChargingConfig& config() const
    {
        return config_;
    }

    const Material::MaterialProperty& materialFor(ReferenceSurfaceRole role,
                                                  const SurfaceModelRuntimeState& state) const
    {
        if (role == ReferenceSurfaceRole::Body)
        {
            return config_.material;
        }
        const SurfaceChargingConfig effective = effectivePatchConfig(config_, state);
        patch_material_cache_ = effective.material;
        return patch_material_cache_;
    }

    static SurfaceModelRuntimeState withPotential(const SurfaceModelRuntimeState& state,
                                                  ReferenceSurfaceRole role,
                                                  double potential_v)
    {
        SurfaceModelRuntimeState updated = state;
        if (role == ReferenceSurfaceRole::Body)
        {
            updated.body_potential_v = potential_v;
        }
        else
        {
            updated.patch_potential_v = potential_v;
        }
        return updated;
    }

    double photoCurrent(ReferenceSurfaceRole role, const SurfaceChargingConfig& effective_config,
                        double potential_v) const
    {
        const double base = role == ReferenceSurfaceRole::Body
                                ? 0.25 * std::max(0.0, effective_config.body_photo_current_density_a_per_m2)
                                : std::max(0.0, effective_config.patch_photo_current_density_a_per_m2) *
                                      std::max(0.0, std::cos(effective_config.patch_incidence_angle_deg *
                                                             kPi / 180.0));
        if (potential_v > 0.0)
        {
            return base * std::exp(-potential_v /
                                   std::max(1.0e-3, 2.0 * effective_config.photoelectron_temperature_ev));
        }
        return base;
    }

    double conductionCurrent(ReferenceSurfaceRole role, const SurfaceModelRuntimeState& state,
                             const SurfaceChargingConfig& effective_config) const
    {
        if (role != ReferenceSurfaceRole::Patch)
        {
            return 0.0;
        }
        return std::max(0.0, state.effective_conductivity_s_per_m) *
               (state.body_potential_v - state.patch_potential_v) /
               std::max(1.0e-9, effective_config.dielectric_thickness_m);
    }

    double ramCurrent(ReferenceSurfaceRole role, const SurfaceChargingConfig& effective_config,
                      double potential_v) const
    {
        if (!effective_config.use_reference_current_balance ||
            effective_config.regime != SurfaceChargingRegime::LeoFlowingPlasma)
        {
            return 0.0;
        }
        const double flow_speed_m_per_s =
            std::max(0.0, effective_config.bulk_flow_velocity_m_per_s);
        if (flow_speed_m_per_s <= 0.0)
        {
            return 0.0;
        }
        return role == ReferenceSurfaceRole::Patch
                   ? legacyRamPatchCurrentDensity(
                         potential_v, effective_config.plasma.ion_density_m3,
                         effective_config.plasma.ion_temperature_ev,
                         effective_config.plasma.ion_mass_amu, flow_speed_m_per_s,
                         effective_config.patch_flow_angle_deg)
                   : legacyRamBodyCurrentDensity(
                         potential_v, effective_config.plasma.ion_density_m3,
                         effective_config.plasma.ion_temperature_ev,
                         effective_config.plasma.ion_mass_amu, flow_speed_m_per_s);
    }

  private:
    SurfaceChargingConfig config_;
    std::string name_;
    mutable Material::MaterialProperty patch_material_cache_{
        2, Mesh::MaterialType::DIELECTRIC, "patch"};
    PicMccCurrentSample latest_sample_{};
};

class LegacyCBenchmarkCurrentModel final : public LegacyBenchmarkCurrentModelBase
{
  public:
    explicit LegacyCBenchmarkCurrentModel(SurfaceChargingConfig config)
        : LegacyBenchmarkCurrentModelBase(std::move(config), "LegacyCBenchmarkModel")
    {
    }

  protected:
    SurfaceCurrents evaluateLegacy(ReferenceSurfaceRole role,
                                   const SurfaceModelRuntimeState& state) const override
    {
        const SurfaceChargingConfig effective_config =
            role == ReferenceSurfaceRole::Patch ? effectivePatchConfig(config(), state) : config();
        const auto& material = materialFor(role, state);
        const double potential = role == ReferenceSurfaceRole::Body ? state.body_potential_v
                                                                    : state.patch_potential_v;

        if (effective_config.regime == SurfaceChargingRegime::LeoFlowingPlasma)
        {
            if (role == ReferenceSurfaceRole::Body)
            {
                return evaluateLegacyLeoBodyCurrents(effective_config, material, potential);
            }

            SurfaceCurrents currents;
            const auto electron_populations = benchmarkElectronPopulations(effective_config);
            const auto ion_populations = benchmarkIonPopulations(effective_config);
            const double electron_characteristic_ev = std::max(
                effective_config.plasma.electron_temperature_ev,
                legacySpectrumCharacteristicEnergyEv(effective_config.electron_spectrum,
                                                     effective_config.plasma.electron_temperature_ev));
            const double ion_characteristic_ev = std::max(
                effective_config.plasma.ion_temperature_ev,
                legacySpectrumCharacteristicEnergyEv(effective_config.ion_spectrum,
                                                     effective_config.plasma.ion_temperature_ev));

            for (const auto& population : electron_populations)
            {
                const double mass_kg =
                    std::max(1.0e-6, population.mass_amu) * 1.66053906660e-27;
                const double base = legacyThermalCurrentDensity(population.density_m3,
                                                                population.temperature_ev, mass_kg);
                const double factor =
                    potential > 0.0
                        ? 1.0
                        : std::exp(potential / std::max(1.0e-6, electron_characteristic_ev));
                const double collection =
                    -base * std::max(0.0, effective_config.electron_collection_coefficient) * factor;
                currents.electron_current_a_per_m2 += collection;

                const double see_integral = integrateLegacyMaxwellian(
                    population.temperature_ev, [&](double energy_ev) {
                        const double shifted = potential > 0.0 ? energy_ev + potential : energy_ev;
                        return legacySeeYield(effective_config, material, shifted);
                    });
                const double backscatter_integral = integrateLegacyMaxwellian(
                    population.temperature_ev, [&](double energy_ev) {
                        const double shifted = potential > 0.0 ? energy_ev + potential : energy_ev;
                        return legacyBackscatterYield(material, shifted);
                    });
                const double emission_escape =
                    potential > 0.0
                        ? std::exp(-potential / std::max(
                                                   1.0e-3,
                                                   material.getScalarProperty(
                                                       "secondary_emission_escape_energy_ev", 2.0)))
                        : 1.0;
                currents.secondary_emission_a_per_m2 +=
                    std::abs(collection) * see_integral * emission_escape;
                currents.backscatter_emission_a_per_m2 +=
                    std::abs(collection) * backscatter_integral * emission_escape;
            }

            for (const auto& population : ion_populations)
            {
                const double mass_kg =
                    std::max(1.0e-6, population.mass_amu) * 1.66053906660e-27;
                const double base = legacyThermalCurrentDensity(population.density_m3,
                                                                population.temperature_ev, mass_kg);
                const double factor =
                    potential > 0.0
                        ? std::exp(-potential / std::max(1.0e-6, ion_characteristic_ev))
                        : 1.0;
                const double collection =
                    base * std::max(0.0, effective_config.ion_collection_coefficient) * factor;
                currents.ion_current_a_per_m2 += collection;

                const double ion_secondary_integral = integrateLegacyMaxwellian(
                    population.temperature_ev, [&](double energy_ev) {
                        const double shifted =
                            potential <= 0.0 ? energy_ev - potential : energy_ev;
                        return legacyIonSecondaryYield(material, shifted);
                    });
                const double emission_escape =
                    potential > 0.0
                        ? std::exp(-potential / std::max(
                                                   1.0e-3,
                                                   material.getScalarProperty(
                                                       "secondary_emission_escape_energy_ev", 2.0)))
                        : 1.0;
                currents.ion_secondary_emission_a_per_m2 +=
                    std::abs(collection) * ion_secondary_integral * emission_escape;
            }

            currents.photo_emission_a_per_m2 = photoCurrent(role, effective_config, potential);
            currents.conduction_current_a_per_m2 = conductionCurrent(role, state, effective_config);
            currents.ram_ion_current_a_per_m2 = ramCurrent(role, effective_config, potential);
            currents.current_derivative_a_per_m2_per_v = 0.0;
            currents.total_current_a_per_m2 =
                currents.electron_current_a_per_m2 + currents.ion_current_a_per_m2 +
                currents.secondary_emission_a_per_m2 + currents.ion_secondary_emission_a_per_m2 +
                currents.backscatter_emission_a_per_m2 + currents.photo_emission_a_per_m2 +
                currents.conduction_current_a_per_m2 + currents.ram_ion_current_a_per_m2;
            return currents;
        }

        const double electron_mass_kg = 9.1093837015e-31;
        const double ion_mass_kg =
            std::max(1.0, effective_config.plasma.ion_mass_amu) * 1.66053906660e-27;
        const double electron_base = legacyThermalCurrentDensity(
            effective_config.plasma.electron_density_m3, effective_config.plasma.electron_temperature_ev,
            electron_mass_kg);
        const double ion_base = legacyThermalCurrentDensity(
            effective_config.plasma.ion_density_m3, effective_config.plasma.ion_temperature_ev,
            ion_mass_kg);

        const double electron_factor =
            potential > 0.0
                ? (1.0 + potential / std::max(1.0e-6, effective_config.plasma.electron_temperature_ev))
                : std::exp(potential /
                           std::max(1.0e-6, effective_config.plasma.electron_temperature_ev));
        const double ion_factor =
            potential > 0.0
                ? std::exp(-potential / std::max(1.0e-6, effective_config.plasma.ion_temperature_ev))
                : (1.0 - potential / std::max(1.0e-6, effective_config.plasma.ion_temperature_ev));

        SurfaceCurrents currents;
        currents.electron_current_a_per_m2 =
            -electron_base * std::max(0.0, effective_config.electron_collection_coefficient) *
            electron_factor;
        currents.ion_current_a_per_m2 =
            ion_base * std::max(0.0, effective_config.ion_collection_coefficient) * ion_factor;

        const double electron_incident_energy =
            std::max(1.0e-6, effective_config.plasma.electron_temperature_ev + std::max(0.0, potential));
        const double ion_incident_energy =
            std::max(1.0e-6, effective_config.plasma.ion_temperature_ev + std::max(0.0, -potential));
        const double secondary_escape =
            potential > 0.0
                ? std::exp(-potential /
                           std::max(1.0e-3, material.getScalarProperty("secondary_emission_escape_energy_ev", 2.0)))
                : 1.0;

        currents.secondary_emission_a_per_m2 =
            std::abs(currents.electron_current_a_per_m2) *
            legacySeeYield(effective_config, material, electron_incident_energy) * secondary_escape;
        currents.ion_secondary_emission_a_per_m2 =
            std::abs(currents.ion_current_a_per_m2) *
            legacyIonSecondaryYield(material, ion_incident_energy) * secondary_escape;
        currents.backscatter_emission_a_per_m2 =
            std::abs(currents.electron_current_a_per_m2) *
            legacyBackscatterYield(material, electron_incident_energy) * secondary_escape;
        currents.photo_emission_a_per_m2 = photoCurrent(role, effective_config, potential);
        currents.conduction_current_a_per_m2 = conductionCurrent(role, state, effective_config);
        currents.ram_ion_current_a_per_m2 = ramCurrent(role, effective_config, potential);
        currents.current_derivative_a_per_m2_per_v = 0.0;
        currents.total_current_a_per_m2 =
            currents.electron_current_a_per_m2 + currents.ion_current_a_per_m2 +
            currents.secondary_emission_a_per_m2 + currents.ion_secondary_emission_a_per_m2 +
            currents.backscatter_emission_a_per_m2 + currents.photo_emission_a_per_m2 +
            currents.conduction_current_a_per_m2 + currents.ram_ion_current_a_per_m2;
        return currents;
    }
};

class LegacyMatlabBenchmarkCurrentModel final : public LegacyBenchmarkCurrentModelBase
{
  public:
    explicit LegacyMatlabBenchmarkCurrentModel(SurfaceChargingConfig config)
        : LegacyBenchmarkCurrentModelBase(std::move(config), "LegacyMatlabBenchmarkModel")
    {
    }

  protected:
    SurfaceCurrents evaluateLegacy(ReferenceSurfaceRole role,
                                   const SurfaceModelRuntimeState& state) const override
    {
        const SurfaceChargingConfig effective_config =
            role == ReferenceSurfaceRole::Patch ? effectivePatchConfig(config(), state) : config();
        const auto& material = materialFor(role, state);
        const double potential = role == ReferenceSurfaceRole::Body ? state.body_potential_v
                                                                    : state.patch_potential_v;

        SurfaceCurrents currents;
        const auto electron_populations = benchmarkElectronPopulations(effective_config);
        const auto ion_populations = benchmarkIonPopulations(effective_config);

        for (const auto& population : electron_populations)
        {
            const double mass_kg = std::max(1.0e-6, population.mass_amu) * 1.66053906660e-27;
            const double base =
                legacyThermalCurrentDensity(population.density_m3, population.temperature_ev, mass_kg);
            const double factor =
                potential > 0.0
                    ? (1.0 + potential / std::max(1.0e-6, population.temperature_ev))
                    : std::exp(potential / std::max(1.0e-6, population.temperature_ev));
            const double collection =
                -base * std::max(0.0, effective_config.electron_collection_coefficient) * factor;
            currents.electron_current_a_per_m2 += collection;

            const double see_integral = integrateLegacyMaxwellian(
                population.temperature_ev, [&](double energy_ev) {
                    const double shifted =
                        potential > 0.0 ? energy_ev + potential : energy_ev;
                    return legacySeeYield(effective_config, material, shifted);
                });
            const double backscatter_integral = integrateLegacyMaxwellian(
                population.temperature_ev, [&](double energy_ev) {
                    const double shifted =
                        potential > 0.0 ? energy_ev + potential : energy_ev;
                    return legacyBackscatterYield(material, shifted);
                });
            const double emission_escape =
                potential > 0.0
                    ? std::exp(-potential /
                               std::max(1.0e-3,
                                        material.getScalarProperty("secondary_emission_escape_energy_ev", 2.0)))
                    : 1.0;
            currents.secondary_emission_a_per_m2 +=
                std::abs(collection) * see_integral * emission_escape;
            currents.backscatter_emission_a_per_m2 +=
                std::abs(collection) * backscatter_integral * emission_escape;
        }

        for (const auto& population : ion_populations)
        {
            const double mass_kg = std::max(1.0e-6, population.mass_amu) * 1.66053906660e-27;
            const double base =
                legacyThermalCurrentDensity(population.density_m3, population.temperature_ev, mass_kg);
            const double factor =
                potential > 0.0
                    ? std::exp(-potential / std::max(1.0e-6, population.temperature_ev))
                    : (1.0 - potential / std::max(1.0e-6, population.temperature_ev));
            const double collection =
                base * std::max(0.0, effective_config.ion_collection_coefficient) * factor;
            currents.ion_current_a_per_m2 += collection;

            const double ion_secondary_integral = integrateLegacyMaxwellian(
                population.temperature_ev, [&](double energy_ev) {
                    const double shifted =
                        potential <= 0.0 ? energy_ev - potential : energy_ev;
                    return legacyIonSecondaryYield(material, shifted);
                });
            const double emission_escape =
                potential > 0.0
                    ? std::exp(-potential /
                               std::max(1.0e-3,
                                        material.getScalarProperty("secondary_emission_escape_energy_ev", 2.0)))
                    : 1.0;
            currents.ion_secondary_emission_a_per_m2 +=
                std::abs(collection) * ion_secondary_integral * emission_escape;
        }

        currents.photo_emission_a_per_m2 = photoCurrent(role, effective_config, potential);
        currents.conduction_current_a_per_m2 = conductionCurrent(role, state, effective_config);
        currents.ram_ion_current_a_per_m2 = ramCurrent(role, effective_config, potential);
        currents.current_derivative_a_per_m2_per_v = 0.0;
        currents.total_current_a_per_m2 =
            currents.electron_current_a_per_m2 + currents.ion_current_a_per_m2 +
            currents.secondary_emission_a_per_m2 + currents.ion_secondary_emission_a_per_m2 +
            currents.backscatter_emission_a_per_m2 + currents.photo_emission_a_per_m2 +
            currents.conduction_current_a_per_m2 + currents.ram_ion_current_a_per_m2;
        return currents;
    }
};

class ReferenceBackedSurfaceCurrentModel final : public SurfaceCurrentModel
{
  public:
    ReferenceBackedSurfaceCurrentModel(SurfaceChargingConfig config, SurfaceCurrentAlgorithmMode mode,
                                       std::string name)
        : config_(std::move(config)), mode_(mode), name_(std::move(name))
    {
    }

    SurfaceCurrents evaluate(ReferenceSurfaceRole role,
                             const SurfaceModelRuntimeState& state) const override
    {
        if (!isReferenceRegime(config_))
        {
            return SurfaceCurrents{};
        }

        const SurfaceChargingConfig effective_config =
            role == ReferenceSurfaceRole::Patch ? effectivePatchConfig(config_, state) : config_;
        ReferenceCurrentBalanceModel reference_model;
        reference_model.configure(makeReferenceConfig(effective_config, decorateState(state)));
        const auto terms =
            reference_model.evaluate(role, state.body_potential_v, state.patch_potential_v);

        SurfaceCurrents currents;
        currents.electron_current_a_per_m2 = terms.electron_collection_a_per_m2;
        currents.ion_current_a_per_m2 = terms.ion_collection_a_per_m2;
        currents.secondary_emission_a_per_m2 = terms.secondary_electron_a_per_m2;
        currents.ion_secondary_emission_a_per_m2 = terms.ion_secondary_electron_a_per_m2;
        currents.backscatter_emission_a_per_m2 = terms.backscatter_electron_a_per_m2;
        currents.photo_emission_a_per_m2 = terms.photoelectron_a_per_m2;
        currents.conduction_current_a_per_m2 = terms.conduction_a_per_m2;
        currents.ram_ion_current_a_per_m2 = terms.ram_ion_a_per_m2;
        if (role == ReferenceSurfaceRole::Patch)
        {
            currents.thermionic_emission_a_per_m2 =
                computePatchThermionicEmission(effective_config, state);
            currents.field_emission_a_per_m2 = computePatchFieldEmission(effective_config, state);
        }
        currents.current_derivative_a_per_m2_per_v = computeCurrentDerivative(role, state);
        currents.total_current_a_per_m2 = currents.electron_current_a_per_m2 +
                                          currents.ion_current_a_per_m2 +
                                          currents.secondary_emission_a_per_m2 +
                                          currents.ion_secondary_emission_a_per_m2 +
                                          currents.backscatter_emission_a_per_m2 +
                                          currents.photo_emission_a_per_m2 +
                                          currents.thermionic_emission_a_per_m2 +
                                          currents.field_emission_a_per_m2 +
                                          currents.conduction_current_a_per_m2 +
                                          currents.ram_ion_current_a_per_m2;
        return currents;
    }

    double solveEquilibriumPotential(ReferenceSurfaceRole role,
                                     const SurfaceModelRuntimeState& state,
                                     double minimum_potential_v,
                                     double maximum_potential_v) const override
    {
        const SurfaceChargingConfig effective_config =
            role == ReferenceSurfaceRole::Patch ? effectivePatchConfig(config_, state) : config_;
        ReferenceCurrentBalanceModel reference_model;
        reference_model.configure(makeReferenceConfig(effective_config, decorateState(state)));
        return reference_model.solveEquilibriumPotential(role, state.body_potential_v,
                                                         state.patch_potential_v,
                                                         minimum_potential_v,
                                                         maximum_potential_v);
    }

    double computeCurrentDerivative(ReferenceSurfaceRole role,
                                    const SurfaceModelRuntimeState& state) const override
    {
        const auto& calibration = calibrationStateFor(state);
        if (role == ReferenceSurfaceRole::Patch &&
            config_.regime == SurfaceChargingRegime::GeoKineticPicLike && calibration.latest_sample.valid &&
            std::isfinite(calibration.latest_sample.current_derivative_a_per_m2_per_v) &&
            std::abs(state.patch_potential_v - calibration.latest_sample.sampled_surface_potential_v) <=
                std::max(0.1, config_.live_pic_probe_delta_v))
        {
            return calibration.latest_sample.current_derivative_a_per_m2_per_v;
        }

        const SurfaceChargingConfig effective_config =
            role == ReferenceSurfaceRole::Patch ? effectivePatchConfig(config_, state) : config_;
        ReferenceCurrentBalanceModel reference_model;
        reference_model.configure(makeReferenceConfig(effective_config, decorateState(state)));
        return reference_model.computeCurrentDerivative(role, state.body_potential_v,
                                                        state.patch_potential_v);
    }

    bool recalibrate(const SurfaceModelRuntimeState& state,
                     std::size_t pic_calibration_samples) override
    {
        if (!config_.enable_pic_calibration && !config_.enable_live_pic_window)
        {
            calibrationStateForWrite(state).latest_sample = PicMccCurrentSample{};
            return false;
        }

        auto& calibration = calibrationStateForWrite(state);

        const double equilibrium = solveEquilibriumPotential(
            ReferenceSurfaceRole::Patch, state,
            -std::max(config_.max_abs_potential_v,
                      20.0 * std::max(1.0, config_.plasma.electron_temperature_ev)),
            std::max(config_.max_abs_potential_v,
                     20.0 * std::max(1.0, config_.plasma.electron_temperature_ev)));
        const double sample_potential =
            std::clamp(state.patch_potential_v + 0.25 * (equilibrium - state.patch_potential_v),
                       -config_.max_abs_potential_v, config_.max_abs_potential_v);

        bool updated = false;
        if (config_.enable_live_pic_window)
        {
            const SurfaceChargingConfig effective_config = effectivePatchConfig(config_, state);
            PicMccSurfaceSamplerConfig sampler_config;
            sampler_config.surface_area_m2 =
                state.surface_area_m2 > 0.0 ? state.surface_area_m2 : effective_config.surface_area_m2;
            sampler_config.gap_distance_m = std::clamp(
                state.effective_sheath_length_m,
                std::max(1.0e-5, effective_config.minimum_sheath_length_m),
                std::min(5.0e-2, std::max(effective_config.minimum_sheath_length_m,
                                          effective_config.maximum_sheath_length_m)));
            sampler_config.surface_potential_v = sample_potential;
            sampler_config.plasma_reference_potential_v = state.reference_plasma_potential_v;
            sampler_config.electron_flow_coupling = effective_config.electron_flow_coupling;
            sampler_config.bulk_flow_velocity_m_per_s = effective_config.bulk_flow_velocity_m_per_s;
            sampler_config.flow_alignment_cosine = effective_config.flow_alignment_cosine;
            sampler_config.flow_angle_deg = effective_config.patch_flow_angle_deg;
            sampler_config.ion_directed_velocity_m_per_s = effective_config.ion_directed_velocity_m_per_s;
            sampler_config.z_layers = std::max<std::size_t>(4, effective_config.live_pic_window_layers);
            sampler_config.particles_per_element =
                std::max<std::size_t>(1, effective_config.live_pic_particles_per_element);
            sampler_config.window_steps = std::max<std::size_t>(4, effective_config.live_pic_window_steps);
            sampler_config.enable_mcc = effective_config.enable_live_pic_mcc;
            sampler_config.plasma = effective_config.plasma;
            sampler_config.electron_spectrum = effective_config.electron_spectrum;
            sampler_config.ion_spectrum = effective_config.ion_spectrum;
            sampler_config.has_electron_spectrum = effective_config.has_electron_spectrum;
            sampler_config.has_ion_spectrum = effective_config.has_ion_spectrum;
            sampler_config.material = effective_config.material;
            calibration.latest_sample =
                sampler_.sampleWithDerivative(sampler_config, effective_config.live_pic_probe_delta_v);

            if (calibration.latest_sample.valid)
            {
                const SurfaceModelRuntimeState sample_state =
                    decorateState(withPatchPotential(state, sample_potential));
                ReferenceCurrentBalanceModel reference_model;
                reference_model.configure(makeReferenceConfig(effective_config, sample_state));
                const auto reference_terms =
                    reference_model.evaluate(ReferenceSurfaceRole::Patch, sample_state.body_potential_v,
                                             sample_state.patch_potential_v);
                updateFactorsFromReference(reference_terms, calibration);
                updated = true;
            }
        }

        if (!updated && config_.enable_pic_calibration)
        {
            calibration.latest_sample = PicMccCurrentSample{};
            const SurfaceChargingConfig effective_config = effectivePatchConfig(config_, state);
            const SurfaceModelRuntimeState sample_state =
                decorateState(withPatchPotential(state, sample_potential));
            ReferenceCurrentBalanceModel reference_model;
            reference_model.configure(makeReferenceConfig(effective_config, sample_state));
            const auto reference_terms =
                reference_model.evaluate(ReferenceSurfaceRole::Patch, sample_state.body_potential_v,
                                         sample_state.patch_potential_v);
            const double electron_scale =
                estimateElectronFluxPicLike(effective_config, sample_potential, pic_calibration_samples);
            const double ion_scale =
                estimateIonFluxPicLike(effective_config, sample_potential, pic_calibration_samples);
            if (std::abs(reference_terms.electron_collection_a_per_m2) > 1.0e-18)
            {
                calibration.electron_calibration_factor = std::clamp(
                    electron_scale / reference_terms.electron_collection_a_per_m2,
                    config_.electron_pic_calibration_min, config_.electron_pic_calibration_max);
                updated = true;
            }
            if (std::abs(reference_terms.ion_collection_a_per_m2) > 1.0e-18)
            {
                calibration.ion_calibration_factor = std::clamp(
                    ion_scale / reference_terms.ion_collection_a_per_m2,
                    config_.ion_pic_calibration_min, config_.ion_pic_calibration_max);
                updated = true;
            }
        }

        last_calibrated_key_ = stateKey(state);

        return updated;
    }

    double electronCalibrationFactor() const override
    {
        return calibrationStateForLastKey().electron_calibration_factor;
    }

    double electronCalibrationFactor(const SurfaceModelRuntimeState& state) const override
    {
        return calibrationStateFor(state).electron_calibration_factor;
    }

    double ionCalibrationFactor() const override
    {
        return calibrationStateForLastKey().ion_calibration_factor;
    }

    double ionCalibrationFactor(const SurfaceModelRuntimeState& state) const override
    {
        return calibrationStateFor(state).ion_calibration_factor;
    }

    const PicMccCurrentSample& latestLivePicSample() const override
    {
        return calibrationStateForLastKey().latest_sample;
    }

    const PicMccCurrentSample& latestLivePicSample(
        const SurfaceModelRuntimeState& state) const override
    {
        return calibrationStateFor(state).latest_sample;
    }

    std::string algorithmName() const override
    {
        return name_;
    }

    SurfaceCurrentAlgorithmMode algorithmMode() const override
    {
        return mode_;
    }

  private:
    struct PatchCalibrationState
    {
        double electron_calibration_factor = 1.0;
        double ion_calibration_factor = 1.0;
        PicMccCurrentSample latest_sample{};
    };

    std::string stateKey(const SurfaceModelRuntimeState& state) const
    {
        if (!state.node_name.empty())
        {
            return state.node_name;
        }
        return "node:" + std::to_string(state.node_index);
    }

    const PatchCalibrationState& calibrationStateFor(const SurfaceModelRuntimeState& state) const
    {
        const auto key = stateKey(state);
        const auto it = patch_calibration_states_.find(key);
        if (it != patch_calibration_states_.end())
        {
            return it->second;
        }
        return default_calibration_state_;
    }

    PatchCalibrationState& calibrationStateForWrite(const SurfaceModelRuntimeState& state)
    {
        const auto key = stateKey(state);
        return patch_calibration_states_[key];
    }

    const PatchCalibrationState& calibrationStateForLastKey() const
    {
        const auto it = patch_calibration_states_.find(last_calibrated_key_);
        if (it != patch_calibration_states_.end())
        {
            return it->second;
        }
        return default_calibration_state_;
    }

    SurfaceModelRuntimeState withPatchPotential(const SurfaceModelRuntimeState& state,
                                                double patch_potential_v) const
    {
        SurfaceModelRuntimeState updated = state;
        const double old_drop =
            std::abs(state.reference_plasma_potential_v - state.patch_potential_v);
        updated.patch_potential_v = patch_potential_v;
        const double new_drop =
            std::abs(state.reference_plasma_potential_v - updated.patch_potential_v);
        if (old_drop > 1.0e-9)
        {
            const double scale = new_drop / old_drop;
            updated.normal_electric_field_v_per_m *= scale;
            updated.local_charge_density_c_per_m3 *= scale;
        }
        return updated;
    }

    SurfaceModelRuntimeState decorateState(const SurfaceModelRuntimeState& state) const
    {
        SurfaceModelRuntimeState updated = state;
        const auto& calibration = calibrationStateFor(state);
        updated.electron_calibration_factor = calibration.electron_calibration_factor;
        updated.ion_calibration_factor = calibration.ion_calibration_factor;
        updated.live_pic_sample = &calibration.latest_sample;
        return updated;
    }

    void updateFactorsFromReference(const ReferenceCurrentComponents& reference_terms,
                                    PatchCalibrationState& calibration)
    {
        if (std::abs(reference_terms.electron_collection_a_per_m2) > 1.0e-18)
        {
            calibration.electron_calibration_factor = std::clamp(
                calibration.latest_sample.electron_collection_current_density_a_per_m2 /
                    reference_terms.electron_collection_a_per_m2,
                config_.electron_pic_calibration_min, config_.electron_pic_calibration_max);
        }
        if (std::abs(reference_terms.ion_collection_a_per_m2) > 1.0e-18)
        {
            calibration.ion_calibration_factor = std::clamp(
                calibration.latest_sample.ion_collection_current_density_a_per_m2 /
                    reference_terms.ion_collection_a_per_m2,
                config_.ion_pic_calibration_min, config_.ion_pic_calibration_max);
        }
    }

    double estimateElectronFluxPicLike(const SurfaceChargingConfig& effective_config,
                                       double surface_potential_v, std::size_t samples) const
    {
        const std::size_t sample_count = std::max<std::size_t>(64, samples);
        const double thermal_velocity =
            std::sqrt(std::max(1.0e-3, effective_config.plasma.electron_temperature_ev) * kElementaryCharge /
                      9.1093837015e-31);
        double modifier = surface_potential_v <= 0.0 ? safeExp(surface_potential_v /
                                                                   std::max(0.25, effective_config.plasma.electron_temperature_ev))
                                                     : (1.0 + std::min(5.0, surface_potential_v /
                                                                               std::max(0.25, effective_config.plasma.electron_temperature_ev)));
        modifier = std::max(0.0, modifier);
        return kElementaryCharge * effective_config.plasma.electron_density_m3 * thermal_velocity *
               modifier * std::max(0.0, effective_config.electron_collection_coefficient) *
               (1.0 - 1.0 / static_cast<double>(sample_count + 1));
    }

    double estimateIonFluxPicLike(const SurfaceChargingConfig& effective_config,
                                  double surface_potential_v, std::size_t samples) const
    {
        const std::size_t sample_count = std::max<std::size_t>(64, samples);
        const double ion_mass =
            std::max(1.0, effective_config.plasma.ion_mass_amu) * 1.66053906660e-27;
        const double thermal_velocity =
            std::sqrt(std::max(1.0e-3, effective_config.plasma.ion_temperature_ev) * kElementaryCharge /
                      ion_mass);
        double modifier = surface_potential_v > 0.0 ? safeExp(-surface_potential_v /
                                                                  std::max(0.25, effective_config.plasma.electron_temperature_ev))
                                                    : std::sqrt(1.0 + std::min(50.0, -surface_potential_v /
                                                                                         std::max(0.25, effective_config.plasma.electron_temperature_ev)));
        modifier = std::max(0.0, modifier);
        return kElementaryCharge * effective_config.plasma.ion_density_m3 * thermal_velocity *
               modifier * std::max(0.0, effective_config.ion_collection_coefficient) *
               (1.0 - 1.0 / static_cast<double>(sample_count + 1));
    }

    SurfaceChargingConfig config_;
    SurfaceCurrentAlgorithmMode mode_;
    std::string name_;
    PatchCalibrationState default_calibration_state_{};
    std::unordered_map<std::string, PatchCalibrationState> patch_calibration_states_;
    std::string last_calibrated_key_ = "patch";
    PicMccSurfaceCurrentSampler sampler_;
};

class LumpedSurfaceCapacitanceModel final : public SurfaceCapacitanceModel
{
  public:
    double computeCapacitancePerArea(const SurfaceChargingConfig& config) const override
    {
        const double thickness = std::max(1.0e-8, config.dielectric_thickness_m);
        const double eps_r = std::max(1.0, config.material.getPermittivity());
        if (config.derive_capacitance_from_material)
        {
            return kEpsilon0 * eps_r / thickness;
        }
        return std::max(1.0e-12, config.capacitance_per_area_f_per_m2);
    }

    double computePatchCapacitanceF(const SurfaceChargingConfig&,
                                    double capacitance_per_area_f_per_m2,
                                    double patch_area_m2) const override
    {
        return std::max(0.0, capacitance_per_area_f_per_m2 * patch_area_m2);
    }

    std::string modelName() const override
    {
        return "LumpedCapacitanceModel";
    }
};

class DenseSurfaceVoltageModel final : public SurfaceVoltageModel
{
  public:
    std::string modelName() const override
    {
        return "DenseSurfaceVoltageModel";
    }
};

class LumpedPotentialReferenceModel final : public PotentialReferenceModel
{
  public:
    double plasmaReferencePotentialV(const SurfaceChargingConfig& config,
                                     const SurfaceModelRuntimeState& state) const override
    {
        return std::isfinite(state.reference_plasma_potential_v) ? state.reference_plasma_potential_v
                                                                 : config.live_pic_reference_potential_v;
    }

    std::string modelName() const override
    {
        return "LumpedPotentialReferenceModel";
    }
};

class SheathElectricFieldProvider final : public ElectricFieldProvider
{
  public:
    double normalFieldVPerM(const SurfaceChargingConfig& config,
                            const SurfaceModelRuntimeState& state) const override
    {
        const SurfaceChargingConfig effective_config =
            isPatchNodeName(state.node_name) ? effectivePatchConfig(config, state) : config;
        const double sheath_length = std::clamp(
            state.effective_sheath_length_m,
            std::max(1.0e-6, effective_config.minimum_sheath_length_m),
            std::max(effective_config.minimum_sheath_length_m, effective_config.maximum_sheath_length_m));
        const double reference_drop_v =
            std::abs(state.reference_plasma_potential_v - state.patch_potential_v);
        const double sheath_field_v_per_m = reference_drop_v / sheath_length;
        const double space_charge_field_v_per_m =
            0.5 * std::abs(state.local_charge_density_c_per_m3) * sheath_length / kEpsilon0;
        const double enhancement =
            std::max(1.0, effective_config.emission.enhancement_factor);
        const double max_field_v_per_m = std::max(
            1.0e6,
            effective_config.material.getScalarProperty("max_provider_field_v_per_m", 1.0e10));
        return std::clamp(enhancement * (sheath_field_v_per_m + space_charge_field_v_per_m), 0.0,
                          max_field_v_per_m);
    }

    std::string modelName() const override
    {
        return "SheathElectricFieldProvider";
    }
};

class BoltzmannVolumeChargeProvider final : public VolumeChargeProvider
{
  public:
    double localChargeDensityCPerM3(const SurfaceChargingConfig& config,
                                    const SurfaceModelRuntimeState& state) const override
    {
        const SurfaceChargingConfig effective_config =
            isPatchNodeName(state.node_name) ? effectivePatchConfig(config, state) : config;
        const double electron_temperature_ev =
            std::max(1.0e-3, effective_config.plasma.electron_temperature_ev);
        const double normalized_potential =
            (state.patch_potential_v - state.reference_plasma_potential_v) / electron_temperature_ev;
        const double electron_density_m3 =
            std::max(0.0, effective_config.plasma.electron_density_m3) *
            safeExp(std::clamp(normalized_potential, -40.0, 8.0));
        double ion_density_m3 = std::max(0.0, effective_config.plasma.ion_density_m3);
        if (normalized_potential <= 0.0)
        {
            ion_density_m3 *= std::sqrt(1.0 + std::min(50.0, -normalized_potential));
        }
        else
        {
            ion_density_m3 *= safeExp(-std::min(20.0, normalized_potential));
        }
        return kElementaryCharge * (ion_density_m3 - electron_density_m3);
    }

    std::string modelName() const override
    {
        return "BoltzmannVolumeChargeProvider";
    }
};

class NullElectricFieldProvider final : public ElectricFieldProvider
{
  public:
    double normalFieldVPerM(const SurfaceChargingConfig&,
                            const SurfaceModelRuntimeState&) const override
    {
        return 0.0;
    }

    std::string modelName() const override
    {
        return "NullElectricFieldProvider";
    }
};

class NullVolumeChargeProvider final : public VolumeChargeProvider
{
  public:
    double localChargeDensityCPerM3(const SurfaceChargingConfig&,
                                    const SurfaceModelRuntimeState&) const override
    {
        return 0.0;
    }

    std::string modelName() const override
    {
        return "NullVolumeChargeProvider";
    }
};

class FieldAwareBubbleCapacitanceEstimator final : public BubbleCapacitanceEstimator
{
  public:
    double estimateCapacitancePerArea(const SurfaceChargingConfig& config,
                                      const SurfaceModelRuntimeState& state,
                                      double fallback_capacitance_per_area_f_per_m2) const override
    {
        const SurfaceChargingConfig effective_config =
            isPatchNodeName(state.node_name) ? effectivePatchConfig(config, state) : config;
        const double base_capacitance =
            std::max(1.0e-12, fallback_capacitance_per_area_f_per_m2);
        const double sheath_length = std::clamp(
            state.effective_sheath_length_m,
            std::max(1.0e-6, effective_config.minimum_sheath_length_m),
            std::max(effective_config.minimum_sheath_length_m, effective_config.maximum_sheath_length_m));
        const double bubble_capacitance =
            kEpsilon0 * std::max(1.0, effective_config.material.getPermittivity()) / sheath_length;
        const double series_capacitance =
            1.0 / (1.0 / base_capacitance + 1.0 / std::max(1.0e-12, bubble_capacitance));
        const double field_scale =
            std::abs(state.normal_electric_field_v_per_m) /
            std::max(1.0, std::abs(state.reference_plasma_potential_v - state.patch_potential_v) /
                               std::max(1.0e-6, sheath_length));
        const double charge_scale =
            std::abs(state.local_charge_density_c_per_m3) * sheath_length /
            std::max(1.0e-12, kEpsilon0 * std::max(1.0, std::abs(state.normal_electric_field_v_per_m)));
        const double base_weight =
            std::clamp(effective_config.material.getScalarProperty("bubble_capacitance_weight", 0.35),
                       0.0, 0.95);
        const double dynamic_weight =
            std::clamp(base_weight + 0.10 * std::log1p(field_scale) + 0.15 * charge_scale, 0.0, 0.95);
        return std::max(1.0e-12,
                        (1.0 - dynamic_weight) * base_capacitance +
                            dynamic_weight * series_capacitance);
    }

    std::string modelName() const override
    {
        return "FieldAwareBubbleCapacitanceEstimator";
    }
};

class ConstantBubbleCapacitanceEstimator final : public BubbleCapacitanceEstimator
{
  public:
    double estimateCapacitancePerArea(const SurfaceChargingConfig&,
                                      const SurfaceModelRuntimeState&,
                                      double fallback_capacitance_per_area_f_per_m2) const override
    {
        return fallback_capacitance_per_area_f_per_m2;
    }

    std::string modelName() const override
    {
        return "ConstantBubbleCapacitanceEstimator";
    }
};

class ConstantSurfaceReferenceStateProvider final : public SurfaceReferenceStateProvider
{
  public:
    double plasmaReferencePotentialV(const SurfaceChargingConfig& config,
                                     const SurfaceModelRuntimeState&) const override
    {
        return config.live_pic_reference_potential_v;
    }

    std::string modelName() const override
    {
        return "ConstantSurfaceReferenceStateProvider";
    }
};

class BranchWeightedReferenceGraphPropagator final : public SurfaceReferenceGraphPropagator
{
  public:
    double propagatedReferencePotentialV(const SurfaceChargingConfig& config,
                                         const SurfaceCircuitModel* circuit_model,
                                         const SurfaceModelRuntimeState& state,
                                         double base_reference_potential_v) const override
    {
        if (circuit_model == nullptr || state.node_index >= circuit_model->nodeCount() ||
            !circuit_model->nodeIsPatch(state.node_index))
        {
            return base_reference_potential_v;
        }

        double total_conductance_s = 0.0;
        double weighted_neighbor_potential_v = 0.0;
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            if (from_node != state.node_index && to_node != state.node_index)
            {
                continue;
            }
            const auto neighbor_node = from_node == state.node_index ? to_node : from_node;
            const double conductance_s = std::max(0.0, circuit_model->branchConductanceS(branch_index));
            if (neighbor_node >= circuit_model->nodeCount() || conductance_s <= 0.0)
            {
                continue;
            }
            total_conductance_s += conductance_s;
            weighted_neighbor_potential_v += conductance_s * circuit_model->nodePotential(neighbor_node);
        }
        if (total_conductance_s <= 0.0)
        {
            return base_reference_potential_v;
        }

        const double neighbor_reference_v = weighted_neighbor_potential_v / total_conductance_s;
        const double propagation_weight = std::clamp(
            config.material.getScalarProperty("graph_reference_propagation_weight", 0.20), 0.0, 0.95);
        return (1.0 - propagation_weight) * base_reference_potential_v +
               propagation_weight * neighbor_reference_v;
    }

    std::string modelName() const override
    {
        return "BranchWeightedReferenceGraphPropagator";
    }
};

class InterfaceAwareGraphCapacitanceMatrixProvider final
    : public SurfaceGraphCapacitanceMatrixProvider
{
  public:
    double diagonalCapacitanceF(const SurfaceChargingConfig& config,
                                const SurfaceCircuitModel* circuit_model,
                                std::size_t node_index) const override
    {
        if (circuit_model == nullptr || node_index >= circuit_model->nodeCount())
        {
            return 0.0;
        }

        const double area_m2 =
            std::max(0.0, circuit_model->nodeAreaM2(node_index));
        const double thickness_m = std::max(1.0e-8, config.dielectric_thickness_m);
        const double self_capacitance_f =
            kEpsilon0 * std::max(1.0, config.material.getPermittivity()) * area_m2 / thickness_m;

        double mutual_sum_f = 0.0;
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            if (from_node != node_index && to_node != node_index)
            {
                continue;
            }
            mutual_sum_f += mutualCapacitanceF(config, circuit_model, branch_index);
        }
        return std::max(0.0, self_capacitance_f + mutual_sum_f);
    }

    double mutualCapacitanceF(const SurfaceChargingConfig& config,
                              const SurfaceCircuitModel* circuit_model,
                              std::size_t branch_index) const override
    {
        if (circuit_model == nullptr || branch_index >= circuit_model->branchCount())
        {
            return 0.0;
        }

        const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
        const auto to_node = circuit_model->branchToNodeIndex(branch_index);
        const double from_area_m2 =
            from_node < circuit_model->nodeCount() ? std::max(0.0, circuit_model->nodeAreaM2(from_node))
                                                   : 0.0;
        const double to_area_m2 =
            to_node < circuit_model->nodeCount() ? std::max(0.0, circuit_model->nodeAreaM2(to_node))
                                                 : 0.0;
        const double coupled_area_m2 = std::min(from_area_m2, to_area_m2);
        const double thickness_m = std::max(1.0e-8, config.dielectric_thickness_m);
        const double branch_scale =
            std::clamp(config.material.getScalarProperty("graph_mutual_capacitance_scale", 0.15), 0.0,
                       1.0);
        return branch_scale * kEpsilon0 * std::max(1.0, config.material.getPermittivity()) *
               coupled_area_m2 / thickness_m;
    }

    double rowSumCapacitanceF(const SurfaceChargingConfig& config,
                              const SurfaceCircuitModel* circuit_model,
                              std::size_t node_index) const override
    {
        if (circuit_model == nullptr || node_index >= circuit_model->nodeCount())
        {
            return 0.0;
        }

        double row_sum_f = std::abs(diagonalCapacitanceF(config, circuit_model, node_index));
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            if (from_node != node_index && to_node != node_index)
            {
                continue;
            }
            row_sum_f += std::abs(mutualCapacitanceF(config, circuit_model, branch_index));
        }
        return row_sum_f;
    }

    double graphCouplingMetric(const SurfaceChargingConfig& config,
                               const SurfaceCircuitModel* circuit_model) const override
    {
        if (circuit_model == nullptr)
        {
            return 0.0;
        }

        double metric = 0.0;
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            metric += std::abs(mutualCapacitanceF(config, circuit_model, branch_index)) *
                      std::max(0.0, circuit_model->branchConductanceS(branch_index));
        }
        return metric;
    }

    std::string matrixFamilyName() const override
    {
        return "InterfaceAwareDenseHeuristic";
    }

    std::string solverAdapterHint() const override
    {
        return "GraphCoupledFieldSolverAdapter";
    }

    bool exposesMutualMatrix() const override
    {
        return true;
    }

    bool exposesDiagonalMatrix() const override
    {
        return true;
    }

    std::string modelName() const override
    {
        return "InterfaceAwareGraphCapacitanceMatrixProvider";
    }
};

class GraphCoupledFieldSolverAdapter final : public SurfaceFieldSolverAdapter
{
  public:
    void adaptFieldState(const SurfaceChargingConfig& config,
                         const SurfaceCircuitModel* circuit_model,
                         SurfaceModelRuntimeState& state) const override
    {
        state.field_solver_reference_potential_v = state.reference_plasma_potential_v;
        state.field_solver_coupling_gain = 0.0;
        if (circuit_model == nullptr || state.node_index >= circuit_model->nodeCount() ||
            !circuit_model->nodeIsPatch(state.node_index))
        {
            return;
        }

        double weighted_neighbor_potential_v = 0.0;
        double weighted_neighbor_reference_v = 0.0;
        double total_weight = 0.0;
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            if (from_node != state.node_index && to_node != state.node_index)
            {
                continue;
            }
            const auto neighbor_node = from_node == state.node_index ? to_node : from_node;
            if (neighbor_node >= circuit_model->nodeCount())
            {
                continue;
            }
            const double conductance_s = std::max(0.0, circuit_model->branchConductanceS(branch_index));
            const double effective_weight = conductance_s > 0.0 ? conductance_s : 1.0e-12;
            total_weight += effective_weight;
            weighted_neighbor_potential_v +=
                effective_weight * circuit_model->nodePotential(neighbor_node);
            weighted_neighbor_reference_v +=
                effective_weight * circuit_model->nodePotential(neighbor_node);
        }
        if (total_weight <= 0.0)
        {
            return;
        }

        const double average_neighbor_potential_v = weighted_neighbor_potential_v / total_weight;
        const double average_neighbor_reference_v = weighted_neighbor_reference_v / total_weight;
        const double graph_weight = std::clamp(
            config.material.getScalarProperty("field_solver_adapter_weight", 0.12), 0.0, 0.45);
        const double field_gain = std::clamp(
            config.material.getScalarProperty("field_solver_field_gain", 0.25), 0.0, 1.0);
        const double reference_shift_v =
            average_neighbor_reference_v - state.reference_plasma_potential_v;
        state.field_solver_coupling_gain =
            graph_weight * std::abs(average_neighbor_potential_v - state.patch_potential_v) /
            std::max(1.0e-6, std::abs(state.reference_plasma_potential_v - state.patch_potential_v) +
                                std::abs(average_neighbor_potential_v - state.patch_potential_v));
        state.field_solver_reference_potential_v =
            state.reference_plasma_potential_v + graph_weight * reference_shift_v;
        const double field_delta_v_per_m =
            (average_neighbor_potential_v - state.patch_potential_v) /
            std::max(1.0e-6, state.effective_sheath_length_m);
        state.normal_electric_field_v_per_m += field_gain * graph_weight * field_delta_v_per_m;
    }

    std::string modelName() const override
    {
        return "GraphCoupledFieldSolverAdapter";
    }
};

class ExternalFileFieldSolverAdapter final : public SurfaceFieldSolverAdapter
{
  public:
    explicit ExternalFileFieldSolverAdapter(SurfaceChargingConfig config)
        : config_(std::move(config))
    {
    }

    void adaptFieldState(const SurfaceChargingConfig& config,
                         const SurfaceCircuitModel* circuit_model,
                         SurfaceModelRuntimeState& state) const override
    {
        fallback_.adaptFieldState(config, circuit_model, state);
        if (!config_.enable_external_field_solver_bridge &&
            !config_.enable_external_volume_solver_bridge)
        {
            return;
        }

        const auto result_map = loadExternalFieldSolveResultMap(config_.external_field_solver_result_path);
        const auto apply_result = [&](const ExternalFieldSolveNodeResult& result) {
            if (std::isfinite(result.reference_potential_v))
            {
                state.field_solver_reference_potential_v = result.reference_potential_v;
            }
            if (std::isfinite(result.normal_field_v_per_m))
            {
                state.normal_electric_field_v_per_m = result.normal_field_v_per_m;
            }
            if (std::isfinite(result.local_charge_density_c_per_m3))
            {
                state.local_charge_density_c_per_m3 = result.local_charge_density_c_per_m3;
            }
            if (std::isfinite(result.capacitance_scale))
            {
                state.field_solver_capacitance_scale =
                    std::clamp(result.capacitance_scale, 0.25, 4.0);
            }
            state.field_solver_coupling_gain =
                std::max(state.field_solver_coupling_gain, 1.0);
        };
        if (!result_map.empty())
        {
            if (const auto it = result_map.find(state.node_name); it != result_map.end())
            {
                apply_result(it->second);
                return;
            }

            const auto logical_id = logicalNodeIdFromRuntimeName(state.node_name);
            if (const auto it = result_map.find(logical_id); it != result_map.end())
            {
                apply_result(it->second);
                return;
            }

            for (const auto& mapping : config_.boundary_mappings)
            {
                if (!mappingTargetsNode(mapping, state.node_name))
                {
                    continue;
                }
                if (const auto it = result_map.find(mapping.boundary_group_id); it != result_map.end())
                {
                    apply_result(it->second);
                    return;
                }
            }
        }

        // External volume bridge results are ingested by SurfaceVolumetricSolverAdapter with
        // stability-gated blending, so field adapter keeps ownership only for field-result files.
    }

    std::string modelName() const override
    {
        return "ExternalFileFieldSolverAdapter";
    }

  private:
    SurfaceChargingConfig config_;
    GraphCoupledFieldSolverAdapter fallback_;
};

class VolumeStubVolumetricSolverAdapter final : public SurfaceVolumetricSolverAdapter
{
  public:
    void adaptVolumeState(const SurfaceChargingConfig& config,
                          const SurfaceCircuitModel* circuit_model,
                          SurfaceModelRuntimeState& state) const override
    {
        const double incoming_projection_weight_sum = std::max(1.0, state.volume_projection_weight_sum);
        const double incoming_mesh_coupling_gain =
            std::clamp(state.volume_mesh_coupling_gain, 0.0, 1.0);
        const double incoming_capacitance_scale =
            std::clamp(state.field_solver_capacitance_scale, 0.25, 4.0);
        const double incoming_normal_field_v_per_m = state.normal_electric_field_v_per_m;
        const double incoming_charge_density_c_per_m3 = state.local_charge_density_c_per_m3;
        const double incoming_reference_potential_v = state.field_solver_reference_potential_v;
        state.external_volume_feedback_blend_factor = 0.0;
        state.external_volume_feedback_mismatch_metric = 0.0;
        state.external_volume_feedback_applied = 0.0;
        std::unordered_map<std::string, ExternalVolumeSolveCellResult> external_volume_results;
        if (config.enable_external_volume_solver_bridge)
        {
            external_volume_results =
                loadExternalVolumeSolveResultMap(config.external_volume_solver_result_path);
            if (external_volume_results.empty())
            {
                const auto& cached_results = loadExternalVolumeResultMapCached(config);
                external_volume_results.insert(cached_results.begin(), cached_results.end());
            }
        }
        const bool has_external_volume_result = !external_volume_results.empty();
        const bool has_external_field_result =
            config.enable_external_field_solver_bridge &&
            !config.external_field_solver_result_path.empty() &&
            std::filesystem::exists(config.external_field_solver_result_path);
        const bool use_internal_volume_poisson =
            config.runtime_route == SurfaceRuntimeRoute::SCDATUnified &&
            (hasStructuredTopologyInput(config) || config.enable_external_volume_solver_bridge ||
             config.enable_external_field_solver_bridge);
        const double current_signature = buildSignature(config, circuit_model, state);

        if (!use_internal_volume_poisson)
        {
            state.pseudo_volume_m3 =
                std::max(1.0e-18, std::max(0.0, state.surface_area_m2) *
                                      std::max(1.0e-6, state.effective_sheath_length_m));
            state.volume_projection_weight_sum =
                std::max(incoming_projection_weight_sum,
                         projectionWeightSumForRuntimeNode(config, state.node_name));
            state.volume_mesh_coupling_gain = incoming_mesh_coupling_gain;
            state.volume_potential_v = incoming_reference_potential_v;
            state.deposited_charge_c =
                state.local_charge_density_c_per_m3 * state.pseudo_volume_m3;
            state.poisson_residual_v_m = 0.0;
            state.volume_solver_mode_id = 0.0;
            state.volume_solver_iterations = 0.0;
            state.volume_solver_linear_iterations = 0.0;
            state.volume_solver_converged = 0.0;
            state.volume_solver_residual_norm = 0.0;
            state.volume_solver_max_delta_v = 0.0;
            state.volume_solver_matrix_nnz = 0.0;
            state.volume_solver_cell_count = 0.0;
            state.external_volume_feedback_blend_factor = 0.0;
            state.external_volume_feedback_mismatch_metric = 0.0;
            state.external_volume_feedback_applied = 0.0;
            return;
        }

        if (!cached_solution_valid_ ||
            !std::isfinite(cached_signature_) ||
            std::abs(cached_signature_ - current_signature) > 1.0e-12)
        {
            cached_solution_ = solvePseudoVolumePoisson(
                config, circuit_model, state, &persisted_volume_potentials_v_,
                &persisted_charge_densities_c_per_m3_);
            cached_signature_ = current_signature;
            cached_solution_valid_ = true;
            persisted_volume_potentials_v_.clear();
            persisted_charge_densities_c_per_m3_.clear();
            for (const auto& cell : cached_solution_.cells)
            {
                persisted_volume_potentials_v_[cell.cell_id] = cell.solved_potential_v;
                persisted_charge_densities_c_per_m3_[cell.cell_id] =
                    cell.reconstructed_charge_density_c_per_m3;
            }
        }

        if (cached_solution_.cells.empty() || state.node_index >= cached_solution_.cells.size())
        {
            state.pseudo_volume_m3 =
                std::max(1.0e-18, std::max(0.0, state.surface_area_m2) *
                                      std::max(1.0e-6, state.effective_sheath_length_m));
            state.volume_projection_weight_sum = incoming_projection_weight_sum;
            state.volume_mesh_coupling_gain = incoming_mesh_coupling_gain;
            state.volume_solver_mode_id = 0.0;
            state.volume_solver_iterations = 0.0;
            state.volume_solver_linear_iterations = 0.0;
            state.volume_solver_converged = 0.0;
            state.volume_solver_residual_norm = 0.0;
            state.volume_solver_max_delta_v = 0.0;
            state.volume_solver_matrix_nnz = 0.0;
            state.volume_solver_cell_count = 0.0;
            state.external_volume_feedback_blend_factor = 0.0;
            state.external_volume_feedback_mismatch_metric = 0.0;
            state.external_volume_feedback_applied = 0.0;
            return;
        }

        const auto& cell = cached_solution_.cells[state.node_index];
        state.pseudo_volume_m3 = cell.pseudo_volume_m3;
        state.volume_projection_weight_sum =
            std::max(incoming_projection_weight_sum, cell.projection_weight_sum);
        state.volume_mesh_coupling_gain =
            std::clamp(std::max(incoming_mesh_coupling_gain, cell.coupling_gain), 0.0, 1.0);
        state.volume_potential_v = cell.solved_potential_v;
        state.deposited_charge_c = cell.deposited_charge_c;
        state.poisson_residual_v_m = cell.poisson_residual_v_m;
        state.volume_solver_mode_id = cached_solution_.solver_mode == "iterative" ? 2.0 : 1.0;
        state.volume_solver_iterations =
            static_cast<double>(cached_solution_.self_consistent_iterations_completed);
        state.volume_solver_linear_iterations =
            static_cast<double>(cached_solution_.linear_iterations_completed);
        state.volume_solver_converged = cached_solution_.converged ? 1.0 : 0.0;
        state.volume_solver_residual_norm = cached_solution_.linear_residual_norm;
        state.volume_solver_max_delta_v = cached_solution_.max_delta_v;
        state.volume_solver_matrix_nnz =
            static_cast<double>(cached_solution_.matrix_nonzeros);
        state.volume_solver_cell_count =
            static_cast<double>(cached_solution_.cells.size());

        const double reference_span_v =
            std::max(1.0, std::abs(cell.surface_potential_v - cell.reference_potential_v));
        const double baseline_field_v_per_m =
            std::abs(cell.surface_potential_v - cell.reference_potential_v) /
            std::max(1.0e-6, 0.5 * cell.sheath_length_m);
        const double internal_capacitance_scale =
            std::clamp(0.75 + 0.55 * (std::abs(cell.normal_field_v_per_m) /
                                      std::max(1.0, baseline_field_v_per_m)) +
                           0.10 * std::max(0.0, cell.projection_weight_sum - 1.0) +
                           0.15 * cell.coupling_gain,
                       0.25, 4.0);

        state.normal_electric_field_v_per_m = cell.normal_field_v_per_m;
        state.local_charge_density_c_per_m3 = cell.reconstructed_charge_density_c_per_m3;
        state.field_solver_coupling_gain =
            std::max(state.field_solver_coupling_gain, 0.4 + 0.6 * cell.coupling_gain);
        state.field_solver_capacitance_scale =
            std::clamp(std::max(incoming_capacitance_scale, internal_capacitance_scale), 0.25, 4.0);
        state.field_solver_reference_potential_v = cell.solved_potential_v;
        state.reference_plasma_potential_v = cell.solved_potential_v;

        auto blend_value = [](double base, double target, double alpha) {
            const double blend = std::clamp(alpha, 0.0, 1.0);
            return (1.0 - blend) * base + blend * target;
        };

        if (has_external_field_result)
        {
            const double field_bridge_blend = 1.0;
            state.field_solver_reference_potential_v =
                blend_value(state.field_solver_reference_potential_v, incoming_reference_potential_v,
                            field_bridge_blend);
            state.reference_plasma_potential_v = state.field_solver_reference_potential_v;
            state.normal_electric_field_v_per_m =
                blend_value(state.normal_electric_field_v_per_m, incoming_normal_field_v_per_m,
                            field_bridge_blend);
            state.local_charge_density_c_per_m3 =
                blend_value(state.local_charge_density_c_per_m3, incoming_charge_density_c_per_m3,
                            field_bridge_blend);
            state.field_solver_capacitance_scale =
                std::clamp(blend_value(state.field_solver_capacitance_scale,
                                       incoming_capacitance_scale, field_bridge_blend),
                           0.25, 4.0);
        }

        if (has_external_volume_result)
        {
            const ExternalVolumeSolveCellResult* external_result = nullptr;
            if (const auto it = external_volume_results.find(cell.cell_id);
                it != external_volume_results.end())
            {
                external_result = &it->second;
            }
            else if (!cell.boundary_group_id.empty())
            {
                if (const auto it = external_volume_results.find(cell.boundary_group_id);
                    it != external_volume_results.end())
                {
                    external_result = &it->second;
                }
                else
                {
                    for (const auto& [_, candidate] : external_volume_results)
                    {
                        if (candidate.boundary_group_id == cell.boundary_group_id)
                        {
                            external_result = &candidate;
                            break;
                        }
                    }
                }
            }
            if (external_result != nullptr)
            {
                const auto& external = *external_result;
                const double potential_tolerance_v =
                    std::max(std::clamp(config.volume_self_consistent_tolerance_v, 1.0e-9, 1.0),
                             0.5 * std::clamp(config.field_volume_outer_tolerance, 1.0e-9, 10.0));
                const double field_tolerance_v_per_m =
                    potential_tolerance_v / std::max(1.0e-6, 0.5 * cell.sheath_length_m);
                const double charge_tolerance_c_per_m3 =
                    kEpsilon0 * field_tolerance_v_per_m /
                    std::max(1.0e-6, 0.5 * cell.sheath_length_m);
                const double scale_tolerance = 0.1;

                double mismatch_metric = 0.0;
                if (std::isfinite(external.potential_v))
                {
                    mismatch_metric = std::max(
                        mismatch_metric,
                        std::abs(external.potential_v - cell.solved_potential_v) /
                            std::max(1.0e-9, potential_tolerance_v));
                }
                if (std::isfinite(external.reference_potential_v))
                {
                    mismatch_metric = std::max(
                        mismatch_metric,
                        std::abs(external.reference_potential_v - cell.solved_potential_v) /
                            std::max(1.0e-9, potential_tolerance_v));
                }
                if (std::isfinite(external.normal_field_v_per_m))
                {
                    mismatch_metric = std::max(
                        mismatch_metric,
                        std::abs(external.normal_field_v_per_m - cell.normal_field_v_per_m) /
                            std::max(1.0e-9, field_tolerance_v_per_m));
                }
                if (std::isfinite(external.local_charge_density_c_per_m3))
                {
                    mismatch_metric = std::max(
                        mismatch_metric,
                        std::abs(external.local_charge_density_c_per_m3 -
                                 cell.reconstructed_charge_density_c_per_m3) /
                            std::max(1.0e-18, charge_tolerance_c_per_m3));
                }
                if (std::isfinite(external.capacitance_scale))
                {
                    mismatch_metric = std::max(
                        mismatch_metric,
                        std::abs(external.capacitance_scale - state.field_solver_capacitance_scale) /
                            scale_tolerance);
                }

                const double residual_gate =
                    1.0 / (1.0 + std::abs(cached_solution_.linear_residual_norm) /
                                     std::max(1.0e-9, potential_tolerance_v));
                const double delta_gate =
                    1.0 / (1.0 + std::abs(cached_solution_.max_delta_v) /
                                     std::max(1.0e-9, potential_tolerance_v));
                const double convergence_gate = cached_solution_.converged ? 1.0 : 0.4;
                const double internal_confidence =
                    convergence_gate * 0.5 * (residual_gate + delta_gate);
                const double stability_gate = 1.0 / (1.0 + mismatch_metric);
                const double external_trust = std::clamp(
                    stability_gate + 0.5 * (1.0 - internal_confidence), 0.05, 1.0);
                const double external_blend = std::clamp(
                    std::clamp(config.volume_linear_relaxation, 0.2, 1.0) * external_trust, 0.05,
                    1.0);

                const double target_volume_potential =
                    std::isfinite(external.potential_v) ? external.potential_v : state.volume_potential_v;
                const double target_reference_potential =
                    std::isfinite(external.reference_potential_v)
                        ? external.reference_potential_v
                        : (std::isfinite(external.potential_v) ? external.potential_v
                                                               : state.field_solver_reference_potential_v);
                const double target_normal_field =
                    std::isfinite(external.normal_field_v_per_m) ? external.normal_field_v_per_m
                                                                  : state.normal_electric_field_v_per_m;
                const double target_charge_density =
                    std::isfinite(external.local_charge_density_c_per_m3)
                        ? external.local_charge_density_c_per_m3
                        : state.local_charge_density_c_per_m3;
                const double target_capacitance_scale =
                    std::isfinite(external.capacitance_scale)
                        ? std::clamp(external.capacitance_scale, 0.25, 4.0)
                        : state.field_solver_capacitance_scale;

                state.volume_potential_v =
                    blend_value(state.volume_potential_v, target_volume_potential, external_blend);
                state.field_solver_reference_potential_v = blend_value(
                    state.field_solver_reference_potential_v, target_reference_potential,
                    external_blend);
                state.reference_plasma_potential_v = state.field_solver_reference_potential_v;
                state.normal_electric_field_v_per_m = blend_value(
                    state.normal_electric_field_v_per_m, target_normal_field, external_blend);
                state.local_charge_density_c_per_m3 = blend_value(
                    state.local_charge_density_c_per_m3, target_charge_density, external_blend);
                state.field_solver_capacitance_scale =
                    std::clamp(blend_value(state.field_solver_capacitance_scale,
                                           target_capacitance_scale, external_blend),
                               0.25, 4.0);

                if (std::isfinite(external.coupling_gain))
                {
                    const double target_coupling_gain =
                        std::clamp(external.coupling_gain, 0.0, 1.0);
                    state.volume_mesh_coupling_gain = std::clamp(
                        blend_value(state.volume_mesh_coupling_gain, target_coupling_gain,
                                    external_blend),
                        0.0, 1.0);
                }
                if (std::isfinite(external.projection_weight_scale))
                {
                    const double target_projection_weight_sum =
                        std::max(1.0, state.volume_projection_weight_sum *
                                          std::clamp(external.projection_weight_scale, 0.5, 4.0));
                    state.volume_projection_weight_sum = std::max(
                        1.0, blend_value(state.volume_projection_weight_sum,
                                         target_projection_weight_sum, external_blend));
                }
                if (std::isfinite(external.sheath_length_scale))
                {
                    const double target_sheath_length_m =
                        std::clamp(state.effective_sheath_length_m *
                                       std::clamp(external.sheath_length_scale, 0.5, 2.0),
                                   std::max(1.0e-6, config.minimum_sheath_length_m),
                                   std::max(config.minimum_sheath_length_m,
                                            config.maximum_sheath_length_m));
                    state.effective_sheath_length_m = std::clamp(
                        blend_value(state.effective_sheath_length_m, target_sheath_length_m,
                                    external_blend),
                        std::max(1.0e-6, config.minimum_sheath_length_m),
                        std::max(config.minimum_sheath_length_m,
                                 config.maximum_sheath_length_m));
                    state.pseudo_volume_m3 = std::max(
                        1.0e-18, std::max(0.0, state.surface_area_m2) *
                                     state.effective_sheath_length_m);
                }

                state.external_volume_feedback_blend_factor = external_blend;
                state.external_volume_feedback_mismatch_metric = mismatch_metric;
                state.external_volume_feedback_applied = 1.0;
                state.field_solver_coupling_gain =
                    std::max(state.field_solver_coupling_gain, 0.6 + 0.4 * external_blend);
            }
        }

        const double surface_volume_reference_drop_v =
            std::abs(cell.surface_potential_v - state.field_solver_reference_potential_v);
        const double self_consistent_gain =
            std::clamp(surface_volume_reference_drop_v / reference_span_v, 0.0, 1.0);
        state.volume_mesh_coupling_gain =
            std::clamp(std::max(state.volume_mesh_coupling_gain, self_consistent_gain), 0.0, 1.0);
    }

    std::string projectionFamilyName() const override
    {
        return "NodeToPseudoCellProjection";
    }

    std::string meshFamilyName() const override
    {
        return "PseudoBoundaryVolumeMesh";
    }

    std::string requestSchemaVersion() const override
    {
        return "scdat.external_volume_request.v1";
    }

    std::string resultSchemaVersion() const override
    {
        return "scdat.external_volume_result.v1";
    }

    std::string modelName() const override
    {
        return "VolumeStubVolumetricSolverAdapter";
    }

    bool supportsBoundaryFaceMapping() const override
    {
        return true;
    }

    bool supportsProjectionWeights() const override
    {
        return true;
    }

  private:
    double externalResultSignature(const std::filesystem::path& path) const
    {
        if (path.empty() || !std::filesystem::exists(path))
        {
            return 0.0;
        }
        double signature = static_cast<double>(path.string().size()) * 1.0e-6;
        std::error_code error;
        const auto file_size = std::filesystem::file_size(path, error);
        if (!error)
        {
            signature += static_cast<double>(file_size) * 1.0e-12;
        }
        error.clear();
        const auto last_write = std::filesystem::last_write_time(path, error);
        if (!error)
        {
            signature += static_cast<double>(last_write.time_since_epoch().count()) * 1.0e-18;
        }
        return signature;
    }

    const std::unordered_map<std::string, ExternalVolumeSolveCellResult>&
    loadExternalVolumeResultMapCached(const SurfaceChargingConfig& config) const
    {
        const double signature = externalResultSignature(config.external_volume_solver_result_path);
        if (!cached_external_volume_result_valid_ ||
            !std::isfinite(cached_external_volume_result_signature_) ||
            std::abs(signature - cached_external_volume_result_signature_) > 1.0e-12)
        {
            cached_external_volume_result_map_ =
                loadExternalVolumeSolveResultMap(config.external_volume_solver_result_path);
            cached_external_volume_result_signature_ = signature;
            cached_external_volume_result_valid_ = true;
        }
        return cached_external_volume_result_map_;
    }

    double buildSignature(const SurfaceChargingConfig& config,
                          const SurfaceCircuitModel* circuit_model,
                          const SurfaceModelRuntimeState& state) const
    {
        auto accumulate_path_signature = [](double& signature, const std::filesystem::path& path,
                                            double path_scale, double size_scale,
                                            double time_scale) {
            if (path.empty())
            {
                return;
            }
            signature += static_cast<double>(path.string().size()) * path_scale;
            if (std::filesystem::exists(path))
            {
                std::error_code error;
                const auto file_size = std::filesystem::file_size(path, error);
                if (!error)
                {
                    signature += static_cast<double>(file_size) * size_scale;
                }
                error.clear();
                const auto last_write = std::filesystem::last_write_time(path, error);
                if (!error)
                {
                    signature += static_cast<double>(last_write.time_since_epoch().count()) * time_scale;
                }
            }
        };

        double signature = 0.0;
        signature += static_cast<double>(circuit_model ? circuit_model->nodeCount() : 0) * 1.0e-3;
        signature += static_cast<double>(circuit_model ? circuit_model->branchCount() : 0) * 1.0e-4;
        signature += state.reference_plasma_potential_v * 1.0e-5;
        signature += state.field_solver_reference_potential_v * 1.0e-5;
        if (circuit_model != nullptr)
        {
            for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
            {
                signature += (static_cast<double>(node_index) + 1.0) *
                             circuit_model->nodePotential(node_index) * 1.0e-6;
            }
        }
        accumulate_path_signature(signature, config.external_volume_mesh_path, 1.0e-7, 1.0e-12,
                                  1.0e-18);
        accumulate_path_signature(signature, config.external_surface_volume_projection_path, 1.0e-7,
                                  1.0e-12, 1.0e-18);
        return signature;
    }

    mutable bool cached_solution_valid_ = false;
    mutable double cached_signature_ = std::numeric_limits<double>::quiet_NaN();
    mutable PseudoVolumePoissonSolution cached_solution_;
    mutable std::unordered_map<std::string, double> persisted_volume_potentials_v_;
    mutable std::unordered_map<std::string, double> persisted_charge_densities_c_per_m3_;
    mutable bool cached_external_volume_result_valid_ = false;
    mutable double cached_external_volume_result_signature_ =
        std::numeric_limits<double>::quiet_NaN();
    mutable std::unordered_map<std::string, ExternalVolumeSolveCellResult>
        cached_external_volume_result_map_;
};

class DenseSurfaceCircuitModel final : public SurfaceCircuitModel
{
  public:
    void configure(const SurfaceChargingConfig& config, double capacitance_per_area_f_per_m2,
                   double body_potential_v, double patch_potential_v,
                   double effective_conductivity_s_per_m) override
    {
        coupling_.clear();
        body_node_index_ = 0;
        primary_patch_index_ = 0;
        primary_branch_index_ = 0;
        body_area_m2_ = 0.0;
        node_areas_m2_.clear();
        node_is_patch_.clear();
        node_names_.clear();
        patch_node_indices_.clear();
        patch_to_body_branch_indices_.clear();
        branch_uses_dynamic_conductance_.clear();
        found_patch_ = false;
        found_branch_ = false;

        if (!config.surface_nodes.empty())
        {
            for (std::size_t i = 0; i < config.surface_nodes.size(); ++i)
            {
                const auto& node_config = config.surface_nodes[i];
                SurfaceCircuitNode node;
                node.name = node_config.name.empty() ? (node_config.is_patch ? "patch" : "body")
                                                     : node_config.name;
                node.potential_v =
                    node_config.is_patch ? node_config.initial_potential_v : body_potential_v;
                if (node_config.is_patch && !found_patch_)
                {
                    node.potential_v = patch_potential_v;
                }
                const std::size_t node_index = coupling_.addNode(node);
                const double node_area =
                    node_config.area_m2 > 0.0 ? node_config.area_m2 : config.surface_area_m2;
                node.capacitance_f = node_config.capacitance_f > 0.0
                                         ? node_config.capacitance_f
                                         : nodeCapacitance(config, node_index, node.name, node_area,
                                                           capacitance_per_area_f_per_m2);
                node.fixed_potential = node_config.fixed_potential;
                node.fixed_value_v = node_config.fixed_value_v;
                coupling_.setNodeCapacitance(node_index, node.capacitance_f);
                coupling_.setNodeFixedPotential(node_index, node.fixed_potential, node.fixed_value_v);
                if (!node_config.is_patch)
                {
                    body_node_index_ = node_index;
                }
                else if (!found_patch_)
                {
                    primary_patch_index_ = node_index;
                    found_patch_ = true;
                }
                registerNodeMetadata(
                    node_index, node.name, node_config.is_patch,
                    node_area);
            }
            if (body_area_m2_ <= 0.0)
            {
                double accumulated_patch_area = 0.0;
                for (const auto& area : node_areas_m2_)
                {
                    accumulated_patch_area += std::max(0.0, area);
                }
                body_area_m2_ = std::max(config.surface_area_m2, accumulated_patch_area);
            }
            for (std::size_t i = 0; i < config.surface_branches.size(); ++i)
            {
                const auto& branch_config = config.surface_branches[i];
                SurfaceCircuitBranch branch;
                branch.from_node = branch_config.from_node;
                branch.to_node = branch_config.to_node;
                branch.conductance_s = branch_config.conductance_s > 0.0
                                           ? branch_config.conductance_s
                                           : (branch_config.resistance_ohm > 0.0
                                                  ? 1.0 / std::max(1.0e-12, branch_config.resistance_ohm)
                                                  : defaultConductanceForBranch(
                                                        branch_config.from_node, branch_config.to_node, config,
                                                        effective_conductivity_s_per_m));
                branch.bias_v = branch_config.bias_v;
                coupling_.addBranch(branch);
                branch_uses_dynamic_conductance_.push_back(
                    branch_config.conductance_s <= 0.0 && branch_config.resistance_ohm <= 0.0);
                updatePatchToBodyBranchIndex(i, branch);
                if (!found_branch_ && connectsPrimaryPatchToBody(branch))
                {
                    primary_branch_index_ = i;
                    found_branch_ = true;
                }
            }
            return;
        }

        SurfaceCircuitNode body_node;
        body_node.name = "body";
        body_node.potential_v = body_potential_v;
        body_node.capacitance_f = std::max(0.0, config.body_capacitance_f);
        body_node.shunt_conductance_s = std::max(0.0, config.body_leakage_conductance_s);
        body_node.fixed_potential = !config.body_floating;
        body_node.fixed_value_v = config.body_initial_potential_v;
        body_node_index_ = coupling_.addNode(body_node);
        registerNodeMetadata(body_node_index_, body_node.name, false, config.surface_area_m2);

        SurfaceCircuitNode patch_node;
        patch_node.name = "patch";
        patch_node.potential_v = patch_potential_v;
        primary_patch_index_ = coupling_.addNode(patch_node);
        patch_node.capacitance_f =
            nodeCapacitance(config, primary_patch_index_, patch_node.name, config.surface_area_m2,
                            capacitance_per_area_f_per_m2);
        coupling_.setNodeCapacitance(primary_patch_index_, patch_node.capacitance_f);
        found_patch_ = true;
        registerNodeMetadata(primary_patch_index_, patch_node.name, true, config.surface_area_m2);

        SurfaceCircuitBranch branch;
        branch.from_node = primary_patch_index_;
        branch.to_node = body_node_index_;
        branch.conductance_s = config.patch_body_resistance_ohm > 0.0
                                   ? 1.0 / std::max(1.0e-12, config.patch_body_resistance_ohm)
                                   : std::max(0.0, effective_conductivity_s_per_m) *
                                         config.surface_area_m2 /
                                         std::max(1.0e-9, config.dielectric_thickness_m);
        coupling_.addBranch(branch);
        branch_uses_dynamic_conductance_.push_back(config.patch_body_resistance_ohm <= 0.0);
        primary_branch_index_ = 0;
        found_branch_ = true;
        updatePatchToBodyBranchIndex(0, branch);
    }

    void updateNodePotentials(double body_potential_v, double patch_potential_v) override
    {
        if (body_node_index_ < coupling_.getNodes().size())
        {
            coupling_.setNodePotential(body_node_index_, body_potential_v);
        }
        if (primary_patch_index_ < coupling_.getNodes().size())
        {
            coupling_.setNodePotential(primary_patch_index_, patch_potential_v);
        }
    }

    void setNodePotential(std::size_t node_index, double potential_v) override
    {
        coupling_.setNodePotential(node_index, potential_v);
    }

    void setNodeCapacitance(std::size_t node_index, double capacitance_f) override
    {
        coupling_.setNodeCapacitance(node_index, std::max(0.0, capacitance_f));
    }

    void setBranchConductance(std::size_t branch_index, double conductance_s) override
    {
        coupling_.setBranchConductance(branch_index, conductance_s);
    }

    SurfaceCircuitAdvanceResult advanceImplicit(
        double dt, const SurfaceCircuitLinearization& linearization,
        double max_delta_potential_v) override
    {
        return coupling_.advanceImplicit(dt, linearization, max_delta_potential_v);
    }

    std::size_t bodyNodeIndex() const override
    {
        return body_node_index_;
    }

    std::size_t primaryPatchNodeIndex() const override
    {
        return primary_patch_index_;
    }

    std::size_t primaryPatchToBodyBranchIndex() const override
    {
        return primary_branch_index_;
    }

    std::size_t patchCount() const override
    {
        return patch_node_indices_.size();
    }

    std::size_t patchNodeIndex(std::size_t patch_ordinal) const override
    {
        return patch_ordinal < patch_node_indices_.size() ? patch_node_indices_[patch_ordinal] : 0;
    }

    std::size_t patchToBodyBranchIndex(std::size_t patch_ordinal) const override
    {
        return patch_ordinal < patch_to_body_branch_indices_.size()
                   ? patch_to_body_branch_indices_[patch_ordinal]
                   : std::numeric_limits<std::size_t>::max();
    }

    bool branchUsesDynamicConductance(std::size_t branch_index) const override
    {
        return branch_index < branch_uses_dynamic_conductance_.size()
                   ? branch_uses_dynamic_conductance_[branch_index]
                   : false;
    }

    std::size_t nodeCount() const override
    {
        return coupling_.getNodes().size();
    }

    std::size_t branchCount() const override
    {
        return coupling_.getBranches().size();
    }

    double nodeAreaM2(std::size_t node_index) const override
    {
        return node_index < node_areas_m2_.size() ? node_areas_m2_[node_index] : 0.0;
    }

    double bodyAreaM2() const override
    {
        return body_area_m2_;
    }

    double nodePotential(std::size_t node_index) const override
    {
        return node_index < coupling_.getNodes().size() ? coupling_.getNodes()[node_index].potential_v : 0.0;
    }

    bool nodeIsPatch(std::size_t node_index) const override
    {
        return node_index < node_is_patch_.size() ? node_is_patch_[node_index] : false;
    }

    std::string nodeName(std::size_t node_index) const override
    {
        return node_index < node_names_.size() ? node_names_[node_index] : std::string{};
    }

    std::size_t branchFromNodeIndex(std::size_t branch_index) const override
    {
        return branch_index < coupling_.getBranches().size()
                   ? coupling_.getBranches()[branch_index].from_node
                   : std::numeric_limits<std::size_t>::max();
    }

    std::size_t branchToNodeIndex(std::size_t branch_index) const override
    {
        return branch_index < coupling_.getBranches().size()
                   ? coupling_.getBranches()[branch_index].to_node
                   : std::numeric_limits<std::size_t>::max();
    }

    double branchConductanceS(std::size_t branch_index) const override
    {
        return branch_index < coupling_.getBranches().size()
                   ? coupling_.getBranches()[branch_index].conductance_s
                   : 0.0;
    }

    double branchResistanceOhm(std::size_t branch_index) const override
    {
        if (branch_index >= coupling_.getBranches().size())
        {
            return 0.0;
        }
        const double conductance_s = coupling_.getBranches()[branch_index].conductance_s;
        return conductance_s > 0.0 ? 1.0 / conductance_s : 0.0;
    }

    double branchBiasV(std::size_t branch_index) const override
    {
        return branch_index < coupling_.getBranches().size()
                   ? coupling_.getBranches()[branch_index].bias_v
                   : 0.0;
    }

    std::string branchName(std::size_t branch_index) const override
    {
        if (branch_index >= coupling_.getBranches().size())
        {
            return {};
        }
        const auto& branch = coupling_.getBranches()[branch_index];
        std::ostringstream label;
        label << nodeName(branch.from_node) << "->" << nodeName(branch.to_node);
        return label.str();
    }

    std::string modelName() const override
    {
        return "DenseSurfaceCircuitModel";
    }

  private:
    double nodeCapacitance(const SurfaceChargingConfig& config, std::size_t node_index,
                           const std::string& node_name, double area_m2,
                           double fallback_capacitance_per_area_f_per_m2) const
    {
        if (area_m2 <= 0.0)
        {
            return 0.0;
        }

        SurfaceModelRuntimeState state;
        state.node_index = node_index;
        state.node_name = node_name;
        state.surface_area_m2 = area_m2;
        const auto effective_config = effectivePatchConfig(config, state);
        const double thickness = std::max(1.0e-8, effective_config.dielectric_thickness_m);
        const double capacitance_per_area =
            effective_config.derive_capacitance_from_material
                ? kEpsilon0 * std::max(1.0, effective_config.material.getPermittivity()) / thickness
                : fallback_capacitance_per_area_f_per_m2;
        return std::max(0.0, capacitance_per_area * area_m2);
    }

    double defaultConductanceForBranch(std::size_t from_node, std::size_t to_node,
                                       const SurfaceChargingConfig& config,
                                       double effective_conductivity_s_per_m) const
    {
        const bool connects_body_patch =
            (from_node == body_node_index_ && to_node < node_is_patch_.size() && node_is_patch_[to_node]) ||
            (to_node == body_node_index_ && from_node < node_is_patch_.size() && node_is_patch_[from_node]);
        if (!connects_body_patch)
        {
            return 0.0;
        }

        const auto patch_node =
            from_node == body_node_index_ ? to_node : from_node;
        return std::max(0.0, effective_conductivity_s_per_m) *
               std::max(0.0, nodeAreaM2(patch_node)) /
               std::max(1.0e-9, config.dielectric_thickness_m);
    }

    void registerNodeMetadata(std::size_t node_index, const std::string& name, bool is_patch,
                              double area_m2)
    {
        if (node_index >= node_areas_m2_.size())
        {
            node_areas_m2_.resize(node_index + 1, 0.0);
            node_is_patch_.resize(node_index + 1, false);
            node_names_.resize(node_index + 1);
        }
        node_areas_m2_[node_index] = std::max(0.0, area_m2);
        node_is_patch_[node_index] = is_patch;
        node_names_[node_index] = name;
        if (is_patch)
        {
            patch_node_indices_.push_back(node_index);
            patch_to_body_branch_indices_.push_back(std::numeric_limits<std::size_t>::max());
        }
        else
        {
            body_area_m2_ = std::max(0.0, area_m2);
        }
    }

    void updatePatchToBodyBranchIndex(std::size_t branch_index, const SurfaceCircuitBranch& branch)
    {
        for (std::size_t patch_ordinal = 0; patch_ordinal < patch_node_indices_.size(); ++patch_ordinal)
        {
            const auto patch_node_index = patch_node_indices_[patch_ordinal];
            const bool connects = (branch.from_node == patch_node_index && branch.to_node == body_node_index_) ||
                                  (branch.from_node == body_node_index_ && branch.to_node == patch_node_index);
            if (connects &&
                patch_to_body_branch_indices_[patch_ordinal] == std::numeric_limits<std::size_t>::max())
            {
                patch_to_body_branch_indices_[patch_ordinal] = branch_index;
                if (patch_node_index == primary_patch_index_)
                {
                    primary_branch_index_ = branch_index;
                    found_branch_ = true;
                }
            }
        }
    }

    bool connectsPrimaryPatchToBody(const SurfaceCircuitBranch& branch) const
    {
        return (branch.from_node == primary_patch_index_ && branch.to_node == body_node_index_) ||
               (branch.from_node == body_node_index_ && branch.to_node == primary_patch_index_);
    }

    SurfaceCircuitCoupling coupling_;
    std::size_t body_node_index_ = 0;
    std::size_t primary_patch_index_ = 0;
    std::size_t primary_branch_index_ = 0;
    double body_area_m2_ = 0.0;
    std::vector<double> node_areas_m2_;
    std::vector<bool> node_is_patch_;
    std::vector<std::string> node_names_;
    std::vector<std::size_t> patch_node_indices_;
    std::vector<std::size_t> patch_to_body_branch_indices_;
    std::vector<bool> branch_uses_dynamic_conductance_;
    bool found_patch_ = false;
    bool found_branch_ = false;
};

class DefaultSurfaceScenarioOrchestrator final : public SurfaceScenarioOrchestrator
{
  public:
    explicit DefaultSurfaceScenarioOrchestrator(SurfaceChargingConfig config)
        : config_(std::move(config))
    {
    }

    std::unique_ptr<SurfaceCurrentModel> createCurrentModel() const override
    {
        const auto mode = config_.current_algorithm_mode;
        if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
        {
            if (config_.benchmark_source == SurfaceBenchmarkSource::MatlabGeo ||
                config_.benchmark_source == SurfaceBenchmarkSource::MatlabLeo)
            {
                return std::make_unique<LegacyMatlabBenchmarkCurrentModel>(config_);
            }
            return std::make_unique<LegacyCBenchmarkCurrentModel>(config_);
        }
        std::string name = "UnifiedReferenceCurrentModel";
        if (mode == SurfaceCurrentAlgorithmMode::LegacyRefCompatible)
        {
            name = "LegacyReferenceCurrentModel";
        }
        return std::make_unique<ReferenceBackedSurfaceCurrentModel>(config_, mode, name);
    }

    std::unique_ptr<SurfaceVoltageModel> createVoltageModel() const override
    {
        return std::make_unique<DenseSurfaceVoltageModel>();
    }

    std::unique_ptr<SurfaceCapacitanceModel> createCapacitanceModel() const override
    {
        return std::make_unique<LumpedSurfaceCapacitanceModel>();
    }

    std::unique_ptr<PotentialReferenceModel> createPotentialReferenceModel() const override
    {
        return std::make_unique<LumpedPotentialReferenceModel>();
    }

    std::unique_ptr<ElectricFieldProvider> createElectricFieldProvider() const override
    {
        if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
        {
            return std::make_unique<NullElectricFieldProvider>();
        }
        return std::make_unique<SheathElectricFieldProvider>();
    }

    std::unique_ptr<VolumeChargeProvider> createVolumeChargeProvider() const override
    {
        if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
        {
            return std::make_unique<NullVolumeChargeProvider>();
        }
        return std::make_unique<BoltzmannVolumeChargeProvider>();
    }

    std::unique_ptr<BubbleCapacitanceEstimator> createBubbleCapacitanceEstimator() const override
    {
        if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
        {
            return std::make_unique<ConstantBubbleCapacitanceEstimator>();
        }
        return std::make_unique<FieldAwareBubbleCapacitanceEstimator>();
    }

    std::unique_ptr<SurfaceReferenceStateProvider> createSurfaceReferenceStateProvider() const override
    {
        return std::make_unique<ConstantSurfaceReferenceStateProvider>();
    }

    std::unique_ptr<SurfaceReferenceGraphPropagator>
    createSurfaceReferenceGraphPropagator() const override
    {
        return std::make_unique<BranchWeightedReferenceGraphPropagator>();
    }

    std::unique_ptr<SurfaceGraphCapacitanceMatrixProvider>
    createSurfaceGraphCapacitanceMatrixProvider() const override
    {
        return std::make_unique<InterfaceAwareGraphCapacitanceMatrixProvider>();
    }

    std::unique_ptr<SurfaceFieldSolverAdapter> createSurfaceFieldSolverAdapter() const override
    {
        if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
        {
            return nullptr;
        }
        if (config_.enable_external_field_solver_bridge ||
            config_.enable_external_volume_solver_bridge)
        {
            return std::make_unique<ExternalFileFieldSolverAdapter>(config_);
        }
        return std::make_unique<GraphCoupledFieldSolverAdapter>();
    }

    std::unique_ptr<SurfaceVolumetricSolverAdapter> createSurfaceVolumetricSolverAdapter() const override
    {
        if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
        {
            return nullptr;
        }
        return std::make_unique<VolumeStubVolumetricSolverAdapter>();
    }

    std::unique_ptr<SurfaceCircuitModel> createCircuitModel() const override
    {
        return std::make_unique<DenseSurfaceCircuitModel>();
    }

  private:
    SurfaceChargingConfig config_;
};
} // namespace

std::unique_ptr<SurfaceScenarioOrchestrator>
makeSurfaceScenarioOrchestrator(const SurfaceChargingConfig& config)
{
    return std::make_unique<DefaultSurfaceScenarioOrchestrator>(normalizeTopologyConfig(config));
}

SurfaceChargingConfig normalizeSurfaceChargingConfig(const SurfaceChargingConfig& config)
{
    return normalizeTopologyConfig(config);
}

bool validateSurfaceChargingConfig(const SurfaceChargingConfig& config,
                                   std::string* error_message)
{
    return validateStructuredTopologyConfig(config, error_message);
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
