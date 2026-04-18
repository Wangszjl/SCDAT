#include "SurfaceBridgeArtifactWriters.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <optional>
#include <sstream>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kEpsilon0 = 8.8541878128e-12;

double latestHistoryValue(const std::vector<std::vector<double>>& history, std::size_t index)
{
    if (index >= history.size() || history[index].empty())
    {
        return 0.0;
    }
    return history[index].back();
}

std::string jsonEscape(const std::string& value)
{
    std::ostringstream escaped;
    for (const char ch : value)
    {
        switch (ch)
        {
        case '\\':
            escaped << "\\\\";
            break;
        case '"':
            escaped << "\\\"";
            break;
        case '\n':
            escaped << "\\n";
            break;
        case '\r':
            escaped << "\\r";
            break;
        case '\t':
            escaped << "\\t";
            break;
        default:
            escaped << ch;
            break;
        }
    }
    return escaped.str();
}

std::string logicalNodeIdFromName(const std::string& node_name)
{
    const auto separator = node_name.find(':');
    return separator == std::string::npos ? node_name : node_name.substr(separator + 1);
}

std::string boundaryGroupIdForNode(const SurfaceChargingConfig& config,
                                   const std::string& node_name)
{
    const auto logical_id = logicalNodeIdFromName(node_name);
    for (const auto& mapping : config.boundary_mappings)
    {
        if (mapping.node_id == node_name || mapping.node_id == logical_id)
        {
            return mapping.boundary_group_id;
        }
    }
    return {};
}

double projectionWeightForNode(const SurfaceChargingConfig& config,
                               const std::string& node_name)
{
    const auto boundary_group_id = boundaryGroupIdForNode(config, node_name);
    if (boundary_group_id.empty())
    {
        return 1.0;
    }

    std::size_t linked_boundary_faces = 0;
    for (const auto& mapping : config.boundary_mappings)
    {
        if (mapping.boundary_group_id == boundary_group_id)
        {
            linked_boundary_faces += 1;
        }
    }
    for (const auto& group : config.body_boundary_groups)
    {
        if (group.id == boundary_group_id)
        {
            linked_boundary_faces += group.boundary_face_ids.size();
        }
    }
    for (const auto& group : config.patch_boundary_groups)
    {
        if (group.id == boundary_group_id)
        {
            linked_boundary_faces += group.boundary_face_ids.size();
        }
    }

    return std::max(1.0, 1.0 + 0.25 * static_cast<double>(linked_boundary_faces));
}

std::string runtimeRouteName(SurfaceRuntimeRoute route)
{
    switch (route)
    {
    case SurfaceRuntimeRoute::SurfacePic:
        return "surface_pic";
    case SurfaceRuntimeRoute::SurfacePicHybrid:
        return "surface_pic_hybrid";
    case SurfaceRuntimeRoute::LegacyBenchmark:
        return "legacy_benchmark";
    case SurfaceRuntimeRoute::SCDATUnified:
    default:
        return "scdat_unified";
    }
}

std::string surfacePicStrategyName(SurfacePicStrategy strategy)
{
    switch (strategy)
    {
    case SurfacePicStrategy::SurfacePicDirect:
        return "surface_pic_direct";
    case SurfacePicStrategy::SurfacePicHybridReference:
        return "surface_pic_hybrid_reference";
    case SurfacePicStrategy::SurfacePicCalibrated:
    default:
        return "surface_pic_calibrated";
    }
}

std::string surfacePicRuntimeName(SurfacePicRuntimeKind runtime)
{
    switch (runtime)
    {
    case SurfacePicRuntimeKind::GraphCoupledSharedSurface:
        return "graph_coupled_shared_surface";
    case SurfacePicRuntimeKind::LocalWindowSampler:
    default:
        return "local_window_sampler";
    }
}

std::string volumeLinearSolverPolicyName(VolumeLinearSolverPolicy policy)
{
    switch (policy)
    {
    case VolumeLinearSolverPolicy::DenseOnly:
        return "dense_only";
    case VolumeLinearSolverPolicy::IterativeOnly:
        return "iterative_only";
    case VolumeLinearSolverPolicy::Auto:
    default:
        return "iterative_or_dense_auto";
    }
}

std::string branchTypeLabel(const SurfaceCircuitModel* circuit_model, std::size_t branch_index)
{
    if (circuit_model == nullptr)
    {
        return "resistive";
    }
    const auto resistance = circuit_model->branchResistanceOhm(branch_index);
    const auto conductance = circuit_model->branchConductanceS(branch_index);
    if (resistance <= 1.0e-12 && conductance > 1.0e9)
    {
        return "ideal_connection";
    }
    if (resistance > 1.0e6 || conductance < 1.0e-6)
    {
        return "weak_coupling";
    }
    return "resistive";
}

const SurfacePatchPhysicsConfig* findPatchPhysicsOverrideForReport(
    const SurfaceChargingConfig& config, std::size_t node_index, const std::string& node_name)
{
    for (const auto& patch_config : config.patch_physics_overrides)
    {
        if (patch_config.match_by_index && patch_config.node_index == node_index)
        {
            return &patch_config;
        }
        if (!patch_config.node_name.empty() && patch_config.node_name == node_name)
        {
            return &patch_config;
        }
    }
    return nullptr;
}

const Material::MaterialProperty& resolvePatchMaterialForReport(const SurfaceChargingConfig& config,
                                                                std::size_t node_index,
                                                                const std::string& node_name)
{
    const auto* patch_config = findPatchPhysicsOverrideForReport(config, node_index, node_name);
    if (patch_config != nullptr && patch_config->override_material)
    {
        return patch_config->material;
    }
    return config.material;
}

std::string graphNodeTypeLabel(bool is_patch)
{
    return is_patch ? "patch" : "body";
}

std::string graphNodeMaterialName(const SurfaceChargingConfig& config,
                                  std::size_t node_index,
                                  const std::string& node_name,
                                  bool is_patch)
{
    if (!is_patch)
    {
        return "conductor-body";
    }
    return resolvePatchMaterialForReport(config, node_index, node_name).getName();
}

std::string graphNodeOwnerLabel(const SurfaceChargingConfig& config,
                                const std::string& node_name,
                                bool is_patch)
{
    if (is_patch && node_name.rfind("patch:", 0) == 0)
    {
        const auto patch_id = node_name.substr(std::string("patch:").size());
        for (const auto& patch : config.patches)
        {
            if (patch.id == patch_id)
            {
                return patch.body_id;
            }
        }
        return patch_id;
    }
    if (node_name.rfind("body:", 0) == 0)
    {
        return node_name.substr(std::string("body:").size());
    }
    return node_name;
}

std::string graphBranchTypeLabel(const SurfaceCircuitModel* circuit_model, std::size_t branch_index)
{
    if (circuit_model == nullptr || branch_index >= circuit_model->branchCount())
    {
        return "unknown";
    }
    const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
    const auto to_node = circuit_model->branchToNodeIndex(branch_index);
    const bool from_patch =
        from_node < circuit_model->nodeCount() && circuit_model->nodeIsPatch(from_node);
    const bool to_patch =
        to_node < circuit_model->nodeCount() && circuit_model->nodeIsPatch(to_node);
    if (from_patch && to_patch)
    {
        return "patch-patch";
    }
    if (!from_patch && !to_patch)
    {
        return "body-body";
    }
    return "body-patch";
}

std::optional<SurfacePatchConfig> findPatchConfigForNode(const SurfaceChargingConfig& config,
                                                         const std::string& node_name)
{
    const auto logical_id = logicalNodeIdFromName(node_name);
    for (const auto& patch : config.patches)
    {
        if (patch.id == logical_id)
        {
            return patch;
        }
    }
    return std::nullopt;
}

std::array<double, 3> pseudoNodeCenterForNode(const SurfaceChargingConfig& config,
                                              const SurfaceCircuitModel* circuit_model,
                                              std::size_t node_index)
{
    const double radius_m = 0.15 + 0.03 * static_cast<double>(node_index);
    if (circuit_model == nullptr || node_index >= circuit_model->nodeCount())
    {
        return {radius_m, 0.0, 0.0};
    }

    const auto node_name = circuit_model->nodeName(node_index);
    if (!circuit_model->nodeIsPatch(node_index))
    {
        const double angle = 0.45 * static_cast<double>(node_index);
        return {radius_m * std::cos(angle), radius_m * std::sin(angle),
                0.02 * static_cast<double>(node_index)};
    }

    const auto patch_config = findPatchConfigForNode(config, node_name);
    const double incidence_deg =
        patch_config && patch_config->patch_incidence_angle_deg ? *patch_config->patch_incidence_angle_deg : 45.0;
    const double flow_deg =
        patch_config && patch_config->patch_flow_angle_deg ? *patch_config->patch_flow_angle_deg : 0.0;
    const double incidence_rad = incidence_deg * kPi / 180.0;
    const double flow_rad = flow_deg * kPi / 180.0;
    return {
        radius_m * std::cos(flow_rad),
        radius_m * std::sin(flow_rad),
        radius_m * 0.35 * std::cos(incidence_rad),
    };
}

std::array<double, 3> pseudoNodeNormalForNode(const SurfaceChargingConfig& config,
                                              const SurfaceCircuitModel* circuit_model,
                                              std::size_t node_index)
{
    if (circuit_model == nullptr || node_index >= circuit_model->nodeCount())
    {
        return {1.0, 0.0, 0.0};
    }

    const auto node_name = circuit_model->nodeName(node_index);
    if (!circuit_model->nodeIsPatch(node_index))
    {
        return {0.0, 0.0, 1.0};
    }

    const auto patch_config = findPatchConfigForNode(config, node_name);
    const double incidence_deg =
        patch_config && patch_config->patch_incidence_angle_deg ? *patch_config->patch_incidence_angle_deg : 45.0;
    const double flow_deg =
        patch_config && patch_config->patch_flow_angle_deg ? *patch_config->patch_flow_angle_deg : 0.0;
    const double incidence_rad = incidence_deg * kPi / 180.0;
    const double flow_rad = flow_deg * kPi / 180.0;
    return {
        std::cos(flow_rad) * std::sin(incidence_rad),
        std::sin(flow_rad) * std::sin(incidence_rad),
        std::cos(incidence_rad),
    };
}

double boundaryFaceAreaForGroup(const SurfaceChargingConfig& config, const SurfaceCircuitModel* circuit_model,
                                const std::string& boundary_group_id)
{
    if (circuit_model == nullptr || boundary_group_id.empty())
    {
        return 0.0;
    }

    auto resolve_node_area = [&](const std::string& runtime_name) {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            if (circuit_model->nodeName(node_index) == runtime_name)
            {
                return std::max(0.0, circuit_model->nodeAreaM2(node_index));
            }
        }
        return 0.0;
    };

    for (const auto& group : config.body_boundary_groups)
    {
        if (group.id == boundary_group_id)
        {
            const auto area_m2 = resolve_node_area("body:" + group.body_id);
            return group.boundary_face_ids.empty() ? area_m2
                                                   : area_m2 / static_cast<double>(group.boundary_face_ids.size());
        }
    }
    for (const auto& group : config.patch_boundary_groups)
    {
        if (group.id == boundary_group_id)
        {
            const auto area_m2 = resolve_node_area("patch:" + group.patch_id);
            return group.boundary_face_ids.empty() ? area_m2
                                                   : area_m2 / static_cast<double>(group.boundary_face_ids.size());
        }
    }
    return 0.0;
}

} // namespace

bool writeExternalFieldSolveResultTemplateJson(const std::filesystem::path& json_path,
                                               const SurfaceCircuitModel* circuit_model)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.external_field_result.v1\",\n";
    output << "  \"nodes\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            output << "    {\"node_id\": \"" << jsonEscape(circuit_model->nodeName(node_index))
                   << "\", \"reference_potential_v\": 0.000000000000"
                   << ", \"normal_field_v_per_m\": 0.000000000000"
                   << ", \"local_charge_density_c_per_m3\": 0.000000000000"
                   << ", \"capacitance_scale\": 1.000000000000}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n}\n";
    return true;
}

bool writeExternalFieldSolveResultJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const std::vector<std::vector<double>>& node_field_solver_reference_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_charge_history,
    const std::vector<std::vector<double>>& node_field_solver_capacitance_scale_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.external_field_result.v1\",\n";
    output << "  \"nodes\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            output << "    {\"node_id\": \"" << jsonEscape(node_name)
                   << "\", \"logical_node_id\": \""
                   << jsonEscape(logicalNodeIdFromName(node_name))
                   << "\", \"node_type\": \""
                   << (circuit_model->nodeIsPatch(node_index) ? "patch" : "body")
                   << "\", \"boundary_group_id\": \""
                   << jsonEscape(boundaryGroupIdForNode(config, node_name))
                   << "\", \"reference_potential_v\": "
                   << latestHistoryValue(node_field_solver_reference_history, node_index)
                   << ", \"normal_field_v_per_m\": "
                   << latestHistoryValue(node_field_history, node_index)
                   << ", \"local_charge_density_c_per_m3\": "
                   << latestHistoryValue(node_charge_history, node_index)
                   << ", \"capacitance_scale\": "
                   << latestHistoryValue(node_field_solver_capacitance_scale_history, node_index)
                   << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n}\n";
    return true;
}

bool writeExternalFieldSolveRequestJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const SurfaceGraphCapacitanceMatrixProvider* matrix_provider,
    const SurfaceFieldSolverAdapter* field_solver_adapter,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_total_current_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_charge_history,
    const std::vector<std::vector<double>>& node_graph_capacitance_history,
    const std::vector<std::vector<double>>& node_graph_capacitance_row_sum_history,
    const std::vector<std::vector<double>>& node_field_solver_capacitance_scale_history,
    const std::vector<std::vector<double>>& branch_current_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.external_field_request.v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_strategy\": \""
           << surfacePicStrategyName(config.surface_pic_strategy) << "\",\n";
    output << "  \"live_pic_collision_cross_section_set_id\": \""
           << jsonEscape(config.live_pic_collision_cross_section_set_id) << "\",\n";
    output << "  \"surface_circuit_model\": \""
           << jsonEscape(circuit_model ? circuit_model->modelName() : "BuiltinCircuit") << "\",\n";
    output << "  \"matrix_family\": \""
           << jsonEscape(matrix_provider ? matrix_provider->matrixFamilyName() : "BuiltinMatrixFamily")
           << "\",\n";
    output << "  \"solver_adapter_hint\": \""
           << jsonEscape(matrix_provider ? matrix_provider->solverAdapterHint()
                                         : (field_solver_adapter ? field_solver_adapter->modelName()
                                                                 : "BuiltinFieldSolverAdapter"))
           << "\",\n";
    output << "  \"nodes\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            output << "    {\"node_id\": \"" << jsonEscape(circuit_model->nodeName(node_index))
                   << "\", \"logical_node_id\": \""
                   << jsonEscape(logicalNodeIdFromName(circuit_model->nodeName(node_index)))
                   << "\", \"node_type\": \""
                   << (circuit_model->nodeIsPatch(node_index) ? "patch" : "body")
                   << "\", \"potential_v\": "
                   << latestHistoryValue(node_potential_history, node_index)
                   << ", \"total_current_density_a_per_m2\": "
                   << latestHistoryValue(node_total_current_history, node_index)
                   << ", \"normal_field_v_per_m\": "
                   << latestHistoryValue(node_field_history, node_index)
                   << ", \"local_charge_density_c_per_m3\": "
                   << latestHistoryValue(node_charge_history, node_index)
                   << ", \"graph_capacitance_diagonal_f\": "
                   << latestHistoryValue(node_graph_capacitance_history, node_index)
                   << ", \"graph_capacitance_row_sum_f\": "
                   << latestHistoryValue(node_graph_capacitance_row_sum_history, node_index)
                   << ", \"field_solver_capacitance_scale\": "
                   << latestHistoryValue(node_field_solver_capacitance_scale_history, node_index)
                   << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"interfaces\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            output << "    {\"interface_id\": \"" << jsonEscape(circuit_model->branchName(branch_index))
                   << "\", \"from_node\": \"" << jsonEscape(circuit_model->nodeName(from_node))
                   << "\", \"to_node\": \"" << jsonEscape(circuit_model->nodeName(to_node))
                   << "\", \"conductance_s\": " << circuit_model->branchConductanceS(branch_index)
                   << ", \"current_a\": " << latestHistoryValue(branch_current_history, branch_index)
                   << ", \"mutual_capacitance_f\": "
                   << (matrix_provider ? matrix_provider->mutualCapacitanceF(config, circuit_model, branch_index)
                                       : 0.0)
                   << "}";
            output << (branch_index + 1 < circuit_model->branchCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"body_boundary_groups\": [\n";
    for (std::size_t i = 0; i < config.body_boundary_groups.size(); ++i)
    {
        const auto& group = config.body_boundary_groups[i];
        output << "    {\"id\": \"" << jsonEscape(group.id) << "\", \"body_id\": \""
               << jsonEscape(group.body_id) << "\", \"external_group_name\": \""
               << jsonEscape(group.external_group_name) << "\", \"boundary_face_count\": "
               << group.boundary_face_ids.size() << "}";
        output << (i + 1 < config.body_boundary_groups.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"patch_boundary_groups\": [\n";
    for (std::size_t i = 0; i < config.patch_boundary_groups.size(); ++i)
    {
        const auto& group = config.patch_boundary_groups[i];
        output << "    {\"id\": \"" << jsonEscape(group.id) << "\", \"patch_id\": \""
               << jsonEscape(group.patch_id) << "\", \"external_group_name\": \""
               << jsonEscape(group.external_group_name) << "\", \"boundary_face_count\": "
               << group.boundary_face_ids.size() << "}";
        output << (i + 1 < config.patch_boundary_groups.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"boundary_mappings\": [\n";
    for (std::size_t i = 0; i < config.boundary_mappings.size(); ++i)
    {
        const auto& mapping = config.boundary_mappings[i];
        output << "    {\"node_id\": \"" << jsonEscape(mapping.node_id)
               << "\", \"boundary_group_id\": \"" << jsonEscape(mapping.boundary_group_id)
               << "\", \"required\": " << (mapping.required ? "true" : "false") << "}";
        output << (i + 1 < config.boundary_mappings.size() ? ",\n" : "\n");
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeSurfaceVolumeProjectionJson(const std::filesystem::path& json_path,
                                      const SurfaceChargingConfig& config,
                                      const SurfaceCircuitModel* circuit_model,
                                      const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_volume_projection.v1\",\n";
    output << "  \"projection_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->projectionFamilyName()
                                                   : "NodeToPseudoCellProjection")
           << "\",\n";
    output << "  \"mesh_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->meshFamilyName()
                                                   : "PseudoBoundaryVolumeMesh")
           << "\",\n";
    output << "  \"supports_projection_weights\": "
           << ((volumetric_solver_adapter && volumetric_solver_adapter->supportsProjectionWeights()) ? "true"
                                                                                                      : "false")
           << ",\n";
    output << "  \"surface_to_volume\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const double projection_weight = projectionWeightForNode(config, node_name);
            output << "    {\"node_id\": \"" << jsonEscape(node_name) << "\", \"cell_id\": \"cell_"
                   << node_index << "\", \"boundary_group_id\": \"" << jsonEscape(boundaryGroupIdForNode(config, node_name)) << "\", \"weight\": "
                   << projection_weight << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"projection_weights\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const double projection_weight = projectionWeightForNode(config, node_name);
            output << "    {\"node_id\": \"" << jsonEscape(node_name)
                   << "\", \"cell_id\": \"cell_" << node_index
                   << "\", \"boundary_group_id\": \"" << jsonEscape(boundaryGroupIdForNode(config, node_name))
                   << "\", \"surface_to_volume_weight\": " << projection_weight
                   << ", \"volume_to_surface_weight\": " << projection_weight << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"volume_to_surface\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const double projection_weight = projectionWeightForNode(config, node_name);
            output << "    {\"cell_id\": \"cell_" << node_index << "\", \"node_id\": \""
                   << jsonEscape(node_name) << "\", \"boundary_group_id\": \"" << jsonEscape(boundaryGroupIdForNode(config, node_name)) << "\", \"weight\": "
                   << projection_weight << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeExternalVolumeMeshStubJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_field_solver_reference_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_charge_history,
    const std::vector<std::vector<double>>& node_field_solver_capacitance_scale_history,
    const std::vector<std::vector<double>>& node_volume_potential_history,
    const std::vector<std::vector<double>>& node_deposited_charge_history,
    const std::vector<std::vector<double>>& node_poisson_residual_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const double sheath_length_m =
        std::clamp(config.sheath_length_m > 0.0 ? config.sheath_length_m : 1.0e-3,
                   std::max(1.0e-6, config.minimum_sheath_length_m),
                   std::max(config.minimum_sheath_length_m, config.maximum_sheath_length_m));

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_volume_stub.v1\",\n";
    output << "  \"pseudo_sheath_length_m\": " << sheath_length_m << ",\n";
    output << "  \"cells\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const auto boundary_group_id = boundaryGroupIdForNode(config, node_name);
            const double area_m2 = std::max(0.0, circuit_model->nodeAreaM2(node_index));
            const double projection_weight = projectionWeightForNode(config, node_name);
            const auto center = pseudoNodeCenterForNode(config, circuit_model, node_index);
            output << "    {\"cell_id\": \"cell_" << node_index << "\", \"node_id\": \""
                   << jsonEscape(node_name) << "\", \"logical_node_id\": \""
                   << jsonEscape(logicalNodeIdFromName(node_name))
                   << "\", \"node_type\": \"" << (circuit_model->nodeIsPatch(node_index) ? "patch" : "body")
                   << "\", \"boundary_group_id\": \"" << jsonEscape(boundary_group_id)
                   << "\", \"pseudo_volume_m3\": " << (area_m2 * sheath_length_m * projection_weight)
                   << ", \"center_x_m\": " << center[0]
                   << ", \"center_y_m\": " << center[1]
                   << ", \"center_z_m\": " << center[2]
                   << ", \"characteristic_length_m\": " << std::sqrt(std::max(0.0, area_m2) / kPi)
                   << ", \"potential_v\": " << latestHistoryValue(node_potential_history, node_index)
                   << ", \"reference_potential_v\": "
                   << latestHistoryValue(node_field_solver_reference_history, node_index)
                   << ", \"normal_field_v_per_m\": "
                   << latestHistoryValue(node_field_history, node_index)
                   << ", \"local_charge_density_c_per_m3\": "
                   << latestHistoryValue(node_charge_history, node_index)
                   << ", \"volume_potential_v\": "
                   << latestHistoryValue(node_volume_potential_history, node_index)
                   << ", \"deposited_charge_c\": "
                   << latestHistoryValue(node_deposited_charge_history, node_index)
                   << ", \"poisson_residual_v_m\": "
                   << latestHistoryValue(node_poisson_residual_history, node_index)
                   << ", \"field_solver_capacitance_scale\": "
                   << latestHistoryValue(node_field_solver_capacitance_scale_history, node_index)
                   << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeVolumeMeshSkeletonJson(const std::filesystem::path& json_path,
                                 const SurfaceChargingConfig& config,
                                 const SurfaceCircuitModel* circuit_model,
                                 const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.volume_mesh_stub.v1\",\n";
    output << "  \"mesh_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->meshFamilyName()
                                                   : "PseudoBoundaryVolumeMesh")
           << "\",\n";
    output << "  \"projection_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->projectionFamilyName()
                                                   : "NodeToPseudoCellProjection")
           << "\",\n";
    output << "  \"supports_boundary_face_mapping\": "
           << ((volumetric_solver_adapter && volumetric_solver_adapter->supportsBoundaryFaceMapping()) ? "true"
                                                                                                        : "false")
           << ",\n";
    output << "  \"boundary_faces\": [\n";

    bool first_face = true;
    for (const auto& group : config.body_boundary_groups)
    {
        const auto runtime_name = "body:" + group.body_id;
        const double face_area_m2 = boundaryFaceAreaForGroup(config, circuit_model, group.id);
        for (std::size_t local_face_index = 0; local_face_index < group.boundary_face_ids.size();
             ++local_face_index)
        {
            const int face_id = group.boundary_face_ids[local_face_index];
            std::size_t node_index = 0;
            if (circuit_model != nullptr)
            {
                for (std::size_t candidate = 0; candidate < circuit_model->nodeCount(); ++candidate)
                {
                    if (circuit_model->nodeName(candidate) == runtime_name)
                    {
                        node_index = candidate;
                        break;
                    }
                }
            }
            const auto center = pseudoNodeCenterForNode(config, circuit_model, node_index);
            const auto normal = pseudoNodeNormalForNode(config, circuit_model, node_index);
            if (!first_face)
            {
                output << ",\n";
            }
            first_face = false;
            output << "    {\"face_id\": \"face_" << jsonEscape(group.id) << "_" << local_face_index
                   << "\", \"boundary_group_id\": \"" << jsonEscape(group.id)
                   << "\", \"node_id\": \"" << jsonEscape(runtime_name)
                   << "\", \"external_face_id\": " << face_id
                   << ", \"area_m2\": " << face_area_m2
                   << ", \"center_x_m\": " << center[0]
                   << ", \"center_y_m\": " << center[1]
                   << ", \"center_z_m\": " << center[2]
                   << ", \"normal_x\": " << normal[0]
                   << ", \"normal_y\": " << normal[1]
                   << ", \"normal_z\": " << normal[2] << "}";
        }
    }
    for (const auto& group : config.patch_boundary_groups)
    {
        const auto runtime_name = "patch:" + group.patch_id;
        const double face_area_m2 = boundaryFaceAreaForGroup(config, circuit_model, group.id);
        for (std::size_t local_face_index = 0; local_face_index < group.boundary_face_ids.size();
             ++local_face_index)
        {
            const int face_id = group.boundary_face_ids[local_face_index];
            std::size_t node_index = 0;
            if (circuit_model != nullptr)
            {
                for (std::size_t candidate = 0; candidate < circuit_model->nodeCount(); ++candidate)
                {
                    if (circuit_model->nodeName(candidate) == runtime_name)
                    {
                        node_index = candidate;
                        break;
                    }
                }
            }
            const auto center = pseudoNodeCenterForNode(config, circuit_model, node_index);
            const auto normal = pseudoNodeNormalForNode(config, circuit_model, node_index);
            if (!first_face)
            {
                output << ",\n";
            }
            first_face = false;
            output << "    {\"face_id\": \"face_" << jsonEscape(group.id) << "_" << local_face_index
                   << "\", \"boundary_group_id\": \"" << jsonEscape(group.id)
                   << "\", \"node_id\": \"" << jsonEscape(runtime_name)
                   << "\", \"external_face_id\": " << face_id
                   << ", \"area_m2\": " << face_area_m2
                   << ", \"center_x_m\": " << center[0]
                   << ", \"center_y_m\": " << center[1]
                   << ", \"center_z_m\": " << center[2]
                   << ", \"normal_x\": " << normal[0]
                   << ", \"normal_y\": " << normal[1]
                   << ", \"normal_z\": " << normal[2] << "}";
        }
    }
    if (!first_face)
    {
        output << "\n";
    }
    output << "  ],\n";
    output << "  \"cells\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const auto boundary_group_id = boundaryGroupIdForNode(config, node_name);
            const auto center = pseudoNodeCenterForNode(config, circuit_model, node_index);
            const auto normal = pseudoNodeNormalForNode(config, circuit_model, node_index);
            std::size_t linked_face_count = 0;
            for (const auto& group : config.body_boundary_groups)
            {
                if (group.id == boundary_group_id)
                {
                    linked_face_count = group.boundary_face_ids.size();
                }
            }
            for (const auto& group : config.patch_boundary_groups)
            {
                if (group.id == boundary_group_id)
                {
                    linked_face_count = group.boundary_face_ids.size();
                }
            }
            output << "    {\"cell_id\": \"cell_" << node_index
                   << "\", \"node_id\": \"" << jsonEscape(node_name)
                   << "\", \"boundary_group_id\": \"" << jsonEscape(boundary_group_id)
                   << "\", \"linked_face_count\": " << linked_face_count
                   << ", \"node_area_m2\": " << std::max(0.0, circuit_model->nodeAreaM2(node_index))
                   << ", \"center_x_m\": " << center[0]
                   << ", \"center_y_m\": " << center[1]
                   << ", \"center_z_m\": " << center[2]
                   << ", \"surface_normal_x\": " << normal[0]
                   << ", \"surface_normal_y\": " << normal[1]
                   << ", \"surface_normal_z\": " << normal[2]
                   << ", \"cell_volume_m3\": "
                   << (std::max(0.0, circuit_model->nodeAreaM2(node_index)) *
                       std::max(1.0e-6, config.sheath_length_m))
                   << ", \"characteristic_length_m\": "
                   << std::sqrt(std::max(0.0, circuit_model->nodeAreaM2(node_index)) / kPi)
                   << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"cell_neighbors\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            const auto from_center = pseudoNodeCenterForNode(config, circuit_model, from_node);
            const auto to_center = pseudoNodeCenterForNode(config, circuit_model, to_node);
            const auto face_center = std::array<double, 3>{
                0.5 * (from_center[0] + to_center[0]),
                0.5 * (from_center[1] + to_center[1]),
                0.5 * (from_center[2] + to_center[2])};
            output << "    {\"source_cell_id\": \"cell_" << from_node
                   << "\", \"target_cell_id\": \"cell_" << to_node
                   << "\", \"interface_id\": \"" << jsonEscape(circuit_model->branchName(branch_index))
                   << "\", \"conductance_s\": " << circuit_model->branchConductanceS(branch_index)
                   << ", \"resistance_ohm\": " << circuit_model->branchResistanceOhm(branch_index)
                   << ", \"shared_face_area_m2\": "
                   << (0.5 * std::min(std::max(0.0, circuit_model->nodeAreaM2(from_node)),
                                      std::max(0.0, circuit_model->nodeAreaM2(to_node))))
                   << ", \"face_distance_m\": "
                   << std::max(1.0e-6,
                               0.5 * (std::sqrt(std::max(0.0, circuit_model->nodeAreaM2(from_node)) / kPi) +
                                      std::sqrt(std::max(0.0, circuit_model->nodeAreaM2(to_node)) / kPi)) +
                                   std::max(1.0e-6, config.sheath_length_m))
                   << ", \"face_center_x_m\": " << face_center[0]
                   << ", \"face_center_y_m\": " << face_center[1]
                   << ", \"face_center_z_m\": " << face_center[2]
                   << ", \"permittivity_scale\": 1.0"
                   << ", \"branch_type\": \"" << jsonEscape(branchTypeLabel(circuit_model, branch_index))
                   << "\"}";
            output << (branch_index + 1 < circuit_model->branchCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"cell_faces\": [\n";
    if (circuit_model != nullptr)
    {
        bool first_cell_face = true;
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const auto boundary_group_id = boundaryGroupIdForNode(config, node_name);
            std::vector<int> boundary_face_ids;
            for (const auto& group : config.body_boundary_groups)
            {
                if (group.id == boundary_group_id)
                {
                    boundary_face_ids = group.boundary_face_ids;
                }
            }
            for (const auto& group : config.patch_boundary_groups)
            {
                if (group.id == boundary_group_id)
                {
                    boundary_face_ids = group.boundary_face_ids;
                }
            }
            const double projection_weight = projectionWeightForNode(config, node_name);
            for (std::size_t local_face_index = 0; local_face_index < boundary_face_ids.size();
                 ++local_face_index)
            {
                if (!first_cell_face)
                {
                    output << ",\n";
                }
                first_cell_face = false;
                output << "    {\"cell_id\": \"cell_" << node_index
                       << "\", \"face_id\": \"face_" << jsonEscape(boundary_group_id) << "_"
                       << local_face_index
                       << "\", \"role\": \"boundary\", \"projection_weight\": "
                       << projection_weight << "}";
            }
        }
        if (!first_cell_face)
        {
            output << "\n";
        }
    }
    output << "  ],\n";
    output << "  \"face_neighbors\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            const auto from_group = boundaryGroupIdForNode(config, circuit_model->nodeName(from_node));
            const auto to_group = boundaryGroupIdForNode(config, circuit_model->nodeName(to_node));
            const auto from_center = pseudoNodeCenterForNode(config, circuit_model, from_node);
            const auto to_center = pseudoNodeCenterForNode(config, circuit_model, to_node);
            output << "    {\"source_boundary_group_id\": \"" << jsonEscape(from_group)
                   << "\", \"target_boundary_group_id\": \"" << jsonEscape(to_group)
                   << "\", \"interface_id\": \"" << jsonEscape(circuit_model->branchName(branch_index))
                   << "\", \"center_distance_m\": "
                   << std::sqrt((from_center[0] - to_center[0]) * (from_center[0] - to_center[0]) +
                                (from_center[1] - to_center[1]) * (from_center[1] - to_center[1]) +
                                (from_center[2] - to_center[2]) * (from_center[2] - to_center[2]))
                   << ", \"branch_type\": \"" << jsonEscape(branchTypeLabel(circuit_model, branch_index))
                   << "\"}";
            output << (branch_index + 1 < circuit_model->branchCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeExternalFieldBridgeManifestJson(const std::filesystem::path& json_path,
                                          const std::filesystem::path& csv_path,
                                          const SurfaceChargingConfig& config,
                                          const SurfaceCircuitModel* circuit_model,
                                          const SurfaceFieldSolverAdapter* field_solver_adapter,
                                          const SurfaceGraphCapacitanceMatrixProvider* matrix_provider)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    auto artifact = [&](const char* extension) {
        auto path = csv_path;
        path.replace_extension(extension);
        return path;
    };

    output << "{\n";
    output << "  \"schema_version\": \"scdat.field_bridge_manifest.v1\",\n";
    output << "  \"bridge_enabled\": " << (config.enable_external_field_solver_bridge ? "true" : "false")
           << ",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_strategy\": \"" << surfacePicStrategyName(config.surface_pic_strategy)
           << "\",\n";
    output << "  \"surface_field_solver_adapter\": \""
           << jsonEscape(field_solver_adapter ? field_solver_adapter->modelName()
                                              : "BuiltinFieldSolverAdapter")
           << "\",\n";
    output << "  \"surface_graph_capacitance_matrix_family\": \""
           << jsonEscape(matrix_provider ? matrix_provider->matrixFamilyName() : "BuiltinMatrixFamily")
           << "\",\n";
    output << "  \"node_count\": " << (circuit_model ? circuit_model->nodeCount() : 0) << ",\n";
    output << "  \"branch_count\": " << (circuit_model ? circuit_model->branchCount() : 0) << ",\n";
    output << "  \"artifacts\": {\n";
    output << "    \"boundary_mapping_json\": \"" << jsonEscape(artifact(".boundary_mapping.json").string())
           << "\",\n";
    output << "    \"field_request_json\": \"" << jsonEscape(artifact(".field_request.json").string())
           << "\",\n";
    output << "    \"field_result_template_json\": \""
           << jsonEscape(artifact(".field_result_template.json").string()) << "\",\n";
    output << "    \"field_result_json\": \"" << jsonEscape(artifact(".field_result.json").string())
           << "\",\n";
    output << "    \"graph_matrix_json\": \"" << jsonEscape(artifact(".graph_matrix.json").string())
           << "\",\n";
    output << "    \"field_adapter_json\": \"" << jsonEscape(artifact(".field_adapter.json").string())
           << "\",\n";
    output << "    \"volume_stub_json\": \"" << jsonEscape(artifact(".volume_stub.json").string())
           << "\",\n";
    output << "    \"volume_mesh_stub_json\": \""
           << jsonEscape(artifact(".volume_mesh_stub.json").string()) << "\"\n";
    output << "  }\n";
    output << "}\n";
    return true;
}

bool writeVolumeHistoryJson(
    const std::filesystem::path& json_path, const SurfaceCircuitModel* circuit_model,
    const std::vector<double>& time_history,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_volume_potential_history,
    const std::vector<std::vector<double>>& node_deposited_charge_history,
    const std::vector<std::vector<double>>& node_poisson_residual_history,
    const std::vector<std::vector<double>>& node_pseudo_volume_history,
    const std::vector<std::vector<double>>& node_projection_weight_history,
    const std::vector<std::vector<double>>& node_mesh_coupling_gain_history,
    const std::vector<std::vector<double>>& node_solver_mode_history,
    const std::vector<std::vector<double>>& node_solver_iterations_history,
    const std::vector<std::vector<double>>& node_solver_linear_iterations_history,
    const std::vector<std::vector<double>>& node_solver_converged_history,
    const std::vector<std::vector<double>>& node_solver_residual_history,
    const std::vector<std::vector<double>>& node_solver_max_delta_history,
    const std::vector<std::vector<double>>& node_solver_matrix_nnz_history,
    const std::vector<std::vector<double>>& node_solver_cell_count_history,
    const std::vector<std::vector<double>>& node_field_volume_coupling_iterations_history,
    const std::vector<std::vector<double>>& node_field_volume_coupling_converged_history,
    const std::vector<std::vector<double>>& node_field_volume_coupling_max_delta_history,
    const std::vector<std::vector<double>>& node_field_volume_coupling_relaxation_used_history,
    const std::vector<std::vector<double>>& node_external_volume_feedback_blend_factor_history,
    const std::vector<std::vector<double>>& node_external_volume_feedback_mismatch_metric_history,
    const std::vector<std::vector<double>>& node_external_volume_feedback_applied_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.volume_history.v1\",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"times_s\": [";
    for (std::size_t i = 0; i < time_history.size(); ++i)
    {
        if (i > 0)
        {
            output << ", ";
        }
        output << time_history[i];
    }
    output << "],\n";
    output << "  \"nodes\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            if (node_index > 0)
            {
                output << ",\n";
            }
            auto write_series = [&](const char* key, const std::vector<std::vector<double>>& history) {
                output << "      \"" << key << "\": [";
                const auto& series = node_index < history.size() ? history[node_index]
                                                                 : std::vector<double>{};
                for (std::size_t i = 0; i < series.size(); ++i)
                {
                    if (i > 0)
                    {
                        output << ", ";
                    }
                    output << series[i];
                }
                output << "]";
            };

            output << "    {\n";
            output << "      \"node_index\": " << node_index << ",\n";
            output << "      \"node_name\": \"" << jsonEscape(circuit_model->nodeName(node_index))
                   << "\",\n";
            output << "      \"is_patch\": "
                   << (circuit_model->nodeIsPatch(node_index) ? "true" : "false") << ",\n";
            write_series("surface_potential_v", node_potential_history);
            output << ",\n";
            write_series("volume_potential_v", node_volume_potential_history);
            output << ",\n";
            write_series("deposited_charge_c", node_deposited_charge_history);
            output << ",\n";
            write_series("poisson_residual_v_m", node_poisson_residual_history);
            output << ",\n";
            write_series("pseudo_volume_m3", node_pseudo_volume_history);
            output << ",\n";
            write_series("projection_weight_sum", node_projection_weight_history);
            output << ",\n";
            write_series("mesh_coupling_gain", node_mesh_coupling_gain_history);
            output << ",\n";
            write_series("solver_mode_id", node_solver_mode_history);
            output << ",\n";
            write_series("solver_iterations", node_solver_iterations_history);
            output << ",\n";
            write_series("solver_linear_iterations", node_solver_linear_iterations_history);
            output << ",\n";
            write_series("solver_converged", node_solver_converged_history);
            output << ",\n";
            write_series("solver_residual_norm", node_solver_residual_history);
            output << ",\n";
            write_series("solver_max_delta_v", node_solver_max_delta_history);
            output << ",\n";
            write_series("solver_matrix_nnz", node_solver_matrix_nnz_history);
            output << ",\n";
            write_series("solver_cell_count", node_solver_cell_count_history);
            output << ",\n";
            write_series("field_volume_coupling_iterations",
                         node_field_volume_coupling_iterations_history);
            output << ",\n";
            write_series("field_volume_coupling_converged",
                         node_field_volume_coupling_converged_history);
            output << ",\n";
            write_series("field_volume_coupling_max_delta",
                         node_field_volume_coupling_max_delta_history);
            output << ",\n";
            write_series("field_volume_coupling_relaxation_used",
                         node_field_volume_coupling_relaxation_used_history);
            output << ",\n";
            write_series("external_volume_feedback_blend_factor",
                         node_external_volume_feedback_blend_factor_history);
            output << ",\n";
            write_series("external_volume_feedback_mismatch_metric",
                         node_external_volume_feedback_mismatch_metric_history);
            output << ",\n";
            write_series("external_volume_feedback_applied",
                         node_external_volume_feedback_applied_history);
            output << "\n    }";
        }
        output << "\n";
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeVolumetricSolverAdapterContractJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter,
    const SurfaceCircuitModel* circuit_model)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << "{\n";
    output << "  \"schema_version\": \"scdat.volumetric_solver_adapter_contract.v1\",\n";
    output << "  \"volumetric_solver_adapter\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->modelName()
                                                   : "BuiltinVolumetricAdapter")
           << "\",\n";
    output << "  \"projection_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->projectionFamilyName()
                                                   : "NodeToPseudoCellProjection")
           << "\",\n";
    output << "  \"mesh_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->meshFamilyName()
                                                   : "PseudoBoundaryVolumeMesh")
           << "\",\n";
    output << "  \"request_schema_version\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->requestSchemaVersion()
                                                   : "scdat.external_volume_request.v1")
           << "\",\n";
    output << "  \"result_schema_version\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->resultSchemaVersion()
                                                   : "scdat.external_volume_result.v1")
           << "\",\n";
    output << "  \"supports_surface_to_volume_projection\": true,\n";
    output << "  \"supports_volume_to_surface_projection\": true,\n";
    output << "  \"supports_boundary_face_mapping\": "
           << ((volumetric_solver_adapter && volumetric_solver_adapter->supportsBoundaryFaceMapping())
                   ? "true"
                   : "false")
           << ",\n";
    output << "  \"supports_projection_weights\": "
           << ((volumetric_solver_adapter && volumetric_solver_adapter->supportsProjectionWeights())
                   ? "true"
                   : "false")
           << ",\n";
    output << "  \"supports_cell_centers\": true,\n";
    output << "  \"supports_face_centers\": true,\n";
    output << "  \"supports_face_normals\": true,\n";
    output << "  \"supports_neighbor_face_geometry\": true,\n";
    output << "  \"supports_cell_face_links\": true,\n";
    output << "  \"supports_cell_volume_override\": true,\n";
    output << "  \"supports_cell_initial_state\": true,\n";
    output << "  \"solver_mode_hint\": \"" << volumeLinearSolverPolicyName(config.volume_linear_solver_policy)
           << "\",\n";
    output << "  \"node_count\": " << (circuit_model ? circuit_model->nodeCount() : 0) << ",\n";
    output << "  \"branch_count\": " << (circuit_model ? circuit_model->branchCount() : 0) << "\n";
    output << "}\n";
    return true;
}

bool writeExternalVolumeSolveRequestJson(const std::filesystem::path& json_path,
                                         const std::filesystem::path& csv_path,
                                         const SurfaceChargingConfig& config,
                                         const SurfaceCircuitModel* circuit_model,
                                         const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    auto artifact = [&](const char* extension) {
        auto path = csv_path;
        path.replace_extension(extension);
        return path;
    };

    output << "{\n";
    output << "  \"schema_version\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->requestSchemaVersion()
                                                   : "scdat.external_volume_request.v1")
           << "\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_strategy\": \"" << surfacePicStrategyName(config.surface_pic_strategy)
           << "\",\n";
    output << "  \"volumetric_solver_adapter\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->modelName()
                                                   : "BuiltinVolumetricAdapter")
           << "\",\n";
    output << "  \"projection_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->projectionFamilyName()
                                                   : "NodeToPseudoCellProjection")
           << "\",\n";
    output << "  \"mesh_family\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->meshFamilyName()
                                                   : "PseudoBoundaryVolumeMesh")
           << "\",\n";
    output << "  \"node_count\": " << (circuit_model ? circuit_model->nodeCount() : 0) << ",\n";
    output << "  \"branch_count\": " << (circuit_model ? circuit_model->branchCount() : 0) << ",\n";
    output << "  \"expects_cell_centers\": true,\n";
    output << "  \"expects_face_centers\": true,\n";
    output << "  \"expects_face_normals\": true,\n";
    output << "  \"expects_neighbor_face_geometry\": true,\n";
    output << "  \"expects_cell_face_links\": true,\n";
    output << "  \"expects_cell_volume_override\": true,\n";
    output << "  \"expects_cell_initial_state\": true,\n";
    output << "  \"solver_mode_hint\": \"" << volumeLinearSolverPolicyName(config.volume_linear_solver_policy)
           << "\",\n";
    output << "  \"linear_max_iterations\": "
           << std::clamp<std::size_t>(config.volume_linear_max_iterations, 8, 4096) << ",\n";
    output << "  \"linear_tolerance_scale\": "
           << std::clamp(config.volume_linear_tolerance_scale, 0.01, 2.0) << ",\n";
    output << "  \"linear_relaxation\": "
           << std::clamp(config.volume_linear_relaxation, 0.2, 1.0) << ",\n";
    output << "  \"field_volume_outer_iterations\": "
           << std::clamp<std::size_t>(config.field_volume_outer_iterations, 1, 8) << ",\n";
    output << "  \"field_volume_outer_tolerance\": "
           << std::clamp(config.field_volume_outer_tolerance, 1.0e-9, 10.0) << ",\n";
    output << "  \"field_volume_outer_relaxation\": "
           << std::clamp(config.field_volume_outer_relaxation, 0.05, 1.0) << ",\n";
    output << "  \"self_consistent_iterations\": "
           << std::clamp<std::size_t>(config.volume_self_consistent_iterations, 1, 8) << ",\n";
    output << "  \"self_consistent_tolerance_v\": "
           << std::clamp(config.volume_self_consistent_tolerance_v, 1.0e-9, 1.0) << ",\n";
    output << "  \"charge_relaxation\": "
           << std::clamp(config.volume_charge_relaxation, 0.05, 1.0) << ",\n";
    output << "  \"potential_relaxation\": "
           << std::clamp(config.volume_potential_relaxation, 0.05, 1.0) << ",\n";
    output << "  \"artifacts\": {\n";
    output << "    \"volume_stub_json\": \"" << jsonEscape(artifact(".volume_stub.json").string())
           << "\",\n";
    output << "    \"volume_mesh_stub_json\": \"" << jsonEscape(artifact(".volume_mesh_stub.json").string())
           << "\",\n";
    output << "    \"surface_volume_projection_json\": \""
           << jsonEscape(artifact(".surface_volume_projection.json").string()) << "\",\n";
    output << "    \"graph_matrix_json\": \"" << jsonEscape(artifact(".graph_matrix.json").string())
           << "\",\n";
    output << "    \"field_request_json\": \"" << jsonEscape(artifact(".field_request.json").string())
           << "\",\n";
    output << "    \"volume_result_json\": \"" << jsonEscape(artifact(".volume_result.json").string())
           << "\"\n";
    output << "  }\n";
    output << "}\n";
    return true;
}

bool writeExternalVolumeSolveResultTemplateJson(const std::filesystem::path& json_path,
                                                const SurfaceCircuitModel* circuit_model,
                                                const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->resultSchemaVersion()
                                                   : "scdat.external_volume_result.v1")
           << "\",\n";
    output << "  \"cells\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            output << "    {\"cell_id\": \"cell_" << node_index
                   << "\", \"potential_v\": 0.000000000000"
                   << ", \"reference_potential_v\": 0.000000000000"
                   << ", \"normal_field_v_per_m\": 0.000000000000"
                   << ", \"local_charge_density_c_per_m3\": 0.000000000000"
                   << ", \"capacitance_scale\": 1.000000000000"
                   << ", \"coupling_gain\": 0.000000000000"
                   << ", \"projection_weight_scale\": 1.000000000000"
                   << ", \"sheath_length_scale\": 1.000000000000}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeExternalVolumeSolveResultJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_volume_potential_history,
    const std::vector<std::vector<double>>& node_field_solver_reference_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_charge_history,
    const std::vector<std::vector<double>>& node_field_solver_capacitance_scale_history,
    const std::vector<std::vector<double>>& node_projection_weight_history,
    const std::vector<std::vector<double>>& node_mesh_coupling_gain_history,
    const std::vector<std::vector<double>>& node_effective_sheath_length_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const double nominal_sheath_length_m =
        std::clamp(config.sheath_length_m > 0.0 ? config.sheath_length_m : 1.0e-3,
                   std::max(1.0e-6, config.minimum_sheath_length_m),
                   std::max(config.minimum_sheath_length_m, config.maximum_sheath_length_m));

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \""
           << jsonEscape(volumetric_solver_adapter ? volumetric_solver_adapter->resultSchemaVersion()
                                                   : "scdat.external_volume_result.v1")
           << "\",\n";
    output << "  \"cells\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            const auto node_name = circuit_model->nodeName(node_index);
            const double effective_sheath_length_m =
                std::max(0.0, latestHistoryValue(node_effective_sheath_length_history, node_index));
            const double sheath_length_scale =
                effective_sheath_length_m > 0.0
                    ? effective_sheath_length_m / std::max(1.0e-12, nominal_sheath_length_m)
                    : 1.0;
            double volume_potential_v = latestHistoryValue(node_volume_potential_history, node_index);
            if (!std::isfinite(volume_potential_v))
            {
                volume_potential_v = latestHistoryValue(node_potential_history, node_index);
            }

            output << "    {\"cell_id\": \"cell_" << node_index
                   << "\", \"node_id\": \"" << jsonEscape(node_name)
                   << "\", \"logical_node_id\": \""
                   << jsonEscape(logicalNodeIdFromName(node_name))
                   << "\", \"node_type\": \"" << (circuit_model->nodeIsPatch(node_index) ? "patch" : "body")
                   << "\", \"boundary_group_id\": \"" << jsonEscape(boundaryGroupIdForNode(config, node_name))
                   << "\", \"potential_v\": " << volume_potential_v
                   << ", \"reference_potential_v\": "
                   << latestHistoryValue(node_field_solver_reference_history, node_index)
                   << ", \"normal_field_v_per_m\": "
                   << latestHistoryValue(node_field_history, node_index)
                   << ", \"local_charge_density_c_per_m3\": "
                   << latestHistoryValue(node_charge_history, node_index)
                   << ", \"capacitance_scale\": "
                   << latestHistoryValue(node_field_solver_capacitance_scale_history, node_index)
                   << ", \"coupling_gain\": "
                   << latestHistoryValue(node_mesh_coupling_gain_history, node_index)
                   << ", \"projection_weight_scale\": "
                   << std::max(1.0, latestHistoryValue(node_projection_weight_history, node_index))
                   << ", \"sheath_length_scale\": " << std::max(1.0e-12, sheath_length_scale)
                   << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeSurfaceGraphMatrixSnapshotCsv(
    const std::filesystem::path& csv_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider)
{
    std::ofstream output(csv_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "entry_type,index,name,from_node,to_node,value_f,aux_value,notes\n";
    if (circuit_model == nullptr || graph_capacitance_matrix_provider == nullptr)
    {
        return true;
    }

    for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
    {
        output << "node_diagonal," << node_index << "," << circuit_model->nodeName(node_index) << ",,,"
               << graph_capacitance_matrix_provider->diagonalCapacitanceF(config, circuit_model, node_index)
               << ","
               << graph_capacitance_matrix_provider->rowSumCapacitanceF(config, circuit_model, node_index)
               << ",row_sum_f\n";
    }
    for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
    {
        const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
        const auto to_node = circuit_model->branchToNodeIndex(branch_index);
        output << "branch_mutual," << branch_index << "," << circuit_model->branchName(branch_index)
               << "," << circuit_model->nodeName(from_node) << "," << circuit_model->nodeName(to_node)
               << ","
               << graph_capacitance_matrix_provider->mutualCapacitanceF(config, circuit_model, branch_index)
               << "," << circuit_model->branchConductanceS(branch_index) << ",conductance_s\n";
    }
    output << "graph_metric,0,matrix_family,,,0.000000000000,0.000000000000,"
           << graph_capacitance_matrix_provider->matrixFamilyName() << "\n";
    output << "graph_metric,1,graph_coupling_metric,,,"
           << graph_capacitance_matrix_provider->graphCouplingMetric(config, circuit_model)
           << ",0.000000000000,graph_coupling_metric\n";
    return true;
}

bool writeSurfaceGraphMatrixSnapshotJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_graph_matrix.v1\",\n";
    output << "  \"matrix_family\": \""
           << jsonEscape(graph_capacitance_matrix_provider
                             ? graph_capacitance_matrix_provider->matrixFamilyName()
                             : "BuiltinMatrixFamily")
           << "\",\n";
    output << "  \"solver_adapter_hint\": \""
           << jsonEscape(graph_capacitance_matrix_provider
                             ? graph_capacitance_matrix_provider->solverAdapterHint()
                             : "BuiltinFieldSolverAdapter")
           << "\",\n";
    output << "  \"graph_coupling_metric\": "
           << (graph_capacitance_matrix_provider && circuit_model
                   ? graph_capacitance_matrix_provider->graphCouplingMetric(config, circuit_model)
                   : 0.0)
           << ",\n";
    output << "  \"nodes\": [\n";
    if (circuit_model && graph_capacitance_matrix_provider)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            output << "    {\"index\": " << node_index << ", \"name\": \""
                   << jsonEscape(circuit_model->nodeName(node_index)) << "\", \"diagonal_f\": "
                   << graph_capacitance_matrix_provider->diagonalCapacitanceF(config, circuit_model, node_index)
                   << ", \"row_sum_f\": "
                   << graph_capacitance_matrix_provider->rowSumCapacitanceF(config, circuit_model, node_index)
                   << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"interfaces\": [\n";
    if (circuit_model && graph_capacitance_matrix_provider)
    {
        for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model->branchToNodeIndex(branch_index);
            output << "    {\"index\": " << branch_index << ", \"name\": \""
                   << jsonEscape(circuit_model->branchName(branch_index)) << "\", \"from_node\": \""
                   << jsonEscape(circuit_model->nodeName(from_node)) << "\", \"to_node\": \""
                   << jsonEscape(circuit_model->nodeName(to_node)) << "\", \"mutual_capacitance_f\": "
                   << graph_capacitance_matrix_provider->mutualCapacitanceF(config, circuit_model, branch_index)
                   << ", \"conductance_s\": " << circuit_model->branchConductanceS(branch_index)
                   << "}";
            output << (branch_index + 1 < circuit_model->branchCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeSurfaceFieldSolverAdapterContractReport(
    const std::filesystem::path& report_path, const SurfaceFieldSolverAdapter* field_solver_adapter,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider,
    const SurfaceCircuitModel* circuit_model)
{
    std::ofstream output(report_path);
    if (!output.is_open())
    {
        return false;
    }

    output << "surface_field_solver_adapter="
           << (field_solver_adapter ? field_solver_adapter->modelName() : "BuiltinFieldSolverAdapter")
           << "\n";
    output << "surface_graph_capacitance_matrix_provider="
           << (graph_capacitance_matrix_provider ? graph_capacitance_matrix_provider->modelName()
                                                 : "BuiltinGraphCapacitance")
           << "\n";
    output << "surface_graph_capacitance_matrix_family="
           << (graph_capacitance_matrix_provider ? graph_capacitance_matrix_provider->matrixFamilyName()
                                                 : "BuiltinMatrixFamily")
           << "\n";
    output << "surface_graph_capacitance_solver_adapter_hint="
           << (graph_capacitance_matrix_provider ? graph_capacitance_matrix_provider->solverAdapterHint()
                                                 : "BuiltinFieldSolverAdapter")
           << "\n";
    output << "expects_node_reference_inputs=1\n";
    output << "expects_node_charge_inputs=1\n";
    output << "expects_node_field_inputs=1\n";
    output << "expects_graph_capacitance_inputs=" << (graph_capacitance_matrix_provider ? 1 : 0) << "\n";
    output << "supports_patch_nodes=1\n";
    output << "supports_body_nodes=1\n";
    output << "supports_mutual_matrix="
           << (graph_capacitance_matrix_provider &&
                       graph_capacitance_matrix_provider->exposesMutualMatrix()
                   ? 1
                   : 0)
           << "\n";
    output << "supports_diagonal_matrix="
           << (graph_capacitance_matrix_provider &&
                       graph_capacitance_matrix_provider->exposesDiagonalMatrix()
                   ? 1
                   : 0)
           << "\n";
    output << "node_count=" << (circuit_model ? circuit_model->nodeCount() : 0) << "\n";
    output << "branch_count=" << (circuit_model ? circuit_model->branchCount() : 0) << "\n";
    return true;
}

bool writeSurfaceFieldSolverAdapterContractJson(
    const std::filesystem::path& json_path, const SurfaceFieldSolverAdapter* field_solver_adapter,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider,
    const SurfaceCircuitModel* circuit_model)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << "{\n";
    output << "  \"schema_version\": \"scdat.field_solver_adapter_contract.v1\",\n";
    output << "  \"surface_field_solver_adapter\": \""
           << jsonEscape(field_solver_adapter ? field_solver_adapter->modelName()
                                              : "BuiltinFieldSolverAdapter")
           << "\",\n";
    output << "  \"surface_graph_capacitance_matrix_provider\": \""
           << jsonEscape(graph_capacitance_matrix_provider
                             ? graph_capacitance_matrix_provider->modelName()
                             : "BuiltinGraphCapacitance")
           << "\",\n";
    output << "  \"surface_graph_capacitance_matrix_family\": \""
           << jsonEscape(graph_capacitance_matrix_provider
                             ? graph_capacitance_matrix_provider->matrixFamilyName()
                             : "BuiltinMatrixFamily")
           << "\",\n";
    output << "  \"surface_graph_capacitance_solver_adapter_hint\": \""
           << jsonEscape(graph_capacitance_matrix_provider
                             ? graph_capacitance_matrix_provider->solverAdapterHint()
                             : "BuiltinFieldSolverAdapter")
           << "\",\n";
    output << "  \"expects_node_reference_inputs\": true,\n";
    output << "  \"expects_node_charge_inputs\": true,\n";
    output << "  \"expects_node_field_inputs\": true,\n";
    output << "  \"expects_graph_capacitance_inputs\": "
           << (graph_capacitance_matrix_provider ? "true" : "false") << ",\n";
    output << "  \"supports_patch_nodes\": true,\n";
    output << "  \"supports_body_nodes\": true,\n";
    output << "  \"supports_mutual_matrix\": "
           << ((graph_capacitance_matrix_provider &&
                graph_capacitance_matrix_provider->exposesMutualMatrix())
                   ? "true"
                   : "false")
           << ",\n";
    output << "  \"supports_diagonal_matrix\": "
           << ((graph_capacitance_matrix_provider &&
                graph_capacitance_matrix_provider->exposesDiagonalMatrix())
                   ? "true"
                   : "false")
           << ",\n";
    output << "  \"node_count\": " << (circuit_model ? circuit_model->nodeCount() : 0) << ",\n";
    output << "  \"branch_count\": " << (circuit_model ? circuit_model->branchCount() : 0) << "\n";
    output << "}\n";
    return true;
}

bool writeSurfaceBoundaryMappingJson(const std::filesystem::path& json_path,
                                     const SurfaceChargingConfig& config,
                                     const SurfaceCircuitModel* circuit_model,
                                     const SurfaceCircuitMappingState& mapping_state)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_boundary_mapping.v1\",\n";
    output << "  \"body_boundary_groups\": [\n";
    for (std::size_t i = 0; i < config.body_boundary_groups.size(); ++i)
    {
        const auto& group = config.body_boundary_groups[i];
        output << "    {\"id\": \"" << jsonEscape(group.id) << "\", \"body_id\": \""
               << jsonEscape(group.body_id) << "\", \"external_group_name\": \""
               << jsonEscape(group.external_group_name) << "\", \"boundary_face_ids\": [";
        for (std::size_t j = 0; j < group.boundary_face_ids.size(); ++j)
        {
            if (j > 0)
            {
                output << ", ";
            }
            output << group.boundary_face_ids[j];
        }
        output << "]}";
        output << (i + 1 < config.body_boundary_groups.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"patch_boundary_groups\": [\n";
    for (std::size_t i = 0; i < config.patch_boundary_groups.size(); ++i)
    {
        const auto& group = config.patch_boundary_groups[i];
        output << "    {\"id\": \"" << jsonEscape(group.id) << "\", \"patch_id\": \""
               << jsonEscape(group.patch_id) << "\", \"external_group_name\": \""
               << jsonEscape(group.external_group_name) << "\", \"boundary_face_ids\": [";
        for (std::size_t j = 0; j < group.boundary_face_ids.size(); ++j)
        {
            if (j > 0)
            {
                output << ", ";
            }
            output << group.boundary_face_ids[j];
        }
        output << "]}";
        output << (i + 1 < config.patch_boundary_groups.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"boundary_mappings\": [\n";
    for (std::size_t i = 0; i < config.boundary_mappings.size(); ++i)
    {
        const auto& mapping = config.boundary_mappings[i];
        output << "    {\"node_id\": \"" << jsonEscape(mapping.node_id)
               << "\", \"boundary_group_id\": \"" << jsonEscape(mapping.boundary_group_id)
               << "\", \"required\": " << (mapping.required ? "true" : "false") << "}";
        output << (i + 1 < config.boundary_mappings.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"surface_to_circuit\": [\n";
    for (std::size_t surface_index = 0;
         surface_index < mapping_state.surface_to_circuit_node_index.size(); ++surface_index)
    {
        const auto circuit_index = mapping_state.surface_to_circuit_node_index[surface_index];
        const auto reduced_index =
            circuit_index < mapping_state.circuit_to_reduced_node_index.size()
                ? mapping_state.circuit_to_reduced_node_index[circuit_index]
                : std::numeric_limits<std::size_t>::max();
        const std::string node_name =
            circuit_model != nullptr && circuit_index < circuit_model->nodeCount()
                ? circuit_model->nodeName(circuit_index)
                : ("node:" + std::to_string(circuit_index));
        output << "    {\"surface_node_index\": " << surface_index
               << ", \"circuit_node_index\": " << circuit_index
               << ", \"node_name\": \"" << jsonEscape(node_name)
               << "\", \"boundary_group_id\": \"" << jsonEscape(boundaryGroupIdForNode(config, node_name))
               << "\", \"reduced_node_index\": ";
        if (reduced_index == std::numeric_limits<std::size_t>::max())
        {
            output << "null";
        }
        else
        {
            output << reduced_index;
        }
        output << "}";
        output << (surface_index + 1 < mapping_state.surface_to_circuit_node_index.size() ? ",\n"
                                                                                           : "\n");
    }
    output << "  ],\n";
    output << "  \"circuit_to_surface\": [\n";
    for (std::size_t circuit_index = 0;
         circuit_index < mapping_state.circuit_to_surface_node_indices.size(); ++circuit_index)
    {
        const auto reduced_index =
            circuit_index < mapping_state.circuit_to_reduced_node_index.size()
                ? mapping_state.circuit_to_reduced_node_index[circuit_index]
                : std::numeric_limits<std::size_t>::max();
        output << "    {\"circuit_node_index\": " << circuit_index
               << ", \"node_name\": \""
               << jsonEscape(circuit_model != nullptr && circuit_index < circuit_model->nodeCount()
                                 ? circuit_model->nodeName(circuit_index)
                                 : ("node:" + std::to_string(circuit_index)))
               << "\", \"surface_node_indices\": [";
        for (std::size_t j = 0; j < mapping_state.circuit_to_surface_node_indices[circuit_index].size();
             ++j)
        {
            if (j > 0)
            {
                output << ", ";
            }
            output << mapping_state.circuit_to_surface_node_indices[circuit_index][j];
        }
        output << "], \"reduced_node_index\": ";
        if (reduced_index == std::numeric_limits<std::size_t>::max())
        {
            output << "null";
        }
        else
        {
            output << reduced_index;
        }
        output << "}";
        output << (circuit_index + 1 < mapping_state.circuit_to_surface_node_indices.size() ? ",\n"
                                                                                            : "\n");
    }
    output << "  ],\n";
    output << "  \"reduced_node_groups\": [\n";
    for (std::size_t i = 0; i < mapping_state.reduced_node_groups.size(); ++i)
    {
        const auto& group = mapping_state.reduced_node_groups[i];
        output << "    {\"reduced_node_index\": " << group.reduced_node_index
               << ", \"group_id\": \"" << jsonEscape(group.group_id)
               << "\", \"member_node_indices\": [";
        for (std::size_t j = 0; j < group.member_node_indices.size(); ++j)
        {
            if (j > 0)
            {
                output << ", ";
            }
            output << group.member_node_indices[j];
        }
        output << "], \"total_area_m2\": " << group.total_area_m2 << "}";
        output << (i + 1 < mapping_state.reduced_node_groups.size() ? ",\n" : "\n");
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeSurfaceGraphReport(
    const std::filesystem::path& report_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model, const std::string& runtime_route_name,
    const std::string& circuit_model_name, const std::string& electric_field_provider_name,
    const std::string& reference_graph_propagator_name,
    const std::string& graph_capacitance_provider_name,
    const std::string& field_solver_adapter_name,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_total_current_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_charge_history,
    const std::vector<std::vector<double>>& node_propagated_reference_history,
    const std::vector<std::vector<double>>& node_field_solver_reference_history,
    const std::vector<std::vector<double>>& node_graph_capacitance_history,
    const std::vector<std::vector<double>>& node_graph_capacitance_row_sum_history,
    const std::vector<std::vector<double>>& node_field_solver_coupling_gain_history,
    const std::vector<std::vector<double>>& node_field_solver_capacitance_scale_history,
    const std::vector<std::vector<double>>& branch_current_history,
    const std::vector<std::vector<double>>& branch_conductance_history,
    const std::vector<std::vector<double>>& branch_voltage_drop_history,
    const std::vector<std::vector<double>>& branch_power_history,
    const std::vector<std::vector<double>>& branch_mutual_capacitance_history)
{
    std::ofstream output(report_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "runtime_route=" << runtime_route_name << "\n";
    output << "surface_circuit_model=" << circuit_model_name << "\n";
    output << "electric_field_provider=" << electric_field_provider_name << "\n";
    output << "surface_reference_graph_propagator=" << reference_graph_propagator_name << "\n";
    output << "surface_graph_capacitance_matrix_provider=" << graph_capacitance_provider_name
           << "\n";
    output << "surface_field_solver_adapter=" << field_solver_adapter_name << "\n";
    if (graph_capacitance_matrix_provider != nullptr)
    {
        output << "surface_graph_capacitance_matrix_family="
               << graph_capacitance_matrix_provider->matrixFamilyName() << "\n";
        output << "surface_graph_capacitance_solver_adapter_hint="
               << graph_capacitance_matrix_provider->solverAdapterHint() << "\n";
        output << "surface_graph_capacitance_exposes_mutual_matrix="
               << (graph_capacitance_matrix_provider->exposesMutualMatrix() ? 1 : 0) << "\n";
        output << "surface_graph_capacitance_exposes_diagonal_matrix="
               << (graph_capacitance_matrix_provider->exposesDiagonalMatrix() ? 1 : 0) << "\n";
    }

    if (circuit_model == nullptr)
    {
        output << "node_count=0\n";
        output << "branch_count=0\n";
        return true;
    }

    const std::size_t node_count = circuit_model->nodeCount();
    const std::size_t branch_count = circuit_model->branchCount();
    std::size_t body_count = 0;
    std::size_t patch_count = 0;
    std::size_t body_body_count = 0;
    std::size_t body_patch_count = 0;
    std::size_t patch_patch_count = 0;
    std::vector<std::size_t> node_degrees(node_count, 0);

    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        if (circuit_model->nodeIsPatch(node_index))
        {
            ++patch_count;
        }
        else
        {
            ++body_count;
        }
    }

    double total_graph_capacitance_diagonal_f = 0.0;
    double max_graph_capacitance_diagonal_f = 0.0;
    double total_graph_capacitance_row_sum_f = 0.0;
    double max_graph_capacitance_row_sum_f = 0.0;
    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        const double diagonal_f = latestHistoryValue(node_graph_capacitance_history, node_index);
        total_graph_capacitance_diagonal_f += diagonal_f;
        max_graph_capacitance_diagonal_f =
            std::max(max_graph_capacitance_diagonal_f, std::abs(diagonal_f));
        const double row_sum_f = latestHistoryValue(node_graph_capacitance_row_sum_history, node_index);
        total_graph_capacitance_row_sum_f += row_sum_f;
        max_graph_capacitance_row_sum_f =
            std::max(max_graph_capacitance_row_sum_f, std::abs(row_sum_f));
    }

    double total_mutual_capacitance_f = 0.0;
    double max_mutual_capacitance_f = 0.0;
    double max_abs_interface_current_a = 0.0;
    double max_abs_interface_voltage_drop_v = 0.0;
    double max_abs_interface_power_w = 0.0;
    double max_interface_conductance_s = 0.0;
    double max_abs_neighbor_potential_delta_v = 0.0;
    double max_abs_neighbor_field_contrast_v_per_m = 0.0;
    double max_abs_node_reference_offset_v = 0.0;
    double max_abs_node_field_solver_reference_offset_v = 0.0;
    double max_field_solver_coupling_gain = 0.0;
    std::size_t strongest_field_node_index = 0;
    double strongest_field_node_magnitude_v_per_m = 0.0;
    std::size_t strongest_coupling_branch_index = 0;
    double strongest_coupling_branch_metric = 0.0;
    for (std::size_t branch_index = 0; branch_index < branch_count; ++branch_index)
    {
        const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
        const auto to_node = circuit_model->branchToNodeIndex(branch_index);
        if (from_node < node_degrees.size())
        {
            ++node_degrees[from_node];
        }
        if (to_node < node_degrees.size())
        {
            ++node_degrees[to_node];
        }

        const auto branch_type = graphBranchTypeLabel(circuit_model, branch_index);
        if (branch_type == "body-body")
        {
            ++body_body_count;
        }
        else if (branch_type == "patch-patch")
        {
            ++patch_patch_count;
        }
        else if (branch_type == "body-patch")
        {
            ++body_patch_count;
        }

        const double mutual_capacitance_f = latestHistoryValue(branch_mutual_capacitance_history, branch_index);
        total_mutual_capacitance_f += mutual_capacitance_f;
        max_mutual_capacitance_f = std::max(max_mutual_capacitance_f, std::abs(mutual_capacitance_f));
        const double abs_current_a = std::abs(latestHistoryValue(branch_current_history, branch_index));
        max_abs_interface_current_a = std::max(max_abs_interface_current_a, abs_current_a);
        const double abs_voltage_drop_v =
            std::abs(latestHistoryValue(branch_voltage_drop_history, branch_index));
        max_abs_interface_voltage_drop_v = std::max(max_abs_interface_voltage_drop_v, abs_voltage_drop_v);
        const double abs_power_w = std::abs(latestHistoryValue(branch_power_history, branch_index));
        max_abs_interface_power_w = std::max(max_abs_interface_power_w, abs_power_w);
        max_interface_conductance_s =
            std::max(max_interface_conductance_s, latestHistoryValue(branch_conductance_history, branch_index));
        const double from_potential_v = latestHistoryValue(node_potential_history, from_node);
        const double to_potential_v = latestHistoryValue(node_potential_history, to_node);
        max_abs_neighbor_potential_delta_v =
            std::max(max_abs_neighbor_potential_delta_v, std::abs(from_potential_v - to_potential_v));
        const double from_field_v_per_m = latestHistoryValue(node_field_history, from_node);
        const double to_field_v_per_m = latestHistoryValue(node_field_history, to_node);
        max_abs_neighbor_field_contrast_v_per_m =
            std::max(max_abs_neighbor_field_contrast_v_per_m, std::abs(from_field_v_per_m - to_field_v_per_m));
        const double coupling_metric =
            std::abs(mutual_capacitance_f) *
                std::max(0.0, latestHistoryValue(branch_conductance_history, branch_index)) +
            abs_voltage_drop_v * std::max(0.0, latestHistoryValue(branch_conductance_history, branch_index));
        if (coupling_metric > strongest_coupling_branch_metric)
        {
            strongest_coupling_branch_metric = coupling_metric;
            strongest_coupling_branch_index = branch_index;
        }
    }

    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        const double field_magnitude_v_per_m = std::abs(latestHistoryValue(node_field_history, node_index));
        if (field_magnitude_v_per_m > strongest_field_node_magnitude_v_per_m)
        {
            strongest_field_node_magnitude_v_per_m = field_magnitude_v_per_m;
            strongest_field_node_index = node_index;
        }
        const double reference_offset_v =
            latestHistoryValue(node_propagated_reference_history, node_index) -
            latestHistoryValue(node_potential_history, node_index);
        max_abs_node_reference_offset_v =
            std::max(max_abs_node_reference_offset_v, std::abs(reference_offset_v));
        const double field_solver_reference_offset_v =
            latestHistoryValue(node_field_solver_reference_history, node_index) -
            latestHistoryValue(node_potential_history, node_index);
        max_abs_node_field_solver_reference_offset_v =
            std::max(max_abs_node_field_solver_reference_offset_v, std::abs(field_solver_reference_offset_v));
        max_field_solver_coupling_gain =
            std::max(max_field_solver_coupling_gain,
                     std::abs(latestHistoryValue(node_field_solver_coupling_gain_history, node_index)));
    }

    const std::size_t max_node_degree =
        node_degrees.empty() ? 0 : *std::max_element(node_degrees.begin(), node_degrees.end());
    double average_node_degree = 0.0;
    if (!node_degrees.empty())
    {
        const double degree_sum =
            std::accumulate(node_degrees.begin(), node_degrees.end(), 0.0,
                            [](double acc, std::size_t value) { return acc + static_cast<double>(value); });
        average_node_degree = degree_sum / static_cast<double>(node_degrees.size());
    }

    output << "node_count=" << node_count << "\n";
    output << "body_count=" << body_count << "\n";
    output << "patch_count=" << patch_count << "\n";
    output << "branch_count=" << branch_count << "\n";
    output << "body_body_interface_count=" << body_body_count << "\n";
    output << "body_patch_interface_count=" << body_patch_count << "\n";
    output << "patch_patch_interface_count=" << patch_patch_count << "\n";
    output << "max_node_degree=" << max_node_degree << "\n";
    output << "average_node_degree=" << average_node_degree << "\n";
    output << "total_graph_capacitance_diagonal_f=" << total_graph_capacitance_diagonal_f << "\n";
    output << "max_graph_capacitance_diagonal_f=" << max_graph_capacitance_diagonal_f << "\n";
    output << "total_graph_capacitance_row_sum_f=" << total_graph_capacitance_row_sum_f << "\n";
    output << "max_graph_capacitance_row_sum_f=" << max_graph_capacitance_row_sum_f << "\n";
    output << "total_interface_mutual_capacitance_f=" << total_mutual_capacitance_f << "\n";
    output << "max_interface_mutual_capacitance_f=" << max_mutual_capacitance_f << "\n";
    output << "max_abs_interface_current_a=" << max_abs_interface_current_a << "\n";
    output << "max_abs_interface_voltage_drop_v=" << max_abs_interface_voltage_drop_v << "\n";
    output << "max_abs_interface_power_w=" << max_abs_interface_power_w << "\n";
    output << "max_interface_conductance_s=" << max_interface_conductance_s << "\n";
    output << "max_abs_neighbor_potential_delta_v=" << max_abs_neighbor_potential_delta_v << "\n";
    output << "max_abs_neighbor_field_contrast_v_per_m=" << max_abs_neighbor_field_contrast_v_per_m << "\n";
    output << "max_abs_node_reference_offset_v=" << max_abs_node_reference_offset_v << "\n";
    output << "max_abs_node_field_solver_reference_offset_v=" << max_abs_node_field_solver_reference_offset_v << "\n";
    output << "max_field_solver_coupling_gain=" << max_field_solver_coupling_gain << "\n";
    output << "graph_coupling_metric="
           << (graph_capacitance_matrix_provider != nullptr
                   ? graph_capacitance_matrix_provider->graphCouplingMetric(config, circuit_model)
                   : 0.0)
           << "\n";
    if (node_count > 0)
    {
        output << "strongest_field_node_index=" << strongest_field_node_index << "\n";
        output << "strongest_field_node_name=" << circuit_model->nodeName(strongest_field_node_index) << "\n";
        output << "strongest_field_node_abs_field_v_per_m=" << strongest_field_node_magnitude_v_per_m << "\n";
    }
    if (branch_count > 0)
    {
        output << "strongest_coupling_interface_index=" << strongest_coupling_branch_index << "\n";
        output << "strongest_coupling_interface_name=" << circuit_model->branchName(strongest_coupling_branch_index) << "\n";
        output << "strongest_coupling_interface_metric=" << strongest_coupling_branch_metric << "\n";
    }

    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        const double reference_offset_v =
            latestHistoryValue(node_propagated_reference_history, node_index) -
            latestHistoryValue(node_potential_history, node_index);
        const double field_solver_reference_offset_v =
            latestHistoryValue(node_field_solver_reference_history, node_index) -
            latestHistoryValue(node_potential_history, node_index);
        const double local_charge_density_c_per_m3 = latestHistoryValue(node_charge_history, node_index);
        const double field_to_charge_ratio =
            std::abs(local_charge_density_c_per_m3) <= 1.0e-30
                ? 0.0
                : latestHistoryValue(node_field_history, node_index) / local_charge_density_c_per_m3;
        output << "node_" << node_index << "="
               << "name:" << circuit_model->nodeName(node_index)
               << ",type:" << graphNodeTypeLabel(circuit_model->nodeIsPatch(node_index))
               << ",owner:" << graphNodeOwnerLabel(config, circuit_model->nodeName(node_index),
                                                   circuit_model->nodeIsPatch(node_index))
               << ",material:" << graphNodeMaterialName(config, node_index, circuit_model->nodeName(node_index),
                                                        circuit_model->nodeIsPatch(node_index))
               << ",degree:" << node_degrees[node_index]
               << ",area_m2:" << circuit_model->nodeAreaM2(node_index)
               << ",last_potential_v:" << latestHistoryValue(node_potential_history, node_index)
               << ",last_total_current_density_a_per_m2:" << latestHistoryValue(node_total_current_history, node_index)
               << ",last_normal_electric_field_v_per_m:" << latestHistoryValue(node_field_history, node_index)
               << ",last_local_charge_density_c_per_m3:" << local_charge_density_c_per_m3
               << ",last_propagated_reference_potential_v:" << latestHistoryValue(node_propagated_reference_history, node_index)
               << ",last_reference_offset_v:" << reference_offset_v
               << ",last_field_solver_reference_potential_v:" << latestHistoryValue(node_field_solver_reference_history, node_index)
               << ",last_field_solver_reference_offset_v:" << field_solver_reference_offset_v
               << ",last_graph_capacitance_diagonal_f:" << latestHistoryValue(node_graph_capacitance_history, node_index)
               << ",last_graph_capacitance_row_sum_f:" << latestHistoryValue(node_graph_capacitance_row_sum_history, node_index)
               << ",last_field_solver_coupling_gain:" << latestHistoryValue(node_field_solver_coupling_gain_history, node_index)
               << ",last_field_solver_capacitance_scale:" << latestHistoryValue(node_field_solver_capacitance_scale_history, node_index)
               << ",field_to_charge_ratio:" << field_to_charge_ratio << "\n";
    }

    for (std::size_t branch_index = 0; branch_index < branch_count; ++branch_index)
    {
        const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
        const auto to_node = circuit_model->branchToNodeIndex(branch_index);
        output << "interface_" << branch_index << "="
               << "name:" << circuit_model->branchName(branch_index)
               << ",type:" << graphBranchTypeLabel(circuit_model, branch_index)
               << ",from_node:" << circuit_model->nodeName(from_node)
               << ",to_node:" << circuit_model->nodeName(to_node)
               << ",conductance_s:" << circuit_model->branchConductanceS(branch_index)
               << ",resistance_ohm:" << circuit_model->branchResistanceOhm(branch_index)
               << ",bias_v:" << circuit_model->branchBiasV(branch_index)
               << ",dynamic_conductance:" << (circuit_model->branchUsesDynamicConductance(branch_index) ? 1 : 0)
               << ",last_current_a:" << latestHistoryValue(branch_current_history, branch_index)
               << ",last_voltage_drop_v:" << latestHistoryValue(branch_voltage_drop_history, branch_index)
               << ",last_power_w:" << latestHistoryValue(branch_power_history, branch_index)
               << ",last_mutual_capacitance_f:" << latestHistoryValue(branch_mutual_capacitance_history, branch_index)
               << "\n";
    }
    return true;
}

bool writeSharedSurfaceRuntimeObserverJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model, const std::vector<double>& time_history,
    const std::vector<double>& shared_patch_potential_history,
    const std::vector<double>& shared_patch_area_history,
    const std::vector<double>& shared_reference_potential_history,
    const std::vector<double>& shared_effective_sheath_length_history,
    const std::vector<double>& shared_sheath_charge_history,
    const std::vector<double>& shared_runtime_enabled_history,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_shared_runtime_enabled_history,
    const std::vector<std::vector<double>>& node_shared_patch_potential_history,
    const std::vector<std::vector<double>>& node_shared_reference_potential_history,
    const std::vector<std::vector<double>>& node_shared_sheath_charge_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_shared_runtime_observer.v1\",\n";
    output << "  \"contract_id\": \"surface-pic-boundary-observer-v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_runtime\": \"" << surfacePicRuntimeName(config.surface_pic_runtime_kind)
           << "\",\n";
    output << "  \"surface_pic_strategy\": \"" << surfacePicStrategyName(config.surface_pic_strategy)
           << "\",\n";
    output << "  \"shared_runtime_enabled\": "
           << (latest_scalar(shared_runtime_enabled_history) >= 0.5 ? "true" : "false") << ",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << latest_scalar(time_history) << ",\n";
    output << "  \"shared_surface\": {\n";
    output << "    \"latest_patch_potential_v\": " << latest_scalar(shared_patch_potential_history)
           << ",\n";
    output << "    \"latest_patch_area_m2\": " << latest_scalar(shared_patch_area_history) << ",\n";
    output << "    \"latest_reference_potential_v\": " << latest_scalar(shared_reference_potential_history)
           << ",\n";
    output << "    \"latest_effective_sheath_length_m\": "
           << latest_scalar(shared_effective_sheath_length_history) << ",\n";
    output << "    \"latest_sheath_charge_c\": " << latest_scalar(shared_sheath_charge_history)
           << "\n";
    output << "  },\n";
    output << "  \"nodes\": [\n";
    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            output << "    {\"index\": " << node_index << ", \"name\": \""
                   << jsonEscape(circuit_model->nodeName(node_index)) << "\", \"is_patch\": "
                   << (circuit_model->nodeIsPatch(node_index) ? "true" : "false")
                   << ", \"latest_potential_v\": "
                   << latestHistoryValue(node_potential_history, node_index)
                   << ", \"shared_runtime_enabled\": "
                   << (latestHistoryValue(node_shared_runtime_enabled_history, node_index) >= 0.5 ? "true"
                                                                                                   : "false")
                   << ", \"latest_shared_patch_potential_v\": "
                   << latestHistoryValue(node_shared_patch_potential_history, node_index)
                   << ", \"latest_shared_reference_potential_v\": "
                   << latestHistoryValue(node_shared_reference_potential_history, node_index)
                   << ", \"latest_shared_sheath_charge_c\": "
                   << latestHistoryValue(node_shared_sheath_charge_history, node_index) << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeSurfaceMonitorJson(const std::filesystem::path& json_path,
                             const SurfaceCircuitModel* circuit_model,
                             const SurfaceGraphMonitorSnapshot& graph_snapshot,
                             const SurfaceFieldMonitorSnapshot& field_snapshot,
                             const SurfaceBenchmarkMonitorSnapshot& benchmark_snapshot)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_monitor.v1\",\n";
    output << "  \"graph\": {\n";
    output << "    \"node_count\": " << graph_snapshot.node_count << ",\n";
    output << "    \"branch_count\": " << graph_snapshot.branch_count << ",\n";
    output << "    \"max_node_degree\": " << graph_snapshot.max_node_degree << ",\n";
    output << "    \"average_node_degree\": " << graph_snapshot.average_node_degree << ",\n";
    output << "    \"graph_coupling_metric\": " << graph_snapshot.graph_coupling_metric << ",\n";
    output << "    \"strongest_field_node_index\": " << graph_snapshot.strongest_field_node_index
           << ",\n";
    output << "    \"strongest_field_node_name\": \""
           << jsonEscape((circuit_model && graph_snapshot.node_count > 0)
                             ? circuit_model->nodeName(graph_snapshot.strongest_field_node_index)
                             : "")
           << "\",\n";
    output << "    \"strongest_field_node_abs_v_per_m\": "
           << graph_snapshot.strongest_field_node_abs_v_per_m << ",\n";
    output << "    \"strongest_coupling_interface_index\": "
           << graph_snapshot.strongest_coupling_interface_index << ",\n";
    output << "    \"strongest_coupling_interface_name\": \""
           << jsonEscape((circuit_model && graph_snapshot.branch_count > 0)
                             ? circuit_model->branchName(
                                   graph_snapshot.strongest_coupling_interface_index)
                             : "")
           << "\",\n";
    output << "    \"strongest_coupling_interface_metric\": "
           << graph_snapshot.strongest_coupling_interface_metric << ",\n";
    output << "    \"max_abs_neighbor_potential_delta_v\": "
           << graph_snapshot.max_abs_neighbor_potential_delta_v << ",\n";
    output << "    \"max_abs_neighbor_field_contrast_v_per_m\": "
           << graph_snapshot.max_abs_neighbor_field_contrast_v_per_m << "\n";
    output << "  },\n";
    output << "  \"field\": {\n";
    output << "    \"strongest_field_node_index\": " << field_snapshot.strongest_field_node_index
           << ",\n";
    output << "    \"strongest_field_node_name\": \""
           << jsonEscape((circuit_model && graph_snapshot.node_count > 0)
                             ? circuit_model->nodeName(field_snapshot.strongest_field_node_index)
                             : "")
           << "\",\n";
    output << "    \"strongest_field_node_abs_v_per_m\": "
           << field_snapshot.strongest_field_node_abs_v_per_m << ",\n";
    output << "    \"max_abs_node_reference_offset_v\": "
           << field_snapshot.max_abs_node_reference_offset_v << ",\n";
    output << "    \"max_abs_node_field_solver_reference_offset_v\": "
           << field_snapshot.max_abs_node_field_solver_reference_offset_v << ",\n";
    output << "    \"max_field_solver_coupling_gain\": "
           << field_snapshot.max_field_solver_coupling_gain << ",\n";
    output << "    \"reference_offset_envelope_v\": "
           << field_snapshot.reference_offset_envelope_v << "\n";
    output << "  },\n";
    output << "  \"benchmark\": {\n";
    output << "    \"runtime_route\": \"" << jsonEscape(benchmark_snapshot.runtime_route)
           << "\",\n";
    output << "    \"benchmark_mode\": \"" << jsonEscape(benchmark_snapshot.benchmark_mode)
           << "\",\n";
    output << "    \"benchmark_source\": \"" << jsonEscape(benchmark_snapshot.benchmark_source)
           << "\",\n";
    output << "    \"benchmark_execution_mode\": \""
           << jsonEscape(benchmark_snapshot.benchmark_execution_mode) << "\",\n";
    output << "    \"consistency_status\": \""
           << jsonEscape(benchmark_snapshot.consistency_status) << "\",\n";
    output << "    \"consistency_authority\": \""
           << jsonEscape(benchmark_snapshot.consistency_authority) << "\",\n";
    output << "    \"patch_rmse_v\": " << benchmark_snapshot.patch_rmse_v << ",\n";
    output << "    \"body_rmse_v\": " << benchmark_snapshot.body_rmse_v << ",\n";
    output << "    \"compared_patch_samples\": " << benchmark_snapshot.compared_patch_samples
           << ",\n";
    output << "    \"compared_body_samples\": " << benchmark_snapshot.compared_body_samples
           << ",\n";
    output << "    \"terminal_potential_v\": " << benchmark_snapshot.terminal_potential_v << ",\n";
    output << "    \"time_to_equilibrium_ms\": " << benchmark_snapshot.time_to_equilibrium_ms
           << "\n";
    output << "  }\n";
    output << "}\n";
    return true;
}

bool writeSharedSurfaceRuntimeConsistencyJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model, const std::vector<double>& time_history,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<double>& pre_global_solve_patch_spread_history,
    const std::vector<double>& patch_spread_reduction_v_history,
    const std::vector<double>& patch_spread_reduction_ratio_history,
    const std::vector<double>& shared_current_matrix_coupling_active_history,
    const std::vector<double>& shared_current_matrix_coupling_offdiag_entry_history,
    const std::vector<double>& shared_global_coupled_solve_active_history,
    const std::vector<double>& shared_global_coupled_solve_iteration_history,
    const std::vector<double>& shared_global_coupled_solve_converged_history,
    const std::vector<double>& shared_global_coupled_solve_max_delta_history,
    const std::vector<double>& shared_live_pic_coupled_refresh_active_history,
    const std::vector<double>& shared_live_pic_coupled_refresh_count_history,
    const std::vector<double>& shared_particle_transport_coupling_active_history,
    const std::vector<double>& shared_particle_transport_offdiag_entry_history,
    const std::vector<double>& shared_particle_transport_total_conductance_history,
    const std::vector<double>& shared_particle_transport_conservation_error_history,
    const std::vector<double>& shared_particle_transport_charge_history,
    const std::vector<double>& shared_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_normal_electric_field_history,
    const std::vector<std::vector<double>>& node_local_charge_density_history,
    const std::vector<std::vector<double>>& node_shared_runtime_enabled_history,
    const std::vector<std::vector<double>>& node_shared_reference_potential_history,
    const std::vector<std::vector<double>>& node_shared_sheath_charge_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_charge_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_net_flux_history,
    const std::vector<std::vector<double>>& edge_particle_transport_charge_matrix_c,
    const std::vector<std::vector<double>>& edge_particle_transport_target_charge_matrix_c,
    const std::vector<std::vector<double>>& edge_particle_transport_operator_drive_matrix_c,
    double shared_particle_transport_edge_graph_operator_iterations,
    bool shared_particle_transport_edge_graph_operator_converged,
    double shared_particle_transport_edge_graph_operator_max_balance_residual_c,
    double shared_particle_transport_edge_graph_operator_branch_graph_edge_count,
    double shared_particle_transport_edge_graph_operator_branch_graph_pair_count,
    double shared_particle_transport_edge_graph_operator_effective_pair_count,
    double shared_particle_transport_edge_graph_operator_total_pair_weight_f,
    double shared_particle_transport_edge_graph_operator_total_conductance_weight_f,
    double shared_particle_transport_edge_graph_operator_min_node_preconditioner,
    double shared_particle_transport_edge_graph_operator_max_node_preconditioner,
    const std::vector<std::vector<double>>& node_live_pic_electron_current_history,
    const std::vector<std::vector<double>>& node_live_pic_ion_current_history,
    const std::vector<std::vector<double>>& node_live_pic_net_current_history,
    const std::vector<std::vector<double>>& node_live_pic_collision_count_history,
    const std::vector<std::vector<double>>& node_electron_calibration_factor_history,
    const std::vector<std::vector<double>>& node_ion_calibration_factor_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    std::vector<std::size_t> shared_patch_nodes;
    double min_reference_potential_v = 0.0;
    double max_reference_potential_v = 0.0;
    double min_sheath_charge_c = 0.0;
    double max_sheath_charge_c = 0.0;
    double min_patch_potential_v = 0.0;
    double max_patch_potential_v = 0.0;
    double min_live_pic_electron_current_a_per_m2 = 0.0;
    double max_live_pic_electron_current_a_per_m2 = 0.0;
    double min_live_pic_ion_current_a_per_m2 = 0.0;
    double max_live_pic_ion_current_a_per_m2 = 0.0;
    double min_live_pic_net_current_a_per_m2 = 0.0;
    double max_live_pic_net_current_a_per_m2 = 0.0;
    double min_live_pic_collision_count = 0.0;
    double max_live_pic_collision_count = 0.0;
    double min_electron_calibration_factor = 0.0;
    double max_electron_calibration_factor = 0.0;
    double min_ion_calibration_factor = 0.0;
    double max_ion_calibration_factor = 0.0;
    double min_normal_electric_field_v_per_m = 0.0;
    double max_normal_electric_field_v_per_m = 0.0;
    double min_local_charge_density_c_per_m3 = 0.0;
    double max_local_charge_density_c_per_m3 = 0.0;
    double min_distributed_particle_transport_charge_c = 0.0;
    double max_distributed_particle_transport_charge_c = 0.0;
    double min_distributed_particle_transport_reference_shift_v = 0.0;
    double max_distributed_particle_transport_reference_shift_v = 0.0;
    double sum_distributed_particle_transport_charge_c = 0.0;
    double min_distributed_particle_transport_net_flux_a = 0.0;
    double max_distributed_particle_transport_net_flux_a = 0.0;
    double sum_distributed_particle_transport_net_flux_a = 0.0;
    double min_node_edge_transport_charge_c = 0.0;
    double max_node_edge_transport_charge_c = 0.0;
    double sum_node_edge_transport_charge_c = 0.0;
    double total_abs_edge_transport_charge_c = 0.0;
    double total_abs_edge_target_charge_c = 0.0;
    double total_abs_edge_operator_drive_charge_c = 0.0;
    double total_node_edge_operator_drive_charge_c = 0.0;
    bool first_shared_node = true;

    if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            if (!circuit_model->nodeIsPatch(node_index))
            {
                continue;
            }
            if (latestHistoryValue(node_shared_runtime_enabled_history, node_index) < 0.5)
            {
                continue;
            }

            const double shared_reference = latestHistoryValue(node_shared_reference_potential_history, node_index);
            const double shared_sheath_charge = latestHistoryValue(node_shared_sheath_charge_history, node_index);
            const double patch_potential = latestHistoryValue(node_potential_history, node_index);
            const double live_pic_electron_current_a_per_m2 =
                latestHistoryValue(node_live_pic_electron_current_history, node_index);
            const double live_pic_ion_current_a_per_m2 =
                latestHistoryValue(node_live_pic_ion_current_history, node_index);
            const double live_pic_net_current_a_per_m2 =
                latestHistoryValue(node_live_pic_net_current_history, node_index);
            const double live_pic_collision_count =
                latestHistoryValue(node_live_pic_collision_count_history, node_index);
            const double electron_calibration_factor =
                latestHistoryValue(node_electron_calibration_factor_history, node_index);
            const double ion_calibration_factor =
                latestHistoryValue(node_ion_calibration_factor_history, node_index);
            const double normal_electric_field_v_per_m =
                latestHistoryValue(node_normal_electric_field_history, node_index);
            const double local_charge_density_c_per_m3 =
                latestHistoryValue(node_local_charge_density_history, node_index);
            const double distributed_particle_transport_charge_c =
                latestHistoryValue(node_distributed_particle_transport_charge_history, node_index);
            const double distributed_particle_transport_reference_shift_v =
                latestHistoryValue(node_distributed_particle_transport_reference_shift_history, node_index);
            const double distributed_particle_transport_net_flux_a =
                latestHistoryValue(node_distributed_particle_transport_net_flux_history, node_index);
            double node_edge_transport_charge_c = 0.0;
            if (node_index < edge_particle_transport_charge_matrix_c.size())
            {
                for (const double edge_charge_c : edge_particle_transport_charge_matrix_c[node_index])
                {
                    node_edge_transport_charge_c += edge_charge_c;
                }
            }
            double node_edge_operator_drive_charge_c = 0.0;
            if (node_index < edge_particle_transport_operator_drive_matrix_c.size())
            {
                for (const double operator_drive_charge_c : edge_particle_transport_operator_drive_matrix_c[node_index])
                {
                    node_edge_operator_drive_charge_c += operator_drive_charge_c;
                }
            }

            shared_patch_nodes.push_back(node_index);
            if (first_shared_node)
            {
                min_reference_potential_v = max_reference_potential_v = shared_reference;
                min_sheath_charge_c = max_sheath_charge_c = shared_sheath_charge;
                min_patch_potential_v = max_patch_potential_v = patch_potential;
                min_live_pic_electron_current_a_per_m2 = max_live_pic_electron_current_a_per_m2 =
                    live_pic_electron_current_a_per_m2;
                min_live_pic_ion_current_a_per_m2 = max_live_pic_ion_current_a_per_m2 =
                    live_pic_ion_current_a_per_m2;
                min_live_pic_net_current_a_per_m2 = max_live_pic_net_current_a_per_m2 =
                    live_pic_net_current_a_per_m2;
                min_live_pic_collision_count = max_live_pic_collision_count = live_pic_collision_count;
                min_electron_calibration_factor = max_electron_calibration_factor = electron_calibration_factor;
                min_ion_calibration_factor = max_ion_calibration_factor = ion_calibration_factor;
                min_normal_electric_field_v_per_m = max_normal_electric_field_v_per_m = normal_electric_field_v_per_m;
                min_local_charge_density_c_per_m3 = max_local_charge_density_c_per_m3 = local_charge_density_c_per_m3;
                min_distributed_particle_transport_charge_c = max_distributed_particle_transport_charge_c =
                    distributed_particle_transport_charge_c;
                min_distributed_particle_transport_reference_shift_v = max_distributed_particle_transport_reference_shift_v =
                    distributed_particle_transport_reference_shift_v;
                sum_distributed_particle_transport_charge_c = distributed_particle_transport_charge_c;
                min_distributed_particle_transport_net_flux_a = max_distributed_particle_transport_net_flux_a =
                    distributed_particle_transport_net_flux_a;
                sum_distributed_particle_transport_net_flux_a = distributed_particle_transport_net_flux_a;
                min_node_edge_transport_charge_c = max_node_edge_transport_charge_c = node_edge_transport_charge_c;
                sum_node_edge_transport_charge_c = node_edge_transport_charge_c;
                total_node_edge_operator_drive_charge_c = node_edge_operator_drive_charge_c;
                first_shared_node = false;
            }
            else
            {
                min_reference_potential_v = std::min(min_reference_potential_v, shared_reference);
                max_reference_potential_v = std::max(max_reference_potential_v, shared_reference);
                min_sheath_charge_c = std::min(min_sheath_charge_c, shared_sheath_charge);
                max_sheath_charge_c = std::max(max_sheath_charge_c, shared_sheath_charge);
                min_patch_potential_v = std::min(min_patch_potential_v, patch_potential);
                max_patch_potential_v = std::max(max_patch_potential_v, patch_potential);
                min_live_pic_electron_current_a_per_m2 =
                    std::min(min_live_pic_electron_current_a_per_m2, live_pic_electron_current_a_per_m2);
                max_live_pic_electron_current_a_per_m2 =
                    std::max(max_live_pic_electron_current_a_per_m2, live_pic_electron_current_a_per_m2);
                min_live_pic_ion_current_a_per_m2 =
                    std::min(min_live_pic_ion_current_a_per_m2, live_pic_ion_current_a_per_m2);
                max_live_pic_ion_current_a_per_m2 =
                    std::max(max_live_pic_ion_current_a_per_m2, live_pic_ion_current_a_per_m2);
                min_live_pic_net_current_a_per_m2 =
                    std::min(min_live_pic_net_current_a_per_m2, live_pic_net_current_a_per_m2);
                max_live_pic_net_current_a_per_m2 =
                    std::max(max_live_pic_net_current_a_per_m2, live_pic_net_current_a_per_m2);
                min_live_pic_collision_count = std::min(min_live_pic_collision_count, live_pic_collision_count);
                max_live_pic_collision_count = std::max(max_live_pic_collision_count, live_pic_collision_count);
                min_electron_calibration_factor = std::min(min_electron_calibration_factor, electron_calibration_factor);
                max_electron_calibration_factor = std::max(max_electron_calibration_factor, electron_calibration_factor);
                min_ion_calibration_factor = std::min(min_ion_calibration_factor, ion_calibration_factor);
                max_ion_calibration_factor = std::max(max_ion_calibration_factor, ion_calibration_factor);
                min_normal_electric_field_v_per_m = std::min(min_normal_electric_field_v_per_m, normal_electric_field_v_per_m);
                max_normal_electric_field_v_per_m = std::max(max_normal_electric_field_v_per_m, normal_electric_field_v_per_m);
                min_local_charge_density_c_per_m3 = std::min(min_local_charge_density_c_per_m3, local_charge_density_c_per_m3);
                max_local_charge_density_c_per_m3 = std::max(max_local_charge_density_c_per_m3, local_charge_density_c_per_m3);
                min_distributed_particle_transport_charge_c =
                    std::min(min_distributed_particle_transport_charge_c, distributed_particle_transport_charge_c);
                max_distributed_particle_transport_charge_c =
                    std::max(max_distributed_particle_transport_charge_c, distributed_particle_transport_charge_c);
                min_distributed_particle_transport_reference_shift_v =
                    std::min(min_distributed_particle_transport_reference_shift_v, distributed_particle_transport_reference_shift_v);
                max_distributed_particle_transport_reference_shift_v =
                    std::max(max_distributed_particle_transport_reference_shift_v, distributed_particle_transport_reference_shift_v);
                sum_distributed_particle_transport_charge_c += distributed_particle_transport_charge_c;
                min_distributed_particle_transport_net_flux_a =
                    std::min(min_distributed_particle_transport_net_flux_a, distributed_particle_transport_net_flux_a);
                max_distributed_particle_transport_net_flux_a =
                    std::max(max_distributed_particle_transport_net_flux_a, distributed_particle_transport_net_flux_a);
                sum_distributed_particle_transport_net_flux_a += distributed_particle_transport_net_flux_a;
                min_node_edge_transport_charge_c = std::min(min_node_edge_transport_charge_c, node_edge_transport_charge_c);
                max_node_edge_transport_charge_c = std::max(max_node_edge_transport_charge_c, node_edge_transport_charge_c);
                sum_node_edge_transport_charge_c += node_edge_transport_charge_c;
                total_node_edge_operator_drive_charge_c += node_edge_operator_drive_charge_c;
            }
        }

        for (std::size_t i = 0; i < shared_patch_nodes.size(); ++i)
        {
            const auto node_i = shared_patch_nodes[i];
            for (std::size_t j = i + 1; j < shared_patch_nodes.size(); ++j)
            {
                const auto node_j = shared_patch_nodes[j];
                if (node_i < edge_particle_transport_charge_matrix_c.size() &&
                    node_j < edge_particle_transport_charge_matrix_c[node_i].size())
                {
                    total_abs_edge_transport_charge_c += std::abs(edge_particle_transport_charge_matrix_c[node_i][node_j]);
                }
                if (node_i < edge_particle_transport_target_charge_matrix_c.size() &&
                    node_j < edge_particle_transport_target_charge_matrix_c[node_i].size())
                {
                    total_abs_edge_target_charge_c += std::abs(edge_particle_transport_target_charge_matrix_c[node_i][node_j]);
                }
                if (node_i < edge_particle_transport_operator_drive_matrix_c.size() &&
                    node_j < edge_particle_transport_operator_drive_matrix_c[node_i].size())
                {
                    total_abs_edge_operator_drive_charge_c += std::abs(edge_particle_transport_operator_drive_matrix_c[node_i][node_j]);
                }
            }
        }
    }

    const double reference_spread_v = first_shared_node ? 0.0 : (max_reference_potential_v - min_reference_potential_v);
    const double sheath_charge_spread_c = first_shared_node ? 0.0 : (max_sheath_charge_c - min_sheath_charge_c);
    const double patch_potential_spread_v = first_shared_node ? 0.0 : (max_patch_potential_v - min_patch_potential_v);
    const double live_pic_electron_current_spread_a_per_m2 =
        first_shared_node ? 0.0 : (max_live_pic_electron_current_a_per_m2 - min_live_pic_electron_current_a_per_m2);
    const double live_pic_ion_current_spread_a_per_m2 =
        first_shared_node ? 0.0 : (max_live_pic_ion_current_a_per_m2 - min_live_pic_ion_current_a_per_m2);
    const double live_pic_net_current_spread_a_per_m2 =
        first_shared_node ? 0.0 : (max_live_pic_net_current_a_per_m2 - min_live_pic_net_current_a_per_m2);
    const double live_pic_collision_count_spread =
        first_shared_node ? 0.0 : (max_live_pic_collision_count - min_live_pic_collision_count);
    const double electron_calibration_factor_spread =
        first_shared_node ? 0.0 : (max_electron_calibration_factor - min_electron_calibration_factor);
    const double ion_calibration_factor_spread =
        first_shared_node ? 0.0 : (max_ion_calibration_factor - min_ion_calibration_factor);
    const double normal_electric_field_spread_v_per_m =
        first_shared_node ? 0.0 : (max_normal_electric_field_v_per_m - min_normal_electric_field_v_per_m);
    const double local_charge_density_spread_c_per_m3 =
        first_shared_node ? 0.0 : (max_local_charge_density_c_per_m3 - min_local_charge_density_c_per_m3);
    const double distributed_particle_transport_charge_spread_c =
        first_shared_node ? 0.0 : (max_distributed_particle_transport_charge_c - min_distributed_particle_transport_charge_c);
    const double distributed_particle_transport_reference_shift_spread_v =
        first_shared_node ? 0.0 : (max_distributed_particle_transport_reference_shift_v - min_distributed_particle_transport_reference_shift_v);
    const double distributed_particle_transport_net_flux_spread_a =
        first_shared_node ? 0.0 : (max_distributed_particle_transport_net_flux_a - min_distributed_particle_transport_net_flux_a);
    const double node_edge_transport_charge_spread_c =
        first_shared_node ? 0.0 : (max_node_edge_transport_charge_c - min_node_edge_transport_charge_c);
    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };
    const double pre_global_solve_patch_potential_spread_v = latest_scalar(pre_global_solve_patch_spread_history);
    const double patch_potential_spread_reduction_v = latest_scalar(patch_spread_reduction_v_history);
    const double patch_potential_spread_reduction_ratio = latest_scalar(patch_spread_reduction_ratio_history);
    const bool shared_current_matrix_coupling_active = latest_scalar(shared_current_matrix_coupling_active_history) >= 0.5;
    const int shared_current_matrix_coupling_offdiag_entry_count =
        static_cast<int>(std::llround(latest_scalar(shared_current_matrix_coupling_offdiag_entry_history)));
    const bool shared_global_coupled_solve_active = latest_scalar(shared_global_coupled_solve_active_history) >= 0.5;
    const int shared_global_coupled_solve_iterations =
        static_cast<int>(std::llround(latest_scalar(shared_global_coupled_solve_iteration_history)));
    const bool shared_global_coupled_solve_converged = latest_scalar(shared_global_coupled_solve_converged_history) >= 0.5;
    const double shared_global_coupled_solve_max_delta_v = latest_scalar(shared_global_coupled_solve_max_delta_history);
    const bool shared_live_pic_coupled_refresh_active = latest_scalar(shared_live_pic_coupled_refresh_active_history) >= 0.5;
    const int shared_live_pic_coupled_refresh_count =
        static_cast<int>(std::llround(latest_scalar(shared_live_pic_coupled_refresh_count_history)));
    const bool shared_particle_transport_coupling_active =
        latest_scalar(shared_particle_transport_coupling_active_history) >= 0.5;
    const int shared_particle_transport_offdiag_entry_count =
        static_cast<int>(std::llround(latest_scalar(shared_particle_transport_offdiag_entry_history)));
    const double shared_particle_transport_total_conductance_s =
        latest_scalar(shared_particle_transport_total_conductance_history);
    const double shared_particle_transport_conservation_error_a_per_v =
        latest_scalar(shared_particle_transport_conservation_error_history);
    const double shared_particle_transport_charge_c = latest_scalar(shared_particle_transport_charge_history);
    const double shared_particle_transport_reference_shift_v =
        latest_scalar(shared_particle_transport_reference_shift_history);
    const bool shared_particle_transport_distribution_active =
        shared_patch_nodes.size() >= 2 && distributed_particle_transport_charge_spread_c > 0.0;
    const double shared_particle_transport_distribution_conservation_error_c =
        std::abs(sum_distributed_particle_transport_charge_c - shared_particle_transport_charge_c);
    const double shared_particle_transport_exchange_flux_conservation_error_a =
        std::abs(sum_distributed_particle_transport_net_flux_a);
    const bool shared_particle_transport_exchange_active = distributed_particle_transport_net_flux_spread_a > 0.0;
    const bool shared_particle_transport_edge_domain_active = total_abs_edge_transport_charge_c > 0.0;
    const double shared_particle_transport_edge_charge_conservation_error_c = std::abs(sum_node_edge_transport_charge_c);
    const bool shared_particle_transport_edge_operator_active = total_abs_edge_operator_drive_charge_c > 0.0;
    const double shared_particle_transport_edge_operator_drive_conservation_error_c =
        std::abs(total_node_edge_operator_drive_charge_c);
    const double shared_global_solve_weight = std::clamp(
        config.material.getScalarProperty("shared_surface_global_solve_weight", 0.0), 0.0, 1.0);
    const bool global_sheath_proxy_solve_active = shared_global_solve_weight > 0.0 && shared_patch_nodes.size() >= 2;
    const bool sheath_consistency_pass = reference_spread_v <= 1.0e-9 && sheath_charge_spread_c <= 1.0e-15;
    const bool shared_particle_pool_consistency_pass =
        live_pic_electron_current_spread_a_per_m2 <= 1.0e-12 &&
        live_pic_ion_current_spread_a_per_m2 <= 1.0e-12 &&
        live_pic_net_current_spread_a_per_m2 <= 1.0e-12 &&
        live_pic_collision_count_spread <= 0.5 &&
        electron_calibration_factor_spread <= 1.0e-12 &&
        ion_calibration_factor_spread <= 1.0e-12;
    const bool consistency_pass = sheath_consistency_pass && shared_particle_pool_consistency_pass;

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_shared_runtime_consistency.v1\",\n";
    output << "  \"contract_id\": \"surface-pic-shared-runtime-consistency-v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_runtime\": \"" << surfacePicRuntimeName(config.surface_pic_runtime_kind)
           << "\",\n";
    output << "  \"surface_pic_strategy\": \"" << surfacePicStrategyName(config.surface_pic_strategy)
           << "\",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << (time_history.empty() ? 0.0 : time_history.back()) << ",\n";
    output << "  \"shared_patch_node_count\": " << shared_patch_nodes.size() << ",\n";
    output << "  \"reference_potential_spread_v\": " << reference_spread_v << ",\n";
    output << "  \"sheath_charge_spread_c\": " << sheath_charge_spread_c << ",\n";
    output << "  \"pre_global_sheath_proxy_patch_potential_spread_v\": "
           << pre_global_solve_patch_potential_spread_v << ",\n";
    output << "  \"patch_potential_spread_v\": " << patch_potential_spread_v << ",\n";
    output << "  \"patch_potential_spread_reduction_v\": " << patch_potential_spread_reduction_v << ",\n";
    output << "  \"patch_potential_spread_reduction_ratio\": " << patch_potential_spread_reduction_ratio << ",\n";
    output << "  \"shared_current_matrix_coupling_active\": "
           << (shared_current_matrix_coupling_active ? "true" : "false") << ",\n";
    output << "  \"shared_current_matrix_coupling_offdiag_entry_count\": "
           << shared_current_matrix_coupling_offdiag_entry_count << ",\n";
    output << "  \"shared_global_coupled_solve_active\": "
           << (shared_global_coupled_solve_active ? "true" : "false") << ",\n";
    output << "  \"shared_global_coupled_solve_iterations\": "
           << shared_global_coupled_solve_iterations << ",\n";
    output << "  \"shared_global_coupled_solve_converged\": "
           << (shared_global_coupled_solve_converged ? "true" : "false") << ",\n";
    output << "  \"shared_global_coupled_solve_max_delta_v\": "
           << shared_global_coupled_solve_max_delta_v << ",\n";
    output << "  \"shared_live_pic_coupled_refresh_active\": "
           << (shared_live_pic_coupled_refresh_active ? "true" : "false") << ",\n";
    output << "  \"shared_live_pic_coupled_refresh_count\": "
           << shared_live_pic_coupled_refresh_count << ",\n";
    output << "  \"shared_particle_transport_coupling_active\": "
           << (shared_particle_transport_coupling_active ? "true" : "false") << ",\n";
    output << "  \"shared_particle_transport_offdiag_entry_count\": "
           << shared_particle_transport_offdiag_entry_count << ",\n";
    output << "  \"shared_particle_transport_total_conductance_s\": "
           << shared_particle_transport_total_conductance_s << ",\n";
    output << "  \"shared_particle_transport_conservation_error_a_per_v\": "
           << shared_particle_transport_conservation_error_a_per_v << ",\n";
    output << "  \"shared_particle_transport_charge_c\": "
           << shared_particle_transport_charge_c << ",\n";
    output << "  \"shared_particle_transport_reference_shift_v\": "
           << shared_particle_transport_reference_shift_v << ",\n";
    output << "  \"shared_particle_transport_distribution_active\": "
           << (shared_particle_transport_distribution_active ? "true" : "false") << ",\n";
    output << "  \"shared_particle_transport_distribution_charge_spread_c\": "
           << distributed_particle_transport_charge_spread_c << ",\n";
    output << "  \"shared_particle_transport_distribution_reference_shift_spread_v\": "
           << distributed_particle_transport_reference_shift_spread_v << ",\n";
    output << "  \"shared_particle_transport_distribution_conservation_error_c\": "
           << shared_particle_transport_distribution_conservation_error_c << ",\n";
    output << "  \"shared_particle_transport_exchange_active\": "
           << (shared_particle_transport_exchange_active ? "true" : "false") << ",\n";
    output << "  \"shared_particle_transport_exchange_flux_spread_a\": "
           << distributed_particle_transport_net_flux_spread_a << ",\n";
    output << "  \"shared_particle_transport_exchange_flux_conservation_error_a\": "
           << shared_particle_transport_exchange_flux_conservation_error_a << ",\n";
    output << "  \"shared_particle_transport_edge_domain_active\": "
           << (shared_particle_transport_edge_domain_active ? "true" : "false") << ",\n";
    output << "  \"shared_particle_transport_edge_charge_total_abs_c\": "
           << total_abs_edge_transport_charge_c << ",\n";
    output << "  \"shared_particle_transport_edge_charge_spread_c\": "
           << node_edge_transport_charge_spread_c << ",\n";
    output << "  \"shared_particle_transport_edge_charge_conservation_error_c\": "
           << shared_particle_transport_edge_charge_conservation_error_c << ",\n";
    output << "  \"shared_particle_transport_edge_operator_active\": "
           << (shared_particle_transport_edge_operator_active ? "true" : "false") << ",\n";
    output << "  \"shared_particle_transport_edge_target_charge_total_abs_c\": "
           << total_abs_edge_target_charge_c << ",\n";
    output << "  \"shared_particle_transport_edge_operator_total_abs_drive_charge_c\": "
           << total_abs_edge_operator_drive_charge_c << ",\n";
    output << "  \"shared_particle_transport_edge_operator_drive_conservation_error_c\": "
           << shared_particle_transport_edge_operator_drive_conservation_error_c << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_iterations\": "
           << shared_particle_transport_edge_graph_operator_iterations << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_converged\": "
           << (shared_particle_transport_edge_graph_operator_converged ? "true" : "false") << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_max_balance_residual_c\": "
           << shared_particle_transport_edge_graph_operator_max_balance_residual_c << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_branch_graph_active\": "
           << (shared_particle_transport_edge_graph_operator_branch_graph_edge_count > 0.0 ? "true" : "false")
           << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_branch_graph_edge_count\": "
           << shared_particle_transport_edge_graph_operator_branch_graph_edge_count << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_branch_graph_pair_count\": "
           << shared_particle_transport_edge_graph_operator_branch_graph_pair_count << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_effective_pair_count\": "
           << shared_particle_transport_edge_graph_operator_effective_pair_count << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_total_pair_weight_f\": "
           << shared_particle_transport_edge_graph_operator_total_pair_weight_f << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_total_conductance_weight_f\": "
           << std::scientific << std::setprecision(12)
           << shared_particle_transport_edge_graph_operator_total_conductance_weight_f << ",\n"
           << std::fixed << std::setprecision(12);
    output << "  \"shared_particle_transport_edge_graph_operator_min_node_preconditioner\": "
           << shared_particle_transport_edge_graph_operator_min_node_preconditioner << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_max_node_preconditioner\": "
           << shared_particle_transport_edge_graph_operator_max_node_preconditioner << ",\n";
    output << "  \"global_sheath_proxy_solve_weight\": " << shared_global_solve_weight << ",\n";
    output << "  \"global_sheath_proxy_solve_active\": "
           << (global_sheath_proxy_solve_active ? "true" : "false") << ",\n";
    output << "  \"live_pic_electron_current_spread_a_per_m2\": "
           << live_pic_electron_current_spread_a_per_m2 << ",\n";
    output << "  \"live_pic_ion_current_spread_a_per_m2\": "
           << live_pic_ion_current_spread_a_per_m2 << ",\n";
    output << "  \"live_pic_net_current_spread_a_per_m2\": "
           << live_pic_net_current_spread_a_per_m2 << ",\n";
    output << "  \"live_pic_collision_count_spread\": " << live_pic_collision_count_spread << ",\n";
    output << "  \"electron_calibration_factor_spread\": " << electron_calibration_factor_spread << ",\n";
    output << "  \"ion_calibration_factor_spread\": " << ion_calibration_factor_spread << ",\n";
    output << "  \"normal_electric_field_spread_v_per_m\": "
           << normal_electric_field_spread_v_per_m << ",\n";
    output << "  \"local_charge_density_spread_c_per_m3\": " << local_charge_density_spread_c_per_m3 << ",\n";
    output << "  \"consistency_threshold_reference_v\": 0.000000001000,\n";
    output << "  \"consistency_threshold_sheath_charge_c\": 0.000000000000001000,\n";
    output << "  \"consistency_threshold_live_pic_current_a_per_m2\": 0.000000000001,\n";
    output << "  \"consistency_threshold_live_pic_collision_count\": 0.500000000000,\n";
    output << "  \"consistency_threshold_calibration_factor\": 0.000000000001,\n";
    output << "  \"sheath_consistency_pass\": " << (sheath_consistency_pass ? "true" : "false") << ",\n";
    output << "  \"shared_particle_pool_consistency_pass\": "
           << (shared_particle_pool_consistency_pass ? "true" : "false") << ",\n";
    output << "  \"consistency_pass\": " << (consistency_pass ? "true" : "false") << ",\n";
    output << "  \"shared_patch_nodes\": [\n";
    for (std::size_t i = 0; i < shared_patch_nodes.size(); ++i)
    {
        const auto node_index = shared_patch_nodes[i];
        output << "    {\"index\": " << node_index << ", \"name\": \""
               << jsonEscape(circuit_model->nodeName(node_index))
               << "\", \"latest_patch_potential_v\": "
               << latestHistoryValue(node_potential_history, node_index)
               << ", \"latest_shared_reference_potential_v\": "
               << latestHistoryValue(node_shared_reference_potential_history, node_index)
               << ", \"latest_shared_sheath_charge_c\": "
               << latestHistoryValue(node_shared_sheath_charge_history, node_index)
               << ", \"latest_live_pic_net_current_a_per_m2\": "
               << latestHistoryValue(node_live_pic_net_current_history, node_index)
               << ", \"latest_live_pic_collision_count\": "
               << latestHistoryValue(node_live_pic_collision_count_history, node_index)
               << ", \"latest_electron_calibration_factor\": "
               << latestHistoryValue(node_electron_calibration_factor_history, node_index)
               << ", \"latest_ion_calibration_factor\": "
               << latestHistoryValue(node_ion_calibration_factor_history, node_index) << "}";
        output << (i + 1 < shared_patch_nodes.size() ? ",\n" : "\n");
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeSharedSurfaceParticleTransportDomainJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const GlobalParticleDomainState& global_particle_domain_state,
    const std::vector<double>& time_history,
    const std::vector<double>& shared_particle_transport_charge_history,
    const std::vector<double>& shared_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_shared_runtime_enabled_history,
    const std::vector<std::vector<double>>& node_shared_reference_potential_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_charge_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_net_flux_history,
    const std::vector<std::vector<double>>& edge_charge_matrix_c,
    const std::vector<std::vector<double>>& edge_target_charge_matrix_c,
    const std::vector<std::vector<double>>& edge_operator_drive_matrix_c,
    double edge_graph_operator_iterations, bool edge_graph_operator_converged,
    double edge_graph_operator_max_balance_residual_c,
    double edge_graph_operator_branch_graph_edge_count,
    double edge_graph_operator_branch_graph_pair_count,
    double edge_graph_operator_effective_pair_count,
    double edge_graph_operator_total_pair_weight_f,
    double edge_graph_operator_total_conductance_weight_f,
    double edge_graph_operator_min_node_preconditioner,
    double edge_graph_operator_max_node_preconditioner,
    const std::vector<std::vector<double>>& exchange_flux_matrix_a)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };

    std::vector<std::size_t> shared_patch_nodes;
    double total_distributed_charge_c = 0.0;
    double total_exchange_flux_a = 0.0;
    double total_abs_edge_charge_c = 0.0;
    double total_node_edge_charge_c = 0.0;
    double total_abs_edge_target_charge_c = 0.0;
    double total_abs_edge_operator_drive_charge_c = 0.0;
    double total_node_edge_operator_drive_charge_c = 0.0;
    std::size_t exchange_edge_count = 0;
    std::vector<std::tuple<std::size_t, std::size_t, double, double, double, double>> exchange_edges;
    const bool use_owned_runtime_state =
        global_particle_domain_state.active && !global_particle_domain_state.nodes.empty();

    if (use_owned_runtime_state)
    {
        shared_patch_nodes.reserve(global_particle_domain_state.nodes.size());
        for (const auto& node_state : global_particle_domain_state.nodes)
        {
            shared_patch_nodes.push_back(node_state.node_index);
            total_distributed_charge_c += node_state.charge_c;
            total_exchange_flux_a += node_state.net_flux_a;
            double node_edge_stored_charge_c = 0.0;
            double node_edge_target_charge_c = 0.0;
            double node_edge_operator_drive_charge_c = 0.0;
            for (const auto& edge_state : global_particle_domain_state.edges)
            {
                if (edge_state.from_node_index == node_state.node_index)
                {
                    node_edge_stored_charge_c += edge_state.stored_charge_c;
                    node_edge_target_charge_c += edge_state.target_charge_c;
                    node_edge_operator_drive_charge_c += edge_state.operator_drive_charge_c;
                }
                else if (edge_state.to_node_index == node_state.node_index)
                {
                    node_edge_stored_charge_c -= edge_state.stored_charge_c;
                    node_edge_target_charge_c -= edge_state.target_charge_c;
                    node_edge_operator_drive_charge_c -= edge_state.operator_drive_charge_c;
                }
            }
            total_node_edge_charge_c += node_edge_stored_charge_c;
            total_node_edge_operator_drive_charge_c += node_edge_operator_drive_charge_c;
        }
        exchange_edges.reserve(global_particle_domain_state.edges.size());
        for (const auto& edge_state : global_particle_domain_state.edges)
        {
            exchange_edges.emplace_back(edge_state.from_node_index, edge_state.to_node_index,
                                        edge_state.exchange_flux_a, edge_state.stored_charge_c,
                                        edge_state.target_charge_c,
                                        edge_state.operator_drive_charge_c);
            total_abs_edge_charge_c += std::abs(edge_state.stored_charge_c);
            total_abs_edge_target_charge_c += std::abs(edge_state.target_charge_c);
            total_abs_edge_operator_drive_charge_c +=
                std::abs(edge_state.operator_drive_charge_c);
        }
    }
    else if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            if (!circuit_model->nodeIsPatch(node_index))
            {
                continue;
            }
            if (latestHistoryValue(node_shared_runtime_enabled_history, node_index) < 0.5)
            {
                continue;
            }

            shared_patch_nodes.push_back(node_index);
            total_distributed_charge_c +=
                latestHistoryValue(node_distributed_particle_transport_charge_history, node_index);
            total_exchange_flux_a +=
                latestHistoryValue(node_distributed_particle_transport_net_flux_history, node_index);
        }

        for (std::size_t row = 0; row < shared_patch_nodes.size(); ++row)
        {
            const auto from_index = shared_patch_nodes[row];
            double node_edge_stored_charge_c = 0.0;
            if (from_index < edge_charge_matrix_c.size())
            {
                for (const double edge_charge_c : edge_charge_matrix_c[from_index])
                {
                    node_edge_stored_charge_c += edge_charge_c;
                }
            }
            total_node_edge_charge_c += node_edge_stored_charge_c;
            double node_edge_target_charge_c = 0.0;
            if (from_index < edge_target_charge_matrix_c.size())
            {
                for (const double target_charge_c : edge_target_charge_matrix_c[from_index])
                {
                    node_edge_target_charge_c += target_charge_c;
                }
            }
            double node_edge_operator_drive_charge_c = 0.0;
            if (from_index < edge_operator_drive_matrix_c.size())
            {
                for (const double operator_drive_charge_c :
                     edge_operator_drive_matrix_c[from_index])
                {
                    node_edge_operator_drive_charge_c += operator_drive_charge_c;
                }
            }
            total_node_edge_operator_drive_charge_c += node_edge_operator_drive_charge_c;
            if (from_index >= exchange_flux_matrix_a.size())
            {
                continue;
            }
            for (std::size_t col = row + 1; col < shared_patch_nodes.size(); ++col)
            {
                const auto to_index = shared_patch_nodes[col];
                if (to_index >= exchange_flux_matrix_a[from_index].size())
                {
                    continue;
                }
                const double exchange_flux_a = exchange_flux_matrix_a[from_index][to_index];
                const double edge_charge_c =
                    (from_index < edge_charge_matrix_c.size() &&
                     to_index < edge_charge_matrix_c[from_index].size())
                        ? edge_charge_matrix_c[from_index][to_index]
                        : 0.0;
                const double target_edge_charge_c =
                    (from_index < edge_target_charge_matrix_c.size() &&
                     to_index < edge_target_charge_matrix_c[from_index].size())
                        ? edge_target_charge_matrix_c[from_index][to_index]
                        : 0.0;
                const double operator_drive_charge_c =
                    (from_index < edge_operator_drive_matrix_c.size() &&
                     to_index < edge_operator_drive_matrix_c[from_index].size())
                        ? edge_operator_drive_matrix_c[from_index][to_index]
                        : 0.0;
                if ((!std::isfinite(exchange_flux_a) || std::abs(exchange_flux_a) <= 0.0) &&
                    (!std::isfinite(edge_charge_c) || std::abs(edge_charge_c) <= 0.0) &&
                    (!std::isfinite(operator_drive_charge_c) ||
                     std::abs(operator_drive_charge_c) <= 0.0))
                {
                    continue;
                }
                exchange_edges.emplace_back(from_index, to_index, exchange_flux_a, edge_charge_c,
                                            target_edge_charge_c, operator_drive_charge_c);
                total_abs_edge_charge_c += std::abs(edge_charge_c);
                total_abs_edge_target_charge_c += std::abs(target_edge_charge_c);
                total_abs_edge_operator_drive_charge_c += std::abs(operator_drive_charge_c);
            }
        }
    }

    exchange_edge_count = exchange_edges.size();
    const double shared_particle_transport_charge_c = use_owned_runtime_state
        ? global_particle_domain_state.total_charge_c
        : latest_scalar(shared_particle_transport_charge_history);
    const double shared_particle_transport_reference_shift_v = use_owned_runtime_state
        ? global_particle_domain_state.total_reference_shift_v
        : latest_scalar(shared_particle_transport_reference_shift_history);
    const double distributed_transport_charge_conservation_error_c = use_owned_runtime_state
        ? global_particle_domain_state.charge_conservation_error_c
        : std::abs(total_distributed_charge_c - shared_particle_transport_charge_c);
    const double exchange_flux_conservation_error_a = std::abs(total_exchange_flux_a);
    const double edge_charge_conservation_error_c = std::abs(total_node_edge_charge_c);
    const double edge_operator_drive_conservation_error_c =
        std::abs(total_node_edge_operator_drive_charge_c);
    const bool domain_active = use_owned_runtime_state
        ? global_particle_domain_state.active
        : (shared_patch_nodes.size() >= 2 &&
           (std::abs(shared_particle_transport_charge_c) > 0.0 ||
            exchange_edge_count > 0 || total_abs_edge_charge_c > 0.0));

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_shared_particle_transport_domain.v1\",\n";
    output << "  \"contract_id\": \"surface-pic-shared-particle-transport-domain-v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_runtime\": \"" << surfacePicRuntimeName(config.surface_pic_runtime_kind)
           << "\",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << (time_history.empty() ? 0.0 : time_history.back()) << ",\n";
    output << "  \"runtime_state_backed\": true,\n";
    output << "  \"shared_particle_transport_bookkeeping_mode\": "
           << "\"" << jsonEscape(global_particle_domain_state.bookkeeping_mode) << "\",\n";
    output << "  \"shared_particle_transport_domain_active\": "
           << (domain_active ? "true" : "false") << ",\n";
    output << "  \"shared_patch_node_count\": " << shared_patch_nodes.size() << ",\n";
    output << "  \"exchange_edge_count\": " << exchange_edge_count << ",\n";
    output << "  \"shared_particle_transport_charge_c\": "
           << shared_particle_transport_charge_c << ",\n";
    output << "  \"shared_particle_transport_reference_shift_v\": "
           << shared_particle_transport_reference_shift_v << ",\n";
    output << "  \"distributed_particle_transport_charge_total_c\": "
           << total_distributed_charge_c << ",\n";
    output << "  \"distributed_particle_transport_charge_conservation_error_c\": "
           << distributed_transport_charge_conservation_error_c << ",\n";
    output << "  \"exchange_flux_conservation_error_a\": "
           << exchange_flux_conservation_error_a << ",\n";
    output << "  \"edge_charge_total_abs_c\": " << total_abs_edge_charge_c << ",\n";
    output << "  \"edge_charge_conservation_error_c\": "
           << edge_charge_conservation_error_c << ",\n";
    output << "  \"edge_target_charge_total_abs_c\": "
           << total_abs_edge_target_charge_c << ",\n";
    output << "  \"edge_operator_drive_total_abs_c\": "
           << total_abs_edge_operator_drive_charge_c << ",\n";
    output << "  \"edge_operator_drive_conservation_error_c\": "
           << edge_operator_drive_conservation_error_c << ",\n";
    output << "  \"edge_graph_operator_iterations\": "
           << edge_graph_operator_iterations << ",\n";
    output << "  \"edge_graph_operator_converged\": "
           << (edge_graph_operator_converged ? "true" : "false") << ",\n";
    output << "  \"edge_graph_operator_max_balance_residual_c\": "
           << edge_graph_operator_max_balance_residual_c << ",\n";
    output << "  \"edge_graph_operator_branch_graph_active\": "
           << (edge_graph_operator_branch_graph_edge_count > 0.0 ? "true" : "false")
           << ",\n";
    output << "  \"edge_graph_operator_branch_graph_edge_count\": "
           << edge_graph_operator_branch_graph_edge_count << ",\n";
    output << "  \"edge_graph_operator_branch_graph_pair_count\": "
           << edge_graph_operator_branch_graph_pair_count << ",\n";
    output << "  \"edge_graph_operator_effective_pair_count\": "
           << edge_graph_operator_effective_pair_count << ",\n";
    output << "  \"edge_graph_operator_total_pair_weight_f\": "
           << edge_graph_operator_total_pair_weight_f << ",\n";
    output << "  \"edge_graph_operator_total_conductance_weight_f\": "
           << std::scientific << std::setprecision(12)
           << edge_graph_operator_total_conductance_weight_f << ",\n"
           << std::fixed << std::setprecision(12);
    output << "  \"edge_graph_operator_min_node_preconditioner\": "
           << edge_graph_operator_min_node_preconditioner << ",\n";
    output << "  \"edge_graph_operator_max_node_preconditioner\": "
           << edge_graph_operator_max_node_preconditioner << ",\n";
    output << "  \"nodes\": [\n";
    for (std::size_t i = 0; i < shared_patch_nodes.size(); ++i)
    {
        const auto node_index = shared_patch_nodes[i];
        double node_edge_stored_charge_c = 0.0;
        double node_edge_target_charge_c = 0.0;
        double node_edge_operator_drive_charge_c = 0.0;
        double latest_patch_potential_v = latestHistoryValue(node_potential_history, node_index);
        double latest_shared_reference_potential_v =
            latestHistoryValue(node_shared_reference_potential_history, node_index);
        double latest_distributed_particle_transport_charge_c =
            latestHistoryValue(node_distributed_particle_transport_charge_history, node_index);
        double latest_distributed_particle_transport_reference_shift_v = latestHistoryValue(
            node_distributed_particle_transport_reference_shift_history, node_index);
        double latest_distributed_particle_transport_net_flux_a =
            latestHistoryValue(node_distributed_particle_transport_net_flux_history, node_index);
        if (use_owned_runtime_state)
        {
            for (const auto& node_state : global_particle_domain_state.nodes)
            {
                if (node_state.node_index != node_index)
                {
                    continue;
                }
                latest_patch_potential_v = node_state.patch_potential_v;
                latest_shared_reference_potential_v =
                    node_state.shared_reference_potential_v;
                latest_distributed_particle_transport_charge_c = node_state.charge_c;
                latest_distributed_particle_transport_reference_shift_v =
                    node_state.reference_shift_v;
                latest_distributed_particle_transport_net_flux_a = node_state.net_flux_a;
                break;
            }
            for (const auto& edge_state : global_particle_domain_state.edges)
            {
                if (edge_state.from_node_index == node_index)
                {
                    node_edge_stored_charge_c += edge_state.stored_charge_c;
                    node_edge_target_charge_c += edge_state.target_charge_c;
                    node_edge_operator_drive_charge_c += edge_state.operator_drive_charge_c;
                }
                else if (edge_state.to_node_index == node_index)
                {
                    node_edge_stored_charge_c -= edge_state.stored_charge_c;
                    node_edge_target_charge_c -= edge_state.target_charge_c;
                    node_edge_operator_drive_charge_c -= edge_state.operator_drive_charge_c;
                }
            }
        }
        else
        {
            if (node_index < edge_charge_matrix_c.size())
            {
                for (const double edge_charge_c : edge_charge_matrix_c[node_index])
                {
                    node_edge_stored_charge_c += edge_charge_c;
                }
            }
            if (node_index < edge_target_charge_matrix_c.size())
            {
                for (const double target_charge_c : edge_target_charge_matrix_c[node_index])
                {
                    node_edge_target_charge_c += target_charge_c;
                }
            }
            if (node_index < edge_operator_drive_matrix_c.size())
            {
                for (const double operator_drive_charge_c :
                     edge_operator_drive_matrix_c[node_index])
                {
                    node_edge_operator_drive_charge_c += operator_drive_charge_c;
                }
            }
        }
        output << "    {\"index\": " << node_index << ", \"name\": \""
               << jsonEscape(circuit_model->nodeName(node_index))
               << "\", \"latest_patch_potential_v\": "
               << latest_patch_potential_v
               << ", \"latest_shared_reference_potential_v\": "
               << latest_shared_reference_potential_v
               << ", \"latest_distributed_particle_transport_charge_c\": "
               << latest_distributed_particle_transport_charge_c
               << ", \"latest_distributed_particle_transport_reference_shift_v\": "
               << latest_distributed_particle_transport_reference_shift_v
               << ", \"latest_distributed_particle_transport_net_flux_a\": "
               << latest_distributed_particle_transport_net_flux_a
               << ", \"net_edge_stored_charge_c\": " << node_edge_stored_charge_c
               << ", \"net_edge_target_charge_c\": " << node_edge_target_charge_c
               << ", \"net_edge_operator_drive_charge_c\": "
               << node_edge_operator_drive_charge_c
               << "}";
        output << (i + 1 < shared_patch_nodes.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"exchange_edges\": [\n";
    for (std::size_t i = 0; i < exchange_edges.size(); ++i)
    {
        const auto [from_index, to_index, exchange_flux_a, edge_charge_c, target_edge_charge_c,
                    operator_drive_charge_c] = exchange_edges[i];
        output << "    {\"from_index\": " << from_index << ", \"from_name\": \""
               << jsonEscape(circuit_model->nodeName(from_index))
               << "\", \"to_index\": " << to_index << ", \"to_name\": \""
               << jsonEscape(circuit_model->nodeName(to_index))
               << "\", \"net_flux_a\": " << exchange_flux_a
               << ", \"abs_flux_a\": " << std::abs(exchange_flux_a)
               << ", \"stored_charge_c\": " << edge_charge_c
               << ", \"abs_stored_charge_c\": " << std::abs(edge_charge_c)
               << ", \"target_charge_c\": " << target_edge_charge_c
               << ", \"abs_target_charge_c\": " << std::abs(target_edge_charge_c)
               << ", \"operator_drive_charge_c\": " << operator_drive_charge_c
               << ", \"abs_operator_drive_charge_c\": "
               << std::abs(operator_drive_charge_c) << "}";
        output << (i + 1 < exchange_edges.size() ? ",\n" : "\n");
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeGlobalSurfaceParticleDomainJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model, const std::vector<double>& time_history,
    const GlobalParticleDomainState& global_particle_domain_state,
    const std::vector<double>& shared_particle_transport_charge_history,
    const std::vector<double>& shared_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_shared_runtime_enabled_history,
    const std::vector<std::vector<double>>& node_shared_reference_potential_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_charge_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_net_flux_history,
    const std::vector<std::vector<double>>& edge_charge_matrix_c,
    const std::vector<std::vector<double>>& edge_target_charge_matrix_c,
    const std::vector<std::vector<double>>& edge_operator_drive_matrix_c,
    double edge_graph_operator_iterations, bool edge_graph_operator_converged,
    double edge_graph_operator_max_balance_residual_c,
    double edge_graph_operator_branch_graph_edge_count,
    double edge_graph_operator_branch_graph_pair_count,
    double edge_graph_operator_effective_pair_count,
    double edge_graph_operator_total_pair_weight_f,
    double edge_graph_operator_total_conductance_weight_f,
    double edge_graph_operator_min_node_preconditioner,
    double edge_graph_operator_max_node_preconditioner,
    const std::vector<std::vector<double>>& exchange_flux_matrix_a)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };

    std::vector<std::size_t> shared_patch_nodes;
    std::vector<std::tuple<std::size_t, std::size_t, double, double, double, double, double>>
        domain_edges;
    double total_node_charge_c = 0.0;
    double total_node_flux_a = 0.0;
    double total_abs_edge_charge_c = 0.0;
    double total_abs_edge_target_charge_c = 0.0;
    double total_abs_edge_operator_drive_charge_c = 0.0;
    double total_edge_conductance_s = 0.0;
    const bool use_owned_runtime_state =
        global_particle_domain_state.active && !global_particle_domain_state.nodes.empty();

    if (use_owned_runtime_state)
    {
        shared_patch_nodes.reserve(global_particle_domain_state.nodes.size());
        for (const auto& node_state : global_particle_domain_state.nodes)
        {
            shared_patch_nodes.push_back(node_state.node_index);
            total_node_charge_c += node_state.charge_c;
            total_node_flux_a += node_state.net_flux_a;
        }
        total_abs_edge_charge_c = global_particle_domain_state.edge_charge_total_abs_c;
        total_abs_edge_target_charge_c =
            global_particle_domain_state.edge_target_charge_total_abs_c;
        total_abs_edge_operator_drive_charge_c =
            global_particle_domain_state.edge_operator_drive_total_abs_c;
        total_edge_conductance_s = global_particle_domain_state.edge_conductance_total_s;
    }
    else if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            if (!circuit_model->nodeIsPatch(node_index))
            {
                continue;
            }
            if (latestHistoryValue(node_shared_runtime_enabled_history, node_index) < 0.5)
            {
                continue;
            }

            shared_patch_nodes.push_back(node_index);
            total_node_charge_c +=
                latestHistoryValue(node_distributed_particle_transport_charge_history, node_index);
            total_node_flux_a +=
                latestHistoryValue(node_distributed_particle_transport_net_flux_history, node_index);
        }

        for (std::size_t row = 0; row < shared_patch_nodes.size(); ++row)
        {
            const auto from_index = shared_patch_nodes[row];
            if (from_index >= exchange_flux_matrix_a.size())
            {
                continue;
            }
            for (std::size_t col = row + 1; col < shared_patch_nodes.size(); ++col)
            {
                const auto to_index = shared_patch_nodes[col];
                if (to_index >= exchange_flux_matrix_a[from_index].size())
                {
                    continue;
                }
                const double exchange_flux_a = exchange_flux_matrix_a[from_index][to_index];
                const double edge_charge_c =
                    (from_index < edge_charge_matrix_c.size() &&
                     to_index < edge_charge_matrix_c[from_index].size())
                        ? edge_charge_matrix_c[from_index][to_index]
                        : 0.0;
                const double target_charge_c =
                    (from_index < edge_target_charge_matrix_c.size() &&
                     to_index < edge_target_charge_matrix_c[from_index].size())
                        ? edge_target_charge_matrix_c[from_index][to_index]
                        : 0.0;
                const double operator_drive_charge_c =
                    (from_index < edge_operator_drive_matrix_c.size() &&
                     to_index < edge_operator_drive_matrix_c[from_index].size())
                        ? edge_operator_drive_matrix_c[from_index][to_index]
                        : 0.0;
                if ((!std::isfinite(exchange_flux_a) || std::abs(exchange_flux_a) <= 0.0) &&
                    (!std::isfinite(edge_charge_c) || std::abs(edge_charge_c) <= 0.0) &&
                    (!std::isfinite(operator_drive_charge_c) ||
                     std::abs(operator_drive_charge_c) <= 0.0))
                {
                    continue;
                }
                domain_edges.emplace_back(from_index, to_index, 0.0, exchange_flux_a,
                                          edge_charge_c, target_charge_c,
                                          operator_drive_charge_c);
                total_abs_edge_charge_c += std::abs(edge_charge_c);
                total_abs_edge_target_charge_c += std::abs(target_charge_c);
                total_abs_edge_operator_drive_charge_c += std::abs(operator_drive_charge_c);
            }
        }
    }

    const double global_particle_charge_c =
        global_particle_domain_state.active ? global_particle_domain_state.total_charge_c
                                            : latest_scalar(shared_particle_transport_charge_history);
    const double global_particle_reference_shift_v =
        global_particle_domain_state.active
            ? global_particle_domain_state.total_reference_shift_v
            : latest_scalar(shared_particle_transport_reference_shift_history);
    const double global_particle_charge_conservation_error_c =
        global_particle_domain_state.active
            ? global_particle_domain_state.charge_conservation_error_c
            : std::abs(total_node_charge_c - global_particle_charge_c);
    const double global_particle_flux_conservation_error_a =
        global_particle_domain_state.active
            ? global_particle_domain_state.flux_conservation_error_a
            : std::abs(total_node_flux_a);
    const bool domain_active = global_particle_domain_state.active ||
                               (shared_patch_nodes.size() >= 2 &&
                                (std::abs(global_particle_charge_c) > 0.0 || !domain_edges.empty() ||
                                 total_abs_edge_charge_c > 0.0));

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_global_particle_domain.v1\",\n";
    output << "  \"contract_id\": \"surface-pic-global-particle-domain-v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_runtime\": \"" << surfacePicRuntimeName(config.surface_pic_runtime_kind)
           << "\",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << latest_scalar(time_history) << ",\n";
    output << "  \"global_particle_domain_active\": "
           << (domain_active ? "true" : "false") << ",\n";
    output << "  \"global_particle_bookkeeping_mode\": "
           << "\"" << jsonEscape(global_particle_domain_state.bookkeeping_mode) << "\",\n";
    output << "  \"runtime_state_backed\": true,\n";
    output << "  \"shared_patch_node_count\": " << shared_patch_nodes.size() << ",\n";
    output << "  \"domain_edge_count\": "
           << (use_owned_runtime_state ? global_particle_domain_state.edges.size()
                                       : domain_edges.size())
           << ",\n";
    output << "  \"global_particle_charge_c\": " << global_particle_charge_c << ",\n";
    output << "  \"global_particle_reference_shift_v\": "
           << global_particle_reference_shift_v << ",\n";
    output << "  \"global_particle_charge_conservation_error_c\": "
           << global_particle_charge_conservation_error_c << ",\n";
    output << "  \"global_particle_flux_conservation_error_a\": "
           << global_particle_flux_conservation_error_a << ",\n";
    output << "  \"global_particle_edge_charge_total_abs_c\": "
           << total_abs_edge_charge_c << ",\n";
    output << "  \"global_particle_edge_target_charge_total_abs_c\": "
           << total_abs_edge_target_charge_c << ",\n";
    output << "  \"global_particle_edge_operator_drive_total_abs_c\": "
           << total_abs_edge_operator_drive_charge_c << ",\n";
    output << "  \"global_particle_edge_conductance_total_s\": "
           << total_edge_conductance_s << ",\n";
    output << "  \"edge_graph_operator_iterations\": " << edge_graph_operator_iterations << ",\n";
    output << "  \"edge_graph_operator_converged\": "
           << (edge_graph_operator_converged ? "true" : "false") << ",\n";
    output << "  \"edge_graph_operator_max_balance_residual_c\": "
           << edge_graph_operator_max_balance_residual_c << ",\n";
    output << "  \"edge_graph_operator_branch_graph_active\": "
           << (edge_graph_operator_branch_graph_edge_count > 0.0 ? "true" : "false")
           << ",\n";
    output << "  \"edge_graph_operator_branch_graph_edge_count\": "
           << edge_graph_operator_branch_graph_edge_count << ",\n";
    output << "  \"edge_graph_operator_branch_graph_pair_count\": "
           << edge_graph_operator_branch_graph_pair_count << ",\n";
    output << "  \"edge_graph_operator_effective_pair_count\": "
           << edge_graph_operator_effective_pair_count << ",\n";
    output << "  \"edge_graph_operator_total_pair_weight_f\": "
           << edge_graph_operator_total_pair_weight_f << ",\n";
    output << "  \"edge_graph_operator_total_conductance_weight_f\": "
           << std::scientific << std::setprecision(12)
           << edge_graph_operator_total_conductance_weight_f << ",\n"
           << std::fixed << std::setprecision(12);
    output << "  \"edge_graph_operator_min_node_preconditioner\": "
           << edge_graph_operator_min_node_preconditioner << ",\n";
    output << "  \"edge_graph_operator_max_node_preconditioner\": "
           << edge_graph_operator_max_node_preconditioner << ",\n";
    output << "  \"nodes\": [\n";
    if (use_owned_runtime_state)
    {
        for (std::size_t i = 0; i < global_particle_domain_state.nodes.size(); ++i)
        {
            const auto& node_state = global_particle_domain_state.nodes[i];
            output << "    {\"index\": " << node_state.node_index << ", \"name\": \""
                   << jsonEscape(node_state.node_name)
                   << "\", \"latest_patch_potential_v\": "
                   << node_state.patch_potential_v
                   << ", \"latest_shared_reference_potential_v\": "
                   << node_state.shared_reference_potential_v
                   << ", \"global_node_charge_c\": " << node_state.charge_c
                   << ", \"global_node_reference_shift_v\": "
                   << node_state.reference_shift_v
                   << ", \"global_node_flux_a\": " << node_state.net_flux_a << "}";
            output << (i + 1 < global_particle_domain_state.nodes.size() ? ",\n" : "\n");
        }
    }
    else
    {
        for (std::size_t i = 0; i < shared_patch_nodes.size(); ++i)
        {
            const auto node_index = shared_patch_nodes[i];
            output << "    {\"index\": " << node_index << ", \"name\": \""
                   << jsonEscape(circuit_model->nodeName(node_index))
                   << "\", \"latest_patch_potential_v\": "
                   << latestHistoryValue(node_potential_history, node_index)
                   << ", \"latest_shared_reference_potential_v\": "
                   << latestHistoryValue(node_shared_reference_potential_history, node_index)
                   << ", \"global_node_charge_c\": "
                   << latestHistoryValue(node_distributed_particle_transport_charge_history,
                                         node_index)
                   << ", \"global_node_reference_shift_v\": "
                   << latestHistoryValue(node_distributed_particle_transport_reference_shift_history,
                                         node_index)
                   << ", \"global_node_flux_a\": "
                   << latestHistoryValue(node_distributed_particle_transport_net_flux_history,
                                         node_index)
                   << "}";
            output << (i + 1 < shared_patch_nodes.size() ? ",\n" : "\n");
        }
    }
    output << "  ],\n";
    output << "  \"domain_edges\": [\n";
    if (use_owned_runtime_state)
    {
        for (std::size_t i = 0; i < global_particle_domain_state.edges.size(); ++i)
        {
            const auto& edge_state = global_particle_domain_state.edges[i];
            output << "    {\"from_index\": " << edge_state.from_node_index
                   << ", \"from_name\": \"" << jsonEscape(edge_state.from_node_name)
                   << "\", \"to_index\": " << edge_state.to_node_index
                   << ", \"to_name\": \"" << jsonEscape(edge_state.to_node_name)
                   << "\", \"conductance_s\": " << edge_state.conductance_s
                   << ", \"net_flux_a\": " << edge_state.exchange_flux_a
                   << ", \"stored_charge_c\": " << edge_state.stored_charge_c
                   << ", \"target_charge_c\": " << edge_state.target_charge_c
                   << ", \"operator_drive_charge_c\": "
                   << edge_state.operator_drive_charge_c << "}";
            output << (i + 1 < global_particle_domain_state.edges.size() ? ",\n" : "\n");
        }
    }
    else
    {
        for (std::size_t i = 0; i < domain_edges.size(); ++i)
        {
            const auto [from_index, to_index, conductance_s, exchange_flux_a, edge_charge_c,
                        target_charge_c, operator_drive_charge_c] = domain_edges[i];
            output << "    {\"from_index\": " << from_index << ", \"from_name\": \""
                   << jsonEscape(circuit_model->nodeName(from_index))
                   << "\", \"to_index\": " << to_index << ", \"to_name\": \""
                   << jsonEscape(circuit_model->nodeName(to_index))
                   << "\", \"conductance_s\": " << conductance_s
                   << ", \"net_flux_a\": " << exchange_flux_a
                   << ", \"stored_charge_c\": " << edge_charge_c
                   << ", \"target_charge_c\": " << target_charge_c
                   << ", \"operator_drive_charge_c\": " << operator_drive_charge_c << "}";
            output << (i + 1 < domain_edges.size() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeGlobalParticleRepositoryJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const std::vector<double>& time_history,
    const GlobalParticleRepositoryState& global_particle_repository_state)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_global_particle_repository.v1\",\n";
    output << "  \"contract_id\": \"surface-pic-global-particle-repository-v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_runtime\": \"" << surfacePicRuntimeName(config.surface_pic_runtime_kind)
           << "\",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << latest_scalar(time_history) << ",\n";
    output << "  \"runtime_state_backed\": true,\n";
    output << "  \"global_particle_repository_active\": "
           << (global_particle_repository_state.active ? "true" : "false") << ",\n";
    output << "  \"global_particle_repository_bookkeeping_mode\": "
           << "\"" << jsonEscape(global_particle_repository_state.bookkeeping_mode) << "\",\n";
    output << "  \"global_particle_repository_lifecycle_mode\": "
           << "\"" << jsonEscape(global_particle_repository_state.lifecycle_mode) << "\",\n";
    output << "  \"shared_patch_node_count\": "
           << global_particle_repository_state.nodes.size() << ",\n";
    output << "  \"repository_edge_count\": "
           << global_particle_repository_state.edges.size() << ",\n";
    output << "  \"total_reservoir_charge_c\": "
           << global_particle_repository_state.total_reservoir_charge_c << ",\n";
    output << "  \"total_target_reservoir_charge_c\": "
           << global_particle_repository_state.total_target_reservoir_charge_c << ",\n";
    output << "  \"total_migration_delta_abs_charge_c\": "
           << global_particle_repository_state.total_migration_delta_abs_charge_c << ",\n";
    output << "  \"total_edge_feedback_abs_charge_c\": "
           << global_particle_repository_state.total_edge_feedback_abs_charge_c << ",\n";
    output << "  \"total_conservation_correction_abs_charge_c\": "
           << global_particle_repository_state.total_conservation_correction_abs_charge_c << ",\n";
    output << "  \"total_migration_edge_abs_charge_c\": "
           << global_particle_repository_state.total_migration_edge_abs_charge_c << ",\n";
    output << "  \"charge_conservation_error_c\": "
           << global_particle_repository_state.charge_conservation_error_c << ",\n";
    output << "  \"migration_charge_conservation_error_c\": "
           << global_particle_repository_state.migration_charge_conservation_error_c << ",\n";
    output << "  \"nodes\": [\n";
    for (std::size_t i = 0; i < global_particle_repository_state.nodes.size(); ++i)
    {
        const auto& node_state = global_particle_repository_state.nodes[i];
        output << "    {\"index\": " << node_state.node_index << ", \"name\": \""
               << jsonEscape(node_state.node_name)
               << "\", \"area_m2\": " << node_state.area_m2
               << ", \"latest_patch_potential_v\": " << node_state.patch_potential_v
               << ", \"latest_shared_reference_potential_v\": "
               << node_state.shared_reference_potential_v
               << ", \"reference_shift_v\": " << node_state.reference_shift_v
               << ", \"reservoir_charge_c\": " << node_state.reservoir_charge_c
               << ", \"target_reservoir_charge_c\": " << node_state.target_reservoir_charge_c
               << ", \"migration_delta_charge_c\": " << node_state.migration_delta_charge_c
               << ", \"edge_feedback_charge_c\": " << node_state.edge_feedback_charge_c
               << ", \"conservation_correction_charge_c\": "
               << node_state.conservation_correction_charge_c
               << ", \"net_flux_a\": " << node_state.net_flux_a << "}";
        output << (i + 1 < global_particle_repository_state.nodes.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"edges\": [\n";
    for (std::size_t i = 0; i < global_particle_repository_state.edges.size(); ++i)
    {
        const auto& edge_state = global_particle_repository_state.edges[i];
        output << "    {\"from_index\": " << edge_state.from_node_index
               << ", \"from_name\": \"" << jsonEscape(edge_state.from_node_name)
               << "\", \"to_index\": " << edge_state.to_node_index
               << ", \"to_name\": \"" << jsonEscape(edge_state.to_node_name)
               << "\", \"conductance_s\": " << edge_state.conductance_s
               << ", \"migration_flux_a\": " << edge_state.migration_flux_a
               << ", \"migration_charge_c\": " << edge_state.migration_charge_c
               << ", \"stored_charge_c\": " << edge_state.stored_charge_c
               << ", \"target_charge_c\": " << edge_state.target_charge_c
               << ", \"operator_drive_charge_c\": "
               << edge_state.operator_drive_charge_c << "}";
        output << (i + 1 < global_particle_repository_state.edges.size() ? ",\n" : "\n");
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

bool writeGlobalSurfaceSheathFieldSolveJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model, const std::vector<double>& time_history,
    const GlobalSheathFieldSolveState& global_sheath_field_solve_state,
    const GlobalParticleDomainState& global_particle_domain_state,
    const std::vector<double>& shared_effective_sheath_length_history,
    const std::vector<double>& shared_particle_transport_charge_history,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_shared_runtime_enabled_history,
    const std::vector<std::vector<double>>& node_shared_reference_potential_history,
    const std::vector<std::vector<double>>& node_normal_electric_field_history,
    const std::vector<std::vector<double>>& node_local_charge_density_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_charge_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_reference_shift_history,
    const std::vector<std::vector<double>>& node_distributed_particle_transport_net_flux_history,
    const std::vector<double>& shared_global_coupled_solve_active_history,
    const std::vector<double>& shared_global_coupled_solve_iteration_history,
    const std::vector<double>& shared_global_coupled_solve_converged_history,
    const std::vector<double>& shared_global_coupled_solve_max_delta_history)
{
    std::ofstream output(json_path);
    if (!output.is_open())
    {
        return false;
    }

    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };

    std::vector<std::size_t> shared_patch_nodes;
    double total_area_m2 = 0.0;
    double weighted_patch_potential_v = 0.0;
    double weighted_reference_potential_v = 0.0;
    double min_patch_potential_v = 0.0;
    double max_patch_potential_v = 0.0;
    double min_reference_potential_v = 0.0;
    double max_reference_potential_v = 0.0;
    double min_normal_field_v_per_m = 0.0;
    double max_normal_field_v_per_m = 0.0;
    double min_charge_density_c_per_m3 = 0.0;
    double max_charge_density_c_per_m3 = 0.0;
    double total_node_charge_c = 0.0;
    double total_node_flux_a = 0.0;
    bool first_shared_node = true;
    const bool use_owned_runtime_state =
        global_sheath_field_solve_state.active && !global_sheath_field_solve_state.nodes.empty();

    if (use_owned_runtime_state)
    {
        shared_patch_nodes.reserve(global_sheath_field_solve_state.nodes.size());
        for (const auto& node_state : global_sheath_field_solve_state.nodes)
        {
            shared_patch_nodes.push_back(node_state.node_index);
        }
    }
    else if (circuit_model != nullptr)
    {
        for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
        {
            if (!circuit_model->nodeIsPatch(node_index))
            {
                continue;
            }
            if (latestHistoryValue(node_shared_runtime_enabled_history, node_index) < 0.5)
            {
                continue;
            }

            shared_patch_nodes.push_back(node_index);
            const double area_m2 = std::max(1.0e-16, circuit_model->nodeAreaM2(node_index));
            const double patch_potential_v = latestHistoryValue(node_potential_history, node_index);
            const double reference_potential_v =
                latestHistoryValue(node_shared_reference_potential_history, node_index);
            const double normal_field_v_per_m =
                latestHistoryValue(node_normal_electric_field_history, node_index);
            const double charge_density_c_per_m3 =
                latestHistoryValue(node_local_charge_density_history, node_index);
            const double distributed_charge_c =
                latestHistoryValue(node_distributed_particle_transport_charge_history, node_index);
            const double distributed_flux_a = latestHistoryValue(
                node_distributed_particle_transport_net_flux_history, node_index);

            total_area_m2 += area_m2;
            weighted_patch_potential_v += area_m2 * patch_potential_v;
            weighted_reference_potential_v += area_m2 * reference_potential_v;
            total_node_charge_c += distributed_charge_c;
            total_node_flux_a += distributed_flux_a;
            if (first_shared_node)
            {
                min_patch_potential_v = max_patch_potential_v = patch_potential_v;
                min_reference_potential_v = max_reference_potential_v = reference_potential_v;
                min_normal_field_v_per_m = max_normal_field_v_per_m = normal_field_v_per_m;
                min_charge_density_c_per_m3 = max_charge_density_c_per_m3 = charge_density_c_per_m3;
                first_shared_node = false;
            }
            else
            {
                min_patch_potential_v = std::min(min_patch_potential_v, patch_potential_v);
                max_patch_potential_v = std::max(max_patch_potential_v, patch_potential_v);
                min_reference_potential_v =
                    std::min(min_reference_potential_v, reference_potential_v);
                max_reference_potential_v =
                    std::max(max_reference_potential_v, reference_potential_v);
                min_normal_field_v_per_m =
                    std::min(min_normal_field_v_per_m, normal_field_v_per_m);
                max_normal_field_v_per_m =
                    std::max(max_normal_field_v_per_m, normal_field_v_per_m);
                min_charge_density_c_per_m3 =
                    std::min(min_charge_density_c_per_m3, charge_density_c_per_m3);
                max_charge_density_c_per_m3 =
                    std::max(max_charge_density_c_per_m3, charge_density_c_per_m3);
            }
        }
    }

    const bool active =
        global_sheath_field_solve_state.active ||
        (latest_scalar(shared_global_coupled_solve_active_history) >= 0.5 &&
         shared_patch_nodes.size() >= 2 && total_area_m2 > 0.0);
    const double global_patch_potential_v =
        global_sheath_field_solve_state.active
            ? global_sheath_field_solve_state.global_patch_potential_v
            : (total_area_m2 > 0.0 ? weighted_patch_potential_v / total_area_m2 : 0.0);
    const double global_reference_potential_v =
        global_sheath_field_solve_state.active
            ? global_sheath_field_solve_state.global_reference_potential_v
            : (total_area_m2 > 0.0 ? weighted_reference_potential_v / total_area_m2 : 0.0);
    const double global_effective_sheath_length_m = std::max(
        1.0e-6, global_sheath_field_solve_state.active
                    ? global_sheath_field_solve_state.effective_sheath_length_m
                    : (latest_scalar(shared_effective_sheath_length_history) > 0.0
                           ? latest_scalar(shared_effective_sheath_length_history)
                           : std::max(config.minimum_sheath_length_m,
                                      std::min(config.maximum_sheath_length_m,
                                               config.sheath_length_m))));
    const double global_normal_electric_field_v_per_m =
        global_sheath_field_solve_state.active
            ? global_sheath_field_solve_state.global_normal_electric_field_v_per_m
            : (global_reference_potential_v - global_patch_potential_v) /
                  global_effective_sheath_length_m;
    const double global_local_charge_density_c_per_m3 =
        global_sheath_field_solve_state.active
            ? global_sheath_field_solve_state.global_local_charge_density_c_per_m3
            : kEpsilon0 * global_normal_electric_field_v_per_m / global_effective_sheath_length_m;
    const double reference_spread_v = max_reference_potential_v - min_reference_potential_v;
    const double patch_spread_v = max_patch_potential_v - min_patch_potential_v;
    const double normal_field_spread_v_per_m =
        max_normal_field_v_per_m - min_normal_field_v_per_m;
    const double local_charge_density_spread_c_per_m3 =
        max_charge_density_c_per_m3 - min_charge_density_c_per_m3;
    const double field_residual_v_per_m =
        global_sheath_field_solve_state.active
            ? global_sheath_field_solve_state.field_residual_v_per_m
            : std::max({normal_field_spread_v_per_m,
                        std::abs(reference_spread_v) / global_effective_sheath_length_m,
                        std::abs(local_charge_density_spread_c_per_m3) *
                            global_effective_sheath_length_m / std::max(1.0e-12, kEpsilon0)});
    const double global_particle_charge_conservation_error_c =
        std::abs(total_node_charge_c - latest_scalar(shared_particle_transport_charge_history));
    const double global_particle_flux_conservation_error_a = std::abs(total_node_flux_a);
    const double particle_field_coupled_residual_v =
        global_sheath_field_solve_state.active
            ? global_sheath_field_solve_state.particle_field_coupled_residual_v
            : std::max({field_residual_v_per_m * global_effective_sheath_length_m,
                        global_particle_charge_conservation_error_c *
                            global_effective_sheath_length_m /
                            std::max(1.0e-18, kEpsilon0 * total_area_m2),
                        global_particle_flux_conservation_error_a *
                            global_effective_sheath_length_m /
                            std::max(1.0e-18, kEpsilon0 * total_area_m2)});
    double multi_step_stability_metric_v = global_sheath_field_solve_state.active
                                               ? global_sheath_field_solve_state.multi_step_stability_metric_v
                                               : latest_scalar(shared_global_coupled_solve_max_delta_history);
    if (shared_global_coupled_solve_max_delta_history.size() >= 2)
    {
        multi_step_stability_metric_v = 0.0;
        const std::size_t count =
            std::min<std::size_t>(3, shared_global_coupled_solve_max_delta_history.size());
        for (std::size_t offset = 1; offset < count; ++offset)
        {
            const auto current =
                shared_global_coupled_solve_max_delta_history[shared_global_coupled_solve_max_delta_history.size() - offset];
            const auto previous =
                shared_global_coupled_solve_max_delta_history[shared_global_coupled_solve_max_delta_history.size() - offset - 1];
            multi_step_stability_metric_v = std::max(
                multi_step_stability_metric_v, std::abs(current - previous));
        }
    }

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.surface_global_sheath_field_solve.v1\",\n";
    output << "  \"contract_id\": \"surface-pic-global-sheath-field-solve-v1\",\n";
    output << "  \"runtime_route\": \"" << runtimeRouteName(config.runtime_route) << "\",\n";
    output << "  \"surface_pic_runtime\": \"" << surfacePicRuntimeName(config.surface_pic_runtime_kind)
           << "\",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << latest_scalar(time_history) << ",\n";
    output << "  \"global_sheath_field_solve_active\": "
           << (active ? "true" : "false") << ",\n";
    output << "  \"global_sheath_field_solve_mode\": "
           << "\"" << jsonEscape(global_sheath_field_solve_state.solve_mode) << "\",\n";
    output << "  \"runtime_state_backed\": true,\n";
    output << "  \"shared_patch_node_count\": " << shared_patch_nodes.size() << ",\n";
    output << "  \"global_reference_potential_v\": " << global_reference_potential_v << ",\n";
    output << "  \"global_patch_potential_v\": " << global_patch_potential_v << ",\n";
    output << "  \"global_effective_sheath_length_m\": "
           << global_effective_sheath_length_m << ",\n";
    output << "  \"global_normal_electric_field_v_per_m\": "
           << global_normal_electric_field_v_per_m << ",\n";
    output << "  \"global_local_charge_density_c_per_m3\": "
           << global_local_charge_density_c_per_m3 << ",\n";
    output << "  \"shared_global_coupled_solve_iterations\": "
           << latest_scalar(shared_global_coupled_solve_iteration_history) << ",\n";
    output << "  \"shared_global_coupled_solve_converged\": "
           << (latest_scalar(shared_global_coupled_solve_converged_history) >= 0.5 ? "true"
                                                                                   : "false")
           << ",\n";
    output << "  \"shared_global_coupled_solve_max_delta_v\": "
           << latest_scalar(shared_global_coupled_solve_max_delta_history) << ",\n";
    output << "  \"global_reference_spread_v\": " << reference_spread_v << ",\n";
    output << "  \"global_patch_potential_spread_v\": " << patch_spread_v << ",\n";
    output << "  \"global_normal_electric_field_spread_v_per_m\": "
           << normal_field_spread_v_per_m << ",\n";
    output << "  \"global_local_charge_density_spread_c_per_m3\": "
           << local_charge_density_spread_c_per_m3 << ",\n";
    output << "  \"field_residual_v_per_m\": " << field_residual_v_per_m << ",\n";
    output << "  \"particle_field_coupled_residual_v\": "
           << particle_field_coupled_residual_v << ",\n";
    output << "  \"multi_step_stability_metric_v\": "
           << multi_step_stability_metric_v << ",\n";
    output << "  \"linear_residual_norm_v\": "
           << (global_sheath_field_solve_state.active
                   ? global_sheath_field_solve_state.linear_residual_norm_v
                   : 0.0)
           << ",\n";
    output << "  \"matrix_row_count\": "
           << (global_sheath_field_solve_state.active
                   ? global_sheath_field_solve_state.matrix_row_count
                   : 0.0)
           << ",\n";
    output << "  \"matrix_nonzeros\": "
           << (global_sheath_field_solve_state.active
                   ? global_sheath_field_solve_state.matrix_nonzeros
                   : 0.0)
           << ",\n";
    output << "  \"global_particle_charge_conservation_error_c\": "
           << global_particle_charge_conservation_error_c << ",\n";
    output << "  \"global_particle_flux_conservation_error_a\": "
           << global_particle_flux_conservation_error_a << ",\n";
    output << "  \"nodes\": [\n";
    if (use_owned_runtime_state)
    {
        const auto find_particle_node = [&global_particle_domain_state](std::size_t node_index)
            -> const GlobalParticleDomainNodeState* {
            for (const auto& node_state : global_particle_domain_state.nodes)
            {
                if (node_state.node_index == node_index)
                {
                    return &node_state;
                }
            }
            return nullptr;
        };

        for (std::size_t i = 0; i < global_sheath_field_solve_state.nodes.size(); ++i)
        {
            const auto& sheath_node = global_sheath_field_solve_state.nodes[i];
            const auto* particle_node = find_particle_node(sheath_node.node_index);
            output << "    {\"index\": " << sheath_node.node_index << ", \"name\": \""
                   << jsonEscape(sheath_node.node_name)
                   << "\", \"latest_patch_potential_v\": "
                   << sheath_node.patch_potential_v
                   << ", \"latest_shared_reference_potential_v\": "
                   << sheath_node.reference_potential_v
                   << ", \"latest_normal_electric_field_v_per_m\": "
                   << sheath_node.normal_electric_field_v_per_m
                   << ", \"latest_local_charge_density_c_per_m3\": "
                   << sheath_node.local_charge_density_c_per_m3
                   << ", \"global_node_charge_c\": "
                   << (particle_node != nullptr ? particle_node->charge_c
                                                : sheath_node.particle_charge_c)
                   << ", \"global_node_reference_shift_v\": "
                   << (particle_node != nullptr ? particle_node->reference_shift_v : 0.0)
                   << ", \"global_node_flux_a\": "
                   << (particle_node != nullptr ? particle_node->net_flux_a
                                                : sheath_node.particle_flux_a)
                   << "}";
            output << (i + 1 < global_sheath_field_solve_state.nodes.size() ? ",\n" : "\n");
        }
    }
    else
    {
        for (std::size_t i = 0; i < shared_patch_nodes.size(); ++i)
        {
            const auto node_index = shared_patch_nodes[i];
            output << "    {\"index\": " << node_index << ", \"name\": \""
                   << jsonEscape(circuit_model->nodeName(node_index))
                   << "\", \"latest_patch_potential_v\": "
                   << latestHistoryValue(node_potential_history, node_index)
                   << ", \"latest_shared_reference_potential_v\": "
                   << latestHistoryValue(node_shared_reference_potential_history, node_index)
                   << ", \"latest_normal_electric_field_v_per_m\": "
                   << latestHistoryValue(node_normal_electric_field_history, node_index)
                   << ", \"latest_local_charge_density_c_per_m3\": "
                   << latestHistoryValue(node_local_charge_density_history, node_index)
                   << ", \"global_node_charge_c\": "
                   << latestHistoryValue(node_distributed_particle_transport_charge_history,
                                         node_index)
                   << ", \"global_node_reference_shift_v\": "
                   << latestHistoryValue(node_distributed_particle_transport_reference_shift_history,
                                         node_index)
                   << ", \"global_node_flux_a\": "
                   << latestHistoryValue(node_distributed_particle_transport_net_flux_history,
                                         node_index)
                   << "}";
            output << (i + 1 < shared_patch_nodes.size() ? ",\n" : "\n");
        }
    }
    output << "  ]\n";
    output << "}\n";
    return true;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
