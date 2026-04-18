#include "SurfaceObserverArtifactExport.h"

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

bool exportSurfaceObserverArtifacts(
    const std::filesystem::path& csv_path,
    const std::function<bool(const std::filesystem::path&)>& write_monitor_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_shared_runtime_observer_artifact,
    std::string& last_error_message)
{
    auto monitor_json_path = csv_path;
    monitor_json_path.replace_extension(".monitor.json");
    if (!write_monitor_artifact(monitor_json_path))
    {
        last_error_message = "Failed to write surface monitor json to: " +
                             monitor_json_path.string();
        return false;
    }

    auto shared_runtime_observer_path = csv_path;
    shared_runtime_observer_path.replace_extension(".shared_runtime_observer.json");
    if (!write_shared_runtime_observer_artifact(shared_runtime_observer_path))
    {
        last_error_message =
            "Failed to write shared surface runtime observer json to: " +
            shared_runtime_observer_path.string();
        return false;
    }

    return true;
}

bool exportSharedRuntimeArtifactSuite(
    const std::filesystem::path& csv_path, bool export_shared_particle_transport_domain,
    const std::function<bool(const std::filesystem::path&)>&
        write_shared_runtime_consistency_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_shared_particle_transport_domain_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_global_particle_domain_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_global_particle_repository_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_global_sheath_field_solve_artifact,
    std::string& last_error_message)
{
    auto shared_runtime_consistency_path = csv_path;
    shared_runtime_consistency_path.replace_extension(".shared_runtime_consistency.json");
    if (!write_shared_runtime_consistency_artifact(shared_runtime_consistency_path))
    {
        last_error_message =
            "Failed to write shared surface runtime consistency json to: " +
            shared_runtime_consistency_path.string();
        return false;
    }

    if (export_shared_particle_transport_domain)
    {
        auto shared_particle_transport_domain_path = csv_path;
        shared_particle_transport_domain_path.replace_extension(
            ".shared_particle_transport_domain.json");
        if (!write_shared_particle_transport_domain_artifact(
                shared_particle_transport_domain_path))
        {
            last_error_message =
                "Failed to write shared particle transport domain json to: " +
                shared_particle_transport_domain_path.string();
            return false;
        }
    }

    auto global_particle_domain_path = csv_path;
    global_particle_domain_path.replace_extension(".global_particle_domain.json");
    if (!write_global_particle_domain_artifact(global_particle_domain_path))
    {
        last_error_message =
            "Failed to write global particle domain json to: " +
            global_particle_domain_path.string();
        return false;
    }

    auto global_particle_repository_path = csv_path;
    global_particle_repository_path.replace_extension(".global_particle_repository.json");
    if (!write_global_particle_repository_artifact(global_particle_repository_path))
    {
        last_error_message =
            "Failed to write global particle repository json to: " +
            global_particle_repository_path.string();
        return false;
    }

    auto global_sheath_field_solve_path = csv_path;
    global_sheath_field_solve_path.replace_extension(".global_sheath_field_solve.json");
    if (!write_global_sheath_field_solve_artifact(global_sheath_field_solve_path))
    {
        last_error_message =
            "Failed to write global sheath-field solve json to: " +
            global_sheath_field_solve_path.string();
        return false;
    }

    return true;
}

bool exportSurfaceBridgeArtifactSuite(
    const std::filesystem::path& csv_path,
    const std::function<bool(const std::filesystem::path&)>& write_graph_report_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_graph_matrix_csv_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_graph_matrix_json_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_field_adapter_report_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_field_adapter_json_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_boundary_mapping_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_field_request_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_field_result_template_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_field_result_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_stub_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_volume_mesh_stub_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_field_bridge_manifest_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_surface_volume_projection_artifact,
    std::string& last_error_message)
{
    auto graph_report_path = csv_path;
    graph_report_path.replace_extension(".graph.txt");
    if (!write_graph_report_artifact(graph_report_path))
    {
        last_error_message =
            "Failed to write surface graph sidecar report to: " +
            graph_report_path.string();
        return false;
    }

    auto graph_matrix_path = csv_path;
    graph_matrix_path.replace_extension(".graph_matrix.csv");
    if (!write_graph_matrix_csv_artifact(graph_matrix_path))
    {
        last_error_message =
            "Failed to write surface graph matrix snapshot to: " +
            graph_matrix_path.string();
        return false;
    }

    auto graph_matrix_json_path = csv_path;
    graph_matrix_json_path.replace_extension(".graph_matrix.json");
    if (!write_graph_matrix_json_artifact(graph_matrix_json_path))
    {
        last_error_message =
            "Failed to write surface graph matrix snapshot json to: " +
            graph_matrix_json_path.string();
        return false;
    }

    auto field_adapter_path = csv_path;
    field_adapter_path.replace_extension(".field_adapter.txt");
    if (!write_field_adapter_report_artifact(field_adapter_path))
    {
        last_error_message =
            "Failed to write surface field-solver adapter report to: " +
            field_adapter_path.string();
        return false;
    }

    auto field_adapter_json_path = csv_path;
    field_adapter_json_path.replace_extension(".field_adapter.json");
    if (!write_field_adapter_json_artifact(field_adapter_json_path))
    {
        last_error_message =
            "Failed to write surface field-solver adapter json to: " +
            field_adapter_json_path.string();
        return false;
    }

    auto boundary_mapping_path = csv_path;
    boundary_mapping_path.replace_extension(".boundary_mapping.json");
    if (!write_boundary_mapping_artifact(boundary_mapping_path))
    {
        last_error_message =
            "Failed to write surface boundary mapping json to: " +
            boundary_mapping_path.string();
        return false;
    }

    auto field_request_path = csv_path;
    field_request_path.replace_extension(".field_request.json");
    if (!write_field_request_artifact(field_request_path))
    {
        last_error_message =
            "Failed to write external field solve request json to: " +
            field_request_path.string();
        return false;
    }

    auto field_result_template_path = csv_path;
    field_result_template_path.replace_extension(".field_result_template.json");
    if (!write_field_result_template_artifact(field_result_template_path))
    {
        last_error_message =
            "Failed to write external field solve result template json to: " +
            field_result_template_path.string();
        return false;
    }

    auto field_result_path = csv_path;
    field_result_path.replace_extension(".field_result.json");
    if (!write_field_result_artifact(field_result_path))
    {
        last_error_message =
            "Failed to write external field solve result json to: " +
            field_result_path.string();
        return false;
    }

    auto volume_stub_path = csv_path;
    volume_stub_path.replace_extension(".volume_stub.json");
    if (!write_volume_stub_artifact(volume_stub_path))
    {
        last_error_message =
            "Failed to write external volume stub json to: " +
            volume_stub_path.string();
        return false;
    }

    auto volume_mesh_stub_path = csv_path;
    volume_mesh_stub_path.replace_extension(".volume_mesh_stub.json");
    if (!write_volume_mesh_stub_artifact(volume_mesh_stub_path))
    {
        last_error_message =
            "Failed to write volume mesh skeleton json to: " +
            volume_mesh_stub_path.string();
        return false;
    }

    auto field_bridge_manifest_path = csv_path;
    field_bridge_manifest_path.replace_extension(".field_bridge_manifest.json");
    if (!write_field_bridge_manifest_artifact(field_bridge_manifest_path))
    {
        last_error_message =
            "Failed to write external field bridge manifest json to: " +
            field_bridge_manifest_path.string();
        return false;
    }

    auto surface_volume_projection_path = csv_path;
    surface_volume_projection_path.replace_extension(".surface_volume_projection.json");
    if (!write_surface_volume_projection_artifact(surface_volume_projection_path))
    {
        last_error_message =
            "Failed to write surface-volume projection json to: " +
            surface_volume_projection_path.string();
        return false;
    }

    return true;
}

bool exportVolumetricArtifactSuite(
    const std::filesystem::path& csv_path,
    const std::function<bool(const std::filesystem::path&)>&
        write_volumetric_adapter_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_request_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_volume_result_template_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_result_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_history_artifact,
    std::string& last_error_message)
{
    auto volumetric_adapter_path = csv_path;
    volumetric_adapter_path.replace_extension(".volumetric_adapter.json");
    if (!write_volumetric_adapter_artifact(volumetric_adapter_path))
    {
        last_error_message =
            "Failed to write volumetric adapter contract json to: " +
            volumetric_adapter_path.string();
        return false;
    }

    auto volume_request_path = csv_path;
    volume_request_path.replace_extension(".volume_request.json");
    if (!write_volume_request_artifact(volume_request_path))
    {
        last_error_message =
            "Failed to write external volume solve request json to: " +
            volume_request_path.string();
        return false;
    }

    auto volume_result_template_path = csv_path;
    volume_result_template_path.replace_extension(".volume_result_template.json");
    if (!write_volume_result_template_artifact(volume_result_template_path))
    {
        last_error_message =
            "Failed to write external volume solve result template json to: " +
            volume_result_template_path.string();
        return false;
    }

    auto volume_result_path = csv_path;
    volume_result_path.replace_extension(".volume_result.json");
    if (!write_volume_result_artifact(volume_result_path))
    {
        last_error_message =
            "Failed to write external volume solve result json to: " +
            volume_result_path.string();
        return false;
    }

    auto volume_history_path = csv_path;
    volume_history_path.replace_extension(".volume_history.json");
    if (!write_volume_history_artifact(volume_history_path))
    {
        last_error_message =
            "Failed to write volume history json to: " +
            volume_history_path.string();
        return false;
    }

    return true;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
