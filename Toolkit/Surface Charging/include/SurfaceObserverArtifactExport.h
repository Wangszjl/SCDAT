#pragma once

#include <filesystem>
#include <functional>
#include <string>

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
    std::string& last_error_message);

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
    std::string& last_error_message);

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
    std::string& last_error_message);

bool exportVolumetricArtifactSuite(
    const std::filesystem::path& csv_path,
    const std::function<bool(const std::filesystem::path&)>&
        write_volumetric_adapter_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_request_artifact,
    const std::function<bool(const std::filesystem::path&)>&
        write_volume_result_template_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_result_artifact,
    const std::function<bool(const std::filesystem::path&)>& write_volume_history_artifact,
    std::string& last_error_message);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
