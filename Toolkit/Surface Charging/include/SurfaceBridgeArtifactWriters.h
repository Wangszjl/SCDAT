#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include <filesystem>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceGraphMonitorSnapshot
{
    std::size_t node_count = 0;
    std::size_t branch_count = 0;
    std::size_t max_node_degree = 0;
    double average_node_degree = 0.0;
    std::size_t strongest_field_node_index = 0;
    double strongest_field_node_abs_v_per_m = 0.0;
    std::size_t strongest_coupling_interface_index = 0;
    double strongest_coupling_interface_metric = 0.0;
    double max_abs_neighbor_potential_delta_v = 0.0;
    double max_abs_neighbor_field_contrast_v_per_m = 0.0;
    double graph_coupling_metric = 0.0;
};

struct SurfaceFieldMonitorSnapshot
{
    std::size_t strongest_field_node_index = 0;
    double strongest_field_node_abs_v_per_m = 0.0;
    double max_abs_node_reference_offset_v = 0.0;
    double max_abs_node_field_solver_reference_offset_v = 0.0;
    double max_field_solver_coupling_gain = 0.0;
    double reference_offset_envelope_v = 0.0;
};

struct SurfaceBenchmarkMonitorSnapshot
{
    std::string runtime_route;
    std::string benchmark_mode;
    std::string benchmark_source;
    std::string benchmark_execution_mode;
    std::string consistency_status;
    std::string consistency_authority;
    double patch_rmse_v = 0.0;
    double body_rmse_v = 0.0;
    double terminal_potential_v = 0.0;
    double time_to_equilibrium_ms = 0.0;
    std::size_t compared_patch_samples = 0;
    std::size_t compared_body_samples = 0;
};

bool writeSurfaceMonitorJson(const std::filesystem::path& json_path,
                             const SurfaceCircuitModel* circuit_model,
                             const SurfaceGraphMonitorSnapshot& graph_snapshot,
                             const SurfaceFieldMonitorSnapshot& field_snapshot,
                             const SurfaceBenchmarkMonitorSnapshot& benchmark_snapshot);
bool writeExternalFieldSolveResultTemplateJson(const std::filesystem::path& json_path,
                                               const SurfaceCircuitModel* circuit_model);
bool writeExternalFieldSolveResultJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const std::vector<std::vector<double>>& node_field_solver_reference_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_charge_history,
    const std::vector<std::vector<double>>& node_field_solver_capacitance_scale_history);
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
    const std::vector<std::vector<double>>& branch_current_history);
bool writeSurfaceVolumeProjectionJson(const std::filesystem::path& json_path,
                                      const SurfaceChargingConfig& config,
                                      const SurfaceCircuitModel* circuit_model,
                                      const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter);
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
    const std::vector<std::vector<double>>& node_poisson_residual_history);
bool writeVolumeMeshSkeletonJson(const std::filesystem::path& json_path,
                                 const SurfaceChargingConfig& config,
                                 const SurfaceCircuitModel* circuit_model,
                                 const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter);
bool writeExternalFieldBridgeManifestJson(const std::filesystem::path& json_path,
                                          const std::filesystem::path& csv_path,
                                          const SurfaceChargingConfig& config,
                                          const SurfaceCircuitModel* circuit_model,
                                          const SurfaceFieldSolverAdapter* field_solver_adapter,
                                          const SurfaceGraphCapacitanceMatrixProvider* matrix_provider);
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
    const std::vector<std::vector<double>>& node_external_volume_feedback_applied_history);
bool writeVolumetricSolverAdapterContractJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter,
    const SurfaceCircuitModel* circuit_model);
bool writeExternalVolumeSolveRequestJson(const std::filesystem::path& json_path,
                                         const std::filesystem::path& csv_path,
                                         const SurfaceChargingConfig& config,
                                         const SurfaceCircuitModel* circuit_model,
                                         const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter);
bool writeExternalVolumeSolveResultTemplateJson(
    const std::filesystem::path& json_path, const SurfaceCircuitModel* circuit_model,
    const SurfaceVolumetricSolverAdapter* volumetric_solver_adapter);
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
    const std::vector<std::vector<double>>& node_effective_sheath_length_history);
bool writeSurfaceGraphMatrixSnapshotCsv(
    const std::filesystem::path& csv_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider);
bool writeSurfaceGraphMatrixSnapshotJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider);
bool writeSurfaceFieldSolverAdapterContractReport(
    const std::filesystem::path& report_path, const SurfaceFieldSolverAdapter* field_solver_adapter,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider,
    const SurfaceCircuitModel* circuit_model);
bool writeSurfaceFieldSolverAdapterContractJson(
    const std::filesystem::path& json_path, const SurfaceFieldSolverAdapter* field_solver_adapter,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider,
    const SurfaceCircuitModel* circuit_model);
bool writeSurfaceBoundaryMappingJson(const std::filesystem::path& json_path,
                                     const SurfaceChargingConfig& config,
                                     const SurfaceCircuitModel* circuit_model,
                                     const SurfaceCircuitMappingState& mapping_state);
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
    const std::vector<std::vector<double>>& branch_mutual_capacitance_history);
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
    const std::vector<std::vector<double>>& node_shared_sheath_charge_history);
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
    const std::vector<std::vector<double>>& node_ion_calibration_factor_history);
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
    const std::vector<std::vector<double>>& exchange_flux_matrix_a);
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
    const std::vector<std::vector<double>>& exchange_flux_matrix_a);
bool writeGlobalParticleRepositoryJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const std::vector<double>& time_history,
    const GlobalParticleRepositoryState& global_particle_repository_state);
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
    const std::vector<double>& shared_global_coupled_solve_max_delta_history);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
