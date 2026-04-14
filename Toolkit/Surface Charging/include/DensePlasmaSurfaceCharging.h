#pragma once

#include "ChargeAccumulationModel.h"
#include "PicMccSurfaceCurrentSampler.h"
#include "ReferenceCurrentBalanceModel.h"

#include "../../Plasma Analysis/include/FluidAlgorithmConfig.h"

#include "../../Tools/Coupling/include/BenchmarkContracts.h"
#include "../../Tools/Coupling/include/SurfaceCircuitCoupling.h"
#include "../../Tools/FieldSolver/include/SurfaceFieldVolumeContracts.h"
#include "../../Tools/Material/include/MaterialDatabase.h"
#include "../../Tools/Material/include/SurfaceInteraction.h"
#include "../../Tools/Output/include/ResultExporter.h"
#include "../../Tools/Particle/include/SurfaceVolumeDistributionContracts.h"
#include "../../Tools/Particle/include/ParticleSource.h"

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceRuntimePlan;

enum class SurfaceChargingRegime
{
    Generic,
    LeoFlowingPlasma,
    GeoKineticPicLike,
    ThrusterPlume
};

enum class SurfaceCurrentAlgorithmMode
{
    UnifiedKernelAligned,
    LegacyRefCompatible
};

enum class SurfaceCapacitanceModelKind
{
    Lumped
};

enum class PotentialReferenceModelKind
{
    Lumped
};

enum class SurfaceBenchmarkMode
{
    UnifiedKernelAligned,
    LegacyRefCompatible
};

enum class SurfaceRuntimeRoute
{
    SCDATUnified,
    SurfacePic,
    SurfacePicHybrid,
    LegacyBenchmark
};

enum class SurfacePicStrategy
{
    SurfacePicDirect,
    SurfacePicCalibrated,
    SurfacePicHybridReference
};

enum class SurfaceLegacyInputAdapterKind
{
    None,
    CTextReferenceDeck,
    MatlabReferenceDeck
};

enum class SurfacePicRuntimeKind
{
    LocalWindowSampler,
    GraphCoupledSharedSurface
};

enum class SurfaceInstrumentSetKind
{
    MetadataOnly,
    SurfacePicObserverSet
};

enum class SurfaceBenchmarkSource
{
    None,
    CGeo,
    CLeoRam,
    CLeoWake,
    MatlabGeo,
    MatlabLeo
};

enum class LegacyBenchmarkExecutionMode
{
    ReplayFromReference,
    ExecuteLegacyAlgorithm
};

enum class VolumeLinearSolverPolicy
{
    Auto,
    DenseOnly,
    IterativeOnly
};

struct SurfaceCurrents
{
    double electron_current_a_per_m2 = 0.0;
    double ion_current_a_per_m2 = 0.0;
    double secondary_emission_a_per_m2 = 0.0;
    double ion_secondary_emission_a_per_m2 = 0.0;
    double backscatter_emission_a_per_m2 = 0.0;
    double photo_emission_a_per_m2 = 0.0;
    double thermionic_emission_a_per_m2 = 0.0;
    double field_emission_a_per_m2 = 0.0;
    double conduction_current_a_per_m2 = 0.0;
    double ram_ion_current_a_per_m2 = 0.0;
    double current_derivative_a_per_m2_per_v = 0.0;
    double total_current_a_per_m2 = 0.0;
};

struct SurfaceDistributionMoments
{
    PlasmaAnalysis::PlasmaDistributionModelKind model =
        PlasmaAnalysis::PlasmaDistributionModelKind::MultiPopulationHybrid;
    double electron_collection_scale = 1.0;
    double ion_collection_scale = 1.0;
    double electron_characteristic_energy_ev = 1.0;
    double ion_characteristic_energy_ev = 1.0;
    double photo_incidence_scale = 1.0;
    double local_reference_shift_v = 0.0;
    bool valid = false;
};

struct SurfaceMaterialInteractionMoments
{
    double electron_absorption_scale = 1.0;
    double secondary_emission_scale = 1.0;
    double ion_secondary_emission_scale = 1.0;
    double backscatter_scale = 1.0;
    double photo_emission_scale = 1.0;
    double conductivity_scale = 1.0;
    double capacitance_scale = 1.0;
    bool valid = false;
};

struct EmissionModelParameters
{
    double surface_temperature_k = 300.0;
    double enhancement_factor = 1.0;
    double photon_flux_m2_s = 0.0;
};

struct SurfaceCircuitNodeConfig
{
    std::string name;
    double area_m2 = 0.0;
    bool is_patch = true;
    double initial_potential_v = 0.0;
    double capacitance_f = 0.0;
    bool fixed_potential = false;
    double fixed_value_v = 0.0;
};

struct SurfaceCircuitBranchConfig
{
    std::size_t from_node = 0;
    std::size_t to_node = 0;
    double conductance_s = 0.0;
    double resistance_ohm = 0.0;
    double bias_v = 0.0;
};

struct SurfacePatchPhysicsConfig
{
    std::string node_name;
    std::size_t node_index = 0;
    bool match_by_index = false;
    bool override_material = false;
    Material::MaterialProperty material{2, Mesh::MaterialType::DIELECTRIC, "kapton"};
    bool override_see_model = false;
    SecondaryElectronEmissionModel reference_see_model =
        SecondaryElectronEmissionModel::Whipple;
    bool override_electron_collection_model = false;
    ElectronCollectionModelKind electron_collection_model =
        ElectronCollectionModelKind::OmlLike;
    bool override_patch_incidence_angle = false;
    double patch_incidence_angle_deg = 0.0;
    bool override_patch_flow_angle = false;
    double patch_flow_angle_deg = 0.0;
    bool override_photoelectron_temperature = false;
    double photoelectron_temperature_ev = 2.0;
    bool override_patch_photo_current_density = false;
    double patch_photo_current_density_a_per_m2 = 0.0;
    bool override_electron_collection_coefficient = false;
    double electron_collection_coefficient = 1.0;
    bool override_ion_collection_coefficient = false;
    double ion_collection_coefficient = 1.0;
    bool override_plasma = false;
    PlasmaAnalysis::PlasmaParameters plasma;
    bool override_electron_spectrum = false;
    Particle::ResolvedSpectrum electron_spectrum;
    bool has_electron_spectrum = false;
    bool override_ion_spectrum = false;
    Particle::ResolvedSpectrum ion_spectrum;
    bool has_ion_spectrum = false;
    bool override_emission = false;
    EmissionModelParameters emission;
};

struct DefaultSurfacePhysicsConfig
{
    Material::MaterialProperty material{2, Mesh::MaterialType::DIELECTRIC, "kapton"};
    SecondaryElectronEmissionModel reference_see_model =
        SecondaryElectronEmissionModel::Whipple;
    ElectronCollectionModelKind electron_collection_model =
        ElectronCollectionModelKind::OmlLike;
    double patch_incidence_angle_deg = 0.0;
    double patch_flow_angle_deg = 0.0;
    double photoelectron_temperature_ev = 2.0;
    double patch_photo_current_density_a_per_m2 = 0.0;
    double electron_collection_coefficient = 1.0;
    double ion_collection_coefficient = 1.0;
    PlasmaAnalysis::PlasmaParameters plasma;
    Particle::ResolvedSpectrum electron_spectrum;
    Particle::ResolvedSpectrum ion_spectrum;
    bool has_electron_spectrum = false;
    bool has_ion_spectrum = false;
    EmissionModelParameters emission;
};

struct StructureBodyConfig
{
    std::string id;
    double area_m2 = 0.0;
    bool floating = false;
    double initial_potential_v = 0.0;
    double capacitance_f = 0.0;
    double leakage_conductance_s = 0.0;
    bool fixed_potential = false;
    double fixed_value_v = 0.0;
};

struct SurfacePatchConfig
{
    std::string id;
    std::string body_id;
    double area_m2 = 0.0;
    double initial_potential_v = 0.0;
    double capacitance_f = 0.0;
    std::optional<Material::MaterialProperty> material;
    std::optional<SecondaryElectronEmissionModel> reference_see_model;
    std::optional<ElectronCollectionModelKind> electron_collection_model;
    std::optional<double> patch_incidence_angle_deg;
    std::optional<double> patch_flow_angle_deg;
    std::optional<double> photoelectron_temperature_ev;
    std::optional<double> patch_photo_current_density_a_per_m2;
    std::optional<double> electron_collection_coefficient;
    std::optional<double> ion_collection_coefficient;
    std::optional<PlasmaAnalysis::PlasmaParameters> plasma;
    std::optional<Particle::ResolvedSpectrum> electron_spectrum;
    std::optional<Particle::ResolvedSpectrum> ion_spectrum;
    std::optional<bool> has_electron_spectrum;
    std::optional<bool> has_ion_spectrum;
    std::optional<EmissionModelParameters> emission;
};

struct PatchInterfaceConfig
{
    std::string id;
    std::string from_id;
    std::string to_id;
    double conductance_s = 0.0;
    double resistance_ohm = 0.0;
    double bias_v = 0.0;
    bool dynamic_conductance = false;
};

struct BodyBoundaryGroup
{
    std::string id;
    std::string body_id;
    std::string external_group_name;
    std::vector<int> boundary_face_ids;
};

struct PatchBoundaryGroup
{
    std::string id;
    std::string patch_id;
    std::string external_group_name;
    std::vector<int> boundary_face_ids;
};

struct SurfaceBoundaryMapping
{
    std::string node_id;
    std::string boundary_group_id;
    bool required = true;
};

using ExternalFieldSolveNodeResult = FieldSolver::ExternalFieldSolveNodeResult;

struct ExternalFieldSolveRequest
{
    std::string runtime_route;
    std::string circuit_model;
    std::string matrix_family;
    std::string solver_adapter_hint;
    std::vector<SurfaceBoundaryMapping> mappings;
};

struct ExternalFieldSolveResult
{
    std::string schema_version = "scdat.external_field_result.v1";
    std::vector<ExternalFieldSolveNodeResult> nodes;
};

struct VolumeMeshCellStub : public Particle::VolumeCellDistributionStubContract
{
};

struct VolumeMeshBoundaryFaceStub : public Particle::VolumeBoundaryFaceStubContract
{
};

struct SurfaceVolumeProjectionEntry : public Particle::SurfaceVolumeProjectionEntryContract
{
};

using ExternalVolumeSolveCellResult = FieldSolver::ExternalVolumeSolveCellResult;

struct ExternalVolumeSolveRequest
{
    std::string schema_version = "scdat.external_volume_request.v1";
    std::string runtime_route;
    std::string volumetric_adapter;
    std::string projection_family;
};

struct ExternalVolumeSolveResult
{
    std::string schema_version = "scdat.external_volume_result.v1";
    std::vector<ExternalVolumeSolveCellResult> cells;
};

struct SurfaceChargingConfig
{
    Coupling::Contracts::SolverConfig solver_config{};
    unsigned int seed = 20260408u;
    std::string sampling_policy = "deterministic";
    SurfaceRuntimeRoute runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    SurfacePicStrategy surface_pic_strategy = SurfacePicStrategy::SurfacePicCalibrated;
    SurfaceLegacyInputAdapterKind legacy_input_adapter_kind =
        SurfaceLegacyInputAdapterKind::None;
    SurfacePicRuntimeKind surface_pic_runtime_kind =
        SurfacePicRuntimeKind::LocalWindowSampler;
    SurfaceInstrumentSetKind surface_instrument_set_kind =
        SurfaceInstrumentSetKind::MetadataOnly;
    SurfaceBenchmarkSource benchmark_source = SurfaceBenchmarkSource::None;
    std::string reference_family = "native";
    std::string reference_case_id;
    std::string reference_matrix_case_id;
    LegacyBenchmarkExecutionMode legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    DefaultSurfacePhysicsConfig default_surface_physics;
    double surface_area_m2 = 1.0e-4;
    bool floating = true;
    bool derive_capacitance_from_material = true;
    bool derive_sheath_length_from_plasma = true;
    double capacitance_per_area_f_per_m2 = 8.0e-10;
    double dielectric_thickness_m = 1.25e-4;
    double sheath_length_m = 1.0e-3;
    double minimum_sheath_length_m = 5.0e-5;
    double maximum_sheath_length_m = 10.0;
    double radiation_conductivity_coefficient = 0.0;
    double radiation_conductivity_exponent = 1.0;
    SurfaceChargingRegime regime = SurfaceChargingRegime::Generic;
    double bulk_flow_velocity_m_per_s = 0.0;
    double flow_alignment_cosine = 1.0;
    double electron_flow_coupling = 0.0;
    ElectronCollectionModelKind electron_collection_model =
        ElectronCollectionModelKind::OmlLike;
    double electron_collection_coefficient = 1.0;
    double ion_collection_coefficient = 1.0;
    double ion_directed_velocity_m_per_s = 0.0;
    SurfaceCurrentAlgorithmMode current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    SurfaceCapacitanceModelKind capacitance_model_kind =
        SurfaceCapacitanceModelKind::Lumped;
    PotentialReferenceModelKind potential_reference_model_kind =
        PotentialReferenceModelKind::Lumped;
    SurfaceBenchmarkMode benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    bool enable_pic_calibration = false;
    std::size_t pic_calibration_samples = 4096;
    std::size_t pic_recalibration_interval_steps = 0;
    double pic_recalibration_trigger_v = 25.0;
    double electron_pic_calibration_min = 0.5;
    double electron_pic_calibration_max = 2.0;
    double ion_pic_calibration_min = 0.5;
    double ion_pic_calibration_max = 2.0;
    bool enable_live_pic_window = false;
    bool enable_live_pic_mcc = true;
    std::string live_pic_collision_cross_section_set_id = "surface_pic_v1";
    std::size_t live_pic_window_steps = 12;
    std::size_t live_pic_window_layers = 8;
    std::size_t live_pic_particles_per_element = 2;
    double live_pic_probe_delta_v = 1.0;
    double live_pic_reference_potential_v = 0.0;
    bool use_reference_current_balance = true;
    bool enable_body_patch_circuit = true;
    bool body_floating = false;
    double body_initial_potential_v = 0.0;
    double body_capacitance_f = 1.0e-10;
    double body_leakage_conductance_s = 0.0;
    double patch_body_resistance_ohm = 0.0;
    double patch_incidence_angle_deg = 0.0;
    double patch_flow_angle_deg = 0.0;
    double photoelectron_temperature_ev = 2.0;
    double body_photo_current_density_a_per_m2 = 0.0;
    double patch_photo_current_density_a_per_m2 = 0.0;
    bool enable_secondary_electron = true;
    bool enable_backscatter = true;
    bool enable_photoelectron = true;
    double max_delta_potential_v_per_step = 250.0;
    SecondaryElectronEmissionModel reference_see_model =
        SecondaryElectronEmissionModel::Whipple;
    double max_abs_potential_v = 5.0e4;
    double max_abs_current_density_a_per_m2 = 5.0e3;
    std::size_t internal_substeps = 16;
    PlasmaAnalysis::PlasmaParameters plasma;
    PlasmaAnalysis::PlasmaDistributionModelKind distribution_model =
        PlasmaAnalysis::PlasmaDistributionModelKind::MultiPopulationHybrid;
    Particle::ResolvedSpectrum electron_spectrum;
    Particle::ResolvedSpectrum ion_spectrum;
    bool has_electron_spectrum = false;
    bool has_ion_spectrum = false;
    Material::MaterialProperty material{2, Mesh::MaterialType::DIELECTRIC, "kapton"};
    std::filesystem::path material_library_path;
    std::string imported_material_name;
    EmissionModelParameters emission;
    std::vector<StructureBodyConfig> bodies;
    std::vector<SurfacePatchConfig> patches;
    std::vector<PatchInterfaceConfig> interfaces;
    std::vector<BodyBoundaryGroup> body_boundary_groups;
    std::vector<PatchBoundaryGroup> patch_boundary_groups;
    std::vector<SurfaceBoundaryMapping> boundary_mappings;
    bool enable_external_field_solver_bridge = false;
    bool enable_external_volume_solver_bridge = false;
    std::filesystem::path external_field_solver_request_path;
    std::filesystem::path external_field_solver_result_path;
    std::filesystem::path external_volume_solver_request_path;
    std::filesystem::path external_volume_solver_result_path;
    std::filesystem::path external_volume_mesh_path;
    std::filesystem::path external_surface_volume_projection_path;
    VolumeLinearSolverPolicy volume_linear_solver_policy = VolumeLinearSolverPolicy::Auto;
    std::size_t volume_linear_max_iterations = 128;
    double volume_linear_tolerance_scale = 0.25;
    double volume_linear_relaxation = 0.85;
    std::size_t field_volume_outer_iterations = 3;
    double field_volume_outer_tolerance = 1.0e-3;
    double field_volume_outer_relaxation = 0.65;
    std::size_t volume_self_consistent_iterations = 4;
    double volume_self_consistent_tolerance_v = 1.0e-3;
    double volume_charge_relaxation = 0.45;
    double volume_potential_relaxation = 0.50;
    std::vector<SurfaceCircuitNodeConfig> surface_nodes;
    std::vector<SurfaceCircuitBranchConfig> surface_branches;
    std::vector<SurfacePatchPhysicsConfig> patch_physics_overrides;
};

struct SurfaceChargingStatus
{
    double time_s = 0.0;
    SurfaceCurrents currents;
    ChargeAccumulationState state;
    double equilibrium_error = 0.0;
    double body_potential_v = 0.0;
    double patch_potential_v = 0.0;
    double circuit_branch_current_a = 0.0;
    bool pic_recalibrated = false;
    bool equilibrium_reached = false;
    std::size_t steps_completed = 0;
    double transition_conductivity_scale = 1.0;
    double transition_source_flux_scale = 1.0;
    double transition_simulation_param_scale = 1.0;
    std::size_t transition_runtime_internal_substeps = 0;
    double transition_runtime_max_delta_potential_v_per_step = 0.0;
    double transition_runtime_solver_relaxation_factor = 0.0;
};

struct SurfaceModelRuntimeState
{
    double body_potential_v = 0.0;
    double patch_potential_v = 0.0;
    double surface_area_m2 = 0.0;
    double effective_conductivity_s_per_m = 0.0;
    double effective_sheath_length_m = 1.0e-3;
    double reference_plasma_potential_v = 0.0;
    double propagated_reference_potential_v = 0.0;
    double field_solver_reference_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double graph_capacitance_diagonal_f = 0.0;
    double graph_capacitance_row_sum_f = 0.0;
    double field_solver_coupling_gain = 0.0;
    double field_solver_capacitance_scale = 1.0;
    double pseudo_volume_m3 = 0.0;
    double volume_projection_weight_sum = 1.0;
    double volume_mesh_coupling_gain = 0.0;
    double volume_potential_v = 0.0;
    double deposited_charge_c = 0.0;
    double poisson_residual_v_m = 0.0;
    double volume_solver_mode_id = 0.0;
    double volume_solver_iterations = 0.0;
    double volume_solver_linear_iterations = 0.0;
    double volume_solver_converged = 0.0;
    double volume_solver_residual_norm = 0.0;
    double volume_solver_max_delta_v = 0.0;
    double volume_solver_matrix_nnz = 0.0;
    double volume_solver_cell_count = 0.0;
    double field_volume_coupling_iterations = 0.0;
    double field_volume_coupling_converged = 0.0;
    double field_volume_coupling_max_delta = 0.0;
    double field_volume_coupling_relaxation_used = 0.0;
    double external_volume_feedback_blend_factor = 0.0;
    double external_volume_feedback_mismatch_metric = 0.0;
    double external_volume_feedback_applied = 0.0;
    double electron_calibration_factor = 1.0;
    double ion_calibration_factor = 1.0;
    double shared_patch_potential_v = 0.0;
    double shared_patch_area_m2 = 0.0;
    double shared_reference_potential_v = 0.0;
    double shared_effective_sheath_length_m = 1.0e-3;
    double shared_particle_transport_charge_c = 0.0;
    double shared_particle_transport_reference_shift_v = 0.0;
    bool shared_particle_transport_domain_active = false;
    bool global_particle_domain_active = false;
    double global_particle_domain_charge_c = 0.0;
    double global_particle_domain_charge_conservation_error_c = 0.0;
    double global_particle_domain_flux_conservation_error_a = 0.0;
    double global_particle_domain_node_charge_c = 0.0;
    double global_particle_domain_node_flux_a = 0.0;
    bool global_sheath_field_solve_active = false;
    double global_sheath_field_reference_potential_v = 0.0;
    double global_sheath_field_residual_v_per_m = 0.0;
    double global_particle_field_coupled_residual_v = 0.0;
    double global_sheath_field_multi_step_stability_metric_v = 0.0;
    double distributed_particle_transport_charge_c = 0.0;
    double distributed_particle_transport_reference_shift_v = 0.0;
    double distributed_particle_transport_net_flux_a = 0.0;
    bool distributed_particle_transport_active = false;
    bool shared_runtime_enabled = false;
    std::size_t node_index = 0;
    std::string node_name;
    const PicMccCurrentSample* live_pic_sample = nullptr;
};

struct SurfaceDidvState
{
    std::size_t node_index = 0;
    std::string node_name;
    double node_area_m2 = 0.0;
    double plasma_current_a = 0.0;
    double plasma_didv_a_per_v = 0.0;
    double conduction_current_a = 0.0;
    double conduction_didv_a_per_v = 0.0;
    double total_current_a = 0.0;
    double total_didv_a_per_v = 0.0;
    bool valid = false;
};

struct SurfaceKernelSnapshot
{
    std::string source_family = "spis_equivalent_surface_kernel_v1";
    SurfaceDistributionMoments distribution{};
    SurfaceMaterialInteractionMoments material_interaction{};
    SurfaceCurrents currents{};
    SurfaceDidvState didv{};
    PicMccCurrentSample live_pic_sample{};
    bool valid = false;
};

struct SurfaceCircuitReducedNodeGroup
{
    std::size_t reduced_node_index = 0;
    std::string group_id;
    std::vector<std::size_t> member_node_indices;
    double total_area_m2 = 0.0;
};

struct SurfaceCircuitMappingState
{
    std::vector<std::size_t> surface_to_circuit_node_index;
    std::vector<std::vector<std::size_t>> circuit_to_surface_node_indices;
    std::vector<std::size_t> circuit_to_reduced_node_index;
    std::vector<SurfaceCircuitReducedNodeGroup> reduced_node_groups;
    bool valid = false;
};

struct GlobalParticleDomainNodeState
{
    std::size_t node_index = 0;
    std::string node_name;
    double area_m2 = 0.0;
    double patch_potential_v = 0.0;
    double shared_reference_potential_v = 0.0;
    double charge_c = 0.0;
    double reference_shift_v = 0.0;
    double net_flux_a = 0.0;
};

struct GlobalParticleDomainEdgeState
{
    std::size_t from_node_index = 0;
    std::size_t to_node_index = 0;
    std::string from_node_name;
    std::string to_node_name;
    double conductance_s = 0.0;
    double exchange_flux_a = 0.0;
    double stored_charge_c = 0.0;
    double target_charge_c = 0.0;
    double operator_drive_charge_c = 0.0;
};

struct GlobalParticleDomainState
{
    bool active = false;
    std::string bookkeeping_mode = "uninitialized";
    double total_charge_c = 0.0;
    double total_reference_shift_v = 0.0;
    double total_node_charge_c = 0.0;
    double total_node_flux_a = 0.0;
    double charge_conservation_error_c = 0.0;
    double flux_conservation_error_a = 0.0;
    double edge_charge_total_abs_c = 0.0;
    double edge_target_charge_total_abs_c = 0.0;
    double edge_operator_drive_total_abs_c = 0.0;
    double edge_conductance_total_s = 0.0;
    double edge_graph_operator_iterations = 0.0;
    bool edge_graph_operator_converged = false;
    double edge_graph_operator_max_balance_residual_c = 0.0;
    double edge_graph_operator_branch_graph_edge_count = 0.0;
    double edge_graph_operator_branch_graph_pair_count = 0.0;
    double edge_graph_operator_effective_pair_count = 0.0;
    double edge_graph_operator_total_pair_weight_f = 0.0;
    double edge_graph_operator_total_conductance_weight_f = 0.0;
    double edge_graph_operator_min_node_preconditioner = 1.0;
    double edge_graph_operator_max_node_preconditioner = 1.0;
    std::vector<GlobalParticleDomainNodeState> nodes;
    std::vector<GlobalParticleDomainEdgeState> edges;
};

struct GlobalParticleRepositoryNodeState
{
    std::size_t node_index = 0;
    std::string node_name;
    double area_m2 = 0.0;
    double patch_potential_v = 0.0;
    double shared_reference_potential_v = 0.0;
    double reference_shift_v = 0.0;
    double reservoir_charge_c = 0.0;
    double target_reservoir_charge_c = 0.0;
    double migration_delta_charge_c = 0.0;
    double edge_feedback_charge_c = 0.0;
    double conservation_correction_charge_c = 0.0;
    double net_flux_a = 0.0;
};

struct GlobalParticleRepositoryEdgeState
{
    std::size_t from_node_index = 0;
    std::size_t to_node_index = 0;
    std::string from_node_name;
    std::string to_node_name;
    double conductance_s = 0.0;
    double migration_flux_a = 0.0;
    double migration_charge_c = 0.0;
    double stored_charge_c = 0.0;
    double target_charge_c = 0.0;
    double operator_drive_charge_c = 0.0;
};

struct GlobalParticleRepositoryState
{
    bool active = false;
    std::string bookkeeping_mode = "uninitialized";
    std::string lifecycle_mode = "uninitialized";
    double total_reservoir_charge_c = 0.0;
    double total_target_reservoir_charge_c = 0.0;
    double total_migration_delta_abs_charge_c = 0.0;
    double total_edge_feedback_abs_charge_c = 0.0;
    double total_conservation_correction_abs_charge_c = 0.0;
    double total_migration_edge_abs_charge_c = 0.0;
    double charge_conservation_error_c = 0.0;
    double migration_charge_conservation_error_c = 0.0;
    std::vector<GlobalParticleRepositoryNodeState> nodes;
    std::vector<GlobalParticleRepositoryEdgeState> edges;
};

struct GlobalSheathFieldNodeState
{
    std::size_t node_index = 0;
    std::string node_name;
    double area_m2 = 0.0;
    double patch_potential_v = 0.0;
    double uncoupled_reference_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double sheath_capacitance_f = 0.0;
    double sheath_charge_c = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double particle_charge_c = 0.0;
    double particle_flux_a = 0.0;
};

struct GlobalSheathFieldSolveState
{
    bool active = false;
    std::string solve_mode = "uninitialized";
    double global_patch_potential_v = 0.0;
    double global_reference_potential_v = 0.0;
    double effective_sheath_length_m = 1.0e-3;
    double global_normal_electric_field_v_per_m = 0.0;
    double global_local_charge_density_c_per_m3 = 0.0;
    double field_residual_v_per_m = 0.0;
    double particle_field_coupled_residual_v = 0.0;
    double multi_step_stability_metric_v = 0.0;
    double linear_residual_norm_v = 0.0;
    double matrix_row_count = 0.0;
    double matrix_nonzeros = 0.0;
    std::vector<GlobalSheathFieldNodeState> nodes;
};

class SurfaceCurrentModel
{
  public:
    virtual ~SurfaceCurrentModel() = default;

    virtual SurfaceCurrents evaluate(ReferenceSurfaceRole role,
                                     const SurfaceModelRuntimeState& state) const = 0;
    virtual double solveEquilibriumPotential(ReferenceSurfaceRole role,
                                             const SurfaceModelRuntimeState& state,
                                             double minimum_potential_v,
                                             double maximum_potential_v) const = 0;
    virtual double computeCurrentDerivative(ReferenceSurfaceRole role,
                                            const SurfaceModelRuntimeState& state) const = 0;
    virtual bool recalibrate(const SurfaceModelRuntimeState& state,
                             std::size_t pic_calibration_samples) = 0;
    virtual double electronCalibrationFactor() const = 0;
    virtual double electronCalibrationFactor(const SurfaceModelRuntimeState& state) const = 0;
    virtual double ionCalibrationFactor() const = 0;
    virtual double ionCalibrationFactor(const SurfaceModelRuntimeState& state) const = 0;
    virtual const PicMccCurrentSample& latestLivePicSample() const = 0;
    virtual const PicMccCurrentSample& latestLivePicSample(
        const SurfaceModelRuntimeState& state) const = 0;
    virtual const SurfaceKernelSnapshot& latestLivePicKernelSnapshot() const = 0;
    virtual const SurfaceKernelSnapshot& latestLivePicKernelSnapshot(
        const SurfaceModelRuntimeState& state) const = 0;
    virtual std::string algorithmName() const = 0;
    virtual SurfaceCurrentAlgorithmMode algorithmMode() const = 0;
};

class SurfaceCapacitanceModel
{
  public:
    virtual ~SurfaceCapacitanceModel() = default;

    virtual double computeCapacitancePerArea(const SurfaceChargingConfig& config) const = 0;
    virtual double computePatchCapacitanceF(const SurfaceChargingConfig& config,
                                            double capacitance_per_area_f_per_m2,
                                            double patch_area_m2) const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceVoltageModel
{
  public:
    virtual ~SurfaceVoltageModel() = default;
    virtual std::string modelName() const = 0;
};

class SurfaceCircuitModel;

class PotentialReferenceModel
{
  public:
    virtual ~PotentialReferenceModel() = default;

    virtual double plasmaReferencePotentialV(const SurfaceChargingConfig& config,
                                             const SurfaceModelRuntimeState& state) const = 0;
    virtual std::string modelName() const = 0;
};

class ElectricFieldProvider
{
  public:
    virtual ~ElectricFieldProvider() = default;
    virtual double normalFieldVPerM(const SurfaceChargingConfig& config,
                                    const SurfaceModelRuntimeState& state) const = 0;
    virtual std::string modelName() const = 0;
};

class VolumeChargeProvider
{
  public:
    virtual ~VolumeChargeProvider() = default;
    virtual double localChargeDensityCPerM3(const SurfaceChargingConfig& config,
                                            const SurfaceModelRuntimeState& state) const = 0;
    virtual std::string modelName() const = 0;
};

class BubbleCapacitanceEstimator
{
  public:
    virtual ~BubbleCapacitanceEstimator() = default;
    virtual double estimateCapacitancePerArea(const SurfaceChargingConfig& config,
                                              const SurfaceModelRuntimeState& state,
                                              double fallback_capacitance_per_area_f_per_m2) const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceReferenceStateProvider
{
  public:
    virtual ~SurfaceReferenceStateProvider() = default;
    virtual double plasmaReferencePotentialV(const SurfaceChargingConfig& config,
                                             const SurfaceModelRuntimeState& state) const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceReferenceGraphPropagator
{
  public:
    virtual ~SurfaceReferenceGraphPropagator() = default;
    virtual double propagatedReferencePotentialV(const SurfaceChargingConfig& config,
                                                 const SurfaceCircuitModel* circuit_model,
                                                 const SurfaceModelRuntimeState& state,
                                                 double base_reference_potential_v) const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceGraphCapacitanceMatrixProvider
{
  public:
    virtual ~SurfaceGraphCapacitanceMatrixProvider() = default;
    virtual double diagonalCapacitanceF(const SurfaceChargingConfig& config,
                                        const SurfaceCircuitModel* circuit_model,
                                        std::size_t node_index) const = 0;
    virtual double mutualCapacitanceF(const SurfaceChargingConfig& config,
                                      const SurfaceCircuitModel* circuit_model,
                                      std::size_t branch_index) const = 0;
    virtual double rowSumCapacitanceF(const SurfaceChargingConfig& config,
                                      const SurfaceCircuitModel* circuit_model,
                                      std::size_t node_index) const = 0;
    virtual double graphCouplingMetric(const SurfaceChargingConfig& config,
                                       const SurfaceCircuitModel* circuit_model) const = 0;
    virtual std::string matrixFamilyName() const = 0;
    virtual std::string solverAdapterHint() const = 0;
    virtual bool exposesMutualMatrix() const = 0;
    virtual bool exposesDiagonalMatrix() const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceFieldSolverAdapter
{
  public:
    virtual ~SurfaceFieldSolverAdapter() = default;
    virtual void adaptFieldState(const SurfaceChargingConfig& config,
                                 const SurfaceCircuitModel* circuit_model,
                                 SurfaceModelRuntimeState& state) const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceVolumetricSolverAdapter
{
  public:
    virtual ~SurfaceVolumetricSolverAdapter() = default;
    virtual void adaptVolumeState(const SurfaceChargingConfig& config,
                                  const SurfaceCircuitModel* circuit_model,
                                  SurfaceModelRuntimeState& state) const = 0;
    virtual std::string projectionFamilyName() const = 0;
    virtual std::string meshFamilyName() const = 0;
    virtual std::string requestSchemaVersion() const = 0;
    virtual std::string resultSchemaVersion() const = 0;
    virtual std::string modelName() const = 0;
    virtual bool supportsBoundaryFaceMapping() const = 0;
    virtual bool supportsProjectionWeights() const = 0;
};

class SurfaceCircuitModel
{
  public:
    virtual ~SurfaceCircuitModel() = default;

    virtual void configure(const SurfaceChargingConfig& config,
                           double capacitance_per_area_f_per_m2, double body_potential_v,
                           double patch_potential_v,
                           double effective_conductivity_s_per_m) = 0;
    virtual void updateNodePotentials(double body_potential_v, double patch_potential_v) = 0;
    virtual void setNodePotential(std::size_t node_index, double potential_v) = 0;
    virtual void setNodeCapacitance(std::size_t node_index, double capacitance_f) = 0;
    virtual void setBranchConductance(std::size_t branch_index, double conductance_s) = 0;
    virtual Coupling::SurfaceCircuitAdvanceResult advanceImplicit(
        double dt, const Coupling::SurfaceCircuitLinearization& linearization,
        double max_delta_potential_v) = 0;
    virtual Coupling::SurfaceCircuitAdvanceResult advanceImplicit(
        double dt, const Coupling::SurfaceCircuitKernelInput& kernel_input,
        double max_delta_potential_v) = 0;
    virtual std::size_t bodyNodeIndex() const = 0;
    virtual std::size_t primaryPatchNodeIndex() const = 0;
    virtual std::size_t primaryPatchToBodyBranchIndex() const = 0;
    virtual std::size_t patchCount() const = 0;
    virtual std::size_t patchNodeIndex(std::size_t patch_ordinal) const = 0;
    virtual std::size_t patchToBodyBranchIndex(std::size_t patch_ordinal) const = 0;
    virtual bool branchUsesDynamicConductance(std::size_t branch_index) const = 0;
    virtual std::size_t nodeCount() const = 0;
    virtual std::size_t branchCount() const = 0;
    virtual double nodeAreaM2(std::size_t node_index) const = 0;
    virtual double bodyAreaM2() const = 0;
    virtual double nodePotential(std::size_t node_index) const = 0;
    virtual bool nodeIsPatch(std::size_t node_index) const = 0;
    virtual std::string nodeName(std::size_t node_index) const = 0;
    virtual std::size_t branchFromNodeIndex(std::size_t branch_index) const = 0;
    virtual std::size_t branchToNodeIndex(std::size_t branch_index) const = 0;
    virtual double branchConductanceS(std::size_t branch_index) const = 0;
    virtual double branchResistanceOhm(std::size_t branch_index) const = 0;
    virtual double branchBiasV(std::size_t branch_index) const = 0;
    virtual std::string branchName(std::size_t branch_index) const = 0;
    virtual std::string modelName() const = 0;
};

class SurfaceScenarioOrchestrator
{
  public:
    virtual ~SurfaceScenarioOrchestrator() = default;

    virtual std::unique_ptr<SurfaceCurrentModel> createCurrentModel() const = 0;
    virtual std::unique_ptr<SurfaceVoltageModel> createVoltageModel() const = 0;
    virtual std::unique_ptr<SurfaceCapacitanceModel> createCapacitanceModel() const = 0;
    virtual std::unique_ptr<PotentialReferenceModel> createPotentialReferenceModel() const = 0;
    virtual std::unique_ptr<ElectricFieldProvider> createElectricFieldProvider() const = 0;
    virtual std::unique_ptr<VolumeChargeProvider> createVolumeChargeProvider() const = 0;
    virtual std::unique_ptr<BubbleCapacitanceEstimator> createBubbleCapacitanceEstimator() const = 0;
    virtual std::unique_ptr<SurfaceReferenceStateProvider> createSurfaceReferenceStateProvider() const = 0;
    virtual std::unique_ptr<SurfaceReferenceGraphPropagator>
    createSurfaceReferenceGraphPropagator() const = 0;
    virtual std::unique_ptr<SurfaceGraphCapacitanceMatrixProvider>
    createSurfaceGraphCapacitanceMatrixProvider() const = 0;
    virtual std::unique_ptr<SurfaceFieldSolverAdapter> createSurfaceFieldSolverAdapter() const = 0;
    virtual std::unique_ptr<SurfaceVolumetricSolverAdapter>
    createSurfaceVolumetricSolverAdapter() const = 0;
    virtual std::unique_ptr<SurfaceCircuitModel> createCircuitModel() const = 0;
};

std::unique_ptr<SurfaceScenarioOrchestrator>
makeSurfaceScenarioOrchestrator(const SurfaceChargingConfig& config);
SurfaceChargingConfig normalizeSurfaceChargingConfig(const SurfaceChargingConfig& config);
bool validateSurfaceChargingConfig(const SurfaceChargingConfig& config,
                                   std::string* error_message = nullptr);

class DensePlasmaSurfaceCharging
{
  public:
    bool initialize(const SurfaceChargingConfig& config);
    bool initialize(const SurfaceRuntimePlan& runtime_plan);
    bool advance(double dt);
    const SurfaceChargingStatus& getStatus() const { return status_; }
    const std::string& lastErrorMessage() const { return last_error_message_; }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

    SurfaceCurrents computeSurfaceCurrents(double surface_potential_v) const;
    double computeFloatingPotential() const;
    double recommendTimeStep(double remaining_time_s, double minimum_dt_s, double maximum_dt_s) const;

  private:
    struct LegacyBenchmarkReplayRow
    {
        double time_s = 0.0;
        double body_potential_v = 0.0;
        double patch_potential_v = 0.0;
        SurfaceCurrents currents;
        double net_current_a_per_m2 = 0.0;
        double branch_current_a = 0.0;
    };

    SurfaceCurrents computeLegacySurfaceCurrents(double surface_potential_v) const;
    ReferenceCurrentBalanceConfig buildReferenceConfig(double conductivity_s_per_m) const;
    SurfaceModelRuntimeState buildRuntimeState(double body_potential_v, double patch_potential_v,
                                               double surface_area_m2,
                                               double effective_conductivity_s_per_m,
                                               double effective_sheath_length_m,
                                               std::size_t node_index = 0,
                                               const std::string& node_name = std::string{}) const;
    double computeCapacitancePerArea() const;
    double computeCapacitancePerArea(const SurfaceModelRuntimeState& state) const;
    double computeBodyNodeCapacitanceF(const SurfaceModelRuntimeState& state) const;
    double computeEffectiveSheathLength() const;
    double computeEffectiveConductivity(double surface_potential_v) const;
    double computeEffectiveConductivity(double surface_potential_v,
                                        const Material::MaterialProperty& material) const;
    bool useSharedSurfacePicRuntime() const;
    double computeSharedSurfacePatchPotentialV() const;
    double computeSharedSurfacePatchAreaM2() const;
    double computeSharedSurfacePatchPotentialSpreadV() const;
    double computeSharedSurfaceEffectiveSheathLengthM(double base_sheath_length_m) const;
    double computeSharedSurfaceReferencePotentialV(double base_reference_potential_v,
                                                  double local_patch_potential_v) const;
    double computeSharedSurfaceGlobalSolveWeight() const;
    bool useSharedSurfaceGlobalCoupledSolve() const;
    std::size_t computeSharedSurfaceGlobalCoupledIterationLimit() const;
    double computeSharedSurfaceGlobalCoupledToleranceV() const;
    bool useSharedSurfaceLivePicCoupledRefresh() const;
    double computeSharedSurfaceLivePicCoupledRefreshThresholdV() const;
    bool useSharedSurfaceParticleTransportCoupling() const;
    double computeSharedSurfaceParticleTransportConductanceSPerM2(
        double shared_live_pic_net_current_density_a_per_m2,
        double shared_live_pic_derivative_a_per_m2_per_v) const;
    void advanceSharedSurfaceParticleTransportState(
        double dt, const PicMccCurrentSample& shared_live_pic_sample, double shared_patch_area_m2,
        double shared_effective_sheath_length_m);
    double computeSharedSurfaceParticleTransportReferenceShiftV(
        double shared_patch_area_m2, double shared_effective_sheath_length_m) const;
    void updateSharedSurfaceDistributedParticleTransportState(
        const std::vector<double>& node_potentials_v, double dt, double shared_patch_area_m2,
        double shared_effective_sheath_length_m,
        double shared_transport_conductance_s_per_m2);
    double sharedSurfaceDistributedParticleTransportChargeC(std::size_t node_index) const;
    double computeSharedSurfaceGlobalParticleDomainChargeConservationErrorC() const;
    double computeSharedSurfaceGlobalParticleDomainFluxConservationErrorA() const;
    double computeSharedSurfaceGlobalParticleDomainEdgeChargeTotalAbsC() const;
    double computeSharedSurfaceDistributedParticleTransportReferenceShiftV(
        std::size_t node_index, double node_area_m2, double shared_effective_sheath_length_m) const;
    double computeSharedSurfaceGlobalSheathFieldResidualVPerM(
        const SurfaceModelRuntimeState& state) const;
    double computeSharedSurfaceGlobalParticleFieldCoupledResidualV(
        const SurfaceModelRuntimeState& state) const;
    double computeSharedSurfaceGlobalSheathFieldMultiStepStabilityMetricV() const;
    double computeTransitionSourceFluxScale() const;
    double computeLeakageCurrentDensity(double surface_potential_v) const;
    double computeNetCurrentDensity(double surface_potential_v) const;
    double computeRadiationInducedConductivity() const;
    bool updateLivePicCalibration(double surface_potential_v);
    void applyGeoPicCalibration();
    void rebuildSurfaceCircuit(double conductivity_s_per_m);
    bool loadLegacyBenchmarkReplay();
    bool loadLegacyBenchmarkExecutionCurve();
    bool advanceLegacyBenchmarkReplay(double dt);
    double estimateElectronFluxPicLike(double surface_potential_v, std::size_t samples) const;
    double estimateIonFluxPicLike(double surface_potential_v, std::size_t samples) const;
    double estimateCurrentDerivative(double surface_potential_v) const;
    double advancePotentialImplicit(double surface_potential_v, double capacitance_per_area_f_per_m2,
                                    double dt) const;
    void appendSurfaceNodeDiagnostics(const std::vector<double>& branch_currents_a,
                                      double effective_sheath_length_m);
    void syncSharedTransportCachesFromOwnedGlobalParticleDomainState();
    void populateOwnedGlobalParticleDomainStateFromTransportBuffers(
        const std::vector<double>& node_potentials_v,
        const std::vector<double>& transport_node_charge_c,
        const std::vector<double>& transport_node_net_flux_a,
        const std::vector<std::vector<double>>& transport_edge_charge_matrix_c,
        const std::vector<std::vector<double>>& transport_edge_target_charge_matrix_c,
        const std::vector<std::vector<double>>& transport_edge_operator_drive_matrix_c,
        const std::vector<std::vector<double>>& transport_exchange_flux_matrix_a,
        const std::vector<std::vector<double>>& transport_conductance_matrix_s);
    void populateOwnedGlobalParticleRepositoryStateFromTransportBuffers(
        const std::vector<double>& node_potentials_v,
        const std::vector<double>& transport_node_charge_c,
        const std::vector<double>& transport_node_net_flux_a,
        const std::vector<double>& target_node_charge_c,
        const std::vector<double>& node_migration_delta_charge_c,
        const std::vector<double>& node_edge_feedback_charge_c,
        const std::vector<double>& node_conservation_correction_charge_c,
        const std::vector<std::vector<double>>& transport_edge_charge_matrix_c,
        const std::vector<std::vector<double>>& transport_edge_target_charge_matrix_c,
        const std::vector<std::vector<double>>& transport_edge_operator_drive_matrix_c,
        const std::vector<std::vector<double>>& transport_exchange_flux_matrix_a,
        const std::vector<std::vector<double>>& transport_conductance_matrix_s,
        double dt);
    void updateOwnedGlobalParticleDomainState(const std::vector<double>& node_potentials_v);
    void rebuildOwnedGlobalParticleDomainEdges();
    void solveOwnedGlobalSheathFieldSystem(const std::vector<double>& node_potentials_v,
                                           double shared_effective_sheath_length_m,
                                           double dt);
    void appendOwnedGlobalCoupledLinearization(
        const std::vector<double>& node_potentials_v, double dt,
        Coupling::SurfaceCircuitLinearization& linearization,
        std::size_t& current_matrix_coupling_offdiag_entries,
        std::size_t& particle_transport_offdiag_entries,
        double& particle_transport_total_conductance_s,
        double& particle_transport_conservation_error_a_per_v) const;
    void refreshGlobalKernelStates(const std::vector<SurfaceModelRuntimeState>& runtime_states);
    const GlobalParticleDomainNodeState* findGlobalParticleDomainNodeState(
        std::size_t node_index) const;
    const GlobalSheathFieldNodeState* findGlobalSheathFieldNodeState(
        std::size_t node_index) const;
    double currentSharedSurfaceParticleTransportChargeC() const;
    double currentSharedSurfaceParticleTransportReferenceShiftV() const;
    SurfaceDidvState assembleSurfaceDidvState(ReferenceSurfaceRole role,
                                              const SurfaceModelRuntimeState& runtime_state,
                                              const SurfaceCurrents& currents) const;
    SurfaceCircuitMappingState buildSurfaceCircuitMappingState() const;

    SurfaceChargingConfig config_;
    SurfaceChargingStatus status_;
    double transition_conductivity_scale_ = 1.0;
    double transition_source_flux_scale_ = 1.0;
    double transition_simulation_param_scale_ = 1.0;
    std::size_t transition_runtime_internal_substeps_ = 1;
    double transition_runtime_max_delta_potential_v_per_step_ = 0.0;
    double transition_runtime_solver_relaxation_factor_ = 1.0;
    ChargeAccumulationModel accumulation_model_;
    Material::SurfaceInteraction surface_interaction_;
    Output::ResultExporter exporter_;
    std::vector<double> history_time_;
    std::vector<double> history_potential_;
    std::vector<double> history_charge_;
    std::vector<double> history_current_;
    std::vector<double> history_electron_current_;
    std::vector<double> history_ion_current_;
    std::vector<double> history_secondary_current_;
    std::vector<double> history_ion_secondary_current_;
    std::vector<double> history_backscatter_current_;
    std::vector<double> history_photo_current_;
    std::vector<double> history_thermionic_current_;
    std::vector<double> history_field_emission_current_;
    std::vector<double> history_leakage_current_;
    std::vector<double> history_ram_current_;
    std::vector<double> history_body_potential_;
    std::vector<double> history_patch_potential_;
    std::vector<double> history_circuit_branch_current_;
    std::vector<double> history_current_derivative_;
    std::vector<double> history_pic_recalibration_marker_;
    std::vector<double> history_live_pic_electron_current_;
    std::vector<double> history_live_pic_ion_current_;
    std::vector<double> history_live_pic_net_current_;
    std::vector<double> history_live_pic_derivative_;
    std::vector<double> history_live_pic_collision_count_;
    std::vector<double> history_live_pic_mcc_enabled_;
    std::vector<double> history_net_current_;
    std::vector<double> history_capacitance_;
    std::vector<double> history_effective_conductivity_;
    std::vector<double> history_effective_sheath_length_;
    std::vector<double> history_shared_surface_patch_potential_;
    std::vector<double> history_shared_surface_patch_area_;
    std::vector<double> history_shared_surface_reference_potential_;
    std::vector<double> history_shared_surface_effective_sheath_length_;
    std::vector<double> history_shared_surface_sheath_charge_;
    std::vector<double> history_shared_surface_runtime_enabled_;
    std::vector<double> history_shared_surface_pre_global_solve_patch_potential_spread_;
    std::vector<double> history_shared_surface_patch_potential_spread_reduction_v_;
    std::vector<double> history_shared_surface_patch_potential_spread_reduction_ratio_;
    std::vector<double> history_shared_surface_current_matrix_coupling_active_;
    std::vector<double> history_shared_surface_current_matrix_coupling_offdiag_entries_;
    std::vector<double> history_shared_surface_global_coupled_solve_active_;
    std::vector<double> history_shared_surface_global_coupled_solve_iterations_;
    std::vector<double> history_shared_surface_global_coupled_solve_converged_;
    std::vector<double> history_shared_surface_global_coupled_solve_max_delta_v_;
    std::vector<double> history_shared_surface_live_pic_coupled_refresh_active_;
    std::vector<double> history_shared_surface_live_pic_coupled_refresh_count_;
    std::vector<double> history_shared_surface_particle_transport_coupling_active_;
    std::vector<double> history_shared_surface_particle_transport_offdiag_entries_;
    std::vector<double> history_shared_surface_particle_transport_total_conductance_s_;
    std::vector<double> history_shared_surface_particle_transport_conservation_error_a_per_v_;
    std::vector<double> history_shared_surface_particle_transport_charge_c_;
    std::vector<double> history_shared_surface_particle_transport_reference_shift_v_;
    std::vector<double> history_normal_electric_field_;
    std::vector<double> history_local_charge_density_;
    std::vector<double> history_adaptive_time_step_;
    std::vector<double> history_internal_substeps_;
    std::vector<double> history_electron_calibration_factor_;
    std::vector<double> history_ion_calibration_factor_;
    std::vector<double> history_equilibrium_potential_;
    std::vector<double> history_equilibrium_error_;
    std::vector<std::vector<double>> history_surface_node_potentials_;
    std::vector<std::vector<double>> history_surface_node_total_currents_;
    std::vector<std::vector<double>> history_surface_node_electron_currents_;
    std::vector<std::vector<double>> history_surface_node_ion_currents_;
    std::vector<std::vector<double>> history_surface_node_secondary_currents_;
    std::vector<std::vector<double>> history_surface_node_ion_secondary_currents_;
    std::vector<std::vector<double>> history_surface_node_backscatter_currents_;
    std::vector<std::vector<double>> history_surface_node_photo_currents_;
    std::vector<std::vector<double>> history_surface_node_thermionic_currents_;
    std::vector<std::vector<double>> history_surface_node_field_emission_currents_;
    std::vector<std::vector<double>> history_surface_node_conduction_currents_;
    std::vector<std::vector<double>> history_surface_node_ram_currents_;
    std::vector<std::vector<double>> history_surface_node_current_derivatives_;
    std::vector<std::vector<double>> history_surface_node_effective_sheath_lengths_;
    std::vector<std::vector<double>> history_surface_node_normal_electric_fields_;
    std::vector<std::vector<double>> history_surface_node_local_charge_densities_;
    std::vector<std::vector<double>> history_surface_node_propagated_reference_potentials_;
    std::vector<std::vector<double>> history_surface_node_field_solver_reference_potentials_;
    std::vector<std::vector<double>> history_surface_node_shared_runtime_enabled_;
    std::vector<std::vector<double>> history_surface_node_shared_patch_potentials_;
    std::vector<std::vector<double>> history_surface_node_shared_patch_areas_;
    std::vector<std::vector<double>> history_surface_node_shared_reference_potentials_;
    std::vector<std::vector<double>> history_surface_node_shared_sheath_charges_;
    std::vector<std::vector<double>> history_surface_node_distributed_particle_transport_charges_;
    std::vector<std::vector<double>> history_surface_node_distributed_particle_transport_reference_shifts_;
    std::vector<std::vector<double>> history_surface_node_distributed_particle_transport_net_fluxes_;
    std::vector<std::vector<double>> history_surface_node_graph_capacitance_diagonals_;
    std::vector<std::vector<double>> history_surface_node_graph_capacitance_row_sums_;
    std::vector<std::vector<double>> history_surface_node_field_solver_coupling_gains_;
    std::vector<std::vector<double>> history_surface_node_field_solver_capacitance_scales_;
    std::vector<std::vector<double>> history_surface_node_pseudo_volumes_;
    std::vector<std::vector<double>> history_surface_node_volume_projection_weight_sums_;
    std::vector<std::vector<double>> history_surface_node_volume_mesh_coupling_gains_;
    std::vector<std::vector<double>> history_surface_node_volume_potentials_;
    std::vector<std::vector<double>> history_surface_node_deposited_charges_;
    std::vector<std::vector<double>> history_surface_node_poisson_residuals_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_mode_ids_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_iterations_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_linear_iterations_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_converged_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_residual_norms_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_max_deltas_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_matrix_nnzs_;
    std::vector<std::vector<double>> history_surface_node_volume_solver_cell_counts_;
    std::vector<std::vector<double>> history_surface_node_field_volume_coupling_iterations_;
    std::vector<std::vector<double>> history_surface_node_field_volume_coupling_converged_;
    std::vector<std::vector<double>> history_surface_node_field_volume_coupling_max_deltas_;
    std::vector<std::vector<double>> history_surface_node_field_volume_coupling_relaxation_used_;
    std::vector<std::vector<double>> history_surface_node_external_volume_feedback_blend_factors_;
    std::vector<std::vector<double>> history_surface_node_external_volume_feedback_mismatch_metrics_;
    std::vector<std::vector<double>> history_surface_node_external_volume_feedback_applied_;
    std::vector<std::vector<double>> history_surface_node_electron_calibration_factors_;
    std::vector<std::vector<double>> history_surface_node_ion_calibration_factors_;
    std::vector<std::vector<double>> history_surface_node_live_pic_electron_currents_;
    std::vector<std::vector<double>> history_surface_node_live_pic_ion_currents_;
    std::vector<std::vector<double>> history_surface_node_live_pic_net_currents_;
    std::vector<std::vector<double>> history_surface_node_live_pic_derivatives_;
    std::vector<std::vector<double>> history_surface_node_live_pic_collision_counts_;
    std::vector<std::vector<double>> history_surface_node_live_pic_mcc_enabled_;
    std::vector<std::vector<double>> history_surface_branch_currents_;
    std::vector<std::vector<double>> history_surface_branch_conductances_;
    std::vector<std::vector<double>> history_surface_branch_voltage_drops_;
    std::vector<std::vector<double>> history_surface_branch_power_w_;
    std::vector<std::vector<double>> history_surface_branch_mutual_capacitances_;
    std::vector<double> circuit_node_potentials_;
    std::unique_ptr<SurfaceScenarioOrchestrator> orchestrator_;
    std::unique_ptr<SurfaceCurrentModel> current_model_;
    std::unique_ptr<SurfaceVoltageModel> voltage_model_;
    std::unique_ptr<SurfaceCapacitanceModel> capacitance_model_;
    std::unique_ptr<PotentialReferenceModel> potential_reference_model_;
    std::unique_ptr<ElectricFieldProvider> electric_field_provider_;
    std::unique_ptr<VolumeChargeProvider> volume_charge_provider_;
    std::unique_ptr<BubbleCapacitanceEstimator> bubble_capacitance_estimator_;
    std::unique_ptr<SurfaceReferenceStateProvider> reference_state_provider_;
    std::unique_ptr<SurfaceReferenceGraphPropagator> reference_graph_propagator_;
    std::unique_ptr<SurfaceGraphCapacitanceMatrixProvider> graph_capacitance_matrix_provider_;
    std::unique_ptr<SurfaceFieldSolverAdapter> field_solver_adapter_;
    std::unique_ptr<SurfaceVolumetricSolverAdapter> volumetric_solver_adapter_;
    std::unique_ptr<SurfaceCircuitModel> circuit_model_;
    std::vector<LegacyBenchmarkReplayRow> legacy_benchmark_replay_;
    std::vector<LegacyBenchmarkReplayRow> legacy_benchmark_body_curve_;
    std::size_t legacy_benchmark_replay_index_ = 0;
    GlobalParticleDomainState global_particle_domain_state_;
    GlobalParticleRepositoryState global_particle_repository_state_;
    GlobalSheathFieldSolveState global_sheath_field_solve_state_;
    double shared_particle_transport_charge_c_ = 0.0;
    double shared_particle_transport_reference_shift_v_ = 0.0;
    std::vector<double> shared_particle_transport_node_charge_c_;
    std::vector<double> shared_particle_transport_node_net_flux_a_;
    std::vector<std::vector<double>> shared_particle_transport_edge_charge_matrix_c_;
    std::vector<std::vector<double>> shared_particle_transport_edge_target_charge_matrix_c_;
    std::vector<std::vector<double>> shared_particle_transport_edge_operator_drive_matrix_c_;
    std::vector<std::vector<double>> shared_particle_transport_exchange_flux_matrix_a_;
    std::vector<std::vector<double>> global_particle_transport_conductance_matrix_s_;
    std::vector<std::vector<double>> global_sheath_field_reference_response_matrix_;
    std::vector<double> global_sheath_field_previous_node_charge_c_;
    double shared_particle_transport_edge_graph_operator_iterations_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_converged_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_effective_pair_count_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_ = 0.0;
    double shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_ = 1.0;
    double shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_ = 1.0;
    mutable double field_volume_adaptive_relaxation_ = 0.65;
    mutable bool field_volume_adaptive_relaxation_initialized_ = false;
    mutable std::string last_error_message_;
    bool initialized_ = false;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
