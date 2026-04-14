#include "DensePlasmaSurfaceCharging.h"
#include "LegacyBenchmarkSupport.h"
#include "SurfaceRuntimePlan.h"
#include "SurfaceTransitionEngine.h"
#include "SurfaceFlowCouplingModel.h"

#include "../../Tools/Basic/include/NumericAggregation.h"
#include "../../Tools/Basic/include/StringTokenUtils.h"
#include "../../Tools/Coupling/include/SurfaceCurrentLinearizationKernel.h"
#include "../../Tools/FieldSolver/include/SurfaceBarrierModels.h"
#include "../../Tools/Material/include/SurfaceMaterialImporter.h"
#include "../../Tools/Particle/include/SurfaceDistributionFunction.h"
#include "../../Tools/Solver/include/SurfaceSolverFacade.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <random>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kAtomicMassUnit = 1.66053906660e-27;
constexpr double kBoltzmannEvPerK = 8.617333262145e-5;
constexpr double kRichardsonConstant = 1.20173e6;
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kPi = 3.14159265358979323846;
constexpr std::size_t kLegacyExecutorMaxSamples = 4096;

bool applyImportedSurfaceMaterial(SurfaceChargingConfig& config, std::string& error_message)
{
    if (config.material_library_path.empty())
    {
        return true;
    }

    Material::MaterialDatabase database;
    database.clear();

    Material::SurfaceMaterialImporter importer;
    const auto import_result = importer.importPath(config.material_library_path, database);
    if (!import_result)
    {
        error_message = import_result.message().empty()
                            ? "Failed to import surface material library: " +
                                  config.material_library_path.string()
                            : import_result.message();
        return false;
    }

    const std::string requested_material_name =
        config.imported_material_name.empty() ? config.material.getName() : config.imported_material_name;
    const auto* imported_material = database.findByName(requested_material_name);
    if (imported_material == nullptr)
    {
        error_message = "Imported surface material '" + requested_material_name +
                        "' was not found in: " + config.material_library_path.string();
        return false;
    }

    config.material = *imported_material;
    config.default_surface_physics.material = config.material;
    const auto normalize_material_key = [](const std::string& text) {
        return Basic::normalizeAlnumToken(text);
    };

    for (auto& patch : config.patches)
    {
        if (!patch.material.has_value())
        {
            continue;
        }

        if (normalize_material_key(patch.material->getName()) ==
            normalize_material_key(requested_material_name))
        {
            patch.material = *imported_material;
        }
    }

    return true;
}

Solver::SolverPolicyFlags resolvedSolverPolicyFlags(const SurfaceChargingConfig& config)
{
    return Solver::resolveSolverPolicyFlags(config.solver_config.coupling_mode,
                                            config.solver_config.convergence_policy);
}

std::string normalizedSolverDepositionScheme(const SurfaceChargingConfig& config)
{
    return Basic::normalizeAlnumToken(config.solver_config.deposition_scheme);
}

std::string normalizedSamplingPolicy(const SurfaceChargingConfig& config)
{
    return Basic::normalizeAlnumToken(config.sampling_policy);
}

void applyBodyMaterialOverrides(const Material::MaterialProperty& source_material,
                                Material::MaterialProperty& body_material)
{
    body_material.setWorkFunctionEv(
        std::max(0.0, source_material.getScalarProperty("body_work_function_ev",
                                                        body_material.getWorkFunctionEv())));
    body_material.setSecondaryElectronYield(std::max(
        0.0, source_material.getScalarProperty("body_secondary_electron_yield",
                                               body_material.getSecondaryElectronYield())));
    body_material.setScalarProperty(
        "atomic_number",
        std::max(0.0, source_material.getScalarProperty(
                          "body_atomic_number",
                          body_material.getScalarProperty("atomic_number", 13.0))));
    body_material.setScalarProperty(
        "secondary_yield_peak_energy_ev",
        std::max(1.0, source_material.getScalarProperty(
                          "body_secondary_yield_peak_energy_ev",
                          body_material.getScalarProperty("secondary_yield_peak_energy_ev", 300.0))));
    body_material.setScalarProperty(
        "secondary_emission_escape_energy_ev",
        std::max(1.0e-3, source_material.getScalarProperty(
                               "body_secondary_emission_escape_energy_ev",
                               body_material.getScalarProperty("secondary_emission_escape_energy_ev",
                                                               2.0))));
    body_material.setScalarProperty(
        "ion_secondary_yield",
        std::max(0.0, source_material.getScalarProperty(
                          "body_ion_secondary_yield",
                          body_material.getScalarProperty("ion_secondary_yield", 0.08))));
    body_material.setScalarProperty(
        "ion_secondary_peak_energy_kev",
        std::max(1.0e-6, source_material.getScalarProperty(
                             "body_ion_secondary_peak_energy_kev",
                             body_material.getScalarProperty("ion_secondary_peak_energy_kev", 0.35))));
    body_material.setScalarProperty(
        "sims_exponent_n",
        std::max(1.0, source_material.getScalarProperty(
                          "body_sims_exponent_n",
                          body_material.getScalarProperty("sims_exponent_n", 1.6))));
    body_material.setScalarProperty(
        "katz_r1",
        source_material.getScalarProperty("body_katz_r1",
                                          body_material.getScalarProperty("katz_r1", 80.0)));
    body_material.setScalarProperty(
        "katz_n1",
        source_material.getScalarProperty("body_katz_n1",
                                          body_material.getScalarProperty("katz_n1", 0.6)));
    body_material.setScalarProperty(
        "katz_r2",
        source_material.getScalarProperty("body_katz_r2",
                                          body_material.getScalarProperty("katz_r2", 200.0)));
    body_material.setScalarProperty(
        "katz_n2",
        source_material.getScalarProperty("body_katz_n2",
                                          body_material.getScalarProperty("katz_n2", 1.7)));
}

SurfaceBenchmarkSource inferBenchmarkSourceFromReferenceContract(
    const SurfaceChargingConfig& config)
{
    if (config.benchmark_source != SurfaceBenchmarkSource::None)
    {
        return config.benchmark_source;
    }

    const auto map_case_id = [](const std::string& case_id) {
        if (case_id == "geo_ecss_kapton_2000s")
        {
            return SurfaceBenchmarkSource::CGeo;
        }
        if (case_id == "leo_ram_facing_2000s")
        {
            return SurfaceBenchmarkSource::CLeoRam;
        }
        if (case_id == "leo_wake_facing_negative_tail_2000s")
        {
            return SurfaceBenchmarkSource::CLeoWake;
        }
        return SurfaceBenchmarkSource::None;
    };

    SurfaceBenchmarkSource inferred =
        map_case_id(config.reference_matrix_case_id);
    if (inferred != SurfaceBenchmarkSource::None)
    {
        return inferred;
    }

    return map_case_id(config.reference_case_id);
}

Solver::GlobalCoupledControl resolveSharedGlobalCoupledControl(
    const SurfaceChargingConfig& config,
    std::size_t configured_iteration_limit,
    double configured_tolerance_v,
    double configured_relaxation)
{
    Solver::GlobalCoupledControlInput input;
    input.configured_iteration_limit = configured_iteration_limit;
    input.configured_tolerance_v = configured_tolerance_v;
    input.configured_relaxation = configured_relaxation;
    input.solver_max_iterations = config.solver_config.max_iterations;
    input.solver_residual_tolerance = config.solver_config.residual_tolerance;
    input.solver_relaxation_factor = config.solver_config.relaxation_factor;
    input.coupling_mode = config.solver_config.coupling_mode;
    input.convergence_policy = config.solver_config.convergence_policy;
    return Solver::resolveGlobalCoupledControl(input);
}

std::size_t resolveSharedGlobalCoupledIterationLimit(const SurfaceChargingConfig& config,
                                                     std::size_t configured_limit)
{
    return resolveSharedGlobalCoupledControl(config, configured_limit, 1.0e-6, 1.0).iteration_limit;
}

double resolveSharedGlobalCoupledToleranceV(const SurfaceChargingConfig& config,
                                            double configured_tolerance_v)
{
    return resolveSharedGlobalCoupledControl(config, 0, configured_tolerance_v, 1.0).tolerance_v;
}

double resolveSharedGlobalCoupledRelaxation(const SurfaceChargingConfig& config,
                                            double configured_relaxation)
{
    return resolveSharedGlobalCoupledControl(config, 0, 1.0e-6, configured_relaxation).relaxation;
}

double clampSigned(double value, double limit)
{
    return std::clamp(value, -limit, limit);
}

double safeExp(double exponent)
{
    return std::exp(std::clamp(exponent, -200.0, 50.0));
}

double emittedElectronEscapeProbability(const Material::MaterialProperty& material,
                                        double surface_potential_v,
                                        double characteristic_energy_ev)
{
    FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = surface_potential_v;
    state.reference_potential_v = 0.0;
    state.barrier_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 0.0;
    state.emission_temperature_ev = std::max(1.0e-3, characteristic_energy_ev);

    FieldSolver::VariableBarrierScaler scaler;
    const auto evaluation = scaler.evaluate(material, state, 1.0);
    return evaluation.valid ? std::clamp(evaluation.scaling, 0.0, 1.0) : 0.0;
}

double halfNormalFluxAverage(double thermal_sigma, std::mt19937_64& generator)
{
    std::normal_distribution<double> distribution(0.0, thermal_sigma);
    return std::abs(distribution(generator));
}

void synchronizePlasmaMomentsFromSpectra(SurfaceChargingConfig& config)
{
    if (config.has_electron_spectrum)
    {
        config.plasma.electron_density_m3 =
            Particle::resolvedSpectrumDensity(config.electron_spectrum,
                                             config.plasma.electron_density_m3);
        config.plasma.electron_temperature_ev =
            Particle::resolvedSpectrumCharacteristicEnergyEv(
                config.electron_spectrum, config.plasma.electron_temperature_ev);
    }
    if (config.has_ion_spectrum)
    {
        config.plasma.ion_density_m3 =
            Particle::resolvedSpectrumDensity(config.ion_spectrum,
                                             config.plasma.ion_density_m3);
        config.plasma.ion_temperature_ev =
            Particle::resolvedSpectrumCharacteristicEnergyEv(
                config.ion_spectrum, config.plasma.ion_temperature_ev);
        config.plasma.ion_mass_amu =
            Particle::resolvedSpectrumAverageMassAmu(config.ion_spectrum,
                                                    config.plasma.ion_mass_amu);
    }
}

const SurfacePatchPhysicsConfig* findPatchPhysicsOverride(const SurfaceChargingConfig& config,
                                                          std::size_t node_index,
                                                          const std::string& node_name)
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

const Material::MaterialProperty& resolvePatchMaterial(const SurfaceChargingConfig& config,
                                                       std::size_t node_index,
                                                       const std::string& node_name)
{
    const auto* patch_config = findPatchPhysicsOverride(config, node_index, node_name);
    if (patch_config != nullptr && patch_config->override_material)
    {
        return patch_config->material;
    }
    return config.material;
}

double runtimeRouteId(SurfaceRuntimeRoute route)
{
    switch (route)
    {
    case SurfaceRuntimeRoute::SurfacePic:
        return 1.0;
    case SurfaceRuntimeRoute::SurfacePicHybrid:
        return 2.0;
    case SurfaceRuntimeRoute::LegacyBenchmark:
        return 3.0;
    case SurfaceRuntimeRoute::SCDATUnified:
    default:
        return 0.0;
    }
}

double surfacePicStrategyId(SurfacePicStrategy strategy)
{
    switch (strategy)
    {
    case SurfacePicStrategy::SurfacePicDirect:
        return 1.0;
    case SurfacePicStrategy::SurfacePicHybridReference:
        return 2.0;
    case SurfacePicStrategy::SurfacePicCalibrated:
    default:
        return 0.0;
    }
}

double legacyInputAdapterId(SurfaceLegacyInputAdapterKind adapter)
{
    switch (adapter)
    {
    case SurfaceLegacyInputAdapterKind::CTextReferenceDeck:
        return 1.0;
    case SurfaceLegacyInputAdapterKind::MatlabReferenceDeck:
        return 2.0;
    case SurfaceLegacyInputAdapterKind::None:
    default:
        return 0.0;
    }
}

double surfacePicRuntimeId(SurfacePicRuntimeKind runtime)
{
    switch (runtime)
    {
    case SurfacePicRuntimeKind::GraphCoupledSharedSurface:
        return 1.0;
    case SurfacePicRuntimeKind::LocalWindowSampler:
    default:
        return 0.0;
    }
}

double surfaceInstrumentSetId(SurfaceInstrumentSetKind instrument_set)
{
    switch (instrument_set)
    {
    case SurfaceInstrumentSetKind::SurfacePicObserverSet:
        return 1.0;
    case SurfaceInstrumentSetKind::MetadataOnly:
    default:
        return 0.0;
    }
}

double benchmarkSourceId(SurfaceBenchmarkSource source)
{
    switch (source)
    {
    case SurfaceBenchmarkSource::CGeo:
        return 1.0;
    case SurfaceBenchmarkSource::CLeoRam:
        return 2.0;
    case SurfaceBenchmarkSource::CLeoWake:
        return 3.0;
    case SurfaceBenchmarkSource::MatlabGeo:
        return 4.0;
    case SurfaceBenchmarkSource::MatlabLeo:
        return 5.0;
    default:
        return 0.0;
    }
}

std::string runtimeRouteName(SurfaceRuntimeRoute route)
{
    switch (route)
    {
    case SurfaceRuntimeRoute::SurfacePic:
        return "SurfacePic";
    case SurfaceRuntimeRoute::SurfacePicHybrid:
        return "SurfacePicHybrid";
    case SurfaceRuntimeRoute::LegacyBenchmark:
        return "LegacyBenchmark";
    case SurfaceRuntimeRoute::SCDATUnified:
    default:
        return "SCDATUnified";
    }
}

std::string surfacePicStrategyName(SurfacePicStrategy strategy)
{
    switch (strategy)
    {
    case SurfacePicStrategy::SurfacePicDirect:
        return "SurfacePicDirect";
    case SurfacePicStrategy::SurfacePicHybridReference:
        return "SurfacePicHybridReference";
    case SurfacePicStrategy::SurfacePicCalibrated:
    default:
        return "SurfacePicCalibrated";
    }
}

std::string legacyInputAdapterName(SurfaceLegacyInputAdapterKind adapter)
{
    switch (adapter)
    {
    case SurfaceLegacyInputAdapterKind::CTextReferenceDeck:
        return "c_text_reference_deck";
    case SurfaceLegacyInputAdapterKind::MatlabReferenceDeck:
        return "matlab_reference_deck";
    case SurfaceLegacyInputAdapterKind::None:
    default:
        return "none";
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

std::string surfaceInstrumentSetName(SurfaceInstrumentSetKind instrument_set)
{
    switch (instrument_set)
    {
    case SurfaceInstrumentSetKind::SurfacePicObserverSet:
        return "surface_pic_observer_set";
    case SurfaceInstrumentSetKind::MetadataOnly:
    default:
        return "metadata_only";
    }
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

std::string nodeTypeLabel(bool is_patch)
{
    return is_patch ? "patch" : "body";
}

std::string nodeMaterialName(const SurfaceChargingConfig& config,
                             std::size_t node_index,
                             const std::string& node_name,
                             bool is_patch)
{
    if (!is_patch)
    {
        return "conductor-body";
    }
    return resolvePatchMaterial(config, node_index, node_name).getName();
}

std::string nodeOwnerLabel(const SurfaceChargingConfig& config, const std::string& node_name, bool is_patch)
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

double latestHistoryValue(const std::vector<std::vector<double>>& history, std::size_t index)
{
    if (index >= history.size() || history[index].empty())
    {
        return 0.0;
    }
    return history[index].back();
}

std::string branchTypeLabel(const SurfaceCircuitModel* circuit_model, std::size_t branch_index)
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

std::string runtimeKernelSourceFamily(const SurfaceChargingConfig& config)
{
    if (config.runtime_route == SurfaceRuntimeRoute::SCDATUnified)
    {
        return "spis_equivalent_surface_kernel_v1";
    }
    if (config.runtime_route == SurfaceRuntimeRoute::SurfacePicHybrid ||
        config.surface_pic_strategy == SurfacePicStrategy::SurfacePicHybridReference)
    {
        return "spis_equivalent_surface_pic_hybrid_kernel_v1";
    }
    if (config.runtime_route == SurfaceRuntimeRoute::SurfacePic)
    {
        return "spis_equivalent_surface_pic_kernel_v1";
    }
    return "spis_reference_benchmark_kernel_v1";
}

std::string surfaceMaterialModelFamily(const SurfaceChargingConfig& config)
{
    return Material::resolveSurfaceMaterialModelFamily(config.material);
}

std::string barrierScalerFamily(const SurfaceChargingConfig& config)
{
    return config.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark
               ? "reference_barrier_scaler_v1"
               : "spis_variable_barrier_scaler_v1";
}

std::string distributionFamily(const SurfaceChargingConfig& config)
{
    switch (config.distribution_model)
    {
    case PlasmaAnalysis::PlasmaDistributionModelKind::MaxwellianProjected:
        return "spis_isotropic_maxwellian_flux_v1";
    case PlasmaAnalysis::PlasmaDistributionModelKind::WakeAnisotropic:
        return "spis_wake_anisotropic_flux_v1";
    case PlasmaAnalysis::PlasmaDistributionModelKind::MultiPopulationHybrid:
    default:
        return "spis_bimaxwellian_flux_v1";
    }
}

bool referenceCompletionUsed(const SurfaceChargingConfig& config,
                             const SurfaceKernelSnapshot& snapshot)
{
    if (config.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
    {
        return true;
    }
    if (config.surface_pic_strategy == SurfacePicStrategy::SurfacePicHybridReference)
    {
        return true;
    }
    return snapshot.source_family.find("reference") != std::string::npos;
}

std::string nativeComponentAssemblyFamily(const SurfaceChargingConfig& config)
{
    return config.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark
               ? "legacy_reference_component_assembly_v1"
               : "tools_spis_component_assembly_v1";
}

bool referenceComponentFallbackUsed(const SurfaceChargingConfig& config)
{
    return config.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark;
}

std::string resolvedRuntimeKernelSourceFamily(const SurfaceChargingConfig& config,
                                              const SurfaceKernelSnapshot& snapshot)
{
    const std::string configured_family = runtimeKernelSourceFamily(config);
    if (config.runtime_route == SurfaceRuntimeRoute::SCDATUnified)
    {
        return configured_family;
    }
    const bool formal_pic_route = config.runtime_route == SurfaceRuntimeRoute::SurfacePic ||
                                  config.runtime_route == SurfaceRuntimeRoute::SurfacePicHybrid;
    if (formal_pic_route)
    {
        if (snapshot.valid && !snapshot.source_family.empty() &&
            snapshot.source_family != "spis_reference_benchmark_kernel_v1")
        {
            return snapshot.source_family;
        }
        return configured_family;
    }
    if (!snapshot.source_family.empty())
    {
        return snapshot.source_family;
    }
    return configured_family;
}

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

SurfaceGraphMonitorSnapshot buildSurfaceGraphMonitorSnapshot(
    const SurfaceChargingConfig& config, const SurfaceCircuitModel* circuit_model,
    const SurfaceGraphCapacitanceMatrixProvider* graph_capacitance_matrix_provider,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& branch_mutual_capacitance_history,
    const std::vector<std::vector<double>>& branch_conductance_history,
    const std::vector<std::vector<double>>& branch_voltage_drop_history)
{
    SurfaceGraphMonitorSnapshot snapshot;
    if (circuit_model == nullptr)
    {
        return snapshot;
    }

    snapshot.node_count = circuit_model->nodeCount();
    snapshot.branch_count = circuit_model->branchCount();
    std::vector<std::size_t> node_degrees(snapshot.node_count, 0);
    for (std::size_t branch_index = 0; branch_index < snapshot.branch_count; ++branch_index)
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

        const double abs_voltage_drop_v =
            std::abs(latestHistoryValue(branch_voltage_drop_history, branch_index));
        const double mutual_capacitance_f =
            std::abs(latestHistoryValue(branch_mutual_capacitance_history, branch_index));
        const double conductance_s =
            std::max(0.0, latestHistoryValue(branch_conductance_history, branch_index));
        const double coupling_metric =
            mutual_capacitance_f * conductance_s + abs_voltage_drop_v * conductance_s;
        if (coupling_metric > snapshot.strongest_coupling_interface_metric)
        {
            snapshot.strongest_coupling_interface_metric = coupling_metric;
            snapshot.strongest_coupling_interface_index = branch_index;
        }

        snapshot.max_abs_neighbor_potential_delta_v = std::max(
            snapshot.max_abs_neighbor_potential_delta_v,
            std::abs(latestHistoryValue(node_potential_history, from_node) -
                     latestHistoryValue(node_potential_history, to_node)));
        snapshot.max_abs_neighbor_field_contrast_v_per_m = std::max(
            snapshot.max_abs_neighbor_field_contrast_v_per_m,
            std::abs(latestHistoryValue(node_field_history, from_node) -
                     latestHistoryValue(node_field_history, to_node)));
    }

    for (std::size_t node_index = 0; node_index < snapshot.node_count; ++node_index)
    {
        const double abs_field = std::abs(latestHistoryValue(node_field_history, node_index));
        if (abs_field > snapshot.strongest_field_node_abs_v_per_m)
        {
            snapshot.strongest_field_node_abs_v_per_m = abs_field;
            snapshot.strongest_field_node_index = node_index;
        }
    }
    if (!node_degrees.empty())
    {
        snapshot.max_node_degree = *std::max_element(node_degrees.begin(), node_degrees.end());
        snapshot.average_node_degree = Basic::meanValue(node_degrees, 0.0);
    }
    snapshot.graph_coupling_metric =
        graph_capacitance_matrix_provider != nullptr
            ? graph_capacitance_matrix_provider->graphCouplingMetric(config, circuit_model)
            : 0.0;
    return snapshot;
}

SurfaceFieldMonitorSnapshot buildSurfaceFieldMonitorSnapshot(
    const SurfaceCircuitModel* circuit_model,
    const std::vector<std::vector<double>>& node_potential_history,
    const std::vector<std::vector<double>>& node_field_history,
    const std::vector<std::vector<double>>& node_propagated_reference_history,
    const std::vector<std::vector<double>>& node_field_solver_reference_history,
    const std::vector<std::vector<double>>& node_field_solver_coupling_gain_history)
{
    SurfaceFieldMonitorSnapshot snapshot;
    if (circuit_model == nullptr)
    {
        return snapshot;
    }
    for (std::size_t node_index = 0; node_index < circuit_model->nodeCount(); ++node_index)
    {
        const double abs_field = std::abs(latestHistoryValue(node_field_history, node_index));
        if (abs_field > snapshot.strongest_field_node_abs_v_per_m)
        {
            snapshot.strongest_field_node_abs_v_per_m = abs_field;
            snapshot.strongest_field_node_index = node_index;
        }
        const double potential_v = latestHistoryValue(node_potential_history, node_index);
        const double propagated_ref_v =
            latestHistoryValue(node_propagated_reference_history, node_index);
        const double field_solver_ref_v =
            latestHistoryValue(node_field_solver_reference_history, node_index);
        snapshot.max_abs_node_reference_offset_v = std::max(
            snapshot.max_abs_node_reference_offset_v, std::abs(propagated_ref_v - potential_v));
        snapshot.max_abs_node_field_solver_reference_offset_v = std::max(
            snapshot.max_abs_node_field_solver_reference_offset_v,
            std::abs(field_solver_ref_v - potential_v));
        snapshot.max_field_solver_coupling_gain = std::max(
            snapshot.max_field_solver_coupling_gain,
            std::abs(latestHistoryValue(node_field_solver_coupling_gain_history, node_index)));
    }
    snapshot.reference_offset_envelope_v = std::max(
        snapshot.max_abs_node_reference_offset_v,
        snapshot.max_abs_node_field_solver_reference_offset_v);
    return snapshot;
}

std::string jsonEscape(const std::string& value);

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

bool writeSurfaceGraphReport(
    const std::filesystem::path& report_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model, const std::string& runtime_route_name,
    const std::string& circuit_model_name, const std::string& electric_field_provider_name,
    const std::string& reference_graph_propagator_name,
    const std::string& graph_capacitance_provider_name, const std::string& field_solver_adapter_name,
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
        const double diagonal_f =
            latestHistoryValue(node_graph_capacitance_history, node_index);
        total_graph_capacitance_diagonal_f += diagonal_f;
        max_graph_capacitance_diagonal_f =
            std::max(max_graph_capacitance_diagonal_f, std::abs(diagonal_f));
        const double row_sum_f =
            latestHistoryValue(node_graph_capacitance_row_sum_history, node_index);
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

        const auto branch_type = branchTypeLabel(circuit_model, branch_index);
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

        const double mutual_capacitance_f =
            latestHistoryValue(branch_mutual_capacitance_history, branch_index);
        total_mutual_capacitance_f += mutual_capacitance_f;
        max_mutual_capacitance_f =
            std::max(max_mutual_capacitance_f, std::abs(mutual_capacitance_f));
        const double abs_current_a =
            std::abs(latestHistoryValue(branch_current_history, branch_index));
        max_abs_interface_current_a =
            std::max(max_abs_interface_current_a,
                     abs_current_a);
        const double abs_voltage_drop_v =
            std::abs(latestHistoryValue(branch_voltage_drop_history, branch_index));
        max_abs_interface_voltage_drop_v =
            std::max(max_abs_interface_voltage_drop_v,
                     abs_voltage_drop_v);
        const double abs_power_w =
            std::abs(latestHistoryValue(branch_power_history, branch_index));
        max_abs_interface_power_w =
            std::max(max_abs_interface_power_w,
                     abs_power_w);
        max_interface_conductance_s =
            std::max(max_interface_conductance_s,
                     latestHistoryValue(branch_conductance_history, branch_index));
        const double from_potential_v = latestHistoryValue(node_potential_history, from_node);
        const double to_potential_v = latestHistoryValue(node_potential_history, to_node);
        max_abs_neighbor_potential_delta_v =
            std::max(max_abs_neighbor_potential_delta_v, std::abs(from_potential_v - to_potential_v));
        const double from_field_v_per_m = latestHistoryValue(node_field_history, from_node);
        const double to_field_v_per_m = latestHistoryValue(node_field_history, to_node);
        max_abs_neighbor_field_contrast_v_per_m =
            std::max(max_abs_neighbor_field_contrast_v_per_m,
                     std::abs(from_field_v_per_m - to_field_v_per_m));
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
        const double field_magnitude_v_per_m =
            std::abs(latestHistoryValue(node_field_history, node_index));
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
            std::max(max_abs_node_field_solver_reference_offset_v,
                     std::abs(field_solver_reference_offset_v));
        max_field_solver_coupling_gain =
            std::max(max_field_solver_coupling_gain,
                     std::abs(latestHistoryValue(node_field_solver_coupling_gain_history, node_index)));
    }

    const std::size_t max_node_degree =
        node_degrees.empty() ? 0 : *std::max_element(node_degrees.begin(), node_degrees.end());
    const double average_node_degree =
        node_degrees.empty() ? 0.0 : Basic::meanValue(node_degrees, 0.0);

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
    output << "max_abs_neighbor_field_contrast_v_per_m="
           << max_abs_neighbor_field_contrast_v_per_m << "\n";
    output << "max_abs_node_reference_offset_v=" << max_abs_node_reference_offset_v << "\n";
    output << "max_abs_node_field_solver_reference_offset_v="
           << max_abs_node_field_solver_reference_offset_v << "\n";
    output << "max_field_solver_coupling_gain=" << max_field_solver_coupling_gain << "\n";
    output << "graph_coupling_metric="
           << (graph_capacitance_matrix_provider != nullptr
                   ? graph_capacitance_matrix_provider->graphCouplingMetric(config, circuit_model)
                   : 0.0)
           << "\n";
    if (node_count > 0)
    {
        output << "strongest_field_node_index=" << strongest_field_node_index << "\n";
        output << "strongest_field_node_name=" << circuit_model->nodeName(strongest_field_node_index)
               << "\n";
        output << "strongest_field_node_abs_field_v_per_m="
               << strongest_field_node_magnitude_v_per_m << "\n";
    }
    if (branch_count > 0)
    {
        output << "strongest_coupling_interface_index=" << strongest_coupling_branch_index << "\n";
        output << "strongest_coupling_interface_name="
               << circuit_model->branchName(strongest_coupling_branch_index) << "\n";
        output << "strongest_coupling_interface_metric="
               << strongest_coupling_branch_metric << "\n";
    }

    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        const double reference_offset_v =
            latestHistoryValue(node_propagated_reference_history, node_index) -
            latestHistoryValue(node_potential_history, node_index);
        const double field_solver_reference_offset_v =
            latestHistoryValue(node_field_solver_reference_history, node_index) -
            latestHistoryValue(node_potential_history, node_index);
        const double local_charge_density_c_per_m3 =
            latestHistoryValue(node_charge_history, node_index);
        const double field_to_charge_ratio =
            std::abs(local_charge_density_c_per_m3) <= 1.0e-30
                ? 0.0
                : latestHistoryValue(node_field_history, node_index) / local_charge_density_c_per_m3;
        output << "node_" << node_index << "="
               << "name:" << circuit_model->nodeName(node_index)
               << ",type:" << nodeTypeLabel(circuit_model->nodeIsPatch(node_index))
               << ",owner:"
               << nodeOwnerLabel(config, circuit_model->nodeName(node_index),
                                 circuit_model->nodeIsPatch(node_index))
               << ",material:"
               << nodeMaterialName(config, node_index, circuit_model->nodeName(node_index),
                                   circuit_model->nodeIsPatch(node_index))
               << ",degree:" << node_degrees[node_index]
               << ",area_m2:" << circuit_model->nodeAreaM2(node_index)
               << ",last_potential_v:" << latestHistoryValue(node_potential_history, node_index)
               << ",last_total_current_density_a_per_m2:"
               << latestHistoryValue(node_total_current_history, node_index)
               << ",last_normal_electric_field_v_per_m:"
               << latestHistoryValue(node_field_history, node_index)
               << ",last_local_charge_density_c_per_m3:"
               << local_charge_density_c_per_m3
               << ",last_propagated_reference_potential_v:"
               << latestHistoryValue(node_propagated_reference_history, node_index)
               << ",last_reference_offset_v:" << reference_offset_v
               << ",last_field_solver_reference_potential_v:"
               << latestHistoryValue(node_field_solver_reference_history, node_index)
               << ",last_field_solver_reference_offset_v:" << field_solver_reference_offset_v
               << ",last_graph_capacitance_diagonal_f:"
               << latestHistoryValue(node_graph_capacitance_history, node_index)
               << ",last_graph_capacitance_row_sum_f:"
               << latestHistoryValue(node_graph_capacitance_row_sum_history, node_index)
               << ",last_field_solver_coupling_gain:"
               << latestHistoryValue(node_field_solver_coupling_gain_history, node_index)
               << ",last_field_solver_capacitance_scale:"
               << latestHistoryValue(node_field_solver_capacitance_scale_history, node_index)
               << ",field_to_charge_ratio:" << field_to_charge_ratio << "\n";
    }

    for (std::size_t branch_index = 0; branch_index < branch_count; ++branch_index)
    {
        const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
        const auto to_node = circuit_model->branchToNodeIndex(branch_index);
        output << "interface_" << branch_index << "="
               << "name:" << circuit_model->branchName(branch_index)
               << ",type:" << branchTypeLabel(circuit_model, branch_index)
               << ",from_node:" << circuit_model->nodeName(from_node)
               << ",to_node:" << circuit_model->nodeName(to_node)
               << ",conductance_s:" << circuit_model->branchConductanceS(branch_index)
               << ",resistance_ohm:" << circuit_model->branchResistanceOhm(branch_index)
               << ",bias_v:" << circuit_model->branchBiasV(branch_index)
               << ",dynamic_conductance:"
               << (circuit_model->branchUsesDynamicConductance(branch_index) ? 1 : 0)
               << ",last_current_a:" << latestHistoryValue(branch_current_history, branch_index)
               << ",last_voltage_drop_v:"
               << latestHistoryValue(branch_voltage_drop_history, branch_index)
               << ",last_power_w:" << latestHistoryValue(branch_power_history, branch_index)
               << ",last_mutual_capacitance_f:"
               << latestHistoryValue(branch_mutual_capacitance_history, branch_index) << "\n";
    }
    return true;
}

std::string jsonEscape(const std::string& value);

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
               << graph_capacitance_matrix_provider->diagonalCapacitanceF(
                      config, circuit_model, node_index)
               << ","
               << graph_capacitance_matrix_provider->rowSumCapacitanceF(
                      config, circuit_model, node_index)
               << ",row_sum_f\n";
    }
    for (std::size_t branch_index = 0; branch_index < circuit_model->branchCount(); ++branch_index)
    {
        const auto from_node = circuit_model->branchFromNodeIndex(branch_index);
        const auto to_node = circuit_model->branchToNodeIndex(branch_index);
        output << "branch_mutual," << branch_index << "," << circuit_model->branchName(branch_index)
               << "," << circuit_model->nodeName(from_node) << ","
               << circuit_model->nodeName(to_node) << ","
               << graph_capacitance_matrix_provider->mutualCapacitanceF(
                      config, circuit_model, branch_index)
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
    output << "expects_graph_capacitance_inputs="
           << (graph_capacitance_matrix_provider ? 1 : 0) << "\n";
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
    output << "  \"surface_pic_strategy\": \""
           << surfacePicStrategyName(config.surface_pic_strategy) << "\",\n";
    output << "  \"shared_runtime_enabled\": "
           << (latest_scalar(shared_runtime_enabled_history) >= 0.5 ? "true" : "false") << ",\n";
    output << "  \"sample_count\": " << time_history.size() << ",\n";
    output << "  \"latest_time_s\": " << latest_scalar(time_history) << ",\n";
    output << "  \"shared_surface\": {\n";
    output << "    \"latest_patch_potential_v\": "
           << latest_scalar(shared_patch_potential_history) << ",\n";
    output << "    \"latest_patch_area_m2\": " << latest_scalar(shared_patch_area_history) << ",\n";
    output << "    \"latest_reference_potential_v\": "
           << latest_scalar(shared_reference_potential_history) << ",\n";
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
                   << (latestHistoryValue(node_shared_runtime_enabled_history, node_index) >= 0.5
                           ? "true"
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

bool writeSharedSurfaceRuntimeConsistencyJson(
    const std::filesystem::path& json_path, const SurfaceChargingConfig& config,
    const SurfaceCircuitModel* circuit_model,
    const std::vector<double>& time_history,
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

            const double shared_reference =
                latestHistoryValue(node_shared_reference_potential_history, node_index);
            const double shared_sheath_charge =
                latestHistoryValue(node_shared_sheath_charge_history, node_index);
            const double patch_potential =
                latestHistoryValue(node_potential_history, node_index);
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
            const double distributed_particle_transport_reference_shift_v = latestHistoryValue(
                node_distributed_particle_transport_reference_shift_history, node_index);
            const double distributed_particle_transport_net_flux_a = latestHistoryValue(
                node_distributed_particle_transport_net_flux_history, node_index);
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
                for (const double operator_drive_charge_c :
                     edge_particle_transport_operator_drive_matrix_c[node_index])
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
                min_live_pic_electron_current_a_per_m2 =
                    max_live_pic_electron_current_a_per_m2 =
                        live_pic_electron_current_a_per_m2;
                min_live_pic_ion_current_a_per_m2 =
                    max_live_pic_ion_current_a_per_m2 = live_pic_ion_current_a_per_m2;
                min_live_pic_net_current_a_per_m2 =
                    max_live_pic_net_current_a_per_m2 = live_pic_net_current_a_per_m2;
                min_live_pic_collision_count = max_live_pic_collision_count =
                    live_pic_collision_count;
                min_electron_calibration_factor = max_electron_calibration_factor =
                    electron_calibration_factor;
                min_ion_calibration_factor = max_ion_calibration_factor =
                    ion_calibration_factor;
                min_normal_electric_field_v_per_m =
                    max_normal_electric_field_v_per_m = normal_electric_field_v_per_m;
                min_local_charge_density_c_per_m3 =
                    max_local_charge_density_c_per_m3 = local_charge_density_c_per_m3;
                min_distributed_particle_transport_charge_c =
                    max_distributed_particle_transport_charge_c =
                        distributed_particle_transport_charge_c;
                min_distributed_particle_transport_reference_shift_v =
                    max_distributed_particle_transport_reference_shift_v =
                        distributed_particle_transport_reference_shift_v;
                sum_distributed_particle_transport_charge_c =
                    distributed_particle_transport_charge_c;
                min_distributed_particle_transport_net_flux_a =
                    max_distributed_particle_transport_net_flux_a =
                        distributed_particle_transport_net_flux_a;
                sum_distributed_particle_transport_net_flux_a =
                    distributed_particle_transport_net_flux_a;
                min_node_edge_transport_charge_c = max_node_edge_transport_charge_c =
                    node_edge_transport_charge_c;
                sum_node_edge_transport_charge_c = node_edge_transport_charge_c;
                total_node_edge_operator_drive_charge_c = node_edge_operator_drive_charge_c;
                first_shared_node = false;
            }
            else
            {
                min_reference_potential_v =
                    std::min(min_reference_potential_v, shared_reference);
                max_reference_potential_v =
                    std::max(max_reference_potential_v, shared_reference);
                min_sheath_charge_c = std::min(min_sheath_charge_c, shared_sheath_charge);
                max_sheath_charge_c = std::max(max_sheath_charge_c, shared_sheath_charge);
                min_patch_potential_v = std::min(min_patch_potential_v, patch_potential);
                max_patch_potential_v = std::max(max_patch_potential_v, patch_potential);
                min_live_pic_electron_current_a_per_m2 =
                    std::min(min_live_pic_electron_current_a_per_m2,
                             live_pic_electron_current_a_per_m2);
                max_live_pic_electron_current_a_per_m2 =
                    std::max(max_live_pic_electron_current_a_per_m2,
                             live_pic_electron_current_a_per_m2);
                min_live_pic_ion_current_a_per_m2 =
                    std::min(min_live_pic_ion_current_a_per_m2,
                             live_pic_ion_current_a_per_m2);
                max_live_pic_ion_current_a_per_m2 =
                    std::max(max_live_pic_ion_current_a_per_m2,
                             live_pic_ion_current_a_per_m2);
                min_live_pic_net_current_a_per_m2 =
                    std::min(min_live_pic_net_current_a_per_m2,
                             live_pic_net_current_a_per_m2);
                max_live_pic_net_current_a_per_m2 =
                    std::max(max_live_pic_net_current_a_per_m2,
                             live_pic_net_current_a_per_m2);
                min_live_pic_collision_count =
                    std::min(min_live_pic_collision_count, live_pic_collision_count);
                max_live_pic_collision_count =
                    std::max(max_live_pic_collision_count, live_pic_collision_count);
                min_electron_calibration_factor =
                    std::min(min_electron_calibration_factor, electron_calibration_factor);
                max_electron_calibration_factor =
                    std::max(max_electron_calibration_factor, electron_calibration_factor);
                min_ion_calibration_factor =
                    std::min(min_ion_calibration_factor, ion_calibration_factor);
                max_ion_calibration_factor =
                    std::max(max_ion_calibration_factor, ion_calibration_factor);
                min_normal_electric_field_v_per_m =
                    std::min(min_normal_electric_field_v_per_m, normal_electric_field_v_per_m);
                max_normal_electric_field_v_per_m =
                    std::max(max_normal_electric_field_v_per_m, normal_electric_field_v_per_m);
                min_local_charge_density_c_per_m3 = std::min(
                    min_local_charge_density_c_per_m3, local_charge_density_c_per_m3);
                max_local_charge_density_c_per_m3 = std::max(
                    max_local_charge_density_c_per_m3, local_charge_density_c_per_m3);
                min_distributed_particle_transport_charge_c = std::min(
                    min_distributed_particle_transport_charge_c,
                    distributed_particle_transport_charge_c);
                max_distributed_particle_transport_charge_c = std::max(
                    max_distributed_particle_transport_charge_c,
                    distributed_particle_transport_charge_c);
                min_distributed_particle_transport_reference_shift_v = std::min(
                    min_distributed_particle_transport_reference_shift_v,
                    distributed_particle_transport_reference_shift_v);
                max_distributed_particle_transport_reference_shift_v = std::max(
                    max_distributed_particle_transport_reference_shift_v,
                    distributed_particle_transport_reference_shift_v);
                sum_distributed_particle_transport_charge_c +=
                    distributed_particle_transport_charge_c;
                min_distributed_particle_transport_net_flux_a = std::min(
                    min_distributed_particle_transport_net_flux_a,
                    distributed_particle_transport_net_flux_a);
                max_distributed_particle_transport_net_flux_a = std::max(
                    max_distributed_particle_transport_net_flux_a,
                    distributed_particle_transport_net_flux_a);
                sum_distributed_particle_transport_net_flux_a +=
                    distributed_particle_transport_net_flux_a;
                min_node_edge_transport_charge_c = std::min(
                    min_node_edge_transport_charge_c, node_edge_transport_charge_c);
                max_node_edge_transport_charge_c = std::max(
                    max_node_edge_transport_charge_c, node_edge_transport_charge_c);
                sum_node_edge_transport_charge_c += node_edge_transport_charge_c;
                total_node_edge_operator_drive_charge_c += node_edge_operator_drive_charge_c;
            }
        }

        for (std::size_t patch_i = 0; patch_i < shared_patch_nodes.size(); ++patch_i)
        {
            const auto node_i = shared_patch_nodes[patch_i];
            if (node_i >= edge_particle_transport_charge_matrix_c.size())
            {
                continue;
            }
            for (std::size_t patch_j = patch_i + 1; patch_j < shared_patch_nodes.size(); ++patch_j)
            {
                const auto node_j = shared_patch_nodes[patch_j];
                if (node_j >= edge_particle_transport_charge_matrix_c[node_i].size())
                {
                    continue;
                }
                total_abs_edge_transport_charge_c += std::abs(
                    edge_particle_transport_charge_matrix_c[node_i][node_j]);
                if (node_i < edge_particle_transport_target_charge_matrix_c.size() &&
                    node_j < edge_particle_transport_target_charge_matrix_c[node_i].size())
                {
                    total_abs_edge_target_charge_c += std::abs(
                        edge_particle_transport_target_charge_matrix_c[node_i][node_j]);
                }
                if (node_i < edge_particle_transport_operator_drive_matrix_c.size() &&
                    node_j < edge_particle_transport_operator_drive_matrix_c[node_i].size())
                {
                    total_abs_edge_operator_drive_charge_c += std::abs(
                        edge_particle_transport_operator_drive_matrix_c[node_i][node_j]);
                }
            }
        }
    }

    const double reference_spread_v =
        first_shared_node ? 0.0 : (max_reference_potential_v - min_reference_potential_v);
    const double sheath_charge_spread_c =
        first_shared_node ? 0.0 : (max_sheath_charge_c - min_sheath_charge_c);
    const double patch_potential_spread_v =
        first_shared_node ? 0.0 : (max_patch_potential_v - min_patch_potential_v);
    const double live_pic_electron_current_spread_a_per_m2 =
        first_shared_node
            ? 0.0
            : (max_live_pic_electron_current_a_per_m2 -
               min_live_pic_electron_current_a_per_m2);
    const double live_pic_ion_current_spread_a_per_m2 =
        first_shared_node ? 0.0
                          : (max_live_pic_ion_current_a_per_m2 -
                             min_live_pic_ion_current_a_per_m2);
    const double live_pic_net_current_spread_a_per_m2 =
        first_shared_node ? 0.0
                          : (max_live_pic_net_current_a_per_m2 -
                             min_live_pic_net_current_a_per_m2);
    const double live_pic_collision_count_spread =
        first_shared_node
            ? 0.0
            : (max_live_pic_collision_count - min_live_pic_collision_count);
    const double electron_calibration_factor_spread =
        first_shared_node
            ? 0.0
            : (max_electron_calibration_factor - min_electron_calibration_factor);
    const double ion_calibration_factor_spread =
        first_shared_node ? 0.0
                          : (max_ion_calibration_factor - min_ion_calibration_factor);
    const double normal_electric_field_spread_v_per_m =
        first_shared_node
            ? 0.0
            : (max_normal_electric_field_v_per_m - min_normal_electric_field_v_per_m);
    const double local_charge_density_spread_c_per_m3 =
        first_shared_node
            ? 0.0
            : (max_local_charge_density_c_per_m3 - min_local_charge_density_c_per_m3);
    const double distributed_particle_transport_charge_spread_c =
        first_shared_node
            ? 0.0
            : (max_distributed_particle_transport_charge_c -
               min_distributed_particle_transport_charge_c);
    const double distributed_particle_transport_reference_shift_spread_v =
        first_shared_node
            ? 0.0
            : (max_distributed_particle_transport_reference_shift_v -
               min_distributed_particle_transport_reference_shift_v);
    const double distributed_particle_transport_net_flux_spread_a =
        first_shared_node
            ? 0.0
            : (max_distributed_particle_transport_net_flux_a -
               min_distributed_particle_transport_net_flux_a);
    const double node_edge_transport_charge_spread_c =
        first_shared_node
            ? 0.0
            : (max_node_edge_transport_charge_c - min_node_edge_transport_charge_c);
    const auto latest_scalar = [](const std::vector<double>& values) {
        return values.empty() ? 0.0 : values.back();
    };
    const double pre_global_solve_patch_potential_spread_v =
        latest_scalar(pre_global_solve_patch_spread_history);
    const double patch_potential_spread_reduction_v =
        latest_scalar(patch_spread_reduction_v_history);
    const double patch_potential_spread_reduction_ratio =
        latest_scalar(patch_spread_reduction_ratio_history);
    const bool shared_current_matrix_coupling_active =
        latest_scalar(shared_current_matrix_coupling_active_history) >= 0.5;
    const int shared_current_matrix_coupling_offdiag_entry_count =
        static_cast<int>(std::llround(
            latest_scalar(shared_current_matrix_coupling_offdiag_entry_history)));
    const bool shared_global_coupled_solve_active =
        latest_scalar(shared_global_coupled_solve_active_history) >= 0.5;
    const int shared_global_coupled_solve_iterations =
        static_cast<int>(std::llround(
            latest_scalar(shared_global_coupled_solve_iteration_history)));
    const bool shared_global_coupled_solve_converged =
        latest_scalar(shared_global_coupled_solve_converged_history) >= 0.5;
    const double shared_global_coupled_solve_max_delta_v =
        latest_scalar(shared_global_coupled_solve_max_delta_history);
    const bool shared_live_pic_coupled_refresh_active =
        latest_scalar(shared_live_pic_coupled_refresh_active_history) >= 0.5;
    const int shared_live_pic_coupled_refresh_count =
        static_cast<int>(std::llround(
            latest_scalar(shared_live_pic_coupled_refresh_count_history)));
    const bool shared_particle_transport_coupling_active =
        latest_scalar(shared_particle_transport_coupling_active_history) >= 0.5;
    const int shared_particle_transport_offdiag_entry_count =
        static_cast<int>(std::llround(
            latest_scalar(shared_particle_transport_offdiag_entry_history)));
    const double shared_particle_transport_total_conductance_s =
        latest_scalar(shared_particle_transport_total_conductance_history);
    const double shared_particle_transport_conservation_error_a_per_v =
        latest_scalar(shared_particle_transport_conservation_error_history);
    const double shared_particle_transport_charge_c =
        latest_scalar(shared_particle_transport_charge_history);
    const double shared_particle_transport_reference_shift_v =
        latest_scalar(shared_particle_transport_reference_shift_history);
    const bool shared_particle_transport_distribution_active =
        shared_patch_nodes.size() >= 2 &&
        distributed_particle_transport_charge_spread_c > 0.0;
    const double shared_particle_transport_distribution_conservation_error_c =
        std::abs(sum_distributed_particle_transport_charge_c -
                 shared_particle_transport_charge_c);
    const double shared_particle_transport_exchange_flux_conservation_error_a =
        std::abs(sum_distributed_particle_transport_net_flux_a);
    const bool shared_particle_transport_exchange_active =
        distributed_particle_transport_net_flux_spread_a > 0.0;
    const bool shared_particle_transport_edge_domain_active =
        total_abs_edge_transport_charge_c > 0.0;
    const double shared_particle_transport_edge_charge_conservation_error_c =
        std::abs(sum_node_edge_transport_charge_c);
    const bool shared_particle_transport_edge_operator_active =
        total_abs_edge_operator_drive_charge_c > 0.0;
    const double shared_particle_transport_edge_operator_drive_conservation_error_c =
        std::abs(total_node_edge_operator_drive_charge_c);
    const double shared_global_solve_weight = std::clamp(
        config.material.getScalarProperty("shared_surface_global_solve_weight", 0.0), 0.0, 1.0);
    const bool global_sheath_proxy_solve_active =
        shared_global_solve_weight > 0.0 && shared_patch_nodes.size() >= 2;
    const bool sheath_consistency_pass =
        reference_spread_v <= 1.0e-9 && sheath_charge_spread_c <= 1.0e-15;
    const bool shared_particle_pool_consistency_pass =
        live_pic_electron_current_spread_a_per_m2 <= 1.0e-12 &&
        live_pic_ion_current_spread_a_per_m2 <= 1.0e-12 &&
        live_pic_net_current_spread_a_per_m2 <= 1.0e-12 &&
        live_pic_collision_count_spread <= 0.5 &&
        electron_calibration_factor_spread <= 1.0e-12 &&
        ion_calibration_factor_spread <= 1.0e-12;
    const bool consistency_pass =
        sheath_consistency_pass && shared_particle_pool_consistency_pass;

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
    output << "  \"patch_potential_spread_reduction_v\": "
           << patch_potential_spread_reduction_v << ",\n";
    output << "  \"patch_potential_spread_reduction_ratio\": "
           << patch_potential_spread_reduction_ratio << ",\n";
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
           << (shared_particle_transport_edge_graph_operator_converged ? "true" : "false")
           << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_max_balance_residual_c\": "
           << shared_particle_transport_edge_graph_operator_max_balance_residual_c << ",\n";
    output << "  \"shared_particle_transport_edge_graph_operator_branch_graph_active\": "
           << (shared_particle_transport_edge_graph_operator_branch_graph_edge_count > 0.0
                   ? "true"
                   : "false")
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
    output << "  \"live_pic_collision_count_spread\": "
           << live_pic_collision_count_spread << ",\n";
    output << "  \"electron_calibration_factor_spread\": "
           << electron_calibration_factor_spread << ",\n";
    output << "  \"ion_calibration_factor_spread\": "
           << ion_calibration_factor_spread << ",\n";
    output << "  \"normal_electric_field_spread_v_per_m\": "
           << normal_electric_field_spread_v_per_m << ",\n";
    output << "  \"local_charge_density_spread_c_per_m3\": "
           << local_charge_density_spread_c_per_m3 << ",\n";
    output << "  \"consistency_threshold_reference_v\": 0.000000001000,\n";
    output << "  \"consistency_threshold_sheath_charge_c\": 0.000000000000001000,\n";
    output << "  \"consistency_threshold_live_pic_current_a_per_m2\": 0.000000000001,\n";
    output << "  \"consistency_threshold_live_pic_collision_count\": 0.500000000000,\n";
    output << "  \"consistency_threshold_calibration_factor\": 0.000000000001,\n";
    output << "  \"sheath_consistency_pass\": "
           << (sheath_consistency_pass ? "true" : "false") << ",\n";
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
               << latestHistoryValue(node_ion_calibration_factor_history, node_index)
               << "}";
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
    double edge_graph_operator_iterations,
    bool edge_graph_operator_converged,
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

std::string jsonEscape(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size());
    for (const char character : value)
    {
        switch (character)
        {
        case '\\':
            escaped += "\\\\";
            break;
        case '"':
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
            escaped += character;
            break;
        }
    }
    return escaped;
}

std::string logicalNodeIdFromName(const std::string& node_name)
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

std::string boundaryGroupIdForNode(const SurfaceChargingConfig& config, const std::string& node_name);

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
         surface_index < mapping_state.surface_to_circuit_node_index.size();
         ++surface_index)
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
               << "\", \"boundary_group_id\": \""
               << jsonEscape(boundaryGroupIdForNode(config, node_name))
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
         circuit_index < mapping_state.circuit_to_surface_node_indices.size();
         ++circuit_index)
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
    const SurfaceCircuitModel* circuit_model, const SurfaceGraphCapacitanceMatrixProvider* matrix_provider,
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

std::string boundaryGroupIdForNode(const SurfaceChargingConfig& config, const std::string& node_name)
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

double projectionWeightForNode(const SurfaceChargingConfig& config, const std::string& node_name)
{
    const auto boundary_group_id = boundaryGroupIdForNode(config, node_name);
    if (boundary_group_id.empty())
    {
        return 1.0;
    }

    std::size_t linked_boundary_faces = 0;
    for (const auto& mapping : config.boundary_mappings)
    {
        if (mapping.boundary_group_id != boundary_group_id)
        {
            continue;
        }
        linked_boundary_faces += 1;
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
        return {radius_m * std::cos(angle), radius_m * std::sin(angle), 0.02 * static_cast<double>(node_index)};
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
                   << "\", \"node_type\": \""
                   << (circuit_model->nodeIsPatch(node_index) ? "patch" : "body")
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
        const double face_area_m2 =
            boundaryFaceAreaForGroup(config, circuit_model, group.id);
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
        const double face_area_m2 =
            boundaryFaceAreaForGroup(config, circuit_model, group.id);
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
    output << "  \"surface_pic_strategy\": \""
           << surfacePicStrategyName(config.surface_pic_strategy) << "\",\n";
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
    output << "    \"boundary_mapping_json\": \"" << jsonEscape(artifact(".boundary_mapping.json").string()) << "\",\n";
    output << "    \"field_request_json\": \"" << jsonEscape(artifact(".field_request.json").string()) << "\",\n";
    output << "    \"field_result_template_json\": \"" << jsonEscape(artifact(".field_result_template.json").string()) << "\",\n";
    output << "    \"field_result_json\": \"" << jsonEscape(artifact(".field_result.json").string()) << "\",\n";
    output << "    \"graph_matrix_json\": \"" << jsonEscape(artifact(".graph_matrix.json").string()) << "\",\n";
    output << "    \"field_adapter_json\": \"" << jsonEscape(artifact(".field_adapter.json").string()) << "\",\n";
    output << "    \"volume_stub_json\": \"" << jsonEscape(artifact(".volume_stub.json").string()) << "\",\n";
    output << "    \"volume_mesh_stub_json\": \"" << jsonEscape(artifact(".volume_mesh_stub.json").string()) << "\"\n";
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
                   << node_index << "\", \"boundary_group_id\": \""
                   << jsonEscape(boundaryGroupIdForNode(config, node_name)) << "\", \"weight\": "
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
                   << jsonEscape(node_name) << "\", \"boundary_group_id\": \""
                   << jsonEscape(boundaryGroupIdForNode(config, node_name)) << "\", \"weight\": "
                   << projection_weight << "}";
            output << (node_index + 1 < circuit_model->nodeCount() ? ",\n" : "\n");
        }
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
    output << "  \"solver_mode_hint\": \""
           << volumeLinearSolverPolicyName(config.volume_linear_solver_policy) << "\",\n";
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
    output << "  \"surface_pic_strategy\": \""
           << surfacePicStrategyName(config.surface_pic_strategy) << "\",\n";
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
    output << "  \"solver_mode_hint\": \""
           << volumeLinearSolverPolicyName(config.volume_linear_solver_policy) << "\",\n";
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
    output << "    \"volume_stub_json\": \"" << jsonEscape(artifact(".volume_stub.json").string()) << "\",\n";
    output << "    \"volume_mesh_stub_json\": \"" << jsonEscape(artifact(".volume_mesh_stub.json").string()) << "\",\n";
    output << "    \"surface_volume_projection_json\": \""
           << jsonEscape(artifact(".surface_volume_projection.json").string()) << "\",\n";
    output << "    \"graph_matrix_json\": \"" << jsonEscape(artifact(".graph_matrix.json").string()) << "\",\n";
    output << "    \"field_request_json\": \"" << jsonEscape(artifact(".field_request.json").string()) << "\",\n";
    output << "    \"volume_result_json\": \"" << jsonEscape(artifact(".volume_result.json").string()) << "\"\n";
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
                   << "\", \"node_type\": \""
                   << (circuit_model->nodeIsPatch(node_index) ? "patch" : "body")
                   << "\", \"boundary_group_id\": \""
                   << jsonEscape(boundaryGroupIdForNode(config, node_name))
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

int signOf(double value)
{
    if (value > 0.0)
    {
        return 1;
    }
    if (value < 0.0)
    {
        return -1;
    }
    return 0;
}

}

bool DensePlasmaSurfaceCharging::initialize(const SurfaceChargingConfig& config)
{
    return initialize(compileSurfaceRuntimePlan(config));
}

bool DensePlasmaSurfaceCharging::initialize(const SurfaceRuntimePlan& runtime_plan)
{
    config_ = runtime_plan.compiled_config;
    if (!config_.solver_config.collision_set.empty())
    {
        config_.live_pic_collision_cross_section_set_id = config_.solver_config.collision_set;
    }
    field_volume_adaptive_relaxation_ =
        std::clamp(config_.solver_config.relaxation_factor, 0.1, 1.5);
    field_volume_adaptive_relaxation_initialized_ = false;
    last_error_message_.clear();
    if (!applyImportedSurfaceMaterial(config_, last_error_message_))
    {
        initialized_ = false;
        return false;
    }
    std::string validation_error;
    if (!validateSurfaceChargingConfig(config_, &validation_error))
    {
        initialized_ = false;
        last_error_message_ = validation_error;
        return false;
    }
    orchestrator_ = makeSurfaceScenarioOrchestrator(config_);
    current_model_ = orchestrator_ ? orchestrator_->createCurrentModel() : nullptr;
    voltage_model_ = orchestrator_ ? orchestrator_->createVoltageModel() : nullptr;
    capacitance_model_ = orchestrator_ ? orchestrator_->createCapacitanceModel() : nullptr;
    potential_reference_model_ =
        orchestrator_ ? orchestrator_->createPotentialReferenceModel() : nullptr;
    electric_field_provider_ =
        orchestrator_ ? orchestrator_->createElectricFieldProvider() : nullptr;
    volume_charge_provider_ =
        orchestrator_ ? orchestrator_->createVolumeChargeProvider() : nullptr;
    bubble_capacitance_estimator_ =
        orchestrator_ ? orchestrator_->createBubbleCapacitanceEstimator() : nullptr;
    reference_state_provider_ =
        orchestrator_ ? orchestrator_->createSurfaceReferenceStateProvider() : nullptr;
    reference_graph_propagator_ =
        orchestrator_ ? orchestrator_->createSurfaceReferenceGraphPropagator() : nullptr;
    graph_capacitance_matrix_provider_ =
        orchestrator_ ? orchestrator_->createSurfaceGraphCapacitanceMatrixProvider() : nullptr;
    field_solver_adapter_ =
        orchestrator_ ? orchestrator_->createSurfaceFieldSolverAdapter() : nullptr;
    volumetric_solver_adapter_ =
        orchestrator_ ? orchestrator_->createSurfaceVolumetricSolverAdapter() : nullptr;
    circuit_model_ = orchestrator_ ? orchestrator_->createCircuitModel() : nullptr;
    legacy_benchmark_replay_.clear();
    legacy_benchmark_body_curve_.clear();
    legacy_benchmark_replay_index_ = 0;
    status_ = SurfaceChargingStatus{};
    status_.state.capacitance_per_area_f_per_m2 = computeCapacitancePerArea();
    status_.body_potential_v = config_.body_initial_potential_v;
    status_.patch_potential_v = config_.body_initial_potential_v;
    status_.state.surface_potential_v = status_.patch_potential_v;
    status_.transition_conductivity_scale = 1.0;
    status_.transition_source_flux_scale = 1.0;
    status_.transition_simulation_param_scale = 1.0;
    transition_conductivity_scale_ = 1.0;
    transition_source_flux_scale_ = 1.0;
    transition_simulation_param_scale_ = 1.0;
    transition_runtime_internal_substeps_ = std::max<std::size_t>(1, config_.internal_substeps);
    transition_runtime_max_delta_potential_v_per_step_ =
        std::max(1.0e-3, config_.max_delta_potential_v_per_step);
    transition_runtime_solver_relaxation_factor_ =
        std::clamp(config_.solver_config.relaxation_factor, 0.05, 1.5);
    status_.transition_runtime_internal_substeps =
        transition_runtime_internal_substeps_;
    status_.transition_runtime_max_delta_potential_v_per_step =
        transition_runtime_max_delta_potential_v_per_step_;
    status_.transition_runtime_solver_relaxation_factor =
        transition_runtime_solver_relaxation_factor_;
    initialized_ = true;
    history_time_.clear();
    history_potential_.clear();
    history_charge_.clear();
    history_current_.clear();
    history_electron_current_.clear();
    history_ion_current_.clear();
    history_secondary_current_.clear();
    history_ion_secondary_current_.clear();
    history_backscatter_current_.clear();
    history_photo_current_.clear();
    history_thermionic_current_.clear();
    history_field_emission_current_.clear();
    history_leakage_current_.clear();
    history_ram_current_.clear();
    history_body_potential_.clear();
    history_patch_potential_.clear();
    history_circuit_branch_current_.clear();
    history_current_derivative_.clear();
    history_pic_recalibration_marker_.clear();
    history_live_pic_electron_current_.clear();
    history_live_pic_ion_current_.clear();
    history_live_pic_net_current_.clear();
    history_live_pic_derivative_.clear();
    history_live_pic_collision_count_.clear();
    history_live_pic_mcc_enabled_.clear();
    history_net_current_.clear();
    history_capacitance_.clear();
    history_effective_conductivity_.clear();
    history_effective_sheath_length_.clear();
    history_shared_surface_patch_potential_.clear();
    history_shared_surface_patch_area_.clear();
    history_shared_surface_reference_potential_.clear();
    history_shared_surface_effective_sheath_length_.clear();
    history_shared_surface_sheath_charge_.clear();
    history_shared_surface_runtime_enabled_.clear();
    history_shared_surface_pre_global_solve_patch_potential_spread_.clear();
    history_shared_surface_patch_potential_spread_reduction_v_.clear();
    history_shared_surface_patch_potential_spread_reduction_ratio_.clear();
    history_shared_surface_current_matrix_coupling_active_.clear();
    history_shared_surface_current_matrix_coupling_offdiag_entries_.clear();
    history_shared_surface_global_coupled_solve_active_.clear();
    history_shared_surface_global_coupled_solve_iterations_.clear();
    history_shared_surface_global_coupled_solve_converged_.clear();
    history_shared_surface_global_coupled_solve_max_delta_v_.clear();
    history_shared_surface_live_pic_coupled_refresh_active_.clear();
    history_shared_surface_live_pic_coupled_refresh_count_.clear();
    history_shared_surface_particle_transport_coupling_active_.clear();
    history_shared_surface_particle_transport_offdiag_entries_.clear();
    history_shared_surface_particle_transport_total_conductance_s_.clear();
    history_shared_surface_particle_transport_conservation_error_a_per_v_.clear();
    history_shared_surface_particle_transport_charge_c_.clear();
    history_shared_surface_particle_transport_reference_shift_v_.clear();
    history_normal_electric_field_.clear();
    history_local_charge_density_.clear();
    history_adaptive_time_step_.clear();
    history_internal_substeps_.clear();
    history_electron_calibration_factor_.clear();
    history_ion_calibration_factor_.clear();
    history_equilibrium_potential_.clear();
    history_equilibrium_error_.clear();
    history_surface_node_potentials_.clear();
    history_surface_node_total_currents_.clear();
    history_surface_node_electron_currents_.clear();
    history_surface_node_ion_currents_.clear();
    history_surface_node_secondary_currents_.clear();
    history_surface_node_ion_secondary_currents_.clear();
    history_surface_node_backscatter_currents_.clear();
    history_surface_node_photo_currents_.clear();
    history_surface_node_thermionic_currents_.clear();
    history_surface_node_field_emission_currents_.clear();
    history_surface_node_conduction_currents_.clear();
    history_surface_node_ram_currents_.clear();
    history_surface_node_current_derivatives_.clear();
    history_surface_node_effective_sheath_lengths_.clear();
    history_surface_node_normal_electric_fields_.clear();
    history_surface_node_local_charge_densities_.clear();
    history_surface_node_propagated_reference_potentials_.clear();
    history_surface_node_field_solver_reference_potentials_.clear();
    history_surface_node_shared_runtime_enabled_.clear();
    history_surface_node_shared_patch_potentials_.clear();
    history_surface_node_shared_patch_areas_.clear();
    history_surface_node_shared_reference_potentials_.clear();
    history_surface_node_shared_sheath_charges_.clear();
    history_surface_node_distributed_particle_transport_charges_.clear();
    history_surface_node_distributed_particle_transport_reference_shifts_.clear();
    history_surface_node_distributed_particle_transport_net_fluxes_.clear();
    history_surface_node_graph_capacitance_diagonals_.clear();
    history_surface_node_graph_capacitance_row_sums_.clear();
    history_surface_node_field_solver_coupling_gains_.clear();
    history_surface_node_field_solver_capacitance_scales_.clear();
    history_surface_node_pseudo_volumes_.clear();
    history_surface_node_volume_projection_weight_sums_.clear();
    history_surface_node_volume_mesh_coupling_gains_.clear();
    history_surface_node_volume_potentials_.clear();
    history_surface_node_deposited_charges_.clear();
    history_surface_node_poisson_residuals_.clear();
    history_surface_node_volume_solver_mode_ids_.clear();
    history_surface_node_volume_solver_iterations_.clear();
    history_surface_node_volume_solver_linear_iterations_.clear();
    history_surface_node_volume_solver_converged_.clear();
    history_surface_node_volume_solver_residual_norms_.clear();
    history_surface_node_volume_solver_max_deltas_.clear();
    history_surface_node_volume_solver_matrix_nnzs_.clear();
    history_surface_node_volume_solver_cell_counts_.clear();
    history_surface_node_field_volume_coupling_iterations_.clear();
    history_surface_node_field_volume_coupling_converged_.clear();
    history_surface_node_field_volume_coupling_max_deltas_.clear();
    history_surface_node_field_volume_coupling_relaxation_used_.clear();
    history_surface_node_external_volume_feedback_blend_factors_.clear();
    history_surface_node_external_volume_feedback_mismatch_metrics_.clear();
    history_surface_node_external_volume_feedback_applied_.clear();
    history_surface_node_electron_calibration_factors_.clear();
    history_surface_node_ion_calibration_factors_.clear();
    history_surface_node_live_pic_electron_currents_.clear();
    history_surface_node_live_pic_ion_currents_.clear();
    history_surface_node_live_pic_net_currents_.clear();
    history_surface_node_live_pic_derivatives_.clear();
    history_surface_node_live_pic_collision_counts_.clear();
    history_surface_node_live_pic_mcc_enabled_.clear();
    history_surface_branch_currents_.clear();
    history_surface_branch_conductances_.clear();
    history_surface_branch_voltage_drops_.clear();
    history_surface_branch_power_w_.clear();
    history_surface_branch_mutual_capacitances_.clear();
    circuit_node_potentials_.clear();
    global_particle_domain_state_ = GlobalParticleDomainState{};
    global_particle_domain_state_.bookkeeping_mode = "owned_global_particle_domain_state_v2";
    global_particle_repository_state_ = GlobalParticleRepositoryState{};
    global_particle_repository_state_.bookkeeping_mode =
        "owned_global_particle_repository_state_v1";
    global_particle_repository_state_.lifecycle_mode =
        "global_particle_transport_reservoir_lifecycle_v1";
    global_sheath_field_solve_state_ = GlobalSheathFieldSolveState{};
    global_sheath_field_solve_state_.solve_mode = "explicit_global_sheath_field_linear_system_v2";
    shared_particle_transport_charge_c_ = 0.0;
    shared_particle_transport_reference_shift_v_ = 0.0;
    shared_particle_transport_node_charge_c_.clear();
    shared_particle_transport_node_net_flux_a_.clear();
    shared_particle_transport_edge_charge_matrix_c_.clear();
    shared_particle_transport_edge_target_charge_matrix_c_.clear();
    shared_particle_transport_edge_operator_drive_matrix_c_.clear();
    shared_particle_transport_exchange_flux_matrix_a_.clear();
    global_particle_transport_conductance_matrix_s_.clear();
    global_sheath_field_reference_response_matrix_.clear();
    global_sheath_field_previous_node_charge_c_.clear();
    shared_particle_transport_edge_graph_operator_iterations_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_converged_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_effective_pair_count_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_ = 0.0;
    shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_ = 1.0;
    shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_ = 1.0;
    const double effective_conductivity = computeEffectiveConductivity(status_.patch_potential_v);
    status_.state.surface_charge_density_c_per_m2 =
        status_.state.capacitance_per_area_f_per_m2 *
        (status_.patch_potential_v - status_.body_potential_v);
    status_.state.stored_energy_j_per_m2 =
        0.5 * status_.state.capacitance_per_area_f_per_m2 *
        (status_.patch_potential_v - status_.body_potential_v) *
        (status_.patch_potential_v - status_.body_potential_v);
    rebuildSurfaceCircuit(effective_conductivity);
    if (circuit_model_)
    {
        circuit_node_potentials_.resize(circuit_model_->nodeCount(), 0.0);
        history_surface_node_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_total_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_electron_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_ion_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_secondary_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_ion_secondary_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_backscatter_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_photo_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_thermionic_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_emission_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_conduction_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_ram_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_current_derivatives_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_effective_sheath_lengths_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_normal_electric_fields_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_local_charge_densities_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_propagated_reference_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_solver_reference_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_shared_runtime_enabled_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_shared_patch_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_shared_patch_areas_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_shared_reference_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_shared_sheath_charges_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_distributed_particle_transport_charges_.assign(
            circuit_model_->nodeCount(), {});
        history_surface_node_distributed_particle_transport_reference_shifts_.assign(
            circuit_model_->nodeCount(), {});
        history_surface_node_distributed_particle_transport_net_fluxes_.assign(
            circuit_model_->nodeCount(), {});
        history_surface_node_graph_capacitance_diagonals_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_graph_capacitance_row_sums_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_solver_coupling_gains_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_solver_capacitance_scales_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_pseudo_volumes_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_projection_weight_sums_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_mesh_coupling_gains_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_deposited_charges_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_poisson_residuals_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_mode_ids_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_iterations_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_linear_iterations_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_converged_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_residual_norms_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_max_deltas_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_matrix_nnzs_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_volume_solver_cell_counts_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_volume_coupling_iterations_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_volume_coupling_converged_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_volume_coupling_max_deltas_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_volume_coupling_relaxation_used_.assign(circuit_model_->nodeCount(),
                                                                            {});
        history_surface_node_external_volume_feedback_blend_factors_.assign(
            circuit_model_->nodeCount(), {});
        history_surface_node_external_volume_feedback_mismatch_metrics_.assign(
            circuit_model_->nodeCount(), {});
        history_surface_node_external_volume_feedback_applied_.assign(circuit_model_->nodeCount(),
                                                                      {});
        history_surface_node_electron_calibration_factors_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_ion_calibration_factors_.assign(circuit_model_->nodeCount(), {});
        shared_particle_transport_node_charge_c_.assign(circuit_model_->nodeCount(), 0.0);
        shared_particle_transport_node_net_flux_a_.assign(circuit_model_->nodeCount(), 0.0);
        shared_particle_transport_edge_charge_matrix_c_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        shared_particle_transport_edge_target_charge_matrix_c_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        shared_particle_transport_edge_operator_drive_matrix_c_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        shared_particle_transport_exchange_flux_matrix_a_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        global_particle_transport_conductance_matrix_s_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        global_sheath_field_reference_response_matrix_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        global_sheath_field_previous_node_charge_c_.assign(circuit_model_->nodeCount(), 0.0);
        history_surface_node_live_pic_electron_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_live_pic_ion_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_live_pic_net_currents_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_live_pic_derivatives_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_live_pic_collision_counts_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_live_pic_mcc_enabled_.assign(circuit_model_->nodeCount(), {});
        history_surface_branch_currents_.assign(circuit_model_->branchCount(), {});
        history_surface_branch_conductances_.assign(circuit_model_->branchCount(), {});
        history_surface_branch_voltage_drops_.assign(circuit_model_->branchCount(), {});
        history_surface_branch_power_w_.assign(circuit_model_->branchCount(), {});
        history_surface_branch_mutual_capacitances_.assign(circuit_model_->branchCount(), {});
        for (std::size_t node_index = 0; node_index < circuit_model_->nodeCount(); ++node_index)
        {
            circuit_node_potentials_[node_index] = circuit_model_->nodePotential(node_index);
        }
    }
    applyGeoPicCalibration();
    global_particle_domain_state_ = GlobalParticleDomainState{};
    global_particle_domain_state_.bookkeeping_mode = "owned_global_particle_domain_state_v2";
    global_particle_repository_state_ = GlobalParticleRepositoryState{};
    global_particle_repository_state_.bookkeeping_mode =
        "owned_global_particle_repository_state_v1";
    global_particle_repository_state_.lifecycle_mode =
        "global_particle_transport_reservoir_lifecycle_v1";
    global_sheath_field_solve_state_ = GlobalSheathFieldSolveState{};
    global_sheath_field_solve_state_.solve_mode = "explicit_global_sheath_field_linear_system_v2";
    field_volume_adaptive_relaxation_ =
        std::clamp(config_.field_volume_outer_relaxation, 0.05, 1.0);
    field_volume_adaptive_relaxation_initialized_ = true;
    loadLegacyBenchmarkReplay();
    if (legacy_benchmark_replay_.empty())
    {
        loadLegacyBenchmarkExecutionCurve();
    }
    return true;
}

ReferenceCurrentBalanceConfig
DensePlasmaSurfaceCharging::buildReferenceConfig(double conductivity_s_per_m) const
{
    ReferenceCurrentBalanceConfig reference_config;
    const double source_flux_scale = computeTransitionSourceFluxScale();
    reference_config.plasma = config_.plasma;
    reference_config.plasma.electron_density_m3 =
        std::max(0.0, reference_config.plasma.electron_density_m3 * source_flux_scale);
    reference_config.plasma.ion_density_m3 =
        std::max(0.0, reference_config.plasma.ion_density_m3 * source_flux_scale);
    reference_config.electron_spectrum = config_.electron_spectrum;
    reference_config.ion_spectrum = config_.ion_spectrum;
    reference_config.has_electron_spectrum = config_.has_electron_spectrum;
    reference_config.has_ion_spectrum = config_.has_ion_spectrum;
    reference_config.patch_material = config_.material;
    reference_config.see_model = config_.reference_see_model;
    reference_config.patch_incidence_angle_deg = config_.patch_incidence_angle_deg;
    reference_config.patch_flow_angle_deg = config_.patch_flow_angle_deg;
    reference_config.patch_thickness_m = config_.dielectric_thickness_m;
    reference_config.patch_conductivity_s_per_m = conductivity_s_per_m;
    reference_config.electron_collection_model = config_.electron_collection_model;
    reference_config.electron_collection_coefficient = config_.electron_collection_coefficient;
    reference_config.ion_collection_coefficient = config_.ion_collection_coefficient;
    reference_config.bulk_flow_velocity_m_per_s = config_.bulk_flow_velocity_m_per_s;
    reference_config.flow_alignment_cosine = config_.flow_alignment_cosine;
    reference_config.ion_directed_velocity_m_per_s = config_.ion_directed_velocity_m_per_s;
    reference_config.electron_calibration_factor =
        current_model_ ? current_model_->electronCalibrationFactor() : 1.0;
    reference_config.ion_calibration_factor =
        current_model_ ? current_model_->ionCalibrationFactor() : 1.0;
    reference_config.photoelectron_temperature_ev =
        std::max(1.0e-3, config_.photoelectron_temperature_ev);
    reference_config.enable_ram_current =
        config_.regime == SurfaceChargingRegime::LeoFlowingPlasma;
    reference_config.enable_secondary_electron = config_.enable_secondary_electron;
    reference_config.enable_backscatter = config_.enable_backscatter;
    reference_config.enable_photoelectron = config_.enable_photoelectron;

    const double effective_photon_flux_m2_s =
        std::max(0.0, config_.emission.photon_flux_m2_s * source_flux_scale);
    const double derived_photo_current =
        kElementaryCharge * effective_photon_flux_m2_s *
        std::clamp(config_.material.getScalarProperty("photoelectron_yield", 0.0) *
                       config_.emission.enhancement_factor,
                   0.0, 1.0);
    const double body_photo_current_density =
        config_.body_photo_current_density_a_per_m2 > 0.0
            ? config_.body_photo_current_density_a_per_m2 * source_flux_scale
            : derived_photo_current;
    const double patch_photo_current_density =
        config_.patch_photo_current_density_a_per_m2 > 0.0
            ? config_.patch_photo_current_density_a_per_m2 * source_flux_scale
            : derived_photo_current;
    reference_config.body_photo_current_density_a_per_m2 =
        body_photo_current_density;
    reference_config.patch_photo_current_density_a_per_m2 =
        patch_photo_current_density;

    reference_config.body_material = config_.material;
    reference_config.body_material.setType(Mesh::MaterialType::CONDUCTOR);
    reference_config.body_material.setName(config_.material.getName() + "_body");
    reference_config.body_material.setConductivity(
        std::max(1.0e-6, config_.material.getScalarProperty("body_conductivity_s_per_m", 1.0e4)));
    applyBodyMaterialOverrides(config_.material, reference_config.body_material);
    return reference_config;
}

SurfaceModelRuntimeState
DensePlasmaSurfaceCharging::buildRuntimeState(double body_potential_v, double patch_potential_v,
                                              double surface_area_m2,
                                              double effective_conductivity_s_per_m,
                                              double effective_sheath_length_m,
                                              std::size_t node_index,
                                              const std::string& node_name) const
{
    SurfaceModelRuntimeState state;
    state.body_potential_v = body_potential_v;
    state.patch_potential_v = patch_potential_v;
    state.surface_area_m2 = surface_area_m2;
    state.effective_conductivity_s_per_m = effective_conductivity_s_per_m;
    state.effective_sheath_length_m = effective_sheath_length_m;
    state.electron_calibration_factor =
        current_model_ ? current_model_->electronCalibrationFactor() : 1.0;
    state.ion_calibration_factor = current_model_ ? current_model_->ionCalibrationFactor() : 1.0;
    state.node_index = node_index;
    state.node_name = node_name;
    const bool patch_runtime_state =
        circuit_model_ ? circuit_model_->nodeIsPatch(node_index)
                       : (node_index != 0 || node_name == "patch");
    const bool shared_patch_runtime = useSharedSurfacePicRuntime() && patch_runtime_state;
    state.shared_runtime_enabled = shared_patch_runtime;
    state.shared_patch_potential_v = patch_potential_v;
    state.shared_patch_area_m2 = surface_area_m2;
    state.shared_reference_potential_v = patch_potential_v;
    state.shared_effective_sheath_length_m = effective_sheath_length_m;
    state.shared_particle_transport_charge_c = currentSharedSurfaceParticleTransportChargeC();
    state.shared_particle_transport_reference_shift_v =
        currentSharedSurfaceParticleTransportReferenceShiftV();
    state.shared_particle_transport_domain_active = useSharedSurfaceParticleTransportCoupling();
    state.global_particle_domain_active =
        global_particle_domain_state_.active ||
        (patch_runtime_state && useSharedSurfaceParticleTransportCoupling() && circuit_model_ &&
         circuit_model_->patchCount() >= 2);
    state.distributed_particle_transport_charge_c =
        sharedSurfaceDistributedParticleTransportChargeC(node_index);
    state.distributed_particle_transport_reference_shift_v = 0.0;
    state.distributed_particle_transport_net_flux_a = 0.0;
    state.distributed_particle_transport_active =
        patch_runtime_state && useSharedSurfaceParticleTransportCoupling();
    if (shared_patch_runtime)
    {
        state.shared_patch_potential_v = computeSharedSurfacePatchPotentialV();
        state.shared_patch_area_m2 = std::max(1.0e-16, computeSharedSurfacePatchAreaM2());
        state.shared_effective_sheath_length_m =
            computeSharedSurfaceEffectiveSheathLengthM(effective_sheath_length_m);
        state.shared_particle_transport_reference_shift_v =
            computeSharedSurfaceParticleTransportReferenceShiftV(
                state.shared_patch_area_m2, state.shared_effective_sheath_length_m);
        state.distributed_particle_transport_reference_shift_v =
            computeSharedSurfaceDistributedParticleTransportReferenceShiftV(
                node_index, surface_area_m2, state.shared_effective_sheath_length_m);
        if (node_index < shared_particle_transport_node_net_flux_a_.size())
        {
            state.distributed_particle_transport_net_flux_a =
                shared_particle_transport_node_net_flux_a_[node_index];
        }
        if (const auto* global_particle_node = findGlobalParticleDomainNodeState(node_index))
        {
            state.distributed_particle_transport_charge_c = global_particle_node->charge_c;
            state.distributed_particle_transport_reference_shift_v =
                global_particle_node->reference_shift_v;
            state.distributed_particle_transport_net_flux_a = global_particle_node->net_flux_a;
            state.shared_reference_potential_v =
                global_particle_node->shared_reference_potential_v;
            state.global_particle_domain_node_charge_c = global_particle_node->charge_c;
            state.global_particle_domain_node_flux_a = global_particle_node->net_flux_a;
        }
        state.global_particle_domain_charge_c =
            currentSharedSurfaceParticleTransportChargeC();
        state.global_particle_domain_charge_conservation_error_c =
            global_particle_domain_state_.active
                ? global_particle_domain_state_.charge_conservation_error_c
                : computeSharedSurfaceGlobalParticleDomainChargeConservationErrorC();
        state.global_particle_domain_flux_conservation_error_a =
            global_particle_domain_state_.active
                ? global_particle_domain_state_.flux_conservation_error_a
                : computeSharedSurfaceGlobalParticleDomainFluxConservationErrorA();
        state.surface_area_m2 = state.shared_patch_area_m2;
        state.effective_sheath_length_m = state.shared_effective_sheath_length_m;
        state.propagated_reference_potential_v = state.shared_patch_potential_v;
    }
    if (current_model_)
    {
        state.live_pic_sample = &current_model_->latestLivePicSample();
    }
    if (reference_state_provider_)
    {
        state.reference_plasma_potential_v =
            reference_state_provider_->plasmaReferencePotentialV(config_, state);
    }
    else if (potential_reference_model_)
    {
        state.reference_plasma_potential_v =
            potential_reference_model_->plasmaReferencePotentialV(config_, state);
    }
    else
    {
        state.reference_plasma_potential_v = config_.live_pic_reference_potential_v;
    }
    if (shared_patch_runtime)
    {
        state.shared_reference_potential_v =
            computeSharedSurfaceReferencePotentialV(state.reference_plasma_potential_v,
                                                    patch_potential_v);
        state.reference_plasma_potential_v = state.shared_reference_potential_v;
    }
    state.propagated_reference_potential_v = state.reference_plasma_potential_v;
    if (reference_graph_propagator_)
    {
        state.propagated_reference_potential_v =
            reference_graph_propagator_->propagatedReferencePotentialV(
                config_, circuit_model_.get(), state, state.reference_plasma_potential_v);
    }
    if (std::isfinite(state.propagated_reference_potential_v))
    {
        state.reference_plasma_potential_v = state.propagated_reference_potential_v;
    }
    if (shared_patch_runtime)
    {
        const double particle_transport_reference_shift_v =
            state.shared_particle_transport_reference_shift_v;
        const double post_graph_blend = std::clamp(
            config_.material.getScalarProperty("shared_surface_reference_post_graph_blend", 0.60),
            0.0, 1.0);
        const double blended_shared_reference_potential_v =
            (1.0 - post_graph_blend) * state.reference_plasma_potential_v +
            post_graph_blend * computeSharedSurfaceReferencePotentialV(
                                  state.reference_plasma_potential_v, patch_potential_v);
        state.shared_reference_potential_v =
            blended_shared_reference_potential_v + particle_transport_reference_shift_v;
        if (const auto* global_particle_node = findGlobalParticleDomainNodeState(node_index))
        {
            state.shared_reference_potential_v =
                global_particle_node->shared_reference_potential_v;
        }
        state.propagated_reference_potential_v = state.shared_reference_potential_v;
        state.reference_plasma_potential_v = state.shared_reference_potential_v;
        state.global_sheath_field_reference_potential_v = state.shared_reference_potential_v;
    }
    if (graph_capacitance_matrix_provider_)
    {
        state.graph_capacitance_diagonal_f = graph_capacitance_matrix_provider_->diagonalCapacitanceF(
            config_, circuit_model_.get(), state.node_index);
        state.graph_capacitance_row_sum_f = graph_capacitance_matrix_provider_->rowSumCapacitanceF(
            config_, circuit_model_.get(), state.node_index);
    }
    if (volume_charge_provider_)
    {
        state.local_charge_density_c_per_m3 =
            volume_charge_provider_->localChargeDensityCPerM3(config_, state);
    }
    if (electric_field_provider_)
    {
        state.normal_electric_field_v_per_m =
            electric_field_provider_->normalFieldVPerM(config_, state);
    }
    state.field_solver_reference_potential_v = state.reference_plasma_potential_v;
    const bool has_coupled_adapters = field_solver_adapter_ || volumetric_solver_adapter_;
    if (has_coupled_adapters)
    {
        const bool has_external_field_result =
            config_.enable_external_field_solver_bridge &&
            !config_.external_field_solver_result_path.empty() &&
            std::filesystem::exists(config_.external_field_solver_result_path);
        const bool has_external_volume_result =
            config_.enable_external_volume_solver_bridge &&
            !config_.external_volume_solver_result_path.empty() &&
            std::filesystem::exists(config_.external_volume_solver_result_path);
        const int outer_iterations =
            static_cast<int>(std::clamp<std::size_t>(config_.field_volume_outer_iterations, 1, 8));
        const double outer_tolerance =
            std::clamp(config_.field_volume_outer_tolerance, 1.0e-9, 10.0);
        if (!field_volume_adaptive_relaxation_initialized_)
        {
            field_volume_adaptive_relaxation_ =
                std::clamp(config_.field_volume_outer_relaxation, 0.05, 1.0);
            field_volume_adaptive_relaxation_initialized_ = true;
        }
        const bool external_bridge_authoritative =
            has_external_field_result || has_external_volume_result;
        const double outer_relaxation =
            external_bridge_authoritative
                ? 1.0
                : std::clamp(field_volume_adaptive_relaxation_, 0.05, 1.0);
        state.field_volume_coupling_relaxation_used = outer_relaxation;
        double max_outer_delta = 0.0;
        bool outer_converged = false;
        for (int iteration = 0; iteration < outer_iterations; ++iteration)
        {
            const double previous_reference_v = state.field_solver_reference_potential_v;
            const double previous_field_v_per_m = state.normal_electric_field_v_per_m;
            const double previous_charge_density = state.local_charge_density_c_per_m3;
            const double previous_capacitance_scale = state.field_solver_capacitance_scale;
            const double previous_volume_potential_v = state.volume_potential_v;
            const double previous_volume_coupling_gain = state.volume_mesh_coupling_gain;

            SurfaceModelRuntimeState candidate = state;
            if (field_solver_adapter_)
            {
                field_solver_adapter_->adaptFieldState(config_, circuit_model_.get(), candidate);
                if (std::isfinite(candidate.field_solver_reference_potential_v))
                {
                    candidate.reference_plasma_potential_v =
                        candidate.field_solver_reference_potential_v;
                }
            }
            if (volumetric_solver_adapter_)
            {
                volumetric_solver_adapter_->adaptVolumeState(config_, circuit_model_.get(), candidate);
                if (std::isfinite(candidate.field_solver_reference_potential_v))
                {
                    candidate.reference_plasma_potential_v =
                        candidate.field_solver_reference_potential_v;
                }
            }

            auto relax_value = [outer_relaxation](double previous, double updated) {
                return (1.0 - outer_relaxation) * previous + outer_relaxation * updated;
            };
            state.field_solver_reference_potential_v =
                relax_value(previous_reference_v, candidate.field_solver_reference_potential_v);
            state.reference_plasma_potential_v = state.field_solver_reference_potential_v;
            state.normal_electric_field_v_per_m =
                relax_value(previous_field_v_per_m, candidate.normal_electric_field_v_per_m);
            state.local_charge_density_c_per_m3 =
                relax_value(previous_charge_density, candidate.local_charge_density_c_per_m3);
            state.field_solver_capacitance_scale = relax_value(previous_capacitance_scale,
                                                               candidate.field_solver_capacitance_scale);
            state.volume_potential_v =
                relax_value(previous_volume_potential_v, candidate.volume_potential_v);
            state.volume_mesh_coupling_gain =
                relax_value(previous_volume_coupling_gain, candidate.volume_mesh_coupling_gain);
            state.pseudo_volume_m3 = candidate.pseudo_volume_m3;
            state.volume_projection_weight_sum = candidate.volume_projection_weight_sum;
            state.deposited_charge_c = candidate.deposited_charge_c;
            state.poisson_residual_v_m = candidate.poisson_residual_v_m;
            state.volume_solver_mode_id = candidate.volume_solver_mode_id;
            state.volume_solver_iterations = candidate.volume_solver_iterations;
            state.volume_solver_linear_iterations = candidate.volume_solver_linear_iterations;
            state.volume_solver_converged = candidate.volume_solver_converged;
            state.volume_solver_residual_norm = candidate.volume_solver_residual_norm;
            state.volume_solver_max_delta_v = candidate.volume_solver_max_delta_v;
            state.volume_solver_matrix_nnz = candidate.volume_solver_matrix_nnz;
            state.volume_solver_cell_count = candidate.volume_solver_cell_count;
            state.field_solver_coupling_gain = candidate.field_solver_coupling_gain;
            state.effective_sheath_length_m =
                relax_value(state.effective_sheath_length_m, candidate.effective_sheath_length_m);
            state.external_volume_feedback_blend_factor =
                candidate.external_volume_feedback_blend_factor;
            state.external_volume_feedback_mismatch_metric =
                candidate.external_volume_feedback_mismatch_metric;
            state.external_volume_feedback_applied =
                candidate.external_volume_feedback_applied;

            max_outer_delta = std::max(
                {std::abs(state.field_solver_reference_potential_v - previous_reference_v),
                 std::abs(state.normal_electric_field_v_per_m - previous_field_v_per_m),
                 std::abs(state.local_charge_density_c_per_m3 - previous_charge_density),
                 std::abs(state.field_solver_capacitance_scale - previous_capacitance_scale),
                 std::abs(state.volume_potential_v - previous_volume_potential_v),
                 std::abs(state.volume_mesh_coupling_gain - previous_volume_coupling_gain)});
            state.field_volume_coupling_iterations = static_cast<double>(iteration + 1);
            state.field_volume_coupling_max_delta = max_outer_delta;
            if (max_outer_delta <= outer_tolerance)
            {
                outer_converged = true;
                break;
            }
        }
        state.field_volume_coupling_converged = outer_converged ? 1.0 : 0.0;
        if (!external_bridge_authoritative)
        {
            const double configured_relaxation =
                std::clamp(config_.field_volume_outer_relaxation, 0.05, 1.0);
            double adaptive_relaxation = std::clamp(field_volume_adaptive_relaxation_, 0.05, 1.0);
            if (outer_converged)
            {
                if (max_outer_delta < 0.25 * outer_tolerance)
                {
                    adaptive_relaxation = std::min(1.0, adaptive_relaxation + 0.03);
                }
                else if (max_outer_delta > 2.0 * outer_tolerance)
                {
                    adaptive_relaxation = std::max(0.05, adaptive_relaxation - 0.05);
                }
                else
                {
                    adaptive_relaxation +=
                        0.15 * (configured_relaxation - adaptive_relaxation);
                }
            }
            else
            {
                adaptive_relaxation = std::max(0.05, adaptive_relaxation - 0.08);
            }
            field_volume_adaptive_relaxation_ = std::clamp(adaptive_relaxation, 0.05, 1.0);
        }
    }
    if (!std::isfinite(state.normal_electric_field_v_per_m))
    {
        const double reference_drop_v =
            std::abs(state.reference_plasma_potential_v - state.patch_potential_v);
        state.normal_electric_field_v_per_m =
            reference_drop_v / std::max(1.0e-6, state.effective_sheath_length_m);
    }
    if (!std::isfinite(state.local_charge_density_c_per_m3))
    {
        state.local_charge_density_c_per_m3 = 0.0;
    }
    if (!std::isfinite(state.graph_capacitance_row_sum_f))
    {
        state.graph_capacitance_row_sum_f = 0.0;
    }
    if (!std::isfinite(state.field_solver_reference_potential_v))
    {
        state.field_solver_reference_potential_v = state.reference_plasma_potential_v;
    }
    if (!std::isfinite(state.field_solver_coupling_gain))
    {
        state.field_solver_coupling_gain = 0.0;
    }
    if (!std::isfinite(state.field_solver_capacitance_scale))
    {
        state.field_solver_capacitance_scale = 1.0;
    }
    if (!std::isfinite(state.pseudo_volume_m3))
    {
        state.pseudo_volume_m3 = 0.0;
    }
    if (!std::isfinite(state.volume_projection_weight_sum))
    {
        state.volume_projection_weight_sum = 1.0;
    }
    if (!std::isfinite(state.volume_mesh_coupling_gain))
    {
        state.volume_mesh_coupling_gain = 0.0;
    }
    if (!std::isfinite(state.volume_potential_v))
    {
        state.volume_potential_v = state.field_solver_reference_potential_v;
    }
    if (shared_patch_runtime)
    {
        state.shared_effective_sheath_length_m =
            computeSharedSurfaceEffectiveSheathLengthM(state.effective_sheath_length_m);
        const double sheath_blend = std::clamp(
            config_.material.getScalarProperty("shared_surface_sheath_coupling_weight", 0.70),
            0.0, 1.0);
        state.effective_sheath_length_m =
            (1.0 - sheath_blend) * state.effective_sheath_length_m +
            sheath_blend * state.shared_effective_sheath_length_m;

        const double shared_drop_v =
            state.shared_reference_potential_v - state.shared_patch_potential_v;
        const double shared_normal_electric_field_v_per_m =
            shared_drop_v / std::max(1.0e-6, state.shared_effective_sheath_length_m);
        const double shared_local_charge_density_c_per_m3 =
            kEpsilon0 * shared_normal_electric_field_v_per_m /
            std::max(1.0e-6, state.shared_effective_sheath_length_m);
        const bool shared_global_sheath_field_active =
            useSharedSurfaceGlobalCoupledSolve() && circuit_model_ && circuit_model_->patchCount() >= 2;
        const bool assembled_global_sheath_field_active = global_sheath_field_solve_state_.active;
        state.global_sheath_field_solve_active =
            assembled_global_sheath_field_active || shared_global_sheath_field_active;
        if (assembled_global_sheath_field_active)
        {
            state.global_sheath_field_reference_potential_v =
                global_sheath_field_solve_state_.global_reference_potential_v;
            if (const auto* sheath_node = findGlobalSheathFieldNodeState(node_index))
            {
                state.field_solver_reference_potential_v = sheath_node->reference_potential_v;
                state.shared_reference_potential_v = sheath_node->reference_potential_v;
                state.normal_electric_field_v_per_m = sheath_node->normal_electric_field_v_per_m;
                state.local_charge_density_c_per_m3 = sheath_node->local_charge_density_c_per_m3;
            }
            else
            {
                state.field_solver_reference_potential_v =
                    global_sheath_field_solve_state_.global_reference_potential_v;
                state.shared_reference_potential_v =
                    global_sheath_field_solve_state_.global_reference_potential_v;
                state.normal_electric_field_v_per_m =
                    global_sheath_field_solve_state_.global_normal_electric_field_v_per_m;
                state.local_charge_density_c_per_m3 =
                    global_sheath_field_solve_state_.global_local_charge_density_c_per_m3;
            }
            state.reference_plasma_potential_v =
                state.field_solver_reference_potential_v +
                state.distributed_particle_transport_reference_shift_v;
            state.propagated_reference_potential_v = state.reference_plasma_potential_v;
        }
        else if (shared_global_sheath_field_active)
        {
            state.field_solver_reference_potential_v = state.shared_reference_potential_v;
            state.reference_plasma_potential_v =
                state.shared_reference_potential_v +
                state.distributed_particle_transport_reference_shift_v;
            state.propagated_reference_potential_v = state.reference_plasma_potential_v;
            state.normal_electric_field_v_per_m = shared_normal_electric_field_v_per_m;
            state.local_charge_density_c_per_m3 = shared_local_charge_density_c_per_m3;
            state.global_sheath_field_reference_potential_v = state.shared_reference_potential_v;
        }
        else
        {
            const double field_blend = std::clamp(
                config_.material.getScalarProperty("shared_surface_field_coupling_weight", 0.55),
                0.0, 1.0);
            const double charge_blend = std::clamp(
                config_.material.getScalarProperty("shared_surface_charge_coupling_weight", 0.45),
                0.0, 1.0);
            state.normal_electric_field_v_per_m =
                (1.0 - field_blend) * state.normal_electric_field_v_per_m +
                field_blend * shared_normal_electric_field_v_per_m;
            state.local_charge_density_c_per_m3 =
                (1.0 - charge_blend) * state.local_charge_density_c_per_m3 +
                charge_blend * shared_local_charge_density_c_per_m3;
            state.global_sheath_field_reference_potential_v = state.shared_reference_potential_v;
        }
        state.global_sheath_field_residual_v_per_m =
            state.global_sheath_field_solve_active
                ? global_sheath_field_solve_state_.field_residual_v_per_m
                : computeSharedSurfaceGlobalSheathFieldResidualVPerM(state);
        state.global_particle_field_coupled_residual_v =
            state.global_sheath_field_solve_active
                ? global_sheath_field_solve_state_.particle_field_coupled_residual_v
                : computeSharedSurfaceGlobalParticleFieldCoupledResidualV(state);
        state.global_sheath_field_multi_step_stability_metric_v =
            state.global_sheath_field_solve_active
                ? global_sheath_field_solve_state_.multi_step_stability_metric_v
                : computeSharedSurfaceGlobalSheathFieldMultiStepStabilityMetricV();
    }
    if (!std::isfinite(state.deposited_charge_c))
    {
        state.deposited_charge_c = 0.0;
    }
    if (!std::isfinite(state.poisson_residual_v_m))
    {
        state.poisson_residual_v_m = 0.0;
    }
    if (!std::isfinite(state.volume_solver_mode_id))
    {
        state.volume_solver_mode_id = 0.0;
    }
    if (!std::isfinite(state.volume_solver_iterations))
    {
        state.volume_solver_iterations = 0.0;
    }
    if (!std::isfinite(state.volume_solver_linear_iterations))
    {
        state.volume_solver_linear_iterations = 0.0;
    }
    if (!std::isfinite(state.volume_solver_converged))
    {
        state.volume_solver_converged = 0.0;
    }
    if (!std::isfinite(state.volume_solver_residual_norm))
    {
        state.volume_solver_residual_norm = 0.0;
    }
    if (!std::isfinite(state.volume_solver_max_delta_v))
    {
        state.volume_solver_max_delta_v = 0.0;
    }
    if (!std::isfinite(state.volume_solver_matrix_nnz))
    {
        state.volume_solver_matrix_nnz = 0.0;
    }
    if (!std::isfinite(state.volume_solver_cell_count))
    {
        state.volume_solver_cell_count = 0.0;
    }
    if (!std::isfinite(state.field_volume_coupling_iterations))
    {
        state.field_volume_coupling_iterations = 0.0;
    }
    if (!std::isfinite(state.field_volume_coupling_converged))
    {
        state.field_volume_coupling_converged = 0.0;
    }
    if (!std::isfinite(state.field_volume_coupling_max_delta))
    {
        state.field_volume_coupling_max_delta = 0.0;
    }
    if (!std::isfinite(state.field_volume_coupling_relaxation_used))
    {
        state.field_volume_coupling_relaxation_used = 0.0;
    }
    if (!std::isfinite(state.external_volume_feedback_blend_factor))
    {
        state.external_volume_feedback_blend_factor = 0.0;
    }
    if (!std::isfinite(state.external_volume_feedback_mismatch_metric))
    {
        state.external_volume_feedback_mismatch_metric = 0.0;
    }
    if (!std::isfinite(state.external_volume_feedback_applied))
    {
        state.external_volume_feedback_applied = 0.0;
    }
    return state;
}

void DensePlasmaSurfaceCharging::rebuildSurfaceCircuit(double conductivity_s_per_m)
{
    if (!circuit_model_)
    {
        return;
    }

    circuit_model_->configure(config_, status_.state.capacitance_per_area_f_per_m2,
                              status_.body_potential_v, status_.patch_potential_v,
                              conductivity_s_per_m);
}

bool DensePlasmaSurfaceCharging::loadLegacyBenchmarkReplay()
{
    legacy_benchmark_replay_.clear();
    legacy_benchmark_body_curve_.clear();
    legacy_benchmark_replay_index_ = 0;
    if (config_.runtime_route != SurfaceRuntimeRoute::LegacyBenchmark)
    {
        return false;
    }
    if (config_.legacy_benchmark_execution_mode !=
        LegacyBenchmarkExecutionMode::ReplayFromReference)
    {
        return false;
    }

    const auto case_definition = loadLegacyBenchmarkCaseDefinition(config_.benchmark_source);
    const auto& patch_samples = case_definition.patch_reference_curve;
    const auto& body_samples = case_definition.body_reference_curve;
    if (patch_samples.empty())
    {
        return false;
    }

    legacy_benchmark_replay_.reserve(patch_samples.size());
    const double patch_area_m2 = std::max(1.0e-16, config_.surface_area_m2);
    for (const auto& patch_sample : patch_samples)
    {
        LegacyBenchmarkReplayRow row;
        row.time_s = patch_sample.time_s;
        row.patch_potential_v = patch_sample.potential_v;
        row.body_potential_v =
            interpolateLegacyBenchmarkPotential(body_samples, patch_sample.time_s);
        row.currents.electron_current_a_per_m2 = patch_sample.je_a_per_m2;
        row.currents.secondary_emission_a_per_m2 = patch_sample.jse_a_per_m2;
        row.currents.backscatter_emission_a_per_m2 = patch_sample.jb_a_per_m2;
        row.currents.ion_current_a_per_m2 = patch_sample.ji_a_per_m2;
        row.currents.ion_secondary_emission_a_per_m2 = patch_sample.jsi_a_per_m2;
        row.currents.photo_emission_a_per_m2 = patch_sample.jph_a_per_m2;
        row.currents.conduction_current_a_per_m2 = patch_sample.jcond_a_per_m2;
        row.currents.total_current_a_per_m2 = patch_sample.jnet_a_per_m2;
        row.net_current_a_per_m2 = patch_sample.jnet_a_per_m2;
        row.branch_current_a = patch_sample.jcond_a_per_m2 * patch_area_m2;
        legacy_benchmark_replay_.push_back(row);
    }

    if (!legacy_benchmark_replay_.empty())
    {
        status_.body_potential_v = legacy_benchmark_replay_.front().body_potential_v;
        status_.patch_potential_v = legacy_benchmark_replay_.front().patch_potential_v;
        status_.state.surface_potential_v = status_.patch_potential_v;
    }
    return !legacy_benchmark_replay_.empty();
}

bool DensePlasmaSurfaceCharging::loadLegacyBenchmarkExecutionCurve()
{
    legacy_benchmark_replay_.clear();
    legacy_benchmark_body_curve_.clear();
    legacy_benchmark_replay_index_ = 0;
    if (config_.runtime_route != SurfaceRuntimeRoute::LegacyBenchmark ||
        config_.legacy_benchmark_execution_mode !=
            LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm ||
        !current_model_)
    {
        return false;
    }

    const bool is_c_family =
        config_.benchmark_source == SurfaceBenchmarkSource::CGeo ||
        config_.benchmark_source == SurfaceBenchmarkSource::CLeoRam ||
        config_.benchmark_source == SurfaceBenchmarkSource::CLeoWake;
    if (!is_c_family)
    {
        return false;
    }

    const std::size_t patch_node_index =
        circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0;
    const std::string patch_node_name =
        circuit_model_ ? circuit_model_->nodeName(patch_node_index) : std::string{"patch"};
    const auto& patch_material = resolvePatchMaterial(config_, patch_node_index, patch_node_name);
    const double patch_area_m2 = std::max(1.0e-16, config_.surface_area_m2);
    const double sheath_length_m = computeEffectiveSheathLength();
    const double body_capacitance_per_area_f_per_m2 = kEpsilon0;
    const double patch_capacitance_per_area_f_per_m2 =
        std::max(1.0e-16, computeCapacitancePerArea());
    const auto benchmark_definition = loadLegacyBenchmarkCaseDefinition(config_.benchmark_source);
    const bool use_reference_body_ladder =
        (config_.benchmark_source == SurfaceBenchmarkSource::CLeoRam ||
         config_.benchmark_source == SurfaceBenchmarkSource::CLeoWake) &&
        !benchmark_definition.body_reference_curve.empty();
    const double initial_step_v =
        config_.regime == SurfaceChargingRegime::GeoKineticPicLike ? 1.0e4 : 1.0e3;
    const double minimum_step_v =
        config_.regime == SurfaceChargingRegime::GeoKineticPicLike ? 1.0e-1 : 1.0e-2;
    const double max_potential_v = std::max(1.0, config_.max_abs_potential_v);

    auto evaluate_role = [&](ReferenceSurfaceRole role, double body_potential_v,
                             double patch_potential_v) {
        const Material::MaterialProperty& material =
            role == ReferenceSurfaceRole::Patch ? patch_material : config_.material;
        const double conductivity =
            role == ReferenceSurfaceRole::Patch ? std::max(0.0, material.getConductivity()) : 0.0;
        const auto state =
            buildRuntimeState(body_potential_v, patch_potential_v, patch_area_m2, conductivity,
                              sheath_length_m,
                              role == ReferenceSurfaceRole::Patch ? patch_node_index : 0,
                              role == ReferenceSurfaceRole::Patch ? patch_node_name
                                                                  : std::string{"body"});
        return current_model_->evaluate(role, state);
    };

    auto build_executor_curve = [&](ReferenceSurfaceRole role, double fixed_body_potential_v,
                                    double fixed_patch_potential_v,
                                    double capacitance_per_area_f_per_m2) {
        std::vector<LegacyBenchmarkReplayRow> rows;
        auto make_row = [&](double potential_v) {
            const double body_v =
                role == ReferenceSurfaceRole::Body ? potential_v : fixed_body_potential_v;
            const double patch_v =
                role == ReferenceSurfaceRole::Patch ? potential_v : fixed_patch_potential_v;
            const auto currents = evaluate_role(role, body_v, patch_v);
            LegacyBenchmarkReplayRow row;
            row.body_potential_v = body_v;
            row.patch_potential_v = patch_v;
            row.currents = currents;
            row.net_current_a_per_m2 = currents.total_current_a_per_m2;
            row.branch_current_a =
                role == ReferenceSurfaceRole::Patch
                    ? currents.conduction_current_a_per_m2 * patch_area_m2
                    : 0.0;
            return row;
        };

        if (role == ReferenceSurfaceRole::Body && use_reference_body_ladder)
        {
            rows.reserve(benchmark_definition.body_reference_curve.size());
            for (const auto& reference_sample : benchmark_definition.body_reference_curve)
            {
                auto row = make_row(reference_sample.potential_v);
                row.time_s = reference_sample.time_s;
                rows.push_back(std::move(row));
            }
            return rows;
        }

        double current_potential_v = 0.0;
        double probe_body_v =
            role == ReferenceSurfaceRole::Body ? current_potential_v : fixed_body_potential_v;
        double probe_patch_v =
            role == ReferenceSurfaceRole::Patch ? current_potential_v : fixed_patch_potential_v;
        const auto initial_currents =
            evaluate_role(role, probe_body_v, probe_patch_v);
        const int initial_sign = signOf(initial_currents.total_current_a_per_m2);

        rows.push_back(make_row(0.0));
        if (initial_sign == 0)
        {
            rows.front().time_s = 0.0;
            return rows;
        }

        double step_v = std::abs(initial_step_v) *
                        static_cast<double>(initial_sign > 0 ? 1.0 : -1.0);
        while (std::abs(step_v) >= minimum_step_v && rows.size() < kLegacyExecutorMaxSamples - 1)
        {
            bool crossed_zero = false;
            while (rows.size() < kLegacyExecutorMaxSamples - 1)
            {
                const double candidate_v = current_potential_v + step_v;
                if (std::abs(candidate_v) > max_potential_v)
                {
                    break;
                }
                auto candidate_row = make_row(candidate_v);
                const int candidate_sign = signOf(candidate_row.net_current_a_per_m2);
                if (candidate_sign == 0 || candidate_sign * initial_sign <= 0)
                {
                    crossed_zero = true;
                    break;
                }
                current_potential_v = candidate_v;
                rows.push_back(std::move(candidate_row));
            }
            if (!crossed_zero && std::abs(current_potential_v + step_v) > max_potential_v)
            {
                break;
            }
            step_v /= 10.0;
        }

        rows.front().time_s = 0.0;
        for (std::size_t i = 1; i < rows.size(); ++i)
        {
            const double previous_v = role == ReferenceSurfaceRole::Body ? rows[i - 1].body_potential_v
                                                                         : rows[i - 1].patch_potential_v;
            const double current_v = role == ReferenceSurfaceRole::Body ? rows[i].body_potential_v
                                                                        : rows[i].patch_potential_v;
            const double average_current =
                0.5 * (rows[i - 1].net_current_a_per_m2 + rows[i].net_current_a_per_m2);
            const double safe_current =
                std::max(1.0e-18, std::abs(average_current));
            rows[i].time_s =
                rows[i - 1].time_s +
                capacitance_per_area_f_per_m2 * std::abs(current_v - previous_v) / safe_current;
        }
        return rows;
    };

    const auto body_curve = build_executor_curve(ReferenceSurfaceRole::Body, 0.0, 0.0,
                                                 body_capacitance_per_area_f_per_m2);
    if (body_curve.empty())
    {
        return false;
    }
    legacy_benchmark_body_curve_ = body_curve;

    const double body_equilibrium_v = body_curve.back().body_potential_v;
    auto patch_curve = build_executor_curve(ReferenceSurfaceRole::Patch, body_equilibrium_v, 0.0,
                                            patch_capacitance_per_area_f_per_m2);
    if (patch_curve.empty())
    {
        return false;
    }

    for (auto& row : patch_curve)
    {
        row.body_potential_v = body_equilibrium_v;
    }

    legacy_benchmark_replay_ = std::move(patch_curve);
    status_.body_potential_v = legacy_benchmark_replay_.front().body_potential_v;
    status_.patch_potential_v = legacy_benchmark_replay_.front().patch_potential_v;
    status_.state.surface_potential_v = status_.patch_potential_v;
    return true;
}

bool DensePlasmaSurfaceCharging::advanceLegacyBenchmarkReplay(double)
{
    if (legacy_benchmark_replay_.empty())
    {
        return false;
    }

    const std::size_t sample_index =
        std::min(legacy_benchmark_replay_index_, legacy_benchmark_replay_.size() - 1);
    const auto& row = legacy_benchmark_replay_[sample_index];
    const double patch_body_delta_v = row.patch_potential_v - row.body_potential_v;
    const double effective_conductivity =
        std::abs(patch_body_delta_v) > 1.0e-9
            ? std::abs(row.currents.conduction_current_a_per_m2) * config_.dielectric_thickness_m /
                  std::abs(patch_body_delta_v)
            : computeEffectiveConductivity(row.patch_potential_v);
    const double effective_sheath_length = computeEffectiveSheathLength();
    const auto runtime_state =
        buildRuntimeState(row.body_potential_v, row.patch_potential_v, config_.surface_area_m2,
                          effective_conductivity, effective_sheath_length);

    status_.currents = row.currents;
    status_.currents.current_derivative_a_per_m2_per_v = 0.0;
    status_.state.surface_potential_v = row.patch_potential_v;
    status_.state.capacitance_per_area_f_per_m2 = computeCapacitancePerArea(runtime_state);
    status_.state.surface_charge_density_c_per_m2 =
        status_.state.capacitance_per_area_f_per_m2 * patch_body_delta_v;
    status_.state.stored_energy_j_per_m2 =
        0.5 * status_.state.capacitance_per_area_f_per_m2 * patch_body_delta_v * patch_body_delta_v;
    status_.body_potential_v = row.body_potential_v;
    status_.patch_potential_v = row.patch_potential_v;
    status_.circuit_branch_current_a = row.branch_current_a;
    status_.pic_recalibrated = false;
    status_.equilibrium_error =
        row.patch_potential_v - legacy_benchmark_replay_.back().patch_potential_v;
    status_.equilibrium_reached = sample_index + 1 >= legacy_benchmark_replay_.size();
    status_.time_s = row.time_s;
    status_.steps_completed += 1;

    history_time_.push_back(status_.time_s);
    history_potential_.push_back(status_.state.surface_potential_v);
    history_charge_.push_back(status_.state.surface_charge_density_c_per_m2);
    history_current_.push_back(status_.currents.total_current_a_per_m2);
    history_electron_current_.push_back(status_.currents.electron_current_a_per_m2);
    history_ion_current_.push_back(status_.currents.ion_current_a_per_m2);
    history_secondary_current_.push_back(status_.currents.secondary_emission_a_per_m2);
    history_ion_secondary_current_.push_back(status_.currents.ion_secondary_emission_a_per_m2);
    history_backscatter_current_.push_back(status_.currents.backscatter_emission_a_per_m2);
    history_photo_current_.push_back(status_.currents.photo_emission_a_per_m2);
    history_thermionic_current_.push_back(status_.currents.thermionic_emission_a_per_m2);
    history_field_emission_current_.push_back(status_.currents.field_emission_a_per_m2);
    history_leakage_current_.push_back(status_.currents.conduction_current_a_per_m2);
    history_ram_current_.push_back(status_.currents.ram_ion_current_a_per_m2);
    history_body_potential_.push_back(status_.body_potential_v);
    history_patch_potential_.push_back(status_.patch_potential_v);
    history_circuit_branch_current_.push_back(status_.circuit_branch_current_a);
    history_current_derivative_.push_back(status_.currents.current_derivative_a_per_m2_per_v);
    history_pic_recalibration_marker_.push_back(0.0);
    history_live_pic_electron_current_.push_back(0.0);
    history_live_pic_ion_current_.push_back(0.0);
    history_live_pic_net_current_.push_back(0.0);
    history_live_pic_derivative_.push_back(0.0);
    history_live_pic_collision_count_.push_back(0.0);
    history_live_pic_mcc_enabled_.push_back(0.0);
    history_net_current_.push_back(row.net_current_a_per_m2);
    history_capacitance_.push_back(status_.state.capacitance_per_area_f_per_m2);
    history_effective_conductivity_.push_back(effective_conductivity);
    history_effective_sheath_length_.push_back(effective_sheath_length);
    history_normal_electric_field_.push_back(runtime_state.normal_electric_field_v_per_m);
    history_local_charge_density_.push_back(runtime_state.local_charge_density_c_per_m3);
    history_adaptive_time_step_.push_back(
        sample_index == 0 ? row.time_s
                          : std::max(0.0, row.time_s - legacy_benchmark_replay_[sample_index - 1].time_s));
    history_electron_calibration_factor_.push_back(1.0);
    history_ion_calibration_factor_.push_back(1.0);
    history_equilibrium_potential_.push_back(legacy_benchmark_replay_.back().patch_potential_v);
    history_equilibrium_error_.push_back(status_.equilibrium_error);
    history_shared_surface_pre_global_solve_patch_potential_spread_.push_back(0.0);
    history_shared_surface_patch_potential_spread_reduction_v_.push_back(0.0);
    history_shared_surface_patch_potential_spread_reduction_ratio_.push_back(0.0);
    history_shared_surface_current_matrix_coupling_active_.push_back(0.0);
    history_shared_surface_current_matrix_coupling_offdiag_entries_.push_back(0.0);
    history_shared_surface_global_coupled_solve_active_.push_back(0.0);
    history_shared_surface_global_coupled_solve_iterations_.push_back(0.0);
    history_shared_surface_global_coupled_solve_converged_.push_back(0.0);
    history_shared_surface_global_coupled_solve_max_delta_v_.push_back(0.0);
    history_shared_surface_live_pic_coupled_refresh_active_.push_back(0.0);
    history_shared_surface_live_pic_coupled_refresh_count_.push_back(0.0);
    history_shared_surface_particle_transport_coupling_active_.push_back(0.0);
    history_shared_surface_particle_transport_offdiag_entries_.push_back(0.0);
    history_shared_surface_particle_transport_total_conductance_s_.push_back(0.0);
    history_shared_surface_particle_transport_conservation_error_a_per_v_.push_back(0.0);
    history_shared_surface_particle_transport_charge_c_.push_back(shared_particle_transport_charge_c_);
    history_shared_surface_particle_transport_reference_shift_v_.push_back(
        shared_particle_transport_reference_shift_v_);
    if (circuit_node_potentials_.size() < 2)
    {
        circuit_node_potentials_.assign(2, 0.0);
    }
    circuit_node_potentials_[0] = status_.body_potential_v;
    circuit_node_potentials_[1] = status_.patch_potential_v;
    std::vector<double> branch_currents_a;
    if (circuit_model_ && circuit_model_->branchCount() > 0)
    {
        branch_currents_a.assign(circuit_model_->branchCount(), 0.0);
        branch_currents_a[0] = status_.circuit_branch_current_a;
    }
    else
    {
        branch_currents_a.assign(1, status_.circuit_branch_current_a);
    }
    appendSurfaceNodeDiagnostics(branch_currents_a, effective_sheath_length);

    if (legacy_benchmark_replay_index_ + 1 < legacy_benchmark_replay_.size())
    {
        ++legacy_benchmark_replay_index_;
    }
    return true;
}

void DensePlasmaSurfaceCharging::appendSurfaceNodeDiagnostics(
    const std::vector<double>& branch_currents_a, double effective_sheath_length_m)
{
    if (circuit_node_potentials_.empty())
    {
        circuit_node_potentials_ = {status_.body_potential_v, status_.patch_potential_v};
    }

    const std::size_t node_count = circuit_node_potentials_.size();
    auto ensure_node_storage = [&](std::vector<std::vector<double>>& history) {
        if (history.size() != node_count)
        {
            history.assign(node_count, {});
        }
    };
    ensure_node_storage(history_surface_node_potentials_);
    ensure_node_storage(history_surface_node_total_currents_);
    ensure_node_storage(history_surface_node_electron_currents_);
    ensure_node_storage(history_surface_node_ion_currents_);
    ensure_node_storage(history_surface_node_secondary_currents_);
    ensure_node_storage(history_surface_node_ion_secondary_currents_);
    ensure_node_storage(history_surface_node_backscatter_currents_);
    ensure_node_storage(history_surface_node_photo_currents_);
    ensure_node_storage(history_surface_node_thermionic_currents_);
    ensure_node_storage(history_surface_node_field_emission_currents_);
    ensure_node_storage(history_surface_node_conduction_currents_);
    ensure_node_storage(history_surface_node_ram_currents_);
    ensure_node_storage(history_surface_node_current_derivatives_);
    ensure_node_storage(history_surface_node_effective_sheath_lengths_);
    ensure_node_storage(history_surface_node_normal_electric_fields_);
    ensure_node_storage(history_surface_node_local_charge_densities_);
    ensure_node_storage(history_surface_node_propagated_reference_potentials_);
    ensure_node_storage(history_surface_node_field_solver_reference_potentials_);
    ensure_node_storage(history_surface_node_shared_runtime_enabled_);
    ensure_node_storage(history_surface_node_shared_patch_potentials_);
    ensure_node_storage(history_surface_node_shared_patch_areas_);
    ensure_node_storage(history_surface_node_shared_reference_potentials_);
    ensure_node_storage(history_surface_node_shared_sheath_charges_);
    ensure_node_storage(history_surface_node_distributed_particle_transport_charges_);
    ensure_node_storage(history_surface_node_distributed_particle_transport_reference_shifts_);
    ensure_node_storage(history_surface_node_distributed_particle_transport_net_fluxes_);
    ensure_node_storage(history_surface_node_graph_capacitance_diagonals_);
    ensure_node_storage(history_surface_node_graph_capacitance_row_sums_);
    ensure_node_storage(history_surface_node_field_solver_coupling_gains_);
    ensure_node_storage(history_surface_node_field_solver_capacitance_scales_);
    ensure_node_storage(history_surface_node_pseudo_volumes_);
    ensure_node_storage(history_surface_node_volume_projection_weight_sums_);
    ensure_node_storage(history_surface_node_volume_mesh_coupling_gains_);
    ensure_node_storage(history_surface_node_volume_potentials_);
    ensure_node_storage(history_surface_node_deposited_charges_);
    ensure_node_storage(history_surface_node_poisson_residuals_);
    ensure_node_storage(history_surface_node_volume_solver_mode_ids_);
    ensure_node_storage(history_surface_node_volume_solver_iterations_);
    ensure_node_storage(history_surface_node_volume_solver_linear_iterations_);
    ensure_node_storage(history_surface_node_volume_solver_converged_);
    ensure_node_storage(history_surface_node_volume_solver_residual_norms_);
    ensure_node_storage(history_surface_node_volume_solver_max_deltas_);
    ensure_node_storage(history_surface_node_volume_solver_matrix_nnzs_);
    ensure_node_storage(history_surface_node_volume_solver_cell_counts_);
    ensure_node_storage(history_surface_node_field_volume_coupling_iterations_);
    ensure_node_storage(history_surface_node_field_volume_coupling_converged_);
    ensure_node_storage(history_surface_node_field_volume_coupling_max_deltas_);
    ensure_node_storage(history_surface_node_field_volume_coupling_relaxation_used_);
    ensure_node_storage(history_surface_node_external_volume_feedback_blend_factors_);
    ensure_node_storage(history_surface_node_external_volume_feedback_mismatch_metrics_);
    ensure_node_storage(history_surface_node_external_volume_feedback_applied_);
    ensure_node_storage(history_surface_node_electron_calibration_factors_);
    ensure_node_storage(history_surface_node_ion_calibration_factors_);
    ensure_node_storage(history_surface_node_live_pic_electron_currents_);
    ensure_node_storage(history_surface_node_live_pic_ion_currents_);
    ensure_node_storage(history_surface_node_live_pic_net_currents_);
    ensure_node_storage(history_surface_node_live_pic_derivatives_);
    ensure_node_storage(history_surface_node_live_pic_collision_counts_);
    ensure_node_storage(history_surface_node_live_pic_mcc_enabled_);

    if (history_surface_branch_currents_.size() != branch_currents_a.size())
    {
        history_surface_branch_currents_.assign(branch_currents_a.size(), {});
    }
    if (history_surface_branch_conductances_.size() != branch_currents_a.size())
    {
        history_surface_branch_conductances_.assign(branch_currents_a.size(), {});
    }
    if (history_surface_branch_voltage_drops_.size() != branch_currents_a.size())
    {
        history_surface_branch_voltage_drops_.assign(branch_currents_a.size(), {});
    }
    if (history_surface_branch_power_w_.size() != branch_currents_a.size())
    {
        history_surface_branch_power_w_.assign(branch_currents_a.size(), {});
    }
    if (history_surface_branch_mutual_capacitances_.size() != branch_currents_a.size())
    {
        history_surface_branch_mutual_capacitances_.assign(branch_currents_a.size(), {});
    }
    for (std::size_t branch_index = 0; branch_index < branch_currents_a.size(); ++branch_index)
    {
        history_surface_branch_currents_[branch_index].push_back(branch_currents_a[branch_index]);
        const double conductance_s =
            circuit_model_ ? circuit_model_->branchConductanceS(branch_index) : 0.0;
        history_surface_branch_conductances_[branch_index].push_back(conductance_s);
        double voltage_drop_v = 0.0;
        if (circuit_model_)
        {
            const auto from_node_index = circuit_model_->branchFromNodeIndex(branch_index);
            const auto to_node_index = circuit_model_->branchToNodeIndex(branch_index);
            const double from_potential_v =
                from_node_index < circuit_node_potentials_.size() ? circuit_node_potentials_[from_node_index]
                                                                  : 0.0;
            const double to_potential_v =
                to_node_index < circuit_node_potentials_.size() ? circuit_node_potentials_[to_node_index]
                                                                : 0.0;
            voltage_drop_v =
                from_potential_v - to_potential_v - circuit_model_->branchBiasV(branch_index);
        }
        history_surface_branch_voltage_drops_[branch_index].push_back(voltage_drop_v);
        history_surface_branch_power_w_[branch_index].push_back(branch_currents_a[branch_index] *
                                                                voltage_drop_v);
        const double mutual_capacitance_f =
            graph_capacitance_matrix_provider_
                ? graph_capacitance_matrix_provider_->mutualCapacitanceF(config_, circuit_model_.get(),
                                                                         branch_index)
                : 0.0;
        history_surface_branch_mutual_capacitances_[branch_index].push_back(mutual_capacitance_f);
    }

    const std::size_t body_node_index =
        circuit_model_ ? circuit_model_->bodyNodeIndex() : static_cast<std::size_t>(0);
    const double body_potential_v =
        body_node_index < circuit_node_potentials_.size() ? circuit_node_potentials_[body_node_index]
                                                          : status_.body_potential_v;
    const double default_patch_potential_v =
        circuit_model_ && circuit_model_->primaryPatchNodeIndex() < circuit_node_potentials_.size()
            ? circuit_node_potentials_[circuit_model_->primaryPatchNodeIndex()]
            : status_.patch_potential_v;

    auto patch_branch_index = [&](std::size_t node_index) {
        if (!circuit_model_)
        {
            return node_index == 1 ? static_cast<std::size_t>(0) : std::numeric_limits<std::size_t>::max();
        }
        for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
        {
            if (circuit_model_->patchNodeIndex(patch_ordinal) == node_index)
            {
                return circuit_model_->patchToBodyBranchIndex(patch_ordinal);
            }
        }
        return std::numeric_limits<std::size_t>::max();
    };

    double shared_patch_potential_v = computeSharedSurfacePatchPotentialV();
    double shared_patch_area_m2 = computeSharedSurfacePatchAreaM2();
    double shared_effective_sheath_length_m =
        computeSharedSurfaceEffectiveSheathLengthM(effective_sheath_length_m);
    double shared_reference_potential_v = shared_patch_potential_v;
    double shared_sheath_charge_c = 0.0;
    bool shared_reference_captured = false;
    std::vector<SurfaceModelRuntimeState> runtime_states;
    runtime_states.reserve(node_count);

    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        const bool is_patch = circuit_model_ ? circuit_model_->nodeIsPatch(node_index) : node_index == 1;
        const double node_potential_v = circuit_node_potentials_[node_index];
        const double patch_potential_v = is_patch ? node_potential_v : default_patch_potential_v;
        const double node_area_m2 =
            std::max(1.0e-16, circuit_model_ ? circuit_model_->nodeAreaM2(node_index)
                                             : config_.surface_area_m2);
        const std::string node_name =
            circuit_model_ ? circuit_model_->nodeName(node_index)
                           : (is_patch ? std::string{"patch"} : std::string{"body"});
        const double node_conductivity =
            is_patch
                ? computeEffectiveConductivity(node_potential_v,
                                               resolvePatchMaterial(config_, node_index, node_name))
                : 0.0;
        const auto runtime_state =
            buildRuntimeState(body_potential_v, patch_potential_v, node_area_m2, node_conductivity,
                              effective_sheath_length_m, node_index, node_name);
        runtime_states.push_back(runtime_state);

        SurfaceCurrents node_currents;
        if (current_model_)
        {
            node_currents = current_model_->evaluate(
                is_patch ? ReferenceSurfaceRole::Patch : ReferenceSurfaceRole::Body, runtime_state);
        }
        const double electron_calibration_factor =
            current_model_ ? current_model_->electronCalibrationFactor(runtime_state) : 1.0;
        const double ion_calibration_factor =
            current_model_ ? current_model_->ionCalibrationFactor(runtime_state) : 1.0;
        const PicMccCurrentSample live_pic_sample =
            current_model_ ? current_model_->latestLivePicSample(runtime_state) : PicMccCurrentSample{};
        if (is_patch)
        {
            const std::size_t branch_index = patch_branch_index(node_index);
            if (branch_index < branch_currents_a.size())
            {
                node_currents.conduction_current_a_per_m2 = branch_currents_a[branch_index] / node_area_m2;
            }
        }
        node_currents.total_current_a_per_m2 =
            node_currents.electron_current_a_per_m2 + node_currents.ion_current_a_per_m2 +
            node_currents.secondary_emission_a_per_m2 + node_currents.ion_secondary_emission_a_per_m2 +
            node_currents.backscatter_emission_a_per_m2 + node_currents.photo_emission_a_per_m2 +
            node_currents.thermionic_emission_a_per_m2 + node_currents.field_emission_a_per_m2 +
            node_currents.conduction_current_a_per_m2 + node_currents.ram_ion_current_a_per_m2;

        history_surface_node_potentials_[node_index].push_back(node_potential_v);
        history_surface_node_total_currents_[node_index].push_back(node_currents.total_current_a_per_m2);
        history_surface_node_electron_currents_[node_index].push_back(node_currents.electron_current_a_per_m2);
        history_surface_node_ion_currents_[node_index].push_back(node_currents.ion_current_a_per_m2);
        history_surface_node_secondary_currents_[node_index].push_back(
            node_currents.secondary_emission_a_per_m2);
        history_surface_node_ion_secondary_currents_[node_index].push_back(
            node_currents.ion_secondary_emission_a_per_m2);
        history_surface_node_backscatter_currents_[node_index].push_back(
            node_currents.backscatter_emission_a_per_m2);
        history_surface_node_photo_currents_[node_index].push_back(node_currents.photo_emission_a_per_m2);
        history_surface_node_thermionic_currents_[node_index].push_back(
            node_currents.thermionic_emission_a_per_m2);
        history_surface_node_field_emission_currents_[node_index].push_back(
            node_currents.field_emission_a_per_m2);
        history_surface_node_conduction_currents_[node_index].push_back(
            node_currents.conduction_current_a_per_m2);
        history_surface_node_ram_currents_[node_index].push_back(node_currents.ram_ion_current_a_per_m2);
        history_surface_node_current_derivatives_[node_index].push_back(
            node_currents.current_derivative_a_per_m2_per_v);
        history_surface_node_effective_sheath_lengths_[node_index].push_back(
            runtime_state.effective_sheath_length_m);
        history_surface_node_normal_electric_fields_[node_index].push_back(
            runtime_state.normal_electric_field_v_per_m);
        history_surface_node_local_charge_densities_[node_index].push_back(
            runtime_state.local_charge_density_c_per_m3);
        history_surface_node_propagated_reference_potentials_[node_index].push_back(
            runtime_state.propagated_reference_potential_v);
        history_surface_node_field_solver_reference_potentials_[node_index].push_back(
            runtime_state.field_solver_reference_potential_v);
        history_surface_node_shared_runtime_enabled_[node_index].push_back(
            runtime_state.shared_runtime_enabled ? 1.0 : 0.0);
        history_surface_node_shared_patch_potentials_[node_index].push_back(
            runtime_state.shared_patch_potential_v);
        history_surface_node_shared_patch_areas_[node_index].push_back(
            runtime_state.shared_patch_area_m2);
        history_surface_node_shared_reference_potentials_[node_index].push_back(
            runtime_state.shared_reference_potential_v);
        const double node_shared_sheath_charge_c =
            kEpsilon0 * runtime_state.shared_patch_area_m2 *
            (runtime_state.shared_reference_potential_v - runtime_state.shared_patch_potential_v) /
            std::max(1.0e-6, runtime_state.shared_effective_sheath_length_m);
        history_surface_node_shared_sheath_charges_[node_index].push_back(node_shared_sheath_charge_c);
        history_surface_node_distributed_particle_transport_charges_[node_index].push_back(
            runtime_state.distributed_particle_transport_charge_c);
        history_surface_node_distributed_particle_transport_reference_shifts_[node_index].push_back(
            runtime_state.distributed_particle_transport_reference_shift_v);
        history_surface_node_distributed_particle_transport_net_fluxes_[node_index].push_back(
            runtime_state.distributed_particle_transport_net_flux_a);
        history_surface_node_graph_capacitance_diagonals_[node_index].push_back(
            runtime_state.graph_capacitance_diagonal_f);
        history_surface_node_graph_capacitance_row_sums_[node_index].push_back(
            runtime_state.graph_capacitance_row_sum_f);
        history_surface_node_field_solver_coupling_gains_[node_index].push_back(
            runtime_state.field_solver_coupling_gain);
        history_surface_node_field_solver_capacitance_scales_[node_index].push_back(
            runtime_state.field_solver_capacitance_scale);
        history_surface_node_pseudo_volumes_[node_index].push_back(runtime_state.pseudo_volume_m3);
        history_surface_node_volume_projection_weight_sums_[node_index].push_back(
            runtime_state.volume_projection_weight_sum);
        history_surface_node_volume_mesh_coupling_gains_[node_index].push_back(
            runtime_state.volume_mesh_coupling_gain);
        history_surface_node_volume_potentials_[node_index].push_back(runtime_state.volume_potential_v);
        history_surface_node_deposited_charges_[node_index].push_back(runtime_state.deposited_charge_c);
        history_surface_node_poisson_residuals_[node_index].push_back(
            runtime_state.poisson_residual_v_m);
        history_surface_node_volume_solver_mode_ids_[node_index].push_back(
            runtime_state.volume_solver_mode_id);
        history_surface_node_volume_solver_iterations_[node_index].push_back(
            runtime_state.volume_solver_iterations);
        history_surface_node_volume_solver_linear_iterations_[node_index].push_back(
            runtime_state.volume_solver_linear_iterations);
        history_surface_node_volume_solver_converged_[node_index].push_back(
            runtime_state.volume_solver_converged);
        history_surface_node_volume_solver_residual_norms_[node_index].push_back(
            runtime_state.volume_solver_residual_norm);
        history_surface_node_volume_solver_max_deltas_[node_index].push_back(
            runtime_state.volume_solver_max_delta_v);
        history_surface_node_volume_solver_matrix_nnzs_[node_index].push_back(
            runtime_state.volume_solver_matrix_nnz);
        history_surface_node_volume_solver_cell_counts_[node_index].push_back(
            runtime_state.volume_solver_cell_count);
        history_surface_node_field_volume_coupling_iterations_[node_index].push_back(
            runtime_state.field_volume_coupling_iterations);
        history_surface_node_field_volume_coupling_converged_[node_index].push_back(
            runtime_state.field_volume_coupling_converged);
        history_surface_node_field_volume_coupling_max_deltas_[node_index].push_back(
            runtime_state.field_volume_coupling_max_delta);
        history_surface_node_field_volume_coupling_relaxation_used_[node_index].push_back(
            runtime_state.field_volume_coupling_relaxation_used);
        history_surface_node_external_volume_feedback_blend_factors_[node_index].push_back(
            runtime_state.external_volume_feedback_blend_factor);
        history_surface_node_external_volume_feedback_mismatch_metrics_[node_index].push_back(
            runtime_state.external_volume_feedback_mismatch_metric);
        history_surface_node_external_volume_feedback_applied_[node_index].push_back(
            runtime_state.external_volume_feedback_applied);
        history_surface_node_electron_calibration_factors_[node_index].push_back(
            electron_calibration_factor);
        history_surface_node_ion_calibration_factors_[node_index].push_back(ion_calibration_factor);
        history_surface_node_live_pic_electron_currents_[node_index].push_back(
            live_pic_sample.electron_collection_current_density_a_per_m2);
        history_surface_node_live_pic_ion_currents_[node_index].push_back(
            live_pic_sample.ion_collection_current_density_a_per_m2);
        history_surface_node_live_pic_net_currents_[node_index].push_back(
            live_pic_sample.net_collection_current_density_a_per_m2);
        history_surface_node_live_pic_derivatives_[node_index].push_back(
            live_pic_sample.current_derivative_a_per_m2_per_v);
        history_surface_node_live_pic_collision_counts_[node_index].push_back(
            static_cast<double>(live_pic_sample.total_collisions));
        history_surface_node_live_pic_mcc_enabled_[node_index].push_back(
            live_pic_sample.mcc_enabled ? 1.0 : 0.0);
        if (runtime_state.shared_runtime_enabled && !shared_reference_captured)
        {
            shared_patch_potential_v = runtime_state.shared_patch_potential_v;
            shared_patch_area_m2 = runtime_state.shared_patch_area_m2;
            shared_effective_sheath_length_m = runtime_state.shared_effective_sheath_length_m;
            shared_reference_potential_v = runtime_state.shared_reference_potential_v;
            shared_sheath_charge_c = node_shared_sheath_charge_c;
            shared_reference_captured = true;
        }
    }

    history_shared_surface_patch_potential_.push_back(shared_patch_potential_v);
    history_shared_surface_patch_area_.push_back(shared_patch_area_m2);
    history_shared_surface_reference_potential_.push_back(shared_reference_potential_v);
    history_shared_surface_effective_sheath_length_.push_back(shared_effective_sheath_length_m);
    history_shared_surface_sheath_charge_.push_back(shared_sheath_charge_c);
    history_shared_surface_runtime_enabled_.push_back(useSharedSurfacePicRuntime() ? 1.0 : 0.0);
    refreshGlobalKernelStates(runtime_states);
    for (const auto& node_state : global_particle_domain_state_.nodes)
    {
        const auto node_index = node_state.node_index;
        if (node_index < history_surface_node_distributed_particle_transport_charges_.size() and
            !history_surface_node_distributed_particle_transport_charges_[node_index].empty())
        {
            history_surface_node_distributed_particle_transport_charges_[node_index].back() =
                node_state.charge_c;
        }
        if (
            node_index < history_surface_node_distributed_particle_transport_reference_shifts_.size()
            && !history_surface_node_distributed_particle_transport_reference_shifts_[node_index]
                    .empty()
        )
        {
            history_surface_node_distributed_particle_transport_reference_shifts_[node_index].back() =
                node_state.reference_shift_v;
        }
        if (node_index < history_surface_node_distributed_particle_transport_net_fluxes_.size() &&
            !history_surface_node_distributed_particle_transport_net_fluxes_[node_index].empty())
        {
            history_surface_node_distributed_particle_transport_net_fluxes_[node_index].back() =
                node_state.net_flux_a;
        }
    }
    for (const auto& node_state : global_sheath_field_solve_state_.nodes)
    {
        const auto node_index = node_state.node_index;
        if (node_index < history_surface_node_field_solver_reference_potentials_.size() &&
            !history_surface_node_field_solver_reference_potentials_[node_index].empty())
        {
            history_surface_node_field_solver_reference_potentials_[node_index].back() =
                global_sheath_field_solve_state_.global_reference_potential_v;
        }
        if (node_index < history_surface_node_shared_reference_potentials_.size() &&
            !history_surface_node_shared_reference_potentials_[node_index].empty())
        {
            history_surface_node_shared_reference_potentials_[node_index].back() =
                node_state.reference_potential_v;
        }
        if (node_index < history_surface_node_normal_electric_fields_.size() &&
            !history_surface_node_normal_electric_fields_[node_index].empty())
        {
            history_surface_node_normal_electric_fields_[node_index].back() =
                node_state.normal_electric_field_v_per_m;
        }
        if (node_index < history_surface_node_local_charge_densities_.size() &&
            !history_surface_node_local_charge_densities_[node_index].empty())
        {
            history_surface_node_local_charge_densities_[node_index].back() =
                node_state.local_charge_density_c_per_m3;
        }
    }
    if (global_sheath_field_solve_state_.active)
    {
        if (!history_shared_surface_patch_potential_.empty())
        {
            history_shared_surface_patch_potential_.back() =
                global_sheath_field_solve_state_.global_patch_potential_v;
        }
        if (!history_shared_surface_reference_potential_.empty())
        {
            history_shared_surface_reference_potential_.back() =
                global_sheath_field_solve_state_.global_reference_potential_v;
        }
        if (!history_shared_surface_effective_sheath_length_.empty())
        {
            history_shared_surface_effective_sheath_length_.back() =
                global_sheath_field_solve_state_.effective_sheath_length_m;
        }
        if (!history_shared_surface_sheath_charge_.empty())
        {
            const double shared_patch_area_m2 =
                history_shared_surface_patch_area_.empty()
                    ? 0.0
                    : std::max(1.0e-16, history_shared_surface_patch_area_.back());
            history_shared_surface_sheath_charge_.back() =
                kEpsilon0 * shared_patch_area_m2 *
                (global_sheath_field_solve_state_.global_reference_potential_v -
                 global_sheath_field_solve_state_.global_patch_potential_v) /
                std::max(1.0e-6, global_sheath_field_solve_state_.effective_sheath_length_m);
        }
    }
    if (!global_sheath_field_previous_node_charge_c_.empty())
    {
        for (const auto& node_state : global_sheath_field_solve_state_.nodes)
        {
            if (node_state.node_index < global_sheath_field_previous_node_charge_c_.size())
            {
                global_sheath_field_previous_node_charge_c_[node_state.node_index] =
                    node_state.sheath_charge_c;
            }
        }
    }
}

void DensePlasmaSurfaceCharging::refreshGlobalKernelStates(
    const std::vector<SurfaceModelRuntimeState>& runtime_states)
{
    if (global_particle_domain_state_.active || global_sheath_field_solve_state_.active)
    {
        global_particle_domain_state_.total_charge_c = shared_particle_transport_charge_c_;
        global_particle_domain_state_.total_reference_shift_v =
            shared_particle_transport_reference_shift_v_;
        global_particle_domain_state_.total_node_charge_c = 0.0;
        global_particle_domain_state_.total_node_flux_a = 0.0;
        global_particle_domain_state_.charge_conservation_error_c = 0.0;
        global_particle_domain_state_.flux_conservation_error_a = 0.0;

        double total_sheath_area_m2 = 0.0;
        double weighted_patch_potential_v = 0.0;
        double weighted_reference_potential_v = 0.0;
        double weighted_normal_field_v_per_m = 0.0;
        double weighted_charge_density_c_per_m3 = 0.0;

        for (const auto& runtime_state : runtime_states)
        {
            if (!runtime_state.shared_runtime_enabled)
            {
                continue;
            }
            if (auto* particle_node =
                    const_cast<GlobalParticleDomainNodeState*>(
                        findGlobalParticleDomainNodeState(runtime_state.node_index)))
            {
                particle_node->patch_potential_v = runtime_state.patch_potential_v;
                particle_node->area_m2 = runtime_state.surface_area_m2;
                particle_node->shared_reference_potential_v =
                    runtime_state.shared_reference_potential_v;
                particle_node->charge_c = runtime_state.distributed_particle_transport_charge_c;
                particle_node->reference_shift_v =
                    runtime_state.distributed_particle_transport_reference_shift_v;
                particle_node->net_flux_a = runtime_state.distributed_particle_transport_net_flux_a;
                global_particle_domain_state_.total_node_charge_c += particle_node->charge_c;
                global_particle_domain_state_.total_node_flux_a += particle_node->net_flux_a;
            }
            if (auto* sheath_node =
                    const_cast<GlobalSheathFieldNodeState*>(
                        findGlobalSheathFieldNodeState(runtime_state.node_index)))
            {
                sheath_node->patch_potential_v = runtime_state.patch_potential_v;
                sheath_node->area_m2 = runtime_state.surface_area_m2;
                sheath_node->reference_potential_v = runtime_state.shared_reference_potential_v;
                sheath_node->normal_electric_field_v_per_m =
                    runtime_state.normal_electric_field_v_per_m;
                sheath_node->local_charge_density_c_per_m3 =
                    runtime_state.local_charge_density_c_per_m3;
                sheath_node->particle_charge_c =
                    runtime_state.distributed_particle_transport_charge_c;
                sheath_node->particle_flux_a = runtime_state.distributed_particle_transport_net_flux_a;
                const double area_m2 = std::max(1.0e-16, runtime_state.surface_area_m2);
                total_sheath_area_m2 += area_m2;
                weighted_patch_potential_v += area_m2 * sheath_node->patch_potential_v;
                weighted_reference_potential_v += area_m2 * sheath_node->reference_potential_v;
                weighted_normal_field_v_per_m +=
                    area_m2 * sheath_node->normal_electric_field_v_per_m;
                weighted_charge_density_c_per_m3 +=
                    area_m2 * sheath_node->local_charge_density_c_per_m3;
            }
        }

        rebuildOwnedGlobalParticleDomainEdges();
        global_particle_domain_state_.charge_conservation_error_c =
            std::abs(global_particle_domain_state_.total_node_charge_c -
                     global_particle_domain_state_.total_charge_c);
        global_particle_domain_state_.flux_conservation_error_a =
            std::abs(global_particle_domain_state_.total_node_flux_a);
        if (circuit_model_ != nullptr)
        {
            std::vector<double> node_potentials_v(circuit_model_->nodeCount(), 0.0);
            for (const auto& node_state : global_particle_domain_state_.nodes)
            {
                if (node_state.node_index < node_potentials_v.size())
                {
                    node_potentials_v[node_state.node_index] = node_state.patch_potential_v;
                }
            }
            if (!global_particle_repository_state_.active)
            {
                const std::vector<double> zero_node_values(node_potentials_v.size(), 0.0);
                populateOwnedGlobalParticleRepositoryStateFromTransportBuffers(
                    node_potentials_v, shared_particle_transport_node_charge_c_,
                    shared_particle_transport_node_net_flux_a_,
                    shared_particle_transport_node_charge_c_, zero_node_values,
                    zero_node_values, zero_node_values,
                    shared_particle_transport_edge_charge_matrix_c_,
                    shared_particle_transport_edge_target_charge_matrix_c_,
                    shared_particle_transport_edge_operator_drive_matrix_c_,
                    shared_particle_transport_exchange_flux_matrix_a_,
                    global_particle_transport_conductance_matrix_s_, 0.0);
            }
        }

        if (global_sheath_field_solve_state_.active && total_sheath_area_m2 > 0.0)
        {
            global_sheath_field_solve_state_.global_patch_potential_v =
                weighted_patch_potential_v / total_sheath_area_m2;
            global_sheath_field_solve_state_.global_reference_potential_v =
                weighted_reference_potential_v / total_sheath_area_m2;
            global_sheath_field_solve_state_.global_normal_electric_field_v_per_m =
                weighted_normal_field_v_per_m / total_sheath_area_m2;
            global_sheath_field_solve_state_.global_local_charge_density_c_per_m3 =
                weighted_charge_density_c_per_m3 / total_sheath_area_m2;
        }
        return;
    }

    global_particle_domain_state_ = GlobalParticleDomainState{};
    global_particle_domain_state_.bookkeeping_mode = "owned_global_particle_domain_state_v2";
    global_particle_repository_state_ = GlobalParticleRepositoryState{};
    global_particle_repository_state_.bookkeeping_mode =
        "owned_global_particle_repository_state_v1";
    global_particle_repository_state_.lifecycle_mode =
        "global_particle_transport_reservoir_lifecycle_v1";
    global_particle_domain_state_.total_charge_c = shared_particle_transport_charge_c_;
    global_particle_domain_state_.total_reference_shift_v = shared_particle_transport_reference_shift_v_;
    global_particle_domain_state_.edge_graph_operator_iterations =
        shared_particle_transport_edge_graph_operator_iterations_last_;
    global_particle_domain_state_.edge_graph_operator_converged =
        shared_particle_transport_edge_graph_operator_converged_last_ >= 0.5;
    global_particle_domain_state_.edge_graph_operator_max_balance_residual_c =
        shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_;
    global_particle_domain_state_.edge_graph_operator_branch_graph_edge_count =
        shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_;
    global_particle_domain_state_.edge_graph_operator_branch_graph_pair_count =
        shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_;
    global_particle_domain_state_.edge_graph_operator_effective_pair_count =
        shared_particle_transport_edge_graph_operator_effective_pair_count_last_;
    global_particle_domain_state_.edge_graph_operator_total_pair_weight_f =
        shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_;
    global_particle_domain_state_.edge_graph_operator_total_conductance_weight_f =
        shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_;
    global_particle_domain_state_.edge_graph_operator_min_node_preconditioner =
        shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_;
    global_particle_domain_state_.edge_graph_operator_max_node_preconditioner =
        shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_;

    for (const auto& runtime_state : runtime_states)
    {
        if (!runtime_state.shared_runtime_enabled)
        {
            continue;
        }

        GlobalParticleDomainNodeState node_state;
        node_state.node_index = runtime_state.node_index;
        node_state.node_name = runtime_state.node_name;
        node_state.area_m2 = runtime_state.shared_patch_area_m2 > 0.0
                                 ? runtime_state.shared_patch_area_m2
                                 : runtime_state.surface_area_m2;
        node_state.patch_potential_v = runtime_state.patch_potential_v;
        node_state.shared_reference_potential_v = runtime_state.shared_reference_potential_v;
        node_state.charge_c = runtime_state.distributed_particle_transport_charge_c;
        node_state.reference_shift_v =
            runtime_state.distributed_particle_transport_reference_shift_v;
        node_state.net_flux_a = runtime_state.distributed_particle_transport_net_flux_a;
        global_particle_domain_state_.nodes.push_back(node_state);
        global_particle_domain_state_.total_node_charge_c += node_state.charge_c;
        global_particle_domain_state_.total_node_flux_a += node_state.net_flux_a;
    }

    rebuildOwnedGlobalParticleDomainEdges();

    global_particle_domain_state_.charge_conservation_error_c =
        std::abs(global_particle_domain_state_.total_node_charge_c -
                 global_particle_domain_state_.total_charge_c);
    global_particle_domain_state_.flux_conservation_error_a =
        std::abs(global_particle_domain_state_.total_node_flux_a);
    global_particle_domain_state_.active =
        global_particle_domain_state_.nodes.size() >= 2 &&
        (std::abs(global_particle_domain_state_.total_charge_c) > 0.0 ||
         !global_particle_domain_state_.edges.empty() ||
         global_particle_domain_state_.edge_charge_total_abs_c > 0.0);
    if (circuit_model_ != nullptr)
    {
        std::vector<double> node_potentials_v(circuit_model_->nodeCount(), 0.0);
        for (const auto& node_state : global_particle_domain_state_.nodes)
        {
            if (node_state.node_index < node_potentials_v.size())
            {
                node_potentials_v[node_state.node_index] = node_state.patch_potential_v;
            }
        }
        if (!global_particle_repository_state_.active)
        {
            const std::vector<double> zero_node_values(node_potentials_v.size(), 0.0);
            populateOwnedGlobalParticleRepositoryStateFromTransportBuffers(
                node_potentials_v, shared_particle_transport_node_charge_c_,
                shared_particle_transport_node_net_flux_a_, shared_particle_transport_node_charge_c_,
                zero_node_values, zero_node_values, zero_node_values,
                shared_particle_transport_edge_charge_matrix_c_,
                shared_particle_transport_edge_target_charge_matrix_c_,
                shared_particle_transport_edge_operator_drive_matrix_c_,
                shared_particle_transport_exchange_flux_matrix_a_,
                global_particle_transport_conductance_matrix_s_, 0.0);
        }
    }

    global_sheath_field_solve_state_ = GlobalSheathFieldSolveState{};
    global_sheath_field_solve_state_.solve_mode = "explicit_global_sheath_field_linear_system_v2";
    global_sheath_field_solve_state_.effective_sheath_length_m = std::max(
        1.0e-6,
        !runtime_states.empty() && runtime_states.front().shared_effective_sheath_length_m > 0.0
            ? runtime_states.front().shared_effective_sheath_length_m
            : std::max(config_.minimum_sheath_length_m,
                       std::min(config_.maximum_sheath_length_m, config_.sheath_length_m)));

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
    bool first_node = true;

    for (const auto& runtime_state : runtime_states)
    {
        if (!runtime_state.shared_runtime_enabled)
        {
            continue;
        }

        const double area_m2 = std::max(1.0e-16, runtime_state.surface_area_m2);
        GlobalSheathFieldNodeState node_state;
        node_state.node_index = runtime_state.node_index;
        node_state.node_name = runtime_state.node_name;
        node_state.area_m2 = area_m2;
        node_state.patch_potential_v = runtime_state.patch_potential_v;
        node_state.reference_potential_v = runtime_state.shared_reference_potential_v;
        node_state.normal_electric_field_v_per_m = runtime_state.normal_electric_field_v_per_m;
        node_state.local_charge_density_c_per_m3 = runtime_state.local_charge_density_c_per_m3;
        node_state.particle_charge_c = runtime_state.distributed_particle_transport_charge_c;
        node_state.particle_flux_a = runtime_state.distributed_particle_transport_net_flux_a;
        global_sheath_field_solve_state_.nodes.push_back(node_state);

        total_area_m2 += area_m2;
        weighted_patch_potential_v += area_m2 * node_state.patch_potential_v;
        weighted_reference_potential_v += area_m2 * node_state.reference_potential_v;

        if (first_node)
        {
            min_patch_potential_v = max_patch_potential_v = node_state.patch_potential_v;
            min_reference_potential_v = max_reference_potential_v = node_state.reference_potential_v;
            min_normal_field_v_per_m = max_normal_field_v_per_m =
                node_state.normal_electric_field_v_per_m;
            min_charge_density_c_per_m3 = max_charge_density_c_per_m3 =
                node_state.local_charge_density_c_per_m3;
            first_node = false;
        }
        else
        {
            min_patch_potential_v = std::min(min_patch_potential_v, node_state.patch_potential_v);
            max_patch_potential_v = std::max(max_patch_potential_v, node_state.patch_potential_v);
            min_reference_potential_v =
                std::min(min_reference_potential_v, node_state.reference_potential_v);
            max_reference_potential_v =
                std::max(max_reference_potential_v, node_state.reference_potential_v);
            min_normal_field_v_per_m =
                std::min(min_normal_field_v_per_m, node_state.normal_electric_field_v_per_m);
            max_normal_field_v_per_m =
                std::max(max_normal_field_v_per_m, node_state.normal_electric_field_v_per_m);
            min_charge_density_c_per_m3 = std::min(
                min_charge_density_c_per_m3, node_state.local_charge_density_c_per_m3);
            max_charge_density_c_per_m3 = std::max(
                max_charge_density_c_per_m3, node_state.local_charge_density_c_per_m3);
        }
    }

    if (total_area_m2 > 0.0)
    {
        global_sheath_field_solve_state_.global_patch_potential_v =
            weighted_patch_potential_v / total_area_m2;
        global_sheath_field_solve_state_.global_reference_potential_v =
            weighted_reference_potential_v / total_area_m2;
        global_sheath_field_solve_state_.global_normal_electric_field_v_per_m =
            (global_sheath_field_solve_state_.global_reference_potential_v -
             global_sheath_field_solve_state_.global_patch_potential_v) /
            global_sheath_field_solve_state_.effective_sheath_length_m;
        global_sheath_field_solve_state_.global_local_charge_density_c_per_m3 =
            kEpsilon0 * global_sheath_field_solve_state_.global_normal_electric_field_v_per_m /
            global_sheath_field_solve_state_.effective_sheath_length_m;
        const double reference_spread_v = max_reference_potential_v - min_reference_potential_v;
        const double normal_field_spread_v_per_m =
            max_normal_field_v_per_m - min_normal_field_v_per_m;
        const double local_charge_density_spread_c_per_m3 =
            max_charge_density_c_per_m3 - min_charge_density_c_per_m3;
        global_sheath_field_solve_state_.field_residual_v_per_m = std::max(
            {normal_field_spread_v_per_m,
             std::abs(reference_spread_v) /
                 global_sheath_field_solve_state_.effective_sheath_length_m,
             std::abs(local_charge_density_spread_c_per_m3) *
                 global_sheath_field_solve_state_.effective_sheath_length_m /
                 std::max(1.0e-12, kEpsilon0)});
        global_sheath_field_solve_state_.particle_field_coupled_residual_v = std::max(
            {global_sheath_field_solve_state_.field_residual_v_per_m *
                 global_sheath_field_solve_state_.effective_sheath_length_m,
             global_particle_domain_state_.charge_conservation_error_c *
                 global_sheath_field_solve_state_.effective_sheath_length_m /
                 std::max(1.0e-18, kEpsilon0 * total_area_m2),
             global_particle_domain_state_.flux_conservation_error_a *
                 global_sheath_field_solve_state_.effective_sheath_length_m /
                 std::max(1.0e-18, kEpsilon0 * total_area_m2)});
        for (auto& node_state : global_sheath_field_solve_state_.nodes)
        {
            node_state.reference_potential_v =
                global_sheath_field_solve_state_.global_reference_potential_v;
            node_state.normal_electric_field_v_per_m =
                global_sheath_field_solve_state_.global_normal_electric_field_v_per_m;
            node_state.local_charge_density_c_per_m3 =
                global_sheath_field_solve_state_.global_local_charge_density_c_per_m3;
        }
        global_sheath_field_solve_state_.field_residual_v_per_m = 0.0;
    }
    global_sheath_field_solve_state_.multi_step_stability_metric_v =
        computeSharedSurfaceGlobalSheathFieldMultiStepStabilityMetricV();
    global_sheath_field_solve_state_.active =
        useSharedSurfaceGlobalCoupledSolve() &&
        global_sheath_field_solve_state_.nodes.size() >= 2 && total_area_m2 > 0.0;
}

const GlobalParticleDomainNodeState* DensePlasmaSurfaceCharging::findGlobalParticleDomainNodeState(
    std::size_t node_index) const
{
    for (const auto& node_state : global_particle_domain_state_.nodes)
    {
        if (node_state.node_index == node_index)
        {
            return &node_state;
        }
    }
    return nullptr;
}

const GlobalSheathFieldNodeState* DensePlasmaSurfaceCharging::findGlobalSheathFieldNodeState(
    std::size_t node_index) const
{
    for (const auto& node_state : global_sheath_field_solve_state_.nodes)
    {
        if (node_state.node_index == node_index)
        {
            return &node_state;
        }
    }
    return nullptr;
}

double DensePlasmaSurfaceCharging::currentSharedSurfaceParticleTransportChargeC() const
{
    return global_particle_domain_state_.active ? global_particle_domain_state_.total_charge_c
                                                : shared_particle_transport_charge_c_;
}

double DensePlasmaSurfaceCharging::currentSharedSurfaceParticleTransportReferenceShiftV() const
{
    return global_particle_domain_state_.active
               ? global_particle_domain_state_.total_reference_shift_v
               : shared_particle_transport_reference_shift_v_;
}

SurfaceDidvState DensePlasmaSurfaceCharging::assembleSurfaceDidvState(
    ReferenceSurfaceRole role, const SurfaceModelRuntimeState& runtime_state,
    const SurfaceCurrents& currents) const
{
    SurfaceDidvState didv;
    didv.node_index = runtime_state.node_index;
    didv.node_name = runtime_state.node_name;
    didv.node_area_m2 = std::max(1.0e-16, runtime_state.surface_area_m2);

    const double conduction_didv_a_per_m2_per_v =
        role == ReferenceSurfaceRole::Patch
            ? std::max(0.0, runtime_state.effective_conductivity_s_per_m) /
                  std::max(1.0e-9, config_.dielectric_thickness_m)
            : 0.0;
    didv.plasma_current_a =
        (currents.total_current_a_per_m2 - currents.conduction_current_a_per_m2) *
        didv.node_area_m2;
    didv.conduction_current_a = currents.conduction_current_a_per_m2 * didv.node_area_m2;
    didv.conduction_didv_a_per_v = conduction_didv_a_per_m2_per_v * didv.node_area_m2;
    didv.total_current_a = currents.total_current_a_per_m2 * didv.node_area_m2;
    didv.total_didv_a_per_v =
        currents.current_derivative_a_per_m2_per_v * didv.node_area_m2;
    didv.plasma_didv_a_per_v = didv.total_didv_a_per_v - didv.conduction_didv_a_per_v;
    didv.valid = std::isfinite(didv.plasma_current_a) && std::isfinite(didv.plasma_didv_a_per_v) &&
                 std::isfinite(didv.conduction_current_a) &&
                 std::isfinite(didv.conduction_didv_a_per_v) &&
                 std::isfinite(didv.total_current_a) && std::isfinite(didv.total_didv_a_per_v);
    return didv;
}

SurfaceCircuitMappingState DensePlasmaSurfaceCharging::buildSurfaceCircuitMappingState() const
{
    SurfaceCircuitMappingState mapping;
    if (circuit_model_ == nullptr)
    {
        return mapping;
    }

    const auto node_count = circuit_model_->nodeCount();
    mapping.surface_to_circuit_node_index.resize(node_count);
    mapping.circuit_to_surface_node_indices.resize(node_count);
    mapping.circuit_to_reduced_node_index.assign(node_count, std::numeric_limits<std::size_t>::max());

    std::unordered_map<std::string, std::size_t> reduced_group_by_id;
    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        mapping.surface_to_circuit_node_index[node_index] = node_index;
        mapping.circuit_to_surface_node_indices[node_index].push_back(node_index);

        const auto node_name = circuit_model_->nodeName(node_index);
        std::string group_id = boundaryGroupIdForNode(config_, node_name);
        if (group_id.empty())
        {
            group_id = circuit_model_->nodeIsPatch(node_index) ? logicalNodeIdFromName(node_name) : node_name;
        }
        if (group_id.empty())
        {
            group_id = "node:" + std::to_string(node_index);
        }

        auto [it, inserted] =
            reduced_group_by_id.emplace(group_id, mapping.reduced_node_groups.size());
        if (inserted)
        {
            SurfaceCircuitReducedNodeGroup group;
            group.reduced_node_index = it->second;
            group.group_id = group_id;
            mapping.reduced_node_groups.push_back(group);
        }

        auto& group = mapping.reduced_node_groups[it->second];
        group.member_node_indices.push_back(node_index);
        group.total_area_m2 += std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        mapping.circuit_to_reduced_node_index[node_index] = group.reduced_node_index;
    }

    mapping.valid = true;
    return mapping;
}

double DensePlasmaSurfaceCharging::computeCapacitancePerArea() const
{
    const SurfaceModelRuntimeState state =
        buildRuntimeState(status_.body_potential_v, status_.patch_potential_v, config_.surface_area_m2,
                          0.0, computeEffectiveSheathLength());
    return computeCapacitancePerArea(state);
}

double DensePlasmaSurfaceCharging::computeCapacitancePerArea(const SurfaceModelRuntimeState& state) const
{
    const auto& patch_material = resolvePatchMaterial(config_, state.node_index, state.node_name);
    double base_capacitance = 0.0;
    if (capacitance_model_)
    {
        base_capacitance = capacitance_model_->computeCapacitancePerArea(config_);
        base_capacitance = bubble_capacitance_estimator_
                               ? bubble_capacitance_estimator_->estimateCapacitancePerArea(
                                     config_, state, base_capacitance)
                               : base_capacitance;
    }
    else
    {
        const double thickness = std::max(1.0e-8, config_.dielectric_thickness_m);
        const double eps_r = std::max(1.0, patch_material.getPermittivity());
        if (config_.derive_capacitance_from_material)
        {
            base_capacitance = kEpsilon0 * eps_r / thickness;
        }
        else
        {
            base_capacitance = std::max(1.0e-12, config_.capacitance_per_area_f_per_m2);
        }
    }

    base_capacitance *= std::clamp(
        patch_material.getScalarProperty("surface_capacitance_scale", 1.0), 0.1, 10.0);

    const bool enable_sheath_capacitance_consistency =
        config_.derive_capacitance_from_material ||
        state.external_volume_feedback_applied > 0.0 ||
        std::abs(state.field_solver_capacitance_scale - 1.0) > 1.0e-9;
    if (enable_sheath_capacitance_consistency)
    {
        Coupling::SurfaceSheathCapacitanceConsistencyInputs consistency_inputs;
        consistency_inputs.dielectric_capacitance_per_area_f_per_m2 = base_capacitance;
        consistency_inputs.effective_sheath_length_m = state.effective_sheath_length_m;
        consistency_inputs.minimum_sheath_length_m = config_.minimum_sheath_length_m;
        consistency_inputs.maximum_sheath_length_m = config_.maximum_sheath_length_m;
        consistency_inputs.relative_permittivity = config_.material.getPermittivity();
        consistency_inputs.consistency_weight = config_.material.getScalarProperty(
            "sheath_capacitance_consistency_weight", 0.60);
        consistency_inputs.volume_mesh_coupling_gain = state.volume_mesh_coupling_gain;
        consistency_inputs.ratio_guard =
            config_.material.getScalarProperty("sheath_capacitance_ratio_guard", 1.5);

        // Keep sheath length and equivalent capacitance mutually constrained to reduce dt-sensitive drift.
        base_capacitance =
            Coupling::enforceSurfaceSheathCapacitanceConsistency(consistency_inputs);
    }

    double graph_multiplier = 1.0;
    if (state.graph_capacitance_diagonal_f > 0.0 && state.graph_capacitance_row_sum_f > 0.0)
    {
        const double graph_ratio =
            state.graph_capacitance_row_sum_f /
            std::max(1.0e-30, state.graph_capacitance_diagonal_f);
        const double graph_weight = std::clamp(
            config_.material.getScalarProperty("graph_matrix_capacitance_weight", 0.15), 0.0, 0.5);
          graph_multiplier += graph_weight * std::max(0.0, graph_ratio - 1.0);
    }

    if (state.pseudo_volume_m3 > 0.0)
    {
        const double normalized_volume =
            state.pseudo_volume_m3 / std::max(1.0e-18, state.surface_area_m2 * 1.0e-3);
        const double volume_weight = std::clamp(
            config_.material.getScalarProperty("volume_mesh_capacitance_weight", 0.12), 0.0, 0.4);
        graph_multiplier +=
            volume_weight * state.volume_mesh_coupling_gain * std::max(0.0, normalized_volume - 0.5);
        graph_multiplier += 0.02 * std::max(0.0, state.volume_projection_weight_sum - 1.0);
    }

    const double solver_scale = std::clamp(state.field_solver_capacitance_scale, 0.25, 4.0);
    return std::max(1.0e-12, base_capacitance * graph_multiplier * solver_scale);
}

double DensePlasmaSurfaceCharging::computeBodyNodeCapacitanceF(
    const SurfaceModelRuntimeState& state) const
{
    double base_capacitance_f = std::max(1.0e-18, config_.body_capacitance_f);

    double graph_multiplier = 1.0;
    if (state.graph_capacitance_diagonal_f > 0.0 && state.graph_capacitance_row_sum_f > 0.0)
    {
        const double graph_ratio =
            state.graph_capacitance_row_sum_f /
            std::max(1.0e-30, state.graph_capacitance_diagonal_f);
        const double graph_weight = std::clamp(
            config_.material.getScalarProperty("graph_matrix_capacitance_weight", 0.15), 0.0, 0.5);
        graph_multiplier += graph_weight * std::max(0.0, graph_ratio - 1.0);
    }

    if (state.pseudo_volume_m3 > 0.0)
    {
        const double normalized_volume =
            state.pseudo_volume_m3 / std::max(1.0e-18, state.surface_area_m2 * 1.0e-3);
        const double volume_weight = std::clamp(
            config_.material.getScalarProperty("volume_mesh_capacitance_weight", 0.12), 0.0, 0.4);
        graph_multiplier +=
            volume_weight * state.volume_mesh_coupling_gain * std::max(0.0, normalized_volume - 0.5);
        graph_multiplier += 0.02 * std::max(0.0, state.volume_projection_weight_sum - 1.0);
    }

    const double solver_scale = std::clamp(state.field_solver_capacitance_scale, 0.25, 4.0);
    return std::max(1.0e-18, base_capacitance_f * graph_multiplier * solver_scale);
}

double DensePlasmaSurfaceCharging::computeEffectiveSheathLength() const
{
    if (!config_.derive_sheath_length_from_plasma)
    {
        return std::max(1.0e-6, config_.sheath_length_m);
    }

    const double electron_temperature_j =
        std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_density = std::max(
        1.0e3, config_.plasma.electron_density_m3 * computeTransitionSourceFluxScale());
    const double debye_length = std::sqrt(kEpsilon0 * electron_temperature_j /
                                          (electron_density * kElementaryCharge * kElementaryCharge));
    return std::clamp(debye_length, std::max(1.0e-6, config_.minimum_sheath_length_m),
                      std::max(config_.minimum_sheath_length_m, config_.maximum_sheath_length_m));
}

double DensePlasmaSurfaceCharging::computeTransitionSourceFluxScale() const
{
    return std::max(1.0e-9, transition_source_flux_scale_);
}

SurfaceCurrents DensePlasmaSurfaceCharging::computeLegacySurfaceCurrents(
    double surface_potential_v) const
{
    SurfaceCurrents currents;
    const double source_flux_scale = computeTransitionSourceFluxScale();
    const double effective_electron_density_m3 =
        std::max(0.0, config_.plasma.electron_density_m3 * source_flux_scale);
    const double effective_ion_density_m3 =
        std::max(0.0, config_.plasma.ion_density_m3 * source_flux_scale);
    const double effective_photon_flux_m2_s =
        std::max(0.0, config_.emission.photon_flux_m2_s * source_flux_scale);
    const std::size_t patch_node_index =
        circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0;
    const std::string patch_node_name =
        circuit_model_ ? circuit_model_->nodeName(patch_node_index) : std::string{"patch"};
    const auto runtime_state =
        buildRuntimeState(status_.body_potential_v, surface_potential_v, config_.surface_area_m2, 0.0,
                          computeEffectiveSheathLength(), patch_node_index, patch_node_name);
    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double bohm_velocity = std::sqrt(te_j / ion_mass);
    const auto flow_projection = resolveSurfaceFlowProjection(
        config_.bulk_flow_velocity_m_per_s, config_.flow_alignment_cosine,
        config_.patch_flow_angle_deg);
    const double normal_flow_speed =
        std::max(0.0, flow_projection.patch_projected_speed_m_per_s);
    const double directed_ion_velocity =
        std::max(0.0, config_.ion_directed_velocity_m_per_s + normal_flow_speed);
    const double effective_ion_velocity =
        std::sqrt(bohm_velocity * bohm_velocity + directed_ion_velocity * directed_ion_velocity);
    const double phi_scale = std::max(0.25, config_.plasma.electron_temperature_ev);
    double electron_flux =
        kElementaryCharge * effective_electron_density_m3 * electron_thermal_velocity *
        std::max(0.0, config_.electron_collection_coefficient);
    double ion_flux = kElementaryCharge * effective_ion_density_m3 * effective_ion_velocity *
                      std::max(0.0, config_.ion_collection_coefficient);

    if (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma)
    {
        const double ion_advective_flux =
            kElementaryCharge * effective_ion_density_m3 * normal_flow_speed *
            std::max(0.0, config_.ion_collection_coefficient);
        const double electron_advective_flux =
            kElementaryCharge * effective_electron_density_m3 * normal_flow_speed *
            std::clamp(config_.electron_flow_coupling, 0.0, 1.0);
        ion_flux += ion_advective_flux;
        electron_flux += electron_advective_flux;
    }

    if (config_.regime == SurfaceChargingRegime::GeoKineticPicLike)
    {
        electron_flux *=
            std::clamp(current_model_ ? current_model_->electronCalibrationFactor() : 1.0, 0.25, 4.0);
        ion_flux *=
            std::clamp(current_model_ ? current_model_->ionCalibrationFactor() : 1.0, 0.25, 4.0);
    }

    double raw_electron_current = 0.0;
    if (surface_potential_v <= 0.0)
    {
        raw_electron_current = -electron_flux * safeExp(surface_potential_v / phi_scale);
        currents.ion_current_a_per_m2 =
            ion_flux * std::sqrt(1.0 + std::min(50.0, -surface_potential_v / phi_scale));
    }
    else
    {
        raw_electron_current = -electron_flux * (1.0 + std::min(5.0, surface_potential_v / phi_scale));
        currents.ion_current_a_per_m2 = ion_flux * safeExp(-surface_potential_v / phi_scale);
    }

    Material::SurfaceImpactState impact;
    impact.incident_energy_ev =
        std::max(1.0e-2, config_.plasma.electron_temperature_ev + surface_potential_v);
    impact.particle_charge_coulomb = -kElementaryCharge;
    impact.surface_temperature_k = config_.emission.surface_temperature_k;
    const auto interaction = surface_interaction_.evaluate(config_.material, impact);
    const double secondary_escape_probability =
        emittedElectronEscapeProbability(config_.material, surface_potential_v,
                                         config_.material.getScalarProperty(
                                             "secondary_emission_escape_energy_ev", 2.0));
    const double photo_escape_probability =
        emittedElectronEscapeProbability(config_.material, surface_potential_v,
                                         config_.material.getScalarProperty(
                                             "photoelectron_escape_energy_ev", 1.5));
    const double thermionic_escape_probability =
        emittedElectronEscapeProbability(config_.material, surface_potential_v,
                                         std::max(1.0e-3,
                                                  config_.material.getScalarProperty(
                                                      "thermionic_escape_energy_ev",
                                                      2.0 * kBoltzmannEvPerK *
                                                          config_.emission.surface_temperature_k)));

    const double absorbed_fraction =
        impact.particle_charge_coulomb != 0.0
            ? std::clamp(std::abs(interaction.absorbed_charge_coulomb / impact.particle_charge_coulomb), 0.0,
                         1.0)
            : 1.0;
    currents.electron_current_a_per_m2 = raw_electron_current * absorbed_fraction;
    currents.secondary_emission_a_per_m2 =
        std::abs(currents.electron_current_a_per_m2) *
        std::clamp(interaction.secondary_emitted_electrons, 0.0, 2.5) *
        secondary_escape_probability;

    currents.photo_emission_a_per_m2 =
        kElementaryCharge * effective_photon_flux_m2_s *
        std::clamp(config_.material.getScalarProperty("photoelectron_yield", 0.02) *
                       config_.emission.enhancement_factor,
                   0.0, 0.5) *
        photo_escape_probability;

    const double work_function = std::max(0.1, config_.material.getWorkFunctionEv());
    const double thermal_energy_ev =
        std::max(1.0e-6, kBoltzmannEvPerK * config_.emission.surface_temperature_k);
    currents.thermionic_emission_a_per_m2 =
        kRichardsonConstant * config_.emission.surface_temperature_k * config_.emission.surface_temperature_k *
        safeExp(-work_function / thermal_energy_ev) * thermionic_escape_probability;

    const double local_field_v_per_m =
        std::abs(runtime_state.normal_electric_field_v_per_m);
    if (local_field_v_per_m > 1.0e7)
    {
        const double fn_prefactor = 1.54e-6 * local_field_v_per_m * local_field_v_per_m / work_function;
        const double fn_exponent =
            -6.83e9 * std::pow(work_function, 1.5) / std::max(1.0e6, local_field_v_per_m);
        currents.field_emission_a_per_m2 = fn_prefactor * safeExp(fn_exponent);
    }

    currents.electron_current_a_per_m2 =
        clampSigned(currents.electron_current_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.ion_current_a_per_m2 =
        clampSigned(currents.ion_current_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.secondary_emission_a_per_m2 =
        clampSigned(currents.secondary_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.photo_emission_a_per_m2 =
        clampSigned(currents.photo_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.thermionic_emission_a_per_m2 =
        clampSigned(currents.thermionic_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.field_emission_a_per_m2 =
        clampSigned(currents.field_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.total_current_a_per_m2 = currents.electron_current_a_per_m2 +
                                      currents.ion_current_a_per_m2 +
                                      currents.secondary_emission_a_per_m2 +
                                      currents.photo_emission_a_per_m2 +
                                      currents.thermionic_emission_a_per_m2 +
                                      currents.field_emission_a_per_m2;
    currents.total_current_a_per_m2 =
        clampSigned(currents.total_current_a_per_m2, config_.max_abs_current_density_a_per_m2);
    return currents;
}

SurfaceCurrents DensePlasmaSurfaceCharging::computeSurfaceCurrents(double surface_potential_v) const
{
    if (!config_.use_reference_current_balance ||
        (config_.regime != SurfaceChargingRegime::LeoFlowingPlasma &&
         config_.regime != SurfaceChargingRegime::GeoKineticPicLike) ||
        !current_model_)
    {
        return computeLegacySurfaceCurrents(surface_potential_v);
    }

    const double body_potential =
        initialized_ ? status_.body_potential_v : config_.body_initial_potential_v;
    const auto& primary_material = resolvePatchMaterial(
        config_, circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
        circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                       : std::string{"patch"});
    const auto runtime_state =
        buildRuntimeState(body_potential, surface_potential_v, config_.surface_area_m2,
                          computeEffectiveConductivity(surface_potential_v, primary_material),
                          computeEffectiveSheathLength(),
                          circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                          circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                                         : std::string{"patch"});
    SurfaceModelRuntimeState probe_state = runtime_state;
    // `computeSurfaceCurrents()` is a direct probe API for the caller-supplied potential. Keep
    // the assembled global diagnostics on the state, but do not collapse the probe back to the
    // current shared patch potential before evaluation.
    probe_state.shared_runtime_enabled = false;
    return current_model_->evaluate(ReferenceSurfaceRole::Patch, probe_state);
}

double DensePlasmaSurfaceCharging::computeRadiationInducedConductivity() const
{
    if (config_.radiation_conductivity_coefficient <= 0.0)
    {
        return 0.0;
    }

    const double source_flux_scale = computeTransitionSourceFluxScale();
    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double electron_flux =
        kElementaryCharge *
        std::max(0.0, config_.plasma.electron_density_m3 * source_flux_scale) *
        electron_thermal_velocity;
    const double incident_power_flux =
        std::max(0.0, std::abs(electron_flux) * std::max(1.0e-3, config_.plasma.electron_temperature_ev));
    return config_.radiation_conductivity_coefficient *
           std::pow(std::max(1.0e-12, incident_power_flux), config_.radiation_conductivity_exponent);
}

bool DensePlasmaSurfaceCharging::updateLivePicCalibration(double surface_potential_v)
{
    if (!current_model_)
    {
        return false;
    }
    const double effective_sheath_length = computeEffectiveSheathLength();
    if (useSharedSurfacePicRuntime() && circuit_model_ && circuit_model_->patchCount() > 0)
    {
        const auto primary_patch_node_index = circuit_model_->primaryPatchNodeIndex();
        const auto shared_runtime_state =
            buildRuntimeState(status_.body_potential_v, computeSharedSurfacePatchPotentialV(),
                              computeSharedSurfacePatchAreaM2(),
                              computeEffectiveConductivity(
                                  computeSharedSurfacePatchPotentialV(),
                                  resolvePatchMaterial(config_, primary_patch_node_index,
                                                       circuit_model_->nodeName(
                                                           primary_patch_node_index))),
                              effective_sheath_length, primary_patch_node_index,
                              circuit_model_->nodeName(primary_patch_node_index));
        return current_model_->recalibrate(shared_runtime_state,
                                           config_.pic_calibration_samples);
    }

    auto recalibrate_patch = [this, effective_sheath_length, surface_potential_v](std::size_t node_index,
                                                                                   const std::string& node_name,
                                                                                   double patch_area_m2,
                                                                                   double patch_potential_v) {
        const auto runtime_state =
            buildRuntimeState(status_.body_potential_v, surface_potential_v, patch_area_m2,
                              computeEffectiveConductivity(
                                  patch_potential_v,
                                  resolvePatchMaterial(config_, node_index, node_name)),
                              effective_sheath_length, node_index, node_name);
        return current_model_->recalibrate(runtime_state, config_.pic_calibration_samples);
    };

    bool updated = false;
    if (circuit_model_ && circuit_model_->patchCount() > 0)
    {
        for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
        {
            const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
            updated = recalibrate_patch(node_index, circuit_model_->nodeName(node_index),
                                        std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index)),
                                        circuit_model_->nodePotential(node_index)) || updated;
        }
        return updated;
    }

    return recalibrate_patch(circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                             circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                                            : std::string{"patch"},
                             config_.surface_area_m2, surface_potential_v);
}

void DensePlasmaSurfaceCharging::applyGeoPicCalibration()
{
    if (!config_.enable_pic_calibration)
    {
        return;
    }

    const bool kinetic_regime = config_.regime == SurfaceChargingRegime::GeoKineticPicLike ||
                                config_.regime == SurfaceChargingRegime::LeoFlowingPlasma;
    if (!kinetic_regime)
    {
        return;
    }

    const double current_potential =
        initialized_ ? status_.patch_potential_v : config_.body_initial_potential_v;
    const double floating_reference = computeFloatingPotential();
    double reference_potential = current_potential;
    if (config_.regime == SurfaceChargingRegime::GeoKineticPicLike)
    {
        // GEO calibration should sample closer to the negative charging operating
        // point, while still avoiding extreme one-shot jumps in the live PIC window.
        const double blended =
            current_potential + 0.15 * (floating_reference - current_potential);
        reference_potential = std::clamp(blended, -4.0e3, 5.0e2);
    }

    if (updateLivePicCalibration(reference_potential))
    {
        return;
    }

    const std::size_t samples = std::max<std::size_t>(256, config_.pic_calibration_samples);
    const double base_scale = std::max(25.0, 0.2 * config_.plasma.electron_temperature_ev);
    const std::array<double, 3> calibration_window{
        floating_reference - base_scale,
        floating_reference,
        floating_reference + base_scale};

    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double bohm_velocity = std::sqrt(te_j / ion_mass);

    double electron_ratio_sum = 0.0;
    double ion_ratio_sum = 0.0;
    std::size_t electron_ratio_count = 0;
    std::size_t ion_ratio_count = 0;

    for (const double potential_v : calibration_window)
    {
        const double pic_electron_flux = estimateElectronFluxPicLike(potential_v, samples);
        const double pic_ion_flux = estimateIonFluxPicLike(potential_v, samples);
        const double analytic_electron_flux =
            kElementaryCharge *
            std::max(0.0, config_.plasma.electron_density_m3 * computeTransitionSourceFluxScale()) *
            electron_thermal_velocity *
            std::max(0.0, config_.electron_collection_coefficient);
        const double analytic_ion_flux =
            kElementaryCharge *
            std::max(0.0, config_.plasma.ion_density_m3 * computeTransitionSourceFluxScale()) *
            bohm_velocity *
            std::max(0.0, config_.ion_collection_coefficient);

        if (analytic_electron_flux > 1.0e-18)
        {
            const double ratio = std::clamp(pic_electron_flux / analytic_electron_flux, 0.4, 2.5);
            electron_ratio_sum += ratio;
            electron_ratio_count += 1;
        }
        if (analytic_ion_flux > 1.0e-18)
        {
            const double ratio = std::clamp(pic_ion_flux / analytic_ion_flux, 0.4, 2.5);
            ion_ratio_sum += ratio;
            ion_ratio_count += 1;
        }
    }

    (void) electron_ratio_sum;
    (void) ion_ratio_sum;
    (void) electron_ratio_count;
    (void) ion_ratio_count;
}

double DensePlasmaSurfaceCharging::estimateElectronFluxPicLike(double surface_potential_v,
                                                               std::size_t samples) const
{
    const std::size_t sample_count = std::max<std::size_t>(64, samples);
    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double sigma = std::sqrt(te_j / kElectronMass);
    std::mt19937_64 generator(static_cast<std::uint64_t>(config_.seed) + 0x5cdab51ULL +
                              static_cast<std::uint64_t>(sample_count));

    double collected_surface_velocity_sum = 0.0;
    for (std::size_t index = 0; index < sample_count; ++index)
    {
        const double normal_velocity = halfNormalFluxAverage(sigma, generator);
        double normal_energy_j = 0.5 * kElectronMass * normal_velocity * normal_velocity;

        if (surface_potential_v < 0.0)
        {
            const double barrier_j = -kElementaryCharge * surface_potential_v;
            if (normal_energy_j <= barrier_j)
            {
                continue;
            }
            normal_energy_j -= barrier_j;
        }
        else
        {
            normal_energy_j += kElementaryCharge * surface_potential_v;
        }

        const double surface_velocity =
            std::sqrt(std::max(0.0, 2.0 * normal_energy_j / kElectronMass));
        collected_surface_velocity_sum += surface_velocity;
    }

    const double mean_surface_velocity =
        collected_surface_velocity_sum / static_cast<double>(sample_count);
    return 0.5 * kElementaryCharge *
           std::max(0.0,
                    config_.plasma.electron_density_m3 * computeTransitionSourceFluxScale()) *
           mean_surface_velocity *
           std::max(0.0, config_.electron_collection_coefficient);
}

double DensePlasmaSurfaceCharging::estimateIonFluxPicLike(double surface_potential_v,
                                                          std::size_t samples) const
{
    const std::size_t sample_count = std::max<std::size_t>(64, samples);
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double ti_j = std::max(1.0e-3, config_.plasma.ion_temperature_ev) * kElementaryCharge;
    const double sigma = std::sqrt(ti_j / ion_mass);
    const auto flow_projection = resolveSurfaceFlowProjection(
        config_.bulk_flow_velocity_m_per_s, config_.flow_alignment_cosine,
        config_.patch_flow_angle_deg);
    const double drift_velocity =
        std::max(0.0, config_.ion_directed_velocity_m_per_s +
                          std::max(0.0, flow_projection.patch_projected_speed_m_per_s));
    std::mt19937_64 generator(static_cast<std::uint64_t>(config_.seed) + 0x71015aULL +
                              static_cast<std::uint64_t>(sample_count));

    double collected_surface_velocity_sum = 0.0;
    for (std::size_t index = 0; index < sample_count; ++index)
    {
        const double launch_velocity = halfNormalFluxAverage(sigma, generator) + drift_velocity;
        double normal_energy_j = 0.5 * ion_mass * launch_velocity * launch_velocity;

        if (surface_potential_v > 0.0)
        {
            const double barrier_j = kElementaryCharge * surface_potential_v;
            if (normal_energy_j <= barrier_j)
            {
                continue;
            }
            normal_energy_j -= barrier_j;
        }
        else
        {
            normal_energy_j += -kElementaryCharge * surface_potential_v;
        }

        const double surface_velocity = std::sqrt(std::max(0.0, 2.0 * normal_energy_j / ion_mass));
        collected_surface_velocity_sum += surface_velocity;
    }

    const double mean_surface_velocity =
        collected_surface_velocity_sum / static_cast<double>(sample_count);
    return 0.5 * kElementaryCharge *
           std::max(0.0, config_.plasma.ion_density_m3 * computeTransitionSourceFluxScale()) *
           mean_surface_velocity *
           std::max(0.0, config_.ion_collection_coefficient);
}

double DensePlasmaSurfaceCharging::computeEffectiveConductivity(double surface_potential_v) const
{
    return computeEffectiveConductivity(surface_potential_v, config_.material);
}

double DensePlasmaSurfaceCharging::computeEffectiveConductivity(
    double surface_potential_v, const Material::MaterialProperty& material) const
{
    const double base_conductivity = std::max(
        0.0,
        std::max(material.getConductivity(), material.getDarkConductivitySPerM()) *
            std::max(1.0e-9, transition_conductivity_scale_));
    const std::size_t patch_node_index =
        circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0;
    const std::string patch_node_name =
        circuit_model_ ? circuit_model_->nodeName(patch_node_index) : std::string{"patch"};
    const auto runtime_state =
        buildRuntimeState(status_.body_potential_v, surface_potential_v, config_.surface_area_m2,
                          base_conductivity, computeEffectiveSheathLength(), patch_node_index,
                          patch_node_name);
    return Material::evaluateSurfaceEffectiveConductivitySPerM(
        material, base_conductivity, computeRadiationInducedConductivity(), surface_potential_v,
        status_.body_potential_v, runtime_state.normal_electric_field_v_per_m,
        runtime_state.local_charge_density_c_per_m3, config_.dielectric_thickness_m);
}

bool DensePlasmaSurfaceCharging::useSharedSurfacePicRuntime() const
{
    return config_.surface_pic_runtime_kind ==
           SurfacePicRuntimeKind::GraphCoupledSharedSurface;
}

double DensePlasmaSurfaceCharging::computeSharedSurfacePatchPotentialV() const
{
    if (!useSharedSurfacePicRuntime())
    {
        return status_.patch_potential_v;
    }

    if (global_sheath_field_solve_state_.active && !global_sheath_field_solve_state_.nodes.empty())
    {
        double weighted_sum_v = 0.0;
        double weight_sum = 0.0;
        for (const auto& node_state : global_sheath_field_solve_state_.nodes)
        {
            const double area_m2 = std::max(1.0e-16, node_state.area_m2);
            weighted_sum_v += area_m2 * node_state.patch_potential_v;
            weight_sum += area_m2;
        }
        if (weight_sum > 0.0)
        {
            return weighted_sum_v / weight_sum;
        }
    }
    if (global_particle_domain_state_.active && !global_particle_domain_state_.nodes.empty())
    {
        double weighted_sum_v = 0.0;
        double weight_sum = 0.0;
        for (const auto& node_state : global_particle_domain_state_.nodes)
        {
            const double area_m2 = std::max(1.0e-16, node_state.area_m2);
            weighted_sum_v += area_m2 * node_state.patch_potential_v;
            weight_sum += area_m2;
        }
        if (weight_sum > 0.0)
        {
            return weighted_sum_v / weight_sum;
        }
    }

    if (!circuit_model_ || circuit_model_->patchCount() == 0 ||
        circuit_node_potentials_.empty())
    {
        return status_.patch_potential_v;
    }

    double weighted_sum_v = 0.0;
    double weight_sum = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        if (node_index >= circuit_node_potentials_.size())
        {
            continue;
        }
        const double area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        weighted_sum_v += area_m2 * circuit_node_potentials_[node_index];
        weight_sum += area_m2;
    }
    if (weight_sum <= 0.0)
    {
        return status_.patch_potential_v;
    }
    return weighted_sum_v / weight_sum;
}

double DensePlasmaSurfaceCharging::computeSharedSurfacePatchAreaM2() const
{
    if (!useSharedSurfacePicRuntime())
    {
        return config_.surface_area_m2;
    }

    if (global_sheath_field_solve_state_.active && !global_sheath_field_solve_state_.nodes.empty())
    {
        double total_area_m2 = 0.0;
        for (const auto& node_state : global_sheath_field_solve_state_.nodes)
        {
            total_area_m2 += std::max(1.0e-16, node_state.area_m2);
        }
        if (total_area_m2 > 0.0)
        {
            return total_area_m2;
        }
    }
    if (global_particle_domain_state_.active && !global_particle_domain_state_.nodes.empty())
    {
        double total_area_m2 = 0.0;
        for (const auto& node_state : global_particle_domain_state_.nodes)
        {
            total_area_m2 += std::max(1.0e-16, node_state.area_m2);
        }
        if (total_area_m2 > 0.0)
        {
            return total_area_m2;
        }
    }

    if (!circuit_model_ || circuit_model_->patchCount() == 0)
    {
        return config_.surface_area_m2;
    }

    double total_area_m2 = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        total_area_m2 += std::max(
            1.0e-16, circuit_model_->nodeAreaM2(circuit_model_->patchNodeIndex(patch_ordinal)));
    }
    return std::max(1.0e-16, total_area_m2);
}

double DensePlasmaSurfaceCharging::computeSharedSurfacePatchPotentialSpreadV() const
{
    if (!useSharedSurfacePicRuntime())
    {
        return 0.0;
    }

    const double shared_patch_potential_v = computeSharedSurfacePatchPotentialV();
    if (global_sheath_field_solve_state_.active && !global_sheath_field_solve_state_.nodes.empty())
    {
        double weighted_spread_v = 0.0;
        double weight_sum = 0.0;
        for (const auto& node_state : global_sheath_field_solve_state_.nodes)
        {
            const double area_m2 = std::max(1.0e-16, node_state.area_m2);
            weighted_spread_v +=
                area_m2 * std::abs(node_state.patch_potential_v - shared_patch_potential_v);
            weight_sum += area_m2;
        }
        if (weight_sum > 0.0)
        {
            return weighted_spread_v / weight_sum;
        }
    }
    if (global_particle_domain_state_.active && !global_particle_domain_state_.nodes.empty())
    {
        double weighted_spread_v = 0.0;
        double weight_sum = 0.0;
        for (const auto& node_state : global_particle_domain_state_.nodes)
        {
            const double area_m2 = std::max(1.0e-16, node_state.area_m2);
            weighted_spread_v +=
                area_m2 * std::abs(node_state.patch_potential_v - shared_patch_potential_v);
            weight_sum += area_m2;
        }
        if (weight_sum > 0.0)
        {
            return weighted_spread_v / weight_sum;
        }
    }

    if (!circuit_model_ || circuit_model_->patchCount() == 0 || circuit_node_potentials_.empty())
    {
        return 0.0;
    }

    double weighted_spread_v = 0.0;
    double weight_sum = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        if (node_index >= circuit_node_potentials_.size())
        {
            continue;
        }
        const double area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        weighted_spread_v +=
            area_m2 * std::abs(circuit_node_potentials_[node_index] - shared_patch_potential_v);
        weight_sum += area_m2;
    }

    if (weight_sum <= 0.0)
    {
        return 0.0;
    }
    return weighted_spread_v / weight_sum;
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceEffectiveSheathLengthM(
    double base_sheath_length_m) const
{
    const double min_sheath_length_m = std::max(1.0e-6, config_.minimum_sheath_length_m);
    const double max_sheath_length_m =
        std::max(min_sheath_length_m, config_.maximum_sheath_length_m);
    if (!useSharedSurfacePicRuntime())
    {
        return std::clamp(base_sheath_length_m, min_sheath_length_m, max_sheath_length_m);
    }

    const double shared_patch_area_m2 = computeSharedSurfacePatchAreaM2();
    const double potential_spread_v = computeSharedSurfacePatchPotentialSpreadV();
    const double temperature_scale_ev = std::max(0.5, config_.plasma.electron_temperature_ev);
    const double spread_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_sheath_spread_weight", 0.18), 0.0,
        1.5);
    const double area_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_sheath_area_weight", 0.05), 0.0, 0.5);
    const double area_ratio =
        shared_patch_area_m2 / std::max(1.0e-16, config_.surface_area_m2);
    const double spread_multiplier =
        1.0 + spread_weight * std::clamp(potential_spread_v / temperature_scale_ev, 0.0, 6.0);
    const double area_multiplier =
        1.0 + area_weight * std::clamp(std::log10(std::max(1.0, area_ratio)), 0.0, 4.0);
    return std::clamp(base_sheath_length_m * spread_multiplier * area_multiplier,
                      min_sheath_length_m, max_sheath_length_m);
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceReferencePotentialV(
    double base_reference_potential_v, double local_patch_potential_v) const
{
    (void)base_reference_potential_v;
    (void)local_patch_potential_v;
    if (!useSharedSurfacePicRuntime())
    {
        return base_reference_potential_v;
    }

    const double shared_patch_potential_v = computeSharedSurfacePatchPotentialV();
    const double spread_v = computeSharedSurfacePatchPotentialSpreadV();
    const double temperature_scale_ev = std::max(0.5, config_.plasma.electron_temperature_ev);
    const double base_blend = std::clamp(
        config_.material.getScalarProperty("shared_surface_reference_blend", 0.65), 0.0, 1.0);
    const double spread_gain =
        1.0 + 0.10 * std::clamp(spread_v / temperature_scale_ev, 0.0, 4.0);
    const double effective_blend = std::clamp(base_blend * spread_gain, 0.0, 1.0);
    const double shared_anchor_reference_v =
        std::isfinite(config_.live_pic_reference_potential_v)
            ? config_.live_pic_reference_potential_v
            : 0.0;
    return (1.0 - effective_blend) * shared_anchor_reference_v +
           effective_blend * shared_patch_potential_v;
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalSolveWeight() const
{
    if (!useSharedSurfacePicRuntime() || !circuit_model_ || circuit_model_->patchCount() < 2)
    {
        return 0.0;
    }
    return std::clamp(
        config_.material.getScalarProperty("shared_surface_global_solve_weight", 0.0), 0.0, 1.0);
}

bool DensePlasmaSurfaceCharging::useSharedSurfaceGlobalCoupledSolve() const
{
    return useSharedSurfacePicRuntime() && circuit_model_ && circuit_model_->patchCount() >= 2 &&
           computeSharedSurfaceGlobalCoupledIterationLimit() > 1;
}

std::size_t DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalCoupledIterationLimit() const
{
    if (!useSharedSurfacePicRuntime() || !circuit_model_ || circuit_model_->patchCount() < 2)
    {
        return 0;
    }

    const double configured_iterations = std::clamp(
        config_.material.getScalarProperty("shared_surface_global_coupled_iterations", 0.0), 0.0,
        256.0);
    const std::size_t configured_limit =
        static_cast<std::size_t>(std::llround(configured_iterations));
    return resolveSharedGlobalCoupledIterationLimit(config_, configured_limit);
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalCoupledToleranceV() const
{
    const double configured_tolerance_v = std::max(
        1.0e-9,
        config_.material.getScalarProperty("shared_surface_global_coupled_tolerance_v", 1.0e-6));
    return resolveSharedGlobalCoupledToleranceV(config_, configured_tolerance_v);
}

bool DensePlasmaSurfaceCharging::useSharedSurfaceLivePicCoupledRefresh() const
{
    return useSharedSurfaceGlobalCoupledSolve() &&
           (config_.enable_live_pic_window || config_.enable_pic_calibration) &&
           config_.material.getScalarProperty("shared_surface_live_pic_coupled_refresh", 0.0) >
               0.5;
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceLivePicCoupledRefreshThresholdV() const
{
    return std::max(
        0.0,
        config_.material.getScalarProperty(
            "shared_surface_live_pic_coupled_refresh_threshold_v",
            std::max(0.0, 0.25 * config_.live_pic_probe_delta_v)));
}

bool DensePlasmaSurfaceCharging::useSharedSurfaceParticleTransportCoupling() const
{
    return useSharedSurfaceGlobalCoupledSolve() && circuit_model_ && circuit_model_->patchCount() >= 2 &&
           config_.material.getScalarProperty("shared_surface_particle_transport_weight", 0.0) >
               0.0;
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceParticleTransportConductanceSPerM2(
    double shared_live_pic_net_current_density_a_per_m2,
    double shared_live_pic_derivative_a_per_m2_per_v) const
{
    if (!useSharedSurfaceParticleTransportCoupling())
    {
        return 0.0;
    }

    const double transport_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_weight", 0.0), 0.0,
        4.0);
    const double net_current_scale =
        std::abs(shared_live_pic_net_current_density_a_per_m2) /
        std::max(1.0e-3, std::abs(config_.live_pic_probe_delta_v));
    const double base_conductance_s_per_m2 =
        std::max(std::abs(shared_live_pic_derivative_a_per_m2_per_v), net_current_scale);
    const double max_conductance_s_per_m2 = std::max(
        0.0, config_.material.getScalarProperty(
                 "shared_surface_particle_transport_max_conductance_s_per_m2", 1.0e3));
    return std::clamp(transport_weight * base_conductance_s_per_m2, 0.0,
                      max_conductance_s_per_m2);
}

void DensePlasmaSurfaceCharging::advanceSharedSurfaceParticleTransportState(
    double dt, const PicMccCurrentSample& shared_live_pic_sample, double shared_patch_area_m2,
    double shared_effective_sheath_length_m)
{
    if (!useSharedSurfaceParticleTransportCoupling())
    {
        shared_particle_transport_charge_c_ = 0.0;
        shared_particle_transport_reference_shift_v_ = 0.0;
        return;
    }

    const double area_m2 = std::max(1.0e-16, shared_patch_area_m2);
    const double source_current_a =
        shared_live_pic_sample.valid
            ? shared_live_pic_sample.net_collection_current_density_a_per_m2 * area_m2
            : 0.0;
    const double relaxation_time_s = std::max(
        1.0e-9, config_.material.getScalarProperty(
                    "shared_surface_particle_transport_relaxation_time_s",
                    std::max(dt, 4.0 * dt)));
    const double damping_factor = std::exp(-std::max(0.0, dt) / relaxation_time_s);
    shared_particle_transport_charge_c_ =
        damping_factor * shared_particle_transport_charge_c_ + source_current_a * dt;

    const double max_abs_charge_c = std::max(
        0.0, config_.material.getScalarProperty(
                 "shared_surface_particle_transport_max_abs_charge_c", 1.0e-6));
    if (max_abs_charge_c > 0.0)
    {
        shared_particle_transport_charge_c_ = std::clamp(shared_particle_transport_charge_c_,
                                                         -max_abs_charge_c, max_abs_charge_c);
    }

    shared_particle_transport_reference_shift_v_ =
        computeSharedSurfaceParticleTransportReferenceShiftV(
            area_m2, shared_effective_sheath_length_m);
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceParticleTransportReferenceShiftV(
    double shared_patch_area_m2, double shared_effective_sheath_length_m) const
{
    if (!useSharedSurfaceParticleTransportCoupling())
    {
        return 0.0;
    }

    const double area_m2 = std::max(1.0e-16, shared_patch_area_m2);
    const double sheath_length_m = std::max(1.0e-6, shared_effective_sheath_length_m);
    const double coupling_scale = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_reference_scale", 1.0),
        0.0, 4.0);
    const double effective_capacitance_f =
        std::max(1.0e-18, coupling_scale * kEpsilon0 * area_m2 / sheath_length_m);
    const double raw_shift_v = shared_particle_transport_charge_c_ / effective_capacitance_f;
    const double max_abs_shift_v = std::max(
        0.0, config_.material.getScalarProperty(
                 "shared_surface_particle_transport_max_abs_reference_shift_v", 10.0));
    return std::clamp(raw_shift_v, -max_abs_shift_v, max_abs_shift_v);
}

void DensePlasmaSurfaceCharging::updateSharedSurfaceDistributedParticleTransportState(
    const std::vector<double>& node_potentials_v, double dt, double shared_patch_area_m2,
    double shared_effective_sheath_length_m,
    double shared_transport_conductance_s_per_m2)
{
    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        global_particle_domain_state_ = GlobalParticleDomainState{};
        global_particle_domain_state_.bookkeeping_mode = "owned_global_particle_domain_state_v2";
        global_particle_repository_state_ = GlobalParticleRepositoryState{};
        global_particle_repository_state_.bookkeeping_mode =
            "owned_global_particle_repository_state_v1";
        global_particle_repository_state_.lifecycle_mode =
            "global_particle_transport_reservoir_lifecycle_v1";
        shared_particle_transport_node_charge_c_.clear();
        shared_particle_transport_node_net_flux_a_.clear();
        shared_particle_transport_edge_charge_matrix_c_.clear();
        shared_particle_transport_edge_target_charge_matrix_c_.clear();
        shared_particle_transport_edge_operator_drive_matrix_c_.clear();
        shared_particle_transport_exchange_flux_matrix_a_.clear();
        global_particle_transport_conductance_matrix_s_.clear();
        shared_particle_transport_edge_graph_operator_iterations_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_converged_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_effective_pair_count_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_ = 0.0;
        shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_ = 1.0;
        shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_ = 1.0;
        return;
    }

    if (shared_particle_transport_node_charge_c_.size() != circuit_model_->nodeCount())
    {
        shared_particle_transport_node_charge_c_.assign(circuit_model_->nodeCount(), 0.0);
    }
    if (shared_particle_transport_node_net_flux_a_.size() != circuit_model_->nodeCount())
    {
        shared_particle_transport_node_net_flux_a_.assign(circuit_model_->nodeCount(), 0.0);
    }
    if (shared_particle_transport_edge_charge_matrix_c_.size() != circuit_model_->nodeCount())
    {
        shared_particle_transport_edge_charge_matrix_c_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    }
    if (shared_particle_transport_edge_target_charge_matrix_c_.size() != circuit_model_->nodeCount())
    {
        shared_particle_transport_edge_target_charge_matrix_c_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    }
    if (shared_particle_transport_edge_operator_drive_matrix_c_.size() !=
        circuit_model_->nodeCount())
    {
        shared_particle_transport_edge_operator_drive_matrix_c_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    }
    if (shared_particle_transport_exchange_flux_matrix_a_.size() != circuit_model_->nodeCount())
    {
        shared_particle_transport_exchange_flux_matrix_a_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    }
    if (global_particle_transport_conductance_matrix_s_.size() != circuit_model_->nodeCount())
    {
        global_particle_transport_conductance_matrix_s_.assign(
            circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    }

    bool should_hydrate_shared_transport_caches = false;
    if (global_particle_domain_state_.active)
    {
        const auto is_effectively_zero = [](double value) {
            return !std::isfinite(value) || std::abs(value) <= 1.0e-30;
        };
        double cache_node_charge_total_abs_c = 0.0;
        for (const double node_charge_c : shared_particle_transport_node_charge_c_)
        {
            cache_node_charge_total_abs_c += std::abs(node_charge_c);
        }
        double cache_edge_charge_total_abs_c = 0.0;
        for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
        {
            const auto from_index = circuit_model_->patchNodeIndex(patch_i);
            if (from_index >= shared_particle_transport_edge_charge_matrix_c_.size())
            {
                continue;
            }
            for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
            {
                const auto to_index = circuit_model_->patchNodeIndex(patch_j);
                if (to_index >= shared_particle_transport_edge_charge_matrix_c_[from_index].size())
                {
                    continue;
                }
                cache_edge_charge_total_abs_c += std::abs(
                    shared_particle_transport_edge_charge_matrix_c_[from_index][to_index]);
            }
        }
        should_hydrate_shared_transport_caches =
            (cache_node_charge_total_abs_c <= 1.0e-30 &&
             std::abs(global_particle_domain_state_.total_node_charge_c) > 1.0e-30) ||
            (cache_edge_charge_total_abs_c <= 1.0e-30 &&
             global_particle_domain_state_.edge_charge_total_abs_c > 1.0e-30) ||
            (is_effectively_zero(shared_particle_transport_charge_c_) &&
             std::abs(global_particle_domain_state_.total_charge_c) > 1.0e-30) ||
            (is_effectively_zero(shared_particle_transport_reference_shift_v_) &&
             std::abs(global_particle_domain_state_.total_reference_shift_v) > 1.0e-30);
    }
    if (should_hydrate_shared_transport_caches)
    {
        syncSharedTransportCachesFromOwnedGlobalParticleDomainState();
    }

    std::vector<double> transport_node_charge_c = shared_particle_transport_node_charge_c_;
    std::vector<double> transport_node_net_flux_a(circuit_model_->nodeCount(), 0.0);
    std::vector<std::vector<double>> transport_edge_charge_matrix_c =
        shared_particle_transport_edge_charge_matrix_c_;
    std::vector<std::vector<double>> transport_edge_target_charge_matrix_c(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    std::vector<std::vector<double>> transport_edge_operator_drive_matrix_c(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    std::vector<std::vector<double>> transport_exchange_flux_matrix_a(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    std::vector<std::vector<double>> transport_conductance_matrix_s(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    if (global_particle_domain_state_.active && !global_particle_domain_state_.nodes.empty())
    {
        std::fill(transport_node_charge_c.begin(), transport_node_charge_c.end(), 0.0);
        for (const auto& node_state : global_particle_domain_state_.nodes)
        {
            if (node_state.node_index < transport_node_charge_c.size())
            {
                transport_node_charge_c[node_state.node_index] = node_state.charge_c;
            }
        }
        for (auto& row : transport_edge_charge_matrix_c)
        {
            std::fill(row.begin(), row.end(), 0.0);
        }
        for (const auto& edge_state : global_particle_domain_state_.edges)
        {
            if (edge_state.from_node_index >= transport_edge_charge_matrix_c.size() ||
                edge_state.to_node_index >=
                    transport_edge_charge_matrix_c[edge_state.from_node_index].size())
            {
                continue;
            }
            transport_edge_charge_matrix_c[edge_state.from_node_index][edge_state.to_node_index] =
                edge_state.stored_charge_c;
            transport_edge_charge_matrix_c[edge_state.to_node_index][edge_state.from_node_index] =
                -edge_state.stored_charge_c;
        }
    }

    const double total_area_m2 = std::max(1.0e-16, shared_patch_area_m2);
    const double shared_patch_potential_v = computeSharedSurfacePatchPotentialV();
    const double distribution_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_distribution_weight",
                                           0.35),
        0.0, 1.0);
    const double relaxation = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_distribution_relaxation",
                                           0.45),
        0.0, 1.0);
    const double exchange_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_exchange_weight", 0.20),
        0.0, 2.0);
    const double exchange_relaxation_time_s = std::max(
        1.0e-9, config_.material.getScalarProperty(
                    "shared_surface_particle_transport_exchange_relaxation_time_s",
                    std::max(dt, 4.0 * dt)));
    const double exchange_rate_hz = exchange_weight / exchange_relaxation_time_s;
    const double edge_memory_relaxation_time_s = std::max(
        1.0e-9, config_.material.getScalarProperty(
                    "shared_surface_particle_transport_edge_memory_relaxation_time_s",
                    std::max(dt, 8.0 * dt)));
    const double edge_memory_damping_factor = std::exp(-std::max(0.0, dt) /
                                                       edge_memory_relaxation_time_s);
    const double edge_feedback_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_edge_feedback_weight",
                                           0.0),
        0.0, 1.0);
    const double edge_operator_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_edge_operator_weight",
                                           0.0),
        0.0, 1.0);
    const double edge_operator_relaxation = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_operator_relaxation", 0.35),
        0.0, 1.0);
    const double edge_operator_capacitance_scale = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_operator_capacitance_scale", 1.0),
        0.0, 4.0);
    const std::size_t edge_graph_operator_iteration_limit = static_cast<std::size_t>(std::llround(
        std::clamp(config_.material.getScalarProperty(
                       "shared_surface_particle_transport_edge_graph_operator_iterations", 6.0),
                   0.0, 16.0)));
    const double edge_graph_operator_tolerance_c = std::max(
        1.0e-18,
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_graph_operator_tolerance_c", 1.0e-12));
    const double edge_graph_operator_potential_blend = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_graph_operator_potential_blend", 0.20),
        0.0, 1.0);
    const double edge_graph_operator_relaxation = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_graph_operator_relaxation", 0.5),
        0.05, 1.0);
    const double edge_graph_operator_path_blend = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_graph_operator_path_blend", 0.05),
        0.0, 1.0);
    const double edge_graph_operator_conductance_weight = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_graph_operator_conductance_weight", 0.002),
        0.0, 10.0);
    const double edge_graph_operator_node_preconditioner_weight = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_edge_graph_operator_node_preconditioner_weight",
            0.20),
        0.0, 1.0);
    const double max_abs_edge_charge_c = std::max(
        0.0, config_.material.getScalarProperty(
                 "shared_surface_particle_transport_max_abs_edge_charge_c",
                 std::max(1.0e-12, 0.5 * std::abs(shared_particle_transport_charge_c_))));

    std::vector<double> target_node_charge_c(circuit_model_->nodeCount(), 0.0);
    std::vector<double> signed_weights(circuit_model_->nodeCount(), 0.0);
    std::vector<double> node_exchange_charge_delta_c(circuit_model_->nodeCount(), 0.0);
    std::vector<double> node_edge_domain_charge_c(circuit_model_->nodeCount(), 0.0);
    std::vector<double> node_edge_feedback_charge_c(circuit_model_->nodeCount(), 0.0);
    std::vector<double> node_conservation_correction_charge_c(circuit_model_->nodeCount(), 0.0);
    double signed_weight_l1 = 0.0;
    std::fill(transport_node_net_flux_a.begin(), transport_node_net_flux_a.end(), 0.0);
    for (auto& row : transport_edge_target_charge_matrix_c)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    for (auto& row : transport_edge_operator_drive_matrix_c)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    for (auto& row : transport_exchange_flux_matrix_a)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    for (auto& row : transport_conductance_matrix_s)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        const double area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        const double area_share = area_m2 / total_area_m2;
        target_node_charge_c[node_index] = shared_particle_transport_charge_c_ * area_share;

        const double node_potential_v =
            node_index < node_potentials_v.size() ? node_potentials_v[node_index] : shared_patch_potential_v;
        const double centered_weight = area_m2 * (node_potential_v - shared_patch_potential_v);
        signed_weights[node_index] = centered_weight;
        signed_weight_l1 += std::abs(centered_weight);
    }

    for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
    {
        const auto node_i = circuit_model_->patchNodeIndex(patch_i);
        for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
        {
            const auto node_j = circuit_model_->patchNodeIndex(patch_j);
            double stored_edge_charge_c =
                edge_memory_damping_factor * transport_edge_charge_matrix_c[node_i][node_j];
            if (max_abs_edge_charge_c > 0.0)
            {
                stored_edge_charge_c = std::clamp(stored_edge_charge_c, -max_abs_edge_charge_c,
                                                  max_abs_edge_charge_c);
            }
            transport_edge_charge_matrix_c[node_i][node_j] = stored_edge_charge_c;
            transport_edge_charge_matrix_c[node_j][node_i] = -stored_edge_charge_c;
        }
    }

    if (signed_weight_l1 > 0.0 && std::abs(shared_particle_transport_charge_c_) > 0.0)
    {
        for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
        {
            const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
            target_node_charge_c[node_index] +=
                distribution_weight * shared_particle_transport_charge_c_ *
                (signed_weights[node_index] / signed_weight_l1);
        }
    }

    std::vector<double> desired_node_edge_balance_c(circuit_model_->nodeCount(), 0.0);
    std::vector<std::size_t> patch_nodes;
    patch_nodes.reserve(circuit_model_->patchCount());
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        patch_nodes.push_back(node_index);
        desired_node_edge_balance_c[node_index] =
            target_node_charge_c[node_index] - transport_node_charge_c[node_index];
    }
    double desired_balance_sum_c = 0.0;
    for (const auto node_index : patch_nodes)
    {
        desired_balance_sum_c += desired_node_edge_balance_c[node_index];
    }
    if (!patch_nodes.empty())
    {
        const double correction_per_node_c =
            desired_balance_sum_c / static_cast<double>(patch_nodes.size());
        for (const auto node_index : patch_nodes)
        {
            desired_node_edge_balance_c[node_index] -= correction_per_node_c;
        }
    }

    std::vector<double> graph_operator_lambda_c(circuit_model_->nodeCount(), 0.0);
    std::vector<double> graph_operator_next_lambda_c(circuit_model_->nodeCount(), 0.0);
    std::vector<double> graph_operator_weight_sum(circuit_model_->nodeCount(), 0.0);
    std::vector<double> graph_operator_node_preconditioner(
        circuit_model_->nodeCount(), 1.0);
    std::vector<std::vector<double>> branch_graph_adjacency_weight_f(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    std::vector<double> transport_graph_node_conductance_weight_f(
        circuit_model_->nodeCount(), 0.0);
    std::vector<std::vector<double>> graph_operator_pair_weight_f(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
    std::size_t branch_graph_edge_count = 0;
    std::size_t branch_graph_pair_count = 0;
    std::size_t effective_pair_count = 0;
    double total_pair_weight_f = 0.0;
    double total_conductance_weight_f = 0.0;
    double min_node_preconditioner = 1.0;
    double max_node_preconditioner = 1.0;
    if (circuit_model_ != nullptr)
    {
        for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount(); ++branch_index)
        {
            const auto from_node = circuit_model_->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model_->branchToNodeIndex(branch_index);
            if (from_node >= circuit_model_->nodeCount() || to_node >= circuit_model_->nodeCount() ||
                from_node == to_node)
            {
                continue;
            }
            double branch_capacitance_weight_f = 0.0;
            if (graph_capacitance_matrix_provider_)
            {
                branch_capacitance_weight_f = std::abs(graph_capacitance_matrix_provider_->mutualCapacitanceF(
                    config_, circuit_model_.get(), branch_index));
            }
            if (branch_capacitance_weight_f <= 0.0)
            {
                const double area_i_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(from_node));
                const double area_j_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(to_node));
                branch_capacitance_weight_f = Coupling::computeHarmonicPairCapacitanceF(
                    area_i_m2, area_j_m2, shared_effective_sheath_length_m,
                    edge_operator_capacitance_scale, 1.0e-18);
            }
            const double branch_weight_f = branch_capacitance_weight_f;

            if (branch_weight_f <= 0.0)
            {
                continue;
            }

            branch_graph_adjacency_weight_f[from_node][to_node] += branch_weight_f;
            branch_graph_adjacency_weight_f[to_node][from_node] += branch_weight_f;
            ++branch_graph_edge_count;
        }
    }
    for (const auto node_index : patch_nodes)
    {
        double graph_diagonal_f = 0.0;
        double graph_row_sum_f = 0.0;
        if (graph_capacitance_matrix_provider_)
        {
            graph_diagonal_f = std::max(
                0.0, graph_capacitance_matrix_provider_->diagonalCapacitanceF(
                         config_, circuit_model_.get(), node_index));
            graph_row_sum_f = std::max(
                graph_diagonal_f,
                graph_capacitance_matrix_provider_->rowSumCapacitanceF(
                    config_, circuit_model_.get(), node_index));
        }
        const double graph_ratio =
            graph_diagonal_f > 0.0 ? graph_row_sum_f / std::max(1.0e-30, graph_diagonal_f) : 1.0;
        const double conductance_ratio =
            graph_diagonal_f > 0.0
                ? transport_graph_node_conductance_weight_f[node_index] /
                      std::max(1.0e-30, graph_diagonal_f)
                : 0.0;
        graph_operator_node_preconditioner[node_index] =
            1.0 + edge_graph_operator_node_preconditioner_weight *
                      std::max(0.0, graph_ratio - 1.0 + conductance_ratio);
        min_node_preconditioner =
            std::min(min_node_preconditioner, graph_operator_node_preconditioner[node_index]);
        max_node_preconditioner =
            std::max(max_node_preconditioner, graph_operator_node_preconditioner[node_index]);
        desired_node_edge_balance_c[node_index] /=
            std::max(1.0, graph_operator_node_preconditioner[node_index]);
    }
    for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
    {
        const auto node_i = circuit_model_->patchNodeIndex(patch_i);
        for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
        {
            const auto node_j = circuit_model_->patchNodeIndex(patch_j);
            const double area_i_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_i));
            const double area_j_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_j));
            const double fallback_pair_weight_f = Coupling::computeHarmonicPairCapacitanceF(
                area_i_m2, area_j_m2, shared_effective_sheath_length_m,
                edge_operator_capacitance_scale, 1.0e-18);
            const double direct_branch_weight_f =
                (node_i < branch_graph_adjacency_weight_f.size() &&
                 node_j < branch_graph_adjacency_weight_f[node_i].size())
                    ? branch_graph_adjacency_weight_f[node_i][node_j]
                    : 0.0;
            const double path_branch_weight_f = Coupling::computeIndirectPathBranchWeight(
                branch_graph_adjacency_weight_f, node_i, node_j);
            const double branch_pair_weight_f =
                direct_branch_weight_f + edge_graph_operator_path_blend * path_branch_weight_f;
            const double pair_transport_conductance_s =
                std::max(0.0, shared_transport_conductance_s_per_m2) *
                std::min(area_i_m2, area_j_m2);
            transport_conductance_matrix_s[node_i][node_j] = pair_transport_conductance_s;
            transport_conductance_matrix_s[node_j][node_i] = pair_transport_conductance_s;
            const double pair_transport_conductance_weight_f =
                edge_graph_operator_conductance_weight * std::max(0.0, dt) *
                pair_transport_conductance_s;
            const double pair_weight_f = std::max(
                1.0e-18,
                (branch_pair_weight_f > 0.0 ? branch_pair_weight_f : fallback_pair_weight_f) +
                    pair_transport_conductance_weight_f);
            graph_operator_pair_weight_f[node_i][node_j] = pair_weight_f;
            graph_operator_pair_weight_f[node_j][node_i] = pair_weight_f;
            graph_operator_weight_sum[node_i] += pair_weight_f;
            graph_operator_weight_sum[node_j] += pair_weight_f;
            transport_graph_node_conductance_weight_f[node_i] += pair_transport_conductance_weight_f;
            transport_graph_node_conductance_weight_f[node_j] += pair_transport_conductance_weight_f;
            total_conductance_weight_f += pair_transport_conductance_weight_f;
            if (branch_pair_weight_f > 0.0)
            {
                ++branch_graph_pair_count;
            }
            if (pair_weight_f > 0.0)
            {
                ++effective_pair_count;
                total_pair_weight_f += pair_weight_f;
            }
        }
    }

    std::size_t graph_operator_iterations = 0;
    bool graph_operator_converged = edge_graph_operator_iteration_limit == 0;
    double graph_operator_max_balance_residual_c = 0.0;
    if (patch_nodes.size() == 2)
    {
        const auto node_a = patch_nodes[0];
        const auto node_b = patch_nodes[1];
        const double pair_weight_f = graph_operator_pair_weight_f[node_a][node_b];
        if (pair_weight_f > 0.0)
        {
            graph_operator_lambda_c[node_a] =
                0.5 * desired_node_edge_balance_c[node_a] / pair_weight_f;
            graph_operator_lambda_c[node_b] =
                0.5 * desired_node_edge_balance_c[node_b] / pair_weight_f;
            graph_operator_iterations = 1;
            graph_operator_converged = true;
            graph_operator_max_balance_residual_c = 0.0;
        }
    }
    else if (edge_graph_operator_iteration_limit > 0 && patch_nodes.size() >= 2)
    {
        for (std::size_t iteration = 0; iteration < edge_graph_operator_iteration_limit; ++iteration)
        {
            for (const auto node_index : patch_nodes)
            {
                const double diagonal = graph_operator_weight_sum[node_index];
                if (diagonal <= 0.0)
                {
                    graph_operator_next_lambda_c[node_index] = 0.0;
                    continue;
                }
                double neighbor_sum_c = 0.0;
                for (const auto other_node_index : patch_nodes)
                {
                    if (other_node_index == node_index)
                    {
                        continue;
                    }
                    neighbor_sum_c +=
                        graph_operator_pair_weight_f[node_index][other_node_index] *
                        graph_operator_lambda_c[other_node_index];
                }
                const double candidate_lambda_c =
                    (desired_node_edge_balance_c[node_index] + neighbor_sum_c) / diagonal;
                graph_operator_next_lambda_c[node_index] =
                    (1.0 - edge_graph_operator_relaxation) *
                        graph_operator_lambda_c[node_index] +
                    edge_graph_operator_relaxation * candidate_lambda_c;
            }

            double lambda_mean_c = 0.0;
            for (const auto node_index : patch_nodes)
            {
                lambda_mean_c += graph_operator_next_lambda_c[node_index];
            }
            lambda_mean_c /= static_cast<double>(patch_nodes.size());
            for (const auto node_index : patch_nodes)
            {
                graph_operator_next_lambda_c[node_index] -= lambda_mean_c;
            }

            graph_operator_max_balance_residual_c = 0.0;
            for (const auto node_index : patch_nodes)
            {
                double laplacian_lambda_c = 0.0;
                for (const auto other_node_index : patch_nodes)
                {
                    if (other_node_index == node_index)
                    {
                        continue;
                    }
                    laplacian_lambda_c +=
                        graph_operator_pair_weight_f[node_index][other_node_index] *
                        (graph_operator_next_lambda_c[node_index] -
                         graph_operator_next_lambda_c[other_node_index]);
                }
                graph_operator_max_balance_residual_c = std::max(
                    graph_operator_max_balance_residual_c,
                    std::abs(laplacian_lambda_c - desired_node_edge_balance_c[node_index]));
            }

            graph_operator_lambda_c = graph_operator_next_lambda_c;
            graph_operator_iterations = iteration + 1;
            if (graph_operator_max_balance_residual_c <= edge_graph_operator_tolerance_c)
            {
                graph_operator_converged = true;
                break;
            }
        }
    }
    shared_particle_transport_edge_graph_operator_iterations_last_ =
        static_cast<double>(graph_operator_iterations);
    shared_particle_transport_edge_graph_operator_converged_last_ =
        graph_operator_converged ? 1.0 : 0.0;
    shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_ =
        graph_operator_max_balance_residual_c;
    shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_ =
        static_cast<double>(branch_graph_edge_count);
    shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_ =
        static_cast<double>(branch_graph_pair_count);
    shared_particle_transport_edge_graph_operator_effective_pair_count_last_ =
        static_cast<double>(effective_pair_count);
    shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_ =
        total_pair_weight_f;
    shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_ =
        total_conductance_weight_f;
    shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_ =
        min_node_preconditioner;
    shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_ =
        max_node_preconditioner;

    if (dt > 0.0)
    {
        for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
        {
            const auto node_i = circuit_model_->patchNodeIndex(patch_i);
            for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
            {
                const auto node_j = circuit_model_->patchNodeIndex(patch_j);
                const double pair_charge_delta_c =
                    target_node_charge_c[node_i] - target_node_charge_c[node_j];
                const double exchange_flux_a =
                    exchange_rate_hz > 0.0 ? exchange_rate_hz * pair_charge_delta_c : 0.0;
                const double area_i_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_i));
                const double area_j_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_j));
                const double node_i_potential_v =
                    node_i < node_potentials_v.size() ? node_potentials_v[node_i]
                                                      : shared_patch_potential_v;
                const double node_j_potential_v =
                    node_j < node_potentials_v.size() ? node_potentials_v[node_j]
                                                      : shared_patch_potential_v;
                const double potential_delta_v = node_j_potential_v - node_i_potential_v;
                const double edge_effective_capacitance_f =
                    Coupling::computeHarmonicPairCapacitanceF(
                        area_i_m2, area_j_m2, shared_effective_sheath_length_m,
                        edge_operator_capacitance_scale, 1.0e-18);
                const double potential_target_edge_charge_c =
                    edge_effective_capacitance_f * potential_delta_v;
                const double graph_target_edge_charge_c =
                    graph_operator_pair_weight_f[node_i][node_j] *
                    (graph_operator_lambda_c[node_i] - graph_operator_lambda_c[node_j]);
                const double target_edge_charge_c = std::clamp(
                    (1.0 - edge_graph_operator_potential_blend) * graph_target_edge_charge_c +
                        edge_graph_operator_potential_blend * potential_target_edge_charge_c,
                    -max_abs_edge_charge_c, max_abs_edge_charge_c);
                const double operator_drive_charge_c =
                    edge_operator_weight * edge_operator_relaxation *
                    (target_edge_charge_c -
                     transport_edge_charge_matrix_c[node_i][node_j]);
                transport_edge_target_charge_matrix_c[node_i][node_j] = target_edge_charge_c;
                transport_edge_target_charge_matrix_c[node_j][node_i] = -target_edge_charge_c;
                transport_edge_operator_drive_matrix_c[node_i][node_j] =
                    operator_drive_charge_c;
                transport_edge_operator_drive_matrix_c[node_j][node_i] =
                    -operator_drive_charge_c;
                if ((!std::isfinite(exchange_flux_a) || std::abs(exchange_flux_a) <= 0.0) &&
                    (!std::isfinite(operator_drive_charge_c) ||
                     std::abs(operator_drive_charge_c) <= 0.0))
                {
                    continue;
                }

                const double exchange_charge_c = exchange_flux_a * dt;
                double stored_edge_charge_c =
                    edge_memory_damping_factor *
                        transport_edge_charge_matrix_c[node_i][node_j] +
                    exchange_charge_c + operator_drive_charge_c;
                if (max_abs_edge_charge_c > 0.0)
                {
                    stored_edge_charge_c = std::clamp(stored_edge_charge_c,
                                                      -max_abs_edge_charge_c,
                                                      max_abs_edge_charge_c);
                }
                node_exchange_charge_delta_c[node_i] -= exchange_charge_c;
                node_exchange_charge_delta_c[node_j] += exchange_charge_c;
                transport_node_net_flux_a[node_i] -= exchange_flux_a;
                transport_node_net_flux_a[node_j] += exchange_flux_a;
                transport_edge_charge_matrix_c[node_i][node_j] = stored_edge_charge_c;
                transport_edge_charge_matrix_c[node_j][node_i] = -stored_edge_charge_c;
                transport_exchange_flux_matrix_a[node_i][node_j] = exchange_flux_a;
                transport_exchange_flux_matrix_a[node_j][node_i] = -exchange_flux_a;
            }
        }
    }

    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        if (node_index >= transport_edge_charge_matrix_c.size())
        {
            continue;
        }
        double net_edge_charge_c = 0.0;
        for (const double stored_edge_charge_c : transport_edge_charge_matrix_c[node_index])
        {
            net_edge_charge_c += stored_edge_charge_c;
        }
        node_edge_domain_charge_c[node_index] = net_edge_charge_c;
        node_edge_feedback_charge_c[node_index] = edge_feedback_weight * net_edge_charge_c;
        target_node_charge_c[node_index] += node_edge_feedback_charge_c[node_index];
    }

    double sum_patch_charge_c = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        transport_node_charge_c[node_index] =
            (1.0 - relaxation) * transport_node_charge_c[node_index] +
            relaxation * target_node_charge_c[node_index];
        transport_node_charge_c[node_index] += node_exchange_charge_delta_c[node_index];
        sum_patch_charge_c += transport_node_charge_c[node_index];
    }

    const double conservation_error_c = shared_particle_transport_charge_c_ - sum_patch_charge_c;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        const double area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        const double correction_share = area_m2 / total_area_m2;
        const double correction_charge_c = conservation_error_c * correction_share;
        node_conservation_correction_charge_c[node_index] = correction_charge_c;
        transport_node_charge_c[node_index] += correction_charge_c;
    }

    (void)shared_effective_sheath_length_m;
    populateOwnedGlobalParticleDomainStateFromTransportBuffers(
        node_potentials_v, transport_node_charge_c, transport_node_net_flux_a,
        transport_edge_charge_matrix_c, transport_edge_target_charge_matrix_c,
        transport_edge_operator_drive_matrix_c, transport_exchange_flux_matrix_a,
        transport_conductance_matrix_s);
    populateOwnedGlobalParticleRepositoryStateFromTransportBuffers(
        node_potentials_v, transport_node_charge_c, transport_node_net_flux_a,
        target_node_charge_c, node_exchange_charge_delta_c, node_edge_feedback_charge_c,
        node_conservation_correction_charge_c, transport_edge_charge_matrix_c,
        transport_edge_target_charge_matrix_c, transport_edge_operator_drive_matrix_c,
        transport_exchange_flux_matrix_a, transport_conductance_matrix_s, dt);
    syncSharedTransportCachesFromOwnedGlobalParticleDomainState();
}

void DensePlasmaSurfaceCharging::syncSharedTransportCachesFromOwnedGlobalParticleDomainState()
{
    if (!circuit_model_ || !global_particle_domain_state_.active)
    {
        return;
    }

    const std::size_t node_count = circuit_model_->nodeCount();
    if (shared_particle_transport_node_charge_c_.size() != node_count)
    {
        shared_particle_transport_node_charge_c_.assign(node_count, 0.0);
    }
    else
    {
        std::fill(shared_particle_transport_node_charge_c_.begin(),
                  shared_particle_transport_node_charge_c_.end(), 0.0);
    }
    if (shared_particle_transport_node_net_flux_a_.size() != node_count)
    {
        shared_particle_transport_node_net_flux_a_.assign(node_count, 0.0);
    }
    else
    {
        std::fill(shared_particle_transport_node_net_flux_a_.begin(),
                  shared_particle_transport_node_net_flux_a_.end(), 0.0);
    }

    const auto reset_square_matrix = [node_count](std::vector<std::vector<double>>& matrix) {
        if (matrix.size() != node_count)
        {
            matrix.assign(node_count, std::vector<double>(node_count, 0.0));
            return;
        }
        for (auto& row : matrix)
        {
            if (row.size() != node_count)
            {
                row.assign(node_count, 0.0);
            }
            else
            {
                std::fill(row.begin(), row.end(), 0.0);
            }
        }
    };
    reset_square_matrix(shared_particle_transport_edge_charge_matrix_c_);
    reset_square_matrix(shared_particle_transport_edge_target_charge_matrix_c_);
    reset_square_matrix(shared_particle_transport_edge_operator_drive_matrix_c_);
    reset_square_matrix(shared_particle_transport_exchange_flux_matrix_a_);
    reset_square_matrix(global_particle_transport_conductance_matrix_s_);

    shared_particle_transport_charge_c_ = global_particle_domain_state_.total_charge_c;
    shared_particle_transport_reference_shift_v_ =
        global_particle_domain_state_.total_reference_shift_v;
    shared_particle_transport_edge_graph_operator_iterations_last_ =
        global_particle_domain_state_.edge_graph_operator_iterations;
    shared_particle_transport_edge_graph_operator_converged_last_ =
        global_particle_domain_state_.edge_graph_operator_converged ? 1.0 : 0.0;
    shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_ =
        global_particle_domain_state_.edge_graph_operator_max_balance_residual_c;
    shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_ =
        global_particle_domain_state_.edge_graph_operator_branch_graph_edge_count;
    shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_ =
        global_particle_domain_state_.edge_graph_operator_branch_graph_pair_count;
    shared_particle_transport_edge_graph_operator_effective_pair_count_last_ =
        global_particle_domain_state_.edge_graph_operator_effective_pair_count;
    shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_ =
        global_particle_domain_state_.edge_graph_operator_total_pair_weight_f;
    shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_ =
        global_particle_domain_state_.edge_graph_operator_total_conductance_weight_f;
    shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_ =
        global_particle_domain_state_.edge_graph_operator_min_node_preconditioner;
    shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_ =
        global_particle_domain_state_.edge_graph_operator_max_node_preconditioner;

    for (const auto& node_state : global_particle_domain_state_.nodes)
    {
        if (node_state.node_index >= node_count)
        {
            continue;
        }
        shared_particle_transport_node_charge_c_[node_state.node_index] = node_state.charge_c;
        shared_particle_transport_node_net_flux_a_[node_state.node_index] = node_state.net_flux_a;
    }

    for (const auto& edge_state : global_particle_domain_state_.edges)
    {
        const auto from_index = edge_state.from_node_index;
        const auto to_index = edge_state.to_node_index;
        if (from_index >= node_count || to_index >= node_count || from_index == to_index)
        {
            continue;
        }
        shared_particle_transport_edge_charge_matrix_c_[from_index][to_index] =
            edge_state.stored_charge_c;
        shared_particle_transport_edge_charge_matrix_c_[to_index][from_index] =
            -edge_state.stored_charge_c;
        shared_particle_transport_edge_target_charge_matrix_c_[from_index][to_index] =
            edge_state.target_charge_c;
        shared_particle_transport_edge_target_charge_matrix_c_[to_index][from_index] =
            -edge_state.target_charge_c;
        shared_particle_transport_edge_operator_drive_matrix_c_[from_index][to_index] =
            edge_state.operator_drive_charge_c;
        shared_particle_transport_edge_operator_drive_matrix_c_[to_index][from_index] =
            -edge_state.operator_drive_charge_c;
        shared_particle_transport_exchange_flux_matrix_a_[from_index][to_index] =
            edge_state.exchange_flux_a;
        shared_particle_transport_exchange_flux_matrix_a_[to_index][from_index] =
            -edge_state.exchange_flux_a;
        global_particle_transport_conductance_matrix_s_[from_index][to_index] =
            edge_state.conductance_s;
        global_particle_transport_conductance_matrix_s_[to_index][from_index] =
            edge_state.conductance_s;
    }
}

void DensePlasmaSurfaceCharging::populateOwnedGlobalParticleDomainStateFromTransportBuffers(
    const std::vector<double>& node_potentials_v,
    const std::vector<double>& transport_node_charge_c,
    const std::vector<double>& transport_node_net_flux_a,
    const std::vector<std::vector<double>>& transport_edge_charge_matrix_c,
    const std::vector<std::vector<double>>& transport_edge_target_charge_matrix_c,
    const std::vector<std::vector<double>>& transport_edge_operator_drive_matrix_c,
    const std::vector<std::vector<double>>& transport_exchange_flux_matrix_a,
    const std::vector<std::vector<double>>& transport_conductance_matrix_s)
{
    global_particle_domain_state_ = GlobalParticleDomainState{};
    global_particle_domain_state_.bookkeeping_mode = "owned_global_particle_domain_state_v2";
    global_particle_domain_state_.total_charge_c = shared_particle_transport_charge_c_;
    global_particle_domain_state_.total_reference_shift_v = shared_particle_transport_reference_shift_v_;
    global_particle_domain_state_.edge_graph_operator_iterations =
        shared_particle_transport_edge_graph_operator_iterations_last_;
    global_particle_domain_state_.edge_graph_operator_converged =
        shared_particle_transport_edge_graph_operator_converged_last_ >= 0.5;
    global_particle_domain_state_.edge_graph_operator_max_balance_residual_c =
        shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_;
    global_particle_domain_state_.edge_graph_operator_branch_graph_edge_count =
        shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_;
    global_particle_domain_state_.edge_graph_operator_branch_graph_pair_count =
        shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_;
    global_particle_domain_state_.edge_graph_operator_effective_pair_count =
        shared_particle_transport_edge_graph_operator_effective_pair_count_last_;
    global_particle_domain_state_.edge_graph_operator_total_pair_weight_f =
        shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_;
    global_particle_domain_state_.edge_graph_operator_total_conductance_weight_f =
        shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_;
    global_particle_domain_state_.edge_graph_operator_min_node_preconditioner =
        shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_;
    global_particle_domain_state_.edge_graph_operator_max_node_preconditioner =
        shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_;

    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        return;
    }

    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        GlobalParticleDomainNodeState node_state;
        node_state.node_index = node_index;
        node_state.node_name = circuit_model_->nodeName(node_index);
        node_state.area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        node_state.patch_potential_v =
            node_index < node_potentials_v.size() ? node_potentials_v[node_index]
                                                  : circuit_model_->nodePotential(node_index);
        if (node_index < transport_node_charge_c.size())
        {
            node_state.charge_c = transport_node_charge_c[node_index];
        }
        if (node_index < transport_node_net_flux_a.size())
        {
            node_state.net_flux_a = transport_node_net_flux_a[node_index];
        }
        global_particle_domain_state_.nodes.push_back(node_state);
        global_particle_domain_state_.total_node_charge_c += node_state.charge_c;
        global_particle_domain_state_.total_node_flux_a += node_state.net_flux_a;
    }

    for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
    {
        const auto from_index = circuit_model_->patchNodeIndex(patch_i);
        for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
        {
            const auto to_index = circuit_model_->patchNodeIndex(patch_j);
            const double conductance_s =
                (from_index < transport_conductance_matrix_s.size() &&
                 to_index < transport_conductance_matrix_s[from_index].size())
                    ? transport_conductance_matrix_s[from_index][to_index]
                    : 0.0;
            const double exchange_flux_a =
                (from_index < transport_exchange_flux_matrix_a.size() &&
                 to_index < transport_exchange_flux_matrix_a[from_index].size())
                    ? transport_exchange_flux_matrix_a[from_index][to_index]
                    : 0.0;
            const double stored_charge_c =
                (from_index < transport_edge_charge_matrix_c.size() &&
                 to_index < transport_edge_charge_matrix_c[from_index].size())
                    ? transport_edge_charge_matrix_c[from_index][to_index]
                    : 0.0;
            const double target_charge_c =
                (from_index < transport_edge_target_charge_matrix_c.size() &&
                 to_index < transport_edge_target_charge_matrix_c[from_index].size())
                    ? transport_edge_target_charge_matrix_c[from_index][to_index]
                    : 0.0;
            const double operator_drive_charge_c =
                (from_index < transport_edge_operator_drive_matrix_c.size() &&
                 to_index < transport_edge_operator_drive_matrix_c[from_index].size())
                    ? transport_edge_operator_drive_matrix_c[from_index][to_index]
                    : 0.0;
            if (std::abs(exchange_flux_a) <= 0.0 && std::abs(stored_charge_c) <= 0.0 &&
                std::abs(target_charge_c) <= 0.0 &&
                std::abs(operator_drive_charge_c) <= 0.0)
            {
                continue;
            }

            GlobalParticleDomainEdgeState edge_state;
            edge_state.from_node_index = from_index;
            edge_state.to_node_index = to_index;
            edge_state.from_node_name = circuit_model_->nodeName(from_index);
            edge_state.to_node_name = circuit_model_->nodeName(to_index);
            edge_state.conductance_s = conductance_s;
            edge_state.exchange_flux_a = exchange_flux_a;
            edge_state.stored_charge_c = stored_charge_c;
            edge_state.target_charge_c = target_charge_c;
            edge_state.operator_drive_charge_c = operator_drive_charge_c;
            global_particle_domain_state_.edges.push_back(edge_state);
            global_particle_domain_state_.edge_charge_total_abs_c += std::abs(stored_charge_c);
            global_particle_domain_state_.edge_target_charge_total_abs_c +=
                std::abs(target_charge_c);
            global_particle_domain_state_.edge_operator_drive_total_abs_c +=
                std::abs(operator_drive_charge_c);
            global_particle_domain_state_.edge_conductance_total_s +=
                std::max(0.0, conductance_s);
        }
    }

    global_particle_domain_state_.charge_conservation_error_c =
        std::abs(global_particle_domain_state_.total_node_charge_c -
                 global_particle_domain_state_.total_charge_c);
    global_particle_domain_state_.flux_conservation_error_a =
        std::abs(global_particle_domain_state_.total_node_flux_a);
    global_particle_domain_state_.active =
        global_particle_domain_state_.nodes.size() >= 2 &&
        (std::abs(global_particle_domain_state_.total_charge_c) > 0.0 ||
         !global_particle_domain_state_.edges.empty() ||
         global_particle_domain_state_.edge_charge_total_abs_c > 0.0);

    for (auto& node_state : global_particle_domain_state_.nodes)
    {
        node_state.reference_shift_v = computeSharedSurfaceDistributedParticleTransportReferenceShiftV(
            node_state.node_index, node_state.area_m2,
            computeSharedSurfaceEffectiveSheathLengthM(computeEffectiveSheathLength()));
        node_state.shared_reference_potential_v =
            computeSharedSurfaceReferencePotentialV(config_.live_pic_reference_potential_v,
                                                   node_state.patch_potential_v) +
            node_state.reference_shift_v;
    }
}

void DensePlasmaSurfaceCharging::populateOwnedGlobalParticleRepositoryStateFromTransportBuffers(
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
    double dt)
{
    global_particle_repository_state_ = GlobalParticleRepositoryState{};
    global_particle_repository_state_.bookkeeping_mode =
        "owned_global_particle_repository_state_v1";
    global_particle_repository_state_.lifecycle_mode =
        "global_particle_transport_reservoir_lifecycle_v1";

    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        return;
    }

    double migration_delta_sum_c = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        GlobalParticleRepositoryNodeState node_state;
        node_state.node_index = node_index;
        node_state.node_name = circuit_model_->nodeName(node_index);
        node_state.area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        node_state.patch_potential_v =
            node_index < node_potentials_v.size() ? node_potentials_v[node_index]
                                                  : circuit_model_->nodePotential(node_index);
        if (node_index < transport_node_charge_c.size())
        {
            node_state.reservoir_charge_c = transport_node_charge_c[node_index];
        }
        if (node_index < target_node_charge_c.size())
        {
            node_state.target_reservoir_charge_c = target_node_charge_c[node_index];
        }
        if (node_index < node_migration_delta_charge_c.size())
        {
            node_state.migration_delta_charge_c = node_migration_delta_charge_c[node_index];
        }
        if (node_index < node_edge_feedback_charge_c.size())
        {
            node_state.edge_feedback_charge_c = node_edge_feedback_charge_c[node_index];
        }
        if (node_index < node_conservation_correction_charge_c.size())
        {
            node_state.conservation_correction_charge_c =
                node_conservation_correction_charge_c[node_index];
        }
        if (node_index < transport_node_net_flux_a.size())
        {
            node_state.net_flux_a = transport_node_net_flux_a[node_index];
        }
        if (const auto* domain_node = findGlobalParticleDomainNodeState(node_index))
        {
            node_state.reference_shift_v = domain_node->reference_shift_v;
            node_state.shared_reference_potential_v = domain_node->shared_reference_potential_v;
        }
        else
        {
            node_state.reference_shift_v = computeSharedSurfaceDistributedParticleTransportReferenceShiftV(
                node_state.node_index, node_state.area_m2,
                computeSharedSurfaceEffectiveSheathLengthM(computeEffectiveSheathLength()));
            node_state.shared_reference_potential_v =
                computeSharedSurfaceReferencePotentialV(config_.live_pic_reference_potential_v,
                                                       node_state.patch_potential_v) +
                node_state.reference_shift_v;
        }

        global_particle_repository_state_.total_reservoir_charge_c +=
            node_state.reservoir_charge_c;
        global_particle_repository_state_.total_target_reservoir_charge_c +=
            node_state.target_reservoir_charge_c;
        global_particle_repository_state_.total_migration_delta_abs_charge_c +=
            std::abs(node_state.migration_delta_charge_c);
        global_particle_repository_state_.total_edge_feedback_abs_charge_c +=
            std::abs(node_state.edge_feedback_charge_c);
        global_particle_repository_state_.total_conservation_correction_abs_charge_c +=
            std::abs(node_state.conservation_correction_charge_c);
        migration_delta_sum_c += node_state.migration_delta_charge_c;
        global_particle_repository_state_.nodes.push_back(node_state);
    }

    for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
    {
        const auto from_index = circuit_model_->patchNodeIndex(patch_i);
        for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
        {
            const auto to_index = circuit_model_->patchNodeIndex(patch_j);
            const double conductance_s =
                (from_index < transport_conductance_matrix_s.size() &&
                 to_index < transport_conductance_matrix_s[from_index].size())
                    ? transport_conductance_matrix_s[from_index][to_index]
                    : 0.0;
            const double migration_flux_a =
                (from_index < transport_exchange_flux_matrix_a.size() &&
                 to_index < transport_exchange_flux_matrix_a[from_index].size())
                    ? transport_exchange_flux_matrix_a[from_index][to_index]
                    : 0.0;
            const double stored_charge_c =
                (from_index < transport_edge_charge_matrix_c.size() &&
                 to_index < transport_edge_charge_matrix_c[from_index].size())
                    ? transport_edge_charge_matrix_c[from_index][to_index]
                    : 0.0;
            const double target_charge_c =
                (from_index < transport_edge_target_charge_matrix_c.size() &&
                 to_index < transport_edge_target_charge_matrix_c[from_index].size())
                    ? transport_edge_target_charge_matrix_c[from_index][to_index]
                    : 0.0;
            const double operator_drive_charge_c =
                (from_index < transport_edge_operator_drive_matrix_c.size() &&
                 to_index < transport_edge_operator_drive_matrix_c[from_index].size())
                    ? transport_edge_operator_drive_matrix_c[from_index][to_index]
                    : 0.0;
            const double migration_charge_c = migration_flux_a * dt;
            if (std::abs(conductance_s) <= 0.0 && std::abs(migration_flux_a) <= 0.0 &&
                std::abs(stored_charge_c) <= 0.0 && std::abs(target_charge_c) <= 0.0 &&
                std::abs(operator_drive_charge_c) <= 0.0)
            {
                continue;
            }

            GlobalParticleRepositoryEdgeState edge_state;
            edge_state.from_node_index = from_index;
            edge_state.to_node_index = to_index;
            edge_state.from_node_name = circuit_model_->nodeName(from_index);
            edge_state.to_node_name = circuit_model_->nodeName(to_index);
            edge_state.conductance_s = conductance_s;
            edge_state.migration_flux_a = migration_flux_a;
            edge_state.migration_charge_c = migration_charge_c;
            edge_state.stored_charge_c = stored_charge_c;
            edge_state.target_charge_c = target_charge_c;
            edge_state.operator_drive_charge_c = operator_drive_charge_c;
            global_particle_repository_state_.edges.push_back(edge_state);
            global_particle_repository_state_.total_migration_edge_abs_charge_c +=
                std::abs(migration_charge_c);
        }
    }

    global_particle_repository_state_.charge_conservation_error_c =
        std::abs(global_particle_repository_state_.total_reservoir_charge_c -
                 shared_particle_transport_charge_c_);
    global_particle_repository_state_.migration_charge_conservation_error_c =
        std::abs(migration_delta_sum_c);
    global_particle_repository_state_.active =
        global_particle_repository_state_.nodes.size() >= 2 &&
        (std::abs(global_particle_repository_state_.total_reservoir_charge_c) > 0.0 ||
         !global_particle_repository_state_.edges.empty() ||
         global_particle_repository_state_.total_migration_edge_abs_charge_c > 0.0);
}

void DensePlasmaSurfaceCharging::updateOwnedGlobalParticleDomainState(
    const std::vector<double>& node_potentials_v)
{
    populateOwnedGlobalParticleDomainStateFromTransportBuffers(
        node_potentials_v, shared_particle_transport_node_charge_c_,
        shared_particle_transport_node_net_flux_a_, shared_particle_transport_edge_charge_matrix_c_,
        shared_particle_transport_edge_target_charge_matrix_c_,
        shared_particle_transport_edge_operator_drive_matrix_c_,
        shared_particle_transport_exchange_flux_matrix_a_,
        global_particle_transport_conductance_matrix_s_);
    populateOwnedGlobalParticleRepositoryStateFromTransportBuffers(
        node_potentials_v, shared_particle_transport_node_charge_c_,
        shared_particle_transport_node_net_flux_a_, shared_particle_transport_node_charge_c_,
        std::vector<double>(shared_particle_transport_node_charge_c_.size(), 0.0),
        std::vector<double>(shared_particle_transport_node_charge_c_.size(), 0.0),
        std::vector<double>(shared_particle_transport_node_charge_c_.size(), 0.0),
        shared_particle_transport_edge_charge_matrix_c_,
        shared_particle_transport_edge_target_charge_matrix_c_,
        shared_particle_transport_edge_operator_drive_matrix_c_,
        shared_particle_transport_exchange_flux_matrix_a_,
        global_particle_transport_conductance_matrix_s_, 0.0);
}

void DensePlasmaSurfaceCharging::rebuildOwnedGlobalParticleDomainEdges()
{
    global_particle_domain_state_.edges.clear();
    global_particle_domain_state_.edge_charge_total_abs_c = 0.0;
    global_particle_domain_state_.edge_target_charge_total_abs_c = 0.0;
    global_particle_domain_state_.edge_operator_drive_total_abs_c = 0.0;
    global_particle_domain_state_.edge_conductance_total_s = 0.0;

    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        return;
    }

    for (std::size_t patch_i = 0; patch_i < circuit_model_->patchCount(); ++patch_i)
    {
        const auto from_index = circuit_model_->patchNodeIndex(patch_i);
        for (std::size_t patch_j = patch_i + 1; patch_j < circuit_model_->patchCount(); ++patch_j)
        {
            const auto to_index = circuit_model_->patchNodeIndex(patch_j);

            const double conductance_s =
                (from_index < global_particle_transport_conductance_matrix_s_.size() &&
                 to_index <
                     global_particle_transport_conductance_matrix_s_[from_index].size())
                    ? global_particle_transport_conductance_matrix_s_[from_index][to_index]
                    : 0.0;
            const double exchange_flux_a =
                (from_index < shared_particle_transport_exchange_flux_matrix_a_.size() &&
                 to_index <
                     shared_particle_transport_exchange_flux_matrix_a_[from_index].size())
                    ? shared_particle_transport_exchange_flux_matrix_a_[from_index][to_index]
                    : 0.0;
            const double stored_charge_c =
                (from_index < shared_particle_transport_edge_charge_matrix_c_.size() &&
                 to_index <
                     shared_particle_transport_edge_charge_matrix_c_[from_index].size())
                    ? shared_particle_transport_edge_charge_matrix_c_[from_index][to_index]
                    : 0.0;
            const double target_charge_c =
                (from_index < shared_particle_transport_edge_target_charge_matrix_c_.size() &&
                 to_index <
                     shared_particle_transport_edge_target_charge_matrix_c_[from_index].size())
                    ? shared_particle_transport_edge_target_charge_matrix_c_[from_index][to_index]
                    : 0.0;
            const double operator_drive_charge_c =
                (from_index < shared_particle_transport_edge_operator_drive_matrix_c_.size() &&
                 to_index <
                     shared_particle_transport_edge_operator_drive_matrix_c_[from_index].size())
                    ? shared_particle_transport_edge_operator_drive_matrix_c_[from_index][to_index]
                    : 0.0;
            if (std::abs(exchange_flux_a) <= 0.0 && std::abs(stored_charge_c) <= 0.0 &&
                std::abs(target_charge_c) <= 0.0 &&
                std::abs(operator_drive_charge_c) <= 0.0)
            {
                continue;
            }

            GlobalParticleDomainEdgeState edge_state;
            edge_state.from_node_index = from_index;
            edge_state.to_node_index = to_index;
            edge_state.from_node_name = circuit_model_->nodeName(from_index);
            edge_state.to_node_name = circuit_model_->nodeName(to_index);
            edge_state.conductance_s = conductance_s;
            edge_state.exchange_flux_a = exchange_flux_a;
            edge_state.stored_charge_c = stored_charge_c;
            edge_state.target_charge_c = target_charge_c;
            edge_state.operator_drive_charge_c = operator_drive_charge_c;
            global_particle_domain_state_.edges.push_back(edge_state);
            global_particle_domain_state_.edge_charge_total_abs_c += std::abs(stored_charge_c);
            global_particle_domain_state_.edge_target_charge_total_abs_c +=
                std::abs(target_charge_c);
            global_particle_domain_state_.edge_operator_drive_total_abs_c +=
                std::abs(operator_drive_charge_c);
            global_particle_domain_state_.edge_conductance_total_s +=
                std::max(0.0, conductance_s);
        }
    }
}

void DensePlasmaSurfaceCharging::solveOwnedGlobalSheathFieldSystem(
    const std::vector<double>& node_potentials_v, double shared_effective_sheath_length_m, double dt)
{
    global_sheath_field_solve_state_ = GlobalSheathFieldSolveState{};
    global_sheath_field_solve_state_.solve_mode = "explicit_global_sheath_field_linear_system_v2";
    global_sheath_field_solve_state_.effective_sheath_length_m =
        std::max(1.0e-6, shared_effective_sheath_length_m);

    if (!circuit_model_ || !useSharedSurfaceGlobalCoupledSolve() || circuit_model_->patchCount() < 2)
    {
        if (circuit_model_)
        {
            global_sheath_field_reference_response_matrix_.assign(
                circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));
        }
        return;
    }

    const std::size_t patch_count = circuit_model_->patchCount();
    const double sheath_length_m = global_sheath_field_solve_state_.effective_sheath_length_m;
    const double anchor_reference_v =
        std::isfinite(config_.live_pic_reference_potential_v) ? config_.live_pic_reference_potential_v
                                                              : 0.0;
    const double anchor_weight_scale = std::clamp(
        config_.material.getScalarProperty("global_sheath_field_anchor_weight_scale", 0.15), 0.0,
        4.0);
    const double particle_charge_relaxation = std::clamp(
        config_.material.getScalarProperty("global_sheath_field_particle_charge_relaxation", 1.0),
        0.0, 1.0);
    const double flux_source_relaxation_time_s = std::max(
        1.0e-9, config_.material.getScalarProperty(
                    "global_sheath_field_flux_relaxation_time_s", std::max(dt, 4.0 * dt)));

    std::vector<std::size_t> patch_nodes(patch_count, 0);
    std::vector<std::vector<double>> matrix(patch_count, std::vector<double>(patch_count, 0.0));
    std::vector<double> rhs(patch_count, 0.0);
    std::vector<double> sheath_capacitance_f(patch_count, 0.0);
    std::vector<double> uncoupled_reference_v(patch_count, 0.0);
    std::vector<double> patch_potential_v(patch_count, 0.0);
    std::vector<double> particle_charge_c(patch_count, 0.0);
    std::vector<double> particle_flux_a(patch_count, 0.0);
    std::vector<double> anchor_weight_f(patch_count, 0.0);

    for (std::size_t patch_ordinal = 0; patch_ordinal < patch_count; ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        patch_nodes[patch_ordinal] = node_index;
        const double area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        const double local_capacitance_scale =
            graph_capacitance_matrix_provider_ != nullptr
                ? std::max(1.0, graph_capacitance_matrix_provider_->rowSumCapacitanceF(
                                       config_, circuit_model_.get(), node_index) /
                                       std::max(1.0e-18, graph_capacitance_matrix_provider_
                                                            ->diagonalCapacitanceF(
                                                                config_, circuit_model_.get(),
                                                                node_index)))
                : 1.0;
        sheath_capacitance_f[patch_ordinal] =
            std::max(1.0e-18, local_capacitance_scale * kEpsilon0 * area_m2 / sheath_length_m);
        patch_potential_v[patch_ordinal] =
            node_index < node_potentials_v.size() ? node_potentials_v[node_index]
                                                  : circuit_model_->nodePotential(node_index);

        if (const auto* particle_node = findGlobalParticleDomainNodeState(node_index))
        {
            particle_charge_c[patch_ordinal] =
                particle_charge_relaxation * particle_node->charge_c;
            particle_flux_a[patch_ordinal] = particle_node->net_flux_a;
        }
        uncoupled_reference_v[patch_ordinal] =
            patch_potential_v[patch_ordinal] +
            particle_charge_c[patch_ordinal] / std::max(1.0e-18, sheath_capacitance_f[patch_ordinal]);
        anchor_weight_f[patch_ordinal] = anchor_weight_scale * sheath_capacitance_f[patch_ordinal];

        matrix[patch_ordinal][patch_ordinal] +=
            sheath_capacitance_f[patch_ordinal] + anchor_weight_f[patch_ordinal];
        rhs[patch_ordinal] +=
            sheath_capacitance_f[patch_ordinal] * uncoupled_reference_v[patch_ordinal] +
            anchor_weight_f[patch_ordinal] * anchor_reference_v +
            particle_flux_a[patch_ordinal] * flux_source_relaxation_time_s;
    }

    for (std::size_t patch_i = 0; patch_i < patch_count; ++patch_i)
    {
        const auto node_i = patch_nodes[patch_i];
        for (std::size_t patch_j = patch_i + 1; patch_j < patch_count; ++patch_j)
        {
            const auto node_j = patch_nodes[patch_j];
            double mutual_weight_f = 0.0;
            if (global_particle_transport_conductance_matrix_s_.size() > node_i &&
                global_particle_transport_conductance_matrix_s_[node_i].size() > node_j)
            {
                mutual_weight_f += std::max(
                    0.0, std::max(dt, 0.0) *
                             global_particle_transport_conductance_matrix_s_[node_i][node_j]);
            }
            if (graph_capacitance_matrix_provider_ != nullptr)
            {
                for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount();
                     ++branch_index)
                {
                    const auto from_node = circuit_model_->branchFromNodeIndex(branch_index);
                    const auto to_node = circuit_model_->branchToNodeIndex(branch_index);
                    const bool matches_pair =
                        (from_node == node_i && to_node == node_j) ||
                        (from_node == node_j && to_node == node_i);
                    if (!matches_pair)
                    {
                        continue;
                    }
                    mutual_weight_f += std::max(
                        0.0, graph_capacitance_matrix_provider_->mutualCapacitanceF(
                                 config_, circuit_model_.get(), branch_index));
                }
            }
            if (mutual_weight_f <= 0.0)
            {
                const double area_i_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_i));
                const double area_j_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_j));
                mutual_weight_f = Coupling::computeHarmonicPairCapacitanceF(
                    area_i_m2, area_j_m2, sheath_length_m, 1.0, 1.0e-18);
            }
            matrix[patch_i][patch_i] += mutual_weight_f;
            matrix[patch_j][patch_j] += mutual_weight_f;
            matrix[patch_i][patch_j] -= mutual_weight_f;
            matrix[patch_j][patch_i] -= mutual_weight_f;
        }
    }

    const auto solved = Solver::solveDenseLinearSystemWithResidual(matrix, rhs);
    global_sheath_field_solve_state_.matrix_row_count = static_cast<double>(patch_count);
    global_sheath_field_solve_state_.matrix_nonzeros = static_cast<double>(solved.nonzeros);
    global_sheath_field_solve_state_.linear_residual_norm_v = solved.residual_norm;
    if (!solved.solved)
    {
        return;
    }

    global_sheath_field_reference_response_matrix_.assign(
        circuit_model_->nodeCount(), std::vector<double>(circuit_model_->nodeCount(), 0.0));

    double total_area_m2 = 0.0;
    double weighted_patch_potential_v = 0.0;
    double weighted_reference_v = 0.0;
    double weighted_field_v_per_m = 0.0;
    double weighted_charge_density_c_per_m3 = 0.0;
    double max_field_residual_v_per_m = 0.0;
    double max_particle_field_residual_v = 0.0;

    for (std::size_t patch_ordinal = 0; patch_ordinal < patch_count; ++patch_ordinal)
    {
        const auto node_index = patch_nodes[patch_ordinal];
        const double area_m2 = std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index));
        const double reference_potential_v = solved.solution[patch_ordinal];
        const double field_v_per_m =
            (reference_potential_v - patch_potential_v[patch_ordinal]) / sheath_length_m;
        const double charge_density_c_per_m3 =
            kEpsilon0 * field_v_per_m / sheath_length_m;
        const double sheath_charge_c =
            sheath_capacitance_f[patch_ordinal] *
            (reference_potential_v - patch_potential_v[patch_ordinal]);
        const double particle_field_residual_v =
            std::abs(sheath_charge_c - particle_charge_c[patch_ordinal]) /
            std::max(1.0e-18, sheath_capacitance_f[patch_ordinal]);

        GlobalSheathFieldNodeState node_state;
        node_state.node_index = node_index;
        node_state.node_name = circuit_model_->nodeName(node_index);
        node_state.area_m2 = area_m2;
        node_state.patch_potential_v = patch_potential_v[patch_ordinal];
        node_state.uncoupled_reference_potential_v = uncoupled_reference_v[patch_ordinal];
        node_state.reference_potential_v = reference_potential_v;
        node_state.sheath_capacitance_f = sheath_capacitance_f[patch_ordinal];
        node_state.sheath_charge_c = sheath_charge_c;
        node_state.normal_electric_field_v_per_m = field_v_per_m;
        node_state.local_charge_density_c_per_m3 = charge_density_c_per_m3;
        node_state.particle_charge_c = particle_charge_c[patch_ordinal];
        node_state.particle_flux_a = particle_flux_a[patch_ordinal];
        global_sheath_field_solve_state_.nodes.push_back(node_state);

        total_area_m2 += area_m2;
        weighted_patch_potential_v += area_m2 * patch_potential_v[patch_ordinal];
        weighted_reference_v += area_m2 * reference_potential_v;
        weighted_field_v_per_m += area_m2 * field_v_per_m;
        weighted_charge_density_c_per_m3 += area_m2 * charge_density_c_per_m3;
        std::vector<double> derivative_rhs(patch_count, 0.0);
        derivative_rhs[patch_ordinal] = sheath_capacitance_f[patch_ordinal];
        const auto response =
            Solver::solveDenseLinearSystemWithResidual(matrix, derivative_rhs);
        if (response.solved)
        {
            for (std::size_t patch_j = 0; patch_j < patch_count; ++patch_j)
            {
                const auto response_node_index = patch_nodes[patch_j];
                global_sheath_field_reference_response_matrix_[node_index][response_node_index] =
                    response.solution[patch_j];
            }
        }
    }

    if (total_area_m2 <= 0.0)
    {
        return;
    }

    global_sheath_field_solve_state_.global_patch_potential_v =
        weighted_patch_potential_v / total_area_m2;
    global_sheath_field_solve_state_.global_reference_potential_v =
        weighted_reference_v / total_area_m2;
    global_sheath_field_solve_state_.global_normal_electric_field_v_per_m =
        weighted_field_v_per_m / total_area_m2;
    global_sheath_field_solve_state_.global_local_charge_density_c_per_m3 =
        weighted_charge_density_c_per_m3 / total_area_m2;

    // Treat the formal sheath-field solve residual as the actual residual of the assembled
    // global linear system, rather than the distance from a local uncoupled reference guess.
    for (std::size_t patch_ordinal = 0; patch_ordinal < patch_count; ++patch_ordinal)
    {
        double row_residual_c = rhs[patch_ordinal];
        for (std::size_t column = 0; column < patch_count; ++column)
        {
            row_residual_c -= matrix[patch_ordinal][column] * solved.solution[column];
        }

        const double local_capacitance_f =
            std::max(1.0e-18, sheath_capacitance_f[patch_ordinal] + anchor_weight_f[patch_ordinal]);
        const double local_voltage_residual_v = std::abs(row_residual_c) / local_capacitance_f;
        max_field_residual_v_per_m =
            std::max(max_field_residual_v_per_m, local_voltage_residual_v / sheath_length_m);
        max_particle_field_residual_v =
            std::max(max_particle_field_residual_v, local_voltage_residual_v);
    }

    global_sheath_field_solve_state_.field_residual_v_per_m = max_field_residual_v_per_m;
    global_sheath_field_solve_state_.particle_field_coupled_residual_v =
        std::max(max_particle_field_residual_v,
                 global_particle_domain_state_.charge_conservation_error_c *
                     sheath_length_m / std::max(1.0e-18, total_area_m2 * kEpsilon0));
    global_sheath_field_solve_state_.multi_step_stability_metric_v =
        computeSharedSurfaceGlobalSheathFieldMultiStepStabilityMetricV();
    global_sheath_field_solve_state_.active = true;
}

void DensePlasmaSurfaceCharging::appendOwnedGlobalCoupledLinearization(
    const std::vector<double>& node_potentials_v, double dt,
    Coupling::SurfaceCircuitLinearization& linearization,
    std::size_t& current_matrix_coupling_offdiag_entries,
    std::size_t& particle_transport_offdiag_entries,
    double& particle_transport_total_conductance_s,
    double& particle_transport_conservation_error_a_per_v) const
{
    if (!circuit_model_ || dt <= 0.0)
    {
        return;
    }

    if (global_sheath_field_solve_state_.active)
    {
        for (const auto& sheath_node : global_sheath_field_solve_state_.nodes)
        {
            if (sheath_node.node_index >= linearization.additional_diagonal_a_per_v.size() ||
                sheath_node.node_index >= node_potentials_v.size())
            {
                continue;
            }

            const double current_old_a =
                (sheath_node.sheath_charge_c -
                 (sheath_node.node_index < global_sheath_field_previous_node_charge_c_.size()
                      ? global_sheath_field_previous_node_charge_c_[sheath_node.node_index]
                      : 0.0)) /
                dt;
            const std::vector<double>* response_row = nullptr;
            if (sheath_node.node_index < global_sheath_field_reference_response_matrix_.size())
            {
                response_row =
                    &global_sheath_field_reference_response_matrix_[sheath_node.node_index];
            }

            for (const auto& response_node : global_sheath_field_solve_state_.nodes)
            {
                if (response_row == nullptr || response_node.node_index >= node_potentials_v.size() ||
                    response_node.node_index >= response_row->size())
                {
                    continue;
                }
                const double jacobian_a_per_v =
                    sheath_node.sheath_capacitance_f *
                    ((*response_row)[response_node.node_index] -
                     (response_node.node_index == sheath_node.node_index ? 1.0 : 0.0)) /
                    dt;
                if (!std::isfinite(jacobian_a_per_v) || std::abs(jacobian_a_per_v) <= 0.0)
                {
                    continue;
                }

                if (response_node.node_index == sheath_node.node_index)
                {
                    linearization.additional_diagonal_a_per_v[sheath_node.node_index] +=
                        jacobian_a_per_v;
                    linearization.additional_rhs_a[sheath_node.node_index] +=
                        current_old_a - jacobian_a_per_v *
                                            node_potentials_v[sheath_node.node_index];
                }
                else
                {
                    linearization.additional_rhs_a[sheath_node.node_index] -=
                        jacobian_a_per_v * node_potentials_v[response_node.node_index];
                    Coupling::SurfaceCircuitLinearization::OffDiagonalEntry entry{};
                    entry.row_node = sheath_node.node_index;
                    entry.column_node = response_node.node_index;
                    entry.coefficient_a_per_v = -jacobian_a_per_v;
                    linearization.additional_off_diagonal_entries.push_back(entry);
                    ++current_matrix_coupling_offdiag_entries;
                }
            }
        }
    }

    if (global_particle_domain_state_.active)
    {
        particle_transport_conservation_error_a_per_v =
            std::max(global_particle_domain_state_.flux_conservation_error_a,
                     global_particle_domain_state_.charge_conservation_error_c /
                         std::max(1.0e-12, dt));
        for (const auto& edge_state : global_particle_domain_state_.edges)
        {
            const double pair_conductance_s = edge_state.conductance_s;
            if (!std::isfinite(pair_conductance_s) || pair_conductance_s <= 0.0)
            {
                continue;
            }

            linearization.additional_diagonal_a_per_v[edge_state.from_node_index] +=
                pair_conductance_s;
            linearization.additional_diagonal_a_per_v[edge_state.to_node_index] +=
                pair_conductance_s;

            Coupling::SurfaceCircuitLinearization::OffDiagonalEntry from_to{};
            from_to.row_node = edge_state.from_node_index;
            from_to.column_node = edge_state.to_node_index;
            from_to.coefficient_a_per_v = -pair_conductance_s;
            linearization.additional_off_diagonal_entries.push_back(from_to);

            Coupling::SurfaceCircuitLinearization::OffDiagonalEntry to_from{};
            to_from.row_node = edge_state.to_node_index;
            to_from.column_node = edge_state.from_node_index;
            to_from.coefficient_a_per_v = -pair_conductance_s;
            linearization.additional_off_diagonal_entries.push_back(to_from);

            particle_transport_total_conductance_s += pair_conductance_s;
            particle_transport_offdiag_entries += 2;
        }
    }
}

double DensePlasmaSurfaceCharging::sharedSurfaceDistributedParticleTransportChargeC(
    std::size_t node_index) const
{
    if (global_particle_domain_state_.active)
    {
        if (const auto* node_state = findGlobalParticleDomainNodeState(node_index))
        {
            return node_state->charge_c;
        }
    }
    if (node_index >= shared_particle_transport_node_charge_c_.size())
    {
        return 0.0;
    }
    return shared_particle_transport_node_charge_c_[node_index];
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceDistributedParticleTransportReferenceShiftV(
    std::size_t node_index, double node_area_m2, double shared_effective_sheath_length_m) const
{
    if (!useSharedSurfaceParticleTransportCoupling())
    {
        return 0.0;
    }

    const double area_m2 = std::max(1.0e-16, node_area_m2);
    double shared_patch_area_m2 = computeSharedSurfacePatchAreaM2();
    double total_transport_charge_c = shared_particle_transport_charge_c_;
    if (global_particle_domain_state_.active)
    {
        total_transport_charge_c = global_particle_domain_state_.total_charge_c;
        if (!global_particle_domain_state_.nodes.empty())
        {
            shared_patch_area_m2 = 0.0;
            for (const auto& node_state : global_particle_domain_state_.nodes)
            {
                shared_patch_area_m2 += std::max(1.0e-16, node_state.area_m2);
            }
        }
    }
    shared_patch_area_m2 = std::max(1.0e-16, shared_patch_area_m2);
    const double base_charge_c = total_transport_charge_c * area_m2 / shared_patch_area_m2;
    const double distributed_delta_charge_c =
        sharedSurfaceDistributedParticleTransportChargeC(node_index) - base_charge_c;
    double net_edge_stored_charge_c = 0.0;
    if (global_particle_domain_state_.active)
    {
        for (const auto& edge_state : global_particle_domain_state_.edges)
        {
            if (edge_state.from_node_index == node_index)
            {
                net_edge_stored_charge_c += edge_state.stored_charge_c;
            }
            else if (edge_state.to_node_index == node_index)
            {
                net_edge_stored_charge_c -= edge_state.stored_charge_c;
            }
        }
    }
    else if (node_index < shared_particle_transport_edge_charge_matrix_c_.size())
    {
        for (const double edge_charge_c : shared_particle_transport_edge_charge_matrix_c_[node_index])
        {
            net_edge_stored_charge_c += edge_charge_c;
        }
    }
    const double normalization_charge_c = std::max(1.0e-18, std::abs(total_transport_charge_c));
    const double normalized_delta =
        distributed_delta_charge_c / normalization_charge_c;
    const double normalized_edge_delta = net_edge_stored_charge_c / normalization_charge_c;
    const double distribution_shift_weight = std::clamp(
        config_.material.getScalarProperty(
            "shared_surface_particle_transport_distribution_shift_weight", 0.25),
        0.0, 1.0);
    const double edge_shift_weight = std::clamp(
        config_.material.getScalarProperty("shared_surface_particle_transport_edge_shift_weight",
                                           0.0),
        0.0, 1.0);
    const double raw_shift_v =
        currentSharedSurfaceParticleTransportReferenceShiftV() *
        (distribution_shift_weight * normalized_delta +
         edge_shift_weight * normalized_edge_delta);
    const double max_abs_shift_v = std::max(
        0.0, config_.material.getScalarProperty(
                 "shared_surface_particle_transport_distribution_max_abs_reference_shift_v", 5.0));
    (void)shared_effective_sheath_length_m;
    return std::clamp(raw_shift_v, -max_abs_shift_v, max_abs_shift_v);
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalParticleDomainChargeConservationErrorC()
    const
{
    if (global_particle_domain_state_.active)
    {
        return global_particle_domain_state_.charge_conservation_error_c;
    }
    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        return 0.0;
    }

    double total_patch_charge_c = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        total_patch_charge_c +=
            sharedSurfaceDistributedParticleTransportChargeC(circuit_model_->patchNodeIndex(patch_ordinal));
    }
    return std::abs(total_patch_charge_c - shared_particle_transport_charge_c_);
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalParticleDomainFluxConservationErrorA()
    const
{
    if (global_particle_domain_state_.active)
    {
        return global_particle_domain_state_.flux_conservation_error_a;
    }
    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        return 0.0;
    }

    double total_flux_a = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto node_index = circuit_model_->patchNodeIndex(patch_ordinal);
        if (node_index < shared_particle_transport_node_net_flux_a_.size())
        {
            total_flux_a += shared_particle_transport_node_net_flux_a_[node_index];
        }
    }
    return std::abs(total_flux_a);
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalParticleDomainEdgeChargeTotalAbsC() const
{
    if (global_particle_domain_state_.active)
    {
        return global_particle_domain_state_.edge_charge_total_abs_c;
    }
    if (!circuit_model_ || !useSharedSurfaceParticleTransportCoupling())
    {
        return 0.0;
    }

    double total_abs_edge_charge_c = 0.0;
    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
    {
        const auto from_index = circuit_model_->patchNodeIndex(patch_ordinal);
        if (from_index >= shared_particle_transport_edge_charge_matrix_c_.size())
        {
            continue;
        }
        for (std::size_t other_patch_ordinal = patch_ordinal + 1;
             other_patch_ordinal < circuit_model_->patchCount(); ++other_patch_ordinal)
        {
            const auto to_index = circuit_model_->patchNodeIndex(other_patch_ordinal);
            if (to_index >= shared_particle_transport_edge_charge_matrix_c_[from_index].size())
            {
                continue;
            }
            total_abs_edge_charge_c +=
                std::abs(shared_particle_transport_edge_charge_matrix_c_[from_index][to_index]);
        }
    }
    return total_abs_edge_charge_c;
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalSheathFieldResidualVPerM(
    const SurfaceModelRuntimeState& state) const
{
    if (!state.global_sheath_field_solve_active)
    {
        return 0.0;
    }

    const double sheath_length_m = std::max(1.0e-6, state.shared_effective_sheath_length_m);
    const double expected_field_v_per_m =
        (state.global_sheath_field_reference_potential_v - state.shared_patch_potential_v) /
        sheath_length_m;
    const double expected_charge_density_c_per_m3 =
        kEpsilon0 * expected_field_v_per_m / sheath_length_m;
    const double field_residual_v_per_m =
        std::abs(state.normal_electric_field_v_per_m - expected_field_v_per_m);
    const double reference_residual_v_per_m =
        std::abs(state.field_solver_reference_potential_v -
                 state.global_sheath_field_reference_potential_v) /
        sheath_length_m;
    const double charge_density_residual_v_per_m =
        std::abs(state.local_charge_density_c_per_m3 - expected_charge_density_c_per_m3) *
        sheath_length_m / std::max(1.0e-12, kEpsilon0);

    return std::max(
        {field_residual_v_per_m, reference_residual_v_per_m, charge_density_residual_v_per_m});
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalParticleFieldCoupledResidualV(
    const SurfaceModelRuntimeState& state) const
{
    if (!state.global_sheath_field_solve_active)
    {
        return 0.0;
    }

    const double sheath_length_m = std::max(1.0e-6, state.shared_effective_sheath_length_m);
    const double patch_area_m2 = std::max(1.0e-16, state.shared_patch_area_m2);
    const double charge_residual_v =
        state.global_particle_domain_charge_conservation_error_c * sheath_length_m /
        std::max(1.0e-18, kEpsilon0 * patch_area_m2);
    const double flux_residual_v =
        state.global_particle_domain_flux_conservation_error_a * sheath_length_m /
        std::max(1.0e-18, kEpsilon0 * patch_area_m2);
    const double field_residual_v =
        computeSharedSurfaceGlobalSheathFieldResidualVPerM(state) * sheath_length_m;

    return std::max(
        {field_residual_v, charge_residual_v, flux_residual_v});
}

double DensePlasmaSurfaceCharging::computeSharedSurfaceGlobalSheathFieldMultiStepStabilityMetricV()
    const
{
    if (history_shared_surface_global_coupled_solve_max_delta_v_.size() < 2)
    {
        return history_shared_surface_global_coupled_solve_max_delta_v_.empty()
                   ? 0.0
                   : history_shared_surface_global_coupled_solve_max_delta_v_.back();
    }

    const std::size_t count = std::min<std::size_t>(3, history_shared_surface_global_coupled_solve_max_delta_v_.size());
    double max_step_change_v = 0.0;
    for (std::size_t offset = 1; offset < count; ++offset)
    {
        const auto current =
            history_shared_surface_global_coupled_solve_max_delta_v_[history_shared_surface_global_coupled_solve_max_delta_v_.size() - offset];
        const auto previous =
            history_shared_surface_global_coupled_solve_max_delta_v_[history_shared_surface_global_coupled_solve_max_delta_v_.size() - offset - 1];
        max_step_change_v = std::max(max_step_change_v, std::abs(current - previous));
    }
    return max_step_change_v;
}

double DensePlasmaSurfaceCharging::computeLeakageCurrentDensity(double surface_potential_v) const
{
    if (config_.use_reference_current_balance &&
        (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
         config_.regime == SurfaceChargingRegime::GeoKineticPicLike))
    {
        const auto currents = computeSurfaceCurrents(surface_potential_v);
        return currents.conduction_current_a_per_m2;
    }

    const double thickness = std::max(1.0e-8, config_.dielectric_thickness_m);
    const double conductivity = computeEffectiveConductivity(surface_potential_v);
    return -conductivity * surface_potential_v / thickness;
}

double DensePlasmaSurfaceCharging::computeNetCurrentDensity(double surface_potential_v) const
{
    const auto currents = computeSurfaceCurrents(surface_potential_v);
    if (config_.use_reference_current_balance &&
        (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
         config_.regime == SurfaceChargingRegime::GeoKineticPicLike))
    {
        return currents.total_current_a_per_m2;
    }
    return currents.total_current_a_per_m2 + computeLeakageCurrentDensity(surface_potential_v);
}

double DensePlasmaSurfaceCharging::estimateCurrentDerivative(double surface_potential_v) const
{
    if (config_.use_reference_current_balance &&
        (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
         config_.regime == SurfaceChargingRegime::GeoKineticPicLike))
    {
        return computeSurfaceCurrents(surface_potential_v).current_derivative_a_per_m2_per_v;
    }

    const double step = std::max(1.0e-3, 1.0e-2 * std::max(1.0, std::abs(surface_potential_v)));
    const double jp = computeNetCurrentDensity(surface_potential_v + step);
    const double jm = computeNetCurrentDensity(surface_potential_v - step);
    return (jp - jm) / (2.0 * step);
}

double DensePlasmaSurfaceCharging::computeFloatingPotential() const
{
    if (config_.use_reference_current_balance &&
        (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
         config_.regime == SurfaceChargingRegime::GeoKineticPicLike) &&
        current_model_)
    {
        const double body_potential =
            initialized_ ? status_.body_potential_v : config_.body_initial_potential_v;
        const auto runtime_state =
            buildRuntimeState(body_potential, status_.patch_potential_v, config_.surface_area_m2,
                              computeEffectiveConductivity(
                                  status_.patch_potential_v,
                                  resolvePatchMaterial(
                                      config_, circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                                      circuit_model_
                                          ? circuit_model_->nodeName(
                                                circuit_model_->primaryPatchNodeIndex())
                                          : std::string{"patch"})),
                              computeEffectiveSheathLength(),
                              circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                              circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                                             : std::string{"patch"});
        return current_model_->solveEquilibriumPotential(
            ReferenceSurfaceRole::Patch, runtime_state,
            -std::max(config_.max_abs_potential_v,
                      20.0 * std::max(1.0, config_.plasma.electron_temperature_ev)),
            std::max(config_.max_abs_potential_v,
                     20.0 * std::max(1.0, config_.plasma.electron_temperature_ev)));
    }

    double left = -std::max(config_.max_abs_potential_v,
                            20.0 * std::max(1.0, config_.plasma.electron_temperature_ev));
    double right = -left;
    double f_left = computeNetCurrentDensity(left);
    double f_right = computeNetCurrentDensity(right);

    for (int expand = 0; expand < 8 && f_left * f_right > 0.0; ++expand)
    {
        left *= 1.5;
        right *= 1.5;
        f_left = computeNetCurrentDensity(left);
        f_right = computeNetCurrentDensity(right);
    }

    if (f_left * f_right > 0.0)
    {
        return std::abs(f_left) < std::abs(f_right) ? left : right;
    }

    for (int iteration = 0; iteration < 80; ++iteration)
    {
        const double mid = 0.5 * (left + right);
        const double f_mid = computeNetCurrentDensity(mid);
        if (std::abs(f_mid) < 1.0e-12 || std::abs(right - left) < 1.0e-3)
        {
            return mid;
        }

        if (f_left * f_mid <= 0.0)
        {
            right = mid;
            f_right = f_mid;
        }
        else
        {
            left = mid;
            f_left = f_mid;
        }
    }
    return 0.5 * (left + right);
}

double DensePlasmaSurfaceCharging::recommendTimeStep(double remaining_time_s, double minimum_dt_s,
                                                     double maximum_dt_s) const
{
    const double min_dt = std::max(1.0e-12, minimum_dt_s);
    const double max_dt = std::max(min_dt, maximum_dt_s);
    const double remaining = std::max(0.0, remaining_time_s);
    if (!initialized_ || remaining <= 0.0)
    {
        return min_dt;
    }
    if (remaining <= min_dt)
    {
        return remaining;
    }

    const double current_potential = status_.state.surface_potential_v;
    const double target_potential =
        config_.floating ? computeFloatingPotential() : current_potential;
    const double potential_error = std::abs(target_potential - current_potential);
    const double relative_error =
        potential_error / std::max(1.0, std::abs(target_potential));
    const double derivative =
        std::abs(estimateCurrentDerivative(current_potential));
    const double capacitance =
        std::max(1.0e-12, status_.state.capacitance_per_area_f_per_m2);
    const double relaxation_time =
        derivative > 1.0e-18 ? capacitance / derivative : max_dt;
    const double net_current = std::abs(computeNetCurrentDensity(current_potential));
    const double steady_current_threshold =
        1.0e-3 * std::max(1.0, config_.max_abs_current_density_a_per_m2);

    double suggested_dt = min_dt;
    if (relative_error < 5.0e-2 || net_current < steady_current_threshold)
    {
        suggested_dt = max_dt;
    }
    else if (relative_error < 2.5e-1)
    {
        suggested_dt = std::min(max_dt, 5.0 * relaxation_time);
    }
    else
    {
        suggested_dt = std::min(max_dt, std::max(min_dt, 0.5 * relaxation_time));
    }

    if (graph_capacitance_matrix_provider_ && circuit_model_)
    {
        const double graph_coupling_metric =
            graph_capacitance_matrix_provider_->graphCouplingMetric(config_,
                                                                    circuit_model_.get());
        const double coupling_penalty =
            std::clamp(graph_coupling_metric * 1.0e12, 0.0, 4.0);
        suggested_dt /= (1.0 + 0.5 * coupling_penalty);
    }

    const SurfaceModelRuntimeState runtime_state =
        buildRuntimeState(status_.body_potential_v, status_.patch_potential_v,
                          config_.surface_area_m2, 0.0,
                          computeEffectiveSheathLength());
    const double volume_penalty =
        1.0 + 0.35 * std::clamp(runtime_state.volume_mesh_coupling_gain, 0.0, 1.0) +
        0.05 * std::max(0.0, runtime_state.volume_projection_weight_sum - 1.0);
    suggested_dt /= volume_penalty;

    if (current_potential < -25.0)
    {
        const double negative_sensitivity =
            std::clamp(std::abs(current_potential) / 250.0, 0.0, 8.0);
        const double stiffness =
            derivative * std::max(min_dt, std::min(remaining, max_dt)) / capacitance;
        const double negative_region_penalty =
            1.0 + 0.35 * negative_sensitivity +
            0.30 * std::clamp(stiffness, 0.0, 8.0);
        suggested_dt /= negative_region_penalty;
    }

    return std::clamp(std::min(remaining, suggested_dt), min_dt, max_dt);
}

double DensePlasmaSurfaceCharging::advancePotentialImplicit(double surface_potential_v,
                                                           double capacitance_per_area_f_per_m2,
                                                           double dt) const
{
    const double dt_safe = std::max(1.0e-12, dt);
    const double capacitance = std::max(1.0e-12, capacitance_per_area_f_per_m2);
    const double search_limit =
        std::max(config_.max_abs_potential_v, 20.0 * std::max(1.0, config_.plasma.electron_temperature_ev));
    const double max_delta = transition_runtime_max_delta_potential_v_per_step_;

    const auto residual = [&](double potential) {
        return capacitance * (potential - surface_potential_v) / dt_safe -
               computeNetCurrentDensity(potential);
    };
    const auto guarded_didv = [&](double potential) {
        const double raw_didv = estimateCurrentDerivative(potential);
        if (!std::isfinite(raw_didv))
        {
            return 0.0;
        }

        const double denominator_guard = 0.95 * capacitance / dt_safe;
        const double slope_guard = std::max(1.0e-9, denominator_guard);
        return std::clamp(raw_didv, -slope_guard, slope_guard);
    };

    double potential = clampSigned(surface_potential_v, search_limit);
    double residual_value = residual(potential);
    double previous_potential = potential;
    double previous_residual = residual_value;

    double bracket_span =
        std::max(max_delta, 0.05 * std::max(1.0, std::abs(surface_potential_v)));
    if (surface_potential_v < -25.0)
    {
        bracket_span *= 1.5;
    }
    bracket_span = std::min(bracket_span, search_limit);

    double lower = clampSigned(surface_potential_v - bracket_span, search_limit);
    double upper = clampSigned(surface_potential_v + bracket_span, search_limit);
    if (lower > upper)
    {
        std::swap(lower, upper);
    }
    double f_lower = residual(lower);
    double f_upper = residual(upper);
    for (int expand = 0; expand < 10 && f_lower * f_upper > 0.0; ++expand)
    {
        bracket_span = std::min(search_limit, bracket_span * 1.8);
        lower = clampSigned(surface_potential_v - bracket_span, search_limit);
        upper = clampSigned(surface_potential_v + bracket_span, search_limit);
        if (lower > upper)
        {
            std::swap(lower, upper);
        }
        f_lower = residual(lower);
        f_upper = residual(upper);
    }

    const bool bracketed = (f_lower * f_upper <= 0.0);
    for (int iteration = 0; iteration < 28; ++iteration)
    {
        if (std::abs(residual_value) < 1.0e-9)
        {
            return potential;
        }

        double candidate = potential;
        bool have_candidate = false;

        const double derivative = capacitance / dt_safe - guarded_didv(potential);
        if (std::isfinite(derivative) && std::abs(derivative) > 1.0e-12)
        {
            candidate = potential - residual_value / derivative;
            have_candidate = std::isfinite(candidate);
        }

        if (!have_candidate &&
            std::abs(residual_value - previous_residual) > 1.0e-12)
        {
            candidate = potential -
                residual_value * (potential - previous_potential) /
                    (residual_value - previous_residual);
            have_candidate = std::isfinite(candidate);
        }

        if (!have_candidate || (bracketed && !(candidate > lower && candidate < upper)))
        {
            if (bracketed)
            {
                candidate = 0.5 * (lower + upper);
            }
            else
            {
                const double explicit_delta =
                    dt_safe * computeNetCurrentDensity(potential) / capacitance;
                candidate = potential + explicit_delta;
            }
        }

        candidate = clampSigned(candidate, search_limit);
        candidate =
            surface_potential_v +
            std::clamp(candidate - surface_potential_v, -max_delta, max_delta);
        candidate = clampSigned(candidate, search_limit);

        const double candidate_residual = residual(candidate);
        if (bracketed)
        {
            if (f_lower * candidate_residual <= 0.0)
            {
                upper = candidate;
                f_upper = candidate_residual;
            }
            else
            {
                lower = candidate;
                f_lower = candidate_residual;
            }
        }

        previous_potential = potential;
        previous_residual = residual_value;
        potential = candidate;
        residual_value = candidate_residual;
        if (std::abs(potential - previous_potential) < 1.0e-6)
        {
            return potential;
        }
    }

    if (bracketed)
    {
        const double bracket_mid = 0.5 * (lower + upper);
        const double bounded =
            surface_potential_v +
            std::clamp(bracket_mid - surface_potential_v, -max_delta, max_delta);
        return clampSigned(bounded, search_limit);
    }

    const double explicit_candidate =
        surface_potential_v +
        dt_safe * computeNetCurrentDensity(surface_potential_v) / capacitance;
    const double bounded_candidate =
        surface_potential_v +
        std::clamp(explicit_candidate - surface_potential_v, -max_delta, max_delta);
    return clampSigned(bounded_candidate, search_limit);
}

bool DensePlasmaSurfaceCharging::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }

    // Prepare stage: evaluate transition decisions before entering solve paths.
    const auto prepare_stage = [&]() {
        SurfaceAdvanceTransitionInput transition_input;
        transition_input.config = config_;
        transition_input.status = status_;
        transition_input.has_legacy_benchmark_replay = !legacy_benchmark_replay_.empty();
        transition_input.has_circuit_model = circuit_model_ != nullptr;
        transition_input.patch_count = circuit_model_ ? circuit_model_->patchCount() : 0;
        return evaluateSurfaceAdvanceTransition(transition_input);
    };
    const auto transition = prepare_stage();
    const auto apply_conductivity_evolution = [&]() {
        transition_conductivity_scale_ = 1.0;
        if (!transition.extended_events.conductivity_evolution_active)
        {
            status_.transition_conductivity_scale = transition_conductivity_scale_;
            return;
        }

        const double base_scale = std::max(
            0.0,
            config_.material.getScalarProperty("transition_conductivity_evolution_base_scale", 1.0));
        const double min_scale = std::max(
            1.0e-9,
            config_.material.getScalarProperty("transition_conductivity_evolution_min_scale",
                                              1.0e-3));
        const double max_scale = std::max(
            min_scale,
            config_.material.getScalarProperty("transition_conductivity_evolution_max_scale",
                                              1.0e3));
        const double time_slope_per_day = config_.material.getScalarProperty(
            "transition_conductivity_evolution_time_slope_per_day", 0.0);
        const double sun_flux_coupling = std::clamp(
            config_.material.getScalarProperty("transition_conductivity_evolution_sun_flux_coupling",
                                              0.0),
            -2.0, 2.0);

        const double elapsed_days = std::max(0.0, status_.time_s / 86400.0);
        const double time_factor = std::max(0.0, 1.0 + time_slope_per_day * elapsed_days);
        const double sun_factor = std::max(
            0.0, 1.0 + sun_flux_coupling * (transition.sun_flux_scale - 1.0));
        transition_conductivity_scale_ = std::clamp(base_scale * time_factor * sun_factor,
                                                    min_scale, max_scale);
        status_.transition_conductivity_scale = transition_conductivity_scale_;
    };
    apply_conductivity_evolution();

    const auto apply_source_flux_updater = [&]() {
        transition_source_flux_scale_ = 1.0;
        if (!transition.extended_events.source_flux_updater_active)
        {
            status_.transition_source_flux_scale = transition_source_flux_scale_;
            return;
        }

        const double base_scale = std::max(
            0.0,
            config_.material.getScalarProperty("transition_source_flux_updater_base_scale", 1.0));
        const double min_scale = std::max(
            1.0e-9,
            config_.material.getScalarProperty("transition_source_flux_updater_min_scale", 1.0e-3));
        const double max_scale = std::max(
            min_scale,
            config_.material.getScalarProperty("transition_source_flux_updater_max_scale", 1.0e3));
        const double time_slope_per_day = config_.material.getScalarProperty(
            "transition_source_flux_updater_time_slope_per_day", 0.0);
        const double sun_flux_coupling = std::clamp(
            config_.material.getScalarProperty("transition_source_flux_updater_sun_flux_coupling",
                                              0.0),
            -2.0, 2.0);

        const double elapsed_days = std::max(0.0, status_.time_s / 86400.0);
        const double time_factor = std::max(0.0, 1.0 + time_slope_per_day * elapsed_days);
        const double sun_factor = std::max(
            0.0, 1.0 + sun_flux_coupling * (transition.sun_flux_scale - 1.0));
        transition_source_flux_scale_ = std::clamp(base_scale * time_factor * sun_factor,
                                                   min_scale, max_scale);
        status_.transition_source_flux_scale = transition_source_flux_scale_;
    };
    apply_source_flux_updater();

    const auto apply_simulation_param_updater = [&]() {
        transition_simulation_param_scale_ = 1.0;
        transition_runtime_internal_substeps_ =
            std::max<std::size_t>(1, config_.internal_substeps);
        transition_runtime_max_delta_potential_v_per_step_ =
            std::max(1.0e-3, config_.max_delta_potential_v_per_step);
        transition_runtime_solver_relaxation_factor_ =
            std::clamp(config_.solver_config.relaxation_factor, 0.05, 1.5);

        if (transition.extended_events.simulation_param_updater_active)
        {
            const double base_scale = std::max(
                0.0,
                config_.material.getScalarProperty("transition_simulation_param_updater_base_scale",
                                                  1.0));
            const double min_scale = std::max(
                1.0e-9,
                config_.material.getScalarProperty("transition_simulation_param_updater_min_scale",
                                                  1.0e-3));
            const double max_scale = std::max(
                min_scale,
                config_.material.getScalarProperty("transition_simulation_param_updater_max_scale",
                                                  1.0e3));
            const double time_slope_per_day = config_.material.getScalarProperty(
                "transition_simulation_param_updater_time_slope_per_day", 0.0);
            const double sun_flux_coupling = std::clamp(
                config_.material.getScalarProperty(
                    "transition_simulation_param_updater_sun_flux_coupling", 0.0),
                -2.0, 2.0);

            const double elapsed_days = std::max(0.0, status_.time_s / 86400.0);
            const double time_factor =
                std::max(0.0, 1.0 + time_slope_per_day * elapsed_days);
            const double sun_factor = std::max(
                0.0, 1.0 + sun_flux_coupling * (transition.sun_flux_scale - 1.0));
            transition_simulation_param_scale_ =
                std::clamp(base_scale * time_factor * sun_factor, min_scale, max_scale);

            const double internal_substep_scale =
                std::max(1.0e-6,
                         config_.material.getScalarProperty(
                             "transition_simulation_param_updater_internal_substep_scale", 1.0));
            const double max_delta_scale =
                std::max(1.0e-6,
                         config_.material.getScalarProperty(
                             "transition_simulation_param_updater_max_delta_scale", 1.0));
            const double solver_relaxation_scale =
                std::max(1.0e-6,
                         config_.material.getScalarProperty(
                             "transition_simulation_param_updater_solver_relaxation_scale", 1.0));

            const double effective_substep_factor =
                transition_simulation_param_scale_ * internal_substep_scale;
            const double effective_max_delta_factor =
                transition_simulation_param_scale_ * max_delta_scale;
            const double effective_relaxation_factor =
                transition_simulation_param_scale_ * solver_relaxation_scale;

            const double base_substeps =
                static_cast<double>(std::max<std::size_t>(1, config_.internal_substeps));
            transition_runtime_internal_substeps_ = std::clamp<std::size_t>(
                static_cast<std::size_t>(std::llround(
                    std::max(1.0, base_substeps * effective_substep_factor))),
                1, 4096);
            transition_runtime_max_delta_potential_v_per_step_ = std::clamp(
                std::max(1.0e-3,
                         config_.max_delta_potential_v_per_step * effective_max_delta_factor),
                1.0e-3, 1.0e6);
            transition_runtime_solver_relaxation_factor_ = std::clamp(
                config_.solver_config.relaxation_factor * effective_relaxation_factor,
                0.05, 1.5);
        }

        field_volume_adaptive_relaxation_ =
            std::clamp(transition_runtime_solver_relaxation_factor_, 0.05, 1.0);
        status_.transition_simulation_param_scale = transition_simulation_param_scale_;
        status_.transition_runtime_internal_substeps =
            transition_runtime_internal_substeps_;
        status_.transition_runtime_max_delta_potential_v_per_step =
            transition_runtime_max_delta_potential_v_per_step_;
        status_.transition_runtime_solver_relaxation_factor =
            transition_runtime_solver_relaxation_factor_;
    };
    apply_simulation_param_updater();

    if (transition.use_legacy_benchmark_replay)
    {
        return advanceLegacyBenchmarkReplay(dt);
    }

    if (transition.use_reference_circuit)
    {
        bool recalibrated = false;
        if (transition.pic_recalibration_requested)
        {
            applyGeoPicCalibration();
            recalibrated = true;
        }

        const double reference_capacitance =
            std::max(1.0e-12, status_.state.capacitance_per_area_f_per_m2);
        const double reference_derivative =
            estimateCurrentDerivative(status_.patch_potential_v);
        Coupling::SurfaceAdaptiveSubstepInputs adaptive_substep_inputs;
        adaptive_substep_inputs.base_substeps = transition_runtime_internal_substeps_;
        adaptive_substep_inputs.reference_potential_v = status_.patch_potential_v;
        adaptive_substep_inputs.dt_s = dt;
        adaptive_substep_inputs.current_derivative_a_per_m2_per_v =
            reference_derivative;
        adaptive_substep_inputs.capacitance_per_area_f_per_m2 =
            reference_capacitance;
        const std::size_t substeps =
            Coupling::computeSurfaceAdaptiveSubstepCount(adaptive_substep_inputs);
        const double sub_dt = dt / static_cast<double>(substeps);
        auto state = status_.state;
        SurfaceCurrents currents;
        double leakage_current = 0.0;
        double net_current = 0.0;
        double branch_current_a = 0.0;
        std::vector<double> latest_branch_currents_a(
            circuit_model_ ? circuit_model_->branchCount() : 0, 0.0);
        double effective_conductivity =
            computeEffectiveConductivity(status_.patch_potential_v,
                                         resolvePatchMaterial(
                                             config_,
                                             circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                                             circuit_model_
                                                 ? circuit_model_->nodeName(
                                                       circuit_model_->primaryPatchNodeIndex())
                                                 : std::string{"patch"}));
        const double effective_sheath_length = computeEffectiveSheathLength();
        double body_potential = status_.body_potential_v;
        double patch_potential = status_.patch_potential_v;
        double latest_pre_global_solve_patch_spread_v = 0.0;
        double latest_patch_spread_reduction_v = 0.0;
        double latest_patch_spread_reduction_ratio = 0.0;
        double latest_shared_current_matrix_coupling_offdiag_entries = 0.0;
        double latest_shared_global_coupled_solve_iterations = 0.0;
        double latest_shared_global_coupled_solve_converged = 0.0;
        double latest_shared_global_coupled_solve_max_delta_v = 0.0;
        double latest_shared_live_pic_coupled_refresh_count = 0.0;
        double latest_shared_particle_transport_offdiag_entries = 0.0;
        double latest_shared_particle_transport_total_conductance_s = 0.0;
        double latest_shared_particle_transport_conservation_error_a_per_v = 0.0;
        if (circuit_model_ && circuit_node_potentials_.size() != circuit_model_->nodeCount())
        {
            circuit_node_potentials_.resize(circuit_model_->nodeCount(), 0.0);
            for (std::size_t node_index = 0; node_index < circuit_model_->nodeCount(); ++node_index)
            {
                circuit_node_potentials_[node_index] = circuit_model_->nodePotential(node_index);
            }
        }

        for (std::size_t step = 0; step < substeps; ++step)
        {
            effective_conductivity =
                computeEffectiveConductivity(patch_potential,
                                             resolvePatchMaterial(
                                                 config_, circuit_model_->primaryPatchNodeIndex(),
                                                 circuit_model_->nodeName(
                                                     circuit_model_->primaryPatchNodeIndex())));
            rebuildSurfaceCircuit(effective_conductivity);
            for (std::size_t node_index = 0; node_index < circuit_node_potentials_.size(); ++node_index)
            {
                circuit_model_->setNodePotential(node_index, circuit_node_potentials_[node_index]);
            }

            std::vector<double> iteration_node_potentials_v = circuit_node_potentials_;
            std::vector<double> solved_node_potentials_v = circuit_node_potentials_;
            std::vector<double> solved_branch_currents_a(
                circuit_model_ ? circuit_model_->branchCount() : 0, 0.0);
            const bool shared_global_coupled_solve_active =
                transition.shared_global_coupled_solve;
            const std::size_t shared_global_coupled_iteration_limit = std::max<std::size_t>(
                1, transition.shared_global_coupled_iteration_limit);
            const double shared_global_coupled_tolerance_v =
                computeSharedSurfaceGlobalCoupledToleranceV();
            const double runtime_shared_global_coupled_relaxation =
                config_.material.getScalarProperty("shared_surface_global_coupled_relaxation",
                                                   1.0) *
                transition_runtime_solver_relaxation_factor_;
            const double shared_global_coupled_relaxation = resolveSharedGlobalCoupledRelaxation(
                config_, runtime_shared_global_coupled_relaxation);
            const bool fixed_iteration_policy = transition.fixed_iteration_policy;
            const bool residual_guarded_policy = transition.residual_guarded_policy;
            const bool enforce_fixed_global_iterations =
                shared_global_coupled_solve_active && fixed_iteration_policy;
            bool shared_global_coupled_converged = !shared_global_coupled_solve_active;
            std::size_t shared_global_coupled_iterations_used = 0;
            double shared_global_coupled_max_delta_v = 0.0;
            std::size_t shared_current_matrix_coupling_offdiag_entries = 0;
            std::size_t shared_particle_transport_offdiag_entries = 0;
            double shared_particle_transport_total_conductance_s = 0.0;
            double shared_particle_transport_conservation_error_a_per_v = 0.0;
            const bool shared_live_pic_coupled_refresh_active =
                transition.shared_live_pic_coupled_refresh;
            const bool shared_particle_transport_coupling_active =
                transition.shared_particle_transport_coupling;
            const double shared_live_pic_coupled_refresh_threshold_v =
                computeSharedSurfaceLivePicCoupledRefreshThresholdV();
            const std::size_t shared_live_pic_coupled_refresh_max_count =
                static_cast<std::size_t>(std::llround(std::clamp(
                    config_.material.getScalarProperty(
                        "shared_surface_live_pic_coupled_refresh_max_count", 1.0),
                    0.0, 8.0)));
            double previous_shared_refresh_potential_v = computeSharedSurfacePatchPotentialV();
            std::size_t shared_live_pic_coupled_refresh_count = 0;

            for (std::size_t coupled_iteration = 0;
                 coupled_iteration < shared_global_coupled_iteration_limit;
                 ++coupled_iteration)
            {
                circuit_node_potentials_ = iteration_node_potentials_v;
                for (std::size_t node_index = 0; node_index < iteration_node_potentials_v.size();
                     ++node_index)
                {
                    circuit_model_->setNodePotential(node_index, iteration_node_potentials_v[node_index]);
                }

                if (shared_live_pic_coupled_refresh_active &&
                    shared_live_pic_coupled_refresh_count <
                        shared_live_pic_coupled_refresh_max_count)
                {
                    const double shared_refresh_potential_v = computeSharedSurfacePatchPotentialV();
                    const bool needs_shared_refresh =
                        coupled_iteration == 0 ||
                        std::abs(shared_refresh_potential_v - previous_shared_refresh_potential_v) >=
                            shared_live_pic_coupled_refresh_threshold_v;
                    if (needs_shared_refresh)
                    {
                        const auto primary_patch_node_index =
                            circuit_model_->primaryPatchNodeIndex();
                        const auto shared_runtime_state =
                            buildRuntimeState(body_potential, shared_refresh_potential_v,
                                              computeSharedSurfacePatchAreaM2(),
                                              computeEffectiveConductivity(
                                                  shared_refresh_potential_v,
                                                  resolvePatchMaterial(
                                                      config_, primary_patch_node_index,
                                                      circuit_model_->nodeName(
                                                          primary_patch_node_index))),
                                              computeEffectiveSheathLength(),
                                              primary_patch_node_index,
                                              circuit_model_->nodeName(primary_patch_node_index));
                        if (current_model_->recalibrate(shared_runtime_state,
                                                        config_.pic_calibration_samples))
                        {
                            ++shared_live_pic_coupled_refresh_count;
                        }
                        previous_shared_refresh_potential_v = shared_refresh_potential_v;
                    }
                }

                if (shared_particle_transport_coupling_active)
                {
                    const auto primary_patch_node_index = circuit_model_->primaryPatchNodeIndex();
                    const double shared_patch_potential_v = computeSharedSurfacePatchPotentialV();
                    const double shared_patch_area_m2 = computeSharedSurfacePatchAreaM2();
                    const double shared_effective_sheath_length_m =
                        computeSharedSurfaceEffectiveSheathLengthM(computeEffectiveSheathLength());
                    const auto shared_runtime_state =
                        buildRuntimeState(body_potential, shared_patch_potential_v,
                                          shared_patch_area_m2,
                                          computeEffectiveConductivity(
                                              shared_patch_potential_v,
                                              resolvePatchMaterial(
                                                  config_, primary_patch_node_index,
                                                  circuit_model_->nodeName(primary_patch_node_index))),
                                          shared_effective_sheath_length_m, primary_patch_node_index,
                                          circuit_model_->nodeName(primary_patch_node_index));
                    const PicMccCurrentSample shared_live_pic_sample =
                        current_model_ ? current_model_->latestLivePicSample(shared_runtime_state)
                                       : PicMccCurrentSample{};
                    const double shared_transport_conductance_s_per_m2 =
                        computeSharedSurfaceParticleTransportConductanceSPerM2(
                            shared_live_pic_sample.net_collection_current_density_a_per_m2,
                            shared_live_pic_sample.current_derivative_a_per_m2_per_v);
                    advanceSharedSurfaceParticleTransportState(
                        sub_dt / static_cast<double>(shared_global_coupled_iteration_limit),
                        shared_live_pic_sample, shared_patch_area_m2,
                        shared_effective_sheath_length_m);
                    updateSharedSurfaceDistributedParticleTransportState(
                        iteration_node_potentials_v,
                        sub_dt / static_cast<double>(shared_global_coupled_iteration_limit),
                        shared_patch_area_m2,
                        shared_effective_sheath_length_m,
                        shared_transport_conductance_s_per_m2);
                }
                else
                {
                    global_particle_domain_state_ = GlobalParticleDomainState{};
                    global_particle_domain_state_.bookkeeping_mode =
                        "owned_global_particle_domain_state_v2";
                    global_particle_repository_state_ = GlobalParticleRepositoryState{};
                    global_particle_repository_state_.bookkeeping_mode =
                        "owned_global_particle_repository_state_v1";
                    global_particle_repository_state_.lifecycle_mode =
                        "global_particle_transport_reservoir_lifecycle_v1";
                }

                solveOwnedGlobalSheathFieldSystem(
                    iteration_node_potentials_v,
                    computeSharedSurfaceEffectiveSheathLengthM(computeEffectiveSheathLength()),
                    sub_dt / static_cast<double>(shared_global_coupled_iteration_limit));

                Coupling::SurfaceCircuitLinearization linearization;
                linearization.plasma_current_a =
                    std::vector<double>(circuit_model_->nodeCount(), 0.0);
                linearization.plasma_didv_a_per_v =
                    std::vector<double>(circuit_model_->nodeCount(), 0.0);
                linearization.additional_rhs_a =
                    std::vector<double>(circuit_model_->nodeCount(), 0.0);
                linearization.additional_diagonal_a_per_v =
                    std::vector<double>(circuit_model_->nodeCount(), 0.0);
                const auto mapping_state = buildSurfaceCircuitMappingState();
                std::vector<double> shared_patch_voltage_weights(circuit_model_->nodeCount(), 0.0);
                std::vector<std::size_t> patch_node_indices;
                patch_node_indices.reserve(circuit_model_->patchCount());
                for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount();
                     ++patch_ordinal)
                {
                    patch_node_indices.push_back(circuit_model_->patchNodeIndex(patch_ordinal));
                }
                double shared_patch_area_m2 = 0.0;
                for (const auto& group : mapping_state.reduced_node_groups)
                {
                    bool patch_group = false;
                    for (const auto member_node_index : group.member_node_indices)
                    {
                        if (member_node_index < circuit_model_->nodeCount() &&
                            circuit_model_->nodeIsPatch(member_node_index))
                        {
                            patch_group = true;
                            break;
                        }
                    }
                    if (patch_group)
                    {
                        shared_patch_area_m2 += std::max(1.0e-16, group.total_area_m2);
                    }
                }
                if (shared_patch_area_m2 > 0.0 && mapping_state.valid)
                {
                    for (const auto& group : mapping_state.reduced_node_groups)
                    {
                        const double reduced_weight =
                            std::max(1.0e-16, group.total_area_m2) / shared_patch_area_m2;
                        for (const auto member_node_index : group.member_node_indices)
                        {
                            if (member_node_index < circuit_model_->nodeCount() &&
                                circuit_model_->nodeIsPatch(member_node_index))
                            {
                                shared_patch_voltage_weights[member_node_index] = reduced_weight;
                            }
                        }
                    }
                }
                shared_current_matrix_coupling_offdiag_entries = 0;
                shared_particle_transport_offdiag_entries = 0;
                shared_particle_transport_total_conductance_s = 0.0;
                shared_particle_transport_conservation_error_a_per_v = 0.0;

                for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount();
                     ++patch_ordinal)
                {
                    const auto patch_node_index = circuit_model_->patchNodeIndex(patch_ordinal);
                    const double local_patch_potential =
                        circuit_model_->nodePotential(patch_node_index);
                    const double patch_area_m2 =
                        std::max(1.0e-16, circuit_model_->nodeAreaM2(patch_node_index));
                    const double local_effective_conductivity =
                        computeEffectiveConductivity(
                            local_patch_potential,
                            resolvePatchMaterial(config_, patch_node_index,
                                                 circuit_model_->nodeName(patch_node_index)));
                    const auto patch_branch_index =
                        circuit_model_->patchToBodyBranchIndex(patch_ordinal);
                    if (patch_branch_index != std::numeric_limits<std::size_t>::max() &&
                        circuit_model_->branchUsesDynamicConductance(patch_branch_index))
                    {
                        circuit_model_->setBranchConductance(
                            patch_branch_index,
                            std::max(0.0, local_effective_conductivity) * patch_area_m2 /
                                std::max(1.0e-9, config_.dielectric_thickness_m));
                    }
                    const auto patch_runtime_state =
                        buildRuntimeState(body_potential, local_patch_potential, patch_area_m2,
                                          local_effective_conductivity, effective_sheath_length,
                                          patch_node_index, circuit_model_->nodeName(patch_node_index));
                    const double local_capacitance_per_area =
                        computeCapacitancePerArea(patch_runtime_state);
                    circuit_model_->setNodeCapacitance(
                        patch_node_index,
                        std::max(0.0, local_capacitance_per_area * patch_area_m2));
                    const auto patch_currents_old =
                        current_model_->evaluate(ReferenceSurfaceRole::Patch, patch_runtime_state);
                    const auto patch_didv = assembleSurfaceDidvState(
                        ReferenceSurfaceRole::Patch, patch_runtime_state, patch_currents_old);
                    linearization.plasma_current_a[patch_node_index] =
                        patch_didv.valid ? patch_didv.plasma_current_a : 0.0;
                    const bool fallback_shared_current_matrix_coupling_enabled =
                        !global_sheath_field_solve_state_.active &&
                        patch_runtime_state.shared_runtime_enabled &&
                        circuit_model_->patchCount() >= 2 && shared_patch_area_m2 > 0.0;
                    const auto shared_matrix_fallback =
                        Coupling::computeSurfaceSharedCurrentMatrixFallback(
                            fallback_shared_current_matrix_coupling_enabled,
                            patch_didv.valid,
                            patch_node_index,
                            patch_didv.plasma_didv_a_per_v,
                            patch_didv.conduction_didv_a_per_v,
                            shared_patch_voltage_weights,
                            patch_node_indices,
                            mapping_state.circuit_to_reduced_node_index);
                    linearization.plasma_didv_a_per_v[patch_node_index] =
                        shared_matrix_fallback.plasma_didv_a_per_v;
                    shared_current_matrix_coupling_offdiag_entries +=
                        static_cast<std::size_t>(
                            shared_matrix_fallback.off_diagonal_entries.size());
                    linearization.additional_off_diagonal_entries.insert(
                        linearization.additional_off_diagonal_entries.end(),
                        shared_matrix_fallback.off_diagonal_entries.begin(),
                        shared_matrix_fallback.off_diagonal_entries.end());
                }

                appendOwnedGlobalCoupledLinearization(
                    iteration_node_potentials_v,
                    sub_dt / static_cast<double>(shared_global_coupled_iteration_limit),
                    linearization, shared_current_matrix_coupling_offdiag_entries,
                    shared_particle_transport_offdiag_entries,
                    shared_particle_transport_total_conductance_s,
                    shared_particle_transport_conservation_error_a_per_v);

                if (config_.body_floating)
                {
                    const auto runtime_state =
                        buildRuntimeState(body_potential, patch_potential,
                                          std::max(1.0e-16, circuit_model_->bodyAreaM2()),
                                          effective_conductivity, effective_sheath_length,
                                          circuit_model_->bodyNodeIndex(),
                                          circuit_model_->nodeName(circuit_model_->bodyNodeIndex()));
                    circuit_model_->setNodeCapacitance(circuit_model_->bodyNodeIndex(),
                                                       computeBodyNodeCapacitanceF(runtime_state));
                    const auto body_terms_old =
                        current_model_->evaluate(ReferenceSurfaceRole::Body, runtime_state);
                    const auto body_didv = assembleSurfaceDidvState(
                        ReferenceSurfaceRole::Body, runtime_state, body_terms_old);
                    linearization.plasma_current_a[circuit_model_->bodyNodeIndex()] =
                        body_didv.valid ? body_didv.total_current_a : 0.0;
                    linearization.plasma_didv_a_per_v[circuit_model_->bodyNodeIndex()] =
                        body_didv.valid ? body_didv.total_didv_a_per_v : 0.0;
                }

                if (graph_capacitance_matrix_provider_ != nullptr && sub_dt > 0.0)
                {
                    const double runtime_matrix_weight = std::clamp(
                        config_.material.getScalarProperty("graph_matrix_runtime_weight", 0.35), 0.0,
                        2.0);
                    if (runtime_matrix_weight > 0.0)
                    {
                        linearization.additional_off_diagonal_entries.reserve(
                            linearization.additional_off_diagonal_entries.size() +
                            2 * circuit_model_->branchCount());
                        for (std::size_t branch_index = 0;
                             branch_index < circuit_model_->branchCount();
                             ++branch_index)
                        {
                            const auto from_node =
                                circuit_model_->branchFromNodeIndex(branch_index);
                            const auto to_node = circuit_model_->branchToNodeIndex(branch_index);
                            if (from_node >= circuit_model_->nodeCount() ||
                                to_node >= circuit_model_->nodeCount())
                            {
                                continue;
                            }
                            const double mutual_capacitance_f = std::max(
                                0.0, runtime_matrix_weight *
                                         std::abs(graph_capacitance_matrix_provider_->mutualCapacitanceF(
                                             config_, circuit_model_.get(), branch_index)));
                            const double mutual_conductance_a_per_v = mutual_capacitance_f / sub_dt;
                            if (!std::isfinite(mutual_conductance_a_per_v) ||
                                mutual_conductance_a_per_v <= 0.0)
                            {
                                continue;
                            }

                            const double from_potential_old_v =
                                circuit_model_->nodePotential(from_node);
                            const double to_potential_old_v =
                                circuit_model_->nodePotential(to_node);
                            const double old_delta_v =
                                from_potential_old_v - to_potential_old_v;

                            linearization.additional_diagonal_a_per_v[from_node] +=
                                mutual_conductance_a_per_v;
                            linearization.additional_diagonal_a_per_v[to_node] +=
                                mutual_conductance_a_per_v;
                            linearization.additional_rhs_a[from_node] +=
                                mutual_conductance_a_per_v * old_delta_v;
                            linearization.additional_rhs_a[to_node] -=
                                mutual_conductance_a_per_v * old_delta_v;

                            Coupling::SurfaceCircuitLinearization::OffDiagonalEntry from_to{};
                            from_to.row_node = from_node;
                            from_to.column_node = to_node;
                            from_to.coefficient_a_per_v = -mutual_conductance_a_per_v;
                            linearization.additional_off_diagonal_entries.push_back(from_to);

                            Coupling::SurfaceCircuitLinearization::OffDiagonalEntry to_from{};
                            to_from.row_node = to_node;
                            to_from.column_node = from_node;
                            to_from.coefficient_a_per_v = -mutual_conductance_a_per_v;
                            linearization.additional_off_diagonal_entries.push_back(to_from);
                        }
                    }
                }

                Coupling::CircuitDidvCompositionInput didv_input;
                didv_input.nodes.resize(circuit_model_->nodeCount());
                didv_input.off_diagonal_entries = linearization.additional_off_diagonal_entries;
                 Coupling::CircuitKernelTopologyState topology_state;
                 topology_state.node_capacitance_f.resize(circuit_model_->nodeCount(), 0.0);
                 topology_state.node_shunt_conductance_s.resize(circuit_model_->nodeCount(), 0.0);
                for (std::size_t node_index = 0; node_index < circuit_model_->nodeCount();
                     ++node_index)
                {
                    auto& didv_node = didv_input.nodes[node_index];
                    didv_node.terms.push_back(Coupling::CircuitDidvTerm{
                        "plasma", linearization.plasma_current_a[node_index],
                        linearization.plasma_didv_a_per_v[node_index]});
                    didv_node.additional_rhs_a = linearization.additional_rhs_a[node_index];
                    didv_node.additional_diagonal_a_per_v =
                        linearization.additional_diagonal_a_per_v[node_index];
                    topology_state.node_capacitance_f[node_index] =
                        circuit_model_->nodeIsPatch(node_index)
                            ? computeCapacitancePerArea(buildRuntimeState(
                                  body_potential, circuit_model_->nodePotential(node_index),
                                  std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index)),
                                  computeEffectiveConductivity(
                                      circuit_model_->nodePotential(node_index),
                                      resolvePatchMaterial(config_, node_index,
                                                           circuit_model_->nodeName(node_index))),
                                  effective_sheath_length, node_index,
                                  circuit_model_->nodeName(node_index))) *
                                  std::max(1.0e-16, circuit_model_->nodeAreaM2(node_index))
                            : computeBodyNodeCapacitanceF(buildRuntimeState(
                                  body_potential, patch_potential,
                                  std::max(1.0e-16, circuit_model_->bodyAreaM2()),
                                  effective_conductivity, effective_sheath_length, node_index,
                                  circuit_model_->nodeName(node_index)));
                    topology_state.node_shunt_conductance_s[node_index] =
                        circuit_model_->nodeIsPatch(node_index) ? 0.0 : config_.body_leakage_conductance_s;
                }

                topology_state.branch_conductance_s.resize(circuit_model_->branchCount(), 0.0);
                topology_state.branch_bias_v.resize(circuit_model_->branchCount(), 0.0);
                topology_state.branch_current_source_a.resize(circuit_model_->branchCount(), 0.0);
                for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount();
                     ++branch_index)
                {
                    topology_state.branch_conductance_s[branch_index] =
                        circuit_model_->branchConductanceS(branch_index);
                    topology_state.branch_bias_v[branch_index] =
                        circuit_model_->branchBiasV(branch_index);
                }

                Coupling::SurfaceCircuitKernelInput kernel_input =
                    Coupling::composeCircuitKernelInput(didv_input, topology_state);

                const auto circuit_result = circuit_model_->advanceImplicit(
                    sub_dt, kernel_input,
                    transition_runtime_max_delta_potential_v_per_step_);
                if (!circuit_result.converged)
                {
                    return false;
                }

                solved_node_potentials_v = circuit_result.node_potentials_v;
                if (shared_global_coupled_solve_active && shared_global_coupled_relaxation < 1.0)
                {
                    double blended_max_delta_v = 0.0;
                    for (std::size_t node_index = 0;
                         node_index < solved_node_potentials_v.size() &&
                         node_index < iteration_node_potentials_v.size();
                         ++node_index)
                    {
                        const double blended_potential_v =
                            (1.0 - shared_global_coupled_relaxation) *
                                iteration_node_potentials_v[node_index] +
                            shared_global_coupled_relaxation *
                                solved_node_potentials_v[node_index];
                        blended_max_delta_v = std::max(
                            blended_max_delta_v,
                            std::abs(blended_potential_v - iteration_node_potentials_v[node_index]));
                        solved_node_potentials_v[node_index] = blended_potential_v;
                    }
                    shared_global_coupled_max_delta_v = blended_max_delta_v;
                    solved_branch_currents_a.assign(circuit_model_->branchCount(), 0.0);
                    for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount();
                         ++branch_index)
                    {
                        const auto from_node_index =
                            circuit_model_->branchFromNodeIndex(branch_index);
                        const auto to_node_index =
                            circuit_model_->branchToNodeIndex(branch_index);
                        const double from_potential_v =
                            from_node_index < solved_node_potentials_v.size()
                                ? solved_node_potentials_v[from_node_index]
                                : 0.0;
                        const double to_potential_v =
                            to_node_index < solved_node_potentials_v.size()
                                ? solved_node_potentials_v[to_node_index]
                                : 0.0;
                        solved_branch_currents_a[branch_index] =
                            circuit_model_->branchConductanceS(branch_index) *
                            (from_potential_v - to_potential_v);
                    }
                }
                else
                {
                    solved_branch_currents_a = circuit_result.branch_currents_a;
                    shared_global_coupled_max_delta_v = circuit_result.max_delta_potential_v;
                }
                shared_global_coupled_iterations_used = coupled_iteration + 1;

                const bool reached_explicit_global_residual_target =
                    global_sheath_field_solve_state_.active &&
                    global_sheath_field_solve_state_.field_residual_v_per_m <=
                        std::max(1.0e-6, shared_global_coupled_tolerance_v /
                                             std::max(1.0e-6, computeSharedSurfaceEffectiveSheathLengthM(
                                                                      computeEffectiveSheathLength()))) &&
                    global_sheath_field_solve_state_.particle_field_coupled_residual_v <=
                        std::max(1.0e-3, shared_global_coupled_tolerance_v) &&
                    global_particle_domain_state_.charge_conservation_error_c <= 1.0e-15 &&
                    global_particle_domain_state_.flux_conservation_error_a <= 1.0e-15 &&
                    coupled_iteration + 1 >= 2;
                const bool reached_coupled_fixed_point =
                    (shared_global_coupled_max_delta_v <= shared_global_coupled_tolerance_v &&
                     (!shared_global_coupled_solve_active || coupled_iteration + 1 >= 2)) ||
                    reached_explicit_global_residual_target;
                if (!shared_global_coupled_solve_active)
                {
                    shared_global_coupled_converged = true;
                    break;
                }
                if (enforce_fixed_global_iterations)
                {
                    if (coupled_iteration + 1 >= shared_global_coupled_iteration_limit)
                    {
                        shared_global_coupled_converged =
                            reached_coupled_fixed_point || !residual_guarded_policy;
                        break;
                    }
                }
                else if (reached_coupled_fixed_point)
                {
                    shared_global_coupled_converged = true;
                    break;
                }

                iteration_node_potentials_v = circuit_result.node_potentials_v;
            }

            const auto compute_patch_spread_from_potentials =
                [this](const std::vector<double>& node_potentials_v) {
                    if (!circuit_model_ || circuit_model_->patchCount() < 2)
                    {
                        return 0.0;
                    }

                    double weighted_sum_v = 0.0;
                    double total_area_m2 = 0.0;
                    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount();
                         ++patch_ordinal)
                    {
                        const auto patch_node_index = circuit_model_->patchNodeIndex(patch_ordinal);
                        if (patch_node_index >= node_potentials_v.size())
                        {
                            continue;
                        }
                        const double area_m2 =
                            std::max(1.0e-16, circuit_model_->nodeAreaM2(patch_node_index));
                        weighted_sum_v += area_m2 * node_potentials_v[patch_node_index];
                        total_area_m2 += area_m2;
                    }

                    if (total_area_m2 <= 0.0)
                    {
                        return 0.0;
                    }

                    const double shared_patch_potential_v = weighted_sum_v / total_area_m2;
                    double weighted_spread_v = 0.0;
                    for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount();
                         ++patch_ordinal)
                    {
                        const auto patch_node_index = circuit_model_->patchNodeIndex(patch_ordinal);
                        if (patch_node_index >= node_potentials_v.size())
                        {
                            continue;
                        }
                        const double area_m2 =
                            std::max(1.0e-16, circuit_model_->nodeAreaM2(patch_node_index));
                        weighted_spread_v +=
                            area_m2 *
                            std::abs(node_potentials_v[patch_node_index] - shared_patch_potential_v);
                    }
                    return weighted_spread_v / total_area_m2;
                };
            const double pre_global_solve_patch_spread_v =
                compute_patch_spread_from_potentials(solved_node_potentials_v);
            const double shared_global_solve_weight =
                global_sheath_field_solve_state_.active ? 0.0 : computeSharedSurfaceGlobalSolveWeight();
            if (shared_global_solve_weight > 0.0 && circuit_model_->patchCount() >= 2)
            {
                double shared_patch_potential_v = 0.0;
                double shared_patch_area_m2 = 0.0;
                for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount();
                     ++patch_ordinal)
                {
                    const auto patch_node_index = circuit_model_->patchNodeIndex(patch_ordinal);
                    if (patch_node_index >= solved_node_potentials_v.size())
                    {
                        continue;
                    }
                    const double area_m2 =
                        std::max(1.0e-16, circuit_model_->nodeAreaM2(patch_node_index));
                    shared_patch_potential_v += area_m2 * solved_node_potentials_v[patch_node_index];
                    shared_patch_area_m2 += area_m2;
                }
                shared_patch_potential_v /=
                    std::max(1.0e-16, shared_patch_area_m2);

                for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount();
                     ++patch_ordinal)
                {
                    const auto patch_node_index = circuit_model_->patchNodeIndex(patch_ordinal);
                    if (patch_node_index >= solved_node_potentials_v.size())
                    {
                        continue;
                    }
                    solved_node_potentials_v[patch_node_index] =
                        (1.0 - shared_global_solve_weight) *
                            solved_node_potentials_v[patch_node_index] +
                        shared_global_solve_weight * shared_patch_potential_v;
                }

                solved_branch_currents_a.assign(circuit_model_->branchCount(), 0.0);
                for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount();
                     ++branch_index)
                {
                    const auto from_node_index =
                        circuit_model_->branchFromNodeIndex(branch_index);
                    const auto to_node_index =
                        circuit_model_->branchToNodeIndex(branch_index);
                    const double from_potential_v =
                        from_node_index < solved_node_potentials_v.size()
                            ? solved_node_potentials_v[from_node_index]
                            : 0.0;
                    const double to_potential_v =
                        to_node_index < solved_node_potentials_v.size()
                            ? solved_node_potentials_v[to_node_index]
                            : 0.0;
                    solved_branch_currents_a[branch_index] =
                        circuit_model_->branchConductanceS(branch_index) *
                        (from_potential_v - to_potential_v);
                }
            }
    const double post_global_solve_patch_spread_v =
        compute_patch_spread_from_potentials(solved_node_potentials_v);
            if (shared_global_coupled_solve_active &&
                post_global_solve_patch_spread_v <=
                    std::max(shared_global_coupled_tolerance_v, 1.0e-6))
            {
                shared_global_coupled_converged = true;
                shared_global_coupled_max_delta_v = std::min(
                    shared_global_coupled_max_delta_v, post_global_solve_patch_spread_v);
            }
            const double patch_spread_reduction_v = std::max(
                0.0, pre_global_solve_patch_spread_v - post_global_solve_patch_spread_v);
            const double patch_spread_reduction_ratio =
                pre_global_solve_patch_spread_v > 1.0e-16
                    ? patch_spread_reduction_v / pre_global_solve_patch_spread_v
                    : 0.0;
            const bool physically_converged_shared_global_solve =
                shared_global_coupled_solve_active &&
                ((global_sheath_field_solve_state_.active &&
                  global_sheath_field_solve_state_.field_residual_v_per_m <=
                      std::max(1.0e-6, shared_global_coupled_tolerance_v /
                                           std::max(1.0e-6, computeSharedSurfaceEffectiveSheathLengthM(
                                                                    computeEffectiveSheathLength()))) &&
                  global_sheath_field_solve_state_.particle_field_coupled_residual_v <=
                      std::max(1.0e-3, shared_global_coupled_tolerance_v) &&
                  global_particle_domain_state_.charge_conservation_error_c <= 1.0e-15 &&
                  global_particle_domain_state_.flux_conservation_error_a <= 1.0e-15) ||
                 post_global_solve_patch_spread_v <=
                     std::max(shared_global_coupled_tolerance_v, 1.0e-6));
            latest_pre_global_solve_patch_spread_v = pre_global_solve_patch_spread_v;
            latest_patch_spread_reduction_v = patch_spread_reduction_v;
            latest_patch_spread_reduction_ratio = patch_spread_reduction_ratio;
            latest_shared_current_matrix_coupling_offdiag_entries =
                static_cast<double>(shared_current_matrix_coupling_offdiag_entries);
            latest_shared_global_coupled_solve_iterations =
                static_cast<double>(shared_global_coupled_iterations_used);
            latest_shared_global_coupled_solve_converged =
                (shared_global_coupled_converged || physically_converged_shared_global_solve)
                    ? 1.0
                    : 0.0;
            latest_shared_global_coupled_solve_max_delta_v =
                physically_converged_shared_global_solve
                    ? std::min(shared_global_coupled_max_delta_v,
                               post_global_solve_patch_spread_v)
                    : shared_global_coupled_max_delta_v;
            latest_shared_live_pic_coupled_refresh_count =
                static_cast<double>(shared_live_pic_coupled_refresh_count);
            latest_shared_particle_transport_offdiag_entries =
                static_cast<double>(shared_particle_transport_offdiag_entries);
            latest_shared_particle_transport_total_conductance_s =
                shared_particle_transport_total_conductance_s;
            latest_shared_particle_transport_conservation_error_a_per_v =
                shared_particle_transport_conservation_error_a_per_v;
            shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_ =
                std::max(
                    shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_,
                    std::clamp(
                        config_.material.getScalarProperty(
                            "shared_surface_particle_transport_edge_graph_operator_conductance_weight",
                            0.002),
                        0.0, 10.0) *
                        (sub_dt / static_cast<double>(std::max<std::size_t>(
                            1, shared_global_coupled_iteration_limit))) *
                        shared_particle_transport_total_conductance_s);

            body_potential = solved_node_potentials_v[circuit_model_->bodyNodeIndex()];
            patch_potential =
                solved_node_potentials_v[circuit_model_->primaryPatchNodeIndex()];
            circuit_node_potentials_ = solved_node_potentials_v;
            if (history_surface_branch_currents_.size() != solved_branch_currents_a.size())
            {
                history_surface_branch_currents_.assign(solved_branch_currents_a.size(), {});
            }
            latest_branch_currents_a = solved_branch_currents_a;
            if (circuit_model_->primaryPatchToBodyBranchIndex() <
                solved_branch_currents_a.size())
            {
                branch_current_a = solved_branch_currents_a
                    [circuit_model_->primaryPatchToBodyBranchIndex()];
            }

            currents = computeSurfaceCurrents(patch_potential);
            currents.conduction_current_a_per_m2 =
                branch_current_a /
                std::max(1.0e-16, circuit_model_->nodeAreaM2(circuit_model_->primaryPatchNodeIndex()));
            currents.total_current_a_per_m2 =
                currents.electron_current_a_per_m2 + currents.ion_current_a_per_m2 +
                currents.secondary_emission_a_per_m2 + currents.ion_secondary_emission_a_per_m2 +
                currents.backscatter_emission_a_per_m2 + currents.photo_emission_a_per_m2 +
                currents.thermionic_emission_a_per_m2 + currents.field_emission_a_per_m2 +
                currents.conduction_current_a_per_m2 + currents.ram_ion_current_a_per_m2;

            leakage_current = currents.conduction_current_a_per_m2;
            net_current = currents.total_current_a_per_m2;

            state.surface_potential_v = patch_potential;
            state.capacitance_per_area_f_per_m2 =
                computeCapacitancePerArea(buildRuntimeState(
                    body_potential, patch_potential,
                    std::max(1.0e-16,
                             circuit_model_->nodeAreaM2(circuit_model_->primaryPatchNodeIndex())),
                    effective_conductivity, effective_sheath_length,
                    circuit_model_->primaryPatchNodeIndex(),
                    circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())));
            state.surface_charge_density_c_per_m2 =
                state.capacitance_per_area_f_per_m2 * (patch_potential - body_potential);
            state.stored_energy_j_per_m2 =
                0.5 * state.capacitance_per_area_f_per_m2 * (patch_potential - body_potential) *
                (patch_potential - body_potential);

            if (!std::isfinite(state.surface_potential_v) ||
                !std::isfinite(state.surface_charge_density_c_per_m2))
            {
                return false;
            }
        }

        const auto commit_reference_stage = [&]() -> bool {
            const double target_potential =
                config_.floating ? computeFloatingPotential() : patch_potential;
            status_.currents = currents;
            status_.state = state;
            status_.body_potential_v = body_potential;
            status_.patch_potential_v = patch_potential;
            status_.circuit_branch_current_a = branch_current_a;
            status_.pic_recalibrated = recalibrated;
            status_.equilibrium_error = status_.patch_potential_v - target_potential;
            status_.equilibrium_reached = std::abs(status_.equilibrium_error) < 1.0e-2;
            status_.time_s += dt;
            status_.steps_completed += 1;
            const auto final_runtime_state =
                buildRuntimeState(status_.body_potential_v, status_.patch_potential_v,
                                  std::max(1.0e-16,
                                           circuit_model_->nodeAreaM2(
                                               circuit_model_->primaryPatchNodeIndex())),
                                  effective_conductivity, effective_sheath_length,
                                  circuit_model_->primaryPatchNodeIndex(),
                                  circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex()));

            history_time_.push_back(status_.time_s);
            history_potential_.push_back(status_.state.surface_potential_v);
            history_charge_.push_back(status_.state.surface_charge_density_c_per_m2);
            history_current_.push_back(status_.currents.total_current_a_per_m2);
            history_electron_current_.push_back(status_.currents.electron_current_a_per_m2);
            history_ion_current_.push_back(status_.currents.ion_current_a_per_m2);
            history_secondary_current_.push_back(status_.currents.secondary_emission_a_per_m2);
            history_ion_secondary_current_.push_back(
                status_.currents.ion_secondary_emission_a_per_m2);
            history_backscatter_current_.push_back(status_.currents.backscatter_emission_a_per_m2);
            history_photo_current_.push_back(status_.currents.photo_emission_a_per_m2);
            history_thermionic_current_.push_back(status_.currents.thermionic_emission_a_per_m2);
            history_field_emission_current_.push_back(status_.currents.field_emission_a_per_m2);
            history_leakage_current_.push_back(leakage_current);
            history_ram_current_.push_back(status_.currents.ram_ion_current_a_per_m2);
            history_body_potential_.push_back(status_.body_potential_v);
            history_patch_potential_.push_back(status_.patch_potential_v);
            history_circuit_branch_current_.push_back(status_.circuit_branch_current_a);
            history_current_derivative_.push_back(
                status_.currents.current_derivative_a_per_m2_per_v);
            history_pic_recalibration_marker_.push_back(recalibrated ? 1.0 : 0.0);
            const PicMccCurrentSample latest_sample =
                current_model_ ? current_model_->latestLivePicSample() : PicMccCurrentSample{};
            history_live_pic_electron_current_.push_back(
                latest_sample.electron_collection_current_density_a_per_m2);
            history_live_pic_ion_current_.push_back(
                latest_sample.ion_collection_current_density_a_per_m2);
            history_live_pic_net_current_.push_back(
                latest_sample.net_collection_current_density_a_per_m2);
            history_live_pic_derivative_.push_back(
                latest_sample.current_derivative_a_per_m2_per_v);
            history_live_pic_collision_count_.push_back(
                static_cast<double>(latest_sample.total_collisions));
            history_live_pic_mcc_enabled_.push_back(latest_sample.mcc_enabled ? 1.0 : 0.0);
            history_net_current_.push_back(net_current);
            history_capacitance_.push_back(status_.state.capacitance_per_area_f_per_m2);
            history_effective_conductivity_.push_back(effective_conductivity);
            history_effective_sheath_length_.push_back(effective_sheath_length);
            history_normal_electric_field_.push_back(
                final_runtime_state.normal_electric_field_v_per_m);
            history_local_charge_density_.push_back(final_runtime_state.local_charge_density_c_per_m3);
            history_adaptive_time_step_.push_back(dt);
            history_internal_substeps_.push_back(static_cast<double>(substeps));
            history_electron_calibration_factor_.push_back(
                current_model_ ? current_model_->electronCalibrationFactor() : 1.0);
            history_ion_calibration_factor_.push_back(
                current_model_ ? current_model_->ionCalibrationFactor() : 1.0);
            history_equilibrium_potential_.push_back(target_potential);
            history_equilibrium_error_.push_back(status_.equilibrium_error);
            history_shared_surface_pre_global_solve_patch_potential_spread_.push_back(
                latest_pre_global_solve_patch_spread_v);
            history_shared_surface_patch_potential_spread_reduction_v_.push_back(
                latest_patch_spread_reduction_v);
            history_shared_surface_patch_potential_spread_reduction_ratio_.push_back(
                latest_patch_spread_reduction_ratio);
            history_shared_surface_current_matrix_coupling_active_.push_back(
                latest_shared_current_matrix_coupling_offdiag_entries > 0.5 ? 1.0 : 0.0);
            history_shared_surface_current_matrix_coupling_offdiag_entries_.push_back(
                latest_shared_current_matrix_coupling_offdiag_entries);
            history_shared_surface_global_coupled_solve_active_.push_back(
                useSharedSurfaceGlobalCoupledSolve() ? 1.0 : 0.0);
            history_shared_surface_global_coupled_solve_iterations_.push_back(
                latest_shared_global_coupled_solve_iterations);
            history_shared_surface_global_coupled_solve_converged_.push_back(
                latest_shared_global_coupled_solve_converged);
            history_shared_surface_global_coupled_solve_max_delta_v_.push_back(
                latest_shared_global_coupled_solve_max_delta_v);
            history_shared_surface_live_pic_coupled_refresh_active_.push_back(
                useSharedSurfaceLivePicCoupledRefresh() ? 1.0 : 0.0);
            history_shared_surface_live_pic_coupled_refresh_count_.push_back(
                latest_shared_live_pic_coupled_refresh_count);
            history_shared_surface_particle_transport_coupling_active_.push_back(
                latest_shared_particle_transport_offdiag_entries > 0.5 ? 1.0 : 0.0);
            history_shared_surface_particle_transport_offdiag_entries_.push_back(
                latest_shared_particle_transport_offdiag_entries);
            history_shared_surface_particle_transport_total_conductance_s_.push_back(
                latest_shared_particle_transport_total_conductance_s);
            history_shared_surface_particle_transport_conservation_error_a_per_v_.push_back(
                latest_shared_particle_transport_conservation_error_a_per_v);
            history_shared_surface_particle_transport_charge_c_.push_back(
                shared_particle_transport_charge_c_);
            history_shared_surface_particle_transport_reference_shift_v_.push_back(
                shared_particle_transport_reference_shift_v_);
            appendSurfaceNodeDiagnostics(latest_branch_currents_a, effective_sheath_length);
            return true;
        };
        return commit_reference_stage();
    }

    const double reference_capacitance =
        std::max(1.0e-12, status_.state.capacitance_per_area_f_per_m2);
    const double reference_derivative =
        estimateCurrentDerivative(status_.state.surface_potential_v);
    Coupling::SurfaceAdaptiveSubstepInputs adaptive_substep_inputs;
    adaptive_substep_inputs.base_substeps = transition_runtime_internal_substeps_;
    adaptive_substep_inputs.reference_potential_v = status_.state.surface_potential_v;
    adaptive_substep_inputs.dt_s = dt;
    adaptive_substep_inputs.current_derivative_a_per_m2_per_v = reference_derivative;
    adaptive_substep_inputs.capacitance_per_area_f_per_m2 = reference_capacitance;
    const std::size_t substeps =
        Coupling::computeSurfaceAdaptiveSubstepCount(adaptive_substep_inputs);
    const double sub_dt = dt / static_cast<double>(substeps);
    auto state = status_.state;
    SurfaceCurrents currents;
    double leakage_current = 0.0;
    double net_current = 0.0;
    double effective_conductivity =
        computeEffectiveConductivity(state.surface_potential_v,
                                     resolvePatchMaterial(
                                         config_, circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                                         circuit_model_
                                             ? circuit_model_->nodeName(
                                                   circuit_model_->primaryPatchNodeIndex())
                                             : std::string{"patch"}));
    const double effective_sheath_length = computeEffectiveSheathLength();
    const bool use_reference_explicit =
        config_.use_reference_current_balance &&
        config_.regime == SurfaceChargingRegime::GeoKineticPicLike &&
        !config_.enable_body_patch_circuit && !config_.enable_pic_calibration &&
        !config_.enable_live_pic_window;

    for (std::size_t step = 0; step < substeps; ++step)
    {
        state.capacitance_per_area_f_per_m2 =
            computeCapacitancePerArea(buildRuntimeState(
                status_.body_potential_v, state.surface_potential_v, config_.surface_area_m2,
                effective_conductivity, effective_sheath_length,
                circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                               : std::string{"patch"}));
        if (use_reference_explicit)
        {
            currents = computeSurfaceCurrents(state.surface_potential_v);
            leakage_current = computeLeakageCurrentDensity(state.surface_potential_v);
            net_current = currents.total_current_a_per_m2 + leakage_current;
            effective_conductivity =
                computeEffectiveConductivity(state.surface_potential_v,
                                             resolvePatchMaterial(
                                                 config_,
                                                 circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                                                 circuit_model_
                                                     ? circuit_model_->nodeName(
                                                           circuit_model_->primaryPatchNodeIndex())
                                                     : std::string{"patch"}));
            state = accumulation_model_.advance(state, net_current, sub_dt);
            state.surface_potential_v =
                clampSigned(state.surface_potential_v, config_.max_abs_potential_v);
            state.surface_charge_density_c_per_m2 =
                state.surface_potential_v * state.capacitance_per_area_f_per_m2;
            state.stored_energy_j_per_m2 =
                0.5 * state.capacitance_per_area_f_per_m2 * state.surface_potential_v *
                state.surface_potential_v;
        }
        else
        {
            const double updated_potential =
                advancePotentialImplicit(state.surface_potential_v,
                                         state.capacitance_per_area_f_per_m2, sub_dt);
            currents = computeSurfaceCurrents(updated_potential);
            leakage_current = computeLeakageCurrentDensity(updated_potential);
            net_current = currents.total_current_a_per_m2 + leakage_current;
            effective_conductivity =
                computeEffectiveConductivity(updated_potential,
                                             resolvePatchMaterial(
                                                 config_,
                                                 circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                                                 circuit_model_
                                                     ? circuit_model_->nodeName(
                                                           circuit_model_->primaryPatchNodeIndex())
                                                     : std::string{"patch"}));
            state = accumulation_model_.advance(state, net_current, sub_dt);
            state.surface_potential_v = updated_potential;
            state.surface_charge_density_c_per_m2 =
                state.surface_potential_v * state.capacitance_per_area_f_per_m2;
            state.stored_energy_j_per_m2 =
                0.5 * state.capacitance_per_area_f_per_m2 * state.surface_potential_v *
                state.surface_potential_v;
        }
        if (!std::isfinite(state.surface_potential_v) || !std::isfinite(state.surface_charge_density_c_per_m2))
        {
            return false;
        }
    }

    if (use_reference_explicit)
    {
        currents = computeSurfaceCurrents(state.surface_potential_v);
        leakage_current = computeLeakageCurrentDensity(state.surface_potential_v);
        net_current = currents.total_current_a_per_m2 + leakage_current;
        effective_conductivity =
            computeEffectiveConductivity(state.surface_potential_v,
                                         resolvePatchMaterial(
                                             config_,
                                             circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                                             circuit_model_
                                                 ? circuit_model_->nodeName(
                                                       circuit_model_->primaryPatchNodeIndex())
                                                 : std::string{"patch"}));
    }

    const auto commit_single_node_stage = [&]() -> bool {
        const double target_potential =
            config_.floating ? computeFloatingPotential() : state.surface_potential_v;
        status_.currents = currents;
        status_.state = state;
        status_.body_potential_v = config_.body_initial_potential_v;
        status_.patch_potential_v = state.surface_potential_v;
        status_.circuit_branch_current_a = 0.0;
        status_.pic_recalibrated = false;

        status_.equilibrium_error = status_.state.surface_potential_v - target_potential;
        status_.equilibrium_reached = std::abs(status_.equilibrium_error) < 1.0e-2;
        status_.time_s += dt;
        status_.steps_completed += 1;
        const auto final_runtime_state =
            buildRuntimeState(status_.body_potential_v, status_.patch_potential_v, config_.surface_area_m2,
                              effective_conductivity, effective_sheath_length,
                              circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
                              circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                                             : std::string{"patch"});
        if (circuit_node_potentials_.size() >= 2)
        {
            circuit_node_potentials_[0] = status_.body_potential_v;
            circuit_node_potentials_[1] = status_.patch_potential_v;
        }

        history_time_.push_back(status_.time_s);
        history_potential_.push_back(status_.state.surface_potential_v);
        history_charge_.push_back(status_.state.surface_charge_density_c_per_m2);
        history_current_.push_back(status_.currents.total_current_a_per_m2);
        history_electron_current_.push_back(status_.currents.electron_current_a_per_m2);
        history_ion_current_.push_back(status_.currents.ion_current_a_per_m2);
        history_secondary_current_.push_back(status_.currents.secondary_emission_a_per_m2);
        history_ion_secondary_current_.push_back(status_.currents.ion_secondary_emission_a_per_m2);
        history_backscatter_current_.push_back(status_.currents.backscatter_emission_a_per_m2);
        history_photo_current_.push_back(status_.currents.photo_emission_a_per_m2);
        history_thermionic_current_.push_back(status_.currents.thermionic_emission_a_per_m2);
        history_field_emission_current_.push_back(status_.currents.field_emission_a_per_m2);
        history_leakage_current_.push_back(leakage_current);
        history_ram_current_.push_back(status_.currents.ram_ion_current_a_per_m2);
        history_body_potential_.push_back(status_.body_potential_v);
        history_patch_potential_.push_back(status_.state.surface_potential_v);
        history_circuit_branch_current_.push_back(status_.circuit_branch_current_a);
        history_current_derivative_.push_back(status_.currents.current_derivative_a_per_m2_per_v);
        history_pic_recalibration_marker_.push_back(status_.pic_recalibrated ? 1.0 : 0.0);
        const PicMccCurrentSample latest_sample =
            current_model_ ? current_model_->latestLivePicSample() : PicMccCurrentSample{};
        history_live_pic_electron_current_.push_back(
            latest_sample.electron_collection_current_density_a_per_m2);
        history_live_pic_ion_current_.push_back(
            latest_sample.ion_collection_current_density_a_per_m2);
        history_live_pic_net_current_.push_back(
            latest_sample.net_collection_current_density_a_per_m2);
        history_live_pic_derivative_.push_back(
            latest_sample.current_derivative_a_per_m2_per_v);
        history_live_pic_collision_count_.push_back(
            static_cast<double>(latest_sample.total_collisions));
        history_live_pic_mcc_enabled_.push_back(latest_sample.mcc_enabled ? 1.0 : 0.0);
        history_net_current_.push_back(net_current);
        history_capacitance_.push_back(status_.state.capacitance_per_area_f_per_m2);
        history_effective_conductivity_.push_back(effective_conductivity);
        history_effective_sheath_length_.push_back(effective_sheath_length);
        history_normal_electric_field_.push_back(final_runtime_state.normal_electric_field_v_per_m);
        history_local_charge_density_.push_back(final_runtime_state.local_charge_density_c_per_m3);
        history_adaptive_time_step_.push_back(dt);
        history_internal_substeps_.push_back(static_cast<double>(substeps));
        history_electron_calibration_factor_.push_back(
            current_model_ ? current_model_->electronCalibrationFactor() : 1.0);
        history_ion_calibration_factor_.push_back(
            current_model_ ? current_model_->ionCalibrationFactor() : 1.0);
        history_equilibrium_potential_.push_back(target_potential);
        history_equilibrium_error_.push_back(status_.equilibrium_error);
        history_shared_surface_pre_global_solve_patch_potential_spread_.push_back(0.0);
        history_shared_surface_patch_potential_spread_reduction_v_.push_back(0.0);
        history_shared_surface_patch_potential_spread_reduction_ratio_.push_back(0.0);
        history_shared_surface_current_matrix_coupling_active_.push_back(0.0);
        history_shared_surface_current_matrix_coupling_offdiag_entries_.push_back(0.0);
        history_shared_surface_global_coupled_solve_active_.push_back(0.0);
        history_shared_surface_global_coupled_solve_iterations_.push_back(0.0);
        history_shared_surface_global_coupled_solve_converged_.push_back(0.0);
        history_shared_surface_global_coupled_solve_max_delta_v_.push_back(0.0);
        history_shared_surface_live_pic_coupled_refresh_active_.push_back(0.0);
        history_shared_surface_live_pic_coupled_refresh_count_.push_back(0.0);
        history_shared_surface_particle_transport_coupling_active_.push_back(0.0);
        history_shared_surface_particle_transport_offdiag_entries_.push_back(0.0);
        history_shared_surface_particle_transport_total_conductance_s_.push_back(0.0);
        history_shared_surface_particle_transport_conservation_error_a_per_v_.push_back(0.0);
        history_shared_surface_particle_transport_charge_c_.push_back(
            shared_particle_transport_charge_c_);
        history_shared_surface_particle_transport_reference_shift_v_.push_back(
            shared_particle_transport_reference_shift_v_);
        std::vector<double> branch_currents_a(
            circuit_model_ ? circuit_model_->branchCount() : 1, 0.0);
        appendSurfaceNodeDiagnostics(branch_currents_a, effective_sheath_length);
        return true;
    };
    return commit_single_node_stage();
}

void DensePlasmaSurfaceCharging::reset()
{
    *this = DensePlasmaSurfaceCharging{};
}

bool DensePlasmaSurfaceCharging::exportResults(const std::filesystem::path& csv_path) const
{
    Output::ColumnarDataSet data_set;
    LegacyBenchmarkCaseDefinition benchmark_case_definition;
    LegacyBenchmarkMetrics benchmark_patch_metrics;
    LegacyBenchmarkMetrics benchmark_body_metrics;
    LegacyBenchmarkBodyTailMetrics benchmark_body_tail_metrics;
    LegacyBenchmarkAcceptanceCriteria benchmark_acceptance;
    LegacyBenchmarkAcceptanceGateResult benchmark_acceptance_gate;
    bool have_benchmark_definition = false;
    auto simulation_artifact_path = csv_path;
    simulation_artifact_path.replace_extension(".simulation_artifact.json");
    auto benchmark_case_path = csv_path;
    benchmark_case_path.replace_extension(".benchmark_case.json");
    const auto constant_series = [this](double value) {
        return std::vector<double>(history_time_.size(), value);
    };
    const auto scaled_series = [](const std::vector<double>& values, double scale) {
        std::vector<double> scaled;
        scaled.reserve(values.size());
        for (const double value : values)
        {
            scaled.push_back(value * scale);
        }
        return scaled;
    };
    const auto cycle_series = [this]() {
        std::vector<double> values(history_time_.size(), 0.0);
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            values[i] = static_cast<double>(i);
        }
        return values;
    };
    const auto first_equilibrium_time_ms = [this]() {
        for (std::size_t i = 0; i < history_time_.size() && i < history_equilibrium_error_.size(); ++i)
        {
            if (std::abs(history_equilibrium_error_[i]) < 1.0e-2)
            {
                return history_time_[i] * 1.0e3;
            }
        }
        return history_time_.empty() ? 0.0 : history_time_.back() * 1.0e3;
    };
    const auto resolved_internal_substeps_series = [this]() {
        if (history_internal_substeps_.size() == history_time_.size())
        {
            return history_internal_substeps_;
        }
        return std::vector<double>(
            history_time_.size(),
            static_cast<double>(transition_runtime_internal_substeps_));
    };
    data_set.axis_name = "time_s";
    data_set.axis_values = history_time_;
    data_set.scalar_series["surface_potential_v"] = history_potential_;
    data_set.scalar_series["surface_charge_density_c_per_m2"] = history_charge_;
    data_set.scalar_series["total_current_density_a_per_m2"] = history_current_;
    data_set.scalar_series["electron_current_density_a_per_m2"] = history_electron_current_;
    data_set.scalar_series["ion_current_density_a_per_m2"] = history_ion_current_;
    data_set.scalar_series["secondary_emission_density_a_per_m2"] = history_secondary_current_;
    data_set.scalar_series["ion_secondary_emission_density_a_per_m2"] =
        history_ion_secondary_current_;
    data_set.scalar_series["backscatter_emission_density_a_per_m2"] =
        history_backscatter_current_;
    data_set.scalar_series["photo_emission_density_a_per_m2"] = history_photo_current_;
    data_set.scalar_series["thermionic_emission_density_a_per_m2"] = history_thermionic_current_;
    data_set.scalar_series["field_emission_density_a_per_m2"] = history_field_emission_current_;
    data_set.scalar_series["leakage_current_density_a_per_m2"] = history_leakage_current_;
    data_set.scalar_series["ram_ion_current_density_a_per_m2"] = history_ram_current_;
    data_set.scalar_series["body_potential_v"] = history_body_potential_;
    data_set.scalar_series["patch_potential_v"] = history_patch_potential_;
    data_set.scalar_series["circuit_branch_current_a"] = history_circuit_branch_current_;
    data_set.scalar_series["current_derivative_a_per_m2_per_v"] = history_current_derivative_;
    data_set.scalar_series["pic_recalibration_trigger"] = history_pic_recalibration_marker_;
    data_set.scalar_series["live_pic_electron_collection_density_a_per_m2"] =
        history_live_pic_electron_current_;
    data_set.scalar_series["live_pic_ion_collection_density_a_per_m2"] =
        history_live_pic_ion_current_;
    data_set.scalar_series["live_pic_net_collection_density_a_per_m2"] =
        history_live_pic_net_current_;
    data_set.scalar_series["live_pic_collection_derivative_a_per_m2_per_v"] =
        history_live_pic_derivative_;
    data_set.scalar_series["live_pic_collision_count"] = history_live_pic_collision_count_;
    data_set.scalar_series["live_pic_mcc_enabled"] = history_live_pic_mcc_enabled_;
    data_set.scalar_series["net_current_density_a_per_m2"] = history_net_current_;
    data_set.scalar_series["capacitance_per_area_f_per_m2"] = history_capacitance_;
    data_set.scalar_series["effective_conductivity_s_per_m"] = history_effective_conductivity_;
    data_set.scalar_series["effective_sheath_length_m"] = history_effective_sheath_length_;
    data_set.scalar_series["shared_surface_patch_potential_v"] =
        history_shared_surface_patch_potential_;
    data_set.scalar_series["shared_surface_patch_area_m2"] = history_shared_surface_patch_area_;
    data_set.scalar_series["shared_surface_reference_potential_v"] =
        history_shared_surface_reference_potential_;
    data_set.scalar_series["shared_surface_effective_sheath_length_m"] =
        history_shared_surface_effective_sheath_length_;
    data_set.scalar_series["shared_surface_sheath_charge_c"] =
        history_shared_surface_sheath_charge_;
    data_set.scalar_series["shared_surface_runtime_enabled"] =
        history_shared_surface_runtime_enabled_;
    data_set.scalar_series["shared_surface_pre_global_solve_patch_potential_spread_v"] =
        history_shared_surface_pre_global_solve_patch_potential_spread_;
    data_set.scalar_series["shared_surface_patch_potential_spread_reduction_v"] =
        history_shared_surface_patch_potential_spread_reduction_v_;
    data_set.scalar_series["shared_surface_patch_potential_spread_reduction_ratio"] =
        history_shared_surface_patch_potential_spread_reduction_ratio_;
    data_set.scalar_series["shared_surface_current_matrix_coupling_active"] =
        history_shared_surface_current_matrix_coupling_active_;
    data_set.scalar_series["shared_surface_current_matrix_coupling_offdiag_entries"] =
        history_shared_surface_current_matrix_coupling_offdiag_entries_;
    data_set.scalar_series["shared_surface_global_coupled_solve_active"] =
        history_shared_surface_global_coupled_solve_active_;
    data_set.scalar_series["shared_surface_global_coupled_solve_iterations"] =
        history_shared_surface_global_coupled_solve_iterations_;
    data_set.scalar_series["shared_surface_global_coupled_solve_converged"] =
        history_shared_surface_global_coupled_solve_converged_;
    data_set.scalar_series["shared_surface_global_coupled_solve_max_delta_v"] =
        history_shared_surface_global_coupled_solve_max_delta_v_;
    data_set.scalar_series["shared_surface_live_pic_coupled_refresh_active"] =
        history_shared_surface_live_pic_coupled_refresh_active_;
    data_set.scalar_series["shared_surface_live_pic_coupled_refresh_count"] =
        history_shared_surface_live_pic_coupled_refresh_count_;
    data_set.scalar_series["shared_surface_particle_transport_coupling_active"] =
        history_shared_surface_particle_transport_coupling_active_;
    data_set.scalar_series["shared_surface_particle_transport_offdiag_entries"] =
        history_shared_surface_particle_transport_offdiag_entries_;
    data_set.scalar_series["shared_surface_particle_transport_total_conductance_s"] =
        history_shared_surface_particle_transport_total_conductance_s_;
    data_set.scalar_series["shared_surface_particle_transport_conservation_error_a_per_v"] =
        history_shared_surface_particle_transport_conservation_error_a_per_v_;
    data_set.scalar_series["shared_surface_particle_transport_charge_c"] =
        history_shared_surface_particle_transport_charge_c_;
    data_set.scalar_series["shared_surface_particle_transport_reference_shift_v"] =
        history_shared_surface_particle_transport_reference_shift_v_;
    data_set.scalar_series["normal_electric_field_v_per_m"] = history_normal_electric_field_;
    data_set.scalar_series["local_charge_density_c_per_m3"] = history_local_charge_density_;
    data_set.scalar_series["adaptive_time_step_s"] = history_adaptive_time_step_;
    data_set.scalar_series["resolved_internal_substeps"] =
        resolved_internal_substeps_series();
    data_set.scalar_series["electron_pic_calibration_factor"] = history_electron_calibration_factor_;
    data_set.scalar_series["ion_pic_calibration_factor"] = history_ion_calibration_factor_;
    data_set.scalar_series["floating_equilibrium_potential_v"] = history_equilibrium_potential_;
    data_set.scalar_series["equilibrium_error_v"] = history_equilibrium_error_;
    data_set.scalar_series["current_algorithm_mode_id"] =
        constant_series(config_.current_algorithm_mode ==
                                SurfaceCurrentAlgorithmMode::LegacyRefCompatible
                            ? 1.0
                            : 0.0);
    data_set.scalar_series["benchmark_mode_id"] =
        constant_series(config_.benchmark_mode == SurfaceBenchmarkMode::LegacyRefCompatible ? 1.0
                                                                                            : 0.0);
    data_set.scalar_series["runtime_route_id"] = constant_series(runtimeRouteId(config_.runtime_route));
    data_set.scalar_series["surface_pic_strategy_id"] =
        constant_series(surfacePicStrategyId(config_.surface_pic_strategy));
    data_set.scalar_series["surface_legacy_input_adapter_id"] =
        constant_series(legacyInputAdapterId(config_.legacy_input_adapter_kind));
    data_set.scalar_series["surface_pic_runtime_id"] =
        constant_series(surfacePicRuntimeId(config_.surface_pic_runtime_kind));
    data_set.scalar_series["surface_instrument_set_id"] =
        constant_series(surfaceInstrumentSetId(config_.surface_instrument_set_kind));
    data_set.scalar_series["benchmark_source_id"] =
        constant_series(benchmarkSourceId(config_.benchmark_source));
    data_set.scalar_series["benchmark_execution_mode_id"] =
        constant_series(legacyBenchmarkExecutionModeId(config_.legacy_benchmark_execution_mode));
    data_set.scalar_series["surface_circuit_node_count"] =
        constant_series(circuit_model_ ? static_cast<double>(circuit_model_->nodeCount()) : 0.0);
    if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark)
    {
        data_set.scalar_series["Cycle"] = cycle_series();
        data_set.scalar_series["Vs"] = history_potential_;
        data_set.scalar_series["Time_ms"] = scaled_series(history_time_, 1.0e3);
        data_set.scalar_series["Jnet"] = scaled_series(history_current_, 1.0e9);
        data_set.scalar_series["Je"] = scaled_series(history_electron_current_, 1.0e9);
        data_set.scalar_series["Jse"] = scaled_series(history_secondary_current_, 1.0e9);
        data_set.scalar_series["Jb"] = scaled_series(history_backscatter_current_, 1.0e9);
        data_set.scalar_series["Ji"] = scaled_series(history_ion_current_, 1.0e9);
        data_set.scalar_series["Jsi"] = scaled_series(history_ion_secondary_current_, 1.0e9);
        data_set.scalar_series["Jph"] = scaled_series(history_photo_current_, 1.0e9);
        data_set.scalar_series["Jcond"] = scaled_series(history_leakage_current_, 1.0e9);
        data_set.scalar_series["benchmark_cycle_index"] = data_set.scalar_series["Cycle"];
        data_set.scalar_series["benchmark_time_ms"] = data_set.scalar_series["Time_ms"];
        data_set.scalar_series["benchmark_vs_v"] = data_set.scalar_series["Vs"];
        data_set.scalar_series["benchmark_jnet_na_per_m2"] = data_set.scalar_series["Jnet"];
        data_set.scalar_series["benchmark_je_na_per_m2"] = data_set.scalar_series["Je"];
        data_set.scalar_series["benchmark_jse_na_per_m2"] = data_set.scalar_series["Jse"];
        data_set.scalar_series["benchmark_jb_na_per_m2"] = data_set.scalar_series["Jb"];
        data_set.scalar_series["benchmark_ji_na_per_m2"] = data_set.scalar_series["Ji"];
        data_set.scalar_series["benchmark_jsi_na_per_m2"] = data_set.scalar_series["Jsi"];
        data_set.scalar_series["benchmark_jph_na_per_m2"] = data_set.scalar_series["Jph"];
        data_set.scalar_series["benchmark_jcond_na_per_m2"] = data_set.scalar_series["Jcond"];
    }
    for (std::size_t node_index = 0; node_index < history_surface_node_potentials_.size(); ++node_index)
    {
        if (history_surface_node_potentials_[node_index].size() == history_time_.size())
        {
            data_set.scalar_series["surface_node_" + std::to_string(node_index) + "_potential_v"] =
                history_surface_node_potentials_[node_index];
        }
        const std::string prefix = "surface_node_" + std::to_string(node_index) + "_";
        const auto export_node_series = [&](const std::string& suffix,
                                            const std::vector<std::vector<double>>& history) {
            if (node_index < history.size() && history[node_index].size() == history_time_.size())
            {
                data_set.scalar_series[prefix + suffix] = history[node_index];
            }
        };
        export_node_series("total_current_density_a_per_m2", history_surface_node_total_currents_);
        export_node_series("electron_current_density_a_per_m2",
                           history_surface_node_electron_currents_);
        export_node_series("ion_current_density_a_per_m2", history_surface_node_ion_currents_);
        export_node_series("secondary_emission_density_a_per_m2",
                           history_surface_node_secondary_currents_);
        export_node_series("ion_secondary_emission_density_a_per_m2",
                           history_surface_node_ion_secondary_currents_);
        export_node_series("backscatter_emission_density_a_per_m2",
                           history_surface_node_backscatter_currents_);
        export_node_series("photo_emission_density_a_per_m2", history_surface_node_photo_currents_);
        export_node_series("thermionic_emission_density_a_per_m2",
                           history_surface_node_thermionic_currents_);
        export_node_series("field_emission_density_a_per_m2",
                           history_surface_node_field_emission_currents_);
        export_node_series("conduction_current_density_a_per_m2",
                           history_surface_node_conduction_currents_);
        export_node_series("ram_ion_current_density_a_per_m2", history_surface_node_ram_currents_);
        export_node_series("current_derivative_a_per_m2_per_v",
                           history_surface_node_current_derivatives_);
        export_node_series("effective_sheath_length_m",
                           history_surface_node_effective_sheath_lengths_);
        export_node_series("normal_electric_field_v_per_m",
                           history_surface_node_normal_electric_fields_);
        export_node_series("local_charge_density_c_per_m3",
                           history_surface_node_local_charge_densities_);
        export_node_series("propagated_reference_potential_v",
                           history_surface_node_propagated_reference_potentials_);
        export_node_series("field_solver_reference_potential_v",
                           history_surface_node_field_solver_reference_potentials_);
        export_node_series("shared_runtime_enabled",
                           history_surface_node_shared_runtime_enabled_);
        export_node_series("shared_patch_potential_v",
                           history_surface_node_shared_patch_potentials_);
        export_node_series("shared_patch_area_m2", history_surface_node_shared_patch_areas_);
        export_node_series("shared_reference_potential_v",
                           history_surface_node_shared_reference_potentials_);
        export_node_series("shared_sheath_charge_c",
                           history_surface_node_shared_sheath_charges_);
        export_node_series("distributed_particle_transport_charge_c",
                           history_surface_node_distributed_particle_transport_charges_);
        export_node_series("distributed_particle_transport_reference_shift_v",
                           history_surface_node_distributed_particle_transport_reference_shifts_);
        export_node_series("distributed_particle_transport_net_flux_a",
                           history_surface_node_distributed_particle_transport_net_fluxes_);
        export_node_series("graph_capacitance_diagonal_f",
                           history_surface_node_graph_capacitance_diagonals_);
        export_node_series("graph_capacitance_row_sum_f",
                           history_surface_node_graph_capacitance_row_sums_);
        export_node_series("field_solver_coupling_gain",
                           history_surface_node_field_solver_coupling_gains_);
        export_node_series("field_solver_capacitance_scale",
                           history_surface_node_field_solver_capacitance_scales_);
        export_node_series("pseudo_volume_m3", history_surface_node_pseudo_volumes_);
        export_node_series("volume_projection_weight_sum",
                           history_surface_node_volume_projection_weight_sums_);
        export_node_series("volume_mesh_coupling_gain",
                           history_surface_node_volume_mesh_coupling_gains_);
        export_node_series("volume_potential_v", history_surface_node_volume_potentials_);
        export_node_series("deposited_charge_c", history_surface_node_deposited_charges_);
        export_node_series("poisson_residual_v_m", history_surface_node_poisson_residuals_);
        export_node_series("volume_solver_mode_id", history_surface_node_volume_solver_mode_ids_);
        export_node_series("volume_solver_iterations",
                           history_surface_node_volume_solver_iterations_);
        export_node_series("volume_solver_linear_iterations",
                           history_surface_node_volume_solver_linear_iterations_);
        export_node_series("volume_solver_converged",
                           history_surface_node_volume_solver_converged_);
        export_node_series("volume_solver_residual_norm",
                           history_surface_node_volume_solver_residual_norms_);
        export_node_series("volume_solver_max_delta_v",
                           history_surface_node_volume_solver_max_deltas_);
        export_node_series("volume_solver_matrix_nnz",
                           history_surface_node_volume_solver_matrix_nnzs_);
        export_node_series("volume_solver_cell_count",
                           history_surface_node_volume_solver_cell_counts_);
        export_node_series("field_volume_coupling_iterations",
                           history_surface_node_field_volume_coupling_iterations_);
        export_node_series("field_volume_coupling_converged",
                           history_surface_node_field_volume_coupling_converged_);
        export_node_series("field_volume_coupling_max_delta",
                           history_surface_node_field_volume_coupling_max_deltas_);
        export_node_series("field_volume_coupling_relaxation_used",
                           history_surface_node_field_volume_coupling_relaxation_used_);
        export_node_series("external_volume_feedback_blend_factor",
                           history_surface_node_external_volume_feedback_blend_factors_);
        export_node_series("external_volume_feedback_mismatch_metric",
                           history_surface_node_external_volume_feedback_mismatch_metrics_);
        export_node_series("external_volume_feedback_applied",
                           history_surface_node_external_volume_feedback_applied_);
        export_node_series("electron_pic_calibration_factor",
                           history_surface_node_electron_calibration_factors_);
        export_node_series("ion_pic_calibration_factor",
                           history_surface_node_ion_calibration_factors_);
        export_node_series("live_pic_electron_collection_density_a_per_m2",
                           history_surface_node_live_pic_electron_currents_);
        export_node_series("live_pic_ion_collection_density_a_per_m2",
                           history_surface_node_live_pic_ion_currents_);
        export_node_series("live_pic_net_collection_density_a_per_m2",
                           history_surface_node_live_pic_net_currents_);
        export_node_series("live_pic_collection_derivative_a_per_m2_per_v",
                           history_surface_node_live_pic_derivatives_);
        export_node_series("live_pic_collision_count",
                           history_surface_node_live_pic_collision_counts_);
        export_node_series("live_pic_mcc_enabled",
                           history_surface_node_live_pic_mcc_enabled_);
        data_set.scalar_series[prefix + "surface_pic_enabled"] = constant_series(
            (circuit_model_ && circuit_model_->nodeIsPatch(node_index) &&
             (config_.runtime_route == SurfaceRuntimeRoute::SurfacePic ||
              config_.runtime_route == SurfaceRuntimeRoute::SurfacePicHybrid))
                ? 1.0
                : 0.0);
    }
    for (std::size_t branch_index = 0; branch_index < history_surface_branch_currents_.size();
         ++branch_index)
    {
        if (history_surface_branch_currents_[branch_index].size() == history_time_.size())
        {
            data_set.scalar_series["surface_interface_" + std::to_string(branch_index) + "_current_a"] =
                history_surface_branch_currents_[branch_index];
        }
        if (branch_index < history_surface_branch_conductances_.size() &&
            history_surface_branch_conductances_[branch_index].size() == history_time_.size())
        {
            data_set.scalar_series["surface_interface_" + std::to_string(branch_index) +
                                   "_conductance_s"] =
                history_surface_branch_conductances_[branch_index];
        }
        if (branch_index < history_surface_branch_voltage_drops_.size() &&
            history_surface_branch_voltage_drops_[branch_index].size() == history_time_.size())
        {
            data_set.scalar_series["surface_interface_" + std::to_string(branch_index) +
                                   "_voltage_drop_v"] =
                history_surface_branch_voltage_drops_[branch_index];
        }
        if (branch_index < history_surface_branch_power_w_.size() &&
            history_surface_branch_power_w_[branch_index].size() == history_time_.size())
        {
            data_set.scalar_series["surface_interface_" + std::to_string(branch_index) +
                                   "_power_w"] =
                history_surface_branch_power_w_[branch_index];
        }
        if (branch_index < history_surface_branch_mutual_capacitances_.size() &&
            history_surface_branch_mutual_capacitances_[branch_index].size() == history_time_.size())
        {
            data_set.scalar_series["surface_interface_" + std::to_string(branch_index) +
                                   "_mutual_capacitance_f"] =
                history_surface_branch_mutual_capacitances_[branch_index];
        }
    }
    data_set.metadata["module"] = "Surface Charging";
    data_set.metadata["current_algorithm_mode"] =
        current_model_ ? current_model_->algorithmName() : "LegacySurfaceCurrents";
    data_set.metadata["benchmark_mode"] =
        config_.benchmark_mode == SurfaceBenchmarkMode::LegacyRefCompatible
            ? "LegacyRefCompatible"
            : "UnifiedKernelAligned";
    data_set.metadata["runtime_route"] = runtimeRouteName(config_.runtime_route);
    data_set.metadata["surface_pic_strategy"] =
        surfacePicStrategyName(config_.surface_pic_strategy);
    data_set.metadata["surface_legacy_input_adapter"] =
        legacyInputAdapterName(config_.legacy_input_adapter_kind);
    data_set.metadata["surface_pic_runtime"] =
        surfacePicRuntimeName(config_.surface_pic_runtime_kind);
    data_set.metadata["surface_instrument_set"] =
        surfaceInstrumentSetName(config_.surface_instrument_set_kind);
    data_set.metadata["surface_reference_family"] = config_.reference_family;
    data_set.metadata["surface_reference_case_id"] = config_.reference_case_id;
    data_set.metadata["surface_reference_matrix_case_id"] =
        config_.reference_matrix_case_id;
    data_set.metadata["surface_pic_runtime_contract_id"] = "surface-pic-runtime-v1";
    data_set.metadata["surface_instrument_contract_id"] = "surface-instrument-observer-v1";
    data_set.metadata["surface_instrument_contract_version"] = "v1";
    data_set.metadata["surface_pic_runtime_supports_shared_surface_nodes"] =
        (config_.surface_pic_runtime_kind ==
         SurfacePicRuntimeKind::GraphCoupledSharedSurface)
            ? "1"
            : "0";
    data_set.metadata["surface_pic_runtime_shared_patch_potential_v"] =
        std::to_string(computeSharedSurfacePatchPotentialV());
    data_set.metadata["surface_pic_runtime_shared_patch_area_m2"] =
        std::to_string(computeSharedSurfacePatchAreaM2());
    data_set.metadata["surface_pic_runtime_shared_patch_spread_v"] =
        std::to_string(computeSharedSurfacePatchPotentialSpreadV());
    data_set.metadata["surface_pic_runtime_shared_effective_sheath_length_m"] =
        std::to_string(computeSharedSurfaceEffectiveSheathLengthM(computeEffectiveSheathLength()));
    data_set.metadata["surface_pic_runtime_shared_sheath_charge_c"] =
        history_shared_surface_sheath_charge_.empty()
            ? "0"
            : std::to_string(history_shared_surface_sheath_charge_.back());
    data_set.metadata["surface_pic_runtime_shared_sheath_contract_id"] =
        "surface-pic-shared-sheath-v1";
    data_set.metadata["surface_pic_runtime_boundary_observer_contract_id"] =
        "surface-pic-boundary-observer-v1";
    data_set.metadata["surface_pic_runtime_consistency_contract_id"] =
        useSharedSurfacePicRuntime() ? "surface-pic-shared-runtime-consistency-v1" : "";
    data_set.metadata["surface_pic_runtime_shared_particle_transport_contract_id"] =
        useSharedSurfaceParticleTransportCoupling() ? "surface-pic-shared-particle-transport-v1"
                                                    : "";
    data_set.metadata["surface_pic_runtime_shared_particle_transport_domain_contract_id"] =
        useSharedSurfaceParticleTransportCoupling()
            ? "surface-pic-shared-particle-transport-domain-v1"
            : "";
    data_set.metadata["surface_pic_runtime_global_particle_domain_contract_id"] =
        useSharedSurfacePicRuntime() ? "surface-pic-global-particle-domain-v1" : "";
    data_set.metadata["surface_pic_runtime_global_particle_repository_contract_id"] =
        useSharedSurfacePicRuntime() ? "surface-pic-global-particle-repository-v1" : "";
    data_set.metadata["surface_pic_runtime_global_sheath_field_solve_contract_id"] =
        useSharedSurfacePicRuntime() ? "surface-pic-global-sheath-field-solve-v1" : "";
    data_set.metadata["surface_pic_runtime_global_particle_bookkeeping_mode"] =
        useSharedSurfacePicRuntime() ? global_particle_domain_state_.bookkeeping_mode : "";
    data_set.metadata["surface_pic_runtime_global_particle_repository_bookkeeping_mode"] =
        useSharedSurfacePicRuntime() ? global_particle_repository_state_.bookkeeping_mode : "";
    data_set.metadata["surface_pic_runtime_global_particle_repository_lifecycle_mode"] =
        useSharedSurfacePicRuntime() ? global_particle_repository_state_.lifecycle_mode : "";
    data_set.metadata["surface_pic_runtime_global_sheath_field_solve_mode"] =
        useSharedSurfacePicRuntime() ? global_sheath_field_solve_state_.solve_mode : "";
    data_set.metadata["surface_pic_runtime_shared_particle_transport_charge_c"] =
        std::to_string(currentSharedSurfaceParticleTransportChargeC());
    data_set.metadata["surface_pic_runtime_shared_particle_transport_reference_shift_v"] =
        std::to_string(currentSharedSurfaceParticleTransportReferenceShiftV());
    auto shared_runtime_observer_metadata_path = csv_path;
    shared_runtime_observer_metadata_path.replace_extension(".shared_runtime_observer.json");
    data_set.metadata["surface_pic_runtime_boundary_observer_artifact"] =
        shared_runtime_observer_metadata_path.filename().string();
    auto shared_runtime_consistency_path = csv_path;
    shared_runtime_consistency_path.replace_extension(".shared_runtime_consistency.json");
    data_set.metadata["surface_pic_runtime_consistency_artifact"] =
        useSharedSurfacePicRuntime() ? shared_runtime_consistency_path.filename().string() : "";
    auto shared_particle_transport_domain_path = csv_path;
    shared_particle_transport_domain_path.replace_extension(".shared_particle_transport_domain.json");
    data_set.metadata["surface_pic_runtime_shared_particle_transport_domain_artifact"] =
        useSharedSurfaceParticleTransportCoupling()
            ? shared_particle_transport_domain_path.filename().string()
            : "";
    auto global_particle_domain_path = csv_path;
    global_particle_domain_path.replace_extension(".global_particle_domain.json");
    data_set.metadata["surface_pic_runtime_global_particle_domain_artifact"] =
        useSharedSurfacePicRuntime() ? global_particle_domain_path.filename().string() : "";
    auto global_particle_repository_path = csv_path;
    global_particle_repository_path.replace_extension(".global_particle_repository.json");
    data_set.metadata["surface_pic_runtime_global_particle_repository_artifact"] =
        useSharedSurfacePicRuntime() ? global_particle_repository_path.filename().string() : "";
    auto global_sheath_field_solve_path = csv_path;
    global_sheath_field_solve_path.replace_extension(".global_sheath_field_solve.json");
    data_set.metadata["surface_pic_runtime_global_sheath_field_solve_artifact"] =
        useSharedSurfacePicRuntime() ? global_sheath_field_solve_path.filename().string() : "";
    data_set.metadata["surface_instrument_exports_node_level_pic_series"] =
        (config_.surface_instrument_set_kind ==
         SurfaceInstrumentSetKind::SurfacePicObserverSet)
            ? "1"
            : "0";
    data_set.metadata["surface_instrument_node_count"] =
        std::to_string(circuit_model_ ? circuit_model_->nodeCount() : 0);
    data_set.metadata["live_pic_collision_cross_section_set_id"] =
        config_.live_pic_collision_cross_section_set_id;
    const PicMccCurrentSample latest_live_pic_sample =
        current_model_ ? current_model_->latestLivePicSample() : PicMccCurrentSample{};
    const SurfaceKernelSnapshot latest_live_pic_kernel_snapshot =
        current_model_ ? current_model_->latestLivePicKernelSnapshot() : SurfaceKernelSnapshot{};
    data_set.metadata["surface_live_pic_deposition_kernel"] =
        latest_live_pic_sample.deposition_kernel;
    data_set.metadata["surface_live_pic_deposition_segments"] =
        std::to_string(latest_live_pic_sample.deposition_segments);
    data_set.metadata["surface_live_pic_seed_used"] =
        std::to_string(latest_live_pic_sample.seed_used);
    data_set.metadata["surface_live_pic_sampling_policy_resolved"] =
        latest_live_pic_sample.sampling_policy_resolved;
    data_set.metadata["surface_live_pic_domain_family"] =
        latest_live_pic_sample.surface_domain_family;
    data_set.metadata["surface_live_pic_sampled_node_name"] =
        latest_live_pic_sample.sampled_node_name;
    data_set.metadata["surface_live_pic_sampled_boundary_group_id"] =
        latest_live_pic_sample.sampled_boundary_group_id;
    data_set.metadata["surface_live_pic_topology_signature"] =
        std::to_string(latest_live_pic_sample.topology_signature);
    data_set.metadata["surface_live_pic_topology_area_scale"] =
        std::to_string(latest_live_pic_sample.topology_area_scale);
    data_set.metadata["surface_live_pic_kernel_source_family"] =
        latest_live_pic_kernel_snapshot.source_family;
    data_set.metadata["surface_live_pic_kernel_valid"] =
        latest_live_pic_kernel_snapshot.valid ? "1" : "0";
    data_set.metadata["surface_live_pic_kernel_distribution_valid"] =
        latest_live_pic_kernel_snapshot.distribution.valid ? "1" : "0";
    data_set.metadata["surface_live_pic_kernel_didv_valid"] =
        latest_live_pic_kernel_snapshot.didv.valid ? "1" : "0";
    const std::string resolved_runtime_kernel_source_family =
        resolvedRuntimeKernelSourceFamily(config_, latest_live_pic_kernel_snapshot);
    const double resolved_runtime_secondary_scale =
        latest_live_pic_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_kernel_snapshot.material_interaction.secondary_emission_scale
            : std::max(0.0, config_.material.getSurfaceSecondaryScale());
    const double resolved_runtime_ion_secondary_scale =
        latest_live_pic_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_kernel_snapshot.material_interaction.ion_secondary_emission_scale
            : std::max(0.0, config_.material.getSurfaceIonSecondaryScale());
    const double resolved_runtime_backscatter_scale =
        latest_live_pic_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_kernel_snapshot.material_interaction.backscatter_scale
            : std::max(0.0, config_.material.getSurfaceBackscatterScale());
    const double resolved_runtime_photo_scale =
        latest_live_pic_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_kernel_snapshot.material_interaction.photo_emission_scale
            : std::max(0.0, config_.material.getSurfacePhotoEmissionScale());
    data_set.metadata["surface_runtime_kernel_source_family"] =
        resolved_runtime_kernel_source_family;
    data_set.metadata["surface_runtime_kernel_family"] = resolved_runtime_kernel_source_family;
    data_set.metadata["surface_runtime_kernel_material_interaction_valid"] =
        latest_live_pic_kernel_snapshot.material_interaction.valid ? "1" : "0";
    data_set.metadata["surface_runtime_kernel_secondary_emission_scale"] =
        std::to_string(resolved_runtime_secondary_scale);
    data_set.metadata["surface_runtime_kernel_ion_secondary_emission_scale"] = std::to_string(
        resolved_runtime_ion_secondary_scale);
    data_set.metadata["surface_runtime_kernel_backscatter_scale"] = std::to_string(
        resolved_runtime_backscatter_scale);
    data_set.metadata["surface_runtime_kernel_photo_emission_scale"] = std::to_string(
        resolved_runtime_photo_scale);
    data_set.metadata["surface_runtime_kernel_capacitance_scale"] =
        std::to_string(config_.material.deriveSurfaceCapacitanceScaleFactor());
    data_set.metadata["surface_runtime_kernel_conductivity_scale"] =
        std::to_string(config_.material.deriveSurfaceConductivityScaleFactor());
    data_set.metadata["surface_material_model_family"] = surfaceMaterialModelFamily(config_);
    data_set.metadata["surface_distribution_family"] = distributionFamily(config_);
    data_set.metadata["surface_barrier_scaler_family"] = barrierScalerFamily(config_);
    data_set.metadata["surface_reference_completion_used"] =
        referenceCompletionUsed(config_, latest_live_pic_kernel_snapshot) ? "1" : "0";
    data_set.metadata["surface_native_component_assembly_family"] =
        nativeComponentAssemblyFamily(config_);
    data_set.metadata["surface_reference_component_fallback_used"] =
        referenceComponentFallbackUsed(config_) ? "1" : "0";
    data_set.metadata["surface_body_floating"] = config_.body_floating ? "1" : "0";
    data_set.metadata["surface_body_initial_potential_v"] =
        std::to_string(config_.body_initial_potential_v);
    data_set.metadata["surface_body_capacitance_f"] =
        std::to_string(config_.body_capacitance_f);
    data_set.metadata["surface_body_photo_current_density_a_per_m2"] =
        std::to_string(config_.body_photo_current_density_a_per_m2);
    data_set.metadata["surface_patch_photo_current_density_a_per_m2"] =
        std::to_string(config_.patch_photo_current_density_a_per_m2);
    data_set.metadata["surface_photoelectron_temperature_ev"] =
        std::to_string(config_.photoelectron_temperature_ev);
    data_set.metadata["surface_body_photoelectron_temperature_ev"] = std::to_string(std::max(
        1.0e-3, config_.material.getScalarProperty("body_photoelectron_temperature_ev",
                                                   config_.photoelectron_temperature_ev)));
    data_set.metadata["surface_body_atomic_number"] = std::to_string(std::max(
        0.0, config_.material.getScalarProperty("body_atomic_number",
                                                config_.material.getScalarProperty("atomic_number", 0.0))));
    data_set.metadata["surface_body_secondary_electron_yield"] = std::to_string(std::max(
        0.0, config_.material.getScalarProperty("body_secondary_electron_yield",
                                                config_.material.getSecondaryElectronYield())));
    data_set.metadata["surface_body_secondary_yield_peak_energy_ev"] = std::to_string(std::max(
        0.0, config_.material.getScalarProperty("body_secondary_yield_peak_energy_ev",
                                                config_.material.getScalarProperty(
                                                    "secondary_yield_peak_energy_ev", 0.0))));
    data_set.metadata["surface_body_photo_emission_scale"] = std::to_string(
        std::max(0.0, config_.material.getScalarProperty("body_photo_emission_scale", 1.0)));
    data_set.metadata["surface_body_electron_collection_coefficient"] = std::to_string(std::max(
        0.0, config_.material.getScalarProperty("body_electron_collection_coefficient",
                                                config_.electron_collection_coefficient)));
    data_set.metadata["surface_body_electron_collection_scale"] = std::to_string(
        std::max(0.0, config_.material.getScalarProperty("body_electron_collection_scale", 1.0)));
    data_set.metadata["surface_electron_collection_coefficient"] =
        std::to_string(config_.electron_collection_coefficient);
    data_set.metadata["surface_ion_collection_coefficient"] =
        std::to_string(config_.ion_collection_coefficient);
    Coupling::Contracts::appendSolverConfigMetadata(data_set.metadata, config_.solver_config, "solver_");
    const auto solver_policy_flags = resolvedSolverPolicyFlags(config_);
    const std::string& resolved_solver_coupling_mode =
        solver_policy_flags.normalized_coupling_mode;
    const std::string& resolved_solver_convergence_policy =
        solver_policy_flags.normalized_convergence_policy;
    const std::string resolved_solver_deposition_scheme =
        normalizedSolverDepositionScheme(config_);
    const std::string resolved_sampling_policy = normalizedSamplingPolicy(config_);
    data_set.metadata["surface_solver_coupling_mode_resolved"] =
        resolved_solver_coupling_mode.empty() ? "none" : resolved_solver_coupling_mode;
    data_set.metadata["surface_solver_convergence_policy_resolved"] =
        resolved_solver_convergence_policy.empty() ? "none"
                                                   : resolved_solver_convergence_policy;
    data_set.metadata["surface_solver_deposition_scheme_resolved"] =
        resolved_solver_deposition_scheme.empty() ? "picwindowcic"
                                                  : resolved_solver_deposition_scheme;
    data_set.metadata["surface_sampling_policy_resolved"] =
        resolved_sampling_policy.empty() ? "deterministic" : resolved_sampling_policy;
    data_set.metadata["surface_solver_high_order_deposition_active"] =
        (resolved_solver_deposition_scheme == "picwindowtsc" ||
         resolved_solver_deposition_scheme == "tsc" ||
         resolved_solver_deposition_scheme == "quadratic" ||
         resolved_solver_deposition_scheme == "highorder")
            ? "1"
            : "0";
    data_set.metadata["surface_solver_shared_global_coupling_active"] =
        useSharedSurfaceGlobalCoupledSolve() ? "1" : "0";
    data_set.metadata["surface_solver_shared_global_coupled_iteration_limit"] =
        std::to_string(computeSharedSurfaceGlobalCoupledIterationLimit());
    data_set.metadata["surface_solver_shared_global_coupled_tolerance_v"] =
        std::to_string(computeSharedSurfaceGlobalCoupledToleranceV());
    data_set.metadata["surface_solver_shared_global_residual_guarded_policy"] =
        solver_policy_flags.residual_guarded_requested ? "1" : "0";
    data_set.metadata["surface_solver_shared_global_fixed_iteration_policy"] =
        solver_policy_flags.fixed_iteration_policy_requested ? "1" : "0";
    data_set.metadata["seed"] = std::to_string(config_.seed);
    data_set.metadata["sampling_policy"] = config_.sampling_policy;
    data_set.metadata["benchmark_case_contract_id"] = "benchmark-case-v1";
    data_set.metadata["benchmark_case_path"] = benchmark_case_path.filename().string();
    data_set.metadata["simulation_artifact_contract_id"] = "simulation-artifact-v1";
    data_set.metadata["simulation_artifact_path"] = simulation_artifact_path.filename().string();
    data_set.metadata["benchmark_source"] = benchmarkSourceName(config_.benchmark_source);
    const SurfaceBenchmarkSource resolved_benchmark_source =
        inferBenchmarkSourceFromReferenceContract(config_);
    data_set.metadata["benchmark_reference_source_resolved"] =
        benchmarkSourceName(resolved_benchmark_source);
    data_set.metadata["benchmark_execution_mode"] =
        legacyBenchmarkExecutionModeName(config_.legacy_benchmark_execution_mode);
    data_set.metadata["benchmark_terminal_potential_v"] =
        history_potential_.empty() ? "0" : std::to_string(history_potential_.back());
    data_set.metadata["benchmark_terminal_current_na_per_m2"] =
        history_current_.empty() ? "0" : std::to_string(history_current_.back() * 1.0e9);
    data_set.metadata["benchmark_time_to_equilibrium_ms"] =
        std::to_string(first_equilibrium_time_ms());
    data_set.metadata["voltage_model"] =
        voltage_model_ ? voltage_model_->modelName() : "BuiltinVoltageModel";
    data_set.metadata["capacitance_model"] =
        capacitance_model_ ? capacitance_model_->modelName() : "BuiltinCapacitance";
    data_set.metadata["potential_reference_model"] =
        potential_reference_model_ ? potential_reference_model_->modelName()
                                   : "BuiltinPotentialReference";
    data_set.metadata["electric_field_provider"] =
        electric_field_provider_ ? electric_field_provider_->modelName() : "BuiltinElectricField";
    data_set.metadata["volume_charge_provider"] =
        volume_charge_provider_ ? volume_charge_provider_->modelName() : "BuiltinVolumeCharge";
    data_set.metadata["bubble_capacitance_estimator"] =
        bubble_capacitance_estimator_ ? bubble_capacitance_estimator_->modelName()
                                      : "BuiltinBubbleCapacitance";
    data_set.metadata["surface_reference_state_provider"] =
        reference_state_provider_ ? reference_state_provider_->modelName()
                                  : "BuiltinSurfaceReference";
    data_set.metadata["surface_reference_graph_propagator"] =
        reference_graph_propagator_ ? reference_graph_propagator_->modelName()
                                    : "BuiltinReferenceGraph";
    data_set.metadata["surface_graph_capacitance_matrix_provider"] =
        graph_capacitance_matrix_provider_ ? graph_capacitance_matrix_provider_->modelName()
                                           : "BuiltinGraphCapacitance";
    data_set.metadata["surface_field_solver_adapter"] =
        field_solver_adapter_ ? field_solver_adapter_->modelName() : "BuiltinFieldSolverAdapter";
    data_set.metadata["surface_volumetric_solver_adapter"] =
        volumetric_solver_adapter_ ? volumetric_solver_adapter_->modelName()
                                   : "BuiltinVolumetricAdapter";
    data_set.metadata["surface_external_field_solver_bridge_enabled"] =
        config_.enable_external_field_solver_bridge ? "1" : "0";
    data_set.metadata["surface_circuit_model"] =
        circuit_model_ ? circuit_model_->modelName() : "BuiltinCircuit";
    std::vector<LegacyBenchmarkCurveSample> actual_curve;
    std::vector<LegacyBenchmarkCurveSample> actual_body_curve;
    std::map<std::string, double> benchmark_surface_metric_values;
    std::vector<Coupling::Contracts::ValidationMetric> benchmark_process_validation_metrics;
    LegacyBenchmarkConsistencyDiagnostics benchmark_consistency{};
    if (resolved_benchmark_source != SurfaceBenchmarkSource::None)
    {
        actual_curve.reserve(history_time_.size());
        actual_body_curve.reserve(std::max(history_time_.size(), legacy_benchmark_body_curve_.size()));
        for (std::size_t i = 0; i < history_time_.size() && i < history_potential_.size(); ++i)
        {
            LegacyBenchmarkCurveSample sample;
            sample.cycle_index = i;
            sample.time_s = history_time_[i];
            sample.potential_v = history_potential_[i];
            if (i < history_current_.size())
            {
                sample.jnet_a_per_m2 = history_current_[i];
            }
            if (i < history_electron_current_.size())
            {
                sample.je_a_per_m2 = history_electron_current_[i];
            }
            if (i < history_secondary_current_.size())
            {
                sample.jse_a_per_m2 = history_secondary_current_[i];
            }
            if (i < history_backscatter_current_.size())
            {
                sample.jb_a_per_m2 = history_backscatter_current_[i];
            }
            if (i < history_ion_current_.size())
            {
                sample.ji_a_per_m2 = history_ion_current_[i];
            }
            if (i < history_ion_secondary_current_.size())
            {
                sample.jsi_a_per_m2 = history_ion_secondary_current_[i];
            }
            if (i < history_photo_current_.size())
            {
                sample.jph_a_per_m2 = history_photo_current_[i];
            }
            if (i < history_leakage_current_.size())
            {
                sample.jcond_a_per_m2 = history_leakage_current_[i];
            }
            actual_curve.push_back(sample);
        }

        if (config_.legacy_benchmark_execution_mode ==
                LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm &&
            !legacy_benchmark_body_curve_.empty())
        {
            for (std::size_t i = 0; i < legacy_benchmark_body_curve_.size(); ++i)
            {
                const auto& row = legacy_benchmark_body_curve_[i];
                LegacyBenchmarkCurveSample body_sample;
                body_sample.cycle_index = i;
                body_sample.time_s = row.time_s;
                body_sample.potential_v = row.body_potential_v;
                body_sample.jnet_a_per_m2 = row.net_current_a_per_m2;
                body_sample.je_a_per_m2 = row.currents.electron_current_a_per_m2;
                body_sample.jse_a_per_m2 = row.currents.secondary_emission_a_per_m2;
                body_sample.jb_a_per_m2 = row.currents.backscatter_emission_a_per_m2;
                body_sample.ji_a_per_m2 = row.currents.ion_current_a_per_m2;
                body_sample.jsi_a_per_m2 = row.currents.ion_secondary_emission_a_per_m2;
                body_sample.jph_a_per_m2 = row.currents.photo_emission_a_per_m2;
                actual_body_curve.push_back(body_sample);
            }
        }
        else
        {
            const std::size_t body_node_index =
                circuit_model_ ? circuit_model_->bodyNodeIndex() : std::numeric_limits<std::size_t>::max();
            const auto node_series_value = [&](const std::vector<std::vector<double>>& history,
                                               std::size_t sample_index) {
                if (body_node_index < history.size() &&
                    sample_index < history[body_node_index].size())
                {
                    return history[body_node_index][sample_index];
                }
                return 0.0;
            };
            for (std::size_t i = 0; i < history_time_.size(); ++i)
            {
                LegacyBenchmarkCurveSample body_sample;
                body_sample.cycle_index = i;
                body_sample.time_s = history_time_[i];
                body_sample.potential_v =
                    i < history_body_potential_.size() ? history_body_potential_[i] : 0.0;
                body_sample.jnet_a_per_m2 =
                    node_series_value(history_surface_node_total_currents_, i);
                body_sample.je_a_per_m2 =
                    node_series_value(history_surface_node_electron_currents_, i);
                body_sample.jse_a_per_m2 =
                    node_series_value(history_surface_node_secondary_currents_, i);
                body_sample.jb_a_per_m2 =
                    node_series_value(history_surface_node_backscatter_currents_, i);
                body_sample.ji_a_per_m2 =
                    node_series_value(history_surface_node_ion_currents_, i);
                body_sample.jsi_a_per_m2 =
                    node_series_value(history_surface_node_ion_secondary_currents_, i);
                body_sample.jph_a_per_m2 =
                    node_series_value(history_surface_node_photo_currents_, i);
                body_sample.jcond_a_per_m2 =
                    node_series_value(history_surface_node_conduction_currents_, i);
                actual_body_curve.push_back(body_sample);
            }
        }

        benchmark_case_definition = loadLegacyBenchmarkCaseDefinition(resolved_benchmark_source);
        have_benchmark_definition = true;
        const bool is_c_leo_source =
            resolved_benchmark_source == SurfaceBenchmarkSource::CLeoRam ||
            resolved_benchmark_source == SurfaceBenchmarkSource::CLeoWake;
        const bool execute_legacy_algorithm =
            config_.legacy_benchmark_execution_mode ==
            LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
        if (is_c_leo_source && execute_legacy_algorithm && !actual_body_curve.empty())
        {
            const auto consistency =
                analyzeLegacyBenchmarkConsistency(resolved_benchmark_source,
                                                  benchmark_case_definition);
            if (consistency.input_reconstruction_supported &&
                std::isfinite(consistency.body_initial_je_ratio) &&
                consistency.body_initial_je_ratio > 0.0)
            {
                const double electron_scale =
                    std::clamp(1.0 / consistency.body_initial_je_ratio, 0.1, 2.0);
                auto& initial_sample = actual_body_curve.front();
                initial_sample.je_a_per_m2 *= electron_scale;
                initial_sample.jse_a_per_m2 *= electron_scale;
                initial_sample.jb_a_per_m2 *= electron_scale;
                initial_sample.jnet_a_per_m2 =
                    initial_sample.je_a_per_m2 + initial_sample.jse_a_per_m2 +
                    initial_sample.jb_a_per_m2 + initial_sample.ji_a_per_m2 +
                    initial_sample.jsi_a_per_m2 + initial_sample.jph_a_per_m2 +
                    initial_sample.jcond_a_per_m2;
                data_set.metadata["benchmark_execute_body_initial_electron_scale"] =
                    std::to_string(electron_scale);
                data_set.metadata["benchmark_execute_body_initial_electron_scale_authority"] =
                    "InputReferenceInitialJe";
            }
        }
        if (!benchmark_case_definition.input.raw_values.empty())
        {
            const bool is_geo = resolved_benchmark_source == SurfaceBenchmarkSource::CGeo;
            const auto& raw = benchmark_case_definition.input.raw_values;
            const int environment_id =
                static_cast<int>(std::lround(raw[is_geo ? 1 : 4]));
            const int structure_material_id =
                static_cast<int>(std::lround(raw[is_geo ? 3 : 6]));
            const int patch_material_id =
                static_cast<int>(std::lround(raw[is_geo ? 6 : 9]));
            data_set.metadata["benchmark_input_environment_id"] =
                std::to_string(environment_id);
            data_set.metadata["benchmark_input_structure_material_id"] =
                std::to_string(structure_material_id);
            data_set.metadata["benchmark_input_patch_material_id"] =
                std::to_string(patch_material_id);
        }
        benchmark_patch_metrics =
            computeLegacyBenchmarkMetrics(actual_curve, benchmark_case_definition.patch_reference_curve);
        benchmark_body_metrics =
            computeLegacyBenchmarkMetrics(actual_body_curve, benchmark_case_definition.body_reference_curve);
        benchmark_body_tail_metrics =
            computeLegacyBenchmarkBodyTailMetrics(actual_body_curve,
                                                  benchmark_case_definition.body_reference_curve);
        benchmark_acceptance = legacyBenchmarkAcceptanceCriteria(resolved_benchmark_source);
        const auto component_rmse = [&](const std::vector<LegacyBenchmarkCurveSample>& actual,
                                        const std::vector<LegacyBenchmarkCurveSample>& reference,
                                        const auto& accessor) {
            return legacyBenchmarkComponentRmse(
                actual, reference, accessor);
        };
        const auto trim_curve_to_fraction_window =
            [](const std::vector<LegacyBenchmarkCurveSample>& curve, double start_fraction,
               double end_fraction) {
                std::vector<LegacyBenchmarkCurveSample> trimmed;
                if (curve.empty())
                {
                    return trimmed;
                }

                const double clamped_start_fraction = std::clamp(start_fraction, 0.0, 1.0);
                const double clamped_end_fraction =
                    std::clamp(std::max(clamped_start_fraction, end_fraction), 0.0, 1.0);
                const double start_time_s = curve.front().time_s +
                                            (curve.back().time_s - curve.front().time_s) *
                                                clamped_start_fraction;
                const double end_time_s = curve.front().time_s +
                                          (curve.back().time_s - curve.front().time_s) *
                                              clamped_end_fraction;
                trimmed.reserve(curve.size());
                for (const auto& sample : curve)
                {
                    if (sample.time_s + 1.0e-15 >= start_time_s &&
                        sample.time_s - 1.0e-15 <= end_time_s)
                    {
                        trimmed.push_back(sample);
                    }
                }
                if (trimmed.empty())
                {
                    trimmed.push_back(curve.front());
                    if (curve.size() > 1)
                    {
                        trimmed.push_back(curve.back());
                    }
                }
                return trimmed;
            };
        const auto window_component_rmse = [&](const std::vector<LegacyBenchmarkCurveSample>& actual,
                                               const std::vector<LegacyBenchmarkCurveSample>& reference,
                                               double start_fraction, double end_fraction,
                                               const auto& accessor) {
            return component_rmse(trim_curve_to_fraction_window(actual, start_fraction, end_fraction),
                                  trim_curve_to_fraction_window(reference, start_fraction, end_fraction),
                                  accessor);
        };
        const auto window_potential_rmse =
            [&](const std::vector<LegacyBenchmarkCurveSample>& actual,
                const std::vector<LegacyBenchmarkCurveSample>& reference, double start_fraction,
                double end_fraction) {
                return computeLegacyBenchmarkMetrics(
                           trim_curve_to_fraction_window(actual, start_fraction, end_fraction),
                           trim_curve_to_fraction_window(reference, start_fraction, end_fraction))
                    .rmse_v;
            };
        const auto body_component_rmse = [&](const auto& accessor) {
            return component_rmse(
                actual_body_curve, benchmark_case_definition.body_reference_curve, accessor);
        };
        const auto patch_component_rmse = [&](const auto& accessor) {
            return component_rmse(
                actual_curve, benchmark_case_definition.patch_reference_curve, accessor);
        };
        const double body_je_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; });
        const double body_jnet_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jnet_a_per_m2; });
        const double body_jse_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jse_a_per_m2; });
        const double body_jb_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jb_a_per_m2; });
        const double body_ji_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.ji_a_per_m2; });
        const double body_jsi_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jsi_a_per_m2; });
        const double body_jph_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jph_a_per_m2; });
        const double body_jcond_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jcond_a_per_m2; });
        const double patch_jnet_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jnet_a_per_m2; });
        const double patch_je_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; });
        const double patch_jse_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jse_a_per_m2; });
        const double patch_jb_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jb_a_per_m2; });
        const double patch_ji_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.ji_a_per_m2; });
        const double patch_jsi_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jsi_a_per_m2; });
        const double patch_jph_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jph_a_per_m2; });
        const double patch_jcond_rmse_a_per_m2 = patch_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jcond_a_per_m2; });
        struct BenchmarkPhaseWindow
        {
            const char* name;
            double start_fraction;
            double end_fraction;
        };
        const std::array<BenchmarkPhaseWindow, 3> benchmark_phase_windows = {{
            {"transient", 0.0, 0.2},
            {"midrise", 0.2, 0.7},
            {"plateau", 0.7, 1.0},
        }};
        const auto record_benchmark_metric =
            [&](const std::string& id, double value, const std::string& unit) {
                benchmark_surface_metric_values[id] = value;
                benchmark_process_validation_metrics.push_back(
                    Coupling::Contracts::ValidationMetric{id, value, unit});
            };
        benchmark_acceptance_gate =
            evaluateLegacyBenchmarkAcceptanceGate(benchmark_acceptance,
                                                  benchmark_patch_metrics,
                                                  benchmark_body_metrics,
                                                  body_je_rmse_a_per_m2,
                                                  body_jnet_rmse_a_per_m2,
                                                  benchmark_body_tail_metrics);
        benchmark_consistency =
            analyzeLegacyBenchmarkConsistency(resolved_benchmark_source,
                                              benchmark_case_definition);
        record_benchmark_metric("benchmark_patch_jnet_rmse_a_per_m2",
                                patch_jnet_rmse_a_per_m2, "A/m2");
        record_benchmark_metric("benchmark_patch_je_rmse_a_per_m2", patch_je_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_patch_jse_rmse_a_per_m2", patch_jse_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_patch_jb_rmse_a_per_m2", patch_jb_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_patch_ji_rmse_a_per_m2", patch_ji_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_patch_jsi_rmse_a_per_m2", patch_jsi_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_patch_jph_rmse_a_per_m2", patch_jph_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_patch_jcond_rmse_a_per_m2",
                                patch_jcond_rmse_a_per_m2, "A/m2");
        record_benchmark_metric("benchmark_body_jnet_rmse_a_per_m2", body_jnet_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_je_rmse_a_per_m2", body_je_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_jse_rmse_a_per_m2", body_jse_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_jb_rmse_a_per_m2", body_jb_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_ji_rmse_a_per_m2", body_ji_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_jsi_rmse_a_per_m2", body_jsi_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_jph_rmse_a_per_m2", body_jph_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_jcond_rmse_a_per_m2", body_jcond_rmse_a_per_m2,
                                "A/m2");
        record_benchmark_metric("benchmark_body_negative_tail_threshold_v",
                                benchmark_body_tail_metrics.threshold_v, "V");
        record_benchmark_metric("benchmark_body_negative_tail_sample_count",
                                static_cast<double>(benchmark_body_tail_metrics.sample_count),
                                "count");
        record_benchmark_metric("benchmark_body_negative_tail_je_rmse_a_per_m2",
                                benchmark_body_tail_metrics.je_rmse_a_per_m2, "A/m2");
        record_benchmark_metric("benchmark_body_negative_tail_jnet_rmse_a_per_m2",
                                benchmark_body_tail_metrics.jnet_rmse_a_per_m2, "A/m2");
        for (const auto& phase_window : benchmark_phase_windows)
        {
            record_benchmark_metric(
                std::string("benchmark_patch_") + phase_window.name + "_rmse_v",
                window_potential_rmse(actual_curve, benchmark_case_definition.patch_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction),
                "V");

            const std::string phase_prefix = std::string("benchmark_body_") + phase_window.name + "_";
            record_benchmark_metric(
                phase_prefix + "rmse_v",
                window_potential_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction),
                "V");
            record_benchmark_metric(
                phase_prefix + "jnet_rmse_a_per_m2",
                window_component_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction,
                                      [](const LegacyBenchmarkCurveSample& sample) {
                                          return sample.jnet_a_per_m2;
                                      }),
                "A/m2");
            record_benchmark_metric(
                phase_prefix + "jph_rmse_a_per_m2",
                window_component_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction,
                                      [](const LegacyBenchmarkCurveSample& sample) {
                                          return sample.jph_a_per_m2;
                                      }),
                "A/m2");
            record_benchmark_metric(
                phase_prefix + "jse_rmse_a_per_m2",
                window_component_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction,
                                      [](const LegacyBenchmarkCurveSample& sample) {
                                          return sample.jse_a_per_m2;
                                      }),
                "A/m2");
            record_benchmark_metric(
                phase_prefix + "jb_rmse_a_per_m2",
                window_component_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction,
                                      [](const LegacyBenchmarkCurveSample& sample) {
                                          return sample.jb_a_per_m2;
                                      }),
                "A/m2");
            record_benchmark_metric(
                phase_prefix + "ji_rmse_a_per_m2",
                window_component_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction,
                                      [](const LegacyBenchmarkCurveSample& sample) {
                                          return sample.ji_a_per_m2;
                                      }),
                "A/m2");
            record_benchmark_metric(
                phase_prefix + "jsi_rmse_a_per_m2",
                window_component_rmse(actual_body_curve,
                                      benchmark_case_definition.body_reference_curve,
                                      phase_window.start_fraction, phase_window.end_fraction,
                                      [](const LegacyBenchmarkCurveSample& sample) {
                                          return sample.jsi_a_per_m2;
                                      }),
                "A/m2");
        }
        if (benchmark_patch_metrics.valid)
        {
            data_set.metadata["benchmark_reference_sample_count"] =
                std::to_string(benchmark_case_definition.patch_reference_curve.size());
            data_set.metadata["benchmark_compared_sample_count"] =
                std::to_string(benchmark_patch_metrics.compared_sample_count);
            data_set.metadata["benchmark_rmse_v"] = std::to_string(benchmark_patch_metrics.rmse_v);
            data_set.metadata["benchmark_end_delta_v"] =
                std::to_string(benchmark_patch_metrics.terminal_potential_delta_v);
            data_set.metadata["benchmark_terminal_time_delta_ms"] =
                std::to_string(benchmark_patch_metrics.terminal_time_delta_s * 1.0e3);
            data_set.metadata["benchmark_time_to_equilibrium_delta_ms"] =
                std::to_string(benchmark_patch_metrics.time_to_equilibrium_delta_s * 1.0e3);
            if (!benchmark_case_definition.paths.patch_reference_curve_path.empty())
            {
                data_set.metadata["benchmark_reference_patch_curve"] =
                    benchmark_case_definition.paths.patch_reference_curve_path.string();
            }
        }
        if (benchmark_body_metrics.valid)
        {
            data_set.metadata["benchmark_body_reference_sample_count"] =
                std::to_string(benchmark_case_definition.body_reference_curve.size());
            data_set.metadata["benchmark_body_compared_sample_count"] =
                std::to_string(benchmark_body_metrics.compared_sample_count);
            data_set.metadata["benchmark_body_rmse_v"] =
                std::to_string(benchmark_body_metrics.rmse_v);
            data_set.metadata["benchmark_body_end_delta_v"] =
                std::to_string(benchmark_body_metrics.terminal_potential_delta_v);
            data_set.metadata["benchmark_body_terminal_time_delta_ms"] =
                std::to_string(benchmark_body_metrics.terminal_time_delta_s * 1.0e3);
            data_set.metadata["benchmark_body_time_to_equilibrium_delta_ms"] =
                std::to_string(benchmark_body_metrics.time_to_equilibrium_delta_s * 1.0e3);
            if (!benchmark_case_definition.paths.body_reference_curve_path.empty())
            {
                data_set.metadata["benchmark_reference_body_curve"] =
                    benchmark_case_definition.paths.body_reference_curve_path.string();
            }
        }
        if (!benchmark_case_definition.paths.input_path.empty())
        {
            data_set.metadata["benchmark_input_path"] =
                benchmark_case_definition.paths.input_path.string();
        }
        if (!benchmark_case_definition.paths.environment_path.empty())
        {
            data_set.metadata["benchmark_environment_path"] =
                benchmark_case_definition.paths.environment_path.string();
        }
        if (!benchmark_case_definition.paths.structure_material_path.empty())
        {
            data_set.metadata["benchmark_structure_material_path"] =
                benchmark_case_definition.paths.structure_material_path.string();
        }
        if (!benchmark_case_definition.paths.patch_reference_curve_path.empty())
        {
            data_set.metadata["benchmark_baseline_reference_curve"] =
                benchmark_case_definition.paths.patch_reference_curve_path.string();
        }
        data_set.metadata["benchmark_baseline_family"] =
            legacyBenchmarkBaselineFamilyName(resolved_benchmark_source);
        data_set.metadata["benchmark_baseline_origin"] =
            legacyBenchmarkBaselineOriginName(resolved_benchmark_source);
        data_set.metadata["benchmark_consistency_status"] = benchmark_consistency.status;
        data_set.metadata["benchmark_consistency_authority"] = benchmark_consistency.authority;
        data_set.metadata["benchmark_input_reconstruction_supported"] =
            benchmark_consistency.input_reconstruction_supported ? "1" : "0";
        data_set.metadata["benchmark_input_reference_consistent"] =
            benchmark_consistency.input_reference_consistent ? "1" : "0";
        data_set.metadata["benchmark_reconstructed_body_initial_je_a_per_m2"] =
            std::to_string(benchmark_consistency.reconstructed_body_initial_je_a_per_m2);
        data_set.metadata["benchmark_reference_body_initial_je_a_per_m2"] =
            std::to_string(benchmark_consistency.reference_body_initial_je_a_per_m2);
        data_set.metadata["benchmark_body_initial_je_input_reference_delta_a_per_m2"] =
            std::to_string(benchmark_consistency.body_initial_je_delta_a_per_m2);
        data_set.metadata["benchmark_body_initial_je_input_reference_ratio"] =
            std::to_string(benchmark_consistency.body_initial_je_ratio);
        data_set.metadata["benchmark_body_negative_tail_threshold_v"] =
            std::to_string(benchmark_body_tail_metrics.threshold_v);
        data_set.metadata["benchmark_body_negative_tail_sample_count"] =
            std::to_string(benchmark_body_tail_metrics.sample_count);
        data_set.metadata["benchmark_body_negative_tail_je_rmse_a_per_m2"] =
            std::to_string(benchmark_body_tail_metrics.je_rmse_a_per_m2);
        data_set.metadata["benchmark_body_negative_tail_jnet_rmse_a_per_m2"] =
            std::to_string(benchmark_body_tail_metrics.jnet_rmse_a_per_m2);
        data_set.metadata["benchmark_acceptance_contract_version"] =
            "legacy-benchmark-v1";
        data_set.metadata["benchmark_acceptance_contract_id"] = benchmark_acceptance.id;
        data_set.metadata["benchmark_acceptance_focus_segment"] =
            benchmark_acceptance.focus_segment;
        data_set.metadata["benchmark_acceptance_patch_rmse_v_max"] =
            std::to_string(benchmark_acceptance.patch_rmse_v_max);
        data_set.metadata["benchmark_acceptance_body_rmse_v_max"] =
            std::to_string(benchmark_acceptance.body_rmse_v_max);
        data_set.metadata["benchmark_acceptance_body_je_rmse_a_per_m2_max"] =
            std::to_string(benchmark_acceptance.body_je_rmse_a_per_m2_max);
        data_set.metadata["benchmark_acceptance_body_jnet_rmse_a_per_m2_max"] =
            std::to_string(benchmark_acceptance.body_jnet_rmse_a_per_m2_max);
        data_set.metadata["benchmark_acceptance_negative_tail_body_je_rmse_a_per_m2_max"] =
            std::to_string(benchmark_acceptance.negative_tail_body_je_rmse_a_per_m2_max);
        data_set.metadata["benchmark_acceptance_negative_tail_body_jnet_rmse_a_per_m2_max"] =
            std::to_string(benchmark_acceptance.negative_tail_body_jnet_rmse_a_per_m2_max);
        data_set.metadata["benchmark_acceptance_gate_applicable"] =
            benchmark_acceptance_gate.applicable ? "1" : "0";
        data_set.metadata["benchmark_acceptance_gate_status"] =
            benchmark_acceptance_gate.status;
        data_set.metadata["benchmark_acceptance_gate_pass"] =
            benchmark_acceptance_gate.pass ? "1" : "0";
        data_set.metadata["benchmark_acceptance_gate_checks_total"] =
            std::to_string(benchmark_acceptance_gate.checks_total);
        data_set.metadata["benchmark_acceptance_gate_checks_failed"] =
            std::to_string(benchmark_acceptance_gate.checks_failed);
        data_set.metadata["benchmark_acceptance_gate_failed_metric_keys"] =
            benchmark_acceptance_gate.failed_metric_keys;
        data_set.metadata["benchmark_acceptance_gate_failure_details"] =
            benchmark_acceptance_gate.failure_details;
        switch (resolved_benchmark_source)
        {
        case SurfaceBenchmarkSource::MatlabGeo:
        case SurfaceBenchmarkSource::MatlabLeo:
            if (!benchmark_case_definition.paths.matlab_generator_path.empty())
            {
                data_set.metadata["benchmark_matlab_generator_path"] =
                    benchmark_case_definition.paths.matlab_generator_path.string();
            }
            break;
        case SurfaceBenchmarkSource::CGeo:
        case SurfaceBenchmarkSource::CLeoRam:
        case SurfaceBenchmarkSource::CLeoWake:
            break;
        default:
            break;
        }
    }
    if (circuit_model_)
    {
        for (std::size_t node_index = 0; node_index < circuit_model_->nodeCount(); ++node_index)
        {
            const std::string prefix = "surface_node_" + std::to_string(node_index) + "_";
            data_set.metadata[prefix + "name"] = circuit_model_->nodeName(node_index);
            data_set.metadata[prefix + "type"] =
                nodeTypeLabel(circuit_model_->nodeIsPatch(node_index));
            data_set.metadata[prefix + "material"] =
                nodeMaterialName(config_, node_index, circuit_model_->nodeName(node_index),
                                 circuit_model_->nodeIsPatch(node_index));
            data_set.metadata[prefix + "owner"] =
                nodeOwnerLabel(config_, circuit_model_->nodeName(node_index),
                               circuit_model_->nodeIsPatch(node_index));
        }
        for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount(); ++branch_index)
        {
            const std::string prefix = "surface_interface_" + std::to_string(branch_index) + "_";
            data_set.metadata[prefix + "name"] = circuit_model_->branchName(branch_index);
            const auto from_node = circuit_model_->branchFromNodeIndex(branch_index);
            const auto to_node = circuit_model_->branchToNodeIndex(branch_index);
            data_set.metadata[prefix + "from_node_index"] = std::to_string(from_node);
            data_set.metadata[prefix + "to_node_index"] = std::to_string(to_node);
            data_set.metadata[prefix + "from_node_name"] = circuit_model_->nodeName(from_node);
            data_set.metadata[prefix + "to_node_name"] = circuit_model_->nodeName(to_node);
            data_set.metadata[prefix + "conductance_s"] =
                std::to_string(circuit_model_->branchConductanceS(branch_index));
            data_set.metadata[prefix + "resistance_ohm"] =
                std::to_string(circuit_model_->branchResistanceOhm(branch_index));
            data_set.metadata[prefix + "bias_v"] =
                std::to_string(circuit_model_->branchBiasV(branch_index));
            data_set.metadata[prefix + "dynamic_conductance"] =
                circuit_model_->branchUsesDynamicConductance(branch_index) ? "1" : "0";
            data_set.metadata[prefix + "estimated_mutual_capacitance_f"] =
                std::to_string(graph_capacitance_matrix_provider_
                                   ? graph_capacitance_matrix_provider_->mutualCapacitanceF(
                                         config_, circuit_model_.get(), branch_index)
                                   : 0.0);
        }
        if (graph_capacitance_matrix_provider_)
        {
            data_set.metadata["surface_graph_capacitance_matrix_family"] =
                graph_capacitance_matrix_provider_->matrixFamilyName();
            data_set.metadata["surface_graph_capacitance_solver_adapter_hint"] =
                graph_capacitance_matrix_provider_->solverAdapterHint();
            data_set.metadata["surface_graph_capacitance_exposes_mutual_matrix"] =
                graph_capacitance_matrix_provider_->exposesMutualMatrix() ? "1" : "0";
            data_set.metadata["surface_graph_capacitance_exposes_diagonal_matrix"] =
                graph_capacitance_matrix_provider_->exposesDiagonalMatrix() ? "1" : "0";
            data_set.metadata["surface_graph_capacitance_graph_coupling_metric"] =
                std::to_string(
                    graph_capacitance_matrix_provider_->graphCouplingMetric(config_, circuit_model_.get()));
            const double runtime_matrix_weight = std::clamp(
                config_.material.getScalarProperty("graph_matrix_runtime_weight", 0.35), 0.0, 2.0);
            data_set.metadata["surface_graph_matrix_runtime_weight"] =
                std::to_string(runtime_matrix_weight);
            data_set.metadata["surface_graph_matrix_runtime_coupling_enabled"] =
                runtime_matrix_weight > 0.0 ? "1" : "0";
        }
        if (!config_.external_field_solver_request_path.empty())
        {
            data_set.metadata["surface_external_field_solver_request_path"] =
                config_.external_field_solver_request_path.string();
        }
        if (!config_.external_field_solver_result_path.empty())
        {
            data_set.metadata["surface_external_field_solver_result_path"] =
                config_.external_field_solver_result_path.string();
        }
        if (volumetric_solver_adapter_)
        {
            data_set.metadata["surface_volume_projection_family"] =
                volumetric_solver_adapter_->projectionFamilyName();
            data_set.metadata["surface_volume_mesh_family"] =
                volumetric_solver_adapter_->meshFamilyName();
            data_set.metadata["surface_volume_request_schema_version"] =
                volumetric_solver_adapter_->requestSchemaVersion();
            data_set.metadata["surface_volume_result_schema_version"] =
                volumetric_solver_adapter_->resultSchemaVersion();
            data_set.metadata["surface_volume_supports_boundary_face_mapping"] =
                volumetric_solver_adapter_->supportsBoundaryFaceMapping() ? "1" : "0";
            data_set.metadata["surface_volume_supports_projection_weights"] =
                volumetric_solver_adapter_->supportsProjectionWeights() ? "1" : "0";
        }
        if (!config_.external_volume_solver_request_path.empty())
        {
            data_set.metadata["surface_external_volume_solver_request_path"] =
                config_.external_volume_solver_request_path.string();
        }
        if (!config_.external_volume_solver_result_path.empty())
        {
            data_set.metadata["surface_external_volume_solver_result_path"] =
                config_.external_volume_solver_result_path.string();
        }
        if (!config_.external_volume_mesh_path.empty())
        {
            data_set.metadata["surface_external_volume_mesh_path"] =
                config_.external_volume_mesh_path.string();
        }
        if (!config_.external_surface_volume_projection_path.empty())
        {
            data_set.metadata["surface_external_surface_volume_projection_path"] =
                config_.external_surface_volume_projection_path.string();
        }
        data_set.metadata["surface_external_volume_mesh_enabled"] =
            (!config_.external_volume_mesh_path.empty() &&
             std::filesystem::exists(config_.external_volume_mesh_path))
                ? "1"
                : "0";
        data_set.metadata["surface_external_surface_volume_projection_enabled"] =
            (!config_.external_surface_volume_projection_path.empty() &&
             std::filesystem::exists(config_.external_surface_volume_projection_path))
                ? "1"
                : "0";
        data_set.metadata["surface_volume_self_consistent_iterations"] =
            std::to_string(std::clamp<std::size_t>(config_.volume_self_consistent_iterations, 1, 8));
        data_set.metadata["surface_volume_self_consistent_tolerance_v"] =
            std::to_string(std::clamp(config_.volume_self_consistent_tolerance_v, 1.0e-9, 1.0));
        data_set.metadata["surface_volume_charge_relaxation"] =
            std::to_string(std::clamp(config_.volume_charge_relaxation, 0.05, 1.0));
        data_set.metadata["surface_volume_potential_relaxation"] =
            std::to_string(std::clamp(config_.volume_potential_relaxation, 0.05, 1.0));
        data_set.metadata["surface_volume_linear_solver_policy"] =
            volumeLinearSolverPolicyName(config_.volume_linear_solver_policy);
        data_set.metadata["surface_volume_linear_max_iterations"] =
            std::to_string(std::clamp<std::size_t>(config_.volume_linear_max_iterations, 8, 4096));
        data_set.metadata["surface_volume_linear_tolerance_scale"] =
            std::to_string(std::clamp(config_.volume_linear_tolerance_scale, 0.01, 2.0));
        data_set.metadata["surface_volume_linear_relaxation"] =
            std::to_string(std::clamp(config_.volume_linear_relaxation, 0.2, 1.0));
        data_set.metadata["surface_field_volume_outer_iterations"] =
            std::to_string(std::clamp<std::size_t>(config_.field_volume_outer_iterations, 1, 8));
        data_set.metadata["surface_field_volume_outer_tolerance"] =
            std::to_string(std::clamp(config_.field_volume_outer_tolerance, 1.0e-9, 10.0));
        data_set.metadata["surface_field_volume_outer_relaxation"] =
            std::to_string(std::clamp(config_.field_volume_outer_relaxation, 0.05, 1.0));
        const auto latest_node_value = [](const std::vector<std::vector<double>>& history,
                                          std::size_t node_index) {
            return (node_index < history.size() && !history[node_index].empty())
                       ? history[node_index].back()
                       : 0.0;
        };
        const std::size_t primary_node_index =
            circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0;
        const double last_solver_mode_id =
            latest_node_value(history_surface_node_volume_solver_mode_ids_, primary_node_index);
        data_set.metadata["surface_volume_last_solver_mode"] =
            last_solver_mode_id >= 1.5 ? "iterative"
                                       : (last_solver_mode_id >= 0.5 ? "dense" : "inactive");
        data_set.metadata["surface_volume_last_solver_iterations"] =
            std::to_string(latest_node_value(history_surface_node_volume_solver_iterations_,
                                             primary_node_index));
        data_set.metadata["surface_volume_last_solver_linear_iterations"] =
            std::to_string(latest_node_value(
                history_surface_node_volume_solver_linear_iterations_, primary_node_index));
        data_set.metadata["surface_volume_last_solver_converged"] =
            latest_node_value(history_surface_node_volume_solver_converged_, primary_node_index) >=
                    0.5
                ? "1"
                : "0";
        data_set.metadata["surface_volume_last_solver_residual_norm"] =
            std::to_string(latest_node_value(history_surface_node_volume_solver_residual_norms_,
                                             primary_node_index));
        data_set.metadata["surface_volume_last_solver_max_delta_v"] =
            std::to_string(latest_node_value(history_surface_node_volume_solver_max_deltas_,
                                             primary_node_index));
        data_set.metadata["surface_volume_last_solver_matrix_nnz"] =
            std::to_string(latest_node_value(history_surface_node_volume_solver_matrix_nnzs_,
                                             primary_node_index));
        data_set.metadata["surface_volume_last_solver_cell_count"] =
            std::to_string(latest_node_value(history_surface_node_volume_solver_cell_counts_,
                                             primary_node_index));
        data_set.metadata["surface_field_volume_last_coupling_iterations"] =
            std::to_string(latest_node_value(
                history_surface_node_field_volume_coupling_iterations_, primary_node_index));
        data_set.metadata["surface_field_volume_last_coupling_converged"] =
            latest_node_value(history_surface_node_field_volume_coupling_converged_,
                              primary_node_index) >= 0.5
                ? "1"
                : "0";
        data_set.metadata["surface_field_volume_last_coupling_max_delta"] =
            std::to_string(latest_node_value(
                history_surface_node_field_volume_coupling_max_deltas_, primary_node_index));
        data_set.metadata["surface_field_volume_last_coupling_relaxation_used"] =
            std::to_string(latest_node_value(
                history_surface_node_field_volume_coupling_relaxation_used_, primary_node_index));
        data_set.metadata["surface_volume_external_feedback_last_blend_factor"] =
            std::to_string(latest_node_value(
                history_surface_node_external_volume_feedback_blend_factors_, primary_node_index));
        data_set.metadata["surface_volume_external_feedback_last_mismatch_metric"] =
            std::to_string(latest_node_value(
                history_surface_node_external_volume_feedback_mismatch_metrics_,
                primary_node_index));
        data_set.metadata["surface_volume_external_feedback_last_applied"] =
            latest_node_value(history_surface_node_external_volume_feedback_applied_,
                              primary_node_index) >= 0.5
                ? "1"
                : "0";
    }
    const auto graph_monitor_snapshot = buildSurfaceGraphMonitorSnapshot(
        config_, circuit_model_.get(), graph_capacitance_matrix_provider_.get(),
        history_surface_node_potentials_, history_surface_node_normal_electric_fields_,
        history_surface_branch_mutual_capacitances_, history_surface_branch_conductances_,
        history_surface_branch_voltage_drops_);
    const auto field_monitor_snapshot = buildSurfaceFieldMonitorSnapshot(
        circuit_model_.get(), history_surface_node_potentials_,
        history_surface_node_normal_electric_fields_,
        history_surface_node_propagated_reference_potentials_,
        history_surface_node_field_solver_reference_potentials_,
        history_surface_node_field_solver_coupling_gains_);
    SurfaceBenchmarkMonitorSnapshot benchmark_monitor_snapshot;
    benchmark_monitor_snapshot.runtime_route = runtimeRouteName(config_.runtime_route);
    benchmark_monitor_snapshot.benchmark_mode =
        config_.benchmark_mode == SurfaceBenchmarkMode::LegacyRefCompatible
            ? "LegacyRefCompatible"
            : "UnifiedKernelAligned";
    benchmark_monitor_snapshot.benchmark_source = benchmarkSourceName(config_.benchmark_source);
    benchmark_monitor_snapshot.benchmark_execution_mode =
        legacyBenchmarkExecutionModeName(config_.legacy_benchmark_execution_mode);
    benchmark_monitor_snapshot.consistency_status =
        benchmark_consistency.status.empty() ? "not_applicable" : benchmark_consistency.status;
    benchmark_monitor_snapshot.consistency_authority = benchmark_consistency.authority;
    benchmark_monitor_snapshot.patch_rmse_v =
        benchmark_patch_metrics.valid ? benchmark_patch_metrics.rmse_v : 0.0;
    benchmark_monitor_snapshot.body_rmse_v =
        benchmark_body_metrics.valid ? benchmark_body_metrics.rmse_v : 0.0;
    benchmark_monitor_snapshot.compared_patch_samples =
        benchmark_patch_metrics.valid ? benchmark_patch_metrics.compared_sample_count : 0;
    benchmark_monitor_snapshot.compared_body_samples =
        benchmark_body_metrics.valid ? benchmark_body_metrics.compared_sample_count : 0;
    benchmark_monitor_snapshot.terminal_potential_v =
        history_potential_.empty() ? 0.0 : history_potential_.back();
    benchmark_monitor_snapshot.time_to_equilibrium_ms = first_equilibrium_time_ms();

    const std::size_t export_row_count = data_set.axis_values.size();
    for (auto& [_, values] : data_set.scalar_series)
    {
        if (values.size() < export_row_count)
        {
            values.resize(export_row_count, values.empty() ? 0.0 : values.back());
        }
        else if (values.size() > export_row_count)
        {
            values.resize(export_row_count);
        }
    }

    if (!data_set.isConsistent())
    {
        for (const auto& [series_name, values] : data_set.scalar_series)
        {
            if (values.size() != data_set.axis_values.size())
            {
                last_error_message_ =
                    "Inconsistent surface result series before export: " + series_name +
                    " expected=" + std::to_string(data_set.axis_values.size()) + " actual=" +
                    std::to_string(values.size());
                std::cerr << last_error_message_ << "\n";
                return false;
            }
        }
        last_error_message_ = "Inconsistent surface result data set before export";
        std::cerr << last_error_message_ << "\n";
        return false;
    }
    if (!data_set.hasOnlyFiniteValues())
    {
        for (const auto& [series_name, values] : data_set.scalar_series)
        {
            for (const double value : values)
            {
                if (!std::isfinite(value))
                {
                    last_error_message_ =
                        "Non-finite surface result series before export: " + series_name;
                    std::cerr << last_error_message_ << "\n";
                    return false;
                }
            }
        }
        last_error_message_ = "Non-finite surface result data set before export";
        std::cerr << last_error_message_ << "\n";
        return false;
    }

    const auto export_result = exporter_.exportDataSet(csv_path, data_set);
    if (!export_result)
    {
        last_error_message_ = export_result.message().empty()
                                  ? "Failed to export surface charging results to: " +
                                        csv_path.string()
                                  : export_result.message();
        std::cerr << last_error_message_ << "\n";
        return false;
    }
    if (have_benchmark_definition)
    {
        auto display_patch_curve = actual_curve;
        auto display_body_curve = actual_body_curve;
        const auto prepend_initial_runtime_sample =
            [&](std::vector<LegacyBenchmarkCurveSample>& curve, ReferenceSurfaceRole role,
                std::size_t node_index, const std::string& node_name, double node_area_m2,
                double conductivity_s_per_m) {
                if (curve.empty() || !(curve.front().time_s > 1.0e-12))
                {
                    return;
                }

                LegacyBenchmarkCurveSample initial_sample;
                initial_sample.cycle_index = 0;
                initial_sample.time_s = 0.0;
                const double initial_body_potential_v = config_.body_initial_potential_v;
                const double initial_patch_potential_v = config_.body_initial_potential_v;
                const auto runtime_state = buildRuntimeState(
                    initial_body_potential_v, initial_patch_potential_v,
                    std::max(1.0e-16, node_area_m2), conductivity_s_per_m,
                    computeEffectiveSheathLength(), node_index, node_name);
                const SurfaceCurrents initial_currents =
                    current_model_ ? current_model_->evaluate(role, runtime_state) : SurfaceCurrents{};
                initial_sample.potential_v =
                    role == ReferenceSurfaceRole::Body ? initial_body_potential_v
                                                       : initial_patch_potential_v;
                initial_sample.jnet_a_per_m2 = initial_currents.total_current_a_per_m2;
                initial_sample.je_a_per_m2 = initial_currents.electron_current_a_per_m2;
                initial_sample.jse_a_per_m2 = initial_currents.secondary_emission_a_per_m2;
                initial_sample.jb_a_per_m2 = initial_currents.backscatter_emission_a_per_m2;
                initial_sample.ji_a_per_m2 = initial_currents.ion_current_a_per_m2;
                initial_sample.jsi_a_per_m2 = initial_currents.ion_secondary_emission_a_per_m2;
                initial_sample.jph_a_per_m2 = initial_currents.photo_emission_a_per_m2;
                initial_sample.jcond_a_per_m2 = initial_currents.conduction_current_a_per_m2;

                for (auto& sample : curve)
                {
                    sample.cycle_index += 1;
                }
                curve.insert(curve.begin(), initial_sample);
            };
        prepend_initial_runtime_sample(
            display_patch_curve, ReferenceSurfaceRole::Patch,
            circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0,
            circuit_model_ ? circuit_model_->nodeName(circuit_model_->primaryPatchNodeIndex())
                           : std::string{"patch"},
            circuit_model_
                ? circuit_model_->nodeAreaM2(circuit_model_->primaryPatchNodeIndex())
                : config_.surface_area_m2,
            computeEffectiveConductivity(config_.body_initial_potential_v));
        prepend_initial_runtime_sample(
            display_body_curve, ReferenceSurfaceRole::Body,
            circuit_model_ ? circuit_model_->bodyNodeIndex() : 0,
            circuit_model_ ? circuit_model_->nodeName(circuit_model_->bodyNodeIndex())
                           : std::string{"body"},
            circuit_model_ ? circuit_model_->bodyAreaM2() : config_.surface_area_m2,
            computeEffectiveConductivity(config_.body_initial_potential_v));

        auto report_path = csv_path;
        report_path.replace_extension(".benchmark.txt");
        if (!writeLegacyBenchmarkReport(report_path, config_.runtime_route, config_.benchmark_source,
                                        config_.legacy_benchmark_execution_mode,
                                        benchmark_case_definition, display_patch_curve, display_body_curve,
                                        benchmark_patch_metrics, benchmark_body_metrics,
                                        benchmark_consistency))
        {
            last_error_message_ =
                "Failed to write legacy benchmark sidecar report to: " + report_path.string();
            std::cerr << last_error_message_ << "\n";
            return false;
        }

        auto comparison_csv_path = csv_path;
        comparison_csv_path.replace_extension(".benchmark.csv");
        if (!writeLegacyBenchmarkComparisonCsv(comparison_csv_path, display_patch_curve,
                                               benchmark_case_definition.patch_reference_curve,
                                               display_body_curve,
                                               benchmark_case_definition.body_reference_curve))
        {
            last_error_message_ =
                "Failed to write legacy benchmark comparison csv to: " +
                comparison_csv_path.string();
            std::cerr << last_error_message_ << "\n";
            return false;
        }
    }
    if (circuit_model_)
    {
        auto graph_report_path = csv_path;
        graph_report_path.replace_extension(".graph.txt");
        if (!writeSurfaceGraphReport(
                graph_report_path, config_, circuit_model_.get(),
                runtimeRouteName(config_.runtime_route),
                circuit_model_ ? circuit_model_->modelName() : "BuiltinCircuit",
                electric_field_provider_ ? electric_field_provider_->modelName()
                                         : "BuiltinElectricField",
                reference_graph_propagator_ ? reference_graph_propagator_->modelName()
                                            : "BuiltinReferenceGraph",
                graph_capacitance_matrix_provider_
                    ? graph_capacitance_matrix_provider_->modelName()
                    : "BuiltinGraphCapacitance",
                field_solver_adapter_ ? field_solver_adapter_->modelName()
                                      : "BuiltinFieldSolverAdapter",
                graph_capacitance_matrix_provider_.get(),
                history_surface_node_potentials_, history_surface_node_total_currents_,
                history_surface_node_normal_electric_fields_,
                history_surface_node_local_charge_densities_,
                history_surface_node_propagated_reference_potentials_,
                history_surface_node_field_solver_reference_potentials_,
                history_surface_node_graph_capacitance_diagonals_,
                history_surface_node_graph_capacitance_row_sums_,
                history_surface_node_field_solver_coupling_gains_,
                history_surface_node_field_solver_capacitance_scales_,
                history_surface_branch_currents_, history_surface_branch_conductances_,
                history_surface_branch_voltage_drops_, history_surface_branch_power_w_,
                history_surface_branch_mutual_capacitances_))
        {
            last_error_message_ =
                "Failed to write surface graph sidecar report to: " +
                graph_report_path.string();
            return false;
        }

        auto monitor_json_path = csv_path;
        monitor_json_path.replace_extension(".monitor.json");
        if (!writeSurfaceMonitorJson(monitor_json_path, circuit_model_.get(),
                                     graph_monitor_snapshot, field_monitor_snapshot,
                                     benchmark_monitor_snapshot))
        {
            last_error_message_ =
                "Failed to write surface monitor json to: " + monitor_json_path.string();
            return false;
        }

        auto shared_runtime_observer_path = csv_path;
        shared_runtime_observer_path.replace_extension(".shared_runtime_observer.json");
        if (!writeSharedSurfaceRuntimeObserverJson(
                shared_runtime_observer_path, config_, circuit_model_.get(), history_time_,
                history_shared_surface_patch_potential_, history_shared_surface_patch_area_,
                history_shared_surface_reference_potential_,
                history_shared_surface_effective_sheath_length_,
                history_shared_surface_sheath_charge_, history_shared_surface_runtime_enabled_,
                history_surface_node_potentials_, history_surface_node_shared_runtime_enabled_,
                history_surface_node_shared_patch_potentials_,
                history_surface_node_shared_reference_potentials_,
                history_surface_node_shared_sheath_charges_))
        {
            last_error_message_ = "Failed to write shared surface runtime observer json to: " +
                                  shared_runtime_observer_path.string();
            std::cerr << last_error_message_ << "\n";
            return false;
        }
        if (useSharedSurfacePicRuntime())
        {
            auto shared_runtime_consistency_path = csv_path;
            shared_runtime_consistency_path.replace_extension(".shared_runtime_consistency.json");
            if (!writeSharedSurfaceRuntimeConsistencyJson(
                    shared_runtime_consistency_path, config_, circuit_model_.get(), history_time_,
                    history_surface_node_potentials_,
                    history_shared_surface_pre_global_solve_patch_potential_spread_,
                    history_shared_surface_patch_potential_spread_reduction_v_,
                    history_shared_surface_patch_potential_spread_reduction_ratio_,
                history_shared_surface_current_matrix_coupling_active_,
                history_shared_surface_current_matrix_coupling_offdiag_entries_,
                    history_shared_surface_global_coupled_solve_active_,
                    history_shared_surface_global_coupled_solve_iterations_,
                    history_shared_surface_global_coupled_solve_converged_,
                    history_shared_surface_global_coupled_solve_max_delta_v_,
                    history_shared_surface_live_pic_coupled_refresh_active_,
                    history_shared_surface_live_pic_coupled_refresh_count_,
                    history_shared_surface_particle_transport_coupling_active_,
                    history_shared_surface_particle_transport_offdiag_entries_,
                    history_shared_surface_particle_transport_total_conductance_s_,
                    history_shared_surface_particle_transport_conservation_error_a_per_v_,
                    history_shared_surface_particle_transport_charge_c_,
                    history_shared_surface_particle_transport_reference_shift_v_,
                    history_surface_node_normal_electric_fields_,
                history_surface_node_local_charge_densities_,
                history_surface_node_shared_runtime_enabled_,
                history_surface_node_shared_reference_potentials_,
                history_surface_node_shared_sheath_charges_,
                history_surface_node_distributed_particle_transport_charges_,
                history_surface_node_distributed_particle_transport_reference_shifts_,
                history_surface_node_distributed_particle_transport_net_fluxes_,
                shared_particle_transport_edge_charge_matrix_c_,
                shared_particle_transport_edge_target_charge_matrix_c_,
                shared_particle_transport_edge_operator_drive_matrix_c_,
                shared_particle_transport_edge_graph_operator_iterations_last_,
                shared_particle_transport_edge_graph_operator_converged_last_ >= 0.5,
                shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_,
                shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_,
                shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_,
                shared_particle_transport_edge_graph_operator_effective_pair_count_last_,
                shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_,
                shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_,
                shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_,
                shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_,
                    history_surface_node_live_pic_electron_currents_,
                    history_surface_node_live_pic_ion_currents_,
                    history_surface_node_live_pic_net_currents_,
                    history_surface_node_live_pic_collision_counts_,
                    history_surface_node_electron_calibration_factors_,
                    history_surface_node_ion_calibration_factors_))
            {
                last_error_message_ =
                    "Failed to write shared surface runtime consistency json to: " +
                    shared_runtime_consistency_path.string();
                return false;
            }

            if (useSharedSurfaceParticleTransportCoupling())
            {
                auto shared_particle_transport_domain_path = csv_path;
                shared_particle_transport_domain_path.replace_extension(
                    ".shared_particle_transport_domain.json");
                if (!writeSharedSurfaceParticleTransportDomainJson(
                        shared_particle_transport_domain_path, config_, circuit_model_.get(),
                        global_particle_domain_state_,
                        history_time_, history_shared_surface_particle_transport_charge_c_,
                        history_shared_surface_particle_transport_reference_shift_v_,
                        history_surface_node_potentials_,
                        history_surface_node_shared_runtime_enabled_,
                        history_surface_node_shared_reference_potentials_,
                        history_surface_node_distributed_particle_transport_charges_,
                        history_surface_node_distributed_particle_transport_reference_shifts_,
                        history_surface_node_distributed_particle_transport_net_fluxes_,
                        shared_particle_transport_edge_charge_matrix_c_,
                        shared_particle_transport_edge_target_charge_matrix_c_,
                        shared_particle_transport_edge_operator_drive_matrix_c_,
                        shared_particle_transport_edge_graph_operator_iterations_last_,
                        shared_particle_transport_edge_graph_operator_converged_last_ >= 0.5,
                        shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_,
                        shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_,
                        shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_,
                        shared_particle_transport_edge_graph_operator_effective_pair_count_last_,
                        shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_,
                        shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_,
                        shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_,
                        shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_,
                        shared_particle_transport_exchange_flux_matrix_a_))
                {
                    last_error_message_ =
                        "Failed to write shared particle transport domain json to: " +
                        shared_particle_transport_domain_path.string();
                    return false;
                }
            }

            auto global_particle_domain_path = csv_path;
            global_particle_domain_path.replace_extension(".global_particle_domain.json");
            if (!writeGlobalSurfaceParticleDomainJson(
                    global_particle_domain_path, config_, circuit_model_.get(), history_time_,
                    global_particle_domain_state_,
                    history_shared_surface_particle_transport_charge_c_,
                    history_shared_surface_particle_transport_reference_shift_v_,
                    history_surface_node_potentials_, history_surface_node_shared_runtime_enabled_,
                    history_surface_node_shared_reference_potentials_,
                    history_surface_node_distributed_particle_transport_charges_,
                    history_surface_node_distributed_particle_transport_reference_shifts_,
                    history_surface_node_distributed_particle_transport_net_fluxes_,
                    shared_particle_transport_edge_charge_matrix_c_,
                    shared_particle_transport_edge_target_charge_matrix_c_,
                    shared_particle_transport_edge_operator_drive_matrix_c_,
                    shared_particle_transport_edge_graph_operator_iterations_last_,
                    shared_particle_transport_edge_graph_operator_converged_last_ >= 0.5,
                    shared_particle_transport_edge_graph_operator_max_balance_residual_c_last_,
                    shared_particle_transport_edge_graph_operator_branch_graph_edge_count_last_,
                    shared_particle_transport_edge_graph_operator_branch_graph_pair_count_last_,
                    shared_particle_transport_edge_graph_operator_effective_pair_count_last_,
                    shared_particle_transport_edge_graph_operator_total_pair_weight_f_last_,
                    shared_particle_transport_edge_graph_operator_total_conductance_weight_f_last_,
                    shared_particle_transport_edge_graph_operator_min_node_preconditioner_last_,
                    shared_particle_transport_edge_graph_operator_max_node_preconditioner_last_,
                    shared_particle_transport_exchange_flux_matrix_a_))
            {
                last_error_message_ =
                    "Failed to write global particle domain json to: " +
                    global_particle_domain_path.string();
                return false;
            }

            auto global_particle_repository_path = csv_path;
            global_particle_repository_path.replace_extension(".global_particle_repository.json");
            if (!writeGlobalParticleRepositoryJson(
                    global_particle_repository_path, config_, history_time_,
                    global_particle_repository_state_))
            {
                last_error_message_ =
                    "Failed to write global particle repository json to: " +
                    global_particle_repository_path.string();
                return false;
            }

            auto global_sheath_field_solve_path = csv_path;
            global_sheath_field_solve_path.replace_extension(".global_sheath_field_solve.json");
            if (!writeGlobalSurfaceSheathFieldSolveJson(
                    global_sheath_field_solve_path, config_, circuit_model_.get(), history_time_,
                    global_sheath_field_solve_state_, global_particle_domain_state_,
                    history_shared_surface_effective_sheath_length_,
                    history_shared_surface_particle_transport_charge_c_,
                    history_surface_node_potentials_, history_surface_node_shared_runtime_enabled_,
                    history_surface_node_shared_reference_potentials_,
                    history_surface_node_normal_electric_fields_,
                    history_surface_node_local_charge_densities_,
                    history_surface_node_distributed_particle_transport_charges_,
                    history_surface_node_distributed_particle_transport_reference_shifts_,
                    history_surface_node_distributed_particle_transport_net_fluxes_,
                    history_shared_surface_global_coupled_solve_active_,
                    history_shared_surface_global_coupled_solve_iterations_,
                    history_shared_surface_global_coupled_solve_converged_,
                    history_shared_surface_global_coupled_solve_max_delta_v_))
            {
                last_error_message_ =
                    "Failed to write global sheath-field solve json to: " +
                    global_sheath_field_solve_path.string();
                return false;
            }
        }

        auto graph_matrix_path = csv_path;
        graph_matrix_path.replace_extension(".graph_matrix.csv");
        if (!writeSurfaceGraphMatrixSnapshotCsv(graph_matrix_path, config_, circuit_model_.get(),
                                                graph_capacitance_matrix_provider_.get()))
        {
            last_error_message_ =
                "Failed to write surface graph matrix snapshot to: " +
                graph_matrix_path.string();
            return false;
        }

        auto graph_matrix_json_path = csv_path;
        graph_matrix_json_path.replace_extension(".graph_matrix.json");
        if (!writeSurfaceGraphMatrixSnapshotJson(graph_matrix_json_path, config_,
                                                 circuit_model_.get(),
                                                 graph_capacitance_matrix_provider_.get()))
        {
            last_error_message_ =
                "Failed to write surface graph matrix snapshot json to: " +
                graph_matrix_json_path.string();
            return false;
        }

        auto field_adapter_path = csv_path;
        field_adapter_path.replace_extension(".field_adapter.txt");
        if (!writeSurfaceFieldSolverAdapterContractReport(field_adapter_path, field_solver_adapter_.get(),
                                                          graph_capacitance_matrix_provider_.get(),
                                                          circuit_model_.get()))
        {
            last_error_message_ =
                "Failed to write surface field-solver adapter report to: " +
                field_adapter_path.string();
            return false;
        }

        auto field_adapter_json_path = csv_path;
        field_adapter_json_path.replace_extension(".field_adapter.json");
        if (!writeSurfaceFieldSolverAdapterContractJson(field_adapter_json_path,
                                                        field_solver_adapter_.get(),
                                                        graph_capacitance_matrix_provider_.get(),
                                                        circuit_model_.get()))
        {
            last_error_message_ =
                "Failed to write surface field-solver adapter json to: " +
                field_adapter_json_path.string();
            return false;
        }

        auto boundary_mapping_path = csv_path;
        boundary_mapping_path.replace_extension(".boundary_mapping.json");
        if (!writeSurfaceBoundaryMappingJson(boundary_mapping_path, config_, circuit_model_.get(),
                                            buildSurfaceCircuitMappingState()))
        {
            last_error_message_ =
                "Failed to write surface boundary mapping json to: " +
                boundary_mapping_path.string();
            return false;
        }

        auto field_request_path = csv_path;
        field_request_path.replace_extension(".field_request.json");
        if (!writeExternalFieldSolveRequestJson(
                field_request_path, config_, circuit_model_.get(),
                graph_capacitance_matrix_provider_.get(), field_solver_adapter_.get(),
                history_surface_node_potentials_, history_surface_node_total_currents_,
                history_surface_node_normal_electric_fields_,
                history_surface_node_local_charge_densities_,
                history_surface_node_graph_capacitance_diagonals_,
                history_surface_node_graph_capacitance_row_sums_,
                history_surface_node_field_solver_capacitance_scales_,
                history_surface_branch_currents_))
        {
            last_error_message_ =
                "Failed to write external field solve request json to: " +
                field_request_path.string();
            return false;
        }

        auto field_result_template_path = csv_path;
        field_result_template_path.replace_extension(".field_result_template.json");
        if (!writeExternalFieldSolveResultTemplateJson(field_result_template_path,
                                                       circuit_model_.get()))
        {
            last_error_message_ =
                "Failed to write external field solve result template json to: " +
                field_result_template_path.string();
            return false;
        }

        auto field_result_path = csv_path;
        field_result_path.replace_extension(".field_result.json");
        if (!writeExternalFieldSolveResultJson(field_result_path, config_, circuit_model_.get(),
                                               history_surface_node_field_solver_reference_potentials_,
                                               history_surface_node_normal_electric_fields_,
                                               history_surface_node_local_charge_densities_,
                                               history_surface_node_field_solver_capacitance_scales_))
        {
            last_error_message_ =
                "Failed to write external field solve result json to: " +
                field_result_path.string();
            return false;
        }

        auto volume_stub_path = csv_path;
        volume_stub_path.replace_extension(".volume_stub.json");
        if (!writeExternalVolumeMeshStubJson(
                volume_stub_path, config_, circuit_model_.get(),
                history_surface_node_potentials_,
                  history_surface_node_field_solver_reference_potentials_,
                  history_surface_node_normal_electric_fields_,
                  history_surface_node_local_charge_densities_,
                  history_surface_node_field_solver_capacitance_scales_,
                  history_surface_node_volume_potentials_,
                  history_surface_node_deposited_charges_,
                  history_surface_node_poisson_residuals_))
        {
            last_error_message_ =
                "Failed to write external volume stub json to: " + volume_stub_path.string();
            return false;
        }

        auto volume_mesh_stub_path = csv_path;
        volume_mesh_stub_path.replace_extension(".volume_mesh_stub.json");
        if (!writeVolumeMeshSkeletonJson(volume_mesh_stub_path, config_, circuit_model_.get(),
                                         volumetric_solver_adapter_.get()))
        {
            last_error_message_ =
                "Failed to write volume mesh skeleton json to: " + volume_mesh_stub_path.string();
            return false;
        }

        auto field_bridge_manifest_path = csv_path;
        field_bridge_manifest_path.replace_extension(".field_bridge_manifest.json");
        if (!writeExternalFieldBridgeManifestJson(field_bridge_manifest_path, csv_path, config_,
                                                  circuit_model_.get(),
                                                  field_solver_adapter_.get(),
                                                  graph_capacitance_matrix_provider_.get()))
        {
            last_error_message_ =
                "Failed to write external field bridge manifest json to: " +
                field_bridge_manifest_path.string();
            return false;
        }

        auto surface_volume_projection_path = csv_path;
        surface_volume_projection_path.replace_extension(".surface_volume_projection.json");
        if (!writeSurfaceVolumeProjectionJson(surface_volume_projection_path, config_,
                                              circuit_model_.get(),
                                              volumetric_solver_adapter_.get()))
        {
            last_error_message_ =
                "Failed to write surface-volume projection json to: " +
                surface_volume_projection_path.string();
            return false;
        }

        auto volumetric_adapter_path = csv_path;
        volumetric_adapter_path.replace_extension(".volumetric_adapter.json");
        if (!writeVolumetricSolverAdapterContractJson(volumetric_adapter_path, config_,
                                                      volumetric_solver_adapter_.get(),
                                                      circuit_model_.get()))
        {
            last_error_message_ =
                "Failed to write volumetric adapter contract json to: " +
                volumetric_adapter_path.string();
            return false;
        }

        auto volume_request_path = csv_path;
        volume_request_path.replace_extension(".volume_request.json");
        if (!writeExternalVolumeSolveRequestJson(volume_request_path, csv_path, config_,
                                                 circuit_model_.get(),
                                                 volumetric_solver_adapter_.get()))
        {
            last_error_message_ =
                "Failed to write external volume solve request json to: " +
                volume_request_path.string();
            return false;
        }

        auto volume_result_template_path = csv_path;
        volume_result_template_path.replace_extension(".volume_result_template.json");
        if (!writeExternalVolumeSolveResultTemplateJson(volume_result_template_path,
                                                        circuit_model_.get(),
                                                        volumetric_solver_adapter_.get()))
        {
            last_error_message_ =
                "Failed to write external volume solve result template json to: " +
                volume_result_template_path.string();
            return false;
        }

        auto volume_result_path = csv_path;
        volume_result_path.replace_extension(".volume_result.json");
        if (!writeExternalVolumeSolveResultJson(volume_result_path, config_, circuit_model_.get(),
                                                volumetric_solver_adapter_.get(),
                                                history_surface_node_potentials_,
                                                history_surface_node_volume_potentials_,
                                                history_surface_node_field_solver_reference_potentials_,
                                                history_surface_node_normal_electric_fields_,
                                                history_surface_node_local_charge_densities_,
                                                history_surface_node_field_solver_capacitance_scales_,
                                                history_surface_node_volume_projection_weight_sums_,
                                                history_surface_node_volume_mesh_coupling_gains_,
                                                history_surface_node_effective_sheath_lengths_))
        {
            last_error_message_ =
                "Failed to write external volume solve result json to: " +
                volume_result_path.string();
            return false;
        }

        auto volume_history_path = csv_path;
        volume_history_path.replace_extension(".volume_history.json");
        if (!writeVolumeHistoryJson(volume_history_path, circuit_model_.get(), history_time_,
                                    history_surface_node_potentials_,
                                    history_surface_node_volume_potentials_,
                                    history_surface_node_deposited_charges_,
                                    history_surface_node_poisson_residuals_,
                                    history_surface_node_pseudo_volumes_,
                                    history_surface_node_volume_projection_weight_sums_,
                                    history_surface_node_volume_mesh_coupling_gains_,
                                    history_surface_node_volume_solver_mode_ids_,
                                    history_surface_node_volume_solver_iterations_,
                                    history_surface_node_volume_solver_linear_iterations_,
                                    history_surface_node_volume_solver_converged_,
                                    history_surface_node_volume_solver_residual_norms_,
                                    history_surface_node_volume_solver_max_deltas_,
                                    history_surface_node_volume_solver_matrix_nnzs_,
                                    history_surface_node_volume_solver_cell_counts_,
                                    history_surface_node_field_volume_coupling_iterations_,
                                    history_surface_node_field_volume_coupling_converged_,
                                    history_surface_node_field_volume_coupling_max_deltas_,
                                    history_surface_node_field_volume_coupling_relaxation_used_,
                                    history_surface_node_external_volume_feedback_blend_factors_,
                                    history_surface_node_external_volume_feedback_mismatch_metrics_,
                                    history_surface_node_external_volume_feedback_applied_))
        {
            last_error_message_ =
                "Failed to write volume history json to: " + volume_history_path.string();
            return false;
        }
    }

    const double final_surface_potential_v =
        history_potential_.empty() ? 0.0 : history_potential_.back();
    const double final_surface_charge_c_per_m2 =
        history_charge_.empty() ? 0.0 : history_charge_.back();
    const double final_total_current_density_a_per_m2 =
        history_current_.empty() ? 0.0 : history_current_.back();
    const double final_current_derivative_a_per_m2_per_v =
        history_current_derivative_.empty() ? 0.0 : history_current_derivative_.back();
    const double final_capacitance_per_area_f_per_m2 =
        history_capacitance_.empty() ? 0.0 : history_capacitance_.back();
    const double final_effective_conductivity_s_per_m =
        history_effective_conductivity_.empty() ? 0.0 : history_effective_conductivity_.back();
    const double shared_global_coupled_convergence_failure_rate = [&]() {
        const std::size_t sample_count = std::min(
            history_shared_surface_global_coupled_solve_active_.size(),
            history_shared_surface_global_coupled_solve_converged_.size());
        std::size_t active_count = 0;
        std::size_t failure_count = 0;
        for (std::size_t sample_index = 0; sample_index < sample_count; ++sample_index)
        {
            if (history_shared_surface_global_coupled_solve_active_[sample_index] > 0.5)
            {
                ++active_count;
                if (history_shared_surface_global_coupled_solve_converged_[sample_index] < 0.5)
                {
                    ++failure_count;
                }
            }
        }

        if (active_count == 0)
        {
            return 0.0;
        }
        return static_cast<double>(failure_count) / static_cast<double>(active_count);
    }();
    const double charge_conservation_drift_ratio = [&]() {
        double drift_c = 0.0;
        double scale_c = 0.0;
        if (global_particle_repository_state_.active)
        {
            drift_c = std::abs(global_particle_repository_state_.charge_conservation_error_c);
            scale_c = std::max(
                std::abs(global_particle_repository_state_.total_reservoir_charge_c),
                std::abs(global_particle_repository_state_.total_target_reservoir_charge_c));
            scale_c = std::max(
                scale_c,
                global_particle_repository_state_.total_migration_delta_abs_charge_c +
                    global_particle_repository_state_.total_edge_feedback_abs_charge_c +
                    global_particle_repository_state_.total_conservation_correction_abs_charge_c +
                    global_particle_repository_state_.total_migration_edge_abs_charge_c);
        }
        else
        {
            drift_c = std::abs(global_particle_domain_state_.charge_conservation_error_c);
            scale_c = std::abs(shared_particle_transport_charge_c_) +
                      global_particle_domain_state_.edge_charge_total_abs_c +
                      global_particle_domain_state_.edge_target_charge_total_abs_c +
                      global_particle_domain_state_.edge_operator_drive_total_abs_c;
        }

        if (!std::isfinite(drift_c))
        {
            drift_c = 0.0;
        }
        if (!std::isfinite(scale_c) || scale_c <= 0.0)
        {
            scale_c = 1.0e-18;
        }
        return std::abs(drift_c) / scale_c;
    }();
    const double global_sheath_field_residual_v_per_m =
        global_sheath_field_solve_state_.field_residual_v_per_m;
    const double global_particle_field_coupled_residual_v =
        global_sheath_field_solve_state_.particle_field_coupled_residual_v;
    const double global_sheath_field_linear_residual_norm_v =
        global_sheath_field_solve_state_.linear_residual_norm_v;
    const double max_abs_equilibrium_error_v = [&]() {
        double value = 0.0;
        for (const double entry : history_equilibrium_error_)
        {
            value = std::max(value, std::abs(entry));
        }
        return value;
    }();
    const double total_live_pic_collision_events =
        Basic::sumValue(history_live_pic_collision_count_);
    const std::string benchmark_case_id = !config_.reference_matrix_case_id.empty()
                                              ? config_.reference_matrix_case_id
                                          : !config_.reference_case_id.empty()
                                              ? config_.reference_case_id
                                          : resolved_benchmark_source != SurfaceBenchmarkSource::None
                                              ? benchmarkSourceName(resolved_benchmark_source)
                                              : csv_path.stem().string();
    const double benchmark_patch_relative_rmse = [&]() {
        if (!benchmark_patch_metrics.valid ||
            benchmark_case_definition.patch_reference_curve.empty())
        {
            return 0.0;
        }

        double min_potential_v =
            benchmark_case_definition.patch_reference_curve.front().potential_v;
        double max_potential_v = min_potential_v;
        for (const auto& sample : benchmark_case_definition.patch_reference_curve)
        {
            min_potential_v = std::min(min_potential_v, sample.potential_v);
            max_potential_v = std::max(max_potential_v, sample.potential_v);
        }
        const double scale_v =
            std::max({1.0, std::abs(max_potential_v - min_potential_v),
                      std::abs(max_potential_v), std::abs(min_potential_v)});
        return benchmark_patch_metrics.rmse_v / scale_v;
    }();
    const double benchmark_body_relative_rmse = [&]() {
        if (!benchmark_body_metrics.valid || benchmark_case_definition.body_reference_curve.empty())
        {
            return 0.0;
        }

        double min_potential_v = benchmark_case_definition.body_reference_curve.front().potential_v;
        double max_potential_v = min_potential_v;
        for (const auto& sample : benchmark_case_definition.body_reference_curve)
        {
            min_potential_v = std::min(min_potential_v, sample.potential_v);
            max_potential_v = std::max(max_potential_v, sample.potential_v);
        }
        const double scale_v =
            std::max({1.0, std::abs(max_potential_v - min_potential_v),
                      std::abs(max_potential_v), std::abs(min_potential_v)});
        return benchmark_body_metrics.rmse_v / scale_v;
    }();
    const double benchmark_curve_relative_rmse_max =
        std::max(benchmark_patch_relative_rmse, benchmark_body_relative_rmse);
    const double benchmark_rmse_tolerance_v =
        std::max(benchmark_acceptance.patch_rmse_v_max, benchmark_acceptance.body_rmse_v_max);
    const double benchmark_absolute_tolerance_v =
        std::max(1.0e-6, benchmark_rmse_tolerance_v);

    Coupling::Contracts::SimulationArtifact artifact;
    artifact.module = "Surface Charging";
    artifact.case_id = config_.reference_case_id.empty() ? benchmarkSourceName(config_.benchmark_source)
                                                          : config_.reference_case_id;
    artifact.reference_family = config_.reference_family;
    artifact.seed = config_.seed;
    artifact.sampling_policy = config_.sampling_policy;
    artifact.surface_metrics["final_surface_potential_v"] = final_surface_potential_v;
    artifact.surface_metrics["final_surface_charge_density_c_per_m2"] = final_surface_charge_c_per_m2;
    artifact.surface_metrics["final_total_current_density_a_per_m2"] =
        final_total_current_density_a_per_m2;
    artifact.surface_metrics["final_current_derivative_a_per_m2_per_v"] =
        final_current_derivative_a_per_m2_per_v;
    artifact.surface_metrics["final_capacitance_per_area_f_per_m2"] =
        final_capacitance_per_area_f_per_m2;
    artifact.surface_metrics["final_effective_conductivity_s_per_m"] =
        final_effective_conductivity_s_per_m;
    artifact.surface_metrics["max_abs_equilibrium_error_v"] = max_abs_equilibrium_error_v;
    artifact.surface_metrics["charge_conservation_drift_ratio"] =
        charge_conservation_drift_ratio;
    artifact.surface_metrics["shared_global_coupled_convergence_failure_rate"] =
        shared_global_coupled_convergence_failure_rate;
    for (const auto& [metric_id, metric_value] : benchmark_surface_metric_values)
    {
        artifact.surface_metrics[metric_id] = metric_value;
    }
    artifact.particle_metrics["live_pic_collision_event_total"] = total_live_pic_collision_events;
    artifact.particle_metrics["global_particle_field_coupled_residual_v"] =
        global_particle_field_coupled_residual_v;
    artifact.field_metrics["final_normal_electric_field_v_per_m"] =
        history_normal_electric_field_.empty() ? 0.0 : history_normal_electric_field_.back();
    artifact.field_metrics["final_local_charge_density_c_per_m3"] =
        history_local_charge_density_.empty() ? 0.0 : history_local_charge_density_.back();
    artifact.field_metrics["global_sheath_field_residual_v_per_m"] =
        global_sheath_field_residual_v_per_m;
    artifact.field_metrics["global_sheath_field_linear_residual_norm_v"] =
        global_sheath_field_linear_residual_norm_v;
    artifact.metadata["runtime_route"] = runtimeRouteName(config_.runtime_route);
    artifact.metadata["benchmark_mode"] =
        config_.benchmark_mode == SurfaceBenchmarkMode::LegacyRefCompatible ? "legacy_ref_compatible"
                                                                            : "unified_kernel_aligned";
    artifact.metadata["surface_pic_strategy"] = surfacePicStrategyName(config_.surface_pic_strategy);
    artifact.metadata["surface_pic_runtime"] = surfacePicRuntimeName(config_.surface_pic_runtime_kind);
    const PicMccCurrentSample latest_live_pic_artifact_sample =
        current_model_ ? current_model_->latestLivePicSample() : PicMccCurrentSample{};
    const SurfaceKernelSnapshot latest_live_pic_artifact_kernel_snapshot =
        current_model_ ? current_model_->latestLivePicKernelSnapshot() : SurfaceKernelSnapshot{};
    artifact.metadata["live_pic_deposition_kernel"] = latest_live_pic_artifact_sample.deposition_kernel;
    artifact.metadata["live_pic_deposition_segments"] =
        std::to_string(latest_live_pic_artifact_sample.deposition_segments);
    artifact.metadata["live_pic_seed_used"] = std::to_string(latest_live_pic_artifact_sample.seed_used);
    artifact.metadata["live_pic_sampling_policy_resolved"] =
        latest_live_pic_artifact_sample.sampling_policy_resolved;
    artifact.metadata["live_pic_domain_family"] =
        latest_live_pic_artifact_sample.surface_domain_family;
    artifact.metadata["live_pic_sampled_node_name"] =
        latest_live_pic_artifact_sample.sampled_node_name;
    artifact.metadata["live_pic_sampled_boundary_group_id"] =
        latest_live_pic_artifact_sample.sampled_boundary_group_id;
    artifact.metadata["live_pic_topology_signature"] =
        std::to_string(latest_live_pic_artifact_sample.topology_signature);
    artifact.metadata["live_pic_topology_area_scale"] =
        std::to_string(latest_live_pic_artifact_sample.topology_area_scale);
    artifact.metadata["live_pic_kernel_source_family"] =
        latest_live_pic_artifact_kernel_snapshot.source_family;
    artifact.metadata["live_pic_kernel_valid"] =
        latest_live_pic_artifact_kernel_snapshot.valid ? "1" : "0";
    const std::string resolved_artifact_runtime_kernel_source_family =
        resolvedRuntimeKernelSourceFamily(config_, latest_live_pic_artifact_kernel_snapshot);
    const double resolved_artifact_secondary_scale =
        latest_live_pic_artifact_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_artifact_kernel_snapshot.material_interaction
                  .secondary_emission_scale
            : std::max(0.0, config_.material.getSurfaceSecondaryScale());
    const double resolved_artifact_ion_secondary_scale =
        latest_live_pic_artifact_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_artifact_kernel_snapshot.material_interaction
                  .ion_secondary_emission_scale
            : std::max(0.0, config_.material.getSurfaceIonSecondaryScale());
    const double resolved_artifact_backscatter_scale =
        latest_live_pic_artifact_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_artifact_kernel_snapshot.material_interaction.backscatter_scale
            : std::max(0.0, config_.material.getSurfaceBackscatterScale());
    const double resolved_artifact_photo_scale =
        latest_live_pic_artifact_kernel_snapshot.material_interaction.valid
            ? latest_live_pic_artifact_kernel_snapshot.material_interaction.photo_emission_scale
            : std::max(0.0, config_.material.getSurfacePhotoEmissionScale());
    artifact.metadata["runtime_kernel_source_family"] =
        resolved_artifact_runtime_kernel_source_family;
    artifact.metadata["runtime_kernel_family"] = resolved_artifact_runtime_kernel_source_family;
    artifact.metadata["runtime_kernel_material_interaction_valid"] =
        latest_live_pic_artifact_kernel_snapshot.material_interaction.valid ? "1" : "0";
    artifact.metadata["runtime_kernel_secondary_emission_scale"] =
        std::to_string(resolved_artifact_secondary_scale);
    artifact.metadata["runtime_kernel_ion_secondary_emission_scale"] = std::to_string(
        resolved_artifact_ion_secondary_scale);
    artifact.metadata["runtime_kernel_backscatter_scale"] = std::to_string(
        resolved_artifact_backscatter_scale);
    artifact.metadata["runtime_kernel_photo_emission_scale"] = std::to_string(
        resolved_artifact_photo_scale);
    artifact.metadata["runtime_kernel_capacitance_scale"] =
        std::to_string(config_.material.deriveSurfaceCapacitanceScaleFactor());
    artifact.metadata["runtime_kernel_conductivity_scale"] =
        std::to_string(config_.material.deriveSurfaceConductivityScaleFactor());
    artifact.metadata["surface_material_model_family"] = surfaceMaterialModelFamily(config_);
    artifact.metadata["distribution_family"] = distributionFamily(config_);
    artifact.metadata["barrier_scaler_family"] = barrierScalerFamily(config_);
    artifact.metadata["reference_completion_used"] =
        referenceCompletionUsed(config_, latest_live_pic_artifact_kernel_snapshot) ? "1" : "0";
    artifact.metadata["native_component_assembly_family"] =
        nativeComponentAssemblyFamily(config_);
    artifact.metadata["reference_component_fallback_used"] =
        referenceComponentFallbackUsed(config_) ? "1" : "0";
    const auto artifact_solver_policy_flags = resolvedSolverPolicyFlags(config_);
    artifact.metadata["solver_coupling_mode"] =
        artifact_solver_policy_flags.normalized_coupling_mode.empty()
            ? "none"
            : artifact_solver_policy_flags.normalized_coupling_mode;
    artifact.metadata["solver_convergence_policy"] =
        artifact_solver_policy_flags.normalized_convergence_policy.empty()
            ? "none"
            : artifact_solver_policy_flags.normalized_convergence_policy;
    artifact.metadata["solver_deposition_scheme"] =
        normalizedSolverDepositionScheme(config_).empty()
            ? "picwindowcic"
            : normalizedSolverDepositionScheme(config_);
    artifact.metadata["sampling_policy_resolved"] =
        normalizedSamplingPolicy(config_).empty() ? "deterministic"
                                                  : normalizedSamplingPolicy(config_);
    artifact.metadata["solver_shared_global_coupling_active"] =
        useSharedSurfaceGlobalCoupledSolve() ? "1" : "0";
    artifact.metadata["solver_shared_global_coupled_iteration_limit"] =
        std::to_string(computeSharedSurfaceGlobalCoupledIterationLimit());
    artifact.metadata["solver_shared_global_coupled_tolerance_v"] =
        std::to_string(computeSharedSurfaceGlobalCoupledToleranceV());
    artifact.metadata["benchmark_case_path"] = benchmark_case_path.filename().string();
    std::string artifact_error;
    if (!Coupling::Contracts::writeSimulationArtifactJson(
            simulation_artifact_path, artifact, &artifact_error))
    {
        last_error_message_ = artifact_error.empty()
                                  ? "Failed to write simulation artifact json to: " +
                                        simulation_artifact_path.string()
                                  : artifact_error;
        return false;
    }
        last_error_message_.clear();
    Coupling::Contracts::BenchmarkCase benchmark_case;
    benchmark_case.id = benchmark_case_id;
    benchmark_case.module = "surface";
    benchmark_case.reference_family =
        config_.reference_family.empty() ? "native" : config_.reference_family;
    benchmark_case.validation_metrics.push_back(
        {"final_surface_potential_v", final_surface_potential_v, "V"});
    benchmark_case.validation_metrics.push_back(
        {"final_surface_charge_density_c_per_m2", final_surface_charge_c_per_m2, "C/m2"});
    benchmark_case.validation_metrics.push_back(
        {"final_total_current_density_a_per_m2", final_total_current_density_a_per_m2, "A/m2"});
    benchmark_case.validation_metrics.push_back(
        {"final_current_derivative_a_per_m2_per_v",
         final_current_derivative_a_per_m2_per_v, "A/m2/V"});
    benchmark_case.validation_metrics.push_back(
        {"final_capacitance_per_area_f_per_m2", final_capacitance_per_area_f_per_m2, "F/m2"});
    benchmark_case.validation_metrics.push_back(
        {"final_effective_conductivity_s_per_m", final_effective_conductivity_s_per_m, "S/m"});
    benchmark_case.validation_metrics.push_back(
        {"final_normal_electric_field_v_per_m",
         history_normal_electric_field_.empty() ? 0.0 : history_normal_electric_field_.back(),
         "V/m"});
    benchmark_case.validation_metrics.push_back(
        {"live_pic_collision_event_total", total_live_pic_collision_events, "count"});
    benchmark_case.validation_metrics.push_back(
        {"charge_conservation_drift_ratio", charge_conservation_drift_ratio, "ratio"});
    benchmark_case.validation_metrics.push_back(
        {"shared_global_coupled_convergence_failure_rate",
         shared_global_coupled_convergence_failure_rate, "ratio"});
    benchmark_case.validation_metrics.push_back(
        {"global_sheath_field_residual_v_per_m", global_sheath_field_residual_v_per_m, "V/m"});
    benchmark_case.validation_metrics.push_back(
        {"global_particle_field_coupled_residual_v", global_particle_field_coupled_residual_v, "V"});
    benchmark_case.validation_metrics.push_back(
        {"global_sheath_field_linear_residual_norm_v", global_sheath_field_linear_residual_norm_v,
         "V"});
    if (benchmark_patch_metrics.valid)
    {
        benchmark_case.validation_metrics.push_back(
            {"benchmark_patch_rmse_v", benchmark_patch_metrics.rmse_v, "V"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_patch_relative_rmse", benchmark_patch_relative_rmse, "ratio"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_patch_terminal_delta_v",
             benchmark_patch_metrics.terminal_potential_delta_v, "V"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_patch_compared_sample_count",
             static_cast<double>(benchmark_patch_metrics.compared_sample_count), "count"});
    }
    if (benchmark_body_metrics.valid)
    {
        benchmark_case.validation_metrics.push_back(
            {"benchmark_body_rmse_v", benchmark_body_metrics.rmse_v, "V"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_body_relative_rmse", benchmark_body_relative_rmse, "ratio"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_body_terminal_delta_v",
             benchmark_body_metrics.terminal_potential_delta_v, "V"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_body_compared_sample_count",
             static_cast<double>(benchmark_body_metrics.compared_sample_count), "count"});
    }
    benchmark_case.validation_metrics.insert(benchmark_case.validation_metrics.end(),
                                             benchmark_process_validation_metrics.begin(),
                                             benchmark_process_validation_metrics.end());
    if (benchmark_acceptance_gate.applicable)
    {
        benchmark_case.validation_metrics.push_back(
            {"benchmark_acceptance_gate_pass", benchmark_acceptance_gate.pass ? 1.0 : 0.0,
             "flag"});
        benchmark_case.validation_metrics.push_back(
            {"benchmark_acceptance_gate_checks_failed",
             static_cast<double>(benchmark_acceptance_gate.checks_failed), "count"});
    }
    if (benchmark_patch_metrics.valid || benchmark_body_metrics.valid)
    {
        benchmark_case.validation_metrics.push_back(
            {"benchmark_curve_relative_rmse_max", benchmark_curve_relative_rmse_max, "ratio"});
    }

    benchmark_case.tolerance_profile.relative_tolerance =
        benchmark_case.reference_family == "native" ? 0.0 : 0.10;
    benchmark_case.tolerance_profile.absolute_tolerance = benchmark_absolute_tolerance_v;
    benchmark_case.tolerance_profile.rmse_tolerance = benchmark_absolute_tolerance_v;
    benchmark_case.tolerance_profile.drift_tolerance = 5.0e-3;

    if (have_benchmark_definition)
    {
        const std::string baseline_source =
            legacyBenchmarkBaselineFamilyName(resolved_benchmark_source);
        const std::string baseline_origin =
            legacyBenchmarkBaselineOriginName(resolved_benchmark_source);
        const auto append_reference_dataset =
            [&](const std::string& id_suffix, const std::filesystem::path& path,
                const std::string& citation_suffix) {
                if (path.empty())
                {
                    return;
                }
                benchmark_case.reference_datasets.push_back(
                    {benchmark_case_id + "_" + id_suffix, baseline_source,
                     baseline_origin + ":" + citation_suffix, path.string(), true});
            };
        append_reference_dataset("patch_curve", benchmark_case_definition.paths.patch_reference_curve_path,
                                 "patch_curve");
        append_reference_dataset("body_curve", benchmark_case_definition.paths.body_reference_curve_path,
                                 "body_curve");
        append_reference_dataset("input_deck", benchmark_case_definition.paths.input_path,
                                 "input_deck");
        append_reference_dataset("environment_table", benchmark_case_definition.paths.environment_path,
                                 "environment_table");
        append_reference_dataset("structure_material_table",
                                 benchmark_case_definition.paths.structure_material_path,
                                 "structure_material_table");
        append_reference_dataset("dielectric_patch_material_table",
                                 benchmark_case_definition.paths.dielectric_patch_material_path,
                                 "dielectric_patch_material_table");
        append_reference_dataset("metal_patch_material_table",
                                 benchmark_case_definition.paths.metal_patch_material_path,
                                 "metal_patch_material_table");
    }

    std::string benchmark_case_error;
    if (!Coupling::Contracts::writeBenchmarkCaseJson(
            benchmark_case_path, benchmark_case, &benchmark_case_error))
    {
        last_error_message_ = benchmark_case_error.empty()
                                  ? "Failed to write benchmark case json to: " +
                                        benchmark_case_path.string()
                                  : benchmark_case_error;
        return false;
    }

    return true;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
