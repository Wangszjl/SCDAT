#include "DensePlasmaSurfaceCharging.h"
#include "LegacyBenchmarkSupport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <random>

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

double clampSigned(double value, double limit)
{
    return std::clamp(value, -limit, limit);
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

double halfNormalFluxAverage(double thermal_sigma, std::mt19937_64& generator)
{
    std::normal_distribution<double> distribution(0.0, thermal_sigma);
    return std::abs(distribution(generator));
}

double spectrumDensity(const Particle::ResolvedSpectrum& spectrum, double fallback_density)
{
    double density = 0.0;
    for (const auto& population : spectrum.populations)
    {
        density += std::max(0.0, population.density_m3);
    }
    return density > 0.0 ? density : std::max(0.0, fallback_density);
}

double spectrumTemperatureEv(const Particle::ResolvedSpectrum& spectrum, double fallback_temperature_ev)
{
    double weighted_temperature = 0.0;
    double total_density = 0.0;
    for (const auto& population : spectrum.populations)
    {
        weighted_temperature += population.density_m3 * population.temperature_ev;
        total_density += population.density_m3;
    }
    if (total_density > 0.0)
    {
        return weighted_temperature / total_density;
    }
    if (!spectrum.energy_grid_ev.empty() &&
        spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size())
    {
        double numerator = 0.0;
        double denominator = 0.0;
        for (size_t i = 1; i < spectrum.energy_grid_ev.size(); ++i)
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
            numerator += 0.5 * (e0 * f0 + e1 * f1) * width;
            denominator += 0.5 * (f0 + f1) * width;
        }
        if (denominator > 0.0)
        {
            return numerator / denominator;
        }
    }
    return std::max(1.0e-3, fallback_temperature_ev);
}

double spectrumAverageMassAmu(const Particle::ResolvedSpectrum& spectrum, double fallback_mass_amu)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (const auto& population : spectrum.populations)
    {
        numerator += std::max(0.0, population.density_m3) * std::max(1.0, population.mass_amu);
        denominator += std::max(0.0, population.density_m3);
    }
    return denominator > 0.0 ? numerator / denominator : std::max(1.0, fallback_mass_amu);
}

void synchronizePlasmaMomentsFromSpectra(SurfaceChargingConfig& config)
{
    if (config.has_electron_spectrum)
    {
        config.plasma.electron_density_m3 =
            spectrumDensity(config.electron_spectrum, config.plasma.electron_density_m3);
        config.plasma.electron_temperature_ev =
            spectrumTemperatureEv(config.electron_spectrum, config.plasma.electron_temperature_ev);
    }
    if (config.has_ion_spectrum)
    {
        config.plasma.ion_density_m3 =
            spectrumDensity(config.ion_spectrum, config.plasma.ion_density_m3);
        config.plasma.ion_temperature_ev =
            spectrumTemperatureEv(config.ion_spectrum, config.plasma.ion_temperature_ev);
        config.plasma.ion_mass_amu =
            spectrumAverageMassAmu(config.ion_spectrum, config.plasma.ion_mass_amu);
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
    return route == SurfaceRuntimeRoute::LegacyBenchmark ? 1.0 : 0.0;
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
    return route == SurfaceRuntimeRoute::LegacyBenchmark ? "LegacyBenchmark"
                                                         : "SCDATUnified";
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
        snapshot.average_node_degree =
            std::accumulate(node_degrees.begin(), node_degrees.end(), 0.0) /
            static_cast<double>(node_degrees.size());
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
        node_degrees.empty()
            ? 0.0
            : std::accumulate(node_degrees.begin(), node_degrees.end(), 0.0) /
                  static_cast<double>(node_degrees.size());

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

bool writeSurfaceBoundaryMappingJson(const std::filesystem::path& json_path,
                                     const SurfaceChargingConfig& config)
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
    output << "    \"field_request_json\": \"" << jsonEscape(artifact(".field_request.json").string()) << "\"\n";
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
    config_ = config;
    last_error_message_.clear();
    std::string validation_error;
    if (!validateSurfaceChargingConfig(config_, &validation_error))
    {
        initialized_ = false;
        last_error_message_ = validation_error;
        return false;
    }
    if (config_.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark &&
        config_.legacy_benchmark_execution_mode ==
            LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm)
    {
        config_ = applyLegacyBenchmarkExecutionConfig(config_);
    }
    config_ = normalizeSurfaceChargingConfig(config_);
    synchronizePlasmaMomentsFromSpectra(config_);
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
    history_normal_electric_field_.clear();
    history_local_charge_density_.clear();
    history_adaptive_time_step_.clear();
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
    history_surface_node_normal_electric_fields_.clear();
    history_surface_node_local_charge_densities_.clear();
    history_surface_node_propagated_reference_potentials_.clear();
    history_surface_node_field_solver_reference_potentials_.clear();
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
        history_surface_node_normal_electric_fields_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_local_charge_densities_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_propagated_reference_potentials_.assign(circuit_model_->nodeCount(), {});
        history_surface_node_field_solver_reference_potentials_.assign(circuit_model_->nodeCount(), {});
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
    reference_config.plasma = config_.plasma;
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

    const double derived_photo_current =
        kElementaryCharge * config_.emission.photon_flux_m2_s *
        std::clamp(config_.material.getScalarProperty("photoelectron_yield", 0.0) *
                       config_.emission.enhancement_factor,
                   0.0, 1.0);
    reference_config.body_photo_current_density_a_per_m2 =
        config_.body_photo_current_density_a_per_m2 > 0.0 ? config_.body_photo_current_density_a_per_m2
                                                          : derived_photo_current;
    reference_config.patch_photo_current_density_a_per_m2 =
        config_.patch_photo_current_density_a_per_m2 > 0.0 ? config_.patch_photo_current_density_a_per_m2
                                                           : derived_photo_current;

    reference_config.body_material = config_.material;
    reference_config.body_material.setType(Mesh::MaterialType::CONDUCTOR);
    reference_config.body_material.setName(config_.material.getName() + "_body");
    reference_config.body_material.setConductivity(
        std::max(1.0e-6, config_.material.getScalarProperty("body_conductivity_s_per_m", 1.0e4)));
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
    ensure_node_storage(history_surface_node_normal_electric_fields_);
    ensure_node_storage(history_surface_node_local_charge_densities_);
    ensure_node_storage(history_surface_node_propagated_reference_potentials_);
    ensure_node_storage(history_surface_node_field_solver_reference_potentials_);
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
        history_surface_node_normal_electric_fields_[node_index].push_back(
            runtime_state.normal_electric_field_v_per_m);
        history_surface_node_local_charge_densities_[node_index].push_back(
            runtime_state.local_charge_density_c_per_m3);
        history_surface_node_propagated_reference_potentials_[node_index].push_back(
            runtime_state.propagated_reference_potential_v);
        history_surface_node_field_solver_reference_potentials_[node_index].push_back(
            runtime_state.field_solver_reference_potential_v);
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
    }
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
        const double eps_r = std::max(1.0, config_.material.getPermittivity());
        if (config_.derive_capacitance_from_material)
        {
            base_capacitance = kEpsilon0 * eps_r / thickness;
        }
        else
        {
            base_capacitance = std::max(1.0e-12, config_.capacitance_per_area_f_per_m2);
        }
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

double DensePlasmaSurfaceCharging::computeEffectiveSheathLength() const
{
    if (!config_.derive_sheath_length_from_plasma)
    {
        return std::max(1.0e-6, config_.sheath_length_m);
    }

    const double electron_temperature_j =
        std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_density = std::max(1.0e3, config_.plasma.electron_density_m3);
    const double debye_length = std::sqrt(kEpsilon0 * electron_temperature_j /
                                          (electron_density * kElementaryCharge * kElementaryCharge));
    return std::clamp(debye_length, std::max(1.0e-6, config_.minimum_sheath_length_m),
                      std::max(config_.minimum_sheath_length_m, config_.maximum_sheath_length_m));
}

SurfaceCurrents DensePlasmaSurfaceCharging::computeLegacySurfaceCurrents(
    double surface_potential_v) const
{
    SurfaceCurrents currents;
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
    const double flow_alignment = std::clamp(config_.flow_alignment_cosine, -1.0, 1.0);
    const double normal_flow_speed =
        std::max(0.0, config_.bulk_flow_velocity_m_per_s * flow_alignment);
    const double directed_ion_velocity =
        std::max(0.0, config_.ion_directed_velocity_m_per_s + normal_flow_speed);
    const double effective_ion_velocity =
        std::sqrt(bohm_velocity * bohm_velocity + directed_ion_velocity * directed_ion_velocity);
    const double phi_scale = std::max(0.25, config_.plasma.electron_temperature_ev);
    double electron_flux =
        kElementaryCharge * config_.plasma.electron_density_m3 * electron_thermal_velocity *
        std::max(0.0, config_.electron_collection_coefficient);
    double ion_flux = kElementaryCharge * config_.plasma.ion_density_m3 * effective_ion_velocity *
                      std::max(0.0, config_.ion_collection_coefficient);

    if (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma)
    {
        const double ion_advective_flux =
            kElementaryCharge * config_.plasma.ion_density_m3 * normal_flow_speed *
            std::max(0.0, config_.ion_collection_coefficient);
        const double electron_advective_flux =
            kElementaryCharge * config_.plasma.electron_density_m3 * normal_flow_speed *
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
        emittedElectronEscapeProbability(surface_potential_v,
                                         config_.material.getScalarProperty(
                                             "secondary_emission_escape_energy_ev", 2.0));
    const double photo_escape_probability =
        emittedElectronEscapeProbability(surface_potential_v,
                                         config_.material.getScalarProperty(
                                             "photoelectron_escape_energy_ev", 1.5));
    const double thermionic_escape_probability =
        emittedElectronEscapeProbability(surface_potential_v,
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
        kElementaryCharge * config_.emission.photon_flux_m2_s *
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
    return current_model_->evaluate(ReferenceSurfaceRole::Patch, runtime_state);
}

double DensePlasmaSurfaceCharging::computeRadiationInducedConductivity() const
{
    if (config_.radiation_conductivity_coefficient <= 0.0)
    {
        return 0.0;
    }

    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double electron_flux =
        kElementaryCharge * config_.plasma.electron_density_m3 * electron_thermal_velocity;
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
            kElementaryCharge * config_.plasma.electron_density_m3 * electron_thermal_velocity *
            std::max(0.0, config_.electron_collection_coefficient);
        const double analytic_ion_flux =
            kElementaryCharge * config_.plasma.ion_density_m3 * bohm_velocity *
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
    std::mt19937_64 generator(0x5cdab51ULL + static_cast<std::uint64_t>(sample_count));

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
    return 0.5 * kElementaryCharge * config_.plasma.electron_density_m3 * mean_surface_velocity *
           std::max(0.0, config_.electron_collection_coefficient);
}

double DensePlasmaSurfaceCharging::estimateIonFluxPicLike(double surface_potential_v,
                                                          std::size_t samples) const
{
    const std::size_t sample_count = std::max<std::size_t>(64, samples);
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double ti_j = std::max(1.0e-3, config_.plasma.ion_temperature_ev) * kElementaryCharge;
    const double sigma = std::sqrt(ti_j / ion_mass);
    const double drift_velocity =
        std::max(0.0, config_.ion_directed_velocity_m_per_s +
                          std::max(0.0, config_.bulk_flow_velocity_m_per_s *
                                            std::clamp(config_.flow_alignment_cosine, -1.0, 1.0)));
    std::mt19937_64 generator(0x71015aULL + static_cast<std::uint64_t>(sample_count));

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
    return 0.5 * kElementaryCharge * config_.plasma.ion_density_m3 * mean_surface_velocity *
           std::max(0.0, config_.ion_collection_coefficient);
}

double DensePlasmaSurfaceCharging::computeEffectiveConductivity(double surface_potential_v) const
{
    return computeEffectiveConductivity(surface_potential_v, config_.material);
}

double DensePlasmaSurfaceCharging::computeEffectiveConductivity(
    double surface_potential_v, const Material::MaterialProperty& material) const
{
    const double base_conductivity = std::max(0.0, material.getConductivity());
    const double thickness = std::max(1.0e-8, config_.dielectric_thickness_m);
    const std::size_t patch_node_index =
        circuit_model_ ? circuit_model_->primaryPatchNodeIndex() : 0;
    const std::string patch_node_name =
        circuit_model_ ? circuit_model_->nodeName(patch_node_index) : std::string{"patch"};
    const auto runtime_state =
        buildRuntimeState(status_.body_potential_v, surface_potential_v, config_.surface_area_m2,
                          base_conductivity, computeEffectiveSheathLength(), patch_node_index,
                          patch_node_name);
    const double electric_field = std::max(
        std::abs(runtime_state.normal_electric_field_v_per_m),
        std::abs(surface_potential_v - status_.body_potential_v) / thickness);
    const double poole_frenkel_beta =
        std::max(0.0, material.getScalarProperty("poole_frenkel_beta", 0.0));
    const double max_enhancement =
        std::max(1.0, material.getScalarProperty("max_field_enhancement_factor", 1.0e8));
    const double field_enhancement =
        std::min(max_enhancement, safeExp(poole_frenkel_beta * std::sqrt(std::max(0.0, electric_field))));
    const double space_charge_ratio =
        std::abs(runtime_state.local_charge_density_c_per_m3) * thickness /
        std::max(1.0e-12, kEpsilon0 * std::max(1.0, electric_field));
    const double space_charge_enhancement = 1.0 + std::min(5.0, space_charge_ratio);
    return base_conductivity * field_enhancement * space_charge_enhancement +
           computeRadiationInducedConductivity();
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
              graph_capacitance_matrix_provider_->graphCouplingMetric(config_, circuit_model_.get());
          const double coupling_penalty =
              std::clamp(graph_coupling_metric * 1.0e12, 0.0, 4.0);
          suggested_dt /= (1.0 + 0.5 * coupling_penalty);
      }

      const SurfaceModelRuntimeState runtime_state =
          buildRuntimeState(status_.body_potential_v, status_.patch_potential_v, config_.surface_area_m2,
                            0.0, computeEffectiveSheathLength());
      const double volume_penalty =
          1.0 + 0.35 * std::clamp(runtime_state.volume_mesh_coupling_gain, 0.0, 1.0) +
          0.05 * std::max(0.0, runtime_state.volume_projection_weight_sum - 1.0);
      suggested_dt /= volume_penalty;

      return std::clamp(std::min(remaining, suggested_dt), min_dt, max_dt);
  }

double DensePlasmaSurfaceCharging::advancePotentialImplicit(double surface_potential_v, double dt) const
{
    const double capacitance = std::max(1.0e-12, status_.state.capacitance_per_area_f_per_m2);
    const double search_limit =
        std::max(config_.max_abs_potential_v, 20.0 * std::max(1.0, config_.plasma.electron_temperature_ev));

    double potential = surface_potential_v;
    for (int iteration = 0; iteration < 20; ++iteration)
    {
        const double residual =
            capacitance * (potential - surface_potential_v) / dt - computeNetCurrentDensity(potential);
        if (std::abs(residual) < 1.0e-9)
        {
            return potential;
        }

        const double derivative =
            capacitance / dt - estimateCurrentDerivative(potential);
        if (std::abs(derivative) < 1.0e-12)
        {
            break;
        }

        const double delta =
            std::clamp(-residual / derivative, -0.2 * search_limit, 0.2 * search_limit);
        potential += delta;
        potential = clampSigned(potential, search_limit);
        if (std::abs(delta) < 1.0e-6)
        {
            return potential;
        }
    }

    const double explicit_candidate =
        surface_potential_v + dt * computeNetCurrentDensity(surface_potential_v) / capacitance;
    return clampSigned(explicit_candidate, search_limit);
}

bool DensePlasmaSurfaceCharging::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }
    if (!legacy_benchmark_replay_.empty())
    {
        return advanceLegacyBenchmarkReplay(dt);
    }

    const bool use_reference_circuit =
        config_.use_reference_current_balance &&
        (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
         config_.regime == SurfaceChargingRegime::GeoKineticPicLike) &&
        config_.enable_body_patch_circuit;

    if (use_reference_circuit)
    {
        bool recalibrated = false;
        if (config_.enable_pic_calibration &&
            (status_.steps_completed == 0 ||
             (config_.pic_recalibration_interval_steps > 0 &&
              status_.steps_completed % config_.pic_recalibration_interval_steps == 0) ||
             std::abs(status_.equilibrium_error) > config_.pic_recalibration_trigger_v))
        {
            applyGeoPicCalibration();
            recalibrated = true;
        }

        const std::size_t substeps = std::max<std::size_t>(1, config_.internal_substeps);
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

            Coupling::SurfaceCircuitLinearization linearization;
            linearization.plasma_current_a = std::vector<double>(circuit_model_->nodeCount(), 0.0);
            linearization.plasma_didv_a_per_v =
                std::vector<double>(circuit_model_->nodeCount(), 0.0);
            linearization.additional_rhs_a =
                std::vector<double>(circuit_model_->nodeCount(), 0.0);
            linearization.additional_diagonal_a_per_v =
                std::vector<double>(circuit_model_->nodeCount(), 0.0);

            for (std::size_t patch_ordinal = 0; patch_ordinal < circuit_model_->patchCount(); ++patch_ordinal)
            {
                const auto patch_node_index = circuit_model_->patchNodeIndex(patch_ordinal);
                const double local_patch_potential = circuit_model_->nodePotential(patch_node_index);
                const double patch_area_m2 =
                    std::max(1.0e-16, circuit_model_->nodeAreaM2(patch_node_index));
                const double local_effective_conductivity =
                    computeEffectiveConductivity(
                        local_patch_potential,
                        resolvePatchMaterial(config_, patch_node_index,
                                             circuit_model_->nodeName(patch_node_index)));
                const auto patch_branch_index = circuit_model_->patchToBodyBranchIndex(patch_ordinal);
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
                linearization.plasma_current_a[patch_node_index] =
                    (patch_currents_old.total_current_a_per_m2 -
                     patch_currents_old.conduction_current_a_per_m2) *
                    patch_area_m2;
                linearization.plasma_didv_a_per_v[patch_node_index] =
                    (patch_currents_old.current_derivative_a_per_m2_per_v +
                     local_effective_conductivity /
                         std::max(1.0e-9, config_.dielectric_thickness_m)) *
                    patch_area_m2;
            }

            if (config_.body_floating)
            {
                const auto runtime_state =
                    buildRuntimeState(body_potential, patch_potential,
                                      std::max(1.0e-16, circuit_model_->bodyAreaM2()),
                                      effective_conductivity, effective_sheath_length,
                                      circuit_model_->bodyNodeIndex(),
                                      circuit_model_->nodeName(circuit_model_->bodyNodeIndex()));
                const auto body_terms_old =
                    current_model_->evaluate(ReferenceSurfaceRole::Body, runtime_state);
                linearization.plasma_current_a[circuit_model_->bodyNodeIndex()] =
                    body_terms_old.total_current_a_per_m2 *
                    std::max(1.0e-16, circuit_model_->bodyAreaM2());
                linearization.plasma_didv_a_per_v[circuit_model_->bodyNodeIndex()] =
                    current_model_->computeCurrentDerivative(ReferenceSurfaceRole::Body,
                                                             runtime_state) *
                    std::max(1.0e-16, circuit_model_->bodyAreaM2());
            }

            if (graph_capacitance_matrix_provider_ != nullptr && sub_dt > 0.0)
            {
                const double runtime_matrix_weight = std::clamp(
                    config_.material.getScalarProperty("graph_matrix_runtime_weight", 0.35), 0.0, 2.0);
                if (runtime_matrix_weight > 0.0)
                {
                    linearization.additional_off_diagonal_entries.clear();
                    linearization.additional_off_diagonal_entries.reserve(
                        2 * circuit_model_->branchCount());
                    for (std::size_t branch_index = 0; branch_index < circuit_model_->branchCount();
                         ++branch_index)
                    {
                        const auto from_node = circuit_model_->branchFromNodeIndex(branch_index);
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

                        const double from_potential_old_v = circuit_model_->nodePotential(from_node);
                        const double to_potential_old_v = circuit_model_->nodePotential(to_node);
                        const double old_delta_v = from_potential_old_v - to_potential_old_v;

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

            const auto circuit_result =
                circuit_model_->advanceImplicit(sub_dt, linearization,
                                                std::max(1.0e-3, config_.max_delta_potential_v_per_step));
            if (!circuit_result.converged)
            {
                return false;
            }

            body_potential = circuit_result.node_potentials_v[circuit_model_->bodyNodeIndex()];
            patch_potential =
                circuit_result.node_potentials_v[circuit_model_->primaryPatchNodeIndex()];
            circuit_node_potentials_ = circuit_result.node_potentials_v;
            if (history_surface_branch_currents_.size() != circuit_result.branch_currents_a.size())
            {
                history_surface_branch_currents_.assign(circuit_result.branch_currents_a.size(), {});
            }
            latest_branch_currents_a = circuit_result.branch_currents_a;
            if (circuit_model_->primaryPatchToBodyBranchIndex() <
                circuit_result.branch_currents_a.size())
            {
                branch_current_a = circuit_result.branch_currents_a
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
        history_ion_secondary_current_.push_back(status_.currents.ion_secondary_emission_a_per_m2);
        history_backscatter_current_.push_back(status_.currents.backscatter_emission_a_per_m2);
        history_photo_current_.push_back(status_.currents.photo_emission_a_per_m2);
        history_thermionic_current_.push_back(status_.currents.thermionic_emission_a_per_m2);
        history_field_emission_current_.push_back(status_.currents.field_emission_a_per_m2);
        history_leakage_current_.push_back(leakage_current);
        history_ram_current_.push_back(status_.currents.ram_ion_current_a_per_m2);
        history_body_potential_.push_back(status_.body_potential_v);
        history_patch_potential_.push_back(status_.patch_potential_v);
        history_circuit_branch_current_.push_back(status_.circuit_branch_current_a);
        history_current_derivative_.push_back(status_.currents.current_derivative_a_per_m2_per_v);
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
        history_normal_electric_field_.push_back(final_runtime_state.normal_electric_field_v_per_m);
        history_local_charge_density_.push_back(final_runtime_state.local_charge_density_c_per_m3);
        history_adaptive_time_step_.push_back(dt);
        history_electron_calibration_factor_.push_back(
            current_model_ ? current_model_->electronCalibrationFactor() : 1.0);
        history_ion_calibration_factor_.push_back(
            current_model_ ? current_model_->ionCalibrationFactor() : 1.0);
        history_equilibrium_potential_.push_back(target_potential);
        history_equilibrium_error_.push_back(status_.equilibrium_error);
        appendSurfaceNodeDiagnostics(latest_branch_currents_a, effective_sheath_length);
        return true;
    }

    const std::size_t substeps = std::max<std::size_t>(1, config_.internal_substeps);
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
            const double updated_potential = advancePotentialImplicit(state.surface_potential_v, sub_dt);
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
    history_electron_calibration_factor_.push_back(
        current_model_ ? current_model_->electronCalibrationFactor() : 1.0);
    history_ion_calibration_factor_.push_back(
        current_model_ ? current_model_->ionCalibrationFactor() : 1.0);
    history_equilibrium_potential_.push_back(target_potential);
    history_equilibrium_error_.push_back(status_.equilibrium_error);
    std::vector<double> branch_currents_a(circuit_model_ ? circuit_model_->branchCount() : 1, 0.0);
    appendSurfaceNodeDiagnostics(branch_currents_a, effective_sheath_length);
    return true;
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
    data_set.scalar_series["normal_electric_field_v_per_m"] = history_normal_electric_field_;
    data_set.scalar_series["local_charge_density_c_per_m3"] = history_local_charge_density_;
    data_set.scalar_series["adaptive_time_step_s"] = history_adaptive_time_step_;
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
        export_node_series("normal_electric_field_v_per_m",
                           history_surface_node_normal_electric_fields_);
        export_node_series("local_charge_density_c_per_m3",
                           history_surface_node_local_charge_densities_);
        export_node_series("propagated_reference_potential_v",
                           history_surface_node_propagated_reference_potentials_);
        export_node_series("field_solver_reference_potential_v",
                           history_surface_node_field_solver_reference_potentials_);
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
            : "UnifiedSpisAligned";
    data_set.metadata["runtime_route"] = runtimeRouteName(config_.runtime_route);
    data_set.metadata["benchmark_source"] = benchmarkSourceName(config_.benchmark_source);
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
    LegacyBenchmarkConsistencyDiagnostics benchmark_consistency{};
    if (config_.benchmark_source != SurfaceBenchmarkSource::None)
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
            for (std::size_t i = 0; i < history_time_.size(); ++i)
            {
                LegacyBenchmarkCurveSample body_sample;
                body_sample.cycle_index = i;
                body_sample.time_s = history_time_[i];
                body_sample.potential_v =
                    i < history_body_potential_.size() ? history_body_potential_[i] : 0.0;
                actual_body_curve.push_back(body_sample);
            }
        }

        benchmark_case_definition = loadLegacyBenchmarkCaseDefinition(config_.benchmark_source);
        have_benchmark_definition = true;
        const bool is_c_leo_source =
            config_.benchmark_source == SurfaceBenchmarkSource::CLeoRam ||
            config_.benchmark_source == SurfaceBenchmarkSource::CLeoWake;
        const bool execute_legacy_algorithm =
            config_.legacy_benchmark_execution_mode ==
            LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
        if (is_c_leo_source && execute_legacy_algorithm && !actual_body_curve.empty())
        {
            const auto consistency =
                analyzeLegacyBenchmarkConsistency(config_.benchmark_source, benchmark_case_definition);
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
            const bool is_geo = config_.benchmark_source == SurfaceBenchmarkSource::CGeo;
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
        benchmark_acceptance = legacyBenchmarkAcceptanceCriteria(config_.benchmark_source);
        const auto body_component_rmse = [&](const auto& accessor) {
            const std::size_t count =
                std::min(actual_body_curve.size(), benchmark_case_definition.body_reference_curve.size());
            if (count == 0)
            {
                return 0.0;
            }

            double error_sum = 0.0;
            for (std::size_t i = 0; i < count; ++i)
            {
                const double delta = accessor(actual_body_curve[i]) -
                                     accessor(benchmark_case_definition.body_reference_curve[i]);
                error_sum += delta * delta;
            }
            return std::sqrt(error_sum / static_cast<double>(count));
        };
        const double body_je_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; });
        const double body_jnet_rmse_a_per_m2 = body_component_rmse(
            [](const LegacyBenchmarkCurveSample& sample) { return sample.jnet_a_per_m2; });
        benchmark_acceptance_gate =
            evaluateLegacyBenchmarkAcceptanceGate(benchmark_acceptance,
                                                  benchmark_patch_metrics,
                                                  benchmark_body_metrics,
                                                  body_je_rmse_a_per_m2,
                                                  body_jnet_rmse_a_per_m2,
                                                  benchmark_body_tail_metrics);
        benchmark_consistency =
            analyzeLegacyBenchmarkConsistency(config_.benchmark_source, benchmark_case_definition);
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
            legacyBenchmarkBaselineFamilyName(config_.benchmark_source);
        data_set.metadata["benchmark_baseline_origin"] =
            legacyBenchmarkBaselineOriginName(config_.benchmark_source);
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
        switch (config_.benchmark_source)
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
            : "UnifiedSpisAligned";
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

    const bool export_ok = static_cast<bool>(exporter_.exportDataSet(csv_path, data_set));
    if (!export_ok)
    {
        last_error_message_ = "Failed to export surface charging results to: " + csv_path.string();
        return false;
    }
    if (have_benchmark_definition)
    {
        auto report_path = csv_path;
        report_path.replace_extension(".benchmark.txt");
        if (!writeLegacyBenchmarkReport(report_path, config_.runtime_route, config_.benchmark_source,
                                        config_.legacy_benchmark_execution_mode,
                                        benchmark_case_definition, actual_curve, actual_body_curve,
                                        benchmark_patch_metrics, benchmark_body_metrics,
                                        benchmark_consistency))
        {
            last_error_message_ =
                "Failed to write legacy benchmark sidecar report to: " + report_path.string();
            return false;
        }

        auto comparison_csv_path = csv_path;
        comparison_csv_path.replace_extension(".benchmark.csv");
        if (!writeLegacyBenchmarkComparisonCsv(comparison_csv_path, actual_curve,
                                               benchmark_case_definition.patch_reference_curve,
                                               actual_body_curve,
                                               benchmark_case_definition.body_reference_curve))
        {
            last_error_message_ =
                "Failed to write legacy benchmark comparison csv to: " +
                comparison_csv_path.string();
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
        if (!writeSurfaceBoundaryMappingJson(boundary_mapping_path, config_))
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
    last_error_message_.clear();
    return true;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
