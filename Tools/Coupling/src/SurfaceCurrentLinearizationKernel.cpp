#include "../include/SurfaceCurrentLinearizationKernel.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SCDAT
{
namespace Coupling
{
namespace
{
constexpr double kEpsilon0 = 8.8541878128e-12;
}

SurfaceCurrentLinearizationResult SurfaceCurrentLinearizationKernel::assemble(
    const std::vector<SurfaceCurrentLinearizationNode>& nodes,
    const std::vector<SurfaceCircuitLinearization::OffDiagonalEntry>& off_diagonal_entries) const
{
    CircuitDidvCompositionInput didv_input;
    didv_input.nodes.resize(nodes.size());
    didv_input.off_diagonal_entries = off_diagonal_entries;

    for (std::size_t index = 0; index < nodes.size(); ++index)
    {
        const auto& node = nodes[index];
        const auto& components = node.components;
        auto& didv_node = didv_input.nodes[index];
        didv_node.terms = {
            CircuitDidvTerm{"electron_collection", components.electron_collection_a,
                            components.electron_collection_didv_a_per_v},
            CircuitDidvTerm{"ion_collection", components.ion_collection_a,
                            components.ion_collection_didv_a_per_v},
            CircuitDidvTerm{"secondary_electron_emission",
                            -components.secondary_electron_emission_a,
                            -components.secondary_electron_emission_didv_a_per_v},
            CircuitDidvTerm{"ion_secondary_electron_emission",
                            -components.ion_secondary_electron_emission_a,
                            -components.ion_secondary_electron_emission_didv_a_per_v},
            CircuitDidvTerm{"backscatter_emission", -components.backscatter_emission_a,
                            -components.backscatter_emission_didv_a_per_v},
            CircuitDidvTerm{"photoelectron_emission", -components.photoelectron_emission_a,
                            -components.photoelectron_emission_didv_a_per_v},
            CircuitDidvTerm{"conduction", components.conduction_a,
                            components.conduction_didv_a_per_v},
        };
        didv_node.additional_rhs_a = node.additional_rhs_a;
        didv_node.additional_diagonal_a_per_v = node.additional_diagonal_a_per_v;
    }

    const auto didv_result = composeCircuitDidv(didv_input);

    SurfaceCurrentLinearizationResult result;
    result.total_current_a.resize(nodes.size(), 0.0);
    result.total_didv_a_per_v.resize(nodes.size(), 0.0);
    result.component_states.reserve(nodes.size());
    result.kernel_input = didv_result.kernel_input;

    for (std::size_t index = 0; index < nodes.size(); ++index)
    {
        const auto& components = nodes[index].components;
        result.total_current_a[index] = didv_result.node_summaries[index].total_current_a;
        result.total_didv_a_per_v[index] = didv_result.node_summaries[index].total_didv_a_per_v;
        result.component_states.push_back(components);
    }

    return result;
}

double estimateCenteredDidv(const std::function<double(double)>& evaluator, double sample_x,
                            double minimum_step, double relative_step)
{
    const double step =
        std::max(minimum_step, relative_step * std::max(1.0, std::abs(sample_x)));
    const double positive = evaluator(sample_x + step);
    const double negative = evaluator(sample_x - step);
    if (!std::isfinite(positive) || !std::isfinite(negative))
    {
        return 0.0;
    }
    return (positive - negative) / (2.0 * step);
}

SurfaceCurrentRootSolveResult solveBisectionCurrentRoot(
    const std::function<double(double)>& evaluator, double minimum_x, double maximum_x,
    double function_tolerance, double bracket_tolerance, int max_iterations)
{
    SurfaceCurrentRootSolveResult result;
    double lo = minimum_x;
    double hi = maximum_x;
    double flo = evaluator(lo);
    double fhi = evaluator(hi);
    if (!std::isfinite(flo) || !std::isfinite(fhi))
    {
        result.root_x = std::clamp(0.0, minimum_x, maximum_x);
        return result;
    }
    if (flo == 0.0)
    {
        result.root_x = lo;
        result.value_at_root = flo;
        result.converged = true;
        return result;
    }
    if (fhi == 0.0)
    {
        result.root_x = hi;
        result.value_at_root = fhi;
        result.converged = true;
        return result;
    }

    if (flo * fhi > 0.0)
    {
        const double mid = 0.0;
        const double fmid = evaluator(mid);
        if (std::isfinite(fmid) && flo * fmid <= 0.0)
        {
            hi = mid;
            fhi = fmid;
        }
        else if (std::isfinite(fmid) && fmid * fhi <= 0.0)
        {
            lo = mid;
            flo = fmid;
        }
        else
        {
            result.root_x = std::clamp(mid, minimum_x, maximum_x);
            result.value_at_root = std::isfinite(fmid) ? fmid : 0.0;
            return result;
        }
    }

    for (int iteration = 0; iteration < max_iterations; ++iteration)
    {
        const double mid = 0.5 * (lo + hi);
        const double fmid = evaluator(mid);
        if (!std::isfinite(fmid))
        {
            break;
        }
        if (std::abs(fmid) <= function_tolerance || std::abs(hi - lo) <= bracket_tolerance)
        {
            result.root_x = std::clamp(mid, minimum_x, maximum_x);
            result.value_at_root = fmid;
            result.converged = true;
            return result;
        }
        if (flo * fmid <= 0.0)
        {
            hi = mid;
            fhi = fmid;
        }
        else
        {
            lo = mid;
            flo = fmid;
        }
    }

    result.root_x = std::clamp(0.5 * (lo + hi), minimum_x, maximum_x);
    result.value_at_root = evaluator(result.root_x);
    result.converged = std::isfinite(result.value_at_root);
    return result;
}

double aggregateSurfaceCurrentDensity(const SurfaceCurrentComponentState& components,
                                      double thermionic_emission_a,
                                      double field_emission_a,
                                      double ram_ion_current_a)
{
    return components.electron_collection_a + components.ion_collection_a -
           components.secondary_electron_emission_a -
           components.ion_secondary_electron_emission_a -
           components.backscatter_emission_a - components.photoelectron_emission_a +
           components.conduction_a + thermionic_emission_a + field_emission_a +
           ram_ion_current_a;
}

SurfaceCurrentDensityLinearizationResult linearizeSurfaceCurrentDensity(
    const SurfaceCurrentDensityLinearizationInputs& inputs)
{
    SurfaceCurrentDensityLinearizationResult result;
    const double area_m2 = std::max(1.0e-12, inputs.area_m2);

    SurfaceCurrentLinearizationNode node;
    node.components.electron_collection_a =
        inputs.component_densities.electron_collection_a * area_m2;
    node.components.ion_collection_a =
        inputs.component_densities.ion_collection_a * area_m2;
    node.components.secondary_electron_emission_a =
        inputs.component_densities.secondary_electron_emission_a * area_m2;
    node.components.ion_secondary_electron_emission_a =
        inputs.component_densities.ion_secondary_electron_emission_a * area_m2;
    node.components.backscatter_emission_a =
        inputs.component_densities.backscatter_emission_a * area_m2;
    node.components.photoelectron_emission_a =
        inputs.component_densities.photoelectron_emission_a * area_m2;
    node.components.conduction_a = inputs.component_densities.conduction_a * area_m2;
    node.components.electron_collection_didv_a_per_v =
        inputs.plasma_didv_density_a_per_m2_per_v * area_m2;
    node.components.conduction_didv_a_per_v =
        inputs.conduction_didv_density_a_per_m2_per_v * area_m2;

    SurfaceCurrentLinearizationKernel kernel;
    result.linearization = kernel.assemble({node});

    result.plasma_current_a =
        aggregateSurfaceCurrentDensity(inputs.component_densities) * area_m2;
    result.plasma_didv_a_per_v = inputs.plasma_didv_density_a_per_m2_per_v * area_m2;
    result.conduction_current_a = inputs.component_densities.conduction_a * area_m2;
    result.conduction_didv_a_per_v = inputs.conduction_didv_density_a_per_m2_per_v * area_m2;
    result.total_current_density_a_per_m2 = aggregateSurfaceCurrentDensity(
        inputs.component_densities, inputs.thermionic_emission_density_a_per_m2,
        inputs.field_emission_density_a_per_m2, inputs.ram_ion_current_density_a_per_m2);
    result.total_current_a = result.total_current_density_a_per_m2 * area_m2;
    result.total_didv_a_per_v = result.linearization.total_didv_a_per_v.empty()
                                    ? (inputs.plasma_didv_density_a_per_m2_per_v +
                                       inputs.conduction_didv_density_a_per_m2_per_v) *
                                          area_m2
                                    : result.linearization.total_didv_a_per_v.front();
    result.total_didv_density_a_per_m2_per_v = result.total_didv_a_per_v / area_m2;
    return result;
}

SurfaceCalibrationFactors estimateSurfaceCalibrationFactors(
    const SurfaceCalibrationFactorInputs& inputs)
{
    SurfaceCalibrationFactors factors;
    factors.electron_factor = estimateCalibrationFactor(
        inputs.measured_electron_current_a, inputs.native_electron_current_a,
        inputs.minimum_factor, inputs.maximum_factor, inputs.default_electron_factor);
    factors.ion_factor = estimateCalibrationFactor(
        inputs.measured_ion_current_a, inputs.native_ion_current_a, inputs.minimum_factor,
        inputs.maximum_factor, inputs.default_ion_factor);
    return factors;
}

SurfaceEmissionBlendResult blendSurfaceEmissionComponents(
    const SurfaceEmissionBlendInputs& inputs)
{
    SurfaceEmissionBlendResult result;
    result.reference_emitted_total_a = inputs.reference_secondary_emission_a +
                                       inputs.reference_ion_secondary_emission_a +
                                       inputs.reference_backscatter_emission_a;

    if (!inputs.hybrid_completion)
    {
        result.secondary_emission_a = inputs.measured_emitted_total_a;
        result.ion_secondary_emission_a = 0.0;
        result.backscatter_emission_a = 0.0;
        result.emission_scale = std::abs(result.reference_emitted_total_a) > 1.0e-18
                                    ? inputs.measured_emitted_total_a /
                                          result.reference_emitted_total_a
                                    : 1.0;
        return result;
    }

    if (std::abs(result.reference_emitted_total_a) > 1.0e-18)
    {
        result.emission_scale = inputs.measured_emitted_total_a / result.reference_emitted_total_a;
        result.secondary_emission_a =
            inputs.reference_secondary_emission_a * result.emission_scale;
        result.ion_secondary_emission_a =
            inputs.reference_ion_secondary_emission_a * result.emission_scale;
        result.backscatter_emission_a =
            inputs.reference_backscatter_emission_a * result.emission_scale;
        return result;
    }

    result.secondary_emission_a = inputs.measured_emitted_total_a;
    result.ion_secondary_emission_a = 0.0;
    result.backscatter_emission_a = 0.0;
    result.emission_scale = 1.0;
    return result;
}

SurfaceDistributionBlendResult blendSurfaceDistributionComponents(
    const SurfaceDistributionBlendInputs& inputs)
{
    SurfaceDistributionBlendResult result;
    result.electron_collection_scale =
        std::abs(inputs.reference_electron_current_a) > 1.0e-18
            ? inputs.measured_electron_current_a / inputs.reference_electron_current_a
            : 1.0;
    result.ion_collection_scale =
        std::abs(inputs.reference_ion_current_a) > 1.0e-18
            ? inputs.measured_ion_current_a / inputs.reference_ion_current_a
            : 1.0;
    result.electron_characteristic_energy_ev =
        std::max(1.0e-6, inputs.electron_characteristic_energy_ev);
    result.ion_characteristic_energy_ev =
        std::max(1.0e-6, inputs.ion_characteristic_energy_ev);
    return result;
}

SurfaceSharedCurrentMatrixFallbackResult computeSurfaceSharedCurrentMatrixFallback(
    bool fallback_shared_current_matrix_coupling_enabled,
    bool patch_didv_valid,
    std::size_t patch_node_index,
    double patch_plasma_didv_a_per_v,
    double patch_conduction_didv_a_per_v,
    const std::vector<double>& shared_patch_voltage_weights,
    const std::vector<std::size_t>& patch_node_indices,
    const std::vector<std::size_t>& circuit_to_reduced_node_index,
    double minimum_coupled_derivative_abs_a_per_v)
{
    SurfaceSharedCurrentMatrixFallbackResult result;

    const double local_shared_weight =
        fallback_shared_current_matrix_coupling_enabled &&
                patch_node_index < shared_patch_voltage_weights.size()
            ? shared_patch_voltage_weights[patch_node_index]
            : 1.0;
    result.plasma_didv_a_per_v =
        patch_didv_valid
            ? (patch_plasma_didv_a_per_v * local_shared_weight +
               patch_conduction_didv_a_per_v)
            : 0.0;

    if (!fallback_shared_current_matrix_coupling_enabled || !patch_didv_valid)
    {
        return result;
    }

    const auto reduced_node_index =
        patch_node_index < circuit_to_reduced_node_index.size()
            ? circuit_to_reduced_node_index[patch_node_index]
            : std::numeric_limits<std::size_t>::max();

    for (const auto other_patch_node_index : patch_node_indices)
    {
        if (other_patch_node_index == patch_node_index)
        {
            continue;
        }

        const auto other_reduced_node_index =
            other_patch_node_index < circuit_to_reduced_node_index.size()
                ? circuit_to_reduced_node_index[other_patch_node_index]
                : std::numeric_limits<std::size_t>::max();
        if (reduced_node_index == other_reduced_node_index)
        {
            continue;
        }
        if (other_patch_node_index >= shared_patch_voltage_weights.size())
        {
            continue;
        }

        const double coupled_derivative_a_per_v =
            patch_plasma_didv_a_per_v *
            shared_patch_voltage_weights[other_patch_node_index];
        if (!std::isfinite(coupled_derivative_a_per_v) ||
            std::abs(coupled_derivative_a_per_v) <=
                std::max(0.0, minimum_coupled_derivative_abs_a_per_v))
        {
            continue;
        }

        SurfaceCircuitLinearization::OffDiagonalEntry entry{};
        entry.row_node = patch_node_index;
        entry.column_node = other_patch_node_index;
        entry.coefficient_a_per_v = -coupled_derivative_a_per_v;
        result.off_diagonal_entries.push_back(entry);
    }

    return result;
}

SurfacePicRouteAdjustmentResult applySurfacePicRouteAdjustment(
    const SurfacePicRouteAdjustmentInputs& inputs,
    const SurfaceCurrentComponentState& current_components,
    double current_total_didv_a_per_v)
{
    SurfacePicRouteAdjustmentResult result;
    result.components = current_components;
    result.total_didv_a_per_v = current_total_didv_a_per_v;

    if (!inputs.formal_surface_pic_route || !inputs.kernel_snapshot_valid)
    {
        return result;
    }

    result.components = inputs.snapshot_components;
    result.total_didv_a_per_v = inputs.snapshot_total_didv_a_per_v;
    result.applied = true;
    return result;
}

std::size_t computeSurfaceAdaptiveSubstepCount(
    const SurfaceAdaptiveSubstepInputs& inputs)
{
    const std::size_t base_substeps = std::max<std::size_t>(1, inputs.base_substeps);
    const std::size_t max_substeps =
        std::max<std::size_t>(base_substeps, inputs.maximum_substeps);
    const double dt_safe =
        std::max(1.0e-12, std::isfinite(inputs.dt_s) ? inputs.dt_s : 0.0);
    const double capacitance_per_area =
        std::max(1.0e-12,
                 std::isfinite(inputs.capacitance_per_area_f_per_m2)
                     ? inputs.capacitance_per_area_f_per_m2
                     : 0.0);
    const double current_derivative =
        std::isfinite(inputs.current_derivative_a_per_m2_per_v)
            ? inputs.current_derivative_a_per_m2_per_v
            : 0.0;

    const double stiffness =
        std::abs(current_derivative) * dt_safe / capacitance_per_area;
    double sensitivity_scale = 1.0;
    if (inputs.reference_potential_v < -25.0)
    {
        const double negative_magnitude = std::abs(inputs.reference_potential_v);
        sensitivity_scale +=
            0.25 * std::clamp(negative_magnitude / 250.0, 0.0, 8.0);
    }

    const double target_scale = std::max(1.0, stiffness * sensitivity_scale);
    const double target_substeps =
        static_cast<double>(base_substeps) * target_scale;
    if (!std::isfinite(target_substeps))
    {
        return max_substeps;
    }

    return std::clamp<std::size_t>(
        static_cast<std::size_t>(std::ceil(target_substeps)), base_substeps,
        max_substeps);
}

double computeSurfaceSheathEquivalentCapacitancePerArea(
    const SurfaceSheathCapacitanceConsistencyInputs& inputs)
{
    const double minimum_capacitance =
        std::max(1.0e-12,
                 std::isfinite(inputs.minimum_capacitance_per_area_f_per_m2)
                     ? inputs.minimum_capacitance_per_area_f_per_m2
                     : 0.0);
    const double minimum_sheath_length_m =
        std::max(1.0e-6,
                 std::isfinite(inputs.minimum_sheath_length_m)
                     ? inputs.minimum_sheath_length_m
                     : 0.0);
    const double maximum_sheath_length_m =
        std::max(minimum_sheath_length_m,
                 std::isfinite(inputs.maximum_sheath_length_m)
                     ? inputs.maximum_sheath_length_m
                     : minimum_sheath_length_m);
    const double effective_sheath_length_m =
        std::isfinite(inputs.effective_sheath_length_m)
            ? inputs.effective_sheath_length_m
            : minimum_sheath_length_m;
    const double sheath_length_m = std::clamp(effective_sheath_length_m,
                                              minimum_sheath_length_m,
                                              maximum_sheath_length_m);

    const double relative_permittivity =
        std::max(1.0,
                 std::isfinite(inputs.relative_permittivity)
                     ? inputs.relative_permittivity
                     : 1.0);
    const double sheath_capacitance_per_area_f_per_m2 =
        kEpsilon0 * relative_permittivity / std::max(1.0e-6, sheath_length_m);
    const double dielectric_capacitance_per_area_f_per_m2 =
        std::max(minimum_capacitance,
                 std::isfinite(inputs.dielectric_capacitance_per_area_f_per_m2)
                     ? inputs.dielectric_capacitance_per_area_f_per_m2
                     : minimum_capacitance);

    return 1.0 /
           (1.0 / dielectric_capacitance_per_area_f_per_m2 +
            1.0 /
                std::max(minimum_capacitance,
                         sheath_capacitance_per_area_f_per_m2));
}

double enforceSurfaceSheathCapacitanceConsistency(
    const SurfaceSheathCapacitanceConsistencyInputs& inputs)
{
    const double minimum_capacitance =
        std::max(1.0e-12,
                 std::isfinite(inputs.minimum_capacitance_per_area_f_per_m2)
                     ? inputs.minimum_capacitance_per_area_f_per_m2
                     : 0.0);
    const double candidate_capacitance =
        std::max(minimum_capacitance,
                 std::isfinite(inputs.dielectric_capacitance_per_area_f_per_m2)
                     ? inputs.dielectric_capacitance_per_area_f_per_m2
                     : minimum_capacitance);

    SurfaceSheathCapacitanceConsistencyInputs equivalent_inputs = inputs;
    equivalent_inputs.dielectric_capacitance_per_area_f_per_m2 =
        candidate_capacitance;
    equivalent_inputs.minimum_capacitance_per_area_f_per_m2 = minimum_capacitance;
    const double consistent_series_capacitance =
        computeSurfaceSheathEquivalentCapacitancePerArea(equivalent_inputs);
    const double mismatch =
        std::abs(candidate_capacitance - consistent_series_capacitance) /
        std::max(minimum_capacitance,
                 std::max(candidate_capacitance,
                          consistent_series_capacitance));

    const double base_weight =
        std::clamp(std::isfinite(inputs.consistency_weight)
                       ? inputs.consistency_weight
                       : 0.60,
                   0.0, 0.95);
    const double coupling_gain =
        std::clamp(std::isfinite(inputs.volume_mesh_coupling_gain)
                       ? inputs.volume_mesh_coupling_gain
                       : 0.0,
                   0.0, 1.0);
    const double coupling_boost = 0.15 * coupling_gain;
    const double effective_weight =
        std::clamp(base_weight * (mismatch + coupling_boost), 0.0, 0.95);

    double blended_capacitance =
        (1.0 - effective_weight) * candidate_capacitance +
        effective_weight * consistent_series_capacitance;

    const double ratio_guard =
        std::clamp(std::isfinite(inputs.ratio_guard) ? inputs.ratio_guard : 1.5,
                   0.1, 8.0);
    const double minimum_allowed = candidate_capacitance / (1.0 + ratio_guard);
    const double maximum_allowed = candidate_capacitance * (1.0 + ratio_guard);
    blended_capacitance =
        std::clamp(blended_capacitance, minimum_allowed, maximum_allowed);

    return std::max(minimum_capacitance, blended_capacitance);
}

double computeIndirectPathBranchWeight(
    const std::vector<std::vector<double>>& branch_graph_adjacency_weight_f,
    std::size_t from_node, std::size_t to_node)
{
    if (from_node >= branch_graph_adjacency_weight_f.size() ||
        to_node >= branch_graph_adjacency_weight_f.size() || from_node == to_node)
    {
        return 0.0;
    }

    double path_branch_weight_f = 0.0;
    const std::size_t node_count = branch_graph_adjacency_weight_f.size();
    for (std::size_t intermediate_node = 0; intermediate_node < node_count;
         ++intermediate_node)
    {
        if (intermediate_node == from_node || intermediate_node == to_node)
        {
            continue;
        }
        if (intermediate_node >= branch_graph_adjacency_weight_f[from_node].size() ||
            intermediate_node >= branch_graph_adjacency_weight_f[to_node].size())
        {
            continue;
        }

        const double forward_weight_f =
            branch_graph_adjacency_weight_f[from_node][intermediate_node];
        const double backward_weight_f =
            branch_graph_adjacency_weight_f[to_node][intermediate_node];
        if (forward_weight_f <= 0.0 || backward_weight_f <= 0.0)
        {
            continue;
        }

        path_branch_weight_f +=
            2.0 * forward_weight_f * backward_weight_f /
            std::max(1.0e-30, forward_weight_f + backward_weight_f);
    }

    return path_branch_weight_f;
}

double computeHarmonicPairCapacitanceF(double area_i_m2, double area_j_m2,
                                       double effective_sheath_length_m,
                                       double capacitance_scale,
                                       double minimum_capacitance_f)
{
    const double safe_area_i_m2 = std::max(1.0e-16, area_i_m2);
    const double safe_area_j_m2 = std::max(1.0e-16, area_j_m2);
    const double harmonic_pair_area_m2 =
        2.0 * safe_area_i_m2 * safe_area_j_m2 /
        std::max(1.0e-16, safe_area_i_m2 + safe_area_j_m2);
    const double safe_sheath_length_m = std::max(1.0e-6, effective_sheath_length_m);
    const double scaled_capacitance_f =
        std::max(0.0, capacitance_scale) * kEpsilon0 * harmonic_pair_area_m2 /
        safe_sheath_length_m;
    return std::max(minimum_capacitance_f, scaled_capacitance_f);
}

double estimateCalibrationFactor(double measured_current_a, double native_current_a,
                                 double minimum_factor, double maximum_factor,
                                 double default_factor)
{
    if (!std::isfinite(measured_current_a) || !std::isfinite(native_current_a) ||
        std::abs(native_current_a) <= 1.0e-18)
    {
        return default_factor;
    }
    return std::clamp(measured_current_a / native_current_a, minimum_factor, maximum_factor);
}

} // namespace Coupling
} // namespace SCDAT
