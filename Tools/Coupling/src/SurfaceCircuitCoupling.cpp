#include "../include/SurfaceCircuitCoupling.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SCDAT
{
namespace Coupling
{
namespace
{
constexpr double kTwoPi = 6.28318530717958647692;

void accumulateDidvTotals(double weight, double current_a, double didv_a_per_v,
                          double& total_current_a, double& total_didv_a_per_v)
{
    if (!std::isfinite(weight))
    {
        return;
    }
    if (std::isfinite(current_a))
    {
        total_current_a += weight * current_a;
    }
    if (std::isfinite(didv_a_per_v))
    {
        total_didv_a_per_v += weight * didv_a_per_v;
    }
}

void accumulateDirectTerms(const std::vector<CircuitDidvTerm>& terms,
                           double& total_current_a, double& total_didv_a_per_v)
{
    for (const auto& term : terms)
    {
        accumulateDidvTotals(1.0, term.current_a, term.didv_a_per_v, total_current_a,
                             total_didv_a_per_v);
    }
}

void accumulateWeightedTerms(const std::vector<CircuitWeightedDidvTerm>& weighted_terms,
                             double& total_current_a, double& total_didv_a_per_v,
                             bool& consumed)
{
    for (const auto& term : weighted_terms)
    {
        consumed = true;
        accumulateDidvTotals(term.weight, term.current_a, term.didv_a_per_v, total_current_a,
                             total_didv_a_per_v);
    }
}

void accumulateSurfaceCompositions(
    const std::vector<CircuitDidvSurfaceComposition>& surface_compositions,
    double& total_current_a, double& total_didv_a_per_v, bool& consumed)
{
    for (const auto& surface : surface_compositions)
    {
        double surface_current_a = 0.0;
        double surface_didv_a_per_v = 0.0;
        for (const auto& term : surface.terms)
        {
            consumed = true;
            accumulateDidvTotals(1.0, term.current_a, term.didv_a_per_v, surface_current_a,
                                 surface_didv_a_per_v);
        }
        accumulateDidvTotals(surface.surface_weight, surface_current_a, surface_didv_a_per_v,
                             total_current_a, total_didv_a_per_v);
    }
}

bool applyCircuitElement(SurfaceCircuitAssemblyResult& result,
                         const CircuitElementDescriptor& element)
{
    switch (element.kind)
    {
    case CircuitElementKind::NodeCapacitance:
        if (element.target_index < result.nodes.size() && std::isfinite(element.value))
        {
            const double capacitance = std::max(0.0, element.value);
            result.nodes[element.target_index].capacitance_f = capacitance;
            auto& node_state = result.kernel_input.nodes[element.target_index];
            node_state.capacitance_f = capacitance;
            node_state.valid = true;
            return true;
        }
        return false;
    case CircuitElementKind::NodeShuntConductance:
        if (element.target_index < result.nodes.size() && std::isfinite(element.value))
        {
            const double conductance = std::max(0.0, element.value);
            result.nodes[element.target_index].shunt_conductance_s = conductance;
            auto& node_state = result.kernel_input.nodes[element.target_index];
            node_state.shunt_conductance_s = conductance;
            node_state.valid = true;
            return true;
        }
        return false;
    case CircuitElementKind::NodeFixedPotential:
        if (element.target_index < result.nodes.size() && std::isfinite(element.value))
        {
            auto& node = result.nodes[element.target_index];
            node.fixed_potential = true;
            node.fixed_value_v = element.value;
            node.potential_v = element.value;
            return true;
        }
        return false;
    case CircuitElementKind::BranchConductance:
        if (element.target_index < result.branches.size() && std::isfinite(element.value))
        {
            result.branches[element.target_index].conductance_s = std::max(0.0, element.value);
            return true;
        }
        return false;
    case CircuitElementKind::BranchBias:
        if (element.target_index < result.branches.size() && std::isfinite(element.value))
        {
            result.branches[element.target_index].bias_v = element.value;
            return true;
        }
        return false;
    case CircuitElementKind::BranchCurrentSource:
        if (element.target_index < result.kernel_input.branches.size() &&
            std::isfinite(element.value))
        {
            auto& branch_state = result.kernel_input.branches[element.target_index];
            branch_state.current_source_a = element.value;
            branch_state.valid = true;
            return true;
        }
        return false;
    }

    return false;
}

bool isExcitationActive(const CircuitExcitationDescriptor& excitation, double time_s)
{
    if (!excitation.enabled || !std::isfinite(time_s))
    {
        return false;
    }
    if (std::isfinite(excitation.start_time_s) && time_s < excitation.start_time_s)
    {
        return false;
    }
    if (std::isfinite(excitation.end_time_s) && time_s > excitation.end_time_s)
    {
        return false;
    }
    return true;
}

double evaluatePiecewiseLinearExcitation(const CircuitExcitationDescriptor& excitation,
                                         double local_time_s)
{
    std::vector<CircuitExcitationPwlPoint> points;
    points.reserve(excitation.pwl_points.size());
    for (const auto& point : excitation.pwl_points)
    {
        if (std::isfinite(point.time_s) && std::isfinite(point.value))
        {
            points.push_back(point);
        }
    }

    if (points.empty())
    {
        return excitation.base_value + excitation.amplitude;
    }

    std::sort(points.begin(), points.end(),
              [](const CircuitExcitationPwlPoint& left,
                 const CircuitExcitationPwlPoint& right) {
                  return left.time_s < right.time_s;
              });

    if (local_time_s <= points.front().time_s)
    {
        return excitation.base_value + points.front().value;
    }
    if (local_time_s >= points.back().time_s)
    {
        return excitation.base_value + points.back().value;
    }

    for (std::size_t index = 0; index + 1 < points.size(); ++index)
    {
        const auto& left = points[index];
        const auto& right = points[index + 1];
        if (local_time_s > right.time_s)
        {
            continue;
        }

        const double span_s = right.time_s - left.time_s;
        if (span_s <= 0.0)
        {
            continue;
        }

        const double ratio = std::clamp((local_time_s - left.time_s) / span_s, 0.0, 1.0);
        return excitation.base_value + left.value + (right.value - left.value) * ratio;
    }

    return excitation.base_value + points.back().value;
}

double evaluateExcitationValue(const CircuitExcitationDescriptor& excitation, double time_s)
{
    const double local_time_s =
        std::isfinite(excitation.start_time_s) ? std::max(0.0, time_s - excitation.start_time_s)
                                               : time_s;
    switch (excitation.waveform)
    {
    case CircuitExcitationWaveformKind::Constant:
        return excitation.base_value + excitation.amplitude;
    case CircuitExcitationWaveformKind::Sinusoidal:
    {
        const double omega = kTwoPi * std::max(0.0, excitation.frequency_hz);
        return excitation.base_value + excitation.amplitude *
                                           std::sin(omega * local_time_s + excitation.phase_rad);
    }
    case CircuitExcitationWaveformKind::Pulse:
    {
        const double frequency_hz = std::max(0.0, excitation.frequency_hz);
        if (frequency_hz <= 0.0)
        {
            return excitation.base_value + excitation.amplitude;
        }
        const double period_s = 1.0 / frequency_hz;
        const double normalized_phase =
            std::fmod(local_time_s / period_s + excitation.phase_rad / kTwoPi, 1.0);
        const double wrapped_phase = normalized_phase < 0.0 ? normalized_phase + 1.0
                                                             : normalized_phase;
        const double duty = std::clamp(excitation.duty_cycle, 0.0, 1.0);
        return excitation.base_value + (wrapped_phase <= duty ? excitation.amplitude : 0.0);
    }
    case CircuitExcitationWaveformKind::PiecewiseLinear:
        return evaluatePiecewiseLinearExcitation(excitation, local_time_s);
    case CircuitExcitationWaveformKind::Exponential:
    {
        if (!std::isfinite(excitation.exponential_tau_s) ||
            excitation.exponential_tau_s <= 0.0)
        {
            return excitation.base_value + excitation.amplitude;
        }
        return excitation.base_value +
               excitation.amplitude *
                   (1.0 - std::exp(-std::max(0.0, local_time_s) / excitation.exponential_tau_s));
    }
    }

    return excitation.base_value;
}
} // namespace

SurfaceCircuitAssemblyResult assembleSurfaceCircuit(const CircuitAssembly& assembly)
{
    return assembleSurfaceCircuit(assembly, {});
}

SurfaceCircuitAssemblyResult assembleSurfaceCircuit(
    const CircuitAssembly& assembly,
    const std::vector<CircuitElementDescriptor>& dynamic_elements)
{
    SurfaceCircuitAssemblyResult result;
    result.nodes.reserve(assembly.nodes.size());
    for (const auto& descriptor : assembly.nodes)
    {
        SurfaceCircuitNode node;
        node.name = descriptor.name;
        node.potential_v = descriptor.initial_potential_v;
        result.nodes.push_back(node);
    }

    result.branches.reserve(assembly.branches.size());
    for (const auto& descriptor : assembly.branches)
    {
        SurfaceCircuitBranch branch;
        branch.from_node = descriptor.from_node;
        branch.to_node = descriptor.to_node;
        result.branches.push_back(branch);
    }

    result.kernel_input.nodes.resize(result.nodes.size());
    result.kernel_input.branches.resize(result.branches.size());

    for (const auto& element : assembly.elements)
    {
        applyCircuitElement(result, element);
    }

    for (const auto& element : dynamic_elements)
    {
        applyCircuitElement(result, element);
    }

    return result;
}

SurfaceCircuitKernelInput composeCircuitKernelInput(
    const CircuitDidvCompositionInput& didv_input,
    const CircuitKernelTopologyState& topology_state)
{
    const auto didv_result = composeCircuitDidv(didv_input);
    SurfaceCircuitKernelInput kernel_input = didv_result.kernel_input;

    for (std::size_t node_index = 0; node_index < kernel_input.nodes.size(); ++node_index)
    {
        auto& node_state = kernel_input.nodes[node_index];
        if (node_index < topology_state.node_capacitance_f.size() &&
            std::isfinite(topology_state.node_capacitance_f[node_index]))
        {
            node_state.capacitance_f =
                std::max(0.0, topology_state.node_capacitance_f[node_index]);
        }
        if (node_index < topology_state.node_shunt_conductance_s.size() &&
            std::isfinite(topology_state.node_shunt_conductance_s[node_index]))
        {
            node_state.shunt_conductance_s =
                std::max(0.0, topology_state.node_shunt_conductance_s[node_index]);
        }
        node_state.valid = true;
    }

    const std::size_t branch_count =
        std::max({topology_state.branch_conductance_s.size(),
                  topology_state.branch_bias_v.size(),
                  topology_state.branch_current_source_a.size()});
    kernel_input.branches.resize(branch_count);
    for (std::size_t branch_index = 0; branch_index < branch_count; ++branch_index)
    {
        auto& branch_state = kernel_input.branches[branch_index];
        if (branch_index < topology_state.branch_conductance_s.size() &&
            std::isfinite(topology_state.branch_conductance_s[branch_index]))
        {
            branch_state.conductance_s =
                std::max(0.0, topology_state.branch_conductance_s[branch_index]);
        }
        if (branch_index < topology_state.branch_bias_v.size() &&
            std::isfinite(topology_state.branch_bias_v[branch_index]))
        {
            branch_state.bias_v = topology_state.branch_bias_v[branch_index];
        }
        if (branch_index < topology_state.branch_current_source_a.size() &&
            std::isfinite(topology_state.branch_current_source_a[branch_index]))
        {
            branch_state.current_source_a =
                topology_state.branch_current_source_a[branch_index];
        }
        branch_state.valid = true;
    }

    return kernel_input;
}

CircuitExcitationSampleResult sampleCircuitExcitations(
    const std::vector<CircuitExcitationDescriptor>& excitations, double time_s)
{
    CircuitExcitationSampleResult result;
    result.elements.reserve(excitations.size());
    for (const auto& excitation : excitations)
    {
        if (!isExcitationActive(excitation, time_s))
        {
            continue;
        }

        const double value = evaluateExcitationValue(excitation, time_s);
        if (!std::isfinite(value))
        {
            continue;
        }

        result.elements.push_back(
            CircuitElementDescriptor{excitation.target_kind, excitation.target_index, value});
    }
    return result;
}

CircuitDidvCompositionResult composeCircuitDidv(const CircuitDidvCompositionInput& input)
{
    CircuitDidvCompositionResult result;
    const std::size_t node_count = input.nodes.size();
    result.node_summaries.resize(node_count);
    result.linearization.plasma_current_a.resize(node_count, 0.0);
    result.linearization.plasma_didv_a_per_v.resize(node_count, 0.0);
    result.linearization.additional_rhs_a.resize(node_count, 0.0);
    result.linearization.additional_diagonal_a_per_v.resize(node_count, 0.0);
    result.linearization.additional_off_diagonal_entries = input.off_diagonal_entries;

    result.kernel_input.nodes.resize(node_count);
    result.kernel_input.off_diagonal_entries = input.off_diagonal_entries;

    for (std::size_t node_index = 0; node_index < node_count; ++node_index)
    {
        const auto& node = input.nodes[node_index];
        double total_current_a = 0.0;
        double total_didv_a_per_v = 0.0;

        bool consumed_semantic_terms = false;
        switch (node.aggregation_kind)
        {
        case CircuitDidvAggregationKind::DirectSum:
            accumulateDirectTerms(node.terms, total_current_a, total_didv_a_per_v);
            consumed_semantic_terms = true;
            break;
        case CircuitDidvAggregationKind::WeightedTerms:
            accumulateWeightedTerms(node.weighted_terms, total_current_a, total_didv_a_per_v,
                                    consumed_semantic_terms);
            break;
        case CircuitDidvAggregationKind::MultiSurface:
            accumulateSurfaceCompositions(node.surface_compositions, total_current_a,
                                          total_didv_a_per_v, consumed_semantic_terms);
            break;
        }

        if (!consumed_semantic_terms)
        {
            accumulateDirectTerms(node.terms, total_current_a, total_didv_a_per_v);
        }

        result.node_summaries[node_index] =
            CircuitDidvNodeSummary{total_current_a, total_didv_a_per_v};
        result.linearization.plasma_current_a[node_index] = total_current_a;
        result.linearization.plasma_didv_a_per_v[node_index] = total_didv_a_per_v;
        result.linearization.additional_rhs_a[node_index] = node.additional_rhs_a;
        result.linearization.additional_diagonal_a_per_v[node_index] =
            node.additional_diagonal_a_per_v;

        auto& kernel_node = result.kernel_input.nodes[node_index];
        kernel_node.plasma_current_a = total_current_a;
        kernel_node.plasma_didv_a_per_v = total_didv_a_per_v;
        kernel_node.additional_rhs_a = node.additional_rhs_a;
        kernel_node.additional_diagonal_a_per_v = node.additional_diagonal_a_per_v;
        kernel_node.valid = true;
    }

    return result;
}

void SurfaceCircuitCoupling::clear()
{
    nodes_.clear();
    branches_.clear();
}

void SurfaceCircuitCoupling::configure(const CircuitAssembly& assembly)
{
    const auto assembled = assembleSurfaceCircuit(assembly);
    nodes_ = assembled.nodes;
    branches_ = assembled.branches;
}

std::size_t SurfaceCircuitCoupling::addNode(const SurfaceCircuitNode& node)
{
    nodes_.push_back(node);
    return nodes_.size() - 1;
}

void SurfaceCircuitCoupling::addBranch(const SurfaceCircuitBranch& branch)
{
    branches_.push_back(branch);
}

void SurfaceCircuitCoupling::setNodePotential(std::size_t node_index, double potential_v)
{
    if (node_index < nodes_.size())
    {
        nodes_[node_index].potential_v = potential_v;
    }
}

void SurfaceCircuitCoupling::setNodeFixedPotential(std::size_t node_index, bool fixed,
                                                   double fixed_value_v)
{
    if (node_index < nodes_.size())
    {
        nodes_[node_index].fixed_potential = fixed;
        nodes_[node_index].fixed_value_v = fixed_value_v;
        if (fixed)
        {
            nodes_[node_index].potential_v = fixed_value_v;
        }
    }
}

void SurfaceCircuitCoupling::setNodeShuntConductance(std::size_t node_index, double conductance_s)
{
    if (node_index < nodes_.size())
    {
        nodes_[node_index].shunt_conductance_s = std::max(0.0, conductance_s);
    }
}

void SurfaceCircuitCoupling::setNodeCapacitance(std::size_t node_index, double capacitance_f)
{
    if (node_index < nodes_.size())
    {
        nodes_[node_index].capacitance_f = std::max(0.0, capacitance_f);
    }
}

void SurfaceCircuitCoupling::setBranchConductance(std::size_t branch_index, double conductance_s)
{
    if (branch_index < branches_.size())
    {
        branches_[branch_index].conductance_s = std::max(0.0, conductance_s);
    }
}

std::vector<double>
SurfaceCircuitCoupling::solveDenseLinearSystem(std::vector<std::vector<double>> matrix,
                                               std::vector<double> rhs)
{
    const std::size_t dimension = rhs.size();
    if (dimension == 0)
    {
        return {};
    }

    for (std::size_t pivot = 0; pivot < dimension; ++pivot)
    {
        std::size_t best_row = pivot;
        double best_value = std::abs(matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < dimension; ++row)
        {
            const double candidate = std::abs(matrix[row][pivot]);
            if (candidate > best_value)
            {
                best_value = candidate;
                best_row = row;
            }
        }

        if (best_value < 1.0e-24)
        {
            throw std::runtime_error("SurfaceCircuitCoupling singular matrix");
        }

        if (best_row != pivot)
        {
            std::swap(matrix[pivot], matrix[best_row]);
            std::swap(rhs[pivot], rhs[best_row]);
        }

        const double diagonal = matrix[pivot][pivot];
        for (std::size_t column = pivot; column < dimension; ++column)
        {
            matrix[pivot][column] /= diagonal;
        }
        rhs[pivot] /= diagonal;

        for (std::size_t row = 0; row < dimension; ++row)
        {
            if (row == pivot)
            {
                continue;
            }

            const double factor = matrix[row][pivot];
            if (std::abs(factor) < 1.0e-24)
            {
                continue;
            }

            for (std::size_t column = pivot; column < dimension; ++column)
            {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    return rhs;
}

SurfaceCircuitAdvanceResult SurfaceCircuitCoupling::advanceImplicit(
    double dt, const SurfaceCircuitLinearization& linearization, double max_delta_potential_v)
{
    SurfaceCircuitKernelInput kernel_input;
    kernel_input.nodes.resize(nodes_.size());
    for (std::size_t index = 0; index < nodes_.size(); ++index)
    {
        auto& node_state = kernel_input.nodes[index];
        if (index < linearization.plasma_current_a.size())
        {
            node_state.plasma_current_a = linearization.plasma_current_a[index];
        }
        if (index < linearization.plasma_didv_a_per_v.size())
        {
            node_state.plasma_didv_a_per_v = linearization.plasma_didv_a_per_v[index];
        }
        if (index < linearization.additional_rhs_a.size())
        {
            node_state.additional_rhs_a = linearization.additional_rhs_a[index];
        }
        if (index < linearization.additional_diagonal_a_per_v.size())
        {
            node_state.additional_diagonal_a_per_v =
                linearization.additional_diagonal_a_per_v[index];
        }
        node_state.valid = true;
    }
    kernel_input.off_diagonal_entries = linearization.additional_off_diagonal_entries;
    return advanceImplicit(dt, kernel_input, max_delta_potential_v);
}

SurfaceCircuitAdvanceResult SurfaceCircuitCoupling::advanceImplicit(
    double dt, const SurfaceCircuitKernelInput& kernel_input, double max_delta_potential_v)
{
    SurfaceCircuitAdvanceResult result;
    result.node_potentials_v.resize(nodes_.size(), 0.0);
    result.plasma_currents_a.resize(nodes_.size(), 0.0);
    result.branch_currents_a.resize(branches_.size(), 0.0);

    if (nodes_.empty())
    {
        result.converged = true;
        return result;
    }

    const std::size_t node_count = nodes_.size();
    std::vector<double> old_potentials(node_count, 0.0);
    std::vector<double> plasma_current(node_count, 0.0);
    std::vector<double> plasma_didv(node_count, 0.0);
    std::vector<double> additional_rhs(node_count, 0.0);
    std::vector<double> additional_diagonal(node_count, 0.0);
    for (std::size_t index = 0; index < node_count; ++index)
    {
        old_potentials[index] = nodes_[index].potential_v;
    }
    for (std::size_t index = 0; index < std::min(node_count, kernel_input.nodes.size()); ++index)
    {
        const auto& node_state = kernel_input.nodes[index];
        plasma_current[index] = node_state.plasma_current_a;
        plasma_didv[index] = node_state.plasma_didv_a_per_v;
        additional_rhs[index] = node_state.additional_rhs_a;
        additional_diagonal[index] = node_state.additional_diagonal_a_per_v;
        if (node_state.valid)
        {
            if (std::isfinite(node_state.capacitance_f))
            {
                nodes_[index].capacitance_f = std::max(0.0, node_state.capacitance_f);
            }
            if (std::isfinite(node_state.shunt_conductance_s))
            {
                nodes_[index].shunt_conductance_s =
                    std::max(0.0, node_state.shunt_conductance_s);
            }
        }
    }

    std::vector<std::size_t> free_nodes;
    free_nodes.reserve(node_count);
    for (std::size_t index = 0; index < node_count; ++index)
    {
        if (!nodes_[index].fixed_potential)
        {
            free_nodes.push_back(index);
        }
        else
        {
            nodes_[index].potential_v = nodes_[index].fixed_value_v;
        }
    }

    if (!free_nodes.empty())
    {
        const std::size_t free_count = free_nodes.size();
        std::vector<std::vector<double>> matrix(free_count, std::vector<double>(free_count, 0.0));
        std::vector<double> rhs(free_count, 0.0);
        std::vector<std::size_t> free_index(node_count, static_cast<std::size_t>(-1));
        for (std::size_t local = 0; local < free_count; ++local)
        {
            free_index[free_nodes[local]] = local;
        }

        for (std::size_t local = 0; local < free_count; ++local)
        {
            const std::size_t node_index = free_nodes[local];
            const auto& node = nodes_[node_index];
            const double capacitance = std::max(0.0, node.capacitance_f);
            const double dt_scale = dt > 0.0 ? capacitance / dt : 0.0;
            matrix[local][local] +=
                dt_scale - plasma_didv[node_index] + node.shunt_conductance_s +
                additional_diagonal[node_index];
            rhs[local] += dt_scale * node.potential_v + plasma_current[node_index] -
                          plasma_didv[node_index] * node.potential_v +
                          additional_rhs[node_index];
        }

        for (const auto& entry : kernel_input.off_diagonal_entries)
        {
            if (entry.row_node >= node_count || entry.column_node >= node_count ||
                !std::isfinite(entry.coefficient_a_per_v) ||
                std::abs(entry.coefficient_a_per_v) <= 0.0)
            {
                continue;
            }

            const auto row_local = free_index[entry.row_node];
            if (row_local == static_cast<std::size_t>(-1))
            {
                continue;
            }

            const auto column_local = free_index[entry.column_node];
            if (column_local != static_cast<std::size_t>(-1))
            {
                matrix[row_local][column_local] += entry.coefficient_a_per_v;
            }
            else
            {
                rhs[row_local] -= entry.coefficient_a_per_v *
                                  nodes_[entry.column_node].potential_v;
            }
        }

        for (std::size_t branch_index = 0; branch_index < branches_.size(); ++branch_index)
        {
            auto branch = branches_[branch_index];
            if (branch_index < kernel_input.branches.size() && kernel_input.branches[branch_index].valid)
            {
                const auto& branch_state = kernel_input.branches[branch_index];
                if (std::isfinite(branch_state.conductance_s))
                {
                    branch.conductance_s = std::max(0.0, branch_state.conductance_s);
                }
                if (std::isfinite(branch_state.bias_v))
                {
                    branch.bias_v = branch_state.bias_v;
                }
            }
            if (branch.from_node >= node_count || branch.to_node >= node_count ||
                branch.conductance_s <= 0.0)
            {
                if (branch_index >= kernel_input.branches.size() ||
                    !kernel_input.branches[branch_index].valid ||
                    std::abs(kernel_input.branches[branch_index].current_source_a) <= 0.0)
                {
                    continue;
                }
            }

            const auto from_local = free_index[branch.from_node];
            const auto to_local = free_index[branch.to_node];
            const bool from_free = from_local != static_cast<std::size_t>(-1);
            const bool to_free = to_local != static_cast<std::size_t>(-1);
            const double conductance = branch.conductance_s;
            const double branch_current_source_a =
                branch_index < kernel_input.branches.size() && kernel_input.branches[branch_index].valid
                    ? kernel_input.branches[branch_index].current_source_a
                    : 0.0;

            if (from_free)
            {
                matrix[from_local][from_local] += conductance;
                rhs[from_local] += conductance * branch.bias_v - branch_current_source_a;
                if (to_free)
                {
                    matrix[from_local][to_local] -= conductance;
                }
                else
                {
                    rhs[from_local] += conductance * nodes_[branch.to_node].potential_v;
                }
            }

            if (to_free)
            {
                matrix[to_local][to_local] += conductance;
                rhs[to_local] -= conductance * branch.bias_v + branch_current_source_a;
                if (from_free)
                {
                    matrix[to_local][from_local] -= conductance;
                }
                else
                {
                    rhs[to_local] += conductance * nodes_[branch.from_node].potential_v;
                }
            }
        }

        const auto solved = solveDenseLinearSystem(matrix, rhs);
        for (std::size_t local = 0; local < free_count; ++local)
        {
            const std::size_t node_index = free_nodes[local];
            double new_potential = solved[local];
            if (std::isfinite(max_delta_potential_v))
            {
                new_potential =
                    std::clamp(new_potential, old_potentials[node_index] - max_delta_potential_v,
                               old_potentials[node_index] + max_delta_potential_v);
            }
            nodes_[node_index].potential_v = new_potential;
            result.max_delta_potential_v =
                std::max(result.max_delta_potential_v,
                         std::abs(new_potential - old_potentials[node_index]));
        }
    }

    for (std::size_t index = 0; index < node_count; ++index)
    {
        result.node_potentials_v[index] = nodes_[index].potential_v;
        result.plasma_currents_a[index] =
            plasma_current[index] +
            plasma_didv[index] * (nodes_[index].potential_v - old_potentials[index]);
    }

    for (std::size_t branch_index = 0; branch_index < branches_.size(); ++branch_index)
    {
        auto branch = branches_[branch_index];
        if (branch_index < kernel_input.branches.size() && kernel_input.branches[branch_index].valid)
        {
            const auto& branch_state = kernel_input.branches[branch_index];
            if (std::isfinite(branch_state.conductance_s))
            {
                branch.conductance_s = std::max(0.0, branch_state.conductance_s);
            }
            if (std::isfinite(branch_state.bias_v))
            {
                branch.bias_v = branch_state.bias_v;
            }
        }
        if (branch.from_node >= node_count || branch.to_node >= node_count)
        {
            continue;
        }
        const double branch_current_source_a =
            branch_index < kernel_input.branches.size() && kernel_input.branches[branch_index].valid
                ? kernel_input.branches[branch_index].current_source_a
                : 0.0;
        result.branch_currents_a[branch_index] =
            branch.conductance_s * (nodes_[branch.from_node].potential_v -
                                    nodes_[branch.to_node].potential_v - branch.bias_v) +
            branch_current_source_a;
    }

    result.converged = true;
    return result;
}

} // namespace Coupling
} // namespace SCDAT
