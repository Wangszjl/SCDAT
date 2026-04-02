#include "../include/SurfaceCircuitCoupling.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SCDAT
{
namespace Coupling
{

void SurfaceCircuitCoupling::clear()
{
    nodes_.clear();
    branches_.clear();
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
    for (std::size_t index = 0; index < std::min(node_count, linearization.plasma_current_a.size());
         ++index)
    {
        plasma_current[index] = linearization.plasma_current_a[index];
    }
    for (std::size_t index = 0;
         index < std::min(node_count, linearization.plasma_didv_a_per_v.size()); ++index)
    {
        plasma_didv[index] = linearization.plasma_didv_a_per_v[index];
    }
    for (std::size_t index = 0;
         index < std::min(node_count, linearization.additional_rhs_a.size()); ++index)
    {
        additional_rhs[index] = linearization.additional_rhs_a[index];
    }
    for (std::size_t index = 0;
         index < std::min(node_count, linearization.additional_diagonal_a_per_v.size()); ++index)
    {
        additional_diagonal[index] = linearization.additional_diagonal_a_per_v[index];
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

        for (const auto& entry : linearization.additional_off_diagonal_entries)
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

        for (const auto& branch : branches_)
        {
            if (branch.from_node >= node_count || branch.to_node >= node_count ||
                branch.conductance_s <= 0.0)
            {
                continue;
            }

            const auto from_local = free_index[branch.from_node];
            const auto to_local = free_index[branch.to_node];
            const bool from_free = from_local != static_cast<std::size_t>(-1);
            const bool to_free = to_local != static_cast<std::size_t>(-1);
            const double conductance = branch.conductance_s;

            if (from_free)
            {
                matrix[from_local][from_local] += conductance;
                rhs[from_local] += conductance * branch.bias_v;
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
                rhs[to_local] -= conductance * branch.bias_v;
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
        const auto& branch = branches_[branch_index];
        if (branch.from_node >= node_count || branch.to_node >= node_count)
        {
            continue;
        }
        result.branch_currents_a[branch_index] =
            branch.conductance_s * (nodes_[branch.from_node].potential_v -
                                    nodes_[branch.to_node].potential_v - branch.bias_v);
    }

    result.converged = true;
    return result;
}

} // namespace Coupling
} // namespace SCDAT
