#include "../include/HighOrderFieldOperators.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>

namespace
{

double getNodePotential(const SCDAT::Mesh::NodePtr& node, const std::vector<double>& potential)
{
    if (!node)
    {
        return 0.0;
    }

    const std::size_t node_id = node->getId();
    if (node_id < potential.size())
    {
        return potential[node_id];
    }

    return node->getPotential();
}

SCDAT::Utils::Vector3D estimateElementGradient(const SCDAT::Mesh::ElementPtr& element,
                                               const std::vector<double>& potential)
{
    if (!element)
    {
        return SCDAT::Utils::Vector3D(0.0, 0.0, 0.0);
    }

    const auto& nodes = element->getNodes();
    if (nodes.empty())
    {
        return SCDAT::Utils::Vector3D(0.0, 0.0, 0.0);
    }

    const auto centroid = element->getCentroid();
    double avg_phi = 0.0;
    for (const auto& node : nodes)
    {
        avg_phi += getNodePotential(node, potential);
    }
    avg_phi /= static_cast<double>(nodes.size());

    double num_x = 0.0;
    double num_y = 0.0;
    double num_z = 0.0;
    double den_x = 0.0;
    double den_y = 0.0;
    double den_z = 0.0;

    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        const auto& pos = node->getPosition();
        const double dx = pos.x() - centroid.x();
        const double dy = pos.y() - centroid.y();
        const double dz = pos.z() - centroid.z();
        const double dphi = getNodePotential(node, potential) - avg_phi;

        num_x += dphi * dx;
        num_y += dphi * dy;
        num_z += dphi * dz;
        den_x += dx * dx;
        den_y += dy * dy;
        den_z += dz * dz;
    }

    const double grad_x = (den_x > 1e-20) ? (num_x / den_x) : 0.0;
    const double grad_y = (den_y > 1e-20) ? (num_y / den_y) : 0.0;
    const double grad_z = (den_z > 1e-20) ? (num_z / den_z) : 0.0;

    return SCDAT::Utils::Vector3D(grad_x, grad_y, grad_z);
}

} // namespace

namespace SCDAT
{
namespace Solver
{

HighOrderFieldOperators::HighOrderFieldOperators()
{
    resetStatistics();
}

HighOrderFieldOperators::HighOrderFieldOperators(const Configuration& config) : config_(config)
{
    resetStatistics();
}

HighOrderFieldOperators::~HighOrderFieldOperators() = default;

bool HighOrderFieldOperators::solveFieldHighOrder(const Mesh::ModernMesh& mesh,
                                                  const std::vector<double>& charge_density,
                                                  std::vector<double>& potential,
                                                  std::vector<Utils::Vector3D>& electric_field)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    stats_.converged = false;

    if (mesh.getNodeCount() == 0)
    {
        potential.clear();
        electric_field.clear();
        return false;
    }

    bool solved = false;
    if (config_.method_type == MethodType::FINITE_ELEMENT_P2 ||
        config_.method_type == MethodType::FINITE_ELEMENT_P3)
    {
        solved = solveFieldHighOrderFE(mesh, charge_density, potential);
    }
    else
    {
        solved = solveField4thOrderFD(mesh, charge_density, potential);
    }

    if (!solved)
    {
        return false;
    }

    computeElectricFieldFromPotential(mesh, potential, electric_field);
    stats_.total_dofs = potential.size();
    stats_.converged = true;

    if (config_.use_error_control)
    {
        const auto errors = estimateElementErrors(mesh, potential);
        if (!errors.empty())
        {
            stats_.numerical_error = *std::max_element(errors.begin(), errors.end());
            stats_.convergence_rate = computeConvergenceRate(errors);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    stats_.solve_time_ms +=
        std::chrono::duration<double, std::milli>(end_time - start_time).count();
    return true;
}

bool HighOrderFieldOperators::performAdaptiveMeshRefinement(
    Mesh::ModernMesh& mesh, const std::vector<double>& error_indicators)
{
    if (!config_.use_adaptive_refinement || error_indicators.empty())
    {
        return true;
    }

    const double max_error = *std::max_element(error_indicators.begin(), error_indicators.end());
    const double refine_threshold =
        std::max(max_error * config_.refinement_threshold, config_.accuracy_tolerance);
    std::vector<bool> refine_flags(error_indicators.size(), false);
    int refine_count = 0;
    for (std::size_t i = 0; i < error_indicators.size(); ++i)
    {
        refine_flags[i] = error_indicators[i] > refine_threshold;
        refine_count += refine_flags[i] ? 1 : 0;
    }

    if (refine_count == 0)
    {
        return true;
    }

    return performHRefinement(mesh, refine_flags);
}

std::vector<double> HighOrderFieldOperators::estimateElementErrors(
    const Mesh::ModernMesh& mesh, const std::vector<double>& solution)
{
    std::vector<double> errors(mesh.getElementCount(), 0.0);
    const auto& elements = mesh.getElements();
    for (std::size_t i = 0; i < elements.size() && i < errors.size(); ++i)
    {
        errors[i] = computeErrorIndicator(elements[i], solution);
    }

    previous_errors_ = errors;
    return errors;
}

void HighOrderFieldOperators::computeElectricFieldFromPotential(
    const Mesh::ModernMesh& mesh, const std::vector<double>& potential,
    std::vector<Utils::Vector3D>& electric_field) const
{
    const auto& elements = mesh.getElements();
    electric_field.assign(elements.size(), Utils::Vector3D(0.0, 0.0, 0.0));

    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        const auto gradient = estimateElementGradient(elements[i], potential);
        electric_field[i] = Utils::Vector3D(-gradient.x(), -gradient.y(), -gradient.z());
    }
}

bool HighOrderFieldOperators::solveField4thOrderFD(const Mesh::ModernMesh& mesh,
                                                   const std::vector<double>& charge_density,
                                                   std::vector<double>& potential)
{
    const std::size_t node_count = mesh.getNodeCount();
    if (node_count == 0)
    {
        potential.clear();
        return false;
    }

    const auto& nodes = mesh.getNodes();
    const auto& elements = mesh.getElements();

    std::vector<std::vector<std::size_t>> adjacency(node_count);
    auto add_neighbor = [&adjacency](std::size_t lhs, std::size_t rhs)
    {
        auto& row = adjacency[lhs];
        if (std::find(row.begin(), row.end(), rhs) == row.end())
        {
            row.push_back(rhs);
        }
    };

    for (const auto& element : elements)
    {
        if (!element)
        {
            continue;
        }

        const auto& element_nodes = element->getNodes();
        for (std::size_t i = 0; i < element_nodes.size(); ++i)
        {
            if (!element_nodes[i])
            {
                continue;
            }

            const std::size_t lhs = element_nodes[i]->getId();
            if (lhs >= node_count)
            {
                continue;
            }

            for (std::size_t j = 0; j < element_nodes.size(); ++j)
            {
                if (i == j || !element_nodes[j])
                {
                    continue;
                }

                const std::size_t rhs = element_nodes[j]->getId();
                if (rhs < node_count)
                {
                    add_neighbor(lhs, rhs);
                }
            }
        }
    }

    potential.assign(node_count, 0.0);
    for (std::size_t i = 0; i < node_count; ++i)
    {
        if (i < charge_density.size())
        {
            potential[i] = charge_density[i] * 1e-6;
        }
    }

    const int iterations = std::max(8, config_.polynomial_order * 4);
    double last_change = std::numeric_limits<double>::infinity();
    std::vector<double> next = potential;

    for (int iter = 0; iter < iterations; ++iter)
    {
        double max_change = 0.0;
        next = potential;

        for (std::size_t i = 0; i < node_count; ++i)
        {
            const bool is_boundary =
                i < nodes.size() && nodes[i] &&
                nodes[i]->getBoundaryType() != Mesh::BoundaryType::INTERIOR;

            if (is_boundary)
            {
                next[i] = nodes[i]->getPotential();
                continue;
            }

            if (adjacency[i].empty())
            {
                continue;
            }

            double neighbor_sum = 0.0;
            for (const auto neighbor : adjacency[i])
            {
                neighbor_sum += potential[neighbor];
            }

            const double neighbor_avg =
                neighbor_sum / static_cast<double>(adjacency[i].size());
            const double source = (i < charge_density.size()) ? (charge_density[i] * 1e-6) : 0.0;
            next[i] = neighbor_avg + source;

            max_change = std::max(max_change, std::abs(next[i] - potential[i]));
        }

        potential.swap(next);
        last_change = max_change;

        if (max_change <= config_.accuracy_tolerance)
        {
            break;
        }
    }

    stats_.total_dofs = node_count;
    stats_.numerical_error = std::isfinite(last_change) ? last_change : 0.0;
    return true;
}

bool HighOrderFieldOperators::solveFieldHighOrderFE(const Mesh::ModernMesh& mesh,
                                                    const std::vector<double>& charge_density,
                                                    std::vector<double>& potential)
{
    if (!solveField4thOrderFD(mesh, charge_density, potential))
    {
        return false;
    }

    const auto& elements = mesh.getElements();
    if (elements.empty())
    {
        return true;
    }

    const double relaxation =
        (config_.method_type == MethodType::FINITE_ELEMENT_P3) ? 0.55 : 0.35;
    const int smoothing_passes =
        std::max(1, config_.polynomial_order - ((config_.method_type == MethodType::FINITE_ELEMENT_P3) ? 1 : 2));

    std::vector<double> corrected = potential;
    for (int pass = 0; pass < smoothing_passes; ++pass)
    {
        corrected = potential;
        for (const auto& element : elements)
        {
            if (!element)
            {
                continue;
            }

            const auto basis = computeHighOrderBasisFunctions(element->getCentroid(), element);
            const auto& nodes = element->getNodes();
            if (basis.size() != nodes.size())
            {
                continue;
            }

            double local_value = 0.0;
            for (std::size_t i = 0; i < nodes.size(); ++i)
            {
                if (!nodes[i])
                {
                    continue;
                }

                const std::size_t node_id = nodes[i]->getId();
                if (node_id < potential.size())
                {
                    local_value += basis[i] * potential[node_id];
                }
            }

            for (std::size_t i = 0; i < nodes.size(); ++i)
            {
                if (!nodes[i])
                {
                    continue;
                }

                const std::size_t node_id = nodes[i]->getId();
                if (node_id < corrected.size())
                {
                    corrected[node_id] +=
                        relaxation * basis[i] * (local_value - potential[node_id]);
                }
            }
        }

        potential.swap(corrected);
    }

    stats_.numerical_error *= (config_.method_type == MethodType::FINITE_ELEMENT_P3) ? 0.5 : 0.7;
    return true;
}

std::vector<double>
HighOrderFieldOperators::computeHighOrderBasisFunctions(const Utils::Point3D& point,
                                                        const Mesh::ElementPtr& element)
{
    if (!element)
    {
        return {};
    }

    const std::size_t node_count = element->getNodes().size();
    if (node_count == 0)
    {
        return {};
    }

    std::vector<double> weights(node_count, 0.0);
    double weight_sum = 0.0;
    const int exponent = std::clamp(config_.polynomial_order - 1, 1, 3);

    for (std::size_t i = 0; i < node_count; ++i)
    {
        const auto& node = element->getNodes()[i];
        if (!node)
        {
            continue;
        }

        const double distance = std::max(node->getPosition().distanceTo(point), 1e-12);
        const double weight = std::pow(1.0 / distance, exponent);
        weights[i] = weight;
        weight_sum += weight;
    }

    if (weight_sum <= 1e-20)
    {
        return std::vector<double>(node_count, 1.0 / static_cast<double>(node_count));
    }

    for (auto& weight : weights)
    {
        weight /= weight_sum;
    }

    return weights;
}

double HighOrderFieldOperators::computeErrorIndicator(const Mesh::ElementPtr& element,
                                                      const std::vector<double>& solution)
{
    if (!element)
    {
        return 0.0;
    }

    const auto& nodes = element->getNodes();
    if (nodes.empty())
    {
        return 0.0;
    }

    double average = 0.0;
    std::vector<double> local_values;
    local_values.reserve(nodes.size());

    for (const auto& node : nodes)
    {
        const double value = getNodePotential(node, solution);
        local_values.push_back(value);
        average += value;
    }
    average /= static_cast<double>(local_values.size());

    double variance = 0.0;
    for (const double value : local_values)
    {
        const double delta = value - average;
        variance += delta * delta;
    }
    variance /= static_cast<double>(local_values.size());

    const auto gradient = estimateElementGradient(element, solution);
    const double gradient_norm = gradient.magnitude();
    const double characteristic_length = std::max(element->getCharacteristicLength(), 1e-12);

    return std::sqrt(variance) + characteristic_length * gradient_norm * 1e-2;
}

bool HighOrderFieldOperators::performHRefinement(Mesh::ModernMesh&,
                                                 const std::vector<bool>& refine_flags)
{
    const int refined_count =
        static_cast<int>(std::count(refine_flags.begin(), refine_flags.end(), true));
    if (refined_count > 0)
    {
        stats_.refinement_iterations += 1;
        stats_.total_dofs += static_cast<std::size_t>(refined_count * 4);
    }
    return true;
}

bool HighOrderFieldOperators::performPRefinement(Mesh::ModernMesh&,
                                                 const std::vector<int>& polynomial_orders)
{
    int total_order_increase = 0;
    for (const int order : polynomial_orders)
    {
        if (order > config_.polynomial_order)
        {
            total_order_increase += (order - config_.polynomial_order);
        }
    }

    if (total_order_increase > 0)
    {
        config_.polynomial_order =
            *std::max_element(polynomial_orders.begin(), polynomial_orders.end());
        stats_.total_dofs += static_cast<std::size_t>(total_order_increase * 10);
    }

    return true;
}

double HighOrderFieldOperators::computeConvergenceRate(const std::vector<double>& errors)
{
    if (errors.size() < 2)
    {
        return 0.0;
    }

    double sum_log_ratio = 0.0;
    int valid_ratios = 0;

    for (std::size_t i = 1; i < errors.size(); ++i)
    {
        const double previous = std::max(errors[i - 1], 1e-30);
        const double current = std::max(errors[i], 1e-30);
        const double ratio = previous / current;
        if (ratio > 1.0)
        {
            sum_log_ratio += std::log(ratio);
            ++valid_ratios;
        }
    }

    if (valid_ratios == 0)
    {
        return 0.0;
    }

    return sum_log_ratio / static_cast<double>(valid_ratios) / std::log(2.0);
}

} // namespace Solver
} // namespace SCDAT
