#include "MultiBoundaryElectricFieldSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace SCDAT
{
namespace FieldSolver
{

namespace
{

constexpr double kVacuumPermittivity = 8.8541878128e-12;
constexpr double kGeometryTolerance = 1.0e-12;
constexpr double kBreakdownTolerance = 1.0e-30;

double resolvePermittivity(double value)
{
    if (value <= 0.0)
    {
        return kVacuumPermittivity;
    }
    return (value < 1.0e-8) ? value : value * kVacuumPermittivity;
}

double dotProduct(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < lhs.size() && i < rhs.size(); ++i)
    {
        sum += lhs[i] * rhs[i];
    }
    return sum;
}

std::vector<double> multiplyMatrix(const std::vector<std::vector<double>>& matrix,
                                   const std::vector<double>& vector)
{
    std::vector<double> result(matrix.size(), 0.0);
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        for (std::size_t j = 0; j < matrix[i].size() && j < vector.size(); ++j)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

bool solve3x3(std::array<std::array<double, 3>, 3> matrix, std::array<double, 3> rhs,
              std::array<double, 3>& solution)
{
    solution = {0.0, 0.0, 0.0};

    for (std::size_t pivot = 0; pivot < 3; ++pivot)
    {
        std::size_t best = pivot;
        for (std::size_t row = pivot + 1; row < 3; ++row)
        {
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot]))
            {
                best = row;
            }
        }

        if (std::abs(matrix[best][pivot]) < 1.0e-18)
        {
            return false;
        }

        if (best != pivot)
        {
            std::swap(matrix[best], matrix[pivot]);
            std::swap(rhs[best], rhs[pivot]);
        }

        for (std::size_t row = pivot + 1; row < 3; ++row)
        {
            const double factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (std::size_t col = pivot; col < 3; ++col)
            {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    for (int row = 2; row >= 0; --row)
    {
        double sum = rhs[static_cast<std::size_t>(row)];
        for (std::size_t col = static_cast<std::size_t>(row) + 1; col < 3; ++col)
        {
            sum -= matrix[static_cast<std::size_t>(row)][col] * solution[col];
        }
        solution[static_cast<std::size_t>(row)] =
            sum / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }

    return true;
}

std::unordered_map<Mesh::NodeId, std::vector<Mesh::NodeId>>
buildNodeNeighbors(const std::vector<Mesh::ElementPtr>& elements)
{
    std::unordered_map<Mesh::NodeId, std::unordered_set<Mesh::NodeId>> unique_neighbors;

    for (const auto& element : elements)
    {
        if (!element)
        {
            continue;
        }

        const auto& nodes = element->getNodes();
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            if (!nodes[i])
            {
                continue;
            }

            for (std::size_t j = 0; j < nodes.size(); ++j)
            {
                if (i == j || !nodes[j])
                {
                    continue;
                }
                unique_neighbors[nodes[i]->getId()].insert(nodes[j]->getId());
            }
        }
    }

    std::unordered_map<Mesh::NodeId, std::vector<Mesh::NodeId>> neighbors;
    for (auto& [node_id, neighbor_set] : unique_neighbors)
    {
        neighbors[node_id] = std::vector<Mesh::NodeId>(neighbor_set.begin(), neighbor_set.end());
    }
    return neighbors;
}

std::unordered_map<Mesh::NodeId, double>
buildNodeControlVolumes(const std::vector<Mesh::ElementPtr>& elements)
{
    std::unordered_map<Mesh::NodeId, double> control_volumes;

    for (const auto& element : elements)
    {
        if (!element)
        {
            continue;
        }

        const auto& nodes = element->getNodes();
        if (nodes.empty())
        {
            continue;
        }

        const double measure = std::max(std::abs(element->getVolume()), kGeometryTolerance);
        const double share = measure / static_cast<double>(nodes.size());
        for (const auto& node : nodes)
        {
            if (node)
            {
                control_volumes[node->getId()] += share;
            }
        }
    }

    return control_volumes;
}

bool solveDenseLinearSystem(std::vector<std::vector<double>> matrix, std::vector<double> rhs,
                            std::vector<double>& solution)
{
    const std::size_t size = matrix.size();
    if (size == 0 || rhs.size() != size)
    {
        return false;
    }

    solution.assign(size, 0.0);
    for (std::size_t pivot = 0; pivot < size; ++pivot)
    {
        std::size_t best = pivot;
        for (std::size_t row = pivot + 1; row < size; ++row)
        {
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot]))
            {
                best = row;
            }
        }

        if (std::abs(matrix[best][pivot]) < kGeometryTolerance)
        {
            return false;
        }

        if (best != pivot)
        {
            std::swap(matrix[best], matrix[pivot]);
            std::swap(rhs[best], rhs[pivot]);
        }

        for (std::size_t row = pivot + 1; row < size; ++row)
        {
            const double factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (std::size_t col = pivot; col < size; ++col)
            {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    for (int row = static_cast<int>(size) - 1; row >= 0; --row)
    {
        double sum = rhs[static_cast<std::size_t>(row)];
        for (std::size_t col = static_cast<std::size_t>(row) + 1; col < size; ++col)
        {
            sum -= matrix[static_cast<std::size_t>(row)][col] * solution[col];
        }
        solution[static_cast<std::size_t>(row)] =
            sum / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }

    return true;
}

Mesh::NodePtr findNearestNode(const std::vector<Mesh::NodePtr>& nodes,
                              const Utils::Point3D& position)
{
    Mesh::NodePtr best_node = nullptr;
    double best_distance = std::numeric_limits<double>::max();

    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        const double distance = node->getPosition().distanceTo(position);
        if (distance < best_distance)
        {
            best_distance = distance;
            best_node = node;
        }
    }

    return best_node;
}

double coordinateByAxis(const Utils::Point3D& point, int axis)
{
    switch (axis)
    {
    case 0:
        return point.x();
    case 1:
        return point.y();
    default:
        return point.z();
    }
}

void setCoordinateByAxis(Utils::Point3D& point, int axis, double value)
{
    switch (axis)
    {
    case 0:
        point.setX(value);
        break;
    case 1:
        point.setY(value);
        break;
    default:
        point.setZ(value);
        break;
    }
}

Mesh::NodePtr findNodeById(const std::vector<Mesh::NodePtr>& nodes, Mesh::NodeId node_id)
{
    if (node_id < nodes.size() && nodes[node_id] && nodes[node_id]->getId() == node_id)
    {
        return nodes[node_id];
    }

    auto it = std::find_if(nodes.begin(), nodes.end(),
                           [node_id](const Mesh::NodePtr& node)
                           { return node && node->getId() == node_id; });
    return (it != nodes.end()) ? *it : nullptr;
}

int boundaryConditionRank(BoundaryConditionType type)
{
    return (type == BoundaryConditionType::DIRICHLET) ? 1 : 0;
}

std::vector<Mesh::NodeId> inferPeriodicPairsWithFallback(
    const Mesh::MeshPtr& mesh, const std::vector<Mesh::NodeId>& primary_nodes,
    const std::vector<Mesh::NodePtr>& primary_node_ptrs,
    const std::vector<Mesh::NodePtr>& opposite_face_nodes, double matching_tolerance);

std::vector<double> buildDistanceSignature(const Mesh::NodePtr& node,
                                           const std::vector<Mesh::NodePtr>& face_nodes)
{
    std::vector<double> signature;
    signature.reserve(face_nodes.size());
    for (const auto& other : face_nodes)
    {
        if (!node || !other || node->getId() == other->getId())
        {
            continue;
        }

        signature.push_back(node->getPosition().distanceTo(other->getPosition()));
    }

    std::sort(signature.begin(), signature.end());
    return signature;
}

double signatureDifference(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size())
    {
        return std::numeric_limits<double>::max();
    }

    double difference = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        difference += std::abs(lhs[i] - rhs[i]);
    }
    return difference;
}

bool inferPeriodicPairsByDistanceSignature(const std::vector<Mesh::NodePtr>& primary_nodes,
                                           const std::vector<Mesh::NodePtr>& opposite_nodes,
                                           double tolerance,
                                           std::vector<Mesh::NodeId>& inferred_pairs)
{
    if (primary_nodes.size() != opposite_nodes.size())
    {
        return false;
    }

    std::vector<std::vector<double>> primary_signatures;
    std::vector<std::vector<double>> opposite_signatures;
    primary_signatures.reserve(primary_nodes.size());
    opposite_signatures.reserve(opposite_nodes.size());

    for (const auto& node : primary_nodes)
    {
        primary_signatures.push_back(buildDistanceSignature(node, primary_nodes));
    }
    for (const auto& node : opposite_nodes)
    {
        opposite_signatures.push_back(buildDistanceSignature(node, opposite_nodes));
    }

    inferred_pairs.clear();
    inferred_pairs.reserve(primary_nodes.size());
    std::vector<bool> opposite_used(opposite_nodes.size(), false);
    std::vector<std::size_t> matched_indices(primary_nodes.size(), 0);

    for (std::size_t i = 0; i < primary_nodes.size(); ++i)
    {
        double best_error = std::numeric_limits<double>::max();
        std::size_t best_index = opposite_nodes.size();

        for (std::size_t j = 0; j < opposite_nodes.size(); ++j)
        {
            if (opposite_used[j])
            {
                continue;
            }

            const double error = signatureDifference(primary_signatures[i], opposite_signatures[j]);
            if (error < best_error)
            {
                best_error = error;
                best_index = j;
            }
        }

        if (best_index >= opposite_nodes.size() || best_error > tolerance)
        {
            inferred_pairs.clear();
            return false;
        }

        opposite_used[best_index] = true;
        matched_indices[i] = best_index;
        inferred_pairs.push_back(opposite_nodes[best_index]->getId());
    }

    for (std::size_t i = 0; i < primary_nodes.size(); ++i)
    {
        for (std::size_t j = i + 1; j < primary_nodes.size(); ++j)
        {
            const double primary_distance =
                primary_nodes[i]->getPosition().distanceTo(primary_nodes[j]->getPosition());
            const double opposite_distance =
                opposite_nodes[matched_indices[i]]->getPosition().distanceTo(
                    opposite_nodes[matched_indices[j]]->getPosition());
            if (std::abs(primary_distance - opposite_distance) > tolerance)
            {
                inferred_pairs.clear();
                return false;
            }
        }
    }

    return true;
}

std::vector<Mesh::NodeId>
inferPeriodicPairsAxisAligned(const Mesh::MeshPtr& mesh,
                              const std::vector<Mesh::NodeId>& primary_nodes)
{
    if (!mesh || primary_nodes.empty())
    {
        return {};
    }

    const auto& nodes = mesh->getNodes();
    std::vector<Mesh::NodePtr> primary_node_ptrs;
    primary_node_ptrs.reserve(primary_nodes.size());
    for (Mesh::NodeId node_id : primary_nodes)
    {
        const auto node = findNodeById(nodes, node_id);
        if (!node)
        {
            return {};
        }
        primary_node_ptrs.push_back(node);
    }

    std::array<double, 3> primary_min = {std::numeric_limits<double>::max(),
                                         std::numeric_limits<double>::max(),
                                         std::numeric_limits<double>::max()};
    std::array<double, 3> primary_max = {-std::numeric_limits<double>::max(),
                                         -std::numeric_limits<double>::max(),
                                         -std::numeric_limits<double>::max()};
    std::array<double, 3> mesh_min = {std::numeric_limits<double>::max(),
                                      std::numeric_limits<double>::max(),
                                      std::numeric_limits<double>::max()};
    std::array<double, 3> mesh_max = {-std::numeric_limits<double>::max(),
                                      -std::numeric_limits<double>::max(),
                                      -std::numeric_limits<double>::max()};

    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        for (int axis = 0; axis < 3; ++axis)
        {
            const double value = coordinateByAxis(node->getPosition(), axis);
            mesh_min[axis] = std::min(mesh_min[axis], value);
            mesh_max[axis] = std::max(mesh_max[axis], value);
        }
    }

    for (const auto& node : primary_node_ptrs)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            const double value = coordinateByAxis(node->getPosition(), axis);
            primary_min[axis] = std::min(primary_min[axis], value);
            primary_max[axis] = std::max(primary_max[axis], value);
        }
    }

    int periodic_axis = 0;
    double min_spread = primary_max[0] - primary_min[0];
    for (int axis = 1; axis < 3; ++axis)
    {
        const double spread = primary_max[axis] - primary_min[axis];
        if (spread < min_spread)
        {
            min_spread = spread;
            periodic_axis = axis;
        }
    }

    const double span = mesh_max[periodic_axis] - mesh_min[periodic_axis];
    if (span <= kGeometryTolerance)
    {
        return {};
    }

    double centroid_coordinate = 0.0;
    for (const auto& node : primary_node_ptrs)
    {
        centroid_coordinate += coordinateByAxis(node->getPosition(), periodic_axis);
    }
    centroid_coordinate /= static_cast<double>(primary_node_ptrs.size());

    const bool primary_is_min_face = std::abs(centroid_coordinate - mesh_min[periodic_axis]) <=
                                     std::abs(centroid_coordinate - mesh_max[periodic_axis]);
    const double target_axis_value =
        primary_is_min_face ? mesh_max[periodic_axis] : mesh_min[periodic_axis];
    const double translation = primary_is_min_face ? span : -span;
    const double axis_tolerance = std::max(1.0e-9, 1.0e-6 * std::abs(span));

    std::unordered_set<Mesh::NodeId> primary_node_set(primary_nodes.begin(), primary_nodes.end());
    std::unordered_set<Mesh::NodeId> used_targets;
    std::vector<Mesh::NodeId> inferred_pairs;
    inferred_pairs.reserve(primary_nodes.size());

    for (const auto& primary_node : primary_node_ptrs)
    {
        Utils::Point3D target_position = primary_node->getPosition();
        setCoordinateByAxis(target_position, periodic_axis,
                            coordinateByAxis(target_position, periodic_axis) + translation);

        Mesh::NodeId best_node_id = static_cast<Mesh::NodeId>(-1);
        double best_score = std::numeric_limits<double>::max();

        for (const auto& candidate : nodes)
        {
            if (!candidate || primary_node_set.count(candidate->getId()) > 0 ||
                used_targets.count(candidate->getId()) > 0)
            {
                continue;
            }

            const double axis_error = std::abs(
                coordinateByAxis(candidate->getPosition(), periodic_axis) - target_axis_value);
            if (axis_error > axis_tolerance)
            {
                continue;
            }

            double tangential_error_sq = 0.0;
            for (int axis = 0; axis < 3; ++axis)
            {
                if (axis == periodic_axis)
                {
                    continue;
                }

                const double delta = coordinateByAxis(candidate->getPosition(), axis) -
                                     coordinateByAxis(target_position, axis);
                tangential_error_sq += delta * delta;
            }

            const double score =
                tangential_error_sq + axis_error * axis_error / std::max(span * span, 1.0);
            if (score < best_score)
            {
                best_score = score;
                best_node_id = candidate->getId();
            }
        }

        if (best_node_id == static_cast<Mesh::NodeId>(-1))
        {
            return {};
        }

        used_targets.insert(best_node_id);
        inferred_pairs.push_back(best_node_id);
    }

    return inferred_pairs;
}

std::vector<Mesh::NodeId> inferPeriodicPairs(const Mesh::MeshPtr& mesh,
                                             const std::vector<Mesh::NodeId>& primary_nodes)
{
    if (!mesh || primary_nodes.empty())
    {
        return {};
    }

    const auto& nodes = mesh->getNodes();
    std::vector<Mesh::NodePtr> primary_node_ptrs;
    primary_node_ptrs.reserve(primary_nodes.size());
    for (Mesh::NodeId node_id : primary_nodes)
    {
        const auto node = findNodeById(nodes, node_id);
        if (!node)
        {
            return {};
        }
        primary_node_ptrs.push_back(node);
    }

    Utils::Point3D primary_centroid(0.0, 0.0, 0.0);
    for (const auto& node : primary_node_ptrs)
    {
        primary_centroid += node->getPosition();
    }
    primary_centroid /= static_cast<double>(primary_node_ptrs.size());

    Utils::Vector3D best_normal(0.0, 0.0, 0.0);
    double best_normal_magnitude_sq = 0.0;
    for (std::size_t i = 0; i < primary_node_ptrs.size(); ++i)
    {
        const Utils::Vector3D vi = primary_node_ptrs[i]->getPosition() - primary_centroid;
        for (std::size_t j = i + 1; j < primary_node_ptrs.size(); ++j)
        {
            const Utils::Vector3D vj = primary_node_ptrs[j]->getPosition() - primary_centroid;
            const Utils::Vector3D normal = vi.cross(vj);
            const double magnitude_sq = normal.magnitudeSquared();
            if (magnitude_sq > best_normal_magnitude_sq)
            {
                best_normal = normal;
                best_normal_magnitude_sq = magnitude_sq;
            }
        }
    }

    if (best_normal_magnitude_sq <= kGeometryTolerance * kGeometryTolerance)
    {
        return inferPeriodicPairsAxisAligned(mesh, primary_nodes);
    }

    best_normal /= std::sqrt(best_normal_magnitude_sq);

    std::unordered_set<Mesh::NodeId> primary_node_set(primary_nodes.begin(), primary_nodes.end());
    double best_distance = 0.0;
    bool found_opposite_plane = false;
    for (const auto& node : nodes)
    {
        if (!node || primary_node_set.count(node->getId()) > 0)
        {
            continue;
        }

        const double signed_distance = (node->getPosition() - primary_centroid).dot(best_normal);
        if (!found_opposite_plane || std::abs(signed_distance) > std::abs(best_distance))
        {
            best_distance = signed_distance;
            found_opposite_plane = true;
        }
    }

    if (!found_opposite_plane || std::abs(best_distance) <= kGeometryTolerance)
    {
        return inferPeriodicPairsAxisAligned(mesh, primary_nodes);
    }

    const double plane_tolerance = std::max(1.0e-8, 1.0e-5 * std::abs(best_distance));
    std::vector<Mesh::NodePtr> opposite_face_nodes;
    opposite_face_nodes.reserve(primary_nodes.size());
    for (const auto& node : nodes)
    {
        if (!node || primary_node_set.count(node->getId()) > 0)
        {
            continue;
        }

        const double signed_distance = (node->getPosition() - primary_centroid).dot(best_normal);
        if (std::abs(signed_distance - best_distance) <= plane_tolerance)
        {
            opposite_face_nodes.push_back(node);
        }
    }

    if (opposite_face_nodes.size() < primary_nodes.size())
    {
        return inferPeriodicPairsAxisAligned(mesh, primary_nodes);
    }

    Utils::Point3D opposite_centroid(0.0, 0.0, 0.0);
    for (const auto& node : opposite_face_nodes)
    {
        opposite_centroid += node->getPosition();
    }
    opposite_centroid /= static_cast<double>(opposite_face_nodes.size());

    const Utils::Vector3D translation = opposite_centroid - primary_centroid;
    double primary_face_diameter = 0.0;
    for (std::size_t i = 0; i < primary_node_ptrs.size(); ++i)
    {
        for (std::size_t j = i + 1; j < primary_node_ptrs.size(); ++j)
        {
            primary_face_diameter =
                std::max(primary_face_diameter, primary_node_ptrs[i]->getPosition().distanceTo(
                                                    primary_node_ptrs[j]->getPosition()));
        }
    }

    const double pairing_tolerance =
        std::max(1.0e-7, 1.0e-4 * std::max(primary_face_diameter, translation.magnitude()));
    std::unordered_set<Mesh::NodeId> used_targets;
    std::vector<Mesh::NodeId> inferred_pairs;
    inferred_pairs.reserve(primary_nodes.size());

    for (const auto& primary_node : primary_node_ptrs)
    {
        const Utils::Point3D target_position = primary_node->getPosition() + translation;
        Mesh::NodeId best_node_id = static_cast<Mesh::NodeId>(-1);
        double best_score = std::numeric_limits<double>::max();

        for (const auto& candidate : opposite_face_nodes)
        {
            if (!candidate || used_targets.count(candidate->getId()) > 0)
            {
                continue;
            }

            const double score = candidate->getPosition().distanceTo(target_position);
            if (score < best_score)
            {
                best_score = score;
                best_node_id = candidate->getId();
            }
        }

        if (best_node_id == static_cast<Mesh::NodeId>(-1) || best_score > pairing_tolerance)
        {
            return inferPeriodicPairsWithFallback(mesh, primary_nodes, primary_node_ptrs,
                                                  opposite_face_nodes, pairing_tolerance);
        }

        used_targets.insert(best_node_id);
        inferred_pairs.push_back(best_node_id);
    }

    return inferred_pairs;
}

bool inferPeriodicPairsNearestNeighbor(const std::vector<Mesh::NodePtr>& primary_nodes,
                                       const std::vector<Mesh::NodePtr>& opposite_nodes,
                                       std::vector<Mesh::NodeId>& inferred_pairs)
{
    if (primary_nodes.size() != opposite_nodes.size())
    {
        return false;
    }

    Utils::Point3D primary_centroid(0.0, 0.0, 0.0);
    for (const auto& n : primary_nodes) primary_centroid += n->getPosition();
    primary_centroid /= static_cast<double>(primary_nodes.size());

    Utils::Point3D opposite_centroid(0.0, 0.0, 0.0);
    for (const auto& n : opposite_nodes) opposite_centroid += n->getPosition();
    opposite_centroid /= static_cast<double>(opposite_nodes.size());

    Utils::Vector3D translation = opposite_centroid - primary_centroid;

    inferred_pairs.clear();
    inferred_pairs.reserve(primary_nodes.size());
    std::vector<bool> opposite_used(opposite_nodes.size(), false);

    for (const auto& p_node : primary_nodes)
    {
        Utils::Point3D target = p_node->getPosition() + translation;
        double best_dist = std::numeric_limits<double>::max();
        std::size_t best_index = opposite_nodes.size();

        for (std::size_t j = 0; j < opposite_nodes.size(); ++j)
        {
            if (opposite_used[j]) continue;
            double d = opposite_nodes[j]->getPosition().distanceTo(target);
            if (d < best_dist)
            {
                best_dist = d;
                best_index = j;
            }
        }

        if (best_index >= opposite_nodes.size()) return false;
        opposite_used[best_index] = true;
        inferred_pairs.push_back(opposite_nodes[best_index]->getId());
    }
    return true;
}

std::vector<Mesh::NodeId> inferPeriodicPairsWithFallback(
    const Mesh::MeshPtr& mesh, const std::vector<Mesh::NodeId>& primary_nodes,
    const std::vector<Mesh::NodePtr>& primary_node_ptrs,
    const std::vector<Mesh::NodePtr>& opposite_face_nodes, double matching_tolerance)
{
    std::vector<Mesh::NodeId> matched_pairs;
    if (inferPeriodicPairsByDistanceSignature(primary_node_ptrs, opposite_face_nodes,
                                              matching_tolerance, matched_pairs))
    {
        return matched_pairs;
    }

    if (inferPeriodicPairsNearestNeighbor(primary_node_ptrs, opposite_face_nodes, matched_pairs))
    {
        return matched_pairs;
    }

    return inferPeriodicPairsAxisAligned(mesh, primary_nodes);
}

} // namespace

MultiBoundaryElectricFieldSolver::MultiBoundaryElectricFieldSolver(Mesh::MeshPtr mesh)
    : mesh_(std::move(mesh)), convergence_tolerance_(1e-12), max_iterations_(10000),
      current_iteration_(0), use_adaptive_refinement_(false), use_multigrid_(false)
{
    initializeSolverConfig();
}

void MultiBoundaryElectricFieldSolver::initializeSolverConfig()
{
    solver_config_ = SolverConfiguration{};
    multigrid_config_ = MultigridConfiguration{};
    coupling_config_ = CouplingIterationConfiguration{};
    convergence_tolerance_ = solver_config_.tolerance;
    max_iterations_ = solver_config_.max_iterations;
}

void MultiBoundaryElectricFieldSolver::addBoundaryCondition(const std::string& boundary_name,
                                                            BoundaryConditionType type,
                                                            double value,
                                                            const std::vector<NodeId>& nodes)
{
    BoundaryCondition bc;
    bc.name = boundary_name;
    bc.type = type;
    bc.value = value;
    bc.nodes = nodes;
    bc.priority = static_cast<int>(boundary_conditions_.size());

    boundary_conditions_[boundary_name] = bc;
    for (NodeId node : nodes)
    {
        node_boundary_map_[node] = boundary_name;
    }
}

void MultiBoundaryElectricFieldSolver::addPeriodicBoundaryCondition(
    const std::string& boundary_name, const std::vector<NodeId>& primary_nodes,
    const std::vector<NodeId>& paired_nodes, double potential_offset)
{
    BoundaryCondition bc;
    bc.name = boundary_name;
    bc.type = BoundaryConditionType::PERIODIC;
    bc.value = potential_offset;
    bc.nodes = primary_nodes;
    bc.paired_nodes = paired_nodes;
    if (bc.paired_nodes.size() != bc.nodes.size())
    {
        const auto inferred_pairs = inferPeriodicPairs(mesh_, bc.nodes);
        if (inferred_pairs.size() == bc.nodes.size())
        {
            bc.paired_nodes = inferred_pairs;
        }
    }
    bc.priority = static_cast<int>(boundary_conditions_.size());

    boundary_conditions_[boundary_name] = bc;
    for (NodeId node : primary_nodes)
    {
        node_boundary_map_[node] = boundary_name;
    }
    for (NodeId node : paired_nodes)
    {
        node_boundary_map_[node] = boundary_name;
    }
}

void MultiBoundaryElectricFieldSolver::setChargeDensity(NodeId node_id, double charge_density)
{
    charge_densities_[node_id] = charge_density;
}

void MultiBoundaryElectricFieldSolver::setMaterialPermittivity(Mesh::MaterialId material_id,
                                                               double permittivity)
{
    material_permittivities_[material_id] = permittivity;
}

void MultiBoundaryElectricFieldSolver::setSolverConfiguration(const SolverConfiguration& config)
{
    solver_config_ = config;
    convergence_tolerance_ = solver_config_.tolerance;
    max_iterations_ = solver_config_.max_iterations;
    use_multigrid_ = (solver_config_.solver_type == SolverType::MULTIGRID);
}

void MultiBoundaryElectricFieldSolver::setMultigridConfiguration(
    const MultigridConfiguration& config)
{
    multigrid_config_ = config;
}

void MultiBoundaryElectricFieldSolver::setCouplingIterationConfiguration(
    const CouplingIterationConfiguration& config)
{
    coupling_config_.max_outer_iterations = std::max(1, config.max_outer_iterations);
    coupling_config_.potential_tolerance_v =
        std::max(1.0e-12, config.potential_tolerance_v);
    coupling_config_.relaxation = std::clamp(config.relaxation, 0.05, 1.0);
}

bool MultiBoundaryElectricFieldSolver::solve()
{
    residual_history_.clear();
    current_iteration_ = 0;

    if (!assembleSystem())
    {
        return false;
    }

    applyBoundaryConditions();

    const bool ok = use_multigrid_ ? solveWithMultigrid() : solveWithDirectMethod();
    if (ok)
    {
        calculateElectricField();
    }

    return ok;
}

bool MultiBoundaryElectricFieldSolver::solveCoupled(const CouplingCallback& callback)
{
    coupling_iteration_count_ = 0;
    coupling_residual_history_.clear();

    const int max_outer_iterations = std::max(1, coupling_config_.max_outer_iterations);
    const double tolerance_v = std::max(1.0e-12, coupling_config_.potential_tolerance_v);
    const double relaxation = std::clamp(coupling_config_.relaxation, 0.05, 1.0);

    auto previous_potentials = snapshotNodePotentials();
    for (int outer = 0; outer < max_outer_iterations; ++outer)
    {
        if (callback)
        {
            callback(outer, previous_potentials);
        }

        if (!solve())
        {
            return false;
        }

        auto current_potentials = snapshotNodePotentials();
        if (current_potentials.size() != previous_potentials.size())
        {
            return false;
        }

        double max_delta = 0.0;
        for (std::size_t i = 0; i < current_potentials.size(); ++i)
        {
            max_delta = std::max(max_delta, std::abs(current_potentials[i] - previous_potentials[i]));
        }

        coupling_iteration_count_ = outer + 1;
        coupling_residual_history_.push_back(max_delta);

        if (max_delta <= tolerance_v)
        {
            return true;
        }

        if (outer + 1 < max_outer_iterations && relaxation < 1.0)
        {
            applyRelaxedNodePotentials(previous_potentials, current_potentials, relaxation);
            current_potentials = snapshotNodePotentials();
        }
        previous_potentials = std::move(current_potentials);
    }

    return !coupling_residual_history_.empty() &&
           coupling_residual_history_.back() <= tolerance_v;
}

std::vector<double> MultiBoundaryElectricFieldSolver::snapshotNodePotentials() const
{
    std::vector<double> values;
    if (!mesh_)
    {
        return values;
    }

    const auto& nodes = mesh_->getNodes();
    values.reserve(nodes.size());
    for (const auto& node : nodes)
    {
        if (node)
        {
            values.push_back(node->getPotential());
        }
    }
    return values;
}

void MultiBoundaryElectricFieldSolver::applyRelaxedNodePotentials(
    const std::vector<double>& previous, const std::vector<double>& current, double relaxation)
{
    if (!mesh_ || previous.size() != current.size())
    {
        return;
    }

    const auto& nodes = mesh_->getNodes();
    std::size_t index = 0;
    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }
        const double relaxed = previous[index] + relaxation * (current[index] - previous[index]);
        node->setPotential(relaxed);
        ++index;
    }
}

bool MultiBoundaryElectricFieldSolver::assembleSystem()
{
    if (!mesh_ || mesh_->getNodes().empty() || mesh_->getElements().empty())
    {
        return false;
    }

    spatial_index_.reset();
    electric_field_cache_.clear();
    node_index_map_.clear();
    index_to_node_.clear();

    const auto& nodes = mesh_->getNodes();
    index_to_node_.reserve(nodes.size());
    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }
        node_index_map_[node->getId()] = index_to_node_.size();
        index_to_node_.push_back(node->getId());
    }

    const size_t n = index_to_node_.size();
    dense_system_matrix_.assign(n, std::vector<double>(n, 0.0));
    rhs_storage_.assign(n, 0.0);
    solution_storage_.assign(n, 0.0);
    preconditioner_diagonal_.assign(n, 1.0);

    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }
        const size_t idx = getGlobalIndex(node->getId());
        if (idx < solution_storage_.size())
        {
            solution_storage_[idx] = node->getPotential();
        }
    }

#ifdef USE_EIGEN
    const Eigen::Index en = static_cast<Eigen::Index>(n);
    system_matrix_.resize(en, en);
    rhs_vector_.resize(en);
    solution_vector_.resize(en);
    rhs_vector_.setZero();
    solution_vector_.setZero();
#else
    system_matrix_.resize(static_cast<int>(n), static_cast<int>(n));
    rhs_vector_.assign(n, 0.0);
    solution_vector_.assign(n, 0.0);
#endif

    for (const auto& element : mesh_->getElements())
    {
        assembleElement(element);
    }

    return true;
}

void MultiBoundaryElectricFieldSolver::assembleElement(Mesh::ElementPtr element)
{
    if (!element)
    {
        return;
    }

    if (element->getType() == Mesh::ElementType::TETRAHEDRON)
    {
        assembleTetrahedronElement(element);
    }
    else
    {
        assembleHexahedronElement(element);
    }
}

void MultiBoundaryElectricFieldSolver::applyBoundaryConditions()
{
    std::vector<std::pair<std::string, BoundaryCondition>> ordered;
    ordered.reserve(boundary_conditions_.size());
    for (const auto& kv : boundary_conditions_)
    {
        ordered.push_back(kv);
    }

    std::sort(ordered.begin(), ordered.end(),
              [](const auto& lhs, const auto& rhs)
              {
                  if (boundaryConditionRank(lhs.second.type) !=
                      boundaryConditionRank(rhs.second.type))
                  {
                      return boundaryConditionRank(lhs.second.type) <
                             boundaryConditionRank(rhs.second.type);
                  }
                  return lhs.second.priority < rhs.second.priority;
              });

    for (const auto& kv : ordered)
    {
        applyBoundaryCondition(kv.second);
    }
}

void MultiBoundaryElectricFieldSolver::applyBoundaryCondition(const BoundaryCondition& bc)
{
    if (bc.type == BoundaryConditionType::PERIODIC)
    {
        const std::size_t pair_count = std::min(bc.nodes.size(), bc.paired_nodes.size());
        for (std::size_t i = 0; i < pair_count; ++i)
        {
            const size_t p_idx = getGlobalIndex(bc.nodes[i]);
            const size_t q_idx = getGlobalIndex(bc.paired_nodes[i]);
            if (p_idx >= dense_system_matrix_.size() || q_idx >= dense_system_matrix_.size() || p_idx == q_idx)
            {
                continue;
            }

            const double offset = bc.value;

            for (size_t row = 0; row < dense_system_matrix_.size(); ++row) {
                rhs_storage_[row] -= dense_system_matrix_[row][q_idx] * offset;
            }
            
            for (size_t col = 0; col < dense_system_matrix_.size(); ++col) {
                dense_system_matrix_[p_idx][col] += dense_system_matrix_[q_idx][col];
            }
            rhs_storage_[p_idx] += rhs_storage_[q_idx]; 
            
            for (size_t row = 0; row < dense_system_matrix_.size(); ++row) {
                if (row != q_idx) {
                    dense_system_matrix_[row][p_idx] += dense_system_matrix_[row][q_idx];
                }
            }
            
            std::fill(dense_system_matrix_[q_idx].begin(), dense_system_matrix_[q_idx].end(), 0.0);
            for (size_t row = 0; row < dense_system_matrix_.size(); ++row) {
                dense_system_matrix_[row][q_idx] = 0.0;
            }
            
            dense_system_matrix_[q_idx][q_idx] = 1.0;
            rhs_storage_[q_idx] = 0.0;
            solution_storage_[q_idx] = 0.0;
        }
        return;
    }

    if (bc.type == BoundaryConditionType::FLOATING)
    {
        if (bc.nodes.empty())
        {
            return;
        }

        const NodeId anchor_node = bc.nodes.front();
        const size_t p_idx = getGlobalIndex(anchor_node);
        if (p_idx >= dense_system_matrix_.size())
        {
            return;
        }

        for (std::size_t i = 1; i < bc.nodes.size(); ++i)
        {
            const size_t q_idx = getGlobalIndex(bc.nodes[i]);
            if (q_idx >= dense_system_matrix_.size() || p_idx == q_idx)
            {
                continue;
            }

            for (size_t col = 0; col < dense_system_matrix_.size(); ++col) {
                dense_system_matrix_[p_idx][col] += dense_system_matrix_[q_idx][col];
            }
            rhs_storage_[p_idx] += rhs_storage_[q_idx];

            for (size_t row = 0; row < dense_system_matrix_.size(); ++row) {
                if (row != q_idx) {
                    dense_system_matrix_[row][p_idx] += dense_system_matrix_[row][q_idx];
                }
            }

            std::fill(dense_system_matrix_[q_idx].begin(), dense_system_matrix_[q_idx].end(), 0.0);
            for (size_t row = 0; row < dense_system_matrix_.size(); ++row) {
                dense_system_matrix_[row][q_idx] = 0.0;
            }

            dense_system_matrix_[q_idx][q_idx] = 1.0;
            rhs_storage_[q_idx] = 0.0;
            solution_storage_[q_idx] = 0.0;
        }

        rhs_storage_[p_idx] += bc.value;
        return;
    }

    for (NodeId node : bc.nodes)
    {
        const size_t idx = getGlobalIndex(node);
        if (idx >= dense_system_matrix_.size())
        {
            continue;
        }

        switch (bc.type)
        {
        case BoundaryConditionType::DIRICHLET:
            for (std::size_t row = 0; row < dense_system_matrix_.size(); ++row)
            {
                if (row != idx)
                {
                    rhs_storage_[row] -= dense_system_matrix_[row][idx] * bc.value;
                    dense_system_matrix_[row][idx] = 0.0;
                }
            }

            std::fill(dense_system_matrix_[idx].begin(), dense_system_matrix_[idx].end(), 0.0);
            dense_system_matrix_[idx][idx] = 1.0;
            rhs_storage_[idx] = bc.value;
            solution_storage_[idx] = bc.value;
            break;

        case BoundaryConditionType::NEUMANN:
            rhs_storage_[idx] += bc.value;
            break;

        case BoundaryConditionType::ROBIN:
        case BoundaryConditionType::MIXED:
            dense_system_matrix_[idx][idx] += std::max(bc.alpha, kGeometryTolerance);
            rhs_storage_[idx] += bc.value + bc.beta;
            break;
        }
    }
}

double MultiBoundaryElectricFieldSolver::getPotential(NodeId node_id) const
{
    const size_t idx = getGlobalIndex(node_id);
    if (idx < solution_storage_.size())
    {
        return solution_storage_[idx];
    }

    if (mesh_)
    {
        const auto& nodes = mesh_->getNodes();
        if (idx < nodes.size() && nodes[idx])
        {
            return nodes[idx]->getPotential();
        }
    }

    return 0.0;
}

Utils::Vector3D MultiBoundaryElectricFieldSolver::getElectricField(NodeId node_id) const
{
    const auto it = electric_field_cache_.find(node_id);
    if (it != electric_field_cache_.end())
    {
        return it->second;
    }

    if (mesh_)
    {
        const auto& nodes = mesh_->getNodes();
        if (const size_t idx = getGlobalIndex(node_id); idx < nodes.size() && nodes[idx])
        {
            return nodes[idx]->getElectricField();
        }
    }

    return Utils::Vector3D(0.0, 0.0, 0.0);
}

Utils::Vector3D
MultiBoundaryElectricFieldSolver::getElectricField(const Utils::Point3D& position) const
{
    if (!mesh_)
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    if (!spatial_index_)
    {
        spatial_index_ = std::make_shared<Mesh::SpatialIndex>();
        spatial_index_->initialize(mesh_->getElements());
    }

    if (spatial_index_)
    {
        const auto candidates = spatial_index_->findElementsContaining(position);
        for (const auto& element : candidates)
        {
            if (!element || !element->contains(position))
            {
                continue;
            }

            Utils::Vector3D accumulated(0.0, 0.0, 0.0);
            double count = 0.0;
            for (const auto& node : element->getNodes())
            {
                if (!node)
                {
                    continue;
                }

                accumulated = accumulated + getElectricField(node->getId());
                count += 1.0;
            }

            if (count > 0.0)
            {
                return accumulated / count;
            }
        }
    }

    if (const auto nearest = findNearestNode(mesh_->getNodes(), position))
    {
        return getElectricField(nearest->getId());
    }

    return Utils::Vector3D(0.0, 0.0, 0.0);
}

double MultiBoundaryElectricFieldSolver::validateSolution() const
{
    return residual_history_.empty() ? 0.0 : residual_history_.back();
}

double MultiBoundaryElectricFieldSolver::calculateElectrostaticEnergy() const
{
    if (!mesh_)
    {
        return 0.0;
    }

    const auto control_volumes = buildNodeControlVolumes(mesh_->getElements());
    double total_energy = 0.0;
    for (const auto& node : mesh_->getNodes())
    {
        if (!node)
        {
            continue;
        }

        const double control_volume = [&, node]()
        {
            const auto it = control_volumes.find(node->getId());
            return (it != control_volumes.end()) ? it->second : 0.0;
        }();

        const Utils::Vector3D electric_field = getElectricField(node->getId());
        total_energy +=
            0.5 * kVacuumPermittivity * electric_field.magnitudeSquared() * control_volume;
    }

    return total_energy;
}

bool MultiBoundaryElectricFieldSolver::solveWithMultigrid()
{
    setupPreconditioner();

    const int max_cycles = std::max(1, std::min(max_iterations_, 200));
    for (current_iteration_ = 0; current_iteration_ < max_cycles; ++current_iteration_)
    {
        performVCycle();
        const double residual = calculateResidualNorm();
        residual_history_.push_back(residual);
        if (residual <= convergence_tolerance_)
        {
            ++current_iteration_;
            return true;
        }
    }

    return !residual_history_.empty() && std::isfinite(residual_history_.back());
}

bool MultiBoundaryElectricFieldSolver::solveWithDirectMethod()
{
    setupPreconditioner();

    bool ok = false;
    if (solver_config_.solver_type == SolverType::DIRECT_LU)
    {
        ok = solveDenseLinearSystem(dense_system_matrix_, rhs_storage_, solution_storage_);
        current_iteration_ = ok ? 1 : 0;
        residual_history_.push_back(calculateResidualNorm());
    }
    else
    {
        ok = conjugateGradientSolve();
        if (residual_history_.empty())
        {
            residual_history_.push_back(calculateResidualNorm());
        }
    }

#ifndef USE_EIGEN
    rhs_vector_ = rhs_storage_;
    solution_vector_ = solution_storage_;
#endif
    return ok;
}

bool MultiBoundaryElectricFieldSolver::conjugateGradientSolve()
{
    const std::size_t size = rhs_storage_.size();
    if (size == 0 || dense_system_matrix_.size() != size)
    {
        return false;
    }

    std::vector<double> x = solution_storage_;
    if (x.size() != size)
    {
        x.assign(size, 0.0);
    }

    std::vector<double> r = rhs_storage_;
    const std::vector<double> ax0 = multiplyMatrix(dense_system_matrix_, x);
    for (std::size_t i = 0; i < size; ++i)
    {
        r[i] -= ax0[i];
    }

    std::vector<double> z(size, 0.0);
    for (std::size_t i = 0; i < size; ++i)
    {
        z[i] = r[i] / std::max(preconditioner_diagonal_[i], kGeometryTolerance);
    }

    std::vector<double> p = z;
    double rz_old = dotProduct(r, z);
    if (std::sqrt(std::max(rz_old, 0.0)) <= convergence_tolerance_)
    {
        solution_storage_ = x;
        return true;
    }

    for (int iteration = 0; iteration < max_iterations_; ++iteration)
    {
        const std::vector<double> ap = multiplyMatrix(dense_system_matrix_, p);
        const double denominator = dotProduct(p, ap);
        if (std::abs(denominator) < kBreakdownTolerance)
        {
            break;
        }

        const double alpha = rz_old / denominator;
        for (std::size_t i = 0; i < size; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }

        const double residual_norm = std::sqrt(std::max(dotProduct(r, r), 0.0));
        residual_history_.push_back(residual_norm);
        current_iteration_ = iteration + 1;
        if (residual_norm <= convergence_tolerance_)
        {
            solution_storage_ = x;
            return true;
        }

        for (std::size_t i = 0; i < size; ++i)
        {
            z[i] = r[i] / std::max(preconditioner_diagonal_[i], kGeometryTolerance);
        }

        const double rz_new = dotProduct(r, z);
        if (std::abs(rz_old) < kBreakdownTolerance)
        {
            break;
        }

        const double beta = rz_new / rz_old;
        for (std::size_t i = 0; i < size; ++i)
        {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }

    solution_storage_ = x;
    return false;
}

void MultiBoundaryElectricFieldSolver::performVCycle()
{
    if (dense_system_matrix_.empty() || rhs_storage_.empty())
    {
        return;
    }

    setupPreconditioner();

    const int sweeps =
        std::max(1, multigrid_config_.pre_smoothing_steps + multigrid_config_.post_smoothing_steps);
    const double omega = 2.0 / 3.0;

    if (solution_storage_.size() != rhs_storage_.size())
    {
        solution_storage_.assign(rhs_storage_.size(), 0.0);
    }

    for (int sweep = 0; sweep < sweeps; ++sweep)
    {
        const std::vector<double> ax = multiplyMatrix(dense_system_matrix_, solution_storage_);
        for (std::size_t i = 0; i < solution_storage_.size(); ++i)
        {
            const double residual = rhs_storage_[i] - ax[i];
            solution_storage_[i] +=
                omega * residual / std::max(preconditioner_diagonal_[i], kGeometryTolerance);
        }
    }
}

void MultiBoundaryElectricFieldSolver::calculateElectricField()
{
    if (!mesh_)
    {
        return;
    }

    for (const auto& bc_pair : boundary_conditions_) {
        const auto& bc = bc_pair.second;
        if (bc.type == BoundaryConditionType::FLOATING && bc.nodes.size() > 1) {
            const size_t p_idx = getGlobalIndex(bc.nodes.front());
            if (p_idx < solution_storage_.size()) {
                const double anchor_pot = solution_storage_[p_idx];
                for (size_t i = 1; i < bc.nodes.size(); ++i) {
                    const size_t q_idx = getGlobalIndex(bc.nodes[i]);
                    if (q_idx < solution_storage_.size()) {
                        solution_storage_[q_idx] = anchor_pot;
                    }
                }
            }
        }
        else if (bc.type == BoundaryConditionType::PERIODIC) {
            const size_t pair_count = std::min(bc.nodes.size(), bc.paired_nodes.size());
            for (size_t i = 0; i < pair_count; ++i) {
                const size_t p_idx = getGlobalIndex(bc.nodes[i]);
                const size_t q_idx = getGlobalIndex(bc.paired_nodes[i]);
                if (p_idx < solution_storage_.size() && q_idx < solution_storage_.size()) {
                    solution_storage_[q_idx] = solution_storage_[p_idx] + bc.value;
                }
            }
        }
    }

    const auto& nodes = mesh_->getNodes();
    const auto neighbors = buildNodeNeighbors(mesh_->getElements());
    electric_field_cache_.clear();

    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        const size_t idx = getGlobalIndex(node->getId());
        if (idx < solution_storage_.size())
        {
            node->setPotential(solution_storage_[idx]);
        }
    }

    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        std::array<std::array<double, 3>, 3> normal_matrix{};
        std::array<double, 3> rhs{};
        bool have_geometry = false;

        if (const auto it = neighbors.find(node->getId()); it != neighbors.end())
        {
            for (NodeId neighbor_id : it->second)
            {
                const size_t center_idx = getGlobalIndex(node->getId());
                const size_t neighbor_idx = getGlobalIndex(neighbor_id);
                if (center_idx >= solution_storage_.size() ||
                    neighbor_idx >= solution_storage_.size())
                {
                    continue;
                }

                const auto& neighbor = mesh_->getNodes()[neighbor_idx];
                if (!neighbor)
                {
                    continue;
                }

                const Utils::Vector3D delta = neighbor->getPosition() - node->getPosition();
                if (delta.magnitudeSquared() <= kGeometryTolerance * kGeometryTolerance)
                {
                    continue;
                }

                const double delta_phi =
                    solution_storage_[neighbor_idx] - solution_storage_[center_idx];

                normal_matrix[0][0] += delta.x() * delta.x();
                normal_matrix[0][1] += delta.x() * delta.y();
                normal_matrix[0][2] += delta.x() * delta.z();
                normal_matrix[1][0] += delta.y() * delta.x();
                normal_matrix[1][1] += delta.y() * delta.y();
                normal_matrix[1][2] += delta.y() * delta.z();
                normal_matrix[2][0] += delta.z() * delta.x();
                normal_matrix[2][1] += delta.z() * delta.y();
                normal_matrix[2][2] += delta.z() * delta.z();

                rhs[0] += delta.x() * delta_phi;
                rhs[1] += delta.y() * delta_phi;
                rhs[2] += delta.z() * delta_phi;
                have_geometry = true;
            }
        }

        std::array<double, 3> gradient{};
        Utils::Vector3D electric_field(0.0, 0.0, 0.0);
        if (have_geometry && solve3x3(normal_matrix, rhs, gradient))
        {
            electric_field = Utils::Vector3D(-gradient[0], -gradient[1], -gradient[2]);
        }

        node->setElectricField(electric_field);
        electric_field_cache_[node->getId()] = electric_field;
    }

    if (!spatial_index_)
    {
        spatial_index_ = std::make_shared<Mesh::SpatialIndex>();
        spatial_index_->initialize(mesh_->getElements());
    }

#ifndef USE_EIGEN
    solution_vector_ = solution_storage_;
#endif
}

double MultiBoundaryElectricFieldSolver::calculateResidualNorm() const
{
    if (dense_system_matrix_.empty() || rhs_storage_.empty() ||
        solution_storage_.size() != rhs_storage_.size())
    {
        return std::numeric_limits<double>::infinity();
    }

    const std::vector<double> ax = multiplyMatrix(dense_system_matrix_, solution_storage_);
    double sum = 0.0;
    for (std::size_t i = 0; i < rhs_storage_.size(); ++i)
    {
        const double residual = rhs_storage_[i] - ax[i];
        sum += residual * residual;
    }
    return std::sqrt(sum);
}

void MultiBoundaryElectricFieldSolver::setupPreconditioner()
{
    preconditioner_diagonal_.assign(dense_system_matrix_.size(), 1.0);
    for (std::size_t i = 0; i < dense_system_matrix_.size(); ++i)
    {
        if (i < dense_system_matrix_[i].size() &&
            std::abs(dense_system_matrix_[i][i]) > kGeometryTolerance)
        {
            preconditioner_diagonal_[i] = dense_system_matrix_[i][i];
        }
    }
}

void MultiBoundaryElectricFieldSolver::assembleTetrahedronElement(Mesh::ElementPtr element)
{
    auto tetra = std::dynamic_pointer_cast<Mesh::TetrahedronElement>(element);
    if (!tetra)
    {
        assembleHexahedronElement(element);
        return;
    }

    const auto geometry = calculateTetrahedronGeometry(tetra);
    if (geometry.volume <= kGeometryTolerance)
    {
        return;
    }

    const auto gradients = calculateShapeFunctionGradients(geometry);
    const auto& nodes = tetra->getNodes();
    const double permittivity = getMaterialPermittivity(element);

    for (std::size_t i = 0; i < nodes.size() && i < gradients.size(); ++i)
    {
        if (!nodes[i])
        {
            continue;
        }

        const size_t global_i = getGlobalIndex(nodes[i]->getId());
        if (global_i >= dense_system_matrix_.size())
        {
            continue;
        }

        for (std::size_t j = 0; j < nodes.size() && j < gradients.size(); ++j)
        {
            if (!nodes[j])
            {
                continue;
            }

            const size_t global_j = getGlobalIndex(nodes[j]->getId());
            if (global_j >= dense_system_matrix_[global_i].size())
            {
                continue;
            }

            const double kij = permittivity * gradients[i].dot(gradients[j]) * geometry.volume;
            dense_system_matrix_[global_i][global_j] += kij;
        }

        if (const auto rho = charge_densities_.find(nodes[i]->getId());
            rho != charge_densities_.end())
        {
            rhs_storage_[global_i] += -rho->second * geometry.volume / 4.0;
        }
    }
}

void MultiBoundaryElectricFieldSolver::assembleHexahedronElement(Mesh::ElementPtr element)
{
    if (!element)
    {
        return;
    }

    const auto& nodes = element->getNodes();
    if (nodes.size() < 2)
    {
        return;
    }

    const double measure = std::max(std::abs(element->getVolume()), kGeometryTolerance);
    const double permittivity = getMaterialPermittivity(element);
    const double pair_count = static_cast<double>((nodes.size() * (nodes.size() - 1)) / 2);
    const double edge_scale = measure / std::max(1.0, pair_count);

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        if (!nodes[i])
        {
            continue;
        }

        const size_t global_i = getGlobalIndex(nodes[i]->getId());
        if (global_i >= dense_system_matrix_.size())
        {
            continue;
        }

        for (std::size_t j = i + 1; j < nodes.size(); ++j)
        {
            if (!nodes[j])
            {
                continue;
            }

            const size_t global_j = getGlobalIndex(nodes[j]->getId());
            if (global_j >= dense_system_matrix_.size())
            {
                continue;
            }

            const double distance = std::max(
                nodes[i]->getPosition().distanceTo(nodes[j]->getPosition()), kGeometryTolerance);
            const double weight = permittivity * edge_scale / (distance * distance);

            dense_system_matrix_[global_i][global_i] += weight;
            dense_system_matrix_[global_j][global_j] += weight;
            dense_system_matrix_[global_i][global_j] -= weight;
            dense_system_matrix_[global_j][global_i] -= weight;
        }

        if (const auto rho = charge_densities_.find(nodes[i]->getId());
            rho != charge_densities_.end())
        {
            rhs_storage_[global_i] += -rho->second * measure / static_cast<double>(nodes.size());
        }
    }
}

ElementGeometry MultiBoundaryElectricFieldSolver::calculateTetrahedronGeometry(
    Mesh::TetrahedronPtr tetrahedron) const
{
    ElementGeometry geometry;

    if (!tetrahedron)
    {
        return geometry;
    }

    const auto& nodes = tetrahedron->getNodes();
    if (nodes.size() != 4 || !nodes[0] || !nodes[1] || !nodes[2] || !nodes[3])
    {
        return geometry;
    }

    for (std::size_t i = 0; i < 4; ++i)
    {
        geometry.tetra_nodes[i] = nodes[i]->getPosition();
    }

    const auto& p0 = geometry.tetra_nodes[0];
    const auto& p1 = geometry.tetra_nodes[1];
    const auto& p2 = geometry.tetra_nodes[2];
    const auto& p3 = geometry.tetra_nodes[3];

    const Utils::Vector3D a = p1 - p0;
    const Utils::Vector3D b = p2 - p0;
    const Utils::Vector3D c = p3 - p0;

    geometry.volume = std::abs(a.dot(b.cross(c))) / 6.0;
    geometry.centroid = Utils::Point3D((p0.x() + p1.x() + p2.x() + p3.x()) / 4.0,
                                       (p0.y() + p1.y() + p2.y() + p3.y()) / 4.0,
                                       (p0.z() + p1.z() + p2.z() + p3.z()) / 4.0);

    geometry.face_normals = {
        (p2 - p1).cross(p3 - p1),
        (p3 - p0).cross(p2 - p0),
        (p1 - p0).cross(p3 - p0),
        (p2 - p0).cross(p1 - p0),
    };

    geometry.face_areas.reserve(geometry.face_normals.size());
    for (const auto& normal : geometry.face_normals)
    {
        geometry.face_areas.push_back(0.5 * normal.magnitude());
    }

    return geometry;
}

std::vector<Utils::Vector3D> MultiBoundaryElectricFieldSolver::calculateShapeFunctionGradients(
    const ElementGeometry& geometry) const
{
    if (geometry.volume <= kGeometryTolerance)
    {
        return std::vector<Utils::Vector3D>(4, Utils::Vector3D(0.0, 0.0, 0.0));
    }

    const auto& p0 = geometry.tetra_nodes[0];
    const auto& p1 = geometry.tetra_nodes[1];
    const auto& p2 = geometry.tetra_nodes[2];
    const auto& p3 = geometry.tetra_nodes[3];
    const double scale = 1.0 / (6.0 * geometry.volume);

    return {
        ((p3 - p1).cross(p2 - p1)) * scale,
        ((p2 - p0).cross(p3 - p0)) * scale,
        ((p3 - p0).cross(p1 - p0)) * scale,
        ((p1 - p0).cross(p2 - p0)) * scale,
    };
}

std::pair<Eigen::Matrix3d, std::vector<Utils::Vector3D>>
MultiBoundaryElectricFieldSolver::calculateHexahedronShapeGradients(Mesh::TetrahedronPtr hexahedron,
                                                                    double /*xi*/, double /*eta*/,
                                                                    double /*zeta*/) const
{
    const auto geometry = calculateTetrahedronGeometry(hexahedron);

#ifdef USE_EIGEN
    Eigen::Matrix3d jacobian = Eigen::Matrix3d::Identity();
#else
    Eigen::Matrix3d jacobian{};
    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            jacobian[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
#endif

    if (geometry.volume > kGeometryTolerance)
    {
        const auto& p0 = geometry.tetra_nodes[0];
        const Utils::Vector3D e1 = geometry.tetra_nodes[1] - p0;
        const Utils::Vector3D e2 = geometry.tetra_nodes[2] - p0;
        const Utils::Vector3D e3 = geometry.tetra_nodes[3] - p0;

#ifdef USE_EIGEN
        jacobian(0, 0) = e1.x();
        jacobian(1, 0) = e1.y();
        jacobian(2, 0) = e1.z();
        jacobian(0, 1) = e2.x();
        jacobian(1, 1) = e2.y();
        jacobian(2, 1) = e2.z();
        jacobian(0, 2) = e3.x();
        jacobian(1, 2) = e3.y();
        jacobian(2, 2) = e3.z();
#else
        jacobian[0][0] = e1.x();
        jacobian[1][0] = e1.y();
        jacobian[2][0] = e1.z();
        jacobian[0][1] = e2.x();
        jacobian[1][1] = e2.y();
        jacobian[2][1] = e2.z();
        jacobian[0][2] = e3.x();
        jacobian[1][2] = e3.y();
        jacobian[2][2] = e3.z();
#endif
    }

    return {jacobian, calculateShapeFunctionGradients(geometry)};
}

size_t MultiBoundaryElectricFieldSolver::getGlobalIndex(NodeId node_id) const
{
    const auto it = node_index_map_.find(node_id);
    if (it != node_index_map_.end())
    {
        return it->second;
    }
    return static_cast<size_t>(node_id);
}

double MultiBoundaryElectricFieldSolver::getMaterialPermittivity(Mesh::ElementPtr element) const
{
    if (!element)
    {
        return kVacuumPermittivity;
    }

    if (const auto it = material_permittivities_.find(element->getMaterialId());
        it != material_permittivities_.end())
    {
        return resolvePermittivity(it->second);
    }

    return kVacuumPermittivity;
}

} // namespace FieldSolver
} // namespace SCDAT
