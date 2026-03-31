#include "../include/PoissonSolver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_set>

namespace SCDAT
{
namespace FieldSolver
{

namespace
{

constexpr double kVacuumPermittivity = 8.8541878128e-12;
constexpr double kGeometryTolerance = 1.0e-12;

double resolvePermittivity(double value)
{
    if (value <= 0.0)
    {
        return kVacuumPermittivity;
    }

    return (value < 1.0e-8) ? value : value * kVacuumPermittivity;
}

double elementMeasure(const Mesh::ElementPtr& element)
{
    if (!element)
    {
        return 0.0;
    }
    return std::max(std::abs(element->getVolume()), kGeometryTolerance);
}

Mesh::NodePtr findNodeById(const std::vector<Mesh::NodePtr>& nodes, Mesh::NodeId node_id)
{
    if (node_id < nodes.size() && nodes[node_id] && nodes[node_id]->getId() == node_id)
    {
        return nodes[node_id];
    }

    auto it = std::find_if(nodes.begin(), nodes.end(), [node_id](const Mesh::NodePtr& node)
                           { return node && node->getId() == node_id; });
    return (it != nodes.end()) ? *it : nullptr;
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
    neighbors.reserve(unique_neighbors.size());
    for (auto& [node_id, neighbor_set] : unique_neighbors)
    {
        neighbors[node_id] =
            std::vector<Mesh::NodeId>(neighbor_set.begin(), neighbor_set.end());
    }
    return neighbors;
}

std::unordered_map<Mesh::NodeId, double>
buildNodalControlVolumes(const std::vector<Mesh::ElementPtr>& elements)
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

        const double nodal_share = elementMeasure(element) / static_cast<double>(nodes.size());
        for (const auto& node : nodes)
        {
            if (node)
            {
                control_volumes[node->getId()] += nodal_share;
            }
        }
    }

    return control_volumes;
}

bool solve3x3(std::array<std::array<double, 3>, 3> matrix, std::array<double, 3> rhs,
              std::array<double, 3>& solution)
{
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

} // namespace

double BoundaryCondition::getValue(const Utils::Point3D& position, double time) const
{
    if (is_time_dependent_ && time_function_)
    {
        return time_function_(position, time);
    }
    return value_;
}

void BoundaryCondition::apply(Mesh::NodePtr node, double time) const
{
    if (!node || !isApplicable(node))
    {
        return;
    }
    if (type_ == BoundaryConditionType::DIRICHLET)
    {
        node->setPotential(getValue(node->getPosition(), time));
    }
}

bool BoundaryCondition::isApplicable(Mesh::NodePtr node) const
{
    return node && node->getBoundaryType() != Mesh::BoundaryType::INTERIOR;
}

std::string BoundaryCondition::toString() const
{
    std::ostringstream stream;
    stream << "BoundaryCondition(type=" << static_cast<int>(type_) << ", value=" << value_
           << ", time_dependent=" << (is_time_dependent_ ? "true" : "false") << ")";
    return stream.str();
}

void PoissonSolver::initialize()
{
    buildNodeIndexMapping();
    system_ = std::make_shared<LinearSystem>(node_to_index_.size());
}

void PoissonSolver::addBoundaryCondition(Mesh::NodeId node_id, BoundaryConditionPtr bc)
{
    boundary_conditions_[node_id] = std::move(bc);
}

void PoissonSolver::clearBoundaryConditions()
{
    boundary_conditions_.clear();
}

void PoissonSolver::addMaterialProperty(Mesh::MaterialId material_id, MaterialPropertyPtr property)
{
    material_properties_[material_id] = std::move(property);
}

void PoissonSolver::clearMaterialProperties()
{
    material_properties_.clear();
}

void PoissonSolver::setChargeDensity(Mesh::NodeId node_id, double density)
{
    charge_densities_[node_id] = density;
}

void PoissonSolver::setChargeDensity(const std::unordered_map<Mesh::NodeId, double>& densities)
{
    charge_densities_ = densities;
}

void PoissonSolver::clearChargeDensities()
{
    charge_densities_.clear();
}

void PoissonSolver::assembleSystem()
{
    if (!hasValidMesh())
    {
        return;
    }

    buildNodeIndexMapping();
    if (!system_ || system_->getSize() != node_to_index_.size())
    {
        system_ = std::make_shared<LinearSystem>(node_to_index_.size());
    }

    system_->zero();
    for (const auto& element : getMeshElements())
    {
        assembleElementMatrix(element);
        assembleElementRHS(element);
    }
    applyBoundaryConditions();
    system_->markAssembled();
}

void PoissonSolver::assembleElementMatrix(Mesh::ElementPtr element)
{
    if (!element || !system_)
    {
        return;
    }

    const auto& nodes = element->getNodes();
    if (nodes.size() < 2)
    {
        return;
    }

    const double permittivity = resolvePermittivity(getPermittivity(element));
    const double measure = elementMeasure(element);
    const double edge_scale =
        measure / static_cast<double>((nodes.size() * (nodes.size() - 1)) / 2);

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        if (!nodes[i])
        {
            continue;
        }

        const auto index_i = node_to_index_.find(nodes[i]->getId());
        if (index_i == node_to_index_.end())
        {
            continue;
        }

        for (std::size_t j = i + 1; j < nodes.size(); ++j)
        {
            if (!nodes[j])
            {
                continue;
            }

            const auto index_j = node_to_index_.find(nodes[j]->getId());
            if (index_j == node_to_index_.end())
            {
                continue;
            }

            const double distance =
                std::max(nodes[i]->getPosition().distanceTo(nodes[j]->getPosition()),
                         kGeometryTolerance);
            const double weight = permittivity * edge_scale / (distance * distance);

            system_->addToMatrix(index_i->second, index_i->second, weight);
            system_->addToMatrix(index_j->second, index_j->second, weight);
            system_->addToMatrix(index_i->second, index_j->second, -weight);
            system_->addToMatrix(index_j->second, index_i->second, -weight);
        }
    }
}

void PoissonSolver::assembleElementRHS(Mesh::ElementPtr element)
{
    if (!element || !system_)
    {
        return;
    }

    const auto& nodes = element->getNodes();
    if (nodes.empty())
    {
        return;
    }

    const double nodal_share = elementMeasure(element) / static_cast<double>(nodes.size());
    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        const auto index = node_to_index_.find(node->getId());
        const auto density = charge_densities_.find(node->getId());
        if (index != node_to_index_.end() && density != charge_densities_.end())
        {
            system_->addToRHS(index->second, -density->second * nodal_share);
        }
    }
}

void PoissonSolver::applyBoundaryConditions()
{
    if (!system_)
    {
        return;
    }

    for (const auto& [node_id, boundary] : boundary_conditions_)
    {
        const auto index = node_to_index_.find(node_id);
        if (index == node_to_index_.end() || !boundary)
        {
            continue;
        }

        if (boundary->getType() == BoundaryConditionType::DIRICHLET)
        {
            system_->applyDirichletBC(index->second, boundary->getValue());
        }
        else if (boundary->getType() == BoundaryConditionType::NEUMANN)
        {
            system_->applyNeumannBC(index->second, boundary->getValue());
        }
    }
}

bool PoissonSolver::solve(SolverType solver_type)
{
    if (!system_ || !system_->isAssembled())
    {
        assembleSystem();
    }

    const bool solved = system_->solve(solver_type);
    if (solved)
    {
        updateNodePotentials();
        updateNodeElectricFields();
    }
    return solved;
}

double PoissonSolver::getPotential(Mesh::NodeId node_id) const
{
    const auto index = node_to_index_.find(node_id);
    if (!system_ || index == node_to_index_.end())
    {
        return 0.0;
    }
    return system_->getSolution(index->second);
}

Utils::Vector3D PoissonSolver::getElectricField(Mesh::NodeId node_id) const
{
    const auto nodes = getMeshNodes();
    if (const auto node = findNodeById(nodes, node_id))
    {
        return node->getElectricField();
    }

    const auto index = node_to_index_.find(node_id);
    if (!system_ || index == node_to_index_.end())
    {
        return Utils::Vector3D(0, 0, 0);
    }
    return computeSimpleFiniteDifferenceGradient(index->second);
}

Utils::Vector3D PoissonSolver::getElectricField(const Utils::Point3D& position) const
{
    initializeSpatialIndex();
    if (spatial_index_)
    {
        const auto candidate_elements = spatial_index_->findElementsContaining(position);
        for (const auto& element : candidate_elements)
        {
            if (element && element->contains(position))
            {
                return interpolateElectricFieldInElement(element, position);
            }
        }
    }

    return getElectricFieldNearestNeighbor(position);
}

Utils::Vector3D PoissonSolver::getElectricFieldAtPosition(const Utils::Point3D& position) const
{
    return getElectricField(position);
}

void PoissonSolver::updateNodePotentials()
{
    if (!system_)
    {
        return;
    }

    const auto nodes = getMeshNodes();
    for (const auto& [index, node_id] : index_to_node_)
    {
        if (const auto node = findNodeById(nodes, node_id))
        {
            node->setPotential(system_->getSolution(index));
        }
    }
}

void PoissonSolver::updateNodeElectricFields()
{
    computeElectricFieldGradient();
}

void PoissonSolver::computeElectricFieldGradient()
{
    const auto nodes = getMeshNodes();
    if (nodes.empty())
    {
        return;
    }

    const auto neighbors = buildNodeNeighbors(getMeshElements());
    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }

        std::vector<Utils::Vector3D> delta_pos;
        std::vector<double> delta_phi;

        if (const auto it = neighbors.find(node->getId()); it != neighbors.end())
        {
            delta_pos.reserve(it->second.size());
            delta_phi.reserve(it->second.size());

            for (Mesh::NodeId neighbor_id : it->second)
            {
                const auto neighbor = findNodeById(nodes, neighbor_id);
                if (!neighbor)
                {
                    continue;
                }

                const auto displacement = neighbor->getPosition() - node->getPosition();
                if (displacement.magnitudeSquared() <= kGeometryTolerance * kGeometryTolerance)
                {
                    continue;
                }

                delta_pos.push_back(displacement);
                delta_phi.push_back(neighbor->getPotential() - node->getPotential());
            }
        }

        Utils::Vector3D gradient(0.0, 0.0, 0.0);
        if (!delta_pos.empty())
        {
            gradient = computeLeastSquaresGradient(delta_pos, delta_phi);
        }
        else if (const auto index = node_to_index_.find(node->getId()); index != node_to_index_.end())
        {
            gradient = computeSimpleGradientFallback(index->second);
        }

        node->setElectricField(Utils::Vector3D(-gradient.x(), -gradient.y(), -gradient.z()));
    }
}

double PoissonSolver::getTotalEnergy() const
{
    const auto nodes = getMeshNodes();
    if (nodes.empty())
    {
        return 0.0;
    }

    const auto control_volumes = buildNodalControlVolumes(getMeshElements());
    double total_energy = 0.0;
    for (const auto& node : nodes)
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

        total_energy += 0.5 * kVacuumPermittivity * node->getElectricField().magnitudeSquared() *
                        control_volume;
    }

    return total_energy;
}

double PoissonSolver::getMaxPotential() const
{
    if (!system_ || system_->getSize() == 0)
    {
        return 0.0;
    }

    double maximum = -std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < system_->getSize(); ++i)
    {
        maximum = std::max(maximum, system_->getSolution(i));
    }
    return maximum;
}

double PoissonSolver::getMinPotential() const
{
    if (!system_ || system_->getSize() == 0)
    {
        return 0.0;
    }

    double minimum = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < system_->getSize(); ++i)
    {
        minimum = std::min(minimum, system_->getSolution(i));
    }
    return minimum;
}

Utils::Vector3D PoissonSolver::getMaxElectricField() const
{
    Utils::Vector3D best(0, 0, 0);
    double best_magnitude = 0.0;
    for (const auto& [node_id, index] : node_to_index_)
    {
        static_cast<void>(index);
        const Utils::Vector3D electric_field = getElectricField(node_id);
        const double magnitude = electric_field.magnitude();
        if (magnitude > best_magnitude)
        {
            best_magnitude = magnitude;
            best = electric_field;
        }
    }
    return best;
}

void PoissonSolver::printSolution() const
{
    if (system_)
    {
        system_->printSolution();
    }
}

std::string PoissonSolver::toString() const
{
    std::ostringstream stream;
    stream << "PoissonSolver(nodes=" << node_to_index_.size() << ", bcs="
           << boundary_conditions_.size() << ", materials=" << material_properties_.size()
           << ", charges=" << charge_densities_.size() << ")";
    return stream.str();
}

void PoissonSolver::buildNodeIndexMapping()
{
    node_to_index_.clear();
    index_to_node_.clear();

    const auto nodes = getMeshNodes();
    std::size_t next_index = 0;
    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }
        node_to_index_[node->getId()] = next_index;
        index_to_node_[next_index] = node->getId();
        ++next_index;
    }
}

MaterialPropertyPtr PoissonSolver::getMaterialProperty(Mesh::MaterialId material_id) const
{
    const auto it = material_properties_.find(material_id);
    if (it != material_properties_.end())
    {
        return it->second;
    }
    return std::make_shared<MaterialProperty>(material_id, Mesh::MaterialType::VACUUM, "default");
}

double PoissonSolver::getPermittivity(Mesh::ElementPtr element) const
{
    if (!element)
    {
        return 1.0;
    }
    return getMaterialProperty(element->getMaterialId())->getPermittivity();
}

void PoissonSolver::assembleTetrahedronElement(Mesh::TetrahedronPtr tetrahedron)
{
    assembleElementMatrix(tetrahedron);
    assembleElementRHS(tetrahedron);
}

void PoissonSolver::assembleTriangleElement(Mesh::TrianglePtr triangle)
{
    assembleElementMatrix(triangle);
    assembleElementRHS(triangle);
}

void PoissonSolver::initializeSpatialIndex() const
{
    if (spatial_index_ || !mesh_)
    {
        return;
    }

    spatial_index_ = std::make_shared<Mesh::SpatialIndex>();
    spatial_index_->initialize(mesh_->getElements());
}

Utils::Vector3D PoissonSolver::interpolateElectricFieldInElement(Mesh::ElementPtr element,
                                                                 const Utils::Point3D&) const
{
    if (!element)
    {
        return Utils::Vector3D(0, 0, 0);
    }

    const auto& nodes = element->getNodes();
    if (nodes.empty())
    {
        return Utils::Vector3D(0, 0, 0);
    }

    Utils::Vector3D accumulated(0.0, 0.0, 0.0);
    double count = 0.0;
    for (const auto& node : nodes)
    {
        if (!node)
        {
            continue;
        }
        accumulated = accumulated + node->getElectricField();
        count += 1.0;
    }

    if (count <= 0.0)
    {
        return Utils::Vector3D(0, 0, 0);
    }

    return accumulated / count;
}

Utils::Vector3D PoissonSolver::getElectricFieldNearestNeighbor(const Utils::Point3D& position) const
{
    const auto nodes = getMeshNodes();
    if (nodes.empty())
    {
        return Utils::Vector3D(0, 0, 0);
    }

    Mesh::NodeId nearest = nodes.front()->getId();
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
            nearest = node->getId();
        }
    }

    return getElectricField(nearest);
}

Utils::Vector3D PoissonSolver::computeLeastSquaresGradient(
    const std::vector<Utils::Vector3D>& delta_pos, const std::vector<double>& delta_phi) const
{
    if (delta_pos.empty() || delta_pos.size() != delta_phi.size())
    {
        return Utils::Vector3D(0, 0, 0);
    }

    std::array<std::array<double, 3>, 3> normal_matrix{};
    std::array<double, 3> rhs{};
    for (std::size_t i = 0; i < delta_pos.size(); ++i)
    {
        const double dx = delta_pos[i].x();
        const double dy = delta_pos[i].y();
        const double dz = delta_pos[i].z();

        normal_matrix[0][0] += dx * dx;
        normal_matrix[0][1] += dx * dy;
        normal_matrix[0][2] += dx * dz;
        normal_matrix[1][0] += dy * dx;
        normal_matrix[1][1] += dy * dy;
        normal_matrix[1][2] += dy * dz;
        normal_matrix[2][0] += dz * dx;
        normal_matrix[2][1] += dz * dy;
        normal_matrix[2][2] += dz * dz;

        rhs[0] += dx * delta_phi[i];
        rhs[1] += dy * delta_phi[i];
        rhs[2] += dz * delta_phi[i];
    }

    std::array<double, 3> solution{};
    if (!solve3x3(normal_matrix, rhs, solution))
    {
        return Utils::Vector3D(0, 0, 0);
    }

    return Utils::Vector3D(solution[0], solution[1], solution[2]);
}

Utils::Vector3D PoissonSolver::computeSimpleFiniteDifferenceGradient(std::size_t index) const
{
    if (!system_ || system_->getSize() == 0 || index >= system_->getSize())
    {
        return Utils::Vector3D(0, 0, 0);
    }

    const auto node_it = index_to_node_.find(index);
    if (node_it == index_to_node_.end())
    {
        return Utils::Vector3D(0, 0, 0);
    }

    const auto nodes = getMeshNodes();
    const auto neighbors = buildNodeNeighbors(getMeshElements());
    const auto center = findNodeById(nodes, node_it->second);
    if (!center)
    {
        return Utils::Vector3D(0, 0, 0);
    }

    const auto neighbor_it = neighbors.find(center->getId());
    if (neighbor_it == neighbors.end() || neighbor_it->second.empty())
    {
        return Utils::Vector3D(0, 0, 0);
    }

    std::vector<Utils::Vector3D> delta_pos;
    std::vector<double> delta_phi;
    delta_pos.reserve(neighbor_it->second.size());
    delta_phi.reserve(neighbor_it->second.size());

    for (Mesh::NodeId neighbor_id : neighbor_it->second)
    {
        const auto neighbor = findNodeById(nodes, neighbor_id);
        if (!neighbor)
        {
            continue;
        }

        const auto displacement = neighbor->getPosition() - center->getPosition();
        if (displacement.magnitudeSquared() <= kGeometryTolerance * kGeometryTolerance)
        {
            continue;
        }

        delta_pos.push_back(displacement);
        delta_phi.push_back(neighbor->getPotential() - center->getPotential());
    }

    return computeLeastSquaresGradient(delta_pos, delta_phi);
}

Utils::Vector3D PoissonSolver::computeSimpleGradientFallback(std::size_t index) const
{
    return computeSimpleFiniteDifferenceGradient(index);
}

bool PoissonSolver::hasValidMesh() const
{
    return mesh_ && !mesh_->getNodes().empty() && !mesh_->getElements().empty();
}

std::vector<Mesh::ElementPtr> PoissonSolver::getMeshElements() const
{
    if (!mesh_)
    {
        return {};
    }
    return mesh_->getElements();
}

std::vector<Mesh::NodePtr> PoissonSolver::getMeshNodes() const
{
    if (!mesh_)
    {
        return {};
    }
    return mesh_->getNodes();
}

} // namespace FieldSolver
} // namespace SCDAT
