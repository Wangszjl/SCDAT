#ifndef SCDAT_FIELD_POISSON_SOLVER_H
#define SCDAT_FIELD_POISSON_SOLVER_H

#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Mesh/include/MeshAlgorithms.h"
#include "../../Mesh/include/MeshPartitioning.h"
#include "../../Mesh/include/MeshParsing.h"
#include "../../Solver/include/LinearSystem.h"
#include "BoundaryCondition.h"
#include "MaterialProperty.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

using ::SCDAT::Solver::LinearSystem;
using ::SCDAT::Solver::SolverType;
using LinearSystemPtr = std::shared_ptr<LinearSystem>;

class PoissonSolver
{
  public:
    explicit PoissonSolver(Mesh::VolMeshPtr mesh) : mesh_(std::move(mesh)), system_(nullptr) {}
    virtual ~PoissonSolver() = default;

    void initialize();
    void setMesh(Mesh::VolMeshPtr mesh) { mesh_ = std::move(mesh); }

    void addBoundaryCondition(Mesh::NodeId node_id, BoundaryConditionPtr bc);
    void clearBoundaryConditions();

    void addMaterialProperty(Mesh::MaterialId material_id, MaterialPropertyPtr property);
    void clearMaterialProperties();

    void setChargeDensity(Mesh::NodeId node_id, double density);
    void setChargeDensity(const std::unordered_map<Mesh::NodeId, double>& densities);
    void clearChargeDensities();

    void assembleSystem();
    void assembleElementMatrix(Mesh::ElementPtr element);
    void assembleElementRHS(Mesh::ElementPtr element);
    void applyBoundaryConditions();

    bool solve(SolverType solver_type = SolverType::DIRECT);

    double getPotential(Mesh::NodeId node_id) const;
    Utils::Vector3D getElectricField(Mesh::NodeId node_id) const;
    Utils::Vector3D getElectricField(const Utils::Point3D& position) const;
    Utils::Vector3D getElectricFieldAtPosition(const Utils::Point3D& position) const;

    void updateNodePotentials();
    void updateNodeElectricFields();
    void computeElectricFieldGradient();

    double getTotalEnergy() const;
    double getMaxPotential() const;
    double getMinPotential() const;
    Utils::Vector3D getMaxElectricField() const;

    void printSolution() const;
    std::string toString() const;

  private:
    Mesh::VolMeshPtr mesh_;
    LinearSystemPtr system_;
    std::unordered_map<Mesh::NodeId, BoundaryConditionPtr> boundary_conditions_;
    std::unordered_map<Mesh::MaterialId, MaterialPropertyPtr> material_properties_;
    std::unordered_map<Mesh::NodeId, double> charge_densities_;
    std::unordered_map<Mesh::NodeId, std::size_t> node_to_index_;
    std::unordered_map<std::size_t, Mesh::NodeId> index_to_node_;
    mutable Mesh::SpatialIndexPtr spatial_index_;

    void buildNodeIndexMapping();
    MaterialPropertyPtr getMaterialProperty(Mesh::MaterialId material_id) const;
    double getPermittivity(Mesh::ElementPtr element) const;

    void assembleTetrahedronElement(Mesh::TetrahedronPtr tetrahedron);
    void assembleTriangleElement(Mesh::TrianglePtr triangle);

    void initializeSpatialIndex() const;
    Utils::Vector3D interpolateElectricFieldInElement(Mesh::ElementPtr element,
                                                      const Utils::Point3D& position) const;
    Utils::Vector3D getElectricFieldNearestNeighbor(const Utils::Point3D& position) const;
    Utils::Vector3D computeLeastSquaresGradient(const std::vector<Utils::Vector3D>& delta_pos,
                                                const std::vector<double>& delta_phi) const;

    Utils::Vector3D computeSimpleFiniteDifferenceGradient(std::size_t index) const;
    Utils::Vector3D computeSimpleGradientFallback(std::size_t index) const;

    bool hasValidMesh() const;
    std::vector<Mesh::ElementPtr> getMeshElements() const;
    std::vector<Mesh::NodePtr> getMeshNodes() const;
};

using PoissonSolverPtr = std::shared_ptr<PoissonSolver>;

} // namespace FieldSolver
} // namespace SCDAT

#endif // SCDAT_FIELD_POISSON_SOLVER_H
