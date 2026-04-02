#pragma once

#include "../../Boundary/include/BoundaryType.h"
#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Mesh/include/MeshAlgorithms.h"
#include "../../Mesh/include/MeshParsing.h"
#include "../../Mesh/include/MeshPartitioning.h"

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#else
namespace Eigen
{
using VectorXd = std::vector<double>;

template <typename T> class SparseMatrix
{
  public:
    using Scalar = T;

    SparseMatrix() = default;
    SparseMatrix(int rows, int cols) : rows_(rows), cols_(cols) {}

    void resize(int rows, int cols)
    {
        rows_ = rows;
        cols_ = cols;
    }

    int rows() const
    {
        return rows_;
    }

    int cols() const
    {
        return cols_;
    }

  private:
    int rows_ = 0;
    int cols_ = 0;
};

using Matrix3d = std::array<std::array<double, 3>, 3>;
} // namespace Eigen
#endif

namespace SCDAT
{
namespace FieldSolver
{

enum class SolverType
{
    CONJUGATE_GRADIENT,
    BICONJUGATE_GRADIENT,
    GMRES,
    DIRECT_LU,
    MULTIGRID
};

enum class PreconditionerType
{
    NONE,
    JACOBI,
    ILU,
    SSOR,
    MULTIGRID
};

enum class MultigridCycleType
{
    V_CYCLE,
    W_CYCLE,
    F_CYCLE
};

struct BoundaryCondition
{
    std::string name;
    BoundaryConditionType type = BoundaryConditionType::DIRICHLET;
    double value = 0.0;
    std::vector<Mesh::NodeId> nodes;
    std::vector<Mesh::NodeId> paired_nodes;
    int priority = 0;
    double alpha = 1.0;
    double beta = 0.0;
};

struct SolverConfiguration
{
    double tolerance = 1e-12;
    int max_iterations = 10000;
    SolverType solver_type = SolverType::CONJUGATE_GRADIENT;
    bool use_preconditioning = true;
    PreconditionerType preconditioner_type = PreconditionerType::ILU;
};

struct MultigridConfiguration
{
    int num_levels = 4;
    int coarsening_factor = 2;
    int pre_smoothing_steps = 2;
    int post_smoothing_steps = 2;
    MultigridCycleType cycle_type = MultigridCycleType::V_CYCLE;
};

struct CouplingIterationConfiguration
{
    int max_outer_iterations = 1;
    double potential_tolerance_v = 1.0e-3;
    double relaxation = 1.0;
};

struct ElementGeometry
{
    double volume = 0.0;
    Utils::Point3D centroid{};
    std::vector<Utils::Vector3D> face_normals;
    std::vector<double> face_areas;
    std::array<Utils::Point3D, 4> tetra_nodes{};
};

class MultiBoundaryElectricFieldSolver
{
  public:
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using DenseVector = Eigen::VectorXd;
    using NodeId = Mesh::NodeId;
        using CouplingCallback = std::function<void(int, const std::vector<double>&)>;

    explicit MultiBoundaryElectricFieldSolver(Mesh::MeshPtr mesh);
    virtual ~MultiBoundaryElectricFieldSolver() = default;

    void addBoundaryCondition(const std::string& boundary_name, BoundaryConditionType type,
                              double value, const std::vector<NodeId>& nodes);
    void addPeriodicBoundaryCondition(const std::string& boundary_name,
                                      const std::vector<NodeId>& primary_nodes,
                                      const std::vector<NodeId>& paired_nodes,
                                      double potential_offset = 0.0);
    void setChargeDensity(NodeId node_id, double charge_density);
    void setMaterialPermittivity(Mesh::MaterialId material_id, double permittivity);
    void setSolverConfiguration(const SolverConfiguration& config);
    void setMultigridConfiguration(const MultigridConfiguration& config);
    void setCouplingIterationConfiguration(const CouplingIterationConfiguration& config);

    bool solve();
    bool solveCoupled(const CouplingCallback& callback = CouplingCallback{});

    double getPotential(NodeId node_id) const;
    Utils::Vector3D getElectricField(NodeId node_id) const;
    Utils::Vector3D getElectricField(const Utils::Point3D& position) const;

    const std::vector<double>& getResidualHistory() const
    {
        return residual_history_;
    }

    int getIterationCount() const
    {
        return current_iteration_;
    }

    const std::vector<double>& getCouplingResidualHistory() const
    {
        return coupling_residual_history_;
    }

    int getCouplingIterationCount() const
    {
        return coupling_iteration_count_;
    }

    double validateSolution() const;
    double calculateElectrostaticEnergy() const;

  private:
    Mesh::MeshPtr mesh_;
    std::unordered_map<std::string, BoundaryCondition> boundary_conditions_;
    std::unordered_map<NodeId, double> charge_densities_;
    std::unordered_map<Mesh::MaterialId, double> material_permittivities_;
    std::unordered_map<NodeId, std::string> node_boundary_map_;
    std::unordered_map<NodeId, size_t> node_index_map_;
    std::vector<NodeId> index_to_node_;

    SparseMatrix system_matrix_;
    DenseVector rhs_vector_;
    DenseVector solution_vector_;
    std::vector<std::vector<double>> dense_system_matrix_;
    std::vector<double> rhs_storage_;
    std::vector<double> solution_storage_;
    std::vector<double> preconditioner_diagonal_;
    std::unordered_map<NodeId, Utils::Vector3D> electric_field_cache_;
    mutable Mesh::SpatialIndexPtr spatial_index_;

    SolverConfiguration solver_config_;
    MultigridConfiguration multigrid_config_;
    CouplingIterationConfiguration coupling_config_;

    double convergence_tolerance_ = 1e-12;
    int max_iterations_ = 10000;
    int current_iteration_ = 0;
    std::vector<double> residual_history_;
    int coupling_iteration_count_ = 0;
    std::vector<double> coupling_residual_history_;

    bool use_adaptive_refinement_ = false;
    bool use_multigrid_ = false;

    void initializeSolverConfig();
    bool assembleSystem();
    void assembleElement(Mesh::ElementPtr element);
    void assembleTetrahedronElement(Mesh::ElementPtr element);
    void assembleHexahedronElement(Mesh::ElementPtr element);
    void applyBoundaryConditions();
    void applyBoundaryCondition(const BoundaryCondition& bc);
    bool solveWithMultigrid();
    bool solveWithDirectMethod();
    bool conjugateGradientSolve();
    std::vector<double> snapshotNodePotentials() const;
    void applyRelaxedNodePotentials(const std::vector<double>& previous,
                                    const std::vector<double>& current,
                                    double relaxation);
    void performVCycle();
    void calculateElectricField();
    double calculateResidualNorm() const;
    void setupPreconditioner();

    ElementGeometry calculateTetrahedronGeometry(Mesh::TetrahedronPtr tetrahedron) const;
    std::vector<Utils::Vector3D>
    calculateShapeFunctionGradients(const ElementGeometry& geometry) const;
    std::pair<Eigen::Matrix3d, std::vector<Utils::Vector3D>>
    calculateHexahedronShapeGradients(Mesh::TetrahedronPtr hexahedron, double xi, double eta,
                                      double zeta) const;

    size_t getGlobalIndex(NodeId node_id) const;
    double getMaterialPermittivity(Mesh::ElementPtr element) const;
};

} // namespace FieldSolver
} // namespace SCDAT
