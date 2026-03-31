#pragma once

#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Mesh/include/ModernMesh.h"
#include <memory>
#include <vector>

namespace SCDAT
{
namespace Solver
{

class HighOrderFieldOperators
{
  public:
    enum class MethodType
    {
        FINITE_DIFFERENCE_4TH,
        FINITE_DIFFERENCE_6TH,
        FINITE_ELEMENT_P2,
        FINITE_ELEMENT_P3,
        SPECTRAL_METHOD,
        DISCONTINUOUS_GALERKIN
    };

    struct Configuration
    {
        MethodType method_type = MethodType::FINITE_DIFFERENCE_4TH;
        int polynomial_order = 4;
        double accuracy_tolerance = 1e-8;
        bool use_adaptive_refinement = true;
        bool use_error_control = true;
        int max_refinement_levels = 5;
        double refinement_threshold = 1e-6;
    };

    struct Statistics
    {
        double numerical_error = 0.0;
        double convergence_rate = 0.0;
        int refinement_iterations = 0;
        double solve_time_ms = 0.0;
        std::size_t total_dofs = 0;
        bool converged = false;
    };

    HighOrderFieldOperators();
    explicit HighOrderFieldOperators(const Configuration& config);
    ~HighOrderFieldOperators();

    void setConfiguration(const Configuration& config) { config_ = config; }
    const Configuration& getConfiguration() const { return config_; }

    bool solveFieldHighOrder(const Mesh::ModernMesh& mesh, const std::vector<double>& charge_density,
                             std::vector<double>& potential,
                             std::vector<Utils::Vector3D>& electric_field);

    bool performAdaptiveMeshRefinement(Mesh::ModernMesh& mesh,
                                       const std::vector<double>& error_indicators);

    std::vector<double> estimateElementErrors(const Mesh::ModernMesh& mesh,
                                              const std::vector<double>& solution);

    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics() { stats_ = Statistics{}; }

  private:
    void computeElectricFieldFromPotential(const Mesh::ModernMesh& mesh,
                                           const std::vector<double>& potential,
                                           std::vector<Utils::Vector3D>& electric_field) const;

    bool solveField4thOrderFD(const Mesh::ModernMesh& mesh,
                              const std::vector<double>& charge_density,
                              std::vector<double>& potential);

    bool solveFieldHighOrderFE(const Mesh::ModernMesh& mesh,
                               const std::vector<double>& charge_density,
                               std::vector<double>& potential);

    std::vector<double> computeHighOrderBasisFunctions(const Utils::Point3D& point,
                                                       const Mesh::ElementPtr& element);

    double computeErrorIndicator(const Mesh::ElementPtr& element,
                                 const std::vector<double>& solution);

    bool performHRefinement(Mesh::ModernMesh& mesh, const std::vector<bool>& refine_flags);
    bool performPRefinement(Mesh::ModernMesh& mesh,
                            const std::vector<int>& polynomial_orders);
    double computeConvergenceRate(const std::vector<double>& errors);

  private:
    Configuration config_;
    Statistics stats_;
    std::vector<double> previous_errors_;
    std::vector<int> refinement_levels_;
    mutable std::vector<double> temp_vector_;
    mutable std::vector<std::vector<double>> basis_cache_;
};

using HighOrderFieldOperatorsPtr = std::shared_ptr<HighOrderFieldOperators>;

} // namespace Solver
} // namespace SCDAT
