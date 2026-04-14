#ifndef SCDAT_SOLVER_SURFACE_SOLVER_FACADE_H
#define SCDAT_SOLVER_SURFACE_SOLVER_FACADE_H

#include <cstddef>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Solver
{

struct SolverPolicyFlags
{
    std::string normalized_coupling_mode;
    std::string normalized_convergence_policy;
    bool implicit_coupling_requested = false;
    bool residual_guarded_requested = false;
    bool fixed_iteration_policy_requested = false;
};

std::string normalizeSolverPolicyToken(std::string text);
bool isImplicitSolverCouplingModeToken(const std::string& text);
bool isResidualGuardedSolverPolicyToken(const std::string& text);
bool isFixedIterationSolverPolicyToken(const std::string& text);
SolverPolicyFlags resolveSolverPolicyFlags(const std::string& coupling_mode,
                                           const std::string& convergence_policy);

struct GlobalCoupledControlInput
{
    std::size_t configured_iteration_limit = 0;
    double configured_tolerance_v = 1.0e-6;
    double configured_relaxation = 1.0;
    std::size_t solver_max_iterations = 0;
    double solver_residual_tolerance = 0.0;
    double solver_relaxation_factor = 1.0;
    std::string coupling_mode;
    std::string convergence_policy;
};

struct GlobalCoupledControl
{
    std::size_t iteration_limit = 0;
    double tolerance_v = 1.0e-6;
    double relaxation = 1.0;
    bool implicit_coupling_requested = false;
    bool residual_guarded_requested = false;
};

GlobalCoupledControl resolveGlobalCoupledControl(const GlobalCoupledControlInput& input);

struct DenseLinearSystemSolveResult
{
    std::vector<double> solution;
    double residual_norm = 0.0;
    bool solved = false;
    std::size_t nonzeros = 0;
};

DenseLinearSystemSolveResult solveDenseLinearSystemWithResidual(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs);

struct IterativeLinearSystemSolveResult
{
    std::vector<double> solution;
    int iterations = 0;
    bool converged = false;
    double residual_norm = 0.0;
    std::size_t nonzeros = 0;
};

IterativeLinearSystemSolveResult solveIterativeLinearSystem(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& rhs,
    const std::vector<double>& initial_guess,
    int max_iterations,
    double tolerance,
    double relaxation);

struct VolumeLinearSolverRoutingInput
{
    std::size_t system_size = 0;
    bool has_external_cells = false;
    bool has_external_cell_faces = false;
    bool iterative_only = false;
    bool dense_only = false;
};

struct VolumeLinearSolverRouting
{
    bool auto_prefers_iterative = false;
    bool use_iterative_solver = false;
    std::string solver_mode = "dense";
};

VolumeLinearSolverRouting resolveVolumeLinearSolverRouting(
    const VolumeLinearSolverRoutingInput& input);

} // namespace Solver
} // namespace SCDAT

#endif // SCDAT_SOLVER_SURFACE_SOLVER_FACADE_H
