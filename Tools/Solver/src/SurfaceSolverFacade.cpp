#include "../include/SurfaceSolverFacade.h"
#include "../../Basic/include/StringTokenUtils.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace SCDAT
{
namespace Solver
{
namespace
{
struct SparseLinearEntry
{
    std::size_t column = 0;
    double value = 0.0;
};

struct SparseLinearRow
{
    double diagonal = 0.0;
    std::vector<SparseLinearEntry> off_diagonal;
};

std::vector<SparseLinearRow> buildSparseLinearRows(const std::vector<std::vector<double>>& matrix,
                                                   std::size_t* nonzero_count)
{
    std::vector<SparseLinearRow> rows(matrix.size());
    std::size_t nnz = 0;
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        rows[i].diagonal = i < matrix[i].size() ? matrix[i][i] : 0.0;
        for (std::size_t j = 0; j < matrix[i].size(); ++j)
        {
            if (j == i)
            {
                if (std::abs(matrix[i][j]) > 1.0e-20)
                {
                    ++nnz;
                }
                continue;
            }
            if (std::abs(matrix[i][j]) <= 1.0e-20)
            {
                continue;
            }
            rows[i].off_diagonal.push_back({j, matrix[i][j]});
            ++nnz;
        }
    }
    if (nonzero_count != nullptr)
    {
        *nonzero_count = nnz;
    }
    return rows;
}

double sparseLinearResidualNorm(const std::vector<SparseLinearRow>& rows,
                                const std::vector<double>& rhs,
                                const std::vector<double>& solution)
{
    if (rows.size() != rhs.size() || rhs.size() != solution.size())
    {
        return 0.0;
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i)
    {
        double residual = rhs[i] - rows[i].diagonal * solution[i];
        for (const auto& entry : rows[i].off_diagonal)
        {
            residual -= entry.value * solution[entry.column];
        }
        sum += residual * residual;
    }
    return std::sqrt(sum / std::max<std::size_t>(1, rhs.size()));
}
} // namespace

std::string normalizeSolverPolicyToken(std::string text)
{
    return Basic::normalizeAlnumToken(text);
}

bool isImplicitSolverCouplingModeToken(const std::string& text)
{
    const std::string token = normalizeSolverPolicyToken(text);
    return token == "implicit" || token == "fieldparticleimplicit" ||
           token == "fieldparticlewallimplicit" || token == "globalcoupled" ||
           token == "fullyimplicit" || token == "iterativecoupled";
}

bool isResidualGuardedSolverPolicyToken(const std::string& text)
{
    const std::string token = normalizeSolverPolicyToken(text);
    return token == "residualguarded" || token == "residualnormguarded" ||
           token == "adaptive";
}

bool isFixedIterationSolverPolicyToken(const std::string& text)
{
    const std::string token = normalizeSolverPolicyToken(text);
    return token == "fixediteration" || token == "fixed";
}

SolverPolicyFlags resolveSolverPolicyFlags(const std::string& coupling_mode,
                                           const std::string& convergence_policy)
{
    SolverPolicyFlags flags;
    flags.normalized_coupling_mode = normalizeSolverPolicyToken(coupling_mode);
    flags.normalized_convergence_policy = normalizeSolverPolicyToken(convergence_policy);
    flags.implicit_coupling_requested =
        isImplicitSolverCouplingModeToken(flags.normalized_coupling_mode);
    flags.residual_guarded_requested =
        isResidualGuardedSolverPolicyToken(flags.normalized_convergence_policy);
    flags.fixed_iteration_policy_requested =
        isFixedIterationSolverPolicyToken(flags.normalized_convergence_policy);
    return flags;
}

GlobalCoupledControl resolveGlobalCoupledControl(const GlobalCoupledControlInput& input)
{
    GlobalCoupledControl control;
    control.iteration_limit = input.configured_iteration_limit;
    control.tolerance_v = std::max(1.0e-9, input.configured_tolerance_v);
    control.relaxation = std::clamp(input.configured_relaxation, 0.05, 1.0);

    const auto policy_flags =
        resolveSolverPolicyFlags(input.coupling_mode, input.convergence_policy);
    control.implicit_coupling_requested = policy_flags.implicit_coupling_requested;
    control.residual_guarded_requested = policy_flags.residual_guarded_requested;

    if ((control.implicit_coupling_requested || control.residual_guarded_requested) &&
        input.solver_max_iterations > 0)
    {
        const std::size_t solver_limit =
            std::clamp<std::size_t>(input.solver_max_iterations, 2, 256);
        control.iteration_limit = std::max(control.iteration_limit, solver_limit);
    }

    if (control.implicit_coupling_requested && control.iteration_limit < 2)
    {
        control.iteration_limit = 2;
    }

    if ((control.implicit_coupling_requested || control.residual_guarded_requested) &&
        std::isfinite(input.solver_residual_tolerance) && input.solver_residual_tolerance > 0.0)
    {
        control.tolerance_v =
            std::min(control.tolerance_v, std::max(1.0e-9, input.solver_residual_tolerance));
    }

    if ((control.implicit_coupling_requested || control.residual_guarded_requested) &&
        std::isfinite(input.solver_relaxation_factor))
    {
        control.relaxation = std::clamp(input.solver_relaxation_factor, 0.05, 1.0);
    }

    return control;
}

DenseLinearSystemSolveResult solveDenseLinearSystemWithResidual(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs)
{
    DenseLinearSystemSolveResult result;
    const std::size_t n = rhs.size();
    result.solution.assign(n, 0.0);
    if (matrix.size() != n)
    {
        return result;
    }

    for (const auto& row : matrix)
    {
        if (row.size() != n)
        {
            return result;
        }

        for (const double value : row)
        {
            if (std::abs(value) > 1.0e-20)
            {
                ++result.nonzeros;
            }
        }
    }

    const auto matrix_original = matrix;
    const auto rhs_original = rhs;
    for (std::size_t pivot = 0; pivot < n; ++pivot)
    {
        std::size_t best_row = pivot;
        double best_value = std::abs(matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < n; ++row)
        {
            const double candidate = std::abs(matrix[row][pivot]);
            if (candidate > best_value)
            {
                best_value = candidate;
                best_row = row;
            }
        }

        if (best_value <= 1.0e-20)
        {
            return result;
        }

        if (best_row != pivot)
        {
            std::swap(matrix[pivot], matrix[best_row]);
            std::swap(rhs[pivot], rhs[best_row]);
        }

        const double diagonal = matrix[pivot][pivot];
        for (std::size_t row = pivot + 1; row < n; ++row)
        {
            const double factor = matrix[row][pivot] / diagonal;
            if (!std::isfinite(factor) || std::abs(factor) <= 0.0)
            {
                continue;
            }

            rhs[row] -= factor * rhs[pivot];
            for (std::size_t column = pivot; column < n; ++column)
            {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
        }
    }

    for (std::size_t back = n; back-- > 0;)
    {
        double value = rhs[back];
        for (std::size_t column = back + 1; column < n; ++column)
        {
            value -= matrix[back][column] * result.solution[column];
        }

        const double diagonal = matrix[back][back];
        if (std::abs(diagonal) <= 1.0e-20)
        {
            return result;
        }
        result.solution[back] = value / diagonal;
    }

    result.solved = true;
    result.residual_norm = 0.0;
    for (std::size_t row = 0; row < n; ++row)
    {
        double residual = rhs_original[row];
        for (std::size_t column = 0; column < n; ++column)
        {
            residual -= matrix_original[row][column] * result.solution[column];
        }
        result.residual_norm = std::max(result.residual_norm, std::abs(residual));
    }

    return result;
}

IterativeLinearSystemSolveResult solveIterativeLinearSystem(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& rhs,
    const std::vector<double>& initial_guess,
    int max_iterations,
    double tolerance,
    double relaxation)
{
    IterativeLinearSystemSolveResult result;
    const std::size_t n = rhs.size();
    if (n == 0 || matrix.size() != n)
    {
        return result;
    }

    for (const auto& row : matrix)
    {
        if (row.size() != n)
        {
            return result;
        }
    }

    std::size_t sparse_nnz = 0;
    const auto sparse_rows = buildSparseLinearRows(matrix, &sparse_nnz);
    result.nonzeros = sparse_nnz;

    const double effective_tolerance = std::max(1.0e-12, tolerance);
    const double relaxation_factor = std::clamp(relaxation, 0.2, 1.0);
    result.solution =
        initial_guess.size() == n ? initial_guess : std::vector<double>(n, 0.0);

    for (int iteration = 0; iteration < std::max(1, max_iterations); ++iteration)
    {
        double max_delta = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const double diagonal = sparse_rows[i].diagonal;
            if (std::abs(diagonal) < 1.0e-20)
            {
                continue;
            }

            double sigma = rhs[i];
            for (const auto& entry : sparse_rows[i].off_diagonal)
            {
                sigma -= entry.value * result.solution[entry.column];
            }

            const double updated = sigma / diagonal;
            max_delta = std::max(max_delta, std::abs(updated - result.solution[i]));
            result.solution[i] =
                relaxation_factor * updated + (1.0 - relaxation_factor) * result.solution[i];
        }

        result.iterations = iteration + 1;
        if (max_delta <= effective_tolerance)
        {
            result.converged = true;
            break;
        }
    }

    result.residual_norm = sparseLinearResidualNorm(sparse_rows, rhs, result.solution);
    return result;
}

VolumeLinearSolverRouting resolveVolumeLinearSolverRouting(
    const VolumeLinearSolverRoutingInput& input)
{
    VolumeLinearSolverRouting routing;
    routing.auto_prefers_iterative =
        input.system_size >= 6 || input.has_external_cells || input.has_external_cell_faces;

    if (input.iterative_only && !input.dense_only)
    {
        routing.use_iterative_solver = true;
    }
    else if (input.dense_only && !input.iterative_only)
    {
        routing.use_iterative_solver = false;
    }
    else
    {
        routing.use_iterative_solver = routing.auto_prefers_iterative;
    }

    routing.solver_mode = routing.use_iterative_solver ? "iterative" : "dense";
    return routing;
}

} // namespace Solver
} // namespace SCDAT
