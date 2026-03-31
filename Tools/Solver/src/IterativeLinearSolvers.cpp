#include "IterativeLinearSolvers.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace SCDAT
{
namespace Solver
{

ILUPreconditioner::ILUPreconditioner(int fill_level, double drop_tolerance)
    : fill_level_(fill_level), drop_tolerance_(drop_tolerance), factorized_(false)
{
}

bool ILUPreconditioner::factorize(const BandMatrix& matrix)
{
    const int n = matrix.getSize();
    if (n <= 0)
    {
        factorized_ = false;
        return false;
    }

    L_ = std::make_shared<BandMatrix>(n);
    U_ = std::make_shared<BandMatrix>(n);
    pivot_.resize(static_cast<size_t>(n));

    std::iota(pivot_.begin(), pivot_.end(), 0);

    try
    {
        iluDecomposition(matrix);
        dropSmallEntries();
        factorized_ = true;
    }
    catch (const std::exception&)
    {
        L_.reset();
        U_.reset();
        pivot_.clear();
        factorized_ = false;
    }

    return factorized_;
}

void ILUPreconditioner::solve(const std::vector<double>& rhs, std::vector<double>& solution) const
{
    if (!factorized_ || !L_ || !U_ || rhs.size() != pivot_.size())
    {
        solution = rhs;
        return;
    }

    std::vector<double> y(rhs.size(), 0.0);
    solveForward(rhs, y);
    solveBackward(y, solution);
}

void ILUPreconditioner::solveForward(const std::vector<double>& rhs, std::vector<double>& y) const
{
    const std::size_t n = rhs.size();
    y.assign(n, 0.0);

    if (!factorized_ || !L_ || pivot_.size() != n)
    {
        y = rhs;
        return;
    }

    for (std::size_t i = 0; i < n; ++i)
    {
        double sum = rhs[static_cast<std::size_t>(pivot_[i])];
        for (std::size_t j = 0; j < i; ++j)
        {
            sum -= (*L_)(static_cast<int>(i), static_cast<int>(j)) * y[j];
        }

        const double diagonal = (*L_)(static_cast<int>(i), static_cast<int>(i));
        y[i] = (std::abs(diagonal) > 1e-30) ? (sum / diagonal) : sum;
    }
}

void ILUPreconditioner::solveBackward(const std::vector<double>& y,
                                      std::vector<double>& solution) const
{
    const std::size_t n = y.size();
    solution.assign(n, 0.0);

    if (!factorized_ || !U_)
    {
        solution = y;
        return;
    }

    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = y[static_cast<std::size_t>(i)];
        for (std::size_t j = static_cast<std::size_t>(i + 1); j < n; ++j)
        {
            sum -= (*U_)(i, static_cast<int>(j)) * solution[j];
        }

        const double diagonal = (*U_)(i, i);
        if (std::abs(diagonal) <= 1e-30)
        {
            solution[static_cast<std::size_t>(i)] = 0.0;
        }
        else
        {
            solution[static_cast<std::size_t>(i)] = sum / diagonal;
        }
    }
}

int ILUPreconditioner::getMemoryUsage() const
{
    if (!L_ || !U_)
    {
        return 0;
    }

    const auto n = static_cast<int>(pivot_.size());
    return static_cast<int>(2 * n * n * sizeof(double) + pivot_.size() * sizeof(int));
}

void ILUPreconditioner::printFactors() const
{
    std::cout << "ILUPreconditioner(factorized=" << (factorized_ ? "true" : "false") << ")"
              << std::endl;
}

void ILUPreconditioner::iluDecomposition(const BandMatrix& matrix)
{
    const int n = matrix.getSize();
    std::vector<std::vector<double>> work(static_cast<std::size_t>(n),
                                          std::vector<double>(static_cast<std::size_t>(n), 0.0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = matrix(i, j);
        }
    }

    std::iota(pivot_.begin(), pivot_.end(), 0);

    for (int k = 0; k < n; ++k)
    {
        int pivot_row = k;
        double pivot_abs = std::abs(work[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)]);

        for (int i = k + 1; i < n; ++i)
        {
            const double candidate =
                std::abs(work[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]);
            if (candidate > pivot_abs)
            {
                pivot_abs = candidate;
                pivot_row = i;
            }
        }

        if (pivot_abs <= 1e-14)
        {
            throw std::runtime_error("ILU factorization failed: singular pivot");
        }

        if (pivot_row != k)
        {
            std::swap(work[static_cast<std::size_t>(k)], work[static_cast<std::size_t>(pivot_row)]);
            std::swap(pivot_[static_cast<std::size_t>(k)], pivot_[static_cast<std::size_t>(pivot_row)]);
        }

        const double pivot = work[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)];
        for (int i = k + 1; i < n; ++i)
        {
            double multiplier =
                work[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] / pivot;
            if (std::abs(multiplier) < drop_tolerance_)
            {
                multiplier = 0.0;
            }

            work[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = multiplier;
            for (int j = k + 1; j < n; ++j)
            {
                work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -=
                    multiplier * work[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];

                if (std::abs(work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) <
                    drop_tolerance_)
                {
                    work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0.0;
                }
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (j < i)
            {
                (*L_)(i, j) = work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                (*U_)(i, j) = 0.0;
            }
            else if (j == i)
            {
                (*L_)(i, j) = 1.0;
                (*U_)(i, j) = work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
            else
            {
                (*L_)(i, j) = 0.0;
                (*U_)(i, j) = work[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
        }
    }
}

void ILUPreconditioner::buildLUMatrices(
    const std::vector<std::vector<std::pair<int, double>>>& work_matrix, int n)
{
    if (n <= 0)
    {
        L_.reset();
        U_.reset();
        return;
    }

    if (!L_ || !U_)
    {
        L_ = std::make_shared<BandMatrix>(n);
        U_ = std::make_shared<BandMatrix>(n);
    }

    for (int i = 0; i < n; ++i)
    {
        (*L_)(i, i) = 1.0;
        if (static_cast<std::size_t>(i) >= work_matrix.size())
        {
            continue;
        }

        for (const auto& [j, value] : work_matrix[static_cast<std::size_t>(i)])
        {
            if (j < 0 || j >= n)
            {
                continue;
            }

            if (j < i)
            {
                (*L_)(i, j) = value;
            }
            else
            {
                (*U_)(i, j) = value;
            }
        }
    }
}

void ILUPreconditioner::dropSmallEntries()
{
    if (!L_ || !U_)
    {
        return;
    }

    const int n = static_cast<int>(pivot_.size());
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i != j && std::abs((*L_)(i, j)) < drop_tolerance_)
            {
                (*L_)(i, j) = 0.0;
            }

            if (std::abs((*U_)(i, j)) < drop_tolerance_)
            {
                (*U_)(i, j) = 0.0;
            }
        }

        (*L_)(i, i) = 1.0;
    }
}

ConjugateGradientSolver::ConjugateGradientSolver() = default;

ConjugateGradientSolver::ConjugateGradientSolver(const Parameters& params) : params_(params) {}

void ConjugateGradientSolver::setPreconditioner(ILUPreconditionerPtr preconditioner)
{
    preconditioner_ = std::move(preconditioner);
    params_.use_preconditioner = static_cast<bool>(preconditioner_);
    params_.preconditioner_type = 0;
}

void ConjugateGradientSolver::setDiagonalPreconditioner(const std::vector<double>& diagonal)
{
    diagonal_preconditioner_ = diagonal;
    params_.use_preconditioner = !diagonal_preconditioner_.empty();
    params_.preconditioner_type = 1;
}

bool ConjugateGradientSolver::solve(const BandMatrix& matrix, const std::vector<double>& rhs,
                                    std::vector<double>& solution)
{
    return solve(matrix, rhs, solution, std::vector<double>(rhs.size(), 0.0));
}

bool ConjugateGradientSolver::solve(const BandMatrix& matrix, const std::vector<double>& rhs,
                                    std::vector<double>& solution,
                                    const std::vector<double>& initial_guess)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    solution = initial_guess;
    bool ok = conjugateGradientIteration(matrix, rhs, solution);
    auto t1 = std::chrono::high_resolution_clock::now();
    stats_.solve_time = std::chrono::duration<double>(t1 - t0).count() * 1000.0;
    stats_.converged = ok;
    return ok;
}

bool ConjugateGradientSolver::conjugateGradientIteration(const BandMatrix& matrix,
                                                         const std::vector<double>& rhs,
                                                         std::vector<double>& solution)
{
    const int n = matrix.getSize();
    if (n <= 0 || rhs.size() != static_cast<size_t>(n))
    {
        return false;
    }

    const auto result = Detail::runConjugateGradient(
        rhs, solution, static_cast<std::size_t>(params_.max_iterations), params_.tolerance,
        [&matrix](const std::vector<double>& input, std::vector<double>& output)
        {
            matrix.multiply(input, output);
        },
        [this](const std::vector<double>& residual, std::vector<double>& preconditioned)
        {
            applyPreconditioner(residual, preconditioned);
        },
        [this](std::size_t iter, double residual_norm, const std::vector<double>&)
        {
            stats_.iterations = static_cast<int>(iter);
            stats_.final_residual = residual_norm;

            if (iter > 0 && params_.verbose &&
                ((static_cast<int>(iter) % params_.print_frequency) == 0 || iter == 1))
            {
                printIteration(static_cast<int>(iter), residual_norm);
            }
        });

    stats_.iterations = static_cast<int>(result.iterations);
    stats_.final_residual = result.final_residual;
    stats_.convergence_rate =
        (result.converged && result.iterations > 0 && result.initial_residual > kBreakdownTolerance)
            ? std::pow(result.final_residual / result.initial_residual,
                       1.0 / static_cast<double>(result.iterations))
            : 0.0;

    return result.converged;
}

void ConjugateGradientSolver::applyPreconditioner(const std::vector<double>& r,
                                                  std::vector<double>& z) const
{
    if (!params_.use_preconditioner)
    {
        z = r;
        return;
    }
    if (params_.preconditioner_type == 0 && preconditioner_ && preconditioner_->isFactorized())
    {
        preconditioner_->solve(r, z);
        return;
    }
    if (params_.preconditioner_type == 1 && diagonal_preconditioner_.size() == r.size())
    {
        z.resize(r.size());
        for (size_t i = 0; i < r.size(); ++i)
        {
            double d = diagonal_preconditioner_[i];
            z[i] = (std::abs(d) > 1e-30) ? (r[i] / d) : r[i];
        }
        return;
    }
    z = r;
}

void ConjugateGradientSolver::printIteration(int iter, double residual) const
{
    std::cout << "CG iter " << iter << ", residual=" << std::scientific << residual << std::endl;
}

double ConjugateGradientSolver::computeResidual(const BandMatrix& matrix,
                                                const std::vector<double>& solution,
                                                const std::vector<double>& rhs) const
{
    std::vector<double> Ax(solution.size(), 0.0);
    matrix.multiply(solution, Ax);
    double s = 0.0;
    for (size_t i = 0; i < rhs.size(); ++i)
    {
        double d = rhs[i] - Ax[i];
        s += d * d;
    }
    return std::sqrt(s);
}

void ConjugateGradientSolver::printStatistics() const
{
    std::cout << "ConjugateGradientSolver: iters=" << stats_.iterations
              << ", residual=" << std::scientific << stats_.final_residual
              << ", time_ms=" << std::fixed << std::setprecision(3) << stats_.solve_time
              << ", converged=" << (stats_.converged ? "true" : "false") << std::endl;
}

} // namespace Solver
} // namespace SCDAT
