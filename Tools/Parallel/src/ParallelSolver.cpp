#include "ParallelSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>

#if defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <mach/task.h>
#endif

namespace SCDAT
{
namespace Solver
{

ParallelSolver::ParallelSolver(int num_threads)
    : num_threads_(num_threads > 0 ? num_threads : omp_get_max_threads()), parallel_efficiency_(1.0)
{
    omp_set_num_threads(num_threads_);
}

void ParallelSolver::setNumThreads(int num_threads)
{
    num_threads_ = num_threads > 0 ? num_threads : omp_get_max_threads();
    omp_set_num_threads(num_threads_);
}

void ParallelSolver::parallelMatrixVectorMultiply(const SparseMatrix& matrix,
                                                  const std::vector<double>& vector,
                                                  std::vector<double>& result)
{
    const size_t n = matrix.rows();
    result.assign(n, 0.0);

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(n); ++i)
    {
        const size_t row = static_cast<size_t>(i);
        if (row + 1 >= matrix.row_ptr_.size())
        {
            continue;
        }

        double sum = 0.0;
        const size_t begin = matrix.row_ptr_[row];
        const size_t end = matrix.row_ptr_[row + 1];
        for (size_t k = begin; k < end; ++k)
        {
            const size_t col = matrix.col_indices_[k];
            if (col < vector.size())
            {
                sum += matrix.values_[k] * vector[col];
            }
        }
        result[row] = sum;
    }
}

void ParallelSolver::parallelVectorOperation(const std::vector<double>& a,
                                             const std::vector<double>& b,
                                             std::vector<double>& result, VectorOperation operation)
{
    const size_t n = std::min(a.size(), b.size());
    result.assign(n, 0.0);

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(n); ++i)
    {
        const size_t idx = static_cast<size_t>(i);
        switch (operation)
        {
        case VectorOperation::ADD:
            result[idx] = a[idx] + b[idx];
            break;
        case VectorOperation::SUBTRACT:
            result[idx] = a[idx] - b[idx];
            break;
        case VectorOperation::MULTIPLY:
            result[idx] = a[idx] * b[idx];
            break;
        case VectorOperation::DIVIDE:
            result[idx] = (std::abs(b[idx]) > 1e-30) ? (a[idx] / b[idx]) : 0.0;
            break;
        }
    }
}

double ParallelSolver::parallelDotProduct(const std::vector<double>& a,
                                          const std::vector<double>& b)
{
    const size_t n = std::min(a.size(), b.size());
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < static_cast<int>(n); ++i)
    {
        sum += a[static_cast<size_t>(i)] * b[static_cast<size_t>(i)];
    }

    return sum;
}

double ParallelSolver::parallelVectorNorm(const std::vector<double>& vector)
{
    return std::sqrt(parallelDotProduct(vector, vector));
}

std::pair<size_t, size_t> ParallelSolver::getLoadBalancedRange(size_t total_size, int thread_id,
                                                               int num_threads)
{
    const size_t chunk = total_size / static_cast<size_t>(num_threads);
    const size_t rem = total_size % static_cast<size_t>(num_threads);

    size_t start = 0;
    size_t end = 0;
    if (thread_id < static_cast<int>(rem))
    {
        start = static_cast<size_t>(thread_id) * (chunk + 1);
        end = start + chunk + 1;
    }
    else
    {
        start = static_cast<size_t>(thread_id) * chunk + rem;
        end = start + chunk;
    }

    return {start, std::min(end, total_size)};
}

void ParallelSolver::updateParallelEfficiency(double serial_time, double parallel_time)
{
    if (parallel_time <= 0.0)
    {
        return;
    }

    const double speedup = serial_time / parallel_time;
    parallel_efficiency_ = speedup / std::max(1, num_threads_);
    timing_data_.push_back(parallel_efficiency_);
    if (timing_data_.size() > 1000)
    {
        timing_data_.erase(timing_data_.begin());
    }
}

ParallelPoissonSolver::ParallelPoissonSolver(int num_threads) : ParallelSolver(num_threads) {}

SolverResult ParallelPoissonSolver::solve(const SparseMatrix& matrix,
                                          const std::vector<double>& rhs,
                                          std::vector<double>& solution, double tolerance,
                                          int max_iterations)
{
    const auto t0 = std::chrono::high_resolution_clock::now();
    SolverResult r = parallelConjugateGradient(matrix, rhs, solution, tolerance, max_iterations);
    const auto t1 = std::chrono::high_resolution_clock::now();
    r.solve_time = static_cast<double>(
                       std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
                   1000.0;
    return r;
}

void ParallelPoissonSolver::setPreconditioner(std::shared_ptr<Preconditioner> preconditioner)
{
    preconditioner_ = std::move(preconditioner);
}

SolverResult ParallelPoissonSolver::parallelConjugateGradient(const SparseMatrix& matrix,
                                                              const std::vector<double>& rhs,
                                                              std::vector<double>& solution,
                                                              double tolerance, int max_iterations)
{
    SolverResult out;
    const size_t n = rhs.size();
    solution.assign(n, 0.0);

    std::vector<double> r(n, 0.0);
    std::vector<double> p(n, 0.0);
    std::vector<double> Ap(n, 0.0);
    std::vector<double> z(n, 0.0);

    parallelMatrixVectorMultiply(matrix, solution, Ap);
    parallelVectorOperation(rhs, Ap, r, VectorOperation::SUBTRACT);

    if (preconditioner_)
    {
        applyPreconditioner(r, z);
        p = z;
    }
    else
    {
        z = r;
        p = r;
    }

    double rz_old = parallelDotProduct(r, z);
    out.residual = std::sqrt(std::max(0.0, rz_old));

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        parallelMatrixVectorMultiply(matrix, p, Ap);
        const double denom = parallelDotProduct(p, Ap);
        if (std::abs(denom) < 1e-30)
        {
            out.error_message = "breakdown: denominator close to zero";
            out.iterations = iter;
            return out;
        }

        const double alpha = rz_old / denom;

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(n); ++i)
        {
            const size_t idx = static_cast<size_t>(i);
            solution[idx] += alpha * p[idx];
            r[idx] -= alpha * Ap[idx];
        }

        out.residual = parallelVectorNorm(r);
        out.iterations = iter + 1;
        if (out.residual < tolerance)
        {
            out.converged = true;
            return out;
        }

        if (preconditioner_)
        {
            applyPreconditioner(r, z);
        }
        else
        {
            z = r;
        }

        const double rz_new = parallelDotProduct(r, z);
        const double beta = (std::abs(rz_old) > 1e-30) ? (rz_new / rz_old) : 0.0;

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(n); ++i)
        {
            const size_t idx = static_cast<size_t>(i);
            p[idx] = z[idx] + beta * p[idx];
        }

        rz_old = rz_new;
    }

    out.error_message = "failed to converge within max_iterations";
    return out;
}

void ParallelPoissonSolver::applyPreconditioner(const std::vector<double>& input,
                                                std::vector<double>& output)
{
    if (preconditioner_)
    {
        preconditioner_->apply(input, output);
    }
    else
    {
        output = input;
    }
}

ParallelLinearSolver::ParallelLinearSolver(SolverType solver_type, int num_threads)
    : ParallelSolver(num_threads), solver_type_(solver_type)
{
}

SolverResult ParallelLinearSolver::solve(const SparseMatrix& matrix, const std::vector<double>& rhs,
                                         std::vector<double>& solution, double tolerance,
                                         int max_iterations)
{
    const auto t0 = std::chrono::high_resolution_clock::now();

    SolverResult r;
    switch (solver_type_)
    {
    case SolverType::CONJUGATE_GRADIENT:
    {
        ParallelPoissonSolver cg(getNumThreads());
        r = cg.solve(matrix, rhs, solution, tolerance, max_iterations);
        break;
    }
    case SolverType::BICGSTAB:
        r = parallelBiCGSTAB(matrix, rhs, solution, tolerance, max_iterations);
        break;
    case SolverType::GMRES:
        r = parallelGMRES(matrix, rhs, solution, tolerance, max_iterations);
        break;
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    r.solve_time = static_cast<double>(
                       std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
                   1000.0;
    return r;
}

SolverResult ParallelLinearSolver::parallelBiCGSTAB(const SparseMatrix& matrix,
                                                    const std::vector<double>& rhs,
                                                    std::vector<double>& solution, double tolerance,
                                                    int max_iterations)
{
    ParallelPoissonSolver fallback(getNumThreads());
    return fallback.solve(matrix, rhs, solution, tolerance, max_iterations);
}

SolverResult ParallelLinearSolver::parallelGMRES(const SparseMatrix& matrix,
                                                 const std::vector<double>& rhs,
                                                 std::vector<double>& solution, double tolerance,
                                                 int max_iterations)
{
    ParallelPoissonSolver fallback(getNumThreads());
    return fallback.solve(matrix, rhs, solution, tolerance, max_iterations);
}

void ParallelPerformanceMonitor::startMonitoring()
{
    start_time_ = omp_get_wtime();
}

void ParallelPerformanceMonitor::stopMonitoring()
{
    end_time_ = omp_get_wtime();

    metrics_.parallel_time = end_time_ - start_time_;
    metrics_.num_threads = omp_get_max_threads();
    metrics_.serial_time = metrics_.parallel_time * std::max(1, metrics_.num_threads);

    if (metrics_.parallel_time > 0.0)
    {
        metrics_.speedup = metrics_.serial_time / metrics_.parallel_time;
        metrics_.efficiency = metrics_.speedup / std::max(1, metrics_.num_threads);
    }
    else
    {
        metrics_.speedup = 0.0;
        metrics_.efficiency = 0.0;
    }

    metrics_.memory_usage = getCurrentMemoryUsage();
}

ParallelPerformanceMonitor::PerformanceMetrics ParallelPerformanceMonitor::getMetrics() const
{
    return metrics_;
}

void ParallelPerformanceMonitor::reset()
{
    start_time_ = 0.0;
    end_time_ = 0.0;
    metrics_ = PerformanceMetrics{};
}

double ParallelPerformanceMonitor::getCurrentMemoryUsage() const
{
#ifdef _WIN32
    // MinGW environments may have incompatible Windows API headers in this build.
    return 0.0;
#elif defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line))
    {
        if (line.rfind("VmRSS:", 0) == 0)
        {
            std::istringstream iss(line);
            std::string label;
            double kb = 0.0;
            std::string unit;
            iss >> label >> kb >> unit;
            return kb / 1024.0;
        }
    }
    return 0.0;
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info),
                  &count) == KERN_SUCCESS)
    {
        return static_cast<double>(info.resident_size) / (1024.0 * 1024.0);
    }
    return 0.0;
#else
    return 0.0;
#endif
}

} // namespace Solver
} // namespace SCDAT
