#pragma once

#include <cmath>
#include <cstddef>
#include <limits>

namespace SCDAT
{
namespace Solver
{
namespace Detail
{

enum class ConjugateGradientStopReason
{
    CONVERGED,
    MAX_ITERATIONS,
    SINGULAR_MATRIX,
    BREAKDOWN,
    INVALID_SYSTEM
};

struct ConjugateGradientResult
{
    bool converged = false;
    std::size_t iterations = 0;
    double initial_residual = 0.0;
    double final_residual = 0.0;
    ConjugateGradientStopReason stop_reason = ConjugateGradientStopReason::INVALID_SYSTEM;
};

template <typename Vector>
auto dotProduct(const Vector& a, const Vector& b) -> typename Vector::value_type
{
    using Scalar = typename Vector::value_type;
    Scalar result{0};
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

template <typename Vector> double norm(const Vector& x)
{
    return std::sqrt(static_cast<double>(dotProduct(x, x)));
}

template <typename Scalar, typename Vector>
void axpy(Scalar alpha, const Vector& x, Vector& y)
{
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        y[i] += alpha * x[i];
    }
}

template <typename Scalar, typename Vector> void scale(Scalar alpha, Vector& x)
{
    for (auto& value : x)
    {
        value *= alpha;
    }
}

template <typename Vector> void addInPlace(const Vector& x, Vector& y)
{
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        y[i] += x[i];
    }
}

template <typename Vector, typename MatrixMultiplyFn, typename PreconditionerFn,
          typename IterationHookFn>
ConjugateGradientResult runConjugateGradient(const Vector& rhs, Vector& solution,
                                             std::size_t max_iterations, double tolerance,
                                             MatrixMultiplyFn&& multiply,
                                             PreconditionerFn&& apply_preconditioner,
                                             IterationHookFn&& on_iteration)
{
    using Scalar = typename Vector::value_type;

    ConjugateGradientResult result;
    if (rhs.size() != solution.size())
    {
        return result;
    }

    Vector r(rhs.size());
    Vector p(rhs.size());
    Vector Ap(rhs.size());
    Vector z(rhs.size());

    multiply(solution, Ap);
    for (std::size_t i = 0; i < rhs.size(); ++i)
    {
        r[i] = rhs[i] - Ap[i];
    }

    apply_preconditioner(r, z);
    p = z;

    Scalar rz_old = dotProduct(r, z);
    result.initial_residual = norm(r);
    result.final_residual = result.initial_residual;

    on_iteration(0, result.initial_residual, solution);

    if (result.initial_residual <= tolerance)
    {
        result.converged = true;
        result.stop_reason = ConjugateGradientStopReason::CONVERGED;
        return result;
    }

    for (std::size_t iter = 1; iter <= max_iterations; ++iter)
    {
        multiply(p, Ap);

        const Scalar pAp = dotProduct(p, Ap);
        if (std::abs(pAp) <= std::numeric_limits<Scalar>::epsilon())
        {
            result.iterations = iter - 1;
            result.stop_reason = ConjugateGradientStopReason::SINGULAR_MATRIX;
            return result;
        }

        const Scalar alpha = rz_old / pAp;
        axpy(alpha, p, solution);
        axpy(-alpha, Ap, r);

        apply_preconditioner(r, z);
        const Scalar rz_new = dotProduct(r, z);

        result.iterations = iter;
        result.final_residual = norm(r);
        on_iteration(iter, result.final_residual, solution);

        if (result.final_residual <= tolerance)
        {
            result.converged = true;
            result.stop_reason = ConjugateGradientStopReason::CONVERGED;
            return result;
        }

        if (std::abs(rz_old) <= std::numeric_limits<Scalar>::epsilon())
        {
            result.stop_reason = ConjugateGradientStopReason::BREAKDOWN;
            return result;
        }

        const Scalar beta = rz_new / rz_old;
        scale(beta, p);
        addInPlace(z, p);
        rz_old = rz_new;
    }

    result.stop_reason = ConjugateGradientStopReason::MAX_ITERATIONS;
    return result;
}

} // namespace Detail
} // namespace Solver
} // namespace SCDAT
