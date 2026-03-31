#pragma once

#include "../../../../Basic/include/ErrorHandling.h"
#include "../../../../Basic/include/Logger.h"
#include "../../ConjugateGradientCommon.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace SCDAT
{
namespace Solver
{
namespace Experimental
{

template <typename T>
concept Numeric = std::is_arithmetic_v<T> && requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
};

template <typename T>
concept VectorLike = requires(T v, std::size_t i) {
    { v.size() } -> std::convertible_to<std::size_t>;
    { v[i] } -> std::convertible_to<typename T::value_type>;
    typename T::value_type;
} && Numeric<typename T::value_type>;

template <typename T>
concept MatrixLike = requires(T m, std::size_t i, std::size_t j) {
    { m.rows() } -> std::convertible_to<std::size_t>;
    { m.cols() } -> std::convertible_to<std::size_t>;
    { m(i, j) } -> std::convertible_to<typename T::value_type>;
    typename T::value_type;
} && Numeric<typename T::value_type>;

template <typename T>
concept SolverLike = requires(T solver) {
    { solver.solve() } -> std::convertible_to<Core::Result<typename T::Vector>>;
    { solver.getIterationCount() } -> std::convertible_to<std::size_t>;
    { solver.getResidualNorm() } -> std::convertible_to<double>;
    { solver.isConverged() } -> std::convertible_to<bool>;
};

enum class SolverType
{
    DIRECT_LU,
    DIRECT_CHOLESKY,
    ITERATIVE_CG,
    ITERATIVE_BICGSTAB,
    ITERATIVE_GMRES,
    MULTIGRID_V,
    MULTIGRID_W,
    ADAPTIVE_HYBRID
};

enum class PreconditionerType
{
    NONE,
    JACOBI,
    GAUSS_SEIDEL,
    ILU,
    MULTIGRID,
    ADAPTIVE
};

struct SolverConfiguration
{
    SolverType solver_type = SolverType::ITERATIVE_CG;
    PreconditionerType preconditioner_type = PreconditionerType::ILU;

    double tolerance = 1e-12;
    std::size_t max_iterations = 10000;
    double relaxation_factor = 1.0;

    bool enable_parallel = true;
    std::size_t num_threads = 0;

    bool enable_monitoring = true;
    std::size_t monitoring_interval = 100;

    bool enable_adaptive = true;
    double adaptive_threshold = 1e-3;
};

struct SolverStatus
{
    bool is_converged = false;
    std::size_t iteration_count = 0;
    double residual_norm = 0.0;
    double solution_norm = 0.0;
    std::chrono::milliseconds solve_time{0};
    std::string convergence_reason;
};

template <VectorLike VectorType, MatrixLike MatrixType> class ModernSolver
{
  public:
    using Vector = VectorType;
    using Matrix = MatrixType;
    using Scalar = typename VectorType::value_type;

    explicit ModernSolver(SolverConfiguration config = {}) : config_(std::move(config)) {}
    virtual ~ModernSolver() = default;

    ModernSolver(const ModernSolver&) = delete;
    ModernSolver& operator=(const ModernSolver&) = delete;
    ModernSolver(ModernSolver&&) = default;
    ModernSolver& operator=(ModernSolver&&) = default;

    virtual Core::VoidResult setSystem(const Matrix& A, const Vector& b) = 0;
    virtual Core::Result<Vector> solve() = 0;
    virtual Core::VoidResult solve(Vector& x) = 0;

    virtual std::vector<SolverStatus> solveIteratively()
    {
        return {};
    }

    [[nodiscard]] const SolverStatus& getStatus() const noexcept
    {
        return status_;
    }

    [[nodiscard]] const SolverConfiguration& getConfiguration() const noexcept
    {
        return config_;
    }

    [[nodiscard]] bool isConverged() const noexcept
    {
        return status_.is_converged;
    }

    [[nodiscard]] std::size_t getIterationCount() const noexcept
    {
        return status_.iteration_count;
    }

    [[nodiscard]] double getResidualNorm() const noexcept
    {
        return status_.residual_norm;
    }

    void setConfiguration(SolverConfiguration config)
    {
        config_ = std::move(config);
    }

    void setTolerance(double tolerance)
    {
        config_.tolerance = tolerance;
    }

    void setMaxIterations(std::size_t max_iterations)
    {
        config_.max_iterations = max_iterations;
    }

    using MonitorCallback = std::function<void(const SolverStatus&)>;
    void setMonitorCallback(MonitorCallback callback)
    {
        monitor_callback_ = std::move(callback);
    }

  protected:
    SolverConfiguration config_;
    SolverStatus status_;
    MonitorCallback monitor_callback_;

    void updateStatus(std::size_t iteration, double residual_norm, double solution_norm = 0.0)
    {
        status_.iteration_count = iteration;
        status_.residual_norm = residual_norm;
        status_.solution_norm = solution_norm;
        status_.is_converged = (residual_norm <= config_.tolerance);

        if (monitor_callback_ && config_.enable_monitoring &&
            (iteration % config_.monitoring_interval == 0))
        {
            monitor_callback_(status_);
        }
    }

    void logProgress(const std::string& message) const
    {
        if (config_.enable_monitoring)
        {
            LOG_DEBUG("Experimental::ModernSolver", "{}", message);
        }
    }
};

template <VectorLike VectorType, MatrixLike MatrixType>
class ConjugateGradientSolver : public ModernSolver<VectorType, MatrixType>
{
  public:
    using Base = ModernSolver<VectorType, MatrixType>;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using Scalar = typename Base::Scalar;

    explicit ConjugateGradientSolver(SolverConfiguration config = {}) : Base(std::move(config)) {}

    Core::VoidResult setSystem(const Matrix& A, const Vector& b) override
    {
        if (A.rows() != A.cols())
        {
            return Core::ErrorHandler::makeError<void>(Core::ErrorCode::SOLVER_MATRIX_SINGULAR);
        }

        if (A.rows() != b.size())
        {
            return Core::ErrorHandler::makeError<void>(Core::ErrorCode::INVALID_ARGUMENT);
        }

        A_ = &A;
        b_ = &b;

        LOG_INFO("Experimental::ConjugateGradientSolver", "system set");
        return Core::ErrorHandler::makeSuccess();
    }

    Core::Result<Vector> solve() override
    {
        if (!A_ || !b_)
        {
            return Core::ErrorHandler::makeError<Vector>(Core::ErrorCode::INVALID_ARGUMENT);
        }

        Vector x(b_->size());
        std::fill(x.begin(), x.end(), Scalar{0});

        auto result = solve(x);
        if (!result.has_value())
        {
            return Core::ErrorHandler::makeError<Vector>(result.error());
        }

        return Core::ErrorHandler::makeSuccess(std::move(x));
    }

    Core::VoidResult solve(Vector& x) override
    {
        if (!A_ || !b_ || x.size() != b_->size())
        {
            return Core::ErrorHandler::makeError<void>(Core::ErrorCode::INVALID_ARGUMENT);
        }

        auto start_time = std::chrono::high_resolution_clock::now();

        const auto result = ::SCDAT::Solver::Detail::runConjugateGradient(
            *b_, x, this->config_.max_iterations, this->config_.tolerance,
            [this](const Vector& input, Vector& output)
            {
                matrixVectorProduct(*A_, input, output);
            },
            [](const Vector& residual, Vector& preconditioned)
            {
                preconditioned = residual;
            },
            [this](std::size_t iter, double residual_norm, const Vector& solution)
            {
                this->updateStatus(iter, residual_norm,
                                   ::SCDAT::Solver::Detail::norm(solution));
            });

        auto end_time = std::chrono::high_resolution_clock::now();
        this->status_.solve_time =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        if (!result.converged)
        {
            this->status_.convergence_reason =
                (result.stop_reason == ::SCDAT::Solver::Detail::ConjugateGradientStopReason::SINGULAR_MATRIX)
                    ? "singular or ill-conditioned matrix"
                    : "reached max iterations";

            return Core::ErrorHandler::makeError<void>(
                result.stop_reason == ::SCDAT::Solver::Detail::ConjugateGradientStopReason::SINGULAR_MATRIX
                    ? Core::ErrorCode::SOLVER_MATRIX_SINGULAR
                    : Core::ErrorCode::SOLVER_CONVERGENCE_FAILED);
        }

        this->status_.convergence_reason = "reached tolerance";
        LOG_INFO("Experimental::ConjugateGradientSolver", "solve finished");
        return Core::ErrorHandler::makeSuccess();
    }

    std::vector<SolverStatus> solveIteratively() override
    {
        std::vector<SolverStatus> statuses;
        if (!A_ || !b_)
        {
            return statuses;
        }

        Vector x(b_->size());
        std::fill(x.begin(), x.end(), Scalar{0});

        const auto result = ::SCDAT::Solver::Detail::runConjugateGradient(
            *b_, x, this->config_.max_iterations, this->config_.tolerance,
            [this](const Vector& input, Vector& output)
            {
                matrixVectorProduct(*A_, input, output);
            },
            [](const Vector& residual, Vector& preconditioned)
            {
                preconditioned = residual;
            },
            [this, &statuses](std::size_t iter, double residual_norm, const Vector& solution)
            {
                this->updateStatus(iter, residual_norm,
                                   ::SCDAT::Solver::Detail::norm(solution));
                statuses.push_back(this->status_);
            });

        if (result.converged)
        {
            this->status_.convergence_reason = "reached tolerance";
        }
        else if (result.stop_reason ==
                 ::SCDAT::Solver::Detail::ConjugateGradientStopReason::SINGULAR_MATRIX)
        {
            this->status_.convergence_reason = "singular or ill-conditioned matrix";
        }
        else if (result.stop_reason ==
                 ::SCDAT::Solver::Detail::ConjugateGradientStopReason::BREAKDOWN)
        {
            this->status_.convergence_reason = "breakdown in residual update";
        }
        else
        {
            this->status_.convergence_reason = "reached max iterations";
        }

        return statuses;
    }

  private:
    const Matrix* A_ = nullptr;
    const Vector* b_ = nullptr;

    void matrixVectorProduct(const Matrix& A, const Vector& x, Vector& y) const
    {
        for (std::size_t i = 0; i < A.rows(); ++i)
        {
            y[i] = Scalar{0};
            for (std::size_t j = 0; j < A.cols(); ++j)
            {
                y[i] += A(i, j) * x[j];
            }
        }
    }
};

} // namespace Experimental
} // namespace Solver
} // namespace SCDAT
