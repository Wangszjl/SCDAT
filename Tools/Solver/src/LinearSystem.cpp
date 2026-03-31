#include "../include/LinearSystem.h"
#include "../include/SparseMatrix.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SCDAT
{
namespace Solver
{

void SparseMatrix::addEntry(std::size_t row, std::size_t col, double value)
{
    if (row >= rows_ || col >= cols_)
    {
        throw std::out_of_range("SparseMatrix index out of range");
    }

    row_ptr_[row + 1] += 1;
    col_indices_.push_back(col);
    values_.push_back(value);
}

void SparseMatrix::finalize()
{
    for (std::size_t i = 1; i < row_ptr_.size(); ++i)
    {
        row_ptr_[i] += row_ptr_[i - 1];
    }
}

void LinearSystem::applyDirichletBC(std::size_t node_index, double value)
{
    if (node_index >= size_)
    {
        throw std::out_of_range("LinearSystem Dirichlet index out of range");
    }

    for (std::size_t i = 0; i < size_; ++i)
    {
        if (i != node_index)
        {
            rhs_[i] -= matrix_[i][node_index] * value;
            matrix_[i][node_index] = 0.0;
        }
    }

    for (std::size_t j = 0; j < size_; ++j)
    {
        matrix_[node_index][j] = (j == node_index) ? 1.0 : 0.0;
    }
    rhs_[node_index] = value;
}

void LinearSystem::applyNeumannBC(std::size_t node_index, double flux)
{
    if (node_index >= size_)
    {
        throw std::out_of_range("LinearSystem Neumann index out of range");
    }

    rhs_[node_index] += flux;
}

bool LinearSystem::solve(SolverType solver_type)
{
    if (solver_type == SolverType::DIRECT)
    {
        return solveDirect();
    }
    return solveIterative();
}

bool LinearSystem::solveDirect()
{
    return gaussianElimination();
}

bool LinearSystem::solveIterative(double tolerance, int max_iterations)
{
    return conjugateGradient(tolerance, max_iterations);
}

double LinearSystem::getResidualNorm() const
{
    double sum = 0.0;
    for (std::size_t i = 0; i < size_; ++i)
    {
        double ri = rhs_[i];
        for (std::size_t j = 0; j < size_; ++j)
        {
            ri -= matrix_[i][j] * solution_[j];
        }
        sum += ri * ri;
    }
    return std::sqrt(sum);
}

double LinearSystem::getConditionNumber() const
{
    if (size_ == 0)
    {
        return 0.0;
    }

    double min_diag = std::numeric_limits<double>::max();
    double max_diag = 0.0;
    for (std::size_t i = 0; i < size_; ++i)
    {
        const double diagonal = std::abs(matrix_[i][i]);
        min_diag = std::min(min_diag, diagonal);
        max_diag = std::max(max_diag, diagonal);
    }

    if (min_diag < 1e-30)
    {
        return std::numeric_limits<double>::infinity();
    }
    return max_diag / min_diag;
}

bool LinearSystem::isSymmetric(double tolerance) const
{
    for (std::size_t i = 0; i < size_; ++i)
    {
        for (std::size_t j = i + 1; j < size_; ++j)
        {
            if (std::abs(matrix_[i][j] - matrix_[j][i]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

bool LinearSystem::isPositiveDefinite() const
{
    for (std::size_t i = 0; i < size_; ++i)
    {
        if (matrix_[i][i] <= 0.0)
        {
            return false;
        }
    }
    return true;
}

void LinearSystem::clear()
{
    for (auto& row : matrix_)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    std::fill(rhs_.begin(), rhs_.end(), 0.0);
    std::fill(solution_.begin(), solution_.end(), 0.0);
    assembled_ = false;
}

void LinearSystem::zero()
{
    clear();
}

void LinearSystem::printMatrix() const
{
    for (std::size_t i = 0; i < size_; ++i)
    {
        for (std::size_t j = 0; j < size_; ++j)
        {
            std::cout << std::setw(10) << matrix_[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void LinearSystem::printRHS() const
{
    for (double value : rhs_)
    {
        std::cout << value << '\n';
    }
}

void LinearSystem::printSolution() const
{
    for (double value : solution_)
    {
        std::cout << value << '\n';
    }
}

std::string LinearSystem::toString() const
{
    std::ostringstream stream;
    stream << "LinearSystem(size=" << size_ << ", assembled=" << (assembled_ ? "true" : "false")
           << ", residual=" << getResidualNorm() << ")";
    return stream.str();
}

bool LinearSystem::gaussianElimination()
{
    std::vector<std::vector<double>> matrix = matrix_;
    std::vector<double> rhs = rhs_;

    for (std::size_t k = 0; k < size_; ++k)
    {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < size_; ++i)
        {
            if (std::abs(matrix[i][k]) > std::abs(matrix[pivot][k]))
            {
                pivot = i;
            }
        }

        if (std::abs(matrix[pivot][k]) < 1e-20)
        {
            return false;
        }

        if (pivot != k)
        {
            std::swap(matrix[pivot], matrix[k]);
            std::swap(rhs[pivot], rhs[k]);
        }

        for (std::size_t i = k + 1; i < size_; ++i)
        {
            const double factor = matrix[i][k] / matrix[k][k];
            for (std::size_t j = k; j < size_; ++j)
            {
                matrix[i][j] -= factor * matrix[k][j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }

    for (int i = static_cast<int>(size_) - 1; i >= 0; --i)
    {
        double sum = rhs[static_cast<std::size_t>(i)];
        for (std::size_t j = static_cast<std::size_t>(i) + 1; j < size_; ++j)
        {
            sum -= matrix[static_cast<std::size_t>(i)][j] * solution_[j];
        }
        solution_[static_cast<std::size_t>(i)] =
            sum / matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }

    return true;
}

bool LinearSystem::conjugateGradient(double tolerance, int max_iterations)
{
    std::vector<double> r(size_, 0.0);
    std::vector<double> p(size_, 0.0);
    std::vector<double> Ap(size_, 0.0);
    std::fill(solution_.begin(), solution_.end(), 0.0);

    for (std::size_t i = 0; i < size_; ++i)
    {
        r[i] = rhs_[i];
        p[i] = r[i];
    }

    auto dot = [](const std::vector<double>& lhs, const std::vector<double>& rhs)
    {
        double sum = 0.0;
        for (std::size_t i = 0; i < lhs.size(); ++i)
        {
            sum += lhs[i] * rhs[i];
        }
        return sum;
    };

    double rsold = dot(r, r);
    if (std::sqrt(rsold) < tolerance)
    {
        return true;
    }

    for (int iteration = 0; iteration < max_iterations; ++iteration)
    {
        for (std::size_t i = 0; i < size_; ++i)
        {
            double sum = 0.0;
            for (std::size_t j = 0; j < size_; ++j)
            {
                sum += matrix_[i][j] * p[j];
            }
            Ap[i] = sum;
        }

        const double pAp = dot(p, Ap);
        if (std::abs(pAp) < 1e-30)
        {
            return false;
        }

        const double alpha = rsold / pAp;
        for (std::size_t i = 0; i < size_; ++i)
        {
            solution_[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        const double rsnew = dot(r, r);
        if (std::sqrt(rsnew) < tolerance)
        {
            return true;
        }

        const double beta = rsnew / rsold;
        for (std::size_t i = 0; i < size_; ++i)
        {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }

    return false;
}

bool LinearSystem::bicgstab(double tolerance, int max_iterations)
{
    return conjugateGradient(tolerance, max_iterations);
}

} // namespace Solver
} // namespace SCDAT
