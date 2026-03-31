#ifndef SCDAT_SOLVER_SPARSE_MATRIX_H
#define SCDAT_SOLVER_SPARSE_MATRIX_H

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace SCDAT
{
namespace Solver
{

enum class VectorOperation
{
    ADD = 1,
    SUBTRACT = 2,
    MULTIPLY = 3,
    DIVIDE = 4
};

/**
 * @brief CSR sparse matrix used by solver-side linear algebra kernels.
 */
class SparseMatrix
{
  public:
    SparseMatrix() : rows_(0), cols_(0) {}

    SparseMatrix(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols)
    {
        row_ptr_.resize(rows + 1, 0);
    }

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t nnz() const { return values_.size(); }

    void addEntry(std::size_t row, std::size_t col, double value);
    void finalize();

    std::vector<std::size_t> row_ptr_;
    std::vector<std::size_t> col_indices_;
    std::vector<double> values_;

  private:
    std::size_t rows_;
    std::size_t cols_;
};

} // namespace Solver
} // namespace SCDAT

#endif // SCDAT_SOLVER_SPARSE_MATRIX_H
