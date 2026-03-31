#ifndef SCDAT_SOLVER_PRECONDITIONER_H
#define SCDAT_SOLVER_PRECONDITIONER_H

#include "SparseMatrix.h"

namespace SCDAT
{
namespace Solver
{

class Preconditioner
{
  public:
    virtual ~Preconditioner() = default;

    virtual void apply(const std::vector<double>& input, std::vector<double>& output) = 0;
    virtual void setup(const SparseMatrix& matrix) = 0;
};

} // namespace Solver
} // namespace SCDAT

#endif // SCDAT_SOLVER_PRECONDITIONER_H
