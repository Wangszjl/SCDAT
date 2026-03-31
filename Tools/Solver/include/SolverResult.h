#ifndef SCDAT_SOLVER_RESULT_H
#define SCDAT_SOLVER_RESULT_H

#include <string>

namespace SCDAT
{
namespace Solver
{

struct SolverResult
{
    bool converged = false;
    int iterations = 0;
    double residual = 0.0;
    double solve_time = 0.0;
    std::string message;
    std::string error_message;
};

} // namespace Solver
} // namespace SCDAT

#endif // SCDAT_SOLVER_RESULT_H
