#ifndef SCDAT_SOLVER_LINEAR_SYSTEM_H
#define SCDAT_SOLVER_LINEAR_SYSTEM_H

#include <cstddef>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Solver
{

enum class SolverType
{
    DIRECT = 1,
    ITERATIVE = 2,
    MULTIGRID = 3
};

/**
 * @brief Dense linear system assembly object used by field-side equation assemblers.
 */
class LinearSystem
{
  public:
    explicit LinearSystem(std::size_t size)
        : size_(size), matrix_(size, std::vector<double>(size, 0.0)), rhs_(size, 0.0),
          solution_(size, 0.0), assembled_(false)
    {
    }

    virtual ~LinearSystem() = default;

    std::size_t getSize() const { return size_; }
    bool isAssembled() const { return assembled_; }

    double& operator()(std::size_t i, std::size_t j) { return matrix_[i][j]; }
    const double& operator()(std::size_t i, std::size_t j) const { return matrix_[i][j]; }

    double& getRHS(std::size_t i) { return rhs_[i]; }
    const double& getRHS(std::size_t i) const { return rhs_[i]; }

    double& getSolution(std::size_t i) { return solution_[i]; }
    const double& getSolution(std::size_t i) const { return solution_[i]; }

    void addToMatrix(std::size_t i, std::size_t j, double value) { matrix_[i][j] += value; }
    void addToRHS(std::size_t i, double value) { rhs_[i] += value; }
    void setMatrixEntry(std::size_t i, std::size_t j, double value) { matrix_[i][j] = value; }
    void setRHSEntry(std::size_t i, double value) { rhs_[i] = value; }

    void applyDirichletBC(std::size_t node_index, double value);
    void applyNeumannBC(std::size_t node_index, double flux);

    virtual bool solve(SolverType solver_type = SolverType::DIRECT);
    virtual bool solveDirect();
    virtual bool solveIterative(double tolerance = 1e-6, int max_iterations = 1000);

    double getResidualNorm() const;
    double getConditionNumber() const;
    bool isSymmetric(double tolerance = 1e-12) const;
    bool isPositiveDefinite() const;

    void clear();
    void zero();
    void markAssembled() { assembled_ = true; }

    void printMatrix() const;
    void printRHS() const;
    void printSolution() const;
    std::string toString() const;

  protected:
    std::size_t size_;
    std::vector<std::vector<double>> matrix_;
    std::vector<double> rhs_;
    std::vector<double> solution_;
    bool assembled_;

    bool gaussianElimination();
    bool conjugateGradient(double tolerance, int max_iterations);
    bool bicgstab(double tolerance, int max_iterations);
};

} // namespace Solver
} // namespace SCDAT

#endif // SCDAT_SOLVER_LINEAR_SYSTEM_H
