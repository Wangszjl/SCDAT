#pragma once

#include "../../Basic/include/VoidResult.h"

#include <functional>
#include <string>
#include <utility>

namespace SCDAT
{
namespace Coupling
{

enum class CouplingType
{
    ONE_WAY = 0,
    TWO_WAY = 1,
    ITERATIVE = 2
};

struct ConvergenceCriteria
{
    double relative_tolerance = 1.0e-6;
    double absolute_tolerance = 1.0e-9;
    int max_iterations = 50;
    int min_iterations = 1;

    bool isConverged(double relative_error, double absolute_error, int iteration) const
    {
        if (iteration < min_iterations)
        {
            return false;
        }
        return relative_error <= relative_tolerance && absolute_error <= absolute_tolerance;
    }
};

struct CouplingStatistics
{
    int total_iterations = 0;
    double final_relative_error = 0.0;
    double final_absolute_error = 0.0;
    bool converged = false;
    std::string convergence_reason;
};

class CouplingBase
{
  public:
    CouplingBase(std::string name, CouplingType type);
    virtual ~CouplingBase() = default;

    const std::string& getName() const;
    CouplingType getType() const;

    void setConvergenceCriteria(const ConvergenceCriteria& criteria);
    const ConvergenceCriteria& getConvergenceCriteria() const;
    const CouplingStatistics& getStatistics() const;
    void updateStatistics(int iteration, double relative_error, double absolute_error,
                          bool converged, const std::string& reason);

    virtual VoidResult execute(double dt) = 0;
    virtual std::pair<double, double> currentErrors() const = 0;

  private:
    std::string name_;
    CouplingType type_ = CouplingType::ONE_WAY;
    ConvergenceCriteria criteria_;
    CouplingStatistics statistics_;
};

class FunctionalCoupling : public CouplingBase
{
  public:
    using StepFunction = std::function<VoidResult(double)>;
    using ErrorFunction = std::function<std::pair<double, double>()>;

    FunctionalCoupling(std::string name, CouplingType type, StepFunction step_function,
                       ErrorFunction error_function);

    VoidResult execute(double dt) override;
    std::pair<double, double> currentErrors() const override;

  private:
    StepFunction step_function_;
    ErrorFunction error_function_;
};

} // namespace Coupling
} // namespace SCDAT
