#include "../include/CouplingBase.h"

namespace SCDAT
{
namespace Coupling
{

CouplingBase::CouplingBase(std::string name, CouplingType type)
    : name_(std::move(name)), type_(type)
{
}

const std::string& CouplingBase::getName() const
{
    return name_;
}

CouplingType CouplingBase::getType() const
{
    return type_;
}

void CouplingBase::setConvergenceCriteria(const ConvergenceCriteria& criteria)
{
    criteria_ = criteria;
}

const ConvergenceCriteria& CouplingBase::getConvergenceCriteria() const
{
    return criteria_;
}

const CouplingStatistics& CouplingBase::getStatistics() const
{
    return statistics_;
}

void CouplingBase::updateStatistics(int iteration, double relative_error, double absolute_error,
                                    bool converged, const std::string& reason)
{
    statistics_.total_iterations = iteration;
    statistics_.final_relative_error = relative_error;
    statistics_.final_absolute_error = absolute_error;
    statistics_.converged = converged;
    statistics_.convergence_reason = reason;
}

FunctionalCoupling::FunctionalCoupling(std::string name, CouplingType type,
                                       StepFunction step_function, ErrorFunction error_function)
    : CouplingBase(std::move(name), type),
      step_function_(std::move(step_function)),
      error_function_(std::move(error_function))
{
}

VoidResult FunctionalCoupling::execute(double dt)
{
    if (!step_function_)
    {
        return VoidResult::failure("FunctionalCoupling step function is not set");
    }
    return step_function_(dt);
}

std::pair<double, double> FunctionalCoupling::currentErrors() const
{
    if (!error_function_)
    {
        return {0.0, 0.0};
    }
    return error_function_();
}

} // namespace Coupling
} // namespace SCDAT
