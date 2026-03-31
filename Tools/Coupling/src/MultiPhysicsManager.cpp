#include "../include/MultiPhysicsManager.h"

namespace SCDAT
{
namespace Coupling
{

void MultiPhysicsManager::addCoupling(std::shared_ptr<CouplingBase> coupling)
{
    if (coupling)
    {
        couplings_.push_back(std::move(coupling));
    }
}

bool MultiPhysicsManager::hasCoupling(const std::string& name) const
{
    return getCoupling(name) != nullptr;
}

std::shared_ptr<CouplingBase> MultiPhysicsManager::getCoupling(const std::string& name) const
{
    for (const auto& coupling : couplings_)
    {
        if (coupling && coupling->getName() == name)
        {
            return coupling;
        }
    }
    return {};
}

VoidResult MultiPhysicsManager::executeAllCouplings(double dt) const
{
    for (const auto& coupling : couplings_)
    {
        if (!coupling)
        {
            continue;
        }

        if (auto result = coupling->execute(dt); !result)
        {
            return result;
        }

        const auto [relative_error, absolute_error] = coupling->currentErrors();
        coupling->updateStatistics(1, relative_error, absolute_error, true, "one-shot");
    }

    return VoidResult::success();
}

VoidResult MultiPhysicsManager::executeIterativeCouplings(double dt) const
{
    for (const auto& coupling : couplings_)
    {
        if (!coupling)
        {
            continue;
        }

        const auto& criteria = coupling->getConvergenceCriteria();
        bool converged = false;
        for (int iteration = 1; iteration <= criteria.max_iterations; ++iteration)
        {
            if (auto result = coupling->execute(dt); !result)
            {
                return result;
            }

            const auto [relative_error, absolute_error] = coupling->currentErrors();
            converged = criteria.isConverged(relative_error, absolute_error, iteration);
            coupling->updateStatistics(iteration, relative_error, absolute_error, converged,
                                       converged ? "criteria-met" : "running");
            if (converged)
            {
                break;
            }
        }
    }

    return VoidResult::success();
}

std::unordered_map<std::string, CouplingStatistics> MultiPhysicsManager::getAllStatistics() const
{
    std::unordered_map<std::string, CouplingStatistics> stats;
    for (const auto& coupling : couplings_)
    {
        if (coupling)
        {
            stats[coupling->getName()] = coupling->getStatistics();
        }
    }
    return stats;
}

} // namespace Coupling
} // namespace SCDAT
