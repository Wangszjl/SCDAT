#include "FluidAlgorithmManager.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

void FluidAlgorithmManager::addRegion(const std::string& name,
                                      std::shared_ptr<FluidAlgorithmAdapter> adapter)
{
    if (adapter)
    {
        regions_[name] = std::move(adapter);
    }
}

bool FluidAlgorithmManager::advanceAll(double dt) const
{
    for (const auto& [_, adapter] : regions_)
    {
        if (adapter && !adapter->advance(dt))
        {
            return false;
        }
    }
    return true;
}

std::shared_ptr<FluidAlgorithmAdapter> FluidAlgorithmManager::getRegion(
    const std::string& name) const
{
    if (const auto it = regions_.find(name); it != regions_.end())
    {
        return it->second;
    }
    return {};
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
