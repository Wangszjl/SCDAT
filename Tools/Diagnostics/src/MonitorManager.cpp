#include "../include/MonitorManager.h"

namespace SCDAT
{
namespace Diagnostics
{

void MonitorManager::addMonitor(std::shared_ptr<Monitor> monitor)
{
    if (monitor)
    {
        monitors_.push_back(std::move(monitor));
    }
}

bool MonitorManager::hasMonitor(const std::string& name) const
{
    return getMonitor(name) != nullptr;
}

std::shared_ptr<Monitor> MonitorManager::getMonitor(const std::string& name) const
{
    for (const auto& monitor : monitors_)
    {
        if (monitor && monitor->getName() == name)
        {
            return monitor;
        }
    }
    return {};
}

VoidResult MonitorManager::sampleAll(double time) const
{
    for (const auto& monitor : monitors_)
    {
        if (monitor)
        {
            if (auto result = monitor->sample(time); !result)
            {
                return result;
            }
        }
    }
    return VoidResult::success();
}

std::unordered_map<std::string, double> MonitorManager::latestScalars() const
{
    std::unordered_map<std::string, double> scalars;
    for (const auto& monitor : monitors_)
    {
        if (!monitor)
        {
            continue;
        }

        if (const auto* sample = monitor->latestSample(); sample != nullptr)
        {
            for (const auto& [key, value] : sample->scalars)
            {
                scalars[monitor->getName() + "." + key] = value;
            }
        }
    }
    return scalars;
}

std::vector<std::string> MonitorManager::getMonitorNames() const
{
    std::vector<std::string> names;
    names.reserve(monitors_.size());
    for (const auto& monitor : monitors_)
    {
        if (monitor)
        {
            names.push_back(monitor->getName());
        }
    }
    return names;
}

} // namespace Diagnostics
} // namespace SCDAT
