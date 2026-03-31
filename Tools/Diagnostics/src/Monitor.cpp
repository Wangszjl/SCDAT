#include "../include/Monitor.h"

namespace SCDAT
{
namespace Diagnostics
{

Monitor::Monitor(std::string name) : name_(std::move(name)) {}

const std::string& Monitor::getName() const
{
    return name_;
}

const std::vector<MonitorSample>& Monitor::getSamples() const
{
    return samples_;
}

const MonitorSample* Monitor::latestSample() const
{
    if (samples_.empty())
    {
        return nullptr;
    }
    return &samples_.back();
}

void Monitor::recordSample(double time, std::unordered_map<std::string, double> scalars)
{
    samples_.push_back(MonitorSample{time, std::move(scalars)});
}

} // namespace Diagnostics
} // namespace SCDAT
