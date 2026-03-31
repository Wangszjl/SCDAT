#pragma once

#include "Monitor.h"

#include "../../Basic/include/VoidResult.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Diagnostics
{

class MonitorManager
{
  public:
    void addMonitor(std::shared_ptr<Monitor> monitor);
    bool hasMonitor(const std::string& name) const;
    std::shared_ptr<Monitor> getMonitor(const std::string& name) const;
    VoidResult sampleAll(double time) const;

    std::unordered_map<std::string, double> latestScalars() const;
    std::vector<std::string> getMonitorNames() const;

  private:
    std::vector<std::shared_ptr<Monitor>> monitors_;
};

} // namespace Diagnostics
} // namespace SCDAT
