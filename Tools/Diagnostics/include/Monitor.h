#pragma once

#include "../../Basic/include/VoidResult.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Diagnostics
{

struct MonitorSample
{
    double time = 0.0;
    std::unordered_map<std::string, double> scalars;
};

class Monitor
{
  public:
    explicit Monitor(std::string name);
    virtual ~Monitor() = default;

    const std::string& getName() const;
    virtual VoidResult sample(double time) = 0;

    const std::vector<MonitorSample>& getSamples() const;
    const MonitorSample* latestSample() const;

  protected:
    void recordSample(double time, std::unordered_map<std::string, double> scalars);

  private:
    std::string name_;
    std::vector<MonitorSample> samples_;
};

} // namespace Diagnostics
} // namespace SCDAT
