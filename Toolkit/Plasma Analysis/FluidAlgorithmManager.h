#pragma once

#include "FluidAlgorithmAdapter.h"

#include <memory>
#include <string>
#include <unordered_map>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

class FluidAlgorithmManager
{
  public:
    void addRegion(const std::string& name, std::shared_ptr<FluidAlgorithmAdapter> adapter);
    bool advanceAll(double dt) const;
    std::shared_ptr<FluidAlgorithmAdapter> getRegion(const std::string& name) const;

  private:
    std::unordered_map<std::string, std::shared_ptr<FluidAlgorithmAdapter>> regions_;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
