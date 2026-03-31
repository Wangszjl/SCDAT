#pragma once

#include "CouplingBase.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Coupling
{

class MultiPhysicsManager
{
  public:
    void addCoupling(std::shared_ptr<CouplingBase> coupling);
    bool hasCoupling(const std::string& name) const;
    std::shared_ptr<CouplingBase> getCoupling(const std::string& name) const;

    VoidResult executeAllCouplings(double dt) const;
    VoidResult executeIterativeCouplings(double dt) const;

    std::unordered_map<std::string, CouplingStatistics> getAllStatistics() const;

  private:
    std::vector<std::shared_ptr<CouplingBase>> couplings_;
};

} // namespace Coupling
} // namespace SCDAT
