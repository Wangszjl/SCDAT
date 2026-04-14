#pragma once

#include "SurfaceChargingCases.h"

#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

class SurfaceScenarioCatalog
{
  public:
    std::vector<std::string> listMainlinePresetNames() const;
    std::vector<std::string> listReplayPresetNames() const;

    bool tryGetMainlinePreset(const std::string& name,
                              SurfaceChargingScenarioPreset& preset) const;
    bool tryGetReplayPreset(const std::string& name,
                            SurfaceChargingScenarioPreset& preset) const;
    bool tryGetPreset(const std::string& name,
                      SurfaceChargingScenarioPreset& preset) const;

    SurfaceChargingScenarioPreset makeDefaultPreset() const;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
