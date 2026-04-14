#include "SurfaceScenarioCatalog.h"

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

std::vector<std::string> SurfaceScenarioCatalog::listMainlinePresetNames() const
{
    return listSurfaceChargingScenarioPresetNames();
}

std::vector<std::string> SurfaceScenarioCatalog::listReplayPresetNames() const
{
    return listSurfaceChargingReplayScenarioPresetNames();
}

bool SurfaceScenarioCatalog::tryGetMainlinePreset(const std::string& name,
                                                  SurfaceChargingScenarioPreset& preset) const
{
    return tryGetSurfaceChargingMainlineScenarioPreset(name, preset);
}

bool SurfaceScenarioCatalog::tryGetReplayPreset(const std::string& name,
                                                SurfaceChargingScenarioPreset& preset) const
{
    return tryGetSurfaceChargingReplayScenarioPreset(name, preset);
}

bool SurfaceScenarioCatalog::tryGetPreset(const std::string& name,
                                          SurfaceChargingScenarioPreset& preset) const
{
    return tryGetSurfaceChargingScenarioPreset(name, preset);
}

SurfaceChargingScenarioPreset SurfaceScenarioCatalog::makeDefaultPreset() const
{
    return makeDefaultSurfaceChargingScenarioPreset();
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
