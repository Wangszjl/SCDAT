#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceChargingScenarioPreset
{
    std::string name;
    std::string description;
    SurfaceChargingConfig config;
    double time_step_s = 1.0e-9;
    std::size_t steps = 10;
    bool adaptive_time_stepping = false;
    double total_duration_s = 0.0;
    double minimum_time_step_s = 0.0;
    double maximum_time_step_s = 0.0;
    std::filesystem::path default_output_csv;
};

std::vector<std::string> listSurfaceChargingScenarioPresetNames();
bool tryGetSurfaceChargingScenarioPreset(const std::string& name,
                                         SurfaceChargingScenarioPreset& preset);
SurfaceChargingScenarioPreset makeDefaultSurfaceChargingScenarioPreset();

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
