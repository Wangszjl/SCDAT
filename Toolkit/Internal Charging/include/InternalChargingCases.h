#pragma once

#include "SpacecraftInternalChargingAlgorithm.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

struct InternalChargingScenarioPreset
{
    std::string name;
    std::string description;
    InternalChargingConfiguration config;
    double time_step_s = 1.0e-6;
    std::size_t steps = 6;
    std::filesystem::path default_output_csv;
};

std::vector<std::string> listInternalChargingScenarioPresetNames();
bool tryGetInternalChargingScenarioPreset(const std::string& name,
                                          InternalChargingScenarioPreset& preset);
InternalChargingScenarioPreset makeDefaultInternalChargingScenarioPreset();

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
