#pragma once

#include "RadiationDoseAlgorithm.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace Radiation
{

struct RadiationScenarioPreset
{
    std::string name;
    std::string description;
    RadiationConfiguration config;
    double time_step_s = 1.0e-3;
    std::size_t steps = 10;
    std::filesystem::path default_output_csv;
};

std::vector<std::string> listRadiationScenarioPresetNames();
bool tryGetRadiationScenarioPreset(const std::string& name, RadiationScenarioPreset& preset);
RadiationScenarioPreset makeDefaultRadiationScenarioPreset();

} // namespace Radiation
} // namespace Toolkit
} // namespace SCDAT
