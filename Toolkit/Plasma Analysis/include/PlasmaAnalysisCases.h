#pragma once

#include "FluidAlgorithmConfig.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct PlasmaScenarioPreset
{
    std::string name;
    std::string description;
    FluidAlgorithmConfig config;
    std::size_t steps = 8;
    std::filesystem::path default_output_csv;
};

std::vector<std::string> listPlasmaScenarioPresetNames();
bool tryGetPlasmaScenarioPreset(const std::string& name, PlasmaScenarioPreset& preset);
PlasmaScenarioPreset makeDefaultPlasmaScenarioPreset();

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
