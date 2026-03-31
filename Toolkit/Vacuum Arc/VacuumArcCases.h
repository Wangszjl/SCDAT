#pragma once

#include "SurfaceDischargeArcAlgorithm.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

struct VacuumArcScenarioPreset
{
    std::string name;
    std::string description;
    DischargeConfiguration config;
    double time_step_s = 1.0e-9;
    std::size_t steps = 8;
    std::filesystem::path default_output_csv;
};

std::vector<std::string> listVacuumArcScenarioPresetNames();
bool tryGetVacuumArcScenarioPreset(const std::string& name, VacuumArcScenarioPreset& preset);
VacuumArcScenarioPreset makeDefaultVacuumArcScenarioPreset();

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
