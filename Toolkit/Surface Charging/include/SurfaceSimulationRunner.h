#pragma once

#include "SurfaceChargingCases.h"

#include <filesystem>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceSimulationRunResult
{
    bool success = false;
    std::string error_message;
};

class SurfaceSimulationRunner
{
  public:
    SurfaceSimulationRunResult run(const SurfaceChargingScenarioPreset& preset,
                                   const std::filesystem::path& output_path) const;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
