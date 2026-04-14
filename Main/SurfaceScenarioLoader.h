#pragma once

#include "SurfaceChargingCases.h"

#include <filesystem>

namespace SCDAT
{
namespace MainEntry
{

class SurfaceScenarioLoader
{
  public:
    Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset
    loadFromJson(const std::filesystem::path& json_path,
                 std::filesystem::path& output_path) const;
};

} // namespace MainEntry
} // namespace SCDAT
