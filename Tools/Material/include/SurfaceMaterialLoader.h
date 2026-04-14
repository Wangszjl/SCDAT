#pragma once

#include "../../Basic/include/VoidResult.h"
#include "MaterialDatabase.h"

#include <filesystem>

namespace SCDAT
{
namespace Material
{

class SurfaceMaterialLoader
{
  public:
    VoidResult loadCsv(const std::filesystem::path& path, MaterialDatabase& database) const;
    VoidResult loadKeyValueFile(const std::filesystem::path& path, MaterialDatabase& database) const;
};

} // namespace Material
} // namespace SCDAT

