#pragma once

#include "../../Basic/include/VoidResult.h"
#include "SurfaceMaterialLoader.h"

#include <filesystem>

namespace SCDAT
{
namespace Material
{

class SurfaceMaterialImporter
{
  public:
    VoidResult importPath(const std::filesystem::path& path, MaterialDatabase& database) const;

  private:
    SurfaceMaterialLoader loader_;
};

} // namespace Material
} // namespace SCDAT

