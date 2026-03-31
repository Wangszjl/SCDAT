#pragma once

#include "../../Basic/include/VoidResult.h"
#include "SPISMaterialLoader.h"

#include <filesystem>

namespace SCDAT
{
namespace Material
{

class SPISMaterialImporter
{
  public:
    VoidResult importPath(const std::filesystem::path& path, MaterialDatabase& database) const;

  private:
    SPISMaterialLoader loader_;
};

} // namespace Material
} // namespace SCDAT
