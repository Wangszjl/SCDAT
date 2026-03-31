#pragma once

#include "../../Basic/include/VoidResult.h"
#include "ResultTypes.h"

#include <filesystem>

namespace SCDAT
{
namespace Output
{

class HDF5Exporter
{
  public:
    VoidResult exportDataSet(const std::filesystem::path& path,
                             const ColumnarDataSet& data_set) const;
};

} // namespace Output
} // namespace SCDAT
