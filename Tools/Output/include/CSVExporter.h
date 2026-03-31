#pragma once

#include "../../Basic/include/VoidResult.h"
#include "ResultTypes.h"

#include <filesystem>

namespace SCDAT
{
namespace Output
{

class CSVExporter
{
  public:
    VoidResult exportDataSet(const std::filesystem::path& path,
                             const ColumnarDataSet& data_set) const;
};

} // namespace Output
} // namespace SCDAT
