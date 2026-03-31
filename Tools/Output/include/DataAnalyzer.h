#pragma once

#include "ResultTypes.h"

#include <filesystem>
#include <optional>
#include <string>

namespace SCDAT
{
namespace Output
{

class DataAnalyzer
{
  public:
    ColumnarDataSet loadCsv(const std::filesystem::path& csv_path) const;
    DataSetSummary summarizeCsv(const std::filesystem::path& csv_path) const;

    std::optional<std::filesystem::path>
    findLatestResultFile(const std::filesystem::path& results_dir, const std::string& extension) const;
};

} // namespace Output
} // namespace SCDAT
