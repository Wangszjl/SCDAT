#include "../include/DataAnalyzer.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>

namespace SCDAT
{
namespace Output
{
namespace
{

std::vector<std::string> splitCsvLine(const std::string& line)
{
    std::vector<std::string> tokens;
    std::stringstream stream(line);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        tokens.push_back(token);
    }
    return tokens;
}

} // namespace

ColumnarDataSet DataAnalyzer::loadCsv(const std::filesystem::path& csv_path) const
{
    std::ifstream input(csv_path);
    if (!input.is_open())
    {
        throw std::runtime_error("Failed to open csv: " + csv_path.string());
    }

    std::string header_line;
    if (!std::getline(input, header_line))
    {
        throw std::runtime_error("CSV is empty: " + csv_path.string());
    }

    const auto headers = splitCsvLine(header_line);
    if (headers.empty())
    {
        throw std::runtime_error("CSV header missing in: " + csv_path.string());
    }

    ColumnarDataSet data_set;
    data_set.axis_name = headers.front();
    for (std::size_t i = 1; i < headers.size(); ++i)
    {
        data_set.scalar_series[headers[i]] = {};
    }

    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty())
        {
            continue;
        }

        const auto tokens = splitCsvLine(line);
        if (tokens.size() != headers.size())
        {
            throw std::runtime_error("Malformed csv row in: " + csv_path.string());
        }

        data_set.axis_values.push_back(std::stod(tokens.front()));
        for (std::size_t i = 1; i < headers.size(); ++i)
        {
            data_set.scalar_series[headers[i]].push_back(std::stod(tokens[i]));
        }
    }

    return data_set;
}

DataSetSummary DataAnalyzer::summarizeCsv(const std::filesystem::path& csv_path) const
{
    const auto data_set = loadCsv(csv_path);
    DataSetSummary summary;
    summary.rows = data_set.axis_values.size();
    summary.axis_name = data_set.axis_name;

    for (const auto& [name, values] : data_set.scalar_series)
    {
        if (values.empty())
        {
            summary.column_summaries[name] = {};
            continue;
        }

        ColumnSummary column_summary;
        column_summary.minimum = std::numeric_limits<double>::max();
        column_summary.maximum = -std::numeric_limits<double>::max();
        for (const double value : values)
        {
            column_summary.minimum = std::min(column_summary.minimum, value);
            column_summary.maximum = std::max(column_summary.maximum, value);
            column_summary.mean += value;
        }
        column_summary.mean /= static_cast<double>(values.size());
        summary.column_summaries[name] = column_summary;
    }

    return summary;
}

std::optional<std::filesystem::path>
DataAnalyzer::findLatestResultFile(const std::filesystem::path& results_dir,
                                   const std::string& extension) const
{
    if (!std::filesystem::exists(results_dir) || !std::filesystem::is_directory(results_dir))
    {
        return std::nullopt;
    }

    std::optional<std::filesystem::directory_entry> latest;
    for (const auto& entry : std::filesystem::directory_iterator(results_dir))
    {
        if (!entry.is_regular_file() || entry.path().extension() != extension)
        {
            continue;
        }

        if (!latest || std::filesystem::last_write_time(entry) >
                           std::filesystem::last_write_time(*latest))
        {
            latest = entry;
        }
    }

    if (!latest)
    {
        return std::nullopt;
    }
    return latest->path();
}

} // namespace Output
} // namespace SCDAT
