#pragma once

#include <map>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Output
{

struct ColumnarDataSet
{
    std::string axis_name = "index";
    std::vector<double> axis_values;
    std::map<std::string, std::vector<double>> scalar_series;
    std::unordered_map<std::string, std::string> metadata;

    bool isConsistent() const
    {
        for (const auto& [_, values] : scalar_series)
        {
            if (values.size() != axis_values.size())
            {
                return false;
            }
        }
        return true;
    }

    bool hasOnlyFiniteValues() const
    {
        for (const double value : axis_values)
        {
            if (!std::isfinite(value))
            {
                return false;
            }
        }

        for (const auto& [_, values] : scalar_series)
        {
            for (const double value : values)
            {
                if (!std::isfinite(value))
                {
                    return false;
                }
            }
        }
        return true;
    }
};

struct ColumnSummary
{
    double minimum = 0.0;
    double maximum = 0.0;
    double mean = 0.0;
};

struct DataSetSummary
{
    std::size_t rows = 0;
    std::string axis_name;
    std::unordered_map<std::string, ColumnSummary> column_summaries;
};

} // namespace Output
} // namespace SCDAT
