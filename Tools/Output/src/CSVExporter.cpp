#include "../include/CSVExporter.h"

#include <fstream>
#include <iomanip>

namespace SCDAT
{
namespace Output
{

VoidResult CSVExporter::exportDataSet(const std::filesystem::path& path,
                                      const ColumnarDataSet& data_set) const
{
    if (!data_set.isConsistent())
    {
        return VoidResult::failure("CSV export rejected inconsistent data set");
    }

    std::filesystem::create_directories(path.parent_path().empty() ? std::filesystem::path(".")
                                                                  : path.parent_path());

    std::ofstream output(path);
    if (!output.is_open())
    {
        return VoidResult::failure("Failed to open csv output: " + path.string());
    }

    output << data_set.axis_name;
    for (const auto& [name, _] : data_set.scalar_series)
    {
        output << ',' << name;
    }
    output << '\n';

    output << std::scientific << std::setprecision(8);
    for (std::size_t row = 0; row < data_set.axis_values.size(); ++row)
    {
        output << data_set.axis_values[row];
        for (const auto& [_, values] : data_set.scalar_series)
        {
            output << ',' << values[row];
        }
        output << '\n';
    }

    return VoidResult::success();
}

} // namespace Output
} // namespace SCDAT
