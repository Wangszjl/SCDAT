#include "../include/HDF5Exporter.h"

#include <fstream>
#include <iomanip>

namespace SCDAT
{
namespace Output
{

VoidResult HDF5Exporter::exportDataSet(const std::filesystem::path& path,
                                       const ColumnarDataSet& data_set) const
{
    if (!data_set.isConsistent())
    {
        return VoidResult::failure("HDF5 export rejected inconsistent data set");
    }

    std::filesystem::create_directories(path.parent_path().empty() ? std::filesystem::path(".")
                                                                  : path.parent_path());

    std::ofstream output(path);
    if (!output.is_open())
    {
        return VoidResult::failure("Failed to open HDF5-style output: " + path.string());
    }

    output << "# SCDAT HDF5-like archive\n";
    output << "[metadata]\n";
    output << "axis_name=" << data_set.axis_name << '\n';
    for (const auto& [key, value] : data_set.metadata)
    {
        output << key << '=' << value << '\n';
    }

    output << "[dataset:" << data_set.axis_name << "]\n";
    output << std::scientific << std::setprecision(8);
    for (const double axis_value : data_set.axis_values)
    {
        output << axis_value << '\n';
    }

    for (const auto& [name, values] : data_set.scalar_series)
    {
        output << "[dataset:" << name << "]\n";
        for (const double value : values)
        {
            output << value << '\n';
        }
    }

    return VoidResult::success();
}

} // namespace Output
} // namespace SCDAT
