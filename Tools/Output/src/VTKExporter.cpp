#include "../include/VTKExporter.h"

#include <fstream>
#include <iomanip>

namespace SCDAT
{
namespace Output
{

VoidResult VTKExporter::exportDataSet(const std::filesystem::path& path,
                                      const ColumnarDataSet& data_set) const
{
    if (!data_set.isConsistent())
    {
        return VoidResult::failure("VTK export rejected inconsistent data set");
    }

    std::filesystem::create_directories(path.parent_path().empty() ? std::filesystem::path(".")
                                                                  : path.parent_path());

    std::ofstream output(path);
    if (!output.is_open())
    {
        return VoidResult::failure("Failed to open vtk output: " + path.string());
    }

    output << "# vtk DataFile Version 3.0\n";
    output << "SCDAT columnar export\n";
    output << "ASCII\n";
    output << "DATASET POLYDATA\n";
    output << "POINTS " << data_set.axis_values.size() << " float\n";
    output << std::scientific << std::setprecision(8);
    for (const double axis_value : data_set.axis_values)
    {
        output << "0 0 " << axis_value << '\n';
    }

    output << "VERTICES " << data_set.axis_values.size() << ' '
           << data_set.axis_values.size() * 2 << '\n';
    for (std::size_t i = 0; i < data_set.axis_values.size(); ++i)
    {
        output << "1 " << i << '\n';
    }

    output << "POINT_DATA " << data_set.axis_values.size() << '\n';
    for (const auto& [name, values] : data_set.scalar_series)
    {
        output << "SCALARS " << name << " float 1\n";
        output << "LOOKUP_TABLE default\n";
        for (const double value : values)
        {
            output << value << '\n';
        }
    }

    return VoidResult::success();
}

} // namespace Output
} // namespace SCDAT
