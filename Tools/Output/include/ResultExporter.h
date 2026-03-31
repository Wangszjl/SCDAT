#pragma once

#include "../../Basic/include/VoidResult.h"
#include "CSVExporter.h"
#include "HDF5Exporter.h"
#include "ResultTypes.h"
#include "VTKExporter.h"

#include <filesystem>

namespace SCDAT
{
namespace Output
{

struct ExportArtifacts
{
    std::filesystem::path csv_path;
    std::filesystem::path vtk_path;
    std::filesystem::path hdf5_path;
};

class ResultExporter
{
  public:
    VoidResult exportDataSet(const std::filesystem::path& csv_path, const ColumnarDataSet& data_set,
                             ExportArtifacts* artifacts = nullptr) const;

  private:
    CSVExporter csv_exporter_;
    VTKExporter vtk_exporter_;
    HDF5Exporter hdf5_exporter_;
};

} // namespace Output
} // namespace SCDAT
