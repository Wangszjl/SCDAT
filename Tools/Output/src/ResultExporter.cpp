#include "../include/ResultExporter.h"

namespace SCDAT
{
namespace Output
{

VoidResult ResultExporter::exportDataSet(const std::filesystem::path& csv_path,
                                         const ColumnarDataSet& data_set,
                                         ExportArtifacts* artifacts) const
{
    if (!data_set.isConsistent())
    {
        return VoidResult::failure("Result export rejected inconsistent data set");
    }
    if (!data_set.hasOnlyFiniteValues())
    {
        return VoidResult::failure("Result export rejected non-finite data set");
    }

    std::filesystem::path stem_path = csv_path;
    stem_path.replace_extension("");
    const auto vtk_path = stem_path.string() + ".vtk";
    const auto hdf5_path = stem_path.string() + ".h5";

    if (auto result = csv_exporter_.exportDataSet(csv_path, data_set); !result)
    {
        return result;
    }
    if (auto result = vtk_exporter_.exportDataSet(vtk_path, data_set); !result)
    {
        return result;
    }
    if (auto result = hdf5_exporter_.exportDataSet(hdf5_path, data_set); !result)
    {
        return result;
    }

    if (artifacts != nullptr)
    {
        artifacts->csv_path = csv_path;
        artifacts->vtk_path = vtk_path;
        artifacts->hdf5_path = hdf5_path;
    }

    return VoidResult::success();
}

} // namespace Output
} // namespace SCDAT
