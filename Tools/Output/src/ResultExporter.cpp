#include "../include/ResultExporter.h"

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace SCDAT
{
namespace Output
{

namespace
{

std::string escapeJson(const std::string& text)
{
    std::string escaped;
    escaped.reserve(text.size() + 8);

    for (const unsigned char ch : text)
    {
        switch (ch)
        {
        case '"':
            escaped += "\\\"";
            break;
        case '\\':
            escaped += "\\\\";
            break;
        case '\b':
            escaped += "\\b";
            break;
        case '\f':
            escaped += "\\f";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            if (ch < 0x20)
            {
                std::ostringstream hex;
                hex << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<int>(ch);
                escaped += hex.str();
            }
            else
            {
                escaped.push_back(static_cast<char>(ch));
            }
            break;
        }
    }

    return escaped;
}

std::string toUtcIso8601()
{
    const auto now = std::chrono::system_clock::now();
    const auto now_time_t = std::chrono::system_clock::to_time_t(now);

    std::tm utc_tm{};
#if defined(_WIN32)
    gmtime_s(&utc_tm, &now_time_t);
#else
    gmtime_r(&now_time_t, &utc_tm);
#endif

    std::ostringstream output;
    output << std::put_time(&utc_tm, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

std::string joinSeriesNames(const ColumnarDataSet& data_set)
{
    std::ostringstream output;
    bool first = true;
    for (const auto& entry : data_set.scalar_series)
    {
        if (!first)
        {
            output << ';';
        }
        output << entry.first;
        first = false;
    }
    return output.str();
}

} // namespace

ColumnarDataSet ResultExporter::buildContractedDataSet(
    const ColumnarDataSet& data_set,
    const std::filesystem::path& csv_path,
    const std::filesystem::path& vtk_path,
    const std::filesystem::path& hdf5_path) const
{
    ColumnarDataSet contracted_data_set = data_set;
    auto& metadata = contracted_data_set.metadata;

    if (metadata.find("module") == metadata.end())
    {
        metadata["module"] = "unknown";
    }

    metadata["output_contract_schema"] = "scdat.output.contract.v1";
    metadata["output_contract_version"] = "v1";
    metadata["output_axis_name"] = contracted_data_set.axis_name;
    metadata["output_row_count"] = std::to_string(contracted_data_set.axis_values.size());
    metadata["output_series_count"] = std::to_string(contracted_data_set.scalar_series.size());
    metadata["output_series_names"] = joinSeriesNames(contracted_data_set);
    metadata["output_artifact_csv"] = std::filesystem::absolute(csv_path).string();
    metadata["output_artifact_vtk"] = std::filesystem::absolute(vtk_path).string();
    metadata["output_artifact_h5"] = std::filesystem::absolute(hdf5_path).string();
    metadata["output_export_utc"] = toUtcIso8601();

    return contracted_data_set;
}

VoidResult ResultExporter::writeMetadataSidecar(const std::filesystem::path& artifact_path,
                                                const char* artifact_format,
                                                const ColumnarDataSet& data_set) const
{
    const auto sidecar_path = std::filesystem::path(artifact_path.string() + ".metadata.json");
    const auto parent_path = sidecar_path.parent_path();
    if (!parent_path.empty())
    {
        std::filesystem::create_directories(parent_path);
    }

    std::ofstream output(sidecar_path);
    if (!output.is_open())
    {
        return VoidResult::failure("Failed to open metadata sidecar output: " +
                                   sidecar_path.string());
    }

    output << "{\n";
    output << "  \"output_contract_schema\": \""
           << escapeJson(data_set.metadata.at("output_contract_schema")) << "\",\n";
    output << "  \"output_contract_version\": \""
           << escapeJson(data_set.metadata.at("output_contract_version")) << "\",\n";
    output << "  \"artifact_format\": \"" << escapeJson(artifact_format) << "\",\n";
    output << "  \"artifact_path\": \""
           << escapeJson(std::filesystem::absolute(artifact_path).string()) << "\",\n";
    output << "  \"metadata\": {\n";

    std::size_t index = 0;
    for (const auto& entry : data_set.metadata)
    {
        output << "    \"" << escapeJson(entry.first) << "\": \""
               << escapeJson(entry.second) << "\"";
        if (index + 1 < data_set.metadata.size())
        {
            output << ',';
        }
        output << "\n";
        ++index;
    }

    output << "  }\n";
    output << "}\n";

    return VoidResult::success();
}

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
    const std::filesystem::path vtk_path = stem_path.string() + ".vtk";
    const std::filesystem::path hdf5_path = stem_path.string() + ".h5";
    const auto contracted_data_set = buildContractedDataSet(data_set, csv_path, vtk_path, hdf5_path);

    if (auto result = csv_exporter_.exportDataSet(csv_path, contracted_data_set); !result)
    {
        return result;
    }
    if (auto result = vtk_exporter_.exportDataSet(vtk_path, contracted_data_set); !result)
    {
        return result;
    }
    if (auto result = hdf5_exporter_.exportDataSet(hdf5_path, contracted_data_set); !result)
    {
        return result;
    }

    if (auto result = writeMetadataSidecar(csv_path, "csv", contracted_data_set); !result)
    {
        return result;
    }
    if (auto result = writeMetadataSidecar(vtk_path, "vtk", contracted_data_set); !result)
    {
        return result;
    }
    if (auto result = writeMetadataSidecar(hdf5_path, "h5", contracted_data_set); !result)
    {
        return result;
    }

    if (artifacts != nullptr)
    {
        artifacts->csv_path = csv_path;
        artifacts->vtk_path = vtk_path;
        artifacts->hdf5_path = hdf5_path;
        artifacts->csv_metadata_json_path = std::filesystem::path(csv_path.string() + ".metadata.json");
        artifacts->vtk_metadata_json_path = std::filesystem::path(vtk_path.string() + ".metadata.json");
        artifacts->hdf5_metadata_json_path = std::filesystem::path(hdf5_path.string() + ".metadata.json");
    }

    return VoidResult::success();
}

} // namespace Output
} // namespace SCDAT
