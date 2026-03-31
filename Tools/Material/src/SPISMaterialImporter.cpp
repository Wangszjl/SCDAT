#include "../include/SPISMaterialImporter.h"

namespace SCDAT
{
namespace Material
{

VoidResult SPISMaterialImporter::importPath(const std::filesystem::path& path,
                                            MaterialDatabase& database) const
{
    if (!std::filesystem::exists(path))
    {
        return VoidResult::failure("Material import path does not exist: " + path.string());
    }

    if (std::filesystem::is_directory(path))
    {
        for (const auto& entry : std::filesystem::directory_iterator(path))
        {
            if (!entry.is_regular_file())
            {
                continue;
            }

            const auto extension = entry.path().extension().string();
            if (extension == ".csv")
            {
                if (auto result = loader_.loadCsv(entry.path(), database); !result)
                {
                    return result;
                }
            }
            else if (extension == ".txt" || extension == ".cfg" || extension == ".dat")
            {
                if (auto result = loader_.loadKeyValueFile(entry.path(), database); !result)
                {
                    return result;
                }
            }
        }

        return VoidResult::success();
    }

    const auto extension = path.extension().string();
    if (extension == ".csv")
    {
        return loader_.loadCsv(path, database);
    }
    return loader_.loadKeyValueFile(path, database);
}

} // namespace Material
} // namespace SCDAT
