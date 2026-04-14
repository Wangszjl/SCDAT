#include "../include/SurfaceMaterialLoader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace SCDAT
{
namespace Material
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

std::string trim(std::string value)
{
    const auto is_space = [](unsigned char ch) { return std::isspace(ch) != 0; };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(),
                             [&](unsigned char ch) { return !is_space(ch); }));
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [&](unsigned char ch) { return !is_space(ch); })
                    .base(),
                value.end());
    return value;
}

Mesh::MaterialType parseMaterialType(const std::string& type_name)
{
    std::string normalized = type_name;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });

    if (normalized == "dielectric")
    {
        return Mesh::MaterialType::DIELECTRIC;
    }
    if (normalized == "conductor")
    {
        return Mesh::MaterialType::CONDUCTOR;
    }
    return Mesh::MaterialType::VACUUM;
}

void applyCommonField(const std::string& key, const std::string& value, MaterialProperty& material)
{
    if (key == "name")
    {
        material.setName(value);
        return;
    }
    if (key == "type")
    {
        material.setType(parseMaterialType(value));
        return;
    }

    const double numeric_value = std::stod(value);
    if (key == "id")
    {
        material.setId(static_cast<Mesh::MaterialId>(std::lround(numeric_value)));
    }
    else if (key == "permittivity")
    {
        material.setPermittivity(numeric_value);
    }
    else if (key == "conductivity")
    {
        material.setConductivity(numeric_value);
    }
    else if (key == "permeability")
    {
        material.setPermeability(numeric_value);
    }
    else if (key == "work_function_ev")
    {
        material.setWorkFunctionEv(numeric_value);
    }
    else if (key == "breakdown_field_v_per_m")
    {
        material.setBreakdownFieldVPerM(numeric_value);
    }
    else if (key == "secondary_electron_yield")
    {
        material.setSecondaryElectronYield(numeric_value);
    }
    else if (key == "mass_density_kg_per_m3")
    {
        material.setMassDensityKgPerM3(numeric_value);
    }
    else
    {
        material.setScalarProperty(key, numeric_value);
    }
}

} // namespace

VoidResult SurfaceMaterialLoader::loadCsv(const std::filesystem::path& path,
                                       MaterialDatabase& database) const
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        return VoidResult::failure("Failed to open material csv: " + path.string());
    }

    std::string header_line;
    if (!std::getline(input, header_line))
    {
        return VoidResult::failure("Material csv is empty: " + path.string());
    }

    const auto headers = splitCsvLine(header_line);
    if (headers.empty())
    {
        return VoidResult::failure("Material csv header is invalid: " + path.string());
    }

    std::string line;
    while (std::getline(input, line))
    {
        if (trim(line).empty())
        {
            continue;
        }

        const auto tokens = splitCsvLine(line);
        MaterialProperty material;
        for (std::size_t i = 0; i < std::min(headers.size(), tokens.size()); ++i)
        {
            const std::string key = trim(headers[i]);
            const std::string value = trim(tokens[i]);
            if (!value.empty())
            {
                applyCommonField(key, value, material);
            }
        }
        database.addMaterial(material);
    }

    return VoidResult::success();
}

VoidResult SurfaceMaterialLoader::loadKeyValueFile(const std::filesystem::path& path,
                                                MaterialDatabase& database) const
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        return VoidResult::failure("Failed to open material file: " + path.string());
    }

    MaterialProperty material;
    bool has_pending_material = false;
    std::string line;
    while (std::getline(input, line))
    {
        line = trim(line);
        if (line.empty())
        {
            if (has_pending_material)
            {
                database.addMaterial(material);
                material = MaterialProperty{};
                has_pending_material = false;
            }
            continue;
        }

        if (line.front() == '#' || line.front() == ';')
        {
            continue;
        }

        const auto delimiter = line.find('=');
        if (delimiter == std::string::npos)
        {
            return VoidResult::failure("Invalid material key-value line: " + line);
        }

        const std::string key = trim(line.substr(0, delimiter));
        const std::string value = trim(line.substr(delimiter + 1));
        applyCommonField(key, value, material);
        has_pending_material = true;
    }

    if (has_pending_material)
    {
        database.addMaterial(material);
    }

    return VoidResult::success();
}

} // namespace Material
} // namespace SCDAT

