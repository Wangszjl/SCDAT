#include "../include/MaterialDatabase.h"

#include <algorithm>
#include <cctype>

namespace SCDAT
{
namespace Material
{

MaterialDatabase::MaterialDatabase()
{
    loadDefaultLibrary();
}

void MaterialDatabase::clear()
{
    materials_by_id_.clear();
    name_to_id_.clear();
}

void MaterialDatabase::addMaterial(const MaterialProperty& property)
{
    materials_by_id_[property.getId()] = property;
    name_to_id_[normalizeName(property.getName())] = property.getId();
}

void MaterialDatabase::loadDefaultLibrary()
{
    clear();

    MaterialProperty vacuum(0, Mesh::MaterialType::VACUUM, "vacuum");
    vacuum.setPermittivity(1.0);
    addMaterial(vacuum);

    MaterialProperty aluminum(1, Mesh::MaterialType::CONDUCTOR, "aluminum");
    aluminum.setConductivity(3.5e7);
    aluminum.setWorkFunctionEv(4.08);
    aluminum.setMassDensityKgPerM3(2700.0);
    addMaterial(aluminum);

    MaterialProperty kapton(2, Mesh::MaterialType::DIELECTRIC, "kapton");
    kapton.setPermittivity(3.4);
    kapton.setConductivity(1.0e-15);
    kapton.setBreakdownFieldVPerM(2.5e8);
    kapton.setSecondaryElectronYield(1.2);
    kapton.setMassDensityKgPerM3(1420.0);
    addMaterial(kapton);

    MaterialProperty ptfe(3, Mesh::MaterialType::DIELECTRIC, "ptfe");
    ptfe.setPermittivity(2.1);
    ptfe.setConductivity(1.0e-18);
    ptfe.setBreakdownFieldVPerM(6.0e7);
    ptfe.setSecondaryElectronYield(1.4);
    ptfe.setMassDensityKgPerM3(2200.0);
    addMaterial(ptfe);
}

bool MaterialDatabase::hasMaterial(Mesh::MaterialId id) const
{
    return materials_by_id_.contains(id);
}

bool MaterialDatabase::hasMaterial(const std::string& name) const
{
    return name_to_id_.contains(normalizeName(name));
}

MaterialProperty* MaterialDatabase::findById(Mesh::MaterialId id)
{
    if (const auto it = materials_by_id_.find(id); it != materials_by_id_.end())
    {
        return &it->second;
    }
    return nullptr;
}

const MaterialProperty* MaterialDatabase::findById(Mesh::MaterialId id) const
{
    if (const auto it = materials_by_id_.find(id); it != materials_by_id_.end())
    {
        return &it->second;
    }
    return nullptr;
}

MaterialProperty* MaterialDatabase::findByName(const std::string& name)
{
    if (const auto it = name_to_id_.find(normalizeName(name)); it != name_to_id_.end())
    {
        return findById(it->second);
    }
    return nullptr;
}

const MaterialProperty* MaterialDatabase::findByName(const std::string& name) const
{
    if (const auto it = name_to_id_.find(normalizeName(name)); it != name_to_id_.end())
    {
        return findById(it->second);
    }
    return nullptr;
}

std::vector<MaterialProperty> MaterialDatabase::getAllMaterials() const
{
    std::vector<MaterialProperty> materials;
    materials.reserve(materials_by_id_.size());
    for (const auto& [_, property] : materials_by_id_)
    {
        materials.push_back(property);
    }
    return materials;
}

std::string MaterialDatabase::normalizeName(std::string_view name)
{
    std::string normalized;
    normalized.reserve(name.size());
    for (const char ch : name)
    {
        if (!std::isspace(static_cast<unsigned char>(ch)))
        {
            normalized.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
        }
    }
    return normalized;
}

void UnifiedMaterialDatabase::addAlias(const std::string& alias, Mesh::MaterialId id)
{
    aliases_[normalizeName(alias)] = id;
}

const MaterialProperty* UnifiedMaterialDatabase::findByAliasOrName(const std::string& name) const
{
    if (const auto it = aliases_.find(normalizeName(name)); it != aliases_.end())
    {
        return findById(it->second);
    }
    return findByName(name);
}

} // namespace Material
} // namespace SCDAT
