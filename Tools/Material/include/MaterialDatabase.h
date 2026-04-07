#pragma once

#include "MaterialProperty.h"

#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Material
{

class MaterialDatabase
{
  public:
    MaterialDatabase();
    virtual ~MaterialDatabase() = default;

    void clear();
    void addMaterial(const MaterialProperty& property);
    void loadDefaultLibrary();

    bool hasMaterial(Mesh::MaterialId id) const;
    bool hasMaterial(const std::string& name) const;

    MaterialProperty* findById(Mesh::MaterialId id);
    const MaterialProperty* findById(Mesh::MaterialId id) const;
    MaterialProperty* findByName(const std::string& name);
    const MaterialProperty* findByName(const std::string& name) const;

    std::vector<MaterialProperty> getAllMaterials() const;

  protected:
    static std::string normalizeName(std::string_view name);

    std::map<Mesh::MaterialId, MaterialProperty> materials_by_id_;
    std::unordered_map<std::string, Mesh::MaterialId> name_to_id_;
};

class UnifiedMaterialDatabase : public MaterialDatabase
{
  public:
    UnifiedMaterialDatabase();

    void addAlias(const std::string& alias, Mesh::MaterialId id);
    const MaterialProperty* findByAliasOrName(const std::string& name) const;

  private:
    std::unordered_map<std::string, Mesh::MaterialId> aliases_;
};

} // namespace Material
} // namespace SCDAT
