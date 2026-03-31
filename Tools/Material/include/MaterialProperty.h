#pragma once

#include "../../Mesh/include/MeshAlgorithms.h"

#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>

namespace SCDAT
{
namespace Material
{

class MaterialProperty
{
  public:
    using ScalarPropertyMap = std::unordered_map<std::string, double>;

    MaterialProperty(Mesh::MaterialId id = 0, Mesh::MaterialType type = Mesh::MaterialType::VACUUM,
                     std::string name = "vacuum")
        : id_(id), type_(type), name_(std::move(name))
    {
    }

    Mesh::MaterialId getId() const { return id_; }
    Mesh::MaterialType getType() const { return type_; }
    const std::string& getName() const { return name_; }

    void setId(Mesh::MaterialId id) { id_ = id; }
    void setType(Mesh::MaterialType type) { type_ = type; }
    void setName(const std::string& name) { name_ = name; }

    double getPermittivity() const { return permittivity_; }
    double getConductivity() const { return conductivity_; }
    double getPermeability() const { return permeability_; }
    double getWorkFunctionEv() const { return work_function_ev_; }
    double getBreakdownFieldVPerM() const { return breakdown_field_v_per_m_; }
    double getSecondaryElectronYield() const { return secondary_electron_yield_; }
    double getMassDensityKgPerM3() const { return mass_density_kg_per_m3_; }

    void setPermittivity(double value) { permittivity_ = value; }
    void setConductivity(double value) { conductivity_ = value; }
    void setPermeability(double value) { permeability_ = value; }
    void setWorkFunctionEv(double value) { work_function_ev_ = value; }
    void setBreakdownFieldVPerM(double value) { breakdown_field_v_per_m_ = value; }
    void setSecondaryElectronYield(double value) { secondary_electron_yield_ = value; }
    void setMassDensityKgPerM3(double value) { mass_density_kg_per_m3_ = value; }

    double getDielectricConstant() const { return permittivity_; }
    bool isConductor() const { return type_ == Mesh::MaterialType::CONDUCTOR; }
    bool isDielectric() const { return type_ == Mesh::MaterialType::DIELECTRIC; }
    bool isVacuum() const { return type_ == Mesh::MaterialType::VACUUM; }

    void setScalarProperty(const std::string& key, double value) { scalar_properties_[key] = value; }

    double getScalarProperty(const std::string& key, double default_value = 0.0) const
    {
        if (const auto it = scalar_properties_.find(key); it != scalar_properties_.end())
        {
            return it->second;
        }
        return default_value;
    }

    const ScalarPropertyMap& getScalarProperties() const { return scalar_properties_; }

    std::string toString() const
    {
        std::ostringstream stream;
        stream << "MaterialProperty{id=" << id_ << ", name=" << name_
               << ", eps=" << permittivity_ << ", sigma=" << conductivity_
               << ", mu=" << permeability_ << ", phi_w=" << work_function_ev_
               << ", Ebd=" << breakdown_field_v_per_m_
               << ", delta=" << secondary_electron_yield_
               << ", rho=" << mass_density_kg_per_m3_ << "}";
        return stream.str();
    }

  private:
    Mesh::MaterialId id_ = 0;
    Mesh::MaterialType type_ = Mesh::MaterialType::VACUUM;
    std::string name_;
    double permittivity_ = 1.0;
    double conductivity_ = 0.0;
    double permeability_ = 1.0;
    double work_function_ev_ = 0.0;
    double breakdown_field_v_per_m_ = 0.0;
    double secondary_electron_yield_ = 0.0;
    double mass_density_kg_per_m3_ = 0.0;
    ScalarPropertyMap scalar_properties_;
};

using MaterialPropertyPtr = std::shared_ptr<MaterialProperty>;

} // namespace Material
} // namespace SCDAT
