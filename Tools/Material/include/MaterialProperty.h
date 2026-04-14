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

    struct SpisSurfaceParameters
    {
        double secondary_electron_yield_max = 0.0;
        double secondary_electron_peak_energy_ev = 300.0;
        double proton_secondary_electron_yield_max = 0.08;
        double proton_secondary_electron_peak_energy_ev = 350.0;
        double electron_backscatter_yield = 0.05;
        double photoelectron_yield = 0.0;
        double photoelectron_temperature_ev = 2.0;
        double surface_absorption_probability = 0.9;
        double surface_reflection_coefficient = 0.05;
        double surface_emission_scaling = 1.0;
        double surface_secondary_scale = 1.0;
        double surface_ion_secondary_scale = 1.0;
        double surface_backscatter_scale = 1.0;
        double surface_photo_emission_scale = 1.0;
        double surface_capacitance_f_per_m2 = 0.0;
        double surface_capacitance_scale = 1.0;
        double surface_conductivity_s_per_m = 0.0;
        double dark_conductivity_s_per_m = 0.0;
        double surface_conductivity_scale = 1.0;
        double relative_dielectric_constant = 1.0;
        double coating_thickness_m = 1.0e-6;
        bool conductor_semantics = false;
    };

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
    const SpisSurfaceParameters& getSpisSurfaceParameters() const { return spis_surface_parameters_; }

    void setPermittivity(double value) { permittivity_ = value; }
    void setConductivity(double value) { conductivity_ = value; }
    void setPermeability(double value) { permeability_ = value; }
    void setWorkFunctionEv(double value) { work_function_ev_ = value; }
    void setBreakdownFieldVPerM(double value) { breakdown_field_v_per_m_ = value; }
    void setSecondaryElectronYield(double value) { secondary_electron_yield_ = value; }
    void setMassDensityKgPerM3(double value) { mass_density_kg_per_m3_ = value; }
    void setSpisSurfaceParameters(const SpisSurfaceParameters& value)
    {
        spis_surface_parameters_ = value;
    }

    double getDielectricConstant() const { return permittivity_; }
    bool isConductor() const
    {
        return type_ == Mesh::MaterialType::CONDUCTOR || spis_surface_parameters_.conductor_semantics;
    }
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

    double getSecondaryElectronYieldMax() const
    {
        return scalarFallback("secondary_electron_yield_max",
                              spis_surface_parameters_.secondary_electron_yield_max,
                              secondary_electron_yield_);
    }

    void setSecondaryElectronYieldMax(double value)
    {
        spis_surface_parameters_.secondary_electron_yield_max = value;
        scalar_properties_["secondary_electron_yield_max"] = value;
    }

    double getSecondaryElectronPeakEnergyEv() const
    {
        return scalarFallback("secondary_electron_peak_energy_ev",
                              spis_surface_parameters_.secondary_electron_peak_energy_ev,
                              getScalarProperty("secondary_yield_peak_energy_ev", 300.0));
    }

    void setSecondaryElectronPeakEnergyEv(double value)
    {
        spis_surface_parameters_.secondary_electron_peak_energy_ev = value;
        scalar_properties_["secondary_electron_peak_energy_ev"] = value;
        scalar_properties_["secondary_yield_peak_energy_ev"] = value;
    }

    double getProtonSecondaryElectronYieldMax() const
    {
        return scalarFallback("proton_secondary_electron_yield_max",
                              spis_surface_parameters_.proton_secondary_electron_yield_max,
                              getScalarProperty("ion_secondary_yield", 0.08));
    }

    void setProtonSecondaryElectronYieldMax(double value)
    {
        spis_surface_parameters_.proton_secondary_electron_yield_max = value;
        scalar_properties_["proton_secondary_electron_yield_max"] = value;
        scalar_properties_["ion_secondary_yield"] = value;
    }

    double getProtonSecondaryElectronPeakEnergyEv() const
    {
        return scalarFallback("proton_secondary_electron_peak_energy_ev",
                              spis_surface_parameters_.proton_secondary_electron_peak_energy_ev,
                              getScalarProperty("ion_secondary_peak_energy_kev", 0.35) * 1.0e3);
    }

    void setProtonSecondaryElectronPeakEnergyEv(double value)
    {
        spis_surface_parameters_.proton_secondary_electron_peak_energy_ev = value;
        scalar_properties_["proton_secondary_electron_peak_energy_ev"] = value;
        scalar_properties_["ion_secondary_peak_energy_kev"] = value / 1.0e3;
    }

    double getElectronBackscatterYield() const
    {
        return scalarFallback("electron_backscatter_yield",
                              spis_surface_parameters_.electron_backscatter_yield,
                              getScalarProperty("surface_backscatter_scale", 0.05));
    }

    void setElectronBackscatterYield(double value)
    {
        spis_surface_parameters_.electron_backscatter_yield = value;
        scalar_properties_["electron_backscatter_yield"] = value;
    }

    double getPhotoelectronYield() const
    {
        return scalarFallback("photoelectron_yield", spis_surface_parameters_.photoelectron_yield,
                              getScalarProperty("photoelectron_yield", 0.0));
    }

    void setPhotoelectronYield(double value)
    {
        spis_surface_parameters_.photoelectron_yield = value;
        scalar_properties_["photoelectron_yield"] = value;
    }

    double getPhotoelectronTemperatureEv() const
    {
        return scalarFallback("photoelectron_temperature_ev",
                              spis_surface_parameters_.photoelectron_temperature_ev,
                              getScalarProperty("photoelectron_temperature_ev", 2.0));
    }

    void setPhotoelectronTemperatureEv(double value)
    {
        spis_surface_parameters_.photoelectron_temperature_ev = value;
        scalar_properties_["photoelectron_temperature_ev"] = value;
    }

    double getSurfaceAbsorptionProbability() const
    {
        return scalarFallback("surface_absorption_probability",
                              spis_surface_parameters_.surface_absorption_probability, 0.9);
    }

    void setSurfaceAbsorptionProbability(double value)
    {
        spis_surface_parameters_.surface_absorption_probability = value;
        scalar_properties_["surface_absorption_probability"] = value;
    }

    double getSurfaceReflectionCoefficient() const
    {
        return scalarFallback("surface_reflection_coefficient",
                              spis_surface_parameters_.surface_reflection_coefficient, 0.05);
    }

    void setSurfaceReflectionCoefficient(double value)
    {
        spis_surface_parameters_.surface_reflection_coefficient = value;
        scalar_properties_["surface_reflection_coefficient"] = value;
    }

    double getSurfaceEmissionScaling() const
    {
        return scalarFallback("surface_emission_scaling",
                              spis_surface_parameters_.surface_emission_scaling, 1.0);
    }

    void setSurfaceEmissionScaling(double value)
    {
        spis_surface_parameters_.surface_emission_scaling = value;
        scalar_properties_["surface_emission_scaling"] = value;
    }

    double getSurfaceSecondaryScale() const
    {
        return scalarFallback("surface_secondary_scale",
                              spis_surface_parameters_.surface_secondary_scale, 1.0);
    }

    void setSurfaceSecondaryScale(double value)
    {
        spis_surface_parameters_.surface_secondary_scale = value;
        scalar_properties_["surface_secondary_scale"] = value;
    }

    double getSurfaceIonSecondaryScale() const
    {
        return scalarFallback("surface_ion_secondary_scale",
                              spis_surface_parameters_.surface_ion_secondary_scale, 1.0);
    }

    void setSurfaceIonSecondaryScale(double value)
    {
        spis_surface_parameters_.surface_ion_secondary_scale = value;
        scalar_properties_["surface_ion_secondary_scale"] = value;
    }

    double getSurfaceBackscatterScale() const
    {
        return scalarFallback("surface_backscatter_scale",
                              spis_surface_parameters_.surface_backscatter_scale, 1.0);
    }

    void setSurfaceBackscatterScale(double value)
    {
        spis_surface_parameters_.surface_backscatter_scale = value;
        scalar_properties_["surface_backscatter_scale"] = value;
    }

    double getSurfacePhotoEmissionScale() const
    {
        return scalarFallback("surface_photo_emission_scale",
                              spis_surface_parameters_.surface_photo_emission_scale, 1.0);
    }

    void setSurfacePhotoEmissionScale(double value)
    {
        spis_surface_parameters_.surface_photo_emission_scale = value;
        scalar_properties_["surface_photo_emission_scale"] = value;
    }

    double getSurfaceCapacitanceFPerM2() const
    {
        return scalarFallback("surface_capacitance_f_per_m2",
                              spis_surface_parameters_.surface_capacitance_f_per_m2,
                              getScalarProperty("surface_capacitance_f_per_m2", 0.0));
    }

    void setSurfaceCapacitanceFPerM2(double value)
    {
        spis_surface_parameters_.surface_capacitance_f_per_m2 = value;
        scalar_properties_["surface_capacitance_f_per_m2"] = value;
    }

    double getSurfaceCapacitanceScale() const
    {
        return scalarFallback("surface_capacitance_scale",
                              spis_surface_parameters_.surface_capacitance_scale, 1.0);
    }

    void setSurfaceCapacitanceScale(double value)
    {
        spis_surface_parameters_.surface_capacitance_scale = value;
        scalar_properties_["surface_capacitance_scale"] = value;
    }

    double deriveSurfaceCapacitanceScaleFactor() const
    {
        return std::clamp(getSurfaceCapacitanceScale() *
                              std::max(1.0, getSurfaceCapacitanceFPerM2() > 0.0 ? 1.0 : 0.0),
                          0.1, 10.0);
    }

    double getSurfaceConductivitySPerM() const
    {
        return scalarFallback("surface_conductivity_s_per_m",
                              spis_surface_parameters_.surface_conductivity_s_per_m,
                              getScalarProperty("surface_conductivity_s_per_m", conductivity_));
    }

    void setSurfaceConductivitySPerM(double value)
    {
        spis_surface_parameters_.surface_conductivity_s_per_m = value;
        scalar_properties_["surface_conductivity_s_per_m"] = value;
    }

    double getSurfaceConductivityScale() const
    {
        return scalarFallback("surface_conductivity_scale",
                              spis_surface_parameters_.surface_conductivity_scale, 1.0);
    }

    void setSurfaceConductivityScale(double value)
    {
        spis_surface_parameters_.surface_conductivity_scale = value;
        scalar_properties_["surface_conductivity_scale"] = value;
    }

    double deriveSurfaceConductivityScaleFactor() const
    {
        return std::clamp(
            getSurfaceConductivityScale() * std::max(0.0, getSurfaceConductivitySPerM()) /
                std::max(1.0e-18, std::max(getConductivity(), getDarkConductivitySPerM())),
            0.0, 1.0e12);
    }

    double getDarkConductivitySPerM() const
    {
        return scalarFallback("dark_conductivity_s_per_m",
                              spis_surface_parameters_.dark_conductivity_s_per_m,
                              getScalarProperty("dark_conductivity_s_per_m", conductivity_));
    }

    void setDarkConductivitySPerM(double value)
    {
        spis_surface_parameters_.dark_conductivity_s_per_m = value;
        scalar_properties_["dark_conductivity_s_per_m"] = value;
    }

    double getRelativeDielectricConstant() const
    {
        return scalarFallback("relative_dielectric_constant",
                              spis_surface_parameters_.relative_dielectric_constant,
                              permittivity_);
    }

    void setRelativeDielectricConstant(double value)
    {
        spis_surface_parameters_.relative_dielectric_constant = value;
        permittivity_ = value;
        scalar_properties_["relative_dielectric_constant"] = value;
    }

    double getCoatingThicknessM() const
    {
        return scalarFallback("coating_thickness_m", spis_surface_parameters_.coating_thickness_m,
                              getScalarProperty("coating_thickness_m", 1.0e-6));
    }

    void setCoatingThicknessM(double value)
    {
        spis_surface_parameters_.coating_thickness_m = value;
        scalar_properties_["coating_thickness_m"] = value;
    }

    void setConductorSemantics(bool value)
    {
        spis_surface_parameters_.conductor_semantics = value;
    }

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
    double scalarFallback(const std::string& key, double typed_value, double fallback_value) const
    {
        if (const auto it = scalar_properties_.find(key); it != scalar_properties_.end())
        {
            return it->second;
        }
        return typed_value != 0.0 ? typed_value : fallback_value;
    }

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
    SpisSurfaceParameters spis_surface_parameters_{};
    ScalarPropertyMap scalar_properties_;
};

using MaterialPropertyPtr = std::shared_ptr<MaterialProperty>;

} // namespace Material
} // namespace SCDAT
