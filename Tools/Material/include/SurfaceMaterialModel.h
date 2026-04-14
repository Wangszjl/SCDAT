#pragma once

#include <functional>
#include <memory>

#include "MaterialProperty.h"

namespace SCDAT
{
namespace Material
{

enum class SurfaceIncidentParticle
{
    Electron,
    Ion,
    Photon
};

struct SurfaceModelContext
{
    SurfaceIncidentParticle incident_particle = SurfaceIncidentParticle::Electron;
    double incident_energy_ev = 0.0;
    double incident_angle_rad = 0.0;
    double source_current_density_a_per_m2 = 0.0;
    double surface_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double counterelectrode_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double emission_temperature_ev = 2.0;
    double transport_length_m = 0.0;
    double exposed_area_m2 = 1.0;
};

struct SurfaceComponentCurrents
{
    double electron_collection_a_per_m2 = 0.0;
    double ion_collection_a_per_m2 = 0.0;
    double secondary_electron_emission_a_per_m2 = 0.0;
    double ion_secondary_electron_emission_a_per_m2 = 0.0;
    double backscatter_emission_a_per_m2 = 0.0;
    double photoelectron_emission_a_per_m2 = 0.0;
    double induced_conduction_a_per_m2 = 0.0;
    double net_current_a_per_m2 = 0.0;
};

enum class SurfaceCurrentRole
{
    Patch,
    Body
};

struct SurfaceRoleCurrentInputs
{
    SurfaceCurrentRole role = SurfaceCurrentRole::Patch;
    double electron_characteristic_energy_ev = 0.0;
    double ion_characteristic_energy_ev = 0.0;
    double patch_incidence_angle_deg = 0.0;
    double patch_flow_angle_deg = 0.0;
    double surface_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double counterelectrode_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double patch_photoelectron_temperature_ev = 2.0;
    double body_photoelectron_temperature_ev = 2.0;
    double patch_photo_current_density_a_per_m2 = 0.0;
    double body_photo_current_density_a_per_m2 = 0.0;
    double body_photo_emission_scale = 1.0;
    double dielectric_thickness_m = 0.0;
    double exposed_area_m2 = 1.0;
    double effective_conductivity_s_per_m = 0.0;
    double conductivity_scale = 1.0;
};

struct SurfaceRoleCurrentContexts
{
    SurfaceModelContext electron_context{};
    SurfaceModelContext ion_context{};
    SurfaceModelContext photo_context{};
};

struct SurfaceRoleCurrentBundle
{
    SurfaceRoleCurrentContexts contexts{};
    SurfaceComponentCurrents electron_currents{};
    SurfaceComponentCurrents ion_currents{};
    SurfaceComponentCurrents photo_currents{};
};

enum class LegacySecondaryYieldModel
{
    Whipple,
    Sims,
    Katz
};

struct LegacyEmissionIntegralBundle
{
    double secondary_integral = 0.0;
    double backscatter_integral = 0.0;
    double ion_secondary_integral = 0.0;
    double emission_escape_probability = 1.0;
};

struct LegacyEmissionResponseInputs
{
    LegacySecondaryYieldModel model = LegacySecondaryYieldModel::Whipple;
    double electron_temperature_ev = 0.0;
    double ion_temperature_ev = 0.0;
    double surface_potential_v = 0.0;
    double electron_collection_current_a_per_m2 = 0.0;
    double ion_collection_current_a_per_m2 = 0.0;
    double sims_exponent_n = 1.6;
};

struct LegacyEmissionResponse
{
    LegacyEmissionIntegralBundle integrals{};
    double secondary_emission_a_per_m2 = 0.0;
    double backscatter_emission_a_per_m2 = 0.0;
    double ion_secondary_emission_a_per_m2 = 0.0;
};

enum class LegacyCollectionParticle
{
    Electron,
    Ion
};

struct LegacyPopulationResponseInputs
{
    LegacyCollectionParticle particle = LegacyCollectionParticle::Electron;
    LegacySecondaryYieldModel model = LegacySecondaryYieldModel::Whipple;
    double density_m3 = 0.0;
    double temperature_ev = 0.0;
    double mass_amu = 1.0;
    double surface_potential_v = 0.0;
    double electron_characteristic_energy_ev = 1.0;
    double ion_characteristic_energy_ev = 1.0;
    double electron_collection_coefficient = 1.0;
    double ion_collection_coefficient = 1.0;
    double sims_exponent_n = 1.6;
};

struct LegacyPopulationResponse
{
    double collection_current_a_per_m2 = 0.0;
    LegacyEmissionResponse emission{};
};

class SurfaceYieldFunction
{
  public:
    virtual ~SurfaceYieldFunction() = default;
    virtual double evaluate(const MaterialProperty& material,
                            const SurfaceModelContext& context) const = 0;
};

class SurfacePhotoEmissionModel
{
  public:
    virtual ~SurfacePhotoEmissionModel() = default;
    virtual double evaluateCurrentDensity(const MaterialProperty& material,
                                          const SurfaceModelContext& context) const = 0;
};

class SurfaceInducedConductModel
{
  public:
    virtual ~SurfaceInducedConductModel() = default;
    virtual double evaluateCurrentDensity(const MaterialProperty& material,
                                          const SurfaceModelContext& context) const = 0;
    virtual double evaluateDidv(const MaterialProperty& material,
                                const SurfaceModelContext& context) const = 0;
};

class SurfaceMaterialModel
{
  public:
    virtual ~SurfaceMaterialModel() = default;

    virtual const char* modelFamily() const = 0;
    virtual SurfaceComponentCurrents evaluateCurrents(const MaterialProperty& material,
                                                      const SurfaceModelContext& context) const = 0;
};

class BasicSurfaceMaterialModel final : public SurfaceMaterialModel
{
  public:
    const char* modelFamily() const override;
    SurfaceComponentCurrents evaluateCurrents(const MaterialProperty& material,
                                              const SurfaceModelContext& context) const override;
};

enum class SurfaceMaterialModelVariant
{
        Basic,
        Erosion
};

class ErosionSurfaceMaterialModel final : public SurfaceMaterialModel
{
    public:
        const char* modelFamily() const override;
        SurfaceComponentCurrents evaluateCurrents(const MaterialProperty& material,
                                                                                            const SurfaceModelContext& context) const override;
};

SurfaceMaterialModelVariant resolveSurfaceMaterialModelVariant(
        const MaterialProperty& material);
const char* surfaceMaterialModelVariantFamily(SurfaceMaterialModelVariant variant);
const char* resolveSurfaceMaterialModelFamily(const MaterialProperty& material);

SurfaceRoleCurrentContexts makeSurfaceRoleCurrentContexts(
    const SurfaceRoleCurrentInputs& inputs);
SurfaceRoleCurrentBundle evaluateSurfaceRoleCurrentBundle(
    const MaterialProperty& material, const SurfaceRoleCurrentInputs& inputs,
    const SurfaceMaterialModel& model);
SurfaceRoleCurrentBundle evaluateSurfaceRoleCurrentBundle(
    const MaterialProperty& material, const SurfaceRoleCurrentInputs& inputs);
double evaluateLegacySecondaryElectronYield(const MaterialProperty& material,
                                            LegacySecondaryYieldModel model,
                                            double incident_energy_ev,
                                            double sims_exponent_n = 1.6);
double evaluateLegacyIonSecondaryElectronYield(const MaterialProperty& material,
                                               double incident_energy_ev);
double evaluateLegacyBackscatterYield(const MaterialProperty& material,
                                      double incident_energy_ev);
double evaluateLegacyRamPatchCurrentDensity(double surface_potential_v,
                                            double ion_density_m3,
                                            double ion_temperature_ev,
                                            double ion_mass_amu,
                                            double projected_flow_speed_m_per_s);
double evaluateLegacyRamBodyCurrentDensity(double surface_potential_v,
                                           double ion_density_m3,
                                           double ion_temperature_ev,
                                           double ion_mass_amu,
                                           double flow_speed_m_per_s);
double integrateMaxwellianYield(double temperature_ev,
                                const std::function<double(double)>& yield_evaluator);
double evaluateLegacyEmissionEscapeProbability(const MaterialProperty& material,
                                               double surface_potential_v);
double evaluateSurfaceFieldEmissionCurrentDensity(const MaterialProperty& material,
                                                  double surface_potential_v,
                                                  double reference_potential_v,
                                                  double normal_electric_field_v_per_m,
                                                  double effective_sheath_length_m);
double evaluateSurfaceThermionicEmissionCurrentDensity(const MaterialProperty& material,
                                                       double surface_temperature_k,
                                                       double surface_potential_v);
double evaluateSurfaceEffectiveConductivitySPerM(const MaterialProperty& material,
                                                 double base_conductivity_s_per_m,
                                                 double radiation_induced_conductivity_s_per_m,
                                                 double surface_potential_v,
                                                 double reference_potential_v,
                                                 double normal_electric_field_v_per_m,
                                                 double local_charge_density_c_per_m3,
                                                 double dielectric_thickness_m);
LegacyEmissionIntegralBundle evaluateLegacyEmissionIntegralBundle(
    const MaterialProperty& material, LegacySecondaryYieldModel model,
    double electron_temperature_ev, double ion_temperature_ev, double surface_potential_v,
    double sims_exponent_n = 1.6);
LegacyEmissionResponse evaluateLegacyEmissionResponse(
    const MaterialProperty& material, const LegacyEmissionResponseInputs& inputs);
LegacyPopulationResponse evaluateLegacyPopulationResponse(
    const MaterialProperty& material, const LegacyPopulationResponseInputs& inputs);

} // namespace Material
} // namespace SCDAT
