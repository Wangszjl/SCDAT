#pragma once

#include "../../FieldSolver/include/SurfaceBarrierModels.h"
#include "MaterialProperty.h"
#include "SurfaceMaterialModel.h"

namespace SCDAT
{
namespace Material
{

struct SurfaceImpactState
{
    double incident_energy_ev = 0.0;
    double incident_angle_rad = 0.0;
    double particle_charge_coulomb = 0.0;
    double surface_temperature_k = 300.0;
    double surface_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double emission_temperature_ev = 0.0;
    double source_current_density_a_per_m2 = 0.0;
    double counterelectrode_potential_v = 0.0;
    double transport_length_m = 0.0;
    double exposed_area_m2 = 1.0;
};

struct SurfaceInteractionResult
{
    double absorbed_charge_coulomb = 0.0;
    double reflected_energy_ev = 0.0;
    double secondary_emitted_electrons = 0.0;
    double ion_secondary_emitted_electrons = 0.0;
    double backscattered_electrons = 0.0;
    double photoelectron_current_a_per_m2 = 0.0;
    double induced_conduction_current_a_per_m2 = 0.0;
    double deposited_heat_j_per_m2 = 0.0;
    SurfaceComponentCurrents component_currents{};
    std::string material_model_family = "legacy_surface_interaction_v1";
};

struct SurfaceInteractionBundleRequest
{
    SurfaceImpactState electron_impact{};
    SurfaceImpactState ion_impact{};
    FieldSolver::SurfaceBarrierState barrier_state{};
    double photo_incidence_scale = 1.0;
};

struct SurfaceInteractionEnvironment
{
    double electron_incident_energy_ev = 0.0;
    double electron_incident_angle_rad = 0.0;
    double ion_incident_energy_ev = 0.0;
    double ion_incident_angle_rad = 0.0;
    double surface_temperature_k = 300.0;
    double surface_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double emission_temperature_ev = 0.0;
    double source_current_density_a_per_m2 = 0.0;
    double counterelectrode_potential_v = 0.0;
    double transport_length_m = 0.0;
    double exposed_area_m2 = 1.0;
    double barrier_potential_v = 0.0;
    double photo_incidence_scale = 1.0;
};

enum class SurfaceInteractionRole
{
    Patch,
    Body
};

struct SurfaceInteractionRoleInputs
{
    SurfaceInteractionRole role = SurfaceInteractionRole::Patch;
    double electron_incident_energy_ev = 0.0;
    double electron_incident_angle_deg = 0.0;
    double ion_incident_energy_ev = 0.0;
    double ion_incident_angle_deg = 0.0;
    double surface_temperature_k = 300.0;
    double surface_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double counterelectrode_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double patch_photoelectron_temperature_ev = 0.0;
    double body_photoelectron_temperature_ev = 0.0;
    double patch_photo_current_density_a_per_m2 = 0.0;
    double body_photo_current_density_a_per_m2 = 0.0;
    double body_photo_emission_scale = 1.0;
    double dielectric_thickness_m = 0.0;
    double exposed_area_m2 = 1.0;
    double photo_incidence_scale = 1.0;
};

struct SurfaceInteractionBundleResult
{
    SurfaceInteractionResult electron_response{};
    SurfaceInteractionResult ion_response{};
    FieldSolver::SurfaceEmissionBarrierOutputs barrier_outputs{};
    bool valid = false;
};

struct SurfaceInteractionMaterialSelection
{
    double absorption_probability = 0.9;
    double reflection_coefficient = 0.05;
    double emission_scaling = 1.0;
    std::string interaction_family = "spis_material_indexed_surface_interaction_v1";
};

class SurfaceInteraction
{
  public:
    SurfaceInteraction();

    void setAbsorptionProbability(double value) { absorption_probability_ = value; }
    void setReflectionCoefficient(double value) { reflection_coefficient_ = value; }
    void setEmissionScaling(double value) { emission_scaling_ = value; }
    void setMaterialModel(const SurfaceMaterialModel* model) { material_model_ = model; }
    const char* materialModelFamily() const;

    SurfaceInteractionResult evaluate(const MaterialProperty& material,
                                      const SurfaceImpactState& impact) const;
    SurfaceInteractionBundleRequest makeBundleRequest(
        const SurfaceInteractionEnvironment& environment) const;
    static SurfaceInteractionMaterialSelection selectMaterialInteraction(
        const MaterialProperty& material);
SurfaceInteractionBundleResult evaluateBundle(const MaterialProperty& material,
                                                  const SurfaceInteractionBundleRequest& request) const;
    SurfaceInteractionBundleResult evaluateBundle(const MaterialProperty& material,
                                                  const SurfaceInteractionEnvironment& environment) const;
    SurfaceInteractionBundleResult evaluateMaterialIndexedBundle(
        const MaterialProperty& material, const SurfaceInteractionEnvironment& environment) const;
    static SurfaceInteractionEnvironment makeInteractionEnvironment(
        const SurfaceInteractionRoleInputs& inputs);

  private:
    double absorption_probability_ = 0.9;
    double reflection_coefficient_ = 0.05;
    double emission_scaling_ = 1.0;
    BasicSurfaceMaterialModel default_material_model_{};
    const SurfaceMaterialModel* material_model_ = nullptr;
};

} // namespace Material
} // namespace SCDAT
