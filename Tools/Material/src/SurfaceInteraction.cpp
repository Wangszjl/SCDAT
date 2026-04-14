#include "../include/SurfaceInteraction.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Material
{

SurfaceInteraction::SurfaceInteraction() : material_model_(&default_material_model_) {}

const char* SurfaceInteraction::materialModelFamily() const
{
    return material_model_ == nullptr ? "legacy_surface_interaction_v1"
                                      : material_model_->modelFamily();
}

SurfaceInteractionResult SurfaceInteraction::evaluate(const MaterialProperty& material,
                                                      const SurfaceImpactState& impact) const
{
    const double clamped_angle = std::clamp(impact.incident_angle_rad, 0.0, 1.5707963267948966);
    const double cosine_factor = std::max(0.05, std::cos(clamped_angle));
    const double reflection =
        std::clamp(reflection_coefficient_ * cosine_factor, 0.0, 0.95);
    const double absorption = std::clamp(absorption_probability_ * (1.0 - 0.5 * reflection), 0.0,
                                         1.0);

    SurfaceModelContext context;
    context.incident_particle = impact.particle_charge_coulomb < 0.0
                                    ? SurfaceIncidentParticle::Electron
                                    : SurfaceIncidentParticle::Ion;
    context.incident_energy_ev = impact.incident_energy_ev;
    context.incident_angle_rad = impact.incident_angle_rad;
    context.source_current_density_a_per_m2 = impact.source_current_density_a_per_m2;
    context.surface_potential_v = impact.surface_potential_v;
    context.reference_potential_v = impact.reference_potential_v;
    context.counterelectrode_potential_v = impact.counterelectrode_potential_v;
    context.normal_electric_field_v_per_m = impact.normal_electric_field_v_per_m;
    context.emission_temperature_ev =
        impact.emission_temperature_ev > 0.0 ? impact.emission_temperature_ev
                                             : material.getPhotoelectronTemperatureEv();
    context.transport_length_m = impact.transport_length_m;
    context.exposed_area_m2 = std::max(1.0e-12, impact.exposed_area_m2);
    const auto component_currents =
        material_model_ == nullptr ? SurfaceComponentCurrents{}
                                   : material_model_->evaluateCurrents(material, context);
    const double incident_collection_a_per_m2 =
        std::max(component_currents.electron_collection_a_per_m2,
                 component_currents.ion_collection_a_per_m2);
    const double emission_scale = std::max(0.0, emission_scaling_);
    const double scaled_secondary_current_a_per_m2 =
        component_currents.secondary_electron_emission_a_per_m2 * emission_scale;
    const double scaled_ion_secondary_current_a_per_m2 =
        component_currents.ion_secondary_electron_emission_a_per_m2 * emission_scale;
    const double scaled_backscatter_current_a_per_m2 =
        component_currents.backscatter_emission_a_per_m2 * emission_scale;

    SurfaceInteractionResult result;
    result.absorbed_charge_coulomb = impact.particle_charge_coulomb * absorption;
    result.reflected_energy_ev = impact.incident_energy_ev * reflection;
    result.secondary_emitted_electrons =
        incident_collection_a_per_m2 > 0.0
            ? std::max(0.0, scaled_secondary_current_a_per_m2 / incident_collection_a_per_m2)
            : 0.0;
    result.ion_secondary_emitted_electrons =
        incident_collection_a_per_m2 > 0.0
            ? std::max(0.0, scaled_ion_secondary_current_a_per_m2 / incident_collection_a_per_m2)
            : 0.0;
    result.backscattered_electrons =
        incident_collection_a_per_m2 > 0.0
            ? std::max(0.0, scaled_backscatter_current_a_per_m2 / incident_collection_a_per_m2)
            : 0.0;
    result.photoelectron_current_a_per_m2 = component_currents.photoelectron_emission_a_per_m2;
    result.induced_conduction_current_a_per_m2 = component_currents.induced_conduction_a_per_m2;
    result.component_currents = component_currents;
    result.component_currents.secondary_electron_emission_a_per_m2 = scaled_secondary_current_a_per_m2;
    result.component_currents.ion_secondary_electron_emission_a_per_m2 =
        scaled_ion_secondary_current_a_per_m2;
    result.component_currents.backscatter_emission_a_per_m2 = scaled_backscatter_current_a_per_m2;
    result.material_model_family = materialModelFamily();
    result.deposited_heat_j_per_m2 = impact.incident_energy_ev * 1.602176634e-19 *
                                     (absorption + 0.25 * result.secondary_emitted_electrons);
    return result;
}

SurfaceInteractionBundleResult SurfaceInteraction::evaluateBundle(
    const MaterialProperty& material, const SurfaceInteractionBundleRequest& request) const
{
    SurfaceInteractionBundleResult result;
    result.electron_response = evaluate(material, request.electron_impact);
    result.ion_response = evaluate(material, request.ion_impact);

    FieldSolver::SurfaceEmissionBarrierComponentInputs components;
    components.electron_collection_a_per_m2 =
        result.electron_response.component_currents.electron_collection_a_per_m2;
    components.ion_collection_a_per_m2 =
        result.ion_response.component_currents.ion_collection_a_per_m2;
    components.secondary_emission_a_per_m2 =
        result.electron_response.component_currents.secondary_electron_emission_a_per_m2;
    components.ion_secondary_emission_a_per_m2 =
        result.ion_response.component_currents.ion_secondary_electron_emission_a_per_m2;
    components.backscatter_emission_a_per_m2 =
        result.electron_response.component_currents.backscatter_emission_a_per_m2;
    components.photo_emission_a_per_m2 =
        result.electron_response.photoelectron_current_a_per_m2;
    components.photo_incidence_scale = request.photo_incidence_scale;
    components.user_secondary_scale = material.getSurfaceSecondaryScale();
    components.user_ion_secondary_scale = material.getSurfaceIonSecondaryScale();
    components.user_backscatter_scale = material.getSurfaceBackscatterScale();
    components.user_photo_scale = material.getSurfacePhotoEmissionScale();
    components.fallback_secondary_yield = result.electron_response.secondary_emitted_electrons;
    components.fallback_ion_secondary_yield = result.ion_response.secondary_emitted_electrons;
    components.fallback_backscatter_yield = result.electron_response.backscattered_electrons;
    result.barrier_outputs =
        FieldSolver::evaluateSurfaceEmissionBarrierBundle(material, request.barrier_state, components);
    result.valid = result.barrier_outputs.valid;
    return result;
}

SurfaceInteractionBundleRequest SurfaceInteraction::makeBundleRequest(
    const SurfaceInteractionEnvironment& environment) const
{
    SurfaceInteractionBundleRequest request;
    request.electron_impact.incident_energy_ev = std::max(1.0e-3, environment.electron_incident_energy_ev);
    request.electron_impact.incident_angle_rad = environment.electron_incident_angle_rad;
    request.electron_impact.particle_charge_coulomb = -1.602176634e-19;
    request.electron_impact.surface_temperature_k = environment.surface_temperature_k;
    request.electron_impact.surface_potential_v = environment.surface_potential_v;
    request.electron_impact.reference_potential_v = environment.reference_potential_v;
    request.electron_impact.normal_electric_field_v_per_m =
        environment.normal_electric_field_v_per_m;
    request.electron_impact.emission_temperature_ev = environment.emission_temperature_ev;
    request.electron_impact.source_current_density_a_per_m2 =
        environment.source_current_density_a_per_m2;
    request.electron_impact.counterelectrode_potential_v =
        environment.counterelectrode_potential_v;
    request.electron_impact.transport_length_m = environment.transport_length_m;
    request.electron_impact.exposed_area_m2 = environment.exposed_area_m2;

    request.ion_impact = request.electron_impact;
    request.ion_impact.incident_energy_ev = std::max(1.0e-3, environment.ion_incident_energy_ev);
    request.ion_impact.incident_angle_rad = environment.ion_incident_angle_rad;
    request.ion_impact.particle_charge_coulomb = 1.602176634e-19;

    request.barrier_state.local_potential_v = environment.surface_potential_v;
    request.barrier_state.reference_potential_v = environment.reference_potential_v;
    request.barrier_state.barrier_potential_v = environment.barrier_potential_v;
    request.barrier_state.normal_electric_field_v_per_m =
        environment.normal_electric_field_v_per_m;
    request.barrier_state.emission_temperature_ev = environment.emission_temperature_ev;
    request.photo_incidence_scale = environment.photo_incidence_scale;
    return request;
}

SurfaceInteractionMaterialSelection SurfaceInteraction::selectMaterialInteraction(
    const MaterialProperty& material)
{
    SurfaceInteractionMaterialSelection selection;
    selection.absorption_probability =
        std::clamp(material.getSurfaceAbsorptionProbability(), 0.0, 1.0);
    selection.reflection_coefficient =
        std::clamp(material.getSurfaceReflectionCoefficient(), 0.0, 0.95);
    selection.emission_scaling = std::max(0.0, material.getSurfaceEmissionScaling());

    if (material.isConductor())
    {
        selection.absorption_probability =
            std::clamp(std::max(0.2, selection.absorption_probability * 0.85), 0.0, 1.0);
        selection.reflection_coefficient =
            std::clamp(std::max(selection.reflection_coefficient, 0.12), 0.0, 0.95);
        selection.interaction_family = "spis_material_indexed_conductor_interaction_v1";
    }
    else if (material.isDielectric())
    {
        selection.absorption_probability =
            std::clamp(std::max(selection.absorption_probability, 0.75), 0.0, 1.0);
        selection.reflection_coefficient =
            std::clamp(std::min(selection.reflection_coefficient, 0.2), 0.0, 0.95);
        selection.interaction_family = "spis_material_indexed_dielectric_interaction_v1";
    }

    return selection;
}

SurfaceInteractionBundleResult SurfaceInteraction::evaluateBundle(
    const MaterialProperty& material, const SurfaceInteractionEnvironment& environment) const
{
    return evaluateBundle(material, makeBundleRequest(environment));
}

SurfaceInteractionBundleResult SurfaceInteraction::evaluateMaterialIndexedBundle(
    const MaterialProperty& material, const SurfaceInteractionEnvironment& environment) const
{
    SurfaceInteraction configured = *this;
    const auto selection = selectMaterialInteraction(material);
    configured.setAbsorptionProbability(selection.absorption_probability);
    configured.setReflectionCoefficient(selection.reflection_coefficient);
    configured.setEmissionScaling(selection.emission_scaling);
    auto result = configured.evaluateBundle(material, environment);
    result.electron_response.material_model_family = selection.interaction_family;
    result.ion_response.material_model_family = selection.interaction_family;
    return result;
}

SurfaceInteractionEnvironment SurfaceInteraction::makeInteractionEnvironment(
    const SurfaceInteractionRoleInputs& inputs)
{
    SurfaceInteractionEnvironment environment;
    environment.electron_incident_energy_ev = std::max(1.0e-3, inputs.electron_incident_energy_ev);
    environment.electron_incident_angle_rad =
        std::clamp(inputs.electron_incident_angle_deg, 0.0, 180.0) * 3.14159265358979323846 / 180.0;
    environment.ion_incident_energy_ev = std::max(1.0e-3, inputs.ion_incident_energy_ev);
    environment.ion_incident_angle_rad =
        std::clamp(inputs.ion_incident_angle_deg, 0.0, 180.0) * 3.14159265358979323846 / 180.0;
    environment.surface_temperature_k = inputs.surface_temperature_k;
    environment.surface_potential_v = inputs.surface_potential_v;
    environment.reference_potential_v = inputs.reference_potential_v;
    environment.normal_electric_field_v_per_m = inputs.normal_electric_field_v_per_m;
    environment.emission_temperature_ev =
        inputs.role == SurfaceInteractionRole::Body
            ? std::max(1.0e-3, inputs.body_photoelectron_temperature_ev)
            : std::max(1.0e-3, inputs.patch_photoelectron_temperature_ev);
    environment.source_current_density_a_per_m2 =
        inputs.role == SurfaceInteractionRole::Body
            ? 0.25 * std::max(0.0, inputs.body_photo_emission_scale) *
                  std::max(0.0, inputs.body_photo_current_density_a_per_m2)
            : std::max(0.0, inputs.patch_photo_current_density_a_per_m2) *
                  std::max(0.0, std::cos(environment.electron_incident_angle_rad));
    environment.counterelectrode_potential_v = inputs.counterelectrode_potential_v;
    environment.transport_length_m =
        inputs.role == SurfaceInteractionRole::Patch ? std::max(1.0e-9, inputs.dielectric_thickness_m)
                                                     : 0.0;
    environment.exposed_area_m2 = std::max(1.0e-12, inputs.exposed_area_m2);
    environment.barrier_potential_v =
        std::max(0.0, inputs.reference_potential_v - inputs.surface_potential_v);
    environment.photo_incidence_scale = inputs.photo_incidence_scale;
    return environment;
}

} // namespace Material
} // namespace SCDAT
