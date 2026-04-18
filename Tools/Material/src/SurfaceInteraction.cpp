#include "../include/SurfaceInteraction.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace SCDAT
{
namespace Material
{

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr std::array<const char*, 19> kSupportedSurfaceInteractorFamilies = {
    "spis_material_model_surface_interactor_v1",
    "spis_material_df_surface_interactor_v1",
    "spis_generic_df_surface_interactor_v1",
    "spis_maxwellian_surface_interactor_v1",
    "spis_maxwellian_surface_interactor_with_recollection_v1",
    "spis_multiple_surface_interactor_v1",
    "spis_multiple_maxwellian_surface_interactor_v1",
    "spis_yield_surface_interactor_v1",
    "spis_reflection_surface_interactor_v1",
    "spis_erosion_surface_interactor_v1",
    "spis_improved_photo_emission_surface_interactor_v1",
    "spis_basic_induced_conduct_surface_interactor_v1",
    "spis_tabulated_sey_surface_interactor_v1",
    "spis_recollected_sey_surface_interactor_v1",
    "spis_default_pee_model_v1",
    "spis_default_seeey_model_v1",
    "spis_default_seep_model_v1",
    "spis_default_erosion_model_v1",
    "spis_device_surface_interactor_v1",
};

enum class ConfiguredInteractorBehavior
{
    MaterialModel,
    Yield,
    MaxwellianWithRecollection,
    Reflection,
    Erosion,
    PhotoEmission,
    InducedConduction,
    TabulatedSey,
    Multiple,
    MultipleMaxwellian,
    RecollectedSey,
};

template <typename Range>
std::string joinFamilyRange(const Range& families)
{
    std::string signature;
    bool first = true;
    for (const auto& family : families)
    {
        if (!first)
        {
            signature += "+";
        }
        first = false;
        signature += family;
    }
    return signature;
}

double clampedCosineFactor(double incident_angle_rad)
{
    const double clamped_angle = std::clamp(incident_angle_rad, 0.0, 0.5 * kPi);
    return std::max(0.05, std::cos(clamped_angle));
}

double incidentCollectionCurrent(const SurfaceComponentCurrents& currents)
{
    return std::max(currents.electron_collection_a_per_m2, currents.ion_collection_a_per_m2);
}

double safeNormalizedCount(double emission_current_a_per_m2, double collection_current_a_per_m2)
{
    return collection_current_a_per_m2 > 0.0
               ? std::max(0.0, emission_current_a_per_m2 / collection_current_a_per_m2)
               : 0.0;
}

double resolveBarrierPotential(const SurfaceImpactState& impact)
{
    return std::max(0.0, impact.reference_potential_v - impact.surface_potential_v);
}

double resolveMaxwellianThermalScale(const MaterialProperty& material,
                                     const SurfaceImpactState& impact)
{
    const double emission_temperature_ev = impact.emission_temperature_ev > 0.0
                                               ? impact.emission_temperature_ev
                                               : material.getPhotoelectronTemperatureEv();
    return std::clamp(0.85 + 0.25 * std::log1p(std::max(1.0e-3, emission_temperature_ev)), 0.85,
                      2.0);
}

double resolveRecollectionFraction(const SurfaceImpactState& impact)
{
    const double barrier_potential_v = resolveBarrierPotential(impact);
    const double barrier_component =
        barrier_potential_v /
        std::max(1.0, std::abs(impact.reference_potential_v) + barrier_potential_v + 1.0);
    const double field_component =
        std::clamp(std::abs(impact.normal_electric_field_v_per_m) / 1.0e5, 0.0, 1.0);
    return std::clamp(0.05 + 0.65 * barrier_component + 0.3 * field_component, 0.0, 0.95);
}

void applyRecollectionAdjustment(SurfaceInteractionResult& result, double collection_current,
                                 double recollection_fraction)
{
    const double clamped_fraction = std::clamp(recollection_fraction, 0.0, 0.95);
    result.component_currents.secondary_electron_emission_a_per_m2 *= (1.0 - clamped_fraction);
    result.component_currents.ion_secondary_electron_emission_a_per_m2 *=
        (1.0 - clamped_fraction);
    result.component_currents.backscatter_emission_a_per_m2 *= (1.0 - clamped_fraction);
    result.secondary_emitted_electrons = safeNormalizedCount(
        result.component_currents.secondary_electron_emission_a_per_m2, collection_current);
    result.ion_secondary_emitted_electrons = safeNormalizedCount(
        result.component_currents.ion_secondary_electron_emission_a_per_m2, collection_current);
    result.backscattered_electrons = safeNormalizedCount(
        result.component_currents.backscatter_emission_a_per_m2, collection_current);
    result.absorbed_charge_coulomb *= (1.0 + 0.1 * clamped_fraction);
    result.deposited_heat_j_per_m2 *= (1.0 + 0.25 * clamped_fraction);
}

double lookupTabulatedSecondaryYield(const MaterialProperty& material, double energy_ev)
{
    const int energy_count = static_cast<int>(
        std::max(0.0, material.getScalarProperty("surface_tabulated_sey_energy_count", 0.0)));
    if (energy_count < 2)
    {
        return std::max(0.0, material.getSecondaryElectronYieldMax());
    }

    std::vector<double> energies;
    std::vector<double> yields;
    for (int i = 0; i < energy_count; ++i)
    {
        energies.push_back(
            material.getScalarProperty("surface_tabulated_sey_energy_" + std::to_string(i) + "_ev",
                                       static_cast<double>(i + 1)));
        yields.push_back(std::max(
            0.0, material.getScalarProperty("surface_tabulated_sey_value_" + std::to_string(i),
                                            material.getSecondaryElectronYieldMax())));
    }

    if (energy_ev <= energies.front())
    {
        return yields.front();
    }
    if (energy_ev >= energies.back())
    {
        return yields.back();
    }

    for (int i = 1; i < energy_count; ++i)
    {
        if (energy_ev <= energies[static_cast<std::size_t>(i)])
        {
            const double e0 = energies[static_cast<std::size_t>(i - 1)];
            const double e1 = energies[static_cast<std::size_t>(i)];
            const double y0 = yields[static_cast<std::size_t>(i - 1)];
            const double y1 = yields[static_cast<std::size_t>(i)];
            const double t = e1 > e0 ? (energy_ev - e0) / (e1 - e0) : 0.0;
            return y0 + (y1 - y0) * std::clamp(t, 0.0, 1.0);
        }
    }
    return yields.back();
}

SurfaceInteractionResult evaluateMaterialModelDrivenInteraction(
    const MaterialProperty& material, const SurfaceImpactState& impact,
    const SurfaceInteractorEvaluationContext& context, const SurfaceMaterialModel* material_model)
{
    const double cosine_factor = clampedCosineFactor(impact.incident_angle_rad);
    const double reflection =
        std::clamp(context.reflection_coefficient * cosine_factor, 0.0, 0.95);
    const double absorption =
        std::clamp(context.absorption_probability * (1.0 - 0.5 * reflection), 0.0, 1.0);

    SurfaceModelContext model_context;
    model_context.incident_particle = impact.particle_charge_coulomb < 0.0
                                          ? SurfaceIncidentParticle::Electron
                                          : SurfaceIncidentParticle::Ion;
    model_context.incident_energy_ev = impact.incident_energy_ev;
    model_context.incident_angle_rad = impact.incident_angle_rad;
    model_context.source_current_density_a_per_m2 = impact.source_current_density_a_per_m2;
    model_context.surface_potential_v = impact.surface_potential_v;
    model_context.reference_potential_v = impact.reference_potential_v;
    model_context.counterelectrode_potential_v = impact.counterelectrode_potential_v;
    model_context.normal_electric_field_v_per_m = impact.normal_electric_field_v_per_m;
    model_context.emission_temperature_ev =
        impact.emission_temperature_ev > 0.0 ? impact.emission_temperature_ev
                                             : material.getPhotoelectronTemperatureEv();
    model_context.transport_length_m = impact.transport_length_m;
    model_context.exposed_area_m2 = std::max(1.0e-12, impact.exposed_area_m2);

    const auto component_currents =
        material_model == nullptr ? SurfaceComponentCurrents{}
                                  : material_model->evaluateCurrents(material, model_context);
    const double incident_collection_a_per_m2 = incidentCollectionCurrent(component_currents);
    const double emission_scale = std::max(0.0, context.emission_scaling);

    SurfaceInteractionResult result;
    result.absorbed_charge_coulomb = impact.particle_charge_coulomb * absorption;
    result.reflected_energy_ev = impact.incident_energy_ev * reflection;
    result.photoelectron_current_a_per_m2 = component_currents.photoelectron_emission_a_per_m2;
    result.induced_conduction_current_a_per_m2 = component_currents.induced_conduction_a_per_m2;
    result.component_currents = component_currents;
    result.component_currents.secondary_electron_emission_a_per_m2 *= emission_scale;
    result.component_currents.ion_secondary_electron_emission_a_per_m2 *= emission_scale;
    result.component_currents.backscatter_emission_a_per_m2 *= emission_scale;
    result.secondary_emitted_electrons = safeNormalizedCount(
        result.component_currents.secondary_electron_emission_a_per_m2, incident_collection_a_per_m2);
    result.ion_secondary_emitted_electrons = safeNormalizedCount(
        result.component_currents.ion_secondary_electron_emission_a_per_m2,
        incident_collection_a_per_m2);
    result.backscattered_electrons = safeNormalizedCount(
        result.component_currents.backscatter_emission_a_per_m2, incident_collection_a_per_m2);
    result.material_model_family =
        material_model == nullptr ? "legacy_surface_interaction_v1" : material_model->modelFamily();
    result.deposited_heat_j_per_m2 = impact.incident_energy_ev * kElementaryCharge *
                                     (absorption + 0.25 * result.secondary_emitted_electrons);
    return result;
}

class ConfiguredSurfaceInteractorPlugin final : public SurfaceInteractorPlugin
{
  public:
    ConfiguredSurfaceInteractorPlugin(std::string family_name, ConfiguredInteractorBehavior behavior)
        : family_name_(std::move(family_name)), behavior_(behavior)
    {
    }

    const char* interactorFamily() const override { return family_name_.c_str(); }

    SurfaceInteractionResult evaluate(
        const MaterialProperty& material, const SurfaceImpactState& impact,
        const SurfaceInteractorEvaluationContext& context) const override
    {
        BasicSurfaceMaterialModel basic_model;
        ErosionSurfaceMaterialModel erosion_model;

        const SurfaceMaterialModel* selected_model =
            behavior_ == ConfiguredInteractorBehavior::Erosion ? &erosion_model
                                                               : context.material_model;
        if (selected_model == nullptr)
        {
            selected_model = &basic_model;
        }

        auto result = evaluateMaterialModelDrivenInteraction(material, impact, context, selected_model);
        const double cosine_factor = clampedCosineFactor(impact.incident_angle_rad);
        const double collection_current =
            std::max(1.0e-18, incidentCollectionCurrent(result.component_currents));

        switch (behavior_)
        {
        case ConfiguredInteractorBehavior::MaterialModel:
            break;
        case ConfiguredInteractorBehavior::Yield:
        {
            const double yield_scale = std::clamp(
                0.35 + material.getSecondaryElectronYieldMax() * (0.4 + 0.6 * cosine_factor), 0.0,
                4.0);
            result.component_currents.secondary_electron_emission_a_per_m2 =
                collection_current * yield_scale;
            result.secondary_emitted_electrons = yield_scale;
            result.deposited_heat_j_per_m2 *= 1.0 + 0.15 * yield_scale;
            break;
        }
        case ConfiguredInteractorBehavior::MaxwellianWithRecollection:
        {
            const double thermal_scale = resolveMaxwellianThermalScale(material, impact);
            const double yield_scale = std::clamp(
                (0.35 + material.getSecondaryElectronYieldMax() * (0.4 + 0.6 * cosine_factor)) *
                    thermal_scale,
                0.0, 4.5);
            result.component_currents.secondary_electron_emission_a_per_m2 =
                collection_current * yield_scale;
            result.secondary_emitted_electrons = yield_scale;
            applyRecollectionAdjustment(result, collection_current,
                                        resolveRecollectionFraction(impact));
            break;
        }
        case ConfiguredInteractorBehavior::Reflection:
        {
            const double reflection_scale =
                std::clamp(1.0 + 4.0 * material.getSurfaceReflectionCoefficient(), 1.0, 3.0);
            result.reflected_energy_ev *= reflection_scale;
            result.absorbed_charge_coulomb *= 1.0 / reflection_scale;
            result.component_currents.backscatter_emission_a_per_m2 *= reflection_scale;
            result.backscattered_electrons = safeNormalizedCount(
                result.component_currents.backscatter_emission_a_per_m2, collection_current);
            break;
        }
        case ConfiguredInteractorBehavior::Erosion:
            result.deposited_heat_j_per_m2 *= 1.25;
            break;
        case ConfiguredInteractorBehavior::PhotoEmission:
        {
            const double photo_boost = std::clamp(
                1.0 + 0.25 * cosine_factor +
                    5.0e4 / std::max(5.0e4, std::abs(impact.normal_electric_field_v_per_m)),
                1.0, 2.5);
            result.photoelectron_current_a_per_m2 *= photo_boost;
            result.component_currents.photoelectron_emission_a_per_m2 =
                result.photoelectron_current_a_per_m2;
            break;
        }
        case ConfiguredInteractorBehavior::InducedConduction:
        {
            const double potential_delta =
                std::abs(impact.counterelectrode_potential_v - impact.surface_potential_v);
            const double conduction_scale = std::clamp(
                material.deriveSurfaceConductivityScaleFactor() *
                    (1.0 + potential_delta /
                               std::max(1.0, std::abs(impact.reference_potential_v) +
                                                 std::abs(impact.surface_potential_v) + 1.0)),
                1.0, 10.0);
            result.induced_conduction_current_a_per_m2 *= conduction_scale;
            result.component_currents.induced_conduction_a_per_m2 =
                result.induced_conduction_current_a_per_m2;
            break;
        }
        case ConfiguredInteractorBehavior::TabulatedSey:
        {
            const double tabulated_yield =
                lookupTabulatedSecondaryYield(material, impact.incident_energy_ev);
            result.component_currents.secondary_electron_emission_a_per_m2 =
                collection_current * std::max(0.0, tabulated_yield);
            result.secondary_emitted_electrons = std::max(0.0, tabulated_yield);
            break;
        }
        case ConfiguredInteractorBehavior::RecollectedSey:
        {
            const double tabulated_yield =
                lookupTabulatedSecondaryYield(material, impact.incident_energy_ev);
            result.component_currents.secondary_electron_emission_a_per_m2 =
                collection_current * std::max(0.0, tabulated_yield);
            result.secondary_emitted_electrons = std::max(0.0, tabulated_yield);
            applyRecollectionAdjustment(result, collection_current,
                                        resolveRecollectionFraction(impact));
            break;
        }
        case ConfiguredInteractorBehavior::Multiple:
        {
            const double yield_weight =
                std::clamp(material.getScalarProperty("surface_multiple_yield_weight", 0.35), 0.0, 1.0);
            const double reflection_weight = std::clamp(
                material.getScalarProperty("surface_multiple_reflection_weight", 0.25), 0.0, 1.0);
            const double photo_weight =
                std::clamp(material.getScalarProperty("surface_multiple_photo_weight", 0.2), 0.0, 1.0);
            const double conduction_weight = std::clamp(
                material.getScalarProperty("surface_multiple_conduction_weight", 0.2), 0.0, 1.0);
            const double total_weight =
                std::max(1.0e-12, yield_weight + reflection_weight + photo_weight + conduction_weight);
            const double yield_scale = std::clamp(
                0.35 + material.getSecondaryElectronYieldMax() * (0.4 + 0.6 * cosine_factor), 0.0,
                4.0);
            const double reflection_scale =
                std::clamp(1.0 + 4.0 * material.getSurfaceReflectionCoefficient(), 1.0, 3.0);
            const double photo_boost = std::clamp(
                1.0 + 0.25 * cosine_factor +
                    5.0e4 / std::max(5.0e4, std::abs(impact.normal_electric_field_v_per_m)),
                1.0, 2.5);
            const double potential_delta =
                std::abs(impact.counterelectrode_potential_v - impact.surface_potential_v);
            const double conduction_scale = std::clamp(
                material.deriveSurfaceConductivityScaleFactor() *
                    (1.0 + potential_delta /
                               std::max(1.0, std::abs(impact.reference_potential_v) +
                                                 std::abs(impact.surface_potential_v) + 1.0)),
                1.0, 10.0);

            result.secondary_emitted_electrons =
                (yield_weight * yield_scale + reflection_weight * result.secondary_emitted_electrons +
                 photo_weight * result.secondary_emitted_electrons +
                 conduction_weight * result.secondary_emitted_electrons) /
                total_weight;
            result.backscattered_electrons =
                (yield_weight * result.backscattered_electrons +
                 reflection_weight * result.backscattered_electrons * reflection_scale +
                 photo_weight * result.backscattered_electrons +
                 conduction_weight * result.backscattered_electrons) /
                total_weight;
            result.photoelectron_current_a_per_m2 =
                (yield_weight * result.photoelectron_current_a_per_m2 +
                 reflection_weight * result.photoelectron_current_a_per_m2 +
                 photo_weight * result.photoelectron_current_a_per_m2 * photo_boost +
                 conduction_weight * result.photoelectron_current_a_per_m2) /
                total_weight;
            result.induced_conduction_current_a_per_m2 =
                (yield_weight * result.induced_conduction_current_a_per_m2 +
                 reflection_weight * result.induced_conduction_current_a_per_m2 +
                 photo_weight * result.induced_conduction_current_a_per_m2 +
                 conduction_weight * result.induced_conduction_current_a_per_m2 * conduction_scale) /
                total_weight;
            result.component_currents.secondary_electron_emission_a_per_m2 =
                result.secondary_emitted_electrons * collection_current;
            result.component_currents.backscatter_emission_a_per_m2 =
                result.backscattered_electrons * collection_current;
            result.component_currents.photoelectron_emission_a_per_m2 =
                result.photoelectron_current_a_per_m2;
            result.component_currents.induced_conduction_a_per_m2 =
                result.induced_conduction_current_a_per_m2;
            break;
        }
        case ConfiguredInteractorBehavior::MultipleMaxwellian:
        {
            const double thermal_scale = resolveMaxwellianThermalScale(material, impact);
            const double yield_weight = std::clamp(
                material.getScalarProperty("surface_multiple_yield_weight", 0.35) *
                    std::sqrt(thermal_scale),
                0.0, 1.0);
            const double reflection_weight = std::clamp(
                material.getScalarProperty("surface_multiple_reflection_weight", 0.25), 0.0,
                1.0);
            const double photo_weight = std::clamp(
                material.getScalarProperty("surface_multiple_photo_weight", 0.2) * thermal_scale,
                0.0, 1.0);
            const double conduction_weight = std::clamp(
                material.getScalarProperty("surface_multiple_conduction_weight", 0.2), 0.0, 1.0);
            const double total_weight =
                std::max(1.0e-12, yield_weight + reflection_weight + photo_weight + conduction_weight);
            const double yield_scale = std::clamp(
                (0.35 + material.getSecondaryElectronYieldMax() * (0.4 + 0.6 * cosine_factor)) *
                    thermal_scale,
                0.0, 4.5);
            const double reflection_scale =
                std::clamp(1.0 + 4.0 * material.getSurfaceReflectionCoefficient(), 1.0, 3.0);
            const double photo_boost = std::clamp(
                (1.0 + 0.25 * cosine_factor +
                 5.0e4 / std::max(5.0e4, std::abs(impact.normal_electric_field_v_per_m))) *
                    std::sqrt(thermal_scale),
                1.0, 3.0);
            const double potential_delta =
                std::abs(impact.counterelectrode_potential_v - impact.surface_potential_v);
            const double conduction_scale = std::clamp(
                material.deriveSurfaceConductivityScaleFactor() *
                    (1.0 + potential_delta /
                               std::max(1.0, std::abs(impact.reference_potential_v) +
                                                 std::abs(impact.surface_potential_v) + 1.0)),
                1.0, 10.0);

            result.secondary_emitted_electrons =
                (yield_weight * yield_scale + reflection_weight * result.secondary_emitted_electrons +
                 photo_weight * result.secondary_emitted_electrons +
                 conduction_weight * result.secondary_emitted_electrons) /
                total_weight;
            result.backscattered_electrons =
                (yield_weight * result.backscattered_electrons +
                 reflection_weight * result.backscattered_electrons * reflection_scale +
                 photo_weight * result.backscattered_electrons +
                 conduction_weight * result.backscattered_electrons) /
                total_weight;
            result.photoelectron_current_a_per_m2 =
                (yield_weight * result.photoelectron_current_a_per_m2 +
                 reflection_weight * result.photoelectron_current_a_per_m2 +
                 photo_weight * result.photoelectron_current_a_per_m2 * photo_boost +
                 conduction_weight * result.photoelectron_current_a_per_m2) /
                total_weight;
            result.induced_conduction_current_a_per_m2 =
                (yield_weight * result.induced_conduction_current_a_per_m2 +
                 reflection_weight * result.induced_conduction_current_a_per_m2 +
                 photo_weight * result.induced_conduction_current_a_per_m2 +
                 conduction_weight * result.induced_conduction_current_a_per_m2 *
                     conduction_scale) /
                total_weight;
            result.component_currents.secondary_electron_emission_a_per_m2 =
                result.secondary_emitted_electrons * collection_current;
            result.component_currents.backscatter_emission_a_per_m2 =
                result.backscattered_electrons * collection_current;
            result.component_currents.photoelectron_emission_a_per_m2 =
                result.photoelectron_current_a_per_m2;
            result.component_currents.induced_conduction_a_per_m2 =
                result.induced_conduction_current_a_per_m2;
            applyRecollectionAdjustment(result, collection_current,
                                        0.5 * resolveRecollectionFraction(impact));
            break;
        }
        }

        result.component_currents.net_current_a_per_m2 =
            result.component_currents.ion_collection_a_per_m2 -
            result.component_currents.electron_collection_a_per_m2 -
            result.component_currents.secondary_electron_emission_a_per_m2 -
            result.component_currents.ion_secondary_electron_emission_a_per_m2 -
            result.component_currents.backscatter_emission_a_per_m2 -
            result.component_currents.photoelectron_emission_a_per_m2 +
            result.component_currents.induced_conduction_a_per_m2;
        return result;
    }

  private:
    std::string family_name_;
    ConfiguredInteractorBehavior behavior_;
};

std::shared_ptr<const SurfaceInteractorPlugin> builtinInteractorPlugin(
    SurfaceInteractorModelVariant variant)
{
    switch (variant)
    {
    case SurfaceInteractorModelVariant::MaterialModel:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_material_model_surface_interactor_v1",
            ConfiguredInteractorBehavior::MaterialModel);
    case SurfaceInteractorModelVariant::MaterialDf:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_material_df_surface_interactor_v1",
            ConfiguredInteractorBehavior::MaterialModel);
    case SurfaceInteractorModelVariant::GenericDf:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_generic_df_surface_interactor_v1", ConfiguredInteractorBehavior::MaterialModel);
    case SurfaceInteractorModelVariant::Maxwellian:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_maxwellian_surface_interactor_v1", ConfiguredInteractorBehavior::Yield);
    case SurfaceInteractorModelVariant::MaxwellianWithRecollection:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_maxwellian_surface_interactor_with_recollection_v1",
            ConfiguredInteractorBehavior::MaxwellianWithRecollection);
    case SurfaceInteractorModelVariant::Multiple:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_multiple_surface_interactor_v1", ConfiguredInteractorBehavior::Multiple);
    case SurfaceInteractorModelVariant::MultipleMaxwellian:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_multiple_maxwellian_surface_interactor_v1",
            ConfiguredInteractorBehavior::MultipleMaxwellian);
    case SurfaceInteractorModelVariant::Yield:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_yield_surface_interactor_v1", ConfiguredInteractorBehavior::Yield);
    case SurfaceInteractorModelVariant::Reflection:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_reflection_surface_interactor_v1", ConfiguredInteractorBehavior::Reflection);
    case SurfaceInteractorModelVariant::Erosion:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_erosion_surface_interactor_v1", ConfiguredInteractorBehavior::Erosion);
    case SurfaceInteractorModelVariant::ImprovedPhotoEmission:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_improved_photo_emission_surface_interactor_v1",
            ConfiguredInteractorBehavior::PhotoEmission);
    case SurfaceInteractorModelVariant::BasicInducedConduction:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_basic_induced_conduct_surface_interactor_v1",
            ConfiguredInteractorBehavior::InducedConduction);
    case SurfaceInteractorModelVariant::TabulatedSey:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_tabulated_sey_surface_interactor_v1",
            ConfiguredInteractorBehavior::TabulatedSey);
    case SurfaceInteractorModelVariant::RecollectedSey:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_recollected_sey_surface_interactor_v1",
            ConfiguredInteractorBehavior::RecollectedSey);
    case SurfaceInteractorModelVariant::DefaultPee:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_default_pee_model_v1", ConfiguredInteractorBehavior::PhotoEmission);
    case SurfaceInteractorModelVariant::DefaultSeee:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_default_seeey_model_v1", ConfiguredInteractorBehavior::Yield);
    case SurfaceInteractorModelVariant::DefaultSeep:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_default_seep_model_v1", ConfiguredInteractorBehavior::Yield);
    case SurfaceInteractorModelVariant::DefaultErosion:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_default_erosion_model_v1", ConfiguredInteractorBehavior::Erosion);
    case SurfaceInteractorModelVariant::Device:
        return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
            "spis_device_surface_interactor_v1", ConfiguredInteractorBehavior::MaterialModel);
    }
    return std::make_shared<ConfiguredSurfaceInteractorPlugin>(
        "spis_material_model_surface_interactor_v1", ConfiguredInteractorBehavior::MaterialModel);
}

} // namespace

const char* MaterialModelSurfaceInteractor::interactorFamily() const
{
    return "spis_material_model_surface_interactor_v1";
}

SurfaceInteractionResult MaterialModelSurfaceInteractor::evaluate(
    const MaterialProperty& material, const SurfaceImpactState& impact,
    const SurfaceInteractorEvaluationContext& context) const
{
    return evaluateMaterialModelDrivenInteraction(material, impact, context, context.material_model);
}

SurfaceInteraction::SurfaceInteraction()
    : material_model_(&default_material_model_),
      default_interactor_plugin_(std::make_shared<MaterialModelSurfaceInteractor>()),
      interactor_plugin_(default_interactor_plugin_)
{
}

void SurfaceInteraction::setInteractorPlugin(std::shared_ptr<const SurfaceInteractorPlugin> plugin)
{
    interactor_plugin_ = std::move(plugin);
    if (!interactor_plugin_)
    {
        interactor_plugin_ = default_interactor_plugin_;
    }
}

void SurfaceInteraction::resetInteractorPlugin()
{
    interactor_plugin_ = default_interactor_plugin_;
}

const char* SurfaceInteraction::materialModelFamily() const
{
    return material_model_ == nullptr ? "legacy_surface_interaction_v1"
                                      : material_model_->modelFamily();
}

const char* SurfaceInteraction::interactorFamily() const
{
    return interactor_plugin_ == nullptr ? "legacy_surface_interaction_v1"
                                         : interactor_plugin_->interactorFamily();
}

const char* SurfaceInteraction::interactorFamily(const MaterialProperty& material) const
{
    if (interactor_plugin_ && interactor_plugin_ != default_interactor_plugin_)
    {
        return interactor_plugin_->interactorFamily();
    }
    return resolveSurfaceInteractorModelFamily(material);
}

SurfaceInteractorFamilyView SurfaceInteraction::interactorFamilyView(
    const MaterialProperty& material) const
{
    SurfaceInteractorFamilyView view;
    view.supported_family_count = kSupportedSurfaceInteractorFamilies.size();
    view.supported_family_signature = joinFamilyRange(kSupportedSurfaceInteractorFamilies);

    const char* active_family = interactorFamily(material);
    if (active_family != nullptr && std::string(active_family).size() > 0)
    {
        view.active_family_count = 1;
        view.active_family_signature = active_family;
    }
    return view;
}

SurfaceInteractionResult SurfaceInteraction::evaluate(const MaterialProperty& material,
                                                      const SurfaceImpactState& impact) const
{
    SurfaceInteractorEvaluationContext context;
    context.material_model = material_model_;
    context.absorption_probability = absorption_probability_;
    context.reflection_coefficient = reflection_coefficient_;
    context.emission_scaling = emission_scaling_;

    if (!interactor_plugin_ || interactor_plugin_ == default_interactor_plugin_)
    {
        return builtinInteractorPlugin(resolveSurfaceInteractorModelVariant(material))
            ->evaluate(material, impact, context);
    }
    return interactor_plugin_->evaluate(material, impact, context);
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
    request.electron_impact.particle_charge_coulomb = -kElementaryCharge;
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
    request.ion_impact.particle_charge_coulomb = kElementaryCharge;

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
        std::clamp(inputs.electron_incident_angle_deg, 0.0, 180.0) * kPi / 180.0;
    environment.ion_incident_energy_ev = std::max(1.0e-3, inputs.ion_incident_energy_ev);
    environment.ion_incident_angle_rad =
        std::clamp(inputs.ion_incident_angle_deg, 0.0, 180.0) * kPi / 180.0;
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

SurfaceInteractorModelVariant resolveSurfaceInteractorModelVariant(
    const MaterialProperty& material)
{
    const int explicit_code = static_cast<int>(
        material.getScalarProperty("surface_interactor_family_code", -1.0));
    if (explicit_code >= 0 && explicit_code <= static_cast<int>(SurfaceInteractorModelVariant::Device))
    {
        return static_cast<SurfaceInteractorModelVariant>(explicit_code);
    }

    if (material.getScalarProperty("surface_interactor_use_multiple_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::Multiple;
    }
    if (material.getScalarProperty("surface_interactor_use_tabulated_sey_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::TabulatedSey;
    }
    if (material.getScalarProperty("surface_interactor_use_erosion_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::Erosion;
    }
    if (material.getScalarProperty("surface_interactor_use_reflection_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::Reflection;
    }
    if (material.getScalarProperty("surface_interactor_use_yield_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::Yield;
    }
    if (material.getScalarProperty("surface_interactor_use_photoemission_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::ImprovedPhotoEmission;
    }
    if (material.getScalarProperty("surface_interactor_use_induced_conduct_model", 0.0) > 0.5)
    {
        return SurfaceInteractorModelVariant::BasicInducedConduction;
    }
    return SurfaceInteractorModelVariant::MaterialModel;
}

const char* resolveSurfaceInteractorModelFamily(const MaterialProperty& material)
{
    switch (resolveSurfaceInteractorModelVariant(material))
    {
    case SurfaceInteractorModelVariant::MaterialModel:
        return "spis_material_model_surface_interactor_v1";
    case SurfaceInteractorModelVariant::MaterialDf:
        return "spis_material_df_surface_interactor_v1";
    case SurfaceInteractorModelVariant::GenericDf:
        return "spis_generic_df_surface_interactor_v1";
    case SurfaceInteractorModelVariant::Maxwellian:
        return "spis_maxwellian_surface_interactor_v1";
    case SurfaceInteractorModelVariant::MaxwellianWithRecollection:
        return "spis_maxwellian_surface_interactor_with_recollection_v1";
    case SurfaceInteractorModelVariant::Multiple:
        return "spis_multiple_surface_interactor_v1";
    case SurfaceInteractorModelVariant::MultipleMaxwellian:
        return "spis_multiple_maxwellian_surface_interactor_v1";
    case SurfaceInteractorModelVariant::Yield:
        return "spis_yield_surface_interactor_v1";
    case SurfaceInteractorModelVariant::Reflection:
        return "spis_reflection_surface_interactor_v1";
    case SurfaceInteractorModelVariant::Erosion:
        return "spis_erosion_surface_interactor_v1";
    case SurfaceInteractorModelVariant::ImprovedPhotoEmission:
        return "spis_improved_photo_emission_surface_interactor_v1";
    case SurfaceInteractorModelVariant::BasicInducedConduction:
        return "spis_basic_induced_conduct_surface_interactor_v1";
    case SurfaceInteractorModelVariant::TabulatedSey:
        return "spis_tabulated_sey_surface_interactor_v1";
    case SurfaceInteractorModelVariant::RecollectedSey:
        return "spis_recollected_sey_surface_interactor_v1";
    case SurfaceInteractorModelVariant::DefaultPee:
        return "spis_default_pee_model_v1";
    case SurfaceInteractorModelVariant::DefaultSeee:
        return "spis_default_seeey_model_v1";
    case SurfaceInteractorModelVariant::DefaultSeep:
        return "spis_default_seep_model_v1";
    case SurfaceInteractorModelVariant::DefaultErosion:
        return "spis_default_erosion_model_v1";
    case SurfaceInteractorModelVariant::Device:
        return "spis_device_surface_interactor_v1";
    }
    return "spis_material_model_surface_interactor_v1";
}

} // namespace Material
} // namespace SCDAT
