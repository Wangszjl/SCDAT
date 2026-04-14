#include "../include/SurfaceMaterialModel.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Material
{

namespace
{

constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kPi = 3.14159265358979323846;
constexpr double kAtomicMassUnitKg = 1.66053906660e-27;
constexpr double kBoltzmannEvPerK = 8.617333262145e-5;
constexpr double kRichardsonConstant = 1.20173e6;
constexpr double kEpsilon0 = 8.8541878128e-12;

double clampAngle(double angle_rad)
{
    return std::clamp(angle_rad, 0.0, 0.5 * kPi);
}

double secondaryYield(const MaterialProperty& material, double energy_ev, double angle_rad)
{
    const double peak_energy_ev = std::max(1.0, material.getSecondaryElectronPeakEnergyEv());
    const double x = std::max(1.0e-6, energy_ev / peak_energy_ev);
    const double cosine_factor = std::max(0.1, std::cos(clampAngle(angle_rad)));
    return std::max(0.0, material.getSecondaryElectronYieldMax() * x * std::exp(1.0 - x) /
                             std::sqrt(cosine_factor));
}

double ionSecondaryYield(const MaterialProperty& material, double energy_ev, double angle_rad)
{
    const double peak_energy_ev = std::max(1.0, material.getProtonSecondaryElectronPeakEnergyEv());
    const double x = std::max(1.0e-6, energy_ev / peak_energy_ev);
    const double angular_boost = 1.0 + 0.15 * std::sin(clampAngle(angle_rad));
    return std::max(0.0, material.getProtonSecondaryElectronYieldMax() * std::sqrt(x) *
                             std::exp(0.5 - 0.5 * x) * angular_boost);
}

double backscatterYield(const MaterialProperty& material, double energy_ev, double angle_rad)
{
    const double dielectric_floor = material.isDielectric() ? 0.5 : 0.15;
    const double conductivity_ratio =
        material.getSurfaceConductivitySPerM() < 1.0e3 ? dielectric_floor : 0.0;
    const double angular_factor = 1.0 + 0.2 * std::sin(clampAngle(angle_rad));
    const double reduced_energy = std::min(1.0, std::max(0.0, energy_ev / 2.0e3));
    return std::clamp((material.getElectronBackscatterYield() + conductivity_ratio) *
                          (0.35 + 0.65 * reduced_energy) * angular_factor,
                      0.0, 0.95);
}

double photoCurrentDensity(const MaterialProperty& material, const SurfaceModelContext& context)
{
    const double base_current_density =
        std::max(0.0, context.source_current_density_a_per_m2) *
        std::max(0.0, material.getPhotoelectronYield());
    const double thermal_tail =
        std::exp(-std::max(0.0, context.surface_potential_v - context.reference_potential_v) /
                 std::max(1.0e-3, material.getPhotoelectronTemperatureEv()));
    return base_current_density * thermal_tail;
}

double inducedConductCurrentDensity(const MaterialProperty& material,
                                    const SurfaceModelContext& context)
{
    if (context.transport_length_m <= 0.0)
    {
        return 0.0;
    }
    const double conductivity = std::max(material.getDarkConductivitySPerM(),
                                         material.getSurfaceConductivitySPerM());
    return conductivity *
           (context.counterelectrode_potential_v - context.surface_potential_v) /
           std::max(1.0e-9, context.transport_length_m);
}

double legacyThermalCurrentDensity(double density_m3, double temperature_ev, double mass_kg)
{
    return kElementaryCharge * std::max(0.0, density_m3) *
           std::sqrt(std::max(1.0e-6, temperature_ev) * kElementaryCharge /
                     (2.0 * kPi * std::max(1.0e-30, mass_kg)));
}

double legacyWhippleYield(double peak_energy_ev, double peak_yield, double incident_energy_ev)
{
    if (incident_energy_ev <= 0.0 || peak_yield <= 0.0)
    {
        return 0.0;
    }
    const double q = 2.28 * std::pow(incident_energy_ev / std::max(1.0, peak_energy_ev), 1.35);
    if (q <= 1.0e-12)
    {
        return 0.0;
    }
    return 2.228 * (q - 1.0 + std::exp(-q)) / q * peak_yield *
           std::pow(std::max(1.0, peak_energy_ev) / incident_energy_ev, 0.35);
}

double legacySimsYield(double peak_energy_ev, double peak_yield, double exponent_n,
                       double incident_energy_ev)
{
    if (incident_energy_ev <= 10.0 || peak_yield <= 0.0)
    {
        return 0.0;
    }
    const double n = std::max(1.05, exponent_n);
    double xm = 2.5;
    for (int i = 1; i < 10001; ++i)
    {
        const double x = 0.5 + static_cast<double>(i) * (2.5 - 0.5) / 10000.0;
        const double y = (1.0 - 1.0 / n) * (std::exp(x) - 1.0);
        if (std::abs(x - y) <= 1.0e-3)
        {
            xm = x;
            break;
        }
    }
    const double x = xm * std::pow(incident_energy_ev / std::max(1.0, peak_energy_ev), n);
    if (x <= 1.0e-12 || std::abs(1.0 - std::exp(-xm)) < 1.0e-12)
    {
        return 0.0;
    }
    return peak_yield / (1.0 - std::exp(-xm)) *
           std::pow(std::max(1.0, peak_energy_ev) / incident_energy_ev, n - 1.0) * 2.0 *
           (x + std::exp(-x) - 1.0) / x;
}

double legacyKatzYield(const MaterialProperty& material, double incident_energy_ev)
{
    if (incident_energy_ev <= 1.0e-9)
    {
        return 0.0;
    }
    const double peak_energy_ev =
        std::max(1.0, material.getScalarProperty("secondary_yield_peak_energy_ev", 400.0));
    const double peak_yield = std::max(0.0, material.getSecondaryElectronYield());
    const double r1 = material.getScalarProperty("katz_r1", 0.85);
    const double n1 = material.getScalarProperty("katz_n1", 0.35);
    const double r2 = material.getScalarProperty("katz_r2", 0.12);
    const double n2 = material.getScalarProperty("katz_n2", 1.05);
    const double energy_kev = incident_energy_ev * 1.0e-3;
    const double peak_energy_kev = peak_energy_ev * 1.0e-3;
    const double range = r1 * std::pow(energy_kev, n1) + r2 * std::pow(energy_kev, n2);
    const double peak_range =
        r1 * std::pow(peak_energy_kev, n1) + r2 * std::pow(peak_energy_kev, n2);
    if (range <= 1.0e-12 || peak_range <= 1.0e-12)
    {
        return 0.0;
    }
    const double attenuation = std::exp(-std::pow(range / peak_range - 1.0, 2.0));
    const double rolloff =
        1.0 / (1.0 + std::pow(incident_energy_ev / std::max(1.0, 3.0 * peak_energy_ev), 0.7));
    return peak_yield * attenuation * (0.55 + 0.45 * rolloff);
}

double safeExp(double exponent)
{
    return std::exp(std::clamp(exponent, -700.0, 700.0));
}

bool erosionModelEnabled(const MaterialProperty& material)
{
    const double explicit_toggle =
        material.getScalarProperty("surface_material_model_use_erosion", -1.0);
    if (explicit_toggle >= 0.0)
    {
        return explicit_toggle > 0.5;
    }

    if (material.getScalarProperty("surface_material_model_variant", 0.0) >= 0.5)
    {
        return true;
    }

    return material.getScalarProperty("erosion_yield_scale", 0.0) > 0.0 ||
           material.getScalarProperty("surface_erosion_yield_scale", 0.0) > 0.0 ||
           material.getScalarProperty("erosion_rate_scale", 0.0) > 0.0;
}

double erosionActivation(const MaterialProperty& material, const SurfaceModelContext& context)
{
    const double yield_scale =
        std::max(material.getScalarProperty("erosion_yield_scale", 0.0),
                 material.getScalarProperty("surface_erosion_yield_scale", 0.0));
    if (yield_scale <= 0.0)
    {
        return 0.0;
    }

    const double threshold_energy_ev =
        std::max(0.0, material.getScalarProperty("erosion_threshold_energy_ev", 80.0));
    const double characteristic_energy_ev =
        std::max(1.0, material.getScalarProperty("erosion_characteristic_energy_ev", 500.0));
    const double energy_excess_ev = std::max(0.0, context.incident_energy_ev - threshold_energy_ev);
    const double energy_factor = std::clamp(energy_excess_ev / characteristic_energy_ev, 0.0, 1.0);
    const double angle_factor = 1.0 + 0.25 * std::sin(clampAngle(context.incident_angle_rad));
    return std::clamp(yield_scale * energy_factor * angle_factor, 0.0, 1.0);
}

void applyErosionResponse(const MaterialProperty& material, const SurfaceModelContext& context,
                         SurfaceComponentCurrents& currents)
{
    const double activation = erosionActivation(material, context);
    if (activation <= 0.0)
    {
        return;
    }

    const double secondary_gain =
        std::max(0.0, material.getScalarProperty("erosion_secondary_gain", 0.35));
    const double ion_secondary_gain =
        std::max(0.0, material.getScalarProperty("erosion_ion_secondary_gain", 0.80));
    const double backscatter_gain =
        std::max(0.0, material.getScalarProperty("erosion_backscatter_gain", 0.25));
    const double photo_suppression =
        std::clamp(material.getScalarProperty("erosion_photo_suppression", 0.25), 0.0, 1.0);
    const double conduction_suppression =
        std::clamp(material.getScalarProperty("erosion_conductivity_suppression", 0.40), 0.0,
                   1.0);

    currents.secondary_electron_emission_a_per_m2 *= 1.0 + activation * secondary_gain;
    currents.ion_secondary_electron_emission_a_per_m2 *=
        1.0 + activation * ion_secondary_gain;
    currents.backscatter_emission_a_per_m2 *= 1.0 + activation * backscatter_gain;
    currents.photoelectron_emission_a_per_m2 *=
        std::max(0.0, 1.0 - activation * photo_suppression);
    currents.induced_conduction_a_per_m2 *=
        std::max(0.0, 1.0 - activation * conduction_suppression);
}

} // namespace

const char* BasicSurfaceMaterialModel::modelFamily() const
{
    return "spis_basic_surface_material_model_v1";
}

SurfaceComponentCurrents BasicSurfaceMaterialModel::evaluateCurrents(
    const MaterialProperty& material, const SurfaceModelContext& context) const
{
    SurfaceComponentCurrents currents;

    const double incident_flux_a_per_m2 = kElementaryCharge * std::max(0.0, context.exposed_area_m2) *
                                          std::max(0.0, context.incident_energy_ev) * 1.0e-6;
    switch (context.incident_particle)
    {
    case SurfaceIncidentParticle::Electron:
        currents.electron_collection_a_per_m2 = incident_flux_a_per_m2;
        currents.secondary_electron_emission_a_per_m2 =
            incident_flux_a_per_m2 *
            secondaryYield(material, context.incident_energy_ev, context.incident_angle_rad);
        currents.backscatter_emission_a_per_m2 =
            incident_flux_a_per_m2 *
            backscatterYield(material, context.incident_energy_ev, context.incident_angle_rad);
        break;
    case SurfaceIncidentParticle::Ion:
        currents.ion_collection_a_per_m2 = incident_flux_a_per_m2;
        currents.ion_secondary_electron_emission_a_per_m2 =
            incident_flux_a_per_m2 *
            ionSecondaryYield(material, context.incident_energy_ev, context.incident_angle_rad);
        break;
    case SurfaceIncidentParticle::Photon:
        break;
    }

    currents.photoelectron_emission_a_per_m2 = photoCurrentDensity(material, context);
    currents.induced_conduction_a_per_m2 = inducedConductCurrentDensity(material, context);
    currents.net_current_a_per_m2 = currents.electron_collection_a_per_m2 +
                                    currents.ion_collection_a_per_m2 -
                                    currents.secondary_electron_emission_a_per_m2 -
                                    currents.ion_secondary_electron_emission_a_per_m2 -
                                    currents.backscatter_emission_a_per_m2 -
                                    currents.photoelectron_emission_a_per_m2 +
                                    currents.induced_conduction_a_per_m2;
    return currents;
}

const char* ErosionSurfaceMaterialModel::modelFamily() const
{
    return "spis_erosion_surface_material_model_v1";
}

SurfaceComponentCurrents ErosionSurfaceMaterialModel::evaluateCurrents(
    const MaterialProperty& material, const SurfaceModelContext& context) const
{
    BasicSurfaceMaterialModel basic_model;
    auto currents = basic_model.evaluateCurrents(material, context);
    applyErosionResponse(material, context, currents);
    currents.net_current_a_per_m2 = currents.electron_collection_a_per_m2 +
                                    currents.ion_collection_a_per_m2 -
                                    currents.secondary_electron_emission_a_per_m2 -
                                    currents.ion_secondary_electron_emission_a_per_m2 -
                                    currents.backscatter_emission_a_per_m2 -
                                    currents.photoelectron_emission_a_per_m2 +
                                    currents.induced_conduction_a_per_m2;
    return currents;
}

SurfaceMaterialModelVariant resolveSurfaceMaterialModelVariant(
    const MaterialProperty& material)
{
    return erosionModelEnabled(material) ? SurfaceMaterialModelVariant::Erosion
                                         : SurfaceMaterialModelVariant::Basic;
}

const char* surfaceMaterialModelVariantFamily(SurfaceMaterialModelVariant variant)
{
    switch (variant)
    {
    case SurfaceMaterialModelVariant::Erosion:
        return "spis_erosion_surface_material_model_v1";
    case SurfaceMaterialModelVariant::Basic:
    default:
        return "spis_basic_surface_material_model_v1";
    }
}

const char* resolveSurfaceMaterialModelFamily(const MaterialProperty& material)
{
    return surfaceMaterialModelVariantFamily(resolveSurfaceMaterialModelVariant(material));
}

SurfaceRoleCurrentContexts makeSurfaceRoleCurrentContexts(
    const SurfaceRoleCurrentInputs& inputs)
{
    SurfaceRoleCurrentContexts contexts;
    contexts.electron_context.incident_particle = SurfaceIncidentParticle::Electron;
    contexts.electron_context.incident_energy_ev =
        std::max(1.0e-3,
                 inputs.electron_characteristic_energy_ev +
                     std::max(0.0, inputs.surface_potential_v));
    contexts.electron_context.incident_angle_rad =
        std::clamp(inputs.patch_incidence_angle_deg, 0.0, 180.0) * kPi / 180.0;
    contexts.electron_context.surface_potential_v = inputs.surface_potential_v;
    contexts.electron_context.reference_potential_v = inputs.reference_potential_v;
    contexts.electron_context.normal_electric_field_v_per_m =
        inputs.normal_electric_field_v_per_m;
    contexts.electron_context.emission_temperature_ev =
        inputs.role == SurfaceCurrentRole::Body
            ? std::max(1.0e-3, inputs.body_photoelectron_temperature_ev)
            : std::max(1.0e-3, inputs.patch_photoelectron_temperature_ev);
    contexts.electron_context.exposed_area_m2 = std::max(1.0e-12, inputs.exposed_area_m2);

    contexts.ion_context = contexts.electron_context;
    contexts.ion_context.incident_particle = SurfaceIncidentParticle::Ion;
    contexts.ion_context.incident_energy_ev =
        std::max(1.0e-3,
                 inputs.ion_characteristic_energy_ev +
                     std::max(0.0, -inputs.surface_potential_v));
    contexts.ion_context.incident_angle_rad =
        std::clamp(inputs.patch_flow_angle_deg, 0.0, 180.0) * kPi / 180.0;

    contexts.photo_context = contexts.electron_context;
    contexts.photo_context.incident_particle = SurfaceIncidentParticle::Photon;
    contexts.photo_context.source_current_density_a_per_m2 =
        inputs.role == SurfaceCurrentRole::Body
            ? 0.25 * std::max(0.0, inputs.body_photo_emission_scale) *
                  std::max(0.0, inputs.body_photo_current_density_a_per_m2)
            : std::max(0.0, inputs.patch_photo_current_density_a_per_m2) *
                  std::max(0.0, std::cos(contexts.electron_context.incident_angle_rad));
    contexts.photo_context.counterelectrode_potential_v =
        inputs.counterelectrode_potential_v;
    contexts.photo_context.transport_length_m =
        inputs.role == SurfaceCurrentRole::Patch
            ? std::max(1.0e-9, inputs.dielectric_thickness_m)
            : 0.0;
    return contexts;
}

SurfaceRoleCurrentBundle evaluateSurfaceRoleCurrentBundle(
    const MaterialProperty& material, const SurfaceRoleCurrentInputs& inputs,
    const SurfaceMaterialModel& model)
{
    SurfaceRoleCurrentBundle bundle;
    bundle.contexts = makeSurfaceRoleCurrentContexts(inputs);
    bundle.electron_currents =
        model.evaluateCurrents(material, bundle.contexts.electron_context);
    bundle.ion_currents = model.evaluateCurrents(material, bundle.contexts.ion_context);

    auto conductive_material = material;
    conductive_material.setSurfaceConductivitySPerM(
        std::max(0.0, inputs.effective_conductivity_s_per_m) *
        std::max(0.0, inputs.conductivity_scale));
    bundle.photo_currents =
        model.evaluateCurrents(conductive_material, bundle.contexts.photo_context);
    return bundle;
}

SurfaceRoleCurrentBundle evaluateSurfaceRoleCurrentBundle(
    const MaterialProperty& material, const SurfaceRoleCurrentInputs& inputs)
{
    BasicSurfaceMaterialModel model;
    return evaluateSurfaceRoleCurrentBundle(material, inputs, model);
}

double evaluateLegacySecondaryElectronYield(const MaterialProperty& material,
                                            LegacySecondaryYieldModel model,
                                            double incident_energy_ev,
                                            double sims_exponent_n)
{
    const double peak_energy_ev =
        std::max(1.0, material.getScalarProperty("secondary_yield_peak_energy_ev", 300.0));
    const double peak_yield = std::max(0.0, material.getSecondaryElectronYield());
    switch (model)
    {
    case LegacySecondaryYieldModel::Whipple:
        return legacyWhippleYield(peak_energy_ev, peak_yield, incident_energy_ev);
    case LegacySecondaryYieldModel::Sims:
        return legacySimsYield(peak_energy_ev, peak_yield, sims_exponent_n, incident_energy_ev);
    case LegacySecondaryYieldModel::Katz:
    default:
        return legacyKatzYield(material, incident_energy_ev);
    }
}

double evaluateLegacyIonSecondaryElectronYield(const MaterialProperty& material,
                                               double incident_energy_ev)
{
    const double peak_energy_kev =
        std::max(1.0e-6, material.getScalarProperty("ion_secondary_peak_energy_kev", 0.35));
    const double energy_kev = std::max(0.0, incident_energy_ev * 1.0e-3);
    const double peak_yield =
        std::max(0.0, material.getScalarProperty("ion_secondary_yield", 0.08));
    if (energy_kev <= 0.0)
    {
        return 0.0;
    }
    if (energy_kev < 0.476)
    {
        return 0.5 * std::pow(2.0, -2.0) * peak_yield * std::sqrt(energy_kev) /
               (1.0 + energy_kev / peak_energy_kev);
    }
    if (energy_kev <= 10.0)
    {
        const double exponent = -(1.0 / energy_kev - 0.1);
        return (2.0 - (1.0 / energy_kev - 0.1) / 2.0) * std::pow(2.0, exponent) * peak_yield *
               std::sqrt(energy_kev) / (1.0 + energy_kev / peak_energy_kev);
    }
    return 2.0 * peak_yield * std::sqrt(energy_kev) / (1.0 + energy_kev / peak_energy_kev);
}

double evaluateLegacyBackscatterYield(const MaterialProperty& material,
                                      double incident_energy_ev)
{
    const double z = material.getScalarProperty("atomic_number", 13.0);
    const double es = incident_energy_ev;
    double ybe = 0.0;
    if (es < 50.0)
    {
        ybe = 0.0;
    }
    else if (es < 1.0e3)
    {
        ybe = 0.3338 * std::log(es / 50.0) *
              (1.0 - std::pow(0.7358, 0.037 * z) + 0.1 * std::exp(-es / 5000.0));
    }
    else if (es < 1.0e4)
    {
        ybe = 1.0 - std::pow(0.7358, 0.037 * z) + 0.1 * std::exp(-es / 5000.0);
    }
    else if (es < 1.0e5)
    {
        ybe = 1.0 - std::pow(0.7358, 0.037 * z);
    }
    if (ybe <= 0.0 || std::abs(std::log(ybe)) < 1.0e-12)
    {
        return 0.0;
    }
    return 2.0 * (1.0 - ybe + ybe * std::log(ybe)) / (std::log(ybe) * std::log(ybe));
}

double evaluateLegacyRamPatchCurrentDensity(double surface_potential_v,
                                            double ion_density_m3,
                                            double ion_temperature_ev,
                                            double ion_mass_amu,
                                            double projected_flow_speed_m_per_s)
{
    const double ni_cm3 = std::max(0.0, ion_density_m3 * 1.0e-6);
    const double ti_ev = std::max(1.0e-6, ion_temperature_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / std::max(1.0, ion_mass_amu));
    const double qd = projected_flow_speed_m_per_s * 1.0e-3 / std::max(1.0e-12, wt);
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    double current_na_per_m2 = 0.08011 * ni_cm3 * wt;
    if (surface_potential_v >= 0.0)
    {
        current_na_per_m2 *=
            0.5642 * std::exp(-(qd - qv) * (qd - qv)) + qd + qd * std::erf(qd - qv);
    }
    else
    {
        current_na_per_m2 *= 0.5642 * std::exp(-qd * qd) + qd + qd * std::erf(qd);
    }
    return std::max(0.0, current_na_per_m2) * 1.0e-9;
}

double evaluateLegacyRamBodyCurrentDensity(double surface_potential_v,
                                           double ion_density_m3,
                                           double ion_temperature_ev,
                                           double ion_mass_amu,
                                           double flow_speed_m_per_s)
{
    const double ni_cm3 = std::max(0.0, ion_density_m3 * 1.0e-6);
    const double ti_ev = std::max(1.0e-6, ion_temperature_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / std::max(1.0, ion_mass_amu));
    const double qd = std::max(1.0e-12, flow_speed_m_per_s * 1.0e-3 / std::max(1.0e-12, wt));
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    double current_na_per_m2 = 0.0;
    if (surface_potential_v >= 0.0)
    {
        current_na_per_m2 =
            0.5642 / qd *
                ((qv + qd) * std::exp(-(qv - qd) * (qv - qd)) -
                 (qv - qd) * std::exp(-(qv + qd) * (qv + qd))) +
            (0.5 / qd + qd - qv * qv / qd) *
                (std::erf(qv + qd) - std::erf(qv - qd));
        current_na_per_m2 *= 0.02003 * ni_cm3 * wt;
    }
    else
    {
        current_na_per_m2 = 0.04005 * ni_cm3 * wt *
                            (0.5642 * std::exp(-qd * qd) +
                             (qd + 0.5 / qd) * std::erf(qd));
        current_na_per_m2 *= 4.0;
    }
    return std::max(0.0, current_na_per_m2) * 1.0e-9;
}

double integrateMaxwellianYield(double temperature_ev,
                                const std::function<double(double)>& yield_evaluator)
{
    if (temperature_ev <= 0.0)
    {
        return 0.0;
    }
    const double max_energy_ev = 20.0 * temperature_ev;
    const std::size_t steps = 400;
    const double step = max_energy_ev / static_cast<double>(steps);
    double weighted_yield = 0.0;
    double weighted_flux = 0.0;
    for (std::size_t index = 0; index <= steps; ++index)
    {
        const double energy_ev = step * static_cast<double>(index);
        const double weight = std::sqrt(std::max(0.0, energy_ev)) *
                              std::exp(-energy_ev / std::max(1.0e-6, temperature_ev));
        const double trapezoid_weight =
            (index == 0 || index == steps) ? 0.5 : 1.0;
        weighted_yield += trapezoid_weight * weight * yield_evaluator(energy_ev);
        weighted_flux += trapezoid_weight * weight;
    }
    return weighted_flux > 0.0 ? weighted_yield / weighted_flux : 0.0;
}

double evaluateLegacyEmissionEscapeProbability(const MaterialProperty& material,
                                               double surface_potential_v)
{
    return surface_potential_v > 0.0
               ? std::exp(-surface_potential_v / std::max(
                                                1.0e-3,
                                                material.getScalarProperty(
                                                    "secondary_emission_escape_energy_ev", 2.0)))
               : 1.0;
}

double evaluateSurfaceFieldEmissionCurrentDensity(const MaterialProperty& material,
                                                  double surface_potential_v,
                                                  double reference_potential_v,
                                                  double normal_electric_field_v_per_m,
                                                  double effective_sheath_length_m)
{
    const double work_function = std::max(0.1, material.getWorkFunctionEv());
    const double sheath_length = std::max(1.0e-6, effective_sheath_length_m);
    const double reference_field_v_per_m =
        std::abs(reference_potential_v - surface_potential_v) / sheath_length;
    const double local_field_v_per_m =
        std::max(std::abs(normal_electric_field_v_per_m), reference_field_v_per_m);
    if (local_field_v_per_m <= 1.0e7)
    {
        return 0.0;
    }

    const double fn_prefactor = 1.54e-6 * local_field_v_per_m * local_field_v_per_m / work_function;
    const double fn_exponent =
        -6.83e9 * std::pow(work_function, 1.5) / std::max(1.0e6, local_field_v_per_m);
    return fn_prefactor * safeExp(fn_exponent);
}

double evaluateSurfaceThermionicEmissionCurrentDensity(const MaterialProperty& material,
                                                       double surface_temperature_k,
                                                       double surface_potential_v)
{
    const double work_function = std::max(0.1, material.getWorkFunctionEv());
    const double temperature_k = std::max(1.0, surface_temperature_k);
    const double thermal_energy_ev = std::max(1.0e-6, kBoltzmannEvPerK * temperature_k);
    return kRichardsonConstant * temperature_k * temperature_k *
           safeExp(-work_function / thermal_energy_ev) *
           evaluateLegacyEmissionEscapeProbability(material, surface_potential_v);
}

double evaluateSurfaceEffectiveConductivitySPerM(const MaterialProperty& material,
                                                 double base_conductivity_s_per_m,
                                                 double radiation_induced_conductivity_s_per_m,
                                                 double surface_potential_v,
                                                 double reference_potential_v,
                                                 double normal_electric_field_v_per_m,
                                                 double local_charge_density_c_per_m3,
                                                 double dielectric_thickness_m)
{
    const double base_conductivity = std::max(0.0, base_conductivity_s_per_m);
    const double conductivity_scale =
        std::clamp(material.deriveSurfaceConductivityScaleFactor(), 0.0, 1.0e12);
    const double thickness_m = std::max(1.0e-8, dielectric_thickness_m);
    const double electric_field_v_per_m =
        std::max(std::abs(normal_electric_field_v_per_m),
                 std::abs(surface_potential_v - reference_potential_v) / thickness_m);
    const double poole_frenkel_beta =
        std::max(0.0, material.getScalarProperty("poole_frenkel_beta", 0.0));
    const double max_enhancement =
        std::max(1.0, material.getScalarProperty("max_field_enhancement_factor", 1.0e8));
    const double field_enhancement = std::min(
        max_enhancement, safeExp(poole_frenkel_beta * std::sqrt(std::max(0.0, electric_field_v_per_m))));
    const double space_charge_ratio =
        std::abs(local_charge_density_c_per_m3) * thickness_m /
        std::max(1.0e-12, kEpsilon0 * std::max(1.0, electric_field_v_per_m));
    const double space_charge_enhancement = 1.0 + std::min(5.0, space_charge_ratio);
    return base_conductivity * conductivity_scale * field_enhancement * space_charge_enhancement +
           std::max(0.0, radiation_induced_conductivity_s_per_m);
}

LegacyEmissionIntegralBundle evaluateLegacyEmissionIntegralBundle(
    const MaterialProperty& material, LegacySecondaryYieldModel model,
    double electron_temperature_ev, double ion_temperature_ev, double surface_potential_v,
    double sims_exponent_n)
{
    LegacyEmissionIntegralBundle bundle;
    bundle.secondary_integral = integrateMaxwellianYield(
        electron_temperature_ev, [&](double energy_ev) {
            const double shifted =
                surface_potential_v > 0.0 ? energy_ev + surface_potential_v : energy_ev;
            return evaluateLegacySecondaryElectronYield(material, model, shifted,
                                                        sims_exponent_n);
        });
    bundle.backscatter_integral = integrateMaxwellianYield(
        electron_temperature_ev, [&](double energy_ev) {
            const double shifted =
                surface_potential_v > 0.0 ? energy_ev + surface_potential_v : energy_ev;
            return evaluateLegacyBackscatterYield(material, shifted);
        });
    bundle.ion_secondary_integral = integrateMaxwellianYield(
        ion_temperature_ev, [&](double energy_ev) {
            const double shifted =
                surface_potential_v <= 0.0 ? energy_ev - surface_potential_v : energy_ev;
            return evaluateLegacyIonSecondaryElectronYield(material, shifted);
        });
    bundle.emission_escape_probability =
        evaluateLegacyEmissionEscapeProbability(material, surface_potential_v);
    return bundle;
}

LegacyEmissionResponse evaluateLegacyEmissionResponse(
    const MaterialProperty& material, const LegacyEmissionResponseInputs& inputs)
{
    LegacyEmissionResponse response;
    response.integrals = evaluateLegacyEmissionIntegralBundle(
        material, inputs.model, inputs.electron_temperature_ev, inputs.ion_temperature_ev,
        inputs.surface_potential_v, inputs.sims_exponent_n);
    response.secondary_emission_a_per_m2 =
        std::abs(inputs.electron_collection_current_a_per_m2) *
        response.integrals.secondary_integral *
        response.integrals.emission_escape_probability;
    response.backscatter_emission_a_per_m2 =
        std::abs(inputs.electron_collection_current_a_per_m2) *
        response.integrals.backscatter_integral *
        response.integrals.emission_escape_probability;
    response.ion_secondary_emission_a_per_m2 =
        std::abs(inputs.ion_collection_current_a_per_m2) *
        response.integrals.ion_secondary_integral *
        response.integrals.emission_escape_probability;
    return response;
}

LegacyPopulationResponse evaluateLegacyPopulationResponse(
    const MaterialProperty& material, const LegacyPopulationResponseInputs& inputs)
{
    LegacyPopulationResponse response;
    const double mass_kg = std::max(1.0e-6, inputs.mass_amu) * kAtomicMassUnitKg;
    const double base =
        legacyThermalCurrentDensity(inputs.density_m3, inputs.temperature_ev, mass_kg);

    if (inputs.particle == LegacyCollectionParticle::Electron)
    {
        const double factor =
            inputs.surface_potential_v > 0.0
                ? 1.0
                : std::exp(inputs.surface_potential_v /
                           std::max(1.0e-6, inputs.electron_characteristic_energy_ev));
        response.collection_current_a_per_m2 =
            -base * std::max(0.0, inputs.electron_collection_coefficient) * factor;
        response.emission = evaluateLegacyEmissionResponse(
            material, LegacyEmissionResponseInputs{
                          inputs.model,
                          inputs.temperature_ev,
                          inputs.temperature_ev,
                          inputs.surface_potential_v,
                          response.collection_current_a_per_m2,
                          0.0,
                          inputs.sims_exponent_n,
                      });
    }
    else
    {
        const double factor =
            inputs.surface_potential_v > 0.0
                ? std::exp(-inputs.surface_potential_v /
                           std::max(1.0e-6, inputs.ion_characteristic_energy_ev))
                : 1.0;
        response.collection_current_a_per_m2 =
            base * std::max(0.0, inputs.ion_collection_coefficient) * factor;
        response.emission = evaluateLegacyEmissionResponse(
            material, LegacyEmissionResponseInputs{
                          inputs.model,
                          inputs.temperature_ev,
                          inputs.temperature_ev,
                          inputs.surface_potential_v,
                          0.0,
                          response.collection_current_a_per_m2,
                          inputs.sims_exponent_n,
                      });
    }

    return response;
}

} // namespace Material
} // namespace SCDAT
