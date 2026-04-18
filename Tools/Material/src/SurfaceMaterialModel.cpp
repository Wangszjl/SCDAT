#include "../include/SurfaceMaterialModel.h"

#include <algorithm>
#include <cmath>
#include <string>

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

std::string normalizedMaterialName(const MaterialProperty& material)
{
    std::string name = material.getName();
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name;
}

double scalarWithAliases(const MaterialProperty& material,
                         std::initializer_list<const char*> keys,
                         double fallback)
{
    for (const char* key : keys)
    {
        if (key == nullptr || key[0] == '\0')
        {
            continue;
        }
        const double value = material.getScalarProperty(key, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(value))
        {
            return value;
        }
    }
    return fallback;
}

double interpolateLinear(double x0, double x1, double y0, double y1, double x)
{
    const double span = x1 - x0;
    if (!std::isfinite(span) || std::abs(span) <= 1.0e-12)
    {
        return y0;
    }
    const double t = std::clamp((x - x0) / span, 0.0, 1.0);
    return y0 + (y1 - y0) * t;
}

bool resolveBackscatterGrid(const MaterialProperty& material,
                            std::vector<double>& energies_ev,
                            std::vector<double>& angles_deg,
                            std::vector<double>& values)
{
    const int energy_count = static_cast<int>(std::llround(
        material.getScalarProperty("backscatter_table_energy_count", 0.0)));
    const int angle_count = static_cast<int>(std::llround(
        material.getScalarProperty("backscatter_table_angle_count", 0.0)));
    if (energy_count < 2 || angle_count < 2)
    {
        return false;
    }

    energies_ev.resize(static_cast<std::size_t>(energy_count));
    angles_deg.resize(static_cast<std::size_t>(angle_count));
    values.resize(static_cast<std::size_t>(energy_count * angle_count));

    for (int i = 0; i < energy_count; ++i)
    {
        const auto value =
            material.getScalarProperty("backscatter_table_energy_" + std::to_string(i) + "_ev",
                                       std::numeric_limits<double>::quiet_NaN());
        if (!std::isfinite(value))
        {
            return false;
        }
        energies_ev[static_cast<std::size_t>(i)] = value;
    }
    for (int j = 0; j < angle_count; ++j)
    {
        const auto value =
            material.getScalarProperty("backscatter_table_angle_" + std::to_string(j) + "_deg",
                                       std::numeric_limits<double>::quiet_NaN());
        if (!std::isfinite(value))
        {
            return false;
        }
        angles_deg[static_cast<std::size_t>(j)] = value;
    }
    for (int i = 0; i < energy_count; ++i)
    {
        for (int j = 0; j < angle_count; ++j)
        {
            const auto value = material.getScalarProperty(
                "backscatter_table_value_" + std::to_string(i) + "_" + std::to_string(j),
                std::numeric_limits<double>::quiet_NaN());
            if (!std::isfinite(value))
            {
                return false;
            }
            values[static_cast<std::size_t>(i * angle_count + j)] = value;
        }
    }
    return true;
}

double evaluateTabulatedBackscatterYield(const MaterialProperty& material,
                                         double incident_energy_ev,
                                         double angle_rad)
{
    std::vector<double> energies_ev;
    std::vector<double> angles_deg;
    std::vector<double> values;
    if (!resolveBackscatterGrid(material, energies_ev, angles_deg, values))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double query_energy = std::clamp(incident_energy_ev, energies_ev.front(), energies_ev.back());
    const double query_angle_deg =
        std::clamp(clampAngle(angle_rad) * 180.0 / kPi, angles_deg.front(), angles_deg.back());

    auto energy_upper =
        std::lower_bound(energies_ev.begin(), energies_ev.end(), query_energy);
    auto angle_upper =
        std::lower_bound(angles_deg.begin(), angles_deg.end(), query_angle_deg);
    const std::size_t e1 = std::min<std::size_t>(energies_ev.size() - 1,
                                                 static_cast<std::size_t>(energy_upper - energies_ev.begin()));
    const std::size_t a1 = std::min<std::size_t>(angles_deg.size() - 1,
                                                 static_cast<std::size_t>(angle_upper - angles_deg.begin()));
    const std::size_t e0 = e1 == 0 ? 0 : e1 - 1;
    const std::size_t a0 = a1 == 0 ? 0 : a1 - 1;
    const std::size_t angle_count = angles_deg.size();

    const double q00 = values[e0 * angle_count + a0];
    const double q01 = values[e0 * angle_count + a1];
    const double q10 = values[e1 * angle_count + a0];
    const double q11 = values[e1 * angle_count + a1];
    const double energy_low_interp = interpolateLinear(
        angles_deg[a0], angles_deg[a1], q00, q01, query_angle_deg);
    const double energy_high_interp = interpolateLinear(
        angles_deg[a0], angles_deg[a1], q10, q11, query_angle_deg);
    return std::clamp(interpolateLinear(energies_ev[e0], energies_ev[e1], energy_low_interp,
                                        energy_high_interp, query_energy),
                      0.0, 0.95);
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
    const double tabulated = evaluateTabulatedBackscatterYield(material, energy_ev, angle_rad);
    if (std::isfinite(tabulated))
    {
        return tabulated;
    }

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

double erosionActivation(const ErosionParamSet& params, const SurfaceModelContext& context)
{
    if (params.yield_scale <= 0.0)
    {
        return 0.0;
    }

    const double energy_excess_ev =
        std::max(0.0, context.incident_energy_ev - params.threshold_energy_ev);
    const double energy_factor =
        std::clamp(energy_excess_ev / std::max(1.0, params.characteristic_energy_ev), 0.0, 1.0);
    const double angle_factor = 1.0 + 0.25 * std::sin(clampAngle(context.incident_angle_rad));
    return std::clamp(params.yield_scale * energy_factor * angle_factor, 0.0, 1.0);
}

void applyErosionResponse(const MaterialProperty& material, const SurfaceModelContext& context,
                         SurfaceComponentCurrents& currents)
{
    const auto params = resolveErosionParamSet(material);
    const double activation = erosionActivation(params, context);
    if (activation <= 0.0)
    {
        return;
    }

    currents.secondary_electron_emission_a_per_m2 *= 1.0 + activation * params.secondary_gain;
    currents.ion_secondary_electron_emission_a_per_m2 *=
        1.0 + activation * params.ion_secondary_gain;
    currents.backscatter_emission_a_per_m2 *= 1.0 + activation * params.backscatter_gain;
    currents.photoelectron_emission_a_per_m2 *=
        std::max(0.0, 1.0 - activation * params.photo_suppression);
    currents.induced_conduction_a_per_m2 *=
        std::max(0.0, 1.0 - activation * params.conductivity_suppression);
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

ErosionParamSet resolveErosionParamSet(const MaterialProperty& material)
{
    ErosionParamSet params;
    const std::string material_name = normalizedMaterialName(material);

    if (material_name.find("kapton") != std::string::npos)
    {
        params.yield_scale = 0.20;
        params.threshold_energy_ev = 60.0;
        params.characteristic_energy_ev = 220.0;
        params.secondary_gain = 0.45;
        params.ion_secondary_gain = 0.90;
        params.backscatter_gain = 0.30;
        params.photo_suppression = 0.30;
        params.conductivity_suppression = 0.45;
    }
    else if (material_name.find("ptfe") != std::string::npos ||
             material_name.find("teflon") != std::string::npos)
    {
        params.yield_scale = 0.16;
        params.threshold_energy_ev = 70.0;
        params.characteristic_energy_ev = 260.0;
        params.secondary_gain = 0.40;
        params.ion_secondary_gain = 0.75;
        params.backscatter_gain = 0.35;
        params.photo_suppression = 0.25;
        params.conductivity_suppression = 0.35;
    }
    else if (material_name.find("al") != std::string::npos)
    {
        params.yield_scale = 0.08;
        params.threshold_energy_ev = 120.0;
        params.characteristic_energy_ev = 400.0;
        params.secondary_gain = 0.20;
        params.ion_secondary_gain = 0.45;
        params.backscatter_gain = 0.18;
        params.photo_suppression = 0.10;
        params.conductivity_suppression = 0.10;
    }

    params.yield_scale = std::max(
        0.0, scalarWithAliases(material, {"erosion_yield_scale", "surface_erosion_yield_scale"},
                               params.yield_scale));
    params.threshold_energy_ev = std::max(
        0.0, scalarWithAliases(material, {"erosion_threshold_energy_ev"}, params.threshold_energy_ev));
    params.characteristic_energy_ev = std::max(
        1.0, scalarWithAliases(material, {"erosion_characteristic_energy_ev"},
                               params.characteristic_energy_ev));
    params.secondary_gain = std::max(
        0.0, scalarWithAliases(material, {"erosion_secondary_gain"}, params.secondary_gain));
    params.ion_secondary_gain = std::max(
        0.0, scalarWithAliases(material, {"erosion_ion_secondary_gain"},
                               params.ion_secondary_gain));
    params.backscatter_gain = std::max(
        0.0, scalarWithAliases(material, {"erosion_backscatter_gain"}, params.backscatter_gain));
    params.photo_suppression = std::clamp(
        scalarWithAliases(material, {"erosion_photo_suppression"}, params.photo_suppression), 0.0,
        1.0);
    params.conductivity_suppression = std::clamp(
        scalarWithAliases(material, {"erosion_conductivity_suppression"},
                          params.conductivity_suppression),
        0.0, 1.0);
    return params;
}

SurfaceFieldEmissionParameters resolveSurfaceFieldEmissionParameters(
    const MaterialProperty& material)
{
    SurfaceFieldEmissionParameters params;
    params.work_function_ev = std::max(0.1, scalarWithAliases(
                                                material, {"surface_work_function_ev"},
                                                material.getWorkFunctionEv() > 0.0
                                                    ? material.getWorkFunctionEv()
                                                    : 4.5));
    params.fn_decay_coefficient_v_per_m_ev_1p5 = std::max(
        1.0, scalarWithAliases(material,
                               {"surface_fn_decay_coefficient_v_per_m_ev_1p5",
                                "fn_decay_coefficient_v_per_m_ev_1p5"},
                               5.0e5));
    params.threshold_field_v_per_m = std::max(
        1.0e3, scalarWithAliases(material,
                                 {"surface_fn_threshold_field_v_per_m",
                                  "fn_threshold_field_v_per_m"},
                                 1.0e7));
    params.effective_sheath_floor_m = std::max(
        1.0e-9, scalarWithAliases(material,
                                  {"surface_fn_effective_sheath_floor_m",
                                   "fn_effective_sheath_floor_m"},
                                  1.0e-6));
    return params;
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
    const double tabulated =
        evaluateTabulatedBackscatterYield(material, incident_energy_ev, 0.0);
    if (std::isfinite(tabulated))
    {
        return tabulated;
    }

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

double evaluateDebyeLengthM(double electron_temperature_ev,
                            double electron_density_m3)
{
    const double density_m3 = std::max(0.0, electron_density_m3);
    if (electron_temperature_ev <= 0.0 || density_m3 <= 0.0)
    {
        return 0.0;
    }

    return std::sqrt(kEpsilon0 * electron_temperature_ev /
                     (density_m3 * kElementaryCharge));
}

double evaluateSimpleSheathFieldVPerM(double surface_potential_v,
                                      double reference_potential_v,
                                      double effective_sheath_length_m,
                                      double sheath_floor_m)
{
    const double sheath_length =
        std::max(std::max(1.0e-12, sheath_floor_m), effective_sheath_length_m);
    return std::abs(reference_potential_v - surface_potential_v) / sheath_length;
}

double evaluateOmlLikeCollectionFactor(double surface_potential_v,
                                       double characteristic_energy_ev,
                                       LegacyCollectionParticle particle)
{
    const double energy_ev = std::max(1.0e-6, characteristic_energy_ev);
    const double normalized_potential = surface_potential_v / energy_ev;
    if (particle == LegacyCollectionParticle::Electron)
    {
        return normalized_potential >= 0.0 ? safeExp(-normalized_potential)
                                           : 1.0 + std::abs(normalized_potential);
    }

    return normalized_potential <= 0.0 ? safeExp(normalized_potential)
                                        : 1.0 + normalized_potential;
}

double evaluateMaxwellianEnergyWeight(double energy_ev,
                                      double temperature_ev)
{
    if (energy_ev < 0.0 || temperature_ev <= 0.0)
    {
        return 0.0;
    }

    return std::sqrt(std::max(0.0, energy_ev)) *
           std::exp(-energy_ev / std::max(1.0e-6, temperature_ev));
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
        const double weight = evaluateMaxwellianEnergyWeight(energy_ev, temperature_ev);
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
    const auto params = resolveSurfaceFieldEmissionParameters(material);
    const double work_function = params.work_function_ev;
    const double reference_field_v_per_m =
        evaluateSimpleSheathFieldVPerM(surface_potential_v, reference_potential_v,
                                       effective_sheath_length_m,
                                       params.effective_sheath_floor_m);
    const double local_field_v_per_m =
        std::max(std::abs(normal_electric_field_v_per_m), reference_field_v_per_m);
    if (local_field_v_per_m <= params.threshold_field_v_per_m)
    {
        return 0.0;
    }

    const double fn_prefactor = 1.54e-6 * local_field_v_per_m * local_field_v_per_m / work_function;
    const double fn_exponent =
        -params.fn_decay_coefficient_v_per_m_ev_1p5 * std::pow(work_function, 1.5) /
        std::max(1.0e6, local_field_v_per_m);
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
