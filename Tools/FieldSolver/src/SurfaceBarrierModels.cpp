#include "../include/SurfaceBarrierModels.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace FieldSolver
{

namespace
{

double finiteOrDefault(double value, double fallback)
{
    return std::isfinite(value) ? value : fallback;
}

SurfaceBarrierScalerKind decodeSurfaceBarrierScalerKind(double raw_kind_value,
                                                        SurfaceBarrierScalerKind fallback_kind)
{
    if (!std::isfinite(raw_kind_value))
    {
        return fallback_kind;
    }

    const int code = static_cast<int>(std::llround(raw_kind_value));
    switch (code)
    {
    case 0:
        return SurfaceBarrierScalerKind::SecondaryRecollection;
    case 1:
        return SurfaceBarrierScalerKind::VariableBarrier;
    case 2:
        return SurfaceBarrierScalerKind::OmlCurrent;
    case 3:
        return SurfaceBarrierScalerKind::LteOmlCurrent;
    case 4:
        return SurfaceBarrierScalerKind::FowlerNordheimCurrent;
    default:
        return fallback_kind;
    }
}

SurfaceBarrierEvaluation evaluateBarrier(const Material::MaterialProperty& material,
                                         const SurfaceBarrierState& state,
                                         double emitted_current_a_per_m2,
                                         double field_weight)
{
    SurfaceBarrierEvaluation evaluation;
    const double material_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(material.getPhotoelectronTemperatureEv(), 2.0));
    const double requested_temperature_ev =
        state.emission_temperature_ev > 0.0 ? state.emission_temperature_ev : material_temperature_ev;
    const double effective_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(requested_temperature_ev, material_temperature_ev));
    const double local_potential_v = finiteOrDefault(state.local_potential_v, 0.0);
    const double reference_potential_v = finiteOrDefault(state.reference_potential_v, 0.0);
    const double barrier_potential_v = finiteOrDefault(state.barrier_potential_v, 0.0);
    const double potential_gap_v = local_potential_v - reference_potential_v + barrier_potential_v;
    const double normal_field_v_per_m =
        finiteOrDefault(state.normal_electric_field_v_per_m, 0.0);
    const double field_term =
        std::clamp(std::abs(normal_field_v_per_m) * field_weight, 0.0, 4.0);

    double exponent = -(potential_gap_v / effective_temperature_ev) + field_term;
    if (!std::isfinite(exponent))
    {
        exponent = 0.0;
    }
    exponent = std::clamp(exponent, -50.0, 50.0);

    double scaling = std::exp(exponent);
    if (!std::isfinite(scaling))
    {
        scaling = exponent >= 0.0 ? 1.0 : 0.0;
    }
    scaling = std::clamp(scaling, 0.0, 1.0);

    const double sanitized_emitted_current =
        std::max(0.0, finiteOrDefault(emitted_current_a_per_m2, 0.0));

    evaluation.scaling = scaling;
    evaluation.escaped_current_a_per_m2 = sanitized_emitted_current * scaling;
    evaluation.recollected_current_a_per_m2 =
        sanitized_emitted_current - evaluation.escaped_current_a_per_m2;
    evaluation.didv_a_per_m2_per_v =
        -evaluation.escaped_current_a_per_m2 / effective_temperature_ev;
    if (!std::isfinite(evaluation.didv_a_per_m2_per_v))
    {
        evaluation.didv_a_per_m2_per_v = 0.0;
    }
    evaluation.valid = std::isfinite(evaluation.scaling) &&
                       std::isfinite(evaluation.escaped_current_a_per_m2) &&
                       std::isfinite(evaluation.recollected_current_a_per_m2) &&
                       std::isfinite(evaluation.didv_a_per_m2_per_v);
    return evaluation;
}

SurfaceBarrierEvaluation evaluateFowlerNordheimBarrier(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2)
{
    SurfaceBarrierEvaluation evaluation;

    const double material_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(material.getPhotoelectronTemperatureEv(), 2.0));
    const double requested_temperature_ev =
        state.emission_temperature_ev > 0.0 ? state.emission_temperature_ev : material_temperature_ev;
    const double effective_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(requested_temperature_ev, material_temperature_ev));
    const double local_potential_v = finiteOrDefault(state.local_potential_v, 0.0);
    const double reference_potential_v = finiteOrDefault(state.reference_potential_v, 0.0);
    const double barrier_potential_v = finiteOrDefault(state.barrier_potential_v, 0.0);
    const double potential_gap_v = local_potential_v - reference_potential_v + barrier_potential_v;

    const double normal_field_v_per_m =
        std::abs(finiteOrDefault(state.normal_electric_field_v_per_m, 0.0));
    const double effective_field_v_per_m = std::max(1.0e3, normal_field_v_per_m);
    const double work_function_ev =
        std::max(1.0e-3, finiteOrDefault(material.getWorkFunctionEv(), 4.5));
    const double fn_decay = std::max(
        1.0,
        finiteOrDefault(
            material.getScalarProperty("surface_fn_decay_coefficient_v_per_m_ev_1p5", 5.0e5),
            5.0e5));

    double field_exponent =
        -fn_decay * std::pow(work_function_ev, 1.5) / effective_field_v_per_m;
    if (!std::isfinite(field_exponent))
    {
        field_exponent = -50.0;
    }
    field_exponent = std::clamp(field_exponent, -50.0, 0.0);
    const double field_scaling = std::exp(field_exponent);

    const double potential_scaling =
        potential_gap_v <= 0.0 ? 1.0 : std::exp(-potential_gap_v / effective_temperature_ev);
    double scaling = field_scaling * potential_scaling;
    if (!std::isfinite(scaling))
    {
        scaling = 0.0;
    }
    scaling = std::clamp(scaling, 0.0, 1.0);

    const double sanitized_emitted_current =
        std::max(0.0, finiteOrDefault(emitted_current_a_per_m2, 0.0));
    const double escaped_current = sanitized_emitted_current * scaling;
    const double recollected_current = sanitized_emitted_current - escaped_current;
    const double didv =
        potential_gap_v > 0.0 ? -escaped_current / effective_temperature_ev : 0.0;

    evaluation.scaling = scaling;
    evaluation.escaped_current_a_per_m2 = escaped_current;
    evaluation.recollected_current_a_per_m2 = recollected_current;
    evaluation.didv_a_per_m2_per_v = std::isfinite(didv) ? didv : 0.0;
    evaluation.valid = std::isfinite(evaluation.scaling) &&
                       std::isfinite(evaluation.escaped_current_a_per_m2) &&
                       std::isfinite(evaluation.recollected_current_a_per_m2) &&
                       std::isfinite(evaluation.didv_a_per_m2_per_v);
    return evaluation;
}

} // namespace

const char* SecondaryRecollectionScaler::family() const
{
    return "spis_secondary_recollection_scaler_v1";
}

SurfaceBarrierEvaluation SecondaryRecollectionScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateBarrier(material, state, emitted_current_a_per_m2, 1.0e-7);
}

const char* VariableBarrierScaler::family() const
{
    return "spis_variable_barrier_scaler_v1";
}

SurfaceBarrierEvaluation VariableBarrierScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateBarrier(material, state, emitted_current_a_per_m2, 2.5e-7);
}

const char* OmlCurrentScaler::family() const
{
    return "spis_oml_current_scaler_v1";
}

SurfaceBarrierEvaluation OmlCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateBarrier(material, state, emitted_current_a_per_m2, 0.0);
}

const char* LteOmlCurrentScaler::family() const
{
    return "spis_lte_oml_current_scaler_v1";
}

SurfaceBarrierEvaluation LteOmlCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    SurfaceBarrierState effective_state = state;
    const double temperature_scale = std::clamp(
        finiteOrDefault(material.getScalarProperty("surface_lte_oml_temperature_scale", 0.85),
                        0.85),
        0.1, 5.0);
    effective_state.emission_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(state.emission_temperature_ev, 2.0) * temperature_scale);
    const double field_weight = std::clamp(
        finiteOrDefault(material.getScalarProperty("surface_lte_oml_field_weight", 5.0e-8),
                        5.0e-8),
        0.0, 1.0e-5);
    return evaluateBarrier(material, effective_state, emitted_current_a_per_m2, field_weight);
}

const char* FowlerNordheimCurrentScaler::family() const
{
    return "spis_fowler_nordheim_current_scaler_v1";
}

SurfaceBarrierEvaluation FowlerNordheimCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateFowlerNordheimBarrier(material, state, emitted_current_a_per_m2);
}

SurfaceBarrierScalerKind resolveSurfaceBarrierScalerKind(
    const Material::MaterialProperty& material, const char* property_key,
    SurfaceBarrierScalerKind fallback_kind)
{
    if (property_key == nullptr || property_key[0] == '\0')
    {
        return fallback_kind;
    }
    const double raw_kind_value = material.getScalarProperty(
        property_key, static_cast<double>(static_cast<int>(fallback_kind)));
    return decodeSurfaceBarrierScalerKind(raw_kind_value, fallback_kind);
}

const BarrierCurrentScaler& selectSurfaceBarrierScaler(SurfaceBarrierScalerKind kind)
{
    static SecondaryRecollectionScaler secondary_recollection_scaler;
    static VariableBarrierScaler variable_barrier_scaler;
    static OmlCurrentScaler oml_current_scaler;
    static LteOmlCurrentScaler lte_oml_current_scaler;
    static FowlerNordheimCurrentScaler fowler_nordheim_current_scaler;

    switch (kind)
    {
    case SurfaceBarrierScalerKind::SecondaryRecollection:
        return secondary_recollection_scaler;
    case SurfaceBarrierScalerKind::VariableBarrier:
        return variable_barrier_scaler;
    case SurfaceBarrierScalerKind::OmlCurrent:
        return oml_current_scaler;
    case SurfaceBarrierScalerKind::LteOmlCurrent:
        return lte_oml_current_scaler;
    case SurfaceBarrierScalerKind::FowlerNordheimCurrent:
        return fowler_nordheim_current_scaler;
    }

    return variable_barrier_scaler;
}

SurfaceEmissionBarrierInputs makeSurfaceEmissionBarrierInputs(
    const SurfaceEmissionBarrierComponentInputs& components)
{
    SurfaceEmissionBarrierInputs inputs;
    inputs.electron_collection_a_per_m2 = components.electron_collection_a_per_m2;
    inputs.ion_collection_a_per_m2 = components.ion_collection_a_per_m2;
    inputs.secondary_emission_a_per_m2 = components.secondary_emission_a_per_m2;
    inputs.ion_secondary_emission_a_per_m2 = components.ion_secondary_emission_a_per_m2;
    inputs.backscatter_emission_a_per_m2 = components.backscatter_emission_a_per_m2;
    inputs.photo_emission_a_per_m2 = components.photo_emission_a_per_m2;
    inputs.photo_incidence_scale = components.photo_incidence_scale;
    inputs.user_secondary_scale = components.user_secondary_scale;
    inputs.user_ion_secondary_scale = components.user_ion_secondary_scale;
    inputs.user_backscatter_scale = components.user_backscatter_scale;
    inputs.user_photo_scale = components.user_photo_scale;
    inputs.fallback_secondary_yield = components.fallback_secondary_yield;
    inputs.fallback_ion_secondary_yield = components.fallback_ion_secondary_yield;
    inputs.fallback_backscatter_yield = components.fallback_backscatter_yield;
    return inputs;
}

SurfaceEmissionBarrierOutputs evaluateSurfaceEmissionBarrierBundle(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    const SurfaceEmissionBarrierComponentInputs& components)
{
    return evaluateSurfaceEmissionBarrierBundle(material, state,
                                                makeSurfaceEmissionBarrierInputs(components));
}

SurfaceEmissionBarrierOutputs evaluateSurfaceEmissionBarrierBundle(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    const SurfaceEmissionBarrierInputs& inputs)
{
    const auto& photo_barrier_scaler = selectSurfaceBarrierScaler(resolveSurfaceBarrierScalerKind(
        material, "surface_photo_barrier_scaler_kind",
        SurfaceBarrierScalerKind::VariableBarrier));
    const auto& secondary_barrier_scaler = selectSurfaceBarrierScaler(
        resolveSurfaceBarrierScalerKind(material, "surface_secondary_barrier_scaler_kind",
                                        SurfaceBarrierScalerKind::SecondaryRecollection));
    const auto& ion_secondary_barrier_scaler = selectSurfaceBarrierScaler(
        resolveSurfaceBarrierScalerKind(material, "surface_ion_secondary_barrier_scaler_kind",
                                        SurfaceBarrierScalerKind::SecondaryRecollection));
    const auto& backscatter_barrier_scaler = selectSurfaceBarrierScaler(
        resolveSurfaceBarrierScalerKind(material, "surface_backscatter_barrier_scaler_kind",
                                        SurfaceBarrierScalerKind::SecondaryRecollection));

    SurfaceEmissionBarrierOutputs outputs;
    outputs.escaped_photo =
        photo_barrier_scaler.evaluate(material, state, inputs.photo_emission_a_per_m2);
    outputs.escaped_secondary =
        secondary_barrier_scaler.evaluate(material, state, inputs.secondary_emission_a_per_m2);
    outputs.escaped_ion_secondary = ion_secondary_barrier_scaler.evaluate(
        material, state, inputs.ion_secondary_emission_a_per_m2);
    outputs.escaped_backscatter = backscatter_barrier_scaler.evaluate(
        material, state, inputs.backscatter_emission_a_per_m2);

    const double electron_collection_abs =
        std::abs(finiteOrDefault(inputs.electron_collection_a_per_m2, 0.0));
    const double ion_collection_abs =
        std::abs(finiteOrDefault(inputs.ion_collection_a_per_m2, 0.0));
    const double photo_emission_abs =
        std::abs(finiteOrDefault(inputs.photo_emission_a_per_m2, 0.0));

    const double user_secondary_scale =
        std::max(0.0, finiteOrDefault(inputs.user_secondary_scale, 0.0));
    const double user_ion_secondary_scale =
        std::max(0.0, finiteOrDefault(inputs.user_ion_secondary_scale, 0.0));
    const double user_backscatter_scale =
        std::max(0.0, finiteOrDefault(inputs.user_backscatter_scale, 0.0));
    const double user_photo_scale = std::max(0.0, finiteOrDefault(inputs.user_photo_scale, 0.0));

    const double fallback_secondary_yield =
        std::max(0.25, finiteOrDefault(inputs.fallback_secondary_yield, 0.25));
    const double fallback_ion_secondary_yield =
        std::max(0.25, finiteOrDefault(inputs.fallback_ion_secondary_yield, 0.25));
    const double fallback_backscatter_yield =
        std::max(0.0, finiteOrDefault(inputs.fallback_backscatter_yield, 0.0));
    const double photo_incidence_scale =
        std::clamp(finiteOrDefault(inputs.photo_incidence_scale, 1.0), 0.0, 2.0);

    const auto saturateScale = [](double value) {
        return std::isfinite(value) ? std::clamp(value, 0.0, 4.0) : 0.0;
    };

    const double secondary_ratio =
        electron_collection_abs > 1.0e-18
            ? outputs.escaped_secondary.escaped_current_a_per_m2 / electron_collection_abs
            : fallback_secondary_yield;
    const double ion_secondary_ratio =
        ion_collection_abs > 1.0e-18
            ? outputs.escaped_ion_secondary.escaped_current_a_per_m2 / ion_collection_abs
            : fallback_ion_secondary_yield;
    const double backscatter_ratio =
        electron_collection_abs > 1.0e-18
            ? outputs.escaped_backscatter.escaped_current_a_per_m2 / electron_collection_abs
            : fallback_backscatter_yield;
    const double photo_ratio =
        photo_emission_abs > 1.0e-18
            ? outputs.escaped_photo.escaped_current_a_per_m2 / photo_emission_abs
            : 1.0;

    outputs.secondary_emission_scale =
        saturateScale(user_secondary_scale * finiteOrDefault(secondary_ratio, 0.0));
    outputs.ion_secondary_emission_scale =
        saturateScale(user_ion_secondary_scale * finiteOrDefault(ion_secondary_ratio, 0.0));
    outputs.backscatter_scale =
        saturateScale(user_backscatter_scale * finiteOrDefault(backscatter_ratio, 0.0));
    outputs.photo_emission_scale = saturateScale(
        user_photo_scale * photo_incidence_scale * finiteOrDefault(photo_ratio, 0.0));
    outputs.valid = outputs.escaped_photo.valid && outputs.escaped_secondary.valid &&
                    outputs.escaped_ion_secondary.valid && outputs.escaped_backscatter.valid &&
                    std::isfinite(outputs.secondary_emission_scale) &&
                    std::isfinite(outputs.ion_secondary_emission_scale) &&
                    std::isfinite(outputs.backscatter_scale) &&
                    std::isfinite(outputs.photo_emission_scale);
    return outputs;
}

double evaluateLegacyThermalIonCollectionNaPerM2(
    const std::vector<Particle::LegacyPopulation>& populations,
    const Material::MaterialProperty& material,
    double surface_potential_v,
    bool include_potential_barrier)
{
    double current_na_per_m2 = 0.0;
    VariableBarrierScaler scaler;
    for (const auto& population : populations)
    {
        if (population.density_m3 <= 0.0)
        {
            continue;
        }

        const double ni_cm3 = population.density_m3 * 1.0e-6;
        const double temperature = std::max(1.0e-6, population.temperature_ev);
        double term =
            0.6255 * ni_cm3 * std::sqrt(temperature / std::max(1.0, population.mass_amu));
        if (include_potential_barrier && surface_potential_v > 0.0)
        {
            SurfaceBarrierState state;
            state.local_potential_v = surface_potential_v;
            state.reference_potential_v = 0.0;
            state.barrier_potential_v = 0.0;
            state.normal_electric_field_v_per_m = 0.0;
            state.emission_temperature_ev = temperature;

            const auto evaluation = scaler.evaluate(material, state, 1.0);
            term *= evaluation.valid ? std::clamp(evaluation.scaling, 0.0, 1.0) : 0.0;
        }
        current_na_per_m2 += term;
    }
    return current_na_per_m2;
}

} // namespace FieldSolver
} // namespace SCDAT
