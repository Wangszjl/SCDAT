#include "../include/SurfaceBarrierModels.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

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

double scalarWithAliases(const Material::MaterialProperty& material,
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

struct SurfaceFieldEmissionConfig
{
    double work_function_ev = 4.5;
    double fn_decay_coefficient_v_per_m_ev_1p5 = 5.0e5;
    double threshold_field_v_per_m = 1.0e7;
    double effective_sheath_floor_m = 1.0e-6;
};

SurfaceFieldEmissionConfig resolveSurfaceFieldEmissionConfig(
    const Material::MaterialProperty& material)
{
    SurfaceFieldEmissionConfig params;
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

std::string propertyScopeFromKey(const char* property_key)
{
    if (property_key == nullptr)
    {
        return {};
    }

    const std::string key(property_key);
    constexpr const char* kSuffix = "_barrier_scaler_kind";
    if (key.size() <= std::char_traits<char>::length(kSuffix))
    {
        return {};
    }
    if (key.rfind(kSuffix) != key.size() - std::char_traits<char>::length(kSuffix))
    {
        return {};
    }
    return key.substr(0, key.size() - std::char_traits<char>::length(kSuffix));
}

std::string scopedPropertyName(const std::string& scope, const char* suffix)
{
    if (suffix == nullptr || suffix[0] == '\0')
    {
        return scope;
    }
    return scope.empty() ? std::string(suffix) : scope + "_" + suffix;
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
    case 5:
        return SurfaceBarrierScalerKind::AutomaticBarrierCurrent;
    case 6:
        return SurfaceBarrierScalerKind::MultipleCurrent;
    case 7:
        return SurfaceBarrierScalerKind::GlobalTempCurrent;
    case 8:
        return SurfaceBarrierScalerKind::SmoothedGlobalTempCurrent;
    case 9:
        return SurfaceBarrierScalerKind::CurrentVariation;
    case 10:
        return SurfaceBarrierScalerKind::LocalVariation;
    default:
        return fallback_kind;
    }
}

const char* familyForSimpleKind(SurfaceBarrierScalerKind kind)
{
    switch (kind)
    {
    case SurfaceBarrierScalerKind::SecondaryRecollection:
        return "spis_secondary_recollection_scaler_v1";
    case SurfaceBarrierScalerKind::VariableBarrier:
        return "spis_variable_barrier_scaler_v1";
    case SurfaceBarrierScalerKind::OmlCurrent:
        return "spis_oml_current_scaler_v1";
    case SurfaceBarrierScalerKind::LteOmlCurrent:
        return "spis_lte_oml_current_scaler_v1";
    case SurfaceBarrierScalerKind::FowlerNordheimCurrent:
        return "spis_fowler_nordheim_current_scaler_v1";
    case SurfaceBarrierScalerKind::AutomaticBarrierCurrent:
        return "spis_automatic_barrier_current_scaler_v1";
    case SurfaceBarrierScalerKind::MultipleCurrent:
        return "spis_multiple_current_scaler_v1";
    case SurfaceBarrierScalerKind::GlobalTempCurrent:
        return "spis_global_temp_current_scaler_v1";
    case SurfaceBarrierScalerKind::SmoothedGlobalTempCurrent:
        return "spis_smoothed_global_temp_current_scaler_v1";
    case SurfaceBarrierScalerKind::CurrentVariation:
        return "spis_current_variation_scaler_v1";
    case SurfaceBarrierScalerKind::LocalVariation:
        return "spis_local_variation_scaler_v1";
    }

    return "spis_variable_barrier_scaler_v1";
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
    const auto fn_params = resolveSurfaceFieldEmissionConfig(material);

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
    const double sheath_floor_m = std::max(1.0e-9, fn_params.effective_sheath_floor_m);
    const double reference_field_v_per_m = std::abs(potential_gap_v) / sheath_floor_m;
    const double effective_field_v_per_m =
        std::max(normal_field_v_per_m, reference_field_v_per_m);
    if (effective_field_v_per_m <= fn_params.threshold_field_v_per_m)
    {
        evaluation.scaling = 0.0;
        evaluation.escaped_current_a_per_m2 = 0.0;
        evaluation.recollected_current_a_per_m2 =
            std::max(0.0, finiteOrDefault(emitted_current_a_per_m2, 0.0));
        evaluation.didv_a_per_m2_per_v = 0.0;
        evaluation.valid = std::isfinite(evaluation.recollected_current_a_per_m2);
        return evaluation;
    }

    double field_exponent =
        -fn_params.fn_decay_coefficient_v_per_m_ev_1p5 *
        std::pow(fn_params.work_function_ev, 1.5) /
        std::max(1.0e6, effective_field_v_per_m);
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

SurfaceBarrierScalerKind sanitizeMultipleInnerKind(SurfaceBarrierScalerKind kind)
{
    switch (kind)
    {
    case SurfaceBarrierScalerKind::MultipleCurrent:
        return SurfaceBarrierScalerKind::VariableBarrier;
    case SurfaceBarrierScalerKind::AutomaticBarrierCurrent:
        return SurfaceBarrierScalerKind::AutomaticBarrierCurrent;
    default:
        return kind;
    }
}

SurfaceBarrierScalerKind resolveAutomaticBarrierInnerKind(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state)
{
    const double local_potential_v = finiteOrDefault(state.local_potential_v, 0.0);
    const double reference_potential_v = finiteOrDefault(state.reference_potential_v, 0.0);
    const double barrier_potential_v = finiteOrDefault(state.barrier_potential_v, 0.0);
    const double potential_gap_v = local_potential_v - reference_potential_v + barrier_potential_v;
    const double abs_field_v_per_m =
        std::abs(finiteOrDefault(state.normal_electric_field_v_per_m, 0.0));
    const double emission_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(state.emission_temperature_ev,
                                         material.getPhotoelectronTemperatureEv()));

    const double fn_field_threshold_v_per_m = std::max(
        1.0e3, finiteOrDefault(
                    material.getScalarProperty("surface_automatic_fn_field_threshold_v_per_m",
                                               7.5e6),
                    7.5e6));
    const double lte_temperature_threshold_ev = std::max(
        1.0e-3, finiteOrDefault(
                    material.getScalarProperty("surface_automatic_lte_temperature_threshold_ev",
                                               2.5),
                    2.5));
    const double oml_gap_threshold_v = std::max(
        0.0, finiteOrDefault(
                 material.getScalarProperty("surface_automatic_oml_gap_threshold_v", 1.0), 1.0));
    const double oml_field_threshold_v_per_m = std::max(
        0.0, finiteOrDefault(
                 material.getScalarProperty("surface_automatic_oml_field_threshold_v_per_m",
                                            2.0e5),
                 2.0e5));

    if (abs_field_v_per_m >= fn_field_threshold_v_per_m)
    {
        return SurfaceBarrierScalerKind::FowlerNordheimCurrent;
    }
    if (potential_gap_v <= 0.0)
    {
        return emission_temperature_ev >= lte_temperature_threshold_ev
                   ? SurfaceBarrierScalerKind::LteOmlCurrent
                   : SurfaceBarrierScalerKind::OmlCurrent;
    }
    if (potential_gap_v <= oml_gap_threshold_v)
    {
        return emission_temperature_ev >= lte_temperature_threshold_ev ||
                       abs_field_v_per_m <= oml_field_threshold_v_per_m
                   ? SurfaceBarrierScalerKind::LteOmlCurrent
                   : SurfaceBarrierScalerKind::OmlCurrent;
    }
    if (abs_field_v_per_m <= oml_field_threshold_v_per_m * 0.5)
    {
        return SurfaceBarrierScalerKind::OmlCurrent;
    }
    return SurfaceBarrierScalerKind::VariableBarrier;
}

SurfaceBarrierEvaluation evaluateSimpleScaler(SurfaceBarrierScalerKind kind,
                                              const Material::MaterialProperty& material,
                                              const SurfaceBarrierState& state,
                                              double emitted_current_a_per_m2)
{
    return selectSurfaceBarrierScaler(kind).evaluate(material, state, emitted_current_a_per_m2);
}

SurfaceBarrierScalerKind sanitizeVariationInnerKind(SurfaceBarrierScalerKind kind)
{
    switch (kind)
    {
    case SurfaceBarrierScalerKind::MultipleCurrent:
    case SurfaceBarrierScalerKind::CurrentVariation:
    case SurfaceBarrierScalerKind::LocalVariation:
        return SurfaceBarrierScalerKind::VariableBarrier;
    default:
        return kind;
    }
}

SurfaceBarrierState buildGlobalTempState(const Material::MaterialProperty& material,
                                         const SurfaceBarrierState& state,
                                         const std::string& property_scope)
{
    SurfaceBarrierState effective_state = state;
    const double global_temp_ev = std::max(
        1.0e-3,
        finiteOrDefault(
            material.getScalarProperty(scopedPropertyName(property_scope, "global_temp_ev"),
                                       material.getPhotoelectronTemperatureEv()),
            material.getPhotoelectronTemperatureEv()));
    const double temperature_scale = std::clamp(
        finiteOrDefault(
            material.getScalarProperty(scopedPropertyName(property_scope, "global_temp_scale"), 1.0),
            1.0),
        0.05, 20.0);
    effective_state.emission_temperature_ev = global_temp_ev * temperature_scale;
    return effective_state;
}

SurfaceBarrierState buildSmoothedGlobalTempState(const Material::MaterialProperty& material,
                                                 const SurfaceBarrierState& state,
                                                 const std::string& property_scope)
{
    const SurfaceBarrierState global_state = buildGlobalTempState(material, state, property_scope);
    SurfaceBarrierState effective_state = state;
    const double smooth_weight = std::clamp(
        finiteOrDefault(
            material.getScalarProperty(scopedPropertyName(property_scope, "smoothed_temp_weight"),
                                       0.1),
            0.1),
        0.0, 1.0);
    effective_state.emission_temperature_ev =
        std::max(1.0e-3, finiteOrDefault(state.emission_temperature_ev, global_state.emission_temperature_ev)) *
            (1.0 - smooth_weight) +
        global_state.emission_temperature_ev * smooth_weight;
    effective_state.normal_electric_field_v_per_m =
        finiteOrDefault(state.normal_electric_field_v_per_m, 0.0) *
        std::clamp(
            finiteOrDefault(
                material.getScalarProperty(scopedPropertyName(property_scope, "smoothed_field_scale"),
                                           1.0),
                1.0),
            0.0, 4.0);
    return effective_state;
}

SurfaceBarrierEvaluation evaluateGlobalTempBarrier(const Material::MaterialProperty& material,
                                                   const SurfaceBarrierState& state,
                                                   double emitted_current_a_per_m2,
                                                   const std::string& property_scope)
{
    auto effective_state = buildGlobalTempState(material, state, property_scope);
    const double field_weight = std::clamp(
        finiteOrDefault(
            material.getScalarProperty(scopedPropertyName(property_scope, "global_temp_field_weight"),
                                       0.0),
            0.0),
        0.0, 1.0e-5);
    return evaluateBarrier(material, effective_state, emitted_current_a_per_m2, field_weight);
}

SurfaceBarrierEvaluation evaluateSmoothedGlobalTempBarrier(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2, const std::string& property_scope)
{
    auto effective_state = buildSmoothedGlobalTempState(material, state, property_scope);
    const double field_weight = std::clamp(
        finiteOrDefault(material.getScalarProperty(
                            scopedPropertyName(property_scope, "smoothed_global_temp_field_weight"),
                            0.0),
                        0.0),
        0.0, 1.0e-5);
    return evaluateBarrier(material, effective_state, emitted_current_a_per_m2, field_weight);
}

SurfaceBarrierState buildVariationReferenceState(const Material::MaterialProperty& material,
                                                 const SurfaceBarrierState& state,
                                                 const std::string& property_scope,
                                                 bool local_reference)
{
    SurfaceBarrierState reference_state = state;
    if (!local_reference)
    {
        reference_state.local_potential_v = finiteOrDefault(
            material.getScalarProperty(
                scopedPropertyName(property_scope, "current_variation_global_potential_v"),
                state.reference_potential_v),
            state.reference_potential_v);
        reference_state.reference_potential_v = finiteOrDefault(
            material.getScalarProperty(
                scopedPropertyName(property_scope, "current_variation_global_reference_potential_v"),
                state.reference_potential_v),
            state.reference_potential_v);
        reference_state.barrier_potential_v = finiteOrDefault(
            material.getScalarProperty(
                scopedPropertyName(property_scope, "current_variation_global_barrier_potential_v"),
                state.barrier_potential_v),
            state.barrier_potential_v);
        reference_state.normal_electric_field_v_per_m = finiteOrDefault(
            material.getScalarProperty(
                scopedPropertyName(property_scope, "current_variation_global_field_v_per_m"),
                0.0),
            0.0);
        reference_state.emission_temperature_ev = std::max(
            1.0e-3,
            finiteOrDefault(
                material.getScalarProperty(
                    scopedPropertyName(property_scope, "current_variation_global_temperature_ev"),
                    state.emission_temperature_ev),
                state.emission_temperature_ev));
        return reference_state;
    }

    reference_state.local_potential_v += finiteOrDefault(
        material.getScalarProperty(
            scopedPropertyName(property_scope, "local_variation_potential_offset_v"), 0.0),
        0.0);
    reference_state.barrier_potential_v += finiteOrDefault(
        material.getScalarProperty(
            scopedPropertyName(property_scope, "local_variation_barrier_offset_v"), 0.0),
        0.0);
    reference_state.normal_electric_field_v_per_m *= std::clamp(
        finiteOrDefault(
            material.getScalarProperty(
                scopedPropertyName(property_scope, "local_variation_field_scale"), 1.0),
            1.0),
        0.0, 4.0);
    reference_state.emission_temperature_ev *= std::clamp(
        finiteOrDefault(
            material.getScalarProperty(
                scopedPropertyName(property_scope, "local_variation_temperature_scale"), 1.0),
            1.0),
        0.05, 20.0);
    return reference_state;
}

SurfaceBarrierEvaluation composeVariationBarrier(const Material::MaterialProperty& material,
                                                 const SurfaceBarrierState& state,
                                                 double emitted_current_a_per_m2,
                                                 const std::string& property_scope,
                                                 bool local_reference)
{
    const auto resolve_kind_property = [&](const char* suffix,
                                           SurfaceBarrierScalerKind fallback_kind) {
        return sanitizeVariationInnerKind(resolveSurfaceBarrierScalerKind(
            material, scopedPropertyName(property_scope, suffix).c_str(), fallback_kind));
    };

    const auto base_kind = resolve_kind_property("current_variation_base_scaler_kind",
                                                 SurfaceBarrierScalerKind::VariableBarrier);
    const auto reference_kind = resolve_kind_property("current_variation_reference_scaler_kind",
                                                      SurfaceBarrierScalerKind::OmlCurrent);
    const double gain = std::clamp(
        finiteOrDefault(material.getScalarProperty(
                            scopedPropertyName(property_scope, local_reference
                                                                   ? "local_variation_gain"
                                                                   : "current_variation_gain"),
                            1.0),
                        1.0),
        0.0, 4.0);

    const auto base_eval = evaluateSimpleScaler(base_kind, material, state, emitted_current_a_per_m2);
    const auto reference_eval = evaluateSimpleScaler(
        reference_kind, material, buildVariationReferenceState(material, state, property_scope,
                                                               local_reference),
        1.0);

    SurfaceBarrierEvaluation evaluation = base_eval;
    const double reference_scaling =
        std::clamp(finiteOrDefault(reference_eval.scaling, 0.0), 1.0e-9, 1.0);
    const double adjusted_scaling =
        std::clamp(base_eval.scaling * std::pow(reference_scaling, gain), 0.0, 1.0);
    const double sanitized_emitted_current =
        std::max(0.0, finiteOrDefault(emitted_current_a_per_m2, 0.0));

    evaluation.scaling = adjusted_scaling;
    evaluation.escaped_current_a_per_m2 = sanitized_emitted_current * adjusted_scaling;
    evaluation.recollected_current_a_per_m2 =
        sanitized_emitted_current - evaluation.escaped_current_a_per_m2;

    double log_didv = 0.0;
    if (base_eval.escaped_current_a_per_m2 > 1.0e-18)
    {
        log_didv += base_eval.didv_a_per_m2_per_v / base_eval.escaped_current_a_per_m2;
    }
    if (reference_eval.escaped_current_a_per_m2 > 1.0e-18)
    {
        log_didv += gain * reference_eval.didv_a_per_m2_per_v /
                    reference_eval.escaped_current_a_per_m2;
    }
    evaluation.didv_a_per_m2_per_v = evaluation.escaped_current_a_per_m2 * log_didv;
    evaluation.valid = base_eval.valid && reference_eval.valid &&
                       std::isfinite(evaluation.scaling) &&
                       std::isfinite(evaluation.escaped_current_a_per_m2) &&
                       std::isfinite(evaluation.recollected_current_a_per_m2) &&
                       std::isfinite(evaluation.didv_a_per_m2_per_v);
    return evaluation;
}

SurfaceBarrierEvaluation evaluateMultipleBarrier(const Material::MaterialProperty& material,
                                                 const SurfaceBarrierState& state,
                                                 double emitted_current_a_per_m2,
                                                 const std::string& property_scope)
{
    const auto resolve_inner_kind = [&](const char* suffix,
                                        SurfaceBarrierScalerKind fallback_kind) {
        const auto property_name = scopedPropertyName(property_scope, suffix);
        return sanitizeMultipleInnerKind(resolveSurfaceBarrierScalerKind(
            material, property_name.c_str(), fallback_kind));
    };
    const auto resolve_weight = [&](const char* suffix, double fallback_weight) {
        const auto property_name = scopedPropertyName(property_scope, suffix);
        return std::max(0.0, finiteOrDefault(material.getScalarProperty(property_name, fallback_weight),
                                             fallback_weight));
    };

    const SurfaceBarrierScalerKind primary_kind = resolve_inner_kind(
        "multiple_current_primary_scaler_kind", SurfaceBarrierScalerKind::VariableBarrier);
    const SurfaceBarrierScalerKind secondary_kind = resolve_inner_kind(
        "multiple_current_secondary_scaler_kind", SurfaceBarrierScalerKind::OmlCurrent);
    const SurfaceBarrierScalerKind tertiary_kind = resolve_inner_kind(
        "multiple_current_tertiary_scaler_kind", SurfaceBarrierScalerKind::LteOmlCurrent);

    const double primary_weight = resolve_weight("multiple_current_primary_weight", 0.5);
    const double secondary_weight = resolve_weight("multiple_current_secondary_weight", 0.3);
    const double tertiary_weight = resolve_weight("multiple_current_tertiary_weight", 0.2);
    const double total_weight =
        std::max(1.0e-12, primary_weight + secondary_weight + tertiary_weight);

    SurfaceBarrierEvaluation evaluation;
    evaluation.scaling = 0.0;
    evaluation.valid = false;
    const auto accumulate = [&](SurfaceBarrierScalerKind kind, double weight) {
        if (weight <= 0.0)
        {
            return;
        }

        const double normalized_weight = weight / total_weight;
        const auto component =
            kind == SurfaceBarrierScalerKind::AutomaticBarrierCurrent
                ? AutomaticBarrierCurrentScaler().evaluate(material, state, emitted_current_a_per_m2)
                : evaluateSimpleScaler(kind, material, state, emitted_current_a_per_m2);

        evaluation.scaling += component.scaling * normalized_weight;
        evaluation.escaped_current_a_per_m2 +=
            component.escaped_current_a_per_m2 * normalized_weight;
        evaluation.recollected_current_a_per_m2 +=
            component.recollected_current_a_per_m2 * normalized_weight;
        evaluation.didv_a_per_m2_per_v += component.didv_a_per_m2_per_v * normalized_weight;
        evaluation.valid = evaluation.valid || component.valid;
    };

    accumulate(primary_kind, primary_weight);
    accumulate(secondary_kind, secondary_weight);
    accumulate(tertiary_kind, tertiary_weight);

    evaluation.scaling = std::clamp(finiteOrDefault(evaluation.scaling, 0.0), 0.0, 1.0);
    evaluation.escaped_current_a_per_m2 =
        std::max(0.0, finiteOrDefault(evaluation.escaped_current_a_per_m2, 0.0));
    evaluation.recollected_current_a_per_m2 =
        std::max(0.0, finiteOrDefault(evaluation.recollected_current_a_per_m2, 0.0));
    evaluation.didv_a_per_m2_per_v = finiteOrDefault(evaluation.didv_a_per_m2_per_v, 0.0);

    const double sanitized_emitted_current =
        std::max(0.0, finiteOrDefault(emitted_current_a_per_m2, 0.0));
    const double reconstructed_current =
        evaluation.escaped_current_a_per_m2 + evaluation.recollected_current_a_per_m2;
    if (std::abs(reconstructed_current - sanitized_emitted_current) >
        std::max(1.0e-18, sanitized_emitted_current * 1.0e-9))
    {
        evaluation.recollected_current_a_per_m2 =
            std::max(0.0, sanitized_emitted_current - evaluation.escaped_current_a_per_m2);
    }

    evaluation.valid = evaluation.valid && std::isfinite(evaluation.scaling) &&
                       std::isfinite(evaluation.escaped_current_a_per_m2) &&
                       std::isfinite(evaluation.recollected_current_a_per_m2) &&
                       std::isfinite(evaluation.didv_a_per_m2_per_v);
    return evaluation;
}

SurfaceBarrierScalerDecision makeDecision(SurfaceBarrierScalerKind configured_kind,
                                          SurfaceBarrierScalerKind resolved_kind,
                                          std::string resolved_family)
{
    SurfaceBarrierScalerDecision decision;
    decision.configured_kind = configured_kind;
    decision.resolved_kind = resolved_kind;
    decision.configured_family = familyForSimpleKind(configured_kind);
    decision.resolved_family =
        resolved_family.empty() ? familyForSimpleKind(resolved_kind) : std::move(resolved_family);
    return decision;
}

std::string describeMultipleBarrierResolvedFamily(const Material::MaterialProperty& material,
                                                  const std::string& property_scope)
{
    const auto resolve_inner_kind = [&](const char* suffix,
                                        SurfaceBarrierScalerKind fallback_kind) {
        const auto property_name = scopedPropertyName(property_scope, suffix);
        return sanitizeMultipleInnerKind(resolveSurfaceBarrierScalerKind(
            material, property_name.c_str(), fallback_kind));
    };

    const auto primary_kind = resolve_inner_kind("multiple_current_primary_scaler_kind",
                                                 SurfaceBarrierScalerKind::VariableBarrier);
    const auto secondary_kind = resolve_inner_kind("multiple_current_secondary_scaler_kind",
                                                   SurfaceBarrierScalerKind::OmlCurrent);
    const auto tertiary_kind = resolve_inner_kind("multiple_current_tertiary_scaler_kind",
                                                  SurfaceBarrierScalerKind::LteOmlCurrent);

    std::ostringstream stream;
    stream << familyForSimpleKind(SurfaceBarrierScalerKind::MultipleCurrent) << "["
           << familyForSimpleKind(primary_kind) << ","
           << familyForSimpleKind(secondary_kind) << ","
           << familyForSimpleKind(tertiary_kind) << "]";
    return stream.str();
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

const char* AutomaticBarrierCurrentScaler::family() const
{
    return "spis_automatic_barrier_current_scaler_v1";
}

SurfaceBarrierEvaluation AutomaticBarrierCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    const auto chosen_kind = resolveAutomaticBarrierInnerKind(material, state);
    return evaluateSimpleScaler(chosen_kind, material, state, emitted_current_a_per_m2);
}

MultipleCurrentScaler::MultipleCurrentScaler(std::string property_scope)
    : property_scope_(std::move(property_scope))
{
}

const char* MultipleCurrentScaler::family() const
{
    return "spis_multiple_current_scaler_v1";
}

SurfaceBarrierEvaluation MultipleCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateMultipleBarrier(material, state, emitted_current_a_per_m2, property_scope_);
}

GlobalTempCurrentScaler::GlobalTempCurrentScaler(std::string property_scope)
    : property_scope_(std::move(property_scope))
{
}

const char* GlobalTempCurrentScaler::family() const
{
    return "spis_global_temp_current_scaler_v1";
}

SurfaceBarrierEvaluation GlobalTempCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateGlobalTempBarrier(material, state, emitted_current_a_per_m2, property_scope_);
}

SmoothedGlobalTempCurrentScaler::SmoothedGlobalTempCurrentScaler(std::string property_scope)
    : property_scope_(std::move(property_scope))
{
}

const char* SmoothedGlobalTempCurrentScaler::family() const
{
    return "spis_smoothed_global_temp_current_scaler_v1";
}

SurfaceBarrierEvaluation SmoothedGlobalTempCurrentScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return evaluateSmoothedGlobalTempBarrier(material, state, emitted_current_a_per_m2,
                                             property_scope_);
}

CurrentVariationScaler::CurrentVariationScaler(std::string property_scope)
    : property_scope_(std::move(property_scope))
{
}

const char* CurrentVariationScaler::family() const
{
    return "spis_current_variation_scaler_v1";
}

SurfaceBarrierEvaluation CurrentVariationScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return composeVariationBarrier(material, state, emitted_current_a_per_m2, property_scope_,
                                   false);
}

LocalVariationScaler::LocalVariationScaler(std::string property_scope)
    : property_scope_(std::move(property_scope))
{
}

const char* LocalVariationScaler::family() const
{
    return "spis_local_variation_scaler_v1";
}

SurfaceBarrierEvaluation LocalVariationScaler::evaluate(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2) const
{
    return composeVariationBarrier(material, state, emitted_current_a_per_m2, property_scope_,
                                   true);
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
    static AutomaticBarrierCurrentScaler automatic_barrier_current_scaler;
    static MultipleCurrentScaler multiple_current_scaler;
    static GlobalTempCurrentScaler global_temp_current_scaler;
    static SmoothedGlobalTempCurrentScaler smoothed_global_temp_current_scaler;
    static CurrentVariationScaler current_variation_scaler;
    static LocalVariationScaler local_variation_scaler;

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
    case SurfaceBarrierScalerKind::AutomaticBarrierCurrent:
        return automatic_barrier_current_scaler;
    case SurfaceBarrierScalerKind::MultipleCurrent:
        return multiple_current_scaler;
    case SurfaceBarrierScalerKind::GlobalTempCurrent:
        return global_temp_current_scaler;
    case SurfaceBarrierScalerKind::SmoothedGlobalTempCurrent:
        return smoothed_global_temp_current_scaler;
    case SurfaceBarrierScalerKind::CurrentVariation:
        return current_variation_scaler;
    case SurfaceBarrierScalerKind::LocalVariation:
        return local_variation_scaler;
    }

    return variable_barrier_scaler;
}

SurfaceBarrierScalerDecision resolveSurfaceBarrierScalerDecision(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    const char* property_key, SurfaceBarrierScalerKind fallback_kind)
{
    const auto configured_kind =
        resolveSurfaceBarrierScalerKind(material, property_key, fallback_kind);
    switch (configured_kind)
    {
    case SurfaceBarrierScalerKind::AutomaticBarrierCurrent:
        return makeDecision(configured_kind, resolveAutomaticBarrierInnerKind(material, state), {});
    case SurfaceBarrierScalerKind::MultipleCurrent:
        return makeDecision(configured_kind, configured_kind,
                            describeMultipleBarrierResolvedFamily(material,
                                                                  propertyScopeFromKey(property_key)));
    case SurfaceBarrierScalerKind::CurrentVariation:
    case SurfaceBarrierScalerKind::LocalVariation:
    case SurfaceBarrierScalerKind::GlobalTempCurrent:
    case SurfaceBarrierScalerKind::SmoothedGlobalTempCurrent:
        return makeDecision(configured_kind, configured_kind, {});
    default:
        return makeDecision(configured_kind, configured_kind, {});
    }
}

SurfaceBarrierEvaluation evaluateConfiguredSurfaceBarrierScaler(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2, const char* property_key,
    SurfaceBarrierScalerKind fallback_kind)
{
    const auto configured_kind =
        resolveSurfaceBarrierScalerKind(material, property_key, fallback_kind);
    switch (configured_kind)
    {
    case SurfaceBarrierScalerKind::AutomaticBarrierCurrent:
        return AutomaticBarrierCurrentScaler().evaluate(material, state, emitted_current_a_per_m2);
    case SurfaceBarrierScalerKind::MultipleCurrent:
        return MultipleCurrentScaler(propertyScopeFromKey(property_key))
            .evaluate(material, state, emitted_current_a_per_m2);
    case SurfaceBarrierScalerKind::GlobalTempCurrent:
        return GlobalTempCurrentScaler(propertyScopeFromKey(property_key))
            .evaluate(material, state, emitted_current_a_per_m2);
    case SurfaceBarrierScalerKind::SmoothedGlobalTempCurrent:
        return SmoothedGlobalTempCurrentScaler(propertyScopeFromKey(property_key))
            .evaluate(material, state, emitted_current_a_per_m2);
    case SurfaceBarrierScalerKind::CurrentVariation:
        return CurrentVariationScaler(propertyScopeFromKey(property_key))
            .evaluate(material, state, emitted_current_a_per_m2);
    case SurfaceBarrierScalerKind::LocalVariation:
        return LocalVariationScaler(propertyScopeFromKey(property_key))
            .evaluate(material, state, emitted_current_a_per_m2);
    default:
        return selectSurfaceBarrierScaler(configured_kind)
            .evaluate(material, state, emitted_current_a_per_m2);
    }
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
    SurfaceEmissionBarrierOutputs outputs;
    outputs.escaped_photo = evaluateConfiguredSurfaceBarrierScaler(
        material, state, inputs.photo_emission_a_per_m2, "surface_photo_barrier_scaler_kind",
        SurfaceBarrierScalerKind::VariableBarrier);
    outputs.escaped_secondary = evaluateConfiguredSurfaceBarrierScaler(
        material, state, inputs.secondary_emission_a_per_m2,
        "surface_secondary_barrier_scaler_kind",
        SurfaceBarrierScalerKind::SecondaryRecollection);
    outputs.escaped_ion_secondary = evaluateConfiguredSurfaceBarrierScaler(
        material, state, inputs.ion_secondary_emission_a_per_m2,
        "surface_ion_secondary_barrier_scaler_kind",
        SurfaceBarrierScalerKind::SecondaryRecollection);
    outputs.escaped_backscatter = evaluateConfiguredSurfaceBarrierScaler(
        material, state, inputs.backscatter_emission_a_per_m2,
        "surface_backscatter_barrier_scaler_kind",
        SurfaceBarrierScalerKind::SecondaryRecollection);

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
