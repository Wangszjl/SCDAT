#include "../include/SurfaceBarrierModels.h"

#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

TEST(SurfaceBarrierModelsTest, SecondaryRecollectionScalerProducesBoundedScaling)
{
    SCDAT::Material::MaterialProperty material(1, SCDAT::Mesh::MaterialType::DIELECTRIC, "kapton");
    material.setPhotoelectronTemperatureEv(2.0);

    SCDAT::FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = -5.0;
    state.reference_potential_v = 0.0;
    state.emission_temperature_ev = 2.0;

    SCDAT::FieldSolver::SecondaryRecollectionScaler scaler;
    const auto evaluation = scaler.evaluate(material, state, 1.0e-6);

    EXPECT_STREQ(scaler.family(), "spis_secondary_recollection_scaler_v1");
    EXPECT_TRUE(evaluation.valid);
    EXPECT_GE(evaluation.scaling, 0.0);
    EXPECT_LE(evaluation.scaling, 1.0);
    EXPECT_NEAR(evaluation.escaped_current_a_per_m2 + evaluation.recollected_current_a_per_m2,
                1.0e-6, 1.0e-12);
}

TEST(SurfaceBarrierModelsTest, VariableBarrierScalerChangesDerivativeWithPotential)
{
    SCDAT::Material::MaterialProperty material(2, SCDAT::Mesh::MaterialType::DIELECTRIC, "teflon");
    material.setPhotoelectronTemperatureEv(1.5);

    SCDAT::FieldSolver::SurfaceBarrierState low_state;
    low_state.local_potential_v = -4.0;
    low_state.reference_potential_v = 0.0;
    low_state.normal_electric_field_v_per_m = 5.0e4;
    low_state.emission_temperature_ev = 1.5;

    auto high_state = low_state;
    high_state.local_potential_v = 4.0;

    SCDAT::FieldSolver::VariableBarrierScaler scaler;
    const auto low = scaler.evaluate(material, low_state, 2.0e-6);
    const auto high = scaler.evaluate(material, high_state, 2.0e-6);

    EXPECT_STREQ(scaler.family(), "spis_variable_barrier_scaler_v1");
    EXPECT_TRUE(low.valid);
    EXPECT_TRUE(high.valid);
    EXPECT_GT(low.escaped_current_a_per_m2, high.escaped_current_a_per_m2);
    EXPECT_TRUE(std::isfinite(low.didv_a_per_m2_per_v));
    EXPECT_TRUE(std::isfinite(high.didv_a_per_m2_per_v));
}

TEST(SurfaceBarrierModelsTest, ScalerKindRoutingResolvesExpectedFamilies)
{
    SCDAT::Material::MaterialProperty material(11, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "scaler-routing");

    const auto default_kind = SCDAT::FieldSolver::resolveSurfaceBarrierScalerKind(
        material, "surface_photo_barrier_scaler_kind",
        SCDAT::FieldSolver::SurfaceBarrierScalerKind::VariableBarrier);
    EXPECT_EQ(default_kind, SCDAT::FieldSolver::SurfaceBarrierScalerKind::VariableBarrier);
    EXPECT_STREQ(SCDAT::FieldSolver::selectSurfaceBarrierScaler(default_kind).family(),
                 "spis_variable_barrier_scaler_v1");

    material.setScalarProperty("surface_photo_barrier_scaler_kind", 2.0);
    const auto oml_kind = SCDAT::FieldSolver::resolveSurfaceBarrierScalerKind(
        material, "surface_photo_barrier_scaler_kind",
        SCDAT::FieldSolver::SurfaceBarrierScalerKind::VariableBarrier);
    EXPECT_EQ(oml_kind, SCDAT::FieldSolver::SurfaceBarrierScalerKind::OmlCurrent);
    EXPECT_STREQ(SCDAT::FieldSolver::selectSurfaceBarrierScaler(oml_kind).family(),
                 "spis_oml_current_scaler_v1");

    material.setScalarProperty("surface_photo_barrier_scaler_kind", 4.0);
    const auto fn_kind = SCDAT::FieldSolver::resolveSurfaceBarrierScalerKind(
        material, "surface_photo_barrier_scaler_kind",
        SCDAT::FieldSolver::SurfaceBarrierScalerKind::VariableBarrier);
    EXPECT_EQ(fn_kind,
              SCDAT::FieldSolver::SurfaceBarrierScalerKind::FowlerNordheimCurrent);
    EXPECT_STREQ(SCDAT::FieldSolver::selectSurfaceBarrierScaler(fn_kind).family(),
                 "spis_fowler_nordheim_current_scaler_v1");

    material.setScalarProperty("surface_photo_barrier_scaler_kind", 99.0);
    const auto fallback_kind = SCDAT::FieldSolver::resolveSurfaceBarrierScalerKind(
        material, "surface_photo_barrier_scaler_kind",
        SCDAT::FieldSolver::SurfaceBarrierScalerKind::SecondaryRecollection);
    EXPECT_EQ(fallback_kind,
              SCDAT::FieldSolver::SurfaceBarrierScalerKind::SecondaryRecollection);
}

TEST(SurfaceBarrierModelsTest, EmissionBarrierBundleRoutesConfiguredScalerFamilies)
{
    SCDAT::Material::MaterialProperty material(12, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "scaler-routing-bundle");
    material.setPhotoelectronTemperatureEv(2.0);
    material.setWorkFunctionEv(4.5);

    SCDAT::FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = 5.0;
    state.reference_potential_v = 0.0;
    state.barrier_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 2.0e5;
    state.emission_temperature_ev = 2.0;

    SCDAT::FieldSolver::SurfaceEmissionBarrierInputs inputs;
    inputs.electron_collection_a_per_m2 = -3.0e-6;
    inputs.ion_collection_a_per_m2 = 2.0e-6;
    inputs.secondary_emission_a_per_m2 = 8.0e-7;
    inputs.ion_secondary_emission_a_per_m2 = 3.0e-7;
    inputs.backscatter_emission_a_per_m2 = 2.0e-7;
    inputs.photo_emission_a_per_m2 = 1.0e-7;
    inputs.photo_incidence_scale = 1.0;
    inputs.user_secondary_scale = 1.0;
    inputs.user_ion_secondary_scale = 1.0;
    inputs.user_backscatter_scale = 1.0;
    inputs.user_photo_scale = 1.0;
    inputs.fallback_secondary_yield = 0.5;
    inputs.fallback_ion_secondary_yield = 0.3;
    inputs.fallback_backscatter_yield = 0.2;

    const auto default_outputs =
        SCDAT::FieldSolver::evaluateSurfaceEmissionBarrierBundle(material, state, inputs);
    ASSERT_TRUE(default_outputs.valid);

    material.setScalarProperty("surface_photo_barrier_scaler_kind", 4.0);
    material.setScalarProperty("surface_secondary_barrier_scaler_kind", 2.0);
    material.setScalarProperty("surface_ion_secondary_barrier_scaler_kind", 3.0);
    material.setScalarProperty("surface_backscatter_barrier_scaler_kind", 1.0);

    const auto routed_outputs =
        SCDAT::FieldSolver::evaluateSurfaceEmissionBarrierBundle(material, state, inputs);
    ASSERT_TRUE(routed_outputs.valid);

    EXPECT_LT(routed_outputs.photo_emission_scale, default_outputs.photo_emission_scale);
    EXPECT_NEAR(routed_outputs.escaped_secondary.scaling,
                default_outputs.escaped_secondary.scaling, 5.0e-2);
    EXPECT_NE(routed_outputs.escaped_ion_secondary.scaling,
              default_outputs.escaped_ion_secondary.scaling);
}

TEST(SurfaceBarrierModelsTest, EmissionBarrierBundleProducesResolvedComponentScales)
{
    SCDAT::Material::MaterialProperty material(3, SCDAT::Mesh::MaterialType::DIELECTRIC, "bundle");
    material.setPhotoelectronTemperatureEv(2.0);

    SCDAT::FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = -3.0;
    state.reference_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 2.0e4;
    state.emission_temperature_ev = 2.0;

    SCDAT::FieldSolver::SurfaceEmissionBarrierInputs inputs;
    inputs.electron_collection_a_per_m2 = -3.0e-6;
    inputs.ion_collection_a_per_m2 = 2.0e-6;
    inputs.secondary_emission_a_per_m2 = 8.0e-7;
    inputs.ion_secondary_emission_a_per_m2 = 3.0e-7;
    inputs.backscatter_emission_a_per_m2 = 2.0e-7;
    inputs.photo_emission_a_per_m2 = 1.0e-7;
    inputs.photo_incidence_scale = 0.8;
    inputs.user_secondary_scale = 1.2;
    inputs.user_ion_secondary_scale = 0.9;
    inputs.user_backscatter_scale = 0.7;
    inputs.user_photo_scale = 1.1;
    inputs.fallback_secondary_yield = 0.5;
    inputs.fallback_ion_secondary_yield = 0.3;
    inputs.fallback_backscatter_yield = 0.2;

    const auto outputs =
        SCDAT::FieldSolver::evaluateSurfaceEmissionBarrierBundle(material, state, inputs);

    EXPECT_TRUE(outputs.valid);
    EXPECT_GT(outputs.secondary_emission_scale, 0.0);
    EXPECT_GT(outputs.ion_secondary_emission_scale, 0.0);
    EXPECT_GT(outputs.backscatter_scale, 0.0);
    EXPECT_GT(outputs.photo_emission_scale, 0.0);
    EXPECT_TRUE(outputs.escaped_photo.valid);
}

TEST(SurfaceBarrierModelsTest, BuildsBarrierInputsFromComponentBundle)
{
    SCDAT::FieldSolver::SurfaceEmissionBarrierComponentInputs components;
    components.electron_collection_a_per_m2 = -2.0e-6;
    components.ion_collection_a_per_m2 = 1.0e-6;
    components.secondary_emission_a_per_m2 = 7.0e-7;
    components.ion_secondary_emission_a_per_m2 = 2.0e-7;
    components.backscatter_emission_a_per_m2 = 1.0e-7;
    components.photo_emission_a_per_m2 = 5.0e-8;
    components.photo_incidence_scale = 0.6;
    components.user_secondary_scale = 1.2;
    components.user_ion_secondary_scale = 0.9;
    components.user_backscatter_scale = 0.7;
    components.user_photo_scale = 1.1;
    components.fallback_secondary_yield = 0.5;
    components.fallback_ion_secondary_yield = 0.3;
    components.fallback_backscatter_yield = 0.2;

    const auto inputs = SCDAT::FieldSolver::makeSurfaceEmissionBarrierInputs(components);
    EXPECT_DOUBLE_EQ(inputs.electron_collection_a_per_m2, components.electron_collection_a_per_m2);
    EXPECT_DOUBLE_EQ(inputs.photo_incidence_scale, components.photo_incidence_scale);
    EXPECT_DOUBLE_EQ(inputs.user_photo_scale, components.user_photo_scale);
    EXPECT_DOUBLE_EQ(inputs.fallback_backscatter_yield, components.fallback_backscatter_yield);
}

TEST(SurfaceBarrierModelsTest, ComponentBundleOverloadMatchesLegacyInputPath)
{
    SCDAT::Material::MaterialProperty material(4, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "kapton-overload");
    material.setPhotoelectronTemperatureEv(1.8);

    SCDAT::FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = -2.5;
    state.reference_potential_v = 0.5;
    state.barrier_potential_v = 0.2;
    state.normal_electric_field_v_per_m = 1.5e4;
    state.emission_temperature_ev = 1.8;

    SCDAT::FieldSolver::SurfaceEmissionBarrierComponentInputs components;
    components.electron_collection_a_per_m2 = -2.5e-6;
    components.ion_collection_a_per_m2 = 1.4e-6;
    components.secondary_emission_a_per_m2 = 8.0e-7;
    components.ion_secondary_emission_a_per_m2 = 3.5e-7;
    components.backscatter_emission_a_per_m2 = 1.5e-7;
    components.photo_emission_a_per_m2 = 9.0e-8;
    components.photo_incidence_scale = 0.75;
    components.user_secondary_scale = 1.15;
    components.user_ion_secondary_scale = 0.85;
    components.user_backscatter_scale = 0.7;
    components.user_photo_scale = 1.05;
    components.fallback_secondary_yield = 0.55;
    components.fallback_ion_secondary_yield = 0.35;
    components.fallback_backscatter_yield = 0.25;

    const auto via_components =
        SCDAT::FieldSolver::evaluateSurfaceEmissionBarrierBundle(material, state, components);
    const auto via_inputs = SCDAT::FieldSolver::evaluateSurfaceEmissionBarrierBundle(
        material, state, SCDAT::FieldSolver::makeSurfaceEmissionBarrierInputs(components));

    EXPECT_TRUE(via_components.valid);
    EXPECT_TRUE(via_inputs.valid);
    EXPECT_NEAR(via_components.secondary_emission_scale, via_inputs.secondary_emission_scale,
                1.0e-12);
    EXPECT_NEAR(via_components.ion_secondary_emission_scale,
                via_inputs.ion_secondary_emission_scale, 1.0e-12);
    EXPECT_NEAR(via_components.backscatter_scale, via_inputs.backscatter_scale, 1.0e-12);
    EXPECT_NEAR(via_components.photo_emission_scale, via_inputs.photo_emission_scale, 1.0e-12);
    EXPECT_NEAR(via_components.escaped_photo.escaped_current_a_per_m2,
                via_inputs.escaped_photo.escaped_current_a_per_m2, 1.0e-12);
}

TEST(SurfaceBarrierModelsTest,
     VariableBarrierScalerMatchesLegacyEscapeLawWhenBarrierStateIsNeutral)
{
    SCDAT::Material::MaterialProperty material(5, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "escape-law-parity");

    SCDAT::FieldSolver::SurfaceBarrierState state;
    state.reference_potential_v = 0.0;
    state.barrier_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 0.0;
    state.emission_temperature_ev = 2.5;

    SCDAT::FieldSolver::VariableBarrierScaler scaler;
    constexpr double kIncomingCurrentApm2 = 4.0e-6;

    const auto assert_escape_law = [&](double local_potential_v) {
        state.local_potential_v = local_potential_v;
        const auto evaluation = scaler.evaluate(material, state, kIncomingCurrentApm2);
        const double expected_scaling =
            local_potential_v <= 0.0
                ? 1.0
                : std::exp(-local_potential_v / state.emission_temperature_ev);
        const double expected_escaped_current = kIncomingCurrentApm2 * expected_scaling;

        EXPECT_TRUE(evaluation.valid);
        EXPECT_NEAR(evaluation.scaling, expected_scaling, 1.0e-12);
        EXPECT_NEAR(evaluation.escaped_current_a_per_m2, expected_escaped_current,
                    1.0e-18 + std::abs(expected_escaped_current) * 1.0e-12);
    };

    assert_escape_law(-5.0);
    assert_escape_law(0.0);
    assert_escape_law(4.0);
}

TEST(SurfaceBarrierModelsTest,
     BarrierAndRecollectionScalersCoverLowHighFieldAndPotentialWindows)
{
    SCDAT::Material::MaterialProperty material(6, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "window-coverage");
    material.setPhotoelectronTemperatureEv(2.0);

    struct WindowCase
    {
        const char* label;
        double local_potential_v;
        double reference_potential_v;
        double normal_field_v_per_m;
    };

    const WindowCase windows[] = {
        {"low_field_negative_delta", -2.0, 2.0, 5.0e3},
        {"low_field_positive_delta", 2.0, -2.0, 5.0e3},
        {"high_field_negative_delta", -4.0, 4.0, 2.0e5},
        {"high_field_positive_delta", 4.0, -4.0, 2.0e5},
    };

    SCDAT::FieldSolver::VariableBarrierScaler barrier_scaler;
    SCDAT::FieldSolver::SecondaryRecollectionScaler recollection_scaler;
    constexpr double kIncomingCurrentApm2 = 2.0e-6;

    for (const auto& window : windows)
    {
        SCDAT::FieldSolver::SurfaceBarrierState state;
        state.local_potential_v = window.local_potential_v;
        state.reference_potential_v = window.reference_potential_v;
        state.barrier_potential_v = 0.0;
        state.normal_electric_field_v_per_m = window.normal_field_v_per_m;
        state.emission_temperature_ev = 2.0;

        const auto barrier = barrier_scaler.evaluate(material, state, kIncomingCurrentApm2);
        const auto recollection =
            recollection_scaler.evaluate(material, state, kIncomingCurrentApm2);

        EXPECT_TRUE(barrier.valid) << window.label;
        EXPECT_TRUE(recollection.valid) << window.label;
        EXPECT_GE(barrier.scaling, 0.0) << window.label;
        EXPECT_LE(barrier.scaling, 1.0) << window.label;
        EXPECT_GE(recollection.scaling, 0.0) << window.label;
        EXPECT_LE(recollection.scaling, 1.0) << window.label;
        EXPECT_NEAR(recollection.escaped_current_a_per_m2 +
                        recollection.recollected_current_a_per_m2,
                    kIncomingCurrentApm2,
                    1.0e-18 + std::abs(kIncomingCurrentApm2) * 1.0e-12)
            << window.label;
    }
}

TEST(SurfaceBarrierModelsTest,
     BarrierScalersAnalyticalDerivativeMatchesFiniteDifference)
{
    SCDAT::Material::MaterialProperty material(7, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "didv-parity");
    material.setPhotoelectronTemperatureEv(2.0);

    SCDAT::FieldSolver::SurfaceBarrierState base_state;
    base_state.local_potential_v = 3.0;
    base_state.reference_potential_v = 0.0;
    base_state.barrier_potential_v = 0.0;
    base_state.normal_electric_field_v_per_m = 4.0e5;
    base_state.emission_temperature_ev = 2.0;

    constexpr double kIncomingCurrentApm2 = 3.0e-6;
    constexpr double kDeltaV = 1.0e-4;

    const auto assert_didv_alignment =
        [&](const SCDAT::FieldSolver::BarrierCurrentScaler& scaler, const char* label) {
            auto minus_state = base_state;
            auto plus_state = base_state;
            minus_state.local_potential_v -= kDeltaV;
            plus_state.local_potential_v += kDeltaV;

            const auto center = scaler.evaluate(material, base_state, kIncomingCurrentApm2);
            const auto minus = scaler.evaluate(material, minus_state, kIncomingCurrentApm2);
            const auto plus = scaler.evaluate(material, plus_state, kIncomingCurrentApm2);
            const double finite_difference =
                (plus.escaped_current_a_per_m2 - minus.escaped_current_a_per_m2) /
                (2.0 * kDeltaV);

            EXPECT_TRUE(center.valid) << label;
            EXPECT_TRUE(minus.valid) << label;
            EXPECT_TRUE(plus.valid) << label;
            EXPECT_TRUE(std::isfinite(center.didv_a_per_m2_per_v)) << label;
            EXPECT_TRUE(std::isfinite(finite_difference)) << label;
            EXPECT_NEAR(center.didv_a_per_m2_per_v, finite_difference,
                        1.0e-12 + std::abs(finite_difference) * 5.0e-5)
                << label;
        };

    SCDAT::FieldSolver::VariableBarrierScaler variable_scaler;
    SCDAT::FieldSolver::SecondaryRecollectionScaler recollection_scaler;
    assert_didv_alignment(variable_scaler, "variable_barrier_scaler");
    assert_didv_alignment(recollection_scaler, "secondary_recollection_scaler");
}

TEST(SurfaceBarrierModelsTest,
     BarrierDidvDeviationGateExportsSummaryAndChecksThresholds)
{
    SCDAT::Material::MaterialProperty material(8, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "didv-gate");
    material.setPhotoelectronTemperatureEv(2.0);

    struct WindowCase
    {
        const char* label;
        double local_potential_v;
        double reference_potential_v;
        double normal_field_v_per_m;
    };

    struct DidvCaseResult
    {
        std::string scaler_family;
        std::string window_label;
        double local_potential_v = 0.0;
        double reference_potential_v = 0.0;
        double normal_field_v_per_m = 0.0;
        double analytic_didv = 0.0;
        double finite_difference_didv = 0.0;
        double absolute_error = 0.0;
        double relative_error = 0.0;
        bool pass = false;
    };

    const WindowCase windows[] = {
        {"low_field_negative_delta", -2.0, 2.0, 5.0e3},
        {"low_field_positive_delta", 2.0, -2.0, 5.0e3},
        {"high_field_negative_delta", -4.0, 4.0, 2.0e5},
        {"high_field_positive_delta", 4.0, -4.0, 2.0e5},
    };

    constexpr double kIncomingCurrentApm2 = 2.0e-6;
    constexpr double kDeltaV = 1.0e-4;
    constexpr double kAbsoluteErrorThreshold = 1.5e-6;
    constexpr double kRelativeErrorThreshold = 1.1;

    std::vector<DidvCaseResult> results;

    const auto collect_results = [&](const SCDAT::FieldSolver::BarrierCurrentScaler& scaler) {
        for (const auto& window : windows)
        {
            SCDAT::FieldSolver::SurfaceBarrierState center_state;
            center_state.local_potential_v = window.local_potential_v;
            center_state.reference_potential_v = window.reference_potential_v;
            center_state.barrier_potential_v = 0.0;
            center_state.normal_electric_field_v_per_m = window.normal_field_v_per_m;
            center_state.emission_temperature_ev = 2.0;

            auto minus_state = center_state;
            auto plus_state = center_state;
            minus_state.local_potential_v -= kDeltaV;
            plus_state.local_potential_v += kDeltaV;

            const auto center = scaler.evaluate(material, center_state, kIncomingCurrentApm2);
            const auto minus = scaler.evaluate(material, minus_state, kIncomingCurrentApm2);
            const auto plus = scaler.evaluate(material, plus_state, kIncomingCurrentApm2);

            ASSERT_TRUE(center.valid);
            ASSERT_TRUE(minus.valid);
            ASSERT_TRUE(plus.valid);

            const double finite_difference =
                (plus.escaped_current_a_per_m2 - minus.escaped_current_a_per_m2) /
                (2.0 * kDeltaV);
            const double absolute_error =
                std::abs(center.didv_a_per_m2_per_v - finite_difference);
            const double relative_error =
                absolute_error / std::max(1.0e-18, std::abs(center.didv_a_per_m2_per_v));

            DidvCaseResult result;
            result.scaler_family = scaler.family();
            result.window_label = window.label;
            result.local_potential_v = window.local_potential_v;
            result.reference_potential_v = window.reference_potential_v;
            result.normal_field_v_per_m = window.normal_field_v_per_m;
            result.analytic_didv = center.didv_a_per_m2_per_v;
            result.finite_difference_didv = finite_difference;
            result.absolute_error = absolute_error;
            result.relative_error = relative_error;
            result.pass = absolute_error <= kAbsoluteErrorThreshold &&
                          relative_error <= kRelativeErrorThreshold;
            results.push_back(result);
        }
    };

    SCDAT::FieldSolver::VariableBarrierScaler variable_scaler;
    SCDAT::FieldSolver::SecondaryRecollectionScaler recollection_scaler;
    collect_results(variable_scaler);
    collect_results(recollection_scaler);

    double max_absolute_error = 0.0;
    double max_relative_error = 0.0;
    bool overall_pass = true;
    for (const auto& result : results)
    {
        max_absolute_error = std::max(max_absolute_error, result.absolute_error);
        max_relative_error = std::max(max_relative_error, result.relative_error);
        overall_pass = overall_pass && result.pass;
    }

    const auto gate_json_path =
        std::filesystem::current_path() / "surface_barrier_didv_gate.json";
    const auto gate_md_path =
        std::filesystem::current_path() / "surface_barrier_didv_gate.md";

    std::ofstream json_output(gate_json_path);
    ASSERT_TRUE(json_output.is_open());
    json_output << std::setprecision(17);
    json_output << "{\n";
    json_output << "  \"schema_version\": \"surface.barrier.didv.gate.v1\",\n";
    json_output << "  \"absolute_error_threshold\": " << kAbsoluteErrorThreshold << ",\n";
    json_output << "  \"relative_error_threshold\": " << kRelativeErrorThreshold << ",\n";
    json_output << "  \"max_absolute_error\": " << max_absolute_error << ",\n";
    json_output << "  \"max_relative_error\": " << max_relative_error << ",\n";
    json_output << "  \"cases\": [\n";
    for (std::size_t i = 0; i < results.size(); ++i)
    {
        const auto& result = results[i];
        json_output << "    {\n";
        json_output << "      \"scaler_family\": \"" << result.scaler_family
                    << "\",\n";
        json_output << "      \"window\": \"" << result.window_label << "\",\n";
        json_output << "      \"local_potential_v\": " << result.local_potential_v << ",\n";
        json_output << "      \"reference_potential_v\": " << result.reference_potential_v
                    << ",\n";
        json_output << "      \"normal_field_v_per_m\": "
                    << result.normal_field_v_per_m << ",\n";
        json_output << "      \"analytic_didv_a_per_m2_per_v\": " << result.analytic_didv
                    << ",\n";
        json_output << "      \"finite_difference_didv_a_per_m2_per_v\": "
                    << result.finite_difference_didv << ",\n";
        json_output << "      \"absolute_error\": " << result.absolute_error << ",\n";
        json_output << "      \"relative_error\": " << result.relative_error << ",\n";
        json_output << "      \"pass\": " << (result.pass ? "true" : "false") << "\n";
        json_output << (i + 1 < results.size() ? "    },\n" : "    }\n");
    }
    json_output << "  ],\n";
    json_output << "  \"overall_pass\": " << (overall_pass ? "true" : "false") << "\n";
    json_output << "}\n";

    std::ofstream md_output(gate_md_path);
    ASSERT_TRUE(md_output.is_open());
    md_output << "# Surface Barrier DIDV Gate Summary\n\n";
    md_output << "- absolute_error_threshold: " << kAbsoluteErrorThreshold << "\n";
    md_output << "- relative_error_threshold: " << kRelativeErrorThreshold << "\n";
    md_output << "- max_absolute_error: " << max_absolute_error << "\n";
    md_output << "- max_relative_error: " << max_relative_error << "\n";
    md_output << "- overall_pass: " << (overall_pass ? "true" : "false") << "\n\n";
    md_output << "| scaler_family | window | analytic_didv | finite_difference_didv | "
                 "absolute_error | relative_error | pass |\n";
    md_output << "| --- | --- | --- | --- | --- | --- | --- |\n";
    for (const auto& result : results)
    {
        md_output << "| " << result.scaler_family << " | " << result.window_label << " | "
                  << result.analytic_didv << " | " << result.finite_difference_didv << " | "
                  << result.absolute_error << " | " << result.relative_error << " | "
                  << (result.pass ? "true" : "false") << " |\n";
    }

    EXPECT_TRUE(std::filesystem::exists(gate_json_path));
    EXPECT_TRUE(std::filesystem::exists(gate_md_path));
    EXPECT_LE(max_absolute_error, kAbsoluteErrorThreshold);
    EXPECT_LE(max_relative_error, kRelativeErrorThreshold);
    EXPECT_TRUE(overall_pass);
}

TEST(SurfaceBarrierModelsTest, BarrierScalersStayBoundedAcrossExtremeInputVectors)
{
    SCDAT::Material::MaterialProperty material(9, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "extreme-vectors");
    material.setPhotoelectronTemperatureEv(2.0);

    struct ExtremeCase
    {
        const char* label;
        double local_potential_v;
        double reference_potential_v;
        double barrier_potential_v;
        double normal_field_v_per_m;
        double emission_temperature_ev;
        double emitted_current_a_per_m2;
    };

    const ExtremeCase cases[] = {
        {"near_zero", 1.0e-12, -1.0e-12, 0.0, 1.0e-6, 1.0e-12, 1.0e-24},
        {"huge_positive", 1.0e6, -1.0e6, 1.0e6, 1.0e9, 1.0e-6, 5.0e-6},
        {"huge_negative", -1.0e6, 1.0e6, -1.0e6, -1.0e9, 1.0e-6, 5.0e-6},
        {"sign_flip", 80.0, -40.0, 5.0, -5.0e7, 0.5, 2.5e-6},
    };

    const auto run_check = [&](const SCDAT::FieldSolver::BarrierCurrentScaler& scaler) {
        for (const auto& item : cases)
        {
            SCDAT::FieldSolver::SurfaceBarrierState state;
            state.local_potential_v = item.local_potential_v;
            state.reference_potential_v = item.reference_potential_v;
            state.barrier_potential_v = item.barrier_potential_v;
            state.normal_electric_field_v_per_m = item.normal_field_v_per_m;
            state.emission_temperature_ev = item.emission_temperature_ev;

            const auto evaluation =
                scaler.evaluate(material, state, item.emitted_current_a_per_m2);
            const double emitted = std::max(0.0, item.emitted_current_a_per_m2);

            EXPECT_TRUE(evaluation.valid) << item.label;
            EXPECT_TRUE(std::isfinite(evaluation.scaling)) << item.label;
            EXPECT_TRUE(std::isfinite(evaluation.escaped_current_a_per_m2)) << item.label;
            EXPECT_TRUE(std::isfinite(evaluation.recollected_current_a_per_m2)) << item.label;
            EXPECT_TRUE(std::isfinite(evaluation.didv_a_per_m2_per_v)) << item.label;
            EXPECT_GE(evaluation.scaling, 0.0) << item.label;
            EXPECT_LE(evaluation.scaling, 1.0) << item.label;
            EXPECT_NEAR(evaluation.escaped_current_a_per_m2 +
                            evaluation.recollected_current_a_per_m2,
                        emitted,
                        1.0e-18 + std::abs(emitted) * 1.0e-12)
                << item.label;
        }
    };

    SCDAT::FieldSolver::VariableBarrierScaler variable_scaler;
    SCDAT::FieldSolver::SecondaryRecollectionScaler recollection_scaler;
    run_check(variable_scaler);
    run_check(recollection_scaler);
}

TEST(SurfaceBarrierModelsTest, BarrierScalersGuardAgainstNaNAndInfInputs)
{
    SCDAT::Material::MaterialProperty material(10, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "nan-inf-guard");
    material.setPhotoelectronTemperatureEv(2.0);

    SCDAT::FieldSolver::SurfaceBarrierState pathological_state;
    pathological_state.local_potential_v = std::numeric_limits<double>::quiet_NaN();
    pathological_state.reference_potential_v = std::numeric_limits<double>::infinity();
    pathological_state.barrier_potential_v = -std::numeric_limits<double>::infinity();
    pathological_state.normal_electric_field_v_per_m = std::numeric_limits<double>::quiet_NaN();
    pathological_state.emission_temperature_ev = std::numeric_limits<double>::quiet_NaN();

    SCDAT::FieldSolver::VariableBarrierScaler variable_scaler;
    SCDAT::FieldSolver::SecondaryRecollectionScaler recollection_scaler;

    const auto variable_eval = variable_scaler.evaluate(
        material, pathological_state, std::numeric_limits<double>::infinity());
    const auto recollection_eval = recollection_scaler.evaluate(
        material, pathological_state, std::numeric_limits<double>::quiet_NaN());

    const auto assert_sanitized = [](const SCDAT::FieldSolver::SurfaceBarrierEvaluation& eval) {
        EXPECT_TRUE(eval.valid);
        EXPECT_TRUE(std::isfinite(eval.scaling));
        EXPECT_TRUE(std::isfinite(eval.escaped_current_a_per_m2));
        EXPECT_TRUE(std::isfinite(eval.recollected_current_a_per_m2));
        EXPECT_TRUE(std::isfinite(eval.didv_a_per_m2_per_v));
        EXPECT_GE(eval.scaling, 0.0);
        EXPECT_LE(eval.scaling, 1.0);
        EXPECT_GE(eval.escaped_current_a_per_m2, 0.0);
        EXPECT_GE(eval.recollected_current_a_per_m2, 0.0);
    };

    assert_sanitized(variable_eval);
    assert_sanitized(recollection_eval);

    SCDAT::FieldSolver::SurfaceEmissionBarrierInputs inputs;
    inputs.electron_collection_a_per_m2 = std::numeric_limits<double>::quiet_NaN();
    inputs.ion_collection_a_per_m2 = std::numeric_limits<double>::infinity();
    inputs.secondary_emission_a_per_m2 = std::numeric_limits<double>::infinity();
    inputs.ion_secondary_emission_a_per_m2 = std::numeric_limits<double>::quiet_NaN();
    inputs.backscatter_emission_a_per_m2 = std::numeric_limits<double>::infinity();
    inputs.photo_emission_a_per_m2 = std::numeric_limits<double>::quiet_NaN();
    inputs.photo_incidence_scale = std::numeric_limits<double>::quiet_NaN();
    inputs.user_secondary_scale = std::numeric_limits<double>::infinity();
    inputs.user_ion_secondary_scale = std::numeric_limits<double>::infinity();
    inputs.user_backscatter_scale = std::numeric_limits<double>::quiet_NaN();
    inputs.user_photo_scale = std::numeric_limits<double>::infinity();
    inputs.fallback_secondary_yield = std::numeric_limits<double>::quiet_NaN();
    inputs.fallback_ion_secondary_yield = std::numeric_limits<double>::quiet_NaN();
    inputs.fallback_backscatter_yield = std::numeric_limits<double>::quiet_NaN();

    const auto outputs =
        SCDAT::FieldSolver::evaluateSurfaceEmissionBarrierBundle(material, pathological_state, inputs);

    EXPECT_TRUE(outputs.valid);
    EXPECT_TRUE(std::isfinite(outputs.secondary_emission_scale));
    EXPECT_TRUE(std::isfinite(outputs.ion_secondary_emission_scale));
    EXPECT_TRUE(std::isfinite(outputs.backscatter_scale));
    EXPECT_TRUE(std::isfinite(outputs.photo_emission_scale));
    EXPECT_GE(outputs.secondary_emission_scale, 0.0);
    EXPECT_LE(outputs.secondary_emission_scale, 4.0);
    EXPECT_GE(outputs.ion_secondary_emission_scale, 0.0);
    EXPECT_LE(outputs.ion_secondary_emission_scale, 4.0);
    EXPECT_GE(outputs.backscatter_scale, 0.0);
    EXPECT_LE(outputs.backscatter_scale, 4.0);
    EXPECT_GE(outputs.photo_emission_scale, 0.0);
    EXPECT_LE(outputs.photo_emission_scale, 4.0);
}
