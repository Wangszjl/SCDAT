#include "../include/SurfaceDistributionFunction.h"
#include "../include/ParticleSource.h"

#include <gtest/gtest.h>

TEST(SurfaceDistributionFunctionTest, MaxwellianFluxMomentsUseExpectedFamily)
{
    SCDAT::Particle::IsotropicMaxwellianFluxDistribution distribution(5.0e11, 12.0);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_isotropic_maxwellian_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_NEAR(moments.mean_energy_ev, 24.0, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest, BiMaxwellianMeanEnergyTracksMixture)
{
    SCDAT::Particle::BiMaxwellianFluxDistribution distribution(4.0e11, 5.0, 1.0e11, 30.0);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_bimaxwellian_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.mean_energy_ev, 10.0);
    EXPECT_LT(moments.mean_energy_ev, 30.0);
}

TEST(SurfaceDistributionFunctionTest, TabulatedDistributionComputesWeightedMoments)
{
    SCDAT::Particle::TabulatedSurfaceDistributionFunction distribution(
        {10.0, 50.0, 200.0}, {1.0, 2.0, 1.0}, {0.2, 0.6, 0.9});
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_tabulated_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_NEAR(moments.total_flux, 4.0, 1.0e-12);
    EXPECT_GT(moments.mean_energy_ev, 50.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
}

TEST(SurfaceDistributionFunctionTest, TabulatedVelocityDistributionComputesWeightedMoments)
{
    SCDAT::Particle::TabulatedVelocitySurfaceDistributionFunction distribution(
        {2.0e4, 4.0e4, 8.0e4}, {0.5, 1.0, 0.5}, 1.0e-10, {0.3, 0.6, 0.9});
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_tabulated_velocity_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_NEAR(moments.total_flux, 2.0, 1.0e-12);
    EXPECT_GT(moments.mean_energy_ev, 0.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
}

TEST(SurfaceDistributionFunctionTest, SamplerTracksDirtyAndCachedState)
{
    SCDAT::Particle::SurfaceFluxSampler sampler;
    EXPECT_TRUE(sampler.dirty());

    SCDAT::Particle::SurfaceFluxMoments moments;
    moments.total_flux = 12.0;
    moments.valid = true;
    sampler.cacheMoments(moments);

    EXPECT_FALSE(sampler.dirty());
    EXPECT_DOUBLE_EQ(sampler.cachedMoments().total_flux, 12.0);

    sampler.markDirty();
    EXPECT_TRUE(sampler.dirty());
}

TEST(SurfaceDistributionFunctionTest, EmissionSamplerTracksDensificationAndEmissionCount)
{
    SCDAT::Particle::SurfaceEmissionSampler sampler;
    sampler.setSuperParticleDensification(250.0);
    sampler.noteEmissionCount(42.5);

    EXPECT_DOUBLE_EQ(sampler.superParticleDensification(), 250.0);
    EXPECT_DOUBLE_EQ(sampler.cachedEmissionCount(), 42.5);
}

TEST(SurfaceDistributionFunctionTest, LegacySpectrumInterpolationInterpolatesInteriorPoint)
{
    const double interpolated = SCDAT::Particle::interpolateLegacySpectrum(
        {1.0, 2.0, 4.0}, {10.0, 20.0, 40.0}, 3.0);
    EXPECT_GT(interpolated, 20.0);
    EXPECT_LT(interpolated, 40.0);
}

TEST(SurfaceDistributionFunctionTest, LegacySpectrumLine3SmoothingPreservesPointCount)
{
    const auto smoothed = SCDAT::Particle::smoothLegacySpectrumLine3(
        {1.0, 2.0, 3.0, 4.0}, {10.0, 40.0, 20.0, 50.0});
    ASSERT_EQ(smoothed.size(), 4U);
    EXPECT_TRUE(std::isfinite(smoothed[0]));
    EXPECT_TRUE(std::isfinite(smoothed[3]));
}

TEST(SurfaceDistributionFunctionTest, LegacyComponentFluxHelpersReturnPositiveFiniteValues)
{
    const double electron_flux =
        SCDAT::Particle::legacyElectronComponentFlux(1.0e11, 10.0, 20.0);
    const double ion_flux =
        SCDAT::Particle::legacyIonComponentFlux(8.0e10, 5.0, 16.0, 12.0);
    EXPECT_GT(electron_flux, 0.0);
    EXPECT_GT(ion_flux, 0.0);
    EXPECT_TRUE(std::isfinite(electron_flux));
    EXPECT_TRUE(std::isfinite(ion_flux));
}

TEST(SurfaceDistributionFunctionTest, LegacyFluxTableBuilderBuildsSharedEnergyGridAndPopulationFlux)
{
    SCDAT::Particle::LegacyFluxTableBuildRequest request;
    request.has_electron_spectrum = true;
    request.has_ion_spectrum = true;
    request.fallback_electron_density_m3 = 1.0e11;
    request.fallback_electron_temperature_ev = 8.0;
    request.fallback_ion_density_m3 = 5.0e10;
    request.fallback_ion_temperature_ev = 3.0;
    request.fallback_ion_mass_amu = 16.0;
    request.sample_count = 32;
    request.max_population_count = 2;
    request.electron_spectrum.energy_grid_ev = {1.0, 5.0, 25.0};
    request.electron_spectrum.differential_number_flux = {2.0, 3.0, 1.0};
    request.electron_spectrum.populations = {
        {6.0e10, 6.0, 0.0, 1.0, -1.0, 0.0},
        {4.0e10, 12.0, 0.0, 1.0, -1.0, 0.0},
    };
    request.ion_spectrum.energy_grid_ev = {0.5, 4.0, 12.0};
    request.ion_spectrum.differential_number_flux = {1.0, 2.0, 0.5};
    request.ion_spectrum.populations = {
        {5.0e10, 3.0, 0.0, 16.0, 1.0, 0.0},
    };

    const auto tables = SCDAT::Particle::buildLegacyFluxTables(request);
    ASSERT_EQ(tables.energy_ev.size(), 32U);
    ASSERT_EQ(tables.width_ev.size(), tables.energy_ev.size());
    ASSERT_EQ(tables.total_electron_flux.size(), tables.energy_ev.size());
    ASSERT_EQ(tables.total_ion_flux.size(), tables.energy_ev.size());
    ASSERT_EQ(tables.electron_component_flux.size(), 2U);
    ASSERT_EQ(tables.ion_component_flux.size(), 1U);
    EXPECT_EQ(tables.spectrum_energy_grid_ev.size(), 3U);
    EXPECT_EQ(tables.smoothed_electron_flux.size(), 3U);
    EXPECT_EQ(tables.smoothed_ion_flux.size(), 3U);
    EXPECT_GT(tables.average_spectrum_electron_temperature_ev, 0.0);
    EXPECT_GT(tables.average_spectrum_ion_temperature_ev, 0.0);
}

TEST(SurfaceDistributionFunctionTest, SynthesisBuildsHybridCollectionScalesFromFluxMoments)
{
    SCDAT::Particle::SurfaceDistributionSynthesisRequest request;
    request.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
    request.electron_flux_moments.total_flux = 4.0;
    request.electron_flux_moments.mean_energy_ev = 18.0;
    request.electron_flux_moments.mean_cos_incidence = 0.6;
    request.electron_flux_moments.valid = true;
    request.ion_flux_moments.total_flux = 2.0;
    request.ion_flux_moments.mean_energy_ev = 7.0;
    request.ion_flux_moments.mean_cos_incidence = 0.8;
    request.ion_flux_moments.valid = true;
    request.fallback_electron_energy_ev = 10.0;
    request.fallback_ion_energy_ev = 5.0;
    request.local_reference_shift_v = 1.5;
    request.incidence_scale = 0.9;
    request.normalized_flow_alignment = 0.7;
    request.local_field_factor = 1.1;
    request.electron_population_count = 2;
    request.ion_population_count = 3;

    const auto result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(request);
    EXPECT_TRUE(result.valid);
    EXPECT_GT(result.electron_collection_scale, 1.0);
    EXPECT_GT(result.ion_collection_scale, 1.0);
    EXPECT_NEAR(result.electron_characteristic_energy_ev, 18.0, 1.0e-12);
    EXPECT_NEAR(result.ion_characteristic_energy_ev, 7.0, 1.0e-12);
    EXPECT_NEAR(result.photo_incidence_scale, 0.9, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest, BuildSurfaceFluxMomentsConstructsWakeAndHybridSources)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest wake_request;
    wake_request.model = SCDAT::Particle::SurfaceDistributionModel::WakeAnisotropic;
    wake_request.electron_temperature_ev = 8.0;
    wake_request.ion_temperature_ev = 3.0;
    const auto wake_result = SCDAT::Particle::buildSurfaceFluxMoments(wake_request);
    EXPECT_TRUE(wake_result.valid);
    EXPECT_TRUE(wake_result.electron_flux_moments.valid);
    EXPECT_TRUE(wake_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionBuildRequest hybrid_request;
    hybrid_request.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
    hybrid_request.electron_density_m3 = 1.0e11;
    hybrid_request.electron_temperature_ev = 10.0;
    hybrid_request.ion_density_m3 = 5.0e10;
    hybrid_request.ion_temperature_ev = 4.0;
    hybrid_request.electron_populations = {
        {3.0e10, 5.0},
        {7.0e10, 20.0},
    };
    hybrid_request.ion_populations = {
        {2.0e10, 2.0},
        {3.0e10, 6.0},
    };
    const auto hybrid_result = SCDAT::Particle::buildSurfaceFluxMoments(hybrid_request);
    EXPECT_TRUE(hybrid_result.valid);
    EXPECT_GT(hybrid_result.electron_flux_moments.mean_energy_ev, 5.0);
    EXPECT_GT(hybrid_result.ion_flux_moments.mean_energy_ev, 2.0);
}

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportTabulatedVelocityAndNonLocalizedModels)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest tabulated_request;
    tabulated_request.model = SCDAT::Particle::SurfaceDistributionModel::TabulatedVelocity;
    tabulated_request.electron_density_m3 = 1.0e11;
    tabulated_request.electron_temperature_ev = 9.0;
    tabulated_request.ion_density_m3 = 6.0e10;
    tabulated_request.ion_temperature_ev = 4.0;
    tabulated_request.ion_directed_velocity_m_per_s = 2.0e3;
    const auto tabulated_result =
        SCDAT::Particle::buildSurfaceFluxMoments(tabulated_request);
    EXPECT_TRUE(tabulated_result.valid);
    EXPECT_TRUE(tabulated_result.electron_flux_moments.valid);
    EXPECT_TRUE(tabulated_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionBuildRequest non_local_request = tabulated_request;
    non_local_request.model = SCDAT::Particle::SurfaceDistributionModel::NonLocalizedHybrid;
    non_local_request.ion_directed_velocity_m_per_s = 9.0e3;
    const auto non_local_result =
        SCDAT::Particle::buildSurfaceFluxMoments(non_local_request);
    EXPECT_TRUE(non_local_result.valid);
    EXPECT_TRUE(non_local_result.electron_flux_moments.valid);
    EXPECT_TRUE(non_local_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest synthesis_request;
    synthesis_request.model = SCDAT::Particle::SurfaceDistributionModel::NonLocalizedHybrid;
    synthesis_request.electron_flux_moments = non_local_result.electron_flux_moments;
    synthesis_request.ion_flux_moments = non_local_result.ion_flux_moments;
    synthesis_request.fallback_electron_energy_ev = 8.0;
    synthesis_request.fallback_ion_energy_ev = 3.0;
    synthesis_request.incidence_scale = 0.9;
    synthesis_request.normalized_flow_alignment = -0.4;
    synthesis_request.wake_factor = 0.7;
    synthesis_request.local_field_factor = 1.1;
    const auto synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(synthesis_request);
    EXPECT_TRUE(synthesized.valid);
    EXPECT_GT(synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(synthesized.ion_collection_scale, 0.1);
    EXPECT_GE(synthesized.photo_incidence_scale, 0.0);
}

TEST(SurfaceDistributionFunctionTest, BundleEvaluationBuildsAndSynthesizesDistributionMoments)
{
    SCDAT::Particle::SurfaceDistributionBundleRequest request;
    request.build_request.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
    request.build_request.electron_density_m3 = 1.0e11;
    request.build_request.electron_temperature_ev = 12.0;
    request.build_request.ion_density_m3 = 8.0e10;
    request.build_request.ion_temperature_ev = 5.0;
    request.build_request.electron_populations = {
        {4.0e10, 6.0},
        {6.0e10, 18.0},
    };
    request.build_request.ion_populations = {
        {3.0e10, 3.0},
        {5.0e10, 7.0},
    };
    request.synthesis_request.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
    request.synthesis_request.fallback_electron_energy_ev = 10.0;
    request.synthesis_request.fallback_ion_energy_ev = 4.0;
    request.synthesis_request.local_reference_shift_v = 1.0;
    request.synthesis_request.incidence_scale = 0.85;
    request.synthesis_request.normalized_flow_alignment = 0.65;
    request.synthesis_request.local_field_factor = 1.1;
    request.synthesis_request.wake_factor = 0.15;
    request.synthesis_request.patch_role = true;
    request.synthesis_request.electron_population_count = 2;
    request.synthesis_request.ion_population_count = 2;

    const auto result = SCDAT::Particle::evaluateSurfaceDistributionBundle(request);
    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.built_flux.valid);
    EXPECT_TRUE(result.synthesized.valid);
    EXPECT_GT(result.synthesized.electron_collection_scale, 1.0);
    EXPECT_GT(result.synthesized.ion_collection_scale, 1.0);
}

TEST(SurfaceDistributionFunctionTest, EnvironmentBundleBuildsRequestAndEvaluatesMoments)
{
    SCDAT::Particle::SurfaceDistributionEnvironment environment;
    environment.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
    environment.electron_density_m3 = 1.2e11;
    environment.electron_temperature_ev = 14.0;
    environment.ion_density_m3 = 6.0e10;
    environment.ion_temperature_ev = 5.5;
    environment.ion_directed_velocity_m_per_s = 1.8e3;
    environment.electron_populations = {
        {5.0e10, 7.0},
        {7.0e10, 19.0},
    };
    environment.ion_populations = {
        {2.0e10, 3.0},
        {4.0e10, 8.0},
    };
    environment.fallback_electron_energy_ev = 11.0;
    environment.fallback_ion_energy_ev = 4.5;
    environment.local_reference_shift_v = 1.5;
    environment.incidence_scale = 0.82;
    environment.normalized_flow_alignment = 0.7;
    environment.local_field_factor = 1.15;
    environment.wake_factor = 0.2;
    environment.patch_role = true;

    const auto request =
        SCDAT::Particle::makeSurfaceDistributionBundleRequest(environment);
    EXPECT_EQ(request.build_request.model, environment.model);
    EXPECT_DOUBLE_EQ(request.build_request.electron_density_m3,
                     environment.electron_density_m3);
    EXPECT_EQ(request.build_request.electron_populations.size(), 2U);
    EXPECT_DOUBLE_EQ(request.synthesis_request.incidence_scale,
                     environment.incidence_scale);
    EXPECT_TRUE(request.synthesis_request.patch_role);

    const auto result =
        SCDAT::Particle::evaluateSurfaceDistributionEnvironment(environment);
    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.built_flux.valid);
    EXPECT_TRUE(result.synthesized.valid);
    EXPECT_GT(result.built_flux.electron_flux_moments.total_flux, 0.0);
    EXPECT_GT(result.synthesized.electron_collection_scale, 1.0);
    EXPECT_GT(result.synthesized.ion_collection_scale, 1.0);
    EXPECT_NEAR(result.synthesized.photo_incidence_scale, 0.82, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest, RoleEnvironmentBuilderComputesIncidenceWakeAndFieldTerms)
{
    SCDAT::Particle::SurfaceDistributionRoleInputs inputs;
    inputs.model = SCDAT::Particle::SurfaceDistributionModel::WakeAnisotropic;
    inputs.electron_density_m3 = 1.0e11;
    inputs.electron_temperature_ev = 10.0;
    inputs.ion_density_m3 = 8.0e10;
    inputs.ion_temperature_ev = 4.0;
    inputs.fallback_electron_energy_ev = 9.0;
    inputs.fallback_ion_energy_ev = 3.5;
    inputs.local_reference_shift_v = 2.0;
    inputs.patch_incidence_angle_deg = 60.0;
    inputs.bulk_flow_velocity_m_per_s = 7500.0;
    inputs.projected_speed_m_per_s = -3750.0;
    inputs.flow_alignment_cosine = -0.5;
    inputs.normal_electric_field_v_per_m = 3.0e4;
    inputs.share_patch_distribution = false;
    inputs.patch_role = true;

    const auto environment = SCDAT::Particle::makeSurfaceDistributionEnvironment(inputs);
    EXPECT_EQ(environment.model, inputs.model);
    EXPECT_NEAR(environment.incidence_scale, 0.5, 1.0e-12);
    EXPECT_NEAR(environment.normalized_flow_alignment, -0.5, 1.0e-12);
    EXPECT_NEAR(environment.wake_factor, 0.75, 1.0e-12);
    EXPECT_NEAR(environment.local_field_factor, 1.15, 1.0e-12);
    EXPECT_TRUE(environment.patch_role);
}

TEST(SurfaceDistributionFunctionTest,
     RoleEnvironmentBuilderCoversBodyPatchIncidenceFlowAndTemperatureCombinations)
{
    struct RoleCase
    {
        const char* label;
        bool patch_role;
        bool share_patch_distribution;
        double patch_incidence_angle_deg;
        double bulk_flow_velocity_m_per_s;
        double projected_speed_m_per_s;
        double flow_alignment_cosine;
        double normal_electric_field_v_per_m;
        double expected_incidence_scale;
        double expected_flow_alignment;
        double expected_field_factor;
        double expected_wake_factor;
    };

    const RoleCase cases[] = {
        {"patch_projected_flow", true, false, 60.0, 7500.0, -3750.0, 0.25, 3.0e4, 0.5, -0.5,
         1.15, 0.75},
        {"body_fallback_alignment", false, false, 15.0, 0.0, 0.0, 0.3, 1.0e4, 0.9659258262890683,
         0.3, 1.05, 0.35},
        {"shared_patch_distribution", true, true, 80.0, 7200.0, 1200.0, -0.6, 8.0e4, 1.0,
         0.16666666666666666, 1.0, 0.4166666666666667},
    };

    for (const auto& role_case : cases)
    {
        SCDAT::Particle::SurfaceDistributionRoleInputs inputs;
        inputs.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
        inputs.electron_density_m3 = 1.5e11;
        inputs.electron_temperature_ev = role_case.patch_role ? 8.0 : 15.0;
        inputs.ion_density_m3 = 7.0e10;
        inputs.ion_temperature_ev = role_case.patch_role ? 3.0 : 6.5;
        inputs.fallback_electron_energy_ev = inputs.electron_temperature_ev;
        inputs.fallback_ion_energy_ev = inputs.ion_temperature_ev;
        inputs.patch_incidence_angle_deg = role_case.patch_incidence_angle_deg;
        inputs.bulk_flow_velocity_m_per_s = role_case.bulk_flow_velocity_m_per_s;
        inputs.projected_speed_m_per_s = role_case.projected_speed_m_per_s;
        inputs.flow_alignment_cosine = role_case.flow_alignment_cosine;
        inputs.normal_electric_field_v_per_m = role_case.normal_electric_field_v_per_m;
        inputs.share_patch_distribution = role_case.share_patch_distribution;
        inputs.patch_role = role_case.patch_role;

        const auto environment = SCDAT::Particle::makeSurfaceDistributionEnvironment(inputs);

        EXPECT_NEAR(environment.incidence_scale, role_case.expected_incidence_scale, 1.0e-12)
            << role_case.label;
        EXPECT_NEAR(environment.normalized_flow_alignment, role_case.expected_flow_alignment,
                    1.0e-12)
            << role_case.label;
        EXPECT_NEAR(environment.local_field_factor, role_case.expected_field_factor, 1.0e-12)
            << role_case.label;
        EXPECT_NEAR(environment.wake_factor, role_case.expected_wake_factor, 1.0e-12)
            << role_case.label;
        EXPECT_EQ(environment.patch_role, role_case.patch_role) << role_case.label;
        EXPECT_NEAR(environment.fallback_electron_energy_ev, inputs.fallback_electron_energy_ev,
                    1.0e-12)
            << role_case.label;
        EXPECT_NEAR(environment.fallback_ion_energy_ev, inputs.fallback_ion_energy_ev, 1.0e-12)
            << role_case.label;
    }
}

TEST(SurfaceDistributionFunctionTest,
     RoleInputBuilderBuildsPopulationsAndFallbackFromResolvedSpectra)
{
    SCDAT::Particle::SurfaceDistributionRoleBuildRequest request;
    request.model = SCDAT::Particle::SurfaceDistributionModel::MultiPopulationHybrid;
    request.electron_density_m3 = 1.0e11;
    request.electron_temperature_ev = 6.0;
    request.ion_density_m3 = 7.0e10;
    request.ion_temperature_ev = 3.0;
    request.ion_directed_velocity_m_per_s = 1200.0;
    request.local_reference_shift_v = 1.2;
    request.patch_incidence_angle_deg = 45.0;
    request.bulk_flow_velocity_m_per_s = 7500.0;
    request.projected_speed_m_per_s = -3200.0;
    request.flow_alignment_cosine = -0.4;
    request.normal_electric_field_v_per_m = 2.0e4;
    request.share_patch_distribution = false;
    request.patch_role = true;
    request.electron_spectrum.energy_grid_ev = {2.0, 10.0, 30.0};
    request.electron_spectrum.differential_number_flux = {1.0, 3.0, 1.0};
    request.electron_spectrum.populations = {
        {6.0e10, 5.0, 0.0, 1.0, -1.0, 0.0},
        {4.0e10, 12.0, 0.0, 1.0, -1.0, 0.0},
    };
    request.ion_spectrum.populations = {
        {7.0e10, 3.0, 0.0, 16.0, 1.0, 0.0},
    };

    const auto inputs = SCDAT::Particle::makeSurfaceDistributionRoleInputs(request);
    EXPECT_EQ(inputs.model, request.model);
    EXPECT_DOUBLE_EQ(inputs.electron_density_m3, request.electron_density_m3);
    EXPECT_DOUBLE_EQ(inputs.ion_density_m3, request.ion_density_m3);
    ASSERT_EQ(inputs.electron_populations.size(), 2U);
    ASSERT_EQ(inputs.ion_populations.size(), 1U);
    EXPECT_DOUBLE_EQ(inputs.electron_populations[0].density_m3, 6.0e10);
    EXPECT_DOUBLE_EQ(inputs.electron_populations[1].temperature_ev, 12.0);
    EXPECT_DOUBLE_EQ(inputs.ion_populations[0].density_m3, 7.0e10);
    EXPECT_GT(inputs.fallback_electron_energy_ev, request.electron_temperature_ev);
    EXPECT_GE(inputs.fallback_ion_energy_ev, request.ion_temperature_ev);
    EXPECT_DOUBLE_EQ(inputs.local_reference_shift_v, request.local_reference_shift_v);
    EXPECT_DOUBLE_EQ(inputs.projected_speed_m_per_s, request.projected_speed_m_per_s);
    EXPECT_TRUE(inputs.patch_role);
}

TEST(SurfaceDistributionFunctionTest, ResolvedSpectrumMomentsPreferPopulationStatistics)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    spectrum.populations = {
        {2.0e11, 4.0, 0.0, 16.0, 1.0, 0.0},
        {1.0e11, 10.0, 0.0, 32.0, 1.0, 0.0},
    };

    const auto moments = SCDAT::Particle::resolveSpectrumMoments(
        spectrum, 5.0e10, 3.0, 20.0);
    EXPECT_TRUE(moments.valid);
    EXPECT_TRUE(moments.used_population_moments);
    EXPECT_NEAR(moments.density_m3, 3.0e11, 1.0);
    EXPECT_NEAR(moments.characteristic_energy_ev, 6.0, 1.0e-12);
    EXPECT_NEAR(moments.average_mass_amu, 21.333333333333332, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest, ResolvedSpectrumMomentsUseDiscreteFluxFallback)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    spectrum.energy_grid_ev = {0.0, 10.0, 20.0};
    spectrum.differential_number_flux = {1.0, 1.0, 1.0};

    const auto moments = SCDAT::Particle::resolveSpectrumMoments(
        spectrum, 7.0e10, 4.0, 28.0);
    EXPECT_TRUE(moments.valid);
    EXPECT_FALSE(moments.used_population_moments);
    EXPECT_NEAR(moments.density_m3, 7.0e10, 1.0);
    EXPECT_NEAR(moments.characteristic_energy_ev, 10.0, 1.0e-12);
    EXPECT_NEAR(moments.average_mass_amu, 28.0, 1.0e-12);

    EXPECT_NEAR(SCDAT::Particle::resolvedSpectrumDensity(spectrum, 3.0), 3.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::Particle::resolvedSpectrumCharacteristicEnergyEv(spectrum, 1.0), 10.0,
                1.0e-12);
    EXPECT_NEAR(SCDAT::Particle::resolvedSpectrumAverageMassAmu(spectrum, 12.0), 12.0,
                1.0e-12);
}

TEST(SurfaceDistributionFunctionTest,
     ResolvedSpectrumDiscreteFluxHelpersDetectIntegrateAndAverageEnergy)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    spectrum.model = SCDAT::Particle::SpatialSamplingModel::TABULATED;
    spectrum.energy_grid_ev = {0.0, 10.0, 20.0};
    spectrum.differential_number_flux = {1.0, 1.0, 1.0};

    EXPECT_TRUE(SCDAT::Particle::usesResolvedSpectrumDiscreteFlux(spectrum));
    const double total_flux = SCDAT::Particle::integrateResolvedSpectrumFlux(
        spectrum, [](double, double flux) { return flux; });
    const double weighted_energy_flux = SCDAT::Particle::integrateResolvedSpectrumFlux(
        spectrum, [](double energy_ev, double flux) { return energy_ev * flux; });
    EXPECT_NEAR(total_flux, 20.0, 1.0e-12);
    EXPECT_NEAR(weighted_energy_flux, 200.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::Particle::resolvedSpectrumAverageEnergyEv(spectrum, 1.0), 10.0,
                1.0e-12);
}

TEST(SurfaceDistributionFunctionTest,
     ResolvedSpectrumDiscreteFluxHelpersFallbackWhenSpectrumIsInvalid)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    spectrum.energy_grid_ev = {5.0, 10.0};
    spectrum.differential_number_flux = {1.0};

    EXPECT_FALSE(SCDAT::Particle::usesResolvedSpectrumDiscreteFlux(spectrum));
    const double total_flux = SCDAT::Particle::integrateResolvedSpectrumFlux(
        spectrum, [](double, double flux) { return flux; });
    EXPECT_NEAR(total_flux, 0.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::Particle::resolvedSpectrumAverageEnergyEv(spectrum, 7.5), 7.5,
                1.0e-12);
}

TEST(SurfaceDistributionFunctionTest,
     PopulationMaxwellianAverageHelperReturnsExpectedConstantAndRampAverage)
{
    const double constant_average = SCDAT::Particle::integratePopulationMaxwellianAverageEv(
        12.0, [](double) { return 0.75; });
    const double ramp_average = SCDAT::Particle::integratePopulationMaxwellianAverageEv(
        5.0, [](double energy_ev) { return energy_ev; });

    EXPECT_NEAR(constant_average, 0.75, 3.0e-3);
    EXPECT_NEAR(ramp_average, 10.0, 1.0e-2);
}

TEST(SurfaceDistributionFunctionTest, SampleResolvedSpectrumEnergyEvUsesFluxWeights)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    spectrum.energy_grid_ev = {5.0, 15.0, 30.0};
    spectrum.differential_number_flux = {0.0, 0.0, 2.0};

    std::mt19937 rng(20260412u);
    const double sampled =
        SCDAT::Particle::sampleResolvedSpectrumEnergyEv(spectrum, rng, 2.0);
    EXPECT_NEAR(sampled, 30.0, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest, SampleResolvedSpectrumEnergyEvFallsBackWhenUnavailable)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    std::mt19937 rng(20260412u);

    const double sampled =
        SCDAT::Particle::sampleResolvedSpectrumEnergyEv(spectrum, rng, 7.5);
    EXPECT_NEAR(sampled, 7.5, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest,
     FluxSamplerRefreshesCachedMomentsWhenDistributionChangesAndSamplingUsesNewSpectrum)
{
    SCDAT::Particle::SurfaceFluxSampler sampler;

    SCDAT::Particle::TabulatedSurfaceDistributionFunction baseline_distribution(
        {5.0, 15.0, 30.0}, {3.0, 1.0, 0.0}, {0.4, 0.7, 0.9});
    const auto baseline_moments = baseline_distribution.computeMoments();
    ASSERT_TRUE(baseline_moments.valid);
    sampler.cacheMoments(baseline_moments);
    ASSERT_FALSE(sampler.dirty());
    EXPECT_NEAR(sampler.cachedMoments().mean_energy_ev, 7.5, 1.0e-12);

    SCDAT::Particle::ResolvedSpectrum baseline_spectrum;
    baseline_spectrum.model = SCDAT::Particle::SpatialSamplingModel::TABULATED;
    baseline_spectrum.energy_grid_ev = {5.0, 15.0, 30.0};
    baseline_spectrum.differential_number_flux = {3.0, 1.0, 0.0};

    std::mt19937 baseline_rng(20260413u);
    const double baseline_sample =
        SCDAT::Particle::sampleResolvedSpectrumEnergyEv(baseline_spectrum, baseline_rng, 1.0);
    EXPECT_NEAR(baseline_sample, 15.0, 1.0e-12);

    sampler.markDirty();
    ASSERT_TRUE(sampler.dirty());

    SCDAT::Particle::TabulatedSurfaceDistributionFunction updated_distribution(
        {5.0, 15.0, 30.0}, {0.0, 0.0, 4.0}, {0.4, 0.7, 0.95});
    const auto updated_moments = updated_distribution.computeMoments();
    ASSERT_TRUE(updated_moments.valid);
    sampler.cacheMoments(updated_moments);

    EXPECT_FALSE(sampler.dirty());
    EXPECT_NEAR(sampler.cachedMoments().mean_energy_ev, 30.0, 1.0e-12);
    EXPECT_GT(sampler.cachedMoments().mean_energy_ev, baseline_moments.mean_energy_ev);

    SCDAT::Particle::ResolvedSpectrum updated_spectrum;
    updated_spectrum.model = SCDAT::Particle::SpatialSamplingModel::TABULATED;
    updated_spectrum.energy_grid_ev = {5.0, 15.0, 30.0};
    updated_spectrum.differential_number_flux = {0.0, 0.0, 4.0};

    std::mt19937 updated_rng(20260413u);
    const double updated_sample =
        SCDAT::Particle::sampleResolvedSpectrumEnergyEv(updated_spectrum, updated_rng, 1.0);
    EXPECT_NEAR(updated_sample, 30.0, 1.0e-12);
    EXPECT_NE(updated_sample, baseline_sample);
}

TEST(SurfaceDistributionFunctionTest, ResolvedSpectrumAverageMassKgUsesPopulationMoments)
{
    SCDAT::Particle::ResolvedSpectrum spectrum;
    spectrum.populations = {
        {2.0e11, 4.0, 0.0, 16.0, 1.0, 0.0},
        {1.0e11, 10.0, 0.0, 32.0, 1.0, 0.0},
    };

    constexpr double kAtomicMassUnit = 1.66053906660e-27;
    const double mass_kg =
        SCDAT::Particle::resolvedSpectrumAverageMassKg(spectrum, 28.0 * kAtomicMassUnit);
    EXPECT_NEAR(mass_kg, 21.333333333333332 * kAtomicMassUnit, 1.0e-39);
}
