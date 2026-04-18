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

TEST(SurfaceDistributionFunctionTest, MultipleSurfDistributionCombinesComponentMoments)
{
    SCDAT::Particle::SurfaceFluxMoments local_component;
    local_component.total_flux = 8.0;
    local_component.mean_energy_ev = 12.0;
    local_component.mean_cos_incidence = 0.45;
    local_component.valid = true;

    SCDAT::Particle::SurfaceFluxMoments remote_component;
    remote_component.total_flux = 4.0;
    remote_component.mean_energy_ev = 36.0;
    remote_component.mean_cos_incidence = 0.85;
    remote_component.valid = true;

    SCDAT::Particle::MultipleSurfSurfaceDistributionFunction distribution(
        {local_component, remote_component}, {1.0, 0.75}, 0.2);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_multiple_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_GT(moments.mean_energy_ev, local_component.mean_energy_ev);
    EXPECT_LT(moments.mean_energy_ev, remote_component.mean_energy_ev);
    EXPECT_GT(moments.mean_cos_incidence, local_component.mean_cos_incidence);
}

TEST(SurfaceDistributionFunctionTest, LocalModifiedPearsonDistributionProducesFiniteShapeAwareMoments)
{
    SCDAT::Particle::LocalModifiedPearsonIVSurfaceDistributionFunction distribution(
        6.5e10, 14.0, 0.8, 4.5, 0.3, 1.2);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_local_modified_pearson_iv_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_GT(moments.mean_energy_ev, 14.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
    EXPECT_LT(moments.mean_cos_incidence, 1.0);
}

TEST(SurfaceDistributionFunctionTest, LocalTabulatedDistributionAppliesLocalShiftAndScale)
{
    SCDAT::Particle::LocalTabulatedSurfaceDistributionFunction distribution(
        {5.0, 15.0, 35.0}, {1.0, 2.0, 1.0}, 4.0, 1.5, {0.25, 0.55, 0.85});
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_local_tabulated_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_NEAR(moments.total_flux, 6.0, 1.0e-12);
    EXPECT_GT(moments.mean_energy_ev, 15.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
}

TEST(SurfaceDistributionFunctionTest,
     TwoAxesTabulatedVelocityDistributionCombinesNormalAndTangentialEnergy)
{
    SCDAT::Particle::TwoAxesTabulatedVelocitySurfaceDistributionFunction distribution(
        {2.0e4, 4.0e4, 6.0e4}, {1.0e4, 2.0e4, 2.5e4}, {0.5, 1.0, 0.75}, 1.0e-10, 5.0e-11,
        {0.35, 0.6, 0.88});
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_two_axes_tabulated_velocity_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_GT(moments.mean_energy_ev, 0.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
}

TEST(SurfaceDistributionFunctionTest,
     AxisymTabulatedVelocityDistributionCombinesAxialAndRadialEnergy)
{
    SCDAT::Particle::AxisymTabulatedVelocitySurfaceDistributionFunction distribution(
        {2.0e4, 4.0e4, 6.0e4}, {0.8e4, 1.6e4, 2.4e4}, {0.5, 1.0, 0.75}, 1.0e-10, 6.0e-11,
        {0.35, 0.6, 0.88});
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_axisym_tabulated_velocity_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_GT(moments.mean_energy_ev, 0.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
}

TEST(SurfaceDistributionFunctionTest, UniformVelocityDistributionProducesSingleModeMoments)
{
    SCDAT::Particle::UniformVelocitySurfaceDistributionFunction distribution(
        5.0e4, 2.5, 1.0e-10, 0.8);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_uniform_velocity_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_NEAR(moments.total_flux, 2.5, 1.0e-12);
    EXPECT_GT(moments.mean_energy_ev, 0.0);
    EXPECT_NEAR(moments.mean_cos_incidence, 0.8, 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest, FluidDistributionTracksDirectedFlowAndCompressibility)
{
    SCDAT::Particle::FluidSurfaceDistributionFunction distribution(
        7.5e10, 9.0, 3500.0, 1.3, 0.2);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_fluid_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_GT(moments.mean_energy_ev, 9.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.0);
}

TEST(SurfaceDistributionFunctionTest, FowlerNordheimDistributionProducesFiniteEmissionMoments)
{
    SCDAT::Particle::FowlerNordheimSurfaceDistributionFunction distribution(
        7.5e7, 4.7, 1.4);
    const auto moments = distribution.computeMoments();

    EXPECT_STREQ(distribution.family(), "spis_fowler_nordheim_surface_flux_v1");
    EXPECT_TRUE(moments.valid);
    EXPECT_GT(moments.total_flux, 0.0);
    EXPECT_GT(moments.mean_energy_ev, 0.0);
    EXPECT_GT(moments.mean_cos_incidence, 0.9);
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

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportMultipleSurfAndLocalModifiedPearsonModels)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest multiple_request;
    multiple_request.model = SCDAT::Particle::SurfaceDistributionModel::MultipleSurf;
    multiple_request.electron_density_m3 = 1.1e11;
    multiple_request.electron_temperature_ev = 11.0;
    multiple_request.ion_density_m3 = 6.0e10;
    multiple_request.ion_temperature_ev = 4.0;
    multiple_request.ion_directed_velocity_m_per_s = 2500.0;
    multiple_request.electron_populations = {
        {7.0e10, 8.0},
        {3.0e10, 18.0},
        {1.0e10, 30.0},
    };
    multiple_request.ion_populations = {
        {4.0e10, 3.0},
        {2.0e10, 7.0},
    };
    const auto multiple_result = SCDAT::Particle::buildSurfaceFluxMoments(multiple_request);
    EXPECT_TRUE(multiple_result.valid);
    EXPECT_TRUE(multiple_result.electron_flux_moments.valid);
    EXPECT_TRUE(multiple_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest multiple_synthesis;
    multiple_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::MultipleSurf;
    multiple_synthesis.electron_flux_moments = multiple_result.electron_flux_moments;
    multiple_synthesis.ion_flux_moments = multiple_result.ion_flux_moments;
    multiple_synthesis.fallback_electron_energy_ev = 10.0;
    multiple_synthesis.fallback_ion_energy_ev = 4.0;
    multiple_synthesis.incidence_scale = 0.8;
    multiple_synthesis.local_field_factor = 1.15;
    multiple_synthesis.electron_population_count = multiple_request.electron_populations.size();
    multiple_synthesis.ion_population_count = multiple_request.ion_populations.size();
    const auto multiple_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(multiple_synthesis);
    EXPECT_TRUE(multiple_synthesized.valid);
    EXPECT_GT(multiple_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(multiple_synthesized.ion_collection_scale, 0.1);

    SCDAT::Particle::SurfaceDistributionBuildRequest pearson_request = multiple_request;
    pearson_request.model = SCDAT::Particle::SurfaceDistributionModel::LocalModifiedPearsonIV;
    const auto pearson_result = SCDAT::Particle::buildSurfaceFluxMoments(pearson_request);
    EXPECT_TRUE(pearson_result.valid);
    EXPECT_TRUE(pearson_result.electron_flux_moments.valid);
    EXPECT_TRUE(pearson_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest pearson_synthesis;
    pearson_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::LocalModifiedPearsonIV;
    pearson_synthesis.electron_flux_moments = pearson_result.electron_flux_moments;
    pearson_synthesis.ion_flux_moments = pearson_result.ion_flux_moments;
    pearson_synthesis.fallback_electron_energy_ev = 9.0;
    pearson_synthesis.fallback_ion_energy_ev = 3.5;
    pearson_synthesis.incidence_scale = 0.75;
    pearson_synthesis.local_field_factor = 1.1;
    pearson_synthesis.local_reference_shift_v = 2.0;
    pearson_synthesis.wake_factor = 0.25;
    const auto pearson_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(pearson_synthesis);
    EXPECT_TRUE(pearson_synthesized.valid);
    EXPECT_GT(pearson_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(pearson_synthesized.ion_collection_scale, 0.1);
    EXPECT_GE(pearson_synthesized.photo_incidence_scale, 0.0);
}

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportLocalTabulatedAndTwoAxesVelocityModels)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest local_tabulated_request;
    local_tabulated_request.model = SCDAT::Particle::SurfaceDistributionModel::LocalTabulated;
    local_tabulated_request.electron_density_m3 = 9.0e10;
    local_tabulated_request.electron_temperature_ev = 10.0;
    local_tabulated_request.ion_density_m3 = 5.5e10;
    local_tabulated_request.ion_temperature_ev = 3.5;
    local_tabulated_request.ion_directed_velocity_m_per_s = 1800.0;
    const auto local_tabulated_result =
        SCDAT::Particle::buildSurfaceFluxMoments(local_tabulated_request);
    EXPECT_TRUE(local_tabulated_result.valid);
    EXPECT_TRUE(local_tabulated_result.electron_flux_moments.valid);
    EXPECT_TRUE(local_tabulated_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest local_tabulated_synthesis;
    local_tabulated_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::LocalTabulated;
    local_tabulated_synthesis.electron_flux_moments = local_tabulated_result.electron_flux_moments;
    local_tabulated_synthesis.ion_flux_moments = local_tabulated_result.ion_flux_moments;
    local_tabulated_synthesis.fallback_electron_energy_ev = 8.0;
    local_tabulated_synthesis.fallback_ion_energy_ev = 3.0;
    local_tabulated_synthesis.incidence_scale = 0.82;
    local_tabulated_synthesis.local_field_factor = 1.12;
    const auto local_tabulated_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(local_tabulated_synthesis);
    EXPECT_TRUE(local_tabulated_synthesized.valid);
    EXPECT_GT(local_tabulated_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(local_tabulated_synthesized.ion_collection_scale, 0.1);

    SCDAT::Particle::SurfaceDistributionBuildRequest two_axes_request = local_tabulated_request;
    two_axes_request.model = SCDAT::Particle::SurfaceDistributionModel::TwoAxesTabulatedVelocity;
    two_axes_request.ion_directed_velocity_m_per_s = 3200.0;
    const auto two_axes_result =
        SCDAT::Particle::buildSurfaceFluxMoments(two_axes_request);
    EXPECT_TRUE(two_axes_result.valid);
    EXPECT_TRUE(two_axes_result.electron_flux_moments.valid);
    EXPECT_TRUE(two_axes_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest two_axes_synthesis;
    two_axes_synthesis.model =
        SCDAT::Particle::SurfaceDistributionModel::TwoAxesTabulatedVelocity;
    two_axes_synthesis.electron_flux_moments = two_axes_result.electron_flux_moments;
    two_axes_synthesis.ion_flux_moments = two_axes_result.ion_flux_moments;
    two_axes_synthesis.fallback_electron_energy_ev = 9.0;
    two_axes_synthesis.fallback_ion_energy_ev = 3.0;
    two_axes_synthesis.incidence_scale = 0.78;
    two_axes_synthesis.normalized_flow_alignment = 0.55;
    two_axes_synthesis.local_field_factor = 1.08;
    const auto two_axes_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(two_axes_synthesis);
    EXPECT_TRUE(two_axes_synthesized.valid);
    EXPECT_GT(two_axes_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(two_axes_synthesized.ion_collection_scale, 0.1);
    EXPECT_GE(two_axes_synthesized.photo_incidence_scale, 0.0);
}

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportUniformFluidAndFowlerNordheimModels)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest uniform_request;
    uniform_request.model = SCDAT::Particle::SurfaceDistributionModel::UniformVelocity;
    uniform_request.electron_density_m3 = 8.0e10;
    uniform_request.electron_temperature_ev = 8.5;
    uniform_request.ion_density_m3 = 4.5e10;
    uniform_request.ion_temperature_ev = 3.0;
    uniform_request.ion_directed_velocity_m_per_s = 2200.0;
    const auto uniform_result = SCDAT::Particle::buildSurfaceFluxMoments(uniform_request);
    EXPECT_TRUE(uniform_result.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest uniform_synthesis;
    uniform_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::UniformVelocity;
    uniform_synthesis.electron_flux_moments = uniform_result.electron_flux_moments;
    uniform_synthesis.ion_flux_moments = uniform_result.ion_flux_moments;
    uniform_synthesis.fallback_electron_energy_ev = 7.5;
    uniform_synthesis.fallback_ion_energy_ev = 2.5;
    uniform_synthesis.local_field_factor = 1.05;
    uniform_synthesis.incidence_scale = 0.85;
    const auto uniform_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(uniform_synthesis);
    EXPECT_TRUE(uniform_synthesized.valid);
    EXPECT_GT(uniform_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(uniform_synthesized.ion_collection_scale, 0.1);

    SCDAT::Particle::SurfaceDistributionBuildRequest fluid_request = uniform_request;
    fluid_request.model = SCDAT::Particle::SurfaceDistributionModel::Fluid;
    fluid_request.ion_directed_velocity_m_per_s = 4200.0;
    const auto fluid_result = SCDAT::Particle::buildSurfaceFluxMoments(fluid_request);
    EXPECT_TRUE(fluid_result.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest fluid_synthesis;
    fluid_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::Fluid;
    fluid_synthesis.electron_flux_moments = fluid_result.electron_flux_moments;
    fluid_synthesis.ion_flux_moments = fluid_result.ion_flux_moments;
    fluid_synthesis.fallback_electron_energy_ev = 7.5;
    fluid_synthesis.fallback_ion_energy_ev = 2.5;
    fluid_synthesis.local_field_factor = 1.1;
    fluid_synthesis.incidence_scale = 0.8;
    fluid_synthesis.normalized_flow_alignment = 0.6;
    const auto fluid_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(fluid_synthesis);
    EXPECT_TRUE(fluid_synthesized.valid);
    EXPECT_GT(fluid_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(fluid_synthesized.ion_collection_scale, 0.1);

    SCDAT::Particle::SurfaceDistributionBuildRequest fn_request = uniform_request;
    fn_request.model = SCDAT::Particle::SurfaceDistributionModel::FowlerNordheim;
    const auto fn_result = SCDAT::Particle::buildSurfaceFluxMoments(fn_request);
    EXPECT_TRUE(fn_result.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest fn_synthesis;
    fn_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::FowlerNordheim;
    fn_synthesis.electron_flux_moments = fn_result.electron_flux_moments;
    fn_synthesis.ion_flux_moments = fn_result.ion_flux_moments;
    fn_synthesis.fallback_electron_energy_ev = 1.0;
    fn_synthesis.fallback_ion_energy_ev = 2.5;
    fn_synthesis.local_field_factor = 1.2;
    fn_synthesis.incidence_scale = 0.7;
    const auto fn_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(fn_synthesis);
    EXPECT_TRUE(fn_synthesized.valid);
    EXPECT_GT(fn_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(fn_synthesized.ion_collection_scale, 0.1);
}

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportAxisymTabulatedVelocityModel)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest axisym_request;
    axisym_request.model = SCDAT::Particle::SurfaceDistributionModel::AxisymTabulatedVelocity;
    axisym_request.electron_density_m3 = 9.5e10;
    axisym_request.electron_temperature_ev = 9.0;
    axisym_request.ion_density_m3 = 5.0e10;
    axisym_request.ion_temperature_ev = 3.2;
    axisym_request.ion_directed_velocity_m_per_s = 2800.0;
    axisym_request.axisym_anisotropy = 2.5;

    const auto axisym_result = SCDAT::Particle::buildSurfaceFluxMoments(axisym_request);
    EXPECT_TRUE(axisym_result.valid);
    EXPECT_TRUE(axisym_result.electron_flux_moments.valid);
    EXPECT_TRUE(axisym_result.ion_flux_moments.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest axisym_synthesis;
    axisym_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::AxisymTabulatedVelocity;
    axisym_synthesis.electron_flux_moments = axisym_result.electron_flux_moments;
    axisym_synthesis.ion_flux_moments = axisym_result.ion_flux_moments;
    axisym_synthesis.fallback_electron_energy_ev = 8.0;
    axisym_synthesis.fallback_ion_energy_ev = 3.0;
    axisym_synthesis.incidence_scale = 0.8;
    axisym_synthesis.normalized_flow_alignment = 0.5;
    axisym_synthesis.axisym_anisotropy = axisym_request.axisym_anisotropy;
    axisym_synthesis.local_field_factor = 1.1;

    const auto axisym_synthesized =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(axisym_synthesis);
    EXPECT_TRUE(axisym_synthesized.valid);
    EXPECT_GT(axisym_synthesized.electron_collection_scale, 0.1);
    EXPECT_GT(axisym_synthesized.ion_collection_scale, 0.1);
    EXPECT_GE(axisym_synthesized.photo_incidence_scale, 0.0);
}

TEST(SurfaceDistributionFunctionTest, BuildAndSynthesisSupportMaxwellFamilyModels)
{
    const std::vector<SCDAT::Particle::SurfaceDistributionModel> models{
        SCDAT::Particle::SurfaceDistributionModel::GlobalMaxwellBoltzmann,
        SCDAT::Particle::SurfaceDistributionModel::GlobalMaxwellBoltzmann2,
        SCDAT::Particle::SurfaceDistributionModel::GlobalMaxwell,
        SCDAT::Particle::SurfaceDistributionModel::LocalMaxwell,
        SCDAT::Particle::SurfaceDistributionModel::LocalMaxwell2,
    };

    for (const auto model : models)
    {
        SCOPED_TRACE(static_cast<int>(model));

        SCDAT::Particle::SurfaceDistributionBuildRequest request;
        request.model = model;
        request.electron_density_m3 = 1.0e11;
        request.electron_temperature_ev = 10.0;
        request.ion_density_m3 = 5.0e10;
        request.ion_temperature_ev = 3.5;
        request.ion_directed_velocity_m_per_s = 3200.0;

        const auto built = SCDAT::Particle::buildSurfaceFluxMoments(request);
        EXPECT_TRUE(built.valid);
        EXPECT_TRUE(built.electron_flux_moments.valid);
        EXPECT_TRUE(built.ion_flux_moments.valid);
        EXPECT_GT(built.electron_flux_moments.total_flux, 0.0);
        EXPECT_GT(built.ion_flux_moments.total_flux, 0.0);

        SCDAT::Particle::SurfaceDistributionSynthesisRequest synthesis;
        synthesis.model = model;
        synthesis.electron_flux_moments = built.electron_flux_moments;
        synthesis.ion_flux_moments = built.ion_flux_moments;
        synthesis.fallback_electron_energy_ev = 8.0;
        synthesis.fallback_ion_energy_ev = 3.0;
        synthesis.local_reference_shift_v = 2.0;
        synthesis.incidence_scale = 0.82;
        synthesis.normalized_flow_alignment = 0.55;
        synthesis.local_field_factor = 1.1;
        synthesis.axisym_anisotropy = 1.8;
        synthesis.patch_role = true;

        const auto result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(synthesis);
        EXPECT_TRUE(result.valid);
        EXPECT_GT(result.electron_collection_scale, 0.1);
        EXPECT_GT(result.ion_collection_scale, 0.1);
        EXPECT_GT(result.electron_characteristic_energy_ev, 0.0);
        EXPECT_GT(result.ion_characteristic_energy_ev, 0.0);
        EXPECT_GE(result.photo_incidence_scale, 0.0);
    }
}

TEST(SurfaceDistributionFunctionTest,
     RecollMaxwellScalesLocalMaxwellFluxAndCollectionByRecollectionRatio)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest local_request;
    local_request.model = SCDAT::Particle::SurfaceDistributionModel::LocalMaxwell;
    local_request.electron_density_m3 = 1.0e11;
    local_request.electron_temperature_ev = 10.0;
    local_request.ion_density_m3 = 5.0e10;
    local_request.ion_temperature_ev = 3.5;
    local_request.ion_directed_velocity_m_per_s = 3200.0;

    auto recoll_request = local_request;
    recoll_request.model = SCDAT::Particle::SurfaceDistributionModel::RecollMaxwell;
    recoll_request.electron_recollection_ratio = 0.35;
    recoll_request.ion_recollection_ratio = 0.60;

    const auto local_built = SCDAT::Particle::buildSurfaceFluxMoments(local_request);
    const auto recoll_built = SCDAT::Particle::buildSurfaceFluxMoments(recoll_request);
    ASSERT_TRUE(local_built.valid);
    ASSERT_TRUE(recoll_built.valid);
    EXPECT_NEAR(recoll_built.electron_flux_moments.total_flux,
                local_built.electron_flux_moments.total_flux * 0.35,
                local_built.electron_flux_moments.total_flux * 1.0e-12 + 1.0e-12);
    EXPECT_NEAR(recoll_built.ion_flux_moments.total_flux,
                local_built.ion_flux_moments.total_flux * 0.60,
                local_built.ion_flux_moments.total_flux * 1.0e-12 + 1.0e-12);
    EXPECT_NEAR(recoll_built.electron_flux_moments.mean_energy_ev,
                local_built.electron_flux_moments.mean_energy_ev, 1.0e-12);
    EXPECT_NEAR(recoll_built.ion_flux_moments.mean_energy_ev,
                local_built.ion_flux_moments.mean_energy_ev, 1.0e-12);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest local_synthesis;
    local_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::LocalMaxwell;
    local_synthesis.electron_flux_moments = local_built.electron_flux_moments;
    local_synthesis.ion_flux_moments = local_built.ion_flux_moments;
    local_synthesis.fallback_electron_energy_ev = 8.0;
    local_synthesis.fallback_ion_energy_ev = 3.0;
    local_synthesis.local_reference_shift_v = 2.0;
    local_synthesis.incidence_scale = 0.82;
    local_synthesis.normalized_flow_alignment = 0.55;
    local_synthesis.local_field_factor = 1.1;
    local_synthesis.patch_role = true;

    auto recoll_synthesis = local_synthesis;
    recoll_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::RecollMaxwell;
    recoll_synthesis.electron_flux_moments = recoll_built.electron_flux_moments;
    recoll_synthesis.ion_flux_moments = recoll_built.ion_flux_moments;
    recoll_synthesis.electron_recollection_ratio = 0.35;
    recoll_synthesis.ion_recollection_ratio = 0.60;

    const auto local_result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(local_synthesis);
    const auto recoll_result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(recoll_synthesis);
    ASSERT_TRUE(local_result.valid);
    ASSERT_TRUE(recoll_result.valid);
    EXPECT_NEAR(recoll_result.electron_collection_scale,
                local_result.electron_collection_scale * 0.35,
                local_result.electron_collection_scale * 1.0e-12 + 1.0e-12);
    EXPECT_NEAR(recoll_result.ion_collection_scale,
                local_result.ion_collection_scale * 0.60,
                local_result.ion_collection_scale * 1.0e-12 + 1.0e-12);
    EXPECT_GT(recoll_result.photo_incidence_scale, 0.0);
    EXPECT_LT(recoll_result.photo_incidence_scale, local_synthesis.incidence_scale + 1.0e-12);
}

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportPicAndGenericSurfaceDistributionFamilies)
{
    const std::vector<SCDAT::Particle::SurfaceDistributionModel> models{
        SCDAT::Particle::SurfaceDistributionModel::PICSurf,
        SCDAT::Particle::SurfaceDistributionModel::NonPICSurf,
        SCDAT::Particle::SurfaceDistributionModel::GenericSurf,
        SCDAT::Particle::SurfaceDistributionModel::GlobalSurf,
        SCDAT::Particle::SurfaceDistributionModel::LocalGenericSurf,
    };

    for (const auto model : models)
    {
        SCOPED_TRACE(static_cast<int>(model));

        SCDAT::Particle::SurfaceDistributionBuildRequest request;
        request.model = model;
        request.electron_density_m3 = 1.2e11;
        request.electron_temperature_ev = 12.0;
        request.ion_density_m3 = 6.5e10;
        request.ion_temperature_ev = 4.5;
        request.ion_directed_velocity_m_per_s = 4200.0;

        const auto built = SCDAT::Particle::buildSurfaceFluxMoments(request);
        ASSERT_TRUE(built.valid);
        EXPECT_TRUE(built.electron_flux_moments.valid);
        EXPECT_TRUE(built.ion_flux_moments.valid);
        EXPECT_GT(built.electron_flux_moments.total_flux, 0.0);
        EXPECT_GT(built.ion_flux_moments.total_flux, 0.0);

        SCDAT::Particle::SurfaceDistributionSynthesisRequest synthesis;
        synthesis.model = model;
        synthesis.electron_flux_moments = built.electron_flux_moments;
        synthesis.ion_flux_moments = built.ion_flux_moments;
        synthesis.fallback_electron_energy_ev = 9.0;
        synthesis.fallback_ion_energy_ev = 4.0;
        synthesis.local_reference_shift_v = 2.5;
        synthesis.incidence_scale = 0.80;
        synthesis.normalized_flow_alignment = 0.35;
        synthesis.local_field_factor = 1.15;
        synthesis.wake_factor = 0.20;
        synthesis.patch_role = true;

        const auto result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(synthesis);
        EXPECT_TRUE(result.valid);
        EXPECT_GT(result.electron_collection_scale, 0.1);
        EXPECT_GT(result.ion_collection_scale, 0.1);
        EXPECT_GT(result.electron_characteristic_energy_ev, 0.0);
        EXPECT_GT(result.ion_characteristic_energy_ev, 0.0);
        EXPECT_GE(result.photo_incidence_scale, 0.0);
    }
}

TEST(SurfaceDistributionFunctionTest, PicSurfaceFamilyExceedsNonPicElectronCollectionStrength)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest build_request;
    build_request.electron_density_m3 = 1.0e11;
    build_request.electron_temperature_ev = 10.0;
    build_request.ion_density_m3 = 5.0e10;
    build_request.ion_temperature_ev = 3.0;
    build_request.ion_directed_velocity_m_per_s = 3600.0;

    auto pic_request = build_request;
    pic_request.model = SCDAT::Particle::SurfaceDistributionModel::PICSurf;
    auto nonpic_request = build_request;
    nonpic_request.model = SCDAT::Particle::SurfaceDistributionModel::NonPICSurf;

    const auto pic_built = SCDAT::Particle::buildSurfaceFluxMoments(pic_request);
    const auto nonpic_built = SCDAT::Particle::buildSurfaceFluxMoments(nonpic_request);
    ASSERT_TRUE(pic_built.valid);
    ASSERT_TRUE(nonpic_built.valid);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest pic_synthesis;
    pic_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::PICSurf;
    pic_synthesis.electron_flux_moments = pic_built.electron_flux_moments;
    pic_synthesis.ion_flux_moments = pic_built.ion_flux_moments;
    pic_synthesis.fallback_electron_energy_ev = 8.0;
    pic_synthesis.fallback_ion_energy_ev = 3.0;
    pic_synthesis.local_reference_shift_v = 2.0;
    pic_synthesis.incidence_scale = 0.82;
    pic_synthesis.normalized_flow_alignment = 0.4;
    pic_synthesis.local_field_factor = 1.2;
    pic_synthesis.wake_factor = 0.25;
    pic_synthesis.patch_role = true;

    auto nonpic_synthesis = pic_synthesis;
    nonpic_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::NonPICSurf;
    nonpic_synthesis.electron_flux_moments = nonpic_built.electron_flux_moments;
    nonpic_synthesis.ion_flux_moments = nonpic_built.ion_flux_moments;

    const auto pic_result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(pic_synthesis);
    const auto nonpic_result =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(nonpic_synthesis);
    ASSERT_TRUE(pic_result.valid);
    ASSERT_TRUE(nonpic_result.valid);
    EXPECT_GT(pic_result.electron_collection_scale, nonpic_result.electron_collection_scale);
    EXPECT_GT(pic_result.photo_incidence_scale, 0.0);
    EXPECT_GT(nonpic_result.photo_incidence_scale, 0.0);
}

TEST(SurfaceDistributionFunctionTest,
     BuildAndSynthesisSupportTestableAndThrusterSurfaceDistributionFamilies)
{
    const std::vector<SCDAT::Particle::SurfaceDistributionModel> models{
        SCDAT::Particle::SurfaceDistributionModel::TestableSurf,
        SCDAT::Particle::SurfaceDistributionModel::TestableForA,
        SCDAT::Particle::SurfaceDistributionModel::MaxwellianThruster,
    };

    for (const auto model : models)
    {
        SCOPED_TRACE(static_cast<int>(model));

        SCDAT::Particle::SurfaceDistributionBuildRequest request;
        request.model = model;
        request.electron_density_m3 = 1.1e11;
        request.electron_temperature_ev = 11.0;
        request.ion_density_m3 = 6.2e10;
        request.ion_temperature_ev = 4.2;
        request.ion_directed_velocity_m_per_s = 4500.0;

        const auto built = SCDAT::Particle::buildSurfaceFluxMoments(request);
        ASSERT_TRUE(built.valid);
        EXPECT_TRUE(built.electron_flux_moments.valid);
        EXPECT_TRUE(built.ion_flux_moments.valid);
        EXPECT_GT(built.electron_flux_moments.total_flux, 0.0);
        EXPECT_GT(built.ion_flux_moments.total_flux, 0.0);

        SCDAT::Particle::SurfaceDistributionSynthesisRequest synthesis;
        synthesis.model = model;
        synthesis.electron_flux_moments = built.electron_flux_moments;
        synthesis.ion_flux_moments = built.ion_flux_moments;
        synthesis.fallback_electron_energy_ev = 9.0;
        synthesis.fallback_ion_energy_ev = 4.0;
        synthesis.local_reference_shift_v = 2.5;
        synthesis.incidence_scale = 0.78;
        synthesis.normalized_flow_alignment = 0.62;
        synthesis.local_field_factor = 1.18;
        synthesis.wake_factor = 0.12;
        synthesis.patch_role = true;

        const auto result = SCDAT::Particle::synthesizeSurfaceDistributionMoments(synthesis);
        EXPECT_TRUE(result.valid);
        EXPECT_GT(result.electron_collection_scale, 0.1);
        EXPECT_GT(result.ion_collection_scale, 0.1);
        EXPECT_GT(result.electron_characteristic_energy_ev, 0.0);
        EXPECT_GT(result.ion_characteristic_energy_ev, 0.0);
        EXPECT_GE(result.photo_incidence_scale, 0.0);
    }
}

TEST(SurfaceDistributionFunctionTest,
     TestableAndThrusterFamiliesProduceDistinctFluxAndCollectionResponses)
{
    SCDAT::Particle::SurfaceDistributionBuildRequest build_request;
    build_request.electron_density_m3 = 1.0e11;
    build_request.electron_temperature_ev = 10.0;
    build_request.ion_density_m3 = 5.0e10;
    build_request.ion_temperature_ev = 3.5;
    build_request.ion_directed_velocity_m_per_s = 5200.0;

    auto testable_request = build_request;
    testable_request.model = SCDAT::Particle::SurfaceDistributionModel::TestableSurf;
    auto testable_for_a_request = build_request;
    testable_for_a_request.model = SCDAT::Particle::SurfaceDistributionModel::TestableForA;
    auto thruster_request = build_request;
    thruster_request.model = SCDAT::Particle::SurfaceDistributionModel::MaxwellianThruster;
    auto local_maxwell_request = build_request;
    local_maxwell_request.model = SCDAT::Particle::SurfaceDistributionModel::LocalMaxwell;

    const auto testable_built = SCDAT::Particle::buildSurfaceFluxMoments(testable_request);
    const auto testable_for_a_built =
        SCDAT::Particle::buildSurfaceFluxMoments(testable_for_a_request);
    const auto thruster_built = SCDAT::Particle::buildSurfaceFluxMoments(thruster_request);
    const auto local_maxwell_built =
        SCDAT::Particle::buildSurfaceFluxMoments(local_maxwell_request);
    ASSERT_TRUE(testable_built.valid);
    ASSERT_TRUE(testable_for_a_built.valid);
    ASSERT_TRUE(thruster_built.valid);
    ASSERT_TRUE(local_maxwell_built.valid);

    EXPECT_NE(testable_built.electron_flux_moments.mean_energy_ev,
              testable_for_a_built.electron_flux_moments.mean_energy_ev);
    EXPECT_GT(thruster_built.ion_flux_moments.total_flux,
              local_maxwell_built.ion_flux_moments.total_flux);

    SCDAT::Particle::SurfaceDistributionSynthesisRequest testable_synthesis;
    testable_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::TestableSurf;
    testable_synthesis.electron_flux_moments = testable_built.electron_flux_moments;
    testable_synthesis.ion_flux_moments = testable_built.ion_flux_moments;
    testable_synthesis.fallback_electron_energy_ev = 8.0;
    testable_synthesis.fallback_ion_energy_ev = 3.0;
    testable_synthesis.local_reference_shift_v = 2.0;
    testable_synthesis.incidence_scale = 0.82;
    testable_synthesis.normalized_flow_alignment = 0.55;
    testable_synthesis.local_field_factor = 1.12;
    testable_synthesis.patch_role = true;

    auto testable_for_a_synthesis = testable_synthesis;
    testable_for_a_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::TestableForA;
    testable_for_a_synthesis.electron_flux_moments = testable_for_a_built.electron_flux_moments;
    testable_for_a_synthesis.ion_flux_moments = testable_for_a_built.ion_flux_moments;

    auto thruster_synthesis = testable_synthesis;
    thruster_synthesis.model = SCDAT::Particle::SurfaceDistributionModel::MaxwellianThruster;
    thruster_synthesis.electron_flux_moments = thruster_built.electron_flux_moments;
    thruster_synthesis.ion_flux_moments = thruster_built.ion_flux_moments;

    const auto testable_result =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(testable_synthesis);
    const auto testable_for_a_result =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(testable_for_a_synthesis);
    const auto thruster_result =
        SCDAT::Particle::synthesizeSurfaceDistributionMoments(thruster_synthesis);
    ASSERT_TRUE(testable_result.valid);
    ASSERT_TRUE(testable_for_a_result.valid);
    ASSERT_TRUE(thruster_result.valid);

    EXPECT_NE(testable_result.electron_collection_scale,
              testable_for_a_result.electron_collection_scale);
    EXPECT_GT(thruster_result.ion_collection_scale, testable_result.ion_collection_scale);
    EXPECT_LT(thruster_result.photo_incidence_scale, testable_result.photo_incidence_scale);
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
