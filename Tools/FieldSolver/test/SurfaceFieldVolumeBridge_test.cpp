#include "../include/SurfaceFieldVolumeBridge.h"

#include "../../Mesh/include/MeshParsing.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace
{

std::shared_ptr<SCDAT::Mesh::VolMesh> makeNativePlateMesh()
{
    SCDAT::Mesh::GeometryDimensions dims;
    dims.shape = SCDAT::Mesh::GeometryShape::PLATE;
    dims.length = 0.01;
    dims.width = 0.01;
    dims.thickness = 0.02;

    SCDAT::Mesh::MeshGenerationOptions options;
    options.nx = 1;
    options.ny = 1;
    options.nz = 6;
    options.tetrahedralize = true;

    auto mesh = SCDAT::Mesh::GeometryMeshGenerator::generateFromDimensions(dims, options);
    EXPECT_TRUE(mesh);
    EXPECT_TRUE(mesh->validate());
    return std::shared_ptr<SCDAT::Mesh::VolMesh>(std::move(mesh));
}

std::vector<SCDAT::Particle::Particle> makeNativeParityParticles()
{
    std::vector<SCDAT::Particle::Particle> particles;
    particles.emplace_back(1, SCDAT::Particle::ParticleType::ELECTRON,
                           SCDAT::Geometry::Point3D(0.002, 0.002, 0.004),
                           SCDAT::Geometry::Vector3D(1.0e5, 0.0, 0.0), 9.1093837015e-31,
                           -1.602176634e-19, 1.0);
    particles.emplace_back(2, SCDAT::Particle::ParticleType::ION,
                           SCDAT::Geometry::Point3D(0.004, 0.004, 0.012),
                           SCDAT::Geometry::Vector3D(-1.0e3, 0.0, 0.0), 1.67262192369e-27,
                           1.602176634e-19, 1.0);
    particles.emplace_back(3, SCDAT::Particle::ParticleType::ELECTRON,
                           SCDAT::Geometry::Point3D(0.006, 0.006, 0.016),
                           SCDAT::Geometry::Vector3D(0.0, 1.5e5, 0.0), 9.1093837015e-31,
                           -1.602176634e-19, 1.0);
    return particles;
}

std::vector<SCDAT::Particle::ParticlePtr> makeParticlePointers(
    std::vector<SCDAT::Particle::Particle>& particles)
{
    std::vector<SCDAT::Particle::ParticlePtr> pointers;
    pointers.reserve(particles.size());
    for (auto& particle : particles)
    {
        pointers.push_back(&particle);
    }
    return pointers;
}

} // namespace

TEST(SurfaceFieldVolumeBridgeTest, BlendFieldVolumeScalarClampsBlendFactor)
{
    EXPECT_NEAR(SCDAT::FieldSolver::blendFieldVolumeScalar(4.0, 10.0, -1.0), 4.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::FieldSolver::blendFieldVolumeScalar(4.0, 10.0, 0.5), 7.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::FieldSolver::blendFieldVolumeScalar(4.0, 10.0, 2.0), 10.0, 1.0e-12);
}

TEST(SurfaceFieldVolumeBridgeTest, BlendFieldVolumeScalarSanitizesNonFiniteInputs)
{
    const double value = SCDAT::FieldSolver::blendFieldVolumeScalar(
        std::numeric_limits<double>::quiet_NaN(), 8.0, 0.5);
    EXPECT_TRUE(std::isfinite(value));
    EXPECT_NEAR(value, 4.0, 1.0e-12);
}

TEST(SurfaceFieldVolumeBridgeTest, UpdateMismatchMetricUsesMaxNormalizedDeviation)
{
    double metric = 0.0;
    metric = SCDAT::FieldSolver::updateFieldVolumeMismatchMetric(metric, 5.0, 4.0, 2.0);
    EXPECT_NEAR(metric, 0.5, 1.0e-12);

    metric = SCDAT::FieldSolver::updateFieldVolumeMismatchMetric(metric, 9.0, 1.0, 2.0);
    EXPECT_NEAR(metric, 4.0, 1.0e-12);

    metric = SCDAT::FieldSolver::updateFieldVolumeMismatchMetric(
        metric, std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0);
    EXPECT_NEAR(metric, 4.0, 1.0e-12);
}

TEST(SurfaceFieldVolumeBridgeTest,
     ComputeExternalBlendFactorDecreasesWithMismatchAndNonConvergence)
{
    const double converged_small_mismatch =
        SCDAT::FieldSolver::computeFieldVolumeExternalBlendFactor(0.1, 1.0e-8, 1.0e-8, true,
                                                                   1.0e-6, 0.8);
    const double unconverged_large_mismatch =
        SCDAT::FieldSolver::computeFieldVolumeExternalBlendFactor(10.0, 1.0e-2, 1.0e-2, false,
                                                                   1.0e-6, 0.8);

    EXPECT_GE(converged_small_mismatch, 0.05);
    EXPECT_LE(converged_small_mismatch, 1.0);
    EXPECT_GE(unconverged_large_mismatch, 0.05);
    EXPECT_LE(unconverged_large_mismatch, 1.0);
    EXPECT_GT(converged_small_mismatch, unconverged_large_mismatch);
}

TEST(SurfaceFieldVolumeBridgeTest, NativeVolumeParitySnapshotBuildsMinimalClosedLoop)
{
    auto mesh = makeNativePlateMesh();
    ASSERT_TRUE(mesh);

    auto particles = makeNativeParityParticles();
    auto particle_ptrs = makeParticlePointers(particles);
    const auto snapshot = SCDAT::FieldSolver::buildNativeVolumeParitySnapshot(
        mesh, particle_ptrs, 10.0, 0.0, 1.0e-9,
        SCDAT::FieldSolver::NativeVolumeParityRoute::NativeMinimal);

    EXPECT_STREQ(snapshot.schema_version.c_str(), "scdat.native_volume_parity.v1");
    EXPECT_STREQ(SCDAT::FieldSolver::nativeVolumeParityRouteName(snapshot.resolved_route),
                 "native_minimal");
    EXPECT_TRUE(snapshot.mesh.valid);
    EXPECT_GT(snapshot.mesh.node_count, 0U);
    EXPECT_GT(snapshot.mesh.element_count, 0U);
    EXPECT_GT(snapshot.mesh.total_volume_m3, 0.0);
    EXPECT_TRUE(snapshot.field.solved);
    EXPECT_EQ(snapshot.field.boundary_condition_count, 2U);
    EXPECT_EQ(snapshot.field.dirichlet_boundary_count, 2U);
    EXPECT_EQ(snapshot.field.boundary_condition_family_count, 2U);
    EXPECT_EQ(snapshot.field.field_family_count, 3U);
    EXPECT_NE(snapshot.field.boundary_condition_family_signature.find("FourierPoissonBC"),
              std::string::npos);
    EXPECT_NE(snapshot.field.boundary_condition_family_signature.find("SurfDistribMatterBC"),
              std::string::npos);
    EXPECT_NE(snapshot.field.field_family_signature.find("UniformBField"),
              std::string::npos);
    EXPECT_NE(snapshot.field.field_family_signature.find("EMField"),
              std::string::npos);
    EXPECT_GT(snapshot.field.max_field_v_per_m, 0.0);
    EXPECT_GE(snapshot.field.max_potential_v, snapshot.field.min_potential_v);
    EXPECT_EQ(snapshot.distribution.active_particle_count, particle_ptrs.size());
    EXPECT_EQ(snapshot.distribution.electron_particle_count, 2U);
    EXPECT_EQ(snapshot.distribution.ion_particle_count, 1U);
    EXPECT_EQ(snapshot.distribution.distribution_channel_count, 3U);
    EXPECT_EQ(snapshot.distribution.family_count, 3U);
    EXPECT_NE(snapshot.distribution.family_signature.find("PICVolDistrib"),
              std::string::npos);
    EXPECT_NE(snapshot.distribution.family_signature.find("CompositeVolDistrib"),
              std::string::npos);
    EXPECT_GT(snapshot.distribution.average_electron_density_m3, 0.0);
    EXPECT_GT(snapshot.distribution.average_ion_density_m3, 0.0);
    EXPECT_EQ(snapshot.interaction.interactor_count, 2U);
    EXPECT_EQ(snapshot.interaction.executed_interaction_count, 2U);
    EXPECT_EQ(snapshot.interaction.family_count, 2U);
    EXPECT_NE(snapshot.interaction.family_signature.find("MCCInteractor"),
              std::string::npos);
    EXPECT_NE(snapshot.interaction.family_signature.find("TrajectoryInteractionFromField"),
              std::string::npos);
    EXPECT_TRUE(snapshot.interaction.executed);
}

TEST(SurfaceFieldVolumeBridgeTest, HybridBlendParityAddsPhase2BoundaryAndInteractionCoverage)
{
    auto mesh = makeNativePlateMesh();
    ASSERT_TRUE(mesh);

    auto particles = makeNativeParityParticles();
    auto particle_ptrs = makeParticlePointers(particles);
    const auto snapshot = SCDAT::FieldSolver::buildNativeVolumeParitySnapshot(
        mesh, particle_ptrs, 12.0, -2.0, 1.0e-9,
        SCDAT::FieldSolver::NativeVolumeParityRoute::HybridBlend);

    EXPECT_STREQ(SCDAT::FieldSolver::nativeVolumeParityRouteName(snapshot.resolved_route),
                 "hybrid_blend");
    EXPECT_TRUE(snapshot.field.solved);
    EXPECT_GE(snapshot.field.boundary_condition_count, 3U);
    EXPECT_EQ(snapshot.field.boundary_condition_family_count, 6U);
    EXPECT_EQ(snapshot.field.field_family_count, 5U);
    EXPECT_NE(snapshot.field.boundary_condition_family_signature.find("VoltageDependentMBC"),
              std::string::npos);
    EXPECT_NE(snapshot.field.boundary_condition_family_signature.find(
                  "MixedDirichletFourierPoissonBC"),
              std::string::npos);
    EXPECT_NE(snapshot.field.field_family_signature.find("SolenoidBField"),
              std::string::npos);
    EXPECT_NE(snapshot.field.field_family_signature.find("DipolarBField"),
              std::string::npos);
    EXPECT_GE(snapshot.distribution.distribution_channel_count, 3U);
    EXPECT_EQ(snapshot.distribution.family_count, 8U);
    EXPECT_NE(snapshot.distribution.family_signature.find("SmartPICVolDistrib"),
              std::string::npos);
    EXPECT_NE(snapshot.distribution.family_signature.find("BackTrackingVolDistrib"),
              std::string::npos);
    EXPECT_EQ(snapshot.interaction.interactor_count, 5U);
    EXPECT_EQ(snapshot.interaction.executed_interaction_count, 5U);
    EXPECT_EQ(snapshot.interaction.family_count, 5U);
    EXPECT_NE(snapshot.interaction.family_signature.find("CEXInteractor"),
              std::string::npos);
    EXPECT_NE(snapshot.interaction.family_signature.find("PhotoIonization"),
              std::string::npos);
    EXPECT_NE(snapshot.interaction.family_signature.find("SpinningSpacecraftTrajectory"),
              std::string::npos);
}

TEST(SurfaceFieldVolumeBridgeTest, ActiveBoundaryAndFieldFamilySignaturesFollowRoute)
{
    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeActiveBoundaryConditionFamilySignature(
            SCDAT::FieldSolver::NativeVolumeParityRoute::NativeMinimal),
        "FourierPoissonBC+SurfDistribMatterBC");
    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeActiveFieldFamilySignature(
            SCDAT::FieldSolver::NativeVolumeParityRoute::NativeMinimal),
        "UniformBField+MultipleEField+EMField");

    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeActiveBoundaryConditionFamilySignature(
            SCDAT::FieldSolver::NativeVolumeParityRoute::HybridBlend),
        "VoltageDependentMBC+MixedDirichletFourierPoissonBC+FourierPoissonBC+"
        "SurfDistribMatterBC+OneSurfDistribTestableMatterBC+CapacitiveVoltageGenerator");
    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeActiveFieldFamilySignature(
            SCDAT::FieldSolver::NativeVolumeParityRoute::HybridBlend),
        "UniformBField+SolenoidBField+DipolarBField+MultipleEField+EMField");
}

TEST(SurfaceFieldVolumeBridgeTest, ExplicitBoundaryAndFieldFamiliesNormalizeDuplicates)
{
    const auto boundary_families =
        SCDAT::FieldSolver::normalizeNativeVolumeBoundaryConditionFamilies(
            {SCDAT::FieldSolver::NativeVolumeBoundaryConditionFamily::VoltageDependentMBC,
             SCDAT::FieldSolver::NativeVolumeBoundaryConditionFamily::FourierPoissonBC,
             SCDAT::FieldSolver::NativeVolumeBoundaryConditionFamily::VoltageDependentMBC,
             SCDAT::FieldSolver::NativeVolumeBoundaryConditionFamily::CapacitiveVoltageGenerator});
    ASSERT_EQ(boundary_families.size(), 3U);
    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeBoundaryConditionFamilySignature(boundary_families),
        "VoltageDependentMBC+FourierPoissonBC+CapacitiveVoltageGenerator");

    const auto field_families = SCDAT::FieldSolver::normalizeNativeVolumeFieldFamilies(
        {SCDAT::FieldSolver::NativeVolumeFieldFamily::UniformBField,
         SCDAT::FieldSolver::NativeVolumeFieldFamily::DipolarBField,
         SCDAT::FieldSolver::NativeVolumeFieldFamily::UniformBField,
         SCDAT::FieldSolver::NativeVolumeFieldFamily::EMField});
    ASSERT_EQ(field_families.size(), 3U);
    EXPECT_EQ(SCDAT::FieldSolver::nativeVolumeFieldFamilySignature(field_families),
              "UniformBField+DipolarBField+EMField");
}

TEST(SurfaceFieldVolumeBridgeTest, ExplicitDistributionFamiliesNormalizeDuplicates)
{
    const auto distribution_families =
        SCDAT::FieldSolver::normalizeNativeVolumeDistributionFamilies(
            {SCDAT::FieldSolver::NativeVolumeDistributionFamily::GlobalMaxwellBoltzmannVolDistrib,
             SCDAT::FieldSolver::NativeVolumeDistributionFamily::PICBoltzmannVolDistrib,
             SCDAT::FieldSolver::NativeVolumeDistributionFamily::Updatable,
             SCDAT::FieldSolver::NativeVolumeDistributionFamily::GlobalMaxwellBoltzmannVolDistrib});
    ASSERT_EQ(distribution_families.size(), 3U);
    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeDistributionFamilySignature(distribution_families),
        "GlobalMaxwellBoltzmannVolDistrib+PICBoltzmannVolDistrib+Updatable");
}

TEST(SurfaceFieldVolumeBridgeTest, ExplicitInteractionFamiliesNormalizeDuplicates)
{
    const auto interaction_families =
        SCDAT::FieldSolver::normalizeNativeVolumeInteractionFamilies(
            {SCDAT::FieldSolver::NativeVolumeInteractionFamily::MCCInteractor,
             SCDAT::FieldSolver::NativeVolumeInteractionFamily::PhotoIonization,
             SCDAT::FieldSolver::NativeVolumeInteractionFamily::MCCInteractor,
             SCDAT::FieldSolver::NativeVolumeInteractionFamily::SpinningSpacecraftTrajectory});
    ASSERT_EQ(interaction_families.size(), 3U);
    EXPECT_EQ(
        SCDAT::FieldSolver::nativeVolumeInteractionFamilySignature(interaction_families),
        "MCCInteractor+PhotoIonization+SpinningSpacecraftTrajectory");
}
