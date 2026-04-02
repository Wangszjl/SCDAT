#include "../include/PICCycle.h"

#include "../../Boundary/include/BoundaryConditions.h"
#include "../../FieldSolver/include/PoissonSolver.h"
#include "../../Mesh/include/MeshParsing.h"
#include "../../Particle/include/ParticleManager.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace
{

using SCDAT::FieldSolver::PoissonSolver;
using SCDAT::Geometry::Point3D;
using SCDAT::Geometry::Vector3D;
using SCDAT::Mesh::BoundaryType;
using SCDAT::Mesh::GeometryDimensions;
using SCDAT::Mesh::GeometryMeshGenerator;
using SCDAT::Mesh::GeometryShape;
using SCDAT::Mesh::MeshGenerationOptions;
using SCDAT::Mesh::VolMesh;
using SCDAT::Particle::AbsorbingBoundaryCondition;
using SCDAT::Particle::ParticleManager;
using SCDAT::Particle::PeriodicBoundaryCondition;
using SCDAT::PICcore::PICCycle;

std::shared_ptr<VolMesh> createPlateMesh()
{
    GeometryDimensions dims;
    dims.shape = GeometryShape::PLATE;
    dims.length = 1.0;
    dims.width = 1.0;
    dims.thickness = 1.0;

    MeshGenerationOptions options;
    options.nx = 1;
    options.ny = 1;
    options.nz = 6;
    options.tetrahedralize = true;

    auto mesh = GeometryMeshGenerator::generateFromDimensions(dims, options);
    EXPECT_TRUE(mesh);
    EXPECT_TRUE(mesh->validate());

    const double z_min = mesh->getBoundingBoxMin().z();
    const double z_max = mesh->getBoundingBoxMax().z();
    for (const auto& node : mesh->getNodes())
    {
        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12 || std::abs(z - z_max) < 1.0e-12)
        {
            node->setBoundaryType(BoundaryType::DIRICHLET);
        }
        else
        {
            node->setBoundaryType(BoundaryType::INTERIOR);
        }
    }

    return std::shared_ptr<VolMesh>(std::move(mesh));
}

PICCycle createConfiguredCycle(const std::shared_ptr<VolMesh>& mesh,
                               const std::shared_ptr<ParticleManager>& particle_manager)
{
    auto poisson_solver = std::make_shared<PoissonSolver>(mesh);

    PICCycle::Parameters params;
    params.time_step = 1.0e-9;
    params.max_iterations = 4;
    params.statistics_frequency = 1;

    PICCycle cycle(params);
    cycle.setMesh(mesh->getNodes(), mesh->getElements());
    cycle.setParticleManager(particle_manager);
    cycle.setPoissonSolver(poisson_solver);

    const Point3D min_corner = mesh->getBoundingBoxMin();
    const Point3D max_corner = mesh->getBoundingBoxMax();
    auto periodic_boundary = std::make_shared<PeriodicBoundaryCondition>(
        std::vector<Vector3D>{Vector3D(max_corner.x() - min_corner.x(), 0.0, 0.0),
                              Vector3D(0.0, max_corner.y() - min_corner.y(), 0.0)});
    cycle.setBoundaryCondition(PICCycle::BoundaryFace::XMin, periodic_boundary);
    cycle.setBoundaryCondition(PICCycle::BoundaryFace::XMax, periodic_boundary);
    cycle.setBoundaryCondition(PICCycle::BoundaryFace::YMin, periodic_boundary);
    cycle.setBoundaryCondition(PICCycle::BoundaryFace::YMax, periodic_boundary);
    cycle.setBoundaryCondition(PICCycle::BoundaryFace::ZMin,
                               std::make_shared<AbsorbingBoundaryCondition>(1.0));
    cycle.setBoundaryCondition(PICCycle::BoundaryFace::ZMax,
                               std::make_shared<AbsorbingBoundaryCondition>(1.0));

    EXPECT_TRUE(cycle.initialize());
    cycle.setPotential(std::vector<double>(mesh->getNodeCount(), 0.0));
    return cycle;
}

std::size_t findElementWithLargestCoordinate(const VolMesh& mesh, int axis)
{
    std::size_t best_index = 0;
    double best_value = -1.0e300;
    for (std::size_t i = 0; i < mesh.getElementCount(); ++i)
    {
        const auto centroid = mesh.getElement(i)->getCentroid();
        const double value = (axis == 0) ? centroid.x() : (axis == 1 ? centroid.y() : centroid.z());
        if (value > best_value)
        {
            best_value = value;
            best_index = i;
        }
    }
    return best_index;
}

} // namespace

TEST(PICCycleBoundaryTest, PeriodicBoundaryKeepsParticleActiveInsideDomain)
{
    auto mesh = createPlateMesh();
    auto particle_manager = std::make_shared<ParticleManager>();
    PICCycle cycle = createConfiguredCycle(mesh, particle_manager);

    const std::size_t boundary_element = findElementWithLargestCoordinate(*mesh, 0);
    const Point3D start = mesh->getElement(boundary_element)->getCentroid();
    const auto particle_id =
        particle_manager->createElectron(start, Vector3D(2.0e9, 0.0, 0.0), 1.0);

    ASSERT_TRUE(cycle.executeTimeStep());

    const auto* particle = particle_manager->getParticle(particle_id);
    ASSERT_NE(particle, nullptr);
    EXPECT_TRUE(particle->isActive());

    const Point3D min_corner = mesh->getBoundingBoxMin();
    const Point3D max_corner = mesh->getBoundingBoxMax();
    EXPECT_GE(particle->getPosition().x(), min_corner.x() - 1.0e-9);
    EXPECT_LE(particle->getPosition().x(), max_corner.x() + 1.0e-9);
}

TEST(PICCycleBoundaryTest, AbsorbingBoundaryMarksParticleInactiveAndAccumulatesCharge)
{
    auto mesh = createPlateMesh();
    auto particle_manager = std::make_shared<ParticleManager>();
    PICCycle cycle = createConfiguredCycle(mesh, particle_manager);

    const std::size_t top_element = findElementWithLargestCoordinate(*mesh, 2);
    const Point3D start = mesh->getElement(top_element)->getCentroid();
    const auto particle_id =
        particle_manager->createElectron(start, Vector3D(0.0, 0.0, 2.0e9), 1.0);

    ASSERT_TRUE(cycle.executeTimeStep());

    const auto* particle = particle_manager->getParticle(particle_id);
    ASSERT_NE(particle, nullptr);
    EXPECT_FALSE(particle->isActive());
    EXPECT_LT(cycle.getStatistics().surface_deposited_charge, 0.0);
    const auto& ledger = cycle.getSurfaceCurrentLedger();
    EXPECT_LT(ledger[static_cast<std::size_t>(PICCycle::BoundaryFace::ZMax)]
                  .absorbed_electron_charge_c,
              0.0);
    EXPECT_EQ(ledger[static_cast<std::size_t>(PICCycle::BoundaryFace::ZMax)]
                  .absorbed_electron_particles,
              1u);
}

TEST(PICCycleBoundaryTest, ManualPotentialProfileReconstructsElectricField)
{
    auto mesh = createPlateMesh();
    auto particle_manager = std::make_shared<ParticleManager>();
    PICCycle cycle = createConfiguredCycle(mesh, particle_manager);

    std::vector<double> potential(mesh->getNodeCount(), 0.0);
    for (const auto& node : mesh->getNodes())
    {
        potential[node->getId()] = 12.0 * node->getPosition().z();
    }

    cycle.setPotential(potential);

    const auto& electric_field = cycle.getElectricField();
    ASSERT_EQ(electric_field.size(), mesh->getNodeCount());

    double mean_ez = 0.0;
    std::size_t finite_samples = 0;
    for (const auto& field : electric_field)
    {
        EXPECT_TRUE(std::isfinite(field.x()));
        EXPECT_TRUE(std::isfinite(field.y()));
        EXPECT_TRUE(std::isfinite(field.z()));
        mean_ez += field.z();
        ++finite_samples;
    }

    ASSERT_GT(finite_samples, 0u);
    mean_ez /= static_cast<double>(finite_samples);
    EXPECT_LT(mean_ez, -1.0);
}

TEST(PICCycleBoundaryTest, LegacyFallbackDoesNotRecoverFarEscapes)
{
    auto mesh = createPlateMesh();
    auto particle_manager = std::make_shared<ParticleManager>();
    PICCycle cycle = createConfiguredCycle(mesh, particle_manager);
    cycle.clearBoundaryConditions();

    const std::size_t top_element = findElementWithLargestCoordinate(*mesh, 2);
    const Point3D start = mesh->getElement(top_element)->getCentroid();
    const auto particle_id =
        particle_manager->createElectron(start, Vector3D(0.0, 0.0, 5.0e12), 1.0);

    ASSERT_TRUE(cycle.executeTimeStep());

    const auto* particle = particle_manager->getParticle(particle_id);
    ASSERT_NE(particle, nullptr);
    EXPECT_FALSE(particle->isActive());
    EXPECT_NEAR(cycle.getStatistics().surface_deposited_charge, 0.0, 1.0e-20);
}
