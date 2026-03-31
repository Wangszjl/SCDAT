#include "../include/PoissonSolver.h"
#include "../../Mesh/include/MeshParsing.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>

namespace
{

using SCDAT::FieldSolver::BoundaryCondition;
using SCDAT::FieldSolver::BoundaryConditionType;
using SCDAT::FieldSolver::PoissonSolver;
using SCDAT::Geometry::Point3D;
using SCDAT::Mesh::BoundaryType;
using SCDAT::Mesh::GeometryDimensions;
using SCDAT::Mesh::GeometryMeshGenerator;
using SCDAT::Mesh::GeometryShape;
using SCDAT::Mesh::MeshGenerationOptions;
using SCDAT::Mesh::VolMesh;
using SCDAT::Solver::SolverType;

std::shared_ptr<VolMesh> makePlateMesh()
{
    GeometryDimensions dims;
    dims.shape = GeometryShape::PLATE;
    dims.length = 0.01;
    dims.width = 0.01;
    dims.thickness = 0.02;

    MeshGenerationOptions options;
    options.nx = 1;
    options.ny = 1;
    options.nz = 6;
    options.tetrahedralize = true;

    auto mesh = GeometryMeshGenerator::generateFromDimensions(dims, options);
    EXPECT_TRUE(mesh);
    EXPECT_TRUE(mesh->validate());

    return std::shared_ptr<VolMesh>(std::move(mesh));
}

} // namespace

TEST(PoissonSolverTest, SolvesFiniteDirichletProfileOnPlateMesh)
{
    auto mesh = makePlateMesh();
    ASSERT_TRUE(mesh);

    const double z_min = mesh->getBoundingBoxMin().z();
    const double z_max = mesh->getBoundingBoxMax().z();

    PoissonSolver solver(mesh);
    for (const auto& node : mesh->getNodes())
    {
        ASSERT_TRUE(node);
        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12)
        {
            node->setBoundaryType(BoundaryType::DIRICHLET);
            solver.addBoundaryCondition(
                node->getId(),
                std::make_shared<BoundaryCondition>(BoundaryConditionType::DIRICHLET, 10.0));
        }
        else if (std::abs(z - z_max) < 1.0e-12)
        {
            node->setBoundaryType(BoundaryType::DIRICHLET);
            solver.addBoundaryCondition(
                node->getId(),
                std::make_shared<BoundaryCondition>(BoundaryConditionType::DIRICHLET, 0.0));
        }
        else
        {
            node->setBoundaryType(BoundaryType::INTERIOR);
        }
    }

    ASSERT_TRUE(solver.solve(SolverType::DIRECT));

    bool found_interior = false;
    double max_field = 0.0;
    for (const auto& node : mesh->getNodes())
    {
        ASSERT_TRUE(node);
        const double phi = node->getPotential();
        EXPECT_TRUE(std::isfinite(phi));

        if (std::abs(node->getPosition().z() - z_min) < 1.0e-12)
        {
            EXPECT_NEAR(phi, 10.0, 1.0e-8);
        }
        else if (std::abs(node->getPosition().z() - z_max) < 1.0e-12)
        {
            EXPECT_NEAR(phi, 0.0, 1.0e-8);
        }
        else
        {
            found_interior = true;
            EXPECT_GE(phi, -1.0e-8);
            EXPECT_LE(phi, 10.0 + 1.0e-8);
        }

        const auto e = solver.getElectricField(node->getId());
        EXPECT_TRUE(std::isfinite(e.x()));
        EXPECT_TRUE(std::isfinite(e.y()));
        EXPECT_TRUE(std::isfinite(e.z()));
        max_field = std::max(max_field, e.magnitude());
    }

    EXPECT_TRUE(found_interior);
    EXPECT_GT(max_field, 0.0);
    EXPECT_GE(solver.getTotalEnergy(), 0.0);
    EXPECT_GE(solver.getMaxPotential(), solver.getMinPotential());
}
