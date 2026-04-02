#include "../include/MultiBoundaryElectricFieldSolver.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>

namespace
{

using SCDAT::FieldSolver::BoundaryConditionType;
using SCDAT::FieldSolver::MultiBoundaryElectricFieldSolver;
using SCDAT::FieldSolver::SolverConfiguration;
using SCDAT::FieldSolver::SolverType;
using SCDAT::Geometry::Point3D;
using SCDAT::Mesh::ElementType;
using SCDAT::Mesh::Mesh;
using SCDAT::Mesh::MeshPtr;

MeshPtr makeSingleTetraMesh()
{
    auto mesh = std::make_shared<Mesh>();
    mesh->addNode(Point3D(0.0, 0.0, 0.0));
    mesh->addNode(Point3D(1.0, 0.0, 0.0));
    mesh->addNode(Point3D(0.0, 1.0, 0.0));
    mesh->addNode(Point3D(0.0, 0.0, 1.0));
    mesh->addElement(ElementType::TETRAHEDRON, {0, 1, 2, 3}, 0);
    return mesh;
}

MeshPtr makeSkewedTetraMesh()
{
    auto mesh = std::make_shared<Mesh>();
    mesh->addNode(Point3D(0.0, 0.0, 0.0));
    mesh->addNode(Point3D(1.25, 0.05, 0.0));
    mesh->addNode(Point3D(0.15, 0.9, 0.2));
    mesh->addNode(Point3D(0.0, 0.25, 1.1));
    mesh->addElement(ElementType::TETRAHEDRON, {0, 1, 2, 3}, 0);
    return mesh;
}

MeshPtr makeUnitCubeTetraMesh()
{
    auto mesh = std::make_shared<Mesh>();
    mesh->addNode(Point3D(0.0, 0.0, 0.0));
    mesh->addNode(Point3D(1.0, 0.0, 0.0));
    mesh->addNode(Point3D(0.0, 1.0, 0.0));
    mesh->addNode(Point3D(1.0, 1.0, 0.0));
    mesh->addNode(Point3D(0.0, 0.0, 1.0));
    mesh->addNode(Point3D(1.0, 0.0, 1.0));
    mesh->addNode(Point3D(0.0, 1.0, 1.0));
    mesh->addNode(Point3D(1.0, 1.0, 1.0));

    mesh->addElement(ElementType::TETRAHEDRON, {0, 1, 2, 4}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {1, 3, 2, 7}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {1, 2, 4, 7}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {1, 4, 5, 7}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {2, 4, 6, 7}, 0);
    return mesh;
}

MeshPtr makeObliquePrismTetraMesh()
{
    auto mesh = std::make_shared<Mesh>();

    const Point3D p0(0.0, 0.0, 0.0);
    const Point3D p1(1.0, 0.1, 0.2);
    const Point3D p2(0.2, 0.9, -0.1);
    const Point3D shift(0.35, 0.4, 1.15);

    mesh->addNode(p0);
    mesh->addNode(p1);
    mesh->addNode(p2);
    mesh->addNode(p0 + shift);
    mesh->addNode(p1 + shift);
    mesh->addNode(p2 + shift);

    mesh->addElement(ElementType::TETRAHEDRON, {0, 1, 2, 3}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {1, 2, 3, 4}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {2, 3, 4, 5}, 0);
    return mesh;
}

MeshPtr makeRotatedPrismTetraMesh()
{
    auto mesh = std::make_shared<Mesh>();

    const Point3D p0(0.0, 0.0, 0.0);
    const Point3D p1(1.2, 0.15, 0.0);
    const Point3D p2(0.25, 0.95, 0.0);
    const Point3D centroid((p0.x() + p1.x() + p2.x()) / 3.0, (p0.y() + p1.y() + p2.y()) / 3.0, 0.0);
    const double angle = 0.9;
    const double cos_theta = std::cos(angle);
    const double sin_theta = std::sin(angle);
    const auto rotateAroundCentroid = [&](const Point3D& point)
    {
        const double dx = point.x() - centroid.x();
        const double dy = point.y() - centroid.y();
        return Point3D(centroid.x() + cos_theta * dx - sin_theta * dy,
                       centroid.y() + sin_theta * dx + cos_theta * dy, 1.1);
    };

    mesh->addNode(p0);
    mesh->addNode(p1);
    mesh->addNode(p2);
    mesh->addNode(rotateAroundCentroid(p0));
    mesh->addNode(rotateAroundCentroid(p1));
    mesh->addNode(rotateAroundCentroid(p2));

    mesh->addElement(ElementType::TETRAHEDRON, {0, 1, 2, 3}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {1, 2, 3, 4}, 0);
    mesh->addElement(ElementType::TETRAHEDRON, {2, 3, 4, 5}, 0);
    return mesh;
}

void addReferenceDirichletBoundary(MultiBoundaryElectricFieldSolver& solver)
{
    solver.addBoundaryCondition("node0", BoundaryConditionType::DIRICHLET, 5.0, {0});
    solver.addBoundaryCondition("node1", BoundaryConditionType::DIRICHLET, 0.0, {1});
    solver.addBoundaryCondition("node2", BoundaryConditionType::DIRICHLET, 0.0, {2});
    solver.addBoundaryCondition("node3", BoundaryConditionType::DIRICHLET, 0.0, {3});
}

} // namespace

TEST(MultiBoundaryElectricFieldSolverTest, DirectModeRespectsDirichletValues)
{
    auto mesh = makeSingleTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);
    addReferenceDirichletBoundary(solver);

    ASSERT_TRUE(solver.solve());

    EXPECT_NEAR(solver.getPotential(0), 5.0, 1.0e-8);
    EXPECT_NEAR(solver.getPotential(1), 0.0, 1.0e-8);
    EXPECT_NEAR(solver.getPotential(2), 0.0, 1.0e-8);
    EXPECT_NEAR(solver.getPotential(3), 0.0, 1.0e-8);

    const auto field = solver.getElectricField(0);
    EXPECT_TRUE(std::isfinite(field.x()));
    EXPECT_TRUE(std::isfinite(field.y()));
    EXPECT_TRUE(std::isfinite(field.z()));
    EXPECT_GT(field.magnitude(), 0.0);
    EXPECT_GE(solver.calculateElectrostaticEnergy(), 0.0);
}

TEST(MultiBoundaryElectricFieldSolverTest, MultigridModeProducesFiniteResidual)
{
    auto mesh = makeSingleTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::MULTIGRID;
    config.max_iterations = 20;
    config.tolerance = 1.0e-10;
    solver.setSolverConfiguration(config);
    addReferenceDirichletBoundary(solver);

    ASSERT_TRUE(solver.solve());
    EXPECT_FALSE(solver.getResidualHistory().empty());
    EXPECT_TRUE(std::isfinite(solver.validateSolution()));

    const auto at_center = solver.getElectricField(Point3D(0.1, 0.1, 0.1));
    EXPECT_TRUE(std::isfinite(at_center.x()));
    EXPECT_TRUE(std::isfinite(at_center.y()));
    EXPECT_TRUE(std::isfinite(at_center.z()));
}

TEST(MultiBoundaryElectricFieldSolverTest, PeriodicPairMaintainsConfiguredPotentialOffset)
{
    auto mesh = makeSkewedTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);
    solver.addBoundaryCondition("node0", BoundaryConditionType::DIRICHLET, 5.0, {0});
    solver.addBoundaryCondition("node3", BoundaryConditionType::DIRICHLET, 0.0, {3});
    solver.addPeriodicBoundaryCondition("periodic_pair", {1}, {2}, 2.0);

    ASSERT_TRUE(solver.solve());

    EXPECT_NEAR(solver.getPotential(2) - solver.getPotential(1), 2.0, 1.0e-6);
    EXPECT_TRUE(std::isfinite(solver.validateSolution()));
}

TEST(MultiBoundaryElectricFieldSolverTest, FloatingBoundaryGroupBecomesEquipotential)
{
    auto mesh = makeSkewedTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);
    solver.addBoundaryCondition("node0", BoundaryConditionType::DIRICHLET, 5.0, {0});
    solver.addBoundaryCondition("node3", BoundaryConditionType::DIRICHLET, 0.0, {3});
    solver.addBoundaryCondition("floating_cap", BoundaryConditionType::FLOATING, 1.0, {1, 2});

    ASSERT_TRUE(solver.solve());

    EXPECT_NEAR(solver.getPotential(1), solver.getPotential(2), 1.0e-6);
    const auto field = solver.getElectricField(Point3D(0.2, 0.2, 0.2));
    EXPECT_TRUE(std::isfinite(field.x()));
    EXPECT_TRUE(std::isfinite(field.y()));
    EXPECT_TRUE(std::isfinite(field.z()));
}

TEST(MultiBoundaryElectricFieldSolverTest, AutoPeriodicPairInferenceMatchesOppositeFaceNodes)
{
    auto mesh = makeUnitCubeTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);
    solver.addPeriodicBoundaryCondition("auto_x_periodic", {0, 2, 4, 6}, {}, 1.5);
    solver.addBoundaryCondition("x0_y0_z0", BoundaryConditionType::DIRICHLET, 5.0, {0});
    solver.addBoundaryCondition("x0_y1_z0", BoundaryConditionType::DIRICHLET, 3.0, {2});
    solver.addBoundaryCondition("x0_y0_z1", BoundaryConditionType::DIRICHLET, 1.0, {4});
    solver.addBoundaryCondition("x0_y1_z1", BoundaryConditionType::DIRICHLET, -1.0, {6});

    ASSERT_TRUE(solver.solve());

    EXPECT_NEAR(solver.getPotential(1) - solver.getPotential(0), 1.5, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(3) - solver.getPotential(2), 1.5, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(5) - solver.getPotential(4), 1.5, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(7) - solver.getPotential(6), 1.5, 1.0e-6);
}

TEST(MultiBoundaryElectricFieldSolverTest, AutoPeriodicPairInferenceSupportsObliqueTranslatedFaces)
{
    auto mesh = makeObliquePrismTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);
    solver.addPeriodicBoundaryCondition("auto_oblique_periodic", {0, 1, 2}, {}, 0.75);
    solver.addBoundaryCondition("base0", BoundaryConditionType::DIRICHLET, 5.0, {0});
    solver.addBoundaryCondition("base1", BoundaryConditionType::DIRICHLET, 2.0, {1});
    solver.addBoundaryCondition("base2", BoundaryConditionType::DIRICHLET, -1.0, {2});

    ASSERT_TRUE(solver.solve());

    EXPECT_NEAR(solver.getPotential(3) - solver.getPotential(0), 0.75, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(4) - solver.getPotential(1), 0.75, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(5) - solver.getPotential(2), 0.75, 1.0e-6);
}

TEST(MultiBoundaryElectricFieldSolverTest, AutoPeriodicPairInferenceSupportsRotatedCongruentFaces)
{
    auto mesh = makeRotatedPrismTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);
    solver.addPeriodicBoundaryCondition("auto_rotated_periodic", {0, 1, 2}, {}, 0.5);
    solver.addBoundaryCondition("base0", BoundaryConditionType::DIRICHLET, 4.0, {0});
    solver.addBoundaryCondition("base1", BoundaryConditionType::DIRICHLET, 1.5, {1});
    solver.addBoundaryCondition("base2", BoundaryConditionType::DIRICHLET, -0.5, {2});

    ASSERT_TRUE(solver.solve());

    EXPECT_NEAR(solver.getPotential(3) - solver.getPotential(0), 0.5, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(4) - solver.getPotential(1), 0.5, 1.0e-6);
    EXPECT_NEAR(solver.getPotential(5) - solver.getPotential(2), 0.5, 1.0e-6);
}

TEST(MultiBoundaryElectricFieldSolverTest, CoupledMultiBodyStressRemainsStable)
{
    auto mesh = makeUnitCubeTetraMesh();
    ASSERT_TRUE(mesh);

    MultiBoundaryElectricFieldSolver solver(mesh);
    SolverConfiguration config;
    config.solver_type = SolverType::DIRECT_LU;
    solver.setSolverConfiguration(config);

    SCDAT::FieldSolver::CouplingIterationConfiguration coupling;
    coupling.max_outer_iterations = 12;
    coupling.potential_tolerance_v = 5.0e-3;
    coupling.relaxation = 0.65;
    solver.setCouplingIterationConfiguration(coupling);

    solver.addBoundaryCondition("coupled_drive", BoundaryConditionType::DIRICHLET, 3.0, {0, 2});
    solver.addBoundaryCondition("coupled_sink", BoundaryConditionType::DIRICHLET, 0.0, {4, 6});
    solver.addPeriodicBoundaryCondition("coupled_periodic", {1, 3}, {5, 7}, 0.25);

    auto coupling_callback = [&solver](int, const std::vector<double>& previous_potentials) {
        double floating_average = 0.0;
        int count = 0;
        for (const std::size_t node_id : {std::size_t{1}, std::size_t{3}})
        {
            if (node_id < previous_potentials.size())
            {
                floating_average += previous_potentials[node_id];
                ++count;
            }
        }
        if (count > 0)
        {
            floating_average /= static_cast<double>(count);
        }

        const double drive_potential = std::clamp(2.8 + 0.08 * floating_average, 1.0, 5.0);
        solver.addBoundaryCondition("coupled_drive", BoundaryConditionType::DIRICHLET,
                                    drive_potential, {0, 2});
        solver.addBoundaryCondition("coupled_sink", BoundaryConditionType::DIRICHLET,
                                    drive_potential - 2.0, {4, 6});
    };

    ASSERT_TRUE(solver.solveCoupled(coupling_callback));
    EXPECT_GE(solver.getCouplingIterationCount(), 1);
    ASSERT_FALSE(solver.getCouplingResidualHistory().empty());
    EXPECT_LE(solver.getCouplingResidualHistory().back(), coupling.potential_tolerance_v * 1.5);

    EXPECT_NEAR(solver.getPotential(5) - solver.getPotential(1), 0.25, 1.0e-5);
    EXPECT_NEAR(solver.getPotential(7) - solver.getPotential(3), 0.25, 1.0e-5);
    EXPECT_GT(solver.getPotential(0), solver.getPotential(4));
    EXPECT_TRUE(std::isfinite(solver.validateSolution()));

    const auto field = solver.getElectricField(Point3D(0.5, 0.5, 0.5));
    EXPECT_TRUE(std::isfinite(field.x()));
    EXPECT_TRUE(std::isfinite(field.y()));
    EXPECT_TRUE(std::isfinite(field.z()));
}
