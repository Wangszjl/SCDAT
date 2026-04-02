#include "BoltzmannSolver.h"
#include "DiffusionDriftSolver.h"
#include "NonlinearPoissonSolver.h"

#include <gtest/gtest.h>

using SCDAT::FieldSolver::BoltzmannSolver;
using SCDAT::FieldSolver::DiffusionDriftParameters;
using SCDAT::FieldSolver::DiffusionDriftSolver;
using SCDAT::FieldSolver::NonlinearPoissonSolver;
using SCDAT::Geometry::Point3D;
using SCDAT::Geometry::Vector3D;

TEST(ExtendedFieldSolverTest, BoltzmannSolverBuildsNormalizedDistribution)
{
    BoltzmannSolver solver;
    const auto& state = solver.solve(4.0);
    double sum = 0.0;
    for (double value : state.distribution)
    {
        sum += value;
    }

    EXPECT_NEAR(sum, 1.0, 1.0e-9);
    EXPECT_GT(state.mean_energy_ev, 0.0);
    EXPECT_GT(state.ionization_rate, 0.0);
}

TEST(ExtendedFieldSolverTest, DiffusionDriftSolverAdvancesStructuredDensityField)
{
    DiffusionDriftSolver solver;
    DiffusionDriftParameters parameters;
    parameters.source_rate = 1.0e18;
    parameters.floor_density = 1.0;
    ASSERT_TRUE(solver.setParameters(parameters));
    ASSERT_TRUE(solver.initializeGrid(Vector3D(1.0, 1.0, 1.0), Vector3D(4.0, 4.0, 4.0)));
    ASSERT_TRUE(solver.setInitialDensity([](const Point3D&) { return 10.0; }));

    std::vector<Vector3D> electric_field(64, Vector3D(0.0, 0.0, 1.0));
    ASSERT_TRUE(solver.advance(electric_field, 1.0e-9));
    EXPECT_GT(solver.getDensity(Point3D(0.5, 0.5, 0.5)), 10.0);
}

TEST(ExtendedFieldSolverTest, NonlinearPoissonSolverRespectsDirichletBoundaries)
{
    NonlinearPoissonSolver solver;
    ASSERT_TRUE(solver.initializeGrid(Vector3D(1.0, 1.0, 1.0), Vector3D(4.0, 4.0, 8.0)));
    solver.setBoundaryPotential([](const Point3D& point) { return point.z() > 0.99 ? 10.0 : 0.0; });
    solver.setChargeDensityFunction([](const Point3D&, double) { return 0.0; });
    solver.setPermittivityFunction([](const Point3D&, double) { return 1.0; });

    ASSERT_TRUE(solver.solve());
    const double mid_potential = solver.getPotential(Point3D(0.5, 0.5, 0.5));
    EXPECT_GT(mid_potential, 0.0);
    EXPECT_LT(mid_potential, 10.0);
    EXPECT_FALSE(solver.getElectricField().empty());
}

TEST(ExtendedFieldSolverTest, NonlinearPoissonQuasiNewtonConvergesUnderStrongNonlinearity)
{
    NonlinearPoissonSolver solver;
    SCDAT::FieldSolver::NonlinearPoissonParameters parameters;
    parameters.max_iterations = 700;
    parameters.tolerance = 2.0e-4;
    parameters.relaxation = 0.9;
    parameters.adaptive_relaxation_floor = 0.2;
    parameters.derivative_step_v = 5.0e-3;
    parameters.minimum_effective_diagonal = 1.0e-8;
    parameters.solve_policy = SCDAT::FieldSolver::NonlinearSolvePolicy::QuasiNewton;

    solver.setParameters(parameters);
    ASSERT_TRUE(solver.initializeGrid(Vector3D(1.0, 1.0, 1.0), Vector3D(7.0, 7.0, 15.0)));
    solver.setBoundaryPotential([](const Point3D& point) { return point.z() > 0.99 ? 25.0 : 0.0; });
    solver.setChargeDensityFunction([](const Point3D& point, double phi) {
        const double bounded_phi = std::clamp(phi / 8.0, -3.0, 3.0);
        return 1.0e-6 * std::tanh(bounded_phi) + 5.0e-7 * (point.z() - 0.5);
    });
    solver.setPermittivityFunction([](const Point3D& point, double) {
        return 2.8 + 0.4 * point.z();
    });

    ASSERT_TRUE(solver.solve());
    EXPECT_TRUE(solver.hasConverged());
    EXPECT_GT(solver.getIterationsUsed(), 0);
    EXPECT_LE(solver.getLastMaxDelta(), parameters.tolerance * 1.5);
    EXPECT_TRUE(std::isfinite(solver.getLastResidualNorm()));

    const double mid_potential = solver.getPotential(Point3D(0.5, 0.5, 0.5));
    EXPECT_GT(mid_potential, 0.0);
    EXPECT_LT(mid_potential, 25.0);
    EXPECT_FALSE(solver.getElectricField().empty());
}
