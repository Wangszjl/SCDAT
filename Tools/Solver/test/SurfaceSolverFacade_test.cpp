#include "../include/SurfaceSolverFacade.h"

#include <gtest/gtest.h>

#include <cmath>

namespace SCDAT
{
namespace Solver
{

TEST(SurfaceSolverFacadeTest, ResolveGlobalCoupledControlAppliesImplicitOverrides)
{
    GlobalCoupledControlInput input;
    input.configured_iteration_limit = 1;
    input.configured_tolerance_v = 1.0e-4;
    input.configured_relaxation = 0.9;
    input.solver_max_iterations = 20;
    input.solver_residual_tolerance = 1.0e-6;
    input.solver_relaxation_factor = 0.35;
    input.coupling_mode = "field_particle_implicit";
    input.convergence_policy = "residual_norm_guarded";

    const auto control = resolveGlobalCoupledControl(input);

    EXPECT_EQ(control.iteration_limit, static_cast<std::size_t>(20));
    EXPECT_NEAR(control.tolerance_v, 1.0e-6, 1.0e-12);
    EXPECT_NEAR(control.relaxation, 0.35, 1.0e-12);
    EXPECT_TRUE(control.implicit_coupling_requested);
    EXPECT_TRUE(control.residual_guarded_requested);
}

TEST(SurfaceSolverFacadeTest, ResolveGlobalCoupledControlKeepsConfiguredValuesForExplicitMode)
{
    GlobalCoupledControlInput input;
    input.configured_iteration_limit = 1;
    input.configured_tolerance_v = 2.0e-7;
    input.configured_relaxation = 1.2;
    input.solver_max_iterations = 64;
    input.solver_residual_tolerance = 1.0e-9;
    input.solver_relaxation_factor = 0.2;
    input.coupling_mode = "explicit";
    input.convergence_policy = "fixed";

    const auto control = resolveGlobalCoupledControl(input);

    EXPECT_EQ(control.iteration_limit, static_cast<std::size_t>(1));
    EXPECT_NEAR(control.tolerance_v, 2.0e-7, 1.0e-12);
    EXPECT_NEAR(control.relaxation, 1.0, 1.0e-12);
    EXPECT_FALSE(control.implicit_coupling_requested);
    EXPECT_FALSE(control.residual_guarded_requested);
}

TEST(SurfaceSolverFacadeTest, ResolveSolverPolicyFlagsNormalizesTokenVariants)
{
    const auto implicit_residual_flags =
        resolveSolverPolicyFlags("field_particle_implicit", "residual_norm_guarded");
    EXPECT_EQ(implicit_residual_flags.normalized_coupling_mode, "fieldparticleimplicit");
    EXPECT_EQ(implicit_residual_flags.normalized_convergence_policy, "residualnormguarded");
    EXPECT_TRUE(implicit_residual_flags.implicit_coupling_requested);
    EXPECT_TRUE(implicit_residual_flags.residual_guarded_requested);
    EXPECT_FALSE(implicit_residual_flags.fixed_iteration_policy_requested);

    const auto fixed_iteration_flags =
        resolveSolverPolicyFlags("explicit", "fixed_iteration");
    EXPECT_FALSE(fixed_iteration_flags.implicit_coupling_requested);
    EXPECT_FALSE(fixed_iteration_flags.residual_guarded_requested);
    EXPECT_TRUE(fixed_iteration_flags.fixed_iteration_policy_requested);

    EXPECT_TRUE(isResidualGuardedSolverPolicyToken("residual_norm_guarded"));
    EXPECT_TRUE(isFixedIterationSolverPolicyToken("fixed_iteration"));
    EXPECT_FALSE(isImplicitSolverCouplingModeToken("explicit"));
}

TEST(SurfaceSolverFacadeTest, SolveDenseLinearSystemWithResidualSolvesWellConditionedSystem)
{
    std::vector<std::vector<double>> matrix{{4.0, -1.0}, {-1.0, 3.0}};
    std::vector<double> rhs{15.0, 10.0};

    const auto solved = solveDenseLinearSystemWithResidual(matrix, rhs);

    ASSERT_TRUE(solved.solved);
    ASSERT_EQ(solved.solution.size(), static_cast<std::size_t>(2));
    EXPECT_NEAR(solved.solution[0], 5.0, 1.0e-12);
    EXPECT_NEAR(solved.solution[1], 5.0, 1.0e-12);
    EXPECT_LE(std::abs(solved.residual_norm), 1.0e-12);
    EXPECT_EQ(solved.nonzeros, static_cast<std::size_t>(4));
}

TEST(SurfaceSolverFacadeTest, SolveDenseLinearSystemWithResidualRejectsSingularSystem)
{
    std::vector<std::vector<double>> matrix{{1.0, 2.0}, {2.0, 4.0}};
    std::vector<double> rhs{3.0, 6.0};

    const auto solved = solveDenseLinearSystemWithResidual(matrix, rhs);

    EXPECT_FALSE(solved.solved);
    EXPECT_EQ(solved.solution.size(), static_cast<std::size_t>(2));
}

TEST(SurfaceSolverFacadeTest, SolveIterativeLinearSystemConvergesForWellConditionedSystem)
{
    const std::vector<std::vector<double>> matrix{{4.0, -1.0}, {-1.0, 3.0}};
    const std::vector<double> rhs{15.0, 10.0};
    const std::vector<double> initial_guess{0.0, 0.0};

    const auto solved =
        solveIterativeLinearSystem(matrix, rhs, initial_guess, 128, 1.0e-10, 1.0);

    ASSERT_EQ(solved.solution.size(), static_cast<std::size_t>(2));
    EXPECT_TRUE(solved.converged);
    EXPECT_GT(solved.iterations, 0);
    EXPECT_NEAR(solved.solution[0], 5.0, 1.0e-8);
    EXPECT_NEAR(solved.solution[1], 5.0, 1.0e-8);
    EXPECT_LE(std::abs(solved.residual_norm), 1.0e-6);
    EXPECT_EQ(solved.nonzeros, static_cast<std::size_t>(4));
}

TEST(SurfaceSolverFacadeTest, ResolveVolumeLinearSolverRoutingHonorsPolicyHints)
{
    VolumeLinearSolverRoutingInput auto_input;
    auto_input.system_size = 8;
    const auto auto_routing = resolveVolumeLinearSolverRouting(auto_input);
    EXPECT_TRUE(auto_routing.auto_prefers_iterative);
    EXPECT_TRUE(auto_routing.use_iterative_solver);
    EXPECT_EQ(auto_routing.solver_mode, "iterative");

    VolumeLinearSolverRoutingInput dense_only_input;
    dense_only_input.system_size = 32;
    dense_only_input.dense_only = true;
    const auto dense_only_routing = resolveVolumeLinearSolverRouting(dense_only_input);
    EXPECT_FALSE(dense_only_routing.use_iterative_solver);
    EXPECT_EQ(dense_only_routing.solver_mode, "dense");

    VolumeLinearSolverRoutingInput iterative_only_input;
    iterative_only_input.system_size = 1;
    iterative_only_input.iterative_only = true;
    const auto iterative_only_routing =
        resolveVolumeLinearSolverRouting(iterative_only_input);
    EXPECT_TRUE(iterative_only_routing.use_iterative_solver);
    EXPECT_EQ(iterative_only_routing.solver_mode, "iterative");
}

} // namespace Solver
} // namespace SCDAT
