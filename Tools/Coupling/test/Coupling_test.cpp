#include "CouplingBase.h"
#include "MultiPhysicsManager.h"
#include "SurfaceCircuitCoupling.h"

#include <gtest/gtest.h>

using SCDAT::Coupling::ConvergenceCriteria;
using SCDAT::Coupling::CouplingType;
using SCDAT::Coupling::FunctionalCoupling;
using SCDAT::Coupling::MultiPhysicsManager;
using SCDAT::Coupling::SurfaceCircuitBranch;
using SCDAT::Coupling::SurfaceCircuitCoupling;
using SCDAT::Coupling::SurfaceCircuitLinearization;
using SCDAT::Coupling::SurfaceCircuitNode;

TEST(CouplingTest, IterativeManagerDrivesFunctionalCouplingToConvergence)
{
    double residual = 1.0;
    auto coupling = std::make_shared<FunctionalCoupling>(
        "test", CouplingType::ITERATIVE,
        [&residual](double) {
            residual *= 0.25;
            return SCDAT::VoidResult::success();
        },
        [&residual] { return std::pair<double, double>{residual, residual}; });

    ConvergenceCriteria criteria;
    criteria.relative_tolerance = 0.01;
    criteria.absolute_tolerance = 0.01;
    criteria.max_iterations = 20;
    criteria.min_iterations = 1;
    coupling->setConvergenceCriteria(criteria);

    MultiPhysicsManager manager;
    manager.addCoupling(coupling);

    ASSERT_TRUE(manager.executeIterativeCouplings(1.0));
    const auto stats = manager.getAllStatistics();
    ASSERT_TRUE(stats.contains("test"));
    EXPECT_TRUE(stats.at("test").converged);
    EXPECT_LT(stats.at("test").final_absolute_error, 0.01);
}

TEST(CouplingTest, SurfaceCircuitCouplingAdvancesTwoNodePatchTowardInjectedPotential)
{
    SurfaceCircuitCoupling circuit;
    SurfaceCircuitNode body;
    body.name = "body";
    body.potential_v = 0.0;
    body.fixed_potential = true;
    body.fixed_value_v = 0.0;
    const auto body_index = circuit.addNode(body);

    SurfaceCircuitNode patch;
    patch.name = "patch";
    patch.potential_v = 0.0;
    patch.capacitance_f = 1.0e-9;
    const auto patch_index = circuit.addNode(patch);

    SurfaceCircuitBranch branch;
    branch.from_node = patch_index;
    branch.to_node = body_index;
    branch.conductance_s = 1.0e-9;
    circuit.addBranch(branch);

    SurfaceCircuitLinearization linearization;
    linearization.plasma_current_a = {0.0, 1.0e-9};
    linearization.plasma_didv_a_per_v = {0.0, -1.0e-11};

    const auto result = circuit.advanceImplicit(1.0, linearization, 50.0);
    ASSERT_TRUE(result.converged);
    ASSERT_EQ(result.node_potentials_v.size(), 2u);
    EXPECT_GT(result.node_potentials_v[patch_index], 0.0);
    EXPECT_LT(result.branch_currents_a.front(), result.plasma_currents_a[patch_index]);
}

TEST(CouplingTest, SurfaceCircuitCouplingRespectsPotentialStepLimit)
{
    SurfaceCircuitCoupling circuit;
    SurfaceCircuitNode body;
    body.name = "body";
    body.fixed_potential = true;
    body.fixed_value_v = 0.0;
    const auto body_index = circuit.addNode(body);

    SurfaceCircuitNode patch;
    patch.name = "patch";
    patch.capacitance_f = 1.0e-12;
    const auto patch_index = circuit.addNode(patch);

    SurfaceCircuitBranch branch;
    branch.from_node = patch_index;
    branch.to_node = body_index;
    branch.conductance_s = 0.0;
    circuit.addBranch(branch);

    SurfaceCircuitLinearization linearization;
    linearization.plasma_current_a = {0.0, 1.0e-6};
    linearization.plasma_didv_a_per_v = {0.0, 0.0};

    const auto result = circuit.advanceImplicit(1.0, linearization, 5.0);
    ASSERT_TRUE(result.converged);
    EXPECT_NEAR(result.node_potentials_v[patch_index], 5.0, 1.0e-9);
}

TEST(CouplingTest, SurfaceCircuitCouplingAppliesAdditionalMatrixEntries)
{
    SurfaceCircuitCoupling circuit;
    SurfaceCircuitNode left;
    left.name = "left";
    left.potential_v = 1.0;
    left.capacitance_f = 1.0;
    const auto left_index = circuit.addNode(left);

    SurfaceCircuitNode right;
    right.name = "right";
    right.potential_v = 0.0;
    right.capacitance_f = 1.0;
    const auto right_index = circuit.addNode(right);

    SurfaceCircuitLinearization linearization;
    linearization.plasma_current_a = {0.0, 0.0};
    linearization.plasma_didv_a_per_v = {0.0, 0.0};
    linearization.additional_rhs_a = {0.0, 0.0};
    linearization.additional_diagonal_a_per_v = {1.0, 1.0};
    linearization.additional_off_diagonal_entries.push_back(
        {left_index, right_index, -1.0});
    linearization.additional_off_diagonal_entries.push_back(
        {right_index, left_index, -1.0});

    const auto result = circuit.advanceImplicit(1.0, linearization);
    ASSERT_TRUE(result.converged);
    ASSERT_EQ(result.node_potentials_v.size(), 2u);
    EXPECT_NEAR(result.node_potentials_v[left_index], 2.0 / 3.0, 1.0e-9);
    EXPECT_NEAR(result.node_potentials_v[right_index], 1.0 / 3.0, 1.0e-9);
}
