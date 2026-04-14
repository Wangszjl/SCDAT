#include "CouplingBase.h"
#include "MultiPhysicsManager.h"
#include "SurfaceCurrentLinearizationKernel.h"
#include "SurfaceCircuitCoupling.h"

#include <gtest/gtest.h>
#include <cmath>
#include <limits>

using SCDAT::Coupling::ConvergenceCriteria;
using SCDAT::Coupling::CircuitAssembly;
using SCDAT::Coupling::CircuitBranchDescriptor;
using SCDAT::Coupling::CircuitDidvCompositionInput;
using SCDAT::Coupling::CircuitDidvAggregationKind;
using SCDAT::Coupling::CircuitDidvSurfaceComposition;
using SCDAT::Coupling::CircuitKernelTopologyState;
using SCDAT::Coupling::CircuitDidvTerm;
using SCDAT::Coupling::CircuitWeightedDidvTerm;
using SCDAT::Coupling::CircuitExcitationDescriptor;
using SCDAT::Coupling::CircuitExcitationWaveformKind;
using SCDAT::Coupling::CircuitElementDescriptor;
using SCDAT::Coupling::CircuitElementKind;
using SCDAT::Coupling::CircuitNodeDescriptor;
using SCDAT::Coupling::CouplingType;
using SCDAT::Coupling::FunctionalCoupling;
using SCDAT::Coupling::MultiPhysicsManager;
using SCDAT::Coupling::SurfaceCircuitBranch;
using SCDAT::Coupling::SurfaceCircuitCoupling;
using SCDAT::Coupling::SurfaceCircuitKernelBranchState;
using SCDAT::Coupling::SurfaceCircuitKernelInput;
using SCDAT::Coupling::SurfaceCircuitKernelNodeState;
using SCDAT::Coupling::SurfaceCircuitLinearization;
using SCDAT::Coupling::SurfaceCurrentLinearizationKernel;
using SCDAT::Coupling::SurfaceCurrentLinearizationNode;
using SCDAT::Coupling::SurfaceCircuitNode;
using SCDAT::Coupling::SurfaceAdaptiveSubstepInputs;
using SCDAT::Coupling::SurfaceSheathCapacitanceConsistencyInputs;

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

TEST(CouplingTest, GenericCircuitAssemblyAppliesNodeBranchAndSourceElements)
{
    CircuitAssembly assembly;
    assembly.nodes = {
        CircuitNodeDescriptor{"body", 0.0},
        CircuitNodeDescriptor{"patch", 1.0},
    };
    assembly.branches = {
        CircuitBranchDescriptor{"patch_to_body", 1, 0},
    };
    assembly.elements = {
        CircuitElementDescriptor{CircuitElementKind::NodeFixedPotential, 0, 0.0},
        CircuitElementDescriptor{CircuitElementKind::NodeCapacitance, 1, 2.0e-9},
        CircuitElementDescriptor{CircuitElementKind::NodeShuntConductance, 1, 3.0e-9},
        CircuitElementDescriptor{CircuitElementKind::BranchConductance, 0, 4.0e-9},
        CircuitElementDescriptor{CircuitElementKind::BranchBias, 0, 5.0},
        CircuitElementDescriptor{CircuitElementKind::BranchCurrentSource, 0, -6.0e-10},
    };

    const auto assembled = SCDAT::Coupling::assembleSurfaceCircuit(assembly);
    ASSERT_EQ(assembled.nodes.size(), 2u);
    ASSERT_EQ(assembled.branches.size(), 1u);
    ASSERT_EQ(assembled.kernel_input.branches.size(), 1u);

    EXPECT_TRUE(assembled.nodes[0].fixed_potential);
    EXPECT_NEAR(assembled.nodes[0].fixed_value_v, 0.0, 1.0e-15);
    EXPECT_NEAR(assembled.nodes[1].capacitance_f, 2.0e-9, 1.0e-21);
    EXPECT_NEAR(assembled.nodes[1].shunt_conductance_s, 3.0e-9, 1.0e-21);
    EXPECT_NEAR(assembled.branches[0].conductance_s, 4.0e-9, 1.0e-21);
    EXPECT_NEAR(assembled.branches[0].bias_v, 5.0, 1.0e-15);
    EXPECT_TRUE(assembled.kernel_input.branches[0].valid);
    EXPECT_NEAR(assembled.kernel_input.branches[0].current_source_a, -6.0e-10, 1.0e-24);
}

TEST(CouplingTest, SurfaceCircuitCouplingConfigureUsesGenericAssembly)
{
    CircuitAssembly assembly;
    assembly.nodes = {
        CircuitNodeDescriptor{"body", 0.0},
        CircuitNodeDescriptor{"patch", 0.0},
    };
    assembly.branches = {
        CircuitBranchDescriptor{"patch_to_body", 1, 0},
    };
    assembly.elements = {
        CircuitElementDescriptor{CircuitElementKind::NodeFixedPotential, 0, 0.0},
        CircuitElementDescriptor{CircuitElementKind::NodeCapacitance, 1, 1.0e-9},
        CircuitElementDescriptor{CircuitElementKind::BranchConductance, 0, 2.0e-9},
        CircuitElementDescriptor{CircuitElementKind::BranchCurrentSource, 0, -5.0e-10},
    };

    SurfaceCircuitCoupling circuit;
    circuit.configure(assembly);
    const auto assembled = SCDAT::Coupling::assembleSurfaceCircuit(assembly);
    auto kernel_input = assembled.kernel_input;
    kernel_input.nodes.resize(2);
    kernel_input.nodes[1].plasma_current_a = 2.0e-9;
    kernel_input.nodes[1].plasma_didv_a_per_v = -1.0e-10;
    kernel_input.nodes[1].valid = true;

    const auto result = circuit.advanceImplicit(1.0, kernel_input, 50.0);
    ASSERT_TRUE(result.converged);
    ASSERT_EQ(result.node_potentials_v.size(), 2u);
    EXPECT_GT(result.node_potentials_v[1], 0.0);
    EXPECT_TRUE(std::isfinite(result.branch_currents_a[0]));
    EXPECT_GT(std::abs(result.branch_currents_a[0]), 1.0e-12);
}

TEST(CouplingTest, CircuitDidvCompositionAggregatesTermsAndBuildsKernelInput)
{
    CircuitDidvCompositionInput input;
    input.nodes.resize(2);
    input.nodes[0].terms = {
        CircuitDidvTerm{"collection", 2.0e-9, 1.0e-11},
        CircuitDidvTerm{"emission", -5.0e-10, -2.0e-12},
    };
    input.nodes[0].additional_rhs_a = 3.0e-10;
    input.nodes[0].additional_diagonal_a_per_v = 4.0e-12;
    input.nodes[1].terms = {
        CircuitDidvTerm{"collection", -1.0e-9, -8.0e-12},
    };
    input.off_diagonal_entries.push_back(
        SCDAT::Coupling::SurfaceCircuitLinearization::OffDiagonalEntry{0, 1, -1.5e-12});

    const auto result = SCDAT::Coupling::composeCircuitDidv(input);
    ASSERT_EQ(result.node_summaries.size(), 2u);
    ASSERT_EQ(result.kernel_input.nodes.size(), 2u);
    ASSERT_EQ(result.kernel_input.off_diagonal_entries.size(), 1u);

    EXPECT_NEAR(result.node_summaries[0].total_current_a, 1.5e-9, 1.0e-24);
    EXPECT_NEAR(result.node_summaries[0].total_didv_a_per_v, 8.0e-12, 1.0e-27);
    EXPECT_NEAR(result.linearization.plasma_current_a[0], 1.5e-9, 1.0e-24);
    EXPECT_NEAR(result.linearization.plasma_didv_a_per_v[0], 8.0e-12, 1.0e-27);
    EXPECT_NEAR(result.kernel_input.nodes[0].plasma_current_a, 1.5e-9, 1.0e-24);
    EXPECT_NEAR(result.kernel_input.nodes[0].plasma_didv_a_per_v, 8.0e-12, 1.0e-27);
    EXPECT_NEAR(result.kernel_input.nodes[0].additional_rhs_a, 3.0e-10, 1.0e-24);
    EXPECT_NEAR(result.kernel_input.nodes[0].additional_diagonal_a_per_v, 4.0e-12, 1.0e-27);
    EXPECT_TRUE(result.kernel_input.nodes[0].valid);
}

TEST(CouplingTest, CircuitKernelInputComposerMergesDidvAndTopologyState)
{
    CircuitDidvCompositionInput didv_input;
    didv_input.nodes.resize(1);
    didv_input.nodes[0].terms = {
        CircuitDidvTerm{"plasma", 2.0e-9, -3.0e-11},
    };
    didv_input.nodes[0].additional_rhs_a = 1.0e-10;
    didv_input.nodes[0].additional_diagonal_a_per_v = 2.0e-12;

    CircuitKernelTopologyState topology;
    topology.node_capacitance_f = {5.0e-9};
    topology.node_shunt_conductance_s = {6.0e-9};
    topology.branch_conductance_s = {7.0e-9};
    topology.branch_bias_v = {4.0};
    topology.branch_current_source_a = {-8.0e-10};

    const auto kernel_input = SCDAT::Coupling::composeCircuitKernelInput(didv_input, topology);
    ASSERT_EQ(kernel_input.nodes.size(), 1u);
    ASSERT_EQ(kernel_input.branches.size(), 1u);
    EXPECT_NEAR(kernel_input.nodes[0].plasma_current_a, 2.0e-9, 1.0e-24);
    EXPECT_NEAR(kernel_input.nodes[0].plasma_didv_a_per_v, -3.0e-11, 1.0e-27);
    EXPECT_NEAR(kernel_input.nodes[0].additional_rhs_a, 1.0e-10, 1.0e-24);
    EXPECT_NEAR(kernel_input.nodes[0].additional_diagonal_a_per_v, 2.0e-12, 1.0e-27);
    EXPECT_NEAR(kernel_input.nodes[0].capacitance_f, 5.0e-9, 1.0e-24);
    EXPECT_NEAR(kernel_input.nodes[0].shunt_conductance_s, 6.0e-9, 1.0e-24);
    EXPECT_NEAR(kernel_input.branches[0].conductance_s, 7.0e-9, 1.0e-24);
    EXPECT_NEAR(kernel_input.branches[0].bias_v, 4.0, 1.0e-15);
    EXPECT_NEAR(kernel_input.branches[0].current_source_a, -8.0e-10, 1.0e-24);
    EXPECT_TRUE(kernel_input.nodes[0].valid);
    EXPECT_TRUE(kernel_input.branches[0].valid);
}

TEST(CouplingTest, CircuitDidvCompositionSupportsWeightedTerms)
{
    CircuitDidvCompositionInput input;
    input.nodes.resize(1);
    auto& node = input.nodes[0];
    node.aggregation_kind = CircuitDidvAggregationKind::WeightedTerms;
    node.weighted_terms = {
        CircuitWeightedDidvTerm{"surface_a", 0.25, 4.0e-9, 2.0e-11},
        CircuitWeightedDidvTerm{"surface_b", 0.75, 8.0e-9, 1.0e-11},
    };

    const auto result = SCDAT::Coupling::composeCircuitDidv(input);
    ASSERT_EQ(result.node_summaries.size(), 1u);
    EXPECT_NEAR(result.node_summaries[0].total_current_a, 7.0e-9, 1.0e-24);
    EXPECT_NEAR(result.node_summaries[0].total_didv_a_per_v, 1.25e-11, 1.0e-26);
    EXPECT_NEAR(result.kernel_input.nodes[0].plasma_current_a, 7.0e-9, 1.0e-24);
    EXPECT_NEAR(result.kernel_input.nodes[0].plasma_didv_a_per_v, 1.25e-11, 1.0e-26);
}

TEST(CouplingTest, CircuitDidvCompositionSupportsMultiSurfaceAggregation)
{
    CircuitDidvCompositionInput input;
    input.nodes.resize(1);
    auto& node = input.nodes[0];
    node.aggregation_kind = CircuitDidvAggregationKind::MultiSurface;

    CircuitDidvSurfaceComposition surface_a;
    surface_a.name = "surface_a";
    surface_a.surface_weight = 1.0;
    surface_a.terms = {
        CircuitDidvTerm{"collection", 2.0e-9, 4.0e-12},
        CircuitDidvTerm{"emission", -1.0e-9, -1.0e-12},
    };

    CircuitDidvSurfaceComposition surface_b;
    surface_b.name = "surface_b";
    surface_b.surface_weight = 0.5;
    surface_b.terms = {
        CircuitDidvTerm{"collection", 4.0e-9, 8.0e-12},
    };

    node.surface_compositions = {surface_a, surface_b};

    const auto result = SCDAT::Coupling::composeCircuitDidv(input);
    ASSERT_EQ(result.node_summaries.size(), 1u);
    EXPECT_NEAR(result.node_summaries[0].total_current_a, 3.0e-9, 1.0e-24);
    EXPECT_NEAR(result.node_summaries[0].total_didv_a_per_v, 7.0e-12, 1.0e-27);
    EXPECT_NEAR(result.linearization.plasma_current_a[0], 3.0e-9, 1.0e-24);
    EXPECT_NEAR(result.linearization.plasma_didv_a_per_v[0], 7.0e-12, 1.0e-27);
}

TEST(CouplingTest, CircuitExcitationSamplingSupportsAllWaveforms)
{
    std::vector<CircuitExcitationDescriptor> excitations;
    CircuitExcitationDescriptor constant_source;
    constant_source.name = "constant_source";
    constant_source.target_kind = CircuitElementKind::BranchCurrentSource;
    constant_source.target_index = 0;
    constant_source.waveform = CircuitExcitationWaveformKind::Constant;
    constant_source.base_value = 1.0e-9;
    excitations.push_back(constant_source);

    CircuitExcitationDescriptor sinusoidal_bias;
    sinusoidal_bias.name = "sinusoidal_bias";
    sinusoidal_bias.target_kind = CircuitElementKind::BranchBias;
    sinusoidal_bias.target_index = 1;
    sinusoidal_bias.waveform = CircuitExcitationWaveformKind::Sinusoidal;
    sinusoidal_bias.base_value = 2.0;
    sinusoidal_bias.amplitude = 1.5;
    sinusoidal_bias.frequency_hz = 1.0;
    excitations.push_back(sinusoidal_bias);

    CircuitExcitationDescriptor pulse_conductance;
    pulse_conductance.name = "pulse_conductance";
    pulse_conductance.target_kind = CircuitElementKind::BranchConductance;
    pulse_conductance.target_index = 2;
    pulse_conductance.waveform = CircuitExcitationWaveformKind::Pulse;
    pulse_conductance.base_value = 2.0e-9;
    pulse_conductance.amplitude = 3.0e-9;
    pulse_conductance.frequency_hz = 2.0;
    pulse_conductance.duty_cycle = 0.25;
    excitations.push_back(pulse_conductance);

    CircuitExcitationDescriptor pwl_bias;
    pwl_bias.name = "pwl_bias";
    pwl_bias.target_kind = CircuitElementKind::BranchBias;
    pwl_bias.target_index = 3;
    pwl_bias.waveform = CircuitExcitationWaveformKind::PiecewiseLinear;
    pwl_bias.base_value = 10.0;
    pwl_bias.pwl_points = {
        {0.0, 0.0},
        {0.5, 4.0},
        {1.0, 2.0},
    };
    excitations.push_back(pwl_bias);

    CircuitExcitationDescriptor exp_source;
    exp_source.name = "exp_source";
    exp_source.target_kind = CircuitElementKind::BranchCurrentSource;
    exp_source.target_index = 4;
    exp_source.waveform = CircuitExcitationWaveformKind::Exponential;
    exp_source.base_value = 1.0;
    exp_source.amplitude = 4.0;
    exp_source.exponential_tau_s = 0.5;
    excitations.push_back(exp_source);

    const auto active_sample = SCDAT::Coupling::sampleCircuitExcitations(excitations, 0.0);
    ASSERT_EQ(active_sample.elements.size(), 5u);
    EXPECT_NEAR(active_sample.elements[0].value, 1.0e-9, 1.0e-24);
    EXPECT_NEAR(active_sample.elements[1].value, 2.0, 1.0e-15);
    EXPECT_NEAR(active_sample.elements[2].value, 5.0e-9, 1.0e-24);
    EXPECT_NEAR(active_sample.elements[3].value, 10.0, 1.0e-15);
    EXPECT_NEAR(active_sample.elements[4].value, 1.0, 1.0e-15);

    const auto mid_sample = SCDAT::Coupling::sampleCircuitExcitations(excitations, 0.25);
    ASSERT_EQ(mid_sample.elements.size(), 5u);
    EXPECT_NEAR(mid_sample.elements[2].value, 2.0e-9, 1.0e-24);
    EXPECT_NEAR(mid_sample.elements[3].value, 12.0, 1.0e-15);
    EXPECT_NEAR(mid_sample.elements[4].value, 1.0 + 4.0 * (1.0 - std::exp(-0.5)), 1.0e-12);

    const auto late_sample = SCDAT::Coupling::sampleCircuitExcitations(excitations, 1.50);
    ASSERT_EQ(late_sample.elements.size(), 5u);
    EXPECT_NEAR(late_sample.elements[3].value, 12.0, 1.0e-15);
    EXPECT_NEAR(late_sample.elements[4].value, 1.0 + 4.0 * (1.0 - std::exp(-3.0)), 1.0e-12);
}

TEST(CouplingTest, SurfaceAssemblyAppliesDynamicExcitationOverrides)
{
    CircuitAssembly assembly;
    assembly.nodes = {
        CircuitNodeDescriptor{"body", 0.0},
        CircuitNodeDescriptor{"patch", 0.0},
    };
    assembly.branches = {
        CircuitBranchDescriptor{"patch_to_body", 1, 0},
    };
    assembly.elements = {
        CircuitElementDescriptor{CircuitElementKind::NodeFixedPotential, 0, 0.0},
        CircuitElementDescriptor{CircuitElementKind::BranchBias, 0, 1.0},
    };

    std::vector<CircuitExcitationDescriptor> excitations;
    excitations.push_back(CircuitExcitationDescriptor{
        "bias_drive",
        CircuitElementKind::BranchBias,
        0,
        CircuitExcitationWaveformKind::Constant,
        3.0,
        0.0,
    });
    excitations.push_back(CircuitExcitationDescriptor{
        "source_drive",
        CircuitElementKind::BranchCurrentSource,
        0,
        CircuitExcitationWaveformKind::Constant,
        -4.0e-10,
        0.0,
    });

    const auto sampled = SCDAT::Coupling::sampleCircuitExcitations(excitations, 0.0);
    const auto assembled =
        SCDAT::Coupling::assembleSurfaceCircuit(assembly, sampled.elements);
    ASSERT_EQ(assembled.branches.size(), 1u);
    ASSERT_EQ(assembled.kernel_input.branches.size(), 1u);
    EXPECT_NEAR(assembled.branches[0].bias_v, 3.0, 1.0e-15);
    EXPECT_TRUE(assembled.kernel_input.branches[0].valid);
    EXPECT_NEAR(assembled.kernel_input.branches[0].current_source_a, -4.0e-10, 1.0e-24);
}

TEST(CouplingTest, SurfaceCircuitKernelInputAppliesBiasCapacitanceAndBranchSource)
{
    SurfaceCircuitCoupling circuit;
    SurfaceCircuitNode body;
    body.name = "body";
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
    branch.conductance_s = 2.0e-9;
    branch.bias_v = 3.0;
    circuit.addBranch(branch);

    SurfaceCircuitKernelInput kernel_input;
    kernel_input.nodes.resize(2);
    kernel_input.branches.resize(1);
    kernel_input.nodes[patch_index] =
        SurfaceCircuitKernelNodeState{5.0e-9, -1.0e-10, 0.0, 0.0, 2.0e-9, 0.0, true};
    kernel_input.branches[0] =
        SurfaceCircuitKernelBranchState{4.0e-9, 6.0, 1.0e-9, true};

    const auto result = circuit.advanceImplicit(1.0, kernel_input, 50.0);
    ASSERT_TRUE(result.converged);
    ASSERT_EQ(result.node_potentials_v.size(), 2u);
    EXPECT_GT(result.node_potentials_v[patch_index], 0.0);
    EXPECT_TRUE(std::isfinite(result.branch_currents_a[0]));
    EXPECT_GT(std::abs(result.branch_currents_a[0]), 1.0e-12);
}

TEST(CouplingTest, SurfaceCurrentLinearizationKernelAggregatesComponentsAndDidv)
{
    SurfaceCurrentLinearizationKernel kernel;
    std::vector<SurfaceCurrentLinearizationNode> nodes(1);
    auto& components = nodes[0].components;
    components.electron_collection_a = 4.0e-9;
    components.ion_collection_a = 1.0e-9;
    components.secondary_electron_emission_a = 2.0e-9;
    components.backscatter_emission_a = 5.0e-10;
    components.photoelectron_emission_a = 2.5e-10;
    components.conduction_a = 1.5e-9;
    components.electron_collection_didv_a_per_v = 1.0e-11;
    components.secondary_electron_emission_didv_a_per_v = 2.0e-12;
    components.conduction_didv_a_per_v = 5.0e-12;
    nodes[0].additional_rhs_a = 7.0e-10;
    nodes[0].additional_diagonal_a_per_v = 3.0e-12;

    const auto result = kernel.assemble(nodes);
    ASSERT_EQ(result.total_current_a.size(), 1u);
    ASSERT_EQ(result.total_didv_a_per_v.size(), 1u);
    EXPECT_NEAR(result.total_current_a[0], 3.75e-9, 1.0e-21);
    EXPECT_NEAR(result.total_didv_a_per_v[0], 1.3e-11, 1.0e-24);
    EXPECT_TRUE(result.kernel_input.nodes[0].valid);
    EXPECT_NEAR(result.kernel_input.nodes[0].additional_rhs_a, 7.0e-10, 1.0e-24);
}

TEST(CouplingTest, CenteredDidvEstimatorMatchesLinearSlope)
{
    const double derivative = SCDAT::Coupling::estimateCenteredDidv(
        [](double potential_v) { return 3.0 * potential_v - 2.0; }, 4.0);
    EXPECT_NEAR(derivative, 3.0, 1.0e-12);
}

TEST(CouplingTest, BisectionRootSolverFindsBracketedCurrentRoot)
{
    const auto result = SCDAT::Coupling::solveBisectionCurrentRoot(
        [](double potential_v) { return potential_v - 1.25; }, -5.0, 5.0);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.root_x, 1.25, 1.0e-6);
    EXPECT_NEAR(result.value_at_root, 0.0, 1.0e-6);
}

TEST(CouplingTest, BisectionRootSolverMidpointFallbackKeepsResultBounded)
{
    const auto result = SCDAT::Coupling::solveBisectionCurrentRoot(
        [](double potential_v) { return potential_v * potential_v - 0.25; }, -5.0, -1.0);
    EXPECT_TRUE(result.converged);
    EXPECT_GE(result.root_x, -5.0);
    EXPECT_LE(result.root_x, -1.0);
    EXPECT_TRUE(std::isfinite(result.value_at_root));
}

TEST(CouplingTest, AggregateSurfaceCurrentDensityIncludesExtraTerms)
{
    SCDAT::Coupling::SurfaceCurrentComponentState components;
    components.electron_collection_a = -4.0e-9;
    components.ion_collection_a = 1.0e-9;
    components.secondary_electron_emission_a = 2.0e-10;
    components.ion_secondary_electron_emission_a = 1.0e-10;
    components.backscatter_emission_a = 5.0e-11;
    components.photoelectron_emission_a = 3.0e-10;
    components.conduction_a = 8.0e-10;

    const double total = SCDAT::Coupling::aggregateSurfaceCurrentDensity(
        components, 4.0e-10, 2.0e-10, 1.5e-9);
    EXPECT_NEAR(total, -7.5e-10, 1.0e-21);
}

TEST(CouplingTest, CalibrationFactorEstimatorClampsAndPreservesFallback)
{
    EXPECT_NEAR(SCDAT::Coupling::estimateCalibrationFactor(6.0, 2.0, 0.5, 2.5, 1.0), 2.5,
                1.0e-12);
    EXPECT_NEAR(SCDAT::Coupling::estimateCalibrationFactor(1.0, 0.0, 0.5, 2.5, 1.3), 1.3,
                1.0e-12);
}

TEST(CouplingTest, AdaptiveSubstepCountRespondsToStiffnessAndNegativeBias)
{
    SurfaceAdaptiveSubstepInputs baseline_inputs;
    baseline_inputs.base_substeps = 4;
    baseline_inputs.reference_potential_v = 0.0;
    baseline_inputs.dt_s = 1.0;
    baseline_inputs.current_derivative_a_per_m2_per_v = 2.0e-12;
    baseline_inputs.capacitance_per_area_f_per_m2 = 1.0e-12;

    const std::size_t baseline =
        SCDAT::Coupling::computeSurfaceAdaptiveSubstepCount(baseline_inputs);
    EXPECT_EQ(baseline, 8u);

    SurfaceAdaptiveSubstepInputs negative_bias_inputs = baseline_inputs;
    negative_bias_inputs.reference_potential_v = -250.0;
    const std::size_t negative_bias =
        SCDAT::Coupling::computeSurfaceAdaptiveSubstepCount(negative_bias_inputs);
    EXPECT_GT(negative_bias, baseline);
    EXPECT_EQ(negative_bias, 10u);
}

TEST(CouplingTest, AdaptiveSubstepCountHonorsConfiguredMaximum)
{
    SurfaceAdaptiveSubstepInputs inputs;
    inputs.base_substeps = 3;
    inputs.maximum_substeps = 40;
    inputs.reference_potential_v = -500.0;
    inputs.dt_s = 10.0;
    inputs.current_derivative_a_per_m2_per_v = 1.0e+2;
    inputs.capacitance_per_area_f_per_m2 = 1.0e-12;

    const std::size_t substeps =
        SCDAT::Coupling::computeSurfaceAdaptiveSubstepCount(inputs);
    EXPECT_EQ(substeps, 40u);
}

TEST(CouplingTest, SheathEquivalentCapacitanceMatchesSeriesCombination)
{
    SurfaceSheathCapacitanceConsistencyInputs inputs;
    inputs.dielectric_capacitance_per_area_f_per_m2 = 1.0e-8;
    inputs.effective_sheath_length_m = 1.0e-4;
    inputs.minimum_sheath_length_m = 1.0e-6;
    inputs.maximum_sheath_length_m = 1.0e-4;
    inputs.relative_permittivity = 2.0;

    const double equivalent =
        SCDAT::Coupling::computeSurfaceSheathEquivalentCapacitancePerArea(inputs);
    const double sheath_capacitance = 8.8541878128e-12 * 2.0 / 1.0e-4;
    const double expected = 1.0 / (1.0 / 1.0e-8 + 1.0 / sheath_capacitance);

    EXPECT_NEAR(equivalent, expected, 1.0e-20);
    EXPECT_LT(equivalent, inputs.dielectric_capacitance_per_area_f_per_m2);
}

TEST(CouplingTest, SheathCapacitanceConsistencyBlendRespectsRatioGuard)
{
    SurfaceSheathCapacitanceConsistencyInputs inputs;
    inputs.dielectric_capacitance_per_area_f_per_m2 = 1.0e-8;
    inputs.effective_sheath_length_m = 0.5;
    inputs.minimum_sheath_length_m = 1.0e-6;
    inputs.maximum_sheath_length_m = 1.0;
    inputs.relative_permittivity = 1.0;
    inputs.consistency_weight = 0.95;
    inputs.volume_mesh_coupling_gain = 1.0;
    inputs.ratio_guard = 0.1;

    const double blended =
        SCDAT::Coupling::enforceSurfaceSheathCapacitanceConsistency(inputs);
    const double min_allowed =
        inputs.dielectric_capacitance_per_area_f_per_m2 / (1.0 + 0.1);
    EXPECT_NEAR(blended, min_allowed, 1.0e-18);
}

TEST(CouplingTest, SurfaceCalibrationFactorsEstimateElectronAndIonTogether)
{
    const auto factors = SCDAT::Coupling::estimateSurfaceCalibrationFactors(
        SCDAT::Coupling::SurfaceCalibrationFactorInputs{
            4.0, 2.0,
            0.0, 0.0,
            0.5, 3.0,
            1.2, 1.4,
        });
    EXPECT_NEAR(factors.electron_factor, 2.0, 1.0e-12);
    EXPECT_NEAR(factors.ion_factor, 1.4, 1.0e-12);
}

TEST(CouplingTest, SurfaceEmissionBlendSupportsDirectAndHybridModes)
{
    const auto direct = SCDAT::Coupling::blendSurfaceEmissionComponents(
        SCDAT::Coupling::SurfaceEmissionBlendInputs{
            2.0e-10, 1.0e-10, 5.0e-11, 7.0e-10, false,
        });
    EXPECT_NEAR(direct.secondary_emission_a, 7.0e-10, 1.0e-24);
    EXPECT_NEAR(direct.ion_secondary_emission_a, 0.0, 1.0e-24);
    EXPECT_NEAR(direct.backscatter_emission_a, 0.0, 1.0e-24);

    const auto hybrid = SCDAT::Coupling::blendSurfaceEmissionComponents(
        SCDAT::Coupling::SurfaceEmissionBlendInputs{
            2.0e-10, 1.0e-10, 5.0e-11, 7.0e-10, true,
        });
    EXPECT_NEAR(hybrid.emission_scale, 2.0, 1.0e-12);
    EXPECT_NEAR(hybrid.secondary_emission_a, 4.0e-10, 1.0e-24);
    EXPECT_NEAR(hybrid.ion_secondary_emission_a, 2.0e-10, 1.0e-24);
    EXPECT_NEAR(hybrid.backscatter_emission_a, 1.0e-10, 1.0e-24);
}

TEST(CouplingTest, SurfaceDistributionBlendComputesCollectionScales)
{
    const auto blend = SCDAT::Coupling::blendSurfaceDistributionComponents(
        SCDAT::Coupling::SurfaceDistributionBlendInputs{
            -2.0e-9, 1.0e-9,
            -5.0e-9, 2.5e-9,
            12.0, 4.0,
        });
    EXPECT_NEAR(blend.electron_collection_scale, 2.5, 1.0e-12);
    EXPECT_NEAR(blend.ion_collection_scale, 2.5, 1.0e-12);
    EXPECT_NEAR(blend.electron_characteristic_energy_ev, 12.0, 1.0e-12);
    EXPECT_NEAR(blend.ion_characteristic_energy_ev, 4.0, 1.0e-12);
}

TEST(CouplingTest, DensityLinearizationAggregatesCurrentAndDidvWithExtraTerms)
{
    SCDAT::Coupling::SurfaceCurrentDensityLinearizationInputs inputs;
    inputs.component_densities.electron_collection_a = -4.0e-9;
    inputs.component_densities.ion_collection_a = 1.2e-9;
    inputs.component_densities.secondary_electron_emission_a = 3.0e-10;
    inputs.component_densities.backscatter_emission_a = 1.0e-10;
    inputs.component_densities.photoelectron_emission_a = 2.0e-10;
    inputs.component_densities.conduction_a = 9.0e-10;
    inputs.area_m2 = 2.0;
    inputs.plasma_didv_density_a_per_m2_per_v = -5.0e-11;
    inputs.conduction_didv_density_a_per_m2_per_v = 2.0e-11;
    inputs.thermionic_emission_density_a_per_m2 = 4.0e-10;
    inputs.field_emission_density_a_per_m2 = 1.0e-10;
    inputs.ram_ion_current_density_a_per_m2 = 3.0e-10;

    const auto result = SCDAT::Coupling::linearizeSurfaceCurrentDensity(inputs);
    EXPECT_NEAR(result.plasma_current_a, -5.0e-9, 1.0e-21);
    EXPECT_NEAR(result.total_current_density_a_per_m2, -1.7e-9, 1.0e-21);
    EXPECT_NEAR(result.total_current_a, -3.4e-9, 1.0e-21);
    EXPECT_NEAR(result.total_didv_a_per_v, -6.0e-11, 1.0e-24);
    EXPECT_NEAR(result.total_didv_density_a_per_m2_per_v, -3.0e-11, 1.0e-24);
}

TEST(CouplingTest, SurfacePicRouteAdjustmentAppliesKernelSnapshotForFormalRoute)
{
    SCDAT::Coupling::SurfaceCurrentComponentState base_components;
    base_components.electron_collection_a = -2.0e-9;
    base_components.ion_collection_a = 1.0e-9;

    SCDAT::Coupling::SurfaceCurrentComponentState snapshot_components;
    snapshot_components.electron_collection_a = -6.0e-9;
    snapshot_components.ion_collection_a = 2.5e-9;
    snapshot_components.secondary_electron_emission_a = 8.0e-10;

    const auto adjusted = SCDAT::Coupling::applySurfacePicRouteAdjustment(
        SCDAT::Coupling::SurfacePicRouteAdjustmentInputs{
            true,
            true,
            snapshot_components,
            -4.2e-11,
        },
        base_components, -1.0e-11);

    EXPECT_TRUE(adjusted.applied);
    EXPECT_NEAR(adjusted.components.electron_collection_a, -6.0e-9, 1.0e-24);
    EXPECT_NEAR(adjusted.components.ion_collection_a, 2.5e-9, 1.0e-24);
    EXPECT_NEAR(adjusted.components.secondary_electron_emission_a, 8.0e-10, 1.0e-24);
    EXPECT_NEAR(adjusted.total_didv_a_per_v, -4.2e-11, 1.0e-24);
}

TEST(CouplingTest, SurfacePicRouteAdjustmentKeepsCurrentStateWhenNotApplicable)
{
    SCDAT::Coupling::SurfaceCurrentComponentState base_components;
    base_components.electron_collection_a = -3.0e-9;
    base_components.ion_collection_a = 1.5e-9;

    const auto adjusted = SCDAT::Coupling::applySurfacePicRouteAdjustment(
        SCDAT::Coupling::SurfacePicRouteAdjustmentInputs{},
        base_components, -1.4e-11);

    EXPECT_FALSE(adjusted.applied);
    EXPECT_NEAR(adjusted.components.electron_collection_a, -3.0e-9, 1.0e-24);
    EXPECT_NEAR(adjusted.components.ion_collection_a, 1.5e-9, 1.0e-24);
    EXPECT_NEAR(adjusted.total_didv_a_per_v, -1.4e-11, 1.0e-24);
}

TEST(CouplingTest, SurfacePicRouteAdjustmentRequiresValidSnapshotParityWithFormalRoute)
{
    SCDAT::Coupling::SurfaceCurrentComponentState base_components;
    base_components.electron_collection_a = -2.2e-9;
    base_components.ion_collection_a = 9.0e-10;
    base_components.secondary_electron_emission_a = 2.0e-10;

    SCDAT::Coupling::SurfaceCurrentComponentState snapshot_components;
    snapshot_components.electron_collection_a = -9.0e-9;
    snapshot_components.ion_collection_a = 4.0e-9;
    snapshot_components.secondary_electron_emission_a = 7.0e-10;

    const auto formal_without_snapshot = SCDAT::Coupling::applySurfacePicRouteAdjustment(
        SCDAT::Coupling::SurfacePicRouteAdjustmentInputs{
            true,
            false,
            snapshot_components,
            -5.0e-11,
        },
        base_components, -1.3e-11);
    EXPECT_FALSE(formal_without_snapshot.applied);
    EXPECT_NEAR(formal_without_snapshot.components.electron_collection_a, -2.2e-9, 1.0e-24);
    EXPECT_NEAR(formal_without_snapshot.total_didv_a_per_v, -1.3e-11, 1.0e-24);

    const auto non_formal_with_snapshot = SCDAT::Coupling::applySurfacePicRouteAdjustment(
        SCDAT::Coupling::SurfacePicRouteAdjustmentInputs{
            false,
            true,
            snapshot_components,
            -5.0e-11,
        },
        base_components, -1.3e-11);
    EXPECT_FALSE(non_formal_with_snapshot.applied);
    EXPECT_NEAR(non_formal_with_snapshot.components.electron_collection_a, -2.2e-9, 1.0e-24);
    EXPECT_NEAR(non_formal_with_snapshot.total_didv_a_per_v, -1.3e-11, 1.0e-24);
}

TEST(CouplingTest, SharedCurrentMatrixFallbackBuildsWeightedDidvAndOffDiagonalEntries)
{
    std::vector<double> shared_patch_voltage_weights(6, 0.0);
    shared_patch_voltage_weights[1] = 0.4;
    shared_patch_voltage_weights[2] = 0.3;
    shared_patch_voltage_weights[3] = 0.2;

    std::vector<std::size_t> patch_node_indices = {1, 2, 3};
    std::vector<std::size_t> circuit_to_reduced_node_index(
        6, std::numeric_limits<std::size_t>::max());
    circuit_to_reduced_node_index[1] = 0;
    circuit_to_reduced_node_index[2] = 1;
    circuit_to_reduced_node_index[3] = 0;

    const auto result = SCDAT::Coupling::computeSurfaceSharedCurrentMatrixFallback(
        true,
        true,
        1,
        -10.0,
        2.0,
        shared_patch_voltage_weights,
        patch_node_indices,
        circuit_to_reduced_node_index,
        0.0);

    EXPECT_NEAR(result.plasma_didv_a_per_v, -2.0, 1.0e-12);
    ASSERT_EQ(result.off_diagonal_entries.size(), 1U);
    EXPECT_EQ(result.off_diagonal_entries[0].row_node, 1U);
    EXPECT_EQ(result.off_diagonal_entries[0].column_node, 2U);
    EXPECT_NEAR(result.off_diagonal_entries[0].coefficient_a_per_v, 3.0, 1.0e-12);
}

TEST(CouplingTest, SharedCurrentMatrixFallbackKeepsCurrentStateWhenDisabledOrInvalid)
{
    std::vector<double> shared_patch_voltage_weights(4, 0.0);
    shared_patch_voltage_weights[1] = 0.25;
    shared_patch_voltage_weights[2] = 0.75;
    std::vector<std::size_t> patch_node_indices = {1, 2};
    std::vector<std::size_t> circuit_to_reduced_node_index(
        4, std::numeric_limits<std::size_t>::max());

    const auto disabled = SCDAT::Coupling::computeSurfaceSharedCurrentMatrixFallback(
        false,
        true,
        1,
        -8.0,
        1.0,
        shared_patch_voltage_weights,
        patch_node_indices,
        circuit_to_reduced_node_index,
        0.0);
    EXPECT_NEAR(disabled.plasma_didv_a_per_v, -7.0, 1.0e-12);
    EXPECT_TRUE(disabled.off_diagonal_entries.empty());

    const auto invalid = SCDAT::Coupling::computeSurfaceSharedCurrentMatrixFallback(
        true,
        false,
        1,
        -8.0,
        1.0,
        shared_patch_voltage_weights,
        patch_node_indices,
        circuit_to_reduced_node_index,
        0.0);
    EXPECT_NEAR(invalid.plasma_didv_a_per_v, 0.0, 1.0e-12);
    EXPECT_TRUE(invalid.off_diagonal_entries.empty());
}

TEST(CouplingTest, IndirectPathBranchWeightAccumulatesTwoHopHarmonicPairs)
{
    std::vector<std::vector<double>> adjacency(4, std::vector<double>(4, 0.0));
    adjacency[0][1] = 2.0;
    adjacency[1][0] = 2.0;
    adjacency[1][3] = 4.0;
    adjacency[3][1] = 4.0;
    adjacency[0][2] = 1.0;
    adjacency[2][0] = 1.0;
    adjacency[2][3] = 1.0;
    adjacency[3][2] = 1.0;

    const double weight =
        SCDAT::Coupling::computeIndirectPathBranchWeight(adjacency, 0, 3);
    EXPECT_NEAR(weight, 11.0 / 3.0, 1.0e-12);
    EXPECT_NEAR(
        SCDAT::Coupling::computeIndirectPathBranchWeight(adjacency, 1, 1), 0.0,
        1.0e-12);
}

TEST(CouplingTest, HarmonicPairCapacitanceHelperMatchesAnalyticalValue)
{
    constexpr double kEpsilon0 = 8.8541878128e-12;
    const double area_i_m2 = 2.0;
    const double area_j_m2 = 3.0;
    const double sheath_length_m = 0.2;
    const double scale = 0.5;
    const double harmonic_pair_area_m2 =
        2.0 * area_i_m2 * area_j_m2 / (area_i_m2 + area_j_m2);
    const double expected = scale * kEpsilon0 * harmonic_pair_area_m2 / sheath_length_m;

    const double actual = SCDAT::Coupling::computeHarmonicPairCapacitanceF(
        area_i_m2, area_j_m2, sheath_length_m, scale, 1.0e-18);
    EXPECT_NEAR(actual, expected, std::abs(expected) * 1.0e-12 + 1.0e-24);
}

TEST(CouplingTest, HarmonicPairCapacitanceHelperAppliesSafetyClamps)
{
    const double clamped = SCDAT::Coupling::computeHarmonicPairCapacitanceF(
        0.0, 0.0, 0.0, -1.0, 7.5e-19);
    EXPECT_NEAR(clamped, 7.5e-19, 1.0e-30);
}
