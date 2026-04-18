#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Coupling
{

struct SurfaceCircuitNode
{
    std::string name;
    double potential_v = 0.0;
    double capacitance_f = 0.0;
    double shunt_conductance_s = 0.0;
    bool fixed_potential = false;
    double fixed_value_v = 0.0;
};

struct SurfaceCircuitBranch
{
    std::size_t from_node = 0;
    std::size_t to_node = 0;
    double conductance_s = 0.0;
    double bias_v = 0.0;
};

struct SurfaceCircuitLinearization
{
    std::vector<double> plasma_current_a;
    std::vector<double> plasma_didv_a_per_v;
    std::vector<double> additional_rhs_a;
    std::vector<double> additional_diagonal_a_per_v;

    struct OffDiagonalEntry
    {
        std::size_t row_node = 0;
        std::size_t column_node = 0;
        double coefficient_a_per_v = 0.0;
    };
    std::vector<OffDiagonalEntry> additional_off_diagonal_entries;
};

enum class CircuitElementKind
{
    NodeCapacitance,
    NodeShuntConductance,
    NodeFixedPotential,
    BranchConductance,
    BranchBias,
    BranchCurrentSource,
};

struct CircuitNodeDescriptor
{
    std::string name;
    double initial_potential_v = 0.0;
};

struct CircuitBranchDescriptor
{
    std::string name;
    std::size_t from_node = 0;
    std::size_t to_node = 0;
};

struct CircuitElementDescriptor
{
    CircuitElementKind kind = CircuitElementKind::NodeCapacitance;
    std::size_t target_index = 0;
    double value = 0.0;
};

struct CircuitAssembly
{
    std::vector<CircuitNodeDescriptor> nodes;
    std::vector<CircuitBranchDescriptor> branches;
    std::vector<CircuitElementDescriptor> elements;
};

struct SurfaceCircuitKernelNodeState
{
    double plasma_current_a = 0.0;
    double plasma_didv_a_per_v = 0.0;
    double additional_rhs_a = 0.0;
    double additional_diagonal_a_per_v = 0.0;
    double capacitance_f = std::numeric_limits<double>::quiet_NaN();
    double shunt_conductance_s = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct SurfaceCircuitKernelBranchState
{
    double conductance_s = std::numeric_limits<double>::quiet_NaN();
    double bias_v = std::numeric_limits<double>::quiet_NaN();
    double current_source_a = 0.0;
    bool valid = false;
};

struct SurfaceCircuitKernelInput
{
    std::vector<SurfaceCircuitKernelNodeState> nodes;
    std::vector<SurfaceCircuitKernelBranchState> branches;
    std::vector<SurfaceCircuitLinearization::OffDiagonalEntry> off_diagonal_entries;
};

struct CircuitDidvTerm
{
    std::string name;
    double current_a = 0.0;
    double didv_a_per_v = 0.0;
};

enum class CircuitDidvAggregationKind
{
    DirectSum,
    WeightedTerms,
    MultiSurface,
};

struct CircuitWeightedDidvTerm
{
    std::string name;
    double weight = 1.0;
    double current_a = 0.0;
    double didv_a_per_v = 0.0;
};

struct CircuitDidvSurfaceComposition
{
    std::string name;
    double surface_weight = 1.0;
    std::vector<CircuitDidvTerm> terms;
};

struct CircuitDidvNodeComposition
{
    std::vector<CircuitDidvTerm> terms;
    CircuitDidvAggregationKind aggregation_kind = CircuitDidvAggregationKind::DirectSum;
    std::vector<CircuitWeightedDidvTerm> weighted_terms;
    std::vector<CircuitDidvSurfaceComposition> surface_compositions;
    double additional_rhs_a = 0.0;
    double additional_diagonal_a_per_v = 0.0;
};

struct CircuitDidvCompositionInput
{
    std::vector<CircuitDidvNodeComposition> nodes;
    std::vector<SurfaceCircuitLinearization::OffDiagonalEntry> off_diagonal_entries;
};

struct CircuitDidvNodeSummary
{
    double total_current_a = 0.0;
    double total_didv_a_per_v = 0.0;
};

enum class CircuitDidvRouteKind
{
    DirectComposition,
    SurfMatrices,
    ReducedSurfDidv,
    ReducedRegularDidv,
};

struct CircuitSurfDidvMatrixTerm
{
    std::string name;
    std::size_t node_index = 0;
    double current_a = 0.0;
    double didv_a_per_v = 0.0;
};

struct CircuitSurfDidvMatrixInput
{
    std::size_t node_count = 0;
    std::vector<CircuitSurfDidvMatrixTerm> terms;
    std::vector<SurfaceCircuitLinearization::OffDiagonalEntry> off_diagonal_entries;
};

struct ReducableCircuitDidv
{
    CircuitDidvRouteKind route_kind = CircuitDidvRouteKind::DirectComposition;
    CircuitDidvCompositionInput composition;
    std::vector<std::size_t> node_reduction_map;
    std::size_t reduced_node_count = 0;
};

struct CircuitDidvCompositionResult
{
    std::vector<CircuitDidvNodeSummary> node_summaries;
    SurfaceCircuitLinearization linearization;
    SurfaceCircuitKernelInput kernel_input;
};

CircuitDidvCompositionResult composeCircuitDidv(const CircuitDidvCompositionInput& input);
CircuitDidvCompositionInput composeSurfDidvFromMatrices(
    const CircuitSurfDidvMatrixInput& input);
ReducableCircuitDidv makeReducableCircuitDidv(const CircuitDidvCompositionInput& composition);
CircuitDidvCompositionInput reduceCircuitDidv(
    const ReducableCircuitDidv& reducable_didv);
CircuitDidvCompositionInput reduceCircuitDidvFromSurfDidv(
    const CircuitSurfDidvMatrixInput& input,
    const std::vector<std::size_t>& node_reduction_map,
    std::size_t reduced_node_count);
CircuitDidvCompositionInput reduceCircuitDidvFromRegularDidv(
    const CircuitDidvCompositionInput& input,
    const std::vector<std::size_t>& node_reduction_map,
    std::size_t reduced_node_count);

struct CircuitKernelTopologyState
{
    std::vector<double> node_capacitance_f;
    std::vector<double> node_shunt_conductance_s;
    std::vector<double> branch_conductance_s;
    std::vector<double> branch_bias_v;
    std::vector<double> branch_current_source_a;
};

SurfaceCircuitKernelInput composeCircuitKernelInput(
    const CircuitDidvCompositionInput& didv_input,
    const CircuitKernelTopologyState& topology_state);

enum class CircuitExcitationWaveformKind
{
    Constant,
    Sinusoidal,
    Pulse,
    PiecewiseLinear,
    Exponential,
};

struct CircuitExcitationPwlPoint
{
    double time_s = 0.0;
    double value = 0.0;
};

struct CircuitExcitationDescriptor
{
    std::string name;
    CircuitElementKind target_kind = CircuitElementKind::BranchCurrentSource;
    std::size_t target_index = 0;
    CircuitExcitationWaveformKind waveform = CircuitExcitationWaveformKind::Constant;
    double base_value = 0.0;
    double amplitude = 0.0;
    double frequency_hz = 0.0;
    double phase_rad = 0.0;
    double duty_cycle = 0.5;
    double start_time_s = -std::numeric_limits<double>::infinity();
    double end_time_s = std::numeric_limits<double>::infinity();
    bool enabled = true;
    std::vector<CircuitExcitationPwlPoint> pwl_points;
    double exponential_tau_s = 0.0;
};

struct CircuitExcitationSampleResult
{
    std::vector<CircuitElementDescriptor> elements;
};

CircuitExcitationSampleResult sampleCircuitExcitations(
    const std::vector<CircuitExcitationDescriptor>& excitations, double time_s);

struct SurfaceCircuitAssemblyResult
{
    std::vector<SurfaceCircuitNode> nodes;
    std::vector<SurfaceCircuitBranch> branches;
    SurfaceCircuitKernelInput kernel_input;
};

struct CircuitFamilyRouteInput
{
    bool has_circuit_model = false;
    bool enable_reference_current_balance = false;
    bool enable_body_patch_circuit = false;
    bool has_dynamic_excitations = false;
    bool has_weighted_didv_terms = false;
    bool has_multi_surface_didv = false;
    bool has_matrix_didv = false;
    bool has_regular_reduction = false;
    bool has_surface_reduction = false;
    bool has_off_diagonal_coupling = false;
    std::size_t node_count = 0;
    std::size_t patch_count = 0;
    std::size_t branch_count = 0;
};

struct CircuitFamilyView
{
    std::size_t circ_family_count = 0;
    std::size_t active_circ_family_count = 0;
    std::string circ_family_signature;
    std::string active_circ_family_signature;
    std::size_t circfield_family_count = 0;
    std::size_t active_circfield_family_count = 0;
    std::string circfield_family_signature;
    std::string active_circfield_family_signature;
    std::size_t didv_family_count = 0;
    std::size_t active_didv_family_count = 0;
    std::string didv_family_signature;
    std::string active_didv_family_signature;
};

CircuitFamilyView resolveCircuitFamilyView(const CircuitFamilyRouteInput& input);

SurfaceCircuitAssemblyResult assembleSurfaceCircuit(const CircuitAssembly& assembly);
SurfaceCircuitAssemblyResult assembleSurfaceCircuit(
    const CircuitAssembly& assembly,
    const std::vector<CircuitElementDescriptor>& dynamic_elements);

struct SurfaceCircuitAdvanceResult
{
    std::vector<double> node_potentials_v;
    std::vector<double> plasma_currents_a;
    std::vector<double> branch_currents_a;
    double max_delta_potential_v = 0.0;
    bool converged = false;
};

class SurfaceCircuitCoupling
{
  public:
    void clear();
    void configure(const CircuitAssembly& assembly);
    std::size_t addNode(const SurfaceCircuitNode& node);
    void addBranch(const SurfaceCircuitBranch& branch);

    void setNodePotential(std::size_t node_index, double potential_v);
    void setNodeFixedPotential(std::size_t node_index, bool fixed, double fixed_value_v);
    void setNodeShuntConductance(std::size_t node_index, double conductance_s);
    void setNodeCapacitance(std::size_t node_index, double capacitance_f);
    void setBranchConductance(std::size_t branch_index, double conductance_s);

    const std::vector<SurfaceCircuitNode>& getNodes() const
    {
        return nodes_;
    }

    const std::vector<SurfaceCircuitBranch>& getBranches() const
    {
        return branches_;
    }

    SurfaceCircuitAdvanceResult advanceImplicit(
        double dt, const SurfaceCircuitLinearization& linearization,
        double max_delta_potential_v = std::numeric_limits<double>::infinity());
    SurfaceCircuitAdvanceResult advanceImplicit(
        double dt, const SurfaceCircuitKernelInput& kernel_input,
        double max_delta_potential_v = std::numeric_limits<double>::infinity());

  private:
    static std::vector<double> solveDenseLinearSystem(std::vector<std::vector<double>> matrix,
                                                      std::vector<double> rhs);

    std::vector<SurfaceCircuitNode> nodes_;
    std::vector<SurfaceCircuitBranch> branches_;
};

} // namespace Coupling
} // namespace SCDAT
