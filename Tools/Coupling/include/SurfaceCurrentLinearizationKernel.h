#pragma once

#include "SurfaceCircuitCoupling.h"

#include <functional>
#include <vector>

namespace SCDAT
{
namespace Coupling
{

struct SurfaceCurrentComponentState
{
    double electron_collection_a = 0.0;
    double ion_collection_a = 0.0;
    double secondary_electron_emission_a = 0.0;
    double ion_secondary_electron_emission_a = 0.0;
    double backscatter_emission_a = 0.0;
    double photoelectron_emission_a = 0.0;
    double conduction_a = 0.0;
    double electron_collection_didv_a_per_v = 0.0;
    double ion_collection_didv_a_per_v = 0.0;
    double secondary_electron_emission_didv_a_per_v = 0.0;
    double ion_secondary_electron_emission_didv_a_per_v = 0.0;
    double backscatter_emission_didv_a_per_v = 0.0;
    double photoelectron_emission_didv_a_per_v = 0.0;
    double conduction_didv_a_per_v = 0.0;
};

struct SurfaceCurrentLinearizationNode
{
    SurfaceCurrentComponentState components;
    double additional_rhs_a = 0.0;
    double additional_diagonal_a_per_v = 0.0;
};

struct SurfaceCurrentLinearizationResult
{
    std::vector<double> total_current_a;
    std::vector<double> total_didv_a_per_v;
    std::vector<SurfaceCurrentComponentState> component_states;
    SurfaceCircuitKernelInput kernel_input;
};

struct SurfaceCurrentRootSolveResult
{
    double root_x = 0.0;
    double value_at_root = 0.0;
    bool converged = false;
};

struct SurfaceCurrentDensityLinearizationInputs
{
    SurfaceCurrentComponentState component_densities;
    double area_m2 = 1.0;
    double plasma_didv_density_a_per_m2_per_v = 0.0;
    double conduction_didv_density_a_per_m2_per_v = 0.0;
    double thermionic_emission_density_a_per_m2 = 0.0;
    double field_emission_density_a_per_m2 = 0.0;
    double ram_ion_current_density_a_per_m2 = 0.0;
};

struct SurfaceCurrentDensityLinearizationResult
{
    double total_current_density_a_per_m2 = 0.0;
    double total_didv_density_a_per_m2_per_v = 0.0;
    double plasma_current_a = 0.0;
    double plasma_didv_a_per_v = 0.0;
    double conduction_current_a = 0.0;
    double conduction_didv_a_per_v = 0.0;
    double total_current_a = 0.0;
    double total_didv_a_per_v = 0.0;
    SurfaceCurrentLinearizationResult linearization;
};

struct SurfaceCalibrationFactorInputs
{
    double measured_electron_current_a = 0.0;
    double native_electron_current_a = 0.0;
    double measured_ion_current_a = 0.0;
    double native_ion_current_a = 0.0;
    double minimum_factor = 0.0;
    double maximum_factor = 1.0;
    double default_electron_factor = 1.0;
    double default_ion_factor = 1.0;
};

struct SurfaceCalibrationFactors
{
    double electron_factor = 1.0;
    double ion_factor = 1.0;
};

struct SurfaceEmissionBlendInputs
{
    double reference_secondary_emission_a = 0.0;
    double reference_ion_secondary_emission_a = 0.0;
    double reference_backscatter_emission_a = 0.0;
    double measured_emitted_total_a = 0.0;
    bool hybrid_completion = false;
};

struct SurfaceEmissionBlendResult
{
    double secondary_emission_a = 0.0;
    double ion_secondary_emission_a = 0.0;
    double backscatter_emission_a = 0.0;
    double reference_emitted_total_a = 0.0;
    double emission_scale = 1.0;
};

struct SurfaceDistributionBlendInputs
{
    double reference_electron_current_a = 0.0;
    double reference_ion_current_a = 0.0;
    double measured_electron_current_a = 0.0;
    double measured_ion_current_a = 0.0;
    double electron_characteristic_energy_ev = 1.0;
    double ion_characteristic_energy_ev = 1.0;
};

struct SurfaceDistributionBlendResult
{
    double electron_collection_scale = 1.0;
    double ion_collection_scale = 1.0;
    double electron_characteristic_energy_ev = 1.0;
    double ion_characteristic_energy_ev = 1.0;
};

struct SurfacePicRouteAdjustmentInputs
{
    bool formal_surface_pic_route = false;
    bool kernel_snapshot_valid = false;
    SurfaceCurrentComponentState snapshot_components;
    double snapshot_total_didv_a_per_v = 0.0;
};

struct SurfacePicRouteAdjustmentResult
{
    SurfaceCurrentComponentState components;
    double total_didv_a_per_v = 0.0;
    bool applied = false;
};

struct SurfaceSharedCurrentMatrixFallbackResult
{
    double plasma_didv_a_per_v = 0.0;
    std::vector<SurfaceCircuitLinearization::OffDiagonalEntry> off_diagonal_entries;
};

struct SurfaceAdaptiveSubstepInputs
{
    std::size_t base_substeps = 1;
    double reference_potential_v = 0.0;
    double dt_s = 0.0;
    double current_derivative_a_per_m2_per_v = 0.0;
    double capacitance_per_area_f_per_m2 = 1.0e-12;
    std::size_t maximum_substeps = 256;
};

struct SurfaceSheathCapacitanceConsistencyInputs
{
    double dielectric_capacitance_per_area_f_per_m2 = 1.0e-12;
    double effective_sheath_length_m = 1.0e-3;
    double minimum_sheath_length_m = 1.0e-6;
    double maximum_sheath_length_m = 1.0e-3;
    double relative_permittivity = 1.0;
    double consistency_weight = 0.60;
    double volume_mesh_coupling_gain = 0.0;
    double ratio_guard = 1.5;
    double minimum_capacitance_per_area_f_per_m2 = 1.0e-12;
};

class SurfaceCurrentLinearizationKernel
{
  public:
    SurfaceCurrentLinearizationResult assemble(
        const std::vector<SurfaceCurrentLinearizationNode>& nodes,
        const std::vector<SurfaceCircuitLinearization::OffDiagonalEntry>& off_diagonal_entries = {}) const;
};

double estimateCenteredDidv(const std::function<double(double)>& evaluator, double sample_x,
                            double minimum_step = 1.0e-3, double relative_step = 1.0e-2);
SurfaceCurrentRootSolveResult solveBisectionCurrentRoot(
    const std::function<double(double)>& evaluator, double minimum_x, double maximum_x,
    double function_tolerance = 1.0e-12, double bracket_tolerance = 1.0e-6,
    int max_iterations = 80);
double aggregateSurfaceCurrentDensity(const SurfaceCurrentComponentState& components,
                                      double thermionic_emission_a = 0.0,
                                      double field_emission_a = 0.0,
                                      double ram_ion_current_a = 0.0);
SurfaceCurrentDensityLinearizationResult linearizeSurfaceCurrentDensity(
    const SurfaceCurrentDensityLinearizationInputs& inputs);
SurfaceCalibrationFactors estimateSurfaceCalibrationFactors(
    const SurfaceCalibrationFactorInputs& inputs);
SurfaceEmissionBlendResult blendSurfaceEmissionComponents(
    const SurfaceEmissionBlendInputs& inputs);
SurfaceDistributionBlendResult blendSurfaceDistributionComponents(
    const SurfaceDistributionBlendInputs& inputs);
SurfaceSharedCurrentMatrixFallbackResult computeSurfaceSharedCurrentMatrixFallback(
    bool fallback_shared_current_matrix_coupling_enabled,
    bool patch_didv_valid,
    std::size_t patch_node_index,
    double patch_plasma_didv_a_per_v,
    double patch_conduction_didv_a_per_v,
    const std::vector<double>& shared_patch_voltage_weights,
    const std::vector<std::size_t>& patch_node_indices,
    const std::vector<std::size_t>& circuit_to_reduced_node_index,
    double minimum_coupled_derivative_abs_a_per_v = 0.0);
SurfacePicRouteAdjustmentResult applySurfacePicRouteAdjustment(
    const SurfacePicRouteAdjustmentInputs& inputs,
    const SurfaceCurrentComponentState& current_components,
    double current_total_didv_a_per_v = 0.0);
std::size_t computeSurfaceAdaptiveSubstepCount(
    const SurfaceAdaptiveSubstepInputs& inputs);
double computeSurfaceSheathEquivalentCapacitancePerArea(
    const SurfaceSheathCapacitanceConsistencyInputs& inputs);
double enforceSurfaceSheathCapacitanceConsistency(
    const SurfaceSheathCapacitanceConsistencyInputs& inputs);
double computeIndirectPathBranchWeight(
    const std::vector<std::vector<double>>& branch_graph_adjacency_weight_f,
    std::size_t from_node, std::size_t to_node);
double computeHarmonicPairCapacitanceF(double area_i_m2, double area_j_m2,
                                       double effective_sheath_length_m,
                                       double capacitance_scale = 1.0,
                                       double minimum_capacitance_f = 1.0e-18);
double estimateCalibrationFactor(double measured_current_a, double native_current_a,
                                 double minimum_factor, double maximum_factor,
                                 double default_factor = 1.0);

} // namespace Coupling
} // namespace SCDAT
