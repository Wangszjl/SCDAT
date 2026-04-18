#pragma once

#include <cstddef>
#include <string>

namespace SCDAT
{
namespace FieldSolver
{

enum class NativeVolumeParityRoute
{
    NativeMinimal = 0,
    ExternalBridge = 1,
    HybridBlend = 2,
};

enum class NativeVolumeBoundaryConditionFamily
{
    VoltageDependentMBC = 0,
    MixedDirichletFourierPoissonBC = 1,
    FourierPoissonBC = 2,
    SurfDistribMatterBC = 3,
    OneSurfDistribTestableMatterBC = 4,
    CapacitiveVoltageGenerator = 5,
};

enum class NativeVolumeFieldFamily
{
    UniformBField = 0,
    SolenoidBField = 1,
    DipolarBField = 2,
    MultipleEField = 3,
    EMField = 4,
};

enum class NativeVolumeDistributionFamily
{
    PICVolDistrib = 0,
    PICVolDistribNoAcc = 1,
    PICVolDistribUpdatable = 2,
    SmartPICVolDistrib = 3,
    CompositeVolDistrib = 4,
    BackTrackingVolDistrib = 5,
    BacktrackingPICCompositeVolDistrib = 6,
    BacktrackingBoltzmannCompositeVolDistrib = 7,
    VolDistrib = 8,
    NonPICVolDistrib = 9,
    AnalyticVolDistrib = 10,
    MultipleVolDistrib = 11,
    LocalMaxwellVolDistrib = 12,
    LocalMaxwellBoltzmannVolDistrib = 13,
    GlobalMaxwellBoltzmannVolDistrib = 14,
    PICBoltzmannVolDistrib = 15,
    SteadyMaxwellBoltzmannVolDistrib = 16,
    UnlimitedGlobalMaxwellBoltzmannVolDistrib = 17,
    SurfaceLimitedGlobalMaxwellBoltzmannVolDistrib = 18,
    TrunckatedGlobalMaxwellBoltzmannVolDistrib = 19,
    ImplicitableVolDistrib = 20,
    Updatable = 21,
};

enum class NativeVolumeInteractionFamily
{
    MCCInteractor = 0,
    CEXInteractor = 1,
    PhotoIonization = 2,
    TrajectoryInteractionFromField = 3,
    SpinningSpacecraftTrajectory = 4,
};

struct ExternalFieldSolveNodeResult
{
    std::string node_id;
    double reference_potential_v = 0.0;
    double normal_field_v_per_m = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double capacitance_scale = 1.0;
};

struct ExternalVolumeSolveCellResult
{
    std::string cell_id;
    std::string boundary_group_id;
    double potential_v = 0.0;
    double reference_potential_v = 0.0;
    double normal_field_v_per_m = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double capacitance_scale = 1.0;
    double coupling_gain = 0.0;
    double projection_weight_scale = 1.0;
    double sheath_length_scale = 1.0;
};

struct NativeVolumeMeshSummary
{
    std::size_t node_count = 0;
    std::size_t element_count = 0;
    double total_volume_m3 = 0.0;
    double bbox_span_x_m = 0.0;
    double bbox_span_y_m = 0.0;
    double bbox_span_z_m = 0.0;
    bool valid = false;
};

struct NativeVolumeFieldSummary
{
    double min_potential_v = 0.0;
    double max_potential_v = 0.0;
    double max_field_v_per_m = 0.0;
    double total_energy_j = 0.0;
    std::size_t boundary_condition_count = 0;
    std::size_t dirichlet_boundary_count = 0;
    std::size_t boundary_condition_family_count = 0;
    std::string boundary_condition_family_signature;
    std::size_t field_family_count = 0;
    std::string field_family_signature;
    bool solved = false;
};

struct NativeVolumeDistributionSummary
{
    std::size_t active_particle_count = 0;
    std::size_t electron_particle_count = 0;
    std::size_t ion_particle_count = 0;
    std::size_t distribution_channel_count = 0;
    std::size_t family_count = 0;
    std::string family_signature;
    double average_electron_density_m3 = 0.0;
    double average_ion_density_m3 = 0.0;
    double average_net_charge_density_c_per_m3 = 0.0;
};

struct NativeVolumeInteractionSummary
{
    std::size_t interactor_count = 0;
    std::size_t active_particle_count = 0;
    std::size_t executed_interaction_count = 0;
    std::size_t family_count = 0;
    std::string family_signature;
    bool executed = false;
};

struct NativeVolumeParitySnapshot
{
    std::string schema_version = "scdat.native_volume_parity.v1";
    NativeVolumeParityRoute resolved_route = NativeVolumeParityRoute::NativeMinimal;
    NativeVolumeMeshSummary mesh{};
    NativeVolumeFieldSummary field{};
    NativeVolumeDistributionSummary distribution{};
    NativeVolumeInteractionSummary interaction{};
};

const char* nativeVolumeParityRouteName(NativeVolumeParityRoute route);
const char* nativeVolumeBoundaryConditionFamilyName(NativeVolumeBoundaryConditionFamily family);
const char* nativeVolumeFieldFamilyName(NativeVolumeFieldFamily family);
const char* nativeVolumeDistributionFamilyName(NativeVolumeDistributionFamily family);
const char* nativeVolumeInteractionFamilyName(NativeVolumeInteractionFamily family);

} // namespace FieldSolver
} // namespace SCDAT
