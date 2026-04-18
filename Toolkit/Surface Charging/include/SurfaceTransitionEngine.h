#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include <cstddef>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceAdvanceTransitionInput
{
    SurfaceChargingConfig config;
    SurfaceChargingStatus status;
    bool has_legacy_benchmark_replay = false;
    bool has_circuit_model = false;
    std::size_t patch_count = 0;
};

struct SurfaceTransitionExtendedEvents
{
    bool conductivity_evolution_active = false;
    bool source_flux_updater_active = false;
    bool simulation_param_updater_active = false;
    bool sheath_or_presheath_poisson_bc_updater_active = false;
    bool rccabs_sc_updater_active = false;
    bool vcross_bfield_updater_active = false;
    bool basic_eclipse_exit_active = false;
    bool transient_artificial_sources_active = false;
    bool langmuir_probe_transition_active = false;
};

struct SurfaceTransitionObserverState
{
    bool active = false;
    std::size_t observed_transition_count = 0;
    double checkpoint_dt_s = 0.0;
    double next_checkpoint_time_s = 0.0;
};

struct SurfaceTransitionFinalizationState
{
    bool active = false;
    double trigger_time_s = 0.0;
    double duration_s = 0.0;
    double validity_renormalization = 1.0;
};

struct SurfaceSunFluxIntensityState
{
    bool active = false;
    double intensity_scale = 1.0;
    double normalized_scale = 1.0;
    double daylight_factor = 1.0;
    double spin_factor = 1.0;
};

struct SurfaceTransitionObjectLayerState
{
    std::size_t interface_layer_family_count = 0;
    std::size_t active_interface_layer_family_count = 0;
    std::string interface_layer_family_signature;
    std::string active_interface_layer_family_signature;
    std::size_t supported_family_count = 0;
    std::size_t active_family_count = 0;
    std::string supported_family_signature;
    std::string active_family_signature;
    std::size_t lifecycle_family_count = 0;
    std::size_t active_lifecycle_family_count = 0;
    std::string lifecycle_family_signature;
    std::string active_lifecycle_family_signature;
};

struct SurfaceAdvanceTransition
{
    SurfaceRuntimeRoute runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    SurfacePicStrategy surface_pic_strategy = SurfacePicStrategy::SurfacePicCalibrated;
    LegacyBenchmarkExecutionMode legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    bool use_legacy_benchmark_replay = false;
    bool use_reference_circuit = false;
    bool shared_surface_pic_runtime = false;
    bool shared_global_coupled_solve = false;
    bool shared_live_pic_coupled_refresh = false;
    bool shared_particle_transport_coupling = false;
    bool pic_recalibration_requested = false;
    bool fixed_iteration_policy = false;
    bool residual_guarded_policy = false;
    std::size_t shared_global_coupled_iteration_limit = 0;
    bool local_time_transition_active = false;
    bool spinning_spacecraft_active = false;
    bool sun_flux_updater_active = false;
    double local_time_hour = 12.0;
    double spinning_phase_rad = 0.0;
    double sun_flux_scale = 1.0;
    SurfaceTransitionObserverState observer{};
    SurfaceTransitionFinalizationState finalization{};
    SurfaceSunFluxIntensityState sun_flux_intensity{};
    SurfaceTransitionExtendedEvents extended_events{};
    SurfaceTransitionObjectLayerState object_layer{};
};

SurfaceAdvanceTransition evaluateSurfaceAdvanceTransition(
    const SurfaceAdvanceTransitionInput& input);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
