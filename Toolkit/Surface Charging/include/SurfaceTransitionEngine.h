#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include <cstddef>

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
    SurfaceTransitionExtendedEvents extended_events{};
};

SurfaceAdvanceTransition evaluateSurfaceAdvanceTransition(
    const SurfaceAdvanceTransitionInput& input);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
