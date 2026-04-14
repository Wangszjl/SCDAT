#include "SurfaceTransitionEngine.h"

#include "../../Tools/Solver/include/SurfaceSolverFacade.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

constexpr double kTwoPi = 6.28318530717958647692;

struct SurfaceTransitionEventEvaluation
{
    bool local_time_transition_active = false;
    bool spinning_spacecraft_active = false;
    bool sun_flux_updater_active = false;
    bool conductivity_evolution_active = false;
    bool source_flux_updater_active = false;
    bool simulation_param_updater_active = false;
    bool sheath_or_presheath_poisson_bc_updater_active = false;
    bool rccabs_sc_updater_active = false;
    bool vcross_bfield_updater_active = false;
    bool basic_eclipse_exit_active = false;
    bool transient_artificial_sources_active = false;
    bool langmuir_probe_transition_active = false;
    double local_time_hour = 12.0;
    double spinning_phase_rad = 0.0;
    double sun_flux_scale = 1.0;
};

double wrapPhaseRad(double phase_rad)
{
    const double wrapped = std::fmod(phase_rad, kTwoPi);
    return wrapped < 0.0 ? wrapped + kTwoPi : wrapped;
}

double wrapHour24(double hour)
{
    const double wrapped = std::fmod(hour, 24.0);
    return wrapped < 0.0 ? wrapped + 24.0 : wrapped;
}

bool hourInsideWindow(double hour, double start_hour, double end_hour)
{
    if (start_hour <= end_hour)
    {
        return hour >= start_hour && hour <= end_hour;
    }
    return hour >= start_hour || hour <= end_hour;
}

bool transitionEventEnabled(const SurfaceChargingConfig& config, const char* property_key)
{
    return config.material.getScalarProperty(property_key, 0.0) > 0.5;
}

SurfaceTransitionEventEvaluation evaluateTransitionEvents(
    const SurfaceAdvanceTransitionInput& input)
{
    SurfaceTransitionEventEvaluation result;

    const bool local_time_enabled =
        input.config.material.getScalarProperty("transition_local_time_enabled", 0.0) > 0.5;
    const double local_time_base_hour =
        input.config.material.getScalarProperty("transition_local_time_base_hour", 12.0);
    const double daylight_start_hour = wrapHour24(
        input.config.material.getScalarProperty("transition_local_time_daylight_start_hour", 6.0));
    const double daylight_end_hour = wrapHour24(
        input.config.material.getScalarProperty("transition_local_time_daylight_end_hour", 18.0));
    result.local_time_hour = wrapHour24(local_time_base_hour + input.status.time_s / 3600.0);
    result.local_time_transition_active =
        local_time_enabled &&
        hourInsideWindow(result.local_time_hour, daylight_start_hour, daylight_end_hour);

    const bool spinning_enabled =
        input.config.material.getScalarProperty("transition_spinning_enabled", 0.0) > 0.5;
    const double spin_period_s = std::max(
        0.0, input.config.material.getScalarProperty("transition_spin_period_s", 0.0));
    const double spin_phase_offset_rad =
        input.config.material.getScalarProperty("transition_spin_phase_offset_rad", 0.0);
    if (spinning_enabled && spin_period_s > 1.0e-9)
    {
        result.spinning_spacecraft_active = true;
        result.spinning_phase_rad = wrapPhaseRad(
            kTwoPi * (input.status.time_s / spin_period_s) + spin_phase_offset_rad);
    }

    const bool sun_flux_enabled =
        input.config.material.getScalarProperty("transition_sun_flux_enabled", 0.0) > 0.5;

    if (sun_flux_enabled)
    {
        result.sun_flux_updater_active = true;
        const double base_scale = std::max(
            0.0, input.config.material.getScalarProperty("transition_sun_flux_base_scale", 1.0));
        const double night_scale = std::clamp(
            input.config.material.getScalarProperty("transition_sun_flux_night_scale", 0.0), 0.0,
            1.0);
        const double spin_amplitude = std::clamp(
            input.config.material.getScalarProperty("transition_sun_flux_spin_amplitude", 0.0), 0.0,
            1.0);

        const double local_time_scale =
            local_time_enabled ? (result.local_time_transition_active ? 1.0 : night_scale) : 1.0;
        const double spin_scale =
            result.spinning_spacecraft_active
                ? std::max(0.0, 1.0 + spin_amplitude * std::cos(result.spinning_phase_rad))
                : 1.0;
        result.sun_flux_scale = std::clamp(base_scale * local_time_scale * spin_scale, 0.0, 10.0);
    }

    result.conductivity_evolution_active = transitionEventEnabled(
        input.config, "transition_conductivity_evolution_enabled");
    result.source_flux_updater_active = transitionEventEnabled(
        input.config, "transition_source_flux_updater_enabled");
    result.simulation_param_updater_active = transitionEventEnabled(
        input.config, "transition_simulation_param_updater_enabled");
    result.sheath_or_presheath_poisson_bc_updater_active = transitionEventEnabled(
        input.config, "transition_sheath_or_presheath_poisson_bc_updater_enabled");
    result.rccabs_sc_updater_active = transitionEventEnabled(
        input.config, "transition_rccabs_sc_updater_enabled");
    result.vcross_bfield_updater_active = transitionEventEnabled(
        input.config, "transition_vcross_bfield_updater_enabled");
    result.basic_eclipse_exit_active = transitionEventEnabled(
        input.config, "transition_basic_eclipse_exit_enabled");
    result.transient_artificial_sources_active = transitionEventEnabled(
        input.config, "transition_transient_artificial_sources_enabled");
    result.langmuir_probe_transition_active = transitionEventEnabled(
        input.config, "transition_langmuir_probe_transition_enabled");

    return result;
}

bool sharedSurfacePicRuntimeEnabled(const SurfaceChargingConfig& config)
{
    return config.surface_pic_runtime_kind == SurfacePicRuntimeKind::GraphCoupledSharedSurface;
}

std::size_t resolveSharedGlobalCoupledIterationLimit(const SurfaceChargingConfig& config,
                                                     std::size_t patch_count)
{
    if (!sharedSurfacePicRuntimeEnabled(config) || patch_count < 2)
    {
        return 0;
    }

    const double configured_iterations = std::clamp(
        config.material.getScalarProperty("shared_surface_global_coupled_iterations", 0.0), 0.0,
        256.0);
    const std::size_t configured_limit =
        static_cast<std::size_t>(std::llround(configured_iterations));
    Solver::GlobalCoupledControlInput input;
    input.configured_iteration_limit = configured_limit;
    input.configured_tolerance_v = std::max(
        1.0e-9,
        config.material.getScalarProperty("shared_surface_global_coupled_tolerance_v", 1.0e-6));
    input.configured_relaxation = 1.0;
    input.solver_max_iterations = config.solver_config.max_iterations;
    input.solver_residual_tolerance = config.solver_config.residual_tolerance;
    input.solver_relaxation_factor = config.solver_config.relaxation_factor;
    input.coupling_mode = config.solver_config.coupling_mode;
    input.convergence_policy = config.solver_config.convergence_policy;
    return Solver::resolveGlobalCoupledControl(input).iteration_limit;
}

} // namespace

SurfaceAdvanceTransition evaluateSurfaceAdvanceTransition(
    const SurfaceAdvanceTransitionInput& input)
{
    SurfaceAdvanceTransition transition;
    transition.runtime_route = input.config.runtime_route;
    transition.surface_pic_strategy = input.config.surface_pic_strategy;
    transition.legacy_benchmark_execution_mode = input.config.legacy_benchmark_execution_mode;
    const bool legacy_route =
        input.config.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark;
    const bool replay_mode =
        input.config.legacy_benchmark_execution_mode ==
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    transition.use_legacy_benchmark_replay =
        input.has_legacy_benchmark_replay && legacy_route && replay_mode;

    transition.use_reference_circuit =
        input.config.use_reference_current_balance &&
        (input.config.regime == SurfaceChargingRegime::LeoFlowingPlasma ||
         input.config.regime == SurfaceChargingRegime::GeoKineticPicLike) &&
        input.config.enable_body_patch_circuit;

    transition.shared_surface_pic_runtime = sharedSurfacePicRuntimeEnabled(input.config);
    transition.shared_global_coupled_iteration_limit =
        resolveSharedGlobalCoupledIterationLimit(input.config, input.patch_count);
    transition.shared_global_coupled_solve =
        transition.shared_surface_pic_runtime && input.has_circuit_model && input.patch_count >= 2 &&
        transition.shared_global_coupled_iteration_limit > 1;

    transition.shared_live_pic_coupled_refresh =
        transition.shared_global_coupled_solve &&
        (input.config.enable_live_pic_window || input.config.enable_pic_calibration) &&
        input.config.material.getScalarProperty("shared_surface_live_pic_coupled_refresh", 0.0) >
            0.5;

    transition.shared_particle_transport_coupling =
        transition.shared_global_coupled_solve && input.has_circuit_model && input.patch_count >= 2 &&
        input.config.material.getScalarProperty("shared_surface_particle_transport_weight", 0.0) >
            0.0;

    const auto event_evaluation = evaluateTransitionEvents(input);
    transition.local_time_transition_active = event_evaluation.local_time_transition_active;
    transition.spinning_spacecraft_active = event_evaluation.spinning_spacecraft_active;
    transition.sun_flux_updater_active = event_evaluation.sun_flux_updater_active;
    transition.local_time_hour = event_evaluation.local_time_hour;
    transition.spinning_phase_rad = event_evaluation.spinning_phase_rad;
    transition.sun_flux_scale = event_evaluation.sun_flux_scale;
    transition.extended_events.conductivity_evolution_active =
        event_evaluation.conductivity_evolution_active;
    transition.extended_events.source_flux_updater_active =
        event_evaluation.source_flux_updater_active;
    transition.extended_events.simulation_param_updater_active =
        event_evaluation.simulation_param_updater_active;
    transition.extended_events.sheath_or_presheath_poisson_bc_updater_active =
        event_evaluation.sheath_or_presheath_poisson_bc_updater_active;
    transition.extended_events.rccabs_sc_updater_active =
        event_evaluation.rccabs_sc_updater_active;
    transition.extended_events.vcross_bfield_updater_active =
        event_evaluation.vcross_bfield_updater_active;
    transition.extended_events.basic_eclipse_exit_active =
        event_evaluation.basic_eclipse_exit_active;
    transition.extended_events.transient_artificial_sources_active =
        event_evaluation.transient_artificial_sources_active;
    transition.extended_events.langmuir_probe_transition_active =
        event_evaluation.langmuir_probe_transition_active;

    transition.pic_recalibration_requested =
        input.config.enable_pic_calibration &&
        (input.status.steps_completed == 0 ||
         (input.config.pic_recalibration_interval_steps > 0 &&
          input.status.steps_completed % input.config.pic_recalibration_interval_steps == 0) ||
         std::abs(input.status.equilibrium_error) > input.config.pic_recalibration_trigger_v);

    if (transition.pic_recalibration_requested && transition.sun_flux_updater_active)
    {
        const double minimum_flux_scale = std::clamp(
            input.config.material.getScalarProperty(
                "transition_sun_flux_pic_recalibration_min_scale", 0.0),
            0.0, 1.0);
        if (transition.sun_flux_scale < minimum_flux_scale)
        {
            transition.pic_recalibration_requested = false;
        }
    }

    const auto solver_policy_flags = Solver::resolveSolverPolicyFlags(
        input.config.solver_config.coupling_mode, input.config.solver_config.convergence_policy);
    transition.fixed_iteration_policy = solver_policy_flags.fixed_iteration_policy_requested;
    transition.residual_guarded_policy = solver_policy_flags.residual_guarded_requested;

    return transition;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
