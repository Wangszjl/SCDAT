#include "SurfaceTransitionEngine.h"

#include "../../Tools/Solver/include/SurfaceSolverFacade.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

constexpr double kTwoPi = 6.28318530717958647692;
constexpr std::array<const char*, 3> kTopTransitionInterfaceFamilies = {
    "TransitionInterface",
    "Transition",
    "SimulationParamUpdater",
};
constexpr std::array<const char*, 15> kTopTransitionSupportedFamilies = {
    "LocalTimeTransition",
    "SpinningSpacecraft",
    "SunFluxUpdater",
    "ConductivityEvolution",
    "SourceFluxUpdater",
    "SimulationParamUpdater",
    "SheathOrPresheathPoissonBCUpdater",
    "RCCabsSCUpdater",
    "VcrossBfieldUpdater",
    "BasicEclipseExit",
    "TransientArtificialSources",
    "LangmuirProbeTransition",
    "TransitionObserver",
    "Finalization",
    "SunFluxIntensityUpdater",
};
constexpr std::array<const char*, 3> kTopTransitionLifecycleFamilies = {
    "TransitionObserver",
    "Finalization",
    "SunFluxIntensityUpdater",
};

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
    double sun_flux_daylight_factor = 1.0;
    double sun_flux_spin_factor = 1.0;
    double sun_flux_base_scale = 1.0;
};

template <typename Range>
std::string joinFamilyRange(const Range& families)
{
    std::string signature;
    bool first = true;
    for (const auto& family : families)
    {
        if (!first)
        {
            signature += "+";
        }
        first = false;
        signature += family;
    }
    return signature;
}

void appendUniqueFamily(std::vector<std::string>& families, const char* family)
{
    if (std::find(families.begin(), families.end(), family) == families.end())
    {
        families.emplace_back(family);
    }
}

std::size_t countActiveTransitionEvents(const SurfaceTransitionEventEvaluation& evaluation)
{
    const bool flags[] = {
        evaluation.local_time_transition_active,
        evaluation.spinning_spacecraft_active,
        evaluation.sun_flux_updater_active,
        evaluation.conductivity_evolution_active,
        evaluation.source_flux_updater_active,
        evaluation.simulation_param_updater_active,
        evaluation.sheath_or_presheath_poisson_bc_updater_active,
        evaluation.rccabs_sc_updater_active,
        evaluation.vcross_bfield_updater_active,
        evaluation.basic_eclipse_exit_active,
        evaluation.transient_artificial_sources_active,
        evaluation.langmuir_probe_transition_active,
    };
    return static_cast<std::size_t>(std::count(std::begin(flags), std::end(flags), true));
}

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

SurfaceTransitionObserverState evaluateTransitionObserver(
    const SurfaceAdvanceTransitionInput& input,
    const SurfaceTransitionEventEvaluation& evaluation)
{
    SurfaceTransitionObserverState observer;
    observer.observed_transition_count = countActiveTransitionEvents(evaluation);
    observer.active =
        input.config.material.getScalarProperty("transition_observer_enabled", 1.0) > 0.5 &&
        observer.observed_transition_count > 0;
    observer.checkpoint_dt_s = std::max(
        0.0, input.config.material.getScalarProperty("transition_observer_checkpoint_dt_s", 0.0));
    if (observer.active && observer.checkpoint_dt_s > 1.0e-12)
    {
        const double current_bucket =
            std::floor(input.status.time_s / observer.checkpoint_dt_s) + 1.0;
        observer.next_checkpoint_time_s = current_bucket * observer.checkpoint_dt_s;
    }
    else
    {
        observer.next_checkpoint_time_s = input.status.time_s;
    }
    return observer;
}

SurfaceTransitionFinalizationState evaluateTransitionFinalization(
    const SurfaceAdvanceTransitionInput& input)
{
    SurfaceTransitionFinalizationState finalization;
    finalization.trigger_time_s = std::max(
        0.0, input.config.material.getScalarProperty("transition_finalization_trigger_time_s", 0.0));
    finalization.duration_s = std::max(
        0.0, input.config.material.getScalarProperty("transition_finalization_duration_s", 0.0));
    finalization.validity_renormalization = std::max(
        0.0,
        input.config.material.getScalarProperty(
            "transition_finalization_validity_renormalization", 1.0));
    finalization.active =
        input.config.material.getScalarProperty("transition_finalization_enabled", 0.0) > 0.5 &&
        input.status.time_s >= finalization.trigger_time_s;
    return finalization;
}

SurfaceSunFluxIntensityState evaluateSunFluxIntensityState(
    const SurfaceTransitionEventEvaluation& evaluation)
{
    SurfaceSunFluxIntensityState intensity;
    intensity.active = evaluation.sun_flux_updater_active;
    intensity.intensity_scale = evaluation.sun_flux_scale;
    intensity.daylight_factor = evaluation.sun_flux_daylight_factor;
    intensity.spin_factor = evaluation.sun_flux_spin_factor;
    const double base_scale = std::max(1.0e-12, evaluation.sun_flux_base_scale);
    intensity.normalized_scale =
        evaluation.sun_flux_updater_active ? evaluation.sun_flux_scale / base_scale : 1.0;
    return intensity;
}

SurfaceTransitionObjectLayerState buildTransitionObjectLayer(
    const SurfaceTransitionEventEvaluation& evaluation,
    const SurfaceTransitionObserverState& observer,
    const SurfaceTransitionFinalizationState& finalization,
    const SurfaceSunFluxIntensityState& sun_flux_intensity)
{
    SurfaceTransitionObjectLayerState object_layer;
    object_layer.interface_layer_family_count = kTopTransitionInterfaceFamilies.size();
    object_layer.interface_layer_family_signature =
        joinFamilyRange(kTopTransitionInterfaceFamilies);
    object_layer.supported_family_count = kTopTransitionSupportedFamilies.size();
    object_layer.supported_family_signature =
        joinFamilyRange(kTopTransitionSupportedFamilies);
    object_layer.lifecycle_family_count = kTopTransitionLifecycleFamilies.size();
    object_layer.lifecycle_family_signature =
        joinFamilyRange(kTopTransitionLifecycleFamilies);

    std::vector<std::string> active_families;
    if (evaluation.local_time_transition_active)
    {
        appendUniqueFamily(active_families, "LocalTimeTransition");
    }
    if (evaluation.spinning_spacecraft_active)
    {
        appendUniqueFamily(active_families, "SpinningSpacecraft");
    }
    if (evaluation.sun_flux_updater_active)
    {
        appendUniqueFamily(active_families, "SunFluxUpdater");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (evaluation.conductivity_evolution_active)
    {
        appendUniqueFamily(active_families, "ConductivityEvolution");
    }
    if (evaluation.source_flux_updater_active)
    {
        appendUniqueFamily(active_families, "SourceFluxUpdater");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (evaluation.simulation_param_updater_active)
    {
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (evaluation.sheath_or_presheath_poisson_bc_updater_active)
    {
        appendUniqueFamily(active_families, "SheathOrPresheathPoissonBCUpdater");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (evaluation.rccabs_sc_updater_active)
    {
        appendUniqueFamily(active_families, "RCCabsSCUpdater");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (evaluation.vcross_bfield_updater_active)
    {
        appendUniqueFamily(active_families, "VcrossBfieldUpdater");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (evaluation.basic_eclipse_exit_active)
    {
        appendUniqueFamily(active_families, "BasicEclipseExit");
    }
    if (evaluation.transient_artificial_sources_active)
    {
        appendUniqueFamily(active_families, "TransientArtificialSources");
    }
    if (evaluation.langmuir_probe_transition_active)
    {
        appendUniqueFamily(active_families, "LangmuirProbeTransition");
    }
    if (observer.active)
    {
        appendUniqueFamily(active_families, "TransitionObserver");
    }
    if (finalization.active)
    {
        appendUniqueFamily(active_families, "Finalization");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }
    if (sun_flux_intensity.active)
    {
        appendUniqueFamily(active_families, "SunFluxIntensityUpdater");
        appendUniqueFamily(active_families, "SimulationParamUpdater");
    }

    std::vector<std::string> active_lifecycle_families;
    if (observer.active)
    {
        appendUniqueFamily(active_lifecycle_families, "TransitionObserver");
    }
    if (finalization.active)
    {
        appendUniqueFamily(active_lifecycle_families, "Finalization");
    }
    if (sun_flux_intensity.active)
    {
        appendUniqueFamily(active_lifecycle_families, "SunFluxIntensityUpdater");
    }

    const bool transition_layer_active =
        evaluation.local_time_transition_active || evaluation.spinning_spacecraft_active ||
        evaluation.conductivity_evolution_active || evaluation.basic_eclipse_exit_active ||
        evaluation.transient_artificial_sources_active ||
        evaluation.langmuir_probe_transition_active;
    const bool simulation_param_updater_layer_active =
        evaluation.sun_flux_updater_active || evaluation.source_flux_updater_active ||
        evaluation.simulation_param_updater_active ||
        evaluation.sheath_or_presheath_poisson_bc_updater_active ||
        evaluation.rccabs_sc_updater_active || evaluation.vcross_bfield_updater_active ||
        finalization.active || sun_flux_intensity.active;

    std::vector<std::string> active_interface_layer_families;
    if (!active_families.empty())
    {
        appendUniqueFamily(active_interface_layer_families, "TransitionInterface");
    }
    if (transition_layer_active)
    {
        appendUniqueFamily(active_interface_layer_families, "Transition");
    }
    if (simulation_param_updater_layer_active)
    {
        appendUniqueFamily(active_interface_layer_families, "SimulationParamUpdater");
    }

    object_layer.active_interface_layer_family_count = active_interface_layer_families.size();
    object_layer.active_interface_layer_family_signature =
        joinFamilyRange(active_interface_layer_families);
    object_layer.active_family_count = active_families.size();
    object_layer.active_family_signature = joinFamilyRange(active_families);
    object_layer.active_lifecycle_family_count = active_lifecycle_families.size();
    object_layer.active_lifecycle_family_signature =
        joinFamilyRange(active_lifecycle_families);
    return object_layer;
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
        result.sun_flux_base_scale = base_scale;
        result.sun_flux_daylight_factor = local_time_scale;
        result.sun_flux_spin_factor = spin_scale;
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
    transition.observer = evaluateTransitionObserver(input, event_evaluation);
    transition.finalization = evaluateTransitionFinalization(input);
    transition.sun_flux_intensity = evaluateSunFluxIntensityState(event_evaluation);
    transition.object_layer = buildTransitionObjectLayer(event_evaluation, transition.observer,
                                                        transition.finalization,
                                                        transition.sun_flux_intensity);
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
