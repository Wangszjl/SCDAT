#include "SurfaceRuntimePlan.h"

#include "LegacyBenchmarkSupport.h"
#include "SurfaceChargingCases.h"

#include "../../Tools/Particle/include/SurfaceDistributionFunction.h"
#include "../../Tools/Solver/include/SurfaceSolverFacade.h"

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

void synchronizePlasmaMomentsFromSpectra(SurfaceChargingConfig& config)
{
    if (config.has_electron_spectrum)
    {
        config.plasma.electron_density_m3 =
            Particle::resolvedSpectrumDensity(config.electron_spectrum,
                                              config.plasma.electron_density_m3);
        config.plasma.electron_temperature_ev =
            Particle::resolvedSpectrumCharacteristicEnergyEv(
                config.electron_spectrum, config.plasma.electron_temperature_ev);
    }
    if (config.has_ion_spectrum)
    {
        config.plasma.ion_density_m3 =
            Particle::resolvedSpectrumDensity(config.ion_spectrum,
                                              config.plasma.ion_density_m3);
        config.plasma.ion_temperature_ev =
            Particle::resolvedSpectrumCharacteristicEnergyEv(
                config.ion_spectrum, config.plasma.ion_temperature_ev);
        config.plasma.ion_mass_amu =
            Particle::resolvedSpectrumAverageMassAmu(config.ion_spectrum,
                                                     config.plasma.ion_mass_amu);
    }
}

SurfaceRuntimePlan buildRuntimePlan(SurfaceChargingConfig config)
{
    SurfaceRuntimePlan plan;

    plan.runtime_route = config.runtime_route;
    plan.surface_pic_strategy = config.surface_pic_strategy;
    plan.legacy_input_adapter_kind = config.legacy_input_adapter_kind;
    plan.surface_pic_runtime_kind = config.surface_pic_runtime_kind;
    plan.surface_instrument_set_kind = config.surface_instrument_set_kind;
    plan.benchmark_source = config.benchmark_source;
    plan.legacy_benchmark_execution_mode = config.legacy_benchmark_execution_mode;
    plan.uses_legacy_benchmark_route =
        config.runtime_route == SurfaceRuntimeRoute::LegacyBenchmark;
    plan.executes_legacy_algorithm =
        plan.uses_legacy_benchmark_route &&
        config.legacy_benchmark_execution_mode ==
            LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
    plan.external_field_solver_bridge_enabled = config.enable_external_field_solver_bridge;
    plan.external_volume_solver_bridge_enabled = config.enable_external_volume_solver_bridge;

    plan.solver_policy_flags = Solver::resolveSolverPolicyFlags(
        config.solver_config.coupling_mode, config.solver_config.convergence_policy);
    if (plan.solver_policy_flags.residual_guarded_requested)
    {
        config.internal_substeps = std::max<std::size_t>(config.internal_substeps, 24);
    }
    if (plan.solver_policy_flags.implicit_coupling_requested)
    {
        const std::size_t solver_iteration_limit =
            std::clamp<std::size_t>(config.solver_config.max_iterations, 2, 256);
        config.material.setScalarProperty("shared_surface_global_coupled_iterations",
                                          static_cast<double>(solver_iteration_limit));
        config.material.setScalarProperty(
            "shared_surface_global_coupled_tolerance_v",
            std::max(1.0e-9, config.solver_config.residual_tolerance));
    }

    if (plan.executes_legacy_algorithm)
    {
        config = applyLegacyBenchmarkExecutionConfig(config);
    }

    config = normalizeSurfaceChargingConfig(config);
    synchronizePlasmaMomentsFromSpectra(config);

    plan.runtime_route = config.runtime_route;
    plan.surface_pic_strategy = config.surface_pic_strategy;
    plan.legacy_input_adapter_kind = config.legacy_input_adapter_kind;
    plan.surface_pic_runtime_kind = config.surface_pic_runtime_kind;
    plan.surface_instrument_set_kind = config.surface_instrument_set_kind;
    plan.benchmark_source = config.benchmark_source;
    plan.legacy_benchmark_execution_mode = config.legacy_benchmark_execution_mode;
    plan.external_field_solver_bridge_enabled = config.enable_external_field_solver_bridge;
    plan.external_volume_solver_bridge_enabled = config.enable_external_volume_solver_bridge;
    plan.compiled_config = std::move(config);
    return plan;
}

} // namespace

SurfaceRuntimePlan compileSurfaceRuntimePlan(const SurfaceChargingConfig& config)
{
    return buildRuntimePlan(config);
}

SurfaceRuntimePlan compileSurfaceRuntimePlan(const SurfaceChargingScenarioPreset& preset)
{
    SurfaceRuntimePlan plan = buildRuntimePlan(preset.config);
    plan.source_name = preset.name;
    plan.source_description = preset.description;
    return plan;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
