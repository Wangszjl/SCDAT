#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include "../../Tools/Solver/include/SurfaceSolverFacade.h"

#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct SurfaceChargingScenarioPreset;

struct SurfaceRuntimePlan
{
    SurfaceChargingConfig compiled_config;
    Solver::SolverPolicyFlags solver_policy_flags;
    SurfaceRuntimeRoute runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    SurfacePicStrategy surface_pic_strategy = SurfacePicStrategy::SurfacePicCalibrated;
    SurfaceLegacyInputAdapterKind legacy_input_adapter_kind =
        SurfaceLegacyInputAdapterKind::None;
    SurfacePicRuntimeKind surface_pic_runtime_kind =
        SurfacePicRuntimeKind::LocalWindowSampler;
    SurfaceInstrumentSetKind surface_instrument_set_kind =
        SurfaceInstrumentSetKind::MetadataOnly;
    SurfaceBenchmarkSource benchmark_source = SurfaceBenchmarkSource::None;
    LegacyBenchmarkExecutionMode legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    bool uses_legacy_benchmark_route = false;
    bool executes_legacy_algorithm = false;
    bool external_field_solver_bridge_enabled = false;
    bool external_volume_solver_bridge_enabled = false;
    std::string source_name;
    std::string source_description;
};

SurfaceRuntimePlan compileSurfaceRuntimePlan(const SurfaceChargingConfig& config);
SurfaceRuntimePlan compileSurfaceRuntimePlan(const SurfaceChargingScenarioPreset& preset);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
