#include "PlasmaAnalysisCases.h"

#include <array>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{

FluidAlgorithmConfig makeCcpArgonConfig()
{
    FluidAlgorithmConfig config;
    config.domain_size = Geometry::Vector3D(2.0e-3, 2.0e-3, 2.0e-2);
    config.resolution = Geometry::Vector3D(8.0, 8.0, 48.0);
    config.time_step_s = 5.0e-10;
    config.dense_plasma_threshold_m3 = 3.0e16;
    config.initial_potential_v = 180.0;
    config.initial_plasma.electron_density_m3 = 3.5e16;
    config.initial_plasma.ion_density_m3 = 3.5e16;
    config.initial_plasma.electron_temperature_ev = 3.5;
    config.initial_plasma.ion_temperature_ev = 0.15;
    config.initial_plasma.neutral_density_m3 = 2.4e21;
    config.initial_plasma.electric_field_v_per_m = 1.5e4;
    config.initial_plasma.pressure_pa = 10.0;
    return config;
}

FluidAlgorithmConfig makeLeoWakeConfig()
{
    FluidAlgorithmConfig config;
    config.domain_size = Geometry::Vector3D(0.1, 0.1, 0.6);
    config.resolution = Geometry::Vector3D(6.0, 6.0, 40.0);
    config.time_step_s = 2.0e-7;
    config.dense_plasma_threshold_m3 = 1.0e12;
    config.initial_potential_v = -8.0;
    config.initial_plasma.electron_density_m3 = 2.5e11;
    config.initial_plasma.ion_density_m3 = 2.0e11;
    config.initial_plasma.electron_temperature_ev = 0.25;
    config.initial_plasma.ion_temperature_ev = 0.12;
    config.initial_plasma.neutral_density_m3 = 5.0e12;
    config.initial_plasma.electric_field_v_per_m = 35.0;
    config.initial_plasma.pressure_pa = 1.0e-5;
    return config;
}

FluidAlgorithmConfig makeThrusterPlumeConfig()
{
    FluidAlgorithmConfig config;
    config.domain_size = Geometry::Vector3D(2.0e-2, 2.0e-2, 0.12);
    config.resolution = Geometry::Vector3D(8.0, 8.0, 36.0);
    config.time_step_s = 2.5e-8;
    config.dense_plasma_threshold_m3 = 5.0e15;
    config.initial_potential_v = 28.0;
    config.initial_plasma.electron_density_m3 = 8.0e15;
    config.initial_plasma.ion_density_m3 = 7.6e15;
    config.initial_plasma.electron_temperature_ev = 10.0;
    config.initial_plasma.ion_temperature_ev = 1.2;
    config.initial_plasma.neutral_density_m3 = 8.0e18;
    config.initial_plasma.electric_field_v_per_m = 2.2e3;
    config.initial_plasma.pressure_pa = 5.0e-2;
    return config;
}

std::array<PlasmaScenarioPreset, 3> buildPresets()
{
    return {PlasmaScenarioPreset{"ccp_argon_10pa", "Low-pressure CCP argon discharge reference case.",
                                 makeCcpArgonConfig(), 20,
                                 "results/plasma_ccp_argon_10pa.csv"},
            PlasmaScenarioPreset{"leo_wake", "Low Earth orbit wake-side plasma relaxation case.",
                                 makeLeoWakeConfig(), 16, "results/plasma_leo_wake.csv"},
            PlasmaScenarioPreset{"hall_thruster_plume",
                                 "Near-field Hall thruster plume transport and sheath case.",
                                 makeThrusterPlumeConfig(), 18,
                                 "results/plasma_hall_thruster_plume.csv"}};
}

const auto& presets()
{
    static const auto kPresets = buildPresets();
    return kPresets;
}

} // namespace

std::vector<std::string> listPlasmaScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(presets().size());
    for (const auto& preset : presets())
    {
        names.push_back(preset.name);
    }
    return names;
}

bool tryGetPlasmaScenarioPreset(const std::string& name, PlasmaScenarioPreset& preset)
{
    for (const auto& candidate : presets())
    {
        if (candidate.name == name)
        {
            preset = candidate;
            return true;
        }
    }
    return false;
}

PlasmaScenarioPreset makeDefaultPlasmaScenarioPreset()
{
    return presets().front();
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
