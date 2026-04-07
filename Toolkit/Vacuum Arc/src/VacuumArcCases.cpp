#include "VacuumArcCases.h"

#include <array>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{
namespace
{

DischargeConfiguration makeMicrogapConfig()
{
    DischargeConfiguration config;
    config.gap_distance_m = 5.0e-4;
    config.applied_field_v_per_m = 4.0e7;
    config.surface_potential_v = 150.0;
    config.surface_charge_density_c_per_m2 = 3.0e-5;
    config.surface_capacitance_f = 3.0e-10;
    config.surface_leakage_conductance_s = 2.5e-8;
    config.surface_reference_potential_v = 0.0;
    config.max_surface_potential_step_v = 10.0;
    config.alignment_mode = ArcPicAlignmentMode::ArcPicAligned;
    config.cathode_temperature_k = 520.0;
    config.anode_temperature_k = 340.0;
    config.channel_radius_m = 2.0e-5;
    config.emission_strategy = ArcEmissionStrategyType::FowlerNordheimLike;
    config.emission_parameters.fn_prefactor = 2.0e-18;
    config.emission_parameters.fn_barrier_v_per_m = 3.5e9;
    config.emission_parameters.fn_temperature_reference_k = 500.0;
    config.emission_parameters.fn_temperature_gain = 8.0e-3;
    return config;
}

DischargeConfiguration makeTripleJunctionConfig()
{
    DischargeConfiguration config;
    config.gap_distance_m = 2.0e-4;
    config.applied_field_v_per_m = 6.0e7;
    config.surface_potential_v = 220.0;
    config.surface_charge_density_c_per_m2 = 8.0e-5;
    config.surface_capacitance_f = 2.0e-10;
    config.surface_leakage_conductance_s = 1.0e-8;
    config.surface_reference_potential_v = 0.0;
    config.max_surface_potential_step_v = 8.0;
    config.alignment_mode = ArcPicAlignmentMode::ArcPicAligned;
    config.cathode_temperature_k = 650.0;
    config.anode_temperature_k = 360.0;
    config.channel_radius_m = 1.2e-5;
    config.emission_strategy = ArcEmissionStrategyType::LinearThreshold;
    return config;
}

DischargeConfiguration makeRestrikeConfig()
{
    DischargeConfiguration config;
    config.gap_distance_m = 7.0e-4;
    config.applied_field_v_per_m = 3.5e7;
    config.surface_potential_v = 80.0;
    config.surface_charge_density_c_per_m2 = 2.5e-5;
    config.surface_capacitance_f = 4.0e-10;
    config.surface_leakage_conductance_s = 3.0e-8;
    config.surface_reference_potential_v = 0.0;
    config.max_surface_potential_step_v = 12.0;
    config.alignment_mode = ArcPicAlignmentMode::ArcPicAligned;
    config.cathode_temperature_k = 560.0;
    config.anode_temperature_k = 345.0;
    config.channel_radius_m = 2.5e-5;
    config.emission_strategy = ArcEmissionStrategyType::FowlerNordheimLike;
    config.emission_parameters.fn_prefactor = 1.0e-18;
    config.emission_parameters.fn_barrier_v_per_m = 4.2e9;
    config.emission_parameters.fn_temperature_reference_k = 520.0;
    config.emission_parameters.fn_temperature_gain = 6.0e-3;
    return config;
}

std::array<VacuumArcScenarioPreset, 3> buildPresets()
{
    return {VacuumArcScenarioPreset{"microgap_flashover",
                                    "Compact micro-gap vacuum flashover development case.",
                                    makeMicrogapConfig(), 1.0e-9, 12,
                                    "results/arc_microgap_flashover.csv"},
            VacuumArcScenarioPreset{"triple_junction_flashover",
                                    "Triple-junction initiated vacuum surface flashover case.",
                                    makeTripleJunctionConfig(), 5.0e-10, 14,
                                    "results/arc_triple_junction_flashover.csv"},
            VacuumArcScenarioPreset{"restrike_recovery",
                                    "Post-breakdown restrike and channel recovery case.",
                                    makeRestrikeConfig(), 1.0e-9, 10,
                                    "results/arc_restrike_recovery.csv"}};
}

const auto& presets()
{
    static const auto kPresets = buildPresets();
    return kPresets;
}

} // namespace

std::vector<std::string> listVacuumArcScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(presets().size());
    for (const auto& preset : presets())
    {
        names.push_back(preset.name);
    }
    return names;
}

bool tryGetVacuumArcScenarioPreset(const std::string& name, VacuumArcScenarioPreset& preset)
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

VacuumArcScenarioPreset makeDefaultVacuumArcScenarioPreset()
{
    return presets().front();
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
