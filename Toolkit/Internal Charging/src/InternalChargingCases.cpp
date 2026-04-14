#include "InternalChargingCases.h"

#include <array>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{
namespace
{

InternalChargingConfiguration makeGeoElectronBeltConfig()
{
    InternalChargingConfiguration config;
    config.material_stack_model = InternalMaterialStackModelKind::SpisLayeredStack;
    config.geometry_model = InternalGeometryModelKind::ShieldedLayerStack1D;
    config.primary_source_model = InternalPrimarySourceModelKind::PresetMonoEnergeticFlux;
    config.physics_process_list = InternalPhysicsProcessListKind::Geant4ShieldingLike;
    config.energy_deposition_model = InternalEnergyDepositionModelKind::Geant4StepRecorderLike;
    config.charge_response_model =
        InternalChargeResponseModelKind::RadiationInducedConductivityRelaxation;
    config.layers = 24;
    config.thickness_m = 4.0e-3;
    config.area_m2 = 5.0e-2;
    config.incident_current_density_a_per_m2 = 1.5e-8;
    config.incident_energy_ev = 1.0e5;
    config.material_name = "kapton";
    return config;
}

InternalChargingConfiguration makeMeoHarnessConfig()
{
    InternalChargingConfiguration config;
    config.material_stack_model = InternalMaterialStackModelKind::SpisHarnessBundle;
    config.geometry_model = InternalGeometryModelKind::LayerStack1D;
    config.primary_source_model = InternalPrimarySourceModelKind::PresetMonoEnergeticFlux;
    config.physics_process_list = InternalPhysicsProcessListKind::Geant4EmStandardLike;
    config.energy_deposition_model = InternalEnergyDepositionModelKind::ContinuousSlabDeposition;
    config.charge_response_model = InternalChargeResponseModelKind::SpisLayeredDielectric;
    config.layers = 20;
    config.thickness_m = 1.5e-3;
    config.area_m2 = 7.0e-3;
    config.incident_current_density_a_per_m2 = 3.0e-9;
    config.incident_energy_ev = 3.0e4;
    config.material_name = "ptfe";
    return config;
}

InternalChargingConfiguration makeBacksheetConfig()
{
    InternalChargingConfiguration config;
    config.material_stack_model = InternalMaterialStackModelKind::SpisBacksheetStack;
    config.geometry_model = InternalGeometryModelKind::LayerStack1D;
    config.primary_source_model = InternalPrimarySourceModelKind::PresetMonoEnergeticFlux;
    config.physics_process_list = InternalPhysicsProcessListKind::Geant4ShieldingLike;
    config.energy_deposition_model = InternalEnergyDepositionModelKind::Geant4StepRecorderLike;
    config.charge_response_model =
        InternalChargeResponseModelKind::RadiationInducedConductivityRelaxation;
    config.layers = 18;
    config.thickness_m = 1.0e-3;
    config.area_m2 = 2.0e-2;
    config.incident_current_density_a_per_m2 = 7.0e-9;
    config.incident_energy_ev = 2.0e4;
    config.material_name = "kapton_hn";
    return config;
}

std::array<InternalChargingScenarioPreset, 3> buildPresets()
{
    return {InternalChargingScenarioPreset{"geo_electron_belt",
                                          "GEO outer-belt dielectric charging case.",
                                          makeGeoElectronBeltConfig(), 5.0e-4, 20,
                                          "results/internal_geo_electron_belt.csv"},
            InternalChargingScenarioPreset{"meo_dielectric_harness",
                                          "MEO harness dielectric volume charging case.",
                                          makeMeoHarnessConfig(), 1.0e-4, 20,
                                          "results/internal_meo_dielectric_harness.csv"},
            InternalChargingScenarioPreset{"solar_array_backsheet",
                                          "Solar-array backsheet charging accumulation case.",
                                          makeBacksheetConfig(), 2.0e-4, 24,
                                          "results/internal_solar_array_backsheet.csv"}};
}

const auto& presets()
{
    static const auto kPresets = buildPresets();
    return kPresets;
}

} // namespace

std::vector<std::string> listInternalChargingScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(presets().size());
    for (const auto& preset : presets())
    {
        names.push_back(preset.name);
    }
    return names;
}

bool tryGetInternalChargingScenarioPreset(const std::string& name,
                                          InternalChargingScenarioPreset& preset)
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

InternalChargingScenarioPreset makeDefaultInternalChargingScenarioPreset()
{
    return presets().front();
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
