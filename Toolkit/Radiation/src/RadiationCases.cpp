#include "RadiationCases.h"

#include <array>

namespace SCDAT
{
namespace Toolkit
{
namespace Radiation
{
namespace
{

RadiationConfiguration makeGeoElectronDoseConfig()
{
    RadiationConfiguration config;
    config.layers = 24;
    config.thickness_m = 4.0e-3;
    config.area_m2 = 4.0e-2;
    config.particle_flux_m2_s = 8.0e10;
    config.mean_energy_ev = 8.5e4;
    config.particle_species = ParticleSpecies::Electron;
    config.physics_list = RadiationPhysicsList::Geant4EmLivermore;
    config.material_name = "kapton";
    config.attenuation_length_fraction = 0.20;
    return config;
}

RadiationConfiguration makeMeoProtonDoseConfig()
{
    RadiationConfiguration config;
    config.layers = 20;
    config.thickness_m = 2.0e-3;
    config.area_m2 = 1.2e-2;
    config.particle_flux_m2_s = 1.8e10;
    config.mean_energy_ev = 1.6e5;
    config.particle_species = ParticleSpecies::Proton;
    config.physics_list = RadiationPhysicsList::Geant4EmStandard;
    config.material_name = "ptfe";
    config.attenuation_length_fraction = 0.28;
    return config;
}

RadiationConfiguration makeHeavyIonStormConfig()
{
    RadiationConfiguration config;
    config.layers = 30;
    config.thickness_m = 5.0e-3;
    config.area_m2 = 3.0e-2;
    config.particle_flux_m2_s = 2.5e9;
    config.mean_energy_ev = 8.0e5;
    config.particle_species = ParticleSpecies::HeavyIon;
    config.physics_list = RadiationPhysicsList::Geant4SpaceShielding;
    config.material_name = "kapton";
    config.attenuation_length_fraction = 0.35;
    return config;
}

RadiationConfiguration makeCouplingZeroFluxBaselineConfig()
{
    RadiationConfiguration config;
    config.layers = 24;
    config.thickness_m = 4.0e-3;
    config.area_m2 = 4.0e-2;
    config.particle_flux_m2_s = 0.0;
    config.mean_energy_ev = 8.5e4;
    config.particle_species = ParticleSpecies::Electron;
    config.physics_list = RadiationPhysicsList::Geant4EmLivermore;
    config.material_name = "kapton";
    config.attenuation_length_fraction = 0.20;
    return config;
}

std::array<RadiationScenarioPreset, 4> buildPresets()
{
    return {
        RadiationScenarioPreset{
            "geo_electron_belt_dose",
            "GEO electron-belt driven dose accumulation benchmark.",
            makeGeoElectronDoseConfig(),
            5.0e-4,
            20,
            "results/radiation_geo_electron_belt_dose.csv",
        },
        RadiationScenarioPreset{
            "meo_proton_harness_dose",
            "MEO proton-driven dielectric harness dose scenario.",
            makeMeoProtonDoseConfig(),
            2.0e-4,
            24,
            "results/radiation_meo_proton_harness_dose.csv",
        },
        RadiationScenarioPreset{
            "heavy_ion_storm_dose",
            "Heavy-ion storm transient dose stress case.",
            makeHeavyIonStormConfig(),
            1.0e-4,
            30,
            "results/radiation_heavy_ion_storm_dose.csv",
        },
        RadiationScenarioPreset{
            "coupling_zero_flux_baseline",
            "Zero-flux radiation preset for coupling-off baseline checks.",
            makeCouplingZeroFluxBaselineConfig(),
            5.0e-4,
            20,
            "results/radiation_coupling_zero_flux_baseline.csv",
        },
    };
}

const auto& presets()
{
    static const auto kPresets = buildPresets();
    return kPresets;
}

} // namespace

std::vector<std::string> listRadiationScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(presets().size());
    for (const auto& preset : presets())
    {
        names.push_back(preset.name);
    }
    return names;
}

bool tryGetRadiationScenarioPreset(const std::string& name, RadiationScenarioPreset& preset)
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

RadiationScenarioPreset makeDefaultRadiationScenarioPreset()
{
    return presets().front();
}

} // namespace Radiation
} // namespace Toolkit
} // namespace SCDAT
