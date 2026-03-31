#include "SurfaceChargingCases.h"

#include <array>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

Material::MaterialProperty makeKaptonSurface()
{
    Material::MaterialProperty material(2, Mesh::MaterialType::DIELECTRIC, "kapton");
    material.setPermittivity(3.4);
    material.setConductivity(1.0e-15);
    material.setWorkFunctionEv(4.7);
    material.setSecondaryElectronYield(1.2);
    material.setBreakdownFieldVPerM(2.5e8);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 300.0);
    material.setScalarProperty("secondary_emission_escape_energy_ev", 2.5);
    material.setScalarProperty("photoelectron_yield", 0.016);
    material.setScalarProperty("photoelectron_escape_energy_ev", 1.8);
    material.setScalarProperty("thermionic_escape_energy_ev", 0.2);
    material.setScalarProperty("poole_frenkel_beta", 3.5e-3);
    material.setScalarProperty("max_field_enhancement_factor", 1.0e7);
    return material;
}

Material::MaterialProperty makePtfeSurface()
{
    Material::MaterialProperty material(3, Mesh::MaterialType::DIELECTRIC, "ptfe");
    material.setPermittivity(2.1);
    material.setConductivity(1.0e-17);
    material.setWorkFunctionEv(5.75);
    material.setSecondaryElectronYield(1.4);
    material.setBreakdownFieldVPerM(6.0e7);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 500.0);
    material.setScalarProperty("photoelectron_yield", 0.015);
    material.setScalarProperty("secondary_emission_escape_energy_ev", 4.0);
    material.setScalarProperty("photoelectron_escape_energy_ev", 1.5);
    material.setScalarProperty("thermionic_escape_energy_ev", 0.15);
    material.setScalarProperty("poole_frenkel_beta", 7.5e-3);
    material.setScalarProperty("max_field_enhancement_factor", 5.0e8);
    return material;
}

SurfaceChargingConfig makeLeoDaylightConfig()
{
    SurfaceChargingConfig config;
    config.regime = SurfaceChargingRegime::LeoFlowingPlasma;
    config.surface_area_m2 = 2.5e-2;
    config.derive_capacitance_from_material = true;
    config.dielectric_thickness_m = 1.25e-4;
    config.floating = true;
    config.bulk_flow_velocity_m_per_s = 7.6e3;
    config.flow_alignment_cosine = 1.0;
    config.electron_flow_coupling = 0.015;
    config.plasma.electron_density_m3 = 8.0e11;
    config.plasma.ion_density_m3 = 8.5e11;
    config.plasma.electron_temperature_ev = 0.28;
    config.plasma.ion_temperature_ev = 0.12;
    config.plasma.ion_mass_amu = 16.0;
    config.plasma.neutral_density_m3 = 1.0e14;
    config.plasma.pressure_pa = 1.0e-5;
    config.material = makeKaptonSurface();
    config.emission.surface_temperature_k = 320.0;
    config.emission.enhancement_factor = 1.05;
    config.emission.photon_flux_m2_s = 3.0e19;
    return config;
}

SurfaceChargingConfig makeLeoWakeConfig()
{
    SurfaceChargingConfig config = makeLeoDaylightConfig();
    config.flow_alignment_cosine = -1.0;
    config.electron_flow_coupling = 0.002;
    config.plasma.electron_density_m3 = 3.2e11;
    config.plasma.ion_density_m3 = 9.0e10;
    config.electron_collection_coefficient = 0.85;
    config.ion_collection_coefficient = 0.35;
    config.material.setScalarProperty("photoelectron_yield", 0.017);
    return config;
}

SurfaceChargingConfig makeGeoEclipseConfig()
{
    SurfaceChargingConfig config;
    config.regime = SurfaceChargingRegime::GeoKineticPicLike;
    config.surface_area_m2 = 4.0e-2;
    config.derive_capacitance_from_material = true;
    config.dielectric_thickness_m = 2.5e-4;
    config.floating = true;
    config.plasma.electron_density_m3 = 8.0e6;
    config.plasma.ion_density_m3 = 4.0e6;
    config.plasma.electron_temperature_ev = 1500.0;
    config.plasma.ion_temperature_ev = 8.0;
    config.plasma.ion_mass_amu = 1.0;
    config.plasma.neutral_density_m3 = 1.0e10;
    config.plasma.pressure_pa = 1.0e-8;
    config.enable_pic_calibration = true;
    config.pic_calibration_samples = 6144;
    config.material = makePtfeSurface();
    config.radiation_conductivity_coefficient = 2.0e-10;
    config.radiation_conductivity_exponent = 1.0;
    config.emission.surface_temperature_k = 250.0;
    config.emission.enhancement_factor = 1.05;
    config.emission.photon_flux_m2_s = 0.0;
    return config;
}

SurfaceChargingConfig makeThrusterPlumeConfig()
{
    SurfaceChargingConfig config;
    config.regime = SurfaceChargingRegime::ThrusterPlume;
    config.surface_area_m2 = 3.0e-3;
    config.derive_capacitance_from_material = true;
    config.dielectric_thickness_m = 1.5e-4;
    config.floating = true;
    config.plasma.electron_density_m3 = 2.5e15;
    config.plasma.ion_density_m3 = 2.0e15;
    config.plasma.electron_temperature_ev = 8.0;
    config.plasma.ion_temperature_ev = 1.5;
    config.plasma.ion_mass_amu = 131.3;
    config.plasma.neutral_density_m3 = 1.0e18;
    config.plasma.pressure_pa = 2.0e-3;
    config.electron_collection_coefficient = 0.30;
    config.ion_collection_coefficient = 1.15;
    config.ion_directed_velocity_m_per_s = 2.2e4;
    config.material = makeKaptonSurface();
    config.emission.surface_temperature_k = 380.0;
    config.emission.enhancement_factor = 1.15;
    config.emission.photon_flux_m2_s = 8.0e17;
    return config;
}

std::array<SurfaceChargingScenarioPreset, 4> buildPresets()
{
    SurfaceChargingScenarioPreset leo_daylight;
    leo_daylight.name = "leo_daylight_kapton";
    leo_daylight.description = "Sunlit LEO ram-facing dielectric panel charging case.";
    leo_daylight.config = makeLeoDaylightConfig();
    leo_daylight.time_step_s = 2.0e-8;
    leo_daylight.steps = 160;
    leo_daylight.adaptive_time_stepping = true;
    leo_daylight.total_duration_s = 2.0e3;
    leo_daylight.minimum_time_step_s = 2.0e-8;
    leo_daylight.maximum_time_step_s = 1.2e2;
    leo_daylight.default_output_csv = "results/surface_leo_daylight_kapton.csv";

    SurfaceChargingScenarioPreset leo_wake;
    leo_wake.name = "leo_daylight_wake_kapton";
    leo_wake.description = "Sunlit LEO wake-facing dielectric panel charging case.";
    leo_wake.config = makeLeoWakeConfig();
    leo_wake.time_step_s = 2.0e-8;
    leo_wake.steps = 160;
    leo_wake.adaptive_time_stepping = true;
    leo_wake.total_duration_s = 2.0e3;
    leo_wake.minimum_time_step_s = 2.0e-8;
    leo_wake.maximum_time_step_s = 1.2e2;
    leo_wake.default_output_csv = "results/surface_leo_daylight_wake_kapton.csv";

    SurfaceChargingScenarioPreset geo_eclipse;
    geo_eclipse.name = "geo_eclipse_dielectric";
    geo_eclipse.description = "GEO eclipse dielectric case with short-window PIC calibration and long-horizon kinetic charging.";
    geo_eclipse.config = makeGeoEclipseConfig();
    geo_eclipse.time_step_s = 2.0e-1;
    geo_eclipse.steps = 200;
    geo_eclipse.adaptive_time_stepping = true;
    geo_eclipse.total_duration_s = 2.0e3;
    geo_eclipse.minimum_time_step_s = 2.0e-1;
    geo_eclipse.maximum_time_step_s = 1.2e2;
    geo_eclipse.default_output_csv = "results/surface_geo_eclipse_dielectric.csv";

    SurfaceChargingScenarioPreset thruster_plume;
    thruster_plume.name = "thruster_plume_dielectric";
    thruster_plume.description = "Near-thruster plume driven surface current balance case.";
    thruster_plume.config = makeThrusterPlumeConfig();
    thruster_plume.time_step_s = 1.0e-9;
    thruster_plume.steps = 80;
    thruster_plume.default_output_csv = "results/surface_thruster_plume_dielectric.csv";

    return {leo_daylight, leo_wake, geo_eclipse, thruster_plume};
}

const auto& presets()
{
    static const auto kPresets = buildPresets();
    return kPresets;
}

} // namespace

std::vector<std::string> listSurfaceChargingScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(presets().size());
    for (const auto& preset : presets())
    {
        names.push_back(preset.name);
    }
    return names;
}

bool tryGetSurfaceChargingScenarioPreset(const std::string& name,
                                         SurfaceChargingScenarioPreset& preset)
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

SurfaceChargingScenarioPreset makeDefaultSurfaceChargingScenarioPreset()
{
    return presets().front();
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
