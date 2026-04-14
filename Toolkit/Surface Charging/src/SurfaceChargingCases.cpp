#include "SurfaceChargingCases.h"

#include "../../Tools/Basic/include/Constants.h"

#include <array>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{
using SCDAT::Basic::Constants::PhysicsConstants;

double thermalSpeedForSourceModel(const Particle::ParticleType& particle_type, double temperature_ev)
{
    double mass = PhysicsConstants::ProtonMass;
    if (particle_type == Particle::ParticleType::ELECTRON ||
        particle_type == Particle::ParticleType::PHOTOELECTRON ||
        particle_type == Particle::ParticleType::SECONDARY_ELECTRON)
    {
        mass = PhysicsConstants::ElectronMass;
    }
    return std::sqrt(2.0 * std::max(1.0e-6, temperature_ev) *
                     PhysicsConstants::ElementaryCharge / mass);
}

Particle::ResolvedSpectrum makeDoubleMaxwellSpectrum(const Particle::ParticleType& particle_type,
                                                     double total_density_m3, double cold_temperature_ev,
                                                     double hot_temperature_ev, double hot_fraction,
                                                     double drift_speed_m_per_s, double mass_amu)
{
    Particle::SamplingParameters params;
    params.density = total_density_m3;
    params.bulk_speed = drift_speed_m_per_s;
    params.thermal_speed = thermalSpeedForSourceModel(particle_type, cold_temperature_ev);
    params.hot_thermal_speed = thermalSpeedForSourceModel(particle_type, hot_temperature_ev);
    params.hot_fraction = hot_fraction;
    auto spectrum = Particle::ParticleSource::buildResolvedSpectrum(
        particle_type, Particle::SpatialSamplingModel::DOUBLE_MAXWELL, params,
        Particle::SpectrumUsage::CurrentBalance);
    for (auto& population : spectrum.populations)
    {
        population.mass_amu = mass_amu;
    }
    return spectrum;
}

Particle::ResolvedSpectrum makeSingleMaxwellSpectrum(const Particle::ParticleType& particle_type,
                                                     double density_m3, double temperature_ev,
                                                     double drift_speed_m_per_s, double mass_amu)
{
    Particle::SamplingParameters params;
    params.density = density_m3;
    params.bulk_speed = drift_speed_m_per_s;
    params.thermal_speed = thermalSpeedForSourceModel(particle_type, temperature_ev);
    auto spectrum = Particle::ParticleSource::buildResolvedSpectrum(
        particle_type, Particle::SpatialSamplingModel::SINGLE_MAXWELL, params,
        Particle::SpectrumUsage::CurrentBalance);
    for (auto& population : spectrum.populations)
    {
        population.mass_amu = mass_amu;
    }
    return spectrum;
}

Material::MaterialProperty makeKaptonSurface()
{
    Material::MaterialProperty material(2, Mesh::MaterialType::DIELECTRIC, "kapton");
    material.setPermittivity(3.4);
    material.setConductivity(1.0e-15);
    material.setWorkFunctionEv(4.7);
    material.setSecondaryElectronYield(2.1);
    material.setBreakdownFieldVPerM(2.5e8);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 150.0);
    material.setScalarProperty("secondary_emission_escape_energy_ev", 2.0);
    material.setScalarProperty("photoelectron_yield", 0.016);
    material.setScalarProperty("photoelectron_escape_energy_ev", 2.0);
    material.setScalarProperty("thermionic_escape_energy_ev", 0.2);
    material.setScalarProperty("sims_exponent_n", 1.6);
    material.setScalarProperty("katz_r1", 0.82);
    material.setScalarProperty("katz_n1", 0.37);
    material.setScalarProperty("katz_r2", 0.11);
    material.setScalarProperty("katz_n2", 1.08);
    material.setScalarProperty("ion_secondary_yield", 0.455);
    // Keep the legacy workbook parameterization here: the reference MATLAB
    // script passes `Emi=140` into a function expecting keV.
    material.setScalarProperty("ion_secondary_peak_energy_kev", 140.0);
    material.setScalarProperty("atomic_number", 5.3);
    material.setScalarProperty("body_conductivity_s_per_m", 3.5e7);
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
    material.setScalarProperty("sims_exponent_n", 1.7);
    material.setScalarProperty("katz_r1", 0.95);
    material.setScalarProperty("katz_n1", 0.32);
    material.setScalarProperty("katz_r2", 0.14);
    material.setScalarProperty("katz_n2", 1.12);
    material.setScalarProperty("ion_secondary_yield", 0.08);
    material.setScalarProperty("ion_secondary_peak_energy_kev", 0.45);
    material.setScalarProperty("atomic_number", 9.0);
    material.setScalarProperty("body_conductivity_s_per_m", 3.5e7);
    material.setScalarProperty("poole_frenkel_beta", 7.5e-3);
    material.setScalarProperty("max_field_enhancement_factor", 5.0e8);
    return material;
}

SurfaceChargingConfig makeLeoBaseConfig()
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
    config.patch_incidence_angle_deg = 0.0;
    config.patch_flow_angle_deg = 0.0;
    // Align the default ram-facing LEO environment with the cold Maxwellian
    // reference deck (ref/LEO case 5): 1e4 cm^-3 and 0.1 eV for both species.
    config.plasma.electron_density_m3 = 1.0e10;
    config.plasma.ion_density_m3 = 1.0e10;
    config.plasma.electron_temperature_ev = 0.10;
    config.plasma.ion_temperature_ev = 0.10;
    config.plasma.ion_mass_amu = 16.0;
    config.plasma.neutral_density_m3 = 1.0e14;
    config.plasma.pressure_pa = 1.0e-5;
    config.material = makeKaptonSurface();
    config.emission.surface_temperature_k = 320.0;
    config.emission.enhancement_factor = 1.05;
    config.emission.photon_flux_m2_s = 0.0;
    config.reference_see_model = SecondaryElectronEmissionModel::Sims;
    config.body_floating = true;
    config.body_initial_potential_v = 0.0;
    config.body_capacitance_f = 2.0e-10;
    config.body_photo_current_density_a_per_m2 = 0.0;
    config.patch_photo_current_density_a_per_m2 = 0.0;
    config.max_delta_potential_v_per_step = 10.0;
    config.has_electron_spectrum = true;
    config.has_ion_spectrum = true;
    config.electron_spectrum = makeSingleMaxwellSpectrum(
        Particle::ParticleType::ELECTRON, config.plasma.electron_density_m3,
        config.plasma.electron_temperature_ev, 0.0,
        PhysicsConstants::ElectronMass / PhysicsConstants::AtomicMassUnit);
    config.ion_spectrum =
        makeSingleMaxwellSpectrum(Particle::ParticleType::ION, config.plasma.ion_density_m3,
                                  config.plasma.ion_temperature_ev, 0.0, 16.0);
    return config;
}

SurfaceChargingConfig makeLeoRamReferenceConfig()
{
    SurfaceChargingConfig config = makeLeoBaseConfig();
    config.enable_pic_calibration = false;
    config.enable_live_pic_window = false;
    config.enable_live_pic_mcc = false;
    config.enable_body_patch_circuit = false;
    config.use_reference_current_balance = true;
    return config;
}

SurfaceChargingConfig makeLeoWakeReferenceConfig()
{
    SurfaceChargingConfig config = makeLeoRamReferenceConfig();
    config.flow_alignment_cosine = -1.0;
    config.electron_flow_coupling = 0.002;
    config.patch_flow_angle_deg = 180.0;
    config.electron_spectrum = makeDoubleMaxwellSpectrum(
        Particle::ParticleType::ELECTRON, 3.2e11, 0.20, 1.2, 0.10, 0.0, PhysicsConstants::ElectronMass /
                                                                      PhysicsConstants::AtomicMassUnit);
    config.ion_spectrum =
        makeDoubleMaxwellSpectrum(Particle::ParticleType::ION, 9.0e10, 0.08, 0.22, 0.08, 0.0, 16.0);
    return config;
}

SurfaceChargingConfig makeLeoRamPicCircuitConfig()
{
    SurfaceChargingConfig config = makeLeoRamReferenceConfig();
    config.enable_pic_calibration = true;
    config.enable_live_pic_window = true;
    config.enable_live_pic_mcc = true;
    config.pic_recalibration_interval_steps = 10;
    config.pic_recalibration_trigger_v = 1.5;
    config.live_pic_window_steps = 12;
    config.live_pic_window_layers = 6;
    config.live_pic_particles_per_element = 4;
    config.live_pic_probe_delta_v = 0.5;
    config.live_pic_reference_potential_v = 0.0;
    config.enable_body_patch_circuit = true;
    return config;
}

SurfaceChargingConfig makeLeoWakePicCircuitConfig()
{
    SurfaceChargingConfig config = makeLeoWakeReferenceConfig();
    config.enable_pic_calibration = true;
    config.enable_live_pic_window = true;
    config.enable_live_pic_mcc = true;
    config.pic_recalibration_interval_steps = 10;
    config.pic_recalibration_trigger_v = 1.5;
    config.live_pic_window_steps = 12;
    config.live_pic_window_layers = 6;
    config.live_pic_particles_per_element = 4;
    config.live_pic_probe_delta_v = 0.5;
    config.live_pic_reference_potential_v = 0.0;
    config.enable_body_patch_circuit = true;
    return config;
}

SurfaceChargingConfig makeGeoEcssKaptonReferenceConfig()
{
    SurfaceChargingConfig config;
    config.regime = SurfaceChargingRegime::GeoKineticPicLike;
    config.surface_area_m2 = 1.0;
    config.derive_capacitance_from_material = false;
    config.capacitance_per_area_f_per_m2 = 3.0e-8;
    config.dielectric_thickness_m = 2.5e-4;
    config.floating = true;
    config.use_reference_current_balance = true;
    config.enable_body_patch_circuit = false;
    config.enable_pic_calibration = false;
    config.enable_live_pic_window = false;
    config.enable_live_pic_mcc = false;
    config.body_initial_potential_v = 0.0;
    config.body_floating = false;
    config.body_capacitance_f = 1.0;
    config.material = makeKaptonSurface();
    config.reference_see_model = SecondaryElectronEmissionModel::Whipple;
    config.radiation_conductivity_coefficient = 0.0;
    config.emission.surface_temperature_k = 300.0;
    config.emission.enhancement_factor = 1.0;
    config.emission.photon_flux_m2_s = 0.0;
    config.body_photo_current_density_a_per_m2 = 0.0;
    config.patch_photo_current_density_a_per_m2 = 0.0;
    config.photoelectron_temperature_ev = 2.0;
    config.max_delta_potential_v_per_step = 500.0;
    config.max_abs_potential_v = 5.0e4;
    config.internal_substeps = 1;
    config.electron_collection_coefficient = 1.0 / 6.0;
    config.ion_collection_coefficient = 1.0 / 6.0;
    config.plasma.neutral_density_m3 = 1.0e10;
    config.plasma.pressure_pa = 1.0e-8;
    config.has_electron_spectrum = true;
    config.has_ion_spectrum = true;
    config.electron_spectrum = makeDoubleMaxwellSpectrum(Particle::ParticleType::ELECTRON, 1.4e6,
                                                         400.0, 27500.0, 1.2e6 / 1.4e6, 0.0,
                                                         PhysicsConstants::ElectronMass /
                                                             PhysicsConstants::AtomicMassUnit);
    config.ion_spectrum =
        makeDoubleMaxwellSpectrum(Particle::ParticleType::ION, 1.9e6, 200.0, 28000.0,
                                  1.3e6 / 1.9e6, 0.0, 1.0);
    config.material.setConductivity(0.0);
    config.material.setScalarProperty("body_conductivity_s_per_m", 0.0);
    return config;
}

SurfaceChargingConfig makeGeoEcssKaptonPicCircuitConfig()
{
    SurfaceChargingConfig config = makeGeoEcssKaptonReferenceConfig();
    config.enable_pic_calibration = true;
    config.enable_live_pic_window = true;
    config.enable_live_pic_mcc = true;
    config.enable_body_patch_circuit = true;
    config.body_floating = true;
    config.body_capacitance_f = 3.0e-10;
    config.capacitance_per_area_f_per_m2 = 5.0e-9;
    config.surface_area_m2 = 4.0e-2;
    config.pic_calibration_samples = 6144;
    config.pic_recalibration_interval_steps = 8;
    config.pic_recalibration_trigger_v = 1.0e6;
    config.live_pic_window_steps = 12;
    config.live_pic_window_layers = 6;
    config.live_pic_particles_per_element = 3;
    config.live_pic_probe_delta_v = 5.0;
    config.material.setScalarProperty("body_atomic_number", 13.0);
    config.material.setScalarProperty("body_ion_secondary_yield", 0.244);
    config.material.setScalarProperty("body_ion_secondary_peak_energy_kev", 230.0);
    config.material.setScalarProperty("body_secondary_electron_yield", 0.97);
    config.material.setScalarProperty("body_secondary_yield_peak_energy_ev", 300.0);
    config.material.setScalarProperty("body_sims_exponent_n", 1.621);
    config.material.setScalarProperty("body_katz_r1", 154.0);
    config.material.setScalarProperty("body_katz_n1", 0.80);
    config.material.setScalarProperty("body_katz_r2", 220.0);
    config.material.setScalarProperty("body_katz_n2", 1.76);
    config.material.setScalarProperty("body_electron_collection_coefficient", 2.0 / 3.0);
    config.body_photo_current_density_a_per_m2 = 1.0e-4;
    config.patch_photo_current_density_a_per_m2 = 6.0e-5;
    config.photoelectron_temperature_ev = 3.75;
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
    config.reference_see_model = SecondaryElectronEmissionModel::Katz;
    config.has_electron_spectrum = true;
    config.has_ion_spectrum = true;
    config.electron_spectrum = makeSingleMaxwellSpectrum(
        Particle::ParticleType::ELECTRON, 2.5e15, 8.0, 0.0,
        PhysicsConstants::ElectronMass / PhysicsConstants::AtomicMassUnit);
    config.ion_spectrum = makeSingleMaxwellSpectrum(Particle::ParticleType::ION, 2.0e15, 1.5,
                                                    2.2e4, 131.3);
    return config;
}

void applyReferenceFamilyContract(SurfaceChargingConfig& config,
                                  const std::string& family,
                                  const std::string& case_id,
                                  const std::string& matrix_case_id)
{
    config.reference_family = family;
    config.reference_case_id = case_id;
    config.reference_matrix_case_id = matrix_case_id;
}

void assignPublicBenchmarkSourceForMatrixCase(SurfaceChargingConfig& config,
                                              const std::string& matrix_case_id)
{
    if (config.benchmark_source != SurfaceBenchmarkSource::None)
    {
        return;
    }

    if (matrix_case_id == "geo_ecss_kapton_2000s")
    {
        config.benchmark_source = SurfaceBenchmarkSource::CGeo;
        return;
    }
    if (matrix_case_id == "leo_ram_facing_2000s")
    {
        config.benchmark_source = SurfaceBenchmarkSource::CLeoRam;
        return;
    }
    if (matrix_case_id == "leo_wake_facing_negative_tail_2000s")
    {
        config.benchmark_source = SurfaceBenchmarkSource::CLeoWake;
    }
}

void applySurfacePicContract(SurfaceChargingConfig& config,
                                 const std::string& case_id,
                                 const std::string& matrix_case_id)
{
    const bool enforce_geo_solver_blueprint = matrix_case_id == "geo_ecss_kapton_2000s";

    if (enforce_geo_solver_blueprint && config.solver_config.coupling_mode.empty())
    {
        config.solver_config.coupling_mode = "field_particle_implicit";
    }
    if (enforce_geo_solver_blueprint && config.solver_config.convergence_policy.empty())
    {
        config.solver_config.convergence_policy = "residual_norm_guarded";
    }
    if (config.solver_config.deposition_scheme.empty())
    {
        // Keep the current CIC kernel as the baseline, but pin it explicitly so
        // surface-route presets always resolve through the same deposition path.
        config.solver_config.deposition_scheme = "pic_window_cic";
    }
    if (config.solver_config.collision_set.empty())
    {
        config.solver_config.collision_set = "surface_pic_v1";
    }
    if (config.solver_config.physics_process_set.empty())
    {
        config.solver_config.physics_process_set = "surface_process_core_v1";
    }
    if (enforce_geo_solver_blueprint)
    {
        config.solver_config.max_iterations =
            std::max<std::size_t>(config.solver_config.max_iterations, 64);
        config.solver_config.residual_tolerance =
            std::max(1.0e-9, config.solver_config.residual_tolerance);
        if (!std::isfinite(config.solver_config.relaxation_factor) ||
            config.solver_config.relaxation_factor <= 0.0)
        {
            config.solver_config.relaxation_factor = 0.65;
        }
        config.solver_config.relaxation_factor =
            std::clamp(config.solver_config.relaxation_factor, 0.05, 0.65);
    }

    config.legacy_input_adapter_kind = SurfaceLegacyInputAdapterKind::None;
    config.surface_pic_runtime_kind = SurfacePicRuntimeKind::GraphCoupledSharedSurface;
    config.surface_instrument_set_kind = SurfaceInstrumentSetKind::SurfacePicObserverSet;
    assignPublicBenchmarkSourceForMatrixCase(config, matrix_case_id);
    applyReferenceFamilyContract(config, "surface_pic_reference_family", case_id,
                                 matrix_case_id);
}

void applyUnifiedReferenceContract(SurfaceChargingConfig& config,
                                   const std::string& case_id,
                                   const std::string& matrix_case_id)
{
    config.legacy_input_adapter_kind = SurfaceLegacyInputAdapterKind::None;
    config.surface_pic_runtime_kind = SurfacePicRuntimeKind::LocalWindowSampler;
    config.surface_instrument_set_kind = SurfaceInstrumentSetKind::MetadataOnly;
    assignPublicBenchmarkSourceForMatrixCase(config, matrix_case_id);
    applyReferenceFamilyContract(config, "geo_leo_mat_reference_family", case_id,
                                 matrix_case_id);
}

void applyLegacyDeckContract(SurfaceChargingConfig& config,
                             SurfaceLegacyInputAdapterKind adapter_kind,
                             const std::string& case_id,
                             const std::string& matrix_case_id)
{
    config.legacy_input_adapter_kind = adapter_kind;
    config.surface_pic_runtime_kind = SurfacePicRuntimeKind::LocalWindowSampler;
    config.surface_instrument_set_kind = SurfaceInstrumentSetKind::MetadataOnly;
    applyReferenceFamilyContract(config, "geo_leo_mat_reference_family", case_id,
                                 matrix_case_id);
}

std::array<SurfaceChargingScenarioPreset, 15> buildPresets()
{
    SurfaceChargingScenarioPreset geo_ref;
    geo_ref.name = "geo_ecss_kapton_ref";
    geo_ref.description = "GEO ECSS Kapton reference charging case aligned with the MATLAB workbook.";
    geo_ref.config = makeGeoEcssKaptonReferenceConfig();
    geo_ref.config.runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    geo_ref.config.benchmark_source = SurfaceBenchmarkSource::None;
    geo_ref.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    geo_ref.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applyUnifiedReferenceContract(geo_ref.config, "geo_ecss_kapton_ref",
                                  "geo_ecss_kapton_2000s");
    geo_ref.time_step_s = 1.0 / 3.0;
    geo_ref.steps = 6006;
    geo_ref.default_output_csv = "results/surface_geo_ecss_kapton_ref.csv";

    SurfaceChargingScenarioPreset geo_pic;
    geo_pic.name = "geo_ecss_kapton_pic_circuit";
    geo_pic.description = "GEO ECSS Kapton case with spectrum-driven current balance and PIC-circuit calibration.";
    geo_pic.config = makeGeoEcssKaptonPicCircuitConfig();
    geo_pic.config.runtime_route = SurfaceRuntimeRoute::SurfacePic;
    geo_pic.config.surface_pic_strategy = SurfacePicStrategy::SurfacePicCalibrated;
    geo_pic.config.live_pic_collision_cross_section_set_id = "surface_pic_v1";
    geo_pic.config.benchmark_source = SurfaceBenchmarkSource::None;
    geo_pic.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    geo_pic.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applySurfacePicContract(geo_pic.config, "geo_ecss_kapton_pic_circuit",
                                "geo_ecss_kapton_2000s");
    geo_pic.time_step_s = 5.0e-5;
    geo_pic.steps = 512;
    geo_pic.adaptive_time_stepping = true;
    geo_pic.total_duration_s = 2.56e-2;
    geo_pic.minimum_time_step_s = 1.0e-6;
    geo_pic.maximum_time_step_s = 5.0e-3;
    geo_pic.default_output_csv = "results/surface_geo_ecss_kapton_pic_circuit.csv";

    SurfaceChargingScenarioPreset leo_ref_ram;
    leo_ref_ram.name = "leo_ref_ram_facing";
    leo_ref_ram.description = "LEO reference ram-facing case driven by resolved spectra.";
    leo_ref_ram.config = makeLeoRamReferenceConfig();
    leo_ref_ram.config.runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    leo_ref_ram.config.benchmark_source = SurfaceBenchmarkSource::None;
    leo_ref_ram.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    leo_ref_ram.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applyUnifiedReferenceContract(leo_ref_ram.config, "leo_ref_ram_facing",
                                  "leo_ram_facing_2000s");
    leo_ref_ram.time_step_s = 2.0e-8;
    leo_ref_ram.steps = 600;
    leo_ref_ram.adaptive_time_stepping = true;
    leo_ref_ram.total_duration_s = 2.0e3;
    leo_ref_ram.minimum_time_step_s = 2.0e-8;
    leo_ref_ram.maximum_time_step_s = 5.0;
    leo_ref_ram.default_output_csv = "results/surface_leo_ref_ram_facing.csv";

    SurfaceChargingScenarioPreset leo_ref_wake;
    leo_ref_wake.name = "leo_ref_wake_facing";
    leo_ref_wake.description = "LEO reference wake-facing case driven by resolved spectra.";
    leo_ref_wake.config = makeLeoWakeReferenceConfig();
    leo_ref_wake.config.runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    leo_ref_wake.config.benchmark_source = SurfaceBenchmarkSource::None;
    leo_ref_wake.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    leo_ref_wake.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applyUnifiedReferenceContract(leo_ref_wake.config, "leo_ref_wake_facing",
                                  "leo_wake_facing_negative_tail_2000s");
    leo_ref_wake.time_step_s = 2.0e-8;
    leo_ref_wake.steps = 600;
    leo_ref_wake.adaptive_time_stepping = true;
    leo_ref_wake.total_duration_s = 2.0e3;
    leo_ref_wake.minimum_time_step_s = 2.0e-8;
    leo_ref_wake.maximum_time_step_s = 5.0;
    leo_ref_wake.default_output_csv = "results/surface_leo_ref_wake_facing.csv";

    SurfaceChargingScenarioPreset leo_pic_ram;
    leo_pic_ram.name = "leo_pic_circuit_ram_facing";
    leo_pic_ram.description = "LEO ram-facing case with spectrum-driven reference current balance and PIC calibration.";
    leo_pic_ram.config = makeLeoRamPicCircuitConfig();
    leo_pic_ram.config.runtime_route = SurfaceRuntimeRoute::SurfacePic;
    leo_pic_ram.config.surface_pic_strategy = SurfacePicStrategy::SurfacePicCalibrated;
    leo_pic_ram.config.live_pic_collision_cross_section_set_id = "surface_pic_v1";
    leo_pic_ram.config.benchmark_source = SurfaceBenchmarkSource::None;
    leo_pic_ram.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    leo_pic_ram.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applySurfacePicContract(leo_pic_ram.config, "leo_pic_circuit_ram_facing",
                                "leo_ram_facing_2000s");
    leo_pic_ram.time_step_s = 2.0e-8;
    leo_pic_ram.steps = 600;
    leo_pic_ram.adaptive_time_stepping = true;
    leo_pic_ram.total_duration_s = 2.0e3;
    leo_pic_ram.minimum_time_step_s = 2.0e-8;
    leo_pic_ram.maximum_time_step_s = 5.0;
    leo_pic_ram.default_output_csv = "results/surface_leo_pic_circuit_ram_facing.csv";

    SurfaceChargingScenarioPreset leo_pic_wake;
    leo_pic_wake.name = "leo_pic_circuit_wake_facing";
    leo_pic_wake.description = "LEO wake-facing case with spectrum-driven reference current balance and PIC calibration.";
    leo_pic_wake.config = makeLeoWakePicCircuitConfig();
    leo_pic_wake.config.runtime_route = SurfaceRuntimeRoute::SurfacePic;
    leo_pic_wake.config.surface_pic_strategy = SurfacePicStrategy::SurfacePicCalibrated;
    leo_pic_wake.config.live_pic_collision_cross_section_set_id = "surface_pic_v1";
    leo_pic_wake.config.benchmark_source = SurfaceBenchmarkSource::None;
    leo_pic_wake.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    leo_pic_wake.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applySurfacePicContract(leo_pic_wake.config, "leo_pic_circuit_wake_facing",
                                "leo_wake_facing_negative_tail_2000s");
    leo_pic_wake.time_step_s = 2.0e-8;
    leo_pic_wake.steps = 600;
    leo_pic_wake.adaptive_time_stepping = true;
    leo_pic_wake.total_duration_s = 2.0e3;
    leo_pic_wake.minimum_time_step_s = 2.0e-8;
    leo_pic_wake.maximum_time_step_s = 5.0;
    leo_pic_wake.default_output_csv = "results/surface_leo_pic_circuit_wake_facing.csv";

    SurfaceChargingScenarioPreset geo_pic_direct = geo_pic;
    geo_pic_direct.name = "geo_ecss_kapton_surface_pic_direct";
    geo_pic_direct.description =
        "GEO ECSS Kapton case using the formal surface-PIC direct current route.";
    geo_pic_direct.config.runtime_route = SurfaceRuntimeRoute::SurfacePic;
    geo_pic_direct.config.surface_pic_strategy = SurfacePicStrategy::SurfacePicDirect;
    geo_pic_direct.config.enable_pic_calibration = false;
    applySurfacePicContract(geo_pic_direct.config, "geo_ecss_kapton_surface_pic_direct",
                                "geo_ecss_kapton_2000s");
    geo_pic_direct.default_output_csv = "results/surface_geo_ecss_kapton_surface_pic_direct.csv";

    SurfaceChargingScenarioPreset leo_pic_ram_hybrid = leo_pic_ram;
    leo_pic_ram_hybrid.name = "leo_pic_circuit_ram_facing_hybrid";
    leo_pic_ram_hybrid.description =
        "LEO ram-facing case using the formal surface-PIC hybrid reference route.";
    leo_pic_ram_hybrid.config.runtime_route = SurfaceRuntimeRoute::SurfacePicHybrid;
    leo_pic_ram_hybrid.config.surface_pic_strategy = SurfacePicStrategy::SurfacePicHybridReference;
    applySurfacePicContract(leo_pic_ram_hybrid.config,
                                "leo_pic_circuit_ram_facing_hybrid",
                                "leo_ram_facing_2000s");
    leo_pic_ram_hybrid.default_output_csv =
        "results/surface_leo_pic_circuit_ram_facing_hybrid.csv";

    SurfaceChargingScenarioPreset thruster_plume;
    thruster_plume.name = "thruster_plume_dielectric";
    thruster_plume.description = "Near-thruster plume driven surface current balance case.";
    thruster_plume.config = makeThrusterPlumeConfig();
    thruster_plume.config.runtime_route = SurfaceRuntimeRoute::SCDATUnified;
    thruster_plume.config.benchmark_source = SurfaceBenchmarkSource::None;
    thruster_plume.config.current_algorithm_mode = SurfaceCurrentAlgorithmMode::UnifiedKernelAligned;
    thruster_plume.config.benchmark_mode = SurfaceBenchmarkMode::UnifiedKernelAligned;
    applyReferenceFamilyContract(thruster_plume.config, "thruster_surface_family",
                                 "thruster_plume_dielectric", "thruster_plume_reference");
    thruster_plume.time_step_s = 1.0e-9;
    thruster_plume.steps = 80;
    thruster_plume.default_output_csv = "results/surface_thruster_plume_dielectric.csv";

    SurfaceChargingScenarioPreset geo_ref_legacy = geo_ref;
    geo_ref_legacy.name = "geo_ecss_kapton_ref_legacy_compatible";
    geo_ref_legacy.description =
        "GEO ECSS Kapton reference case running through the legacy-compatible benchmark path.";
    geo_ref_legacy.config.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    geo_ref_legacy.config.benchmark_source = SurfaceBenchmarkSource::CGeo;
    geo_ref_legacy.config.legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    geo_ref_legacy.config.current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    geo_ref_legacy.config.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    applyLegacyDeckContract(geo_ref_legacy.config,
                            SurfaceLegacyInputAdapterKind::CTextReferenceDeck,
                            "geo_ecss_kapton_ref_legacy_compatible",
                            "geo_ecss_kapton_2000s");
    geo_ref_legacy.time_step_s = 1.0e-3;
    geo_ref_legacy.steps = 12;
    geo_ref_legacy.default_output_csv =
        "results/surface_geo_ecss_kapton_ref_legacy_compatible.csv";

    SurfaceChargingScenarioPreset leo_ref_ram_legacy = leo_ref_ram;
    leo_ref_ram_legacy.name = "leo_ref_ram_facing_legacy_compatible";
    leo_ref_ram_legacy.description =
        "LEO ram-facing reference case running through the legacy-compatible benchmark path.";
    leo_ref_ram_legacy.config.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    leo_ref_ram_legacy.config.benchmark_source = SurfaceBenchmarkSource::CLeoRam;
    leo_ref_ram_legacy.config.legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    leo_ref_ram_legacy.config.current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    leo_ref_ram_legacy.config.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    applyLegacyDeckContract(leo_ref_ram_legacy.config,
                            SurfaceLegacyInputAdapterKind::CTextReferenceDeck,
                            "leo_ref_ram_facing_legacy_compatible",
                            "leo_ram_facing_2000s");
    leo_ref_ram_legacy.time_step_s = 1.0;
    leo_ref_ram_legacy.steps = 68;
    leo_ref_ram_legacy.default_output_csv =
        "results/surface_leo_ref_ram_facing_legacy_compatible.csv";

    SurfaceChargingScenarioPreset leo_ref_wake_legacy = leo_ref_wake;
    leo_ref_wake_legacy.name = "leo_ref_wake_facing_legacy_compatible";
    leo_ref_wake_legacy.description =
        "LEO wake-facing reference case running through the legacy-compatible benchmark path.";
    leo_ref_wake_legacy.config.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    leo_ref_wake_legacy.config.benchmark_source = SurfaceBenchmarkSource::CLeoWake;
    leo_ref_wake_legacy.config.legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    leo_ref_wake_legacy.config.current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    leo_ref_wake_legacy.config.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    applyLegacyDeckContract(leo_ref_wake_legacy.config,
                            SurfaceLegacyInputAdapterKind::CTextReferenceDeck,
                            "leo_ref_wake_facing_legacy_compatible",
                            "leo_wake_facing_negative_tail_2000s");
    leo_ref_wake_legacy.time_step_s = 1.0;
    leo_ref_wake_legacy.steps = 68;
    leo_ref_wake_legacy.default_output_csv =
        "results/surface_leo_ref_wake_facing_legacy_compatible.csv";

    SurfaceChargingScenarioPreset geo_ref_matlab = geo_ref;
    geo_ref_matlab.name = "geo_ecss_kapton_ref_matlab_compatible";
    geo_ref_matlab.description =
        "GEO ECSS Kapton reference case running through the MATLAB-compatible benchmark path.";
    geo_ref_matlab.config.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    geo_ref_matlab.config.benchmark_source = SurfaceBenchmarkSource::MatlabGeo;
    geo_ref_matlab.config.legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    geo_ref_matlab.config.current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    geo_ref_matlab.config.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    applyLegacyDeckContract(geo_ref_matlab.config,
                            SurfaceLegacyInputAdapterKind::MatlabReferenceDeck,
                            "geo_ecss_kapton_ref_matlab_compatible",
                            "geo_ecss_kapton_2000s");
    geo_ref_matlab.time_step_s = 3.285188e-2;
    geo_ref_matlab.steps = 60880;
    geo_ref_matlab.default_output_csv =
        "results/surface_geo_ecss_kapton_ref_matlab_compatible.csv";

    SurfaceChargingScenarioPreset leo_ref_ram_matlab = leo_ref_ram;
    leo_ref_ram_matlab.name = "leo_ref_ram_facing_matlab_compatible";
    leo_ref_ram_matlab.description =
        "LEO ram-facing reference case running through the MATLAB-compatible benchmark path.";
    leo_ref_ram_matlab.config.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    leo_ref_ram_matlab.config.benchmark_source = SurfaceBenchmarkSource::MatlabLeo;
    leo_ref_ram_matlab.config.legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    leo_ref_ram_matlab.config.current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    leo_ref_ram_matlab.config.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    applyLegacyDeckContract(leo_ref_ram_matlab.config,
                            SurfaceLegacyInputAdapterKind::MatlabReferenceDeck,
                            "leo_ref_ram_facing_matlab_compatible",
                            "leo_ram_facing_2000s");
    leo_ref_ram_matlab.time_step_s = 1.0;
    leo_ref_ram_matlab.steps = 68;
    leo_ref_ram_matlab.default_output_csv =
        "results/surface_leo_ref_ram_facing_matlab_compatible.csv";

    SurfaceChargingScenarioPreset leo_ref_wake_matlab = leo_ref_wake;
    leo_ref_wake_matlab.name = "leo_ref_wake_facing_matlab_compatible";
    leo_ref_wake_matlab.description =
        "LEO wake-facing reference case running through the MATLAB-compatible benchmark path.";
    leo_ref_wake_matlab.config.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    leo_ref_wake_matlab.config.benchmark_source = SurfaceBenchmarkSource::MatlabLeo;
    leo_ref_wake_matlab.config.legacy_benchmark_execution_mode =
        LegacyBenchmarkExecutionMode::ReplayFromReference;
    leo_ref_wake_matlab.config.current_algorithm_mode =
        SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    leo_ref_wake_matlab.config.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    applyLegacyDeckContract(leo_ref_wake_matlab.config,
                            SurfaceLegacyInputAdapterKind::MatlabReferenceDeck,
                            "leo_ref_wake_facing_matlab_compatible",
                            "leo_wake_facing_negative_tail_2000s");
    leo_ref_wake_matlab.time_step_s = 1.0;
    leo_ref_wake_matlab.steps = 68;
    leo_ref_wake_matlab.default_output_csv =
        "results/surface_leo_ref_wake_facing_matlab_compatible.csv";

    return {geo_ref,
            geo_pic,
            geo_pic_direct,
            leo_ref_ram,
            leo_ref_wake,
            leo_pic_ram,
            leo_pic_ram_hybrid,
            leo_pic_wake,
            thruster_plume,
            geo_ref_legacy,
            leo_ref_ram_legacy,
            leo_ref_wake_legacy,
            geo_ref_matlab,
            leo_ref_ram_matlab,
            leo_ref_wake_matlab};
}

const auto& presets()
{
    static const auto kPresets = buildPresets();
    return kPresets;
}

bool hasPresetNameSuffix(const std::string& value, const std::string& suffix)
{
    if (suffix.size() > value.size())
    {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isLegacyReplayPresetName(const std::string& name)
{
    return hasPresetNameSuffix(name, "_legacy_compatible") ||
           hasPresetNameSuffix(name, "_matlab_compatible");
}

const auto& mainlinePresets()
{
    static const auto kMainlinePresets = []() {
        std::vector<SurfaceChargingScenarioPreset> filtered;
        filtered.reserve(presets().size());
        for (const auto& preset : presets())
        {
            if (!isLegacyReplayPresetName(preset.name))
            {
                filtered.push_back(preset);
            }
        }
        return filtered;
    }();
    return kMainlinePresets;
}

const auto& replayPresets()
{
    static const auto kReplayPresets = []() {
        std::vector<SurfaceChargingScenarioPreset> filtered;
        filtered.reserve(presets().size());
        for (const auto& preset : presets())
        {
            if (isLegacyReplayPresetName(preset.name))
            {
                filtered.push_back(preset);
            }
        }
        return filtered;
    }();
    return kReplayPresets;
}

bool tryGetPresetFromCollection(const std::string& canonical_name,
                                const std::vector<SurfaceChargingScenarioPreset>& candidates,
                                SurfaceChargingScenarioPreset& preset)
{
    for (const auto& candidate : candidates)
    {
        if (candidate.name == canonical_name)
        {
            preset = candidate;
            return true;
        }
    }
    return false;
}

bool tryResolveAlias(const std::string& name, std::string& canonical_name)
{
    if (name == "geo_eclipse_dielectric")
    {
        canonical_name = "geo_ecss_kapton_pic_circuit";
        return true;
    }
    if (name == "leo_daylight_kapton")
    {
        canonical_name = "leo_pic_circuit_ram_facing";
        return true;
    }
    if (name == "leo_daylight_wake_kapton")
    {
        canonical_name = "leo_pic_circuit_wake_facing";
        return true;
    }
    return false;
}

} // namespace

std::vector<std::string> listSurfaceChargingScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(mainlinePresets().size() + 3);
    for (const auto& preset : mainlinePresets())
    {
        names.push_back(preset.name);
    }
    names.push_back("geo_eclipse_dielectric");
    names.push_back("leo_daylight_kapton");
    names.push_back("leo_daylight_wake_kapton");
    return names;
}

std::vector<std::string> listSurfaceChargingReplayScenarioPresetNames()
{
    std::vector<std::string> names;
    names.reserve(replayPresets().size());
    for (const auto& preset : replayPresets())
    {
        names.push_back(preset.name);
    }
    return names;
}

bool tryGetSurfaceChargingMainlineScenarioPreset(const std::string& name,
                                                 SurfaceChargingScenarioPreset& preset)
{
    std::string canonical_name = name;
    tryResolveAlias(name, canonical_name);
    return tryGetPresetFromCollection(canonical_name, mainlinePresets(), preset);
}

bool tryGetSurfaceChargingReplayScenarioPreset(const std::string& name,
                                               SurfaceChargingScenarioPreset& preset)
{
    std::string canonical_name = name;
    return tryGetPresetFromCollection(canonical_name, replayPresets(), preset);
}

bool tryGetSurfaceChargingScenarioPreset(const std::string& name,
                                         SurfaceChargingScenarioPreset& preset)
{
    if (tryGetSurfaceChargingMainlineScenarioPreset(name, preset))
    {
        return true;
    }
    return tryGetSurfaceChargingReplayScenarioPreset(name, preset);
}

SurfaceChargingScenarioPreset makeDefaultSurfaceChargingScenarioPreset()
{
    SurfaceChargingScenarioPreset preset;
    tryGetSurfaceChargingMainlineScenarioPreset("leo_daylight_kapton", preset);
    return preset;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT


