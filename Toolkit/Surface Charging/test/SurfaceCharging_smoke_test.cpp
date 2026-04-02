#include "DensePlasmaSurfaceCharging.h"
#include "LegacyBenchmarkSupport.h"
#include "ReferenceCurrentBalanceModel.h"
#include "SurfaceChargingCases.h"

#include <gtest/gtest.h>

#include <fstream>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <limits>
#include <sstream>
#include <map>
#include <vector>

using SCDAT::Toolkit::SurfaceCharging::DensePlasmaSurfaceCharging;
using SCDAT::Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset;

namespace
{

struct CurveSample
{
    double time_ms = 0.0;
    double potential_v = 0.0;
};

std::filesystem::path resolveFromRepoRoot(const std::filesystem::path& path)
{
    if (path.empty() || path.is_absolute())
    {
        return path;
    }

    const auto source_path = std::filesystem::path(__FILE__);
    const auto repo_root = source_path.parent_path().parent_path().parent_path().parent_path();
    return std::filesystem::weakly_canonical(repo_root / path);
}

std::vector<std::string> splitCsvLine(const std::string& line)
{
    std::vector<std::string> values;
    std::stringstream stream(line);
    std::string token;
    while (std::getline(stream, token, ','))
    {
        values.push_back(token);
    }
    return values;
}

std::size_t findColumnIndex(const std::vector<std::string>& header, const std::string& name)
{
    for (std::size_t i = 0; i < header.size(); ++i)
    {
        if (header[i] == name)
        {
            return i;
        }
    }
    return header.size();
}

std::vector<CurveSample> readReferenceDatCurve(const std::filesystem::path& path)
{
    std::ifstream input(resolveFromRepoRoot(path));
    std::vector<CurveSample> curve;
    std::string line;
    while (std::getline(input, line))
    {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
        {
            continue;
        }
        const char c = line[first];
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+'))
        {
            continue;
        }

        std::istringstream stream(line);
        double cycle = 0.0;
        CurveSample sample;
        if (stream >> cycle >> sample.potential_v >> sample.time_ms)
        {
            curve.push_back(sample);
        }
    }
    return curve;
}

std::vector<CurveSample> readBenchmarkCsvCurve(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::vector<CurveSample> curve;
    std::string header;
    std::getline(input, header);
    const auto header_columns = splitCsvLine(header);
    const auto time_index = findColumnIndex(header_columns, "Time_ms");
    const auto potential_index = findColumnIndex(header_columns, "Vs");
    if (time_index >= header_columns.size() || potential_index >= header_columns.size())
    {
        return curve;
    }

    std::string row;
    while (std::getline(input, row))
    {
        if (row.empty())
        {
            continue;
        }
        const auto values = splitCsvLine(row);
        if (values.size() <= std::max(time_index, potential_index))
        {
            continue;
        }
        curve.push_back({std::stod(values[time_index]), std::stod(values[potential_index])});
    }
    return curve;
}

double curveRmse(const std::vector<CurveSample>& lhs, const std::vector<CurveSample>& rhs)
{
    const std::size_t count = std::min(lhs.size(), rhs.size());
    if (count == 0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < count; ++i)
    {
        const double delta = lhs[i].potential_v - rhs[i].potential_v;
        sum += delta * delta;
    }
    return std::sqrt(sum / static_cast<double>(count));
}

std::map<std::string, std::string> readKeyValueReport(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(input, line))
    {
        const auto pos = line.find('=');
        if (pos == std::string::npos)
        {
            continue;
        }
        values.emplace(line.substr(0, pos), line.substr(pos + 1));
    }
    return values;
}

double readReportDoubleOrDefault(const std::map<std::string, std::string>& values,
                                 const std::string& key,
                                 double fallback)
{
    const auto it = values.find(key);
    if (it == values.end() || it->second.empty())
    {
        return fallback;
    }
    try
    {
        return std::stod(it->second);
    }
    catch (const std::exception&)
    {
        return fallback;
    }
}

SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceConfig
buildOneDimensionalReferenceConfig(
    const SCDAT::Toolkit::SurfaceCharging::SurfaceChargingConfig& config,
    double effective_conductivity_s_per_m)
{
    constexpr double kElementaryCharge = 1.602176634e-19;

    SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceConfig reference_config;
    reference_config.plasma = config.plasma;
    reference_config.electron_spectrum = config.electron_spectrum;
    reference_config.ion_spectrum = config.ion_spectrum;
    reference_config.has_electron_spectrum = config.has_electron_spectrum;
    reference_config.has_ion_spectrum = config.has_ion_spectrum;
    reference_config.patch_material = config.material;
    reference_config.see_model = config.reference_see_model;
    reference_config.patch_incidence_angle_deg = config.patch_incidence_angle_deg;
    reference_config.patch_flow_angle_deg = config.patch_flow_angle_deg;
    reference_config.patch_thickness_m = config.dielectric_thickness_m;
    reference_config.patch_conductivity_s_per_m = effective_conductivity_s_per_m;
    reference_config.electron_collection_coefficient = config.electron_collection_coefficient;
    reference_config.ion_collection_coefficient = config.ion_collection_coefficient;
    reference_config.bulk_flow_velocity_m_per_s = config.bulk_flow_velocity_m_per_s;
    reference_config.flow_alignment_cosine = config.flow_alignment_cosine;
    reference_config.ion_directed_velocity_m_per_s = config.ion_directed_velocity_m_per_s;
    reference_config.electron_calibration_factor = 1.0;
    reference_config.ion_calibration_factor = 1.0;
    reference_config.photoelectron_temperature_ev =
        std::max(1.0e-3, config.photoelectron_temperature_ev);
    reference_config.enable_ram_current =
        config.regime == SCDAT::Toolkit::SurfaceCharging::SurfaceChargingRegime::LeoFlowingPlasma;

    const double derived_photo_current =
        kElementaryCharge * config.emission.photon_flux_m2_s *
        std::clamp(config.material.getScalarProperty("photoelectron_yield", 0.0) *
                       config.emission.enhancement_factor,
                   0.0, 1.0);
    reference_config.body_photo_current_density_a_per_m2 =
        config.body_photo_current_density_a_per_m2 > 0.0
            ? config.body_photo_current_density_a_per_m2
            : derived_photo_current;
    reference_config.patch_photo_current_density_a_per_m2 =
        config.patch_photo_current_density_a_per_m2 > 0.0
            ? config.patch_photo_current_density_a_per_m2
            : derived_photo_current;

    reference_config.body_material = config.material;
    reference_config.body_material.setType(SCDAT::Mesh::MaterialType::CONDUCTOR);
    reference_config.body_material.setName(config.material.getName() + "_body");
    reference_config.body_material.setConductivity(std::max(
        1.0e-6, config.material.getScalarProperty("body_conductivity_s_per_m", 1.0e4)));
    return reference_config;
}

} // namespace

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkSupportParsesCGeoDefinition)
{
    const auto definition =
        SCDAT::Toolkit::SurfaceCharging::loadLegacyBenchmarkCaseDefinition(
            SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::CGeo);

    ASSERT_EQ(definition.input.raw_values.size(), 10u);
    EXPECT_DOUBLE_EQ(definition.input.raw_values.front(), 1.0);
    EXPECT_DOUBLE_EQ(definition.input.raw_values[2], 30.0);
    ASSERT_FALSE(definition.environments.empty());
    EXPECT_EQ(definition.environments.front().case_id, 1);
    EXPECT_NE(definition.environments.front().name.find("ECSS"), std::string::npos);
    ASSERT_FALSE(definition.structure_materials.empty());
    EXPECT_EQ(definition.structure_materials.front().name, "Aluminium");
    ASSERT_FALSE(definition.dielectric_patch_materials.empty());
    EXPECT_EQ(definition.dielectric_patch_materials.front().name, "Kapton");
    EXPECT_FALSE(definition.patch_reference_curve.empty());
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkSupportParsesCLeoDefinition)
{
    const auto definition =
        SCDAT::Toolkit::SurfaceCharging::loadLegacyBenchmarkCaseDefinition(
            SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::CLeoRam);

    ASSERT_EQ(definition.input.raw_values.size(), 13u);
    EXPECT_DOUBLE_EQ(definition.input.raw_values[2], 800.0);
    ASSERT_FALSE(definition.environments.empty());
    EXPECT_EQ(definition.environments.front().case_id, 1);
    EXPECT_NE(definition.environments.front().name.find("IRI"), std::string::npos);
    ASSERT_FALSE(definition.metal_patch_materials.empty());
    EXPECT_EQ(definition.metal_patch_materials.front().material_id, 6);
    EXPECT_EQ(definition.metal_patch_materials.front().name, "Conductive white paint PCB-Z");
    EXPECT_FALSE(definition.patch_reference_curve.empty());
    EXPECT_FALSE(definition.body_reference_curve.empty());
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkSupportResolvesMatlabLeoPaths)
{
    SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkPaths paths;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::resolveLegacyBenchmarkPaths(
        SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::MatlabLeo, paths));

    EXPECT_TRUE(std::filesystem::exists(paths.input_path));
    EXPECT_TRUE(std::filesystem::exists(paths.environment_path));
    EXPECT_TRUE(std::filesystem::exists(paths.patch_reference_curve_path));
    EXPECT_TRUE(std::filesystem::exists(paths.body_reference_curve_path));
    EXPECT_TRUE(std::filesystem::exists(paths.matlab_generator_path));
    EXPECT_NE(paths.matlab_generator_path.string().find(
                  "reconstruct_legacy_leo_initial_body_current.py"),
              std::string::npos);
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkSupportComputesZeroMetricsForSelfComparison)
{
    const auto reference_curve =
        SCDAT::Toolkit::SurfaceCharging::loadLegacyBenchmarkCurve(
            "ref/GEO/reference_result_patch.dat", true);
    const auto metrics = SCDAT::Toolkit::SurfaceCharging::computeLegacyBenchmarkMetrics(
        reference_curve, reference_curve);

    ASSERT_TRUE(metrics.valid);
    EXPECT_EQ(metrics.compared_sample_count, reference_curve.size());
    EXPECT_NEAR(metrics.rmse_v, 0.0, 1.0e-12);
    EXPECT_NEAR(metrics.terminal_potential_delta_v, 0.0, 1.0e-12);
    EXPECT_NEAR(metrics.terminal_time_delta_s, 0.0, 1.0e-12);
    EXPECT_NEAR(metrics.time_to_equilibrium_delta_s, 0.0, 1.0e-12);
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkExecutionConfigAppliesGeoInputs)
{
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;

    const auto translated =
        SCDAT::Toolkit::SurfaceCharging::applyLegacyBenchmarkExecutionConfig(preset.config);

    EXPECT_EQ(translated.regime,
              SCDAT::Toolkit::SurfaceCharging::SurfaceChargingRegime::GeoKineticPicLike);
    EXPECT_NEAR(translated.dielectric_thickness_m, 1.0e-3, 1.0e-12);
    EXPECT_NEAR(translated.patch_incidence_angle_deg, 30.0, 1.0e-12);
    EXPECT_EQ(translated.reference_see_model,
              SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel::Whipple);
    EXPECT_EQ(translated.material.getName(), "Aluminium");
    ASSERT_FALSE(translated.patch_physics_overrides.empty());
    EXPECT_EQ(translated.patch_physics_overrides.front().material.getName(), "Aluminium");
    EXPECT_NEAR(translated.plasma.electron_density_m3, 1.4e6, 1.0);
    EXPECT_GT(translated.body_photo_current_density_a_per_m2, 0.0);
    EXPECT_GT(translated.patch_photo_current_density_a_per_m2, 0.0);
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkExecutionConfigAppliesLeoInputs)
{
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;

    const auto translated =
        SCDAT::Toolkit::SurfaceCharging::applyLegacyBenchmarkExecutionConfig(preset.config);

    EXPECT_EQ(translated.regime,
              SCDAT::Toolkit::SurfaceCharging::SurfaceChargingRegime::LeoFlowingPlasma);
    EXPECT_NEAR(translated.dielectric_thickness_m, 1.0e-3, 1.0e-12);
    EXPECT_NEAR(translated.patch_incidence_angle_deg, 90.0, 1.0e-12);
    EXPECT_NEAR(translated.bulk_flow_velocity_m_per_s, 0.0, 1.0e-12);
    EXPECT_EQ(translated.material.getName(), "Aluminium");
    ASSERT_FALSE(translated.patch_physics_overrides.empty());
    EXPECT_EQ(translated.patch_physics_overrides.front().material.getName(), "Aluminium");
    EXPECT_GT(translated.plasma.ion_density_m3, 1.0e9);
    EXPECT_NEAR(translated.body_photo_current_density_a_per_m2, 0.0, 1.0e-18);
    EXPECT_NEAR(translated.patch_photo_current_density_a_per_m2, 0.0, 1.0e-18);
}

TEST(SurfaceChargingSmokeTest, ToolkitInitializesAdvancesAndExports)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_pic_circuit", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    for (int i = 0; i < 8; ++i)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path = std::filesystem::temp_directory_path() / "surface_charging_smoke.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));
    EXPECT_TRUE(std::filesystem::exists(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    EXPECT_NE(header.find("floating_equilibrium_potential_v"), std::string::npos);
    EXPECT_NE(header.find("net_current_density_a_per_m2"), std::string::npos);
    EXPECT_NE(header.find("effective_conductivity_s_per_m"), std::string::npos);
    EXPECT_NE(header.find("adaptive_time_step_s"), std::string::npos);
    EXPECT_NE(header.find("normal_electric_field_v_per_m"), std::string::npos);
    EXPECT_NE(header.find("local_charge_density_c_per_m3"), std::string::npos);
    EXPECT_NE(header.find("electron_pic_calibration_factor"), std::string::npos);
    EXPECT_NE(header.find("body_potential_v"), std::string::npos);
    EXPECT_NE(header.find("circuit_branch_current_a"), std::string::npos);
    EXPECT_NE(header.find("current_derivative_a_per_m2_per_v"), std::string::npos);
    EXPECT_NE(header.find("live_pic_mcc_enabled"), std::string::npos);
    EXPECT_NE(header.find("current_algorithm_mode_id"), std::string::npos);
    EXPECT_NE(header.find("benchmark_mode_id"), std::string::npos);
    EXPECT_NE(header.find("runtime_route_id"), std::string::npos);
    EXPECT_NE(header.find("benchmark_source_id"), std::string::npos);
    EXPECT_NE(header.find("benchmark_execution_mode_id"), std::string::npos);
    EXPECT_NE(header.find("surface_circuit_node_count"), std::string::npos);

    EXPECT_GE(charging.getStatus().steps_completed, 8u);
    EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
    EXPECT_TRUE(std::isfinite(charging.getStatus().currents.total_current_a_per_m2));
}

TEST(SurfaceChargingSmokeTest, ThrusterPlumePresetRemainsFinite)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "thruster_plume_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    for (int i = 0; i < 12; ++i)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
        EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
        EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_charge_density_c_per_m2));
        EXPECT_TRUE(std::isfinite(charging.getStatus().currents.total_current_a_per_m2));
    }
}

TEST(SurfaceChargingSmokeTest, GeoEclipseFloatingPotentialIsFiniteAndNegative)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const double floating_potential = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(floating_potential));
    EXPECT_LT(floating_potential, -1.0e3);
    EXPECT_GT(floating_potential, -5.0e4);
}

TEST(SurfaceChargingSmokeTest, LeoDaylightFloatingPotentialIsModeratelyPositive)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const double floating_potential = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(floating_potential));
    EXPECT_GT(floating_potential, 0.0);
    EXPECT_LT(floating_potential, 20.0);
}

TEST(SurfaceChargingSmokeTest, LeoWakePotentialExceedsRamFacingPotential)
{
    DensePlasmaSurfaceCharging ram_facing;
    DensePlasmaSurfaceCharging wake_facing;
    SurfaceChargingScenarioPreset ram_preset;
    SurfaceChargingScenarioPreset wake_preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", ram_preset));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_wake_facing", wake_preset));

    ASSERT_TRUE(ram_facing.initialize(ram_preset.config));
    ASSERT_TRUE(wake_facing.initialize(wake_preset.config));

    const double ram_potential = ram_facing.computeFloatingPotential();
    const double wake_potential = wake_facing.computeFloatingPotential();

    EXPECT_TRUE(std::isfinite(ram_potential));
    EXPECT_TRUE(std::isfinite(wake_potential));
    EXPECT_GT(wake_potential, ram_potential);
}

TEST(SurfaceChargingSmokeTest,
     ThickKaptonPotentialMatchesOneDimensionalReferenceModel)
{
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref", preset));

    auto config = preset.config;
    config.runtime_route =
        SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::SCDATUnified;
    config.benchmark_source =
        SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::None;
    config.current_algorithm_mode =
        SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode::UnifiedSpisAligned;
    config.benchmark_mode =
        SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode::UnifiedSpisAligned;
    config.use_reference_current_balance = true;
    config.enable_pic_calibration = false;
    config.enable_live_pic_window = false;
    config.enable_live_pic_mcc = false;
    config.enable_body_patch_circuit = false;

    config.surface_area_m2 = 2.0e-2;
    config.derive_capacitance_from_material = true;
    config.dielectric_thickness_m = 7.5e-4;
    config.material.setPermittivity(3.4);
    config.material.setConductivity(5.0e-9);
    config.material.setScalarProperty("poole_frenkel_beta", 0.0);
    config.material.setScalarProperty("max_field_enhancement_factor", 1.0);
    config.radiation_conductivity_coefficient = 0.0;
    config.patch_photo_current_density_a_per_m2 = 0.0;
    config.body_photo_current_density_a_per_m2 = 0.0;
    config.emission.photon_flux_m2_s = 0.0;
    config.body_floating = false;
    config.body_initial_potential_v = 0.0;
    config.max_abs_potential_v = 2.0e4;

    DensePlasmaSurfaceCharging charging;
    ASSERT_TRUE(charging.initialize(config));

    const double advanced_surface_potential_v = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(advanced_surface_potential_v));

    SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceModel one_dimensional_model;
    const auto reference_config =
        buildOneDimensionalReferenceConfig(config, config.material.getConductivity());
    ASSERT_TRUE(one_dimensional_model.configure(reference_config));

    const double one_dimensional_surface_potential_v =
        one_dimensional_model.solveEquilibriumPotential(
            SCDAT::Toolkit::SurfaceCharging::ReferenceSurfaceRole::Patch,
            config.body_initial_potential_v, config.body_initial_potential_v,
            -config.max_abs_potential_v, config.max_abs_potential_v);
    EXPECT_TRUE(std::isfinite(one_dimensional_surface_potential_v));

    const double tolerance_v =
        std::max(25.0, 0.12 * std::abs(one_dimensional_surface_potential_v));
    EXPECT_NEAR(advanced_surface_potential_v,
                one_dimensional_surface_potential_v,
                tolerance_v);
}

TEST(SurfaceChargingSmokeTest, LegacyAliasesResolveToNewSpectrumPresets)
{
    SurfaceChargingScenarioPreset geo_alias;
    SurfaceChargingScenarioPreset leo_alias;
    SurfaceChargingScenarioPreset wake_alias;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_eclipse_dielectric", geo_alias));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_daylight_kapton", leo_alias));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_daylight_wake_kapton", wake_alias));

    EXPECT_TRUE(geo_alias.config.has_electron_spectrum);
    EXPECT_TRUE(geo_alias.config.enable_pic_calibration);
    EXPECT_TRUE(leo_alias.config.has_electron_spectrum);
    EXPECT_TRUE(wake_alias.config.has_ion_spectrum);
}

TEST(SurfaceChargingSmokeTest, LegacyCompatibleBenchmarkPresetsResolve)
{
    SurfaceChargingScenarioPreset geo_legacy;
    SurfaceChargingScenarioPreset ram_legacy;
    SurfaceChargingScenarioPreset wake_legacy;
    SurfaceChargingScenarioPreset matlab_legacy;
    SurfaceChargingScenarioPreset ram_matlab_legacy;
    SurfaceChargingScenarioPreset wake_matlab_legacy;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", geo_legacy));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing_legacy_compatible", ram_legacy));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_wake_facing_legacy_compatible", wake_legacy));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_matlab_compatible", matlab_legacy));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing_matlab_compatible", ram_matlab_legacy));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_wake_facing_matlab_compatible", wake_matlab_legacy));

    EXPECT_EQ(geo_legacy.config.current_algorithm_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode::LegacyRefCompatible);
    EXPECT_EQ(ram_legacy.config.current_algorithm_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode::LegacyRefCompatible);
    EXPECT_EQ(wake_legacy.config.current_algorithm_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode::LegacyRefCompatible);
    EXPECT_EQ(geo_legacy.config.benchmark_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode::LegacyRefCompatible);
    EXPECT_EQ(ram_legacy.config.benchmark_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode::LegacyRefCompatible);
    EXPECT_EQ(wake_legacy.config.benchmark_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode::LegacyRefCompatible);
    EXPECT_EQ(geo_legacy.config.runtime_route,
              SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark);
    EXPECT_EQ(geo_legacy.config.benchmark_source,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::CGeo);
    EXPECT_EQ(ram_legacy.config.benchmark_source,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::CLeoRam);
    EXPECT_EQ(wake_legacy.config.benchmark_source,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::CLeoWake);
    EXPECT_EQ(matlab_legacy.config.runtime_route,
              SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark);
    EXPECT_EQ(matlab_legacy.config.benchmark_source,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::MatlabGeo);
    EXPECT_EQ(ram_matlab_legacy.config.runtime_route,
              SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark);
    EXPECT_EQ(wake_matlab_legacy.config.runtime_route,
              SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute::LegacyBenchmark);
    EXPECT_EQ(ram_matlab_legacy.config.current_algorithm_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode::LegacyRefCompatible);
    EXPECT_EQ(wake_matlab_legacy.config.current_algorithm_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode::LegacyRefCompatible);
    EXPECT_EQ(ram_matlab_legacy.config.benchmark_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode::LegacyRefCompatible);
    EXPECT_EQ(wake_matlab_legacy.config.benchmark_mode,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode::LegacyRefCompatible);
    EXPECT_EQ(ram_matlab_legacy.config.benchmark_source,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::MatlabLeo);
    EXPECT_EQ(wake_matlab_legacy.config.benchmark_source,
              SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkSource::MatlabLeo);
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkRoutesExportDistinctBenchmarkSourceIds)
{
    auto export_and_read_source_id = [](
                                         const std::string& preset_name,
                                         const std::filesystem::path& csv_path) {
        DensePlasmaSurfaceCharging charging;
        SurfaceChargingScenarioPreset preset;
        EXPECT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
            preset_name, preset));
        EXPECT_TRUE(charging.initialize(preset.config));
        EXPECT_TRUE(charging.advance(preset.time_step_s));
        EXPECT_TRUE(charging.exportResults(csv_path));

        std::ifstream input(csv_path);
        EXPECT_TRUE(input.is_open());
        std::string header;
        std::getline(input, header);
        const auto header_columns = splitCsvLine(header);
        const std::size_t source_index =
            findColumnIndex(header_columns, "benchmark_source_id");
        EXPECT_LT(source_index, header_columns.size());

        std::string row;
        std::getline(input, row);
        const auto values = splitCsvLine(row);
        EXPECT_GT(values.size(), source_index);
        return std::stod(values[source_index]);
    };

    const auto c_csv =
        std::filesystem::temp_directory_path() / "surface_charging_legacy_c_benchmark.csv";
    const auto matlab_csv =
        std::filesystem::temp_directory_path() / "surface_charging_legacy_matlab_benchmark.csv";

    const double c_source_id =
        export_and_read_source_id("geo_ecss_kapton_ref_legacy_compatible", c_csv);
    const double matlab_source_id =
        export_and_read_source_id("geo_ecss_kapton_ref_matlab_compatible", matlab_csv);

    EXPECT_DOUBLE_EQ(c_source_id, 1.0);
    EXPECT_DOUBLE_EQ(matlab_source_id, 4.0);
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkExecuteModeAdvancesWithoutReplay)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
    preset.steps = 4;

    ASSERT_TRUE(charging.initialize(preset.config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
    EXPECT_TRUE(std::isfinite(charging.getStatus().currents.total_current_a_per_m2));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_legacy_execute.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    const auto header_columns = splitCsvLine(header);
    const auto mode_index =
        findColumnIndex(header_columns, "benchmark_execution_mode_id");
    ASSERT_LT(mode_index, header_columns.size());

    std::string row;
    std::getline(input, row);
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), mode_index);
    EXPECT_DOUBLE_EQ(std::stod(values[mode_index]), 1.0);
}

TEST(SurfaceChargingSmokeTest, LegacyBenchmarkExportsRefStyleColumns)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_legacy_ref_style.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);

    EXPECT_NE(header.find("Cycle"), std::string::npos);
    EXPECT_NE(header.find("Vs"), std::string::npos);
    EXPECT_NE(header.find("Time_ms"), std::string::npos);
    EXPECT_NE(header.find("Jnet"), std::string::npos);
    EXPECT_NE(header.find("Je"), std::string::npos);
    EXPECT_NE(header.find("Jse"), std::string::npos);
    EXPECT_NE(header.find("Jb"), std::string::npos);
    EXPECT_NE(header.find("Ji"), std::string::npos);
    EXPECT_NE(header.find("Jsi"), std::string::npos);
    EXPECT_NE(header.find("Jph"), std::string::npos);
    EXPECT_NE(header.find("Jcond"), std::string::npos);
    EXPECT_NE(header.find("benchmark_time_ms"), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, MatlabBenchmarkExportWritesBaselineMetadata)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_matlab_compatible", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_matlab_metadata.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::filesystem::path h5_path = csv_path;
    h5_path.replace_extension(".h5");
    std::ifstream input(h5_path);
    ASSERT_TRUE(input.is_open());
    const std::string archive((std::istreambuf_iterator<char>(input)),
                              std::istreambuf_iterator<char>());

    EXPECT_NE(archive.find("benchmark_baseline_family=MATLAB"), std::string::npos);
    EXPECT_NE(archive.find("benchmark_baseline_origin=GeneratedFromMatlabReferenceScript"),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_matlab_generator_path="), std::string::npos);
    EXPECT_NE(archive.find("generate_reference_result_geo.py"),
              std::string::npos);
}

TEST(SurfaceChargingSmokeTest, MatlabLeoBenchmarkExportWritesBaselineMetadata)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing_matlab_compatible", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_matlab_leo_metadata.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::filesystem::path h5_path = csv_path;
    h5_path.replace_extension(".h5");
    std::ifstream h5_input(h5_path);
    ASSERT_TRUE(h5_input.is_open());
    const std::string archive((std::istreambuf_iterator<char>(h5_input)),
                              std::istreambuf_iterator<char>());

    EXPECT_NE(archive.find("benchmark_baseline_family=MATLAB"), std::string::npos);
    EXPECT_NE(archive.find("benchmark_baseline_origin=LegacyLeoReferenceWithMatlabDiagnostics"),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_matlab_generator_path="), std::string::npos);
    EXPECT_NE(archive.find("benchmark_acceptance_contract_id=matlab-leo-v1"),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_body_negative_tail_je_rmse_a_per_m2="),
              std::string::npos);
    EXPECT_NE(archive.find("reconstruct_legacy_leo_initial_body_current.py"),
              std::string::npos);

    std::filesystem::path report_path = csv_path;
    report_path.replace_extension(".benchmark.txt");
    std::ifstream report_input(report_path);
    ASSERT_TRUE(report_input.is_open());
    const std::string report((std::istreambuf_iterator<char>(report_input)),
                             std::istreambuf_iterator<char>());

    EXPECT_NE(report.find("benchmark_source=MATLAB-LEO"), std::string::npos);
    EXPECT_NE(report.find("baseline_origin=LegacyLeoReferenceWithMatlabDiagnostics"),
              std::string::npos);
    EXPECT_NE(report.find("acceptance_contract_id=matlab-leo-v1"), std::string::npos);
    EXPECT_NE(report.find("acceptance_focus_segment=negative_potential_tail"),
              std::string::npos);
}

TEST(SurfaceChargingSmokeTest, LegacyCBenchmarkExportWritesComparisonMetadata)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_c_metadata.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::filesystem::path h5_path = csv_path;
    h5_path.replace_extension(".h5");
    std::ifstream input(h5_path);
    ASSERT_TRUE(input.is_open());
    const std::string archive((std::istreambuf_iterator<char>(input)),
                              std::istreambuf_iterator<char>());

    EXPECT_NE(archive.find("benchmark_baseline_family=LegacyC"), std::string::npos);
    EXPECT_NE(archive.find("benchmark_baseline_origin=ReferenceDatTable"), std::string::npos);
    EXPECT_NE(archive.find("benchmark_reference_patch_curve="), std::string::npos);
    EXPECT_NE(archive.find("reference_result_patch.dat"),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_acceptance_contract_id=legacy-c-geo-v1"),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_acceptance_patch_rmse_v_max="),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_acceptance_gate_status="),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_acceptance_gate_pass="),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_body_negative_tail_jnet_rmse_a_per_m2="),
              std::string::npos);
    EXPECT_NE(archive.find("benchmark_rmse_v="), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, BenchmarkExportWritesReadableSidecarReport)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_benchmark_report.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    auto report_path = csv_path;
    report_path.replace_extension(".benchmark.txt");
    std::ifstream input(report_path);
    ASSERT_TRUE(input.is_open());
    const std::string report((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());

    EXPECT_NE(report.find("runtime_route=LegacyBenchmark"), std::string::npos);
    EXPECT_NE(report.find("benchmark_source=C-GEO"), std::string::npos);
    EXPECT_NE(report.find("execution_mode=ReplayFromReference"), std::string::npos);
    EXPECT_NE(report.find("baseline_family=LegacyC"), std::string::npos);
    EXPECT_NE(report.find("patch_reference_curve="), std::string::npos);
    EXPECT_NE(report.find("reference_result_patch.dat"),
              std::string::npos);
    EXPECT_NE(report.find("patch_rmse_v="), std::string::npos);
    EXPECT_NE(report.find("body_rmse_v="), std::string::npos);
    EXPECT_NE(report.find("body_je_rmse_a_per_m2="), std::string::npos);
    EXPECT_NE(report.find("consistency_status="), std::string::npos);
    EXPECT_NE(report.find("consistency_authority="), std::string::npos);
    EXPECT_NE(report.find("body_je_mean_signed_delta_a_per_m2="), std::string::npos);
    EXPECT_NE(report.find("body_je_initial_efolding_energy_ev_actual="), std::string::npos);
    EXPECT_NE(report.find("body_je_initial_delta_a_per_m2="), std::string::npos);
    EXPECT_NE(report.find("body_je_initial_ratio_actual_to_reference="), std::string::npos);
    EXPECT_NE(report.find("acceptance_contract_version=legacy-benchmark-v1"),
              std::string::npos);
    EXPECT_NE(report.find("acceptance_contract_id=legacy-c-geo-v1"),
              std::string::npos);
    EXPECT_NE(report.find("acceptance_patch_rmse_v_max="), std::string::npos);
    EXPECT_NE(report.find("acceptance_body_rmse_v_max="), std::string::npos);
    EXPECT_NE(report.find("acceptance_body_je_rmse_a_per_m2_max="),
              std::string::npos);
    EXPECT_NE(report.find("acceptance_gate_status="), std::string::npos);
    EXPECT_NE(report.find("acceptance_gate_checks_failed="), std::string::npos);
    EXPECT_NE(report.find("acceptance_gate_failure_details="), std::string::npos);
    EXPECT_NE(report.find("body_negative_tail_threshold_v="), std::string::npos);
    EXPECT_NE(report.find("body_negative_tail_sample_count="), std::string::npos);
    EXPECT_NE(report.find("body_negative_tail_je_rmse_a_per_m2="), std::string::npos);
    EXPECT_NE(report.find("body_negative_tail_jnet_rmse_a_per_m2="), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, BenchmarkReportAppliesAcceptanceThresholdGate)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
    preset.steps = 4;

    ASSERT_TRUE(charging.initialize(preset.config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_acceptance_gate.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::filesystem::path report_path = csv_path;
    report_path.replace_extension(".benchmark.txt");
    const auto report = readKeyValueReport(report_path);
    ASSERT_FALSE(report.empty());

    const double patch_rmse = readReportDoubleOrDefault(report, "patch_rmse_v", 0.0);
    const double patch_rmse_max =
        readReportDoubleOrDefault(report, "acceptance_patch_rmse_v_max", 0.0);
    const double body_rmse = readReportDoubleOrDefault(report, "body_rmse_v", 0.0);
    const double body_rmse_max =
        readReportDoubleOrDefault(report, "acceptance_body_rmse_v_max", 0.0);
    const double body_je_rmse =
        readReportDoubleOrDefault(report, "body_je_rmse_a_per_m2", 0.0);
    const double body_je_rmse_max =
        readReportDoubleOrDefault(report, "acceptance_body_je_rmse_a_per_m2_max", 0.0);
    const double body_jnet_rmse =
        readReportDoubleOrDefault(report, "body_jnet_rmse_a_per_m2", 0.0);
    const double body_jnet_rmse_max =
        readReportDoubleOrDefault(report, "acceptance_body_jnet_rmse_a_per_m2_max", 0.0);
    const double tail_samples =
        readReportDoubleOrDefault(report, "body_negative_tail_sample_count", 0.0);
    const double tail_je_rmse =
        readReportDoubleOrDefault(report, "body_negative_tail_je_rmse_a_per_m2", 0.0);
    const double tail_je_rmse_max = readReportDoubleOrDefault(
        report, "acceptance_negative_tail_body_je_rmse_a_per_m2_max", 0.0);
    const double tail_jnet_rmse =
        readReportDoubleOrDefault(report, "body_negative_tail_jnet_rmse_a_per_m2", 0.0);
    const double tail_jnet_rmse_max = readReportDoubleOrDefault(
        report, "acceptance_negative_tail_body_jnet_rmse_a_per_m2_max", 0.0);
    const double gate_pass =
        readReportDoubleOrDefault(report, "acceptance_gate_pass", 0.0);
    const double gate_checks_failed =
        readReportDoubleOrDefault(report, "acceptance_gate_checks_failed", -1.0);
    std::string gate_status = "";
    const auto gate_status_it = report.find("acceptance_gate_status");
    if (gate_status_it != report.end())
    {
        gate_status = gate_status_it->second;
    }

    EXPECT_LE(patch_rmse, patch_rmse_max)
        << "patch_rmse_v exceeds threshold: rmse=" << patch_rmse
        << ", max=" << patch_rmse_max;
    EXPECT_LE(body_rmse, body_rmse_max)
        << "body_rmse_v exceeds threshold: rmse=" << body_rmse
        << ", max=" << body_rmse_max;
    EXPECT_LE(body_je_rmse, body_je_rmse_max)
        << "body_je_rmse_a_per_m2 exceeds threshold: rmse=" << body_je_rmse
        << ", max=" << body_je_rmse_max;
    EXPECT_LE(body_jnet_rmse, body_jnet_rmse_max)
        << "body_jnet_rmse_a_per_m2 exceeds threshold: rmse=" << body_jnet_rmse
        << ", max=" << body_jnet_rmse_max;

    if (tail_samples > 0.0 && tail_je_rmse_max > 0.0)
    {
        EXPECT_LE(tail_je_rmse, tail_je_rmse_max)
            << "body_negative_tail_je_rmse_a_per_m2 exceeds threshold: rmse="
            << tail_je_rmse << ", max=" << tail_je_rmse_max;
    }
    if (tail_samples > 0.0 && tail_jnet_rmse_max > 0.0)
    {
        EXPECT_LE(tail_jnet_rmse, tail_jnet_rmse_max)
            << "body_negative_tail_jnet_rmse_a_per_m2 exceeds threshold: rmse="
            << tail_jnet_rmse << ", max=" << tail_jnet_rmse_max;
    }

    EXPECT_EQ(gate_status, "PASS");
    EXPECT_EQ(gate_pass, 1.0);
    EXPECT_EQ(gate_checks_failed, 0.0);
}

TEST(SurfaceChargingSmokeTest, LegacyExecuteModeWritesLowErrorGeoReport)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
    preset.steps = 4;

    ASSERT_TRUE(charging.initialize(preset.config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_execute_geo_report.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    auto report_path = csv_path;
    report_path.replace_extension(".benchmark.txt");
    std::ifstream input(report_path);
    ASSERT_TRUE(input.is_open());

    std::string line;
    double rmse_v = std::numeric_limits<double>::quiet_NaN();
    while (std::getline(input, line))
    {
        if (line.rfind("patch_rmse_v=", 0) == 0)
        {
            rmse_v = std::stod(line.substr(std::string("patch_rmse_v=").size()));
            break;
        }
    }

    EXPECT_TRUE(std::isfinite(rmse_v));
    EXPECT_LT(rmse_v, 1.0);
}

TEST(SurfaceChargingSmokeTest, LegacyExecuteModeWritesFiniteLeoReport)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
    preset.steps = 8;

    ASSERT_TRUE(charging.initialize(preset.config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_execute_leo_report.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    auto report_path = csv_path;
    report_path.replace_extension(".benchmark.txt");
    std::ifstream input(report_path);
    ASSERT_TRUE(input.is_open());

    std::string line;
    double patch_rmse_v = std::numeric_limits<double>::quiet_NaN();
    double body_rmse_v = std::numeric_limits<double>::quiet_NaN();
    double body_je_rmse_a_per_m2 = std::numeric_limits<double>::quiet_NaN();
    double body_je_efolding_actual_ev = std::numeric_limits<double>::quiet_NaN();
    double body_je_efolding_reference_ev = std::numeric_limits<double>::quiet_NaN();
    double body_je_initial_ratio = std::numeric_limits<double>::quiet_NaN();
    std::string consistency_status;
    std::string consistency_authority;
    std::string acceptance_gate_status;
    while (std::getline(input, line))
    {
        if (line.rfind("patch_rmse_v=", 0) == 0)
        {
            patch_rmse_v = std::stod(line.substr(std::string("patch_rmse_v=").size()));
        }
        else if (line.rfind("body_rmse_v=", 0) == 0)
        {
            body_rmse_v = std::stod(line.substr(std::string("body_rmse_v=").size()));
        }
        else if (line.rfind("body_je_rmse_a_per_m2=", 0) == 0)
        {
            body_je_rmse_a_per_m2 =
                std::stod(line.substr(std::string("body_je_rmse_a_per_m2=").size()));
        }
        else if (line.rfind("body_je_initial_efolding_energy_ev_actual=", 0) == 0)
        {
            body_je_efolding_actual_ev =
                std::stod(line.substr(std::string("body_je_initial_efolding_energy_ev_actual=").size()));
        }
        else if (line.rfind("body_je_initial_efolding_energy_ev_reference=", 0) == 0)
        {
            body_je_efolding_reference_ev = std::stod(
                line.substr(std::string("body_je_initial_efolding_energy_ev_reference=").size()));
        }
        else if (line.rfind("body_je_initial_ratio_actual_to_reference=", 0) == 0)
        {
            body_je_initial_ratio = std::stod(
                line.substr(std::string("body_je_initial_ratio_actual_to_reference=").size()));
        }
        else if (line.rfind("consistency_status=", 0) == 0)
        {
            consistency_status = line.substr(std::string("consistency_status=").size());
        }
        else if (line.rfind("consistency_authority=", 0) == 0)
        {
            consistency_authority = line.substr(std::string("consistency_authority=").size());
        }
        else if (line.rfind("acceptance_gate_status=", 0) == 0)
        {
            acceptance_gate_status = line.substr(std::string("acceptance_gate_status=").size());
        }
    }

    EXPECT_TRUE(std::isfinite(patch_rmse_v));
    EXPECT_TRUE(std::isfinite(body_rmse_v));
    EXPECT_TRUE(std::isfinite(body_je_rmse_a_per_m2));
    EXPECT_TRUE(std::isfinite(body_je_efolding_actual_ev));
    EXPECT_TRUE(std::isfinite(body_je_efolding_reference_ev));
    EXPECT_TRUE(std::isfinite(body_je_initial_ratio));
    // The legacy LEO execute path is still benchmark-driven and not yet a
    // strict formula-for-formula replay of the historical C body solver. Keep
    // a bounded regression guard while the remaining body-current mismatch is
    // being closed.
    EXPECT_LT(body_rmse_v, 5.5e3);
    EXPECT_LT(body_je_rmse_a_per_m2, 1.0e-6);
    EXPECT_GT(body_je_efolding_actual_ev, 0.0);
    EXPECT_GT(body_je_efolding_reference_ev, 0.0);
    EXPECT_LT(std::abs(body_je_initial_ratio - 1.0), 0.2);
    EXPECT_EQ(acceptance_gate_status, "PASS");
    EXPECT_EQ(consistency_status, "InputReferenceMismatch");
    EXPECT_EQ(consistency_authority, "ReferenceCurve");
}

TEST(SurfaceChargingSmokeTest, LegacyExecuteModeWritesFiniteLeoWakeReport)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_wake_facing_legacy_compatible", preset));
    preset.config.legacy_benchmark_execution_mode =
        SCDAT::Toolkit::SurfaceCharging::LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm;
    preset.steps = 8;

    ASSERT_TRUE(charging.initialize(preset.config));
    for (std::size_t step = 0; step < preset.steps; ++step)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_execute_leo_wake_report.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    auto report_path = csv_path;
    report_path.replace_extension(".benchmark.txt");
    std::ifstream input(report_path);
    ASSERT_TRUE(input.is_open());

    std::string line;
    double patch_rmse_v = std::numeric_limits<double>::quiet_NaN();
    double body_rmse_v = std::numeric_limits<double>::quiet_NaN();
    double body_je_rmse_a_per_m2 = std::numeric_limits<double>::quiet_NaN();
    double body_je_efolding_actual_ev = std::numeric_limits<double>::quiet_NaN();
    double body_je_efolding_reference_ev = std::numeric_limits<double>::quiet_NaN();
    double body_je_initial_ratio = std::numeric_limits<double>::quiet_NaN();
    std::string consistency_status;
    std::string consistency_authority;
    std::string acceptance_gate_status;
    while (std::getline(input, line))
    {
        if (line.rfind("patch_rmse_v=", 0) == 0)
        {
            patch_rmse_v = std::stod(line.substr(std::string("patch_rmse_v=").size()));
        }
        else if (line.rfind("body_rmse_v=", 0) == 0)
        {
            body_rmse_v = std::stod(line.substr(std::string("body_rmse_v=").size()));
        }
        else if (line.rfind("body_je_rmse_a_per_m2=", 0) == 0)
        {
            body_je_rmse_a_per_m2 =
                std::stod(line.substr(std::string("body_je_rmse_a_per_m2=").size()));
        }
        else if (line.rfind("body_je_initial_efolding_energy_ev_actual=", 0) == 0)
        {
            body_je_efolding_actual_ev =
                std::stod(line.substr(std::string("body_je_initial_efolding_energy_ev_actual=").size()));
        }
        else if (line.rfind("body_je_initial_efolding_energy_ev_reference=", 0) == 0)
        {
            body_je_efolding_reference_ev = std::stod(
                line.substr(std::string("body_je_initial_efolding_energy_ev_reference=").size()));
        }
        else if (line.rfind("body_je_initial_ratio_actual_to_reference=", 0) == 0)
        {
            body_je_initial_ratio = std::stod(
                line.substr(std::string("body_je_initial_ratio_actual_to_reference=").size()));
        }
        else if (line.rfind("consistency_status=", 0) == 0)
        {
            consistency_status = line.substr(std::string("consistency_status=").size());
        }
        else if (line.rfind("consistency_authority=", 0) == 0)
        {
            consistency_authority = line.substr(std::string("consistency_authority=").size());
        }
        else if (line.rfind("acceptance_gate_status=", 0) == 0)
        {
            acceptance_gate_status = line.substr(std::string("acceptance_gate_status=").size());
        }
    }

    EXPECT_TRUE(std::isfinite(patch_rmse_v));
    EXPECT_TRUE(std::isfinite(body_rmse_v));
    EXPECT_TRUE(std::isfinite(body_je_rmse_a_per_m2));
    EXPECT_TRUE(std::isfinite(body_je_efolding_actual_ev));
    EXPECT_TRUE(std::isfinite(body_je_efolding_reference_ev));
    EXPECT_TRUE(std::isfinite(body_je_initial_ratio));
    EXPECT_LT(body_rmse_v, 6.5e3);
    EXPECT_LT(body_je_rmse_a_per_m2, 1.0e-6);
    EXPECT_GT(body_je_efolding_actual_ev, 0.0);
    EXPECT_GT(body_je_efolding_reference_ev, 0.0);
    EXPECT_LT(std::abs(body_je_initial_ratio - 1.0), 0.2);
    EXPECT_EQ(acceptance_gate_status, "PASS");
    EXPECT_EQ(consistency_status, "InputReferenceMismatch");
    EXPECT_EQ(consistency_authority, "ReferenceCurve");
}

TEST(SurfaceChargingSmokeTest, UnifiedRouteExportsFiniteFieldAndChargeDiagnostics)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_pic_circuit", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_field_charge.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    const auto header_columns = splitCsvLine(header);
    const auto field_index =
        findColumnIndex(header_columns, "normal_electric_field_v_per_m");
    const auto charge_index =
        findColumnIndex(header_columns, "local_charge_density_c_per_m3");
    const auto cap_index =
        findColumnIndex(header_columns, "capacitance_per_area_f_per_m2");
    ASSERT_LT(field_index, header_columns.size());
    ASSERT_LT(charge_index, header_columns.size());
    ASSERT_LT(cap_index, header_columns.size());

    bool saw_nonzero_field = false;
    bool saw_nonzero_charge = false;
    bool saw_finite_capacitance = false;
    std::string row;
    while (std::getline(input, row))
    {
        if (row.empty())
        {
            continue;
        }
        const auto values = splitCsvLine(row);
        ASSERT_GT(values.size(), cap_index);
        const double field_v_per_m = std::stod(values[field_index]);
        const double charge_c_per_m3 = std::stod(values[charge_index]);
        const double capacitance = std::stod(values[cap_index]);
        EXPECT_TRUE(std::isfinite(field_v_per_m));
        EXPECT_TRUE(std::isfinite(charge_c_per_m3));
        EXPECT_TRUE(std::isfinite(capacitance));
        saw_nonzero_field = saw_nonzero_field || std::abs(field_v_per_m) > 0.0;
        saw_nonzero_charge = saw_nonzero_charge || std::abs(charge_c_per_m3) > 0.0;
        saw_finite_capacitance = saw_finite_capacitance || capacitance > 0.0;
    }

    EXPECT_TRUE(saw_nonzero_field);
    EXPECT_TRUE(saw_nonzero_charge);
    EXPECT_TRUE(saw_finite_capacitance);
}

TEST(SurfaceChargingSmokeTest, LegacyCBenchmarkCurvesCanBeComparedAgainstReferenceFiles)
{
    auto run_and_compare = [](const std::string& preset_name,
                              const std::filesystem::path& csv_path,
                              const std::filesystem::path& ref_path,
                              std::size_t max_steps = std::numeric_limits<std::size_t>::max()) {
        DensePlasmaSurfaceCharging charging;
        SurfaceChargingScenarioPreset preset;
        EXPECT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
            preset_name, preset));
        EXPECT_TRUE(charging.initialize(preset.config));
        const std::size_t steps = std::min(preset.steps, max_steps);
        for (std::size_t step = 0; step < steps; ++step)
        {
            EXPECT_TRUE(charging.advance(preset.time_step_s));
        }
        EXPECT_TRUE(charging.exportResults(csv_path));

        const auto benchmark_curve = readBenchmarkCsvCurve(csv_path);
        const auto reference_curve = readReferenceDatCurve(ref_path);
        EXPECT_GE(benchmark_curve.size(), 8u);
        EXPECT_GE(reference_curve.size(), 8u);

        const double rmse = curveRmse(benchmark_curve, reference_curve);
        const double end_delta =
            benchmark_curve.empty() || reference_curve.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : benchmark_curve.back().potential_v - reference_curve.back().potential_v;

        EXPECT_TRUE(std::isfinite(rmse));
        EXPECT_TRUE(std::isfinite(end_delta));
    };

    run_and_compare("geo_ecss_kapton_ref_legacy_compatible",
                    std::filesystem::temp_directory_path() /
                        "surface_charging_legacy_geo_compare.csv",
                    "ref/GEO/reference_result_patch.dat");
    run_and_compare("leo_ref_ram_facing_legacy_compatible",
                    std::filesystem::temp_directory_path() /
                        "surface_charging_legacy_leo_ram_compare.csv",
                    "ref/LEO/reference_result_patch.dat");
}

TEST(SurfaceChargingSmokeTest, LegacyMatlabBenchmarkCurveCanBeComparedAgainstGeneratedBaseline)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref_matlab_compatible", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    constexpr std::size_t kPrefixSteps = 1024;
    for (std::size_t step = 0; step < kPrefixSteps; ++step)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_legacy_matlab_compare.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    const auto benchmark_curve = readBenchmarkCsvCurve(csv_path);
    const auto reference_curve = readReferenceDatCurve("ref/mat/reference_result_geo_patch.dat");
    ASSERT_GE(benchmark_curve.size(), kPrefixSteps);
    ASSERT_GE(reference_curve.size(), kPrefixSteps);

    double sum = 0.0;
    for (std::size_t i = 0; i < kPrefixSteps; ++i)
    {
        const double delta = benchmark_curve[i].potential_v - reference_curve[i].potential_v;
        sum += delta * delta;
    }

    const double rmse = std::sqrt(sum / static_cast<double>(kPrefixSteps));
    const double end_delta = benchmark_curve[kPrefixSteps - 1].potential_v -
                             reference_curve[kPrefixSteps - 1].potential_v;
    EXPECT_TRUE(std::isfinite(rmse));
    EXPECT_TRUE(std::isfinite(end_delta));
    EXPECT_NEAR(rmse, 0.0, 1.0e-9);
    EXPECT_NEAR(end_delta, 0.0, 1.0e-9);
}

TEST(SurfaceChargingSmokeTest, StructuredBodiesPatchesAndInterfacesExportNodeAndInterfaceSeries)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_pic_calibration = true;
    config.enable_live_pic_window = true;
    config.enable_live_pic_mcc = true;
    config.internal_substeps = 2;
    config.bodies = {
        {"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0},
        {"appendage", 2.0e-2, true, 0.2, 8.0e-11, 0.0, false, 0.0},
    };
    config.patches.clear();
    config.patches.push_back({"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0});
    config.patches.push_back({"ptfe_wake", "appendage", 9.0e-3, 2.2, 0.0});
    config.interfaces = {
        {"body_link", "chassis", "appendage", 6.0e-11, 0.0, 0.0, false},
        {"patch_link", "kapton_ram", "ptfe_wake", 1.0e-11, 0.0, 0.0, false},
    };

    auto ptfe_material = preset.config.material;
    ptfe_material.setName("ptfe_patch");
    ptfe_material.setPermittivity(2.1);
    ptfe_material.setConductivity(1.0e-17);
    ptfe_material.setSecondaryElectronYield(1.4);
    ptfe_material.setScalarProperty("photoelectron_yield", 0.005);
    config.patches[1].material = ptfe_material;
    config.patches[1].patch_incidence_angle_deg = 78.0;
    config.patches[1].electron_collection_coefficient = 1.25;

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_structured_topology.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    EXPECT_NE(header.find("runtime_route_id"), std::string::npos);
    EXPECT_NE(header.find("benchmark_source_id"), std::string::npos);
    EXPECT_NE(header.find("benchmark_execution_mode_id"), std::string::npos);
    EXPECT_NE(header.find("surface_node_0_potential_v"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_total_current_density_a_per_m2"), std::string::npos);
    EXPECT_NE(header.find("surface_node_3_electron_current_density_a_per_m2"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_normal_electric_field_v_per_m"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_propagated_reference_potential_v"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_field_solver_reference_potential_v"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_graph_capacitance_diagonal_f"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_graph_capacitance_row_sum_f"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_field_solver_coupling_gain"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_field_solver_capacitance_scale"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_pseudo_volume_m3"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_volume_projection_weight_sum"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_volume_mesh_coupling_gain"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_electron_pic_calibration_factor"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_live_pic_mcc_enabled"), std::string::npos);
    EXPECT_NE(header.find("surface_interface_0_current_a"), std::string::npos);
    EXPECT_NE(header.find("surface_interface_0_conductance_s"), std::string::npos);
    EXPECT_NE(header.find("surface_interface_0_voltage_drop_v"), std::string::npos);
    EXPECT_NE(header.find("surface_interface_0_power_w"), std::string::npos);
    EXPECT_NE(header.find("surface_interface_0_mutual_capacitance_f"), std::string::npos);
    EXPECT_NE(header.find("surface_interface_1_current_a"), std::string::npos);

    const auto header_columns = splitCsvLine(header);
    const auto node2_index = findColumnIndex(header_columns, "surface_node_2_potential_v");
    const auto node3_index = findColumnIndex(header_columns, "surface_node_3_potential_v");
    const auto node2_total_current_index =
        findColumnIndex(header_columns, "surface_node_2_total_current_density_a_per_m2");
    const auto node3_electron_current_index =
        findColumnIndex(header_columns, "surface_node_3_electron_current_density_a_per_m2");
    const auto node2_field_index =
        findColumnIndex(header_columns, "surface_node_2_normal_electric_field_v_per_m");
    const auto node2_propagated_reference_index =
        findColumnIndex(header_columns, "surface_node_2_propagated_reference_potential_v");
    const auto node2_field_solver_reference_index =
        findColumnIndex(header_columns, "surface_node_2_field_solver_reference_potential_v");
    const auto node2_graph_capacitance_index =
        findColumnIndex(header_columns, "surface_node_2_graph_capacitance_diagonal_f");
    const auto node2_graph_capacitance_row_sum_index =
        findColumnIndex(header_columns, "surface_node_2_graph_capacitance_row_sum_f");
    const auto node2_field_solver_coupling_gain_index =
        findColumnIndex(header_columns, "surface_node_2_field_solver_coupling_gain");
    const auto node2_field_solver_capacitance_scale_index =
        findColumnIndex(header_columns, "surface_node_2_field_solver_capacitance_scale");
    const auto node2_pseudo_volume_index =
        findColumnIndex(header_columns, "surface_node_2_pseudo_volume_m3");
    const auto node2_volume_projection_weight_sum_index =
        findColumnIndex(header_columns, "surface_node_2_volume_projection_weight_sum");
    const auto node2_volume_mesh_coupling_gain_index =
        findColumnIndex(header_columns, "surface_node_2_volume_mesh_coupling_gain");
    const auto node2_volume_potential_index =
        findColumnIndex(header_columns, "surface_node_2_volume_potential_v");
    const auto node2_deposited_charge_index =
        findColumnIndex(header_columns, "surface_node_2_deposited_charge_c");
    const auto node2_poisson_residual_index =
        findColumnIndex(header_columns, "surface_node_2_poisson_residual_v_m");
    const auto node2_electron_calibration_index =
        findColumnIndex(header_columns, "surface_node_2_electron_pic_calibration_factor");
    const auto node2_live_pic_mcc_enabled_index =
        findColumnIndex(header_columns, "surface_node_2_live_pic_mcc_enabled");
    const auto interface0_index =
        findColumnIndex(header_columns, "surface_interface_0_current_a");
    const auto interface0_conductance_index =
        findColumnIndex(header_columns, "surface_interface_0_conductance_s");
    const auto interface0_voltage_drop_index =
        findColumnIndex(header_columns, "surface_interface_0_voltage_drop_v");
    const auto interface0_power_index =
        findColumnIndex(header_columns, "surface_interface_0_power_w");
    const auto interface0_mutual_capacitance_index =
        findColumnIndex(header_columns, "surface_interface_0_mutual_capacitance_f");
    ASSERT_LT(node2_index, header_columns.size());
    ASSERT_LT(node3_index, header_columns.size());
    ASSERT_LT(node2_total_current_index, header_columns.size());
    ASSERT_LT(node3_electron_current_index, header_columns.size());
    ASSERT_LT(node2_field_index, header_columns.size());
    ASSERT_LT(node2_propagated_reference_index, header_columns.size());
    ASSERT_LT(node2_field_solver_reference_index, header_columns.size());
    ASSERT_LT(node2_graph_capacitance_index, header_columns.size());
    ASSERT_LT(node2_graph_capacitance_row_sum_index, header_columns.size());
    ASSERT_LT(node2_field_solver_coupling_gain_index, header_columns.size());
    ASSERT_LT(node2_field_solver_capacitance_scale_index, header_columns.size());
    ASSERT_LT(node2_pseudo_volume_index, header_columns.size());
    ASSERT_LT(node2_volume_projection_weight_sum_index, header_columns.size());
    ASSERT_LT(node2_volume_mesh_coupling_gain_index, header_columns.size());
    ASSERT_LT(node2_volume_potential_index, header_columns.size());
    ASSERT_LT(node2_deposited_charge_index, header_columns.size());
    ASSERT_LT(node2_poisson_residual_index, header_columns.size());
    ASSERT_LT(node2_electron_calibration_index, header_columns.size());
    ASSERT_LT(node2_live_pic_mcc_enabled_index, header_columns.size());
    ASSERT_LT(interface0_index, header_columns.size());
    ASSERT_LT(interface0_conductance_index, header_columns.size());
    ASSERT_LT(interface0_voltage_drop_index, header_columns.size());
    ASSERT_LT(interface0_power_index, header_columns.size());
    ASSERT_LT(interface0_mutual_capacitance_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(std::getline(input, row));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node3_index);
    ASSERT_GT(values.size(), node2_total_current_index);
    ASSERT_GT(values.size(), node3_electron_current_index);
    ASSERT_GT(values.size(), node2_field_index);
    ASSERT_GT(values.size(), node2_propagated_reference_index);
    ASSERT_GT(values.size(), node2_field_solver_reference_index);
    ASSERT_GT(values.size(), node2_graph_capacitance_index);
    ASSERT_GT(values.size(), node2_graph_capacitance_row_sum_index);
    ASSERT_GT(values.size(), node2_field_solver_coupling_gain_index);
    ASSERT_GT(values.size(), node2_field_solver_capacitance_scale_index);
    ASSERT_GT(values.size(), node2_pseudo_volume_index);
    ASSERT_GT(values.size(), node2_volume_projection_weight_sum_index);
    ASSERT_GT(values.size(), node2_volume_mesh_coupling_gain_index);
    ASSERT_GT(values.size(), node2_volume_potential_index);
    ASSERT_GT(values.size(), node2_deposited_charge_index);
    ASSERT_GT(values.size(), node2_poisson_residual_index);
    ASSERT_GT(values.size(), node2_electron_calibration_index);
    ASSERT_GT(values.size(), node2_live_pic_mcc_enabled_index);
    ASSERT_GT(values.size(), interface0_index);
    ASSERT_GT(values.size(), interface0_conductance_index);
    ASSERT_GT(values.size(), interface0_voltage_drop_index);
    ASSERT_GT(values.size(), interface0_power_index);
    ASSERT_GT(values.size(), interface0_mutual_capacitance_index);

    const double patch_a = std::stod(values[node2_index]);
    const double patch_b = std::stod(values[node3_index]);
    const double node2_total_current = std::stod(values[node2_total_current_index]);
    const double node3_electron_current = std::stod(values[node3_electron_current_index]);
    const double node2_field = std::stod(values[node2_field_index]);
    const double node2_propagated_reference = std::stod(values[node2_propagated_reference_index]);
    const double node2_field_solver_reference =
        std::stod(values[node2_field_solver_reference_index]);
    const double node2_graph_capacitance = std::stod(values[node2_graph_capacitance_index]);
    const double node2_graph_capacitance_row_sum =
        std::stod(values[node2_graph_capacitance_row_sum_index]);
    const double node2_field_solver_coupling_gain =
        std::stod(values[node2_field_solver_coupling_gain_index]);
    const double node2_field_solver_capacitance_scale =
        std::stod(values[node2_field_solver_capacitance_scale_index]);
    const double node2_pseudo_volume = std::stod(values[node2_pseudo_volume_index]);
    const double node2_volume_projection_weight_sum =
        std::stod(values[node2_volume_projection_weight_sum_index]);
    const double node2_volume_mesh_coupling_gain =
        std::stod(values[node2_volume_mesh_coupling_gain_index]);
    const double node2_volume_potential = std::stod(values[node2_volume_potential_index]);
    const double node2_deposited_charge = std::stod(values[node2_deposited_charge_index]);
    const double node2_poisson_residual = std::stod(values[node2_poisson_residual_index]);
    const double node2_electron_calibration = std::stod(values[node2_electron_calibration_index]);
    const double node2_live_pic_mcc_enabled = std::stod(values[node2_live_pic_mcc_enabled_index]);
    const double branch_current = std::stod(values[interface0_index]);
    const double branch_conductance = std::stod(values[interface0_conductance_index]);
    const double branch_voltage_drop = std::stod(values[interface0_voltage_drop_index]);
    const double branch_power = std::stod(values[interface0_power_index]);
    const double branch_mutual_capacitance = std::stod(values[interface0_mutual_capacitance_index]);
    EXPECT_TRUE(std::isfinite(patch_a));
    EXPECT_TRUE(std::isfinite(patch_b));
    EXPECT_TRUE(std::isfinite(node2_total_current));
    EXPECT_TRUE(std::isfinite(node3_electron_current));
    EXPECT_TRUE(std::isfinite(node2_field));
    EXPECT_TRUE(std::isfinite(node2_propagated_reference));
    EXPECT_TRUE(std::isfinite(node2_field_solver_reference));
    EXPECT_TRUE(std::isfinite(node2_graph_capacitance));
    EXPECT_TRUE(std::isfinite(node2_graph_capacitance_row_sum));
    EXPECT_TRUE(std::isfinite(node2_field_solver_coupling_gain));
    EXPECT_TRUE(std::isfinite(node2_field_solver_capacitance_scale));
    EXPECT_TRUE(std::isfinite(node2_pseudo_volume));
    EXPECT_TRUE(std::isfinite(node2_volume_projection_weight_sum));
    EXPECT_TRUE(std::isfinite(node2_volume_mesh_coupling_gain));
    EXPECT_TRUE(std::isfinite(node2_volume_potential));
    EXPECT_TRUE(std::isfinite(node2_deposited_charge));
    EXPECT_TRUE(std::isfinite(node2_poisson_residual));
    EXPECT_TRUE(std::isfinite(node2_electron_calibration));
    EXPECT_TRUE(std::isfinite(node2_live_pic_mcc_enabled));
    EXPECT_TRUE(std::isfinite(branch_current));
    EXPECT_TRUE(std::isfinite(branch_conductance));
    EXPECT_TRUE(std::isfinite(branch_voltage_drop));
    EXPECT_TRUE(std::isfinite(branch_power));
    EXPECT_TRUE(std::isfinite(branch_mutual_capacitance));
    EXPECT_GT(node2_pseudo_volume, 0.0);
    EXPECT_GE(node2_volume_projection_weight_sum, 1.0);
    EXPECT_GE(node2_volume_mesh_coupling_gain, 0.0);
    EXPECT_GT(std::abs(node2_deposited_charge), 0.0);
    EXPECT_LT(std::abs(node2_poisson_residual), 1.0);
    EXPECT_GT(std::abs(node2_volume_potential - patch_a), 1.0e-6);
    EXPECT_GT(std::abs(patch_a - patch_b), 1.0e-6);
}

TEST(SurfaceChargingSmokeTest, StructuredTopologyRejectsPatchWithoutValidBody)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.bodies = {{"body_a", 1.0e-2, true, 0.0, 1.0e-10, 0.0, false, 0.0}};
    config.patches = {{"patch_a", "missing_body", 5.0e-3, 0.0, 0.0}};
    config.interfaces.clear();
    config.surface_nodes.clear();
    config.surface_branches.clear();

    EXPECT_FALSE(charging.initialize(config));
    EXPECT_NE(charging.lastErrorMessage().find("unknown body"), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, StructuredTopologyRejectsBoundaryMappingWithUnknownGroup)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.bodies = {{"body_a", 1.0e-2, true, 0.0, 1.0e-10, 0.0, false, 0.0}};
    config.patches = {{"patch_a", "body_a", 5.0e-3, 0.0, 0.0}};
    config.boundary_mappings = {{"patch_a", "missing_group", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    EXPECT_FALSE(charging.initialize(config));
    EXPECT_NE(charging.lastErrorMessage().find("unknown boundary group"), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, StructuredTopologyWritesReadableGraphSidecarReport)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_pic_calibration = true;
    config.enable_live_pic_window = true;
    config.enable_live_pic_mcc = true;
    config.internal_substeps = 2;
    config.bodies = {
        {"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0},
        {"appendage", 2.0e-2, true, 0.2, 8.0e-11, 0.0, false, 0.0},
    };
    config.patches.clear();
    config.patches.push_back({"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0});
    config.patches.push_back({"ptfe_wake", "appendage", 9.0e-3, 2.2, 0.0});
    config.interfaces = {
        {"body_link", "chassis", "appendage", 6.0e-11, 0.0, 0.0, false},
        {"patch_link", "kapton_ram", "ptfe_wake", 1.0e-11, 0.0, 0.0, false},
    };

    auto ptfe_material = preset.config.material;
    ptfe_material.setName("ptfe_patch");
    ptfe_material.setPermittivity(2.1);
    ptfe_material.setConductivity(1.0e-17);
    ptfe_material.setSecondaryElectronYield(1.4);
    ptfe_material.setScalarProperty("photoelectron_yield", 0.005);
    config.patches[1].material = ptfe_material;
    config.patches[1].patch_incidence_angle_deg = 78.0;
    config.patches[1].electron_collection_coefficient = 1.25;

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_structured_topology_graph.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    auto report_path = csv_path;
    report_path.replace_extension(".graph.txt");
    const auto report = readKeyValueReport(report_path);
    ASSERT_FALSE(report.empty());

    EXPECT_EQ(report.at("runtime_route"), "SCDATUnified");
    EXPECT_EQ(report.at("surface_reference_graph_propagator"),
              "BranchWeightedReferenceGraphPropagator");
    EXPECT_EQ(report.at("surface_graph_capacitance_matrix_provider"),
              "InterfaceAwareGraphCapacitanceMatrixProvider");
    EXPECT_EQ(report.at("surface_field_solver_adapter"), "GraphCoupledFieldSolverAdapter");
    EXPECT_EQ(report.at("surface_graph_capacitance_matrix_family"),
              "InterfaceAwareDenseHeuristic");
    EXPECT_EQ(report.at("node_count"), "4");
    EXPECT_EQ(report.at("body_count"), "2");
    EXPECT_EQ(report.at("patch_count"), "2");
    EXPECT_EQ(report.at("body_body_interface_count"), "1");
    EXPECT_EQ(report.at("body_patch_interface_count"), "2");
    EXPECT_EQ(report.at("patch_patch_interface_count"), "1");
    EXPECT_TRUE(std::stod(report.at("max_abs_neighbor_potential_delta_v")) > 0.0);
    EXPECT_TRUE(std::stod(report.at("max_abs_neighbor_field_contrast_v_per_m")) >= 0.0);
    EXPECT_TRUE(std::stod(report.at("max_abs_node_field_solver_reference_offset_v")) >= 0.0);
    EXPECT_TRUE(std::stod(report.at("max_field_solver_coupling_gain")) >= 0.0);
    EXPECT_TRUE(std::stod(report.at("graph_coupling_metric")) >= 0.0);
    EXPECT_FALSE(report.at("strongest_field_node_name").empty());
    EXPECT_NE(report.at("strongest_coupling_interface_name").find("->"), std::string::npos);
    EXPECT_NE(report.at("node_2").find("material:kapton"), std::string::npos);
    EXPECT_NE(report.at("node_2").find("last_reference_offset_v:"), std::string::npos);
    EXPECT_NE(report.at("node_2").find("last_field_solver_reference_offset_v:"), std::string::npos);
    EXPECT_NE(report.at("node_2").find("last_graph_capacitance_row_sum_f:"), std::string::npos);
    EXPECT_NE(report.at("node_2").find("last_field_solver_capacitance_scale:"), std::string::npos);
    EXPECT_NE(report.at("node_3").find("material:ptfe_patch"), std::string::npos);
    EXPECT_NE(report.at("interface_0").find("type:body-patch"), std::string::npos);

    auto monitor_json_path = csv_path;
    monitor_json_path.replace_extension(".monitor.json");
    std::ifstream monitor_json_input(monitor_json_path);
    ASSERT_TRUE(monitor_json_input.is_open());
    std::string monitor_json_text((std::istreambuf_iterator<char>(monitor_json_input)),
                                  std::istreambuf_iterator<char>());
    EXPECT_NE(monitor_json_text.find("\"schema_version\": \"scdat.surface_monitor.v1\""),
              std::string::npos);
    EXPECT_NE(monitor_json_text.find("\"graph_coupling_metric\""), std::string::npos);
    EXPECT_NE(monitor_json_text.find("\"reference_offset_envelope_v\""), std::string::npos);
    EXPECT_NE(monitor_json_text.find("\"consistency_status\""), std::string::npos);
    EXPECT_NE(monitor_json_text.find("\"benchmark_mode\": \"UnifiedSpisAligned\""),
              std::string::npos);

    auto h5_path = csv_path;
    h5_path.replace_extension(".h5");
    std::ifstream h5_input(h5_path);
    ASSERT_TRUE(h5_input.is_open());
    std::string h5_text((std::istreambuf_iterator<char>(h5_input)),
                        std::istreambuf_iterator<char>());
    EXPECT_NE(h5_text.find("surface_graph_matrix_runtime_weight="), std::string::npos);
    EXPECT_NE(h5_text.find("surface_graph_matrix_runtime_coupling_enabled="), std::string::npos);

    auto matrix_path = csv_path;
    matrix_path.replace_extension(".graph_matrix.csv");
    std::ifstream matrix_input(matrix_path);
    ASSERT_TRUE(matrix_input.is_open());
    std::string matrix_header;
    std::getline(matrix_input, matrix_header);
    EXPECT_NE(matrix_header.find("entry_type"), std::string::npos);
    EXPECT_NE(matrix_header.find("value_f"), std::string::npos);
    std::string matrix_row;
    ASSERT_TRUE(std::getline(matrix_input, matrix_row));
    EXPECT_NE(matrix_row.find("node_diagonal"), std::string::npos);

    auto matrix_json_path = csv_path;
    matrix_json_path.replace_extension(".graph_matrix.json");
    std::ifstream matrix_json_input(matrix_json_path);
    ASSERT_TRUE(matrix_json_input.is_open());
    std::string matrix_json_text((std::istreambuf_iterator<char>(matrix_json_input)),
                                 std::istreambuf_iterator<char>());
    EXPECT_NE(matrix_json_text.find("\"schema_version\": \"scdat.surface_graph_matrix.v1\""),
              std::string::npos);
    EXPECT_NE(matrix_json_text.find("\"interfaces\""), std::string::npos);

    auto field_adapter_path = csv_path;
    field_adapter_path.replace_extension(".field_adapter.txt");
    const auto field_adapter_report = readKeyValueReport(field_adapter_path);
    ASSERT_FALSE(field_adapter_report.empty());
    EXPECT_EQ(field_adapter_report.at("surface_field_solver_adapter"),
              "GraphCoupledFieldSolverAdapter");
    EXPECT_EQ(field_adapter_report.at("surface_graph_capacitance_matrix_family"),
              "InterfaceAwareDenseHeuristic");
    EXPECT_EQ(field_adapter_report.at("supports_mutual_matrix"), "1");
    EXPECT_EQ(field_adapter_report.at("supports_diagonal_matrix"), "1");
    EXPECT_EQ(field_adapter_report.at("node_count"), "4");
    EXPECT_EQ(field_adapter_report.at("branch_count"), "4");

    auto field_adapter_json_path = csv_path;
    field_adapter_json_path.replace_extension(".field_adapter.json");
    std::ifstream field_adapter_json_input(field_adapter_json_path);
    ASSERT_TRUE(field_adapter_json_input.is_open());
    std::string field_adapter_json_text((std::istreambuf_iterator<char>(field_adapter_json_input)),
                                        std::istreambuf_iterator<char>());
    EXPECT_NE(field_adapter_json_text.find("\"schema_version\": \"scdat.field_solver_adapter_contract.v1\""),
              std::string::npos);
    EXPECT_NE(field_adapter_json_text.find("\"supports_mutual_matrix\": true"), std::string::npos);

    auto boundary_mapping_path = csv_path;
    boundary_mapping_path.replace_extension(".boundary_mapping.json");
    std::ifstream boundary_mapping_input(boundary_mapping_path);
    ASSERT_TRUE(boundary_mapping_input.is_open());
    std::string boundary_mapping_text((std::istreambuf_iterator<char>(boundary_mapping_input)),
                                      std::istreambuf_iterator<char>());
    EXPECT_NE(boundary_mapping_text.find("\"boundary_mappings\""), std::string::npos);

    auto field_request_path = csv_path;
    field_request_path.replace_extension(".field_request.json");
    std::ifstream field_request_input(field_request_path);
    ASSERT_TRUE(field_request_input.is_open());
    std::string field_request_text((std::istreambuf_iterator<char>(field_request_input)),
                                   std::istreambuf_iterator<char>());
    EXPECT_NE(field_request_text.find("\"schema_version\": \"scdat.external_field_request.v1\""),
              std::string::npos);
    EXPECT_NE(field_request_text.find("\"nodes\""), std::string::npos);
    EXPECT_NE(field_request_text.find("\"interfaces\""), std::string::npos);
    EXPECT_NE(field_request_text.find("\"body_boundary_groups\""), std::string::npos);
    EXPECT_NE(field_request_text.find("\"patch_boundary_groups\""), std::string::npos);

    auto field_result_template_path = csv_path;
    field_result_template_path.replace_extension(".field_result_template.json");
    std::ifstream field_result_template_input(field_result_template_path);
    ASSERT_TRUE(field_result_template_input.is_open());
    std::string field_result_template_text(
        (std::istreambuf_iterator<char>(field_result_template_input)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(field_result_template_text.find("\"schema_version\": \"scdat.external_field_result.v1\""),
              std::string::npos);
    EXPECT_NE(field_result_template_text.find("\"capacitance_scale\""), std::string::npos);

    auto volume_stub_path = csv_path;
    volume_stub_path.replace_extension(".volume_stub.json");
    std::ifstream volume_stub_input(volume_stub_path);
    ASSERT_TRUE(volume_stub_input.is_open());
    std::string volume_stub_text((std::istreambuf_iterator<char>(volume_stub_input)),
                                 std::istreambuf_iterator<char>());
    EXPECT_NE(volume_stub_text.find("\"schema_version\": \"scdat.surface_volume_stub.v1\""),
              std::string::npos);
    EXPECT_NE(volume_stub_text.find("\"pseudo_volume_m3\""), std::string::npos);
    EXPECT_NE(volume_stub_text.find("\"volume_potential_v\""), std::string::npos);
    EXPECT_NE(volume_stub_text.find("\"deposited_charge_c\""), std::string::npos);
    EXPECT_NE(volume_stub_text.find("\"poisson_residual_v_m\""), std::string::npos);
    EXPECT_NE(volume_stub_text.find("\"field_solver_capacitance_scale\""), std::string::npos);

    auto volume_mesh_stub_path = csv_path;
    volume_mesh_stub_path.replace_extension(".volume_mesh_stub.json");
    std::ifstream volume_mesh_stub_input(volume_mesh_stub_path);
    ASSERT_TRUE(volume_mesh_stub_input.is_open());
    std::string volume_mesh_stub_text((std::istreambuf_iterator<char>(volume_mesh_stub_input)),
                                      std::istreambuf_iterator<char>());
    EXPECT_NE(volume_mesh_stub_text.find("\"schema_version\": \"scdat.volume_mesh_stub.v1\""),
              std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"boundary_faces\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"mesh_family\": \"PseudoBoundaryVolumeMesh\""),
              std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"center_x_m\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"characteristic_length_m\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"cell_neighbors\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"face_neighbors\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"shared_face_area_m2\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"face_distance_m\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"surface_normal_x\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"face_center_x_m\""), std::string::npos);
    EXPECT_NE(volume_mesh_stub_text.find("\"cell_faces\""), std::string::npos);

    auto field_bridge_manifest_path = csv_path;
    field_bridge_manifest_path.replace_extension(".field_bridge_manifest.json");
    std::ifstream field_bridge_manifest_input(field_bridge_manifest_path);
    ASSERT_TRUE(field_bridge_manifest_input.is_open());
    std::string field_bridge_manifest_text(
        (std::istreambuf_iterator<char>(field_bridge_manifest_input)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(field_bridge_manifest_text.find("\"schema_version\": \"scdat.field_bridge_manifest.v1\""),
              std::string::npos);
    EXPECT_NE(field_bridge_manifest_text.find("\"artifacts\""), std::string::npos);
    EXPECT_NE(field_bridge_manifest_text.find("\"volume_mesh_stub_json\""), std::string::npos);

    auto surface_volume_projection_path = csv_path;
    surface_volume_projection_path.replace_extension(".surface_volume_projection.json");
    std::ifstream surface_volume_projection_input(surface_volume_projection_path);
    ASSERT_TRUE(surface_volume_projection_input.is_open());
    std::string surface_volume_projection_text(
        (std::istreambuf_iterator<char>(surface_volume_projection_input)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(surface_volume_projection_text.find("\"schema_version\": \"scdat.surface_volume_projection.v1\""),
              std::string::npos);
    EXPECT_NE(surface_volume_projection_text.find("\"mesh_family\": \"PseudoBoundaryVolumeMesh\""),
              std::string::npos);
    EXPECT_NE(surface_volume_projection_text.find("\"supports_projection_weights\": true"),
              std::string::npos);
    EXPECT_NE(surface_volume_projection_text.find("\"projection_weights\""), std::string::npos);
    EXPECT_NE(surface_volume_projection_text.find("\"surface_to_volume\""), std::string::npos);
    EXPECT_NE(surface_volume_projection_text.find("\"volume_to_surface\""), std::string::npos);

    auto volumetric_adapter_path = csv_path;
    volumetric_adapter_path.replace_extension(".volumetric_adapter.json");
    std::ifstream volumetric_adapter_input(volumetric_adapter_path);
    ASSERT_TRUE(volumetric_adapter_input.is_open());
    std::string volumetric_adapter_text((std::istreambuf_iterator<char>(volumetric_adapter_input)),
                                        std::istreambuf_iterator<char>());
    EXPECT_NE(volumetric_adapter_text.find("\"schema_version\": \"scdat.volumetric_solver_adapter_contract.v1\""),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"projection_family\": \"NodeToPseudoCellProjection\""),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"mesh_family\": \"PseudoBoundaryVolumeMesh\""),
              std::string::npos);

    auto volume_history_path = csv_path;
    volume_history_path.replace_extension(".volume_history.json");
    std::ifstream volume_history_input(volume_history_path);
    ASSERT_TRUE(volume_history_input.is_open());
    std::string volume_history_text((std::istreambuf_iterator<char>(volume_history_input)),
                                    std::istreambuf_iterator<char>());
    EXPECT_NE(volume_history_text.find("\"schema_version\": \"scdat.volume_history.v1\""),
              std::string::npos);
    EXPECT_NE(volume_history_text.find("\"volume_potential_v\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"deposited_charge_c\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"poisson_residual_v_m\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"solver_mode_id\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"solver_iterations\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"solver_residual_norm\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"field_volume_coupling_relaxation_used\""),
              std::string::npos);
    EXPECT_NE(volume_history_text.find("\"external_volume_feedback_blend_factor\""),
              std::string::npos);
    EXPECT_NE(volume_history_text.find("\"external_volume_feedback_mismatch_metric\""),
              std::string::npos);
    EXPECT_NE(volume_history_text.find("\"external_volume_feedback_applied\""),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_boundary_face_mapping\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_projection_weights\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_cell_centers\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_face_centers\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_face_normals\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_neighbor_face_geometry\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_cell_face_links\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_cell_volume_override\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"supports_cell_initial_state\": true"),
              std::string::npos);
    EXPECT_NE(volumetric_adapter_text.find("\"solver_mode_hint\": \"iterative_or_dense_auto\""),
              std::string::npos);

    auto volume_request_path = csv_path;
    volume_request_path.replace_extension(".volume_request.json");
    std::ifstream volume_request_input(volume_request_path);
    ASSERT_TRUE(volume_request_input.is_open());
    std::string volume_request_text((std::istreambuf_iterator<char>(volume_request_input)),
                                    std::istreambuf_iterator<char>());
    EXPECT_NE(volume_request_text.find("\"schema_version\": \"scdat.external_volume_request.v1\""),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"volume_stub_json\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"volume_mesh_stub_json\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"expects_cell_face_links\": true"),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"expects_cell_volume_override\": true"),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"expects_cell_initial_state\": true"),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"solver_mode_hint\": \"iterative_or_dense_auto\""),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"linear_max_iterations\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"linear_tolerance_scale\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"linear_relaxation\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"self_consistent_iterations\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"self_consistent_tolerance_v\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"charge_relaxation\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"potential_relaxation\""), std::string::npos);

    auto volume_result_template_path = csv_path;
    volume_result_template_path.replace_extension(".volume_result_template.json");
    std::ifstream volume_result_template_input(volume_result_template_path);
    ASSERT_TRUE(volume_result_template_input.is_open());
    std::string volume_result_template_text(
        (std::istreambuf_iterator<char>(volume_result_template_input)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(volume_result_template_text.find("\"schema_version\": \"scdat.external_volume_result.v1\""),
              std::string::npos);
    EXPECT_NE(volume_result_template_text.find("\"cells\""), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, ExternalFieldSolverBridgeIngestsNodeResult)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_external_field_solver_bridge = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_field_solver_result_path =
        temp_dir / "surface_charging_external_bridge_result.json";

    {
        std::ofstream output(config.external_field_solver_result_path);
        ASSERT_TRUE(output.is_open());
        output << "{\n"
               << "  \"nodes\": [\n"
               << "    {\"node_id\": \"patch:kapton_ram\", \"reference_potential_v\": 7.5, "
                  "\"normal_field_v_per_m\": 123.0, \"local_charge_density_c_per_m3\": 4.5e-8, "
                  "\"capacitance_scale\": 1.7}\n"
               << "  ]\n"
               << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_external_bridge.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto node1_reference_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_reference_potential_v");
    const auto node1_scale_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_capacitance_scale");
    ASSERT_LT(node1_reference_index, header_columns.size());
    ASSERT_LT(node1_scale_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node1_reference_index);
    ASSERT_GT(values.size(), node1_scale_index);
    EXPECT_NEAR(std::stod(values[node1_reference_index]), 7.5, 1.0e-9);
    EXPECT_GE(std::stod(values[node1_scale_index]), 1.7);

    auto field_adapter_path = csv_path;
    field_adapter_path.replace_extension(".field_adapter.txt");
    const auto field_adapter_report = readKeyValueReport(field_adapter_path);
    ASSERT_FALSE(field_adapter_report.empty());
    EXPECT_EQ(field_adapter_report.at("surface_field_solver_adapter"),
              "ExternalFileFieldSolverAdapter");
}

TEST(SurfaceChargingSmokeTest, ExternalFieldSolverBridgeIngestsBoundaryGroupResult)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_external_field_solver_bridge = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_field_solver_result_path =
        temp_dir / "surface_charging_external_bridge_boundary_group_result.json";

    {
        std::ofstream output(config.external_field_solver_result_path);
        ASSERT_TRUE(output.is_open());
        output << "{\n"
               << "  \"schema_version\": \"scdat.external_field_result.v1\",\n"
               << "  \"nodes\": [\n"
               << "    {\"boundary_group_id\": \"bg_patch\", \"capacitance_scale\": 1.65, "
                  "\"local_charge_density_c_per_m3\": 3.9e-8, "
                  "\"normal_field_v_per_m\": 222.0, \"reference_potential_v\": 5.25}\n"
               << "  ]\n"
               << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path = std::filesystem::temp_directory_path() /
                          "surface_charging_external_bridge_boundary_group.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto node1_reference_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_reference_potential_v");
    const auto node1_scale_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_capacitance_scale");
    ASSERT_LT(node1_reference_index, header_columns.size());
    ASSERT_LT(node1_scale_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node1_reference_index);
    ASSERT_GT(values.size(), node1_scale_index);
    EXPECT_NEAR(std::stod(values[node1_reference_index]), 5.25, 1.0e-9);
    EXPECT_GE(std::stod(values[node1_scale_index]), 1.65);
}

TEST(SurfaceChargingSmokeTest, ExternalFieldSolverBridgeRejectsSchemaMismatch)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_external_field_solver_bridge = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_field_solver_result_path =
        temp_dir / "surface_charging_external_bridge_bad_schema_result.json";

    {
        std::ofstream output(config.external_field_solver_result_path);
        ASSERT_TRUE(output.is_open());
        output << "{\n"
               << "  \"schema_version\": \"scdat.external_field_result.v999\",\n"
               << "  \"nodes\": [\n"
               << "    {\"node_id\": \"patch:kapton_ram\", \"reference_potential_v\": 9999.0, "
                  "\"normal_field_v_per_m\": 8888.0, \"local_charge_density_c_per_m3\": 1.2e-3, "
                  "\"capacitance_scale\": 3.7}\n"
               << "  ]\n"
               << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path = std::filesystem::temp_directory_path() /
                          "surface_charging_external_bridge_bad_schema.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto node1_reference_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_reference_potential_v");
    const auto node1_scale_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_capacitance_scale");
    ASSERT_LT(node1_reference_index, header_columns.size());
    ASSERT_LT(node1_scale_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node1_reference_index);
    ASSERT_GT(values.size(), node1_scale_index);
    const double node1_reference_v = std::stod(values[node1_reference_index]);
    const double node1_scale = std::stod(values[node1_scale_index]);
    EXPECT_TRUE(std::isfinite(node1_reference_v));
    EXPECT_TRUE(std::isfinite(node1_scale));
    EXPECT_LT(node1_reference_v, 1000.0);
    EXPECT_GT(std::abs(node1_reference_v - 9999.0), 1.0);
    EXPECT_GT(std::abs(node1_scale - 3.7), 1.0e-3);
}

TEST(SurfaceChargingSmokeTest, ExternalVolumeSolverBridgeIngestsCellResult)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_external_volume_solver_bridge = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_volume_solver_result_path =
        temp_dir / "surface_charging_external_volume_bridge_result.json";

    {
        std::ofstream output(config.external_volume_solver_result_path);
        ASSERT_TRUE(output.is_open());
        output << "{\n"
               << "  \"cells\": [\n"
               << "    {\"cell_id\": \"cell_1\", \"potential_v\": 6.25, "
                  "\"reference_potential_v\": 6.75, "
                  "\"normal_field_v_per_m\": 456.0, \"local_charge_density_c_per_m3\": 8.0e-8, "
                  "\"capacitance_scale\": 1.9, \"coupling_gain\": 0.8, "
                  "\"projection_weight_scale\": 1.5, \"sheath_length_scale\": 1.1}\n"
               << "  ]\n"
               << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_external_volume_bridge.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto node1_reference_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_reference_potential_v");
    const auto node1_scale_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_capacitance_scale");
    const auto node1_projection_weight_index =
        findColumnIndex(header_columns, "surface_node_1_volume_projection_weight_sum");
    const auto node1_mesh_coupling_gain_index =
        findColumnIndex(header_columns, "surface_node_1_volume_mesh_coupling_gain");
    const auto node1_external_blend_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_blend_factor");
    const auto node1_external_mismatch_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_mismatch_metric");
    const auto node1_external_applied_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_applied");
    ASSERT_LT(node1_reference_index, header_columns.size());
    ASSERT_LT(node1_scale_index, header_columns.size());
    ASSERT_LT(node1_projection_weight_index, header_columns.size());
    ASSERT_LT(node1_mesh_coupling_gain_index, header_columns.size());
    ASSERT_LT(node1_external_blend_index, header_columns.size());
    ASSERT_LT(node1_external_mismatch_index, header_columns.size());
    ASSERT_LT(node1_external_applied_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node1_reference_index);
    ASSERT_GT(values.size(), node1_scale_index);
    ASSERT_GT(values.size(), node1_projection_weight_index);
    ASSERT_GT(values.size(), node1_mesh_coupling_gain_index);
    ASSERT_GT(values.size(), node1_external_blend_index);
    ASSERT_GT(values.size(), node1_external_mismatch_index);
    ASSERT_GT(values.size(), node1_external_applied_index);
    const double node1_reference_v = std::stod(values[node1_reference_index]);
    const double node1_capacitance_scale = std::stod(values[node1_scale_index]);
    const double node1_projection_weight = std::stod(values[node1_projection_weight_index]);
    const double node1_mesh_coupling = std::stod(values[node1_mesh_coupling_gain_index]);
    const double node1_external_blend = std::stod(values[node1_external_blend_index]);
    const double node1_external_mismatch = std::stod(values[node1_external_mismatch_index]);
    const double node1_external_applied = std::stod(values[node1_external_applied_index]);
    EXPECT_TRUE(std::isfinite(node1_reference_v));
    EXPECT_TRUE(std::isfinite(node1_capacitance_scale));
    EXPECT_TRUE(std::isfinite(node1_projection_weight));
    EXPECT_TRUE(std::isfinite(node1_mesh_coupling));
    EXPECT_TRUE(std::isfinite(node1_external_blend));
    EXPECT_TRUE(std::isfinite(node1_external_mismatch));
    EXPECT_TRUE(std::isfinite(node1_external_applied));
    EXPECT_GT(node1_reference_v, 0.0);
    EXPECT_LE(node1_reference_v, 6.75);
    EXPECT_GE(node1_capacitance_scale, 1.0);
    EXPECT_GE(node1_projection_weight, 1.0);
    EXPECT_GE(node1_mesh_coupling, 0.0);
    EXPECT_GE(node1_external_blend, 0.05);
    EXPECT_LE(node1_external_blend, 1.0);
    EXPECT_GE(node1_external_mismatch, 0.0);
    EXPECT_GE(node1_external_applied, 0.5);

    auto volumetric_adapter_path = csv_path;
    volumetric_adapter_path.replace_extension(".volumetric_adapter.json");
    std::ifstream volumetric_adapter_input(volumetric_adapter_path);
    ASSERT_TRUE(volumetric_adapter_input.is_open());
    std::string volumetric_adapter_text((std::istreambuf_iterator<char>(volumetric_adapter_input)),
                                        std::istreambuf_iterator<char>());
    EXPECT_NE(volumetric_adapter_text.find("\"volumetric_solver_adapter\": \"VolumeStubVolumetricSolverAdapter\""),
              std::string::npos);

    auto volume_result_template_path = csv_path;
    volume_result_template_path.replace_extension(".volume_result_template.json");
    std::ifstream volume_result_template_input(volume_result_template_path);
    ASSERT_TRUE(volume_result_template_input.is_open());
    std::string volume_result_template_text(
        (std::istreambuf_iterator<char>(volume_result_template_input)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(volume_result_template_text.find("\"capacitance_scale\": 1.000000000000"),
              std::string::npos);
    EXPECT_NE(volume_result_template_text.find("\"reference_potential_v\": 0.000000000000"),
              std::string::npos);
    EXPECT_NE(volume_result_template_text.find("\"projection_weight_scale\": 1.000000000000"),
              std::string::npos);
}

TEST(SurfaceChargingSmokeTest, ExternalVolumeSolverBridgeIngestsBoundaryGroupResult)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_external_volume_solver_bridge = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_volume_solver_result_path =
        temp_dir / "surface_charging_external_volume_bridge_bg_result.json";

    {
        std::ofstream output(config.external_volume_solver_result_path);
        ASSERT_TRUE(output.is_open());
        output << "{\n"
               << "  \"schema_version\": \"scdat.external_volume_result.v1\",\n"
               << "  \"cells\": [\n"
               << "    {\"boundary_group_id\": \"bg_patch\", \"potential_v\": 5.75, "
                  "\"reference_potential_v\": 6.1, \"normal_field_v_per_m\": 350.0, "
                  "\"local_charge_density_c_per_m3\": 6.5e-8, \"capacitance_scale\": 1.8, "
                  "\"coupling_gain\": 0.6, \"projection_weight_scale\": 1.3, "
                  "\"sheath_length_scale\": 1.05}\n"
               << "  ]\n"
               << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path = std::filesystem::temp_directory_path() /
                          "surface_charging_external_volume_bridge_bg.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto node1_reference_index =
        findColumnIndex(header_columns, "surface_node_1_field_solver_reference_potential_v");
    const auto node1_external_blend_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_blend_factor");
    const auto node1_external_applied_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_applied");
    ASSERT_LT(node1_reference_index, header_columns.size());
    ASSERT_LT(node1_external_blend_index, header_columns.size());
    ASSERT_LT(node1_external_applied_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node1_reference_index);
    ASSERT_GT(values.size(), node1_external_blend_index);
    ASSERT_GT(values.size(), node1_external_applied_index);
    const double node1_reference_v = std::stod(values[node1_reference_index]);
    const double node1_external_blend = std::stod(values[node1_external_blend_index]);
    const double node1_external_applied = std::stod(values[node1_external_applied_index]);
    EXPECT_TRUE(std::isfinite(node1_reference_v));
    EXPECT_GT(node1_reference_v, 0.0);
    EXPECT_LE(node1_reference_v, 6.1);
    EXPECT_GE(node1_external_blend, 0.05);
    EXPECT_GE(node1_external_applied, 0.5);
}

TEST(SurfaceChargingSmokeTest, ExternalVolumeSolverBridgeRejectsSchemaMismatch)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_external_volume_solver_bridge = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_volume_solver_result_path =
        temp_dir / "surface_charging_external_volume_bridge_bad_schema_result.json";

    {
        std::ofstream output(config.external_volume_solver_result_path);
        ASSERT_TRUE(output.is_open());
        output << "{\n"
               << "  \"schema_version\": \"scdat.external_volume_result.v999\",\n"
               << "  \"cells\": [\n"
               << "    {\"cell_id\": \"cell_1\", \"potential_v\": 9999.0, "
                  "\"reference_potential_v\": 9999.0, \"normal_field_v_per_m\": 9999.0, "
                  "\"local_charge_density_c_per_m3\": 9.9e-4, \"capacitance_scale\": 4.0, "
                  "\"coupling_gain\": 1.0, \"projection_weight_scale\": 4.0, "
                  "\"sheath_length_scale\": 2.0}\n"
               << "  ]\n"
               << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path = std::filesystem::temp_directory_path() /
                          "surface_charging_external_volume_bridge_bad_schema.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto node1_external_blend_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_blend_factor");
    const auto node1_external_applied_index =
        findColumnIndex(header_columns, "surface_node_1_external_volume_feedback_applied");
    ASSERT_LT(node1_external_blend_index, header_columns.size());
    ASSERT_LT(node1_external_applied_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node1_external_blend_index);
    ASSERT_GT(values.size(), node1_external_applied_index);
    EXPECT_NEAR(std::stod(values[node1_external_blend_index]), 0.0, 1.0e-12);
    EXPECT_NEAR(std::stod(values[node1_external_applied_index]), 0.0, 1.0e-12);
}

TEST(SurfaceChargingSmokeTest, VolumetricRuntimeUsesExternalMeshAndProjectionInputs)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_volume_mesh_path =
        temp_dir / "surface_charging_external_mesh_runtime.volume_mesh.json";
    config.external_surface_volume_projection_path =
        temp_dir / "surface_charging_external_mesh_runtime.projection.json";
    config.volume_linear_solver_policy =
        SCDAT::Toolkit::SurfaceCharging::VolumeLinearSolverPolicy::IterativeOnly;

    {
        std::ofstream mesh_output(config.external_volume_mesh_path);
        ASSERT_TRUE(mesh_output.is_open());
        mesh_output << "{\n"
                    << "  \"cells\": [\n"
                    << "    {\"cell_id\": \"ext_body_cell\", \"node_id\": \"body:chassis\", "
                       "\"boundary_group_id\": \"bg_body\", \"node_area_m2\": 0.05, "
                       "\"cell_volume_m3\": 2.5e-5, "
                       "\"initial_potential_v\": 0.1, "
                       "\"initial_charge_density_c_per_m3\": 1.0e-8, "
                       "\"characteristic_length_m\": 0.11, \"center_x_m\": -0.04, "
                       "\"center_y_m\": 0.00, \"center_z_m\": 0.01},\n"
                    << "    {\"cell_id\": \"ext_patch_cell\", \"node_id\": \"patch:kapton_ram\", "
                       "\"boundary_group_id\": \"bg_patch\", \"node_area_m2\": 0.011, "
                       "\"cell_volume_m3\": 1.1e-5, "
                       "\"initial_potential_v\": 0.4, "
                       "\"initial_charge_density_c_per_m3\": 2.0e-8, "
                       "\"characteristic_length_m\": 0.07, \"center_x_m\": 0.21, "
                       "\"center_y_m\": 0.03, \"center_z_m\": 0.06}\n"
                    << "  ],\n"
                    << "  \"cell_neighbors\": [\n"
                    << "    {\"source_cell_id\": \"ext_patch_cell\", "
                       "\"target_cell_id\": \"ext_body_cell\", \"conductance_s\": 3.0e-10, "
                       "\"shared_face_area_m2\": 0.0075, \"face_distance_m\": 0.24, "
                       "\"permittivity_scale\": 1.3, "
                       "\"face_center_x_m\": 0.08, \"face_center_y_m\": 0.01, \"face_center_z_m\": 0.02}\n"
                    << "  ]\n"
                    << "}\n";
    }

    {
        std::ofstream projection_output(config.external_surface_volume_projection_path);
        ASSERT_TRUE(projection_output.is_open());
        projection_output << "{\n"
                          << "  \"projection_weights\": [\n"
                          << "    {\"node_id\": \"body:chassis\", \"cell_id\": \"ext_body_cell\", "
                             "\"surface_to_volume_weight\": 1.4, \"volume_to_surface_weight\": 1.2},\n"
                          << "    {\"node_id\": \"patch:kapton_ram\", \"cell_id\": \"ext_patch_cell\", "
                             "\"surface_to_volume_weight\": 3.5, \"volume_to_surface_weight\": 2.5}\n"
                          << "  ]\n"
                          << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_external_mesh_runtime.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto projection_weight_index =
        findColumnIndex(header_columns, "surface_node_1_volume_projection_weight_sum");
    const auto pseudo_volume_index =
        findColumnIndex(header_columns, "surface_node_1_pseudo_volume_m3");
    const auto volume_potential_index =
        findColumnIndex(header_columns, "surface_node_1_volume_potential_v");
    const auto solver_mode_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_mode_id");
    const auto solver_iterations_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_iterations");
    const auto solver_linear_iterations_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_linear_iterations");
    const auto solver_converged_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_converged");
    const auto solver_residual_norm_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_residual_norm");
    const auto solver_matrix_nnz_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_matrix_nnz");
    const auto solver_cell_count_index =
        findColumnIndex(header_columns, "surface_node_1_volume_solver_cell_count");
    const auto coupling_iterations_index =
        findColumnIndex(header_columns, "surface_node_1_field_volume_coupling_iterations");
    const auto coupling_converged_index =
        findColumnIndex(header_columns, "surface_node_1_field_volume_coupling_converged");
    const auto coupling_max_delta_index =
        findColumnIndex(header_columns, "surface_node_1_field_volume_coupling_max_delta");
    const auto coupling_relaxation_used_index =
        findColumnIndex(header_columns, "surface_node_1_field_volume_coupling_relaxation_used");
    ASSERT_LT(projection_weight_index, header_columns.size());
    ASSERT_LT(pseudo_volume_index, header_columns.size());
    ASSERT_LT(volume_potential_index, header_columns.size());
    ASSERT_LT(solver_mode_index, header_columns.size());
    ASSERT_LT(solver_iterations_index, header_columns.size());
    ASSERT_LT(solver_linear_iterations_index, header_columns.size());
    ASSERT_LT(solver_converged_index, header_columns.size());
    ASSERT_LT(solver_residual_norm_index, header_columns.size());
    ASSERT_LT(solver_matrix_nnz_index, header_columns.size());
    ASSERT_LT(solver_cell_count_index, header_columns.size());
    ASSERT_LT(coupling_iterations_index, header_columns.size());
    ASSERT_LT(coupling_converged_index, header_columns.size());
    ASSERT_LT(coupling_max_delta_index, header_columns.size());
    ASSERT_LT(coupling_relaxation_used_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), projection_weight_index);
    ASSERT_GT(values.size(), pseudo_volume_index);
    ASSERT_GT(values.size(), volume_potential_index);
    ASSERT_GT(values.size(), solver_mode_index);
    ASSERT_GT(values.size(), solver_iterations_index);
    ASSERT_GT(values.size(), solver_linear_iterations_index);
    ASSERT_GT(values.size(), solver_converged_index);
    ASSERT_GT(values.size(), solver_residual_norm_index);
    ASSERT_GT(values.size(), solver_matrix_nnz_index);
    ASSERT_GT(values.size(), solver_cell_count_index);
    ASSERT_GT(values.size(), coupling_iterations_index);
    ASSERT_GT(values.size(), coupling_converged_index);
    ASSERT_GT(values.size(), coupling_max_delta_index);
    ASSERT_GT(values.size(), coupling_relaxation_used_index);
    EXPECT_GE(std::stod(values[projection_weight_index]), 3.5);
    EXPECT_GT(std::stod(values[pseudo_volume_index]), 0.0);
    EXPECT_TRUE(std::isfinite(std::stod(values[volume_potential_index])));
    EXPECT_GE(std::stod(values[solver_mode_index]), 1.5);
    EXPECT_GE(std::stod(values[solver_iterations_index]), 1.0);
    EXPECT_GE(std::stod(values[solver_linear_iterations_index]), 1.0);
    EXPECT_GE(std::stod(values[solver_converged_index]), 0.0);
    EXPECT_TRUE(std::isfinite(std::stod(values[solver_residual_norm_index])));
    EXPECT_GE(std::stod(values[solver_matrix_nnz_index]), 1.0);
    EXPECT_GE(std::stod(values[solver_cell_count_index]), 2.0);
    EXPECT_GE(std::stod(values[coupling_iterations_index]), 1.0);
    EXPECT_GE(std::stod(values[coupling_converged_index]), 0.0);
    EXPECT_TRUE(std::isfinite(std::stod(values[coupling_max_delta_index])));
    EXPECT_TRUE(std::isfinite(std::stod(values[coupling_relaxation_used_index])));
    EXPECT_GE(std::stod(values[coupling_relaxation_used_index]), 0.05);
    EXPECT_LE(std::stod(values[coupling_relaxation_used_index]), 1.0);

    auto volume_request_path = csv_path;
    volume_request_path.replace_extension(".volume_request.json");
    std::ifstream volume_request_input(volume_request_path);
    ASSERT_TRUE(volume_request_input.is_open());
    std::string volume_request_text((std::istreambuf_iterator<char>(volume_request_input)),
                                    std::istreambuf_iterator<char>());
    EXPECT_NE(volume_request_text.find("\"solver_mode_hint\": \"iterative_only\""),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"field_volume_outer_iterations\""),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"field_volume_outer_relaxation\""),
              std::string::npos);
    EXPECT_NE(volume_request_text.find("\"linear_max_iterations\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"linear_tolerance_scale\""), std::string::npos);
    EXPECT_NE(volume_request_text.find("\"linear_relaxation\""), std::string::npos);

    auto volume_history_path = csv_path;
    volume_history_path.replace_extension(".volume_history.json");
    std::ifstream volume_history_input(volume_history_path);
    ASSERT_TRUE(volume_history_input.is_open());
    std::string volume_history_text((std::istreambuf_iterator<char>(volume_history_input)),
                                    std::istreambuf_iterator<char>());
    EXPECT_NE(volume_history_text.find("\"solver_matrix_nnz\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"solver_cell_count\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"field_volume_coupling_iterations\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"field_volume_coupling_max_delta\""), std::string::npos);
    EXPECT_NE(volume_history_text.find("\"field_volume_coupling_relaxation_used\""),
              std::string::npos);
    EXPECT_NE(volume_history_text.find("\"external_volume_feedback_blend_factor\""),
              std::string::npos);
    EXPECT_NE(volume_history_text.find("\"external_volume_feedback_applied\""),
              std::string::npos);
}

TEST(SurfaceChargingSmokeTest, VolumetricRuntimeDerivesGeometryFromExternalBoundaryFaces)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11, 12}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_volume_mesh_path =
        temp_dir / "surface_charging_external_face_geometry.volume_mesh.json";
    config.external_surface_volume_projection_path =
        temp_dir / "surface_charging_external_face_geometry.projection.json";

    {
        std::ofstream mesh_output(config.external_volume_mesh_path);
        ASSERT_TRUE(mesh_output.is_open());
        mesh_output << "{\n"
                    << "  \"boundary_faces\": [\n"
                    << "    {\"face_id\": \"f_body\", \"boundary_group_id\": \"bg_body\", "
                       "\"node_id\": \"body:chassis\", \"area_m2\": 0.015, "
                       "\"center_x_m\": -0.02, \"center_y_m\": 0.0, \"center_z_m\": 0.0, "
                       "\"normal_x\": 0.0, \"normal_y\": 0.0, \"normal_z\": 1.0},\n"
                    << "    {\"face_id\": \"f_patch_0\", \"boundary_group_id\": \"bg_patch\", "
                       "\"node_id\": \"patch:kapton_ram\", \"area_m2\": 0.010, "
                       "\"center_x_m\": 0.21, \"center_y_m\": 0.00, \"center_z_m\": 0.02, "
                       "\"normal_x\": 1.0, \"normal_y\": 0.0, \"normal_z\": 0.0},\n"
                    << "    {\"face_id\": \"f_patch_1\", \"boundary_group_id\": \"bg_patch\", "
                       "\"node_id\": \"patch:kapton_ram\", \"area_m2\": 0.009, "
                       "\"center_x_m\": 0.22, \"center_y_m\": 0.01, \"center_z_m\": 0.03, "
                       "\"normal_x\": 0.9, \"normal_y\": 0.1, \"normal_z\": 0.1}\n"
                    << "  ],\n"
                    << "  \"cells\": [\n"
                    << "    {\"cell_id\": \"ext_body_cell\", \"node_id\": \"body:chassis\", "
                       "\"boundary_group_id\": \"bg_body\"},\n"
                    << "    {\"cell_id\": \"ext_patch_cell\", \"node_id\": \"patch:kapton_ram\", "
                       "\"boundary_group_id\": \"bg_patch\"}\n"
                    << "  ],\n"
                    << "  \"cell_neighbors\": [\n"
                    << "    {\"source_cell_id\": \"ext_patch_cell\", "
                       "\"target_cell_id\": \"ext_body_cell\", \"conductance_s\": 2.5e-10}\n"
                    << "  ]\n"
                    << "}\n";
    }

    {
        std::ofstream projection_output(config.external_surface_volume_projection_path);
        ASSERT_TRUE(projection_output.is_open());
        projection_output << "{\n"
                          << "  \"projection_weights\": [\n"
                          << "    {\"node_id\": \"body:chassis\", \"cell_id\": \"ext_body_cell\", "
                             "\"surface_to_volume_weight\": 1.2, \"volume_to_surface_weight\": 1.1},\n"
                          << "    {\"node_id\": \"patch:kapton_ram\", \"cell_id\": \"ext_patch_cell\", "
                             "\"surface_to_volume_weight\": 2.0, \"volume_to_surface_weight\": 1.8}\n"
                          << "  ]\n"
                          << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_external_face_geometry.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto pseudo_volume_index =
        findColumnIndex(header_columns, "surface_node_1_pseudo_volume_m3");
    ASSERT_LT(pseudo_volume_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), pseudo_volume_index);
    const double pseudo_volume = std::stod(values[pseudo_volume_index]);
    EXPECT_GT(pseudo_volume, 3.0e-5);
}

TEST(SurfaceChargingSmokeTest, VolumetricRuntimeRespondsToExternalCellCenterGeometry)
{
    auto run_case = [](const std::filesystem::path& mesh_path,
                       const std::filesystem::path& projection_path,
                       const std::filesystem::path& csv_path,
                       double patch_center_x_m) {
        DensePlasmaSurfaceCharging charging;
        SurfaceChargingScenarioPreset preset;
        EXPECT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
            "leo_ref_ram_facing", preset));

        auto config = preset.config;
        config.enable_body_patch_circuit = true;
        config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
        config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
        config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1}}};
        config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10}}};
        config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
        config.surface_nodes.clear();
        config.surface_branches.clear();
        config.external_volume_mesh_path = mesh_path;
        config.external_surface_volume_projection_path = projection_path;

        {
            std::ofstream mesh_output(config.external_volume_mesh_path);
            EXPECT_TRUE(mesh_output.is_open());
            mesh_output << "{\n"
                        << "  \"cells\": [\n"
                        << "    {\"cell_id\": \"ext_body_cell\", \"node_id\": \"body:chassis\", "
                           "\"boundary_group_id\": \"bg_body\", \"node_area_m2\": 0.05, "
                           "\"characteristic_length_m\": 0.09, \"center_x_m\": 0.0, "
                           "\"center_y_m\": 0.0, \"center_z_m\": 0.0},\n"
                        << "    {\"cell_id\": \"ext_patch_cell\", \"node_id\": \"patch:kapton_ram\", "
                           "\"boundary_group_id\": \"bg_patch\", \"node_area_m2\": 0.011, "
                           "\"characteristic_length_m\": 0.05, \"center_x_m\": " << patch_center_x_m
                        << ", \"center_y_m\": 0.0, \"center_z_m\": 0.02}\n"
                        << "  ],\n"
                        << "  \"cell_neighbors\": [\n"
                        << "    {\"source_cell_id\": \"ext_patch_cell\", "
                           "\"target_cell_id\": \"ext_body_cell\", \"conductance_s\": 2.0e-10}\n"
                        << "  ]\n"
                        << "}\n";
        }

        {
            std::ofstream projection_output(config.external_surface_volume_projection_path);
            EXPECT_TRUE(projection_output.is_open());
            projection_output << "{\n"
                              << "  \"projection_weights\": [\n"
                              << "    {\"node_id\": \"body:chassis\", \"cell_id\": \"ext_body_cell\", "
                                 "\"surface_to_volume_weight\": 1.0, \"volume_to_surface_weight\": 1.0},\n"
                              << "    {\"node_id\": \"patch:kapton_ram\", \"cell_id\": \"ext_patch_cell\", "
                                 "\"surface_to_volume_weight\": 2.0, \"volume_to_surface_weight\": 2.0}\n"
                              << "  ]\n"
                              << "}\n";
        }

        EXPECT_TRUE(charging.initialize(config));
        EXPECT_TRUE(charging.advance(1.0e-5));
        EXPECT_TRUE(charging.exportResults(csv_path));

        std::ifstream input(csv_path);
        EXPECT_TRUE(input.is_open());
        std::string header;
        EXPECT_TRUE(static_cast<bool>(std::getline(input, header)));
        const auto header_columns = splitCsvLine(header);
        const auto volume_potential_index =
            findColumnIndex(header_columns, "surface_node_1_volume_potential_v");
        const auto coupling_gain_index =
            findColumnIndex(header_columns, "surface_node_1_volume_mesh_coupling_gain");
        EXPECT_LT(volume_potential_index, header_columns.size());
        EXPECT_LT(coupling_gain_index, header_columns.size());

        std::string row;
        EXPECT_TRUE(static_cast<bool>(std::getline(input, row)));
        const auto values = splitCsvLine(row);
        EXPECT_GT(values.size(), volume_potential_index);
        EXPECT_GT(values.size(), coupling_gain_index);
        return std::make_pair(std::stod(values[volume_potential_index]),
                              std::stod(values[coupling_gain_index]));
    };

    const auto temp_dir = std::filesystem::temp_directory_path();
    const auto near_mesh = temp_dir / "surface_charging_external_center_near.volume_mesh.json";
    const auto far_mesh = temp_dir / "surface_charging_external_center_far.volume_mesh.json";
    const auto projection = temp_dir / "surface_charging_external_center.projection.json";
    const auto near_csv = temp_dir / "surface_charging_external_center_near.csv";
    const auto far_csv = temp_dir / "surface_charging_external_center_far.csv";

    const auto near_result = run_case(near_mesh, projection, near_csv, 0.08);
    const auto far_result = run_case(far_mesh, projection, far_csv, 0.45);

    EXPECT_TRUE(std::isfinite(near_result.first));
    EXPECT_TRUE(std::isfinite(far_result.first));
    EXPECT_TRUE(std::isfinite(near_result.second));
    EXPECT_TRUE(std::isfinite(far_result.second));
    const double combined_difference =
        std::abs(near_result.first - far_result.first) +
        std::abs(near_result.second - far_result.second);
    EXPECT_GT(combined_difference, 1.0e-6);
}

TEST(SurfaceChargingSmokeTest, VolumetricRuntimeDerivesCellGeometryFromExternalCellFaces)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.bodies = {{"chassis", 5.0e-2, false, 0.0, config.body_capacitance_f, 0.0, true, 0.0}};
    config.patches = {{"kapton_ram", "chassis", 1.1e-2, 0.8, 0.0}};
    config.body_boundary_groups = {{"bg_body", "chassis", "body_group", {1, 2}}};
    config.patch_boundary_groups = {{"bg_patch", "kapton_ram", "patch_group", {10, 11}}};
    config.boundary_mappings = {{"patch:kapton_ram", "bg_patch", true}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    const auto temp_dir = std::filesystem::temp_directory_path();
    config.external_volume_mesh_path =
        temp_dir / "surface_charging_external_cell_faces.volume_mesh.json";
    config.external_surface_volume_projection_path =
        temp_dir / "surface_charging_external_cell_faces.projection.json";

    {
        std::ofstream mesh_output(config.external_volume_mesh_path);
        ASSERT_TRUE(mesh_output.is_open());
        mesh_output << "{\n"
                    << "  \"boundary_faces\": [\n"
                    << "    {\"face_id\": \"body_face_0\", \"boundary_group_id\": \"bg_body\", "
                       "\"node_id\": \"body:chassis\", \"area_m2\": 0.020, "
                       "\"center_x_m\": -0.03, \"center_y_m\": 0.00, \"center_z_m\": 0.00, "
                       "\"normal_x\": 0.0, \"normal_y\": 0.0, \"normal_z\": 1.0},\n"
                    << "    {\"face_id\": \"body_face_1\", \"boundary_group_id\": \"bg_body\", "
                       "\"node_id\": \"body:chassis\", \"area_m2\": 0.018, "
                       "\"center_x_m\": -0.02, \"center_y_m\": 0.01, \"center_z_m\": 0.00, "
                       "\"normal_x\": 0.0, \"normal_y\": 0.0, \"normal_z\": 1.0},\n"
                    << "    {\"face_id\": \"patch_face_0\", \"boundary_group_id\": \"bg_patch\", "
                       "\"node_id\": \"patch:kapton_ram\", \"area_m2\": 0.010, "
                       "\"center_x_m\": 0.18, \"center_y_m\": 0.00, \"center_z_m\": 0.02, "
                       "\"normal_x\": 1.0, \"normal_y\": 0.0, \"normal_z\": 0.0},\n"
                    << "    {\"face_id\": \"patch_face_1\", \"boundary_group_id\": \"bg_patch\", "
                       "\"node_id\": \"patch:kapton_ram\", \"area_m2\": 0.009, "
                       "\"center_x_m\": 0.20, \"center_y_m\": 0.02, \"center_z_m\": 0.04, "
                       "\"normal_x\": 0.9, \"normal_y\": 0.2, \"normal_z\": 0.1}\n"
                    << "  ],\n"
                    << "  \"cells\": [\n"
                    << "    {\"cell_id\": \"ext_body_cell\", \"node_id\": \"body:chassis\", "
                       "\"boundary_group_id\": \"bg_body\"},\n"
                    << "    {\"cell_id\": \"ext_patch_cell\", \"node_id\": \"patch:kapton_ram\", "
                       "\"boundary_group_id\": \"bg_patch\"}\n"
                    << "  ],\n"
                    << "  \"cell_faces\": [\n"
                    << "    {\"cell_id\": \"ext_body_cell\", \"face_id\": \"body_face_0\", "
                       "\"role\": \"boundary\", \"projection_weight\": 1.0},\n"
                    << "    {\"cell_id\": \"ext_body_cell\", \"face_id\": \"body_face_1\", "
                       "\"role\": \"boundary\", \"projection_weight\": 1.0},\n"
                    << "    {\"cell_id\": \"ext_patch_cell\", \"face_id\": \"patch_face_0\", "
                       "\"role\": \"boundary\", \"projection_weight\": 2.0},\n"
                    << "    {\"cell_id\": \"ext_patch_cell\", \"face_id\": \"patch_face_1\", "
                       "\"role\": \"boundary\", \"projection_weight\": 2.0}\n"
                    << "  ],\n"
                    << "  \"cell_neighbors\": [\n"
                    << "    {\"source_cell_id\": \"ext_patch_cell\", \"target_cell_id\": "
                       "\"ext_body_cell\", \"conductance_s\": 2.5e-10}\n"
                    << "  ]\n"
                    << "}\n";
    }

    {
        std::ofstream projection_output(config.external_surface_volume_projection_path);
        ASSERT_TRUE(projection_output.is_open());
        projection_output << "{\n"
                          << "  \"projection_weights\": [\n"
                          << "    {\"node_id\": \"body:chassis\", \"cell_id\": \"ext_body_cell\", "
                             "\"surface_to_volume_weight\": 1.3, \"volume_to_surface_weight\": 1.1},\n"
                          << "    {\"node_id\": \"patch:kapton_ram\", \"cell_id\": \"ext_patch_cell\", "
                             "\"surface_to_volume_weight\": 2.4, \"volume_to_surface_weight\": 2.2}\n"
                          << "  ]\n"
                          << "}\n";
    }

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_external_cell_faces.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, header)));
    const auto header_columns = splitCsvLine(header);
    const auto pseudo_volume_index =
        findColumnIndex(header_columns, "surface_node_1_pseudo_volume_m3");
    const auto volume_potential_index =
        findColumnIndex(header_columns, "surface_node_1_volume_potential_v");
    ASSERT_LT(pseudo_volume_index, header_columns.size());
    ASSERT_LT(volume_potential_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, row)));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), pseudo_volume_index);
    ASSERT_GT(values.size(), volume_potential_index);
    EXPECT_GT(std::stod(values[pseudo_volume_index]), 0.0);
    EXPECT_TRUE(std::isfinite(std::stod(values[volume_potential_index])));
}

TEST(SurfaceChargingSmokeTest, StructuredTopologyRejectsInterfaceWithUnknownEndpoint)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.bodies = {{"body_a", 1.0e-2, true, 0.0, 1.0e-10, 0.0, false, 0.0}};
    config.patches = {{"patch_a", "body_a", 5.0e-3, 0.0, 0.0}};
    config.interfaces = {{"bad_link", "patch_a", "ghost_patch", 1.0e-9, 0.0, 0.0, false}};
    config.surface_nodes.clear();
    config.surface_branches.clear();

    EXPECT_FALSE(charging.initialize(config));
    EXPECT_NE(charging.lastErrorMessage().find("unknown to_id"), std::string::npos);
}

TEST(SurfaceChargingSmokeTest, StructuredTopologyRejectsDuplicateNodeIds)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.bodies = {
        {"body_a", 1.0e-2, true, 0.0, 1.0e-10, 0.0, false, 0.0},
        {"body_a", 2.0e-2, true, 0.0, 2.0e-10, 0.0, false, 0.0},
    };
    config.patches.clear();
    config.interfaces.clear();
    config.surface_nodes.clear();
    config.surface_branches.clear();

    EXPECT_FALSE(charging.initialize(config));
    EXPECT_NE(charging.lastErrorMessage().find("Duplicate structured topology body id"),
              std::string::npos);
}

TEST(SurfaceChargingSmokeTest, MultiPatchCircuitExportsIndependentNodePotentials)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_pic_calibration = false;
    config.enable_live_pic_window = false;
    config.enable_live_pic_mcc = false;
    config.internal_substeps = 2;
    config.body_floating = false;
    config.body_initial_potential_v = 0.0;
    config.surface_nodes = {
        {"body", 3.0e-2, false, 0.0, config.body_capacitance_f, true, 0.0},
        {"patch_a", 1.0e-2, true, 0.5, 0.0, false, 0.0},
        {"patch_b", 1.5e-2, true, 2.5, 0.0, false, 0.0},
    };
    config.surface_branches = {
        {1, 0, 8.0e-11, 0.0, 0.0},
        {2, 0, 4.0e-11, 0.0, 0.0},
        {1, 2, 1.0e-11, 0.0, 0.0},
    };

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_multi_patch.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    EXPECT_NE(header.find("surface_node_0_potential_v"), std::string::npos);
    EXPECT_NE(header.find("surface_node_1_potential_v"), std::string::npos);
    EXPECT_NE(header.find("surface_node_2_potential_v"), std::string::npos);

    const auto header_columns = splitCsvLine(header);
    const std::size_t node1_index = findColumnIndex(header_columns, "surface_node_1_potential_v");
    const std::size_t node2_index = findColumnIndex(header_columns, "surface_node_2_potential_v");
    ASSERT_LT(node1_index, header_columns.size());
    ASSERT_LT(node2_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(std::getline(input, row));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node2_index);

    const double patch_a = std::stod(values[node1_index]);
    const double patch_b = std::stod(values[node2_index]);
    EXPECT_TRUE(std::isfinite(patch_a));
    EXPECT_TRUE(std::isfinite(patch_b));
    EXPECT_NEAR(charging.getStatus().body_potential_v, 0.0, 1.0e-12);
    EXPECT_GT(std::abs(patch_a - patch_b), 1.0e-6);
}

TEST(SurfaceChargingSmokeTest, MultiPatchPhysicsOverridesDriveDistinctPatchPotentials)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    auto config = preset.config;
    config.enable_body_patch_circuit = true;
    config.enable_pic_calibration = false;
    config.enable_live_pic_window = false;
    config.enable_live_pic_mcc = false;
    config.internal_substeps = 2;
    config.body_floating = false;
    config.body_initial_potential_v = 0.0;
    config.surface_nodes = {
        {"body", 4.0e-2, false, 0.0, config.body_capacitance_f, true, 0.0},
        {"patch_a", 1.0e-2, true, 1.0, 0.0, false, 0.0},
        {"patch_b", 1.0e-2, true, 1.0, 0.0, false, 0.0},
    };
    config.surface_branches = {
        {1, 0, 8.0e-11, 0.0, 0.0},
        {2, 0, 8.0e-11, 0.0, 0.0},
    };

    auto patch_b_spectrum = config.electron_spectrum;
    for (auto& population : patch_b_spectrum.populations)
    {
        population.density_m3 *= 2.0;
        population.temperature_ev *= 1.5;
    }

    SCDAT::Material::MaterialProperty patch_b_material(
        9, SCDAT::Mesh::MaterialType::DIELECTRIC, "ptfe_override");
    patch_b_material.setPermittivity(2.1);
    patch_b_material.setConductivity(1.0e-17);
    patch_b_material.setWorkFunctionEv(5.75);
    patch_b_material.setSecondaryElectronYield(1.4);
    patch_b_material.setScalarProperty("photoelectron_yield", 0.004);
    patch_b_material.setScalarProperty("photoelectron_escape_energy_ev", 1.2);
    patch_b_material.setScalarProperty("secondary_emission_escape_energy_ev", 4.0);
    patch_b_material.setScalarProperty("poole_frenkel_beta", 7.5e-3);
    patch_b_material.setScalarProperty("max_field_enhancement_factor", 5.0e8);

    SCDAT::Toolkit::SurfaceCharging::SurfacePatchPhysicsConfig patch_b_override;
    patch_b_override.node_name = "patch_b";
    patch_b_override.override_material = true;
    patch_b_override.material = patch_b_material;
    patch_b_override.override_patch_incidence_angle = true;
    patch_b_override.patch_incidence_angle_deg = 82.0;
    patch_b_override.override_electron_spectrum = true;
    patch_b_override.electron_spectrum = patch_b_spectrum;
    patch_b_override.has_electron_spectrum = true;
    patch_b_override.override_patch_photo_current_density = true;
    patch_b_override.patch_photo_current_density_a_per_m2 = 0.0;
    patch_b_override.override_electron_collection_coefficient = true;
    patch_b_override.electron_collection_coefficient = 1.3;
    config.patch_physics_overrides.push_back(patch_b_override);

    ASSERT_TRUE(charging.initialize(config));
    ASSERT_TRUE(charging.advance(1.0e-5));

    const auto csv_path =
        std::filesystem::temp_directory_path() / "surface_charging_multi_patch_override.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    const auto header_columns = splitCsvLine(header);
    const std::size_t node1_index = findColumnIndex(header_columns, "surface_node_1_potential_v");
    const std::size_t node2_index = findColumnIndex(header_columns, "surface_node_2_potential_v");
    ASSERT_LT(node1_index, header_columns.size());
    ASSERT_LT(node2_index, header_columns.size());

    std::string row;
    ASSERT_TRUE(std::getline(input, row));
    const auto values = splitCsvLine(row);
    ASSERT_GT(values.size(), node2_index);
    const double patch_a = std::stod(values[node1_index]);
    const double patch_b = std::stod(values[node2_index]);

    EXPECT_TRUE(std::isfinite(patch_a));
    EXPECT_TRUE(std::isfinite(patch_b));
    EXPECT_GT(std::abs(patch_a - patch_b), 1.0e-4);
    EXPECT_GT(patch_a, patch_b);
}

TEST(SurfaceChargingSmokeTest, PositiveBiasSuppressesPhotoelectronEscape)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_pic_circuit_ram_facing", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const auto near_zero = charging.computeSurfaceCurrents(0.0);
    const auto positive_bias = charging.computeSurfaceCurrents(2.0);

    EXPECT_LT(positive_bias.photo_emission_a_per_m2, near_zero.photo_emission_a_per_m2);
}

TEST(SurfaceChargingSmokeTest, GeoPicCalibrationProducesFiniteTrajectory)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_pic_circuit", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    ASSERT_TRUE(charging.advance(preset.time_step_s));
    ASSERT_TRUE(charging.advance(preset.time_step_s));

    const auto csv_path = std::filesystem::temp_directory_path() / "surface_charging_geo_pic.csv";
    ASSERT_TRUE(charging.exportResults(csv_path));

    std::ifstream input(csv_path);
    ASSERT_TRUE(input.is_open());
    std::string header;
    std::getline(input, header);
    EXPECT_NE(header.find("ion_pic_calibration_factor"), std::string::npos);
    EXPECT_NE(header.find("pic_recalibration_trigger"), std::string::npos);
    EXPECT_NE(header.find("live_pic_electron_collection_density_a_per_m2"), std::string::npos);
    EXPECT_NE(header.find("live_pic_collection_derivative_a_per_m2_per_v"), std::string::npos);
    EXPECT_NE(header.find("live_pic_mcc_enabled"), std::string::npos);

    const auto header_columns = splitCsvLine(header);
    const std::size_t live_pic_index =
        findColumnIndex(header_columns, "live_pic_electron_collection_density_a_per_m2");
    const std::size_t live_pic_mcc_enabled_index =
        findColumnIndex(header_columns, "live_pic_mcc_enabled");
    ASSERT_LT(live_pic_index, header_columns.size());
    ASSERT_LT(live_pic_mcc_enabled_index, header_columns.size());

    bool found_nonzero_live_pic = false;
    bool found_mcc_enabled = false;
    std::string row;
    while (std::getline(input, row))
    {
        if (row.empty())
        {
            continue;
        }
        const auto row_values = splitCsvLine(row);
        if (row_values.size() <= live_pic_index)
        {
            continue;
        }
        if (std::abs(std::stod(row_values[live_pic_index])) > 0.0)
        {
            found_nonzero_live_pic = true;
        }
        if (std::stod(row_values[live_pic_mcc_enabled_index]) > 0.5)
        {
            found_mcc_enabled = true;
        }
        if (found_nonzero_live_pic && found_mcc_enabled)
        {
            break;
        }
    }

    EXPECT_TRUE(found_nonzero_live_pic);
    EXPECT_TRUE(found_mcc_enabled);
    EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
}

TEST(SurfaceChargingSmokeTest, PicCircuitPresetsKeepLivePicMccEnabled)
{
    SurfaceChargingScenarioPreset geo_preset;
    SurfaceChargingScenarioPreset leo_preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_pic_circuit", geo_preset));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_pic_circuit_ram_facing", leo_preset));

    EXPECT_TRUE(geo_preset.config.enable_live_pic_window);
    EXPECT_TRUE(geo_preset.config.enable_live_pic_mcc);
    EXPECT_TRUE(geo_preset.config.enable_pic_calibration);

    EXPECT_TRUE(leo_preset.config.enable_live_pic_window);
    EXPECT_TRUE(leo_preset.config.enable_live_pic_mcc);
    EXPECT_TRUE(leo_preset.config.enable_pic_calibration);
}

TEST(SurfaceChargingSmokeTest, LeoPicCircuitWakeRemainsAboveRamFacing)
{
    DensePlasmaSurfaceCharging ram_charging;
    DensePlasmaSurfaceCharging wake_charging;
    SurfaceChargingScenarioPreset ram_preset;
    SurfaceChargingScenarioPreset wake_preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_pic_circuit_ram_facing", ram_preset));
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_pic_circuit_wake_facing", wake_preset));

    ASSERT_TRUE(ram_charging.initialize(ram_preset.config));
    ASSERT_TRUE(wake_charging.initialize(wake_preset.config));
    for (int i = 0; i < 4; ++i)
    {
        ASSERT_TRUE(ram_charging.advance(ram_preset.time_step_s));
        ASSERT_TRUE(wake_charging.advance(wake_preset.time_step_s));
    }

    EXPECT_GT(wake_charging.getStatus().state.surface_potential_v,
              ram_charging.getStatus().state.surface_potential_v);
    EXPECT_TRUE(ram_preset.config.enable_live_pic_mcc);
    EXPECT_TRUE(wake_preset.config.enable_live_pic_mcc);
}

TEST(SurfaceChargingSmokeTest, ReferenceModelRamCurrentExceedsWakeCurrent)
{
    SCDAT::Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "leo_ref_ram_facing", preset));

    SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceModel model;
    SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceConfig config;
    config.plasma = preset.config.plasma;
    config.electron_spectrum = preset.config.electron_spectrum;
    config.ion_spectrum = preset.config.ion_spectrum;
    config.has_electron_spectrum = preset.config.has_electron_spectrum;
    config.has_ion_spectrum = preset.config.has_ion_spectrum;
    config.patch_material = preset.config.material;
    config.body_material = preset.config.material;
    config.body_material.setType(SCDAT::Mesh::MaterialType::CONDUCTOR);
    config.body_photo_current_density_a_per_m2 = 0.0;
    config.patch_photo_current_density_a_per_m2 = 0.0;
    config.enable_ram_current = true;
    config.bulk_flow_velocity_m_per_s = preset.config.bulk_flow_velocity_m_per_s;
    config.flow_alignment_cosine = 1.0;
    config.patch_flow_angle_deg = 0.0;
    config.patch_conductivity_s_per_m = preset.config.material.getConductivity();
    config.patch_thickness_m = preset.config.dielectric_thickness_m;
    config.electron_collection_coefficient = preset.config.electron_collection_coefficient;
    config.ion_collection_coefficient = preset.config.ion_collection_coefficient;

    ASSERT_TRUE(model.configure(config));
    const auto ram_terms =
        model.evaluate(SCDAT::Toolkit::SurfaceCharging::ReferenceSurfaceRole::Patch, 0.0, 0.0);

    config.patch_flow_angle_deg = 180.0;
    ASSERT_TRUE(model.configure(config));
    const auto wake_terms =
        model.evaluate(SCDAT::Toolkit::SurfaceCharging::ReferenceSurfaceRole::Patch, 0.0, 0.0);

    EXPECT_GT(ram_terms.ram_ion_a_per_m2, wake_terms.ram_ion_a_per_m2);
}

TEST(SurfaceChargingSmokeTest, ReferenceModelSupportsMultipleSeeModels)
{
    SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceConfig config;
    config.plasma.electron_density_m3 = 8.0e6;
    config.plasma.ion_density_m3 = 4.0e6;
    config.plasma.electron_temperature_ev = 1500.0;
    config.plasma.ion_temperature_ev = 8.0;
    config.plasma.ion_mass_amu = 1.0;
    config.patch_material = SCDAT::Material::MaterialProperty(
        2, SCDAT::Mesh::MaterialType::DIELECTRIC, "ptfe");
    config.patch_material.setSecondaryElectronYield(1.4);
    config.patch_material.setScalarProperty("secondary_yield_peak_energy_ev", 500.0);
    config.patch_material.setScalarProperty("sims_exponent_n", 1.7);
    config.patch_material.setScalarProperty("katz_r1", 0.95);
    config.patch_material.setScalarProperty("katz_n1", 0.32);
    config.patch_material.setScalarProperty("katz_r2", 0.14);
    config.patch_material.setScalarProperty("katz_n2", 1.12);
    config.patch_material.setScalarProperty("atomic_number", 9.0);
    config.patch_material.setConductivity(1.0e-17);
    config.body_material = config.patch_material;
    config.body_material.setType(SCDAT::Mesh::MaterialType::CONDUCTOR);

    SCDAT::Toolkit::SurfaceCharging::ReferenceCurrentBalanceModel model;
    config.see_model = SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel::Whipple;
    ASSERT_TRUE(model.configure(config));
    const auto whipple =
        model.evaluate(SCDAT::Toolkit::SurfaceCharging::ReferenceSurfaceRole::Patch, 0.0, -50.0);

    config.see_model = SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel::Sims;
    ASSERT_TRUE(model.configure(config));
    const auto sims =
        model.evaluate(SCDAT::Toolkit::SurfaceCharging::ReferenceSurfaceRole::Patch, 0.0, -50.0);

    config.see_model = SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel::Katz;
    ASSERT_TRUE(model.configure(config));
    const auto katz =
        model.evaluate(SCDAT::Toolkit::SurfaceCharging::ReferenceSurfaceRole::Patch, 0.0, -50.0);

    EXPECT_TRUE(std::isfinite(whipple.secondary_electron_a_per_m2));
    EXPECT_TRUE(std::isfinite(sims.secondary_electron_a_per_m2));
    EXPECT_TRUE(std::isfinite(katz.secondary_electron_a_per_m2));
    EXPECT_NE(whipple.secondary_electron_a_per_m2, sims.secondary_electron_a_per_m2);
}

TEST(SurfaceChargingSmokeTest, ThrusterPlumeFloatingPotentialIsModeratelyNegative)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "thruster_plume_dielectric", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    const double floating_potential = charging.computeFloatingPotential();
    EXPECT_TRUE(std::isfinite(floating_potential));
    EXPECT_LT(floating_potential, -1.0);
    EXPECT_GT(floating_potential, -25.0);
}

TEST(SurfaceChargingSmokeTest, GeoReferenceModeTracksEarlyWorkbookTransient)
{
    DensePlasmaSurfaceCharging charging;
    SurfaceChargingScenarioPreset preset;
    ASSERT_TRUE(SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
        "geo_ecss_kapton_ref", preset));

    ASSERT_TRUE(charging.initialize(preset.config));
    for (int i = 0; i < 32; ++i)
    {
        ASSERT_TRUE(charging.advance(preset.time_step_s));
    }

    EXPECT_LT(charging.getStatus().state.surface_potential_v, -1.0e2);
    EXPECT_GT(charging.getStatus().state.surface_potential_v, -3.0e2);
    EXPECT_TRUE(std::isfinite(charging.getStatus().state.surface_potential_v));
}
