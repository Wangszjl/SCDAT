#include "SurfaceScenarioLoader.h"

#include "SurfaceScenarioCatalog.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace SCDAT
{
namespace MainEntry
{

namespace detail
{

std::string readTextFile(const std::filesystem::path& path);
std::optional<std::string> extractObjectField(const std::string& text,
                                              const std::string& key);
std::optional<std::string> extractArrayField(const std::string& text,
                                             const std::string& key);
std::optional<std::string> extractStringField(const std::string& text,
                                              const std::string& key);
std::optional<double> extractNumberField(const std::string& text,
                                         const std::string& key);
std::optional<bool> extractBoolField(const std::string& text,
                                     const std::string& key);
std::vector<std::string> parseStringArray(const std::string& array_text);

void applyUnifiedSolverConfigFromJsonObject(
    const std::string& text,
    Coupling::Contracts::SolverConfig& solver_config);
void applyReproducibilityConfigFromJsonObject(const std::string& text,
                                              unsigned int& seed,
                                              std::string& sampling_policy);

std::optional<Toolkit::SurfaceCharging::SurfaceRuntimeRoute>
parseSurfaceRuntimeRoute(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SurfacePicStrategy>
parseSurfacePicStrategy(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SurfaceLegacyInputAdapterKind>
parseSurfaceLegacyInputAdapterKind(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SurfacePicRuntimeKind>
parseSurfacePicRuntimeKind(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SurfaceInstrumentSetKind>
parseSurfaceInstrumentSetKind(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode>
parseSurfaceCurrentAlgorithmMode(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SurfaceBenchmarkMode>
parseSurfaceBenchmarkMode(const std::string& text);
std::optional<Toolkit::SurfaceCharging::VolumeLinearSolverPolicy>
parseVolumeLinearSolverPolicy(const std::string& text);
std::optional<FieldSolver::NativeVolumeBoundaryConditionFamily>
parseNativeVolumeBoundaryConditionFamily(const std::string& text);
std::optional<FieldSolver::NativeVolumeFieldFamily>
parseNativeVolumeFieldFamily(const std::string& text);
std::optional<FieldSolver::NativeVolumeDistributionFamily>
parseNativeVolumeDistributionFamily(const std::string& text);
std::optional<FieldSolver::NativeVolumeInteractionFamily>
parseNativeVolumeInteractionFamily(const std::string& text);
std::optional<Toolkit::SurfaceCharging::SecondaryElectronEmissionModel>
parseSecondaryElectronEmissionModel(const std::string& text);
std::optional<Toolkit::SurfaceCharging::ElectronCollectionModelKind>
parseElectronCollectionModelKind(const std::string& text);
std::optional<Toolkit::PlasmaAnalysis::PlasmaDistributionModelKind>
parsePlasmaDistributionModelKind(const std::string& text);

void applyPlasmaParametersFromJsonObject(const std::string& object_text,
                                         Toolkit::PlasmaAnalysis::PlasmaParameters& target);
void applyEmissionParametersFromJsonObject(
    const std::string& object_text,
    Toolkit::SurfaceCharging::EmissionModelParameters& target);
const Material::MaterialProperty*
resolveMaterialAliasOrName(const std::string& text);
void applyMaterialFromJsonObject(const std::string& object_text,
                                 Material::MaterialProperty& material);
bool applySurfaceInteractorFamilySelection(const std::string& text,
                                           Material::MaterialProperty& material);
void applyStructuredTopologyFromJson(const std::string& config_scope,
                                     Toolkit::SurfaceCharging::SurfaceChargingConfig& config);
bool applySpectrumFromJsonObject(const std::string& spectrum_text,
                                 Particle::ParticleType particle_type,
                                 Particle::ResolvedSpectrum& target);

} // namespace detail

namespace
{

namespace fs = std::filesystem;

} // namespace

namespace detail
{

std::string normalizeInteractorFamilyToken(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    std::string normalized;
    normalized.reserve(text.size());
    for (const unsigned char ch : text)
    {
        if (std::isalnum(ch))
        {
            normalized.push_back(static_cast<char>(ch));
        }
    }
    return normalized;
}

bool applySurfaceInteractorFamilySelection(const std::string& text,
                                           Material::MaterialProperty& material)
{
    const auto token = normalizeInteractorFamilyToken(text);
    if (token.empty())
    {
        return false;
    }

    const auto set_family_code = [&](double family_code) {
        material.setScalarProperty("surface_interactor_family_code", family_code);
        material.setScalarProperty("surface_interactor_use_multiple_model", 0.0);
        material.setScalarProperty("surface_interactor_use_tabulated_sey_model", 0.0);
        material.setScalarProperty("surface_interactor_use_erosion_model", 0.0);
        material.setScalarProperty("surface_interactor_use_reflection_model", 0.0);
        material.setScalarProperty("surface_interactor_use_yield_model", 0.0);
        material.setScalarProperty("surface_interactor_use_photoemission_model", 0.0);
        material.setScalarProperty("surface_interactor_use_induced_conduct_model", 0.0);
    };

    if (token == "spismaterialmodelsurfaceinteractorv1" || token == "materialmodel" ||
        token == "material")
    {
        set_family_code(0.0);
        return true;
    }
    if (token == "spismaterialdfsurfaceinteractorv1" || token == "materialdf")
    {
        set_family_code(1.0);
        return true;
    }
    if (token == "spisgenericdfsurfaceinteractorv1" || token == "genericdf" ||
        token == "generic")
    {
        set_family_code(2.0);
        return true;
    }
    if (token == "spismaxwelliansurfaceinteractorv1" || token == "maxwellian")
    {
        set_family_code(3.0);
        return true;
    }
    if (token == "spismaxwelliansurfaceinteractorwithrecollectionv1" ||
        token == "maxwellianwithrecollection" || token == "maxwellianrecollection")
    {
        set_family_code(4.0);
        return true;
    }
    if (token == "spismultiplesurfaceinteractorv1" || token == "multiple")
    {
        set_family_code(5.0);
        return true;
    }
    if (token == "spismultiplemaxwelliansurfaceinteractorv1" ||
        token == "multiplemaxwellian")
    {
        set_family_code(6.0);
        return true;
    }
    if (token == "spisyieldsurfaceinteractorv1" || token == "yield")
    {
        set_family_code(7.0);
        return true;
    }
    if (token == "spisreflectionsurfaceinteractorv1" || token == "reflection")
    {
        set_family_code(8.0);
        return true;
    }
    if (token == "spiserosionsurfaceinteractorv1" || token == "erosion")
    {
        set_family_code(9.0);
        return true;
    }
    if (token == "spisimprovedphotoemissionsurfaceinteractorv1" ||
        token == "improvedphotoemission" || token == "photoemission" || token == "photoem")
    {
        set_family_code(10.0);
        return true;
    }
    if (token == "spisbasicinducedconductsurfaceinteractorv1" ||
        token == "basicinducedconduction" || token == "inducedconduction" ||
        token == "inducedconduct")
    {
        set_family_code(11.0);
        return true;
    }
    if (token == "spistabulatedseysurfaceinteractorv1" || token == "tabulatedsey" ||
        token == "tabulatedsecondaryyield")
    {
        set_family_code(12.0);
        return true;
    }
    if (token == "spisrecollectedseysurfaceinteractorv1" || token == "recollectedsey")
    {
        set_family_code(13.0);
        return true;
    }
    if (token == "spisdefaultpeesurfaceinteractorv1" || token == "spisdefaultpeemodelv1" ||
        token == "defaultpee" || token == "defaultpeemodel")
    {
        set_family_code(14.0);
        return true;
    }
    if (token == "spisdefaultseeesurfaceinteractorv1" || token == "spisdefaultseeeymodelv1" ||
        token == "defaultseee" || token == "defaultseeeymodel")
    {
        set_family_code(15.0);
        return true;
    }
    if (token == "spisdefaultseepsurfaceinteractorv1" || token == "spisdefaultseepmodelv1" ||
        token == "defaultseep" || token == "defaultseepmodel")
    {
        set_family_code(16.0);
        return true;
    }
    if (token == "spisdefaulterosionsurfaceinteractorv1" ||
        token == "spisdefaulterosionmodelv1" || token == "defaulterosion" ||
        token == "defaulterosionmodel")
    {
        set_family_code(17.0);
        return true;
    }
    if (token == "spisdevicesurfaceinteractorv1" || token == "device")
    {
        set_family_code(18.0);
        return true;
    }
    return false;
}

} // namespace detail

Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset
SurfaceScenarioLoader::loadFromJson(const std::filesystem::path& json_path,
                                    std::filesystem::path& output_path) const
{
    const std::string content = detail::readTextFile(json_path);
    const auto run_scope = detail::extractObjectField(content, "run").value_or(content);
    const auto config_scope = detail::extractObjectField(content, "config").value_or(content);

    const Toolkit::SurfaceCharging::SurfaceScenarioCatalog scenario_catalog;
    auto preset = scenario_catalog.makeDefaultPreset();
    const auto base_preset = detail::extractStringField(content, "base_preset")
                                 .value_or(detail::extractStringField(content, "preset").value_or(""));
    if (!base_preset.empty())
    {
        if (!scenario_catalog.tryGetPreset(base_preset, preset))
        {
            throw std::runtime_error("Unknown surface base preset in json: " + base_preset);
        }
    }

    if (const auto name = detail::extractStringField(content, "name"); name)
    {
        preset.name = *name;
    }

    if (const auto time_step_s = detail::extractNumberField(run_scope, "time_step_s"); time_step_s)
    {
        preset.time_step_s = *time_step_s;
    }
    if (const auto steps = detail::extractNumberField(run_scope, "steps"); steps)
    {
        preset.steps = static_cast<std::size_t>(std::max(0.0, std::floor(*steps + 0.5)));
    }
    if (const auto adaptive = detail::extractBoolField(run_scope, "adaptive_time_stepping"); adaptive)
    {
        preset.adaptive_time_stepping = *adaptive;
    }
    if (const auto total_duration = detail::extractNumberField(run_scope, "total_duration_s");
        total_duration)
    {
        preset.total_duration_s = *total_duration;
    }
    if (const auto minimum_dt = detail::extractNumberField(run_scope, "minimum_time_step_s"); minimum_dt)
    {
        preset.minimum_time_step_s = *minimum_dt;
    }
    if (const auto maximum_dt = detail::extractNumberField(run_scope, "maximum_time_step_s"); maximum_dt)
    {
        preset.maximum_time_step_s = *maximum_dt;
    }

    if (const auto output_csv = detail::extractStringField(content, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else if (const auto output_csv = detail::extractStringField(run_scope, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else
    {
        output_path = preset.default_output_csv;
    }

    auto apply_bool = [&](const std::string& key, bool& target) {
        if (const auto value = detail::extractBoolField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = detail::extractNumberField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_size_t = [&](const std::string& key, std::size_t& target) {
        if (const auto value = detail::extractNumberField(config_scope, key); value)
        {
            target = static_cast<std::size_t>(std::max(0.0, std::floor(*value + 0.5)));
        }
    };
    auto apply_path = [&](const std::string& key, fs::path& target) {
        if (const auto value = detail::extractStringField(config_scope, key); value)
        {
            target = *value;
        }
    };

    detail::applyUnifiedSolverConfigFromJsonObject(config_scope, preset.config.solver_config);
    detail::applyReproducibilityConfigFromJsonObject(config_scope, preset.config.seed,
                                                     preset.config.sampling_policy);

    apply_number("surface_area_m2", preset.config.surface_area_m2);
    apply_number("capacitance_per_area_f_per_m2", preset.config.capacitance_per_area_f_per_m2);
    apply_number("dielectric_thickness_m", preset.config.dielectric_thickness_m);
    apply_number("bulk_flow_velocity_m_per_s", preset.config.bulk_flow_velocity_m_per_s);
    apply_number("flow_alignment_cosine", preset.config.flow_alignment_cosine);
    apply_number("electron_flow_coupling", preset.config.electron_flow_coupling);
    apply_number("ion_directed_velocity_m_per_s", preset.config.ion_directed_velocity_m_per_s);
    apply_number("live_pic_probe_delta_v", preset.config.live_pic_probe_delta_v);
    apply_number("live_pic_reference_potential_v", preset.config.live_pic_reference_potential_v);
    apply_number("pic_recalibration_trigger_v", preset.config.pic_recalibration_trigger_v);
    apply_number("body_initial_potential_v", preset.config.body_initial_potential_v);
    apply_number("body_capacitance_f", preset.config.body_capacitance_f);
    apply_number("body_photo_current_density_a_per_m2",
                 preset.config.body_photo_current_density_a_per_m2);
    apply_number("patch_photo_current_density_a_per_m2",
                 preset.config.patch_photo_current_density_a_per_m2);
    apply_number("photoelectron_temperature_ev", preset.config.photoelectron_temperature_ev);
    if (const auto value = detail::extractNumberField(config_scope, "body_photoelectron_temperature_ev");
        value)
    {
        preset.config.material.setScalarProperty("body_photoelectron_temperature_ev",
                                                 std::max(1.0e-3, *value));
    }
    if (const auto value = detail::extractNumberField(config_scope, "body_atomic_number"); value)
    {
        preset.config.material.setScalarProperty("body_atomic_number", std::max(0.0, *value));
    }
    if (const auto value = detail::extractNumberField(config_scope, "body_secondary_electron_yield");
        value)
    {
        preset.config.material.setScalarProperty("body_secondary_electron_yield",
                                                 std::max(0.0, *value));
    }
    if (const auto value =
            detail::extractNumberField(config_scope, "body_secondary_yield_peak_energy_ev");
        value)
    {
        preset.config.material.setScalarProperty("body_secondary_yield_peak_energy_ev",
                                                 std::max(0.0, *value));
    }
    if (const auto value = detail::extractNumberField(config_scope, "body_photo_emission_scale"); value)
    {
        preset.config.material.setScalarProperty("body_photo_emission_scale",
                                                 std::max(0.0, *value));
    }
    if (const auto value = detail::extractNumberField(config_scope, "body_electron_collection_coefficient");
        value)
    {
        preset.config.material.setScalarProperty("body_electron_collection_coefficient",
                                                 std::max(0.0, *value));
    }
    if (const auto value = detail::extractNumberField(config_scope, "body_electron_collection_scale");
        value)
    {
        preset.config.material.setScalarProperty("body_electron_collection_scale",
                                                 std::max(0.0, *value));
    }
    apply_number("electron_collection_coefficient",
                 preset.config.electron_collection_coefficient);
    apply_number("ion_collection_coefficient", preset.config.ion_collection_coefficient);
    apply_number("max_delta_potential_v_per_step", preset.config.max_delta_potential_v_per_step);
    apply_number("plasma_electron_density_m3", preset.config.plasma.electron_density_m3);
    apply_number("plasma_ion_density_m3", preset.config.plasma.ion_density_m3);
    apply_number("plasma_electron_temperature_ev", preset.config.plasma.electron_temperature_ev);
    apply_number("plasma_ion_temperature_ev", preset.config.plasma.ion_temperature_ev);
    apply_number("plasma_ion_mass_amu", preset.config.plasma.ion_mass_amu);
    apply_size_t("live_pic_window_steps", preset.config.live_pic_window_steps);
    apply_size_t("live_pic_window_layers", preset.config.live_pic_window_layers);
    apply_size_t("live_pic_particles_per_element", preset.config.live_pic_particles_per_element);
    apply_size_t("pic_recalibration_interval_steps", preset.config.pic_recalibration_interval_steps);
    apply_size_t("pic_calibration_samples", preset.config.pic_calibration_samples);
    apply_size_t("internal_substeps", preset.config.internal_substeps);

    apply_bool("floating", preset.config.floating);
    apply_bool("derive_capacitance_from_material", preset.config.derive_capacitance_from_material);
    apply_bool("use_reference_current_balance", preset.config.use_reference_current_balance);
    apply_bool("enable_pic_calibration", preset.config.enable_pic_calibration);
    apply_bool("enable_live_pic_window", preset.config.enable_live_pic_window);
    apply_bool("enable_live_pic_mcc", preset.config.enable_live_pic_mcc);
    apply_bool("enable_body_patch_circuit", preset.config.enable_body_patch_circuit);
    apply_bool("body_floating", preset.config.body_floating);
    apply_bool("enable_secondary_electron", preset.config.enable_secondary_electron);
    apply_bool("enable_backscatter", preset.config.enable_backscatter);
    apply_bool("enable_photoelectron", preset.config.enable_photoelectron);
    apply_bool("enable_external_field_solver_bridge",
               preset.config.enable_external_field_solver_bridge);
    apply_bool("enable_external_volume_solver_bridge",
               preset.config.enable_external_volume_solver_bridge);
    apply_bool("has_electron_spectrum", preset.config.has_electron_spectrum);
    apply_bool("has_ion_spectrum", preset.config.has_ion_spectrum);

    apply_path("external_field_solver_request_path", preset.config.external_field_solver_request_path);
    apply_path("external_field_solver_result_path", preset.config.external_field_solver_result_path);
    apply_path("external_volume_solver_request_path", preset.config.external_volume_solver_request_path);
    apply_path("external_volume_solver_result_path", preset.config.external_volume_solver_result_path);
    apply_path("external_volume_mesh_path", preset.config.external_volume_mesh_path);
    apply_path("external_surface_volume_projection_path",
               preset.config.external_surface_volume_projection_path);

    auto native_volume_boundary_condition_families =
        detail::extractArrayField(config_scope,
                                  "surface_native_volume_boundary_condition_families");
    if (!native_volume_boundary_condition_families)
    {
        native_volume_boundary_condition_families =
            detail::extractArrayField(config_scope,
                                      "native_volume_boundary_condition_families");
    }
    if (native_volume_boundary_condition_families)
    {
        preset.config.native_volume_boundary_condition_families.clear();
        for (const auto& family_text :
             detail::parseStringArray(*native_volume_boundary_condition_families))
        {
            const auto parsed =
                detail::parseNativeVolumeBoundaryConditionFamily(family_text);
            if (!parsed)
            {
                throw std::runtime_error(
                    "Unsupported native_volume_boundary_condition_families entry in "
                    "surface config json: " +
                    family_text);
            }
            preset.config.native_volume_boundary_condition_families.push_back(*parsed);
        }
    }

    auto native_volume_field_families =
        detail::extractArrayField(config_scope, "surface_native_volume_field_families");
    if (!native_volume_field_families)
    {
        native_volume_field_families =
            detail::extractArrayField(config_scope, "native_volume_field_families");
    }
    if (native_volume_field_families)
    {
        preset.config.native_volume_field_families.clear();
        for (const auto& family_text : detail::parseStringArray(*native_volume_field_families))
        {
            const auto parsed = detail::parseNativeVolumeFieldFamily(family_text);
            if (!parsed)
            {
                throw std::runtime_error(
                    "Unsupported native_volume_field_families entry in surface config json: " +
                    family_text);
            }
            preset.config.native_volume_field_families.push_back(*parsed);
        }
    }

    auto native_volume_distribution_families =
        detail::extractArrayField(config_scope,
                                  "surface_native_volume_distribution_families");
    if (!native_volume_distribution_families)
    {
        native_volume_distribution_families =
            detail::extractArrayField(config_scope,
                                      "native_volume_distribution_families");
    }
    if (native_volume_distribution_families)
    {
        preset.config.native_volume_distribution_families.clear();
        for (const auto& family_text :
             detail::parseStringArray(*native_volume_distribution_families))
        {
            const auto parsed =
                detail::parseNativeVolumeDistributionFamily(family_text);
            if (!parsed)
            {
                throw std::runtime_error(
                    "Unsupported native_volume_distribution_families entry in "
                    "surface config json: " +
                    family_text);
            }
            preset.config.native_volume_distribution_families.push_back(*parsed);
        }
    }

    auto native_volume_interaction_families =
        detail::extractArrayField(config_scope,
                                  "surface_native_volume_interaction_families");
    if (!native_volume_interaction_families)
    {
        native_volume_interaction_families =
            detail::extractArrayField(config_scope,
                                      "native_volume_interaction_families");
    }
    if (native_volume_interaction_families)
    {
        preset.config.native_volume_interaction_families.clear();
        for (const auto& family_text :
             detail::parseStringArray(*native_volume_interaction_families))
        {
            const auto parsed =
                detail::parseNativeVolumeInteractionFamily(family_text);
            if (!parsed)
            {
                throw std::runtime_error(
                    "Unsupported native_volume_interaction_families entry in "
                    "surface config json: " +
                    family_text);
            }
            preset.config.native_volume_interaction_families.push_back(*parsed);
        }
    }

    auto runtime_route = detail::extractStringField(config_scope, "surface_runtime_route");
    if (!runtime_route)
    {
        runtime_route = detail::extractStringField(config_scope, "runtime_route");
    }
    if (runtime_route)
    {
        const auto parsed = detail::parseSurfaceRuntimeRoute(*runtime_route);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported runtime_route in surface config json: " +
                                     *runtime_route);
        }
        preset.config.runtime_route = *parsed;
    }

    if (const auto surface_pic_strategy = detail::extractStringField(config_scope, "surface_pic_strategy");
        surface_pic_strategy)
    {
        const auto parsed = detail::parseSurfacePicStrategy(*surface_pic_strategy);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported surface_pic_strategy in surface config json: " +
                *surface_pic_strategy);
        }
        preset.config.surface_pic_strategy = *parsed;
    }

    if (const auto legacy_input_adapter =
            detail::extractStringField(config_scope, "legacy_input_adapter");
        legacy_input_adapter)
    {
        const auto parsed = detail::parseSurfaceLegacyInputAdapterKind(*legacy_input_adapter);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported legacy_input_adapter in surface config json: " +
                *legacy_input_adapter);
        }
        preset.config.legacy_input_adapter_kind = *parsed;
    }

    if (const auto surface_pic_runtime = detail::extractStringField(config_scope, "surface_pic_runtime");
        surface_pic_runtime)
    {
        const auto parsed = detail::parseSurfacePicRuntimeKind(*surface_pic_runtime);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported surface_pic_runtime in surface config json: " +
                *surface_pic_runtime);
        }
        preset.config.surface_pic_runtime_kind = *parsed;
    }

    if (const auto surface_instrument_set =
            detail::extractStringField(config_scope, "surface_instrument_set");
        surface_instrument_set)
    {
        const auto parsed = detail::parseSurfaceInstrumentSetKind(*surface_instrument_set);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported surface_instrument_set in surface config json: " +
                *surface_instrument_set);
        }
        preset.config.surface_instrument_set_kind = *parsed;
    }

    if (const auto reference_family = detail::extractStringField(config_scope, "reference_family");
        reference_family)
    {
        preset.config.reference_family = *reference_family;
    }
    if (const auto reference_case_id = detail::extractStringField(config_scope, "reference_case_id");
        reference_case_id)
    {
        preset.config.reference_case_id = *reference_case_id;
    }
    if (const auto reference_matrix_case_id =
            detail::extractStringField(config_scope, "reference_matrix_case_id");
        reference_matrix_case_id)
    {
        preset.config.reference_matrix_case_id = *reference_matrix_case_id;
    }

    if (const auto algorithm_mode = detail::extractStringField(config_scope, "current_algorithm_mode");
        algorithm_mode)
    {
        const auto parsed = detail::parseSurfaceCurrentAlgorithmMode(*algorithm_mode);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported current_algorithm_mode in surface config json: " +
                *algorithm_mode);
        }
        preset.config.current_algorithm_mode = *parsed;
    }

    if (const auto benchmark_mode = detail::extractStringField(config_scope, "benchmark_mode");
        benchmark_mode)
    {
        const auto parsed = detail::parseSurfaceBenchmarkMode(*benchmark_mode);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported benchmark_mode in surface config json: " +
                                     *benchmark_mode);
        }
        preset.config.benchmark_mode = *parsed;
    }

    if (const auto cross_section_set_id =
            detail::extractStringField(config_scope, "live_pic_collision_cross_section_set_id");
        cross_section_set_id)
    {
        preset.config.live_pic_collision_cross_section_set_id = *cross_section_set_id;
    }

    if (const auto solver_policy = detail::extractStringField(config_scope, "volume_linear_solver_policy");
        solver_policy)
    {
        const auto parsed = detail::parseVolumeLinearSolverPolicy(*solver_policy);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported volume_linear_solver_policy in surface config json: " +
                *solver_policy);
        }
        preset.config.volume_linear_solver_policy = *parsed;
    }

    if (const auto see_model = detail::extractStringField(config_scope, "reference_see_model");
        see_model)
    {
        const auto parsed_see_model = detail::parseSecondaryElectronEmissionModel(*see_model);
        if (!parsed_see_model)
        {
            throw std::runtime_error("Unsupported reference_see_model in surface config json: " +
                                     *see_model);
        }
        preset.config.reference_see_model = *parsed_see_model;
    }

    if (const auto model_text = detail::extractStringField(config_scope, "electron_collection_model");
        model_text)
    {
        const auto parsed_model = detail::parseElectronCollectionModelKind(*model_text);
        if (!parsed_model)
        {
            throw std::runtime_error(
                "Unsupported electron_collection_model in surface config json: " + *model_text);
        }
        preset.config.electron_collection_model = *parsed_model;
    }

    auto distribution_model = detail::extractStringField(config_scope, "distribution_model");
    if (!distribution_model)
    {
        distribution_model = detail::extractStringField(config_scope, "plasma_distribution_model");
    }
    if (distribution_model)
    {
        const auto parsed_distribution_model =
            detail::parsePlasmaDistributionModelKind(*distribution_model);
        if (!parsed_distribution_model)
        {
            throw std::runtime_error(
                "Unsupported distribution_model in surface config json: " +
                *distribution_model);
        }
        preset.config.distribution_model = *parsed_distribution_model;
    }

    auto plasma_text = detail::extractObjectField(config_scope, "plasma_model");
    if (!plasma_text)
    {
        plasma_text = detail::extractObjectField(config_scope, "plasma");
    }
    if (plasma_text)
    {
        detail::applyPlasmaParametersFromJsonObject(*plasma_text, preset.config.plasma);
    }

    auto emission_text = detail::extractObjectField(config_scope, "surface_emission");
    if (!emission_text)
    {
        emission_text = detail::extractObjectField(config_scope, "emission");
    }
    if (emission_text)
    {
        detail::applyEmissionParametersFromJsonObject(*emission_text, preset.config.emission);
    }

    if (const auto material_library_path =
            detail::extractStringField(config_scope, "material_library_path");
        material_library_path)
    {
        preset.config.material_library_path = *material_library_path;
    }
    if (const auto imported_material_name =
            detail::extractStringField(config_scope, "imported_material_name");
        imported_material_name)
    {
        preset.config.imported_material_name = *imported_material_name;
    }

    if (const auto material_alias = detail::extractStringField(config_scope, "material_alias");
        material_alias)
    {
        const auto* resolved = detail::resolveMaterialAliasOrName(*material_alias);
        if (resolved == nullptr)
        {
            throw std::runtime_error("Unknown material_alias in surface config json: " +
                                     *material_alias);
        }
        preset.config.material = *resolved;
    }
    if (const auto material_name = detail::extractStringField(config_scope, "material_name");
        material_name)
    {
        const auto* resolved = detail::resolveMaterialAliasOrName(*material_name);
        if (resolved == nullptr)
        {
            if (!preset.config.material_library_path.empty())
            {
                preset.config.imported_material_name = *material_name;
            }
            else
            {
                throw std::runtime_error("Unknown material_name in surface config json: " +
                                         *material_name);
            }
        }
        else
        {
            preset.config.material = *resolved;
        }
    }
    if (const auto interactor_family =
            detail::extractStringField(config_scope, "surface_interactor_family");
        interactor_family)
    {
        if (!detail::applySurfaceInteractorFamilySelection(*interactor_family,
                                                           preset.config.material))
        {
            throw std::runtime_error(
                "Unsupported surface_interactor_family in surface config json: " +
                *interactor_family);
        }
    }

    auto material_text = detail::extractObjectField(config_scope, "surface_material");
    if (!material_text)
    {
        material_text = detail::extractObjectField(config_scope, "material");
    }
    if (material_text)
    {
        detail::applyMaterialFromJsonObject(*material_text, preset.config.material);
    }

    detail::applyStructuredTopologyFromJson(config_scope, preset.config);

    if (const auto spectrum_text = detail::extractObjectField(config_scope, "electron_spectrum");
        spectrum_text)
    {
        if (!detail::applySpectrumFromJsonObject(*spectrum_text, Particle::ParticleType::ELECTRON,
                                                 preset.config.electron_spectrum))
        {
            throw std::runtime_error("Failed to parse electron_spectrum in surface config json.");
        }
        preset.config.has_electron_spectrum = true;
    }
    if (const auto spectrum_text = detail::extractObjectField(config_scope, "ion_spectrum");
        spectrum_text)
    {
        if (!detail::applySpectrumFromJsonObject(*spectrum_text, Particle::ParticleType::ION,
                                                 preset.config.ion_spectrum))
        {
            throw std::runtime_error("Failed to parse ion_spectrum in surface config json.");
        }
        preset.config.has_ion_spectrum = true;
    }

    return preset;
}

} // namespace MainEntry
} // namespace SCDAT
