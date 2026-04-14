#include "LegacyBenchmarkSupport.h"
#include "../../Tools/Basic/include/NumericAggregation.h"
#include "../../Tools/Basic/include/StringTokenUtils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

std::filesystem::path legacyBenchmarkRepoRoot()
{
    static const std::filesystem::path root = []() {
        auto cursor = std::filesystem::weakly_canonical(std::filesystem::path(__FILE__)).parent_path();
        for (int depth = 0; depth < 10; ++depth)
        {
            if (std::filesystem::exists(cursor / "ref") &&
                std::filesystem::exists(cursor / "Toolkit") &&
                std::filesystem::exists(cursor / "Tools"))
            {
                return cursor;
            }
            if (!cursor.has_parent_path())
            {
                break;
            }
            cursor = cursor.parent_path();
        }

        return std::filesystem::weakly_canonical(std::filesystem::path(__FILE__))
            .parent_path()
            .parent_path()
            .parent_path()
            .parent_path();
    }();
    return root;
}

std::filesystem::path legacyBenchmarkAbsolutePath(const std::filesystem::path& relative_path)
{
    if (relative_path.empty() || relative_path.is_absolute())
    {
        return relative_path;
    }
    return legacyBenchmarkRepoRoot() / relative_path;
}

bool startsWith(const std::string& value, const std::string& prefix)
{
    return value.rfind(prefix, 0) == 0;
}

std::string sanitizedMaterialName(const std::string& raw_name)
{
    std::string value = Basic::trimAscii(raw_name);
    for (char& ch : value)
    {
        if (std::isspace(static_cast<unsigned char>(ch)))
        {
            ch = '_';
        }
        else if (ch == ',' || ch == '/' || ch == '(' || ch == ')' || ch == '-')
        {
            ch = '_';
        }
        else
        {
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
    }
    return value;
}

std::vector<std::string> splitWhitespace(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> parts;
    std::string token;
    while (stream >> token)
    {
        parts.push_back(token);
    }
    return parts;
}

bool isNumericLine(const std::string& line)
{
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return false;
    }
    const char c = line[first];
    return std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+';
}

std::vector<double> parseNumbersFromLine(const std::string& line)
{
    std::vector<double> values;
    std::istringstream stream(line);
    double value = 0.0;
    while (stream >> value)
    {
        values.push_back(value);
    }
    return values;
}

double weightedAverage(const std::vector<double>& weights, const std::vector<double>& values,
                       double fallback)
{
    double numerator = 0.0;
    double denominator = 0.0;
    const std::size_t count = std::min(weights.size(), values.size());
    for (std::size_t i = 0; i < count; ++i)
    {
        const double weight = std::max(0.0, weights[i]);
        numerator += weight * values[i];
        denominator += weight;
    }
    if (denominator > 0.0)
    {
        return numerator / denominator;
    }
    return fallback;
}

double legacyLeoOrbitalSpeedMPerS(double altitude_km)
{
    return std::sqrt(39.8634 / std::max(1.0, 6371.004 + altitude_km)) * 1000.0 * 100.0;
}

SecondaryElectronEmissionModel mapLegacySeeModel(int see_model_id)
{
    switch (see_model_id)
    {
    case 2:
        return SecondaryElectronEmissionModel::Sims;
    case 3:
        return SecondaryElectronEmissionModel::Katz;
    default:
        return SecondaryElectronEmissionModel::Whipple;
    }
}

std::optional<LegacyEnvironmentRecord> findLegacyEnvironmentRecord(
    const LegacyBenchmarkCaseDefinition& definition, int environment_id)
{
    for (const auto& record : definition.environments)
    {
        if (record.case_id == environment_id)
        {
            return record;
        }
    }
    if (!definition.environments.empty())
    {
        return definition.environments.front();
    }
    return std::nullopt;
}

std::optional<LegacyMaterialRecord> findLegacyMaterialRecord(
    const std::vector<LegacyMaterialRecord>& records, int material_id)
{
    for (const auto& record : records)
    {
        if (record.material_id == material_id)
        {
            return record;
        }
    }
    if (!records.empty())
    {
        return records.front();
    }
    return std::nullopt;
}

std::optional<LegacyMaterialRecord> findExactLegacyMaterialRecord(
    const std::vector<LegacyMaterialRecord>& records, int material_id)
{
    for (const auto& record : records)
    {
        if (record.material_id == material_id)
        {
            return record;
        }
    }
    return std::nullopt;
}

Material::MaterialProperty buildLegacyMaterial(const LegacyMaterialRecord& record,
                                               bool dielectric_patch,
                                               bool isolated_metal_patch,
                                               double explicit_permittivity,
                                               double explicit_conductivity)
{
    const auto material_type =
        dielectric_patch ? Mesh::MaterialType::DIELECTRIC : Mesh::MaterialType::CONDUCTOR;
    Material::MaterialProperty material(record.material_id, material_type,
                                        sanitizedMaterialName(record.name));
    material.setName(record.name);

    std::size_t cursor = 0;
    if (dielectric_patch && record.numeric_values.size() >= 13)
    {
        material.setPermittivity(std::max(1.0, record.numeric_values[cursor++]));
        material.setConductivity(std::max(0.0, record.numeric_values[cursor++]));
    }
    else
    {
        material.setPermittivity(std::max(1.0, explicit_permittivity));
        material.setConductivity(std::max(0.0, explicit_conductivity));
    }

    if (record.numeric_values.size() < cursor + 11)
    {
        return material;
    }

    const double atomic_number = record.numeric_values[cursor++];
    const double photo_current_density_na_per_m2 = record.numeric_values[cursor++];
    const double ion_secondary_yield = record.numeric_values[cursor++];
    const double ion_secondary_peak_energy_kev = record.numeric_values[cursor++];
    const double secondary_peak_yield = record.numeric_values[cursor++];
    const double secondary_peak_energy_kev = record.numeric_values[cursor++];
    const double sims_exponent_n = record.numeric_values[cursor++];
    const double katz_r1 = record.numeric_values[cursor++];
    const double katz_n1 = record.numeric_values[cursor++];
    const double katz_r2 = record.numeric_values[cursor++];
    const double katz_n2 = record.numeric_values[cursor++];

    material.setSecondaryElectronYield(std::max(0.0, secondary_peak_yield));
    material.setWorkFunctionEv(isolated_metal_patch ? 4.5 : 4.7);
    material.setScalarProperty("atomic_number", atomic_number);
    material.setScalarProperty("legacy_photo_current_density_a_per_m2",
                               photo_current_density_na_per_m2 * 1.0e-9);
    material.setScalarProperty("ion_secondary_yield", std::max(0.0, ion_secondary_yield));
    material.setScalarProperty("ion_secondary_peak_energy_kev",
                               std::max(1.0e-6, ion_secondary_peak_energy_kev));
    material.setScalarProperty("secondary_yield_peak_energy_ev",
                               std::max(1.0, secondary_peak_energy_kev * 1.0e3));
    material.setScalarProperty("sims_exponent_n", std::max(1.0, sims_exponent_n));
    material.setScalarProperty("katz_r1", katz_r1);
    material.setScalarProperty("katz_n1", katz_n1);
    material.setScalarProperty("katz_r2", katz_r2);
    material.setScalarProperty("katz_n2", katz_n2);
    material.setScalarProperty("secondary_emission_escape_energy_ev", 2.0);
    material.setScalarProperty("photoelectron_escape_energy_ev", 2.0);
    material.setScalarProperty("thermionic_escape_energy_ev", 0.2);
    material.setScalarProperty("body_conductivity_s_per_m",
                               isolated_metal_patch ? 3.5e7 : material.getConductivity());
    material.setScalarProperty("poole_frenkel_beta", 0.0);
    material.setScalarProperty("max_field_enhancement_factor", 1.0);
    return material;
}

void applyLegacyEnvironmentToConfig(const LegacyEnvironmentRecord& environment,
                                    SurfaceChargingConfig& config)
{
    std::vector<double> electron_densities_m3;
    std::vector<double> electron_temperatures_ev;
    std::vector<double> ion_densities_m3;
    std::vector<double> ion_temperatures_ev;
    std::vector<double> ion_masses_amu;

    if (environment.thermal_values.size() >= 15)
    {
        for (std::size_t i = 0; i < 3; ++i)
        {
            electron_densities_m3.push_back(
                std::max(0.0, environment.thermal_values[2 * i]) * 1.0e6);
            electron_temperatures_ev.push_back(
                std::max(1.0e-6, environment.thermal_values[2 * i + 1]));
        }
        for (std::size_t i = 0; i < 3; ++i)
        {
            const std::size_t base = 6 + 3 * i;
            ion_densities_m3.push_back(std::max(0.0, environment.thermal_values[base]) * 1.0e6);
            ion_temperatures_ev.push_back(
                std::max(1.0e-6, environment.thermal_values[base + 1]));
            ion_masses_amu.push_back(std::max(1.0, environment.thermal_values[base + 2]));
        }
    }

    config.plasma.electron_density_m3 =
        Basic::sumValue(electron_densities_m3);
    config.plasma.ion_density_m3 =
        Basic::sumValue(ion_densities_m3);
    config.plasma.electron_temperature_ev =
        weightedAverage(electron_densities_m3, electron_temperatures_ev, 1.0);
    config.plasma.ion_temperature_ev =
        weightedAverage(ion_densities_m3, ion_temperatures_ev, 1.0);
    config.plasma.ion_mass_amu = weightedAverage(ion_densities_m3, ion_masses_amu, 1.0);

    config.electron_spectrum = Particle::ResolvedSpectrum{};
    config.electron_spectrum.model = Particle::SpatialSamplingModel::TABULATED;
    config.electron_spectrum.energy_grid_ev = environment.electron_energy_ev;
    config.electron_spectrum.differential_number_flux = environment.electron_flux_per_m2_s_sr_ev;
    for (std::size_t i = 0; i < electron_densities_m3.size(); ++i)
    {
        if (electron_densities_m3[i] <= 0.0)
        {
            continue;
        }
        Particle::SpectrumPopulation population;
        population.density_m3 = electron_densities_m3[i];
        population.temperature_ev = electron_temperatures_ev[i];
        population.mass_amu = 5.48579909065e-4;
        population.charge_number = -1.0;
        config.electron_spectrum.populations.push_back(population);
    }

    config.ion_spectrum = Particle::ResolvedSpectrum{};
    config.ion_spectrum.model = Particle::SpatialSamplingModel::TABULATED;
    config.ion_spectrum.energy_grid_ev = environment.ion_energy_ev;
    config.ion_spectrum.differential_number_flux = environment.ion_flux_per_m2_s_sr_ev;
    for (std::size_t i = 0; i < ion_densities_m3.size(); ++i)
    {
        if (ion_densities_m3[i] <= 0.0)
        {
            continue;
        }
        Particle::SpectrumPopulation population;
        population.density_m3 = ion_densities_m3[i];
        population.temperature_ev = ion_temperatures_ev[i];
        population.mass_amu = ion_masses_amu[i];
        population.charge_number = 1.0;
        config.ion_spectrum.populations.push_back(population);
    }

    config.has_electron_spectrum = !config.electron_spectrum.populations.empty() ||
                                   !config.electron_spectrum.energy_grid_ev.empty();
    config.has_ion_spectrum =
        !config.ion_spectrum.populations.empty() || !config.ion_spectrum.energy_grid_ev.empty();
}

std::string joinNameTokens(const std::vector<std::string>& tokens,
                           std::size_t first_index,
                           std::size_t last_index_exclusive)
{
    std::string name;
    for (std::size_t i = first_index; i < last_index_exclusive; ++i)
    {
        if (!name.empty())
        {
            name += ' ';
        }
        name += tokens[i];
    }
    return name;
}

std::string joinStrings(const std::vector<std::string>& values, const char* delimiter)
{
    std::ostringstream output;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
        {
            output << delimiter;
        }
        output << values[i];
    }
    return output.str();
}

} // namespace

double legacyBenchmarkExecutionModeId(LegacyBenchmarkExecutionMode mode)
{
    return mode == LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm ? 1.0 : 0.0;
}

std::string legacyBenchmarkExecutionModeName(LegacyBenchmarkExecutionMode mode)
{
    return mode == LegacyBenchmarkExecutionMode::ExecuteLegacyAlgorithm
               ? "ExecuteLegacyAlgorithm"
               : "ReplayFromReference";
}

std::string legacyBenchmarkBaselineFamilyName(SurfaceBenchmarkSource source)
{
    switch (source)
    {
    case SurfaceBenchmarkSource::MatlabGeo:
    case SurfaceBenchmarkSource::MatlabLeo:
        return "MATLAB";
    case SurfaceBenchmarkSource::CGeo:
    case SurfaceBenchmarkSource::CLeoRam:
    case SurfaceBenchmarkSource::CLeoWake:
        return "LegacyC";
    default:
        return "None";
    }
}

std::string legacyBenchmarkBaselineOriginName(SurfaceBenchmarkSource source)
{
    switch (source)
    {
    case SurfaceBenchmarkSource::MatlabGeo:
        return "GeneratedFromMatlabReferenceScript";
    case SurfaceBenchmarkSource::MatlabLeo:
        return "LegacyLeoReferenceWithMatlabDiagnostics";
    case SurfaceBenchmarkSource::CGeo:
    case SurfaceBenchmarkSource::CLeoRam:
    case SurfaceBenchmarkSource::CLeoWake:
        return "ReferenceDatTable";
    default:
        return "None";
    }
}

LegacyBenchmarkAcceptanceCriteria legacyBenchmarkAcceptanceCriteria(
    SurfaceBenchmarkSource source)
{
    LegacyBenchmarkAcceptanceCriteria criteria;
    switch (source)
    {
    case SurfaceBenchmarkSource::CGeo:
        criteria.id = "legacy-c-geo-v1";
        criteria.note = "Legacy C GEO strict replay envelope.";
        criteria.focus_segment = "global";
        criteria.patch_rmse_v_max = 1.0;
        criteria.body_rmse_v_max = 1.0;
        criteria.body_je_rmse_a_per_m2_max = 2.0e-6;
        criteria.body_jnet_rmse_a_per_m2_max = 2.0e-6;
        break;
    case SurfaceBenchmarkSource::CLeoRam:
        criteria.id = "legacy-c-leo-ram-v1";
        criteria.note = "Legacy C LEO RAM with known body-current gap; tail monitored.";
        criteria.focus_segment = "negative_potential_tail";
        criteria.patch_rmse_v_max = 3.0e2;
        criteria.body_rmse_v_max = 5.5e3;
        criteria.body_je_rmse_a_per_m2_max = 1.0e-5;
        criteria.body_jnet_rmse_a_per_m2_max = 1.0e-5;
        criteria.negative_tail_body_je_rmse_a_per_m2_max = 8.0e-6;
        criteria.negative_tail_body_jnet_rmse_a_per_m2_max = 8.0e-6;
        break;
    case SurfaceBenchmarkSource::CLeoWake:
        criteria.id = "legacy-c-leo-wake-v1";
        criteria.note = "Legacy C LEO WAKE with known body-current gap; tail monitored.";
        criteria.focus_segment = "negative_potential_tail";
        criteria.patch_rmse_v_max = 3.0e2;
        criteria.body_rmse_v_max = 6.5e3;
        criteria.body_je_rmse_a_per_m2_max = 1.0e-5;
        criteria.body_jnet_rmse_a_per_m2_max = 1.0e-5;
        criteria.negative_tail_body_je_rmse_a_per_m2_max = 8.0e-6;
        criteria.negative_tail_body_jnet_rmse_a_per_m2_max = 8.0e-6;
        break;
    case SurfaceBenchmarkSource::MatlabGeo:
        criteria.id = "matlab-geo-v1";
        criteria.note = "MATLAB GEO replay target; prefix segment should be near machine-zero.";
        criteria.focus_segment = "global";
        criteria.patch_rmse_v_max = 1.0e-6;
        criteria.body_rmse_v_max = 1.0e-6;
        criteria.body_je_rmse_a_per_m2_max = 1.0e-7;
        criteria.body_jnet_rmse_a_per_m2_max = 1.0e-7;
        break;
    case SurfaceBenchmarkSource::MatlabLeo:
        criteria.id = "matlab-leo-v1";
        criteria.note = "MATLAB-compatible LEO route on top of legacy LEO references.";
        criteria.focus_segment = "negative_potential_tail";
        criteria.patch_rmse_v_max = 3.0e2;
        criteria.body_rmse_v_max = 5.5e3;
        criteria.body_je_rmse_a_per_m2_max = 1.0e-5;
        criteria.body_jnet_rmse_a_per_m2_max = 1.0e-5;
        criteria.negative_tail_body_je_rmse_a_per_m2_max = 8.0e-6;
        criteria.negative_tail_body_jnet_rmse_a_per_m2_max = 8.0e-6;
        break;
    default:
        criteria.id = "none";
        criteria.note = "No acceptance contract for this benchmark source.";
        criteria.focus_segment = "global";
        break;
    }
    return criteria;
}

bool resolveLegacyBenchmarkPaths(SurfaceBenchmarkSource source, LegacyBenchmarkPaths& paths)
{
    paths = LegacyBenchmarkPaths{};
    const auto absolute = [](const char* relative_path) {
        return legacyBenchmarkAbsolutePath(relative_path);
    };
    switch (source)
    {
    case SurfaceBenchmarkSource::CGeo:
        paths.input_path = absolute("ref/GEO/input_GEO_Charging.dat");
        paths.environment_path = absolute("ref/GEO/GEO_Environment.dat");
        paths.structure_material_path = absolute("ref/GEO/Parameter_Structure.dat");
        paths.dielectric_patch_material_path = absolute("ref/GEO/Parameter_Dielectric_Patch.dat");
        paths.metal_patch_material_path = absolute("ref/GEO/Parameter_Metal_Patch.dat");
        paths.patch_reference_curve_path = absolute("ref/GEO/reference_result_patch.dat");
        paths.body_reference_curve_path = absolute("ref/GEO/reference_result_body.dat");
        return true;
    case SurfaceBenchmarkSource::CLeoRam:
    case SurfaceBenchmarkSource::CLeoWake:
        paths.input_path = absolute("ref/LEO/input_LEO_Charging.dat");
        paths.environment_path = absolute("ref/LEO/LEO_Environment.dat");
        paths.structure_material_path = absolute("ref/LEO/Parameter_Structure.dat");
        paths.dielectric_patch_material_path = absolute("ref/LEO/Parameter_Dielectric_Patch.dat");
        paths.metal_patch_material_path = absolute("ref/LEO/Parameter_Metal_Patch.dat");
        paths.patch_reference_curve_path = absolute("ref/LEO/reference_result_patch.dat");
        paths.body_reference_curve_path = absolute("ref/LEO/reference_result_body.dat");
        return true;
    case SurfaceBenchmarkSource::MatlabGeo:
        paths.patch_reference_curve_path = absolute("ref/mat/reference_result_geo_patch.dat");
        paths.matlab_generator_path = absolute("ref/mat/generate_reference_result_geo.py");
        return true;
    case SurfaceBenchmarkSource::MatlabLeo:
        paths.input_path = absolute("ref/LEO/input_LEO_Charging.dat");
        paths.environment_path = absolute("ref/LEO/LEO_Environment.dat");
        paths.structure_material_path = absolute("ref/LEO/Parameter_Structure.dat");
        paths.dielectric_patch_material_path = absolute("ref/LEO/Parameter_Dielectric_Patch.dat");
        paths.metal_patch_material_path = absolute("ref/LEO/Parameter_Metal_Patch.dat");
        paths.patch_reference_curve_path = absolute("ref/LEO/reference_result_patch.dat");
        paths.body_reference_curve_path = absolute("ref/LEO/reference_result_body.dat");
        paths.matlab_generator_path =
            absolute("Toolkit/Surface Charging/scripts/reconstruct_legacy_leo_initial_body_current.py");
        return true;
    default:
        return false;
    }
}

LegacyBenchmarkInputConfig parseLegacyBenchmarkInputFile(const std::filesystem::path& path)
{
    LegacyBenchmarkInputConfig config;
    if (path.empty())
    {
        return config;
    }

    std::ifstream input(path);
    std::string line;
    while (std::getline(input, line))
    {
        if (!isNumericLine(line))
        {
            continue;
        }

        std::istringstream stream(line);
        double value = 0.0;
        if (stream >> value)
        {
            config.raw_values.push_back(value);
        }
    }
    return config;
}

std::vector<LegacyEnvironmentRecord>
parseLegacyEnvironmentFile(const std::filesystem::path& path)
{
    std::ifstream input(path);
    std::vector<LegacyEnvironmentRecord> records;
    std::string line;
    LegacyEnvironmentRecord* current = nullptr;
    bool expect_thermal_values = false;
    bool expect_electron_flux = false;
    bool expect_ion_flux = false;

    while (std::getline(input, line))
    {
        const std::string trimmed = Basic::trimAscii(line);
        if (trimmed.empty())
        {
            continue;
        }

        if (startsWith(trimmed, "No."))
        {
            const auto parts = splitWhitespace(trimmed);
            if (parts.size() >= 2)
            {
                try
                {
                    records.push_back(LegacyEnvironmentRecord{});
                    current = &records.back();
                    current->case_id =
                        static_cast<int>(std::max(0.0, std::stod(parts[1])));
                    const std::size_t id_pos = trimmed.find(parts[1]);
                    current->name = Basic::trimAscii(trimmed.substr(id_pos + parts[1].size()));
                    expect_thermal_values = false;
                    expect_electron_flux = false;
                    expect_ion_flux = false;
                }
                catch (const std::exception&)
                {
                    current = nullptr;
                }
            }
            continue;
        }

        if (current == nullptr)
        {
            continue;
        }

        if (trimmed.find("Ne1(1/cm^3)") != std::string::npos)
        {
            expect_thermal_values = true;
            continue;
        }

        if (expect_thermal_values && isNumericLine(trimmed))
        {
            current->thermal_values = parseNumbersFromLine(trimmed);
            expect_thermal_values = false;
            continue;
        }

        if (startsWith(trimmed, "Electron energy(eV):"))
        {
            current->electron_energy_ev =
                parseNumbersFromLine(trimmed.substr(std::string("Electron energy(eV):").size()));
            expect_electron_flux = true;
            continue;
        }

        if (expect_electron_flux && startsWith(trimmed, "Flux("))
        {
            const auto colon = trimmed.find(':');
            current->electron_flux_per_m2_s_sr_ev =
                colon == std::string::npos ? std::vector<double>{}
                                           : parseNumbersFromLine(trimmed.substr(colon + 1));
            expect_electron_flux = false;
            continue;
        }

        if (startsWith(trimmed, "Ion energy(eV):"))
        {
            current->ion_energy_ev =
                parseNumbersFromLine(trimmed.substr(std::string("Ion energy(eV):").size()));
            expect_ion_flux = true;
            continue;
        }

        if (expect_ion_flux && startsWith(trimmed, "Flux("))
        {
            const auto colon = trimmed.find(':');
            current->ion_flux_per_m2_s_sr_ev =
                colon == std::string::npos ? std::vector<double>{}
                                           : parseNumbersFromLine(trimmed.substr(colon + 1));
            expect_ion_flux = false;
            continue;
        }
    }

    return records;
}

std::vector<LegacyMaterialRecord>
parseLegacyMaterialTable(const std::filesystem::path& path, std::size_t trailing_numeric_count)
{
    std::ifstream input(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        lines.push_back(line);
    }

    std::vector<LegacyMaterialRecord> records;
    const std::string filename = path.filename().string();
    std::string section_marker;
    std::string alternate_section_marker;
    if (filename.find("Parameter_Structure") != std::string::npos)
    {
        section_marker = "Parameter Structure";
        alternate_section_marker = "Parameter_Structure";
    }
    else if (filename.find("Parameter_Metal") != std::string::npos)
    {
        section_marker = "Parameter Metal Patch";
        alternate_section_marker = "Parameter_Metal_Patch";
    }
    else if (filename.find("Parameter_Dielectric") != std::string::npos)
    {
        section_marker = "Parameter Dielectric Patch";
        alternate_section_marker = "Parameter_Dielectric_Patch";
    }

    std::size_t start_index = 0;
    bool found_section_marker = section_marker.empty();
    if (!section_marker.empty())
    {
        for (std::size_t i = 0; i < lines.size(); ++i)
        {
            const std::string trimmed = Basic::trimAscii(lines[i]);
            if (trimmed.find(section_marker) != std::string::npos ||
                (!alternate_section_marker.empty() &&
                 trimmed.find(alternate_section_marker) != std::string::npos))
            {
                found_section_marker = true;
                start_index = i + 1;
                break;
            }
        }
    }

    for (std::size_t line_index = start_index; line_index < lines.size(); ++line_index)
    {
        const std::string& current_line = lines[line_index];
        const std::string trimmed = Basic::trimAscii(current_line);
        if (trimmed.rfind("Result ", 0) == 0)
        {
            break;
        }
        if (found_section_marker &&
            trimmed.rfind("Parameter ", 0) == 0 &&
            trimmed.find(section_marker) == std::string::npos)
        {
            break;
        }
        if (!isNumericLine(current_line))
        {
            continue;
        }

        const auto parts = splitWhitespace(current_line);
        if (parts.size() < trailing_numeric_count + 2)
        {
            continue;
        }

        try
        {
            LegacyMaterialRecord record;
            record.material_id = static_cast<int>(std::max(0.0, std::stod(parts.front())));
            const std::size_t numeric_begin = parts.size() - trailing_numeric_count;
            record.name = joinNameTokens(parts, 1, numeric_begin);
            record.numeric_values.reserve(trailing_numeric_count);
            for (std::size_t i = numeric_begin; i < parts.size(); ++i)
            {
                record.numeric_values.push_back(std::stod(parts[i]));
            }
            records.push_back(record);
        }
        catch (const std::exception&)
        {
            continue;
        }
    }
    return records;
}

std::vector<LegacyBenchmarkCurveSample>
loadLegacyBenchmarkCurve(const std::filesystem::path& path, bool has_conduction_column)
{
    const auto resolved_path = legacyBenchmarkAbsolutePath(path);
    std::ifstream input(resolved_path);
    std::vector<LegacyBenchmarkCurveSample> samples;
    std::string line;
    while (std::getline(input, line))
    {
        if (!isNumericLine(line))
        {
            continue;
        }

        const auto parts = splitWhitespace(line);
        const std::size_t expected_count = has_conduction_column ? 11u : 10u;
        if (parts.size() < expected_count)
        {
            continue;
        }

        try
        {
            LegacyBenchmarkCurveSample sample;
            sample.cycle_index =
                static_cast<std::size_t>(std::max(0.0, std::stod(parts[0])));
            sample.potential_v = std::stod(parts[1]);
            sample.time_s = 1.0e-3 * std::stod(parts[2]);
            sample.jnet_a_per_m2 = 1.0e-9 * std::stod(parts[3]);
            sample.je_a_per_m2 = 1.0e-9 * std::stod(parts[4]);
            sample.jse_a_per_m2 = 1.0e-9 * std::stod(parts[5]);
            sample.jb_a_per_m2 = 1.0e-9 * std::stod(parts[6]);
            sample.ji_a_per_m2 = 1.0e-9 * std::stod(parts[7]);
            sample.jsi_a_per_m2 = 1.0e-9 * std::stod(parts[8]);
            sample.jph_a_per_m2 = 1.0e-9 * std::stod(parts[9]);
            if (has_conduction_column && parts.size() > 10)
            {
                sample.jcond_a_per_m2 = 1.0e-9 * std::stod(parts[10]);
            }
            samples.push_back(sample);
        }
        catch (const std::exception&)
        {
            continue;
        }
    }
    return samples;
}

double interpolateLegacyBenchmarkPotential(const std::vector<LegacyBenchmarkCurveSample>& samples,
                                          double time_s)
{
    if (samples.empty())
    {
        return 0.0;
    }
    if (time_s <= samples.front().time_s)
    {
        return samples.front().potential_v;
    }
    if (time_s >= samples.back().time_s)
    {
        return samples.back().potential_v;
    }

    for (std::size_t i = 1; i < samples.size(); ++i)
    {
        if (time_s <= samples[i].time_s)
        {
            const auto& a = samples[i - 1];
            const auto& b = samples[i];
            const double dt = std::max(1.0e-12, b.time_s - a.time_s);
            const double alpha = std::clamp((time_s - a.time_s) / dt, 0.0, 1.0);
            return a.potential_v + alpha * (b.potential_v - a.potential_v);
        }
    }

    return samples.back().potential_v;
}

template <typename Accessor>
double interpolateLegacyBenchmarkComponent(
    const std::vector<LegacyBenchmarkCurveSample>& samples,
    double time_s,
    const Accessor& accessor)
{
    if (samples.empty())
    {
        return 0.0;
    }
    if (time_s <= samples.front().time_s)
    {
        return accessor(samples.front());
    }
    if (time_s >= samples.back().time_s)
    {
        return accessor(samples.back());
    }

    for (std::size_t i = 1; i < samples.size(); ++i)
    {
        if (time_s <= samples[i].time_s)
        {
            const auto& a = samples[i - 1];
            const auto& b = samples[i];
            const double dt = std::max(1.0e-12, b.time_s - a.time_s);
            const double alpha = std::clamp((time_s - a.time_s) / dt, 0.0, 1.0);
            return accessor(a) + alpha * (accessor(b) - accessor(a));
        }
    }

    return accessor(samples.back());
}

double legacyBenchmarkTimeAlignmentToleranceS(double overlap_start_s, double overlap_end_s)
{
    const double span_s = std::abs(overlap_end_s - overlap_start_s);
    return std::max(1.0e-12, 1.0e-9 * std::max(1.0, span_s));
}

bool legacyBenchmarkHasMeaningfulTimeOverlap(const std::vector<LegacyBenchmarkCurveSample>& actual,
                                             const std::vector<LegacyBenchmarkCurveSample>& reference,
                                             double overlap_start_s,
                                             double overlap_end_s)
{
    if (actual.empty() || reference.empty())
    {
        return false;
    }

    const double overlap_span_s = overlap_end_s - overlap_start_s;
    if (!(overlap_span_s > 0.0))
    {
        return false;
    }

    const double actual_span_s = std::max(1.0e-12, actual.back().time_s - actual.front().time_s);
    const double reference_span_s =
        std::max(1.0e-12, reference.back().time_s - reference.front().time_s);
    const double actual_overlap_ratio = overlap_span_s / actual_span_s;
    const double reference_overlap_ratio = overlap_span_s / reference_span_s;

    return actual_overlap_ratio >= 0.25 && reference_overlap_ratio >= 0.25;
}

struct LegacyBenchmarkAlignedCoordinate
{
    double actual_time_s = 0.0;
    double reference_time_s = 0.0;
};

double legacyBenchmarkNormalizedProgress(double time_s, double start_s, double end_s)
{
    const double span_s = end_s - start_s;
    if (std::abs(span_s) <= 1.0e-12)
    {
        return 0.0;
    }
    return std::clamp((time_s - start_s) / span_s, 0.0, 1.0);
}

std::vector<LegacyBenchmarkAlignedCoordinate> buildLegacyBenchmarkAlignedCoordinates(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference)
{
    if (actual.empty() || reference.empty())
    {
        return {};
    }

    const double overlap_start_s = std::max(actual.front().time_s, reference.front().time_s);
    const double overlap_end_s = std::min(actual.back().time_s, reference.back().time_s);
    const double tolerance_s =
        legacyBenchmarkTimeAlignmentToleranceS(overlap_start_s, overlap_end_s);
    if (overlap_start_s <= overlap_end_s + tolerance_s &&
        legacyBenchmarkHasMeaningfulTimeOverlap(actual, reference, overlap_start_s, overlap_end_s))
    {
        std::vector<double> times;
        times.reserve(actual.size() + reference.size() + 2);
        const auto append_time = [&](double time_s) {
            if (time_s < overlap_start_s - tolerance_s || time_s > overlap_end_s + tolerance_s)
            {
                return;
            }
            times.push_back(std::clamp(time_s, overlap_start_s, overlap_end_s));
        };

        append_time(overlap_start_s);
        append_time(overlap_end_s);
        for (const auto& sample : actual)
        {
            append_time(sample.time_s);
        }
        for (const auto& sample : reference)
        {
            append_time(sample.time_s);
        }

        std::sort(times.begin(), times.end());
        std::vector<LegacyBenchmarkAlignedCoordinate> coordinates;
        coordinates.reserve(times.size());
        for (const double time_s : times)
        {
            if (coordinates.empty() ||
                std::abs(time_s - coordinates.back().actual_time_s) > tolerance_s)
            {
                coordinates.push_back({time_s, time_s});
            }
        }
        return coordinates;
    }

    const double actual_start_s = actual.front().time_s;
    const double actual_end_s = actual.back().time_s;
    const double reference_start_s = reference.front().time_s;
    const double reference_end_s = reference.back().time_s;
    constexpr double kNormalizedTolerance = 1.0e-12;

    std::vector<double> normalized_progresses;
    normalized_progresses.reserve(actual.size() + reference.size() + 2);
    normalized_progresses.push_back(0.0);
    normalized_progresses.push_back(1.0);
    for (const auto& sample : actual)
    {
        normalized_progresses.push_back(
            legacyBenchmarkNormalizedProgress(sample.time_s, actual_start_s, actual_end_s));
    }
    for (const auto& sample : reference)
    {
        normalized_progresses.push_back(
            legacyBenchmarkNormalizedProgress(sample.time_s, reference_start_s, reference_end_s));
    }

    std::sort(normalized_progresses.begin(), normalized_progresses.end());
    std::vector<LegacyBenchmarkAlignedCoordinate> coordinates;
    coordinates.reserve(normalized_progresses.size());
    double last_progress = std::numeric_limits<double>::quiet_NaN();
    for (const double progress : normalized_progresses)
    {
        if (std::isfinite(last_progress) &&
            std::abs(progress - last_progress) <= kNormalizedTolerance)
        {
            continue;
        }
        last_progress = progress;

        const double actual_time_s = actual_start_s + progress * (actual_end_s - actual_start_s);
        const double reference_time_s =
            reference_start_s + progress * (reference_end_s - reference_start_s);
        coordinates.push_back({actual_time_s, reference_time_s});
    }
    return coordinates;
}

LegacyBenchmarkCaseDefinition loadLegacyBenchmarkCaseDefinition(SurfaceBenchmarkSource source)
{
    LegacyBenchmarkCaseDefinition definition;
    if (!resolveLegacyBenchmarkPaths(source, definition.paths))
    {
        return definition;
    }

    definition.input = parseLegacyBenchmarkInputFile(definition.paths.input_path);
    definition.environments = parseLegacyEnvironmentFile(definition.paths.environment_path);
    definition.structure_materials =
        parseLegacyMaterialTable(definition.paths.structure_material_path, 11);
    definition.dielectric_patch_materials =
        parseLegacyMaterialTable(definition.paths.dielectric_patch_material_path, 13);
    definition.metal_patch_materials =
        parseLegacyMaterialTable(definition.paths.metal_patch_material_path, 11);
    definition.patch_reference_curve =
        loadLegacyBenchmarkCurve(definition.paths.patch_reference_curve_path, true);
    definition.body_reference_curve =
        loadLegacyBenchmarkCurve(definition.paths.body_reference_curve_path, false);
    return definition;
}

SurfaceChargingConfig applyLegacyBenchmarkExecutionConfig(
    const SurfaceChargingConfig& fallback_config,
    LegacyBenchmarkCaseDefinition* resolved_definition)
{
    SurfaceChargingConfig translated = fallback_config;
    const auto definition = loadLegacyBenchmarkCaseDefinition(fallback_config.benchmark_source);
    if (resolved_definition != nullptr)
    {
        *resolved_definition = definition;
    }
    if (definition.input.raw_values.empty())
    {
        return translated;
    }

    translated.runtime_route = SurfaceRuntimeRoute::LegacyBenchmark;
    translated.current_algorithm_mode = SurfaceCurrentAlgorithmMode::LegacyRefCompatible;
    translated.benchmark_mode = SurfaceBenchmarkMode::LegacyRefCompatible;
    translated.enable_pic_calibration = false;
    translated.enable_live_pic_window = false;
    translated.enable_live_pic_mcc = false;
    translated.enable_body_patch_circuit = false;
    translated.use_reference_current_balance = true;
    translated.surface_area_m2 = 1.0;
    translated.floating = true;
    translated.body_floating = false;
    translated.body_initial_potential_v = 0.0;
    translated.body_capacitance_f = 1.0;
    translated.max_abs_potential_v = 1.0e5;
    translated.max_delta_potential_v_per_step = 1.0e3;
    translated.internal_substeps = 1;
    translated.photoelectron_temperature_ev = 3.0;
    translated.electron_collection_coefficient = 1.0;
    translated.ion_collection_coefficient = 1.0;
    translated.patch_physics_overrides.clear();

    const bool is_geo = fallback_config.benchmark_source == SurfaceBenchmarkSource::CGeo;
    const auto& raw = definition.input.raw_values;
    const int jsun = raw.empty() ? 0 : static_cast<int>(std::lround(raw[0]));
    int jenv = 1;
    double theta_deg = 90.0;
    int jmb = 1;
    int jsee = 1;
    int jp = 1;
    int jmp = 1;
    double thickness_m = std::max(1.0e-6, fallback_config.dielectric_thickness_m);
    double layer_permittivity = std::max(1.0, fallback_config.material.getPermittivity());
    double layer_conductivity = std::max(0.0, fallback_config.material.getConductivity());
    int jram = 1;
    double altitude_km = 0.0;
    double flow_angle_deg = 0.0;

    if (is_geo)
    {
        if (raw.size() >= 8)
        {
            jenv = static_cast<int>(std::lround(raw[1]));
            theta_deg = raw[2];
            jmb = static_cast<int>(std::lround(raw[3]));
            jsee = static_cast<int>(std::lround(raw[4]));
            jp = static_cast<int>(std::lround(raw[5]));
            jmp = static_cast<int>(std::lround(raw[6]));
            thickness_m = std::max(1.0e-9, raw[7]);
        }
        if (jp == 2 && raw.size() >= 10)
        {
            layer_permittivity = std::max(1.0, raw[8]);
            layer_conductivity = std::max(0.0, raw[9]);
        }
        translated.regime = SurfaceChargingRegime::GeoKineticPicLike;
        translated.bulk_flow_velocity_m_per_s = 0.0;
        translated.flow_alignment_cosine = 0.0;
    }
    else
    {
        if (raw.size() >= 11)
        {
            jram = static_cast<int>(std::lround(raw[1]));
            altitude_km = raw[2];
            flow_angle_deg = raw[3];
            jenv = static_cast<int>(std::lround(raw[4]));
            theta_deg = raw[5];
            jmb = static_cast<int>(std::lround(raw[6]));
            jsee = static_cast<int>(std::lround(raw[7]));
            jp = static_cast<int>(std::lround(raw[8]));
            jmp = static_cast<int>(std::lround(raw[9]));
            thickness_m = std::max(1.0e-9, raw[10]);
        }
        if (jp == 2 && raw.size() >= 13)
        {
            layer_permittivity = std::max(1.0, raw[11]);
            layer_conductivity = std::max(0.0, raw[12]);
        }
        translated.regime = SurfaceChargingRegime::LeoFlowingPlasma;
        if (jram == 1)
        {
            translated.bulk_flow_velocity_m_per_s = 0.0;
            translated.flow_alignment_cosine = 0.0;
        }
        else
        {
            translated.bulk_flow_velocity_m_per_s = legacyLeoOrbitalSpeedMPerS(altitude_km);
            const double cosine = std::cos(flow_angle_deg * 3.14159265358979323846 / 180.0);
            translated.flow_alignment_cosine =
                jram == 3 ? -std::abs(cosine) : std::max(0.0, cosine);
        }
        translated.patch_flow_angle_deg = flow_angle_deg;
    }

    translated.patch_incidence_angle_deg = jsun == 0 ? 90.0 : theta_deg;
    translated.reference_see_model = mapLegacySeeModel(jsee);
    translated.dielectric_thickness_m = thickness_m;
    translated.derive_capacitance_from_material = false;

    const auto environment = findLegacyEnvironmentRecord(definition, jenv);
    if (environment.has_value())
    {
        applyLegacyEnvironmentToConfig(*environment, translated);
    }

    auto structure_material = findExactLegacyMaterialRecord(definition.structure_materials, jmb);
    if (!structure_material.has_value())
    {
        structure_material = findExactLegacyMaterialRecord(definition.metal_patch_materials, jmb);
    }
    if (!structure_material.has_value())
    {
        structure_material =
            findExactLegacyMaterialRecord(definition.dielectric_patch_materials, jmb);
    }
    if (!structure_material.has_value())
    {
        structure_material = findLegacyMaterialRecord(definition.structure_materials, jmb);
    }
    const auto patch_material_record =
        findLegacyMaterialRecord(jp == 1 ? definition.dielectric_patch_materials
                                         : definition.metal_patch_materials,
                                 jmp);

    if (structure_material.has_value())
    {
        translated.material = buildLegacyMaterial(*structure_material, false, false,
                                                  translated.material.getPermittivity(),
                                                  translated.material.getConductivity());
    }

    Material::MaterialProperty patch_material = translated.material;
    if (patch_material_record.has_value())
    {
        patch_material = buildLegacyMaterial(*patch_material_record, jp == 1, jp == 2,
                                             layer_permittivity, layer_conductivity);
    }

    if (jp == 1)
    {
        translated.capacitance_per_area_f_per_m2 =
            8.8541878128e-12 * std::max(1.0, patch_material.getPermittivity()) /
            std::max(1.0e-9, thickness_m);
    }
    else
    {
        translated.capacitance_per_area_f_per_m2 =
            8.8541878128e-12 * std::max(1.0, layer_permittivity) /
            std::max(1.0e-9, thickness_m);
        patch_material.setPermittivity(layer_permittivity);
        patch_material.setConductivity(layer_conductivity);
    }

    const double body_photo_current_density =
        jsun == 0 ? 0.0
                  : translated.material.getScalarProperty("legacy_photo_current_density_a_per_m2",
                                                          0.0);
    const double patch_photo_current_density =
        jsun == 0 ? 0.0
                  : patch_material.getScalarProperty("legacy_photo_current_density_a_per_m2", 0.0);
    translated.body_photo_current_density_a_per_m2 = body_photo_current_density;
    translated.patch_photo_current_density_a_per_m2 = patch_photo_current_density;

    SurfacePatchPhysicsConfig patch_override;
    patch_override.match_by_index = true;
    patch_override.node_index = 0;
    patch_override.override_material = true;
    patch_override.material = patch_material;
    patch_override.override_see_model = true;
    patch_override.reference_see_model = translated.reference_see_model;
    patch_override.override_patch_incidence_angle = true;
    patch_override.patch_incidence_angle_deg = translated.patch_incidence_angle_deg;
    patch_override.override_patch_flow_angle = true;
    patch_override.patch_flow_angle_deg = translated.patch_flow_angle_deg;
    patch_override.override_photoelectron_temperature = true;
    patch_override.photoelectron_temperature_ev = translated.photoelectron_temperature_ev;
    patch_override.override_patch_photo_current_density = true;
    patch_override.patch_photo_current_density_a_per_m2 = patch_photo_current_density;
    patch_override.override_electron_collection_coefficient = true;
    patch_override.electron_collection_coefficient = translated.electron_collection_coefficient;
    patch_override.override_ion_collection_coefficient = true;
    patch_override.ion_collection_coefficient = translated.ion_collection_coefficient;
    patch_override.override_plasma = true;
    patch_override.plasma = translated.plasma;
    patch_override.override_electron_spectrum = true;
    patch_override.electron_spectrum = translated.electron_spectrum;
    patch_override.has_electron_spectrum = translated.has_electron_spectrum;
    patch_override.override_ion_spectrum = true;
    patch_override.ion_spectrum = translated.ion_spectrum;
    patch_override.has_ion_spectrum = translated.has_ion_spectrum;
    patch_override.override_emission = true;
    patch_override.emission = translated.emission;
    translated.patch_physics_overrides.push_back(patch_override);

    translated.default_surface_physics.material = patch_material;
    translated.default_surface_physics.reference_see_model = translated.reference_see_model;
    translated.default_surface_physics.patch_incidence_angle_deg =
        translated.patch_incidence_angle_deg;
    translated.default_surface_physics.patch_flow_angle_deg = translated.patch_flow_angle_deg;
    translated.default_surface_physics.photoelectron_temperature_ev =
        translated.photoelectron_temperature_ev;
    translated.default_surface_physics.patch_photo_current_density_a_per_m2 =
        patch_photo_current_density;
    translated.default_surface_physics.electron_collection_coefficient =
        translated.electron_collection_coefficient;
    translated.default_surface_physics.ion_collection_coefficient =
        translated.ion_collection_coefficient;
    translated.default_surface_physics.plasma = translated.plasma;
    translated.default_surface_physics.electron_spectrum = translated.electron_spectrum;
    translated.default_surface_physics.ion_spectrum = translated.ion_spectrum;
    translated.default_surface_physics.has_electron_spectrum = translated.has_electron_spectrum;
    translated.default_surface_physics.has_ion_spectrum = translated.has_ion_spectrum;
    translated.default_surface_physics.emission = translated.emission;
    return translated;
}

double firstLegacyEquilibriumTimeS(const std::vector<LegacyBenchmarkCurveSample>& curve,
                                   double equilibrium_tolerance_v)
{
    if (curve.empty())
    {
        return 0.0;
    }

    const double terminal_potential = curve.back().potential_v;
    for (const auto& sample : curve)
    {
        if (std::abs(sample.potential_v - terminal_potential) <= equilibrium_tolerance_v)
        {
            return sample.time_s;
        }
    }
    return curve.back().time_s;
}

double legacyBenchmarkComponentRmse(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    const std::function<double(const LegacyBenchmarkCurveSample&)>& accessor)
{
    const auto aligned_coordinates =
        buildLegacyBenchmarkAlignedCoordinates(actual, reference);
    if (aligned_coordinates.empty())
    {
        return 0.0;
    }

    double error_sum = 0.0;
    for (const auto& coordinate : aligned_coordinates)
    {
        const double delta =
            interpolateLegacyBenchmarkComponent(actual, coordinate.actual_time_s, accessor) -
            interpolateLegacyBenchmarkComponent(reference, coordinate.reference_time_s, accessor);
        error_sum += delta * delta;
    }
    return std::sqrt(error_sum / static_cast<double>(aligned_coordinates.size()));
}

double legacyBenchmarkComponentMeanSignedDelta(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    const std::function<double(const LegacyBenchmarkCurveSample&)>& accessor)
{
    const auto aligned_coordinates =
        buildLegacyBenchmarkAlignedCoordinates(actual, reference);
    if (aligned_coordinates.empty())
    {
        return 0.0;
    }

    double delta_sum = 0.0;
    for (const auto& coordinate : aligned_coordinates)
    {
        delta_sum +=
            interpolateLegacyBenchmarkComponent(actual, coordinate.actual_time_s, accessor) -
            interpolateLegacyBenchmarkComponent(reference, coordinate.reference_time_s, accessor);
    }
    return delta_sum / static_cast<double>(aligned_coordinates.size());
}

double estimateInitialEFoldingEnergyEv(
    const std::vector<LegacyBenchmarkCurveSample>& samples,
    const std::function<double(const LegacyBenchmarkCurveSample&)>& accessor)
{
    if (samples.size() < 2)
    {
        return 0.0;
    }

    const auto& base = samples.front();
    const double base_current = std::abs(accessor(base));
    if (base_current <= 1.0e-30)
    {
        return 0.0;
    }

    for (std::size_t i = 1; i < samples.size(); ++i)
    {
        const auto& sample = samples[i];
        const double delta_v = std::abs(sample.potential_v - base.potential_v);
        if (delta_v <= 0.0 || sample.potential_v >= base.potential_v)
        {
            continue;
        }

        const double current = std::abs(accessor(sample));
        if (current <= 1.0e-30 || current >= base_current)
        {
            continue;
        }

        const double ratio = current / base_current;
        if (ratio <= 0.0 || ratio >= 1.0)
        {
            continue;
        }
        return delta_v / std::abs(std::log(ratio));
    }

    return 0.0;
}

double initialComponentDelta(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    const std::function<double(const LegacyBenchmarkCurveSample&)>& accessor)
{
    if (actual.empty() || reference.empty())
    {
        return 0.0;
    }
    return accessor(actual.front()) - accessor(reference.front());
}

double initialComponentRatio(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    const std::function<double(const LegacyBenchmarkCurveSample&)>& accessor)
{
    if (actual.empty() || reference.empty())
    {
        return 0.0;
    }
    const double denominator = accessor(reference.front());
    if (std::abs(denominator) <= 1.0e-30)
    {
        return 0.0;
    }
    return accessor(actual.front()) / denominator;
}

struct PotentialTailMetric
{
    bool valid = false;
    std::size_t sample_count = 0;
    double threshold_v = 0.0;
    double rmse = 0.0;
};

PotentialTailMetric legacyBenchmarkNegativeTailRmse(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    const std::function<double(const LegacyBenchmarkCurveSample&)>& accessor)
{
    PotentialTailMetric metric;
    const auto aligned_coordinates =
        buildLegacyBenchmarkAlignedCoordinates(actual, reference);
    if (aligned_coordinates.empty())
    {
        return metric;
    }

    double min_reference_potential_v = std::numeric_limits<double>::infinity();
    for (const auto& coordinate : aligned_coordinates)
    {
        min_reference_potential_v =
            std::min(min_reference_potential_v,
                     interpolateLegacyBenchmarkPotential(reference,
                                                         coordinate.reference_time_s));
    }
    if (!(min_reference_potential_v < 0.0))
    {
        return metric;
    }

    metric.threshold_v = std::min(-1.0, 0.5 * min_reference_potential_v);
    double error_sum = 0.0;
    for (const auto& coordinate : aligned_coordinates)
    {
        if (interpolateLegacyBenchmarkPotential(reference, coordinate.reference_time_s) >
            metric.threshold_v)
        {
            continue;
        }
        const double delta =
            interpolateLegacyBenchmarkComponent(actual, coordinate.actual_time_s, accessor) -
            interpolateLegacyBenchmarkComponent(reference, coordinate.reference_time_s, accessor);
        error_sum += delta * delta;
        metric.sample_count += 1;
    }
    if (metric.sample_count == 0)
    {
        return metric;
    }

    metric.valid = true;
    metric.rmse = std::sqrt(error_sum / static_cast<double>(metric.sample_count));
    return metric;
}

std::vector<double> legacySmoothLine3Values(const std::vector<double>& x, const std::vector<double>& y)
{
    const std::size_t count = std::min(x.size(), y.size());
    if (count <= 2)
    {
        return std::vector<double>(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(count));
    }

    std::vector<double> smoothed(count, 0.0);
    for (std::size_t i = 0; i < count; ++i)
    {
        std::size_t first = 0;
        if (i == 0)
        {
            first = 0;
        }
        else if (i >= count - 1)
        {
            first = count - 3;
        }
        else
        {
            first = i - 1;
        }

        const auto i0 = first;
        const auto i1 = first + 1;
        const auto i2 = first + 2;
        const double t0 = x[i0];
        const double t1 = x[i1];
        const double t2 = x[i2];
        const double p0 = y[i0];
        const double p1 = y[i1];
        const double p2 = y[i2];
        const double det = (t0 - t1) * (t0 - t2) * (t1 - t2);
        if (std::abs(det) <= 1.0e-30)
        {
            smoothed[i] = y[i];
            continue;
        }
        const double a0 =
            (p0 * t1 * t2 * (t1 - t2) + p1 * t2 * t0 * (t2 - t0) + p2 * t0 * t1 * (t0 - t1)) /
            det;
        const double a1 =
            (p0 * (t2 * t2 - t1 * t1) + p1 * (t0 * t0 - t2 * t2) + p2 * (t1 * t1 - t0 * t0)) /
            det;
        const double a2 =
            (p0 * (t1 - t2) + p1 * (t2 - t0) + p2 * (t0 - t1)) / det;
        smoothed[i] = a0 + a1 * x[i] + a2 * x[i] * x[i];
    }
    return smoothed;
}

double legacyInterpolateThreePoint(const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   double target)
{
    const std::size_t count = std::min(x.size(), y.size());
    if (count == 0)
    {
        return 0.0;
    }
    if (count == 1)
    {
        return y.front();
    }
    if (count == 2)
    {
        const double x0 = x[0];
        const double x1 = x[1];
        if (std::abs(x1 - x0) <= 1.0e-30)
        {
            return y[0];
        }
        return (y[0] * (target - x1) - y[1] * (target - x0)) / (x0 - x1);
    }

    std::ptrdiff_t k = 0;
    std::ptrdiff_t m = 0;
    if (target <= x[1])
    {
        k = 0;
        m = 2;
    }
    else if (target >= x[count - 2])
    {
        k = static_cast<std::ptrdiff_t>(count) - 3;
        m = static_cast<std::ptrdiff_t>(count) - 1;
    }
    else
    {
        k = 1;
        m = static_cast<std::ptrdiff_t>(count);
        while (m - k != 1)
        {
            const auto i = (k + m) / 2;
            if (target < x[static_cast<std::size_t>(i - 1)])
            {
                m = i;
            }
            else
            {
                k = i;
            }
        }
        k -= 1;
        m -= 1;
        if (std::abs(target - x[static_cast<std::size_t>(k)]) <
            std::abs(target - x[static_cast<std::size_t>(m)]))
        {
            k -= 1;
        }
        else
        {
            m += 1;
        }
    }

    double value = 0.0;
    for (auto i = k; i <= m; ++i)
    {
        double weight = 1.0;
        for (auto j = k; j <= m; ++j)
        {
            if (j != i)
            {
                weight *=
                    (target - x[static_cast<std::size_t>(j)]) /
                    (x[static_cast<std::size_t>(i)] - x[static_cast<std::size_t>(j)]);
            }
        }
        value += weight * y[static_cast<std::size_t>(i)];
    }
    return value;
}

double reconstructLegacyLeoInitialBodyJeApm2(const LegacyEnvironmentRecord& environment)
{
    constexpr double kLegacyPi = 3.14159265358979323846;
    constexpr double kLegacyElementaryCharge = 1.60217733e-19;
    if (environment.thermal_values.size() < 15 || environment.electron_energy_ev.empty() ||
        environment.electron_flux_per_m2_s_sr_ev.empty() || environment.ion_energy_ev.empty())
    {
        return 0.0;
    }

    const std::vector<double> ne = {environment.thermal_values[0], environment.thermal_values[2],
                                    environment.thermal_values[4]};
    const std::vector<double> te = {environment.thermal_values[1], environment.thermal_values[3],
                                    environment.thermal_values[5]};
    const std::vector<double> ti = {environment.thermal_values[7], environment.thermal_values[10],
                                    environment.thermal_values[13]};

    std::vector<double> electron_energy_ev;
    std::vector<double> ion_energy_ev;
    std::vector<double> electron_flux_per_m2_s_sr_ev;
    for (std::size_t i = 0; i < environment.electron_energy_ev.size() &&
                            i < environment.electron_flux_per_m2_s_sr_ev.size();
         ++i)
    {
        if (environment.electron_energy_ev[i] > 0.0)
        {
            electron_energy_ev.push_back(environment.electron_energy_ev[i]);
            electron_flux_per_m2_s_sr_ev.push_back(environment.electron_flux_per_m2_s_sr_ev[i]);
        }
    }
    for (double energy_ev : environment.ion_energy_ev)
    {
        if (energy_ev > 0.0)
        {
            ion_energy_ev.push_back(energy_ev);
        }
    }
    if (electron_energy_ev.empty() || ion_energy_ev.empty())
    {
        return 0.0;
    }

    const auto smoothed_flux =
        legacySmoothLine3Values(electron_energy_ev, electron_flux_per_m2_s_sr_ev);
    constexpr int kLegacyLeoSamples = 1000;
    double minimum_energy_ev = std::min(electron_energy_ev.front(), ion_energy_ev.front());
    double maximum_energy_ev = std::max(electron_energy_ev.back(), ion_energy_ev.back());
    for (std::size_t i = 0; i < te.size(); ++i)
    {
        if (te[i] > 0.0 && ti[i] > 0.0)
        {
            minimum_energy_ev = std::min(minimum_energy_ev, 0.05 * std::min(te[i], ti[i]));
            maximum_energy_ev = std::max(maximum_energy_ev, 10.0 * std::max(te[i], ti[i]));
        }
    }
    const double nmin = std::log10(std::max(1.0e-6, minimum_energy_ev));
    const double nmax = std::log10(std::max(10.0 * minimum_energy_ev, maximum_energy_ev));
    const double dn = (nmax - nmin) / static_cast<double>(kLegacyLeoSamples);

    double xe = 0.0;
    for (int i = 0; i < kLegacyLeoSamples; ++i)
    {
        const double energy_ev = std::pow(10.0, nmin + static_cast<double>(i) * dn);
        const double next_energy_ev =
            std::pow(10.0, nmin + static_cast<double>(i + 1) * dn);
        const double width_ev = next_energy_ev - energy_ev;
        double fluxe0 = 0.0;
        if (energy_ev < electron_energy_ev.front())
        {
            for (std::size_t k = 0; k < te.size(); ++k)
            {
                if (te[k] > 0.0)
                {
                    fluxe0 += 5.325e10 * ne[k] * energy_ev * std::pow(te[k], -1.5) *
                              std::exp(-energy_ev / te[k]);
                }
            }
        }
        else if (energy_ev <= electron_energy_ev.back())
        {
            fluxe0 = std::max(0.0,
                              legacyInterpolateThreePoint(electron_energy_ev, smoothed_flux,
                                                          energy_ev));
        }
        else
        {
            for (std::size_t k = 0; k < te.size(); ++k)
            {
                if (te[k] > 0.0)
                {
                    fluxe0 += 5.325e10 * ne[k] * energy_ev * std::pow(te[k], -1.5) *
                              std::exp(-energy_ev / te[k]);
                }
            }
        }
        xe += fluxe0 * width_ev;
    }
    return -1.0e9 * kLegacyPi * kLegacyElementaryCharge * xe * 1.0e-9;
}

LegacyBenchmarkMetrics computeLegacyBenchmarkMetrics(
    const std::vector<LegacyBenchmarkCurveSample>& actual,
    const std::vector<LegacyBenchmarkCurveSample>& reference,
    double equilibrium_tolerance_v)
{
    LegacyBenchmarkMetrics metrics;
    const auto aligned_coordinates =
        buildLegacyBenchmarkAlignedCoordinates(actual, reference);
    if (aligned_coordinates.empty())
    {
        return metrics;
    }

    double error_sum = 0.0;
    for (const auto& coordinate : aligned_coordinates)
    {
        const double delta =
            interpolateLegacyBenchmarkPotential(actual, coordinate.actual_time_s) -
            interpolateLegacyBenchmarkPotential(reference, coordinate.reference_time_s);
        error_sum += delta * delta;
    }

    const auto& terminal_coordinate = aligned_coordinates.back();
    metrics.valid = true;
    metrics.compared_sample_count = aligned_coordinates.size();
    metrics.rmse_v = std::sqrt(error_sum / static_cast<double>(aligned_coordinates.size()));
    metrics.terminal_potential_delta_v =
        interpolateLegacyBenchmarkPotential(actual, terminal_coordinate.actual_time_s) -
        interpolateLegacyBenchmarkPotential(reference, terminal_coordinate.reference_time_s);
    metrics.terminal_time_delta_s = actual.back().time_s - reference.back().time_s;
    metrics.time_to_equilibrium_delta_s =
        firstLegacyEquilibriumTimeS(actual, equilibrium_tolerance_v) -
        firstLegacyEquilibriumTimeS(reference, equilibrium_tolerance_v);
    return metrics;
}

LegacyBenchmarkBodyTailMetrics computeLegacyBenchmarkBodyTailMetrics(
    const std::vector<LegacyBenchmarkCurveSample>& actual_body,
    const std::vector<LegacyBenchmarkCurveSample>& reference_body)
{
    LegacyBenchmarkBodyTailMetrics metrics;
    const auto je_tail = legacyBenchmarkNegativeTailRmse(
        actual_body, reference_body,
        [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; });
    const auto jnet_tail = legacyBenchmarkNegativeTailRmse(
        actual_body, reference_body,
        [](const LegacyBenchmarkCurveSample& sample) { return sample.jnet_a_per_m2; });

    metrics.valid = je_tail.valid || jnet_tail.valid;
    metrics.sample_count = std::max(je_tail.sample_count, jnet_tail.sample_count);
    metrics.threshold_v = je_tail.valid ? je_tail.threshold_v : jnet_tail.threshold_v;
    metrics.je_rmse_a_per_m2 = je_tail.rmse;
    metrics.jnet_rmse_a_per_m2 = jnet_tail.rmse;
    return metrics;
}

LegacyBenchmarkAcceptanceGateResult evaluateLegacyBenchmarkAcceptanceGate(
    const LegacyBenchmarkAcceptanceCriteria& criteria,
    const LegacyBenchmarkMetrics& patch_metrics,
    const LegacyBenchmarkMetrics& body_metrics,
    double body_je_rmse_a_per_m2,
    double body_jnet_rmse_a_per_m2,
    const LegacyBenchmarkBodyTailMetrics& body_tail_metrics)
{
    LegacyBenchmarkAcceptanceGateResult result;
    result.applicable = criteria.id != "none";
    if (!result.applicable)
    {
        result.status = "NOT_APPLICABLE";
        return result;
    }

    std::vector<std::string> failed_keys;
    std::vector<std::string> failure_details;
    const auto evaluate = [&](const char* key, double actual, double max_allowed) {
        if (!(max_allowed > 0.0) || !std::isfinite(max_allowed))
        {
            return;
        }
        result.checks_total += 1;
        if (std::isfinite(actual) && actual <= max_allowed)
        {
            return;
        }

        failed_keys.emplace_back(key);
        std::ostringstream detail;
        detail << key << "(actual=" << actual << ", max=" << max_allowed << ")";
        failure_details.push_back(detail.str());
    };

    evaluate("patch_rmse_v", patch_metrics.rmse_v, criteria.patch_rmse_v_max);
    evaluate("body_rmse_v", body_metrics.rmse_v, criteria.body_rmse_v_max);
    evaluate("body_je_rmse_a_per_m2", body_je_rmse_a_per_m2,
             criteria.body_je_rmse_a_per_m2_max);
    evaluate("body_jnet_rmse_a_per_m2", body_jnet_rmse_a_per_m2,
             criteria.body_jnet_rmse_a_per_m2_max);

    if (body_tail_metrics.sample_count > 0)
    {
        evaluate("body_negative_tail_je_rmse_a_per_m2", body_tail_metrics.je_rmse_a_per_m2,
                 criteria.negative_tail_body_je_rmse_a_per_m2_max);
        evaluate("body_negative_tail_jnet_rmse_a_per_m2", body_tail_metrics.jnet_rmse_a_per_m2,
                 criteria.negative_tail_body_jnet_rmse_a_per_m2_max);
    }

    result.checks_failed = failed_keys.size();
    result.pass = result.checks_failed == 0;
    result.status = result.pass ? "PASS" : "FAIL";
    result.failed_metric_keys = failed_keys.empty() ? "none" : joinStrings(failed_keys, ",");
    result.failure_details =
        failure_details.empty() ? "none" : joinStrings(failure_details, "; ");
    return result;
}

LegacyBenchmarkConsistencyDiagnostics analyzeLegacyBenchmarkConsistency(
    SurfaceBenchmarkSource source, const LegacyBenchmarkCaseDefinition& definition)
{
    LegacyBenchmarkConsistencyDiagnostics diagnostics;
    diagnostics.valid = true;
    diagnostics.authority = "ReferenceCurve";
    diagnostics.status = "NotApplicable";

    if (source != SurfaceBenchmarkSource::CLeoRam && source != SurfaceBenchmarkSource::CLeoWake)
    {
        return diagnostics;
    }

    diagnostics.input_reconstruction_supported = true;
    diagnostics.status = "InputReconstructionAvailable";
    if (definition.input.raw_values.size() <= 4 || definition.body_reference_curve.empty())
    {
        diagnostics.status = "MissingInputOrReference";
        return diagnostics;
    }

    const int environment_id = static_cast<int>(std::lround(definition.input.raw_values[4]));
    const auto environment_it =
        std::find_if(definition.environments.begin(), definition.environments.end(),
                     [&](const LegacyEnvironmentRecord& record) {
                         return record.case_id == environment_id;
                     });
    if (environment_it == definition.environments.end())
    {
        diagnostics.status = "EnvironmentRecordMissing";
        return diagnostics;
    }

    diagnostics.reconstructed_body_initial_je_a_per_m2 =
        reconstructLegacyLeoInitialBodyJeApm2(*environment_it);
    diagnostics.reference_body_initial_je_a_per_m2 =
        definition.body_reference_curve.front().je_a_per_m2;
    diagnostics.body_initial_je_delta_a_per_m2 =
        diagnostics.reconstructed_body_initial_je_a_per_m2 -
        diagnostics.reference_body_initial_je_a_per_m2;
    const double denominator = diagnostics.reference_body_initial_je_a_per_m2;
    diagnostics.body_initial_je_ratio =
        std::abs(denominator) <= 1.0e-30
            ? 0.0
            : diagnostics.reconstructed_body_initial_je_a_per_m2 / denominator;

    const bool ratio_consistent =
        diagnostics.body_initial_je_ratio > 0.8 && diagnostics.body_initial_je_ratio < 1.25;
    const bool delta_consistent =
        std::abs(diagnostics.body_initial_je_delta_a_per_m2) <= 1.0e-6;
    diagnostics.input_reference_consistent = ratio_consistent || delta_consistent;
    diagnostics.status = diagnostics.input_reference_consistent ? "InputReferenceConsistent"
                                                                : "InputReferenceMismatch";
    diagnostics.authority =
        diagnostics.input_reference_consistent ? "InputAndReferenceAgree" : "ReferenceCurve";
    return diagnostics;
}

bool writeLegacyBenchmarkReport(const std::filesystem::path& report_path,
                                SurfaceRuntimeRoute runtime_route,
                                SurfaceBenchmarkSource source,
                                LegacyBenchmarkExecutionMode execution_mode,
                                const LegacyBenchmarkCaseDefinition& definition,
                                const std::vector<LegacyBenchmarkCurveSample>& actual_patch,
                                const std::vector<LegacyBenchmarkCurveSample>& actual_body,
                                const LegacyBenchmarkMetrics& patch_metrics,
                                const LegacyBenchmarkMetrics& body_metrics,
                                const LegacyBenchmarkConsistencyDiagnostics& consistency)
{
    std::ofstream output(report_path);
    if (!output.is_open())
    {
        return false;
    }

    const auto acceptance = legacyBenchmarkAcceptanceCriteria(source);
    const auto body_tail_metrics =
        computeLegacyBenchmarkBodyTailMetrics(actual_body,
                                              definition.body_reference_curve);
    const double body_je_rmse_a_per_m2 =
        legacyBenchmarkComponentRmse(actual_body, definition.body_reference_curve,
                                     [](const LegacyBenchmarkCurveSample& sample) {
                                         return sample.je_a_per_m2;
                                     });
    const double body_jnet_rmse_a_per_m2 =
        legacyBenchmarkComponentRmse(actual_body, definition.body_reference_curve,
                                     [](const LegacyBenchmarkCurveSample& sample) {
                                         return sample.jnet_a_per_m2;
                                     });
    const auto acceptance_gate =
        evaluateLegacyBenchmarkAcceptanceGate(acceptance, patch_metrics, body_metrics,
                                              body_je_rmse_a_per_m2,
                                              body_jnet_rmse_a_per_m2,
                                              body_tail_metrics);

    output << std::fixed << std::setprecision(12);
    output << "runtime_route="
           << (runtime_route == SurfaceRuntimeRoute::LegacyBenchmark ? "LegacyBenchmark"
                                                                    : "SCDATUnified")
           << "\n";
    output << "benchmark_source=";
    switch (source)
    {
    case SurfaceBenchmarkSource::CGeo:
        output << "C-GEO";
        break;
    case SurfaceBenchmarkSource::CLeoRam:
        output << "C-LEO-RAM";
        break;
    case SurfaceBenchmarkSource::CLeoWake:
        output << "C-LEO-WAKE";
        break;
    case SurfaceBenchmarkSource::MatlabGeo:
        output << "MATLAB-GEO";
        break;
    case SurfaceBenchmarkSource::MatlabLeo:
        output << "MATLAB-LEO";
        break;
    default:
        output << "None";
        break;
    }
    output << "\n";
    output << "execution_mode=" << legacyBenchmarkExecutionModeName(execution_mode) << "\n";
    output << "baseline_family=" << legacyBenchmarkBaselineFamilyName(source) << "\n";
    output << "baseline_origin=" << legacyBenchmarkBaselineOriginName(source) << "\n";
        output << "acceptance_contract_version=legacy-benchmark-v1\n";
        output << "acceptance_contract_id=" << acceptance.id << "\n";
        output << "acceptance_focus_segment=" << acceptance.focus_segment << "\n";
        output << "acceptance_note=" << acceptance.note << "\n";
        output << "acceptance_patch_rmse_v_max=" << acceptance.patch_rmse_v_max << "\n";
        output << "acceptance_body_rmse_v_max=" << acceptance.body_rmse_v_max << "\n";
        output << "acceptance_body_je_rmse_a_per_m2_max="
            << acceptance.body_je_rmse_a_per_m2_max << "\n";
        output << "acceptance_body_jnet_rmse_a_per_m2_max="
            << acceptance.body_jnet_rmse_a_per_m2_max << "\n";
        output << "acceptance_negative_tail_body_je_rmse_a_per_m2_max="
            << acceptance.negative_tail_body_je_rmse_a_per_m2_max << "\n";
        output << "acceptance_negative_tail_body_jnet_rmse_a_per_m2_max="
            << acceptance.negative_tail_body_jnet_rmse_a_per_m2_max << "\n";
            output << "acceptance_gate_applicable=" << (acceptance_gate.applicable ? 1 : 0) << "\n";
            output << "acceptance_gate_status=" << acceptance_gate.status << "\n";
            output << "acceptance_gate_pass=" << (acceptance_gate.pass ? 1 : 0) << "\n";
            output << "acceptance_gate_checks_total=" << acceptance_gate.checks_total << "\n";
            output << "acceptance_gate_checks_failed=" << acceptance_gate.checks_failed << "\n";
            output << "acceptance_gate_failed_metric_keys=" << acceptance_gate.failed_metric_keys
               << "\n";
            output << "acceptance_gate_failure_details=" << acceptance_gate.failure_details << "\n";
    output << "input_path=" << definition.paths.input_path.string() << "\n";
    output << "environment_path=" << definition.paths.environment_path.string() << "\n";
    output << "structure_material_path=" << definition.paths.structure_material_path.string()
           << "\n";
    output << "patch_reference_curve=" << definition.paths.patch_reference_curve_path.string()
           << "\n";
    output << "body_reference_curve=" << definition.paths.body_reference_curve_path.string()
           << "\n";
    if (!definition.paths.matlab_generator_path.empty())
    {
        output << "matlab_generator_path=" << definition.paths.matlab_generator_path.string()
               << "\n";
    }
    output << "reference_patch_sample_count=" << definition.patch_reference_curve.size() << "\n";
    output << "reference_body_sample_count=" << definition.body_reference_curve.size() << "\n";
    output << "consistency_status=" << consistency.status << "\n";
    output << "consistency_authority=" << consistency.authority << "\n";
    output << "input_reconstruction_supported="
           << (consistency.input_reconstruction_supported ? 1 : 0) << "\n";
    output << "input_reference_consistent=" << (consistency.input_reference_consistent ? 1 : 0)
           << "\n";
    output << "reconstructed_body_initial_je_a_per_m2="
           << consistency.reconstructed_body_initial_je_a_per_m2 << "\n";
    output << "reference_body_initial_je_a_per_m2="
           << consistency.reference_body_initial_je_a_per_m2 << "\n";
    output << "body_initial_je_input_reference_delta_a_per_m2="
           << consistency.body_initial_je_delta_a_per_m2 << "\n";
    output << "body_initial_je_input_reference_ratio="
           << consistency.body_initial_je_ratio << "\n";
    output << "patch_metrics_valid=" << (patch_metrics.valid ? 1 : 0) << "\n";
    output << "patch_compared_sample_count=" << patch_metrics.compared_sample_count << "\n";
    output << "patch_rmse_v=" << patch_metrics.rmse_v << "\n";
    output << "patch_terminal_potential_delta_v=" << patch_metrics.terminal_potential_delta_v
           << "\n";
    output << "patch_terminal_time_delta_s=" << patch_metrics.terminal_time_delta_s << "\n";
    output << "patch_time_to_equilibrium_delta_s=" << patch_metrics.time_to_equilibrium_delta_s
           << "\n";
    output << "body_metrics_valid=" << (body_metrics.valid ? 1 : 0) << "\n";
    output << "body_compared_sample_count=" << body_metrics.compared_sample_count << "\n";
    output << "body_rmse_v=" << body_metrics.rmse_v << "\n";
    output << "body_terminal_potential_delta_v=" << body_metrics.terminal_potential_delta_v
           << "\n";
    output << "body_terminal_time_delta_s=" << body_metrics.terminal_time_delta_s << "\n";
    output << "body_time_to_equilibrium_delta_s=" << body_metrics.time_to_equilibrium_delta_s
           << "\n";
            output << "body_negative_tail_threshold_v=" << body_tail_metrics.threshold_v << "\n";
            output << "body_negative_tail_sample_count=" << body_tail_metrics.sample_count << "\n";
            output << "body_negative_tail_je_rmse_a_per_m2=" << body_tail_metrics.je_rmse_a_per_m2
                << "\n";
            output << "body_negative_tail_jnet_rmse_a_per_m2="
                << body_tail_metrics.jnet_rmse_a_per_m2
            << "\n";
    const auto write_component_rmse = [&](const char* prefix,
                                          const std::vector<LegacyBenchmarkCurveSample>& actual,
                                          const std::vector<LegacyBenchmarkCurveSample>& reference) {
        output << prefix << "_jnet_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.jnet_a_per_m2; })
               << "\n";
        output << prefix << "_je_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; })
               << "\n";
        output << prefix << "_jse_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.jse_a_per_m2; })
               << "\n";
        output << prefix << "_jb_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.jb_a_per_m2; })
               << "\n";
        output << prefix << "_ji_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.ji_a_per_m2; })
               << "\n";
        output << prefix << "_jsi_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.jsi_a_per_m2; })
               << "\n";
        output << prefix << "_jph_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.jph_a_per_m2; })
               << "\n";
        output << prefix << "_jcond_rmse_a_per_m2="
               << legacyBenchmarkComponentRmse(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.jcond_a_per_m2; })
               << "\n";
        output << prefix << "_je_mean_signed_delta_a_per_m2="
               << legacyBenchmarkComponentMeanSignedDelta(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; })
               << "\n";
        output << prefix << "_ji_mean_signed_delta_a_per_m2="
               << legacyBenchmarkComponentMeanSignedDelta(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.ji_a_per_m2; })
               << "\n";
        output << prefix << "_je_initial_efolding_energy_ev_actual="
               << estimateInitialEFoldingEnergyEv(
                      actual,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; })
               << "\n";
        output << prefix << "_je_initial_efolding_energy_ev_reference="
               << estimateInitialEFoldingEnergyEv(
                      reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; })
               << "\n";
        output << prefix << "_je_initial_delta_a_per_m2="
               << initialComponentDelta(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; })
               << "\n";
        output << prefix << "_je_initial_ratio_actual_to_reference="
               << initialComponentRatio(
                      actual, reference,
                      [](const LegacyBenchmarkCurveSample& sample) { return sample.je_a_per_m2; })
               << "\n";
    };
    write_component_rmse("patch", actual_patch, definition.patch_reference_curve);
    write_component_rmse("body", actual_body, definition.body_reference_curve);
    const auto write_preview = [&](const char* label,
                                   const std::vector<LegacyBenchmarkCurveSample>& actual,
                                   const std::vector<LegacyBenchmarkCurveSample>& reference) {
        const std::size_t preview_count =
            std::min<std::size_t>(12, std::max(actual.size(), reference.size()));
        output << label << "_preview_count=" << preview_count << "\n";
        for (std::size_t i = 0; i < preview_count; ++i)
        {
            const LegacyBenchmarkCurveSample empty{};
            const auto& a = i < actual.size() ? actual[i] : empty;
            const auto& r = i < reference.size() ? reference[i] : empty;
            output << label << "_sample_" << i << "="
                   << "actual_cycle:" << a.cycle_index
                   << ",reference_cycle:" << r.cycle_index
                   << ",actual_time_s:" << a.time_s
                   << ",reference_time_s:" << r.time_s
                   << ",actual_potential_v:" << a.potential_v
                   << ",reference_potential_v:" << r.potential_v
                   << ",actual_jnet_a_per_m2:" << a.jnet_a_per_m2
                   << ",reference_jnet_a_per_m2:" << r.jnet_a_per_m2
                   << ",actual_je_a_per_m2:" << a.je_a_per_m2
                   << ",reference_je_a_per_m2:" << r.je_a_per_m2
                   << ",actual_ji_a_per_m2:" << a.ji_a_per_m2
                   << ",reference_ji_a_per_m2:" << r.ji_a_per_m2
                   << "\n";
        }
    };
    write_preview("patch", actual_patch, definition.patch_reference_curve);
    write_preview("body", actual_body, definition.body_reference_curve);
    return true;
}

bool writeLegacyBenchmarkComparisonCsv(
    const std::filesystem::path& csv_path,
    const std::vector<LegacyBenchmarkCurveSample>& actual_patch,
    const std::vector<LegacyBenchmarkCurveSample>& reference_patch,
    const std::vector<LegacyBenchmarkCurveSample>& actual_body,
    const std::vector<LegacyBenchmarkCurveSample>& reference_body)
{
    std::ofstream output(csv_path);
    if (!output.is_open())
    {
        return false;
    }

    output << std::fixed << std::setprecision(12);
    output << "role,index,actual_cycle,reference_cycle,actual_time_s,reference_time_s,"
              "actual_potential_v,reference_potential_v,delta_potential_v,"
              "actual_jnet_a_per_m2,reference_jnet_a_per_m2,"
              "actual_je_a_per_m2,reference_je_a_per_m2,"
              "actual_jse_a_per_m2,reference_jse_a_per_m2,"
              "actual_jb_a_per_m2,reference_jb_a_per_m2,"
              "actual_ji_a_per_m2,reference_ji_a_per_m2,"
              "actual_jsi_a_per_m2,reference_jsi_a_per_m2,"
              "actual_jph_a_per_m2,reference_jph_a_per_m2,"
              "actual_jcond_a_per_m2,reference_jcond_a_per_m2\n";

    const auto write_rows = [&](const char* role,
                                const std::vector<LegacyBenchmarkCurveSample>& actual,
                                const std::vector<LegacyBenchmarkCurveSample>& reference) {
        const std::size_t count = std::max(actual.size(), reference.size());
        for (std::size_t i = 0; i < count; ++i)
        {
            const LegacyBenchmarkCurveSample empty{};
            const auto& a = i < actual.size() ? actual[i] : empty;
            const auto& r = i < reference.size() ? reference[i] : empty;
            output << role << ',' << i << ',' << a.cycle_index << ',' << r.cycle_index << ','
                   << a.time_s << ',' << r.time_s << ',' << a.potential_v << ',' << r.potential_v
                   << ',' << (a.potential_v - r.potential_v) << ',' << a.jnet_a_per_m2 << ','
                   << r.jnet_a_per_m2 << ',' << a.je_a_per_m2 << ',' << r.je_a_per_m2 << ','
                   << a.jse_a_per_m2 << ',' << r.jse_a_per_m2 << ',' << a.jb_a_per_m2 << ','
                   << r.jb_a_per_m2 << ',' << a.ji_a_per_m2 << ',' << r.ji_a_per_m2 << ','
                   << a.jsi_a_per_m2 << ',' << r.jsi_a_per_m2 << ',' << a.jph_a_per_m2 << ','
                   << r.jph_a_per_m2 << ',' << a.jcond_a_per_m2 << ',' << r.jcond_a_per_m2
                   << '\n';
        }
    };

    write_rows("patch", actual_patch, reference_patch);
    write_rows("body", actual_body, reference_body);
    return true;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
