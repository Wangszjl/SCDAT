#include "DataAnalyzer.h"
#include "DensePlasmaSurfaceCharging.h"
#include "CouplingBase.h"
#include "InternalChargingCases.h"
#include "MultiPhysicsManager.h"
#include "PICFluidIntegration.h"
#include "PlasmaAnalysisCases.h"
#include "RadiationCases.h"
#include "MaterialDatabase.h"
#include "ModelRegistry.h"
#include "SpacecraftInternalChargingAlgorithm.h"
#include "SurfaceChargingCases.h"
#include "SurfaceDischargeArcAlgorithm.h"
#include "VacuumArcCases.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

namespace fs = std::filesystem;
namespace Particle = SCDAT::Particle;

using SCDAT::Output::DataAnalyzer;
using SCDAT::Output::DataSetSummary;
using SCDAT::Coupling::CouplingType;
using SCDAT::Coupling::FunctionalCoupling;
using SCDAT::Coupling::MultiPhysicsManager;
using SCDAT::Toolkit::InternalCharging::InternalChargingRadiationDrive;
using SCDAT::Toolkit::InternalCharging::InternalChargingScenarioPreset;
using SCDAT::Toolkit::InternalCharging::InternalChargingSourceMode;
using SCDAT::Toolkit::InternalCharging::SpacecraftInternalChargingAlgorithm;
using SCDAT::Toolkit::PlasmaAnalysis::PlasmaScenarioPreset;
using SCDAT::Toolkit::PlasmaAnalysis::PICFluidIntegration;
using SCDAT::Toolkit::Radiation::RadiationDoseAlgorithm;
using SCDAT::Toolkit::Radiation::RadiationScenarioPreset;
using SCDAT::Toolkit::SurfaceCharging::DensePlasmaSurfaceCharging;
using SCDAT::Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset;
using SCDAT::Toolkit::VacuumArc::VacuumArcScenarioPreset;
using SCDAT::Toolkit::VacuumArc::SurfaceDischargeArcAlgorithm;

struct RunSelection
{
    std::string preset_name;
    fs::path output_path;
};

struct InternalRadiationRunSelection
{
    std::string internal_preset_name;
    std::string radiation_preset_name;
    fs::path output_path;
};

struct InternalRadiationCouplingHistoryRow
{
    std::size_t macro_step = 0;
    std::size_t substep = 0;
    int iteration = 0;
    double relative_error = 0.0;
    double absolute_error = 0.0;
    double deposited_energy_delta_j_per_m2 = 0.0;
    double incident_current_density_a_per_m2 = 0.0;
    double feedback_scale = 1.0;
    double internal_max_electric_field_v_per_m = 0.0;
    double internal_average_dose_gy = 0.0;
    bool converged = false;
};

bool looksLikeOutputPath(const std::string& argument)
{
    const auto extension = fs::path(argument).extension().string();
    return extension == ".csv" || extension == ".vtk" || extension == ".h5" || extension == ".png";
}

std::string readTextFile(const fs::path& path)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        throw std::runtime_error("Unable to open surface config json: " + path.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::size_t skipWhitespace(const std::string& text, std::size_t index)
{
    while (index < text.size() && std::isspace(static_cast<unsigned char>(text[index])))
    {
        ++index;
    }
    return index;
}

std::optional<std::size_t> findValueStart(const std::string& text, const std::string& key)
{
    const auto key_pos = text.find("\"" + key + "\"");
    if (key_pos == std::string::npos)
    {
        return std::nullopt;
    }
    const auto colon_pos = text.find(':', key_pos);
    if (colon_pos == std::string::npos)
    {
        return std::nullopt;
    }
    const auto value_pos = skipWhitespace(text, colon_pos + 1);
    if (value_pos >= text.size())
    {
        return std::nullopt;
    }
    return value_pos;
}

std::optional<std::size_t> findMatchingDelimiter(const std::string& text, std::size_t start_index,
                                                 char open_char, char close_char)
{
    if (start_index >= text.size() || text[start_index] != open_char)
    {
        return std::nullopt;
    }

    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t i = start_index; i < text.size(); ++i)
    {
        const char ch = text[i];
        if (in_string)
        {
            if (escaped)
            {
                escaped = false;
                continue;
            }
            if (ch == '\\')
            {
                escaped = true;
                continue;
            }
            if (ch == '"')
            {
                in_string = false;
            }
            continue;
        }

        if (ch == '"')
        {
            in_string = true;
            continue;
        }
        if (ch == open_char)
        {
            ++depth;
            continue;
        }
        if (ch == close_char)
        {
            --depth;
            if (depth == 0)
            {
                return i;
            }
        }
    }
    return std::nullopt;
}

std::optional<std::string> extractStringField(const std::string& text, const std::string& key)
{
    const auto value_pos = findValueStart(text, key);
    if (!value_pos || text[*value_pos] != '"')
    {
        return std::nullopt;
    }

    std::string value;
    bool escaped = false;
    for (std::size_t i = *value_pos + 1; i < text.size(); ++i)
    {
        const char ch = text[i];
        if (escaped)
        {
            switch (ch)
            {
            case '"':
                value.push_back('"');
                break;
            case '\\':
                value.push_back('\\');
                break;
            case '/':
                value.push_back('/');
                break;
            case 'b':
                value.push_back('\b');
                break;
            case 'f':
                value.push_back('\f');
                break;
            case 'n':
                value.push_back('\n');
                break;
            case 'r':
                value.push_back('\r');
                break;
            case 't':
                value.push_back('\t');
                break;
            default:
                value.push_back(ch);
                break;
            }
            escaped = false;
            continue;
        }
        if (ch == '\\')
        {
            escaped = true;
            continue;
        }
        if (ch == '"')
        {
            return value;
        }
        value.push_back(ch);
    }
    return std::nullopt;
}

std::optional<double> extractNumberField(const std::string& text, const std::string& key)
{
    const auto value_pos = findValueStart(text, key);
    if (!value_pos)
    {
        return std::nullopt;
    }
    const auto number_start = text.find_first_of("+-0123456789.", *value_pos);
    if (number_start == std::string::npos)
    {
        return std::nullopt;
    }
    const auto number_end = text.find_first_not_of("+-0123456789.eE", number_start);
    try
    {
        return std::stod(text.substr(number_start, number_end - number_start));
    }
    catch (const std::exception&)
    {
        return std::nullopt;
    }
}

std::optional<bool> extractBoolField(const std::string& text, const std::string& key)
{
    const auto value_pos = findValueStart(text, key);
    if (!value_pos)
    {
        return std::nullopt;
    }
    if (text.compare(*value_pos, 4, "true") == 0)
    {
        return true;
    }
    if (text.compare(*value_pos, 5, "false") == 0)
    {
        return false;
    }
    return std::nullopt;
}

std::optional<std::string> extractObjectField(const std::string& text, const std::string& key)
{
    const auto value_pos = findValueStart(text, key);
    if (!value_pos || text[*value_pos] != '{')
    {
        return std::nullopt;
    }
    const auto end_pos = findMatchingDelimiter(text, *value_pos, '{', '}');
    if (!end_pos)
    {
        return std::nullopt;
    }
    return text.substr(*value_pos, *end_pos - *value_pos + 1);
}

std::optional<std::string> extractArrayField(const std::string& text, const std::string& key)
{
    const auto value_pos = findValueStart(text, key);
    if (!value_pos || text[*value_pos] != '[')
    {
        return std::nullopt;
    }
    const auto end_pos = findMatchingDelimiter(text, *value_pos, '[', ']');
    if (!end_pos)
    {
        return std::nullopt;
    }
    return text.substr(*value_pos, *end_pos - *value_pos + 1);
}

std::vector<double> parseNumberArray(const std::string& array_text)
{
    std::vector<double> values;
    std::size_t index = 0;
    while (index < array_text.size())
    {
        const auto token_start = array_text.find_first_of("+-0123456789.", index);
        if (token_start == std::string::npos)
        {
            break;
        }
        const auto token_end = array_text.find_first_not_of("+-0123456789.eE", token_start);
        try
        {
            values.push_back(std::stod(array_text.substr(token_start, token_end - token_start)));
        }
        catch (const std::exception&)
        {
            // Keep parsing later numbers even if one token is malformed.
        }
        index = (token_end == std::string::npos) ? array_text.size() : token_end;
    }
    return values;
}

std::vector<int> parseIntegerArray(const std::string& array_text)
{
    std::vector<int> values;
    for (const double number : parseNumberArray(array_text))
    {
        values.push_back(static_cast<int>(std::llround(number)));
    }
    return values;
}

std::vector<std::string> parseObjectArray(const std::string& array_text)
{
    std::vector<std::string> objects;
    if (array_text.empty() || array_text.front() != '[')
    {
        return objects;
    }

    std::size_t index = 1;
    while (index < array_text.size())
    {
        index = skipWhitespace(array_text, index);
        if (index >= array_text.size() || array_text[index] == ']')
        {
            break;
        }
        if (array_text[index] == ',')
        {
            ++index;
            continue;
        }
        if (array_text[index] != '{')
        {
            ++index;
            continue;
        }

        const auto end_pos = findMatchingDelimiter(array_text, index, '{', '}');
        if (!end_pos)
        {
            break;
        }
        objects.push_back(array_text.substr(index, *end_pos - index + 1));
        index = *end_pos + 1;
    }
    return objects;
}

std::string lowerCase(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return text;
}

std::string normalizeToken(std::string text)
{
    text = lowerCase(std::move(text));
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

const SCDAT::Material::MaterialProperty*
resolveMaterialAliasOrName(const std::string& text)
{
    static const SCDAT::Material::UnifiedMaterialDatabase kMaterialDatabase;
    return kMaterialDatabase.findByAliasOrName(text);
}

bool hasExplicitMaterialDefinition(const std::string& object_text)
{
    return extractNumberField(object_text, "id").has_value() ||
           extractStringField(object_text, "type").has_value() ||
           extractNumberField(object_text, "permittivity").has_value() ||
           extractNumberField(object_text, "conductivity").has_value() ||
           extractNumberField(object_text, "permeability").has_value() ||
           extractNumberField(object_text, "work_function_ev").has_value() ||
           extractNumberField(object_text, "breakdown_field_v_per_m").has_value() ||
           extractNumberField(object_text, "secondary_electron_yield").has_value() ||
           extractNumberField(object_text, "mass_density_kg_per_m3").has_value() ||
           extractObjectField(object_text, "scalar_properties").has_value();
}

std::optional<SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel>
parseSecondaryElectronEmissionModel(const std::string& text)
{
    using SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel;
    using Registry = SCDAT::Basic::StringModelRegistry<SecondaryElectronEmissionModel>;
    static const Registry kRegistry{
        {"whipple", SecondaryElectronEmissionModel::Whipple},
        {"sims", SecondaryElectronEmissionModel::Sims},
        {"katz", SecondaryElectronEmissionModel::Katz},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::SurfaceCharging::ElectronCollectionModelKind>
parseElectronCollectionModelKind(const std::string& text)
{
    using SCDAT::Toolkit::SurfaceCharging::ElectronCollectionModelKind;
    using Registry = SCDAT::Basic::StringModelRegistry<ElectronCollectionModelKind>;
    static const Registry kRegistry{
        {"omllike", ElectronCollectionModelKind::OmlLike},
        {"oml", ElectronCollectionModelKind::OmlLike},
        {"legacy", ElectronCollectionModelKind::OmlLike},
        {"shiftedenergy", ElectronCollectionModelKind::ShiftedEnergy},
        {"shifted", ElectronCollectionModelKind::ShiftedEnergy},
        {"barrierlimited", ElectronCollectionModelKind::BarrierLimited},
        {"barrier", ElectronCollectionModelKind::BarrierLimited},
        {"limited", ElectronCollectionModelKind::BarrierLimited},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Mesh::MaterialType> parseMaterialType(const std::string& text)
{
    const auto token = normalizeToken(text);
    if (token == "vacuum" || token == "vac")
    {
        return SCDAT::Mesh::MaterialType::VACUUM;
    }
    if (token == "dielectric" || token == "insulator")
    {
        return SCDAT::Mesh::MaterialType::DIELECTRIC;
    }
    if (token == "conductor" || token == "metal")
    {
        return SCDAT::Mesh::MaterialType::CONDUCTOR;
    }
    return std::nullopt;
}

std::optional<SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute>
parseSurfaceRuntimeRoute(const std::string& text)
{
    using SCDAT::Toolkit::SurfaceCharging::SurfaceRuntimeRoute;
    using Registry = SCDAT::Basic::StringModelRegistry<SurfaceRuntimeRoute>;
    static const Registry kRegistry{
        {"scdatunified", SurfaceRuntimeRoute::SCDATUnified},
        {"unified", SurfaceRuntimeRoute::SCDATUnified},
        {"scdat", SurfaceRuntimeRoute::SCDATUnified},
        {"legacybenchmark", SurfaceRuntimeRoute::LegacyBenchmark},
        {"legacy", SurfaceRuntimeRoute::LegacyBenchmark},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode>
parseSurfaceCurrentAlgorithmMode(const std::string& text)
{
    using SCDAT::Toolkit::SurfaceCharging::SurfaceCurrentAlgorithmMode;
    using Registry = SCDAT::Basic::StringModelRegistry<SurfaceCurrentAlgorithmMode>;
    static const Registry kRegistry{
        {"unifiedspisaligned", SurfaceCurrentAlgorithmMode::UnifiedSpisAligned},
        {"unified", SurfaceCurrentAlgorithmMode::UnifiedSpisAligned},
        {"legacyrefcompatible", SurfaceCurrentAlgorithmMode::LegacyRefCompatible},
        {"legacy", SurfaceCurrentAlgorithmMode::LegacyRefCompatible},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode>
parseSurfaceBenchmarkMode(const std::string& text)
{
    using SCDAT::Toolkit::SurfaceCharging::SurfaceBenchmarkMode;
    using Registry = SCDAT::Basic::StringModelRegistry<SurfaceBenchmarkMode>;
    static const Registry kRegistry{
        {"unifiedspisaligned", SurfaceBenchmarkMode::UnifiedSpisAligned},
        {"unified", SurfaceBenchmarkMode::UnifiedSpisAligned},
        {"legacyrefcompatible", SurfaceBenchmarkMode::LegacyRefCompatible},
        {"legacy", SurfaceBenchmarkMode::LegacyRefCompatible},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::SurfaceCharging::VolumeLinearSolverPolicy>
parseVolumeLinearSolverPolicy(const std::string& text)
{
    using SCDAT::Toolkit::SurfaceCharging::VolumeLinearSolverPolicy;
    using Registry = SCDAT::Basic::StringModelRegistry<VolumeLinearSolverPolicy>;
    static const Registry kRegistry{
        {"auto", VolumeLinearSolverPolicy::Auto},
        {"denseonly", VolumeLinearSolverPolicy::DenseOnly},
        {"dense", VolumeLinearSolverPolicy::DenseOnly},
        {"iterativeonly", VolumeLinearSolverPolicy::IterativeOnly},
        {"iterative", VolumeLinearSolverPolicy::IterativeOnly},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::Radiation::ParticleSpecies>
parseRadiationParticleSpecies(const std::string& text)
{
    using SCDAT::Toolkit::Radiation::ParticleSpecies;
    using Registry = SCDAT::Basic::StringModelRegistry<ParticleSpecies>;
    static const Registry kRegistry{
        {"electron", ParticleSpecies::Electron},
        {"e", ParticleSpecies::Electron},
        {"proton", ParticleSpecies::Proton},
        {"p", ParticleSpecies::Proton},
        {"heavyion", ParticleSpecies::HeavyIon},
        {"ion", ParticleSpecies::HeavyIon},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::Radiation::RadiationPhysicsList>
parseRadiationPhysicsList(const std::string& text)
{
    using SCDAT::Toolkit::Radiation::RadiationPhysicsList;
    using Registry = SCDAT::Basic::StringModelRegistry<RadiationPhysicsList>;
    static const Registry kRegistry{
        {"geant4emstandard", RadiationPhysicsList::Geant4EmStandard},
        {"emstandard", RadiationPhysicsList::Geant4EmStandard},
        {"standard", RadiationPhysicsList::Geant4EmStandard},
        {"geant4emlivermore", RadiationPhysicsList::Geant4EmLivermore},
        {"livermore", RadiationPhysicsList::Geant4EmLivermore},
        {"geant4empenelope", RadiationPhysicsList::Geant4EmPenelope},
        {"penelope", RadiationPhysicsList::Geant4EmPenelope},
        {"geant4spaceshielding", RadiationPhysicsList::Geant4SpaceShielding},
        {"spaceshielding", RadiationPhysicsList::Geant4SpaceShielding},
        {"shielding", RadiationPhysicsList::Geant4SpaceShielding},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::PlasmaAnalysis::NonEquilibriumClosureModel>
parseNonEquilibriumClosureModel(const std::string& text)
{
    using SCDAT::Toolkit::PlasmaAnalysis::NonEquilibriumClosureModel;
    using Registry = SCDAT::Basic::StringModelRegistry<NonEquilibriumClosureModel>;
    static const Registry kRegistry{
        {"disabled", NonEquilibriumClosureModel::Disabled},
        {"off", NonEquilibriumClosureModel::Disabled},
        {"none", NonEquilibriumClosureModel::Disabled},
        {"twotemperaturerelaxation", NonEquilibriumClosureModel::TwoTemperatureRelaxation},
        {"twotemperature", NonEquilibriumClosureModel::TwoTemperatureRelaxation},
        {"relaxation", NonEquilibriumClosureModel::TwoTemperatureRelaxation},
        {"twotemp", NonEquilibriumClosureModel::TwoTemperatureRelaxation},
        {"collisionalthermalization", NonEquilibriumClosureModel::CollisionalThermalization},
        {"thermalization", NonEquilibriumClosureModel::CollisionalThermalization},
        {"collisional", NonEquilibriumClosureModel::CollisionalThermalization},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::PlasmaAnalysis::TurbulenceClosureModel>
parseTurbulenceClosureModel(const std::string& text)
{
    using SCDAT::Toolkit::PlasmaAnalysis::TurbulenceClosureModel;
    using Registry = SCDAT::Basic::StringModelRegistry<TurbulenceClosureModel>;
    static const Registry kRegistry{
        {"disabled", TurbulenceClosureModel::Disabled},
        {"off", TurbulenceClosureModel::Disabled},
        {"none", TurbulenceClosureModel::Disabled},
        {"mixinglengtheddydiffusivity", TurbulenceClosureModel::MixingLengthEddyDiffusivity},
        {"mixinglength", TurbulenceClosureModel::MixingLengthEddyDiffusivity},
        {"mixing", TurbulenceClosureModel::MixingLengthEddyDiffusivity},
        {"ml", TurbulenceClosureModel::MixingLengthEddyDiffusivity},
        {"collisionaldamping", TurbulenceClosureModel::CollisionalDamping},
        {"damping", TurbulenceClosureModel::CollisionalDamping},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::VacuumArc::ArcPicAlignmentMode>
parseArcPicAlignmentMode(const std::string& text)
{
    using SCDAT::Toolkit::VacuumArc::ArcPicAlignmentMode;
    using Registry = SCDAT::Basic::StringModelRegistry<ArcPicAlignmentMode>;
    static const Registry kRegistry{
        {"arcpicaligned", ArcPicAlignmentMode::ArcPicAligned},
        {"aligned", ArcPicAlignmentMode::ArcPicAligned},
        {"legacybaseline", ArcPicAlignmentMode::LegacyBaseline},
        {"legacy", ArcPicAlignmentMode::LegacyBaseline},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::VacuumArc::ArcIntegratorMode>
parseArcIntegratorMode(const std::string& text)
{
    using SCDAT::Toolkit::VacuumArc::ArcIntegratorMode;
    using Registry = SCDAT::Basic::StringModelRegistry<ArcIntegratorMode>;
    static const Registry kRegistry{
        {"explicitgrowth", ArcIntegratorMode::ExplicitGrowth},
        {"explicit", ArcIntegratorMode::ExplicitGrowth},
        {"boundedrelaxation", ArcIntegratorMode::BoundedRelaxation},
        {"bounded", ArcIntegratorMode::BoundedRelaxation},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::VacuumArc::ArcEmissionStrategyType>
parseArcEmissionStrategyType(const std::string& text)
{
    using SCDAT::Toolkit::VacuumArc::ArcEmissionStrategyType;
    using Registry = SCDAT::Basic::StringModelRegistry<ArcEmissionStrategyType>;
    static const Registry kRegistry{
        {"linearthreshold", ArcEmissionStrategyType::LinearThreshold},
        {"linear", ArcEmissionStrategyType::LinearThreshold},
        {"fowlernordheimlike", ArcEmissionStrategyType::FowlerNordheimLike},
        {"fowlernordheim", ArcEmissionStrategyType::FowlerNordheimLike},
        {"fn", ArcEmissionStrategyType::FowlerNordheimLike},
    };
    return kRegistry.tryParse(text);
}

std::optional<SCDAT::Toolkit::VacuumArc::SurfaceCircuitLoadModel>
parseArcSurfaceLoadModel(const std::string& text)
{
    using SCDAT::Toolkit::VacuumArc::SurfaceCircuitLoadModel;
    using Registry = SCDAT::Basic::StringModelRegistry<SurfaceCircuitLoadModel>;
    static const Registry kRegistry{
        {"legacyleakagecapacitor", SurfaceCircuitLoadModel::LegacyLeakageCapacitor},
        {"legacy", SurfaceCircuitLoadModel::LegacyLeakageCapacitor},
        {"rc", SurfaceCircuitLoadModel::LegacyLeakageCapacitor},
        {"resistiveshunt", SurfaceCircuitLoadModel::ResistiveShunt},
        {"resistive", SurfaceCircuitLoadModel::ResistiveShunt},
        {"resistor", SurfaceCircuitLoadModel::ResistiveShunt},
        {"rlbranch", SurfaceCircuitLoadModel::RlBranch},
        {"rl", SurfaceCircuitLoadModel::RlBranch},
    };
    return kRegistry.tryParse(text);
}

void applyPlasmaParametersFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters& plasma)
{
    if (const auto value = extractNumberField(object_text, "electron_density_m3"); value)
    {
        plasma.electron_density_m3 = *value;
    }
    if (const auto value = extractNumberField(object_text, "ion_density_m3"); value)
    {
        plasma.ion_density_m3 = *value;
    }
    if (const auto value = extractNumberField(object_text, "electron_temperature_ev"); value)
    {
        plasma.electron_temperature_ev = *value;
    }
    if (const auto value = extractNumberField(object_text, "ion_temperature_ev"); value)
    {
        plasma.ion_temperature_ev = *value;
    }
    if (const auto value = extractNumberField(object_text, "ion_mass_amu"); value)
    {
        plasma.ion_mass_amu = *value;
    }
    if (const auto value = extractNumberField(object_text, "neutral_density_m3"); value)
    {
        plasma.neutral_density_m3 = *value;
    }
    if (const auto value = extractNumberField(object_text, "electric_field_v_per_m"); value)
    {
        plasma.electric_field_v_per_m = *value;
    }
    if (const auto value = extractNumberField(object_text, "pressure_pa"); value)
    {
        plasma.pressure_pa = *value;
    }
}

void applyReactionCollisionConfigFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::PlasmaAnalysis::ReactionCollisionConfig& config)
{
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(object_text, key); value)
        {
            target = *value;
        }
    };
    auto apply_bool = [&](const std::string& key, bool& target) {
        if (const auto value = extractBoolField(object_text, key); value)
        {
            target = *value;
        }
    };

    apply_bool("enable_electron_impact_ionization", config.enable_electron_impact_ionization);
    apply_bool("enable_radiative_recombination", config.enable_radiative_recombination);
    apply_bool("enable_charge_exchange", config.enable_charge_exchange);
    apply_bool("enable_elastic_momentum_transfer", config.enable_elastic_momentum_transfer);

    apply_number("ionization_threshold_ev", config.ionization_threshold_ev);
    apply_number("ionization_cross_section_m2", config.ionization_cross_section_m2);
    apply_number("charge_exchange_cross_section_m2", config.charge_exchange_cross_section_m2);
    apply_number("elastic_momentum_cross_section_m2", config.elastic_momentum_cross_section_m2);
    apply_number("max_relative_density_change_per_step",
                 config.max_relative_density_change_per_step);
}

void applyAdvancedClosureConfigFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::PlasmaAnalysis::AdvancedClosureConfig& config)
{
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(object_text, key); value)
        {
            target = *value;
        }
    };
    auto apply_bool = [&](const std::string& key, bool& target) {
        if (const auto value = extractBoolField(object_text, key); value)
        {
            target = *value;
        }
    };

    apply_bool("enable_non_equilibrium_closure", config.enable_non_equilibrium_closure);
    apply_bool("enable_turbulence_closure", config.enable_turbulence_closure);
    apply_number("non_equilibrium_relaxation_gain", config.non_equilibrium_relaxation_gain);
    apply_number("non_equilibrium_ratio_cap", config.non_equilibrium_ratio_cap);
    apply_number("turbulence_mixing_length_m", config.turbulence_mixing_length_m);
    apply_number("turbulence_gain", config.turbulence_gain);
    apply_number("max_temperature_correction_ev_per_step",
                 config.max_temperature_correction_ev_per_step);
    apply_number("max_relative_density_correction_per_step",
                 config.max_relative_density_correction_per_step);

    if (const auto model = extractStringField(object_text, "non_equilibrium_model"); model)
    {
        const auto parsed = parseNonEquilibriumClosureModel(*model);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported non_equilibrium_model in plasma config json: " + *model);
        }
        config.non_equilibrium_model = *parsed;
    }

    if (const auto model = extractStringField(object_text, "turbulence_model"); model)
    {
        const auto parsed = parseTurbulenceClosureModel(*model);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported turbulence_model in plasma config json: " +
                                     *model);
        }
        config.turbulence_model = *parsed;
    }
}

void applyVector3DFromJsonObject(const std::string& parent_text, const std::string& key,
                                 SCDAT::Geometry::Vector3D& target)
{
    if (const auto object_text = extractObjectField(parent_text, key); object_text)
    {
        const double x = extractNumberField(*object_text, "x").value_or(target.x());
        const double y = extractNumberField(*object_text, "y").value_or(target.y());
        const double z = extractNumberField(*object_text, "z").value_or(target.z());
        target = SCDAT::Geometry::Vector3D(x, y, z);
        return;
    }

    if (const auto array_text = extractArrayField(parent_text, key); array_text)
    {
        const auto values = parseNumberArray(*array_text);
        if (values.size() >= 3)
        {
            target = SCDAT::Geometry::Vector3D(values[0], values[1], values[2]);
        }
    }
}

void applyEmissionParametersFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::SurfaceCharging::EmissionModelParameters& emission)
{
    if (const auto value = extractNumberField(object_text, "surface_temperature_k"); value)
    {
        emission.surface_temperature_k = *value;
    }
    if (const auto value = extractNumberField(object_text, "enhancement_factor"); value)
    {
        emission.enhancement_factor = *value;
    }
    if (const auto value = extractNumberField(object_text, "photon_flux_m2_s"); value)
    {
        emission.photon_flux_m2_s = *value;
    }
}

void applyMaterialScalarPropertiesFromJsonObject(
    const std::string& object_text,
    SCDAT::Material::MaterialProperty& material)
{
    std::size_t index = 0;
    while (index < object_text.size())
    {
        const auto key_start = object_text.find('"', index);
        if (key_start == std::string::npos)
        {
            break;
        }
        const auto key_end = object_text.find('"', key_start + 1);
        if (key_end == std::string::npos)
        {
            break;
        }
        const std::string key = object_text.substr(key_start + 1, key_end - key_start - 1);
        const auto colon_pos = object_text.find(':', key_end + 1);
        if (colon_pos == std::string::npos)
        {
            break;
        }
        const auto number_start = object_text.find_first_of("+-0123456789.", colon_pos + 1);
        if (number_start == std::string::npos)
        {
            index = colon_pos + 1;
            continue;
        }
        const auto number_end = object_text.find_first_not_of("+-0123456789.eE", number_start);
        try
        {
            const double value =
                std::stod(object_text.substr(number_start, number_end - number_start));
            material.setScalarProperty(key, value);
        }
        catch (const std::exception&)
        {
            // Ignore malformed scalar value and continue parsing remaining keys.
        }
        index = (number_end == std::string::npos) ? object_text.size() : number_end;
    }
}

void applyMaterialFromJsonObject(const std::string& object_text,
                                 SCDAT::Material::MaterialProperty& material)
{
    bool resolved_material_reference = false;
    if (const auto alias = extractStringField(object_text, "alias"); alias)
    {
        const auto* resolved = resolveMaterialAliasOrName(*alias);
        if (resolved == nullptr)
        {
            throw std::runtime_error("Unknown material alias in surface config json: " +
                                     *alias);
        }
        material = *resolved;
        resolved_material_reference = true;
    }

    std::optional<std::string> requested_name =
        extractStringField(object_text, "name");
    bool requested_name_is_alias = false;
    if (requested_name)
    {
        if (const auto* resolved = resolveMaterialAliasOrName(*requested_name); resolved)
        {
            material = *resolved;
            requested_name_is_alias =
                normalizeToken(*requested_name) != normalizeToken(resolved->getName());
            resolved_material_reference = true;
        }
        else if (!hasExplicitMaterialDefinition(object_text))
        {
            throw std::runtime_error(
                "Unknown material name without explicit material fields in surface config json: " +
                *requested_name);
        }
    }

    if (const auto id_value = extractNumberField(object_text, "id"); id_value)
    {
        material.setId(static_cast<SCDAT::Mesh::MaterialId>(
            std::max(0.0, std::floor(*id_value + 0.5))));
    }
    if (const auto type_text = extractStringField(object_text, "type"); type_text)
    {
        if (const auto type = parseMaterialType(*type_text); type)
        {
            material.setType(*type);
        }
    }
    if (requested_name && !requested_name_is_alias)
    {
        material.setName(*requested_name);
    }
    else if (resolved_material_reference)
    {
        // Keep canonical material name after alias resolution.
        material.setName(material.getName());
    }
    if (const auto value = extractNumberField(object_text, "permittivity"); value)
    {
        material.setPermittivity(*value);
    }
    if (const auto value = extractNumberField(object_text, "conductivity"); value)
    {
        material.setConductivity(*value);
    }
    if (const auto value = extractNumberField(object_text, "permeability"); value)
    {
        material.setPermeability(*value);
    }
    if (const auto value = extractNumberField(object_text, "work_function_ev"); value)
    {
        material.setWorkFunctionEv(*value);
    }
    if (const auto value = extractNumberField(object_text, "breakdown_field_v_per_m"); value)
    {
        material.setBreakdownFieldVPerM(*value);
    }
    if (const auto value = extractNumberField(object_text, "secondary_electron_yield"); value)
    {
        material.setSecondaryElectronYield(*value);
    }
    if (const auto value = extractNumberField(object_text, "mass_density_kg_per_m3"); value)
    {
        material.setMassDensityKgPerM3(*value);
    }

    if (const auto scalar_properties = extractObjectField(object_text, "scalar_properties");
        scalar_properties)
    {
        applyMaterialScalarPropertiesFromJsonObject(*scalar_properties, material);
    }
}

Particle::ResolvedSpectrum buildSingleMaxwellSpectrum(const Particle::ParticleType& particle_type,
                                                      double density_m3, double temperature_ev,
                                                      double drift_speed_m_per_s, double mass_amu)
{
    Particle::SamplingParameters params;
    params.density = density_m3;
    params.bulk_speed = drift_speed_m_per_s;
    params.thermal_speed = std::sqrt(2.0 * std::max(1.0e-6, temperature_ev) *
                                     1.602176634e-19 /
                                     (particle_type == Particle::ParticleType::ELECTRON
                                          ? 9.1093837015e-31
                                          : 1.67262192369e-27 * std::max(1.0, mass_amu)));
    auto spectrum = Particle::ParticleSource::buildResolvedSpectrum(
        particle_type, Particle::SpatialSamplingModel::SINGLE_MAXWELL, params,
        Particle::SpectrumUsage::CurrentBalance);
    for (auto& population : spectrum.populations)
    {
        population.mass_amu = mass_amu;
    }
    return spectrum;
}

Particle::ResolvedSpectrum buildDoubleMaxwellSpectrum(const Particle::ParticleType& particle_type,
                                                      double total_density_m3,
                                                      double cold_temperature_ev,
                                                      double hot_temperature_ev,
                                                      double hot_fraction,
                                                      double drift_speed_m_per_s,
                                                      double mass_amu)
{
    Particle::SamplingParameters params;
    params.density = total_density_m3;
    params.bulk_speed = drift_speed_m_per_s;
    params.hot_fraction = hot_fraction;
    params.thermal_speed = std::sqrt(2.0 * std::max(1.0e-6, cold_temperature_ev) *
                                     1.602176634e-19 /
                                     (particle_type == Particle::ParticleType::ELECTRON
                                          ? 9.1093837015e-31
                                          : 1.67262192369e-27 * std::max(1.0, mass_amu)));
    params.hot_thermal_speed = std::sqrt(2.0 * std::max(1.0e-6, hot_temperature_ev) *
                                         1.602176634e-19 /
                                         (particle_type == Particle::ParticleType::ELECTRON
                                              ? 9.1093837015e-31
                                              : 1.67262192369e-27 * std::max(1.0, mass_amu)));
    auto spectrum = Particle::ParticleSource::buildResolvedSpectrum(
        particle_type, Particle::SpatialSamplingModel::DOUBLE_MAXWELL, params,
        Particle::SpectrumUsage::CurrentBalance);
    for (auto& population : spectrum.populations)
    {
        population.mass_amu = mass_amu;
    }
    return spectrum;
}

bool applySpectrumFromJsonObject(const std::string& spectrum_text,
                                 const Particle::ParticleType& particle_type,
                                 Particle::ResolvedSpectrum& destination)
{
    const std::string mode = lowerCase(extractStringField(spectrum_text, "mode").value_or(""));
    const double default_mass_amu =
        (particle_type == Particle::ParticleType::ELECTRON) ? 5.48579909070e-4 : 1.0;
    const double mass_amu = extractNumberField(spectrum_text, "mass_amu").value_or(default_mass_amu);

    if (mode == "single_maxwell")
    {
        const double density_m3 = extractNumberField(spectrum_text, "density_m3").value_or(1.0e12);
        const double temperature_ev =
            extractNumberField(spectrum_text, "temperature_ev").value_or(1.0);
        const double drift_speed_m_per_s =
            extractNumberField(spectrum_text, "drift_speed_m_per_s").value_or(0.0);
        destination = buildSingleMaxwellSpectrum(
            particle_type, density_m3, temperature_ev, drift_speed_m_per_s, mass_amu);
        return true;
    }

    if (mode == "double_maxwell")
    {
        const double density_m3 = extractNumberField(spectrum_text, "density_m3").value_or(1.0e12);
        const double cold_temperature_ev =
            extractNumberField(spectrum_text, "cold_temperature_ev").value_or(1.0);
        const double hot_temperature_ev =
            extractNumberField(spectrum_text, "hot_temperature_ev").value_or(10.0);
        const double hot_fraction =
            std::clamp(extractNumberField(spectrum_text, "hot_fraction").value_or(0.1), 0.0, 1.0);
        const double drift_speed_m_per_s =
            extractNumberField(spectrum_text, "drift_speed_m_per_s").value_or(0.0);
        destination = buildDoubleMaxwellSpectrum(particle_type, density_m3, cold_temperature_ev,
                                                 hot_temperature_ev, hot_fraction,
                                                 drift_speed_m_per_s, mass_amu);
        return true;
    }

    if (mode == "tabulated" || mode == "table")
    {
        const auto energy_array_text = extractArrayField(spectrum_text, "energy_grid_ev");
        const auto flux_array_text = extractArrayField(spectrum_text, "differential_number_flux");
        if (!energy_array_text || !flux_array_text)
        {
            return false;
        }

        const auto energy = parseNumberArray(*energy_array_text);
        const auto flux = parseNumberArray(*flux_array_text);
        if (energy.size() < 2 || energy.size() != flux.size())
        {
            return false;
        }

        Particle::ResolvedSpectrum spectrum;
        spectrum.model = Particle::SpatialSamplingModel::TABULATED;
        spectrum.energy_grid_ev = energy;
        spectrum.differential_number_flux = flux;

        Particle::SpectrumPopulation primary_population;
        primary_population.mass_amu = mass_amu;
        primary_population.density_m3 = extractNumberField(spectrum_text, "density_m3").value_or(0.0);
        primary_population.temperature_ev =
            extractNumberField(spectrum_text, "temperature_ev").value_or(0.0);
        primary_population.drift_speed_m_per_s =
            extractNumberField(spectrum_text, "drift_speed_m_per_s").value_or(0.0);
        primary_population.weight = 1.0;
        spectrum.populations.push_back(primary_population);

        destination = std::move(spectrum);
        return true;
    }

    return false;
}

void applyStructuredTopologyFromJson(
    const std::string& config_scope,
    SCDAT::Toolkit::SurfaceCharging::SurfaceChargingConfig& config)
{
    using SCDAT::Toolkit::SurfaceCharging::BodyBoundaryGroup;
    using SCDAT::Toolkit::SurfaceCharging::ElectronCollectionModelKind;
    using SCDAT::Toolkit::SurfaceCharging::EmissionModelParameters;
    using SCDAT::Toolkit::SurfaceCharging::PatchBoundaryGroup;
    using SCDAT::Toolkit::SurfaceCharging::PatchInterfaceConfig;
    using SCDAT::Toolkit::SurfaceCharging::SecondaryElectronEmissionModel;
    using SCDAT::Toolkit::SurfaceCharging::StructureBodyConfig;
    using SCDAT::Toolkit::SurfaceCharging::SurfaceBoundaryMapping;
    using SCDAT::Toolkit::SurfaceCharging::SurfacePatchConfig;
    using SCDAT::Toolkit::PlasmaAnalysis::PlasmaParameters;

    if (const auto array_text = extractArrayField(config_scope, "bodies"); array_text)
    {
        std::vector<StructureBodyConfig> bodies;
        for (const auto& object_text : parseObjectArray(*array_text))
        {
            StructureBodyConfig body;
            body.id = extractStringField(object_text, "id").value_or("");
            if (body.id.empty())
            {
                continue;
            }
            body.area_m2 = extractNumberField(object_text, "area_m2").value_or(body.area_m2);
            body.floating = extractBoolField(object_text, "floating").value_or(body.floating);
            body.initial_potential_v =
                extractNumberField(object_text, "initial_potential_v").value_or(
                    body.initial_potential_v);
            body.capacitance_f =
                extractNumberField(object_text, "capacitance_f").value_or(body.capacitance_f);
            body.leakage_conductance_s =
                extractNumberField(object_text, "leakage_conductance_s")
                    .value_or(body.leakage_conductance_s);
            body.fixed_potential =
                extractBoolField(object_text, "fixed_potential").value_or(body.fixed_potential);
            body.fixed_value_v =
                extractNumberField(object_text, "fixed_value_v").value_or(body.fixed_value_v);
            bodies.push_back(std::move(body));
        }
        config.bodies = std::move(bodies);
    }

    if (const auto array_text = extractArrayField(config_scope, "patches"); array_text)
    {
        std::vector<SurfacePatchConfig> patches;
        for (const auto& object_text : parseObjectArray(*array_text))
        {
            SurfacePatchConfig patch;
            patch.id = extractStringField(object_text, "id").value_or("");
            patch.body_id = extractStringField(object_text, "body_id").value_or("");
            if (patch.id.empty() || patch.body_id.empty())
            {
                continue;
            }

            patch.area_m2 = extractNumberField(object_text, "area_m2").value_or(patch.area_m2);
            patch.initial_potential_v =
                extractNumberField(object_text, "initial_potential_v")
                    .value_or(patch.initial_potential_v);
            patch.capacitance_f =
                extractNumberField(object_text, "capacitance_f").value_or(patch.capacitance_f);

            if (const auto value = extractNumberField(object_text, "patch_incidence_angle_deg"); value)
            {
                patch.patch_incidence_angle_deg = *value;
            }
            if (const auto value = extractNumberField(object_text, "patch_flow_angle_deg"); value)
            {
                patch.patch_flow_angle_deg = *value;
            }
            if (const auto value = extractNumberField(object_text, "photoelectron_temperature_ev"); value)
            {
                patch.photoelectron_temperature_ev = *value;
            }
            if (const auto value = extractNumberField(object_text, "patch_photo_current_density_a_per_m2");
                value)
            {
                patch.patch_photo_current_density_a_per_m2 = *value;
            }
            if (const auto value = extractNumberField(object_text, "electron_collection_coefficient"); value)
            {
                patch.electron_collection_coefficient = *value;
            }
            if (const auto value = extractNumberField(object_text, "ion_collection_coefficient"); value)
            {
                patch.ion_collection_coefficient = *value;
            }
            if (const auto value = extractBoolField(object_text, "has_electron_spectrum"); value)
            {
                patch.has_electron_spectrum = *value;
            }
            if (const auto value = extractBoolField(object_text, "has_ion_spectrum"); value)
            {
                patch.has_ion_spectrum = *value;
            }

            if (const auto see_model = extractStringField(object_text, "reference_see_model");
                see_model)
            {
                const auto parsed_see_model = parseSecondaryElectronEmissionModel(*see_model);
                if (!parsed_see_model)
                {
                    throw std::runtime_error("Unsupported patch reference_see_model: " +
                                             *see_model);
                }
                patch.reference_see_model = *parsed_see_model;
            }
            if (const auto model_text = extractStringField(object_text, "electron_collection_model");
                model_text)
            {
                const auto parsed_model = parseElectronCollectionModelKind(*model_text);
                if (!parsed_model)
                {
                    throw std::runtime_error("Unsupported patch electron_collection_model: " +
                                             *model_text);
                }
                patch.electron_collection_model = *parsed_model;
            }

            if (const auto material_text = extractObjectField(object_text, "material");
                material_text)
            {
                SCDAT::Material::MaterialProperty material = config.material;
                applyMaterialFromJsonObject(*material_text, material);
                patch.material = std::move(material);
            }

            if (const auto plasma_text = extractObjectField(object_text, "plasma"); plasma_text)
            {
                PlasmaParameters plasma = config.plasma;
                applyPlasmaParametersFromJsonObject(*plasma_text, plasma);
                patch.plasma = std::move(plasma);
            }

            if (const auto emission_text = extractObjectField(object_text, "emission");
                emission_text)
            {
                EmissionModelParameters emission = config.emission;
                applyEmissionParametersFromJsonObject(*emission_text, emission);
                patch.emission = std::move(emission);
            }

            if (const auto spectrum_text = extractObjectField(object_text, "electron_spectrum");
                spectrum_text)
            {
                Particle::ResolvedSpectrum spectrum;
                if (!applySpectrumFromJsonObject(*spectrum_text, Particle::ParticleType::ELECTRON,
                                                 spectrum))
                {
                    throw std::runtime_error(
                        "Failed to parse patch electron_spectrum in surface config json.");
                }
                patch.electron_spectrum = std::move(spectrum);
                patch.has_electron_spectrum = true;
            }
            if (const auto spectrum_text = extractObjectField(object_text, "ion_spectrum");
                spectrum_text)
            {
                Particle::ResolvedSpectrum spectrum;
                if (!applySpectrumFromJsonObject(*spectrum_text, Particle::ParticleType::ION,
                                                 spectrum))
                {
                    throw std::runtime_error(
                        "Failed to parse patch ion_spectrum in surface config json.");
                }
                patch.ion_spectrum = std::move(spectrum);
                patch.has_ion_spectrum = true;
            }

            patches.push_back(std::move(patch));
        }
        config.patches = std::move(patches);
    }

    if (const auto array_text = extractArrayField(config_scope, "interfaces"); array_text)
    {
        std::vector<PatchInterfaceConfig> interfaces;
        for (const auto& object_text : parseObjectArray(*array_text))
        {
            PatchInterfaceConfig interface_config;
            interface_config.id = extractStringField(object_text, "id").value_or("");
            interface_config.from_id = extractStringField(object_text, "from_id").value_or("");
            interface_config.to_id = extractStringField(object_text, "to_id").value_or("");
            if (interface_config.from_id.empty() || interface_config.to_id.empty())
            {
                continue;
            }
            interface_config.conductance_s =
                extractNumberField(object_text, "conductance_s").value_or(
                    interface_config.conductance_s);
            interface_config.resistance_ohm =
                extractNumberField(object_text, "resistance_ohm").value_or(
                    interface_config.resistance_ohm);
            interface_config.bias_v =
                extractNumberField(object_text, "bias_v").value_or(interface_config.bias_v);
            interface_config.dynamic_conductance =
                extractBoolField(object_text, "dynamic_conductance")
                    .value_or(interface_config.dynamic_conductance);
            interfaces.push_back(std::move(interface_config));
        }
        config.interfaces = std::move(interfaces);
    }

    if (const auto array_text = extractArrayField(config_scope, "body_boundary_groups"); array_text)
    {
        std::vector<BodyBoundaryGroup> groups;
        for (const auto& object_text : parseObjectArray(*array_text))
        {
            BodyBoundaryGroup group;
            group.id = extractStringField(object_text, "id").value_or("");
            group.body_id = extractStringField(object_text, "body_id").value_or("");
            if (group.id.empty() || group.body_id.empty())
            {
                continue;
            }
            group.external_group_name =
                extractStringField(object_text, "external_group_name").value_or("");
            if (const auto faces_text = extractArrayField(object_text, "boundary_face_ids"); faces_text)
            {
                group.boundary_face_ids = parseIntegerArray(*faces_text);
            }
            groups.push_back(std::move(group));
        }
        config.body_boundary_groups = std::move(groups);
    }

    if (const auto array_text = extractArrayField(config_scope, "patch_boundary_groups"); array_text)
    {
        std::vector<PatchBoundaryGroup> groups;
        for (const auto& object_text : parseObjectArray(*array_text))
        {
            PatchBoundaryGroup group;
            group.id = extractStringField(object_text, "id").value_or("");
            group.patch_id = extractStringField(object_text, "patch_id").value_or("");
            if (group.id.empty() || group.patch_id.empty())
            {
                continue;
            }
            group.external_group_name =
                extractStringField(object_text, "external_group_name").value_or("");
            if (const auto faces_text = extractArrayField(object_text, "boundary_face_ids"); faces_text)
            {
                group.boundary_face_ids = parseIntegerArray(*faces_text);
            }
            groups.push_back(std::move(group));
        }
        config.patch_boundary_groups = std::move(groups);
    }

    if (const auto array_text = extractArrayField(config_scope, "boundary_mappings"); array_text)
    {
        std::vector<SurfaceBoundaryMapping> mappings;
        for (const auto& object_text : parseObjectArray(*array_text))
        {
            SurfaceBoundaryMapping mapping;
            mapping.node_id = extractStringField(object_text, "node_id").value_or("");
            mapping.boundary_group_id =
                extractStringField(object_text, "boundary_group_id").value_or("");
            if (mapping.node_id.empty() || mapping.boundary_group_id.empty())
            {
                continue;
            }
            mapping.required = extractBoolField(object_text, "required").value_or(true);
            mappings.push_back(std::move(mapping));
        }
        config.boundary_mappings = std::move(mappings);
    }
}

PlasmaScenarioPreset loadPlasmaScenarioPresetFromJson(const fs::path& json_path,
                                                      fs::path& output_path)
{
    const std::string content = readTextFile(json_path);
    const auto run_scope = extractObjectField(content, "run").value_or(content);
    const auto config_scope = extractObjectField(content, "config").value_or(content);

    auto preset = SCDAT::Toolkit::PlasmaAnalysis::makeDefaultPlasmaScenarioPreset();
    const auto base_preset = extractStringField(content, "base_preset")
                                 .value_or(extractStringField(content, "preset").value_or(""));
    if (!base_preset.empty())
    {
        if (!SCDAT::Toolkit::PlasmaAnalysis::tryGetPlasmaScenarioPreset(base_preset, preset))
        {
            throw std::runtime_error("Unknown plasma base preset in json: " + base_preset);
        }
    }

    if (const auto name = extractStringField(content, "name"); name)
    {
        preset.name = *name;
    }

    if (const auto time_step_s = extractNumberField(run_scope, "time_step_s"); time_step_s)
    {
        preset.config.time_step_s = *time_step_s;
    }
    if (const auto steps = extractNumberField(run_scope, "steps"); steps)
    {
        preset.steps = static_cast<std::size_t>(std::max(0.0, std::floor(*steps + 0.5)));
    }

    if (const auto output_csv = extractStringField(content, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else if (const auto output_csv = extractStringField(run_scope, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else
    {
        output_path = preset.default_output_csv;
    }

    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = *value;
        }
    };

    applyVector3DFromJsonObject(config_scope, "domain_size", preset.config.domain_size);
    applyVector3DFromJsonObject(config_scope, "domain_size_m", preset.config.domain_size);
    applyVector3DFromJsonObject(config_scope, "resolution", preset.config.resolution);

    apply_number("time_step_s", preset.config.time_step_s);
    apply_number("dense_plasma_threshold_m3", preset.config.dense_plasma_threshold_m3);
    apply_number("initial_potential_v", preset.config.initial_potential_v);

    auto plasma_text = extractObjectField(config_scope, "initial_plasma");
    if (!plasma_text)
    {
        plasma_text = extractObjectField(config_scope, "plasma_model");
    }
    if (!plasma_text)
    {
        plasma_text = extractObjectField(config_scope, "plasma");
    }
    if (plasma_text)
    {
        applyPlasmaParametersFromJsonObject(*plasma_text, preset.config.initial_plasma);
    }

    if (const auto reaction_text = extractObjectField(config_scope, "reaction_collision");
        reaction_text)
    {
        applyReactionCollisionConfigFromJsonObject(*reaction_text,
                                                   preset.config.reaction_collision);
    }

    if (const auto closure_text = extractObjectField(config_scope, "advanced_closure");
        closure_text)
    {
        applyAdvancedClosureConfigFromJsonObject(*closure_text,
                                                 preset.config.advanced_closure);
    }

    return preset;
}

SurfaceChargingScenarioPreset loadSurfaceScenarioPresetFromJson(const fs::path& json_path,
                                                                fs::path& output_path)
{
    const std::string content = readTextFile(json_path);
    const auto run_scope = extractObjectField(content, "run").value_or(content);
    const auto config_scope = extractObjectField(content, "config").value_or(content);

    auto preset = SCDAT::Toolkit::SurfaceCharging::makeDefaultSurfaceChargingScenarioPreset();
    const auto base_preset = extractStringField(content, "base_preset")
                                 .value_or(extractStringField(content, "preset").value_or(""));
    if (!base_preset.empty())
    {
        if (!SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(base_preset, preset))
        {
            throw std::runtime_error("Unknown surface base preset in json: " + base_preset);
        }
    }

    if (const auto name = extractStringField(content, "name"); name)
    {
        preset.name = *name;
    }

    if (const auto time_step_s = extractNumberField(run_scope, "time_step_s"); time_step_s)
    {
        preset.time_step_s = *time_step_s;
    }
    if (const auto steps = extractNumberField(run_scope, "steps"); steps)
    {
        preset.steps = static_cast<std::size_t>(std::max(0.0, std::floor(*steps + 0.5)));
    }
    if (const auto adaptive = extractBoolField(run_scope, "adaptive_time_stepping"); adaptive)
    {
        preset.adaptive_time_stepping = *adaptive;
    }
    if (const auto total_duration = extractNumberField(run_scope, "total_duration_s"); total_duration)
    {
        preset.total_duration_s = *total_duration;
    }
    if (const auto minimum_dt = extractNumberField(run_scope, "minimum_time_step_s"); minimum_dt)
    {
        preset.minimum_time_step_s = *minimum_dt;
    }
    if (const auto maximum_dt = extractNumberField(run_scope, "maximum_time_step_s"); maximum_dt)
    {
        preset.maximum_time_step_s = *maximum_dt;
    }

    if (const auto output_csv = extractStringField(content, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else if (const auto output_csv = extractStringField(run_scope, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else
    {
        output_path = preset.default_output_csv;
    }

    auto apply_bool = [&](const std::string& key, bool& target) {
        if (const auto value = extractBoolField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_size_t = [&](const std::string& key, std::size_t& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = static_cast<std::size_t>(std::max(0.0, std::floor(*value + 0.5)));
        }
    };
    auto apply_path = [&](const std::string& key, fs::path& target) {
        if (const auto value = extractStringField(config_scope, key); value)
        {
            target = *value;
        }
    };

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
    apply_number("body_capacitance_f", preset.config.body_capacitance_f);
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

    apply_bool("floating", preset.config.floating);
    apply_bool("derive_capacitance_from_material", preset.config.derive_capacitance_from_material);
    apply_bool("use_reference_current_balance", preset.config.use_reference_current_balance);
    apply_bool("enable_pic_calibration", preset.config.enable_pic_calibration);
    apply_bool("enable_live_pic_window", preset.config.enable_live_pic_window);
    apply_bool("enable_live_pic_mcc", preset.config.enable_live_pic_mcc);
    apply_bool("enable_body_patch_circuit", preset.config.enable_body_patch_circuit);
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

    if (const auto runtime_route = extractStringField(config_scope, "runtime_route");
        runtime_route)
    {
        const auto parsed = parseSurfaceRuntimeRoute(*runtime_route);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported runtime_route in surface config json: " +
                                     *runtime_route);
        }
        preset.config.runtime_route = *parsed;
    }

    if (const auto algorithm_mode = extractStringField(config_scope, "current_algorithm_mode");
        algorithm_mode)
    {
        const auto parsed = parseSurfaceCurrentAlgorithmMode(*algorithm_mode);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported current_algorithm_mode in surface config json: " +
                *algorithm_mode);
        }
        preset.config.current_algorithm_mode = *parsed;
    }

    if (const auto benchmark_mode = extractStringField(config_scope, "benchmark_mode");
        benchmark_mode)
    {
        const auto parsed = parseSurfaceBenchmarkMode(*benchmark_mode);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported benchmark_mode in surface config json: " +
                                     *benchmark_mode);
        }
        preset.config.benchmark_mode = *parsed;
    }

    if (const auto solver_policy = extractStringField(config_scope, "volume_linear_solver_policy");
        solver_policy)
    {
        const auto parsed = parseVolumeLinearSolverPolicy(*solver_policy);
        if (!parsed)
        {
            throw std::runtime_error(
                "Unsupported volume_linear_solver_policy in surface config json: " +
                *solver_policy);
        }
        preset.config.volume_linear_solver_policy = *parsed;
    }

    if (const auto see_model = extractStringField(config_scope, "reference_see_model"); see_model)
    {
        const auto parsed_see_model = parseSecondaryElectronEmissionModel(*see_model);
        if (!parsed_see_model)
        {
            throw std::runtime_error("Unsupported reference_see_model in surface config json: " +
                                     *see_model);
        }
        preset.config.reference_see_model = *parsed_see_model;
    }

    if (const auto model_text = extractStringField(config_scope, "electron_collection_model");
        model_text)
    {
        const auto parsed_model = parseElectronCollectionModelKind(*model_text);
        if (!parsed_model)
        {
            throw std::runtime_error(
                "Unsupported electron_collection_model in surface config json: " + *model_text);
        }
        preset.config.electron_collection_model = *parsed_model;
    }

    auto plasma_text = extractObjectField(config_scope, "plasma_model");
    if (!plasma_text)
    {
        plasma_text = extractObjectField(config_scope, "plasma");
    }
    if (plasma_text)
    {
        applyPlasmaParametersFromJsonObject(*plasma_text, preset.config.plasma);
    }

    auto emission_text = extractObjectField(config_scope, "surface_emission");
    if (!emission_text)
    {
        emission_text = extractObjectField(config_scope, "emission");
    }
    if (emission_text)
    {
        applyEmissionParametersFromJsonObject(*emission_text, preset.config.emission);
    }

    if (const auto material_alias = extractStringField(config_scope, "material_alias");
        material_alias)
    {
        const auto* resolved = resolveMaterialAliasOrName(*material_alias);
        if (resolved == nullptr)
        {
            throw std::runtime_error("Unknown material_alias in surface config json: " +
                                     *material_alias);
        }
        preset.config.material = *resolved;
    }
    if (const auto material_name = extractStringField(config_scope, "material_name");
        material_name)
    {
        const auto* resolved = resolveMaterialAliasOrName(*material_name);
        if (resolved == nullptr)
        {
            throw std::runtime_error("Unknown material_name in surface config json: " +
                                     *material_name);
        }
        preset.config.material = *resolved;
    }

    auto material_text = extractObjectField(config_scope, "surface_material");
    if (!material_text)
    {
        material_text = extractObjectField(config_scope, "material");
    }
    if (material_text)
    {
        applyMaterialFromJsonObject(*material_text, preset.config.material);
    }

    applyStructuredTopologyFromJson(config_scope, preset.config);

    if (const auto spectrum_text = extractObjectField(config_scope, "electron_spectrum"); spectrum_text)
    {
        if (!applySpectrumFromJsonObject(*spectrum_text, Particle::ParticleType::ELECTRON,
                                         preset.config.electron_spectrum))
        {
            throw std::runtime_error("Failed to parse electron_spectrum in surface config json.");
        }
        preset.config.has_electron_spectrum = true;
    }
    if (const auto spectrum_text = extractObjectField(config_scope, "ion_spectrum"); spectrum_text)
    {
        if (!applySpectrumFromJsonObject(*spectrum_text, Particle::ParticleType::ION,
                                         preset.config.ion_spectrum))
        {
            throw std::runtime_error("Failed to parse ion_spectrum in surface config json.");
        }
        preset.config.has_ion_spectrum = true;
    }

    return preset;
}

RadiationScenarioPreset loadRadiationScenarioPresetFromJson(const fs::path& json_path,
                                                            fs::path& output_path)
{
    const std::string content = readTextFile(json_path);
    const auto run_scope = extractObjectField(content, "run").value_or(content);
    const auto config_scope = extractObjectField(content, "config").value_or(content);

    auto preset = SCDAT::Toolkit::Radiation::makeDefaultRadiationScenarioPreset();
    const auto base_preset = extractStringField(content, "base_preset")
                                 .value_or(extractStringField(content, "preset").value_or(""));
    if (!base_preset.empty())
    {
        if (!SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset(base_preset, preset))
        {
            throw std::runtime_error("Unknown radiation base preset in json: " + base_preset);
        }
    }

    if (const auto name = extractStringField(content, "name"); name)
    {
        preset.name = *name;
    }

    if (const auto time_step_s = extractNumberField(run_scope, "time_step_s"); time_step_s)
    {
        preset.time_step_s = *time_step_s;
    }
    if (const auto steps = extractNumberField(run_scope, "steps"); steps)
    {
        preset.steps = static_cast<std::size_t>(std::max(0.0, std::floor(*steps + 0.5)));
    }

    if (const auto output_csv = extractStringField(content, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else if (const auto output_csv = extractStringField(run_scope, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else
    {
        output_path = preset.default_output_csv;
    }

    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_size_t = [&](const std::string& key, std::size_t& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = static_cast<std::size_t>(std::max(0.0, std::floor(*value + 0.5)));
        }
    };
    auto apply_bool = [&](const std::string& key, bool& target) {
        if (const auto value = extractBoolField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_unsigned = [&](const std::string& key, unsigned int& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = static_cast<unsigned int>(std::max(0.0, std::floor(*value + 0.5)));
        }
    };

    apply_size_t("layers", preset.config.layers);
    apply_number("thickness_m", preset.config.thickness_m);
    apply_number("area_m2", preset.config.area_m2);
    apply_number("particle_flux_m2_s", preset.config.particle_flux_m2_s);
    apply_number("mean_energy_ev", preset.config.mean_energy_ev);
    apply_number("attenuation_length_fraction", preset.config.attenuation_length_fraction);
    apply_bool("enable_monte_carlo_transport", preset.config.enable_monte_carlo_transport);
    apply_bool("stopping_data_allow_rollback", preset.config.stopping_data_allow_rollback);
    apply_size_t("monte_carlo_histories_per_step", preset.config.monte_carlo_histories_per_step);
    apply_size_t("monte_carlo_max_steps_per_track", preset.config.monte_carlo_max_steps_per_track);
    apply_size_t("monte_carlo_max_recorded_points", preset.config.monte_carlo_max_recorded_points);
    apply_number("monte_carlo_lateral_span_m", preset.config.monte_carlo_lateral_span_m);
    apply_number("monte_carlo_angular_sigma_deg", preset.config.monte_carlo_angular_sigma_deg);
    apply_unsigned("monte_carlo_seed", preset.config.monte_carlo_seed);

    if (const auto stopping_data_version =
            extractStringField(config_scope, "stopping_data_version");
        stopping_data_version)
    {
        preset.config.stopping_data_version = *stopping_data_version;
    }

    if (const auto material_name = extractStringField(config_scope, "material_name"); material_name)
    {
        const auto* resolved = resolveMaterialAliasOrName(*material_name);
        if (resolved == nullptr)
        {
            throw std::runtime_error(
                "Unsupported material_name in radiation config json: " + *material_name);
        }
        preset.config.material_name = resolved->getName();
    }
    else if (const auto material_name = extractStringField(config_scope, "material"); material_name)
    {
        const auto* resolved = resolveMaterialAliasOrName(*material_name);
        if (resolved == nullptr)
        {
            throw std::runtime_error(
                "Unsupported material in radiation config json: " + *material_name);
        }
        preset.config.material_name = resolved->getName();
    }

    if (const auto species_text = extractStringField(config_scope, "particle_species"); species_text)
    {
        const auto parsed_species = parseRadiationParticleSpecies(*species_text);
        if (!parsed_species)
        {
            throw std::runtime_error("Unsupported particle_species in radiation config json: " +
                                     *species_text);
        }
        preset.config.particle_species = *parsed_species;
    }

    if (const auto physics_list_text = extractStringField(config_scope, "physics_list");
        physics_list_text)
    {
        const auto parsed_physics_list = parseRadiationPhysicsList(*physics_list_text);
        if (!parsed_physics_list)
        {
            throw std::runtime_error("Unsupported physics_list in radiation config json: " +
                                     *physics_list_text);
        }
        preset.config.physics_list = *parsed_physics_list;
    }

    return preset;
}

void applyArcChannelParametersFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::VacuumArc::PlasmaChannelParameters& parameters)
{
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(object_text, key); value)
        {
            target = *value;
        }
    };

    apply_number("field_conductivity_gain", parameters.field_conductivity_gain);
    apply_number("emission_conductivity_gain", parameters.emission_conductivity_gain);
    apply_number("field_current_gain", parameters.field_current_gain);
    apply_number("resistance_scale", parameters.resistance_scale);
    apply_number("min_conductivity_s_per_m", parameters.min_conductivity_s_per_m);
    apply_number("max_conductivity_s_per_m", parameters.max_conductivity_s_per_m);
    apply_number("min_current_density_a_per_m2", parameters.min_current_density_a_per_m2);
    apply_number("max_current_density_a_per_m2", parameters.max_current_density_a_per_m2);
    apply_number("min_resistance_ohm", parameters.min_resistance_ohm);
    apply_number("max_resistance_ohm", parameters.max_resistance_ohm);
}

void applyArcIntegratorParametersFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::VacuumArc::ArcIntegratorParameters& parameters)
{
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(object_text, key); value)
        {
            target = *value;
        }
    };

    if (const auto mode_text = extractStringField(object_text, "mode"); mode_text)
    {
        const auto parsed = parseArcIntegratorMode(*mode_text);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported arc integrator mode in config json: " +
                                     *mode_text);
        }
        parameters.mode = *parsed;
    }

    apply_number("conductivity_growth_per_ns", parameters.conductivity_growth_per_ns);
    apply_number("current_density_growth_per_ns", parameters.current_density_growth_per_ns);
    apply_number("resistance_relaxation_per_ns", parameters.resistance_relaxation_per_ns);
    apply_number("min_conductivity_s_per_m", parameters.min_conductivity_s_per_m);
    apply_number("max_conductivity_s_per_m", parameters.max_conductivity_s_per_m);
    apply_number("min_current_density_a_per_m2", parameters.min_current_density_a_per_m2);
    apply_number("max_current_density_a_per_m2", parameters.max_current_density_a_per_m2);
    apply_number("min_resistance_ohm", parameters.min_resistance_ohm);
    apply_number("max_resistance_ohm", parameters.max_resistance_ohm);
}

void applyArcEmissionParametersFromJsonObject(
    const std::string& object_text,
    SCDAT::Toolkit::VacuumArc::ArcEmissionStrategyParameters& parameters)
{
    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(object_text, key); value)
        {
            target = *value;
        }
    };

    apply_number("field_enhancement_factor", parameters.field_enhancement_factor);
    apply_number("regional_gain", parameters.regional_gain);
    apply_number("thermal_activation_exponent", parameters.thermal_activation_exponent);

    apply_number("linear_field_threshold_v_per_m", parameters.linear_field_threshold_v_per_m);
    apply_number("linear_field_gain", parameters.linear_field_gain);
    apply_number("linear_temperature_threshold_k", parameters.linear_temperature_threshold_k);
    apply_number("linear_temperature_gain", parameters.linear_temperature_gain);

    apply_number("fn_prefactor", parameters.fn_prefactor);
    apply_number("fn_barrier_v_per_m", parameters.fn_barrier_v_per_m);
    apply_number("fn_temperature_reference_k", parameters.fn_temperature_reference_k);
    apply_number("fn_temperature_gain", parameters.fn_temperature_gain);
}

VacuumArcScenarioPreset loadArcScenarioPresetFromJson(const fs::path& json_path,
                                                      fs::path& output_path)
{
    const std::string content = readTextFile(json_path);
    const auto run_scope = extractObjectField(content, "run").value_or(content);
    const auto config_scope = extractObjectField(content, "config").value_or(content);

    auto preset = SCDAT::Toolkit::VacuumArc::makeDefaultVacuumArcScenarioPreset();
    const auto base_preset = extractStringField(content, "base_preset")
                                 .value_or(extractStringField(content, "preset").value_or(""));
    if (!base_preset.empty())
    {
        if (!SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset(base_preset, preset))
        {
            throw std::runtime_error("Unknown arc base preset in json: " + base_preset);
        }
    }

    if (const auto name = extractStringField(content, "name"); name)
    {
        preset.name = *name;
    }

    if (const auto time_step_s = extractNumberField(run_scope, "time_step_s"); time_step_s)
    {
        preset.time_step_s = *time_step_s;
    }
    if (const auto steps = extractNumberField(run_scope, "steps"); steps)
    {
        preset.steps = static_cast<std::size_t>(std::max(0.0, std::floor(*steps + 0.5)));
    }

    if (const auto output_csv = extractStringField(content, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else if (const auto output_csv = extractStringField(run_scope, "output_csv"); output_csv)
    {
        output_path = *output_csv;
    }
    else
    {
        output_path = preset.default_output_csv;
    }

    auto apply_number = [&](const std::string& key, double& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = *value;
        }
    };
    auto apply_size_t = [&](const std::string& key, std::size_t& target) {
        if (const auto value = extractNumberField(config_scope, key); value)
        {
            target = static_cast<std::size_t>(std::max(0.0, std::floor(*value + 0.5)));
        }
    };
    auto apply_bool = [&](const std::string& key, bool& target) {
        if (const auto value = extractBoolField(config_scope, key); value)
        {
            target = *value;
        }
    };

    apply_number("gap_distance_m", preset.config.gap_distance_m);
    apply_number("applied_field_v_per_m", preset.config.applied_field_v_per_m);
    apply_number("surface_potential_v", preset.config.surface_potential_v);
    apply_number("surface_charge_density_c_per_m2", preset.config.surface_charge_density_c_per_m2);
    apply_number("cathode_temperature_k", preset.config.cathode_temperature_k);
    apply_number("anode_temperature_k", preset.config.anode_temperature_k);
    apply_number("channel_radius_m", preset.config.channel_radius_m);
    apply_number("surface_capacitance_f", preset.config.surface_capacitance_f);
    apply_number("surface_leakage_conductance_s", preset.config.surface_leakage_conductance_s);
    apply_number("surface_reference_potential_v", preset.config.surface_reference_potential_v);
    apply_number("max_surface_potential_step_v", preset.config.max_surface_potential_step_v);

    apply_bool("enable_pic_mcc_collisions", preset.config.enable_pic_mcc_collisions);
    apply_size_t("collision_reconfigure_interval_steps",
                 preset.config.collision_reconfigure_interval_steps);
    apply_number("collision_neutral_density_floor_m3", preset.config.collision_neutral_density_floor_m3);
    apply_number("collision_ion_density_floor_m3", preset.config.collision_ion_density_floor_m3);
    apply_number("collision_electron_density_floor_m3", preset.config.collision_electron_density_floor_m3);
    apply_bool("enable_collision_reaction_breakdown",
               preset.config.enable_collision_reaction_breakdown);
    if (const auto cross_section_set_id =
            extractStringField(config_scope, "collision_cross_section_set_id");
        cross_section_set_id)
    {
        preset.config.collision_cross_section_set_id = *cross_section_set_id;
    }
    apply_number("collision_anomaly_event_rate_limit_per_s",
                 preset.config.collision_anomaly_event_rate_limit_per_s);
    apply_bool("collision_fallback_to_background_mcc_on_error",
               preset.config.collision_fallback_to_background_mcc_on_error);
    apply_number("collision_emission_feedback_gain", preset.config.collision_emission_feedback_gain);
    apply_number("collision_channel_feedback_gain", preset.config.collision_channel_feedback_gain);

    apply_number("surface_load_resistance_ohm", preset.config.surface_load_resistance_ohm);
    apply_number("surface_load_inductance_h", preset.config.surface_load_inductance_h);
    apply_number("surface_load_current_limit_a", preset.config.surface_load_current_limit_a);

    apply_bool("enable_neutral_outgassing_feedback",
               preset.config.enable_neutral_outgassing_feedback);
    apply_number("neutral_outgassing_gain_m3_per_a", preset.config.neutral_outgassing_gain_m3_per_a);
    apply_number("neutral_outgassing_relaxation_per_s",
                 preset.config.neutral_outgassing_relaxation_per_s);
    apply_number("neutral_outgassing_max_density_boost_m3",
                 preset.config.neutral_outgassing_max_density_boost_m3);

    apply_bool("enable_secondary_electron_emission", preset.config.enable_secondary_electron_emission);
    apply_number("secondary_electron_yield_per_ion", preset.config.secondary_electron_yield_per_ion);
    apply_number("secondary_electron_speed_m_per_s", preset.config.secondary_electron_speed_m_per_s);

    apply_number("min_internal_substep_s", preset.config.min_internal_substep_s);
    apply_number("max_internal_substep_s", preset.config.max_internal_substep_s);
    apply_bool("enable_stability_rollback", preset.config.enable_stability_rollback);
    apply_size_t("max_stability_rollbacks", preset.config.max_stability_rollbacks);
    apply_number("anomaly_current_density_limit_a_per_m2",
                 preset.config.anomaly_current_density_limit_a_per_m2);
    apply_number("anomaly_surface_potential_limit_v",
                 preset.config.anomaly_surface_potential_limit_v);

    if (const auto alignment = extractStringField(config_scope, "alignment_mode"); alignment)
    {
        const auto parsed = parseArcPicAlignmentMode(*alignment);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported alignment_mode in arc config json: " + *alignment);
        }
        preset.config.alignment_mode = *parsed;
    }

    if (const auto strategy = extractStringField(config_scope, "emission_strategy"); strategy)
    {
        const auto parsed = parseArcEmissionStrategyType(*strategy);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported emission_strategy in arc config json: " + *strategy);
        }
        preset.config.emission_strategy = *parsed;
    }

    if (const auto load_model = extractStringField(config_scope, "surface_load_model"); load_model)
    {
        const auto parsed = parseArcSurfaceLoadModel(*load_model);
        if (!parsed)
        {
            throw std::runtime_error("Unsupported surface_load_model in arc config json: " + *load_model);
        }
        preset.config.surface_load_model = *parsed;
    }

    if (const auto channel_object = extractObjectField(config_scope, "channel_parameters");
        channel_object)
    {
        applyArcChannelParametersFromJsonObject(*channel_object, preset.config.channel_parameters);
    }

    if (const auto integrator_object = extractObjectField(config_scope, "integrator_parameters");
        integrator_object)
    {
        applyArcIntegratorParametersFromJsonObject(*integrator_object,
                                                   preset.config.integrator_parameters);
    }

    if (const auto emission_object = extractObjectField(config_scope, "emission_parameters");
        emission_object)
    {
        applyArcEmissionParametersFromJsonObject(*emission_object,
                                                 preset.config.emission_parameters);
    }

    return preset;
}

RunSelection resolveRunSelection(int argc, char* argv[], int first_optional_index,
                                 const std::string& default_preset,
                                 const fs::path& default_output_path,
                                 const std::string& output_prefix)
{
    RunSelection selection{default_preset, default_output_path};
    if (argc <= first_optional_index)
    {
        return selection;
    }

    const std::string first_argument = argv[first_optional_index];
    if (looksLikeOutputPath(first_argument))
    {
        selection.output_path = first_argument;
        return selection;
    }

    selection.preset_name = first_argument;
    selection.output_path = fs::path("results") / (output_prefix + "_" + selection.preset_name + ".csv");
    if (argc > first_optional_index + 1)
    {
        selection.output_path = argv[first_optional_index + 1];
    }
    return selection;
}

fs::path makeDefaultInternalRadiationOutputPath(const std::string& internal_preset_name,
                                                const std::string& radiation_preset_name)
{
    return fs::path("results") /
           ("internal_radiation_" + internal_preset_name + "__" + radiation_preset_name + ".csv");
}

fs::path makeInternalRadiationCouplingHistoryOutputPath(const fs::path& output_path)
{
    fs::path stem_path = output_path;
    stem_path.replace_extension();
    return fs::path(stem_path.string() + ".coupling_history.csv");
}

InternalRadiationRunSelection resolveInternalRadiationRunSelection(
    int argc, char* argv[], int first_optional_index, const std::string& default_internal_preset,
    const std::string& default_radiation_preset)
{
    InternalRadiationRunSelection selection{
        default_internal_preset,
        default_radiation_preset,
        makeDefaultInternalRadiationOutputPath(default_internal_preset, default_radiation_preset),
    };

    if (argc <= first_optional_index)
    {
        return selection;
    }

    const std::string first_argument = argv[first_optional_index];
    if (looksLikeOutputPath(first_argument))
    {
        selection.output_path = first_argument;
        return selection;
    }

    selection.internal_preset_name = first_argument;
    selection.output_path =
        makeDefaultInternalRadiationOutputPath(selection.internal_preset_name,
                                               selection.radiation_preset_name);

    if (argc <= first_optional_index + 1)
    {
        return selection;
    }

    const std::string second_argument = argv[first_optional_index + 1];
    if (looksLikeOutputPath(second_argument))
    {
        selection.output_path = second_argument;
        return selection;
    }

    selection.radiation_preset_name = second_argument;
    selection.output_path =
        makeDefaultInternalRadiationOutputPath(selection.internal_preset_name,
                                               selection.radiation_preset_name);

    if (argc > first_optional_index + 2)
    {
        selection.output_path = argv[first_optional_index + 2];
    }
    return selection;
}

double mapDepositedEnergyToCurrentDensity(double deposited_energy_j_per_m2, double dt,
                                          double mean_energy_ev,
                                          double incident_charge_state_abs)
{
    if (deposited_energy_j_per_m2 <= 0.0 || dt <= 0.0 || mean_energy_ev <= 0.0 ||
        incident_charge_state_abs <= 0.0)
    {
        return 0.0;
    }

    constexpr double kEvToJ = 1.602176634e-19;
    constexpr double kElementaryChargeC = 1.602176634e-19;
    const double particle_energy_j = mean_energy_ev * kEvToJ;
    const double particle_flux_m2_s = deposited_energy_j_per_m2 / (dt * particle_energy_j);
    const double particle_charge_c = incident_charge_state_abs * kElementaryChargeC;
    return particle_flux_m2_s * particle_charge_c;
}

double mapRadiationSpeciesToChargeStateAbs(SCDAT::Toolkit::Radiation::ParticleSpecies species)
{
    using SCDAT::Toolkit::Radiation::ParticleSpecies;
    switch (species)
    {
    case ParticleSpecies::Electron:
        return 1.0;
    case ParticleSpecies::Proton:
        return 1.0;
    case ParticleSpecies::HeavyIon:
        return 2.0;
    }
    return 1.0;
}

double computeInternalToRadiationFeedbackScale(
    const SCDAT::Toolkit::InternalCharging::InternalChargingStatus& internal_status)
{
    const double field_shielding_term =
        1.0 / (1.0 + std::abs(internal_status.max_electric_field_v_per_m) / 1.0e8);
    const double conductivity_relief_term =
        1.0 + std::clamp(internal_status.effective_conductivity_s_per_m / 5.0e-8, 0.0, 0.15);
    return std::clamp(field_shielding_term * conductivity_relief_term, 0.85, 1.15);
}

bool exportInternalRadiationCouplingHistory(
    const fs::path& csv_path,
    const std::vector<InternalRadiationCouplingHistoryRow>& rows)
{
    std::ofstream output(csv_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        return false;
    }

    output << "macro_step,substep,iteration,relative_error,absolute_error,deposited_energy_delta_j_per_m2,incident_current_density_a_per_m2,feedback_scale,internal_max_electric_field_v_per_m,internal_average_dose_gy,converged\n";
    for (const auto& row : rows)
    {
        output << row.macro_step << ',' << row.substep << ',' << row.iteration << ','
               << row.relative_error << ',' << row.absolute_error << ','
               << row.deposited_energy_delta_j_per_m2 << ','
               << row.incident_current_density_a_per_m2 << ',' << row.feedback_scale << ','
               << row.internal_max_electric_field_v_per_m << ','
               << row.internal_average_dose_gy << ',' << (row.converged ? 1 : 0) << '\n';
    }

    return static_cast<bool>(output);
}

void printPresetNames(const std::string& module, const std::vector<std::string>& names)
{
    std::cout << module << " presets\n";
    for (const auto& name : names)
    {
        std::cout << "  " << name << '\n';
    }
}

void printUsage(const std::string& exe_name)
{
    std::cout << "Spacecraft Charging and Discharging Analysis Toolkit\n\n";
    std::cout << "Usage:\n";
    std::cout << "  " << exe_name << " help\n";
    std::cout << "  " << exe_name << " presets [plasma|surface|internal|radiation|arc]\n";
    std::cout << "  " << exe_name << " status\n";
    std::cout << "  " << exe_name << " summary [result_csv]\n";
    std::cout << "  " << exe_name << " plasma [preset] [output_csv]\n";
    std::cout << "  " << exe_name << " plasma-config <config_json> [output_csv]\n";
    std::cout << "  " << exe_name << " surface [preset] [output_csv]\n";
    std::cout << "  " << exe_name << " surface-config <config_json> [output_csv]\n";
    std::cout << "  " << exe_name << " internal [preset] [output_csv]\n";
    std::cout << "  " << exe_name
              << " internal-radiation [internal_preset] [radiation_preset] [output_csv]\n";
    std::cout << "  " << exe_name << " radiation [preset] [output_csv]\n";
    std::cout << "  " << exe_name << " radiation-config <config_json> [output_csv]\n";
    std::cout << "  " << exe_name << " arc [preset] [output_csv]\n\n";
    std::cout << "  " << exe_name << " arc-config <config_json> [output_csv]\n\n";
    std::cout << "Default presets:\n";
    std::cout << "  plasma  : ccp_argon_10pa\n";
    std::cout << "  surface : leo_daylight_kapton\n";
    std::cout << "  internal: geo_electron_belt\n";
    std::cout << "  radiation: geo_electron_belt_dose\n";
    std::cout << "  arc     : microgap_flashover\n\n";
    std::cout << "Benchmark entry points:\n";
    std::cout << "  SCDAT_PICcore                 Run the PIC-MCC reference benchmark\n";
    std::cout << "  SCDAT_mesh_export            Export mesh csv files for inspection\n";
}

void printSummary(const DataSetSummary& summary, const fs::path& csv_path)
{
    std::cout << "Result summary: " << csv_path.string() << '\n';
    std::cout << "  rows              : " << summary.rows << '\n';
    std::cout << "  axis              : " << summary.axis_name << '\n';

    for (const auto& [name, column] : summary.column_summaries)
    {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  " << std::left << std::setw(18) << name << ": " << column.minimum << " .. "
                  << column.maximum << "  mean=" << column.mean << '\n';
    }
}

void printStatus()
{
    const fs::path results_dir("results");
    DataAnalyzer analyzer;

    std::cout << "Project status\n";
    std::cout << "  Results directory: " << fs::absolute(results_dir).string() << '\n';
    std::cout << "  Results exists   : " << (fs::exists(results_dir) ? "yes" : "no") << '\n';

    const auto latest_csv = analyzer.findLatestResultFile(results_dir, ".csv");
    const auto latest_png = analyzer.findLatestResultFile(results_dir, ".png");
    const auto latest_vtk = analyzer.findLatestResultFile(results_dir, ".vtk");
    const auto latest_h5 = analyzer.findLatestResultFile(results_dir, ".h5");

    std::cout << "  Latest csv       : "
              << (latest_csv ? latest_csv->string() : std::string("<none>")) << '\n';
    std::cout << "  Latest png       : "
              << (latest_png ? latest_png->string() : std::string("<none>")) << '\n';
    std::cout << "  Latest vtk       : "
              << (latest_vtk ? latest_vtk->string() : std::string("<none>")) << '\n';
    std::cout << "  Latest h5        : "
              << (latest_h5 ? latest_h5->string() : std::string("<none>")) << '\n';

    if (latest_csv)
    {
        const auto summary = analyzer.summarizeCsv(*latest_csv);
        std::cout << "  Latest csv rows  : " << summary.rows << '\n';
        std::cout << "  Axis name        : " << summary.axis_name << '\n';
    }
}

int runPlasma(const PlasmaScenarioPreset& preset, const fs::path& output_path)
{
    PICFluidIntegration integration;
    if (!integration.initialize(preset.config))
    {
        throw std::runtime_error("Failed to initialize plasma analysis toolkit");
    }

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        if (!integration.advance(preset.config.time_step_s))
        {
            throw std::runtime_error("Plasma analysis advance failed");
        }
    }

    if (!integration.exportResults(output_path))
    {
        throw std::runtime_error("Failed to export plasma analysis results");
    }

    std::cout << "Plasma preset '" << preset.name << "' wrote results to " << output_path.string()
              << '\n';
    return 0;
}

int runSurface(const SurfaceChargingScenarioPreset& preset, const fs::path& output_path)
{
    DensePlasmaSurfaceCharging charging;
    if (!charging.initialize(preset.config))
    {
        throw std::runtime_error("Failed to initialize surface charging toolkit");
    }

    if (preset.adaptive_time_stepping && preset.total_duration_s > 0.0)
    {
        double elapsed_time_s = 0.0;
        std::size_t sample_count = 0;
        const double minimum_dt_s = std::max(1.0e-12, preset.minimum_time_step_s);
        const double maximum_dt_s =
            std::max(minimum_dt_s, preset.maximum_time_step_s);
        while (elapsed_time_s + 1.0e-12 < preset.total_duration_s)
        {
            const double remaining_time_s = preset.total_duration_s - elapsed_time_s;
            double dt = 0.0;
            if (preset.steps > 0 && sample_count >= preset.steps)
            {
                // Avoid a final oversized tail step after the adaptive budget is
                // exhausted; finish the horizon with bounded chunks instead.
                dt = std::min(remaining_time_s, maximum_dt_s);
            }
            else
            {
                dt = charging.recommendTimeStep(remaining_time_s, minimum_dt_s, maximum_dt_s);
            }
            if (!charging.advance(dt))
            {
                throw std::runtime_error("Surface charging adaptive advance failed");
            }
            elapsed_time_s += dt;
            sample_count += 1;
        }
    }
    else
    {
        for (std::size_t i = 0; i < preset.steps; ++i)
        {
            if (!charging.advance(preset.time_step_s))
            {
                throw std::runtime_error("Surface charging advance failed");
            }
        }
    }

    if (!charging.exportResults(output_path))
    {
        throw std::runtime_error("Failed to export surface charging results");
    }

    std::cout << "Surface preset '" << preset.name << "' wrote results to " << output_path.string()
              << '\n';
    return 0;
}

int runInternal(const InternalChargingScenarioPreset& preset, const fs::path& output_path)
{
    SpacecraftInternalChargingAlgorithm algorithm;
    if (!algorithm.initialize(preset.config))
    {
        throw std::runtime_error("Failed to initialize internal charging toolkit");
    }

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        if (!algorithm.advance(preset.time_step_s))
        {
            throw std::runtime_error("Internal charging advance failed");
        }
    }

    if (!algorithm.exportResults(output_path))
    {
        throw std::runtime_error("Failed to export internal charging results");
    }

    std::cout << "Internal preset '" << preset.name << "' wrote results to " << output_path.string()
              << '\n';
    return 0;
}

int runInternalRadiationCoupled(const InternalChargingScenarioPreset& internal_preset,
                                const RadiationScenarioPreset& radiation_preset,
                                const fs::path& output_path)
{
    constexpr std::size_t kCouplingSubstepsPerStep = 4;
    constexpr int kCouplingIterationsPerSubstep = 4;
    constexpr double kCouplingRelativeTolerance = 2.0;
    constexpr double kCouplingAbsoluteTolerance = 1.0e-3;

    RadiationDoseAlgorithm radiation_algorithm;
    if (!radiation_algorithm.initialize(radiation_preset.config))
    {
        throw std::runtime_error("Failed to initialize radiation toolkit for internal-radiation coupling");
    }

    auto internal_config = internal_preset.config;
    internal_config.source_mode = InternalChargingSourceMode::Radiation;

    SpacecraftInternalChargingAlgorithm internal_algorithm;
    if (!internal_algorithm.initialize(internal_config))
    {
        throw std::runtime_error(
            "Failed to initialize internal charging toolkit for internal-radiation coupling");
    }

    const double dt = internal_preset.time_step_s;
    if (dt <= 0.0)
    {
        throw std::runtime_error("internal-radiation coupling requires positive internal preset time_step_s");
    }

    const double substep_dt = dt / static_cast<double>(kCouplingSubstepsPerStep);
    const double iteration_dt =
        substep_dt / static_cast<double>(std::max(1, kCouplingIterationsPerSubstep));
    if (iteration_dt <= 0.0)
    {
        throw std::runtime_error("internal-radiation coupling requires positive iterative substep dt");
    }

    double deposited_energy_previous_j_per_m2 = 0.0;
    double previous_drive_current_density_a_per_m2 = 0.0;
    double current_relative_error = 0.0;
    double current_absolute_error = 0.0;
    double current_deposited_energy_delta_j_per_m2 = 0.0;
    double current_drive_current_density_a_per_m2 = 0.0;
    double current_feedback_scale = 1.0;
    std::size_t active_macro_step = 0;
    std::size_t active_substep = 0;
    int active_iteration = 0;
    std::vector<InternalRadiationCouplingHistoryRow> coupling_history;
    coupling_history.reserve(internal_preset.steps * kCouplingSubstepsPerStep *
                             static_cast<std::size_t>(kCouplingIterationsPerSubstep));

    auto coupling = std::make_shared<FunctionalCoupling>(
        "radiation_internal_iterative", CouplingType::ITERATIVE,
        [&](double step_dt) {
            ++active_iteration;

            if (!radiation_algorithm.advance(step_dt))
            {
                return SCDAT::VoidResult::failure("Radiation advance failed in one-way coupling step");
            }

            const auto& radiation_status = radiation_algorithm.getStatus();
            current_deposited_energy_delta_j_per_m2 =
                std::max(0.0, radiation_status.deposited_energy_j_per_m2 - deposited_energy_previous_j_per_m2);
            deposited_energy_previous_j_per_m2 = radiation_status.deposited_energy_j_per_m2;

            const auto& internal_status_before = internal_algorithm.getStatus();
            current_feedback_scale = computeInternalToRadiationFeedbackScale(internal_status_before);

            InternalChargingRadiationDrive drive;
            drive.incident_energy_ev = std::max(1.0, radiation_preset.config.mean_energy_ev);
            drive.incident_charge_state_abs =
                mapRadiationSpeciesToChargeStateAbs(radiation_preset.config.particle_species);
            const double base_drive_current_density_a_per_m2 =
                mapDepositedEnergyToCurrentDensity(current_deposited_energy_delta_j_per_m2, step_dt,
                                                   drive.incident_energy_ev,
                                                   drive.incident_charge_state_abs);
            current_drive_current_density_a_per_m2 =
                base_drive_current_density_a_per_m2 * current_feedback_scale;
            drive.incident_current_density_a_per_m2 = current_drive_current_density_a_per_m2;
            internal_algorithm.setRadiationDrive(drive);

            if (!internal_algorithm.advance(step_dt))
            {
                return SCDAT::VoidResult::failure(
                    "Internal charging advance failed in one-way coupling step");
            }

            current_absolute_error = std::abs(
                current_drive_current_density_a_per_m2 - previous_drive_current_density_a_per_m2);
            const double relative_denominator = std::max(
                std::max(std::abs(previous_drive_current_density_a_per_m2),
                         std::abs(current_drive_current_density_a_per_m2)),
                1.0e-12);
            current_relative_error = current_absolute_error / relative_denominator;
            previous_drive_current_density_a_per_m2 = current_drive_current_density_a_per_m2;

            const bool converged =
                current_relative_error <= kCouplingRelativeTolerance &&
                current_absolute_error <= kCouplingAbsoluteTolerance;
            const auto& internal_status = internal_algorithm.getStatus();
            coupling_history.push_back(InternalRadiationCouplingHistoryRow{
                active_macro_step,
                active_substep,
                active_iteration,
                current_relative_error,
                current_absolute_error,
                current_deposited_energy_delta_j_per_m2,
                current_drive_current_density_a_per_m2,
                current_feedback_scale,
                internal_status.max_electric_field_v_per_m,
                internal_status.average_dose_gy,
                converged,
            });

            return SCDAT::VoidResult::success();
        },
        [&]() { return std::pair<double, double>{current_relative_error, current_absolute_error}; });

    SCDAT::Coupling::ConvergenceCriteria criteria;
    criteria.relative_tolerance = kCouplingRelativeTolerance;
    criteria.absolute_tolerance = kCouplingAbsoluteTolerance;
    criteria.max_iterations = kCouplingIterationsPerSubstep;
    criteria.min_iterations = kCouplingIterationsPerSubstep;
    coupling->setConvergenceCriteria(criteria);

    MultiPhysicsManager manager;
    manager.addCoupling(coupling);
    for (std::size_t i = 0; i < internal_preset.steps; ++i)
    {
        for (std::size_t substep = 0; substep < kCouplingSubstepsPerStep; ++substep)
        {
            active_macro_step = i;
            active_substep = substep;
            active_iteration = 0;

            if (const auto result = manager.executeIterativeCouplings(iteration_dt); !result)
            {
                throw std::runtime_error("internal-radiation coupling failed: " + result.message());
            }

            const auto stats = manager.getAllStatistics();
            const auto stats_it = stats.find("radiation_internal_iterative");
            if (stats_it == stats.end())
            {
                throw std::runtime_error("internal-radiation coupling failed: missing iterative statistics");
            }
            if (!stats_it->second.converged)
            {
                throw std::runtime_error(
                    "internal-radiation coupling failed: iterative substep did not converge");
            }
        }
    }

    if (!internal_algorithm.exportResults(output_path))
    {
        throw std::runtime_error("Failed to export internal-radiation coupled results");
    }

    const auto coupling_history_path = makeInternalRadiationCouplingHistoryOutputPath(output_path);
    if (!exportInternalRadiationCouplingHistory(coupling_history_path, coupling_history))
    {
        throw std::runtime_error("Failed to export internal-radiation coupling history");
    }

    std::cout << "Internal-radiation presets '" << internal_preset.name << "' + '"
              << radiation_preset.name << "' wrote results to " << output_path.string()
              << " and coupling history to " << coupling_history_path.string() << '\n';
    return 0;
}

int runRadiation(const RadiationScenarioPreset& preset, const fs::path& output_path)
{
    RadiationDoseAlgorithm algorithm;
    if (!algorithm.initialize(preset.config))
    {
        throw std::runtime_error("Failed to initialize radiation toolkit");
    }

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        if (!algorithm.advance(preset.time_step_s))
        {
            throw std::runtime_error("Radiation advance failed");
        }
    }

    if (!algorithm.exportResults(output_path))
    {
        throw std::runtime_error("Failed to export radiation results");
    }

    std::cout << "Radiation preset '" << preset.name << "' wrote results to " << output_path.string()
              << '\n';
    return 0;
}

int runArc(const VacuumArcScenarioPreset& preset, const fs::path& output_path)
{
    SurfaceDischargeArcAlgorithm algorithm;
    if (!algorithm.initialize(preset.config))
    {
        throw std::runtime_error("Failed to initialize vacuum arc toolkit");
    }

    for (std::size_t i = 0; i < preset.steps; ++i)
    {
        if (!algorithm.advance(preset.time_step_s))
        {
            throw std::runtime_error("Vacuum arc advance failed");
        }
    }

    if (!algorithm.exportResults(output_path))
    {
        throw std::runtime_error("Failed to export vacuum arc results");
    }

    std::cout << "Vacuum arc preset '" << preset.name << "' wrote results to "
              << output_path.string() << '\n';
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    try
    {
        const std::string exe_name = (argc > 0) ? fs::path(argv[0]).filename().string() : "SCDAT";
        const std::string command = (argc > 1) ? argv[1] : "help";
        DataAnalyzer analyzer;

        if (command == "help" || command == "--help" || command == "-h")
        {
            printUsage(exe_name);
            return 0;
        }

        if (command == "status")
        {
            printStatus();
            return 0;
        }

        if (command == "presets")
        {
            const std::string module = (argc > 2) ? argv[2] : "all";
            if (module == "all" || module == "plasma")
            {
                printPresetNames("plasma", SCDAT::Toolkit::PlasmaAnalysis::listPlasmaScenarioPresetNames());
            }
            if (module == "all" || module == "surface")
            {
                printPresetNames("surface",
                                 SCDAT::Toolkit::SurfaceCharging::listSurfaceChargingScenarioPresetNames());
            }
            if (module == "all" || module == "internal")
            {
                printPresetNames(
                    "internal", SCDAT::Toolkit::InternalCharging::listInternalChargingScenarioPresetNames());
            }
            if (module == "all" || module == "radiation")
            {
                printPresetNames("radiation", SCDAT::Toolkit::Radiation::listRadiationScenarioPresetNames());
            }
            if (module == "all" || module == "arc")
            {
                printPresetNames("arc", SCDAT::Toolkit::VacuumArc::listVacuumArcScenarioPresetNames());
            }

            if (module != "all" && module != "plasma" && module != "surface" &&
                module != "internal" && module != "radiation" && module != "arc")
            {
                throw std::runtime_error("Unknown module for presets command: " + module);
            }
            return 0;
        }

        if (command == "summary")
        {
            fs::path csv_path;
            if (argc > 2)
            {
                csv_path = argv[2];
            }
            else
            {
                const auto latest_csv = analyzer.findLatestResultFile("results", ".csv");
                if (!latest_csv)
                {
                    throw std::runtime_error("No csv file found under results/.");
                }
                csv_path = *latest_csv;
            }

            printSummary(analyzer.summarizeCsv(csv_path), csv_path);
            return 0;
        }

        std::filesystem::create_directories("results");

        if (command == "plasma")
        {
            auto preset = SCDAT::Toolkit::PlasmaAnalysis::makeDefaultPlasmaScenarioPreset();
            const auto selection =
                resolveRunSelection(argc, argv, 2, preset.name, preset.default_output_csv, "plasma");
            if (!SCDAT::Toolkit::PlasmaAnalysis::tryGetPlasmaScenarioPreset(selection.preset_name, preset))
            {
                throw std::runtime_error("Unknown plasma preset: " + selection.preset_name);
            }
            return runPlasma(preset, selection.output_path);
        }
        if (command == "plasma-config")
        {
            if (argc < 3)
            {
                throw std::runtime_error(
                    "plasma-config command requires a JSON file path: plasma-config <config_json> [output_csv]");
            }

            fs::path output_path;
            auto preset = loadPlasmaScenarioPresetFromJson(argv[2], output_path);
            if (argc > 3)
            {
                output_path = argv[3];
            }
            return runPlasma(preset, output_path);
        }
        if (command == "surface")
        {
            auto preset = SCDAT::Toolkit::SurfaceCharging::makeDefaultSurfaceChargingScenarioPreset();
            const auto selection =
                resolveRunSelection(argc, argv, 2, preset.name, preset.default_output_csv, "surface");
            if (!SCDAT::Toolkit::SurfaceCharging::tryGetSurfaceChargingScenarioPreset(
                    selection.preset_name, preset))
            {
                throw std::runtime_error("Unknown surface preset: " + selection.preset_name);
            }
            return runSurface(preset, selection.output_path);
        }
        if (command == "surface-config")
        {
            if (argc < 3)
            {
                throw std::runtime_error(
                    "surface-config command requires a JSON file path: surface-config <config_json> [output_csv]");
            }

            fs::path output_path;
            auto preset = loadSurfaceScenarioPresetFromJson(argv[2], output_path);
            if (argc > 3)
            {
                output_path = argv[3];
            }
            return runSurface(preset, output_path);
        }
        if (command == "internal")
        {
            auto preset = SCDAT::Toolkit::InternalCharging::makeDefaultInternalChargingScenarioPreset();
            const auto selection =
                resolveRunSelection(argc, argv, 2, preset.name, preset.default_output_csv, "internal");
            if (!SCDAT::Toolkit::InternalCharging::tryGetInternalChargingScenarioPreset(
                    selection.preset_name, preset))
            {
                throw std::runtime_error("Unknown internal preset: " + selection.preset_name);
            }
            return runInternal(preset, selection.output_path);
        }
        if (command == "internal-radiation")
        {
            auto internal_preset =
                SCDAT::Toolkit::InternalCharging::makeDefaultInternalChargingScenarioPreset();
            auto radiation_preset = SCDAT::Toolkit::Radiation::makeDefaultRadiationScenarioPreset();
            const auto selection =
                resolveInternalRadiationRunSelection(argc, argv, 2, internal_preset.name,
                                                     radiation_preset.name);
            if (!SCDAT::Toolkit::InternalCharging::tryGetInternalChargingScenarioPreset(
                    selection.internal_preset_name, internal_preset))
            {
                throw std::runtime_error("Unknown internal preset: " + selection.internal_preset_name);
            }
            if (!SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset(
                    selection.radiation_preset_name, radiation_preset))
            {
                throw std::runtime_error("Unknown radiation preset: " + selection.radiation_preset_name);
            }
            return runInternalRadiationCoupled(internal_preset, radiation_preset,
                                               selection.output_path);
        }
        if (command == "radiation")
        {
            auto preset = SCDAT::Toolkit::Radiation::makeDefaultRadiationScenarioPreset();
            const auto selection =
                resolveRunSelection(argc, argv, 2, preset.name, preset.default_output_csv, "radiation");
            if (!SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset(selection.preset_name, preset))
            {
                throw std::runtime_error("Unknown radiation preset: " + selection.preset_name);
            }
            return runRadiation(preset, selection.output_path);
        }
        if (command == "radiation-config")
        {
            if (argc < 3)
            {
                throw std::runtime_error(
                    "radiation-config command requires a JSON file path: radiation-config <config_json> [output_csv]");
            }

            fs::path output_path;
            auto preset = loadRadiationScenarioPresetFromJson(argv[2], output_path);
            if (argc > 3)
            {
                output_path = argv[3];
            }
            return runRadiation(preset, output_path);
        }
        if (command == "arc")
        {
            auto preset = SCDAT::Toolkit::VacuumArc::makeDefaultVacuumArcScenarioPreset();
            const auto selection =
                resolveRunSelection(argc, argv, 2, preset.name, preset.default_output_csv, "arc");
            if (!SCDAT::Toolkit::VacuumArc::tryGetVacuumArcScenarioPreset(selection.preset_name, preset))
            {
                throw std::runtime_error("Unknown arc preset: " + selection.preset_name);
            }
            return runArc(preset, selection.output_path);
        }
        if (command == "arc-config")
        {
            if (argc < 3)
            {
                throw std::runtime_error(
                    "arc-config command requires a JSON file path: arc-config <config_json> [output_csv]");
            }

            fs::path output_path;
            auto preset = loadArcScenarioPresetFromJson(argv[2], output_path);
            if (argc > 3)
            {
                output_path = argv[3];
            }
            return runArc(preset, output_path);
        }

        std::cerr << "Unknown command: " << command << "\n\n";
        printUsage(exe_name);
        return 1;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "SCDAT failed: " << ex.what() << '\n';
        return 1;
    }
}
