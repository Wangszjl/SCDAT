#include "DataAnalyzer.h"
#include "DensePlasmaSurfaceCharging.h"
#include "InternalChargingCases.h"
#include "PICFluidIntegration.h"
#include "PlasmaAnalysisCases.h"
#include "SpacecraftInternalChargingAlgorithm.h"
#include "SurfaceChargingCases.h"
#include "SurfaceDischargeArcAlgorithm.h"
#include "VacuumArcCases.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

namespace fs = std::filesystem;

using SCDAT::Output::DataAnalyzer;
using SCDAT::Output::DataSetSummary;
using SCDAT::Toolkit::InternalCharging::InternalChargingScenarioPreset;
using SCDAT::Toolkit::InternalCharging::SpacecraftInternalChargingAlgorithm;
using SCDAT::Toolkit::PlasmaAnalysis::PlasmaScenarioPreset;
using SCDAT::Toolkit::PlasmaAnalysis::PICFluidIntegration;
using SCDAT::Toolkit::SurfaceCharging::DensePlasmaSurfaceCharging;
using SCDAT::Toolkit::SurfaceCharging::SurfaceChargingScenarioPreset;
using SCDAT::Toolkit::VacuumArc::VacuumArcScenarioPreset;
using SCDAT::Toolkit::VacuumArc::SurfaceDischargeArcAlgorithm;

struct RunSelection
{
    std::string preset_name;
    fs::path output_path;
};

bool looksLikeOutputPath(const std::string& argument)
{
    const auto extension = fs::path(argument).extension().string();
    return extension == ".csv" || extension == ".vtk" || extension == ".h5" || extension == ".png";
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
    std::cout << "  " << exe_name << " presets [plasma|surface|internal|arc]\n";
    std::cout << "  " << exe_name << " status\n";
    std::cout << "  " << exe_name << " summary [result_csv]\n";
    std::cout << "  " << exe_name << " plasma [preset] [output_csv]\n";
    std::cout << "  " << exe_name << " surface [preset] [output_csv]\n";
    std::cout << "  " << exe_name << " internal [preset] [output_csv]\n";
    std::cout << "  " << exe_name << " arc [preset] [output_csv]\n\n";
    std::cout << "Default presets:\n";
    std::cout << "  plasma  : ccp_argon_10pa\n";
    std::cout << "  surface : leo_daylight_kapton\n";
    std::cout << "  internal: geo_electron_belt\n";
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
        while (elapsed_time_s < preset.total_duration_s && sample_count < preset.steps)
        {
            const double remaining_time_s = preset.total_duration_s - elapsed_time_s;
            const double dt = charging.recommendTimeStep(
                remaining_time_s, std::max(1.0e-12, preset.minimum_time_step_s),
                std::max(preset.minimum_time_step_s, preset.maximum_time_step_s));
            if (!charging.advance(dt))
            {
                throw std::runtime_error("Surface charging adaptive advance failed");
            }
            elapsed_time_s += dt;
            sample_count += 1;
        }

        if (elapsed_time_s + 1.0e-12 < preset.total_duration_s)
        {
            const double tail_dt = preset.total_duration_s - elapsed_time_s;
            if (!charging.advance(tail_dt))
            {
                throw std::runtime_error("Surface charging final adaptive advance failed");
            }
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
            if (module == "all" || module == "arc")
            {
                printPresetNames("arc", SCDAT::Toolkit::VacuumArc::listVacuumArcScenarioPresetNames());
            }

            if (module != "all" && module != "plasma" && module != "surface" &&
                module != "internal" && module != "arc")
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
