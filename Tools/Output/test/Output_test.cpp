#include "DataAnalyzer.h"
#include "ResultExporter.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <limits>

using SCDAT::Output::ColumnarDataSet;
using SCDAT::Output::DataAnalyzer;
using SCDAT::Output::ExportArtifacts;
using SCDAT::Output::ResultExporter;

TEST(OutputTest, ResultExporterWritesCsvVtkAndHdf5Artifacts)
{
    ColumnarDataSet data_set;
    data_set.axis_name = "z_m";
    data_set.axis_values = {0.0, 0.5, 1.0};
    data_set.scalar_series["potential_v"] = {0.0, 1.0, 2.0};
    data_set.scalar_series["electric_field_z_v_per_m"] = {2.0, 1.0, 0.0};
    data_set.scalar_series["electron_density_m3"] = {1.0e15, 1.1e15, 1.2e15};
    data_set.scalar_series["ion_density_m3"] = {1.0e15, 1.0e15, 1.0e15};
    data_set.scalar_series["electron_temperature_ev"] = {3.0, 3.0, 3.0};

    const auto csv_path = std::filesystem::temp_directory_path() / "scdat_output_test.csv";
    ResultExporter exporter;
    ExportArtifacts artifacts;
    ASSERT_TRUE(exporter.exportDataSet(csv_path, data_set, &artifacts));

    EXPECT_TRUE(std::filesystem::exists(artifacts.csv_path));
    EXPECT_TRUE(std::filesystem::exists(artifacts.vtk_path));
    EXPECT_TRUE(std::filesystem::exists(artifacts.hdf5_path));
}

TEST(OutputTest, DataAnalyzerSummarizesCsvColumns)
{
    ColumnarDataSet data_set;
    data_set.axis_name = "z_m";
    data_set.axis_values = {0.0, 1.0};
    data_set.scalar_series["potential_v"] = {-1.0, 3.0};
    data_set.scalar_series["electric_field_z_v_per_m"] = {5.0, 7.0};
    data_set.scalar_series["electron_density_m3"] = {1.0, 2.0};
    data_set.scalar_series["ion_density_m3"] = {2.0, 2.0};
    data_set.scalar_series["electron_temperature_ev"] = {4.0, 6.0};

    const auto csv_path = std::filesystem::temp_directory_path() / "scdat_output_summary.csv";
    ResultExporter exporter;
    ASSERT_TRUE(exporter.exportDataSet(csv_path, data_set));

    DataAnalyzer analyzer;
    const auto summary = analyzer.summarizeCsv(csv_path);
    EXPECT_EQ(summary.rows, 2u);
    EXPECT_EQ(summary.axis_name, "z_m");
    EXPECT_NEAR(summary.column_summaries.at("potential_v").minimum, -1.0, 1.0e-12);
    EXPECT_NEAR(summary.column_summaries.at("potential_v").maximum, 3.0, 1.0e-12);
    EXPECT_NEAR(summary.column_summaries.at("electron_temperature_ev").mean, 5.0, 1.0e-12);
}

TEST(OutputTest, ResultExporterRejectsNonFiniteData)
{
    ColumnarDataSet data_set;
    data_set.axis_name = "time_s";
    data_set.axis_values = {0.0, 1.0};
    data_set.scalar_series["surface_potential_v"] = {0.0, std::numeric_limits<double>::infinity()};

    const auto csv_path = std::filesystem::temp_directory_path() / "scdat_output_bad.csv";
    std::filesystem::remove(csv_path);
    ResultExporter exporter;
    const auto result = exporter.exportDataSet(csv_path, data_set);
    EXPECT_FALSE(result);
    EXPECT_FALSE(std::filesystem::exists(csv_path));
}
