#include "MaterialDatabase.h"
#include "SPISMaterialLoader.h"
#include "SurfaceInteraction.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

using SCDAT::Material::MaterialDatabase;
using SCDAT::Material::SPISMaterialLoader;
using SCDAT::Material::SurfaceImpactState;
using SCDAT::Material::SurfaceInteraction;

TEST(MaterialDatabaseTest, DefaultLibraryContainsKapton)
{
    MaterialDatabase database;
    const auto* kapton = database.findByName("kapton");
    ASSERT_NE(kapton, nullptr);
    EXPECT_TRUE(kapton->isDielectric());
    EXPECT_GT(kapton->getBreakdownFieldVPerM(), 0.0);
}

TEST(MaterialDatabaseTest, UnifiedDatabaseResolvesCommonAliases)
{
    SCDAT::Material::UnifiedMaterialDatabase database;

    const auto* kapton = database.findByAliasOrName("kapton_hn");
    ASSERT_NE(kapton, nullptr);
    EXPECT_EQ(kapton->getName(), "kapton");

    const auto* aluminum = database.findByAliasOrName("al6061");
    ASSERT_NE(aluminum, nullptr);
    EXPECT_EQ(aluminum->getName(), "aluminum");

    const auto* ptfe = database.findByAliasOrName("teflon");
    ASSERT_NE(ptfe, nullptr);
    EXPECT_EQ(ptfe->getName(), "ptfe");
}

TEST(MaterialDatabaseTest, LoaderImportsCsvMaterialRecords)
{
    const auto path = std::filesystem::temp_directory_path() / "scdat_material_test.csv";
    std::ofstream output(path);
    output << "id,name,type,permittivity,conductivity,breakdown_field_v_per_m\n";
    output << "10,Quartz,dielectric,3.8,1e-18,9.0e7\n";
    output.close();

    MaterialDatabase database;
    SPISMaterialLoader loader;
    ASSERT_TRUE(loader.loadCsv(path, database));

    const auto* quartz = database.findByName("quartz");
    ASSERT_NE(quartz, nullptr);
    EXPECT_NEAR(quartz->getPermittivity(), 3.8, 1.0e-12);
    EXPECT_NEAR(quartz->getBreakdownFieldVPerM(), 9.0e7, 1.0);
}

TEST(SurfaceInteractionTest, InteractionDepositsChargeAndHeat)
{
    MaterialDatabase database;
    const auto* kapton = database.findByName("kapton");
    ASSERT_NE(kapton, nullptr);

    SurfaceInteraction interaction;
    SurfaceImpactState impact;
    impact.incident_energy_ev = 50.0;
    impact.particle_charge_coulomb = -1.602176634e-19;
    impact.incident_angle_rad = 0.2;

    const auto result = interaction.evaluate(*kapton, impact);
    EXPECT_LT(result.absorbed_charge_coulomb, 0.0);
    EXPECT_GT(result.deposited_heat_j_per_m2, 0.0);
    EXPECT_GE(result.secondary_emitted_electrons, 0.0);
}

TEST(SurfaceInteractionTest, SecondaryYieldPeaksNearConfiguredEnergy)
{
    SCDAT::Material::MaterialProperty material(9, SCDAT::Mesh::MaterialType::DIELECTRIC, "test");
    material.setSecondaryElectronYield(1.5);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 400.0);

    SurfaceInteraction interaction;
    SurfaceImpactState low_energy;
    low_energy.incident_energy_ev = 50.0;
    low_energy.particle_charge_coulomb = -1.602176634e-19;

    SurfaceImpactState peak_energy = low_energy;
    peak_energy.incident_energy_ev = 400.0;

    SurfaceImpactState high_energy = low_energy;
    high_energy.incident_energy_ev = 1600.0;

    const auto low_result = interaction.evaluate(material, low_energy);
    const auto peak_result = interaction.evaluate(material, peak_energy);
    const auto high_result = interaction.evaluate(material, high_energy);

    EXPECT_GT(peak_result.secondary_emitted_electrons, low_result.secondary_emitted_electrons);
    EXPECT_GT(peak_result.secondary_emitted_electrons, high_result.secondary_emitted_electrons);
}
