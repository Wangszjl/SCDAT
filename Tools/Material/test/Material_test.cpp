#include "MaterialDatabase.h"
#include "SurfaceMaterialLoader.h"
#include "SurfaceInteraction.h"
#include "SurfaceMaterialModel.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cmath>

using SCDAT::Material::MaterialDatabase;
using SCDAT::Material::BasicSurfaceMaterialModel;
using SCDAT::Material::ErosionSurfaceMaterialModel;
using SCDAT::Material::SurfaceMaterialModelVariant;
using SCDAT::Material::SurfaceMaterialLoader;
using SCDAT::Material::SurfaceImpactState;
using SCDAT::Material::SurfaceInteraction;
using SCDAT::Material::SurfaceIncidentParticle;
using SCDAT::Material::SurfaceModelContext;

namespace
{

struct SurfaceRolePhotoConductFixture
{
    SCDAT::Material::MaterialProperty material;
    SCDAT::Material::SurfaceRoleCurrentInputs patch_inputs;
    SCDAT::Material::SurfaceRoleCurrentInputs body_inputs;
};

SurfaceRolePhotoConductFixture makeSurfaceRolePhotoConductFixture()
{
    SurfaceRolePhotoConductFixture fixture{
        SCDAT::Material::MaterialProperty(
            23, SCDAT::Mesh::MaterialType::DIELECTRIC, "role_patch_body_split"),
        {},
        {},
    };
    fixture.material.setPhotoelectronYield(2.0);
    fixture.material.setPhotoelectronTemperatureEv(2.0);

    fixture.patch_inputs.role = SCDAT::Material::SurfaceCurrentRole::Patch;
    fixture.patch_inputs.electron_characteristic_energy_ev = 12.0;
    fixture.patch_inputs.ion_characteristic_energy_ev = 4.0;
    fixture.patch_inputs.patch_incidence_angle_deg = 60.0;
    fixture.patch_inputs.patch_flow_angle_deg = 20.0;
    fixture.patch_inputs.surface_potential_v = -5.0;
    fixture.patch_inputs.reference_potential_v = 0.0;
    fixture.patch_inputs.counterelectrode_potential_v = 10.0;
    fixture.patch_inputs.patch_photoelectron_temperature_ev = 1.5;
    fixture.patch_inputs.patch_photo_current_density_a_per_m2 = 3.0e-6;
    fixture.patch_inputs.dielectric_thickness_m = 1.0e-4;
    fixture.patch_inputs.exposed_area_m2 = 1.0;
    fixture.patch_inputs.effective_conductivity_s_per_m = 5.0e-10;
    fixture.patch_inputs.conductivity_scale = 2.0;

    fixture.body_inputs = fixture.patch_inputs;
    fixture.body_inputs.role = SCDAT::Material::SurfaceCurrentRole::Body;
    fixture.body_inputs.body_photoelectron_temperature_ev = 2.5;
    fixture.body_inputs.body_photo_current_density_a_per_m2 = 1.0e-5;
    fixture.body_inputs.body_photo_emission_scale = 1.5;

    return fixture;
}

} // namespace

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
    SurfaceMaterialLoader loader;
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

TEST(MaterialPropertyTest, TypedSpisAccessorsMirrorScalarSurfaceFields)
{
    SCDAT::Material::MaterialProperty material(7, SCDAT::Mesh::MaterialType::DIELECTRIC, "typed");
    material.setSecondaryElectronYieldMax(1.9);
    material.setSecondaryElectronPeakEnergyEv(420.0);
    material.setProtonSecondaryElectronYieldMax(0.12);
    material.setPhotoelectronYield(2.5e-5);
    material.setPhotoelectronTemperatureEv(1.8);
    material.setSurfaceAbsorptionProbability(0.82);
    material.setSurfaceReflectionCoefficient(0.12);
    material.setSurfaceEmissionScaling(1.7);
    material.setSurfaceSecondaryScale(0.9);
    material.setSurfaceIonSecondaryScale(1.1);
    material.setSurfaceBackscatterScale(0.6);
    material.setSurfacePhotoEmissionScale(2.2);
    material.setSurfaceCapacitanceFPerM2(4.0e-7);
    material.setSurfaceCapacitanceScale(2.4);
    material.setSurfaceConductivitySPerM(8.0e-11);
    material.setSurfaceConductivityScale(5.5);
    material.setCoatingThicknessM(6.0e-6);

    EXPECT_NEAR(material.getSecondaryElectronYieldMax(), 1.9, 1.0e-12);
    EXPECT_NEAR(material.getSecondaryElectronPeakEnergyEv(), 420.0, 1.0e-12);
    EXPECT_NEAR(material.getProtonSecondaryElectronYieldMax(), 0.12, 1.0e-12);
    EXPECT_NEAR(material.getPhotoelectronYield(), 2.5e-5, 1.0e-18);
    EXPECT_NEAR(material.getPhotoelectronTemperatureEv(), 1.8, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceAbsorptionProbability(), 0.82, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceReflectionCoefficient(), 0.12, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceEmissionScaling(), 1.7, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceSecondaryScale(), 0.9, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceIonSecondaryScale(), 1.1, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceBackscatterScale(), 0.6, 1.0e-12);
    EXPECT_NEAR(material.getSurfacePhotoEmissionScale(), 2.2, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceCapacitanceFPerM2(), 4.0e-7, 1.0e-18);
    EXPECT_NEAR(material.getSurfaceCapacitanceScale(), 2.4, 1.0e-12);
    EXPECT_NEAR(material.getSurfaceConductivitySPerM(), 8.0e-11, 1.0e-24);
    EXPECT_NEAR(material.getSurfaceConductivityScale(), 5.5, 1.0e-12);
    EXPECT_NEAR(material.getCoatingThicknessM(), 6.0e-6, 1.0e-18);
}

TEST(MaterialPropertyTest, DerivesSurfaceCapacitanceAndConductivityScaleFactors)
{
    SCDAT::Material::MaterialProperty material(12, SCDAT::Mesh::MaterialType::DIELECTRIC, "derived");
    material.setConductivity(2.0e-12);
    material.setDarkConductivitySPerM(5.0e-12);
    material.setSurfaceConductivitySPerM(4.0e-10);
    material.setSurfaceConductivityScale(6.0);
    material.setSurfaceCapacitanceFPerM2(4.5e-7);
    material.setSurfaceCapacitanceScale(3.2);

    EXPECT_NEAR(material.deriveSurfaceCapacitanceScaleFactor(), 3.2, 1.0e-12);
    EXPECT_NEAR(material.deriveSurfaceConductivityScaleFactor(),
                6.0 * 4.0e-10 / 5.0e-12, 1.0e-6);
}

TEST(BasicSurfaceMaterialModelTest, ProducesComponentResolvedCurrents)
{
    SCDAT::Material::MaterialProperty material(8, SCDAT::Mesh::MaterialType::DIELECTRIC, "kapton_like");
    material.setSecondaryElectronYieldMax(1.8);
    material.setSecondaryElectronPeakEnergyEv(300.0);
    material.setProtonSecondaryElectronYieldMax(0.10);
    material.setProtonSecondaryElectronPeakEnergyEv(350.0);
    material.setElectronBackscatterYield(0.2);
    material.setPhotoelectronYield(3.0e-5);
    material.setDarkConductivitySPerM(5.0e-10);

    SurfaceModelContext context;
    context.incident_particle = SurfaceIncidentParticle::Electron;
    context.incident_energy_ev = 300.0;
    context.incident_angle_rad = 0.2;
    context.surface_potential_v = -5.0;
    context.reference_potential_v = 0.0;
    context.source_current_density_a_per_m2 = 2.0e-6;

    BasicSurfaceMaterialModel model;
    const auto currents = model.evaluateCurrents(material, context);
    EXPECT_GT(currents.electron_collection_a_per_m2, 0.0);
    EXPECT_GT(currents.secondary_electron_emission_a_per_m2, 0.0);
    EXPECT_GT(currents.backscatter_emission_a_per_m2, 0.0);
    EXPECT_GT(currents.photoelectron_emission_a_per_m2, 0.0);
    EXPECT_TRUE(std::isfinite(currents.net_current_a_per_m2));
}

TEST(BasicSurfaceMaterialModelTest, UsesContextSourceAndCounterelectrodeForPhotoAndConductiveTerms)
{
    SCDAT::Material::MaterialProperty material(10, SCDAT::Mesh::MaterialType::DIELECTRIC, "photo_cond");
    material.setPhotoelectronYield(2.0);
    material.setPhotoelectronTemperatureEv(2.0);
    material.setDarkConductivitySPerM(4.0e-10);

    SurfaceModelContext context;
    context.incident_particle = SurfaceIncidentParticle::Photon;
    context.source_current_density_a_per_m2 = 3.0e-6;
    context.surface_potential_v = -4.0;
    context.reference_potential_v = 0.0;
    context.counterelectrode_potential_v = 8.0;
    context.transport_length_m = 2.0e-4;

    BasicSurfaceMaterialModel model;
    const auto currents = model.evaluateCurrents(material, context);

    EXPECT_GT(currents.photoelectron_emission_a_per_m2, 0.0);
    EXPECT_GT(currents.induced_conduction_a_per_m2, 0.0);
    EXPECT_TRUE(std::isfinite(currents.net_current_a_per_m2));
}

TEST(ErosionSurfaceMaterialModelTest, ModulatesEmissionAndConductionWithErosionActivation)
{
    SCDAT::Material::MaterialProperty material(31, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "erosion_material");
    material.setSecondaryElectronYieldMax(1.6);
    material.setSecondaryElectronPeakEnergyEv(320.0);
    material.setElectronBackscatterYield(0.2);
    material.setPhotoelectronYield(1.2);
    material.setPhotoelectronTemperatureEv(2.0);
    material.setDarkConductivitySPerM(3.0e-7);
    material.setSurfaceConductivitySPerM(3.0e-7);
    material.setScalarProperty("surface_material_model_use_erosion", 1.0);
    material.setScalarProperty("erosion_yield_scale", 0.8);
    material.setScalarProperty("erosion_threshold_energy_ev", 40.0);
    material.setScalarProperty("erosion_characteristic_energy_ev", 180.0);
    material.setScalarProperty("erosion_secondary_gain", 0.6);
    material.setScalarProperty("erosion_backscatter_gain", 0.5);
    material.setScalarProperty("erosion_photo_suppression", 0.4);
    material.setScalarProperty("erosion_conductivity_suppression", 0.5);

    SurfaceModelContext context;
    context.incident_particle = SurfaceIncidentParticle::Electron;
    context.incident_energy_ev = 700.0;
    context.incident_angle_rad = 0.7;
    context.source_current_density_a_per_m2 = 3.0e-6;
    context.surface_potential_v = -4.0;
    context.reference_potential_v = 0.0;
    context.counterelectrode_potential_v = 9.0;
    context.transport_length_m = 1.0e-4;

    BasicSurfaceMaterialModel basic_model;
    ErosionSurfaceMaterialModel erosion_model;

    const auto baseline = basic_model.evaluateCurrents(material, context);
    const auto erosion = erosion_model.evaluateCurrents(material, context);

    EXPECT_STREQ(erosion_model.modelFamily(), "spis_erosion_surface_material_model_v1");
    EXPECT_GT(erosion.secondary_electron_emission_a_per_m2,
              baseline.secondary_electron_emission_a_per_m2);
    EXPECT_GT(erosion.backscatter_emission_a_per_m2,
              baseline.backscatter_emission_a_per_m2);
    EXPECT_LT(erosion.photoelectron_emission_a_per_m2,
              baseline.photoelectron_emission_a_per_m2);
    EXPECT_LT(erosion.induced_conduction_a_per_m2,
              baseline.induced_conduction_a_per_m2);
    EXPECT_TRUE(std::isfinite(erosion.net_current_a_per_m2));
}

TEST(SurfaceMaterialModelVariantTest, ResolvesModelFamilyFromMaterialProperties)
{
    SCDAT::Material::MaterialProperty material(32, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "variant_selector");

    EXPECT_EQ(SCDAT::Material::resolveSurfaceMaterialModelVariant(material),
              SurfaceMaterialModelVariant::Basic);
    EXPECT_STREQ(SCDAT::Material::resolveSurfaceMaterialModelFamily(material),
                 "spis_basic_surface_material_model_v1");

    material.setScalarProperty("erosion_yield_scale", 0.2);
    EXPECT_EQ(SCDAT::Material::resolveSurfaceMaterialModelVariant(material),
              SurfaceMaterialModelVariant::Erosion);
    EXPECT_STREQ(SCDAT::Material::resolveSurfaceMaterialModelFamily(material),
                 "spis_erosion_surface_material_model_v1");

    material.setScalarProperty("surface_material_model_use_erosion", 0.0);
    EXPECT_EQ(SCDAT::Material::resolveSurfaceMaterialModelVariant(material),
              SurfaceMaterialModelVariant::Basic);

    material.setScalarProperty("surface_material_model_use_erosion", 1.0);
    EXPECT_EQ(SCDAT::Material::resolveSurfaceMaterialModelVariant(material),
              SurfaceMaterialModelVariant::Erosion);
}

TEST(MaterialEmissionHelperTest, FieldEmissionUsesFieldThreshold)
{
    SCDAT::Material::MaterialProperty material(
        13, SCDAT::Mesh::MaterialType::CONDUCTOR, "field_emission");
    material.setWorkFunctionEv(4.5);

    const double below_threshold = SCDAT::Material::evaluateSurfaceFieldEmissionCurrentDensity(
        material, -10.0, 0.0, 5.0e6, 1.0e-4);
    const double above_threshold = SCDAT::Material::evaluateSurfaceFieldEmissionCurrentDensity(
        material, -10.0, 0.0, 2.0e7, 1.0e-4);

    EXPECT_DOUBLE_EQ(below_threshold, 0.0);
    EXPECT_GT(above_threshold, 0.0);
}

TEST(MaterialEmissionHelperTest, ThermionicEmissionIncreasesWithTemperature)
{
    SCDAT::Material::MaterialProperty material(
        14, SCDAT::Mesh::MaterialType::CONDUCTOR, "thermionic");
    material.setWorkFunctionEv(4.2);

    const double low_temperature =
        SCDAT::Material::evaluateSurfaceThermionicEmissionCurrentDensity(material, 300.0, -2.0);
    const double high_temperature =
        SCDAT::Material::evaluateSurfaceThermionicEmissionCurrentDensity(material, 1200.0, -2.0);

    EXPECT_GE(low_temperature, 0.0);
    EXPECT_GT(high_temperature, low_temperature);
}

TEST(MaterialEmissionHelperTest, EffectiveConductivityTracksFieldAndSpaceChargeEnhancement)
{
    SCDAT::Material::MaterialProperty material(
        15, SCDAT::Mesh::MaterialType::DIELECTRIC, "conductivity_helper");
    material.setDarkConductivitySPerM(2.0e-12);
    material.setSurfaceConductivitySPerM(5.0e-10);
    material.setSurfaceConductivityScale(3.0);
    material.setScalarProperty("poole_frenkel_beta", 2.0e-4);
    material.setScalarProperty("max_field_enhancement_factor", 50.0);

    const double weak_state = SCDAT::Material::evaluateSurfaceEffectiveConductivitySPerM(
        material, 2.0e-12, 0.0, -5.0, 0.0, 5.0e3, 0.0, 1.0e-4);
    const double strong_state = SCDAT::Material::evaluateSurfaceEffectiveConductivitySPerM(
        material, 2.0e-12, 1.0e-12, -50.0, 0.0, 5.0e6, 3.0e-6, 1.0e-4);

    EXPECT_GT(weak_state, 0.0);
    EXPECT_GT(strong_state, weak_state);
}

TEST(BasicSurfaceMaterialModelTest, RoleCurrentBundleBuildsElectronIonPhotoAndConductionTerms)
{
    SCDAT::Material::MaterialProperty material(18, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "role_bundle");
    material.setPhotoelectronYield(1.5);
    material.setPhotoelectronTemperatureEv(1.9);
    material.setDarkConductivitySPerM(3.0e-10);

    SCDAT::Material::SurfaceRoleCurrentInputs inputs;
    inputs.role = SCDAT::Material::SurfaceCurrentRole::Patch;
    inputs.electron_characteristic_energy_ev = 14.0;
    inputs.ion_characteristic_energy_ev = 5.0;
    inputs.patch_incidence_angle_deg = 45.0;
    inputs.patch_flow_angle_deg = 20.0;
    inputs.surface_potential_v = -4.0;
    inputs.reference_potential_v = 0.0;
    inputs.counterelectrode_potential_v = 8.0;
    inputs.normal_electric_field_v_per_m = 1.5e4;
    inputs.patch_photoelectron_temperature_ev = 1.7;
    inputs.patch_photo_current_density_a_per_m2 = 3.0e-6;
    inputs.dielectric_thickness_m = 2.0e-4;
    inputs.exposed_area_m2 = 0.5;
    inputs.effective_conductivity_s_per_m = 4.0e-10;
    inputs.conductivity_scale = 2.0;

    const auto bundle = SCDAT::Material::evaluateSurfaceRoleCurrentBundle(material, inputs);
    EXPECT_GT(bundle.electron_currents.electron_collection_a_per_m2, 0.0);
    EXPECT_GT(bundle.ion_currents.ion_collection_a_per_m2, 0.0);
    EXPECT_GT(bundle.photo_currents.photoelectron_emission_a_per_m2, 0.0);
    EXPECT_GT(bundle.photo_currents.induced_conduction_a_per_m2, 0.0);
    EXPECT_TRUE(std::isfinite(bundle.photo_currents.net_current_a_per_m2));
}

TEST(BasicSurfaceMaterialModelTest, RoleCurrentBundleSeparatesPatchAndBodyPhotoConductiveResponses)
{
    const auto fixture = makeSurfaceRolePhotoConductFixture();

    const auto patch_bundle =
        SCDAT::Material::evaluateSurfaceRoleCurrentBundle(fixture.material, fixture.patch_inputs);
    const auto body_bundle =
        SCDAT::Material::evaluateSurfaceRoleCurrentBundle(fixture.material, fixture.body_inputs);

    EXPECT_GT(patch_bundle.photo_currents.photoelectron_emission_a_per_m2, 0.0);
    EXPECT_GT(body_bundle.photo_currents.photoelectron_emission_a_per_m2,
              patch_bundle.photo_currents.photoelectron_emission_a_per_m2);
    EXPECT_GT(patch_bundle.photo_currents.induced_conduction_a_per_m2, 0.0);
    EXPECT_DOUBLE_EQ(body_bundle.photo_currents.induced_conduction_a_per_m2, 0.0);
}

TEST(BasicSurfaceMaterialModelTest, RolePhotoConductInputFixtureBuildsPatchAndBodyVariants)
{
    const auto fixture = makeSurfaceRolePhotoConductFixture();

    EXPECT_EQ(fixture.patch_inputs.role, SCDAT::Material::SurfaceCurrentRole::Patch);
    EXPECT_EQ(fixture.body_inputs.role, SCDAT::Material::SurfaceCurrentRole::Body);
    EXPECT_GT(fixture.patch_inputs.patch_photo_current_density_a_per_m2, 0.0);
    EXPECT_GT(fixture.body_inputs.body_photo_current_density_a_per_m2,
              fixture.patch_inputs.patch_photo_current_density_a_per_m2);
    EXPECT_GT(fixture.patch_inputs.dielectric_thickness_m, 0.0);
    EXPECT_DOUBLE_EQ(fixture.body_inputs.dielectric_thickness_m,
                     fixture.patch_inputs.dielectric_thickness_m);
    EXPECT_GT(fixture.patch_inputs.effective_conductivity_s_per_m, 0.0);
}

TEST(BasicSurfaceMaterialModelTest, RoleSwitchKeepsCollectionConsistencyAndExpectedPhotoConductDelta)
{
    const auto fixture = makeSurfaceRolePhotoConductFixture();
    const auto patch_bundle =
        SCDAT::Material::evaluateSurfaceRoleCurrentBundle(fixture.material, fixture.patch_inputs);
    const auto body_bundle =
        SCDAT::Material::evaluateSurfaceRoleCurrentBundle(fixture.material, fixture.body_inputs);

    const double electron_abs_delta = std::abs(
        patch_bundle.electron_currents.electron_collection_a_per_m2 -
        body_bundle.electron_currents.electron_collection_a_per_m2);
    const double ion_abs_delta = std::abs(
        patch_bundle.ion_currents.ion_collection_a_per_m2 -
        body_bundle.ion_currents.ion_collection_a_per_m2);

    // Body/Patch role switching should preserve collection terms while photo and
    // conduction terms follow role-specific semantics.
    EXPECT_LT(electron_abs_delta, 1.0e-18);
    EXPECT_LT(ion_abs_delta, 1.0e-18);

    const double photo_delta =
        body_bundle.photo_currents.photoelectron_emission_a_per_m2 -
        patch_bundle.photo_currents.photoelectron_emission_a_per_m2;
    EXPECT_GT(photo_delta, 0.0);
    EXPECT_GT(body_bundle.photo_currents.photoelectron_emission_a_per_m2,
              1.1 * patch_bundle.photo_currents.photoelectron_emission_a_per_m2);

    EXPECT_GT(patch_bundle.photo_currents.induced_conduction_a_per_m2, 0.0);
    EXPECT_DOUBLE_EQ(body_bundle.photo_currents.induced_conduction_a_per_m2, 0.0);
}

TEST(BasicSurfaceMaterialModelTest, LegacyYieldHelpersProduceFinitePositiveResponses)
{
    SCDAT::Material::MaterialProperty material(19, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "legacy_yields");
    material.setSecondaryElectronYield(1.6);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 350.0);
    material.setScalarProperty("sims_exponent_n", 1.7);
    material.setScalarProperty("ion_secondary_peak_energy_kev", 0.4);
    material.setScalarProperty("ion_secondary_yield", 0.09);
    material.setScalarProperty("atomic_number", 13.0);

    const double whipple = SCDAT::Material::evaluateLegacySecondaryElectronYield(
        material, SCDAT::Material::LegacySecondaryYieldModel::Whipple, 350.0);
    const double sims = SCDAT::Material::evaluateLegacySecondaryElectronYield(
        material, SCDAT::Material::LegacySecondaryYieldModel::Sims, 350.0, 1.7);
    const double katz = SCDAT::Material::evaluateLegacySecondaryElectronYield(
        material, SCDAT::Material::LegacySecondaryYieldModel::Katz, 350.0);
    const double ion_secondary =
        SCDAT::Material::evaluateLegacyIonSecondaryElectronYield(material, 800.0);
    const double backscatter =
        SCDAT::Material::evaluateLegacyBackscatterYield(material, 1200.0);

    EXPECT_GT(whipple, 0.0);
    EXPECT_GT(sims, 0.0);
    EXPECT_GT(katz, 0.0);
    EXPECT_GT(ion_secondary, 0.0);
    EXPECT_GT(backscatter, 0.0);
}

TEST(BasicSurfaceMaterialModelTest, LegacyRamCurrentHelpersMatchReferenceFormulas)
{
    const double surface_potential_v = 5.0;
    const double ion_density_m3 = 8.5e10;
    const double ion_temperature_ev = 0.22;
    const double ion_mass_amu = 16.0;
    const double projected_flow_speed_m_per_s = 7800.0;

    const auto expected_patch = [&]() {
        const double ni_cm3 = std::max(0.0, ion_density_m3 * 1.0e-6);
        const double ti_ev = std::max(1.0e-6, ion_temperature_ev);
        const double wt = 13.84 * std::sqrt(ti_ev / std::max(1.0, ion_mass_amu));
        const double qd = projected_flow_speed_m_per_s * 1.0e-3 / std::max(1.0e-12, wt);
        const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
        double current_na_per_m2 = 0.08011 * ni_cm3 * wt;
        current_na_per_m2 *=
            0.5642 * std::exp(-(qd - qv) * (qd - qv)) + qd + qd * std::erf(qd - qv);
        return std::max(0.0, current_na_per_m2) * 1.0e-9;
    }();

    const double actual_patch = SCDAT::Material::evaluateLegacyRamPatchCurrentDensity(
        surface_potential_v, ion_density_m3, ion_temperature_ev, ion_mass_amu,
        projected_flow_speed_m_per_s);
    EXPECT_NEAR(actual_patch, expected_patch,
                1.0e-12 * std::max(1.0, std::abs(expected_patch)));

    const auto expected_body = [&](double body_surface_potential_v) {
        const double ni_cm3 = std::max(0.0, ion_density_m3 * 1.0e-6);
        const double ti_ev = std::max(1.0e-6, ion_temperature_ev);
        const double wt = 13.84 * std::sqrt(ti_ev / std::max(1.0, ion_mass_amu));
        const double qd = std::max(1.0e-12,
                                   projected_flow_speed_m_per_s * 1.0e-3 /
                                       std::max(1.0e-12, wt));
        const double qv = std::sqrt(std::abs(body_surface_potential_v) / ti_ev);
        double current_na_per_m2 = 0.0;
        if (body_surface_potential_v >= 0.0)
        {
            current_na_per_m2 =
                0.5642 / qd *
                    ((qv + qd) * std::exp(-(qv - qd) * (qv - qd)) -
                     (qv - qd) * std::exp(-(qv + qd) * (qv + qd))) +
                (0.5 / qd + qd - qv * qv / qd) *
                    (std::erf(qv + qd) - std::erf(qv - qd));
            current_na_per_m2 *= 0.02003 * ni_cm3 * wt;
        }
        else
        {
            current_na_per_m2 = 0.04005 * ni_cm3 * wt *
                                (0.5642 * std::exp(-qd * qd) +
                                 (qd + 0.5 / qd) * std::erf(qd));
            current_na_per_m2 *= 4.0;
        }
        return std::max(0.0, current_na_per_m2) * 1.0e-9;
    };

    const double expected_body_positive = expected_body(5.0);
    const double expected_body_negative = expected_body(-3.0);
    const double actual_body_positive = SCDAT::Material::evaluateLegacyRamBodyCurrentDensity(
        5.0, ion_density_m3, ion_temperature_ev, ion_mass_amu,
        projected_flow_speed_m_per_s);
    const double actual_body_negative = SCDAT::Material::evaluateLegacyRamBodyCurrentDensity(
        -3.0, ion_density_m3, ion_temperature_ev, ion_mass_amu,
        projected_flow_speed_m_per_s);

    EXPECT_NEAR(actual_body_positive, expected_body_positive,
                1.0e-12 * std::max(1.0, std::abs(expected_body_positive)));
    EXPECT_NEAR(actual_body_negative, expected_body_negative,
                1.0e-12 * std::max(1.0, std::abs(expected_body_negative)));
    EXPECT_GE(actual_body_positive, 0.0);
    EXPECT_GE(actual_body_negative, 0.0);
}

TEST(BasicSurfaceMaterialModelTest, MaxwellianYieldIntegratorReturnsWeightedAverage)
{
    const double constant_yield =
        SCDAT::Material::integrateMaxwellianYield(12.0, [](double) { return 0.75; });
    const double ramp_yield = SCDAT::Material::integrateMaxwellianYield(
        12.0, [](double energy_ev) { return energy_ev > 8.0 ? 1.0 : 0.25; });

    EXPECT_NEAR(constant_yield, 0.75, 1.0e-6);
    EXPECT_GT(ramp_yield, 0.25);
    EXPECT_LT(ramp_yield, 1.0);
}

TEST(BasicSurfaceMaterialModelTest, LegacyEmissionIntegralBundleResolvesAllComponents)
{
    SCDAT::Material::MaterialProperty material(20, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "legacy_bundle");
    material.setSecondaryElectronYield(1.5);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 300.0);
    material.setScalarProperty("sims_exponent_n", 1.6);
    material.setScalarProperty("ion_secondary_peak_energy_kev", 0.35);
    material.setScalarProperty("ion_secondary_yield", 0.08);
    material.setScalarProperty("atomic_number", 13.0);
    material.setScalarProperty("secondary_emission_escape_energy_ev", 2.0);

    const auto bundle = SCDAT::Material::evaluateLegacyEmissionIntegralBundle(
        material, SCDAT::Material::LegacySecondaryYieldModel::Sims, 12.0, 4.0, 3.0, 1.6);
    EXPECT_GT(bundle.secondary_integral, 0.0);
    EXPECT_GT(bundle.backscatter_integral, 0.0);
    EXPECT_GT(bundle.ion_secondary_integral, 0.0);
    EXPECT_GT(bundle.emission_escape_probability, 0.0);
    EXPECT_LT(bundle.emission_escape_probability, 1.0);
}

TEST(BasicSurfaceMaterialModelTest, LegacyEmissionResponseScalesWithCollectionCurrents)
{
    SCDAT::Material::MaterialProperty material(21, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "legacy_response");
    material.setSecondaryElectronYield(1.4);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 280.0);
    material.setScalarProperty("sims_exponent_n", 1.6);
    material.setScalarProperty("ion_secondary_peak_energy_kev", 0.35);
    material.setScalarProperty("ion_secondary_yield", 0.08);
    material.setScalarProperty("atomic_number", 13.0);
    material.setScalarProperty("secondary_emission_escape_energy_ev", 2.0);

    SCDAT::Material::LegacyEmissionResponseInputs inputs;
    inputs.model = SCDAT::Material::LegacySecondaryYieldModel::Whipple;
    inputs.electron_temperature_ev = 10.0;
    inputs.ion_temperature_ev = 4.0;
    inputs.surface_potential_v = 2.0;
    inputs.electron_collection_current_a_per_m2 = -3.0e-6;
    inputs.ion_collection_current_a_per_m2 = 1.5e-6;
    inputs.sims_exponent_n = 1.6;

    const auto response = SCDAT::Material::evaluateLegacyEmissionResponse(material, inputs);
    EXPECT_GT(response.secondary_emission_a_per_m2, 0.0);
    EXPECT_GT(response.backscatter_emission_a_per_m2, 0.0);
    EXPECT_GT(response.ion_secondary_emission_a_per_m2, 0.0);
    EXPECT_GT(response.integrals.emission_escape_probability, 0.0);
    EXPECT_LT(response.integrals.emission_escape_probability, 1.0);
}

TEST(BasicSurfaceMaterialModelTest, LegacyPopulationResponseBuildsCollectionAndEmission)
{
    SCDAT::Material::MaterialProperty material(22, SCDAT::Mesh::MaterialType::DIELECTRIC,
                                               "legacy_population");
    material.setSecondaryElectronYield(1.3);
    material.setScalarProperty("secondary_yield_peak_energy_ev", 250.0);
    material.setScalarProperty("sims_exponent_n", 1.6);
    material.setScalarProperty("ion_secondary_peak_energy_kev", 0.3);
    material.setScalarProperty("ion_secondary_yield", 0.07);
    material.setScalarProperty("atomic_number", 13.0);

    const auto electron_response = SCDAT::Material::evaluateLegacyPopulationResponse(
        material,
        SCDAT::Material::LegacyPopulationResponseInputs{
            SCDAT::Material::LegacyCollectionParticle::Electron,
            SCDAT::Material::LegacySecondaryYieldModel::Whipple,
            1.0e11,
            8.0,
            5.48579909065e-4,
            -2.0,
            8.0,
            4.0,
            1.1,
            0.8,
            1.6,
        });
    EXPECT_LT(electron_response.collection_current_a_per_m2, 0.0);
    EXPECT_GT(electron_response.emission.secondary_emission_a_per_m2, 0.0);
    EXPECT_GT(electron_response.emission.backscatter_emission_a_per_m2, 0.0);

    const auto ion_response = SCDAT::Material::evaluateLegacyPopulationResponse(
        material,
        SCDAT::Material::LegacyPopulationResponseInputs{
            SCDAT::Material::LegacyCollectionParticle::Ion,
            SCDAT::Material::LegacySecondaryYieldModel::Whipple,
            8.0e10,
            4.0,
            16.0,
            2.0,
            8.0,
            4.0,
            1.1,
            0.8,
            1.6,
        });
    EXPECT_GT(ion_response.collection_current_a_per_m2, 0.0);
    EXPECT_GT(ion_response.emission.ion_secondary_emission_a_per_m2, 0.0);
}

TEST(SurfaceInteractionTest, ExportsMaterialModelFamilyAndComponentBundle)
{
    SCDAT::Material::MaterialProperty material(9, SCDAT::Mesh::MaterialType::DIELECTRIC, "test");
    material.setSecondaryElectronYieldMax(1.4);
    material.setPhotoelectronYield(2.0e-5);

    SurfaceInteraction interaction;
    SurfaceImpactState impact;
    impact.incident_energy_ev = 200.0;
    impact.particle_charge_coulomb = -1.602176634e-19;
    impact.incident_angle_rad = 0.4;

    const auto result = interaction.evaluate(material, impact);
    EXPECT_EQ(result.material_model_family, "spis_basic_surface_material_model_v1");
    EXPECT_TRUE(std::isfinite(result.component_currents.net_current_a_per_m2));
    EXPECT_GE(result.backscattered_electrons, 0.0);
}

TEST(SurfaceInteractionTest, EmittedCountsTrackComponentCurrents)
{
    SCDAT::Material::MaterialProperty material(11, SCDAT::Mesh::MaterialType::DIELECTRIC, "counts");
    material.setSecondaryElectronYieldMax(1.7);
    material.setSecondaryElectronPeakEnergyEv(250.0);
    material.setElectronBackscatterYield(0.25);

    SurfaceInteraction interaction;
    interaction.setEmissionScaling(0.75);

    SurfaceImpactState impact;
    impact.incident_energy_ev = 250.0;
    impact.particle_charge_coulomb = -1.602176634e-19;
    impact.incident_angle_rad = 0.3;

    const auto result = interaction.evaluate(material, impact);
    const double incident_collection = result.component_currents.electron_collection_a_per_m2;

    ASSERT_GT(incident_collection, 0.0);
    EXPECT_NEAR(result.secondary_emitted_electrons,
                result.component_currents.secondary_electron_emission_a_per_m2 / incident_collection,
                1.0e-12);
    EXPECT_NEAR(result.backscattered_electrons,
                result.component_currents.backscatter_emission_a_per_m2 / incident_collection,
                1.0e-12);
}

TEST(SurfaceInteractionTest, BundleEvaluationProducesBarrierResolvedOutputs)
{
    SCDAT::Material::MaterialProperty material(13, SCDAT::Mesh::MaterialType::DIELECTRIC, "bundle");
    material.setSecondaryElectronYieldMax(1.6);
    material.setProtonSecondaryElectronYieldMax(0.15);
    material.setElectronBackscatterYield(0.2);
    material.setPhotoelectronYield(1.0e-5);
    material.setPhotoelectronTemperatureEv(2.0);
    material.setSurfaceSecondaryScale(1.1);
    material.setSurfaceIonSecondaryScale(0.9);
    material.setSurfaceBackscatterScale(0.8);
    material.setSurfacePhotoEmissionScale(1.2);

    SurfaceInteraction interaction;

    SCDAT::Material::SurfaceInteractionBundleRequest request;
    request.electron_impact.incident_energy_ev = 300.0;
    request.electron_impact.incident_angle_rad = 0.2;
    request.electron_impact.particle_charge_coulomb = -1.602176634e-19;
    request.electron_impact.surface_potential_v = -3.0;
    request.electron_impact.reference_potential_v = 0.0;
    request.electron_impact.normal_electric_field_v_per_m = 2.0e4;
    request.electron_impact.emission_temperature_ev = 2.0;
    request.electron_impact.source_current_density_a_per_m2 = 1.0e-6;
    request.electron_impact.counterelectrode_potential_v = 2.0;
    request.electron_impact.transport_length_m = 1.0e-4;

    request.ion_impact = request.electron_impact;
    request.ion_impact.incident_energy_ev = 40.0;
    request.ion_impact.incident_angle_rad = 0.1;
    request.ion_impact.particle_charge_coulomb = 1.602176634e-19;

    request.barrier_state.local_potential_v = -3.0;
    request.barrier_state.reference_potential_v = 0.0;
    request.barrier_state.normal_electric_field_v_per_m = 2.0e4;
    request.barrier_state.emission_temperature_ev = 2.0;
    request.photo_incidence_scale = 0.75;

    const auto result = interaction.evaluateBundle(material, request);
    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.barrier_outputs.valid);
    EXPECT_GT(result.barrier_outputs.secondary_emission_scale, 0.0);
    EXPECT_GT(result.barrier_outputs.photo_emission_scale, 0.0);
}

TEST(SurfaceInteractionTest, EnvironmentBuilderProducesBundleRequestAndValidOutputs)
{
    SCDAT::Material::MaterialProperty material(14, SCDAT::Mesh::MaterialType::DIELECTRIC, "env");
    material.setSecondaryElectronYieldMax(1.4);
    material.setProtonSecondaryElectronYieldMax(0.12);
    material.setElectronBackscatterYield(0.18);
    material.setPhotoelectronYield(8.0e-6);
    material.setPhotoelectronTemperatureEv(1.8);

    SurfaceInteraction interaction;
    SCDAT::Material::SurfaceInteractionEnvironment environment;
    environment.electron_incident_energy_ev = 250.0;
    environment.electron_incident_angle_rad = 0.2;
    environment.ion_incident_energy_ev = 35.0;
    environment.ion_incident_angle_rad = 0.1;
    environment.surface_potential_v = -2.5;
    environment.reference_potential_v = 0.0;
    environment.normal_electric_field_v_per_m = 1.5e4;
    environment.emission_temperature_ev = 1.8;
    environment.source_current_density_a_per_m2 = 8.0e-7;
    environment.counterelectrode_potential_v = 3.0;
    environment.transport_length_m = 1.0e-4;
    environment.photo_incidence_scale = 0.7;

    const auto request = interaction.makeBundleRequest(environment);
    EXPECT_LT(request.electron_impact.particle_charge_coulomb, 0.0);
    EXPECT_GT(request.ion_impact.particle_charge_coulomb, 0.0);
    EXPECT_NEAR(request.barrier_state.local_potential_v, environment.surface_potential_v, 1.0e-12);

    const auto result = interaction.evaluateBundle(material, environment);
    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.barrier_outputs.valid);
}

TEST(SurfaceInteractionTest, MaterialIndexedBundleSelectsDielectricAndConductorFamilies)
{
    SCDAT::Material::MaterialProperty dielectric(
        15, SCDAT::Mesh::MaterialType::DIELECTRIC, "kapton_like");
    dielectric.setSurfaceAbsorptionProbability(0.6);
    dielectric.setSurfaceReflectionCoefficient(0.08);
    dielectric.setSurfaceEmissionScaling(1.1);
    dielectric.setPhotoelectronYield(1.0e-5);
    dielectric.setPhotoelectronTemperatureEv(2.0);

    SCDAT::Material::MaterialProperty conductor(
        16, SCDAT::Mesh::MaterialType::CONDUCTOR, "al_like");
    conductor.setSurfaceAbsorptionProbability(0.6);
    conductor.setSurfaceReflectionCoefficient(0.02);
    conductor.setSurfaceEmissionScaling(0.9);
    conductor.setPhotoelectronYield(1.0e-5);
    conductor.setPhotoelectronTemperatureEv(2.0);

    SurfaceInteraction interaction;
    SCDAT::Material::SurfaceInteractionEnvironment environment;
    environment.electron_incident_energy_ev = 200.0;
    environment.electron_incident_angle_rad = 0.2;
    environment.ion_incident_energy_ev = 30.0;
    environment.ion_incident_angle_rad = 0.1;
    environment.surface_potential_v = -2.0;
    environment.reference_potential_v = 0.0;
    environment.normal_electric_field_v_per_m = 1.0e4;
    environment.emission_temperature_ev = 2.0;
    environment.source_current_density_a_per_m2 = 1.0e-6;
    environment.counterelectrode_potential_v = 2.0;
    environment.transport_length_m = 1.0e-4;

    const auto dielectric_result = interaction.evaluateMaterialIndexedBundle(dielectric, environment);
    const auto conductor_result = interaction.evaluateMaterialIndexedBundle(conductor, environment);

    EXPECT_EQ(dielectric_result.electron_response.material_model_family,
              "spis_material_indexed_dielectric_interaction_v1");
    EXPECT_EQ(conductor_result.electron_response.material_model_family,
              "spis_material_indexed_conductor_interaction_v1");
}

TEST(SurfaceInteractionTest, MaterialIndexedBundleCouplesYieldBackscatterAndBarrierScales)
{
    SCDAT::Material::MaterialProperty dielectric(
        17, SCDAT::Mesh::MaterialType::DIELECTRIC, "indexed_barrier");
    dielectric.setSurfaceAbsorptionProbability(0.78);
    dielectric.setSurfaceReflectionCoefficient(0.10);
    dielectric.setSurfaceEmissionScaling(1.35);
    dielectric.setSurfaceSecondaryScale(1.4);
    dielectric.setSurfaceIonSecondaryScale(1.2);
    dielectric.setSurfaceBackscatterScale(0.9);
    dielectric.setPhotoelectronYield(1.2e-5);
    dielectric.setPhotoelectronTemperatureEv(2.1);
    dielectric.setSecondaryElectronYieldMax(1.7);
    dielectric.setSecondaryElectronPeakEnergyEv(260.0);
    dielectric.setElectronBackscatterYield(0.22);

    SurfaceInteraction interaction;
    SCDAT::Material::SurfaceInteractionEnvironment low_barrier_environment;
    low_barrier_environment.electron_incident_energy_ev = 260.0;
    low_barrier_environment.electron_incident_angle_rad = 0.25;
    low_barrier_environment.ion_incident_energy_ev = 45.0;
    low_barrier_environment.ion_incident_angle_rad = 0.1;
    low_barrier_environment.surface_potential_v = -2.0;
    low_barrier_environment.reference_potential_v = 0.0;
    low_barrier_environment.barrier_potential_v = 0.0;
    low_barrier_environment.normal_electric_field_v_per_m = 1.2e4;
    low_barrier_environment.emission_temperature_ev = 2.1;
    low_barrier_environment.source_current_density_a_per_m2 = 1.5e-6;
    low_barrier_environment.counterelectrode_potential_v = 2.0;
    low_barrier_environment.transport_length_m = 1.0e-4;
    low_barrier_environment.photo_incidence_scale = 0.85;

    auto high_barrier_environment = low_barrier_environment;
    high_barrier_environment.barrier_potential_v = 22.0;

    const auto low_barrier_result =
        interaction.evaluateMaterialIndexedBundle(dielectric, low_barrier_environment);
    const auto high_barrier_result =
        interaction.evaluateMaterialIndexedBundle(dielectric, high_barrier_environment);

    EXPECT_TRUE(low_barrier_result.valid);
    EXPECT_TRUE(high_barrier_result.valid);
    EXPECT_TRUE(low_barrier_result.barrier_outputs.valid);
    EXPECT_TRUE(high_barrier_result.barrier_outputs.valid);

    EXPECT_GT(low_barrier_result.electron_response.component_currents
                  .secondary_electron_emission_a_per_m2,
              0.0);
    EXPECT_GT(low_barrier_result.electron_response.component_currents
                  .backscatter_emission_a_per_m2,
              0.0);

    EXPECT_LE(high_barrier_result.barrier_outputs.escaped_secondary.escaped_current_a_per_m2,
              low_barrier_result.barrier_outputs.escaped_secondary.escaped_current_a_per_m2);
    EXPECT_LE(high_barrier_result.barrier_outputs.escaped_backscatter.escaped_current_a_per_m2,
              low_barrier_result.barrier_outputs.escaped_backscatter.escaped_current_a_per_m2);
}

TEST(SurfaceInteractionTest, RoleEnvironmentBuilderResolvesBodyAndPatchSemantics)
{
    SCDAT::Material::SurfaceInteractionRoleInputs patch_inputs;
    patch_inputs.role = SCDAT::Material::SurfaceInteractionRole::Patch;
    patch_inputs.electron_incident_energy_ev = 18.0;
    patch_inputs.electron_incident_angle_deg = 60.0;
    patch_inputs.ion_incident_energy_ev = 7.0;
    patch_inputs.ion_incident_angle_deg = 25.0;
    patch_inputs.surface_potential_v = -5.0;
    patch_inputs.reference_potential_v = 0.0;
    patch_inputs.counterelectrode_potential_v = 8.0;
    patch_inputs.normal_electric_field_v_per_m = 2.0e4;
    patch_inputs.patch_photoelectron_temperature_ev = 1.7;
    patch_inputs.patch_photo_current_density_a_per_m2 = 4.0e-6;
    patch_inputs.dielectric_thickness_m = 2.5e-4;
    patch_inputs.exposed_area_m2 = 0.25;
    patch_inputs.photo_incidence_scale = 0.8;

    const auto patch_environment =
        SurfaceInteraction::makeInteractionEnvironment(patch_inputs);
    EXPECT_NEAR(patch_environment.electron_incident_angle_rad, 3.14159265358979323846 / 3.0,
                1.0e-12);
    EXPECT_NEAR(patch_environment.source_current_density_a_per_m2, 2.0e-6, 1.0e-12);
    EXPECT_NEAR(patch_environment.transport_length_m, 2.5e-4, 1.0e-12);
    EXPECT_NEAR(patch_environment.barrier_potential_v, 5.0, 1.0e-12);
    EXPECT_NEAR(patch_environment.emission_temperature_ev, 1.7, 1.0e-12);

    SCDAT::Material::SurfaceInteractionRoleInputs body_inputs = patch_inputs;
    body_inputs.role = SCDAT::Material::SurfaceInteractionRole::Body;
    body_inputs.body_photoelectron_temperature_ev = 2.3;
    body_inputs.body_photo_current_density_a_per_m2 = 6.0e-6;
    body_inputs.body_photo_emission_scale = 1.4;

    const auto body_environment =
        SurfaceInteraction::makeInteractionEnvironment(body_inputs);
    EXPECT_NEAR(body_environment.source_current_density_a_per_m2, 0.25 * 1.4 * 6.0e-6, 1.0e-12);
    EXPECT_DOUBLE_EQ(body_environment.transport_length_m, 0.0);
    EXPECT_NEAR(body_environment.emission_temperature_ev, 2.3, 1.0e-12);
    EXPECT_NEAR(body_environment.barrier_potential_v, 5.0, 1.0e-12);
}
