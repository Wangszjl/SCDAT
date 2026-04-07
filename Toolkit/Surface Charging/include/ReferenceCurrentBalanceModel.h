#pragma once

#include "../../Plasma Analysis/include/FluidAlgorithmConfig.h"
#include "SurfaceElectronCollectionModel.h"

#include "../../Tools/Material/include/MaterialProperty.h"
#include "../../Tools/Particle/include/ParticleSource.h"

#include <limits>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

enum class SecondaryElectronEmissionModel
{
    Whipple,
    Sims,
    Katz
};

enum class ReferenceSurfaceRole
{
    Body,
    Patch
};

struct ReferenceCurrentComponents
{
    double electron_collection_a_per_m2 = 0.0;
    double ion_collection_a_per_m2 = 0.0;
    double secondary_electron_a_per_m2 = 0.0;
    double ion_secondary_electron_a_per_m2 = 0.0;
    double backscatter_electron_a_per_m2 = 0.0;
    double photoelectron_a_per_m2 = 0.0;
    double conduction_a_per_m2 = 0.0;
    double ram_ion_a_per_m2 = 0.0;
    double net_a_per_m2 = 0.0;
};

struct ReferenceCurrentBalanceConfig
{
    PlasmaAnalysis::PlasmaParameters plasma;
    Particle::ResolvedSpectrum electron_spectrum;
    Particle::ResolvedSpectrum ion_spectrum;
    bool has_electron_spectrum = false;
    bool has_ion_spectrum = false;
    Material::MaterialProperty body_material{1, Mesh::MaterialType::CONDUCTOR, "body"};
    Material::MaterialProperty patch_material{2, Mesh::MaterialType::DIELECTRIC, "patch"};
    SecondaryElectronEmissionModel see_model = SecondaryElectronEmissionModel::Whipple;
    double body_photo_current_density_a_per_m2 = 0.0;
    double patch_photo_current_density_a_per_m2 = 0.0;
    double photoelectron_temperature_ev = 2.0;
    double patch_incidence_angle_deg = 0.0;
    double patch_flow_angle_deg = 0.0;
    double patch_thickness_m = 1.0e-4;
    double patch_conductivity_s_per_m = 0.0;
    ElectronCollectionModelKind electron_collection_model =
        ElectronCollectionModelKind::OmlLike;
    double electron_collection_coefficient = 1.0;
    double ion_collection_coefficient = 1.0;
    double bulk_flow_velocity_m_per_s = 0.0;
    double flow_alignment_cosine = 1.0;
    double ion_directed_velocity_m_per_s = 0.0;
    double electron_calibration_factor = 1.0;
    double ion_calibration_factor = 1.0;
    double ion_secondary_yield = 0.08;
    double backscatter_atomic_number = 13.0;
    bool enable_ram_current = false;
    bool use_photoelectron_suppression = true;
    bool enable_secondary_electron = true;
    bool enable_backscatter = true;
    bool enable_photoelectron = true;
};

class ReferenceCurrentBalanceModel
{
  public:
    bool configure(const ReferenceCurrentBalanceConfig& config);

    ReferenceCurrentComponents evaluate(ReferenceSurfaceRole role, double body_potential_v,
                                        double patch_potential_v) const;

    double computeNetCurrentDensity(ReferenceSurfaceRole role, double body_potential_v,
                                    double patch_potential_v) const;

    double computeCurrentDerivative(ReferenceSurfaceRole role, double body_potential_v,
                                    double patch_potential_v) const;

    double solveEquilibriumPotential(ReferenceSurfaceRole role, double body_potential_v,
                                     double patch_potential_v,
                                     double minimum_potential_v =
                                         -std::numeric_limits<double>::infinity(),
                                     double maximum_potential_v =
                                         std::numeric_limits<double>::infinity()) const;

    const ReferenceCurrentBalanceConfig& getConfig() const
    {
        return config_;
    }

  private:
    static double whippleYield(double peak_energy_ev, double peak_yield,
                               double incident_energy_ev, double surface_potential_v);
    static double simsYield(double peak_energy_ev, double peak_yield, double exponent_n,
                            double incident_energy_ev, double surface_potential_v);
    static double katzYield(const Material::MaterialProperty& material, double incident_energy_ev,
                            double surface_potential_v);
    static double ionSecondaryYield(double peak_energy_kev, double peak_yield,
                                    double incident_energy_ev, double surface_potential_v);
    static double backscatterYield(double incident_energy_ev, double surface_potential_v,
                                   double atomic_number);
    static double bodyPhotoCurrentDensity(double base_photo_current_density_a_per_m2);
    static double patchPhotoCurrentDensity(double base_photo_current_density_a_per_m2,
                                           double incidence_angle_deg);
    static double conductionCurrentDensity(double body_potential_v, double patch_potential_v,
                                           double conductivity_s_per_m, double thickness_m);
    static double ramBodyCurrentDensity(double surface_potential_v, double ion_density_m3,
                                        double ion_temperature_ev, double ion_mass_amu,
                                        double flow_speed_m_per_s);
    static double ramPatchCurrentDensity(double surface_potential_v, double ion_density_m3,
                                         double ion_temperature_ev, double ion_mass_amu,
                                         double projected_flow_speed_m_per_s);
    static double safeExp(double exponent);
    static double erfApprox(double x);

    double computeElectronCollection(double surface_potential_v) const;
    double computeIonCollection(double surface_potential_v) const;
    double computeSeeYield(const Material::MaterialProperty& material, double incident_energy_ev,
                           double surface_potential_v) const;
    double emissionEscapeProbability(double surface_potential_v, double characteristic_energy_ev) const;
    double rolePotential(ReferenceSurfaceRole role, double body_potential_v,
                         double patch_potential_v) const;

    ReferenceCurrentBalanceConfig config_{};
    bool configured_ = false;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
