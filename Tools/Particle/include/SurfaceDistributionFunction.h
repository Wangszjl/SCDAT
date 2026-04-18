#pragma once

#include "ParticleSource.h"

#include <cstddef>
#include <functional>
#include <random>
#include <vector>

namespace SCDAT
{
namespace Particle
{

struct SurfaceFluxMoments
{
    double total_flux = 0.0;
    double mean_energy_ev = 0.0;
    double mean_cos_incidence = 1.0;
    bool valid = false;
};

enum class SurfaceDistributionModel
{
    MaxwellianProjected,
    WakeAnisotropic,
    MultiPopulationHybrid,
    TabulatedVelocity,
    NonLocalizedHybrid,
    MultipleSurf,
    LocalModifiedPearsonIV,
    LocalTabulated,
    TwoAxesTabulatedVelocity,
    UniformVelocity,
    Fluid,
    FowlerNordheim,
    AxisymTabulatedVelocity,
    GlobalMaxwellBoltzmann,
    GlobalMaxwellBoltzmann2,
    GlobalMaxwell,
    LocalMaxwell,
    LocalMaxwell2,
    RecollMaxwell,
    PICSurf,
    NonPICSurf,
    GenericSurf,
    GlobalSurf,
    LocalGenericSurf,
    TestableSurf,
    TestableForA,
    MaxwellianThruster
};

struct SurfaceDistributionSynthesisRequest
{
    SurfaceDistributionModel model = SurfaceDistributionModel::MultiPopulationHybrid;
    SurfaceFluxMoments electron_flux_moments{};
    SurfaceFluxMoments ion_flux_moments{};
    double fallback_electron_energy_ev = 1.0;
    double fallback_ion_energy_ev = 1.0;
    double local_reference_shift_v = 0.0;
    double incidence_scale = 1.0;
    double normalized_flow_alignment = 1.0;
    double local_field_factor = 1.0;
    double wake_factor = 0.0;
    double axisym_anisotropy = 1.0;
    double electron_recollection_ratio = 1.0;
    double ion_recollection_ratio = 1.0;
    bool patch_role = false;
    std::size_t electron_population_count = 1;
    std::size_t ion_population_count = 1;
};

struct SurfaceDistributionSynthesisResult
{
    double electron_collection_scale = 1.0;
    double ion_collection_scale = 1.0;
    double electron_characteristic_energy_ev = 1.0;
    double ion_characteristic_energy_ev = 1.0;
    double photo_incidence_scale = 1.0;
    bool valid = false;
};

struct SurfacePopulationDescriptor
{
    double density_m3 = 0.0;
    double temperature_ev = 0.0;
};

struct SurfaceDistributionBuildRequest
{
    SurfaceDistributionModel model = SurfaceDistributionModel::MultiPopulationHybrid;
    double electron_density_m3 = 0.0;
    double electron_temperature_ev = 0.0;
    double ion_density_m3 = 0.0;
    double ion_temperature_ev = 0.0;
    double ion_directed_velocity_m_per_s = 0.0;
    double axisym_anisotropy = 1.0;
    double electron_recollection_ratio = 1.0;
    double ion_recollection_ratio = 1.0;
    std::vector<SurfacePopulationDescriptor> electron_populations;
    std::vector<SurfacePopulationDescriptor> ion_populations;
};

struct SurfaceDistributionBuildResult
{
    SurfaceFluxMoments electron_flux_moments{};
    SurfaceFluxMoments ion_flux_moments{};
    bool valid = false;
};

struct SurfaceDistributionBundleRequest
{
    SurfaceDistributionBuildRequest build_request{};
    SurfaceDistributionSynthesisRequest synthesis_request{};
};

struct SurfaceDistributionEnvironment
{
    SurfaceDistributionModel model = SurfaceDistributionModel::MultiPopulationHybrid;
    double electron_density_m3 = 0.0;
    double electron_temperature_ev = 0.0;
    double ion_density_m3 = 0.0;
    double ion_temperature_ev = 0.0;
    double ion_directed_velocity_m_per_s = 0.0;
    std::vector<SurfacePopulationDescriptor> electron_populations;
    std::vector<SurfacePopulationDescriptor> ion_populations;
    double fallback_electron_energy_ev = 1.0;
    double fallback_ion_energy_ev = 1.0;
    double local_reference_shift_v = 0.0;
    double incidence_scale = 1.0;
    double normalized_flow_alignment = 1.0;
    double local_field_factor = 1.0;
    double wake_factor = 0.0;
    double axisym_anisotropy = 1.0;
    double electron_recollection_ratio = 1.0;
    double ion_recollection_ratio = 1.0;
    bool patch_role = false;
};

struct SurfaceDistributionRoleInputs
{
    SurfaceDistributionModel model = SurfaceDistributionModel::MultiPopulationHybrid;
    double electron_density_m3 = 0.0;
    double electron_temperature_ev = 0.0;
    double ion_density_m3 = 0.0;
    double ion_temperature_ev = 0.0;
    double ion_directed_velocity_m_per_s = 0.0;
    std::vector<SurfacePopulationDescriptor> electron_populations;
    std::vector<SurfacePopulationDescriptor> ion_populations;
    double fallback_electron_energy_ev = 1.0;
    double fallback_ion_energy_ev = 1.0;
    double local_reference_shift_v = 0.0;
    double patch_incidence_angle_deg = 0.0;
    double bulk_flow_velocity_m_per_s = 0.0;
    double projected_speed_m_per_s = 0.0;
    double flow_alignment_cosine = 1.0;
    double normal_electric_field_v_per_m = 0.0;
    double axisym_anisotropy = 1.0;
    double electron_recollection_ratio = 1.0;
    double ion_recollection_ratio = 1.0;
    bool share_patch_distribution = false;
    bool patch_role = false;
};

  struct SurfaceDistributionRoleBuildRequest
  {
    SurfaceDistributionModel model = SurfaceDistributionModel::MultiPopulationHybrid;
    double electron_density_m3 = 0.0;
    double electron_temperature_ev = 0.0;
    double ion_density_m3 = 0.0;
    double ion_temperature_ev = 0.0;
    double ion_directed_velocity_m_per_s = 0.0;
    ResolvedSpectrum electron_spectrum{};
    ResolvedSpectrum ion_spectrum{};
    double local_reference_shift_v = 0.0;
    double patch_incidence_angle_deg = 0.0;
    double bulk_flow_velocity_m_per_s = 0.0;
    double projected_speed_m_per_s = 0.0;
    double flow_alignment_cosine = 1.0;
    double normal_electric_field_v_per_m = 0.0;
    double axisym_anisotropy = 1.0;
    double electron_recollection_ratio = 1.0;
    double ion_recollection_ratio = 1.0;
    bool share_patch_distribution = false;
    bool patch_role = false;
  };

struct SurfaceDistributionBundleResult
{
    SurfaceDistributionBuildResult built_flux{};
    SurfaceDistributionSynthesisResult synthesized{};
    bool valid = false;
};

struct ResolvedSpectrumMoments
{
  double density_m3 = 0.0;
  double characteristic_energy_ev = 0.0;
  double average_mass_amu = 1.0;
  bool used_population_moments = false;
  bool valid = false;
};

struct LegacyPopulation
{
    double density_m3 = 0.0;
    double temperature_ev = 0.0;
    double mass_amu = 1.0;
};

struct LegacyFluxTableBuildRequest
{
    ResolvedSpectrum electron_spectrum{};
    ResolvedSpectrum ion_spectrum{};
    bool has_electron_spectrum = false;
    bool has_ion_spectrum = false;
    double fallback_electron_density_m3 = 0.0;
    double fallback_electron_temperature_ev = 0.0;
    double fallback_electron_mass_amu = 9.1093837015e-31 / 1.66053906660e-27;
    double fallback_ion_density_m3 = 0.0;
    double fallback_ion_temperature_ev = 0.0;
    double fallback_ion_mass_amu = 1.0;
    std::size_t max_population_count = 3;
    int sample_count = 1000;
};

struct LegacyFluxTables
{
    std::vector<double> energy_ev;
    std::vector<double> width_ev;
    std::vector<double> spectrum_energy_grid_ev;
    std::vector<double> smoothed_electron_flux;
    std::vector<double> smoothed_ion_flux;
    std::vector<double> total_electron_flux;
    std::vector<double> spectrum_electron_flux;
    std::vector<double> total_ion_flux;
    std::vector<double> spectrum_ion_flux;
    std::vector<std::vector<double>> electron_component_flux;
    std::vector<std::vector<double>> ion_component_flux;
    std::vector<LegacyPopulation> electron_populations;
    std::vector<LegacyPopulation> ion_populations;
    double average_spectrum_electron_temperature_ev = 0.0;
    double average_spectrum_ion_temperature_ev = 0.0;
};

class SurfaceDistributionFunction
{
  public:
    virtual ~SurfaceDistributionFunction() = default;
    virtual const char* family() const = 0;
    virtual SurfaceFluxMoments computeMoments() const = 0;
};

class IsotropicMaxwellianFluxDistribution final : public SurfaceDistributionFunction
{
  public:
    IsotropicMaxwellianFluxDistribution(double density_m3, double temperature_ev,
                                        double drift_speed_m_per_s = 0.0);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    double density_m3_ = 0.0;
    double temperature_ev_ = 0.0;
    double drift_speed_m_per_s_ = 0.0;
};

class BiMaxwellianFluxDistribution final : public SurfaceDistributionFunction
{
  public:
    BiMaxwellianFluxDistribution(double density_primary_m3, double temperature_primary_ev,
                                 double density_secondary_m3, double temperature_secondary_ev);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    double density_primary_m3_ = 0.0;
    double temperature_primary_ev_ = 0.0;
    double density_secondary_m3_ = 0.0;
    double temperature_secondary_ev_ = 0.0;
};

class TabulatedSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    TabulatedSurfaceDistributionFunction(std::vector<double> energies_ev,
                                         std::vector<double> flux_weights,
                                         std::vector<double> cosine_weights = {});

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

    std::size_t sampleCount() const { return energies_ev_.size(); }

  private:
    std::vector<double> energies_ev_;
    std::vector<double> flux_weights_;
    std::vector<double> cosine_weights_;
};

class TabulatedVelocitySurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    TabulatedVelocitySurfaceDistributionFunction(std::vector<double> normal_velocity_m_per_s,
                                                 std::vector<double> flux_weights,
                                                 double energy_scale_ev_per_velocity2,
                                                 std::vector<double> cosine_weights = {});

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    std::vector<double> normal_velocity_m_per_s_;
    std::vector<double> flux_weights_;
    std::vector<double> cosine_weights_;
    double energy_scale_ev_per_velocity2_ = 1.0e-10;
};

class NonLocalizedHybridSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    NonLocalizedHybridSurfaceDistributionFunction(const SurfaceFluxMoments& local_moments,
                                                  const SurfaceFluxMoments& non_local_moments,
                                                  double non_local_weight);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    SurfaceFluxMoments local_moments_{};
    SurfaceFluxMoments non_local_moments_{};
    double non_local_weight_ = 0.0;
};

class MultipleSurfSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    MultipleSurfSurfaceDistributionFunction(std::vector<SurfaceFluxMoments> component_moments,
                                            std::vector<double> component_weights,
                                            double localization_bias = 0.0);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    std::vector<SurfaceFluxMoments> component_moments_;
    std::vector<double> component_weights_;
    double localization_bias_ = 0.0;
};

class LocalModifiedPearsonIVSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    LocalModifiedPearsonIVSurfaceDistributionFunction(double density_m3,
                                                      double characteristic_energy_ev,
                                                      double skewness,
                                                      double kurtosis,
                                                      double incidence_bias = 0.5,
                                                      double localization_gain = 1.0);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    double density_m3_ = 0.0;
    double characteristic_energy_ev_ = 0.0;
    double skewness_ = 0.0;
    double kurtosis_ = 3.0;
    double incidence_bias_ = 0.5;
    double localization_gain_ = 1.0;
};

class LocalTabulatedSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    LocalTabulatedSurfaceDistributionFunction(std::vector<double> energies_ev,
                                             std::vector<double> flux_weights,
                                             double local_energy_shift_ev = 0.0,
                                             double local_flux_scale = 1.0,
                                             std::vector<double> cosine_weights = {});

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    std::vector<double> energies_ev_;
    std::vector<double> flux_weights_;
    std::vector<double> cosine_weights_;
    double local_energy_shift_ev_ = 0.0;
    double local_flux_scale_ = 1.0;
};

class TwoAxesTabulatedVelocitySurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    TwoAxesTabulatedVelocitySurfaceDistributionFunction(
        std::vector<double> normal_velocity_m_per_s,
        std::vector<double> tangential_velocity_m_per_s,
        std::vector<double> flux_weights,
        double normal_energy_scale_ev_per_velocity2,
        double tangential_energy_scale_ev_per_velocity2,
        std::vector<double> cosine_weights = {});

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    std::vector<double> normal_velocity_m_per_s_;
    std::vector<double> tangential_velocity_m_per_s_;
    std::vector<double> flux_weights_;
    std::vector<double> cosine_weights_;
    double normal_energy_scale_ev_per_velocity2_ = 1.0e-10;
    double tangential_energy_scale_ev_per_velocity2_ = 1.0e-10;
};

class AxisymTabulatedVelocitySurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    AxisymTabulatedVelocitySurfaceDistributionFunction(
        std::vector<double> axial_velocity_m_per_s,
        std::vector<double> radial_velocity_m_per_s,
        std::vector<double> flux_weights,
        double axial_energy_scale_ev_per_velocity2,
        double radial_energy_scale_ev_per_velocity2,
        std::vector<double> cosine_weights = {});

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    std::vector<double> axial_velocity_m_per_s_;
    std::vector<double> radial_velocity_m_per_s_;
    std::vector<double> flux_weights_;
    std::vector<double> cosine_weights_;
    double axial_energy_scale_ev_per_velocity2_ = 1.0e-10;
    double radial_energy_scale_ev_per_velocity2_ = 1.0e-10;
};

class UniformVelocitySurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    UniformVelocitySurfaceDistributionFunction(double velocity_m_per_s,
                                               double flux_weight,
                                               double energy_scale_ev_per_velocity2,
                                               double cosine_weight = 1.0);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    double velocity_m_per_s_ = 0.0;
    double flux_weight_ = 0.0;
    double energy_scale_ev_per_velocity2_ = 1.0e-10;
    double cosine_weight_ = 1.0;
};

class FluidSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    FluidSurfaceDistributionFunction(double density_m3,
                                     double characteristic_energy_ev,
                                     double directed_velocity_m_per_s,
                                     double compressibility = 1.0,
                                     double incidence_bias = 0.5);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    double density_m3_ = 0.0;
    double characteristic_energy_ev_ = 0.0;
    double directed_velocity_m_per_s_ = 0.0;
    double compressibility_ = 1.0;
    double incidence_bias_ = 0.5;
};

class FowlerNordheimSurfaceDistributionFunction final : public SurfaceDistributionFunction
{
  public:
    FowlerNordheimSurfaceDistributionFunction(double field_strength_v_per_m,
                                              double work_function_ev,
                                              double emission_area_scale = 1.0);

    const char* family() const override;
    SurfaceFluxMoments computeMoments() const override;

  private:
    double field_strength_v_per_m_ = 0.0;
    double work_function_ev_ = 4.5;
    double emission_area_scale_ = 1.0;
};

class SurfaceFluxSampler
{
  public:
    void markDirty() { dirty_ = true; }
    void cacheMoments(const SurfaceFluxMoments& moments)
    {
        cached_moments_ = moments;
        dirty_ = false;
    }

    bool dirty() const { return dirty_; }
    const SurfaceFluxMoments& cachedMoments() const { return cached_moments_; }

  private:
    bool dirty_ = true;
    SurfaceFluxMoments cached_moments_{};
};

class SurfaceEmissionSampler
{
  public:
    void setSuperParticleDensification(double value) { densification_ = value > 0.0 ? value : 1.0; }
    double superParticleDensification() const { return densification_; }

    void noteEmissionCount(double emitted_particles) { cached_emission_count_ = emitted_particles; }
    double cachedEmissionCount() const { return cached_emission_count_; }

  private:
    double densification_ = 1.0;
    double cached_emission_count_ = 0.0;
};

SurfaceDistributionSynthesisResult synthesizeSurfaceDistributionMoments(
    const SurfaceDistributionSynthesisRequest& request);
SurfaceDistributionBuildResult buildSurfaceFluxMoments(
    const SurfaceDistributionBuildRequest& request);
SurfaceDistributionBundleResult evaluateSurfaceDistributionBundle(
    const SurfaceDistributionBundleRequest& request);
SurfaceDistributionBundleRequest makeSurfaceDistributionBundleRequest(
    const SurfaceDistributionEnvironment& environment);
SurfaceDistributionBundleResult evaluateSurfaceDistributionEnvironment(
    const SurfaceDistributionEnvironment& environment);
SurfaceDistributionEnvironment makeSurfaceDistributionEnvironment(
    const SurfaceDistributionRoleInputs& inputs);
SurfaceDistributionRoleInputs makeSurfaceDistributionRoleInputs(
  const SurfaceDistributionRoleBuildRequest& request);
bool usesResolvedSpectrumDiscreteFlux(const ResolvedSpectrum& spectrum);
double integrateResolvedSpectrumFlux(
  const ResolvedSpectrum& spectrum,
  const std::function<double(double, double)>& integrand);
double integratePopulationMaxwellianAverageEv(
  double temperature_ev,
  const std::function<double(double)>& evaluator);
double resolvedSpectrumAverageEnergyEv(const ResolvedSpectrum& spectrum,
                     double fallback_energy_ev = 1.0);
ResolvedSpectrumMoments resolveSpectrumMoments(const ResolvedSpectrum& spectrum,
                         double fallback_density_m3,
                         double fallback_characteristic_energy_ev,
                         double fallback_mass_amu = 1.0);
double resolvedSpectrumDensity(const ResolvedSpectrum& spectrum,
                double fallback_density_m3);
double resolvedSpectrumCharacteristicEnergyEv(
  const ResolvedSpectrum& spectrum, double fallback_characteristic_energy_ev);
double resolvedSpectrumAverageMassAmu(const ResolvedSpectrum& spectrum,
                    double fallback_mass_amu);
double resolvedSpectrumAverageMassKg(const ResolvedSpectrum& spectrum,
                                     double fallback_mass_kg);
double sampleResolvedSpectrumEnergyEv(const ResolvedSpectrum& spectrum,
                                      std::mt19937& rng,
                                      double fallback_energy_ev = 1.0);
double interpolateLegacySpectrum(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 double target);
std::vector<double> smoothLegacySpectrumLine3(const std::vector<double>& x,
                                              const std::vector<double>& y);
double legacyElectronComponentFlux(double density_m3,
                                   double temperature_ev,
                                   double energy_ev);
double legacyIonComponentFlux(double density_m3,
                              double temperature_ev,
                              double mass_amu,
                              double energy_ev);
LegacyFluxTables buildLegacyFluxTables(const LegacyFluxTableBuildRequest& request);

} // namespace Particle
} // namespace SCDAT
