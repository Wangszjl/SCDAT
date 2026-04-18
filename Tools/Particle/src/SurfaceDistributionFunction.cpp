#include "../include/SurfaceDistributionFunction.h"
#include "../include/ParticleSource.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <random>

namespace SCDAT
{
namespace Particle
{

namespace
{

SurfaceFluxMoments computeWeightedMoments(const std::vector<double>& energies_ev,
                                          const std::vector<double>& weights,
                                          const std::vector<double>& cosine_weights)
{
    SurfaceFluxMoments moments;
    if (energies_ev.empty() || energies_ev.size() != weights.size())
    {
        return moments;
    }

    const double total_weight =
        std::accumulate(weights.begin(), weights.end(), 0.0, [](double lhs, double rhs) {
            return lhs + std::max(0.0, rhs);
        });
    if (total_weight <= 0.0)
    {
        return moments;
    }

    double weighted_energy = 0.0;
    double weighted_cosine = 0.0;
    for (std::size_t index = 0; index < energies_ev.size(); ++index)
    {
        const double weight = std::max(0.0, weights[index]);
        weighted_energy += weight * std::max(0.0, energies_ev[index]);
        const double cosine = cosine_weights.empty()
                                  ? 1.0
                                  : std::clamp(cosine_weights[index], 0.0, 1.0);
        weighted_cosine += weight * cosine;
    }

    moments.total_flux = total_weight;
    moments.mean_energy_ev = weighted_energy / total_weight;
    moments.mean_cos_incidence = weighted_cosine / total_weight;
    moments.valid = true;
    return moments;
}

SurfaceFluxMoments blendFluxMoments(const SurfaceFluxMoments& local_moments,
                                    const SurfaceFluxMoments& non_local_moments,
                                    double non_local_weight)
{
    SurfaceFluxMoments moments;
    const double weight = std::clamp(non_local_weight, 0.0, 1.0);
    if (!local_moments.valid && !non_local_moments.valid)
    {
        return moments;
    }
    if (!local_moments.valid)
    {
        moments = non_local_moments;
        return moments;
    }
    if (!non_local_moments.valid)
    {
        moments = local_moments;
        return moments;
    }

    const double local_flux = std::max(0.0, local_moments.total_flux);
    const double non_local_flux = std::max(0.0, non_local_moments.total_flux);
    const double blended_local_flux = (1.0 - weight) * local_flux;
    const double blended_non_local_flux = weight * non_local_flux;
    const double total_flux = blended_local_flux + blended_non_local_flux;
    if (total_flux <= 0.0)
    {
        return moments;
    }

    moments.total_flux = total_flux;
    moments.mean_energy_ev =
        (blended_local_flux * std::max(0.0, local_moments.mean_energy_ev) +
         blended_non_local_flux * std::max(0.0, non_local_moments.mean_energy_ev)) /
        total_flux;
    moments.mean_cos_incidence =
        std::clamp((blended_local_flux * std::clamp(local_moments.mean_cos_incidence, 0.0, 1.0) +
                    blended_non_local_flux *
                        std::clamp(non_local_moments.mean_cos_incidence, 0.0, 1.0)) /
                       total_flux,
                   0.0, 1.0);
    moments.valid = true;
    return moments;
}

SurfaceFluxMoments computeParametricMaxwellMoments(double density_m3, double temperature_ev,
                                                   double directed_velocity_m_per_s,
                                                   double energy_multiplier,
                                                   double flux_multiplier,
                                                   double boltzmann_bias,
                                                   double cosine_floor,
                                                   double cosine_gain)
{
    SurfaceFluxMoments moments;
    const double density = std::max(0.0, density_m3);
    const double temperature = std::max(1.0e-6, temperature_ev);
    if (!(density > 0.0))
    {
        return moments;
    }

    const double drift_factor =
        1.0 + 0.25 * std::clamp(std::abs(directed_velocity_m_per_s) / 5000.0, 0.0, 4.0);
    const double thermal_factor =
        std::clamp(1.0 + boltzmann_bias * std::sqrt(temperature) / 8.0, 0.25, 4.0);
    moments.total_flux =
        density * std::max(0.05, flux_multiplier) * drift_factor * thermal_factor;
    moments.mean_energy_ev = std::max(1.0e-6, energy_multiplier * temperature) *
                             std::clamp(1.0 + 0.08 * boltzmann_bias, 0.4, 2.5);
    moments.mean_cos_incidence = std::clamp(
        cosine_floor + cosine_gain *
                           std::clamp(std::abs(directed_velocity_m_per_s) /
                                          std::max(1.0, std::abs(directed_velocity_m_per_s) +
                                                            2.5e4 * std::sqrt(temperature)),
                                      0.0, 1.0),
        0.05, 0.99);
    moments.valid = true;
    return moments;
}

double clampRecollectionRatio(double ratio)
{
    if (!std::isfinite(ratio))
    {
        return 1.0;
    }
    return std::clamp(ratio, 0.0, 1.0);
}

SurfaceFluxMoments scaleFluxMomentsForRecollection(const SurfaceFluxMoments& moments, double ratio)
{
    if (!moments.valid)
    {
        return moments;
    }

    SurfaceFluxMoments scaled = moments;
    scaled.total_flux *= clampRecollectionRatio(ratio);
    scaled.valid = true;
    return scaled;
}

SurfaceFluxMoments combineMultipleSurfaceMoments(
    const std::vector<SurfaceFluxMoments>& component_moments,
    const std::vector<double>& component_weights,
    double localization_bias)
{
    SurfaceFluxMoments moments;
    if (component_moments.empty())
    {
        return moments;
    }

    double total_weight = 0.0;
    double weighted_flux = 0.0;
    double weighted_energy = 0.0;
    double weighted_cosine = 0.0;
    const double clamped_bias = std::clamp(localization_bias, -1.0, 1.0);

    for (std::size_t i = 0; i < component_moments.size(); ++i)
    {
        const auto& component = component_moments[i];
        if (!component.valid)
        {
            continue;
        }

        const double base_weight =
            i < component_weights.size() ? std::max(0.0, component_weights[i]) : 1.0;
        const double locality_scale = 1.0 + clamped_bias * (static_cast<double>(i) + 1.0) /
                                                static_cast<double>(component_moments.size());
        const double component_weight = std::max(0.0, base_weight * locality_scale);
        if (!(component_weight > 0.0))
        {
            continue;
        }

        const double component_flux = std::max(0.0, component.total_flux);
        total_weight += component_weight;
        weighted_flux += component_weight * component_flux;
        weighted_energy += component_weight * component_flux *
                           std::max(0.0, component.mean_energy_ev);
        weighted_cosine += component_weight * component_flux *
                           std::clamp(component.mean_cos_incidence, 0.0, 1.0);
    }

    if (!(weighted_flux > 0.0) || !(total_weight > 0.0))
    {
        return moments;
    }

    moments.total_flux = weighted_flux / total_weight;
    moments.mean_energy_ev = weighted_energy / weighted_flux;
    moments.mean_cos_incidence = std::clamp(weighted_cosine / weighted_flux, 0.0, 1.0);
    moments.valid = true;
    return moments;
}

std::vector<LegacyPopulation> collectLegacyPopulations(const ResolvedSpectrum& spectrum,
                                                       bool has_spectrum,
                                                       double fallback_density_m3,
                                                       double fallback_temperature_ev,
                                                       double fallback_mass_amu,
                                                       std::size_t max_population_count)
{
    std::vector<LegacyPopulation> populations;
    if (has_spectrum && !spectrum.populations.empty())
    {
        for (const auto& population : spectrum.populations)
        {
            populations.push_back(
                {population.density_m3, population.temperature_ev, population.mass_amu});
            if (populations.size() >= max_population_count)
            {
                break;
            }
        }
    }

    if (populations.empty())
    {
        populations.push_back(
            {fallback_density_m3, fallback_temperature_ev, fallback_mass_amu});
    }
    return populations;
}

} // namespace

IsotropicMaxwellianFluxDistribution::IsotropicMaxwellianFluxDistribution(
    double density_m3, double temperature_ev, double drift_speed_m_per_s)
    : density_m3_(density_m3), temperature_ev_(temperature_ev),
      drift_speed_m_per_s_(drift_speed_m_per_s)
{
}

const char* IsotropicMaxwellianFluxDistribution::family() const
{
    return "spis_isotropic_maxwellian_flux_v1";
}

SurfaceFluxMoments IsotropicMaxwellianFluxDistribution::computeMoments() const
{
    SurfaceFluxMoments moments;
    moments.total_flux = std::max(0.0, density_m3_) * (1.0 + std::max(0.0, drift_speed_m_per_s_) * 1.0e-6);
    moments.mean_energy_ev = std::max(0.0, 2.0 * temperature_ev_);
    moments.mean_cos_incidence = 0.5;
    moments.valid = moments.total_flux > 0.0;
    return moments;
}

BiMaxwellianFluxDistribution::BiMaxwellianFluxDistribution(double density_primary_m3,
                                                           double temperature_primary_ev,
                                                           double density_secondary_m3,
                                                           double temperature_secondary_ev)
    : density_primary_m3_(density_primary_m3), temperature_primary_ev_(temperature_primary_ev),
      density_secondary_m3_(density_secondary_m3),
      temperature_secondary_ev_(temperature_secondary_ev)
{
}

const char* BiMaxwellianFluxDistribution::family() const
{
    return "spis_bimaxwellian_flux_v1";
}

SurfaceFluxMoments BiMaxwellianFluxDistribution::computeMoments() const
{
    const double total_density = std::max(0.0, density_primary_m3_) + std::max(0.0, density_secondary_m3_);
    SurfaceFluxMoments moments;
    if (total_density <= 0.0)
    {
        return moments;
    }
    moments.total_flux = total_density;
    moments.mean_energy_ev =
        (std::max(0.0, density_primary_m3_) * std::max(0.0, 2.0 * temperature_primary_ev_) +
         std::max(0.0, density_secondary_m3_) * std::max(0.0, 2.0 * temperature_secondary_ev_)) /
        total_density;
    moments.mean_cos_incidence = 0.5;
    moments.valid = true;
    return moments;
}

TabulatedSurfaceDistributionFunction::TabulatedSurfaceDistributionFunction(
    std::vector<double> energies_ev, std::vector<double> flux_weights,
    std::vector<double> cosine_weights)
    : energies_ev_(std::move(energies_ev)), flux_weights_(std::move(flux_weights)),
      cosine_weights_(std::move(cosine_weights))
{
    if (!cosine_weights_.empty() && cosine_weights_.size() != energies_ev_.size())
    {
        cosine_weights_.clear();
    }
}

const char* TabulatedSurfaceDistributionFunction::family() const
{
    return "spis_tabulated_surface_flux_v1";
}

SurfaceFluxMoments TabulatedSurfaceDistributionFunction::computeMoments() const
{
    return computeWeightedMoments(energies_ev_, flux_weights_, cosine_weights_);
}

TabulatedVelocitySurfaceDistributionFunction::TabulatedVelocitySurfaceDistributionFunction(
    std::vector<double> normal_velocity_m_per_s, std::vector<double> flux_weights,
    double energy_scale_ev_per_velocity2, std::vector<double> cosine_weights)
    : normal_velocity_m_per_s_(std::move(normal_velocity_m_per_s)),
      flux_weights_(std::move(flux_weights)),
      cosine_weights_(std::move(cosine_weights)),
      energy_scale_ev_per_velocity2_(energy_scale_ev_per_velocity2)
{
    if (!cosine_weights_.empty() && cosine_weights_.size() != normal_velocity_m_per_s_.size())
    {
        cosine_weights_.clear();
    }
}

const char* TabulatedVelocitySurfaceDistributionFunction::family() const
{
    return "spis_tabulated_velocity_surface_flux_v1";
}

SurfaceFluxMoments TabulatedVelocitySurfaceDistributionFunction::computeMoments() const
{
    std::vector<double> energies_ev;
    energies_ev.reserve(normal_velocity_m_per_s_.size());
    const double energy_scale = std::max(1.0e-16, std::abs(energy_scale_ev_per_velocity2_));
    for (const double velocity_m_per_s : normal_velocity_m_per_s_)
    {
        const double clamped_velocity = std::max(0.0, velocity_m_per_s);
        energies_ev.push_back(energy_scale * clamped_velocity * clamped_velocity);
    }
    return computeWeightedMoments(energies_ev, flux_weights_, cosine_weights_);
}

NonLocalizedHybridSurfaceDistributionFunction::NonLocalizedHybridSurfaceDistributionFunction(
    const SurfaceFluxMoments& local_moments, const SurfaceFluxMoments& non_local_moments,
    double non_local_weight)
    : local_moments_(local_moments),
      non_local_moments_(non_local_moments),
      non_local_weight_(non_local_weight)
{
}

const char* NonLocalizedHybridSurfaceDistributionFunction::family() const
{
    return "spis_non_localized_surface_flux_v1";
}

SurfaceFluxMoments NonLocalizedHybridSurfaceDistributionFunction::computeMoments() const
{
    return blendFluxMoments(local_moments_, non_local_moments_, non_local_weight_);
}

MultipleSurfSurfaceDistributionFunction::MultipleSurfSurfaceDistributionFunction(
    std::vector<SurfaceFluxMoments> component_moments, std::vector<double> component_weights,
    double localization_bias)
    : component_moments_(std::move(component_moments)),
      component_weights_(std::move(component_weights)),
      localization_bias_(localization_bias)
{
}

const char* MultipleSurfSurfaceDistributionFunction::family() const
{
    return "spis_multiple_surface_flux_v1";
}

SurfaceFluxMoments MultipleSurfSurfaceDistributionFunction::computeMoments() const
{
    return combineMultipleSurfaceMoments(component_moments_, component_weights_,
                                         localization_bias_);
}

LocalModifiedPearsonIVSurfaceDistributionFunction::LocalModifiedPearsonIVSurfaceDistributionFunction(
    double density_m3, double characteristic_energy_ev, double skewness, double kurtosis,
    double incidence_bias, double localization_gain)
    : density_m3_(density_m3), characteristic_energy_ev_(characteristic_energy_ev),
      skewness_(skewness), kurtosis_(kurtosis), incidence_bias_(incidence_bias),
      localization_gain_(localization_gain)
{
}

const char* LocalModifiedPearsonIVSurfaceDistributionFunction::family() const
{
    return "spis_local_modified_pearson_iv_surface_flux_v1";
}

SurfaceFluxMoments LocalModifiedPearsonIVSurfaceDistributionFunction::computeMoments() const
{
    SurfaceFluxMoments moments;
    const double density = std::max(0.0, density_m3_);
    if (!(density > 0.0))
    {
        return moments;
    }

    const double energy = std::max(1.0e-6, characteristic_energy_ev_);
    const double skew_term = std::clamp(skewness_, -3.0, 3.0);
    const double kurtosis_term = std::clamp(kurtosis_, 1.5, 12.0);
    const double localization_term = std::clamp(localization_gain_, 0.25, 4.0);

    moments.total_flux = density * localization_term;
    moments.mean_energy_ev = energy *
                             std::clamp(1.0 + 0.18 * skew_term +
                                            0.04 * (kurtosis_term - 3.0),
                                        0.15, 6.0);
    moments.mean_cos_incidence =
        std::clamp(0.5 + 0.18 * std::clamp(incidence_bias_, -1.0, 1.0) -
                       0.03 * std::max(0.0, kurtosis_term - 3.0),
                   0.05, 0.98);
    moments.valid = true;
    return moments;
}

LocalTabulatedSurfaceDistributionFunction::LocalTabulatedSurfaceDistributionFunction(
    std::vector<double> energies_ev, std::vector<double> flux_weights,
    double local_energy_shift_ev, double local_flux_scale, std::vector<double> cosine_weights)
    : energies_ev_(std::move(energies_ev)), flux_weights_(std::move(flux_weights)),
      cosine_weights_(std::move(cosine_weights)),
      local_energy_shift_ev_(local_energy_shift_ev), local_flux_scale_(local_flux_scale)
{
    if (!cosine_weights_.empty() && cosine_weights_.size() != energies_ev_.size())
    {
        cosine_weights_.clear();
    }
}

const char* LocalTabulatedSurfaceDistributionFunction::family() const
{
    return "spis_local_tabulated_surface_flux_v1";
}

SurfaceFluxMoments LocalTabulatedSurfaceDistributionFunction::computeMoments() const
{
    std::vector<double> shifted_energies;
    shifted_energies.reserve(energies_ev_.size());
    for (const double energy_ev : energies_ev_)
    {
        shifted_energies.push_back(std::max(0.0, energy_ev + local_energy_shift_ev_));
    }

    std::vector<double> scaled_flux_weights;
    scaled_flux_weights.reserve(flux_weights_.size());
    const double flux_scale = std::clamp(local_flux_scale_, 1.0e-6, 1.0e3);
    for (const double weight : flux_weights_)
    {
        scaled_flux_weights.push_back(std::max(0.0, weight) * flux_scale);
    }

    return computeWeightedMoments(shifted_energies, scaled_flux_weights, cosine_weights_);
}

TwoAxesTabulatedVelocitySurfaceDistributionFunction::
    TwoAxesTabulatedVelocitySurfaceDistributionFunction(
        std::vector<double> normal_velocity_m_per_s,
        std::vector<double> tangential_velocity_m_per_s,
        std::vector<double> flux_weights,
        double normal_energy_scale_ev_per_velocity2,
        double tangential_energy_scale_ev_per_velocity2,
        std::vector<double> cosine_weights)
    : normal_velocity_m_per_s_(std::move(normal_velocity_m_per_s)),
      tangential_velocity_m_per_s_(std::move(tangential_velocity_m_per_s)),
      flux_weights_(std::move(flux_weights)),
      cosine_weights_(std::move(cosine_weights)),
      normal_energy_scale_ev_per_velocity2_(normal_energy_scale_ev_per_velocity2),
      tangential_energy_scale_ev_per_velocity2_(tangential_energy_scale_ev_per_velocity2)
{
    if (tangential_velocity_m_per_s_.size() != normal_velocity_m_per_s_.size())
    {
        tangential_velocity_m_per_s_.assign(normal_velocity_m_per_s_.size(), 0.0);
    }
    if (!cosine_weights_.empty() && cosine_weights_.size() != normal_velocity_m_per_s_.size())
    {
        cosine_weights_.clear();
    }
}

const char* TwoAxesTabulatedVelocitySurfaceDistributionFunction::family() const
{
    return "spis_two_axes_tabulated_velocity_surface_flux_v1";
}

SurfaceFluxMoments TwoAxesTabulatedVelocitySurfaceDistributionFunction::computeMoments() const
{
    std::vector<double> energies_ev;
    energies_ev.reserve(normal_velocity_m_per_s_.size());
    const double normal_scale =
        std::max(1.0e-16, std::abs(normal_energy_scale_ev_per_velocity2_));
    const double tangential_scale =
        std::max(1.0e-16, std::abs(tangential_energy_scale_ev_per_velocity2_));
    for (std::size_t i = 0; i < normal_velocity_m_per_s_.size(); ++i)
    {
        const double normal_velocity = std::max(0.0, normal_velocity_m_per_s_[i]);
        const double tangential_velocity = std::max(0.0, tangential_velocity_m_per_s_[i]);
        energies_ev.push_back(normal_scale * normal_velocity * normal_velocity +
                              tangential_scale * tangential_velocity * tangential_velocity);
    }
    return computeWeightedMoments(energies_ev, flux_weights_, cosine_weights_);
}

AxisymTabulatedVelocitySurfaceDistributionFunction::
    AxisymTabulatedVelocitySurfaceDistributionFunction(
        std::vector<double> axial_velocity_m_per_s,
        std::vector<double> radial_velocity_m_per_s,
        std::vector<double> flux_weights,
        double axial_energy_scale_ev_per_velocity2,
        double radial_energy_scale_ev_per_velocity2,
        std::vector<double> cosine_weights)
    : axial_velocity_m_per_s_(std::move(axial_velocity_m_per_s)),
      radial_velocity_m_per_s_(std::move(radial_velocity_m_per_s)),
      flux_weights_(std::move(flux_weights)),
      cosine_weights_(std::move(cosine_weights)),
      axial_energy_scale_ev_per_velocity2_(axial_energy_scale_ev_per_velocity2),
      radial_energy_scale_ev_per_velocity2_(radial_energy_scale_ev_per_velocity2)
{
    if (radial_velocity_m_per_s_.size() != axial_velocity_m_per_s_.size())
    {
        radial_velocity_m_per_s_.assign(axial_velocity_m_per_s_.size(), 0.0);
    }
    if (!cosine_weights_.empty() && cosine_weights_.size() != axial_velocity_m_per_s_.size())
    {
        cosine_weights_.clear();
    }
}

const char* AxisymTabulatedVelocitySurfaceDistributionFunction::family() const
{
    return "spis_axisym_tabulated_velocity_surface_flux_v1";
}

SurfaceFluxMoments AxisymTabulatedVelocitySurfaceDistributionFunction::computeMoments() const
{
    std::vector<double> energies_ev;
    energies_ev.reserve(axial_velocity_m_per_s_.size());
    const double axial_scale =
        std::max(1.0e-16, std::abs(axial_energy_scale_ev_per_velocity2_));
    const double radial_scale =
        std::max(1.0e-16, std::abs(radial_energy_scale_ev_per_velocity2_));
    for (std::size_t i = 0; i < axial_velocity_m_per_s_.size(); ++i)
    {
        const double axial_velocity = std::max(0.0, axial_velocity_m_per_s_[i]);
        const double radial_velocity = std::max(0.0, radial_velocity_m_per_s_[i]);
        energies_ev.push_back(axial_scale * axial_velocity * axial_velocity +
                              radial_scale * radial_velocity * radial_velocity);
    }
    return computeWeightedMoments(energies_ev, flux_weights_, cosine_weights_);
}

UniformVelocitySurfaceDistributionFunction::UniformVelocitySurfaceDistributionFunction(
    double velocity_m_per_s, double flux_weight, double energy_scale_ev_per_velocity2,
    double cosine_weight)
    : velocity_m_per_s_(velocity_m_per_s), flux_weight_(flux_weight),
      energy_scale_ev_per_velocity2_(energy_scale_ev_per_velocity2),
      cosine_weight_(cosine_weight)
{
}

const char* UniformVelocitySurfaceDistributionFunction::family() const
{
    return "spis_uniform_velocity_surface_flux_v1";
}

SurfaceFluxMoments UniformVelocitySurfaceDistributionFunction::computeMoments() const
{
    return computeWeightedMoments(
        {std::max(1.0e-16, std::abs(energy_scale_ev_per_velocity2_)) *
         std::max(0.0, velocity_m_per_s_) * std::max(0.0, velocity_m_per_s_)},
        {std::max(0.0, flux_weight_)}, {std::clamp(cosine_weight_, 0.0, 1.0)});
}

FluidSurfaceDistributionFunction::FluidSurfaceDistributionFunction(
    double density_m3, double characteristic_energy_ev, double directed_velocity_m_per_s,
    double compressibility, double incidence_bias)
    : density_m3_(density_m3), characteristic_energy_ev_(characteristic_energy_ev),
      directed_velocity_m_per_s_(directed_velocity_m_per_s), compressibility_(compressibility),
      incidence_bias_(incidence_bias)
{
}

const char* FluidSurfaceDistributionFunction::family() const
{
    return "spis_fluid_surface_flux_v1";
}

SurfaceFluxMoments FluidSurfaceDistributionFunction::computeMoments() const
{
    SurfaceFluxMoments moments;
    const double density = std::max(0.0, density_m3_);
    if (!(density > 0.0))
    {
        return moments;
    }

    const double compressibility = std::clamp(compressibility_, 0.1, 10.0);
    const double velocity_gain =
        1.0 + 0.2 * std::clamp(std::abs(directed_velocity_m_per_s_) / 5000.0, 0.0, 5.0);
    moments.total_flux = density * compressibility * velocity_gain;
    moments.mean_energy_ev = std::max(1.0e-6, characteristic_energy_ev_) * velocity_gain;
    moments.mean_cos_incidence = std::clamp(
        0.5 + 0.25 * std::clamp(incidence_bias_, -1.0, 1.0) +
            0.10 * std::clamp(directed_velocity_m_per_s_ / 5000.0, -1.0, 1.0),
        0.05, 0.99);
    moments.valid = true;
    return moments;
}

FowlerNordheimSurfaceDistributionFunction::FowlerNordheimSurfaceDistributionFunction(
    double field_strength_v_per_m, double work_function_ev, double emission_area_scale)
    : field_strength_v_per_m_(field_strength_v_per_m), work_function_ev_(work_function_ev),
      emission_area_scale_(emission_area_scale)
{
}

const char* FowlerNordheimSurfaceDistributionFunction::family() const
{
    return "spis_fowler_nordheim_surface_flux_v1";
}

SurfaceFluxMoments FowlerNordheimSurfaceDistributionFunction::computeMoments() const
{
    SurfaceFluxMoments moments;
    const double field_strength = std::max(0.0, field_strength_v_per_m_);
    if (!(field_strength > 0.0))
    {
        return moments;
    }

    const double work_function = std::max(1.0, work_function_ev_);
    const double area_scale = std::clamp(emission_area_scale_, 1.0e-6, 1.0e6);
    const double reduced_field = field_strength / 1.0e7;
    const double transmission =
        std::exp(-std::clamp(work_function / std::max(1.0e-6, reduced_field), 0.0, 40.0));

    moments.total_flux = area_scale * reduced_field * reduced_field * transmission;
    moments.mean_energy_ev =
        std::max(1.0e-3, 0.45 * std::sqrt(std::max(1.0, work_function) * reduced_field));
    moments.mean_cos_incidence = 0.98;
    moments.valid = moments.total_flux > 0.0;
    return moments;
}

SurfaceDistributionSynthesisResult synthesizeSurfaceDistributionMoments(
    const SurfaceDistributionSynthesisRequest& request)
{
    SurfaceDistributionSynthesisResult result;
    result.electron_characteristic_energy_ev =
        std::max(std::max(1.0e-6, request.fallback_electron_energy_ev),
                 request.electron_flux_moments.valid ? request.electron_flux_moments.mean_energy_ev : 0.0);
    result.ion_characteristic_energy_ev =
        std::max(std::max(1.0e-6, request.fallback_ion_energy_ev),
                 request.ion_flux_moments.valid ? request.ion_flux_moments.mean_energy_ev : 0.0);

    switch (request.model)
    {
    case SurfaceDistributionModel::MaxwellianProjected:
        result.electron_collection_scale =
            request.electron_flux_moments.valid
                ? 0.75 + 0.5 * request.electron_flux_moments.mean_cos_incidence
                : 1.0;
        result.ion_collection_scale =
            request.patch_role
                ? std::max(0.2, (0.35 + 0.65 * request.incidence_scale) *
                                      (request.ion_flux_moments.valid
                                           ? 0.75 + 0.5 * request.ion_flux_moments.mean_cos_incidence
                                           : 1.0))
                : 1.0;
        result.photo_incidence_scale = request.incidence_scale;
        break;
    case SurfaceDistributionModel::WakeAnisotropic:
        result.electron_collection_scale =
            (0.85 + 0.15 * request.wake_factor) *
            (request.electron_flux_moments.valid
                 ? 0.5 + request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            (0.20 + 0.80 * (1.0 - request.wake_factor)) *
            (request.ion_flux_moments.valid ? 0.5 + request.ion_flux_moments.mean_cos_incidence
                                            : 1.0);
        result.photo_incidence_scale = request.incidence_scale * (0.8 + 0.2 * request.wake_factor);
        break;
    case SurfaceDistributionModel::TabulatedVelocity:
        result.electron_collection_scale =
            request.local_field_factor *
            (request.electron_flux_moments.valid
                 ? 0.65 + 0.7 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            (0.9 + 0.1 * std::max(0.0, request.normalized_flow_alignment)) *
            (request.ion_flux_moments.valid
                 ? 0.60 + 0.8 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * (0.85 + 0.15 * request.local_field_factor);
        break;
    case SurfaceDistributionModel::NonLocalizedHybrid:
    {
        const double non_local_weight = std::clamp(
            0.5 * request.wake_factor +
                0.5 * (1.0 - std::clamp(request.normalized_flow_alignment, -1.0, 1.0)),
            0.0, 1.0);
        result.electron_collection_scale =
            (1.0 - non_local_weight) *
                (request.electron_flux_moments.valid
                     ? 0.75 + 0.5 * request.electron_flux_moments.mean_cos_incidence
                     : 1.0) +
            non_local_weight * (0.60 + 0.4 * request.incidence_scale);
        result.ion_collection_scale =
            (1.0 - non_local_weight) *
                (request.ion_flux_moments.valid
                     ? 0.75 + 0.5 * request.ion_flux_moments.mean_cos_incidence
                     : 1.0) +
            non_local_weight * (0.70 + 0.3 * (1.0 - request.wake_factor));
        result.photo_incidence_scale =
            request.incidence_scale * (1.0 - 0.3 * non_local_weight);
        break;
    }
    case SurfaceDistributionModel::MultipleSurf:
    {
        const double population_factor =
            1.0 + 0.08 * std::max<std::ptrdiff_t>(
                              0, static_cast<std::ptrdiff_t>(request.electron_population_count +
                                                             request.ion_population_count) -
                                     2);
        result.electron_collection_scale =
            population_factor * request.local_field_factor *
            (request.electron_flux_moments.valid
                 ? 0.70 + 0.65 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            population_factor * (0.85 + 0.15 * request.incidence_scale) *
            (request.ion_flux_moments.valid
                 ? 0.72 + 0.55 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * (0.95 + 0.05 * request.local_field_factor);
        break;
    }
    case SurfaceDistributionModel::LocalModifiedPearsonIV:
    {
        const double energy_contrast = std::clamp(
            result.electron_characteristic_energy_ev /
                std::max(1.0e-6, result.ion_characteristic_energy_ev),
            0.25, 8.0);
        const double shape_factor =
            std::clamp(0.85 + 0.10 * request.local_field_factor +
                           0.05 * std::abs(request.local_reference_shift_v) /
                               std::max(1.0, result.electron_characteristic_energy_ev),
                       0.5, 2.0);
        result.electron_collection_scale =
            shape_factor * (0.65 + 0.35 * std::sqrt(energy_contrast)) *
            (request.electron_flux_moments.valid
                 ? 0.60 + 0.70 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            shape_factor * (0.70 + 0.25 / std::sqrt(energy_contrast)) *
            (request.ion_flux_moments.valid
                 ? 0.65 + 0.60 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(1.05 - 0.15 * request.wake_factor, 0.5, 1.2);
        break;
    }
    case SurfaceDistributionModel::LocalTabulated:
        result.electron_collection_scale =
            request.local_field_factor *
            (request.electron_flux_moments.valid
                 ? 0.70 + 0.60 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            (0.90 + 0.10 * request.incidence_scale) *
            (request.ion_flux_moments.valid
                 ? 0.70 + 0.55 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.9 + 0.1 * request.local_field_factor, 0.5, 1.5);
        break;
    case SurfaceDistributionModel::TwoAxesTabulatedVelocity:
        result.electron_collection_scale =
            request.local_field_factor *
            (request.electron_flux_moments.valid
                 ? 0.60 + 0.75 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            (0.85 + 0.15 * std::max(0.0, request.normalized_flow_alignment)) *
            (request.ion_flux_moments.valid
                 ? 0.65 + 0.70 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.85 + 0.2 * request.local_field_factor, 0.5, 1.5);
        break;
    case SurfaceDistributionModel::AxisymTabulatedVelocity:
    {
        const double anisotropy = std::clamp(request.axisym_anisotropy, 0.25, 4.0);
        result.electron_collection_scale =
            request.local_field_factor *
            std::clamp(0.90 + 0.12 * std::sqrt(anisotropy), 0.7, 1.3) *
            (request.electron_flux_moments.valid
                 ? 0.62 + 0.72 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            std::clamp(0.90 + 0.10 / std::sqrt(anisotropy), 0.7, 1.3) *
            (0.85 + 0.15 * std::max(0.0, request.normalized_flow_alignment)) *
            (request.ion_flux_moments.valid
                 ? 0.66 + 0.68 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.88 + 0.08 * request.local_field_factor, 0.5, 1.5);
        break;
    }
    case SurfaceDistributionModel::GlobalMaxwellBoltzmann:
        result.electron_collection_scale =
            (request.electron_flux_moments.valid
                 ? 0.78 + 0.42 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            (request.ion_flux_moments.valid
                 ? 0.82 + 0.38 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale = request.incidence_scale * 0.95;
        break;
    case SurfaceDistributionModel::GlobalMaxwellBoltzmann2:
        result.electron_collection_scale =
            (request.electron_flux_moments.valid
                 ? 0.74 + 0.50 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            (request.ion_flux_moments.valid
                 ? 0.80 + 0.42 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.electron_characteristic_energy_ev *= 1.12;
        result.ion_characteristic_energy_ev *= 1.08;
        result.photo_incidence_scale = request.incidence_scale * 0.92;
        break;
    case SurfaceDistributionModel::GlobalMaxwell:
        result.electron_collection_scale =
            (request.electron_flux_moments.valid
                 ? 0.80 + 0.36 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            (0.90 + 0.10 * std::max(0.0, request.normalized_flow_alignment)) *
            (request.ion_flux_moments.valid
                 ? 0.78 + 0.34 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale = request.incidence_scale;
        break;
    case SurfaceDistributionModel::LocalMaxwell:
        result.electron_collection_scale =
            request.local_field_factor *
            std::clamp(0.76 + 0.18 * std::abs(request.local_reference_shift_v) /
                                  std::max(1.0, result.electron_characteristic_energy_ev),
                       0.6, 1.5) *
            (request.electron_flux_moments.valid
                 ? 0.68 + 0.54 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            (request.ion_flux_moments.valid
                 ? 0.72 + 0.48 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.90 + 0.08 * request.local_field_factor, 0.5, 1.5);
        break;
    case SurfaceDistributionModel::LocalMaxwell2:
        result.electron_collection_scale =
            request.local_field_factor *
            std::clamp(0.82 + 0.12 * request.axisym_anisotropy +
                                  0.20 * std::abs(request.local_reference_shift_v) /
                                      std::max(1.0, result.electron_characteristic_energy_ev),
                       0.6, 1.8) *
            (request.electron_flux_moments.valid
                 ? 0.66 + 0.56 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            (0.88 + 0.12 * std::max(0.0, request.normalized_flow_alignment)) *
            (request.ion_flux_moments.valid
                 ? 0.70 + 0.50 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale *
            std::clamp((request.patch_role ? 1.05 : 0.95) + 0.05 * request.local_field_factor,
                       0.5, 1.6);
        break;
    case SurfaceDistributionModel::RecollMaxwell:
    {
        const double electron_recollection_ratio =
            clampRecollectionRatio(request.electron_recollection_ratio);
        const double ion_recollection_ratio =
            clampRecollectionRatio(request.ion_recollection_ratio);
        result.electron_collection_scale =
            electron_recollection_ratio * request.local_field_factor *
            std::clamp(0.76 + 0.18 * std::abs(request.local_reference_shift_v) /
                                  std::max(1.0, result.electron_characteristic_energy_ev),
                       0.6, 1.5) *
            (request.electron_flux_moments.valid
                 ? 0.68 + 0.54 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            ion_recollection_ratio * request.local_field_factor *
            (request.ion_flux_moments.valid
                 ? 0.72 + 0.48 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.88 + 0.06 * request.local_field_factor, 0.5, 1.4);
        break;
    }
    case SurfaceDistributionModel::PICSurf:
        result.electron_collection_scale =
            std::clamp(0.95 + 0.35 * request.local_field_factor + 0.20 * request.wake_factor,
                       0.2, 4.0) *
            (request.electron_flux_moments.valid
                 ? 0.65 + 0.70 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            std::clamp(0.90 + 0.25 * request.local_field_factor +
                           0.15 * std::max(0.0, request.normalized_flow_alignment),
                       0.2, 4.0) *
            (request.ion_flux_moments.valid
                 ? 0.75 + 0.40 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.75 + 0.25 * request.local_field_factor, 0.0, 1.8);
        break;
    case SurfaceDistributionModel::NonPICSurf:
        result.electron_collection_scale =
            std::clamp(0.85 + 0.10 * request.local_field_factor - 0.08 * request.wake_factor,
                       0.2, 3.0) *
            (request.electron_flux_moments.valid
                 ? 0.72 + 0.35 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            std::clamp(0.88 + 0.10 * std::max(0.0, request.normalized_flow_alignment), 0.2, 3.0) *
            (request.ion_flux_moments.valid
                 ? 0.80 + 0.20 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale = request.incidence_scale * 0.95;
        break;
    case SurfaceDistributionModel::GenericSurf:
        result.electron_collection_scale =
            std::clamp(0.90 + 0.15 * request.local_field_factor + 0.05 * request.wake_factor,
                       0.2, 3.5) *
            (request.electron_flux_moments.valid
                 ? 0.70 + 0.45 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            std::clamp(0.92 + 0.08 * std::abs(request.normalized_flow_alignment), 0.2, 3.5) *
            (request.ion_flux_moments.valid
                 ? 0.78 + 0.28 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.85 + 0.10 * request.local_field_factor, 0.0, 1.6);
        break;
    case SurfaceDistributionModel::GlobalSurf:
        result.electron_collection_scale =
            std::clamp(0.88 + 0.18 * std::max(0.0, request.normalized_flow_alignment), 0.2, 3.5) *
            (request.electron_flux_moments.valid
                 ? 0.76 + 0.30 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            std::clamp(0.94 + 0.12 * std::max(0.0, request.normalized_flow_alignment), 0.2, 3.5) *
            (request.ion_flux_moments.valid
                 ? 0.82 + 0.20 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale = request.incidence_scale;
        break;
    case SurfaceDistributionModel::LocalGenericSurf:
        result.electron_collection_scale =
            request.local_field_factor *
            std::clamp(0.82 + 0.25 * std::abs(request.local_reference_shift_v) /
                                  std::max(1.0, result.electron_characteristic_energy_ev),
                       0.4, 2.5) *
            (request.electron_flux_moments.valid
                 ? 0.72 + 0.42 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            std::clamp(0.86 + 0.15 * std::abs(request.local_reference_shift_v) /
                                  std::max(1.0, result.ion_characteristic_energy_ev),
                       0.4, 2.5) *
            (request.ion_flux_moments.valid
                 ? 0.80 + 0.22 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.88 + 0.12 * request.local_field_factor, 0.0, 1.8);
        break;
    case SurfaceDistributionModel::TestableSurf:
    {
        const double reference_shift_factor =
            std::clamp(1.0 / (1.0 + 0.03 * std::abs(request.local_reference_shift_v) /
                                        std::max(1.0, result.electron_characteristic_energy_ev)),
                       0.6, 1.2);
        result.electron_collection_scale =
            reference_shift_factor *
            std::clamp(0.94 + 0.08 * request.local_field_factor, 0.5, 2.0) *
            (request.electron_flux_moments.valid
                 ? 0.74 + 0.30 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            reference_shift_factor *
            std::clamp(0.96 + 0.06 * std::max(0.0, request.normalized_flow_alignment), 0.5, 2.0) *
            (request.ion_flux_moments.valid
                 ? 0.82 + 0.18 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.96 + 0.04 * request.local_field_factor, 0.0, 1.4);
        break;
    }
    case SurfaceDistributionModel::TestableForA:
    {
        const double localization_factor =
            std::clamp(0.88 + 0.20 * std::abs(request.local_reference_shift_v) /
                                  std::max(1.0, result.electron_characteristic_energy_ev),
                       0.6, 1.8);
        result.electron_collection_scale =
            request.local_field_factor * localization_factor *
            (request.electron_flux_moments.valid
                 ? 0.70 + 0.36 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            std::clamp(0.90 + 0.10 * request.incidence_scale, 0.4, 1.8) *
            (request.ion_flux_moments.valid
                 ? 0.78 + 0.22 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.92 + 0.06 * request.local_field_factor, 0.0, 1.6);
        break;
    }
    case SurfaceDistributionModel::MaxwellianThruster:
    {
        const double plume_alignment =
            std::clamp(0.5 + 0.5 * std::max(0.0, request.normalized_flow_alignment), 0.25, 1.0);
        const double thruster_bias =
            std::clamp(1.0 + 0.25 * plume_alignment + 0.10 * request.local_field_factor, 0.7, 2.5);
        result.electron_collection_scale =
            std::clamp(0.78 + 0.18 * request.local_field_factor - 0.08 * request.wake_factor,
                       0.25, 2.0) *
            (request.electron_flux_moments.valid
                 ? 0.66 + 0.42 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            thruster_bias * (request.patch_role ? 1.05 : 0.95) *
            (request.ion_flux_moments.valid
                 ? 0.85 + 0.50 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_characteristic_energy_ev *= 1.10;
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.70 + 0.20 * plume_alignment, 0.0, 1.2);
        break;
    }
    case SurfaceDistributionModel::UniformVelocity:
        result.electron_collection_scale =
            request.local_field_factor *
            (request.electron_flux_moments.valid
                 ? 0.75 + 0.55 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            request.local_field_factor *
            (request.ion_flux_moments.valid
                 ? 0.72 + 0.60 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale = request.incidence_scale;
        break;
    case SurfaceDistributionModel::Fluid:
        result.electron_collection_scale =
            (0.85 + 0.15 * request.local_field_factor) *
            (request.electron_flux_moments.valid
                 ? 0.65 + 0.60 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            (0.95 + 0.10 * std::max(0.0, request.normalized_flow_alignment)) *
            (request.ion_flux_moments.valid
                 ? 0.70 + 0.55 * request.ion_flux_moments.mean_cos_incidence
                 : 1.0);
        result.photo_incidence_scale =
            request.incidence_scale * std::clamp(0.9 + 0.1 * request.local_field_factor, 0.5, 1.5);
        break;
    case SurfaceDistributionModel::FowlerNordheim:
        result.electron_collection_scale =
            std::clamp(0.25 + 1.2 * request.local_field_factor, 0.1, 4.0);
        result.ion_collection_scale =
            std::clamp(0.25 + 0.5 * request.incidence_scale, 0.1, 4.0);
        result.photo_incidence_scale = std::clamp(0.2 * request.incidence_scale, 0.0, 2.0);
        break;
    case SurfaceDistributionModel::MultiPopulationHybrid:
    default:
    {
        const double electron_population_factor =
            1.0 + 0.12 * std::max<std::ptrdiff_t>(
                              0, static_cast<std::ptrdiff_t>(request.electron_population_count) - 1);
        const double ion_population_factor =
            1.0 + 0.10 * std::max<std::ptrdiff_t>(
                              0, static_cast<std::ptrdiff_t>(request.ion_population_count) - 1);
        const double reference_shift_factor =
            1.0 / (1.0 + 0.05 * std::abs(request.local_reference_shift_v) /
                             std::max(1.0, result.electron_characteristic_energy_ev));
        result.electron_collection_scale =
            electron_population_factor * reference_shift_factor * request.local_field_factor *
            (request.electron_flux_moments.valid
                 ? 0.75 + 0.5 * request.electron_flux_moments.mean_cos_incidence
                 : 1.0);
        result.ion_collection_scale =
            ion_population_factor * (0.9 + 0.1 * std::max(0.0, request.normalized_flow_alignment)) *
            request.local_field_factor *
            (request.ion_flux_moments.valid ? 0.75 + 0.5 * request.ion_flux_moments.mean_cos_incidence
                                            : 1.0);
        result.photo_incidence_scale = request.incidence_scale;
        break;
    }
    }

    result.electron_collection_scale = std::clamp(result.electron_collection_scale, 0.1, 4.0);
    result.ion_collection_scale = std::clamp(result.ion_collection_scale, 0.1, 4.0);
    result.photo_incidence_scale = std::clamp(result.photo_incidence_scale, 0.0, 2.0);
    result.valid = true;
    return result;
}

SurfaceDistributionBuildResult buildSurfaceFluxMoments(
    const SurfaceDistributionBuildRequest& request)
{
    SurfaceDistributionBuildResult result;
    std::unique_ptr<SurfaceDistributionFunction> electron_distribution;
    std::unique_ptr<SurfaceDistributionFunction> ion_distribution;

    switch (request.model)
    {
    case SurfaceDistributionModel::MaxwellianProjected:
        electron_distribution = std::make_unique<IsotropicMaxwellianFluxDistribution>(
            request.electron_density_m3, request.electron_temperature_ev);
        ion_distribution = std::make_unique<IsotropicMaxwellianFluxDistribution>(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s);
        break;
    case SurfaceDistributionModel::WakeAnisotropic:
        electron_distribution = std::make_unique<TabulatedSurfaceDistributionFunction>(
            std::vector<double>{std::max(1.0e-3, 0.5 * request.electron_temperature_ev),
                                std::max(1.0e-3, request.electron_temperature_ev),
                                std::max(1.0e-3, 2.0 * request.electron_temperature_ev)},
            std::vector<double>{1.0, 1.5, 0.8}, std::vector<double>{0.15, 0.35, 0.65});
        ion_distribution = std::make_unique<TabulatedSurfaceDistributionFunction>(
            std::vector<double>{std::max(1.0e-3, 0.5 * request.ion_temperature_ev),
                                std::max(1.0e-3, request.ion_temperature_ev),
                                std::max(1.0e-3, 2.0 * request.ion_temperature_ev)},
            std::vector<double>{0.6, 1.2, 1.8}, std::vector<double>{0.65, 0.8, 0.95});
        break;
    case SurfaceDistributionModel::TabulatedVelocity:
    {
        const double electron_velocity_scale =
            std::max(1.0e3, 1.0e5 * std::sqrt(std::max(1.0e-6, request.electron_temperature_ev)));
        const double ion_velocity_scale = std::max(
            1.0e3,
            std::abs(request.ion_directed_velocity_m_per_s) +
                2.5e4 * std::sqrt(std::max(1.0e-6, request.ion_temperature_ev)));

        electron_distribution = std::make_unique<TabulatedVelocitySurfaceDistributionFunction>(
            std::vector<double>{0.5 * electron_velocity_scale, electron_velocity_scale,
                                1.5 * electron_velocity_scale, 2.0 * electron_velocity_scale},
            std::vector<double>{0.5, 1.0, 0.7, 0.3}, 1.0e-10,
            std::vector<double>{0.20, 0.40, 0.70, 0.90});
        ion_distribution = std::make_unique<TabulatedVelocitySurfaceDistributionFunction>(
            std::vector<double>{0.5 * ion_velocity_scale, ion_velocity_scale,
                                1.5 * ion_velocity_scale, 2.0 * ion_velocity_scale},
            std::vector<double>{0.4, 0.8, 1.2, 0.6}, 2.5e-11,
            std::vector<double>{0.60, 0.75, 0.90, 0.98});
        break;
    }
    case SurfaceDistributionModel::NonLocalizedHybrid:
    {
        auto local_electron_distribution =
            IsotropicMaxwellianFluxDistribution(request.electron_density_m3,
                                                request.electron_temperature_ev);
        auto local_ion_distribution = IsotropicMaxwellianFluxDistribution(
            request.ion_density_m3, request.ion_temperature_ev,
            request.ion_directed_velocity_m_per_s);
        auto non_local_electron_distribution = TabulatedSurfaceDistributionFunction(
            std::vector<double>{std::max(1.0e-3, 0.5 * request.electron_temperature_ev),
                                std::max(1.0e-3, request.electron_temperature_ev),
                                std::max(1.0e-3, 2.5 * request.electron_temperature_ev)},
            std::vector<double>{0.7, 1.1, 0.9}, std::vector<double>{0.25, 0.50, 0.80});
        auto non_local_ion_distribution = TabulatedSurfaceDistributionFunction(
            std::vector<double>{std::max(1.0e-3, 0.4 * request.ion_temperature_ev),
                                std::max(1.0e-3, request.ion_temperature_ev),
                                std::max(1.0e-3, 3.0 * request.ion_temperature_ev)},
            std::vector<double>{0.5, 1.0, 1.5}, std::vector<double>{0.55, 0.75, 0.95});

        const double non_local_weight =
            std::clamp(std::abs(request.ion_directed_velocity_m_per_s) / 1.0e4, 0.0, 1.0) * 0.6;
        electron_distribution =
            std::make_unique<NonLocalizedHybridSurfaceDistributionFunction>(
                local_electron_distribution.computeMoments(),
                non_local_electron_distribution.computeMoments(), non_local_weight);
        ion_distribution = std::make_unique<NonLocalizedHybridSurfaceDistributionFunction>(
            local_ion_distribution.computeMoments(),
            non_local_ion_distribution.computeMoments(), non_local_weight);
        break;
    }
    case SurfaceDistributionModel::MultipleSurf:
    {
        std::vector<SurfaceFluxMoments> electron_components;
        std::vector<double> electron_weights;
        if (!request.electron_populations.empty())
        {
            for (std::size_t i = 0; i < request.electron_populations.size(); ++i)
            {
                const auto& population = request.electron_populations[i];
                electron_components.push_back(
                    IsotropicMaxwellianFluxDistribution(population.density_m3,
                                                        population.temperature_ev)
                        .computeMoments());
                electron_weights.push_back(1.0 /
                                           static_cast<double>(i + 1));
            }
        }
        else
        {
            electron_components.push_back(
                IsotropicMaxwellianFluxDistribution(request.electron_density_m3,
                                                    request.electron_temperature_ev)
                    .computeMoments());
            electron_components.push_back(
                TabulatedSurfaceDistributionFunction(
                    std::vector<double>{std::max(1.0e-3, 0.6 * request.electron_temperature_ev),
                                        std::max(1.0e-3, 1.4 * request.electron_temperature_ev),
                                        std::max(1.0e-3, 2.8 * request.electron_temperature_ev)},
                    std::vector<double>{1.0, 0.8, 0.35},
                    std::vector<double>{0.35, 0.60, 0.82})
                    .computeMoments());
            electron_weights = {1.0, 0.65};
        }

        std::vector<SurfaceFluxMoments> ion_components;
        std::vector<double> ion_weights;
        if (!request.ion_populations.empty())
        {
            for (std::size_t i = 0; i < request.ion_populations.size(); ++i)
            {
                const auto& population = request.ion_populations[i];
                ion_components.push_back(
                    IsotropicMaxwellianFluxDistribution(
                        population.density_m3, population.temperature_ev,
                        i == 0 ? request.ion_directed_velocity_m_per_s : 0.0)
                        .computeMoments());
                ion_weights.push_back(1.0 + 0.2 * static_cast<double>(i));
            }
        }
        else
        {
            ion_components.push_back(
                IsotropicMaxwellianFluxDistribution(
                    request.ion_density_m3, request.ion_temperature_ev,
                    request.ion_directed_velocity_m_per_s)
                    .computeMoments());
            ion_components.push_back(
                TabulatedSurfaceDistributionFunction(
                    std::vector<double>{std::max(1.0e-3, 0.5 * request.ion_temperature_ev),
                                        std::max(1.0e-3, 1.2 * request.ion_temperature_ev),
                                        std::max(1.0e-3, 2.6 * request.ion_temperature_ev)},
                    std::vector<double>{0.45, 0.85, 1.10},
                    std::vector<double>{0.60, 0.78, 0.93})
                    .computeMoments());
            ion_weights = {0.85, 1.0};
        }

        electron_distribution = std::make_unique<MultipleSurfSurfaceDistributionFunction>(
            std::move(electron_components), std::move(electron_weights), -0.2);
        ion_distribution = std::make_unique<MultipleSurfSurfaceDistributionFunction>(
            std::move(ion_components), std::move(ion_weights), 0.3);
        break;
    }
    case SurfaceDistributionModel::LocalModifiedPearsonIV:
    {
        const double electron_skewness = std::clamp(
            0.2 + 0.6 * std::tanh(request.electron_temperature_ev / 20.0), -2.5, 2.5);
        const double ion_skewness = std::clamp(
            std::abs(request.ion_directed_velocity_m_per_s) / 5000.0 - 0.35, -2.5, 2.5);
        const double electron_kurtosis =
            std::clamp(3.5 + 0.8 * std::max(0.0, request.electron_temperature_ev / 10.0 - 1.0),
                       2.0, 8.0);
        const double ion_kurtosis =
            std::clamp(3.2 + 0.5 * std::abs(request.ion_directed_velocity_m_per_s) / 3000.0,
                       2.0, 8.0);
        electron_distribution =
            std::make_unique<LocalModifiedPearsonIVSurfaceDistributionFunction>(
                request.electron_density_m3, 2.0 * request.electron_temperature_ev,
                electron_skewness, electron_kurtosis, -0.15, 1.05);
        ion_distribution =
            std::make_unique<LocalModifiedPearsonIVSurfaceDistributionFunction>(
                request.ion_density_m3, 2.0 * request.ion_temperature_ev, ion_skewness,
                ion_kurtosis, 0.25, 1.0 + 0.2 * std::clamp(
                                               std::abs(request.ion_directed_velocity_m_per_s) /
                                                   5000.0,
                                               0.0, 1.0));
        break;
    }
    case SurfaceDistributionModel::LocalTabulated:
    {
        electron_distribution = std::make_unique<LocalTabulatedSurfaceDistributionFunction>(
            std::vector<double>{std::max(1.0e-3, 0.45 * request.electron_temperature_ev),
                                std::max(1.0e-3, request.electron_temperature_ev),
                                std::max(1.0e-3, 2.2 * request.electron_temperature_ev),
                                std::max(1.0e-3, 3.8 * request.electron_temperature_ev)},
            std::vector<double>{0.8, 1.4, 0.9, 0.35},
            0.1 * request.electron_temperature_ev, 1.0,
            std::vector<double>{0.25, 0.45, 0.68, 0.88});
        ion_distribution = std::make_unique<LocalTabulatedSurfaceDistributionFunction>(
            std::vector<double>{std::max(1.0e-3, 0.35 * request.ion_temperature_ev),
                                std::max(1.0e-3, request.ion_temperature_ev),
                                std::max(1.0e-3, 1.8 * request.ion_temperature_ev),
                                std::max(1.0e-3, 3.2 * request.ion_temperature_ev)},
            std::vector<double>{0.5, 0.9, 1.1, 0.7},
            0.05 * request.ion_temperature_ev,
            1.0 + 0.1 * std::clamp(std::abs(request.ion_directed_velocity_m_per_s) / 4000.0, 0.0, 1.0),
            std::vector<double>{0.55, 0.72, 0.86, 0.95});
        break;
    }
    case SurfaceDistributionModel::TwoAxesTabulatedVelocity:
    {
        const double electron_normal_scale =
            std::max(1.0e3, 8.5e4 * std::sqrt(std::max(1.0e-6, request.electron_temperature_ev)));
        const double electron_tangential_scale = 0.65 * electron_normal_scale;
        const double ion_normal_scale = std::max(
            1.0e3,
            std::abs(request.ion_directed_velocity_m_per_s) +
                2.0e4 * std::sqrt(std::max(1.0e-6, request.ion_temperature_ev)));
        const double ion_tangential_scale = 0.55 * ion_normal_scale;

        electron_distribution =
            std::make_unique<TwoAxesTabulatedVelocitySurfaceDistributionFunction>(
                std::vector<double>{0.5 * electron_normal_scale, 1.0 * electron_normal_scale,
                                    1.6 * electron_normal_scale, 2.1 * electron_normal_scale},
                std::vector<double>{0.2 * electron_tangential_scale, 0.7 * electron_tangential_scale,
                                    1.1 * electron_tangential_scale, 1.4 * electron_tangential_scale},
                std::vector<double>{0.55, 1.0, 0.75, 0.28}, 1.0e-10, 6.0e-11,
                std::vector<double>{0.28, 0.46, 0.72, 0.9});
        ion_distribution =
            std::make_unique<TwoAxesTabulatedVelocitySurfaceDistributionFunction>(
                std::vector<double>{0.55 * ion_normal_scale, 1.0 * ion_normal_scale,
                                    1.45 * ion_normal_scale, 1.9 * ion_normal_scale},
                std::vector<double>{0.25 * ion_tangential_scale, 0.5 * ion_tangential_scale,
                                    0.85 * ion_tangential_scale, 1.15 * ion_tangential_scale},
                std::vector<double>{0.4, 0.85, 1.15, 0.65}, 2.5e-11, 1.5e-11,
                std::vector<double>{0.58, 0.74, 0.88, 0.96});
        break;
    }
    case SurfaceDistributionModel::AxisymTabulatedVelocity:
    {
        const double anisotropy = std::clamp(request.axisym_anisotropy, 0.25, 4.0);
        const double electron_velocity_scale =
            std::max(1.0e3, 8.0e4 * std::sqrt(std::max(1.0e-6, request.electron_temperature_ev)));
        const double ion_velocity_scale = std::max(
            1.0e3,
            std::abs(request.ion_directed_velocity_m_per_s) +
                1.8e4 * std::sqrt(std::max(1.0e-6, request.ion_temperature_ev)));
        const double axial_scale = std::sqrt(anisotropy);
        const double radial_scale = 1.0 / axial_scale;

        electron_distribution =
            std::make_unique<AxisymTabulatedVelocitySurfaceDistributionFunction>(
                std::vector<double>{0.6 * electron_velocity_scale * axial_scale,
                                    1.0 * electron_velocity_scale * axial_scale,
                                    1.5 * electron_velocity_scale * axial_scale,
                                    2.0 * electron_velocity_scale * axial_scale},
                std::vector<double>{0.35 * electron_velocity_scale * radial_scale,
                                    0.65 * electron_velocity_scale * radial_scale,
                                    0.95 * electron_velocity_scale * radial_scale,
                                    1.20 * electron_velocity_scale * radial_scale},
                std::vector<double>{0.55, 1.0, 0.75, 0.32}, 1.0e-10, 7.0e-11,
                std::vector<double>{0.30, 0.48, 0.74, 0.90});
        ion_distribution = std::make_unique<AxisymTabulatedVelocitySurfaceDistributionFunction>(
            std::vector<double>{0.6 * ion_velocity_scale * axial_scale,
                                1.0 * ion_velocity_scale * axial_scale,
                                1.45 * ion_velocity_scale * axial_scale,
                                1.9 * ion_velocity_scale * axial_scale},
            std::vector<double>{0.30 * ion_velocity_scale * radial_scale,
                                0.55 * ion_velocity_scale * radial_scale,
                                0.80 * ion_velocity_scale * radial_scale,
                                1.05 * ion_velocity_scale * radial_scale},
            std::vector<double>{0.42, 0.86, 1.10, 0.62}, 2.5e-11, 1.6e-11,
            std::vector<double>{0.58, 0.74, 0.88, 0.96});
        break;
    }
    case SurfaceDistributionModel::GlobalMaxwellBoltzmann:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.15, 1.0, 0.55,
            0.42, 0.18);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            2.05, 1.0, 0.40, 0.55, 0.22);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::GlobalMaxwellBoltzmann2:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.45, 1.08, 0.85,
            0.40, 0.20);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            2.25, 1.05, 0.70, 0.52, 0.24);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::GlobalMaxwell:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.0, 0.96, 0.25,
            0.48, 0.16);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            2.0, 0.94, 0.15, 0.60, 0.20);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::LocalMaxwell:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.05, 1.02, 0.35,
            0.50, 0.14);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            1.95, 1.00, 0.20, 0.62, 0.16);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::LocalMaxwell2:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.20, 1.06, 0.55,
            0.46, 0.18);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            2.05, 1.03, 0.35, 0.60, 0.22);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::RecollMaxwell:
        result.electron_flux_moments = scaleFluxMomentsForRecollection(
            computeParametricMaxwellMoments(request.electron_density_m3,
                                            request.electron_temperature_ev,
                                            0.0,
                                            2.05,
                                            1.02,
                                            0.35,
                                            0.50,
                                            0.14),
            request.electron_recollection_ratio);
        result.ion_flux_moments = scaleFluxMomentsForRecollection(
            computeParametricMaxwellMoments(request.ion_density_m3,
                                            request.ion_temperature_ev,
                                            request.ion_directed_velocity_m_per_s,
                                            1.95,
                                            1.00,
                                            0.20,
                                            0.62,
                                            0.16),
            request.ion_recollection_ratio);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::PICSurf:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.30, 1.12, 0.65,
            0.46, 0.22);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            2.15, 1.10, 0.50, 0.60, 0.26);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::NonPICSurf:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 1.90, 0.92, 0.10,
            0.52, 0.12);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            1.90, 0.90, 0.05, 0.64, 0.14);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::GenericSurf:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.00, 0.98, 0.25,
            0.50, 0.15);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            1.95, 0.97, 0.18, 0.62, 0.16);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::GlobalSurf:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.05, 1.00, 0.20,
            0.54, 0.14);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            2.00, 0.99, 0.12, 0.66, 0.18);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::LocalGenericSurf:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.08, 1.00, 0.30,
            0.50, 0.16);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s,
            1.98, 1.00, 0.22, 0.64, 0.18);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::TestableSurf:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.02, 0.98, 0.18,
            0.52, 0.12);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev,
            request.ion_directed_velocity_m_per_s, 1.98, 0.98, 0.14, 0.64, 0.14);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::TestableForA:
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.08, 1.01, 0.28,
            0.50, 0.15);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev,
            request.ion_directed_velocity_m_per_s, 2.00, 1.00, 0.20, 0.62, 0.16);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    case SurfaceDistributionModel::MaxwellianThruster:
    {
        const double ion_drift_speed_m_per_s =
            std::max(std::abs(request.ion_directed_velocity_m_per_s), 3500.0);
        result.electron_flux_moments = computeParametricMaxwellMoments(
            request.electron_density_m3, request.electron_temperature_ev, 0.0, 2.05, 1.00, 0.25,
            0.48, 0.16);
        result.ion_flux_moments = computeParametricMaxwellMoments(
            request.ion_density_m3, request.ion_temperature_ev, ion_drift_speed_m_per_s, 2.35,
            1.20, 0.85, 0.78, 0.30);
        result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
        return result;
    }
    case SurfaceDistributionModel::UniformVelocity:
        electron_distribution = std::make_unique<UniformVelocitySurfaceDistributionFunction>(
            std::max(1.0e3, 7.5e4 * std::sqrt(std::max(1.0e-6, request.electron_temperature_ev))),
            std::max(1.0, request.electron_density_m3 * 1.0e-11), 1.0e-10, 0.7);
        ion_distribution = std::make_unique<UniformVelocitySurfaceDistributionFunction>(
            std::max(1.0e3,
                     std::abs(request.ion_directed_velocity_m_per_s) +
                         1.5e4 * std::sqrt(std::max(1.0e-6, request.ion_temperature_ev))),
            std::max(1.0, request.ion_density_m3 * 1.0e-11), 2.5e-11, 0.82);
        break;
    case SurfaceDistributionModel::Fluid:
        electron_distribution = std::make_unique<FluidSurfaceDistributionFunction>(
            request.electron_density_m3, 2.0 * request.electron_temperature_ev, 0.0, 1.05, -0.1);
        ion_distribution = std::make_unique<FluidSurfaceDistributionFunction>(
            request.ion_density_m3, 2.0 * request.ion_temperature_ev,
            request.ion_directed_velocity_m_per_s,
            1.0 + 0.2 * std::clamp(std::abs(request.ion_directed_velocity_m_per_s) / 4000.0, 0.0, 1.0),
            0.25);
        break;
    case SurfaceDistributionModel::FowlerNordheim:
        electron_distribution = std::make_unique<FowlerNordheimSurfaceDistributionFunction>(
            std::max(1.0e6, 5.0e6 + 1.0e5 * request.electron_temperature_ev), 4.5, 1.0);
        ion_distribution = std::make_unique<UniformVelocitySurfaceDistributionFunction>(
            std::max(1.0e3,
                     std::abs(request.ion_directed_velocity_m_per_s) +
                         1.0e4 * std::sqrt(std::max(1.0e-6, request.ion_temperature_ev))),
            std::max(1.0, request.ion_density_m3 * 5.0e-12), 2.5e-11, 0.65);
        break;
    case SurfaceDistributionModel::MultiPopulationHybrid:
    default:
        if (request.electron_populations.size() >= 2)
        {
            electron_distribution = std::make_unique<BiMaxwellianFluxDistribution>(
                request.electron_populations[0].density_m3,
                request.electron_populations[0].temperature_ev,
                request.electron_populations[1].density_m3,
                request.electron_populations[1].temperature_ev);
        }
        else
        {
            electron_distribution = std::make_unique<IsotropicMaxwellianFluxDistribution>(
                request.electron_density_m3, request.electron_temperature_ev);
        }
        if (request.ion_populations.size() >= 2)
        {
            ion_distribution = std::make_unique<BiMaxwellianFluxDistribution>(
                request.ion_populations[0].density_m3,
                request.ion_populations[0].temperature_ev,
                request.ion_populations[1].density_m3,
                request.ion_populations[1].temperature_ev);
        }
        else
        {
            ion_distribution = std::make_unique<IsotropicMaxwellianFluxDistribution>(
                request.ion_density_m3, request.ion_temperature_ev, request.ion_directed_velocity_m_per_s);
        }
        break;
    }

    result.electron_flux_moments =
        electron_distribution ? electron_distribution->computeMoments() : SurfaceFluxMoments{};
    result.ion_flux_moments =
        ion_distribution ? ion_distribution->computeMoments() : SurfaceFluxMoments{};
    result.valid = result.electron_flux_moments.valid || result.ion_flux_moments.valid;
    return result;
}

SurfaceDistributionBundleResult evaluateSurfaceDistributionBundle(
    const SurfaceDistributionBundleRequest& request)
{
    SurfaceDistributionBundleResult result;
    result.built_flux = buildSurfaceFluxMoments(request.build_request);

    auto synthesis_request = request.synthesis_request;
    synthesis_request.electron_flux_moments = result.built_flux.electron_flux_moments;
    synthesis_request.ion_flux_moments = result.built_flux.ion_flux_moments;
    result.synthesized = synthesizeSurfaceDistributionMoments(synthesis_request);
    result.valid = result.built_flux.valid && result.synthesized.valid;
    return result;
}

SurfaceDistributionBundleRequest makeSurfaceDistributionBundleRequest(
    const SurfaceDistributionEnvironment& environment)
{
    SurfaceDistributionBundleRequest request;
    request.build_request.model = environment.model;
    request.build_request.electron_density_m3 = environment.electron_density_m3;
    request.build_request.electron_temperature_ev = environment.electron_temperature_ev;
    request.build_request.ion_density_m3 = environment.ion_density_m3;
    request.build_request.ion_temperature_ev = environment.ion_temperature_ev;
    request.build_request.ion_directed_velocity_m_per_s =
        environment.ion_directed_velocity_m_per_s;
    request.build_request.axisym_anisotropy = environment.axisym_anisotropy;
    request.build_request.electron_recollection_ratio = environment.electron_recollection_ratio;
    request.build_request.ion_recollection_ratio = environment.ion_recollection_ratio;
    request.build_request.electron_populations = environment.electron_populations;
    request.build_request.ion_populations = environment.ion_populations;

    request.synthesis_request.model = environment.model;
    request.synthesis_request.fallback_electron_energy_ev =
        environment.fallback_electron_energy_ev;
    request.synthesis_request.fallback_ion_energy_ev = environment.fallback_ion_energy_ev;
    request.synthesis_request.local_reference_shift_v = environment.local_reference_shift_v;
    request.synthesis_request.incidence_scale = environment.incidence_scale;
    request.synthesis_request.normalized_flow_alignment =
        environment.normalized_flow_alignment;
    request.synthesis_request.local_field_factor = environment.local_field_factor;
    request.synthesis_request.wake_factor = environment.wake_factor;
    request.synthesis_request.axisym_anisotropy = environment.axisym_anisotropy;
    request.synthesis_request.electron_recollection_ratio =
        environment.electron_recollection_ratio;
    request.synthesis_request.ion_recollection_ratio = environment.ion_recollection_ratio;
    request.synthesis_request.patch_role = environment.patch_role;
    request.synthesis_request.electron_population_count =
        environment.electron_populations.size();
    request.synthesis_request.ion_population_count = environment.ion_populations.size();
    return request;
}

SurfaceDistributionBundleResult evaluateSurfaceDistributionEnvironment(
    const SurfaceDistributionEnvironment& environment)
{
    return evaluateSurfaceDistributionBundle(makeSurfaceDistributionBundleRequest(environment));
}

SurfaceDistributionEnvironment makeSurfaceDistributionEnvironment(
    const SurfaceDistributionRoleInputs& inputs)
{
    SurfaceDistributionEnvironment environment;
    environment.model = inputs.model;
    environment.electron_density_m3 = inputs.electron_density_m3;
    environment.electron_temperature_ev = inputs.electron_temperature_ev;
    environment.ion_density_m3 = inputs.ion_density_m3;
    environment.ion_temperature_ev = inputs.ion_temperature_ev;
    environment.ion_directed_velocity_m_per_s = inputs.ion_directed_velocity_m_per_s;
    environment.electron_populations = inputs.electron_populations;
    environment.ion_populations = inputs.ion_populations;
    environment.fallback_electron_energy_ev = inputs.fallback_electron_energy_ev;
    environment.fallback_ion_energy_ev = inputs.fallback_ion_energy_ev;
    environment.local_reference_shift_v = inputs.local_reference_shift_v;
    environment.incidence_scale = inputs.share_patch_distribution
                                      ? 1.0
                                      : std::max(0.0, std::cos(inputs.patch_incidence_angle_deg *
                                                               3.14159265358979323846 / 180.0));
    environment.normalized_flow_alignment =
        inputs.bulk_flow_velocity_m_per_s > 0.0
            ? std::clamp(inputs.projected_speed_m_per_s /
                             std::max(1.0, std::abs(inputs.bulk_flow_velocity_m_per_s)),
                         -1.0, 1.0)
            : std::clamp(inputs.flow_alignment_cosine, -1.0, 1.0);
    environment.local_field_factor =
        inputs.share_patch_distribution
            ? 1.0
            : (1.0 + 0.05 *
                         std::clamp(std::abs(inputs.normal_electric_field_v_per_m) / 1.0e4, 0.0,
                                    4.0));
    environment.wake_factor = 0.5 * (1.0 - environment.normalized_flow_alignment);
    environment.axisym_anisotropy = std::clamp(inputs.axisym_anisotropy, 0.25, 4.0);
    environment.electron_recollection_ratio =
        clampRecollectionRatio(inputs.electron_recollection_ratio);
    environment.ion_recollection_ratio = clampRecollectionRatio(inputs.ion_recollection_ratio);
    environment.patch_role = inputs.patch_role;
    return environment;
}

SurfaceDistributionRoleInputs makeSurfaceDistributionRoleInputs(
    const SurfaceDistributionRoleBuildRequest& request)
{
    SurfaceDistributionRoleInputs inputs;
    inputs.model = request.model;
    inputs.electron_density_m3 = request.electron_density_m3;
    inputs.electron_temperature_ev = request.electron_temperature_ev;
    inputs.ion_density_m3 = request.ion_density_m3;
    inputs.ion_temperature_ev = request.ion_temperature_ev;
    inputs.ion_directed_velocity_m_per_s = request.ion_directed_velocity_m_per_s;
    for (const auto& population : request.electron_spectrum.populations)
    {
        inputs.electron_populations.push_back(
            {population.density_m3, population.temperature_ev});
    }
    for (const auto& population : request.ion_spectrum.populations)
    {
        inputs.ion_populations.push_back(
            {population.density_m3, population.temperature_ev});
    }
    inputs.fallback_electron_energy_ev = std::max(
        request.electron_temperature_ev,
        resolvedSpectrumCharacteristicEnergyEv(request.electron_spectrum,
                                              request.electron_temperature_ev));
    inputs.fallback_ion_energy_ev = std::max(
        request.ion_temperature_ev,
        resolvedSpectrumCharacteristicEnergyEv(request.ion_spectrum,
                                              request.ion_temperature_ev));
    inputs.local_reference_shift_v = request.local_reference_shift_v;
    inputs.patch_incidence_angle_deg = request.patch_incidence_angle_deg;
    inputs.bulk_flow_velocity_m_per_s = request.bulk_flow_velocity_m_per_s;
    inputs.projected_speed_m_per_s = request.projected_speed_m_per_s;
    inputs.flow_alignment_cosine = request.flow_alignment_cosine;
    inputs.normal_electric_field_v_per_m = request.normal_electric_field_v_per_m;
    inputs.axisym_anisotropy = request.axisym_anisotropy;
    inputs.electron_recollection_ratio = request.electron_recollection_ratio;
    inputs.ion_recollection_ratio = request.ion_recollection_ratio;
    inputs.share_patch_distribution = request.share_patch_distribution;
    inputs.patch_role = request.patch_role;
    return inputs;
}

double interpolateLegacySpectrum(const std::vector<double>& x,
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
        const double dx = x1 - x0;
        if (std::abs(dx) <= 1.0e-20)
        {
            return y[0];
        }
        const double alpha = (target - x0) / dx;
        return y[0] + alpha * (y[1] - y[0]);
    }

    std::ptrdiff_t k = 0;
    std::ptrdiff_t m = 2;
    if (target >= x[count - 2])
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

std::vector<double> smoothLegacySpectrumLine3(const std::vector<double>& x,
                                              const std::vector<double>& y)
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
        const std::size_t last = std::min(count - 1, first + 2);
        const double x0 = x[first];
        const double x1 = x[first + 1];
        const double x2 = x[last];
        const double y0 = y[first];
        const double y1 = y[first + 1];
        const double y2 = y[last];
        const double t0 = 3.0;
        const double t1 = x0 + x1 + x2;
        const double t2 = x0 * x0 + x1 * x1 + x2 * x2;
        const double p1 = y0 + y1 + y2;
        const double p2 = y0 * x0 + y1 * x1 + y2 * x2;
        const double det = t0 * t2 - t1 * t1;
        if (std::abs(det) <= 1.0e-20)
        {
            smoothed[i] = y[i];
            continue;
        }
        const double a0 = (p1 * t2 - p2 * t1) / det;
        const double a1 = (t0 * p2 - t1 * p1) / det;
        smoothed[i] = a0 + a1 * x[i];
    }
    return smoothed;
}

double legacyElectronComponentFlux(double density_m3,
                                   double temperature_ev,
                                   double energy_ev)
{
    const double density_cm3 = std::max(0.0, density_m3) * 1.0e-6;
    const double temperature = std::max(1.0e-6, temperature_ev);
    return 5.325e10 * density_cm3 * energy_ev * std::pow(temperature, -1.5) *
           std::exp(-energy_ev / temperature);
}

double legacyIonComponentFlux(double density_m3,
                              double temperature_ev,
                              double mass_amu,
                              double energy_ev)
{
    const double density_cm3 = std::max(0.0, density_m3) * 1.0e-6;
    const double temperature = std::max(1.0e-6, temperature_ev);
    return 1.244e9 * density_cm3 * energy_ev * std::pow(temperature, -1.5) *
           std::exp(-energy_ev / temperature) / std::sqrt(std::max(1.0, mass_amu));
}

LegacyFluxTables buildLegacyFluxTables(const LegacyFluxTableBuildRequest& request)
{
    const int sample_count = std::max(2, request.sample_count);

    LegacyFluxTables tables;
    tables.electron_populations = collectLegacyPopulations(
        request.electron_spectrum, request.has_electron_spectrum,
        request.fallback_electron_density_m3, request.fallback_electron_temperature_ev,
        request.fallback_electron_mass_amu, request.max_population_count);
    tables.ion_populations = collectLegacyPopulations(
        request.ion_spectrum, request.has_ion_spectrum, request.fallback_ion_density_m3,
        request.fallback_ion_temperature_ev, request.fallback_ion_mass_amu,
        request.max_population_count);

    double spectrum_min_ev = std::numeric_limits<double>::max();
    double spectrum_max_ev = 0.0;
    auto extend_spectrum_range = [&](const std::vector<double>& energies) {
        for (double energy_ev : energies)
        {
            if (energy_ev > 0.0)
            {
                spectrum_min_ev = std::min(spectrum_min_ev, energy_ev);
                spectrum_max_ev = std::max(spectrum_max_ev, energy_ev);
            }
        }
    };
    extend_spectrum_range(request.electron_spectrum.energy_grid_ev);
    extend_spectrum_range(request.ion_spectrum.energy_grid_ev);

    double min_thermal_ev = std::numeric_limits<double>::max();
    double max_thermal_ev = 0.0;
    auto extend_thermal_range = [&](const std::vector<LegacyPopulation>& populations) {
        for (const auto& population : populations)
        {
            if (population.temperature_ev > 0.0)
            {
                min_thermal_ev = std::min(min_thermal_ev, population.temperature_ev);
                max_thermal_ev = std::max(max_thermal_ev, population.temperature_ev);
            }
        }
    };
    extend_thermal_range(tables.electron_populations);
    extend_thermal_range(tables.ion_populations);

    double minimum_energy_ev = std::isfinite(spectrum_min_ev) ? spectrum_min_ev : 0.0;
    double maximum_energy_ev = spectrum_max_ev;
    if (std::isfinite(min_thermal_ev))
    {
        minimum_energy_ev = minimum_energy_ev > 0.0
                                ? std::min(minimum_energy_ev, 0.05 * min_thermal_ev)
                                : 0.05 * min_thermal_ev;
    }
    if (max_thermal_ev > 0.0)
    {
        maximum_energy_ev = std::max(maximum_energy_ev, 10.0 * max_thermal_ev);
    }
    minimum_energy_ev = std::max(1.0e-4, minimum_energy_ev > 0.0 ? minimum_energy_ev : 1.0e-2);
    maximum_energy_ev = std::max(10.0 * minimum_energy_ev, maximum_energy_ev);

    const double nmin = std::log10(minimum_energy_ev);
    const double nmax = std::log10(maximum_energy_ev);
    const double dn = (nmax - nmin) / static_cast<double>(sample_count);

    tables.energy_ev.reserve(static_cast<std::size_t>(sample_count));
    tables.width_ev.reserve(static_cast<std::size_t>(sample_count));
    tables.total_electron_flux.reserve(static_cast<std::size_t>(sample_count));
    tables.spectrum_electron_flux.reserve(static_cast<std::size_t>(sample_count));
    tables.total_ion_flux.reserve(static_cast<std::size_t>(sample_count));
    tables.spectrum_ion_flux.reserve(static_cast<std::size_t>(sample_count));
    tables.electron_component_flux.assign(tables.electron_populations.size(),
                                          std::vector<double>{});
    tables.ion_component_flux.assign(tables.ion_populations.size(), std::vector<double>{});

    double spectrum_electron_energy_moment = 0.0;
    double spectrum_electron_flux_sum = 0.0;
    double spectrum_ion_energy_moment = 0.0;
    double spectrum_ion_flux_sum = 0.0;

    const bool has_electron_spectrum =
        request.has_electron_spectrum && !request.electron_spectrum.energy_grid_ev.empty() &&
        !request.electron_spectrum.differential_number_flux.empty();
    const bool has_ion_spectrum =
        request.has_ion_spectrum && !request.ion_spectrum.energy_grid_ev.empty() &&
        !request.ion_spectrum.differential_number_flux.empty();

    const double electron_spectrum_min_ev =
        has_electron_spectrum
            ? *std::min_element(request.electron_spectrum.energy_grid_ev.begin(),
                                request.electron_spectrum.energy_grid_ev.end())
            : 0.0;
    const double electron_spectrum_max_ev =
        has_electron_spectrum
            ? *std::max_element(request.electron_spectrum.energy_grid_ev.begin(),
                                request.electron_spectrum.energy_grid_ev.end())
            : 0.0;
    const double ion_spectrum_min_ev =
        has_ion_spectrum ? *std::min_element(request.ion_spectrum.energy_grid_ev.begin(),
                                             request.ion_spectrum.energy_grid_ev.end())
                         : 0.0;
    const double ion_spectrum_max_ev =
        has_ion_spectrum ? *std::max_element(request.ion_spectrum.energy_grid_ev.begin(),
                                             request.ion_spectrum.energy_grid_ev.end())
                         : 0.0;

    tables.smoothed_electron_flux =
        has_electron_spectrum
            ? smoothLegacySpectrumLine3(request.electron_spectrum.energy_grid_ev,
                                        request.electron_spectrum.differential_number_flux)
            : std::vector<double>{};
    tables.smoothed_ion_flux =
        has_ion_spectrum
            ? smoothLegacySpectrumLine3(request.ion_spectrum.energy_grid_ev,
                                        request.ion_spectrum.differential_number_flux)
            : std::vector<double>{};
    tables.spectrum_energy_grid_ev = request.electron_spectrum.energy_grid_ev;

    for (int i = 0; i < sample_count; ++i)
    {
        const double energy_ev = std::pow(10.0, nmin + dn * static_cast<double>(i));
        const double next_energy_ev = std::pow(10.0, nmin + dn * static_cast<double>(i + 1));
        const double width_ev = next_energy_ev - energy_ev;
        tables.energy_ev.push_back(energy_ev);
        tables.width_ev.push_back(width_ev);

        double electron_component_sum = 0.0;
        for (std::size_t j = 0; j < tables.electron_populations.size(); ++j)
        {
            const auto& population = tables.electron_populations[j];
            const double value =
                legacyElectronComponentFlux(population.density_m3, population.temperature_ev, energy_ev);
            tables.electron_component_flux[j].push_back(value);
            electron_component_sum += value;
        }

        double ion_component_sum = 0.0;
        for (std::size_t j = 0; j < tables.ion_populations.size(); ++j)
        {
            const auto& population = tables.ion_populations[j];
            const double value = legacyIonComponentFlux(population.density_m3,
                                                        population.temperature_ev,
                                                        population.mass_amu, energy_ev);
            tables.ion_component_flux[j].push_back(value);
            ion_component_sum += value;
        }

        double spectrum_electron_flux = 0.0;
        if (has_electron_spectrum && energy_ev >= electron_spectrum_min_ev &&
            energy_ev <= electron_spectrum_max_ev)
        {
            spectrum_electron_flux =
                std::max(0.0, interpolateLegacySpectrum(request.electron_spectrum.energy_grid_ev,
                                                        tables.smoothed_electron_flux, energy_ev));
        }

        double spectrum_ion_flux = 0.0;
        if (has_ion_spectrum && energy_ev >= ion_spectrum_min_ev && energy_ev <= ion_spectrum_max_ev)
        {
            spectrum_ion_flux = std::max(
                0.0,
                interpolateLegacySpectrum(request.ion_spectrum.energy_grid_ev,
                                          tables.smoothed_ion_flux, energy_ev));
        }

        tables.total_electron_flux.push_back(spectrum_electron_flux > 0.0 ? spectrum_electron_flux
                                                                           : electron_component_sum);
        tables.spectrum_electron_flux.push_back(spectrum_electron_flux);
        tables.total_ion_flux.push_back(spectrum_ion_flux > 0.0 ? spectrum_ion_flux : ion_component_sum);
        tables.spectrum_ion_flux.push_back(spectrum_ion_flux);

        spectrum_electron_energy_moment += spectrum_electron_flux * energy_ev * width_ev;
        spectrum_electron_flux_sum += spectrum_electron_flux * width_ev;
        spectrum_ion_energy_moment += spectrum_ion_flux * energy_ev * width_ev;
        spectrum_ion_flux_sum += spectrum_ion_flux * width_ev;
    }

    if (spectrum_electron_flux_sum > 0.0)
    {
        tables.average_spectrum_electron_temperature_ev =
            2.0 * spectrum_electron_energy_moment / (3.0 * spectrum_electron_flux_sum);
    }
    if (spectrum_ion_flux_sum > 0.0)
    {
        tables.average_spectrum_ion_temperature_ev =
            2.0 * spectrum_ion_energy_moment / (3.0 * spectrum_ion_flux_sum);
    }

    return tables;
}

bool usesResolvedSpectrumDiscreteFlux(const ResolvedSpectrum& spectrum)
{
    switch (spectrum.model)
    {
    case SpatialSamplingModel::KAPPA:
    case SpatialSamplingModel::POWER_LAW:
    case SpatialSamplingModel::TABULATED:
        return !spectrum.energy_grid_ev.empty() &&
               spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size();
    default:
        return spectrum.populations.empty() && !spectrum.energy_grid_ev.empty() &&
               spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size();
    }
}

double integrateResolvedSpectrumFlux(
    const ResolvedSpectrum& spectrum,
    const std::function<double(double, double)>& integrand)
{
    if (spectrum.energy_grid_ev.size() < 2 ||
        spectrum.energy_grid_ev.size() != spectrum.differential_number_flux.size())
    {
        return 0.0;
    }

    double integral = 0.0;
    for (std::size_t i = 1; i < spectrum.energy_grid_ev.size(); ++i)
    {
        const double e0 = std::max(0.0, spectrum.energy_grid_ev[i - 1]);
        const double e1 = std::max(e0, spectrum.energy_grid_ev[i]);
        const double f0 = std::max(0.0, spectrum.differential_number_flux[i - 1]);
        const double f1 = std::max(0.0, spectrum.differential_number_flux[i]);
        const double width = e1 - e0;
        if (width <= 0.0)
        {
            continue;
        }
        integral += 0.5 * (integrand(e0, f0) + integrand(e1, f1)) * width;
    }
    return integral;
}

double integratePopulationMaxwellianAverageEv(
    double temperature_ev, const std::function<double(double)>& evaluator)
{
    constexpr int kSteps = 128;
    constexpr double kUpper = 32.0;
    const double thermal_ev = std::max(temperature_ev, 1.0e-6);
    const double dy = kUpper / static_cast<double>(kSteps);
    double integral = 0.0;
    for (int i = 0; i < kSteps; ++i)
    {
        const double y = (static_cast<double>(i) + 0.5) * dy;
        integral += evaluator(thermal_ev * y) * y * std::exp(-y);
    }
    return integral * dy;
}

double resolvedSpectrumAverageEnergyEv(const ResolvedSpectrum& spectrum,
                                       double fallback_energy_ev)
{
    const double denom = integrateResolvedSpectrumFlux(
        spectrum, [](double, double flux) { return flux; });
    if (!(denom > 0.0))
    {
        return std::max(1.0e-6, fallback_energy_ev);
    }

    const double numer = integrateResolvedSpectrumFlux(
        spectrum, [](double energy_ev, double flux) { return energy_ev * flux; });
    return numer / denom;
}

ResolvedSpectrumMoments resolveSpectrumMoments(const ResolvedSpectrum& spectrum,
                                               double fallback_density_m3,
                                               double fallback_characteristic_energy_ev,
                                               double fallback_mass_amu)
{
    ResolvedSpectrumMoments moments;

    double population_density_m3 = 0.0;
    double weighted_energy_ev = 0.0;
    double weighted_mass_amu = 0.0;
    for (const auto& population : spectrum.populations)
    {
        const double density_m3 = std::max(0.0, population.density_m3);
        if (!(density_m3 > 0.0))
        {
            continue;
        }

        population_density_m3 += density_m3;
        weighted_energy_ev += density_m3 * std::max(1.0e-6, population.temperature_ev);
        weighted_mass_amu += density_m3 * std::max(1.0, population.mass_amu);
    }

    if (population_density_m3 > 0.0)
    {
        moments.density_m3 = population_density_m3;
        moments.characteristic_energy_ev = weighted_energy_ev / population_density_m3;
        moments.average_mass_amu = weighted_mass_amu / population_density_m3;
        moments.used_population_moments = true;
        moments.valid = true;
        return moments;
    }

    double flux_integral = 0.0;
    double weighted_energy_flux = 0.0;
    if (spectrum.energy_grid_ev.size() >= 2 &&
        spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size())
    {
        for (std::size_t i = 1; i < spectrum.energy_grid_ev.size(); ++i)
        {
            const double e0 = std::max(0.0, spectrum.energy_grid_ev[i - 1]);
            const double e1 = std::max(e0, spectrum.energy_grid_ev[i]);
            const double f0 = std::max(0.0, spectrum.differential_number_flux[i - 1]);
            const double f1 = std::max(0.0, spectrum.differential_number_flux[i]);
            const double width = e1 - e0;
            if (width <= 0.0)
            {
                continue;
            }
            const double mean_flux = 0.5 * (f0 + f1);
            const double mean_energy_ev = 0.5 * (e0 + e1);
            flux_integral += mean_flux * width;
            weighted_energy_flux += mean_energy_ev * mean_flux * width;
        }
    }

    moments.density_m3 = std::max(0.0, fallback_density_m3);
    moments.characteristic_energy_ev =
        flux_integral > 0.0
            ? weighted_energy_flux / flux_integral
            : std::max(1.0e-6, fallback_characteristic_energy_ev);
    moments.average_mass_amu = std::max(1.0, fallback_mass_amu);
    moments.used_population_moments = false;
    moments.valid = true;
    return moments;
}

double sampleResolvedSpectrumEnergyEv(const ResolvedSpectrum& spectrum,
                                      std::mt19937& rng,
                                      double fallback_energy_ev)
{
    if (spectrum.energy_grid_ev.empty() ||
        spectrum.energy_grid_ev.size() != spectrum.differential_number_flux.size())
    {
        return std::max(1.0e-3, fallback_energy_ev);
    }

    std::vector<double> cumulative(spectrum.differential_number_flux.size(), 0.0);
    double total_flux = 0.0;
    for (std::size_t i = 0; i < spectrum.differential_number_flux.size(); ++i)
    {
        total_flux += std::max(0.0, spectrum.differential_number_flux[i]);
        cumulative[i] = total_flux;
    }

    if (!(total_flux > 0.0))
    {
        return std::max(1.0e-3, std::max(1.0e-3, spectrum.energy_grid_ev.front()));
    }

    std::uniform_real_distribution<double> dist(0.0, total_flux);
    const double pick = dist(rng);
    const auto it = std::lower_bound(cumulative.begin(), cumulative.end(), pick);
    const std::size_t index =
        static_cast<std::size_t>(std::min<std::ptrdiff_t>(
            std::distance(cumulative.begin(), it),
            static_cast<std::ptrdiff_t>(cumulative.size() - 1)));
    return std::max(1.0e-3, spectrum.energy_grid_ev[index]);
}

double resolvedSpectrumDensity(const ResolvedSpectrum& spectrum,
                              double fallback_density_m3)
{
    return resolveSpectrumMoments(spectrum, fallback_density_m3, 1.0).density_m3;
}

double resolvedSpectrumCharacteristicEnergyEv(
    const ResolvedSpectrum& spectrum, double fallback_characteristic_energy_ev)
{
    return resolveSpectrumMoments(spectrum, 0.0, fallback_characteristic_energy_ev)
        .characteristic_energy_ev;
}

double resolvedSpectrumAverageMassAmu(const ResolvedSpectrum& spectrum,
                                      double fallback_mass_amu)
{
    return resolveSpectrumMoments(spectrum, 0.0, 1.0, fallback_mass_amu).average_mass_amu;
}

double resolvedSpectrumAverageMassKg(const ResolvedSpectrum& spectrum,
                                     double fallback_mass_kg)
{
    constexpr double kAtomicMassUnit = 1.66053906660e-27;
    const double fallback_mass_amu = std::max(1.0, fallback_mass_kg / kAtomicMassUnit);
    const double mass_amu = resolvedSpectrumAverageMassAmu(spectrum, fallback_mass_amu);
    return std::max(1.0e-32, mass_amu * kAtomicMassUnit);
}

} // namespace Particle
} // namespace SCDAT
