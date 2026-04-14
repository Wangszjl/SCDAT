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
