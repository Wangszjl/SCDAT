#include "ReferenceCurrentBalanceModel.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kAtomicMassUnit = 1.66053906660e-27;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kPi = 3.14159265358979323846;

double clampPositive(double value, double minimum)
{
    return std::max(value, minimum);
}

double safeExp(double exponent)
{
    return std::exp(std::clamp(exponent, -200.0, 60.0));
}

double thermalCurrentDensity(double density_m3, double temperature_ev, double mass_kg)
{
    const double temperature_j = clampPositive(temperature_ev, 1.0e-9) * kElementaryCharge;
    return kElementaryCharge * clampPositive(density_m3, 0.0) *
           std::sqrt(temperature_j / (2.0 * kPi * std::max(1.0e-32, mass_kg)));
}

double integrateGammaLike(double temperature_ev, const std::function<double(double)>& fn)
{
    constexpr int kSteps = 128;
    constexpr double kUpper = 32.0;
    const double thermal_ev = clampPositive(temperature_ev, 1.0e-6);
    const double dy = kUpper / static_cast<double>(kSteps);
    double integral = 0.0;
    for (int i = 0; i < kSteps; ++i)
    {
        const double y = (static_cast<double>(i) + 0.5) * dy;
        integral += fn(thermal_ev * y) * y * std::exp(-y);
    }
    return integral * dy;
}

bool useDiscreteSpectrum(const Particle::ResolvedSpectrum& spectrum)
{
    switch (spectrum.model)
    {
    case Particle::SpatialSamplingModel::KAPPA:
    case Particle::SpatialSamplingModel::POWER_LAW:
    case Particle::SpatialSamplingModel::TABULATED:
        return !spectrum.energy_grid_ev.empty() &&
               spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size();
    default:
        return spectrum.populations.empty() && !spectrum.energy_grid_ev.empty() &&
               spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size();
    }
}

double integrateDiscreteFlux(const Particle::ResolvedSpectrum& spectrum,
                             const std::function<double(double, double)>& fn)
{
    if (spectrum.energy_grid_ev.size() < 2 ||
        spectrum.energy_grid_ev.size() != spectrum.differential_number_flux.size())
    {
        return 0.0;
    }

    double integral = 0.0;
    for (size_t i = 1; i < spectrum.energy_grid_ev.size(); ++i)
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
        integral += 0.5 * (fn(e0, f0) + fn(e1, f1)) * width;
    }
    return integral;
}

double averageDiscreteEnergy(const Particle::ResolvedSpectrum& spectrum)
{
    const double denom = integrateDiscreteFlux(
        spectrum, [](double, double flux) { return flux; });
    if (!(denom > 0.0))
    {
        return 1.0;
    }
    const double numer = integrateDiscreteFlux(
        spectrum, [](double energy_ev, double flux) { return energy_ev * flux; });
    return numer / denom;
}

double spectrumTemperatureScaleEv(const Particle::ResolvedSpectrum& spectrum, double fallback_ev)
{
    double weighted_temperature = 0.0;
    double total_density = 0.0;
    for (const auto& population : spectrum.populations)
    {
        weighted_temperature += population.density_m3 * population.temperature_ev;
        total_density += population.density_m3;
    }
    if (total_density > 0.0)
    {
        return weighted_temperature / total_density;
    }
    if (!spectrum.energy_grid_ev.empty())
    {
        return averageDiscreteEnergy(spectrum);
    }
    return clampPositive(fallback_ev, 1.0);
}

double spectrumDensity(const Particle::ResolvedSpectrum& spectrum, double fallback_density)
{
    double density = 0.0;
    for (const auto& population : spectrum.populations)
    {
        density += std::max(0.0, population.density_m3);
    }
    return density > 0.0 ? density : std::max(0.0, fallback_density);
}

double spectrumIonMassAmu(const Particle::ResolvedSpectrum& spectrum, double fallback_mass_amu)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (const auto& population : spectrum.populations)
    {
        numerator += std::max(0.0, population.density_m3) *
                     clampPositive(population.mass_amu, 1.0);
        denominator += std::max(0.0, population.density_m3);
    }
    return denominator > 0.0 ? numerator / denominator : clampPositive(fallback_mass_amu, 1.0);
}

template <typename PopulationFn, typename FluxFn>
double evaluateSpectrum(const Particle::ResolvedSpectrum& spectrum, PopulationFn population_fn,
                        FluxFn flux_fn)
{
    if (useDiscreteSpectrum(spectrum))
    {
        return flux_fn(spectrum);
    }

    double value = 0.0;
    for (const auto& population : spectrum.populations)
    {
        value += population_fn(population);
    }
    return value;
}

} // namespace

bool ReferenceCurrentBalanceModel::configure(const ReferenceCurrentBalanceConfig& config)
{
    config_ = config;
    configured_ = true;
    return true;
}

double ReferenceCurrentBalanceModel::safeExp(double exponent)
{
    return std::exp(std::clamp(exponent, -200.0, 60.0));
}

double ReferenceCurrentBalanceModel::erfApprox(double x)
{
    return std::erf(x);
}

double ReferenceCurrentBalanceModel::whippleYield(double peak_energy_ev, double peak_yield,
                                                  double incident_energy_ev,
                                                  double surface_potential_v)
{
    const double shifted_energy = incident_energy_ev + surface_potential_v;
    if (shifted_energy <= 1.0e-6)
    {
        return 0.0;
    }
    const double q = 2.28 * std::pow(shifted_energy / clampPositive(peak_energy_ev, 1.0), 1.35);
    if (q <= 1.0e-12)
    {
        return 0.0;
    }
    return 1.114 * (1.0 - std::exp(-q)) * peak_yield *
           std::pow(clampPositive(peak_energy_ev, 1.0) / shifted_energy, 0.35);
}

double ReferenceCurrentBalanceModel::simsYield(double peak_energy_ev, double peak_yield,
                                               double exponent_n, double incident_energy_ev,
                                               double surface_potential_v)
{
    const double shifted_energy = incident_energy_ev + surface_potential_v;
    if (shifted_energy <= 1.0e-6)
    {
        return 0.0;
    }
    const double n = clampPositive(exponent_n, 1.05);
    double xm = 1.0;
    for (int iteration = 0; iteration < 32; ++iteration)
    {
        const double residual = xm - (1.0 - 1.0 / n) * (std::exp(xm) - 1.0);
        const double derivative = 1.0 - (1.0 - 1.0 / n) * std::exp(xm);
        if (std::abs(derivative) < 1.0e-10)
        {
            break;
        }
        xm = std::clamp(xm - residual / derivative, 0.1, 5.0);
        if (std::abs(residual) < 1.0e-6)
        {
            break;
        }
    }
    const double x = xm * std::pow(shifted_energy / clampPositive(peak_energy_ev, 1.0), n);
    if (x <= 1.0e-12 || std::abs(1.0 - std::exp(-xm)) < 1.0e-12)
    {
        return 0.0;
    }
    return peak_yield / (1.0 - std::exp(-xm)) *
           std::pow(clampPositive(peak_energy_ev, 1.0) / shifted_energy, n - 1.0) * 2.0 *
           (x + std::exp(-x) - 1.0) / x;
}

double ReferenceCurrentBalanceModel::katzYield(const Material::MaterialProperty& material,
                                               double incident_energy_ev,
                                               double surface_potential_v)
{
    const double shifted_energy = incident_energy_ev + surface_potential_v;
    if (shifted_energy <= 1.0e-6)
    {
        return 0.0;
    }
    const double peak_energy_ev =
        clampPositive(material.getScalarProperty("secondary_yield_peak_energy_ev", 400.0), 1.0);
    const double peak_yield = clampPositive(material.getSecondaryElectronYield(), 0.0);
    const double r1 = material.getScalarProperty("katz_r1", 0.85);
    const double n1 = material.getScalarProperty("katz_n1", 0.35);
    const double r2 = material.getScalarProperty("katz_r2", 0.12);
    const double n2 = material.getScalarProperty("katz_n2", 1.05);
    const double energy_kev = shifted_energy * 1.0e-3;
    const double peak_energy_kev = peak_energy_ev * 1.0e-3;
    const double range = r1 * std::pow(energy_kev, n1) + r2 * std::pow(energy_kev, n2);
    const double peak_range =
        r1 * std::pow(peak_energy_kev, n1) + r2 * std::pow(peak_energy_kev, n2);
    if (range <= 1.0e-12 || peak_range <= 1.0e-12)
    {
        return 0.0;
    }
    const double attenuation = std::exp(-std::pow(range / peak_range - 1.0, 2.0));
    const double high_energy_rolloff =
        1.0 / (1.0 + std::pow(shifted_energy / clampPositive(3.0 * peak_energy_ev, 1.0), 0.7));
    return peak_yield * attenuation * (0.55 + 0.45 * high_energy_rolloff);
}

double ReferenceCurrentBalanceModel::ionSecondaryYield(double peak_energy_kev, double peak_yield,
                                                       double incident_energy_ev,
                                                       double surface_potential_v)
{
    const double incident_energy_kev =
        std::max(0.0, (incident_energy_ev - surface_potential_v) * 1.0e-3);
    if (incident_energy_kev <= 0.0)
    {
        return 0.0;
    }
    return peak_yield * std::sqrt(incident_energy_kev) *
           (1.0 + 1.0 / clampPositive(peak_energy_kev, 1.0e-6)) /
           (1.0 + incident_energy_kev / clampPositive(peak_energy_kev, 1.0e-6));
}

double ReferenceCurrentBalanceModel::backscatterYield(double incident_energy_ev,
                                                      double surface_potential_v,
                                                      double atomic_number)
{
    const double shifted_energy = incident_energy_ev + surface_potential_v;
    double ybe = 0.0;
    if (shifted_energy >= 50.0 && shifted_energy < 1.0e3)
    {
        ybe = 0.3338 * std::log(shifted_energy / 50.0) *
              (1.0 - std::pow(0.7358, 0.037 * atomic_number) +
               0.1 * std::exp(-shifted_energy / 50.0));
    }
    else if (shifted_energy >= 1.0e3 && shifted_energy < 1.0e4)
    {
        ybe = 1.0 - std::pow(0.7358, 0.037 * atomic_number) +
              0.1 * std::exp(-shifted_energy / 5000.0);
    }
    else if (shifted_energy >= 1.0e4 && shifted_energy < 1.0e5)
    {
        ybe = 1.0 - std::pow(0.7358, 0.037 * atomic_number);
    }
    return std::max(0.0, ybe);
}

double ReferenceCurrentBalanceModel::bodyPhotoCurrentDensity(double base_photo_current_density_a_per_m2)
{
    return 0.25 * std::max(0.0, base_photo_current_density_a_per_m2);
}

double ReferenceCurrentBalanceModel::patchPhotoCurrentDensity(
    double base_photo_current_density_a_per_m2, double incidence_angle_deg)
{
    return std::max(0.0, base_photo_current_density_a_per_m2) *
           std::max(0.0, std::cos(incidence_angle_deg * kPi / 180.0));
}

double ReferenceCurrentBalanceModel::conductionCurrentDensity(double body_potential_v,
                                                              double patch_potential_v,
                                                              double conductivity_s_per_m,
                                                              double thickness_m)
{
    return clampPositive(conductivity_s_per_m, 0.0) * (body_potential_v - patch_potential_v) /
           clampPositive(thickness_m, 1.0e-9);
}

double ReferenceCurrentBalanceModel::emissionEscapeProbability(double surface_potential_v,
                                                               double characteristic_energy_ev) const
{
    if (surface_potential_v <= 0.0)
    {
        return 1.0;
    }
    return std::clamp(
        safeExp(-surface_potential_v / clampPositive(characteristic_energy_ev, 1.0e-3)), 0.0, 1.0);
}

double ReferenceCurrentBalanceModel::rolePotential(ReferenceSurfaceRole role, double body_potential_v,
                                                   double patch_potential_v) const
{
    return role == ReferenceSurfaceRole::Body ? body_potential_v : patch_potential_v;
}

double ReferenceCurrentBalanceModel::computeSeeYield(const Material::MaterialProperty& material,
                                                     double incident_energy_ev,
                                                     double surface_potential_v) const
{
    const double peak_energy =
        clampPositive(material.getScalarProperty("secondary_yield_peak_energy_ev", 300.0), 1.0);
    const double peak_yield = clampPositive(material.getSecondaryElectronYield(), 0.0);
    switch (config_.see_model)
    {
    case SecondaryElectronEmissionModel::Whipple:
        return whippleYield(peak_energy, peak_yield, incident_energy_ev, surface_potential_v);
    case SecondaryElectronEmissionModel::Sims:
        return simsYield(peak_energy, peak_yield,
                         material.getScalarProperty("sims_exponent_n", 1.6), incident_energy_ev,
                         surface_potential_v);
    case SecondaryElectronEmissionModel::Katz:
        return katzYield(material, incident_energy_ev, surface_potential_v);
    }
    return 0.0;
}

double ReferenceCurrentBalanceModel::ramBodyCurrentDensity(double surface_potential_v,
                                                           double ion_density_m3,
                                                           double ion_temperature_ev,
                                                           double ion_mass_amu,
                                                           double flow_speed_m_per_s)
{
    const double ni_cm3 = clampPositive(ion_density_m3 * 1.0e-6, 0.0);
    const double ti_ev = clampPositive(ion_temperature_ev, 1.0e-6);
    const double flow_speed_km_s = flow_speed_m_per_s * 1.0e-3;
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / clampPositive(ion_mass_amu, 1.0));
    const double qd = flow_speed_km_s / clampPositive(wt, 1.0e-12);
    double current_na_per_m2 = 0.0;
    if (surface_potential_v >= 0.0)
    {
        current_na_per_m2 =
            0.5642 / clampPositive(qd, 1.0e-12) *
            ((qv + qd) * std::exp(-(qv - qd) * (qv - qd)) -
             (qv - qd) * std::exp(-(qv + qd) * (qv + qd)));
        current_na_per_m2 +=
            (0.5 / clampPositive(qd, 1.0e-12) + qd - qv * qv / clampPositive(qd, 1.0e-12)) *
            (erfApprox(qv + qd) - erfApprox(qv - qd));
        current_na_per_m2 *= 0.02003 * ni_cm3 * wt;
    }
    else
    {
        current_na_per_m2 =
            0.04005 * ni_cm3 * wt *
            (0.5642 * std::exp(-qd * qd) +
             (qd + 0.5 / clampPositive(qd, 1.0e-12)) * erfApprox(qd));
        current_na_per_m2 *= 4.0;
    }
    return current_na_per_m2 * 1.0e-9;
}

double ReferenceCurrentBalanceModel::ramPatchCurrentDensity(double surface_potential_v,
                                                            double ion_density_m3,
                                                            double ion_temperature_ev,
                                                            double ion_mass_amu,
                                                            double flow_speed_m_per_s,
                                                            double patch_flow_angle_deg)
{
    const double ni_cm3 = clampPositive(ion_density_m3 * 1.0e-6, 0.0);
    const double ti_ev = clampPositive(ion_temperature_ev, 1.0e-6);
    const double flow_speed_km_s = flow_speed_m_per_s * 1.0e-3;
    const double alpha = patch_flow_angle_deg * kPi / 180.0;
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / clampPositive(ion_mass_amu, 1.0));
    const double qd = flow_speed_km_s * std::cos(alpha) / clampPositive(wt, 1.0e-12);
    double current_na_per_m2 = 0.08011 * ni_cm3 * wt;
    if (surface_potential_v >= 0.0)
    {
        current_na_per_m2 *=
            0.5642 * std::exp(-(qd - qv) * (qd - qv)) + qd + qd * erfApprox(qd - qv);
    }
    else
    {
        current_na_per_m2 *= 0.5642 * std::exp(-qd * qd) + qd + qd * erfApprox(qd);
    }
    return std::max(0.0, current_na_per_m2) * 1.0e-9;
}

double ReferenceCurrentBalanceModel::computeElectronCollection(double surface_potential_v) const
{
    const auto population_current = [&](const Particle::SpectrumPopulation& population) {
        const double base_flux = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = surface_potential_v <= 0.0
                                  ? safeExp(surface_potential_v /
                                            clampPositive(population.temperature_ev, 1.0e-3))
                                  : (1.0 + surface_potential_v /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return -base_flux * clampPositive(config_.electron_collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * factor;
    };
    const auto flux_current = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double mean_energy_ev = averageDiscreteEnergy(spectrum);
        const double current = integrateDiscreteFlux(spectrum, [&](double energy_ev, double flux) {
            if (surface_potential_v <= 0.0 && energy_ev + surface_potential_v <= 0.0)
            {
                return 0.0;
            }
            const double enhancement =
                surface_potential_v > 0.0
                    ? (1.0 + surface_potential_v / clampPositive(mean_energy_ev, 1.0e-3))
                    : 1.0;
            return flux * enhancement;
        });
        return -kElementaryCharge * clampPositive(config_.electron_collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * current;
    };
    if (config_.has_electron_spectrum)
    {
        return evaluateSpectrum(config_.electron_spectrum, population_current, flux_current);
    }
    const double base_flux = thermalCurrentDensity(
        config_.plasma.electron_density_m3, config_.plasma.electron_temperature_ev, kElectronMass);
    const double factor = surface_potential_v <= 0.0
                              ? safeExp(surface_potential_v /
                                        clampPositive(config_.plasma.electron_temperature_ev, 1.0e-3))
                              : (1.0 + surface_potential_v /
                                           clampPositive(config_.plasma.electron_temperature_ev, 1.0e-3));
    return -base_flux * clampPositive(config_.electron_collection_coefficient, 0.0) *
           clampPositive(config_.electron_calibration_factor, 0.1) * factor;
}

double ReferenceCurrentBalanceModel::computeIonCollection(double surface_potential_v) const
{
    const auto population_current = [&](const Particle::SpectrumPopulation& population) {
        const double base_flux = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = surface_potential_v >= 0.0
                                  ? safeExp(-surface_potential_v /
                                            clampPositive(population.temperature_ev, 1.0e-3))
                                  : (1.0 - surface_potential_v /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return base_flux * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * factor;
    };
    const auto flux_current = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double mean_energy_ev = averageDiscreteEnergy(spectrum);
        const double current = integrateDiscreteFlux(spectrum, [&](double energy_ev, double flux) {
            if (surface_potential_v > 0.0 && energy_ev - surface_potential_v <= 0.0)
            {
                return 0.0;
            }
            const double enhancement =
                surface_potential_v < 0.0
                    ? (1.0 + (-surface_potential_v) / clampPositive(mean_energy_ev, 1.0e-3))
                    : 1.0;
            return flux * enhancement;
        });
        return kElementaryCharge * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * current;
    };
    if (config_.has_ion_spectrum)
    {
        return evaluateSpectrum(config_.ion_spectrum, population_current, flux_current);
    }
    const double ion_mass = clampPositive(config_.plasma.ion_mass_amu, 1.0) * kAtomicMassUnit;
    const double base_flux = thermalCurrentDensity(
        config_.plasma.ion_density_m3, config_.plasma.ion_temperature_ev, ion_mass);
    const double factor = surface_potential_v >= 0.0
                              ? safeExp(-surface_potential_v /
                                        clampPositive(config_.plasma.ion_temperature_ev, 1.0e-3))
                              : (1.0 - surface_potential_v /
                                           clampPositive(config_.plasma.ion_temperature_ev, 1.0e-3));
    return base_flux * clampPositive(config_.ion_collection_coefficient, 0.0) *
           clampPositive(config_.ion_calibration_factor, 0.1) * factor;
}

ReferenceCurrentComponents
ReferenceCurrentBalanceModel::evaluate(ReferenceSurfaceRole role, double body_potential_v,
                                       double patch_potential_v) const
{
    ReferenceCurrentComponents components;
    if (!configured_)
    {
        return components;
    }

    const auto& material =
        role == ReferenceSurfaceRole::Body ? config_.body_material : config_.patch_material;
    const double potential = rolePotential(role, body_potential_v, patch_potential_v);
    components.electron_collection_a_per_m2 = computeElectronCollection(potential);
    components.ion_collection_a_per_m2 = computeIonCollection(potential);

    const double secondary_escape = emissionEscapeProbability(
        potential, material.getScalarProperty("secondary_emission_escape_energy_ev", 2.0));
    const double photo_escape = config_.use_photoelectron_suppression
                                    ? emissionEscapeProbability(potential,
                                                                config_.photoelectron_temperature_ev)
                                    : 1.0;

    const auto electron_secondary_from_pop = [&](const Particle::SpectrumPopulation& population) {
        const double average_yield = integrateGammaLike(population.temperature_ev, [&](double energy_ev) {
            return computeSeeYield(material, energy_ev + std::max(0.0, potential), 0.0);
        });
        const double current = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = potential <= 0.0
                                  ? safeExp(potential /
                                            clampPositive(population.temperature_ev, 1.0e-3))
                                  : (1.0 + potential /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return current * clampPositive(config_.electron_collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * factor * average_yield *
               secondary_escape;
    };
    const auto electron_secondary_from_flux = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double flux = integrateDiscreteFlux(spectrum, [&](double energy_ev, double value) {
            if (potential <= 0.0 && energy_ev + potential <= 0.0)
            {
                return 0.0;
            }
            return value * computeSeeYield(material, energy_ev + std::max(0.0, potential), 0.0);
        });
        return flux * kElementaryCharge *
               clampPositive(config_.electron_collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * secondary_escape;
    };
    const auto ion_secondary_from_pop = [&](const Particle::SpectrumPopulation& population) {
        const double average_yield = integrateGammaLike(population.temperature_ev, [&](double energy_ev) {
            return ionSecondaryYield(
                material.getScalarProperty("ion_secondary_peak_energy_kev", 0.35),
                clampPositive(
                    material.getScalarProperty("ion_secondary_yield", config_.ion_secondary_yield), 0.0),
                energy_ev + std::max(0.0, -potential), 0.0);
        });
        const double current = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = potential >= 0.0
                                  ? safeExp(-potential /
                                            clampPositive(population.temperature_ev, 1.0e-3))
                                  : (1.0 - potential /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return current * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * factor * average_yield *
               secondary_escape;
    };
    const auto ion_secondary_from_flux = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double flux = integrateDiscreteFlux(spectrum, [&](double energy_ev, double value) {
            if (potential > 0.0 && energy_ev - potential <= 0.0)
            {
                return 0.0;
            }
            return value * ionSecondaryYield(
                               material.getScalarProperty("ion_secondary_peak_energy_kev", 0.35),
                               clampPositive(material.getScalarProperty("ion_secondary_yield",
                                                                       config_.ion_secondary_yield),
                                             0.0),
                               energy_ev + std::max(0.0, -potential), 0.0);
        });
        return flux * kElementaryCharge * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * secondary_escape;
    };
    const auto backscatter_from_pop = [&](const Particle::SpectrumPopulation& population) {
        const double average_yield = integrateGammaLike(population.temperature_ev, [&](double energy_ev) {
            return backscatterYield(
                energy_ev + std::max(0.0, potential), 0.0,
                material.getScalarProperty("atomic_number", config_.backscatter_atomic_number));
        });
        const double current = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = potential <= 0.0
                                  ? safeExp(potential /
                                            clampPositive(population.temperature_ev, 1.0e-3))
                                  : (1.0 + potential /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return current * clampPositive(config_.electron_collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * factor * average_yield *
               secondary_escape;
    };
    const auto backscatter_from_flux = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double flux = integrateDiscreteFlux(spectrum, [&](double energy_ev, double value) {
            if (potential <= 0.0 && energy_ev + potential <= 0.0)
            {
                return 0.0;
            }
            return value * backscatterYield(
                               energy_ev + std::max(0.0, potential), 0.0,
                               material.getScalarProperty("atomic_number",
                                                          config_.backscatter_atomic_number));
        });
        return flux * kElementaryCharge *
               clampPositive(config_.electron_collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * secondary_escape;
    };

    if (config_.has_electron_spectrum)
    {
        components.secondary_electron_a_per_m2 =
            evaluateSpectrum(config_.electron_spectrum, electron_secondary_from_pop,
                             electron_secondary_from_flux);
        components.backscatter_electron_a_per_m2 =
            evaluateSpectrum(config_.electron_spectrum, backscatter_from_pop, backscatter_from_flux);
    }
    else
    {
        const double see_yield = integrateGammaLike(config_.plasma.electron_temperature_ev, [&](double energy_ev) {
            return computeSeeYield(material, energy_ev + std::max(0.0, potential), 0.0);
        });
        const double backscatter_yield = integrateGammaLike(
            config_.plasma.electron_temperature_ev, [&](double energy_ev) {
                return backscatterYield(
                    energy_ev + std::max(0.0, potential), 0.0,
                    material.getScalarProperty("atomic_number", config_.backscatter_atomic_number));
            });
        components.secondary_electron_a_per_m2 =
            std::abs(components.electron_collection_a_per_m2) * see_yield * secondary_escape;
        components.backscatter_electron_a_per_m2 =
            std::abs(components.electron_collection_a_per_m2) * backscatter_yield * secondary_escape;
    }

    if (config_.has_ion_spectrum)
    {
        components.ion_secondary_electron_a_per_m2 =
            evaluateSpectrum(config_.ion_spectrum, ion_secondary_from_pop, ion_secondary_from_flux);
    }
    else
    {
        const double ion_secondary_yield = integrateGammaLike(
            config_.plasma.ion_temperature_ev, [&](double energy_ev) {
                return ionSecondaryYield(
                    material.getScalarProperty("ion_secondary_peak_energy_kev", 0.35),
                    clampPositive(material.getScalarProperty("ion_secondary_yield",
                                                             config_.ion_secondary_yield),
                                  0.0),
                    energy_ev + std::max(0.0, -potential), 0.0);
            });
        components.ion_secondary_electron_a_per_m2 =
            std::abs(components.ion_collection_a_per_m2) * ion_secondary_yield * secondary_escape;
    }

    const double base_photo_current =
        role == ReferenceSurfaceRole::Body
            ? bodyPhotoCurrentDensity(config_.body_photo_current_density_a_per_m2)
            : patchPhotoCurrentDensity(config_.patch_photo_current_density_a_per_m2,
                                       config_.patch_incidence_angle_deg);
    components.photoelectron_a_per_m2 = base_photo_current * photo_escape;

    if (role == ReferenceSurfaceRole::Patch)
    {
        components.conduction_a_per_m2 = conductionCurrentDensity(
            body_potential_v, patch_potential_v, config_.patch_conductivity_s_per_m,
            config_.patch_thickness_m);
    }

    if (config_.enable_ram_current)
    {
        const double ion_density = config_.has_ion_spectrum
                                       ? spectrumDensity(config_.ion_spectrum,
                                                         config_.plasma.ion_density_m3)
                                       : config_.plasma.ion_density_m3;
        const double ion_temperature =
            config_.has_ion_spectrum
                ? spectrumTemperatureScaleEv(config_.ion_spectrum, config_.plasma.ion_temperature_ev)
                : config_.plasma.ion_temperature_ev;
        const double ion_mass_amu =
            config_.has_ion_spectrum
                ? spectrumIonMassAmu(config_.ion_spectrum, config_.plasma.ion_mass_amu)
                : config_.plasma.ion_mass_amu;
        if (role == ReferenceSurfaceRole::Body)
        {
            components.ram_ion_a_per_m2 = ramBodyCurrentDensity(
                potential, ion_density, ion_temperature, ion_mass_amu,
                std::max(0.0, config_.bulk_flow_velocity_m_per_s *
                                  std::clamp(config_.flow_alignment_cosine, -1.0, 1.0)));
        }
        else
        {
            components.ram_ion_a_per_m2 = ramPatchCurrentDensity(
                potential, ion_density, ion_temperature, ion_mass_amu,
                std::max(0.0, config_.bulk_flow_velocity_m_per_s), config_.patch_flow_angle_deg);
        }
    }

    components.net_a_per_m2 = components.electron_collection_a_per_m2 +
                              components.ion_collection_a_per_m2 +
                              components.secondary_electron_a_per_m2 +
                              components.ion_secondary_electron_a_per_m2 +
                              components.backscatter_electron_a_per_m2 +
                              components.photoelectron_a_per_m2 +
                              components.conduction_a_per_m2 + components.ram_ion_a_per_m2;
    return components;
}

double ReferenceCurrentBalanceModel::computeNetCurrentDensity(ReferenceSurfaceRole role,
                                                              double body_potential_v,
                                                              double patch_potential_v) const
{
    return evaluate(role, body_potential_v, patch_potential_v).net_a_per_m2;
}

double ReferenceCurrentBalanceModel::computeCurrentDerivative(ReferenceSurfaceRole role,
                                                              double body_potential_v,
                                                              double patch_potential_v) const
{
    const double role_potential = rolePotential(role, body_potential_v, patch_potential_v);
    const double temperature_scale = config_.has_electron_spectrum
                                         ? spectrumTemperatureScaleEv(config_.electron_spectrum,
                                                                      config_.plasma.electron_temperature_ev)
                                         : clampPositive(config_.plasma.electron_temperature_ev, 1.0);
    const double step = std::max(1.0e-3, 1.0e-2 * std::max(1.0, std::abs(role_potential)) +
                                           1.0e-2 * temperature_scale);
    if (role == ReferenceSurfaceRole::Body)
    {
        const double jp = computeNetCurrentDensity(role, body_potential_v + step, patch_potential_v);
        const double jm = computeNetCurrentDensity(role, body_potential_v - step, patch_potential_v);
        return (jp - jm) / (2.0 * step);
    }
    const double jp = computeNetCurrentDensity(role, body_potential_v, patch_potential_v + step);
    const double jm = computeNetCurrentDensity(role, body_potential_v, patch_potential_v - step);
    return (jp - jm) / (2.0 * step);
}

double ReferenceCurrentBalanceModel::solveEquilibriumPotential(ReferenceSurfaceRole role,
                                                               double body_potential_v,
                                                               double patch_potential_v,
                                                               double minimum_potential_v,
                                                               double maximum_potential_v) const
{
    const double temperature_scale = config_.has_electron_spectrum
                                         ? spectrumTemperatureScaleEv(config_.electron_spectrum,
                                                                      config_.plasma.electron_temperature_ev)
                                         : clampPositive(config_.plasma.electron_temperature_ev, 1.0);
    double lower = std::isfinite(minimum_potential_v)
                       ? minimum_potential_v
                       : -std::max(5.0e4, 40.0 * temperature_scale);
    double upper = std::isfinite(maximum_potential_v)
                       ? maximum_potential_v
                       : std::max(5.0e4, 40.0 * temperature_scale);
    auto evaluate_at = [&](double candidate) {
        return role == ReferenceSurfaceRole::Body
                   ? computeNetCurrentDensity(role, candidate, patch_potential_v)
                   : computeNetCurrentDensity(role, body_potential_v, candidate);
    };
    double f_lower = evaluate_at(lower);
    double f_upper = evaluate_at(upper);
    for (int expand = 0; expand < 12 && f_lower * f_upper > 0.0; ++expand)
    {
        lower *= 1.5;
        upper *= 1.5;
        f_lower = evaluate_at(lower);
        f_upper = evaluate_at(upper);
    }
    if (f_lower * f_upper > 0.0)
    {
        return std::abs(f_lower) < std::abs(f_upper) ? lower : upper;
    }
    for (int iteration = 0; iteration < 128; ++iteration)
    {
        const double middle = 0.5 * (lower + upper);
        const double f_middle = evaluate_at(middle);
        if (std::abs(f_middle) < 1.0e-12 || std::abs(upper - lower) < 1.0e-6)
        {
            return middle;
        }
        if (f_lower * f_middle <= 0.0)
        {
            upper = middle;
            f_upper = f_middle;
        }
        else
        {
            lower = middle;
            f_lower = f_middle;
        }
    }
    return 0.5 * (lower + upper);
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
