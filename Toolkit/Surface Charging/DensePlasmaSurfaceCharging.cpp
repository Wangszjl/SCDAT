#include "DensePlasmaSurfaceCharging.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kAtomicMassUnit = 1.66053906660e-27;
constexpr double kBoltzmannEvPerK = 8.617333262145e-5;
constexpr double kRichardsonConstant = 1.20173e6;
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kPi = 3.14159265358979323846;

double clampSigned(double value, double limit)
{
    return std::clamp(value, -limit, limit);
}

double safeExp(double exponent)
{
    return std::exp(std::clamp(exponent, -200.0, 50.0));
}

double emittedElectronEscapeProbability(double surface_potential_v, double characteristic_energy_ev)
{
    if (surface_potential_v <= 0.0)
    {
        return 1.0;
    }

    const double escape_energy = std::max(1.0e-3, characteristic_energy_ev);
    return std::clamp(safeExp(-surface_potential_v / escape_energy), 0.0, 1.0);
}

double halfNormalFluxAverage(double thermal_sigma, std::mt19937_64& generator)
{
    std::normal_distribution<double> distribution(0.0, thermal_sigma);
    return std::abs(distribution(generator));
}
}

bool DensePlasmaSurfaceCharging::initialize(const SurfaceChargingConfig& config)
{
    config_ = config;
    status_ = SurfaceChargingStatus{};
    status_.state.capacitance_per_area_f_per_m2 = computeCapacitancePerArea();
    initialized_ = true;
    electron_pic_calibration_factor_ = 1.0;
    ion_pic_calibration_factor_ = 1.0;
    history_time_.clear();
    history_potential_.clear();
    history_charge_.clear();
    history_current_.clear();
    history_electron_current_.clear();
    history_ion_current_.clear();
    history_secondary_current_.clear();
    history_photo_current_.clear();
    history_thermionic_current_.clear();
    history_field_emission_current_.clear();
    history_leakage_current_.clear();
    history_net_current_.clear();
    history_capacitance_.clear();
    history_effective_conductivity_.clear();
    history_effective_sheath_length_.clear();
    history_adaptive_time_step_.clear();
    history_electron_calibration_factor_.clear();
    history_ion_calibration_factor_.clear();
    history_equilibrium_potential_.clear();
    history_equilibrium_error_.clear();
    applyGeoPicCalibration();
    return true;
}

double DensePlasmaSurfaceCharging::computeCapacitancePerArea() const
{
    const double thickness = std::max(1.0e-8, config_.dielectric_thickness_m);
    const double eps_r = std::max(1.0, config_.material.getPermittivity());
    if (config_.derive_capacitance_from_material)
    {
        return kEpsilon0 * eps_r / thickness;
    }
    return std::max(1.0e-12, config_.capacitance_per_area_f_per_m2);
}

double DensePlasmaSurfaceCharging::computeEffectiveSheathLength() const
{
    if (!config_.derive_sheath_length_from_plasma)
    {
        return std::max(1.0e-6, config_.sheath_length_m);
    }

    const double electron_temperature_j =
        std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_density = std::max(1.0e3, config_.plasma.electron_density_m3);
    const double debye_length = std::sqrt(kEpsilon0 * electron_temperature_j /
                                          (electron_density * kElementaryCharge * kElementaryCharge));
    return std::clamp(debye_length, std::max(1.0e-6, config_.minimum_sheath_length_m),
                      std::max(config_.minimum_sheath_length_m, config_.maximum_sheath_length_m));
}

SurfaceCurrents DensePlasmaSurfaceCharging::computeSurfaceCurrents(double surface_potential_v) const
{
    SurfaceCurrents currents;
    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double bohm_velocity = std::sqrt(te_j / ion_mass);
    const double flow_alignment = std::clamp(config_.flow_alignment_cosine, -1.0, 1.0);
    const double normal_flow_speed =
        std::max(0.0, config_.bulk_flow_velocity_m_per_s * flow_alignment);
    const double directed_ion_velocity =
        std::max(0.0, config_.ion_directed_velocity_m_per_s + normal_flow_speed);
    const double effective_ion_velocity =
        std::sqrt(bohm_velocity * bohm_velocity + directed_ion_velocity * directed_ion_velocity);
    const double phi_scale = std::max(0.25, config_.plasma.electron_temperature_ev);
    double electron_flux =
        kElementaryCharge * config_.plasma.electron_density_m3 * electron_thermal_velocity *
        std::max(0.0, config_.electron_collection_coefficient);
    double ion_flux = kElementaryCharge * config_.plasma.ion_density_m3 * effective_ion_velocity *
                      std::max(0.0, config_.ion_collection_coefficient);

    if (config_.regime == SurfaceChargingRegime::LeoFlowingPlasma)
    {
        const double ion_advective_flux =
            kElementaryCharge * config_.plasma.ion_density_m3 * normal_flow_speed *
            std::max(0.0, config_.ion_collection_coefficient);
        const double electron_advective_flux =
            kElementaryCharge * config_.plasma.electron_density_m3 * normal_flow_speed *
            std::clamp(config_.electron_flow_coupling, 0.0, 1.0);
        ion_flux += ion_advective_flux;
        electron_flux += electron_advective_flux;
    }

    if (config_.regime == SurfaceChargingRegime::GeoKineticPicLike)
    {
        electron_flux *= std::clamp(electron_pic_calibration_factor_, 0.25, 4.0);
        ion_flux *= std::clamp(ion_pic_calibration_factor_, 0.25, 4.0);
    }

    double raw_electron_current = 0.0;
    if (surface_potential_v <= 0.0)
    {
        raw_electron_current = -electron_flux * safeExp(surface_potential_v / phi_scale);
        currents.ion_current_a_per_m2 =
            ion_flux * std::sqrt(1.0 + std::min(50.0, -surface_potential_v / phi_scale));
    }
    else
    {
        raw_electron_current = -electron_flux * (1.0 + std::min(5.0, surface_potential_v / phi_scale));
        currents.ion_current_a_per_m2 = ion_flux * safeExp(-surface_potential_v / phi_scale);
    }

    Material::SurfaceImpactState impact;
    impact.incident_energy_ev =
        std::max(1.0e-2, config_.plasma.electron_temperature_ev + surface_potential_v);
    impact.particle_charge_coulomb = -kElementaryCharge;
    impact.surface_temperature_k = config_.emission.surface_temperature_k;
    const auto interaction = surface_interaction_.evaluate(config_.material, impact);
    const double secondary_escape_probability =
        emittedElectronEscapeProbability(surface_potential_v,
                                         config_.material.getScalarProperty(
                                             "secondary_emission_escape_energy_ev", 2.0));
    const double photo_escape_probability =
        emittedElectronEscapeProbability(surface_potential_v,
                                         config_.material.getScalarProperty(
                                             "photoelectron_escape_energy_ev", 1.5));
    const double thermionic_escape_probability =
        emittedElectronEscapeProbability(surface_potential_v,
                                         std::max(1.0e-3,
                                                  config_.material.getScalarProperty(
                                                      "thermionic_escape_energy_ev",
                                                      2.0 * kBoltzmannEvPerK *
                                                          config_.emission.surface_temperature_k)));

    const double absorbed_fraction =
        impact.particle_charge_coulomb != 0.0
            ? std::clamp(std::abs(interaction.absorbed_charge_coulomb / impact.particle_charge_coulomb), 0.0,
                         1.0)
            : 1.0;
    currents.electron_current_a_per_m2 = raw_electron_current * absorbed_fraction;
    currents.secondary_emission_a_per_m2 =
        std::abs(currents.electron_current_a_per_m2) *
        std::clamp(interaction.secondary_emitted_electrons, 0.0, 2.5) *
        secondary_escape_probability;

    currents.photo_emission_a_per_m2 =
        kElementaryCharge * config_.emission.photon_flux_m2_s *
        std::clamp(config_.material.getScalarProperty("photoelectron_yield", 0.02) *
                       config_.emission.enhancement_factor,
                   0.0, 0.5) *
        photo_escape_probability;

    const double work_function = std::max(0.1, config_.material.getWorkFunctionEv());
    const double thermal_energy_ev =
        std::max(1.0e-6, kBoltzmannEvPerK * config_.emission.surface_temperature_k);
    currents.thermionic_emission_a_per_m2 =
        kRichardsonConstant * config_.emission.surface_temperature_k * config_.emission.surface_temperature_k *
        safeExp(-work_function / thermal_energy_ev) * thermionic_escape_probability;

    const double sheath_length = computeEffectiveSheathLength();
    const double local_field_v_per_m =
        std::abs(surface_potential_v) * std::max(1.0, config_.emission.enhancement_factor) / sheath_length;
    if (local_field_v_per_m > 1.0e7)
    {
        const double fn_prefactor = 1.54e-6 * local_field_v_per_m * local_field_v_per_m / work_function;
        const double fn_exponent =
            -6.83e9 * std::pow(work_function, 1.5) / std::max(1.0e6, local_field_v_per_m);
        currents.field_emission_a_per_m2 = fn_prefactor * safeExp(fn_exponent);
    }

    currents.electron_current_a_per_m2 =
        clampSigned(currents.electron_current_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.ion_current_a_per_m2 =
        clampSigned(currents.ion_current_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.secondary_emission_a_per_m2 =
        clampSigned(currents.secondary_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.photo_emission_a_per_m2 =
        clampSigned(currents.photo_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.thermionic_emission_a_per_m2 =
        clampSigned(currents.thermionic_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.field_emission_a_per_m2 =
        clampSigned(currents.field_emission_a_per_m2, config_.max_abs_current_density_a_per_m2);
    currents.total_current_a_per_m2 = currents.electron_current_a_per_m2 +
                                      currents.ion_current_a_per_m2 +
                                      currents.secondary_emission_a_per_m2 +
                                      currents.photo_emission_a_per_m2 +
                                      currents.thermionic_emission_a_per_m2 +
                                      currents.field_emission_a_per_m2;
    currents.total_current_a_per_m2 =
        clampSigned(currents.total_current_a_per_m2, config_.max_abs_current_density_a_per_m2);
    return currents;
}

double DensePlasmaSurfaceCharging::computeRadiationInducedConductivity() const
{
    if (config_.radiation_conductivity_coefficient <= 0.0)
    {
        return 0.0;
    }

    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double electron_flux =
        kElementaryCharge * config_.plasma.electron_density_m3 * electron_thermal_velocity;
    const double incident_power_flux =
        std::max(0.0, std::abs(electron_flux) * std::max(1.0e-3, config_.plasma.electron_temperature_ev));
    return config_.radiation_conductivity_coefficient *
           std::pow(std::max(1.0e-12, incident_power_flux), config_.radiation_conductivity_exponent);
}

void DensePlasmaSurfaceCharging::applyGeoPicCalibration()
{
    if (config_.regime != SurfaceChargingRegime::GeoKineticPicLike || !config_.enable_pic_calibration)
    {
        return;
    }

    const std::size_t samples = std::max<std::size_t>(256, config_.pic_calibration_samples);
    const double base_scale = std::max(25.0, 0.2 * config_.plasma.electron_temperature_ev);
    const double reference_potential = computeFloatingPotential();
    const std::array<double, 3> calibration_window{
        reference_potential - base_scale,
        reference_potential,
        reference_potential + base_scale};

    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double electron_thermal_velocity = std::sqrt(te_j / (2.0 * kPi * kElectronMass));
    const double bohm_velocity = std::sqrt(te_j / ion_mass);

    double electron_ratio_sum = 0.0;
    double ion_ratio_sum = 0.0;
    std::size_t electron_ratio_count = 0;
    std::size_t ion_ratio_count = 0;

    for (const double potential_v : calibration_window)
    {
        const double pic_electron_flux = estimateElectronFluxPicLike(potential_v, samples);
        const double pic_ion_flux = estimateIonFluxPicLike(potential_v, samples);
        const double analytic_electron_flux =
            kElementaryCharge * config_.plasma.electron_density_m3 * electron_thermal_velocity *
            std::max(0.0, config_.electron_collection_coefficient);
        const double analytic_ion_flux =
            kElementaryCharge * config_.plasma.ion_density_m3 * bohm_velocity *
            std::max(0.0, config_.ion_collection_coefficient);

        if (analytic_electron_flux > 1.0e-18)
        {
            const double ratio = std::clamp(pic_electron_flux / analytic_electron_flux, 0.4, 2.5);
            electron_ratio_sum += ratio;
            electron_ratio_count += 1;
        }
        if (analytic_ion_flux > 1.0e-18)
        {
            const double ratio = std::clamp(pic_ion_flux / analytic_ion_flux, 0.4, 2.5);
            ion_ratio_sum += ratio;
            ion_ratio_count += 1;
        }
    }

    if (electron_ratio_count > 0)
    {
        electron_pic_calibration_factor_ = std::clamp(electron_ratio_sum / electron_ratio_count, 0.5, 2.0);
    }
    if (ion_ratio_count > 0)
    {
        ion_pic_calibration_factor_ = std::clamp(ion_ratio_sum / ion_ratio_count, 0.5, 2.0);
    }
}

double DensePlasmaSurfaceCharging::estimateElectronFluxPicLike(double surface_potential_v,
                                                               std::size_t samples) const
{
    const std::size_t sample_count = std::max<std::size_t>(64, samples);
    const double te_j = std::max(1.0e-3, config_.plasma.electron_temperature_ev) * kElementaryCharge;
    const double sigma = std::sqrt(te_j / kElectronMass);
    std::mt19937_64 generator(0x5cdab51ULL + static_cast<std::uint64_t>(sample_count));

    double collected_surface_velocity_sum = 0.0;
    for (std::size_t index = 0; index < sample_count; ++index)
    {
        const double normal_velocity = halfNormalFluxAverage(sigma, generator);
        double normal_energy_j = 0.5 * kElectronMass * normal_velocity * normal_velocity;

        if (surface_potential_v < 0.0)
        {
            const double barrier_j = -kElementaryCharge * surface_potential_v;
            if (normal_energy_j <= barrier_j)
            {
                continue;
            }
            normal_energy_j -= barrier_j;
        }
        else
        {
            normal_energy_j += kElementaryCharge * surface_potential_v;
        }

        const double surface_velocity =
            std::sqrt(std::max(0.0, 2.0 * normal_energy_j / kElectronMass));
        collected_surface_velocity_sum += surface_velocity;
    }

    const double mean_surface_velocity =
        collected_surface_velocity_sum / static_cast<double>(sample_count);
    return 0.5 * kElementaryCharge * config_.plasma.electron_density_m3 * mean_surface_velocity *
           std::max(0.0, config_.electron_collection_coefficient);
}

double DensePlasmaSurfaceCharging::estimateIonFluxPicLike(double surface_potential_v,
                                                          std::size_t samples) const
{
    const std::size_t sample_count = std::max<std::size_t>(64, samples);
    const double ion_mass = std::max(1.0, config_.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double ti_j = std::max(1.0e-3, config_.plasma.ion_temperature_ev) * kElementaryCharge;
    const double sigma = std::sqrt(ti_j / ion_mass);
    const double drift_velocity =
        std::max(0.0, config_.ion_directed_velocity_m_per_s +
                          std::max(0.0, config_.bulk_flow_velocity_m_per_s *
                                            std::clamp(config_.flow_alignment_cosine, -1.0, 1.0)));
    std::mt19937_64 generator(0x71015aULL + static_cast<std::uint64_t>(sample_count));

    double collected_surface_velocity_sum = 0.0;
    for (std::size_t index = 0; index < sample_count; ++index)
    {
        const double launch_velocity = halfNormalFluxAverage(sigma, generator) + drift_velocity;
        double normal_energy_j = 0.5 * ion_mass * launch_velocity * launch_velocity;

        if (surface_potential_v > 0.0)
        {
            const double barrier_j = kElementaryCharge * surface_potential_v;
            if (normal_energy_j <= barrier_j)
            {
                continue;
            }
            normal_energy_j -= barrier_j;
        }
        else
        {
            normal_energy_j += -kElementaryCharge * surface_potential_v;
        }

        const double surface_velocity = std::sqrt(std::max(0.0, 2.0 * normal_energy_j / ion_mass));
        collected_surface_velocity_sum += surface_velocity;
    }

    const double mean_surface_velocity =
        collected_surface_velocity_sum / static_cast<double>(sample_count);
    return 0.5 * kElementaryCharge * config_.plasma.ion_density_m3 * mean_surface_velocity *
           std::max(0.0, config_.ion_collection_coefficient);
}

double DensePlasmaSurfaceCharging::computeEffectiveConductivity(double surface_potential_v) const
{
    const double base_conductivity = std::max(0.0, config_.material.getConductivity());
    const double thickness = std::max(1.0e-8, config_.dielectric_thickness_m);
    const double electric_field = std::abs(surface_potential_v) / thickness;
    const double poole_frenkel_beta =
        std::max(0.0, config_.material.getScalarProperty("poole_frenkel_beta", 0.0));
    const double max_enhancement =
        std::max(1.0, config_.material.getScalarProperty("max_field_enhancement_factor", 1.0e8));
    const double field_enhancement =
        std::min(max_enhancement, safeExp(poole_frenkel_beta * std::sqrt(std::max(0.0, electric_field))));
    return base_conductivity * field_enhancement + computeRadiationInducedConductivity();
}

double DensePlasmaSurfaceCharging::computeLeakageCurrentDensity(double surface_potential_v) const
{
    const double thickness = std::max(1.0e-8, config_.dielectric_thickness_m);
    const double conductivity = computeEffectiveConductivity(surface_potential_v);
    return -conductivity * surface_potential_v / thickness;
}

double DensePlasmaSurfaceCharging::computeNetCurrentDensity(double surface_potential_v) const
{
    const auto currents = computeSurfaceCurrents(surface_potential_v);
    return currents.total_current_a_per_m2 + computeLeakageCurrentDensity(surface_potential_v);
}

double DensePlasmaSurfaceCharging::estimateCurrentDerivative(double surface_potential_v) const
{
    const double step = std::max(1.0e-3, 1.0e-2 * std::max(1.0, std::abs(surface_potential_v)));
    const double jp = computeNetCurrentDensity(surface_potential_v + step);
    const double jm = computeNetCurrentDensity(surface_potential_v - step);
    return (jp - jm) / (2.0 * step);
}

double DensePlasmaSurfaceCharging::computeFloatingPotential() const
{
    double left = -std::max(config_.max_abs_potential_v,
                            20.0 * std::max(1.0, config_.plasma.electron_temperature_ev));
    double right = -left;
    double f_left = computeNetCurrentDensity(left);
    double f_right = computeNetCurrentDensity(right);

    for (int expand = 0; expand < 8 && f_left * f_right > 0.0; ++expand)
    {
        left *= 1.5;
        right *= 1.5;
        f_left = computeNetCurrentDensity(left);
        f_right = computeNetCurrentDensity(right);
    }

    if (f_left * f_right > 0.0)
    {
        return std::abs(f_left) < std::abs(f_right) ? left : right;
    }

    for (int iteration = 0; iteration < 80; ++iteration)
    {
        const double mid = 0.5 * (left + right);
        const double f_mid = computeNetCurrentDensity(mid);
        if (std::abs(f_mid) < 1.0e-12 || std::abs(right - left) < 1.0e-3)
        {
            return mid;
        }

        if (f_left * f_mid <= 0.0)
        {
            right = mid;
            f_right = f_mid;
        }
        else
        {
            left = mid;
            f_left = f_mid;
        }
    }
    return 0.5 * (left + right);
}

double DensePlasmaSurfaceCharging::recommendTimeStep(double remaining_time_s, double minimum_dt_s,
                                                     double maximum_dt_s) const
{
    const double min_dt = std::max(1.0e-12, minimum_dt_s);
    const double max_dt = std::max(min_dt, maximum_dt_s);
    const double remaining = std::max(0.0, remaining_time_s);
    if (!initialized_ || remaining <= 0.0)
    {
        return min_dt;
    }
    if (remaining <= min_dt)
    {
        return remaining;
    }

    const double current_potential = status_.state.surface_potential_v;
    const double target_potential =
        config_.floating ? computeFloatingPotential() : current_potential;
    const double potential_error = std::abs(target_potential - current_potential);
    const double relative_error =
        potential_error / std::max(1.0, std::abs(target_potential));
    const double derivative =
        std::abs(estimateCurrentDerivative(current_potential));
    const double capacitance =
        std::max(1.0e-12, status_.state.capacitance_per_area_f_per_m2);
    const double relaxation_time =
        derivative > 1.0e-18 ? capacitance / derivative : max_dt;
    const double net_current = std::abs(computeNetCurrentDensity(current_potential));
    const double steady_current_threshold =
        1.0e-3 * std::max(1.0, config_.max_abs_current_density_a_per_m2);

    double suggested_dt = min_dt;
    if (relative_error < 5.0e-2 || net_current < steady_current_threshold)
    {
        suggested_dt = max_dt;
    }
    else if (relative_error < 2.5e-1)
    {
        suggested_dt = std::min(max_dt, 5.0 * relaxation_time);
    }
    else
    {
        suggested_dt = std::min(max_dt, std::max(min_dt, 0.5 * relaxation_time));
    }

    return std::clamp(std::min(remaining, suggested_dt), min_dt, max_dt);
}

double DensePlasmaSurfaceCharging::advancePotentialImplicit(double surface_potential_v, double dt) const
{
    const double capacitance = std::max(1.0e-12, status_.state.capacitance_per_area_f_per_m2);
    const double search_limit =
        std::max(config_.max_abs_potential_v, 20.0 * std::max(1.0, config_.plasma.electron_temperature_ev));

    double potential = surface_potential_v;
    for (int iteration = 0; iteration < 20; ++iteration)
    {
        const double residual =
            capacitance * (potential - surface_potential_v) / dt - computeNetCurrentDensity(potential);
        if (std::abs(residual) < 1.0e-9)
        {
            return potential;
        }

        const double derivative =
            capacitance / dt - estimateCurrentDerivative(potential);
        if (std::abs(derivative) < 1.0e-12)
        {
            break;
        }

        const double delta =
            std::clamp(-residual / derivative, -0.2 * search_limit, 0.2 * search_limit);
        potential += delta;
        potential = clampSigned(potential, search_limit);
        if (std::abs(delta) < 1.0e-6)
        {
            return potential;
        }
    }

    const double explicit_candidate =
        surface_potential_v + dt * computeNetCurrentDensity(surface_potential_v) / capacitance;
    return clampSigned(explicit_candidate, search_limit);
}

bool DensePlasmaSurfaceCharging::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }

    const std::size_t substeps = std::max<std::size_t>(1, config_.internal_substeps);
    const double sub_dt = dt / static_cast<double>(substeps);
    auto state = status_.state;
    SurfaceCurrents currents;
    double leakage_current = 0.0;
    double net_current = 0.0;
    double effective_conductivity = computeEffectiveConductivity(state.surface_potential_v);
    const double effective_sheath_length = computeEffectiveSheathLength();

    for (std::size_t step = 0; step < substeps; ++step)
    {
        const double updated_potential = advancePotentialImplicit(state.surface_potential_v, sub_dt);
        currents = computeSurfaceCurrents(updated_potential);
        leakage_current = computeLeakageCurrentDensity(updated_potential);
        net_current = currents.total_current_a_per_m2 + leakage_current;
        effective_conductivity = computeEffectiveConductivity(updated_potential);
        state = accumulation_model_.advance(state, net_current, sub_dt);
        state.surface_potential_v = updated_potential;
        state.surface_charge_density_c_per_m2 =
            state.surface_potential_v * state.capacitance_per_area_f_per_m2;
        state.stored_energy_j_per_m2 =
            0.5 * state.capacitance_per_area_f_per_m2 * state.surface_potential_v *
            state.surface_potential_v;
        if (!std::isfinite(state.surface_potential_v) || !std::isfinite(state.surface_charge_density_c_per_m2))
        {
            return false;
        }
    }

    const double target_potential =
        config_.floating ? computeFloatingPotential() : state.surface_potential_v;
    status_.currents = currents;
    status_.state = state;

    status_.equilibrium_error = status_.state.surface_potential_v - target_potential;
    status_.equilibrium_reached = std::abs(status_.equilibrium_error) < 1.0e-2;
    status_.time_s += dt;
    status_.steps_completed += 1;

    history_time_.push_back(status_.time_s);
    history_potential_.push_back(status_.state.surface_potential_v);
    history_charge_.push_back(status_.state.surface_charge_density_c_per_m2);
    history_current_.push_back(status_.currents.total_current_a_per_m2);
    history_electron_current_.push_back(status_.currents.electron_current_a_per_m2);
    history_ion_current_.push_back(status_.currents.ion_current_a_per_m2);
    history_secondary_current_.push_back(status_.currents.secondary_emission_a_per_m2);
    history_photo_current_.push_back(status_.currents.photo_emission_a_per_m2);
    history_thermionic_current_.push_back(status_.currents.thermionic_emission_a_per_m2);
    history_field_emission_current_.push_back(status_.currents.field_emission_a_per_m2);
    history_leakage_current_.push_back(leakage_current);
    history_net_current_.push_back(net_current);
    history_capacitance_.push_back(status_.state.capacitance_per_area_f_per_m2);
    history_effective_conductivity_.push_back(effective_conductivity);
    history_effective_sheath_length_.push_back(effective_sheath_length);
    history_adaptive_time_step_.push_back(dt);
    history_electron_calibration_factor_.push_back(electron_pic_calibration_factor_);
    history_ion_calibration_factor_.push_back(ion_pic_calibration_factor_);
    history_equilibrium_potential_.push_back(target_potential);
    history_equilibrium_error_.push_back(status_.equilibrium_error);
    return true;
}

void DensePlasmaSurfaceCharging::reset()
{
    *this = DensePlasmaSurfaceCharging{};
}

bool DensePlasmaSurfaceCharging::exportResults(const std::filesystem::path& csv_path) const
{
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "time_s";
    data_set.axis_values = history_time_;
    data_set.scalar_series["surface_potential_v"] = history_potential_;
    data_set.scalar_series["surface_charge_density_c_per_m2"] = history_charge_;
    data_set.scalar_series["total_current_density_a_per_m2"] = history_current_;
    data_set.scalar_series["electron_current_density_a_per_m2"] = history_electron_current_;
    data_set.scalar_series["ion_current_density_a_per_m2"] = history_ion_current_;
    data_set.scalar_series["secondary_emission_density_a_per_m2"] = history_secondary_current_;
    data_set.scalar_series["photo_emission_density_a_per_m2"] = history_photo_current_;
    data_set.scalar_series["thermionic_emission_density_a_per_m2"] = history_thermionic_current_;
    data_set.scalar_series["field_emission_density_a_per_m2"] = history_field_emission_current_;
    data_set.scalar_series["leakage_current_density_a_per_m2"] = history_leakage_current_;
    data_set.scalar_series["net_current_density_a_per_m2"] = history_net_current_;
    data_set.scalar_series["capacitance_per_area_f_per_m2"] = history_capacitance_;
    data_set.scalar_series["effective_conductivity_s_per_m"] = history_effective_conductivity_;
    data_set.scalar_series["effective_sheath_length_m"] = history_effective_sheath_length_;
    data_set.scalar_series["adaptive_time_step_s"] = history_adaptive_time_step_;
    data_set.scalar_series["electron_pic_calibration_factor"] = history_electron_calibration_factor_;
    data_set.scalar_series["ion_pic_calibration_factor"] = history_ion_calibration_factor_;
    data_set.scalar_series["floating_equilibrium_potential_v"] = history_equilibrium_potential_;
    data_set.scalar_series["equilibrium_error_v"] = history_equilibrium_error_;
    data_set.metadata["module"] = "Surface Charging";
    return static_cast<bool>(exporter_.exportDataSet(csv_path, data_set));
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
