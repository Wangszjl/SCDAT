#include "ReferenceCurrentBalanceModel.h"
#include "SurfaceFlowCouplingModel.h"
#include "../../Tools/FieldSolver/include/SurfaceBarrierModels.h"
#include "../../Tools/Material/include/SurfaceMaterialModel.h"
#include "../../Tools/Particle/include/SurfaceDistributionFunction.h"

#include <algorithm>
#include <cmath>
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

double bodyElectronCollectionScale(const Material::MaterialProperty& material)
{
    return clampPositive(material.getScalarProperty("body_electron_collection_scale", 1.0), 0.0);
}

double resolveBodyElectronCollectionCoefficient(double base_coefficient,
                                                const Material::MaterialProperty& material)
{
    return clampPositive(
        material.getScalarProperty("body_electron_collection_coefficient", base_coefficient), 0.0) *
           bodyElectronCollectionScale(material);
}

double bodyPhotoEmissionScale(const Material::MaterialProperty& material)
{
    return clampPositive(material.getScalarProperty("body_photo_emission_scale", 1.0), 0.0);
}

double bodyPhotoelectronTemperatureEv(const Material::MaterialProperty& material,
                                      double fallback_ev)
{
    return clampPositive(
        material.getScalarProperty("body_photoelectron_temperature_ev", fallback_ev), 1.0e-3);
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

template <typename PopulationFn, typename FluxFn>
double evaluateSpectrum(const Particle::ResolvedSpectrum& spectrum, PopulationFn population_fn,
                        FluxFn flux_fn)
{
    if (Particle::usesResolvedSpectrumDiscreteFlux(spectrum))
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

double ReferenceCurrentBalanceModel::ionSecondaryYield(
    const Material::MaterialProperty& material, double incident_energy_ev,
    double surface_potential_v)
{
    return Material::evaluateLegacyIonSecondaryElectronYield(
        material, std::max(0.0, incident_energy_ev - surface_potential_v));
}

double ReferenceCurrentBalanceModel::backscatterYield(
    const Material::MaterialProperty& material, double incident_energy_ev,
    double surface_potential_v)
{
    return Material::evaluateLegacyBackscatterYield(
        material, std::max(0.0, incident_energy_ev + surface_potential_v));
}

double ReferenceCurrentBalanceModel::bodyPhotoCurrentDensity(double base_photo_current_density_a_per_m2,
                                                             double emission_scale)
{
    return 0.25 * clampPositive(emission_scale, 0.0) *
           std::max(0.0, base_photo_current_density_a_per_m2);
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
    FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = surface_potential_v;
    state.reference_potential_v = 0.0;
    state.barrier_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 0.0;
    state.emission_temperature_ev = clampPositive(characteristic_energy_ev, 1.0e-3);

    FieldSolver::VariableBarrierScaler scaler;
    const auto evaluation = scaler.evaluate(config_.patch_material, state, 1.0);
    return evaluation.valid ? std::clamp(evaluation.scaling, 0.0, 1.0) : 0.0;
}

double ReferenceCurrentBalanceModel::positivePotentialIonBarrierScale(
    const Material::MaterialProperty& material, double surface_potential_v,
    double characteristic_energy_ev) const
{
    if (surface_potential_v <= 0.0)
    {
        return 1.0;
    }

    FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = surface_potential_v;
    state.reference_potential_v = 0.0;
    state.barrier_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 0.0;
    state.emission_temperature_ev = clampPositive(characteristic_energy_ev, 1.0e-3);

    FieldSolver::VariableBarrierScaler scaler;
    const auto evaluation = scaler.evaluate(material, state, 1.0);
    return evaluation.valid ? std::clamp(evaluation.scaling, 0.0, 1.0) : 0.0;
}

double ReferenceCurrentBalanceModel::ionSecondaryRecollectionScale(
    const Material::MaterialProperty& material, double surface_potential_v,
    double characteristic_energy_ev) const
{
    if (surface_potential_v <= 0.0)
    {
        return 1.0;
    }

    FieldSolver::SurfaceBarrierState state;
    state.local_potential_v = surface_potential_v;
    state.reference_potential_v = 0.0;
    state.barrier_potential_v = 0.0;
    state.normal_electric_field_v_per_m = 0.0;
    state.emission_temperature_ev = clampPositive(characteristic_energy_ev, 1.0e-3);

    FieldSolver::SecondaryRecollectionScaler scaler;
    const auto evaluation = scaler.evaluate(material, state, 1.0);
    return evaluation.valid ? std::clamp(evaluation.scaling, 0.0, 1.0) : 0.0;
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
    Material::LegacySecondaryYieldModel model =
        Material::LegacySecondaryYieldModel::Whipple;
    switch (config_.see_model)
    {
    case SecondaryElectronEmissionModel::Whipple:
        model = Material::LegacySecondaryYieldModel::Whipple;
        break;
    case SecondaryElectronEmissionModel::Sims:
        model = Material::LegacySecondaryYieldModel::Sims;
        break;
    case SecondaryElectronEmissionModel::Katz:
        model = Material::LegacySecondaryYieldModel::Katz;
        break;
    }
    return Material::evaluateLegacySecondaryElectronYield(
        material, model, std::max(0.0, incident_energy_ev + surface_potential_v),
        material.getScalarProperty("sims_exponent_n", 1.6));
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
                                                            double projected_flow_speed_m_per_s)
{
    const double ni_cm3 = clampPositive(ion_density_m3 * 1.0e-6, 0.0);
    const double ti_ev = clampPositive(ion_temperature_ev, 1.0e-6);
    const double flow_speed_km_s = projected_flow_speed_m_per_s * 1.0e-3;
    const double qv = std::sqrt(std::abs(surface_potential_v) / ti_ev);
    const double wt = 13.84 * std::sqrt(ti_ev / clampPositive(ion_mass_amu, 1.0));
    const double qd = flow_speed_km_s / clampPositive(wt, 1.0e-12);
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

double ReferenceCurrentBalanceModel::computeElectronCollection(double surface_potential_v,
                                                              double collection_coefficient) const
{
    const auto& collection_model =
        selectSurfaceElectronCollectionModel(config_.electron_collection_model);
    const auto population_current = [&](const Particle::SpectrumPopulation& population) {
        const double base_flux = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor =
            collection_model.populationFactor(surface_potential_v, population.temperature_ev);
        return -base_flux * clampPositive(collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * factor;
    };
    const auto flux_current = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double mean_energy_ev = Particle::resolvedSpectrumAverageEnergyEv(spectrum, 1.0);
        const double current = Particle::integrateResolvedSpectrumFlux(
            spectrum, [&](double energy_ev, double flux) {
            return flux * collection_model.fluxCollectionWeight(surface_potential_v, energy_ev,
                                                                mean_energy_ev);
        });
        return -kElementaryCharge * clampPositive(collection_coefficient, 0.0) *
               clampPositive(config_.electron_calibration_factor, 0.1) * current;
    };
    if (config_.has_electron_spectrum)
    {
        return evaluateSpectrum(config_.electron_spectrum, population_current, flux_current);
    }
    const double base_flux = thermalCurrentDensity(
        config_.plasma.electron_density_m3, config_.plasma.electron_temperature_ev, kElectronMass);
    const double factor = collection_model.populationFactor(
        surface_potential_v, config_.plasma.electron_temperature_ev);
    return -base_flux * clampPositive(collection_coefficient, 0.0) *
           clampPositive(config_.electron_calibration_factor, 0.1) * factor;
}

double ReferenceCurrentBalanceModel::computeIonCollection(double surface_potential_v,
                                                          const Material::MaterialProperty& material) const
{
    const auto population_current = [&](const Particle::SpectrumPopulation& population) {
        const double base_flux = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = surface_potential_v >= 0.0
                                  ? positivePotentialIonBarrierScale(
                                        material, surface_potential_v, population.temperature_ev)
                                  : (1.0 - surface_potential_v /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return base_flux * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * factor;
    };
    const auto flux_current = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double mean_energy_ev = Particle::resolvedSpectrumAverageEnergyEv(spectrum, 1.0);
        const double current = Particle::integrateResolvedSpectrumFlux(
            spectrum, [&](double energy_ev, double flux) {
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
                              ? positivePotentialIonBarrierScale(
                                    material, surface_potential_v, config_.plasma.ion_temperature_ev)
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
    const auto& collection_model =
        selectSurfaceElectronCollectionModel(config_.electron_collection_model);
    const double electron_collection_coefficient =
        role == ReferenceSurfaceRole::Body
            ? resolveBodyElectronCollectionCoefficient(
                  clampPositive(config_.electron_collection_coefficient, 0.0), material)
            : clampPositive(config_.electron_collection_coefficient, 0.0);
    const double photo_emission_scale =
        role == ReferenceSurfaceRole::Body ? bodyPhotoEmissionScale(material) : 1.0;
    const double photoelectron_temperature_ev =
        role == ReferenceSurfaceRole::Body
            ? bodyPhotoelectronTemperatureEv(material, config_.photoelectron_temperature_ev)
            : clampPositive(config_.photoelectron_temperature_ev, 1.0e-3);
    components.electron_collection_a_per_m2 =
        computeElectronCollection(potential, electron_collection_coefficient);
    components.ion_collection_a_per_m2 = computeIonCollection(potential, material);

    const double secondary_escape = emissionEscapeProbability(
        potential, material.getScalarProperty("secondary_emission_escape_energy_ev", 2.0));
    const double photo_escape = config_.use_photoelectron_suppression
                                    ? emissionEscapeProbability(potential,
                                                                photoelectron_temperature_ev)
                                    : 1.0;

    const auto electron_secondary_from_pop = [&](const Particle::SpectrumPopulation& population) {
        const double average_yield = Particle::integratePopulationMaxwellianAverageEv(
            population.temperature_ev, [&](double energy_ev) {
                return computeSeeYield(material, energy_ev + std::max(0.0, potential), 0.0);
            });
        const double current = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = collection_model.populationFactor(potential, population.temperature_ev);
        return current * electron_collection_coefficient *
               clampPositive(config_.electron_calibration_factor, 0.1) * factor * average_yield *
               secondary_escape;
    };
    const auto electron_secondary_from_flux = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double mean_energy_ev = Particle::resolvedSpectrumAverageEnergyEv(spectrum, 1.0);
        const double flux = Particle::integrateResolvedSpectrumFlux(
            spectrum, [&](double energy_ev, double value) {
            const double weight =
                collection_model.fluxEmissionWeight(potential, energy_ev, mean_energy_ev);
            return value * weight *
                   computeSeeYield(material, energy_ev + std::max(0.0, potential), 0.0);
        });
        return flux * kElementaryCharge * electron_collection_coefficient *
               clampPositive(config_.electron_calibration_factor, 0.1) * secondary_escape;
    };
    const auto ion_secondary_from_pop = [&](const Particle::SpectrumPopulation& population) {
        const double average_yield = Particle::integratePopulationMaxwellianAverageEv(
            population.temperature_ev, [&](double energy_ev) {
                return ionSecondaryYield(material, energy_ev + std::max(0.0, -potential), 0.0);
            });
        const double current = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = potential >= 0.0
                                  ? ionSecondaryRecollectionScale(
                                        material, potential, population.temperature_ev)
                                  : (1.0 - potential /
                                               clampPositive(population.temperature_ev, 1.0e-3));
        return current * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * factor * average_yield *
               secondary_escape;
    };
    const auto ion_secondary_from_flux = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double flux = Particle::integrateResolvedSpectrumFlux(
            spectrum, [&](double energy_ev, double value) {
            if (potential > 0.0 && energy_ev - potential <= 0.0)
            {
                return 0.0;
            }
            return value * ionSecondaryYield(material,
                                             energy_ev + std::max(0.0, -potential), 0.0);
        });
        return flux * kElementaryCharge * clampPositive(config_.ion_collection_coefficient, 0.0) *
               clampPositive(config_.ion_calibration_factor, 0.1) * secondary_escape;
    };
    const auto backscatter_from_pop = [&](const Particle::SpectrumPopulation& population) {
        const double average_yield = Particle::integratePopulationMaxwellianAverageEv(
            population.temperature_ev, [&](double energy_ev) {
                return backscatterYield(material, energy_ev + std::max(0.0, potential), 0.0);
            });
        const double current = thermalCurrentDensity(
            population.density_m3, population.temperature_ev,
            clampPositive(population.mass_amu, 1.0e-6) * kAtomicMassUnit);
        const double factor = collection_model.populationFactor(potential, population.temperature_ev);
        return current * electron_collection_coefficient *
               clampPositive(config_.electron_calibration_factor, 0.1) * factor * average_yield *
               secondary_escape;
    };
    const auto backscatter_from_flux = [&](const Particle::ResolvedSpectrum& spectrum) {
        const double mean_energy_ev = Particle::resolvedSpectrumAverageEnergyEv(spectrum, 1.0);
        const double flux = Particle::integrateResolvedSpectrumFlux(
            spectrum, [&](double energy_ev, double value) {
            const double weight =
                collection_model.fluxEmissionWeight(potential, energy_ev, mean_energy_ev);
            return value * weight * backscatterYield(
                               material, energy_ev + std::max(0.0, potential), 0.0);
        });
        return flux * kElementaryCharge * electron_collection_coefficient *
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
        const double see_yield = Particle::integratePopulationMaxwellianAverageEv(
            config_.plasma.electron_temperature_ev, [&](double energy_ev) {
                return computeSeeYield(material, energy_ev + std::max(0.0, potential), 0.0);
            });
        const double backscatter_yield = Particle::integratePopulationMaxwellianAverageEv(
            config_.plasma.electron_temperature_ev, [&](double energy_ev) {
                return backscatterYield(material, energy_ev + std::max(0.0, potential), 0.0);
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
        const double ion_secondary_yield = Particle::integratePopulationMaxwellianAverageEv(
            config_.plasma.ion_temperature_ev, [&](double energy_ev) {
                return ionSecondaryYield(material, energy_ev + std::max(0.0, -potential),
                                         0.0);
            });
        components.ion_secondary_electron_a_per_m2 =
            std::abs(components.ion_collection_a_per_m2) * ion_secondary_yield * secondary_escape;
    }

    const double base_photo_current =
        role == ReferenceSurfaceRole::Body
            ? bodyPhotoCurrentDensity(config_.body_photo_current_density_a_per_m2,
                                      photo_emission_scale)
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
        const auto flow_projection = resolveSurfaceFlowProjection(
            config_.bulk_flow_velocity_m_per_s, config_.flow_alignment_cosine,
            config_.patch_flow_angle_deg);
        const double ion_density = config_.has_ion_spectrum
                                       ? Particle::resolvedSpectrumDensity(
                                             config_.ion_spectrum,
                                             config_.plasma.ion_density_m3)
                                       : config_.plasma.ion_density_m3;
        const double ion_temperature =
            config_.has_ion_spectrum
                ? Particle::resolvedSpectrumCharacteristicEnergyEv(
                      config_.ion_spectrum,
                      clampPositive(config_.plasma.ion_temperature_ev, 1.0))
                : config_.plasma.ion_temperature_ev;
        const double ion_mass_amu =
            config_.has_ion_spectrum
                ? Particle::resolvedSpectrumAverageMassAmu(
                      config_.ion_spectrum,
                      clampPositive(config_.plasma.ion_mass_amu, 1.0))
                : config_.plasma.ion_mass_amu;
        if (role == ReferenceSurfaceRole::Body)
        {
            components.ram_ion_a_per_m2 = ramBodyCurrentDensity(
                potential, ion_density, ion_temperature, ion_mass_amu,
                std::max(0.0, flow_projection.body_projected_speed_m_per_s));
        }
        else
        {
            components.ram_ion_a_per_m2 = ramPatchCurrentDensity(
                potential, ion_density, ion_temperature, ion_mass_amu,
                flow_projection.patch_projected_speed_m_per_s);
        }
    }

    if (!config_.enable_secondary_electron)
    {
        components.secondary_electron_a_per_m2 = 0.0;
        components.ion_secondary_electron_a_per_m2 = 0.0;
    }
    if (!config_.enable_backscatter)
    {
        components.backscatter_electron_a_per_m2 = 0.0;
    }
    if (!config_.enable_photoelectron)
    {
        components.photoelectron_a_per_m2 = 0.0;
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
                                         ? Particle::resolvedSpectrumCharacteristicEnergyEv(
                                               config_.electron_spectrum,
                                               clampPositive(config_.plasma
                                                                 .electron_temperature_ev,
                                                             1.0))
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
                                         ? Particle::resolvedSpectrumCharacteristicEnergyEv(
                                               config_.electron_spectrum,
                                               clampPositive(config_.plasma
                                                                 .electron_temperature_ev,
                                                             1.0))
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
