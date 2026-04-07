#include "PlasmaReactionCollisionLibrary.h"

#include <algorithm>
#include <cmath>
#include <memory>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kAtomicMassUnit = 1.66053906660e-27;
constexpr double kPi = 3.14159265358979323846;

inline double electronThermalVelocity(double electron_temperature_ev)
{
    const double temperature_j = std::max(0.05, electron_temperature_ev) * kElementaryCharge;
    return std::sqrt(8.0 * temperature_j / (kPi * kElectronMass));
}

inline double ionThermalVelocity(double ion_temperature_ev, double ion_mass_amu)
{
    const double ion_mass = std::max(1.0, ion_mass_amu) * kAtomicMassUnit;
    const double temperature_j = std::max(0.02, ion_temperature_ev) * kElementaryCharge;
    return std::sqrt(8.0 * temperature_j / (kPi * ion_mass));
}

class ElectronImpactIonizationProcess final : public IReactionCollisionProcess
{
  public:
    ElectronImpactIonizationProcess(double threshold_ev, double cross_section_m2)
        : threshold_ev_(std::max(1.0, threshold_ev)),
          cross_section_m2_(std::max(1.0e-22, cross_section_m2))
    {
    }

    void evaluate(const PlasmaParameters& parameters,
                  const DensePlasmaAssessment&,
                  ReactionCollisionState& state) const override
    {
        const double ne = std::max(1.0e6, parameters.electron_density_m3);
        const double nn = std::max(1.0e6, parameters.neutral_density_m3);
        const double te = std::max(0.05, parameters.electron_temperature_ev);

        const double velocity = electronThermalVelocity(te);
        const double activation = std::exp(-threshold_ev_ / te);
        const double rate_coefficient = cross_section_m2_ * velocity * activation;
        const double source = ne * nn * rate_coefficient;

        state.ionization_source_m3_per_s += source;
        state.effective_collision_frequency_hz += nn * cross_section_m2_ * velocity;
        state.electron_energy_loss_ev_per_s += (source / ne) * (threshold_ev_ + 1.5) * 0.05;
    }

  private:
    double threshold_ev_ = 15.8;
    double cross_section_m2_ = 2.5e-20;
};

class RadiativeRecombinationProcess final : public IReactionCollisionProcess
{
  public:
    void evaluate(const PlasmaParameters& parameters,
                  const DensePlasmaAssessment&,
                  ReactionCollisionState& state) const override
    {
        const double ne = std::max(1.0e6, parameters.electron_density_m3);
        const double ni = std::max(1.0e6, parameters.ion_density_m3);
        const double te = std::max(0.05, parameters.electron_temperature_ev);

        const double rate_coefficient = 2.5e-19 / std::sqrt(te);
        const double sink = ne * ni * rate_coefficient;

        state.recombination_sink_m3_per_s += sink;
        state.electron_energy_loss_ev_per_s += (sink / ne) * 0.5;
    }
};

class ChargeExchangeProcess final : public IReactionCollisionProcess
{
  public:
    explicit ChargeExchangeProcess(double cross_section_m2)
        : cross_section_m2_(std::max(1.0e-22, cross_section_m2))
    {
    }

    void evaluate(const PlasmaParameters& parameters,
                  const DensePlasmaAssessment& assessment,
                  ReactionCollisionState& state) const override
    {
        const double nn = std::max(1.0e6, parameters.neutral_density_m3);
        const double ion_velocity =
            ionThermalVelocity(parameters.ion_temperature_ev, parameters.ion_mass_amu);

        const double nu_cx = nn * cross_section_m2_ * ion_velocity;
        state.charge_exchange_frequency_hz += nu_cx;
        state.effective_collision_frequency_hz += nu_cx;

        const double coupling =
            nu_cx / (nu_cx + std::max(1.0, assessment.plasma_frequency_hz));
        state.momentum_transfer_ratio += 0.65 * coupling;
    }

  private:
    double cross_section_m2_ = 2.0e-19;
};

class ElasticMomentumTransferProcess final : public IReactionCollisionProcess
{
  public:
    explicit ElasticMomentumTransferProcess(double cross_section_m2)
        : cross_section_m2_(std::max(1.0e-22, cross_section_m2))
    {
    }

    void evaluate(const PlasmaParameters& parameters,
                  const DensePlasmaAssessment& assessment,
                  ReactionCollisionState& state) const override
    {
        const double nn = std::max(1.0e6, parameters.neutral_density_m3);
        const double effective_ion_temperature = std::max(
            parameters.ion_temperature_ev,
            0.2 * std::max(0.05, parameters.electron_temperature_ev));
        const double ion_velocity =
            ionThermalVelocity(effective_ion_temperature, parameters.ion_mass_amu);

        const double nu_momentum = nn * cross_section_m2_ * ion_velocity;
        state.effective_collision_frequency_hz += nu_momentum;

        const double collisionality_weight =
            assessment.collisionality / (1.0 + std::max(0.0, assessment.collisionality));
        state.momentum_transfer_ratio += 0.25 * std::clamp(collisionality_weight, 0.0, 1.0);
    }

  private:
    double cross_section_m2_ = 1.5e-19;
};

} // namespace

PlasmaReactionCollisionLibrary::PlasmaReactionCollisionLibrary()
{
    rebuildBuiltinProcesses();
}

void PlasmaReactionCollisionLibrary::configure(const ReactionCollisionConfig& config)
{
    config_ = config;
    rebuildBuiltinProcesses();
}

void PlasmaReactionCollisionLibrary::registerProcess(
    std::shared_ptr<const IReactionCollisionProcess> process)
{
    if (process)
    {
        custom_processes_.push_back(std::move(process));
    }
}

void PlasmaReactionCollisionLibrary::clearCustomProcesses()
{
    custom_processes_.clear();
}

ReactionCollisionState PlasmaReactionCollisionLibrary::evaluate(
    const PlasmaParameters& parameters,
    const DensePlasmaAssessment& assessment) const
{
    ReactionCollisionState state;

    for (const auto& process : builtin_processes_)
    {
        process->evaluate(parameters, assessment, state);
        state.active_processes += 1;
    }

    for (const auto& process : custom_processes_)
    {
        process->evaluate(parameters, assessment, state);
        state.active_processes += 1;
    }

    state.ionization_source_m3_per_s = std::max(0.0, state.ionization_source_m3_per_s);
    state.recombination_sink_m3_per_s = std::max(0.0, state.recombination_sink_m3_per_s);
    state.effective_collision_frequency_hz =
        std::max(0.0, state.effective_collision_frequency_hz);
    state.charge_exchange_frequency_hz = std::max(0.0, state.charge_exchange_frequency_hz);
    state.electron_energy_loss_ev_per_s = std::max(0.0, state.electron_energy_loss_ev_per_s);
    state.momentum_transfer_ratio = std::clamp(state.momentum_transfer_ratio, 0.0, 1.0);
    return state;
}

void PlasmaReactionCollisionLibrary::rebuildBuiltinProcesses()
{
    builtin_processes_.clear();

    if (config_.enable_electron_impact_ionization)
    {
        builtin_processes_.push_back(std::make_shared<ElectronImpactIonizationProcess>(
            config_.ionization_threshold_ev,
            config_.ionization_cross_section_m2));
    }

    if (config_.enable_radiative_recombination)
    {
        builtin_processes_.push_back(std::make_shared<RadiativeRecombinationProcess>());
    }

    if (config_.enable_charge_exchange)
    {
        builtin_processes_.push_back(std::make_shared<ChargeExchangeProcess>(
            config_.charge_exchange_cross_section_m2));
    }

    if (config_.enable_elastic_momentum_transfer)
    {
        builtin_processes_.push_back(std::make_shared<ElasticMomentumTransferProcess>(
            config_.elastic_momentum_cross_section_m2));
    }
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
