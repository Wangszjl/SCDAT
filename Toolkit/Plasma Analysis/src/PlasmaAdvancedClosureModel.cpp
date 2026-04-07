#include "PlasmaAdvancedClosureModel.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{
constexpr double kMinDensity = 1.0e6;
constexpr double kMinTemperatureEv = 0.05;

double safeNonEquilibriumRatio(double electron_temperature_ev, double ion_temperature_ev)
{
    const double te = std::max(kMinTemperatureEv, electron_temperature_ev);
    const double ti = std::max(kMinTemperatureEv, ion_temperature_ev);
    return std::max(te / ti, ti / te);
}

double computeVelocityScale(double electric_field_v_per_m)
{
    return 0.3 * std::sqrt(std::max(0.0, std::abs(electric_field_v_per_m)));
}

} // namespace

void PlasmaAdvancedClosureModel::configure(const AdvancedClosureConfig& config)
{
    config_ = config;
}

AdvancedClosureState PlasmaAdvancedClosureModel::evaluate(
    const PlasmaParameters& parameters,
    const DensePlasmaAssessment& assessment,
    const ReactionCollisionState& reaction_state,
    double characteristic_length_m,
    double dt) const
{
    AdvancedClosureState state;
    state.non_equilibrium_ratio =
        safeNonEquilibriumRatio(parameters.electron_temperature_ev, parameters.ion_temperature_ev);

    const bool non_equilibrium_enabled =
        config_.enable_non_equilibrium_closure &&
        config_.non_equilibrium_model != NonEquilibriumClosureModel::Disabled;
    const bool turbulence_enabled =
        config_.enable_turbulence_closure &&
        config_.turbulence_model != TurbulenceClosureModel::Disabled;

    if (non_equilibrium_enabled)
    {
        const double ratio_excess =
            std::clamp(state.non_equilibrium_ratio - 1.0, 0.0, config_.non_equilibrium_ratio_cap);
        const double excess_factor = ratio_excess / std::max(1.0, config_.non_equilibrium_ratio_cap);

        const double collisional_base_rate =
            2.0e5 * std::max(0.0, assessment.collisionality) +
            1.0e-2 * std::max(0.0, reaction_state.effective_collision_frequency_hz);

        if (config_.non_equilibrium_model == NonEquilibriumClosureModel::TwoTemperatureRelaxation)
        {
            state.non_equilibrium_relaxation_rate_hz =
                config_.non_equilibrium_relaxation_gain * collisional_base_rate * excess_factor;
        }
        else if (config_.non_equilibrium_model ==
                 NonEquilibriumClosureModel::CollisionalThermalization)
        {
            state.non_equilibrium_relaxation_rate_hz =
                config_.non_equilibrium_relaxation_gain * std::sqrt(std::max(0.0, collisional_base_rate)) *
                (2.0e3 * excess_factor);
        }

        const double temperature_gap_ev =
            parameters.electron_temperature_ev - parameters.ion_temperature_ev;
        const double proposed_delta_ev =
            -temperature_gap_ev * state.non_equilibrium_relaxation_rate_hz * dt;
        state.electron_temperature_delta_ev =
            std::clamp(proposed_delta_ev,
                       -std::max(0.0, config_.max_temperature_correction_ev_per_step),
                       std::max(0.0, config_.max_temperature_correction_ev_per_step));

        if (std::abs(state.electron_temperature_delta_ev) > 0.0)
        {
            state.active_terms += 1;
        }
    }

    if (turbulence_enabled)
    {
        const double characteristic_length = std::max(1.0e-6, characteristic_length_m);
        const double velocity_scale = computeVelocityScale(parameters.electric_field_v_per_m);
        const double collisionality_weight =
            std::clamp(assessment.collisionality / (1.0 + std::max(0.0, assessment.collisionality)),
                       0.0,
                       1.0);

        state.turbulence_intensity =
            std::clamp(config_.turbulence_gain * (0.4 + 0.6 * collisionality_weight) *
                           velocity_scale / 8.0,
                       0.0,
                       2.0);

        if (config_.turbulence_model == TurbulenceClosureModel::MixingLengthEddyDiffusivity)
        {
            state.turbulence_eddy_diffusivity_m2_per_s =
                config_.turbulence_mixing_length_m * velocity_scale * state.turbulence_intensity;
        }
        else if (config_.turbulence_model == TurbulenceClosureModel::CollisionalDamping)
        {
            state.turbulence_eddy_diffusivity_m2_per_s =
                0.8 * config_.turbulence_mixing_length_m * velocity_scale *
                state.turbulence_intensity * (1.0 + 0.5 * collisionality_weight);
        }

        const double electric_energy_density =
            8.8541878128e-12 * parameters.electric_field_v_per_m * parameters.electric_field_v_per_m;
        state.turbulence_dissipation_rate_w_per_m3 =
            state.turbulence_eddy_diffusivity_m2_per_s * electric_energy_density /
            std::max(1.0e-6, characteristic_length * characteristic_length);

        const double density_fraction = std::clamp(
            state.turbulence_intensity * dt * 2.0e5,
            0.0,
            std::max(0.0, config_.max_relative_density_correction_per_step));
        state.density_correction_m3 =
            -density_fraction * std::max(kMinDensity, parameters.electron_density_m3);

        if (state.turbulence_eddy_diffusivity_m2_per_s > 0.0)
        {
            state.active_terms += 1;
        }
    }

    return state;
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
