#pragma once

#include "../../Tools/Geometry/include/Vector3D.h"

#include <cstddef>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct PlasmaParameters
{
    double electron_density_m3 = 1.0e17;
    double ion_density_m3 = 1.0e17;
    double electron_temperature_ev = 5.0;
    double ion_temperature_ev = 0.2;
    double ion_mass_amu = 40.0;
    double neutral_density_m3 = 2.5e20;
    double electric_field_v_per_m = 0.0;
    double pressure_pa = 1.0;
};

struct ReactionCollisionConfig
{
    bool enable_electron_impact_ionization = true;
    bool enable_radiative_recombination = true;
    bool enable_charge_exchange = true;
    bool enable_elastic_momentum_transfer = true;
    double ionization_threshold_ev = 15.8;
    double ionization_cross_section_m2 = 2.5e-20;
    double charge_exchange_cross_section_m2 = 2.0e-19;
    double elastic_momentum_cross_section_m2 = 1.5e-19;
    double max_relative_density_change_per_step = 0.08;
};

enum class NonEquilibriumClosureModel
{
    Disabled = 0,
    TwoTemperatureRelaxation = 1,
    CollisionalThermalization = 2,
};

enum class TurbulenceClosureModel
{
    Disabled = 0,
    MixingLengthEddyDiffusivity = 1,
    CollisionalDamping = 2,
};

struct AdvancedClosureConfig
{
    bool enable_non_equilibrium_closure = false;
    bool enable_turbulence_closure = false;
    NonEquilibriumClosureModel non_equilibrium_model =
        NonEquilibriumClosureModel::TwoTemperatureRelaxation;
    TurbulenceClosureModel turbulence_model =
        TurbulenceClosureModel::MixingLengthEddyDiffusivity;
    double non_equilibrium_relaxation_gain = 0.35;
    double non_equilibrium_ratio_cap = 3.5;
    double turbulence_mixing_length_m = 2.0e-4;
    double turbulence_gain = 0.6;
    double max_temperature_correction_ev_per_step = 0.4;
    double max_relative_density_correction_per_step = 0.05;
};

struct FluidAlgorithmConfig
{
    Geometry::Vector3D domain_size{1.0e-3, 1.0e-3, 5.0e-3};
    Geometry::Vector3D resolution{6.0, 6.0, 24.0};
    double time_step_s = 1.0e-9;
    double dense_plasma_threshold_m3 = 1.0e18;
    double initial_potential_v = 25.0;
    PlasmaParameters initial_plasma;
    ReactionCollisionConfig reaction_collision;
    AdvancedClosureConfig advanced_closure;
};

struct DensePlasmaAssessment
{
    bool is_dense = false;
    double debye_length_m = 0.0;
    double plasma_frequency_hz = 0.0;
    double collisionality = 0.0;
};

struct FluidAlgorithmStatus
{
    double time_s = 0.0;
    bool dense_plasma_detected = false;
    double average_potential_v = 0.0;
    double average_density_m3 = 0.0;
    double sheath_thickness_m = 0.0;
    double presheath_thickness_m = 0.0;
    double sheath_potential_drop_v = 0.0;
    double presheath_potential_drop_v = 0.0;
    double boundary_wall_field_v_per_m = 0.0;
    double ion_bohm_velocity_m_per_s = 0.0;
    double ion_mach_at_sheath_edge = 1.0;
    double ionization_source_m3_per_s = 0.0;
    double recombination_sink_m3_per_s = 0.0;
    double effective_collision_frequency_hz = 0.0;
    double charge_exchange_frequency_hz = 0.0;
    double reaction_energy_loss_ev_per_s = 0.0;
    double reaction_momentum_transfer_ratio = 0.0;
    std::size_t reaction_active_processes = 0;
    bool advanced_closure_enabled = false;
    double non_equilibrium_ratio = 1.0;
    double non_equilibrium_relaxation_rate_hz = 0.0;
    double electron_temperature_closure_delta_ev = 0.0;
    double turbulence_intensity = 0.0;
    double turbulence_eddy_diffusivity_m2_per_s = 0.0;
    double turbulence_dissipation_rate_w_per_m3 = 0.0;
    double closure_density_correction_m3 = 0.0;
    std::size_t closure_active_terms = 0;
    double debye_length_m = 0.0;
    std::size_t steps_completed = 0;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
