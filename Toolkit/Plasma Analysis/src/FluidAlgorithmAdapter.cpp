#include "FluidAlgorithmAdapter.h"

#include <algorithm>
#include <numeric>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{

const char* environmentModelName(PlasmaEnvironmentModelKind model)
{
    switch (model)
    {
    case PlasmaEnvironmentModelKind::SpisCcpReference:
        return "spis_ccp_reference";
    case PlasmaEnvironmentModelKind::SpisOrbitalWake:
        return "spis_orbital_wake";
    case PlasmaEnvironmentModelKind::SpisThrusterPlume:
    default:
        return "spis_thruster_plume";
    }
}

const char* distributionModelName(PlasmaDistributionModelKind model)
{
    switch (model)
    {
    case PlasmaDistributionModelKind::MaxwellianProjected:
        return "maxwellian_projected";
    case PlasmaDistributionModelKind::WakeAnisotropic:
        return "wake_anisotropic";
    case PlasmaDistributionModelKind::MultiPopulationHybrid:
    default:
        return "multi_population_hybrid";
    }
}

const char* reactionRegistryName(PlasmaReactionRegistryKind registry)
{
    switch (registry)
    {
    case PlasmaReactionRegistryKind::SpisDensePlasma:
        return "spis_dense_plasma";
    case PlasmaReactionRegistryKind::SpisCore:
    default:
        return "spis_core";
    }
}

const char* diagnosticSetName(PlasmaDiagnosticSetKind diagnostic_set)
{
    switch (diagnostic_set)
    {
    case PlasmaDiagnosticSetKind::SpisCoreDiagnostics:
        return "spis_core_diagnostics";
    case PlasmaDiagnosticSetKind::SheathMultiscaleDiagnostics:
        return "sheath_multiscale_diagnostics";
    case PlasmaDiagnosticSetKind::FullPhysicsDiagnostics:
    default:
        return "full_physics_diagnostics";
    }
}

} // namespace

bool FluidAlgorithmAdapter::initialize(const FluidAlgorithmConfig& config)
{
    config_ = config;
    plasma_state_ = config.initial_plasma;
    detector_ = DensePlasmaDetector(config.dense_plasma_threshold_m3);
    reaction_collision_library_.configure(config.reaction_collision);
    advanced_closure_model_.configure(config.advanced_closure);

    FieldSolver::BoltzmannParameters boltzmann_parameters;
    boltzmann_parameters.reduced_field_td =
        std::max(1.0, plasma_state_.electric_field_v_per_m / plasma_state_.neutral_density_m3 * 1.0e21);
    boltzmann_solver_.setParameters(boltzmann_parameters);

    if (!diffusion_solver_.initialize(config_, plasma_state_.electron_density_m3))
    {
        return false;
    }

    poisson_solver_.setParameters(FieldSolver::NonlinearPoissonParameters{});
    if (!poisson_solver_.initializeGrid(config.domain_size, config.resolution))
    {
        return false;
    }

    poisson_solver_.setBoundaryPotential([config](const Geometry::Point3D& point) {
        return point.z() >= config.domain_size.z() ? config.initial_potential_v : 0.0;
    });
    poisson_solver_.setChargeDensityFunction([this](const Geometry::Point3D&, double) {
        return 1.602176634e-19 * (plasma_state_.ion_density_m3 - plasma_state_.electron_density_m3);
    });
    poisson_solver_.setPermittivityFunction([](const Geometry::Point3D&, double) { return 1.0; });
    poisson_solver_.solve();

    assessment_ = detector_.assess(plasma_state_);
    boundary_state_ = boundary_layer_.analyze(plasma_state_, assessment_, config.initial_potential_v);
    reaction_state_ = reaction_collision_library_.evaluate(plasma_state_, assessment_);
    const double characteristic_length_m =
        std::max(1.0e-9,
                 config.domain_size.z() /
                     std::max(1.0, config.resolution.z()));
    closure_state_ = advanced_closure_model_.evaluate(plasma_state_, assessment_, reaction_state_,
                                                      characteristic_length_m,
                                                      config.time_step_s);
    status_.average_density_m3 = plasma_state_.electron_density_m3;
    status_.average_potential_v = 0.5 * config.initial_potential_v;
    status_.debye_length_m = assessment_.debye_length_m;
    status_.sheath_thickness_m = boundary_state_.sheath_thickness_m;
    status_.presheath_thickness_m = boundary_state_.presheath_thickness_m;
    status_.sheath_potential_drop_v = boundary_state_.sheath_potential_drop_v;
    status_.presheath_potential_drop_v = boundary_state_.presheath_potential_drop_v;
    status_.boundary_wall_field_v_per_m = boundary_state_.wall_field_v_per_m;
    status_.ion_bohm_velocity_m_per_s = boundary_state_.bohm_velocity_m_per_s;
    status_.ion_mach_at_sheath_edge = boundary_state_.ion_mach_at_sheath_edge;
    status_.ionization_source_m3_per_s = reaction_state_.ionization_source_m3_per_s;
    status_.recombination_sink_m3_per_s = reaction_state_.recombination_sink_m3_per_s;
    status_.effective_collision_frequency_hz = reaction_state_.effective_collision_frequency_hz;
    status_.charge_exchange_frequency_hz = reaction_state_.charge_exchange_frequency_hz;
    status_.reaction_energy_loss_ev_per_s = reaction_state_.electron_energy_loss_ev_per_s;
    status_.reaction_momentum_transfer_ratio = reaction_state_.momentum_transfer_ratio;
    status_.reaction_active_processes = reaction_state_.active_processes;
    status_.advanced_closure_enabled =
        config.advanced_closure.enable_non_equilibrium_closure ||
        config.advanced_closure.enable_turbulence_closure;
    status_.non_equilibrium_ratio = closure_state_.non_equilibrium_ratio;
    status_.non_equilibrium_relaxation_rate_hz =
        closure_state_.non_equilibrium_relaxation_rate_hz;
    status_.electron_temperature_closure_delta_ev =
        closure_state_.electron_temperature_delta_ev;
    status_.turbulence_intensity = closure_state_.turbulence_intensity;
    status_.turbulence_eddy_diffusivity_m2_per_s =
        closure_state_.turbulence_eddy_diffusivity_m2_per_s;
    status_.turbulence_dissipation_rate_w_per_m3 =
        closure_state_.turbulence_dissipation_rate_w_per_m3;
    status_.closure_density_correction_m3 = closure_state_.density_correction_m3;
    status_.closure_active_terms = closure_state_.active_terms;
    status_.dense_plasma_detected = assessment_.is_dense;
    initialized_ = true;
    return true;
}

bool FluidAlgorithmAdapter::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }

    assessment_ = detector_.assess(plasma_state_);
    reaction_state_ = reaction_collision_library_.evaluate(plasma_state_, assessment_);

    const double net_reaction_source_m3_per_s =
        reaction_state_.ionization_source_m3_per_s - reaction_state_.recombination_sink_m3_per_s;
    const double density_baseline = std::max(1.0e6, plasma_state_.electron_density_m3);
    const double max_density_delta =
        std::max(0.0, config_.reaction_collision.max_relative_density_change_per_step) *
        density_baseline;
    const double reaction_density_delta =
        std::clamp(net_reaction_source_m3_per_s * dt, -max_density_delta, max_density_delta);

    const double characteristic_length_m =
        std::max(1.0e-9,
                 config_.domain_size.z() /
                     std::max(1.0, config_.resolution.z()));
    closure_state_ = advanced_closure_model_.evaluate(plasma_state_, assessment_, reaction_state_,
                                                      characteristic_length_m,
                                                      dt);
    const double closure_density_delta = std::clamp(
        closure_state_.density_correction_m3,
        -0.5 * max_density_delta,
        0.5 * max_density_delta);

    plasma_state_.electron_temperature_ev =
        std::max(0.05,
                 plasma_state_.electron_temperature_ev -
                     reaction_state_.electron_energy_loss_ev_per_s * dt);
    plasma_state_.electron_temperature_ev =
        std::max(0.05,
                 plasma_state_.electron_temperature_ev +
                     closure_state_.electron_temperature_delta_ev);

    const auto& boltzmann_state = boltzmann_solver_.solve(plasma_state_.electron_temperature_ev);
    plasma_state_.electric_field_v_per_m = boltzmann_state.mean_energy_ev * 500.0;

    const auto& electric_field = poisson_solver_.getElectricField();
    if (!diffusion_solver_.advance(electric_field.empty()
                                       ? std::vector<Geometry::Vector3D>(
                                             diffusion_solver_.getState().density.size())
                                       : electric_field,
                                   dt))
    {
        return false;
    }

    const auto& density = diffusion_solver_.getDensityDistribution();
    const double average_density =
        density.empty() ? 0.0
                        : std::accumulate(density.begin(), density.end(), 0.0) /
                              static_cast<double>(density.size());
    const double reaction_adjusted_density =
        std::max(1.0e6, average_density + reaction_density_delta + closure_density_delta);
    plasma_state_.electron_density_m3 = reaction_adjusted_density;
    plasma_state_.ion_density_m3 =
        std::max(1.0e6, plasma_state_.ion_density_m3 + reaction_density_delta +
                           0.5 * closure_density_delta);

    poisson_solver_.setChargeDensityFunction([this](const Geometry::Point3D& point, double phi) {
        const double sheath_boost =
            1.0 + 0.05 * std::abs(phi) / std::max(1.0e-3, config_.initial_potential_v + 1.0);
        const double ne = plasma_state_.electron_density_m3 / sheath_boost;
        (void)point;
        return 1.602176634e-19 * (plasma_state_.ion_density_m3 - ne);
    });
    poisson_solver_.solve();

    const auto& potential = poisson_solver_.getPotentialField();
    status_.average_potential_v =
        potential.empty()
            ? 0.0
            : std::accumulate(potential.begin(), potential.end(), 0.0) /
                  static_cast<double>(potential.size());
    status_.average_density_m3 = reaction_adjusted_density;
    status_.time_s += dt;
    status_.steps_completed += 1;
    status_.dense_plasma_detected = assessment_.is_dense;
    status_.debye_length_m = assessment_.debye_length_m;
    boundary_state_ = boundary_layer_.analyze(plasma_state_, assessment_, status_.average_potential_v);
    status_.sheath_thickness_m = boundary_state_.sheath_thickness_m;
    status_.presheath_thickness_m = boundary_state_.presheath_thickness_m;
    status_.sheath_potential_drop_v = boundary_state_.sheath_potential_drop_v;
    status_.presheath_potential_drop_v = boundary_state_.presheath_potential_drop_v;
    status_.boundary_wall_field_v_per_m = boundary_state_.wall_field_v_per_m;
    status_.ion_bohm_velocity_m_per_s = boundary_state_.bohm_velocity_m_per_s;
    status_.ion_mach_at_sheath_edge = boundary_state_.ion_mach_at_sheath_edge;
    status_.ionization_source_m3_per_s = reaction_state_.ionization_source_m3_per_s;
    status_.recombination_sink_m3_per_s = reaction_state_.recombination_sink_m3_per_s;
    status_.effective_collision_frequency_hz = reaction_state_.effective_collision_frequency_hz;
    status_.charge_exchange_frequency_hz = reaction_state_.charge_exchange_frequency_hz;
    status_.reaction_energy_loss_ev_per_s = reaction_state_.electron_energy_loss_ev_per_s;
    status_.reaction_momentum_transfer_ratio = reaction_state_.momentum_transfer_ratio;
    status_.reaction_active_processes = reaction_state_.active_processes;
    status_.advanced_closure_enabled =
        config_.advanced_closure.enable_non_equilibrium_closure ||
        config_.advanced_closure.enable_turbulence_closure;
    status_.non_equilibrium_ratio = closure_state_.non_equilibrium_ratio;
    status_.non_equilibrium_relaxation_rate_hz =
        closure_state_.non_equilibrium_relaxation_rate_hz;
    status_.electron_temperature_closure_delta_ev =
        closure_state_.electron_temperature_delta_ev;
    status_.turbulence_intensity = closure_state_.turbulence_intensity;
    status_.turbulence_eddy_diffusivity_m2_per_s =
        closure_state_.turbulence_eddy_diffusivity_m2_per_s;
    status_.turbulence_dissipation_rate_w_per_m3 =
        closure_state_.turbulence_dissipation_rate_w_per_m3;
    status_.closure_density_correction_m3 = closure_state_.density_correction_m3;
    status_.closure_active_terms = closure_state_.active_terms;
    return true;
}

void FluidAlgorithmAdapter::reset()
{
    *this = FluidAlgorithmAdapter{};
}

Output::ColumnarDataSet FluidAlgorithmAdapter::buildProfileDataSet() const
{
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "z_m";

    const auto& diffusion_state = diffusion_solver_.getState();
    const auto& density = diffusion_state.density;
    const auto& potential = poisson_solver_.getPotentialField();
    const auto& electric_field = poisson_solver_.getElectricField();

    const std::size_t nx = static_cast<std::size_t>(diffusion_state.resolution.x());
    const std::size_t ny = static_cast<std::size_t>(diffusion_state.resolution.y());
    const std::size_t nz = static_cast<std::size_t>(diffusion_state.resolution.z());
    data_set.scalar_series["potential_v"] = std::vector<double>(nz, 0.0);
    data_set.scalar_series["electric_field_z_v_per_m"] = std::vector<double>(nz, 0.0);
    data_set.scalar_series["electron_density_m3"] = std::vector<double>(nz, 0.0);
    data_set.scalar_series["ion_density_m3"] = std::vector<double>(nz, plasma_state_.ion_density_m3);
    data_set.scalar_series["electron_temperature_ev"] =
        std::vector<double>(nz, plasma_state_.electron_temperature_ev);

    for (std::size_t k = 0; k < nz; ++k)
    {
        data_set.axis_values.push_back(static_cast<double>(k) * diffusion_state.spacing.z());
        double potential_sum = 0.0;
        double field_sum = 0.0;
        double density_sum = 0.0;

        for (std::size_t j = 0; j < ny; ++j)
        {
            for (std::size_t i = 0; i < nx; ++i)
            {
                const std::size_t index = i + nx * (j + ny * k);
                if (index < potential.size())
                {
                    potential_sum += potential[index];
                    field_sum += electric_field[index].z();
                }
                if (index < density.size())
                {
                    density_sum += density[index];
                }
            }
        }

        const double plane_scale = 1.0 / static_cast<double>(nx * ny);
        data_set.scalar_series["potential_v"][k] = potential_sum * plane_scale;
        data_set.scalar_series["electric_field_z_v_per_m"][k] = field_sum * plane_scale;
        data_set.scalar_series["electron_density_m3"][k] = density_sum * plane_scale;
    }

    data_set.metadata["module"] = "Plasma Analysis";
    data_set.metadata["organization_family"] =
        config_.enable_spis_style_organization ? "spis_numeric_v1" : "native";
    data_set.metadata["environment_model"] = environmentModelName(config_.environment_model);
    data_set.metadata["distribution_model"] = distributionModelName(config_.distribution_model);
    data_set.metadata["reaction_registry"] = reactionRegistryName(config_.reaction_registry);
    data_set.metadata["diagnostic_set"] = diagnosticSetName(config_.diagnostic_set);
    data_set.metadata["reaction_contract_id"] = "plasma-reaction-balance-v1";
    data_set.metadata["diagnostic_contract_id"] = "plasma-physics-diagnostics-v1";
    data_set.metadata["diagnostic_contract_family"] = "spis_physics_diagnostics";
    data_set.metadata["diagnostic_contract_supports_boundary_layer"] = "true";
    data_set.metadata["diagnostic_contract_supports_non_equilibrium"] = "true";
    data_set.metadata["domain_size_xyz_m"] =
        std::to_string(config_.domain_size.x()) + "," + std::to_string(config_.domain_size.y()) +
        "," + std::to_string(config_.domain_size.z());
    data_set.metadata["resolution_xyz"] =
        std::to_string(static_cast<long long>(config_.resolution.x())) + "," +
        std::to_string(static_cast<long long>(config_.resolution.y())) + "," +
        std::to_string(static_cast<long long>(config_.resolution.z()));
    data_set.metadata["sheath_thickness_m"] = std::to_string(status_.sheath_thickness_m);
    data_set.metadata["presheath_thickness_m"] = std::to_string(status_.presheath_thickness_m);
    data_set.metadata["sheath_potential_drop_v"] = std::to_string(status_.sheath_potential_drop_v);
    data_set.metadata["presheath_potential_drop_v"] = std::to_string(status_.presheath_potential_drop_v);
    data_set.metadata["boundary_wall_field_v_per_m"] = std::to_string(status_.boundary_wall_field_v_per_m);
    data_set.metadata["ion_bohm_velocity_m_per_s"] = std::to_string(status_.ion_bohm_velocity_m_per_s);
    data_set.metadata["ion_mach_at_sheath_edge"] = std::to_string(status_.ion_mach_at_sheath_edge);
    data_set.metadata["ionization_source_m3_per_s"] = std::to_string(status_.ionization_source_m3_per_s);
    data_set.metadata["recombination_sink_m3_per_s"] = std::to_string(status_.recombination_sink_m3_per_s);
    data_set.metadata["effective_collision_frequency_hz"] =
        std::to_string(status_.effective_collision_frequency_hz);
    data_set.metadata["charge_exchange_frequency_hz"] =
        std::to_string(status_.charge_exchange_frequency_hz);
    data_set.metadata["reaction_energy_loss_ev_per_s"] =
        std::to_string(status_.reaction_energy_loss_ev_per_s);
    data_set.metadata["reaction_momentum_transfer_ratio"] =
        std::to_string(status_.reaction_momentum_transfer_ratio);
    data_set.metadata["reaction_active_processes"] =
        std::to_string(status_.reaction_active_processes);
    data_set.metadata["advanced_closure_enabled"] =
        status_.advanced_closure_enabled ? "1" : "0";
    data_set.metadata["non_equilibrium_ratio"] =
        std::to_string(status_.non_equilibrium_ratio);
    data_set.metadata["non_equilibrium_relaxation_rate_hz"] =
        std::to_string(status_.non_equilibrium_relaxation_rate_hz);
    data_set.metadata["electron_temperature_closure_delta_ev"] =
        std::to_string(status_.electron_temperature_closure_delta_ev);
    data_set.metadata["turbulence_intensity"] =
        std::to_string(status_.turbulence_intensity);
    data_set.metadata["turbulence_eddy_diffusivity_m2_per_s"] =
        std::to_string(status_.turbulence_eddy_diffusivity_m2_per_s);
    data_set.metadata["turbulence_dissipation_rate_w_per_m3"] =
        std::to_string(status_.turbulence_dissipation_rate_w_per_m3);
    data_set.metadata["closure_density_correction_m3"] =
        std::to_string(status_.closure_density_correction_m3);
    data_set.metadata["closure_active_terms"] =
        std::to_string(status_.closure_active_terms);
    return data_set;
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
