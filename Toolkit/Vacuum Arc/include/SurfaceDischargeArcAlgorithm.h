#pragma once

#include "AnodeHeatingModel.h"
#include "ArcDischargeDetector.h"
#include "ArcFieldEmissionBoundaryCondition.h"
#include "ArcEmissionStrategy.h"
#include "ArcPICEmissionModel.h"
#include "ArcPICIntegrator.h"
#include "CathodeSpotModel.h"
#include "CylindricalPICSolver.h"
#include "PlasmaChannelModel.h"
#include "SurfaceArcCoupling.h"
#include "TownsendAvalancheModel.h"

#include "../../Tools/Coupling/include/BenchmarkContracts.h"
#include "../../Tools/FieldSolver/include/PoissonSolver.h"
#include "../../Tools/Interactions/Collisions/include/CollisionPicAdapter.h"
#include "../../Tools/Mesh/include/MeshParsing.h"
#include "../../Tools/Output/include/ResultExporter.h"
#include "../../Tools/PICcore/include/PICCycle.h"
#include "../../Tools/Particle/include/ParticleManager.h"

#include <filesystem>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

enum class ArcPicAlignmentMode
{
  ArcPicAligned,
  LegacyBaseline,
};

enum class ArcPipelineStage : std::uint8_t
{
  Initialize = 0,
  FieldSolve = 1,
  Advance = 2,
  Collisions = 3,
  Circuit = 4,
  Output = 5,
};

enum class SurfaceCircuitLoadModel : std::uint8_t
{
    LegacyLeakageCapacitor = 0,
    ResistiveShunt = 1,
    RlBranch = 2,
};

struct DischargeConfiguration
{
    Coupling::Contracts::SolverConfig solver_config{};
    unsigned int seed = 20260408u;
    std::string sampling_policy = "deterministic";
    double gap_distance_m = 5.0e-4;
    double applied_field_v_per_m = 2.5e7;
    double surface_potential_v = 0.0;
    double surface_charge_density_c_per_m2 = 0.0;
    double cathode_temperature_k = 450.0;
    double anode_temperature_k = 320.0;
    double channel_radius_m = 2.0e-5;
    double surface_capacitance_f = 2.0e-10;
    double surface_leakage_conductance_s = 0.0;
    double surface_reference_potential_v = 0.0;
    double max_surface_potential_step_v = 25.0;
    SurfaceCircuitLoadModel surface_load_model = SurfaceCircuitLoadModel::LegacyLeakageCapacitor;
    double surface_load_resistance_ohm = 0.0;
    double surface_load_inductance_h = 0.0;
    double surface_load_current_limit_a = 1.0e6;
    ArcPicAlignmentMode alignment_mode = ArcPicAlignmentMode::ArcPicAligned;
    PlasmaChannelParameters channel_parameters{};
    ArcIntegratorParameters integrator_parameters{};
    ArcEmissionStrategyType emission_strategy = ArcEmissionStrategyType::LinearThreshold;
    ArcEmissionStrategyParameters emission_parameters{};
    bool enable_pic_mcc_collisions = true;
    std::size_t collision_reconfigure_interval_steps = 16;
    double collision_neutral_density_floor_m3 = 1.0e20;
    double collision_ion_density_floor_m3 = 1.0e16;
    double collision_electron_density_floor_m3 = 1.0e16;
    bool enable_collision_reaction_breakdown = false;
    std::string collision_cross_section_set_id = "background_mcc_v1";
    double collision_anomaly_event_rate_limit_per_s = 0.0;
    bool collision_fallback_to_background_mcc_on_error = true;
    double collision_emission_feedback_gain = 0.0;
    double collision_channel_feedback_gain = 0.0;
    double collision_ionization_emission_feedback_gain = 0.0;
    double collision_excitation_emission_feedback_gain = 0.0;
    double collision_charge_exchange_emission_feedback_gain = 0.0;
    double collision_ionization_channel_feedback_gain = 0.0;
    double collision_excitation_channel_feedback_gain = 0.0;
    double collision_charge_exchange_channel_feedback_gain = 0.0;
    bool enable_neutral_outgassing_feedback = false;
    double neutral_outgassing_gain_m3_per_a = 0.0;
    double neutral_outgassing_relaxation_per_s = 0.0;
    double neutral_outgassing_max_density_boost_m3 = 0.0;
    bool enable_neutral_outgassing_reinjection = false;
    double neutral_outgassing_reinjection_gain = 1.0;
    double neutral_outgassing_reinjection_macro_weight = 1.0e6;
    std::size_t neutral_outgassing_reinjection_max_particles_per_step = 16;
    double neutral_outgassing_reinjection_speed_m_per_s = 1.0e3;
    double neutral_outgassing_reinjection_mass_amu = 28.0;
    bool enable_residual_monitoring = true;
    double residual_field_tolerance_v_per_m = 1.0e6;
    double residual_particle_current_tolerance_a = 1.0e-6;
    double residual_circuit_charge_tolerance_c = 1.0e-12;
    std::size_t residual_alarm_consecutive_steps = 2;
    std::size_t residual_alarm_clear_consecutive_steps = 2;
    bool enable_residual_ema_filter = true;
    double residual_ema_alpha = 0.25;
    bool enable_secondary_electron_emission = false;
    double secondary_electron_yield_per_ion = 0.0;
    double secondary_electron_speed_m_per_s = 8.0e5;
    double min_internal_substep_s = 1.0e-12;
    double max_internal_substep_s = 1.0;
    bool enable_adaptive_internal_timestep = false;
    double adaptive_substep_shrink_factor = 0.5;
    double adaptive_substep_recovery_factor = 1.25;
    double adaptive_substep_min_scale = 0.125;
    bool enable_stability_rollback = true;
    std::size_t max_stability_rollbacks = 2;
    double anomaly_current_density_limit_a_per_m2 = 1.0e9;
    double anomaly_surface_potential_limit_v = 1.0e6;
};

struct DischargeStatus
{
    bool discharge_active = false;
    double current_time_s = 0.0;
    double total_discharge_current_a = 0.0;
    double peak_current_density_a_per_m2 = 0.0;
    std::size_t active_arc_channels = 0;
    double cathode_temperature_k = 0.0;
    double anode_temperature_k = 0.0;
    double channel_conductivity_s_per_m = 0.0;
    double pic_absorbed_electron_current_a = 0.0;
    double pic_absorbed_ion_current_a = 0.0;
    double pic_emitted_electron_current_a = 0.0;
    double pic_boundary_net_current_a = 0.0;
    double pic_surface_charge_delta_c = 0.0;
    double surface_potential_v = 0.0;
    double surface_charge_density_c_per_m2 = 0.0;
    double surface_circuit_drive_current_a = 0.0;
    double surface_circuit_leak_current_a = 0.0;
    double surface_load_branch_current_a = 0.0;
    double surface_relative_potential_v = 0.0;
    double surface_circuit_charge_c = 0.0;
    double surface_circuit_capacitor_energy_j = 0.0;
    double surface_circuit_drive_power_w = 0.0;
    double surface_circuit_leak_power_w = 0.0;
    double surface_load_branch_power_w = 0.0;
    double surface_load_resistive_power_w = 0.0;
    double surface_load_inductive_energy_j = 0.0;
    double surface_circuit_power_balance_error_w = 0.0;
    double surface_circuit_drive_energy_j = 0.0;
    double surface_circuit_leak_dissipated_energy_j = 0.0;
    double surface_load_dissipated_energy_j = 0.0;
    double surface_circuit_energy_balance_error_j = 0.0;
    double charge_conservation_error_c = 0.0;
    double collision_events_step = 0.0;
    double collision_particles_created_step = 0.0;
    double collision_ionization_events_step = 0.0;
    double collision_excitation_events_step = 0.0;
    double collision_charge_exchange_events_step = 0.0;
    double collision_ionization_fraction_step = 0.0;
    double collision_excitation_fraction_step = 0.0;
    double collision_charge_exchange_fraction_step = 0.0;
    double collision_reaction_weighted_emission_feedback_step = 1.0;
    double collision_reaction_weighted_channel_feedback_step = 1.0;
    double collision_event_rate_per_s = 0.0;
    bool collision_stage_fallback_triggered = false;
    double collision_effective_neutral_density_m3 = 0.0;
    double collision_effective_ion_density_m3 = 0.0;
    double collision_effective_electron_density_m3 = 0.0;
    double neutral_outgassing_density_boost_m3 = 0.0;
    double neutral_outgassing_reinjected_particles_step = 0.0;
    double neutral_outgassing_reinjection_rate_per_s = 0.0;
    double field_residual_v_per_m = 0.0;
    double particle_residual_a = 0.0;
    double circuit_residual_c = 0.0;
    double field_residual_filtered_v_per_m = 0.0;
    double particle_residual_filtered_a = 0.0;
    double circuit_residual_filtered_c = 0.0;
    bool residual_alarm_active = false;
    std::size_t residual_alarm_counter = 0;
    std::size_t residual_alarm_clear_counter = 0;
    std::size_t stability_substeps_used = 0;
    std::size_t stability_rollbacks = 0;
    double stability_effective_substep_s = 0.0;
    double stability_adaptive_scale = 1.0;
    std::size_t stability_adaptive_reductions = 0;
    bool stability_anomaly_isolation_triggered = false;
    std::size_t stability_anomaly_isolation_events = 0;
    ArcPicAlignmentMode alignment_mode = ArcPicAlignmentMode::ArcPicAligned;
};

class SurfaceDischargeArcAlgorithm
{
  public:
    bool initialize(const DischargeConfiguration& config);
    bool advance(double dt);
    const DischargeStatus& getStatus() const { return status_; }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

  private:
    double computeEffectiveField() const;
    bool updateDischargeActivation(double effective_field, double avalanche_gain);
    double computeEmissionCurrentDensity(double effective_field) const;
    void updateChannelState(double effective_field, double emitted_current_density, double dt);
    void updateThermalState(double emitted_current_density, double dt);
    bool initializePiccoreBoundaryCoupling();
    void applyPiccoreBoundaryPotentials();
    void seedBoundaryProbeParticles(double emitted_current_density, double dt);
    void resetPicBoundaryFeedback();
    void advancePiccoreBoundaryCoupling(double emitted_current_density, double dt);
    double estimateCollisionVolumeM3() const;
    void runCollisionStage(double dt);
    void injectOutgassingNeutralParticles(double volume_m3, double source_density_boost_m3,
                        double dt);
    void updateSurfaceCircuitState(double dt);
    void updateResidualMonitoring(double effective_field, double dt);
    bool hasAnomalousState() const;
    void appendHistory();
    static const char* alignmentModeName(ArcPicAlignmentMode mode);
    static const char* surfaceLoadModelName(SurfaceCircuitLoadModel mode);
    static const char* pipelineStageOrder();
    static const char* pipelineContractId();

    DischargeConfiguration config_;
    DischargeStatus status_;
    TownsendAvalancheModel avalanche_model_;
    CathodeSpotModel cathode_spot_model_;
    AnodeHeatingModel anode_heating_model_;
    PlasmaChannelModel plasma_channel_model_;
    ArcDischargeDetector detector_;
    ArcPICEmissionModel emission_model_;
    ArcPICIntegrator integrator_;
    CylindricalPICSolver cylindrical_solver_;
    SurfaceArcCoupling coupling_;
    Output::ResultExporter exporter_;
    ArcChannelState channel_state_;
    CylindricalProfile radial_profile_;
    std::vector<double> history_time_;
    std::vector<double> history_current_;
    std::vector<double> history_current_density_;
    std::vector<double> history_cathode_temperature_;
    std::vector<double> history_anode_temperature_;
    std::vector<double> history_conductivity_;
    std::vector<double> history_pic_boundary_net_current_;
    std::vector<double> history_pic_surface_charge_delta_;
    std::vector<double> history_surface_potential_;
    std::vector<double> history_surface_charge_density_;
    std::vector<double> history_surface_circuit_drive_current_;
    std::vector<double> history_surface_circuit_leak_current_;
    std::vector<double> history_surface_load_branch_current_;
    std::vector<double> history_surface_relative_potential_;
    std::vector<double> history_surface_circuit_capacitor_energy_;
    std::vector<double> history_surface_circuit_drive_power_;
    std::vector<double> history_surface_circuit_leak_power_;
    std::vector<double> history_surface_load_branch_power_;
    std::vector<double> history_surface_load_resistive_power_;
    std::vector<double> history_surface_load_inductive_energy_;
    std::vector<double> history_surface_circuit_power_balance_error_;
    std::vector<double> history_surface_circuit_drive_energy_;
    std::vector<double> history_surface_circuit_leak_dissipated_energy_;
    std::vector<double> history_surface_load_dissipated_energy_;
    std::vector<double> history_surface_circuit_energy_balance_error_;
    std::vector<double> history_charge_conservation_error_;
    std::vector<double> history_collision_events_;
    std::vector<double> history_collision_particles_created_;
    std::vector<double> history_collision_ionization_events_;
    std::vector<double> history_collision_excitation_events_;
    std::vector<double> history_collision_charge_exchange_events_;
    std::vector<double> history_collision_ionization_fraction_;
    std::vector<double> history_collision_excitation_fraction_;
    std::vector<double> history_collision_charge_exchange_fraction_;
    std::vector<double> history_collision_emission_feedback_multiplier_;
    std::vector<double> history_collision_channel_feedback_multiplier_;
    std::vector<double> history_collision_event_rate_;
    std::vector<double> history_collision_fallback_triggered_;
    std::vector<double> history_collision_neutral_density_;
    std::vector<double> history_collision_ion_density_;
    std::vector<double> history_collision_electron_density_;
    std::vector<double> history_neutral_outgassing_boost_;
    std::vector<double> history_neutral_outgassing_reinjected_particles_;
    std::vector<double> history_neutral_outgassing_reinjection_rate_;
    std::vector<double> history_field_residual_;
    std::vector<double> history_particle_residual_;
    std::vector<double> history_circuit_residual_;
    std::vector<double> history_field_residual_filtered_;
    std::vector<double> history_particle_residual_filtered_;
    std::vector<double> history_circuit_residual_filtered_;
    std::vector<double> history_residual_alarm_active_;
    std::vector<double> history_residual_alarm_counter_;
    std::vector<double> history_residual_alarm_clear_counter_;
    std::vector<double> history_stability_substeps_;
    std::vector<double> history_stability_rollbacks_;
    std::vector<double> history_stability_effective_substep_;
    std::vector<double> history_stability_adaptive_scale_;
    std::vector<double> history_stability_adaptive_reductions_;
    std::vector<double> history_stability_anomaly_isolation_triggered_;
    std::vector<double> history_stability_anomaly_isolation_events_;

    Mesh::VolMeshPtr pic_boundary_mesh_;
    std::shared_ptr<Particle::ParticleManager> pic_boundary_particle_manager_;
    std::shared_ptr<FieldSolver::PoissonSolver> pic_boundary_poisson_solver_;
    std::unique_ptr<PICcore::PICCycle> pic_boundary_cycle_;
    std::shared_ptr<ArcFieldEmissionBoundaryCondition> pic_cathode_boundary_;
    bool piccore_boundary_coupling_initialized_ = false;
    double pic_surface_area_m2_ = 1.0;
    std::size_t pic_top_element_id_ = 0;
    double pic_last_absorbed_electron_charge_c_ = 0.0;
    double pic_last_absorbed_ion_charge_c_ = 0.0;
    double pic_last_emitted_electron_charge_c_ = 0.0;
    double surface_circuit_charge_c_ = 0.0;
    double surface_load_branch_current_a_ = 0.0;
    double accumulated_drive_charge_c_ = 0.0;
    double accumulated_leak_charge_c_ = 0.0;
    double accumulated_load_charge_c_ = 0.0;
    double accumulated_drive_energy_j_ = 0.0;
    double accumulated_leak_dissipated_energy_j_ = 0.0;
    double accumulated_load_dissipated_energy_j_ = 0.0;
    double accumulated_power_balance_error_j_ = 0.0;
    std::size_t residual_alarm_consecutive_counter_ = 0;
    std::size_t residual_alarm_clear_consecutive_counter_ = 0;
    bool residual_filter_initialized_ = false;
    double filtered_field_residual_state_ = 0.0;
    double filtered_particle_residual_state_ = 0.0;
    double filtered_circuit_residual_state_ = 0.0;
    double pending_outgassing_reinjection_particles_ = 0.0;
    double dynamic_neutral_density_boost_m3_ = 0.0;
    double collision_emission_feedback_multiplier_ = 1.0;
    double collision_channel_feedback_multiplier_ = 1.0;
    std::mt19937 neutral_outgassing_reinjection_rng_{20260403u};
    Collision::MonteCarloCollisionHandler collision_handler_{12345};
    bool collision_handler_initialized_ = false;
    std::size_t collision_stage_step_counter_ = 0;
    std::size_t collision_last_total_events_ = 0;
    std::size_t collision_last_total_particles_created_ = 0;

    bool initialized_ = false;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
