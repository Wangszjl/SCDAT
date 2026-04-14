#include "SurfaceDischargeArcAlgorithm.h"

#include "../../Tools/Boundary/include/BoundaryConditions.h"
#include "../../Tools/Basic/include/Constants.h"
#include "../../Tools/Solver/include/SurfaceSolverFacade.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

namespace
{
using BoundaryFace = PICcore::PICCycle::BoundaryFace;

constexpr const char* kArcPipelineContractId = "arcpic-6-stage-v1";
constexpr const char* kArcPipelineStageOrder =
    "initialize>field_solve>advance>collisions>circuit>output";

const char* integratorModeName(ArcIntegratorMode mode)
{
    switch (mode)
    {
    case ArcIntegratorMode::ExplicitGrowth:
        return "explicit_growth";
    case ArcIntegratorMode::BoundedRelaxation:
        return "bounded_relaxation";
    default:
        return "unknown";
    }
}

bool isVacuumArcDebugEnabled()
{
    return std::getenv("VACUUM_ARC_DEBUG") != nullptr;
}

std::string normalizeToken(std::string text)
{
    std::string token;
    token.reserve(text.size());
    for (const unsigned char ch : text)
    {
        if (std::isalnum(ch))
        {
            token.push_back(static_cast<char>(std::tolower(ch)));
        }
    }
    return token;
}

std::size_t boundaryIndex(BoundaryFace face)
{
    return static_cast<std::size_t>(face);
}

std::size_t findTopElementId(const Mesh::VolMesh& mesh)
{
    std::size_t best_id = 0;
    double best_z = -1.0e300;
    for (std::size_t i = 0; i < mesh.getElementCount(); ++i)
    {
        const auto& element = mesh.getElement(i);
        if (!element)
        {
            continue;
        }

        const double z = element->getCentroid().z();
        if (z > best_z)
        {
            best_z = z;
            best_id = i;
        }
    }

    return best_id;
}

double estimateMeshVolume(const Mesh::VolMesh& mesh)
{
    const auto min_corner = mesh.getBoundingBoxMin();
    const auto max_corner = mesh.getBoundingBoxMax();
    const double dx = std::max(0.0, max_corner.x() - min_corner.x());
    const double dy = std::max(0.0, max_corner.y() - min_corner.y());
    const double dz = std::max(0.0, max_corner.z() - min_corner.z());
    const double bounding_volume = dx * dy * dz;
    if (bounding_volume > 0.0)
    {
        return bounding_volume;
    }

    double element_volume_sum = 0.0;
    for (std::size_t i = 0; i < mesh.getElementCount(); ++i)
    {
        const auto& element = mesh.getElement(i);
        if (!element)
        {
            continue;
        }
        element_volume_sum += std::abs(element->getVolume());
    }

    return element_volume_sum;
}

enum class CollisionStageAnomalyReason
{
    None = 0,
    NonFiniteDiagnostics,
    EventRateExceeded,
    BreakdownInconsistency,
};

bool hasFiniteCollisionDiagnostics(const Collision::CollisionStepResult& step_result)
{
    return std::isfinite(step_result.effective_neutral_density_m3) &&
           std::isfinite(step_result.effective_ion_density_m3) &&
           std::isfinite(step_result.effective_electron_density_m3);
}

double sanitizeNonNegativeFinite(double value)
{
    if (!std::isfinite(value))
    {
        return 0.0;
    }
    return std::max(0.0, value);
}

struct CollisionReactionFractions
{
    double ionization_fraction = 0.0;
    double excitation_fraction = 0.0;
    double charge_exchange_fraction = 0.0;
};

CollisionReactionFractions
computeCollisionReactionFractions(const Collision::CollisionStepResult& step_result,
                                  bool enable_collision_reaction_breakdown)
{
    CollisionReactionFractions fractions{};

    if (enable_collision_reaction_breakdown)
    {
        const double ionization_events =
            sanitizeNonNegativeFinite(static_cast<double>(step_result.ionization_events));
        const double excitation_events =
            sanitizeNonNegativeFinite(static_cast<double>(step_result.excitation_events));
        const double charge_exchange_events =
            sanitizeNonNegativeFinite(static_cast<double>(step_result.charge_exchange_events));
        const double total_reaction_events =
            ionization_events + excitation_events + charge_exchange_events;
        if (total_reaction_events > 0.0)
        {
            const double inv_total = 1.0 / total_reaction_events;
            fractions.ionization_fraction = ionization_events * inv_total;
            fractions.excitation_fraction = excitation_events * inv_total;
            fractions.charge_exchange_fraction = charge_exchange_events * inv_total;
            return fractions;
        }
    }

    // Fallback proxy when explicit reaction counters are unavailable: use effective species
    // densities to keep stage-C structured feedback active and finite.
    const double effective_electron_density =
        sanitizeNonNegativeFinite(step_result.effective_electron_density_m3);
    const double effective_neutral_density =
        sanitizeNonNegativeFinite(step_result.effective_neutral_density_m3);
    const double effective_ion_density =
        sanitizeNonNegativeFinite(step_result.effective_ion_density_m3);
    const double density_sum =
        effective_electron_density + effective_neutral_density + effective_ion_density;
    if (density_sum <= 0.0)
    {
        return fractions;
    }

    const double inv_sum = 1.0 / density_sum;
    fractions.ionization_fraction = effective_electron_density * inv_sum;
    fractions.excitation_fraction = effective_neutral_density * inv_sum;
    fractions.charge_exchange_fraction = effective_ion_density * inv_sum;
    return fractions;
}

double computeReactionGainDelta(const CollisionReactionFractions& fractions,
                                double ionization_gain,
                                double excitation_gain,
                                double charge_exchange_gain)
{
    return fractions.ionization_fraction * ionization_gain +
           fractions.excitation_fraction * excitation_gain +
           fractions.charge_exchange_fraction * charge_exchange_gain;
}

CollisionStageAnomalyReason detectCollisionStageAnomaly(
    const Collision::CollisionStepResult& step_result, double dt,
    double collision_anomaly_event_rate_limit_per_s, bool enable_collision_reaction_breakdown)
{
    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        return CollisionStageAnomalyReason::NonFiniteDiagnostics;
    }

    if (!hasFiniteCollisionDiagnostics(step_result))
    {
        return CollisionStageAnomalyReason::NonFiniteDiagnostics;
    }

    const double event_rate_per_s =
        static_cast<double>(step_result.collision_events) / std::max(dt, 1.0e-18);
    if (!std::isfinite(event_rate_per_s))
    {
        return CollisionStageAnomalyReason::NonFiniteDiagnostics;
    }

    if (collision_anomaly_event_rate_limit_per_s > 0.0 &&
        event_rate_per_s > collision_anomaly_event_rate_limit_per_s)
    {
        return CollisionStageAnomalyReason::EventRateExceeded;
    }

    if (enable_collision_reaction_breakdown)
    {
        const std::size_t breakdown_sum =
            step_result.ionization_events + step_result.excitation_events +
            step_result.charge_exchange_events;
        if (breakdown_sum > step_result.collision_events)
        {
            return CollisionStageAnomalyReason::BreakdownInconsistency;
        }
    }

    return CollisionStageAnomalyReason::None;
}
} // namespace

const char* SurfaceDischargeArcAlgorithm::alignmentModeName(ArcPicAlignmentMode mode)
{
    switch (mode)
    {
    case ArcPicAlignmentMode::ArcPicAligned:
        return "arcpic_aligned";
    case ArcPicAlignmentMode::LegacyBaseline:
        return "legacy_baseline";
    default:
        return "unknown";
    }
}

const char* SurfaceDischargeArcAlgorithm::surfaceLoadModelName(SurfaceCircuitLoadModel mode)
{
    switch (mode)
    {
    case SurfaceCircuitLoadModel::LegacyLeakageCapacitor:
        return "legacy_leakage_capacitor";
    case SurfaceCircuitLoadModel::ResistiveShunt:
        return "resistive_shunt";
    case SurfaceCircuitLoadModel::RlBranch:
        return "rl_branch";
    default:
        return "unknown";
    }
}

const char* SurfaceDischargeArcAlgorithm::pipelineStageOrder()
{
    return kArcPipelineStageOrder;
}

const char* SurfaceDischargeArcAlgorithm::pipelineContractId()
{
    return kArcPipelineContractId;
}

double SurfaceDischargeArcAlgorithm::computeEffectiveField() const
{
    return coupling_.computeEffectiveField(config_.applied_field_v_per_m, config_.surface_potential_v,
                                           config_.surface_charge_density_c_per_m2);
}

bool SurfaceDischargeArcAlgorithm::updateDischargeActivation(double effective_field,
                                                             double avalanche_gain)
{
    if (!status_.discharge_active)
    {
        status_.discharge_active = detector_.shouldTrigger(effective_field, avalanche_gain);
        status_.active_arc_channels = status_.discharge_active ? 1u : 0u;
    }

    return status_.discharge_active;
}

double SurfaceDischargeArcAlgorithm::computeEmissionCurrentDensity(double effective_field) const
{
    const double cathode_emission =
        cathode_spot_model_.emitCurrentDensity(effective_field, status_.cathode_temperature_k);
    const double strategy_emission =
        emission_model_.computeEmissionCurrentDensity(effective_field, status_.cathode_temperature_k);

    return std::max(0.0, (strategy_emission + cathode_emission) *
                             std::max(0.0, collision_emission_feedback_multiplier_));
}

void SurfaceDischargeArcAlgorithm::updateChannelState(double effective_field,
                                                      double emitted_current_density, double dt)
{
    const double effective_emitted_current_density =
        emitted_current_density * std::max(0.0, collision_channel_feedback_multiplier_);

    channel_state_ =
        plasma_channel_model_.update(effective_field, effective_emitted_current_density,
                                     config_.gap_distance_m);
    channel_state_ = integrator_.advance(channel_state_, dt);
    radial_profile_ =
        cylindrical_solver_.solve(config_.channel_radius_m, channel_state_.current_density_a_per_m2, 24);

    status_.peak_current_density_a_per_m2 =
        std::max(status_.peak_current_density_a_per_m2, channel_state_.current_density_a_per_m2);
    status_.total_discharge_current_a =
        channel_state_.current_density_a_per_m2 * 3.14159265358979323846 * config_.channel_radius_m *
        config_.channel_radius_m;
    status_.channel_conductivity_s_per_m = channel_state_.conductivity_s_per_m;
}

void SurfaceDischargeArcAlgorithm::updateThermalState(double emitted_current_density, double dt)
{
    status_.cathode_temperature_k += emitted_current_density * dt * 2.0e2;
    status_.anode_temperature_k = anode_heating_model_.advanceTemperature(
        status_.anode_temperature_k, channel_state_.current_density_a_per_m2, dt);
}

bool SurfaceDischargeArcAlgorithm::initializePiccoreBoundaryCoupling()
{
    if (isVacuumArcDebugEnabled())
    {
        std::cerr << "[VACUUM_ARC_DEBUG] initializePiccoreBoundaryCoupling: begin" << std::endl;
    }

    Mesh::GeometryDimensions dims;
    dims.shape = Mesh::GeometryShape::PLATE;
    dims.length = std::max(2.0 * config_.channel_radius_m, 1.0e-4);
    dims.width = std::max(2.0 * config_.channel_radius_m, 1.0e-4);
    dims.thickness = std::max(config_.gap_distance_m, 1.0e-4);

    Mesh::MeshGenerationOptions options;
    options.nx = 1;
    options.ny = 1;
    options.nz = 6;
    options.tetrahedralize = true;

    auto generated_mesh = Mesh::GeometryMeshGenerator::generateFromDimensions(dims, options);
    if (!generated_mesh || !generated_mesh->validate())
    {
        if (isVacuumArcDebugEnabled())
        {
            std::cerr << "[VACUUM_ARC_DEBUG] initializePiccoreBoundaryCoupling: mesh generation/validation failed"
                      << std::endl;
        }
        return false;
    }

    if (isVacuumArcDebugEnabled())
    {
        std::cerr << "[VACUUM_ARC_DEBUG] mesh nodes=" << generated_mesh->getNodeCount()
                  << " elements=" << generated_mesh->getElementCount() << std::endl;
    }

    pic_surface_area_m2_ = std::max(dims.length * dims.width, 1.0e-12);
    pic_boundary_mesh_ = std::shared_ptr<Mesh::VolMesh>(std::move(generated_mesh));

    const double z_min = pic_boundary_mesh_->getBoundingBoxMin().z();
    const double z_max = pic_boundary_mesh_->getBoundingBoxMax().z();
    for (const auto& node : pic_boundary_mesh_->getNodes())
    {
        if (!node)
        {
            continue;
        }

        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12 || std::abs(z - z_max) < 1.0e-12)
        {
            node->setBoundaryType(Mesh::BoundaryType::DIRICHLET);
        }
        else
        {
            node->setBoundaryType(Mesh::BoundaryType::INTERIOR);
        }
    }

    pic_top_element_id_ = findTopElementId(*pic_boundary_mesh_);
    pic_boundary_particle_manager_ = std::make_shared<Particle::ParticleManager>();
    pic_boundary_poisson_solver_ = std::make_shared<FieldSolver::PoissonSolver>(pic_boundary_mesh_);

    PICcore::PICCycle::Parameters cycle_params;
    cycle_params.time_step = 1.0e-10;
    cycle_params.max_iterations = 128;
    cycle_params.statistics_frequency = 16;
    cycle_params.charge_conservation = true;
    cycle_params.verbose = false;

    pic_boundary_cycle_ = std::make_unique<PICcore::PICCycle>(cycle_params);
    pic_boundary_cycle_->setMesh(pic_boundary_mesh_->getNodes(), pic_boundary_mesh_->getElements());
    pic_boundary_cycle_->setParticleManager(pic_boundary_particle_manager_);
    pic_boundary_cycle_->setPoissonSolver(pic_boundary_poisson_solver_);

    const auto domain_min = pic_boundary_mesh_->getBoundingBoxMin();
    const auto domain_max = pic_boundary_mesh_->getBoundingBoxMax();
    std::vector<Geometry::Vector3D> periodic_vectors;
    if (domain_max.x() - domain_min.x() > 1.0e-12)
    {
        periodic_vectors.emplace_back(domain_max.x() - domain_min.x(), 0.0, 0.0);
    }
    if (domain_max.y() - domain_min.y() > 1.0e-12)
    {
        periodic_vectors.emplace_back(0.0, domain_max.y() - domain_min.y(), 0.0);
    }
    if (!periodic_vectors.empty())
    {
        auto periodic_boundary = std::make_shared<Particle::PeriodicBoundaryCondition>(periodic_vectors);
        pic_boundary_cycle_->setBoundaryCondition(BoundaryFace::XMin, periodic_boundary);
        pic_boundary_cycle_->setBoundaryCondition(BoundaryFace::XMax, periodic_boundary);
        pic_boundary_cycle_->setBoundaryCondition(BoundaryFace::YMin, periodic_boundary);
        pic_boundary_cycle_->setBoundaryCondition(BoundaryFace::YMax, periodic_boundary);
    }

    ArcFieldEmissionBoundaryCondition::Parameters emission_params;
    emission_params.emission_area_m2 = pic_surface_area_m2_;
    emission_params.max_emitted_particles_per_event = 8;
    emission_params.enable_secondary_electron_emission =
        config_.enable_secondary_electron_emission;
    emission_params.secondary_electron_yield_per_ion =
        std::max(0.0, config_.secondary_electron_yield_per_ion);
    emission_params.secondary_electron_speed_m_per_s =
        std::max(0.0, config_.secondary_electron_speed_m_per_s);
    pic_cathode_boundary_ = std::make_shared<ArcFieldEmissionBoundaryCondition>(emission_params);

    pic_boundary_cycle_->setBoundaryCondition(BoundaryFace::ZMin, pic_cathode_boundary_);
    pic_boundary_cycle_->setBoundaryCondition(
        BoundaryFace::ZMax, std::make_shared<Particle::AbsorbingBoundaryCondition>(1.0));

    PICcore::PICCycle::SurfaceBoundaryMetadata cathode_metadata;
    cathode_metadata.surface_id = 0;
    cathode_metadata.material_id = 0;
    cathode_metadata.circuit_node_id = 0;
    pic_boundary_cycle_->setBoundaryMetadata(BoundaryFace::ZMin, cathode_metadata);

    PICcore::PICCycle::SurfaceBoundaryMetadata anode_metadata;
    anode_metadata.surface_id = 1;
    anode_metadata.material_id = 0;
    anode_metadata.circuit_node_id = 1;
    pic_boundary_cycle_->setBoundaryMetadata(BoundaryFace::ZMax, anode_metadata);

    if (!pic_boundary_cycle_->initialize())
    {
        if (isVacuumArcDebugEnabled())
        {
            std::cerr << "[VACUUM_ARC_DEBUG] initializePiccoreBoundaryCoupling: PICCycle initialize failed"
                      << std::endl;
        }
        return false;
    }

    applyPiccoreBoundaryPotentials();
    pic_last_absorbed_electron_charge_c_ = 0.0;
    pic_last_absorbed_ion_charge_c_ = 0.0;
    pic_last_emitted_electron_charge_c_ = 0.0;

    if (isVacuumArcDebugEnabled())
    {
        std::cerr << "[VACUUM_ARC_DEBUG] initializePiccoreBoundaryCoupling: ready" << std::endl;
    }

    return true;
}

void SurfaceDischargeArcAlgorithm::applyPiccoreBoundaryPotentials()
{
    if (!pic_boundary_cycle_ || !pic_boundary_mesh_)
    {
        return;
    }

    std::vector<double> potential(pic_boundary_mesh_->getNodeCount(), 0.0);
    const double z_min = pic_boundary_mesh_->getBoundingBoxMin().z();
    const double z_max = pic_boundary_mesh_->getBoundingBoxMax().z();
    for (const auto& node : pic_boundary_mesh_->getNodes())
    {
        if (!node)
        {
            continue;
        }

        const std::size_t node_id = node->getId();
        if (node_id >= potential.size())
        {
            continue;
        }

        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12)
        {
            potential[node_id] = config_.surface_potential_v;
        }
        else if (std::abs(z - z_max) < 1.0e-12)
        {
            potential[node_id] = 0.0;
        }
    }

    pic_boundary_cycle_->setPotential(potential);
}

void SurfaceDischargeArcAlgorithm::seedBoundaryProbeParticles(double emitted_current_density,
                                                              double dt)
{
    if (!pic_boundary_mesh_ || !pic_boundary_particle_manager_ || !pic_boundary_cycle_ ||
        pic_boundary_mesh_->getElementCount() == 0)
    {
        return;
    }

    const auto& top_element = pic_boundary_mesh_->getElement(pic_top_element_id_);
    if (!top_element)
    {
        return;
    }

    const double elementary_charge = Basic::Constants::PhysicsConstants::ElementaryCharge;
    const double expected_from_emission =
        std::max(0.0, emitted_current_density) * pic_surface_area_m2_ * std::max(dt, 1.0e-12) /
        elementary_charge;
    const double scaled_probe = std::max(1.0, std::ceil(expected_from_emission * 0.002));
    const std::size_t probe_count =
        static_cast<std::size_t>(std::min(4.0, std::max(1.0, scaled_probe)));

    const auto seed_position = top_element->getCentroid();
    for (std::size_t i = 0; i < probe_count; ++i)
    {
        const Geometry::Point3D position(seed_position.x(), seed_position.y(),
                                         seed_position.z() - 1.0e-9 * static_cast<double>(i + 1));
        const auto particle_id = pic_boundary_particle_manager_->createIon(
            position, Geometry::Vector3D(0.0, 0.0, -2.0e6), 1, 1, 1.0);
        auto* particle = pic_boundary_particle_manager_->getParticle(particle_id);
        if (particle)
        {
            pic_boundary_cycle_->addParticleLocation(particle, pic_top_element_id_);
        }
    }
}

void SurfaceDischargeArcAlgorithm::resetPicBoundaryFeedback()
{
    status_.pic_absorbed_electron_current_a = 0.0;
    status_.pic_absorbed_ion_current_a = 0.0;
    status_.pic_emitted_electron_current_a = 0.0;
    status_.pic_boundary_net_current_a = 0.0;
    status_.pic_surface_charge_delta_c = 0.0;
}

void SurfaceDischargeArcAlgorithm::advancePiccoreBoundaryCoupling(double emitted_current_density,
                                                                   double dt)
{
    resetPicBoundaryFeedback();

    if (!piccore_boundary_coupling_initialized_ || !pic_boundary_cycle_ || !pic_cathode_boundary_ ||
        dt <= 0.0)
    {
        return;
    }

    if (isVacuumArcDebugEnabled())
    {
        std::cerr << "[VACUUM_ARC_DEBUG] advancePiccoreBoundaryCoupling: dt=" << dt
                  << " J_emit=" << emitted_current_density
                  << " particle_count_before=" << pic_boundary_particle_manager_->getContainer().size()
                  << std::endl;
    }

    pic_boundary_cycle_->setTimeStep(dt);
    pic_cathode_boundary_->setEmissionCurrentDensity(emitted_current_density);
    applyPiccoreBoundaryPotentials();
    seedBoundaryProbeParticles(emitted_current_density, dt);
    pic_boundary_cycle_->executeTimeStep();

    if (isVacuumArcDebugEnabled())
    {
        std::cerr << "[VACUUM_ARC_DEBUG] advancePiccoreBoundaryCoupling: particle_count_after_step="
                  << pic_boundary_particle_manager_->getContainer().size() << std::endl;
    }

    const auto& ledger = pic_boundary_cycle_->getSurfaceCurrentLedger();
    const auto& cathode_ledger = ledger[boundaryIndex(BoundaryFace::ZMin)];
    const auto& anode_ledger = ledger[boundaryIndex(BoundaryFace::ZMax)];

    const double absorbed_electron_charge_c =
        cathode_ledger.absorbed_electron_charge_c + anode_ledger.absorbed_electron_charge_c;
    const double absorbed_ion_charge_c =
        cathode_ledger.absorbed_ion_charge_c + anode_ledger.absorbed_ion_charge_c;
    const double emitted_electron_charge_c =
        cathode_ledger.emitted_electron_charge_c + anode_ledger.emitted_electron_charge_c;

    const double delta_absorbed_electron_charge_c =
        absorbed_electron_charge_c - pic_last_absorbed_electron_charge_c_;
    const double delta_absorbed_ion_charge_c =
        absorbed_ion_charge_c - pic_last_absorbed_ion_charge_c_;
    const double delta_emitted_electron_charge_c =
        emitted_electron_charge_c - pic_last_emitted_electron_charge_c_;

    pic_last_absorbed_electron_charge_c_ = absorbed_electron_charge_c;
    pic_last_absorbed_ion_charge_c_ = absorbed_ion_charge_c;
    pic_last_emitted_electron_charge_c_ = emitted_electron_charge_c;

    const double inv_dt = 1.0 / dt;
    status_.pic_absorbed_electron_current_a = delta_absorbed_electron_charge_c * inv_dt;
    status_.pic_absorbed_ion_current_a = delta_absorbed_ion_charge_c * inv_dt;
    status_.pic_emitted_electron_current_a = -delta_emitted_electron_charge_c * inv_dt;
    status_.pic_boundary_net_current_a = status_.pic_absorbed_electron_current_a +
                                         status_.pic_absorbed_ion_current_a +
                                         status_.pic_emitted_electron_current_a;

    status_.pic_surface_charge_delta_c = pic_boundary_cycle_->getAndResetSurfaceCharge();
    config_.surface_charge_density_c_per_m2 +=
        status_.pic_surface_charge_delta_c / std::max(pic_surface_area_m2_, 1.0e-12);

    pic_boundary_particle_manager_->removeInactiveParticles();

    if (isVacuumArcDebugEnabled())
    {
        std::cerr << "[VACUUM_ARC_DEBUG] advancePiccoreBoundaryCoupling: net_current="
                  << status_.pic_boundary_net_current_a
                  << " particle_count_after_cleanup="
                  << pic_boundary_particle_manager_->getContainer().size() << std::endl;
    }
}

double SurfaceDischargeArcAlgorithm::estimateCollisionVolumeM3() const
{
    constexpr double kPi = 3.14159265358979323846;
    const double cylindrical_volume = kPi * std::max(0.0, config_.channel_radius_m) *
                                      std::max(0.0, config_.channel_radius_m) *
                                      std::max(0.0, config_.gap_distance_m);

    if (!pic_boundary_mesh_)
    {
        return std::max(1.0e-18, cylindrical_volume);
    }

    return std::max({1.0e-18, cylindrical_volume, estimateMeshVolume(*pic_boundary_mesh_)});
}

void SurfaceDischargeArcAlgorithm::runCollisionStage(double dt)
{
    status_.collision_events_step = 0.0;
    status_.collision_particles_created_step = 0.0;
    status_.collision_ionization_events_step = 0.0;
    status_.collision_excitation_events_step = 0.0;
    status_.collision_charge_exchange_events_step = 0.0;
    status_.collision_ionization_fraction_step = 0.0;
    status_.collision_excitation_fraction_step = 0.0;
    status_.collision_charge_exchange_fraction_step = 0.0;
    status_.collision_reaction_weighted_emission_feedback_step = 1.0;
    status_.collision_reaction_weighted_channel_feedback_step = 1.0;
    status_.collision_event_rate_per_s = 0.0;
    status_.collision_stage_fallback_triggered = false;
    status_.collision_effective_neutral_density_m3 = 0.0;
    status_.collision_effective_ion_density_m3 = 0.0;
    status_.collision_effective_electron_density_m3 = 0.0;
    status_.neutral_outgassing_density_boost_m3 = dynamic_neutral_density_boost_m3_;
    status_.neutral_outgassing_reinjected_particles_step = 0.0;
    status_.neutral_outgassing_reinjection_rate_per_s = 0.0;
    collision_emission_feedback_multiplier_ = 1.0;
    collision_channel_feedback_multiplier_ = 1.0;

    if (!config_.enable_pic_mcc_collisions || dt <= 0.0 || !pic_boundary_particle_manager_)
    {
        return;
    }

    double source_density_boost_m3 = 0.0;
    if (config_.enable_neutral_outgassing_feedback)
    {
        const double gain = std::max(0.0, config_.neutral_outgassing_gain_m3_per_a);
        const double relaxation = std::max(0.0, config_.neutral_outgassing_relaxation_per_s);
        const double max_boost = std::max(0.0, config_.neutral_outgassing_max_density_boost_m3);

        source_density_boost_m3 = std::abs(status_.total_discharge_current_a) * gain * dt;
        dynamic_neutral_density_boost_m3_ += source_density_boost_m3;
        dynamic_neutral_density_boost_m3_ = std::max(0.0, dynamic_neutral_density_boost_m3_);

        if (relaxation > 0.0)
        {
            dynamic_neutral_density_boost_m3_ *= std::exp(-relaxation * dt);
        }
        if (max_boost > 0.0)
        {
            dynamic_neutral_density_boost_m3_ =
                std::min(dynamic_neutral_density_boost_m3_, max_boost);
        }
    }
    else
    {
        dynamic_neutral_density_boost_m3_ = 0.0;
    }

    status_.neutral_outgassing_density_boost_m3 = dynamic_neutral_density_boost_m3_;

    const double volume_m3 = estimateCollisionVolumeM3();
    injectOutgassingNeutralParticles(volume_m3, source_density_boost_m3, dt);
    const auto reconfigure_interval = std::max<std::size_t>(1, config_.collision_reconfigure_interval_steps);
    const bool reinitialize_configuration =
        !collision_handler_initialized_ || (collision_stage_step_counter_ % reconfigure_interval == 0);

    Collision::CollisionBackgroundDensities floor_densities{};
    floor_densities.neutral_density_m3 =
        std::max(0.0, config_.collision_neutral_density_floor_m3 + dynamic_neutral_density_boost_m3_);
    floor_densities.ion_density_m3 = std::max(0.0, config_.collision_ion_density_floor_m3);
    floor_densities.electron_density_m3 =
        std::max(0.0, config_.collision_electron_density_floor_m3);
    const std::string collision_set_id = config_.collision_cross_section_set_id.empty()
                                             ? config_.solver_config.collision_set
                                             : config_.collision_cross_section_set_id;
    const auto cross_section_set_id =
        Collision::parseCollisionCrossSectionSetId(collision_set_id)
            .value_or(Collision::CollisionCrossSectionSetId::BackgroundMccV1);
    Collision::BackgroundMccRuntimeConfig runtime_config;
    runtime_config.floor_densities = floor_densities;
    runtime_config.cross_section_set_id = cross_section_set_id;
    runtime_config.reinitialize_configuration = reinitialize_configuration;

    auto& collision_container = pic_boundary_particle_manager_->getContainer();
    std::vector<Particle::Particle> collision_container_snapshot;
    if (config_.collision_fallback_to_background_mcc_on_error)
    {
        collision_container_snapshot.assign(collision_container.begin(), collision_container.end());
    }

    auto step_result = Collision::executeBackgroundMccStep(
        collision_container, collision_handler_, dt, volume_m3, runtime_config);

    const auto primary_anomaly = detectCollisionStageAnomaly(
        step_result, dt, config_.collision_anomaly_event_rate_limit_per_s,
        config_.enable_collision_reaction_breakdown);
    if (primary_anomaly != CollisionStageAnomalyReason::None &&
        config_.collision_fallback_to_background_mcc_on_error)
    {
        status_.collision_stage_fallback_triggered = true;

        collision_container.clear();
        if (!collision_container_snapshot.empty())
        {
            collision_container.addParticles(collision_container_snapshot);
        }

        collision_handler_.resetStatistics();

        Collision::CollisionBackgroundDensities conservative_fallback_densities{};
        conservative_fallback_densities.neutral_density_m3 =
            std::max(0.0, config_.collision_neutral_density_floor_m3);
        conservative_fallback_densities.ion_density_m3 =
            std::max(0.0, config_.collision_ion_density_floor_m3);
        conservative_fallback_densities.electron_density_m3 =
            std::max(0.0, config_.collision_electron_density_floor_m3);
        Collision::BackgroundMccRuntimeConfig fallback_runtime_config;
        fallback_runtime_config.floor_densities = conservative_fallback_densities;
        fallback_runtime_config.cross_section_set_id =
            Collision::CollisionCrossSectionSetId::BackgroundMccV1;
        fallback_runtime_config.reinitialize_configuration = true;

        step_result = Collision::executeBackgroundMccStep(
            collision_container, collision_handler_, dt, volume_m3, fallback_runtime_config);

        const auto fallback_anomaly =
            detectCollisionStageAnomaly(step_result, dt, std::numeric_limits<double>::infinity(),
                                        config_.enable_collision_reaction_breakdown);
        if (fallback_anomaly != CollisionStageAnomalyReason::None)
        {
            step_result = {};
            step_result.effective_neutral_density_m3 =
                conservative_fallback_densities.neutral_density_m3;
            step_result.effective_ion_density_m3 = conservative_fallback_densities.ion_density_m3;
            step_result.effective_electron_density_m3 =
                conservative_fallback_densities.electron_density_m3;
        }
    }

    collision_handler_initialized_ = true;
    ++collision_stage_step_counter_;

    const auto& stats = collision_handler_.getStatistics();
    collision_last_total_events_ = static_cast<std::size_t>(stats.total_collisions);
    collision_last_total_particles_created_ =
        static_cast<std::size_t>(stats.particles_created);

    status_.collision_events_step = static_cast<double>(step_result.collision_events);
    status_.collision_particles_created_step =
        static_cast<double>(step_result.generated_particles);

    if (config_.enable_collision_reaction_breakdown)
    {
        status_.collision_ionization_events_step =
            static_cast<double>(step_result.ionization_events);
        status_.collision_excitation_events_step =
            static_cast<double>(step_result.excitation_events);
        status_.collision_charge_exchange_events_step =
            static_cast<double>(step_result.charge_exchange_events);
    }

    const double safe_dt = std::max(dt, 1.0e-18);
    status_.collision_event_rate_per_s = status_.collision_events_step / safe_dt;

    status_.collision_effective_neutral_density_m3 = step_result.effective_neutral_density_m3;
    status_.collision_effective_ion_density_m3 = step_result.effective_ion_density_m3;
    status_.collision_effective_electron_density_m3 = step_result.effective_electron_density_m3;

    const auto reaction_fractions =
        computeCollisionReactionFractions(step_result, config_.enable_collision_reaction_breakdown);
    status_.collision_ionization_fraction_step = reaction_fractions.ionization_fraction;
    status_.collision_excitation_fraction_step = reaction_fractions.excitation_fraction;
    status_.collision_charge_exchange_fraction_step = reaction_fractions.charge_exchange_fraction;

    const double electron_reference_density =
        std::max(1.0, config_.collision_electron_density_floor_m3);
    const double ion_reference_density = std::max(1.0, config_.collision_ion_density_floor_m3);
    const double emission_density_ratio = std::max(
        0.0, status_.collision_effective_electron_density_m3 / electron_reference_density - 1.0);
    const double channel_density_ratio =
        std::max(0.0, status_.collision_effective_ion_density_m3 / ion_reference_density - 1.0);

    const double base_emission_feedback_delta =
        config_.collision_emission_feedback_gain * emission_density_ratio;
    const double base_channel_feedback_delta =
        config_.collision_channel_feedback_gain * channel_density_ratio;
    const double reaction_emission_feedback_delta = computeReactionGainDelta(
        reaction_fractions, config_.collision_ionization_emission_feedback_gain,
        config_.collision_excitation_emission_feedback_gain,
        config_.collision_charge_exchange_emission_feedback_gain);
    const double reaction_channel_feedback_delta = computeReactionGainDelta(
        reaction_fractions, config_.collision_ionization_channel_feedback_gain,
        config_.collision_excitation_channel_feedback_gain,
        config_.collision_charge_exchange_channel_feedback_gain);

    collision_emission_feedback_multiplier_ =
        std::max(0.0, 1.0 + base_emission_feedback_delta + reaction_emission_feedback_delta);
    collision_channel_feedback_multiplier_ =
        std::max(0.0, 1.0 + base_channel_feedback_delta + reaction_channel_feedback_delta);
    status_.collision_reaction_weighted_emission_feedback_step =
        collision_emission_feedback_multiplier_;
    status_.collision_reaction_weighted_channel_feedback_step =
        collision_channel_feedback_multiplier_;
}

void SurfaceDischargeArcAlgorithm::injectOutgassingNeutralParticles(double volume_m3,
                                                                    double source_density_boost_m3,
                                                                    double dt)
{
    if (!config_.enable_neutral_outgassing_feedback ||
        !config_.enable_neutral_outgassing_reinjection ||
        !pic_boundary_particle_manager_ || !pic_boundary_mesh_ ||
        !(volume_m3 > 0.0) || !(dt > 0.0) || !(source_density_boost_m3 > 0.0))
    {
        return;
    }

    const double reinjection_budget_particles =
        source_density_boost_m3 * std::max(0.0, config_.neutral_outgassing_reinjection_gain);
    if (!(reinjection_budget_particles > 0.0) || !std::isfinite(reinjection_budget_particles))
    {
        return;
    }

    const double macro_weight = std::max(1.0, config_.neutral_outgassing_reinjection_macro_weight);
    pending_outgassing_reinjection_particles_ += reinjection_budget_particles / macro_weight;

    std::size_t macro_particles = static_cast<std::size_t>(
        std::floor(std::max(0.0, pending_outgassing_reinjection_particles_)));
    macro_particles = std::min(macro_particles,
        std::max<std::size_t>(1, config_.neutral_outgassing_reinjection_max_particles_per_step));
    if (macro_particles == 0)
    {
        return;
    }
    pending_outgassing_reinjection_particles_ = std::max(
        0.0, pending_outgassing_reinjection_particles_ - static_cast<double>(macro_particles));

    const auto min_corner = pic_boundary_mesh_->getBoundingBoxMin();
    const auto max_corner = pic_boundary_mesh_->getBoundingBoxMax();
    std::uniform_real_distribution<double> random_x(min_corner.x(), max_corner.x());
    std::uniform_real_distribution<double> random_y(min_corner.y(), max_corner.y());
    std::uniform_real_distribution<double> random01(0.0, 1.0);
    constexpr double kTwoPi = 6.28318530717958647692;
    const double speed_m_per_s =
        std::max(0.0, config_.neutral_outgassing_reinjection_speed_m_per_s);
    const double neutral_mass_kg =
        std::max(1.0, config_.neutral_outgassing_reinjection_mass_amu) *
        Basic::Constants::PhysicsConstants::ProtonMass;

    std::vector<Particle::Particle> reinjected_particles;
    reinjected_particles.reserve(macro_particles);

    for (std::size_t i = 0; i < macro_particles; ++i)
    {
        const double phi = kTwoPi * random01(neutral_outgassing_reinjection_rng_);
        const double cos_theta = random01(neutral_outgassing_reinjection_rng_);
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

        const SCDAT::Geometry::Vector3D velocity(
            speed_m_per_s * sin_theta * std::cos(phi),
            speed_m_per_s * sin_theta * std::sin(phi),
            speed_m_per_s * cos_theta);

        const SCDAT::Geometry::Point3D position(
            random_x(neutral_outgassing_reinjection_rng_),
            random_y(neutral_outgassing_reinjection_rng_),
            min_corner.z() + 1.0e-9 * static_cast<double>(i + 1));

        Particle::Particle neutral(
            Particle::ParticleFactory::getNextId(),
            Particle::ParticleType::ATOM,
            position,
            velocity,
            neutral_mass_kg,
            0.0,
            macro_weight);
        reinjected_particles.push_back(neutral);
    }

    if (!reinjected_particles.empty())
    {
        pic_boundary_particle_manager_->getContainer().addParticles(reinjected_particles);
        status_.neutral_outgassing_reinjected_particles_step =
            static_cast<double>(reinjected_particles.size());
        status_.neutral_outgassing_reinjection_rate_per_s =
            status_.neutral_outgassing_reinjected_particles_step / dt;
    }
}

void SurfaceDischargeArcAlgorithm::updateSurfaceCircuitState(double dt)
{
    status_.surface_circuit_drive_current_a = status_.pic_boundary_net_current_a;
    status_.surface_circuit_leak_current_a = 0.0;
    status_.surface_load_branch_current_a = 0.0;
    status_.surface_circuit_drive_power_w = 0.0;
    status_.surface_circuit_leak_power_w = 0.0;
    status_.surface_load_branch_power_w = 0.0;
    status_.surface_load_resistive_power_w = 0.0;
    status_.surface_load_inductive_energy_j = 0.0;
    status_.surface_circuit_power_balance_error_w = 0.0;
    status_.charge_conservation_error_c = 0.0;

    const double reference_potential_v = config_.surface_reference_potential_v;

    if (dt > 0.0)
    {
        const double previous_potential_v = config_.surface_potential_v;
        const double previous_charge_c = surface_circuit_charge_c_;
        const double previous_branch_current_a = surface_load_branch_current_a_;
        const double previous_relative_potential_v = previous_potential_v - reference_potential_v;

        const double capacitance_f = std::max(config_.surface_capacitance_f, 1.0e-18);
        const double leakage_conductance_s = std::max(0.0, config_.surface_leakage_conductance_s);
        const double leakage_current_a = leakage_conductance_s * previous_relative_potential_v;

        double load_branch_current_a = 0.0;
        double load_resistive_power_w = 0.0;
        double load_inductive_energy_j = 0.0;

        if (config_.surface_load_model == SurfaceCircuitLoadModel::ResistiveShunt)
        {
            surface_load_branch_current_a_ = 0.0;
            const double resistance_ohm = std::max(1.0e-12, config_.surface_load_resistance_ohm);
            load_branch_current_a = previous_relative_potential_v / resistance_ohm;
            load_resistive_power_w =
                load_branch_current_a * load_branch_current_a * resistance_ohm;
        }
        else if (config_.surface_load_model == SurfaceCircuitLoadModel::RlBranch)
        {
            const double resistance_ohm = std::max(1.0e-12, config_.surface_load_resistance_ohm);
            const double inductance_h = std::max(1.0e-18, config_.surface_load_inductance_h);
            const double drive_voltage_v = previous_relative_potential_v;
            const double di_dt =
                (drive_voltage_v - resistance_ohm * surface_load_branch_current_a_) / inductance_h;
            surface_load_branch_current_a_ += di_dt * dt;

            const double current_limit_a = std::max(0.0, config_.surface_load_current_limit_a);
            if (current_limit_a > 0.0)
            {
                surface_load_branch_current_a_ =
                    std::clamp(surface_load_branch_current_a_, -current_limit_a, current_limit_a);
            }
            load_branch_current_a = surface_load_branch_current_a_;
            load_resistive_power_w =
                load_branch_current_a * load_branch_current_a * resistance_ohm;
            load_inductive_energy_j =
                0.5 * inductance_h * load_branch_current_a * load_branch_current_a;
        }
        else
        {
            surface_load_branch_current_a_ = 0.0;
        }

        status_.surface_load_branch_current_a = load_branch_current_a;
        status_.surface_load_resistive_power_w = load_resistive_power_w;
        status_.surface_load_inductive_energy_j = load_inductive_energy_j;

        const double net_circuit_current_a =
            status_.surface_circuit_drive_current_a - leakage_current_a - load_branch_current_a;

        const double potential_delta_v = (net_circuit_current_a * dt) / capacitance_f;
        const double max_step_v = std::max(0.0, config_.max_surface_potential_step_v);
        const double clamped_potential_delta_v =
            std::clamp(potential_delta_v, -max_step_v, max_step_v);

        config_.surface_potential_v += clamped_potential_delta_v;
        status_.surface_circuit_leak_current_a = leakage_current_a;

        const double expected_charge_delta_c = net_circuit_current_a * dt;
        const double expected_charge_c = previous_charge_c + expected_charge_delta_c;
        const double updated_charge_c = capacitance_f *
            (config_.surface_potential_v - reference_potential_v);
        status_.charge_conservation_error_c = updated_charge_c - expected_charge_c;

        bool rolled_back = false;

        if (!std::isfinite(config_.surface_potential_v) ||
            std::abs(config_.surface_potential_v) >
                std::max(1.0, config_.anomaly_surface_potential_limit_v))
        {
            config_.surface_potential_v = previous_potential_v;
            surface_circuit_charge_c_ = previous_charge_c;
            surface_load_branch_current_a_ = previous_branch_current_a;
            status_.surface_load_branch_current_a = previous_branch_current_a;
            status_.charge_conservation_error_c = 0.0;
            ++status_.stability_rollbacks;
            rolled_back = true;
        }

        const double updated_relative_potential_v = config_.surface_potential_v - reference_potential_v;
        status_.surface_circuit_drive_power_w =
            previous_relative_potential_v * status_.surface_circuit_drive_current_a;
        status_.surface_circuit_leak_power_w =
            previous_relative_potential_v * leakage_current_a;
        status_.surface_load_branch_power_w =
            previous_relative_potential_v * status_.surface_load_branch_current_a;

        const double capacitor_power_w =
            0.5 * (previous_relative_potential_v + updated_relative_potential_v) *
            net_circuit_current_a;
        status_.surface_circuit_power_balance_error_w =
            status_.surface_circuit_drive_power_w -
            status_.surface_circuit_leak_power_w -
            status_.surface_load_branch_power_w -
            capacitor_power_w;

        if (rolled_back)
        {
            status_.surface_circuit_power_balance_error_w = 0.0;
        }

        accumulated_drive_charge_c_ += status_.surface_circuit_drive_current_a * dt;
        accumulated_leak_charge_c_ += leakage_current_a * dt;
        accumulated_load_charge_c_ += load_branch_current_a * dt;
        accumulated_drive_energy_j_ += status_.surface_circuit_drive_power_w * dt;
        accumulated_leak_dissipated_energy_j_ +=
            std::max(0.0, status_.surface_circuit_leak_power_w) * dt;
        accumulated_load_dissipated_energy_j_ +=
            std::max(0.0, status_.surface_load_resistive_power_w) * dt;
        accumulated_power_balance_error_j_ +=
            status_.surface_circuit_power_balance_error_w * dt;
    }

    const double effective_capacitance_f = std::max(config_.surface_capacitance_f, 1.0e-18);
    const double relative_potential_v = config_.surface_potential_v - reference_potential_v;
    surface_circuit_charge_c_ =
        effective_capacitance_f * relative_potential_v;

    status_.surface_relative_potential_v = relative_potential_v;
    status_.surface_circuit_charge_c = surface_circuit_charge_c_;
    status_.surface_circuit_capacitor_energy_j =
        0.5 * effective_capacitance_f * relative_potential_v * relative_potential_v;
    if (config_.surface_load_model == SurfaceCircuitLoadModel::RlBranch)
    {
        const double inductance_h = std::max(1.0e-18, config_.surface_load_inductance_h);
        status_.surface_load_inductive_energy_j =
            0.5 * inductance_h * surface_load_branch_current_a_ * surface_load_branch_current_a_;
    }
    else
    {
        status_.surface_load_inductive_energy_j = 0.0;
    }
    status_.surface_circuit_drive_energy_j = accumulated_drive_energy_j_;
    status_.surface_circuit_leak_dissipated_energy_j = accumulated_leak_dissipated_energy_j_;
    status_.surface_load_dissipated_energy_j = accumulated_load_dissipated_energy_j_;
    status_.surface_circuit_energy_balance_error_j = accumulated_power_balance_error_j_;
    status_.surface_potential_v = config_.surface_potential_v;
    status_.surface_charge_density_c_per_m2 = config_.surface_charge_density_c_per_m2;
}

void SurfaceDischargeArcAlgorithm::updateResidualMonitoring(double effective_field, double dt)
{
    status_.field_residual_v_per_m = 0.0;
    status_.particle_residual_a = 0.0;
    status_.circuit_residual_c = status_.charge_conservation_error_c;
    status_.field_residual_filtered_v_per_m = 0.0;
    status_.particle_residual_filtered_a = 0.0;
    status_.circuit_residual_filtered_c = status_.circuit_residual_c;
    status_.residual_alarm_counter = residual_alarm_consecutive_counter_;
    status_.residual_alarm_clear_counter = residual_alarm_clear_consecutive_counter_;

    if (!config_.enable_residual_monitoring || dt <= 0.0)
    {
        residual_alarm_consecutive_counter_ = 0;
        residual_alarm_clear_consecutive_counter_ = 0;
        residual_filter_initialized_ = false;
        filtered_field_residual_state_ = 0.0;
        filtered_particle_residual_state_ = 0.0;
        filtered_circuit_residual_state_ = 0.0;
        status_.residual_alarm_active = false;
        status_.residual_alarm_counter = 0;
        status_.residual_alarm_clear_counter = 0;
        return;
    }

    const double conductivity = std::max(1.0e-18, channel_state_.conductivity_s_per_m);
    const double ohmic_field_v_per_m = channel_state_.current_density_a_per_m2 / conductivity;
    status_.field_residual_v_per_m = effective_field - ohmic_field_v_per_m;

    const double generated_particle_current_proxy_a =
        status_.collision_particles_created_step *
        Basic::Constants::PhysicsConstants::ElementaryCharge / std::max(dt, 1.0e-18);
    status_.particle_residual_a =
        status_.pic_boundary_net_current_a - generated_particle_current_proxy_a;

    if (config_.enable_residual_ema_filter)
    {
        const double alpha = std::clamp(config_.residual_ema_alpha, 0.0, 1.0);
        if (!residual_filter_initialized_)
        {
            filtered_field_residual_state_ = status_.field_residual_v_per_m;
            filtered_particle_residual_state_ = status_.particle_residual_a;
            filtered_circuit_residual_state_ = status_.circuit_residual_c;
            residual_filter_initialized_ = true;
        }
        else
        {
            const double keep = 1.0 - alpha;
            filtered_field_residual_state_ =
                alpha * status_.field_residual_v_per_m + keep * filtered_field_residual_state_;
            filtered_particle_residual_state_ =
                alpha * status_.particle_residual_a + keep * filtered_particle_residual_state_;
            filtered_circuit_residual_state_ =
                alpha * status_.circuit_residual_c + keep * filtered_circuit_residual_state_;
        }
    }
    else
    {
        residual_filter_initialized_ = false;
        filtered_field_residual_state_ = status_.field_residual_v_per_m;
        filtered_particle_residual_state_ = status_.particle_residual_a;
        filtered_circuit_residual_state_ = status_.circuit_residual_c;
    }

    status_.field_residual_filtered_v_per_m = filtered_field_residual_state_;
    status_.particle_residual_filtered_a = filtered_particle_residual_state_;
    status_.circuit_residual_filtered_c = filtered_circuit_residual_state_;

    const bool field_alarm =
        std::abs(status_.field_residual_filtered_v_per_m) >
        std::max(0.0, config_.residual_field_tolerance_v_per_m);
    const bool particle_alarm =
        std::abs(status_.particle_residual_filtered_a) >
        std::max(0.0, config_.residual_particle_current_tolerance_a);
    const bool circuit_alarm =
        std::abs(status_.circuit_residual_filtered_c) >
        std::max(0.0, config_.residual_circuit_charge_tolerance_c);

    const bool residual_out_of_band = field_alarm || particle_alarm || circuit_alarm;
    if (residual_out_of_band)
    {
        ++residual_alarm_consecutive_counter_;
        residual_alarm_clear_consecutive_counter_ = 0;
    }
    else
    {
        residual_alarm_consecutive_counter_ = 0;
        ++residual_alarm_clear_consecutive_counter_;
    }

    const std::size_t trigger_steps =
        std::max<std::size_t>(1, config_.residual_alarm_consecutive_steps);
    const std::size_t clear_steps =
        std::max<std::size_t>(1, config_.residual_alarm_clear_consecutive_steps);

    if (!status_.residual_alarm_active)
    {
        status_.residual_alarm_active = residual_alarm_consecutive_counter_ >= trigger_steps;
    }
    else if (!residual_out_of_band &&
             residual_alarm_clear_consecutive_counter_ >= clear_steps)
    {
        status_.residual_alarm_active = false;
    }

    status_.residual_alarm_counter = residual_alarm_consecutive_counter_;
    status_.residual_alarm_clear_counter = residual_alarm_clear_consecutive_counter_;
}

bool SurfaceDischargeArcAlgorithm::hasAnomalousState() const
{
    if (!std::isfinite(status_.surface_potential_v) ||
        !std::isfinite(status_.total_discharge_current_a) ||
        !std::isfinite(channel_state_.current_density_a_per_m2))
    {
        return true;
    }

    if (std::abs(channel_state_.current_density_a_per_m2) >
        std::max(1.0, config_.anomaly_current_density_limit_a_per_m2))
    {
        return true;
    }

    if (std::abs(status_.surface_potential_v) >
        std::max(1.0, config_.anomaly_surface_potential_limit_v))
    {
        return true;
    }

    return false;
}

void SurfaceDischargeArcAlgorithm::appendHistory()
{
    history_time_.push_back(status_.current_time_s);
    history_current_.push_back(status_.total_discharge_current_a);
    history_current_density_.push_back(channel_state_.current_density_a_per_m2);
    history_cathode_temperature_.push_back(status_.cathode_temperature_k);
    history_anode_temperature_.push_back(status_.anode_temperature_k);
    history_conductivity_.push_back(status_.channel_conductivity_s_per_m);
    history_pic_boundary_net_current_.push_back(status_.pic_boundary_net_current_a);
    history_pic_surface_charge_delta_.push_back(status_.pic_surface_charge_delta_c);
    history_surface_potential_.push_back(status_.surface_potential_v);
    history_surface_charge_density_.push_back(status_.surface_charge_density_c_per_m2);
    history_surface_circuit_drive_current_.push_back(status_.surface_circuit_drive_current_a);
    history_surface_circuit_leak_current_.push_back(status_.surface_circuit_leak_current_a);
    history_surface_load_branch_current_.push_back(status_.surface_load_branch_current_a);
    history_surface_relative_potential_.push_back(status_.surface_relative_potential_v);
    history_surface_circuit_capacitor_energy_.push_back(status_.surface_circuit_capacitor_energy_j);
    history_surface_circuit_drive_power_.push_back(status_.surface_circuit_drive_power_w);
    history_surface_circuit_leak_power_.push_back(status_.surface_circuit_leak_power_w);
    history_surface_load_branch_power_.push_back(status_.surface_load_branch_power_w);
    history_surface_load_resistive_power_.push_back(status_.surface_load_resistive_power_w);
    history_surface_load_inductive_energy_.push_back(status_.surface_load_inductive_energy_j);
    history_surface_circuit_power_balance_error_.push_back(
        status_.surface_circuit_power_balance_error_w);
    history_surface_circuit_drive_energy_.push_back(status_.surface_circuit_drive_energy_j);
    history_surface_circuit_leak_dissipated_energy_.push_back(
        status_.surface_circuit_leak_dissipated_energy_j);
    history_surface_load_dissipated_energy_.push_back(status_.surface_load_dissipated_energy_j);
    history_surface_circuit_energy_balance_error_.push_back(
        status_.surface_circuit_energy_balance_error_j);
    history_charge_conservation_error_.push_back(status_.charge_conservation_error_c);
    history_collision_events_.push_back(status_.collision_events_step);
    history_collision_particles_created_.push_back(status_.collision_particles_created_step);
    history_collision_ionization_events_.push_back(status_.collision_ionization_events_step);
    history_collision_excitation_events_.push_back(status_.collision_excitation_events_step);
    history_collision_charge_exchange_events_.push_back(status_.collision_charge_exchange_events_step);
    history_collision_ionization_fraction_.push_back(status_.collision_ionization_fraction_step);
    history_collision_excitation_fraction_.push_back(status_.collision_excitation_fraction_step);
    history_collision_charge_exchange_fraction_.push_back(
        status_.collision_charge_exchange_fraction_step);
    history_collision_emission_feedback_multiplier_.push_back(
        status_.collision_reaction_weighted_emission_feedback_step);
    history_collision_channel_feedback_multiplier_.push_back(
        status_.collision_reaction_weighted_channel_feedback_step);
    history_collision_event_rate_.push_back(status_.collision_event_rate_per_s);
    history_collision_fallback_triggered_.push_back(
        status_.collision_stage_fallback_triggered ? 1.0 : 0.0);
    history_collision_neutral_density_.push_back(status_.collision_effective_neutral_density_m3);
    history_collision_ion_density_.push_back(status_.collision_effective_ion_density_m3);
    history_collision_electron_density_.push_back(status_.collision_effective_electron_density_m3);
    history_neutral_outgassing_boost_.push_back(status_.neutral_outgassing_density_boost_m3);
    history_neutral_outgassing_reinjected_particles_.push_back(
        status_.neutral_outgassing_reinjected_particles_step);
    history_neutral_outgassing_reinjection_rate_.push_back(
        status_.neutral_outgassing_reinjection_rate_per_s);
    history_field_residual_.push_back(status_.field_residual_v_per_m);
    history_particle_residual_.push_back(status_.particle_residual_a);
    history_circuit_residual_.push_back(status_.circuit_residual_c);
    history_field_residual_filtered_.push_back(status_.field_residual_filtered_v_per_m);
    history_particle_residual_filtered_.push_back(status_.particle_residual_filtered_a);
    history_circuit_residual_filtered_.push_back(status_.circuit_residual_filtered_c);
    history_residual_alarm_active_.push_back(status_.residual_alarm_active ? 1.0 : 0.0);
    history_residual_alarm_counter_.push_back(static_cast<double>(status_.residual_alarm_counter));
    history_residual_alarm_clear_counter_.push_back(
        static_cast<double>(status_.residual_alarm_clear_counter));
    history_stability_substeps_.push_back(static_cast<double>(status_.stability_substeps_used));
    history_stability_rollbacks_.push_back(static_cast<double>(status_.stability_rollbacks));
    history_stability_effective_substep_.push_back(status_.stability_effective_substep_s);
    history_stability_adaptive_scale_.push_back(status_.stability_adaptive_scale);
    history_stability_adaptive_reductions_.push_back(
        static_cast<double>(status_.stability_adaptive_reductions));
    history_stability_anomaly_isolation_triggered_.push_back(
        status_.stability_anomaly_isolation_triggered ? 1.0 : 0.0);
    history_stability_anomaly_isolation_events_.push_back(
        static_cast<double>(status_.stability_anomaly_isolation_events));
}

bool SurfaceDischargeArcAlgorithm::initialize(const DischargeConfiguration& config)
{
    config_ = config;
    neutral_outgassing_reinjection_rng_.seed(config_.seed);
    if (!config_.solver_config.collision_set.empty())
    {
        config_.collision_cross_section_set_id = config_.solver_config.collision_set;
    }
    const bool residual_guarded_policy_requested =
        Solver::isResidualGuardedSolverPolicyToken(config_.solver_config.convergence_policy);
    if (residual_guarded_policy_requested &&
        !config_.enable_adaptive_internal_timestep)
    {
        config_.enable_adaptive_internal_timestep = true;
    }
    plasma_channel_model_.setParameters(config_.channel_parameters);
    integrator_.setParameters(config_.integrator_parameters);
    emission_model_.setStrategyType(config_.emission_strategy);
    emission_model_.configureStrategy(config_.emission_parameters);
    status_ = DischargeStatus{};
    status_.cathode_temperature_k = config.cathode_temperature_k;
    status_.anode_temperature_k = config.anode_temperature_k;
    initialized_ = true;
    history_time_.clear();
    history_current_.clear();
    history_current_density_.clear();
    history_cathode_temperature_.clear();
    history_anode_temperature_.clear();
    history_conductivity_.clear();
    history_pic_boundary_net_current_.clear();
    history_pic_surface_charge_delta_.clear();
    history_surface_potential_.clear();
    history_surface_charge_density_.clear();
    history_surface_circuit_drive_current_.clear();
    history_surface_circuit_leak_current_.clear();
    history_surface_load_branch_current_.clear();
    history_surface_relative_potential_.clear();
    history_surface_circuit_capacitor_energy_.clear();
    history_surface_circuit_drive_power_.clear();
    history_surface_circuit_leak_power_.clear();
    history_surface_load_branch_power_.clear();
    history_surface_load_resistive_power_.clear();
    history_surface_load_inductive_energy_.clear();
    history_surface_circuit_power_balance_error_.clear();
    history_surface_circuit_drive_energy_.clear();
    history_surface_circuit_leak_dissipated_energy_.clear();
    history_surface_load_dissipated_energy_.clear();
    history_surface_circuit_energy_balance_error_.clear();
    history_charge_conservation_error_.clear();
    history_collision_events_.clear();
    history_collision_particles_created_.clear();
    history_collision_ionization_events_.clear();
    history_collision_excitation_events_.clear();
    history_collision_charge_exchange_events_.clear();
    history_collision_ionization_fraction_.clear();
    history_collision_excitation_fraction_.clear();
    history_collision_charge_exchange_fraction_.clear();
    history_collision_emission_feedback_multiplier_.clear();
    history_collision_channel_feedback_multiplier_.clear();
    history_collision_event_rate_.clear();
    history_collision_fallback_triggered_.clear();
    history_collision_neutral_density_.clear();
    history_collision_ion_density_.clear();
    history_collision_electron_density_.clear();
    history_neutral_outgassing_boost_.clear();
    history_neutral_outgassing_reinjected_particles_.clear();
    history_neutral_outgassing_reinjection_rate_.clear();
    history_field_residual_.clear();
    history_particle_residual_.clear();
    history_circuit_residual_.clear();
    history_field_residual_filtered_.clear();
    history_particle_residual_filtered_.clear();
    history_circuit_residual_filtered_.clear();
    history_residual_alarm_active_.clear();
    history_residual_alarm_counter_.clear();
    history_residual_alarm_clear_counter_.clear();
    history_stability_substeps_.clear();
    history_stability_rollbacks_.clear();
    history_stability_effective_substep_.clear();
    history_stability_adaptive_scale_.clear();
    history_stability_adaptive_reductions_.clear();
    history_stability_anomaly_isolation_triggered_.clear();
    history_stability_anomaly_isolation_events_.clear();

    collision_handler_.resetStatistics();
    collision_handler_initialized_ = false;
    collision_stage_step_counter_ = 0;
    collision_last_total_events_ = 0;
    collision_last_total_particles_created_ = 0;

    status_.surface_potential_v = config_.surface_potential_v;
    status_.surface_relative_potential_v =
        config_.surface_potential_v - config_.surface_reference_potential_v;
    status_.surface_charge_density_c_per_m2 = config_.surface_charge_density_c_per_m2;
    status_.alignment_mode = config_.alignment_mode;
    surface_circuit_charge_c_ = std::max(config_.surface_capacitance_f, 1.0e-18) *
                                (config_.surface_potential_v - config_.surface_reference_potential_v);
    status_.surface_circuit_charge_c = surface_circuit_charge_c_;
    status_.surface_circuit_capacitor_energy_j =
        0.5 * std::max(config_.surface_capacitance_f, 1.0e-18) *
        status_.surface_relative_potential_v * status_.surface_relative_potential_v;
    status_.collision_events_step = 0.0;
    status_.collision_particles_created_step = 0.0;
    status_.collision_ionization_events_step = 0.0;
    status_.collision_excitation_events_step = 0.0;
    status_.collision_charge_exchange_events_step = 0.0;
    status_.collision_ionization_fraction_step = 0.0;
    status_.collision_excitation_fraction_step = 0.0;
    status_.collision_charge_exchange_fraction_step = 0.0;
    status_.collision_reaction_weighted_emission_feedback_step = 1.0;
    status_.collision_reaction_weighted_channel_feedback_step = 1.0;
    status_.collision_event_rate_per_s = 0.0;
    status_.collision_stage_fallback_triggered = false;
    status_.collision_effective_neutral_density_m3 = 0.0;
    status_.collision_effective_ion_density_m3 = 0.0;
    status_.collision_effective_electron_density_m3 = 0.0;
    status_.neutral_outgassing_density_boost_m3 = 0.0;
    status_.neutral_outgassing_reinjected_particles_step = 0.0;
    status_.neutral_outgassing_reinjection_rate_per_s = 0.0;
    status_.field_residual_v_per_m = 0.0;
    status_.particle_residual_a = 0.0;
    status_.circuit_residual_c = 0.0;
    status_.field_residual_filtered_v_per_m = 0.0;
    status_.particle_residual_filtered_a = 0.0;
    status_.circuit_residual_filtered_c = 0.0;
    status_.residual_alarm_active = false;
    status_.residual_alarm_counter = 0;
    status_.residual_alarm_clear_counter = 0;
    status_.surface_load_branch_current_a = 0.0;
    status_.surface_circuit_drive_power_w = 0.0;
    status_.surface_circuit_leak_power_w = 0.0;
    status_.surface_load_branch_power_w = 0.0;
    status_.surface_load_resistive_power_w = 0.0;
    status_.surface_load_inductive_energy_j = 0.0;
    status_.surface_circuit_power_balance_error_w = 0.0;
    status_.surface_circuit_drive_energy_j = 0.0;
    status_.surface_circuit_leak_dissipated_energy_j = 0.0;
    status_.surface_load_dissipated_energy_j = 0.0;
    status_.surface_circuit_energy_balance_error_j = 0.0;
    status_.charge_conservation_error_c = 0.0;
    status_.stability_substeps_used = 0;
    status_.stability_rollbacks = 0;
    status_.stability_effective_substep_s = 0.0;
    status_.stability_adaptive_scale = 1.0;
    status_.stability_adaptive_reductions = 0;
    status_.stability_anomaly_isolation_triggered = false;
    status_.stability_anomaly_isolation_events = 0;

    surface_load_branch_current_a_ = 0.0;
    accumulated_drive_charge_c_ = 0.0;
    accumulated_leak_charge_c_ = 0.0;
    accumulated_load_charge_c_ = 0.0;
    accumulated_drive_energy_j_ = 0.0;
    accumulated_leak_dissipated_energy_j_ = 0.0;
    accumulated_load_dissipated_energy_j_ = 0.0;
    accumulated_power_balance_error_j_ = 0.0;
    residual_alarm_consecutive_counter_ = 0;
    residual_alarm_clear_consecutive_counter_ = 0;
    residual_filter_initialized_ = false;
    filtered_field_residual_state_ = 0.0;
    filtered_particle_residual_state_ = 0.0;
    filtered_circuit_residual_state_ = 0.0;
    pending_outgassing_reinjection_particles_ = 0.0;
    dynamic_neutral_density_boost_m3_ = 0.0;
    collision_emission_feedback_multiplier_ = 1.0;
    collision_channel_feedback_multiplier_ = 1.0;

    piccore_boundary_coupling_initialized_ = false;
    if (config_.alignment_mode == ArcPicAlignmentMode::ArcPicAligned)
    {
        piccore_boundary_coupling_initialized_ = initializePiccoreBoundaryCoupling();
        if (!piccore_boundary_coupling_initialized_)
        {
            std::cerr << "VacuumArc: PICcore boundary coupling initialization failed, continuing "
                         "without PIC boundary feedback."
                      << std::endl;
        }
    }

    resetPicBoundaryFeedback();
    return true;
}

bool SurfaceDischargeArcAlgorithm::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }

    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        return false;
    }

    const double min_substep_s = std::max(1.0e-15, config_.min_internal_substep_s);
    const double max_substep_s = std::max(min_substep_s, config_.max_internal_substep_s);
    const bool adaptive_timestep_enabled = config_.enable_adaptive_internal_timestep;
    const double adaptive_shrink_factor =
        std::clamp(config_.adaptive_substep_shrink_factor, 0.1, 0.95);
    const double adaptive_recovery_factor =
        std::max(1.0, config_.adaptive_substep_recovery_factor);
    const double adaptive_min_scale =
        std::clamp(config_.adaptive_substep_min_scale, 0.01, 1.0);
    const std::size_t max_stability_rollbacks =
        std::max<std::size_t>(1, config_.max_stability_rollbacks);

    status_.stability_substeps_used = 0;
    status_.stability_rollbacks = 0;
    status_.stability_effective_substep_s = 0.0;
    status_.stability_adaptive_scale = 1.0;
    status_.stability_adaptive_reductions = 0;
    status_.stability_anomaly_isolation_triggered = false;
    status_.stability_anomaly_isolation_events = 0;

    double remaining_time_s = dt;
    double adaptive_scale = 1.0;
    const std::size_t max_internal_substeps = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(dt / min_substep_s)) + 8);
    std::size_t guard_substeps = 0;

    while (remaining_time_s > 1.0e-18)
    {
        ++guard_substeps;
        if (guard_substeps > max_internal_substeps)
        {
            return false;
        }

        const double candidate_dt = std::min(max_substep_s, remaining_time_s);
        double local_dt = candidate_dt;
        if (adaptive_timestep_enabled)
        {
            local_dt = std::max(min_substep_s, candidate_dt * adaptive_scale);
        }
        local_dt = std::clamp(local_dt, min_substep_s, candidate_dt);

        if (!(local_dt > 0.0) || !std::isfinite(local_dt))
        {
            return false;
        }

        const double previous_surface_potential_v = config_.surface_potential_v;
        const double previous_total_current_a = status_.total_discharge_current_a;
        const ArcChannelState previous_channel_state = channel_state_;

        const double effective_field = computeEffectiveField();
        const double avalanche_gain =
            avalanche_model_.computeGain(effective_field, config_.gap_distance_m);

        const bool discharge_active = updateDischargeActivation(effective_field, avalanche_gain);

        if (discharge_active)
        {
            const double emitted_current_density = computeEmissionCurrentDensity(effective_field);
            updateChannelState(effective_field, emitted_current_density, local_dt);
            updateThermalState(emitted_current_density, local_dt);

            if (config_.alignment_mode == ArcPicAlignmentMode::ArcPicAligned)
            {
                advancePiccoreBoundaryCoupling(emitted_current_density, local_dt);
            }
            else
            {
                resetPicBoundaryFeedback();
            }

            status_.total_discharge_current_a += status_.pic_boundary_net_current_a;
        }
        else
        {
            resetPicBoundaryFeedback();
        }

        runCollisionStage(local_dt);
        updateSurfaceCircuitState(local_dt);
        updateResidualMonitoring(effective_field, local_dt);

        const bool anomalous_state = hasAnomalousState();
        if (anomalous_state)
        {
            if (!config_.enable_stability_rollback)
            {
                return false;
            }

            config_.surface_potential_v = previous_surface_potential_v;
            status_.surface_potential_v = previous_surface_potential_v;
            status_.total_discharge_current_a = previous_total_current_a;
            channel_state_ = previous_channel_state;
            status_.channel_conductivity_s_per_m = channel_state_.conductivity_s_per_m;
            ++status_.stability_rollbacks;
            status_.stability_anomaly_isolation_triggered = true;
            ++status_.stability_anomaly_isolation_events;

            if (status_.stability_rollbacks > max_stability_rollbacks &&
                !adaptive_timestep_enabled)
            {
                return false;
            }
        }

        if (adaptive_timestep_enabled)
        {
            const bool shrink_requested = anomalous_state || status_.residual_alarm_active;
            double next_scale = adaptive_scale;
            if (shrink_requested)
            {
                next_scale = std::max(adaptive_min_scale, adaptive_scale * adaptive_shrink_factor);
                if (next_scale + 1.0e-15 < adaptive_scale)
                {
                    ++status_.stability_adaptive_reductions;
                }
            }
            else
            {
                next_scale = std::min(1.0, adaptive_scale * adaptive_recovery_factor);
            }
            adaptive_scale = next_scale;
        }

        status_.stability_effective_substep_s = local_dt;
        status_.stability_adaptive_scale = adaptive_scale;
        status_.current_time_s += local_dt;
        remaining_time_s = std::max(0.0, remaining_time_s - local_dt);
        if (remaining_time_s <= 1.0e-15)
        {
            remaining_time_s = 0.0;
        }
        ++status_.stability_substeps_used;
        appendHistory();
    }

    return true;
}

void SurfaceDischargeArcAlgorithm::reset()
{
    *this = SurfaceDischargeArcAlgorithm{};
}

bool SurfaceDischargeArcAlgorithm::exportResults(const std::filesystem::path& csv_path) const
{
    const auto max_abs_from_history = [](const std::vector<double>& values) {
        double maximum = 0.0;
        for (const double value : values)
        {
            maximum = std::max(maximum, std::abs(value));
        }
        return maximum;
    };
    const auto sum_history = [](const std::vector<double>& values) {
        double sum = 0.0;
        for (const double value : values)
        {
            sum += value;
        }
        return sum;
    };
    const auto any_triggered = [](const std::vector<double>& values) {
        for (const double value : values)
        {
            if (value >= 0.5)
            {
                return true;
            }
        }
        return false;
    };
    const auto ignition_delay_s = [&]() {
        const std::size_t count = std::min(history_time_.size(), history_current_.size());
        for (std::size_t i = 0; i < count; ++i)
        {
            if (std::abs(history_current_[i]) > 1.0e-12)
            {
                return history_time_[i];
            }
        }
        return history_time_.empty() ? 0.0 : history_time_.front();
    }();
    const double total_ionization_events = sum_history(history_collision_ionization_events_);
    const double total_excitation_events = sum_history(history_collision_excitation_events_);
    const double total_charge_exchange_events =
        sum_history(history_collision_charge_exchange_events_);
    const double total_reaction_events =
        std::max(1.0e-18, total_ionization_events + total_excitation_events +
                              total_charge_exchange_events);
    std::filesystem::path benchmark_metrics_path = csv_path;
    benchmark_metrics_path.replace_extension(".benchmark_metrics.json");
    std::filesystem::path simulation_artifact_path = csv_path;
    simulation_artifact_path.replace_extension(".simulation_artifact.json");

    Output::ColumnarDataSet data_set;
    data_set.axis_name = "time_s";
    data_set.axis_values = history_time_;
    data_set.scalar_series["discharge_current_a"] = history_current_;
    data_set.scalar_series["current_density_a_per_m2"] = history_current_density_;
    data_set.scalar_series["cathode_temperature_k"] = history_cathode_temperature_;
    data_set.scalar_series["anode_temperature_k"] = history_anode_temperature_;
    data_set.scalar_series["channel_conductivity_s_per_m"] = history_conductivity_;
    data_set.scalar_series["pic_boundary_net_current_a"] = history_pic_boundary_net_current_;
    data_set.scalar_series["pic_surface_charge_delta_c"] = history_pic_surface_charge_delta_;
    data_set.scalar_series["surface_potential_v"] = history_surface_potential_;
    data_set.scalar_series["surface_charge_density_c_per_m2"] = history_surface_charge_density_;
    data_set.scalar_series["surface_circuit_drive_current_a"] = history_surface_circuit_drive_current_;
    data_set.scalar_series["surface_circuit_leak_current_a"] = history_surface_circuit_leak_current_;
    data_set.scalar_series["surface_load_branch_current_a"] = history_surface_load_branch_current_;
    data_set.scalar_series["surface_relative_potential_v"] = history_surface_relative_potential_;
    data_set.scalar_series["surface_circuit_capacitor_energy_j"] =
        history_surface_circuit_capacitor_energy_;
    data_set.scalar_series["surface_circuit_drive_power_w"] =
        history_surface_circuit_drive_power_;
    data_set.scalar_series["surface_circuit_leak_power_w"] = history_surface_circuit_leak_power_;
    data_set.scalar_series["surface_load_branch_power_w"] = history_surface_load_branch_power_;
    data_set.scalar_series["surface_load_resistive_power_w"] =
        history_surface_load_resistive_power_;
    data_set.scalar_series["surface_load_inductive_energy_j"] =
        history_surface_load_inductive_energy_;
    data_set.scalar_series["surface_circuit_power_balance_error_w"] =
        history_surface_circuit_power_balance_error_;
    data_set.scalar_series["surface_circuit_drive_energy_j"] =
        history_surface_circuit_drive_energy_;
    data_set.scalar_series["surface_circuit_leak_dissipated_energy_j"] =
        history_surface_circuit_leak_dissipated_energy_;
    data_set.scalar_series["surface_load_dissipated_energy_j"] =
        history_surface_load_dissipated_energy_;
    data_set.scalar_series["surface_circuit_energy_balance_error_j"] =
        history_surface_circuit_energy_balance_error_;
    data_set.scalar_series["charge_conservation_error_c"] = history_charge_conservation_error_;
    data_set.scalar_series["collision_events_step"] = history_collision_events_;
    data_set.scalar_series["collision_particles_created_step"] = history_collision_particles_created_;
    data_set.scalar_series["collision_ionization_events_step"] = history_collision_ionization_events_;
    data_set.scalar_series["collision_excitation_events_step"] = history_collision_excitation_events_;
    data_set.scalar_series["collision_charge_exchange_events_step"] =
        history_collision_charge_exchange_events_;
    data_set.scalar_series["collision_ionization_fraction_step"] =
        history_collision_ionization_fraction_;
    data_set.scalar_series["collision_excitation_fraction_step"] =
        history_collision_excitation_fraction_;
    data_set.scalar_series["collision_charge_exchange_fraction_step"] =
        history_collision_charge_exchange_fraction_;
    data_set.scalar_series["collision_reaction_weighted_emission_feedback_step"] =
        history_collision_emission_feedback_multiplier_;
    data_set.scalar_series["collision_reaction_weighted_channel_feedback_step"] =
        history_collision_channel_feedback_multiplier_;
    data_set.scalar_series["collision_event_rate_per_s"] = history_collision_event_rate_;
    data_set.scalar_series["collision_stage_fallback_triggered"] = history_collision_fallback_triggered_;
    data_set.scalar_series["collision_effective_neutral_density_m3"] = history_collision_neutral_density_;
    data_set.scalar_series["collision_effective_ion_density_m3"] = history_collision_ion_density_;
    data_set.scalar_series["collision_effective_electron_density_m3"] = history_collision_electron_density_;
    data_set.scalar_series["neutral_outgassing_density_boost_m3"] = history_neutral_outgassing_boost_;
    data_set.scalar_series["neutral_outgassing_reinjected_particles_step"] =
        history_neutral_outgassing_reinjected_particles_;
    data_set.scalar_series["neutral_outgassing_reinjection_rate_per_s"] =
        history_neutral_outgassing_reinjection_rate_;
    data_set.scalar_series["field_residual_v_per_m"] = history_field_residual_;
    data_set.scalar_series["particle_residual_a"] = history_particle_residual_;
    data_set.scalar_series["circuit_residual_c"] = history_circuit_residual_;
    data_set.scalar_series["field_residual_filtered_v_per_m"] = history_field_residual_filtered_;
    data_set.scalar_series["particle_residual_filtered_a"] = history_particle_residual_filtered_;
    data_set.scalar_series["circuit_residual_filtered_c"] = history_circuit_residual_filtered_;
    data_set.scalar_series["residual_alarm_active"] = history_residual_alarm_active_;
    data_set.scalar_series["residual_alarm_counter"] = history_residual_alarm_counter_;
    data_set.scalar_series["residual_alarm_clear_counter"] = history_residual_alarm_clear_counter_;
    data_set.scalar_series["stability_substeps_used"] = history_stability_substeps_;
    data_set.scalar_series["stability_rollbacks"] = history_stability_rollbacks_;
    data_set.scalar_series["stability_effective_substep_s"] =
        history_stability_effective_substep_;
    data_set.scalar_series["stability_adaptive_scale"] = history_stability_adaptive_scale_;
    data_set.scalar_series["stability_adaptive_reductions"] =
        history_stability_adaptive_reductions_;
    data_set.scalar_series["stability_anomaly_isolation_triggered"] =
        history_stability_anomaly_isolation_triggered_;
    data_set.scalar_series["stability_anomaly_isolation_events"] =
        history_stability_anomaly_isolation_events_;
    data_set.metadata["module"] = "Vacuum Arc";
    data_set.metadata["emission_strategy"] = emission_model_.getStrategyName();
    data_set.metadata["alignment_mode"] = alignmentModeName(config_.alignment_mode);
    data_set.metadata["pipeline_contract_id"] = pipelineContractId();
    data_set.metadata["pipeline_stage_order"] = pipelineStageOrder();
    data_set.metadata["benchmark_contract_family"] = "arcpic_core_metrics_v1";
    data_set.metadata["benchmark_focus_metrics"] =
        "ignition_delay,peak_current,channel_conductivity,surface_potential,reaction_split";
    data_set.metadata["benchmark_metrics_contract_id"] =
        "vacuum-arc-benchmark-metrics-v1";
    data_set.metadata["benchmark_metrics_artifact_path"] =
        benchmark_metrics_path.filename().string();
    data_set.metadata["collision_diagnostic_contract_id"] =
        "vacuum-arc-collision-emission-channel-v1";
    data_set.metadata["piccore_boundary_coupling"] =
        piccore_boundary_coupling_initialized_ ? "enabled" : "disabled";
    data_set.metadata["surface_capacitance_f"] = std::to_string(config_.surface_capacitance_f);
    data_set.metadata["surface_leakage_conductance_s"] =
        std::to_string(config_.surface_leakage_conductance_s);
    data_set.metadata["surface_reference_potential_v"] =
        std::to_string(config_.surface_reference_potential_v);
    data_set.metadata["surface_load_model"] = surfaceLoadModelName(config_.surface_load_model);
    data_set.metadata["surface_load_resistance_ohm"] =
        std::to_string(config_.surface_load_resistance_ohm);
    data_set.metadata["surface_load_inductance_h"] =
        std::to_string(config_.surface_load_inductance_h);
    data_set.metadata["integrator_mode"] =
        integratorModeName(config_.integrator_parameters.mode);
    data_set.metadata["channel_max_current_density_a_per_m2"] =
        std::to_string(config_.channel_parameters.max_current_density_a_per_m2);
    data_set.metadata["integrator_max_current_density_a_per_m2"] =
        std::to_string(config_.integrator_parameters.max_current_density_a_per_m2);
    data_set.metadata["collision_stage_enabled"] = config_.enable_pic_mcc_collisions ? "true" : "false";
    data_set.metadata["collision_reconfigure_interval_steps"] =
        std::to_string(config_.collision_reconfigure_interval_steps);
    data_set.metadata["collision_neutral_density_floor_m3"] =
        std::to_string(config_.collision_neutral_density_floor_m3);
    data_set.metadata["collision_reaction_breakdown_enabled"] =
        config_.enable_collision_reaction_breakdown ? "true" : "false";
    data_set.metadata["collision_cross_section_set_id"] = config_.collision_cross_section_set_id;
    data_set.metadata["collision_anomaly_event_rate_limit_per_s"] =
        std::to_string(config_.collision_anomaly_event_rate_limit_per_s);
    data_set.metadata["collision_fallback_to_background_mcc_on_error"] =
        config_.collision_fallback_to_background_mcc_on_error ? "true" : "false";
    data_set.metadata["collision_emission_feedback_gain"] =
        std::to_string(config_.collision_emission_feedback_gain);
    data_set.metadata["collision_channel_feedback_gain"] =
        std::to_string(config_.collision_channel_feedback_gain);
    data_set.metadata["collision_ionization_emission_feedback_gain"] =
        std::to_string(config_.collision_ionization_emission_feedback_gain);
    data_set.metadata["collision_excitation_emission_feedback_gain"] =
        std::to_string(config_.collision_excitation_emission_feedback_gain);
    data_set.metadata["collision_charge_exchange_emission_feedback_gain"] =
        std::to_string(config_.collision_charge_exchange_emission_feedback_gain);
    data_set.metadata["collision_ionization_channel_feedback_gain"] =
        std::to_string(config_.collision_ionization_channel_feedback_gain);
    data_set.metadata["collision_excitation_channel_feedback_gain"] =
        std::to_string(config_.collision_excitation_channel_feedback_gain);
    data_set.metadata["collision_charge_exchange_channel_feedback_gain"] =
        std::to_string(config_.collision_charge_exchange_channel_feedback_gain);
    data_set.metadata["secondary_electron_emission_enabled"] =
        config_.enable_secondary_electron_emission ? "true" : "false";
    data_set.metadata["secondary_electron_yield_per_ion"] =
        std::to_string(config_.secondary_electron_yield_per_ion);
    data_set.metadata["neutral_outgassing_feedback_enabled"] =
        config_.enable_neutral_outgassing_feedback ? "true" : "false";
    data_set.metadata["neutral_outgassing_gain_m3_per_a"] =
        std::to_string(config_.neutral_outgassing_gain_m3_per_a);
    data_set.metadata["neutral_outgassing_reinjection_enabled"] =
        config_.enable_neutral_outgassing_reinjection ? "true" : "false";
    data_set.metadata["neutral_outgassing_reinjection_gain"] =
        std::to_string(config_.neutral_outgassing_reinjection_gain);
    data_set.metadata["neutral_outgassing_reinjection_macro_weight"] =
        std::to_string(config_.neutral_outgassing_reinjection_macro_weight);
    data_set.metadata["neutral_outgassing_reinjection_max_particles_per_step"] =
        std::to_string(config_.neutral_outgassing_reinjection_max_particles_per_step);
    data_set.metadata["neutral_outgassing_reinjection_speed_m_per_s"] =
        std::to_string(config_.neutral_outgassing_reinjection_speed_m_per_s);
    data_set.metadata["neutral_outgassing_reinjection_mass_amu"] =
        std::to_string(config_.neutral_outgassing_reinjection_mass_amu);
    data_set.metadata["residual_monitoring_enabled"] =
        config_.enable_residual_monitoring ? "true" : "false";
    data_set.metadata["residual_field_tolerance_v_per_m"] =
        std::to_string(config_.residual_field_tolerance_v_per_m);
    data_set.metadata["residual_particle_current_tolerance_a"] =
        std::to_string(config_.residual_particle_current_tolerance_a);
    data_set.metadata["residual_circuit_charge_tolerance_c"] =
        std::to_string(config_.residual_circuit_charge_tolerance_c);
    data_set.metadata["residual_alarm_consecutive_steps"] =
        std::to_string(config_.residual_alarm_consecutive_steps);
    data_set.metadata["residual_alarm_clear_consecutive_steps"] =
        std::to_string(config_.residual_alarm_clear_consecutive_steps);
    data_set.metadata["residual_ema_filter_enabled"] =
        config_.enable_residual_ema_filter ? "true" : "false";
    data_set.metadata["residual_ema_alpha"] =
        std::to_string(config_.residual_ema_alpha);
    data_set.metadata["stability_min_internal_substep_s"] =
        std::to_string(config_.min_internal_substep_s);
    data_set.metadata["stability_max_internal_substep_s"] =
        std::to_string(config_.max_internal_substep_s);
    data_set.metadata["stability_adaptive_timestep_enabled"] =
        config_.enable_adaptive_internal_timestep ? "true" : "false";
    data_set.metadata["stability_adaptive_substep_shrink_factor"] =
        std::to_string(config_.adaptive_substep_shrink_factor);
    data_set.metadata["stability_adaptive_substep_recovery_factor"] =
        std::to_string(config_.adaptive_substep_recovery_factor);
    data_set.metadata["stability_adaptive_substep_min_scale"] =
        std::to_string(config_.adaptive_substep_min_scale);
    data_set.metadata["stability_rollback_enabled"] =
        config_.enable_stability_rollback ? "true" : "false";
    data_set.metadata["stability_max_rollbacks"] =
        std::to_string(config_.max_stability_rollbacks);
    Coupling::Contracts::appendSolverConfigMetadata(data_set.metadata, config_.solver_config, "solver_");
    data_set.metadata["seed"] = std::to_string(config_.seed);
    data_set.metadata["sampling_policy"] = config_.sampling_policy;
    data_set.metadata["simulation_artifact_contract_id"] = "simulation-artifact-v1";
    data_set.metadata["simulation_artifact_path"] = simulation_artifact_path.filename().string();
    if (!static_cast<bool>(exporter_.exportDataSet(csv_path, data_set)))
    {
        return false;
    }

    std::ofstream benchmark_output(benchmark_metrics_path, std::ios::out | std::ios::trunc);
    if (!benchmark_output.is_open())
    {
        return false;
    }

    benchmark_output << std::fixed << std::setprecision(12);
    benchmark_output << "{\n";
    benchmark_output << "  \"schema_version\": \"scdat.vacuum_arc.benchmark_metrics.v1\",\n";
    benchmark_output << "  \"contract_id\": \"vacuum-arc-benchmark-metrics-v1\",\n";
    benchmark_output << "  \"module\": \"Vacuum Arc\",\n";
    benchmark_output << "  \"alignment_mode\": \"" << alignmentModeName(config_.alignment_mode)
                     << "\",\n";
    benchmark_output << "  \"ignition_delay_s\": " << ignition_delay_s << ",\n";
    benchmark_output << "  \"peak_current_a\": " << max_abs_from_history(history_current_) << ",\n";
    benchmark_output << "  \"peak_current_density_a_per_m2\": "
                     << max_abs_from_history(history_current_density_) << ",\n";
    benchmark_output << "  \"max_channel_conductivity_s_per_m\": "
                     << max_abs_from_history(history_conductivity_) << ",\n";
    benchmark_output << "  \"final_surface_potential_v\": "
                     << (history_surface_potential_.empty() ? 0.0 : history_surface_potential_.back())
                     << ",\n";
    benchmark_output << "  \"max_abs_charge_conservation_error_c\": "
                     << max_abs_from_history(history_charge_conservation_error_) << ",\n";
    benchmark_output << "  \"collision_fallback_triggered\": "
                     << (any_triggered(history_collision_fallback_triggered_) ? "true" : "false")
                     << ",\n";
    benchmark_output << "  \"reaction_split\": {\n";
    benchmark_output << "    \"ionization_fraction\": "
                     << (total_ionization_events / total_reaction_events) << ",\n";
    benchmark_output << "    \"excitation_fraction\": "
                     << (total_excitation_events / total_reaction_events) << ",\n";
    benchmark_output << "    \"charge_exchange_fraction\": "
                     << (total_charge_exchange_events / total_reaction_events) << "\n";
    benchmark_output << "  }\n";
    benchmark_output << "}\n";
    if (!static_cast<bool>(benchmark_output))
    {
        return false;
    }

    Coupling::Contracts::SimulationArtifact artifact;
    artifact.module = "Vacuum Arc";
    artifact.case_id = alignmentModeName(config_.alignment_mode);
    artifact.reference_family = "arcpic";
    artifact.seed = config_.seed;
    artifact.sampling_policy = config_.sampling_policy;
    artifact.particle_metrics["collision_events_total"] = total_reaction_events;
    artifact.particle_metrics["collision_ionization_events_total"] = total_ionization_events;
    artifact.particle_metrics["collision_excitation_events_total"] = total_excitation_events;
    artifact.particle_metrics["collision_charge_exchange_events_total"] = total_charge_exchange_events;
    artifact.surface_metrics["ignition_delay_s"] = ignition_delay_s;
    artifact.surface_metrics["peak_current_a"] = max_abs_from_history(history_current_);
    artifact.surface_metrics["peak_current_density_a_per_m2"] =
        max_abs_from_history(history_current_density_);
    artifact.surface_metrics["max_channel_conductivity_s_per_m"] =
        max_abs_from_history(history_conductivity_);
    artifact.surface_metrics["final_surface_potential_v"] =
        history_surface_potential_.empty() ? 0.0 : history_surface_potential_.back();
    artifact.surface_metrics["max_abs_charge_conservation_error_c"] =
        max_abs_from_history(history_charge_conservation_error_);
    artifact.metadata["pipeline_contract_id"] = pipelineContractId();
    artifact.metadata["benchmark_metrics_contract_id"] = "vacuum-arc-benchmark-metrics-v1";
    artifact.metadata["collision_set"] = config_.collision_cross_section_set_id;
    std::string artifact_error;
    if (!Coupling::Contracts::writeSimulationArtifactJson(
            simulation_artifact_path, artifact, &artifact_error))
    {
        return false;
    }
    return true;
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
