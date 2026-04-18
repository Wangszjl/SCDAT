#include "../include/SurfaceFieldVolumeBridge.h"

#include "../include/BoundaryCondition.h"
#include "../include/PoissonSolver.h"
#include "../include/VolDistrib.h"
#include "../include/VolInteract.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace SCDAT
{
namespace FieldSolver
{
namespace
{

double finiteOr(double value, double fallback)
{
    return std::isfinite(value) ? value : fallback;
}

bool isElectronLike(Particle::ParticleType type)
{
    using Particle::ParticleType;
    switch (type)
    {
    case ParticleType::ELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::BACKSCATTERED_ELECTRON:
    case ParticleType::AUGER_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
        return true;
    default:
        return false;
    }
}

bool isIonLike(Particle::ParticleType type)
{
    using Particle::ParticleType;
    switch (type)
    {
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
    case ParticleType::ALPHA:
    case ParticleType::BETA_PARTICLE:
        return true;
    default:
        return false;
    }
}

template <typename Family>
std::string buildFamilySignature(const std::vector<Family>& families,
                                 const char* (*name_of)(Family))
{
    std::ostringstream stream;
    bool first = true;
    for (const auto family : families)
    {
        if (!first)
        {
            stream << "+";
        }
        stream << name_of(family);
        first = false;
    }
    return stream.str();
}

template <typename Family>
std::vector<Family> normalizeFamilyList(const std::vector<Family>& families)
{
    std::vector<Family> normalized;
    normalized.reserve(families.size());
    for (const auto family : families)
    {
        if (std::find(normalized.begin(), normalized.end(), family) == normalized.end())
        {
            normalized.push_back(family);
        }
    }
    return normalized;
}

} // namespace

const char* nativeVolumeParityRouteName(NativeVolumeParityRoute route)
{
    switch (route)
    {
    case NativeVolumeParityRoute::NativeMinimal:
        return "native_minimal";
    case NativeVolumeParityRoute::ExternalBridge:
        return "external_bridge";
    case NativeVolumeParityRoute::HybridBlend:
        return "hybrid_blend";
    }
    return "native_minimal";
}

const char* nativeVolumeBoundaryConditionFamilyName(NativeVolumeBoundaryConditionFamily family)
{
    switch (family)
    {
    case NativeVolumeBoundaryConditionFamily::VoltageDependentMBC:
        return "VoltageDependentMBC";
    case NativeVolumeBoundaryConditionFamily::MixedDirichletFourierPoissonBC:
        return "MixedDirichletFourierPoissonBC";
    case NativeVolumeBoundaryConditionFamily::FourierPoissonBC:
        return "FourierPoissonBC";
    case NativeVolumeBoundaryConditionFamily::SurfDistribMatterBC:
        return "SurfDistribMatterBC";
    case NativeVolumeBoundaryConditionFamily::OneSurfDistribTestableMatterBC:
        return "OneSurfDistribTestableMatterBC";
    case NativeVolumeBoundaryConditionFamily::CapacitiveVoltageGenerator:
        return "CapacitiveVoltageGenerator";
    }
    return "FourierPoissonBC";
}

const char* nativeVolumeFieldFamilyName(NativeVolumeFieldFamily family)
{
    switch (family)
    {
    case NativeVolumeFieldFamily::UniformBField:
        return "UniformBField";
    case NativeVolumeFieldFamily::SolenoidBField:
        return "SolenoidBField";
    case NativeVolumeFieldFamily::DipolarBField:
        return "DipolarBField";
    case NativeVolumeFieldFamily::MultipleEField:
        return "MultipleEField";
    case NativeVolumeFieldFamily::EMField:
        return "EMField";
    }
    return "EMField";
}

const char* nativeVolumeDistributionFamilyName(NativeVolumeDistributionFamily family)
{
    switch (family)
    {
    case NativeVolumeDistributionFamily::PICVolDistrib:
        return "PICVolDistrib";
    case NativeVolumeDistributionFamily::PICVolDistribNoAcc:
        return "PICVolDistribNoAcc";
    case NativeVolumeDistributionFamily::PICVolDistribUpdatable:
        return "PICVolDistribUpdatable";
    case NativeVolumeDistributionFamily::SmartPICVolDistrib:
        return "SmartPICVolDistrib";
    case NativeVolumeDistributionFamily::CompositeVolDistrib:
        return "CompositeVolDistrib";
    case NativeVolumeDistributionFamily::BackTrackingVolDistrib:
        return "BackTrackingVolDistrib";
    case NativeVolumeDistributionFamily::BacktrackingPICCompositeVolDistrib:
        return "BacktrackingPICCompositeVolDistrib";
    case NativeVolumeDistributionFamily::BacktrackingBoltzmannCompositeVolDistrib:
        return "BacktrackingBoltzmannCompositeVolDistrib";
    case NativeVolumeDistributionFamily::VolDistrib:
        return "VolDistrib";
    case NativeVolumeDistributionFamily::NonPICVolDistrib:
        return "NonPICVolDistrib";
    case NativeVolumeDistributionFamily::AnalyticVolDistrib:
        return "AnalyticVolDistrib";
    case NativeVolumeDistributionFamily::MultipleVolDistrib:
        return "MultipleVolDistrib";
    case NativeVolumeDistributionFamily::LocalMaxwellVolDistrib:
        return "LocalMaxwellVolDistrib";
    case NativeVolumeDistributionFamily::LocalMaxwellBoltzmannVolDistrib:
        return "LocalMaxwellBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::GlobalMaxwellBoltzmannVolDistrib:
        return "GlobalMaxwellBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::PICBoltzmannVolDistrib:
        return "PICBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::SteadyMaxwellBoltzmannVolDistrib:
        return "SteadyMaxwellBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::UnlimitedGlobalMaxwellBoltzmannVolDistrib:
        return "UnlimitedGlobalMaxwellBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::SurfaceLimitedGlobalMaxwellBoltzmannVolDistrib:
        return "SurfaceLimitedGlobalMaxwellBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::TrunckatedGlobalMaxwellBoltzmannVolDistrib:
        return "TrunckatedGlobalMaxwellBoltzmannVolDistrib";
    case NativeVolumeDistributionFamily::ImplicitableVolDistrib:
        return "ImplicitableVolDistrib";
    case NativeVolumeDistributionFamily::Updatable:
        return "Updatable";
    }
    return "PICVolDistrib";
}

const char* nativeVolumeInteractionFamilyName(NativeVolumeInteractionFamily family)
{
    switch (family)
    {
    case NativeVolumeInteractionFamily::MCCInteractor:
        return "MCCInteractor";
    case NativeVolumeInteractionFamily::CEXInteractor:
        return "CEXInteractor";
    case NativeVolumeInteractionFamily::PhotoIonization:
        return "PhotoIonization";
    case NativeVolumeInteractionFamily::TrajectoryInteractionFromField:
        return "TrajectoryInteractionFromField";
    case NativeVolumeInteractionFamily::SpinningSpacecraftTrajectory:
        return "SpinningSpacecraftTrajectory";
    }
    return "MCCInteractor";
}

std::vector<NativeVolumeBoundaryConditionFamily> resolveNativeVolumeBoundaryConditionFamilies(
    NativeVolumeParityRoute route)
{
    if (route == NativeVolumeParityRoute::NativeMinimal)
    {
        return {
            NativeVolumeBoundaryConditionFamily::FourierPoissonBC,
            NativeVolumeBoundaryConditionFamily::SurfDistribMatterBC,
        };
    }
    return {
        NativeVolumeBoundaryConditionFamily::VoltageDependentMBC,
        NativeVolumeBoundaryConditionFamily::MixedDirichletFourierPoissonBC,
        NativeVolumeBoundaryConditionFamily::FourierPoissonBC,
        NativeVolumeBoundaryConditionFamily::SurfDistribMatterBC,
        NativeVolumeBoundaryConditionFamily::OneSurfDistribTestableMatterBC,
        NativeVolumeBoundaryConditionFamily::CapacitiveVoltageGenerator,
    };
}

std::vector<NativeVolumeBoundaryConditionFamily> resolveActiveNativeVolumeBoundaryConditionFamilies(
    NativeVolumeParityRoute route)
{
    return resolveNativeVolumeBoundaryConditionFamilies(route);
}

std::vector<NativeVolumeBoundaryConditionFamily> normalizeNativeVolumeBoundaryConditionFamilies(
    const std::vector<NativeVolumeBoundaryConditionFamily>& families)
{
    return normalizeFamilyList(families);
}

std::vector<NativeVolumeFieldFamily> resolveNativeVolumeFieldFamilies(
    NativeVolumeParityRoute route)
{
    if (route == NativeVolumeParityRoute::NativeMinimal)
    {
        return {
            NativeVolumeFieldFamily::UniformBField,
            NativeVolumeFieldFamily::MultipleEField,
            NativeVolumeFieldFamily::EMField,
        };
    }
    return {
        NativeVolumeFieldFamily::UniformBField,
        NativeVolumeFieldFamily::SolenoidBField,
        NativeVolumeFieldFamily::DipolarBField,
        NativeVolumeFieldFamily::MultipleEField,
        NativeVolumeFieldFamily::EMField,
    };
}

std::vector<NativeVolumeFieldFamily> resolveActiveNativeVolumeFieldFamilies(
    NativeVolumeParityRoute route)
{
    return resolveNativeVolumeFieldFamilies(route);
}

std::vector<NativeVolumeFieldFamily> normalizeNativeVolumeFieldFamilies(
    const std::vector<NativeVolumeFieldFamily>& families)
{
    return normalizeFamilyList(families);
}

std::vector<NativeVolumeDistributionFamily> resolveNativeVolumeDistributionFamilies(
    NativeVolumeParityRoute route)
{
    if (route == NativeVolumeParityRoute::NativeMinimal)
    {
        return {
            NativeVolumeDistributionFamily::PICVolDistrib,
            NativeVolumeDistributionFamily::PICVolDistribUpdatable,
            NativeVolumeDistributionFamily::CompositeVolDistrib,
        };
    }
    return {
        NativeVolumeDistributionFamily::PICVolDistrib,
        NativeVolumeDistributionFamily::PICVolDistribNoAcc,
        NativeVolumeDistributionFamily::PICVolDistribUpdatable,
        NativeVolumeDistributionFamily::SmartPICVolDistrib,
        NativeVolumeDistributionFamily::CompositeVolDistrib,
        NativeVolumeDistributionFamily::BackTrackingVolDistrib,
        NativeVolumeDistributionFamily::BacktrackingPICCompositeVolDistrib,
        NativeVolumeDistributionFamily::BacktrackingBoltzmannCompositeVolDistrib,
    };
}

std::vector<NativeVolumeDistributionFamily> normalizeNativeVolumeDistributionFamilies(
    const std::vector<NativeVolumeDistributionFamily>& families)
{
    return normalizeFamilyList(families);
}

std::vector<NativeVolumeInteractionFamily> resolveNativeVolumeInteractionFamilies(
    NativeVolumeParityRoute route)
{
    if (route == NativeVolumeParityRoute::NativeMinimal)
    {
        return {
            NativeVolumeInteractionFamily::MCCInteractor,
            NativeVolumeInteractionFamily::TrajectoryInteractionFromField,
        };
    }
    return {
        NativeVolumeInteractionFamily::MCCInteractor,
        NativeVolumeInteractionFamily::CEXInteractor,
        NativeVolumeInteractionFamily::PhotoIonization,
        NativeVolumeInteractionFamily::TrajectoryInteractionFromField,
        NativeVolumeInteractionFamily::SpinningSpacecraftTrajectory,
    };
}

std::vector<NativeVolumeInteractionFamily> resolveActiveNativeVolumeInteractionFamilies(
    NativeVolumeParityRoute route)
{
    return resolveNativeVolumeInteractionFamilies(route);
}

std::vector<NativeVolumeInteractionFamily> normalizeNativeVolumeInteractionFamilies(
    const std::vector<NativeVolumeInteractionFamily>& families)
{
    return normalizeFamilyList(families);
}

std::string nativeVolumeBoundaryConditionFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeBoundaryConditionFamilySignature(
        resolveNativeVolumeBoundaryConditionFamilies(route));
}

std::string nativeVolumeActiveBoundaryConditionFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeBoundaryConditionFamilySignature(
        resolveActiveNativeVolumeBoundaryConditionFamilies(route));
}

std::string nativeVolumeBoundaryConditionFamilySignature(
    const std::vector<NativeVolumeBoundaryConditionFamily>& families)
{
    return buildFamilySignature(normalizeNativeVolumeBoundaryConditionFamilies(families),
                                nativeVolumeBoundaryConditionFamilyName);
}

std::string nativeVolumeFieldFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeFieldFamilySignature(resolveNativeVolumeFieldFamilies(route));
}

std::string nativeVolumeActiveFieldFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeFieldFamilySignature(resolveActiveNativeVolumeFieldFamilies(route));
}

std::string nativeVolumeFieldFamilySignature(
    const std::vector<NativeVolumeFieldFamily>& families)
{
    return buildFamilySignature(normalizeNativeVolumeFieldFamilies(families),
                                nativeVolumeFieldFamilyName);
}

std::string nativeVolumeDistributionFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeDistributionFamilySignature(resolveNativeVolumeDistributionFamilies(route));
}

std::string nativeVolumeDistributionFamilySignature(
    const std::vector<NativeVolumeDistributionFamily>& families)
{
    return buildFamilySignature(normalizeNativeVolumeDistributionFamilies(families),
                                nativeVolumeDistributionFamilyName);
}

std::string nativeVolumeInteractionFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeInteractionFamilySignature(resolveNativeVolumeInteractionFamilies(route));
}

std::string nativeVolumeActiveInteractionFamilySignature(NativeVolumeParityRoute route)
{
    return nativeVolumeInteractionFamilySignature(
        resolveActiveNativeVolumeInteractionFamilies(route));
}

std::string nativeVolumeInteractionFamilySignature(
    const std::vector<NativeVolumeInteractionFamily>& families)
{
    return buildFamilySignature(normalizeNativeVolumeInteractionFamilies(families),
                                nativeVolumeInteractionFamilyName);
}

double blendFieldVolumeScalar(double base_value, double target_value, double blend_factor)
{
    const double base = finiteOr(base_value, 0.0);
    const double target = finiteOr(target_value, base);
    const double blend = std::isfinite(blend_factor)
                             ? std::clamp(blend_factor, 0.0, 1.0)
                             : 0.0;
    return (1.0 - blend) * base + blend * target;
}

double updateFieldVolumeMismatchMetric(double current_metric, double candidate_value,
                                       double baseline_value, double tolerance)
{
    double metric = finiteOr(current_metric, 0.0);
    if (!std::isfinite(candidate_value) || !std::isfinite(baseline_value))
    {
        return metric;
    }

    const double safe_tolerance = std::max(1.0e-12, std::abs(finiteOr(tolerance, 0.0)));
    const double normalized = std::abs(candidate_value - baseline_value) / safe_tolerance;
    if (!std::isfinite(normalized))
    {
        return metric;
    }
    return std::max(metric, normalized);
}

double computeFieldVolumeExternalBlendFactor(double mismatch_metric,
                                             double linear_residual_norm,
                                             double max_delta_v,
                                             bool linear_converged,
                                             double potential_tolerance_v,
                                             double volume_linear_relaxation)
{
    const double safe_tolerance = std::max(1.0e-9, std::abs(finiteOr(potential_tolerance_v, 0.0)));
    const double safe_mismatch = std::max(0.0, finiteOr(mismatch_metric, 0.0));
    const double safe_residual = std::abs(finiteOr(linear_residual_norm, 0.0));
    const double safe_delta = std::abs(finiteOr(max_delta_v, 0.0));

    const double residual_gate = 1.0 / (1.0 + safe_residual / safe_tolerance);
    const double delta_gate = 1.0 / (1.0 + safe_delta / safe_tolerance);
    const double convergence_gate = linear_converged ? 1.0 : 0.4;
    const double internal_confidence = convergence_gate * 0.5 * (residual_gate + delta_gate);
    const double stability_gate = 1.0 / (1.0 + safe_mismatch);
    const double external_trust =
        std::clamp(stability_gate + 0.5 * (1.0 - internal_confidence), 0.05, 1.0);
    const double relaxed_weight =
        std::clamp(finiteOr(volume_linear_relaxation, 0.2), 0.2, 1.0) * external_trust;
    return std::clamp(relaxed_weight, 0.05, 1.0);
}

NativeVolumeMeshSummary summarizeNativeVolumeMesh(const Mesh::VolMesh& mesh)
{
    NativeVolumeMeshSummary summary;
    summary.node_count = mesh.getNodeCount();
    summary.element_count = mesh.getElementCount();
    summary.valid = mesh.validate();

    const auto bbox_min = mesh.getBoundingBoxMin();
    const auto bbox_max = mesh.getBoundingBoxMax();
    summary.bbox_span_x_m = std::max(0.0, bbox_max.x() - bbox_min.x());
    summary.bbox_span_y_m = std::max(0.0, bbox_max.y() - bbox_min.y());
    summary.bbox_span_z_m = std::max(0.0, bbox_max.z() - bbox_min.z());

    double total_volume = 0.0;
    for (std::size_t i = 0; i < mesh.getElementCount(); ++i)
    {
        total_volume += std::max(0.0, finiteOr(mesh.getElementVolume(i), 0.0));
    }
    summary.total_volume_m3 = total_volume;
    return summary;
}

NativeVolumeFieldSummary solveNativeVolumePotential(const Mesh::VolMeshPtr& mesh,
                                                    double lower_dirichlet_v,
                                                    double upper_dirichlet_v,
                                                    NativeVolumeParityRoute route)
{
    NativeVolumeFieldSummary summary;
    const auto boundary_families = resolveNativeVolumeBoundaryConditionFamilies(route);
    summary.boundary_condition_family_count = boundary_families.size();
    summary.boundary_condition_family_signature =
        nativeVolumeBoundaryConditionFamilySignature(route);
    const auto field_families = resolveNativeVolumeFieldFamilies(route);
    summary.field_family_count = field_families.size();
    summary.field_family_signature = nativeVolumeFieldFamilySignature(route);
    if (!mesh || mesh->getNodeCount() == 0)
    {
        return summary;
    }

    const double z_min = mesh->getBoundingBoxMin().z();
    const double z_max = mesh->getBoundingBoxMax().z();
    const double z_mid = 0.5 * (z_min + z_max);
    const bool use_midplane_constraint = route != NativeVolumeParityRoute::NativeMinimal;
    bool lower_boundary_registered = false;
    bool upper_boundary_registered = false;
    bool midplane_boundary_registered = false;

    PoissonSolver solver(mesh);
    for (const auto& node : mesh->getNodes())
    {
        if (!node)
        {
            continue;
        }

        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12)
        {
            node->setBoundaryType(Mesh::BoundaryType::DIRICHLET);
            solver.addBoundaryCondition(
                node->getId(),
                std::make_shared<BoundaryCondition>(BoundaryConditionType::DIRICHLET,
                                                    lower_dirichlet_v));
            if (!lower_boundary_registered)
            {
                lower_boundary_registered = true;
                ++summary.boundary_condition_count;
                ++summary.dirichlet_boundary_count;
            }
        }
        else if (std::abs(z - z_max) < 1.0e-12)
        {
            node->setBoundaryType(Mesh::BoundaryType::DIRICHLET);
            solver.addBoundaryCondition(
                node->getId(),
                std::make_shared<BoundaryCondition>(BoundaryConditionType::DIRICHLET,
                                                    upper_dirichlet_v));
            if (!upper_boundary_registered)
            {
                upper_boundary_registered = true;
                ++summary.boundary_condition_count;
                ++summary.dirichlet_boundary_count;
            }
        }
        else if (use_midplane_constraint && std::abs(z - z_mid) < 1.0e-12)
        {
            node->setBoundaryType(Mesh::BoundaryType::DIRICHLET);
            solver.addBoundaryCondition(
                node->getId(),
                std::make_shared<BoundaryCondition>(
                    BoundaryConditionType::DIRICHLET,
                    0.5 * (lower_dirichlet_v + upper_dirichlet_v)));
            if (!midplane_boundary_registered)
            {
                midplane_boundary_registered = true;
                ++summary.boundary_condition_count;
                ++summary.dirichlet_boundary_count;
            }
        }
        else
        {
            node->setBoundaryType(Mesh::BoundaryType::INTERIOR);
        }
    }

    summary.solved = solver.solve(Solver::SolverType::DIRECT);
    if (!summary.solved)
    {
        return summary;
    }

    summary.min_potential_v = finiteOr(solver.getMinPotential(), 0.0);
    summary.max_potential_v = finiteOr(solver.getMaxPotential(), 0.0);
    summary.total_energy_j = std::max(0.0, finiteOr(solver.getTotalEnergy(), 0.0));

    double max_field = 0.0;
    for (const auto& node : mesh->getNodes())
    {
        if (!node)
        {
            continue;
        }
        max_field = std::max(max_field, finiteOr(node->getElectricField().magnitude(), 0.0));
    }
    summary.max_field_v_per_m = max_field;
    return summary;
}

NativeVolumeDistributionSummary evaluateNativeVolumeDistributions(
    std::size_t node_count, const std::vector<Particle::ParticlePtr>& particles,
    NativeVolumeParityRoute route)
{
    NativeVolumeDistributionSummary summary;
    const auto distribution_families = resolveNativeVolumeDistributionFamilies(route);
    summary.family_count = distribution_families.size();
    summary.family_signature = nativeVolumeDistributionFamilySignature(route);
    if (node_count == 0)
    {
        return summary;
    }

    std::vector<Particle::ParticlePtr> electron_particles;
    std::vector<Particle::ParticlePtr> ion_particles;
    for (auto* particle : particles)
    {
        if (!particle || !particle->isActive())
        {
            continue;
        }
        ++summary.active_particle_count;
        if (isElectronLike(particle->getType()))
        {
            electron_particles.push_back(particle);
        }
        if (isIonLike(particle->getType()))
        {
            ion_particles.push_back(particle);
        }
    }

    summary.electron_particle_count = electron_particles.size();
    summary.ion_particle_count = ion_particles.size();
    summary.distribution_channel_count = 3;

    Field::ParticleDensityDistrib electron_density(node_count, Particle::ParticleType::ELECTRON);
    electron_density.updateFromParticles(electron_particles);
    Field::ChargeDensityDistrib electron_charge(node_count, Particle::ParticleType::ELECTRON);
    electron_charge.computeFromParticleDensity(electron_density);

    Field::ParticleDensityDistrib ion_density(node_count, Particle::ParticleType::ION);
    ion_density.updateFromParticles(ion_particles);
    Field::ChargeDensityDistrib ion_charge(node_count, Particle::ParticleType::ION);
    ion_charge.computeFromParticleDensity(ion_density);

    double electron_sum = 0.0;
    double ion_sum = 0.0;
    double net_charge_sum = 0.0;
    for (std::size_t i = 0; i < node_count; ++i)
    {
        electron_sum += electron_density.getDensity(i);
        ion_sum += ion_density.getDensity(i);
        net_charge_sum += electron_charge.getChargeDensity(i) + ion_charge.getChargeDensity(i);
    }

    const double denom = static_cast<double>(node_count);
    summary.average_electron_density_m3 = electron_sum / denom;
    summary.average_ion_density_m3 = ion_sum / denom;
    summary.average_net_charge_density_c_per_m3 = net_charge_sum / denom;
    return summary;
}

NativeVolumeInteractionSummary executeNativeVolumeInteractions(
    std::vector<Particle::ParticlePtr>& particles, double dt, NativeVolumeParityRoute route)
{
    const auto interaction_families = resolveNativeVolumeInteractionFamilies(route);
    return executeNativeVolumeInteractions(particles, dt, interaction_families);
}

NativeVolumeInteractionSummary executeNativeVolumeInteractions(
    std::vector<Particle::ParticlePtr>& particles, double dt,
    const std::vector<NativeVolumeInteractionFamily>& active_families)
{
    NativeVolumeInteractionSummary summary;
    const auto normalized_families = normalizeNativeVolumeInteractionFamilies(active_families);
    summary.family_count = normalized_families.size();
    summary.family_signature = nativeVolumeInteractionFamilySignature(normalized_families);
    if (particles.empty() || dt <= 0.0)
    {
        return summary;
    }

    for (auto* particle : particles)
    {
        if (particle && particle->isActive())
        {
            ++summary.active_particle_count;
        }
    }

    auto manager = Field::buildNativeVolumeInteractionManager(normalized_families);
    summary.interactor_count = manager.getInteractorCount();
    manager.executeAllInteractions(particles, dt);
    summary.executed = summary.interactor_count > 0 && summary.active_particle_count > 0;
    summary.executed_interaction_count = summary.executed ? summary.interactor_count : 0;
    return summary;
}

NativeVolumeParitySnapshot buildNativeVolumeParitySnapshot(
    const Mesh::VolMeshPtr& mesh, std::vector<Particle::ParticlePtr>& particles,
    double lower_dirichlet_v, double upper_dirichlet_v, double dt, NativeVolumeParityRoute route)
{
    NativeVolumeParitySnapshot snapshot;
    snapshot.resolved_route = route;
    if (!mesh)
    {
        return snapshot;
    }

    snapshot.mesh = summarizeNativeVolumeMesh(*mesh);
    snapshot.field = solveNativeVolumePotential(mesh, lower_dirichlet_v, upper_dirichlet_v, route);
    snapshot.distribution =
        evaluateNativeVolumeDistributions(mesh->getNodeCount(), particles, route);
    snapshot.interaction = executeNativeVolumeInteractions(particles, dt, route);
    return snapshot;
}

} // namespace FieldSolver
} // namespace SCDAT
