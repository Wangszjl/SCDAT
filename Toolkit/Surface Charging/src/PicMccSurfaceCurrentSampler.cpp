#include "PicMccSurfaceCurrentSampler.h"
#include "SurfaceFlowCouplingModel.h"

#include "../../Tools/Basic/include/StringTokenUtils.h"
#include "../../Tools/Boundary/include/BoundaryConditions.h"
#include "../../Tools/FieldSolver/include/PoissonSolver.h"
#include "../../Tools/Interactions/Collisions/include/CollisionAlgorithm.h"
#include "../../Tools/Interactions/Collisions/include/CollisionPicAdapter.h"
#include "../../Tools/Mesh/include/MeshParsing.h"
#include "../../Tools/PICcore/include/PICCycle.h"
#include "../../Tools/Particle/include/ParticleManager.h"
#include "../../Tools/Particle/include/SurfaceDistributionFunction.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{
using SCDAT::Collision::CollisionParameters;
using SCDAT::Collision::CollisionProcessFactory;
using SCDAT::Collision::MonteCarloCollisionHandler;
using SCDAT::Collision::ParticleSpecies;
using SCDAT::Collision::CollisionCrossSectionSetId;
using SCDAT::Geometry::Point3D;
using SCDAT::Geometry::Vector3D;
using SCDAT::Mesh::GeometryDimensions;
using SCDAT::Mesh::GeometryMeshGenerator;
using SCDAT::Mesh::GeometryShape;
using SCDAT::Mesh::MeshGenerationOptions;
using SCDAT::Mesh::VolMesh;
using SCDAT::PICcore::PICCycle;
using BoundaryFace = SCDAT::PICcore::PICCycle::BoundaryFace;

constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kAtomicMassUnit = 1.66053906660e-27;
constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kPi = 3.14159265358979323846;

bool useDeterministicSamplingPolicy(const std::string& sampling_policy)
{
    const std::string token = Basic::normalizeAlnumToken(sampling_policy);
    return token.empty() || token == "deterministic" || token == "fixed" ||
           token == "replay" || token == "reproducible";
}

unsigned int hashCombine(unsigned int state, unsigned int value)
{
    constexpr unsigned int kOffset = 2166136261u;
    constexpr unsigned int kPrime = 16777619u;
    state = (state == 0u) ? kOffset : state;
    state ^= value;
    state *= kPrime;
    return state;
}

unsigned int resolveSamplerSeed(const PicMccSurfaceSamplerConfig& config)
{
    if (!useDeterministicSamplingPolicy(config.sampling_policy))
    {
        std::random_device rd;
        return rd();
    }

    unsigned int seed = config.seed == 0u ? 20260408u : config.seed;
    seed = hashCombine(seed, static_cast<unsigned int>(config.window_steps));
    seed = hashCombine(seed, static_cast<unsigned int>(config.z_layers));
    seed = hashCombine(seed, static_cast<unsigned int>(config.particles_per_element));
    seed = hashCombine(seed, static_cast<unsigned int>(config.enable_mcc ? 1 : 0));
    seed = hashCombine(seed, static_cast<unsigned int>(config.linked_boundary_face_count));
    seed = hashCombine(seed, static_cast<unsigned int>(std::lround(
                                1000.0 * std::clamp(config.projection_weight_sum, 0.0, 1000.0))));
    seed = hashCombine(seed, static_cast<unsigned int>(std::lround(
                                1000.0 * std::clamp(config.surface_aspect_ratio, 0.1, 10.0))));
    for (const char ch : config.node_name)
    {
        seed = hashCombine(seed, static_cast<unsigned int>(static_cast<unsigned char>(ch)));
    }
    for (const char ch : config.boundary_group_id)
    {
        seed = hashCombine(seed, static_cast<unsigned int>(static_cast<unsigned char>(ch)));
    }
    return seed;
}

std::vector<std::string> resolveOrderedSourceKeys(const PicMccSurfaceSamplerConfig& config)
{
    std::vector<std::string> resolved;
    resolved.reserve(config.source_keys.size());
    for (const auto& raw_key : config.source_keys)
    {
        const auto source_key = Basic::trimAscii(raw_key);
        if (source_key.empty())
        {
            continue;
        }
        if (std::find(resolved.begin(), resolved.end(), source_key) == resolved.end())
        {
            resolved.push_back(source_key);
        }
    }
    return resolved;
}

double topologySignature(const PicMccSurfaceSamplerConfig& config)
{
    double signature = static_cast<double>(config.linked_boundary_face_count);
    signature += 0.01 * std::clamp(config.projection_weight_sum, 0.0, 1000.0);
    signature += 0.001 * std::clamp(config.surface_aspect_ratio, 0.1, 10.0);
    for (const char ch : config.node_name)
    {
        signature += static_cast<double>(static_cast<unsigned char>(ch)) * 1.0e-6;
    }
    for (const char ch : config.boundary_group_id)
    {
        signature += static_cast<double>(static_cast<unsigned char>(ch)) * 1.0e-7;
    }
    return signature;
}

struct DepositionSchemeResolution
{
    Particle::TrajectoryChargeDepositionKernel kernel =
        Particle::TrajectoryChargeDepositionKernel::LinearSegmentCloud;
    std::size_t segment_count = 10;
    std::string kernel_name = "linear_segment_cloud";
};

DepositionSchemeResolution resolveDepositionScheme(const PicMccSurfaceSamplerConfig& config)
{
    DepositionSchemeResolution resolved;
    resolved.segment_count =
        std::clamp<std::size_t>(std::max<std::size_t>(4, config.window_steps), 1, 64);

    const std::string deposition = Basic::normalizeAlnumToken(config.deposition_scheme);
    if (deposition == "picwindowngp" || deposition == "ngp" || deposition == "nearest")
    {
        resolved.kernel = Particle::TrajectoryChargeDepositionKernel::NearestPoint;
        resolved.segment_count = 1;
        resolved.kernel_name = "nearest_point";
        return resolved;
    }

    if (deposition == "picwindowtsc" || deposition == "tsc" || deposition == "quadratic" ||
        deposition == "highorder")
    {
        resolved.kernel = Particle::TrajectoryChargeDepositionKernel::QuadratureGauss3;
        resolved.segment_count =
            std::clamp<std::size_t>(std::max<std::size_t>(8, config.window_steps), 2, 64);
        resolved.kernel_name = "quadrature_gauss3";
        return resolved;
    }

    resolved.kernel = Particle::TrajectoryChargeDepositionKernel::LinearSegmentCloud;
    resolved.kernel_name = "linear_segment_cloud";
    return resolved;
}

double computeDebyeLength(const PicMccSurfaceSamplerConfig& config)
{
    const double electron_temperature_j =
        std::max(1.0e-3, config.plasma.electron_temperature_ev) * kElementaryCharge;
    const double electron_density = std::max(1.0e3, config.plasma.electron_density_m3);
    return std::sqrt(kEpsilon0 * electron_temperature_j /
                     (electron_density * kElementaryCharge * kElementaryCharge));
}

double computeStableTimeStep(const PicMccSurfaceSamplerConfig& config)
{
    const double gap = std::max(1.0e-6, config.gap_distance_m);
    const double electron_temperature_j =
        std::max(1.0e-3, config.plasma.electron_temperature_ev) * kElementaryCharge;
    const double ion_mass = std::max(1.0, config.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double electron_speed = std::sqrt(electron_temperature_j / std::max(1.0e-32, kElectronMass));
    const double ion_speed = std::sqrt(
        std::max(1.0e-3, config.plasma.ion_temperature_ev) * kElementaryCharge /
        std::max(1.0e-32, ion_mass));
    const double drift_speed =
        std::max(0.0, config.bulk_flow_velocity_m_per_s) + std::max(0.0, config.ion_directed_velocity_m_per_s);
    const double max_speed = std::max({1.0, electron_speed, ion_speed + drift_speed});
    return std::clamp(gap / (40.0 * max_speed), 5.0e-12, 5.0e-9);
}

std::shared_ptr<VolMesh> createSamplingMesh(const PicMccSurfaceSamplerConfig& config)
{
    GeometryDimensions dims;
    dims.shape = GeometryShape::PLATE;
    const double topology_area_scale =
        std::max(1.0, std::sqrt(std::max(1.0, static_cast<double>(config.linked_boundary_face_count))) *
                          std::max(1.0, config.projection_weight_sum));
    const double effective_area_m2 = std::max(1.0e-12, config.surface_area_m2) * topology_area_scale;
    const double aspect_ratio = std::clamp(config.surface_aspect_ratio, 0.25, 4.0);
    dims.length = std::sqrt(effective_area_m2 * aspect_ratio);
    dims.width = std::sqrt(effective_area_m2 / aspect_ratio);
    dims.thickness = std::max(1.0e-6, config.gap_distance_m);

    MeshGenerationOptions options;
    options.nx = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(std::sqrt(
               std::max(1.0, static_cast<double>(config.linked_boundary_face_count))))));
    options.ny = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(static_cast<double>(options.nx) /
                                              std::max(0.5, aspect_ratio))));
    options.nz = std::max<std::size_t>(
        4, config.z_layers + std::min<std::size_t>(8, config.linked_boundary_face_count / 2));
    options.tetrahedralize = true;

    auto mesh = GeometryMeshGenerator::generateFromDimensions(dims, options);
    if (!mesh || !mesh->validate())
    {
        return nullptr;
    }

    const double z_min = mesh->getBoundingBoxMin().z();
    const double z_max = mesh->getBoundingBoxMax().z();
    for (const auto& node : mesh->getNodes())
    {
        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12 || std::abs(z - z_max) < 1.0e-12)
        {
            node->setBoundaryType(SCDAT::Mesh::BoundaryType::DIRICHLET);
        }
        else
        {
            node->setBoundaryType(SCDAT::Mesh::BoundaryType::INTERIOR);
        }
    }

    return std::shared_ptr<VolMesh>(std::move(mesh));
}

Vector3D sampleVelocity(double sigma, double drift_z, bool toward_surface, std::mt19937& rng)
{
    std::normal_distribution<double> normal(0.0, sigma);
    double vz = toward_surface ? -std::abs(normal(rng)) : std::abs(normal(rng));
    vz += drift_z;
    if (toward_surface && vz > 0.0)
    {
        vz = -0.25 * std::abs(vz);
    }
    return Vector3D(normal(rng), normal(rng), vz);
}

Point3D samplePointInTriangle(const Point3D& a, const Point3D& b, const Point3D& c,
                              std::mt19937& rng)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double r1 = std::sqrt(dist(rng));
    const double r2 = dist(rng);

    const double w0 = 1.0 - r1;
    const double w1 = r1 * (1.0 - r2);
    const double w2 = r1 * r2;
    return Point3D(w0 * a.x() + w1 * b.x() + w2 * c.x(),
                   w0 * a.y() + w1 * b.y() + w2 * c.y(),
                   w0 * a.z() + w1 * b.z() + w2 * c.z());
}

Point3D samplePointInTetra(const Point3D& a, const Point3D& b, const Point3D& c,
                           const Point3D& d, std::mt19937& rng)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double e0 = -std::log(std::max(1.0e-12, dist(rng)));
    const double e1 = -std::log(std::max(1.0e-12, dist(rng)));
    const double e2 = -std::log(std::max(1.0e-12, dist(rng)));
    const double e3 = -std::log(std::max(1.0e-12, dist(rng)));
    const double s = e0 + e1 + e2 + e3;

    const double w0 = e0 / s;
    const double w1 = e1 / s;
    const double w2 = e2 / s;
    const double w3 = e3 / s;
    return Point3D(w0 * a.x() + w1 * b.x() + w2 * c.x() + w3 * d.x(),
                   w0 * a.y() + w1 * b.y() + w2 * c.y() + w3 * d.y(),
                   w0 * a.z() + w1 * b.z() + w2 * c.z() + w3 * d.z());
}

Point3D samplePointInElement(const VolMesh& mesh, std::size_t element_index, std::mt19937& rng)
{
    const auto& element = mesh.getElement(element_index);
    const auto& nodes = element->getNodes();
    if (nodes.size() == 4 && nodes[0] && nodes[1] && nodes[2] && nodes[3])
    {
        return samplePointInTetra(nodes[0]->getPosition(), nodes[1]->getPosition(),
                                  nodes[2]->getPosition(), nodes[3]->getPosition(), rng);
    }
    if (nodes.size() == 3 && nodes[0] && nodes[1] && nodes[2])
    {
        return samplePointInTriangle(nodes[0]->getPosition(), nodes[1]->getPosition(),
                                     nodes[2]->getPosition(), rng);
    }
    return mesh.getRandomPointInElement(element_index);
}

Vector3D sampleVelocityFromSpectrum(const Particle::ResolvedSpectrum& spectrum,
                                    double fallback_temperature_ev, double fallback_mass_kg,
                                    double directed_drift_z, bool toward_surface, std::mt19937& rng)
{
    if (!spectrum.populations.empty())
    {
        double total_density = 0.0;
        for (const auto& population : spectrum.populations)
        {
            total_density += std::max(0.0, population.density_m3);
        }

        std::uniform_real_distribution<double> u(0.0, std::max(1.0e-18, total_density));
        const double pick = u(rng);
        double running = 0.0;
        const Particle::SpectrumPopulation* selected = &spectrum.populations.front();
        for (const auto& population : spectrum.populations)
        {
            running += std::max(0.0, population.density_m3);
            if (pick <= running)
            {
                selected = &population;
                break;
            }
        }

        const double mass_kg = std::max(1.0e-32, selected->mass_amu * kAtomicMassUnit);
        const double sigma =
            std::sqrt(std::max(1.0e-3, selected->temperature_ev) * kElementaryCharge / mass_kg);
        return sampleVelocity(sigma, directed_drift_z + selected->drift_speed_m_per_s, toward_surface, rng);
    }

    if (!spectrum.energy_grid_ev.empty() &&
        spectrum.energy_grid_ev.size() == spectrum.differential_number_flux.size())
    {
        const double energy_ev = Particle::sampleResolvedSpectrumEnergyEv(
            spectrum, rng, fallback_temperature_ev);
        const double mass_kg =
            Particle::resolvedSpectrumAverageMassKg(spectrum, fallback_mass_kg);
        const double speed =
            std::sqrt(2.0 * std::max(1.0e-3, energy_ev) * kElementaryCharge / mass_kg);
        const double sigma = speed / std::sqrt(3.0);
        return sampleVelocity(sigma, directed_drift_z, toward_surface, rng);
    }

    const double sigma =
        std::sqrt(std::max(1.0e-3, fallback_temperature_ev) * kElementaryCharge / fallback_mass_kg);
    return sampleVelocity(sigma, directed_drift_z, toward_surface, rng);
}

void initializeParticles(SCDAT::Particle::ParticleManager& particle_manager, const VolMesh& mesh,
                         const PicMccSurfaceSamplerConfig& config, std::mt19937& rng)
{
    const double electron_sigma = std::sqrt(
        std::max(1.0e-3, config.plasma.electron_temperature_ev) * kElementaryCharge / kElectronMass);
    const double ion_mass = std::max(1.0, config.plasma.ion_mass_amu) * kAtomicMassUnit;
    const double ion_sigma = std::sqrt(
        std::max(1.0e-3, config.plasma.ion_temperature_ev) * kElementaryCharge / ion_mass);
    const double volume = std::max(1.0e-18, config.surface_area_m2 * config.gap_distance_m);
    const double macro_particles_per_species =
        static_cast<double>(mesh.getElementCount() * std::max<std::size_t>(1, config.particles_per_element));
    const double electron_density =
        config.has_electron_spectrum
            ? Particle::resolvedSpectrumDensity(config.electron_spectrum,
                                               config.plasma.electron_density_m3)
            : config.plasma.electron_density_m3;
    const double ion_density =
        config.has_ion_spectrum ? Particle::resolvedSpectrumDensity(
                                      config.ion_spectrum, config.plasma.ion_density_m3)
                                : config.plasma.ion_density_m3;
    const double electron_weight =
        electron_density * volume / std::max(1.0, macro_particles_per_species);
    const double ion_weight =
        ion_density * volume / std::max(1.0, macro_particles_per_species);

    const auto flow_projection = resolveSurfaceFlowProjection(
        config.bulk_flow_velocity_m_per_s, config.flow_alignment_cosine, config.flow_angle_deg);
    const double toward_surface_flow_speed =
        std::max(0.0, flow_projection.patch_projected_speed_m_per_s);
    const double ion_drift_z =
        -(toward_surface_flow_speed + std::max(0.0, config.ion_directed_velocity_m_per_s));
    const double electron_drift_z =
        -toward_surface_flow_speed * std::clamp(config.electron_flow_coupling, 0.0, 1.0);

    for (std::size_t element_index = 0; element_index < mesh.getElementCount(); ++element_index)
    {
        for (std::size_t n = 0; n < std::max<std::size_t>(1, config.particles_per_element); ++n)
        {
            const Vector3D electron_velocity =
                config.has_electron_spectrum
                    ? sampleVelocityFromSpectrum(config.electron_spectrum,
                                                config.plasma.electron_temperature_ev, kElectronMass,
                                                electron_drift_z, true, rng)
                    : sampleVelocity(electron_sigma, electron_drift_z, true, rng);
            const Vector3D ion_velocity =
                config.has_ion_spectrum
                    ? sampleVelocityFromSpectrum(config.ion_spectrum, config.plasma.ion_temperature_ev,
                                                ion_mass, ion_drift_z, true, rng)
                    : sampleVelocity(ion_sigma, ion_drift_z, true, rng);

            particle_manager.createElectron(samplePointInElement(mesh, element_index, rng),
                                            electron_velocity, electron_weight);
            particle_manager.createIon(samplePointInElement(mesh, element_index, rng), ion_velocity,
                                       static_cast<int>(std::max(
                                           1.0, config.has_ion_spectrum
                                                ? Particle::resolvedSpectrumAverageMassKg(
                                                    config.ion_spectrum, ion_mass) /
                                                          kAtomicMassUnit
                                                    : config.plasma.ion_mass_amu)),
                                       1, ion_weight);
        }
    }
}

void applyBoundaryPotentials(PICCycle& cycle, const VolMesh& mesh, double surface_potential_v,
                             double plasma_reference_potential_v)
{
    const double z_min = mesh.getBoundingBoxMin().z();
    const double z_max = mesh.getBoundingBoxMax().z();
    for (const auto& node : mesh.getNodes())
    {
        const double z = node->getPosition().z();
        if (std::abs(z - z_min) < 1.0e-12)
        {
            cycle.setNodePotential(node->getId(), surface_potential_v);
        }
        else if (std::abs(z - z_max) < 1.0e-12)
        {
            cycle.setNodePotential(node->getId(), plasma_reference_potential_v);
        }
    }
}

void applyCollisionStep(SCDAT::Particle::ParticleManager& particle_manager,
                        MonteCarloCollisionHandler& collision_handler, double dt)
{
    std::vector<SCDAT::Collision::ParticleObject> working_particles;
    const auto live_particles = particle_manager.getAllParticles();
    working_particles.reserve(live_particles.size());
    for (const auto* particle : live_particles)
    {
        if (particle)
        {
            working_particles.push_back(*particle);
        }
    }

    auto created_particles = collision_handler.processCollisions(working_particles, dt);
    for (const auto& updated_particle : working_particles)
    {
        if (auto* particle = particle_manager.getParticle(updated_particle.getId()))
        {
            *particle = updated_particle;
        }
    }
    if (!created_particles.empty())
    {
        particle_manager.getContainer().addParticles(created_particles);
    }
    particle_manager.removeInactiveParticles();
}

MonteCarloCollisionHandler makeCollisionHandler(const PicMccSurfaceSamplerConfig& config,
                                                unsigned int seed)
{
    MonteCarloCollisionHandler collision_handler(seed);
    collision_handler.clearCollisionProcesses();

    if (!config.enable_mcc || config.plasma.neutral_density_m3 <= 1.0e6)
    {
        return collision_handler;
    }

    auto cross_section_set_id =
        SCDAT::Collision::parseCollisionCrossSectionSetId(config.collision_cross_section_set_id)
            .value_or(CollisionCrossSectionSetId::SurfacePicV1);
    SCDAT::Collision::initializeConfiguredCollisionHandler(
        collision_handler, cross_section_set_id, config.plasma.neutral_density_m3,
        config.plasma.ion_density_m3, config.plasma.electron_density_m3);

    return collision_handler;
}

PicMccCurrentSample sampleSingleWindow(const PicMccSurfaceSamplerConfig& config)
{
    PicMccCurrentSample sample;
    sample.sampled_node_name = config.node_name;
    sample.sampled_boundary_group_id = config.boundary_group_id;
    sample.topology_signature = topologySignature(config);
    sample.topology_area_scale =
        std::max(1.0, std::sqrt(std::max(1.0, static_cast<double>(config.linked_boundary_face_count))) *
                          std::max(1.0, config.projection_weight_sum));
    sample.sampled_surface_potential_v = config.surface_potential_v;
    sample.mcc_enabled = config.enable_mcc && config.plasma.neutral_density_m3 > 1.0e6;
    sample.sampling_policy_resolved =
        useDeterministicSamplingPolicy(config.sampling_policy) ? "deterministic" : "stochastic";
    sample.seed_used = resolveSamplerSeed(config);

    const DepositionSchemeResolution deposition = resolveDepositionScheme(config);
    sample.deposition_segments = deposition.segment_count;
    sample.deposition_kernel = deposition.kernel_name;

    auto mesh = createSamplingMesh(config);
    if (!mesh)
    {
        return sample;
    }

    auto particle_manager = std::make_shared<SCDAT::Particle::ParticleManager>();
    std::mt19937 rng(sample.seed_used);
    initializeParticles(*particle_manager, *mesh, config, rng);

    auto poisson_solver = std::make_shared<SCDAT::FieldSolver::PoissonSolver>(mesh);
    PICCycle::Parameters cycle_params;
    cycle_params.time_step = computeStableTimeStep(config);
    cycle_params.max_iterations = static_cast<int>(std::max<std::size_t>(1, config.window_steps));
    cycle_params.statistics_frequency = std::max(1, static_cast<int>(config.window_steps / 4));
    cycle_params.charge_conservation = true;
    cycle_params.verbose = false;
    cycle_params.trajectory_charge_deposition_segments = deposition.segment_count;
    cycle_params.trajectory_charge_deposition_kernel = deposition.kernel;

    PICCycle cycle(cycle_params);
    cycle.setMesh(mesh->getNodes(), mesh->getElements());
    cycle.setParticleManager(particle_manager);
    cycle.setPoissonSolver(poisson_solver);

    const Point3D domain_min = mesh->getBoundingBoxMin();
    const Point3D domain_max = mesh->getBoundingBoxMax();
    std::vector<Vector3D> periodic_vectors;
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
        auto periodic_boundary =
            std::make_shared<SCDAT::Particle::PeriodicBoundaryCondition>(periodic_vectors);
        cycle.setBoundaryCondition(BoundaryFace::XMin, periodic_boundary);
        cycle.setBoundaryCondition(BoundaryFace::XMax, periodic_boundary);
        cycle.setBoundaryCondition(BoundaryFace::YMin, periodic_boundary);
        cycle.setBoundaryCondition(BoundaryFace::YMax, periodic_boundary);
    }

    cycle.setBoundaryCondition(BoundaryFace::ZMin,
                               std::make_shared<SCDAT::Particle::AbsorbingBoundaryCondition>(1.0));
    cycle.setBoundaryCondition(BoundaryFace::ZMax,
                               std::make_shared<SCDAT::Particle::AbsorbingBoundaryCondition>(1.0));
    PICCycle::SurfaceBoundaryMetadata metadata;
    metadata.surface_id = 0;
    metadata.material_id = static_cast<int>(config.material.getId());
    metadata.circuit_node_id = 0;
    cycle.setBoundaryMetadata(BoundaryFace::ZMin, metadata);

    if (!cycle.initialize())
    {
        return sample;
    }

    applyBoundaryPotentials(cycle, *mesh, config.surface_potential_v, config.plasma_reference_potential_v);
    MonteCarloCollisionHandler collision_handler =
        makeCollisionHandler(config, hashCombine(sample.seed_used, 0x9e3779b9u));

    std::size_t executed_steps = 0;
    for (std::size_t step = 0; step < std::max<std::size_t>(1, config.window_steps); ++step)
    {
        applyBoundaryPotentials(cycle, *mesh, config.surface_potential_v, config.plasma_reference_potential_v);
        cycle.executeTimeStep();
        ++executed_steps;

        if (config.enable_mcc)
        {
            applyCollisionStep(*particle_manager, collision_handler, cycle_params.time_step);
        }
        particle_manager->updateParticleAges(cycle_params.time_step);
    }

    const auto& ledger =
        cycle.getSurfaceCurrentLedger()[static_cast<std::size_t>(BoundaryFace::ZMin)];
    const double total_time =
        cycle_params.time_step * static_cast<double>(std::max<std::size_t>(1, executed_steps));
    const double area = std::max(1.0e-12, config.surface_area_m2);

    sample.electron_collection_current_density_a_per_m2 =
        ledger.absorbed_electron_charge_c / (total_time * area);
    sample.ion_collection_current_density_a_per_m2 =
        ledger.absorbed_ion_charge_c / (total_time * area);
    sample.emitted_electron_current_density_a_per_m2 =
        -ledger.emitted_electron_charge_c / (total_time * area);
    sample.net_collection_current_density_a_per_m2 =
        sample.electron_collection_current_density_a_per_m2 +
        sample.ion_collection_current_density_a_per_m2 +
        sample.emitted_electron_current_density_a_per_m2;
    sample.total_collisions = collision_handler.getStatistics().total_collisions;

    const auto source_keys = resolveOrderedSourceKeys(config);
    if (!source_keys.empty())
    {
        sample.source_resolved_samples.reserve(source_keys.size());
        for (const auto& source_key : source_keys)
        {
            PicMccSourceResolvedSample source_sample;
            source_sample.source_key = source_key;
            sample.source_resolved_samples.push_back(std::move(source_sample));
        }

        if (sample.source_resolved_samples.size() == 1)
        {
            auto& source_sample = sample.source_resolved_samples.front();
            source_sample.collected_current_density_a_per_m2 =
                sample.electron_collection_current_density_a_per_m2 +
                sample.ion_collection_current_density_a_per_m2;
            source_sample.superparticle_count =
                static_cast<double>(particle_manager->getAllParticles().size());
            source_sample.strictly_attributed = true;
        }
    }

    sample.valid =
        executed_steps > 0 && std::isfinite(sample.net_collection_current_density_a_per_m2);
    return sample;
}

} // namespace

PicMccCurrentSample PicMccSurfaceCurrentSampler::sample(
    const PicMccSurfaceSamplerConfig& config) const
{
    return sampleSingleWindow(config);
}

PicMccCurrentSample PicMccSurfaceCurrentSampler::sampleWithDerivative(
    const PicMccSurfaceSamplerConfig& config, double probe_delta_v) const
{
    PicMccCurrentSample center = sampleSingleWindow(config);
    if (!center.valid)
    {
        return center;
    }

    const double delta = std::max(0.1, std::abs(probe_delta_v));
    PicMccSurfaceSamplerConfig plus = config;
    PicMccSurfaceSamplerConfig minus = config;
    plus.surface_potential_v += delta;
    minus.surface_potential_v -= delta;

    const auto sample_plus = sampleSingleWindow(plus);
    const auto sample_minus = sampleSingleWindow(minus);
    if (sample_plus.valid && sample_minus.valid)
    {
        center.current_derivative_a_per_m2_per_v =
            (sample_plus.net_collection_current_density_a_per_m2 -
             sample_minus.net_collection_current_density_a_per_m2) /
            (2.0 * delta);
        center.total_collisions =
            std::max({center.total_collisions, sample_plus.total_collisions,
                      sample_minus.total_collisions});
    }
    center.mcc_enabled = center.mcc_enabled || sample_plus.mcc_enabled || sample_minus.mcc_enabled;
    return center;
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
