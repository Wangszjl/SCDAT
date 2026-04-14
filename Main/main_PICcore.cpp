#include "Basic/include/Constants.h"
#include "FieldSolver/include/PoissonSolver.h"
#include "Interactions/Collisions/include/CollisionPicAdapter.h"
#include "Mesh/include/MeshParsing.h"
#include "Output/include/ResultExporter.h"
#include "PICcore/include/PICCycle.h"
#include "Particle/include/ParticleManager.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using SCDAT::Basic::Constants::MathConstants;
using SCDAT::Basic::Constants::PhysicsConstants;
using SCDAT::Collision::BackgroundMccRuntimeConfig;
using SCDAT::Collision::CollisionCrossSectionSetId;
using SCDAT::Collision::MonteCarloCollisionHandler;
using SCDAT::Geometry::Point3D;
using SCDAT::Geometry::Vector3D;
using SCDAT::Particle::AbsorbingBoundaryCondition;
using SCDAT::Mesh::ElementPtr;
using SCDAT::Mesh::GeometryDimensions;
using SCDAT::Mesh::GeometryMeshGenerator;
using SCDAT::Mesh::GeometryShape;
using SCDAT::Mesh::MeshGenerationOptions;
using SCDAT::Mesh::NodePtr;
using SCDAT::Mesh::VolMesh;
using SCDAT::Particle::PeriodicBoundaryCondition;
using SCDAT::FieldSolver::PoissonSolver;
using SCDAT::Output::ColumnarDataSet;
using SCDAT::Output::ResultExporter;
using SCDAT::Particle::ParticleManager;
using SCDAT::PICcore::PICCycle;

using LayerMap = std::map<double, std::vector<std::size_t>>;
using BoundaryFace = SCDAT::PICcore::PICCycle::BoundaryFace;

struct BenchmarkConfig
{
    double gap_distance = 0.067;
    double electrode_area = 0.0314;
    double pressure = 10.0;
    double neutral_density = 2.4e20;
    double rf_voltage = 450.0;
    double rf_frequency = 13.56e6;
    int steps_per_cycle = 200;
    int rf_cycles = 2;
    std::size_t nz = 40;
    std::size_t particles_per_element = 4;
    double electron_temperature_ev = 9.36;
    double ion_temperature_ev = 0.026;
    double initial_plasma_density = 1.4e14;
    std::string output_csv = "results/pic_mcc_turner_minimal.csv";
    std::string collision_cross_section_set_id = "background_mcc_v1";
    bool collision_reinitialize_each_step = false;

    [[nodiscard]] double electrodeSideLength() const
    {
        return std::sqrt(electrode_area);
    }

    [[nodiscard]] double timeStep() const
    {
        return 1.0 / (rf_frequency * static_cast<double>(steps_per_cycle));
    }

    [[nodiscard]] int totalSteps() const
    {
        return rf_cycles * steps_per_cycle;
    }
};

struct LayerProfile
{
    std::vector<double> z_positions;
    std::vector<double> potentials;
    std::vector<double> electric_field_z;
    std::vector<double> electron_density;
    std::vector<double> ion_density;
    std::vector<double> electron_temperature_ev;
};

double roundLayerKey(double value)
{
    return static_cast<double>(std::llround(value * 1.0e9)) * 1.0e-9;
}

LayerMap buildLayerMap(const std::vector<NodePtr>& nodes)
{
    LayerMap layers;
    for (const auto& node : nodes)
    {
        layers[roundLayerKey(node->getPosition().z())].push_back(node->getId());
    }
    return layers;
}

std::vector<double> makeLinearPotentialProfile(const std::vector<NodePtr>& nodes, double left_value,
                                               double right_value, double gap_distance)
{
    std::vector<double> potential(nodes.size(), 0.0);
    if (gap_distance <= 0.0)
    {
        return potential;
    }

    for (const auto& node : nodes)
    {
        const double z_fraction = std::clamp(node->getPosition().z() / gap_distance, 0.0, 1.0);
        potential[node->getId()] = left_value + (right_value - left_value) * z_fraction;
    }

    return potential;
}

void applyBoundaryPotentials(PICCycle& cycle, const LayerMap& layers, double left_value,
                             double right_value)
{
    if (layers.empty())
    {
        return;
    }

    const auto& left_nodes = layers.begin()->second;
    const auto& right_nodes = layers.rbegin()->second;

    for (std::size_t node_id : left_nodes)
    {
        cycle.setNodePotential(node_id, left_value);
    }
    for (std::size_t node_id : right_nodes)
    {
        cycle.setNodePotential(node_id, right_value);
    }
}

std::shared_ptr<VolMesh> createBenchmarkMesh(const BenchmarkConfig& config)
{
    GeometryDimensions dims;
    dims.shape = GeometryShape::PLATE;
    dims.length = config.electrodeSideLength();
    dims.width = config.electrodeSideLength();
    dims.thickness = config.gap_distance;

    MeshGenerationOptions options;
    options.nx = 1;
    options.ny = 1;
    options.nz = config.nz;
    options.tetrahedralize = true;

    auto mesh = GeometryMeshGenerator::generateFromDimensions(dims, options);
    if (!mesh || !mesh->validate())
    {
        throw std::runtime_error("Failed to generate a valid PIC benchmark mesh.");
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

Vector3D sampleMaxwellVelocity(double component_sigma, std::mt19937& rng)
{
    std::normal_distribution<double> normal(0.0, component_sigma);
    return Vector3D(normal(rng), normal(rng), normal(rng));
}

void initializeBenchmarkParticles(ParticleManager& particle_manager, const VolMesh& mesh,
                                  const BenchmarkConfig& config)
{
    std::mt19937 rng(42);

    const double electron_component_sigma = std::sqrt(config.electron_temperature_ev *
                                                      PhysicsConstants::EVToJoule /
                                                      PhysicsConstants::ElectronMass);
    const double argon_ion_mass = 40.0 * PhysicsConstants::AtomicMassUnit;
    const double ion_component_sigma =
        std::sqrt(config.ion_temperature_ev * PhysicsConstants::EVToJoule / argon_ion_mass);

    const std::size_t macro_particles_per_species =
        mesh.getElementCount() * config.particles_per_element;
    const double discharge_volume = config.electrode_area * config.gap_distance;
    const double macro_weight = macro_particles_per_species > 0
                                    ? (config.initial_plasma_density * discharge_volume /
                                       static_cast<double>(macro_particles_per_species))
                                    : 1.0;

    for (std::size_t element_index = 0; element_index < mesh.getElementCount(); ++element_index)
    {
        for (std::size_t n = 0; n < config.particles_per_element; ++n)
        {
            const Point3D electron_position = mesh.getRandomPointInElement(element_index);
            const Point3D ion_position = mesh.getRandomPointInElement(element_index);

            particle_manager.createElectron(electron_position,
                                            sampleMaxwellVelocity(electron_component_sigma, rng),
                                            macro_weight);
            particle_manager.createIon(ion_position, sampleMaxwellVelocity(ion_component_sigma, rng),
                                       40, 1, macro_weight);
        }
    }
}

LayerProfile collectLayerProfile(const BenchmarkConfig& config, const LayerMap& layers,
                                 const PICCycle& cycle, const ParticleManager& particle_manager)
{
    LayerProfile profile;
    const auto& node_potential = cycle.getPotential();
    const auto& node_electric_field = cycle.getElectricField();

    profile.z_positions.reserve(layers.size());
    profile.potentials.reserve(layers.size());
    profile.electric_field_z.reserve(layers.size());
    profile.electron_density.assign(layers.size(), 0.0);
    profile.ion_density.assign(layers.size(), 0.0);
    profile.electron_temperature_ev.assign(layers.size(), 0.0);

    std::vector<double> electron_energy_accumulator(layers.size(), 0.0);

    std::size_t layer_index = 0;
    for (const auto& [z, node_ids] : layers)
    {
        profile.z_positions.push_back(z);

        double potential_sum = 0.0;
        double field_sum = 0.0;
        for (std::size_t node_id : node_ids)
        {
            if (node_id < node_potential.size())
            {
                potential_sum += node_potential[node_id];
            }
            if (node_id < node_electric_field.size())
            {
                field_sum += node_electric_field[node_id].z();
            }
        }

        profile.potentials.push_back(potential_sum / static_cast<double>(node_ids.size()));
        profile.electric_field_z.push_back(field_sum / static_cast<double>(node_ids.size()));
        ++layer_index;
    }

    const double dz = config.gap_distance / static_cast<double>(std::max<std::size_t>(1, config.nz));
    const double layer_volume = config.electrode_area * dz;
    const std::size_t max_layer_index = profile.z_positions.empty() ? 0 : profile.z_positions.size() - 1;

    for (const auto* particle : particle_manager.getAllParticles())
    {
        if (!particle)
        {
            continue;
        }

        const double z = std::clamp(particle->getPosition().z(), 0.0, config.gap_distance);
        const auto bucket = static_cast<std::size_t>(
            std::clamp<long long>(static_cast<long long>(std::llround(z / dz)), 0,
                                  static_cast<long long>(max_layer_index)));

        if (particle->getCharge() < 0.0)
        {
            profile.electron_density[bucket] += particle->getWeight() / layer_volume;
            electron_energy_accumulator[bucket] +=
                particle->getWeight() * particle->getKineticEnergy() * PhysicsConstants::JouleToEV;
        }
        else if (particle->getCharge() > 0.0)
        {
            profile.ion_density[bucket] += particle->getWeight() / layer_volume;
        }
    }

    for (std::size_t i = 0; i < profile.electron_temperature_ev.size(); ++i)
    {
        if (profile.electron_density[i] > 0.0)
        {
            const double macro_particle_count = profile.electron_density[i] * layer_volume;
            profile.electron_temperature_ev[i] = electron_energy_accumulator[i] /
                                                 std::max(1.0, macro_particle_count) / 1.5;
        }
    }

    return profile;
}

void writeLayerProfileResults(const std::string& filename, const LayerProfile& profile,
                              const BenchmarkConfig& config)
{
    ColumnarDataSet data_set;
    data_set.axis_name = "z_m";
    data_set.axis_values = profile.z_positions;
    data_set.scalar_series["potential_v"] = profile.potentials;
    data_set.scalar_series["electric_field_z_v_per_m"] = profile.electric_field_z;
    data_set.scalar_series["electron_density_m3"] = profile.electron_density;
    data_set.scalar_series["ion_density_m3"] = profile.ion_density;
    data_set.scalar_series["electron_temperature_ev"] = profile.electron_temperature_ev;
    data_set.metadata["module"] = "PIC-MCC";
    data_set.metadata["collision_cross_section_set_id"] = config.collision_cross_section_set_id;
    data_set.metadata["collision_runtime_source"] = "tools.interactions.collisions";
    data_set.metadata["collision_reinitialize_each_step"] =
        config.collision_reinitialize_each_step ? "true" : "false";

    ResultExporter exporter;
    if (auto result = exporter.exportDataSet(filename, data_set); !result)
    {
        throw std::runtime_error(result.message());
    }
}

void printRuntimeSummary(const BenchmarkConfig& config, const VolMesh& mesh,
                         const ParticleManager::Statistics& particle_stats,
                         const MonteCarloCollisionHandler::CollisionStatistics& collision_stats,
                         double wall_time_seconds, const std::string& output_csv,
                         double deposited_surface_charge)
{
    std::cout << "\n=== PIC-MCC Minimal Benchmark Summary ===\n";
    std::cout << "Pressure: " << config.pressure << " Pa\n";
    std::cout << "RF frequency: " << config.rf_frequency / 1.0e6 << " MHz\n";
    std::cout << "RF voltage amplitude: " << config.rf_voltage << " V\n";
    std::cout << "Collision cross-section set: " << config.collision_cross_section_set_id << '\n';
    std::cout << "Collision runtime source: tools.interactions.collisions\n";
    std::cout << "Mesh nodes/elements: " << mesh.getNodeCount() << " / " << mesh.getElementCount()
              << '\n';
    std::cout << "Time step: " << std::scientific << config.timeStep() << " s\n";
    std::cout << "Total steps: " << config.totalSteps() << '\n';
    std::cout << "Active particles: " << particle_stats.active_particles << '\n';
    std::cout << "Total collisions: " << collision_stats.total_collisions << '\n';
    std::cout << "Particles created by MCC: " << collision_stats.particles_created << '\n';
    std::cout << "Deposited surface charge: " << std::scientific << deposited_surface_charge
              << " C\n";
    std::cout << "Wall time: " << std::fixed << std::setprecision(3) << wall_time_seconds
              << " s\n";
    std::cout << "Profile csv: " << output_csv << '\n';
    std::cout << "=========================================\n";
}

} // namespace

int main(int argc, char* argv[])
{
    try
    {
        BenchmarkConfig config;
        if (argc > 1)
        {
            config.rf_cycles = std::max(1, std::stoi(argv[1]));
        }
        if (argc > 2)
        {
            config.nz = std::max<std::size_t>(4, static_cast<std::size_t>(std::stoul(argv[2])));
        }
        if (argc > 3)
        {
            config.particles_per_element =
                std::max<std::size_t>(1, static_cast<std::size_t>(std::stoul(argv[3])));
        }
        if (argc > 4)
        {
            config.output_csv = argv[4];
        }
        if (argc > 5)
        {
            config.collision_cross_section_set_id = argv[5];
        }
        if (argc > 6)
        {
            config.collision_reinitialize_each_step = std::stoi(argv[6]) != 0;
        }

        const auto cross_section_set_id =
            SCDAT::Collision::parseCollisionCrossSectionSetId(config.collision_cross_section_set_id);
        if (!cross_section_set_id)
        {
            throw std::runtime_error("Unsupported collision cross-section set id: " +
                                     config.collision_cross_section_set_id);
        }

        std::filesystem::create_directories(
            std::filesystem::path(config.output_csv).parent_path().empty()
                ? std::filesystem::path(".")
                : std::filesystem::path(config.output_csv).parent_path());

        auto mesh = createBenchmarkMesh(config);
        auto particle_manager = std::make_shared<ParticleManager>();
        initializeBenchmarkParticles(*particle_manager, *mesh, config);

        auto poisson_solver = std::make_shared<PoissonSolver>(mesh);

        PICCycle::Parameters cycle_params;
        cycle_params.time_step = config.timeStep();
        cycle_params.max_iterations = config.totalSteps();
        cycle_params.statistics_frequency = std::max(1, config.steps_per_cycle / 4);
        cycle_params.charge_conservation = true;
        cycle_params.verbose = false;

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
            auto periodic_boundary = std::make_shared<PeriodicBoundaryCondition>(periodic_vectors);
            cycle.setBoundaryCondition(BoundaryFace::XMin, periodic_boundary);
            cycle.setBoundaryCondition(BoundaryFace::XMax, periodic_boundary);
            cycle.setBoundaryCondition(BoundaryFace::YMin, periodic_boundary);
            cycle.setBoundaryCondition(BoundaryFace::YMax, periodic_boundary);
        }

        cycle.setBoundaryCondition(BoundaryFace::ZMin,
                                   std::make_shared<AbsorbingBoundaryCondition>(1.0));
        cycle.setBoundaryCondition(BoundaryFace::ZMax,
                                   std::make_shared<AbsorbingBoundaryCondition>(1.0));

        if (!cycle.initialize())
        {
            std::cerr << "PIC cycle initialization failed.\n";
            return 1;
        }

        const LayerMap layers = buildLayerMap(mesh->getNodes());
        cycle.setPotential(
            makeLinearPotentialProfile(mesh->getNodes(), 0.0, 0.0, config.gap_distance));

        MonteCarloCollisionHandler collision_handler(12345);
        SCDAT::Collision::initializeConfiguredCollisionHandler(
            collision_handler, *cross_section_set_id, config.neutral_density,
            config.initial_plasma_density, config.initial_plasma_density);

        BackgroundMccRuntimeConfig collision_runtime_config;
        collision_runtime_config.floor_densities.neutral_density_m3 = config.neutral_density;
        collision_runtime_config.floor_densities.ion_density_m3 = config.initial_plasma_density;
        collision_runtime_config.floor_densities.electron_density_m3 = config.initial_plasma_density;
        collision_runtime_config.cross_section_set_id = *cross_section_set_id;
        collision_runtime_config.reinitialize_configuration =
            config.collision_reinitialize_each_step;

        const double discharge_volume_m3 = config.electrode_area * config.gap_distance;

        std::cout << "Running minimal PIC-MCC benchmark\n";
        std::cout << "  cycles           : " << config.rf_cycles << '\n';
        std::cout << "  z layers         : " << config.nz << '\n';
        std::cout << "  particles/element: " << config.particles_per_element << '\n';
        std::cout << "  dt               : " << std::scientific << config.timeStep() << " s\n";
        std::cout << "  cross sections   : " << config.collision_cross_section_set_id << '\n';
        std::cout << "  total steps      : " << config.totalSteps() << "\n\n";

        auto wall_start = std::chrono::steady_clock::now();

        for (int step = 0; step < config.totalSteps(); ++step)
        {
            const double current_time = static_cast<double>(step) * config.timeStep();
            const double electrode_potential = config.rf_voltage *
                                               std::sin(MathConstants::TwoPi * config.rf_frequency *
                                                        current_time);

            applyBoundaryPotentials(cycle, layers, electrode_potential, 0.0);
            if (!cycle.executeTimeStep())
            {
                std::cerr << "PIC cycle terminated early at step " << step << ".\n";
                break;
            }

            const auto collision_result = SCDAT::Collision::executeBackgroundMccStep(
                particle_manager->getContainer(), collision_handler, config.timeStep(),
                discharge_volume_m3, collision_runtime_config);
            particle_manager->updateParticleAges(config.timeStep());

            if ((step + 1) % std::max(1, config.steps_per_cycle / 2) == 0 ||
                step + 1 == config.totalSteps())
            {
                const auto particle_stats = particle_manager->getStatistics();
                std::cout << "Step " << std::setw(4) << (step + 1) << '/' << config.totalSteps()
                          << "  active=" << particle_stats.active_particles
                          << "  total_collisions=" << collision_handler.getStatistics().total_collisions
                          << "  generated=" << collision_result.generated_particles << '\n';
            }
        }

        const auto wall_end = std::chrono::steady_clock::now();
        const double wall_time_seconds =
            std::chrono::duration<double>(wall_end - wall_start).count();

        const LayerProfile profile =
            collectLayerProfile(config, layers, cycle, *particle_manager);
        writeLayerProfileResults(config.output_csv, profile, config);

        const auto particle_stats = particle_manager->getStatistics();
        const auto collision_stats = collision_handler.getStatistics();
        printRuntimeSummary(config, *mesh, particle_stats, collision_stats, wall_time_seconds,
                            config.output_csv,
                            cycle.getStatistics().surface_deposited_charge);

        return 0;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "SCDAT_PICcore failed: " << ex.what() << '\n';
        return 1;
    }
}
