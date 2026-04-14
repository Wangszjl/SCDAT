#include "PICCycle.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <unordered_set>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace SCDAT
{
namespace PICcore
{

namespace
{
constexpr ElementId kInvalidElementId = static_cast<ElementId>(-1);

std::size_t toBoundaryIndex(PICCycle::BoundaryFace face)
{
    return static_cast<std::size_t>(face);
}

Utils::Vector3D inwardBoundaryNormal(PICCycle::BoundaryFace face)
{
    switch (face)
    {
    case PICCycle::BoundaryFace::XMin:
        return Utils::Vector3D(1.0, 0.0, 0.0);
    case PICCycle::BoundaryFace::XMax:
        return Utils::Vector3D(-1.0, 0.0, 0.0);
    case PICCycle::BoundaryFace::YMin:
        return Utils::Vector3D(0.0, 1.0, 0.0);
    case PICCycle::BoundaryFace::YMax:
        return Utils::Vector3D(0.0, -1.0, 0.0);
    case PICCycle::BoundaryFace::ZMin:
        return Utils::Vector3D(0.0, 0.0, 1.0);
    case PICCycle::BoundaryFace::ZMax:
        return Utils::Vector3D(0.0, 0.0, -1.0);
    }

    return Utils::Vector3D(0.0, 0.0, 0.0);
}

bool isElectronLike(const ParticleTypeDef& particle)
{
    switch (particle.getType())
    {
    case ::SCDAT::Particle::ParticleType::ELECTRON:
    case ::SCDAT::Particle::ParticleType::PHOTOELECTRON:
    case ::SCDAT::Particle::ParticleType::SECONDARY_ELECTRON:
    case ::SCDAT::Particle::ParticleType::BACKSCATTERED_ELECTRON:
    case ::SCDAT::Particle::ParticleType::AUGER_ELECTRON:
    case ::SCDAT::Particle::ParticleType::THERMAL_ELECTRON:
    case ::SCDAT::Particle::ParticleType::FIELD_EMISSION_ELECTRON:
    case ::SCDAT::Particle::ParticleType::BETA_PARTICLE:
        return true;
    default:
        return particle.getCharge() < 0.0;
    }
}

bool isIonLike(const ParticleTypeDef& particle)
{
    switch (particle.getType())
    {
    case ::SCDAT::Particle::ParticleType::ION:
    case ::SCDAT::Particle::ParticleType::POSITIVE_ION:
    case ::SCDAT::Particle::ParticleType::MOLECULAR_ION:
    case ::SCDAT::Particle::ParticleType::CLUSTER_ION:
    case ::SCDAT::Particle::ParticleType::PROTON:
    case ::SCDAT::Particle::ParticleType::ALPHA:
    case ::SCDAT::Particle::ParticleType::HEAVY_ION:
        return true;
    default:
        return particle.getCharge() > 0.0;
    }
}

std::unordered_map<Mesh::NodeId, std::vector<Mesh::NodeId>>
buildNodeNeighbors(const std::vector<Mesh::ElementPtr>& elements)
{
    std::unordered_map<Mesh::NodeId, std::unordered_set<Mesh::NodeId>> unique_neighbors;

    for (const auto& element : elements)
    {
        if (!element)
        {
            continue;
        }

        const auto& nodes = element->getNodes();
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            if (!nodes[i])
            {
                continue;
            }

            for (std::size_t j = 0; j < nodes.size(); ++j)
            {
                if (i == j || !nodes[j])
                {
                    continue;
                }

                unique_neighbors[nodes[i]->getId()].insert(nodes[j]->getId());
            }
        }
    }

    std::unordered_map<Mesh::NodeId, std::vector<Mesh::NodeId>> neighbors;
    neighbors.reserve(unique_neighbors.size());
    for (auto& [node_id, neighbor_set] : unique_neighbors)
    {
        neighbors[node_id] = std::vector<Mesh::NodeId>(neighbor_set.begin(), neighbor_set.end());
    }

    return neighbors;
}

bool solve3x3(std::array<std::array<double, 3>, 3> matrix, std::array<double, 3> rhs,
              std::array<double, 3>& solution)
{
    for (std::size_t pivot = 0; pivot < 3; ++pivot)
    {
        std::size_t best = pivot;
        for (std::size_t row = pivot + 1; row < 3; ++row)
        {
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot]))
            {
                best = row;
            }
        }

        if (std::abs(matrix[best][pivot]) < 1.0e-18)
        {
            return false;
        }

        if (best != pivot)
        {
            std::swap(matrix[best], matrix[pivot]);
            std::swap(rhs[best], rhs[pivot]);
        }

        for (std::size_t row = pivot + 1; row < 3; ++row)
        {
            const double factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (std::size_t col = pivot; col < 3; ++col)
            {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    for (int row = 2; row >= 0; --row)
    {
        double sum = rhs[static_cast<std::size_t>(row)];
        for (std::size_t col = static_cast<std::size_t>(row) + 1; col < 3; ++col)
        {
            sum -= matrix[static_cast<std::size_t>(row)][col] * solution[col];
        }

        solution[static_cast<std::size_t>(row)] =
            sum / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }

    return true;
}

Utils::Vector3D computeLeastSquaresGradient(const std::vector<Utils::Vector3D>& delta_pos,
                                            const std::vector<double>& delta_phi)
{
    std::array<std::array<double, 3>, 3> normal_matrix{};
    std::array<double, 3> rhs{};

    for (std::size_t i = 0; i < delta_pos.size(); ++i)
    {
        const double dx = delta_pos[i].x();
        const double dy = delta_pos[i].y();
        const double dz = delta_pos[i].z();
        const double dphi = delta_phi[i];

        normal_matrix[0][0] += dx * dx;
        normal_matrix[0][1] += dx * dy;
        normal_matrix[0][2] += dx * dz;
        normal_matrix[1][0] += dy * dx;
        normal_matrix[1][1] += dy * dy;
        normal_matrix[1][2] += dy * dz;
        normal_matrix[2][0] += dz * dx;
        normal_matrix[2][1] += dz * dy;
        normal_matrix[2][2] += dz * dz;

        rhs[0] += dx * dphi;
        rhs[1] += dy * dphi;
        rhs[2] += dz * dphi;
    }

    constexpr double kRegularization = 1.0e-18;
    normal_matrix[0][0] += kRegularization;
    normal_matrix[1][1] += kRegularization;
    normal_matrix[2][2] += kRegularization;

    std::array<double, 3> solution{0.0, 0.0, 0.0};
    if (!solve3x3(normal_matrix, rhs, solution))
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    return Utils::Vector3D(solution[0], solution[1], solution[2]);
}
} // namespace

PICCycle::PICCycle() {}

PICCycle::PICCycle(const Parameters& params) : params_(params) {}

void PICCycle::setMesh(const std::vector<Mesh::NodePtr>& nodes,
                       const std::vector<Mesh::ElementPtr>& elements)
{
    nodes_ = nodes;
    elements_ = elements;
    updateDomainBounds();
}

void PICCycle::setParticleManager(std::shared_ptr<ParticleManager> particle_manager)
{
    particle_manager_ = particle_manager;
}

void PICCycle::setPoissonSolver(std::shared_ptr<PoissonSolver> poisson_solver)
{
    field_poisson_solver_ = poisson_solver;
}

void PICCycle::setBoundaryCondition(BoundaryFace face, BoundaryConditionPtr condition)
{
    boundary_conditions_[toBoundaryIndex(face)] = std::move(condition);
}

void PICCycle::clearBoundaryConditions()
{
    boundary_conditions_.fill(nullptr);
}

bool PICCycle::initialize()
{
    if (nodes_.empty() || elements_.empty() || !particle_manager_ || !field_poisson_solver_)
    {
        std::cerr << "PICCycle initialization failed: missing required components" << std::endl;
        return false;
    }

    // Allocate field storage.
    size_t num_nodes = nodes_.size();
    charge_density_.resize(num_nodes, 0.0);
    potential_.resize(num_nodes, 0.0);
    electric_field_.resize(num_nodes, Utils::Vector3D(0.0, 0.0, 0.0));
    magnetic_field_.resize(num_nodes, Utils::Vector3D(0.0, 0.0, 0.0));

    // Create core algorithm helpers.
    cross_tetrahedron_ = std::make_shared<CrossTetrahedron>();
    field_interpolator_ = std::make_shared<FieldInterpolator>();
    charge_depositor_ = std::make_shared<ChargeDepositor>();

    // Wire helper dependencies.
    cross_tetrahedron_->setMesh(nodes_, elements_);
    cross_tetrahedron_->setFieldInterpolator(field_interpolator_);
    cross_tetrahedron_->setChargeDepositor(charge_depositor_);

    field_interpolator_->setMesh(nodes_, elements_);
    charge_depositor_->setMesh(nodes_, elements_);
    charge_depositor_->setTrajectoryDepositionScheme(
        params_.trajectory_charge_deposition_kernel,
        params_.trajectory_charge_deposition_segments);

    // Reset runtime statistics.
    cycle_stats_ = CycleStatistics();
    surface_current_ledger_.fill(SurfaceCurrentLedgerEntry{});
    for (std::size_t face_index = 0; face_index < surface_boundary_metadata_.size(); ++face_index)
    {
        if (surface_boundary_metadata_[face_index].surface_id < 0)
        {
            surface_boundary_metadata_[face_index].surface_id = static_cast<int>(face_index);
        }
        surface_current_ledger_[face_index].metadata = surface_boundary_metadata_[face_index];
    }
    field_solver_solution_valid_ = false;
    updateDomainBounds();

    if (params_.verbose)
    {
        std::cout << "PIC cycle initialized" << std::endl;
        std::cout << "Mesh nodes: " << num_nodes << std::endl;
        std::cout << "Mesh elements: " << elements_.size() << std::endl;
        std::cout << "Active particles: " << particle_manager_->getContainer().size() << std::endl;
    }

    return true;
}

bool PICCycle::executeTimeStep()
{
    auto start_time = std::chrono::high_resolution_clock::now();
    using Clock = std::chrono::high_resolution_clock;

    const bool profile_cycle = (std::getenv("PIC_PROFILE_CYCLE") != nullptr);
    static double t_clear_acc = 0.0;
    static double t_push_acc = 0.0;
    static double t_charge_acc = 0.0;
    static double t_poisson_acc = 0.0;
    static double t_field_acc = 0.0;
    static double t_stats_acc = 0.0;

    auto toSec = [](const Clock::time_point& a, const Clock::time_point& b)
    { return std::chrono::duration<double>(b - a).count(); };
    auto stage_t0 = start_time;
    auto stage_t1 = start_time;

    // 1. Clear the previous charge density.
    stage_t0 = Clock::now();
    charge_depositor_->clearChargeDensity(charge_density_);
    stage_t1 = Clock::now();
    if (profile_cycle)
    {
        t_clear_acc += toSec(stage_t0, stage_t1);
    }

    // 2. Push particles and accumulate charge.
    stage_t0 = Clock::now();
    pushParticles();
    stage_t1 = Clock::now();
    if (profile_cycle)
    {
        t_push_acc += toSec(stage_t0, stage_t1);
    }

    stage_t0 = Clock::now();
    computeChargeDensity();
    stage_t1 = Clock::now();
    if (profile_cycle)
    {
        t_charge_acc += toSec(stage_t0, stage_t1);
    }

    // 3. Solve the Poisson equation.
    stage_t0 = Clock::now();
    solvePoissonEquation();
    stage_t1 = Clock::now();
    if (profile_cycle)
    {
        t_poisson_acc += toSec(stage_t0, stage_t1);
    }

    // 4. Update the electric field.
    stage_t0 = Clock::now();
    updateElectricField();
    stage_t1 = Clock::now();
    if (profile_cycle)
    {
        t_field_acc += toSec(stage_t0, stage_t1);
    }

    // 5. Refresh cycle statistics.
    stage_t0 = Clock::now();
    updateStatistics();
    stage_t1 = Clock::now();
    if (profile_cycle)
    {
        t_stats_acc += toSec(stage_t0, stage_t1);
    }

    // 6. Check convergence status.
    bool converged = checkConvergence();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    cycle_stats_.completed_steps++;
    cycle_stats_.simulation_time += params_.time_step;
    cycle_stats_.average_step_time =
        (cycle_stats_.average_step_time * (cycle_stats_.completed_steps - 1) +
         static_cast<double>(duration.count()) / 1000.0) /
        cycle_stats_.completed_steps;

    if (profile_cycle && (cycle_stats_.completed_steps % 100 == 0))
    {
        const double total = std::max(t_clear_acc + t_push_acc + t_charge_acc + t_poisson_acc +
                                          t_field_acc + t_stats_acc,
                                      1.0e-12);
        auto pct = [total](double x) { return 100.0 * x / total; };
        std::cout << "[PIC_PROFILE] steps=" << cycle_stats_.completed_steps
                  << " clear=" << pct(t_clear_acc) << "%"
                  << " push=" << pct(t_push_acc) << "%"
                  << " charge=" << pct(t_charge_acc) << "%"
                  << " poisson=" << pct(t_poisson_acc) << "%"
                  << " field=" << pct(t_field_acc) << "%"
                  << " stats=" << pct(t_stats_acc) << "%" << std::endl;
    }

    if (params_.verbose && (cycle_stats_.completed_steps % params_.statistics_frequency == 0))
    {
        outputStatistics();
    }

    return !converged; // Return false when convergence has been reached.
}

bool PICCycle::executeSimulation(double total_time)
{
    if (!initialize())
    {
        return false;
    }

    int total_steps = static_cast<int>(total_time / params_.time_step);

    if (params_.verbose)
    {
        std::cout << "\n=== Starting PIC simulation ===" << std::endl;
        std::cout << "Total simulation time: " << total_time << " s" << std::endl;
        std::cout << "Time step: " << params_.time_step << " s" << std::endl;
        std::cout << "Total steps: " << total_steps << std::endl;
    }

    for (int step = 0; step < total_steps; ++step)
    {
        if (!executeTimeStep())
        {
            if (params_.verbose)
            {
                std::cout << "Simulation failed at step " << step << std::endl;
            }
            break;
        }

        if (cycle_stats_.completed_steps >= params_.max_iterations)
        {
            if (params_.verbose)
            {
                std::cout << "Reached the maximum iteration count, stopping simulation"
                          << std::endl;
            }
            break;
        }
    }

    if (params_.verbose)
    {
        std::cout << "\n=== Simulation complete ===" << std::endl;
        outputStatistics();
    }

    return true;
}

void PICCycle::pushParticles()
{
    if (!particle_manager_ || !cross_tetrahedron_)
        return;

    // Push the latest field state into the particle-side helpers.
    field_interpolator_->setElectricField(electric_field_);
    field_interpolator_->setMagneticField(magnetic_field_);
    field_interpolator_->setPotential(potential_);

    // Traverse all active particles.
    auto& container = particle_manager_->getContainer();
    size_t num_particles = container.size();
    auto it_begin = container.begin();
    size_t num_nodes = nodes_.size();
    const bool debug_push = (std::getenv("PIC_DEBUG_PUSH") != nullptr);

#ifdef USE_OPENMP
    int max_threads = omp_get_max_threads();
    // Thread-local charge buffers.
    std::vector<std::vector<double>> thread_charge_densities(max_threads,
                                                             std::vector<double>(num_nodes, 0.0));

    // Thread-local location updates: {index, element_id}.
    struct LocUpdate
    {
        size_t index;
        ElementId eid;
        LocUpdate(size_t i, ElementId e) : index(i), eid(e) {}
    };
    std::vector<std::vector<LocUpdate>> thread_loc_updates(max_threads);
    for (auto& vec : thread_loc_updates)
        vec.reserve(num_particles / max_threads + 100);

    std::vector<std::vector<RuntimeBoundaryEvent>> thread_boundary_events(max_threads);
    for (auto& vec : thread_boundary_events)
        vec.reserve(num_particles / max_threads + 32);

    // Debug counters for short profiling runs.
    std::atomic<int> debug_lost_count{0};
    std::atomic<int> debug_survived_count{0};

#pragma omp parallel for schedule(dynamic, 1000)
    for (int i = 0; i < static_cast<int>(num_particles); ++i)
    {
        int thread_id = omp_get_thread_num();
        auto& particle = *(it_begin + i);
        const ElementId hint_elem = findParticleElement(&particle);

        // Integrate the particle with thread-local charge accumulation.
        auto result = cross_tetrahedron_->crossTetrahedron(
            &particle, params_.time_step, thread_charge_densities[thread_id], true, hint_elem);

        // Update the particle state.
        if (!result.lost_particle)
        {
            if (debug_push)
            {
                debug_survived_count.fetch_add(1, std::memory_order_relaxed);
            }
            particle.setPosition(result.final_position);
            particle.setVelocity(result.final_velocity);
            // Defer the location-map write until after the parallel loop.
            thread_loc_updates[thread_id].emplace_back(i, result.final_element);
        }
        else
        {
            if (debug_push)
            {
                debug_lost_count.fetch_add(1, std::memory_order_relaxed);
            }

            particle.setPosition(result.final_position);
            particle.setVelocity(result.final_velocity);
            thread_boundary_events[thread_id].push_back(
                RuntimeBoundaryEvent{static_cast<std::size_t>(i), result.final_element});
        }
    }

    // Print a short debug summary for the first few steps.
    static int step_counter = 0;
    step_counter++;
    if (debug_push && step_counter <= 3)
    {
        std::cout << "DEBUG pushParticles step=" << step_counter
                  << " survived=" << debug_survived_count.load()
                  << " lost=" << debug_lost_count.load() << std::endl;
    }

// Reduce thread-local charge buffers into the global array.
#pragma omp parallel for
    for (int n = 0; n < static_cast<int>(num_nodes); ++n)
    {
        for (int t = 0; t < max_threads; ++t)
        {
            charge_density_[n] += thread_charge_densities[t][n];
        }
    }

    // Apply location-map updates serially.
    for (int t = 0; t < max_threads; ++t)
    {
        for (const auto& update : thread_loc_updates[t])
        {
            auto& particle = *(it_begin + update.index);
            updateParticleLocation(&particle, update.eid);
        }
    }

    std::vector<RuntimeBoundaryEvent> boundary_events;
    for (int t = 0; t < max_threads; ++t)
    {
        boundary_events.insert(boundary_events.end(), thread_boundary_events[t].begin(),
                               thread_boundary_events[t].end());
    }

    std::vector<ParticleTypeDef> emitted_particles;
    processRuntimeBoundaryEvents(boundary_events, emitted_particles);
    if (!emitted_particles.empty())
    {
        particle_manager_->getContainer().addParticles(emitted_particles);
    }
#else
    std::vector<RuntimeBoundaryEvent> boundary_events;
    boundary_events.reserve(num_particles / 8 + 8);

    for (auto& particle : container)
    {
        const ElementId hint_elem = findParticleElement(&particle);

        // Integrate the particle trajectory.
        auto result = cross_tetrahedron_->crossTetrahedron(&particle, params_.time_step,
                                                           charge_density_, true, hint_elem);

        // Update the particle state.
        if (!result.lost_particle)
        {
            particle.setPosition(result.final_position);
            particle.setVelocity(result.final_velocity);
            // Update the cached particle location.
            updateParticleLocation(&particle, result.final_element);
        }
        else
        {
            particle.setPosition(result.final_position);
            particle.setVelocity(result.final_velocity);
            boundary_events.push_back(RuntimeBoundaryEvent{
                static_cast<std::size_t>(&particle - &(*it_begin)), result.final_element});

            if (params_.verbose)
            {
                std::cout << "Particle lost: " << result.loss_reason << std::endl;
            }
        }
    }

    std::vector<ParticleTypeDef> emitted_particles;
    processRuntimeBoundaryEvents(boundary_events, emitted_particles);
    if (!emitted_particles.empty())
    {
        particle_manager_->getContainer().addParticles(emitted_particles);
    }
#endif
}

void PICCycle::computeChargeDensity()
{
    // Charge deposition is performed during pushParticles().
    // This hook is kept for optional post-processing such as smoothing.
    double total_charge = charge_depositor_->computeTotalCharge(charge_density_);
    cycle_stats_.total_charge = total_charge;

    if (params_.charge_conservation && params_.verbose &&
        (cycle_stats_.completed_steps % params_.statistics_frequency == 0))
    {
        std::cout << "Total charge: " << total_charge << " C" << std::endl;
    }
}

void PICCycle::solvePoissonEquation()
{
    if (nodes_.empty())
        return;

    const size_t num_nodes = nodes_.size();
    if (potential_.size() != num_nodes)
    {
        potential_.resize(num_nodes, 0.0);
    }

    const auto tryFieldSolver = [&](SCDAT::Solver::SolverType solver_type) -> bool
    {
        if (!field_poisson_solver_)
        {
            return false;
        }

        field_poisson_solver_->clearBoundaryConditions();
        field_poisson_solver_->clearChargeDensities();

        bool has_dirichlet_boundary = false;
        for (size_t i = 0; i < num_nodes; ++i)
        {
            const auto& node = nodes_[i];
            if (!node)
            {
                continue;
            }

            const double rho = (i < charge_density_.size()) ? charge_density_[i] : 0.0;
            node->setChargeDensity(rho);
            node->setPotential(potential_[i]);
            field_poisson_solver_->setChargeDensity(node->getId(), rho);

            switch (node->getBoundaryType())
            {
            case Mesh::BoundaryType::DIRICHLET:
                field_poisson_solver_->addBoundaryCondition(
                    node->getId(),
                    std::make_shared<FieldSolver::BoundaryCondition>(
                        FieldSolver::BoundaryConditionType::DIRICHLET, potential_[i]));
                has_dirichlet_boundary = true;
                break;
            case Mesh::BoundaryType::NEUMANN:
                field_poisson_solver_->addBoundaryCondition(
                    node->getId(), std::make_shared<FieldSolver::BoundaryCondition>(
                                       FieldSolver::BoundaryConditionType::NEUMANN, 0.0));
                break;
            default:
                break;
            }
        }

        if (!has_dirichlet_boundary || !field_poisson_solver_->solve(solver_type))
        {
            return false;
        }

        bool all_values_finite = true;
        for (size_t i = 0; i < num_nodes; ++i)
        {
            const auto& node = nodes_[i];
            if (!node)
            {
                continue;
            }

            potential_[i] = node->getPotential();
            electric_field_[i] = node->getElectricField();

            all_values_finite = all_values_finite && std::isfinite(potential_[i]) &&
                                std::isfinite(electric_field_[i].x()) &&
                                std::isfinite(electric_field_[i].y()) &&
                                std::isfinite(electric_field_[i].z());
        }

        if (all_values_finite)
        {
            field_solver_solution_valid_ = true;
            return true;
        }

        return false;
    };

    if (tryFieldSolver(SCDAT::Solver::SolverType::DIRECT))
    {
        return;
    }

    static bool warned_iterative_retry = false;
    if (!warned_iterative_retry)
    {
        std::cerr << "PICCycle: direct Poisson solve failed, retrying with iterative solver."
                  << std::endl;
        warned_iterative_retry = true;
    }

    if (tryFieldSolver(SCDAT::Solver::SolverType::ITERATIVE))
    {
        return;
    }

    field_solver_solution_valid_ = false;

    static bool warned_solver_failure = false;
    if (!warned_solver_failure)
    {
        std::cerr << "PICCycle: generic Poisson solver failed; retaining the previous potential "
                     "field."
                  << std::endl;
        warned_solver_failure = true;
    }

    for (size_t i = 0; i < num_nodes; ++i)
    {
        if (nodes_[i])
        {
            nodes_[i]->setPotential(potential_[i]);
        }
    }
}

void PICCycle::updateElectricField()
{
    if (nodes_.empty() || electric_field_.empty())
        return;

    // Emit extra diagnostics only when field debugging is explicitly enabled.
    static int update_count = 0;
    update_count++;
    const bool debug_field = (std::getenv("PIC_DEBUG_FIELD") != nullptr);
    const bool debug_first_call = debug_field && (update_count <= 2);

    if (field_poisson_solver_ && field_solver_solution_valid_ &&
        nodes_.size() == electric_field_.size())
    {
        for (size_t i = 0; i < nodes_.size(); ++i)
        {
            if (!nodes_[i])
            {
                electric_field_[i] = Utils::Vector3D(0.0, 0.0, 0.0);
                continue;
            }

            potential_[i] = nodes_[i]->getPotential();
            electric_field_[i] = nodes_[i]->getElectricField();
        }

        if (debug_first_call)
        {
            std::cout << "DEBUG updateElectricField: using FieldSolver node gradients" << std::endl;
        }
        return;
    }

    if (potential_.size() != nodes_.size())
    {
        std::fill(electric_field_.begin(), electric_field_.end(), Utils::Vector3D(0.0, 0.0, 0.0));
        if (debug_first_call)
        {
            std::cout << "DEBUG updateElectricField: no valid potential profile available"
                      << std::endl;
        }
        return;
    }

    const auto neighbors = buildNodeNeighbors(elements_);
    bool all_values_finite = true;
    std::size_t reconstructed_nodes = 0;

    for (size_t i = 0; i < nodes_.size(); ++i)
    {
        const auto& node = nodes_[i];
        if (!node)
        {
            electric_field_[i] = Utils::Vector3D(0.0, 0.0, 0.0);
            continue;
        }

        std::vector<Utils::Vector3D> delta_pos;
        std::vector<double> delta_phi;
        if (const auto it = neighbors.find(node->getId()); it != neighbors.end())
        {
            delta_pos.reserve(it->second.size());
            delta_phi.reserve(it->second.size());

            for (Mesh::NodeId neighbor_id : it->second)
            {
                if (neighbor_id >= nodes_.size())
                {
                    continue;
                }

                const auto& neighbor = nodes_[neighbor_id];
                if (!neighbor)
                {
                    continue;
                }

                const auto displacement = neighbor->getPosition() - node->getPosition();
                if (displacement.magnitudeSquared() <= 1.0e-24)
                {
                    continue;
                }

                delta_pos.push_back(displacement);
                delta_phi.push_back(potential_[neighbor_id] - potential_[i]);
            }
        }

        Utils::Vector3D gradient(0.0, 0.0, 0.0);
        if (!delta_pos.empty())
        {
            gradient = computeLeastSquaresGradient(delta_pos, delta_phi);
            ++reconstructed_nodes;
        }

        electric_field_[i] = Utils::Vector3D(-gradient.x(), -gradient.y(), -gradient.z());
        node->setPotential(potential_[i]);
        node->setElectricField(electric_field_[i]);

        all_values_finite = all_values_finite && std::isfinite(electric_field_[i].x()) &&
                            std::isfinite(electric_field_[i].y()) &&
                            std::isfinite(electric_field_[i].z());
    }

    if (debug_first_call)
    {
        const auto [min_it, max_it] = std::minmax_element(potential_.begin(), potential_.end());

        std::cout << "DEBUG updateElectricField call #" << update_count
                  << ": using reconstructed nodal gradients for " << reconstructed_nodes << " nodes"
                  << std::endl;
        std::cout << "  Potential range: [" << *min_it << ", " << *max_it << "] V" << std::endl;
    }

    if (!all_values_finite)
    {
        std::fill(electric_field_.begin(), electric_field_.end(), Utils::Vector3D(0.0, 0.0, 0.0));
        for (size_t i = 0; i < nodes_.size(); ++i)
        {
            if (nodes_[i])
            {
                nodes_[i]->setElectricField(electric_field_[i]);
            }
        }
    }
}

Utils::Vector3D PICCycle::computeGradient(size_t node_id)
{
    if (node_id >= nodes_.size() || potential_.size() != nodes_.size() || !nodes_[node_id])
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    const auto neighbors = buildNodeNeighbors(elements_);
    const auto it = neighbors.find(nodes_[node_id]->getId());
    if (it == neighbors.end())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    std::vector<Utils::Vector3D> delta_pos;
    std::vector<double> delta_phi;
    delta_pos.reserve(it->second.size());
    delta_phi.reserve(it->second.size());

    for (Mesh::NodeId neighbor_id : it->second)
    {
        if (neighbor_id >= nodes_.size() || !nodes_[neighbor_id])
        {
            continue;
        }

        const auto displacement =
            nodes_[neighbor_id]->getPosition() - nodes_[node_id]->getPosition();
        if (displacement.magnitudeSquared() <= 1.0e-24)
        {
            continue;
        }

        delta_pos.push_back(displacement);
        delta_phi.push_back(potential_[neighbor_id] - potential_[node_id]);
    }

    if (delta_pos.empty())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    return computeLeastSquaresGradient(delta_pos, delta_phi);
}

void PICCycle::updateStatistics()
{
    if (!particle_manager_)
        return;

    cycle_stats_.active_particles = static_cast<int>(particle_manager_->getContainer().size());

    // Compute the total kinetic energy.
    double total_energy = 0.0;
    const auto& container = particle_manager_->getContainer();

    for (const auto& particle : container)
    {
        double mass = particle.getMass();
        double vel_sq = particle.getVelocity().magnitudeSquared();

        // Skip invalid floating-point values.
        if (!std::isnan(mass) && !std::isinf(mass) && !std::isnan(vel_sq) && !std::isinf(vel_sq))
        {
            double kinetic_energy = 0.5 * mass * vel_sq;
            total_energy += kinetic_energy;
        }
    }

    cycle_stats_.total_energy = total_energy;
}

bool PICCycle::checkConvergence()
{
    // For the Turner benchmark we intentionally avoid early convergence termination.
    return false;
}

void PICCycle::outputStatistics()
{
    std::cout << "\n=== PIC cycle statistics (step " << cycle_stats_.completed_steps
              << ") ===" << std::endl;
    std::cout << "Simulation time: " << std::scientific << cycle_stats_.simulation_time << " s"
              << std::endl;
    std::cout << "Active particles: " << cycle_stats_.active_particles << std::endl;
    std::cout << "Total charge: " << std::scientific << cycle_stats_.total_charge << " C"
              << std::endl;
    std::cout << "Total energy: " << std::scientific << cycle_stats_.total_energy << " J"
              << std::endl;
    std::cout << "Surface deposited charge: " << std::scientific
              << cycle_stats_.surface_deposited_charge << " C" << std::endl;
    std::cout << "Average step time: " << std::fixed << std::setprecision(2)
              << cycle_stats_.average_step_time << " ms" << std::endl;

    // Print particle-tracing statistics.
    const auto& cross_stats = cross_tetrahedron_->getStatistics();
    std::cout << "Trajectory crossings: " << cross_stats.total_crossings << std::endl;
    std::cout << "Boundary crossings: " << cross_stats.boundary_crossings << std::endl;
    std::cout << "Lost particles: " << cross_stats.lost_particles << std::endl;
    std::cout << "Reflection events: " << cross_stats.reflection_events << std::endl;
    std::cout << "==============================\n" << std::endl;
}

void PICCycle::addParticleLocation(ParticlePtr particle, ElementId element_id)
{
    if (!particle)
        return;

    ParticleLocation location(particle, element_id);
    particle_locations_[particle->getId()] = location;
}

void PICCycle::updateParticleLocation(ParticlePtr particle, ElementId new_element_id)
{
    if (!particle)
        return;

    auto it = particle_locations_.find(particle->getId());
    if (it != particle_locations_.end())
    {
        it->second.current_element = new_element_id;
        it->second.is_valid = true;
    }
    else
    {
        addParticleLocation(particle, new_element_id);
    }
}

ElementId PICCycle::findParticleElement(ParticlePtr particle) const
{
    if (!particle)
        return kInvalidElementId;

    auto it = particle_locations_.find(particle->getId());
    if (it != particle_locations_.end() && it->second.is_valid)
    {
        return it->second.current_element;
    }

    return kInvalidElementId;
}

void PICCycle::updateDomainBounds()
{
    domain_bounds_.valid = false;
    if (nodes_.empty())
    {
        return;
    }

    domain_bounds_.min_corner = nodes_.front()->getPosition();
    domain_bounds_.max_corner = nodes_.front()->getPosition();
    for (const auto& node : nodes_)
    {
        if (!node)
        {
            continue;
        }

        const auto& position = node->getPosition();
        domain_bounds_.min_corner.setX(std::min(domain_bounds_.min_corner.x(), position.x()));
        domain_bounds_.min_corner.setY(std::min(domain_bounds_.min_corner.y(), position.y()));
        domain_bounds_.min_corner.setZ(std::min(domain_bounds_.min_corner.z(), position.z()));
        domain_bounds_.max_corner.setX(std::max(domain_bounds_.max_corner.x(), position.x()));
        domain_bounds_.max_corner.setY(std::max(domain_bounds_.max_corner.y(), position.y()));
        domain_bounds_.max_corner.setZ(std::max(domain_bounds_.max_corner.z(), position.z()));
    }

    domain_bounds_.valid = true;
}

bool PICCycle::hasRuntimeBoundaryConditions() const
{
    return std::any_of(boundary_conditions_.begin(), boundary_conditions_.end(),
                       [](const auto& condition) { return static_cast<bool>(condition); });
}

std::optional<PICCycle::BoundaryFace> PICCycle::detectBoundaryFace(const Utils::Point3D& position,
                                                                   double outside_tolerance,
                                                                   bool allow_near_match,
                                                                   double near_tolerance) const
{
    if (!domain_bounds_.valid)
    {
        return std::nullopt;
    }

    if (position.x() < domain_bounds_.min_corner.x() - outside_tolerance)
    {
        return BoundaryFace::XMin;
    }
    if (position.x() > domain_bounds_.max_corner.x() + outside_tolerance)
    {
        return BoundaryFace::XMax;
    }
    if (position.y() < domain_bounds_.min_corner.y() - outside_tolerance)
    {
        return BoundaryFace::YMin;
    }
    if (position.y() > domain_bounds_.max_corner.y() + outside_tolerance)
    {
        return BoundaryFace::YMax;
    }
    if (position.z() < domain_bounds_.min_corner.z() - outside_tolerance)
    {
        return BoundaryFace::ZMin;
    }
    if (position.z() > domain_bounds_.max_corner.z() + outside_tolerance)
    {
        return BoundaryFace::ZMax;
    }

    if (!allow_near_match)
    {
        return std::nullopt;
    }

    std::optional<BoundaryFace> nearest_face;
    double nearest_distance = std::numeric_limits<double>::max();
    const auto considerFace = [&](BoundaryFace face, double distance)
    {
        if (distance <= near_tolerance && distance < nearest_distance)
        {
            nearest_face = face;
            nearest_distance = distance;
        }
    };

    considerFace(BoundaryFace::XMin, std::abs(position.x() - domain_bounds_.min_corner.x()));
    considerFace(BoundaryFace::XMax, std::abs(position.x() - domain_bounds_.max_corner.x()));
    considerFace(BoundaryFace::YMin, std::abs(position.y() - domain_bounds_.min_corner.y()));
    considerFace(BoundaryFace::YMax, std::abs(position.y() - domain_bounds_.max_corner.y()));
    considerFace(BoundaryFace::ZMin, std::abs(position.z() - domain_bounds_.min_corner.z()));
    considerFace(BoundaryFace::ZMax, std::abs(position.z() - domain_bounds_.max_corner.z()));

    return nearest_face;
}

Utils::Point3D PICCycle::projectToBoundaryFace(BoundaryFace face, const Utils::Point3D& position,
                                               double inward_offset) const
{
    Utils::Point3D projected = position;
    projected.setX(
        std::clamp(projected.x(), domain_bounds_.min_corner.x(), domain_bounds_.max_corner.x()));
    projected.setY(
        std::clamp(projected.y(), domain_bounds_.min_corner.y(), domain_bounds_.max_corner.y()));
    projected.setZ(
        std::clamp(projected.z(), domain_bounds_.min_corner.z(), domain_bounds_.max_corner.z()));

    switch (face)
    {
    case BoundaryFace::XMin:
        projected.setX(domain_bounds_.min_corner.x());
        break;
    case BoundaryFace::XMax:
        projected.setX(domain_bounds_.max_corner.x());
        break;
    case BoundaryFace::YMin:
        projected.setY(domain_bounds_.min_corner.y());
        break;
    case BoundaryFace::YMax:
        projected.setY(domain_bounds_.max_corner.y());
        break;
    case BoundaryFace::ZMin:
        projected.setZ(domain_bounds_.min_corner.z());
        break;
    case BoundaryFace::ZMax:
        projected.setZ(domain_bounds_.max_corner.z());
        break;
    }

    projected += inwardBoundaryNormal(face) * inward_offset;
    projected.setX(
        std::clamp(projected.x(), domain_bounds_.min_corner.x(), domain_bounds_.max_corner.x()));
    projected.setY(
        std::clamp(projected.y(), domain_bounds_.min_corner.y(), domain_bounds_.max_corner.y()));
    projected.setZ(
        std::clamp(projected.z(), domain_bounds_.min_corner.z(), domain_bounds_.max_corner.z()));
    return projected;
}

ElementId PICCycle::locateParticleElement(const Utils::Point3D& position,
                                          ElementId hint_element) const
{
    if (!field_interpolator_ || elements_.empty())
    {
        return kInvalidElementId;
    }

    std::vector<double> barycentric(4, 0.0);
    const auto canLocate = [&](ElementId element_id) -> bool
    {
        return element_id != kInvalidElementId && element_id < elements_.size() &&
               field_interpolator_->computeBarycentricCoordinates(position, element_id,
                                                                  barycentric);
    };

    if (canLocate(hint_element))
    {
        return hint_element;
    }

    if (hint_element != kInvalidElementId && hint_element < elements_.size())
    {
        for (const auto& neighbor : elements_[hint_element]->getNeighbors())
        {
            if (neighbor && canLocate(neighbor->getId()))
            {
                return neighbor->getId();
            }
        }
    }

    for (ElementId element_id = 0; element_id < elements_.size(); ++element_id)
    {
        if (canLocate(element_id))
        {
            return element_id;
        }
    }

    return kInvalidElementId;
}

bool PICCycle::applyRuntimeBoundaryCondition(ParticlePtr particle, ElementId hint_element,
                                             std::vector<ParticleTypeDef>& emitted_particles,
                                             ElementId& resolved_element)
{
    resolved_element = kInvalidElementId;
    if (!particle || !domain_bounds_.valid)
    {
        return false;
    }

    const double span_x =
        std::max(domain_bounds_.max_corner.x() - domain_bounds_.min_corner.x(), 1.0e-12);
    const double span_y =
        std::max(domain_bounds_.max_corner.y() - domain_bounds_.min_corner.y(), 1.0e-12);
    const double span_z =
        std::max(domain_bounds_.max_corner.z() - domain_bounds_.min_corner.z(), 1.0e-12);
    const double max_span = std::max({span_x, span_y, span_z});
    const double outside_tol = std::max(1.0e-10, 1.0e-8 * max_span);
    const double near_tol = std::max(5.0 * outside_tol, 1.0e-9 * max_span);

    for (int pass = 0; pass < 6; ++pass)
    {
        if (!particle->isActive())
        {
            return true;
        }

        const auto face =
            detectBoundaryFace(particle->getPosition(), outside_tol, pass == 0, near_tol);
        if (!face)
        {
            resolved_element = locateParticleElement(particle->getPosition(), hint_element);
            return resolved_element != kInvalidElementId;
        }

        const auto& condition = boundary_conditions_[toBoundaryIndex(*face)];
        if (!condition)
        {
            return false;
        }

        auto generated_particles = condition->processParticle(
            *particle, projectToBoundaryFace(*face, particle->getPosition()),
            inwardBoundaryNormal(*face), params_.time_step);
        emitted_particles.insert(emitted_particles.end(), generated_particles.begin(),
                                 generated_particles.end());
        recordSurfaceBoundaryInteraction(*face, *particle, generated_particles);

        if (particle->getStatus() == ::SCDAT::Particle::ParticleStatus::ABSORBED)
        {
            cycle_stats_.surface_deposited_charge += particle->getCharge() * particle->getWeight();
            return true;
        }
        if (particle->getStatus() == ::SCDAT::Particle::ParticleStatus::ESCAPED)
        {
            return true;
        }

        hint_element = locateParticleElement(particle->getPosition(), hint_element);
        if (hint_element != kInvalidElementId)
        {
            resolved_element = hint_element;
        }
    }

    if (particle->isActive())
    {
        resolved_element = locateParticleElement(particle->getPosition(), hint_element);
        return resolved_element != kInvalidElementId;
    }

    return true;
}

void PICCycle::recordSurfaceBoundaryInteraction(BoundaryFace face, const ParticleTypeDef& particle,
                                                const std::vector<ParticleTypeDef>& emitted_particles)
{
    auto& ledger = surface_current_ledger_[toBoundaryIndex(face)];
    ledger.metadata = surface_boundary_metadata_[toBoundaryIndex(face)];

    const double particle_charge = particle.getCharge() * particle.getWeight();
    const double particle_energy =
        particle.getKineticEnergy() * std::max(1.0, particle.getWeight());

    if (particle.getStatus() == ::SCDAT::Particle::ParticleStatus::ABSORBED)
    {
        ledger.incident_energy_j += particle_energy;
        if (isElectronLike(particle))
        {
            ledger.absorbed_electron_charge_c += particle_charge;
            ledger.absorbed_electron_particles += 1;
        }
        else if (isIonLike(particle))
        {
            ledger.absorbed_ion_charge_c += particle_charge;
            ledger.absorbed_ion_particles += 1;
        }
    }

    for (const auto& emitted_particle : emitted_particles)
    {
        if (isElectronLike(emitted_particle))
        {
            ledger.emitted_electron_charge_c +=
                emitted_particle.getCharge() * emitted_particle.getWeight();
            ledger.emitted_electron_particles += 1;
        }
    }
}

void PICCycle::processRuntimeBoundaryEvents(const std::vector<RuntimeBoundaryEvent>& events,
                                            std::vector<ParticleTypeDef>& emitted_particles)
{
    if (!particle_manager_ || events.empty())
    {
        return;
    }

    auto& container = particle_manager_->getContainer();
    auto it_begin = container.begin();
    for (const auto& event : events)
    {
        if (event.particle_index >= container.size())
        {
            continue;
        }

        auto& particle = *(it_begin + event.particle_index);
        if (!particle.isActive())
        {
            particle_locations_.erase(particle.getId());
            continue;
        }

        ElementId resolved_element = kInvalidElementId;
        if (!applyRuntimeBoundaryCondition(&particle, event.hint_element, emitted_particles,
                                           resolved_element))
        {
            particle.markAsEscaped();
        }

        if (particle.isActive())
        {
            if (resolved_element != kInvalidElementId)
            {
                updateParticleLocation(&particle, resolved_element);
            }
            else
            {
                particle.markAsEscaped();
                particle_locations_.erase(particle.getId());
            }
        }
        else
        {
            particle_locations_.erase(particle.getId());
        }
    }
}

} // namespace PICcore
} // namespace SCDAT
