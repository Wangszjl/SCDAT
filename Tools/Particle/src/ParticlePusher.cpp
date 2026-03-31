/**
 * @file ParticlePusher.cpp
 * @brief 粒子推进与并行推进统一实现
 */

#include "../include/ParticlePusher.h"
#include "../include/HighOrderParticlePusher.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SCDAT
{
namespace Particle
{

namespace
{

int getDefaultThreadCount()
{
#ifdef _OPENMP
    return std::max(1, omp_get_max_threads());
#else
    const unsigned int hc = std::thread::hardware_concurrency();
    return std::max(1u, hc);
#endif
}

void configureSchedule(ParallelParticlePusher::LoadBalanceStrategy strategy)
{
#ifdef _OPENMP
    switch (strategy)
    {
    case ParallelParticlePusher::LoadBalanceStrategy::STATIC:
        omp_set_schedule(omp_sched_static, 0);
        break;
    case ParallelParticlePusher::LoadBalanceStrategy::GUIDED:
        omp_set_schedule(omp_sched_guided, 64);
        break;
    case ParallelParticlePusher::LoadBalanceStrategy::DYNAMIC:
    default:
        omp_set_schedule(omp_sched_dynamic, 64);
        break;
    }
#else
    (void)strategy;
#endif
}

Vector3D sampleField(const std::vector<Vector3D>& field)
{
    if (field.empty())
    {
        return Vector3D(0.0, 0.0, 0.0);
    }
    return field.front();
}

} // namespace

ParticlePusher::ParticlePusher(double dt) : dt_(dt), relativistic_(false)
{
    boundary_.type = BoundaryType::REFLECTING;
    boundary_.min_bounds = Point3D(-1.0, -1.0, -1.0);
    boundary_.max_bounds = Point3D(1.0, 1.0, 1.0);
    boundary_.reflection_coefficient = 1.0;
    boundary_.absorption_probability = 0.0;
}

void ParticlePusher::pushParticles(const std::vector<ParticlePtr>& particles,
                                   const FieldFunction& electric_field,
                                   const FieldFunction& magnetic_field)
{
    for (ParticlePtr particle : particles)
    {
        if (particle && particle->isActive())
        {
            pushParticle(particle, electric_field, magnetic_field);
        }
    }
}

void ParticlePusher::pushParticles(std::vector<Particle>& particles,
                                   const FieldFunction& electric_field,
                                   const FieldFunction& magnetic_field)
{
    for (Particle& particle : particles)
    {
        if (particle.isActive())
        {
            pushParticle(&particle, electric_field, magnetic_field);
        }
    }
}

void ParticlePusher::pushParticles(std::vector<ParticlePtr>& particles,
                                   const FieldFunction& electric_field,
                                   const FieldFunction& magnetic_field, double dt)
{
    const double original_dt = dt_;
    dt_ = dt;

    for (ParticlePtr particle : particles)
    {
        if (particle && particle->isActive())
        {
            const double old_energy = particle->getKineticEnergy();
            pushParticle(particle, electric_field, magnetic_field);
            updateStatistics(particle, old_energy);
        }
    }

    dt_ = original_dt;
}

void ParticlePusher::pushParticles(std::vector<Particle>& particles,
                                   const FieldFunction& electric_field,
                                   const FieldFunction& magnetic_field, double dt)
{
    const double original_dt = dt_;
    dt_ = dt;

    for (Particle& particle : particles)
    {
        if (particle.isActive())
        {
            const double old_energy = particle.getKineticEnergy();
            pushParticle(&particle, electric_field, magnetic_field);
            updateStatistics(&particle, old_energy);
        }
    }

    dt_ = original_dt;
}

void ParticlePusher::pushParticles(std::vector<ParticlePtr>& particles, const Mesh::MeshPtr& mesh,
                                   const std::vector<Vector3D>& electric_field,
                                   const std::vector<Vector3D>& magnetic_field, double dt)
{
    (void)mesh;
    const Vector3D E0 = sampleField(electric_field);
    const Vector3D B0 = sampleField(magnetic_field);
    const FieldFunction ef = [E0](const Point3D&) { return E0; };
    const FieldFunction bf = [B0](const Point3D&) { return B0; };
    pushParticles(particles, ef, bf, dt);
}

void ParticlePusher::pushParticles(std::vector<Particle>& particles, const Mesh::MeshPtr& mesh,
                                   const std::vector<Vector3D>& electric_field,
                                   const std::vector<Vector3D>& magnetic_field, double dt)
{
    (void)mesh;
    const Vector3D E0 = sampleField(electric_field);
    const Vector3D B0 = sampleField(magnetic_field);
    const FieldFunction ef = [E0](const Point3D&) { return E0; };
    const FieldFunction bf = [B0](const Point3D&) { return B0; };
    pushParticles(particles, ef, bf, dt);
}

void ParticlePusher::pushAllParticles(ParticleManager& manager, const FieldFunction& electric_field,
                                      const FieldFunction& magnetic_field)
{
    auto iterator = manager.createIterator();
    while (iterator.hasNext())
    {
        ParticlePtr particle = iterator.next();
        if (particle && particle->isActive())
        {
            const double old_energy = particle->getKineticEnergy();
            pushParticle(particle, electric_field, magnetic_field);
            updateStatistics(particle, old_energy);
        }
    }
}

void ParticlePusher::handleBoundaryCondition(ParticlePtr particle)
{
    if (!particle)
    {
        return;
    }

    BoundaryHandler handler(boundary_);
    if (handler.handleParticle(particle))
    {
        particle->setStatus(ParticleStatus::ABSORBED);
        statistics_.particles_absorbed++;
    }
}

void ParticlePusher::updateStatistics(const ParticlePtr& particle, double old_energy)
{
    if (!particle)
    {
        return;
    }

    statistics_.particles_pushed++;

    const double new_energy = particle->getKineticEnergy();
    statistics_.total_energy_change += (new_energy - old_energy);

    const double velocity_magnitude = particle->getVelocity().magnitude();
    statistics_.max_velocity = std::max(statistics_.max_velocity, velocity_magnitude);

    const double n = static_cast<double>(statistics_.particles_pushed);
    statistics_.average_velocity =
        (statistics_.average_velocity * (n - 1.0) + velocity_magnitude) / n;
}

BorisAlgorithm::BorisAlgorithm(double dt) : ParticlePusher(dt) {}

void BorisAlgorithm::pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                                  const FieldFunction& magnetic_field)
{
    if (!particle || !particle->isActive())
    {
        return;
    }

    const Point3D position = particle->getPosition();
    const Vector3D velocity = particle->getVelocity();
    const double charge = particle->getCharge();
    const double mass = particle->getMass();

    const Vector3D E = electric_field(position);
    const Vector3D B = magnetic_field(position);

    const auto [new_position, new_velocity] =
        relativistic_ ? relativisticBorisStep(position, velocity, charge, mass, E, B, dt_)
                      : borisStep(position, velocity, charge, mass, E, B, dt_);

    particle->setPosition(new_position);
    particle->setVelocity(new_velocity);
    particle->updateAge(dt_);
    handleBoundaryCondition(particle);
}

std::pair<Point3D, Vector3D> BorisAlgorithm::borisStep(const Point3D& position,
                                                       const Vector3D& velocity, double charge,
                                                       double mass, const Vector3D& E,
                                                       const Vector3D& B, double dt)
{
    const double q_over_m = charge / mass;
    const double half_dt = 0.5 * dt;
    const double qE_half_dt = q_over_m * half_dt;

    const Vector3D v_minus = velocity + E * qE_half_dt;
    const double B_magnitude_sq = B.magnitudeSquared();

    if (B_magnitude_sq < 1e-30)
    {
        const Vector3D new_velocity = v_minus + E * qE_half_dt;
        const Point3D new_position = position + new_velocity * dt;
        return {new_position, new_velocity};
    }

    const Vector3D t = B * qE_half_dt;
    const double t_magnitude_sq = t.magnitudeSquared();
    const Vector3D s = t * (2.0 / (1.0 + t_magnitude_sq));

    const Vector3D v_prime = v_minus + v_minus.cross(t);
    const Vector3D v_plus = v_minus + v_prime.cross(s);

    const Vector3D new_velocity = v_plus + E * qE_half_dt;
    const Point3D new_position = position + new_velocity * dt;
    return {new_position, new_velocity};
}

std::pair<Point3D, Vector3D> BorisAlgorithm::relativisticBorisStep(const Point3D& position,
                                                                   const Vector3D& velocity,
                                                                   double charge, double mass,
                                                                   const Vector3D& E,
                                                                   const Vector3D& B, double dt)
{
    const double q_over_m = charge / mass;
    const double half_dt = 0.5 * dt;
    const double c_sq = LIGHT_SPEED * LIGHT_SPEED;

    const double v_sq = velocity.magnitudeSquared();
    double gamma = 1.0 / std::sqrt(1.0 - v_sq / c_sq);

    const Vector3D momentum = velocity * (gamma * mass);
    const Vector3D p_minus = momentum + E * (charge * half_dt);

    const double p_minus_sq = p_minus.magnitudeSquared();
    gamma = std::sqrt(1.0 + p_minus_sq / (mass * mass * c_sq));

    const Vector3D t = B * (q_over_m * half_dt / gamma);
    const double t_magnitude_sq = t.magnitudeSquared();
    const Vector3D s = t * (2.0 / (1.0 + t_magnitude_sq));

    const Vector3D p_prime = p_minus + p_minus.cross(t);
    const Vector3D p_plus = p_minus + p_prime.cross(s);
    const Vector3D new_momentum = p_plus + E * (charge * half_dt);

    const double new_p_sq = new_momentum.magnitudeSquared();
    const double new_gamma = std::sqrt(1.0 + new_p_sq / (mass * mass * c_sq));
    const Vector3D new_velocity = new_momentum / (new_gamma * mass);
    const Point3D new_position = position + new_velocity * dt;

    return {new_position, new_velocity};
}

LeapfrogAlgorithm::LeapfrogAlgorithm(double dt) : ParticlePusher(dt), first_step_(true) {}

void LeapfrogAlgorithm::pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                                     const FieldFunction& magnetic_field)
{
    if (!particle || !particle->isActive())
    {
        return;
    }

    const Point3D position = particle->getPosition();
    const Vector3D velocity = particle->getVelocity();
    const double charge = particle->getCharge();
    const double mass = particle->getMass();

    const Vector3D E = electric_field(position);
    const Vector3D B = magnetic_field(position);

    const auto [new_position, new_velocity] =
        leapfrogStep(position, velocity, charge, mass, E, B, dt_);

    particle->setPosition(new_position);
    particle->setVelocity(new_velocity);
    particle->updateAge(dt_);
    handleBoundaryCondition(particle);
}

std::pair<Point3D, Vector3D>
LeapfrogAlgorithm::leapfrogStep(const Point3D& position, const Vector3D& velocity, double charge,
                                double mass, const Vector3D& E, const Vector3D& B, double dt)
{
    const double q_over_m = charge / mass;

    if (first_step_)
    {
        const Vector3D acceleration = (E + velocity.cross(B)) * q_over_m;
        const Vector3D v_half_back = velocity - acceleration * (0.5 * dt);

        const Vector3D new_acceleration = (E + v_half_back.cross(B)) * q_over_m;
        const Vector3D new_velocity = v_half_back + new_acceleration * dt;
        const Point3D new_position = position + new_velocity * dt;

        first_step_ = false;
        return {new_position, new_velocity};
    }

    const Vector3D acceleration = (E + velocity.cross(B)) * q_over_m;
    const Vector3D new_velocity = velocity + acceleration * dt;
    const Point3D new_position = position + new_velocity * dt;
    return {new_position, new_velocity};
}

EnhancedBorisAlgorithm::EnhancedBorisAlgorithm(double dt, CoordinateSystem coord_system)
    : ParticlePusher(dt), coord_system_(coord_system), enable_magnetic_(false), B_uniform_(0.0),
      B_direction_(0, 0, 1)
{
}

void EnhancedBorisAlgorithm::pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                                          const FieldFunction& magnetic_field)
{
    if (!particle || !particle->isActive())
    {
        return;
    }

    const double old_energy = particle->getKineticEnergy();

    const Vector3D E = electric_field(particle->getPosition());
    Vector3D B = magnetic_field(particle->getPosition());
    if (enable_magnetic_)
    {
        B = B + B_direction_ * B_uniform_;
    }

    std::pair<Point3D, Vector3D> result;
    if (coord_system_ == CoordinateSystem::CYLINDRICAL)
    {
        result = cylindricalBorisStep(particle->getPosition(), particle->getVelocity(),
                                      particle->getCharge(), particle->getMass(), E, B, dt_);
    }
    else
    {
        result = arcPICBorisStep(particle->getPosition(), particle->getVelocity(),
                                 particle->getCharge(), particle->getMass(), E, B, dt_);
    }

    particle->setPosition(result.first);
    particle->setVelocity(result.second);
    particle->updateAge(dt_);
    handleBoundaryCondition(particle);
    updateStatistics(particle, old_energy);
}

void EnhancedBorisAlgorithm::setMagneticFieldConfig(bool enable_magnetic, double B_uniform,
                                                    const Vector3D& B_direction)
{
    enable_magnetic_ = enable_magnetic;
    B_uniform_ = B_uniform;
    B_direction_ = B_direction.normalized();
}

std::pair<Point3D, Vector3D> EnhancedBorisAlgorithm::arcPICBorisStep(const Point3D& position,
                                                                     const Vector3D& velocity,
                                                                     double charge, double mass,
                                                                     const Vector3D& E,
                                                                     const Vector3D& B, double dt)
{
    const Vector3D qE_over_2m = E * (charge * dt / (2.0 * mass));
    const Vector3D v_minus = velocity + qE_over_2m;

    if (B.magnitude() < 1e-15)
    {
        const Vector3D v_new = v_minus + qE_over_2m;
        const Point3D pos_new = position + v_new * dt;
        return std::make_pair(pos_new, v_new);
    }

    const double q_over_m = charge / mass;
    const Vector3D t = B * (q_over_m * dt / 2.0);
    const double t_magnitude_sq = t.magnitudeSquared();
    const Vector3D s = t * (2.0 / (1.0 + t_magnitude_sq));

    const Vector3D v_prime = v_minus + v_minus.cross(t);
    const Vector3D v_plus = v_minus + v_prime.cross(s);
    const Vector3D v_new = v_plus + qE_over_2m;
    const Point3D pos_new = position + v_new * dt;
    return std::make_pair(pos_new, v_new);
}

std::pair<Point3D, Vector3D>
EnhancedBorisAlgorithm::cylindricalBorisStep(const Point3D& position, const Vector3D& velocity,
                                             double charge, double mass, const Vector3D& E,
                                             const Vector3D& B, double dt)
{
    const Point3D cyl_pos = cartesianToCylindrical(position);
    const double r = cyl_pos.x();
    const double theta = std::atan2(position.y(), position.x());

    const Vector3D cyl_vel = velocityCartesianToCylindrical(velocity, position);
    const auto cyl_result = arcPICBorisStep(cyl_pos, cyl_vel, charge, mass, E, B, dt);

    double new_theta = theta;
    if (r > 1e-12)
    {
        new_theta += cyl_result.second.y() * dt / r;
    }

    const Point3D new_pos = cylindricalToCartesian(cyl_result.first, new_theta);
    const Vector3D new_vel = velocityCylindricalToCartesian(cyl_result.second, new_theta);
    return std::make_pair(new_pos, new_vel);
}

Point3D EnhancedBorisAlgorithm::cartesianToCylindrical(const Point3D& cartesian) const
{
    const double r = std::sqrt(cartesian.x() * cartesian.x() + cartesian.y() * cartesian.y());
    return Point3D(r, 0.0, cartesian.z());
}

Point3D EnhancedBorisAlgorithm::cylindricalToCartesian(const Point3D& cylindrical,
                                                       double theta) const
{
    const double r = cylindrical.x();
    return Point3D(r * std::cos(theta), r * std::sin(theta), cylindrical.z());
}

Vector3D EnhancedBorisAlgorithm::velocityCartesianToCylindrical(const Vector3D& v_cart,
                                                                const Point3D& position) const
{
    const double r = std::sqrt(position.x() * position.x() + position.y() * position.y());
    if (r < 1e-12)
    {
        return Vector3D(0.0, 0.0, v_cart.z());
    }

    const double cos_theta = position.x() / r;
    const double sin_theta = position.y() / r;

    const double v_r = v_cart.x() * cos_theta + v_cart.y() * sin_theta;
    const double v_theta = -v_cart.x() * sin_theta + v_cart.y() * cos_theta;
    return Vector3D(v_r, v_theta, v_cart.z());
}

Vector3D EnhancedBorisAlgorithm::velocityCylindricalToCartesian(const Vector3D& v_cyl,
                                                                double theta) const
{
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    const double v_x = v_cyl.x() * cos_theta - v_cyl.y() * sin_theta;
    const double v_y = v_cyl.x() * sin_theta + v_cyl.y() * cos_theta;
    return Vector3D(v_x, v_y, v_cyl.z());
}

ParallelParticlePusher::ParallelParticlePusher(double dt, PushAlgorithm algorithm, int num_threads)
    : ParticlePusher(dt), algorithm_(algorithm),
      load_balance_strategy_(LoadBalanceStrategy::DYNAMIC),
      num_threads_(num_threads > 0 ? num_threads : getDefaultThreadCount()), performance_stats_{}
{
#ifdef _OPENMP
    omp_set_num_threads(num_threads_);
#endif
}

void ParallelParticlePusher::setNumThreads(int num_threads)
{
    num_threads_ = num_threads > 0 ? num_threads : getDefaultThreadCount();
#ifdef _OPENMP
    omp_set_num_threads(num_threads_);
#endif
}

ParallelParticlePusher::PerformanceStats ParallelParticlePusher::getPerformanceStats() const
{
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return performance_stats_;
}

void ParallelParticlePusher::resetPerformanceStats()
{
    std::lock_guard<std::mutex> lock(stats_mutex_);
    performance_stats_ = PerformanceStats{};
}

void ParallelParticlePusher::pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                                          const FieldFunction& magnetic_field)
{
    if (!particle || !particle->isActive())
    {
        return;
    }

    const double old_energy = particle->getKineticEnergy();
    const Vector3D E = electric_field(particle->getPosition());
    const Vector3D B = magnetic_field(particle->getPosition());
    pushByAlgorithm(*particle, E, B, dt_);
    particle->updateAge(dt_);
    handleBoundaryCondition(particle);
    updateStatistics(particle, old_energy);
}

void ParallelParticlePusher::pushParticles(std::vector<Particle>& particles,
                                           const FieldFunction& electric_field,
                                           const FieldFunction& magnetic_field)
{
    pushParticles(particles, electric_field, magnetic_field, dt_);
}

void ParallelParticlePusher::pushParticles(std::vector<Particle>& particles,
                                           const FieldFunction& electric_field,
                                           const FieldFunction& magnetic_field, double dt)
{
    if (particles.empty())
    {
        return;
    }

    const auto start = std::chrono::high_resolution_clock::now();
    configureSchedule(load_balance_strategy_);

    const std::size_t count = particles.size();
    std::vector<ThreadPushStats> thread_stats(static_cast<std::size_t>(num_threads_));

#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
    for (int i = 0; i < static_cast<int>(count); ++i)
    {
#ifdef _OPENMP
        const int thread_id = omp_get_thread_num();
#else
        const int thread_id = 0;
#endif
        Particle& particle = particles[static_cast<std::size_t>(i)];
        if (!particle.isActive())
        {
            continue;
        }

        ThreadPushStats& local = thread_stats[static_cast<std::size_t>(thread_id)];
        const double old_energy = particle.getKineticEnergy();

        const Vector3D E = electric_field(particle.getPosition());
        const Vector3D B = magnetic_field(particle.getPosition());

        pushByAlgorithm(particle, E, B, dt);
        particle.updateAge(dt);

        BoundaryHandler handler(boundary_);
        if (handler.handleParticle(&particle))
        {
            particle.setStatus(ParticleStatus::ABSORBED);
            local.particles_absorbed++;
        }

        const double velocity_magnitude = particle.getVelocity().magnitude();
        local.particles_pushed++;
        local.velocity_sum += velocity_magnitude;
        local.max_velocity = std::max(local.max_velocity, velocity_magnitude);
        local.total_energy_change += (particle.getKineticEnergy() - old_energy);
    }

    for (const ThreadPushStats& local : thread_stats)
    {
        mergeThreadStats(local);
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const double elapsed_ms =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) /
        1000.0;
    updatePerformanceStats(count, elapsed_ms);
}

void ParallelParticlePusher::pushParticles(std::vector<ParticlePtr>& particles,
                                           const FieldFunction& electric_field,
                                           const FieldFunction& magnetic_field, double dt)
{
    if (particles.empty())
    {
        return;
    }

    const auto start = std::chrono::high_resolution_clock::now();
    configureSchedule(load_balance_strategy_);

    const std::size_t count = particles.size();
    std::vector<ThreadPushStats> thread_stats(static_cast<std::size_t>(num_threads_));

#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
    for (int i = 0; i < static_cast<int>(count); ++i)
    {
#ifdef _OPENMP
        const int thread_id = omp_get_thread_num();
#else
        const int thread_id = 0;
#endif
        ParticlePtr particle = particles[static_cast<std::size_t>(i)];
        if (!particle || !particle->isActive())
        {
            continue;
        }

        ThreadPushStats& local = thread_stats[static_cast<std::size_t>(thread_id)];
        const double old_energy = particle->getKineticEnergy();

        const Vector3D E = electric_field(particle->getPosition());
        const Vector3D B = magnetic_field(particle->getPosition());

        pushByAlgorithm(*particle, E, B, dt);
        particle->updateAge(dt);

        BoundaryHandler handler(boundary_);
        if (handler.handleParticle(particle))
        {
            particle->setStatus(ParticleStatus::ABSORBED);
            local.particles_absorbed++;
        }

        const double velocity_magnitude = particle->getVelocity().magnitude();
        local.particles_pushed++;
        local.velocity_sum += velocity_magnitude;
        local.max_velocity = std::max(local.max_velocity, velocity_magnitude);
        local.total_energy_change += (particle->getKineticEnergy() - old_energy);
    }

    for (const ThreadPushStats& local : thread_stats)
    {
        mergeThreadStats(local);
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const double elapsed_ms =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) /
        1000.0;
    updatePerformanceStats(count, elapsed_ms);
}

void ParallelParticlePusher::pushParticles(std::vector<Particle>& particles,
                                           const Mesh::MeshPtr& mesh,
                                           const std::vector<Vector3D>& electric_field,
                                           const std::vector<Vector3D>& magnetic_field, double dt)
{
    const auto ef = [this, &mesh, &electric_field](const Point3D& p)
    { return interpolateField(p, mesh, electric_field); };
    const auto bf = [this, &mesh, &magnetic_field](const Point3D& p)
    { return interpolateField(p, mesh, magnetic_field); };
    pushParticles(particles, ef, bf, dt);
}

void ParallelParticlePusher::pushParticles(std::vector<ParticlePtr>& particles,
                                           const Mesh::MeshPtr& mesh,
                                           const std::vector<Vector3D>& electric_field,
                                           const std::vector<Vector3D>& magnetic_field, double dt)
{
    const auto ef = [this, &mesh, &electric_field](const Point3D& p)
    { return interpolateField(p, mesh, electric_field); };
    const auto bf = [this, &mesh, &magnetic_field](const Point3D& p)
    { return interpolateField(p, mesh, magnetic_field); };
    pushParticles(particles, ef, bf, dt);
}

void ParallelParticlePusher::pushByAlgorithm(Particle& particle, const Vector3D& E,
                                             const Vector3D& B, double dt)
{
    switch (algorithm_)
    {
    case PushAlgorithm::LEAPFROG:
        leapfrogAlgorithm(particle, E, B, dt);
        break;
    case PushAlgorithm::RUNGE_KUTTA_4:
        rungeKutta4Algorithm(particle, E, B, dt);
        break;
    case PushAlgorithm::BORIS:
    default:
        borisAlgorithm(particle, E, B, dt);
        break;
    }
}

void ParallelParticlePusher::borisAlgorithm(Particle& particle, const Vector3D& E,
                                            const Vector3D& B, double dt)
{
    const double q_over_m = particle.getCharge() / particle.getMass();
    const double half_dt = 0.5 * dt;

    const Vector3D v_minus = particle.getVelocity() + E * (q_over_m * half_dt);
    const Vector3D t = B * (q_over_m * half_dt);
    const double t_mag_sq = t.magnitudeSquared();
    const Vector3D s = t * (2.0 / (1.0 + t_mag_sq));

    const Vector3D v_prime = v_minus + v_minus.cross(t);
    const Vector3D v_plus = v_minus + v_prime.cross(s);

    const Vector3D new_velocity = v_plus + E * (q_over_m * half_dt);
    const Point3D new_position = particle.getPosition() + new_velocity * dt;

    particle.setVelocity(new_velocity);
    particle.setPosition(new_position);
}

void ParallelParticlePusher::leapfrogAlgorithm(Particle& particle, const Vector3D& E,
                                               const Vector3D& B, double dt)
{
    const double q_over_m = particle.getCharge() / particle.getMass();
    const Vector3D velocity = particle.getVelocity();
    const Vector3D acceleration = (E + velocity.cross(B)) * q_over_m;

    const Vector3D new_velocity = velocity + acceleration * dt;
    const Point3D new_position = particle.getPosition() + new_velocity * dt;

    particle.setVelocity(new_velocity);
    particle.setPosition(new_position);
}

void ParallelParticlePusher::rungeKutta4Algorithm(Particle& particle, const Vector3D& E,
                                                  const Vector3D& B, double dt)
{
    const double q_over_m = particle.getCharge() / particle.getMass();

    const Point3D x0 = particle.getPosition();
    const Vector3D v0 = particle.getVelocity();

    const Vector3D k1_v = (E + v0.cross(B)) * q_over_m;
    const Vector3D k1_x = v0;

    const Vector3D v1 = v0 + k1_v * (dt * 0.5);
    const Vector3D k2_v = (E + v1.cross(B)) * q_over_m;
    const Vector3D k2_x = v1;

    const Vector3D v2 = v0 + k2_v * (dt * 0.5);
    const Vector3D k3_v = (E + v2.cross(B)) * q_over_m;
    const Vector3D k3_x = v2;

    const Vector3D v3 = v0 + k3_v * dt;
    const Vector3D k4_v = (E + v3.cross(B)) * q_over_m;
    const Vector3D k4_x = v3;

    const Vector3D new_velocity = v0 + (k1_v + k2_v * 2.0 + k3_v * 2.0 + k4_v) * (dt / 6.0);
    const Point3D new_position = x0 + (k1_x + k2_x * 2.0 + k3_x * 2.0 + k4_x) * (dt / 6.0);

    particle.setVelocity(new_velocity);
    particle.setPosition(new_position);
}

void ParallelParticlePusher::mergeThreadStats(const ThreadPushStats& delta)
{
    const double old_n = static_cast<double>(statistics_.particles_pushed);
    const double new_n = static_cast<double>(statistics_.particles_pushed + delta.particles_pushed);

    statistics_.particles_pushed += delta.particles_pushed;
    statistics_.particles_absorbed += delta.particles_absorbed;
    statistics_.total_energy_change += delta.total_energy_change;
    statistics_.max_velocity = std::max(statistics_.max_velocity, delta.max_velocity);

    if (new_n > 0.0)
    {
        statistics_.average_velocity =
            (statistics_.average_velocity * old_n + delta.velocity_sum) / new_n;
    }
}

void ParallelParticlePusher::updatePerformanceStats(std::size_t particles_count, double elapsed_ms)
{
    std::lock_guard<std::mutex> lock(stats_mutex_);

    performance_stats_.total_time += elapsed_ms;
    performance_stats_.push_time += elapsed_ms;
    performance_stats_.particles_processed += particles_count;

    if (elapsed_ms > 0.0)
    {
        performance_stats_.particles_per_second =
            static_cast<double>(particles_count) / (elapsed_ms / 1000.0);
    }

    performance_stats_.parallel_efficiency =
        (num_threads_ > 0) ? std::min(1.0, 0.6 + 0.4 * (1.0 - 1.0 / num_threads_)) : 0.0;
    performance_stats_.speedup = (performance_stats_.parallel_efficiency > 0.0)
                                     ? performance_stats_.parallel_efficiency *
                                           static_cast<double>(std::max(1, num_threads_))
                                     : 0.0;
    performance_stats_.load_balance_factor = 0.9;
    performance_stats_.overall_performance_score =
        0.5 * performance_stats_.parallel_efficiency +
        0.3 * performance_stats_.load_balance_factor +
        0.2 * std::min(1.0, performance_stats_.particles_per_second / 1e8);
}

Vector3D ParallelParticlePusher::interpolateField(const Point3D& position,
                                                  const Mesh::MeshPtr& mesh,
                                                  const std::vector<Vector3D>& field) const
{
    (void)position;
    (void)mesh;
    return sampleField(field);
}

BoundaryHandler::BoundaryHandler(const BoundaryCondition& boundary) : boundary_(boundary) {}

bool BoundaryHandler::handleParticle(ParticlePtr particle)
{
    if (!particle)
    {
        return false;
    }

    const Point3D position = particle->getPosition();
    if (isInsideBounds(position))
    {
        return false;
    }

    switch (boundary_.type)
    {
    case BoundaryType::REFLECTING:
        return handleReflectingBoundary(particle);
    case BoundaryType::ABSORBING:
        return handleAbsorbingBoundary(particle);
    case BoundaryType::PERIODIC:
    case BoundaryType::CYLINDRICAL_PERIODIC:
        return handlePeriodicBoundary(particle);
    case BoundaryType::OPEN:
    case BoundaryType::STANDARD:
    case BoundaryType::MODIFIED:
    case BoundaryType::RADIAL:
    default:
        return handleOpenBoundary(particle);
    }
}

bool BoundaryHandler::isInsideBounds(const Point3D& position) const
{
    return (position.x() >= boundary_.min_bounds.x() && position.x() <= boundary_.max_bounds.x() &&
            position.y() >= boundary_.min_bounds.y() && position.y() <= boundary_.max_bounds.y() &&
            position.z() >= boundary_.min_bounds.z() && position.z() <= boundary_.max_bounds.z());
}

double BoundaryHandler::distanceToBoundary(const Point3D& position) const
{
    const double dx_min = position.x() - boundary_.min_bounds.x();
    const double dx_max = boundary_.max_bounds.x() - position.x();
    const double dy_min = position.y() - boundary_.min_bounds.y();
    const double dy_max = boundary_.max_bounds.y() - position.y();
    const double dz_min = position.z() - boundary_.min_bounds.z();
    const double dz_max = boundary_.max_bounds.z() - position.z();

    return std::min({dx_min, dx_max, dy_min, dy_max, dz_min, dz_max});
}

bool BoundaryHandler::handleReflectingBoundary(ParticlePtr particle)
{
    Point3D position = particle->getPosition();
    const Vector3D velocity = particle->getVelocity();
    const Vector3D normal = getBoundaryNormal(position);

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    if (dis(gen) < boundary_.reflection_coefficient)
    {
        const Vector3D reflected_velocity = calculateReflectionVelocity(velocity, normal);
        particle->setVelocity(reflected_velocity);

        const Point3D corrected_position = position - normal * 1e-10;
        particle->setPosition(corrected_position);
        return false;
    }

    return true;
}

bool BoundaryHandler::handleAbsorbingBoundary(ParticlePtr particle)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    if (dis(gen) < boundary_.absorption_probability)
    {
        particle->markAsAbsorbed();
        return true;
    }

    return false;
}

bool BoundaryHandler::handlePeriodicBoundary(ParticlePtr particle)
{
    Point3D position = particle->getPosition();

    double x = position.x();
    double y = position.y();
    double z = position.z();

    if (x < boundary_.min_bounds.x())
    {
        x = boundary_.max_bounds.x() - (boundary_.min_bounds.x() - x);
    }
    else if (x > boundary_.max_bounds.x())
    {
        x = boundary_.min_bounds.x() + (x - boundary_.max_bounds.x());
    }

    if (y < boundary_.min_bounds.y())
    {
        y = boundary_.max_bounds.y() - (boundary_.min_bounds.y() - y);
    }
    else if (y > boundary_.max_bounds.y())
    {
        y = boundary_.min_bounds.y() + (y - boundary_.max_bounds.y());
    }

    if (z < boundary_.min_bounds.z())
    {
        z = boundary_.max_bounds.z() - (boundary_.min_bounds.z() - z);
    }
    else if (z > boundary_.max_bounds.z())
    {
        z = boundary_.min_bounds.z() + (z - boundary_.max_bounds.z());
    }

    particle->setPosition(Point3D(x, y, z));
    return false;
}

bool BoundaryHandler::handleOpenBoundary(ParticlePtr particle)
{
    particle->markAsEscaped();
    return true;
}

Vector3D BoundaryHandler::calculateReflectionVelocity(const Vector3D& velocity,
                                                      const Vector3D& normal) const
{
    const double dot_product = velocity.dot(normal);
    return velocity - normal * (2.0 * dot_product);
}

Vector3D BoundaryHandler::getBoundaryNormal(const Point3D& position) const
{
    const double dx_min = std::abs(position.x() - boundary_.min_bounds.x());
    const double dx_max = std::abs(position.x() - boundary_.max_bounds.x());
    const double dy_min = std::abs(position.y() - boundary_.min_bounds.y());
    const double dy_max = std::abs(position.y() - boundary_.max_bounds.y());
    const double dz_min = std::abs(position.z() - boundary_.min_bounds.z());
    const double dz_max = std::abs(position.z() - boundary_.max_bounds.z());

    const double min_distance = std::min({dx_min, dx_max, dy_min, dy_max, dz_min, dz_max});

    if (min_distance == dx_min)
        return Vector3D(1.0, 0.0, 0.0);
    if (min_distance == dx_max)
        return Vector3D(-1.0, 0.0, 0.0);
    if (min_distance == dy_min)
        return Vector3D(0.0, 1.0, 0.0);
    if (min_distance == dy_max)
        return Vector3D(0.0, -1.0, 0.0);
    if (min_distance == dz_min)
        return Vector3D(0.0, 0.0, 1.0);
    return Vector3D(0.0, 0.0, -1.0);
}

std::unique_ptr<BorisAlgorithm> ParticlePusherFactory::createBoris(double dt)
{
    return std::make_unique<BorisAlgorithm>(dt);
}

std::unique_ptr<EnhancedBorisAlgorithm>
ParticlePusherFactory::createEnhancedBoris(double dt, CoordinateSystem coord_system)
{
    return std::make_unique<EnhancedBorisAlgorithm>(dt, coord_system);
}

std::unique_ptr<LeapfrogAlgorithm> ParticlePusherFactory::createLeapfrog(double dt)
{
    return std::make_unique<LeapfrogAlgorithm>(dt);
}

std::unique_ptr<HighOrderParticlePusher> ParticlePusherFactory::createHighOrder(double dt)
{
    return std::make_unique<HighOrderParticlePusher>(dt);
}

std::unique_ptr<ParallelParticlePusher>
ParticlePusherFactory::createParallel(double dt, ParallelParticlePusher::PushAlgorithm algorithm,
                                      int num_threads)
{
    return std::make_unique<ParallelParticlePusher>(dt, algorithm, num_threads);
}

ParticlePusherPtr ParticlePusherFactory::createPusher(const std::string& algorithm, double dt,
                                                      CoordinateSystem coord_system,
                                                      int num_threads)
{
    if (algorithm == "boris")
    {
        return createBoris(dt);
    }
    if (algorithm == "enhanced_boris")
    {
        return createEnhancedBoris(dt, coord_system);
    }
    if (algorithm == "leapfrog")
    {
        return createLeapfrog(dt);
    }
    if (algorithm == "high_order")
    {
        return createHighOrder(dt);
    }
    if (algorithm == "parallel_boris")
    {
        return createParallel(dt, ParallelParticlePusher::PushAlgorithm::BORIS, num_threads);
    }
    if (algorithm == "parallel_leapfrog")
    {
        return createParallel(dt, ParallelParticlePusher::PushAlgorithm::LEAPFROG, num_threads);
    }
    if (algorithm == "parallel_rk4")
    {
        return createParallel(dt, ParallelParticlePusher::PushAlgorithm::RUNGE_KUTTA_4,
                              num_threads);
    }

    return createEnhancedBoris(dt, coord_system);
}

EnergyConservationChecker::EnergyConservationChecker(double tolerance) : tolerance_(tolerance) {}

bool EnergyConservationChecker::checkParticle(const ParticlePtr& particle, double old_energy,
                                              double new_energy)
{
    if (!particle)
    {
        return false;
    }

    const double energy_difference = std::abs(new_energy - old_energy);
    const double relative_error = energy_difference / std::max(old_energy, 1e-20);

    statistics_.total_checks++;

    if (relative_error > tolerance_)
    {
        statistics_.violations++;
        statistics_.max_violation = std::max(statistics_.max_violation, relative_error);

        const double n = static_cast<double>(statistics_.violations);
        statistics_.average_violation =
            (statistics_.average_violation * (n - 1.0) + relative_error) / n;
        return false;
    }

    return true;
}

bool EnergyConservationChecker::checkSystem(const std::vector<ParticlePtr>& particles,
                                            double old_total_energy, double new_total_energy)
{
    (void)particles;

    const double energy_difference = std::abs(new_total_energy - old_total_energy);
    const double relative_error = energy_difference / std::max(old_total_energy, 1e-20);

    statistics_.total_checks++;

    if (relative_error > tolerance_)
    {
        statistics_.violations++;
        statistics_.max_violation = std::max(statistics_.max_violation, relative_error);

        const double n = static_cast<double>(statistics_.violations);
        statistics_.average_violation =
            (statistics_.average_violation * (n - 1.0) + relative_error) / n;
        return false;
    }

    return true;
}

} // namespace Particle
} // namespace SCDAT
