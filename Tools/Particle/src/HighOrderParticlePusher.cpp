#include "HighOrderParticlePusher.h"

#include <chrono>

namespace SCDAT
{
namespace Particle
{

HighOrderParticlePusher::HighOrderParticlePusher(double dt, MethodType method_type)
    : ParticlePusher(dt), method_type_(method_type)
{
}

void HighOrderParticlePusher::pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                                           const FieldFunction& magnetic_field)
{
    if (!particle || !particle->isActive())
    {
        return;
    }

    const Point3D position = particle->getPosition();
    const Vector3D electric = electric_field(position);
    const Vector3D magnetic = magnetic_field(position);
    rungeKutta4Push(particle, electric, magnetic, dt_);
    stats_.pushed_particles += 1;
}

bool HighOrderParticlePusher::pushParticlesHighOrder(std::vector<ParticlePtr>& particles,
                                                     const std::vector<Vector3D>& electric_field,
                                                     const std::vector<Vector3D>& magnetic_field,
                                                     double dt)
{
    const auto start = std::chrono::high_resolution_clock::now();
    dt_ = dt;

    for (auto& particle : particles)
    {
        if (!particle || !particle->isActive())
        {
            continue;
        }

        const Vector3D electric = interpolateFieldAtParticle(particle, electric_field);
        const Vector3D magnetic = interpolateFieldAtParticle(particle, magnetic_field);
        rungeKutta4Push(particle, electric, magnetic, dt);
        stats_.pushed_particles += 1;
    }

    const auto end = std::chrono::high_resolution_clock::now();
    stats_.solve_time_ms += std::chrono::duration<double, std::milli>(end - start).count();
    return true;
}

void HighOrderParticlePusher::pushParticlesWithHistory(
    std::vector<ParticlePtr>& particles,
    const std::vector<std::vector<Vector3D>>& electric_field_history,
    const std::vector<std::vector<Vector3D>>& magnetic_field_history, double dt)
{
    dt_ = dt;
    for (auto& particle : particles)
    {
        if (!particle || !particle->isActive())
        {
            continue;
        }

        if (method_type_ == MethodType::ADAMS_BASHFORTH_6)
        {
            adamsBashforth6Push(particle, electric_field_history, magnetic_field_history, dt);
        }
        else
        {
            const Vector3D electric = electric_field_history.empty()
                                          ? Vector3D(0, 0, 0)
                                          : interpolateFieldAtParticle(particle, electric_field_history.back());
            const Vector3D magnetic = magnetic_field_history.empty()
                                          ? Vector3D(0, 0, 0)
                                          : interpolateFieldAtParticle(particle, magnetic_field_history.back());
            rungeKutta4Push(particle, electric, magnetic, dt);
        }
        stats_.pushed_particles += 1;
    }
}

void HighOrderParticlePusher::rungeKutta4Push(ParticlePtr particle, const Vector3D& electric_field,
                                              const Vector3D& magnetic_field, double dt)
{
    if (!particle || particle->getMass() <= 0.0)
    {
        return;
    }

    const Point3D position = particle->getPosition();
    const Vector3D velocity = particle->getVelocity();
    const double q_over_m = particle->getCharge() / particle->getMass();

    const Vector3D k1_v = (electric_field + velocity.cross(magnetic_field)) * q_over_m;
    const Vector3D k1_x = velocity;

    const Vector3D v1 = velocity + k1_v * (dt * 0.5);
    const Vector3D k2_v = (electric_field + v1.cross(magnetic_field)) * q_over_m;
    const Vector3D k2_x = v1;

    const Vector3D v2 = velocity + k2_v * (dt * 0.5);
    const Vector3D k3_v = (electric_field + v2.cross(magnetic_field)) * q_over_m;
    const Vector3D k3_x = v2;

    const Vector3D v3 = velocity + k3_v * dt;
    const Vector3D k4_v = (electric_field + v3.cross(magnetic_field)) * q_over_m;
    const Vector3D k4_x = v3;

    const Vector3D new_velocity =
        velocity + (k1_v + k2_v * 2.0 + k3_v * 2.0 + k4_v) * (dt / 6.0);
    const Point3D new_position =
        position + (k1_x + k2_x * 2.0 + k3_x * 2.0 + k4_x) * (dt / 6.0);

    particle->setVelocity(new_velocity);
    particle->setPosition(new_position);
}

void HighOrderParticlePusher::adamsBashforth6Push(
    ParticlePtr particle, const std::vector<std::vector<Vector3D>>& electric_field_history,
    const std::vector<std::vector<Vector3D>>& magnetic_field_history, double dt)
{
    if (!particle || electric_field_history.empty() || magnetic_field_history.empty())
    {
        return;
    }

    const Vector3D electric = interpolateFieldAtParticle(particle, electric_field_history.back());
    const Vector3D magnetic = interpolateFieldAtParticle(particle, magnetic_field_history.back());
    rungeKutta4Push(particle, electric, magnetic, dt);
}

Vector3D HighOrderParticlePusher::interpolateFieldAtParticle(const ParticlePtr& particle,
                                                             const std::vector<Vector3D>& field) const
{
    if (!particle || field.empty())
    {
        return Vector3D(0, 0, 0);
    }

    const Point3D position = particle->getPosition();
    const std::size_t index =
        static_cast<std::size_t>(std::abs(position.x() + position.y() + position.z())) %
        field.size();
    return field[index];
}

} // namespace Particle
} // namespace SCDAT
