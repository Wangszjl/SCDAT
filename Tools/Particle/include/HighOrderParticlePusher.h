#ifndef SCDAT_HIGH_ORDER_PARTICLE_PUSHER_H
#define SCDAT_HIGH_ORDER_PARTICLE_PUSHER_H

#include "ParticlePusher.h"
#include <memory>
#include <vector>

namespace SCDAT
{
namespace Particle
{

class HighOrderParticlePusher : public ParticlePusher
{
  public:
    enum class MethodType
    {
        RUNGE_KUTTA_4,
        ADAMS_BASHFORTH_6
    };

    struct Statistics
    {
        double solve_time_ms = 0.0;
        std::size_t pushed_particles = 0;
    };

    explicit HighOrderParticlePusher(double dt, MethodType method_type = MethodType::RUNGE_KUTTA_4);

    void pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                      const FieldFunction& magnetic_field) override;

    bool pushParticlesHighOrder(std::vector<ParticlePtr>& particles,
                                const std::vector<Vector3D>& electric_field,
                                const std::vector<Vector3D>& magnetic_field, double dt);

    void pushParticlesWithHistory(
        std::vector<ParticlePtr>& particles,
        const std::vector<std::vector<Vector3D>>& electric_field_history,
        const std::vector<std::vector<Vector3D>>& magnetic_field_history, double dt);

    void setMethodType(MethodType method_type) { method_type_ = method_type; }
    MethodType getMethodType() const { return method_type_; }

    const Statistics& getStatistics() const { return stats_; }
    void resetStatistics() { stats_ = Statistics{}; }

  private:
    void rungeKutta4Push(ParticlePtr particle, const Vector3D& electric_field,
                         const Vector3D& magnetic_field, double dt);

    void adamsBashforth6Push(
        ParticlePtr particle, const std::vector<std::vector<Vector3D>>& electric_field_history,
        const std::vector<std::vector<Vector3D>>& magnetic_field_history, double dt);

    Vector3D interpolateFieldAtParticle(const ParticlePtr& particle,
                                        const std::vector<Vector3D>& field) const;

  private:
    MethodType method_type_;
    Statistics stats_;
};

using HighOrderParticlePusherPtr = std::shared_ptr<HighOrderParticlePusher>;

} // namespace Particle
} // namespace SCDAT

#endif // SCDAT_HIGH_ORDER_PARTICLE_PUSHER_H
