#ifndef SCDAT_PARTICLE_PUSHER_H
#define SCDAT_PARTICLE_PUSHER_H

/**
 * @file ParticlePusher.h
 * @ingroup ParticleModule
 */

#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
#include "ParticleManager.h"
#include <cstddef>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Mesh
{
class Mesh;
using MeshPtr = std::shared_ptr<Mesh>;
} // namespace Mesh
} // namespace SCDAT

namespace SCDAT
{
namespace Particle
{

// 前向声明
class ParticlePusher;
class BorisAlgorithm;
class EnhancedBorisAlgorithm;
class LeapfrogAlgorithm;
class ParallelParticlePusher;
class BoundaryHandler;

using Point3D = SCDAT::Geometry::Point3D;
using Vector3D = SCDAT::Geometry::Vector3D;
using ParticlePusherPtr = std::shared_ptr<ParticlePusher>;
using FieldFunction = std::function<Vector3D(const Point3D&)>;

/**
 * @brief 边界条件类型枚举
 * @details 多种边界条件类型 
 */
enum class BoundaryType
{
    REFLECTING = 0, ///< 反射边界 
    ABSORBING = 1,  ///< 吸收边界
    PERIODIC = 2,   ///< 周期边界
    OPEN = 3,       ///< 开放边界
    STANDARD = 4,   ///< 标准边界条件
    MODIFIED =5,    ///< 修改的边界条件
    RADIAL = 6,     ///< 径向边界条件
    CYLINDRICAL_PERIODIC = 7 ///< 柱坐标周期边界
};

/**
 * @brief 坐标系类型枚举
 * @details 支持笛卡尔和柱坐标系
 */
enum class CoordinateSystem
{
    CARTESIAN = 0,  ///< 笛卡尔坐标系
    CYLINDRICAL = 1 ///< 柱坐标系 (r, θ, z)
};

/**
 * @brief 边界条件结构
 */
struct BoundaryCondition
{
    BoundaryType type = BoundaryType::REFLECTING;
    Point3D min_bounds = Point3D(-1.0, -1.0, -1.0);
    Point3D max_bounds = Point3D(1.0, 1.0, 1.0);
    double reflection_coefficient = 1.0;
    double absorption_probability = 0.0;
};

class ParticlePusher
{
  public:
    /**
     * @brief 构造函数
     * @param dt 时间步长 (s)
     */
    explicit ParticlePusher(double dt);

    /**
     * @brief 虚析构函数
     */
    virtual ~ParticlePusher() = default;

    /**
     * @brief 推进单个粒子
     * @param particle 要推进的粒子
     * @param electric_field 电场函数
     * @param magnetic_field 磁场函数
     */
    virtual void pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                              const FieldFunction& magnetic_field) = 0;

    /**
     * @brief 批量推进粒子
     * @param particles 粒子列表
     * @param electric_field 电场函数
     * @param magnetic_field 磁场函数
     */
    virtual void pushParticles(const std::vector<ParticlePtr>& particles,
                               const FieldFunction& electric_field,
                               const FieldFunction& magnetic_field);

    /**
     * @brief 批量推进粒子 (值类型版本)
     */
    virtual void pushParticles(std::vector<Particle>& particles,
                               const FieldFunction& electric_field,
                               const FieldFunction& magnetic_field);

    /**
     * @brief 推进粒子（函数场版本）
     * @param particles 粒子列表（非const引用，允许修改）
     * @param electric_field 电场函数
     * @param magnetic_field 磁场函数
     * @param dt 时间步长
     */
    virtual void pushParticles(std::vector<ParticlePtr>& particles,
                               const FieldFunction& electric_field,
                               const FieldFunction& magnetic_field, double dt);

    /**
     * @brief 推进粒子（函数场版本，值类型）
     */
    virtual void pushParticles(std::vector<Particle>& particles,
                               const FieldFunction& electric_field,
                               const FieldFunction& magnetic_field, double dt);

    /**
     * @brief 推进粒子（网格场版本）
     * @param particles 粒子列表（非const引用，允许修改）
     * @param mesh 网格指针
     * @param electric_field 电场数据
     * @param magnetic_field 磁场数据
     * @param dt 时间步长
     */
    virtual void pushParticles(std::vector<ParticlePtr>& particles, const Mesh::MeshPtr& mesh,
                               const std::vector<Vector3D>& electric_field,
                               const std::vector<Vector3D>& magnetic_field, double dt);

    /**
     * @brief 推进粒子（网格场版本，值类型）
     */
    virtual void pushParticles(std::vector<Particle>& particles, const Mesh::MeshPtr& mesh,
                               const std::vector<Vector3D>& electric_field,
                               const std::vector<Vector3D>& magnetic_field, double dt);

    /**
     * @brief 推进粒子管理器中的所有粒子
     * @param manager 粒子管理器
     * @param electric_field 电场函数
     * @param magnetic_field 磁场函数
     */
    void pushAllParticles(ParticleManager& manager, const FieldFunction& electric_field,
                          const FieldFunction& magnetic_field);

    // 参数设置
    void setTimeStep(double dt)
    {
        dt_ = dt;
    }
    double getTimeStep() const
    {
        return dt_;
    }

    void setRelativistic(bool relativistic)
    {
        relativistic_ = relativistic;
    }
    bool isRelativistic() const
    {
        return relativistic_;
    }

    void setBoundaryCondition(const BoundaryCondition& boundary)
    {
        boundary_ = boundary;
    }
    const BoundaryCondition& getBoundaryCondition() const
    {
        return boundary_;
    }

    // 统计信息
    struct PushStatistics
    {
        std::size_t particles_pushed = 0;
        std::size_t particles_absorbed = 0;
        std::size_t particles_reflected = 0;
        std::size_t particles_escaped = 0;
        double total_energy_change = 0.0;
        double max_velocity = 0.0;
        double average_velocity = 0.0;
    };

    const PushStatistics& getStatistics() const
    {
        return statistics_;
    }
    void resetStatistics()
    {
        statistics_ = PushStatistics{};
    }

  protected:
    double dt_;                  ///< 时间步长 (s)
    bool relativistic_;          ///< 是否使用相对论修正
    BoundaryCondition boundary_; ///< 边界条件
    PushStatistics statistics_;  ///< 统计信息

    // 辅助方法
    void handleBoundaryCondition(ParticlePtr particle);
    void updateStatistics(const ParticlePtr& particle, double old_energy);

    // 物理常数
    static constexpr double LIGHT_SPEED = 299792458.0; ///< 光速 (m/s)
};

/**
 * @brief Boris算法推进器
 *
 * 实现Boris算法进行粒子推进，适用于强磁场环境
 */
class BorisAlgorithm : public ParticlePusher
{
  public:
    /**
     * @brief 构造函数
     * @param dt 时间步长 (s) / Time step (s)
     */
    explicit BorisAlgorithm(double dt);

    /**
     * @brief 推进单个粒子
     */
    void pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                      const FieldFunction& magnetic_field) override;

  private:
    std::pair<Point3D, Vector3D> borisStep(const Point3D& position, const Vector3D& velocity,
                                           double charge, double mass, const Vector3D& E,
                                           const Vector3D& B, double dt);

    std::pair<Point3D, Vector3D>
    relativisticBorisStep(const Point3D& position, const Vector3D& velocity, double charge,
                          double mass, const Vector3D& E, const Vector3D& B, double dt);
};

/**
 * @brief 增强版Boris算法推进器
 * @details 支持柱坐标系和多边界条件
 */
class EnhancedBorisAlgorithm : public ParticlePusher
{
  public:
    /**
     * @brief 构造函数 
     * @param dt 时间步长 (s)
     * @param coord_system 坐标系类型
     */
    explicit EnhancedBorisAlgorithm(double dt,
                                    CoordinateSystem coord_system = CoordinateSystem::CARTESIAN);

    /**
     * @brief 推进单个粒子 / Push single particle
     */
    void pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                      const FieldFunction& magnetic_field) override;

    /**
     * @brief 设置坐标系类型 / Set coordinate system type
     */
    void setCoordinateSystem(CoordinateSystem coord_system)
    {
        coord_system_ = coord_system;
    }

    /**
     * @brief 获取坐标系类型 / Get coordinate system type
     */
    CoordinateSystem getCoordinateSystem() const
    {
        return coord_system_;
    }

    /**
     * @brief 设置磁场配置 / Set magnetic field configuration
     * @param enable_magnetic 是否启用磁场 / Enable magnetic field
     * @param B_uniform 均匀磁场强度 / Uniform magnetic field strength
     * @param B_direction 磁场方向 / Magnetic field direction
     */
    void setMagneticFieldConfig(bool enable_magnetic, double B_uniform = 0.0,
                                const Vector3D& B_direction = Vector3D(0, 0, 1));

  private:
    CoordinateSystem coord_system_;
    bool enable_magnetic_;
    double B_uniform_;
    Vector3D B_direction_;

    std::pair<Point3D, Vector3D> arcPICBorisStep(const Point3D& position, const Vector3D& velocity,
                                                 double charge, double mass, const Vector3D& E,
                                                 const Vector3D& B, double dt);

    std::pair<Point3D, Vector3D>
    cylindricalBorisStep(const Point3D& position, const Vector3D& velocity, double charge,
                         double mass, const Vector3D& E, const Vector3D& B, double dt);

    Point3D cartesianToCylindrical(const Point3D& cartesian) const;
    Point3D cylindricalToCartesian(const Point3D& cylindrical, double theta) const;
    Vector3D velocityCartesianToCylindrical(const Vector3D& v_cart,
                                            const Point3D& position) const;
    Vector3D velocityCylindricalToCartesian(const Vector3D& v_cyl, double theta) const;
};

/**
 * @brief Leapfrog算法推进器
 *
 * 实现Leapfrog算法进行粒子推进，具有良好的能量守恒性质
 */
class LeapfrogAlgorithm : public ParticlePusher
{
  public:
    /**
     * @brief 构造函数
     * @param dt 时间步长 (s)
     */
    explicit LeapfrogAlgorithm(double dt);

    /**
     * @brief 推进单个粒子
     */
    void pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                      const FieldFunction& magnetic_field) override;

  private:
    std::pair<Point3D, Vector3D> leapfrogStep(const Point3D& position, const Vector3D& velocity,
                                              double charge, double mass, const Vector3D& E,
                                              const Vector3D& B, double dt);

    bool first_step_;
};

class ParallelParticlePusher : public ParticlePusher
{
  public:
    enum class PushAlgorithm
    {
        BORIS,
        LEAPFROG,
        RUNGE_KUTTA_4
    };

    enum class LoadBalanceStrategy
    {
        STATIC,
        DYNAMIC,
        GUIDED
    };

    struct PerformanceStats
    {
        double total_time = 0.0;
        double push_time = 0.0;
        double field_interpolation_time = 0.0;
        std::size_t particles_processed = 0;
        double particles_per_second = 0.0;
        double parallel_efficiency = 0.0;
        double speedup = 0.0;
        double load_balance_factor = 0.0;
        double overall_performance_score = 0.0;
    };

    explicit ParallelParticlePusher(double dt, PushAlgorithm algorithm = PushAlgorithm::BORIS,
                                    int num_threads = 0);

    void pushParticle(ParticlePtr particle, const FieldFunction& electric_field,
                      const FieldFunction& magnetic_field) override;

    void pushParticles(std::vector<Particle>& particles, const FieldFunction& electric_field,
                       const FieldFunction& magnetic_field) override;

    void pushParticles(std::vector<Particle>& particles, const FieldFunction& electric_field,
                       const FieldFunction& magnetic_field, double dt) override;

    void pushParticles(std::vector<ParticlePtr>& particles, const FieldFunction& electric_field,
                       const FieldFunction& magnetic_field, double dt) override;

    void pushParticles(std::vector<Particle>& particles, const Mesh::MeshPtr& mesh,
                       const std::vector<Vector3D>& electric_field,
                       const std::vector<Vector3D>& magnetic_field, double dt) override;

    void pushParticles(std::vector<ParticlePtr>& particles, const Mesh::MeshPtr& mesh,
                       const std::vector<Vector3D>& electric_field,
                       const std::vector<Vector3D>& magnetic_field, double dt) override;

    void setPushAlgorithm(PushAlgorithm algorithm)
    {
        algorithm_ = algorithm;
    }

    void setLoadBalanceStrategy(LoadBalanceStrategy strategy)
    {
        load_balance_strategy_ = strategy;
    }

    void setNumThreads(int num_threads);

    int getNumThreads() const
    {
        return num_threads_;
    }

    PerformanceStats getPerformanceStats() const;
    void resetPerformanceStats();

  private:
    struct ThreadPushStats
    {
        std::size_t particles_pushed = 0;
        std::size_t particles_absorbed = 0;
        double total_energy_change = 0.0;
        double velocity_sum = 0.0;
        double max_velocity = 0.0;
    };

    PushAlgorithm algorithm_;
    LoadBalanceStrategy load_balance_strategy_;
    int num_threads_;

    mutable std::mutex stats_mutex_;
    PerformanceStats performance_stats_;

    void pushByAlgorithm(Particle& particle, const Vector3D& E, const Vector3D& B, double dt);
    void borisAlgorithm(Particle& particle, const Vector3D& E, const Vector3D& B, double dt);
    void leapfrogAlgorithm(Particle& particle, const Vector3D& E, const Vector3D& B, double dt);
    void rungeKutta4Algorithm(Particle& particle, const Vector3D& E, const Vector3D& B,
                              double dt);

    void mergeThreadStats(const ThreadPushStats& delta);
    void updatePerformanceStats(std::size_t particles_count, double elapsed_ms);

    Vector3D interpolateField(const Point3D& position, const Mesh::MeshPtr& mesh,
                              const std::vector<Vector3D>& field) const;
};

class BoundaryHandler
{
  public:
    /**
     * @brief 构造函数
     * @param boundary 边界条件
     */
    explicit BoundaryHandler(const BoundaryCondition& boundary);

    /**
     * @brief 处理粒子边界条件
     * @param particle 粒子
     * @return 是否需要移除粒子
     */
    bool handleParticle(ParticlePtr particle);
    bool isInsideBounds(const Point3D& position) const;
    double distanceToBoundary(const Point3D& position) const;

  private:
    BoundaryCondition boundary_;

    // 边界处理方法
    bool handleReflectingBoundary(ParticlePtr particle);
    bool handleAbsorbingBoundary(ParticlePtr particle);
    bool handlePeriodicBoundary(ParticlePtr particle);
    bool handleOpenBoundary(ParticlePtr particle);

    Vector3D calculateReflectionVelocity(const Vector3D& velocity, const Vector3D& normal) const;
    Vector3D getBoundaryNormal(const Point3D& position) const;
};

/**
 * @brief 粒子推进器工厂类 / Particle Pusher Factory Class
 * @details 增强工厂支持Arc-PIC算法 / Enhanced factory supporting Arc-PIC algorithms
 */
class ParticlePusherFactory
{
  public:
    /**
     * @brief 创建Boris算法推进器 / Create Boris algorithm pusher
     */
    static std::unique_ptr<BorisAlgorithm> createBoris(double dt);

    /**
     * @brief 创建增强版Boris算法推进器 / Create enhanced Boris algorithm pusher
     * @param dt 时间步长 / Time step
     * @param coord_system 坐标系类型 / Coordinate system type
     */
    static std::unique_ptr<EnhancedBorisAlgorithm>
    createEnhancedBoris(double dt, CoordinateSystem coord_system = CoordinateSystem::CARTESIAN);

    /**
     * @brief 创建Leapfrog算法推进器 / Create Leapfrog algorithm pusher
     */
    static std::unique_ptr<LeapfrogAlgorithm> createLeapfrog(double dt);
    static std::unique_ptr<ParallelParticlePusher>
    createParallel(double dt,
                   ParallelParticlePusher::PushAlgorithm algorithm =
                       ParallelParticlePusher::PushAlgorithm::BORIS,
                   int num_threads = 0);

    static ParticlePusherPtr createPusher(
        const std::string& algorithm, double dt,
        CoordinateSystem coord_system = CoordinateSystem::CARTESIAN, int num_threads = 0);
};

class EnergyConservationChecker
{
  public:
    /**
     * @brief 构造函数
     * @param tolerance 能量守恒容差
     */
    explicit EnergyConservationChecker(double tolerance = 1e-12);

    /**
     * @brief 检查单个粒子的能量守恒
     * @param particle 粒子
     * @param old_energy 推进前的能量
     * @param new_energy 推进后的能量
     * @return 是否满足能量守恒
     */
    bool checkParticle(const ParticlePtr& particle, double old_energy, double new_energy);

    /**
     * @brief 检查粒子系统的总能量守恒
     * @param particles 粒子列表
     * @param old_total_energy 推进前的总能量
     * @param new_total_energy 推进后的总能量
     * @return 是否满足能量守恒
     */
    bool checkSystem(const std::vector<ParticlePtr>& particles, double old_total_energy,
                     double new_total_energy);

    /**
     * @brief 获取统计信息
     */
    struct ConservationStatistics
    {
        std::size_t total_checks = 0;
        std::size_t violations = 0;
        double max_violation = 0.0;
        double average_violation = 0.0;
    };

    const ConservationStatistics& getStatistics() const
    {
        return statistics_;
    }
    void resetStatistics()
    {
        statistics_ = ConservationStatistics{};
    }

  private:
    double tolerance_;
    ConservationStatistics statistics_;
};

} // namespace Particle
} // namespace SCDAT

#endif // SCDAT_PARTICLE_PUSHER_H
