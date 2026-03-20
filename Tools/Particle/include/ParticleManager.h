#ifndef SCDAT_PARTICLE_MANAGER_H
#define SCDAT_PARTICLE_MANAGER_H

/**
 * @file ParticleManager.h
 * @ingroup ParticleModule
 */

#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
#include "ParticleDefinitions.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <execution>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Particle
{

// 前向声明
class ParticleContainer;
class ParticleDataManager;
class ParticleManager;
class ParticleIterator;
class CollisionDetector;

// 类型别名
using ParticleContainerPtr = std::shared_ptr<ParticleContainer>;
using ParticleManagerPtr = std::shared_ptr<ParticleManager>;

// 兼容旧命名
using ParticleManagerContainer = ParticleContainer;
using ParticleDataContainer = ParticleDataManager;

/**
 * @brief 粒子容器类
 *
 * 高效存储和管理粒子的容器，支持快速添加、删除和查找
 */
class ParticleContainer
{
    friend class ParticleIterator; // 允许ParticleIterator访问私有成员

  public:
    /**
     * @brief 构造函数
     * @param initial_capacity 初始容量
     */
    explicit ParticleContainer(std::size_t initial_capacity = 10000);

    /**
     * @brief 析构函数
     */
    ~ParticleContainer() = default;

    // 基本操作
    void addParticle(const Particle& particle);
    bool removeParticle(ParticleId id);
    ParticlePtr getParticle(ParticleId id) const;

    // 批量操作
    void addParticles(const std::vector<Particle>& particles);
    std::size_t removeInactiveParticles();
    std::vector<ParticlePtr> getParticlesByType(ParticleType type) const; // 获取特定类型的粒子列表

    // 容器信息
    /**
     * @brief 获取粒子数量、容量和是否为空
     * @return 粒子数量、容量和是否为空
     */
    std::size_t size() const
    {
        return particles_.size();
    }
    std::size_t capacity() const
    {
        return particles_.capacity();
    }
    bool empty() const
    {
        return particles_.empty();
    }

    // 统计信息
    std::size_t getActiveCount() const;
    std::size_t getTypeCount(ParticleType type) const;
    std::unordered_map<ParticleType, std::size_t> getTypeCounts() const;

    // 内存管理
    void reserve(std::size_t capacity);
    void shrinkToFit(); // 将粒子容器的内存容量收缩到正好能容纳当前元素数量，释放多余的内存
    void clear();

    // 迭代器支持
    std::vector<Particle>::iterator begin()
    {
        return particles_.begin();
    }
    std::vector<Particle>::iterator end()
    {
        return particles_.end();
    }
    std::vector<Particle>::const_iterator begin() const
    {
        return particles_.begin();
    }
    std::vector<Particle>::const_iterator end() const
    {
        return particles_.end();
    }

    // 线程安全访问
    void lockForRead() const
    {
        mutex_.lock_shared();
    }
    void unlockForRead() const
    {
        mutex_.unlock_shared();
    }
    void lockForWrite()
    {
        mutex_.lock();
    }
    void unlockForWrite()
    {
        mutex_.unlock();
    }

  private:
    std::vector<Particle> particles_;                    ///< 粒子存储
    std::unordered_map<ParticleId, std::size_t> id_map_; ///< ID到索引的映射
    mutable std::shared_mutex mutex_;                    ///< 读写锁

    // 内部辅助方法
    void
    rebuildIdMap(); // 重新构建
                    // id_map_，即根据当前粒子列表，把每个粒子的唯一ID映射到其在容器中的索引。常用于批量增删粒子后，保证通过ID能快速查找粒子
    void
    compactStorage(); // 整理和压缩底层存储结构，通常用于移除无效/已删除粒子后，减少碎片、优化内存布局，提高遍历和访问效率。
};

/**
 * @brief 粒子数据容器类（SoA + 批量推进）
 *
 */
class ParticleDataManager
{
  public:
    using ParticleData = SCDAT::Particle::ParticleData;
    using Container = std::pmr::vector<ParticleData>;

    struct Statistics
    {
        size_t total_particles = 0;
        size_t electrons = 0;
        size_t ions = 0;
        size_t capacity = 0;
        double memory_mb = 0.0;

        void update(size_t n, size_t cap, size_t n_electrons, size_t n_ions)
        {
            total_particles = n;
            capacity = cap;
            electrons = n_electrons;
            ions = n_ions;
            memory_mb = cap * 9 * sizeof(double) / (1024.0 * 1024.0);
        }
    };

    explicit ParticleDataManager(
        std::pmr::memory_resource* memory_resource = std::pmr::get_default_resource())
        : particles_(memory_resource), memory_resource_(memory_resource)
    {
    }

    ParticleDataManager(const ParticleDataManager&) = delete;
    ParticleDataManager& operator=(const ParticleDataManager&) = delete;
    ParticleDataManager(ParticleDataManager&&) noexcept = default;
    ParticleDataManager& operator=(ParticleDataManager&&) noexcept = default;

    ~ParticleDataManager() = default;

    // -----------------------------------------------------------------
    // Core Container API（推荐优先使用）
    // -----------------------------------------------------------------
    void reserve(size_t capacity)
    {
        syncSoAToParticlesIfNeeded();
        particles_.reserve(capacity);
        invalidateDerivedCaches();
    }
    void resize(size_t size)
    {
        syncSoAToParticlesIfNeeded();
        particles_.resize(size);
        invalidateDerivedCaches();
    }
    void clear()
    {
        syncSoAToParticlesIfNeeded();
        particles_.clear();
        invalidateDerivedCaches();
    }

    // -----------------------------------------------------------------
    // SoA Compatibility API（为兼容旧 ParticleSoA 接口保留）
    // -----------------------------------------------------------------
    void push(double x, double y, double z, double vx, double vy, double vz, double charge,
              double mass, double weight = 1.0)
    {
        syncSoAToParticlesIfNeeded();

        ParticleData particle{};
        particle.position[0] = x;
        particle.position[1] = y;
        particle.position[2] = z;
        particle.velocity[0] = vx;
        particle.velocity[1] = vy;
        particle.velocity[2] = vz;
        particle.charge = charge;
        particle.mass = mass;
        particle.weight = weight;
        particle.id = particles_.size();
        particle.status = ParticleStatus::ACTIVE;
        particle.type = charge < 0.0   ? ParticleType::ELECTRON
                        : charge > 0.0 ? ParticleType::ION
                                       : ParticleType::UNKNOWN;
        particle.updateEnergy();

        particles_.push_back(particle);
        invalidateDerivedCaches();
    }

    void erase(size_t index)
    {
        syncSoAToParticlesIfNeeded();
        if (index >= particles_.size())
        {
            throw std::out_of_range("Index out of range");
        }

        if (index != particles_.size() - 1)
        {
            std::swap(particles_[index], particles_.back());
        }
        particles_.pop_back();
        invalidateDerivedCaches();
    }

    void shrink_to_fit()
    {
        syncSoAToParticlesIfNeeded();
        particles_.shrink_to_fit();
        invalidateDerivedCaches();
    }

    [[nodiscard]] size_t size() const noexcept
    {
        return particles_.size();
    }
    [[nodiscard]] size_t capacity() const noexcept
    {
        return particles_.capacity();
    }
    [[nodiscard]] bool empty() const noexcept
    {
        return particles_.empty();
    }

    [[nodiscard]] ParticleData& operator[](size_t index)
    {
        syncSoAToParticlesIfNeeded();
        invalidateDerivedCaches();
        return particles_[index];
    }
    [[nodiscard]] const ParticleData& operator[](size_t index) const
    {
        return particles_[index];
    }

    [[nodiscard]] ParticleData& at(size_t index)
    {
        syncSoAToParticlesIfNeeded();
        invalidateDerivedCaches();
        return particles_.at(index);
    }
    [[nodiscard]] const ParticleData& at(size_t index) const
    {
        return particles_.at(index);
    }

    [[nodiscard]] std::span<ParticleData> getParticles()
    {
        syncSoAToParticlesIfNeeded();
        invalidateDerivedCaches();
        return std::span(particles_);
    }
    [[nodiscard]] std::span<const ParticleData> getParticles() const
    {
        return std::span(particles_);
    }

    [[nodiscard]] std::span<ParticleData> data()
    {
        return getParticles();
    }
    [[nodiscard]] std::span<const ParticleData> data() const
    {
        return getParticles();
    }

    auto begin()
    {
        syncSoAToParticlesIfNeeded();
        invalidateDerivedCaches();
        return particles_.begin();
    }
    auto end()
    {
        syncSoAToParticlesIfNeeded();
        invalidateDerivedCaches();
        return particles_.end();
    }
    auto begin() const
    {
        return particles_.begin();
    }
    auto end() const
    {
        return particles_.end();
    }

    template <typename Func>
    void forEachParticle(Func&& func, std::execution::parallel_policy policy = std::execution::par)
    {
        syncSoAToParticlesIfNeeded();
        std::for_each(policy, particles_.begin(), particles_.end(), std::forward<Func>(func));
        invalidateDerivedCaches();
    }

    template <typename Func>
    void forEachParticle(Func&& func,
                         std::execution::parallel_policy policy = std::execution::par) const
    {
        syncSoAToParticlesIfNeeded();
        std::for_each(policy, particles_.begin(), particles_.end(), std::forward<Func>(func));
    }

    double getX(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].position[0];
    }
    double getY(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].position[1];
    }
    double getZ(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].position[2];
    }
    double getVX(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].velocity[0];
    }
    double getVY(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].velocity[1];
    }
    double getVZ(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].velocity[2];
    }
    double getCharge(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].charge;
    }
    double getMass(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].mass;
    }
    double getWeight(size_t i) const
    {
        syncSoAToParticlesIfNeeded();
        return particles_[i].weight;
    }

    void setX(size_t i, double val)
    {
        syncSoAToParticlesIfNeeded();
        particles_[i].position[0] = val;
        invalidateDerivedCaches();
    }
    void setY(size_t i, double val)
    {
        syncSoAToParticlesIfNeeded();
        particles_[i].position[1] = val;
        invalidateDerivedCaches();
    }
    void setZ(size_t i, double val)
    {
        syncSoAToParticlesIfNeeded();
        particles_[i].position[2] = val;
        invalidateDerivedCaches();
    }
    void setVX(size_t i, double val)
    {
        syncSoAToParticlesIfNeeded();
        particles_[i].velocity[0] = val;
        particles_[i].updateEnergy();
        invalidateDerivedCaches();
    }
    void setVY(size_t i, double val)
    {
        syncSoAToParticlesIfNeeded();
        particles_[i].velocity[1] = val;
        particles_[i].updateEnergy();
        invalidateDerivedCaches();
    }
    void setVZ(size_t i, double val)
    {
        syncSoAToParticlesIfNeeded();
        particles_[i].velocity[2] = val;
        particles_[i].updateEnergy();
        invalidateDerivedCaches();
    }

    double* x_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_x_.data();
    }
    double* y_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_y_.data();
    }
    double* z_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_z_.data();
    }
    double* vx_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_vx_.data();
    }
    double* vy_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_vy_.data();
    }
    double* vz_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_vz_.data();
    }
    double* charge_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_charge_.data();
    }
    double* mass_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_mass_.data();
    }
    double* weight_data()
    {
        ensureSoACache();
        soa_cache_modified_ = true;
        stats_dirty_ = true;
        return soa_weight_.data();
    }

    const double* x_data() const
    {
        ensureSoACache();
        return soa_x_.data();
    }
    const double* y_data() const
    {
        ensureSoACache();
        return soa_y_.data();
    }
    const double* z_data() const
    {
        ensureSoACache();
        return soa_z_.data();
    }
    const double* vx_data() const
    {
        ensureSoACache();
        return soa_vx_.data();
    }
    const double* vy_data() const
    {
        ensureSoACache();
        return soa_vy_.data();
    }
    const double* vz_data() const
    {
        ensureSoACache();
        return soa_vz_.data();
    }
    const double* charge_data() const
    {
        ensureSoACache();
        return soa_charge_.data();
    }
    const double* mass_data() const
    {
        ensureSoACache();
        return soa_mass_.data();
    }
    const double* weight_data() const
    {
        ensureSoACache();
        return soa_weight_.data();
    }

    size_t n_electrons() const
    {
        updateStatistics();
        return soa_stats_.electrons;
    }
    size_t n_ions() const
    {
        updateStatistics();
        return soa_stats_.ions;
    }
    const Statistics& getStatistics() const
    {
        updateStatistics();
        return soa_stats_;
    }

    size_t getAlignment() const
    {
        ensureSoACache();
        if (soa_x_.empty())
        {
            return 0;
        }
        return reinterpret_cast<std::uintptr_t>(soa_x_.data()) % 32;
    }

    bool isAligned() const
    {
        ensureSoACache();
        if (soa_x_.empty())
        {
            return true;
        }
        return (reinterpret_cast<std::uintptr_t>(soa_x_.data()) % 32 == 0);
    }

    void borisPushAll(const SCDAT::Geometry::Vector3D& E_field,
                      const SCDAT::Geometry::Vector3D& B_field, double dt);
    void updatePositionsAll(double dt);
    void updateEnergiesAll();

    [[nodiscard]] std::vector<size_t> findActiveParticles() const;
    [[nodiscard]] std::vector<size_t> findParticlesByType(ParticleType type) const;
    [[nodiscard]] std::vector<size_t>
    findParticlesInRegion(const SCDAT::Geometry::Point3D& min_corner,
                          const SCDAT::Geometry::Point3D& max_corner) const;

    [[nodiscard]] size_t getActiveParticleCount() const;
    [[nodiscard]] double getTotalKineticEnergy() const;
    [[nodiscard]] SCDAT::Geometry::Vector3D getTotalMomentum() const;
    [[nodiscard]] SCDAT::Geometry::Point3D getCenterOfMass() const;

    [[nodiscard]] size_t getMemoryUsage() const;
    void compactMemory();

  private:
    Container particles_;
    std::pmr::memory_resource* memory_resource_;

    mutable std::vector<double> soa_x_;
    mutable std::vector<double> soa_y_;
    mutable std::vector<double> soa_z_;
    mutable std::vector<double> soa_vx_;
    mutable std::vector<double> soa_vy_;
    mutable std::vector<double> soa_vz_;
    mutable std::vector<double> soa_charge_;
    mutable std::vector<double> soa_mass_;
    mutable std::vector<double> soa_weight_;
    mutable bool soa_cache_valid_ = false;
    mutable bool soa_cache_modified_ = false;

    mutable Statistics soa_stats_;
    mutable bool stats_dirty_ = true;

    void invalidateDerivedCaches() noexcept
    {
        soa_cache_valid_ = false;
        stats_dirty_ = true;
    }

    void syncSoAToParticlesIfNeeded() const
    {
        if (!soa_cache_modified_)
        {
            return;
        }

        auto& mutable_particles = const_cast<Container&>(particles_);
        const size_t n = mutable_particles.size();
        if (soa_x_.size() != n || soa_y_.size() != n || soa_z_.size() != n || soa_vx_.size() != n ||
            soa_vy_.size() != n || soa_vz_.size() != n || soa_charge_.size() != n ||
            soa_mass_.size() != n || soa_weight_.size() != n)
        {
            soa_cache_modified_ = false;
            soa_cache_valid_ = false;
            stats_dirty_ = true;
            return;
        }

        for (size_t i = 0; i < n; ++i)
        {
            mutable_particles[i].position[0] = soa_x_[i];
            mutable_particles[i].position[1] = soa_y_[i];
            mutable_particles[i].position[2] = soa_z_[i];
            mutable_particles[i].velocity[0] = soa_vx_[i];
            mutable_particles[i].velocity[1] = soa_vy_[i];
            mutable_particles[i].velocity[2] = soa_vz_[i];
            mutable_particles[i].charge = soa_charge_[i];
            mutable_particles[i].mass = soa_mass_[i];
            mutable_particles[i].weight = soa_weight_[i];
            mutable_particles[i].updateEnergy();
        }

        soa_cache_modified_ = false;
        stats_dirty_ = true;
    }

    void ensureSoACache() const
    {
        if (soa_cache_valid_)
        {
            return;
        }

        syncSoAToParticlesIfNeeded();

        const size_t n = particles_.size();
        soa_x_.resize(n);
        soa_y_.resize(n);
        soa_z_.resize(n);
        soa_vx_.resize(n);
        soa_vy_.resize(n);
        soa_vz_.resize(n);
        soa_charge_.resize(n);
        soa_mass_.resize(n);
        soa_weight_.resize(n);

        for (size_t i = 0; i < n; ++i)
        {
            const auto& p = particles_[i];
            soa_x_[i] = p.position[0];
            soa_y_[i] = p.position[1];
            soa_z_[i] = p.position[2];
            soa_vx_[i] = p.velocity[0];
            soa_vy_[i] = p.velocity[1];
            soa_vz_[i] = p.velocity[2];
            soa_charge_[i] = p.charge;
            soa_mass_[i] = p.mass;
            soa_weight_[i] = p.weight;
        }

        soa_cache_valid_ = true;
    }

    void updateStatistics() const
    {
        if (!stats_dirty_)
        {
            return;
        }

        syncSoAToParticlesIfNeeded();

        size_t electrons = 0;
        size_t ions = 0;
        for (const auto& particle : particles_)
        {
            if (particle.charge < 0.0)
            {
                ++electrons;
            }
            else if (particle.charge > 0.0)
            {
                ++ions;
            }
        }

        soa_stats_.update(particles_.size(), particles_.capacity(), electrons, ions);
        stats_dirty_ = false;
    }

    void simdBorisPush(std::span<ParticleData> particles, const SCDAT::Geometry::Vector3D& E_field,
                       const SCDAT::Geometry::Vector3D& B_field, double dt);
};

template <typename ParticleAoS>
ParticleDataManager convertToSoA(const std::vector<ParticleAoS>& particles)
{
    ParticleDataManager soa;
    soa.reserve(particles.size());

    for (const auto& p : particles)
    {
        soa.push(p.x, p.y, p.z, p.vx, p.vy, p.vz, p.charge, p.mass, p.weight);
    }

    return soa;
}

inline void printStatistics(const ParticleDataManager& particles, const std::string& label = "")
{
    const auto& stats = particles.getStatistics();

    if (!label.empty())
    {
        std::cout << "\n[" << label << " - SoA统计]\n";
    }
    else
    {
        std::cout << "\n[SoA统计]\n";
    }

    std::cout << "  总粒子数: " << stats.total_particles << "\n";
    std::cout << "  电子数:   " << stats.electrons << "\n";
    std::cout << "  离子数:   " << stats.ions << "\n";
    std::cout << "  容量:     " << stats.capacity << "\n";
    std::cout << "  内存:     " << std::fixed << std::setprecision(2) << stats.memory_mb << " MB\n";
    std::cout << "  对齐:     " << (particles.isAligned() ? "是(32B)" : "否") << "\n";
}

/**
 * @brief 粒子迭代器类
 *
 * 支持类型过滤和空间区域过滤的粒子迭代器
 */
class ParticleIterator
{
  public:
    /**
     * @brief 构造函数
     * @param container 粒子容器
     * @param type_filter 类型过滤器 (可选)
     * @param region_filter 空间区域过滤器 (可选)
     */
    ParticleIterator(ParticleContainer& container,
                     std::optional<ParticleType> type_filter = std::nullopt,
                     std::function<bool(const SCDAT::Geometry::Point3D&)> region_filter = nullptr);

    // 迭代器操作
    bool hasNext() const;
    ParticlePtr next();
    void reset();

    // 统计信息
    std::size_t count() const;

  private:
    ParticleContainer& container_;
    std::optional<ParticleType> type_filter_;
    std::function<bool(const SCDAT::Geometry::Point3D&)> region_filter_;
    std::size_t current_index_;

    bool passesFilter(const ParticlePtr& particle) const;
};

/**
 * @brief 粒子管理器类
 *
 * 高级粒子管理功能，包括生命周期管理、统计收集等
 */
class ParticleManager
{
  public:
    /**
     * @brief 构造函数
     * @param max_particles 最大粒子数量
     */
    explicit ParticleManager(std::size_t max_particles = 1000000);

    /**
     * @brief 析构函数
     */
    ~ParticleManager() = default;

    // 粒子创建和管理
    ParticleId createElectron(const SCDAT::Geometry::Point3D& position,
                              const SCDAT::Geometry::Vector3D& velocity, double weight = 1.0);

    ParticleId createIon(const SCDAT::Geometry::Point3D& position,
                         const SCDAT::Geometry::Vector3D& velocity, int mass_number,
                         int charge_number, double weight = 1.0);

    ParticleId createPhotoelectron(const SCDAT::Geometry::Point3D& position,
                                   const SCDAT::Geometry::Vector3D& velocity, double photon_energy,
                                   double work_function, double weight = 1.0);

    // 粒子访问
    ParticlePtr getParticle(ParticleId id) const;
    std::vector<ParticlePtr> getParticlesByType(ParticleType type) const;

    // 批量操作
    void updateParticleAges(double dt);
    std::size_t removeOldParticles(double max_age);
    std::size_t removeInactiveParticles();

    // 统计信息
    struct Statistics
    {
        std::size_t total_particles = 0;
        std::size_t active_particles = 0;
        std::unordered_map<ParticleType, std::size_t> type_counts;
        std::unordered_map<ParticleStatus, std::size_t> status_counts;
        double total_kinetic_energy = 0.0;
        double total_charge = 0.0;
        double average_age = 0.0;
    };

    Statistics getStatistics() const;
    void printStatistics() const;

    // 容器访问
    const ParticleContainer& getContainer() const
    {
        return *container_;
    }
    ParticleContainer& getContainer()
    {
        return *container_;
    }

    // 迭代器创建
    ParticleIterator createIterator(
        std::optional<ParticleType> type_filter = std::nullopt,
        std::function<bool(const SCDAT::Geometry::Point3D&)> region_filter = nullptr) const;

    // 内存管理
    void reserveCapacity(std::size_t capacity);
    void optimizeMemory();

    // 线程安全
    void setThreadSafe(bool thread_safe)
    {
        thread_safe_ = thread_safe;
    }
    bool isThreadSafe() const
    {
        return thread_safe_;
    }

    // 碰撞检测
    void enableCollisionDetection(bool enable = true);
    bool isCollisionDetectionEnabled() const
    {
        return collision_detection_enabled_;
    }
    void processCollisions(double dt);

    // 获取所有粒子（用于碰撞检测）
    std::vector<ParticlePtr> getAllParticles() const;

  private:
    ParticleContainerPtr container_;   ///< 粒子容器
    std::size_t max_particles_;        ///< 最大粒子数量
    std::atomic<ParticleId> next_id_;  ///< 下一个可用ID
    bool thread_safe_;                 ///< 线程安全标志
    bool collision_detection_enabled_; ///< 碰撞检测启用标志

    mutable std::mutex stats_mutex_; ///< 统计信息互斥锁

    // 内部辅助方法
    ParticleId getNextId()
    {
        return ++next_id_;
    }
    void enforceParticleLimit();
};

/**
 * @brief 粒子内存池类
 *
 * 高效的粒子内存管理，减少内存分配开销
 */
class ParticleMemoryPool
{
  public:
    /**
     * @brief 构造函数
     * @param block_size 内存块大小
     * @param initial_blocks 初始块数量
     */
    ParticleMemoryPool(std::size_t block_size = 1024, std::size_t initial_blocks = 10);

    /**
     * @brief 析构函数
     */
    ~ParticleMemoryPool();

    // 内存分配和释放
    template <typename T, typename... Args> std::shared_ptr<T> allocate(Args&&... args);

    void deallocate(void* ptr);

    // 统计信息
    std::size_t getTotalAllocated() const
    {
        return total_allocated_;
    }
    std::size_t getTotalDeallocated() const
    {
        return total_deallocated_;
    }
    std::size_t getActiveAllocations() const
    {
        return total_allocated_ - total_deallocated_;
    }

  private:
    struct MemoryBlock
    {
        void* data;
        std::size_t size;
        std::size_t used;
    };

    std::vector<MemoryBlock> blocks_;
    std::size_t block_size_;
    std::size_t current_block_;

    std::atomic<std::size_t> total_allocated_;
    std::atomic<std::size_t> total_deallocated_;

    std::mutex pool_mutex_;

    // 内部方法
    void allocateNewBlock();
    void* allocateFromBlock(std::size_t size);
};

// 模板方法实现
template <typename T, typename... Args>
std::shared_ptr<T> ParticleMemoryPool::allocate(Args&&... args)
{
    std::lock_guard<std::mutex> lock(pool_mutex_);

    void* ptr = allocateFromBlock(sizeof(T));
    if (!ptr)
    {
        allocateNewBlock();
        ptr = allocateFromBlock(sizeof(T));
    }

    if (ptr)
    {
        T* obj = new (ptr) T(std::forward<Args>(args)...);
        ++total_allocated_;

        return std::shared_ptr<T>(obj,
                                  [this](T* p)
                                  {
                                      p->~T();
                                      this->deallocate(p);
                                  });
    }

    return nullptr;
}

} // namespace Particle
} // namespace SCDAT

#endif // SCDAT_PARTICLE_MANAGER_H
