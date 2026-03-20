/**
 * @file ParticleManager.cpp
 * @brief 粒子管理器实现
 *
 * 实现了高效的粒子容器和管理器，支持百万级粒子管理
 *
 * @author Wang Sizhan
 * @date 2026.3.19 16:02:02
 * @version 0.0.1
 */

#include "../include/ParticleManager.h"
#include "../include/ParticlePusher.h"
#include <algorithm>
#include <execution>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <shared_mutex>

namespace SCDAT
{
namespace Particle
{

namespace
{

BoundaryCondition makeUnboundedBoundary()
{
    BoundaryCondition boundary;
    boundary.min_bounds = SCDAT::Geometry::Point3D(std::numeric_limits<double>::lowest(),
                                                   std::numeric_limits<double>::lowest(),
                                                   std::numeric_limits<double>::lowest());
    boundary.max_bounds = SCDAT::Geometry::Point3D(std::numeric_limits<double>::max(),
                                                   std::numeric_limits<double>::max(),
                                                   std::numeric_limits<double>::max());
    return boundary;
}

void applyBorisPushToData(ParticleData& particle_data, const SCDAT::Geometry::Vector3D& E_field,
                          const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    if (particle_data.status != ParticleStatus::ACTIVE || particle_data.mass <= 0.0)
    {
        return;
    }

    const double old_age = particle_data.age;
    const double old_time = particle_data.last_update_time;
    const ParticleStatus old_status = particle_data.status;

    Particle particle(particle_data);
    BorisAlgorithm boris_pusher(dt);
    boris_pusher.setBoundaryCondition(makeUnboundedBoundary());

    const FieldFunction electric = [&E_field](const SCDAT::Geometry::Point3D&) { return E_field; };
    const FieldFunction magnetic = [&B_field](const SCDAT::Geometry::Point3D&) { return B_field; };
    boris_pusher.pushParticle(&particle, electric, magnetic);

    particle_data = particle.getData();
    particle_data.age = old_age;
    particle_data.last_update_time = old_time;
    particle_data.status = old_status;
}

} // namespace

//==============================================================================
// ParticleContainer 实现
//==============================================================================

ParticleContainer::ParticleContainer(std::size_t initial_capacity)
{
    particles_.reserve(initial_capacity);
    id_map_.reserve(initial_capacity);
}

void ParticleContainer::addParticle(const Particle& particle)
{
    std::unique_lock<std::shared_mutex> lock(mutex_);

    ParticleId id = particle.getId();

    // 检查是否已存在
    if (id_map_.find(id) != id_map_.end())
    {
        return; // 粒子已存在
    }

    // 添加粒子
    std::size_t index = particles_.size();
    particles_.push_back(particle);
    id_map_[id] = index;
}

bool ParticleContainer::removeParticle(ParticleId id)
{
    std::unique_lock<std::shared_mutex> lock(mutex_);

    auto it = id_map_.find(id);
    if (it == id_map_.end())
    {
        return false; // 粒子不存在
    }

    std::size_t index = it->second;

    // 使用swap-and-pop技术快速删除
    if (index < particles_.size() - 1)
    {
        std::swap(particles_[index], particles_.back());
        // 更新被交换粒子的索引
        ParticleId swapped_id = particles_[index].getId();
        id_map_[swapped_id] = index;
    }

    particles_.pop_back();
    id_map_.erase(it);

    return true;
}

ParticlePtr ParticleContainer::getParticle(ParticleId id) const
{
    std::shared_lock<std::shared_mutex> lock(mutex_);

    auto it = id_map_.find(id);
    if (it != id_map_.end() && it->second < particles_.size())
    {

        return const_cast<Particle*>(&particles_[it->second]);
    }
    return nullptr;
}

void ParticleContainer::addParticles(const std::vector<Particle>& particles)
{
    std::unique_lock<std::shared_mutex> lock(mutex_);

    particles_.reserve(particles_.size() + particles.size());

    for (const auto& particle : particles)
    {
        ParticleId id = particle.getId();
        if (id_map_.find(id) == id_map_.end())
        {
            std::size_t index = particles_.size();
            particles_.push_back(particle);
            id_map_[id] = index;
        }
    }
}

std::size_t ParticleContainer::removeInactiveParticles()
{
    std::unique_lock<std::shared_mutex> lock(mutex_);

    std::size_t removed_count = 0;

    // 使用partition算法将活跃粒子移到前面
    auto partition_point = std::partition(particles_.begin(), particles_.end(),
                                          [](const Particle& p) { return p.isActive(); });

    // 计算移除的粒子数量
    removed_count = static_cast<std::size_t>(std::distance(partition_point, particles_.end()));

    // 移除非活跃粒子
    particles_.erase(partition_point, particles_.end());

    // 重建ID映射
    rebuildIdMap();

    return removed_count;
}

std::vector<ParticlePtr> ParticleContainer::getParticlesByType(ParticleType type) const
{
    std::shared_lock<std::shared_mutex> lock(mutex_);

    std::vector<ParticlePtr> result;
    result.reserve(particles_.size() / 4); // 预估容量

    for (const auto& particle : particles_)
    {
        if (particle.getType() == type)
        {
            result.push_back(const_cast<Particle*>(&particle));
        }
    }

    return result;
}

std::size_t ParticleContainer::getActiveCount() const
{
    std::shared_lock<std::shared_mutex> lock(mutex_);

    return static_cast<std::size_t>(std::count_if(particles_.begin(), particles_.end(),
                                                  [](const Particle& p) { return p.isActive(); }));
}

std::size_t ParticleContainer::getTypeCount(ParticleType type) const
{
    std::shared_lock<std::shared_mutex> lock(mutex_);

    return static_cast<std::size_t>(std::count_if(particles_.begin(), particles_.end(),
                                                  [type](const Particle& p)
                                                  { return p.getType() == type; }));
}

std::unordered_map<ParticleType, std::size_t> ParticleContainer::getTypeCounts() const
{
    std::shared_lock<std::shared_mutex> lock(mutex_);

    std::unordered_map<ParticleType, std::size_t> counts;

    for (const auto& particle : particles_)
    {
        counts[particle.getType()]++;
    }

    return counts;
}

void ParticleContainer::reserve(std::size_t capacity)
{
    std::unique_lock<std::shared_mutex> lock(mutex_);
    particles_.reserve(capacity);
    id_map_.reserve(capacity);
}

void ParticleContainer::shrinkToFit()
{
    std::unique_lock<std::shared_mutex> lock(mutex_);
    particles_.shrink_to_fit();
}

void ParticleContainer::clear()
{
    std::unique_lock<std::shared_mutex> lock(mutex_);
    particles_.clear();
    id_map_.clear();
}

void ParticleContainer::rebuildIdMap()
{
    id_map_.clear();
    id_map_.reserve(particles_.size());

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        id_map_[particles_[i].getId()] = i;
    }
}

void ParticleContainer::compactStorage()
{
    // 对于值类型存储，不需要移除空指针
    // 但可以移除无效粒子
    removeInactiveParticles();
}

//==============================================================================
// ParticleDataManager 实现
//==============================================================================

void ParticleDataManager::borisPushAll(const SCDAT::Geometry::Vector3D& E_field,
                                       const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    syncSoAToParticlesIfNeeded();

    if (particles_.empty())
    {
        return;
    }

#ifdef __AVX2__
    if (particles_.size() >= 4)
    {
        simdBorisPush(std::span(particles_), E_field, B_field, dt);
        invalidateDerivedCaches();
        return;
    }
#endif

    for (auto& particle : particles_)
    {
        applyBorisPushToData(particle, E_field, B_field, dt);
    }

    invalidateDerivedCaches();
}

void ParticleDataManager::updatePositionsAll(double dt)
{
    syncSoAToParticlesIfNeeded();

    std::for_each(std::execution::par_unseq, particles_.begin(), particles_.end(),
                  [dt](ParticleData& particle)
                  {
                      if (particle.status == ParticleStatus::ACTIVE)
                      {
                          particle.position[0] += particle.velocity[0] * dt;
                          particle.position[1] += particle.velocity[1] * dt;
                          particle.position[2] += particle.velocity[2] * dt;
                      }
                  });

    invalidateDerivedCaches();
}

void ParticleDataManager::updateEnergiesAll()
{
    syncSoAToParticlesIfNeeded();

    std::for_each(std::execution::par_unseq, particles_.begin(), particles_.end(),
                  [](ParticleData& particle) { particle.updateEnergy(); });

    invalidateDerivedCaches();
}

std::vector<size_t> ParticleDataManager::findActiveParticles() const
{
    syncSoAToParticlesIfNeeded();

    std::vector<size_t> active_indices;

    for (size_t i = 0; i < particles_.size(); ++i)
    {
        if (particles_[i].status == ParticleStatus::ACTIVE)
        {
            active_indices.push_back(i);
        }
    }

    return active_indices;
}

std::vector<size_t> ParticleDataManager::findParticlesByType(ParticleType type) const
{
    syncSoAToParticlesIfNeeded();

    std::vector<size_t> type_indices;

    for (size_t i = 0; i < particles_.size(); ++i)
    {
        if (particles_[i].type == type)
        {
            type_indices.push_back(i);
        }
    }

    return type_indices;
}

size_t ParticleDataManager::getActiveParticleCount() const
{
    syncSoAToParticlesIfNeeded();

    return static_cast<size_t>(std::count_if(std::execution::par_unseq, particles_.begin(),
                                             particles_.end(), [](const ParticleData& p)
                                             { return p.status == ParticleStatus::ACTIVE; }));
}

double ParticleDataManager::getTotalKineticEnergy() const
{
    syncSoAToParticlesIfNeeded();

    return std::transform_reduce(std::execution::par_unseq, particles_.begin(), particles_.end(),
                                 0.0, std::plus<>{}, [](const ParticleData& p)
                                 { return (p.status == ParticleStatus::ACTIVE) ? p.energy : 0.0; });
}

size_t ParticleDataManager::getMemoryUsage() const
{
    syncSoAToParticlesIfNeeded();

    return particles_.size() * sizeof(ParticleData) + particles_.capacity() * sizeof(ParticleData);
}

std::vector<size_t>
ParticleDataManager::findParticlesInRegion(const SCDAT::Geometry::Point3D& min_corner,
                                           const SCDAT::Geometry::Point3D& max_corner) const
{
    syncSoAToParticlesIfNeeded();

    std::vector<size_t> region_indices;

    for (size_t i = 0; i < particles_.size(); ++i)
    {
        const auto& p = particles_[i];
        const bool inside_x = p.position[0] >= min_corner.x() && p.position[0] <= max_corner.x();
        const bool inside_y = p.position[1] >= min_corner.y() && p.position[1] <= max_corner.y();
        const bool inside_z = p.position[2] >= min_corner.z() && p.position[2] <= max_corner.z();

        if (inside_x && inside_y && inside_z)
        {
            region_indices.push_back(i);
        }
    }

    return region_indices;
}

SCDAT::Geometry::Vector3D ParticleDataManager::getTotalMomentum() const
{
    syncSoAToParticlesIfNeeded();

    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;

    for (const auto& p : particles_)
    {
        if (p.status != ParticleStatus::ACTIVE || p.mass <= 0.0)
        {
            continue;
        }
        px += p.mass * p.velocity[0];
        py += p.mass * p.velocity[1];
        pz += p.mass * p.velocity[2];
    }

    return SCDAT::Geometry::Vector3D(px, py, pz);
}

SCDAT::Geometry::Point3D ParticleDataManager::getCenterOfMass() const
{
    syncSoAToParticlesIfNeeded();

    double weighted_x = 0.0;
    double weighted_y = 0.0;
    double weighted_z = 0.0;
    double total_mass = 0.0;

    for (const auto& p : particles_)
    {
        if (p.status != ParticleStatus::ACTIVE || p.mass <= 0.0)
        {
            continue;
        }
        weighted_x += p.mass * p.position[0];
        weighted_y += p.mass * p.position[1];
        weighted_z += p.mass * p.position[2];
        total_mass += p.mass;
    }

    if (total_mass <= 0.0)
    {
        return SCDAT::Geometry::Point3D(0.0, 0.0, 0.0);
    }

    return SCDAT::Geometry::Point3D(weighted_x / total_mass, weighted_y / total_mass,
                                    weighted_z / total_mass);
}

void ParticleDataManager::compactMemory()
{
    syncSoAToParticlesIfNeeded();

    auto new_end = std::remove_if(particles_.begin(), particles_.end(), [](const ParticleData& p)
                                  { return p.status != ParticleStatus::ACTIVE; });
    particles_.erase(new_end, particles_.end());
    particles_.shrink_to_fit();

    invalidateDerivedCaches();
}

#ifdef __AVX2__
#include <immintrin.h>

void ParticleDataManager::simdBorisPush(std::span<ParticleData> particles,
                                        const SCDAT::Geometry::Vector3D& E_field,
                                        const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    const size_t simd_width = 4;
    const size_t aligned_size = (particles.size() / simd_width) * simd_width;

    const double half_dt = 0.5 * dt;

    __m256d E_x = _mm256_set1_pd(E_field.x());
    __m256d E_y = _mm256_set1_pd(E_field.y());
    __m256d E_z = _mm256_set1_pd(E_field.z());

    __m256d B_x = _mm256_set1_pd(B_field.x());
    __m256d B_y = _mm256_set1_pd(B_field.y());
    __m256d B_z = _mm256_set1_pd(B_field.z());

    __m256d half_dt_vec = _mm256_set1_pd(half_dt);
    __m256d dt_vec = _mm256_set1_pd(dt);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d two = _mm256_set1_pd(2.0);

    for (size_t i = 0; i < aligned_size; i += simd_width)
    {
        __m256i active_mask = _mm256_setzero_si256();
        __m256d q_over_m_vec = _mm256_setzero_pd();

        __m256d pos_x = _mm256_setzero_pd();
        __m256d pos_y = _mm256_setzero_pd();
        __m256d pos_z = _mm256_setzero_pd();

        __m256d vel_x = _mm256_setzero_pd();
        __m256d vel_y = _mm256_setzero_pd();
        __m256d vel_z = _mm256_setzero_pd();

        for (size_t j = 0; j < simd_width; ++j)
        {
            const auto& particle = particles[i + j];

            if (particle.status == ParticleStatus::ACTIVE && particle.mass > 0.0)
            {
                reinterpret_cast<int64_t*>(&active_mask)[j] = -1;
                reinterpret_cast<double*>(&q_over_m_vec)[j] = particle.charge / particle.mass;

                reinterpret_cast<double*>(&pos_x)[j] = particle.position[0];
                reinterpret_cast<double*>(&pos_y)[j] = particle.position[1];
                reinterpret_cast<double*>(&pos_z)[j] = particle.position[2];

                reinterpret_cast<double*>(&vel_x)[j] = particle.velocity[0];
                reinterpret_cast<double*>(&vel_y)[j] = particle.velocity[1];
                reinterpret_cast<double*>(&vel_z)[j] = particle.velocity[2];
            }
        }

        __m256d E_accel_x = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, E_x), half_dt_vec);
        __m256d E_accel_y = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, E_y), half_dt_vec);
        __m256d E_accel_z = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, E_z), half_dt_vec);

        __m256d v_minus_x = _mm256_add_pd(vel_x, E_accel_x);
        __m256d v_minus_y = _mm256_add_pd(vel_y, E_accel_y);
        __m256d v_minus_z = _mm256_add_pd(vel_z, E_accel_z);

        __m256d t_x = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, B_x), half_dt_vec);
        __m256d t_y = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, B_y), half_dt_vec);
        __m256d t_z = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, B_z), half_dt_vec);

        __m256d t_mag_sq =
            _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(t_x, t_x), _mm256_mul_pd(t_y, t_y)),
                          _mm256_mul_pd(t_z, t_z));
        __m256d denominator = _mm256_add_pd(one, t_mag_sq);
        __m256d s_factor = _mm256_div_pd(two, denominator);

        __m256d s_x = _mm256_mul_pd(t_x, s_factor);
        __m256d s_y = _mm256_mul_pd(t_y, s_factor);
        __m256d s_z = _mm256_mul_pd(t_z, s_factor);

        __m256d cross1_x =
            _mm256_sub_pd(_mm256_mul_pd(v_minus_y, t_z), _mm256_mul_pd(v_minus_z, t_y));
        __m256d cross1_y =
            _mm256_sub_pd(_mm256_mul_pd(v_minus_z, t_x), _mm256_mul_pd(v_minus_x, t_z));
        __m256d cross1_z =
            _mm256_sub_pd(_mm256_mul_pd(v_minus_x, t_y), _mm256_mul_pd(v_minus_y, t_x));

        __m256d v_prime_x = _mm256_add_pd(v_minus_x, cross1_x);
        __m256d v_prime_y = _mm256_add_pd(v_minus_y, cross1_y);
        __m256d v_prime_z = _mm256_add_pd(v_minus_z, cross1_z);

        __m256d cross2_x =
            _mm256_sub_pd(_mm256_mul_pd(v_prime_y, s_z), _mm256_mul_pd(v_prime_z, s_y));
        __m256d cross2_y =
            _mm256_sub_pd(_mm256_mul_pd(v_prime_z, s_x), _mm256_mul_pd(v_prime_x, s_z));
        __m256d cross2_z =
            _mm256_sub_pd(_mm256_mul_pd(v_prime_x, s_y), _mm256_mul_pd(v_prime_y, s_x));

        __m256d v_plus_x = _mm256_add_pd(v_minus_x, cross2_x);
        __m256d v_plus_y = _mm256_add_pd(v_minus_y, cross2_y);
        __m256d v_plus_z = _mm256_add_pd(v_minus_z, cross2_z);

        __m256d new_vel_x = _mm256_add_pd(v_plus_x, E_accel_x);
        __m256d new_vel_y = _mm256_add_pd(v_plus_y, E_accel_y);
        __m256d new_vel_z = _mm256_add_pd(v_plus_z, E_accel_z);

        __m256d new_pos_x = _mm256_add_pd(pos_x, _mm256_mul_pd(new_vel_x, dt_vec));
        __m256d new_pos_y = _mm256_add_pd(pos_y, _mm256_mul_pd(new_vel_y, dt_vec));
        __m256d new_pos_z = _mm256_add_pd(pos_z, _mm256_mul_pd(new_vel_z, dt_vec));

        for (size_t j = 0; j < simd_width; ++j)
        {
            if (reinterpret_cast<int64_t*>(&active_mask)[j] != 0)
            {
                auto& particle = particles[i + j];

                particle.velocity[0] = reinterpret_cast<double*>(&new_vel_x)[j];
                particle.velocity[1] = reinterpret_cast<double*>(&new_vel_y)[j];
                particle.velocity[2] = reinterpret_cast<double*>(&new_vel_z)[j];

                particle.position[0] = reinterpret_cast<double*>(&new_pos_x)[j];
                particle.position[1] = reinterpret_cast<double*>(&new_pos_y)[j];
                particle.position[2] = reinterpret_cast<double*>(&new_pos_z)[j];
                particle.updateEnergy();
            }
        }
    }

    for (size_t i = aligned_size; i < particles.size(); ++i)
    {
        applyBorisPushToData(particles[i], E_field, B_field, dt);
    }
}
#endif

//==============================================================================
// ParticleIterator 实现
//==============================================================================

ParticleIterator::ParticleIterator(
    ParticleContainer& container, std::optional<ParticleType> type_filter,
    std::function<bool(const SCDAT::Geometry::Point3D&)> region_filter)
    : container_(container), type_filter_(type_filter), region_filter_(region_filter),
      current_index_(0)
{
}

bool ParticleIterator::hasNext() const
{
    container_.lockForRead();

    for (std::size_t i = current_index_; i < container_.size(); ++i)
    {
        auto& particle = container_.particles_[i]; // Non-const reference
        if (passesFilter(&particle))
        {
            container_.unlockForRead();
            return true;
        }
    }

    container_.unlockForRead();
    return false;
}

ParticlePtr ParticleIterator::next()
{
    container_.lockForRead();

    while (current_index_ < container_.size())
    {
        auto& particle = container_.particles_[current_index_++];
        if (passesFilter(&particle))
        {
            container_.unlockForRead();
            return &particle;
        }
    }

    container_.unlockForRead();
    return nullptr;
}

void ParticleIterator::reset()
{
    current_index_ = 0;
}

std::size_t ParticleIterator::count() const
{
    container_.lockForRead();

    std::size_t count = 0;
    for (auto& particle : container_.particles_)
    {
        if (passesFilter(&particle))
        {
            count++;
        }
    }

    container_.unlockForRead();
    return count;
}

bool ParticleIterator::passesFilter(const ParticlePtr& particle) const
{
    if (!particle)
        return false;

    // 类型过滤
    if (type_filter_ && particle->getType() != *type_filter_)
    {
        return false;
    }

    // 空间区域过滤
    if (region_filter_ && !region_filter_(particle->getPosition()))
    {
        return false;
    }

    return true;
}

//==============================================================================
// ParticleManager 实现
//==============================================================================

ParticleManager::ParticleManager(std::size_t max_particles)
    : container_(std::make_shared<ParticleContainer>(max_particles / 10)),
      max_particles_(max_particles), next_id_(0), thread_safe_(false),
      collision_detection_enabled_(false)
{
}

ParticleId ParticleManager::createElectron(const SCDAT::Geometry::Point3D& position,
                                           const SCDAT::Geometry::Vector3D& velocity, double weight)
{
    ParticleId id = getNextId();
    auto electron = ParticleFactory::createElectron(id, position, velocity, weight);
    container_->addParticle(electron);

    enforceParticleLimit();
    return id;
}

ParticleId ParticleManager::createIon(const SCDAT::Geometry::Point3D& position,
                                      const SCDAT::Geometry::Vector3D& velocity, int mass_number,
                                      int charge_number, double weight)
{
    ParticleId id = getNextId();
    auto ion =
        ParticleFactory::createIon(id, position, velocity, mass_number, charge_number, weight);
    container_->addParticle(ion);

    enforceParticleLimit();
    return id;
}

ParticleId ParticleManager::createPhotoelectron(const SCDAT::Geometry::Point3D& position,
                                                const SCDAT::Geometry::Vector3D& velocity,
                                                double photon_energy, double work_function,
                                                double weight)
{
    ParticleId id = getNextId();
    auto photoelectron = ParticleFactory::createPhotoelectron(id, position, velocity, photon_energy,
                                                              work_function, weight);
    container_->addParticle(photoelectron);

    enforceParticleLimit();
    return id;
}

ParticlePtr ParticleManager::getParticle(ParticleId id) const
{
    return container_->getParticle(id);
}

std::vector<ParticlePtr> ParticleManager::getParticlesByType(ParticleType type) const
{
    return container_->getParticlesByType(type);
}

void ParticleManager::updateParticleAges(double dt)
{
    if (thread_safe_)
        container_->lockForWrite();

    for (auto& particle : *container_)
    {
        if (particle.isActive())
        {
            particle.updateAge(dt);
        }
    }

    if (thread_safe_)
        container_->unlockForWrite();
}

std::size_t ParticleManager::removeOldParticles(double max_age)
{
    if (thread_safe_)
        container_->lockForWrite();

    std::size_t removed_count = 0;

    for (auto& particle : *container_)
    {
        if (particle.getAge() > max_age)
        {
            particle.setStatus(ParticleStatus::INVALID);
            removed_count++;
        }
    }

    // 移除标记为无效的粒子
    removed_count = container_->removeInactiveParticles();

    if (thread_safe_)
        container_->unlockForWrite();

    return removed_count;
}

std::size_t ParticleManager::removeInactiveParticles()
{
    return container_->removeInactiveParticles();
}

void ParticleManager::enforceParticleLimit()
{
    if (container_->size() > max_particles_)
    {
        // 移除最老的粒子
        std::vector<std::pair<double, ParticlePtr>> aged_particles;

        for (auto& particle : *container_)
        {
            aged_particles.emplace_back(particle.getAge(), &particle);
        }

        // 按年龄排序
        std::sort(aged_particles.begin(), aged_particles.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        // 移除超出限制的粒子
        std::size_t to_remove = container_->size() - max_particles_;
        for (std::size_t i = 0; i < to_remove && i < aged_particles.size(); ++i)
        {
            aged_particles[i].second->setStatus(ParticleStatus::INVALID);
        }

        container_->removeInactiveParticles();
    }
}

ParticleManager::Statistics ParticleManager::getStatistics() const
{
    std::lock_guard<std::mutex> lock(stats_mutex_);

    Statistics stats;
    stats.total_particles = container_->size();
    stats.active_particles = container_->getActiveCount();
    stats.type_counts = container_->getTypeCounts();

    // 计算状态统计
    for (const auto& particle : *container_)
    {
        stats.status_counts[particle.getStatus()]++;
        stats.total_kinetic_energy += particle.getKineticEnergy();
        stats.total_charge += particle.getCharge();
    }

    // 计算平均年龄
    if (stats.active_particles > 0)
    {
        double total_age = 0.0;
        std::size_t active_count = 0;

        for (const auto& particle : *container_)
        {
            if (particle.isActive())
            {
                total_age += particle.getAge();
                active_count++;
            }
        }

        stats.average_age = total_age / static_cast<double>(active_count);
    }

    return stats;
}

void ParticleManager::printStatistics() const
{
    auto stats = getStatistics();

    std::cout << "\n=== 粒子系统统计信息 ===\n";
    std::cout << "总粒子数: " << stats.total_particles << "\n";
    std::cout << "活跃粒子数: " << stats.active_particles << "\n";
    std::cout << "平均年龄: " << std::fixed << std::setprecision(6) << stats.average_age << " s\n";
    std::cout << "总动能: " << std::scientific << std::setprecision(3) << stats.total_kinetic_energy
              << " J\n";
    std::cout << "总电荷: " << std::scientific << std::setprecision(3) << stats.total_charge
              << " C\n";

    std::cout << "\n类型分布:\n";
    for (const auto& [type, count] : stats.type_counts)
    {
        std::cout << "  类型 " << static_cast<int>(type) << ": " << count << "\n";
    }

    std::cout << "\n状态分布:\n";
    for (const auto& [status, count] : stats.status_counts)
    {
        std::cout << "  状态 " << static_cast<int>(status) << ": " << count << "\n";
    }
    std::cout << "========================\n\n";
}

ParticleIterator ParticleManager::createIterator(
    std::optional<ParticleType> type_filter,
    std::function<bool(const SCDAT::Geometry::Point3D&)> region_filter) const
{
    return ParticleIterator(*container_, type_filter, region_filter);
}

void ParticleManager::reserveCapacity(std::size_t capacity)
{
    container_->reserve(capacity);
}

void ParticleManager::optimizeMemory()
{
    container_->removeInactiveParticles();
    container_->shrinkToFit();
}

void ParticleManager::enableCollisionDetection(bool enable)
{
    collision_detection_enabled_ = enable;
}

void ParticleManager::processCollisions(double dt)
{
    (void)dt; // 避免未使用参数警告
    if (!collision_detection_enabled_)
        return;

    // 简化的碰撞检测实现
    // 实际实现需要更复杂的空间分割和碰撞检测算法
    auto particles = getAllParticles();

    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        for (std::size_t j = i + 1; j < particles.size(); ++j)
        {
            if (particles[i] && particles[j] && particles[i]->isActive() &&
                particles[j]->isActive())
            {

                // 简单的距离检测
                auto pos1 = particles[i]->getPosition();
                auto pos2 = particles[j]->getPosition();
                double distance = (pos1 - pos2).magnitude();

                // 如果距离小于某个阈值，认为发生碰撞
                double collision_threshold = 1e-9; // 1 nm
                if (distance < collision_threshold)
                {
                    // 处理碰撞（简化实现）
                    // 实际应该根据物理模型计算碰撞结果
                }
            }
        }
    }
}

std::vector<ParticlePtr> ParticleManager::getAllParticles() const
{
    std::vector<ParticlePtr> particles;
    particles.reserve(container_->size());

    for (auto& particle : *container_)
    {
        if (particle.isActive())
        {
            particles.push_back(&particle);
        }
    }

    return particles;
}

//==============================================================================
// ParticleMemoryPool 实现
//==============================================================================

ParticleMemoryPool::ParticleMemoryPool(std::size_t block_size, std::size_t initial_blocks)
    : block_size_(block_size), current_block_(0), total_allocated_(0), total_deallocated_(0)
{

    blocks_.reserve(initial_blocks * 2);

    // 分配初始内存块
    for (std::size_t i = 0; i < initial_blocks; ++i)
    {
        allocateNewBlock();
    }
}

ParticleMemoryPool::~ParticleMemoryPool()
{
    for (auto& block : blocks_)
    {
        std::free(block.data);
    }
}

void ParticleMemoryPool::deallocate(void* ptr)
{
    if (ptr)
    {
        ++total_deallocated_;
        // 注意：这里简化实现，实际应该将内存返回到池中
        // 完整实现需要跟踪每个分配的内存块
    }
}

void ParticleMemoryPool::allocateNewBlock()
{
    MemoryBlock block;
    block.data = std::malloc(block_size_);
    block.size = block_size_;
    block.used = 0;

    if (block.data)
    {
        blocks_.push_back(block);
    }
}

void* ParticleMemoryPool::allocateFromBlock(std::size_t size)
{
    // 对齐到8字节边界
    size = (size + 7) & ~static_cast<std::size_t>(7);

    for (std::size_t i = current_block_; i < blocks_.size(); ++i)
    {
        auto& block = blocks_[i];
        if (block.used + size <= block.size)
        {
            void* ptr = static_cast<char*>(block.data) + block.used;
            block.used += size;
            current_block_ = i;
            return ptr;
        }
    }

    return nullptr;
}

} // namespace Particle
} // namespace SCDAT
