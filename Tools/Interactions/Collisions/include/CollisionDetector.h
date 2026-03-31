/**
 * @file CollisionDetector.h
 * @brief 粒子碰撞检测（头文件实现）
 * @details
 * 该文件提供“检测 + 处理 + 统计”一体化碰撞器：
 * - 检测阶段：基于粒子对遍历、距离门限、截面模型计算概率；
 * - 处理阶段：按事件类型分派到默认处理器或用户自定义处理器；
 * - 统计阶段：累计各类碰撞次数与平均碰撞率。
 */

#pragma once

#include "CollisonTypes.h"
#include "../../../Particle/include/ParticleManager.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <tuple>
#include <vector>

namespace SCDAT {
namespace Particle {

using CollisionType = SCDAT::Collision::CollisionType;
using PairCrossSectionCalculator = SCDAT::Collision::PairCrossSectionCalculator;
using CollisionParameters = SCDAT::Collision::CollisionParameters;

struct CollisionEvent {
    /// 第一个粒子 ID
    ParticleId particle1_id;
    /// 第二个粒子 ID
    ParticleId particle2_id;
    /// 碰撞类型
    CollisionType type;
    /// 本事件对应的截面
    double cross_section;
    /// 抽样概率
    double probability;
    /// 事件位置（通常取两粒子中点）
    SCDAT::Geometry::Point3D position;
    /// 事件发生时刻（当前实现保留字段）
    double time;

    CollisionEvent(ParticleId p1, ParticleId p2, CollisionType t, double cs, double prob,
                   const SCDAT::Geometry::Point3D& pos, double tm)
        : particle1_id(p1), particle2_id(p2), type(t), cross_section(cs), probability(prob),
          position(pos), time(tm)
    {
    }
};

class CollisionDetector {
  public:
        /**
         * @brief 构造碰撞检测器并绑定粒子管理器。
         * @details 同时注册默认处理器（弹性/非弹性/电荷交换/电离）。
         */
    explicit CollisionDetector(ParticleManager& particle_manager) : particle_manager_(particle_manager)
    {
        resetStatistics();

        collision_handlers_[CollisionType::ELASTIC] =
            [this](const CollisionEvent& event, ParticleManager& manager) {
                handleElasticCollision(event, manager);
            };

        collision_handlers_[CollisionType::INELASTIC] =
            [this](const CollisionEvent& event, ParticleManager& manager) {
                handleInelasticCollision(event, manager);
            };

        collision_handlers_[CollisionType::CHARGE_EXCHANGE] =
            [this](const CollisionEvent& event, ParticleManager& manager) {
                handleChargeExchange(event, manager);
            };

        collision_handlers_[CollisionType::IONIZATION] =
            [this](const CollisionEvent& event, ParticleManager& manager) {
                handleIonization(event, manager);
            };

        collision_handlers_[CollisionType::EXCITATION] =
            [this](const CollisionEvent& event, ParticleManager& manager) {
                handleExcitation(event, manager);
            };

        collision_handlers_[CollisionType::RECOMBINATION] =
            [this](const CollisionEvent& event, ParticleManager& manager) {
                handleRecombination(event, manager);
            };

        initializeDefaultConfiguration();
    }

    ~CollisionDetector() = default;

    /**
     * @brief 检测时间步内的候选碰撞事件。
     * @param dt 时间步长（s）。
     * @return 检测到的碰撞事件列表。
     * @details
     * 使用两两粒子扫描（O(N^2)）与概率抽样。
     * 后续若粒子数较大可替换为空间网格/邻域搜索以降低复杂度。
     */
    std::vector<CollisionEvent> detectCollisions(double dt)
    {
        std::vector<CollisionEvent> events;
        const auto particles = particle_manager_.getAllParticles();

        for (size_t i = 0; i < particles.size(); ++i)
        {
            for (size_t j = i + 1; j < particles.size(); ++j)
            {
                auto& p1 = particles[i];
                auto& p2 = particles[j];
                if (!p1 || !p2)
                {
                    continue;
                }
                if (p1->getStatus() != ParticleStatus::ACTIVE ||
                    p2->getStatus() != ParticleStatus::ACTIVE)
                {
                    continue;
                }

                const SCDAT::Geometry::Vector3D rel_vel = p2->getVelocity() - p1->getVelocity();
                const double rel_speed = rel_vel.magnitude();
                if (rel_speed < 1e-10)
                {
                    continue;
                }

                const double max_collision_distance = 1e-6;
                if (!isInCollisionRange(*p1, *p2, max_collision_distance))
                {
                    continue;
                }

                const std::vector<CollisionType> candidates = {
                    CollisionType::ELASTIC, CollisionType::INELASTIC, CollisionType::CHARGE_EXCHANGE,
                    CollisionType::IONIZATION, CollisionType::EXCITATION,
                    CollisionType::RECOMBINATION};

                for (CollisionType type : candidates)
                {
                    const double prob = calculateCollisionProbability(*p1, *p2, dt, type);
                    if (prob <= 0.0)
                    {
                        continue;
                    }

                    static std::mt19937 rng(std::random_device{}());
                    std::uniform_real_distribution<double> u(0.0, 1.0);
                    if (u(rng) >= prob)
                    {
                        continue;
                    }

                    const SCDAT::Geometry::Point3D pos =
                        (p1->getPosition() + p2->getPosition()) * 0.5;

                    const auto key = std::make_tuple(p1->getType(), p2->getType(), type);
                    double sigma = 0.0;
                    auto it = cross_section_calculators_.find(key);
                    if (it != cross_section_calculators_.end() && it->second)
                    {
                        sigma = it->second->calculateCrossSection(*p1, *p2, rel_speed);
                    }

                    events.emplace_back(p1->getId(), p2->getId(), type, sigma, prob, pos, 0.0);
                    break;
                }
            }
        }

        return events;
    }

    /**
     * @brief 执行碰撞事件并更新统计。
     * @param events 待处理事件列表。
     */
    void processCollisions(const std::vector<CollisionEvent>& events)
    {
        for (const auto& event : events)
        {
            ++statistics_.total_collisions;
            statistics_.total_cross_section += event.cross_section;

            switch (event.type)
            {
            case CollisionType::ELASTIC:
                ++statistics_.elastic_collisions;
                break;
            case CollisionType::INELASTIC:
                ++statistics_.inelastic_collisions;
                break;
            case CollisionType::CHARGE_EXCHANGE:
                ++statistics_.charge_exchange_collisions;
                break;
            case CollisionType::IONIZATION:
                ++statistics_.ionization_collisions;
                break;
            case CollisionType::EXCITATION:
                ++statistics_.excitation_collisions;
                break;
            case CollisionType::RECOMBINATION:
                ++statistics_.recombination_collisions;
                break;
            default:
                break;
            }

            auto it = collision_handlers_.find(event.type);
            if (it != collision_handlers_.end())
            {
                it->second(event, particle_manager_);
            }
        }

        const auto stats = particle_manager_.getStatistics();
        if (statistics_.total_collisions > 0 && stats.total_particles > 0)
        {
            statistics_.average_collision_rate =
                static_cast<double>(statistics_.total_collisions) /
                static_cast<double>(stats.total_particles);
        }
    }

    /**
     * @brief 注册粒子对截面计算器。
     * @details 会自动注册(type2,type1)的反向映射。
     */
    void addCrossSectionCalculator(ParticleType type1, ParticleType type2,
                                   CollisionType collision_type,
                                   std::shared_ptr<PairCrossSectionCalculator> calculator)
    {
        const auto key = std::make_tuple(type1, type2, collision_type);
        const auto reverse_key = std::make_tuple(type2, type1, collision_type);
        cross_section_calculators_[key] = calculator;
        cross_section_calculators_[reverse_key] = std::move(calculator);
    }

    /**
     * @brief 设置指定类型碰撞处理器。
     * @details 可覆盖默认处理逻辑，便于扩展定制化物理过程。
     */
    void setCollisionHandler(CollisionType collision_type,
                             std::function<void(const CollisionEvent&, ParticleManager&)> handler)
    {
        collision_handlers_[collision_type] = std::move(handler);
    }

    /**
     * @brief 获取指定碰撞类型的默认参数。
     * @throws std::out_of_range 当该类型尚未配置默认参数。
     */
    const CollisionParameters& getDefaultCollisionParameters(CollisionType collision_type) const
    {
        return default_collision_parameters_.at(collision_type);
    }

    /**
     * @brief 覆盖指定碰撞类型的默认参数。
     */
    void setDefaultCollisionParameters(CollisionType collision_type,
                                       const CollisionParameters& params)
    {
        default_collision_parameters_[collision_type] = params;
    }

    /**
     * @brief 重新执行默认初始化配置。
     * @details
     * 用于用户希望回到框架默认参数与默认截面注册场景。
     */
    void initializeDefaultConfiguration()
    {
        initializeDefaultCollisionParameters();
        registerDefaultPairCrossSectionCalculators();
    }

    /**
     * @brief 碰撞统计结构。
     */
    struct CollisionStatistics {
        size_t total_collisions = 0;
        size_t elastic_collisions = 0;
        size_t inelastic_collisions = 0;
        size_t charge_exchange_collisions = 0;
        size_t ionization_collisions = 0;
        size_t excitation_collisions = 0;
        size_t recombination_collisions = 0;
        double total_cross_section = 0.0;
        double average_collision_rate = 0.0;
    };

    /**
     * @brief 获取统计快照。
     */
    CollisionStatistics getStatistics() const
    {
        return statistics_;
    }

    /**
     * @brief 重置统计计数。
     */
    void resetStatistics()
    {
        statistics_ = CollisionStatistics{};
    }

  private:
    /**
     * @brief 初始化默认碰撞参数。
     */
    void initializeDefaultCollisionParameters()
    {
        CollisionParameters elastic;
        elastic.energy_threshold = 0.0;
        elastic.cross_section_max = 2.0e-19;
        elastic.energy_loss = 0.0;
        elastic.reaction_equation = "A + B -> A + B";
        elastic.additional_params["reference_energy"] = 10.0;
        elastic.additional_params["power_index"] = 0.3;
        default_collision_parameters_[CollisionType::ELASTIC] = elastic;

        CollisionParameters inelastic;
        inelastic.energy_threshold = 1.0;
        inelastic.cross_section_max = 8.0e-20;
        inelastic.energy_loss = 1.0;
        inelastic.reaction_equation = "A + B -> A* + B";
        default_collision_parameters_[CollisionType::INELASTIC] = inelastic;

        CollisionParameters ionization;
        ionization.energy_threshold = 15.76;
        ionization.cross_section_max = 3.0e-20;
        ionization.energy_loss = 15.76;
        ionization.reaction_equation = "e + A -> 2e + A+";
        default_collision_parameters_[CollisionType::IONIZATION] = ionization;

        CollisionParameters excitation;
        excitation.energy_threshold = 11.5;
        excitation.cross_section_max = 2.5e-20;
        excitation.energy_loss = 11.5;
        excitation.reaction_equation = "e + A -> e + A*";
        excitation.additional_params["alpha"] = 1.2;
        excitation.additional_params["beta"] = 0.35;
        default_collision_parameters_[CollisionType::EXCITATION] = excitation;

        CollisionParameters charge_exchange;
        charge_exchange.energy_threshold = 0.0;
        charge_exchange.cross_section_max = 1.2e-19;
        charge_exchange.energy_loss = 0.2;
        charge_exchange.reaction_equation = "A+ + B -> A + B+";
        default_collision_parameters_[CollisionType::CHARGE_EXCHANGE] = charge_exchange;

        CollisionParameters recombination;
        recombination.energy_threshold = 0.0;
        recombination.cross_section_max = 5.0e-19;
        recombination.energy_loss = 0.0;
        recombination.reaction_equation = "e + A+ -> A";
        recombination.additional_params["reference_energy"] = 1.0;
        recombination.additional_params["power_index"] = 0.7;
        default_collision_parameters_[CollisionType::RECOMBINATION] = recombination;

        CollisionParameters coulomb;
        coulomb.energy_threshold = 0.0;
        coulomb.cross_section_max = 1.0e-18;
        coulomb.reaction_equation = "q1 + q2 -> q1 + q2";
        coulomb.additional_params["electron_density"] = 1.0e16;
        coulomb.additional_params["electron_temperature"] = 1000.0;
        default_collision_parameters_[CollisionType::COULOMB] = coulomb;
    }

    /**
     * @brief 注册默认粒子对截面计算器。
     * @details
     * 使用硬球近似给出常见碰撞通道，确保检测层开箱即用。
     */
    void registerDefaultPairCrossSectionCalculators()
    {
        auto register_hs = [this](ParticleType t1, ParticleType t2, CollisionType ctype,
                                  double r1, double r2) {
            addCrossSectionCalculator(t1, t2, ctype,
                                      std::make_shared<SCDAT::Collision::HardSpherePairCalculator>(
                                          r1, r2));
        };

        const double r_e = 2.0e-10;
        const double r_ion = 1.5e-10;
        const double r_atom = 1.4e-10;

        register_hs(ParticleType::ELECTRON, ParticleType::ATOM, CollisionType::ELASTIC, r_e,
                    r_atom);
        register_hs(ParticleType::ELECTRON, ParticleType::ATOM, CollisionType::INELASTIC, r_e,
                    r_atom);
        register_hs(ParticleType::ELECTRON, ParticleType::ATOM, CollisionType::IONIZATION, r_e,
                    r_atom);
        register_hs(ParticleType::ELECTRON, ParticleType::ATOM, CollisionType::EXCITATION, r_e,
                    r_atom);

        register_hs(ParticleType::ION, ParticleType::ATOM, CollisionType::ELASTIC, r_ion, r_atom);
        register_hs(ParticleType::ION, ParticleType::ATOM, CollisionType::CHARGE_EXCHANGE, r_ion,
                    r_atom);

        register_hs(ParticleType::ELECTRON, ParticleType::ION, CollisionType::ELASTIC, r_e, r_ion);
        register_hs(ParticleType::ELECTRON, ParticleType::ION, CollisionType::RECOMBINATION, r_e,
                    r_ion);

        addCrossSectionCalculator(ParticleType::ELECTRON, ParticleType::ION,
                                  CollisionType::COULOMB,
                                  std::make_shared<SCDAT::Collision::CoulombPairCalculator>(
                                      2.0e-6));
    }

    /**
     * @brief 计算粒子对碰撞概率。
     * @details
     * 当前实现采用：P = sigma * v_rel * n_ref * dt
     * 其中 n_ref 为参考密度。后续可替换为局部网格密度。
     */
    double calculateCollisionProbability(const Particle& particle1, const Particle& particle2,
                                         double dt, CollisionType collision_type) const
    {
        const auto key = std::make_tuple(particle1.getType(), particle2.getType(), collision_type);
        auto it = cross_section_calculators_.find(key);
        if (it == cross_section_calculators_.end() || !it->second)
        {
            return 0.0;
        }

        const SCDAT::Geometry::Vector3D rel_vel = particle2.getVelocity() - particle1.getVelocity();
        const double rel_speed = rel_vel.magnitude();
        if (rel_speed < 1e-10)
        {
            return 0.0;
        }

        const double sigma = it->second->calculateCrossSection(particle1, particle2, rel_speed);
        const double probability = sigma * rel_speed * 1e15 * dt;
        return std::min(1.0, std::max(0.0, probability));
    }

    /**
     * @brief 距离门限判据。
     */
    bool isInCollisionRange(const Particle& particle1, const Particle& particle2,
                            double max_distance) const
    {
        const SCDAT::Geometry::Vector3D d = particle2.getPosition() - particle1.getPosition();
        return d.magnitude() <= max_distance;
    }

    /**
     * @brief 默认弹性碰撞处理。
     * @details 保持两体体系动量守恒，近似处理散射后速度更新。
     */
    void handleElasticCollision(const CollisionEvent& event, ParticleManager& manager)
    {
        auto p1 = manager.getParticle(event.particle1_id);
        auto p2 = manager.getParticle(event.particle2_id);
        if (!p1 || !p2)
        {
            return;
        }

        const double m1 = p1->getMass();
        const double m2 = p2->getMass();

        const SCDAT::Geometry::Vector3D v1 = p1->getVelocity();
        const SCDAT::Geometry::Vector3D v2 = p2->getVelocity();
        const SCDAT::Geometry::Vector3D vcm = (v1 * m1 + v2 * m2) / (m1 + m2);
        const SCDAT::Geometry::Vector3D vrel = v1 - v2;

        p1->setVelocity(vcm + vrel * (m2 / (m1 + m2)));
        p2->setVelocity(vcm - vrel * (m1 / (m1 + m2)));
    }

    /**
     * @brief 默认非弹性碰撞处理。
     * @details 在弹性模型基础上衰减相对速度，模拟动能损失。
     */
    void handleInelasticCollision(const CollisionEvent& event, ParticleManager& manager)
    {
        auto p1 = manager.getParticle(event.particle1_id);
        auto p2 = manager.getParticle(event.particle2_id);
        if (!p1 || !p2)
        {
            return;
        }

        const double m1 = p1->getMass();
        const double m2 = p2->getMass();
        const SCDAT::Geometry::Vector3D v1 = p1->getVelocity();
        const SCDAT::Geometry::Vector3D v2 = p2->getVelocity();
        const SCDAT::Geometry::Vector3D vcm = (v1 * m1 + v2 * m2) / (m1 + m2);
        const SCDAT::Geometry::Vector3D vrel = (v1 - v2) * std::sqrt(0.5);

        p1->setVelocity(vcm + vrel * (m2 / (m1 + m2)));
        p2->setVelocity(vcm - vrel * (m1 / (m1 + m2)));
    }

    /**
     * @brief 默认电荷交换处理。
     */
    void handleChargeExchange(const CollisionEvent& event, ParticleManager& manager)
    {
        auto p1 = manager.getParticle(event.particle1_id);
        auto p2 = manager.getParticle(event.particle2_id);
        if (!p1 || !p2)
        {
            return;
        }

        const double q1 = p1->getCharge();
        const double q2 = p2->getCharge();
        p1->setCharge(q2);
        p2->setCharge(q1);
    }

    /**
     * @brief 默认电离处理。
     * @details 创建次级电子并将中性粒子转为离子，附加反冲速度。
     */
    void handleIonization(const CollisionEvent& event, ParticleManager& manager)
    {
        auto p1 = manager.getParticle(event.particle1_id);
        auto p2 = manager.getParticle(event.particle2_id);
        if (!p1 || !p2)
        {
            return;
        }

        if (std::abs(p2->getCharge()) > 1e-20)
        {
            return;
        }

        static std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> angle(0.0, 2.0 * 3.14159265358979323846);
        std::uniform_real_distribution<double> cos_dist(-1.0, 1.0);

        const double e_ev = 1.0;
        const double speed = std::sqrt(2.0 * e_ev * 1.602176634e-19 / 9.10938356e-31);
        const double theta = std::acos(cos_dist(rng));
        const double phi = angle(rng);

        const SCDAT::Geometry::Vector3D ev(speed * std::sin(theta) * std::cos(phi),
                            speed * std::sin(theta) * std::sin(phi),
                            speed * std::cos(theta));

        manager.createElectron(event.position, ev, 1.0);
        p2->setCharge(1.602176634e-19);

        const SCDAT::Geometry::Vector3D recoil =
            ev * (-9.10938356e-31 / std::max(1e-30, p2->getMass()));
        p2->setVelocity(p2->getVelocity() + recoil);
    }

    /**
     * @brief 默认激发处理。
     * @details
     * 对两体相对动能施加衰减，并保持系统中心质点速度不变。
     */
    void handleExcitation(const CollisionEvent& event, ParticleManager& manager)
    {
        auto p1 = manager.getParticle(event.particle1_id);
        auto p2 = manager.getParticle(event.particle2_id);
        if (!p1 || !p2)
        {
            return;
        }

        const double m1 = p1->getMass();
        const double m2 = p2->getMass();
        const SCDAT::Geometry::Vector3D v1 = p1->getVelocity();
        const SCDAT::Geometry::Vector3D v2 = p2->getVelocity();
        const SCDAT::Geometry::Vector3D vcm = (v1 * m1 + v2 * m2) / std::max(1e-30, m1 + m2);
        const SCDAT::Geometry::Vector3D vrel = v1 - v2;

        const double excitation_loss_factor = std::sqrt(0.8);
        const SCDAT::Geometry::Vector3D vrel_after = vrel * excitation_loss_factor;

        p1->setVelocity(vcm + vrel_after * (m2 / std::max(1e-30, m1 + m2)));
        p2->setVelocity(vcm - vrel_after * (m1 / std::max(1e-30, m1 + m2)));
    }

    /**
     * @brief 默认复合处理。
     * @details
     * 对异号带电粒子执行复合：电荷归零、速度并拢并标记为复合态。
     */
    void handleRecombination(const CollisionEvent& event, ParticleManager& manager)
    {
        auto p1 = manager.getParticle(event.particle1_id);
        auto p2 = manager.getParticle(event.particle2_id);
        if (!p1 || !p2)
        {
            return;
        }

        if (p1->getCharge() * p2->getCharge() >= 0.0)
        {
            return;
        }

        const double m1 = p1->getMass();
        const double m2 = p2->getMass();
        const SCDAT::Geometry::Vector3D drift =
            (p1->getVelocity() * m1 + p2->getVelocity() * m2) / std::max(1e-30, m1 + m2);

        p1->setCharge(0.0);
        p2->setCharge(0.0);
        p1->setVelocity(drift);
        p2->setVelocity(drift);
        p1->setStatus(ParticleStatus::RECOMBINED);
        p2->setStatus(ParticleStatus::RECOMBINED);
    }

    ParticleManager& particle_manager_;

    std::map<std::tuple<ParticleType, ParticleType, CollisionType>,
             std::shared_ptr<PairCrossSectionCalculator>>
        cross_section_calculators_;

    std::map<CollisionType, std::function<void(const CollisionEvent&, ParticleManager&)>>
        collision_handlers_;

    std::map<CollisionType, CollisionParameters> default_collision_parameters_;

    CollisionStatistics statistics_;
};

using CollisionDetectorPtr = std::shared_ptr<CollisionDetector>;

} // namespace Particle
} // namespace SCDAT
