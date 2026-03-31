/**
 * @file CollisionAlgorithm.h
 * @brief 粒子碰撞算法（Monte Carlo 与 DSMC）
 * @details
 * 本文件定义算法层接口与流程控制：
 * - CollisionProcess: 单一碰撞过程的抽象执行单元；
 * - MonteCarloCollisionHandler: 经典 MCC 主循环；
 * - DSMCCollisionHandler: 单元内随机配对的 DSMC 框架入口；
 * - CollisionProcessFactory: 依据参数快速生成过程对象。
 */

#pragma once

#include "CollisonTypes.h"

#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Collision
{

using CollisionCrossSectionPtr = std::shared_ptr<CollisionCrossSection>;

/**
 * @brief 统一调度已启用的碰撞模型并执行碰撞。
 * @details
 * 将“遍历模型 + enabled 判定 + 累加碰撞计数”的算法逻辑集中在算法层，
 * 便于管理器与不同算法入口复用同一实现。
 */
size_t processEnabledCollisionModels(std::vector<std::shared_ptr<CollisionModel>>& models,
                                     std::vector<ParticleObject>& particles,
                                     std::vector<ParticleObject>& new_particles, double dt,
                                     double volume);

/**
 * @brief 通过模型管理器统一执行已启用模型的碰撞计算。
 * @details
 * 该入口将碰撞调度主循环保留在算法层，实现“Types 负责管理、Algorithm 负责计算”的分层边界。
 */
size_t processAllCollisions(CollisionModelManager& manager,
                            std::vector<ParticleObject>& particles,
                            std::vector<ParticleObject>& new_particles, double dt,
                            double volume);

/**
 * @brief 单一碰撞过程抽象基类。
 * @details
 * 将“碰撞概率计算”与“碰撞后状态更新”封装到同一对象中，
 * 便于 MonteCarloCollisionHandler 进行统一调度。
 */
class CollisionProcess
{
  public:
    CollisionProcess(CollisionCrossSectionPtr cross_section, ParticleSpecies incident_species,
                     ParticleSpecies target_species, double target_density);
    virtual ~CollisionProcess() = default;

    double calculateCollisionFrequency(const ParticleObject& particle) const;
    bool isApplicableTo(const ParticleObject& particle) const;

    virtual std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                         std::mt19937& rng) = 0;

    ParticleSpecies getIncidentSpecies() const;
    CollisionCrossSectionPtr getCrossSection() const;
    ParticleSpecies getTargetSpecies() const;

    void setTargetDensity(double density);
    double getTargetDensity() const;

  protected:
    double calculateKineticEnergyEv(const ParticleObject& particle) const;
    double generateIsotropicScatteringAngle(std::mt19937& rng) const;
    double generateRandomAzimuthalAngle(std::mt19937& rng) const;
    static ParticleSpecies classifyParticleSpecies(const ParticleObject& particle);

    CollisionCrossSectionPtr cross_section_;
    ParticleSpecies incident_species_;
    ParticleSpecies target_species_;
    double target_density_;
};

class ElasticCollisionProcess : public CollisionProcess
{
  public:
    ElasticCollisionProcess(CollisionCrossSectionPtr cross_section, ParticleSpecies incident_species,
                            ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

class InelasticCollisionProcess : public CollisionProcess
{
  public:
    InelasticCollisionProcess(CollisionCrossSectionPtr cross_section,
                              ParticleSpecies incident_species,
                              ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

/**
 * @brief 电离碰撞过程。
 * @details
 * 在判定碰撞发生后，可生成二次电子和离子等新粒子。
 */
class IonizationCollisionProcess : public CollisionProcess
{
  public:
    IonizationCollisionProcess(CollisionCrossSectionPtr cross_section,
                               ParticleSpecies incident_species,
                               ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

class ExcitationCollisionProcess : public CollisionProcess
{
  public:
    ExcitationCollisionProcess(CollisionCrossSectionPtr cross_section,
                               ParticleSpecies incident_species,
                               ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

class RecombinationCollisionProcess : public CollisionProcess
{
  public:
    RecombinationCollisionProcess(CollisionCrossSectionPtr cross_section,
                                  ParticleSpecies incident_species,
                                  ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

class ChargeExchangeCollisionProcess : public CollisionProcess
{
  public:
    ChargeExchangeCollisionProcess(CollisionCrossSectionPtr cross_section,
                                   ParticleSpecies incident_species,
                                   ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

class CoulombCollisionProcess : public CollisionProcess
{
  public:
    CoulombCollisionProcess(CollisionCrossSectionPtr cross_section, ParticleSpecies incident_species,
                            ParticleSpecies target_species, double target_density);

    std::vector<ParticleObject> executeCollision(ParticleObject& particle, double dt,
                                                 std::mt19937& rng) override;
};

/**
 * @brief Monte Carlo 碰撞处理器。
 * @details
 * 对粒子列表逐个执行：
 * 1) 计算每个过程的碰撞概率；
 * 2) 随机抽样决定是否发生；
 * 3) 执行碰撞并累积统计信息。
 */
class MonteCarloCollisionHandler
{
  public:
    explicit MonteCarloCollisionHandler(unsigned int random_seed = 12345);

    void addCollisionProcess(const std::string& process_name,
                             std::shared_ptr<CollisionProcess> process);
    void removeCollisionProcess(const std::string& process_name);
    void clearCollisionProcesses();

    std::vector<ParticleObject> processCollisions(std::vector<ParticleObject>& particles,
                                                  double dt);

    /**
     * @brief Monte Carlo 统计结构。
     */
    struct CollisionStatistics
    {
        std::map<std::string, int> collision_counts;
        std::map<std::string, double> total_rates;
        int total_collisions = 0;
        int particles_created = 0;
    };

    const CollisionStatistics& getStatistics() const;
    void resetStatistics();

    void setCollisionEnabled(bool enabled);
    bool isCollisionEnabled() const;

    const CollisionParameters& getDefaultCollisionParameters(CollisionType collision_type) const;
    void setDefaultCollisionParameters(CollisionType collision_type,
                                       const CollisionParameters& params);

    void initializeDefaultConfiguration(double neutral_density = 1.0e19,
                                        double ion_density = 1.0e16,
                                        double electron_density = 1.0e16);

  private:
    bool shouldCollisionOccur(double collision_probability);
    void updateStatistics(const std::string& process_name, int new_particles_count);
    void initializeDefaultCollisionParameters();
    void registerDefaultCollisionProcesses(double neutral_density, double ion_density,
                                           double electron_density);

    std::map<std::string, std::shared_ptr<CollisionProcess>> collision_processes_;
    std::map<CollisionType, CollisionParameters> default_collision_parameters_;
    std::mt19937 rng_;
    CollisionStatistics statistics_;
    bool collision_enabled_;
};

class DSMCCollisionHandler
{
  public:
    explicit DSMCCollisionHandler(unsigned int random_seed = 24680);

    void setModel(std::shared_ptr<CollisionModel> model);

    size_t processCollisions(std::vector<ParticleObject>& particles,
                             std::vector<ParticleObject>& new_particles, double dt,
                             double cell_volume);

  private:
    std::shared_ptr<CollisionModel> model_;
    std::mt19937 rng_;
};

/**
 * @brief 碰撞过程工厂。
 * @details
 * 将“参数组装 -> 具体过程对象”标准化，减少上层业务重复代码。
 */
class CollisionProcessFactory
{
  public:
    static std::shared_ptr<CollisionProcess>
    createElasticProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                         double target_density,
                         const CollisionParameters& params);

    static std::shared_ptr<CollisionProcess>
    createIonizationProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                            double target_density,
                            const CollisionParameters& params);

    static std::shared_ptr<CollisionProcess>
    createInelasticProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                           double target_density,
                           const CollisionParameters& params);

    static std::shared_ptr<CollisionProcess>
    createCoulombProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                         double target_density,
                         const CollisionParameters& params);

    static std::shared_ptr<CollisionProcess>
    createChargeExchangeProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                                double target_density,
                                const CollisionParameters& params);

    static std::shared_ptr<CollisionProcess>
    createExcitationProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                            double target_density,
                            const CollisionParameters& params);

    static std::shared_ptr<CollisionProcess>
    createRecombinationProcess(ParticleSpecies incident_species, ParticleSpecies target_species,
                               double target_density,
                               const CollisionParameters& params);
};

} // namespace Collision
} // namespace SCDAT
