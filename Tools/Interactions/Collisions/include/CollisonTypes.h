/**
 * @file CollisonTypes.h
 * @brief 粒子-粒子相互作用模型与碰撞截面定义
 * @details
 * 本文件集中定义：
 * 1) 碰撞类型与粒子种类枚举；
 * 2) 各类碰撞截面模型（弹性/电离/库仑）；
 * 3) 粒子对截面计算器（硬球/库仑上界）；
 * 4) 多种粒子相互作用模型与统一管理器。
 *
 * 设计目标是将“物理模型与截面计算”与“算法调度（MCC/DSMC）”解耦，
 * 以便后续替换截面公式或加入实验拟合参数时无需改动算法框架。
 */

#pragma once

#include "../../../Geometry/include/Point3D.h"
#include "../../../Geometry/include/Vector3D.h"
#include "../../../Particle/include/ParticleDefinitions.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Collision
{

using ParticleObject = SCDAT::Particle::Particle;
using ParticleType = SCDAT::Particle::ParticleType;
using Vector3D = SCDAT::Geometry::Vector3D;

enum class CollisionType
{
    ELASTIC = 0,
    INELASTIC = 1,
    IONIZATION = 2,
    EXCITATION = 3,
    ATTACHMENT = 4,
    CHARGE_EXCHANGE = 5,
    COULOMB = 6,
    RECOMBINATION = 7,
    DSMC = 8
};

/**
 * @brief 宏观算法中的目标粒子类别。
 * @details
 * 该枚举用于 CollisionProcess 中描述“背景气体/背景等离子体”的种类，
 * 与 ParticleType 的细粒度枚举形成“算法层 vs 物理层”的分层。
 */
enum class ParticleSpecies
{
    ELECTRON = 0,
    ION = 1,
    NEUTRAL = 2,
    PHOTON = 3
};

/**
 * @brief 粒子类型物理属性集合。
 * @param name 可读名称（用于日志、模型命名）。
 * @param mass_amu 质量（amu）。
 * @param charge_in_e 电荷数（单位为 e）。
 * @param type 工程枚举类型（ParticleType）。
 */
struct ParticleTypeInfo
{
    std::string name;
    double mass_amu;
    double charge_in_e;
    ParticleType type;

    ParticleTypeInfo(const std::string& n, double m, double c, ParticleType t)
        : name(n), mass_amu(m), charge_in_e(c), type(t)
    {
    }
};

/**
 * @brief 碰撞模型参数。
 * @details
 * - energy_threshold: 反应阈值能量（eV）
 * - cross_section_max: 截面上限（m^2）
 * - energy_loss: 单次碰撞能损（eV）
 * - mass_ratio: 质量比（模型可选项）
 * - reaction_equation: 反应描述文本
 * - additional_params: 扩展参数字典（用于经验模型拟合）
 */
struct CollisionParameters
{
    double energy_threshold = 0.0;
    double cross_section_max = 1e-20;
    double energy_loss = 0.0;
    double mass_ratio = 1.0;
    std::string reaction_equation;
    std::map<std::string, double> additional_params;
};

/**
 * @brief 碰撞截面基类。
 * @details
 * 统一定义 “sigma(E)” 接口，输入粒子能量（eV），输出截面（m^2）。
 * 派生类可采用理论公式、经验公式或查表插值实现。
 */
class CollisionCrossSection
{
  public:
    CollisionCrossSection(CollisionType type, const CollisionParameters& params);
    virtual ~CollisionCrossSection() = default;

    virtual double calculateCrossSection(double energy_ev) const = 0;

    bool isEnergyAboveThreshold(double energy_ev) const;
    CollisionType getCollisionType() const;
    const CollisionParameters& getParameters() const;

  protected:
    CollisionType collision_type_;
    CollisionParameters params_;
};

/**
 * @brief 弹性碰撞截面模型。
 * @details 采用简单的常数截面作为示例，实际应用中可替换为更复杂的能量依赖模型或查表数据。
 */
class ElasticCollisionCrossSection : public CollisionCrossSection
{
  public:
    explicit ElasticCollisionCrossSection(const CollisionParameters& params);
    double calculateCrossSection(double energy_ev) const override;
};

/**
 * @brief 电离截面模型（阈值型）。
 * @details
 * 采用近似 Bethe-Born 形式，满足 E <= I 时截面为 0。
 */
class IonizationCollisionCrossSection : public CollisionCrossSection
{
  public:
    explicit IonizationCollisionCrossSection(const CollisionParameters& params);
    double calculateCrossSection(double energy_ev) const override;
};

/**
 * @brief 激发截面模型（阈值型）。
 * @details
 * 当 E <= Eth 时截面为 0；超过阈值后先升高再随能量缓慢衰减，
 * 适用于缺少查表数据时的经验近似。
 */
class ExcitationCollisionCrossSection : public CollisionCrossSection
{
  public:
    explicit ExcitationCollisionCrossSection(const CollisionParameters& params);
    double calculateCrossSection(double energy_ev) const override;
};

/**
 * @brief 复合截面模型。
 * @details
 * 对电子-离子复合采用低能更强、高能更弱的经验形式。
 */
class RecombinationCollisionCrossSection : public CollisionCrossSection
{
  public:
    explicit RecombinationCollisionCrossSection(const CollisionParameters& params);
    double calculateCrossSection(double energy_ev) const override;
};
/**
 * @brief 库仑碰撞截面模型。
 * @details
 * 采用经典库仑散射理论，截面随能量降低而增加，
 * 但在数值实现中引入 Debye 长度和 Coulomb 对数以避免发散。
 */
class CoulombCollisionCrossSection : public CollisionCrossSection
{
  public:
    explicit CoulombCollisionCrossSection(const CollisionParameters& params);
    double calculateCrossSection(double energy_ev) const override;

    static double calculateDebyeLength(double temperature_k, double density_m3);
    static double calculateCoulombLogarithm(double temperature_k, double density_m3);

  private:
    double electron_density_;
    double electron_temperature_;
    double coulomb_logarithm_;
};

/**
 * @brief 粒子对截面计算器接口。
 * @details
 * 与能量截面不同，该接口显式接收两粒子状态与相对速度，
 * 便于在检测器中使用“几何判据 + 物理判据”混合计算。
 */
class PairCrossSectionCalculator
{
  public:
    virtual ~PairCrossSectionCalculator() = default;

    virtual double calculateCrossSection(const ParticleObject& p1, const ParticleObject& p2,
                                         double relative_velocity) const = 0;
};

class HardSpherePairCalculator : public PairCrossSectionCalculator
{
  public:
    HardSpherePairCalculator(double radius1, double radius2);

    double calculateCrossSection(const ParticleObject& p1, const ParticleObject& p2,
                                 double relative_velocity) const override;

  private:
    double radius1_;
    double radius2_;
};

/**
 * @brief 库仑粒子对截面计算器。
 * @details
 * 对长程库仑作用引入 impact_parameter_limit 作为数值上界，
 * 防止极低能情况下截面数值发散。
 */
class CoulombPairCalculator : public PairCrossSectionCalculator
{
  public:
    explicit CoulombPairCalculator(double impact_parameter_limit);

    double calculateCrossSection(const ParticleObject& p1, const ParticleObject& p2,
                                 double relative_velocity) const override;

  private:
    double impact_parameter_limit_;
};

class CollisionModel
{
  public:
    explicit CollisionModel(const std::string& name);
    virtual ~CollisionModel();

    void setEnabled(bool enabled);
    bool isEnabled() const;
    const std::string& getName() const;
    size_t getTotalCollisions() const;
    void resetStatistics();

    virtual size_t processCollisions(std::vector<ParticleObject>& particles,
                                     std::vector<ParticleObject>& new_particles, double dt,
                                     double volume) = 0;

  protected:
    std::string name_;
    bool enabled_;
    size_t total_collisions_;
};

/**
 * @brief 弹性碰撞模型。
 * @details
 * 对目标粒子集合进行候选配对并执行散射，保持两体系统动量守恒。
 */
class ElasticCollisionModel : public CollisionModel
{
  public:
    ElasticCollisionModel(const ParticleTypeInfo& type1, const ParticleTypeInfo& type2,
                          double cross_section);

    size_t processCollisions(std::vector<ParticleObject>& particles,
                             std::vector<ParticleObject>& new_particles, double dt,
                             double volume) override;

  private:
    void performElasticCollision(ParticleObject& p1, ParticleObject& p2);

    ParticleTypeInfo type1_;
    ParticleTypeInfo type2_;
    double cross_section_;
};

class IonizationCollisionModel : public CollisionModel
{
  public:
    IonizationCollisionModel(const ParticleTypeInfo& electron_type,
                             const ParticleTypeInfo& neutral_type, const ParticleTypeInfo& ion_type,
                             double cross_section, double ionization_energy_ev);

    size_t processCollisions(std::vector<ParticleObject>& particles,
                             std::vector<ParticleObject>& new_particles, double dt,
                             double volume) override;

  private:
    bool performIonization(ParticleObject& electron, const ParticleObject& neutral,
                           std::vector<ParticleObject>& new_particles);

    ParticleTypeInfo electron_type_;
    ParticleTypeInfo neutral_type_;
    ParticleTypeInfo ion_type_;
    double cross_section_;
    double ionization_energy_ev_;
};

/**
 * @brief 电荷交换模型。
 * @details
 * 用于离子-中性粒子反应，典型处理为交换荷态并近似保持速度量级。
 */
class ChargeExchangeCollisionModel : public CollisionModel
{
  public:
    ChargeExchangeCollisionModel(const ParticleTypeInfo& ion_type,
                                 const ParticleTypeInfo& neutral_type, double cross_section);

    size_t processCollisions(std::vector<ParticleObject>& particles,
                             std::vector<ParticleObject>& new_particles, double dt,
                             double volume) override;

  private:
    void performChargeExchange(ParticleObject& ion, ParticleObject& neutral);

    ParticleTypeInfo ion_type_;
    ParticleTypeInfo neutral_type_;
    double cross_section_;
};

class RecombinationCollisionModel : public CollisionModel
{
  public:
    RecombinationCollisionModel(const ParticleTypeInfo& electron_type,
                                const ParticleTypeInfo& ion_type,
                                const ParticleTypeInfo& neutral_type, double cross_section);

    size_t processCollisions(std::vector<ParticleObject>& particles,
                             std::vector<ParticleObject>& new_particles, double dt,
                             double volume) override;

  private:
    ParticleTypeInfo electron_type_;
    ParticleTypeInfo ion_type_;
    ParticleTypeInfo neutral_type_;
    double cross_section_;
};

/**
 * @brief 碰撞模型管理器。
 * @details
 * 负责模型生命周期管理与统计归并。
 */
class CollisionModelManager
{
  public:
    CollisionModelManager();
    ~CollisionModelManager();

    void addModel(std::shared_ptr<CollisionModel> model);
    void removeModel(const std::string& name);
    std::shared_ptr<CollisionModel> getModel(const std::string& name);

    const std::vector<std::shared_ptr<CollisionModel>>& getModels() const;

    void setModelEnabled(const std::string& name, bool enabled);
    std::vector<std::string> getModelNames() const;
    size_t getTotalCollisions() const;
    void resetStatistics();
    void clearModels();

  private:
    std::vector<std::shared_ptr<CollisionModel>> models_;
};

} // namespace Collision
} // namespace SCDAT
