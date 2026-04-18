/**
 * @file VolInteract.h
 * @brief 体积交互处理模块头文件
 * @details 定义粒子-粒子碰撞、粒子-中性气体相互作用、电离和复合过程
 *
 * @author Wang Sizhan
 * @date 2026.3.26 16:39:53
 * @version V0.0.1
 * @ingroup FieldSolverModule
 */

#ifndef SCDAT_FIELD_VOL_INTERACT_H
#define SCDAT_FIELD_VOL_INTERACT_H

#include "SurfaceFieldVolumeContracts.h"

#include "../Particle/include/ParticleDefinitions.h"
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Field
{

// 类型别名，简化使用
using ParticleType = SCDAT::Particle::ParticleType;
using ParticlePtr = SCDAT::Particle::ParticlePtr;

class VolInteractModel
{
  public:
    explicit VolInteractModel(std::string name);
    virtual ~VolInteractModel();

    const std::string& getName() const;

  private:
    std::string name_;
};

/**
 * @brief 体积交互器基类
 * @details 定义所有体积交互的通用接口
 */
class VolInteractor
{
  public:
    /**
     * @brief 构造函数
     * @param name 交互器名称
     * @param family SPIS 对标的体域交互 family
     */
    explicit VolInteractor(const std::string& name,
                           FieldSolver::NativeVolumeInteractionFamily family);

    /**
     * @brief 虚析构函数
     */
    virtual ~VolInteractor();

    /**
     * @brief 设置是否启用
     * @param enabled 是否启用
     */
    void setEnabled(bool enabled);

    /**
     * @brief 检查是否启用
     * @return 是否启用
     */
    bool isEnabled() const;

    /**
     * @brief 获取交互器名称
     * @return 名称
     */
    const std::string& getName() const;

    /**
     * @brief 获取交互器所属 family
     * @return family
     */
    FieldSolver::NativeVolumeInteractionFamily getFamily() const;

    /**
     * @brief 获取交互器绑定的交互模型
     * @return 若该 family 使用模型对象，则返回其指针；否则返回空
     */
    const VolInteractModel* getModel() const;

    /**
     * @brief 对当前粒子集合执行交互
     * @param particles 全量活跃粒子列表
     * @param dt 时间步长
     * @return 本次真正执行的 family 数，通常为 0 或 1
     */
    virtual std::size_t execute(std::vector<ParticlePtr>& particles, double dt) = 0;

  protected:
    void setModel(std::shared_ptr<VolInteractModel> model);

    std::string name_; ///< 交互器名称
    bool enabled_;     ///< 是否启用
    FieldSolver::NativeVolumeInteractionFamily family_;
    std::shared_ptr<VolInteractModel> model_;
};

/**
 * @brief 弹性碰撞交互器
 * @details 处理两种粒子之间的弹性碰撞
 */
class ElasticCollisionInteractor : public VolInteractor
{
  public:
    /**
     * @brief 构造函数
     * @param name 交互器名称
     * @param family SPIS 对标的体域交互 family
     * @param type1 第一种粒子类型
     * @param type2 第二种粒子类型
     * @param crossSection 碰撞截面 [m²]
     */
    ElasticCollisionInteractor(const std::string& name,
                               FieldSolver::NativeVolumeInteractionFamily family,
                               const ParticleType& type1,
                               const ParticleType& type2, double crossSection);

    std::size_t execute(std::vector<ParticlePtr>& particles, double dt) override;

    /**
     * @brief 计算碰撞交互
     * @param particles1 第一种粒子列表
     * @param particles2 第二种粒子列表
     * @param dt 时间步长
     * @return 是否完成一次显式 family 级执行
     */
    std::size_t computeInteraction(std::vector<ParticlePtr>& particles1,
                                   std::vector<ParticlePtr>& particles2, double dt);

  private:
    /**
     * @brief 执行弹性碰撞
     * @param p1 第一个粒子
     * @param p2 第二个粒子
     */
    void performElasticCollision(ParticlePtr& p1, ParticlePtr& p2);

    ParticleType type1_;  ///< 第一种粒子类型
    ParticleType type2_;  ///< 第二种粒子类型
    double crossSection_; ///< 碰撞截面
    double reducedMass_;  ///< 约化质量
};

class ElasticCollisionsInteractor : public ElasticCollisionInteractor
{
  public:
    ElasticCollisionsInteractor(const std::string& name,
                                FieldSolver::NativeVolumeInteractionFamily family,
                                const ParticleType& type1,
                                const ParticleType& type2, double crossSection);
};

class MCCInteractor final : public ElasticCollisionsInteractor
{
  public:
    MCCInteractor();
};

class CEXInteractor final : public ElasticCollisionsInteractor
{
  public:
    CEXInteractor();
};

class ConstantIonizationInteractor : public VolInteractor
{
  public:
    ConstantIonizationInteractor(const std::string& name,
                                 FieldSolver::NativeVolumeInteractionFamily family,
                                 double ionization_yield_scale = 1.0);

    std::size_t execute(std::vector<ParticlePtr>& particles, double dt) override;

    double ionizationYieldScale() const;

  protected:
    double ionization_yield_scale_;
};

class PhotoIonizationInteractor final : public ConstantIonizationInteractor
{
  public:
    explicit PhotoIonizationInteractor(double ionization_yield_scale = 1.0);
};

class VolInteractionAlongTrajectoryInteractor : public VolInteractor
{
  public:
    VolInteractionAlongTrajectoryInteractor(const std::string& name,
                                            FieldSolver::NativeVolumeInteractionFamily family);

    double integrationTimeSeconds() const;
    void resetIntegrationTime();

  protected:
    void recordIntegrationTime(double dt);

  private:
    double integration_time_s_ = 0.0;
};

class TrajectoryInteractionFromFieldInteractor final : public VolInteractionAlongTrajectoryInteractor
{
  public:
    explicit TrajectoryInteractionFromFieldInteractor(double field_coupling_scale = 1.0);

    std::size_t execute(std::vector<ParticlePtr>& particles, double dt) override;

  private:
    double field_coupling_scale_;
};

class SpinningSpacecraftTrajectoryInteractor final : public VolInteractionAlongTrajectoryInteractor
{
  public:
    explicit SpinningSpacecraftTrajectoryInteractor(double angular_velocity_rad_per_s = 1.0);

    std::size_t execute(std::vector<ParticlePtr>& particles, double dt) override;

  private:
    double angular_velocity_rad_per_s_;
};

/**
 * @brief 体积交互管理器
 * @details 管理所有体积交互过程
 */
class VolInteractionManager
{
  public:
    /**
     * @brief 构造函数
     */
    VolInteractionManager();

    /**
     * @brief 析构函数
     */
    ~VolInteractionManager();

    /**
     * @brief 添加交互器
     * @param interactor 交互器
     */
    void addInteractor(std::shared_ptr<VolInteractor> interactor);

    /**
     * @brief 移除交互器
     * @param name 交互器名称
     */
    void removeInteractor(const std::string& name);

    /**
     * @brief 获取交互器
     * @param name 交互器名称
     * @return 交互器指针
     */
    std::shared_ptr<VolInteractor> getInteractor(const std::string& name);

    /**
     * @brief 执行所有交互
     * @param particles 粒子列表
     * @param dt 时间步长
     * @return 实际执行到的显式 family 数
     */
    std::size_t executeAllInteractions(std::vector<ParticlePtr>& particles, double dt);

    /**
     * @brief 设置交互器启用状态
     * @param name 交互器名称
     * @param enabled 是否启用
     */
    void setInteractorEnabled(const std::string& name, bool enabled);

    /**
     * @brief 获取交互器数量
     * @return 交互器数量
     */
    size_t getInteractorCount() const;

    /**
     * @brief 获取当前启用的交互 family 列表
     * @return 启用 family 列表
     */
    std::vector<FieldSolver::NativeVolumeInteractionFamily> activeFamilies() const;

    /**
     * @brief 清空所有交互器
     */
    void clear();

  private:
    std::vector<std::shared_ptr<VolInteractor>> interactors_; ///< 交互器列表
};

std::shared_ptr<VolInteractor> makeNativeVolumeInteractor(
    FieldSolver::NativeVolumeInteractionFamily family);

VolInteractionManager buildNativeVolumeInteractionManager(
    const std::vector<FieldSolver::NativeVolumeInteractionFamily>& families);

} // namespace Field
} // namespace SCDAT

#endif // SCDAT_FIELD_VOL_INTERACT_H
