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

#include "../Particle/include/ParticleDefinitions.h"
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
     */
    explicit VolInteractor(const std::string& name);

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

  protected:
    std::string name_; ///< 交互器名称
    bool enabled_;     ///< 是否启用
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
     * @param type1 第一种粒子类型
     * @param type2 第二种粒子类型
     * @param crossSection 碰撞截面 [m²]
     */
    ElasticCollisionInteractor(const ParticleType& type1, const ParticleType& type2,
                               double crossSection);

    /**
     * @brief 计算碰撞交互
     * @param particles1 第一种粒子列表
     * @param particles2 第二种粒子列表
     * @param dt 时间步长
     */
    void computeInteraction(std::vector<ParticlePtr>& particles1,
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
     */
    void executeAllInteractions(std::vector<ParticlePtr>& particles, double dt);

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
     * @brief 清空所有交互器
     */
    void clear();

  private:
    std::vector<std::shared_ptr<VolInteractor>> interactors_; ///< 交互器列表
};

} // namespace Field
} // namespace SCDAT

#endif // SCDAT_FIELD_VOL_INTERACT_H
