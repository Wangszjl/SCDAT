#ifndef SCDAT_ADVANCED_BOUNDARY_CONDITIONS_H
#define SCDAT_ADVANCED_BOUNDARY_CONDITIONS_H

/**
 * @file AdvancedBoundaryConditions.h
 * @ingroup ParticleModule
 * @note 该头文件已迁移为 SCDAT 命名空间体系。
 */

#include "../../Geometry/include/Vector3D.h"
#include "../../Geometry/include/Matrix3x3.h"
#include "../../Geometry/include/Point3D.h"
#include "../../Particle/include/ParticleDefinitions.h"
#include "BoundaryType.h"
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{

namespace Particle
{

// 前向声明
class AdvancedBoundaryHandler;
class PeriodicBoundaryCondition;
class ReflectiveBoundaryCondition;
class AbsorbingBoundaryCondition;
class EmissionBoundaryCondition;

// 类型别名
using AdvancedBoundaryHandlerPtr = std::shared_ptr<AdvancedBoundaryHandler>;
using BoundaryConditionPtr = std::shared_ptr<class BoundaryConditionBase>;
using ParticleTypeDef = Particle;

/**
 * @brief 边界统计信息
 */
struct BoundaryStatistics
{
    size_t particles_reflected = 0;     ///< 反射粒子数
    size_t particles_absorbed = 0;      ///< 吸收粒子数
    size_t particles_emitted = 0;       ///< 发射粒子数
    size_t particles_transmitted = 0;   ///< 透射粒子数
    double total_energy_flux = 0.0;     ///< 总能量通量 (J/s)
    double total_charge_flux = 0.0;     ///< 总电荷通量 (A)
    double average_angle = 0.0;         ///< 平均入射角 (rad)
    double surface_temperature = 300.0; ///< 表面温度 (K)

    /**
     * @brief 重置统计信息
     */
    void reset()
    {
        particles_reflected = 0;
        particles_absorbed = 0;
        particles_emitted = 0;
        particles_transmitted = 0;
        total_energy_flux = 0.0;
        total_charge_flux = 0.0;
        average_angle = 0.0;
        surface_temperature = 300.0;
    }

    /**
     * @brief 获取总粒子数
     */
    size_t getTotalParticles() const
    {
        return particles_reflected + particles_absorbed + particles_emitted + particles_transmitted;
    }

    /**
     * @brief 获取反射系数
     */
    double getReflectionCoefficient() const
    {
        size_t total = getTotalParticles();
        return total > 0 ? static_cast<double>(particles_reflected) / static_cast<double>(total)
                         : 0.0;
    }
};

/**
 * @brief 边界条件基类
 */
class BoundaryConditionBase
{
  public:
    /**
     * @brief 构造函数
     * @param type 边界类型
     * @param geometry 几何类型
     */
    BoundaryConditionBase(AdvancedBoundaryType type, BoundaryGeometry geometry)
        : type_(type), geometry_(geometry), enabled_(true)
    {
    }

    /**
     * @brief 虚析构函数
     */
    virtual ~BoundaryConditionBase() = default;

    /**
     * @brief 处理粒子边界相互作用
     * @param particle 入射粒子
     * @param intersection_point 边界交点
     * @param normal 边界法向量
     * @param dt 时间步长
     * @return 发射粒子列表
     */
    virtual std::vector<ParticleTypeDef>
    processParticle(ParticleTypeDef& particle,
                    const SCDAT::Geometry::Point3D& intersection_point,
                    const SCDAT::Geometry::Vector3D& normal, double dt) = 0;

    /**
     * @brief 更新边界统计信息
     * @param particle 粒子
     * @param action 动作类型
     */
    virtual void updateStatistics(const ParticleTypeDef& particle, const std::string& action);

    /**
     * @brief 获取边界类型
     */
    AdvancedBoundaryType getType() const
    {
        return type_;
    }

    /**
     * @brief 获取几何类型
     */
    BoundaryGeometry getGeometry() const
    {
        return geometry_;
    }

    /**
     * @brief 设置启用状态
     */
    void setEnabled(bool enabled)
    {
        enabled_ = enabled;
    }

    /**
     * @brief 检查是否启用
     */
    bool isEnabled() const
    {
        return enabled_;
    }

    /**
     * @brief 获取统计信息
     */
    const BoundaryStatistics& getStatistics() const
    {
        return statistics_;
    }

    /**
     * @brief 重置统计信息
     */
    void resetStatistics()
    {
        statistics_.reset();
    }

  protected:
    AdvancedBoundaryType type_;     ///< 边界类型
    BoundaryGeometry geometry_;     ///< 几何类型
    bool enabled_;                  ///< 是否启用
    BoundaryStatistics statistics_; ///< 统计信息
};

/**
 * @brief 3D周期边界条件
 *
 * 支持非正交周期边界和复杂几何
 */
class PeriodicBoundaryCondition : public BoundaryConditionBase
{
  public:
    /**
     * @brief 构造函数
     * @param period_vectors 周期向量 (最多3个)
     */
    explicit PeriodicBoundaryCondition(const std::vector<SCDAT::Geometry::Vector3D>& period_vectors);

    /**
     * @brief 析构函数
     */
    ~PeriodicBoundaryCondition() override = default;

    /**
     * @brief 处理粒子边界相互作用
     */
    std::vector<ParticleTypeDef>
    processParticle(ParticleTypeDef& particle, const SCDAT::Geometry::Point3D& intersection_point,
                    const SCDAT::Geometry::Vector3D& normal, double dt) override;

    /**
     * @brief 设置周期向量
     * @param period_vectors 周期向量列表
     */
    void setPeriodVectors(const std::vector<SCDAT::Geometry::Vector3D>& period_vectors);

    /**
     * @brief 获取周期向量
     */
    const std::vector<SCDAT::Geometry::Vector3D>& getPeriodVectors() const
    {
        return period_vectors_;
    }

    /**
     * @brief 计算周期映射位置
     * @param position 原始位置
     * @return 映射后位置
     */
    SCDAT::Geometry::Point3D mapPeriodicPosition(const SCDAT::Geometry::Point3D& position) const;

    /**
     * @brief 检查是否为有效周期配置
     */
    bool isValidPeriodicConfiguration() const;

  private:
    std::vector<SCDAT::Geometry::Vector3D> period_vectors_;            ///< 周期向量
    std::vector<std::vector<double>> period_matrix_;         ///< 周期矩阵
    std::vector<std::vector<double>> inverse_period_matrix_; ///< 逆周期矩阵
    std::vector<SCDAT::Geometry::Vector3D> reciprocal_vectors_;        ///< 倒格子向量
    double period_determinant_;                              ///< 周期矩阵行列式
    bool is_orthogonal_;                                     ///< 是否正交周期

    /**
     * @brief 更新周期矩阵
     */
    void updatePeriodMatrix();

    /**
     * @brief 计算最近周期像
     * @param position 位置
     * @return 最近像位置
     */
    SCDAT::Geometry::Point3D
    findNearestPeriodicImage(const SCDAT::Geometry::Point3D& position) const;

    /**
     * @brief 构建非正交周期矩阵
     */
    void buildNonOrthogonalPeriodMatrix();

    /**
     * @brief 计算倒格子
     */
    void calculateReciprocalLattice();

    /**
     * @brief 设置晶格变换
     */
    void setupLatticeTransformations();

    /**
     * @brief 计算矩阵行列式
     */
    double calculateMatrixDeterminant(const std::vector<std::vector<double>>& matrix) const;

    /**
     * @brief 计算矩阵逆
     */
    std::vector<std::vector<double>>
    calculateMatrixInverse(const std::vector<std::vector<double>>& matrix) const;

    /**
     * @brief 矩阵乘法
     */
    std::vector<std::vector<double>>
    multiplyMatrices(const std::vector<std::vector<double>>& A,
                     const std::vector<std::vector<double>>& B) const;
};

/**
 * @brief 曲面反射边界条件
 *
 * 支持任意曲面的精确法向量计算和反射
 */
class ReflectiveBoundaryCondition : public BoundaryConditionBase
{
  public:
    /**
     * @brief 构造函数
     * @param reflection_coefficient 反射系数 (0-1)
     * @param energy_loss_factor 能量损失因子 (0-1)
     */
    ReflectiveBoundaryCondition(double reflection_coefficient = 1.0,
                                double energy_loss_factor = 0.0);

    /**
     * @brief 析构函数
     */
    ~ReflectiveBoundaryCondition() override = default;

    /**
     * @brief 处理粒子边界相互作用
     */
    std::vector<ParticleTypeDef>
    processParticle(ParticleTypeDef& particle,
                    const SCDAT::Geometry::Point3D& intersection_point,
                    const SCDAT::Geometry::Vector3D& normal, double dt) override;

    /**
     * @brief 设置表面函数
     * @param surface_function 表面函数 f(x,y,z) = 0
     * @param gradient_function 梯度函数 ∇f
     */
    void
    setSurfaceFunction(std::function<double(const SCDAT::Geometry::Point3D&)> surface_function,
                       std::function<SCDAT::Geometry::Vector3D(const SCDAT::Geometry::Point3D&)> gradient_function);

    /**
     * @brief 设置反射参数
     */
    void setReflectionCoefficient(double coefficient)
    {
        reflection_coefficient_ = coefficient;
    }
    void setEnergyLossFactor(double factor)
    {
        energy_loss_factor_ = factor;
    }

    /**
     * @brief 计算精确法向量
     * @param point 表面点
     * @return 单位法向量
     */
    SCDAT::Geometry::Vector3D calculateSurfaceNormal(const SCDAT::Geometry::Point3D& point) const;

  private:
    double reflection_coefficient_; ///< 反射系数
    double energy_loss_factor_;     ///< 能量损失因子

    std::function<double(const SCDAT::Geometry::Point3D&)> surface_function_;           ///< 表面函数
    std::function<SCDAT::Geometry::Vector3D(const SCDAT::Geometry::Point3D&)> gradient_function_; ///< 梯度函数

    /**
     * @brief 计算反射速度
     * @param incident_velocity 入射速度
     * @param normal 表面法向量
     * @return 反射速度
     */
    SCDAT::Geometry::Vector3D calculateReflectedVelocity(const SCDAT::Geometry::Vector3D& incident_velocity,
                                                          const SCDAT::Geometry::Vector3D& normal) const;

    /**
     * @brief 应用能量损失
     * @param velocity 速度
     * @return 损失后速度
     */
    SCDAT::Geometry::Vector3D
    applyEnergyLoss(const SCDAT::Geometry::Vector3D& velocity) const;
};

/**
 * @brief 统计吸收边界条件
 *
 * 提供详细的吸收统计分析功能
 */
class AbsorbingBoundaryCondition : public BoundaryConditionBase
{
  public:
    /**
     * @brief 构造函数
     * @param absorption_probability 吸收概率 (0-1)
     */
    explicit AbsorbingBoundaryCondition(double absorption_probability = 1.0);

    /**
     * @brief 析构函数
     */
    ~AbsorbingBoundaryCondition() override = default;

    /**
     * @brief 处理粒子边界相互作用
     */
    std::vector<ParticleTypeDef>
    processParticle(ParticleTypeDef& particle,
                    const SCDAT::Geometry::Point3D& intersection_point,
                    const SCDAT::Geometry::Vector3D& normal, double dt) override;

    /**
     * @brief 设置吸收概率
     */
    void setAbsorptionProbability(double probability);

    /**
     * @brief 获取吸收概率
     */
    double getAbsorptionProbability() const
    {
        return absorption_probability_;
    }

    /**
     * @brief 获取平均吸收能量
     */
    double getAverageAbsorptionEnergy() const;

    /**
     * @brief 获取平均入射角
     */
    double getAverageIncidenceAngle() const;

    /**
     * @brief 获取总吸收能量
     */
    double getTotalAbsorbedEnergy() const
    {
        return total_absorbed_energy_;
    }

    /**
     * @brief 获取总吸收电荷
     */
    double getTotalAbsorbedCharge() const
    {
        return total_absorbed_charge_;
    }

    /**
     * @brief 获取吸收计数
     */
    size_t getAbsorptionCount() const
    {
        return absorption_count_;
    }

    /**
     * @brief 重置吸收统计
     */
    void resetAbsorptionStatistics();

  private:
    double absorption_probability_; ///< 吸收概率
    double total_absorbed_energy_;  ///< 总吸收能量 (J)
    double total_absorbed_charge_;  ///< 总吸收电荷 (C)
    size_t absorption_count_;       ///< 吸收计数

    /**
     * @brief 计算入射角
     * @param velocity 粒子速度
     * @param normal 表面法向量
     * @return 入射角 (弧度)
     */
    double calculateIncidenceAngle(const SCDAT::Geometry::Vector3D& velocity,
                                   const SCDAT::Geometry::Vector3D& normal) const;
};

/**
 * @brief 热发射边界条件
 *
 * 基于Richardson-Dushman方程的热电子发射
 */
class ThermalEmissionBoundaryCondition : public BoundaryConditionBase
{
  public:
    /**
     * @brief 构造函数
     * @param temperature 表面温度 (K)
     * @param work_function 功函数 (eV)
     * @param richardson_constant Richardson常数 (A/(m²·K²))
     */
    ThermalEmissionBoundaryCondition(double temperature = 300.0, double work_function = 4.5,
                                     double richardson_constant = 1.2e6);

    /**
     * @brief 析构函数
     */
    ~ThermalEmissionBoundaryCondition() override = default;

    /**
     * @brief 处理粒子边界相互作用
     */
    std::vector<ParticleTypeDef>
    processParticle(ParticleTypeDef& particle,
                    const SCDAT::Geometry::Point3D& intersection_point,
                    const SCDAT::Geometry::Vector3D& normal, double dt) override;

    /**
     * @brief 设置表面温度
     */
    void setTemperature(double temperature);

    /**
     * @brief 设置发射面积
     */
    void setEmissionArea(double area);

    /**
     * @brief 获取当前发射电流密度
     */
    double getCurrentDensity() const;

    /**
     * @brief 获取表面温度
     */
    double getTemperature() const
    {
        return temperature_;
    }

    /**
     * @brief 获取功函数
     */
    double getWorkFunction() const
    {
        return work_function_;
    }

  private:
    double temperature_;         ///< 表面温度 (K)
    double work_function_;       ///< 功函数 (eV)
    double richardson_constant_; ///< Richardson常数 (A/(m²·K²))
    double emission_area_;       ///< 发射面积 (m²)

    /**
     * @brief 生成热发射电子
     * @param position 发射位置
     * @param normal 表面法向量
     * @return 发射电子
     */
    ParticleTypeDef generateThermalElectron(const SCDAT::Geometry::Point3D& position,
                                            const SCDAT::Geometry::Vector3D& normal);

    /**
     * @brief 生成Maxwell-Boltzmann速度分布
     * @param normal 表面法向量
     * @return 发射速度
     */
    SCDAT::Geometry::Vector3D generateMaxwellBoltzmannVelocity(const SCDAT::Geometry::Vector3D& normal);
};

} // namespace Particle
} // namespace SCDAT

#endif // SCDAT_ADVANCED_BOUNDARY_CONDITIONS_H
