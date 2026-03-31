/**
 * @file VolDistrib.h
 * @brief 体积分布模块头文件
 * @details 定义粒子密度分布、电荷密度分布、能量密度分布的计算方法
 *
 * @author Wang Sizhan
 * @date 2026.3.26
 * @version V0.0.1
 * @ingroup FieldSolverModule
 */

#ifndef SCDAT_FIELD_VOL_DISTRIB_H
#define SCDAT_FIELD_VOL_DISTRIB_H

#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
#include "../Particle/include/ParticleDefinitions.h"
#include <functional>
#include <memory>
#include <vector>

namespace SCDAT
{
namespace Field
{

// 类型别名，简化使用
using Point3D = Geometry::Point3D;
using Vector3D = Geometry::Vector3D;
using ParticleType = Particle::ParticleType;
using ParticlePtr = Particle::ParticlePtr;

/**
 * @brief 简化的标量场类
 * @details 专门为VolDistrib模块设计的轻量级标量场
 */
class ScalarField
{
  public:
    explicit ScalarField(size_t size) : values_(size, 0.0) {}

    void setValue(size_t index, double value)
    {
        if (index < values_.size())
            values_[index] = value;
    }

    double getValue(size_t index) const
    {
        return index < values_.size() ? values_[index] : 0.0;
    }

    const std::vector<double>& getValues() const
    {
        return values_;
    }
    void setValues(const std::vector<double>& values)
    {
        values_ = values;
    }
    void reset()
    {
        std::fill(values_.begin(), values_.end(), 0.0);
    }
    size_t size() const
    {
        return values_.size();
    }

  private:
    std::vector<double> values_;
};

/**
 * @brief 简化的向量场类
 * @details 专门为VolDistrib模块设计的轻量级向量场
 */
class VectorField
{
  public:
    explicit VectorField(size_t size) : values_(size, Vector3D(0.0, 0.0, 0.0)) {}

    void setValue(size_t index, const Vector3D& value)
    {
        if (index < values_.size())
            values_[index] = value;
    }

    Vector3D getValue(size_t index) const
    {
        return index < values_.size() ? values_[index] : Vector3D(0.0, 0.0, 0.0);
    }

    const std::vector<Vector3D>& getValues() const
    {
        return values_;
    }
    void setValues(const std::vector<Vector3D>& values)
    {
        values_ = values;
    }
    void reset()
    {
        std::fill(values_.begin(), values_.end(), Vector3D(0.0, 0.0, 0.0));
    }
    size_t size() const
    {
        return values_.size();
    }

  private:
    std::vector<Vector3D> values_;
};

/**
 * @brief 体积分布基类
 * @details 定义所有体积分布的通用接口
 */
class VolDistrib
{
  public:
    /**
     * @brief 构造函数
     * @param nodeCount 节点数量
     * @param particleType 粒子类型
     */
    VolDistrib(size_t nodeCount, ParticleType particleType);

    /**
     * @brief 虚析构函数
     */
    virtual ~VolDistrib();

    /**
     * @brief 获取节点数量
     * @return 节点数量
     */
    size_t getNodeCount() const;

    /**
     * @brief 获取粒子类型
     * @return 粒子类型
     */
    ParticleType getParticleType() const;

  protected:
    size_t nodeCount_;          ///< 节点数量
    ParticleType particleType_; ///< 粒子类型
};

/**
 * @brief 粒子密度分布类
 * @details 处理粒子数密度的空间分布
 */
class ParticleDensityDistrib : public VolDistrib
{
  public:
    /**
     * @brief 构造函数
     * @param nodeCount 节点数量
     * @param particleType 粒子类型
     */
    ParticleDensityDistrib(size_t nodeCount, ParticleType particleType);

    /**
     * @brief 设置均匀密度
     * @param density 密度值 [m^-3]
     */
    void setUniformDensity(double density);

    /**
     * @brief 设置密度分布函数
     * @param profile 密度分布函数
     */
    void setDensityProfile(std::function<double(size_t)> profile);

    /**
     * @brief 添加高斯分布
     * @param centerIndex 中心节点索引
     * @param amplitude 幅度
     * @param sigma 标准差
     */
    void addGaussianDistribution(size_t centerIndex, double amplitude, double sigma);

    /**
     * @brief 获取指定节点的密度
     * @param nodeIndex 节点索引
     * @return 密度值
     */
    double getDensity(size_t nodeIndex) const;

    /**
     * @brief 获取密度场
     * @return 密度场的常引用
     */
    const ScalarField& getDensityField() const;

    /**
     * @brief 计算矩
     * @param result 结果场
     * @param order 矩的阶数
     */
    void computeMoment(ScalarField& result, int order) const;

    /**
     * @brief 从粒子更新密度分布
     * @param particles 粒子列表
     */
    void updateFromParticles(const std::vector<ParticlePtr>& particles);

  private:
    ScalarField density_; ///< 密度场
};

/**
 * @brief 电荷密度分布类
 * @details 处理电荷密度的空间分布
 */
class ChargeDensityDistrib : public VolDistrib
{
  public:
    /**
     * @brief 构造函数
     * @param nodeCount 节点数量
     * @param particleType 粒子类型
     */
    ChargeDensityDistrib(size_t nodeCount, ParticleType particleType);

    /**
     * @brief 从粒子密度计算电荷密度
     * @param particleDensity 粒子密度分布
     */
    void computeFromParticleDensity(const ParticleDensityDistrib& particleDensity);

    /**
     * @brief 添加电荷密度贡献
     * @param particleDensity 粒子密度分布
     * @param chargeMultiplier 电荷倍数
     */
    void addContribution(const ParticleDensityDistrib& particleDensity,
                         double chargeMultiplier = 1.0);

    /**
     * @brief 获取指定节点的电荷密度
     * @param nodeIndex 节点索引
     * @return 电荷密度值
     */
    double getChargeDensity(size_t nodeIndex) const;

    /**
     * @brief 获取电荷密度场
     * @return 电荷密度场的常引用
     */
    const ScalarField& getChargeDensityField() const;

  private:
    ScalarField chargeDensity_; ///< 电荷密度场
};

/**
 * @brief 能量密度分布类
 * @details 处理能量密度的空间分布
 */
class EnergyDensityDistrib : public VolDistrib
{
  public:
    /**
     * @brief 构造函数
     * @param nodeCount 节点数量
     * @param particleType 粒子类型
     */
    EnergyDensityDistrib(size_t nodeCount, ParticleType particleType);

    /**
     * @brief 计算动能密度
     * @param particleDensity 粒子密度分布
     * @param velocityField 速度场
     */
    void computeKineticEnergy(const ParticleDensityDistrib& particleDensity,
                              const VectorField& velocityField);

    /**
     * @brief 获取指定节点的能量密度
     * @param nodeIndex 节点索引
     * @return 能量密度值
     */
    double getEnergyDensity(size_t nodeIndex) const;

    /**
     * @brief 获取能量密度场
     * @return 能量密度场的常引用
     */
    const ScalarField& getEnergyDensityField() const;

    /**
     * @brief 从粒子更新能量密度分布
     * @param particles 粒子列表
     */
    void updateFromParticles(const std::vector<ParticlePtr>& particles);

  private:
    ScalarField energyDensity_; ///< 能量密度场
};

} // namespace Field
} // namespace SCDAT

#endif // SCDAT_FIELD_VOL_DISTRIB_H
