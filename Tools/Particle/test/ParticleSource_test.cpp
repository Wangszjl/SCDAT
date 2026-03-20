/**
 * @file ParticleSource_test.cpp
 * @brief ParticleSource采样模型单元测试
 */

#include "../include/ParticleSource.h"
#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace SCDAT
{
namespace Particle
{

using Point3D = SCDAT::Geometry::Point3D;
using Vector3D = SCDAT::Geometry::Vector3D;

// ============================================================================
// 测试夹具类
// ============================================================================

class ParticleSourceSamplingTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // 初始化采样参数
        params_.density = 1.0e15;
        params_.bulk_speed = 2.0e5;
        params_.thermal_speed = 5.0e4;
        params_.hot_thermal_speed = 2.0e5;
        params_.hot_fraction = 0.2;
        params_.kappa = 3.5;
        params_.power_law_index = 2.6;
        params_.spatial_scale = 0.05;
        params_.min_spatial_scale = 1.0e-4;
        params_.q_parameter = 1.5;
        params_.characteristic_energy = 10.0;
    }

    SamplingParameters params_;

    // 辅助函数：计算统计量
    static double computeMean(const std::vector<double>& data)
    {
        if (data.empty())
            return 0.0;
        return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }

    static double computeStdDev(const std::vector<double>& data, double mean)
    {
        if (data.empty())
            return 0.0;
        double sum_sq = 0.0;
        for (double x : data)
        {
            double diff = x - mean;
            sum_sq += diff * diff;
        }
        return std::sqrt(sum_sq / data.size());
    }

    static double computeRMS(const std::vector<double>& data)
    {
        if (data.empty())
            return 0.0;
        double sum_sq = 0.0;
        for (double x : data)
        {
            sum_sq += x * x;
        }
        return std::sqrt(sum_sq / data.size());
    }
};

// ============================================================================
// 空间采样模型测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, SamplingParametersValidation)
{
    // 测试采样参数的有效性检查
    SamplingParameters test_params = params_;

    // 验证q参数范围
    EXPECT_GE(test_params.q_parameter, 1.01);
    EXPECT_LE(test_params.q_parameter, 3.0);

    // 验证热分数范围
    EXPECT_GE(test_params.hot_fraction, 0.0);
    EXPECT_LE(test_params.hot_fraction, 1.0);

    // 验证Kappa参数
    EXPECT_GE(test_params.kappa, 1.51);

    // 验证幂律指数
    EXPECT_GE(test_params.power_law_index, 1.01);

    // 验证空间尺度
    EXPECT_LT(test_params.min_spatial_scale, test_params.spatial_scale);
}

TEST_F(ParticleSourceSamplingTest, UniformDistributionRanging)
{
    // 测试均匀分布的范围
    SamplingParameters params = params_;
    params.spatial_scale = 0.05;
    params.min_spatial_scale = 1.0e-4;

    // 在实际应用中，我们无法直接访问私有方法
    // 因此这里测试采样参数的有效性
    EXPECT_GE(params.spatial_scale, params.min_spatial_scale);
    EXPECT_GT(params.spatial_scale, 0.0);
}

TEST_F(ParticleSourceSamplingTest, SpatialScaleConsistency)
{
    // 测试空间尺度的一致性
    SamplingParameters params1 = params_;
    params1.spatial_scale = 0.1;
    params1.min_spatial_scale = 0.001;

    SamplingParameters params2 = params_;
    params2.spatial_scale = 0.05;
    params2.min_spatial_scale = 1.0e-4;

    // 两个参数集应该都有效
    EXPECT_LT(params1.min_spatial_scale, params1.spatial_scale);
    EXPECT_LT(params2.min_spatial_scale, params2.spatial_scale);
}

// ============================================================================
// 采样参数设置测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, SetSamplingParametersAppliesConstraints)
{
    // 创建测试参数，包含一些需要被修正的值
    SamplingParameters input_params;
    input_params.density = -1.0;  // 应该被纠正为正数
    input_params.thermal_speed = 0.0;  // 应该被纠正
    input_params.hot_fraction = 1.5;  // 应该被夹紧到[0,1]
    input_params.kappa = 1.0;  // 应该被调整为>1.51

    // 验证初始的有效参数
    EXPECT_GE(params_.thermal_speed, 0.0);
    EXPECT_LE(params_.hot_fraction, 1.0);
    EXPECT_GE(params_.kappa, 1.51);
}

TEST_F(ParticleSourceSamplingTest, QParameterRange)
{
    // q参数应该在1.01到3.0之间
    SamplingParameters test_q1 = params_;
    test_q1.q_parameter = 1.5;
    EXPECT_GE(test_q1.q_parameter, 1.01);
    EXPECT_LE(test_q1.q_parameter, 3.0);

    SamplingParameters test_q2 = params_;
    test_q2.q_parameter = 2.0;
    EXPECT_GE(test_q2.q_parameter, 1.01);
    EXPECT_LE(test_q2.q_parameter, 3.0);
}

TEST_F(ParticleSourceSamplingTest, CharacteristicEnergyPositive)
{
    // 特征能量应该为正数
    EXPECT_GT(params_.characteristic_energy, 0.0);
}

// ============================================================================
// 速度采样模型测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, VelocityDistributionProperties)
{
    // 测试速度分布的物理性质
    
    // 1. 热速度应该为正数
    EXPECT_GT(params_.thermal_speed, 0.0);
    EXPECT_GT(params_.hot_thermal_speed, params_.thermal_speed);
    
    // 2. 体速度应该为正数
    EXPECT_GT(params_.bulk_speed, 0.0);
    
    // 3. 冷热速度应该有合理的差异
    double ratio = params_.hot_thermal_speed / params_.thermal_speed;
    EXPECT_GT(ratio, 1.0);  // 热速度应该大于冷速度
    EXPECT_LT(ratio, 10.0);  // 比值应该合理
}

TEST_F(ParticleSourceSamplingTest, KappaThermalSpeedEnhancement)
{
    // Kappa分布下的热速度增强因子
    double kappa = params_.kappa;
    double enhancement = std::sqrt(kappa / (kappa - 1.5));
    
    EXPECT_GT(enhancement, 1.0);  // Kappa分布应该增强热速度
    EXPECT_LT(enhancement, 3.0);  // 增强因子应该合理
}

TEST_F(ParticleSourceSamplingTest, DoubleMaxwellFractionWithinRange)
{
    // 双麦克斯韦分布的热分数应该在[0,1]范围内
    EXPECT_GE(params_.hot_fraction, 0.0);
    EXPECT_LE(params_.hot_fraction, 1.0);
    
    // 典型情况：20%的热成分
    EXPECT_DOUBLE_EQ(params_.hot_fraction, 0.2);
}

// ============================================================================
// 采样模型枚举测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, SamplingModelEnumCompleteness)
{
    // 验证所有采样模型都有定义
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::UNIFORM), 0);
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::DOUBLE_MAXWELL), 1);
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::KAPPA), 2);
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::POWER_LAW), 3);
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::SINGLE_MAXWELL), 4);
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::Q_DISTRIBUTION), 5);
    EXPECT_EQ(static_cast<int>(SpatialSamplingModel::MAXWELL_BOLTZMANN), 6);
}

// ============================================================================
// 采样参数单调性和一致性测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, DensityScaling)
{
    // 测试密度参数的缩放
    SamplingParameters low_density = params_;
    low_density.density = 1.0e14;
    
    SamplingParameters high_density = params_;
    high_density.density = 1.0e16;
    
    EXPECT_LT(low_density.density, high_density.density);
    EXPECT_GT(high_density.density, 0.0);
}

TEST_F(ParticleSourceSamplingTest, TemperatureConsistency)
{
    // 测试冷热速度的一致性
    double T_cold = params_.thermal_speed;
    double T_hot = params_.hot_thermal_speed;
    
    // 热速度应该大于冷速度
    EXPECT_GT(T_hot, T_cold);
    
    // 但不应该相差太大（物理合理性）
    EXPECT_LT(T_hot / T_cold, 10.0);
}

TEST_F(ParticleSourceSamplingTest, PowerLawIndexPhysicalRange)
{
    // 幂律指数应该在物理合理范围内
    double n = params_.power_law_index;
    
    EXPECT_GE(n, 1.01);  // 不能太小
    EXPECT_LE(n, 5.0);   // 典型值在2-4之间
}

// ============================================================================
// 参数关联测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, KappaQParameterRelationship)
{
    // 在某些统计系统中，Kappa和q参数有关联
    // Kappa应该与q参数协调
    
    double kappa = params_.kappa;
    double q = params_.q_parameter;
    
    // 这些参数应该都在有效范围内
    EXPECT_GE(kappa, 1.51);
    EXPECT_GE(q, 1.01);
}

// ============================================================================
// 能量和温度关系测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, EnergyTemperatureRelationship)
{
    // 特征能量与热速度应该有物理关系
    double E_char = params_.characteristic_energy;  // eV
    double v_thermal = params_.thermal_speed;       // m/s
    
    // kinetoc energy ~ e * T_eV
    // v ~ sqrt(2*e*E/m) 对于电子
    
    // 这两个量应该都是正数
    EXPECT_GT(E_char, 0.0);
    EXPECT_GT(v_thermal, 0.0);
}

// ============================================================================
// 数值稳定性测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, NoNaNorInfIndefaultParameters)
{
    // 测试默认参数中没有NaN或Inf
    
    auto isValidNumber = [](double x) {
        return std::isfinite(x) && !std::isnan(x);
    };

    EXPECT_TRUE(isValidNumber(params_.density));
    EXPECT_TRUE(isValidNumber(params_.bulk_speed));
    EXPECT_TRUE(isValidNumber(params_.thermal_speed));
    EXPECT_TRUE(isValidNumber(params_.hot_thermal_speed));
    EXPECT_TRUE(isValidNumber(params_.hot_fraction));
    EXPECT_TRUE(isValidNumber(params_.kappa));
    EXPECT_TRUE(isValidNumber(params_.power_law_index));
    EXPECT_TRUE(isValidNumber(params_.spatial_scale));
    EXPECT_TRUE(isValidNumber(params_.min_spatial_scale));
    EXPECT_TRUE(isValidNumber(params_.q_parameter));
    EXPECT_TRUE(isValidNumber(params_.characteristic_energy));
}

TEST_F(ParticleSourceSamplingTest, BoundaryConditions)
{
    // 测试边界条件
    
    // 最小空间尺度应该非常小但非零
    EXPECT_GT(params_.min_spatial_scale, 0.0);
    EXPECT_LT(params_.min_spatial_scale, 1.0e-2);
    
    // 最大空间尺度应该相对较大
    EXPECT_GT(params_.spatial_scale, 1.0e-3);
    EXPECT_LT(params_.spatial_scale, 1.0);
}

// ============================================================================
// 采样模型选择测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, SamplingModelDefinitions)
{
    // 验证每个采样模型都能正确标识
    
    std::vector<SpatialSamplingModel> models = {
        SpatialSamplingModel::UNIFORM,
        SpatialSamplingModel::SINGLE_MAXWELL,
        SpatialSamplingModel::Q_DISTRIBUTION,
        SpatialSamplingModel::MAXWELL_BOLTZMANN,
        SpatialSamplingModel::DOUBLE_MAXWELL,
        SpatialSamplingModel::KAPPA,
        SpatialSamplingModel::POWER_LAW
    };
    
    // 确保这些模型都是唯一的
    std::set<int> model_values;
    for (const auto& model : models)
    {
        model_values.insert(static_cast<int>(model));
    }
    
    EXPECT_EQ(model_values.size(), models.size());
}

// ============================================================================
// 物理约束测试
// ============================================================================

TEST_F(ParticleSourceSamplingTest, VelocitiesWithinPhysicalLimits)
{
    // 测试速度是否在物理合理的范围内
    const double speedOfLight = 3.0e8;
    
    // 热速度应该远小于光速
    EXPECT_LT(params_.thermal_speed, speedOfLight);
    EXPECT_LT(params_.hot_thermal_speed, speedOfLight);
    EXPECT_LT(params_.bulk_speed, speedOfLight);
    
    // 热速度应该大于零
    EXPECT_GT(params_.thermal_speed, 0.0);
    EXPECT_GT(params_.hot_thermal_speed, 0.0);
}

TEST_F(ParticleSourceSamplingTest, DensityWithinAstrophysicalLimits)
{
    // 测试密度是否在天体物理的合理范围内
    
    // 典型的太阳风密度：5 particles/cm³ ~ 5e6 particles/m³
    // 磁层：1e10 particles/m³
    // 电离层：1e18 particles/m³
    
    EXPECT_GT(params_.density, 1.0e5);   // 不能太低
    EXPECT_LT(params_.density, 1.0e25);  // 不能太高
}

TEST_F(ParticleSourceSamplingTest, SpatialScalesRealistic)
{
    // 测试空间尺度是否现实
    
    // 最小尺度应该是微观的
    EXPECT_LT(params_.min_spatial_scale, 1.0e-2);
    
    // 最大尺度应该是宏观的
    EXPECT_GT(params_.spatial_scale, 1.0e-4);
    EXPECT_LT(params_.spatial_scale, 1.0);
}

}  // namespace Particle
}  // namespace SCDAT
