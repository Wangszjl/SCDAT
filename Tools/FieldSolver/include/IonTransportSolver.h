/**
 * @file IonTransportSolver.h
 * @brief 离子输运求解器头文件 / Ion Transport Solver Header
 * @details 定义离子输运方程的数值求解算法接口，包括漂移-扩散方程和碰撞过程
 *          Defines numerical solving algorithm interfaces for ion transport equations, including
 * drift-diffusion equations and collision processes
 *
 * @author PIC-SPIS开发团队 / PIC-SPIS Development Team
 * @date 2025年6月10日 13:45:00
 * @version V1.0.0-enhanced
 * @ingroup FieldSolverModule
 */

#pragma once

#include "../../Basic/include/VoidResult.h"
#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

/**
 * @brief 离子输运参数结构
 */
struct IonTransportParameters
{
    double mobility = 1e-4;              ///< 离子迁移率 (m²/V·s)
    double diffusion_coefficient = 1e-6; ///< 扩散系数 (m²/s)
    double thermal_velocity = 1e3;       ///< 热速度 (m/s)
    double collision_frequency = 1e6;    ///< 碰撞频率 (Hz)
    double temperature = 300.0;          ///< 离子温度 (K)
    double mass = 6.634e-26;             ///< 离子质量 (kg) - 氩离子
    double charge = 1.602e-19;           ///< 离子电荷 (C)
    double atomic_number = 18;           ///< 原子序数
    double ionization_energy = 15.76;    ///< 电离能 (eV)
};

/**
 * @brief 离子状态结构
 */
struct IonState
{
    double density;           ///< 离子密度 (m⁻³)
    Utils::Vector3D velocity; ///< 平均速度 (m/s)
    Utils::Vector3D flux;     ///< 粒子通量 (m⁻²·s⁻¹)
    double temperature;       ///< 离子温度 (K)
    double pressure;          ///< 离子压强 (Pa)
    double charge_state;      ///< 电荷态 (平均)

    IonState() : density(0.0), temperature(300.0), pressure(0.0), charge_state(1.0) {}
};

/**
 * @brief 离子输运系数结构
 */
struct IonTransportCoefficients
{
    double mobility;              ///< 迁移率 (m²/V·s)
    double diffusion_coefficient; ///< 扩散系数 (m²/s)
    double thermal_conductivity;  ///< 热导率 (W/m·K)
    double viscosity;             ///< 粘滞系数 (Pa·s)

    IonTransportCoefficients()
        : mobility(0.0), diffusion_coefficient(0.0), thermal_conductivity(0.0), viscosity(0.0)
    {
    }
};

/**
 * @brief 离子输运统计信息结构
 */
struct IonTransportStatistics
{
    double total_ions = 0.0;         ///< 总离子数
    double average_density = 0.0;    ///< 平均密度 (m⁻³)
    double peak_density = 0.0;       ///< 峰值密度 (m⁻³)
    double average_velocity = 0.0;   ///< 平均速度 (m/s)
    double total_current = 0.0;      ///< 总电流 (A)
    double diffusion_loss = 0.0;     ///< 扩散损失 (s⁻¹)
    double recombination_loss = 0.0; ///< 复合损失 (s⁻¹)
    size_t grid_points = 0;          ///< 网格点数
    int iterations_used = 0;         ///< 使用的迭代次数
};

/**
 * @brief 离子输运求解器类
 * @details 实现离子输运方程的数值求解，包括：
 *          - 离子漂移-扩散方程求解
 *          - 离子-中性粒子碰撞
 *          - 离子-电子复合过程
 *          - 多电荷态离子处理
 */
class IonTransportSolver
{
  public:
    /**
     * @brief 构造函数
     * @details 初始化离子输运求解器，设置默认参数
     */
    IonTransportSolver();

    /**
     * @brief 析构函数
     */
    ~IonTransportSolver();

    /**
     * @brief 设置求解器参数
     * @param params 离子输运参数
     * @return 设置结果
     */
    VoidResult setParameters(const IonTransportParameters& params);

    /**
     * @brief 初始化求解网格
     * @param grid_size 网格尺寸
     * @param resolution 网格分辨率
     * @return 初始化结果
     */
    VoidResult initializeGrid(const Utils::Vector3D& grid_size, const Utils::Vector3D& resolution);

    /**
     * @brief 设置初始条件
     * @param initial_state 初始离子状态
     * @return 设置结果
     */
    VoidResult
    setInitialConditions(const std::function<IonState(const Utils::Point3D&)>& initial_state);

    /**
     * @brief 计算离子输运系数
     * @param state 离子状态
     * @param electric_field 电场 (V/m)
     * @return 输运系数
     * @details 基于动力学理论计算离子迁移率、扩散系数等
     */
    IonTransportCoefficients calculateTransportCoefficients(const IonState& state,
                                                            const Utils::Vector3D& electric_field);

    /**
     * @brief 求解离子漂移-扩散方程
     * @param electric_field 电场分布
     * @param electron_density 电子密度分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 求解连续性方程: ∂n/∂t + ∇·Γ = S
     */
    VoidResult solveDriftDiffusion(const std::vector<Utils::Vector3D>& electric_field,
                                   const std::vector<double>& electron_density, double dt);

    /**
     * @brief 计算离子-电子复合速率
     * @param ion_state 离子状态
     * @param electron_density 电子密度
     * @return 复合速率 (m⁻³·s⁻¹)
     * @details 计算辐射复合和三体复合
     */
    double calculateRecombinationRate(const IonState& ion_state, double electron_density);

    /**
     * @brief 计算离子-中性粒子碰撞
     * @param ion_state 离子状态
     * @param neutral_density 中性粒子密度
     * @return 碰撞频率 (Hz)
     * @details 计算弹性和电荷交换碰撞
     */
    double calculateIonNeutralCollision(const IonState& ion_state, double neutral_density);

    /**
     * @brief 执行一个时间步
     * @param electric_field 电场分布
     * @param electron_density 电子密度分布
     * @param neutral_density 中性粒子密度分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     */
    VoidResult timeStep(const std::vector<Utils::Vector3D>& electric_field,
                        const std::vector<double>& electron_density,
                        const std::vector<double>& neutral_density, double dt);

    /**
     * @brief 获取离子状态
     * @param position 位置坐标
     * @return 离子状态
     */
    IonState getIonState(const Utils::Point3D& position) const;

    /**
     * @brief 获取所有网格点的离子状态
     * @return 离子状态向量
     */
    const std::vector<IonState>& getAllIonStates() const;

    /**
     * @brief 获取离子密度分布
     * @return 密度分布向量
     */
    std::vector<double> getDensityDistribution() const;

    /**
     * @brief 获取离子速度分布
     * @return 速度分布向量
     */
    std::vector<Utils::Vector3D> getVelocityDistribution() const;

    /**
     * @brief 获取离子电流密度分布
     * @return 电流密度分布向量
     */
    std::vector<Utils::Vector3D> getCurrentDensityDistribution() const;

    /**
     * @brief 获取求解器统计信息
     * @return 统计信息字符串
     */
    std::string getStatistics() const;

    /**
     * @brief 检查求解器是否已初始化
     * @return 是否已初始化
     */
    bool isInitialized() const
    {
        return initialized_;
    }

    /**
     * @brief 重置求解器状态
     */
    void reset();

  private:
    // 求解器状态
    bool initialized_; ///< 是否已初始化

    // 求解器参数
    IonTransportParameters params_; ///< 离子输运参数

    // 网格信息
    Utils::Vector3D grid_size_;  ///< 网格尺寸
    Utils::Vector3D resolution_; ///< 网格分辨率
    Utils::Vector3D spacing_;    ///< 网格间距
    size_t total_points_;        ///< 总网格点数

    // 离子状态数据
    std::vector<IonState> ion_states_;     ///< 离子状态
    std::vector<IonState> ion_states_old_; ///< 上一时刻的离子状态

    // 统计信息
    IonTransportStatistics statistics_; ///< 统计信息

    // 私有方法 - 网格操作
    size_t getGridIndex(size_t i, size_t j, size_t k) const;
    Utils::Point3D getGridPosition(size_t i, size_t j, size_t k) const;
    void getGridIndices(size_t index, size_t& i, size_t& j, size_t& k) const;

    // 私有方法 - 数值方法
    double calculateLaplacian(size_t index, const std::vector<double>& field) const;
    Utils::Vector3D calculateGradient(size_t index, const std::vector<double>& field) const;
    double calculateDivergence(size_t index, const std::vector<Utils::Vector3D>& field) const;

    // 私有方法 - 物理计算
    double calculateIonMobility(const IonState& state, double electric_field_magnitude);
    double calculateIonDiffusionCoefficient(const IonState& state);
    double calculateIonCollisionFrequency(const IonState& state, double neutral_density);
    double calculateChargeExchangeRate(const IonState& state, double neutral_density);

    // 私有方法 - 边界条件
    void applyBoundaryConditions();
    bool isBoundaryPoint(size_t i, size_t j, size_t k) const;

    // 私有方法 - 辅助函数
    void updateStatistics();
    double interpolateField(const Utils::Point3D& position, const std::vector<double>& field) const;
};

} // namespace FieldSolver
} // namespace SCDAT
