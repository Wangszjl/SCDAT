/**
 * @file ElectronTransportSolver.h
 * @brief 电子输运求解器头文件
 * @details 定义电子输运方程的数值求解算法接口，包括漂移-扩散方程、能量输运和碰撞积分
 *
 * @author Wang Sizhan
 * @date 2026.3.27 10:52:49
 * @version V0.0.1
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
 * @brief 电子输运参数结构
 */
struct ElectronTransportParameters
{
    double mobility = 1e-2;                  ///< 电子迁移率 (m²/V·s)
    double diffusion_coefficient = 1e-4;     ///< 扩散系数 (m²/s)
    double thermal_velocity = 1e5;           ///< 热速度 (m/s)
    double collision_frequency = 1e8;        ///< 碰撞频率 (Hz)
    double energy_relaxation_time = 1e-9;    ///< 能量弛豫时间 (s)
    double momentum_relaxation_time = 1e-12; ///< 动量弛豫时间 (s)
    double temperature = 300.0;              ///< 电子温度 (K)
    double mass = 9.109e-31;                 ///< 电子质量 (kg)
    double charge = -1.602e-19;              ///< 电子电荷 (C)
};

/**
 * @brief 电子状态结构
 */
struct ElectronState
{
    double density;              ///< 电子密度 (m⁻³)
    double energy_density;       ///< 能量密度 (J/m³)
    SCDAT::Geometry::Vector3D velocity;    ///< 平均速度 (m/s)
    SCDAT::Geometry::Vector3D flux;        ///< 粒子通量 (m⁻²·s⁻¹)
    SCDAT::Geometry::Vector3D energy_flux; ///< 能量通量 (W/m²)
    double temperature;          ///< 电子温度 (K)
    double pressure;             ///< 电子压强 (Pa)

    ElectronState() : density(0.0), energy_density(0.0), temperature(300.0), pressure(0.0) {}
};

/**
 * @brief 输运系数结构
 */
struct TransportCoefficients
{
    double mobility;              ///< 迁移率 (m²/V·s)
    double diffusion_coefficient; ///< 扩散系数 (m²/s)
    double thermal_conductivity;  ///< 热导率 (W/m·K)
    double energy_mobility;       ///< 能量迁移率 (m²/V·s)
    double energy_diffusion;      ///< 能量扩散系数 (m²/s)

    TransportCoefficients()
        : mobility(0.0), diffusion_coefficient(0.0), thermal_conductivity(0.0),
          energy_mobility(0.0), energy_diffusion(0.0)
    {
    }
};

/**
 * @brief 求解器配置结构
 */
struct SolverConfiguration
{
    double time_step = 1e-12;                     ///< 时间步长 (s)
    double spatial_step = 1e-6;                   ///< 空间步长 (m)
    double tolerance = 1e-10;                     ///< 收敛精度
    int max_iterations = 1000;                    ///< 最大迭代次数
    bool use_implicit_scheme = true;              ///< 是否使用隐式格式
    bool enable_energy_equation = true;           ///< 是否求解能量方程
    bool enable_collision_terms = true;           ///< 是否包含碰撞项
    std::string boundary_condition = "dirichlet"; ///< 边界条件类型
};

/**
 * @brief 电子输运统计信息结构
 */
struct ElectronTransportStatistics
{
    double total_electrons = 0.0;   ///< 总电子数
    double average_density = 0.0;   ///< 平均密度 (m⁻³)
    double peak_density = 0.0;      ///< 峰值密度 (m⁻³)
    double average_energy = 0.0;    ///< 平均能量 (eV)
    double total_current = 0.0;     ///< 总电流 (A)
    double power_dissipation = 0.0; ///< 功率耗散 (W)
    double diffusion_loss = 0.0;    ///< 扩散损失 (s⁻¹)
    double collision_loss = 0.0;    ///< 碰撞损失 (s⁻¹)
    size_t grid_points = 0;         ///< 网格点数
    int iterations_used = 0;        ///< 使用的迭代次数
};

/**
 * @brief 电子输运求解器类
 * @details 实现电子输运方程的数值求解，包括：
 *          - 漂移-扩散方程求解 (连续性方程 + 动量方程)
 *          - 能量输运方程求解 (能量平衡方程)
 *          - 碰撞积分计算 (弹性和非弹性碰撞)
 *          - 输运系数计算 (迁移率、扩散系数等)
 */
class ElectronTransportSolver
{
  public:
    /**
     * @brief 构造函数
     * @details 初始化电子输运求解器，设置默认参数
     */
    ElectronTransportSolver();

    /**
     * @brief 析构函数
     */
    ~ElectronTransportSolver();

    /**
     * @brief 设置求解器参数
     * @param params 电子输运参数
     * @return 设置结果
     */
    VoidResult setParameters(const ElectronTransportParameters& params);

    /**
     * @brief 设置求解器配置
     * @param config 求解器配置
     * @return 设置结果
     */
    VoidResult setConfiguration(const SolverConfiguration& config);

    /**
     * @brief 初始化求解网格
     * @param grid_size 网格尺寸
     * @param resolution 网格分辨率
     * @return 初始化结果
     */
    VoidResult initializeGrid(const Geometry::Vector3D& grid_size, const Geometry::Vector3D& resolution);

    /**
     * @brief 设置初始条件
     * @param initial_state 初始电子状态
     * @return 设置结果
     */
    VoidResult
    setInitialConditions(const std::function<ElectronState(const Geometry::Point3D&)>& initial_state);

    /**
     * @brief 设置边界条件
     * @param boundary_state 边界电子状态函数
     * @return 设置结果
     */
    VoidResult setBoundaryConditions(
        const std::function<ElectronState(const Geometry::Point3D&, double)>& boundary_state);

    /**
     * @brief 计算输运系数
     * @param state 电子状态
     * @param electric_field 电场 (V/m)
     * @return 输运系数
     * @details 基于动力学理论计算迁移率、扩散系数等
     */
    TransportCoefficients calculateTransportCoefficients(const ElectronState& state,
                                                         const Geometry::Vector3D& electric_field);

    /**
     * @brief 求解漂移-扩散方程
     * @param electric_field 电场分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 求解连续性方程: ∂n/∂t + ∇·Γ = S
     */
    VoidResult solveDriftDiffusion(const std::vector<Geometry::Vector3D>& electric_field, double dt);

    /**
     * @brief 求解能量输运方程
     * @param electric_field 电场分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 求解能量方程: ∂(nε)/∂t + ∇·Γε = Sε
     */
    VoidResult solveEnergyTransport(const std::vector<Geometry::Vector3D>& electric_field, double dt);

    /**
     * @brief 计算碰撞积分
     * @param state 电子状态
     * @return 碰撞源项 (m⁻³·s⁻¹)
     * @details 计算弹性和非弹性碰撞的源项
     */
    double calculateCollisionIntegral(const ElectronState& state);

    /**
     * @brief 执行一个时间步
     * @param electric_field 电场分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 综合求解漂移-扩散和能量输运方程
     */
    VoidResult timeStep(const std::vector<Geometry::Vector3D>& electric_field, double dt);

    /**
     * @brief 获取电子状态
     * @param position 位置坐标
     * @return 电子状态
     */
    ElectronState getElectronState(const Geometry::Point3D& position) const;

    /**
     * @brief 获取所有网格点的电子状态
     * @return 电子状态向量
     */
    const std::vector<ElectronState>& getAllElectronStates() const;

    /**
     * @brief 获取电子密度分布
     * @return 密度分布向量
     */
    std::vector<double> getDensityDistribution() const;

    /**
     * @brief 获取电子温度分布
     * @return 温度分布向量
     */
    std::vector<double> getTemperatureDistribution() const;

    /**
     * @brief 获取电流密度分布
     * @param electric_field 电场分布
     * @return 电流密度分布向量
     */
    std::vector<Geometry::Vector3D>
    getCurrentDensityDistribution(const std::vector<Geometry::Vector3D>& electric_field) const;

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

    // 求解器参数和配置
    ElectronTransportParameters params_; ///< 电子输运参数
    SolverConfiguration config_;         ///< 求解器配置

    // 网格信息
    Geometry::Vector3D grid_size_;  ///< 网格尺寸
    Geometry::Vector3D resolution_; ///< 网格分辨率
    Geometry::Vector3D spacing_;    ///< 网格间距
    size_t total_points_;        ///< 总网格点数

    // 电子状态数据
    std::vector<ElectronState> electron_states_;     ///< 电子状态
    std::vector<ElectronState> electron_states_old_; ///< 上一时刻的电子状态

    // 边界条件函数
    std::function<ElectronState(const Geometry::Point3D&, double)> boundary_function_;

    // 统计信息
    ElectronTransportStatistics statistics_; ///< 统计信息

    // 私有方法 - 网格操作
    size_t getGridIndex(size_t i, size_t j, size_t k) const;
    Geometry::Point3D getGridPosition(size_t i, size_t j, size_t k) const;
    void getGridIndices(size_t index, size_t& i, size_t& j, size_t& k) const;

    // 私有方法 - 数值方法
    VoidResult solveExplicitScheme(const std::vector<Geometry::Vector3D>& electric_field, double dt);
    VoidResult solveImplicitScheme(const std::vector<Geometry::Vector3D>& electric_field, double dt);
    double calculateLaplacian(size_t index, const std::vector<double>& field) const;
    Geometry::Vector3D calculateGradient(size_t index, const std::vector<double>& field) const;
    double calculateDivergence(size_t index, const std::vector<Geometry::Vector3D>& field) const;

    // 私有方法 - 物理计算
    double calculateMobility(const ElectronState& state, double electric_field_magnitude);
    double calculateDiffusionCoefficient(const ElectronState& state);
    double calculateThermalConductivity(const ElectronState& state);
    double calculateCollisionFrequency(const ElectronState& state);

    // 私有方法 - 边界条件
    void applyBoundaryConditions(double time);
    bool isBoundaryPoint(size_t i, size_t j, size_t k) const;

    // 私有方法 - 辅助函数
    void updateStatistics();
    double interpolateField(const Geometry::Point3D& position, const std::vector<double>& field) const;
    Geometry::Vector3D interpolateVectorField(const Geometry::Point3D& position,
                                           const std::vector<Geometry::Vector3D>& field) const;
};

} // namespace FieldSolver
} // namespace SCDAT
