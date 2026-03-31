/**
 * @file EnergyBalanceSolver.h
 * @brief 能量平衡求解器头文件 / Energy Balance Solver Header
 * @details 定义能量平衡方程的数值求解算法接口，包括电子和离子的能量输运
 *          Defines numerical solving algorithm interfaces for energy balance equations, including
 * electron and ion energy transport
 *
 * @author PIC-SPIS开发团队 / PIC-SPIS Development Team
 * @date 2025年6月10日 13:50:00
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
 * @brief 能量平衡参数结构
 */
struct EnergyBalanceParameters
{
    double electron_heat_capacity = 1.5;         ///< 电子热容 (kB)
    double ion_heat_capacity = 1.5;              ///< 离子热容 (kB)
    double electron_thermal_conductivity = 1e-2; ///< 电子热导率 (W/m·K)
    double ion_thermal_conductivity = 1e-4;      ///< 离子热导率 (W/m·K)
    double elastic_collision_frequency = 1e8;    ///< 弹性碰撞频率 (Hz)
    double inelastic_collision_frequency = 1e6;  ///< 非弹性碰撞频率 (Hz)
    double energy_exchange_coefficient = 1e-20;  ///< 能量交换系数 (m³/s)
    double radiation_coefficient = 1e-40;        ///< 辐射系数 (W·m³)
};

/**
 * @brief 能量状态结构
 */
struct EnergyState
{
    double electron_temperature;        ///< 电子温度 (K)
    double ion_temperature;             ///< 离子温度 (K)
    double electron_energy_density;     ///< 电子能量密度 (J/m³)
    double ion_energy_density;          ///< 离子能量密度 (J/m³)
    Utils::Vector3D electron_heat_flux; ///< 电子热流 (W/m²)
    Utils::Vector3D ion_heat_flux;      ///< 离子热流 (W/m²)
    double joule_heating_rate;          ///< 焦耳加热率 (W/m³)
    double elastic_cooling_rate;        ///< 弹性冷却率 (W/m³)
    double inelastic_cooling_rate;      ///< 非弹性冷却率 (W/m³)
    double radiation_loss_rate;         ///< 辐射损失率 (W/m³)

    EnergyState()
        : electron_temperature(300.0), ion_temperature(300.0), electron_energy_density(0.0),
          ion_energy_density(0.0), joule_heating_rate(0.0), elastic_cooling_rate(0.0),
          inelastic_cooling_rate(0.0), radiation_loss_rate(0.0)
    {
    }
};

/**
 * @brief 能量源项结构
 */
struct EnergySourceTerms
{
    double electron_joule_heating;  ///< 电子焦耳加热 (W/m³)
    double ion_joule_heating;       ///< 离子焦耳加热 (W/m³)
    double elastic_energy_exchange; ///< 弹性能量交换 (W/m³)
    double inelastic_energy_loss;   ///< 非弹性能量损失 (W/m³)
    double radiation_loss;          ///< 辐射损失 (W/m³)
    double thermal_conduction_loss; ///< 热传导损失 (W/m³)

    EnergySourceTerms()
        : electron_joule_heating(0.0), ion_joule_heating(0.0), elastic_energy_exchange(0.0),
          inelastic_energy_loss(0.0), radiation_loss(0.0), thermal_conduction_loss(0.0)
    {
    }
};

/**
 * @brief 能量平衡统计信息结构
 */
struct EnergyBalanceStatistics
{
    double total_electron_energy = 0.0; ///< 总电子能量 (J)
    double total_ion_energy = 0.0;      ///< 总离子能量 (J)
    double average_electron_temp = 0.0; ///< 平均电子温度 (K)
    double average_ion_temp = 0.0;      ///< 平均离子温度 (K)
    double peak_electron_temp = 0.0;    ///< 峰值电子温度 (K)
    double peak_ion_temp = 0.0;         ///< 峰值离子温度 (K)
    double total_joule_heating = 0.0;   ///< 总焦耳加热 (W)
    double total_radiation_loss = 0.0;  ///< 总辐射损失 (W)
    double total_conduction_loss = 0.0; ///< 总传导损失 (W)
    double energy_balance_error = 0.0;  ///< 能量平衡误差 (%)
    size_t grid_points = 0;             ///< 网格点数
    int iterations_used = 0;            ///< 使用的迭代次数
};

/**
 * @brief 能量平衡求解器类
 * @details 实现能量平衡方程的数值求解，包括：
 *          - 电子能量平衡方程求解
 *          - 离子能量平衡方程求解
 *          - 能量源项计算 (焦耳加热、碰撞、辐射)
 *          - 热传导方程求解
 */
class EnergyBalanceSolver
{
  public:
    /**
     * @brief 构造函数
     * @details 初始化能量平衡求解器，设置默认参数
     */
    EnergyBalanceSolver();

    /**
     * @brief 析构函数
     */
    ~EnergyBalanceSolver();

    /**
     * @brief 设置求解器参数
     * @param params 能量平衡参数
     * @return 设置结果
     */
    VoidResult setParameters(const EnergyBalanceParameters& params);

    /**
     * @brief 初始化求解网格
     * @param grid_size 网格尺寸
     * @param resolution 网格分辨率
     * @return 初始化结果
     */
    VoidResult initializeGrid(const Utils::Vector3D& grid_size, const Utils::Vector3D& resolution);

    /**
     * @brief 设置初始条件
     * @param initial_state 初始能量状态
     * @return 设置结果
     */
    VoidResult
    setInitialConditions(const std::function<EnergyState(const Utils::Point3D&)>& initial_state);

    /**
     * @brief 计算能量源项
     * @param position 位置坐标
     * @param electron_density 电子密度 (m⁻³)
     * @param ion_density 离子密度 (m⁻³)
     * @param electric_field 电场 (V/m)
     * @param current_density 电流密度 (A/m²)
     * @return 能量源项
     * @details 计算焦耳加热、碰撞冷却、辐射损失等源项
     */
    EnergySourceTerms calculateEnergySourceTerms(const Utils::Point3D& position,
                                                 double electron_density, double ion_density,
                                                 const Utils::Vector3D& electric_field,
                                                 const Utils::Vector3D& current_density);

    /**
     * @brief 求解电子能量平衡方程
     * @param electron_density 电子密度分布
     * @param electric_field 电场分布
     * @param current_density 电流密度分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 求解: ∂(nεe)/∂t + ∇·qe = Se
     */
    VoidResult solveElectronEnergyBalance(const std::vector<double>& electron_density,
                                          const std::vector<Utils::Vector3D>& electric_field,
                                          const std::vector<Utils::Vector3D>& current_density,
                                          double dt);

    /**
     * @brief 求解离子能量平衡方程
     * @param ion_density 离子密度分布
     * @param electric_field 电场分布
     * @param current_density 电流密度分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 求解: ∂(nεi)/∂t + ∇·qi = Si
     */
    VoidResult solveIonEnergyBalance(const std::vector<double>& ion_density,
                                     const std::vector<Utils::Vector3D>& electric_field,
                                     const std::vector<Utils::Vector3D>& current_density,
                                     double dt);

    /**
     * @brief 计算热传导
     * @param temperature_field 温度场
     * @param thermal_conductivity 热导率
     * @return 热传导项 (W/m³)
     * @details 计算: -∇·(κ∇T)
     */
    std::vector<double> calculateThermalConduction(const std::vector<double>& temperature_field,
                                                   double thermal_conductivity);

    /**
     * @brief 执行一个时间步
     * @param electron_density 电子密度分布
     * @param ion_density 离子密度分布
     * @param electric_field 电场分布
     * @param current_density 电流密度分布
     * @param dt 时间步长 (s)
     * @return 求解结果
     * @details 综合求解电子和离子能量平衡方程
     */
    VoidResult timeStep(const std::vector<double>& electron_density,
                        const std::vector<double>& ion_density,
                        const std::vector<Utils::Vector3D>& electric_field,
                        const std::vector<Utils::Vector3D>& current_density, double dt);

    /**
     * @brief 获取能量状态
     * @param position 位置坐标
     * @return 能量状态
     */
    EnergyState getEnergyState(const Utils::Point3D& position) const;

    /**
     * @brief 获取所有网格点的能量状态
     * @return 能量状态向量
     */
    const std::vector<EnergyState>& getAllEnergyStates() const;

    /**
     * @brief 获取电子温度分布
     * @return 电子温度分布向量
     */
    std::vector<double> getElectronTemperatureDistribution() const;

    /**
     * @brief 获取离子温度分布
     * @return 离子温度分布向量
     */
    std::vector<double> getIonTemperatureDistribution() const;

    /**
     * @brief 获取能量密度分布
     * @return 能量密度分布向量 (电子+离子)
     */
    std::vector<double> getEnergyDensityDistribution() const;

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
    EnergyBalanceParameters params_; ///< 能量平衡参数

    // 网格信息
    Utils::Vector3D grid_size_;  ///< 网格尺寸
    Utils::Vector3D resolution_; ///< 网格分辨率
    Utils::Vector3D spacing_;    ///< 网格间距
    size_t total_points_;        ///< 总网格点数

    // 能量状态数据
    std::vector<EnergyState> energy_states_;     ///< 能量状态
    std::vector<EnergyState> energy_states_old_; ///< 上一时刻的能量状态

    // 统计信息
    EnergyBalanceStatistics statistics_; ///< 统计信息

    // 私有方法 - 网格操作
    size_t getGridIndex(size_t i, size_t j, size_t k) const;
    Utils::Point3D getGridPosition(size_t i, size_t j, size_t k) const;
    void getGridIndices(size_t index, size_t& i, size_t& j, size_t& k) const;

    // 私有方法 - 数值方法
    double calculateLaplacian(size_t index, const std::vector<double>& field) const;
    Utils::Vector3D calculateGradient(size_t index, const std::vector<double>& field) const;
    double calculateDivergence(size_t index, const std::vector<Utils::Vector3D>& field) const;

    // 私有方法 - 物理计算
    double calculateJouleHeating(double current_density_magnitude, double conductivity);
    double calculateElasticCooling(double electron_temp, double ion_temp, double density);
    double calculateInelasticCooling(double electron_temp, double density);
    double calculateRadiationLoss(double electron_temp, double ion_temp, double density);
    double calculateEnergyExchange(double electron_temp, double ion_temp, double density);

    // 私有方法 - 边界条件
    void applyBoundaryConditions();
    bool isBoundaryPoint(size_t i, size_t j, size_t k) const;

    // 私有方法 - 辅助函数
    void updateStatistics();
    double interpolateField(const Utils::Point3D& position, const std::vector<double>& field) const;
    void checkEnergyConservation();
};

} // namespace FieldSolver
} // namespace SCDAT
