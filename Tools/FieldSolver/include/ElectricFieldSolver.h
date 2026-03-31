/**
 * @file EnhancedElectricFieldSolver.h
 * @brief 增强版电场求解器头文件
 * @details 基于Arc-PIC实现的多边界条件电场求解器
 *
 * @author Wang Sizhan
 * @date 2026.3.26 7:46:01
 * @version V0.0.1
 *
 * @ingroup FieldSolverModule
 */

#ifndef SCDAT_FIELD_ENHANCED_ELECTRIC_FIELD_SOLVER_H
#define SCDAT_FIELD_ENHANCED_ELECTRIC_FIELD_SOLVER_H

#include "../../Boundary/include/BoundaryType.h"
#include "../../Basic/include/Constants.h"
#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
#include "VolField.h"
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{

namespace Field
{

/**
 * @brief 求解器类型枚举
 */
enum class SolverType
{
    FINITE_DIFFERENCE = 0, ///< 有限差分方法
    FINITE_ELEMENT = 1,    ///< 有限元方法
    SPECTRAL = 2,          ///< 谱方法 / Spectral method
    MULTIGRID = 3          ///< 多重网格方法 / Multigrid method
};

/**
 * @brief 边界条件配置结构体
 */
struct BoundaryConditionConfig
{   
    SCDAT::BoundaryType type = SCDAT::BoundaryType::STANDARD;              ///< 边界条件类型
    double value = 0.0;                                                    ///< 边界值
    SCDAT::Geometry::Vector3D normal = SCDAT::Geometry::Vector3D(1, 0, 0); ///< 边界法向量
    std::function<double(const SCDAT::Geometry::Point3D&)> function;       ///< 边界函数
    bool is_time_dependent = false;                                        ///< 是否时间相关
    double time_coefficient = 1.0;                                         ///< 时间系数
};

/**
 * @brief 求解器配置结构体
 */
struct SolverConfig
{
    SolverType solver_type = SolverType::FINITE_DIFFERENCE; ///< 求解器类型
    double tolerance = 1e-12;                               ///< 收敛容差
    int max_iterations = 10000;                             ///< 最大迭代次数
    double relaxation_factor = 1.0;                         ///< 松弛因子
    bool enable_superlu = false;                            ///< 启用SuperLU求解器
    bool enable_openmp = false;                             ///< 启用OpenMP并行
    bool enable_preconditioner = true;                      ///< 启用预条件器
    std::string preconditioner_type = "jacobi";             ///< 预条件器类型
};

/**
 * @brief 网格配置结构体
 */
struct GridConfig
{
    SCDAT::Geometry::Point3D domain_min = SCDAT::Geometry::Point3D(-1, -1, -1); ///< 计算域最小值
    SCDAT::Geometry::Point3D domain_max = SCDAT::Geometry::Point3D(1, 1, 1);    ///< 计算域最大值
    std::vector<int> grid_size = {64, 64, 64};                                  ///< 网格尺寸，支持1D/2D：{Nx,1,1} 或 {Nx,Ny,1}
    bool is_uniform = true;                                                     ///< 是否均匀网格
    double grid_spacing_x = 0.03125;                                            ///< X方向网格间距
    double grid_spacing_y = 0.03125;                                            ///< Y方向网格间距
    double grid_spacing_z = 0.03125;                                            ///< Z方向网格间距
};

/**
 * @brief 电场求解器类
 * @details 高性能电场求解器
 */
class ElectricFieldSolver
{
  public:
    /**
     * @brief 构造函数
     * @param grid_config 网格配置
     * @param solver_config 求解器配置
     */
    ElectricFieldSolver(const GridConfig& grid_config, const SolverConfig& solver_config,
                           double eps0 = SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity);

    /**
     * @brief 析构函数
     */
    ~ElectricFieldSolver();

    /**
     * @brief 初始化求解器
     * @return 是否成功
     */
    bool initialize();

    /**
     * @brief 设置边界条件
     * @param boundary_id 边界ID
     * @param config 边界条件配置
     */
    void setBoundaryCondition(const std::string& boundary_id,
                              const BoundaryConditionConfig& config);

    /**
     * @brief 设置电荷密度分布
     * @param charge_density 电荷密度场
     */
    void setChargeDensity(const std::vector<double>& charge_density);

    /**
     * @brief 求解泊松方程
     * @details 求解 ∇²φ = -ρ/ε₀
     * @return 是否收敛
     */
    bool solvePoissonEquation();

    /**
     * @brief 计算电场
     * @details 计算 E = -∇φ
     */
    void calculateElectricField();

    /**
     * @brief 获取电势场
     * @return 电势场数据
     */
    const std::vector<double>& getPotential() const
    {
        return potential_;
    }

    /**
     * @brief 获取电场
     * @return 电场数据
     */
    const std::vector<SCDAT::Geometry::Vector3D>& getElectricField() const
    {
        return electric_field_;
    }

    /**
     * @brief 获取求解器统计信息
     */
    struct SolverStatistics
    {
        int iterations_used = 0;      ///< 使用的迭代次数
        double final_residual = 0.0;  ///< 最终残差
        double solve_time = 0.0;      ///< 求解时间
        bool converged = false;       ///< 是否收敛
        std::string convergence_info; ///< 收敛信息
    };

    /**
     * @brief 获取统计信息
     */
    const SolverStatistics& getStatistics() const
    {
        return statistics_;
    }

    /**
     * @brief 重置求解器
     */
    void reset();

    /**
     * @brief 更新时间相关边界条件
     * @param current_time 当前时间
     */
    void updateTimeDependentBoundaries(double current_time);

  private:
    GridConfig grid_config_;                                    ///< 网格配置
    SolverConfig solver_config_;                                ///< 求解器配置
    std::map<std::string, BoundaryConditionConfig> boundaries_; ///< 边界条件映射

    std::vector<double> potential_;                         ///< 电势场
    std::vector<SCDAT::Geometry::Vector3D> electric_field_; ///< 电场
    std::vector<double> charge_density_;                    ///< 电荷密度

    SolverStatistics statistics_; ///< 统计信息
    bool is_initialized_;         ///< 初始化标志

    double eps0_;                ///< 真空介电常数

    bool solveWithFiniteDifference();     ///< 有限差分求解
    bool solveWithSuperLU();              ///< SuperLU求解
    bool solveWithIterativeMethod();      ///< 迭代方法求解
    bool solveWithFiniteElement();        ///< 有限元方法求解
    bool solveWithIterativeFEM();         ///< 迭代有限元求解
    bool solveWithDirectFEM();            ///< 直接有限元求解
    bool solveWithPreconditionedMethod(); ///< 预条件方法求解

    // 预条件器相关方法
    void buildPreconditioner(std::vector<double>& preconditioner);       ///< 构建预条件器
    void buildJacobiPreconditioner(std::vector<double>& preconditioner); ///< 构建Jacobi预条件器
    void buildSSORPreconditioner(std::vector<double>& preconditioner);   ///< 构建SSOR预条件器
    void applyPreconditioner(const std::vector<double>& input, std::vector<double>& output,
                             const std::vector<double>& preconditioner); ///< 应用预条件器
    void calculateResidualVector(std::vector<double>& residual);         ///< 计算残差向量
    void calculateMatrixVectorProduct(const std::vector<double>& input,
                                      std::vector<double>& output); ///< 计算矩阵向量乘积

    // 边界条件应用方法
    void applyAllBoundaryConditions(); ///< 应用所有边界条件

    // 边界条件处理方法
    void applyStandardBoundary(const std::string& boundary_id);  ///< 应用标准边界条件
    void applyModifiedBoundary(const std::string& boundary_id);  ///< 应用修改边界条件
    void applyRadialBoundary(const std::string& boundary_id);    ///< 应用径向边界条件
    void applyPeriodicBoundary(const std::string& boundary_id);  ///< 应用周期边界条件
    void applyDirichletBoundary(const std::string& boundary_id); ///< 应用狄利克雷边界
    void applyNeumannBoundary(const std::string& boundary_id);   ///< 应用诺伊曼边界
    void applyRobinBoundary(const std::string& boundary_id);     ///< 应用罗宾边界

    // 辅助方法
    size_t getLinearIndex(int i, int j, int k) const;                ///< 获取线性索引
    std::tuple<int, int, int> get3DIndex(size_t linear_index) const; ///< 获取3D索引
    bool isValidIndex(int i, int j, int k) const;                    ///< 检查索引有效性
    double calculateResidual() const;                                ///< 计算残差
    void updateStatistics(int iterations, double residual, double solve_time,
                          bool converged); ///< 更新统计
};

/**
 * @brief 电场求解器工厂类
 */
class ElectricFieldSolverFactory
{
  public:
    /**
     * @brief 创建求解器
     */
    static std::unique_ptr<ElectricFieldSolver>
    createArcPICSolver(const GridConfig& grid_config, const SolverConfig& solver_config);

    /**
     * @brief 创建标准求解器 / Create standard solver
     */
    static std::unique_ptr<ElectricFieldSolver>
    createStandardSolver(const SCDAT::Geometry::Point3D& domain_min,
                         const SCDAT::Geometry::Point3D& domain_max,
                         const std::vector<int>& grid_size);

    /**
     * @brief 从配置文件创建求解器
     */
    static std::unique_ptr<ElectricFieldSolver> createFromConfig(const std::string& config_file);
};

} // namespace Field
} // namespace SCDAT

#endif // SCDAT_FIELD_ELECTRIC_FIELD_SOLVER_H
