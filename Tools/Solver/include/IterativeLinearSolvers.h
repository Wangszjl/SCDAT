/**
 * @file PoissonSolver.h
 * @brief 泊松方程求解器模块头文件 / Poisson Equation Solver Module Header
 * @details 定义三维非结构化网格上的泊松方程求解器，包括带状稀疏矩阵、
 *          ILU预处理器和预处理共轭梯度法求解器
 *
 * 主要功能 / Main Features:
 * - 带状稀疏矩阵存储和运算 / Band sparse matrix storage and operations
 * - 不完全LU分解预处理 / Incomplete LU decomposition preconditioning
 * - 预处理共轭梯度法求解 / Preconditioned conjugate gradient solving
 * - 多边界条件支持 / Multi-boundary condition support
 * - 高精度数值求解 / High-precision numerical solving
 *
 * 算法基础 / Algorithm Foundation:
 * 基于SPIS 6.1.2原版ConjGrad3DUnstructPoissonSolver的C++17现代化实现
 *
 * @author PIC-SPIS开发团队 / PIC-SPIS Development Team
 * @date 2025年1月27日 16:35:00
 * @version V2.0.0-Enhanced
 * @copyright 开源许可证 / Open Source License
 * @ingroup FieldSolverModule
 */

#ifndef SCDAT_ITERATIVE_LINEAR_SOLVERS_H
#define SCDAT_ITERATIVE_LINEAR_SOLVERS_H

#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Mesh/include/MeshAlgorithms.h"
#include "../../Mesh/include/MeshPartitioning.h"
#include "../../Mesh/include/MeshParsing.h"
#include "ConjugateGradientCommon.h"
#include "BandMatrix.h"
#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>
#include <string>

namespace SCDAT {
namespace Solver {

// 前向声明
class ConjugateGradientSolver;
class ILUPreconditioner;

// 引入 Core 中的 BandMatrix
using BandMatrix = SCDAT::Core::BandMatrix;

// 类型别名
using BandMatrixPtr = std::shared_ptr<BandMatrix>;
using ConjugateGradientSolverPtr = std::shared_ptr<ConjugateGradientSolver>;
using ILUPreconditionerPtr = std::shared_ptr<ILUPreconditioner>;

using NodeId = Mesh::NodeId;
using ElementId = Mesh::ElementId;

/**
 * @brief 带状稀疏矩阵类 / Band Sparse Matrix Class
 * @details 基于SPIS原版BandMatrix的C++17现代化实现，使用压缩行存储格式(CSR)
 *          优化了带状矩阵的存储和运算效率，专门用于有限元泊松方程求解
 *
 * 技术特点 / Technical Features:
 * - CSR压缩存储格式，内存效率高 / CSR compressed storage format, high memory efficiency
 * - 带宽优化算法，减少填充元素 / Bandwidth optimization algorithm, reduce fill-in elements
 * - Cuthill-McKee重排序支持 / Cuthill-McKee reordering support
 * - 高效矩阵向量乘法 / Efficient matrix-vector multiplication
 * - 对角元素快速访问 / Fast diagonal element access
 *
 * 数学基础 / Mathematical Foundation:
 * 用于求解线性系统 Ax = b，其中A为带状稀疏矩阵
 *
 * 性能指标 / Performance Metrics:
 * - 内存使用：O(nnz) 其中nnz为非零元素数量
 * - 矩阵向量乘法：O(nnz) 时间复杂度
 * - 带宽：通常 < 0.1% 的矩阵大小
 *
 * @note 矩阵索引从0开始，使用前需要调用finalize()完成组装
 * @warning 大型矩阵操作可能消耗大量内存，建议监控内存使用
 */
using BandMatrix = SCDAT::Core::BandMatrix;

/**
 * @brief 不完全LU分解预处理器
 * 
 * 基于SPIS原版预处理器的C++实现
 * 提供ILU(0)和ILU(k)分解
 */
class ILUPreconditioner {
public:
    // 构造函数
    ILUPreconditioner(int fill_level = 0, double drop_tolerance = 1e-6);
    
    // 析构函数
    ~ILUPreconditioner() = default;
    
    // 预处理器设置
    void setFillLevel(int level) { fill_level_ = level; }
    void setDropTolerance(double tol) { drop_tolerance_ = tol; }
    
    // 分解和求解
    bool factorize(const BandMatrix& matrix);
    void solve(const std::vector<double>& rhs, std::vector<double>& solution) const;
    void solveForward(const std::vector<double>& rhs, std::vector<double>& y) const;
    void solveBackward(const std::vector<double>& y, std::vector<double>& solution) const;
    
    // 状态查询
    bool isFactorized() const { return factorized_; }
    int getMemoryUsage() const;
    
    // 调试
    void printFactors() const;

private:
    int fill_level_;                    ///< 填充级别
    double drop_tolerance_;             ///< 丢弃容差
    bool factorized_;                   ///< 是否已分解
    
    BandMatrixPtr L_;                   ///< 下三角矩阵
    BandMatrixPtr U_;                   ///< 上三角矩阵
    std::vector<int> pivot_;            ///< 主元置换
    
    // 分解算法
    void iluDecomposition(const BandMatrix& matrix);
    void buildLUMatrices(const std::vector<std::vector<std::pair<int, double>>>& work_matrix, int n);
    void dropSmallEntries();
};

/**
 * @brief 共轭梯度求解器
 * 
 * 基于SPIS原版ConjGrad3DUnstructPoissonSolver的C++实现
 * 支持预处理共轭梯度法求解线性系统
 */
class ConjugateGradientSolver {
public:
    // 求解器参数结构
    struct Parameters {
        double tolerance = 1e-8;        ///< 收敛容差
        int max_iterations = 1000;      ///< 最大迭代次数
        bool use_preconditioner = true; ///< 是否使用预处理器
        int preconditioner_type = 0;    ///< 预处理器类型 (0=ILU, 1=Diagonal)
        bool verbose = false;           ///< 是否输出详细信息
        int print_frequency = 10;       ///< 打印频率
    };
    
    // 求解器统计信息
    struct Statistics {
        int iterations = 0;             ///< 实际迭代次数
        double final_residual = 0.0;    ///< 最终残差
        double convergence_rate = 0.0;   ///< 收敛率
        double solve_time = 0.0;         ///< 求解时间
        bool converged = false;          ///< 是否收敛
    };
    
    // 构造函数
    ConjugateGradientSolver();
    ConjugateGradientSolver(const Parameters& params);
    
    // 析构函数
    ~ConjugateGradientSolver() = default;
    
    // 参数设置
    void setParameters(const Parameters& params) { params_ = params; }
    const Parameters& getParameters() const { return params_; }
    
    // 预处理器设置
    void setPreconditioner(ILUPreconditionerPtr preconditioner);
    void setDiagonalPreconditioner(const std::vector<double>& diagonal);
    
    // 求解方法
    bool solve(const BandMatrix& matrix, 
               const std::vector<double>& rhs,
               std::vector<double>& solution);
    
    bool solve(const BandMatrix& matrix,
               const std::vector<double>& rhs,
               std::vector<double>& solution,
               const std::vector<double>& initial_guess);
    
    // 统计信息
    const Statistics& getStatistics() const { return stats_; }
    void printStatistics() const;
    
    // 残差计算
    double computeResidual(const BandMatrix& matrix,
                          const std::vector<double>& solution,
                          const std::vector<double>& rhs) const;

private:
    Parameters params_;                 ///< 求解器参数
    Statistics stats_;                  ///< 统计信息
    ILUPreconditionerPtr preconditioner_; ///< ILU预处理器
    std::vector<double> diagonal_preconditioner_; ///< 对角预处理器
    
    // 核心算法
    bool conjugateGradientIteration(const BandMatrix& matrix,
                                   const std::vector<double>& rhs,
                                   std::vector<double>& solution);
    
    // 预处理器应用
    void applyPreconditioner(const std::vector<double>& r, 
                           std::vector<double>& z) const;
    
    // 辅助方法
    void printIteration(int iter, double residual) const;

    static constexpr double kBreakdownTolerance = 1e-30;
};

} // namespace Solver
} // namespace SCDAT

#endif // SCDAT_ITERATIVE_LINEAR_SOLVERS_H
