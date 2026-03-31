#pragma once

/**
 * @file ParallelSolver.h
 * @ingroup FieldSolverModule
 */

#include <vector>
#include <memory>
#include <functional>
#include <omp.h>
#include "../../Solver/include/Preconditioner.h"
#include "../../Solver/include/SolverResult.h"
#include "../../Solver/include/SparseMatrix.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Geometry/include/Matrix3x3.h"

namespace SCDAT {
namespace Solver {

/**
 * @brief 并行求解器基类
 * 
 * 提供OpenMP并行化的求解器基础功能，包括负载均衡、内存访问优化等
 */
class ParallelSolver {
public:
    /**
     * @brief 构造函数
     * @param num_threads 线程数量，0表示使用系统默认
     */
    explicit ParallelSolver(int num_threads = 0);
    
    /**
     * @brief 析构函数
     */
    virtual ~ParallelSolver() = default;

    /**
     * @brief 设置线程数量
     * @param num_threads 线程数量
     */
    void setNumThreads(int num_threads);
    
    /**
     * @brief 获取当前线程数量
     * @return 线程数量
     */
    int getNumThreads() const { return num_threads_; }
    
    /**
     * @brief 获取并行效率
     * @return 并行效率 (0-1)
     */
    double getParallelEfficiency() const { return parallel_efficiency_; }

protected:
    /**
     * @brief 并行矩阵向量乘法
     * @param matrix 稀疏矩阵
     * @param vector 输入向量
     * @param result 结果向量
     */
    void parallelMatrixVectorMultiply(
        const SparseMatrix& matrix,
        const std::vector<double>& vector,
        std::vector<double>& result
    );
    
    /**
     * @brief 并行向量运算
     * @param a 向量a
     * @param b 向量b
     * @param result 结果向量
     * @param operation 运算类型
     */
    void parallelVectorOperation(
        const std::vector<double>& a,
        const std::vector<double>& b,
        std::vector<double>& result,
        VectorOperation operation
    );
    
    /**
     * @brief 并行点积计算
     * @param a 向量a
     * @param b 向量b
     * @return 点积结果
     */
    double parallelDotProduct(
        const std::vector<double>& a,
        const std::vector<double>& b
    );
    
    /**
     * @brief 并行向量范数计算
     * @param vector 输入向量
     * @return 向量范数
     */
    double parallelVectorNorm(const std::vector<double>& vector);
    
    /**
     * @brief 负载均衡分配
     * @param total_size 总大小
     * @param thread_id 线程ID
     * @param num_threads 线程数量
     * @return 分配的起始和结束索引
     */
    std::pair<size_t, size_t> getLoadBalancedRange(
        size_t total_size,
        int thread_id,
        int num_threads
    );
    
    /**
     * @brief 更新并行效率统计
     * @param serial_time 串行时间
     * @param parallel_time 并行时间
     */
    void updateParallelEfficiency(double serial_time, double parallel_time);

private:
    int num_threads_;                    ///< 线程数量
    double parallel_efficiency_;         ///< 并行效率
    mutable std::vector<double> timing_data_;  ///< 性能统计数据
};

/**
 * @brief 并行泊松求解器
 * 
 * 基于共轭梯度法的并行泊松方程求解器
 */
class ParallelPoissonSolver : public ParallelSolver {
public:
    /**
     * @brief 构造函数
     * @param num_threads 线程数量
     */
    explicit ParallelPoissonSolver(int num_threads = 0);
    
    /**
     * @brief 求解泊松方程
     * @param matrix 系数矩阵
     * @param rhs 右端项
     * @param solution 解向量
     * @param tolerance 收敛容差
     * @param max_iterations 最大迭代次数
     * @return 求解结果
     */
    SolverResult solve(
        const SparseMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        double tolerance = 1e-10,
        int max_iterations = 1000
    );
    
    /**
     * @brief 设置预处理器
     * @param preconditioner 预处理器
     */
    void setPreconditioner(std::shared_ptr<Preconditioner> preconditioner);

private:
    /**
     * @brief 并行共轭梯度求解
     */
    SolverResult parallelConjugateGradient(
        const SparseMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        double tolerance,
        int max_iterations
    );
    
    /**
     * @brief 应用预处理器
     */
    void applyPreconditioner(
        const std::vector<double>& input,
        std::vector<double>& output
    );

private:
    std::shared_ptr<Preconditioner> preconditioner_;  ///< 预处理器
    std::vector<double> temp_vectors_[4];             ///< 临时向量
};

/**
 * @brief 并行线性系统求解器
 * 
 * 支持多种并行求解算法的通用线性系统求解器
 */
class ParallelLinearSolver : public ParallelSolver {
public:
    /**
     * @brief 求解器类型
     */
    enum class SolverType {
        CONJUGATE_GRADIENT,    ///< 共轭梯度法
        BICGSTAB,             ///< BiCGSTAB方法
        GMRES                 ///< GMRES方法
    };
    
    /**
     * @brief 构造函数
     * @param solver_type 求解器类型
     * @param num_threads 线程数量
     */
    explicit ParallelLinearSolver(
        SolverType solver_type = SolverType::CONJUGATE_GRADIENT,
        int num_threads = 0
    );
    
    /**
     * @brief 求解线性系统
     * @param matrix 系数矩阵
     * @param rhs 右端项
     * @param solution 解向量
     * @param tolerance 收敛容差
     * @param max_iterations 最大迭代次数
     * @return 求解结果
     */
    SolverResult solve(
        const SparseMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        double tolerance = 1e-10,
        int max_iterations = 1000
    );
    
    /**
     * @brief 设置求解器类型
     * @param solver_type 求解器类型
     */
    void setSolverType(SolverType solver_type) { solver_type_ = solver_type; }
    
    /**
     * @brief 获取求解器类型
     * @return 求解器类型
     */
    SolverType getSolverType() const { return solver_type_; }

private:
    /**
     * @brief BiCGSTAB求解器
     */
    SolverResult parallelBiCGSTAB(
        const SparseMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        double tolerance,
        int max_iterations
    );
    
    /**
     * @brief GMRES求解器
     */
    SolverResult parallelGMRES(
        const SparseMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        double tolerance,
        int max_iterations
    );

private:
    SolverType solver_type_;              ///< 求解器类型
    std::vector<std::vector<double>> workspace_;  ///< 工作空间
};

/**
 * @brief 性能监控器
 *
 * 监控并行求解器的性能指标
 */
class ParallelPerformanceMonitor {
public:
    /**
     * @brief 性能指标结构
     */
    struct PerformanceMetrics {
        double serial_time;           ///< 串行时间
        double parallel_time;         ///< 并行时间
        double speedup;              ///< 加速比
        double efficiency;           ///< 并行效率
        double memory_usage;         ///< 内存使用量
        int num_threads;             ///< 线程数量
    };

    /**
     * @brief 开始性能监控
     */
    void startMonitoring();

    /**
     * @brief 结束性能监控
     */
    void stopMonitoring();

    /**
     * @brief 获取性能指标
     * @return 性能指标
     */
    PerformanceMetrics getMetrics() const;

    /**
     * @brief 重置统计数据
     */
    void reset();

private:
    /**
     * @brief 获取当前内存使用情况
     * @return 内存使用量（MB）
     */
    double getCurrentMemoryUsage() const;

private:
    double start_time_;
    double end_time_;
    PerformanceMetrics metrics_;
};

} // namespace Solver
} // namespace SCDAT
