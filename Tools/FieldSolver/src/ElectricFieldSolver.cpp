#include "ElectricFieldSolver.h"

#include "../Basic/include/Constants.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace SCDAT
{
namespace Field
{

namespace Constants = SCDAT::Basic::Constants;

/**
 * @brief 电场求解器构造函数。
 * @details
 * 根据网格与求解器配置初始化内部存储，包括电势场、电场和电荷密度数组。
 * 若输入网格维度配置非法（非3个分量），将回退为默认网格尺寸 {64,64,64}。
 *
 * @param grid_config 网格配置。
 * @param solver_config 求解器配置。
 * @param eps0 真空介电常数，默认值见头文件。
 */
ElectricFieldSolver::ElectricFieldSolver(const GridConfig& grid_config,
                                         const SolverConfig& solver_config, double eps0)
    : grid_config_(grid_config), solver_config_(solver_config), is_initialized_(false), eps0_(eps0)
{
    if (grid_config_.grid_size.size() != 3)
    {
        grid_config_.grid_size = {64, 64, 64};
    }

    const size_t total_points = static_cast<size_t>(std::max(1, grid_config_.grid_size[0])) *
                                static_cast<size_t>(std::max(1, grid_config_.grid_size[1])) *
                                static_cast<size_t>(std::max(1, grid_config_.grid_size[2]));

    potential_.assign(total_points, 0.0);
    electric_field_.assign(total_points, SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0));
    charge_density_.assign(total_points, 0.0);
}

/**
 * @brief 电场求解器析构函数。
 * @details 当前实现使用默认析构，无额外资源释放逻辑。
 */
ElectricFieldSolver::~ElectricFieldSolver() = default;

/**
 * @brief 初始化求解器。
 * @details
 * 执行网格有效性检查、网格步长计算（均匀网格）以及默认边界条件构建。
 * 支持降维场景：1D={Nx,1,1}，2D={Nx,Ny,1}。
 *
 * @return 初始化成功返回 true，否则返回 false。
 */
bool ElectricFieldSolver::initialize()
{
    if (is_initialized_)
    {
        return true;
    }

    if (grid_config_.grid_size.size() != 3)
    {
        return false;
    }

    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];
    if (nx < 1 || ny < 1 || nz < 1)
    {
        return false;
    }

    // 支持降维求解：1D={Nx,1,1}，2D={Nx,Ny,1}
    if (nx <= 1 && ny <= 1 && nz <= 1)
    {
        return false;
    }

    if (grid_config_.is_uniform)
    {
        grid_config_.grid_spacing_x =
            (nx > 1) ? (grid_config_.domain_max.x() - grid_config_.domain_min.x()) /
                           static_cast<double>(nx - 1)
                     : 1.0;
        grid_config_.grid_spacing_y =
            (ny > 1) ? (grid_config_.domain_max.y() - grid_config_.domain_min.y()) /
                           static_cast<double>(ny - 1)
                     : 1.0;
        grid_config_.grid_spacing_z =
            (nz > 1) ? (grid_config_.domain_max.z() - grid_config_.domain_min.z()) /
                           static_cast<double>(nz - 1)
                     : 1.0;
    }

    if (boundaries_.empty())
    {
        BoundaryConditionConfig d;
        d.type = SCDAT::BoundaryType::STANDARD;
        d.value = 0.0;

        setBoundaryCondition("x_min", d);
        setBoundaryCondition("x_max", d);
        setBoundaryCondition("y_min", d);
        setBoundaryCondition("y_max", d);
        setBoundaryCondition("z_min", d);
        setBoundaryCondition("z_max", d);
    }

    is_initialized_ = true;
    return true;
}

/**
 * @brief 设置指定边界的边界条件配置。
 * @param boundary_id 边界标识，例如 x_min/x_max/y_min/y_max/z_min/z_max。
 * @param config 边界条件配置对象。
 */
void ElectricFieldSolver::setBoundaryCondition(const std::string& boundary_id,
                                               const BoundaryConditionConfig& config)
{
    boundaries_[boundary_id] = config;
}

/**
 * @brief 设置电荷密度场。
 * @details 输入数组长度必须与网格总节点数一致。
 *
 * @param charge_density 电荷密度数组。
 * @throws std::invalid_argument 当输入长度与内部场长度不一致时抛出。
 */
void ElectricFieldSolver::setChargeDensity(const std::vector<double>& charge_density)
{
    if (charge_density.size() != charge_density_.size())
    {
        throw std::invalid_argument("charge density size mismatch");
    }
    charge_density_ = charge_density;
}

/**
 * @brief 求解泊松方程。
 * @details
 * 根据配置选择数值方法求解，并统计求解耗时、残差和收敛信息。
 * 方程形式为 @f$\nabla^2\phi=-\rho/\varepsilon_0@f$。
 *
 * @return 求解成功并达到对应流程返回条件时为 true，否则为 false。
 */
bool ElectricFieldSolver::solvePoissonEquation()
{
    if (!is_initialized_ && !initialize())
    {
        return false;
    }

    const auto t0 = std::chrono::high_resolution_clock::now();

    bool ok = false;
    switch (solver_config_.solver_type)
    {
    case SolverType::FINITE_DIFFERENCE:
        ok = solver_config_.enable_superlu ? solveWithSuperLU() : solveWithIterativeMethod();
        break;
    case SolverType::FINITE_ELEMENT:
        ok = solveWithFiniteElement();
        break;
    case SolverType::SPECTRAL:
        ok = solveWithFiniteDifference();
        break;
    case SolverType::MULTIGRID:
        ok = solveWithPreconditionedMethod();
        break;
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    const double dt = static_cast<double>(
                          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
                      1e6;

    updateStatistics(statistics_.iterations_used, calculateResidual(), dt, ok);
    return ok;
}

/**
 * @brief 由电势场计算电场。
 * @details
 * 基于中心差分近似空间梯度，按照关系式 @f$\mathbf{E}=-\nabla\phi@f$ 计算每个网格点电场。
 * 边界点采用索引夹取方式避免越界访问。
 */
void ElectricFieldSolver::calculateElectricField()
{
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];

    const double dx = std::max(1e-30, grid_config_.grid_spacing_x);
    const double dy = std::max(1e-30, grid_config_.grid_spacing_y);
    const double dz = std::max(1e-30, grid_config_.grid_spacing_z);

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                const size_t idx = getLinearIndex(i, j, k);

                auto phi = [&](int ii, int jj, int kk) -> double
                {
                    ii = std::max(0, std::min(ii, nx - 1));
                    jj = std::max(0, std::min(jj, ny - 1));
                    kk = std::max(0, std::min(kk, nz - 1));
                    return potential_[getLinearIndex(ii, jj, kk)];
                };

                const double ex = -(phi(i + 1, j, k) - phi(i - 1, j, k)) / (2.0 * dx);
                const double ey = -(phi(i, j + 1, k) - phi(i, j - 1, k)) / (2.0 * dy);
                const double ez = -(phi(i, j, k + 1) - phi(i, j, k - 1)) / (2.0 * dz);

                electric_field_[idx] = SCDAT::Geometry::Vector3D(ex, ey, ez);
            }
        }
    }
}

/**
 * @brief 重置求解器状态。
 * @details 将电势、电场、电荷密度和统计信息恢复到初始状态。
 */
void ElectricFieldSolver::reset()
{
    std::fill(potential_.begin(), potential_.end(), 0.0);
    std::fill(electric_field_.begin(), electric_field_.end(),
              SCDAT::Geometry::Vector3D(0.0, 0.0, 0.0));
    std::fill(charge_density_.begin(), charge_density_.end(), 0.0);
    statistics_ = SolverStatistics{};
}

/**
 * @brief 更新时间相关边界条件。
 * @details
 * 对启用时间依赖的边界进行更新：
 * - 若提供函数，则使用函数值乘以时间系数与当前时间；
 * - 否则直接将边界值按时间系数递推更新。
 *
 * @param current_time 当前物理时间。
 */
void ElectricFieldSolver::updateTimeDependentBoundaries(double current_time)
{
    for (auto& kv : boundaries_)
    {
        if (!kv.second.is_time_dependent)
        {
            continue;
        }

        if (kv.second.function)
        {
            kv.second.value = kv.second.function(SCDAT::Geometry::Point3D()) *
                              kv.second.time_coefficient * current_time;
        }
        else
        {
            kv.second.value *= kv.second.time_coefficient;
        }
    }
}

/**
 * @brief 有限差分求解入口。
 * @details 当前实现转发到迭代求解器。
 * @return 求解是否达到收敛判据。
 */
bool ElectricFieldSolver::solveWithFiniteDifference()
{
    return solveWithIterativeMethod();
}

/**
 * @brief SuperLU 直接法求解入口。
 * @details 当前版本在无直接后端时回退到迭代求解。
 * @return 求解是否达到收敛判据。
 */
bool ElectricFieldSolver::solveWithSuperLU()
{
    // Fallback to iterative solve when direct backend is unavailable.
    return solveWithIterativeMethod();
}

/**
 * @brief 迭代法求解泊松方程。
 * @details
 * 使用基于离散拉普拉斯算子的迭代更新，并结合松弛因子执行加权更新。
 * 每轮迭代后应用边界条件，并使用最大改变量判断收敛。
 *
 * @return 收敛返回 true，否则返回 false。
 */
bool ElectricFieldSolver::solveWithIterativeMethod()
{
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];
    const double dx = std::max(1e-30, grid_config_.grid_spacing_x);
    const double dy = std::max(1e-30, grid_config_.grid_spacing_y);
    const double dz = std::max(1e-30, grid_config_.grid_spacing_z);
    const double eps0 = SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity;

    std::vector<double> new_phi = potential_;

    const bool use_x = nx > 1;
    const bool use_y = ny > 1;
    const bool use_z = nz > 1;
    const double idx2 = use_x ? (1.0 / (dx * dx)) : 0.0;
    const double idy2 = use_y ? (1.0 / (dy * dy)) : 0.0;
    const double idz2 = use_z ? (1.0 / (dz * dz)) : 0.0;
    const double denom = 2.0 * (idx2 + idy2 + idz2);
    if (denom <= 0.0)
    {
        return false;
    }

    const int i_begin = use_x ? 1 : 0;
    const int i_end = use_x ? (nx - 1) : 1;
    const int j_begin = use_y ? 1 : 0;
    const int j_end = use_y ? (ny - 1) : 1;
    const int k_begin = use_z ? 1 : 0;
    const int k_end = use_z ? (nz - 1) : 1;

    bool converged = false;
    int iter = 0;

    for (iter = 0; iter < solver_config_.max_iterations; ++iter)
    {
        double max_change = 0.0;

        for (int i = i_begin; i < i_end; ++i)
        {
            for (int j = j_begin; j < j_end; ++j)
            {
                for (int k = k_begin; k < k_end; ++k)
                {
                    const size_t c = getLinearIndex(i, j, k);
                    const size_t xm = use_x ? getLinearIndex(i - 1, j, k) : c;
                    const size_t xp = use_x ? getLinearIndex(i + 1, j, k) : c;
                    const size_t ym = use_y ? getLinearIndex(i, j - 1, k) : c;
                    const size_t yp = use_y ? getLinearIndex(i, j + 1, k) : c;
                    const size_t zm = use_z ? getLinearIndex(i, j, k - 1) : c;
                    const size_t zp = use_z ? getLinearIndex(i, j, k + 1) : c;

                    const double rhs = -charge_density_[c] / eps0;
                    const double updated = ((potential_[xm] + potential_[xp]) * idx2 +
                                            (potential_[ym] + potential_[yp]) * idy2 +
                                            (potential_[zm] + potential_[zp]) * idz2 - rhs) /
                                           denom;

                    const double blended = potential_[c] + solver_config_.relaxation_factor *
                                                               (updated - potential_[c]);
                    new_phi[c] = blended;
                    max_change = std::max(max_change, std::abs(blended - potential_[c]));
                }
            }
        }

        potential_.swap(new_phi);
        applyAllBoundaryConditions();

        if (max_change < solver_config_.tolerance)
        {
            converged = true;
            break;
        }
    }

    statistics_.iterations_used = iter + 1;
    statistics_.final_residual = calculateResidual();
    statistics_.converged = converged;
    return converged;
}

/**
 * @brief 有限元求解入口。
 * @details 根据配置决定使用迭代有限元或直接有限元路径。
 * @return 求解是否成功。
 */
bool ElectricFieldSolver::solveWithFiniteElement()
{
    return solver_config_.enable_preconditioner ? solveWithIterativeFEM() : solveWithDirectFEM();
}

/**
 * @brief 迭代有限元求解。
 * @details 当前实现复用统一迭代求解流程。
 * @return 求解是否达到收敛判据。
 */
bool ElectricFieldSolver::solveWithIterativeFEM()
{
    return solveWithIterativeMethod();
}

/**
 * @brief 直接有限元求解。
 * @details 当前实现复用统一迭代求解流程。
 * @return 求解是否达到收敛判据。
 */
bool ElectricFieldSolver::solveWithDirectFEM()
{
    return solveWithIterativeMethod();
}

/**
 * @brief 预条件校正求解流程。
 * @details
 * 构建预条件器后，对残差向量进行预条件校正并回写到电势场。
 * 执行一次校正后更新边界和统计信息。
 *
 * @return 校正后残差满足容差时返回 true，否则返回 false。
 */
bool ElectricFieldSolver::solveWithPreconditionedMethod()
{
    std::vector<double> preconditioner;
    buildPreconditioner(preconditioner);

    std::vector<double> residual;
    calculateResidualVector(residual);

    std::vector<double> corr;
    applyPreconditioner(residual, corr, preconditioner);

    for (size_t i = 0; i < potential_.size() && i < corr.size(); ++i)
    {
        potential_[i] += corr[i];
    }

    applyAllBoundaryConditions();
    statistics_.iterations_used = 1;
    statistics_.final_residual = calculateResidual();
    statistics_.converged = statistics_.final_residual < solver_config_.tolerance;
    return statistics_.converged;
}

/**
 * @brief 构建预条件器。
 * @details 根据配置选择 Jacobi 或 SSOR 预条件器。
 * @param preconditioner 输出预条件器对角（或等效）向量。
 */
void ElectricFieldSolver::buildPreconditioner(std::vector<double>& preconditioner)
{
    if (solver_config_.preconditioner_type == "ssor")
    {
        buildSSORPreconditioner(preconditioner);
    }
    else
    {
        buildJacobiPreconditioner(preconditioner);
    }
}

/**
 * @brief 构建 Jacobi 预条件器。
 * @details 当前使用单位对角近似。
 * @param preconditioner 输出预条件器向量。
 */
void ElectricFieldSolver::buildJacobiPreconditioner(std::vector<double>& preconditioner)
{
    preconditioner.assign(potential_.size(), 1.0);
}

/**
 * @brief 构建 SSOR 预条件器。
 * @details
 * 采用与松弛因子相关的简化缩放形式，避免出现除零风险时使用下限保护。
 *
 * @param preconditioner 输出预条件器向量。
 */
void ElectricFieldSolver::buildSSORPreconditioner(std::vector<double>& preconditioner)
{
    preconditioner.assign(potential_.size(),
                          1.0 / std::max(1e-30, 2.0 - solver_config_.relaxation_factor));
}

/**
 * @brief 应用预条件器。
 * @details 对输入向量执行逐点缩放，生成预条件后的输出向量。
 * @param input 输入向量。
 * @param output 输出向量。
 * @param preconditioner 预条件器系数向量。
 */
void ElectricFieldSolver::applyPreconditioner(const std::vector<double>& input,
                                              std::vector<double>& output,
                                              const std::vector<double>& preconditioner)
{
    output.assign(input.size(), 0.0);
    for (size_t i = 0; i < input.size(); ++i)
    {
        const double p = (i < preconditioner.size()) ? preconditioner[i] : 1.0;
        output[i] = input[i] * p;
    }
}

/**
 * @brief 计算残差向量。
 * @details 残差定义为 @f$r=b-Ax@f$，其中 @f$b@f$ 由电荷密度给出。
 * @param residual 输出残差向量。
 */
void ElectricFieldSolver::calculateResidualVector(std::vector<double>& residual)
{
    std::vector<double> Ax;
    calculateMatrixVectorProduct(potential_, Ax);

    residual.assign(potential_.size(), 0.0);
    for (size_t i = 0; i < residual.size(); ++i)
    {
        residual[i] = charge_density_[i] - Ax[i];
    }
}

/**
 * @brief 计算离散算子与向量乘积。
 * @details
 * 对输入电势向量应用离散拉普拉斯算子，得到 @f$A\phi@f$。
 * 可自动适配 1D/2D/3D 维度场景。
 *
 * @param input 输入向量（通常为电势场）。
 * @param output 输出向量（算子作用结果）。
 */
void ElectricFieldSolver::calculateMatrixVectorProduct(const std::vector<double>& input,
                                                       std::vector<double>& output)
{
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];
    const double dx = std::max(1e-30, grid_config_.grid_spacing_x);
    const double dy = std::max(1e-30, grid_config_.grid_spacing_y);
    const double dz = std::max(1e-30, grid_config_.grid_spacing_z);

    output.assign(input.size(), 0.0);

    const bool use_x = nx > 1;
    const bool use_y = ny > 1;
    const bool use_z = nz > 1;
    const double idx2 = use_x ? (1.0 / (dx * dx)) : 0.0;
    const double idy2 = use_y ? (1.0 / (dy * dy)) : 0.0;
    const double idz2 = use_z ? (1.0 / (dz * dz)) : 0.0;

    const int i_begin = use_x ? 1 : 0;
    const int i_end = use_x ? (nx - 1) : 1;
    const int j_begin = use_y ? 1 : 0;
    const int j_end = use_y ? (ny - 1) : 1;
    const int k_begin = use_z ? 1 : 0;
    const int k_end = use_z ? (nz - 1) : 1;

    for (int i = i_begin; i < i_end; ++i)
    {
        for (int j = j_begin; j < j_end; ++j)
        {
            for (int k = k_begin; k < k_end; ++k)
            {
                const size_t c = getLinearIndex(i, j, k);
                const size_t xm = use_x ? getLinearIndex(i - 1, j, k) : c;
                const size_t xp = use_x ? getLinearIndex(i + 1, j, k) : c;
                const size_t ym = use_y ? getLinearIndex(i, j - 1, k) : c;
                const size_t yp = use_y ? getLinearIndex(i, j + 1, k) : c;
                const size_t zm = use_z ? getLinearIndex(i, j, k - 1) : c;
                const size_t zp = use_z ? getLinearIndex(i, j, k + 1) : c;

                output[c] = (input[xm] - 2.0 * input[c] + input[xp]) * idx2 +
                            (input[ym] - 2.0 * input[c] + input[yp]) * idy2 +
                            (input[zm] - 2.0 * input[c] + input[zp]) * idz2;
            }
        }
    }
}

/**
 * @brief 应用全部边界条件。
 * @details 按边界类型分发到对应边界处理函数。
 */
void ElectricFieldSolver::applyAllBoundaryConditions()
{
    for (const auto& kv : boundaries_)
    {
        switch (kv.second.type)
        {
        case SCDAT::BoundaryType::STANDARD:
        case SCDAT::BoundaryType::MODIFIED_1:
        case SCDAT::BoundaryType::MODIFIED_4:
            applyStandardBoundary(kv.first);
            break;
        case SCDAT::BoundaryType::RADIAL:
            applyRadialBoundary(kv.first);
            break;
        case SCDAT::BoundaryType::PERIODIC:
            applyPeriodicBoundary(kv.first);
            break;
        case SCDAT::BoundaryType::DIRICHLET:
            applyDirichletBoundary(kv.first);
            break;
        case SCDAT::BoundaryType::NEUMANN:
        case SCDAT::BoundaryType::MIXED:
            applyNeumannBoundary(kv.first);
            break;
        }
    }
}

/**
 * @brief 应用标准边界条件。
 * @details 当前映射为 Dirichlet 处理。
 * @param boundary_id 边界标识。
 */
void ElectricFieldSolver::applyStandardBoundary(const std::string& boundary_id)
{
    applyDirichletBoundary(boundary_id);
}

/**
 * @brief 应用改进边界条件。
 * @details 当前映射为 Dirichlet 处理。
 * @param boundary_id 边界标识。
 */
void ElectricFieldSolver::applyModifiedBoundary(const std::string& boundary_id)
{
    applyDirichletBoundary(boundary_id);
}

/**
 * @brief 应用径向边界条件。
 * @details 当前映射为 Dirichlet 处理。
 * @param boundary_id 边界标识。
 */
void ElectricFieldSolver::applyRadialBoundary(const std::string& boundary_id)
{
    applyDirichletBoundary(boundary_id);
}

/**
 * @brief 应用周期边界条件。
 * @details
 * 当前实现对 x 方向执行周期复制：将 x=0 面的电势复制到 x=nx-1 面。
 *
 * @param boundary_id 边界标识（当前未使用）。
 */
void ElectricFieldSolver::applyPeriodicBoundary(const std::string& boundary_id)
{
    (void)boundary_id;
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];

    for (int j = 0; j < ny; ++j)
    {
        for (int k = 0; k < nz; ++k)
        {
            potential_[getLinearIndex(nx - 1, j, k)] = potential_[getLinearIndex(0, j, k)];
        }
    }
}

/**
 * @brief 应用 Dirichlet 边界条件。
 * @details 按边界标识将对应边界面电势直接设置为常值。
 * @param boundary_id 边界标识。
 */
void ElectricFieldSolver::applyDirichletBoundary(const std::string& boundary_id)
{
    auto it = boundaries_.find(boundary_id);
    if (it == boundaries_.end())
    {
        return;
    }

    const double v = it->second.value;
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];

    if (boundary_id == "x_min")
    {
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(0, j, k)] = v;
    }
    else if (boundary_id == "x_max")
    {
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(nx - 1, j, k)] = v;
    }
    else if (boundary_id == "y_min")
    {
        for (int i = 0; i < nx; ++i)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(i, 0, k)] = v;
    }
    else if (boundary_id == "y_max")
    {
        for (int i = 0; i < nx; ++i)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(i, ny - 1, k)] = v;
    }
    else if (boundary_id == "z_min")
    {
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                potential_[getLinearIndex(i, j, 0)] = v;
    }
    else if (boundary_id == "z_max")
    {
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                potential_[getLinearIndex(i, j, nz - 1)] = v;
    }
}

/**
 * @brief 应用 Neumann 边界条件。
 * @details
 * 使用最简一阶近似：边界值取相邻内部点，从而近似法向导数为 0。
 *
 * @param boundary_id 边界标识。
 */
void ElectricFieldSolver::applyNeumannBoundary(const std::string& boundary_id)
{
    // Minimal Neumann approximation: copy adjacent interior values.
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];

    if (boundary_id == "x_min" && nx > 1)
    {
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(0, j, k)] = potential_[getLinearIndex(1, j, k)];
    }
    else if (boundary_id == "x_max" && nx > 1)
    {
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(nx - 1, j, k)] = potential_[getLinearIndex(nx - 2, j, k)];
    }
    else if (boundary_id == "y_min" && ny > 1)
    {
        for (int i = 0; i < nx; ++i)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(i, 0, k)] = potential_[getLinearIndex(i, 1, k)];
    }
    else if (boundary_id == "y_max" && ny > 1)
    {
        for (int i = 0; i < nx; ++i)
            for (int k = 0; k < nz; ++k)
                potential_[getLinearIndex(i, ny - 1, k)] = potential_[getLinearIndex(i, ny - 2, k)];
    }
    else if (boundary_id == "z_min" && nz > 1)
    {
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                potential_[getLinearIndex(i, j, 0)] = potential_[getLinearIndex(i, j, 1)];
    }
    else if (boundary_id == "z_max" && nz > 1)
    {
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                potential_[getLinearIndex(i, j, nz - 1)] = potential_[getLinearIndex(i, j, nz - 2)];
    }
}

/**
 * @brief 将三维网格索引映射为一维线性索引。
 * @param i x 方向索引。
 * @param j y 方向索引。
 * @param k z 方向索引。
 * @return 线性存储下标。
 */
size_t ElectricFieldSolver::getLinearIndex(int i, int j, int k) const
{
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    return static_cast<size_t>(k * nx * ny + j * nx + i);
}

/**
 * @brief 将线性索引反解为三维网格索引。
 * @param linear_index 线性下标。
 * @return 返回 (i,j,k) 三元组。
 */
std::tuple<int, int, int> ElectricFieldSolver::get3DIndex(size_t linear_index) const
{
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];

    const int k = static_cast<int>(linear_index / static_cast<size_t>(nx * ny));
    const size_t rem = linear_index % static_cast<size_t>(nx * ny);
    const int j = static_cast<int>(rem / static_cast<size_t>(nx));
    const int i = static_cast<int>(rem % static_cast<size_t>(nx));
    return {i, j, k};
}

/**
 * @brief 检查三维索引是否有效。
 * @param i x 方向索引。
 * @param j y 方向索引。
 * @param k z 方向索引。
 * @return 位于网格范围内返回 true，否则返回 false。
 */
bool ElectricFieldSolver::isValidIndex(int i, int j, int k) const
{
    return i >= 0 && j >= 0 && k >= 0 && i < grid_config_.grid_size[0] &&
           j < grid_config_.grid_size[1] && k < grid_config_.grid_size[2];
}

/**
 * @brief 计算当前电势解的残差。
 * @details
 * 在内部节点上计算离散方程残差绝对值，并返回其最大值作为收敛指标。
 *
 * @return 最大残差范数（无穷范数近似）。
 */
double ElectricFieldSolver::calculateResidual() const
{
    const int nx = grid_config_.grid_size[0];
    const int ny = grid_config_.grid_size[1];
    const int nz = grid_config_.grid_size[2];

    const double dx = std::max(1e-30, grid_config_.grid_spacing_x);
    const double dy = std::max(1e-30, grid_config_.grid_spacing_y);
    const double dz = std::max(1e-30, grid_config_.grid_spacing_z);
    const double eps0 = SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity;

    const bool use_x = nx > 1;
    const bool use_y = ny > 1;
    const bool use_z = nz > 1;
    const double idx2 = use_x ? (1.0 / (dx * dx)) : 0.0;
    const double idy2 = use_y ? (1.0 / (dy * dy)) : 0.0;
    const double idz2 = use_z ? (1.0 / (dz * dz)) : 0.0;

    const int i_begin = use_x ? 1 : 0;
    const int i_end = use_x ? (nx - 1) : 1;
    const int j_begin = use_y ? 1 : 0;
    const int j_end = use_y ? (ny - 1) : 1;
    const int k_begin = use_z ? 1 : 0;
    const int k_end = use_z ? (nz - 1) : 1;

    double max_r = 0.0;
    for (int i = i_begin; i < i_end; ++i)
    {
        for (int j = j_begin; j < j_end; ++j)
        {
            for (int k = k_begin; k < k_end; ++k)
            {
                const size_t c = getLinearIndex(i, j, k);
                const size_t xm = use_x ? getLinearIndex(i - 1, j, k) : c;
                const size_t xp = use_x ? getLinearIndex(i + 1, j, k) : c;
                const size_t ym = use_y ? getLinearIndex(i, j - 1, k) : c;
                const size_t yp = use_y ? getLinearIndex(i, j + 1, k) : c;
                const size_t zm = use_z ? getLinearIndex(i, j, k - 1) : c;
                const size_t zp = use_z ? getLinearIndex(i, j, k + 1) : c;

                const double lap = (potential_[xm] - 2.0 * potential_[c] + potential_[xp]) * idx2 +
                                   (potential_[ym] - 2.0 * potential_[c] + potential_[yp]) * idy2 +
                                   (potential_[zm] - 2.0 * potential_[c] + potential_[zp]) * idz2;
                const double r = std::abs(lap + charge_density_[c] / eps0);
                max_r = std::max(max_r, r);
            }
        }
    }
    return max_r;
}

/**
 * @brief 更新求解统计信息。
 * @param iterations 实际迭代次数。
 * @param residual 最终残差。
 * @param solve_time 求解耗时（秒）。
 * @param converged 是否收敛。
 */
void ElectricFieldSolver::updateStatistics(int iterations, double residual, double solve_time,
                                           bool converged)
{
    statistics_.iterations_used = iterations;
    statistics_.final_residual = residual;
    statistics_.solve_time = solve_time;
    statistics_.converged = converged;
    statistics_.convergence_info = converged ? "converged" : "not converged";
}

/**
 * @brief 创建 Arc-PIC 风格电场求解器。
 * @param grid_config 网格配置。
 * @param solver_config 求解器配置。
 * @return 求解器智能指针。
 */
std::unique_ptr<ElectricFieldSolver>
ElectricFieldSolverFactory::createArcPICSolver(const GridConfig& grid_config,
                                               const SolverConfig& solver_config)
{
    return std::make_unique<ElectricFieldSolver>(grid_config, solver_config);
}

/**
 * @brief 创建标准求解器。
 * @details
 * 使用给定区域和网格尺寸生成默认有限差分求解器，并启用预条件设置。
 *
 * @param domain_min 计算域最小坐标。
 * @param domain_max 计算域最大坐标。
 * @param grid_size 网格尺寸向量，期望长度为 3。
 * @return 求解器智能指针。
 */
std::unique_ptr<ElectricFieldSolver>
ElectricFieldSolverFactory::createStandardSolver(const SCDAT::Geometry::Point3D& domain_min,
                                                 const SCDAT::Geometry::Point3D& domain_max,
                                                 const std::vector<int>& grid_size)
{
    GridConfig grid;
    grid.domain_min = domain_min;
    grid.domain_max = domain_max;
    if (grid_size.size() == 3)
    {
        grid.grid_size = grid_size;
    }

    SolverConfig solver;
    solver.solver_type = SolverType::FINITE_DIFFERENCE;
    solver.enable_preconditioner = true;
    return std::make_unique<ElectricFieldSolver>(grid, solver);
}

/**
 * @brief 从配置文件创建求解器。
 * @details
 * 读取简单键值配置（如 nx/ny/nz/tol/max_iter），覆盖默认配置后构造求解器。
 * 若文件不可读，则返回基于默认配置构建的求解器。
 *
 * @param config_file 配置文件路径。
 * @return 求解器智能指针。
 */
std::unique_ptr<ElectricFieldSolver>
ElectricFieldSolverFactory::createFromConfig(const std::string& config_file)
{
    GridConfig grid;
    SolverConfig solver;

    std::ifstream fin(config_file);
    if (!fin)
    {
        return std::make_unique<ElectricFieldSolver>(grid, solver);
    }

    std::string line;
    while (std::getline(fin, line))
    {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key))
        {
            continue;
        }

        if (key == "nx")
        {
            iss >> grid.grid_size[0];
        }
        else if (key == "ny")
        {
            iss >> grid.grid_size[1];
        }
        else if (key == "nz")
        {
            iss >> grid.grid_size[2];
        }
        else if (key == "tol")
        {
            iss >> solver.tolerance;
        }
        else if (key == "max_iter")
        {
            iss >> solver.max_iterations;
        }
    }

    return std::make_unique<ElectricFieldSolver>(grid, solver);
}

} // namespace Field
} // namespace SCDAT
