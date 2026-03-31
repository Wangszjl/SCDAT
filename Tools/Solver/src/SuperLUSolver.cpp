/**
 * @file SuperLUSolver.cpp
 * @brief SuperLU稀疏直接求解器实现
 * @details 提供回退实现（Gauss-Seidel迭代）和SuperLU接口
 * @author PIC-SPIS Development Team
 * @version 1.0.0
 * @date 2025-12-17
 */

#include "SuperLUSolver.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <iostream>

namespace SCDAT {
namespace Solver {

//=============================================================================
// SuperLUSolver 实现
//=============================================================================

#ifdef HAVE_SUPERLU

bool SuperLUSolver::factorizeSuperLU(const CSRMatrix& matrix) {
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        // 验证矩阵
        if (!matrix.validate()) {
            throw std::runtime_error("Invalid CSR matrix");
        }
        
        int n = matrix.nrows;
        int nnz = matrix.nnz;
        
        if (params_.verbose) {
            std::cout << "SuperLU factorization starting...\n";
            std::cout << "  Matrix size: " << n << "x" << n << "\n";
            std::cout << "  Non-zeros: " << nnz << "\n";
        }
        
        // 1. 转换为SuperLU格式 (CompCol)
        double* values = const_cast<double*>(matrix.values.data());
        int* col_ind = const_cast<int*>(matrix.col_indices.data());
        int* row_ptr = const_cast<int*>(matrix.row_pointers.data());
        
        // 转换 CSR -> CSC (列压缩格�?
        std::vector<double> values_csc(nnz);
        std::vector<int> row_ind_csc(nnz);
        std::vector<int> col_ptr_csc(n + 1, 0);
        
        // 计算每列的非零元素数
        for (int i = 0; i < nnz; ++i) {
            col_ptr_csc[col_ind[i] + 1]++;
        }
        
        // 计算列指�?
        for (int i = 0; i < n; ++i) {
            col_ptr_csc[i + 1] += col_ptr_csc[i];
        }
        
        // 填充CSC数据
        std::vector<int> col_counts(n, 0);
        for (int row = 0; row < n; ++row) {
            for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx) {
                int col = col_ind[idx];
                int pos = col_ptr_csc[col] + col_counts[col];
                values_csc[pos] = values[idx];
                row_ind_csc[pos] = row;
                col_counts[col]++;
            }
        }
        
        // 2. 创建SuperLU矩阵
        SuperMatrix A;
        dCreate_CompCol_Matrix(&A, n, n, nnz,
                              values_csc.data(),
                              row_ind_csc.data(),
                              col_ptr_csc.data(),
                              SLU_NC, SLU_D, SLU_GE);
        
        // 3. 分配置换向量
        perm_c_ = new int[n];
        perm_r_ = new int[n];
        
        // 4. 设置选项
        superlu_options_t options;
        set_default_options(&options);
        options.ColPerm = COLAMD;  // 列排序算�?
        options.Equil = params_.equilibrate ? YES : NO;
        options.IterRefine = params_.use_iterative_refine ? SLU_DOUBLE : NOREFINE;
        
        // 5. 获取列置�?
        SuperMatrix AC;
        int* etree = new int[n];
        SuperLUStat_t stat;
        StatInit(&stat);
        
        get_perm_c(options.ColPerm, &A, perm_c_);
        sp_preorder(&options, &A, perm_c_, etree, &AC);
        
        // 6. 执行数值分�?
        GlobalLU_t Glu;
        int panel_size = sp_ienv(1);
        int relax = sp_ienv(2);
        int info;
        
        dgstrf(&options, &AC, relax, panel_size, etree, NULL, 0,
               perm_c_, perm_r_, &L_, &U_, &Glu, &stat, &info);
        
        if (info != 0) {
            delete[] etree;
            StatFree(&stat);
            Destroy_SuperMatrix_Store(&A);
            Destroy_CompCol_Permuted(&AC);
            throw std::runtime_error("SuperLU factorization failed, info=" + std::to_string(info));
        }
        
        // 7. 收集统计信息
        stats_.nnz_L = ((SCformat*)L_.Store)->nnz;
        stats_.nnz_U = ((NCformat*)U_.Store)->nnz;
        stats_.fill_ratio = static_cast<double>(stats_.nnz_L + stats_.nnz_U) / nnz;
        stats_.factorized = true;
        
        // 清理
        delete[] etree;
        StatFree(&stat);
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
        
        auto end = std::chrono::high_resolution_clock::now();
        stats_.factorization_time = std::chrono::duration<double>(end - start).count();
        
        if (params_.verbose) {
            std::cout << "SuperLU factorization complete!\n";
            std::cout << stats_.toString();
        }
        
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "SuperLU factorization error: " << e.what() << "\n";
        return false;
    }
}

bool SuperLUSolver::solveSuperLU(const std::vector<double>& rhs, std::vector<double>& solution) {
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        int n = static_cast<int>(rhs.size());
        solution.resize(n);
        
        // 创建右端项矩�?
        SuperMatrix B;
        std::vector<double> rhs_copy = rhs;
        dCreate_Dense_Matrix(&B, n, 1, rhs_copy.data(), n, SLU_DN, SLU_D, SLU_GE);
        
        // 求解
        SuperLUStat_t stat;
        StatInit(&stat);
        
        int info;
        trans_t trans = NOTRANS;
        dgstrs(trans, &L_, &U_, perm_c_, perm_r_, &B, &stat, &info);
        
        if (info != 0) {
            StatFree(&stat);
            Destroy_SuperMatrix_Store(&B);
            throw std::runtime_error("SuperLU solve failed, info=" + std::to_string(info));
        }
        
        // 复制�?
        std::copy(rhs_copy.begin(), rhs_copy.end(), solution.begin());
        
        // 可选：迭代精化
        if (params_.use_iterative_refine) {
            // TODO: 实现迭代精化
            stats_.refine_steps = 0;
        }
        
        StatFree(&stat);
        Destroy_SuperMatrix_Store(&B);
        
        auto end = std::chrono::high_resolution_clock::now();
        stats_.solve_time = std::chrono::duration<double>(end - start).count();
        
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "SuperLU solve error: " << e.what() << "\n";
        return false;
    }
}

void SuperLUSolver::cleanupSuperLU() {
    if (stats_.factorized) {
        Destroy_SuperNode_Matrix(&L_);
        Destroy_CompCol_Matrix(&U_);
        delete[] perm_c_;
        delete[] perm_r_;
        perm_c_ = nullptr;
        perm_r_ = nullptr;
        stats_.factorized = false;
    }
}

#endif // HAVE_SUPERLU

//=============================================================================
// 回退实现（Gauss-Seidel迭代�?
//=============================================================================

bool SuperLUSolver::factorizeFallback(const CSRMatrix& matrix) {
    auto start = std::chrono::high_resolution_clock::now();
    
    // 验证矩阵
    if (!matrix.validate()) {
        std::cerr << "Invalid CSR matrix\n";
        return false;
    }
    
    // 保存矩阵副本
    matrix_copy_ = std::make_shared<CSRMatrix>(matrix);
    
    stats_.factorized = true;
    
    auto end = std::chrono::high_resolution_clock::now();
    stats_.factorization_time = std::chrono::duration<double>(end - start).count();
    
    if (params_.verbose) {
        std::cout << "Fallback solver initialized (Gauss-Seidel)\n";
        std::cout << "  Matrix: " << matrix_copy_->info() << "\n";
        std::cout << "  WARNING: Not using SuperLU. Install SuperLU for better performance!\n";
    }
    
    return true;
}

bool SuperLUSolver::solveFallback(const std::vector<double>& rhs, std::vector<double>& solution) {
    auto start = std::chrono::high_resolution_clock::now();
    
    if (!matrix_copy_) {
        std::cerr << "Matrix not factorized\n";
        return false;
    }
    
    int n = matrix_copy_->nrows;
    solution.resize(n, 0.0);
    
    // Gauss-Seidel迭代
    const int max_iter = 1000;
    const double tol = 1e-6;
    
    bool converged = false;
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_change = 0.0;
        
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            double diag = 0.0;
            
            int row_start = matrix_copy_->row_pointers[i];
            int row_end = matrix_copy_->row_pointers[i + 1];
            
            for (int idx = row_start; idx < row_end; ++idx) {
                int j = matrix_copy_->col_indices[idx];
                double aij = matrix_copy_->values[idx];
                
                if (i == j) {
                    diag = aij;
                } else {
                    sum += aij * solution[j];
                }
            }
            
            if (std::abs(diag) < 1e-14) {
                std::cerr << "Zero diagonal at row " << i << "\n";
                return false;
            }
            
            double new_val = (rhs[i] - sum) / diag;
            double change = std::abs(new_val - solution[i]);
            max_change = std::max(max_change, change);
            solution[i] = new_val;
        }
        
        if (max_change < tol) {
            stats_.final_residual = max_change;
            converged = true;
            if (params_.verbose) {
                std::cout << "Gauss-Seidel converged in " << iter + 1 << " iterations\n";
            }
            break;
        }
        
        if (iter == max_iter - 1) {
            std::cerr << "WARNING: Gauss-Seidel did not converge (residual=" << max_change << ")\n";
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    stats_.solve_time = std::chrono::duration<double>(end - start).count();
    
    if (!converged && !solution.empty()) {
        double residual = 0.0;
        for (int i = 0; i < n; ++i) {
            double ax = 0.0;
            for (int idx = matrix_copy_->row_pointers[i]; idx < matrix_copy_->row_pointers[i + 1];
                 ++idx) {
                ax += matrix_copy_->values[idx] *
                      solution[static_cast<std::size_t>(matrix_copy_->col_indices[idx])];
            }
            const double diff = rhs[static_cast<std::size_t>(i)] - ax;
            residual += diff * diff;
        }
        stats_.final_residual = std::sqrt(residual);
    }

    return converged;
}

//=============================================================================
// PoissonSuperLUSolver 实现
//=============================================================================

PoissonSuperLUSolver::PoissonSuperLUSolver(int nx, int ny, int nz,
                                           double dx, double dy, double dz)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz),
      boundary_values_(static_cast<std::size_t>(std::max(0, nx * ny * nz)), 0.0) {
    buildLaplacianMatrix();
}

void PoissonSuperLUSolver::buildLaplacianMatrix() {
    // 3D拉普拉斯算子：∇² = ∂�?∂x² + ∂�?∂y² + ∂�?∂z²
    // 7点差分格�?
    
    int n = nx_ * ny_ * nz_;
    int nnz_estimate = 7 * n;  // 每个内部节点7个非零元�?
    
    laplacian_ = CSRMatrix(n, n, 0);
    laplacian_.values.reserve(nnz_estimate);
    laplacian_.col_indices.reserve(nnz_estimate);
    laplacian_.row_pointers.resize(n + 1, 0);
    
    double cx = 1.0 / (dx_ * dx_);
    double cy = 1.0 / (dy_ * dy_);
    double cz = 1.0 / (dz_ * dz_);
    double c_center = -2.0 * (cx + cy + cz);
    
    int row = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i, ++row) {
                laplacian_.row_pointers[row] = static_cast<int>(laplacian_.values.size());
                
                // 边界节点：Dirichlet B.C. (φ = 0)
                if (i == 0 || i == nx_-1 ||
                    j == 0 || j == ny_-1 ||
                    k == 0 || k == nz_-1) {
                    laplacian_.values.push_back(1.0);
                    laplacian_.col_indices.push_back(row);
                    continue;
                }
                
                // 内部节点�?点差�?
                int idx = i + j * nx_ + k * nx_ * ny_;
                
                // -x方向
                laplacian_.values.push_back(cx);
                laplacian_.col_indices.push_back(idx - 1);
                
                // -y方向
                laplacian_.values.push_back(cy);
                laplacian_.col_indices.push_back(idx - nx_);
                
                // -z方向
                laplacian_.values.push_back(cz);
                laplacian_.col_indices.push_back(idx - nx_ * ny_);
                
                // 中心
                laplacian_.values.push_back(c_center);
                laplacian_.col_indices.push_back(idx);
                
                // +z方向
                laplacian_.values.push_back(cz);
                laplacian_.col_indices.push_back(idx + nx_ * ny_);
                
                // +y方向
                laplacian_.values.push_back(cy);
                laplacian_.col_indices.push_back(idx + nx_);
                
                // +x方向
                laplacian_.values.push_back(cx);
                laplacian_.col_indices.push_back(idx + 1);
            }
        }
    }
    
    laplacian_.row_pointers[n] = static_cast<int>(laplacian_.values.size());
    laplacian_.nnz = static_cast<int>(laplacian_.values.size());
}

void PoissonSuperLUSolver::setBoundaryConditions(const std::vector<double>& boundary_values) {
    boundary_values_.assign(static_cast<std::size_t>(nx_ * ny_ * nz_), 0.0);
    const std::size_t count = std::min(boundary_values_.size(), boundary_values.size());
    std::copy_n(boundary_values.begin(), count, boundary_values_.begin());
}

bool PoissonSuperLUSolver::solve(const std::vector<double>& charge_density,
                                 std::vector<double>& potential) {
    // 1. 构造右端项: b = -ρ/ε₀
    const double eps0 = 8.854187817e-12;  // 真空介电常数
    const std::size_t system_size = static_cast<std::size_t>(nx_ * ny_ * nz_);
    std::vector<double> rhs(system_size, 0.0);
    for (std::size_t i = 0; i < system_size && i < charge_density.size(); ++i) {
        rhs[i] = -charge_density[i] / eps0;
    }

    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                if (i == 0 || i == nx_ - 1 ||
                    j == 0 || j == ny_ - 1 ||
                    k == 0 || k == nz_ - 1) {
                    const std::size_t idx =
                        static_cast<std::size_t>(i + j * nx_ + k * nx_ * ny_);
                    rhs[idx] = (idx < boundary_values_.size()) ? boundary_values_[idx] : 0.0;
                }
            }
        }
    }
    
    // 2. 首次求解需要分�?
    if (!solver_.isFactorized()) {
        if (!solver_.factorize(laplacian_)) {
            return false;
        }
    }
    
    // 3. 求解
    return solver_.solve(rhs, potential);
}

} // namespace Solver
} // namespace SCDAT
