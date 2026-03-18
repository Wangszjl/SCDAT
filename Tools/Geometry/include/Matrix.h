/**
 * @file Matrix.h
 * @brief 大规模矩阵计算头文件
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026年3月18日
 *
 * 提供动态尺寸矩阵（m×n）运算能力，并支持基于Matrix3x3构建分块大矩阵。
 */

#ifndef SCDAT_MATRIX_H
#define SCDAT_MATRIX_H

#include "Matrix3x3.h"
#include "Vector3D.h"
#include <cstddef>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Geometry
{

class Matrix
{
  public:
    Matrix() noexcept;
    Matrix(std::size_t rows, std::size_t cols, double initValue = 0.0);
    Matrix(std::size_t rows, std::size_t cols, const std::vector<double>& rowMajorData);

    static Matrix identity(std::size_t size);
    static Matrix blockDiagonalFrom3x3(const std::vector<Matrix3x3>& blocks);

    std::size_t rows() const noexcept;
    std::size_t cols() const noexcept;
    bool empty() const noexcept;

    double& operator()(std::size_t row, std::size_t col);
    const double& operator()(std::size_t row, std::size_t col) const;

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;

    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix& operator*=(double scalar);

    Matrix transpose() const;
    double frobeniusNorm() const;

    std::vector<double> multiply(const std::vector<double>& vec) const;
    std::vector<Vector3D> multiplyBlocks3(const std::vector<Vector3D>& vecBlocks) const;

    std::vector<double> rowMajorData() const;
    std::string toString() const;

  private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;

    std::size_t index(std::size_t row, std::size_t col) const;
    void ensureSameShape(const Matrix& other, const char* opName) const;
};

Matrix operator*(double scalar, const Matrix& matrix);

} // namespace Geometry
} // namespace SCDAT

#endif // SCDAT_MATRIX_H
