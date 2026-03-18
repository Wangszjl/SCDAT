/**
 * @file Matrix.cpp
 * @brief 大规模矩阵计算实现
 */

#include "../include/Matrix.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace SCDAT
{
namespace Geometry
{

Matrix::Matrix() noexcept : rows_(0), cols_(0), data_() {}

Matrix::Matrix(std::size_t rows, std::size_t cols, double initValue)
    : rows_(rows), cols_(cols), data_(rows * cols, initValue)
{
    if (rows == 0 || cols == 0)
    {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
}

Matrix::Matrix(std::size_t rows, std::size_t cols, const std::vector<double>& rowMajorData)
    : rows_(rows), cols_(cols), data_(rowMajorData)
{
    if (rows == 0 || cols == 0)
    {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }

    if (rowMajorData.size() != rows * cols)
    {
        throw std::invalid_argument("Matrix data size does not match dimensions");
    }
}

Matrix Matrix::identity(std::size_t size)
{
    Matrix m(size, size, 0.0);
    for (std::size_t i = 0; i < size; ++i)
    {
        m(i, i) = 1.0;
    }
    return m;
}

Matrix Matrix::blockDiagonalFrom3x3(const std::vector<Matrix3x3>& blocks)
{
    if (blocks.empty())
    {
        throw std::invalid_argument("At least one Matrix3x3 block is required");
    }

    const std::size_t dimension = blocks.size() * 3;
    Matrix result(dimension, dimension, 0.0);

    for (std::size_t blockIndex = 0; blockIndex < blocks.size(); ++blockIndex)
    {
        const Matrix3x3& block = blocks[blockIndex];
        const std::size_t base = blockIndex * 3;

        for (std::size_t r = 0; r < 3; ++r)
        {
            for (std::size_t c = 0; c < 3; ++c)
            {
                result(base + r, base + c) = block(static_cast<int>(r), static_cast<int>(c));
            }
        }
    }

    return result;
}

std::size_t Matrix::rows() const noexcept
{
    return rows_;
}

std::size_t Matrix::cols() const noexcept
{
    return cols_;
}

bool Matrix::empty() const noexcept
{
    return rows_ == 0 || cols_ == 0;
}

std::size_t Matrix::index(std::size_t row, std::size_t col) const
{
    if (row >= rows_ || col >= cols_)
    {
        throw std::out_of_range("Matrix index out of range");
    }
    return row * cols_ + col;
}

double& Matrix::operator()(std::size_t row, std::size_t col)
{
    return data_[index(row, col)];
}

const double& Matrix::operator()(std::size_t row, std::size_t col) const
{
    return data_[index(row, col)];
}

void Matrix::ensureSameShape(const Matrix& other, const char* opName) const
{
    if (rows_ != other.rows_ || cols_ != other.cols_)
    {
        throw std::invalid_argument(std::string(opName) + ": matrix shape mismatch");
    }
}

Matrix Matrix::operator+(const Matrix& other) const
{
    ensureSameShape(other, "addition");

    Matrix result(rows_, cols_, 0.0);
    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        result.data_[i] = data_[i] + other.data_[i];
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const
{
    ensureSameShape(other, "subtraction");

    Matrix result(rows_, cols_, 0.0);
    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        result.data_[i] = data_[i] - other.data_[i];
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const
{
    if (cols_ != other.rows_)
    {
        throw std::invalid_argument("matrix multiplication shape mismatch");
    }

    Matrix result(rows_, other.cols_, 0.0);

    for (std::size_t i = 0; i < rows_; ++i)
    {
        for (std::size_t k = 0; k < cols_; ++k)
        {
            const double aik = (*this)(i, k);
            if (aik == 0.0)
            {
                continue;
            }
            for (std::size_t j = 0; j < other.cols_; ++j)
            {
                result(i, j) += aik * other(k, j);
            }
        }
    }

    return result;
}

Matrix Matrix::operator*(double scalar) const
{
    Matrix result(rows_, cols_, 0.0);
    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        result.data_[i] = data_[i] * scalar;
    }
    return result;
}

Matrix& Matrix::operator+=(const Matrix& other)
{
    ensureSameShape(other, "addition assignment");

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        data_[i] += other.data_[i];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& other)
{
    ensureSameShape(other, "subtraction assignment");

    for (std::size_t i = 0; i < data_.size(); ++i)
    {
        data_[i] -= other.data_[i];
    }
    return *this;
}

Matrix& Matrix::operator*=(double scalar)
{
    for (double& value : data_)
    {
        value *= scalar;
    }
    return *this;
}

Matrix Matrix::transpose() const
{
    Matrix result(cols_, rows_, 0.0);
    for (std::size_t i = 0; i < rows_; ++i)
    {
        for (std::size_t j = 0; j < cols_; ++j)
        {
            result(j, i) = (*this)(i, j);
        }
    }
    return result;
}

double Matrix::frobeniusNorm() const
{
    double sum = 0.0;
    for (double value : data_)
    {
        sum += value * value;
    }
    return std::sqrt(sum);
}

std::vector<double> Matrix::multiply(const std::vector<double>& vec) const
{
    if (vec.size() != cols_)
    {
        throw std::invalid_argument("matrix-vector multiplication shape mismatch");
    }

    std::vector<double> result(rows_, 0.0);
    for (std::size_t i = 0; i < rows_; ++i)
    {
        double acc = 0.0;
        for (std::size_t j = 0; j < cols_; ++j)
        {
            acc += (*this)(i, j) * vec[j];
        }
        result[i] = acc;
    }

    return result;
}

std::vector<Vector3D> Matrix::multiplyBlocks3(const std::vector<Vector3D>& vecBlocks) const
{
    if (rows_ % 3 != 0 || cols_ % 3 != 0)
    {
        throw std::invalid_argument("matrix dimensions must be multiples of 3 for block-vector multiplication");
    }

    const std::size_t expectedBlocks = cols_ / 3;
    if (vecBlocks.size() != expectedBlocks)
    {
        throw std::invalid_argument("block-vector size mismatch");
    }

    std::vector<double> flatVector(cols_, 0.0);
    for (std::size_t i = 0; i < vecBlocks.size(); ++i)
    {
        flatVector[i * 3] = vecBlocks[i].x();
        flatVector[i * 3 + 1] = vecBlocks[i].y();
        flatVector[i * 3 + 2] = vecBlocks[i].z();
    }

    const std::vector<double> flatResult = multiply(flatVector);

    std::vector<Vector3D> result(rows_ / 3, Vector3D::zero());
    for (std::size_t i = 0; i < result.size(); ++i)
    {
        result[i] = Vector3D(flatResult[i * 3], flatResult[i * 3 + 1], flatResult[i * 3 + 2]);
    }

    return result;
}

std::vector<double> Matrix::rowMajorData() const
{
    return data_;
}

std::string Matrix::toString() const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "Matrix(" << rows_ << "x" << cols_ << ")\n";
    for (std::size_t i = 0; i < rows_; ++i)
    {
        oss << "  [";
        for (std::size_t j = 0; j < cols_; ++j)
        {
            oss << std::setw(10) << (*this)(i, j);
            if (j + 1 < cols_)
            {
                oss << ", ";
            }
        }
        oss << "]";
        if (i + 1 < rows_)
        {
            oss << "\n";
        }
    }
    return oss.str();
}

Matrix operator*(double scalar, const Matrix& matrix)
{
    return matrix * scalar;
}

} // namespace Geometry
} // namespace SCDAT
