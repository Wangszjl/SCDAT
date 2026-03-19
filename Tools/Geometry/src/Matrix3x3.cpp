/**
 * @file Matrix3x3.cpp
 * @brief 3×3矩阵数学运算实现
 * @author PIC-SPIS开发团队
 * @version V0.0.1
 * @date 2026年3月18日 14:30:17
 *
 * 实现3×3矩阵的基本数学运算、线性代数操作、矩阵分解和几何变换功能。
 * 支持矩阵运算、求逆、特征值计算、旋转变换等核心功能。
 *
 *
 * @note 矩阵采用行主序存储，索引从0开始
 * @note 所有旋转角度使用弧度制
 */

#include "../include/Matrix3x3.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SCDAT
{
namespace Geometry
{

namespace
{
constexpr double kEigenEpsilon = 1e-12;
constexpr double kPi = 3.14159265358979323846;

struct CharacteristicCoefficients
{
    double c1; // tr(A)
    double c2; // 1/2 * (tr(A)^2 - tr(A^2))
    double c3; // det(A)
};

CharacteristicCoefficients computeCharacteristicCoefficients(const Matrix3x3& A)
{
    const double c1 = A.trace();

    double traceA2 = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            traceA2 += A(i, j) * A(j, i);
        }
    }

    const double c2 = 0.5 * (c1 * c1 - traceA2);
    const double c3 = A.determinant();
    return {c1, c2, c3};
}

double evaluateCharacteristicPolynomial(const CharacteristicCoefficients& coeff, double lambda)
{
    // p(lambda) = lambda^3 - c1*lambda^2 + c2*lambda - c3
    return ((lambda - coeff.c1) * lambda + coeff.c2) * lambda - coeff.c3;
}

Vector3D dominantNullVector(const Matrix3x3& A, double lambda, double tol)
{
    const Matrix3x3 shifted = A - Matrix3x3::identity() * lambda;

    const Vector3D r0 = shifted.getRow(0);
    const Vector3D r1 = shifted.getRow(1);
    const Vector3D r2 = shifted.getRow(2);

    const Vector3D c01 = Vector3D(r0.cross(r1));
    const Vector3D c02 = Vector3D(r0.cross(r2));
    const Vector3D c12 = Vector3D(r1.cross(r2));

    const double n01 = c01.magnitudeSquared();
    const double n02 = c02.magnitudeSquared();
    const double n12 = c12.magnitudeSquared();

    Vector3D best = c01;
    double bestNorm = n01;
    if (n02 > bestNorm)
    {
        best = c02;
        bestNorm = n02;
    }
    if (n12 > bestNorm)
    {
        best = c12;
        bestNorm = n12;
    }

    if (bestNorm > tol * tol)
    {
        return best.safeNormalized();
    }

    const std::array<Vector3D, 3> rows = {r0, r1, r2};
    int bestRowIndex = 0;
    double bestRowNorm = rows[0].magnitudeSquared();
    for (int k = 1; k < 3; ++k)
    {
        const double n = rows[k].magnitudeSquared();
        if (n > bestRowNorm)
        {
            bestRowNorm = n;
            bestRowIndex = k;
        }
    }

    if (bestRowNorm > tol * tol)
    {
        const Vector3D row = rows[bestRowIndex];
        Vector3D candidate;
        if (std::abs(row.x()) <= std::abs(row.y()) && std::abs(row.x()) <= std::abs(row.z()))
        {
            candidate = Vector3D::unitX();
        }
        else if (std::abs(row.y()) <= std::abs(row.z()))
        {
            candidate = Vector3D::unitY();
        }
        else
        {
            candidate = Vector3D::unitZ();
        }

        const Vector3D v = Vector3D(row.cross(candidate)).safeNormalized();
        if (v.magnitudeSquared() > tol * tol)
        {
            return v;
        }
    }

    // A - lambda I 接近零矩阵时，任意非零向量都在零空间中。
    return Vector3D::unitX();
}

} // namespace

// ============================================================================
// 构造函数实现
// ============================================================================

/**
 * @brief 默认构造函数，创建零矩阵
 *
 * 初始化所有元素为0.0的3×3矩阵。
 */
Matrix3x3::Matrix3x3() noexcept
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            data_[i][j] = 0.0;
        }
    }
}

/**
 * @brief 从二维数组构造矩阵
 *
 * @param data 3×3二维数组数据
 */
Matrix3x3::Matrix3x3(const std::array<std::array<double, 3>, 3>& data) noexcept : data_(data) {}

/**
 * @brief 从9个标量值构造矩阵
 *
 * 按行主序排列：
 * Row-major order:
 * [m00 m01 m02]
 * [m10 m11 m12]
 * [m20 m21 m22]
 *
 * @param m00,m01,m02 第一行元素 / First row elements
 * @param m10,m11,m12 第二行元素 / Second row elements
 * @param m20,m21,m22 第三行元素 / Third row elements
 */
Matrix3x3::Matrix3x3(double m00, double m01, double m02, double m10, double m11, double m12,
                     double m20, double m21, double m22) noexcept
{
    data_[0] = {m00, m01, m02};
    data_[1] = {m10, m11, m12};
    data_[2] = {m20, m21, m22};
}

/**
 * @brief 从三个行向量构造矩阵
 * Construct matrix from three row vectors
 *
 * @param row0 第一行向量 / First row vector
 * @param row1 第二行向量 / Second row vector
 * @param row2 第三行向量 / Third row vector
 */
Matrix3x3::Matrix3x3(const Vector3D& row0, const Vector3D& row1, const Vector3D& row2) noexcept
{
    data_[0] = {row0.x(), row0.y(), row0.z()};
    data_[1] = {row1.x(), row1.y(), row1.z()};
    data_[2] = {row2.x(), row2.y(), row2.z()};
}

// 元素访问
double& Matrix3x3::operator()(int row, int col)
{
    if (row < 0 || row >= 3 || col < 0 || col >= 3)
    {
        throw std::out_of_range("Matrix3x3 index out of range");
    }
    return data_[static_cast<size_t>(row)][static_cast<size_t>(col)];
}

const double& Matrix3x3::operator()(int row, int col) const
{
    if (row < 0 || row >= 3 || col < 0 || col >= 3)
    {
        throw std::out_of_range("Matrix3x3 index out of range");
    }
    return data_[static_cast<size_t>(row)][static_cast<size_t>(col)];
}

// 行和列访问
Vector3D Matrix3x3::getRow(int row) const
{
    if (row < 0 || row >= 3)
    {
        throw std::out_of_range("Matrix3x3 row index out of range");
    }
    size_t r = static_cast<size_t>(row);
    return Vector3D(data_[r][0], data_[r][1], data_[r][2]);
}

Vector3D Matrix3x3::getColumn(int col) const
{
    if (col < 0 || col >= 3)
    {
        throw std::out_of_range("Matrix3x3 column index out of range");
    }
    size_t c = static_cast<size_t>(col);
    return Vector3D(data_[0][c], data_[1][c], data_[2][c]);
}

void Matrix3x3::setRow(int row, const Vector3D& vec)
{
    if (row < 0 || row >= 3)
    {
        throw std::out_of_range("Matrix3x3 row index out of range");
    }
    size_t r = static_cast<size_t>(row);
    data_[r][0] = vec.x();
    data_[r][1] = vec.y();
    data_[r][2] = vec.z();
}

void Matrix3x3::setColumn(int col, const Vector3D& vec)
{
    if (col < 0 || col >= 3)
    {
        throw std::out_of_range("Matrix3x3 column index out of range");
    }
    size_t c = static_cast<size_t>(col);
    data_[0][c] = vec.x();
    data_[1][c] = vec.y();
    data_[2][c] = vec.z();
}

// 矩阵运算
Matrix3x3 Matrix3x3::operator+(const Matrix3x3& other) const
{
    Matrix3x3 result;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            result.data_[i][j] = data_[i][j] + other.data_[i][j];
        }
    }
    return result;
}

Matrix3x3 Matrix3x3::operator-(const Matrix3x3& other) const
{
    Matrix3x3 result;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            result.data_[i][j] = data_[i][j] - other.data_[i][j];
        }
    }
    return result;
}

Matrix3x3 Matrix3x3::operator*(const Matrix3x3& other) const
{
    Matrix3x3 result;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            result.data_[i][j] = 0.0;
            for (size_t k = 0; k < 3; ++k)
            {
                result.data_[i][j] += data_[i][k] * other.data_[k][j];
            }
        }
    }
    return result;
}

Matrix3x3 Matrix3x3::operator*(double scalar) const
{
    Matrix3x3 result;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            result.data_[i][j] = data_[i][j] * scalar;
        }
    }
    return result;
}

Matrix3x3 Matrix3x3::operator/(double scalar) const
{
    if (std::abs(scalar) < EPSILON)
    {
        throw std::runtime_error("Division by zero in Matrix3x3");
    }
    return *this * (1.0 / scalar);
}

// 复合赋值运算符
Matrix3x3& Matrix3x3::operator+=(const Matrix3x3& other)
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            data_[i][j] += other.data_[i][j];
        }
    }
    return *this;
}

Matrix3x3& Matrix3x3::operator-=(const Matrix3x3& other)
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            data_[i][j] -= other.data_[i][j];
        }
    }
    return *this;
}

Matrix3x3& Matrix3x3::operator*=(const Matrix3x3& other)
{
    *this = *this * other;
    return *this;
}

Matrix3x3& Matrix3x3::operator*=(double scalar)
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            data_[i][j] *= scalar;
        }
    }
    return *this;
}

Matrix3x3& Matrix3x3::operator/=(double scalar)
{
    if (std::abs(scalar) < EPSILON)
    {
        throw std::runtime_error("Division by zero in Matrix3x3");
    }
    return *this *= (1.0 / scalar);
}

// 矩阵-向量运算
Vector3D Matrix3x3::operator*(const Vector3D& vec) const
{
    return Vector3D(data_[0][0] * vec.x() + data_[0][1] * vec.y() + data_[0][2] * vec.z(),
                    data_[1][0] * vec.x() + data_[1][1] * vec.y() + data_[1][2] * vec.z(),
                    data_[2][0] * vec.x() + data_[2][1] * vec.y() + data_[2][2] * vec.z());
}

// 比较运算符
bool Matrix3x3::operator==(const Matrix3x3& other) const
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            if (std::abs(data_[i][j] - other.data_[i][j]) > EPSILON)
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix3x3::operator!=(const Matrix3x3& other) const
{
    return !(*this == other);
}

// 矩阵属性
double Matrix3x3::determinant() const
{
    return data_[0][0] * (data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1]) -
           data_[0][1] * (data_[1][0] * data_[2][2] - data_[1][2] * data_[2][0]) +
           data_[0][2] * (data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0]);
}

double Matrix3x3::trace() const
{
    return data_[0][0] + data_[1][1] + data_[2][2];
}

Matrix3x3 Matrix3x3::transpose() const
{
    Matrix3x3 result;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            result.data_[i][j] = data_[j][i];
        }
    }
    return result;
}

Matrix3x3 Matrix3x3::inverse() const
{
    double det = determinant();
    if (std::abs(det) < EPSILON)
    {
        throw std::runtime_error("Matrix is singular and cannot be inverted");
    }

    Matrix3x3 result;
    double inv_det = 1.0 / det;

    // 计算伴随矩阵的转置
    result.data_[0][0] = (data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1]) * inv_det;
    result.data_[0][1] = (data_[0][2] * data_[2][1] - data_[0][1] * data_[2][2]) * inv_det;
    result.data_[0][2] = (data_[0][1] * data_[1][2] - data_[0][2] * data_[1][1]) * inv_det;

    result.data_[1][0] = (data_[1][2] * data_[2][0] - data_[1][0] * data_[2][2]) * inv_det;
    result.data_[1][1] = (data_[0][0] * data_[2][2] - data_[0][2] * data_[2][0]) * inv_det;
    result.data_[1][2] = (data_[0][2] * data_[1][0] - data_[0][0] * data_[1][2]) * inv_det;

    result.data_[2][0] = (data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0]) * inv_det;
    result.data_[2][1] = (data_[0][1] * data_[2][0] - data_[0][0] * data_[2][1]) * inv_det;
    result.data_[2][2] = (data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0]) * inv_det;

    return result;
}

// 矩阵属性检查
bool Matrix3x3::isSymmetric(double tolerance) const
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            if (std::abs(data_[i][j] - data_[j][i]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix3x3::isOrthogonal(double tolerance) const
{
    Matrix3x3 product = (*this) * transpose();
    return product.isIdentity(tolerance);
}

bool Matrix3x3::isIdentity(double tolerance) const
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(data_[i][j] - expected) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix3x3::isSingular(double tolerance) const
{
    return std::abs(determinant()) < tolerance;
}

// 特征值计算（基于特征多项式的解析解，返回实部）
std::array<double, 3> Matrix3x3::eigenvalues() const
{
    const CharacteristicCoefficients coeff = computeCharacteristicCoefficients(*this);

    // lambda^3 - c1*lambda^2 + c2*lambda - c3 = 0
    // 变换 lambda = t + c1/3 -> t^3 + pt + q = 0
    const double p = coeff.c2 - coeff.c1 * coeff.c1 / 3.0;
    const double q =
        -2.0 * coeff.c1 * coeff.c1 * coeff.c1 / 27.0 + coeff.c1 * coeff.c2 / 3.0 - coeff.c3;

    const double disc = (q * q) / 4.0 + (p * p * p) / 27.0;
    const double discScale = std::abs((q * q) / 4.0) + std::abs((p * p * p) / 27.0) + 1.0;
    const double discTol = 64.0 * kEigenEpsilon * discScale;

    const double shift = coeff.c1 / 3.0;
    std::array<double, 3> eigenvals{};

    if (disc > discTol)
    {
        // 1个实根 + 1对共轭复根（接口为实数数组，后两项返回复根实部）
        const double sqrtDisc = std::sqrt(disc);
        const double u = std::cbrt(-q / 2.0 + sqrtDisc);
        const double v = std::cbrt(-q / 2.0 - sqrtDisc);

        const double tReal = u + v;
        const double realPart = -0.5 * tReal;

        eigenvals[0] = tReal + shift;
        eigenvals[1] = realPart + shift;
        eigenvals[2] = realPart + shift;
    }
    else if (std::abs(disc) <= discTol)
    {
        if (std::abs(p) <= kEigenEpsilon && std::abs(q) <= kEigenEpsilon)
        {
            eigenvals[0] = shift;
            eigenvals[1] = shift;
            eigenvals[2] = shift;
        }
        else
        {
            const double u = std::cbrt(-q / 2.0);
            eigenvals[0] = 2.0 * u + shift;
            eigenvals[1] = -u + shift;
            eigenvals[2] = -u + shift;
        }
    }
    else
    {
        // 3个实根
        const double phiArgDen = std::sqrt(-(p * p * p) / 27.0);
        double phiArg = 0.0;
        if (phiArgDen > kEigenEpsilon)
        {
            phiArg = -q / (2.0 * phiArgDen);
        }
        phiArg = std::max(-1.0, std::min(1.0, phiArg));

        const double phi = std::acos(phiArg);
        const double t = 2.0 * std::sqrt(-p / 3.0);

        eigenvals[0] = t * std::cos(phi / 3.0) + shift;
        eigenvals[1] = t * std::cos((phi + 2.0 * kPi) / 3.0) + shift;
        eigenvals[2] = t * std::cos((phi + 4.0 * kPi) / 3.0) + shift;
    }

    for (double& val : eigenvals)
    {
        if (std::abs(val) < kEigenEpsilon)
        {
            val = 0.0;
        }
    }

    std::sort(eigenvals.begin(), eigenvals.end(), std::greater<double>());

    return eigenvals;
}

std::array<Vector3D, 3> Matrix3x3::eigenvectors() const
{
    const auto lambdas = eigenvalues();
    const CharacteristicCoefficients coeff = computeCharacteristicCoefficients(*this);

    std::array<Vector3D, 3> vectors = {Vector3D::zero(), Vector3D::zero(), Vector3D::zero()};

    const double normA = std::max(1.0, frobeniusNorm());
    const double lambdaTol = 128.0 * kEigenEpsilon * normA;

    for (size_t k = 0; k < 3; ++k)
    {
        const double lambda = lambdas[k];

        // 复根场景中，eigenvalues()会返回复根实部；通过特征多项式残差做过滤。
        const double polyResidual = std::abs(evaluateCharacteristicPolynomial(coeff, lambda));
        const double polyScale = 1.0 + std::abs(coeff.c1) + std::abs(coeff.c2) +
                                 std::abs(coeff.c3) + std::abs(lambda * lambda * lambda);

        if (polyResidual > 256.0 * kEigenEpsilon * polyScale)
        {
            continue;
        }

        const Vector3D v = dominantNullVector(*this, lambda, lambdaTol);
        if (v.magnitudeSquared() <= lambdaTol * lambdaTol)
        {
            continue;
        }

        const Vector3D residual = (*this) * v - v * lambda;
        const double residualTol = 256.0 * kEigenEpsilon * (normA + std::abs(lambda));
        if (residual.magnitude() <= residualTol)
        {
            vectors[k] = v.safeNormalized();
        }
    }

    return vectors;
}

// 矩阵范数
double Matrix3x3::frobeniusNorm() const
{
    double sum = 0.0;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            sum += data_[i][j] * data_[i][j];
        }
    }
    return std::sqrt(sum);
}

double Matrix3x3::infinityNorm() const
{
    double max_row_sum = 0.0;
    for (size_t i = 0; i < 3; ++i)
    {
        double row_sum = 0.0;
        for (size_t j = 0; j < 3; ++j)
        {
            row_sum += std::abs(data_[i][j]);
        }
        max_row_sum = std::max(max_row_sum, row_sum);
    }
    return max_row_sum;
}

double Matrix3x3::oneNorm() const
{
    double max_col_sum = 0.0;
    for (size_t j = 0; j < 3; ++j)
    {
        double col_sum = 0.0;
        for (size_t i = 0; i < 3; ++i)
        {
            col_sum += std::abs(data_[i][j]);
        }
        max_col_sum = std::max(max_col_sum, col_sum);
    }
    return max_col_sum;
}

// 就地操作
void Matrix3x3::transposeInPlace()
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = i + 1; j < 3; ++j)
        {
            std::swap(data_[i][j], data_[j][i]);
        }
    }
}

void Matrix3x3::invertInPlace()
{
    *this = inverse();
}

// 静态工厂方法
Matrix3x3 Matrix3x3::identity()
{
    return Matrix3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
}

Matrix3x3 Matrix3x3::zero()
{
    return Matrix3x3();
}

Matrix3x3 Matrix3x3::diagonal(double d0, double d1, double d2)
{
    return Matrix3x3(d0, 0.0, 0.0, 0.0, d1, 0.0, 0.0, 0.0, d2);
}

Matrix3x3 Matrix3x3::diagonal(const Vector3D& diag)
{
    return diagonal(diag.x(), diag.y(), diag.z());
}

// 旋转矩阵
Matrix3x3 Matrix3x3::rotationX(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);
    return Matrix3x3(1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c);
}

Matrix3x3 Matrix3x3::rotationY(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);
    return Matrix3x3(c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c);
}

Matrix3x3 Matrix3x3::rotationZ(double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);
    return Matrix3x3(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
}

// 字符串表示
std::string Matrix3x3::toString() const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "Matrix3x3(\n";
    for (size_t i = 0; i < 3; ++i)
    {
        oss << "  [";
        for (size_t j = 0; j < 3; ++j)
        {
            oss << std::setw(10) << data_[i][j];
            if (j < 2)
                oss << ", ";
        }
        oss << "]";
        if (i < 2)
            oss << ",";
        oss << "\n";
    }
    oss << ")";
    return oss.str();
}

// 全局运算符
/**
 * @brief 一元负号运算符
 */
Matrix3x3 Matrix3x3::operator-() const
{
    Matrix3x3 result;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            result.data_[i][j] = -data_[i][j];
        }
    }
    return result;
}

/**
 * @brief 缩放矩阵
 */
Matrix3x3 Matrix3x3::scaling(double sx, double sy, double sz)
{
    return Matrix3x3(sx, 0.0, 0.0, 0.0, sy, 0.0, 0.0, 0.0, sz);
}

/**
 * @brief 均匀缩放矩阵
 */
Matrix3x3 Matrix3x3::uniformScaling(double scale)
{
    return scaling(scale, scale, scale);
}

/**
 * @brief 反射矩阵
 */
Matrix3x3 Matrix3x3::reflection(const Vector3D& normal)
{
    Vector3D n = normal.normalized();
    return Matrix3x3(1.0 - 2.0 * n.x() * n.x(), -2.0 * n.x() * n.y(), -2.0 * n.x() * n.z(),
                     -2.0 * n.y() * n.x(), 1.0 - 2.0 * n.y() * n.y(), -2.0 * n.y() * n.z(),
                     -2.0 * n.z() * n.x(), -2.0 * n.z() * n.y(), 1.0 - 2.0 * n.z() * n.z());
}

Matrix3x3 operator*(double scalar, const Matrix3x3& matrix)
{
    return matrix * scalar;
}

/**
 * @brief 从三个列向量构造矩阵
 */
Matrix3x3 Matrix3x3::fromColumns(const Vector3D& col0, const Vector3D& col1,
                                 const Vector3D& col2) noexcept
{
    return Matrix3x3(col0.x(), col1.x(), col2.x(), col0.y(), col1.y(), col2.y(), col0.z(), col1.z(),
                     col2.z());
}

/**
 * @brief 检查是否为零矩阵
 */
bool Matrix3x3::isZero(double tolerance) const
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            if (std::abs(data_[i][j]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 检查是否为对角矩阵
 */
bool Matrix3x3::isDiagonal(double tolerance) const
{
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            if (i != j && std::abs(data_[i][j]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 检查是否为上三角矩阵
 */
bool Matrix3x3::isUpperTriangular(double tolerance) const
{
    for (size_t i = 1; i < 3; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            if (std::abs(data_[i][j]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 检查是否为下三角矩阵
 */
bool Matrix3x3::isLowerTriangular(double tolerance) const
{
    for (size_t i = 0; i < 2; ++i)
    {
        for (size_t j = i + 1; j < 3; ++j)
        {
            if (std::abs(data_[i][j]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief 检查是否可逆
 */
bool Matrix3x3::isInvertible(double tolerance) const
{
    return !isSingular(tolerance);
}

/**
 * @brief 计算谱范数（最大奇异值）
 */
double Matrix3x3::spectralNorm() const
{
    // 对于任意实矩阵，||A||_2 = sqrt(lambda_max(A^T A))。
    const Matrix3x3 AtA = transpose() * (*this);
    const auto evals = AtA.eigenvalues();
    const double maxEval = std::max({0.0, evals[0], evals[1], evals[2]});
    return std::sqrt(maxEval);
}

/**
 * @brief 剪切矩阵XY
 */
Matrix3x3 Matrix3x3::shearXY(double shx, double shy)
{
    return Matrix3x3(1.0, shx, 0.0, shy, 1.0, 0.0, 0.0, 0.0, 1.0);
}

/**
 * @brief 剪切矩阵XZ
 */
Matrix3x3 Matrix3x3::shearXZ(double shx, double shz)
{
    return Matrix3x3(1.0, 0.0, shx, 0.0, 1.0, 0.0, shz, 0.0, 1.0);
}

/**
 * @brief 剪切矩阵YZ
 */
Matrix3x3 Matrix3x3::shearYZ(double shy, double shz)
{
    return Matrix3x3(1.0, 0.0, 0.0, 0.0, 1.0, shy, 0.0, shz, 1.0);
}

/**
 * @brief 条件数计算
 */
double Matrix3x3::conditionNumber() const
{
    try
    {
        Matrix3x3 inv = inverse();
        return frobeniusNorm() * inv.frobeniusNorm();
    }
    catch (...)
    {
        return std::numeric_limits<double>::infinity();
    }
}

std::ostream& operator<<(std::ostream& os, const Matrix3x3& matrix)
{
    os << matrix.toString();
    return os;
}

} // namespace Geometry
} // namespace SCDAT
