/**
 * @file Matrix3x3.h
 * @brief 3×3矩阵数学运算头文件
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026年3月18日 14:25:19
 *
 * 定义3×3矩阵类，提供完整的矩阵数学运算和线性变换功能。
 * 支持矩阵代数、线性变换、特征值计算、矩阵分解等核心功能。
 *
 * @note 矩阵采用行主序存储，索引从0开始
 * @note 旋转角度使用弧度制，遵循右手定则
 *
 * 主要功能
 * - 矩阵基本运算
 * - 线性变换
 * - 旋转、缩放、反射
 * - 矩阵分解和特征值
 * - 数值稳定性优化
 * @ingroup ToolsModule
 */

#ifndef SCDAT_MATRIX3X3_H
#define SCDAT_MATRIX3X3_H

#include "Vector3D.h"
#include <array>
#include <string>

namespace SCDAT
{
namespace Geometry
{

/**
 * @brief 3×3矩阵类
 *
 * 用于三维空间的线性变换，包括旋转、缩放、反射等。
 * 采用行主序存储，支持与Vector3D的高效运算。
 *
 * @note 矩阵元素按行主序存储：data_[row][col]
 * @note 支持链式运算和就地操作优化
 */
class Matrix3x3
{
  private:
    std::array<std::array<double, 3>, 3> data_;

  public:
    // 构造函数
    Matrix3x3() noexcept; // 零矩阵
    Matrix3x3(const std::array<std::array<double, 3>, 3>& data) noexcept;
    Matrix3x3(double m00, double m01, double m02, double m10, double m11, double m12, double m20,
              double m21, double m22) noexcept;

    // 从三个行向量构造
    Matrix3x3(const Vector3D& row0, const Vector3D& row1, const Vector3D& row2) noexcept;

    // 从三个列向量构造
    static Matrix3x3 fromColumns(const Vector3D& col0, const Vector3D& col1,
                                 const Vector3D& col2) noexcept;

    // 模板化构造函数支持不同数值类型
    template <typename T>
    constexpr Matrix3x3(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) noexcept
        : data_{{{{static_cast<double>(m00), static_cast<double>(m01), static_cast<double>(m02)}},
                 {{static_cast<double>(m10), static_cast<double>(m11), static_cast<double>(m12)}},
                 {{static_cast<double>(m20), static_cast<double>(m21), static_cast<double>(m22)}}}}
    {
    }

    // 拷贝构造和赋值
    Matrix3x3(const Matrix3x3& other) = default;
    Matrix3x3& operator=(const Matrix3x3& other) = default;
    Matrix3x3(Matrix3x3&& other) noexcept = default;
    Matrix3x3& operator=(Matrix3x3&& other) noexcept = default;

    // 析构函数
    ~Matrix3x3() = default;

    // 元素访问
    double& operator()(int row, int col);
    const double& operator()(int row, int col) const;

    // 行和列访问
    Vector3D getRow(int row) const;
    Vector3D getColumn(int col) const;

    void setRow(int row, const Vector3D& vec);
    void setColumn(int col, const Vector3D& vec);

    // 矩阵运算
    Matrix3x3 operator+(const Matrix3x3& other) const;
    Matrix3x3 operator-(const Matrix3x3& other) const;
    Matrix3x3 operator*(const Matrix3x3& other) const;
    Matrix3x3 operator*(double scalar) const;
    Matrix3x3 operator/(double scalar) const;
    Matrix3x3 operator-() const; // 一元负号

    // 复合赋值运算符
    Matrix3x3& operator+=(const Matrix3x3& other);
    Matrix3x3& operator-=(const Matrix3x3& other);
    Matrix3x3& operator*=(const Matrix3x3& other);
    Matrix3x3& operator*=(double scalar);
    Matrix3x3& operator/=(double scalar);

    // 矩阵-向量运算
    Vector3D operator*(const Vector3D& vec) const;

    // 比较运算符
    bool operator==(const Matrix3x3& other) const;
    bool operator!=(const Matrix3x3& other) const;

    // 矩阵属性
    double determinant() const;
    double trace() const;
    Matrix3x3 transpose() const;
    Matrix3x3 inverse() const;

    // 矩阵检查
    bool isSymmetric(double tolerance = 1e-12) const;
    bool isOrthogonal(double tolerance = 1e-12) const;
    bool isIdentity(double tolerance = 1e-12) const;
    bool isSingular(double tolerance = 1e-12) const;
    bool isZero(double tolerance = 1e-12) const;
    bool isDiagonal(double tolerance = 1e-12) const;
    bool isUpperTriangular(double tolerance = 1e-12) const;
    bool isLowerTriangular(double tolerance = 1e-12) const;
    bool isInvertible(double tolerance = 1e-12) const;

    // 特征值和特征向量（实数域）
    // 若存在复特征值，eigenvalues()中对应项返回复根实部，eigenvectors()对应项返回零向量。
    std::array<double, 3> eigenvalues() const;
    std::array<Vector3D, 3> eigenvectors() const;

    // 矩阵范数
    double frobeniusNorm() const;
    double infinityNorm() const;
    double oneNorm() const;
    double spectralNorm() const; // 最大奇异值

    // 就地操作
    void transposeInPlace();
    void invertInPlace();

    // 静态工厂方法
    static Matrix3x3 identity();
    static Matrix3x3 zero();
    static Matrix3x3 diagonal(double d0, double d1, double d2);
    static Matrix3x3 diagonal(const Vector3D& diag);

    // 旋转矩阵
    static Matrix3x3 rotationX(double angle);
    static Matrix3x3 rotationY(double angle);
    static Matrix3x3 rotationZ(double angle);
    static Matrix3x3 rotationAroundAxis(const Vector3D& axis, double angle);
    static Matrix3x3 rotationFromEuler(double roll, double pitch, double yaw);

    // 缩放矩阵
    static Matrix3x3 scaling(double sx, double sy, double sz);
    static Matrix3x3 scaling(const Vector3D& scale);
    static Matrix3x3 uniformScaling(double scale);

    // 反射矩阵
    static Matrix3x3 reflection(const Vector3D& normal);

    // 从向量构造反对称矩阵（叉积矩阵）
    static Matrix3x3 crossProductMatrix(const Vector3D& vec);

    // 从两个向量构造旋转矩阵
    static Matrix3x3 rotationFromTo(const Vector3D& from, const Vector3D& to);

    // 剪切矩阵
    static Matrix3x3 shearXY(double shx, double shy);
    static Matrix3x3 shearXZ(double shx, double shz);
    static Matrix3x3 shearYZ(double shy, double shz);

    // 矩阵分解
    bool luDecomposition(Matrix3x3& L, Matrix3x3& U) const;
    bool qrDecomposition(Matrix3x3& Q, Matrix3x3& R) const;

    // 矩阵幂运算
    Matrix3x3 power(int n) const;
    Matrix3x3 sqrt() const; // 矩阵平方根

    // 矩阵指数和对数
    Matrix3x3 exp() const; // 矩阵指数
    Matrix3x3 log() const; // 矩阵对数

    // 条件数
    double conditionNumber() const;

    // 数据访问
    const std::array<std::array<double, 3>, 3>& data() const
    {
        return data_;
    }
    std::array<std::array<double, 3>, 3>& data()
    {
        return data_;
    }

    // 字符串表示
    std::string toString() const;

  private:
    static constexpr double EPSILON = 1e-12;
};

// 全局运算符
Matrix3x3 operator*(double scalar, const Matrix3x3& matrix);

// 输出流运算符
std::ostream& operator<<(std::ostream& os, const Matrix3x3& matrix);

} // namespace Geometry
} // namespace SCDAT

#endif // SCDAT_MATRIX3X3_H
