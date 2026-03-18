/**
 * @file Vector3D.cpp
 * @brief 三维向量数学运算实现
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026年3月18日 8:11:13
 *
 * 实现三维向量的基本数学运算、几何变换、坐标系转换和物理计算功能。
 * 支持向量代数、旋转变换、投影运算和坐标系转换等核心功能。
 *
 * @note 所有角度计算使用弧度制，坐标系遵循右手定则
 */

#include "../include/Vector3D.h"
#include <algorithm> // std::max需要
#include <cmath>
#include <cstdlib> // rand函数需要
#include <iomanip> // 格式化输出需要
#include <random>  // 随机数生成需要
#include <sstream>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace SCDAT
{
namespace Geometry
{

/**
 * @brief 计算两个向量之间的夹角（弧度制）
 *
 * 使用向量点积公式计算夹角：cos(θ) = (a·b) / (|a||b|)
 *
 * @param other 另一个向量
 * @return 夹角值（弧度，范围[0, π]）
 * @throws std::runtime_error 当任一向量为零向量时
 *
 * @note 返回值范围为[0, π]，符合数学定义
 */
double Vector3D::angleTo(const Vector3D& other) const
{
    double dot_product = this->dot(other);
    double mag_product = this->magnitude() * other.magnitude();

    if (mag_product < EPSILON)
    {
        throw std::runtime_error("Cannot compute angle with zero-length vector");
    }

    double cos_angle = dot_product / mag_product;
    // 确保cos_angle在[-1, 1]范围内，避免数值误差
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

    return std::acos(cos_angle);
}

/**
 * @brief 计算两个向量之间的夹角（角度制）
 *
 * @param other 另一个向量
 * @return 夹角值（角度，范围[0°, 180°]）
 * @throws std::runtime_error 当任一向量为零向量时
 */
double Vector3D::angleToInDegrees(const Vector3D& other) const
{
    return angleTo(other) * 180.0 / M_PI;
}

/**
 * @brief 计算标量三重积 a·(b×c)
 *
 * 标量三重积的几何意义是由三个向量构成的平行六面体的体积。
 *
 * @param b 第二个向量
 * @param c 第三个向量
 * @return 标量三重积值
 *
 * @note 标量三重积满足循环性质：a·(b×c) = b·(c×a) = c·(a×b)
 */
double Vector3D::scalarTripleProduct(const Vector3D& b, const Vector3D& c) const
{
    return this->dot(b.cross(c));
}

/**
 * @brief 计算向量三重积 a×(b×c)
 *
 * 向量三重积可以用公式展开：a×(b×c) = b(a·c) - c(a·b)
 *
 * @param b 第二个向量
 * @param c 第三个向量
 * @return 向量三重积结果
 */
Vector3D Vector3D::vectorTripleProduct(const Vector3D& b, const Vector3D& c) const
{
    return Vector3D(this->cross(b.cross(c)));
}

/**
 * @brief 检查两个向量是否平行
 *
 * 通过计算叉积的模长来判断平行性，平行向量的叉积为零向量。
 *
 * @param other 另一个向量
 * @param tolerance 容差值
 * @return true表示平行，false表示不平行
 */
bool Vector3D::isParallelTo(const Vector3D& other, double tolerance) const noexcept
{
    Vector3D cross_product = this->cross(other);
    return cross_product.magnitude() < tolerance;
}

/**
 * @brief 检查两个向量是否垂直
 *
 * 通过计算点积来判断垂直性，垂直向量的点积为零。
 *
 * @param other 另一个向量
 * @param tolerance 容差值
 * @return true表示垂直，false表示不垂直
 */
bool Vector3D::isPerpendicularTo(const Vector3D& other, double tolerance) const noexcept
{
    return std::abs(this->dot(other)) < tolerance;
}

/**
 * @brief 检查是否为单位向量
 *
 * 单位向量的模长为1，在数值计算中允许一定的容差。
 *
 * @param tolerance 容差值
 * @return true表示是单位向量，false表示不是
 */
bool Vector3D::isUnit(double tolerance) const noexcept
{
    return std::abs(magnitude() - 1.0) < tolerance;
}

/**
 * @brief 计算向量在另一向量上的投影
 *
 * 投影公式：proj_b(a) = (a·b / |b|²) * b
 *
 * @param other 投影目标向量
 * @return 投影向量
 * @throws std::runtime_error 当目标向量为零向量时
 *
 * @note 投影向量与目标向量平行
 */
Vector3D Vector3D::projectionOnto(const Vector3D& other) const
{
    double other_mag_sq = other.magnitudeSquared();
    if (other_mag_sq < EPSILON)
    {
        throw std::runtime_error("Cannot project onto zero-length vector");
    }

    double dot_product = this->dot(other);
    return Vector3D(other * (dot_product / other_mag_sq));
}

/**
 * @brief 计算向量在另一向量上的投影长度（标量投影）
 *
 * 标量投影公式：comp_b(a) = a·b / |b|
 *
 * @param other 投影目标向量
 * @return 投影长度（可为负值）
 * @throws std::runtime_error 当目标向量为零向量时
 */
double Vector3D::projectionLengthOnto(const Vector3D& other) const
{
    double other_mag = other.magnitude();
    if (other_mag < EPSILON)
    {
        throw std::runtime_error("Cannot project onto zero-length vector");
    }

    return this->dot(other) / other_mag;
}

/**
 * @brief 计算向量垂直于另一向量的分量
 *
 * 垂直分量 = 原向量 - 投影分量
 *
 * @param other 参考向量
 * @return 垂直分量向量
 *
 * @note 结果向量与参考向量垂直
 */
Vector3D Vector3D::perpendicularComponent(const Vector3D& other) const
{
    return Vector3D(*this - projectionOnto(other));
}

// 绕轴旋转（罗德里格斯公式）
Vector3D Vector3D::rotateAroundAxis(const Vector3D& axis, double angle) const
{
    Vector3D k = axis.unit();
    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);

    Vector3D v_rot =
        (*this) * cos_angle + k.cross(*this) * sin_angle + k * (k.dot(*this) * (1.0 - cos_angle));

    return v_rot;
}

void Vector3D::rotateAroundAxisInPlace(const Vector3D& axis, double angle)
{
    *this = rotateAroundAxis(axis, angle);
}

// 获取垂直向量
/** @note 在三维空间中，给定一个非零向量A，自动生成一个与A垂直（正交）的单位向量B，可以用A和B构造一个正交基，或者用于旋转、投影、几何变换等场景 
 */
Vector3D Vector3D::getPerpendicularVector() const
{
    // 选择与当前向量最不平行的坐标轴
    Vector3D candidate;
    if (std::abs(x()) < std::abs(y()) && std::abs(x()) < std::abs(z()))
    {
        candidate = Vector3D::unitX(); 
    }
    else if (std::abs(y()) < std::abs(z()))
    {
        candidate = Vector3D::unitY();
    }
    else
    {
        candidate = Vector3D::unitZ();
    }

    return Vector3D(this->cross(candidate).normalized());
}

Vector3D Vector3D::getPerpendicularVector(const Vector3D& other) const
{
    return Vector3D(this->cross(other).normalized());
}

/** @brief 从球面坐标系转换到笛卡尔坐标系
 *
 * @param radius 半径
 * @param theta 极角（与z轴的夹角，范围[0, π]）
 * @param phi 方位角（在xy平面上的投影与x轴的夹角，范围[0, 2π)）
 * @return 笛卡尔坐标系下的向量
 */
Vector3D Vector3D::fromSpherical(double radius, double theta, double phi) noexcept
{
    double sin_theta = std::sin(theta);
    return Vector3D(radius * sin_theta * std::cos(phi), radius * sin_theta * std::sin(phi),
                    radius * std::cos(theta));
}

/** @brief 从柱面坐标系转换到笛卡尔坐标系
 *
 * @param rho 半径
 * @param phi 角度
 * @param z z坐标
 * @return 笛卡尔坐标系下的向量
 */
Vector3D Vector3D::fromCylindrical(double rho, double phi, double z) noexcept
{
    return Vector3D(rho * std::cos(phi), rho * std::sin(phi), z);
}
/** @brief 将向量从笛卡尔坐标系转换到球面坐标系
 *
 * @return 包含半径、极角和方位角的数组
 */
std::array<double, 3> Vector3D::toSpherical() const noexcept
{
    double radius = magnitude();
    if (radius < EPSILON)
    {
        return {0.0, 0.0, 0.0};
    }

    double theta = std::acos(z() / radius);
    double phi = std::atan2(y(), x());

    return {radius, theta, phi};
}
/** @brief 将向量从笛卡尔坐标系转换到柱面坐标系
 *
 * @return 包含半径、角度和z坐标的数组
 */
std::array<double, 3> Vector3D::toCylindrical() const noexcept
{
    double rho = std::sqrt(x() * x() + y() * y());
    double phi = std::atan2(y(), x());

    return {rho, phi, z()};
}

/**
 * @brief 安全归一化，零向量返回零向量
 */
Vector3D Vector3D::safeNormalized() const noexcept
{
    double mag = magnitude();
    if (mag < EPSILON)
    {
        return Vector3D::zero();
    }
    return Vector3D(*this / mag);
}

/**
 * @brief 安全就地归一化，零向量保持不变
 */
void Vector3D::safeNormalize() noexcept
{
    double mag = magnitude();
    if (mag >= EPSILON)
    {
        *this /= mag;
    }
}

/**
 * @brief 计算相对于另一向量的反射
 */
Vector3D Vector3D::reflect(const Vector3D& normal) const
{
    Vector3D normalized_normal = normal.unit();
    return Vector3D(*this - normalized_normal * (2.0 * this->dot(normalized_normal)));
}

/**
 * @brief 计算两点的中点（作为位置向量）
 */
Vector3D Vector3D::midpoint(const Vector3D& other) const noexcept
{
    return Vector3D((x() + other.x()) * 0.5, (y() + other.y()) * 0.5, (z() + other.z()) * 0.5);
}

/**
 * @brief 分量级乘法
 */
Vector3D Vector3D::componentMultiply(const Vector3D& other) const noexcept
{
    return Vector3D(x() * other.x(), y() * other.y(), z() * other.z());
}

/**
 * @brief 分量级除法
 */
Vector3D Vector3D::componentDivide(const Vector3D& other) const
{
    if (std::abs(other.x()) < EPSILON || std::abs(other.y()) < EPSILON ||
        std::abs(other.z()) < EPSILON)
    {
        throw std::runtime_error("Division by zero component in Vector3D");
    }
    return Vector3D(x() / other.x(), y() / other.y(), z() / other.z());
}

/**
 * @brief 创建随机单位向量
 */
Vector3D Vector3D::randomUnit() noexcept
{
    // 使用改进的Marsaglia方法生成均匀分布的单位向量
    static thread_local std::random_device rd;
    static thread_local std::mt19937 gen(rd());
    static thread_local std::uniform_real_distribution<double> dis(-1.0, 1.0);

    double x1, x2, w;
    do
    {
        x1 = dis(gen);
        x2 = dis(gen);
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0 || w == 0.0);

    double w_sqrt = std::sqrt(1.0 - w);
    return Vector3D(2.0 * x1 * w_sqrt, 2.0 * x2 * w_sqrt, 1.0 - 2.0 * w);
}

/**
 * @brief 计算向量的方向余弦
 */
std::array<double, 3> Vector3D::directionCosines() const
{
    double mag = magnitude();
    if (mag < EPSILON)
    {
        throw std::runtime_error("Cannot compute direction cosines for zero vector");
    }

    return {x() / mag, y() / mag, z() / mag};
}

/**
 * @brief 快速距离比较，避免开方运算
 */
bool Vector3D::isCloserThan(const Vector3D& other, double distance) const noexcept
{
    return distanceSquaredTo(other) < (distance * distance);
}

bool Vector3D::isFartherThan(const Vector3D& other, double distance) const noexcept
{
    return distanceSquaredTo(other) > (distance * distance);
}

/**
 * @brief 计算曼哈顿距离
 */
double Vector3D::manhattanDistance() const noexcept
{
    return std::abs(x()) + std::abs(y()) + std::abs(z());
}

double Vector3D::manhattanDistanceTo(const Vector3D& other) const noexcept
{
    return std::abs(x() - other.x()) + std::abs(y() - other.y()) + std::abs(z() - other.z());
}

/**
 * @brief 计算切比雪夫距离
 */
double Vector3D::chebyshevDistance() const noexcept
{
    return std::max({std::abs(x()), std::abs(y()), std::abs(z())});
}

double Vector3D::chebyshevDistanceTo(const Vector3D& other) const noexcept
{
    return std::max(
        {std::abs(x() - other.x()), std::abs(y() - other.y()), std::abs(z() - other.z())});
}

/**
 * @brief 计算向量的最大分量值
 */
double Vector3D::maxComponent() const noexcept
{
    return std::max({x(), y(), z()});
}

/**
 * @brief 计算向量的最小分量值
 */
double Vector3D::minComponent() const noexcept
{
    return std::min({x(), y(), z()});
}

/**
 * @brief 获取最大分量的轴索引
 */
int Vector3D::maxAxis() const noexcept
{
    if (std::abs(x()) >= std::abs(y()) && std::abs(x()) >= std::abs(z()))
        return 0;
    if (std::abs(y()) >= std::abs(z()))
        return 1;
    return 2;
}

/**
 * @brief 获取最小分量的轴索引
 */
int Vector3D::minAxis() const noexcept
{
    if (std::abs(x()) <= std::abs(y()) && std::abs(x()) <= std::abs(z()))
        return 0;
    if (std::abs(y()) <= std::abs(z()))
        return 1;
    return 2;
}

/**
 * @brief 计算向量的绝对值（各分量取绝对值）
 */
Vector3D Vector3D::abs() const noexcept
{
    return Vector3D(std::abs(x()), std::abs(y()), std::abs(z()));
}

/**
 * @brief 向量分量的符号函数
 */
Vector3D Vector3D::sign() const noexcept
{
    auto sign_func = [](double val) -> double
    {
        if (val > 0.0)
            return 1.0;
        if (val < 0.0)
            return -1.0;
        return 0.0;
    };
    return Vector3D(sign_func(x()), sign_func(y()), sign_func(z()));
}

/**
 * @brief 向量分量的向下取整
 */
Vector3D Vector3D::floor() const noexcept
{
    return Vector3D(std::floor(x()), std::floor(y()), std::floor(z()));
}

/**
 * @brief 向量分量的向上取整
 */
Vector3D Vector3D::ceil() const noexcept
{
    return Vector3D(std::ceil(x()), std::ceil(y()), std::ceil(z()));
}

/**
 * @brief 向量分量的四舍五入
 */
Vector3D Vector3D::round() const noexcept
{
    return Vector3D(std::round(x()), std::round(y()), std::round(z()));
}

/**
 * @brief 向量的分量级最小值
 */
Vector3D Vector3D::min(const Vector3D& other) const noexcept
{
    return Vector3D(std::min(x(), other.x()), std::min(y(), other.y()), std::min(z(), other.z()));
}

/**
 * @brief 向量的分量级最大值
 */
Vector3D Vector3D::max(const Vector3D& other) const noexcept
{
    return Vector3D(std::max(x(), other.x()), std::max(y(), other.y()), std::max(z(), other.z()));
}

/**
 * @brief 向量的分量级钳制（限制在范围内）
 */
Vector3D Vector3D::clamp(const Vector3D& min_vec, const Vector3D& max_vec) const noexcept
{
    return Vector3D(std::max(min_vec.x(), std::min(max_vec.x(), x())),
                    std::max(min_vec.y(), std::min(max_vec.y(), y())),
                    std::max(min_vec.z(), std::min(max_vec.z(), z())));
}

/**
 * @brief 向量的分量级钳制（限制在标量范围内）
 */
Vector3D Vector3D::clamp(double min_val, double max_val) const noexcept
{
    return Vector3D(std::max(min_val, std::min(max_val, x())),
                    std::max(min_val, std::min(max_val, y())),
                    std::max(min_val, std::min(max_val, z())));
}

/**
 * @brief 计算向量的p范数
 */
double Vector3D::pNorm(double p) const
{
    if (p <= 0.0)
    {
        throw std::invalid_argument("p must be positive for p-norm calculation");
    }

    if (std::isinf(p))
    {
        return chebyshevDistance(); // 无穷范数
    }

    return std::pow(std::pow(std::abs(x()), p) + std::pow(std::abs(y()), p) +
                        std::pow(std::abs(z()), p),
                    1.0 / p);
}

/**
 * @brief 检查向量是否有限（无NaN或无穷大）
 */
bool Vector3D::isFinite() const noexcept
{
    return std::isfinite(x()) && std::isfinite(y()) && std::isfinite(z());
}

/**
 * @brief 检查向量是否包含NaN
 */
bool Vector3D::hasNaN() const noexcept
{
    return std::isnan(x()) || std::isnan(y()) || std::isnan(z());
}

/**
 * @brief 检查向量是否包含无穷大
 */
bool Vector3D::hasInf() const noexcept
{
    return std::isinf(x()) || std::isinf(y()) || std::isinf(z());
}

// 字符串表示
std::string Vector3D::toString() const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6); // 控制精度
    oss << "Vector3D(" << x() << ", " << y() << ", " << z() << ")";
    return oss.str();
}

} // namespace Geometry
} // namespace SCDAT
