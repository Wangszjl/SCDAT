/**
 * @file Point3D.h
 * @brief 点定义头文件
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026年3月17日 16:30:25
 * @ingroup ToolsModule
 */

#ifndef SCDAT_POINT3D_H
#define SCDAT_POINT3D_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace SCDAT
{
namespace Geometry
{

/**
 * @brief 坐标轴枚举，提供类型安全的轴访问
 */
enum class Axis : int
{
    X = 0,
    Y = 1,
    Z = 2
};

/**
 * @brief 三维点类
 *
 * 表示三维空间中的一个点，提供基本的几何运算功能。
 * 设计目标：高性能、类型安全、易于使用。
 */
class Point3D
{
  public:
    // 构造函数
    constexpr Point3D() : x_(0.0), y_(0.0), z_(0.0) {}
    constexpr Point3D(double x, double y, double z) : x_(x), y_(y), z_(z) {}
    constexpr Point3D(const std::array<double, 3>& coords)
        : x_(coords[0]), y_(coords[1]), z_(coords[2])
    {
    }
    constexpr Point3D(std::array<double, 3>&& coords) noexcept
        : x_(coords[0]), y_(coords[1]), z_(coords[2])
    {
    }

    // 模板化构造函数支持不同数值类型
    template <typename T>
    constexpr Point3D(T x, T y, T z) noexcept
        : x_(static_cast<double>(x)), y_(static_cast<double>(y)), z_(static_cast<double>(z))
    {
    }

    // 拷贝构造和赋值
    Point3D(const Point3D& other) = default;
    Point3D& operator=(const Point3D& other) = default;
    Point3D(Point3D&& other) noexcept = default;
    Point3D& operator=(Point3D&& other) noexcept = default;

    // 析构函数
    ~Point3D() = default;

    // 访问器
    constexpr double x() const
    {
        return x_;
    }
    constexpr double y() const
    {
        return y_;
    }
    constexpr double z() const
    {
        return z_;
    }

    double& x()
    {
        return x_;
    }
    double& y()
    {
        return y_;
    }
    double& z()
    {
        return z_;
    }

    // 下标访问
    double operator[](int index) const
    {
        if (index < 0 || index > 2)
        {
            throw std::out_of_range("Point3D index out of range: " + std::to_string(index));
        }
        return (&x_)[index]; // 更高效的实现
    }

    double& operator[](int index)
    {
        if (index < 0 || index > 2)
        {
            throw std::out_of_range("Point3D index out of range: " + std::to_string(index));
        }
        return (&x_)[index];
    }

    // 类型安全的轴访问
    constexpr double operator[](Axis axis) const noexcept
    {
        return (&x_)[static_cast<int>(axis)];
    }

    double& operator[](Axis axis) noexcept
    {
        return (&x_)[static_cast<int>(axis)];
    }

    // 设置器
    void set(double x, double y, double z) noexcept
    {
        x_ = x;
        y_ = y;
        z_ = z;
    }

    void setX(double x) noexcept
    {
        x_ = x;
    }
    void setY(double y) noexcept
    {
        y_ = y;
    }
    void setZ(double z) noexcept
    {
        z_ = z;
    }

    // 算术运算符
    constexpr Point3D operator+(const Point3D& other) const noexcept
    {
        return Point3D(x_ + other.x_, y_ + other.y_, z_ + other.z_);
    }

    constexpr Point3D operator-(const Point3D& other) const noexcept
    {
        return Point3D(x_ - other.x_, y_ - other.y_, z_ - other.z_);
    }

    constexpr Point3D operator*(double scalar) const noexcept
    {
        return Point3D(x_ * scalar, y_ * scalar, z_ * scalar);
    }

    Point3D operator/(double scalar) const
    {
        if (std::abs(scalar) < EPSILON)
        {
            throw std::runtime_error("Division by zero in Point3D");
        }
        return Point3D(x_ / scalar, y_ / scalar, z_ / scalar);
    }

    // 复合赋值运算符
    Point3D& operator+=(const Point3D& other) noexcept
    {
        x_ += other.x_;
        y_ += other.y_;
        z_ += other.z_;
        return *this;
    }

    Point3D& operator-=(const Point3D& other) noexcept
    {
        x_ -= other.x_;
        y_ -= other.y_;
        z_ -= other.z_;
        return *this;
    }

    Point3D& operator*=(double scalar) noexcept
    {
        x_ *= scalar;
        y_ *= scalar;
        z_ *= scalar;
        return *this;
    }

    Point3D& operator/=(double scalar)
    {
        if (std::abs(scalar) < EPSILON)
        {
            throw std::runtime_error("Division by zero in Point3D");
        }
        x_ /= scalar;
        y_ /= scalar;
        z_ /= scalar;
        return *this;
    }

    // 比较运算符
    bool operator==(const Point3D& other) const noexcept
    {
        return std::abs(x_ - other.x_) < EPSILON && std::abs(y_ - other.y_) < EPSILON &&
               std::abs(z_ - other.z_) < EPSILON;
    }

    bool operator!=(const Point3D& other) const noexcept
    {
        return !(*this == other);
    }

    // 添加三路比较支持 (C++20)
#if __cplusplus >= 202002L
    auto operator<=>(const Point3D& other) const noexcept
    {
        if (auto cmp = x_ <=> other.x_; cmp != 0)
            return cmp;
        if (auto cmp = y_ <=> other.y_; cmp != 0)
            return cmp;
        return z_ <=> other.z_;
    }
#endif

    // 几何运算
    /**
     * @brief 计算到另一点的距离
     */
    double distanceTo(const Point3D& other) const noexcept
    {
        double dx = x_ - other.x_;
        double dy = y_ - other.y_;
        double dz = z_ - other.z_;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    /**
     * @brief 快速距离比较，避免开方运算
     * @param other 另一个点
     * @param distance 比较距离
     * @return 是否距离小于指定值
     */
    bool isCloserThan(const Point3D& other, double distance) const noexcept
    {
        return distanceSquaredTo(other) < (distance * distance);
    }

    /**
     * @brief 快速距离比较，避免开方运算
     * @param other 另一个点
     * @param distance 比较距离
     * @return 是否距离大于指定值
     */
    bool isFartherThan(const Point3D& other, double distance) const noexcept
    {
        return distanceSquaredTo(other) > (distance * distance);
    }

    /**
     * @brief 计算到另一点的距离的平方（避免开方运算）
     */
    constexpr double distanceSquaredTo(const Point3D& other) const noexcept
    {
        double dx = x_ - other.x_;
        double dy = y_ - other.y_;
        double dz = z_ - other.z_;
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * @brief 计算到原点的距离
     */
    double magnitude() const noexcept
    {
        return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    }

    /**
     * @brief 计算到原点的距离的平方
     */
    constexpr double magnitudeSquared() const noexcept
    {
        return x_ * x_ + y_ * y_ + z_ * z_;
    }

    /**
     * @brief 点积（Vector3D需要此功能）
     */
    constexpr double dot(const Point3D& other) const noexcept
    {
        return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
    }

    /**
     * @brief 叉积（Vector3D需要此功能）
     */
    constexpr Point3D cross(const Point3D& other) const noexcept
    {
        return Point3D(y_ * other.z_ - z_ * other.y_, z_ * other.x_ - x_ * other.z_,
                       x_ * other.y_ - y_ * other.x_);
    }

    /**
     * @brief 线性插值
     */
    constexpr Point3D lerp(const Point3D& other, double t) const noexcept
    {
        return Point3D(x_ + t * (other.x_ - x_), y_ + t * (other.y_ - y_),
                       z_ + t * (other.z_ - z_));
    }

    /**
     * @brief 归一化（Vector3D需要此功能）
     */
    Point3D normalized() const
    {
        double mag = magnitude();
        if (mag < EPSILON)
        {
            throw std::runtime_error("Cannot normalize zero-length vector");
        }
        return *this / mag;
    }

    /**
     * @brief 就地归一化（Vector3D需要此功能）
     */
    void normalize()
    {
        double mag = magnitude();
        if (mag < EPSILON)
        {
            throw std::runtime_error("Cannot normalize zero-length vector");
        }
        *this /= mag;
    }

    /**
     * @brief 计算两点的中点
     */
    constexpr Point3D midpoint(const Point3D& other) const noexcept
    {
        return Point3D((x_ + other.x_) * 0.5, (y_ + other.y_) * 0.5, (z_ + other.z_) * 0.5);
    }

    /**
     * @brief 计算到原点的曼哈顿距离
     */
    inline double manhattanDistance() const noexcept
    {
        return std::abs(x_) + std::abs(y_) + std::abs(z_);
    }

    /**
     * @brief 计算到另一点的曼哈顿距离
     */
    inline double manhattanDistanceTo(const Point3D& other) const noexcept
    {
        return std::abs(x_ - other.x_) + std::abs(y_ - other.y_) + std::abs(z_ - other.z_);
    }

    /**
     * @brief 计算到原点的切比雪夫距离（无穷范数）
     */
    inline double chebyshevDistance() const noexcept
    {
        return std::max({std::abs(x_), std::abs(y_), std::abs(z_)});
    }

    /**
     * @brief 计算到另一点的切比雪夫距离
     */
    inline double chebyshevDistanceTo(const Point3D& other) const noexcept
    {
        return std::max(
            {std::abs(x_ - other.x_), std::abs(y_ - other.y_), std::abs(z_ - other.z_)});
    }

    /**
     * @brief 转换为数组
     */
    constexpr std::array<double, 3> toArray() const noexcept
    {
        return {x_, y_, z_};
    }

    /**
     * @brief 检查是否为零向量
     */
    bool isZero(double epsilon = 1e-12) const noexcept
    {
        return std::abs(x_) < epsilon && std::abs(y_) < epsilon && std::abs(z_) < epsilon;
    }

    /**
     * @brief 检查是否包含无效值（NaN或无穷大）
     */
    bool isValid() const noexcept
    {
        return std::isfinite(x_) && std::isfinite(y_) && std::isfinite(z_);
    }

    /**
     * @brief 获取最大分量的值
     */
    constexpr double maxComponent() const noexcept
    {
        return std::max({x_, y_, z_});
    }

    /**
     * @brief 获取最小分量的值
     */
    constexpr double minComponent() const noexcept
    {
        return std::min({x_, y_, z_});
    }

    /**
     * @brief 获取最大分量的轴
     */
    Axis maxAxis() const noexcept
    {
        if (std::abs(x_) >= std::abs(y_) && std::abs(x_) >= std::abs(z_))
            return Axis::X;
        if (std::abs(y_) >= std::abs(z_))
            return Axis::Y;
        return Axis::Z;
    }

    /**
     * @brief 获取最小分量的轴
     */
    Axis minAxis() const noexcept
    {
        if (std::abs(x_) <= std::abs(y_) && std::abs(x_) <= std::abs(z_))
            return Axis::X;
        if (std::abs(y_) <= std::abs(z_))
            return Axis::Y;
        return Axis::Z;
    }

    /**
     * @brief 分量级乘法
     */
    constexpr Point3D componentMultiply(const Point3D& other) const noexcept
    {
        return Point3D(x_ * other.x_, y_ * other.y_, z_ * other.z_);
    }

    /**
     * @brief 分量级除法
     */
    Point3D componentDivide(const Point3D& other) const
    {
        if (std::abs(other.x_) < EPSILON || std::abs(other.y_) < EPSILON ||
            std::abs(other.z_) < EPSILON)
        {
            throw std::runtime_error("Division by zero component in Point3D");
        }
        return Point3D(x_ / other.x_, y_ / other.y_, z_ / other.z_);
    }

    /**
     * @brief 哈希函数支持
     */
    std::size_t hash() const noexcept;

    /**
     * @brief 字符串表示
     */
    virtual std::string toString() const;

    // 静态工厂方法
    static constexpr Point3D zero()
    {
        return Point3D(0.0, 0.0, 0.0);
    }
    static constexpr Point3D unitX()
    {
        return Point3D(1.0, 0.0, 0.0);
    }
    static constexpr Point3D unitY()
    {
        return Point3D(0.0, 1.0, 0.0);
    }
    static constexpr Point3D unitZ()
    {
        return Point3D(0.0, 0.0, 1.0);
    }
    static constexpr Point3D one()
    {
        return Point3D(1.0, 1.0, 1.0);
    }

  private:
    double x_, y_, z_;

  public:
    static const double EPSILON;
};

// 全局运算符
inline constexpr Point3D operator*(double scalar, const Point3D& point) noexcept
{
    return point * scalar;
}

// 一元负号运算符
inline constexpr Point3D operator-(const Point3D& point) noexcept
{
    return Point3D(-point.x(), -point.y(), -point.z());
}

// 输出流运算符
inline std::ostream& operator<<(std::ostream& os, const Point3D& point)
{
    os << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")";
    return os;
}

// 输入流运算符
inline std::istream& operator>>(std::istream& is, Point3D& point)
{
    double x, y, z;
    is >> x >> y >> z;
    point.set(x, y, z);
    return is;
}

// 注意：Vector3D是独立的类，继承自Point3D

} // namespace Geometry
} // namespace SCDAT

// 添加std::hash特化支持
namespace std
{
template <> struct hash<SCDAT::Geometry::Point3D>
{
    std::size_t operator()(const SCDAT::Geometry::Point3D& point) const noexcept
    {
        return point.hash();
    }
};
} // namespace std

#endif // SCDAT_POINT3D_H
