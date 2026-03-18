/**
 * @file Vector3D.h
 * @brief 三维向量数学运算头文件
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026年3月18日 7:39:41
 *
 * 定义三维向量类，提供完整的向量数学运算功能。
 * 继承自Point3D类，扩展向量特有的运算如叉积、投影、旋转等。
 *
 * @note 向量运算遵循右手定则，角度计算使用弧度制
 * @note 支持球坐标和柱坐标系转换
 *
 * 主要功能：
 * - 向量代数运算
 * - 几何变换
 * - 坐标系转换
 * - 投影和分解
 * - 旋转变换
 * @ingroup ToolsModule
 */

#ifndef SCDAT_VECTOR3D_H
#define SCDAT_VECTOR3D_H

#include "Point3D.h"
#include <array>
#include <string>

namespace SCDAT
{
namespace Geometry
{

/**
 * @brief 三维向量类
 *
 * 继承自Point3D，专门用于表示向量运算。
 * 提供向量特有的运算功能，如叉积、标量三重积等。
 */
class Vector3D : public Point3D
{
  public:
    // 构造函数
    constexpr Vector3D() noexcept : Point3D() {}
    constexpr Vector3D(double x, double y, double z) noexcept : Point3D(x, y, z) {}
    constexpr Vector3D(const Point3D& point) noexcept : Point3D(point) {}
    constexpr Vector3D(const std::array<double, 3>& coords) : Point3D(coords) {}
    constexpr Vector3D(std::array<double, 3>&& coords) noexcept : Point3D(std::move(coords)) {}

    // 从两点构造向量
    constexpr Vector3D(const Point3D& from, const Point3D& to) noexcept
        : Point3D(to.x() - from.x(), to.y() - from.y(), to.z() - from.z())
    {
    }

    // 模板化构造函数支持不同数值类型
    template <typename T>
    constexpr Vector3D(T x, T y, T z) noexcept
        : Point3D(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z))
    {
    }

    // 拷贝构造和赋值
    Vector3D(const Vector3D& other) = default;
    Vector3D& operator=(const Vector3D& other) = default;
    Vector3D(Vector3D&& other) noexcept = default;
    Vector3D& operator=(Vector3D&& other) noexcept = default;

    // 析构函数
    virtual ~Vector3D() = default;

    // 向量特有的运算

    /**
     * @brief 向量长度（模）
     */
    double length() const noexcept
    {
        return magnitude();
    }

    /**
     * @brief 向量长度的平方
     */
    constexpr double lengthSquared() const noexcept
    {
        return magnitudeSquared();
    }

    /**
     * @brief 单位向量
     */
    Vector3D unit() const
    {
        return Vector3D(normalized());
    }

    /**
     * @brief 就地单位化
     */
    void makeUnit()
    {
        normalize();
    }

    /**
     * @brief 向量夹角（弧度）
     */
    double angleTo(const Vector3D& other) const;

    /**
     * @brief 向量夹角（度）
     */
    double angleToInDegrees(const Vector3D& other) const;

    /**
     * @brief 标量三重积 [a, b, c] = a · (b × c)
     */
    double scalarTripleProduct(const Vector3D& b, const Vector3D& c) const;

    /**
     * @brief 向量三重积 a × (b × c)
     */
    Vector3D vectorTripleProduct(const Vector3D& b, const Vector3D& c) const;

    /**
     * @brief 检查是否与另一向量平行
     */
    bool isParallelTo(const Vector3D& other, double tolerance = 1e-6) const noexcept;

    /**
     * @brief 检查是否与另一向量垂直
     */
    bool isPerpendicularTo(const Vector3D& other, double tolerance = 1e-6) const noexcept;

    /**
     * @brief 检查是否为单位向量
     */
    bool isUnit(double tolerance = 1e-6) const noexcept;

    /**
     * @brief 向量在另一向量上的投影
     */
    Vector3D projectionOnto(const Vector3D& other) const;

    /**
     * @brief 向量在另一向量上的投影长度
     */
    double projectionLengthOnto(const Vector3D& other) const;

    /**
     * @brief 向量的垂直分量（相对于另一向量）
     */
    Vector3D perpendicularComponent(const Vector3D& other) const;

    /**
     * @brief 向量绕轴旋转
     */
    Vector3D rotateAroundAxis(const Vector3D& axis, double angle) const;

    /**
     * @brief 就地绕轴旋转
     */
    void rotateAroundAxisInPlace(const Vector3D& axis, double angle);

    /**
     * @brief 获取与当前向量垂直的任意向量
     */
    Vector3D getPerpendicularVector() const;

    /**
     * @brief 获取与当前向量和给定向量都垂直的向量
     */
    Vector3D getPerpendicularVector(const Vector3D& other) const;

    /**
     * @brief 安全归一化，零向量返回零向量
     */
    Vector3D safeNormalized() const noexcept;

    /**
     * @brief 安全就地归一化，零向量保持不变
     */
    void safeNormalize() noexcept;

    /**
     * @brief 计算相对于另一向量的反射
     */
    Vector3D reflect(const Vector3D& normal) const;

    /**
     * @brief 计算两点的中点（作为位置向量）
     */
    Vector3D midpoint(const Vector3D& other) const noexcept;

    /**
     * @brief 分量级乘法
     */
    Vector3D componentMultiply(const Vector3D& other) const noexcept;

    /**
     * @brief 分量级除法
     */
    Vector3D componentDivide(const Vector3D& other) const;

    // 静态工厂方法
    static Vector3D unitX() noexcept
    {
        return Vector3D(1.0, 0.0, 0.0);
    }
    static Vector3D unitY() noexcept
    {
        return Vector3D(0.0, 1.0, 0.0);
    }
    static Vector3D unitZ() noexcept
    {
        return Vector3D(0.0, 0.0, 1.0);
    }
    static Vector3D zero() noexcept
    {
        return Vector3D(0.0, 0.0, 0.0);
    }
    static Vector3D one() noexcept
    {
        return Vector3D(1.0, 1.0, 1.0);
    }

    /**
     * @brief 创建随机单位向量
     */
    static Vector3D randomUnit() noexcept;

    /**
     * @brief 从球坐标创建向量
     * @param radius 半径
     * @param theta 极角（与z轴夹角，弧度）
     * @param phi 方位角（在xy平面内与x轴夹角，弧度）
     */
    static Vector3D fromSpherical(double radius, double theta, double phi) noexcept;

    /**
     * @brief 从柱坐标创建向量
     * @param rho 径向距离
     * @param phi 方位角（弧度）
     * @param z z坐标
     */
    static Vector3D fromCylindrical(double rho, double phi, double z) noexcept;

    /**
     * @brief 转换为球坐标
     * @return {radius, theta, phi}
     */
    std::array<double, 3> toSpherical() const noexcept;

    /**
     * @brief 转换为柱坐标
     * @return {rho, phi, z}
     */
    std::array<double, 3> toCylindrical() const noexcept;

    /**
     * @brief 计算向量的方向余弦
     * @return {cos(alpha), cos(beta), cos(gamma)}
     */
    std::array<double, 3> directionCosines() const;

    /**
     * @brief 快速距离比较，避免开方运算
     */
    bool isCloserThan(const Vector3D& other, double distance) const noexcept;
    bool isFartherThan(const Vector3D& other, double distance) const noexcept;

    /**
     * @brief 计算曼哈顿距离
     */
    double manhattanDistance() const noexcept;
    double manhattanDistanceTo(const Vector3D& other) const noexcept;

    /**
     * @brief 计算切比雪夫距离
     */
    double chebyshevDistance() const noexcept;
    double chebyshevDistanceTo(const Vector3D& other) const noexcept;

    /**
     * @brief 计算向量的最大分量值
     */
    double maxComponent() const noexcept;

    /**
     * @brief 计算向量的最小分量值
     */
    double minComponent() const noexcept;

    /**
     * @brief 获取最大分量的轴索引
     */
    int maxAxis() const noexcept;

    /**
     * @brief 获取最小分量的轴索引
     */
    int minAxis() const noexcept;

    /**
     * @brief 计算向量的绝对值（各分量取绝对值）
     */
    Vector3D abs() const noexcept;

    /**
     * @brief 向量分量的符号函数
     */
    Vector3D sign() const noexcept;

    /**
     * @brief 向量分量的向下取整
     */
    Vector3D floor() const noexcept;

    /**
     * @brief 向量分量的向上取整
     */
    Vector3D ceil() const noexcept;

    /**
     * @brief 向量分量的四舍五入
     */
    Vector3D round() const noexcept;

    /**
     * @brief 向量的分量级最小值
     */
    Vector3D min(const Vector3D& other) const noexcept;

    /**
     * @brief 向量的分量级最大值
     */
    Vector3D max(const Vector3D& other) const noexcept;

    /**
     * @brief 向量的分量级钳制（限制在范围内）
     */
    Vector3D clamp(const Vector3D& min_vec, const Vector3D& max_vec) const noexcept;

    /**
     * @brief 向量的分量级钳制（限制在标量范围内）
     */
    Vector3D clamp(double min_val, double max_val) const noexcept;

    /**
     * @brief 计算向量的p范数
     */
    double pNorm(double p) const;

    /**
     * @brief 检查向量是否有限（无NaN或无穷大）
     */
    bool isFinite() const noexcept;

    /**
     * @brief 检查向量是否包含NaN
     */
    bool hasNaN() const noexcept;

    /**
     * @brief 检查向量是否包含无穷大
     */
    bool hasInf() const noexcept;

    // 重载运算符以返回Vector3D类型
    Vector3D operator+(const Vector3D& other) const noexcept
    {
        return Vector3D(Point3D::operator+(other));
    }

    Vector3D operator-(const Vector3D& other) const noexcept
    {
        return Vector3D(Point3D::operator-(other));
    }

    Vector3D operator*(double scalar) const noexcept
    {
        return Vector3D(Point3D::operator*(scalar));
    }

    Vector3D operator/(double scalar) const
    {
        return Vector3D(Point3D::operator/(scalar));
    }

    Vector3D operator-() const noexcept
    {
        return Vector3D(-x(), -y(), -z());
    }

    // 字符串表示
    std::string toString() const override;
};

// 全局运算符
inline Vector3D operator*(double scalar, const Vector3D& vector) noexcept
{
    return vector * scalar;
}

// 流输出运算符
inline std::ostream& operator<<(std::ostream& os, const Vector3D& vec)
{
    return os << vec.toString();
}

} // namespace Geometry
} // namespace SCDAT

#endif // SCDAT_VECTOR3D_H
