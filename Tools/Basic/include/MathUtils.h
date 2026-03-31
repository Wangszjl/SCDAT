/**
 * @file MathUtils.h
 * @brief 数学工具函数头文件
 * @details 提供各种数学计算工具函数，包括插值、积分、微分等
 *
 * @author Wang Sizhan
 * @date 2026年3月18日 15:38:55
 * @version V0.0.1
 * @ingroup ToolsModule
 */

#ifndef SCDAT_BASIC_MATH_UTILS_H
#define SCDAT_BASIC_MATH_UTILS_H

#include "../Geometry/include/Matrix3x3.h"
#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
#include "Constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <vector>

namespace SCDAT
{
namespace Basic
{
namespace MathUtils
{

using Vector3D = SCDAT::Geometry::Vector3D;
using Point3D = SCDAT::Geometry::Point3D;
using Matrix3x3 = SCDAT::Geometry::Matrix3x3;

/**
 * @brief 插值方法枚举
 */
enum class InterpolationMethod
{
    NEAREST = 0,
    LINEAR = 1,
    QUADRATIC = 2,
    CUBIC = 3,
    SPLINE = 4
};

/**
 * @brief 积分方法枚举
 */
enum class IntegrationMethod
{
    RECTANGLE = 0,
    TRAPEZOIDAL = 1,
    SIMPSON = 2,
    GAUSS_LEGENDRE = 3,
    ADAPTIVE = 4
};

// ============================================================================
// 基础数学函数
// ============================================================================

/**
 * @brief 安全除法，避免除零错误
 * @param numerator 分子
 * @param denominator 分母
 * @param default_value 分母为零时的默认值
 * @return 除法结果
 */
template <typename T> inline T safeDivide(T numerator, T denominator, T default_value = T{0})
{
    return (std::abs(denominator) > SCDAT::Basic::Constants::PhysicsConstants::ZeroTolerance)
               ? (numerator / denominator)
               : default_value;
}

/**
 * @brief 计算两个数的平方和的平方根
 * @param a 第一个数
 * @param b 第二个数
 * @return sqrt(a² + b²)
 */
template <typename T> inline T hypot(T a, T b)
{
    return std::sqrt(a * a + b * b);
}

/**
 * @brief 计算三个数的平方和的平方根
 * @param a 第一个数
 * @param b 第二个数
 * @param c 第三个数
 * @return sqrt(a² + b² + c²)
 */
template <typename T> inline T hypot(T a, T b, T c)
{
    return std::sqrt(a * a + b * b + c * c);
}

/**
 * @brief 限制数值在指定范围内
 * @param value 输入值
 * @param min_val 最小值
 * @param max_val 最大值
 * @return 限制后的值
 */
template <typename T> inline T clamp(T value, T min_val, T max_val)
{
    return std::max(min_val, std::min(value, max_val));
}

/**
 * @brief 线性插值
 * @param t 插值参数 [0, 1]
 * @param a 起始值
 * @param b 结束值
 * @return 插值结果
 */
template <typename T> inline T lerp(double t, const T& a, const T& b)
{
    return a + t * (b - a);
}

/**
 * @brief 平滑步函数 (smoothstep)
 * @param t 输入参数
 * @return 平滑插值结果
 */
inline double smoothstep(double t)
{
    t = std::clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

/**
 * @brief 更平滑的步函数 (smootherstep)
 * @param t 输入参数
 * @return 更平滑的插值结果
 */
inline double smootherstep(double t)
{
    t = std::clamp(t, 0.0, 1.0);
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

// ============================================================================
// 向量和矩阵运算
// ============================================================================

/**
 * @brief 计算向量的单位向量
 * @param v 输入向量
 * @return 单位向量
 */
Vector3D normalize(const Vector3D& v);

/**
 * @brief 计算两个向量的夹角
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 夹角（弧度）
 */
double angleBetween(const Vector3D& a, const Vector3D& b);

/**
 * @brief 计算向量在另一个向量上的投影
 * @param v 被投影向量
 * @param onto 投影目标向量
 * @return 投影向量
 */
Vector3D project(const Vector3D& v, const Vector3D& onto);

/**
 * @brief 计算向量的反射
 * @param incident 入射向量
 * @param normal 法向量
 * @return 反射向量
 */
Vector3D reflect(const Vector3D& incident, const Vector3D& normal);

/**
 * @brief 计算向量的折射
 * @param incident 入射向量
 * @param normal 法向量
 * @param eta 折射率比
 * @return 折射向量
 */
Vector3D refract(const Vector3D& incident, const Vector3D& normal, double eta);

/**
 * @brief 生成正交基
 * @param normal 法向量
 * @param tangent1 第一个切向量（输出）
 * @param tangent2 第二个切向量（输出）
 */
void generateOrthonormalBasis(const Vector3D& normal, Vector3D& tangent1, Vector3D& tangent2);

// ============================================================================
// 插值函数
// ============================================================================

/**
 * @brief 一维线性插值
 * @param x 插值点
 * @param x_data x坐标数组
 * @param y_data y坐标数组
 * @return 插值结果
 */
double interpolate1D(double x, const std::vector<double>& x_data,
                     const std::vector<double>& y_data);

/**
 * @brief 二维双线性插值
 * @param x x坐标
 * @param y y坐标
 * @param x_data x坐标数组
 * @param y_data y坐标数组
 * @param z_data z值二维数组
 * @return 插值结果
 */
double interpolate2D(double x, double y, const std::vector<double>& x_data,
                     const std::vector<double>& y_data,
                     const std::vector<std::vector<double>>& z_data);

/**
 * @brief 三维三线性插值
 * @param point 插值点
 * @param grid_origin 网格原点
 * @param grid_spacing 网格间距
 * @param values 三维值数组
 * @return 插值结果
 */
double interpolate3D(const Point3D& point, const Point3D& grid_origin, const Vector3D& grid_spacing,
                     const std::vector<std::vector<std::vector<double>>>& values);

/**
 * @brief 样条插值
 * @param x 插值点
 * @param x_data x坐标数组
 * @param y_data y坐标数组
 * @return 插值结果
 */
double splineInterpolate(double x, const std::vector<double>& x_data,
                         const std::vector<double>& y_data);

// ============================================================================
// 数值积分
// ============================================================================

/**
 * @brief 一维数值积分
 * @param func 被积函数
 * @param a 积分下限
 * @param b 积分上限
 * @param method 积分方法
 * @param n 分割数
 * @return 积分结果
 */
double integrate1D(std::function<double(double)> func, double a, double b,
                   IntegrationMethod method = IntegrationMethod::SIMPSON, int n = 1000);

/**
 * @brief 高斯-勒让德积分
 * @param func 被积函数
 * @param a 积分下限
 * @param b 积分上限
 * @param n 积分点数
 * @return 积分结果
 */
double gaussLegendreIntegrate(std::function<double(double)> func, double a, double b, int n = 5);

/**
 * @brief 自适应积分
 * @param func 被积函数
 * @param a 积分下限
 * @param b 积分上限
 * @param tolerance 容差
 * @return 积分结果
 */
double adaptiveIntegrate(std::function<double(double)> func, double a, double b,
                         double tolerance = 1e-8);

// ============================================================================
// 数值微分
// ============================================================================

/**
 * @brief 数值求导（中心差分）
 * @param func 函数
 * @param x 求导点
 * @param h 步长
 * @return 导数值
 */
double derivative(std::function<double(double)> func, double x, double h = 1e-6);

/**
 * @brief 数值求二阶导数
 * @param func 函数
 * @param x 求导点
 * @param h 步长
 * @return 二阶导数值
 */
double secondDerivative(std::function<double(double)> func, double x, double h = 1e-6);

/**
 * @brief 计算梯度（数值方法）
 * @param func 三维标量函数
 * @param point 计算点
 * @param h 步长
 * @return 梯度向量
 */
Vector3D gradient(std::function<double(const Point3D&)> func, const Point3D& point,
                  double h = 1e-6);

// ============================================================================
// 统计函数
// ============================================================================

/**
 * @brief 计算平均值
 * @param data 数据数组
 * @return 平均值
 */
template <typename T> double mean(const std::vector<T>& data)
{
    if (data.empty())
        return 0.0;
    double sum = 0.0;
    for (const auto& value : data)
    {
        sum += static_cast<double>(value);
    }
    return sum / data.size();
}

/**
 * @brief 计算标准差
 * @param data 数据数组
 * @return 标准差
 */
template <typename T> double standardDeviation(const std::vector<T>& data)
{
    if (data.size() < 2)
        return 0.0;
    double avg = mean(data);
    double sum_sq_diff = 0.0;
    for (const auto& value : data)
    {
        double diff = static_cast<double>(value) - avg;
        sum_sq_diff += diff * diff;
    }
    return std::sqrt(sum_sq_diff / (data.size() - 1));
}

/**
 * @brief 计算方差
 * @param data 数据数组
 * @return 方差
 */
template <typename T> double variance(const std::vector<T>& data)
{
    double std_dev = standardDeviation(data);
    return std_dev * std_dev;
}

// ============================================================================
// 随机数生成
// ============================================================================

/**
 * @brief 生成均匀分布随机数
 * @param min 最小值
 * @param max 最大值
 * @return 随机数
 */
double uniformRandom(double min = 0.0, double max = 1.0);

/**
 * @brief 生成正态分布随机数
 * @param mean 均值
 * @param stddev 标准差
 * @return 随机数
 */
double normalRandom(double mean = 0.0, double stddev = 1.0);

/**
 * @brief 生成随机单位向量
 * @return 随机单位向量
 */
Vector3D randomUnitVector();

/**
 * @brief 生成球面上的随机点
 * @param radius 球半径
 * @return 随机点
 */
Point3D randomPointOnSphere(double radius = 1.0);

// ============================================================================
// 几何计算
// ============================================================================

/**
 * @brief 计算点到直线的距离
 * @param point 点
 * @param line_start 直线起点
 * @param line_end 直线终点
 * @return 距离
 */
double pointToLineDistance(const Point3D& point, const Point3D& line_start,
                           const Point3D& line_end);

/**
 * @brief 计算点到平面的距离
 * @param point 点
 * @param plane_point 平面上一点
 * @param plane_normal 平面法向量
 * @return 距离（带符号）
 */
double pointToPlaneDistance(const Point3D& point, const Point3D& plane_point,
                            const Vector3D& plane_normal);

/**
 * @brief 计算三角形面积
 * @param a 顶点A
 * @param b 顶点B
 * @param c 顶点C
 * @return 面积
 */
double triangleArea(const Point3D& a, const Point3D& b, const Point3D& c);

/**
 * @brief 计算四面体体积
 * @param a 顶点A
 * @param b 顶点B
 * @param c 顶点C
 * @param d 顶点D
 * @return 体积
 */
double tetrahedronVolume(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

/**
 * @brief 计算重心坐标
 * @param point 目标点
 * @param a 顶点A
 * @param b 顶点B
 * @param c 顶点C
 * @return 重心坐标 [u, v, w]，满足 u + v + w = 1
 */
std::array<double, 3> barycentricCoordinates(const Point3D& point, const Point3D& a,
                                             const Point3D& b, const Point3D& c);

// ============================================================================
// 特殊函数
// ============================================================================

/**
 * @brief 误差函数
 * @param x 输入值
 * @return erf(x)
 */
double errorFunction(double x);

/**
 * @brief 互补误差函数
 * @param x 输入值
 * @return erfc(x) = 1 - erf(x)
 */
double complementaryErrorFunction(double x);

/**
 * @brief 伽马函数
 * @param x 输入值
 * @return Γ(x)
 */
double gammaFunction(double x);

/**
 * @brief 贝塞尔函数 J₀
 * @param x 输入值
 * @return J₀(x)
 */
double besselJ0(double x);

/**
 * @brief 贝塞尔函数 J₁
 * @param x 输入值
 * @return J₁(x)
 */
double besselJ1(double x);

} // namespace MathUtils
} // namespace Basic
} // namespace SCDAT

#endif // SCDAT_UTILS_MATH_UTILS_H
