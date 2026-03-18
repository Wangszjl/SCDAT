/**
 * @file Point3D.cpp
 * @brief 三维点几何运算实现
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026年3月18日 7:46:08
 *
 * 实现三维点的基本几何运算、距离计算、坐标变换和空间关系判断功能。
 * 支持点与点、点与线、点与面的几何关系计算。
 *
 * @note 坐标系遵循右手定则，距离计算使用欧几里得度量
 */

#include "../include/Point3D.h"
#include <functional>
#include <iomanip> // 格式化输出需要
#include <sstream>

// ============================================================================
// 静态常量定义
// ============================================================================

/**
 * @brief 数值计算精度常量
 */
const double SCDAT::Geometry::Point3D::EPSILON = 1e-12;

// ============================================================================
// 工具函数实现
// ============================================================================

/**
 * @brief 计算点的哈希值，用于哈希表存储
 *
 * 使用坐标分量的哈希值组合生成唯一哈希值。
 *
 * @return 哈希值
 *
 * @note 相同坐标的点具有相同的哈希值
 */
std::size_t SCDAT::Geometry::Point3D::hash() const noexcept
{
    // 使用更好的哈希组合算法，减少冲突
    std::size_t h1 = std::hash<double>{}(x_);
    std::size_t h2 = std::hash<double>{}(y_);
    std::size_t h3 = std::hash<double>{}(z_);

    // 使用boost::hash_combine的算法
    std::size_t seed = h1;
    seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

/**
 * @brief 将点转换为字符串表示
 *
 * @return 格式化的字符串
 *
 * @note 输出格式：Point3D(x, y, z)
 */
std::string SCDAT::Geometry::Point3D::toString() const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6); // 控制精度
    oss << "Point3D(" << x_ << ", " << y_ << ", " << z_ << ")";
    return oss.str();
}
