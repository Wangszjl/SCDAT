/**
 * @file Vector3D_test.cpp
 * @brief Vector3D 单元测试
 * @details
 * 本文件验证向量代数与几何运算，包括：
 * 点积/叉积、长度、夹角、投影、安全归一化、绕轴旋转。
 */

#include "../include/Vector3D.h"
#include <cmath>
#include <gtest/gtest.h>

namespace SCDAT
{
namespace Geometry
{

TEST(Vector3DTest, DotCrossAndLength)
{
    // 标准正交基测试：x × y = z。
    Vector3D x = Vector3D::unitX();
    Vector3D y = Vector3D::unitY();

    Vector3D z = Vector3D(x.cross(y));

    EXPECT_DOUBLE_EQ(x.dot(y), 0.0);
    EXPECT_DOUBLE_EQ(z.x(), 0.0);
    EXPECT_DOUBLE_EQ(z.y(), 0.0);
    EXPECT_DOUBLE_EQ(z.z(), 1.0);
    EXPECT_DOUBLE_EQ(z.length(), 1.0);
}

TEST(Vector3DTest, AngleAndProjection)
{
    // 45° 夹角场景，便于检验反三角函数结果。
    Vector3D a(1.0, 1.0, 0.0);
    Vector3D b = Vector3D::unitX();

    const double angle = a.angleTo(b);
    EXPECT_NEAR(angle, std::acos(1.0 / std::sqrt(2.0)), 1e-12);

    Vector3D proj = a.projectionOnto(b);
    EXPECT_DOUBLE_EQ(proj.x(), 1.0);
    EXPECT_DOUBLE_EQ(proj.y(), 0.0);
    EXPECT_DOUBLE_EQ(proj.z(), 0.0);
}

TEST(Vector3DTest, SafeNormalizeZeroVector)
{
    // safeNormalized 对零向量应返回零向量而非抛异常。
    Vector3D zero = Vector3D::zero();
    Vector3D normalized = zero.safeNormalized();

    EXPECT_DOUBLE_EQ(normalized.x(), 0.0);
    EXPECT_DOUBLE_EQ(normalized.y(), 0.0);
    EXPECT_DOUBLE_EQ(normalized.z(), 0.0);
}

TEST(Vector3DTest, RotateAroundAxis)
{
    // 绕 Z 轴旋转 90°：X 轴应旋到 Y 轴。
    Vector3D v = Vector3D::unitX();
    Vector3D axis = Vector3D::unitZ();

    Vector3D rotated = v.rotateAroundAxis(axis, M_PI / 2.0);

    EXPECT_NEAR(rotated.x(), 0.0, 1e-12);
    EXPECT_NEAR(rotated.y(), 1.0, 1e-12);
    EXPECT_NEAR(rotated.z(), 0.0, 1e-12);
}

} // namespace Geometry
} // namespace SCDAT
