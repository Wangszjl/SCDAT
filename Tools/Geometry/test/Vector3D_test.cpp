#include "../include/Vector3D.h"
#include <cmath>
#include <gtest/gtest.h>

namespace SCDAT
{
namespace Geometry
{

TEST(Vector3DTest, DotCrossAndLength)
{
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
    Vector3D zero = Vector3D::zero();
    Vector3D normalized = zero.safeNormalized();

    EXPECT_DOUBLE_EQ(normalized.x(), 0.0);
    EXPECT_DOUBLE_EQ(normalized.y(), 0.0);
    EXPECT_DOUBLE_EQ(normalized.z(), 0.0);
}

TEST(Vector3DTest, RotateAroundAxis)
{
    Vector3D v = Vector3D::unitX();
    Vector3D axis = Vector3D::unitZ();

    Vector3D rotated = v.rotateAroundAxis(axis, M_PI / 2.0);

    EXPECT_NEAR(rotated.x(), 0.0, 1e-12);
    EXPECT_NEAR(rotated.y(), 1.0, 1e-12);
    EXPECT_NEAR(rotated.z(), 0.0, 1e-12);
}

} // namespace Geometry
} // namespace SCDAT
