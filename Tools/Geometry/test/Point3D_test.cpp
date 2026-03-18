#include "../include/Point3D.h"
#include <gtest/gtest.h>

namespace SCDAT
{
namespace Geometry
{

TEST(Point3DTest, ArithmeticOperators)
{
    Point3D p1(1.0, 2.0, 3.0);
    Point3D p2(4.0, 5.0, 6.0);

    Point3D sum = p1 + p2;
    Point3D diff = p2 - p1;

    EXPECT_DOUBLE_EQ(sum.x(), 5.0);
    EXPECT_DOUBLE_EQ(sum.y(), 7.0);
    EXPECT_DOUBLE_EQ(sum.z(), 9.0);

    EXPECT_DOUBLE_EQ(diff.x(), 3.0);
    EXPECT_DOUBLE_EQ(diff.y(), 3.0);
    EXPECT_DOUBLE_EQ(diff.z(), 3.0);
}

TEST(Point3DTest, DistanceAndMagnitude)
{
    Point3D origin(0.0, 0.0, 0.0);
    Point3D p(3.0, 4.0, 12.0);

    EXPECT_DOUBLE_EQ(origin.distanceTo(p), 13.0);
    EXPECT_DOUBLE_EQ(p.magnitude(), 13.0);
    EXPECT_DOUBLE_EQ(p.distanceSquaredTo(origin), 169.0);
}

TEST(Point3DTest, IndexOutOfRangeThrows)
{
    Point3D p(1.0, 2.0, 3.0);
    EXPECT_THROW((void)p[3], std::out_of_range);
    EXPECT_THROW((void)p[-1], std::out_of_range);
}

TEST(Point3DTest, NormalizeZeroThrows)
{
    Point3D p = Point3D::zero();
    EXPECT_THROW((void)p.normalized(), std::runtime_error);
}

} // namespace Geometry
} // namespace SCDAT
