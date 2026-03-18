#include "../include/Matrix.h"
#include "../include/Matrix3x3.h"
#include <gtest/gtest.h>

namespace SCDAT
{
namespace Geometry
{

TEST(MatrixTest, ConstructAndAccess)
{
    Matrix m(2, 3, 0.0);
    m(0, 1) = 2.5;
    m(1, 2) = -4.0;

    EXPECT_DOUBLE_EQ(m(0, 1), 2.5);
    EXPECT_DOUBLE_EQ(m(1, 2), -4.0);
    EXPECT_EQ(m.rows(), 2u);
    EXPECT_EQ(m.cols(), 3u);
}

TEST(MatrixTest, MatrixMultiplication)
{
    Matrix a(2, 3, std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    Matrix b(3, 2, std::vector<double>{7.0, 8.0, 9.0, 10.0, 11.0, 12.0});

    Matrix c = a * b;

    EXPECT_EQ(c.rows(), 2u);
    EXPECT_EQ(c.cols(), 2u);
    EXPECT_DOUBLE_EQ(c(0, 0), 58.0);
    EXPECT_DOUBLE_EQ(c(0, 1), 64.0);
    EXPECT_DOUBLE_EQ(c(1, 0), 139.0);
    EXPECT_DOUBLE_EQ(c(1, 1), 154.0);
}

TEST(MatrixTest, BuildFromMatrix3x3Blocks)
{
    Matrix3x3 b0 = Matrix3x3::identity();
    Matrix3x3 b1 = Matrix3x3::diagonal(2.0, 3.0, 4.0);

    Matrix blockDiag = Matrix::blockDiagonalFrom3x3({b0, b1});

    EXPECT_EQ(blockDiag.rows(), 6u);
    EXPECT_EQ(blockDiag.cols(), 6u);

    EXPECT_DOUBLE_EQ(blockDiag(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(blockDiag(2, 2), 1.0);
    EXPECT_DOUBLE_EQ(blockDiag(3, 3), 2.0);
    EXPECT_DOUBLE_EQ(blockDiag(4, 4), 3.0);
    EXPECT_DOUBLE_EQ(blockDiag(5, 5), 4.0);

    EXPECT_DOUBLE_EQ(blockDiag(0, 3), 0.0);
    EXPECT_DOUBLE_EQ(blockDiag(4, 1), 0.0);
}

TEST(MatrixTest, MultiplyBlockVectors)
{
    Matrix3x3 b0 = Matrix3x3::identity();
    Matrix3x3 b1 = Matrix3x3::diagonal(2.0, 3.0, 4.0);

    Matrix blockDiag = Matrix::blockDiagonalFrom3x3({b0, b1});

    std::vector<Vector3D> input = {Vector3D(1.0, 2.0, 3.0), Vector3D(1.0, 1.0, 1.0)};
    std::vector<Vector3D> output = blockDiag.multiplyBlocks3(input);

    ASSERT_EQ(output.size(), 2u);

    EXPECT_DOUBLE_EQ(output[0].x(), 1.0);
    EXPECT_DOUBLE_EQ(output[0].y(), 2.0);
    EXPECT_DOUBLE_EQ(output[0].z(), 3.0);

    EXPECT_DOUBLE_EQ(output[1].x(), 2.0);
    EXPECT_DOUBLE_EQ(output[1].y(), 3.0);
    EXPECT_DOUBLE_EQ(output[1].z(), 4.0);
}

} // namespace Geometry
} // namespace SCDAT
