/**
 * @file Matrix_test.cpp
 * @brief Matrix/Matrix3x3 单元测试
 * @details
 * 覆盖内容：
 * 1) 基本构造与元素访问；
 * 2) 矩阵乘法数值正确性；
 * 3) 由 3x3 块拼装块对角矩阵；
 * 4) 块矩阵对 Vector3D 组的批量乘法。
 */

#include "../include/Matrix.h"
#include "../include/Matrix3x3.h"
#include <gtest/gtest.h>

namespace SCDAT
{
namespace Geometry
{

TEST(MatrixTest, ConstructAndAccess)
{
    // 验证行列规模、下标访问和写回逻辑。
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
    // 使用手算可验证的样例检查矩阵乘法结果。
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
    // 检查块对角构造是否保持块内值并正确填充非对角零元。
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
    // 检查 blockDiagonalFrom3x3 与 multiplyBlocks3 的组合一致性。
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
