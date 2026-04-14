#include "../include/SurfaceFieldVolumeBridge.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

TEST(SurfaceFieldVolumeBridgeTest, BlendFieldVolumeScalarClampsBlendFactor)
{
    EXPECT_NEAR(SCDAT::FieldSolver::blendFieldVolumeScalar(4.0, 10.0, -1.0), 4.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::FieldSolver::blendFieldVolumeScalar(4.0, 10.0, 0.5), 7.0, 1.0e-12);
    EXPECT_NEAR(SCDAT::FieldSolver::blendFieldVolumeScalar(4.0, 10.0, 2.0), 10.0, 1.0e-12);
}

TEST(SurfaceFieldVolumeBridgeTest, BlendFieldVolumeScalarSanitizesNonFiniteInputs)
{
    const double value = SCDAT::FieldSolver::blendFieldVolumeScalar(
        std::numeric_limits<double>::quiet_NaN(), 8.0, 0.5);
    EXPECT_TRUE(std::isfinite(value));
    EXPECT_NEAR(value, 4.0, 1.0e-12);
}

TEST(SurfaceFieldVolumeBridgeTest, UpdateMismatchMetricUsesMaxNormalizedDeviation)
{
    double metric = 0.0;
    metric = SCDAT::FieldSolver::updateFieldVolumeMismatchMetric(metric, 5.0, 4.0, 2.0);
    EXPECT_NEAR(metric, 0.5, 1.0e-12);

    metric = SCDAT::FieldSolver::updateFieldVolumeMismatchMetric(metric, 9.0, 1.0, 2.0);
    EXPECT_NEAR(metric, 4.0, 1.0e-12);

    metric = SCDAT::FieldSolver::updateFieldVolumeMismatchMetric(
        metric, std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0);
    EXPECT_NEAR(metric, 4.0, 1.0e-12);
}

TEST(SurfaceFieldVolumeBridgeTest,
     ComputeExternalBlendFactorDecreasesWithMismatchAndNonConvergence)
{
    const double converged_small_mismatch =
        SCDAT::FieldSolver::computeFieldVolumeExternalBlendFactor(0.1, 1.0e-8, 1.0e-8, true,
                                                                   1.0e-6, 0.8);
    const double unconverged_large_mismatch =
        SCDAT::FieldSolver::computeFieldVolumeExternalBlendFactor(10.0, 1.0e-2, 1.0e-2, false,
                                                                   1.0e-6, 0.8);

    EXPECT_GE(converged_small_mismatch, 0.05);
    EXPECT_LE(converged_small_mismatch, 1.0);
    EXPECT_GE(unconverged_large_mismatch, 0.05);
    EXPECT_LE(unconverged_large_mismatch, 1.0);
    EXPECT_GT(converged_small_mismatch, unconverged_large_mismatch);
}
