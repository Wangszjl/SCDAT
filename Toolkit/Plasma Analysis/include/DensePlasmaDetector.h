#pragma once

#include "FluidAlgorithmConfig.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

class DensePlasmaDetector
{
  public:
    explicit DensePlasmaDetector(double dense_threshold_m3 = 1.0e18)
        : dense_threshold_m3_(dense_threshold_m3)
    {
    }

    DensePlasmaAssessment assess(const PlasmaParameters& parameters) const;

  private:
    double dense_threshold_m3_ = 1.0e18;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
