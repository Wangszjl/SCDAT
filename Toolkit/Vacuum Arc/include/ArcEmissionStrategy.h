#pragma once

#include <memory>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

enum class ArcEmissionStrategyType
{
    LinearThreshold,
    FowlerNordheimLike,
};

struct ArcEmissionStrategyParameters
{
  // Shared modulation parameters.
  double field_enhancement_factor = 1.0;
  double regional_gain = 1.0;
  double thermal_activation_exponent = 1.0;

    // Linear-threshold model parameters.
    double linear_field_threshold_v_per_m = 1.5e7;
    double linear_field_gain = 1.0e-6;
    double linear_temperature_threshold_k = 500.0;
    double linear_temperature_gain = 1.0e-2;

    // Fowler-Nordheim-like model parameters.
    double fn_prefactor = 1.0e-18;
    double fn_barrier_v_per_m = 4.0e9;
    double fn_temperature_reference_k = 500.0;
    double fn_temperature_gain = 5.0e-3;
};

class ArcEmissionStrategy
{
  public:
    virtual ~ArcEmissionStrategy() = default;

    virtual double computeCurrentDensity(double electric_field_v_per_m,
                                         double cathode_temperature_k) const = 0;
    virtual std::string name() const = 0;
};

class LinearThresholdEmissionStrategy : public ArcEmissionStrategy
{
  public:
    explicit LinearThresholdEmissionStrategy(const ArcEmissionStrategyParameters& parameters)
        : parameters_(parameters)
    {
    }

    double computeCurrentDensity(double electric_field_v_per_m,
                                 double cathode_temperature_k) const override;

    std::string name() const override
    {
        return "linear_threshold";
    }

  private:
    ArcEmissionStrategyParameters parameters_{};
};

class FowlerNordheimLikeEmissionStrategy : public ArcEmissionStrategy
{
  public:
    explicit FowlerNordheimLikeEmissionStrategy(const ArcEmissionStrategyParameters& parameters)
        : parameters_(parameters)
    {
    }

    double computeCurrentDensity(double electric_field_v_per_m,
                                 double cathode_temperature_k) const override;

    std::string name() const override
    {
        return "fowler_nordheim_like";
    }

  private:
    ArcEmissionStrategyParameters parameters_{};
};

std::unique_ptr<ArcEmissionStrategy>
createArcEmissionStrategy(ArcEmissionStrategyType type,
                          const ArcEmissionStrategyParameters& parameters);

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
