#include "ArcEmissionStrategy.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

double LinearThresholdEmissionStrategy::computeCurrentDensity(double electric_field_v_per_m,
                                                              double cathode_temperature_k) const
{
    const double enhanced_field =
        std::max(0.0, electric_field_v_per_m) *
        std::max(0.0, parameters_.field_enhancement_factor);
    const double field_term =
        std::max(0.0, enhanced_field - parameters_.linear_field_threshold_v_per_m) *
        parameters_.linear_field_gain;
    const double thermal_excess =
        std::max(0.0, cathode_temperature_k - parameters_.linear_temperature_threshold_k);
    const double temperature_term =
        std::pow(thermal_excess, std::max(0.0, parameters_.thermal_activation_exponent)) *
        parameters_.linear_temperature_gain;
    return (field_term + temperature_term) * std::max(0.0, parameters_.regional_gain);
}

double FowlerNordheimLikeEmissionStrategy::computeCurrentDensity(double electric_field_v_per_m,
                                                                  double cathode_temperature_k) const
{
    const double enhanced_field =
        std::max(0.0, electric_field_v_per_m) *
        std::max(0.0, parameters_.field_enhancement_factor);
    const double bounded_field = std::max(1.0, enhanced_field);
    const double field_component =
        parameters_.fn_prefactor * bounded_field * bounded_field *
        std::exp(-parameters_.fn_barrier_v_per_m / bounded_field);
    const double thermal_excess =
        std::max(0.0, cathode_temperature_k - parameters_.fn_temperature_reference_k);
    const double thermal_component =
        std::pow(thermal_excess, std::max(0.0, parameters_.thermal_activation_exponent)) *
        parameters_.fn_temperature_gain;
    return (field_component + thermal_component) * std::max(0.0, parameters_.regional_gain);
}

std::unique_ptr<ArcEmissionStrategy>
createArcEmissionStrategy(ArcEmissionStrategyType type,
                          const ArcEmissionStrategyParameters& parameters)
{
    switch (type)
    {
    case ArcEmissionStrategyType::LinearThreshold:
        return std::make_unique<LinearThresholdEmissionStrategy>(parameters);
    case ArcEmissionStrategyType::FowlerNordheimLike:
        return std::make_unique<FowlerNordheimLikeEmissionStrategy>(parameters);
    default:
        throw std::invalid_argument("Unknown ArcEmissionStrategyType");
    }
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
