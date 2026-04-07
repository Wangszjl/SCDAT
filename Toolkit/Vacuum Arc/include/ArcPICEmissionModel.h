#pragma once

#include "ArcEmissionStrategy.h"

#include <memory>
#include <string>

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class ArcPICEmissionModel
{
  public:
    ArcPICEmissionModel();

    void setStrategyType(ArcEmissionStrategyType type);
    ArcEmissionStrategyType getStrategyType() const
    {
        return strategy_type_;
    }

    void configureStrategy(const ArcEmissionStrategyParameters& parameters);
    const ArcEmissionStrategyParameters& getStrategyParameters() const
    {
        return strategy_parameters_;
    }

    std::string getStrategyName() const;

    double computeEmissionCurrentDensity(double electric_field_v_per_m,
                                         double cathode_temperature_k) const;

  private:
    void rebuildStrategy();

    ArcEmissionStrategyType strategy_type_ = ArcEmissionStrategyType::LinearThreshold;
    ArcEmissionStrategyParameters strategy_parameters_{};
    std::unique_ptr<ArcEmissionStrategy> strategy_;
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
