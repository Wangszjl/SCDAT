#include "ArcPICEmissionModel.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

ArcPICEmissionModel::ArcPICEmissionModel()
{
    rebuildStrategy();
}

void ArcPICEmissionModel::setStrategyType(ArcEmissionStrategyType type)
{
    if (strategy_type_ == type)
    {
        return;
    }

    strategy_type_ = type;
    rebuildStrategy();
}

void ArcPICEmissionModel::configureStrategy(const ArcEmissionStrategyParameters& parameters)
{
    strategy_parameters_ = parameters;
    rebuildStrategy();
}

std::string ArcPICEmissionModel::getStrategyName() const
{
    if (!strategy_)
    {
        return "uninitialized";
    }

    return strategy_->name();
}

void ArcPICEmissionModel::rebuildStrategy()
{
    strategy_ = createArcEmissionStrategy(strategy_type_, strategy_parameters_);
}

double ArcPICEmissionModel::computeEmissionCurrentDensity(double electric_field_v_per_m,
                                                          double cathode_temperature_k) const
{
    if (!strategy_)
    {
        return 0.0;
    }

    return strategy_->computeCurrentDensity(electric_field_v_per_m, cathode_temperature_k);
}

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
