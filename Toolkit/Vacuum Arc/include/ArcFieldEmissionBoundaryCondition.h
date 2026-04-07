#pragma once

#include "../../Tools/Boundary/include/FieldEmissionBoundaryCondition.h"

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

class ArcFieldEmissionBoundaryCondition : public Particle::FieldEmissionBoundaryCondition
{
  public:
    using Parameters = Particle::FieldEmissionBoundaryCondition::Parameters;

    ArcFieldEmissionBoundaryCondition() = default;
    explicit ArcFieldEmissionBoundaryCondition(const Parameters& parameters)
        : Particle::FieldEmissionBoundaryCondition(parameters)
    {
    }
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
