#pragma once

#include "FluidAlgorithmConfig.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct BoundaryLayerState
{
    double sheath_thickness_m = 0.0;
  double presheath_thickness_m = 0.0;
    double bohm_velocity_m_per_s = 0.0;
  double ion_mach_at_sheath_edge = 1.0;
  double sheath_potential_drop_v = 0.0;
  double presheath_potential_drop_v = 0.0;
  double debye_to_sheath_ratio = 0.0;
    double wall_field_v_per_m = 0.0;
};

class DensePlasmaBoundaryLayer
{
  public:
    BoundaryLayerState analyze(const PlasmaParameters& parameters,
                               const DensePlasmaAssessment& assessment,
                               double surface_potential_v) const;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
