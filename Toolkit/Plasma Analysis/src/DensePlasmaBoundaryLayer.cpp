#include "DensePlasmaBoundaryLayer.h"

#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kAtomicMassUnit = 1.66053906660e-27;
}

BoundaryLayerState DensePlasmaBoundaryLayer::analyze(const PlasmaParameters& parameters,
                                                     const DensePlasmaAssessment& assessment,
                                                     double surface_potential_v) const
{
    BoundaryLayerState state;
    const double ion_mass = 40.0 * kAtomicMassUnit;
    const double electron_temperature_j =
        std::max(1.0e-3, parameters.electron_temperature_ev) * kElementaryCharge;

    state.bohm_velocity_m_per_s = std::sqrt(electron_temperature_j / ion_mass);
    state.sheath_thickness_m =
        std::max(assessment.debye_length_m * (assessment.is_dense ? 3.0 : 8.0), 1.0e-6);
    state.wall_field_v_per_m = surface_potential_v / std::max(state.sheath_thickness_m, 1.0e-9);
    return state;
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
