#include "DensePlasmaBoundaryLayer.h"

#include <algorithm>
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
    const double ion_mass_amu = std::clamp(parameters.ion_mass_amu, 1.0, 240.0);
    const double ion_mass = ion_mass_amu * kAtomicMassUnit;
    const double electron_temperature_j =
        std::max(1.0e-3, parameters.electron_temperature_ev) * kElementaryCharge;
    const double debye_length_m = std::max(assessment.debye_length_m, 1.0e-9);
    const double collisionality = std::max(0.0, assessment.collisionality);

    const double dense_base = assessment.is_dense ? 2.8 : 4.6;
    const double collisionality_boost = 1.0 + std::min(4.0, 0.75 * collisionality);
    const double sheath_scale = dense_base * collisionality_boost;

    const double presheath_base_scale = 3.0 + std::min(5.0, 1.4 * collisionality);
    const double presheath_scale = std::max(presheath_base_scale, 5.0);

    const double sheath_drop_fraction = std::clamp(
        0.68 + (assessment.is_dense ? 0.06 : 0.0) + 0.03 * std::min(2.0, collisionality),
        0.60,
        0.88);

    state.bohm_velocity_m_per_s = std::sqrt(electron_temperature_j / ion_mass);
    state.ion_mach_at_sheath_edge = 1.0 + 0.08 * std::min(3.0, collisionality);
    state.sheath_thickness_m = std::max(debye_length_m * sheath_scale, 2.0e-7);
    state.presheath_thickness_m =
        std::max(debye_length_m * presheath_scale, 2.5 * state.sheath_thickness_m);
    state.sheath_potential_drop_v = surface_potential_v * sheath_drop_fraction;
    state.presheath_potential_drop_v = surface_potential_v - state.sheath_potential_drop_v;
    state.debye_to_sheath_ratio = debye_length_m / std::max(state.sheath_thickness_m, 1.0e-9);
    state.wall_field_v_per_m =
        state.sheath_potential_drop_v / std::max(state.sheath_thickness_m, 1.0e-9);
    return state;
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
