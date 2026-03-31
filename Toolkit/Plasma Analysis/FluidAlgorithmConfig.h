#pragma once

#include "../../Tools/Geometry/include/Vector3D.h"

#include <cstddef>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct PlasmaParameters
{
    double electron_density_m3 = 1.0e17;
    double ion_density_m3 = 1.0e17;
    double electron_temperature_ev = 5.0;
    double ion_temperature_ev = 0.2;
    double ion_mass_amu = 40.0;
    double neutral_density_m3 = 2.5e20;
    double electric_field_v_per_m = 0.0;
    double pressure_pa = 1.0;
};

struct FluidAlgorithmConfig
{
    Geometry::Vector3D domain_size{1.0e-3, 1.0e-3, 5.0e-3};
    Geometry::Vector3D resolution{6.0, 6.0, 24.0};
    double time_step_s = 1.0e-9;
    double dense_plasma_threshold_m3 = 1.0e18;
    double initial_potential_v = 25.0;
    PlasmaParameters initial_plasma;
};

struct DensePlasmaAssessment
{
    bool is_dense = false;
    double debye_length_m = 0.0;
    double plasma_frequency_hz = 0.0;
    double collisionality = 0.0;
};

struct FluidAlgorithmStatus
{
    double time_s = 0.0;
    bool dense_plasma_detected = false;
    double average_potential_v = 0.0;
    double average_density_m3 = 0.0;
    double sheath_thickness_m = 0.0;
    double debye_length_m = 0.0;
    std::size_t steps_completed = 0;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
