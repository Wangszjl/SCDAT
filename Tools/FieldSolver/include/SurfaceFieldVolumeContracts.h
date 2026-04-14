#pragma once

#include <string>

namespace SCDAT
{
namespace FieldSolver
{

struct ExternalFieldSolveNodeResult
{
    std::string node_id;
    double reference_potential_v = 0.0;
    double normal_field_v_per_m = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double capacitance_scale = 1.0;
};

struct ExternalVolumeSolveCellResult
{
    std::string cell_id;
    std::string boundary_group_id;
    double potential_v = 0.0;
    double reference_potential_v = 0.0;
    double normal_field_v_per_m = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double capacitance_scale = 1.0;
    double coupling_gain = 0.0;
    double projection_weight_scale = 1.0;
    double sheath_length_scale = 1.0;
};

} // namespace FieldSolver
} // namespace SCDAT
