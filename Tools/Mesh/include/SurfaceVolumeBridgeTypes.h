#pragma once

#include <string>

namespace SCDAT
{
namespace Mesh
{

struct ExternalVolumeMeshCellInput
{
    std::string cell_id;
    std::string node_id;
    std::string boundary_group_id;
    double node_area_m2 = 0.0;
    double characteristic_length_m = 0.0;
    double center_x_m = 0.0;
    double center_y_m = 0.0;
    double center_z_m = 0.0;
    double cell_volume_m3 = 0.0;
    double initial_potential_v = 0.0;
    double initial_charge_density_c_per_m3 = 0.0;
};

struct ExternalVolumeMeshNeighborInput
{
    std::string source_cell_id;
    std::string target_cell_id;
    double conductance_s = 0.0;
    double resistance_ohm = 0.0;
    double shared_face_area_m2 = 0.0;
    double face_distance_m = 0.0;
    double permittivity_scale = 1.0;
    double face_center_x_m = 0.0;
    double face_center_y_m = 0.0;
    double face_center_z_m = 0.0;
};

struct ExternalVolumeMeshFaceInput
{
    std::string face_id;
    std::string boundary_group_id;
    std::string node_id;
    double area_m2 = 0.0;
    double center_x_m = 0.0;
    double center_y_m = 0.0;
    double center_z_m = 0.0;
    double normal_x = 1.0;
    double normal_y = 0.0;
    double normal_z = 0.0;
};

struct ExternalVolumeMeshCellFaceInput
{
    std::string cell_id;
    std::string face_id;
    std::string role;
    double projection_weight = 1.0;
};

struct ExternalSurfaceProjectionInput
{
    std::string node_id;
    std::string cell_id;
    double surface_to_volume_weight = 1.0;
    double volume_to_surface_weight = 1.0;
};

} // namespace Mesh
} // namespace SCDAT
