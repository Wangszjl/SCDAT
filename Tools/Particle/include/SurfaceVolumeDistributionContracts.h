#pragma once

#include <string>

namespace SCDAT
{
namespace Particle
{

struct VolumeCellDistributionStubContract
{
    std::string cell_id;
    std::string node_id;
    std::string boundary_group_id;
    double pseudo_volume_m3 = 0.0;
    double reference_potential_v = 0.0;
    double local_charge_density_c_per_m3 = 0.0;
    double center_x_m = 0.0;
    double center_y_m = 0.0;
    double center_z_m = 0.0;
    double characteristic_length_m = 0.0;
};

struct VolumeBoundaryFaceStubContract
{
    std::string face_id;
    std::string boundary_group_id;
    std::string node_id;
    int external_face_id = -1;
    double area_m2 = 0.0;
    double center_x_m = 0.0;
    double center_y_m = 0.0;
    double center_z_m = 0.0;
    double normal_x = 1.0;
    double normal_y = 0.0;
    double normal_z = 0.0;
};

struct SurfaceVolumeProjectionEntryContract
{
    std::string node_id;
    std::string cell_id;
    std::string boundary_group_id;
    double weight = 1.0;
    double volume_to_surface_weight = 1.0;
};

} // namespace Particle
} // namespace SCDAT
