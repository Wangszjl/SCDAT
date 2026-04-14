#pragma once

#include <string>

namespace SCDAT
{
namespace Collision
{

struct PicMccControlContract
{
    bool enable_mcc = true;
    std::string collision_cross_section_set_id = "surface_pic_v1";
};

struct PicMccTelemetryContract
{
    bool mcc_enabled = false;
    int total_collisions = 0;
};

} // namespace Collision
} // namespace SCDAT
