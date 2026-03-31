#ifndef SCDAT_FIELD_BOUNDARY_TYPE_H
#define SCDAT_FIELD_BOUNDARY_TYPE_H

namespace SCDAT
{

/**
 * @brief 电场求解边界条件类型
 */
enum class BoundaryType
{
    STANDARD = 0,   ///< 标准边界条件 (BC=0)
    MODIFIED_1 = 1, ///< 修改的边界条件1 (BC=1)
    RADIAL = 2,     ///< 径向边界条件 (BC=2)
    PERIODIC = 3,   ///< 周期性边界条件 (BC=3)
    MODIFIED_4 = 4, ///< 修改的边界条件4 (BC=4)
    DIRICHLET = 5,  ///< 狄利克雷边界条件
    NEUMANN = 6,    ///< 诺伊曼边界条件
    MIXED = 7       ///< 混合边界条件
};

/**
 * @brief 粒子推进模块边界类型
 */
enum class ParticleBoundaryType
{
    REFLECTING = 0,          ///< 反射边界
    ABSORBING = 1,           ///< 吸收边界
    PERIODIC = 2,            ///< 周期边界
    OPEN = 3,                ///< 开放边界
    STANDARD = 4,            ///< 标准边界条件
    MODIFIED = 5,            ///< 修改的边界条件
    RADIAL = 6,              ///< 径向边界条件
    CYLINDRICAL_PERIODIC = 7 ///< 柱坐标周期边界
};

namespace Mesh
{

/**
 * @brief 网格模块边界类型
 */
enum class BoundaryType
{
    INTERIOR = 0,
    DIRICHLET = 1,
    NEUMANN = 2,
    MIXED = 3,
    PERIODIC = 4
};

} // namespace Mesh

namespace Particle
{

/**
 * @brief 高级边界条件类型枚举
 */
enum class AdvancedBoundaryType
{
    PERIODIC_3D = 0,           ///< 3D周期边界条件
    REFLECTIVE_CURVED = 1,     ///< 曲面反射边界条件
    ABSORBING_STATISTICAL = 2, ///< 统计吸收边界条件
    EMISSION_THERMAL = 3,      ///< 热发射边界条件
    EMISSION_FIELD = 4,        ///< 场发射边界条件
    EMISSION_SECONDARY = 5     ///< 二次发射边界条件
};

/**
 * @brief 边界几何类型
 */
enum class BoundaryGeometry
{
    PLANAR = 0,      ///< 平面边界
    CYLINDRICAL = 1, ///< 圆柱面边界
    SPHERICAL = 2,   ///< 球面边界
    ARBITRARY = 3    ///< 任意曲面边界
};

} // namespace Particle

namespace FieldSolver
{

/**
 * @brief 解器边界条件类型
 */
enum class BoundaryConditionType
{
    DIRICHLET = 1, ///< 固定电位
    NEUMANN = 2,   ///< 固定电场
    ROBIN = 3,     ///< Robin 边界条件
    PERIODIC = 4,  ///< 周期边界
    MIXED = 5,     ///< 混合边界
    FLOATING = 6   ///< 浮动电位
};

} // namespace FieldSolver

} // namespace SCDAT




#endif // SCDAT_FIELD_BOUNDARY_TYPE_H
