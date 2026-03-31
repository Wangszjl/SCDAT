/**
 * @file MeshPartitioning.h
 * @brief 网格空间划分与查询结构体定义，包括包围盒、八叉树、空间索引。
 * @author Wang Sizhan
 * @version V0.0.1
 * @date 2026.3.24 13:30:42
 */
#ifndef SCDAT_MESH_PARTITIONING_H
#define SCDAT_MESH_PARTITIONING_H

#include "Geometry/include/Point3D.h"
#include "Geometry/include/Vector3D.h"
#include "MeshAlgorithms.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Mesh
{

/**
 * @class BoundingBox
 * @brief 三维空间包围盒，用于空间划分和碰撞检测。
 */
class BoundingBox
{
  public:
    /// @brief 默认构造函数
    BoundingBox() = default;
    /// @brief 构造函数
    BoundingBox(const Utils::Point3D& min_pt, const Utils::Point3D& max_pt)
        : min_point_(min_pt), max_point_(max_pt)
    {
    }

    /// @brief 获取最小点
    const Utils::Point3D& getMinPoint() const
    {
        return min_point_;
    }
    /// @brief 获取最大点
    const Utils::Point3D& getMaxPoint() const
    {
        return max_point_;
    }

    /// @brief 获取中心点
    Utils::Point3D getCenter() const;
    /// @brief 获取尺寸向量
    Utils::Vector3D getSize() const;

    /// @brief 判断点是否在包围盒内
    bool contains(const Utils::Point3D& point) const;
    /// @brief 判断与另一个包围盒是否相交
    bool intersects(const BoundingBox& other) const;

    /// @brief 扩展包围盒以包含指定点
    void expand(const Utils::Point3D& point);
    /// @brief 八分割子盒
    std::array<BoundingBox, 8> subdivide() const;

    /// @brief 转字符串
    std::string toString() const;

  private:
    Utils::Point3D min_point_{}; ///< 最小点
    Utils::Point3D max_point_{}; ///< 最大点
};

class OctreeNode;
using OctreeNodePtr = std::shared_ptr<OctreeNode>;

/**
 * @class OctreeNode
 * @brief 八叉树节点，用于空间层次划分和快速查询。
 */
class OctreeNode
{
  public:
    /**
     * @brief 构造函数
     * @param bounds 空间包围盒
     * @param depth 当前深度
     * @param max_depth 最大深度
     * @param max_elements 每节点最大元素数
     */
    OctreeNode(const BoundingBox& bounds, int depth, int max_depth, int max_elements);

    /// @brief 添加单元到节点
    void addElement(ElementPtr element);
    /// @brief 清空节点
    void clear();

    /// @brief 查询包含指定点的所有单元
    std::vector<ElementPtr> query(const Utils::Point3D& point) const;
    /// @brief 查询与区域相交的所有单元
    std::vector<ElementPtr> query(const BoundingBox& region) const;
    /// @brief 查询半径范围内的所有单元
    std::vector<ElementPtr> queryRadius(const Utils::Point3D& center, double radius) const;

    /// @brief 获取当前节点元素数
    int getElementCount() const;
    /// @brief 获取节点总数
    int getNodeCount() const;
    /// @brief 获取最大深度
    int getMaxDepth() const;

    /// @brief 获取节点包围盒
    const BoundingBox& getBounds() const
    {
        return bounds_;
    }

  private:
    void subdivide();
    int getChildIndex(const Utils::Point3D& point) const;

    BoundingBox bounds_;    ///< 空间包围盒
    int depth_ = 0;         ///< 当前深度
    int max_depth_ = 8;     ///< 最大深度
    int max_elements_ = 16; ///< 每节点最大元素数
    bool is_leaf_ = true;   ///< 是否叶子节点

    std::vector<ElementPtr> elements_;        ///< 元素集合
    std::array<OctreeNodePtr, 8> children_{}; ///< 子节点
};

/**
 * @class SpatialIndex
 * @brief 空间索引结构，基于八叉树实现高效空间查询。
 */
class SpatialIndex
{
  public:
    /**
     * @brief 构造函数
     * @param max_depth 最大深度
     * @param max_elements_per_node 每节点最大元素数
     */
    SpatialIndex(int max_depth = 8, int max_elements_per_node = 16);

    /// @brief 用元素集合初始化索引
    void initialize(const std::vector<ElementPtr>& elements);
    /// @brief 用指定包围盒初始化索引
    void initialize(const BoundingBox& bounds);

    /// @brief 添加单元到索引
    void addElement(ElementPtr element);
    /// @brief 清空索引
    void clear();

    /// @brief 查询包含点的所有单元
    std::vector<ElementPtr> findElementsContaining(const Utils::Point3D& point) const;
    /// @brief 查询区域内所有单元
    std::vector<ElementPtr> findElementsInRegion(const BoundingBox& region) const;
    /// @brief 查询半径范围内所有单元
    std::vector<ElementPtr> findElementsInRadius(const Utils::Point3D& center, double radius) const;

    /// @brief 获取总元素数
    int getElementCount() const
    {
        return total_elements_;
    }
    /// @brief 获取节点总数
    int getNodeCount() const;
    /// @brief 获取最大深度
    int getMaxDepth() const;

  private:
    BoundingBox calculateBounds(const std::vector<ElementPtr>& elements) const;

    int max_depth_ = 8;              ///< 最大深度
    int max_elements_per_node_ = 16; ///< 每节点最大元素数
    OctreeNodePtr root_;             ///< 根节点
    int total_elements_ = 0;         ///< 总元素数
};

using SpatialIndexPtr = std::shared_ptr<SpatialIndex>;

} // namespace Mesh
} // namespace SCDAT

#endif // SCDAT_MESH_PARTITIONING_H
