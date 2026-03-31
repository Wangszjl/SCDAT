
/**
 * @file MeshAlgorithms.h
 * @brief 网格基础数据结构与元素算法定义，包含节点、单元、基本几何操作。
 * @author 自动生成
 */
#ifndef SCDAT_MESH_ALGORITHMS_H
#define SCDAT_MESH_ALGORITHMS_H

#include "../../Boundary/include/BoundaryType.h"
#include "Geometry/include/Matrix3x3.h"
#include "Geometry/include/Point3D.h"
#include "Geometry/include/Vector3D.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
/// @brief Geometry 工具别名
namespace Utils = Geometry;

/// @brief 网格相关命名空间
namespace Mesh
{

/// @brief 节点编号类型
using NodeId = std::size_t;
/// @brief 单元编号类型
using ElementId = std::size_t;
/// @brief 材料编号类型
using MaterialId = std::size_t;

enum class MaterialType
{
    VACUUM = 0,
    DIELECTRIC = 1,
    CONDUCTOR = 2
};

enum class ElementType
{
    TETRAHEDRON = 4,
    TRIANGLE = 3,
    LINE = 2,
    POINT = 1
};

class Node;
class Element;
class TetrahedronElement;
class TriangleElement;
class Mesh;

using NodePtr = std::shared_ptr<Node>;
using ElementPtr = std::shared_ptr<Element>;
using TetrahedronPtr = std::shared_ptr<TetrahedronElement>;
using TrianglePtr = std::shared_ptr<TriangleElement>;
using MeshPtr = std::shared_ptr<Mesh>;

/**
 * @class Node
 * @brief 网格节点，包含空间位置、电势、电场等属性。
 */
class Node
{
  public:
    /**
     * @brief 默认构造函数
     */
    Node() = default;
    /**
     * @brief 构造函数
     * @param id 节点编号
     * @param position 空间坐标
     */
    Node(NodeId id, const Utils::Point3D& position) : id_(id), position_(position) {}

    /// @brief 获取节点编号
    NodeId getId() const
    {
        return id_;
    }
    /// @brief 获取节点位置（常量）
    const Utils::Point3D& getPosition() const
    {
        return position_;
    }
    /// @brief 获取节点位置（可修改）
    Utils::Point3D& getPosition()
    {
        return position_;
    }

    /// @brief 获取边界类型
    BoundaryType getBoundaryType() const
    {
        return boundary_type_;
    }
    /// @brief 设置边界类型
    void setBoundaryType(BoundaryType type)
    {
        boundary_type_ = type;
    }

    /// @brief 获取节点电势
    double getPotential() const
    {
        return potential_;
    }
    /// @brief 设置节点电势
    void setPotential(double potential)
    {
        potential_ = potential;
    }

    /// @brief 获取节点电场
    const Utils::Vector3D& getElectricField() const
    {
        return electric_field_;
    }
    /// @brief 设置节点电场
    void setElectricField(const Utils::Vector3D& field)
    {
        electric_field_ = field;
    }

    /// @brief 获取节点电荷密度
    double getChargeDensity() const
    {
        return charge_density_;
    }
    /// @brief 设置节点电荷密度
    void setChargeDensity(double density)
    {
        charge_density_ = density;
    }

    /// @brief 节点信息转字符串
    std::string toString() const;

  private:
    NodeId id_ = 0;                                       ///< 节点编号
    Utils::Point3D position_{};                           ///< 空间坐标
    BoundaryType boundary_type_ = BoundaryType::INTERIOR; ///< 边界类型
    double potential_ = 0.0;                              ///< 电势
    Utils::Vector3D electric_field_{};                    ///< 电场
    double charge_density_ = 0.0;                         ///< 电荷密度
};

/**
 * @class Element
 * @brief 网格单元基类，支持多种类型元素。
 */
class Element
{
  public:
    /**
     * @brief 构造函数
     * @param id 单元编号
     * @param type 单元类型
     * @param nodes 节点指针数组
     */
    Element(ElementId id, ElementType type, std::vector<NodePtr> nodes,
            MaterialId material_id = 0)
        : id_(id), type_(type), nodes_(std::move(nodes)), material_id_(material_id)
    {
    }
    virtual ~Element() = default;

    /// @brief 获取单元编号
    ElementId getId() const
    {
        return id_;
    }
    /// @brief 获取单元类型
    ElementType getType() const
    {
        return type_;
    }
    /// @brief 获取单元节点
    const std::vector<NodePtr>& getNodes() const
    {
        return nodes_;
    }
    const std::vector<ElementPtr>& getNeighbors() const
    {
        return neighbors_;
    }
    void addNeighbor(const ElementPtr& neighbor)
    {
        if (neighbor)
        {
            neighbors_.push_back(neighbor);
        }
    }
    void clearNeighbors()
    {
        neighbors_.clear();
    }
    MaterialId getMaterialId() const
    {
        return material_id_;
    }
    void setMaterialId(MaterialId material_id)
    {
        material_id_ = material_id;
    }

    /// @brief 获取单元体积（或面积）
    virtual double getVolume() const = 0;
    /// @brief 获取单元重心
    virtual Utils::Point3D getCentroid() const = 0;
    /// @brief 获取特征长度
    virtual double getCharacteristicLength() const = 0;
    /// @brief 判断点是否在单元内
    virtual bool contains(const Utils::Point3D& point) const = 0;

  protected:
    ElementId id_ = 0;                      ///< 单元编号
    ElementType type_ = ElementType::POINT; ///< 单元类型
    std::vector<NodePtr> nodes_;            ///< 节点指针数组
    std::vector<ElementPtr> neighbors_;     ///< 相邻单元
    MaterialId material_id_ = 0;            ///< 材料编号
};

/**
 * @class TetrahedronElement
 * @brief 四面体单元
 */
class TetrahedronElement : public Element
{
  public:
    /**
     * @brief 构造函数
     * @param id 单元编号
     * @param nodes 节点指针数组（4个）
     */
    TetrahedronElement(ElementId id, const std::vector<NodePtr>& nodes,
                       MaterialId material_id = 0);

    double getVolume() const override;
    Utils::Point3D getCentroid() const override;
    double getCharacteristicLength() const override;
    bool contains(const Utils::Point3D& point) const override;
};

/**
 * @class TriangleElement
 * @brief 三角形单元
 */
class TriangleElement : public Element
{
  public:
    /**
     * @brief 构造函数
     * @param id 单元编号
     * @param nodes 节点指针数组（3个）
     */
    TriangleElement(ElementId id, const std::vector<NodePtr>& nodes,
                    MaterialId material_id = 0);

    double getVolume() const override;
    Utils::Point3D getCentroid() const override;
    double getCharacteristicLength() const override;
    bool contains(const Utils::Point3D& point) const override;
};

/**
 * @class Mesh
 * @brief 网格对象，包含节点和单元集合。
 */
class Mesh
{
  public:
    /// @brief 添加节点
    NodePtr addNode(const Utils::Point3D& position);
    /// @brief 添加单元
    ElementPtr addElement(ElementType type, const std::vector<NodeId>& node_ids,
                          MaterialId material_id = 0);

    /// @brief 获取所有节点
    const std::vector<NodePtr>& getNodes() const
    {
        return nodes_;
    }
    /// @brief 获取所有单元
    const std::vector<ElementPtr>& getElements() const
    {
        return elements_;
    }

  private:
    std::vector<NodePtr> nodes_;       ///< 节点集合
    std::vector<ElementPtr> elements_; ///< 单元集合
};

} // namespace Mesh
} // namespace SCDAT

#endif // SCDAT_MESH_ALGORITHMS_H
