#include "MeshAlgorithms.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SCDAT
{
namespace Mesh
{

std::string Node::toString() const
{
    std::ostringstream oss;
    oss << "Node(" << id_ << ")";
    return oss.str();
}

TetrahedronElement::TetrahedronElement(ElementId id, const std::vector<NodePtr>& nodes,
                                       MaterialId material_id)
    : Element(id, ElementType::TETRAHEDRON, nodes, material_id)
{
    if (nodes_.size() != 4)
    {
        throw std::invalid_argument("TetrahedronElement requires 4 nodes");
    }
}

double TetrahedronElement::getVolume() const
{
    if (nodes_.size() != 4)
        return 0.0;

    const auto& a = nodes_[0]->getPosition();
    const auto& b = nodes_[1]->getPosition();
    const auto& c = nodes_[2]->getPosition();
    const auto& d = nodes_[3]->getPosition();

    const double ax = b.x() - a.x();
    const double ay = b.y() - a.y();
    const double az = b.z() - a.z();

    const double bx = c.x() - a.x();
    const double by = c.y() - a.y();
    const double bz = c.z() - a.z();

    const double cx = d.x() - a.x();
    const double cy = d.y() - a.y();
    const double cz = d.z() - a.z();

    const double det =
        ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
    return std::abs(det) / 6.0;
}

Utils::Point3D TetrahedronElement::getCentroid() const
{
    if (nodes_.size() != 4)
        return Utils::Point3D();
    const auto& a = nodes_[0]->getPosition();
    const auto& b = nodes_[1]->getPosition();
    const auto& c = nodes_[2]->getPosition();
    const auto& d = nodes_[3]->getPosition();
    return Utils::Point3D((a.x() + b.x() + c.x() + d.x()) / 4.0,
                          (a.y() + b.y() + c.y() + d.y()) / 4.0,
                          (a.z() + b.z() + c.z() + d.z()) / 4.0);
}

double TetrahedronElement::getCharacteristicLength() const
{
    return std::cbrt(getVolume());
}

bool TetrahedronElement::contains(const Utils::Point3D& /*point*/) const
{
    return false;
}

TriangleElement::TriangleElement(ElementId id, const std::vector<NodePtr>& nodes,
                                 MaterialId material_id)
    : Element(id, ElementType::TRIANGLE, nodes, material_id)
{
    if (nodes_.size() != 3)
    {
        throw std::invalid_argument("TriangleElement requires 3 nodes");
    }
}

double TriangleElement::getVolume() const
{
    if (nodes_.size() != 3)
        return 0.0;
    const auto& a = nodes_[0]->getPosition();
    const auto& b = nodes_[1]->getPosition();
    const auto& c = nodes_[2]->getPosition();

    const double abx = b.x() - a.x();
    const double aby = b.y() - a.y();
    const double abz = b.z() - a.z();

    const double acx = c.x() - a.x();
    const double acy = c.y() - a.y();
    const double acz = c.z() - a.z();

    const double cxp = aby * acz - abz * acy;
    const double cyp = abz * acx - abx * acz;
    const double czp = abx * acy - aby * acx;
    return 0.5 * std::sqrt(cxp * cxp + cyp * cyp + czp * czp);
}

Utils::Point3D TriangleElement::getCentroid() const
{
    if (nodes_.size() != 3)
        return Utils::Point3D();
    const auto& a = nodes_[0]->getPosition();
    const auto& b = nodes_[1]->getPosition();
    const auto& c = nodes_[2]->getPosition();
    return Utils::Point3D((a.x() + b.x() + c.x()) / 3.0, (a.y() + b.y() + c.y()) / 3.0,
                          (a.z() + b.z() + c.z()) / 3.0);
}

double TriangleElement::getCharacteristicLength() const
{
    return std::sqrt(getVolume());
}

bool TriangleElement::contains(const Utils::Point3D& /*point*/) const
{
    return false;
}

NodePtr Mesh::addNode(const Utils::Point3D& position)
{
    auto node = std::make_shared<Node>(nodes_.size(), position);
    nodes_.push_back(node);
    return node;
}

ElementPtr Mesh::addElement(ElementType type, const std::vector<NodeId>& node_ids,
                            MaterialId material_id)
{
    std::vector<NodePtr> local_nodes;
    local_nodes.reserve(node_ids.size());
    for (auto id : node_ids)
    {
        if (id < nodes_.size())
        {
            local_nodes.push_back(nodes_[id]);
        }
    }

    ElementPtr element;
    if (type == ElementType::TETRAHEDRON && local_nodes.size() == 4)
    {
        element = std::make_shared<TetrahedronElement>(elements_.size(), local_nodes, material_id);
    }
    else if (type == ElementType::TRIANGLE && local_nodes.size() == 3)
    {
        element = std::make_shared<TriangleElement>(elements_.size(), local_nodes, material_id);
    }

    if (element)
    {
        elements_.push_back(element);
    }
    return element;
}

} // namespace Mesh
} // namespace SCDAT
