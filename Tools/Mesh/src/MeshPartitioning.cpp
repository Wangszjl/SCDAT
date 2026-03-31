#include "MeshPartitioning.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace SCDAT {
namespace Mesh {

Utils::Point3D BoundingBox::getCenter() const {
    return Utils::Point3D((min_point_.x() + max_point_.x()) * 0.5,
                          (min_point_.y() + max_point_.y()) * 0.5,
                          (min_point_.z() + max_point_.z()) * 0.5);
}

Utils::Vector3D BoundingBox::getSize() const {
    return Utils::Vector3D(max_point_.x() - min_point_.x(),
                           max_point_.y() - min_point_.y(),
                           max_point_.z() - min_point_.z());
}

bool BoundingBox::contains(const Utils::Point3D& point) const {
    return point.x() >= min_point_.x() && point.x() <= max_point_.x() &&
           point.y() >= min_point_.y() && point.y() <= max_point_.y() &&
           point.z() >= min_point_.z() && point.z() <= max_point_.z();
}

bool BoundingBox::intersects(const BoundingBox& other) const {
    return !(max_point_.x() < other.min_point_.x() || min_point_.x() > other.max_point_.x() ||
             max_point_.y() < other.min_point_.y() || min_point_.y() > other.max_point_.y() ||
             max_point_.z() < other.min_point_.z() || min_point_.z() > other.max_point_.z());
}

void BoundingBox::expand(const Utils::Point3D& point) {
    min_point_.x() = std::min(min_point_.x(), point.x());
    min_point_.y() = std::min(min_point_.y(), point.y());
    min_point_.z() = std::min(min_point_.z(), point.z());

    max_point_.x() = std::max(max_point_.x(), point.x());
    max_point_.y() = std::max(max_point_.y(), point.y());
    max_point_.z() = std::max(max_point_.z(), point.z());
}

std::array<BoundingBox, 8> BoundingBox::subdivide() const {
    const auto c = getCenter();
    return {
        BoundingBox(min_point_, c),
        BoundingBox(Utils::Point3D(c.x(), min_point_.y(), min_point_.z()), Utils::Point3D(max_point_.x(), c.y(), c.z())),
        BoundingBox(Utils::Point3D(min_point_.x(), c.y(), min_point_.z()), Utils::Point3D(c.x(), max_point_.y(), c.z())),
        BoundingBox(Utils::Point3D(c.x(), c.y(), min_point_.z()), Utils::Point3D(max_point_.x(), max_point_.y(), c.z())),
        BoundingBox(Utils::Point3D(min_point_.x(), min_point_.y(), c.z()), Utils::Point3D(c.x(), c.y(), max_point_.z())),
        BoundingBox(Utils::Point3D(c.x(), min_point_.y(), c.z()), Utils::Point3D(max_point_.x(), c.y(), max_point_.z())),
        BoundingBox(Utils::Point3D(min_point_.x(), c.y(), c.z()), Utils::Point3D(c.x(), max_point_.y(), max_point_.z())),
        BoundingBox(c, max_point_)
    };
}

std::string BoundingBox::toString() const {
    std::ostringstream oss;
    oss << "BoundingBox[(" << min_point_.x() << "," << min_point_.y() << "," << min_point_.z()
        << ")->(" << max_point_.x() << "," << max_point_.y() << "," << max_point_.z() << ")]";
    return oss.str();
}

OctreeNode::OctreeNode(const BoundingBox& bounds, int depth, int max_depth, int max_elements)
    : bounds_(bounds), depth_(depth), max_depth_(max_depth), max_elements_(max_elements) {}

void OctreeNode::addElement(ElementPtr element) {
    if (!element) return;

    if (is_leaf_) {
        elements_.push_back(element);
        if (static_cast<int>(elements_.size()) > max_elements_ && depth_ < max_depth_) {
            subdivide();
        }
        return;
    }

    const int idx = getChildIndex(element->getCentroid());
    if (children_[static_cast<std::size_t>(idx)]) {
        children_[static_cast<std::size_t>(idx)]->addElement(element);
    }
}

void OctreeNode::subdivide() {
    if (!is_leaf_) return;
    auto boxes = bounds_.subdivide();
    for (std::size_t i = 0; i < 8; ++i) {
        children_[i] = std::make_shared<OctreeNode>(boxes[i], depth_ + 1, max_depth_, max_elements_);
    }

    auto current = std::move(elements_);
    elements_.clear();
    is_leaf_ = false;

    for (const auto& e : current) {
        addElement(e);
    }
}

int OctreeNode::getChildIndex(const Utils::Point3D& point) const {
    const auto c = bounds_.getCenter();
    int idx = 0;
    if (point.x() >= c.x()) idx |= 1;
    if (point.y() >= c.y()) idx |= 2;
    if (point.z() >= c.z()) idx |= 4;
    return idx;
}

void OctreeNode::clear() {
    elements_.clear();
    for (auto& child : children_) {
        if (child) child->clear();
        child.reset();
    }
    is_leaf_ = true;
}

std::vector<ElementPtr> OctreeNode::query(const Utils::Point3D& point) const {
    if (!bounds_.contains(point)) return {};
    if (is_leaf_) return elements_;

    std::vector<ElementPtr> out;
    for (const auto& child : children_) {
        if (!child) continue;
        auto sub = child->query(point);
        out.insert(out.end(), sub.begin(), sub.end());
    }
    return out;
}

std::vector<ElementPtr> OctreeNode::query(const BoundingBox& region) const {
    if (!bounds_.intersects(region)) return {};
    if (is_leaf_) return elements_;

    std::vector<ElementPtr> out;
    for (const auto& child : children_) {
        if (!child) continue;
        auto sub = child->query(region);
        out.insert(out.end(), sub.begin(), sub.end());
    }
    return out;
}

std::vector<ElementPtr> OctreeNode::queryRadius(const Utils::Point3D& center, double radius) const {
    const auto candidates = query(center);
    std::vector<ElementPtr> out;
    out.reserve(candidates.size());
    for (const auto& e : candidates) {
        if (center.distanceTo(e->getCentroid()) <= radius) {
            out.push_back(e);
        }
    }
    return out;
}

int OctreeNode::getElementCount() const {
    if (is_leaf_) return static_cast<int>(elements_.size());
    int count = 0;
    for (const auto& child : children_) {
        if (child) count += child->getElementCount();
    }
    return count;
}

int OctreeNode::getNodeCount() const {
    int count = 1;
    if (!is_leaf_) {
        for (const auto& child : children_) {
            if (child) count += child->getNodeCount();
        }
    }
    return count;
}

int OctreeNode::getMaxDepth() const {
    int d = depth_;
    if (!is_leaf_) {
        for (const auto& child : children_) {
            if (child) d = std::max(d, child->getMaxDepth());
        }
    }
    return d;
}

SpatialIndex::SpatialIndex(int max_depth, int max_elements_per_node)
    : max_depth_(max_depth), max_elements_per_node_(max_elements_per_node) {}

void SpatialIndex::initialize(const std::vector<ElementPtr>& elements) {
    if (elements.empty()) {
        clear();
        return;
    }
    initialize(calculateBounds(elements));
    for (const auto& e : elements) {
        addElement(e);
    }
}

void SpatialIndex::initialize(const BoundingBox& bounds) {
    root_ = std::make_shared<OctreeNode>(bounds, 0, max_depth_, max_elements_per_node_);
    total_elements_ = 0;
}

void SpatialIndex::addElement(ElementPtr element) {
    if (!element) return;
    if (!root_) {
        const auto c = element->getCentroid();
        initialize(BoundingBox(c, c));
    }
    root_->addElement(element);
    ++total_elements_;
}

void SpatialIndex::clear() {
    if (root_) root_->clear();
    root_.reset();
    total_elements_ = 0;
}

std::vector<ElementPtr> SpatialIndex::findElementsContaining(const Utils::Point3D& point) const {
    if (!root_) return {};
    return root_->query(point);
}

std::vector<ElementPtr> SpatialIndex::findElementsInRegion(const BoundingBox& region) const {
    if (!root_) return {};
    return root_->query(region);
}

std::vector<ElementPtr> SpatialIndex::findElementsInRadius(const Utils::Point3D& center, double radius) const {
    if (!root_) return {};
    return root_->queryRadius(center, radius);
}

int SpatialIndex::getNodeCount() const {
    return root_ ? root_->getNodeCount() : 0;
}

int SpatialIndex::getMaxDepth() const {
    return root_ ? root_->getMaxDepth() : 0;
}

BoundingBox SpatialIndex::calculateBounds(const std::vector<ElementPtr>& elements) const {
    const auto first = elements.front()->getCentroid();
    BoundingBox bounds(first, first);
    for (const auto& e : elements) {
        bounds.expand(e->getCentroid());
    }
    return bounds;
}

} // namespace Mesh
} // namespace SCDAT
