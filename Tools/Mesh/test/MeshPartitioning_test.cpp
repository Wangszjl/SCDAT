#include "../include/MeshAlgorithms.h"
#include "../include/MeshPartitioning.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <vector>

namespace SCDAT
{
namespace Mesh
{
namespace
{

ElementPtr makeTetrahedronElement(ElementId id, const Utils::Point3D& center, double half_extent)
{
    const auto make_node = [&](NodeId node_id, double dx, double dy, double dz) {
        return std::make_shared<Node>(
            node_id,
            Utils::Point3D(center.x() + dx, center.y() + dy, center.z() + dz));
    };

    std::vector<NodePtr> nodes;
    nodes.reserve(4);
    const NodeId base_id = id * 4;
    nodes.push_back(make_node(base_id + 0, -half_extent, -half_extent, -half_extent));
    nodes.push_back(make_node(base_id + 1, half_extent, -half_extent, -half_extent));
    nodes.push_back(make_node(base_id + 2, -half_extent, half_extent, -half_extent));
    nodes.push_back(make_node(base_id + 3, -half_extent, -half_extent, half_extent));
    return std::make_shared<TetrahedronElement>(id, nodes);
}

std::vector<ElementId> collectIds(const std::vector<ElementPtr>& elements)
{
    std::vector<ElementId> ids;
    ids.reserve(elements.size());
    for (const auto& element : elements)
    {
        if (element)
        {
            ids.push_back(element->getId());
        }
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

} // namespace

TEST(MeshPartitioningTest, BoundingBoxSubdivisionCoversMeshOctants)
{
    const BoundingBox bounds(Utils::Point3D(0.0, 0.0, 0.0), Utils::Point3D(8.0, 8.0, 8.0));
    const auto boxes = bounds.subdivide();

    ASSERT_EQ(boxes.size(), 8u);
    EXPECT_TRUE(boxes[0].contains(Utils::Point3D(1.0, 1.0, 1.0)));
    EXPECT_TRUE(boxes[3].contains(Utils::Point3D(7.0, 7.0, 1.0)));
    EXPECT_TRUE(boxes[4].contains(Utils::Point3D(1.0, 1.0, 7.0)));
    EXPECT_TRUE(boxes[7].contains(Utils::Point3D(7.0, 7.0, 7.0)));
    EXPECT_TRUE(bounds.contains(boxes[0].getMinPoint()));
    EXPECT_TRUE(bounds.contains(boxes[7].getMaxPoint()));
}

TEST(MeshPartitioningTest, OctreeNodeSubdividesAndQueriesByPoint)
{
    OctreeNode root(
        BoundingBox(Utils::Point3D(0.0, 0.0, 0.0), Utils::Point3D(8.0, 8.0, 8.0)),
        0,
        4,
        1);

    const auto element_a = makeTetrahedronElement(1, Utils::Point3D(1.0, 1.0, 1.0), 0.2);
    const auto element_b = makeTetrahedronElement(2, Utils::Point3D(7.0, 1.0, 1.0), 0.2);
    const auto element_c = makeTetrahedronElement(3, Utils::Point3D(7.0, 7.0, 7.0), 0.2);

    root.addElement(element_a);
    root.addElement(element_b);
    root.addElement(element_c);

    EXPECT_GT(root.getNodeCount(), 1);
    EXPECT_GT(root.getMaxDepth(), 0);

    EXPECT_EQ(collectIds(root.query(Utils::Point3D(1.0, 1.0, 1.0))),
              std::vector<ElementId>{1});
    EXPECT_EQ(collectIds(root.query(Utils::Point3D(7.0, 1.0, 1.0))),
              std::vector<ElementId>{2});
    EXPECT_EQ(collectIds(root.query(Utils::Point3D(7.0, 7.0, 7.0))),
              std::vector<ElementId>{3});
    EXPECT_TRUE(root.query(Utils::Point3D(9.0, 9.0, 9.0)).empty());
}

TEST(MeshPartitioningTest, SpatialIndexTracksDepthAndRadiusQueries)
{
    SpatialIndex index(4, 1);
    const std::vector<ElementPtr> elements{
        makeTetrahedronElement(10, Utils::Point3D(1.0, 1.0, 1.0), 0.1),
        makeTetrahedronElement(11, Utils::Point3D(5.0, 5.0, 5.0), 0.1),
        makeTetrahedronElement(12, Utils::Point3D(7.0, 7.0, 7.0), 0.1),
    };

    index.initialize(elements);

    EXPECT_EQ(index.getElementCount(), 3);
    EXPECT_GT(index.getNodeCount(), 1);
    EXPECT_GT(index.getMaxDepth(), 0);
    EXPECT_EQ(collectIds(index.findElementsInRadius(Utils::Point3D(1.0, 1.0, 1.0), 0.5)),
              std::vector<ElementId>{10});

    index.clear();
    EXPECT_EQ(index.getElementCount(), 0);
    EXPECT_EQ(index.getNodeCount(), 0);
}

TEST(MeshPartitioningTest, SpatialIndexBuildsFromMeshElementsAndRegionQueriesMatch)
{
    Mesh mesh;
    mesh.addNode(Utils::Point3D(0.0, 0.0, 0.0));
    mesh.addNode(Utils::Point3D(1.0, 0.0, 0.0));
    mesh.addNode(Utils::Point3D(0.0, 1.0, 0.0));
    mesh.addNode(Utils::Point3D(0.0, 0.0, 1.0));
    mesh.addNode(Utils::Point3D(5.0, 5.0, 0.0));
    mesh.addNode(Utils::Point3D(6.0, 5.0, 0.0));
    mesh.addNode(Utils::Point3D(5.0, 6.0, 0.0));

    ASSERT_TRUE(mesh.addElement(ElementType::TETRAHEDRON, {0, 1, 2, 3}));
    ASSERT_TRUE(mesh.addElement(ElementType::TRIANGLE, {4, 5, 6}));

    SpatialIndex index(3, 1);
    index.initialize(mesh.getElements());

    const BoundingBox local_region(
        Utils::Point3D(-0.5, -0.5, -0.5),
        Utils::Point3D(1.5, 1.5, 1.5));
    const BoundingBox remote_region(
        Utils::Point3D(4.5, 4.5, -0.5),
        Utils::Point3D(6.5, 6.5, 0.5));

    EXPECT_EQ(collectIds(index.findElementsInRegion(local_region)),
              std::vector<ElementId>{0});
    EXPECT_EQ(collectIds(index.findElementsInRegion(remote_region)),
              std::vector<ElementId>{1});
}

} // namespace Mesh
} // namespace SCDAT
