
/**
 * @file MeshParsing.h
 * @brief 网格解析、加载、生成相关接口与数据结构。
 * @author Wang Sizhan
 * @date 2026.3.25 8:48:10
 */
#ifndef SCDAT_MESH_PARSING_H
#define SCDAT_MESH_PARSING_H

#include "Geometry/include/Matrix3x3.h"
#include "Geometry/include/Point3D.h"
#include "Geometry/include/Vector3D.h"
#include "MeshAlgorithms.h"
#include "MeshPartitioning.h"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Mesh
{

using FaceNodeIndex = std::size_t;

/**
 * @struct MeshQuality
 * @brief 网格质量统计信息。
 */
struct MeshQuality
{
    double min_volume = 0.0;             ///< 最小单元体积
    double max_volume = 0.0;             ///< 最大单元体积
    double avg_volume = 0.0;             ///< 平均单元体积
    std::size_t degenerate_elements = 0; ///< 退化单元数

    /// @brief 获取质量报告字符串
    std::string getReport() const;
};

/**
 * @struct ElementRecord
 * @brief 单元记录，包含类型、节点索引、材料编号。
 */
struct ElementRecord
{
    ElementType type = ElementType::POINT; ///< 单元类型
    std::vector<std::size_t> node_indices; ///< 节点索引
    int material_id = 0;                   ///< 材料编号
};

/**
 * @class Face
 * @brief 网格面片，包含节点索引和边界编号。
 */
class Face
{
  public:
    Face() = default;
    explicit Face(std::vector<FaceNodeIndex> nodes, int boundary_id = 0)
        : node_indices_(std::move(nodes)), boundary_id_(boundary_id)
    {
    }

    /// @brief 获取面片节点索引
    const std::vector<FaceNodeIndex>& nodeIndices() const
    {
        return node_indices_;
    }
    /// @brief 获取边界编号
    int boundaryId() const
    {
        return boundary_id_;
    }

  private:
    std::vector<FaceNodeIndex> node_indices_; ///< 节点索引
    int boundary_id_ = 0;                     ///< 边界编号
};

using FacePtr = std::shared_ptr<Face>;

/**
 * @class VolMesh
 * @brief 体网格对象，包含节点、单元、面片集合及相关操作。
 */
class VolMesh
{
  public:
    VolMesh() = default;
    virtual ~VolMesh() = default;

    VolMesh(const VolMesh&) = delete;
    VolMesh& operator=(const VolMesh&) = delete;

    /// @brief 获取节点数
    std::size_t getNodeCount() const
    {
        return nodes_.size();
    }
    /// @brief 获取单元数
    std::size_t getElementCount() const
    {
        return elements_.size();
    }
    /// @brief 获取边界面数
    std::size_t getBoundaryFaceCount() const
    {
        return boundary_faces_.size();
    }

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
    /// @brief 获取所有单元记录
    const std::vector<ElementRecord>& getElementRecords() const
    {
        return element_records_;
    }
    /// @brief 获取所有边界面
    const std::vector<FacePtr>& getBoundaryFaces() const
    {
        return boundary_faces_;
    }

    /// @brief 获取指定节点
    const NodePtr& getNode(std::size_t index) const
    {
        return nodes_.at(index);
    }
    /// @brief 获取指定单元
    const ElementPtr& getElement(std::size_t index) const
    {
        return elements_.at(index);
    }

    /// @brief 添加节点
    std::size_t addNode(const Geometry::Point3D& position);
    /// @brief 添加单元
    std::size_t addElement(ElementType type, const std::vector<std::size_t>& node_indices,
                           int material_id = 0);
    /// @brief 添加边界面
    std::size_t addBoundaryFace(const std::vector<std::size_t>& node_indices, int boundary_id = 0);

    /// @brief 获取单元面积
    double getElementArea(std::size_t element_index) const;
    /// @brief 获取单元体积
    double getElementVolume(std::size_t element_index) const;
    /// @brief 获取单元法向量
    Geometry::Vector3D getElementNormal(std::size_t element_index) const;
    /// @brief 获取单元内随机点
    Geometry::Point3D getRandomPointInElement(std::size_t element_index) const;

    /// @brief 校验网格有效性
    bool validate() const;
    /// @brief 清空网格
    void clear();

    /// @brief 获取包围盒最小点
    Geometry::Point3D getBoundingBoxMin() const;
    /// @brief 获取包围盒最大点
    Geometry::Point3D getBoundingBoxMax() const;

  private:
    std::vector<NodePtr> nodes_;                 ///< 节点集合
    std::vector<ElementPtr> elements_;           ///< 单元集合
    std::vector<ElementRecord> element_records_; ///< 单元记录
    std::vector<FacePtr> boundary_faces_;        ///< 边界面集合
};

using VolMeshPtr = std::shared_ptr<VolMesh>;
using UnstructuredMesh = VolMesh;
using UnstructuredMeshPtr = VolMeshPtr;

enum class MeshFormat
{
    AUTO = 0,
    GMSH = 1,
    VTK = 2,
    STL = 3,
    OBJ = 4,
    UNKNOWN = 99
};

/**
 * @struct MeshLoadOptions
 * @brief 网格加载选项。
 */
struct MeshLoadOptions
{
    bool validate_mesh = true;          ///< 是否校验网格
    bool optimize_memory = true;        ///< 是否优化内存
    bool merge_duplicate_nodes = false; ///< 是否合并重复节点
    double merge_tolerance = 1e-12;     ///< 合并容差
    double scale_factor = 1.0;          ///< 缩放因子
    std::vector<int> material_filter;   ///< 材料过滤
    bool force_format = false;          ///< 是否强制指定格式
};

/**
 * @struct MeshLoadStatistics
 * @brief 网格加载统计信息。
 */
struct MeshLoadStatistics
{
    std::size_t node_count = 0;     ///< 节点数
    std::size_t element_count = 0;  ///< 单元数
    double load_time_seconds = 0.0; ///< 加载耗时

    /// @brief 转字符串
    std::string toString() const;
};

/**
 * @class MeshLoader
 * @brief 网格加载器基类，支持多格式解析。
 */
class MeshLoader
{
  public:
    virtual ~MeshLoader() = default;

    /// @brief 加载网格（自动识别格式）
    static std::unique_ptr<VolMesh> loadMesh(const std::string& filename,
                                             const MeshLoadOptions& options = MeshLoadOptions{});
    /// @brief 加载网格（指定格式）
    static std::unique_ptr<VolMesh> loadMesh(const std::string& filename, MeshFormat format,
                                             const MeshLoadOptions& options = MeshLoadOptions{});
    /// @brief 检测文件格式
    static MeshFormat detectFormat(const std::string& filename);
    /// @brief 创建指定格式加载器
    static std::unique_ptr<MeshLoader> createLoader(MeshFormat format);
    /// @brief 获取支持格式列表
    static std::vector<MeshFormat> getSupportedFormats();
    /// @brief 获取格式名称
    static std::string getFormatName(MeshFormat format);
    /// @brief 获取格式扩展名
    static std::vector<std::string> getFormatExtensions(MeshFormat format);

    /// @brief 加载网格文件（虚接口）
    virtual std::unique_ptr<VolMesh> loadMeshFile(const std::string& filename,
                                                  const MeshLoadOptions& options) = 0;

    /// @brief 校验网格
    virtual bool validateMesh(const VolMesh& mesh) const;
    /// @brief 优化内存
    virtual void optimizeMemory(VolMesh& mesh) const;
    /// @brief 合并重复节点
    virtual std::size_t mergeDuplicateNodes(VolMesh& mesh, double tolerance) const;
    /// @brief 应用缩放
    virtual void applyScaling(VolMesh& mesh, double scale_factor) const;
    /// @brief 材料过滤
    virtual void filterMaterials(VolMesh& mesh, const std::vector<int>& material_filter) const;
};

/**
 * @class GmshLoader
 * @brief Gmsh 格式网格加载器。
 */
class GmshLoader : public MeshLoader
{
  public:
    std::unique_ptr<VolMesh> loadMeshFile(const std::string& filename,
                                          const MeshLoadOptions& options) override;
};

/**
 * @class VTKLoader
 * @brief VTK 格式网格加载器。
 */
class VTKLoader : public MeshLoader
{
  public:
    std::unique_ptr<VolMesh> loadMeshFile(const std::string& filename,
                                          const MeshLoadOptions& options) override;
};

/**
 * @class STLLoader
 * @brief STL 格式网格加载器。
 */
class STLLoader : public MeshLoader
{
  public:
    std::unique_ptr<VolMesh> loadMeshFile(const std::string& filename,
                                          const MeshLoadOptions& options) override;
};

/**
 * @class AdvancedMeshLoader
 * @brief 高级网格加载器，支持批量加载与统计。
 */
class AdvancedMeshLoader
{
  public:
    explicit AdvancedMeshLoader(const MeshLoadOptions& options = MeshLoadOptions{});

    /// @brief 加载单个网格
    std::unique_ptr<VolMesh> loadMesh(const std::string& filename);
    /// @brief 按指定格式加载网格
    std::unique_ptr<VolMesh> loadMeshFromSoftware(const std::string& filename, MeshFormat format);
    /// @brief 批量加载网格
    std::vector<std::unique_ptr<VolMesh>> loadMeshes(const std::vector<std::string>& filenames);

    /// @brief 设置加载选项
    void setOptions(const MeshLoadOptions& options);
    /// @brief 获取加载选项
    const MeshLoadOptions& getOptions() const;

    /// @brief 获取统计信息
    const MeshLoadStatistics& getStatistics() const;
    /// @brief 获取统计字符串
    std::string getStatisticsString() const;

  private:
    MeshLoadOptions options_;       ///< 加载选项
    MeshLoadStatistics statistics_; ///< 统计信息
};

enum class GeometryShape
{
    BOX = 0,
    CYLINDER = 1,
    SPHERE = 2,
    PLATE = 3
};

/**
 * @struct GeometryDimensions
 * @brief 几何体尺寸参数。
 */
struct GeometryDimensions
{
    GeometryShape shape = GeometryShape::BOX; ///< 形状类型
    double length = 1.0;                      ///< 长度
    double width = 1.0;                       ///< 宽度
    double height = 1.0;                      ///< 高度
    double radius = 0.5;                      ///< 半径
    double thickness = 0.1;                   ///< 厚度
};

/**
 * @struct MeshGenerationOptions
 * @brief 网格生成参数。
 */
struct MeshGenerationOptions
{
    std::size_t nx = 10;              ///< X方向分段数
    std::size_t ny = 10;              ///< Y方向分段数
    std::size_t nz = 10;              ///< Z方向分段数
    std::size_t radial_segments = 24; ///< 径向分段数
    std::size_t height_segments = 10; ///< 高度分段数
    bool tetrahedralize = true;       ///< 是否四面体化
};

/**
 * @class GeometryMeshGenerator
 * @brief 几何体网格生成器，支持多种形状。
 */
class GeometryMeshGenerator
{
  public:
    /// @brief 根据尺寸参数生成网格
    static std::unique_ptr<VolMesh>
    generateFromDimensions(const GeometryDimensions& dims,
                           const MeshGenerationOptions& options = MeshGenerationOptions{});

  private:
    static std::unique_ptr<VolMesh> generateBoxMesh(const GeometryDimensions& dims,
                                                    const MeshGenerationOptions& options);
    static std::unique_ptr<VolMesh> generateCylinderMesh(const GeometryDimensions& dims,
                                                         const MeshGenerationOptions& options);
    static std::unique_ptr<VolMesh> generateSphereMesh(const GeometryDimensions& dims,
                                                       const MeshGenerationOptions& options);
    static std::unique_ptr<VolMesh> generatePlateMesh(const GeometryDimensions& dims,
                                                      const MeshGenerationOptions& options);
};

} // namespace Mesh

/**
 * @class SurfaceElement
 * @brief 表面单元，包含中心、法向、面积。
 */
class SurfaceElement
{
  public:
    SurfaceElement() = default;

    /// @brief 设置中心点
    void setCenter(const Geometry::Point3D& center)
    {
        center_ = center;
    }
    /// @brief 设置法向量
    void setNormal(const Geometry::Vector3D& normal)
    {
        normal_ = normal;
    }
    /// @brief 设置面积
    void setArea(double area)
    {
        area_ = area;
    }

    /// @brief 获取中心点
    const Geometry::Point3D& getCenter() const
    {
        return center_;
    }
    /// @brief 获取法向量
    const Geometry::Vector3D& getNormal() const
    {
        return normal_;
    }
    /// @brief 获取面积
    double getArea() const
    {
        return area_;
    }

  private:
    Geometry::Point3D center_;  ///< 中心点
    Geometry::Vector3D normal_; ///< 法向量
    double area_ = 0.0;         ///< 面积
};

} // namespace SCDAT

#endif // SCDAT_MESH_PARSING_H
