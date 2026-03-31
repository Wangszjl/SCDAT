#include "MeshParsing.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace SCDAT
{
namespace Mesh
{

namespace
{

constexpr double kEps = 1e-12;

Geometry::Vector3D toVec(const Geometry::Point3D& p)
{
    return Geometry::Vector3D(p.x(), p.y(), p.z());
}

Geometry::Vector3D cross3(const Geometry::Vector3D& a, const Geometry::Vector3D& b)
{
    return Geometry::Vector3D(a.y() * b.z() - a.z() * b.y(), a.z() * b.x() - a.x() * b.z(),
                              a.x() * b.y() - a.y() * b.x());
}

double dot3(const Geometry::Vector3D& a, const Geometry::Vector3D& b)
{
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

double norm3(const Geometry::Vector3D& v)
{
    return std::sqrt(dot3(v, v));
}

Geometry::Vector3D normalize3(const Geometry::Vector3D& v)
{
    const double n = norm3(v);
    if (n <= kEps)
    {
        return Geometry::Vector3D(0.0, 0.0, 1.0);
    }
    return Geometry::Vector3D(v.x() / n, v.y() / n, v.z() / n);
}

double triangleArea(const Geometry::Point3D& a, const Geometry::Point3D& b,
                    const Geometry::Point3D& c)
{
    const Geometry::Vector3D ab(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
    const Geometry::Vector3D ac(c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
    return 0.5 * norm3(cross3(ab, ac));
}

double tetraVolume(const Geometry::Point3D& a, const Geometry::Point3D& b,
                   const Geometry::Point3D& c, const Geometry::Point3D& d)
{
    const Geometry::Vector3D ab(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
    const Geometry::Vector3D ac(c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
    const Geometry::Vector3D ad(d.x() - a.x(), d.y() - a.y(), d.z() - a.z());
    return std::abs(dot3(ab, cross3(ac, ad))) / 6.0;
}

Geometry::Point3D triangleCentroid(const Geometry::Point3D& a, const Geometry::Point3D& b,
                                   const Geometry::Point3D& c)
{
    return Geometry::Point3D((a.x() + b.x() + c.x()) / 3.0, (a.y() + b.y() + c.y()) / 3.0,
                             (a.z() + b.z() + c.z()) / 3.0);
}

Geometry::Point3D tetraCentroid(const Geometry::Point3D& a, const Geometry::Point3D& b,
                                const Geometry::Point3D& c, const Geometry::Point3D& d)
{
    return Geometry::Point3D((a.x() + b.x() + c.x() + d.x()) / 4.0,
                             (a.y() + b.y() + c.y() + d.y()) / 4.0,
                             (a.z() + b.z() + c.z() + d.z()) / 4.0);
}

std::mt19937& rng()
{
    thread_local std::mt19937 gen(std::random_device{}());
    return gen;
}

Geometry::Point3D randomPointInTriangle(const Geometry::Point3D& a, const Geometry::Point3D& b,
                                        const Geometry::Point3D& c)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double r1 = std::sqrt(dist(rng()));
    const double r2 = dist(rng());

    const double w0 = 1.0 - r1;
    const double w1 = r1 * (1.0 - r2);
    const double w2 = r1 * r2;

    return Geometry::Point3D(w0 * a.x() + w1 * b.x() + w2 * c.x(),
                             w0 * a.y() + w1 * b.y() + w2 * c.y(),
                             w0 * a.z() + w1 * b.z() + w2 * c.z());
}

Geometry::Point3D randomPointInTetra(const Geometry::Point3D& a, const Geometry::Point3D& b,
                                     const Geometry::Point3D& c, const Geometry::Point3D& d)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double e0 = -std::log(std::max(kEps, dist(rng())));
    const double e1 = -std::log(std::max(kEps, dist(rng())));
    const double e2 = -std::log(std::max(kEps, dist(rng())));
    const double e3 = -std::log(std::max(kEps, dist(rng())));
    const double s = e0 + e1 + e2 + e3;

    const double w0 = e0 / s;
    const double w1 = e1 / s;
    const double w2 = e2 / s;
    const double w3 = e3 / s;

    return Geometry::Point3D(w0 * a.x() + w1 * b.x() + w2 * c.x() + w3 * d.x(),
                             w0 * a.y() + w1 * b.y() + w2 * c.y() + w3 * d.y(),
                             w0 * a.z() + w1 * b.z() + w2 * c.z() + w3 * d.z());
}
// 把三维点按 tol 精度映射为唯一的字符串 key，用于空间哈希或去重
std::string keyFromPoint(const Geometry::Point3D& p, double tol)
{
    const auto qx = static_cast<long long>(std::llround(p.x() / tol));
    const auto qy = static_cast<long long>(std::llround(p.y() / tol));
    const auto qz = static_cast<long long>(std::llround(p.z() / tol));
    return std::to_string(qx) + ":" + std::to_string(qy) + ":" + std::to_string(qz);
}

std::vector<std::size_t> toSizeIndices(const std::vector<int>& ids)
{
    std::vector<std::size_t> out;
    out.reserve(ids.size());
    for (int id : ids)
    {
        if (id < 0)
        {
            throw std::runtime_error("Negative node index in mesh element");
        }
        out.push_back(static_cast<std::size_t>(id));
    }
    return out;
}

bool startsWith(const std::string& s, const std::string& p)
{
    return s.rfind(p, 0) == 0;
}

struct PointKey
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    bool operator==(const PointKey& other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct PointKeyHash
{
    std::size_t operator()(const PointKey& key) const
    {
        const std::size_t hx = std::hash<double>{}(key.x);
        const std::size_t hy = std::hash<double>{}(key.y);
        const std::size_t hz = std::hash<double>{}(key.z);
        return hx ^ (hy + 0x9e3779b9 + (hx << 6U) + (hx >> 2U)) ^
               (hz + 0x9e3779b9 + (hy << 6U) + (hy >> 2U));
    }
};

} // namespace

std::string MeshQuality::getReport() const
{
    std::ostringstream oss;
    oss << "min_volume=" << min_volume << ", max_volume=" << max_volume
        << ", avg_volume=" << avg_volume << ", degenerate=" << degenerate_elements;
    return oss.str();
}

std::size_t VolMesh::addNode(const Geometry::Point3D& position)
{
    const std::size_t id = nodes_.size();
    nodes_.push_back(std::make_shared<Node>(id, position));
    return id;
}

std::size_t VolMesh::addElement(ElementType type, const std::vector<std::size_t>& node_indices,
                                int material_id)
{
    std::vector<NodePtr> local_nodes;
    local_nodes.reserve(node_indices.size());

    for (std::size_t idx : node_indices)
    {
        if (idx >= nodes_.size())
        {
            throw std::out_of_range("VolMesh::addElement node index out of range");
        }
        local_nodes.push_back(nodes_[idx]);
    }

    ElementPtr element;
    if (type == ElementType::TRIANGLE)
    {
        if (local_nodes.size() != 3)
        {
            throw std::invalid_argument("Triangle element needs exactly 3 nodes");
        }
        element = std::make_shared<TriangleElement>(elements_.size(), local_nodes,
                                                    static_cast<MaterialId>(material_id));
    }
    else if (type == ElementType::TETRAHEDRON)
    {
        if (local_nodes.size() != 4)
        {
            throw std::invalid_argument("Tetrahedron element needs exactly 4 nodes");
        }
        element = std::make_shared<TetrahedronElement>(elements_.size(), local_nodes,
                                                       static_cast<MaterialId>(material_id));
    }
    else
    {
        throw std::invalid_argument("Unsupported element type in VolMesh::addElement");
    }

    elements_.push_back(element);
    element_records_.push_back(ElementRecord{type, node_indices, material_id});
    return elements_.size() - 1;
}

std::size_t VolMesh::addBoundaryFace(const std::vector<std::size_t>& node_indices, int boundary_id)
{
    for (std::size_t idx : node_indices)
    {
        if (idx >= nodes_.size())
        {
            throw std::out_of_range("VolMesh::addBoundaryFace node index out of range");
        }
    }
    boundary_faces_.push_back(std::make_shared<Face>(node_indices, boundary_id));
    return boundary_faces_.size() - 1;
}

double VolMesh::getElementArea(std::size_t element_index) const
{
    if (element_index >= element_records_.size())
    {
        throw std::out_of_range("VolMesh::getElementArea index out of range");
    }

    const auto& rec = element_records_[element_index];
    if (rec.type == ElementType::TRIANGLE && rec.node_indices.size() == 3)
    {
        const auto& a = nodes_[rec.node_indices[0]]->getPosition();
        const auto& b = nodes_[rec.node_indices[1]]->getPosition();
        const auto& c = nodes_[rec.node_indices[2]]->getPosition();
        return triangleArea(a, b, c);
    }

    if (rec.type == ElementType::TETRAHEDRON && rec.node_indices.size() == 4)
    {
        const auto& a = nodes_[rec.node_indices[0]]->getPosition();
        const auto& b = nodes_[rec.node_indices[1]]->getPosition();
        const auto& c = nodes_[rec.node_indices[2]]->getPosition();
        const auto& d = nodes_[rec.node_indices[3]]->getPosition();

        return triangleArea(a, b, c) + triangleArea(a, b, d) + triangleArea(a, c, d) +
               triangleArea(b, c, d);
    }

    return 0.0;
}

double VolMesh::getElementVolume(std::size_t element_index) const
{
    if (element_index >= element_records_.size())
    {
        throw std::out_of_range("VolMesh::getElementVolume index out of range");
    }

    const auto& rec = element_records_[element_index];
    if (rec.type == ElementType::TETRAHEDRON && rec.node_indices.size() == 4)
    {
        const auto& a = nodes_[rec.node_indices[0]]->getPosition();
        const auto& b = nodes_[rec.node_indices[1]]->getPosition();
        const auto& c = nodes_[rec.node_indices[2]]->getPosition();
        const auto& d = nodes_[rec.node_indices[3]]->getPosition();
        return tetraVolume(a, b, c, d);
    }

    return 0.0;
}

Geometry::Vector3D VolMesh::getElementNormal(std::size_t element_index) const
{
    if (element_index >= element_records_.size())
    {
        throw std::out_of_range("VolMesh::getElementNormal index out of range");
    }

    const auto& rec = element_records_[element_index];
    if (rec.node_indices.size() < 3)
    {
        return Geometry::Vector3D(0.0, 0.0, 1.0);
    }

    const auto& a = nodes_[rec.node_indices[0]]->getPosition();
    const auto& b = nodes_[rec.node_indices[1]]->getPosition();
    const auto& c = nodes_[rec.node_indices[2]]->getPosition();

    const Geometry::Vector3D ab(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
    const Geometry::Vector3D ac(c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
    return normalize3(cross3(ab, ac));
}

Geometry::Point3D VolMesh::getRandomPointInElement(std::size_t element_index) const
{
    if (element_index >= element_records_.size())
    {
        throw std::out_of_range("VolMesh::getRandomPointInElement index out of range");
    }

    const auto& rec = element_records_[element_index];
    if (rec.type == ElementType::TRIANGLE && rec.node_indices.size() == 3)
    {
        const auto& a = nodes_[rec.node_indices[0]]->getPosition();
        const auto& b = nodes_[rec.node_indices[1]]->getPosition();
        const auto& c = nodes_[rec.node_indices[2]]->getPosition();
        return randomPointInTriangle(a, b, c);
    }

    if (rec.type == ElementType::TETRAHEDRON && rec.node_indices.size() == 4)
    {
        const auto& a = nodes_[rec.node_indices[0]]->getPosition();
        const auto& b = nodes_[rec.node_indices[1]]->getPosition();
        const auto& c = nodes_[rec.node_indices[2]]->getPosition();
        const auto& d = nodes_[rec.node_indices[3]]->getPosition();
        return randomPointInTetra(a, b, c, d);
    }

    return Geometry::Point3D(0.0, 0.0, 0.0);
}

bool VolMesh::validate() const
{
    if (nodes_.empty())
    {
        return false;
    }

    if (elements_.size() != element_records_.size())
    {
        return false;
    }

    for (const auto& rec : element_records_)
    {
        if (rec.type == ElementType::TRIANGLE && rec.node_indices.size() != 3)
        {
            return false;
        }
        if (rec.type == ElementType::TETRAHEDRON && rec.node_indices.size() != 4)
        {
            return false;
        }
        for (std::size_t idx : rec.node_indices)
        {
            if (idx >= nodes_.size())
            {
                return false;
            }
        }
    }

    return true;
}

void VolMesh::clear()
{
    nodes_.clear();
    elements_.clear();
    element_records_.clear();
    boundary_faces_.clear();
}

Geometry::Point3D VolMesh::getBoundingBoxMin() const
{
    if (nodes_.empty())
    {
        return Geometry::Point3D(0.0, 0.0, 0.0);
    }

    Geometry::Point3D p = nodes_.front()->getPosition();
    for (const auto& n : nodes_)
    {
        p.x() = std::min(p.x(), n->getPosition().x());
        p.y() = std::min(p.y(), n->getPosition().y());
        p.z() = std::min(p.z(), n->getPosition().z());
    }
    return p;
}

Geometry::Point3D VolMesh::getBoundingBoxMax() const
{
    if (nodes_.empty())
    {
        return Geometry::Point3D(0.0, 0.0, 0.0);
    }

    Geometry::Point3D p = nodes_.front()->getPosition();
    for (const auto& n : nodes_)
    {
        p.x() = std::max(p.x(), n->getPosition().x());
        p.y() = std::max(p.y(), n->getPosition().y());
        p.z() = std::max(p.z(), n->getPosition().z());
    }
    return p;
}

std::string MeshLoadStatistics::toString() const
{
    std::ostringstream oss;
    oss << "nodes=" << node_count << ", elements=" << element_count
        << ", load_time_seconds=" << load_time_seconds;
    return oss.str();
}

std::unique_ptr<VolMesh> MeshLoader::loadMesh(const std::string& filename,
                                              const MeshLoadOptions& options)
{
    const MeshFormat format = detectFormat(filename);
    return loadMesh(filename, format, options);
}

std::unique_ptr<VolMesh> MeshLoader::loadMesh(const std::string& filename, MeshFormat format,
                                              const MeshLoadOptions& options)
{
    auto loader = createLoader(format);
    auto mesh = loader->loadMeshFile(filename, options);

    if (!mesh)
    {
        throw std::runtime_error("MeshLoader returned null mesh for file: " + filename);
    }

    if (options.validate_mesh && !loader->validateMesh(*mesh))
    {
        throw std::runtime_error("Loaded mesh validation failed: " + filename);
    }

    return mesh;
}

MeshFormat MeshLoader::detectFormat(const std::string& filename)
{
    const auto dot_pos = filename.find_last_of('.');
    if (dot_pos == std::string::npos)
    {
        return MeshFormat::UNKNOWN;
    }

    std::string ext = filename.substr(dot_pos);
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (ext == ".msh")
        return MeshFormat::GMSH;
    if (ext == ".vtk")
        return MeshFormat::VTK;
    if (ext == ".stl")
        return MeshFormat::STL;
    if (ext == ".obj")
        return MeshFormat::OBJ;
    return MeshFormat::UNKNOWN;
}

std::unique_ptr<MeshLoader> MeshLoader::createLoader(MeshFormat format)
{
    switch (format)
    {
    case MeshFormat::GMSH:
        return std::make_unique<GmshLoader>();
    case MeshFormat::VTK:
        return std::make_unique<VTKLoader>();
    case MeshFormat::STL:
        return std::make_unique<STLLoader>();
    case MeshFormat::OBJ:
    case MeshFormat::AUTO:
    case MeshFormat::UNKNOWN:
    default:
        throw std::runtime_error("Unsupported mesh format for parser creation");
    }
}

std::vector<MeshFormat> MeshLoader::getSupportedFormats()
{
    return {MeshFormat::GMSH, MeshFormat::VTK, MeshFormat::STL};
}

std::string MeshLoader::getFormatName(MeshFormat format)
{
    switch (format)
    {
    case MeshFormat::GMSH:
        return "Gmsh";
    case MeshFormat::VTK:
        return "VTK";
    case MeshFormat::STL:
        return "STL";
    case MeshFormat::OBJ:
        return "OBJ";
    case MeshFormat::AUTO:
        return "AUTO";
    default:
        return "Unknown";
    }
}

std::vector<std::string> MeshLoader::getFormatExtensions(MeshFormat format)
{
    switch (format)
    {
    case MeshFormat::GMSH:
        return {".msh"};
    case MeshFormat::VTK:
        return {".vtk"};
    case MeshFormat::STL:
        return {".stl"};
    case MeshFormat::OBJ:
        return {".obj"};
    default:
        return {};
    }
}

bool MeshLoader::validateMesh(const VolMesh& mesh) const
{
    return mesh.validate();
}
void MeshLoader::optimizeMemory(VolMesh& mesh) const
{
    const auto& source_nodes = mesh.getNodes();
    const auto& source_elements = mesh.getElementRecords();
    const auto& source_faces = mesh.getBoundaryFaces();

    if (source_nodes.empty())
    {
        return;
    }

    std::vector<std::size_t> index_map(source_nodes.size(), 0);
    std::unordered_map<PointKey, std::size_t, PointKeyHash> unique_nodes;
    unique_nodes.reserve(source_nodes.size());

    VolMesh optimized_mesh;

    for (std::size_t i = 0; i < source_nodes.size(); ++i)
    {
        const auto& source_node = source_nodes[i];
        const auto& position = source_node->getPosition();
        const PointKey key{position.x(), position.y(), position.z()};

        const auto [it, inserted] = unique_nodes.emplace(key, optimized_mesh.getNodeCount());
        if (inserted)
        {
            const std::size_t new_index = optimized_mesh.addNode(position);
            auto& new_node = optimized_mesh.getNode(new_index);
            new_node->setBoundaryType(source_node->getBoundaryType());
            new_node->setPotential(source_node->getPotential());
            new_node->setElectricField(source_node->getElectricField());
            new_node->setChargeDensity(source_node->getChargeDensity());
            it->second = new_index;
        }

        index_map[i] = it->second;
    }

    for (const auto& element : source_elements)
    {
        std::vector<std::size_t> remapped_indices;
        remapped_indices.reserve(element.node_indices.size());
        for (std::size_t node_index : element.node_indices)
        {
            remapped_indices.push_back(index_map[node_index]);
        }

        optimized_mesh.addElement(element.type, remapped_indices, element.material_id);
    }

    for (const auto& face : source_faces)
    {
        std::vector<std::size_t> remapped_indices;
        remapped_indices.reserve(face->nodeIndices().size());
        for (std::size_t node_index : face->nodeIndices())
        {
            remapped_indices.push_back(index_map[node_index]);
        }

        optimized_mesh.addBoundaryFace(remapped_indices, face->boundaryId());
    }

    mesh.clear();

    for (const auto& node : optimized_mesh.getNodes())
    {
        const std::size_t new_index = mesh.addNode(node->getPosition());
        auto& target_node = mesh.getNode(new_index);
        target_node->setBoundaryType(node->getBoundaryType());
        target_node->setPotential(node->getPotential());
        target_node->setElectricField(node->getElectricField());
        target_node->setChargeDensity(node->getChargeDensity());
    }

    for (const auto& element : optimized_mesh.getElementRecords())
    {
        mesh.addElement(element.type, element.node_indices, element.material_id);
    }

    for (const auto& face : optimized_mesh.getBoundaryFaces())
    {
        mesh.addBoundaryFace(face->nodeIndices(), face->boundaryId());
    }
}

std::size_t MeshLoader::mergeDuplicateNodes(VolMesh& /*mesh*/, double /*tolerance*/) const
{
    return 0;
}

void MeshLoader::applyScaling(VolMesh& /*mesh*/, double /*scale_factor*/) const {}

void MeshLoader::filterMaterials(VolMesh& /*mesh*/,
                                 const std::vector<int>& /*material_filter*/) const
{
}

namespace
{

std::unique_ptr<VolMesh> parseGmshAscii(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("Failed to open Gmsh file: " + filename);
    }

    auto mesh = std::make_unique<VolMesh>();
    std::unordered_map<int, std::size_t> node_id_to_index;
    std::unordered_map<int, int> surface_entity_to_physical;
    std::unordered_map<int, int> volume_entity_to_physical;
    std::string line;

    auto parseIntList = [](const std::string& text)
    {
        std::istringstream iss(text);
        std::vector<long long> values;
        long long v = 0;
        while (iss >> v)
        {
            values.push_back(v);
        }
        return values;
    };

    while (std::getline(in, line))
    {
        if (line == "$Entities")
        {
            std::getline(in, line);
            const auto entity_header = parseIntList(line);
            if (entity_header.size() < 4)
            {
                throw std::runtime_error("Invalid Gmsh $Entities header");
            }

            const int num_points = static_cast<int>(entity_header[0]);
            const int num_curves = static_cast<int>(entity_header[1]);
            const int num_surfaces = static_cast<int>(entity_header[2]);
            const int num_volumes = static_cast<int>(entity_header[3]);

            for (int i = 0; i < num_points; ++i)
            {
                std::getline(in, line);
            }
            for (int i = 0; i < num_curves; ++i)
            {
                std::getline(in, line);
            }

            for (int i = 0; i < num_surfaces; ++i)
            {
                std::getline(in, line);
                std::istringstream iss(line);

                int tag = 0;
                double min_x = 0.0, min_y = 0.0, min_z = 0.0;
                double max_x = 0.0, max_y = 0.0, max_z = 0.0;
                int num_physical_tags = 0;
                iss >> tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z >>
                    num_physical_tags;

                int first_physical_tag = 0;
                for (int p = 0; p < num_physical_tags; ++p)
                {
                    int phys = 0;
                    iss >> phys;
                    if (p == 0)
                    {
                        first_physical_tag = phys;
                    }
                }

                if (num_physical_tags > 0)
                {
                    surface_entity_to_physical[tag] = first_physical_tag;
                }
            }

            for (int i = 0; i < num_volumes; ++i)
            {
                std::getline(in, line);
                std::istringstream iss(line);

                int tag = 0;
                double min_x = 0.0, min_y = 0.0, min_z = 0.0;
                double max_x = 0.0, max_y = 0.0, max_z = 0.0;
                int num_physical_tags = 0;
                iss >> tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z >>
                    num_physical_tags;

                int first_physical_tag = 0;
                for (int p = 0; p < num_physical_tags; ++p)
                {
                    int phys = 0;
                    iss >> phys;
                    if (p == 0)
                    {
                        first_physical_tag = phys;
                    }
                }

                if (num_physical_tags > 0)
                {
                    volume_entity_to_physical[tag] = first_physical_tag;
                }
            }
        }
        else if (line == "$Nodes")
        {
            std::getline(in, line);
            const auto header_vals = parseIntList(line);

            if (header_vals.size() == 1)
            {
                // Gmsh 2.x: single node count line.
                const int n = static_cast<int>(header_vals[0]);
                for (int i = 0; i < n; ++i)
                {
                    std::getline(in, line);
                    std::istringstream iss(line);
                    int id = 0;
                    double x = 0.0, y = 0.0, z = 0.0;
                    iss >> id >> x >> y >> z;
                    node_id_to_index[id] = mesh->addNode(Geometry::Point3D(x, y, z));
                }
            }
            else if (header_vals.size() >= 4)
            {
                // Gmsh 4.x: block-based node sections.
                const int num_entity_blocks = static_cast<int>(header_vals[0]);

                for (int b = 0; b < num_entity_blocks; ++b)
                {
                    std::getline(in, line);
                    const auto block_header = parseIntList(line);
                    if (block_header.size() < 4)
                    {
                        throw std::runtime_error("Invalid Gmsh 4.x node block header");
                    }

                    const int parametric = static_cast<int>(block_header[2]);
                    const int num_nodes_in_block = static_cast<int>(block_header[3]);

                    std::vector<int> node_tags;
                    node_tags.reserve(static_cast<std::size_t>(num_nodes_in_block));
                    while (static_cast<int>(node_tags.size()) < num_nodes_in_block)
                    {
                        std::getline(in, line);
                        std::istringstream id_stream(line);
                        int node_tag = 0;
                        while (id_stream >> node_tag)
                        {
                            node_tags.push_back(node_tag);
                            if (static_cast<int>(node_tags.size()) == num_nodes_in_block)
                            {
                                break;
                            }
                        }
                    }

                    for (int i = 0; i < num_nodes_in_block; ++i)
                    {
                        std::getline(in, line);
                        std::istringstream coord_stream(line);
                        double x = 0.0, y = 0.0, z = 0.0;
                        coord_stream >> x >> y >> z;

                        if (parametric == 1)
                        {
                            // Ignore optional parametric coordinates if present.
                            double dummy = 0.0;
                            while (coord_stream >> dummy)
                            {
                            }
                        }

                        node_id_to_index[node_tags[static_cast<std::size_t>(i)]] =
                            mesh->addNode(Geometry::Point3D(x, y, z));
                    }
                }
            }
            else
            {
                throw std::runtime_error("Unsupported Gmsh $Nodes header");
            }
        }
        else if (line == "$Elements")
        {
            std::getline(in, line);
            const auto header_vals = parseIntList(line);

            if (header_vals.size() == 1)
            {
                // Gmsh 2.x element section.
                const int m = static_cast<int>(header_vals[0]);
                for (int i = 0; i < m; ++i)
                {
                    std::getline(in, line);
                    std::istringstream iss(line);

                    int id = 0;
                    int gmsh_type = 0;
                    int ntags = 0;
                    iss >> id >> gmsh_type >> ntags;

                    int physical_tag = 0;
                    for (int t = 0; t < ntags; ++t)
                    {
                        int dummy = 0;
                        iss >> dummy;
                        if (t == 0)
                        {
                            physical_tag = dummy;
                        }
                    }

                    if (gmsh_type == 2)
                    {
                        int n1 = 0, n2 = 0, n3 = 0;
                        iss >> n1 >> n2 >> n3;
                        mesh->addElement(ElementType::TRIANGLE,
                                         {node_id_to_index.at(n1), node_id_to_index.at(n2),
                                          node_id_to_index.at(n3)},
                                         physical_tag);
                    }
                    else if (gmsh_type == 4)
                    {
                        int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
                        iss >> n1 >> n2 >> n3 >> n4;
                        mesh->addElement(ElementType::TETRAHEDRON,
                                         {node_id_to_index.at(n1), node_id_to_index.at(n2),
                                          node_id_to_index.at(n3), node_id_to_index.at(n4)},
                                         physical_tag);
                    }
                }
            }
            else if (header_vals.size() >= 4)
            {
                // Gmsh 4.x element section.
                const int num_entity_blocks = static_cast<int>(header_vals[0]);

                for (int b = 0; b < num_entity_blocks; ++b)
                {
                    std::getline(in, line);
                    const auto block_header = parseIntList(line);
                    if (block_header.size() < 4)
                    {
                        throw std::runtime_error("Invalid Gmsh 4.x element block header");
                    }

                    const int entity_dim = static_cast<int>(block_header[0]);
                    const int entity_tag = static_cast<int>(block_header[1]);
                    const int gmsh_type = static_cast<int>(block_header[2]);
                    const int num_elements_in_block = static_cast<int>(block_header[3]);

                    int physical_tag = 0;
                    if (entity_dim == 2)
                    {
                        const auto it = surface_entity_to_physical.find(entity_tag);
                        if (it != surface_entity_to_physical.end())
                        {
                            physical_tag = it->second;
                        }
                    }
                    else if (entity_dim == 3)
                    {
                        const auto it = volume_entity_to_physical.find(entity_tag);
                        if (it != volume_entity_to_physical.end())
                        {
                            physical_tag = it->second;
                        }
                    }

                    for (int i = 0; i < num_elements_in_block; ++i)
                    {
                        std::getline(in, line);
                        std::istringstream iss(line);

                        int element_tag = 0;
                        iss >> element_tag;

                        if (gmsh_type == 2)
                        {
                            int n1 = 0, n2 = 0, n3 = 0;
                            iss >> n1 >> n2 >> n3;
                            mesh->addElement(ElementType::TRIANGLE,
                                             {node_id_to_index.at(n1), node_id_to_index.at(n2),
                                              node_id_to_index.at(n3)},
                                             physical_tag);
                        }
                        else if (gmsh_type == 4)
                        {
                            int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
                            iss >> n1 >> n2 >> n3 >> n4;
                            mesh->addElement(ElementType::TETRAHEDRON,
                                             {node_id_to_index.at(n1), node_id_to_index.at(n2),
                                              node_id_to_index.at(n3), node_id_to_index.at(n4)},
                                             physical_tag);
                        }
                    }
                }
            }
            else
            {
                throw std::runtime_error("Unsupported Gmsh $Elements header");
            }
        }
    }

    return mesh;
}

std::unique_ptr<VolMesh> parseVtkLegacy(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("Failed to open VTK file: " + filename);
    }

    auto mesh = std::make_unique<VolMesh>();
    std::string token;
    std::vector<std::vector<std::size_t>> cells;

    while (in >> token)
    {
        if (token == "POINTS")
        {
            int n = 0;
            std::string type;
            in >> n >> type;
            for (int i = 0; i < n; ++i)
            {
                double x = 0.0, y = 0.0, z = 0.0;
                in >> x >> y >> z;
                mesh->addNode(Geometry::Point3D(x, y, z));
            }
        }
        else if (token == "CELLS")
        {
            int n_cells = 0;
            int total = 0;
            in >> n_cells >> total;
            cells.resize(static_cast<std::size_t>(n_cells));
            for (int i = 0; i < n_cells; ++i)
            {
                int csize = 0;
                in >> csize;
                cells[static_cast<std::size_t>(i)].resize(static_cast<std::size_t>(csize));
                for (int j = 0; j < csize; ++j)
                {
                    int idx = 0;
                    in >> idx;
                    cells[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                        static_cast<std::size_t>(idx);
                }
            }
        }
        else if (token == "CELL_TYPES")
        {
            int n_types = 0;
            in >> n_types;
            for (int i = 0; i < n_types; ++i)
            {
                int vtk_type = 0;
                in >> vtk_type;
                const auto& cell = cells[static_cast<std::size_t>(i)];
                if (vtk_type == 5 && cell.size() == 3)
                {
                    mesh->addElement(ElementType::TRIANGLE, cell);
                }
                else if (vtk_type == 10 && cell.size() == 4)
                {
                    mesh->addElement(ElementType::TETRAHEDRON, cell);
                }
            }
        }
    }

    return mesh;
}

std::unique_ptr<VolMesh> parseStlAscii(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        throw std::runtime_error("Failed to open STL file: " + filename);
    }

    std::string first;
    in >> first;
    if (first != "solid")
    {
        throw std::runtime_error("Only ASCII STL is currently supported");
    }

    in.clear();
    in.seekg(0, std::ios::beg);

    auto mesh = std::make_unique<VolMesh>();
    std::unordered_map<std::string, std::size_t> node_map;
    std::string line;
    std::vector<std::size_t> tri;

    while (std::getline(in, line))
    {
        std::istringstream iss(line);
        std::string word;
        iss >> word;
        if (word == "vertex")
        {
            double x = 0.0, y = 0.0, z = 0.0;
            iss >> x >> y >> z;
            Geometry::Point3D p(x, y, z);
            const std::string key = keyFromPoint(p, 1e-12);
            auto it = node_map.find(key);
            std::size_t idx = 0;
            if (it == node_map.end())
            {
                idx = mesh->addNode(p);
                node_map.emplace(key, idx);
            }
            else
            {
                idx = it->second;
            }

            tri.push_back(idx);
            if (tri.size() == 3)
            {
                mesh->addElement(ElementType::TRIANGLE, tri);
                mesh->addBoundaryFace(tri, 1);
                tri.clear();
            }
        }
    }

    return mesh;
}

} // namespace

std::unique_ptr<VolMesh> GmshLoader::loadMeshFile(const std::string& filename,
                                                  const MeshLoadOptions& options)
{
    auto mesh = parseGmshAscii(filename);
    if (options.validate_mesh && !mesh->validate())
    {
        throw std::runtime_error("Invalid Gmsh mesh: " + filename);
    }
    return mesh;
}

std::unique_ptr<VolMesh> VTKLoader::loadMeshFile(const std::string& filename,
                                                 const MeshLoadOptions& options)
{
    auto mesh = parseVtkLegacy(filename);
    if (options.validate_mesh && !mesh->validate())
    {
        throw std::runtime_error("Invalid VTK mesh: " + filename);
    }
    return mesh;
}

std::unique_ptr<VolMesh> STLLoader::loadMeshFile(const std::string& filename,
                                                 const MeshLoadOptions& options)
{
    auto mesh = parseStlAscii(filename);
    if (options.validate_mesh && !mesh->validate())
    {
        throw std::runtime_error("Invalid STL mesh: " + filename);
    }
    return mesh;
}

AdvancedMeshLoader::AdvancedMeshLoader(const MeshLoadOptions& options) : options_(options) {}

std::unique_ptr<VolMesh> AdvancedMeshLoader::loadMesh(const std::string& filename)
{
    auto start = std::chrono::steady_clock::now();
    auto mesh = MeshLoader::loadMesh(filename, options_);
    auto end = std::chrono::steady_clock::now();

    statistics_.node_count = mesh ? mesh->getNodeCount() : 0;
    statistics_.element_count = mesh ? mesh->getElementCount() : 0;
    statistics_.load_time_seconds = std::chrono::duration<double>(end - start).count();

    return mesh;
}

std::unique_ptr<VolMesh> AdvancedMeshLoader::loadMeshFromSoftware(const std::string& filename,
                                                                  MeshFormat format)
{
    auto start = std::chrono::steady_clock::now();
    auto mesh = MeshLoader::loadMesh(filename, format, options_);
    auto end = std::chrono::steady_clock::now();

    statistics_.node_count = mesh ? mesh->getNodeCount() : 0;
    statistics_.element_count = mesh ? mesh->getElementCount() : 0;
    statistics_.load_time_seconds = std::chrono::duration<double>(end - start).count();

    return mesh;
}

std::vector<std::unique_ptr<VolMesh>>
AdvancedMeshLoader::loadMeshes(const std::vector<std::string>& filenames)
{
    std::vector<std::unique_ptr<VolMesh>> result;
    result.reserve(filenames.size());
    for (const auto& filename : filenames)
    {
        result.push_back(loadMesh(filename));
    }
    return result;
}

void AdvancedMeshLoader::setOptions(const MeshLoadOptions& options)
{
    options_ = options;
}

const MeshLoadOptions& AdvancedMeshLoader::getOptions() const
{
    return options_;
}

const MeshLoadStatistics& AdvancedMeshLoader::getStatistics() const
{
    return statistics_;
}

std::string AdvancedMeshLoader::getStatisticsString() const
{
    return statistics_.toString();
}

std::unique_ptr<VolMesh>
GeometryMeshGenerator::generateFromDimensions(const GeometryDimensions& dims,
                                              const MeshGenerationOptions& options)
{
    switch (dims.shape)
    {
    case GeometryShape::BOX:
        return generateBoxMesh(dims, options);
    case GeometryShape::CYLINDER:
        return generateCylinderMesh(dims, options);
    case GeometryShape::PLATE:
        return generatePlateMesh(dims, options);
    case GeometryShape::SPHERE:
    default:
        throw std::runtime_error("GeometryShape currently not supported for generation");
    }
}

std::unique_ptr<VolMesh>
GeometryMeshGenerator::generateBoxMesh(const GeometryDimensions& dims,
                                       const MeshGenerationOptions& options)
{
    auto mesh = std::make_unique<VolMesh>();

    const std::size_t nx = std::max<std::size_t>(1, options.nx);
    const std::size_t ny = std::max<std::size_t>(1, options.ny);
    const std::size_t nz = std::max<std::size_t>(1, options.nz);

    const double dx = dims.length / static_cast<double>(nx);
    const double dy = dims.width / static_cast<double>(ny);
    const double dz = dims.height / static_cast<double>(nz);

    auto nodeIndex = [ny, nz](std::size_t i, std::size_t j, std::size_t k)
    { return i * (ny + 1) * (nz + 1) + j * (nz + 1) + k; };

    for (std::size_t i = 0; i <= nx; ++i)
    {
        for (std::size_t j = 0; j <= ny; ++j)
        {
            for (std::size_t k = 0; k <= nz; ++k)
            {
                mesh->addNode(Geometry::Point3D(i * dx, j * dy, k * dz));
            }
        }
    }

    for (std::size_t i = 0; i < nx; ++i)
    {
        for (std::size_t j = 0; j < ny; ++j)
        {
            for (std::size_t k = 0; k < nz; ++k)
            {
                const std::size_t v000 = nodeIndex(i, j, k);
                const std::size_t v100 = nodeIndex(i + 1, j, k);
                const std::size_t v010 = nodeIndex(i, j + 1, k);
                const std::size_t v110 = nodeIndex(i + 1, j + 1, k);
                const std::size_t v001 = nodeIndex(i, j, k + 1);
                const std::size_t v101 = nodeIndex(i + 1, j, k + 1);
                const std::size_t v011 = nodeIndex(i, j + 1, k + 1);
                const std::size_t v111 = nodeIndex(i + 1, j + 1, k + 1);

                mesh->addElement(ElementType::TETRAHEDRON, {v000, v100, v010, v001});
                mesh->addElement(ElementType::TETRAHEDRON, {v100, v110, v010, v111});
                mesh->addElement(ElementType::TETRAHEDRON, {v100, v010, v001, v111});
                mesh->addElement(ElementType::TETRAHEDRON, {v100, v001, v101, v111});
                mesh->addElement(ElementType::TETRAHEDRON, {v010, v001, v011, v111});
            }
        }
    }

    return mesh;
}

std::unique_ptr<VolMesh>
GeometryMeshGenerator::generateCylinderMesh(const GeometryDimensions& dims,
                                            const MeshGenerationOptions& options)
{
    auto mesh = std::make_unique<VolMesh>();

    const std::size_t nr = std::max<std::size_t>(6, options.radial_segments);
    const std::size_t nz = std::max<std::size_t>(1, options.height_segments);
    const double dz = dims.height / static_cast<double>(nz);

    std::vector<std::size_t> center_nodes(nz + 1);
    std::vector<std::vector<std::size_t>> rings(nz + 1, std::vector<std::size_t>(nr));

    for (std::size_t k = 0; k <= nz; ++k)
    {
        const double z = k * dz;
        center_nodes[k] = mesh->addNode(Geometry::Point3D(0.0, 0.0, z));
        for (std::size_t i = 0; i < nr; ++i)
        {
            const double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(nr);
            const double x = dims.radius * std::cos(theta);
            const double y = dims.radius * std::sin(theta);
            rings[k][i] = mesh->addNode(Geometry::Point3D(x, y, z));
        }
    }

    for (std::size_t k = 0; k < nz; ++k)
    {
        for (std::size_t i = 0; i < nr; ++i)
        {
            const std::size_t in = (i + 1) % nr;

            const auto c0 = center_nodes[k];
            const auto c1 = center_nodes[k + 1];
            const auto a0 = rings[k][i];
            const auto b0 = rings[k][in];
            const auto a1 = rings[k + 1][i];
            const auto b1 = rings[k + 1][in];

            mesh->addElement(ElementType::TETRAHEDRON, {c0, a0, b0, c1});
            mesh->addElement(ElementType::TETRAHEDRON, {c1, a1, b1, a0});
            mesh->addElement(ElementType::TETRAHEDRON, {c1, b1, b0, a0});
        }
    }

    return mesh;
}

std::unique_ptr<VolMesh>
GeometryMeshGenerator::generatePlateMesh(const GeometryDimensions& dims,
                                         const MeshGenerationOptions& options)
{
    GeometryDimensions box = dims;
    box.shape = GeometryShape::BOX;
    box.height = std::max(dims.thickness, 1e-6);

    MeshGenerationOptions opts = options;
    opts.nz = std::max<std::size_t>(1, options.nz);

    return generateBoxMesh(box, opts);
}

} // namespace Mesh
} // namespace SCDAT
