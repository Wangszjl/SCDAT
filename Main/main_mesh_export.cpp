#include "Mesh/include/MeshParsing.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace
{

std::unordered_map<int, std::string> parseSurfacePhysicalNames(const std::string& msh_file)
{
    std::ifstream in(msh_file);
    if (!in)
    {
        return {};
    }

    std::unordered_map<int, std::string> names;
    std::string line;
    while (std::getline(in, line))
    {
        if (line != "$PhysicalNames")
        {
            continue;
        }

        if (!std::getline(in, line))
        {
            break;
        }

        int count = 0;
        {
            std::istringstream iss(line);
            iss >> count;
        }

        for (int i = 0; i < count; ++i)
        {
            if (!std::getline(in, line))
            {
                break;
            }

            std::istringstream iss(line);
            int dim = 0;
            int id = 0;
            iss >> dim >> id;

            const auto first_quote = line.find('"');
            const auto last_quote = line.find_last_of('"');
            if (dim == 2 && first_quote != std::string::npos && last_quote != std::string::npos &&
                last_quote > first_quote)
            {
                names[id] = line.substr(first_quote + 1, last_quote - first_quote - 1);
            }
        }

        break;
    }

    return names;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: SCDAT_mesh_export <input.msh> <output_dir>\n";
        return 1;
    }

    const std::string input_file = argv[1];
    const std::filesystem::path output_dir = argv[2];

    try
    {
        std::filesystem::create_directories(output_dir);

        auto mesh = SCDAT::Mesh::MeshLoader::loadMesh(input_file);
        if (!mesh)
        {
            throw std::runtime_error("Failed to load mesh.");
        }

        const auto physical_names = parseSurfacePhysicalNames(input_file);

        const auto nodes_path = output_dir / "nodes.csv";
        const auto triangles_path = output_dir / "triangles.csv";
        const auto tetra_path = output_dir / "tetrahedra.csv";

        std::ofstream nodes_out(nodes_path);
        std::ofstream tri_out(triangles_path);
        std::ofstream tet_out(tetra_path);

        if (!nodes_out || !tri_out || !tet_out)
        {
            throw std::runtime_error("Failed to create output csv files.");
        }

        nodes_out << "id,x,y,z\n";
        const auto& nodes = mesh->getNodes();
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            const auto& p = nodes[i]->getPosition();
            nodes_out << i << ',' << p.x() << ',' << p.y() << ',' << p.z() << '\n';
        }

        tri_out << "eid,n1,n2,n3,boundary_id,boundary_name\n";
        tet_out << "eid,n1,n2,n3,n4\n";

        const auto& records = mesh->getElementRecords();
        std::size_t tri_count = 0;
        std::size_t tet_count = 0;

        for (std::size_t i = 0; i < records.size(); ++i)
        {
            const auto& rec = records[i];
            if (rec.type == SCDAT::Mesh::ElementType::TRIANGLE && rec.node_indices.size() == 3)
            {
                std::string boundary_name;
                const auto it = physical_names.find(rec.material_id);
                if (it != physical_names.end())
                {
                    boundary_name = it->second;
                }

                tri_out << i << ',' << rec.node_indices[0] << ',' << rec.node_indices[1] << ','
                        << rec.node_indices[2] << ',' << rec.material_id << ',' << boundary_name
                        << '\n';
                ++tri_count;
            }
            else if (rec.type == SCDAT::Mesh::ElementType::TETRAHEDRON &&
                     rec.node_indices.size() == 4)
            {
                tet_out << i << ',' << rec.node_indices[0] << ',' << rec.node_indices[1] << ','
                        << rec.node_indices[2] << ',' << rec.node_indices[3] << '\n';
                ++tet_count;
            }
        }

        std::cout << "Mesh loaded successfully.\n";
        std::cout << "Nodes: " << mesh->getNodeCount() << "\n";
        std::cout << "Elements: " << mesh->getElementCount() << "\n";
        std::cout << "Triangle elements: " << tri_count << "\n";
        std::cout << "Tetrahedron elements: " << tet_count << "\n";
        std::cout << "Exported: " << nodes_path.string() << "\n";
        std::cout << "Exported: " << triangles_path.string() << "\n";
        std::cout << "Exported: " << tetra_path.string() << "\n";

        return 0;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 2;
    }
}
