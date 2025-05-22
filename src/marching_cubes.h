#pragma once

#include <vector>
#include <unordered_map>
#include <array>
#include <unordered_map>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

namespace std
{
    template <>
    struct hash<std::array<int, 6>>
    {
        std::size_t operator()(const std::array<int, 6> &arr) const
        {
            std::size_t seed = 0;
            for (int val : arr)
            {
                seed ^= std::hash<int>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
struct GridCell
{
    int i, j, k; // Cell indices

    bool operator==(const GridCell &other) const
    {
        return i == other.i && j == other.j && k == other.k;
    }
};

// Hash function for GridCell
struct GridCellHash
{
    std::size_t operator()(const GridCell &cell) const
    {
        return std::hash<int>()(cell.i) ^
               std::hash<int>()(cell.j) ^
               std::hash<int>()(cell.k);
    }
};
typedef std::unordered_set<GridCell, GridCellHash> active_cells_set;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef std::unordered_map<int, std::vector<std::pair<int, double>>> TruncationRules;
typedef std::vector<std::pair<int, double>> SurfaceIsoValues;
typedef std::vector<std::vector<std::vector<std::vector<double>>>> ScalarFields;

class MarchingCubes
{
public:
    MarchingCubes(const ScalarFields &scalar_fields,
                  const SurfaceIsoValues &iso_values,
                  const TruncationRules &truncation_rules,
                  const Point &grid_origin,
                  double grid_spacing);

    std::vector<Mesh> generate_mesh();

private:
    const ScalarFields &scalar_fields_;
    int num_scalar_fields_;
    int num_vertices_x_;
    int num_vertices_y_;
    int num_vertices_z_;
    const SurfaceIsoValues &iso_values_;
    const TruncationRules &truncation_rules_;
    Point grid_origin_;
    double grid_spacing_;

    std::unordered_map<std::array<int, 6>, typename Mesh::Vertex_index> edge_vertex_map_;

    int compute_cube_index(const std::array<double, 8> &cube_values, double isovalue) const;
    Point interpolate_vertex(const Point &p1, const Point &p2, double val1, double val2, double isovalue) const;
    int process_cube(int x, int y, int z, std::vector<Mesh> &meshes);
};