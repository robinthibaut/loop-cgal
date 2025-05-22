#include "marching_cubes.h"
#include "edge_table.h"
MarchingCubes::MarchingCubes(const ScalarFields &scalar_fields,
                             const SurfaceIsoValues &iso_values,
                                const TruncationRules &truncation_rules,
                             const Point &grid_origin,
                             double grid_spacing)
    : scalar_fields_(scalar_fields), iso_values_(iso_values), grid_origin_(grid_origin), grid_spacing_(grid_spacing), truncation_rules_(truncation_rules)   
{
    num_scalar_fields_ = scalar_fields.size();
    num_vertices_x_ = scalar_fields[0].size();
    num_vertices_y_ = scalar_fields[0][0].size();
    num_vertices_z_ = scalar_fields[0][0][0].size();
}

std::vector<Mesh> MarchingCubes::generate_mesh()
{
    std::vector<Mesh> meshes;
    for (size_t i = 0; i < num_scalar_fields_; ++i)
    {
        meshes.emplace_back();
    }
    for (size_t i = 0; i < num_vertices_x_ - 1; ++i)
    {
        for (size_t j = 0; j < num_vertices_y_ - 1; ++j)
        {
            for (size_t k = 0; k < num_vertices_z_ - 1; ++k)
            {
                process_cube(i, j, k, meshes);
            }
        }
    }
    
    return meshes;
}

int MarchingCubes::compute_cube_index(const std::array<double, 8> &cube_values, double isovalue) const
{
    // Compute cube index based on scalar values
    int index = 0;
    for (int i = 0; i < 8; ++i)
    {
        if (cube_values[i] < isovalue)
        {
            index |= (1 << i);
        }
    }
    return index;
}

Point MarchingCubes::interpolate_vertex(const Point &p1, const Point &p2, double val1, double val2, double isovalue) const
{
    // Linear interpolation of vertex position
    double t = (isovalue - val1) / (val2 - val1);
    return Point(p1.x() + t * (p2.x() - p1.x()),
                 p1.y() + t * (p2.y() - p1.y()),
                 p1.z() + t * (p2.z() - p1.z()));
}

int MarchingCubes::process_cube(int x, int y, int z, std::vector<Mesh> &meshes)
{
    // Define the 8 corners of the cube in grid coordinates
    std::array<Point, 8> cube_corners = {
        grid_origin_ + Vector(x * grid_spacing_, y * grid_spacing_, z * grid_spacing_),
        grid_origin_ + Vector((x + 1) * grid_spacing_, y * grid_spacing_, z * grid_spacing_),
        grid_origin_ + Vector((x + 1) * grid_spacing_, (y + 1) * grid_spacing_, z * grid_spacing_),
        grid_origin_ + Vector(x * grid_spacing_, (y + 1) * grid_spacing_, z * grid_spacing_),
        grid_origin_ + Vector(x * grid_spacing_, y * grid_spacing_, (z + 1) * grid_spacing_),
        grid_origin_ + Vector((x + 1) * grid_spacing_, y * grid_spacing_, (z + 1) * grid_spacing_),
        grid_origin_ + Vector((x + 1) * grid_spacing_, (y + 1) * grid_spacing_, (z + 1) * grid_spacing_),
        grid_origin_ + Vector(x * grid_spacing_, (y + 1) * grid_spacing_, (z + 1) * grid_spacing_)};
    for (int i = 0; i < iso_values_.size(); ++i)
    {
        double isovalue = iso_values_[i].second;
        int scalar_field_index = iso_values_[i].first;
        // Check if the cube intersects with the isosurface
        // Get the scalar values at the 8 corners of the cube
        std::array<double, 8> cube_values = {
            scalar_fields_[scalar_field_index][x][y][z],
            scalar_fields_[scalar_field_index][x + 1][y][z],
            scalar_fields_[scalar_field_index][x + 1][y + 1][z],
            scalar_fields_[scalar_field_index][x][y + 1][z],
            scalar_fields_[scalar_field_index][x][y][z + 1],
            scalar_fields_[scalar_field_index][x + 1][y][z + 1],
            scalar_fields_[scalar_field_index][x + 1][y + 1][z + 1],
            scalar_fields_[scalar_field_index][x][y + 1][z + 1]};

        // Compute the cube index
        int cube_index = compute_cube_index(cube_values, isovalue);

        // If the cube is entirely inside or outside the isosurface, skip it
        if (cube_index == 0 || cube_index == 255)
        {
            continue;
        }
        else {
            const int edges = edgeTable[cube_index]; // edgeTable is a predefined lookup table

            // Interpolate the vertices on the intersected edges
            std::array<Point, 12> edge_vertices;
            for (int i = 0; i < 12; ++i)
            {
                if (edges & (1<<i))
                {
                    int v1 = edgeConnection[i][0]; // edgeConnection maps edges to corner indices
                    int v2 = edgeConnection[i][1];
                    edge_vertices[i] = interpolate_vertex(cube_corners[v1], cube_corners[v2], cube_values[v1], cube_values[v2],isovalue);
                }
            }

            // Generate triangles for the current cube
            const auto &triangles = triTable[cube_index]; // triTable is a predefined lookup table
            for (int i = 0; triangles[i] != -1; i += 3)
            {
                // Add a triangle to the mesh
                auto v1 = meshes[i].add_vertex(edge_vertices[triangles[i]]);
                auto v2 = meshes[i].add_vertex(edge_vertices[triangles[i + 1]]);
                auto v3 = meshes[i].add_vertex(edge_vertices[triangles[i + 2]]);
                meshes[i].add_face(v1, v2, v3);
            }
        }
           
    }


}