#include <pybind11/numpy.h>
#include "marching_cubes.h"
#include "api.h"
NumpyMesh generate_mesh_from_numpy(
    pybind11::array_t<double> scalar_field,
    pybind11::array_t<double> origin,
    pybind11::array_t<double> step_vector,
    pybind11::array_t<int> num_steps,
    double iso_value)
{
    // Convert NumPy arrays to C++ types
    auto scalar_field_buf = scalar_field.unchecked<3>(); // Assume 3D scalar field
    auto origin_buf = origin.unchecked<1>();
    auto step_vector_buf = step_vector.unchecked<1>();
    auto num_steps_buf = num_steps.unchecked<1>();

    // Extract origin, step vector, and number of steps
    Point grid_origin(origin_buf(0), origin_buf(1), origin_buf(2));
    double grid_spacing = step_vector_buf(0); // Assume uniform spacing for simplicity

    // Convert scalar field to std::vector
    std::vector<std::vector<std::vector<double>>> scalar_field_vec(
        num_steps_buf(0),
        std::vector<std::vector<double>>(num_steps_buf(1),
                                         std::vector<double>(num_steps_buf(2))));

    for (int x = 0; x < num_steps_buf(0); ++x)
    {
        for (int y = 0; y < num_steps_buf(1); ++y)
        {
            for (int z = 0; z < num_steps_buf(2); ++z)
            {
                scalar_field_vec[x][y][z] = scalar_field_buf(x, y, z);
            }
        }
    }
    // Create MarchingCubes instance and generate the mesh
    // MarchingCubes mc(scalar_field_vec, iso_value, grid_origin, grid_spacing);
    // Mesh mesh;
    // active_cells_set active_cells;
    // std::tie(mesh, active_cells) = mc.generate_mesh();
    

    // // Extract vertices and triangles from the mesh
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;

    // for (const auto &vertex : mesh.vertices())
    // {
    //     const auto &point = mesh.point(vertex);
    //     vertices.push_back({point.x(), point.y(), point.z()});
    // }

    // for (const auto &face : mesh.faces())
    // {
    //     std::array<int, 3> triangle;
    //     int i = 0;
    //     for (auto halfedge : halfedges_around_face(mesh.halfedge(face), mesh))
    //     {
    //         triangle[i++] = target(halfedge, mesh).idx(); // Assuming `idx()` gives the vertex index
    //     }
    //     triangles.push_back(triangle);
    // }

    // Convert vertices to a NumPy array
    pybind11::array_t<double> vertices_array({static_cast<int>(vertices.size()), 3});
    // auto vertices_buf = vertices_array.mutable_unchecked<2>();
    // for (size_t i = 0; i < vertices.size(); ++i)
    // {
    //     vertices_buf(i, 0) = vertices[i][0];
    //     vertices_buf(i, 1) = vertices[i][1];
    //     vertices_buf(i, 2) = vertices[i][2];
    // }

    // Convert triangles to a NumPy array
    pybind11::array_t<int> triangles_array({static_cast<int>(triangles.size()), 3});
    // auto triangles_buf = triangles_array.mutable_unchecked<2>();
    // for (size_t i = 0; i < triangles.size(); ++i)
    // {
        // triangles_buf(i, 0) = triangles[i][0];
        // triangles_buf(i, 1) = triangles[i][1];
        // triangles_buf(i, 2) = triangles[i][2];
    // }

    return {vertices_array, triangles_array};
}

