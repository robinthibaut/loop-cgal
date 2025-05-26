#include "clip.h"
#include "numpymesh.h"
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/does_intersect.h>

NumpyMesh clip_surface(pybind11::array_t<ssize_t> surface_1_triangles,
                       pybind11::array_t<double> surface_1_vertices,
                       pybind11::array_t<ssize_t> surface_2_triangles,
                       pybind11::array_t<double> surface_2_vertices)
{
    auto vertices1_buf = surface_1_vertices.unchecked<2>();
    auto triangles1_buf = surface_1_triangles.unchecked<2>();
    auto vertices2_buf = surface_2_vertices.unchecked<2>();
    auto triangles2_buf = surface_2_triangles.unchecked<2>();

    TriangleMesh mesh1;
    TriangleMesh mesh2;
    std::vector<TriangleMesh::Vertex_index> vertex_indices;
    for (ssize_t i = 0; i < vertices1_buf.shape(0); ++i)
    {
        vertex_indices.push_back(mesh1.add_vertex(Point(vertices1_buf(i, 0), vertices1_buf(i, 1), vertices1_buf(i, 2))));
    }
    for (ssize_t i = 0; i < triangles1_buf.shape(0); ++i)
    {
        mesh1.add_face(vertex_indices[triangles1_buf(i, 0)], vertex_indices[triangles1_buf(i, 1)], vertex_indices[triangles1_buf(i, 2)]);
    }

    vertex_indices.clear();
    for (ssize_t i = 0; i < vertices2_buf.shape(0); ++i)
    {
        vertex_indices.push_back(mesh2.add_vertex(Point(vertices2_buf(i, 0), vertices2_buf(i, 1), vertices2_buf(i, 2))));
    }
    for (ssize_t i = 0; i < triangles2_buf.shape(0); ++i)
    {
        mesh2.add_face(vertex_indices[triangles2_buf(i, 0)], vertex_indices[triangles2_buf(i, 1)], vertex_indices[triangles2_buf(i, 2)]);
    }
    for (const auto &vertex : mesh1.vertices())
    {
        const auto &point = mesh1.point(vertex);
    }
    for (const auto &vertex : mesh2.vertices())
    {
        const auto &point = mesh2.point(vertex);
    }
    if (!CGAL::is_valid_polygon_mesh(mesh1))
    {
        std::cerr << "Mesh1 is invalid!" << std::endl;
    }
    if (!CGAL::is_valid_polygon_mesh(mesh2))
    {
        std::cerr << "Mesh2 is invalid!" << std::endl;
    }
    bool intersection = CGAL::Polygon_mesh_processing::do_intersect(mesh1, mesh2);
    if (intersection)
    {
        bool flag = CGAL::Polygon_mesh_processing::clip(mesh1, mesh2, CGAL::parameters::do_not_modify(true));
        std::cerr << "Meshes intersect." << std::endl;
        if (!flag)
        {
            std::cerr << "Clipping failed." << std::endl;
            return {};
        }
    }

    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::map<TriangleMesh::Vertex_index, int> vertex_index_map;
    int index = 0;
    for (const auto &vertex : mesh1.vertices())

    {
        const auto &point = mesh1.point(vertex);
        vertices.push_back({point.x(), point.y(), point.z()});
        vertex_index_map[vertex] = index++;
    }
    for (const auto &face : mesh1.faces())
    {
        std::array<int, 3> triangle;
        int i = 0;
        for (auto halfedge : CGAL::halfedges_around_face(mesh1.halfedge(face), mesh1))
        {
            triangle[i++] = vertex_index_map[CGAL::target(halfedge, mesh1)]; // Assuming `idx()` gives the vertex index
        }
        triangles.push_back(triangle);
    }
    // Convert vertices to a NumPy array
    pybind11::array_t<double> vertices_array({static_cast<int>(vertices.size()), 3});
    auto vertices_buf = vertices_array.mutable_unchecked<2>();
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices_buf(i, 0) = vertices[i][0];
        vertices_buf(i, 1) = vertices[i][1];
        vertices_buf(i, 2) = vertices[i][2];
    }
    // Convert triangles to a NumPy array
    pybind11::array_t<int> triangles_array({static_cast<int>(triangles.size()), 3});
    auto triangles_buf = triangles_array.mutable_unchecked<2>();
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        triangles_buf(i, 0) = triangles[i][0];
        triangles_buf(i, 1) = triangles[i][1];
        triangles_buf(i, 2) = triangles[i][2];
    }
    // Create NumpyMesh object
    NumpyMesh result;
    result.vertices = vertices_array;
    result.triangles = triangles_array;
    // Return the NumpyMesh object
    for (const auto &vertex : mesh1.vertices())
    {
        const auto &point = mesh1.point(vertex);
    }
    // Return the NumpyMesh object
    return result;
}