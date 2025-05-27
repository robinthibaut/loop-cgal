#include "clip.h"
#include "numpymesh.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>

NumpyMesh clip_surface(NumpyMesh tm, NumpyMesh clipper)
{
    auto vertices1_buf = tm.vertices.unchecked<2>();
    auto triangles1_buf = tm.triangles.unchecked<2>();
    auto vertices2_buf = clipper.vertices.unchecked<2>();
    auto triangles2_buf = clipper.triangles.unchecked<2>();

    TriangleMesh _tm;
    TriangleMesh _clipper;
    std::vector<TriangleMesh::Vertex_index> vertex_indices;
    // assemble cgal mesh objects from numpy/pybind11 arrays
    for (ssize_t i = 0; i < vertices1_buf.shape(0); ++i)
    {
        vertex_indices.push_back(_tm.add_vertex(Point(vertices1_buf(i, 0), vertices1_buf(i, 1), vertices1_buf(i, 2))));
    }
    for (ssize_t i = 0; i < triangles1_buf.shape(0); ++i)
    {
        _tm.add_face(vertex_indices[triangles1_buf(i, 0)], vertex_indices[triangles1_buf(i, 1)], vertex_indices[triangles1_buf(i, 2)]);
    }

    vertex_indices.clear();
    for (ssize_t i = 0; i < vertices2_buf.shape(0); ++i)
    {
        vertex_indices.push_back(_clipper.add_vertex(Point(vertices2_buf(i, 0), vertices2_buf(i, 1), vertices2_buf(i, 2))));
    }
    for (ssize_t i = 0; i < triangles2_buf.shape(0); ++i)
    {
        _clipper.add_face(vertex_indices[triangles2_buf(i, 0)], vertex_indices[triangles2_buf(i, 1)], vertex_indices[triangles2_buf(i, 2)]);
    }
    for (const auto &vertex : _tm.vertices())
    {
        const auto &point = _tm.point(vertex);
    }
    for (const auto &vertex : _clipper.vertices())
    {
        const auto &point = _clipper.point(vertex);
    }
    if (!CGAL::is_valid_polygon_mesh(_tm))
    {
        std::cerr << "tm is invalid!" << std::endl;
    }
    if (!CGAL::is_valid_polygon_mesh(_clipper))
    {
        std::cerr << "clipper is invalid!" << std::endl;
    }
    // make sure the meshes actually intersect. If they don't, just return mesh 1
    bool intersection = CGAL::Polygon_mesh_processing::do_intersect(_tm, _clipper);
    if (intersection)
    {
        // Clip tm with clipper
        bool flag = CGAL::Polygon_mesh_processing::clip(_tm, _clipper);
        if (!flag)
        {
            std::cerr << "Clipping failed." << std::endl;
            return {};
        }
    }
    else
    {
        std::cout << "Meshes do not intersect. Returning tm." << std::endl;
    }
    // store the result in a numpymesh object for sending back to Python
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::map<TriangleMesh::Vertex_index, int> vertex_index_map;
    int index = 0;
    for (const auto &vertex : _tm.vertices())

    {
        const auto &point = _tm.point(vertex);
        vertices.push_back({point.x(), point.y(), point.z()});
        vertex_index_map[vertex] = index++;
    }
    for (const auto &face : _tm.faces())
    {
        std::array<int, 3> triangle;
        int i = 0;
        for (auto halfedge : CGAL::halfedges_around_face(_tm.halfedge(face), _tm))
        {
            triangle[i++] = vertex_index_map[CGAL::target(halfedge, _tm)]; // Assuming `idx()` gives the vertex index
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
    for (const auto &vertex : _tm.vertices())
    {
        const auto &point = _tm.point(vertex);
    }
    // Return the NumpyMesh object
    return result;
}