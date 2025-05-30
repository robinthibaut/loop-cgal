#include "clip.h"
#include "numpymesh.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>

std::set<TriangleMesh::Edge_index> collect_border_edges(const TriangleMesh &tm)
{
    std::set<TriangleMesh::Edge_index> border_edges;
    for (const auto &halfedge : tm.halfedges())
    {
        if (tm.is_border(halfedge))
        {
            border_edges.insert(CGAL::edge(halfedge, tm)); // Convert halfedge to edge
        }
    }
    return border_edges;
}
double calculate_triangle_area(const std::array<double, 3> &v1, const std::array<double, 3> &v2, const std::array<double, 3> &v3)
{
    // Compute vectors
    std::array<double, 3> vec1 = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
    std::array<double, 3> vec2 = {v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]};

    // Compute cross product
    std::array<double, 3> cross_product = {
        vec1[1] * vec2[2] - vec1[2] * vec2[1],
        vec1[2] * vec2[0] - vec1[0] * vec2[2],
        vec1[0] * vec2[1] - vec1[1] * vec2[0]};

    // Compute magnitude of cross product
    double magnitude = std::sqrt(
        cross_product[0] * cross_product[0] +
        cross_product[1] * cross_product[1] +
        cross_product[2] * cross_product[2]);

    // Area is half the magnitude of the cross product
    return 0.5 * magnitude;
}
NumpyMesh clip_surface(NumpyMesh tm, NumpyMesh clipper, double target_edge_length, bool remesh_before_clipping, bool remesh_after_clipping, bool remove_degenerate_faces, double duplicate_vertex_threshold, double area_threshold)
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
    // Parameters for isotropic remeshing
    const unsigned int number_of_iterations = 3; // Number of remeshing iterations
    if (remesh_before_clipping)
    {
        // Split long edges before remeshing
        CGAL::Polygon_mesh_processing::split_long_edges(edges(_tm), target_edge_length, _tm);
        CGAL::Polygon_mesh_processing::split_long_edges(edges(_clipper), target_edge_length, _clipper);
        // Perform isotropic remeshing on _tm
        std::set<TriangleMesh::Edge_index> protected_edges = collect_border_edges(_tm);
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
            faces(_tm), // Range of faces to remesh
            target_edge_length,
            _tm,
            CGAL::parameters::number_of_iterations(number_of_iterations).protect_constraints(true).edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
        // Perform isotropic remeshing on _clipper
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
            faces(_clipper),
            target_edge_length,
            _clipper,
            CGAL::parameters::number_of_iterations(number_of_iterations));
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
        else
        {
            if (remesh_after_clipping)
            {

                CGAL::Polygon_mesh_processing::stitch_borders(_tm);
                CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycles(_tm);
                CGAL::Polygon_mesh_processing::split_long_edges(edges(_tm), target_edge_length, _tm);

                std::set<TriangleMesh::Edge_index> protected_edges = collect_border_edges(_tm);
                // Perform isotropic remeshing with border edge protection

                CGAL::Polygon_mesh_processing::isotropic_remeshing(
                    faces(_tm), // Range of faces to remesh
                    target_edge_length,
                    _tm,
                    CGAL::parameters::number_of_iterations(number_of_iterations).protect_constraints(true).edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
            }
            if (remove_degenerate_faces)
            {
                std::set<TriangleMesh::Edge_index> protected_edges = collect_border_edges(_tm);
                bool beautify_flag = CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces(faces(_tm), _tm, CGAL::parameters::edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
                if (!beautify_flag)
                {
                    std::cout << "removing degenrate faces failed." << std::endl;
                }
            }
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
    std::set<std::pair<TriangleMesh::Vertex_index, TriangleMesh::Vertex_index>> duplicate_vertex_indices;
    for (const auto &v1 : _tm.vertices())
    {
        const auto &point = _tm.point(v1);
        for (const auto &v2 : _tm.vertices())
        {
            if (v1 == v2)
                continue; // Skip self-comparison
            const auto &point2 = _tm.point(v2);
            // Check if the points are close enough to be considered duplicates.
            // use l2 norm to compare points
            if (CGAL::squared_distance(point, point2) < duplicate_vertex_threshold * duplicate_vertex_threshold)
            {

                if (v1 < v2)
                { // Only add the pair once, maintaining a consistent order
                    duplicate_vertex_indices.insert({v1, v2});
                }
            }
        }
    }
    int index = 0;
    for (const auto &vertex : _tm.vertices())

    {
        // Check if the vertex is a duplicate
        bool is_duplicate = false;
        int first_occurrence_index = -1;
        for (const auto &pair : duplicate_vertex_indices)
        {
            if (pair.second == vertex)
            {
                is_duplicate = true;
                first_occurrence_index = vertex_index_map[pair.first];
                break;
            }
        }
        if (is_duplicate)
        {
            // If the vertex is a duplicate, map it to the first occurrence
            vertex_index_map[vertex] = first_occurrence_index;
            first_occurrence_index = -1; // Reset for next iteration
            continue;
        }

        const auto &point = _tm.point(vertex);
        vertices.push_back({point.x(), point.y(), point.z()});
        vertex_index_map[vertex] = index++;
    }

    std::cout << "Found " << duplicate_vertex_indices.size() << " duplicate vertices." << std::endl;

    std::cout << "After removing duplicates, we have " << vertices.size() << " unique vertices. from" << _tm.vertices().size() << std::endl;
    for (const auto &face : _tm.faces())
    {
        std::array<int, 3> triangle;
        int i = 0;
        for (auto halfedge : CGAL::halfedges_around_face(_tm.halfedge(face), _tm))
        {
            triangle[i++] = vertex_index_map[CGAL::target(halfedge, _tm)]; // Assuming `idx()` gives the vertex index
        }
        // calculate triangle area
        double area = calculate_triangle_area(
            vertices[triangle[0]],
            vertices[triangle[1]],
            vertices[triangle[2]]);
        if (area < area_threshold)
        {
            std::cout << "Found degenerate triangle with area " << area << " for triangle: "
                      << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << std::endl;
            continue; // Skip degenerate triangles
        }
        triangles.push_back(triangle);
    }
    std::cout << "Found " << triangles.size() << " triangles." << std::endl;
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
    return result;
}