#include "clip.h"
#include "numpymesh.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/version.h>


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
TriangleMesh load_mesh(NumpyMesh mesh, bool verbose)
{
    TriangleMesh tm;
    auto vertices_buf = mesh.vertices.unchecked<2>();
    auto triangles_buf = mesh.triangles.unchecked<2>();
    std::vector<TriangleMesh::Vertex_index> vertex_indices;

    if (verbose)
    {
        std::cout << "Loading mesh with " << vertices_buf.shape(0) << " vertices and " << triangles_buf.shape(0) << " triangles." << std::endl;
    }

    // Assemble CGAL mesh objects from numpy/pybind11 arrays
    for (ssize_t i = 0; i < vertices_buf.shape(0); ++i)
    {
        vertex_indices.push_back(tm.add_vertex(Point(vertices_buf(i, 0), vertices_buf(i, 1), vertices_buf(i, 2))));
    }
    for (ssize_t i = 0; i < triangles_buf.shape(0); ++i)
    {
        tm.add_face(vertex_indices[triangles_buf(i, 0)], vertex_indices[triangles_buf(i, 1)], vertex_indices[triangles_buf(i, 2)]);
    }

    if (verbose)
    {
        std::cout << "Loaded mesh with " << tm.number_of_vertices() << " vertices and " << tm.number_of_faces() << " faces." << std::endl;
    }

    return tm;
}

Plane load_plane(NumpyPlane plane, bool verbose)
{
    auto normal_buf = plane.normal.unchecked<1>();
    auto point_buf = plane.origin.unchecked<1>();
    if (verbose)
    {
        std::cout << "Loading plane with normal (" << normal_buf(0) << ", " << normal_buf(1) << ", " << normal_buf(2) << ") and point (" << point_buf(0) << ", " << point_buf(1) << ", " << point_buf(2) << ")." << std::endl;
    }
    return Plane(Point(point_buf(0), point_buf(1), point_buf(2)), Vector(normal_buf(0), normal_buf(1), normal_buf(2)));
}
void refine_mesh(TriangleMesh &mesh, bool split_long_edges , bool verbose , double target_edge_length , int number_of_iterations, bool protect_constraints, bool relax_constraints)
{
    // Split long edges before remeshing
    CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length, mesh);
    // Perform isotropic remeshing on _tm
    std::set<TriangleMesh::Edge_index> protected_edges = collect_border_edges(mesh);
    CGAL::Polygon_mesh_processing::isotropic_remeshing(
        faces(mesh), // Range of faces to remesh
        target_edge_length,
        mesh,
        CGAL::parameters::number_of_iterations(number_of_iterations).protect_constraints(protect_constraints).relax_constraints(relax_constraints).edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
    // Perform isotropic remeshing on _clipper
    if (verbose)
    {
        std::cout << "Refined mesh with " << mesh.number_of_vertices() << " vertices and " << mesh.number_of_faces() << " faces." << std::endl;
    }
    return;
}
NumpyMesh export_mesh(const TriangleMesh &tm, double area_threshold, double duplicate_vertex_threshold, bool verbose)
{
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::map<TriangleMesh::Vertex_index, int> vertex_index_map;
    std::set<std::pair<TriangleMesh::Vertex_index, TriangleMesh::Vertex_index>> duplicate_vertex_indices;
    for (const auto &v1 : tm.vertices())
    {
        const auto &point = tm.point(v1);
        for (const auto &v2 : tm.vertices())
        {
            if (v1 == v2)
                continue; // Skip self-comparison
            const auto &point2 = tm.point(v2);
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
    for (const auto &vertex : tm.vertices())

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

        const auto &point = tm.point(vertex);
        vertices.push_back({point.x(), point.y(), point.z()});
        vertex_index_map[vertex] = index++;
    }
    if (verbose)
    {
        std::cout << "Vertices after remeshing: " << vertices.size() << std::endl;

        std::cout << "Found " << duplicate_vertex_indices.size() << " duplicate vertices." << std::endl;

        std::cout << "After removing duplicates, we have " << vertices.size() << " unique vertices from " << tm.vertices().size() << std::endl;
    }
    for (const auto &face : tm.faces())
    {
        std::array<int, 3> triangle;
        int i = 0;
        for (auto halfedge : CGAL::halfedges_around_face(tm.halfedge(face), tm))
        {
            triangle[i++] = vertex_index_map[CGAL::target(halfedge, tm)]; // Assuming `idx()` gives the vertex index
        }
        // calculate triangle area
        double area = calculate_triangle_area(
            vertices[triangle[0]],
            vertices[triangle[1]],
            vertices[triangle[2]]);
        if (area < area_threshold)
        {
            if (verbose)
            {
                std::cout << "Found degenerate triangle with area " << area << " for triangle: "
                          << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << std::endl;
            }
            continue; // Skip degenerate triangles
        }
        triangles.push_back(triangle);
    }
    if (verbose)
    {
        std::cout << "Found " << triangles.size() << " triangles." << std::endl;
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
    return result;
}

NumpyMesh clip_plane(NumpyMesh tm, NumpyPlane clipper, double target_edge_length, bool remesh_before_clipping, bool remesh_after_clipping, bool remove_degenerate_faces, double duplicate_vertex_threshold, double area_threshold, bool verbose)
{
    int number_of_iterations = 3; // Number of remeshing iterations
    if (verbose)
    {
        std::cout << "Starting clipping process." << std::endl;
        std::cout << "Loading data from NumpyMesh." << std::endl;
    }
    TriangleMesh _tm = load_mesh(tm, verbose);
    if (verbose)
    {
        std::cout << "Loaded mesh." << std::endl;
    }
    Plane _clipper = load_plane(clipper, verbose);
    if (verbose)
    {
        std::cout << "Loaded plane." << std::endl;
    }
    if (remesh_before_clipping)
    {
        if (verbose)
        {
            std::cout << "Remeshing before clipping." << std::endl;
        }
        refine_mesh(_tm, true, verbose, target_edge_length, number_of_iterations,true,false);
        // refine_mesh(_clipper, true, verbose, target_edge_length, number_of_iterations);

        if (verbose)
        {
            std::cout << "Remeshing before clipping done." << std::endl;
        }
    }

    // make sure the meshes actually intersect. If they don't, just return mesh 1
    // bool intersection = CGAL::Polygon_mesh_processing::do_intersect(_tm, _clipper);
    bool intersection = true;
    if (intersection)
    {
        // Clip tm with clipper
        if (verbose)
        {
            std::cout << "Clipping tm with clipper." << std::endl;
        }
        bool flag = CGAL::Polygon_mesh_processing::clip(_tm, _clipper);
        if (verbose)
        {
            std::cout << "Clipping done." << std::endl;
        }
        if (!flag)
        {
            std::cerr << "Clipping failed." << std::endl;
            return {};
        }
        else
        {
            if (remesh_after_clipping)
            {

                if (verbose)
                {
                    std::cout << "Remeshing after clipping." << std::endl;
                }
                CGAL::Polygon_mesh_processing::stitch_borders(_tm);
                CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycles(_tm);
                refine_mesh(_tm, true, verbose, target_edge_length, number_of_iterations);

                if (verbose)
                {
                    std::cout << "Remeshing after clipping done." << std::endl;
                }
            }
            if (remove_degenerate_faces)
            {
                if (verbose)
                {
                    std::cout << "Removing degenerate faces." << std::endl;
                }
                std::set<TriangleMesh::Edge_index> protected_edges = collect_border_edges(_tm);
                #if CGAL_VERSION_NR >= 1060000000
                bool beautify_flag = CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces(faces(_tm), _tm, CGAL::parameters::edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
#else
                bool beautify_flag = CGAL::Polygon_mesh_processing::remove_degenerate_faces(faces(_tm), _tm, CGAL::parameters::edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
#endif                if (!beautify_flag)
                {
                    std::cout << "removing degenrate faces failed." << std::endl;
                }
                if (verbose)
                {
                    std::cout << "Removing degenerate faces done." << std::endl;
                }
            }
        }
    }
    else
    {
        std::cout << "Meshes do not intersect. Returning tm." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Clipping done." << std::endl;
    }

    // store the result in a numpymesh object for sending back to Python

    NumpyMesh result = export_mesh(_tm, area_threshold, duplicate_vertex_threshold, verbose);
    if (verbose)
    {
        std::cout << "Exported clipped mesh with " << result.vertices.shape(0) << " vertices and " << result.triangles.shape(0) << " triangles." << std::endl;
    }
    return result;
}
NumpyMesh clip_surface(NumpyMesh tm, NumpyMesh clipper, double target_edge_length, bool remesh_before_clipping, bool remesh_after_clipping, bool remove_degenerate_faces, double duplicate_vertex_threshold, double area_threshold, bool verbose)
{
    if (verbose)
    {
        std::cout << "Starting clipping process." << std::endl;
        std::cout << "Loading data from NumpyMesh." << std::endl;
    }
    TriangleMesh _tm = load_mesh(tm, verbose);
    TriangleMesh _clipper = load_mesh(clipper, verbose);
    if (verbose)
    {
        std::cout << "Loaded meshes." << std::endl;
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
        if (verbose)
        {
            std::cout << "Remeshing before clipping." << std::endl;
        }
        refine_mesh(_tm, true, verbose, target_edge_length, number_of_iterations,true,true);
        refine_mesh(_clipper, true, verbose, target_edge_length, number_of_iterations,true,true);

        if (verbose)
        {
            std::cout << "Remeshing before clipping done." << std::endl;
        }
    }

    // make sure the meshes actually intersect. If they don't, just return mesh 1
    bool intersection = CGAL::Polygon_mesh_processing::do_intersect(_tm, _clipper);
    if (intersection)
    {
        // Clip tm with clipper
        if (verbose)
        {
            std::cout << "Clipping tm with clipper." << std::endl;
        }
        bool flag = CGAL::Polygon_mesh_processing::clip(_tm, _clipper);
        if (verbose)
        {
            std::cout << "Clipping done." << std::endl;
        }
        if (!flag)
        {
            std::cerr << "Clipping failed." << std::endl;
            return {};
        }
        else
        {
            if (remesh_after_clipping)
            {

                if (verbose)
                {
                    std::cout << "Remeshing after clipping." << std::endl;
                }
                CGAL::Polygon_mesh_processing::stitch_borders(_tm);
                CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycles(_tm);
                refine_mesh(_tm, true, verbose, target_edge_length, number_of_iterations);

                if (verbose)
                {
                    std::cout << "Remeshing after clipping done." << std::endl;
                }
            }
            if (remove_degenerate_faces)
            {
                if (verbose)
                {
                    std::cout << "Removing degenerate faces." << std::endl;
                }
                std::set<TriangleMesh::Edge_index> protected_edges = collect_border_edges(_tm);

#if CGAL_VERSION_NR >= 1060000000
                bool beautify_flag = CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces(faces(_tm), _tm, CGAL::parameters::edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
#else
                bool beautify_flag = CGAL::Polygon_mesh_processing::remove_degenerate_faces(faces(_tm), _tm, CGAL::parameters::edge_is_constrained_map(CGAL::make_boolean_property_map(protected_edges)));
#endif
                if (!beautify_flag)
                {
                    std::cout << "removing degenrate faces failed." << std::endl;
                }
                if (verbose)
                {
                    std::cout << "Removing degenerate faces done." << std::endl;
                }
            }
        }
    }
    else
    {
        std::cout << "Meshes do not intersect. Returning tm." << std::endl;
    }
    if (verbose)
    {
        std::cout << "Clipping done." << std::endl;
    }

    // store the result in a numpymesh object for sending back to Python

    NumpyMesh result = export_mesh(_tm, area_threshold, duplicate_vertex_threshold, verbose);
    if (verbose)
    {
        std::cout << "Exported clipped mesh with " << result.vertices.shape(0) << " vertices and " << result.triangles.shape(0) << " triangles." << std::endl;
    }
    return result;
}

std::vector<NumpyMesh> corefine_mesh(NumpyMesh tm1, NumpyMesh tm2,double target_edge_length, double duplicate_vertex_threshold, double area_threshold, int number_of_iterations, bool relax_constraints, bool protect_constraints,bool verbose)
{
    // Load the meshes
    TriangleMesh _tm1 = load_mesh(tm1, false);
    TriangleMesh _tm2 = load_mesh(tm2, false);
    CGAL::Polygon_mesh_processing::split_long_edges(edges(_tm1), target_edge_length, _tm1);
    CGAL::Polygon_mesh_processing::split_long_edges(edges(_tm2), target_edge_length, _tm2);

    // Perform corefinement
    CGAL::Polygon_mesh_processing::corefine(_tm1, _tm2);
    // Find shared edges
    std::set<TriangleMesh::Edge_index> tm_1_shared_edges;
    std::set<TriangleMesh::Edge_index> tm_2_shared_edges;
    for (const auto &edge1 : _tm1.edges())
    {
        Point p1 = _tm1.point(CGAL::source(edge1, _tm1));
        Point p2 = _tm1.point(CGAL::target(edge1, _tm1));

        for (const auto &edge2 : _tm2.edges())
        {
            Point q1 = _tm2.point(CGAL::source(edge2, _tm2));
            Point q2 = _tm2.point(CGAL::target(edge2, _tm2));

            // Check if the edges are identical (considering both orientations)
            if ((p1 == q1 && p2 == q2) || (p1 == q2 && p2 == q1))
            {
                tm_1_shared_edges.insert(edge1);
                tm_2_shared_edges.insert(edge2);
                break;
            }
        }
    }
    std::cout << "Found " << tm_1_shared_edges.size() << " shared edges in tm1 and " << tm_2_shared_edges.size() << " shared edges in tm2." << std::endl;
    
    // std::set<TriangleMesh::Edge_index> constrained_edges;
    
    std::set<TriangleMesh::Edge_index> boundary_edges = collect_border_edges(_tm1);
    std::set<TriangleMesh::Edge_index> boundary_edges2 = collect_border_edges(_tm2);

    tm_1_shared_edges.insert(boundary_edges.begin(), boundary_edges.end());
    tm_2_shared_edges.insert(boundary_edges2.begin(), boundary_edges2.end());
    // Refine the meshes
    // Perform isotropic remeshing on _tm
    CGAL::Polygon_mesh_processing::isotropic_remeshing(
        faces(_tm1), // Range of faces to remesh
        target_edge_length,
        _tm1,
        CGAL::parameters::number_of_iterations(number_of_iterations).edge_is_constrained_map(CGAL::make_boolean_property_map(tm_1_shared_edges)).relax_constraints(relax_constraints).protect_constraints(protect_constraints));
    CGAL::Polygon_mesh_processing::isotropic_remeshing(
        faces(_tm2), // Range of faces to remesh
        target_edge_length,
        _tm2,
        CGAL::parameters::number_of_iterations(number_of_iterations).edge_is_constrained_map(CGAL::make_boolean_property_map(tm_2_shared_edges)).relax_constraints(relax_constraints).protect_constraints(protect_constraints));

    std::cout << "Corefinement done." << std::endl;

    return {export_mesh(_tm1, area_threshold, duplicate_vertex_threshold, verbose), export_mesh(_tm2, area_threshold, duplicate_vertex_threshold, verbose)};
}

