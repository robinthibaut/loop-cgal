#ifndef MESH_H
#define MESH_H

#include <CGAL/Plane_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Vector_3.h>
#include <CGAL/property_map.h>
#include <numpymesh.h>
#include <pybind11/numpy.h>
#include <set>
#include <utility> // For std::pair
#include <vector>
#include "meshenums.h"
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_K;
typedef CGAL::Surface_mesh<Exact_K::Point_3> Exact_Mesh;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> TriangleMesh;
typedef CGAL::Plane_3<Kernel> Plane;
typedef CGAL::Vector_3<Kernel> Vector;
class TriMesh
{
public:
        // Constructor
        TriMesh(const std::vector<std::vector<int>> &triangles,
                const std::vector<std::pair<double, double>> &vertices);
        TriMesh(const pybind11::array_t<double> &vertices,
                const pybind11::array_t<int> &triangles);

        // Method to cut the mesh with another surface object
        // Returns true if clipping occurred, false if meshes don't intersect or operation failed
        bool cutWithSurface(TriMesh &surface,
                            bool preserve_intersection = false,
                            bool use_exact_kernel = true);

        // Method to remesh the triangle mesh
        void remesh(bool split_long_edges,  double target_edge_length,
                    int number_of_iterations, bool protect_constraints,
                    bool relax_constraints);
        void init();
        void cut_with_implicit_function(const std::vector<double>& property, double value, ImplicitCutMode cutmode = ImplicitCutMode::KEEP_POSITIVE_SIDE);
        // Getters for mesh properties
        void reverseFaceOrientation();
        NumpyMesh save(double area_threshold, double duplicate_vertex_threshold);
        void add_fixed_edges(const pybind11::array_t<int> &pairs);
        pybind11::array_t<int> orig_vertex_map() const;
        const TriangleMesh& get_mesh() const { return _mesh; }
        void set_mesh(const TriangleMesh& mesh) { _mesh = mesh; }
private:
        std::set<TriangleMesh::Edge_index> _fixedEdges;
        TriangleMesh _mesh; // The underlying CGAL surface mesh
        CGAL::Boolean_property_map<std::set<TriangleMesh::Edge_index>>
            _edge_is_constrained_map;
        TriangleMesh::Property_map<TriangleMesh::Vertex_index, int> _orig_vertex_id;
};

#endif // MESH_HANDLER_H
