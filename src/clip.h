#ifndef CLIP_H
#define CLIP_H
#include <pybind11/numpy.h>
#include "numpymesh.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> TriangleMesh;
std::set<TriangleMesh::Edge_index> collect_border_edges(const TriangleMesh& tm);
double calculate_triangle_area(const std::array<double, 3>& v1, const std::array<double, 3>& v2, const std::array<double, 3>& v3);
NumpyMesh clip_surface(
    NumpyMesh tm,
    NumpyMesh clipper,
    double target_edge_length = 10.0
    , bool remesh_before_clipping = true,
    bool remesh_after_clipping = true,
    bool remove_degenerate_faces = true,
    double duplicate_vertex_threshold = 1e-6,
    double area_threshold = 1e-6
);
#endif