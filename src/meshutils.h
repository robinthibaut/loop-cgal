#ifndef MESHUTILS_H
#define MESHUTILS_H
#include "mesh.h"
#include <set>
#include <string>

std::set<TriangleMesh::Edge_index> collect_border_edges(const TriangleMesh &tm);
NumpyMesh export_mesh(const TriangleMesh &tm, double area_threshold,
                      double duplicate_vertex_threshold);
double calculate_triangle_area(const std::array<double, 3> &v1,
                               const std::array<double, 3> &v2,
                               const std::array<double, 3> &v3);
Exact_Mesh convert_to_exact(const TriMesh& input);
TriangleMesh convert_to_double_mesh(const Exact_Mesh& input);

// Attempt to repair a mesh so that subsequent CGAL operations do not abort.
// Returns true if CGAL considers the mesh a valid polygon mesh after repair.
bool ensure_valid_mesh(TriangleMesh &tm, const std::string &label,
                       bool verbose);

#endif // MESHUTILS_H
