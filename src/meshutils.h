#ifndef MESHUTILS_H
#define MESHUTILS_H
#include "mesh.h"

std::set<TriangleMesh::Edge_index> collect_border_edges(const TriangleMesh &tm);
NumpyMesh export_mesh(const TriangleMesh &tm, double area_threshold,
                      double duplicate_vertex_threshold, bool verbose = false);
double calculate_triangle_area(const std::array<double, 3> &v1,
                               const std::array<double, 3> &v2,
                               const std::array<double, 3> &v3);
#endif // MESHUTILS_H
