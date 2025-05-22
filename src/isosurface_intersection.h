#ifndef ISOSURFACE_INTERSECTION_H
#define ISOSURFACE_INTERSECTION_H
#include <marching_cubes.h>

Mesh compute_optimized_intersection(
    const std::vector<std::vector<std::vector<double>>> &scalar_field1,
    const std::vector<std::vector<std::vector<double>>> &scalar_field2,
    double iso_value1,
    double iso_value2,
    double origin_x, double origin_y, double origin_z,
    double cell_size_x, double cell_size_y, double cell_size_z);

Mesh extract_submesh(
    const Mesh &mesh,
    const Point &min_corner,
    const Point &max_corner);


void merge_meshes(
    const Mesh &mesh1,
    const Mesh &mesh2,
    Mesh &result);

#endif // ISOSURFACE_INTERSECTION_H
