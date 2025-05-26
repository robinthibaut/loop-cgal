#ifndef API_H
#define API_H

#include <pybind11/numpy.h>
#include <utility>
struct NumpyMesh
{
    pybind11::array_t<double> vertices;
    pybind11::array_t<int> triangles;
};
NumpyMesh generate_mesh_from_numpy(
    pybind11::array_t<double> scalar_field,
    pybind11::array_t<double> origin,
    pybind11::array_t<double> step_vector,
    pybind11::array_t<int> num_steps,
    double iso_value);
typedef CGAL::Surface_mesh<Point> TriangleMesh;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
std::vector<NumpyMesh> calculate_mesh_intersection(
                                pybind11::array_t<double> scalar_field1,
                                pybind11::array_t<double> scalar_field2,
                                pybind11::array_t<double> origin,
                                pybind11::array_t<double> step_vector,
                                pybind11::array_t<int> num_steps,
                                double iso_value1,
                                double iso_value2);
NumpyMesh clip_surface(
                                pybind11::array_t<ssize_t> surface_1_triangles,
                                pybind11::array_t<double> surface_1_vertices,
                                pybind11::array_t<ssize_t> surface_2_triangles,
                                pybind11::array_t<double> surface_2_vertices
                                );


#endif // API_H