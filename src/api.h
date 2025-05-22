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

std::vector<NumpyMesh> calculate_mesh_intersection(
                                pybind11::array_t<double> scalar_field1,
                                pybind11::array_t<double> scalar_field2,
                                pybind11::array_t<double> origin,
                                pybind11::array_t<double> step_vector,
                                pybind11::array_t<int> num_steps,
                                double iso_value1,
                                double iso_value2);

#endif // API_H