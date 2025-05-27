#ifndef NUMPYMESH_H
#define NUMPYMESH_H
#include <pybind11/numpy.h>
struct NumpyMesh
{
    pybind11::array_t<double> vertices;
    pybind11::array_t<int> triangles;
};
#endif // NUMPYMESH_H