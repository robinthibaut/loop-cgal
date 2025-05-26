#ifndef CLIP_H
#define CLIP_H
#include <pybind11/numpy.h>
#include "numpymesh.h"

NumpyMesh clip_surface(
    pybind11::array_t<ssize_t> surface_1_triangles,
    pybind11::array_t<double> surface_1_vertices,
    pybind11::array_t<ssize_t> surface_2_triangles,
    pybind11::array_t<double> surface_2_vertices);
#endif