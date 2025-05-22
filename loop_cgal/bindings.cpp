#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "marching_cubes.h"
#include "isosurface_intersection.h"
#include "api.h" // Include the API implementation

namespace py = pybind11;

PYBIND11_MODULE(loop_cgal, m)
{
    // Bind the new API function
    m.def("generate_mesh_from_numpy", &generate_mesh_from_numpy,
          py::arg("scalar_field"),
          py::arg("origin"),
          py::arg("step_vector"),
          py::arg("num_steps"),
          py::arg("iso_value"),
          "Generate a mesh from a scalar field using Marching Cubes.");
    m.def("compute_optimized_intersection", &compute_optimized_intersection,
            py::arg("scalar_field1"),
            py::arg("scalar_field2"),
            py::arg("iso_value1"),
            py::arg("iso_value2"),
            py::arg("origin_x"),
            py::arg("origin_y"),
            py::arg("origin_z"),
            py::arg("cell_size_x"),
            py::arg("cell_size_y"),
            py::arg("cell_size_z"),
            "Compute the optimized intersection of two isosurfaces.");
            


}