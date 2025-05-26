#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "clip.h" // Include the API implementation

namespace py = pybind11;

PYBIND11_MODULE(loop_cgal, m)
{
    
      m.def("clip_surface", &clip_surface,
          py::arg("surface_1_triangles"),
          py::arg("surface_1_vertices"),
          py::arg("surface_2_triangles"),
          py::arg("surface_2_vertices"),
          "Clip one surface with another.");
      py::class_<NumpyMesh>(m, "NumpyMesh")
          .def(py::init<>())
          .def_readwrite("vertices", &NumpyMesh::vertices)
          .def_readwrite("triangles", &NumpyMesh::triangles);

}