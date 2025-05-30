#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "clip.h" // Include the API implementation

namespace py = pybind11;

PYBIND11_MODULE(loop_cgal, m)
{
    
      m.def("clip_surface", &clip_surface,
          py::arg("tm"),
          py::arg("clipper"),
          py::arg("target_edge_length") = 10.0,
            py::arg("remesh_before_clipping") = true,
            py::arg("remesh_after_clipping") = true,
            py::arg("remove_degenerate_faces") = true,
            py::arg("duplicate_vertex_threshold") = 1e-6,
            py::arg("area_threshold") = 1e-6,
          "Clip one surface with another.");
      py::class_<NumpyMesh>(m, "NumpyMesh")
          .def(py::init<>())
          .def_readwrite("vertices", &NumpyMesh::vertices)
          .def_readwrite("triangles", &NumpyMesh::triangles);

}