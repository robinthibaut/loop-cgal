#ifndef CLIP_H
#define CLIP_H
#include <pybind11/numpy.h>
#include "numpymesh.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> TriangleMesh;
NumpyMesh clip_surface(
    NumpyMesh tm,
    NumpyMesh clipper
);
#endif