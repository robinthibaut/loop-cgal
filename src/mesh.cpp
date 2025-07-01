#include "mesh.h"
#include "meshutils.h"
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/version.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace PMP = CGAL::Polygon_mesh_processing;

TriMesh::TriMesh(const std::vector<std::vector<int>> &triangles,
                 const std::vector<std::pair<double, double>> &vertices) {
  std::vector<TriangleMesh::Vertex_index> vertex_indices;
  bool verbose = true;
  if (verbose) {
    std::cout << "Loading mesh with " << vertices.size() << " vertices and "
              << triangles.size() << " triangles." << std::endl;
  }

  // Assemble CGAL mesh objects from numpy/pybind11 arrays
  for (ssize_t i = 0; i < vertices.size(); ++i) {
    vertex_indices.push_back(
        _mesh.add_vertex(Point(vertices[i].first, vertices[i].second, 0.0)));
  }
  for (ssize_t i = 0; i < triangles.size(); ++i) {
    _mesh.add_face(vertex_indices[triangles[i][0]],
                   vertex_indices[triangles[i][1]],
                   vertex_indices[triangles[i][2]]);
  }

  if (verbose) {
    std::cout << "Loaded mesh with " << _mesh.number_of_vertices()
              << " vertices and " << _mesh.number_of_faces() << " faces."
              << std::endl;
  }
  init();
}

TriMesh::TriMesh(const pybind11::array_t<double> &vertices,
                 const pybind11::array_t<int> &triangles) {
  auto verts = vertices.unchecked<2>();
  auto tris = triangles.unchecked<2>();

  std::vector<TriangleMesh::Vertex_index> vertex_indices;

  for (ssize_t i = 0; i < verts.shape(0); ++i) {
    vertex_indices.push_back(
        _mesh.add_vertex(Point(verts(i, 0), verts(i, 1), verts(i, 2))));
  }

  for (ssize_t i = 0; i < tris.shape(0); ++i) {
    _mesh.add_face(vertex_indices[tris(i, 0)], vertex_indices[tris(i, 1)],
                   vertex_indices[tris(i, 2)]);
  }

  std::cout << "Loaded mesh with " << _mesh.number_of_vertices()
            << " vertices and " << _mesh.number_of_faces() << " faces."
            << std::endl;
  init();
}
void TriMesh::init() {
  _fixedEdges = collect_border_edges(_mesh);
  _edge_is_constrained_map = CGAL::make_boolean_property_map(_fixedEdges);
}
void TriMesh::remesh(bool split_long_edges, bool verbose,
                     double target_edge_length, int number_of_iterations,
                     bool protect_constraints, bool relax_constraints)

{

  // ------------------------------------------------------------------
  // 0.  Guard‑rail: sensible target length w.r.t. bbox
  // ------------------------------------------------------------------
  CGAL::Bbox_3 bb = PMP::bbox(_mesh);
  const double bbox_diag = std::sqrt(CGAL::square(bb.xmax() - bb.xmin()) +
                                     CGAL::square(bb.ymax() - bb.ymin()) +
                                     CGAL::square(bb.zmax() - bb.zmin()));
  PMP::remove_isolated_vertices(_mesh);
  if (target_edge_length < 1e-4 * bbox_diag) {
    if (verbose)
      std::cout << "  ! target_edge_length (" << target_edge_length
                << ") too small – skipping remesh\n";
    return;
  }

  // ------------------------------------------------------------------
  // 1.  Quick diagnostics
  // ------------------------------------------------------------------
  double min_e = std::numeric_limits<double>::max(), max_e = 0.0;
  for (auto e : _mesh.edges()) {
    const double l = PMP::edge_length(e, _mesh);
    min_e = std::min(min_e, l);
    max_e = std::max(max_e, l);
  }
  if (verbose)
    std::cout << "      edge length range: [" << min_e << ", " << max_e
              << "]  target = " << target_edge_length << '\n';

  if (!CGAL::is_valid_polygon_mesh(_mesh, verbose) && verbose)
    std::cout << "      ! mesh is not a valid polygon mesh\n";

  // ------------------------------------------------------------------
  // 2.  Normal isotropic remeshing loop
  // ------------------------------------------------------------------
  // Convert _fixedEdges to a compatible property map

  // Update remeshing calls to use the property map
  if (split_long_edges)
    if (verbose)
      std::cout << "Splitting long edges before remeshing.\n";
  PMP::split_long_edges(
      edges(_mesh), target_edge_length, _mesh,
      CGAL::parameters::edge_is_constrained_map(_edge_is_constrained_map));

  for (int iter = 0; iter < number_of_iterations; ++iter) {
    if (split_long_edges)
      if (verbose)
        std::cout << "Splitting long edges in iteration " << iter + 1 << ".\n";
    PMP::split_long_edges(
        edges(_mesh), target_edge_length, _mesh,
        CGAL::parameters::edge_is_constrained_map(_edge_is_constrained_map));
    if (verbose)
      std::cout << "Remeshing iteration " << iter + 1 << " of "
                << number_of_iterations << ".\n";
    PMP::isotropic_remeshing(
        faces(_mesh), target_edge_length, _mesh,
        CGAL::parameters::number_of_iterations(1) // one sub‑iteration per loop
            .edge_is_constrained_map(_edge_is_constrained_map)
            .protect_constraints(protect_constraints)
            .relax_constraints(relax_constraints));
  }

  if (verbose)
    std::cout << "Refined mesh → " << _mesh.number_of_vertices() << " V, "
              << _mesh.number_of_faces() << " F\n";
  if (!CGAL::is_valid_polygon_mesh(_mesh, verbose) && verbose)
    std::cout << "      ! mesh is not a valid polygon mesh after remeshing\n";
}

void TriMesh::reverseFaceOrientation() {
  // Reverse the face orientation of the mesh
  PMP::reverse_face_orientations(_mesh);
  if (!CGAL::is_valid_polygon_mesh(_mesh, true)) {
    std::cerr << "Error: Mesh is not valid after reversing face orientations."
              << std::endl;
  }
}

void TriMesh::cutWithSurface(TriMesh &clipper, bool verbose,
                             bool preserve_intersection,
                             bool preserve_intersection_clipper) {
  if (verbose) {
    std::cout << "Cutting mesh with surface." << std::endl;
  }
  bool intersection = PMP::do_intersect(_mesh, clipper._mesh);
  if (intersection) {
    // Clip tm with clipper
    if (verbose) {
      std::cout << "Clipping tm with clipper." << std::endl;
    }
    bool flag =
        PMP::clip(_mesh, clipper._mesh, CGAL::parameters::clip_volume(false));
  }
}

NumpyMesh TriMesh::save(double area_threshold,
                        double duplicate_vertex_threshold, bool verbose) {
  return export_mesh(_mesh, area_threshold, duplicate_vertex_threshold,
                     verbose);
}