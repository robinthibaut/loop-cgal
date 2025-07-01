#include "weld.h"

#include "clip.h"          // load_mesh(), export_mesh(), refine helpers …
#include "meshutils.h"     // collect_border_edges()

// ➊  CGAL / PMP headers ------------------------------------------------------
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>      // <— isotropic_remeshing
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>                // is_closed()
#include <utility>                                   // std::move

namespace PMP = CGAL::Polygon_mesh_processing;
using SurfaceMesh = TriangleMesh;   // alias for brevity

// ---------------------------------------------------------------------------
//   Helper ‑‑ append src faces & stitch borders
// ---------------------------------------------------------------------------
static void append_mesh(SurfaceMesh& dst, const SurfaceMesh& src)
{
  CGAL::copy_face_graph(src, dst);
  PMP::stitch_borders(dst);
  PMP::merge_duplicated_vertices_in_boundary_cycles(dst);
}

// ---------------------------------------------------------------------------
//   Helper ‑‑ constrained isotropic remesh
// ---------------------------------------------------------------------------
static void robust_remesh(SurfaceMesh& m,
                          double  L,
                          int     n_iter,
                          bool    protect,
                          bool    relax,
                          bool    verbose)
{
  std::set<SurfaceMesh::Edge_index> border = collect_border_edges(m);
  PMP::isotropic_remeshing(
          faces(m), L, m,
          CGAL::parameters::number_of_iterations(n_iter)
              .edge_is_constrained_map(
                  CGAL::make_boolean_property_map(border))
              .protect_constraints(protect)
              .relax_constraints(relax));
  if(verbose)
      std::cout << "    ↳ remeshed to " << m.number_of_faces() << " faces\n";
}

// ---------------------------------------------------------------------------
//   weld_meshes  (public API)
// ---------------------------------------------------------------------------
NumpyMesh weld_meshes(const std::vector<NumpyMesh>& meshes,
                      double   target_edge_length,
                      double   duplicate_vertex_threshold,
                      double   area_threshold,
                      int      remesh_iterations,
                      bool     protect_constraints,
                      bool     relax_constraints,
                      bool     verbose)
{
  if(meshes.empty())
      throw std::runtime_error("[weld] empty input list.");

  // 1. seed
  SurfaceMesh out = load_mesh(meshes.front(), verbose);
  PMP::remove_isolated_vertices(out);

  // 2. loop over the remaining meshes
  for(std::size_t i = 1; i < meshes.size(); ++i)
  {
      if(verbose) std::cout << "[weld] === mesh " << i+1 << " / "
                            << meshes.size() << " ===\n";
      SurfaceMesh nxt = load_mesh(meshes[i], verbose);
      PMP::remove_isolated_vertices(nxt);

      // (a) optional pre‑remesh
      robust_remesh(out, target_edge_length, remesh_iterations,
                    protect_constraints, relax_constraints, verbose);
      robust_remesh(nxt, target_edge_length, remesh_iterations,
                    protect_constraints, relax_constraints, verbose);

      // (b) corefine so intersection curves coincide
      PMP::corefine(out, nxt);

      // (c)  decide: Boolean union or simple glue
      bool closed_out  = CGAL::is_closed(out);
      bool closed_nxt  = CGAL::is_closed(nxt);

      if(closed_out && closed_nxt)
      {
          SurfaceMesh tmp;
          if(!PMP::corefine_and_compute_union(out, nxt, tmp))
              throw std::runtime_error("[weld] union failed.");
          out = std::move(tmp);   // Surface_mesh has move‑assignment
      }
      else
      {
          append_mesh(out, nxt);
      }

      // (d) deduplicate border vertices again
      PMP::merge_duplicated_vertices_in_boundary_cycles(out);
  }

      // --- 3.  final seam‑welding pass  -------------------------------
    PMP::stitch_borders(out);                                                // join identical halfedges
    PMP::merge_duplicated_vertices_in_boundary_cycles(out);                  // collapse equal points

    // 3. post‑clean
    PMP::remove_almost_degenerate_faces(faces(out), out);
    PMP::remove_isolated_vertices(out);

  // 4. export
  return export_mesh(out, area_threshold, duplicate_vertex_threshold, verbose);
}
