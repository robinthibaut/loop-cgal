#pragma once
#include "numpymesh.h"          // your existing struct
#include <vector>

/**  Corefine + weld an arbitrary list of triangulated surface meshes.
 *
 *  – Intersection curves are corefined so every mesh shares the very
 *    same vertices/edges along them.
 *  – Closed components that meet another closed component are Boolean‑unioned;
 *    open patches remain open.
 *  – Duplicate vertices closer than `duplicate_vertex_threshold`
 *    are collapsed.
 *  – Facet orientation is preserved.
 */
NumpyMesh weld_meshes(const std::vector<NumpyMesh>& meshes,
                      double   target_edge_length         = 10.0,
                      double   duplicate_vertex_threshold = 1e-6,
                      double   area_threshold             = 1e-6,
                      int      remesh_iterations          = 3,
                      bool     protect_constraints        = false,
                      bool     relax_constraints          = true,
                      bool     verbose                    = false);
