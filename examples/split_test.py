from __future__ import annotations

import loop_cgal
import numpy as np
from LoopStructural.datatypes import BoundingBox


def test():
    bb = BoundingBox(np.zeros(3), np.ones(3))
    grid = bb.structured_grid().vtk()
    print(grid)
    grid["scalars"] = grid.points[:,0]
    surface = grid.contour([.5])
    surface_1_tri = surface.faces.reshape(-1, 4)[:, 1:].copy()
    surface_1_verts = surface.points.copy()
    grid["scalars"] = grid.points[:,1]
    surface_2 = grid.contour([.5])
    surface_2_tri = surface_2.faces.reshape(-1, 4)[:, 1:].copy()
    surface_2_verts = surface_2.points.copy()

    refined1, refined2  = loop_cgal.corefine_pyvista_polydata(surface, surface_2,
    target_edge_length=0.10,
    duplicate_vertex_threshold=0.000001,
    area_threshold=0.0001,
    number_of_iterations=10,
    protect_constraints = True,
    relax_constraints = True,
    verbose = True,)
# ) -> Tuple[pv.PolyData, pv.PolyData]:
#     surface_clipped = loop_cgal.corefine(
#         surface, surface_2,target_edge_length=0.1,
#         remesh_before_clipping=True,
#         remesh_after_clipping=False,
#         remove_degenerate_faces=False,
#         duplicate_vertex_threshold=0.000001,
#         area_threshold=0.0001,
#         protect_constraints=True,
#         relax_constraints=False,
#         verbose=True
#     )  # surface_1_tri, surface_1_verts, surface_2_tri, surface_2_verts
#     # )
    refined1.save('refined1.vtk')
    refined2.save('refined2.vtk')
    surface.save('surface_1.vtk')
    surface_2.save('surface_2.vtk')
    # print(surface_clipped)
    # print(mesh.vertices.shape)
    # print(surface_1_verts.shape)
    # print(mesh.triangles)
    # print(mesh.triangles, mesh.vertices)
    # surface_clipped = pv.PolyData.from_regular_faces(mesh.vertices, mesh.triangles)
    # p = pv.Plotter()
    # # p.add_mesh(surface)
    # p.add_mesh(surface_2, color="blue", show_edges=True)
    # p.add_mesh(surface_clipped, color="red", show_edges=True)
    # # p.add_points(mesh.vertices, render_points_as_spheres=True, point_size=5)
    # p.show()
if __name__ == "__main__":
    test()
