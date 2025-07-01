import loop_cgal
import numpy as np
import pyvista as pv
from LoopStructural.datatypes import BoundingBox
def test():
    bb = BoundingBox(np.zeros(3), np.ones(3))
    grid = bb.structured_grid().vtk()
    print(grid)
    grid["scalars"] = grid.points[:,0]
    surface = grid.contour([.5])
    surface_1_tri = surface.faces.reshape(-1, 4)[:, 1:].copy()
    surface_1_verts = surface.points.copy()
    grid["scalars"] = grid.points[:,1] * 0.5 + np.cos(grid.points[:,0]) * 0.5
    surface_2 = grid.contour([.45])
    surface_2_tri = surface_2.faces.reshape(-1, 4)[:, 1:].copy()
    surface_2_verts = surface_2.points.copy()
    surface_cgal = loop_cgal.TriMesh(surface)
    # surface_cgal.remesh(target_edge_length=0.02, verbose=True,protect_constraints=True, relax_constraints=False, number_of_iterations=1,split_long_edges=False)
    surface_cgal_2 = loop_cgal.TriMesh(surface_2)
    # surface_cgal_2.remesh(target_edge_length=0.02, verbose=True,protect_constraints=True, relax_constraints=False, number_of_iterations=1,split_long_edges=False)
    print("clipping 1")
    surface_cgal.cut_with_surface(surface_cgal_2,verbose=True, preserve_intersection=False, preserve_intersection_clipper=False)

    print("reverse face orientation")
    surface_cgal_2.reverse_face_orientation()
    print("clipping 3")
    surface_cgal_3 = loop_cgal.TriMesh(surface)
    # surface_cgal_3.remesh(target_edge_length=0.02, verbose=True,protect_constraints=True, relax_constraints=False, number_of_iterations=1,split_long_edges=False)
    surface_cgal_3.cut_with_surface(surface_cgal_2)


    print(surface_cgal.to_pyvista() )
    # mesh = surface_cgal.save(area_threshold=0.0001, duplicate_vertex_threshold=0.000001, verbose=True)
    # print(mesh.vertices.shape, mesh.triangles.shape)
    # surface_clipped = loop_cgal.clip_pyvista_polydata(
    #     surface, surface_2,
    #     target_edge_length=0.2,
    #     remesh_before_clipping=True,
    #     remesh_after_clipping=True,
    #     remove_degenerate_faces=False,
    #     duplicate_vertex_threshold=0.000001,
    #     area_threshold=0.0001,
    #     protect_constraints=True,
    #     relax_constraints=False,
    #     verbose=True
    # )  # surface_1_tri, surface_1_verts, surface_2_tri, surface_2_verts
    # # )
    # print(surface_clipped)
    # # print(mesh.vertices.shape)
    # # print(surface_1_verts.shape)
    # # print(mesh.triangles)
    # # print(mesh.triangles, mesh.vertices)
    # # surface_clipped = pv.PolyData.from_regular_faces(mesh.vertices, mesh.triangles)
    surface_cgal.to_pyvista().save('f1.ply')
    surface_cgal_2.to_pyvista().save('f2.ply')
    surface_cgal_3.to_pyvista().save('f3.ply')
    # p = pv.Plotter()
    # # p.add_mesh(surface)
    # p.add_mesh(surface_cgal.to_pyvista(), color="blue", show_edges=True)
    # p.add_mesh(surface_cgal_2.to_pyvista(), color="red", show_edges=True)
    # # p.add_points(mesh.vertices, render_points_as_spheres=True, point_size=5)
    # p.show()
if __name__ == "__main__":
    test()
