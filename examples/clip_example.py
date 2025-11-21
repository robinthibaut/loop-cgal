import loop_cgal
import numpy as np
import pyvista as pv
from LoopStructural.datatypes import BoundingBox


def test():
    """Example demonstrating surface clipping using the TriMesh API."""
    # Create bounding box and generate test surfaces
    bb = BoundingBox(np.zeros(3), np.ones(3))
    grid = bb.structured_grid().vtk()

    # Create first surface (X = 0.5 plane)
    grid["scalars"] = grid.points[:, 0]
    surface_1 = grid.contour([0.5])

    # Create second surface (Y = 0.5 plane) - this will clip the first surface
    grid["scalars"] = grid.points[:, 1]
    surface_2 = grid.contour([0.5])

    # Convert PyVista surfaces to TriMesh objects
    print("Creating TriMesh objects...")
    mesh_1 = loop_cgal.TriMesh(surface_1)
    mesh_2 = loop_cgal.TriMesh(surface_2)

    # Remesh before clipping for better quality
    print("Remeshing surfaces...")
    mesh_1.remesh(
        target_edge_length=0.2,
        remove_degenerate_faces=False,
        protect_constraints=True,
        relax_constraints=False
    )
    mesh_2.remesh(
        target_edge_length=0.2,
        remove_degenerate_faces=False,
        protect_constraints=True,
        relax_constraints=False
    )

    # Perform clipping: mesh_1 will be clipped by mesh_2
    print("Clipping surface_1 with surface_2...")
    mesh_1.cutWithSurface(
        mesh_2,
        preserve_intersection=True,
        preserve_intersection_clipper=False,
        use_exact_kernel=False
    )

    # Remesh after clipping to clean up
    print("Remeshing after clipping...")
    mesh_1.remesh(
        target_edge_length=0.2,
        remove_degenerate_faces=False,
        protect_constraints=True,
        relax_constraints=False
    )

    # Convert back to PyVista for visualization
    surface_clipped = mesh_1.to_pyvista(
        area_threshold=0.0001,
        duplicate_vertex_threshold=0.000001
    )

    print(f"Original surface_1: {surface_1.n_points} points, {surface_1.n_faces} faces")
    print(f"Original surface_2: {surface_2.n_points} points, {surface_2.n_faces} faces")
    print(f"Clipped surface: {surface_clipped.n_points} points, {surface_clipped.n_faces} faces")

    # Visualize results
    p = pv.Plotter()
    p.add_mesh(surface_2, color="blue", opacity=0.5, show_edges=True, label="Clipper (surface_2)")
    p.add_mesh(surface_clipped, color="red", show_edges=True, label="Clipped (surface_1)")
    p.add_legend()
    p.show()


if __name__ == "__main__":
    test()
