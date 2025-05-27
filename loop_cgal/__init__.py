from .loop_cgal import clip_surface, NumpyMesh
import pyvista as pv
import numpy as np
def clip_pyvista_polydata(
    surface_1: pv.PolyData,
    surface_2: pv.PolyData,
) -> pv.PolyData:
    """
    Clip two pyvista PolyData objects using the CGAL library.

    Parameters
    ----------
    surface_1 : pyvista.PolyData
        The first surface to be clipped.
    surface_2 : pyvista.PolyData
        The second surface to be used for clipping.

    Returns
    -------
    pyvista.PolyData
        The resulting clipped surface.
    """
    surface_1 = surface_1.triangulate()
    surface_2 = surface_2.triangulate()
    tm = NumpyMesh()
    tm.vertices = np.array(surface_1.points).copy()
    tm.triangles = surface_1.faces.reshape(-1, 4)[:, 1:].copy()
    clipper = NumpyMesh()
    clipper.vertices = np.array(surface_2.points).copy()
    clipper.triangles = surface_2.faces.reshape(-1, 4)[:, 1:].copy()
    mesh = clip_surface(
        tm,
        clipper,
    )
    return pv.PolyData.from_regular_faces(mesh.vertices, mesh.triangles)
