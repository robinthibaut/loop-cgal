from __future__ import annotations

from typing import Tuple

import numpy as np
import pyvista as pv

from ._loop_cgal import NumpyMesh, NumpyPlane, clip_plane, clip_surface, corefine_mesh
from ._loop_cgal import TriMesh as _TriMesh
from ._loop_cgal import verbose
from ._loop_cgal import set_verbose as set_verbose



def validate_numpy_mesh(mesh: NumpyMesh, mesh_name: str = "mesh") -> None:
    """Validate a NumpyMesh object to prevent crashes in C++ code.
    
    Parameters
    ----------
    mesh : NumpyMesh
        The mesh to validate
    mesh_name : str
        Name of the mesh for error messages
        
    Raises
    ------
    ValueError
        If the mesh is invalid
    """
    if not isinstance(mesh, NumpyMesh):
        raise ValueError(f"{mesh_name} must be a NumpyMesh object")
    
    # Check vertices
    if not hasattr(mesh, 'vertices') or mesh.vertices is None:
        raise ValueError(f"{mesh_name} vertices are None")
    
    vertices = np.asarray(mesh.vertices)
    if vertices.size == 0:
        raise ValueError(f"{mesh_name} has no vertices")
    
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError(f"{mesh_name} vertices must be Nx3 array, got shape {vertices.shape}")
    
    if not np.isfinite(vertices).all():
        raise ValueError(f"{mesh_name} vertices contain NaN or infinite values")
    
    # Check triangles
    if not hasattr(mesh, 'triangles') or mesh.triangles is None:
        raise ValueError(f"{mesh_name} triangles are None")
    
    triangles = np.asarray(mesh.triangles)
    if triangles.size == 0:
        raise ValueError(f"{mesh_name} has no triangles")
    
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError(f"{mesh_name} triangles must be Nx3 array, got shape {triangles.shape}")
    
    if triangles.dtype.kind not in ['i', 'u']:  # integer types
        raise ValueError(f"{mesh_name} triangle indices must be integers, got {triangles.dtype}")
    
    # Check triangle indices are within bounds
    max_vertex_index = vertices.shape[0] - 1
    if triangles.min() < 0:
        raise ValueError(f"{mesh_name} has negative triangle indices")
    
    if triangles.max() > max_vertex_index:
        raise ValueError(f"{mesh_name} triangle indices exceed vertex count (max index: {triangles.max()}, vertex count: {vertices.shape[0]})")
    
    # Check for degenerate triangles (same vertex used multiple times)
    for i, tri in enumerate(triangles):
        if len(set(tri)) != 3:
            raise ValueError(f"{mesh_name} triangle {i} is degenerate: vertices {tri}")


def validate_pyvista_polydata(surface: pv.PolyData, surface_name: str = "surface") -> None:
    """Validate a PyVista PolyData object.
    
    Parameters
    ----------
    surface : pv.PolyData
        The surface to validate
    surface_name : str
        Name of the surface for error messages
        
    Raises
    ------
    ValueError
        If the surface is invalid
    """
    if not isinstance(surface, pv.PolyData):
        raise ValueError(f"{surface_name} must be a pyvista.PolyData object")
    
    if surface.n_points == 0:
        raise ValueError(f"{surface_name} has no points")
    
    if surface.n_cells == 0:
        raise ValueError(f"{surface_name} has no cells")
    
    points = np.asarray(surface.points)
    if not np.isfinite(points).all():
        raise ValueError(f"{surface_name} points contain NaN or infinite values")


def validate_plane(plane_origin: np.ndarray, plane_normal: np.ndarray) -> None:
    """Validate plane parameters.
    
    Parameters
    ----------
    plane_origin : np.ndarray
        The plane origin point
    plane_normal : np.ndarray
        The plane normal vector
        
    Raises
    ------
    ValueError
        If the plane parameters are invalid
    """
    origin = np.asarray(plane_origin, dtype=np.float64)
    normal = np.asarray(plane_normal, dtype=np.float64)
    
    if origin.shape != (3,):
        raise ValueError(f"Plane origin must be a 3-element array, got shape {origin.shape}")
    
    if normal.shape != (3,):
        raise ValueError(f"Plane normal must be a 3-element array, got shape {normal.shape}")
    
    if not np.isfinite(origin).all():
        raise ValueError("Plane origin contains NaN or infinite values")
    
    if not np.isfinite(normal).all():
        raise ValueError("Plane normal contains NaN or infinite values")
    
    normal_length = np.linalg.norm(normal)
    if normal_length < 1e-10:
        raise ValueError(f"Plane normal vector is too small (length: {normal_length})")



class TriMesh(_TriMesh):
    """
    A class for handling triangular meshes using CGAL.
    
    Inherits from the base TriMesh class and provides additional functionality.
    """
    def __init__(self, surface: pv.PolyData):
        # Validate input surface
        validate_pyvista_polydata(surface, "input surface")
        
        # Triangulate to ensure we have triangular faces
        surface = surface.triangulate()
        
        # Extract vertices and triangles
        verts = np.array(surface.points, dtype=np.float64).copy()
        faces = surface.faces.reshape(-1, 4)[:, 1:].copy().astype(np.int32)
        
        # Additional validation on extracted data
        if verts.size == 0:
            raise ValueError("Surface has no vertices after triangulation")
        
        if faces.size == 0:
            raise ValueError("Surface has no triangular faces after triangulation")
        
        if not np.isfinite(verts).all():
            raise ValueError("Surface vertices contain NaN or infinite values")
        
        # Check triangle indices
        max_vertex_index = verts.shape[0] - 1
        if faces.min() < 0:
            raise ValueError("Surface has negative triangle indices")
        
        if faces.max() > max_vertex_index:
            raise ValueError(f"Surface triangle indices exceed vertex count (max index: {faces.max()}, vertex count: {verts.shape[0]})")
        
        # Check for degenerate triangles
        for i, tri in enumerate(faces):
            if len(set(tri)) != 3:
                raise ValueError(f"Surface triangle {i} is degenerate: vertices {tri}")
        
        super().__init__(verts, faces)
        
    def to_pyvista(self, area_threshold: float = 1e-6,  # this is the area threshold for the faces, if the area is smaller than this it will be removed
            duplicate_vertex_threshold: float = 1e-4,  # this is the threshold for duplicate vertices
            ) -> pv.PolyData:
        """
        Convert the TriMesh to a pyvista PolyData object.
        
        Returns
        -------
        pyvista.PolyData
            The converted PolyData object.
        """
        np_mesh = self.save(area_threshold, duplicate_vertex_threshold)
        vertices = np.array(np_mesh.vertices).copy()
        triangles = np.array(np_mesh.triangles).copy()
        return pv.PolyData.from_regular_faces(vertices, triangles)

def clip_pyvista_polydata_with_plane(
    surface: pv.PolyData,
    plane_origin: np.ndarray,
    plane_normal: np.ndarray,
    target_edge_length: float = 10.0,
    remesh_before_clipping: bool = True,
    remesh_after_clipping: bool = True,
    remove_degenerate_faces: bool = True,
    duplicate_vertex_threshold: float = 0.001,
    area_threshold: float = 0.0001,
    protect_constraints: bool = False,
    relax_constraints: bool = True,
) -> pv.PolyData:
    """
    Clip a pyvista PolyData object with a plane using the CGAL library.
    Parameters
    ----------
    surface : pyvista.PolyData

        The surface to be clipped.
    plane_origin : np.ndarray
        The origin point of the clipping plane.
    plane_normal : np.ndarray
        The normal vector of the clipping plane.
    target_edge_length : float, optional
        The target edge length for the remeshing process, by default 10.0
    remesh_before_clipping : bool, optional
        Whether to remesh the surface before clipping, by default True
    remesh_after_clipping : bool, optional
        Whether to remesh the surface after clipping, by default True
    remove_degenerate_faces : bool, optional
        Whether to remove degenerate faces from the resulting mesh, by default True
    duplicate_vertex_threshold : float, optional
        The threshold for merging duplicate vertices, by default 0.001
    area_threshold : float, optional
        The area threshold for removing small faces, by default 0.0001
    
    Returns
    -------
    pyvista.PolyData
        The resulting clipped surface.
    """
    # Validate inputs
    validate_pyvista_polydata(surface, "input surface")
    validate_plane(plane_origin, plane_normal)
    
    if target_edge_length <= 0:
        raise ValueError(f"target_edge_length must be positive, got {target_edge_length}")
    
    surface = surface.triangulate()
    tm = NumpyMesh()
    tm.vertices = np.array(surface.points, dtype=np.float64).copy()
    tm.triangles = surface.faces.reshape(-1, 4)[:, 1:].copy().astype(np.int32)
    
    # Validate the created NumpyMesh
    validate_numpy_mesh(tm, "surface mesh")
    
    plane = NumpyPlane()
    plane.origin = np.asarray(plane_origin, dtype=np.float64)
    plane.normal = np.asarray(plane_normal, dtype=np.float64)

    mesh = clip_plane(
        tm,
        plane,
        target_edge_length=target_edge_length,
        remesh_before_clipping=remesh_before_clipping,
        remesh_after_clipping=remesh_after_clipping,
        remove_degenerate_faces=remove_degenerate_faces,
        duplicate_vertex_threshold=duplicate_vertex_threshold,
        area_threshold=area_threshold,
        protect_constraints=protect_constraints,
        relax_constraints=relax_constraints,
    )
    return pv.PolyData.from_regular_faces(mesh.vertices, mesh.triangles)


def clip_pyvista_polydata(
    surface_1: pv.PolyData,
    surface_2: pv.PolyData,
    target_edge_length: float = 10.0,
    remesh_before_clipping: bool = True,
    remesh_after_clipping: bool = True,
    remove_degenerate_faces: bool = True,
    duplicate_vertex_threshold: float = 0.001,
    area_threshold: float = 0.0001,
    protect_constraints: bool = False,
    relax_constraints: bool = True,
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
    # Validate inputs
    validate_pyvista_polydata(surface_1, "surface_1")
    validate_pyvista_polydata(surface_2, "surface_2")
    
    if target_edge_length <= 0:
        raise ValueError(f"target_edge_length must be positive, got {target_edge_length}")
    
    surface_1 = surface_1.triangulate()
    surface_2 = surface_2.triangulate()
    tm = NumpyMesh()
    tm.vertices = np.array(surface_1.points, dtype=np.float64).copy()
    tm.triangles = surface_1.faces.reshape(-1, 4)[:, 1:].copy().astype(np.int32)
    clipper = NumpyMesh()
    clipper.vertices = np.array(surface_2.points, dtype=np.float64).copy()
    clipper.triangles = surface_2.faces.reshape(-1, 4)[:, 1:].copy().astype(np.int32)
    
    # Validate the created NumpyMesh objects
    validate_numpy_mesh(tm, "surface_1 mesh")
    validate_numpy_mesh(clipper, "surface_2 mesh")
    mesh = clip_surface(
        tm,
        clipper,
        target_edge_length=target_edge_length,
        remesh_before_clipping=remesh_before_clipping,
        remesh_after_clipping=remesh_after_clipping,
        remove_degenerate_faces=remove_degenerate_faces,
        duplicate_vertex_threshold=duplicate_vertex_threshold,
        area_threshold=area_threshold,
        protect_constraints=protect_constraints,
        relax_constraints=relax_constraints,
    )
    return pv.PolyData.from_regular_faces(mesh.vertices, mesh.triangles)


def corefine_pyvista_polydata(
    surface_1: pv.PolyData,
    surface_2: pv.PolyData,
    target_edge_length: float = 10.0,
    duplicate_vertex_threshold: float = 0.001,
    area_threshold: float = 0.0001,
    number_of_iterations: int = 10,
    protect_constraints: bool = True,
    relax_constraints: bool = True,
) -> Tuple[pv.PolyData, pv.PolyData]:
    """
    Corefine two pyvista PolyData objects using the CGAL library.

    Parameters
    ----------
    surface_1 : pyvista.PolyData
        The first surface to be cored.
    surface_2 : pyvista.PolyData
        The second surface to be used for cording.

    Returns
    -------
    pyvista.PolyData
        The resulting cored surface.
    """
    # Validate inputs
    validate_pyvista_polydata(surface_1, "surface_1")
    validate_pyvista_polydata(surface_2, "surface_2")
    
    if target_edge_length <= 0:
        raise ValueError(f"target_edge_length must be positive, got {target_edge_length}")
    
    if number_of_iterations <= 0:
        raise ValueError(f"number_of_iterations must be positive, got {number_of_iterations}")
    
    surface_1 = surface_1.triangulate()
    surface_2 = surface_2.triangulate()
    tm1 = NumpyMesh()
    tm1.vertices = np.array(surface_1.points, dtype=np.float64).copy()
    tm1.triangles = surface_1.faces.reshape(-1, 4)[:, 1:].copy().astype(np.int32)
    tm2 = NumpyMesh()
    tm2.vertices = np.array(surface_2.points, dtype=np.float64).copy()
    tm2.triangles = surface_2.faces.reshape(-1, 4)[:, 1:].copy().astype(np.int32)
    
    # Validate the created NumpyMesh objects
    validate_numpy_mesh(tm1, "surface_1 mesh")
    validate_numpy_mesh(tm2, "surface_2 mesh")

    tm1, tm2 = corefine_mesh(
        tm1,
        tm2,
        target_edge_length=target_edge_length,
        duplicate_vertex_threshold=duplicate_vertex_threshold,
        area_threshold=area_threshold,
        number_of_iterations=number_of_iterations,
        relax_constraints=relax_constraints,
        protect_constraints=protect_constraints,
    )
    return (
        pv.PolyData.from_regular_faces(tm1.vertices, tm1.triangles),
        pv.PolyData.from_regular_faces(tm2.vertices, tm2.triangles),
    )
