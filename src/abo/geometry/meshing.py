"""Mesh generation for reflector geometries using Gmsh.

Generates surface meshes suitable for BEM computation from parametric
reflector definitions. Mesh density is controlled by elements-per-wavelength.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import gmsh
import numpy as np
from bempp_cl.api import Grid

if TYPE_CHECKING:
    from abo.geometry.reflectors import (
        ChebyshevProfile,
        EllipsoidalCap,
        SphericalCap,
    )


def mesh_spherical_cap(
    cap: SphericalCap,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> Grid:
    """Generate a surface mesh for a spherical cap reflector.

    The cap is centred at the origin with its axis along +z. The pole
    is at (0, 0, r) and the concave side faces -z (toward the
    approaching bat). The rim lies in the plane z = r*cos(theta).

    Parameters
    ----------
    cap : SphericalCap
        Reflector geometry.
    elements_per_wavelength : float
        Number of mesh elements per acoustic wavelength.
    frequency : float
        Design frequency in Hz (determines element size).
    speed_of_sound : float
        Speed of sound in m/s.

    Returns
    -------
    Grid
        Surface mesh compatible with bempp-cl.
    """
    wavelength = speed_of_sound / frequency
    element_size = wavelength / elements_per_wavelength

    r = cap.radius
    theta = cap.half_angle
    cut_z = r * np.cos(theta)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)

        # Sphere centred at origin
        gmsh.model.occ.addSphere(0, 0, 0, r, tag=1)

        # Box covering everything below z = cut_z
        box_extent = 4 * r
        gmsh.model.occ.addBox(
            -box_extent / 2,
            -box_extent / 2,
            -box_extent / 2,
            box_extent,
            box_extent,
            box_extent / 2 + cut_z,
            tag=2,
        )

        # Boolean cut: sphere minus box = cap above cut plane
        gmsh.model.occ.cut([(3, 1)], [(3, 2)])
        gmsh.model.occ.synchronize()

        # Mesh size control
        gmsh.option.setNumber(
            "Mesh.CharacteristicLengthMax", element_size,
        )
        gmsh.option.setNumber(
            "Mesh.CharacteristicLengthMin", element_size * 0.5,
        )
        gmsh.model.mesh.generate(2)

        # Identify the spherical surface (not the flat disc at the rim).
        # The flat disc has zero z-extent in its bounding box.
        spherical_tag = _find_spherical_surface()

        vertices, elements = _extract_surface_mesh(spherical_tag)
    finally:
        gmsh.finalize()

    return Grid(vertices, elements)


def _find_spherical_surface() -> int:
    """Find the Gmsh surface tag of the spherical cap (not the flat disc).

    Returns
    -------
    int
        Surface tag of the curved spherical cap surface.

    Raises
    ------
    RuntimeError
        If no curved surface is found.
    """
    for _, tag in gmsh.model.getEntities(dim=2):
        bb = gmsh.model.getBoundingBox(2, tag)
        z_range = bb[5] - bb[2]
        if z_range > 1e-10:
            return tag
    msg = "No curved surface found after boolean cut"
    raise RuntimeError(msg)


def _extract_surface_mesh(
    surface_tag: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract vertices and triangles from a specific Gmsh surface.

    Parameters
    ----------
    surface_tag : int
        Gmsh surface entity tag.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        vertices : float64 array of shape (3, N)
        elements : int32 array of shape (3, M)
    """
    # Get all nodes from the entire model (includes shared boundary nodes)
    all_tags, all_coords, _ = gmsh.model.mesh.getNodes()
    all_coords = all_coords.reshape(-1, 3)
    global_tag_to_coord = {
        int(t): all_coords[i] for i, t in enumerate(all_tags)
    }

    # Get triangles for this surface
    _, _, elem_node_tags = gmsh.model.mesh.getElements(
        dim=2, tag=surface_tag,
    )
    triangles = elem_node_tags[0].reshape(-1, 3).astype(int)

    # Collect unique node tags referenced by the surface triangles
    unique_tags = np.unique(triangles)

    # Build compact vertex array and remapping
    coords = np.array(
        [global_tag_to_coord[int(t)] for t in unique_tags],
        dtype=np.float64,
    )
    tag_to_idx = {int(t): i for i, t in enumerate(unique_tags)}

    # Remap element connectivity to 0-based indices
    elements = np.array(
        [[tag_to_idx[int(t)] for t in tri] for tri in triangles],
        dtype=np.int32,
    )

    return coords.T, elements.T


def mesh_ellipsoidal_cap(
    cap: EllipsoidalCap,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> Grid:
    """Generate a surface mesh for an ellipsoidal cap reflector."""
    raise NotImplementedError


def mesh_chebyshev_profile(
    profile: ChebyshevProfile,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> Grid:
    """Generate a surface mesh for a Chebyshev-profile reflector."""
    raise NotImplementedError
