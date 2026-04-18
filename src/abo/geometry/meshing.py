"""Mesh generation for reflector geometries using Gmsh.

Generates surface meshes suitable for BEM computation from parametric
reflector definitions. Mesh density is controlled by elements-per-wavelength.

Spherical and ellipsoidal caps are meshed via Gmsh OCC boolean operations.
Chebyshev-profile reflectors are meshed by direct structured triangulation
of the surface of revolution (Gmsh OCC has no primitive for arbitrary
profile curves).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import gmsh
import numpy as np
from bempp_cl.api import Grid

from abo.geometry.profiles import chebyshev_profile

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
        spherical_tag = _find_curved_surface()

        vertices, elements = _extract_surface_mesh(spherical_tag)
    finally:
        gmsh.finalize()

    return Grid(vertices, elements)


def mesh_ellipsoidal_cap(
    cap: EllipsoidalCap,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> Grid:
    """Generate a surface mesh for an ellipsoidal cap reflector.

    Built by dilating a sphere to an ellipsoid (Gmsh OCC has no direct
    ellipsoid primitive), then cutting below z = c*cos(half_angle).

    Parameters
    ----------
    cap : EllipsoidalCap
        Reflector geometry.
    elements_per_wavelength : float
        Number of mesh elements per acoustic wavelength.
    frequency : float
        Design frequency in Hz.
    speed_of_sound : float
        Speed of sound in m/s.

    Returns
    -------
    Grid
        Surface mesh compatible with bempp-cl.
    """
    wavelength = speed_of_sound / frequency
    element_size = wavelength / elements_per_wavelength

    a = cap.semi_major
    c = cap.semi_minor
    theta = cap.half_angle
    cut_z = c * np.cos(theta)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)

        # Start with a sphere of radius a, then dilate in z by c/a to
        # produce an ellipsoid with semi-axes (a, a, c).
        gmsh.model.occ.addSphere(0, 0, 0, a, tag=1)
        gmsh.model.occ.dilate([(3, 1)], 0, 0, 0, 1.0, 1.0, c / a)

        # Box covering everything below z = cut_z
        box_extent = 4.0 * max(a, c)
        gmsh.model.occ.addBox(
            -box_extent / 2,
            -box_extent / 2,
            -box_extent / 2,
            box_extent,
            box_extent,
            box_extent / 2 + cut_z,
            tag=2,
        )

        gmsh.model.occ.cut([(3, 1)], [(3, 2)])
        gmsh.model.occ.synchronize()

        gmsh.option.setNumber(
            "Mesh.CharacteristicLengthMax", element_size,
        )
        gmsh.option.setNumber(
            "Mesh.CharacteristicLengthMin", element_size * 0.5,
        )
        gmsh.model.mesh.generate(2)

        curved_tag = _find_curved_surface()
        vertices, elements = _extract_surface_mesh(curved_tag)
    finally:
        gmsh.finalize()

    return Grid(vertices, elements)


def mesh_chebyshev_profile(
    profile: ChebyshevProfile,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> Grid:
    """Generate a surface mesh for a Chebyshev-profile reflector.

    Uses direct structured triangulation of the surface of revolution.
    A polar grid of (rho, phi) samples produces quadrilateral patches
    split into triangles, with a triangle fan at the axial pole.

    Parameters
    ----------
    profile : ChebyshevProfile
        Reflector geometry.
    elements_per_wavelength : float
        Number of mesh elements per acoustic wavelength.
    frequency : float
        Design frequency in Hz.
    speed_of_sound : float
        Speed of sound in m/s.

    Returns
    -------
    Grid
        Surface mesh compatible with bempp-cl.
    """
    wavelength = speed_of_sound / frequency
    element_size = wavelength / elements_per_wavelength

    a = profile.aperture_radius
    coeffs = profile.coefficients

    # Determine grid resolution from target element size
    n_radial = max(int(np.ceil(a / element_size)), 4)
    n_azim = max(int(np.ceil(2.0 * np.pi * a / element_size)), 8)

    # Radial positions excluding the pole (which is a single apex vertex)
    rho = np.linspace(a / n_radial, a, n_radial)
    z_rho = chebyshev_profile(rho, coeffs, a)

    # Azimuthal positions (non-closing; connectivity wraps)
    phi = np.linspace(0.0, 2.0 * np.pi, n_azim, endpoint=False)

    # Build vertex array: apex first, then ring-by-ring
    z_apex = float(chebyshev_profile(np.asarray([0.0]), coeffs, a)[0])
    n_verts = 1 + n_radial * n_azim
    vertices = np.empty((3, n_verts), dtype=np.float64)
    vertices[:, 0] = [0.0, 0.0, z_apex]

    for i, r_i in enumerate(rho):
        for j, phi_j in enumerate(phi):
            idx = 1 + i * n_azim + j
            vertices[0, idx] = r_i * np.cos(phi_j)
            vertices[1, idx] = r_i * np.sin(phi_j)
            vertices[2, idx] = z_rho[i]

    def ring_index(i: int, j: int) -> int:
        """Vertex index on radial ring i, azimuthal step j (mod n_azim)."""
        return 1 + i * n_azim + (j % n_azim)

    # Triangle fan at apex
    fan_tris = []
    for j in range(n_azim):
        fan_tris.append([0, ring_index(0, j), ring_index(0, j + 1)])

    # Quad strips between rings, each split into two triangles
    strip_tris = []
    for i in range(n_radial - 1):
        for j in range(n_azim):
            v00 = ring_index(i, j)
            v01 = ring_index(i, j + 1)
            v10 = ring_index(i + 1, j)
            v11 = ring_index(i + 1, j + 1)
            strip_tris.append([v00, v10, v11])
            strip_tris.append([v00, v11, v01])

    elements = np.array(fan_tris + strip_tris, dtype=np.int32).T

    return Grid(vertices, elements)


def _find_curved_surface() -> int:
    """Find the Gmsh surface tag with non-zero z-extent.

    The boolean cut of a closed shape with a cutting box leaves
    typically two surfaces: the curved outer surface and the flat
    disc at the cut plane. The flat disc has zero z-range.

    Returns
    -------
    int
        Surface tag of the curved surface.

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
