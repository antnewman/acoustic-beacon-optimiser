"""Mesh generation for reflector geometries using Gmsh.

Generates surface meshes suitable for BEM computation from parametric
reflector definitions. Mesh density is controlled by elements-per-wavelength.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from abo.geometry.reflectors import ChebyshevProfile, EllipsoidalCap, SphericalCap


def mesh_spherical_cap(
    cap: SphericalCap,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> object:
    """Generate a surface mesh for a spherical cap reflector.

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
    bempp.api.Grid
        Surface mesh compatible with bempp-cl.
    """
    raise NotImplementedError


def mesh_ellipsoidal_cap(
    cap: EllipsoidalCap,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> object:
    """Generate a surface mesh for an ellipsoidal cap reflector."""
    raise NotImplementedError


def mesh_chebyshev_profile(
    profile: ChebyshevProfile,
    elements_per_wavelength: float = 10.0,
    frequency: float = 60_000.0,
    speed_of_sound: float = 343.0,
) -> object:
    """Generate a surface mesh for a Chebyshev-profile reflector."""
    raise NotImplementedError
