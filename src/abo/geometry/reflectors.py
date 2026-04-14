"""Parametric reflector definitions.

Defines families of axially symmetric concave reflectors:
- Spherical cap (2 parameters: radius, half-angle)
- Ellipsoidal cap (3 parameters: semi-major, semi-minor, half-angle)
- General Chebyshev profile (N parameters)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class SphericalCap:
    """Spherical cap reflector defined by radius and half-angle.

    Parameters
    ----------
    radius : float
        Sphere radius in metres.
    half_angle : float
        Half-angle of the cap in radians.
    """

    radius: float
    half_angle: float

    @property
    def depth(self) -> float:
        """Cap depth d = r(1 - cos(theta))."""
        return self.radius * (1.0 - np.cos(self.half_angle))

    @property
    def aperture(self) -> float:
        """Aperture radius a = r*sin(theta)."""
        return self.radius * np.sin(self.half_angle)

    @property
    def surface_area(self) -> float:
        """Surface area A = 2*pi*r^2*(1 - cos(theta))."""
        return 2.0 * np.pi * self.radius**2 * (1.0 - np.cos(self.half_angle))


@dataclass(frozen=True)
class EllipsoidalCap:
    """Ellipsoidal cap reflector.

    Parameters
    ----------
    semi_major : float
        Semi-major axis in metres.
    semi_minor : float
        Semi-minor axis in metres.
    half_angle : float
        Half-angle of the cap in radians.
    """

    semi_major: float
    semi_minor: float
    half_angle: float


@dataclass(frozen=True)
class ChebyshevProfile:
    """General axial profile defined by truncated Chebyshev series.

    z(rho) = sum_{n=0}^{N-1} c_n * T_n(2*rho/a - 1)

    Parameters
    ----------
    coefficients : NDArray[np.float64]
        Chebyshev coefficients c_n.
    aperture_radius : float
        Aperture radius a in metres.
    """

    coefficients: NDArray[np.float64]
    aperture_radius: float

    def profile(self, rho: NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluate the axial profile z(rho).

        Parameters
        ----------
        rho : NDArray[np.float64]
            Radial distances from axis, in [0, aperture_radius].

        Returns
        -------
        NDArray[np.float64]
            Axial heights z(rho).
        """
        raise NotImplementedError
