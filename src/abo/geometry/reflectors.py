"""Parametric reflector definitions.

Defines families of axially symmetric concave reflectors:
- Spherical cap (2 parameters: radius, half-angle)
- Ellipsoidal cap (3 parameters: semi-major, semi-minor, half-angle)
- General Chebyshev profile (N parameters)

All geometries are oriented with axis along +z, pole at positive z,
concave side facing -z.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy import integrate

from abo.geometry.profiles import (
    chebyshev_profile,
    chebyshev_profile_derivative,
    is_concave,
)


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
    """Oblate/prolate ellipsoidal cap reflector.

    The underlying surface is the ellipsoid (rho/a)^2 + (z/c)^2 = 1,
    where a = semi_major is the radial semi-axis and c = semi_minor is
    the axial semi-axis. The cap is the portion with polar angle in
    [0, half_angle] (measured from the +z pole), parameterised by

        rho(t) = a * sin(t)
        z(t)   = c * cos(t)

    for t in [0, half_angle]. When a = c this reduces exactly to a
    spherical cap of radius a.

    Parameters
    ----------
    semi_major : float
        Radial semi-axis in metres (aperture-direction).
    semi_minor : float
        Axial semi-axis in metres (depth-direction).
    half_angle : float
        Half-angle of the cap in radians, measured from the pole.
    """

    semi_major: float
    semi_minor: float
    half_angle: float

    @property
    def depth(self) -> float:
        """Cap depth d = c * (1 - cos(theta))."""
        return self.semi_minor * (1.0 - np.cos(self.half_angle))

    @property
    def aperture(self) -> float:
        """Aperture radius a_rim = semi_major * sin(theta)."""
        return self.semi_major * np.sin(self.half_angle)

    @property
    def surface_area(self) -> float:
        """Surface area via numerical integration.

        A = integral_0^theta 2*pi*rho(t) * sqrt((d rho/dt)^2 + (dz/dt)^2) dt
          = integral_0^theta 2*pi*a*sin(t) *
                sqrt(a^2*cos^2(t) + c^2*sin^2(t)) dt
        """
        a = self.semi_major
        c = self.semi_minor

        def integrand(t: float) -> float:
            return (
                2.0 * np.pi * a * np.sin(t)
                * np.sqrt(a**2 * np.cos(t) ** 2 + c**2 * np.sin(t) ** 2)
            )

        area, _ = integrate.quad(integrand, 0.0, self.half_angle)
        return float(area)


@dataclass(frozen=True)
class ChebyshevProfile:
    """Axisymmetric reflector defined by a truncated Chebyshev series.

    The axial profile is

        z(rho) = sum_{n=0}^{N-1} c_n * T_n(2*rho/a - 1)

    where T_n is the Chebyshev polynomial of the first kind. The mapping
    (2*rho/a - 1) takes rho in [0, a] to u in [-1, 1], so T_n terms are
    orthogonal over the aperture. The pole is at rho = 0 (z = z(0)) and
    the rim is at rho = aperture_radius (z = z(a)).

    Parameters
    ----------
    coefficients : NDArray[np.float64]
        Chebyshev coefficients c_0, ..., c_{N-1}.
    aperture_radius : float
        Aperture radius a in metres.
    """

    coefficients: NDArray[np.float64]
    aperture_radius: float

    @property
    def depth(self) -> float:
        """Depth at the pole: z(0) - z(a)."""
        return float(
            chebyshev_profile(
                np.asarray([0.0]), self.coefficients, self.aperture_radius,
            )[0]
            - chebyshev_profile(
                np.asarray([self.aperture_radius]),
                self.coefficients,
                self.aperture_radius,
            )[0]
        )

    @property
    def aperture(self) -> float:
        """Aperture radius (identity)."""
        return self.aperture_radius

    def is_concave(self, n_samples: int = 100) -> bool:
        """Whether the profile satisfies dz/drho < 0 on (0, a)."""
        return is_concave(
            self.coefficients, self.aperture_radius, n_samples,
        )

    def profile(
        self, rho: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Evaluate the axial profile z(rho)."""
        return chebyshev_profile(
            rho, self.coefficients, self.aperture_radius,
        )

    @property
    def surface_area(self) -> float:
        """Surface area of revolution.

        A = integral_0^a 2*pi*rho * sqrt(1 + (dz/drho)^2) drho
        """

        def integrand(rho: float) -> float:
            dzdrho = chebyshev_profile_derivative(
                np.asarray([rho]),
                self.coefficients,
                self.aperture_radius,
            )[0]
            return 2.0 * np.pi * rho * np.sqrt(1.0 + dzdrho**2)

        area, _ = integrate.quad(
            integrand, 0.0, self.aperture_radius,
        )
        return float(area)
