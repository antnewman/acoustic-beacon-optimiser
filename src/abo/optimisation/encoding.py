"""Parameter-vector encoding for reflector families.

Converts between dense numerical parameter vectors (as consumed by
CMA-ES and NSGA-II) and the typed reflector dataclasses. Each encoder
provides:

- `dim`: dimensionality of the parameter vector
- `default_bounds`: a sensible default box-bound list for the family
- `decode(x)`: parameter vector to reflector object
- `x0`: a default starting point
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from abo.geometry.reflectors import (
    ChebyshevProfile,
    EllipsoidalCap,
    SphericalCap,
)


@dataclass(frozen=True)
class SphericalCapEncoding:
    """Two-parameter encoding for a SphericalCap.

    Parameters are (radius, half_angle). Default bounds span
    radii 10-100 mm and half-angles 0.2 rad (~11 deg) to near-pi/2.
    """

    radius_min: float = 0.010
    radius_max: float = 0.100
    half_angle_min: float = 0.2
    half_angle_max: float = np.pi / 2 - 0.05

    @property
    def dim(self) -> int:
        return 2

    @property
    def default_bounds(self) -> list[tuple[float, float]]:
        return [
            (self.radius_min, self.radius_max),
            (self.half_angle_min, self.half_angle_max),
        ]

    @property
    def x0(self) -> NDArray[np.float64]:
        return np.array(
            [
                0.5 * (self.radius_min + self.radius_max),
                0.5 * (self.half_angle_min + self.half_angle_max),
            ],
            dtype=np.float64,
        )

    def decode(self, x: NDArray[np.float64]) -> SphericalCap:
        return SphericalCap(radius=float(x[0]), half_angle=float(x[1]))


@dataclass(frozen=True)
class EllipsoidalCapEncoding:
    """Three-parameter encoding for an EllipsoidalCap.

    Parameters are (semi_major, semi_minor, half_angle).
    """

    semi_major_min: float = 0.010
    semi_major_max: float = 0.100
    semi_minor_min: float = 0.005
    semi_minor_max: float = 0.100
    half_angle_min: float = 0.2
    half_angle_max: float = np.pi / 2 - 0.05

    @property
    def dim(self) -> int:
        return 3

    @property
    def default_bounds(self) -> list[tuple[float, float]]:
        return [
            (self.semi_major_min, self.semi_major_max),
            (self.semi_minor_min, self.semi_minor_max),
            (self.half_angle_min, self.half_angle_max),
        ]

    @property
    def x0(self) -> NDArray[np.float64]:
        return np.array(
            [
                0.5 * (self.semi_major_min + self.semi_major_max),
                0.5 * (self.semi_minor_min + self.semi_minor_max),
                0.5 * (self.half_angle_min + self.half_angle_max),
            ],
            dtype=np.float64,
        )

    def decode(self, x: NDArray[np.float64]) -> EllipsoidalCap:
        return EllipsoidalCap(
            semi_major=float(x[0]),
            semi_minor=float(x[1]),
            half_angle=float(x[2]),
        )


@dataclass(frozen=True)
class ChebyshevProfileEncoding:
    """N-parameter encoding for a ChebyshevProfile with fixed aperture.

    Parameters are the N Chebyshev coefficients (c_0, ..., c_{N-1}).
    Aperture radius is fixed at construction to decouple the acoustic
    shape optimisation from the aperture-size degree of freedom
    (which is typically handled separately as a size constraint).

    Default coefficient bounds (+/- coeff_bound) are generous enough
    to include most concave axisymmetric profiles once the concavity
    penalty is applied.
    """

    n_coefficients: int
    aperture_radius: float
    coeff_bound: float = 0.050

    @property
    def dim(self) -> int:
        return self.n_coefficients

    @property
    def default_bounds(self) -> list[tuple[float, float]]:
        return [(-self.coeff_bound, self.coeff_bound)] * self.n_coefficients

    @property
    def x0(self) -> NDArray[np.float64]:
        # Default: a linearly decreasing profile (depth 0.020 at pole,
        # 0 at rim), expressed in Chebyshev coefficients:
        # z(rho) = 0.010 - 0.010 * T_1(u). Higher-order coeffs zero.
        coeffs = np.zeros(self.n_coefficients, dtype=np.float64)
        coeffs[0] = 0.010
        if self.n_coefficients >= 2:
            coeffs[1] = -0.010
        return coeffs

    def decode(self, x: NDArray[np.float64]) -> ChebyshevProfile:
        return ChebyshevProfile(
            coefficients=np.asarray(x, dtype=np.float64),
            aperture_radius=self.aperture_radius,
        )
