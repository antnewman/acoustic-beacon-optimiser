"""Optimisation constraints for reflector shape parameters.

Penalty functions for:
- Concavity (dz/drho <= 0 for rho in (0, a))
- Surface-area bound (area <= max_area)
- Physical dimension bounds via bound-clamping (handled by the optimiser)
"""

from __future__ import annotations

import numpy as np
from bempp_cl.api import Grid
from numpy.typing import NDArray

from abo.acoustics.objectives import surface_area
from abo.geometry.profiles import chebyshev_profile_derivative


def concavity_penalty(
    coefficients: NDArray[np.float64],
    aperture_radius: float,
    n_samples: int = 100,
) -> float:
    """Penalty for concavity-constraint violation on a Chebyshev profile.

    The penalty is the integral of max(0, dz/drho) along the profile,
    normalised by the aperture radius. Returns 0 if the profile is
    concave (dz/drho <= 0 everywhere), positive otherwise.

    Parameters
    ----------
    coefficients : NDArray[np.float64]
        Chebyshev profile coefficients.
    aperture_radius : float
        Aperture radius in metres.
    n_samples : int
        Number of sample points for numerical integration.

    Returns
    -------
    float
        Non-negative penalty; 0 iff feasible.
    """
    eps = 1e-9
    rho = np.linspace(eps, aperture_radius - eps, n_samples)
    slope = chebyshev_profile_derivative(rho, coefficients, aperture_radius)
    violation = np.maximum(slope, 0.0)
    return float(np.trapezoid(violation, rho) / aperture_radius)


def area_constraint(grid: Grid, max_area: float) -> float:
    """Penalty for surface-area-constraint violation.

    Returns 0 if surface_area(grid) <= max_area, otherwise the excess
    area in m^2.

    Parameters
    ----------
    grid : Grid
        Reflector mesh.
    max_area : float
        Maximum allowed surface area in m^2.

    Returns
    -------
    float
        Non-negative penalty; 0 iff feasible.
    """
    excess = surface_area(grid) - max_area
    return float(max(0.0, excess))
