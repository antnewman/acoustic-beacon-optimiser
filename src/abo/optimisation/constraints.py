"""Optimisation constraints for reflector shape parameters.

Enforces geometric constraints:
- Concavity (dz/drho < 0 for rho in (0, a))
- Surface area bounds
- Physical dimension bounds
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def concavity_penalty(
    coefficients: NDArray[np.float64],
    aperture_radius: float,
    n_samples: int = 100,
) -> float:
    """Compute a penalty for concavity constraint violation.

    Returns 0 if the profile is concave, positive otherwise.

    Parameters
    ----------
    coefficients : NDArray[np.float64]
        Chebyshev profile coefficients.
    aperture_radius : float
        Aperture radius.
    n_samples : int
        Number of sample points for evaluation.

    Returns
    -------
    float
        Penalty value (0 = feasible).
    """
    raise NotImplementedError


def area_constraint(
    grid: object,
    max_area: float,
) -> float:
    """Compute surface area constraint violation.

    Returns 0 if area <= max_area, positive otherwise.

    Parameters
    ----------
    grid : bempp.api.Grid
        Reflector mesh.
    max_area : float
        Maximum allowed surface area in m^2.

    Returns
    -------
    float
        Constraint violation (0 = feasible).
    """
    raise NotImplementedError
