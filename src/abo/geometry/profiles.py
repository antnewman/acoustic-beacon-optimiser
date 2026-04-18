"""Chebyshev profile curve utilities.

Provides functions for evaluating, differentiating, and constraining
axial profile curves used in the general reflector parameterisation:

    z(rho) = sum_{n=0}^{N-1} c_n * T_n(u(rho)),   u(rho) = 2*rho/a - 1

where T_n is the Chebyshev polynomial of the first kind and a is the
aperture radius.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial import chebyshev as npc
from numpy.typing import NDArray


def _normalise(
    rho: NDArray[np.float64], aperture_radius: float,
) -> NDArray[np.float64]:
    """Map rho in [0, a] to u in [-1, 1]."""
    return 2.0 * rho / aperture_radius - 1.0


def chebyshev_profile(
    rho: NDArray[np.float64],
    coefficients: NDArray[np.float64],
    aperture_radius: float,
) -> NDArray[np.float64]:
    """Evaluate a Chebyshev-series axial profile.

    z(rho) = sum_{n=0}^{N-1} c_n * T_n(2*rho/a - 1)

    Parameters
    ----------
    rho : NDArray[np.float64]
        Radial coordinates in [0, aperture_radius].
    coefficients : NDArray[np.float64]
        Chebyshev coefficients.
    aperture_radius : float
        Aperture radius.

    Returns
    -------
    NDArray[np.float64]
        Profile heights z(rho).
    """
    u = _normalise(np.asarray(rho, dtype=np.float64), aperture_radius)
    result = npc.chebval(u, np.asarray(coefficients, dtype=np.float64))
    return np.asarray(result, dtype=np.float64)


def chebyshev_profile_derivative(
    rho: NDArray[np.float64],
    coefficients: NDArray[np.float64],
    aperture_radius: float,
) -> NDArray[np.float64]:
    """Evaluate dz/drho for a Chebyshev profile.

    Parameters
    ----------
    rho : NDArray[np.float64]
        Radial coordinates in [0, aperture_radius].
    coefficients : NDArray[np.float64]
        Chebyshev coefficients.
    aperture_radius : float
        Aperture radius.

    Returns
    -------
    NDArray[np.float64]
        Slopes dz/drho at each rho.
    """
    coeffs = np.asarray(coefficients, dtype=np.float64)
    deriv_coeffs = npc.chebder(coeffs, m=1)
    u = _normalise(np.asarray(rho, dtype=np.float64), aperture_radius)
    # Chain rule: dz/drho = (dz/du) * (du/drho) = (dz/du) * (2/a)
    dzdu = npc.chebval(u, deriv_coeffs)
    return np.asarray(dzdu, dtype=np.float64) * (2.0 / aperture_radius)


def is_concave(
    coefficients: NDArray[np.float64],
    aperture_radius: float,
    n_samples: int = 100,
) -> bool:
    """Check whether a Chebyshev profile satisfies the concavity constraint.

    Concavity (in the reflector sense, pole above rim) requires
    dz/drho <= 0 for all rho in (0, a).

    Parameters
    ----------
    coefficients : NDArray[np.float64]
        Chebyshev coefficients.
    aperture_radius : float
        Aperture radius.
    n_samples : int
        Number of sample points for the check (uniform in rho).

    Returns
    -------
    bool
        True if dz/drho <= 0 everywhere tested.
    """
    # Exclude the strict endpoints to avoid numerical edge effects
    eps = 1e-9
    rho = np.linspace(eps, aperture_radius - eps, n_samples)
    slopes = chebyshev_profile_derivative(
        rho, coefficients, aperture_radius,
    )
    return bool(np.all(slopes <= 0.0))
