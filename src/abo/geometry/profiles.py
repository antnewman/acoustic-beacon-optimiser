"""Chebyshev and Fourier profile curve utilities.

Provides functions for evaluating, differentiating, and constraining
axial profile curves used in the general reflector parameterisation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


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
    raise NotImplementedError


def is_concave(
    coefficients: NDArray[np.float64],
    aperture_radius: float,
    n_samples: int = 100,
) -> bool:
    """Check whether a Chebyshev profile satisfies the concavity constraint.

    Requires dz/drho < 0 for rho in (0, a).

    Parameters
    ----------
    coefficients : NDArray[np.float64]
        Chebyshev coefficients.
    aperture_radius : float
        Aperture radius.
    n_samples : int
        Number of sample points for the check.

    Returns
    -------
    bool
        True if the profile is concave.
    """
    raise NotImplementedError
