"""Helmholtz BEM solver using bempp-cl.

Solves the exterior Helmholtz problem for a rigid (sound-hard) scatterer
using the Burton-Miller combined field integral equation (CFIE):

    (1/2 + D_k + alpha * H_k) phi = -(1/2 + D_k + alpha * H_k) p_inc

with alpha = i/k (standard coupling parameter).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def solve_scattering(
    grid: object,
    frequency: float,
    incidence_direction: NDArray[np.float64],
    speed_of_sound: float = 343.0,
    density: float = 1.225,
) -> object:
    """Solve the Helmholtz scattering problem for a rigid reflector.

    Parameters
    ----------
    grid : bempp.api.Grid
        Surface mesh of the reflector.
    frequency : float
        Excitation frequency in Hz.
    incidence_direction : NDArray[np.float64]
        Unit vector for the incident plane wave direction.
    speed_of_sound : float
        Speed of sound in m/s.
    density : float
        Air density in kg/m^3.

    Returns
    -------
    bempp.api.GridFunction
        Total field on the boundary surface.
    """
    raise NotImplementedError


def scattered_field(
    solution: object,
    evaluation_points: NDArray[np.float64],
    frequency: float,
    speed_of_sound: float = 343.0,
) -> NDArray[np.complex128]:
    """Evaluate the scattered pressure field at exterior points.

    Parameters
    ----------
    solution : bempp.api.GridFunction
        BEM solution (total surface field).
    evaluation_points : NDArray[np.float64]
        Points at which to evaluate, shape (3, N).
    frequency : float
        Frequency in Hz.
    speed_of_sound : float
        Speed of sound in m/s.

    Returns
    -------
    NDArray[np.complex128]
        Scattered pressure at each evaluation point.
    """
    raise NotImplementedError
