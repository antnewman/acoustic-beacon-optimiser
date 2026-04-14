"""Biologically meaningful objective functions for reflector optimisation.

- Integrated Conspicuousness (IC): frequency- and angle-weighted backscatter TS
- Catchment Volume (CV): spatial volume exceeding detection threshold
- Surface Area (SA): structural/photosynthetic cost proxy
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def integrated_conspicuousness(
    ts_matrix: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    angles: NDArray[np.float64],
    call_spectrum: NDArray[np.float64],
) -> float:
    """Compute integrated conspicuousness.

    IC = integral W(f) * TS(f, theta) * sin(theta) d(theta) df

    Parameters
    ----------
    ts_matrix : NDArray[np.float64]
        Target strength in dB, shape (n_freq, n_angle).
    frequencies : NDArray[np.float64]
        Frequency array in Hz.
    angles : NDArray[np.float64]
        Incidence angle array in radians.
    call_spectrum : NDArray[np.float64]
        Weighting function W(f) matching bat call power spectrum.

    Returns
    -------
    float
        Integrated conspicuousness value.
    """
    raise NotImplementedError


def surface_area(grid: object) -> float:
    """Compute the surface area of a reflector mesh.

    Parameters
    ----------
    grid : bempp.api.Grid
        Reflector surface mesh.

    Returns
    -------
    float
        Surface area in m^2.
    """
    raise NotImplementedError
