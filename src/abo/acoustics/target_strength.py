"""Target strength computation from BEM solutions.

TS(f, theta) = 10 * log10( |p_scat(r, theta)|^2 * r^2 / |p_inc|^2 )

evaluated in the far field, in the backscatter (monostatic) direction.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def monostatic_target_strength(
    grid: object,
    frequencies: NDArray[np.float64],
    incidence_angles: NDArray[np.float64],
    speed_of_sound: float = 343.0,
    far_field_range: float = 1.0,
) -> NDArray[np.float64]:
    """Compute monostatic target strength over frequencies and angles.

    Parameters
    ----------
    grid : bempp.api.Grid
        Reflector surface mesh.
    frequencies : NDArray[np.float64]
        Array of frequencies in Hz.
    incidence_angles : NDArray[np.float64]
        Array of incidence angles in radians (0 = on-axis).
    speed_of_sound : float
        Speed of sound in m/s.
    far_field_range : float
        Range for far-field evaluation in metres.

    Returns
    -------
    NDArray[np.float64]
        Target strength in dB, shape (len(frequencies), len(incidence_angles)).
    """
    raise NotImplementedError
