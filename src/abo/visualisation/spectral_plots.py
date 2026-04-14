"""Spectral directional pattern plots.

Visualises target strength as a function of frequency and incidence angle,
using polar and heatmap representations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def plot_ts_heatmap(
    ts_matrix: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    angles: NDArray[np.float64],
) -> Figure:
    """Plot target strength as a frequency-angle heatmap.

    Parameters
    ----------
    ts_matrix : NDArray[np.float64]
        Target strength in dB, shape (n_freq, n_angle).
    frequencies : NDArray[np.float64]
        Frequency array in Hz.
    angles : NDArray[np.float64]
        Angle array in radians.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    raise NotImplementedError


def plot_ts_polar(
    ts_at_freq: NDArray[np.float64],
    angles: NDArray[np.float64],
    frequency: float,
) -> Figure:
    """Plot target strength as a polar directivity pattern at one frequency.

    Parameters
    ----------
    ts_at_freq : NDArray[np.float64]
        TS values in dB at the given frequency.
    angles : NDArray[np.float64]
        Angle array in radians.
    frequency : float
        Frequency label in Hz.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    raise NotImplementedError
