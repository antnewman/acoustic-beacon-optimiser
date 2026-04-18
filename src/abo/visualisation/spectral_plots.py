"""Spectral directional pattern plots.

Visualises target strength as a function of frequency and incidence
angle, using both heatmap (frequency-angle) and polar directivity
representations.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray


def plot_ts_heatmap(
    ts_matrix: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    angles: NDArray[np.float64],
    title: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Figure:
    """Plot target strength as a frequency-angle heatmap.

    Parameters
    ----------
    ts_matrix : NDArray[np.float64]
        Target strength in dB, shape (n_freq, n_angle).
    frequencies : NDArray[np.float64]
        Frequency array in Hz.
    angles : NDArray[np.float64]
        Angle array in radians (0 = on-axis).
    title : str, optional
        Plot title.
    vmin, vmax : float, optional
        Colour-scale limits in dB.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    ts = np.asarray(ts_matrix, dtype=np.float64)
    f_khz = np.asarray(frequencies, dtype=np.float64) / 1000.0
    theta_deg = np.degrees(np.asarray(angles, dtype=np.float64))

    fig, ax = plt.subplots(figsize=(8, 5))
    extent = (
        float(theta_deg.min()),
        float(theta_deg.max()),
        float(f_khz.min()),
        float(f_khz.max()),
    )
    im = ax.imshow(
        ts,
        aspect="auto",
        origin="lower",
        extent=extent,
        vmin=vmin,
        vmax=vmax,
        cmap="viridis",
    )
    ax.set_xlabel("Incidence angle (deg)")
    ax.set_ylabel("Frequency (kHz)")
    if title:
        ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Target strength (dB)")
    fig.tight_layout()
    return fig


def plot_ts_polar(
    ts_at_freq: NDArray[np.float64],
    angles: NDArray[np.float64],
    frequency: float,
    title: str | None = None,
) -> Figure:
    """Plot target strength as a polar directivity pattern at one frequency.

    The pattern is mirrored around the axis (0 degrees) to show a full
    symmetric lobe, a common convention in transducer directivity plots.

    Parameters
    ----------
    ts_at_freq : NDArray[np.float64]
        TS values in dB at the given frequency.
    angles : NDArray[np.float64]
        Angle array in radians (0 = on-axis).
    frequency : float
        Frequency in Hz (for the label).
    title : str, optional
        Plot title override.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    ts = np.asarray(ts_at_freq, dtype=np.float64)
    theta = np.asarray(angles, dtype=np.float64)

    # Mirror for full symmetric polar plot
    theta_full = np.concatenate([-theta[::-1], theta])
    ts_full = np.concatenate([ts[::-1], ts])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="polar")
    ax.plot(theta_full, ts_full, linewidth=1.5)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(-90, 91, 30))
    ax.set_rlabel_position(135.0)
    ax.grid(visible=True, linestyle=":", linewidth=0.5)
    default_title = f"Directivity at {frequency / 1000:.1f} kHz"
    ax.set_title(title if title else default_title)
    fig.tight_layout()
    return fig


def plot_on_axis_vs_frequency(
    ts_matrix: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    angles: NDArray[np.float64],
    title: str | None = None,
) -> Figure:
    """Plot on-axis TS as a function of frequency.

    Convenience plot for comparing the spectral fingerprint of different
    reflector geometries.

    Parameters
    ----------
    ts_matrix : NDArray[np.float64]
        TS matrix in dB, shape (n_freq, n_angle).
    frequencies : NDArray[np.float64]
        Frequency array in Hz.
    angles : NDArray[np.float64]
        Angle array; the on-axis index is found as argmin(abs(angles)).
    title : str, optional
        Plot title.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    ts = np.asarray(ts_matrix, dtype=np.float64)
    f_khz = np.asarray(frequencies, dtype=np.float64) / 1000.0
    theta = np.asarray(angles, dtype=np.float64)
    i_axis = int(np.argmin(np.abs(theta)))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(f_khz, ts[:, i_axis], linewidth=1.5)
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("On-axis TS (dB)")
    ax.grid(visible=True, linestyle=":", linewidth=0.5)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig
