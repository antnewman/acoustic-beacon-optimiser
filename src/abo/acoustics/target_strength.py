"""Target strength computation from BEM solutions.

TS(f, theta) = 10 * log10( |p_scat(r, theta)|^2 * r^2 / |p_inc|^2 )

evaluated in the far field, in the backscatter (monostatic) direction.
For a unit-amplitude incident plane wave, |p_inc|^2 = 1.
"""

from __future__ import annotations

import numpy as np
from bempp_cl.api import Grid as BemppGrid
from numpy.typing import NDArray

from abo.acoustics.bem_solver import scattered_field, solve_scattering


def monostatic_target_strength(
    grid: BemppGrid,
    frequencies: NDArray[np.float64],
    incidence_angles: NDArray[np.float64],
    speed_of_sound: float = 343.0,
    far_field_range: float = 1.0,
    gmres_tol: float = 1e-5,
) -> NDArray[np.float64]:
    """Compute monostatic target strength over frequencies and angles.

    For each (frequency, angle) pair, solves the BEM scattering problem
    with the incident wave arriving from the given angle, then evaluates
    the backscattered pressure at the far-field range.

    The incidence direction lies in the xz-plane: the wave travels
    toward the reflector, and theta = 0 corresponds to on-axis
    incidence along -z (head-on to a reflector with pole at +z).

    Parameters
    ----------
    grid : BemppGrid
        Reflector surface mesh.
    frequencies : NDArray[np.float64]
        Array of frequencies in Hz.
    incidence_angles : NDArray[np.float64]
        Array of incidence angles in radians (0 = on-axis).
    speed_of_sound : float
        Speed of sound in m/s.
    far_field_range : float
        Range for far-field evaluation in metres.
    gmres_tol : float
        GMRES convergence tolerance.

    Returns
    -------
    NDArray[np.float64]
        Target strength in dB, shape (len(frequencies), len(angles)).
    """
    freqs = np.asarray(frequencies, dtype=np.float64).ravel()
    angles = np.asarray(incidence_angles, dtype=np.float64).ravel()
    ts = np.empty((len(freqs), len(angles)), dtype=np.float64)

    for i, f in enumerate(freqs):
        for j, theta in enumerate(angles):
            # Incidence direction: wave travels toward reflector
            # theta=0 is along -z (on-axis for cap with pole at +z)
            d = np.array([
                -np.sin(theta), 0.0, -np.cos(theta),
            ])

            # Backscatter point: far-field range in the direction
            # the wave came from (opposite to incidence)
            backscatter_pt = -far_field_range * d
            eval_pts = backscatter_pt.reshape(3, 1)

            total_field = solve_scattering(
                grid, f, d, speed_of_sound, gmres_tol,
            )
            p_scat = scattered_field(
                grid, total_field, eval_pts, f, d, speed_of_sound,
            )

            # TS = 10*log10(|p_scat|^2 * r^2) for unit incident amplitude
            ts[i, j] = 10.0 * np.log10(
                np.abs(p_scat[0]) ** 2 * far_field_range**2,
            )

    return ts
