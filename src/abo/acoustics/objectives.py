"""Biologically meaningful objective functions for reflector optimisation.

- Integrated Conspicuousness (IC): frequency- and angle-weighted backscatter
- Surface Area (SA): structural/photosynthetic cost proxy
- Catchment Volume (CV): spatial volume exceeding detection threshold
  (placeholder; requires atmospheric attenuation and source-level model)

All IC/SA computations take either a pre-computed target strength matrix
(for offline post-processing) or an optimised pipeline that runs the BEM
solver internally (for use as optimisation objectives).
"""

from __future__ import annotations

import numpy as np
from bempp_cl.api import Grid
from numpy.typing import NDArray

from abo.acoustics.target_strength import monostatic_target_strength


def integrated_conspicuousness(
    ts_matrix: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    angles: NDArray[np.float64],
    call_spectrum: NDArray[np.float64],
) -> float:
    """Compute integrated conspicuousness.

    IC = int_f int_theta W(f) * (|p_scat|/|p_inc|)^2 * sin(theta)
             d(theta) df

    where (|p_scat|/|p_inc|)^2 = 10^(TS/10) converts the TS (dB) matrix
    back to a linear scattering intensity. The result is then converted
    back to dB.

    Parameters
    ----------
    ts_matrix : NDArray[np.float64]
        Target strength in dB, shape (n_freq, n_angle).
    frequencies : NDArray[np.float64]
        Frequency array in Hz.
    angles : NDArray[np.float64]
        Incidence angle array in radians.
    call_spectrum : NDArray[np.float64]
        Weighting function W(f) at each frequency (same length as
        `frequencies`). Expected to be normalised so that integral over
        frequency is 1.

    Returns
    -------
    float
        Integrated conspicuousness in dB.
    """
    ts = np.asarray(ts_matrix, dtype=np.float64)
    f = np.asarray(frequencies, dtype=np.float64)
    theta = np.asarray(angles, dtype=np.float64)
    w = np.asarray(call_spectrum, dtype=np.float64)

    intensity = 10.0 ** (ts / 10.0)
    # Weight angular integral by sin(theta) (solid-angle element)
    angular_weight = np.sin(theta)
    # I(f) = int_theta intensity * sin(theta) d(theta)
    i_of_f = np.trapezoid(intensity * angular_weight[np.newaxis, :], theta)
    # Weighted frequency integral
    ic_linear = np.trapezoid(w * i_of_f, f)

    if ic_linear <= 0.0:
        return -np.inf
    return float(10.0 * np.log10(ic_linear))


def surface_area(grid: Grid) -> float:
    """Compute the total surface area of a BEM mesh.

    Sums the triangle areas from the vertex coordinates and element
    connectivity of a bempp-cl Grid.

    Parameters
    ----------
    grid : Grid
        Reflector surface mesh.

    Returns
    -------
    float
        Total surface area in m^2.
    """
    vertices = np.asarray(grid.vertices, dtype=np.float64)
    elements = np.asarray(grid.elements, dtype=np.int64)

    v0 = vertices[:, elements[0, :]]
    v1 = vertices[:, elements[1, :]]
    v2 = vertices[:, elements[2, :]]

    # Cross product magnitudes give twice the triangle areas
    edge1 = v1 - v0
    edge2 = v2 - v0
    cross = np.cross(edge1, edge2, axis=0)
    areas = 0.5 * np.linalg.norm(cross, axis=0)

    return float(np.sum(areas))


def compute_ic_for_grid(
    grid: Grid,
    frequencies: NDArray[np.float64],
    angles: NDArray[np.float64],
    call_spectrum: NDArray[np.float64],
    speed_of_sound: float = 343.0,
    far_field_range: float = 1.0,
    gmres_tol: float = 1e-5,
) -> float:
    """Compute integrated conspicuousness directly from a mesh.

    Runs the BEM solver across the given frequency and angle grid,
    then integrates against the call spectrum. Convenience wrapper for
    use as an optimisation objective.

    Parameters
    ----------
    grid : Grid
        Reflector surface mesh.
    frequencies : NDArray[np.float64]
        Frequencies in Hz.
    angles : NDArray[np.float64]
        Incidence angles in radians.
    call_spectrum : NDArray[np.float64]
        Weighting function matched to frequencies.
    speed_of_sound : float
        Speed of sound in m/s.
    far_field_range : float
        Far-field evaluation range in metres.
    gmres_tol : float
        GMRES tolerance.

    Returns
    -------
    float
        Integrated conspicuousness in dB.
    """
    ts_matrix = monostatic_target_strength(
        grid,
        frequencies,
        angles,
        speed_of_sound=speed_of_sound,
        far_field_range=far_field_range,
        gmres_tol=gmres_tol,
    )
    return integrated_conspicuousness(
        ts_matrix, frequencies, angles, call_spectrum,
    )
