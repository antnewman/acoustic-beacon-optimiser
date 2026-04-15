"""Analytical Mie series solution for scattering from a rigid sphere.

Provides the exact backscattered target strength for validating the BEM
solver. The Mie series expansion uses spherical Bessel and Hankel functions
to compute the scattering coefficients for a sound-hard (Neumann) sphere.
"""

from __future__ import annotations

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def _sph_hankel1_deriv(n: int, x: float) -> complex:
    """Derivative of spherical Hankel function of the first kind."""
    return (
        spherical_jn(n, x, derivative=True)
        + 1j * spherical_yn(n, x, derivative=True)
    )


def mie_backscatter_ts(
    sphere_radius: float,
    frequency: float,
    speed_of_sound: float = 343.0,
    far_field_range: float = 1.0,
) -> float:
    """Compute the exact monostatic target strength of a rigid sphere.

    Uses the Mie series expansion for the far-field backscattering
    amplitude of a sound-hard sphere.

    Parameters
    ----------
    sphere_radius : float
        Radius of the sphere in metres.
    frequency : float
        Frequency in Hz.
    speed_of_sound : float
        Speed of sound in m/s.
    far_field_range : float
        Far-field evaluation range in metres.

    Returns
    -------
    float
        Target strength in dB.
    """
    k = 2.0 * np.pi * frequency / speed_of_sound
    ka = k * sphere_radius

    # Truncation: standard rule for Mie series convergence
    n_terms = int(ka + 10.0 * ka ** (1.0 / 3.0) + 10)

    # Scattering coefficients: a_n = -j_n'(ka) / h_n'^(1)(ka)
    f_back = 0.0 + 0j
    for n in range(n_terms + 1):
        jn_d = spherical_jn(n, ka, derivative=True)
        hn_d = _sph_hankel1_deriv(n, ka)
        a_n = -jn_d / hn_d
        f_back += (2 * n + 1) * (-1) ** n * a_n

    f_back /= 1j * k

    # Far-field scattered pressure
    p_scat_mag = np.abs(f_back) / far_field_range

    return float(10.0 * np.log10(p_scat_mag**2 * far_field_range**2))
