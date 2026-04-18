"""Glossophagine bat FM sweep call spectrum models.

Models broadband frequency-modulated echolocation calls of nectar-feeding
bats based on published recordings:

- Glossophaga soricina: two-harmonic FM sweep, fundamental 56-70 kHz,
  second harmonic 112-140 kHz, total bandwidth roughly 55-140 kHz.
  (Simon et al., 2006, J Exp Biol 209: 3599-3609.)
- Leptonycteris yerbabuenae: FM sweep with bandwidth around 50-115 kHz,
  dominant energy 70-90 kHz. Call parameters adapt across search,
  approach, and terminal phases.
  (Gonzalez-Terrazas et al., 2016, PLoS ONE 11: e0163492.)

The spectra returned by the functions below are power spectral density
profiles W(f) used as weighting functions in the Integrated
Conspicuousness (IC) objective. They are normalised so that the integral
over the relevant band equals 1 (unit total energy).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _gaussian_bump(
    frequencies: NDArray[np.float64],
    centre: float,
    bandwidth: float,
) -> NDArray[np.float64]:
    """Un-normalised Gaussian centred at `centre` with full-width half-max
    approximately equal to `bandwidth`."""
    sigma = bandwidth / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return np.exp(-0.5 * ((frequencies - centre) / sigma) ** 2)


def _normalise(
    frequencies: NDArray[np.float64], spectrum: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Scale the spectrum so that the trapezoidal integral equals 1."""
    total = np.trapezoid(spectrum, frequencies)
    if total <= 0.0:
        return spectrum
    return spectrum / total


def glossophaga_call_spectrum(
    frequencies: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Power spectrum weighting for Glossophaga soricina FM call.

    Two-harmonic model based on Simon et al. (2006). The fundamental is
    centred at 63 kHz (average of 56-70 kHz band), and the second
    harmonic is centred at 126 kHz. Second-harmonic energy is about
    half that of the fundamental in typical recordings.

    Parameters
    ----------
    frequencies : NDArray[np.float64]
        Frequency array in Hz.

    Returns
    -------
    NDArray[np.float64]
        Normalised power spectrum W(f) integrating to 1 over the input
        range.
    """
    f = np.asarray(frequencies, dtype=np.float64)
    fundamental = _gaussian_bump(f, centre=63_000.0, bandwidth=14_000.0)
    second = 0.5 * _gaussian_bump(f, centre=126_000.0, bandwidth=28_000.0)
    raw = fundamental + second
    return _normalise(f, raw)


def leptonycteris_call_spectrum(
    frequencies: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Power spectrum weighting for Leptonycteris yerbabuenae FM call.

    Single dominant sweep centred near 80 kHz with about 45 kHz of
    bandwidth (covers 50-115 kHz at approximately half-maximum), based
    on Gonzalez-Terrazas et al. (2016).

    Parameters
    ----------
    frequencies : NDArray[np.float64]
        Frequency array in Hz.

    Returns
    -------
    NDArray[np.float64]
        Normalised power spectrum W(f).
    """
    f = np.asarray(frequencies, dtype=np.float64)
    raw = _gaussian_bump(f, centre=80_000.0, bandwidth=45_000.0)
    return _normalise(f, raw)


def flat_band_spectrum(
    frequencies: NDArray[np.float64],
    f_min: float = 30_000.0,
    f_max: float = 120_000.0,
) -> NDArray[np.float64]:
    """Uniform (flat) spectrum over [f_min, f_max], zero elsewhere.

    Useful as a neutral weighting that treats all frequencies equally,
    for sensitivity analyses and sanity checks.

    Parameters
    ----------
    frequencies : NDArray[np.float64]
        Frequency array in Hz.
    f_min : float
        Lower band edge in Hz.
    f_max : float
        Upper band edge in Hz.

    Returns
    -------
    NDArray[np.float64]
        Normalised flat spectrum W(f).
    """
    f = np.asarray(frequencies, dtype=np.float64)
    raw = np.where((f >= f_min) & (f <= f_max), 1.0, 0.0)
    return _normalise(f, raw)
