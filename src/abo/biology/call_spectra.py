"""Glossophagine bat FM sweep call spectrum models.

Models broadband frequency-modulated echolocation calls spanning
approximately 55-140 kHz, based on published recordings of
Glossophaga soricina and Leptonycteris yerbabuenae.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def glossophaga_call_spectrum(
    frequencies: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Power spectrum weighting for Glossophaga soricina FM call.

    Parameters
    ----------
    frequencies : NDArray[np.float64]
        Frequency array in Hz.

    Returns
    -------
    NDArray[np.float64]
        Normalised power spectrum W(f) at each frequency.
    """
    raise NotImplementedError


def leptonycteris_call_spectrum(
    frequencies: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Power spectrum weighting for Leptonycteris yerbabuenae FM call.

    Parameters
    ----------
    frequencies : NDArray[np.float64]
        Frequency array in Hz.

    Returns
    -------
    NDArray[np.float64]
        Normalised power spectrum W(f) at each frequency.
    """
    raise NotImplementedError
