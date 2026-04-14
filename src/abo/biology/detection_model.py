"""Atmospheric attenuation and bat detection threshold models.

Used for computing catchment volume (CV): the spatial volume within
which a reflector echo exceeds the bat's detection threshold.
"""

from __future__ import annotations


def atmospheric_attenuation(
    frequency: float,
    distance: float,
    temperature: float = 20.0,
    humidity: float = 50.0,
) -> float:
    """Compute atmospheric attenuation in dB.

    Parameters
    ----------
    frequency : float
        Frequency in Hz.
    distance : float
        Propagation distance in metres.
    temperature : float
        Temperature in degrees Celsius.
    humidity : float
        Relative humidity in percent.

    Returns
    -------
    float
        Attenuation in dB (positive value).
    """
    raise NotImplementedError


def detection_threshold_db() -> float:
    """Return the estimated bat detection threshold.

    Returns
    -------
    float
        Detection threshold in dB re 1 Pa (typically -60 to -70 dB).
    """
    return -65.0
