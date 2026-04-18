"""Integration test: validate BEM against Simon et al. (2020) FEM results.

Simon et al. (2020, PNAS 117: 1367-1374) computed FEM impulse responses for
spherical cap reflectors with r in {35, 50, 70} mm and d in {25, 30, 49} mm.

Proper validation requires the raw FEM target strength values from the
published work, which are not included in the paper or its supplementary
materials. These tests are skipped until the original data can be obtained
(planned: direct request to the corresponding author, Marc Holderied,
University of Bristol).

The helper functions `_cap_from_radius_depth` and `_on_axis_ts` are
retained for future use when the reference data is available.
"""

from __future__ import annotations

import numpy as np
import pytest

from abo.acoustics.bem_solver import scattered_field, solve_scattering
from abo.geometry.meshing import mesh_spherical_cap
from abo.geometry.reflectors import SphericalCap

SPEED_OF_SOUND = 343.0
FAR_FIELD_RANGE = 1.0


def _cap_from_radius_depth(radius_mm: float, depth_mm: float) -> SphericalCap:
    """Create a SphericalCap from radius and depth in mm."""
    r = radius_mm * 1e-3
    d = depth_mm * 1e-3
    theta = np.arccos(1.0 - d / r)
    return SphericalCap(radius=r, half_angle=theta)


def _on_axis_ts(cap: SphericalCap, frequency: float) -> float:
    """Compute on-axis backscatter TS for a spherical cap."""
    grid = mesh_spherical_cap(cap, frequency=frequency)
    direction = np.array([0.0, 0.0, -1.0])
    total = solve_scattering(grid, frequency, direction)
    backscatter_pt = np.array([[0.0], [0.0], [FAR_FIELD_RANGE]])
    p_scat = scattered_field(
        grid, total, backscatter_pt, frequency, direction,
    )
    return float(
        10.0 * np.log10(np.abs(p_scat[0]) ** 2 * FAR_FIELD_RANGE**2),
    )


@pytest.mark.integration
@pytest.mark.skip(
    reason="Awaiting raw FEM data from Simon et al. (2020) authors.",
)
class TestSimon2020Validation:
    def test_spectral_notch_positions(self) -> None:
        """Spectral notch positions should match published FEM (within 1 kHz)."""

    def test_target_strength_agreement(self) -> None:
        """BEM TS should agree with published FEM TS within 2 dB across
        frequencies and angles for all nine configurations."""

    def test_impulse_response_peaks(self) -> None:
        """IR peak positions should match published FEM to within
        one sample at 500 kHz sampling."""
