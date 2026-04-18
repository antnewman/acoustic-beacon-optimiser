"""Unit tests for objective functions."""

from __future__ import annotations

import numpy as np
import pytest

from abo.acoustics.objectives import (
    integrated_conspicuousness,
    surface_area,
)
from abo.biology.call_spectra import flat_band_spectrum
from abo.geometry.meshing import mesh_spherical_cap
from abo.geometry.reflectors import SphericalCap


class TestIntegratedConspicuousness:
    def test_constant_ts_reproduces_value(self) -> None:
        """With constant TS and a flat spectrum, IC should equal that TS."""
        freqs = np.linspace(30_000.0, 120_000.0, 50)
        angles = np.linspace(0.0, np.pi / 2, 30)
        # Constant TS = -20 dB everywhere
        ts = -20.0 * np.ones((len(freqs), len(angles)))
        spectrum = flat_band_spectrum(freqs, f_min=30_000.0, f_max=120_000.0)

        ic = integrated_conspicuousness(ts, freqs, angles, spectrum)

        # Expected value: IC integrates intensity * sin(theta) over angles
        # and then over frequency. The angular integral of sin(theta)
        # from 0 to pi/2 equals 1, so IC = 10*log10(10^(-2) * 1) = -20 dB.
        assert ic == pytest.approx(-20.0, abs=0.5)

    def test_higher_ts_gives_higher_ic(self) -> None:
        freqs = np.linspace(30_000.0, 120_000.0, 50)
        angles = np.linspace(0.0, np.pi / 2, 30)
        spectrum = flat_band_spectrum(freqs)
        ts_low = -30.0 * np.ones((len(freqs), len(angles)))
        ts_high = -10.0 * np.ones((len(freqs), len(angles)))
        ic_low = integrated_conspicuousness(ts_low, freqs, angles, spectrum)
        ic_high = integrated_conspicuousness(
            ts_high, freqs, angles, spectrum,
        )
        assert ic_high > ic_low

    def test_shape_mismatch_raises(self) -> None:
        freqs = np.linspace(30_000.0, 120_000.0, 50)
        angles = np.linspace(0.0, np.pi / 2, 30)
        spectrum = flat_band_spectrum(freqs)
        ts = np.zeros((10, 10))  # wrong shape
        with pytest.raises((ValueError, IndexError)):
            integrated_conspicuousness(ts, freqs, angles, spectrum)


class TestSurfaceArea:
    def test_spherical_cap_area_matches_analytical(self) -> None:
        """Mesh-computed area should match the analytical formula."""
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 4)
        grid = mesh_spherical_cap(cap, frequency=30_000.0)
        area_mesh = surface_area(grid)
        area_analytical = cap.surface_area
        rel_err = abs(area_mesh - area_analytical) / area_analytical
        assert rel_err < 0.05

    def test_area_positive(self) -> None:
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 4)
        grid = mesh_spherical_cap(cap, frequency=30_000.0)
        assert surface_area(grid) > 0.0
