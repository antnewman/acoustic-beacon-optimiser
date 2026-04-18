"""Unit tests for glossophagine call spectrum models."""

from __future__ import annotations

import numpy as np
import pytest

from abo.biology.call_spectra import (
    flat_band_spectrum,
    glossophaga_call_spectrum,
    leptonycteris_call_spectrum,
)


@pytest.fixture
def frequencies() -> np.ndarray:
    """Dense frequency grid spanning the glossophagine call range."""
    return np.linspace(20_000.0, 150_000.0, 2000)


class TestGlossophaga:
    def test_normalised(self, frequencies: np.ndarray) -> None:
        w = glossophaga_call_spectrum(frequencies)
        integral = np.trapezoid(w, frequencies)
        assert integral == pytest.approx(1.0, abs=1e-6)

    def test_positive(self, frequencies: np.ndarray) -> None:
        w = glossophaga_call_spectrum(frequencies)
        assert np.all(w >= 0.0)

    def test_fundamental_peak_in_range(self, frequencies: np.ndarray) -> None:
        """Maximum energy should sit in the fundamental band (55-75 kHz)."""
        w = glossophaga_call_spectrum(frequencies)
        f_peak = frequencies[np.argmax(w)]
        assert 55_000.0 <= f_peak <= 75_000.0


class TestLeptonycteris:
    def test_normalised(self, frequencies: np.ndarray) -> None:
        w = leptonycteris_call_spectrum(frequencies)
        integral = np.trapezoid(w, frequencies)
        assert integral == pytest.approx(1.0, abs=1e-6)

    def test_peak_near_80khz(self, frequencies: np.ndarray) -> None:
        w = leptonycteris_call_spectrum(frequencies)
        f_peak = frequencies[np.argmax(w)]
        assert 70_000.0 <= f_peak <= 90_000.0


class TestFlatBand:
    def test_normalised(self, frequencies: np.ndarray) -> None:
        w = flat_band_spectrum(frequencies)
        integral = np.trapezoid(w, frequencies)
        assert integral == pytest.approx(1.0, abs=1e-6)

    def test_zero_outside_band(self, frequencies: np.ndarray) -> None:
        w = flat_band_spectrum(frequencies, f_min=40_000.0, f_max=100_000.0)
        below = w[frequencies < 40_000.0]
        above = w[frequencies > 100_000.0]
        assert np.all(below == 0.0)
        assert np.all(above == 0.0)

    def test_nonzero_inside_band(self, frequencies: np.ndarray) -> None:
        w = flat_band_spectrum(frequencies, f_min=40_000.0, f_max=100_000.0)
        inside = w[(frequencies >= 40_000.0) & (frequencies <= 100_000.0)]
        assert np.all(inside > 0.0)
