"""Unit tests for parametric reflector definitions."""

from __future__ import annotations

import numpy as np
import pytest

from abo.geometry.reflectors import SphericalCap


class TestSphericalCap:
    def test_depth(self) -> None:
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 4)
        expected = 0.050 * (1.0 - np.cos(np.pi / 4))
        assert cap.depth == pytest.approx(expected, rel=1e-10)

    def test_aperture(self) -> None:
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 4)
        expected = 0.050 * np.sin(np.pi / 4)
        assert cap.aperture == pytest.approx(expected, rel=1e-10)

    def test_surface_area(self) -> None:
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 4)
        expected = 2.0 * np.pi * 0.050**2 * (1.0 - np.cos(np.pi / 4))
        assert cap.surface_area == pytest.approx(expected, rel=1e-10)

    def test_hemisphere(self) -> None:
        """A hemisphere has depth = radius and area = 2*pi*r^2."""
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 2)
        assert cap.depth == pytest.approx(0.050, rel=1e-10)
        assert cap.surface_area == pytest.approx(2.0 * np.pi * 0.050**2, rel=1e-10)

    def test_frozen(self) -> None:
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 4)
        with pytest.raises(AttributeError):
            cap.radius = 0.1  # type: ignore[misc]
