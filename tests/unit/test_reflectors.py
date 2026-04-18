"""Unit tests for parametric reflector definitions."""

from __future__ import annotations

import numpy as np
import pytest

from abo.geometry.reflectors import (
    ChebyshevProfile,
    EllipsoidalCap,
    SphericalCap,
)


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


class TestEllipsoidalCap:
    def test_depth(self) -> None:
        cap = EllipsoidalCap(
            semi_major=0.050, semi_minor=0.025, half_angle=np.pi / 3,
        )
        expected = 0.025 * (1.0 - np.cos(np.pi / 3))
        assert cap.depth == pytest.approx(expected, rel=1e-10)

    def test_aperture(self) -> None:
        cap = EllipsoidalCap(
            semi_major=0.050, semi_minor=0.025, half_angle=np.pi / 3,
        )
        expected = 0.050 * np.sin(np.pi / 3)
        assert cap.aperture == pytest.approx(expected, rel=1e-10)

    def test_reduces_to_sphere(self) -> None:
        """Ellipsoidal cap with semi_major == semi_minor equals a sphere."""
        r = 0.050
        theta = np.pi / 4
        ellipsoid = EllipsoidalCap(
            semi_major=r, semi_minor=r, half_angle=theta,
        )
        sphere = SphericalCap(radius=r, half_angle=theta)
        assert ellipsoid.depth == pytest.approx(sphere.depth, rel=1e-10)
        assert ellipsoid.aperture == pytest.approx(sphere.aperture, rel=1e-10)
        # Surface area from quadrature must match the analytical spherical
        # cap area.
        assert ellipsoid.surface_area == pytest.approx(
            sphere.surface_area, rel=1e-6,
        )

    def test_surface_area_positive(self) -> None:
        cap = EllipsoidalCap(
            semi_major=0.050, semi_minor=0.025, half_angle=np.pi / 3,
        )
        assert cap.surface_area > 0.0


class TestChebyshevProfile:
    def test_constant_profile_flat(self) -> None:
        """Coefficients [c0] produce a flat disc at z = c0."""
        profile = ChebyshevProfile(
            coefficients=np.array([0.010]),
            aperture_radius=0.050,
        )
        rho = np.linspace(0.0, 0.050, 5)
        z = profile.profile(rho)
        np.testing.assert_allclose(z, 0.010, atol=1e-12)

    def test_linear_profile(self) -> None:
        """Coefficients [c0, c1] give z(rho) = c0 + c1 * T_1(u) where
        u = 2*rho/a - 1 is linear in rho."""
        profile = ChebyshevProfile(
            coefficients=np.array([0.020, -0.010]),
            aperture_radius=0.050,
        )
        # At rho=0: u=-1, T_1(-1)=-1 -> z = 0.020 + (-0.010)*(-1) = 0.030
        # At rho=a: u=1,  T_1(1)=1   -> z = 0.020 + (-0.010)*1    = 0.010
        assert profile.profile(np.array([0.0]))[0] == pytest.approx(0.030)
        assert profile.profile(np.array([0.050]))[0] == pytest.approx(0.010)

    def test_depth_property(self) -> None:
        profile = ChebyshevProfile(
            coefficients=np.array([0.020, -0.010]),
            aperture_radius=0.050,
        )
        assert profile.depth == pytest.approx(0.020, abs=1e-12)

    def test_is_concave_for_monotone_decreasing(self) -> None:
        """A linear profile with negative slope is concave."""
        profile = ChebyshevProfile(
            coefficients=np.array([0.020, -0.010]),
            aperture_radius=0.050,
        )
        assert profile.is_concave() is True

    def test_not_concave_for_monotone_increasing(self) -> None:
        profile = ChebyshevProfile(
            coefficients=np.array([0.020, 0.010]),
            aperture_radius=0.050,
        )
        assert profile.is_concave() is False

    def test_surface_area_flat_disc(self) -> None:
        """A flat profile should have area = pi * a^2."""
        profile = ChebyshevProfile(
            coefficients=np.array([0.010]),
            aperture_radius=0.050,
        )
        expected = np.pi * 0.050**2
        assert profile.surface_area == pytest.approx(expected, rel=1e-6)
