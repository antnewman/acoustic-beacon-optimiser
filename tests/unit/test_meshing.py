"""Unit tests for mesh generation."""

from __future__ import annotations

import numpy as np
import pytest

from abo.geometry.meshing import mesh_spherical_cap
from abo.geometry.reflectors import SphericalCap


@pytest.fixture
def cap_45deg() -> SphericalCap:
    """Spherical cap with 50 mm radius and 45 degree half-angle."""
    return SphericalCap(radius=0.050, half_angle=np.pi / 4)


class TestMeshSphericalCap:
    def test_returns_grid(self, cap_45deg: SphericalCap) -> None:
        grid = mesh_spherical_cap(cap_45deg, frequency=60_000.0)
        assert grid.number_of_elements > 0
        assert grid.number_of_vertices > 0

    def test_vertices_shape(self, cap_45deg: SphericalCap) -> None:
        grid = mesh_spherical_cap(cap_45deg, frequency=60_000.0)
        assert grid.vertices.shape[0] == 3
        assert grid.elements.shape[0] == 3

    def test_vertices_on_sphere(self, cap_45deg: SphericalCap) -> None:
        """All mesh vertices should lie on the sphere surface."""
        grid = mesh_spherical_cap(cap_45deg, frequency=60_000.0)
        r = cap_45deg.radius
        distances = np.linalg.norm(grid.vertices, axis=0)
        np.testing.assert_allclose(distances, r, atol=1e-6)

    def test_cap_above_cut_plane(self, cap_45deg: SphericalCap) -> None:
        """No vertex should be below the cut plane z = r*cos(theta)."""
        grid = mesh_spherical_cap(cap_45deg, frequency=60_000.0)
        cut_z = cap_45deg.radius * np.cos(cap_45deg.half_angle)
        z_coords = grid.vertices[2, :]
        assert np.all(z_coords >= cut_z - 1e-10)

    def test_pole_present(self, cap_45deg: SphericalCap) -> None:
        """The mesh should include a vertex near the pole (0, 0, r)."""
        grid = mesh_spherical_cap(cap_45deg, frequency=60_000.0)
        r = cap_45deg.radius
        z_max = grid.vertices[2, :].max()
        assert z_max == pytest.approx(r, abs=1e-6)

    def test_more_elements_at_higher_frequency(
        self, cap_45deg: SphericalCap,
    ) -> None:
        """Higher frequency requires finer mesh; more elements."""
        grid_lo = mesh_spherical_cap(cap_45deg, frequency=30_000.0)
        grid_hi = mesh_spherical_cap(cap_45deg, frequency=90_000.0)
        assert grid_hi.number_of_elements > grid_lo.number_of_elements

    def test_hemisphere(self) -> None:
        """A hemisphere (half_angle = pi/2) should have cut plane at z = 0."""
        cap = SphericalCap(radius=0.050, half_angle=np.pi / 2)
        grid = mesh_spherical_cap(cap, frequency=60_000.0)
        z_min = grid.vertices[2, :].min()
        assert z_min >= -1e-10
