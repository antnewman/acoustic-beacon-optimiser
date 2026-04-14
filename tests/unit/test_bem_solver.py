"""Unit tests for BEM solver."""

from __future__ import annotations

import sys

import numpy as np
import pytest

sys.modules["bempp"] = __import__("bempp_cl")
sys.modules["bempp.api"] = __import__("bempp_cl.api", fromlist=["api"])
import bempp.api  # noqa: E402

from abo.acoustics.bem_solver import scattered_field, solve_scattering  # noqa: E402


@pytest.fixture
def sphere_grid() -> bempp.api.Grid:
    """Coarse regular sphere mesh for fast tests."""
    return bempp.api.shapes.regular_sphere(3)


class TestSolveScattering:
    def test_returns_grid_function(
        self, sphere_grid: bempp.api.Grid,
    ) -> None:
        direction = np.array([0.0, 0.0, -1.0])
        result = solve_scattering(
            sphere_grid, frequency=1000.0, incidence_direction=direction,
        )
        assert isinstance(result, bempp.api.GridFunction)

    def test_coefficients_not_zero(
        self, sphere_grid: bempp.api.Grid,
    ) -> None:
        direction = np.array([0.0, 0.0, -1.0])
        result = solve_scattering(
            sphere_grid, frequency=1000.0, incidence_direction=direction,
        )
        assert np.any(np.abs(result.coefficients) > 0)


class TestScatteredField:
    def test_returns_complex_array(
        self, sphere_grid: bempp.api.Grid,
    ) -> None:
        direction = np.array([0.0, 0.0, -1.0])
        total = solve_scattering(
            sphere_grid, frequency=1000.0, incidence_direction=direction,
        )
        pts = np.array([[0.0], [0.0], [2.0]])
        p_scat = scattered_field(
            sphere_grid, total, pts, 1000.0, direction,
        )
        assert p_scat.shape == (1,)
        assert np.issubdtype(p_scat.dtype, np.complexfloating)

    def test_scattered_field_nonzero(
        self, sphere_grid: bempp.api.Grid,
    ) -> None:
        direction = np.array([0.0, 0.0, -1.0])
        total = solve_scattering(
            sphere_grid, frequency=1000.0, incidence_direction=direction,
        )
        pts = np.array([[0.0], [0.0], [2.0]])
        p_scat = scattered_field(
            sphere_grid, total, pts, 1000.0, direction,
        )
        assert np.abs(p_scat[0]) > 0

    def test_multiple_points(
        self, sphere_grid: bempp.api.Grid,
    ) -> None:
        direction = np.array([0.0, 0.0, -1.0])
        total = solve_scattering(
            sphere_grid, frequency=1000.0, incidence_direction=direction,
        )
        pts = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 2.0], [2.0, 2.0, 2.0]])
        p_scat = scattered_field(
            sphere_grid, total, pts, 1000.0, direction,
        )
        assert p_scat.shape == (3,)
