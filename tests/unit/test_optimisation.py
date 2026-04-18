"""Unit tests for optimisation wrappers.

These tests use fast synthetic objectives (sphere, Rastrigin, dual
quadratics) rather than BEM. Full BEM-backed optimisation runs are
performed by the CLI runners and are not part of the unit test suite.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

from abo.optimisation.constraints import concavity_penalty
from abo.optimisation.encoding import (
    ChebyshevProfileEncoding,
    EllipsoidalCapEncoding,
    SphericalCapEncoding,
)
from abo.optimisation.multi_objective import pareto_frontier
from abo.optimisation.single_objective import optimise_cmaes


def _sphere(x: NDArray[np.float64]) -> float:
    """f(x) = sum x_i^2. Minimum at origin, f=0."""
    return float(np.sum(x**2))


class TestOptimiseCmaes:
    def test_converges_on_sphere(self) -> None:
        result = optimise_cmaes(
            objective=_sphere,
            x0=np.array([2.0, -3.0]),
            sigma0=1.0,
            max_evaluations=500,
            seed=42,
        )
        assert result.f_best < 1e-6
        assert np.linalg.norm(result.x_best) < 1e-3
        assert result.n_evals > 0

    def test_respects_bounds(self) -> None:
        # Minimum at origin, but bounds force x[0] >= 1
        result = optimise_cmaes(
            objective=_sphere,
            x0=np.array([2.0, 2.0]),
            sigma0=1.0,
            bounds=[(1.0, 5.0), (-5.0, 5.0)],
            max_evaluations=500,
            seed=42,
        )
        assert result.x_best[0] >= 1.0 - 1e-6

    def test_max_evaluations_budget(self) -> None:
        result = optimise_cmaes(
            objective=_sphere,
            x0=np.array([2.0, -3.0]),
            sigma0=1.0,
            max_evaluations=50,
            seed=0,
        )
        assert result.n_evals <= 60  # may slightly overshoot within a gen

    def test_history_monotone_non_increasing(self) -> None:
        result = optimise_cmaes(
            objective=_sphere,
            x0=np.array([2.0, -3.0]),
            sigma0=1.0,
            max_evaluations=200,
            seed=42,
        )
        hist = np.array(result.history)
        assert np.all(np.diff(hist) <= 1e-12)


class TestParetoFrontier:
    def test_two_objective_frontier(self) -> None:
        """Two competing quadratics with minima at 0 and at 1."""

        def f1(x: NDArray[np.float64]) -> float:
            return float((x[0] - 0.0) ** 2)

        def f2(x: NDArray[np.float64]) -> float:
            return float((x[0] - 1.0) ** 2)

        result = pareto_frontier(
            objectives=[f1, f2],
            bounds=[(-2.0, 2.0)],
            population_size=20,
            n_generations=30,
            seed=42,
        )
        # Pareto-optimal set must lie within [0, 1] on x[0]
        assert np.all(result.pareto_set[:, 0] >= -1e-2)
        assert np.all(result.pareto_set[:, 0] <= 1.0 + 1e-2)
        assert result.pareto_front.shape[1] == 2

    def test_n_evals_counted(self) -> None:
        def f1(x: NDArray[np.float64]) -> float:
            return float(x[0] ** 2)

        def f2(x: NDArray[np.float64]) -> float:
            return float((x[0] - 1.0) ** 2)

        result = pareto_frontier(
            objectives=[f1, f2],
            bounds=[(-2.0, 2.0)],
            population_size=10,
            n_generations=5,
            seed=0,
        )
        assert result.n_evals > 0


class TestEncodings:
    def test_spherical_cap_roundtrip(self) -> None:
        enc = SphericalCapEncoding()
        assert enc.dim == 2
        x = enc.x0
        cap = enc.decode(x)
        assert cap.radius == pytest.approx(x[0])
        assert cap.half_angle == pytest.approx(x[1])

    def test_ellipsoidal_cap_roundtrip(self) -> None:
        enc = EllipsoidalCapEncoding()
        assert enc.dim == 3
        cap = enc.decode(enc.x0)
        assert cap.semi_major > 0.0
        assert cap.semi_minor > 0.0

    def test_chebyshev_default_is_concave(self) -> None:
        enc = ChebyshevProfileEncoding(
            n_coefficients=4, aperture_radius=0.050,
        )
        assert enc.dim == 4
        profile = enc.decode(enc.x0)
        assert profile.is_concave() is True

    def test_bounds_consistent_with_dim(self) -> None:
        for enc in [
            SphericalCapEncoding(),
            EllipsoidalCapEncoding(),
            ChebyshevProfileEncoding(
                n_coefficients=5, aperture_radius=0.050,
            ),
        ]:
            bounds = enc.default_bounds
            assert len(bounds) == enc.dim
            assert all(lo < hi for lo, hi in bounds)


class TestConcavityPenalty:
    def test_zero_for_decreasing_profile(self) -> None:
        coeffs = np.array([0.020, -0.010])  # z(0)=0.030, z(a)=0.010
        p = concavity_penalty(coeffs, aperture_radius=0.050)
        assert p == pytest.approx(0.0, abs=1e-12)

    def test_positive_for_increasing_profile(self) -> None:
        coeffs = np.array([0.010, 0.010])  # z(0)=0, z(a)=0.020
        p = concavity_penalty(coeffs, aperture_radius=0.050)
        assert p > 0.0
