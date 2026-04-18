"""Tests for optimisation runners.

Exercises the encoder factory and the structure of the generated cost
functions. Does not run BEM; that would require many seconds per call
and is exercised end-to-end by the CLI in integration runs.
"""

from __future__ import annotations

import numpy as np
import pytest

from abo.biology.call_spectra import flat_band_spectrum
from abo.optimisation.encoding import (
    ChebyshevProfileEncoding,
    EllipsoidalCapEncoding,
    SphericalCapEncoding,
)
from abo.optimisation.runners import (
    EvaluationConfig,
    build_encoder,
    make_ic_cost,
    make_pareto_objectives,
)


@pytest.fixture
def eval_config() -> EvaluationConfig:
    freqs = np.linspace(40_000.0, 100_000.0, 5)
    angles = np.linspace(0.0, np.pi / 3, 5)
    spectrum = flat_band_spectrum(freqs)
    return EvaluationConfig(
        frequencies=freqs, angles=angles, call_spectrum=spectrum,
    )


class TestBuildEncoder:
    def test_spherical(self) -> None:
        enc = build_encoder("spherical-cap")
        assert isinstance(enc, SphericalCapEncoding)

    def test_ellipsoidal(self) -> None:
        enc = build_encoder("ellipsoidal")
        assert isinstance(enc, EllipsoidalCapEncoding)

    def test_chebyshev(self) -> None:
        enc = build_encoder(
            "chebyshev", n_coefficients=3, aperture_radius=0.040,
        )
        assert isinstance(enc, ChebyshevProfileEncoding)
        assert enc.n_coefficients == 3
        assert enc.aperture_radius == 0.040

    def test_unknown_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown family"):
            build_encoder("not-a-family")


class TestMakeIcCost:
    def test_returns_callable(self, eval_config: EvaluationConfig) -> None:
        enc = SphericalCapEncoding()
        cost = make_ic_cost(enc, eval_config, area_max=0.010)
        assert callable(cost)


class TestMakeParetoObjectives:
    def test_returns_two_callables(
        self, eval_config: EvaluationConfig,
    ) -> None:
        enc = SphericalCapEncoding()
        objs = make_pareto_objectives(enc, eval_config)
        assert len(objs) == 2
        for obj in objs:
            assert callable(obj)
