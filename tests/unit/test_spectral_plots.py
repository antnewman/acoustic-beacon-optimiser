"""Unit tests for spectral directional pattern plots."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for tests

import numpy as np
import pytest
from matplotlib.figure import Figure

from abo.visualisation.spectral_plots import (
    plot_on_axis_vs_frequency,
    plot_ts_heatmap,
    plot_ts_polar,
)


@pytest.fixture
def ts_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frequencies = np.linspace(30_000.0, 120_000.0, 40)
    angles = np.linspace(0.0, np.pi / 2, 30)
    ts = np.random.default_rng(seed=0).uniform(
        -40.0, -10.0, size=(len(frequencies), len(angles)),
    )
    return ts, frequencies, angles


class TestHeatmap:
    def test_returns_figure(
        self, ts_data: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        ts, f, theta = ts_data
        fig = plot_ts_heatmap(ts, f, theta)
        assert isinstance(fig, Figure)
        assert len(fig.axes) >= 1

    def test_title_and_clim(
        self, ts_data: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        ts, f, theta = ts_data
        fig = plot_ts_heatmap(
            ts, f, theta, title="Test heatmap", vmin=-40.0, vmax=0.0,
        )
        ax = fig.axes[0]
        assert ax.get_title() == "Test heatmap"


class TestPolar:
    def test_returns_figure(
        self, ts_data: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        ts, f, theta = ts_data
        fig = plot_ts_polar(ts[0, :], theta, f[0])
        assert isinstance(fig, Figure)

    def test_default_title_contains_khz(
        self, ts_data: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        ts, f, theta = ts_data
        fig = plot_ts_polar(ts[0, :], theta, f[0])
        assert "kHz" in fig.axes[0].get_title()


class TestOnAxis:
    def test_returns_figure(
        self, ts_data: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        ts, f, theta = ts_data
        fig = plot_on_axis_vs_frequency(ts, f, theta)
        assert isinstance(fig, Figure)

    def test_xaxis_is_frequency(
        self, ts_data: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> None:
        ts, f, theta = ts_data
        fig = plot_on_axis_vs_frequency(ts, f, theta)
        assert "Frequency" in fig.axes[0].get_xlabel()
