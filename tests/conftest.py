"""Shared pytest fixtures for acoustic-beacon-optimiser tests."""

from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture
def sample_frequencies() -> np.ndarray:
    """Frequency array spanning glossophagine call range."""
    return np.linspace(30_000.0, 120_000.0, 50)


@pytest.fixture
def sample_angles() -> np.ndarray:
    """Incidence angle array from 0 to 90 degrees."""
    return np.linspace(0.0, np.pi / 2, 37)
