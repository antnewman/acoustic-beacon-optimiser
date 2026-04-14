"""Pareto frontier visualisation.

Plots the multi-objective Pareto frontier (IC vs SA) with natural
reflector geometries overlaid for comparison.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def plot_pareto_frontier(
    pareto_front: NDArray[np.float64],
    natural_points: dict[str, tuple[float, float]] | None = None,
) -> Figure:
    """Plot a 2D Pareto frontier with optional natural reflector overlay.

    Parameters
    ----------
    pareto_front : NDArray[np.float64]
        Pareto front points, shape (N, 2) with columns [IC, SA].
    natural_points : dict, optional
        Named natural reflector positions {name: (IC, SA)}.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    raise NotImplementedError
