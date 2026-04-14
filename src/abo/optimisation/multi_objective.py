"""Multi-objective optimisation for Pareto frontier computation.

Uses pymoo NSGA-II to find the Pareto frontier of integrated
conspicuousness (IC) vs surface area (SA).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import NDArray


def pareto_frontier(
    objectives: list[Callable[[NDArray[np.float64]], float]],
    x0: NDArray[np.float64],
    bounds: list[tuple[float, float]],
    population_size: int = 50,
    n_generations: int = 100,
) -> dict[str, Any]:
    """Compute the Pareto frontier using NSGA-II.

    Parameters
    ----------
    objectives : list[Callable]
        List of objective functions (all minimised).
    x0 : NDArray[np.float64]
        Reference initial point.
    bounds : list[tuple[float, float]]
        Parameter bounds.
    population_size : int
        NSGA-II population size.
    n_generations : int
        Number of generations.

    Returns
    -------
    dict
        Result dict with keys 'pareto_set', 'pareto_front', 'n_evals'.
    """
    raise NotImplementedError
