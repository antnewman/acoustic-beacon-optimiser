"""CMA-ES single-objective optimisation for reflector shapes.

Maximises integrated conspicuousness (IC) subject to a surface area
constraint, using the pycma ask-and-tell interface.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import NDArray


def optimise_cmaes(
    objective: Callable[[NDArray[np.float64]], float],
    x0: NDArray[np.float64],
    sigma0: float,
    bounds: list[tuple[float, float]],
    max_evaluations: int = 1000,
    options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run CMA-ES optimisation.

    Parameters
    ----------
    objective : Callable
        Function mapping parameter vector to scalar cost (minimised).
    x0 : NDArray[np.float64]
        Initial parameter vector.
    sigma0 : float
        Initial step size.
    bounds : list[tuple[float, float]]
        Parameter bounds [(lo, hi), ...].
    max_evaluations : int
        Maximum number of objective evaluations.
    options : dict, optional
        Additional pycma options.

    Returns
    -------
    dict
        Result dict with keys 'x_best', 'f_best', 'n_evals', 'history'.
    """
    raise NotImplementedError
