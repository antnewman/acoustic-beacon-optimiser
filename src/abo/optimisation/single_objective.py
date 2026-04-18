"""CMA-ES single-objective optimisation for reflector shapes.

Thin wrapper around the `cma` package (pycma) providing an ask-and-tell
driver that returns the best-found parameters and fitness history for
an arbitrary scalar objective.

The scientific use case is: given a parameter-to-objective mapping
`x -> cost`, find `x*` that minimises cost subject to box bounds. By
convention the objective is a cost to minimise; to maximise
integrated conspicuousness, pass `lambda x: -IC(x)` as the objective.

Constraint handling is via penalty: callers add the concavity or
surface-area penalties to the returned cost before calling `tell`.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import cma
import numpy as np
from numpy.typing import NDArray


@dataclass
class OptimisationResult:
    """Result of a single-objective CMA-ES run.

    Parameters
    ----------
    x_best : NDArray[np.float64]
        Best parameter vector found.
    f_best : float
        Best (lowest) objective value found.
    n_evals : int
        Total number of objective evaluations.
    history : list[float]
        Best-so-far cost at each generation.
    """

    x_best: NDArray[np.float64]
    f_best: float
    n_evals: int
    history: list[float]


def optimise_cmaes(
    objective: Callable[[NDArray[np.float64]], float],
    x0: NDArray[np.float64],
    sigma0: float,
    bounds: list[tuple[float, float]] | None = None,
    max_evaluations: int = 1000,
    population_size: int | None = None,
    options: dict[str, Any] | None = None,
    seed: int | None = None,
) -> OptimisationResult:
    """Run CMA-ES optimisation using the ask-and-tell API.

    Parameters
    ----------
    objective : Callable
        Function mapping parameter vector to scalar cost (minimised).
    x0 : NDArray[np.float64]
        Initial parameter vector.
    sigma0 : float
        Initial step size (standard deviation). Should be roughly 1/4
        of the expected parameter range.
    bounds : list[tuple[float, float]] or None
        Per-parameter box bounds [(lo, hi), ...]. Use None for
        unbounded optimisation.
    max_evaluations : int
        Maximum number of objective evaluations.
    population_size : int or None
        Lambda (population size). None uses the pycma default
        `4 + floor(3 * log(dim))`.
    options : dict or None
        Additional pycma options; merged with defaults.
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    OptimisationResult
        Best parameters, best cost, evaluation count, and history.
    """
    x0_arr = np.asarray(x0, dtype=np.float64)

    opts: dict[str, Any] = {
        "maxfevals": max_evaluations,
        "verbose": -9,
    }
    if bounds is not None:
        lower = [b[0] for b in bounds]
        upper = [b[1] for b in bounds]
        opts["bounds"] = [lower, upper]
    if population_size is not None:
        opts["popsize"] = population_size
    if seed is not None:
        opts["seed"] = seed
    if options is not None:
        opts.update(options)

    es = cma.CMAEvolutionStrategy(x0_arr.tolist(), sigma0, opts)

    history: list[float] = []
    n_evals = 0
    best_x: NDArray[np.float64] = x0_arr.copy()
    best_f = float("inf")

    while not es.stop():
        xs = es.ask()
        fitnesses = [float(objective(np.asarray(x, dtype=np.float64)))
                     for x in xs]
        es.tell(xs, fitnesses)
        n_evals += len(xs)

        gen_best_idx = int(np.argmin(fitnesses))
        gen_best = fitnesses[gen_best_idx]
        if gen_best < best_f:
            best_f = gen_best
            best_x = np.asarray(xs[gen_best_idx], dtype=np.float64)
        history.append(best_f)

        if n_evals >= max_evaluations:
            break

    return OptimisationResult(
        x_best=best_x,
        f_best=best_f,
        n_evals=n_evals,
        history=history,
    )
