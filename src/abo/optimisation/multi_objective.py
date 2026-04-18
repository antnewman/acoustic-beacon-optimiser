"""Multi-objective optimisation for Pareto frontier computation.

Wraps pymoo's NSGA-II to find the Pareto frontier of the
(negative integrated conspicuousness, surface area) pair, yielding the
acoustic-cost / structural-cost trade-off for a reflector family.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pymoo.termination import get_termination


@dataclass
class ParetoResult:
    """Result of a multi-objective NSGA-II run.

    Parameters
    ----------
    pareto_set : NDArray[np.float64]
        Non-dominated parameter vectors, shape (N, n_params).
    pareto_front : NDArray[np.float64]
        Objective values at each Pareto point, shape (N, n_objectives).
        All objectives are minimised; callers that pass `-IC` receive
        it in the same convention.
    n_evals : int
        Total number of objective evaluations.
    """

    pareto_set: NDArray[np.float64]
    pareto_front: NDArray[np.float64]
    n_evals: int


class _ElementwiseObjective(ElementwiseProblem):
    """pymoo problem adapter evaluating a list of scalar objectives
    per individual."""

    def __init__(
        self,
        objectives: list[Callable[[NDArray[np.float64]], float]],
        bounds: list[tuple[float, float]],
    ) -> None:
        self._objectives = objectives
        xl = np.array([b[0] for b in bounds], dtype=np.float64)
        xu = np.array([b[1] for b in bounds], dtype=np.float64)
        super().__init__(
            n_var=len(bounds),
            n_obj=len(objectives),
            n_ieq_constr=0,
            xl=xl,
            xu=xu,
        )

    def _evaluate(
        self,
        x: NDArray[np.float64],
        out: dict[str, NDArray[np.float64]],
        *args: object,
        **kwargs: object,
    ) -> None:
        out["F"] = np.array(
            [float(obj(x)) for obj in self._objectives],
            dtype=np.float64,
        )


def pareto_frontier(
    objectives: list[Callable[[NDArray[np.float64]], float]],
    bounds: list[tuple[float, float]],
    population_size: int = 50,
    n_generations: int = 100,
    seed: int | None = None,
) -> ParetoResult:
    """Compute the Pareto frontier using NSGA-II.

    All objectives are minimised. To maximise integrated conspicuousness,
    pass `lambda x: -IC(x)` as the first objective; the resulting Pareto
    front will report `-IC` values which the caller can negate for
    display.

    Parameters
    ----------
    objectives : list[Callable]
        Each objective maps a parameter vector to a scalar to minimise.
    bounds : list[tuple[float, float]]
        Per-parameter box bounds.
    population_size : int
        NSGA-II population size.
    n_generations : int
        Number of generations.
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    ParetoResult
        Pareto-optimal parameter set, objective values, evaluation count.
    """
    problem = _ElementwiseObjective(objectives, bounds)
    algorithm = NSGA2(pop_size=population_size)
    termination = get_termination("n_gen", n_generations)

    result = minimize(
        problem,
        algorithm,
        termination,
        seed=seed,
        verbose=False,
    )

    pareto_set = np.asarray(result.X, dtype=np.float64)
    pareto_front = np.asarray(result.F, dtype=np.float64)
    # pymoo returns 1D for single-solution edge cases; normalise to 2D
    if pareto_set.ndim == 1:
        pareto_set = pareto_set.reshape(1, -1)
        pareto_front = pareto_front.reshape(1, -1)

    n_evals = int(result.algorithm.evaluator.n_eval)

    return ParetoResult(
        pareto_set=pareto_set,
        pareto_front=pareto_front,
        n_evals=n_evals,
    )
