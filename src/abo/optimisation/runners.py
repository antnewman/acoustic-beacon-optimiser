"""Optimisation runners that wire encoders, meshers, BEM, and objectives.

These functions turn a reflector family plus a call spectrum into a
callable cost function consumable by CMA-ES or NSGA-II.

The actual optimisation runs are expensive (many BEM solves) and
produce scientific figures rather than test output. They are invoked
from the CLI (`abo optimise`, `abo pareto`) and from notebooks.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Protocol

import numpy as np
from bempp_cl.api import Grid
from numpy.typing import NDArray

from abo.acoustics.objectives import (
    integrated_conspicuousness,
    surface_area,
)
from abo.acoustics.target_strength import monostatic_target_strength
from abo.geometry.meshing import (
    mesh_chebyshev_profile,
    mesh_ellipsoidal_cap,
    mesh_spherical_cap,
)
from abo.geometry.reflectors import (
    ChebyshevProfile,
    EllipsoidalCap,
    SphericalCap,
)
from abo.optimisation.constraints import (
    area_constraint,
    concavity_penalty,
)
from abo.optimisation.encoding import (
    ChebyshevProfileEncoding,
    EllipsoidalCapEncoding,
    SphericalCapEncoding,
)


class _ReflectorEncoder(Protocol):
    @property
    def dim(self) -> int: ...

    @property
    def default_bounds(self) -> list[tuple[float, float]]: ...

    @property
    def x0(self) -> NDArray[np.float64]: ...

    def decode(self, x: NDArray[np.float64]) -> object: ...


@dataclass(frozen=True)
class EvaluationConfig:
    """Configuration for the BEM evaluation grid used in objectives.

    Keep the frequency and angle grids coarse for optimisation; refine
    for final scoring of optimised shapes.
    """

    frequencies: NDArray[np.float64]
    angles: NDArray[np.float64]
    call_spectrum: NDArray[np.float64]
    speed_of_sound: float = 343.0
    far_field_range: float = 1.0
    gmres_tol: float = 1e-4
    elements_per_wavelength: float = 6.0


def _mesh_from_reflector(
    reflector: object,
    frequency: float,
    config: EvaluationConfig,
) -> Grid:
    """Dispatch to the appropriate mesher for the reflector type."""
    kwargs = {
        "elements_per_wavelength": config.elements_per_wavelength,
        "frequency": frequency,
        "speed_of_sound": config.speed_of_sound,
    }
    if isinstance(reflector, SphericalCap):
        return mesh_spherical_cap(reflector, **kwargs)
    if isinstance(reflector, EllipsoidalCap):
        return mesh_ellipsoidal_cap(reflector, **kwargs)
    if isinstance(reflector, ChebyshevProfile):
        return mesh_chebyshev_profile(reflector, **kwargs)
    msg = f"Unknown reflector type: {type(reflector)!r}"
    raise TypeError(msg)


def make_ic_cost(
    encoder: _ReflectorEncoder,
    config: EvaluationConfig,
    area_max: float | None = None,
    area_penalty_weight: float = 10.0,
    concavity_penalty_weight: float = 1000.0,
    infeasible_cost: float = 1e6,
) -> Callable[[NDArray[np.float64]], float]:
    """Construct a cost function: `-IC + penalties`, suitable for CMA-ES.

    The cost is the negative integrated conspicuousness (so lower is
    better from the optimiser's perspective) plus weighted penalty
    terms for area excess and concavity violation. Meshing or solve
    failures return `infeasible_cost`.

    Parameters
    ----------
    encoder : _ReflectorEncoder
        Encoder mapping parameter vectors to reflector objects.
    config : EvaluationConfig
        BEM evaluation grid.
    area_max : float or None
        Maximum allowed surface area in m^2. None disables the area
        penalty.
    area_penalty_weight : float
        Multiplier applied to the excess-area penalty (in m^2).
    concavity_penalty_weight : float
        Multiplier applied to the concavity penalty (for Chebyshev
        profiles; ignored for spherical/ellipsoidal caps).
    infeasible_cost : float
        Fallback cost returned when meshing or solving raises.

    Returns
    -------
    Callable
        Cost function `x -> float` ready for `optimise_cmaes`.
    """

    def cost(x: NDArray[np.float64]) -> float:
        try:
            reflector = encoder.decode(x)
            design_f = float(np.max(config.frequencies))
            grid = _mesh_from_reflector(reflector, design_f, config)

            ts = monostatic_target_strength(
                grid,
                config.frequencies,
                config.angles,
                speed_of_sound=config.speed_of_sound,
                far_field_range=config.far_field_range,
                gmres_tol=config.gmres_tol,
            )
            ic_db = integrated_conspicuousness(
                ts, config.frequencies, config.angles, config.call_spectrum,
            )

            penalty = 0.0
            if area_max is not None:
                penalty += area_penalty_weight * area_constraint(
                    grid, area_max,
                )
            if isinstance(encoder, ChebyshevProfileEncoding):
                penalty += concavity_penalty_weight * concavity_penalty(
                    x, encoder.aperture_radius,
                )

            return -ic_db + penalty
        except (RuntimeError, ValueError):
            return infeasible_cost

    return cost


def make_pareto_objectives(
    encoder: _ReflectorEncoder,
    config: EvaluationConfig,
    concavity_penalty_weight: float = 1000.0,
    infeasible_cost: float = 1e6,
) -> list[Callable[[NDArray[np.float64]], float]]:
    """Construct (negative IC, surface area) objectives for NSGA-II.

    Both objectives are minimised. The first is negative IC in dB
    (to maximise conspicuousness); the second is surface area in m^2
    (to minimise structural cost). Concavity is enforced via a penalty
    added to both objectives.

    Parameters
    ----------
    encoder : _ReflectorEncoder
        Encoder mapping parameter vectors to reflector objects.
    config : EvaluationConfig
        BEM evaluation grid.
    concavity_penalty_weight : float
        Multiplier applied to the concavity penalty (for Chebyshev
        profiles; ignored otherwise).
    infeasible_cost : float
        Fallback cost returned when meshing or solving raises.

    Returns
    -------
    list[Callable]
        Two-element list [neg_ic, sa] of cost functions.
    """

    def _evaluate(
        x: NDArray[np.float64],
    ) -> tuple[float, float]:
        try:
            reflector = encoder.decode(x)
            design_f = float(np.max(config.frequencies))
            grid = _mesh_from_reflector(reflector, design_f, config)

            ts = monostatic_target_strength(
                grid,
                config.frequencies,
                config.angles,
                speed_of_sound=config.speed_of_sound,
                far_field_range=config.far_field_range,
                gmres_tol=config.gmres_tol,
            )
            ic_db = integrated_conspicuousness(
                ts, config.frequencies, config.angles, config.call_spectrum,
            )
            sa = surface_area(grid)

            concav_pen = 0.0
            if isinstance(encoder, ChebyshevProfileEncoding):
                concav_pen = concavity_penalty_weight * concavity_penalty(
                    x, encoder.aperture_radius,
                )

            return -ic_db + concav_pen, sa + concav_pen
        except (RuntimeError, ValueError):
            return infeasible_cost, infeasible_cost

    def neg_ic(x: NDArray[np.float64]) -> float:
        return _evaluate(x)[0]

    def sa(x: NDArray[np.float64]) -> float:
        return _evaluate(x)[1]

    return [neg_ic, sa]


def build_encoder(family: str, **kwargs: object) -> _ReflectorEncoder:
    """Factory for the built-in reflector encoders.

    Parameters
    ----------
    family : str
        One of {"spherical-cap", "ellipsoidal", "chebyshev"}.
    **kwargs
        Passed to the encoder constructor (e.g. `n_coefficients` and
        `aperture_radius` for Chebyshev).

    Returns
    -------
    _ReflectorEncoder
        The encoder object.
    """
    if family == "spherical-cap":
        return SphericalCapEncoding(**kwargs)  # type: ignore[arg-type]
    if family == "ellipsoidal":
        return EllipsoidalCapEncoding(**kwargs)  # type: ignore[arg-type]
    if family == "chebyshev":
        return ChebyshevProfileEncoding(**kwargs)  # type: ignore[arg-type]
    msg = f"Unknown family: {family!r}"
    raise ValueError(msg)
