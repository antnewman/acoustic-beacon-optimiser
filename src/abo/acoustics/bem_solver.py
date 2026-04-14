"""Helmholtz BEM solver using bempp-cl.

Solves the exterior Helmholtz problem for a rigid (sound-hard) scatterer
using the Burton-Miller combined field integral equation (CFIE):

    (1/2 + D_k + alpha * H_k) phi = -(p_inc + alpha * dp_inc/dn)

with alpha = i/k (standard coupling parameter).
"""

from __future__ import annotations

import numpy as np
from bempp_cl.api import Grid as BemppGrid
from bempp_cl.api import (
    GridFunction,
    complex_callable,
    function_space,
    linalg,
    operators,
)
from numpy.typing import NDArray


def solve_scattering(
    grid: BemppGrid,
    frequency: float,
    incidence_direction: NDArray[np.float64],
    speed_of_sound: float = 343.0,
    gmres_tol: float = 1e-5,
) -> GridFunction:
    """Solve the Helmholtz scattering problem for a rigid reflector.

    Uses the Burton-Miller CFIE to avoid spurious interior resonances.
    The incident field is a unit-amplitude plane wave.

    Parameters
    ----------
    grid : BemppGrid
        Surface mesh of the reflector.
    frequency : float
        Excitation frequency in Hz.
    incidence_direction : NDArray[np.float64]
        Unit vector for the incident plane wave direction.
    speed_of_sound : float
        Speed of sound in m/s.
    gmres_tol : float
        GMRES convergence tolerance.

    Returns
    -------
    GridFunction
        Total field on the boundary surface.
    """
    k = 2.0 * np.pi * frequency / speed_of_sound
    d = np.asarray(incidence_direction, dtype=np.float64)
    d = d / np.linalg.norm(d)

    space = function_space(grid, "P", 1)

    # Boundary operators
    identity = operators.boundary.sparse.identity(space, space, space)
    dlp = operators.boundary.helmholtz.double_layer(
        space, space, space, k,
    )
    hyp = operators.boundary.helmholtz.hypersingular(
        space, space, space, k,
    )

    # Burton-Miller coupling parameter
    alpha = 1j / k

    # LHS: (1/2 I + D_k + alpha * H_k)
    lhs = 0.5 * identity + dlp + alpha * hyp

    # Incident field grid functions
    @complex_callable
    def _p_inc(x, n, domain_index, result):  # type: ignore[no-untyped-def]
        result[0] = np.exp(1j * k * np.dot(d, x))

    @complex_callable
    def _dp_inc_dn(x, n, domain_index, result):  # type: ignore[no-untyped-def]
        result[0] = (
            1j * k * np.dot(d, n) * np.exp(1j * k * np.dot(d, x))
        )

    p_inc_gf = GridFunction(space, fun=_p_inc)
    dp_inc_gf = GridFunction(space, fun=_dp_inc_dn)

    # RHS: -(p_inc + alpha * dp_inc/dn)
    rhs = -(p_inc_gf + alpha * dp_inc_gf)

    # Solve
    phi, info, *_ = linalg.gmres(
        lhs, rhs, tol=gmres_tol, return_iteration_count=True,
    )
    if info != 0:
        msg = f"GMRES did not converge (info={info})"
        raise RuntimeError(msg)

    return phi


def scattered_field(
    grid: BemppGrid,
    total_field: GridFunction,
    evaluation_points: NDArray[np.float64],
    frequency: float,
    incidence_direction: NDArray[np.float64],
    speed_of_sound: float = 343.0,
) -> NDArray[np.complex128]:
    """Evaluate the scattered pressure field at exterior points.

    Uses the Kirchhoff-Helmholtz representation formula:
        p_scat(x) = DLP[phi_scat](x) + SLP[dp_inc/dn](x)

    where phi_scat = phi_total - p_inc on the boundary, and
    dp_scat/dn = -dp_inc/dn (since dp_total/dn = 0 for Neumann BC).

    Parameters
    ----------
    grid : BemppGrid
        Surface mesh of the reflector.
    total_field : GridFunction
        BEM solution (total surface field from solve_scattering).
    evaluation_points : NDArray[np.float64]
        Points at which to evaluate, shape (3, N).
    frequency : float
        Frequency in Hz.
    incidence_direction : NDArray[np.float64]
        Unit vector for the incident plane wave direction.
    speed_of_sound : float
        Speed of sound in m/s.

    Returns
    -------
    NDArray[np.complex128]
        Scattered pressure at each evaluation point, shape (N,).
    """
    k = 2.0 * np.pi * frequency / speed_of_sound
    d = np.asarray(incidence_direction, dtype=np.float64)
    d = d / np.linalg.norm(d)
    pts = np.asarray(evaluation_points, dtype=np.float64)

    space = total_field.space

    # Incident field on the boundary
    @complex_callable
    def _p_inc(x, n, domain_index, result):  # type: ignore[no-untyped-def]
        result[0] = np.exp(1j * k * np.dot(d, x))

    @complex_callable
    def _dp_inc_dn(x, n, domain_index, result):  # type: ignore[no-untyped-def]
        result[0] = (
            1j * k * np.dot(d, n) * np.exp(1j * k * np.dot(d, x))
        )

    p_inc_gf = GridFunction(space, fun=_p_inc)
    dp_inc_gf = GridFunction(space, fun=_dp_inc_dn)

    # Scattered surface field
    phi_scat = total_field - p_inc_gf

    # Potential operators for evaluation at exterior points
    slp_pot = operators.potential.helmholtz.single_layer(
        space, pts, k,
    )
    dlp_pot = operators.potential.helmholtz.double_layer(
        space, pts, k,
    )

    # Kirchhoff-Helmholtz: p_scat = DLP[phi_scat] + SLP[dp_inc/dn]
    p_scat = dlp_pot @ phi_scat + slp_pot @ dp_inc_gf

    return p_scat.ravel().astype(np.complex128)
