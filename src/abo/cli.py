"""Command-line interface for acoustic-beacon-optimiser."""

from __future__ import annotations

import json
from pathlib import Path

import click
import numpy as np

from abo.acoustics.bem_solver import scattered_field, solve_scattering
from abo.acoustics.target_strength import monostatic_target_strength
from abo.biology.call_spectra import (
    flat_band_spectrum,
    glossophaga_call_spectrum,
    leptonycteris_call_spectrum,
)
from abo.geometry.meshing import mesh_spherical_cap
from abo.geometry.reflectors import SphericalCap
from abo.optimisation.multi_objective import pareto_frontier
from abo.optimisation.runners import (
    EvaluationConfig,
    build_encoder,
    make_ic_cost,
    make_pareto_objectives,
)
from abo.optimisation.single_objective import optimise_cmaes


def _call_spectrum(
    name: str, frequencies: np.ndarray,
) -> np.ndarray:
    if name == "glossophaga":
        return glossophaga_call_spectrum(frequencies)
    if name == "leptonycteris":
        return leptonycteris_call_spectrum(frequencies)
    if name == "flat":
        return flat_band_spectrum(frequencies)
    msg = f"Unknown call spectrum: {name!r}"
    raise click.ClickException(msg)


@click.group()
@click.version_option(package_name="acoustic-beacon-optimiser")
def main() -> None:
    """BEM solver and shape optimiser for floral acoustic reflectors."""


@main.command()
@click.option(
    "--radius", type=float, required=True, help="Sphere radius in mm.",
)
@click.option(
    "--depth", type=float, required=True, help="Cap depth in mm.",
)
@click.option(
    "--freq-min", type=float, default=30_000.0, help="Min frequency in Hz.",
)
@click.option(
    "--freq-max", type=float, default=120_000.0, help="Max frequency in Hz.",
)
@click.option(
    "--n-freq", type=int, default=20, help="Number of frequency points.",
)
@click.option(
    "--n-angles", type=int, default=19, help="Number of incidence angles.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output path for TS matrix (NumPy .npz).",
)
def solve(
    radius: float,
    depth: float,
    freq_min: float,
    freq_max: float,
    n_freq: int,
    n_angles: int,
    output: Path | None,
) -> None:
    """Compute target strength for a spherical cap reflector.

    Radius and depth are given in mm. Depth must be <= radius.
    """
    r = radius * 1e-3
    d = depth * 1e-3
    if d > r:
        msg = (
            f"Depth {depth} mm exceeds radius {radius} mm; "
            "cap is geometrically impossible."
        )
        raise click.ClickException(msg)
    theta = float(np.arccos(1.0 - d / r))
    cap = SphericalCap(radius=r, half_angle=theta)

    frequencies = np.linspace(freq_min, freq_max, n_freq)
    angles = np.linspace(0.0, np.pi / 2, n_angles)

    click.echo(
        f"Solving: r={radius} mm, d={depth} mm, theta={np.degrees(theta):.1f} deg",
    )
    click.echo(
        f"Frequency: {freq_min / 1000:.0f}-{freq_max / 1000:.0f} kHz"
        f" ({n_freq} points)",
    )
    click.echo(f"Angles: {n_angles} points from 0 to 90 deg")

    # Mesh at the highest frequency for adequate resolution
    grid = mesh_spherical_cap(cap, frequency=float(freq_max))
    click.echo(f"Mesh: {grid.number_of_elements} elements")

    ts = monostatic_target_strength(
        grid,
        frequencies,
        angles,
        far_field_range=1.0,
    )

    click.echo(
        f"TS range: [{ts.min():.2f}, {ts.max():.2f}] dB",
    )
    click.echo(f"TS on-axis mean: {ts[:, 0].mean():.2f} dB")

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            output,
            target_strength_db=ts,
            frequencies_hz=frequencies,
            angles_rad=angles,
            radius_m=r,
            half_angle_rad=theta,
        )
        click.echo(f"Saved to {output}")


@main.command()
@click.option(
    "--family",
    type=click.Choice(["spherical-cap", "ellipsoidal", "chebyshev"]),
    default="spherical-cap",
    help="Reflector parameterisation family.",
)
@click.option(
    "--call",
    type=click.Choice(["glossophaga", "leptonycteris", "flat"]),
    default="glossophaga",
    help="Bat call spectrum weighting.",
)
@click.option(
    "--area-max", type=float, default=0.010,
    help="Maximum surface area in m^2.",
)
@click.option(
    "--max-evals", type=int, default=200,
    help="Maximum objective evaluations.",
)
@click.option(
    "--n-coeffs", type=int, default=4,
    help="Chebyshev coefficient count (family=chebyshev only).",
)
@click.option(
    "--aperture", type=float, default=0.040,
    help="Chebyshev aperture radius in m (family=chebyshev only).",
)
@click.option(
    "--seed", type=int, default=None, help="Random seed.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output path for result JSON.",
)
def optimise(
    family: str,
    call: str,
    area_max: float,
    max_evals: int,
    n_coeffs: int,
    aperture: float,
    seed: int | None,
    output: Path | None,
) -> None:
    """Run single-objective CMA-ES shape optimisation.

    Maximises integrated conspicuousness subject to a surface-area
    constraint, using the selected bat-call spectrum as frequency
    weighting.
    """
    click.echo(
        f"CMA-ES: family={family}, call={call}, "
        f"A_max={area_max} m^2, budget={max_evals}",
    )

    if family == "chebyshev":
        encoder = build_encoder(
            family, n_coefficients=n_coeffs, aperture_radius=aperture,
        )
    else:
        encoder = build_encoder(family)

    frequencies = np.linspace(40_000.0, 100_000.0, 8)
    angles = np.linspace(0.0, np.pi / 3, 7)
    spectrum = _call_spectrum(call, frequencies)
    config = EvaluationConfig(
        frequencies=frequencies,
        angles=angles,
        call_spectrum=spectrum,
    )

    cost = make_ic_cost(encoder, config, area_max=area_max)

    result = optimise_cmaes(
        objective=cost,
        x0=encoder.x0,
        sigma0=0.25 * np.mean(
            [hi - lo for lo, hi in encoder.default_bounds],
        ),
        bounds=encoder.default_bounds,
        max_evaluations=max_evals,
        seed=seed,
    )

    click.echo(f"Best x:      {result.x_best.tolist()}")
    click.echo(f"Best cost:   {result.f_best:.4f} (IC = {-result.f_best:.2f} dB)")
    click.echo(f"Evaluations: {result.n_evals}")

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(
                {
                    "family": family,
                    "call": call,
                    "area_max": area_max,
                    "x_best": result.x_best.tolist(),
                    "f_best": result.f_best,
                    "n_evals": result.n_evals,
                    "history": result.history,
                },
                indent=2,
            ),
        )
        click.echo(f"Saved to {output}")


@main.command()
@click.option(
    "--family",
    type=click.Choice(["spherical-cap", "ellipsoidal", "chebyshev"]),
    default="spherical-cap",
    help="Reflector parameterisation family.",
)
@click.option(
    "--call",
    type=click.Choice(["glossophaga", "leptonycteris", "flat"]),
    default="glossophaga",
    help="Bat call spectrum weighting.",
)
@click.option(
    "--pop-size", type=int, default=20, help="NSGA-II population size.",
)
@click.option(
    "--n-gen", type=int, default=20, help="Number of generations.",
)
@click.option(
    "--n-coeffs", type=int, default=4,
    help="Chebyshev coefficient count (family=chebyshev only).",
)
@click.option(
    "--aperture", type=float, default=0.040,
    help="Chebyshev aperture radius in m (family=chebyshev only).",
)
@click.option(
    "--seed", type=int, default=None, help="Random seed.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output path for Pareto data (NumPy .npz).",
)
def pareto(
    family: str,
    call: str,
    pop_size: int,
    n_gen: int,
    n_coeffs: int,
    aperture: float,
    seed: int | None,
    output: Path | None,
) -> None:
    """Compute the (IC, SA) Pareto frontier via NSGA-II.

    Returns the set of non-dominated (integrated conspicuousness,
    surface area) trade-offs.
    """
    click.echo(
        f"NSGA-II: family={family}, call={call}, "
        f"pop={pop_size}, gen={n_gen}",
    )

    if family == "chebyshev":
        encoder = build_encoder(
            family, n_coefficients=n_coeffs, aperture_radius=aperture,
        )
    else:
        encoder = build_encoder(family)

    frequencies = np.linspace(40_000.0, 100_000.0, 8)
    angles = np.linspace(0.0, np.pi / 3, 7)
    spectrum = _call_spectrum(call, frequencies)
    config = EvaluationConfig(
        frequencies=frequencies,
        angles=angles,
        call_spectrum=spectrum,
    )

    objectives = make_pareto_objectives(encoder, config)

    result = pareto_frontier(
        objectives=objectives,
        bounds=encoder.default_bounds,
        population_size=pop_size,
        n_generations=n_gen,
        seed=seed,
    )

    click.echo(f"Pareto points: {result.pareto_set.shape[0]}")
    click.echo(f"Evaluations:   {result.n_evals}")
    # Front columns: (-IC, SA); convert for display
    ic_values = -result.pareto_front[:, 0]
    sa_values = result.pareto_front[:, 1]
    click.echo(
        f"IC range: [{ic_values.min():.2f}, {ic_values.max():.2f}] dB",
    )
    click.echo(
        f"SA range: [{sa_values.min():.4f}, {sa_values.max():.4f}] m^2",
    )

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            output,
            pareto_set=result.pareto_set,
            pareto_front=result.pareto_front,
            n_evals=result.n_evals,
        )
        click.echo(f"Saved to {output}")


__all__ = ["main", "scattered_field", "solve_scattering"]
