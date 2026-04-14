"""Command-line interface for acoustic-beacon-optimiser."""

from __future__ import annotations

import click


@click.group()
@click.version_option(package_name="acoustic-beacon-optimiser")
def main() -> None:
    """BEM solver and shape optimiser for floral acoustic reflectors."""


@main.command()
@click.option("--radius", type=float, required=True, help="Sphere radius in mm.")
@click.option("--depth", type=float, required=True, help="Cap depth in mm.")
@click.option(
    "--freq-min", type=float, default=30_000.0, help="Min frequency in Hz.",
)
@click.option(
    "--freq-max", type=float, default=120_000.0, help="Max frequency in Hz.",
)
@click.option("--n-freq", type=int, default=50, help="Number of frequency points.")
@click.option(
    "--n-angles", type=int, default=37, help="Number of incidence angles.",
)
def solve(
    radius: float,
    depth: float,
    freq_min: float,
    freq_max: float,
    n_freq: int,
    n_angles: int,
) -> None:
    """Compute target strength for a spherical cap reflector."""
    click.echo(f"Solving: radius={radius} mm, depth={depth} mm")
    click.echo(
        f"Frequency: {freq_min / 1000:.0f}-{freq_max / 1000:.0f} kHz"
        f" ({n_freq} points)",
    )
    click.echo(f"Angles: {n_angles} points from 0 to 90 deg")
    raise NotImplementedError("BEM solver not yet implemented.")


@main.command()
@click.option(
    "--family",
    type=click.Choice(["spherical-cap", "ellipsoidal", "chebyshev"]),
    default="spherical-cap",
    help="Reflector parameterisation family.",
)
@click.option(
    "--objective",
    type=click.Choice(["ic", "ts"]),
    default="ic",
    help="Objective function.",
)
@click.option(
    "--area-max", type=float, default=3000.0, help="Max surface area in mm^2.",
)
@click.option(
    "--max-evals", type=int, default=500, help="Max objective evaluations.",
)
def optimise(
    family: str,
    objective: str,
    area_max: float,
    max_evals: int,
) -> None:
    """Run single-objective CMA-ES shape optimisation."""
    click.echo(f"Optimising: family={family}, objective={objective}")
    click.echo(f"Area constraint: {area_max} mm^2, max evals: {max_evals}")
    raise NotImplementedError("Optimisation not yet implemented.")


@main.command()
@click.option(
    "--family",
    type=click.Choice(["spherical-cap", "ellipsoidal", "chebyshev"]),
    default="spherical-cap",
    help="Reflector parameterisation family.",
)
@click.option("--pop-size", type=int, default=50, help="NSGA-II population size.")
@click.option("--n-gen", type=int, default=100, help="Number of generations.")
def pareto(
    family: str,
    pop_size: int,
    n_gen: int,
) -> None:
    """Compute Pareto frontier (IC vs SA) using NSGA-II."""
    click.echo(f"Pareto: family={family}, pop={pop_size}, gen={n_gen}")
    raise NotImplementedError("Multi-objective optimisation not yet implemented.")
