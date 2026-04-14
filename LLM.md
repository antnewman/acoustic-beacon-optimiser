# acoustic-beacon-optimiser -- LLM context file

This file provides a concise machine-readable summary of the project for use by large language models.

## Identity

- Name: acoustic-beacon-optimiser
- Version: 0.1.0 (early development)
- Author: Ant Newman
- Licence: MIT
- Repository: https://github.com/antnewman/acoustic-beacon-optimiser
- Language: Python 3.11+

## Purpose

Computes acoustic scattering from parameterised concave reflectors using the boundary element method (BEM) and optimises reflector shapes against biologically meaningful objective functions. The biological context is bat-pollinated flowers that have evolved concave structures to reflect echolocation signals back toward approaching bats.

## Core capabilities (planned)

1. Helmholtz BEM solver (bempp-cl, Burton-Miller CFIE) for rigid scatterers
2. Parametric reflector families: spherical cap, ellipsoidal cap, Chebyshev profile
3. Objective functions: target strength (TS), integrated conspicuousness (IC), catchment volume (CV), surface area (SA)
4. Single-objective optimisation via CMA-ES
5. Multi-objective Pareto frontier via NSGA-II
6. Validation against Mie series (analytical) and Simon et al. 2020 FEM results

## Package structure

```
src/abo/
  geometry/       reflectors.py, meshing.py, profiles.py
  acoustics/      bem_solver.py, target_strength.py, objectives.py
  optimisation/   single_objective.py, multi_objective.py, constraints.py
  biology/        call_spectra.py, natural_shapes.py, detection_model.py
  visualisation/  spectral_plots.py, pareto_plots.py, mesh_plots.py
  cli.py          Click-based CLI entry point (command: abo)
```

## Key dependencies

bempp-cl (BEM), gmsh (meshing), cma (CMA-ES), pymoo (NSGA-II), numpy, scipy, matplotlib, pyvista, click

## Commands

```
pip install -e ".[dev]"     # Install
abo --help                  # CLI
pytest                      # Tests
python -m ruff check src/   # Lint
mypy src/                   # Type check
```

## Scattering regime

Helmholtz number He = 2*pi*f*d/c in range 3-30 (resonance-to-optical transition). Reflector diameters 15-50 mm, frequencies 30-120 kHz. Full-wave BEM required; ray optics and Rayleigh approximation are not valid.

## Biological systems

- Marcgravia evenia: dish-leaf, ~35 mm diameter, ~25 mm depth
- Mucuna holtonii: vexillum petal, ~25 mm diameter, ~15 mm depth
- Bat pollinators: Glossophaga soricina, Leptonycteris yerbabuenae (FM calls 55-140 kHz)

## Documentation

- docs/SPEC.md: full technical specification and mathematical formulation
- docs/RESEARCH.md: literature review and gap analysis
