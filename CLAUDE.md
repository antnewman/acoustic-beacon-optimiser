# Claude Code Instructions - acoustic-beacon-optimiser

## Project Overview
Python library and CLI tool for computing acoustic scattering from parameterised concave reflectors using BEM, and optimising those shapes against biologically meaningful objective functions. Companion to a forthcoming academic paper.

## Tech Stack
- **Language:** Python 3.11+
- **BEM solver:** bempp-cl
- **Mesh generation:** Gmsh (Python API)
- **Optimiser:** cma (CMA-ES), pymoo (NSGA-II)
- **Numerics:** NumPy, SciPy
- **Visualisation:** Matplotlib, PyVista
- **CLI:** Click
- **Testing:** pytest
- **Linting:** Ruff
- **Type checking:** mypy (strict)

## Key Commands
- **Install:** `pip install -e ".[dev]"`
- **Tests:** `pytest` (skip slow: `pytest -m "not slow"`)
- **Lint:** `ruff check src/ tests/`
- **Type check:** `mypy src/`
- **CLI:** `abo --help`

## Working with Claude
1. For ALL git commits, use author "antnewman <antjsnewman@outlook.com>" -- NEVER include a "Co-Authored-By: Claude" line
2. Branch workflow: work on feature branches off `main`. Never push directly to `main`.
3. Always run `ruff check` after changes
4. Follow existing patterns: type annotations, frozen dataclasses, NumPy docstrings
5. Use British English throughout all prose (optimisation, behaviour, parameterised)
6. No em dashes. Use hyphens, semicolons, or parentheses instead.
7. Use conventional commit format for commit messages

## Package Structure
- `src/abo/geometry/` -- Parametric reflector definitions and mesh generation
- `src/abo/acoustics/` -- Helmholtz BEM solver, target strength, objectives
- `src/abo/optimisation/` -- CMA-ES and NSGA-II wrappers
- `src/abo/biology/` -- Bat call spectra, natural shapes, detection models
- `src/abo/visualisation/` -- Plotting utilities
- `src/abo/cli.py` -- Click CLI entry point

## Documentation
- `docs/SPEC.md` -- Full technical specification and mathematical formulation
- `docs/RESEARCH.md` -- Literature review and gap analysis
- `LLM.md` -- Machine-readable project summary for LLM context
