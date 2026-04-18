# acoustic-beacon-optimiser

[![Licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENCE)
[![Python: 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![Status: early development](https://img.shields.io/badge/status-early%20development-orange.svg)]()

## Summary

`acoustic-beacon-optimiser` computes acoustic scattering from parameterised concave reflectors using the boundary element method (BEM) and optimises reflector shapes against biologically weighted objective functions. The biological motivation comes from bat-pollinated flowers, which have evolved concave acoustic reflectors that attract echolocating pollinators; despite two decades of empirical work, no formal scattering model or systematic shape optimisation for these structures has been published. Key species of interest include the dish-shaped leaf of *Marcgravia evenia* and the vexillum petal of *Mucuna holtonii*, both of which produce strongly directional echo returns. The tool serves as a computational companion to a forthcoming academic paper on optimal floral acoustic reflector geometry.

## Why this matters

Empirical biology on bat-pollinated floral reflectors is over twenty years mature, and the mathematical machinery required to model them (boundary element solvers, evolutionary optimisation) is well established; yet nobody has connected the two disciplines in a rigorous computational framework. The relevant Helmholtz number regime (He approximately 3 to 30) places these reflectors squarely in the resonance-to-optical transition, where geometric or Rayleigh approximations break down and full-wave methods are required. A validated scattering model coupled with shape optimisation would allow quantitative testing of evolutionary hypotheses about reflector form and function.

## Key features

- Helmholtz BEM scattering computation via `bempp-cl` (Burton-Miller CFIE formulation)
- Parametric reflector families: spherical cap, ellipsoidal cap, general Chebyshev profile
- Biologically weighted objective functions: integrated conspicuousness, surface area
- Single-objective optimisation via CMA-ES (`cma`)
- Multi-objective Pareto frontier computation via NSGA-II (`pymoo`)
- Validation against the analytical Mie series (agreement to 0.03 dB at He = 5 and 10)
- Spectral directional pattern visualisation
- Command-line interface for end-to-end runs

## Roadmap

- **Phase 1** -- BEM solver and validation (complete; Mie agreement confirmed)
- **Phase 2** -- Parameterisation, call spectra, objectives, plots (complete)
- **Phase 3** -- CMA-ES, NSGA-II, CLI runners (complete)
- **Phase 4** -- Notebooks, paper manuscript, production optimisation runs (in progress)

## Installation

This project is in early development and is not yet available on PyPI.

```bash
git clone https://github.com/antnewman/acoustic-beacon-optimiser.git
cd acoustic-beacon-optimiser
pip install -e ".[dev]"
```

System prerequisites:

- Python 3.11 or later
- A working C/C++ compiler (for `numba` and `bempp-cl` kernels)
- OpenCL drivers are optional; `bempp-cl` falls back to Numba JIT on systems without them.

## Quickstart

### 1. Compute target strength for a Marcgravia-approximating cap

```bash
abo solve --radius 12 --depth 25 --freq-min 45000 --freq-max 100000 \
          --output results/marcgravia.npz
```

This meshes a spherical cap (sphere radius 12 mm, depth 25 mm), runs the BEM solver across the supplied frequency range and a default 19-angle grid, and saves the target-strength matrix to a NumPy `.npz`.

### 2. Optimise a spherical cap for Glossophaga conspicuousness

```bash
abo optimise --family spherical-cap --call glossophaga \
             --area-max 0.008 --max-evals 200 --seed 42 \
             --output results/sc_optimum.json
```

Runs CMA-ES with a 200-evaluation budget, maximising integrated conspicuousness subject to a 8000 mm^2 surface-area cap. Writes the best parameters, history, and evaluation count to JSON.

### 3. Compute a Pareto frontier

```bash
abo pareto --family spherical-cap --call glossophaga \
           --pop-size 20 --n-gen 20 --seed 42 \
           --output results/sc_pareto.npz
```

Runs NSGA-II over the `(-IC, SA)` objective pair and saves the Pareto set and front as a NumPy archive.

### 4. Use the Python API directly

```python
import numpy as np
from abo.geometry.meshing import mesh_spherical_cap
from abo.geometry.reflectors import SphericalCap
from abo.acoustics.target_strength import monostatic_target_strength
from abo.biology.call_spectra import glossophaga_call_spectrum
from abo.acoustics.objectives import integrated_conspicuousness

cap = SphericalCap(radius=0.012, half_angle=np.radians(55.0))
frequencies = np.linspace(45_000.0, 100_000.0, 8)
angles = np.linspace(0.0, np.pi / 3, 7)

grid = mesh_spherical_cap(cap, frequency=float(frequencies.max()))
ts = monostatic_target_strength(grid, frequencies, angles)
ic = integrated_conspicuousness(
    ts, frequencies, angles, glossophaga_call_spectrum(frequencies),
)
print(f"IC = {ic:.2f} dB")
```

## Notebooks

The `notebooks/` directory contains the scientific pipeline that feeds the paper figures:

| Notebook | Purpose | Approximate runtime |
|---|---|---|
| `01_sphere_validation.ipynb` | BEM vs Mie benchmark | 10 min |
| `02_marcgravia_reconstruction.ipynb` | Marcgravia dish-leaf directivity | 20-30 min |
| `03_shape_optimisation.ipynb` | CMA-ES single-objective run | 30 min |
| `04_pareto_analysis.ipynb` | NSGA-II Pareto frontier | 1-3 h |

Notebook outputs are not committed to the repository. Reviewers and readers can rerun them locally, or in a free cloud environment (Google Colab, Kaggle). A follow-up PR will add "Open in Colab" badges directly to each notebook header.

## Development

```bash
# Full test suite (excludes slow integration tests)
pytest -m "not slow"

# Lint
ruff check src/ tests/

# Type check
mypy src/
```

As of Phase 3 completion, the test suite has 67 passing unit tests covering the BEM solver, meshing, call spectra, objectives, plots, optimisation wrappers, and encoders.

## Technology stack

| Library | Purpose |
|---|---|
| `bempp-cl` | Helmholtz boundary element solver |
| `Gmsh` | Surface mesh generation |
| `cma` | CMA-ES single-objective optimisation |
| `pymoo` | NSGA-II multi-objective optimisation |
| NumPy / SciPy | Numerical computation |
| Matplotlib / PyVista | Visualisation and 3D rendering |
| Click | Command-line interface |

## Documentation

- [Technical specification](docs/SPEC.md) -- mathematical formulation, architecture, validation strategy
- [Research background](docs/RESEARCH.md) -- literature review and gap analysis (28 references)
- [API reference](docs/API.md) -- index of public functions and classes
- [LLM context summary](LLM.md) -- machine-readable project summary

## Citation

Citation guidance is provided in [CITATION.cff](CITATION.cff). This will be updated when the companion paper is accepted; the repository itself will be archived on Zenodo to mint a code DOI for citation.

## Licence

This project is released under the MIT licence. See [LICENCE](LICENCE) for details.

## Author

Ant Newman, CEng MIET MIEEE -- [tortoiseai.co.uk/about/ant-newman](https://tortoiseai.co.uk/about/ant-newman)
