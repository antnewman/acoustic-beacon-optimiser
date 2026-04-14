# acoustic-beacon-optimiser

[![Licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENCE)
[![Python: 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![Status: early development](https://img.shields.io/badge/status-early%20development-orange.svg)]()

## Summary

`acoustic-beacon-optimiser` computes acoustic scattering from parameterised concave reflectors using the boundary element method (BEM) and optimises reflector shapes against biologically weighted objective functions. The biological motivation comes from bat-pollinated flowers, which have evolved concave acoustic reflectors that attract echolocating pollinators; despite two decades of empirical work, no formal scattering model or systematic shape optimisation for these structures has been published. Key species of interest include the dish-shaped leaf of *Marcgravia evenia* and the vexillum petal of *Mucuna holtonii*, both of which produce strongly directional echo returns. The tool serves as a computational companion to a forthcoming academic paper on optimal floral acoustic reflector geometry.

## Why this matters

Empirical biology on bat-pollinated floral reflectors is over twenty years mature, and the mathematical machinery required to model them (boundary element solvers, evolutionary optimisation) is well established; yet nobody has connected the two disciplines in a rigorous computational framework. The relevant Helmholtz number regime (He approximately 3 to 30) places these reflectors squarely in the resonance-to-optical transition, where geometric or Rayleigh approximations break down and full-wave methods are required. A validated scattering model coupled with shape optimisation would allow quantitative testing of evolutionary hypotheses about reflector form and function.

## Key features (planned)

- Helmholtz BEM scattering computation via bempp-cl (Burton-Miller CFIE formulation)
- Parametric reflector families: spherical cap, ellipsoidal cap, general Chebyshev profile
- Biologically weighted objective functions: integrated conspicuousness, catchment volume, surface area
- Single-objective optimisation via CMA-ES (pycma)
- Multi-objective Pareto frontier computation via NSGA-II (pymoo)
- Validation against analytical Mie series and published FEM results (Simon et al., 2020)
- Spectral directional pattern visualisation

## Roadmap

- **Phase 1** -- BEM solver and validation
- **Phase 2** -- Parameterisation and objectives
- **Phase 3** -- Shape optimisation and analysis

## Installation

This project is in early development and is not yet available on PyPI.

```bash
git clone https://github.com/antnewman/acoustic-beacon-optimiser.git
cd acoustic-beacon-optimiser
pip install -e ".[dev]"
```

## Technology stack

| Library | Purpose |
|---------|---------|
| bempp-cl | Helmholtz boundary element solver |
| Gmsh | Surface mesh generation |
| pycma | CMA-ES single-objective optimisation |
| pymoo | NSGA-II multi-objective optimisation |
| NumPy / SciPy | Numerical computation |
| Matplotlib / PyVista | Visualisation and 3D rendering |
| Click | Command-line interface |

## Documentation

- [Technical specification](docs/SPEC.md)
- [Research background](docs/RESEARCH.md)

## Citation

Citation guidance is provided in [CITATION.cff](CITATION.cff). This will be updated when the companion paper is published.

## Licence

This project is released under the MIT licence. See [LICENCE](LICENCE) for details.

## Author

Ant Newman, CEng MIET MIEEE -- [tortoiseai.co.uk/about/ant-newman](https://tortoiseai.co.uk/about/ant-newman)
