# Project Specification

| Field   | Value            |
|---------|------------------|
| Version | 0.1.0-draft      |
| Author  | Ant Newman       |
| Date    | April 2026       |

---

## 1. Problem Statement

### 1.1 Physical Setting

Glossophagine (nectar-feeding) bats emit broadband frequency-modulated echolocation calls spanning 55--140 kHz, corresponding to wavelengths of approximately 2.4--6.2 mm in air. A number of bat-pollinated plants in the Neotropics have evolved concave structures that function as acoustic retroreflectors or beacons, returning a conspicuous echo signature to approaching pollinators.

Two canonical systems dominate the literature:

1. **Mucuna holtonii vexillum.** A concave petal (approximately 20--30 mm diameter) that acts as a cat's-eye retroreflector, returning a strong echo along the axis of incidence. First described by von Helversen and von Helversen (1999).

2. **Marcgravia evenia dish-leaf.** A concave modified leaf (approximately 35--50 mm diameter, 25--49 mm depth) positioned above the inflorescence that acts as a multidirectional acoustic beacon, enhancing detectability from a wide angular range. Described and experimentally verified by Simon et al. (2011).

### 1.2 Scattering Regime

The relevant non-dimensional frequency parameter (Helmholtz number) is:

```
He = 2 * pi * f * d / c
```

where `f` is frequency, `d` is reflector diameter, and `c` is the speed of sound. For the biological systems of interest, `He` spans the range 3--30, placing the scattering firmly in the resonance-to-optical transition regime where neither Rayleigh approximations nor geometric optics suffice.

### 1.3 Prior Work and Gaps

Simon et al. (2020) applied finite element method (FEM) modelling to spherical caps, exploring a 3x3 parameter grid of three radii (35, 50, 70 mm) and three depths (25, 30, 49 mm). This remains the only published numerical simulation of bat-flower acoustic reflectors.

The following gaps remain unaddressed:

- No boundary element method (BEM) computation has been applied to these structures.
- No parametric study beyond the 3x3 grid has been conducted.
- No shape optimisation of any kind has been attempted.
- No Pareto analysis trading off echo conspicuousness against reflector size has been performed.
- No open-source computational tool exists for studying these systems.

This project aims to fill all five gaps.

---

## 2. Mathematical Formulation

### 2.1 Forward Problem

We solve the exterior Helmholtz equation for acoustic scattering from a rigid body. The incident field is a unit-amplitude plane wave:

```
p_inc(x) = exp(i k n . x)
```

where `k = 2 pi f / c` is the wavenumber and `n` is the direction of incidence.

The total pressure field `p = p_inc + p_scat` satisfies:

1. **Helmholtz equation** in the exterior domain:

```
nabla^2 p + k^2 p = 0,  x in Omega_ext
```

2. **Neumann boundary condition** on the rigid scatterer surface (sound-hard):

```
dp/dn = 0  on Gamma
```

3. **Sommerfeld radiation condition** at infinity:

```
lim_{r -> inf} r (dp_scat/dr - i k p_scat) = 0
```

ensuring outgoing waves only.

#### Burton-Miller Combined Field Integral Equation

The standard boundary integral formulation for the exterior Neumann problem suffers from spurious solutions at interior resonance frequencies of the scatterer cavity. At these frequencies, the boundary integral operator becomes singular and the solution is non-unique. Since the concave reflectors studied here have cavity-like geometry, interior resonances will fall within the frequency band of interest and must be suppressed.

The Burton-Miller combined field integral equation (CFIE) resolves this by coupling the standard boundary integral equation with its normal derivative, yielding a uniquely solvable system at all frequencies:

```
(1/2 + D_k + alpha * H_k) phi = -(1/2 + D_k + alpha * H_k) p_inc
```

where:

- `D_k` is the double-layer boundary integral operator
- `H_k` is the hypersingular boundary integral operator
- `alpha = i / k` is the coupling parameter
- `phi` is the unknown surface potential

The coupling parameter `alpha = i / k` is the standard choice that balances conditioning across frequencies.

### 2.2 Reflector Parameterisation

Three families of axisymmetric reflector profiles are defined, offering increasing geometric flexibility.

#### Spherical Cap (2 parameters)

Parameterised by radius of curvature `r` and half-angle `theta`:

```
depth:    d = r (1 - cos(theta))
aperture: a = r sin(theta)
area:     A = 2 pi r^2 (1 - cos(theta))
```

This is the geometry studied by Simon et al. (2020) and serves as the baseline.

#### Ellipsoidal Cap (3 parameters)

Parameterised by semi-major axis `a_1`, semi-minor axis `a_2`, and half-angle `theta`. This family allows independent control of aperture width and depth curvature, capturing reflectors that are oblate or prolate relative to the spherical case.

#### General Chebyshev Profile (N parameters)

The meridional profile is represented as a Chebyshev series:

```
z(rho) = sum_{n=0}^{N-1} c_n T_n(2 rho / a - 1)
```

where `rho` is the radial coordinate, `a` is the aperture radius, `T_n` is the Chebyshev polynomial of the first kind of degree `n`, and `c_n` are the free coefficients.

The following constraints are imposed:

- `z(0) = d` (prescribed depth at centre)
- `z(a) = 0` (profile meets the aperture plane)
- `dz/drho < 0` for all `rho in (0, a)` (monotonically decreasing; no overhangs)

This family can approximate arbitrary smooth concave profiles, including biological shapes extracted from micro-CT data.

### 2.3 Objective Functions

#### Target Strength

The monostatic target strength at frequency `f` and incidence angle `theta_inc` is:

```
TS(f, theta_inc) = 10 * log10( |p_scat(r, theta_inc)|^2 * r^2 / |p_inc|^2 )
```

expressed in decibels (dB). This is the standard sonar/acoustics measure of echo strength.

#### Integrated Conspicuousness

The integrated conspicuousness captures the overall echo return weighted by the bat's call spectrum and integrated over solid angle:

```
IC = integral W(f) * TS(f, theta) * sin(theta) d(theta) df
```

where `W(f)` is a frequency weighting function derived from measured glossophagine call spectra (see Section 3, biology module). The `sin(theta)` factor accounts for the solid-angle element in spherical coordinates.

#### Catchment Volume

The catchment volume `V_c` is defined as the spatial volume within which the returning echo exceeds the bat's detection threshold:

```
V_c = { x in R^3 : echo_level(x) >= L_thresh }
```

The detection threshold `L_thresh` is estimated at -60 to -70 dB re 1 Pa, following psychoacoustic measurements in the literature (Yovel et al., 2008).

#### Surface Area

The surface area `SA` of the reflector serves as a cost proxy, representing the metabolic investment the plant must make to grow the structure.

#### Optimisation Problems

Two formulations are considered:

1. **Single-objective:** Maximise `IC` subject to `SA <= SA_max`, where `SA_max` is a prescribed upper bound on reflector area.

2. **Multi-objective:** Compute the Pareto frontier of `(IC, SA)`, revealing the trade-off between echo conspicuousness and reflector cost. Natural geometries are plotted on this frontier to assess their optimality.

### 2.4 Optimiser Selection

**Single-objective: CMA-ES** (Covariance Matrix Adaptation Evolution Strategy), via the `pycma` library. CMA-ES is a derivative-free global optimiser well-suited to moderate-dimensional (2--20 parameters), noisy, non-convex landscapes. It adapts a full covariance matrix, making it effective for correlated parameters such as Chebyshev coefficients. It requires only objective function evaluations, which aligns well with the BEM forward solver acting as a black box.

**Multi-objective: NSGA-II** (Non-dominated Sorting Genetic Algorithm II), via the `pymoo` library. NSGA-II is the standard algorithm for generating Pareto frontiers in engineering optimisation. It maintains a population of non-dominated solutions and uses crowding distance to ensure uniform coverage of the frontier. The `pymoo` implementation provides a well-tested, extensible framework with built-in convergence metrics.

---

## 3. Architecture

### 3.1 Directory Structure

```
acoustic-beacon-optimiser/
|-- src/abo/
|   |-- __init__.py
|   |-- geometry/
|   |   |-- __init__.py
|   |   |-- reflectors.py       # Reflector surface classes (SphericalCap, EllipsoidalCap, ChebyshevProfile)
|   |   |-- meshing.py          # Gmsh-based surface mesh generation
|   |   |-- profiles.py         # Meridional profile utilities and constraints
|   |-- acoustics/
|   |   |-- __init__.py
|   |   |-- bem_solver.py       # Burton-Miller CFIE solver wrapping bempp-cl
|   |   |-- target_strength.py  # TS computation from far-field patterns
|   |   |-- objectives.py       # IC, catchment volume, and composite objectives
|   |-- optimisation/
|   |   |-- __init__.py
|   |   |-- single_objective.py # CMA-ES wrapper via pycma
|   |   |-- multi_objective.py  # NSGA-II wrapper via pymoo
|   |   |-- constraints.py      # Bound and nonlinear constraint definitions
|   |-- biology/
|   |   |-- __init__.py
|   |   |-- call_spectra.py     # Glossophagine call frequency weighting W(f)
|   |   |-- natural_shapes.py   # Digitised natural reflector profiles
|   |   |-- detection_model.py  # Echo detection threshold and catchment volume
|   |-- visualisation/
|   |   |-- __init__.py
|   |   |-- spectral_plots.py   # TS vs frequency plots
|   |   |-- pareto_plots.py     # Pareto frontier visualisation
|   |   |-- mesh_plots.py       # 3D mesh rendering via PyVista
|   |-- cli.py                  # Click-based command-line interface
|-- tests/
|   |-- unit/                   # Unit tests (per-module)
|   |-- integration/            # End-to-end solver and optimisation tests
|-- notebooks/
|   |-- 01_mie_validation.ipynb
|   |-- 02_spherical_cap_study.ipynb
|   |-- 03_optimisation_runs.ipynb
|   |-- 04_pareto_analysis.ipynb
|-- data/
|   |-- natural_reflectors/     # Digitised profiles and micro-CT meshes
|   |-- call_parameters/        # Bat call spectra and audiogram data
|-- paper/
|   |-- figures/
|   |-- main.tex
|   |-- references.bib
|-- docs/
|   |-- SPEC.md
|-- pyproject.toml
|-- README.md
|-- LICENSE
```

---

## 4. Technology Stack

| Component              | Choice         | Rationale                                                                                          |
|------------------------|----------------|----------------------------------------------------------------------------------------------------|
| Language               | Python 3.11+   | Ecosystem maturity for scientific computing; type hints and performance improvements in 3.11+       |
| BEM solver             | bempp-cl       | Open-source boundary element library with OpenCL acceleration; native Burton-Miller CFIE support    |
| Mesh generation        | Gmsh           | Scriptable, parametric meshing with Python API; fine control over element size and curvature        |
| Single-obj. optimiser  | pycma (CMA-ES) | Derivative-free global optimisation; handles noisy, non-convex landscapes with correlated parameters|
| Multi-obj. optimiser   | pymoo (NSGA-II)| Standard Pareto frontier algorithm; extensible framework with convergence diagnostics               |
| Numerical computing    | NumPy          | Foundational array operations; universal dependency across the scientific Python stack              |
| Scientific routines    | SciPy          | Special functions (spherical Bessel/Hankel for Mie series), integration, interpolation              |
| 2D plotting            | Matplotlib     | Publication-quality figures; extensive customisation for journal submission                         |
| 3D visualisation       | PyVista        | VTK-based mesh rendering; interactive and scriptable 3D surface plots                              |
| CLI framework          | Click          | Composable command-line interface with automatic help generation and parameter validation           |
| Testing                | pytest         | Fixtures, parametrised tests, and coverage reporting; de facto standard for Python testing          |
| Linting                | Ruff           | Fast, comprehensive linter and formatter replacing flake8, isort, and black                        |
| Type checking          | mypy           | Static type analysis; catches interface errors before runtime                                      |

---

## 5. Validation Strategy

Validation proceeds at three levels of increasing physical fidelity.

### 5.1 Level 1: Analytical (Mie Series)

The BEM solver is validated against the analytical Mie series solution for scattering from a rigid sphere. Tests are conducted at four Helmholtz numbers:

- `He = 5`
- `He = 10`
- `He = 20`
- `He = 30`

These span the full range relevant to the biological reflectors. For each case, the far-field scattering pattern and target strength are compared. The acceptance criterion is agreement within **0.5 dB** across all angles.

The Mie series coefficients are computed using spherical Bessel and Hankel functions from SciPy, providing a fully independent reference.

### 5.2 Level 2: Numerical Cross-Validation (Simon et al., 2020)

The BEM solver results are compared against the FEM simulations of Simon et al. (2020) for their nine spherical-cap configurations:

| Configuration | Radius (mm) | Depth (mm) |
|---------------|-------------|------------|
| 1             | 35          | 25         |
| 2             | 35          | 30         |
| 3             | 35          | 49         |
| 4             | 50          | 25         |
| 5             | 50          | 30         |
| 6             | 50          | 49         |
| 7             | 70          | 25         |
| 8             | 70          | 30         |
| 9             | 70          | 49         |

Key quantities compared: impulse response (IR) peak amplitudes, broadband target strength, and spectral notch positions. These comparisons verify that the BEM and FEM approaches yield consistent results for the same geometries.

### 5.3 Level 3: Empirical (Simon et al., 2011)

The BEM predictions for the Marcgravia evenia dish-leaf geometry are compared against the experimental ensonification measurements reported by Simon et al. (2011). The acceptance criterion is agreement within **3 dB**, reflecting the combined uncertainty of the biological specimen geometry, measurement setup, and material property assumptions.

---

## 6. Development Phases

### Phase 1: BEM Solver and Validation (Weeks 1--3)

1. Set up project repository, CI pipeline, and development environment.
2. Implement Gmsh-based parametric meshing for spherical caps with adaptive element sizing.
3. Implement Burton-Miller CFIE solver wrapping bempp-cl, including far-field evaluation.
4. Implement target strength computation from far-field patterns.
5. Implement Mie series reference solution for rigid spheres (Level 1 validation).
6. Run Level 1 validation at `He = 5, 10, 20, 30`; confirm agreement within 0.5 dB.

### Phase 2: Parameterisation and Objectives (Weeks 4--5)

7. Implement ellipsoidal cap parameterisation and meshing.
8. Implement general Chebyshev profile parameterisation with monotonicity constraints.
9. Implement integrated conspicuousness (IC) objective with configurable call spectra.
10. Implement catchment volume computation and detection threshold model.

### Phase 3: Optimisation (Weeks 6--8)

11. Implement CMA-ES single-objective optimisation loop (max IC subject to SA constraint).
12. Implement NSGA-II multi-objective optimisation loop (IC vs SA Pareto frontier).
13. Run Level 2 validation against Simon et al. (2020) nine configurations.
14. Run Level 3 validation against Simon et al. (2011) Marcgravia measurements.
15. Execute optimisation campaigns across the three parameterisation families.

### Phase 4: Analysis and Paper (Weeks 9--12)

16. Analyse Pareto frontiers; plot natural geometries on the frontier to assess biological optimality.
17. Generate publication-quality figures (spectral plots, Pareto diagrams, 3D mesh renderings).
18. Prepare Jupyter notebooks documenting the full computational workflow.
19. Draft manuscript following the outline in Section 8.
20. Submit to target journal; release code and data as open-source.

---

## 7. Risks and Mitigations

### 7.1 BEM Computational Cost at High Frequencies

**Risk.** At the upper end of the frequency range (140 kHz, `He ~ 30`), the BEM mesh must resolve the surface to at least 6 elements per wavelength, leading to large dense matrices and long solve times. A single frequency solve could take minutes, and optimisation requires hundreds of evaluations.

**Mitigation.** Use bempp-cl's OpenCL-accelerated assembly for GPU offloading. Employ hierarchical matrix (H-matrix) compression where supported. Limit the Chebyshev profile order to keep the optimisation dimension manageable (N <= 8). If necessary, use a coarse-frequency sweep for optimisation and refine only near the optimum.

### 7.2 bempp-cl Installation Complexity

**Risk.** bempp-cl depends on OpenCL drivers, Numba, and several compiled dependencies, making installation non-trivial across platforms.

**Mitigation.** Provide a Docker container with all dependencies pre-installed. Document tested platform configurations (Ubuntu 22.04 with Intel/NVIDIA OpenCL). Include a Makefile target for environment setup.

### 7.3 Optimisation Convergence

**Risk.** The objective landscape may be multimodal or exhibit flat regions, causing CMA-ES to converge slowly or to local optima. NSGA-II may produce sparse Pareto fronts.

**Mitigation.** Run CMA-ES with multiple restarts (BIPOP-CMA-ES variant). For NSGA-II, use a population size of at least 100 and run for a minimum of 200 generations. Validate convergence by comparing independent runs.

### 7.4 Insufficient Validation Data

**Risk.** Only one published numerical dataset (Simon et al., 2020) and one experimental dataset (Simon et al., 2011) are available for comparison, limiting the strength of validation.

**Mitigation.** Supplement with the analytical Mie series (which admits arbitrary precision). Where possible, compare against additional published TS measurements for simple concave reflectors from the underwater acoustics or radar cross-section literature.

---

## 8. Paper Outline

**Target journal.** The Journal of the Acoustical Society of America (JASA). Alternative: Bioinspiration and Biomimetics.

**Working title.** "Shape optimisation of acoustic beacons in bat-pollinated flowers: a boundary element method study"

### Proposed Structure

1. **Introduction.** Bat-flower acoustic communication; prior experimental and modelling work; statement of the gap (no BEM, no optimisation, no Pareto analysis).

2. **Mathematical formulation.** Exterior Helmholtz problem; Burton-Miller CFIE; reflector parameterisation (spherical, ellipsoidal, Chebyshev); objective function definitions.

3. **Numerical methods and validation.** BEM implementation details; mesh convergence study; Mie series validation; cross-validation with Simon et al. (2020); comparison with Simon et al. (2011) experimental data.

4. **Optimisation methodology.** CMA-ES for single-objective; NSGA-II for multi-objective; constraint handling; convergence criteria.

5. **Results and discussion.** Parametric study of spherical caps (extending Simon et al.); Pareto frontiers for each parameterisation family; natural geometries plotted on the Pareto frontier; biological interpretation of optimality (or near-optimality).

6. **Conclusions.** Summary of findings; implications for understanding plant-bat coevolution; potential biomimetic applications; open-source tool availability.

**Key result.** The Pareto frontier of integrated conspicuousness vs surface area, with the natural Mucuna holtonii and Marcgravia evenia geometries plotted, demonstrating whether evolution has converged on Pareto-optimal or near-optimal acoustic designs.

---

## 9. References

1. von Helversen, D. and von Helversen, O. (1999). Acoustic guide in bat-pollinated flower. *Nature*, 398, 759--760. DOI: [10.1038/19648](https://doi.org/10.1038/19648)

2. von Helversen, O. and von Helversen, D. (2003). Acoustic guide in bat-pollinated flower. *Journal of Comparative Physiology A*, 189, 327--336. DOI: [10.1007/s00359-003-0405-3](https://doi.org/10.1007/s00359-003-0405-3)

3. von Helversen, D., Holderied, M.W. and von Helversen, O. (2003). Echoes of bat-pollinated bell-shaped flowers: conspicuous for nectar-feeding bats? *Journal of Experimental Biology*, 206, 1025--1034. DOI: [10.1242/jeb.00203](https://doi.org/10.1242/jeb.00203)

4. Von Helversen, D. (2004). Object classification by echolocation in nectar feeding bats: size-independent generalisation of shape. *Journal of Comparative Physiology A*, 190, 515--521. DOI: [10.1007/s00359-004-0492-9](https://doi.org/10.1007/s00359-004-0492-9)

5. Simon, R., Holderied, M.W. and von Helversen, O. (2006). Size discrimination of hollow hemispheres by echolocation in a nectar feeding bat. *Journal of Experimental Biology*, 209, 3599--3609. DOI: [10.1242/jeb.02398](https://doi.org/10.1242/jeb.02398)

6. Yovel, Y., Franz, M.O., Stilz, P. and Schnitzler, H.-U. (2008). Plant classification from bat-like echolocation signals. *PLoS Computational Biology*, 4, e1000032. DOI: [10.1371/journal.pcbi.1000032](https://doi.org/10.1371/journal.pcbi.1000032)

7. Simon, R., Holderied, M.W., Koch, C.U. and von Helversen, O. (2011). Floral acoustics: conspicuous echoes of a dish-shaped leaf attract bat pollinators. *Science*, 333, 631--633. DOI: [10.1126/science.1204210](https://doi.org/10.1126/science.1204210)

8. Schoener, C.R., Schoener, M.G. and Kerth, G. (2015). Similar is not the same: social calls of conspecifics are more effective in attracting wild bats to day roosts than those of other bat species. *Current Biology*, 25, 1911--1916. DOI: [10.1016/j.cub.2015.05.054](https://doi.org/10.1016/j.cub.2015.05.054)

9. Gonzalez-Terrazas, T.P., Koblitz, J.C., Fleming, T.H., Medellín, R.A., Kalko, E.K.V., Schnitzler, H.-U. and Tschapka, M. (2016). How nectar-feeding bats localize their food: echolocation behaviour and search strategies of *Glossophaga soricina*. *PLoS ONE*, 11, e0163492. DOI: [10.1371/journal.pone.0163492](https://doi.org/10.1371/journal.pone.0163492)

10. Simon, R., Holderied, M.W. and Koch, C.U. (2020). Echo-acoustic scanning with noseleaf and ears in phyllostomid bats. *Proceedings of the National Academy of Sciences*, 117, 1367--1374. DOI: [10.1073/pnas.1909890117](https://doi.org/10.1073/pnas.1909890117)

11. Simon, R., Rupitsch, S.J. and Guenther, F. (2021). Modelling the acoustic reflector of *Marcgravia evenia*. *PLoS Computational Biology*, 17, e1009706. DOI: [10.1371/journal.pcbi.1009706](https://doi.org/10.1371/journal.pcbi.1009706)

12. Simon, R., Guenther, F. and Rupitsch, S.J. (2023). The acoustic properties of bat-pollinated flowers: a modelling approach. *Journal of Experimental Biology*, 226, jeb245263. DOI: [10.1242/jeb.245263](https://doi.org/10.1242/jeb.245263)

13. Jansen, R. and Steckel, J. (2024). Bio-inspired sonar reflector design for enhanced target detection. *IEEE Access*, 12, 69371--69382. DOI: [10.1109/ACCESS.2024.3401226](https://doi.org/10.1109/ACCESS.2024.3401226)
