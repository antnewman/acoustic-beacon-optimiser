# API reference

This document indexes the public API of `acoustic-beacon-optimiser` for
humans reading the code. It mirrors the NumPy-style docstrings on each
function and class.

## `abo.geometry`

### Reflector classes -- `abo.geometry.reflectors`

| Class | Parameters | Purpose |
|---|---|---|
| `SphericalCap` | `radius`, `half_angle` | Spherical cap with two-parameter family |
| `EllipsoidalCap` | `semi_major`, `semi_minor`, `half_angle` | Oblate/prolate ellipsoidal cap |
| `ChebyshevProfile` | `coefficients`, `aperture_radius` | Axisymmetric profile defined by truncated Chebyshev series |

Every class exposes the properties `depth`, `aperture`, `surface_area`.
`ChebyshevProfile` additionally exposes `is_concave()` and `profile(rho)`.

### Meshing -- `abo.geometry.meshing`

| Function | Input | Returns |
|---|---|---|
| `mesh_spherical_cap(cap, ...)` | `SphericalCap` | `bempp_cl.api.Grid` |
| `mesh_ellipsoidal_cap(cap, ...)` | `EllipsoidalCap` | `bempp_cl.api.Grid` |
| `mesh_chebyshev_profile(profile, ...)` | `ChebyshevProfile` | `bempp_cl.api.Grid` |

Mesh density is controlled by `elements_per_wavelength`, `frequency`,
`speed_of_sound`. Spherical and ellipsoidal caps use Gmsh OCC Boolean
operations; the Chebyshev profile uses direct structured triangulation
of the surface of revolution.

### Profile utilities -- `abo.geometry.profiles`

- `chebyshev_profile(rho, coefficients, aperture_radius)` -- evaluate
  the profile z(rho).
- `chebyshev_profile_derivative(rho, coefficients, aperture_radius)`
  -- evaluate dz/drho.
- `is_concave(coefficients, aperture_radius, n_samples=100)` -- check
  dz/drho <= 0 on (0, a).

## `abo.acoustics`

### BEM solver -- `abo.acoustics.bem_solver`

- `solve_scattering(grid, frequency, incidence_direction, ...)` --
  solve the Burton-Miller CFIE for a rigid scatterer; returns the
  total field on the boundary as a `bempp_cl.api.GridFunction`.
- `scattered_field(grid, total_field, evaluation_points, frequency,
  incidence_direction, ...)` -- evaluate the scattered pressure at
  arbitrary exterior points using the Kirchhoff-Helmholtz
  representation formula.

### Target strength -- `abo.acoustics.target_strength`

- `monostatic_target_strength(grid, frequencies, incidence_angles,
  ...)` -- compute TS in dB over a frequency and angle grid. For
  each pair, solves the BEM problem and evaluates backscatter at
  `far_field_range`.

### Objectives -- `abo.acoustics.objectives`

- `integrated_conspicuousness(ts_matrix, frequencies, angles,
  call_spectrum)` -- IC in dB given a pre-computed TS matrix and a
  call-spectrum weighting.
- `surface_area(grid)` -- total triangle-mesh area in m^2.
- `compute_ic_for_grid(grid, frequencies, angles, call_spectrum, ...)`
  -- end-to-end convenience that runs the BEM pipeline and returns IC.

## `abo.biology`

### Call spectra -- `abo.biology.call_spectra`

- `glossophaga_call_spectrum(frequencies)` -- two-harmonic FM model
  based on Simon et al. (2006).
- `leptonycteris_call_spectrum(frequencies)` -- single-sweep model
  based on Gonzalez-Terrazas et al. (2016).
- `flat_band_spectrum(frequencies, f_min, f_max)` -- neutral uniform
  weighting for sensitivity analyses.

All spectra are normalised so that the integral over the input
frequency range equals 1.

### Natural shapes -- `abo.biology.natural_shapes`

- Constants: `MARCGRAVIA_EVENIA`, `MUCUNA_HOLTONII` -- published
  reflector dimensions.
- `SIMON_2020_CONFIGS` -- list of spherical cap configurations from
  Simon et al. (2020); geometrically impossible combinations excluded.

### Detection model -- `abo.biology.detection_model`

- `atmospheric_attenuation(frequency, distance, temperature, humidity)`
  -- placeholder (not yet implemented).
- `detection_threshold_db()` -- returns -65 dB by default.

## `abo.optimisation`

### Single-objective -- `abo.optimisation.single_objective`

- `optimise_cmaes(objective, x0, sigma0, bounds=None, ...)` --
  CMA-ES ask-and-tell driver. Returns `OptimisationResult` with
  best parameters, cost, evaluation count, and convergence history.

### Multi-objective -- `abo.optimisation.multi_objective`

- `pareto_frontier(objectives, bounds, population_size=50,
  n_generations=100, seed=None)` -- NSGA-II wrapper over a list
  of elementwise scalar objectives. Returns `ParetoResult` with
  `pareto_set`, `pareto_front`, `n_evals`.

### Encoding -- `abo.optimisation.encoding`

- `SphericalCapEncoding` -- (radius, half_angle); 2 parameters.
- `EllipsoidalCapEncoding` -- (semi_major, semi_minor, half_angle);
  3 parameters.
- `ChebyshevProfileEncoding(n_coefficients, aperture_radius)` -- N
  Chebyshev coefficients at fixed aperture.

Each encoder provides `dim`, `default_bounds`, `x0`, and `decode(x)`.

### Runners -- `abo.optimisation.runners`

- `EvaluationConfig` -- dataclass bundling frequencies, angles,
  call spectrum, and BEM parameters.
- `make_ic_cost(encoder, config, area_max=None, ...)` -- build a
  scalar cost (negative IC plus area penalty) consumable by CMA-ES.
- `make_pareto_objectives(encoder, config, ...)` -- build a
  `[neg_ic, sa]` pair consumable by NSGA-II.
- `build_encoder(family, **kwargs)` -- factory for the three
  built-in encoder families.

### Constraints -- `abo.optimisation.constraints`

- `concavity_penalty(coefficients, aperture_radius, n_samples=100)` --
  non-negative penalty; zero iff the profile is concave.
- `area_constraint(grid, max_area)` -- non-negative excess area;
  zero iff `surface_area(grid) <= max_area`.

## `abo.visualisation`

### Spectral plots -- `abo.visualisation.spectral_plots`

- `plot_ts_heatmap(ts_matrix, frequencies, angles, ...)` --
  frequency-angle heatmap.
- `plot_ts_polar(ts_at_freq, angles, frequency, ...)` -- polar
  directivity at one frequency.
- `plot_on_axis_vs_frequency(ts_matrix, frequencies, angles, ...)` --
  on-axis TS as a function of frequency.

### Pareto plots -- `abo.visualisation.pareto_plots`

Placeholders; to be filled in during production Phase 4 runs.

### Mesh plots -- `abo.visualisation.mesh_plots`

Placeholders; intended to wrap PyVista for 3D rendering.

## CLI -- `abo`

Commands (see `abo --help`):

- `abo solve --radius R --depth D [...]` -- target strength for a
  spherical cap.
- `abo optimise --family F --call C [...]` -- CMA-ES run.
- `abo pareto --family F --call C [...]` -- NSGA-II Pareto frontier.

See the `notebooks/` directory for worked examples of the non-CLI
API.
