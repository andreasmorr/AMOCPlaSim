# AMOCPlaSim

Edge-state and equilibrium analysis of the AMOC in the general circulation model PlaSim, using Gaussian covariance ellipsoids in EOF space as state-space boundaries.

Bisection trajectories from the saddle (edge) state to either the AMOC-on or AMOC-off attractor are available at pre-industrial (285 ppm) and current (360 ppm) CO₂ levels. Long equilibrium runs at all three states (on, off, edge) are used to fit Gaussian distributions in EOF space and to compute local stability metrics.

## File structure

```
AMOCPlaSim/
├── scripts/
│   └── plasim_edge_analysis.jl        # Main analysis script
├── src/
│   └── plasim_utils.jl                # NetCDF loading, EOF projection, ellipsoid fitting
├── data/
│   └── plasim/
│       ├── resilience_metrics.csv     # Key metrics for all states and CO2 levels
│       └── resilience_summaries.jld2  # Full cached results
├── Project.toml
└── Manifest.toml
```

## Analysis

### `plasim_edge_analysis.jl`
For each CO₂ level (285 ppm, 360 ppm) the script:
1. Loads all NetCDF edge-track files and classifies trajectories as converging to AMOC-on or AMOC-off.
2. Loads attractor and edge-state positions from converged equilibrium files.
3. Fits Gaussian covariance ellipsoids (at a chosen nσ level, default 3σ) to each state in EOF space.
4. Computes two primary resilience metrics:
   - **Convergence time**: transit time from last visit inside the edge ellipse to first entry into the target attractor ellipse (in the EOF1–EOF2 plane).
   - **Edge-to-attractor distance**: gap between the surfaces of the edge and attractor ellipsoids (zero if overlapping).
5. Computes local stability metrics from equilibrium runs: variance, dominant variance, lag-1 autocorrelation, integrated autocorrelation time per EOF, and mean AMOC strength from NetCDF output.
6. Produces scatter plots in EOF space, bar charts of convergence times and edge distances, and saves all key metrics to `data/plasim/resilience_metrics.csv`.

### `plasim_utils.jl`
Utility functions for loading and pre-processing PlaSim NetCDF files, projecting fields onto EOFs, and fitting/evaluating Gaussian ellipsoids.

## Usage

Edit the `DATA_DIR` and `N_FILES_*` constants at the top of `plasim_edge_analysis.jl` to point to your NetCDF data, then run from the project root:

```bash
julia --project scripts/plasim_edge_analysis.jl
```

Results are cached in `data/plasim/` and figures are written to `plots/`.

## Dependencies

Julia with [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/), [NCDatasets.jl](https://alexander-barth.github.io/NCDatasets.jl/stable/), [DataFrames.jl](https://dataframes.juliadata.org/stable/), and [CairoMakie](https://makie.org/). Install with:

```julia
using Pkg; Pkg.instantiate()
```
