# AMOCPlaSim

Edge-state and equilibrium analysis of the AMOC in the general circulation model PlaSim, using Gaussian covariance ellipsoids in EOF space as state-space boundaries.

Bisection trajectories from the saddle (edge) state to either the AMOC-on or AMOC-off attractor are available at pre-industrial (285 ppm) and present-day (360 ppm) CO₂ levels. Long equilibrium runs at all three states (on, off, edge) are used to fit Gaussian distributions in EOF space and to compute local stability metrics.

---

## File structure

```
AMOCPlaSim/
├── scripts/
│   ├── plasim_edge_analysis.jl        # Main analysis script
│   ├── plasim_export_paper_data.jl    # Export CSVs for paper figures
│   └── plotting_paper.py              # Paper figure (reads exported CSVs)
├── src/
│   └── plasim_utils.jl                # NetCDF loading, EOF projection, ellipsoid fitting
├── data/
│   └── plasim/
│       ├── resilience_metrics.csv     # Key metrics for all states and CO₂ levels
│       ├── resilience_summaries.jld2  # Full cached results
│       └── paper/                     # CSV exports for plotting
│           ├── trajectories_{285,360}ppm.csv  # Filtered converged trajectory time series
│           ├── equilibria_{285,360}ppm.csv    # Equilibrium run time series (on/off/edge)
│           ├── ellipses_{285,360}ppm.csv      # Gaussian ellipse boundary coordinates
│           └── state_means_{285,360}ppm.csv   # Mean EOF positions for each state
├── plots/
│   └── plasim_paper.pdf               # Output paper figure
├── Project.toml
└── Manifest.toml
```

---

## Scripts

### `plasim_edge_analysis.jl`

For each CO₂ level (285 ppm, 360 ppm) the script:
1. Loads all NetCDF edge-track files and classifies trajectories as converging to AMOC-on or AMOC-off.
2. Loads attractor and edge-state positions from converged equilibrium files.
3. Fits Gaussian covariance ellipsoids (at a chosen nσ level, default 3σ) to each state in EOF space.
4. Computes two primary resilience metrics:
   - **Convergence time**: transit time from last visit inside the edge ellipse to first entry into the target attractor ellipse (in the EOF1–EOF2 plane).
   - **Edge-to-attractor distance**: gap between the surfaces of the edge and attractor ellipsoids (zero if overlapping).
5. Computes local stability metrics from equilibrium runs: variance, dominant variance, lag-1 autocorrelation, integrated autocorrelation time per EOF, mean AMOC strength, and the **1σ semi-major axis of the Gaussian ellipse** fitted in the (EOF1, EOF2) plane (`ellipse_long_axis_1sigma` = `sqrt(max eigenvalue of C[1:2,1:2])`). This quantity is used as a stand-in for characteristic return time in the synthesis figure.
6. Saves all key metrics to `data/plasim/resilience_metrics.csv`.

### `plasim_export_paper_data.jl`

Exports all data required by `plotting_paper.py` to `data/plasim/paper/`. For each CO₂ level: filtered converged trajectories in EOF space, equilibrium run time series (on/off/edge states), Gaussian ellipse boundary coordinates, and state mean positions. Run this once after `plasim_edge_analysis.jl` has cached results.

### `src/plasim_utils.jl`

Utility functions for loading and pre-processing PlaSim NetCDF files, projecting fields onto EOFs, and fitting and evaluating Gaussian ellipsoids.

---

## Paper figure

`plotting_paper.py` produces a publication-quality 4-panel figure using the shared design language from `../amoc_plot_style.py`.

**Figure layout:**
- **Top row** (shorter): AMOC strength vs time for all converged trajectories at 285 ppm (left) and 360 ppm (right).
- **Bottom row** (square): 2D EOF phase portrait. Background shows faint equilibrium run lines for on/off/edge states. Gaussian covariance ellipsoids are overlaid and annotated. Trajectories are time-shaded (alpha increases with time). Attractor mean positions marked with stars.

Output: `plots/plasim_paper.pdf`

---

## Usage

Edit the `DATA_DIR` and `N_FILES_*` constants at the top of `plasim_edge_analysis.jl` to point to your NetCDF data, then run from the project root:

```bash
# 1. Run main analysis (fits ellipsoids, computes metrics)
julia --project scripts/plasim_edge_analysis.jl

# 2. Export paper data (trajectories, equilibria, ellipses)
julia --project scripts/plasim_export_paper_data.jl

# 3. Generate paper figure
python scripts/plotting_paper.py
```

Results are cached in `data/plasim/` and figures are written to `plots/`.

---

## Dependencies

Julia with [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/), [NCDatasets.jl](https://alexander-barth.github.io/NCDatasets.jl/stable/), [DataFrames.jl](https://dataframes.juliadata.org/stable/), and [CairoMakie](https://makie.org/). Python requires `numpy` and `matplotlib`. Install Julia packages with:

```julia
using Pkg; Pkg.instantiate()
```
