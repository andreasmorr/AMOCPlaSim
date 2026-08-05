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
│   ├── plotting_paper.py              # Paper figure (reads exported CSVs)
│   ├── compute_box_salinity.py        # Box-mean salinity time series (Wood/CLIMBER-X boxes)
│   └── regions.py                     # Börner's PlaRegion basin masks (vendored)
├── src/
│   └── plasim_utils.jl                # NetCDF loading, EOF projection, ellipsoid fitting
├── data/
│   ├── plasim/
│   │   ├── resilience_metrics.csv     # Key metrics for all states and CO₂ levels
│   │   ├── resilience_summaries.jld2  # Full cached results
│   │   └── paper/                     # CSV exports for plotting
│   │       ├── trajectories_{285,360}ppm.csv  # Filtered converged trajectory time series
│   │       ├── equilibria_{285,360}ppm.csv    # Equilibrium run time series (on/off/edge)
│   │       ├── ellipses_{285,360}ppm.csv      # Gaussian ellipse boundary coordinates
│   │       └── state_means_{285,360}ppm.csv   # Mean EOF positions for each state
│   └── plasim_boxsalt/               # Box-mean salinity files (compute_box_salinity.py output)
│       ├── plasimelancholia_{285ppm,360ppm}_edgetrack_{itx,iter}NNN.etc.nc  # Edge tracks (track dim)
│       └── plasimelancholia_{285ppm,360ppm}_{on,of,ed}.etc.nc               # Equilibria
├── plots/
│   └── plasim_paper.png               # Output paper figure (200 dpi PNG)
├── Project.toml
└── Manifest.toml
```

---

## Scripts

### `plasim_edge_analysis.jl`

For each CO₂ level (285 ppm, 360 ppm) the script:
1. Loads all NetCDF edge-track files and classifies trajectories as converging to AMOC-on or AMOC-off.
2. Loads attractor and edge-state positions from converged equilibrium files.
3. Fits Gaussian covariance ellipsoids (at a chosen nσ level, `ELLIPSE_SIGMA`, currently 4σ) to each state in EOF space. The start of each edge (saddle) equilibrium run is trimmed before its mean/covariance/ellipse are computed (transient spin-up; 285 ppm −4 yr, 360 ppm −60 yr, set by `EDGE_CUT` in the analysis and export scripts — not in the data preprocessing).
4. Computes two primary resilience metrics:
   - **Convergence time**: transit time from last visit inside the edge ellipse to first entry into the target attractor ellipse (in the EOF1–EOF2 plane).
   - **Edge-to-attractor distance**: gap between the surfaces of the edge and attractor ellipsoids (zero if overlapping).
5. Computes local stability metrics from equilibrium runs: variance, dominant variance, lag-1 autocorrelation, integrated autocorrelation time per EOF, mean AMOC strength, and **local resilience** as the inverse of the long axis of the 1σ Gaussian ellipse in the (EOF1, EOF2) plane: `local_resilience = 1 / ellipse_long_axis_1sigma = 1 / (2000 × sqrt(λ_max(C[1:2,1:2])))`, where `λ_max` is the largest eigenvalue of the 2×2 marginal covariance. A smaller (less elongated) attractor ellipse implies higher local resilience. This quantity is the stand-in for local resilience in the synthesis figure.
6. Saves all key metrics to `data/plasim/resilience_metrics.csv`.

### `plasim_export_paper_data.jl`

Exports all data required by `plotting_paper.py` to `data/plasim/paper/`. For each CO₂ level: filtered converged trajectories in EOF space, equilibrium run time series (on/off/edge states), Gaussian ellipse boundary coordinates, and state mean positions. Run this once after `plasim_edge_analysis.jl` has cached results.

### `src/plasim_utils.jl`

Utility functions for loading and pre-processing PlaSim NetCDF files, projecting fields onto EOFs, and fitting and evaluating Gaussian ellipsoids.

### `compute_box_salinity.py`

Builds an alternative, box-based state-space reduction from the raw PlaSim-LSG salinity fields (Börner et al. 2025), replacing the three EOF-reduced coordinates (`redu1/2/3`) with three **box-mean salinity** time series. The boxes are the Wood/CLIMBER-X perturbation regions re-created on the PlaSim grid (see `../plotting_perturbations.py`), each averaged over the top 100 m (top 2 depth levels), with no tapering:

| variable | box | latitude band |
|----------|-----|---------------|
| `salt_na`    | North Atlantic | 35°N – 80°N |
| `salt_trop`  | Tropical       | 35°S – 35°N |
| `salt_south` | Southern       | 90°S – 35°S |

Each box mean is a volume-weighted average over its (lat × depth) cells, weighting by `cos(lat) × layer_thickness` and skipping land cells (NaN).

**Inputs** (external, on the data drive — see the constants at the top of the script):
- Edge tracks: `<co2>/{lower,upper}/plasimedge_<co2>_edgetrack_itxNNN_<branch>.s_zonav.nc` — already zonally-averaged salinity `s(time, depth, lat)`. For each index the `lower` (→ AMOC-off) and `upper` (→ AMOC-on) branches are combined into one output file with a `track = ['upper', 'lower']` dimension; unequal track lengths are padded with NaN.
- Equilibria: `<co2>/plasimedge_<co2>_{on,of,ed}.tsv.nc` — full 3D salinity `s(time, depth, lat, lon)`. These are first reduced to the Atlantic zonal mean using Börner's `PlaRegion` mask (`Atlantic3D`, southern border −34°, i.e. `s.where(Atlantic3D()).mean(dim='lon')`), then averaged into the same boxes. Written as single-trajectory `(time,)` files. `PlaRegion` lives in the vendored `scripts/regions.py` and reads its grid files (`wet.nc`, `basins_scalar.nc`, `basins_vector.nc`) from the `GRID_DIR` path set at the top of the script.

**Output** → `data/plasim_boxsalt/`, formatted to match the example `.etc.nc` files in `data/plasim/` so the existing loaders in `src/plasim_utils.jl` can read them directly. Edge-track files are re-indexed contiguously with the existing `itx` (285 ppm) / `iter` (360 ppm) naming; global/variable attributes record the source file and box geometry.

> **Note:** the source is an **Atlantic-only** zonal mean (≈ 35°S – 80°N), so the global Southern Ocean box has no data and `salt_south` is all-NaN in every file — effectively a two-box (NA, Trop) reduction with an explicit no-coverage marker for the third box.

To run the existing Julia analysis on these files, point `DATA_DIR` at `data/plasim_boxsalt`, set `VARIABLE_NAMES = ["salt_na", "salt_trop", "salt_south"]` in `plasim_edge_analysis.jl`, and note that the AMOC-on branch is now the *saltier* NA box (so the on/off split may need its sign flipped).

Run it with:

```bash
python scripts/compute_box_salinity.py                 # edge tracks + equilibria
python scripts/compute_box_salinity.py --skip-equilibria   # edge-track files only
python scripts/compute_box_salinity.py --skip-edgetrack    # equilibrium files only
```

---

## Paper figure

`plotting_paper.py` produces a publication-quality 4-panel figure using the shared design language from `../amoc_plot_style.py`.

**Figure layout:**
- **Top row** (shorter): AMOC strength vs time for all converged on-state and off-state edge-track trajectories at 285 ppm (left) and 360 ppm (right). AMOC-on and AMOC-off equilibria shown as dashed horizontal lines. If the `amoc_strength` column is absent from the trajectory CSV, a placeholder message is shown.
- **Bottom row** (square): 2D EOF phase portrait. Gaussian covariance ellipsoids (thick dashed lines) outline each state (on, off, edge). Converged trajectories are plotted at fixed opacity.

Output: `plots/plasim_paper.png`

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

### Choosing the state-space reduction (EOF vs. box salinity)

The three analysis scripts can run on any of three reduced state-space representations, selected by a single `MODE` constant near the top of each:

| `MODE` | Coordinates | Data dir | Dims |
|--------|-------------|----------|------|
| `:boxsalt_deep` / `"boxsalt_deep"` (default) | `salt_na_deep` (full column), `salt_trop_deep` (0–500 m) | `data/plasim_boxsalt` | 2 |
| `:boxsalt` / `"boxsalt"` | `salt_na`, `salt_trop` (0–100 m box-mean salinities) | `data/plasim_boxsalt` | 2 |
| `:eof` / `"eof"` | `redu1`, `redu2`, `redu3` (Börner EOFs) | `data/plasim` | 3 |

Set `MODE` **consistently** in `plasim_edge_analysis.jl`, `plasim_export_paper_data.jl`, and `plotting_paper.py`, then run the three steps as above. In the box-salinity modes:

- The box files are produced by `compute_box_salinity.py` (run that first); each file carries both the shallow and deep variables.
- Each is a **2-D** reduction — the Southern box has no Atlantic data — so the 3-D scatter figures are skipped and the AMOC-on branch is the *saltier* North Atlantic box (handled via `on_is_low=false` in `classify_trajectories`). Note the deep NA box (full column) barely distinguishes AMOC-on from AMOC-off, since it is dominated by the stable deep ocean.
- Outputs are written with a mode suffix so they never overwrite each other or the EOF results: `resilience_metrics_<suffix>.csv`, `resilience_summaries_<suffix>.jld2`, `plots/plasim_*_<suffix>.png`, and paper CSVs in `data/plasim/paper_<suffix>/`, where `<suffix>` is `boxsalt` or `boxsalt_deep`.
- The box files carry no `amoc_strength`, so it is read from the EOF equilibrium/edge-track files (`AMOC_DIR = data/plasim`) for the paper-figure top row and the metrics CSV.

The umbrella `synthesis_figure.py` reads the deep box-salinity metrics (`data/plasim/resilience_metrics_boxsalt_deep.csv`, matching the pipeline default); point `PLASIM_CSV` at `resilience_metrics_boxsalt.csv` (shallow) or `resilience_metrics.csv` (EOF) for the other reductions.

---

## Dependencies

Julia with [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/), [NCDatasets.jl](https://alexander-barth.github.io/NCDatasets.jl/stable/), [DataFrames.jl](https://dataframes.juliadata.org/stable/), and [CairoMakie](https://makie.org/). Python requires `numpy` and `matplotlib` (plus `xarray` + `netCDF4` for `compute_box_salinity.py`; its Atlantic mask uses the vendored `scripts/regions.py`, which needs the grid files at `GRID_DIR`). Install Julia packages with:

```julia
using Pkg; Pkg.instantiate()
```
