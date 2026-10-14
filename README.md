# AMOCPlaSim

Edge-state and equilibrium analysis of the AMOC in PlaSim-LSG using box-mean
salinity coordinates. The module does not run new PlaSim experiments. It reduces
existing PlaSim-LSG edge-track and equilibrium salinity output to two
paper-facing salinity coordinates, then computes convergence and covariance
metrics from those reduced trajectories.

Bisection trajectories from the saddle state toward either the AMOC-on or
AMOC-off attractor are available at 285 ppm and 360 ppm CO2. Long equilibrium
runs at the AMOC-on, AMOC-off, and edge states provide the state clouds used for
Gaussian ellipses, convergence times, edge-to-attractor distances, local
variability, and the PlaSim rows of the synthesis figure.

## File Structure

```text
AMOCPlaSim/
├── scripts/
│   ├── compute_box_salinity.py        # Build box-mean salinity NetCDFs
│   ├── extract_amoc_strength.py       # Build AMOC-strength sidecar from Zenodo data
│   ├── plasim_edge_analysis.jl        # Compute resilience metrics and diagnostics
│   ├── plasim_export_paper_data.jl    # Export CSVs for paper plotting
│   ├── plotting_paper.py              # Paper figure from exported CSVs
│   └── regions.py                     # Vendored PlaRegion basin masks
├── src/
│   └── plasim_utils.jl                # NetCDF loading, classification, ellipse metrics
├── data/
│   ├── results/
│   │   ├── resilience_metrics.csv
│   │   ├── resilience_summaries.jld2
│   │   └── paper/
│   └── custom_readouts/
│       ├── amoc_strength_timeseries.nc
│       └── plasimelancholia_*.etc.nc  # Generated deep-box salinity NetCDFs
└── plots/
    ├── plasim_paper.png
    └── plasim_*.png
```

## Current Reduction

The paper-facing reduction uses two deep-box salinity coordinates:

| variable | box | depth range |
|----------|-----|-------------|
| `salt_na` | North Atlantic, 35N-80N | 0-1000 m |
| `salt_trop` | Tropical Atlantic, 35S-35N | 0-500 m |

The edge-track source salinity fields are Atlantic zonal means, so both boxes
are latitude-depth boxes (latitude centres 36.25-78.75N and 33.75S-33.75N on
the 2.5-degree grid). The analysis works in the two-dimensional plane
`(salt_na, salt_trop)`.

## Scripts

### `compute_box_salinity.py`

Builds `data/custom_readouts/` from the external PlaSim-LSG salinity data set.
Edge tracks are read from zonally averaged files
`<co2>/{lower,upper}/plasimedge_<co2>_edgetrack_itxNNN_<branch>.s_zonav.nc`.
For each shared index, the `upper` branch (AMOC-on) and `lower` branch
(AMOC-off) are combined into one output file with a `track` dimension.
Output files are numbered contiguously from 000. At 285 ppm the source indices
itx014 and itx025 do not exist in either branch, so output `itx014` is source
`itx015`, and so on; each file stores its source index in the
`source_itx_index` attribute.

Equilibrium files are read from `<co2>/plasimedge_<co2>_{on,of,ed}.tsv.nc`.
They contain full salinity fields and are first reduced to the Atlantic zonal
mean with the vendored `PlaRegion` mask before the same box averages are
computed.

The source root defaults to `/Volumes/KINGSTON/BoernerEtAl` (override with
`--src`). The grid files for the mask (`wet.nc`, `basins_scalar.nc`,
`basins_vector.nc`) are read from the hard-coded `GRID_DIR`. On the data drive
they are in `BoernerEtAl/part1/grid`, but `GRID_DIR` currently points to
`/Volumes/KINGSTON/part1/grid`, so adjust it before running the equilibrium
step.

### `extract_amoc_strength.py`

Builds `data/custom_readouts/amoc_strength_timeseries.nc`, a compact sidecar
with the AMOC-strength time series (streamfunction maximum, 46-66N, Sv) used
for the top row of the paper figure and the mean-AMOC-strength column of the
metrics CSV. It reads the public Zenodo Part I data (default root
`/Volumes/KINGSTON/BoernerEtAl/part1`, override with `--src`):

- equilibrium runs from `<co2>/plasimedge_<co2>_{on,of,ed}.lsg.nc`,
- edge tracks from the `.etc.nc` members of
  `<co2>/plasimedge_<co2>_ed-iterations.zip` (read directly from the archive).

Variables are keyed like the matching salinity readout without `.etc.nc`
(e.g. `plasimelancholia_285ppm_edgetrack_itx000`). The archive members use the
same contiguous numbering as the salinity readouts. When the readouts exist,
the script checks that every AMOC series has the same number of valid years as
the matching salinity track and stops otherwise.

### `plasim_edge_analysis.jl`

Computes the PlaSim resilience diagnostics for 285 ppm and 360 ppm:

- trajectory labels: AMOC-on versus AMOC-off, from the mean North Atlantic
  salinity over the last 10% of each trajectory, split at the median over all
  trajectories of that CO2 level (so exactly half are labelled AMOC-on),
- convergence time: from the last visit inside the edge-state 4-sigma ellipse
  to the first later entry into the target attractor's 4-sigma ellipse;
  trajectories that never get there are excluded from the mean,
- edge-to-attractor distance: gap between the edge and attractor 4-sigma
  ellipses, measured along the line joining their centres,
- local variability and local resilience from the AMOC-on and AMOC-off
  equilibrium covariance matrices,
- metrics CSV, `resilience_summaries.jld2`, and five diagnostic plots.

The start of each edge-state equilibrium run is trimmed before computing the
edge mean/covariance: 4 years at 285 ppm and 60 years at 360 ppm.

### `plasim_export_paper_data.jl`

Exports the current paper data to `data/results/paper/`:

- `trajectories_{285ppm,360ppm}.csv`
- `equilibria_{285ppm,360ppm}.csv`
- `ellipses_{285ppm,360ppm}.csv`
- `state_means_{285ppm,360ppm}.csv`

`trajectories_*.csv` contains only the trajectories that converged in the
sense above.

### `plotting_paper.py`

Reads the exported CSVs and writes `plots/plasim_paper.png`.

## Usage

Run from the AMOCPlaSim directory:

```bash
python scripts/compute_box_salinity.py
python scripts/extract_amoc_strength.py
julia --project scripts/plasim_edge_analysis.jl
julia --project scripts/plasim_export_paper_data.jl
python scripts/plotting_paper.py
```

All analysis and plotting scripts use the two deep-box variables directly; no
alternate shallow or EOF mode is retained.

## Dependencies

Julia dependencies are listed in `Project.toml`; instantiate with:

```julia
using Pkg; Pkg.instantiate()
```

Both Julia analysis scripts also load `CSV`, which is not yet listed in
`Project.toml`; add it with `Pkg.add("CSV")` unless it is in your default
environment.

Python requires `numpy`, `matplotlib`, `pandas`, `xarray`, and `netCDF4`.
`plotting_paper.py` also imports `../amoc_plot_style.py` from the parent
AMOCResilience repository.
`compute_box_salinity.py` also needs the PlaSim grid files referenced by
`GRID_DIR` for the vendored `scripts/regions.py` mask helper.
