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
│   ├── extract_amoc_strength.jl       # Build compact AMOC-strength sidecar
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

The source salinity fields are Atlantic zonal means, so the Southern Ocean box
has no data coverage in this reduction. The effective state space used by the
analysis is therefore two-dimensional.

## Scripts

### `compute_box_salinity.py`

Builds `data/custom_readouts/` from the external PlaSim-LSG salinity data set.
Edge tracks are read from zonally averaged files
`<co2>/{lower,upper}/plasimedge_<co2>_edgetrack_itxNNN_<branch>.s_zonav.nc`.
For each shared index, the `upper` branch (AMOC-on) and `lower` branch
(AMOC-off) are combined into one output file with a `track` dimension.

Equilibrium files are read from `<co2>/plasimedge_<co2>_{on,of,ed}.tsv.nc`.
They contain full salinity fields and are first reduced to the Atlantic zonal
mean with the vendored `PlaRegion` mask before the same box averages are
computed.

### `extract_amoc_strength.jl`

Builds `data/custom_readouts/amoc_strength_timeseries.nc`, a compact sidecar containing
only the AMOC-strength time series needed for the top row of the paper figure
and for the mean-AMOC-strength column of the metrics CSV. Variables in the
sidecar are keyed by the original source filename without `.etc.nc`. This script
is only needed when regenerating the sidecar from external legacy
reduced-coordinate NetCDF files.

### `plasim_edge_analysis.jl`

Computes the PlaSim resilience diagnostics for 285 ppm and 360 ppm:

- trajectory labels: AMOC-on versus AMOC-off, classified by final North
  Atlantic salinity,
- convergence time: from the last visit inside the edge ellipse to first entry
  into the target attractor ellipse,
- edge-to-attractor distance: distance between the edge and attractor Gaussian
  ellipses,
- local variability and local resilience from the AMOC-on equilibrium
  covariance ellipse,
- metrics CSV and diagnostic plots.

The start of each edge-state equilibrium run is trimmed before computing the
edge mean/covariance: 4 years at 285 ppm and 60 years at 360 ppm.

### `plasim_export_paper_data.jl`

Exports the current paper data to `data/results/paper/`:

- `trajectories_{285ppm,360ppm}.csv`
- `equilibria_{285ppm,360ppm}.csv`
- `ellipses_{285ppm,360ppm}.csv`
- `state_means_{285ppm,360ppm}.csv`

### `plotting_paper.py`

Reads the exported CSVs and writes `plots/plasim_paper.png`.

## Usage

Run from the AMOCPlaSim directory:

```bash
python scripts/compute_box_salinity.py
julia --project scripts/plasim_edge_analysis.jl
julia --project scripts/plasim_export_paper_data.jl
python scripts/plotting_paper.py
```

If `data/custom_readouts/amoc_strength_timeseries.nc` needs to be regenerated,
place the legacy reduced-coordinate `.etc.nc` source files in `data/results/`
and run:

```bash
julia --project scripts/extract_amoc_strength.jl
```

All analysis and plotting scripts use the two deep-box variables directly; no
alternate shallow or EOF mode is retained.

## Dependencies

Julia dependencies are listed in `Project.toml`; instantiate with:

```julia
using Pkg; Pkg.instantiate()
```

Python requires `numpy`, `matplotlib`, `pandas`, `xarray`, and `netCDF4`.
`compute_box_salinity.py` also needs the PlaSim grid files referenced by
`GRID_DIR` for the vendored `scripts/regions.py` mask helper.
