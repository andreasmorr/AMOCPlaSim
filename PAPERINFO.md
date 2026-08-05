# PAPERINFO: AMOCPlaSim

## Purpose

This file is the paper-facing summary for the AMOCPlaSim submodule. Use it as
the local source of truth when updating `../main.tex`, especially the Methods,
Appendix, synthesis table, and figure captions. It intentionally omits general
usage instructions and keeps only information needed to describe the scientific
data reduction and resilience analysis.

If this file conflicts with `../main.tex`, prefer this file for the AMOCPlaSim
analysis setup.

## Module Role In The Paper

AMOCPlaSim provides the PlaSim-LSG member of the AMOC resilience hierarchy. The
submodule does not run new PlaSim experiments. It reduces existing PlaSim-LSG
edge-track and equilibrium salinity output to two custom deep-box salinity
coordinates, then computes convergence, edge-distance, AMOC-strength, and local
variability metrics from those reduced trajectories.

The submodule contributes:

- the PlaSim rows of the cross-model resilience synthesis,
- a PlaSim paper figure showing AMOC-strength time series and salinity-space
  edge-track trajectories,
- diagnostics for convergence times, edge-to-attractor distances, and local
  covariance-based resilience.

Do not describe the current PlaSim analysis as EOF-based. The old EOF time
series and EOF-analysis outputs have been removed.

## Source Files

The implemented analysis is defined by:

- `scripts/compute_box_salinity.py`
- `scripts/extract_amoc_strength.jl`
- `scripts/plasim_edge_analysis.jl`
- `scripts/plasim_export_paper_data.jl`
- `scripts/plotting_paper.py`
- `src/plasim_utils.jl`
- `scripts/regions.py`

Primary paper-facing outputs are in:

- `data/custom_readouts/plasimelancholia_*.etc.nc`
- `data/custom_readouts/amoc_strength_timeseries.nc`
- `data/results/resilience_metrics.csv`
- `data/results/resilience_summaries.jld2`
- `data/results/paper/`
- `plots/plasim_paper.png`

The umbrella synthesis figure reads the PlaSim metrics from
`AMOCPlaSim/data/results/resilience_metrics.csv`.

## External PlaSim Data

The raw PlaSim-LSG data are external to this repository. The source data set is
identified in the generated NetCDF metadata as:

```text
https://doi.org/10.5281/zenodo.17053348
```

Edge-track raw files are expected under external `285ppm/` and `360ppm/`
folders:

```text
<co2>/{lower,upper}/plasimedge_<co2>_edgetrack_itxNNN_<branch>.s_zonav.nc
```

Equilibrium raw files are expected as full 3D salinity fields:

```text
<co2>/plasimedge_<co2>_{on,of,ed}.tsv.nc
```

## Custom Salinity Readouts

The tracked custom readouts in `data/custom_readouts/` contain only two
paper-facing salinity variables:

| variable | box | depth range |
|----------|-----|-------------|
| `salt_na` | North Atlantic, 35N-80N | 0-1000 m |
| `salt_trop` | Tropical Atlantic, 35S-35N | 0-500 m |

The names `salt_na` and `salt_trop` now refer to these deep boxes. There is no
retained shallow-box mode and no retained `_deep` variable suffix.

The readouts are volume-weighted box means. The weight is
`cos(latitude) * layer_thickness`, and land cells represented by NaNs are
excluded from both the weighted sum and the normalization. Depth selection uses
PlaSim layer centers:

- `salt_na`: top 13 layers, deepest included center at 950 m, bottom face at
  1025 m.
- `salt_trop`: top 8 layers, deepest included center at 450 m, bottom face at
  500 m.

Because the edge-track salinity source fields are already Atlantic zonal means,
the North Atlantic and Tropical Atlantic boxes reduce to latitude-depth boxes.
For the equilibrium runs, `scripts/compute_box_salinity.py` first applies the
vendored `PlaRegion` Atlantic mask (`Atlantic3D`, southern border -34 degrees),
takes the Atlantic zonal mean, and then applies the same box averaging.

## Retained Files

The generated salinity readouts consist of:

- 37 paired 285 ppm edge-track files:
  `plasimelancholia_285ppm_edgetrack_itx000.etc.nc` through
  `plasimelancholia_285ppm_edgetrack_itx036.etc.nc`.
- 38 paired 360 ppm edge-track files:
  `plasimelancholia_360ppm_edgetrack_iter000.etc.nc` through
  `plasimelancholia_360ppm_edgetrack_iter037.etc.nc`.
- Six equilibrium files:
  `plasimelancholia_{285ppm,360ppm}_{on,of,ed}.etc.nc`.

Each edge-track file contains two tracks with coordinate labels `upper` and
`lower`. The source convention is:

- `upper`: branch converging to AMOC-on,
- `lower`: branch converging to AMOC-off.

The analysis still classifies trajectories from their final North Atlantic
salinity rather than from the branch label. In the retained coordinates,
AMOC-on has higher final `salt_na`, so `ON_IS_LOW = false`.

## AMOC Strength

The generated salinity readouts do not carry AMOC strength. AMOC strength is
stored separately in:

```text
data/custom_readouts/amoc_strength_timeseries.nc
```

This sidecar was cut from the legacy reduced-coordinate PlaSim NetCDF files.
Each variable is keyed by the original source filename without the `.etc.nc`
suffix, for example:

```text
plasimelancholia_285ppm_on
plasimelancholia_285ppm_edgetrack_itx000
```

Equilibrium AMOC-strength variables are one-dimensional. Edge-track
AMOC-strength variables are two-dimensional over time and track. The sidecar is
used for the top row of `plots/plasim_paper.png` and for the
`mean_amoc_strength_Sv` column in `data/results/resilience_metrics.csv`.

## Resilience Analysis

The current analysis is two-dimensional in `(salt_na, salt_trop)`.

Constants in `scripts/plasim_edge_analysis.jl` and
`scripts/plasim_export_paper_data.jl`:

| setting | value |
|---------|-------|
| `CO2_LABEL_PREINDUSTRIAL` | `285ppm` |
| `CO2_LABEL_CURRENT` | `360ppm` |
| `N_FILES_285` | `37` |
| `N_FILES_360` | `38` |
| `VARIABLE_NAMES` | `["salt_na", "salt_trop"]` |
| `FINAL_FRACTION` | `0.1` |
| `ELLIPSE_SIGMA` | `4` |
| `ON_IS_LOW` | `false` |
| 285 ppm edge trim | first 4 years removed |
| 360 ppm edge trim | first 60 years removed |

Attractor positions are not estimated from the final states of the edge-track
trajectories. They are time means of the dedicated equilibrium files:

- AMOC-on attractor: `plasimelancholia_<co2>_on.etc.nc`
- AMOC-off attractor: `plasimelancholia_<co2>_of.etc.nc`
- edge state: `plasimelancholia_<co2>_ed.etc.nc`, after the edge-trim above

Local covariance matrices are computed from the dedicated AMOC-on and AMOC-off
equilibrium files. Edge covariance is computed from the trimmed edge
equilibrium file.

Convergence time is computed with the ellipse-based method in
`src/plasim_utils.jl`: for each converged trajectory, it is the time from the
last visit inside the edge-state ellipse to the first entry into the target
attractor ellipse. The script supplies the covariance matrices, so the legacy
epsilon-ball fallback is not the paper-facing method.

Edge-to-attractor distance is the gap between the edge-state Gaussian ellipse
and the target-attractor Gaussian ellipse in the retained two-dimensional
salinity plane. The paper's PlaSim critical-shock proxy is the AMOC-on row of
`mean_edge_dist`.

Local resilience is stored as:

```text
local_resilience = 1 / (2000 * sqrt(lambda_max(C[1:2, 1:2])))
```

where `C` is the local equilibrium covariance matrix in the retained salinity
coordinates. This is a covariance-ellipse proxy, not a linear eigenvalue of the
underlying PlaSim-LSG dynamics.

No basin volume is estimated for PlaSim in the current analysis.

## Current Metrics

Current `data/results/resilience_metrics.csv` values:

| CO2 ppm | state | mean convergence time yr | mean edge distance | local resilience | mean AMOC strength Sv |
|---------|-------|--------------------------|--------------------|------------------|-----------------------|
| 285 | AMOC-on | 105.25 | 0.0498135 | 0.0356814 | 14.7227 |
| 285 | AMOC-off | 116.00 | 0.2040272 | 0.0115056 | 1.1161 |
| 360 | AMOC-on | 150.6316 | 0.1101618 | 0.0368075 | 15.0263 |
| 360 | AMOC-off | 74.8182 | 0.0932569 | 0.0086736 | 1.6823 |

The paper CSV bundle in `data/results/paper/` contains:

- `trajectories_{285ppm,360ppm}.csv`
- `equilibria_{285ppm,360ppm}.csv`
- `ellipses_{285ppm,360ppm}.csv`
- `state_means_{285ppm,360ppm}.csv`

The current row counts, including headers, are:

| file group | 285 ppm | 360 ppm |
|------------|---------|---------|
| `trajectories` | 1859 | 17384 |
| `equilibria` | 4495 | 5344 |
| `ellipses` | 364 | 364 |
| `state_means` | 4 | 4 |

## Paper Wording Notes

Use "custom deep-box salinity readouts" or "deep-box salinity coordinates" for
the PlaSim state space. Avoid "EOF coordinates", "3D EOF time series",
"shallow-box mode", or "PlaSim experiments" for this submodule.

The 285 ppm case should be described as pre-industrial and the 360 ppm case as
current CO2. The current PlaSim contribution has only these two CO2 levels.
