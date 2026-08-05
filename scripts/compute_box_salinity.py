#!/usr/bin/env python3
"""
compute_box_salinity.py — box-mean salinity time series from PlaSim edge tracks.

Reads the zonally-averaged salinity fields ``s(time, depth, lat)`` produced by
Börner et al. (files ``plasimedge_<co2>_edgetrack_itx<NNN>_<branch>.s_zonav.nc``
in the ``285ppm`` / ``360ppm`` source folders) and reduces each field to box-mean
salinity time series using the box geometry re-created for the PlaSim grid in
``../plotting_perturbations.py`` (CLIMBER-X boxes, no tapering).  Only the
paper-facing deep boxes are written:

    salt_na     North Atlantic, 35N-80N , 0-1000 m
    salt_trop   Tropical Atlantic, 35S-35N , 0-500 m

500 m is the cell edge below the 450 m layer (top 8 levels). 1000 m selects
the top 13 levels with bottom cell edge 1025 m.

Because the input field is already zonally averaged, the CLIMBER-X Atlantic-basin
restriction (NA/Trop) reduces to a plain latitude band; each box mean is a
volume-weighted average over its (lat × depth) cells, weighting by
``cos(lat) × layer_thickness`` and skipping land cells (NaN).

For every trajectory index that exists in both the ``lower`` and ``upper``
branch folders, the two branches are combined into a single output file with a
``track`` dimension (``['upper', 'lower']`` — upper → AMOC-on, lower → AMOC-off),
matching the format expected by ``src/plasim_utils.jl`` so
the existing analysis pipeline (``src/plasim_utils.jl``) can read them directly.

The three converged-state equilibrium files per CO2 level
(``plasimedge_<co2>_<state>.tsv.nc``, state ∈ {on, of, ed}) sit directly in the
CO2 folders as full 3D salinity fields ``s(time, depth, lat, lon)``.  They are
first reduced to the Atlantic zonal mean using Börner's ``PlaRegion`` mask
(``Atlantic3D``, southern border -34° — the same Atlantic definition behind the
edge-track s_zonav files), i.e. ``s.where(Atlantic3D()).mean(dim='lon')``, and
then the identical box averaging is applied.  Each is written as a single
``(time,)`` file matching the analysis loader's expected naming convention.

Output: ``data/custom_readouts/plasimelancholia_<co2>_edgetrack_<tag><NNN>.etc.nc``
(``tag`` = ``itx`` for 285 ppm, ``iter`` for 360 ppm, contiguously indexed) and
``.../plasimelancholia_<co2>_<state>.etc.nc`` for the equilibria.

To run the existing Julia analysis on these files, use the default
``VARIABLE_NAMES = ["salt_na", "salt_trop"]`` in
``scripts/plasim_edge_analysis.jl``.

Usage:
    python scripts/compute_box_salinity.py [--src SRC_ROOT] [--co2 285ppm 360ppm]
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SRC_ROOT = Path("/Volumes/KINGSTON/BoernerEtAl")
OUT_ROOT = SCRIPT_DIR.parent / "data" / "custom_readouts"

# Box depths. The North Atlantic box reaches ~1000 m; the Tropical box reaches
# the 500 m cell edge (the bottom face of the 450 m layer = top 8 layers).
TROP_DEEP_DEPTH = 500.0     # top 8 layers (25..450 m); 500 m is an exact cell edge
NA_DEEP_DEPTH   = 1000.0    # top 13 layers (25..950 m); bottom cell edge 1025 m (~1000 m)
# depth_max = None  ->  full water column

# Box latitude bands (cell-centre test, identical to plotting_perturbations.py).
BOXES = {
    "salt_na":   dict(lat_min=35.0,  lat_max=80.0, depth_max=NA_DEEP_DEPTH,
                      long_name="North Atlantic box mean salinity (35-80N, 0-1000 m)"),
    "salt_trop": dict(lat_min=-35.0, lat_max=35.0, depth_max=TROP_DEEP_DEPTH,
                      long_name="Tropical box mean salinity (35S-35N, 0-500 m)"),
}


def _depth_attr(box: dict):
    """NetCDF attribute value for a box's depth extent."""
    return "full water column" if box["depth_max"] is None else box["depth_max"]

# Track order (matches the example .etc.nc files).
TRACKS = ["upper", "lower"]          # upper -> AMOC-on, lower -> AMOC-off

# Filename tag per CO2 level (matches the existing analysis naming).
CO2_TAG = {"285ppm": "itx", "360ppm": "iter"}

# Equilibrium (converged-state) files live directly in the CO2 folders as full
# 3D salinity fields ``plasimedge_<co2>_<state>.tsv.nc``.  They are reduced to
# the Atlantic zonal mean with Börner's PlaRegion mask (same Atlantic definition,
# southern border -34°, that produced the edge-track s_zonav files) before the
# identical box averaging is applied.
EQ_STATES = ["on", "of", "ed"]       # AMOC-on, AMOC-off, edge (saddle)
# regions.py (Börner's PlaRegion) is vendored next to this script; the grid files
# it needs (wet.nc, basins_scalar.nc, basins_vector.nc) live on the data drive.
GRID_DIR = Path("/Volumes/KINGSTON/part1/grid")               # wet.nc, basins_*.nc


# ---------------------------------------------------------------------------
# Box averaging
# ---------------------------------------------------------------------------
def layer_thickness(depth: np.ndarray) -> np.ndarray:
    """Layer thicknesses (m) from layer-centre depths: faces at 0 m and the
    midpoints between successive centres (bottom face extrapolated)."""
    edges = np.concatenate([
        [0.0],
        0.5 * (depth[:-1] + depth[1:]),
        [depth[-1] + 0.5 * (depth[-1] - depth[-2])],
    ])
    return np.diff(edges)


def box_means(s: np.ndarray, lat: np.ndarray, depth: np.ndarray) -> dict[str, np.ndarray]:
    """Volume-weighted box-mean salinity time series.

    ``s`` has shape (time, depth, lat); returns {box_name: (time,)}.  Weights are
    ``cos(lat) * layer_thickness``; NaN (land) cells are excluded from both the
    weighted sum and the weight normalisation.  Each box is averaged over its own
    depth range (``depth_max`` = None means the full water column), so the full
    depth axis must be supplied.
    """
    thick = layer_thickness(depth)                 # (depth,)
    coslat = np.cos(np.deg2rad(lat))               # (lat,)
    weight = thick[:, None] * coslat[None, :]      # (depth, lat)

    out: dict[str, np.ndarray] = {}
    for name, box in BOXES.items():
        dmax = box["depth_max"]
        dsel = np.ones(len(depth), dtype=bool) if dmax is None else (depth <= dmax)
        lsel = (lat >= box["lat_min"]) & (lat <= box["lat_max"])
        sub = s[:, dsel][:, :, lsel]               # (time, nd, nl)
        w = weight[dsel][:, lsel]                   # (nd, nl)
        valid = ~np.isnan(sub)
        wb = np.broadcast_to(w[None, :, :], sub.shape)
        num = np.sum(np.where(valid, sub * wb, 0.0), axis=(1, 2))
        den = np.sum(np.where(valid, wb, 0.0), axis=(1, 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            out[name] = np.where(den > 0, num / den, np.nan)
    return out


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------
def index_map(folder: Path) -> dict[int, Path]:
    """Map source itx index -> file path for one branch folder."""
    out = {}
    for f in glob.glob(str(folder / "*.nc")):
        m = re.search(r"itx(\d+)", os.path.basename(f))
        if m:
            out[int(m.group(1))] = Path(f)
    return out


def compute_track(path: Path) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Load one branch file and return (box means, lat, depth)."""
    with xr.open_dataset(path, decode_times=False) as ds:
        s = ds["s"].values.astype(np.float64)      # (time, depth, lat)
        lat = ds["lat"].values.astype(np.float64)
        depth = ds["depth"].values.astype(np.float64)
    return box_means(s, lat, depth), lat, depth


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------
def build_dataset(per_track: dict[str, dict[str, np.ndarray]],
                  co2: str, out_index: int, src_index: int,
                  src_files: dict[str, str]) -> xr.Dataset:
    """Assemble a (track, time) dataset from the per-track box means."""
    ntime = max(len(next(iter(m.values()))) for m in per_track.values())

    data_vars = {}
    for name, box in BOXES.items():
        arr = np.full((len(TRACKS), ntime), np.nan, dtype=np.float64)
        for it, tr in enumerate(TRACKS):
            series = per_track[tr][name]
            arr[it, :len(series)] = series
        data_vars[name] = (
            ("track", "time"), arr,
            {"long_name": box["long_name"], "units": "g/kg",
             "box_lat_min": box["lat_min"], "box_lat_max": box["lat_max"],
             "box_depth_max": _depth_attr(box)},
        )

    ds = xr.Dataset(
        data_vars,
        coords={
            "track": ("track", np.array(TRACKS, dtype=object),
                      {"standard_name": "Edge track branch "
                                        "(upper: going to ON, lower: going to OFF)"}),
            "time": ("time", np.arange(ntime, dtype="int32"),
                     {"standard_name": "Time", "units": "years (since edge tracking start)"}),
        },
        attrs={
            "info": "Box-mean salinity time series derived from Börner et al. (2025) "
                    "PlaSim-LSG zonally-averaged salinity (s_zonav) fields.",
            "source": "Data repository: https://doi.org/10.5281/zenodo.17053348",
            "processing": "compute_box_salinity.py — volume-weighted (cos(lat) x layer "
                          "thickness) mean over each box's (lat x depth) cells (see per-"
                          "variable box_depth_max), land cells excluded.",
            "co2_level": co2,
            "source_itx_index": src_index,
            "source_file_upper": src_files["upper"],
            "source_file_lower": src_files["lower"],
        },
    )
    return ds


def process_co2(src_root: Path, co2: str) -> int:
    lower = index_map(src_root / co2 / "lower")
    upper = index_map(src_root / co2 / "upper")
    common = sorted(set(lower) & set(upper))
    only_lower = sorted(set(lower) - set(upper))
    only_upper = sorted(set(upper) - set(lower))
    if only_lower or only_upper:
        print(f"  [warn] {co2}: unpaired indices skipped — "
              f"lower-only {only_lower}, upper-only {only_upper}")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    tag = CO2_TAG[co2]
    written = 0
    for out_index, src_index in enumerate(common):
        paths = {"lower": lower[src_index], "upper": upper[src_index]}
        per_track = {}
        lens = {}
        for tr in TRACKS:
            means, _, _ = compute_track(paths[tr])
            per_track[tr] = means
            lens[tr] = len(next(iter(means.values())))
        ds = build_dataset(
            per_track, co2, out_index, src_index,
            {tr: os.path.basename(paths[tr]) for tr in TRACKS})

        out_name = f"plasimelancholia_{co2}_edgetrack_{tag}{out_index:03d}.etc.nc"
        out_path = OUT_ROOT / out_name
        encoding = {n: {"_FillValue": np.nan} for n in BOXES}
        ds.to_netcdf(out_path, encoding=encoding)
        written += 1
        note = "" if lens["lower"] == lens["upper"] else \
            f"  (padded: upper={lens['upper']}, lower={lens['lower']})"
        print(f"  {out_name}  <- itx{src_index:03d}  "
              f"[{max(lens.values())} yr]{note}")
    return written


# ---------------------------------------------------------------------------
# Equilibrium (converged-state) files: 3D salinity -> Atlantic zonal mean -> boxes
# ---------------------------------------------------------------------------
_ATL_MASK = None


def atlantic_mask():
    """Boolean Atlantic mask (depth, lat, lon) from Börner's PlaRegion.

    Uses the same Atlantic definition (southern border -34°, Mediterranean /
    Hudson Bay / Arctic excluded) that produced the edge-track s_zonav fields,
    so the equilibrium zonal mean is consistent with the edge-track processing.
    """
    global _ATL_MASK
    if _ATL_MASK is None:
        import sys
        sys.path.insert(0, str(SCRIPT_DIR))      # vendored regions.py
        from regions import PlaRegion
        m = PlaRegion(path=str(GRID_DIR) + "/")
        _ATL_MASK = m.Atlantic3D() > 0           # (depth, lat, lon) bool
    return _ATL_MASK


def equilibrium_box_means(path: Path) -> dict[str, np.ndarray]:
    """Atlantic zonal mean of the full 3D salinity field, then the box means.

    The deep boxes span the full water column, so all depth levels are read.  To
    keep memory modest for the large (multi-GB) equilibrium files, the Atlantic
    zonal mean is accumulated in time chunks.
    """
    mask = np.asarray(atlantic_mask().values)             # (depth, lat, lon) bool
    with xr.open_dataset(path, decode_times=False) as ds:
        depth = ds["depth"].values.astype(np.float64)
        lat = ds["lat"].values.astype(np.float64)
        ntime = ds.dims["time"]
        svar = ds["s"]
        zonav = np.empty((ntime, len(depth), len(lat)), dtype=np.float64)
        chunk = 200
        for t0 in range(0, ntime, chunk):
            t1 = min(t0 + chunk, ntime)
            block = svar.isel(time=slice(t0, t1)).values      # (nt, depth, lat, lon)
            block = np.where(mask[None, :, :, :], block, np.nan)
            with warnings.catch_warnings():                    # all-land lats -> NaN
                warnings.simplefilter("ignore", category=RuntimeWarning)
                zonav[t0:t1] = np.nanmean(block, axis=3)       # Atlantic zonal mean
    return box_means(zonav, lat, depth)


def build_equilibrium_dataset(means: dict[str, np.ndarray], co2: str,
                              state: str, src_file: str) -> xr.Dataset:
    ntime = len(next(iter(means.values())))
    data_vars = {}
    for name, box in BOXES.items():
        data_vars[name] = (
            ("time",), means[name],
            {"long_name": box["long_name"], "units": "g/kg",
             "box_lat_min": box["lat_min"], "box_lat_max": box["lat_max"],
             "box_depth_max": _depth_attr(box)},
        )
    return xr.Dataset(
        data_vars,
        coords={"time": ("time", np.arange(ntime, dtype="int32"),
                         {"standard_name": "Time", "units": "years"})},
        attrs={
            "info": "Box-mean salinity time series derived from Börner et al. (2025) "
                    "PlaSim-LSG 3D salinity, reduced via the Atlantic zonal mean.",
            "source": "Data repository: https://doi.org/10.5281/zenodo.17053348",
            "processing": "compute_box_salinity.py — Atlantic zonal mean (PlaRegion "
                          "Atlantic3D, southern border -34) then volume-weighted "
                          "(cos(lat) x layer thickness) box mean over each box's depth "
                          "range (see per-variable box_depth_max), land cells excluded.",
            "co2_level": co2,
            "state": {"on": "AMOC-on", "of": "AMOC-off", "ed": "edge (saddle)"}[state],
            "source_file": src_file,
        },
    )


def process_equilibria(src_root: Path, co2: str) -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    written = 0
    for state in EQ_STATES:
        src = src_root / co2 / f"plasimedge_{co2}_{state}.tsv.nc"
        if not src.exists():
            print(f"  [warn] missing equilibrium file, skipping: {src.name}")
            continue
        means = equilibrium_box_means(src)
        ds = build_equilibrium_dataset(means, co2, state, src.name)
        out_path = OUT_ROOT / f"plasimelancholia_{co2}_{state}.etc.nc"
        ds.to_netcdf(out_path, encoding={n: {"_FillValue": np.nan} for n in BOXES})
        written += 1
        summ = "  ".join(
            f"{n.split('_')[1]}={np.nanmean(means[n]):.3f}" if np.isfinite(means[n]).any()
            else f"{n.split('_')[1]}=nan" for n in BOXES)
        print(f"  {out_path.name}  <- {src.name}  "
              f"[{len(means['salt_na'])} yr]  mean: {summ}")
    return written


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", type=Path, default=DEFAULT_SRC_ROOT,
                    help=f"source root (default: {DEFAULT_SRC_ROOT})")
    ap.add_argument("--co2", nargs="+", default=["285ppm", "360ppm"],
                    help="CO2 levels to process (default: 285ppm 360ppm)")
    ap.add_argument("--skip-edgetrack", action="store_true",
                    help="do not process the edge-track (itx/iter) files")
    ap.add_argument("--skip-equilibria", action="store_true",
                    help="do not process the on/of/ed equilibrium files")
    args = ap.parse_args()

    if not args.src.exists():
        raise SystemExit(f"Source root not found: {args.src}")

    print(f"Source : {args.src}")
    print(f"Output : {OUT_ROOT}")
    total = 0
    for co2 in args.co2:
        print(f"\n[{co2}]")
        if not args.skip_edgetrack:
            total += process_co2(args.src, co2)
        if not args.skip_equilibria:
            total += process_equilibria(args.src, co2)
    print(f"\nDone — wrote {total} files to {OUT_ROOT}")


if __name__ == "__main__":
    main()
