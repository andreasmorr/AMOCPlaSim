#!/usr/bin/env python3
"""
extract_amoc_strength.py — AMOC-strength sidecar from the public PlaSim-LSG data.

Reads the AMOC strength (maximum of the Atlantic overturning streamfunction,
46-66N, Sv) from Part I of the Börner et al. PlaSim-LSG data set
(https://doi.org/10.5281/zenodo.17053348):

    <src>/<co2>/plasimedge_<co2>_{on,of,ed}.lsg.nc        equilibrium runs
    <src>/<co2>/plasimedge_<co2>_ed-iterations.zip        edge tracks, members
        plasimedge_<co2>_edgetrack_<tag>NNN.etc.nc        (tag = itx / iter)

and writes the compact sidecar

    data/custom_readouts/amoc_strength_timeseries.nc

Each variable is keyed like the matching salinity readout written by
``compute_box_salinity.py`` without ``.etc.nc``, e.g.

    plasimelancholia_285ppm_on
    plasimelancholia_285ppm_edgetrack_itx000

Equilibrium variables are one-dimensional over their own time dimension
``time__<key>``. Edge-track variables are ``(track, time__<key>)`` with
track 1 = upper (converging to AMOC-on) and track 2 = lower (AMOC-off).

The edge-track members of the zip archives are numbered contiguously in the
same way as the salinity readouts (at 285 ppm the s_zonav source indices
itx014 and itx025 do not exist, and both data sets skip them). If the salinity
readouts are present, the script checks that every series has the same number
of valid years as the matching salinity track and stops if not.

Usage (from the AMOCPlaSim directory):
    python scripts/extract_amoc_strength.py [--src SRC] [--out OUT]
"""
from __future__ import annotations

import argparse
import datetime as dt
import io
import re
import zipfile
from pathlib import Path

import netCDF4
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SRC = Path("/Volumes/KINGSTON/BoernerEtAl/part1")
READOUT_DIR = SCRIPT_DIR.parent / "data" / "custom_readouts"
DEFAULT_OUT = READOUT_DIR / "amoc_strength_timeseries.nc"

CO2_LEVELS = ["285ppm", "360ppm"]
CO2_TAG = {"285ppm": "itx", "360ppm": "iter"}
EQ_STATES = ["on", "of", "ed"]
DOI = "https://doi.org/10.5281/zenodo.17053348"


def _read(ds: netCDF4.Dataset) -> tuple[np.ndarray, np.ndarray | None, dict]:
    """Return (amoc_strength, time values or None, attributes) from a dataset."""
    var = ds["amoc_strength"]
    var.set_auto_mask(True)
    data = np.ma.filled(var[:].astype(np.float64), np.nan)
    attrs = {k: var.getncattr(k) for k in var.ncattrs() if k in ("long_name", "units")}
    time = None
    if "time" in ds.variables:
        time = np.ma.filled(ds["time"][:].astype(np.float64), np.nan)
    return data, time, attrs


def _write_var(out: netCDF4.Dataset, key: str, data: np.ndarray,
               time: np.ndarray | None, attrs: dict, source: str) -> None:
    dim_time = f"time__{key}"
    n_time = data.shape[-1]
    out.createDimension(dim_time, n_time)
    tv = out.createVariable(dim_time, "f8", (dim_time,))
    tv.long_name = f"Time coordinate copied from {source}"
    tv.units = "years"
    tv[:] = time if time is not None and len(time) == n_time else np.arange(n_time, dtype=float)

    dims = (dim_time,) if data.ndim == 1 else ("track", dim_time)
    v = out.createVariable(key, "f8", dims, fill_value=np.nan)
    v.long_name = attrs.get("long_name", "AMOC strength")
    v.units = attrs.get("units", "Sv")
    v.source_file = source
    v[:] = data


def _valid_track_lengths(path: Path) -> list[int] | None:
    """Valid (non-NaN) years per track of a salinity readout, if it exists."""
    if not path.exists():
        return None
    with netCDF4.Dataset(path) as ds:
        s = np.ma.filled(ds["salt_na"][:].astype(np.float64), np.nan)
    s = np.atleast_2d(s)
    return [int(np.isfinite(row).sum()) for row in s]


def _check(key: str, data: np.ndarray) -> None:
    expected = _valid_track_lengths(READOUT_DIR / f"{key}.etc.nc")
    if expected is None:
        return
    got = [int(np.isfinite(row).sum()) for row in np.atleast_2d(data)]
    if got != expected:
        raise SystemExit(f"{key}: AMOC series has {got} valid years per track, "
                         f"salinity readout has {expected}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", type=Path, default=DEFAULT_SRC,
                    help=f"Zenodo Part I root (default: {DEFAULT_SRC})")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help=f"output file (default: {DEFAULT_OUT})")
    args = ap.parse_args()
    if not args.src.exists():
        raise SystemExit(f"Source root not found: {args.src}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with netCDF4.Dataset(args.out, "w", format="NETCDF4") as out:
        out.title = "PlaSim-LSG AMOC-strength time series"
        out.created_by = "scripts/extract_amoc_strength.py"
        out.created_on = dt.datetime.now().isoformat(timespec="seconds")
        out.source = f"Börner et al. PlaSim-LSG data, Part I: {DOI}"
        out.source_directory = str(args.src)
        out.variable_key = ("Variables are keyed by the matching salinity readout "
                            "filename without .etc.nc.")
        out.createDimension("track", 2)
        tr = out.createVariable("track", "i4", ("track",))
        tr.long_name = "Edge-track branch index"
        tr.track_1 = "upper branch, converging to AMOC-on"
        tr.track_2 = "lower branch, converging to AMOC-off"
        tr[:] = [1, 2]

        for co2 in CO2_LEVELS:
            # Equilibrium runs: full-length .lsg.nc files.
            for state in EQ_STATES:
                src = args.src / co2 / f"plasimedge_{co2}_{state}.lsg.nc"
                key = f"plasimelancholia_{co2}_{state}"
                with netCDF4.Dataset(src) as ds:
                    data, time, attrs = _read(ds)
                _check(key, data)
                _write_var(out, key, data, time, attrs, src.name)
                n_written += 1
                print(f"  {key:45s} <- {src.name}  [{data.shape[-1]} yr]")

            # Edge tracks: .etc.nc members of the iterations archive.
            zpath = args.src / co2 / f"plasimedge_{co2}_ed-iterations.zip"
            tag = CO2_TAG[co2]
            pattern = re.compile(rf"plasimedge_{co2}_edgetrack_{tag}(\d{{3}})\.etc\.nc$")
            with zipfile.ZipFile(zpath) as zf:
                members = sorted((int(m.group(1)), name) for name in zf.namelist()
                                 if (m := pattern.search(name)))
                indices = [i for i, _ in members]
                if indices != list(range(len(indices))):
                    raise SystemExit(f"{zpath.name}: edge-track members not numbered "
                                     f"contiguously: {indices}")
                for idx, name in members:
                    key = f"plasimelancholia_{co2}_edgetrack_{tag}{idx:03d}"
                    raw = zf.read(name)
                    with netCDF4.Dataset("inmemory.nc", memory=raw) as ds:
                        data, time, attrs = _read(ds)
                    if data.ndim != 2 or data.shape[0] != 2:
                        raise SystemExit(f"{name}: expected (track=2, time), got {data.shape}")
                    _check(key, data)
                    _write_var(out, key, data, time, attrs, f"{zpath.name}:{Path(name).name}")
                    n_written += 1
                print(f"  {len(members)} edge tracks <- {zpath.name}")

    print(f"Done — wrote {n_written} AMOC-strength series to {args.out}")


if __name__ == "__main__":
    main()
