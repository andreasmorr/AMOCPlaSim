"""
plotting_paper.py  –  AMOCPlaSim 4-panel paper figure.

Reads CSV files exported by plasim_export_paper_data.jl:
  data/plasim/paper/trajectories_{285ppm,360ppm}.csv
  data/plasim/paper/equilibria_{285ppm,360ppm}.csv
  data/plasim/paper/ellipses_{285ppm,360ppm}.csv
  data/plasim/paper/state_means_{285ppm,360ppm}.csv

Output: plots/plasim_paper.pdf

Run from the AMOCPlaSim directory or the project root:
    python scripts/plotting_paper.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent          # AMOCPlaSim/scripts/
PLASIM_DIR = SCRIPT_DIR.parent                        # AMOCPlaSim/
UMBRELLA   = PLASIM_DIR.parent                        # AMOCResilience/
DATA_DIR   = PLASIM_DIR / "data" / "plasim" / "paper"
PLOTS_DIR  = PLASIM_DIR / "plots"

sys.path.insert(0, str(UMBRELLA))
from amoc_plot_style import (
    COL_ON, COL_OFF, COL_EDGE,
    TRAJ_COLORS,
    make_paper_figure, add_panel_label, savefig_pdf,
)

# ---------------------------------------------------------------------------
# Trajectory alpha shading helper
# ---------------------------------------------------------------------------

def plot_traj_shaded(ax, x, y, color, alpha_start=0.10, alpha_end=0.85, lw=0.8):
    """Plot trajectory as a LineCollection with alpha increasing with time."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # Drop NaNs
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 2:
        return
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segs   = np.concatenate([points[:-1], points[1:]], axis=1)
    n      = len(segs)
    alphas = np.linspace(alpha_start, alpha_end, n)
    colors_with_alpha = [(*matplotlib.colors.to_rgb(color), a) for a in alphas]
    lc = LineCollection(segs, colors=colors_with_alpha, linewidth=lw)
    ax.add_collection(lc)
    ax.autoscale()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

SCENARIOS = [
    ("285ppm", "285 ppm (pre-industrial)"),
    ("360ppm", "360 ppm (current CO\u2082)"),
]


def load_scenario(co2_label: str) -> dict | None:
    """Load all CSV data for one CO2 scenario. Returns None if files are missing."""
    files = {
        "trajs":  DATA_DIR / f"trajectories_{co2_label}.csv",
        "equil":  DATA_DIR / f"equilibria_{co2_label}.csv",
        "ellip":  DATA_DIR / f"ellipses_{co2_label}.csv",
        "means":  DATA_DIR / f"state_means_{co2_label}.csv",
    }
    missing = [str(p) for p in files.values() if not p.exists()]
    if missing:
        print(f"[{co2_label}] Missing CSV files:")
        for p in missing:
            print(f"  {p}")
        return None
    return {k: pd.read_csv(v) for k, v in files.items()}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    data = {}
    for co2_label, _ in SCENARIOS:
        d = load_scenario(co2_label)
        if d is None:
            print(f"\nRun plasim_export_paper_data.jl first to generate CSVs in {DATA_DIR}.")
            sys.exit(1)
        data[co2_label] = d

    fig, axes_top, axes_bottom = make_paper_figure()
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for col, (co2_label, title) in enumerate(SCENARIOS):
        d         = data[co2_label]
        df_trajs  = d["trajs"]
        df_equil  = d["equil"]
        df_ellip  = d["ellip"]
        df_means  = d["means"]

        ax_top = axes_top[col]
        ax_bot = axes_bottom[col]

        # ── TOP panel: AMOC strength vs time ─────────────────────────────────
        has_amoc = "amoc_strength" in df_trajs.columns

        if has_amoc:
            # Separate on/off trajectories
            on_tids  = df_trajs[df_trajs["label"] == "on"]["traj_id"].unique()
            off_tids = df_trajs[df_trajs["label"] == "off"]["traj_id"].unique()

            for tid in on_tids:
                sub = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
                amoc_arr = sub["amoc_strength"].values
                if np.all(np.isnan(amoc_arr)):
                    continue
                ax_top.plot(sub["time"].values, amoc_arr,
                            color=COL_ON, lw=0.7, alpha=0.5)

            for tid in off_tids:
                sub = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
                amoc_arr = sub["amoc_strength"].values
                if np.all(np.isnan(amoc_arr)):
                    continue
                ax_top.plot(sub["time"].values, amoc_arr,
                            color=COL_OFF, lw=0.7, alpha=0.5)

            # Proxy legend entries
            ax_top.plot([], [], color=COL_ON,  lw=1.2, label="AMOC-on")
            ax_top.plot([], [], color=COL_OFF, lw=1.2, label="AMOC-off")
            if col == 0:
                ax_top.legend(loc="upper right", framealpha=0.8, fontsize=7)
        else:
            # No AMOC strength variable available
            ax_top.text(0.5, 0.5,
                        "AMOC strength not available\nfor this model output",
                        ha="center", va="center", transform=ax_top.transAxes,
                        fontsize=8, color="gray")

        ax_top.set_title(title, fontsize=9)
        if col == 0:
            ax_top.set_ylabel("AMOC strength (Sv)")
        else:
            ax_top.tick_params(labelleft=False)
        ax_top.set_xlabel("Time (model years)")
        add_panel_label(ax_top, panel_labels[col])

        # ── BOTTOM panel: EOF1 vs EOF2 phase portrait ─────────────────────────
        # 1. Gaussian ellipses (thick dashed)
        for state_name, color in [
            ("on",   COL_ON  ),
            ("off",  COL_OFF ),
            ("edge", COL_EDGE),
        ]:
            sub = df_ellip[df_ellip["state"] == state_name]
            if not sub.empty:
                ax_bot.plot(sub["x"].values, sub["y"].values,
                            color=color, lw=1.8, ls="--", alpha=0.85, zorder=2)

        # 2. Mean position stars
        for _, row in df_means.iterrows():
            color = (COL_ON   if row["state"] == "on"   else
                     COL_OFF  if row["state"] == "off"  else
                     COL_EDGE)
            ax_bot.scatter([row["x1"]], [row["x2"]],
                           marker="*", s=100, color=color, zorder=5,
                           edgecolors="white", linewidths=0.4)

        # 3. Edge-track trajectories with time shading
        on_tids  = df_trajs[df_trajs["label"] == "on"]["traj_id"].unique()
        off_tids = df_trajs[df_trajs["label"] == "off"]["traj_id"].unique()

        for tid in on_tids:
            sub = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
            plot_traj_shaded(ax_bot, sub["x1"].values, sub["x2"].values,
                             color=COL_ON, alpha_start=0.08, alpha_end=0.7, lw=0.7)

        for tid in off_tids:
            sub = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
            plot_traj_shaded(ax_bot, sub["x1"].values, sub["x2"].values,
                             color=COL_OFF, alpha_start=0.08, alpha_end=0.7, lw=0.7)

        ax_bot.set_xlabel("EOF 1")
        if col == 0:
            ax_bot.set_ylabel("EOF 2")
        else:
            ax_bot.tick_params(labelleft=False)
        add_panel_label(ax_bot, panel_labels[col + 2])

    out_path = PLOTS_DIR / "plasim_paper.pdf"
    savefig_pdf(fig, out_path)
    plt.close(fig)


if __name__ == "__main__":
    main()
