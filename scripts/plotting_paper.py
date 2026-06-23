"""
plotting_paper.py  –  AMOCPlaSim 4-panel paper figure.

Reads CSV files exported by plasim_export_paper_data.jl:
  data/plasim/paper/trajectories_{285ppm,360ppm}.csv
  data/plasim/paper/equilibria_{285ppm,360ppm}.csv
  data/plasim/paper/ellipses_{285ppm,360ppm}.csv
  data/plasim/paper/state_means_{285ppm,360ppm}.csv

Output: plots/plasim_paper.png

Run from the AMOCPlaSim directory or the project root:
    python scripts/plotting_paper.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    for ax in axes_bottom:
        ax.set_aspect("auto")
    axes_bottom[1].sharex(axes_bottom[0])
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

        else:
            # No AMOC strength variable available
            ax_top.text(0.5, 0.5,
                        "AMOC strength not available\nfor this model output",
                        ha="center", va="center", transform=ax_top.transAxes,
                        fontsize=8, color="gray")

        # Reference lines from equilibrium data
        for state, color in [("on", COL_ON), ("off", COL_OFF)]:
            sub = df_equil[df_equil["state"] == state]["amoc_strength"]
            if not sub.empty:
                ref = float(sub.iloc[-min(100, len(sub)):].mean())
                ax_top.axhline(ref, color=color, lw=1.0, ls="--", alpha=0.8)

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

        # 2. Mean position markers removed

        # 3. Edge-track trajectories with time shading
        on_tids  = df_trajs[df_trajs["label"] == "on"]["traj_id"].unique()
        off_tids = df_trajs[df_trajs["label"] == "off"]["traj_id"].unique()

        for tid in on_tids:
            sub = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
            ax_bot.plot(sub["x1"].values, sub["x2"].values,
                        color=COL_ON, alpha=0.4, lw=0.7)

        for tid in off_tids:
            sub = df_trajs[df_trajs["traj_id"] == tid].sort_values("time")
            ax_bot.plot(sub["x1"].values, sub["x2"].values,
                        color=COL_OFF, alpha=0.4, lw=0.7)

        ax_bot.set_xlabel("EOF 1")
        if col == 0:
            ax_bot.set_ylabel("EOF 2")
        else:
            ax_bot.tick_params(labelleft=False)
        add_panel_label(ax_bot, panel_labels[col + 2])

    x_lo = min(ax.get_xlim()[0] for ax in axes_top)
    x_hi = max(ax.get_xlim()[1] for ax in axes_top)
    for ax in axes_top:
        ax.set_xlim(x_lo, x_hi)

    out_path = PLOTS_DIR / "plasim_paper.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
