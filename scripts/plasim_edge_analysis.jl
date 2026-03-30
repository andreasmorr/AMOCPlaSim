"""
    plasim_edge_analysis.jl

Analyze PlaSim high-dimensional edge-state trajectories for current (360 ppm)
and heightened (720 ppm) CO2.

For each CO2 level the script:
  1. Loads all NetCDF edge-track files into a DataFrame.
  2. Classifies each trajectory as converging to AMOC-on or AMOC-off.
  3. Estimates attractor positions from the final states.
  4. Computes convergence times and edge-to-attractor distances.
  5. Produces three figures:
       Fig 1 — Scatter in (redu1, redu2) space, colored by label, edge marked
       Fig 2 — Bar chart of mean convergence times (on vs off, CO2 comparison)
       Fig 3 — Bar chart of mean edge-to-attractor distances

Run from the project root:
    julia --project scripts/plasim_edge_analysis.jl

Edit the DATA_DIR and N_FILES_* constants below to point to your NetCDF data.
"""

using DrWatson
@quickactivate "AMOCResilience"

include(srcdir("plasim_utils.jl"))

using DataFrames
using CairoMakie
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# Configuration — adjust paths and counts to your data layout
# ─────────────────────────────────────────────────────────────────────────────

# Directory containing the PlaSim NetCDF files
const DATA_DIR = datadir("plasim")

# CO2 level labels (used in filenames)
const CO2_LABEL_CURRENT  = "360ppm"
const CO2_LABEL_HEIGHTEN = "720ppm"

# Number of NetCDF files per CO2 level (files that don't exist are skipped)
const N_FILES_360 = 38
const N_FILES_720 = 38

# EOF dimensions to use from each file
const VARIABLE_NAMES = ["redu1", "redu2", "redu3"]
const N_DIMS = length(VARIABLE_NAMES)

# Convergence threshold in EOF units (adjust if needed)
const EPSILON_CONV = 0.05

# Fraction of trajectory end used to estimate attractor positions
const FINAL_FRACTION = 0.1

# ─────────────────────────────────────────────────────────────────────────────
# File pattern helpers
# ─────────────────────────────────────────────────────────────────────────────

file_pattern_360(i) = "plasimelancholia_$(CO2_LABEL_CURRENT)_edgetrack_iter$(i).etc.nc"
file_pattern_720(i) = "plasimelancholia_$(CO2_LABEL_HEIGHTEN)_edgetrack_iter$(i).etc.nc"

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

@info "Loading 360 ppm trajectories..."
df_360 = load_plasim_trajectories(;
    co2_label      = CO2_LABEL_CURRENT,
    data_dir       = DATA_DIR,
    n_files        = N_FILES_360,
    variable_names = VARIABLE_NAMES,
    file_pattern   = file_pattern_360,
)

@info "Loading 720 ppm trajectories..."
df_720 = load_plasim_trajectories(;
    co2_label      = CO2_LABEL_HEIGHTEN,
    data_dir       = DATA_DIR,
    n_files        = N_FILES_720,
    variable_names = VARIABLE_NAMES,
    file_pattern   = file_pattern_720,
)

# ─────────────────────────────────────────────────────────────────────────────
# Run resilience summary for both CO2 levels
# ─────────────────────────────────────────────────────────────────────────────

@info "Computing resilience summary for 360 ppm..."
summary_360 = plasim_resilience_summary(df_360, N_DIMS;
    final_fraction = FINAL_FRACTION,
    ε              = EPSILON_CONV,
)

@info "Computing resilience summary for 720 ppm..."
summary_720 = plasim_resilience_summary(df_720, N_DIMS;
    final_fraction = FINAL_FRACTION,
    ε              = EPSILON_CONV,
)

# ─────────────────────────────────────────────────────────────────────────────
# Print summary statistics
# ─────────────────────────────────────────────────────────────────────────────

println("\n=== PlaSim Edge State Analysis ===\n")

for (label, s) in [("360 ppm", summary_360), ("720 ppm", summary_720)]
    println("--- $label ---")
    println("  Trajectories AMOC-on  : $(s.n_on)")
    println("  Trajectories AMOC-off : $(s.n_off)")
    println("  Mean convergence time  (AMOC-on ) : $(round(s.mean_conv_time_on;  digits=1)) steps")
    println("  Mean convergence time  (AMOC-off) : $(round(s.mean_conv_time_off; digits=1)) steps")
    println("  Mean edge→attractor dist (on ) : $(round(s.mean_dist_on;  digits=4))")
    println("  Mean edge→attractor dist (off) : $(round(s.mean_dist_off; digits=4))")
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# Helper: collect x1, x2 per trajectory at initial time step
# ─────────────────────────────────────────────────────────────────────────────

"""
    trajectory_initial_states(df, n_dims) → (ids, Matrix)

Return trajectory IDs and a (n_traj × n_dims) matrix of initial states.
"""
function trajectory_initial_states(df::DataFrame, n_dims::Int)
    ids = sort(unique(df.trajectory_id))
    states = zeros(Float64, length(ids), n_dims)
    for (j, tid) in enumerate(ids)
        traj = sort(filter(r -> r.trajectory_id == tid, df), :time)
        for k in 1:n_dims
            states[j, k] = traj[1, Symbol("x$k")]
        end
    end
    return ids, states
end

ids_360, init_360 = trajectory_initial_states(df_360, N_DIMS)
ids_720, init_720 = trajectory_initial_states(df_720, N_DIMS)

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: Trajectory scatter in (redu1, redu2) space
# ─────────────────────────────────────────────────────────────────────────────

# Build per-trajectory (x1, x2) time series for plotting
function traj_x1x2(df::DataFrame, tid::Int)
    sub = sort(filter(r -> r.trajectory_id == tid, df), :time)
    return sub.x1, sub.x2
end

col_on  = :steelblue
col_off = :firebrick
col_edge = :black

fig1 = Figure(size = (1400, 600))

for (col_idx, (label, df_plot, summary, init_states, ids)) in enumerate([
    ("360 ppm (current CO₂)", df_360, summary_360, init_360, ids_360),
    ("720 ppm (2× CO₂)",      df_720, summary_720, init_720, ids_720),
])
    ax = Axis(fig1[1, col_idx];
        xlabel    = "redu1 (1st EOF)",
        ylabel    = "redu2 (2nd EOF)",
        title     = label,
        titlesize = 13,
    )

    # Plot each trajectory
    for (j, tid) in enumerate(ids)
        lbl = summary.attractor_labels[j]
        c   = lbl == 1 ? col_on : col_off
        x1, x2 = traj_x1x2(df_plot, tid)
        lines!(ax, x1, x2; color = (c, 0.35), linewidth = 0.8)
    end

    # Overlay attractor mean positions
    att_on  = summary.attractors[1]
    att_off = summary.attractors[2]
    scatter!(ax, [att_on[1]],  [att_on[2]];  color = col_on,  marker = :star5, markersize = 18,
             label = "AMOC-on attractor")
    scatter!(ax, [att_off[1]], [att_off[2]]; color = col_off, marker = :star5, markersize = 18,
             label = "AMOC-off attractor")

    # Overlay edge states (initial positions of all trajectories)
    scatter!(ax, init_states[:, 1], init_states[:, 2];
             color = (col_edge, 0.6), marker = :circle, markersize = 5,
             label = "Edge state")

    axislegend(ax; position = :rt, labelsize = 10)
end

# Add a manual legend for trajectory colors
Legend(fig1[1, 3],
    [LineElement(color = col_on), LineElement(color = col_off)],
    ["AMOC-on trajectory", "AMOC-off trajectory"],
    "Trajectory type",
    labelsize = 11,
)

fig1_path = plotsdir("plasim_trajectories_scatter.png")
wsave(fig1_path, fig1)
@info "Figure 1 saved to: $fig1_path"

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Bar chart — mean convergence times
# ─────────────────────────────────────────────────────────────────────────────

fig2 = Figure(size = (700, 500))
ax2  = Axis(fig2[1, 1];
    xlabel    = "",
    ylabel    = "Mean convergence time (steps)",
    title     = "Mean convergence times by CO₂ level",
    titlesize = 14,
)

categories  = ["AMOC-on\n360 ppm", "AMOC-off\n360 ppm",
               "AMOC-on\n720 ppm", "AMOC-off\n720 ppm"]
conv_values = [
    summary_360.mean_conv_time_on,
    summary_360.mean_conv_time_off,
    summary_720.mean_conv_time_on,
    summary_720.mean_conv_time_off,
]
bar_colors = [col_on, col_off, col_on, col_off]
bar_alpha   = [1.0, 1.0, 0.55, 0.55]  # lighter for 720 ppm

barplot!(ax2, 1:4, conv_values;
    color = [(c, a) for (c, a) in zip(bar_colors, bar_alpha)],
    width = 0.6,
)

# x-axis tick labels
ax2.xticks = (1:4, categories)

# Legend: shade for CO2 level
Legend(fig2[1, 2],
    [PolyElement(color = :gray80), PolyElement(color = :gray50)],
    ["360 ppm", "720 ppm"],
    "CO₂ level",
)

fig2_path = plotsdir("plasim_convergence_times_bar.png")
wsave(fig2_path, fig2)
@info "Figure 2 saved to: $fig2_path"

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: Bar chart — mean edge-to-attractor distances
# ─────────────────────────────────────────────────────────────────────────────

fig3 = Figure(size = (700, 500))
ax3  = Axis(fig3[1, 1];
    xlabel    = "",
    ylabel    = "Mean edge→attractor distance (EOF units)",
    title     = "Edge-state distance to attractor by CO₂ level",
    titlesize = 14,
)

dist_values = [
    summary_360.mean_dist_on,
    summary_360.mean_dist_off,
    summary_720.mean_dist_on,
    summary_720.mean_dist_off,
]

barplot!(ax3, 1:4, dist_values;
    color = [(c, a) for (c, a) in zip(bar_colors, bar_alpha)],
    width = 0.6,
)

ax3.xticks = (1:4, categories)

Legend(fig3[1, 2],
    [PolyElement(color = :gray80), PolyElement(color = :gray50)],
    ["360 ppm", "720 ppm"],
    "CO₂ level",
)

fig3_path = plotsdir("plasim_edge_distances_bar.png")
wsave(fig3_path, fig3)
@info "Figure 3 saved to: $fig3_path"

# ─────────────────────────────────────────────────────────────────────────────
# Save processed results
# ─────────────────────────────────────────────────────────────────────────────

results_to_save = Dict(
    "summary_360" => summary_360,
    "summary_720" => summary_720,
)

wsave(datadir("plasim", "resilience_summaries.jld2"), results_to_save)
@info "Results saved to: $(datadir("plasim", "resilience_summaries.jld2"))"

# Display figures if running interactively
display(fig1)
display(fig2)
display(fig3)
