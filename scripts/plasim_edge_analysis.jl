"""
    plasim_edge_analysis.jl

Analyze PlaSim high-dimensional edge-state trajectories for current (360 ppm)
and preindustrial (285 ppm) CO2.

For each CO2 level the script:
  1. Loads all NetCDF edge-track files into a DataFrame.
  2. Classifies each trajectory as converging to AMOC-on or AMOC-off.
  3. Loads attractor positions from converged `_on.etc.nc` / `_off.etc.nc` files.
  4. Loads the edge-state position from the converged `_ed.etc.nc` file.
  5. Computes convergence times and edge-to-attractor distances.
  6. Produces three figures:
       Fig 1 — Scatter in (EOF1, EOF2) space, colored by label, edge marked
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
using LinearAlgebra
using CSV

# ─────────────────────────────────────────────────────────────────────────────
# Configuration — adjust paths and counts to your data layout
# ─────────────────────────────────────────────────────────────────────────────

# Directory containing the PlaSim NetCDF files
const DATA_DIR = datadir("plasim")

# CO2 level labels (used in filenames)
const CO2_LABEL_CURRENT  = "360ppm"
const CO2_LABEL_PREINDUSTRIAL = "285ppm"

# Number of NetCDF files per CO2 level (files that don't exist are skipped)
const N_FILES_360 = 38
const N_FILES_285 = 37

# EOF dimensions to use from each file
const VARIABLE_NAMES = ["redu1", "redu2", "redu3"]
const N_DIMS = length(VARIABLE_NAMES)

# Convergence thresholds in EOF units — separate for each attractor (adjust if needed)
const EPSILON_ON  = 0.1   # threshold for AMOC-on trajectories
const EPSILON_OFF = 0.2   # threshold for AMOC-off trajectories

# Fraction of trajectory end used to classify trajectories (not for attractor
# estimation — attractor positions come from the dedicated equilibrium files)
const FINAL_FRACTION = 0.1

# Number of standard deviations for Gaussian ellipse plots and ellipse-based metrics
const ELLIPSE_SIGMA = 3

# ─────────────────────────────────────────────────────────────────────────────
# File pattern helpers
# ─────────────────────────────────────────────────────────────────────────────

file_pattern_360(i) = "plasimelancholia_$(CO2_LABEL_CURRENT)_edgetrack_iter$(lpad(i-1, 3, '0')).etc.nc"
file_pattern_285(i) = "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_edgetrack_itx$(lpad(i-1, 3, '0')).etc.nc"

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

@info "Loading 285 ppm trajectories..."
df_285 = load_plasim_trajectories(;
    co2_label      = CO2_LABEL_PREINDUSTRIAL,
    data_dir       = DATA_DIR,
    n_files        = N_FILES_285,
    variable_names = VARIABLE_NAMES,
    file_pattern   = file_pattern_285,
)

@info "Loading 360 ppm trajectories..."
df_360 = load_plasim_trajectories(;
    co2_label      = CO2_LABEL_CURRENT,
    data_dir       = DATA_DIR,
    n_files        = N_FILES_360,
    variable_names = VARIABLE_NAMES,
    file_pattern   = file_pattern_360,
)

# ─────────────────────────────────────────────────────────────────────────────
# Load converged equilibrium states (attractor & edge positions)
#
# These files contain trajectories that have already converged to the vicinity
# of the AMOC-on, AMOC-off, and edge (saddle) equilibria.  Their time-mean
# gives better attractor/edge-state estimates than the final states of the
# bisection trajectories.
# ─────────────────────────────────────────────────────────────────────────────

@info "Loading converged equilibrium states for 285 ppm..."
att_on_285  = load_plasim_state_mean(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_on.etc.nc"),
    VARIABLE_NAMES)
att_off_285 = load_plasim_state_mean(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_of.etc.nc"),
    VARIABLE_NAMES)
edge_285    = load_plasim_state_mean(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_ed.etc.nc"),
    VARIABLE_NAMES)
attractors_285 = Dict{Int, Vector{Float64}}(1 => att_on_285, 2 => att_off_285)

# Full time series for the equilibrium runs (used for trajectory plots)
ts_on_285  = load_plasim_state_timeseries(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_on.etc.nc"),
    VARIABLE_NAMES)
ts_off_285 = load_plasim_state_timeseries(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_of.etc.nc"),
    VARIABLE_NAMES)
ts_ed_285  = load_plasim_state_timeseries(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_ed.etc.nc"),
    VARIABLE_NAMES)

@info "Loading converged equilibrium states for 360 ppm..."
att_on_360  = load_plasim_state_mean(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_on.etc.nc"),
    VARIABLE_NAMES)
att_off_360 = load_plasim_state_mean(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_of.etc.nc"),
    VARIABLE_NAMES)
edge_360    = load_plasim_state_mean(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_ed.etc.nc"),
    VARIABLE_NAMES)
attractors_360 = Dict{Int, Vector{Float64}}(1 => att_on_360, 2 => att_off_360)

# Full time series for the equilibrium runs (used for trajectory plots)
ts_on_360  = load_plasim_state_timeseries(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_on.etc.nc"),
    VARIABLE_NAMES)
ts_off_360 = load_plasim_state_timeseries(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_of.etc.nc"),
    VARIABLE_NAMES)
ts_ed_360  = load_plasim_state_timeseries(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_ed.etc.nc"),
    VARIABLE_NAMES)

# ─────────────────────────────────────────────────────────────────────────────
# Local variability (covariance matrices needed for ellipse-based metrics)
# ─────────────────────────────────────────────────────────────────────────────

@info "Computing local variability for 285 ppm..."
var_on_285  = compute_local_variability(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_on.etc.nc"), VARIABLE_NAMES)
var_off_285 = compute_local_variability(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_of.etc.nc"), VARIABLE_NAMES)
var_ed_285  = compute_local_variability(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_ed.etc.nc"), VARIABLE_NAMES)

@info "Computing local variability for 360 ppm..."
var_on_360  = compute_local_variability(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_on.etc.nc"), VARIABLE_NAMES)
var_off_360 = compute_local_variability(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_of.etc.nc"), VARIABLE_NAMES)
var_ed_360  = compute_local_variability(
    joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_ed.etc.nc"), VARIABLE_NAMES)

# ─────────────────────────────────────────────────────────────────────────────
# Run resilience summary for both CO2 levels
# ─────────────────────────────────────────────────────────────────────────────

@info "Computing resilience summary for 285 ppm..."
summary_285 = plasim_resilience_summary(df_285, N_DIMS;
    final_fraction = FINAL_FRACTION,
    ε_on           = EPSILON_ON,
    ε_off          = EPSILON_OFF,
    attractors     = attractors_285,
    edge_state     = edge_285,
    cov_on         = var_on_285.covariance,
    cov_off        = var_off_285.covariance,
    edge_cov       = var_ed_285.covariance,
    sigma          = ELLIPSE_SIGMA,
    check_dims     = 2,
)

@info "Computing resilience summary for 360 ppm..."
summary_360 = plasim_resilience_summary(df_360, N_DIMS;
    final_fraction = FINAL_FRACTION,
    ε_on           = EPSILON_ON,
    ε_off          = EPSILON_OFF,
    attractors     = attractors_360,
    edge_state     = edge_360,
    cov_on         = var_on_360.covariance,
    cov_off        = var_off_360.covariance,
    edge_cov       = var_ed_360.covariance,
    sigma          = ELLIPSE_SIGMA,
    check_dims     = 2,
)

# ─────────────────────────────────────────────────────────────────────────────
# Print summary statistics
# ─────────────────────────────────────────────────────────────────────────────

println("\n=== PlaSim Edge State Analysis ===\n")

for (label, s) in [("285 ppm", summary_285), ("360 ppm", summary_360)]
    println("--- $label ---")
    println("  Trajectories AMOC-on  : $(s.n_on)  (converged: $(s.n_converged_on))")
    println("  Trajectories AMOC-off : $(s.n_off)  (converged: $(s.n_converged_off))")
    println("  Mean convergence time  (AMOC-on ) : $(round(s.mean_conv_time_on;  digits=1)) yr  [converged only]")
    println("  Mean convergence time  (AMOC-off) : $(round(s.mean_conv_time_off; digits=1)) yr  [converged only]")
    println("  Mean edge→attractor dist (on ) : $(round(s.mean_dist_on;  digits=4))")
    println("  Mean edge→attractor dist (off) : $(round(s.mean_dist_off; digits=4))")
    println()
end

ids_285 = summary_285.trajectory_ids
ids_360 = summary_360.trajectory_ids

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: Trajectory scatter in (EOF1, EOF2) space
# ─────────────────────────────────────────────────────────────────────────────

# Build per-trajectory time series for plotting
function traj_x1x2(df::DataFrame, tid::Int)
    sub = sort(filter(r -> r.trajectory_id == tid, df), :time)
    return sub.x1, sub.x2
end

function traj_x1x2x3(df::DataFrame, tid::Int)
    sub = sort(filter(r -> r.trajectory_id == tid, df), :time)
    return sub.x1, sub.x2, sub.x3
end

# Return a Bool mask selecting rows of a time-series matrix that have no NaN
valid_rows(ts::Matrix) = [!any(isnan, ts[t, :]) for t in axes(ts, 1)]

# Return (x, y) points tracing a circle of given radius around (cx, cy)
function circle_points(cx, cy, r; n = 120)
    θ = LinRange(0, 2π, n + 1)
    return cx .+ r .* cos.(θ), cy .+ r .* sin.(θ)
end

# Return (X, Y, Z) surface matrices for a sphere of radius r centered at (cx, cy, cz)
function sphere_surface(cx, cy, cz, r; n = 40)
    θ = LinRange(0,  π, n)   # polar
    φ = LinRange(0, 2π, n)   # azimuthal
    X = [cx + r * sin(t) * cos(p) for t in θ, p in φ]
    Y = [cy + r * sin(t) * sin(p) for t in θ, p in φ]
    Z = [cz + r * cos(t)          for t in θ, p in φ]
    return X, Y, Z
end

# Return (x, y) points tracing the n-sigma contour of a 2D Gaussian fitted to ts[:,dim1] and ts[:,dim2]
function gaussian_ellipse_points(ts::Matrix, dim1::Int = 1, dim2::Int = 2; n = 120, sigma = 1)
    valid = [!any(isnan, ts[t, :]) for t in axes(ts, 1)]
    data  = ts[valid, [dim1, dim2]]
    μ     = vec(mean(data, dims = 1))
    C     = cov(data)
    vals, vecs = eigen(Symmetric(C))
    vals  = max.(vals, 0.0)   # guard against tiny negative floats
    θ     = LinRange(0, 2π, n + 1)
    # unit circle scaled by sigma * sqrt(eigenvalues), then rotated to data frame
    pts   = vecs * (sigma .* sqrt.(vals) .* [cos.(θ)'; sin.(θ)'])
    return μ[1] .+ pts[1, :], μ[2] .+ pts[2, :]
end

# Return (X, Y, Z) surface matrices for the n-sigma ellipsoid of a 3D Gaussian fitted to ts[:,1:3]
function gaussian_ellipsoid_surface(ts::Matrix; n = 40, sigma = 1)
    valid = [!any(isnan, ts[t, :]) for t in axes(ts, 1)]
    data  = ts[valid, 1:3]
    μ     = vec(mean(data, dims = 1))
    C     = cov(data)
    vals, vecs = eigen(Symmetric(C))
    vals  = max.(vals, 0.0)
    sv    = sigma .* sqrt.(vals)
    θ = LinRange(0,  π, n)   # polar
    φ = LinRange(0, 2π, n)   # azimuthal
    X = zeros(n, n); Y = zeros(n, n); Z = zeros(n, n)
    for i in 1:n, j in 1:n
        v = vecs * (sv .* [sin(θ[i]) * cos(φ[j]), sin(θ[i]) * sin(φ[j]), cos(θ[i])])
        X[i, j] = μ[1] + v[1]
        Y[i, j] = μ[2] + v[2]
        Z[i, j] = μ[3] + v[3]
    end
    return X, Y, Z
end


col_on  = :steelblue
col_off = :firebrick
col_edge = :black

fig1 = Figure(size = (1400, 600))

for (col_idx, (label, df_plot, summary, edge_st, ids, ts_on, ts_off, ts_ed)) in enumerate([
    ("285 ppm (pre-industrial)", df_285, summary_285, edge_285, ids_285, ts_on_285, ts_off_285, ts_ed_285),
    ("360 ppm (current CO₂)",   df_360, summary_360, edge_360, ids_360, ts_on_360, ts_off_360, ts_ed_360),
])
    ax = Axis(fig1[1, col_idx];
        xlabel    = "1st EOF",
        ylabel    = "2nd EOF",
        title     = label,
        titlesize = 13,
    )

    # Plot each trajectory — solid for converged, dashed+faint for unconverged
    for (j, tid) in enumerate(ids)
        lbl       = summary.attractor_labels[j]
        c         = lbl == 1 ? col_on : col_off
        did_conv  = summary.converged[tid]
        x1, x2   = traj_x1x2(df_plot, tid)
        if did_conv
            lines!(ax, x1, x2; color = (c, 0.6), linewidth = 1.2, linestyle = :solid)
        else
            lines!(ax, x1, x2; color = (c, 0.2), linewidth = 0.7, linestyle = :dash)
        end
    end

    # Overlay attractor mean positions and 1σ Gaussian ellipses
    att_on  = summary.attractors[1]
    att_off = summary.attractors[2]
    scatter!(ax, [att_on[1]],  [att_on[2]];  color = col_on,  marker = :star5, markersize = 18,
             label = "AMOC-on attractor")
    scatter!(ax, [att_off[1]], [att_off[2]]; color = col_off, marker = :star5, markersize = 18,
             label = "AMOC-off attractor")

    ex_on,  ey_on  = gaussian_ellipse_points(ts_on,  1, 2; sigma = ELLIPSE_SIGMA)
    ex_off, ey_off = gaussian_ellipse_points(ts_off, 1, 2; sigma = ELLIPSE_SIGMA)
    ex_ed,  ey_ed  = gaussian_ellipse_points(ts_ed,  1, 2; sigma = ELLIPSE_SIGMA)
    lines!(ax, ex_on,  ey_on;  color = (col_on,  0.8), linewidth = 1.4, linestyle = :dot,
           label = "$(ELLIPSE_SIGMA)σ AMOC-on")
    lines!(ax, ex_off, ey_off; color = (col_off, 0.8), linewidth = 1.4, linestyle = :dot,
           label = "$(ELLIPSE_SIGMA)σ AMOC-off")
    lines!(ax, ex_ed,  ey_ed;  color = (col_edge, 0.8), linewidth = 1.4, linestyle = :dot,
           label = "$(ELLIPSE_SIGMA)σ edge")

    # Overlay edge state (single converged position from _ed.etc.nc)
    scatter!(ax, [edge_st[1]], [edge_st[2]];
             color = col_edge, marker = :diamond, markersize = 14,
             label = "Edge state")

    axislegend(ax; position = :rt, labelsize = 10)
end

# Add a manual legend for trajectory colors and convergence status
Legend(fig1[1, 3],
    [
        LineElement(color = col_on,  linestyle = :solid, linewidth = 1.2),
        LineElement(color = col_off, linestyle = :solid, linewidth = 1.2),
        LineElement(color = :gray40, linestyle = :dash,  linewidth = 0.7),
    ],
    ["AMOC-on (converged)", "AMOC-off (converged)", "Not converged"],
    "Trajectory type",
    labelsize = 11,
)

fig1_path = plotsdir("plasim_trajectories_scatter.png")
wsave(fig1_path, fig1)
@info "Figure 1 saved to: $fig1_path"

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1b: Trajectory scatter in full (EOF1, EOF2, EOF3) space
# ─────────────────────────────────────────────────────────────────────────────

fig1b = Figure(size = (1400, 700))

for (col_idx, (label, df_plot, summary, edge_st, ids, ts_on, ts_off, ts_ed)) in enumerate([
    ("285 ppm (pre-industrial)", df_285, summary_285, edge_285, ids_285, ts_on_285, ts_off_285, ts_ed_285),
    ("360 ppm (current CO₂)",   df_360, summary_360, edge_360, ids_360, ts_on_360, ts_off_360, ts_ed_360),
])
    ax3d = Axis3(fig1b[1, col_idx];
        xlabel    = "EOF 1",
        ylabel    = "EOF 2",
        zlabel    = "EOF 3",
        title     = label,
        titlesize = 13,
    )

    # Plot each trajectory — solid for converged, dashed+faint for unconverged
    for (j, tid) in enumerate(ids)
        lbl      = summary.attractor_labels[j]
        c        = lbl == 1 ? col_on : col_off
        did_conv = summary.converged[tid]
        x1, x2, x3 = traj_x1x2x3(df_plot, tid)
        if did_conv
            lines!(ax3d, x1, x2, x3; color = (c, 0.6), linewidth = 1.2, linestyle = :solid)
        else
            lines!(ax3d, x1, x2, x3; color = (c, 0.2), linewidth = 0.7, linestyle = :dash)
        end
    end

    att_on  = summary.attractors[1]
    att_off = summary.attractors[2]
    scatter!(ax3d, [att_on[1]],  [att_on[2]],  [att_on[3]];
             color = col_on,  marker = :star5, markersize = 18, label = "AMOC-on attractor")
    scatter!(ax3d, [att_off[1]], [att_off[2]], [att_off[3]];
             color = col_off, marker = :star5, markersize = 18, label = "AMOC-off attractor")
    scatter!(ax3d, [edge_st[1]], [edge_st[2]], [edge_st[3]];
             color = col_edge, marker = :diamond, markersize = 14, label = "Edge state")

    # 1σ Gaussian ellipsoids for each state
    Xon, Yon, Zon = gaussian_ellipsoid_surface(ts_on;  sigma = ELLIPSE_SIGMA)
    surface!(ax3d, Xon, Yon, Zon; color = (col_on,  0.12), shading = NoShading)
    Xof, Yof, Zof = gaussian_ellipsoid_surface(ts_off; sigma = ELLIPSE_SIGMA)
    surface!(ax3d, Xof, Yof, Zof; color = (col_off, 0.12), shading = NoShading)
    Xed, Yed, Zed = gaussian_ellipsoid_surface(ts_ed;  sigma = ELLIPSE_SIGMA)
    surface!(ax3d, Xed, Yed, Zed; color = (col_edge, 0.12), shading = NoShading)

    axislegend(ax3d; position = :rt, labelsize = 10)
end

Legend(fig1b[1, 3],
    [
        LineElement(color = col_on,  linestyle = :solid, linewidth = 1.2),
        LineElement(color = col_off, linestyle = :solid, linewidth = 1.2),
        LineElement(color = :gray40, linestyle = :dash,  linewidth = 0.7),
    ],
    ["AMOC-on (converged)", "AMOC-off (converged)", "Not converged"],
    "Trajectory type",
    labelsize = 11,
)

fig1b_path = plotsdir("plasim_trajectories_scatter_3d.png")
wsave(fig1b_path, fig1b)
@info "Figure 1b saved to: $fig1b_path"

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1c: Equilibrium runs in (EOF1, EOF2) space
# ─────────────────────────────────────────────────────────────────────────────

fig1c = Figure(size = (1400, 600))

for (col_idx, (label, ts_on, ts_off, ts_ed, summary, edge_st)) in enumerate([
    ("285 ppm (pre-industrial)", ts_on_285, ts_off_285, ts_ed_285, summary_285, edge_285),
    ("360 ppm (current CO₂)",   ts_on_360, ts_off_360, ts_ed_360, summary_360, edge_360),
])
    ax = Axis(fig1c[1, col_idx];
        xlabel    = "1st EOF",
        ylabel    = "2nd EOF",
        title     = label,
        titlesize = 13,
    )

    vm_on  = valid_rows(ts_on)
    vm_off = valid_rows(ts_off)
    vm_ed  = valid_rows(ts_ed)
    lines!(ax, ts_on[vm_on, 1],  ts_on[vm_on, 2];
           color = (col_on,  0.7), linewidth = 1.2, linestyle = :solid, label = "AMOC-on eq. run")
    lines!(ax, ts_off[vm_off, 1], ts_off[vm_off, 2];
           color = (col_off, 0.7), linewidth = 1.2, linestyle = :solid, label = "AMOC-off eq. run")
    lines!(ax, ts_ed[vm_ed, 1],  ts_ed[vm_ed, 2];
           color = (col_edge, 0.6), linewidth = 1.0, linestyle = :solid, label = "Edge eq. run")

    att_on  = summary.attractors[1]
    att_off = summary.attractors[2]
    scatter!(ax, [att_on[1]],  [att_on[2]];  color = col_on,  marker = :star5, markersize = 18,
             label = "AMOC-on attractor")
    scatter!(ax, [att_off[1]], [att_off[2]]; color = col_off, marker = :star5, markersize = 18,
             label = "AMOC-off attractor")
    scatter!(ax, [edge_st[1]], [edge_st[2]];
             color = col_edge, marker = :diamond, markersize = 14, label = "Edge state")

    ex_on,  ey_on  = gaussian_ellipse_points(ts_on,  1, 2; sigma = ELLIPSE_SIGMA)
    ex_off, ey_off = gaussian_ellipse_points(ts_off, 1, 2; sigma = ELLIPSE_SIGMA)
    ex_ed,  ey_ed  = gaussian_ellipse_points(ts_ed,  1, 2; sigma = ELLIPSE_SIGMA)
    lines!(ax, ex_on,  ey_on;  color = (col_on,  0.8), linewidth = 1.4, linestyle = :dot,
           label = "$(ELLIPSE_SIGMA)σ AMOC-on")
    lines!(ax, ex_off, ey_off; color = (col_off, 0.8), linewidth = 1.4, linestyle = :dot,
           label = "$(ELLIPSE_SIGMA)σ AMOC-off")
    lines!(ax, ex_ed,  ey_ed;  color = (col_edge, 0.8), linewidth = 1.4, linestyle = :dot,
           label = "$(ELLIPSE_SIGMA)σ edge")

    axislegend(ax; position = :rt, labelsize = 10)
end

Legend(fig1c[1, 3],
    [
        LineElement(color = col_on,   linestyle = :solid, linewidth = 1.2),
        LineElement(color = col_off,  linestyle = :solid, linewidth = 1.2),
        LineElement(color = col_edge, linestyle = :solid, linewidth = 1.0),
    ],
    ["AMOC-on eq. run", "AMOC-off eq. run", "Edge eq. run"],
    "Equilibrium run",
    labelsize = 11,
)

fig1c_path = plotsdir("plasim_equilibrium_scatter.png")
wsave(fig1c_path, fig1c)
@info "Figure 1c saved to: $fig1c_path"

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1d: Equilibrium runs in (EOF1, EOF2, EOF3) space
# ─────────────────────────────────────────────────────────────────────────────

fig1d = Figure(size = (1400, 700))

for (col_idx, (label, ts_on, ts_off, ts_ed, summary, edge_st)) in enumerate([
    ("285 ppm (pre-industrial)", ts_on_285, ts_off_285, ts_ed_285, summary_285, edge_285),
    ("360 ppm (current CO₂)",   ts_on_360, ts_off_360, ts_ed_360, summary_360, edge_360),
])
    ax3d = Axis3(fig1d[1, col_idx];
        xlabel    = "EOF 1",
        ylabel    = "EOF 2",
        zlabel    = "EOF 3",
        title     = label,
        titlesize = 13,
    )

    vm_on  = valid_rows(ts_on)
    vm_off = valid_rows(ts_off)
    vm_ed  = valid_rows(ts_ed)
    lines!(ax3d, ts_on[vm_on, 1],  ts_on[vm_on, 2],  ts_on[vm_on, 3];
           color = (col_on,  0.7), linewidth = 1.2, linestyle = :solid, label = "AMOC-on eq. run")
    lines!(ax3d, ts_off[vm_off, 1], ts_off[vm_off, 2], ts_off[vm_off, 3];
           color = (col_off, 0.7), linewidth = 1.2, linestyle = :solid, label = "AMOC-off eq. run")
    lines!(ax3d, ts_ed[vm_ed, 1],  ts_ed[vm_ed, 2],   ts_ed[vm_ed, 3];
           color = (col_edge, 0.6), linewidth = 1.0, linestyle = :solid, label = "Edge eq. run")

    att_on  = summary.attractors[1]
    att_off = summary.attractors[2]
    scatter!(ax3d, [att_on[1]],  [att_on[2]],  [att_on[3]];
             color = col_on,  marker = :star5, markersize = 18, label = "AMOC-on attractor")
    scatter!(ax3d, [att_off[1]], [att_off[2]], [att_off[3]];
             color = col_off, marker = :star5, markersize = 18, label = "AMOC-off attractor")
    scatter!(ax3d, [edge_st[1]], [edge_st[2]], [edge_st[3]];
             color = col_edge, marker = :diamond, markersize = 14, label = "Edge state")

    Xon, Yon, Zon = gaussian_ellipsoid_surface(ts_on;  sigma = ELLIPSE_SIGMA)
    surface!(ax3d, Xon, Yon, Zon; color = (col_on,  0.12), shading = NoShading)
    Xof, Yof, Zof = gaussian_ellipsoid_surface(ts_off; sigma = ELLIPSE_SIGMA)
    surface!(ax3d, Xof, Yof, Zof; color = (col_off, 0.12), shading = NoShading)
    Xed, Yed, Zed = gaussian_ellipsoid_surface(ts_ed;  sigma = ELLIPSE_SIGMA)
    surface!(ax3d, Xed, Yed, Zed; color = (col_edge, 0.12), shading = NoShading)

    axislegend(ax3d; position = :rt, labelsize = 10)
end

Legend(fig1d[1, 3],
    [
        LineElement(color = col_on,   linestyle = :solid, linewidth = 1.2),
        LineElement(color = col_off,  linestyle = :solid, linewidth = 1.2),
        LineElement(color = col_edge, linestyle = :solid, linewidth = 1.0),
    ],
    ["AMOC-on eq. run", "AMOC-off eq. run", "Edge eq. run"],
    "Equilibrium run",
    labelsize = 11,
)

fig1d_path = plotsdir("plasim_equilibrium_scatter_3d.png")
wsave(fig1d_path, fig1d)
@info "Figure 1d saved to: $fig1d_path"

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Bar chart — mean convergence times
# ─────────────────────────────────────────────────────────────────────────────

fig2 = Figure(size = (700, 500))
ax2  = Axis(fig2[1, 1];
    xlabel    = "",
    ylabel    = "Mean convergence time (yr)",
    title     = "Mean convergence times by CO₂ level",
    titlesize = 14,
)

categories  = ["AMOC-on\n285 ppm", "AMOC-off\n285 ppm",
               "AMOC-on\n360 ppm", "AMOC-off\n360 ppm"]
conv_values = [
    summary_285.mean_conv_time_on,
    summary_285.mean_conv_time_off,
    summary_360.mean_conv_time_on,
    summary_360.mean_conv_time_off,
]
bar_colors = [col_on, col_off, col_on, col_off]
bar_alpha   = [1.0, 1.0, 0.55, 0.55]  # darker for 285 ppm

barplot!(ax2, 1:4, conv_values;
    color = [(c, a) for (c, a) in zip(bar_colors, bar_alpha)],
    width = 0.6,
)

# x-axis tick labels
ax2.xticks = (1:4, categories)

# Legend: shade for CO2 level
Legend(fig2[1, 2],
    [PolyElement(color = :gray50), PolyElement(color = :gray80)],
    ["285 ppm", "360 ppm"],
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
    summary_285.mean_dist_on,
    summary_285.mean_dist_off,
    summary_360.mean_dist_on,
    summary_360.mean_dist_off,
]

barplot!(ax3, 1:4, dist_values;
    color = [(c, a) for (c, a) in zip(bar_colors, bar_alpha)],
    width = 0.6,
)

ax3.xticks = (1:4, categories)

Legend(fig3[1, 2],
    [PolyElement(color = :gray50), PolyElement(color = :gray80)],
    ["285 ppm", "360 ppm"],
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
    "summary_285" => summary_285,
)

wsave(datadir("plasim", "resilience_summaries.jld2"), results_to_save)
@info "Results saved to: $(datadir("plasim", "resilience_summaries.jld2"))"

# ─────────────────────────────────────────────────────────────────────────────
# Save key metrics to CSV
# ─────────────────────────────────────────────────────────────────────────────

# Volume of the 1σ ellipsoid in 3D EOF space: V = (4π/3) · √det(C)
ellipsoid_volume(C::Matrix) = (4π / 3) * sqrt(det(Symmetric(C)))

# Mean AMOC strength from an equilibrium NetCDF file
function mean_amoc_strength(filepath::String)
    NCDataset(filepath, "r") do ds
        mean(skipmissing(ds["amoc_strength"][:]))
    end
end

amoc_on_285  = mean_amoc_strength(joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_on.etc.nc"))
amoc_off_285 = mean_amoc_strength(joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_of.etc.nc"))
amoc_on_360  = mean_amoc_strength(joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_on.etc.nc"))
amoc_off_360 = mean_amoc_strength(joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_of.etc.nc"))

ellipse_long_axis_1sigma(C::Matrix) = 2000 * sqrt(maximum(eigvals(Symmetric(C[1:2, 1:2]))))

# local_resilience = 1 / (π × √det(C[1:2,1:2])):
# area of the 1σ Gaussian ellipse in (EOF1, EOF2) space = π × a × b = π × √(λ₁ × λ₂) = π × √det(C₁₂)
local_resilience_ellipse(C::Matrix) = 1.0 / (π * sqrt(det(Symmetric(C[1:2, 1:2]))))

metrics_df = DataFrame(
    co2_ppm            = [285,                        285,                         360,                        360                        ],
    state              = ["AMOC-on",                  "AMOC-off",                  "AMOC-on",                  "AMOC-off"                 ],
    mean_conv_time_yr  = [summary_285.mean_conv_time_on,  summary_285.mean_conv_time_off,  summary_360.mean_conv_time_on,  summary_360.mean_conv_time_off ],
    mean_edge_dist     = [summary_285.mean_dist_on,       summary_285.mean_dist_off,       summary_360.mean_dist_on,       summary_360.mean_dist_off      ],
    ellipsoid_volume_1sigma = [
        ellipsoid_volume(var_on_285.covariance),
        ellipsoid_volume(var_off_285.covariance),
        ellipsoid_volume(var_on_360.covariance),
        ellipsoid_volume(var_off_360.covariance),
    ],
    ellipse_long_axis_1sigma = [
        ellipse_long_axis_1sigma(var_on_285.covariance),
        ellipse_long_axis_1sigma(var_off_285.covariance),
        ellipse_long_axis_1sigma(var_on_360.covariance),
        ellipse_long_axis_1sigma(var_off_360.covariance),
    ],
    local_resilience = [
        local_resilience_ellipse(var_on_285.covariance),
        local_resilience_ellipse(var_off_285.covariance),
        local_resilience_ellipse(var_on_360.covariance),
        local_resilience_ellipse(var_off_360.covariance),
    ],
    mean_amoc_strength_Sv = [amoc_on_285, amoc_off_285, amoc_on_360, amoc_off_360],
)

csv_path = datadir("plasim", "resilience_metrics.csv")
CSV.write(csv_path, metrics_df)
@info "Metrics CSV saved to: $csv_path"

# ─────────────────────────────────────────────────────────────────────────────
# Local variability summary print
# ─────────────────────────────────────────────────────────────────────────────

# Print summary
println("\n=== Local Variability around Equilibria ===\n")
for (lbl, von, voff) in [
    ("285 ppm", var_on_285, var_off_285),
    ("360 ppm", var_on_360, var_off_360),
]
    println("--- $lbl ---")
    for (state, v) in [("AMOC-on", von), ("AMOC-off", voff)]
        println("  $state (n=$(v.n_samples) yr):")
        println("    Std per dim        : $(round.(v.std_per_dim;    digits=4))")
        println("    Total variance     : $(round(v.total_variance;  digits=6))")
        println("    Dominant variance  : $(round(v.dominant_variance; digits=6))")
        println("    Mean dist from μ   : $(round(v.mean_dist;       digits=4))")
        println("    Lag-1 autocorr     : $(round.(v.lag1_autocorr;  digits=3))")
        println("    Autocorr time (yr) : $(round.(v.ac_times;       digits=1))")
    end
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Local variability comparison (on vs off, 285 vs 360 ppm)
# ─────────────────────────────────────────────────────────────────────────────

var_measure_labels = [
    "Total variance",
    "Dominant variance",
    "Mean dist from μ",
    "Mean lag-1 AC",
    "Mean AC time (yr)",
]
var_measure_getters = [
    v -> v.total_variance,
    v -> v.dominant_variance,
    v -> v.mean_dist,
    v -> v.mean_lag1_autocorr,
    v -> v.mean_ac_time,
]

fig4 = Figure(size = (1000, 300 * length(var_measure_labels)))

for (row, (mlabel, getter)) in enumerate(zip(var_measure_labels, var_measure_getters))
    ax4 = Axis(fig4[row, 1];
        ylabel             = mlabel,
        xticklabelsvisible = row == length(var_measure_labels),
        xlabel             = row == length(var_measure_labels) ? "" : "",
    )

    vals = [getter(var_on_285), getter(var_off_285), getter(var_on_360), getter(var_off_360)]
    barplot!(ax4, 1:4, vals;
        color = [(c, a) for (c, a) in zip(bar_colors, bar_alpha)],
        width = 0.6,
    )
    ax4.xticks = (1:4, categories)
end

Legend(fig4[1, 2],
    [PolyElement(color = col_on), PolyElement(color = col_off)],
    ["AMOC-on", "AMOC-off"],
    "State",
)
Legend(fig4[2, 2],
    [PolyElement(color = (col_on, 1.0)), PolyElement(color = (col_on, 0.55))],
    ["285 ppm", "360 ppm"],
    "CO₂ level",
)

fig4_path = plotsdir("plasim_local_variability.png")
wsave(fig4_path, fig4)
@info "Figure 4 saved to: $fig4_path"

# Display figures if running interactively
display(fig1)
display(fig1b)
display(fig1c)
display(fig1d)
display(fig2)
display(fig3)
display(fig4)
