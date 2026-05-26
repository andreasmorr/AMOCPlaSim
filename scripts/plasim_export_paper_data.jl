"""
    plasim_export_paper_data.jl

Export PlaSim edge-track data to CSV files for Python plotting.

Outputs (all in data/plasim/paper/):
    trajectories_{co2_label}.csv   — traj_id, time, x1, x2, amoc_strength, label
    equilibria_{co2_label}.csv     — state, time, x1, x2, amoc_strength
    ellipses_{co2_label}.csv       — state, x, y
    state_means_{co2_label}.csv    — state, x1, x2

Run from the project root:
    julia --project scripts/plasim_export_paper_data.jl
"""

using DrWatson
@quickactivate "AMOCResilience"

include(srcdir("plasim_utils.jl"))

using NCDatasets
using DataFrames
using Statistics
using LinearAlgebra
using CSV

# ─────────────────────────────────────────────────────────────────────────────
# Configuration (mirrors plasim_edge_analysis.jl exactly)
# ─────────────────────────────────────────────────────────────────────────────

const DATA_DIR                = datadir("plasim")
const CO2_LABEL_CURRENT       = "360ppm"
const CO2_LABEL_PREINDUSTRIAL = "285ppm"
const N_FILES_360             = 38
const N_FILES_285             = 37
const VARIABLE_NAMES          = ["redu1", "redu2", "redu3"]
const N_DIMS                  = length(VARIABLE_NAMES)
const EPSILON_ON              = 0.1
const EPSILON_OFF             = 0.2
const FINAL_FRACTION          = 0.1
const ELLIPSE_SIGMA           = 3

# ─────────────────────────────────────────────────────────────────────────────
# File pattern helpers
# ─────────────────────────────────────────────────────────────────────────────

file_pattern_360(i) = "plasimelancholia_$(CO2_LABEL_CURRENT)_edgetrack_iter$(lpad(i-1, 3, '0')).etc.nc"
file_pattern_285(i) = "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_edgetrack_itx$(lpad(i-1, 3, '0')).etc.nc"

# ─────────────────────────────────────────────────────────────────────────────
# Helper: try to read AMOC strength from a NetCDF file
# ─────────────────────────────────────────────────────────────────────────────

function try_load_amoc_strength(filepath::String)
    try
        NCDataset(filepath, "r") do ds
            haskey(ds, "amoc_strength") || return nothing
            return Float64.(coalesce.(ds["amoc_strength"][:], NaN))
        end
    catch
        return nothing
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Export trajectories for one CO2 level
# ─────────────────────────────────────────────────────────────────────────────

function export_trajectories(
    df::DataFrame,
    summary,
    co2_label::String,
    n_files::Int,
    file_pattern::Function,
)
    ids    = summary.trajectory_ids
    labels = summary.attractor_labels
    conv   = summary.converged

    # Build per-trajectory time-series rows
    rows = DataFrame(
        traj_id       = Int[],
        time          = Float64[],
        x1            = Float64[],
        x2            = Float64[],
        amoc_strength = Float64[],
        label         = String[],
    )

    for (j, tid) in enumerate(ids)
        conv[tid] || continue  # only converged trajectories

        sub = sort(filter(r -> r.trajectory_id == tid, df), :time)
        state_label = labels[j] == 1 ? "on" : "off"

        # Try to load AMOC strength from the originating NetCDF file
        # file_id and track_id are stored in the DataFrame
        file_id  = sub[1, :file_id]
        track_id = sub[1, :track_id]
        fname    = joinpath(DATA_DIR, file_pattern(file_id))
        amoc_raw = try_load_amoc_strength(fname)

        for t in 1:nrow(sub)
            # Determine the amoc_strength for this row's time index
            amoc_val = NaN
            if amoc_raw !== nothing
                # amoc_raw may be shaped (n_time,) or (n_tracks, n_time) etc.
                # Use track_id offset if applicable; guard with try/catch
                try
                    if ndims(amoc_raw) == 1
                        amoc_val = t <= length(amoc_raw) ? amoc_raw[t] : NaN
                    elseif ndims(amoc_raw) == 2
                        # Try (track, time) layout
                        amoc_val = amoc_raw[track_id, t]
                    end
                catch
                    amoc_val = NaN
                end
            end

            push!(rows, (
                tid,
                sub[t, :time],
                sub[t, :x1],
                sub[t, :x2],
                amoc_val,
                state_label,
            ))
        end
    end

    out_path = datadir("plasim", "paper", "trajectories_$(co2_label).csv")
    CSV.write(out_path, rows)
    @info "  Trajectories saved: $out_path  ($(nrow(rows)) rows)"
end

# ─────────────────────────────────────────────────────────────────────────────
# Export equilibrium run time-series
# ─────────────────────────────────────────────────────────────────────────────

function export_equilibria(
    ts_on::Matrix,
    ts_off::Matrix,
    ts_ed::Matrix,
    co2_label::String,
)
    # Try loading AMOC strength from the equilibrium NetCDF files
    amoc_on_arr  = try_load_amoc_strength(
        joinpath(DATA_DIR, "plasimelancholia_$(co2_label)_on.etc.nc"))
    amoc_off_arr = try_load_amoc_strength(
        joinpath(DATA_DIR, "plasimelancholia_$(co2_label)_of.etc.nc"))
    amoc_ed_arr  = try_load_amoc_strength(
        joinpath(DATA_DIR, "plasimelancholia_$(co2_label)_ed.etc.nc"))

    rows = DataFrame(
        state         = String[],
        time          = Float64[],
        x1            = Float64[],
        x2            = Float64[],
        amoc_strength = Float64[],
    )

    for (state_name, ts, amoc_arr) in [
            ("on",  ts_on,  amoc_on_arr),
            ("off", ts_off, amoc_off_arr),
            ("edge", ts_ed, amoc_ed_arr),
    ]
        n_time = size(ts, 1)
        for t in 1:n_time
            any(isnan, ts[t, :]) && continue
            amoc_val = (amoc_arr !== nothing && t <= length(amoc_arr)) ?
                       amoc_arr[t] : NaN
            push!(rows, (state_name, Float64(t), ts[t, 1], ts[t, 2], amoc_val))
        end
    end

    out_path = datadir("plasim", "paper", "equilibria_$(co2_label).csv")
    CSV.write(out_path, rows)
    @info "  Equilibria saved: $out_path  ($(nrow(rows)) rows)"
end

# ─────────────────────────────────────────────────────────────────────────────
# Export Gaussian ellipse points and state means
# ─────────────────────────────────────────────────────────────────────────────

function gaussian_ellipse_points(ts::Matrix, dim1::Int = 1, dim2::Int = 2;
                                  n::Int = 120, sigma::Real = 1)
    valid = [!any(isnan, ts[t, :]) for t in axes(ts, 1)]
    data  = ts[valid, [dim1, dim2]]
    μ     = vec(mean(data, dims = 1))
    C     = cov(data)
    vals, vecs = eigen(Symmetric(C))
    vals  = max.(vals, 0.0)
    θ     = LinRange(0, 2π, n + 1)
    pts   = vecs * (sigma .* sqrt.(vals) .* [cos.(θ)'; sin.(θ)'])
    return μ[1] .+ pts[1, :], μ[2] .+ pts[2, :]
end

function export_ellipses(
    ts_on::Matrix,
    ts_off::Matrix,
    ts_ed::Matrix,
    co2_label::String,
)
    ellipse_rows = DataFrame(state = String[], x = Float64[], y = Float64[])

    for (state_name, ts) in [("on", ts_on), ("off", ts_off), ("edge", ts_ed)]
        ex, ey = gaussian_ellipse_points(ts, 1, 2; sigma = ELLIPSE_SIGMA)
        for (xi, yi) in zip(ex, ey)
            push!(ellipse_rows, (state_name, xi, yi))
        end
    end

    ell_path = datadir("plasim", "paper", "ellipses_$(co2_label).csv")
    CSV.write(ell_path, ellipse_rows)
    @info "  Ellipses saved: $ell_path"
end

function export_state_means(
    ts_on::Matrix,
    ts_off::Matrix,
    ts_ed::Matrix,
    co2_label::String,
)
    function ts_mean2d(ts)
        valid = [!any(isnan, ts[t, :]) for t in axes(ts, 1)]
        data  = ts[valid, :]
        return vec(mean(data, dims = 1))
    end

    mu_on  = ts_mean2d(ts_on)
    mu_off = ts_mean2d(ts_off)
    mu_ed  = ts_mean2d(ts_ed)

    means_df = DataFrame(
        state = ["on",     "off",      "edge"   ],
        x1    = [mu_on[1], mu_off[1],  mu_ed[1] ],
        x2    = [mu_on[2], mu_off[2],  mu_ed[2] ],
    )

    means_path = datadir("plasim", "paper", "state_means_$(co2_label).csv")
    CSV.write(means_path, means_df)
    @info "  State means saved: $means_path"
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

function main()
    mkpath(datadir("plasim", "paper"))

    # ── Load trajectories ─────────────────────────────────────────────────────
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

    # ── Load equilibrium states ───────────────────────────────────────────────
    @info "Loading equilibrium states for 285 ppm..."
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

    ts_on_285  = load_plasim_state_timeseries(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_on.etc.nc"),
        VARIABLE_NAMES)
    ts_off_285 = load_plasim_state_timeseries(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_of.etc.nc"),
        VARIABLE_NAMES)
    ts_ed_285  = load_plasim_state_timeseries(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_ed.etc.nc"),
        VARIABLE_NAMES)

    @info "Loading equilibrium states for 360 ppm..."
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

    ts_on_360  = load_plasim_state_timeseries(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_on.etc.nc"),
        VARIABLE_NAMES)
    ts_off_360 = load_plasim_state_timeseries(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_of.etc.nc"),
        VARIABLE_NAMES)
    ts_ed_360  = load_plasim_state_timeseries(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_ed.etc.nc"),
        VARIABLE_NAMES)

    # ── Local variability ─────────────────────────────────────────────────────
    @info "Computing local variability for 285 ppm..."
    var_on_285  = compute_local_variability(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_on.etc.nc"),
        VARIABLE_NAMES)
    var_off_285 = compute_local_variability(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_of.etc.nc"),
        VARIABLE_NAMES)
    var_ed_285  = compute_local_variability(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_PREINDUSTRIAL)_ed.etc.nc"),
        VARIABLE_NAMES)

    @info "Computing local variability for 360 ppm..."
    var_on_360  = compute_local_variability(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_on.etc.nc"),
        VARIABLE_NAMES)
    var_off_360 = compute_local_variability(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_of.etc.nc"),
        VARIABLE_NAMES)
    var_ed_360  = compute_local_variability(
        joinpath(DATA_DIR, "plasimelancholia_$(CO2_LABEL_CURRENT)_ed.etc.nc"),
        VARIABLE_NAMES)

    # ── Resilience summaries ──────────────────────────────────────────────────
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

    # ── Export CSVs ───────────────────────────────────────────────────────────
    @info "Exporting 285 ppm data..."
    export_trajectories(df_285, summary_285, CO2_LABEL_PREINDUSTRIAL,
                        N_FILES_285, file_pattern_285)
    export_equilibria(ts_on_285, ts_off_285, ts_ed_285, CO2_LABEL_PREINDUSTRIAL)
    export_ellipses(ts_on_285, ts_off_285, ts_ed_285, CO2_LABEL_PREINDUSTRIAL)
    export_state_means(ts_on_285, ts_off_285, ts_ed_285, CO2_LABEL_PREINDUSTRIAL)

    @info "Exporting 360 ppm data..."
    export_trajectories(df_360, summary_360, CO2_LABEL_CURRENT,
                        N_FILES_360, file_pattern_360)
    export_equilibria(ts_on_360, ts_off_360, ts_ed_360, CO2_LABEL_CURRENT)
    export_ellipses(ts_on_360, ts_off_360, ts_ed_360, CO2_LABEL_CURRENT)
    export_state_means(ts_on_360, ts_off_360, ts_ed_360, CO2_LABEL_CURRENT)

    @info "All PlaSim paper data exported to $(datadir("plasim", "paper"))"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
