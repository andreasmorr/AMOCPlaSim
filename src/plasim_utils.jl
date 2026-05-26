"""
    plasim_utils.jl

Utilities for loading and analyzing PlaSim edge-state trajectory data.

PlaSim (Planet Simulator) is a general circulation model. The edge-state
(saddle) trajectories are computed by bisecting between initial conditions
that converge to AMOC-on vs. AMOC-off attractors. The data is stored as
NetCDF files with EOF-reduced dimensions (redu1, redu2, redu3).

Each NetCDF file contains 2 trajectories (tracks), one converging to AMOC-on
and one to AMOC-off.

Key workflow:
  1. `load_plasim_trajectories` — load all NetCDF files into a DataFrame
  2. `classify_trajectories`    — label each trajectory AMOC-on or AMOC-off
  3. `estimate_attractors`      — estimate attractor positions from final states
  4. `compute_convergence_times`           — time to converge to attractor
  5. `compute_edge_to_attractor_distances` — distance from edge state to attractor
  6. `plasim_resilience_summary`           — convenience wrapper for all of the above

Local variability around equilibria:
  7. `load_plasim_state_timeseries` — load full time series from an equilibrium file
  8. `compute_local_variability`    — variance, autocorrelation, and distance measures
"""

using NCDatasets
using DataFrames
using Statistics
using LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# Gaussian-ellipse helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _in_ellipse(x, μ, C_inv; sigma=1) → Bool

Return `true` if point `x` lies inside the nσ Gaussian ellipse/ellipsoid
defined by mean `μ` and precision matrix `C_inv` (inverse covariance).

The criterion is the Mahalanobis distance: (x-μ)ᵀ C⁻¹ (x-μ) ≤ sigma².
"""
function _in_ellipse(x::AbstractVector, μ::AbstractVector, C_inv::AbstractMatrix;
                     sigma::Real = 1)
    d = x .- μ
    return dot(d, C_inv * d) ≤ Float64(sigma)^2
end

"""
    ellipse_to_ellipse_distance(μ1, C1, μ2, C2; sigma=1) → Float64

Distance between the surfaces of two nσ Gaussian ellipsoids.

Computed as: max(0, ‖μ2-μ1‖ − σ·r1(v) − σ·r2(v)), where v is the unit
vector from μ1 to μ2 and r(v) = 1/√(vᵀ C⁻¹ v) is the 1σ radius along v.
Returns 0 when the ellipsoids overlap.
"""
function ellipse_to_ellipse_distance(μ1::Vector{Float64}, C1::Matrix{Float64},
                                      μ2::Vector{Float64}, C2::Matrix{Float64};
                                      sigma::Real = 1)
    dir       = μ2 .- μ1
    d_centers = norm(dir)
    d_centers == 0.0 && return 0.0
    v  = dir ./ d_centers
    r1 = Float64(sigma) / sqrt(dot(v, Symmetric(C1) \ v))
    r2 = Float64(sigma) / sqrt(dot(v, Symmetric(C2) \ v))
    return max(0.0, d_centers - r1 - r2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_plasim_state_mean(filepath, variable_names) → Vector{Float64}

Load a converged-state NetCDF file (e.g. `plasimelancholia_285ppm_on.etc.nc`)
and return the time-mean of each variable as a state vector.

These files contain trajectories that have already converged to the vicinity of
an equilibrium (AMOC-on, AMOC-off, or edge state), so their time-mean is a
better estimate of the attractor/edge position than the final states of
bisection trajectories.

# Arguments
- `filepath`:       full path to the NetCDF file
- `variable_names`: variable names to read (same order as used elsewhere, e.g.
                    `["redu1", "redu2", "redu3"]`)

# Returns
`Vector{Float64}` of length `length(variable_names)` with the time-mean state.
"""
function load_plasim_state_mean(filepath::String, variable_names::Vector{String})
    means = zeros(Float64, length(variable_names))
    NCDataset(filepath, "r") do ds
        for (k, vname) in enumerate(variable_names)
            data = Float64.(coalesce.(ds[vname][:], NaN))
            valid = filter(!isnan, data)
            means[k] = isempty(valid) ? NaN : mean(valid)
        end
    end
    return means
end

"""
    load_plasim_trajectories(;
        co2_label::String,
        data_dir::String,
        n_files::Int,
        variable_names::Vector{String} = ["redu1", "redu2", "redu3"],
        file_pattern::Function = i -> "plasimelancholia_\$(co2_label)_edgetrack_iter\$(i).etc.nc"
    ) → DataFrame

Load all PlaSim edge-track NetCDF files for a given CO2 level into a DataFrame.

Each file contains two tracks (trajectories). The function reads `variable_names`
variables from each file and concatenates them. The resulting DataFrame has columns:
- `trajectory_id::Int`   — unique integer across all files and tracks
- `file_id::Int`         — which file (1-indexed)
- `track_id::Int`        — which track within the file (1 or 2)
- `time::Float64`        — time step index (1-based)
- `x1::Float64`          — first EOF component (redu1)
- `x2::Float64`          — second EOF component (redu2)
- `x3::Float64`          — third EOF component (redu3)  [if present]
- (additional xi columns for each variable)

# Arguments
- `co2_label`: string label used in file names, e.g. `"360ppm"` or `"720ppm"`
- `data_dir`:  directory containing the NetCDF files
- `n_files`:   number of files to attempt loading (files that don't exist are skipped)
- `variable_names`: names of NetCDF variables to read (default: redu1, redu2, redu3)
- `file_pattern`: function `i -> filename` mapping file index to filename

# Returns
DataFrame with one row per (trajectory, time) pair.
"""
function load_plasim_trajectories(;
    co2_label::String,
    data_dir::String,
    n_files::Int,
    variable_names::Vector{String} = ["redu1", "redu2", "redu3"],
    file_pattern::Function = i -> "plasimelancholia_$(co2_label)_edgetrack_iter$(i).etc.nc"
)
    n_dims = length(variable_names)
    rows = Vector{NamedTuple}()
    trajectory_counter = 0

    for i in 1:n_files
        fname = joinpath(data_dir, file_pattern(i))
        if !isfile(fname)
            @warn "File not found, skipping: $fname"
            continue
        end

        NCDataset(fname, "r") do ds
            # Determine number of time steps and tracks
            # Convention: dims are (time, track) or (track, time) — check both
            dim_names = collect(keys(ds.dim))

            # Read first variable to determine shape
            first_var = ds[variable_names[1]][:, :]  # (time, track) or (track, time)
            n_time, n_tracks = size(first_var)

            # If n_tracks > n_time it's likely transposed; swap
            # Heuristic: tracks are expected to be 2, time should be many steps
            if n_tracks > n_time
                # data is (track, time) — transpose
                for track in 1:n_time
                    trajectory_counter += 1
                    tid = trajectory_counter
                    data_cols = Dict{Symbol, Vector{Float64}}()
                    for (k, vname) in enumerate(variable_names)
                        raw = ds[vname][:, :]   # (track_dim, time_dim)
                        data_cols[Symbol("x$k")] = Float64.(coalesce.(raw[track, :], NaN))
                    end
                    t_steps = 1:n_tracks  # n_tracks is the time dimension here
                    for t in t_steps
                        push!(rows, (
                            trajectory_id = tid,
                            file_id       = i,
                            track_id      = track,
                            time          = Float64(t),
                            (Symbol("x$k") => data_cols[Symbol("x$k")][t]
                             for k in 1:n_dims)...
                        ))
                    end
                end
            else
                # data is (time, track) — standard layout
                for track in 1:n_tracks
                    trajectory_counter += 1
                    tid = trajectory_counter
                    data_cols = Dict{Symbol, Vector{Float64}}()
                    for (k, vname) in enumerate(variable_names)
                        raw = ds[vname][:, :]   # (time_dim, track_dim)
                        data_cols[Symbol("x$k")] = Float64.(coalesce.(raw[:, track], NaN))
                    end
                    for t in 1:n_time
                        push!(rows, (
                            trajectory_id = tid,
                            file_id       = i,
                            track_id      = track,
                            time          = Float64(t),
                            (Symbol("x$k") => data_cols[Symbol("x$k")][t]
                             for k in 1:n_dims)...
                        ))
                    end
                end
            end
        end
    end

    if isempty(rows)
        error("No data loaded. Check data_dir and file_pattern.")
    end

    df = DataFrame(rows)
    @info "Loaded $(trajectory_counter) trajectories from $(data_dir)"
    return df
end

# ─────────────────────────────────────────────────────────────────────────────
# Trajectory classification
# ─────────────────────────────────────────────────────────────────────────────

"""
    _get_final_states(df, n_dims; final_fraction=0.1) → (ids, Matrix)

Internal helper. Returns a vector of trajectory IDs and a matrix of mean
final states (n_trajectories × n_dims). The final state is averaged over the
last `final_fraction` of each trajectory.
"""
function _get_final_states(df::DataFrame, n_dims::Int; final_fraction::Float64 = 0.1)
    ids = sort(unique(df.trajectory_id))
    final_states = zeros(Float64, length(ids), n_dims)

    for (j, tid) in enumerate(ids)
        traj = filter(r -> r.trajectory_id == tid, df)
        sort!(traj, :time)
        # Drop rows where any EOF dimension is NaN (NetCDF fill values)
        valid_mask = [!any(isnan(traj[t, Symbol("x$k")]) for k in 1:n_dims) for t in 1:nrow(traj)]
        traj = traj[valid_mask, :]
        n_pts = nrow(traj)
        n_final = max(1, round(Int, final_fraction * n_pts))
        tail = traj[(end - n_final + 1):end, :]
        for k in 1:n_dims
            final_states[j, k] = mean(tail[:, Symbol("x$k")])
        end
    end

    return ids, final_states
end

"""
    classify_trajectories(df, n_dims; final_fraction=0.1) → Vector{Int}

Classify each trajectory as ending at the AMOC-on (label=1) or AMOC-off
(label=2) attractor, based on the first EOF component (x1 = redu1).

The AMOC-on state has a lower value of redu1, and the AMOC-off state has a
higher value. We split by the median of final redu1 values.

Returns a vector of length `n_trajectories` with values 1 (AMOC-on) or 2
(AMOC-off), in the same order as `sort(unique(df.trajectory_id))`.
"""
function classify_trajectories(df::DataFrame, n_dims::Int;
                                final_fraction::Float64 = 0.1)
    ids, final_states = _get_final_states(df, n_dims; final_fraction)

    # Use first EOF (redu1) to separate AMOC-on from AMOC-off
    # Convention: lower redu1 = AMOC-on
    redu1_vals = final_states[:, 1]
    threshold  = median(redu1_vals)

    labels = [v < threshold ? 1 : 2 for v in redu1_vals]
    return labels
end

# ─────────────────────────────────────────────────────────────────────────────
# Attractor estimation
# ─────────────────────────────────────────────────────────────────────────────

"""
    estimate_attractors(df, n_dims; final_fraction=0.1) → Dict{Int, Vector{Float64}}

Estimate the positions of the AMOC-on and AMOC-off attractors by averaging
the final states of trajectories that converge to each one.

The split between on/off is done by the median of final redu1 values (first EOF).

Returns a Dict:
- `1 => mean_state_on`   (AMOC-on attractor, length n_dims)
- `2 => mean_state_off`  (AMOC-off attractor, length n_dims)
"""
function estimate_attractors(df::DataFrame, n_dims::Int;
                              final_fraction::Float64 = 0.1)
    ids, final_states = _get_final_states(df, n_dims; final_fraction)
    labels = classify_trajectories(df, n_dims; final_fraction)

    on_idx  = findall(==(1), labels)
    off_idx = findall(==(2), labels)

    mean_on  = vec(mean(final_states[on_idx,  :]; dims = 1))
    mean_off = vec(mean(final_states[off_idx, :]; dims = 1))

    attractors = Dict{Int, Vector{Float64}}(
        1 => mean_on,
        2 => mean_off
    )

    @info "Estimated attractors: $(length(on_idx)) on-trajectories, $(length(off_idx)) off-trajectories"
    return attractors
end

# ─────────────────────────────────────────────────────────────────────────────
# Convergence times
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_convergence_times(df, n_dims, attractors; ...) → Dict{Int, Float64}

For each trajectory, compute the convergence time as the span between:
  - the first step when the trajectory exits the edge-state 1σ ellipsoid, and
  - the first step when it enters (and stays inside) the target attractor's 1σ ellipsoid.

This ellipse-based definition requires `cov_on`, `cov_off`, `edge_state`, and
`edge_cov` to be supplied.  When any of these is `nothing` the function falls
back to the legacy behaviour: time to first step where Euclidean distance to
the attractor mean stays ≤ ε for all remaining steps.

# Arguments
- `df`, `n_dims`, `attractors`: as before
- `ε_on`, `ε_off`:   fallback ε-ball radii (used only when covariances are absent)
- `cov_on`:    covariance matrix of the AMOC-on equilibrium run (n_dims × n_dims)
- `cov_off`:   covariance matrix of the AMOC-off equilibrium run
- `edge_state`: mean position of the edge equilibrium run
- `edge_cov`:  covariance matrix of the edge equilibrium run

# Returns
Dict mapping `trajectory_id => convergence_time` (NaN if not converged).
When ellipse mode is active, the time is measured in the same units as
`df.time` (steps from edge-exit to attractor-entry).
"""
function compute_convergence_times(df::DataFrame, n_dims::Int,
                                    attractors::Dict{Int, Vector{Float64}};
                                    ε_on::Float64 = 0.05, ε_off::Float64 = 0.05,
                                    cov_on::Union{Nothing, Matrix{Float64}}    = nothing,
                                    cov_off::Union{Nothing, Matrix{Float64}}   = nothing,
                                    edge_state::Union{Nothing, Vector{Float64}} = nothing,
                                    edge_cov::Union{Nothing, Matrix{Float64}}  = nothing,
                                    sigma::Real = 1,
                                    check_dims::Int = 2)
    use_ellipses = (cov_on !== nothing && cov_off !== nothing &&
                    edge_state !== nothing && edge_cov !== nothing)

    labels = classify_trajectories(df, n_dims)
    ids    = sort(unique(df.trajectory_id))
    result = Dict{Int, Float64}()

    # Slice covariances and means to the first check_dims dimensions
    d = check_dims
    C_on_inv  = use_ellipses ? Matrix(inv(Symmetric(cov_on[1:d, 1:d])))  : nothing
    C_off_inv = use_ellipses ? Matrix(inv(Symmetric(cov_off[1:d, 1:d]))) : nothing
    C_ed_inv  = use_ellipses ? Matrix(inv(Symmetric(edge_cov[1:d, 1:d]))) : nothing

    for (j, tid) in enumerate(ids)
        lbl              = labels[j]
        target_attractor = attractors[lbl][1:d]

        traj   = sort(filter(r -> r.trajectory_id == tid, df), :time)
        times  = traj.time
        n_pts  = nrow(traj)
        states = hcat([traj[:, Symbol("x$k")] for k in 1:d]...)

        if use_ellipses
            C_att_inv = lbl == 1 ? C_on_inv : C_off_inv
            ed_mean   = edge_state[1:d]

            # Last step inside the edge ellipse (final departure point)
            # If the trajectory never entered the edge ellipse, exclude it (NaN)
            t_exit_idx = findlast(t -> _in_ellipse(states[t, :], ed_mean, C_ed_inv; sigma), 1:n_pts)
            if t_exit_idx === nothing
                result[tid] = NaN
                continue
            end

            # First step inside the attractor ellipse after t_exit
            conv_time = NaN
            for t in t_exit_idx:n_pts
                if _in_ellipse(states[t, :], target_attractor, C_att_inv; sigma)
                    conv_time = times[t] - times[t_exit_idx]
                    break
                end
            end
            result[tid] = conv_time
        else
            ε     = lbl == 1 ? ε_on : ε_off
            dists = [norm(states[t, :] .- target_attractor) for t in 1:n_pts]
            conv_time = NaN
            for t in 1:n_pts
                if all(dists[t:end] .<= ε)
                    conv_time = times[t]
                    break
                end
            end
            result[tid] = conv_time
        end
    end

    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Edge-to-attractor distances
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_edge_to_attractor_distances(df, n_dims, attractors; ...) → Dict{Int, Float64}

Compute the distance from the edge state to each trajectory's target attractor.

When `cov_on`, `cov_off`, `edge_state`, and `edge_cov` are all provided the
distance is the **ellipse-to-ellipse gap**: the distance between the surfaces
of the two 1σ Gaussian ellipsoids (zero when they overlap).

When covariances are absent the function falls back to the legacy Euclidean
distance between the edge-state mean (or trajectory initial state) and the
attractor mean.

# Returns
Dict mapping `trajectory_id => distance`.  In ellipse mode all trajectories
converging to the same attractor share the same value.
"""
function compute_edge_to_attractor_distances(df::DataFrame, n_dims::Int,
                                              attractors::Dict{Int, Vector{Float64}};
                                              edge_state::Union{Nothing, Vector{Float64}} = nothing,
                                              cov_on::Union{Nothing, Matrix{Float64}}    = nothing,
                                              cov_off::Union{Nothing, Matrix{Float64}}   = nothing,
                                              edge_cov::Union{Nothing, Matrix{Float64}}  = nothing,
                                              sigma::Real = 1,
                                              check_dims::Int = 2)
    use_ellipses = (cov_on !== nothing && cov_off !== nothing &&
                    edge_state !== nothing && edge_cov !== nothing)

    d      = check_dims
    labels = classify_trajectories(df, n_dims)
    ids    = sort(unique(df.trajectory_id))
    result = Dict{Int, Float64}()

    for (j, tid) in enumerate(ids)
        lbl              = labels[j]
        target_attractor = attractors[lbl]

        if use_ellipses
            C_att = lbl == 1 ? cov_on[1:d, 1:d] : cov_off[1:d, 1:d]
            result[tid] = ellipse_to_ellipse_distance(
                edge_state[1:d], edge_cov[1:d, 1:d],
                target_attractor[1:d], C_att; sigma)
        else
            es = if edge_state !== nothing
                edge_state
            else
                traj = sort(filter(r -> r.trajectory_id == tid, df), :time)
                [traj[1, Symbol("x$k")] for k in 1:n_dims]
            end
            result[tid] = norm(es .- target_attractor)
        end
    end

    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

"""
    plasim_resilience_summary(df, n_dims;
                               final_fraction=0.1, ε_on=0.05, ε_off=0.05,
                               attractors=nothing,
                               edge_state=nothing) → NamedTuple

Convenience wrapper that computes all PlaSim resilience diagnostics.

Calls `classify_trajectories`, `estimate_attractors`, `compute_convergence_times`,
and `compute_edge_to_attractor_distances`, then aggregates the results.

When `attractors` is supplied (a `Dict{Int,Vector{Float64}}` with keys 1=on,
2=off) it is used directly instead of being estimated from the final states of
edgetrack trajectories.  Supply this from the converged `_on.etc.nc` /
`_off.etc.nc` files for more accurate attractor positions.

When `edge_state` is supplied (a `Vector{Float64}`) it is used as the shared
edge-state position for all distance calculations instead of the per-trajectory
initial states.  Supply this from the converged `_ed.etc.nc` file.

# Returns
A NamedTuple with fields:
- `attractor_labels::Vector{Int}`      — per-trajectory label (1=on, 2=off)
- `trajectory_ids::Vector{Int}`        — corresponding trajectory IDs
- `attractors::Dict{Int,Vector{Float64}}` — attractor positions used
- `convergence_times::Dict{Int,Float64}`  — per-trajectory convergence times (NaN = not converged)
- `converged::Dict{Int,Bool}`             — whether each trajectory converged
- `edge_distances::Dict{Int,Float64}`     — edge-to-attractor distances
- `n_on::Int`                           — number of AMOC-on trajectories
- `n_off::Int`                          — number of AMOC-off trajectories
- `mean_conv_time_on::Float64`          — mean convergence time for AMOC-on (NaN if none converged)
- `mean_conv_time_off::Float64`         — mean convergence time for AMOC-off (NaN if none converged)
- `n_converged_on::Int`                 — number of AMOC-on trajectories that converged within ε_on
- `n_converged_off::Int`                — number of AMOC-off trajectories that converged within ε_off
- `mean_dist_on::Float64`               — mean edge distance for AMOC-on
- `mean_dist_off::Float64`              — mean edge distance for AMOC-off
- `edge_state::Union{Nothing,Vector{Float64}}` — edge-state position used
"""
function plasim_resilience_summary(df::DataFrame, n_dims::Int;
                                    final_fraction::Float64 = 0.1,
                                    ε_on::Float64  = 0.05,
                                    ε_off::Float64 = 0.05,
                                    attractors::Union{Nothing, Dict{Int, Vector{Float64}}} = nothing,
                                    edge_state::Union{Nothing, Vector{Float64}} = nothing,
                                    cov_on::Union{Nothing, Matrix{Float64}}    = nothing,
                                    cov_off::Union{Nothing, Matrix{Float64}}   = nothing,
                                    edge_cov::Union{Nothing, Matrix{Float64}}  = nothing,
                                    sigma::Real = 1,
                                    check_dims::Int = 2)
    ids    = sort(unique(df.trajectory_id))
    labels = classify_trajectories(df, n_dims; final_fraction)
    attrs  = attractors !== nothing ? attractors : estimate_attractors(df, n_dims; final_fraction)
    conv_t = compute_convergence_times(df, n_dims, attrs;
                 ε_on, ε_off, cov_on, cov_off, edge_state, edge_cov, sigma, check_dims)
    edist  = compute_edge_to_attractor_distances(df, n_dims, attrs;
                 edge_state, cov_on, cov_off, edge_cov, sigma, check_dims)

    on_idx  = findall(==(1), labels)
    off_idx = findall(==(2), labels)

    on_ids  = ids[on_idx]
    off_ids = ids[off_idx]

    n_on  = length(on_idx)
    n_off = length(off_idx)

    # Bool lookup: did each trajectory converge?
    converged = Dict{Int, Bool}(tid => !isnan(conv_t[tid]) for tid in ids)

    # Only include trajectories that actually converged (conv_time != NaN)
    on_conv_times  = filter(!isnan, [conv_t[tid] for tid in on_ids])
    off_conv_times = filter(!isnan, [conv_t[tid] for tid in off_ids])

    mean_conv_on  = isempty(on_conv_times)  ? NaN : mean(on_conv_times)
    mean_conv_off = isempty(off_conv_times) ? NaN : mean(off_conv_times)

    mean_dist_on  = n_on  > 0 ? mean(edist[tid] for tid in on_ids)  : NaN
    mean_dist_off = n_off > 0 ? mean(edist[tid] for tid in off_ids) : NaN

    return (
        attractor_labels   = labels,
        trajectory_ids     = ids,
        attractors         = attrs,
        convergence_times  = conv_t,
        converged          = converged,
        edge_distances     = edist,
        n_on               = n_on,
        n_off              = n_off,
        n_converged_on     = length(on_conv_times),
        n_converged_off    = length(off_conv_times),
        mean_conv_time_on  = mean_conv_on,
        mean_conv_time_off = mean_conv_off,
        mean_dist_on       = mean_dist_on,
        mean_dist_off      = mean_dist_off,
        edge_state         = edge_state,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Local variability around equilibria
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_plasim_state_timeseries(filepath, variable_names) → Matrix{Float64}

Load a converged-state NetCDF file and return the full time series as a
matrix of shape `(n_time × n_dims)`.  NaN rows (NetCDF fill values) are
kept so that the caller can decide how to handle them.
"""
function load_plasim_state_timeseries(filepath::String, variable_names::Vector{String})
    NCDataset(filepath, "r") do ds
        n_time = length(ds["time"])
        mat = Matrix{Float64}(undef, n_time, length(variable_names))
        for (k, vname) in enumerate(variable_names)
            mat[:, k] = Float64.(coalesce.(ds[vname][:], NaN))
        end
        return mat
    end
end

"""
    _integrated_ac_time(x) → Float64

Estimate the integrated autocorrelation time of a zero-mean time series `x`
as  τ = 0.5 + Σ_{l=1}^{L} ρ(l), where the sum is truncated at the first
lag where the autocorrelation function drops below zero.  This is the
standard estimator used in MCMC diagnostics and CSD analysis.
"""
function _integrated_ac_time(x::AbstractVector{<:Real})
    n   = length(x)
    τ   = 0.5
    for l in 1:min(n ÷ 2, 500)
        ρ = cor(x[1:n-l], x[l+1:n])
        ρ <= 0 && break
        τ += ρ
    end
    return τ
end

"""
    compute_local_variability(filepath, variable_names) → NamedTuple

Compute several measures of local variability from a converged-equilibrium
NetCDF file (`_on.etc.nc` or `_of.etc.nc`).

The time series is detrended by subtracting its mean before computing
second-order statistics.

# Measures returned
- `mean_state::Vector{Float64}`         — time-mean position (attractor estimate)
- `std_per_dim::Vector{Float64}`        — standard deviation along each EOF axis
- `total_variance::Float64`             — trace of the covariance matrix
- `dominant_variance::Float64`          — largest eigenvalue of the covariance matrix
                                          (variance in the most variable direction)
- `dominant_direction::Vector{Float64}` — corresponding eigenvector
- `mean_dist::Float64`                  — mean Euclidean distance from the attractor mean
- `lag1_autocorr::Vector{Float64}`      — lag-1 autocorrelation per EOF dimension
                                          (classical CSD indicator; → 1 near tipping)
- `ac_times::Vector{Float64}`           — integrated autocorrelation time per EOF dim (yr)
- `mean_lag1_autocorr::Float64`         — mean lag-1 autocorrelation across dims
- `mean_ac_time::Float64`               — mean integrated autocorrelation time (yr)
- `n_samples::Int`                      — number of valid (non-NaN) time steps used
"""
function compute_local_variability(filepath::String, variable_names::Vector{String})
    ts_raw = load_plasim_state_timeseries(filepath, variable_names)

    # Drop rows with any NaN
    valid = [!any(isnan, ts_raw[t, :]) for t in axes(ts_raw, 1)]
    ts    = ts_raw[valid, :]
    n, d  = size(ts)

    μ        = vec(mean(ts; dims = 1))
    centered = ts .- μ'

    # Per-dimension standard deviation
    std_per_dim = vec(std(ts; dims = 1))

    # Covariance matrix and its eigendecomposition
    C      = cov(ts)
    evals  = eigvals(Symmetric(C))        # sorted ascending
    evecs  = eigvecs(Symmetric(C))
    dom_idx = argmax(evals)
    dominant_variance  = evals[dom_idx]
    dominant_direction = evecs[:, dom_idx]

    total_variance = tr(C)

    # Mean distance from attractor
    mean_dist = mean(norm(centered[t, :]) for t in 1:n)

    # Lag-1 autocorrelation and integrated autocorrelation time per dimension
    lag1_autocorr = [cor(centered[1:n-1, k], centered[2:n, k]) for k in 1:d]
    ac_times      = [_integrated_ac_time(centered[:, k])       for k in 1:d]

    return (
        mean_state          = μ,
        covariance          = C,
        std_per_dim         = std_per_dim,
        total_variance      = total_variance,
        dominant_variance   = dominant_variance,
        dominant_direction  = dominant_direction,
        mean_dist           = mean_dist,
        lag1_autocorr       = lag1_autocorr,
        ac_times            = ac_times,
        mean_lag1_autocorr  = mean(lag1_autocorr),
        mean_ac_time        = mean(ac_times),
        n_samples           = n,
    )
end
