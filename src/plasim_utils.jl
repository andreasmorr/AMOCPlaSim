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
  4. `compute_convergence_times`       — time to converge to attractor
  5. `compute_edge_to_attractor_distances` — distance from edge state to attractor
  6. `plasim_resilience_summary` — convenience wrapper for all of the above
"""

using NCDatasets
using DataFrames
using Statistics
using LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

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
                        data_cols[Symbol("x$k")] = Float64.(raw[track, :])
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
                        data_cols[Symbol("x$k")] = Float64.(raw[:, track])
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

The AMOC-on state has a higher value of redu1 (stronger overturning), and the
AMOC-off state has a lower value. We split by the median of final redu1 values.

Returns a vector of length `n_trajectories` with values 1 (AMOC-on) or 2
(AMOC-off), in the same order as `sort(unique(df.trajectory_id))`.
"""
function classify_trajectories(df::DataFrame, n_dims::Int;
                                final_fraction::Float64 = 0.1)
    ids, final_states = _get_final_states(df, n_dims; final_fraction)

    # Use first EOF (redu1) to separate AMOC-on from AMOC-off
    # Convention: higher redu1 = AMOC-on
    redu1_vals = final_states[:, 1]
    threshold  = median(redu1_vals)

    labels = [v >= threshold ? 1 : 2 for v in redu1_vals]
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
    compute_convergence_times(df, n_dims, attractors; ε=0.05) → Dict{Int, Float64}

For each trajectory, find the first time step at which the trajectory enters
and stays within ε (Euclidean distance) of its target attractor.

"Stays within" means all subsequent time steps also remain within ε (to avoid
transient close passes).

If a trajectory never converges within ε, the convergence time is set to the
final time step of that trajectory.

# Arguments
- `df`:         DataFrame from `load_plasim_trajectories`
- `n_dims`:     number of EOF dimensions
- `attractors`: Dict from `estimate_attractors`
- `ε`:          convergence threshold (in EOF units)

# Returns
Dict mapping `trajectory_id => convergence_time` (as a Float64 time-step index).
"""
function compute_convergence_times(df::DataFrame, n_dims::Int,
                                    attractors::Dict{Int, Vector{Float64}};
                                    ε::Float64 = 0.05)
    labels = classify_trajectories(df, n_dims)
    ids    = sort(unique(df.trajectory_id))
    result = Dict{Int, Float64}()

    for (j, tid) in enumerate(ids)
        target_attractor = attractors[labels[j]]
        traj = sort(filter(r -> r.trajectory_id == tid, df), :time)

        times = traj.time
        n_pts = nrow(traj)

        # Build matrix of states: (n_pts × n_dims)
        states = hcat([traj[:, Symbol("x$k")] for k in 1:n_dims]...)

        # Compute distance to target attractor at each time step
        dists = [norm(states[t, :] .- target_attractor) for t in 1:n_pts]

        # Find first time step where distance stays ≤ ε for all remaining steps
        conv_time = times[end]  # default: last time step
        for t in 1:n_pts
            if all(dists[t:end] .<= ε)
                conv_time = times[t]
                break
            end
        end

        result[tid] = conv_time
    end

    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Edge-to-attractor distances
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_edge_to_attractor_distances(df, n_dims, attractors) → Dict{Int, Float64}

For each trajectory, compute the Euclidean distance from the initial state
(the edge/saddle state at t=1) to the trajectory's target attractor.

This measures how far the edge state is from each attractor in EOF space.

# Arguments
- `df`:         DataFrame from `load_plasim_trajectories`
- `n_dims`:     number of EOF dimensions
- `attractors`: Dict from `estimate_attractors`

# Returns
Dict mapping `trajectory_id => distance_from_edge_to_attractor`.
"""
function compute_edge_to_attractor_distances(df::DataFrame, n_dims::Int,
                                              attractors::Dict{Int, Vector{Float64}})
    labels = classify_trajectories(df, n_dims)
    ids    = sort(unique(df.trajectory_id))
    result = Dict{Int, Float64}()

    for (j, tid) in enumerate(ids)
        target_attractor = attractors[labels[j]]
        traj = sort(filter(r -> r.trajectory_id == tid, df), :time)

        # Edge state = initial state (first time step)
        edge_state = [traj[1, Symbol("x$k")] for k in 1:n_dims]

        result[tid] = norm(edge_state .- target_attractor)
    end

    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

"""
    plasim_resilience_summary(df, n_dims; final_fraction=0.1, ε=0.05) → NamedTuple

Convenience wrapper that computes all PlaSim resilience diagnostics.

Calls `classify_trajectories`, `estimate_attractors`, `compute_convergence_times`,
and `compute_edge_to_attractor_distances`, then aggregates the results.

# Returns
A NamedTuple with fields:
- `attractor_labels::Vector{Int}`      — per-trajectory label (1=on, 2=off)
- `trajectory_ids::Vector{Int}`        — corresponding trajectory IDs
- `attractors::Dict{Int,Vector{Float64}}` — estimated attractor positions
- `convergence_times::Dict{Int,Float64}`  — per-trajectory convergence times
- `edge_distances::Dict{Int,Float64}`     — edge-to-attractor distances
- `n_on::Int`                           — number of AMOC-on trajectories
- `n_off::Int`                          — number of AMOC-off trajectories
- `mean_conv_time_on::Float64`          — mean convergence time for AMOC-on
- `mean_conv_time_off::Float64`         — mean convergence time for AMOC-off
- `mean_dist_on::Float64`               — mean edge distance for AMOC-on
- `mean_dist_off::Float64`              — mean edge distance for AMOC-off
"""
function plasim_resilience_summary(df::DataFrame, n_dims::Int;
                                    final_fraction::Float64 = 0.1,
                                    ε::Float64 = 0.05)
    ids    = sort(unique(df.trajectory_id))
    labels = classify_trajectories(df, n_dims; final_fraction)
    attrs  = estimate_attractors(df, n_dims; final_fraction)
    conv_t = compute_convergence_times(df, n_dims, attrs; ε)
    edist  = compute_edge_to_attractor_distances(df, n_dims, attrs)

    on_idx  = findall(==(1), labels)
    off_idx = findall(==(2), labels)

    on_ids  = ids[on_idx]
    off_ids = ids[off_idx]

    n_on  = length(on_idx)
    n_off = length(off_idx)

    mean_conv_on  = n_on  > 0 ? mean(conv_t[tid] for tid in on_ids)  : NaN
    mean_conv_off = n_off > 0 ? mean(conv_t[tid] for tid in off_ids) : NaN

    mean_dist_on  = n_on  > 0 ? mean(edist[tid] for tid in on_ids)  : NaN
    mean_dist_off = n_off > 0 ? mean(edist[tid] for tid in off_ids) : NaN

    return (
        attractor_labels   = labels,
        trajectory_ids     = ids,
        attractors         = attrs,
        convergence_times  = conv_t,
        edge_distances     = edist,
        n_on               = n_on,
        n_off              = n_off,
        mean_conv_time_on  = mean_conv_on,
        mean_conv_time_off = mean_conv_off,
        mean_dist_on       = mean_dist_on,
        mean_dist_off      = mean_dist_off,
    )
end
