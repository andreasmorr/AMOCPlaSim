"""
    extract_amoc_strength.jl

Extract AMOC-strength time series from the original PlaSim reduced-coordinate
NetCDF files into a compact standalone sidecar file.

Output:
    data/custom_readouts/amoc_strength_timeseries.nc

Each output variable is named after the original source filename without the
`.etc.nc` suffix. For example:

    plasimelancholia_285ppm_on
    plasimelancholia_285ppm_edgetrack_itx000

Equilibrium variables are one-dimensional over their own time coordinate.
Edge-track variables are two-dimensional `(time, track)` arrays.

Run from the AMOCPlaSim directory:
    julia --project scripts/extract_amoc_strength.jl
"""

using DrWatson
@quickactivate "AMOCResilience"

using Dates
using NCDatasets

const SOURCE_DIR = datadir("results")
const OUT_PATH = datadir("custom_readouts", "amoc_strength_timeseries.nc")

source_files() = sort(filter(
    path -> endswith(basename(path), ".etc.nc"),
    readdir(SOURCE_DIR; join = true),
))

variable_key(path::AbstractString) = replace(basename(path), r"\.etc\.nc$" => "")
time_dim_name(key::AbstractString) = "time__$(key)"

function read_time_coordinate(src::NCDataset, n_time::Int)
    if haskey(src, "time")
        time = Float64.(coalesce.(src["time"][:], NaN))
        length(time) == n_time && return time
    end
    return Float64.(0:(n_time - 1))
end

function variable_attributes(src_path::AbstractString, source_var)
    attrs = Dict{String, Any}(
        "source_file" => basename(src_path),
        "long_name" => "AMOC strength",
        "units" => "Sv",
    )
    for attr in ("long_name", "units")
        if attr in keys(source_var.attrib)
            attrs[attr] = source_var.attrib[attr]
        end
    end
    return attrs
end

function copy_amoc_strength!(out::NCDataset, src_path::AbstractString)
    key = variable_key(src_path)
    dim_time = time_dim_name(key)

    NCDataset(src_path, "r") do src
        haskey(src, "amoc_strength") ||
            error("Missing amoc_strength variable in $(src_path)")

        source_var = src["amoc_strength"]
        n_dims = length(dimnames(source_var))
        attrs = variable_attributes(src_path, source_var)

        if n_dims == 1
            data = Float64.(coalesce.(source_var[:], NaN))
            defDim(out, dim_time, length(data))

            time_var = defVar(out, dim_time, Float64, (dim_time,),
                attrib = Dict(
                    "long_name" => "Time coordinate copied from $(basename(src_path))",
                    "units" => haskey(src, "time") && "units" in keys(src["time"].attrib) ?
                               src["time"].attrib["units"] : "years",
                    "source_file" => basename(src_path),
                ),
            )
            time_var[:] = read_time_coordinate(src, length(data))

            out_var = defVar(out, key, Float64, (dim_time,), attrib = attrs)
            out_var[:] = data
        elseif n_dims == 2
            data = Float64.(coalesce.(source_var[:, :], NaN))
            size(data, 2) == 2 ||
                error("Expected two tracks in $(src_path), found $(size(data, 2))")
            defDim(out, dim_time, size(data, 1))

            time_var = defVar(out, dim_time, Float64, (dim_time,),
                attrib = Dict(
                    "long_name" => "Time coordinate copied from $(basename(src_path))",
                    "units" => haskey(src, "time") && "units" in keys(src["time"].attrib) ?
                               src["time"].attrib["units"] : "years",
                    "source_file" => basename(src_path),
                ),
            )
            time_var[:] = read_time_coordinate(src, size(data, 1))

            out_var = defVar(out, key, Float64, (dim_time, "track"), attrib = attrs)
            out_var[:, :] = data
        else
            error("Unsupported amoc_strength dimensionality in $(src_path): $(n_dims)")
        end
    end
end

function main()
    files = source_files()
    isempty(files) && error("No .etc.nc files found in $(SOURCE_DIR)")

    isfile(OUT_PATH) && rm(OUT_PATH)

    mkpath(dirname(OUT_PATH))
    NCDataset(OUT_PATH, "c") do out
        out.attrib["title"] = "PlaSim-LSG AMOC-strength time series"
        out.attrib["created_by"] = "scripts/extract_amoc_strength.jl"
        out.attrib["created_on"] = string(Dates.now())
        out.attrib["source_directory"] = SOURCE_DIR
        out.attrib["variable_key"] = "Variables are keyed by source filename without .etc.nc."

        defDim(out, "track", 2)
        track_var = defVar(out, "track", Int32, ("track",),
            attrib = Dict(
                "long_name" => "Edge-track branch index",
                "track_1" => "upper branch, converging to AMOC-on",
                "track_2" => "lower branch, converging to AMOC-off",
            ),
        )
        track_var[:] = Int32[1, 2]

        for src_path in files
            copy_amoc_strength!(out, src_path)
            @info "Extracted AMOC strength" source = basename(src_path)
        end
    end

    @info "AMOC-strength sidecar written" path = OUT_PATH n_sources = length(files)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
