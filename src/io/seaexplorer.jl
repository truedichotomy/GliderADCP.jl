# Layer 1 — SeaExplorer navigation (.gli.sub/.raw) and payload (.pld1/.legato/...) parsers.
#
# ALSEAMAR file conventions (verified on sea064 M38):
#   - semicolon-separated with a trailing `;` (empty last column), optionally gzipped
#   - segment-numbered names: `sea064.38.gli.sub.<N>[.gz]`, `sea064.38.pld1.raw.<N>[.gz]`,
#     `sea064.38.legato.raw.<N>[.gz]`, `sea064.38.ad2cp.raw.<N>[.gz]` — N needs natural sort
#   - nav `Timestamp` = "dd/mm/yyyy HH:MM:SS"; payload `PLD_REALTIMECLOCK` adds ".sss"
#   - Lat/Lon in NMEA degrees·100 + decimal minutes (7001.296 = 70° 01.296′)
#   - nav `DeadReckoning`: 1 subsurface (dead-reckoned), 0 surface GPS fix; the position
#     jump at the 1→0 transition carries the depth-averaged current displacement.

const _SEAEXPLORER_NAV_FMT = dateformat"dd/mm/yyyy HH:MM:SS"
const _SEAEXPLORER_PLD_FMT = dateformat"dd/mm/yyyy HH:MM:SS.sss"

"""
    nmea2deg(x) -> Float64

Convert an NMEA-style coordinate (degrees·100 + decimal minutes, sign = hemisphere) to
decimal degrees. `NaN` passes through. Example: `nmea2deg(7001.296) ≈ 70.0216`.
"""
function nmea2deg(x::Real)
    isfinite(x) || return Float64(x)
    a = abs(x)
    d = floor(a / 100)
    copysign(d + (a - 100d) / 60, x)
end

"""
    seaexplorer_files(dir, stream) -> Vector{String}

List segment files of one SeaExplorer stream (e.g. `"gli.sub"`, `"pld1.raw"`,
`"legato.raw"`, `"ad2cp.raw"`) in `dir`, naturally sorted by segment number.
"""
function seaexplorer_files(dir::AbstractString, stream::AbstractString)
    pat = Regex("\\." * replace(stream, "." => "\\.") * "\\.(\\d+)(\\.gz)?\$")
    hits = Tuple{Int,String}[]
    for f in readdir(dir)
        m = match(pat, f)
        m === nothing || push!(hits, (parse(Int, m.captures[1]), joinpath(dir, f)))
    end
    return last.(sort(hits))
end

# read one (possibly gzipped) semicolon-separated segment file into a DataFrame
function _read_segment(path::AbstractString; timestamp_col::String)
    bytes = endswith(path, ".gz") ?
        open(io -> read(GzipDecompressorStream(io)), path) : read(path)
    df = CSV.read(bytes, DataFrame;
        delim=';', missingstring=["", "NaN"], types=Dict(timestamp_col => String),
        silencewarnings=true, strict=false)
    # drop the phantom column created by the trailing semicolon
    if !isempty(names(df)) && all(ismissing, df[!, end])
        select!(df, Not(names(df)[end]))
    end
    return df
end

function _read_stream(files::AbstractVector{<:AbstractString}; timestamp_col::String,
                      fmt::DateFormat, mintime::Union{DateTime,Nothing})
    isempty(files) && error("no SeaExplorer files to read")
    dfs = [_read_segment(f; timestamp_col) for f in files]
    df = reduce((a, b) -> vcat(a, b; cols=:union), dfs)
    ts = [t === missing ? missing : tryparse(DateTime, t, fmt) for t in df[!, timestamp_col]]
    df.time = ts
    # drop unparseable timestamps and (by default) boot-time records logged before the
    # glider's clock is set (epoch-1970 stamps, typically with Lat=Lon=0)
    good = findall(t -> t !== missing && (mintime === nothing || t >= mintime), ts)
    ndropped = nrow(df) - length(good)
    ndropped > 0 && @debug "SeaExplorer stream: dropped $ndropped invalid-timestamp records"
    df = df[good, :]
    sort!(df, :time)
    return df
end

_float_col(df, name) = hasproperty(df, name) ?
    Float64.(coalesce.(df[!, name], NaN)) : fill(NaN, nrow(df))
_int_col(df, name, ::Type{T}, default) where {T} = hasproperty(df, name) ?
    T.(coalesce.(df[!, name], default)) : fill(T(default), nrow(df))

"""
    load_seaexplorer_nav(src; stream="gli.sub", mintime=DateTime(2000)) -> GliderNav

Read SeaExplorer navigation files (`.gli`) into a [`GliderNav`](@ref). `src` is a
directory (all segments of `stream`, naturally sorted) or an explicit vector of file
paths. Positions are converted from NMEA to decimal degrees. Records timestamped before
`mintime` (boot records logged before the clock is set, stamped 1970) are dropped;
pass `mintime=nothing` to keep everything.

```julia
nav = load_seaexplorer_nav("…/delayed/nav/logs")
```
"""
function load_seaexplorer_nav(src; stream::AbstractString="gli.sub",
                              mintime::Union{DateTime,Nothing}=DateTime(2000))
    files = src isa AbstractVector ? String.(src) : seaexplorer_files(src, stream)
    df = _read_stream(files; timestamp_col="Timestamp", fmt=_SEAEXPLORER_NAV_FMT, mintime)
    time = Vector{DateTime}(df.time)
    GliderNav(time, datetime2unix.(time),
        nmea2deg.(_float_col(df, "Lon")),
        nmea2deg.(_float_col(df, "Lat")),
        _float_col(df, "Heading"),
        _float_col(df, "Declination"),
        _float_col(df, "Pitch"),
        _float_col(df, "Roll"),
        _float_col(df, "Depth"),
        _int_col(df, "NavState", Int16, -1),
        _int_col(df, "DeadReckoning", Int8, -1),
        _float_col(df, "Altitude"),
        df)
end

"""
    load_seaexplorer_pld(src; stream="pld1.raw", mintime=DateTime(2000)) -> DataFrame

Read SeaExplorer payload science files (`pld1`, `legato`, …) into a time-sorted
`DataFrame` (column set is payload-configuration dependent, so no strong typing here).
A `time::DateTime` column is added from `PLD_REALTIMECLOCK`. `src` is a directory or a
vector of file paths — pass a subset of segments to limit memory.

```julia
pld = load_seaexplorer_pld("…/delayed/pld1/logs"; stream="legato.raw")
```
"""
function load_seaexplorer_pld(src; stream::AbstractString="pld1.raw",
                              mintime::Union{DateTime,Nothing}=DateTime(2000))
    files = src isa AbstractVector ? String.(src) : seaexplorer_files(src, stream)
    _read_stream(files; timestamp_col="PLD_REALTIMECLOCK", fmt=_SEAEXPLORER_PLD_FMT, mintime)
end
