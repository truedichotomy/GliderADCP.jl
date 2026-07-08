# Layer 1 — SeaExplorer navigation (.gli) and payload (.pld1/.legato/…) loading.
#
# The file parsing itself lives in SeaExplorerIO.jl (shared with ATOMIXjulia.jl,
# so loader fixes and new sensors land once): stream discovery, missing-segment
# detection, corrupt-file skipping, NMEA and timestamp normalization. This file
# only adapts its GliderTable into the types the ADCP pipeline uses (GliderNav,
# DataFrame) and re-exports the bookkeeping helpers.
#
# ALSEAMAR file-format notes live in SeaExplorerIO; the DeadReckoning semantics
# the DAC computation relies on: 1 subsurface (dead-reckoned), 0 surface GPS
# fix; the position jump at the 1→0 transition carries the depth-averaged
# current displacement.

"""
    nmea2deg(x) -> Float64

Convert an NMEA-style coordinate (degrees·100 + decimal minutes, sign = hemisphere) to
decimal degrees. `NaN` passes through. Example: `nmea2deg(7001.296) ≈ 70.0216`.
(Alias of `SeaExplorerIO.nmea_to_deg`.)
"""
const nmea2deg = SeaExplorerIO.nmea_to_deg

# DataFrame view of a GliderTable: NaN → missing to keep the historical
# `skipmissing`/`ismissing` idioms working on payload science columns.
function _gt_dataframe(t::SeaExplorerIO.GliderTable)
    df = DataFrame()
    df.time = copy(t.time)
    for (k, v) in t.cols
        df[!, k] = [isnan(x) ? missing : x for x in v]
    end
    return df
end

"""
    load_seaexplorer_nav(src; stream="gli.sub", mintime=DateTime(2000)) -> GliderNav

Read SeaExplorer navigation files (`.gli`) into a [`GliderNav`](@ref). `src` is a
directory (all segments of `stream`, naturally sorted, with a warning listing any
missing segment numbers), an explicit vector of file paths, or a vector of directories
covering several download routes — e.g. glider-computer files plus a GLIMPSE-server
export dir — which are merged and deduplicated by timestamp (earlier directories win;
see `SeaExplorerIO.merge_tables`). Positions are converted from NMEA to decimal
degrees. Records timestamped before `mintime` (boot records logged before the clock is
set, stamped 1970) are dropped; pass `mintime=nothing` to keep everything.

```julia
nav = load_seaexplorer_nav("…/delayed/nav/logs")
nav = load_seaexplorer_nav(["…/delayed/nav/logs", "…/glimpse"])   # both routes, deduped
```
"""
function load_seaexplorer_nav(src; stream::AbstractString="gli.sub",
                              mintime::Union{DateTime,Nothing}=DateTime(2000))
    t = SeaExplorerIO.read_gli(src; stream, epoch_min=mintime)
    n = length(t)
    col(name) = haskey(t, name) ? t[name] : fill(NaN, n)
    icol(name, ::Type{T}, default) where {T} =
        T[isfinite(v) ? T(v) : T(default) for v in col(name)]
    GliderNav(t.time, datetime2unix.(t.time),
        col("Lon"), col("Lat"),
        col("Heading"), col("Declination"), col("Pitch"), col("Roll"), col("Depth"),
        icol("NavState", Int16, -1), icol("DeadReckoning", Int8, -1),
        col("Altitude"), _gt_dataframe(t))
end

"""
    load_seaexplorer_pld(src; stream="pld1.raw", mintime=DateTime(2000)) -> DataFrame

Read SeaExplorer payload science files (`pld1`, `legato`, …) into a time-sorted
`DataFrame` with a `time::DateTime` column; science columns are `Float64` with
`missing` for empty cells (±9999 instrument-off sentinels included). Known coordinate
columns (`NAV_LATITUDE`/`NAV_LONGITUDE`) are converted from NMEA to decimal degrees.
`src` is a directory, a vector of file paths, or a vector of directories covering
several download routes (glider computer + GLIMPSE server) — merged and deduplicated
by timestamp with the highest resolution preserved. `stream` may be a priority-ranked
vector, e.g. `["pld1.raw", "pld1.sub"]` to let telemetered rows fill raw-file gaps.
Pass a subset of segments (or use `SeaExplorerIO.read_pld` with a column selection) to
limit memory on full-resolution streams.

```julia
pld = load_seaexplorer_pld("…/delayed/pld1/logs"; stream="legato.raw")
pld = load_seaexplorer_pld(["…/delayed/pld1/logs", "…/glimpse"];
                           stream=["pld1.raw", "pld1.sub"])   # all routes, deduped
```
"""
function load_seaexplorer_pld(src;
                              stream::Union{AbstractString,AbstractVector{<:AbstractString}}="pld1.raw",
                              mintime::Union{DateTime,Nothing}=DateTime(2000))
    _gt_dataframe(SeaExplorerIO.read_stream(src, stream;
        skip_empty=false, epoch_min=mintime))
end
