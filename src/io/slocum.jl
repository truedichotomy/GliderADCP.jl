# Layer 1 — Slocum glider data ingestion (platform-agnostic tables).
#
# Works from any Slocum-derived table (ERDDAP tabledap export, Python dbdreader, or the
# user's pure-Julia JLDBDReader.jl — https://github.com/truedichotomy/JLDBDReader.jl):
#
#   using JLDBDReader
#   d = MultiDBD(dir="…", eng=true, sci=true)
#   df = DataFrame(get_sync(d, "m_water_vx", "m_water_vy", "m_gps_mag_var",
#                              "m_lat", "m_lon", "m_depth", ...))
#
# then `slocum_nav(df)` / `dac_from_slocum(df)` feed the standard pipeline.

_col(df, names...) = begin
    for n in names
        hasproperty(df, n) && return df[!, n]
    end
    nothing
end

_to_datetime(t::AbstractVector{<:DateTime}) = collect(t)
_to_datetime(t::AbstractVector) = unix2datetime.(Float64.(coalesce.(t, NaN)))

"""
    slocum_nav(df) -> GliderNav

Build a [`GliderNav`](@ref) from a Slocum table. Recognized columns (first match wins):
time (`time`), position (`latitude`/`m_gps_lat`, `longitude`/`m_gps_lon`, decimal deg),
depth (`depth`/`m_depth`), attitude (`m_heading`/`m_pitch`/`m_roll`, radians),
declination (`m_gps_mag_var`, radians). Missing columns become NaN;
`deadreckoning`/`navstate` are set to unknown (Slocum DAC comes from
[`dac_from_slocum`](@ref) instead of DR/GPS jumps).
"""
function slocum_nav(df::DataFrame)
    time = _to_datetime(_col(df, :time))
    n = length(time)
    f(names...; scale=1.0) = begin
        c = _col(df, names...)
        c === nothing ? fill(NaN, n) : Float64.(coalesce.(c, NaN)) .* scale
    end
    GliderNav(time, datetime2unix.(time),
        f(:longitude, :m_gps_lon, :lon), f(:latitude, :m_gps_lat, :lat),
        f(:m_heading, :heading; scale=180 / π),
        f(:m_gps_mag_var; scale=180 / π),
        f(:m_pitch, :pitch; scale=180 / π), f(:m_roll, :roll; scale=180 / π),
        f(:depth, :m_depth), fill(Int16(-1), n), fill(Int8(-1), n), fill(NaN, n), df)
end

"""
    dac_from_slocum(df; by=:source_file, min_depth=10.0) -> DataFrame

Per-segment depth-averaged current from the glider's own dead-reckoned estimate
(`m_water_vx/vy`, magnetic frame), rotated to true east/north by `m_gps_mag_var`
(Gradone et al. 2023 recipe: last non-missing value per segment, mean declination).
Segments come from the `by` column (Slocum `source_file` ≈ surfacing-to-surfacing).
Output matches the [`compute_dac`](@ref) schema used by the solvers
(`yo, t_start, t_end, t_mid, duration, u, v`).
"""
function dac_from_slocum(df::DataFrame; by::Symbol=:source_file, min_depth::Real=10.0)
    hasproperty(df, by) || error("dac_from_slocum: no `$by` column")
    time = _to_datetime(_col(df, :time))
    rows = NamedTuple[]
    yo = 0
    for g in groupby(DataFrame(df; copycols=false), by)
        idx = parentindices(g)[1]
        t = time[idx]
        vx = _col(g, :m_water_vx); vy = _col(g, :m_water_vy)
        (vx === nothing || vy === nothing) && continue
        iv = findlast(i -> !ismissing(vx[i]) && !ismissing(vy[i]), 1:nrow(g))
        iv === nothing && continue
        dep = _col(g, :depth, :m_depth)
        if dep !== nothing
            dmax = maximum(skipmissing(dep); init=0.0)
            dmax < min_depth && continue
        end
        mv = _col(g, :m_gps_mag_var)
        mvdeg = mv === nothing ? 0.0 :
                rad2deg(mean(skipmissing(mv)) isa Number ? mean(skipmissing(mv)) : 0.0)
        u0, v0 = Float64(vx[iv]), Float64(vy[iv])
        c, s = cosd(mvdeg), sind(mvdeg)
        yo += 1
        t1, t2 = extrema(t)
        push!(rows, (yo=yo, t_start=t1, t_end=t2,
            t_mid=t1 + Millisecond(round(Int, (t2 - t1).value / 2)),
            duration=(t2 - t1).value / 1000,
            u=u0 * c - v0 * s, v=u0 * s + v0 * c))
    end
    return sort!(DataFrame(rows), :t_start)
end
