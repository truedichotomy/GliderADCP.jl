# Cross-cutting diagnostics — data coverage and gap reporting.
#
# Loaders never crash on missing/corrupt inputs (they skip with warnings); these
# functions turn "what actually loaded" into a structured, inspectable report so gaps
# are surfaced rather than silently swallowed.

"""
    data_gaps(t::AbstractVector; max_gap=:auto) -> DataFrame

Gaps in a sorted unix-seconds time base: rows `(start, stop, duration)` (DateTimes and
seconds) wherever consecutive samples are more than `max_gap` seconds apart
(`:auto` = 10× the median spacing, at least 60 s). Methods exist for
[`AD2CPData`](@ref), [`GliderNav`](@ref) and [`ProcessedPings`](@ref).
"""
function data_gaps(t::AbstractVector{<:Real}; max_gap::Union{Real,Symbol}=:auto)
    out = DataFrame(start=DateTime[], stop=DateTime[], duration=Float64[])
    n = length(t)
    n < 2 && return out
    dts = diff(t)
    thr = max_gap === :auto ? max(10 * median(dts), 60.0) : Float64(max_gap)
    for i in 1:n-1
        dts[i] > thr &&
            push!(out, (unix2datetime(t[i]), unix2datetime(t[i+1]), dts[i]))
    end
    return out
end
data_gaps(a::AD2CPData; kwargs...) = data_gaps(a.t; kwargs...)
data_gaps(n::GliderNav; kwargs...) = data_gaps(n.t; kwargs...)
data_gaps(p::ProcessedPings; kwargs...) = data_gaps(p.t; kwargs...)

_finitefrac(A) = isempty(A) ? NaN : count(isfinite, A) / length(A)

"""
    coverage(x; max_gap=:auto) -> NamedTuple

Structured coverage report for a loaded dataset — time span, record count, median
sampling interval, the [`data_gaps`](@ref) table with total gap duration, and
finite-data fractions per key field. Methods for [`AD2CPData`](@ref) (adds beam-sample
and sensor fractions plus bottom-track record count), [`GliderNav`](@ref) (adds GPS-fix
count and position coverage) and [`ProcessedPings`](@ref) (adds ENU sample fraction and
the count of pings with no usable cells).
"""
function coverage(a::AD2CPData; max_gap::Union{Real,Symbol}=:auto)
    g = data_gaps(a; max_gap)
    (n=length(a),
        t_start=isempty(a.time) ? missing : first(a.time),
        t_end=isempty(a.time) ? missing : last(a.time),
        median_dt=length(a.t) > 1 ? median(diff(a.t)) : NaN,
        gaps=g, gap_total=sum(g.duration; init=0.0),
        finite_vel=_finitefrac(a.vel), finite_amp=_finitefrac(a.amp),
        finite_corr=_finitefrac(a.corr),
        finite_heading=_finitefrac(a.heading), finite_pitch=_finitefrac(a.pitch),
        finite_pressure=_finitefrac(a.pressure),
        n_bt=a.bt === nothing ? 0 : length(a.bt))
end

function coverage(nav::GliderNav; max_gap::Union{Real,Symbol}=:auto)
    g = data_gaps(nav; max_gap)
    (n=length(nav),
        t_start=isempty(nav.time) ? missing : first(nav.time),
        t_end=isempty(nav.time) ? missing : last(nav.time),
        median_dt=length(nav.t) > 1 ? median(diff(nav.t)) : NaN,
        gaps=g, gap_total=sum(g.duration; init=0.0),
        n_gps_fixes=count(==(0), nav.deadreckoning),
        dr_unknown_frac=count(==(-1), nav.deadreckoning) / max(1, length(nav)),
        finite_lon=_finitefrac(nav.lon), finite_lat=_finitefrac(nav.lat),
        finite_depth=_finitefrac(nav.depth))
end

function coverage(p::ProcessedPings; max_gap::Union{Real,Symbol}=:auto)
    g = data_gaps(p; max_gap)
    empty_pings = count(i -> !any(isfinite, @view p.E[:, i]), 1:length(p))
    (n=length(p),
        t_start=isempty(p.time) ? missing : first(p.time),
        t_end=isempty(p.time) ? missing : last(p.time),
        median_dt=length(p.t) > 1 ? median(diff(p.t)) : NaN,
        gaps=g, gap_total=sum(g.duration; init=0.0),
        finite_E=_finitefrac(p.E),
        empty_pings=empty_pings)
end
