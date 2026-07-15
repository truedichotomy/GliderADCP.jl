# Layer 3 — depth-averaged current (DAC) and surface drift from navigation.
#
# SeaExplorer principle: while submerged (`DeadReckoning == 1`) the vehicle dead-reckons
# its position from heading + a through-water flight model, ignoring currents. At
# surfacing, the first GPS fix (`DeadReckoning == 0`) snaps the position; the jump
# between the last dead-reckoned position and that fix is the current-induced
# displacement accumulated over the submerged interval:
#
#     DAC = (pos_fix − pos_DR_end) / (t_fix − t_submerged_start)
#
# This is the standard glider DAC (Rudnick et al. 2018; Gradone et al. 2023 report
# 1–2 cm/s RMS accuracy for the equivalent Slocum estimate).

const _EARTH_R = 6.371e6  # m

"""
    lonlat_to_dxdy(lon0, lat0, lon1, lat1) -> (dx, dy)

Local-tangent displacement in meters (east, north) from position 0 to position 1
(decimal degrees; spherical earth, cosine of the mean latitude).
"""
function lonlat_to_dxdy(lon0, lat0, lon1, lat1)
    latm = (lat0 + lat1) / 2
    dx = deg2rad(lon1 - lon0) * _EARTH_R * cosd(latm)
    dy = deg2rad(lat1 - lat0) * _EARTH_R
    return dx, dy
end

"""
    compute_dac(nav::GliderNav; min_duration=600.0, min_depth=10.0,
                max_speed=1.5, max_fix_delay=900.0) -> DataFrame

Depth-averaged current per submerged segment from SeaExplorer navigation.

Each maximal `DeadReckoning == 1` block bounded by finite-position GPS fixes yields one
estimate. Quality control drops segments that are too short (`min_duration`, s), too
shallow (`min_depth`, m — surface drift intervals), whose first fix arrives too long
after the last DR record (`max_fix_delay`, s), or that imply unphysical speeds
(`max_speed`, m/s).

Returns a `DataFrame` with one row per accepted segment:
`yo, t_start, t_end, t_mid, duration, lon0, lat0, lon_dr, lat_dr, lon_fix, lat_fix,
maxdepth, u, v` — where `(u, v)` is the DAC (m/s, east/north), `(lon0, lat0)` the last
fix before diving, `(lon_dr, lat_dr)` the final dead-reckoned position and
`(lon_fix, lat_fix)` the first fix after surfacing.
"""
function compute_dac(nav::GliderNav; min_duration::Real=600.0, min_depth::Real=10.0,
                     max_speed::Real=1.5, max_fix_delay::Real=900.0)
    rows = NamedTuple[]
    for (yo, seg) in enumerate(_dr_segments(nav))
        core = _dac_segment(nav, seg; min_duration, min_depth, max_speed, max_fix_delay)
        core === nothing && continue
        push!(rows, (yo = yo, core...))
    end
    return DataFrame(rows)
end

# Candidate submerged segments: maximal DeadReckoning == 1 blocks with a finite
# final dead-reckoned position, bracketed by finite-position GPS fixes. Returns
# (ib, i1, i2, ia): last fix before, block bounds, first fix after. Candidates are
# numbered in order by the compute_dac methods, so yo ids agree across methods.
function _dr_segments(nav::GliderNav)
    n = length(nav)
    dr = nav.deadreckoning
    isfix(i) = dr[i] == 0 && isfinite(nav.lon[i]) && isfinite(nav.lat[i])
    segs = NTuple{4,Int}[]
    i = 1
    while i <= n
        if dr[i] != 1
            i += 1
            continue
        end
        i1 = i
        i2 = i
        while i2 < n && dr[i2+1] == 1
            i2 += 1
        end
        i = i2 + 1
        ib = 0
        for j in i1-1:-1:1
            isfix(j) && (ib = j; break)
        end
        ia = 0
        for j in i2+1:n
            isfix(j) && (ia = j; break)
        end
        (ib == 0 || ia == 0) && continue
        isfinite(nav.lon[i2]) && isfinite(nav.lat[i2]) || continue
        push!(segs, (ib, i1, i2, ia))
    end
    return segs
end

# Onboard-DR DAC core for one candidate segment; `nothing` when quality control
# rejects it. DR error (and current drift) accrues from the last pre-dive fix to
# the first post-surfacing fix — the displacement window is fix-to-fix.
function _dac_segment(nav::GliderNav, (ib, i1, i2, ia)::NTuple{4,Int};
                      min_duration, min_depth, max_speed, max_fix_delay)
    duration = nav.t[ia] - nav.t[ib]
    fix_delay = nav.t[ia] - nav.t[i2]
    depths = filter(isfinite, nav.depth[i1:i2])
    maxdepth = isempty(depths) ? NaN : maximum(depths)

    dx, dy = lonlat_to_dxdy(nav.lon[i2], nav.lat[i2], nav.lon[ia], nav.lat[ia])
    u = dx / duration
    v = dy / duration

    ok = duration >= min_duration && fix_delay <= max_fix_delay &&
         isfinite(maxdepth) && maxdepth >= min_depth &&
         isfinite(u) && isfinite(v) && abs(u) <= max_speed && abs(v) <= max_speed
    ok || return nothing

    return (
        t_start = nav.time[ib], t_end = nav.time[ia],
        t_mid = nav.time[ib] + Millisecond(round(Int, 500duration)),
        duration = duration,
        lon0 = nav.lon[ib], lat0 = nav.lat[ib],
        lon_dr = nav.lon[i2], lat_dr = nav.lat[i2],
        lon_fix = nav.lon[ia], lat_fix = nav.lat[ia],
        maxdepth = maxdepth,
        u = u, v = v,
    )
end

"""
    compute_dac(nav::GliderNav, pings::ProcessedPings;
                fallback=nothing,
                o_min=4.0, o_max=16.0, min_cells=3, coverage_min=0.85,
                dt_grid=10.0, max_gap=30.0, fallback_max_gap=120.0,
                surface_depth=5.0, min_duration=600.0, min_depth=10.0,
                max_speed=1.5, max_fix_delay=900.0) -> DataFrame

**ADCP-referenced (water-track) DAC** — the onboard dead-reckoning is replaced by
the time integral of the ADCP-measured through-water velocity
([`throughwater_velocity`](@ref)) over the same fix-to-fix window:

    DAC = (pos_fix_after − pos_fix_before − ∫ u_tw dt) / T

This removes ALSEAMAR's onboard flight model from the estimate entirely (cf. the
water-track referencing of Todd et al. 2017). On the validated missions the onboard
model runs 5–15 % fast, biasing the standard [`compute_dac`](@ref) 2–4 cm/s against
the direction of travel; the water-track form is limited instead by the mean shear
across the cell offset (≲ 1 cm/s with the default 4–16 m window).

The integral is evaluated on a `dt_grid`-second grid: through-water velocity is
interpolated across gaps up to `max_gap` seconds; grid samples with the glider
shallower than `surface_depth` (from nav) count as zero through-water motion;
remaining unsampled flying time is filled with the segment's mean flying velocity.
Keep `surface_depth` below `min_depth` (the defaults are), so a segment can never
be all-surface yet pass quality control.

Segments whose flying-time coverage falls below `coverage_min` (instrument off —
duty cycling) descend the fallback ladder: with `fallback = flight_model(nav)`
(a [`GliderFlight`](@ref)) the segment is dead-reckoned from the **flight model**
instead (~1.4 cm/s median from the ADCP water track, no systematic bias, on the
validated missions — `fallback_max_gap` spans its apogee masking); segments no
source can cover keep the **onboard estimate**. The `method` column
records the rung per row (`:adcp`, `:flight`, `:onboard`), `coverage` the
flying-time coverage of the series used, and `u_ob, v_ob` always carry the
onboard estimate; QC, the remaining schema and yo numbering match the
one-argument method.

Pass the same post-QC `pings` the solvers consume. Works on any route with enough
near cells — including realtime-telemetered pings (`load_pld_adcp`), where cells
within the 4–16 m window exist for typical 2 m cell configurations.
"""
function compute_dac(nav::GliderNav, pings::ProcessedPings;
                     fallback::Union{Nothing,GliderFlight}=nothing,
                     o_min::Real=4.0, o_max::Real=16.0, min_cells::Int=3,
                     coverage_min::Real=0.85, dt_grid::Real=10.0,
                     max_gap::Real=30.0, fallback_max_gap::Real=120.0,
                     surface_depth::Real=5.0,
                     min_duration::Real=600.0, min_depth::Real=10.0,
                     max_speed::Real=1.5, max_fix_delay::Real=900.0)
    tw = throughwater_velocity(pings; o_min, o_max, min_cells)
    sources = [(t=pings.t, u=tw.u, v=tw.v, label=:adcp, max_gap=Float64(max_gap))]
    fallback === nothing ||
        push!(sources, (; _flight_tw(fallback, nav)..., label=:flight,
                          max_gap=Float64(fallback_max_gap)))
    return _dac_watertrack(nav, sources; coverage_min, dt_grid, surface_depth,
                           min_duration, min_depth, max_speed, max_fix_delay)
end

"""
    compute_dac(nav::GliderNav, fl::GliderFlight;
                max_gap=120.0, coverage_min=0.85, dt_grid=10.0,
                surface_depth=5.0, min_duration=600.0, min_depth=10.0,
                max_speed=1.5, max_fix_delay=900.0) -> DataFrame

**Flight-model water-track DAC** — for deployments without a (running) ADCP:
the onboard dead-reckoning is replaced by the integral of the flight-model
through-water velocity (`U cos γ` from [`flight_model`](@ref), projected along
the true nav heading). On the validated missions this lands within ~1.4 cm/s
median of the ADCP water track with **no systematic along-track bias**, versus
4–6.5 cm/s systematically anti-track for the onboard model — so it is the
preferred DAC wherever ADCP pings are unavailable. Assumes travel along heading
(no sideslip) and zero vertical water velocity (strong convection inflates the
per-yo scatter), and its accuracy tracks the polar's provenance: a same-glider
calibrated polar gives the ~1.4 cm/s figure; published airframe presets sit
between that and the onboard error. Segments the flight model cannot cover
(`coverage_min` of flying time; `max_gap` spans apogee masking) keep the
onboard estimate. Columns as in the pings method, with `method` ∈
`(:flight, :onboard)`.

```julia
fl  = flight_model(nav)                        # or with 1 Hz CTD pressure
dac = compute_dac(nav, fl)
```
"""
function compute_dac(nav::GliderNav, fl::GliderFlight;
                     max_gap::Real=120.0, coverage_min::Real=0.85,
                     dt_grid::Real=10.0, surface_depth::Real=5.0,
                     min_duration::Real=600.0, min_depth::Real=10.0,
                     max_speed::Real=1.5, max_fix_delay::Real=900.0)
    src = (; _flight_tw(fl, nav)..., label=:flight, max_gap=Float64(max_gap))
    return _dac_watertrack(nav, [src]; coverage_min, dt_grid, surface_depth,
                           min_duration, min_depth, max_speed, max_fix_delay)
end

# Water-track DAC over an ordered list of through-water sources — NamedTuples
# (t [unix s], u, v, label::Symbol, max_gap) — tried per segment in order; the
# first whose flying-time coverage reaches coverage_min supplies the estimate,
# else the onboard estimate stands (method = :onboard, coverage = best seen).
function _dac_watertrack(nav::GliderNav, sources;
                         coverage_min::Real, dt_grid::Real, surface_depth::Real,
                         min_duration::Real, min_depth::Real, max_speed::Real,
                         max_fix_delay::Real)
    rows = NamedTuple[]
    counts = Dict{Symbol,Int}()
    for (yo, seg) in enumerate(_dr_segments(nav))
        core = _dac_segment(nav, seg; min_duration, min_depth,
                            max_speed=Inf, max_fix_delay)   # speed-check the final estimate instead
        core === nothing && continue
        ib, _, _, ia = seg
        u, v = core.u, core.v
        method = :onboard
        cov = 0.0
        for src in sources
            sx, sy, c = _tw_displacement(src.t, src.u, src.v, nav,
                                         nav.t[ib], nav.t[ia];
                                         dt_grid, max_gap=src.max_gap, surface_depth)
            if c >= coverage_min
                dxf, dyf = lonlat_to_dxdy(nav.lon[ib], nav.lat[ib], nav.lon[ia], nav.lat[ia])
                u = (dxf - sx) / core.duration
                v = (dyf - sy) / core.duration
                method = src.label
                cov = c
                break
            end
            cov = max(cov, c)
        end
        (isfinite(u) && isfinite(v) && abs(u) <= max_speed && abs(v) <= max_speed) || continue
        counts[method] = get(counts, method, 0) + 1
        push!(rows, merge((yo = yo,), core,
                          (u = u, v = v, u_ob = core.u, v_ob = core.v,
                           coverage = cov, method = method)))
    end
    df = DataFrame(rows)
    labelname = Dict(:adcp => "ADCP water-track", :flight => "flight-model")
    parts = ["$(get(counts, s.label, 0)) $(get(labelname, s.label, String(s.label)))"
             for s in sources]
    push!(parts, "$(get(counts, :onboard, 0)) onboard fallback")
    @info "compute_dac: $(nrow(df)) segments — " * join(parts, ", ")
    return df
end

# Through-water displacement (m east/north) over [t1, t2] and its flying-time
# coverage. Sampled on a dt_grid-second grid: interpolated u_tw where available
# (gaps ≤ max_gap bridged), zero where the glider is at the surface, and the
# residual unsampled flying time mean-filled (displacement = mean flying velocity
# × flying time). coverage = sampled fraction of the flying time.
function _tw_displacement(tp::Vector{Float64}, u::Vector{Float64}, v::Vector{Float64},
                          nav::GliderNav, t1::Real, t2::Real;
                          dt_grid::Real, max_gap::Real, surface_depth::Real)
    tg = collect(t1:dt_grid:t2)
    ug = _interp_capped(tg, tp, u; max_gap)
    vg = _interp_capped(tg, tp, v; max_gap)
    dg = _interp_capped(tg, nav.t, nav.depth; max_gap=300.0)
    n = length(tg)
    nsurf = nval = 0
    su = sv = 0.0
    for k in 1:n
        if isfinite(ug[k]) && isfinite(vg[k])
            su += ug[k]; sv += vg[k]; nval += 1
        elseif isfinite(dg[k]) && dg[k] < surface_depth
            nsurf += 1
        end
    end
    nfly = n - nsurf
    nfly <= 0 && return (0.0, 0.0, 1.0)          # never flying: nothing to integrate
    cov = nval / nfly
    T_fly = (t2 - t1) * nfly / n
    nval == 0 && return (0.0, 0.0, 0.0)
    return (su / nval * T_fly, sv / nval * T_fly, cov)
end

# Linear interpolation that never bridges source gaps wider than `max_gap` seconds
# (unlike `_interp1`, which interpolates across any interior gap).
function _interp_capped(tq::AbstractVector{<:Real}, ts::Vector{Float64},
                        vs::Vector{Float64}; max_gap::Real)
    keep = findall(isfinite, vs)
    out = fill(NaN, length(tq))
    isempty(keep) && return out
    tk = ts[keep]; vk = vs[keep]
    for (i, t) in enumerate(tq)
        (t < tk[1] || t > tk[end]) && continue
        j = clamp(searchsortedlast(tk, t), 1, length(tk) - 1)
        tk[j+1] - tk[j] <= max_gap || continue
        f = tk[j+1] == tk[j] ? 0.0 : (t - tk[j]) / (tk[j+1] - tk[j])
        out[i] = (1 - f) * vk[j] + f * vk[j+1]
    end
    return out
end

"""
    surface_drift(nav::GliderNav; min_gap=30.0, max_gap=1200.0, max_speed=1.5) -> DataFrame

Near-surface drift velocities from consecutive GPS fixes within the same surface
interval (no submerged records between them). Returns `t_mid, duration, lon, lat, u, v`
per fix pair — a near-surface velocity constraint for the inverse solution.
"""
function surface_drift(nav::GliderNav; min_gap::Real=30.0, max_gap::Real=1200.0,
                       max_speed::Real=1.5)
    n = length(nav)
    rows = NamedTuple[]
    prev = 0
    for i in 1:n
        if nav.deadreckoning[i] == 1
            prev = 0
        elseif nav.deadreckoning[i] == 0 && isfinite(nav.lon[i]) && isfinite(nav.lat[i])
            if prev > 0
                dt = nav.t[i] - nav.t[prev]
                if min_gap <= dt <= max_gap
                    dx, dy = lonlat_to_dxdy(nav.lon[prev], nav.lat[prev], nav.lon[i], nav.lat[i])
                    u, v = dx / dt, dy / dt
                    if abs(u) <= max_speed && abs(v) <= max_speed
                        push!(rows, (
                            t_mid = nav.time[prev] + Millisecond(round(Int, 500dt)),
                            duration = dt,
                            lon = (nav.lon[prev] + nav.lon[i]) / 2,
                            lat = (nav.lat[prev] + nav.lat[i]) / 2,
                            u = u, v = v,
                        ))
                    end
                end
            end
            prev = i
        end
    end
    return DataFrame(rows)
end
