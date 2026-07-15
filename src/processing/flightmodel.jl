# Layer 2 — steady-state glider flight model: angle of attack and speed through
# water from pitch + pressure alone.
#
# TWIN COPY: GliderTurbulence.jl carries the same flight model (its ε estimates
# scale as U⁴). The duplication is deliberate — each package stands alone — so
# physics fixes must land in BOTH files (this one and
# GliderTurbulence.jl/src/flightmodel.jl). Loading both packages requires
# qualifying the shared names (`GliderADCP.flight_model`, …).
#
# For a glider in steady, non-accelerating flight the lift/drag force balance
# ties the glide-path angle γ (from horizontal, positive up) to the angle of
# attack α (angle between body axis and flow, α = θ − γ with θ = pitch):
#
#     tan(γ) = −(C_D0 + C_D1·α²) / (a·α)                 [drag/lift polar]
#
# (Merckelbach, Smeed & Griffiths 2010, JTECH 27; Fer, Peterson & Ullgren 2014.)
# Given pitch, this is solved for α; the speed through water along the flight
# path then follows from the MEASURED pressure fall rate:
#
#     U = w_up / sin(γ),   w_up = −dP/dt   (dbar ≈ m)
#
# Using measured dP/dt avoids calibrating ballast volume, hull compressibility
# and displacement — at the cost of assuming zero vertical water velocity
# (adequate away from strong convection/internal-wave events). Only the polar
# SHAPE (C_D0/a, C_D1/a) enters, and only through the α correction (~3–8°).
#
# Sign conventions: θ, γ positive nose-up/path-up; on descent both < 0 and
# α > 0 (nose above path); on ascent α < 0. U ≥ 0 along-path. Results are
# masked (NaN) near apogees/inflections (|γ| or |w| below thresholds) and
# across time gaps.
#
# In this package the flight model is the middle rung of the DAC ladder
# (ADCP water track → flight model → onboard DR): on the four validated
# missions its dead-reckoning lands within ~1.4 cm/s median of the ADCP water
# track with no systematic along-track bias, while the onboard model is
# 4–6.5 cm/s off and systematically anti-track (validation doc 2026-07-15).

"""
    FlightParams(; C_D0 = 0.15, C_D1 = 3.2, a = 4.0)

Lift/drag polar coefficients for the steady-state flight model, all per
radian (`C_D1` per rad², `a` = total lift-curve slope per rad). Defaults
are the pooled AD2CP calibration of SEA064 over three missions (2022 Jan
Mayen, 2022 Lofoten, 2024 NESMA; ~360k steady-flight pings, equal
mission weighting, Tanaka-method AOA regression): C_D0 = 0.150,
C_D1 = 3.18, a = 4.0. With this single polar the flight model matches the
ADCP-measured through-water speed at ×1.000–×1.025 across all four
validated missions (M48 out of sample). Presets: `FLIGHT_SEA064`
(= default), `FLIGHT_SEAEXPLORER_TANAKA22` (0.20/5.0/4.0, the published
SeaExplorer calibration — a different airframe configuration),
`FLIGHT_SLOCUM_TANAKA22` (0.18/5.92/5.4), `FLIGHT_SLOCUM_MEA10`
(0.10/2.88/6.1, Merckelbach et al. 2010). Only the ratios `C_D0/a` and
`C_D1/a` affect the result; prefer a per-mission fit whenever an ADCP is
aboard: `fit_flightparams(measure_aoa(pings, nav)...)`.
"""
Base.@kwdef struct FlightParams
    C_D0::Float64 = 0.15
    C_D1::Float64 = 3.2
    a::Float64    = 4.0
end

const FLIGHT_SEA064 = FlightParams()
const FLIGHT_SEAEXPLORER = FLIGHT_SEA064
const FLIGHT_SEAEXPLORER_TANAKA22 = FlightParams(C_D0 = 0.20, C_D1 = 5.0, a = 4.0)
const FLIGHT_SLOCUM_TANAKA22 = FlightParams(C_D0 = 0.18, C_D1 = 5.92, a = 5.4)
const FLIGHT_SLOCUM_MEA10   = FlightParams(C_D0 = 0.10, C_D1 = 2.88, a = 6.1)

"""
    solve_aoa(pitch_rad::Real; params = FlightParams(), tol = 1e-10,
              maxiter = 50, aoa_max = deg2rad(15)) -> Float64

Angle of attack α (radians) for a given pitch θ (radians) from the
steady-state polar `tan(θ − α) = −(C_D0 + C_D1 α²)/(a α)`, by fixed-point
iteration with a bisection fallback. Returns NaN when no physical solution
exists: the iteration diverges, or it lands on a root with |α| > `aoa_max`
— at small pitch the polar's only roots are stall angles where the
quadratic-drag model is invalid.
"""
function solve_aoa(pitch_rad::Real; params::FlightParams = FlightParams(),
                   tol::Real = 1e-10, maxiter::Int = 50,
                   aoa_max::Real = deg2rad(15))
    isnan(pitch_rad) && return NaN
    θ = float(pitch_rad)
    abs(tan(θ)) < 1e-3 && return NaN              # near-horizontal: no steady flight
    α = 0.0
    for _ in 1:maxiter
        tanγ = tan(θ - α)
        abs(tanγ) < 1e-3 && break                 # fixed point failed — try bisection
        α′ = -(params.C_D0 + params.C_D1 * α^2) / (params.a * tanγ)
        abs(α′) > π / 4 && break                  # fixed point diverging — try bisection
        if abs(α′ - α) < tol
            return abs(α′) <= aoa_max ? α′ : NaN  # reject post-stall roots
        end
        α = α′
    end
    # The fixed-point map is only contracting for gentle polars; steeper
    # C_D1/a makes it diverge at moderate pitch even though a physical root
    # exists. Solve the residual f(α) = a·α·tan(θ − α) + C_D0 + C_D1·α² by
    # bisection: by odd symmetry work on a dive (θ < 0, α > 0), where
    # f(0⁺) = C_D0 > 0, and take the first sign change below aoa_max
    # (the root nearest α = 0 is the pre-stall branch).
    θn = -abs(θ)
    f(x) = params.a * x * tan(θn - x) + params.C_D0 + params.C_D1 * x^2
    lo = 0.0; flo = params.C_D0
    hi = NaN
    x = 0.0
    step = aoa_max / 256
    for _ in 1:256
        x += step
        fx = f(x)
        isfinite(fx) || return NaN
        if sign(fx) != sign(flo)
            hi = x
            break
        end
        lo = x; flo = fx
    end
    isnan(hi) && return NaN                        # no root ≤ aoa_max: stalled
    for _ in 1:60
        mid = (lo + hi) / 2
        (hi - lo) < tol && break
        (sign(f(mid)) == sign(flo)) ? (lo = mid; flo = f(mid)) : (hi = mid)
    end
    αp = (lo + hi) / 2
    return θ < 0 ? αp : -αp
end

"""
    GliderFlight

Flight-model solution on the pressure time base. All angles in degrees,
speeds in m/s, NaN where masked/invalid.

- `time`      : sample times
- `pressure`  : smoothed pressure/depth (dbar ≈ m, positive down)
- `w`         : vertical velocity, positive up (−dP/dt)
- `pitch`     : pitch θ interpolated & smoothed onto `time`
- `aoa`       : angle of attack α (positive on descent)
- `glide`     : glide-path angle γ = θ − α
- `U`         : speed through water along the flight path, `w / sin(γ)`
"""
struct GliderFlight
    time::Vector{DateTime}
    pressure::Vector{Float64}
    w::Vector{Float64}
    pitch::Vector{Float64}
    aoa::Vector{Float64}
    glide::Vector{Float64}
    U::Vector{Float64}
end

Base.length(f::GliderFlight) = length(f.time)

function Base.show(io::IO, f::GliderFlight)
    ok = count(!isnan, f.U)
    print(io, "GliderFlight(", length(f.time), " samples, ", ok, " valid U")
    if ok > 0
        v = filter(!isnan, f.U)
        @printf(io, ", median U=%.3f m/s", sort(v)[(length(v) + 1) ÷ 2])
    end
    print(io, ")")
end

_secs(t::Vector{DateTime}) = [Dates.value(x - t[1]) / 1000 for x in t]

# Windowed local linear fit over an irregular time series (two-pointer bounds,
# per-window least squares), skipping NaNs. Returns (value, slope) at each
# sample. Unlike a moving mean followed by differencing, the regression slope
# is unbiased where the window is truncated one-sidedly (record edges, gaps).
# Samples whose window holds fewer than `min_pts` finite points or spans less
# than `halfwidth` seconds return NaN.
function _local_linfit(ts::Vector{Float64}, v::Vector{Float64}, halfwidth::Real;
                       min_pts::Int = 3)
    n = length(v)
    val = fill(NaN, n)
    slp = fill(NaN, n)
    lo = 1; hi = 0
    for i in 1:n
        while hi < n && ts[hi + 1] <= ts[i] + halfwidth
            hi += 1
        end
        while ts[lo] < ts[i] - halfwidth
            lo += 1
        end
        st = 0.0; stt = 0.0; sv = 0.0; stv = 0.0; cnt = 0
        tmin = Inf; tmax = -Inf
        for k in lo:hi
            isnan(v[k]) && continue
            t = ts[k] - ts[i]                     # center for conditioning
            st += t; stt += t * t; sv += v[k]; stv += t * v[k]
            cnt += 1
            tmin = min(tmin, ts[k]); tmax = max(tmax, ts[k])
        end
        (cnt < min_pts || tmax - tmin < halfwidth) && continue
        det = cnt * stt - st * st
        det <= 0 && continue
        b = (cnt * stv - st * sv) / det
        a = (sv - b * st) / cnt
        val[i] = a                                # fit evaluated at ts[i]
        slp[i] = b
    end
    return val, slp
end

"""
    flight_model(t_pres, pressure, t_pitch, pitch_deg;
                 params = FlightParams(),
                 smooth = 20.0,        # smoothing half-width [s]
                 max_gap = 60.0,       # NaN across data gaps wider than this [s]
                 gamma_min = 10.0,     # mask |glide angle| < this [deg]
                 w_min = 0.02,
                 aoa_scale = 1.0) -> GliderFlight

Steady-state flight solution from a pressure record (dbar, positive down)
and a pitch record (degrees), which may be on different time bases; pitch is
interpolated onto the pressure times. Both series are smoothed with a centered
`smooth`-second half-width local linear fit before differentiating/solving,
suppressing surge oscillations and sensor noise.

Masking (NaN in `U`, `aoa`, `glide`): glide angles below `gamma_min` degrees,
|w| below `w_min` m/s (apogees, surface), unconverged α, and samples adjacent
to time gaps wider than `max_gap` seconds.

`aoa_scale` multiplies the solved angle of attack — a per-mission tuning knob;
1.0 uses the polar solution as is.
"""
function flight_model(t_pres::Vector{DateTime}, pressure::Vector{Float64},
                      t_pitch::Vector{DateTime}, pitch_deg::Vector{Float64};
                      params::FlightParams = FlightParams(),
                      smooth::Real = 20.0,
                      max_gap::Real = 60.0,
                      gamma_min::Real = 10.0,
                      w_min::Real = 0.02,
                      aoa_scale::Real = 1.0)
    n = length(t_pres)
    length(pressure) == n || throw(ArgumentError("pressure/time length mismatch"))
    length(t_pitch) == length(pitch_deg) ||
        throw(ArgumentError("pitch/time length mismatch"))
    issorted(t_pres) || throw(ArgumentError("pressure times must be sorted"))
    issorted(t_pitch) || throw(ArgumentError("pitch times must be sorted"))
    n == 0 && return GliderFlight(DateTime[], Float64[], Float64[],
                                  Float64[], Float64[], Float64[], Float64[])

    tp = _secs(t_pres)
    # Pitch → pressure time base (offset by shared origin), then a windowed
    # local linear fit of both series. The fit's slope gives w = −dP/dt
    # directly and stays unbiased where the window is one-sided.
    tθ = [Dates.value(x - t_pres[1]) / 1000 for x in t_pitch]
    θ = _interp_capped(tp, tθ, pitch_deg; max_gap)
    θfit, _ = _local_linfit(tp, θ, smooth)
    # Re-impose the interpolation's gap mask: the fit window would otherwise
    # refill up to `smooth` seconds into a masked pitch dropout.
    θs = [isnan(θ[i]) ? NaN : θfit[i] for i in 1:n]
    Ps, dPdt = _local_linfit(tp, pressure, smooth)
    w = .-dPdt

    aoa   = fill(NaN, n)
    glide = fill(NaN, n)
    U     = fill(NaN, n)
    for i in 1:n
        (isnan(θs[i]) || isnan(w[i])) && continue
        α = solve_aoa(deg2rad(θs[i]); params) * aoa_scale
        isnan(α) && continue
        γ = deg2rad(θs[i]) - α
        (abs(γ) < deg2rad(gamma_min) || abs(w[i]) < w_min) && continue
        u = w[i] / sin(γ)
        u < 0 && continue        # pitch and motion disagree (inflection): not steady
        aoa[i]   = rad2deg(α)
        glide[i] = rad2deg(γ)
        U[i]     = u
    end

    return GliderFlight(copy(t_pres), Ps, w, θs, aoa, glide, U)
end

"""
    flight_model(nav::GliderNav; kwargs...) -> GliderFlight

Convenience method from navigation alone: nav depth (m, positive down) as the
pressure record and nav pitch, both on the nav cadence (~10–20 s — coarse but
sufficient for the water-track DAC fallback; prefer the two-record method with
1 Hz CTD pressure for turbulence-grade speeds).
"""
flight_model(nav::GliderNav; kwargs...) =
    flight_model(nav.time, nav.depth, nav.time, nav.pitch; kwargs...)

# True where the buoyancy pump or battery-translation actuator is moving
# (|Δ| > thresh between consecutive finite nav records of BallastPos/LinPos in
# nav.df), dilated pad s before and pad_after s after. All-false (no screen)
# when the columns are absent. Mirrors GliderTurbulence.actuator_mask.
function _actuator_mask(nav::GliderNav; thresh::Real = 5.0, pad::Real = 30.0,
                        pad_after::Real = pad)
    n = length(nav)
    raw = falses(n)
    any_col = false
    cols = isempty(nav.df) ? String[] : names(nav.df)
    for name in ("BallastPos", "LinPos")
        name in cols || continue
        any_col = true
        col = nav.df[!, name]
        lastv = NaN
        lasti = 0
        for i in 1:n
            x = col[i]
            (x === missing || !isfinite(x)) && continue
            if lasti > 0 && abs(x - lastv) > thresh
                raw[i] = true
                raw[lasti] = true
            end
            lastv = x
            lasti = i
        end
    end
    mask = falses(n)
    any_col || return mask
    for i in findall(raw)
        lo = searchsortedfirst(nav.t, nav.t[i] - pad)
        hi = searchsortedlast(nav.t, nav.t[i] + pad_after)
        mask[lo:hi] .= true
    end
    return mask
end

"""
    measure_aoa(p::ProcessedPings, nav::GliderNav;
                o_min = 4.0, o_max = 16.0, min_cells = 3,
                pitch_min = deg2rad(8), pitch_max = deg2rad(45),
                aoa_max = deg2rad(15), act_pad = 30.0, act_settle = 60.0)
        -> NamedTuple

Measure the glider's angle of attack per steady-flight ping from the AD2CP
relative velocities, after Tanaka et al. (2022, JTECH 39:1331): the near-cell
mean relative flow is the through-water velocity (currents cancel to first
order — [`throughwater_velocity`](@ref)), γ = atan(w_tw, h_tw), α = θ_nav − γ.
Samples during/after actuator strokes are excluded when the nav table carries
`BallastPos`/`LinPos` columns. Returns steady-ping vectors
`(time, pitch, aoa, U_adcp)` (radians / m/s) ready for
[`fit_flightparams`](@ref) — the per-mission polar calibration chain:

```julia
m   = measure_aoa(pings, nav)
fit = fit_flightparams(m.pitch, m.aoa)
fl  = flight_model(nav; params = fit.params)
```
"""
function measure_aoa(p::ProcessedPings, nav::GliderNav;
                     o_min::Real = 4.0, o_max::Real = 16.0, min_cells::Int = 3,
                     pitch_min::Real = deg2rad(8), pitch_max::Real = deg2rad(45),
                     aoa_max::Real = deg2rad(15),
                     act_pad::Real = 30.0, act_settle::Real = 60.0)
    tw = throughwater_velocity(p; o_min, o_max, min_cells)
    n = length(p)
    h_tw = [isfinite(tw.u[i]) && isfinite(tw.v[i]) ? hypot(tw.u[i], tw.v[i]) : NaN
            for i in 1:n]
    w_tw = tw.w
    θ = deg2rad.(_interp_capped(p.t, nav.t, nav.pitch; max_gap = 30.0))
    am = _actuator_mask(nav; pad = act_pad, pad_after = act_settle)
    moving = _interp_capped(p.t, nav.t, Float64.(am); max_gap = 30.0)

    γ = atan.(w_tw, h_tw)
    α = θ .- γ
    ok = [i for i in 1:n if isfinite(α[i]) && isfinite(θ[i]) &&
          !(moving[i] > 0) &&
          pitch_min <= abs(θ[i]) <= pitch_max &&
          h_tw[i] > 0.05 && abs(w_tw[i]) > 0.02 && abs(α[i]) < aoa_max]
    (time = p.time[ok], pitch = θ[ok], aoa = α[ok],
     U_adcp = hypot.(h_tw[ok], w_tw[ok]))
end

"""
    fit_flightparams(pitch, aoa; a = 4.0, nsigma = 3.0, niter = 3,
                     binwidth = deg2rad(1), min_bin = 20)
        -> (params, k1, k2, n, rms, bins)

Calibrate the steady-flight polar from measured angle-of-attack samples —
e.g. AD2CP relative-flow AOA from [`measure_aoa`](@ref). With γ = θ − α, the
polar `tan γ = −(C_D0 + C_D1 α²)/(a α)` is LINEAR in the ratios k1 = C_D0/a
and k2 = C_D1/a:

    −α·tan(θ − α) = k1 + k2·α²

so ordinary least squares of y = −α·tan(θ−α) on [1, α²] (with `nsigma`·MAD
outlier rejection, `niter` rounds) seeds the ratios. That seed is then
REFINED in α-space: the samples are collapsed into `binwidth`-wide pitch bins
(median α, √n weights) and (k1, k2) are grid-searched to minimize the robust
misfit of `solve_aoa`-predicted α across the bins, with bins where the
candidate polar has no steady solution charged a stall penalty.

`pitch`/`aoa` are radians, positive up (descending flight: θ < 0, α > 0).
Only the ratios are identifiable from flight geometry — the returned
`FlightParams` expresses them at the reference lift slope `a`. `rms` is the
bin-weighted rms of the α misfit [rad]; `bins` carries the binned
(pitch, aoa, n) medians for inspection.
"""
function fit_flightparams(pitch::AbstractVector{<:Real},
                          aoa::AbstractVector{<:Real};
                          a::Real = 4.0, nsigma::Real = 3.0, niter::Int = 3,
                          binwidth::Real = deg2rad(1), min_bin::Int = 20)
    length(pitch) == length(aoa) ||
        throw(ArgumentError("pitch/aoa length mismatch"))
    ok = [i for i in eachindex(pitch)
          if isfinite(pitch[i]) && isfinite(aoa[i]) &&
             abs(aoa[i]) > 1e-4 && abs(tan(pitch[i] - aoa[i])) > 1e-3]
    length(ok) >= 10 || error("fit_flightparams: need ≥10 valid samples, got $(length(ok))")

    # ── stage 1: linear seed in y = −α·tan(θ−α) space ────────────────────────
    y0 = [-aoa[i] * tan(pitch[i] - aoa[i]) for i in ok]
    x0 = [aoa[i]^2 for i in ok]
    keep = trues(length(ok))
    k1 = NaN; k2 = NaN
    for _ in 1:niter
        X = hcat(ones(count(keep)), x0[keep])
        k1, k2 = X \ y0[keep]
        resid = y0 .- (k1 .+ k2 .* x0)
        m = median(abs.(resid[keep]))
        m > 0 || break
        newkeep = abs.(resid) .<= nsigma * 1.4826 * m
        (count(newkeep) >= 10 && newkeep != keep) || break
        keep = newkeep
    end

    # ── stage 2: α-space refinement on pitch-binned medians ─────────────────
    groups = Dict{Int, Vector{Float64}}()
    cents = Dict{Int, Vector{Float64}}()
    for (j, i) in enumerate(ok)
        keep[j] || continue
        b = floor(Int, pitch[i] / binwidth)
        push!(get!(groups, b, Float64[]), aoa[i])
        push!(get!(cents, b, Float64[]), pitch[i])
    end
    bkeys = [b for (b, v) in groups if length(v) >= min_bin]
    θb = [median(cents[b]) for b in bkeys]
    αb = [median(groups[b]) for b in bkeys]
    nb = [length(groups[b]) for b in bkeys]
    wb = sqrt.(nb)
    length(bkeys) >= 3 || error("fit_flightparams: need ≥3 populated pitch bins")
    stall_pen = deg2rad(20)                       # misfit charged when model stalls
    function cost(c1, c2)
        par = FlightParams(C_D0 = c1 * a, C_D1 = c2 * a, a = float(a))
        s = 0.0
        for i in eachindex(θb)
            αm = solve_aoa(θb[i]; params = par)
            s += wb[i] * (isnan(αm) ? stall_pen :
                          min(abs(αm - αb[i]), stall_pen))
        end
        s
    end
    c1 = max(k1, 1e-4); c2 = max(k2, 1e-3)
    span = 1.2                                    # ± e^1.2 ≈ ×3.3 first round
    for _ in 1:4
        f1 = exp.(range(-span, span, length = 13))
        f2 = exp.(range(-span, span, length = 13))
        best = (cost(c1, c2), c1, c2)
        for a1 in f1, a2 in f2
            c = cost(c1 * a1, c2 * a2)
            c < best[1] && (best = (c, c1 * a1, c2 * a2))
        end
        c1, c2 = best[2], best[3]
        span /= 3
    end
    k1, k2 = c1, c2
    par = FlightParams(C_D0 = k1 * a, C_D1 = k2 * a, a = float(a))
    res = [(αm = solve_aoa(t; params = par); isnan(αm) ? NaN : αm - x)
           for (t, x) in zip(θb, αb)]
    fin = findall(!isnan, res)
    rms = isempty(fin) ? NaN :
          sqrt(sum(wb[fin] .* res[fin] .^ 2) / sum(wb[fin]))
    (params = par, k1 = k1, k2 = k2, n = count(keep), rms = rms,
     bins = (pitch = θb, aoa = αb, n = nb))
end

# Flight-model through-water velocity series (t [unix s], u, v [m/s, E/N]) for
# the water-track DAC: horizontal speed U·cos γ projected along the TRUE nav
# heading (the gli Heading column has declination applied onboard). Assumes
# travel along heading (no sideslip) and zero vertical water velocity.
function _flight_tw(fl::GliderFlight, nav::GliderNav; max_gap::Real = 120.0)
    t = datetime2unix.(fl.time)
    s = [isfinite(h) ? sind(h) : NaN for h in nav.heading]
    c = [isfinite(h) ? cosd(h) : NaN for h in nav.heading]
    si = _interp_capped(t, nav.t, s; max_gap)
    ci = _interp_capped(t, nav.t, c; max_gap)
    u = fill(NaN, length(t))
    v = fill(NaN, length(t))
    for i in eachindex(t)
        (isfinite(fl.U[i]) && isfinite(si[i]) && isfinite(ci[i])) || continue
        nrm = hypot(si[i], ci[i])
        nrm > 0.1 || continue
        h = fl.U[i] * cosd(fl.glide[i])
        u[i] = h * si[i] / nrm
        v[i] = h * ci[i] / nrm
    end
    return (t = t, u = u, v = v)
end
