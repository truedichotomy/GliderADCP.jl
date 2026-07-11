# Layer 2 — sound speed correction.
#
# Along-beam velocities scale linearly with the sound speed at the transducer
# (Fischer & Visbeck 1993): v_true = v_recorded · c_true/c_used. The AD2CP computes
# c_used from its own temperature and a CONFIGURED salinity (M38: 38 vs actual ≈35.05
# → ≈0.3 % velocity bias). c_true comes from CTD via TEOS-10.

"""
    _interp1(x, y, xi) -> Vector

NaN-tolerant 1-D linear interpolation (finite pairs only; NaN outside the data range).
`x` must be sorted ascending.
"""
function _interp1(x::AbstractVector, y::AbstractVector, xi::AbstractVector)
    good = findall(i -> isfinite(x[i]) && isfinite(y[i]), eachindex(x, y))
    out = fill(NaN, length(xi))
    length(good) < 2 && return out
    xg = Float64.(x[good]); yg = Float64.(y[good])
    # Gridded interpolation needs strictly increasing knots
    keep = [1; findall(i -> xg[i] > xg[i-1], 2:length(xg)) .+ 1]
    length(keep) < 2 && return out
    xg = xg[keep]; yg = yg[keep]
    itp = interpolate((xg,), yg, Gridded(Linear()))
    for (k, v) in pairs(xi)
        (isfinite(v) && xg[1] <= v <= xg[end]) || continue
        out[k] = itp(v)
    end
    return out
end

"""
    soundspeed_from_ctd(SP, t, p, lon, lat) -> c

TEOS-10 sound speed (m/s) from practical salinity, in-situ temperature (°C),
pressure (dbar) and position (broadcastable).
"""
function soundspeed_from_ctd(SP, t, p, lon, lat)
    SA = gsw_sa_from_sp.(SP, p, lon, lat)
    CT = gsw_ct_from_t.(SA, t, p)
    return gsw_sound_speed.(SA, CT, p)
end

"""
    onboard_soundspeed!(adcp, ctd_t, ctd_temperature; salinity, lat, lon=0.0) -> Vector

Reconstruct the per-ping sound speed the instrument used **onboard** and write it into
`adcp.soundspeed` (returned). Needed for the realtime-telemetered route
([`load_pld_adcp`](@ref)), which does not transmit the onboard value: the AD2CP
computes it from its own temperature and the **configured** `salinity` (from the
deployment plan; e.g. 38.0 on the reference mission). The instrument temperature is
approximated by the payload CTD temperature (`ctd_t` unix seconds,
`ctd_temperature` °C) interpolated to the pings — a ≲0.1 % residual — with the
instrument's own transmitted pressure. After this, the standard
[`soundspeed_correction`](@ref) → [`apply_soundspeed!`](@ref) chain applies unchanged:

```julia
tele = load_pld_adcp(srcs; stream="38.pld1.sub", cellsize=2.0, blanking=0.7)
onboard_soundspeed!(tele, ctd_t, ctd_T; salinity=38.0, lat=69.5)
apply_soundspeed!(tele, soundspeed_correction(tele, ctd_t, ctd_c_true))
```
"""
function onboard_soundspeed!(a::AD2CPData, ctd_t::AbstractVector, ctd_T::AbstractVector;
                             salinity::Real, lat::Real, lon::Real=0.0)
    Tping = _interp1(ctd_t, ctd_T, a.t)
    c = soundspeed_from_ctd(salinity, Tping, a.pressure, lon, lat)
    a.soundspeed .= Float32.(c)
    return a.soundspeed
end

"""
    soundspeed_correction(adcp, ctd_t, ctd_c) -> Vector

Per-ping velocity scale factor `c_true/c_recorded`, with `c_true` linearly interpolated
from a CTD sound-speed time series (`ctd_t` in unix seconds, `ctd_c` in m/s) onto the
ADCP pings. NaN outside CTD coverage.
"""
function soundspeed_correction(a::AD2CPData, ctd_t::AbstractVector, ctd_c::AbstractVector)
    c_true = _interp1(ctd_t, ctd_c, a.t)
    return c_true ./ Float64.(a.soundspeed)
end

"""
    apply_soundspeed!(adcp, scale::AbstractVector; scale_bt=true)

Scale beam velocities (and bottom-track velocities, matched by nearest ping) by the
per-ping factor `scale`. Pings with non-finite `scale` are left untouched (use QC to
mask them if desired). Mutates `adcp.vel` (and `adcp.bt.vel`) in place; returns the
number of pings scaled.
"""
function apply_soundspeed!(a::AD2CPData, scale::AbstractVector; scale_bt::Bool=true)
    n = 0
    for i in eachindex(scale)
        s = scale[i]
        isfinite(s) || continue
        @views a.vel[:, :, i] .*= Float32(s)
        n += 1
    end
    if scale_bt && a.bt !== nothing
        sbt = _interp1(a.t, collect(scale), a.bt.t)
        for i in eachindex(sbt)
            isfinite(sbt[i]) || continue
            @views a.bt.vel[:, i] .*= Float32(sbt[i])
        end
    end
    return n
end
