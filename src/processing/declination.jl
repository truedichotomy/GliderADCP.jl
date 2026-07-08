# Layer 2/3 — magnetic declination from the IGRF model.
#
# The AD2CP compass reports magnetic heading; ENU velocities need true heading. The
# reference mission (and typical deployments) leave the instrument's declination unset,
# so processing must supply it: user value or IGRF (SatelliteToolboxGeomagneticField).

"""
    magnetic_declination(lat, lon, time::DateTime; alt=0.0) -> Float64

IGRF magnetic declination (degrees East of true north) at a position and time.
"""
function magnetic_declination(lat::Real, lon::Real, time::DateTime; alt::Real=0.0)
    y = year(time)
    doy = dayofyear(time)
    ndays = daysinyear(y)
    decyear = y + (doy - 1 + (hour(time) / 24)) / ndays
    b = igrfd(decyear, alt, lat, lon, Val(:geodetic))   # NED, nT
    return atand(b[2], b[1])
end

"""
    magnetic_declination(nav::GliderNav, t::AbstractVector; every=3600.0) -> Vector

Per-ping declination for a mission: glider positions interpolated onto times `t`
(unix s), IGRF evaluated at `every`-second knots and linearly interpolated between
(declination varies slowly). NaN where the track is unknown.
"""
function magnetic_declination(nav::GliderNav, t::AbstractVector; every::Real=3600.0)
    lat = _interp1(nav.t, nav.lat, t)
    lon = _interp1(nav.t, nav.lon, t)
    isempty(t) && return Float64[]
    knots = collect(range(minimum(t), maximum(t); step=every))
    knots[end] < maximum(t) && push!(knots, maximum(t))
    latk = _interp1(nav.t, nav.lat, knots)
    lonk = _interp1(nav.t, nav.lon, knots)
    dk = fill(NaN, length(knots))
    for i in eachindex(knots)
        (isfinite(latk[i]) && isfinite(lonk[i])) || continue
        dk[i] = magnetic_declination(latk[i], lonk[i], unix2datetime(knots[i]))
    end
    out = _interp1(knots, dk, collect(Float64, t))
    # constant-extrapolate at the edges: declination varies slowly, and NaN declination
    # would silently drop every affected ping from the ENU transform downstream
    fin = findall(isfinite, out)
    if !isempty(fin) && length(fin) < length(out)
        nfix = 0
        first_v = out[fin[1]]; last_v = out[fin[end]]
        for i in 1:fin[1]-1
            out[i] = first_v; nfix += 1
        end
        for i in fin[end]+1:length(out)
            out[i] = last_v; nfix += 1
        end
        nfix > 0 && @warn "magnetic_declination: constant-extrapolated $nfix pings " *
                          "outside navigation coverage"
    end
    return out
end
