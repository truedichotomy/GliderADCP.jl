# Layers 2–3 orchestration: from a QC'd AD2CPData to solver-ready per-ping ENU
# relative velocities with absolute cell depths.

"""
    ProcessedPings

Solver-ready product of the Layer-2 chain: earth-frame water velocities **relative to
the glider** on a common vertical-offset grid, with absolute cell depths.

Fields: `time`, `t` (unix s), `depth` (glider depth per ping, m +down), `heading`
(true heading per ping, deg, declination applied), `offsets` (m, +down, length ngrid),
`E`, `N`, `U` (ngrid × ntime, m/s), `celldepth` (ngrid × ntime, m = glider depth +
offset), `look`, `beams` (per-ping 3-beam selection).
"""
struct ProcessedPings
    time::Vector{DateTime}
    t::Vector{Float64}
    depth::Vector{Float64}
    heading::Vector{Float64}
    offsets::Vector{Float64}
    E::Matrix{Float64}
    N::Matrix{Float64}
    U::Matrix{Float64}
    celldepth::Matrix{Float64}
    look::Symbol
    beams::Vector{NTuple{3,Int}}
end

Base.length(p::ProcessedPings) = length(p.t)

function Base.show(io::IO, p::ProcessedPings)
    print(io, "ProcessedPings: $(length(p)) pings × $(length(p.offsets)) offsets, ",
        "$(count(isfinite, p.E)) finite samples, look=$(p.look)")
end

"""
    glider_depth(adcp; lat=45.0) -> Vector{Float64}

Glider depth (m, positive down) from the AD2CP pressure via TEOS-10 `z_from_p`.
`lat` is a scalar or per-ping vector of latitudes.
"""
function glider_depth(a::AD2CPData; lat=45.0)
    latv = lat isa Number ? fill(Float64(lat), length(a)) : Float64.(lat)
    z = similar(a.pressure)
    for i in eachindex(z)
        z[i] = isfinite(a.pressure[i]) ? -gsw_z_from_p(a.pressure[i], latv[i], 0.0, 0.0) : NaN
    end
    return z
end

"""
    process_pings(adcp; lat=45.0, look=:auto, declination=0.0,
                  range_gating=:nominal25, method=:exact, offsets=nothing)
        -> ProcessedPings

Layer-2 chain (isobaric regrid → 3-beam ENU transform) plus absolute cell depths.
Apply `qc!` and `apply_soundspeed!` to `adcp` **before** calling. `declination` is a
scalar or per-ping vector in degrees, added to heading.
"""
function process_pings(a::AD2CPData; lat=45.0, look::Symbol=:auto, declination=0.0,
                       range_gating::Symbol=:nominal25, method::Symbol=:exact,
                       offsets=nothing)
    lk = look === :auto ? detect_look_direction(a) : look
    E, N, U, offs, used = enu_on_isobars(a; look=lk, declination, range_gating, method,
        offsets)
    gd = glider_depth(a; lat)
    celldepth = offs .+ gd'                      # ngrid × ntime
    decl(i) = declination isa Number ? Float64(declination) : Float64(declination[i])
    hdg = [isfinite(a.heading[i]) ? Float64(a.heading[i]) + decl(i) : NaN
           for i in 1:length(a)]
    ProcessedPings(copy(a.time), copy(a.t), gd, hdg, offs, E, N, U, celldepth, lk, used)
end

"""
    segment_indices(p::ProcessedPings, t_start::DateTime, t_end::DateTime) -> Vector{Int}

Ping indices within a time window (e.g. one `compute_dac` row's fix-to-fix segment).
"""
function segment_indices(p::ProcessedPings, t_start::DateTime, t_end::DateTime)
    t1, t2 = datetime2unix(t_start), datetime2unix(t_end)
    findall(t -> t1 <= t <= t2, p.t)
end

"""
    throughwater_velocity(p::ProcessedPings; o_min=4.0, o_max=16.0, min_cells=3)
        -> (u, v, w)

Glider through-water velocity per ping: **minus** the mean relative water velocity
over the near cells (`o_min`–`o_max` m along-beam offset — close enough that the
water there moves with the water at the glider, so currents cancel to first order;
Tanaka et al. 2022). Components are east/north/up (m/s), referenced like the pings
(true north when `process_pings` was given a declination); `NaN` where fewer than
`min_cells` cells are finite. This is a direct measurement — no flight model — and
is what the water-track [`compute_dac`](@ref) method integrates. The residual is
the mean shear across the cell offset (≲ 1 cm/s for typical stratification over
the default 4–16 m window).
"""
function throughwater_velocity(p::ProcessedPings; o_min::Real=4.0, o_max::Real=16.0,
                               min_cells::Int=3)
    sel = findall(o -> o_min <= abs(o) <= o_max, p.offsets)
    n = length(p)
    u = fill(NaN, n); v = fill(NaN, n); w = fill(NaN, n)
    for i in 1:n
        se = sn = su = 0.0
        ne = nn = nu = 0
        for j in sel
            isfinite(p.E[j, i]) && (se += p.E[j, i]; ne += 1)
            isfinite(p.N[j, i]) && (sn += p.N[j, i]; nn += 1)
            isfinite(p.U[j, i]) && (su += p.U[j, i]; nu += 1)
        end
        if ne >= min_cells && nn >= min_cells
            u[i] = -se / ne
            v[i] = -sn / nn
        end
        nu >= min_cells && (w[i] = -su / nu)
    end
    return (u=u, v=v, w=w)
end
