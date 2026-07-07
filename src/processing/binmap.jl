# Layer 2 — cell geometry: along-beam ranges, vertical offsets, isobaric regridding.
#
# Nortek firmware fact (gliderad2cp docs, pers. comm. S. Nylund 2024; see research docs):
# with >2 beams enabled the instrument time-gates ALL beams assuming a nominal 25° slant,
# so the along-beam cell distance is `range/cos(25°)` for every beam (range_gating =
# :nominal25, default). :alongbeam treats the reported ranges as already along-beam
# (Slocum-AD2CP parity).
#
# Vertical offsets are positive DOWN (below the glider): off = −r · (R·F·e_beam)_z.
# Cell depth = glider depth + offset. Regridding interpolates each beam's velocities
# onto a common vertical-offset grid per ping (gliderad2cp's isobaric regrid), so that
# the subsequent 3-beam transform combines measurements from a common depth.

_alongbeam_ranges(a::AD2CPData, range_gating::Symbol) =
    range_gating === :nominal25 ? a.range ./ cosd(25.0) :
    range_gating === :alongbeam ? copy(a.range) :
    error("range_gating must be :nominal25 or :alongbeam")

"""
    vertical_cosines(adcp; look, declination=0.0) -> 4 × ntime Matrix

Per-ping, per-beam vertical direction cosines: the downward component of each beam's
unit vector in the earth frame (positive when the beam points down).
"""
function vertical_cosines(a::AD2CPData; look::Symbol, declination=0.0)
    e = beam_unit_vectors(a.config)
    F = head2vehicle(look)
    nt = length(a)
    vc = fill(NaN, 4, nt)
    decl(i) = declination isa Number ? Float64(declination) : Float64(declination[i])
    for i in 1:nt
        (isfinite(a.pitch[i]) && isfinite(a.roll[i]) && isfinite(a.heading[i])) || continue
        R = rotmat_xyz2enu(a.heading[i], a.pitch[i], a.roll[i]; declination=decl(i)) * F
        for b in 1:4
            vc[b, i] = -(R*e[b])[3]
        end
    end
    return vc
end

"""
    offset_grid(cfg; spacing=cfg.cellsize/2, look=:down) -> Vector

Common vertical-offset grid (m, signed positive down): covers blanking to
`blanking + (ncells + 0.5)·cellsize`, at `spacing` resolution; negated for up-looking.
"""
function offset_grid(cfg::AD2CPConfig; spacing::Real=cfg.cellsize / 2, look::Symbol=:down)
    maxd = cfg.blanking + (cfg.ncells + 0.5) * cfg.cellsize
    g = collect(0.0:spacing:maxd)
    return look === :down ? g : -g
end

"""
    regrid_beams(adcp; look, declination=0.0, range_gating=:nominal25,
                 offsets=offset_grid(adcp.config; look)) -> (V, offsets, vc)

Interpolate each beam's velocities from their per-ping tilt-dependent vertical offsets
onto the common `offsets` grid. Returns `V` (`length(offsets) × 4 × ntime`, NaN outside
each beam's span), the grid, and the vertical cosines used.
"""
function regrid_beams(a::AD2CPData; look::Symbol, declination=0.0,
                      range_gating::Symbol=:nominal25,
                      offsets::AbstractVector=offset_grid(a.config; look))
    r = _alongbeam_ranges(a, range_gating)
    vc = vertical_cosines(a; look, declination)
    nt = length(a)
    ng = length(offsets)
    V = fill(NaN32, ng, 4, nt)
    x = similar(r)
    for i in 1:nt, b in 1:4
        # skip near-horizontal beams (degenerate vertical mapping)
        (isfinite(vc[b, i]) && abs(vc[b, i]) > 0.05) || continue
        @. x = r * vc[b, i]                      # this beam's cell offsets (signed)
        y = @view a.vel[:, b, i]
        # x is monotone (vc const per ping); ensure ascending for interpolation
        xi, yi = vc[b, i] >= 0 ? (x, y) : (reverse(x), reverse(y))
        V[:, b, i] = _interp1(xi, Float64.(yi), Float64.(offsets))
    end
    return V, collect(Float64, offsets), vc
end

"""
    enu_on_isobars(adcp; look=:auto, declination=0.0, range_gating=:nominal25,
                   method=:exact, offsets=nothing)
        -> (E, N, U, offsets, used)

Full Layer-2 chain: isobaric regrid of beam velocities, then per-ping 3-beam
transformation to earth-frame relative velocities on the common offset grid.
Returns `noffsets × ntime` matrices (relative water velocity, m/s), the offset grid
(m, positive down), and the per-ping beam selection.
"""
function enu_on_isobars(a::AD2CPData; look::Symbol=:auto, declination=0.0,
                        range_gating::Symbol=:nominal25, method::Symbol=:exact,
                        offsets=nothing)
    lk = look === :auto ? detect_look_direction(a) : look
    offs = offsets === nothing ? offset_grid(a.config; look=lk) : collect(Float64, offsets)
    V, offs, _ = regrid_beams(a; look=lk, declination, range_gating, offsets=offs)
    e = beam_unit_vectors(a.config)
    F = head2vehicle(lk)
    ng, _, nt = size(V)
    E = fill(NaN, ng, nt); N = fill(NaN, ng, nt); U = fill(NaN, ng, nt)
    used = Vector{NTuple{3,Int}}(undef, nt)
    decl(i) = declination isa Number ? Float64(declination) : Float64(declination[i])
    S = Dict(sel => inv(beams_matrix(e, sel)) for sel in ((1, 2, 4), (2, 3, 4)))
    for i in 1:nt
        h, p, r = a.heading[i], a.pitch[i], a.roll[i]
        if !(isfinite(h) && isfinite(p) && isfinite(r))
            used[i] = (0, 0, 0); continue
        end
        sel = select_beams(p; look=lk)
        used[i] = sel
        R = rotmat_xyz2enu(h, p, r; declination=decl(i)) * F
        for k in 1:ng
            bk = (V[k, sel[1], i], V[k, sel[2], i], V[k, sel[3], i])
            all(isfinite, bk) || continue
            v = method === :exact ? S[sel] * collect(Float64, bk) :
                xyz_from_beams(collect(bk), sel, e; method)
            w = R * v
            E[k, i], N[k, i], U[k, i] = w[1], w[2], w[3]
        end
    end
    return E, N, U, offs, used
end
