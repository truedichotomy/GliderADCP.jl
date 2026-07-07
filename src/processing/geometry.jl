# Layer 2 — beam geometry and coordinate transforms (the heart of the package).
#
# Frames and conventions (derived from first principles; consistent with the factory
# beam→XYZ matrix, the Nortek transform convention, and the paper-validated behavior of
# gliderad2cp / Slocum-AD2CP — see docs/research/*):
#
#   * Head frame (factory): right-handed, X toward beam-1 azimuth, Z up out of the
#     transducer face, Y = Z × X. Beam *away* unit vectors:
#       e_i = (sinθ_i cosφ_i, sinθ_i sinφ_i, cosθ_i)
#     Nortek beam velocities are positive AWAY from the transducer: b_i = e_i ⋅ v.
#     (Verified: with these conventions the least-squares X,Y rows equal the factory
#     beam2xyz X,Y rows, e.g. 1/(2 sin47.5°)=0.6782.)
#   * Vehicle frame: X forward, Z up, Y = Z × X. A head mounted looking DOWN is rotated
#     180° about X relative to the vehicle: (x, y, z)_vehicle = (x, −y, −z)_head — the
#     origin of the "flip rows Y and Z for downward-looking" rule in both Python packages.
#   * Earth frame: ENU. v_enu = R(heading, pitch, roll) ⋅ v_vehicle with R = H·P:
#       hh = deg2rad(heading + declination − 90);  pp, rr = deg2rad(pitch), deg2rad(roll)
#       H = [ cos hh  sin hh  0;  −sin hh  cos hh  0;  0  0  1 ]
#       P = [ cos pp  −sin pp sin rr  −cos rr sin pp;
#             0        cos rr         −sin rr;
#             sin pp   sin rr cos pp   cos pp cos rr ]
#     (The validated form; Slocum-AD2CP v2.0.0 accidentally applies Rᵀ — see research doc.)
#
# 3-beam solution: the fore/aft beam nearest horizontal is dropped
# (down-looking: dive keeps {1,2,4}, climb keeps {2,3,4}; mirrored for up-looking) and
# the remaining 3×3 system E v = b is solved EXACTLY. This equals gliderad2cp's
# "synthesize the dropped beam with zero error velocity" trick; Slocum-AD2CP's
# row-dropped factory submatrix is *not* the exact inverse (halves X and Z and leaves
# vertical-velocity leakage in X) — available here as `method=:rowdrop` for parity.

"""
    beam_unit_vectors(θ, φ) -> Vector{SVector-like}

Beam *away* unit vectors in the factory head frame (Z up out of the face), from beam
inclinations `θ` and azimuths `φ` in degrees (4-tuples). Returns a vector of four
3-element `Vector{Float64}`.
"""
function beam_unit_vectors(θ::NTuple{4,<:Real}, φ::NTuple{4,<:Real})
    [Float64[sind(θ[i]) * cosd(φ[i]), sind(θ[i]) * sind(φ[i]), cosd(θ[i])] for i in 1:4]
end

beam_unit_vectors(cfg::AD2CPConfig) = beam_unit_vectors(cfg.beam_theta, cfg.beam_phi)

"""
    beams_matrix(e, beams) -> Matrix

Rows `e[i]'` for the selected beams: the forward model `b = E ⋅ v_head`.
"""
beams_matrix(e::AbstractVector, beams) = permutedims(hcat((e[i] for i in beams)...))

"""
    xyz_from_beams(b, beams, e; method=:exact) -> Vector (3)

Head-frame velocity from along-beam velocities `b` (values for `beams`, positive away).
`method=:exact` solves the selected 3×3 (or least-squares 4×3) system `E v = b`.
`method=:rowdrop` reproduces Slocum-AD2CP's factory-submatrix behavior (parity only).
"""
function xyz_from_beams(b::AbstractVector, beams, e::AbstractVector; method::Symbol=:exact)
    E = beams_matrix(e, beams)
    if method === :exact
        return length(beams) == 3 ? E \ collect(Float64, b) :
               (E' * E) \ (E' * collect(Float64, b))
    elseif method === :rowdrop
        # factory rows X,Y,Z1 with dropped-beam column removed (Slocum-AD2CP v2 behavior)
        θs = (47.5, 25.0, 47.5, 25.0)
        a(t) = 1 / (2 * sind(t)); c(t) = 1 / (2 * cosd(t))
        M = [a(θs[1]) 0 -a(θs[3]) 0; 0 -a(θs[2]) 0 a(θs[4]); c(θs[1]) 0 c(θs[3]) 0]
        return M[:, collect(beams)] * collect(Float64, b)
    else
        error("xyz_from_beams: unknown method $method")
    end
end

"""
    select_beams(pitch; look=:down) -> NTuple{3,Int}

Which three beams to use for a ping: drop the fore/aft (47.5°) beam that is furthest
from vertical. Down-looking glider mount: dive (pitch < 0) keeps (1,2,4), climb keeps
(2,3,4). Up-looking mount is mirrored.
"""
function select_beams(pitch::Real; look::Symbol=:down)
    dive = pitch < 0
    if look === :down
        return dive ? (1, 2, 4) : (2, 3, 4)
    elseif look === :up
        return dive ? (2, 3, 4) : (1, 2, 4)
    end
    error("select_beams: look must be :down or :up")
end

"""
    head2vehicle(look) -> Diagonal-like 3×3

Head→vehicle frame rotation: identity for an up-looking head, 180° about X
(flip Y and Z) for a down-looking head.
"""
head2vehicle(look::Symbol) =
    look === :up ? Float64[1 0 0; 0 1 0; 0 0 1] :
    look === :down ? Float64[1 0 0; 0 -1 0; 0 0 -1] :
    error("head2vehicle: look must be :down or :up")

"""
    detect_look_direction(adcp) -> :down | :up

Mounting direction from the median Z accelerometer (gravity): negative ⇒ head faces
down (same heuristic as gliderad2cp, inverted sign convention checked on data:
M38 down-looking has median(AccelerometerZ) < 0).
"""
function detect_look_direction(a::AD2CPData)
    az = nanmedian(a.accel[3, :])
    isnan(az) && error("detect_look_direction: no AccelerometerZ data — pass look explicitly")
    return az < 0 ? :down : :up
end

"""
    rotmat_xyz2enu(heading, pitch, roll; declination=0.0) -> 3×3 Matrix

Vehicle→ENU rotation `R = H·P` (angles in degrees; declination added to heading).
"""
function rotmat_xyz2enu(heading::Real, pitch::Real, roll::Real; declination::Real=0.0)
    hh = deg2rad(heading + declination - 90)
    pp = deg2rad(pitch)
    rr = deg2rad(roll)
    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1]
    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);
         0 cos(rr) -sin(rr);
         sin(pp) sin(rr)*cos(pp) cos(pp)*cos(rr)]
    return H * P
end

"""
    beams_to_enu(adcp; look=:auto, declination=0.0, method=:exact)
        -> (E, N, U, beams_used)

Transform beam velocities to earth-frame relative velocities per ping.
Returns `ncells × ntime` matrices `E, N, U` (water velocity relative to glider, m/s)
and the per-ping beam selection (`Vector{NTuple{3,Int}}`).

`declination` may be a scalar or a per-ping vector (deg, added to heading).
"""
function beams_to_enu(a::AD2CPData; look::Symbol=:auto, declination=0.0,
                      method::Symbol=:exact)
    a.config.coordsystem === :beam ||
        error("beams_to_enu: data are in $(a.config.coordsystem) coordinates, expected :beam")
    lk = look === :auto ? detect_look_direction(a) : look
    e = beam_unit_vectors(a.config)
    F = head2vehicle(lk)
    nc, _, nt = size(a.vel)
    E = fill(NaN, nc, nt); N = fill(NaN, nc, nt); U = fill(NaN, nc, nt)
    used = Vector{NTuple{3,Int}}(undef, nt)
    decl(i) = declination isa Number ? Float64(declination) : Float64(declination[i])
    # precompute per-selection head-frame solve matrices: v_head = S ⋅ b
    sel_dive, sel_climb = (1, 2, 4), (2, 3, 4)
    S = Dict(sel => inv(beams_matrix(e, sel)) for sel in (sel_dive, sel_climb))
    for i in 1:nt
        p = a.pitch[i]
        h = a.heading[i]
        r = a.roll[i]
        if !isfinite(p) || !isfinite(h) || !isfinite(r)
            used[i] = (0, 0, 0)
            continue
        end
        sel = select_beams(p; look=lk)
        used[i] = sel
        R = rotmat_xyz2enu(h, p, r; declination=decl(i)) * F
        b1, b2, b3 = a.vel[:, sel[1], i], a.vel[:, sel[2], i], a.vel[:, sel[3], i]
        for k in 1:nc
            bk = (b1[k], b2[k], b3[k])
            all(isfinite, bk) || continue
            v = method === :exact ? S[sel] * collect(Float64, bk) :
                xyz_from_beams(collect(bk), sel, e; method)
            w = R * v
            E[k, i], N[k, i], U[k, i] = w[1], w[2], w[3]
        end
    end
    return E, N, U, used
end
