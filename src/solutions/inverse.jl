# Layer 4b — least-squares inverse method (Visbeck 2002; Todd et al. 2017;
# Gradone et al. 2023).
#
# Per segment, unknowns m = [u_glider(1..nt) ; u_ocean(1..nz)] with the measurement
# model  u_ocean(z_cell) − u_glider(t_ping) = v_relative.  U and V decouple for a real
# design matrix, so the complex problem is solved as two real least-squares problems
# sharing one sparse QR factorization. Constraint rows (all weighted, all optional):
#
#   DAC (:ocean form, Gradone):    (1/H) Σ dz·u_ocean(k) = DAC   over glider-depth bins
#   DAC (:platform form, Visbeck): (1/T) Σ Δt·u_glider(j) = DAC  (GPS displacement / T)
#   bottom track:                  u_glider(j) = v_bt(j)          (over-ground velocity)
#   smoothness:                    interior second differences on the ocean (and
#                                  optionally glider) block
#
# Differences from Slocum-AD2CP v2.0.0 (see docs/research/slocum_ad2cp_analysis.md):
# exact 3-beam transform upstream, correct rotation, no off-by-one bin drop, clean
# interior-only D2 rows, and bottom-track support.

Base.@kwdef struct InverseOptions
    dz::Float64 = 10.0              # ocean bin size (m); ≥ cell size
    wdac::Float64 = 5.0             # DAC constraint weight (Todd et al. 2017)
    wsmooth_ocean::Float64 = 1.0    # ocean-profile curvature weight (Todd et al. 2017)
    wsmooth_glider::Float64 = 0.0   # glider-velocity curvature weight (off by default)
    wbt::Float64 = 5.0              # bottom-track row weight
    bt_max_dt::Float64 = 6.0        # s, BT fix → ping association tolerance
    dac_form::Symbol = :ocean       # :ocean | :ocean_timeweighted | :platform
    min_pings::Int = 30             # skip segments with fewer usable pings
    min_bin_obs::Int = 1            # trim edge bins with fewer observations
end

"""
    invert_segment(E, N, celldepth, tping, maxgliderdepth;
                   dacu=NaN, dacv=NaN, bt=nothing, opts=InverseOptions())

Solve one segment. `E`, `N`, `celldepth` are `ngrid × nt` (relative velocities and
absolute cell depths); `tping` unix seconds; `maxgliderdepth` limits the DAC-covered
bins. `bt` is an optional table with columns `t, u, v` (glider over-ground velocity,
from [`bt_velocity`](@ref)). Returns a NamedTuple
`(z, u, v, nobs, ug, vg, tping, nbt, resid)` or `nothing` if the segment is unusable.
"""
function invert_segment(E::AbstractMatrix, N::AbstractMatrix, celldepth::AbstractMatrix,
                        tping::AbstractVector, maxgliderdepth::Real;
                        dacu::Real=NaN, dacv::Real=NaN, bt=nothing, gliderdepth=nothing,
                        opts::InverseOptions=InverseOptions())
    ngrid, nt = size(E)
    dz = opts.dz

    # --- pass 1: bin census -----------------------------------------------------------
    maxk = 0
    for i in 1:nt, k in 1:ngrid
        z = celldepth[k, i]
        (isfinite(E[k, i]) && isfinite(N[k, i]) && isfinite(z) && z >= 0) || continue
        maxk = max(maxk, floor(Int, z / dz) + 1)
    end
    maxk == 0 && return nothing
    counts = zeros(Int, maxk)
    for i in 1:nt, k in 1:ngrid
        z = celldepth[k, i]
        (isfinite(E[k, i]) && isfinite(N[k, i]) && isfinite(z) && z >= 0) || continue
        counts[floor(Int, z / dz)+1] += 1
    end
    kmin = findfirst(>=(opts.min_bin_obs), counts)
    kmax = findlast(>=(opts.min_bin_obs), counts)
    kmin === nothing && return nothing
    nz = kmax - kmin + 1
    nz < 2 && return nothing
    nobs = counts[kmin:kmax]
    if any(==(0), nobs) && opts.wsmooth_ocean <= 0
        error("invert_segment: interior ocean bins without observations require " *
              "wsmooth_ocean > 0 (or coarser dz)")
    end

    # --- assemble sparse triplets -----------------------------------------------------
    Is = Int[]; Js = Int[]; Vs = Float64[]
    du = Float64[]; dv = Float64[]
    r = 0
    addentry(i, j, v) = (push!(Is, i); push!(Js, j); push!(Vs, v))

    nused = 0
    for i in 1:nt, k in 1:ngrid
        z = celldepth[k, i]
        (isfinite(E[k, i]) && isfinite(N[k, i]) && isfinite(z) && z >= 0) || continue
        kb = floor(Int, z / dz) + 1
        (kmin <= kb <= kmax) || continue
        r += 1
        addentry(r, i, -1.0)                    # glider velocity column
        addentry(r, nt + (kb - kmin + 1), 1.0)  # ocean bin column
        push!(du, E[k, i]); push!(dv, N[k, i])
        nused += 1
    end
    r == 0 && return nothing

    # --- DAC constraint ---------------------------------------------------------------
    if isfinite(dacu) && isfinite(dacv) && opts.wdac > 0
        if opts.dac_form === :ocean || opts.dac_form === :ocean_timeweighted
            covered = [k for k in 1:nz if (kmin + k - 2 + 0.5) * dz <= maxgliderdepth &&
                       nobs[k] > 0]
            weights = if opts.dac_form === :ocean_timeweighted && gliderdepth !== nothing
                # DAC is a time average: weight each bin by glider residence time
                zcent = [(kmin + k - 2 + 0.5) * dz for k in 1:nz]
                w = time_in_bin(gliderdepth, tping, zcent; dz)[covered]
                sum(w) > 0 ? w ./ sum(w) : fill(1.0 / length(covered), length(covered))
            else
                fill(1.0 / max(1, length(covered)), length(covered))   # plain depth mean
            end
            if !isempty(covered)
                Cn = 1 / norm(weights)                   # unit-norm row scaling
                r += 1
                for (j, k) in enumerate(covered)
                    addentry(r, nt + k, opts.wdac * Cn * weights[j])
                end
                push!(du, opts.wdac * Cn * dacu)
                push!(dv, opts.wdac * Cn * dacv)
            end
        elseif opts.dac_form === :platform
            dt = diff(tping)
            mdt = median(dt)
            dts = [clamp(i == nt ? mdt : dt[i], 0, 10mdt) for i in 1:nt]
            T = sum(dts)
            Cn = 1 / sqrt(sum((dts ./ T) .^ 2))
            r += 1
            for j in 1:nt
                addentry(r, j, opts.wdac * Cn * dts[j] / T)
            end
            push!(du, opts.wdac * Cn * dacu)
            push!(dv, opts.wdac * Cn * dacv)
        else
            error("invert_segment: dac_form must be :ocean or :platform")
        end
    end

    # --- bottom-track constraints -----------------------------------------------------
    nbt = 0
    if bt !== nothing && opts.wbt > 0
        for row in eachrow(bt)
            j = searchsortedfirst(tping, row.t)
            best = 0; bestdt = opts.bt_max_dt
            for jj in (j - 1, j)
                if 1 <= jj <= nt && abs(tping[jj] - row.t) <= bestdt
                    best = jj; bestdt = abs(tping[jj] - row.t)
                end
            end
            best == 0 && continue
            (isfinite(row.u) && isfinite(row.v)) || continue
            r += 1
            addentry(r, best, opts.wbt)
            push!(du, opts.wbt * row.u)
            push!(dv, opts.wbt * row.v)
            nbt += 1
        end
    end

    # --- smoothness (interior second differences) --------------------------------------
    if opts.wsmooth_ocean > 0
        for k in 2:nz-1
            r += 1
            addentry(r, nt + k - 1, -opts.wsmooth_ocean)
            addentry(r, nt + k, 2opts.wsmooth_ocean)
            addentry(r, nt + k + 1, -opts.wsmooth_ocean)
            push!(du, 0.0); push!(dv, 0.0)
        end
    end
    if opts.wsmooth_glider > 0
        for j in 2:nt-1
            r += 1
            addentry(r, j - 1, -opts.wsmooth_glider)
            addentry(r, j, 2opts.wsmooth_glider)
            addentry(r, j + 1, -opts.wsmooth_glider)
            push!(du, 0.0); push!(dv, 0.0)
        end
    end

    A = sparse(Is, Js, Vs, r, nt + nz)
    F = qr(A)
    xu = F \ du
    xv = F \ dv

    z = [(kmin + k - 2 + 0.5) * dz for k in 1:nz]
    resid = (norm(A * xu - du) + norm(A * xv - dv)) / max(1, nused)
    return (z=z, u=xu[nt+1:end], v=xv[nt+1:end], nobs=nobs,
        ug=xu[1:nt], vg=xv[1:nt], tping=collect(tping), nbt=nbt, resid=resid)
end

"""
    solve_inverse(pings::ProcessedPings, dac::DataFrame;
                  bt=nothing, opts=InverseOptions()) -> DataFrame

Run the inverse per DAC segment (fix-to-fix yo windows from [`compute_dac`](@ref)).
`bt` is the output of [`bt_velocity`](@ref) (optional). Returns a long-format table:
`yo, t_mid, z, u, v, nobs, nbt` — one row per ocean depth bin per segment
(the same shape as the reference `absolute_ocean_vel.csv`).
"""
function solve_inverse(p::ProcessedPings, dac::DataFrame;
                       bt=nothing, opts::InverseOptions=InverseOptions())
    out = DataFrame(yo=Int[], t_mid=DateTime[], z=Float64[], u=Float64[], v=Float64[],
        nobs=Int[], nbt=Int[])
    for row in eachrow(dac)
        idx = segment_indices(p, row.t_start, row.t_end)
        length(idx) >= opts.min_pings || continue
        btseg = nothing
        if bt !== nothing
            t1, t2 = datetime2unix(row.t_start), datetime2unix(row.t_end)
            btseg = bt[(bt.t .>= t1) .& (bt.t .<= t2), :]
            isempty(btseg) && (btseg = nothing)
        end
        gd = filter(isfinite, p.depth[idx])
        isempty(gd) && continue
        sol = invert_segment(view(p.E, :, idx), view(p.N, :, idx),
            view(p.celldepth, :, idx), p.t[idx], maximum(gd);
            dacu=row.u, dacv=row.v, bt=btseg, gliderdepth=p.depth[idx], opts)
        sol === nothing && continue
        for k in eachindex(sol.z)
            push!(out, (row.yo, row.t_mid, sol.z[k], sol.u[k], sol.v[k],
                sol.nobs[k], sol.nbt))
        end
    end
    return out
end
