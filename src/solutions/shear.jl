# Layer 4a — shear method (lADCP tradition: Fischer & Visbeck 1993; Thurnherr et al.
# 2015; gliderad2cp).
#
#   1. per-ping vertical shear of the ENU relative velocities across adjacent offsets
#      (glider motion is constant within a ping, so relative shear = ocean shear);
#      computed against the *physical* signed offset spacing (this deliberately differs
#      from gliderad2cp's per-grid-index derivative — see research doc §E1)
#   2. bin shear samples into depth bins (median by default) per segment
#   3. integrate top-down → baroclinic profile (gaps hold velocity constant)
#   4. reference: depth-mean over glider-sampled bins matched to the segment DAC

Base.@kwdef struct ShearOptions
    dz::Float64 = 10.0              # depth bin size (m)
    stat::Symbol = :median          # :median | :mean bin statistic
    min_bin_obs::Int = 4            # bins with fewer shear samples are treated as gaps
    min_pings::Int = 30
    referencing::Symbol = :timeweighted   # :timeweighted | :simple
end

"""
    time_in_bin(depth, t, z; dz) -> Vector

Seconds the glider spent in each depth bin (centers `z`, width `dz`), from the per-ping
glider depths and times. Ping intervals are capped at 10× the median (surfacing gaps).
"""
function time_in_bin(depth::AbstractVector, t::AbstractVector, z::AbstractVector; dz::Real)
    w = zeros(length(z))
    n = length(t)
    n < 2 && return w
    dts = diff(t)
    mdt = median(dts)
    for i in 1:n-1
        isfinite(depth[i]) || continue
        kb = floor(Int, depth[i] / dz) + 1
        k = findfirst(zz -> abs(zz - (kb - 0.5) * dz) < dz / 4, z)
        k === nothing && continue
        w[k] += clamp(dts[i], 0, 10mdt)
    end
    return w
end

"""
    shear_segment(E, N, celldepth, offsets; opts=ShearOptions())

Bin-averaged vertical shear profile for one segment. Returns
`(z, sh_u, sh_v, nobs)` on bin centers `z` (NaN shear where `nobs < min_bin_obs`),
or `nothing` if no usable samples.
"""
function shear_segment(E::AbstractMatrix, N::AbstractMatrix, celldepth::AbstractMatrix,
                       offsets::AbstractVector; opts::ShearOptions=ShearOptions())
    ngrid, nt = size(E)
    dz = opts.dz
    zs = Float64[]; su = Float64[]; sv = Float64[]
    for i in 1:nt, k in 1:ngrid-1
        doff = offsets[k+1] - offsets[k]
        doff == 0 && continue
        (isfinite(E[k, i]) && isfinite(E[k+1, i]) &&
         isfinite(N[k, i]) && isfinite(N[k+1, i])) || continue
        zmid = (celldepth[k, i] + celldepth[k+1, i]) / 2
        (isfinite(zmid) && zmid >= 0) || continue
        push!(zs, zmid)
        push!(su, (E[k+1, i] - E[k, i]) / doff)   # ∂u/∂z, z positive down
        push!(sv, (N[k+1, i] - N[k, i]) / doff)
    end
    isempty(zs) && return nothing
    maxk = floor(Int, maximum(zs) / dz) + 1
    acc_u = [Float64[] for _ in 1:maxk]
    acc_v = [Float64[] for _ in 1:maxk]
    for m in eachindex(zs)
        kb = floor(Int, zs[m] / dz) + 1
        push!(acc_u[kb], su[m])
        push!(acc_v[kb], sv[m])
    end
    stat = opts.stat === :median ? median : mean
    z = [(k - 0.5) * dz for k in 1:maxk]
    sh_u = fill(NaN, maxk); sh_v = fill(NaN, maxk)
    nobs = length.(acc_u)
    for k in 1:maxk
        nobs[k] >= opts.min_bin_obs || continue
        sh_u[k] = stat(acc_u[k])
        sh_v[k] = stat(acc_v[k])
    end
    return (z=z, sh_u=sh_u, sh_v=sh_v, nobs=nobs)
end

"""
    integrate_shear(z, sh, dz) -> Vector

Top-down integration of a binned shear profile to velocity (arbitrary constant):
gaps (NaN shear) contribute zero (velocity held constant across them).
"""
function integrate_shear(z::AbstractVector, sh::AbstractVector, dz::Real)
    v = zeros(length(sh))
    acc = 0.0
    for k in eachindex(sh)
        isfinite(sh[k]) && (acc += sh[k] * dz)
        v[k] = acc
    end
    return v
end

"""
    solve_shear(pings::ProcessedPings, dac::DataFrame; opts=ShearOptions()) -> DataFrame

Shear-method absolute velocity per DAC segment: bin shear → integrate → reference the
mean over glider-sampled bins to the segment DAC. Output schema matches
[`solve_inverse`](@ref): `yo, t_mid, z, u, v, nobs`. Bins without shear samples carry
bridged (integrated-through) values flagged by `nobs = 0`.
"""
function solve_shear(p::ProcessedPings, dac::DataFrame; opts::ShearOptions=ShearOptions())
    out = DataFrame(yo=Int[], t_mid=DateTime[], z=Float64[], u=Float64[], v=Float64[],
        nobs=Int[])
    for row in eachrow(dac)
        idx = segment_indices(p, row.t_start, row.t_end)
        length(idx) >= opts.min_pings || continue
        seg = shear_segment(view(p.E, :, idx), view(p.N, :, idx),
            view(p.celldepth, :, idx), p.offsets; opts)
        seg === nothing && continue
        u_bc = integrate_shear(seg.z, seg.sh_u, opts.dz)
        v_bc = integrate_shear(seg.z, seg.sh_v, opts.dz)
        gd = filter(isfinite, p.depth[idx])
        isempty(gd) && continue
        covered = [k for k in eachindex(seg.z)
                   if seg.nobs[k] >= opts.min_bin_obs && seg.z[k] <= maximum(gd)]
        isempty(covered) && continue
        if opts.referencing === :timeweighted
            # DAC is a TIME average over the yo — weight bins by glider residence time
            w = time_in_bin(p.depth[idx], p.t[idx], seg.z; dz=opts.dz)[covered]
            sum(w) > 0 || (w = ones(length(covered)))
            uref = row.u - sum(w .* u_bc[covered]) / sum(w)
            vref = row.v - sum(w .* v_bc[covered]) / sum(w)
        else
            uref = row.u - mean(u_bc[covered])
            vref = row.v - mean(v_bc[covered])
        end
        for k in eachindex(seg.z)
            push!(out, (row.yo, row.t_mid, seg.z[k], u_bc[k] + uref, v_bc[k] + vref,
                seg.nobs[k]))
        end
    end
    return out
end
