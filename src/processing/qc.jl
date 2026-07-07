# Layer 2 — quality control (masks beam velocities to NaN; never destroys amp/corr).
#
# Defaults follow Gradone et al. 2023 / gliderad2cp / Todd et al. 2017 (PLAN.md §6).

Base.@kwdef struct QCThresholds
    correlation::Float64 = 50.0      # min %, Gradone 2023 (80 = gliderad2cp strict)
    amplitude_max::Float64 = 75.0    # max dB, Gradone 2023
    snr_db::Float64 = 3.0            # min dB above pooled noise floor (gliderad2cp); NaN disables
    velocity_max::Float64 = 0.8      # max |along-beam| m/s (gliderad2cp)
    ambiguity_frac::Float64 = 0.9    # flag |v| > frac × configured velocity range; NaN disables
    surface_depth::Float64 = 5.0     # mask pings shallower than this (m); NaN disables
    first_cells::Int = 1             # cells nearest the transducer to drop
    drop_error_pings::Bool = true    # instrument Error ≠ 0 → drop whole ping
end

"""
    qc!(adcp; thr=QCThresholds(), depth=nothing) -> NamedTuple

Apply QC masks in place to `adcp.vel` (NaN-ing rejected samples) and return per-screen
rejection statistics (fractions of previously-finite samples). `depth` is an optional
per-ping glider depth vector (m, positive down) for the surface mask; when `nothing`,
`adcp.pressure` is used as a proxy (1 dbar ≈ 1 m).
"""
function qc!(a::AD2CPData; thr::QCThresholds=QCThresholds(), depth=nothing)
    v = a.vel
    n0 = count(isfinite, v)
    n0 == 0 && return (; total=0.0)
    z = depth === nothing ? a.pressure : depth

    reject(mask) = (k = 0;
        @inbounds for idx in eachindex(v)
            if mask[idx] && isfinite(v[idx])
                v[idx] = NaN32; k += 1
            end
        end; k / n0)

    stats = Dict{Symbol,Float64}()
    stats[:correlation] = reject(a.corr .< thr.correlation)
    stats[:amplitude_max] = reject(a.amp .> thr.amplitude_max)
    if isfinite(thr.snr_db)
        floor_db = nanpctile(a.amp, 0.5)
        stats[:snr] = reject(a.amp .< floor_db + thr.snr_db)
    end
    stats[:velocity_max] = reject(abs.(v) .> thr.velocity_max)
    if isfinite(thr.ambiguity_frac) && isfinite(a.config.velocity_range)
        stats[:ambiguity] = reject(abs.(v) .> thr.ambiguity_frac * a.config.velocity_range)
    end

    nc, nb, nt = size(v)
    if thr.first_cells > 0
        m = falses(nc, nb, nt); m[1:min(thr.first_cells, nc), :, :] .= true
        stats[:first_cells] = reject(m)
    end
    if isfinite(thr.surface_depth)
        m = falses(nc, nb, nt)
        for i in 1:nt
            isfinite(z[i]) && z[i] <= thr.surface_depth && (m[:, :, i] .= true)
        end
        stats[:surface] = reject(m)
    end
    if thr.drop_error_pings
        m = falses(nc, nb, nt)
        for i in 1:nt
            isfinite(a.error[i]) && a.error[i] != 0 && (m[:, :, i] .= true)
        end
        stats[:error_pings] = reject(m)
    end

    stats[:total] = (n0 - count(isfinite, v)) / n0
    return NamedTuple(stats)
end

"""
    bt_valid(bt; max_range=Inf) -> BitMatrix (4 × ntime)

Bottom-track validity: figure of merit ≠ 65535 (Nortek invalid marker), finite velocity,
and detection distance within `(0, max_range]`.
"""
function bt_valid(bt::BottomTrackData; max_range::Real=Inf)
    m = falses(4, length(bt))
    for i in 1:length(bt), b in 1:4
        m[b, i] = isfinite(bt.vel[b, i]) && bt.fom[b, i] != 65535 &&
                  isfinite(bt.distance[b, i]) && 0 < bt.distance[b, i] <= max_range
    end
    return m
end
