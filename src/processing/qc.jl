# Layer 2 — quality control (masks beam velocities to NaN; never destroys amp/corr).
#
# Defaults follow Gradone et al. 2023 / gliderad2cp / Todd et al. 2017 (PLAN.md §6).

"""
    QCThresholds(; correlation=50, amplitude_max=75, snr_db=3, velocity_max=0.8,
                 ambiguity_frac=0.9, surface_depth=5, first_cells=0,
                 drop_error_pings=true)

Thresholds for [`qc!`](@ref). Defaults follow Gradone et al. (2023), gliderad2cp and
Todd et al. (2017); set a numeric field to `NaN` (or `first_cells = 0`) to disable that
screen. See the tutorial for tuning guidance.

`first_cells = 0` (keep the first cell) is validated for large-blanking
configurations: with ≥ 0.5 m blanking, cell 1 shows full correlation, on-curve
amplitude and unbiased velocities on four reference missions (~1.5× per-sample noise;
see the QA/QC guide) — dropping it discards good data. Deployments with small
blanking (Nortek default ~0.1 m) should set `first_cells = 1`: [`qc!`](@ref) warns
when the configured blanking is below 0.5 m and the first cell is being kept.
"""
Base.@kwdef struct QCThresholds
    correlation::Float64 = 50.0      # min %, Gradone 2023 (80 = gliderad2cp strict)
    amplitude_max::Float64 = 75.0    # max dB, Gradone 2023
    snr_db::Float64 = 3.0            # min dB above pooled noise floor (gliderad2cp); NaN disables
    velocity_max::Float64 = 0.8      # max |along-beam| m/s (gliderad2cp)
    ambiguity_frac::Float64 = 0.9    # flag |v| > frac × configured velocity range; NaN disables
    surface_depth::Float64 = 5.0     # mask pings shallower than this (m); NaN disables
    first_cells::Int = 0             # cells nearest the transducer to drop (see docstring)
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
    elseif isfinite(a.config.blanking) && a.config.blanking < 0.5
        @warn "qc!: keeping the first cell with blanking = $(a.config.blanking) m — " *
              "cell 1 is likely ringing-contaminated at small blanking; consider " *
              "QCThresholds(first_cells=1) (the 0-default is validated for ≥ 0.5 m)"
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

Bottom-track validity. Per-beam screens: figure of merit ≠ 65535 (Nortek invalid
marker), finite velocity, and detection distance within `[min_range, max_range]`.
Per-record screen (`bathymetry_check=true`): the implied bottom depth
(instrument depth + vertical range) must not be contradicted by the platform itself
diving deeper than it within ±`window` seconds (+`margin` m).

Both extra screens target **false bottom locks on near-field/water-borne targets**
(wake, scattering layers): on the reference mission, 99.7 % of "locks" were a
persistent target 0.6–2.8 m below the transducer, moving with the water — anchoring
the inverse to such targets contradicts the (earth-frame) DAC and injects strong
spurious shear into the solution (see docs/research/m38_validation.md, Task 3).
Genuine seafloor approaches have ranges of order 10 m or more, so `min_range = 5`
rejects the near-field cluster while keeping real locks.
"""
function bt_valid(bt::BottomTrackData; max_range::Real=Inf, min_range::Real=5.0,
                  bathymetry_check::Bool=true, window::Real=7200.0, margin::Real=20.0)
    n = length(bt)
    m = falses(4, n)
    for i in 1:n, b in 1:4
        m[b, i] = isfinite(bt.vel[b, i]) && bt.fom[b, i] != 65535 &&
                  isfinite(bt.distance[b, i]) &&
                  min_range <= bt.distance[b, i] <= max_range
    end
    if bathymetry_check
        for i in 1:n
            any(@view m[:, i]) || continue
            isfinite(bt.pressure[i]) || continue
            rng = [bt.distance[b, i] for b in 1:4 if m[b, i]]
            bottom = bt.pressure[i] + mean(rng) * cosd(30)   # ≈ vertical range in flight
            lo = searchsortedfirst(bt.t, bt.t[i] - window)
            hi = searchsortedlast(bt.t, bt.t[i] + window)
            deepest = -Inf
            for j in lo:hi
                isfinite(bt.pressure[j]) && bt.pressure[j] > deepest &&
                    (deepest = bt.pressure[j])
            end
            deepest > bottom + margin && (m[:, i] .= false)
        end
    end
    return m
end

"""
    cell_quality(adcp; thr=QCThresholds()) -> DataFrame

Per-cell, per-beam data-quality summary: median correlation and amplitude, and the
fraction of samples that would survive each QC screen (`keep_corr`, `keep_amp`,
`keep_all` — correlation, amplitude window, and their combination with the velocity
cap). Run on **unmasked** data (before `qc!`) to see how ping quality varies with range
and which cells actually contribute; after `qc!` the same call reports the surviving
fractions. The noise floor is the pooled 0.5th amplitude percentile, as in [`qc!`](@ref).
"""
function cell_quality(a::AD2CPData; thr::QCThresholds=QCThresholds())
    nc, nb, nt = size(a.vel)
    floor_db = isfinite(thr.snr_db) ? nanpctile(a.amp, 0.5) : -Inf
    rows = NamedTuple[]
    for b in 1:nb, k in 1:nc
        c = @view a.corr[k, b, :]
        m = @view a.amp[k, b, :]
        v = @view a.vel[k, b, :]
        fin = findall(i -> isfinite(c[i]) && isfinite(m[i]) && isfinite(v[i]), 1:nt)
        isempty(fin) && continue
        n = length(fin)
        kc = count(i -> c[i] >= thr.correlation, fin) / n
        ka = count(i -> m[i] <= thr.amplitude_max &&
                        (!isfinite(thr.snr_db) || m[i] >= floor_db + thr.snr_db), fin) / n
        kall = count(i -> c[i] >= thr.correlation &&
                          m[i] <= thr.amplitude_max &&
                          (!isfinite(thr.snr_db) || m[i] >= floor_db + thr.snr_db) &&
                          abs(v[i]) <= thr.velocity_max, fin) / n
        push!(rows, (cell=k, range=a.range[k], beam=b, n=n,
            med_corr=median(c[fin]), med_amp=median(m[fin]),
            keep_corr=kc, keep_amp=ka, keep_all=kall))
    end
    return DataFrame(rows)
end
