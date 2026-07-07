# Layer 0 — core types.
#
# Conventions (PLAN.md §4.2): depth positive down; velocities m/s; angles in degrees at
# API boundaries; missing data = NaN inside numeric arrays; time = DateTime (UTC) plus a
# derived Float64 Unix-seconds vector `t` for interpolation.

"""
    AD2CPConfig

Instrument configuration extracted from the netCDF `Config` group (or, later, from the
`.ad2cp` binary string record / `\$PNORI`). `beam2xyz` is the factory beam→XYZ matrix with
rows X, Y, Z1, Z2 and columns beams 1–4. `attrs` retains every raw Config attribute for
provenance.
"""
struct AD2CPConfig
    serial::Int
    frequency::Float64              # kHz
    beam_theta::NTuple{4,Float64}   # deg from instrument Z (1 fwd, 2 stbd, 3 aft, 4 port)
    beam_phi::NTuple{4,Float64}     # deg azimuth
    beam2xyz::Matrix{Float64}       # 4×4, rows X,Y,Z1,Z2
    ncells::Int
    cellsize::Float64               # m
    blanking::Float64               # m
    coordsystem::Symbol             # :beam | :xyz | :enu
    velocity_range::Float64         # ambiguity velocity (m/s)
    pressure_offset::Float64        # dbar (applied onboard)
    declination::Float64            # deg, as configured onboard (often 0 = unset)
    salinity_setting::Float64       # PSU used onboard for sound speed
    attrs::Dict{String,Any}
end

"""
    BottomTrackData

Bottom-track records (`Data/AverageBT`): per-beam over-ground velocity (beam coordinates,
m/s), detected range along beam (m) and figure of merit (`fom == 65535` ⇒ invalid), plus
the attitude/pressure needed to rotate BT velocities. Arrays are `4 × ntime`.
"""
struct BottomTrackData
    time::Vector{DateTime}
    t::Vector{Float64}              # unix seconds
    vel::Matrix{Float32}            # 4 × ntime
    distance::Matrix{Float64}       # 4 × ntime
    fom::Matrix{Float32}            # 4 × ntime
    pressure::Vector{Float64}
    heading::Vector{Float32}
    pitch::Vector{Float32}
    roll::Vector{Float32}
    soundspeed::Vector{Float32}
end

Base.length(bt::BottomTrackData) = length(bt.time)

"""
    AD2CPData

Profiling data from one plan (`Data/Average` or `Data/Burst`) of one or more
MIDAS-exported files, concatenated and time-sorted. Beam arrays are
`ncells × 4 beams × ntime`; `range` is the nominal cell-center range reported by the
instrument ("Velocity Range", m — see PLAN.md §2.2 on the 25° range-gating convention).
Velocities are in the coordinate system given by `config.coordsystem` (BEAM for the
standard glider configuration) and are **uncorrected** (sound speed as recorded).
"""
struct AD2CPData
    time::Vector{DateTime}
    t::Vector{Float64}              # unix seconds
    range::Vector{Float64}          # ncells
    vel::Array{Float32,3}           # ncells × 4 × ntime
    amp::Array{Float32,3}           # dB
    corr::Array{Float32,3}          # %
    heading::Vector{Float32}        # deg
    pitch::Vector{Float32}          # deg
    roll::Vector{Float32}           # deg
    pressure::Vector{Float64}       # dbar
    temperature::Vector{Float32}    # °C (instrument)
    soundspeed::Vector{Float32}     # m/s (as used onboard)
    accel::Matrix{Float32}          # 3 × ntime (g)
    mag::Matrix{Float32}            # 3 × ntime
    error::Vector{Float64}
    status::Vector{Float64}
    ensemble::Vector{Float64}
    config::AD2CPConfig
    bt::Union{BottomTrackData,Nothing}
end

Base.length(a::AD2CPData) = length(a.time)
ncells(a::AD2CPData) = length(a.range)

function Base.show(io::IO, a::AD2CPData)
    span = isempty(a.time) ? "empty" : "$(first(a.time)) → $(last(a.time))"
    print(io, "AD2CPData: $(length(a)) ensembles × $(ncells(a)) cells, ",
        "$(round(Int, a.config.frequency)) kHz SN$(a.config.serial), ",
        "coord=$(a.config.coordsystem), BT=$(isnothing(a.bt) ? "no" : "$(length(a.bt)) recs"), ",
        span)
end

"""
    GliderNav

Glider navigation time series (SeaExplorer `.gli` files or equivalent). Positions are in
decimal degrees; `deadreckoning` is 1 while subsurface (dead-reckoned position) and 0 when
GPS-fixed at the surface; `navstate` retains the platform state code. The full parsed
table is kept in `df` for anything not promoted to a typed field.
"""
struct GliderNav
    time::Vector{DateTime}
    t::Vector{Float64}
    lon::Vector{Float64}
    lat::Vector{Float64}
    heading::Vector{Float64}        # deg
    declination::Vector{Float64}    # deg (as logged; often 0)
    pitch::Vector{Float64}          # deg
    roll::Vector{Float64}           # deg
    depth::Vector{Float64}          # m, positive down
    navstate::Vector{Int16}
    deadreckoning::Vector{Int8}     # 1 = DR, 0 = GPS fix, -1 = unknown
    altitude::Vector{Float64}       # m above bottom (-1 = no lock)
    df::DataFrame
end

Base.length(n::GliderNav) = length(n.time)

function Base.show(io::IO, n::GliderNav)
    span = isempty(n.time) ? "empty" : "$(first(n.time)) → $(last(n.time))"
    nfix = count(==(0), n.deadreckoning)
    print(io, "GliderNav: $(length(n)) records, $nfix GPS-fixed, ", span)
end

"""
    a[idx] -> AD2CPData

Time-subset of the dataset by ping indices (bottom track restricted to the same time
window), e.g. `adcp[1:5000]` or `adcp[findall(p -> p > 100, adcp.pressure)]`.
"""
function Base.getindex(a::AD2CPData, idx::AbstractVector)
    time = a.time[idx]
    bt = a.bt
    if bt !== nothing && !isempty(time)
        t1, t2 = extrema(datetime2unix.(time))
        keep = findall(t -> t1 - 5 <= t <= t2 + 5, bt.t)
        bt = BottomTrackData(bt.time[keep], bt.t[keep], bt.vel[:, keep],
            bt.distance[:, keep], bt.fom[:, keep], bt.pressure[keep], bt.heading[keep],
            bt.pitch[keep], bt.roll[keep], bt.soundspeed[keep])
    end
    AD2CPData(time, a.t[idx], a.range, a.vel[:, :, idx], a.amp[:, :, idx],
        a.corr[:, :, idx], a.heading[idx], a.pitch[idx], a.roll[idx], a.pressure[idx],
        a.temperature[idx], a.soundspeed[idx], a.accel[:, idx], a.mag[:, idx],
        a.error[idx], a.status[idx], a.ensemble[idx], a.config, bt)
end
