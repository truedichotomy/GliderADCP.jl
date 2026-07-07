# Layer 1 — Nortek $PNOR ASCII stream parser (SeaExplorer real-time payload feed).
#
# The AD2CP streams NMEA-style sentences to the payload computer in real time
# (logged as `sea<glider>.<mission>.ad2cp.raw.<N>[.gz]`):
#
#   $PNORI,itype,headid,nbeams,ncells,blanking,cellsize,coord*cs      (configuration)
#   $PNORS,MMDDYY,HHMMSS,errhex,statushex,battery,soundspeed,
#          heading,pitch,roll,pressure,temperature,a1,a2*cs           (per-ensemble)
#   $PNORC,MMDDYY,HHMMSS,cell,v1,v2,v3,v4,speed,dir,unit,
#          amp1..amp4,corr1..corr4*cs                                 (per cell)
#
# Amplitudes arrive in counts (≈0.5 dB/count — converted to dB here for consistency
# with the MIDAS netCDF variables); velocities in m/s at 0.01 resolution. This enables
# real-time/piloting processing and cross-checks of the delayed-mode pipeline.

function _nmea_checksum_ok(line::AbstractString)
    startswith(line, '$') || return false
    star = findlast('*', line)
    (star === nothing || star + 2 > lastindex(line)) && return false
    cs = 0x00
    for c in codeunits(SubString(line, 2, star - 1))
        cs ⊻= c
    end
    parsed = tryparse(UInt8, line[star+1:star+2]; base=16)
    return parsed !== nothing && parsed == cs
end

_pnor_fields(line) = split(SubString(line, 1, something(findlast('*', line), lastindex(line) + 1) - 1), ',')
_pnor_f(s) = isempty(s) ? NaN : something(tryparse(Float64, s), NaN)

function _pnor_datetime(datestr, timestr)
    (length(datestr) == 6 && length(timestr) == 6) || return nothing
    mm = tryparse(Int, datestr[1:2]); dd = tryparse(Int, datestr[3:4])
    yy = tryparse(Int, datestr[5:6])
    H = tryparse(Int, timestr[1:2]); M = tryparse(Int, timestr[3:4])
    S = tryparse(Int, timestr[5:6])
    any(isnothing, (mm, dd, yy, H, M, S)) && return nothing
    try
        return DateTime(2000 + yy, mm, dd, H, M, S)
    catch
        return nothing
    end
end

"""
    load_pnor(src; stream="ad2cp.raw", validate_checksum=true,
              mintime=DateTime(2000)) -> AD2CPData

Parse the AD2CP real-time NMEA stream into an [`AD2CPData`](@ref) (same structure as
the MIDAS netCDF path, minus bottom track and accelerometer/magnetometer — pass
`look=:down`/`:up` explicitly to downstream functions). `src` is a directory (all
segments of `stream`) or a vector of files (gzip transparent).
"""
function load_pnor(src; stream::AbstractString="ad2cp.raw", validate_checksum::Bool=true,
                   mintime::Union{DateTime,Nothing}=DateTime(2000))
    files = src isa AbstractVector ? String.(src) : seaexplorer_files(src, stream)
    isempty(files) && error("load_pnor: no input files")

    nbeams = 4; ncells = 0; blanking = NaN; cellsize = NaN; coord = 2; serial = 0
    ens_time = DateTime[]
    sens = NamedTuple[]              # per-ensemble sensor tuple
    cells = Vector{Matrix{Float32}}[]  # per-ensemble [vel amp corr] (ncells × 4 each)
    cur = nothing                    # current ensemble datetime

    for f in files
        raw = endswith(f, ".gz") ? open(io -> read(GzipDecompressorStream(io)), f) : read(f)
        for line in eachline(IOBuffer(raw))
            line = strip(line)
            startswith(line, "\$PNOR") || continue
            validate_checksum && !_nmea_checksum_ok(line) && continue
            p = _pnor_fields(line)
            if p[1] == "\$PNORI" && length(p) >= 8
                nbeams = something(tryparse(Int, p[4]), 4)
                ncells = something(tryparse(Int, p[5]), ncells)
                blanking = _pnor_f(p[6])
                cellsize = _pnor_f(p[7])
                coord = something(tryparse(Int, p[8]), 2)
                m = match(r"(\d+)", p[3])
                m === nothing || (serial = parse(Int, m.captures[1]))
            elseif p[1] == "\$PNORS" && length(p) >= 12
                dt = _pnor_datetime(p[2], p[3])
                (dt === nothing || (mintime !== nothing && dt < mintime)) && (cur = nothing; continue)
                cur = dt
                push!(ens_time, dt)
                push!(sens, (
                    error = something(tryparse(UInt32, p[4]; base=16), UInt32(0)),
                    status = something(tryparse(UInt32, p[5]; base=16), UInt32(0)),
                    soundspeed = _pnor_f(p[7]), heading = _pnor_f(p[8]),
                    pitch = _pnor_f(p[9]), roll = _pnor_f(p[10]),
                    pressure = _pnor_f(p[11]), temperature = _pnor_f(p[12])))
                push!(cells, [fill(NaN32, max(ncells, 1), 4) for _ in 1:3])
            elseif p[1] == "\$PNORC" && length(p) >= 19 && cur !== nothing
                dt = _pnor_datetime(p[2], p[3])
                dt == cur || continue
                k = tryparse(Int, p[4])
                (k === nothing || k < 1) && continue
                V, A, C = cells[end]
                if k > size(V, 1)   # grow if PNORI undersold the cell count
                    for (m, M) in enumerate(cells[end])
                        cells[end][m] = vcat(M, fill(NaN32, k - size(M, 1), 4))
                    end
                    V, A, C = cells[end]
                end
                for b in 1:4
                    V[k, b] = Float32(_pnor_f(p[4+b]))
                    A[k, b] = Float32(_pnor_f(p[11+b]) * 0.5)   # counts → dB
                    C[k, b] = Float32(_pnor_f(p[15+b]))
                end
            end
        end
    end
    isempty(ens_time) && error("load_pnor: no ensembles parsed")

    nc = maximum(size(c[1], 1) for c in cells)
    nt = length(ens_time)
    vel = fill(NaN32, nc, 4, nt); amp = fill(NaN32, nc, 4, nt); corr = fill(NaN32, nc, 4, nt)
    for (i, c) in enumerate(cells)
        r = size(c[1], 1)
        vel[1:r, :, i] = c[1]; amp[1:r, :, i] = c[2]; corr[1:r, :, i] = c[3]
    end
    perm = sortperm(ens_time)
    time = ens_time[perm]
    g(f) = [Float32(getfield(s, f)) for s in sens][perm]

    cfg = AD2CPConfig(serial, 1000.0, (47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0),
        fill(NaN, 4, 4), nc, cellsize, blanking,
        coord == 2 ? :beam : coord == 1 ? :xyz : :enu,
        NaN, NaN, NaN, NaN, Dict{String,Any}("source" => "\$PNOR stream"))
    range = blanking .+ cellsize .* (1:nc)
    AD2CPData(time, datetime2unix.(time), collect(range),
        vel[:, :, perm], amp[:, :, perm], corr[:, :, perm],
        g(:heading), g(:pitch), g(:roll), Float64.(g(:pressure)),
        g(:temperature), g(:soundspeed),
        fill(NaN32, 3, nt), fill(NaN32, 3, nt),
        [Float64(s.error) for s in sens][perm], [Float64(s.status) for s in sens][perm],
        collect(1.0:nt), cfg, nothing)
end
