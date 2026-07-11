# Layer 1 — telemetered AD2CP pings from the SeaExplorer payload record (pld1.sub).
#
# The AD2CP driver broadcasts one instrument ensemble every ~30 s into the payload
# record; the subsampled `pld1.sub` files carrying it are what the glider actually
# transmits over Iridium mid-mission. Per broadcast row:
#
#   AD2CP_TIME                    instrument timestamp, "MMDDYY HH:MM:SS"
#   AD2CP_HEADING/PITCH/ROLL      deg (0.1 resolution)
#   AD2CP_PRESSURE                dbar
#   AD2CP_V<b>_CN<k>              beam velocity, beam b ∈ 1:4, cell k (m/s, 0.01
#                                 quantization; ±9999 while the instrument is off)
#
# Verified against the M38 binary: each row is a single subsampled ensemble
# (values = instrument cells 1..K rounded to 0.01 m/s; NOT a 30-s average), in beam
# coordinates, with attitude matching the same ensemble. No amplitude, correlation,
# or bottom track is transmitted — the corresponding QC screens become no-ops, and
# `first_cells` (on by default) removes the ringing-contaminated cell 1.
#
# This is the true shore-side real-time route (the `$PNOR` stream in `ad2cp.raw`
# is payload-logged only and comes back with the glider). ALSEAMAR's own GLIMPSE
# processing of these rows appears as `AD2CP_*_c` columns in server exports; this
# reader ingests the raw beam data so the product can be computed openly.

const _PLDADCP_TIME_RE = r"^\d{6} \d{2}:\d{2}:\d{2}$"

# AD2CP_TIME is MMDDYY HH:MM:SS (Nortek convention, as in $PNORS)
function _pldadcp_time(s::AbstractString)
    occursin(_PLDADCP_TIME_RE, s) || return nothing
    mm = parse(Int, s[1:2]); dd = parse(Int, s[3:4]); yy = parse(Int, s[5:6])
    try
        return DateTime(2000 + yy, mm, dd,
            parse(Int, s[8:9]), parse(Int, s[11:12]), parse(Int, s[14:15]))
    catch
        return nothing
    end
end

_pldadcp_f(s) = begin
    v = isempty(s) ? NaN : something(tryparse(Float64, s), NaN)
    (v == 9999.0 || v == -9999.0) ? NaN : v
end

_pldadcp_open(path) = endswith(path, ".gz") ?
    IOBuffer(open(io -> read(CodecZlib.GzipDecompressorStream(io)), path)) : open(path, "r")

# parse one pld1.sub file (segment log or GLIMPSE csv) into per-ping tuples
function _read_pldadcp_file!(out::Dict{DateTime,Any}, path::AbstractString)
    io = _pldadcp_open(path)
    try
        hdr = String.(split(strip(readline(io), [';', ' ', '\r']), ';'))
        col = Dict(n => i for (i, n) in enumerate(hdr))
        haskey(col, "AD2CP_TIME") || return 0     # stream without the AD2CP subset
        it = col["AD2CP_TIME"]
        ih = get(col, "AD2CP_HEADING", 0); ip = get(col, "AD2CP_PITCH", 0)
        ir = get(col, "AD2CP_ROLL", 0); ipr = get(col, "AD2CP_PRESSURE", 0)
        # velocity columns present in this file, as (cell, beam) => column
        vcols = Tuple{Int,Int,Int}[]
        for (n, i) in col
            m = match(r"^AD2CP_V([1-4])_CN(\d+)$", n)
            m === nothing || push!(vcols, (parse(Int, m.captures[2]), parse(Int, m.captures[1]), i))
        end
        isempty(vcols) && return 0
        nc = maximum(first.(vcols))
        nread = 0
        for ln in eachline(io)
            isempty(ln) && continue
            p = split(rstrip(ln, '\r'), ';')
            length(p) < it && continue
            tt = _pldadcp_time(p[it])
            tt === nothing && continue
            haskey(out, tt) && continue           # multi-route duplicate
            h = ih > 0 && ih <= length(p) ? _pldadcp_f(p[ih]) : NaN
            isfinite(h) || continue               # instrument off / fill row
            V = fill(NaN32, nc, 4)
            for (k, b, i) in vcols
                i <= length(p) && (V[k, b] = Float32(_pldadcp_f(p[i])))
            end
            all(isnan, V) && continue
            out[tt] = (h,
                ip > 0 && ip <= length(p) ? _pldadcp_f(p[ip]) : NaN,
                ir > 0 && ir <= length(p) ? _pldadcp_f(p[ir]) : NaN,
                ipr > 0 && ipr <= length(p) ? _pldadcp_f(p[ipr]) : NaN,
                V)
            nread += 1
        end
        return nread
    finally
        close(io)
    end
end

"""
    load_pld_adcp(src; stream="pld1.sub", cellsize, blanking, serial=0,
                  soundspeed=NaN, mintime=DateTime(2000)) -> AD2CPData

Read the **telemetered** AD2CP pings from SeaExplorer payload records — the
`AD2CP_V<beam>_CN<cell>` beam velocities plus attitude/pressure that the glider
transmits over Iridium inside `pld1.sub` (one subsampled instrument ensemble every
~30 s; verified single-ensemble, beam coordinates, 0.01 m/s quantization). This is
the data available *shore-side during the mission*, unlike the `\$PNOR` stream
(`load_pnor`), which is payload-logged only and recovered with the glider.

`src` is a directory, a vector of directories (glider-computer segment logs +
GLIMPSE-server exports; duplicate instrument timestamps are deduplicated, earlier
sources first), a single file, or a file vector. Segment-numbered `.gz` logs,
GLIMPSE `.all.csv` and per-cycle `.NNN.csv` exports are all recognized.

`cellsize` and `blanking` (m) are **required** — they are not transmitted, and come
from the deployment configuration. `soundspeed` (scalar or per-ping vector, m/s)
records the sound speed the instrument used onboard so that
[`soundspeed_correction`](@ref) can rescale; it is likewise not transmitted
(reconstruct it from the configured salinity and the payload CTD temperature, or
leave `NaN` to skip the correction). No amplitude/correlation/bottom-track data
exist on this route; the corresponding QC screens pass everything, `first_cells`
still removes the ringing cell, and `process_pings` needs an explicit `look=`.

```julia
tele = load_pld_adcp(["delayed/pld1/logs", "glimpse"]; stream="38.pld1.sub",
                     cellsize=2.0, blanking=0.7, serial=102381)
```
"""
function load_pld_adcp(src; stream::AbstractString="pld1.sub",
                       cellsize::Real, blanking::Real, serial::Integer=0,
                       soundspeed=NaN,
                       mintime::Union{DateTime,Nothing}=DateTime(2000))
    dirs = src isa AbstractString ? [String(src)] : String.(src)
    files = String[]
    if all(isdir, dirs)
        for d in dirs
            append!(files, seaexplorer_files(d, stream))
            append!(files, SeaExplorerIO.glimpse_files(d, stream))
            # GLIMPSE per-cycle exports: <glider>.<mission>.<stream>.NNN.csv
            pat = Regex("\\." * replace(stream, "." => "\\.") * "\\.(\\d+)\\.csv\$", "i")
            append!(files, sort([joinpath(d, f) for f in readdir(d) if occursin(pat, f)]))
        end
    else
        files = dirs                              # explicit file list
    end
    isempty(files) && error("load_pld_adcp: no files matching stream `$stream`")

    out = Dict{DateTime,Any}()
    nbad = 0
    for f in files
        try
            _read_pldadcp_file!(out, f)
        catch err
            nbad += 1
            @warn "load_pld_adcp: skipping unreadable file" file = basename(f) error = sprint(showerror, err)
        end
    end
    nbad == length(files) && error("load_pld_adcp: all $(length(files)) files unreadable")
    mintime !== nothing && filter!(kv -> kv.first >= mintime, out)
    isempty(out) && error("load_pld_adcp: no telemetered AD2CP pings found (instrument off, or wrong stream?)")

    time = sort!(collect(keys(out)))
    nt = length(time)
    nc = maximum(size(out[t][5], 1) for t in time)
    vel = fill(NaN32, nc, 4, nt)
    heading = fill(NaN32, nt); pitch = fill(NaN32, nt); roll = fill(NaN32, nt)
    pressure = fill(NaN, nt)
    for (i, tt) in enumerate(time)
        h, p, r, pr, V = out[tt]
        heading[i] = h; pitch[i] = p; roll[i] = r; pressure[i] = pr
        vel[1:size(V, 1), :, i] = V
    end
    ss = soundspeed isa Real ? fill(Float32(soundspeed), nt) : Float32.(soundspeed)
    length(ss) == nt || error("load_pld_adcp: soundspeed vector length $(length(ss)) ≠ $nt pings")

    cfg = AD2CPConfig(Int(serial), 1000.0, (47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0),
        fill(NaN, 4, 4), nc, Float64(cellsize), Float64(blanking), :beam,
        NaN, NaN, NaN, NaN,
        Dict{String,Any}("source" => "telemetered pld1.sub AD2CP subset"))
    range = Float64(blanking) .+ Float64(cellsize) .* (1:nc)
    AD2CPData(time, datetime2unix.(time), collect(range),
        vel, fill(NaN32, nc, 4, nt), fill(NaN32, nc, 4, nt),
        heading, pitch, roll, pressure,
        fill(NaN32, nt), ss,
        fill(NaN32, 3, nt), fill(NaN32, 3, nt),
        zeros(nt), zeros(nt),
        collect(1.0:nt), cfg, nothing)
end
