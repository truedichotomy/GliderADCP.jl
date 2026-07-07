# Layer 1 (native) — .ad2cp binary reader, removing the MIDAS/Windows dependency.
#
# Format: Nortek "Integrators Guide — AD2CP" (N3015-007), §Data formats; cross-checked
# against DOLfYN and validated bit-for-bit against the MIDAS netCDF export of the same
# file (see tests). Full spec notes: docs/research/formats_and_ecosystem.md.
#
#   File = sequence of (10-byte header, data record), all little-endian.
#   Header: 0xA5, headerSize=10, ID, family=0x10, dataSize u16, dataChecksum u16,
#           headerChecksum u16 (over header bytes 0–7).
#   Checksum: 16-bit wraparound sum of u16 LE words, init 0xB58C; odd trailing byte
#             added as (byte << 8).
#   IDs: 0x15 burst, 0x16 average, 0x17 bottom track, 0xA0 string record
#        (string ID 0x10 = the full instrument configuration dump).
#   DF3 (burst/average) records are self-describing: a config bitmask says which data
#   blocks follow the 76-byte fixed prefix; beams/coordinate-system/cells are packed in
#   one u16; velocities scale as 10^velocityScaling.

# ---- little-endian byte readers (1-based `o` relative to a record's first byte) -----
@inline _u8(b, o) = b[o]
@inline _u16(b, o) = UInt16(b[o]) | (UInt16(b[o+1]) << 8)
@inline _i16(b, o) = reinterpret(Int16, _u16(b, o))
@inline _u32(b, o) = UInt32(_u16(b, o)) | (UInt32(_u16(b, o + 2)) << 16)
@inline _i32(b, o) = reinterpret(Int32, _u32(b, o))
@inline _i8(b, o) = reinterpret(Int8, b[o])

function _ad2cp_checksum(b::AbstractVector{UInt8})
    cs::UInt16 = 0xB58C
    n = length(b)
    @inbounds for i in 1:2:n-1
        cs += UInt16(b[i]) | (UInt16(b[i+1]) << 8)
    end
    isodd(n) && (cs += UInt16(b[n]) << 8)
    return cs
end

"""
    _scan_ad2cp(buf; validate=true) -> (records, nresync)

Scan a raw .ad2cp buffer for records: returns `(id, start, size)` (start = 1-based index
of the first data byte) for every record whose header (and, when `validate`, data)
checksum verifies; resynchronizes on corruption by searching for the next valid header.
"""
function _scan_ad2cp(buf::Vector{UInt8}; validate::Bool=true)
    recs = @NamedTuple{id::UInt8, start::Int, size::Int}[]
    pos = 1
    n = length(buf)
    nresync = 0
    while pos + 9 <= n
        if !(buf[pos] == 0xA5 && buf[pos+1] == 0x0A)
            pos += 1
            nresync += 1
            continue
        end
        if _ad2cp_checksum(view(buf, pos:pos+7)) != _u16(buf, pos + 8)
            pos += 1
            nresync += 1
            continue
        end
        dsz = Int(_u16(buf, pos + 4))
        if pos + 9 + dsz > n
            break
        end
        ok = !validate ||
             _ad2cp_checksum(view(buf, pos+10:pos+9+dsz)) == _u16(buf, pos + 6)
        ok && push!(recs, (id=buf[pos+2], start=pos + 10, size=dsz))
        pos += 10 + dsz
    end
    return recs, nresync
end

_bit(x, k) = (x >> k) & 0x1 == 1

# clock fields at record offset 8 (0-based): year-1900, month 0-based, day, h, m, s + 100µs
function _df_datetime(d, o1)
    y = Int(_u8(d, o1)) + 1900
    mo = Int(_u8(d, o1 + 1)) + 1
    dy = Int(_u8(d, o1 + 2))
    (1 <= mo <= 12 && 1 <= dy <= 31) || return nothing
    try
        return DateTime(y, mo, dy, Int(_u8(d, o1 + 3)), Int(_u8(d, o1 + 4)),
            Int(_u8(d, o1 + 5))) + Millisecond(round(Int, Int(_u16(d, o1 + 6)) / 10))
    catch
        return nothing
    end
end

# ---- instrument configuration from the 0xA0/0x10 string record ----------------------
function _parse_config_string(txt::AbstractString)
    kv(line) = Dict(m.captures[1] => m.captures[2]
                    for m in eachmatch(r"(\w+)=(\"[^\"]*\"|[^,\s]+)", line))
    unq(s) = strip(s, '"')
    getl(cmd) = begin
        m = match(Regex("(?:^|\\n)(?:GET)?" * cmd * ",([^\\r\\n]*)"), txt)
        m === nothing ? nothing : kv(m.captures[1])
    end
    fnum(d, k, def) = d === nothing || !haskey(d, k) ? def :
                      something(tryparse(Float64, unq(d[k])), def)

    avg = getl("AVG")
    plan = getl("PLAN")
    user = getl("USER")
    idl = getl("ID")

    θ = fill(NaN, 4); φ = fill(NaN, 4)
    for m in eachmatch(r"BEAMCFGLIST,BEAM=(\d+),THETA=([-\d.]+),PHI=([-\d.]+)", txt)
        b = parse(Int, m.captures[1])
        1 <= b <= 4 || continue
        θ[b] = parse(Float64, m.captures[2])
        φ[b] = parse(Float64, m.captures[3])
    end
    b2x = fill(NaN, 4, 4)
    mx = match(r"GETXFAVG,ROWS=4,COLS=4,([^\r\n]*)", txt)
    if mx !== nothing
        for m in eachmatch(r"M(\d)(\d)=([-\d.]+)", mx.captures[1])
            b2x[parse(Int, m.captures[1]), parse(Int, m.captures[2])] =
                parse(Float64, m.captures[3])
        end
    end
    cy = avg === nothing || !haskey(avg, "CY") ? "BEAM" : unq(avg["CY"])
    return (
        serial = idl === nothing ? 0 : round(Int, fnum(idl, "SN", 0.0)),
        frequency = fnum(plan, "FREQ", NaN),
        θ = (θ[1], θ[2], θ[3], θ[4]), φ = (φ[1], φ[2], φ[3], φ[4]),
        beam2xyz = b2x,
        ncells = round(Int, fnum(avg, "NC", 0.0)),
        cellsize = fnum(avg, "CS", NaN),
        blanking = fnum(avg, "BD", NaN),
        coordsystem = cy == "BEAM" ? :beam : cy == "XYZ" ? :xyz : :enu,
        velocity_range = fnum(avg, "VR", NaN),
        pressure_offset = fnum(user, "POFF", NaN),
        declination = fnum(user, "DECL", NaN),
        salinity = fnum(plan, "SA", NaN),
    )
end

"""
    read_ad2cp(path; plan=:average, validate_checksums=true) -> AD2CPData

Native reader for Nortek `.ad2cp` binary files — same output structure as
[`load_ad2cp`](@ref) (MIDAS netCDF path), with no MIDAS dependency. Reads the selected
plan's profiling records (`:average` → ID 0x16, `:burst` → 0x15), bottom-track records
(0x17) when present, and the instrument configuration from the embedded string record.

Velocities/amplitudes/correlations and all per-ping sensors reproduce the MIDAS export
bit-for-bit (validated on the reference mission). Records failing their checksums are
skipped (counted in a `@warn`).
"""
function read_ad2cp(path::AbstractString; plan::Symbol=:average,
                    validate_checksums::Bool=true)
    buf = read(path)
    recs, nresync = _scan_ad2cp(buf; validate=validate_checksums)
    isempty(recs) && error("read_ad2cp: no valid records in $path")
    nresync > 0 && @warn "read_ad2cp: resynchronized past $nresync corrupt bytes"

    want = plan === :burst ? 0x15 : 0x16
    prof = [r for r in recs if r.id == want]
    btrecs = [r for r in recs if r.id == 0x17]
    isempty(prof) &&
        error("read_ad2cp: no $(plan) records (IDs present: $(unique(getfield.(recs, :id))))")

    # instrument configuration string (string record, string-ID 0x10)
    cfgtxt = ""
    for r in recs
        r.id == 0xA0 || continue
        d = view(buf, r.start:r.start+r.size-1)
        d[1] == 0x10 || continue
        z = findfirst(==(0x00), d)
        cfgtxt = String(UInt8.(d[2:(z === nothing ? length(d) : z - 1)]))
        break
    end
    sc = _parse_config_string(cfgtxt)

    # ---- profiling records (DF3) ------------------------------------------------------
    # geometry must be uniform across the file (multi-plan files are out of scope here)
    d0 = view(buf, prof[1].start:prof[1].start+prof[1].size-1)
    Int(_u8(d0, 1)) == 3 || error("read_ad2cp: unsupported record version $(_u8(d0, 1))")
    serial_rec = Int(_u32(d0, 5))
    beams_cy = _u16(d0, 31)
    nc = Int(beams_cy & 0x03FF)
    csys = Int((beams_cy >> 10) & 0x3)
    nb = Int((beams_cy >> 12) & 0xF)
    nb == 4 || @warn "read_ad2cp: $nb beams (expected 4 for the glider AD2CP)"

    nt = length(prof)
    vel = fill(NaN32, nc, 4, nt)
    amp = fill(NaN32, nc, 4, nt)
    corr = fill(NaN32, nc, 4, nt)
    time = Vector{DateTime}(undef, nt)
    heading = fill(NaN32, nt); pitch = fill(NaN32, nt); roll = fill(NaN32, nt)
    pressure = fill(NaN, nt); temperature = fill(NaN32, nt); soundspeed = fill(NaN32, nt)
    accel = fill(NaN32, 3, nt); mag = fill(NaN32, 3, nt)
    errorv = fill(NaN, nt); status = fill(NaN, nt); ensemble = fill(NaN, nt)
    ambig = fill(NaN, nt); cellsize_rec = NaN; blank_rec = NaN
    keep = trues(nt)

    for (i, r) in enumerate(prof)
        d = view(buf, r.start:r.start+r.size-1)
        if _u16(d, 31) != beams_cy || Int(_u8(d, 1)) != 3
            keep[i] = false
            continue
        end
        dt = _df_datetime(d, 9)
        dt === nothing && (keep[i] = false; continue)
        time[i] = dt
        cfg = _u16(d, 3)
        st = _u32(d, 69)
        vscale = 10.0^Int(_i8(d, 59))
        soundspeed[i] = _u16(d, 17) * 0.1f0
        temperature[i] = _i16(d, 19) * 0.01f0
        pressure[i] = _u32(d, 21) * 0.001
        heading[i] = _u16(d, 25) * 0.01f0
        pitch[i] = _i16(d, 27) * 0.01f0
        roll[i] = _i16(d, 29) * 0.01f0
        cellsize_rec = _u16(d, 33) * 0.001              # mm → m
        blank_rec = _u16(d, 35) * (_bit(st, 1) ? 0.01 : 0.001)
        for a in 1:3
            mag[a, i] = _i16(d, 41 + 2(a - 1))
            accel[a, i] = _i16(d, 47 + 2(a - 1)) / 16384.0f0
        end
        ambig[i] = _u16(d, 53) * vscale
        errorv[i] = Float64(_u16(d, 65))
        status[i] = Float64(st)
        ensemble[i] = Float64(_u32(d, 73))

        o = Int(_u8(d, 2)) + 1                           # first data byte (1-based)
        if _bit(cfg, 5)                                  # velocity i16[nb][nc]
            for b in 1:nb, k in 1:nc
                vel[k, b, i] = Float32(_i16(d, o + 2 * ((b - 1) * nc + (k - 1))) * vscale)
            end
            o += 2 * nb * nc
        end
        if _bit(cfg, 6)                                  # amplitude u8, 0.5 dB/count
            for b in 1:nb, k in 1:nc
                amp[k, b, i] = _u8(d, o + (b - 1) * nc + (k - 1)) * 0.5f0
            end
            o += nb * nc
        end
        if _bit(cfg, 7)                                  # correlation u8, %
            for b in 1:nb, k in 1:nc
                corr[k, b, i] = Float32(_u8(d, o + (b - 1) * nc + (k - 1)))
            end
        end
    end
    kept = findall(keep)
    n_dropped = nt - length(kept)
    n_dropped > 0 && @warn "read_ad2cp: skipped $n_dropped inconsistent profiling records"

    # ---- bottom-track records (DF20) ---------------------------------------------------
    bt = nothing
    if !isempty(btrecs)
        nbt = length(btrecs)
        btvel = fill(NaN32, 4, nbt); btdist = fill(NaN, 4, nbt); btfom = fill(NaN32, 4, nbt)
        bttime = Vector{DateTime}(undef, nbt)
        btp = fill(NaN, nbt)
        bth = fill(NaN32, nbt); btpi = fill(NaN32, nbt); btr = fill(NaN32, nbt)
        btss = fill(NaN32, nbt)
        btkeep = trues(nbt)
        for (i, r) in enumerate(btrecs)
            d = view(buf, r.start:r.start+r.size-1)
            dt = _df_datetime(d, 9)
            dt === nothing && (btkeep[i] = false; continue)
            bttime[i] = dt
            cfg = _u16(d, 3)
            bnb = Int((_u16(d, 31) >> 12) & 0xF)
            vscale = 10.0^Int(_i8(d, 61))
            btss[i] = _u16(d, 17) * 0.1f0
            btp[i] = _u32(d, 21) * 0.001
            bth[i] = _u16(d, 25) * 0.01f0
            btpi[i] = _i16(d, 27) * 0.01f0
            btr[i] = _i16(d, 29) * 0.01f0
            o = Int(_u8(d, 2)) + 1
            if _bit(cfg, 5)
                for b in 1:min(bnb, 4)
                    btvel[b, i] = Float32(_i32(d, o + 4(b - 1)) * vscale)
                end
                o += 4bnb
            end
            if _bit(cfg, 8)
                for b in 1:min(bnb, 4)
                    btdist[b, i] = _i32(d, o + 4(b - 1)) * 0.001
                end
                o += 4bnb
            end
            if _bit(cfg, 9)
                for b in 1:min(bnb, 4)
                    btfom[b, i] = Float32(_u16(d, o + 2(b - 1)))
                end
            end
        end
        bk = findall(btkeep)
        bt = BottomTrackData(bttime[bk], datetime2unix.(bttime[bk]), btvel[:, bk],
            btdist[:, bk], btfom[:, bk], btp[bk], bth[bk], btpi[bk], btr[bk], btss[bk])
    end

    config = AD2CPConfig(
        sc.serial != 0 ? sc.serial : serial_rec,
        sc.frequency, sc.θ, sc.φ, sc.beam2xyz,
        sc.ncells != 0 ? sc.ncells : nc,
        isfinite(sc.cellsize) ? sc.cellsize : cellsize_rec,
        isfinite(sc.blanking) ? sc.blanking : blank_rec,
        csys == 2 ? :beam : csys == 1 ? :xyz : :enu,
        isfinite(sc.velocity_range) ? sc.velocity_range : nanmedian(ambig),
        sc.pressure_offset, sc.declination, sc.salinity,
        Dict{String,Any}("rawConfiguration" => cfgtxt, "source" => "native .ad2cp reader"))

    range = config.blanking .+ config.cellsize .* (1:nc)
    perm = kept[sortperm(time[kept])]
    AD2CPData(time[perm], datetime2unix.(time[perm]), collect(range),
        vel[:, :, perm], amp[:, :, perm], corr[:, :, perm],
        heading[perm], pitch[perm], roll[perm], pressure[perm], temperature[perm],
        soundspeed[perm], accel[:, perm], mag[:, perm], errorv[perm], status[perm],
        ensemble[perm], config, bt)
end
