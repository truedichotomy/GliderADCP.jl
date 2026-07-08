# Layer 1 — Nortek MIDAS netCDF reader (primary input).
#
# MIDAS layout: groups `Config` (≈390 attrs incl. avg_beam2xyz, cell size, blanking,
# coord system, beam θ/φ) and `Data/Average` (+ `Data/AverageBT` when bottom track is
# enabled, `Data/Burst` + `Data/BurstBT` for burst plans). Beam variables are stored as
# (cells × time) per beam. Multi-file missions are split as `*.ad2cp.00000.nc`, `00001`, …
# and are concatenated on time here; Config is taken from the first file
# (same convention as gliderad2cp and Slocum-AD2CP).

# missing→NaN conversions (MIDAS files usually carry no _FillValue, but be defensive)
_f32(A) = Float32.(coalesce.(Array(A), NaN32))
_f64(A) = Float64.(coalesce.(Array(A), NaN))

"""
    _getvec(grp, name, n) -> Vector{Float32}

Read a per-ensemble variable, tolerating its absence (all-NaN fallback).
"""
_getvec(grp, name, n) =
    haskey(grp, name) ? vec(_f32(grp[name][:])) : fill(NaN32, n)

_getvec64(grp, name, n) =
    haskey(grp, name) ? vec(_f64(grp[name][:])) : fill(NaN, n)

"""
    _beamstack(grp, prefix, ncells, ntime) -> Array{Float32,3}

Stack `<prefix>Beam1..4` (each cells × time) into a `ncells × 4 × ntime` array.
"""
function _beamstack(grp, prefix::String, nc::Int, nt::Int)
    A = fill(NaN32, nc, 4, nt)
    for b in 1:4
        name = "$(prefix)Beam$(b)"
        haskey(grp, name) || continue
        A[:, b, :] = _f32(grp[name][:, :])
    end
    return A
end

_attr(attrs, key, default) = haskey(attrs, key) ? attrs[key] : default
_fattr(attrs, key, default) = Float64(_attr(attrs, key, default))

"""
    _parse_config(attrs::Dict, plan::Symbol) -> AD2CPConfig

Build an `AD2CPConfig` from the raw `Config`-group attributes. `plan` selects the
`avg_*` or `burst_*` attribute family.
"""
function _parse_config(attrs::Dict{String,Any}, plan::Symbol)
    p = plan === :burst ? "burst" : "avg"

    # factory beam→XYZ: flat 16 floats, row-major rows X,Y,Z1,Z2
    b2x_raw = _attr(attrs, "$(p)_beam2xyz", _attr(attrs, "beam2xyz", nothing))
    beam2xyz = if b2x_raw === nothing
        fill(NaN, 4, 4)
    else
        Matrix(transpose(reshape(Float64.(vec(b2x_raw)), 4, 4)))
    end

    θ = ntuple(i -> _fattr(attrs, "beamConfiguration$(i)_theta", NaN), 4)
    φ = ntuple(i -> _fattr(attrs, "beamConfiguration$(i)_phi", NaN), 4)

    cs = uppercase(string(_attr(attrs, "$(p)_coordSystem", "BEAM")))
    coordsystem = cs == "BEAM" ? :beam : cs == "XYZ" ? :xyz : :enu

    AD2CPConfig(
        0,  # serial filled by caller from the data group
        _fattr(attrs, "plan_frequency", _fattr(attrs, "headFrequency", NaN)),
        θ, φ, beam2xyz,
        Int(_attr(attrs, "$(p)_nCells", 0)),
        _fattr(attrs, "$(p)_cellSize", NaN),
        _fattr(attrs, "$(p)_blankingDistance", NaN),
        coordsystem,
        _fattr(attrs, "$(p)_velocityRange", NaN),
        _fattr(attrs, "pressureOffset", NaN),
        _fattr(attrs, "user_decl", NaN),
        _fattr(attrs, "salinity", NaN),
        attrs,
    )
end

"""
    _load_bt(grp) -> BottomTrackData

Read a `Data/AverageBT` (or `BurstBT`) group.
"""
function _load_bt(grp)
    rawtime = grp["time"][:]
    good = findall(!ismissing, rawtime)
    time = DateTime.(rawtime[good])
    nraw = length(rawtime)
    nt = length(good)
    vel = fill(NaN32, 4, nt)
    dist = fill(NaN, 4, nt)
    fom = fill(NaN32, 4, nt)
    for b in 1:4
        vel[b, :] = _getvec(grp, "VelocityBeam$(b)", nraw)[good]
        dist[b, :] = _getvec64(grp, "DistanceBeam$(b)", nraw)[good]
        fom[b, :] = _getvec(grp, "FOMBeam$(b)", nraw)[good]
    end
    BottomTrackData(time, datetime2unix.(time), vel, dist, fom,
        _getvec64(grp, "Pressure", nraw)[good],
        _getvec(grp, "Heading", nraw)[good],
        _getvec(grp, "Pitch", nraw)[good],
        _getvec(grp, "Roll", nraw)[good],
        _getvec(grp, "SpeedOfSound", nraw)[good])
end

# per-file payload prior to concatenation
function _load_one(path::AbstractString, plan::Symbol)
    NCDataset(path) do ds
        haskey(ds.group, "Data") || error("$path: no `Data` group — not a MIDAS AD2CP export")
        data = ds.group["Data"]
        gname = plan === :burst ? "Burst" : "Average"
        if !haskey(data.group, gname)
            avail = collect(keys(data.group))
            error("$path: no `Data/$gname` group (found: $(join(avail, ", "))). " *
                  "Pass plan=:burst for burst-only files.")
        end
        grp = data.group[gname]

        rawtime = grp["time"][:]
        good = findall(!ismissing, rawtime)
        time = DateTime.(rawtime[good])
        nraw = length(rawtime)

        range = vec(_f64(grp["Velocity Range"][:]))
        nc = length(range)

        sub3(A) = A[:, :, good]
        payload = (
            time = time,
            range = range,
            vel = sub3(_beamstack(grp, "Velocity", nc, nraw)),
            amp = sub3(_beamstack(grp, "Amplitude", nc, nraw)),
            corr = sub3(_beamstack(grp, "Correlation", nc, nraw)),
            heading = _getvec(grp, "Heading", nraw)[good],
            pitch = _getvec(grp, "Pitch", nraw)[good],
            roll = _getvec(grp, "Roll", nraw)[good],
            pressure = _getvec64(grp, "Pressure", nraw)[good],
            temperature = _getvec(grp, "WaterTemperature", nraw)[good],
            soundspeed = _getvec(grp, "SpeedOfSound", nraw)[good],
            accel = permutedims(hcat((_getvec(grp, "Accelerometer$a", nraw)[good] for a in ("X", "Y", "Z"))...)),
            mag = permutedims(hcat((_getvec(grp, "Magnetometer$a", nraw)[good] for a in ("X", "Y", "Z"))...)),
            error = _getvec64(grp, "Error", nraw)[good],
            status = _getvec64(grp, "Status", nraw)[good],
            ensemble = _getvec64(grp, "EnsembleCount", nraw)[good],
            serial = haskey(grp, "SerialNumber") ? Int(coalesce(grp["SerialNumber"][1], 0)) : 0,
        )

        btname = gname * "BT"
        bt = haskey(data.group, btname) ? _load_bt(data.group[btname]) : nothing

        attrs = haskey(ds.group, "Config") ? Dict{String,Any}(ds.group["Config"].attrib) :
                Dict{String,Any}()

        return payload, bt, attrs
    end
end

"""
    load_ad2cp(paths; plan=:average) -> AD2CPData

Read one or more Nortek MIDAS-exported AD2CP netCDF files (`*.ad2cp.NNNNN.nc`),
concatenate along time, sort, and attach the instrument configuration and bottom-track
records (when present).

`paths` may be a single file, a vector of files, or a directory (all `*.nc` files whose
name contains `.ad2cp.` are used, in name order). `plan` selects the `Data/Average`
(default) or `Data/Burst` group.

```julia
adcp = load_ad2cp("sea064_M38.ad2cp.00000.nc")
adcp.vel        # ncells × 4 beams × ntime, beam coordinates
adcp.bt         # BottomTrackData (M38 has bottom track enabled)
```
"""
function load_ad2cp(paths; plan::Symbol=:average)
    # native binary files dispatch to the pure-Julia reader (no MIDAS required)
    if paths isa AbstractString && isfile(paths) && endswith(lowercase(paths), ".ad2cp")
        return read_ad2cp(paths; plan)
    end
    files = _ad2cp_files(paths)
    isempty(files) && error("load_ad2cp: no input files found for $paths")

    loaded = []
    okfiles = String[]
    for f in files
        try
            push!(loaded, _load_one(f, plan))
            push!(okfiles, f)
        catch err
            @warn "load_ad2cp: skipping unreadable file" file = basename(f) error = sprint(showerror, err)
        end
    end
    isempty(loaded) &&
        error("load_ad2cp: no readable files (all $(length(files)) failed)")
    files = okfiles
    payloads = first.(loaded)

    # consistency across files
    range = payloads[1].range
    for (f, p) in zip(files[2:end], payloads[2:end])
        p.range == range ||
            error("load_ad2cp: cell ranges differ between files ($(files[1]) vs $f)")
    end

    cat3(field) = cat((getfield(p, field) for p in payloads)...; dims=3)
    cat1(field) = reduce(vcat, (getfield(p, field) for p in payloads))
    cat2(field) = hcat((getfield(p, field) for p in payloads)...)

    time = cat1(:time)
    perm = sortperm(time)

    bts = [bt for (_, bt, _) in loaded if bt !== nothing]
    bt = if isempty(bts)
        nothing
    else
        btime = reduce(vcat, (b.time for b in bts))
        bperm = sortperm(btime)
        BottomTrackData(btime[bperm], datetime2unix.(btime[bperm]),
            hcat((b.vel for b in bts)...)[:, bperm],
            hcat((b.distance for b in bts)...)[:, bperm],
            hcat((b.fom for b in bts)...)[:, bperm],
            reduce(vcat, (b.pressure for b in bts))[bperm],
            reduce(vcat, (b.heading for b in bts))[bperm],
            reduce(vcat, (b.pitch for b in bts))[bperm],
            reduce(vcat, (b.roll for b in bts))[bperm],
            reduce(vcat, (b.soundspeed for b in bts))[bperm])
    end

    attrs = loaded[1][3]
    config0 = _parse_config(attrs, plan)
    config = AD2CPConfig(payloads[1].serial, config0.frequency, config0.beam_theta,
        config0.beam_phi, config0.beam2xyz,
        config0.ncells == 0 ? length(range) : config0.ncells,
        config0.cellsize, config0.blanking, config0.coordsystem, config0.velocity_range,
        config0.pressure_offset, config0.declination, config0.salinity_setting, attrs)

    AD2CPData(time[perm], datetime2unix.(time[perm]), range,
        cat3(:vel)[:, :, perm], cat3(:amp)[:, :, perm], cat3(:corr)[:, :, perm],
        cat1(:heading)[perm], cat1(:pitch)[perm], cat1(:roll)[perm],
        cat1(:pressure)[perm], cat1(:temperature)[perm], cat1(:soundspeed)[perm],
        cat2(:accel)[:, perm], cat2(:mag)[:, perm],
        cat1(:error)[perm], cat1(:status)[perm], cat1(:ensemble)[perm],
        config, bt)
end

# path expansion: file | vector of files | directory containing *.ad2cp.*.nc
function _ad2cp_files(paths::AbstractString)
    if isdir(paths)
        fs = filter(f -> occursin(".ad2cp.", f) && endswith(f, ".nc"), readdir(paths))
        return joinpath.(paths, sort(fs))
    end
    isfile(paths) && return [String(paths)]
    error("load_ad2cp: $paths is neither a file nor a directory")
end
_ad2cp_files(paths::AbstractVector) = String.(paths)
