using GliderADCP
using Test
using Dates
using DataFrames
using CSV
using NCDatasets
using CodecZlib
using Statistics

# Reference data (present on the development machine only; testsets skip elsewhere)
const M38_DIR = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete"
const M38_NC = joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp.00000.nc")
const M38_NAV = joinpath(M38_DIR, "delayed/nav/logs")
const M38_PLD = joinpath(M38_DIR, "delayed/pld1/logs")
const M48_NC = "/Users/gong/GitHub/jlglider/ad2cp/data/sea064_M48.ad2cp.00000.nc"

"Write a miniature MIDAS-style AD2CP netCDF (Config + Data/Average + Data/AverageBT)."
function write_synthetic_midas(path; nt=6, nc=3)
    t0 = DateTime(2022, 11, 4)
    NCDataset(path, "c") do ds
        cfg = defGroup(ds, "Config")
        cfg.attrib["avg_beam2xyz"] = Float32[0.6782, 0, -0.6782, 0,
                                             0, -1.1831, 0, 1.1831,
                                             0.74, 0, 0.74, 0,
                                             0, 0.5518, 0, 0.5518]
        cfg.attrib["avg_nCells"] = Int32(nc)
        cfg.attrib["avg_cellSize"] = 2.0f0
        cfg.attrib["avg_blankingDistance"] = 0.7
        cfg.attrib["avg_coordSystem"] = "BEAM"
        cfg.attrib["avg_velocityRange"] = 2.5f0
        cfg.attrib["pressureOffset"] = 9.5f0
        cfg.attrib["user_decl"] = 0.0f0
        cfg.attrib["salinity"] = 38.0f0
        cfg.attrib["plan_frequency"] = Int32(1000)
        for (i, (th, ph)) in enumerate(zip((47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0)))
            cfg.attrib["beamConfiguration$(i)_theta"] = Float32(th)
            cfg.attrib["beamConfiguration$(i)_phi"] = Float32(ph)
        end
        data = defGroup(ds, "Data")
        avg = defGroup(data, "Average")
        defDim(avg, "time", nt)
        defDim(avg, "Velocity Range", nc)
        defVar(avg, "time", t0 .+ Second.(10 .* (0:nt-1)), ("time",))
        defVar(avg, "Velocity Range", Float32.(0.7 .+ 2.0 .* (1:nc)), ("Velocity Range",))
        for b in 1:4
            defVar(avg, "VelocityBeam$b", Float32.(b .* reshape(1:nc*nt, nc, nt) ./ 100),
                ("Velocity Range", "time"))
            defVar(avg, "CorrelationBeam$b", fill(Float32(90 + b), nc, nt),
                ("Velocity Range", "time"))
            defVar(avg, "AmplitudeBeam$b", fill(Float32(70 + b), nc, nt),
                ("Velocity Range", "time"))
        end
        defVar(avg, "Heading", Float32.(collect(1.0:nt)), ("time",))
        defVar(avg, "Pitch", fill(-17.5f0, nt), ("time",))
        defVar(avg, "Roll", zeros(Float32, nt), ("time",))
        defVar(avg, "Pressure", Float64.(collect(10.0:10.0+nt-1)), ("time",))
        defVar(avg, "SpeedOfSound", fill(1485.0f0, nt), ("time",))
        defVar(avg, "WaterTemperature", fill(7.5f0, nt), ("time",))
        defVar(avg, "SerialNumber", fill(102381.0, nt), ("time",))
        defVar(avg, "AccelerometerZ", fill(-0.95f0, nt), ("time",))
        bt = defGroup(data, "AverageBT")
        defDim(bt, "time", nt - 1)
        defVar(bt, "time", t0 .+ Second.(10 .* (0:nt-2)) .+ Millisecond(558), ("time",))
        for b in 1:4
            defVar(bt, "VelocityBeam$b", fill(Float32(0.1b), nt - 1), ("time",))
            defVar(bt, "DistanceBeam$b", fill(5.0 + b, nt - 1), ("time",))
            defVar(bt, "FOMBeam$b", fill(Float32(100b), nt - 1), ("time",))
        end
        defVar(bt, "Pressure", fill(50.0, nt - 1), ("time",))
        defVar(bt, "Heading", fill(90.0f0, nt - 1), ("time",))
        defVar(bt, "Pitch", fill(-17.0f0, nt - 1), ("time",))
        defVar(bt, "Roll", fill(0.5f0, nt - 1), ("time",))
        defVar(bt, "SpeedOfSound", fill(1485.0f0, nt - 1), ("time",))
    end
    return path
end

"Write a miniature gzipped SeaExplorer gli segment."
function write_synthetic_gli(path; nrows=4)
    hdr = "Timestamp;NavState;SecurityLevel;Heading;Declination;Pitch;Roll;Depth;Temperature;Pa;Lat;Lon;DeadReckoning;DesiredH;BallastCmd;BallastPos;LinCmd;LinPos;AngCmd;AngPos;Voltage;Altitude;"
    rows = String[]
    for i in 1:nrows
        dr = i < nrows ? 1 : 0
        push!(rows,
            "04/11/2022 06:23:$(20+i);117;0;128.06;0;21.32;8.05;$(5.0-i);6.90;72434;7001.296;201.504;$dr;125;275;277.4;36.4;36.2;29;29.8;28.9;-1.0;")
    end
    open(path, "w") do io
        gz = GzipCompressorStream(io)
        write(gz, join([hdr; rows], "\n") * "\n")
        close(gz)
    end
    return path
end

"Synthetic .ad2cp records: (string-config record, one DF3 average record)."
function synthetic_ad2cp_records()
    cs(b) = GliderADCP._ad2cp_checksum(b)
    le16(x) = reinterpret(UInt8, [UInt16(x)])
    le32(x) = reinterpret(UInt8, [UInt32(x)])
    function mkrec(id, payload)
        hdr = UInt8[0xA5, 0x0A, id, 0x10, le16(length(payload))...,
            le16(cs(payload))..., 0x00, 0x00]
        hdr[9:10] = le16(cs(hdr[1:8]))
        vcat(hdr, payload)
    end
    cfgtxt = "ID,STR=\"Glider\",SN=42\nGETPLAN,SA=35.0,FREQ=1000\n" *
             "GETAVG,NC=2,CS=1.00,BD=0.50,CY=\"BEAM\",VR=2.50\n" *
             "GETUSER,POFF=1.50,DECL=3.00\n" *
             "BEAMCFGLIST,BEAM=1,THETA=47.50,PHI=0.00\n" *
             "BEAMCFGLIST,BEAM=2,THETA=25.00,PHI=-90.00\n"
    strrec = mkrec(0xA0, vcat(UInt8[0x10], Vector{UInt8}(cfgtxt), UInt8[0x00]))
    p = zeros(UInt8, 76)
    p[1] = 3; p[2] = 76                              # version, offsetOfData
    p[3:4] = le16(0x00E0)                            # bits 5,6,7
    p[5:8] = le32(42)                                # serial
    p[9:14] = UInt8[122, 10, 4, 6, 0, 0]             # 2022-11-04T06:00:00
    p[15:16] = le16(4380)                            # .438 s
    p[17:18] = le16(15000)                           # sound speed 1500.0
    p[19:20] = le16(750)                             # 7.50 °C
    p[21:24] = le32(12345)                           # 12.345 dbar
    p[25:26] = le16(9000)                            # heading 90.00°
    p[27:28] = reinterpret(UInt8, [Int16(-1750)])    # pitch −17.50°
    p[29:30] = reinterpret(UInt8, [Int16(100)])      # roll 1.00°
    p[31:32] = le16(UInt16(2) | (UInt16(2) << 10) | (UInt16(4) << 12))
    p[33:34] = le16(1000)                            # cell size 1.000 m
    p[35:36] = le16(500)                             # blanking 0.500 m (mm scaling)
    p[59] = reinterpret(UInt8, Int8(-3))             # velocity scaling 10^-3
    p[73:76] = le32(7)                               # ensemble counter
    velraw = Int16[100, 200, -300, -400, 500, 600, -700, -800]
    amps = UInt8[140, 142, 144, 146, 148, 150, 152, 154]
    corrs = UInt8[90, 91, 92, 93, 94, 95, 96, 97]
    payload = vcat(p, reinterpret(UInt8, velraw), amps, corrs)
    return strrec, mkrec(0x16, payload)
end

@testset "GliderADCP.jl" begin

    @testset "nmea2deg" begin
        @test nmea2deg(7001.296) ≈ 70 + 1.296 / 60
        @test nmea2deg(201.504) ≈ 2 + 1.504 / 60
        @test nmea2deg(-3830.5) ≈ -(38 + 30.5 / 60)
        @test nmea2deg(0.0) == 0.0
        @test isnan(nmea2deg(NaN))
    end

    @testset "seaexplorer_files natural sort" begin
        mktempdir() do d
            for n in (1, 2, 10, 100)
                touch(joinpath(d, "sea064.38.gli.sub.$n.gz"))
            end
            touch(joinpath(d, "sea.0.gli.evt.1"))            # must be ignored
            touch(joinpath(d, "sea064.38.pld1.raw.3.gz"))    # different stream
            fs = seaexplorer_files(d, "gli.sub")
            @test length(fs) == 4
            @test [match(r"\.(\d+)\.gz$", f).captures[1] for f in fs] == ["1", "2", "10", "100"]
        end
    end

    @testset "synthetic MIDAS netCDF round-trip" begin
        mktempdir() do d
            f = write_synthetic_midas(joinpath(d, "syn.ad2cp.00000.nc"))
            a = load_ad2cp(f)
            @test length(a) == 6
            @test ncells(a) == 3
            @test a.range ≈ [2.7, 4.7, 6.7] atol = 1e-5
            @test a.config.serial == 102381
            @test a.config.cellsize == 2.0
            @test a.config.blanking ≈ 0.7
            @test a.config.coordsystem === :beam
            @test a.config.velocity_range ≈ 2.5
            @test a.config.beam_theta == (47.5, 25.0, 47.5, 25.0)
            @test a.config.beam2xyz[1, :] ≈ [0.6782, 0, -0.6782, 0] atol = 1e-6
            @test a.config.beam2xyz[2, 2] ≈ -1.1831 atol = 1e-6
            @test size(a.vel) == (3, 4, 6)
            @test a.vel[1, 2, 1] ≈ 2 * 1 / 100          # beam 2, first cell, first ping
            @test a.pitch[1] ≈ -17.5f0
            @test all(isnan, a.mag)                     # absent variable → NaN fill
            @test a.bt !== nothing
            @test length(a.bt) == 5
            @test a.bt.vel[3, 1] ≈ 0.3f0
            @test a.bt.distance[4, 2] ≈ 9.0
            # directory form finds the file too
            a2 = load_ad2cp(d)
            @test length(a2) == 6
            # multi-file concatenation + sorting
            f2 = write_synthetic_midas(joinpath(d, "syn.ad2cp.00001.nc"))
            a3 = load_ad2cp([f, f2])
            @test length(a3) == 12
            @test issorted(a3.time)
        end
    end

    @testset "synthetic SeaExplorer nav" begin
        mktempdir() do d
            write_synthetic_gli(joinpath(d, "sea064.38.gli.sub.1.gz"); nrows=4)
            write_synthetic_gli(joinpath(d, "sea064.38.gli.sub.2.gz"); nrows=3)
            nav = load_seaexplorer_nav(d)
            @test length(nav) == 7
            @test issorted(nav.time)
            @test nav.lat[1] ≈ 70 + 1.296 / 60
            @test nav.lon[1] ≈ 2 + 1.504 / 60
            @test count(==(0), nav.deadreckoning) == 2
            @test nav.navstate[1] == 117
            @test hasproperty(nav.df, :Voltage)  # full table retained
        end
    end

    @testset "geometry: factory-matrix consistency" begin
        e = beam_unit_vectors((47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0))
        E4 = permutedims(hcat(e...))              # 4×3 forward model rows
        P = (E4' * E4) \ E4'                      # 3×4 least-squares inverse
        @test P[1, :] ≈ [0.6782, 0, -0.6782, 0] atol = 2e-4   # factory X row
        @test P[2, :] ≈ [0, -1.1831, 0, 1.1831] atol = 2e-4   # factory Y row
    end

    @testset "geometry: synthetic round-trip (3-beam exact)" begin
        e = beam_unit_vectors((47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0))
        F = head2vehicle(:down)
        for (h, p, r) in ((0.0, -17.5, 0.0), (123.0, -15.0, 8.0), (280.0, 20.0, -5.0), (45.0, 17.0, 3.0))
            R = rotmat_xyz2enu(h, p, r)
            v_enu = [0.3, -0.1, 0.05]
            b = [(R * F * e[i])' * v_enu for i in 1:4]
            sel = select_beams(p; look=:down)
            @test sel == (p < 0 ? (1, 2, 4) : (2, 3, 4))
            v = xyz_from_beams(b[collect(sel)], sel, e)
            @test R * F * v ≈ v_enu atol = 1e-10
            # 4-beam least-squares solve recovers exactly for consistent data
            v4 = xyz_from_beams(b, (1, 2, 3, 4), e)
            @test R * F * v4 ≈ v_enu atol = 1e-10
        end
        # declination is equivalent to a heading offset
        @test rotmat_xyz2enu(100.0, -17.0, 2.0; declination=7.5) ≈
              rotmat_xyz2enu(107.5, -17.0, 2.0) atol = 1e-12
        # Slocum-AD2CP row-drop parity mode: halves X for a pure-X flow (documented defect)
        vh = [0.4, 0.0, 0.0]
        b = [e[i]' * vh for i in 1:4]
        vr = xyz_from_beams(b[[1, 2, 4]], (1, 2, 4), e; method=:rowdrop)
        @test vr[1] ≈ 0.5 * vh[1] atol = 1e-6
        @test xyz_from_beams(b[[1, 2, 4]], (1, 2, 4), e)[1] ≈ vh[1] atol = 1e-10
    end

    @testset "geometry: vertical cosines & offsets" begin
        cfg = AD2CPConfig(1, 1000.0, (47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0),
            fill(NaN, 4, 4), 15, 2.0, 0.7, :beam, 2.5, 0.0, 0.0, 38.0, Dict{String,Any}())
        t0 = DateTime(2022, 11, 4)
        mk(h, p, r) = AD2CPData([t0], [datetime2unix(t0)], collect(2.7:2.0:30.7),
            fill(0.1f0, 15, 4, 1), fill(80.0f0, 15, 4, 1), fill(90.0f0, 15, 4, 1),
            [Float32(h)], [Float32(p)], [Float32(r)], [100.0], [7.0f0], [1485.0f0],
            reshape(Float32[0, 0, -0.95], 3, 1), fill(NaN32, 3, 1),
            [0.0], [0.0], [1.0], cfg, nothing)
        a = mk(90.0, -17.5, 0.0)
        vc = vertical_cosines(a; look=:down)
        @test vc[1, 1] ≈ cosd(30.0) atol = 1e-6      # fore beam on a dive: 47.5−17.5
        @test vc[3, 1] ≈ cosd(65.0) atol = 1e-6      # aft beam: 47.5+17.5
        @test vc[2, 1] ≈ cosd(17.5) * cosd(25.0) atol = 1e-6   # side beams
        @test vc[4, 1] ≈ cosd(17.5) * cosd(25.0) atol = 1e-6
        g = offset_grid(cfg)
        @test g[1] == 0.0
        @test g[2] - g[1] == 1.0                     # cellsize/2
        @test g[end] == 31.0
        @test detect_look_direction(a) === :down
    end

    @testset "regrid + ENU on isobars (uniform flow)" begin
        cfg = AD2CPConfig(1, 1000.0, (47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0),
            fill(NaN, 4, 4), 15, 2.0, 0.7, :beam, 2.5, 0.0, 0.0, 38.0, Dict{String,Any}())
        e = beam_unit_vectors(cfg)
        F = head2vehicle(:down)
        t0 = DateTime(2022, 11, 4)
        nt = 5
        v_enu = [0.25, -0.15, 0.03]
        h, p, r = 240.0, -16.0, 4.0
        R = rotmat_xyz2enu(h, p, r)
        vel = Array{Float32}(undef, 15, 4, nt)
        for b in 1:4
            vel[:, b, :] .= Float32((R * F * e[b])' * v_enu)
        end
        tv = t0 .+ Second.(10 .* (0:nt-1))
        a = AD2CPData(tv, datetime2unix.(tv),
            collect(2.7:2.0:30.7), vel, fill(80.0f0, 15, 4, nt), fill(90.0f0, 15, 4, nt),
            fill(Float32(h), nt), fill(Float32(p), nt), fill(Float32(r), nt),
            fill(100.0, nt), fill(7.0f0, nt), fill(1485.0f0, nt),
            repeat(Float32[0, 0, -0.95], 1, nt), fill(NaN32, 3, nt),
            zeros(nt), zeros(nt), collect(1.0:nt), cfg, nothing)
        E, N, U, offs, used = enu_on_isobars(a)
        @test used[1] == (1, 2, 4)
        finiteE = filter(isfinite, E)
        @test !isempty(finiteE)
        @test all(x -> isapprox(x, v_enu[1]; atol=1e-5), finiteE)
        @test all(x -> isapprox(x, v_enu[2]; atol=1e-5), filter(isfinite, N))
        @test all(x -> isapprox(x, v_enu[3]; atol=1e-5), filter(isfinite, U))
        # beams_to_enu (native cell grid) agrees
        E2, N2, U2, _ = beams_to_enu(a)
        @test all(x -> isapprox(x, v_enu[1]; atol=1e-5), filter(isfinite, E2))
    end

    @testset "qc! masks and stats" begin
        mktempdir() do d
            f = write_synthetic_midas(joinpath(d, "syn.ad2cp.00000.nc"); nt=6, nc=3)
            a = load_ad2cp(f)
            a.corr[1, 1, 1] = 30.0f0            # below 50 → rejected
            a.amp[2, 2, 2] = 90.0f0             # above 75 → rejected
            a.vel[3, 3, 3] = 1.5f0              # |v| > 0.8 → rejected
            stats = qc!(a; thr=QCThresholds(snr_db=NaN, surface_depth=NaN, first_cells=0))
            @test stats.correlation > 0
            @test stats.amplitude_max > 0
            @test stats.velocity_max > 0
            @test isnan(a.vel[1, 1, 1]) && isnan(a.vel[2, 2, 2]) && isnan(a.vel[3, 3, 3])
            @test stats.total ≥ stats.correlation + stats.amplitude_max
        end
    end

    @testset "soundspeed correction" begin
        mktempdir() do d
            f = write_synthetic_midas(joinpath(d, "syn.ad2cp.00000.nc"))
            a = load_ad2cp(f)
            v0 = copy(a.vel)
            scale = soundspeed_correction(a, a.t, fill(1485.0 * 1.01, length(a)))
            @test all(s -> isapprox(s, 1.01; atol=1e-6), scale)
            apply_soundspeed!(a, scale)
            @test a.vel ≈ 1.01f0 .* v0
            # TEOS-10 sanity: warmer/saltier → faster
            @test soundspeed_from_ctd(35.0, 10.0, 100.0, 5.0, 70.0) >
                  soundspeed_from_ctd(34.0, 8.0, 100.0, 5.0, 70.0)
        end
    end

    @testset "compute_dac: synthetic current recovery" begin
        # Build a nav track: fix → 1 h submerged DR (currents ignored) → fix displaced
        # by current × duration → next dive. DAC must recover the prescribed current.
        u_true, v_true = 0.12, -0.07
        lat0, lon0 = 70.0, 2.0
        t0 = DateTime(2022, 11, 4)
        times = DateTime[]; lons = Float64[]; lats = Float64[]
        dr = Int8[]; depth = Float64[]
        # pre-dive fix
        push!(times, t0); push!(lons, lon0); push!(lats, lat0); push!(dr, 0); push!(depth, 0.0)
        # submerged: DR walks east at 0.3 m/s through water, 3600 s, 60 records
        T = 3600.0
        for k in 1:60
            t = k * T / 60
            push!(times, t0 + Second(round(Int, t)))
            push!(lons, lon0 + rad2deg(0.3t / (6.371e6 * cosd(lat0))))
            push!(lats, lat0)
            push!(dr, 1)
            push!(depth, 100.0)
        end
        # surfacing fix: DR endpoint + current drift over T (+30 s fix delay)
        dx, dy = u_true * (T + 30), v_true * (T + 30)
        lon_dr, lat_dr = lons[end], lats[end]
        push!(times, t0 + Second(round(Int, T + 30)))
        push!(lons, lon_dr + rad2deg(dx / (6.371e6 * cosd(lat0))))
        push!(lats, lat_dr + rad2deg(dy / 6.371e6))
        push!(dr, 0); push!(depth, 0.0)
        nav = GliderNav(times, datetime2unix.(times), lons, lats,
            fill(90.0, length(times)), zeros(length(times)), zeros(length(times)),
            zeros(length(times)), depth, fill(Int16(110), length(times)), dr,
            fill(-1.0, length(times)), DataFrame())
        dac = compute_dac(nav; min_duration=100.0)
        @test nrow(dac) == 1
        @test dac.u[1] ≈ u_true atol = 1e-3
        @test dac.v[1] ≈ v_true atol = 1e-3
        @test dac.duration[1] ≈ T + 30 atol = 1.0
        @test dac.maxdepth[1] == 100.0
        # too-short segments are rejected by default QC
        @test nrow(compute_dac(nav; min_duration=7200.0)) == 0
    end

    @testset "end-to-end structure orientation (depth-varying flow through geometry)" begin
        # A pure baroclinic sign flip preserves depth means, dive/climb consistency and
        # DAC closure — so this dedicated test drives a depth-VARYING flow through the
        # full beam forward model and requires the recovered structure to be upright.
        θ = (47.5, 25.0, 47.5, 25.0)
        φ = (0.0, -90.0, 180.0, 90.0)
        e = beam_unit_vectors(θ, φ)
        F = head2vehicle(:down)
        cfg = AD2CPConfig(1, 1000.0, θ, φ, fill(NaN, 4, 4), 15, 2.0, 0.7, :beam, 2.5,
            0.0, 0.0, 38.0, Dict{String,Any}())
        u_o(z) = 0.001 * z
        ug = 0.25
        nt = 400
        t0 = DateTime(2022, 11, 4)
        depth = [10 + 180 * (1 - abs(1 - 2i / nt)) for i in 1:nt]
        pitchv = [i <= nt ÷ 2 ? -17.0 : 17.0 for i in 1:nt]
        rng = collect(2.7:2.0:30.7)
        vel = fill(NaN32, 15, 4, nt)
        for i in 1:nt
            R = rotmat_xyz2enu(90.0, pitchv[i], 0.0)
            for b in 1:4
                eb = R * F * e[b]
                for k in 1:15
                    z = depth[i] + (rng[k] / cosd(25.0)) * (-eb[3])
                    vel[k, b, i] = Float32(eb[1] * (u_o(z) - ug))
                end
            end
        end
        tv = t0 .+ Second.(10 .* (0:nt-1))
        a = AD2CPData(tv, datetime2unix.(tv), rng, vel, fill(80f0, 15, 4, nt),
            fill(90f0, 15, 4, nt), fill(90f0, nt), Float32.(pitchv), fill(0f0, nt),
            Float64.(depth), fill(7f0, nt), fill(1500f0, nt),
            repeat(Float32[0, 0, -0.95], 1, nt), fill(NaN32, 3, nt),
            zeros(nt), zeros(nt), collect(1.0:nt), cfg, nothing)
        pp = process_pings(a; lat=0.0)
        errs = Float64[]
        for i in 1:nt, k in 1:length(pp.offsets)
            isfinite(pp.E[k, i]) || continue
            push!(errs, pp.E[k, i] - (u_o(pp.celldepth[k, i]) - ug))
        end
        @test length(errs) > 5000
        @test maximum(abs.(errs)) < 0.005
        dacdf = DataFrame(yo=[1], t_start=[tv[1]], t_end=[tv[end]], t_mid=[tv[nt÷2]],
            u=[mean(u_o.(5:10:185.0))], v=[0.0])
        inv = solve_inverse(pp, dacdf)
        g = inv.nobs .> 20
        slope = cov(inv.z[g], inv.u[g]) / var(inv.z[g])
        @test 0.0007 < slope < 0.0013            # upright, right magnitude
        @test maximum(abs.(inv.u[g] .- u_o.(inv.z[g]))) < 0.03
    end

    @testset "solvers: synthetic truth recovery" begin
        # One dive+climb segment with a prescribed ocean profile and glider velocity.
        u_o(z) = 0.2 * sin(2π * z / 150) + 0.05
        v_o(z) = -0.1 + 0.0008 * z
        nt = 900
        t0 = DateTime(2022, 11, 4)
        times = t0 .+ Second.(10 .* (0:nt-1))
        tunix = datetime2unix.(times)
        u_g = [0.3 * cos(2π * i / nt) - 0.1 for i in 1:nt]
        v_g = [0.25 * sin(4π * i / nt) for i in 1:nt]
        depth = [5 + 195 * (1 - abs(1 - 2i / nt)) for i in 1:nt]   # 5 → 200 → 5 m
        offsets = collect(0.0:1.0:30.0)
        celldepth = offsets .+ depth'
        E = [u_o(celldepth[k, i]) - u_g[i] for k in 1:31, i in 1:nt]
        N = [v_o(celldepth[k, i]) - v_g[i] for k in 1:31, i in 1:nt]
        pp = ProcessedPings(times, tunix, depth, fill(90.0, nt), offsets, E, N,
            zeros(31, nt), celldepth, :down, fill((1, 2, 4), nt))

        centers = 5.0:10.0:195.0                       # DAC over glider-covered bins
        dacu, dacv = mean(u_o.(centers)), mean(v_o.(centers))
        dacdf = DataFrame(yo=[1], t_start=[times[1]], t_end=[times[end]],
            t_mid=[times[nt÷2]], u=[dacu], v=[dacv])

        # --- inverse, DAC-constrained ---
        inv1 = solve_inverse(pp, dacdf)
        @test !isempty(inv1)
        good = inv1.nobs .> 20
        # residual is 10-m bin discretization of the 150-m-wavelength profile
        @test maximum(abs.(inv1.u[good] .- u_o.(inv1.z[good]))) < 0.03
        @test maximum(abs.(inv1.v[good] .- v_o.(inv1.z[good]))) < 0.03

        # --- inverse, platform-form DAC ---
        inv2 = solve_inverse(pp, dacdf; opts=InverseOptions(dac_form=:platform))
        # platform form constrains mean glider velocity = DAC; with our synthetic
        # (mean u_g ≠ dac) the profile shifts by that difference — check structure only
        g2 = inv2.nobs .> 20
        du = inv2.u[g2] .- u_o.(inv2.z[g2])
        @test std(du) < 0.02                            # shape right, offset allowed

        # --- inverse, bottom-track only (no DAC): absolute without GPS reference ---
        btdf = DataFrame(t=tunix[1:40:end], u=u_g[1:40:end], v=v_g[1:40:end])
        inv3 = solve_inverse(pp, dacdf; bt=btdf,
            opts=InverseOptions(wdac=0.0, wbt=5.0))
        g3 = inv3.nobs .> 20
        @test maximum(abs.(inv3.u[g3] .- u_o.(inv3.z[g3]))) < 0.03
        @test maximum(abs.(inv3.v[g3] .- v_o.(inv3.z[g3]))) < 0.03
        @test inv3.nbt[1] == length(tunix[1:40:end])

        # --- glider velocity recovery ---
        sol = invert_segment(E, N, celldepth, tunix, maximum(depth);
            dacu, dacv, opts=InverseOptions())
        @test cor(sol.ug, u_g) > 0.995
        @test maximum(abs.(sol.ug .- u_g)) < 0.03

        # --- shear method ---
        sh = solve_shear(pp, dacdf)
        @test !isempty(sh)
        gs = (sh.nobs .>= 4) .&& (sh.z .< 200)
        @test maximum(abs.(sh.u[gs] .- u_o.(sh.z[gs]))) < 0.06
        @test maximum(abs.(sh.v[gs] .- v_o.(sh.z[gs]))) < 0.06

        # shear and inverse agree with each other
        both = intersect(inv1.z[good], sh.z[gs])
        iu = Dict(zip(inv1.z, inv1.u)); su = Dict(zip(sh.z, sh.u))
        @test maximum(abs(iu[z] - su[z]) for z in both) < 0.06

        # --- shear-content products (Task 1) ---
        dudz(z) = 0.2 * (2π / 150) * cos(2π * z / 150)
        shp = solve_shear_profile(pp, dacdf)
        gp = (shp.nobs .>= 4) .&& isfinite.(shp.sh_u) .&& (shp.z .< 200)
        @test count(gp) > 10
        @test maximum(abs.(shp.sh_u[gp] .- dudz.(shp.z[gp]))) < 1.5e-3
        ish = inverse_shear(inv1)
        gi = (ish.nobs .> 20) .&& isfinite.(ish.sh_u)
        @test count(gi) > 10
        # smoothness + 2·dz centered differencing attenuate; shape must match
        @test cor(ish.sh_u[gi], dudz.(ish.z[gi])) > 0.98
        @test maximum(abs.(ish.sh_u[gi] .- dudz.(ish.z[gi]))) < 3e-3
        # and the two shear products agree on this noiseless synthetic
        jj = innerjoin(shp[gp, :], ish[gi, :]; on=[:yo, :z], makeunique=true)
        @test nrow(jj) > 8
        @test cor(jj.sh_u, jj.sh_u_1) > 0.97
    end

    @testset "declination, gridding, export" begin
        # IGRF sanity at two known sites
        @test 0.5 < magnetic_declination(69.7, 5.0, DateTime(2022, 11, 15)) < 4.5
        @test abs(magnetic_declination(40.7, -74.0, DateTime(2022, 6, 1)) - (-12.7)) < 3

        prof = DataFrame(
            yo=[1, 1, 2, 2],
            t_mid=DateTime(2022, 11, 4) .+ Hour.([0, 0, 2, 2]),
            z=[5.0, 15.0, 15.0, 25.0],
            u=[0.1, 0.2, 0.3, 0.4], v=[-0.1, -0.2, -0.3, -0.4],
            nobs=[5, 6, 7, 8])
        sec = grid_profiles(prof)
        @test sec.z == [5.0, 15.0, 25.0]
        @test length(sec.t) == 2
        @test sec.U[1, 1] == 0.1 && sec.U[2, 2] == 0.3
        @test isnan(sec.U[3, 1])                    # depth-matched, not row-index
        @test sec.Nobs[3, 2] == 8

        mktempdir() do dd
            f = joinpath(dd, "sec.nc")
            export_sections(f, sec; attrs=Dict{String,Any}("mission" => "test"))
            NCDataset(f) do ds
                @test ds.attrib["mission"] == "test"
                @test occursin("GliderADCP.jl", ds.attrib["source"])
                @test size(ds["u"]) == (3, 2)
                @test ds["u"][1, 1] ≈ 0.1
                @test ismissing(ds["u"][3, 1])
            end
        end
    end

    @testset "time-weighted DAC referencing" begin
        # Glider dwells shallow, then sweeps deep: DAC (a time average along the track)
        # differs from the depth average — referencing must weight by residence time.
        u_o(z) = 0.3 * exp(-z / 60)                  # smooth, surface-intensified
        nt = 1000
        t0 = DateTime(2022, 11, 4)
        times = t0 .+ Second.(10 .* (0:nt-1))
        tunix = datetime2unix.(times)
        nd = round(Int, 0.6nt)
        depth = [i <= nd ? 20.0 + 5 * sin(i / 9) :
                 5 + 190 * (i - nd) / (nt - nd) for i in 1:nt]
        offsets = collect(0.0:1.0:30.0)
        celldepth = offsets .+ depth'
        ug = 0.2
        E = [u_o(celldepth[k, i]) - ug for k in 1:31, i in 1:nt]
        pp = ProcessedPings(times, tunix, depth, fill(90.0, nt), offsets, E,
            zeros(31, nt), zeros(31, nt), celldepth, :down, fill((1, 2, 4), nt))
        dacu = mean(u_o.(depth))                     # what the nav DAC measures
        dacdf = DataFrame(yo=[1], t_start=[times[1]], t_end=[times[end]],
            t_mid=[times[nt÷2]], u=[dacu], v=[0.0])

        sh_tw = solve_shear(pp, dacdf; opts=ShearOptions(referencing=:timeweighted))
        sh_si = solve_shear(pp, dacdf; opts=ShearOptions(referencing=:simple))
        g = (sh_tw.nobs .>= 4) .&& (sh_tw.z .< 200)
        err_tw = maximum(abs.(sh_tw.u[g] .- u_o.(sh_tw.z[g])))
        err_si = maximum(abs.(sh_si.u[g] .- u_o.(sh_si.z[g])))
        @test err_tw < 0.04
        @test err_si > err_tw + 0.02                 # simple referencing is biased here

        inv_tw = solve_inverse(pp, dacdf; opts=InverseOptions(dac_form=:ocean_timeweighted))
        inv_pl = solve_inverse(pp, dacdf)            # plain :ocean form
        gi = inv_tw.nobs .> 20
        err_itw = maximum(abs.(inv_tw.u[gi] .- u_o.(inv_tw.z[gi])))
        err_ipl = maximum(abs.(inv_pl.u[inv_pl.nobs .> 20] .-
                               u_o.(inv_pl.z[inv_pl.nobs .> 20])))
        @test err_itw < 0.04
        @test err_ipl > err_itw                      # plain depth-mean form is biased here
    end

    @testset "slocum helpers" begin
        t0 = DateTime(2023, 5, 1)
        df = DataFrame(
            time=t0 .+ Minute.(0:9),
            source_file=[fill("seg1", 5); fill("seg2", 5)],
            m_water_vx=[missing, missing, 0.1, missing, 0.1,
                        missing, 0.05, missing, missing, 0.05],
            m_water_vy=[missing, missing, 0.0, missing, 0.0,
                        missing, 0.0, missing, missing, 0.0],
            m_gps_mag_var=fill(deg2rad(90.0), 10),   # extreme declination for clarity
            depth=fill(50.0, 10),
            latitude=fill(18.0, 10), longitude=fill(-64.0, 10),
            m_heading=fill(π / 2, 10), m_pitch=fill(-0.3, 10), m_roll=zeros(10))
        dacs = dac_from_slocum(df)
        @test nrow(dacs) == 2
        # (vx=0.1, vy=0) rotated CCW by +90°: u→0, v→0.1
        @test dacs.u[1] ≈ 0.0 atol = 1e-12
        @test dacs.v[1] ≈ 0.1 atol = 1e-12
        nav = slocum_nav(df)
        @test length(nav) == 10
        @test nav.heading[1] ≈ 90.0
        @test nav.declination[1] ≈ 90.0
        @test nav.lat[1] == 18.0
    end

    if isdir(M38_PLD) && isfile(M38_NC)
        @testset "M38 acceptance: \$PNOR real-time stream" begin
            files = seaexplorer_files(M38_PLD, "ad2cp.raw")[10:11]
            a = load_pnor(files)
            @test ncells(a) == 15
            @test a.config.coordsystem === :beam
            @test a.config.cellsize ≈ 2.0 atol = 1e-6
            @test a.config.blanking ≈ 0.7 atol = 1e-6
            @test a.config.serial == 102381
            @test a.range ≈ collect(2.7:2.0:30.7) atol = 1e-3
            # cross-check against the full-resolution netCDF (match by time)
            nc = load_ad2cp(M38_NC)
            i1 = searchsortedfirst(nc.t, a.t[1] - 2)
            i2 = searchsortedlast(nc.t, a.t[end] + 2)
            hs = Float64[]; hn = Float64[]; vs = Float64[]; vn = Float64[]
            for i in i1:i2
                j = searchsortedfirst(a.t, nc.t[i] - 0.6)
                (1 <= j <= length(a)) || continue
                abs(a.t[j] - nc.t[i]) < 1.0 || continue
                if isfinite(a.heading[j]) && isfinite(nc.heading[i])
                    push!(hs, a.heading[j]); push!(hn, nc.heading[i])
                end
                if isfinite(a.vel[5, 1, j]) && isfinite(nc.vel[5, 1, i])
                    push!(vs, a.vel[5, 1, j]); push!(vn, nc.vel[5, 1, i])
                end
            end
            @test length(hs) > 200
            @test cor(hs, hn) > 0.999
            @test length(vs) > 100
            @test cor(vs, vn) > 0.98
            @test median(abs.(vs .- vn)) < 0.02       # 0.01 m/s stream quantization
            @info "PNOR vs netCDF: $(length(hs)) matched ensembles, " *
                  "r_heading=$(round(cor(hs, hn), digits=4)), r_vel=$(round(cor(vs, vn), digits=3))"
        end
    end

    @testset "native .ad2cp reader: synthetic binary" begin
        strrec, avgrec = synthetic_ad2cp_records()
        mktempdir() do d
            f = joinpath(d, "syn.ad2cp")
            # junk bytes between records exercise the resynchronization path
            write(f, vcat(strrec, UInt8[0xFF, 0xA5, 0x00], avgrec))
            a = @test_logs (:warn, r"resynchronized") read_ad2cp(f)
            @test length(a) == 1
            @test ncells(a) == 2
            @test a.config.serial == 42
            @test a.config.cellsize ≈ 1.0
            @test a.config.blanking ≈ 0.5
            @test a.config.coordsystem === :beam
            @test a.config.declination ≈ 3.0
            @test a.config.salinity_setting ≈ 35.0
            @test a.config.beam_theta[1] == 47.5 && a.config.beam_phi[2] == -90.0
            @test a.time[1] == DateTime(2022, 11, 4, 6, 0, 0, 438)
            @test a.vel[1, 1, 1] ≈ 0.1f0 && a.vel[2, 1, 1] ≈ 0.2f0
            @test a.vel[1, 2, 1] ≈ -0.3f0 && a.vel[2, 4, 1] ≈ -0.8f0
            @test a.amp[1, 1, 1] ≈ 70.0f0                # counts × 0.5 dB
            @test a.corr[2, 4, 1] ≈ 97.0f0
            @test a.heading[1] ≈ 90.0f0 && a.pitch[1] ≈ -17.5f0
            @test a.pressure[1] ≈ 12.345
            @test a.ensemble[1] == 7.0
            # corrupt the data checksum → record is dropped
            buf = read(f)
            buf[end] ⊻= 0xFF
            f2 = joinpath(d, "bad.ad2cp")
            write(f2, buf)
            @test_throws Exception read_ad2cp(f2)       # no average records survive
        end
    end

    M38_BIN = joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp")
    if isfile(M38_BIN) && isfile(M38_NC)
        @testset "M38 acceptance: native binary ≡ MIDAS netCDF" begin
            a = read_ad2cp(M38_BIN)
            b = load_ad2cp(M38_NC)
            @test length(a) == length(b) == 124_752
            @test a.time == b.time
            same32(x, y) = (d = filter(isfinite, vec(x .- y)); isempty(d) ? 0.0 : maximum(abs.(d)))
            @test same32(a.vel, b.vel) == 0
            @test same32(a.amp, b.amp) == 0
            @test same32(a.corr, b.corr) == 0
            @test same32(a.heading, b.heading) == 0
            @test same32(a.pitch, b.pitch) == 0
            @test same32(a.roll, b.roll) == 0
            @test same32(a.soundspeed, b.soundspeed) == 0
            @test same32(a.temperature, b.temperature) == 0
            @test same32(a.accel, b.accel) == 0
            @test same32(a.mag, b.mag) == 0
            @test maximum(abs.(a.pressure .- b.pressure)) < 1e-4   # MIDAS float32 rounding
            @test a.ensemble == b.ensemble
            @test a.config.cellsize == b.config.cellsize
            @test a.config.blanking ≈ b.config.blanking atol = 1e-6
            @test a.config.beam2xyz[1, 1] ≈ b.config.beam2xyz[1, 1] atol = 1e-6
            @test a.range ≈ b.range atol = 1e-3
            @test length(a.bt) == length(b.bt) == 124_751
            @test same32(a.bt.vel, b.bt.vel) == 0
            @test same32(a.bt.distance, b.bt.distance) == 0
            @test same32(a.bt.fom, b.bt.fom) == 0
            # dispatch through load_ad2cp
            a2 = load_ad2cp(M38_BIN)
            @test length(a2) == length(a)
            @info "native .ad2cp reader ≡ MIDAS netCDF: bit-identical on M38 " *
                  "($(length(a)) ensembles + $(length(a.bt)) BT records)"
        end
    end

    @testset "shear-bias calibration (synthetic injection)" begin
        u_o(z) = 0.15 * sin(2π * z / 160) + 0.05
        nyo, npy = 3, 400
        nt = nyo * npy
        headings = [30.0, 150.0, 270.0]
        t0 = DateTime(2022, 11, 4)
        times = t0 .+ Second.(10 .* (0:nt-1))
        tunix = datetime2unix.(times)
        offsets = collect(0.0:1.0:30.0)
        depth = Float64[]; hdg = Float64[]
        for y in 1:nyo, i in 1:npy
            push!(depth, 10 + 180 * (1 - abs(1 - 2i / npy)))
            push!(hdg, headings[y])
        end
        celldepth = offsets .+ depth'
        strue = -4e-4                                # injected track-frame bias slope
        B = strue .* (offsets .- mean(offsets))
        E = Matrix{Float64}(undef, 31, nt); N = similar(E)
        for i in 1:nt
            h = hdg[i]
            ugE, ugN = 0.3 * sind(h), 0.3 * cosd(h)  # glider moves along heading
            for k in 1:31
                E[k, i] = u_o(celldepth[k, i]) - ugE + B[k] * sind(h)
                N[k, i] = -ugN + B[k] * cosd(h)
            end
        end
        pp = ProcessedPings(times, tunix, depth, hdg, offsets, E, N, zeros(31, nt),
            celldepth, :down, fill((1, 2, 4), nt))
        pingmean0 = [mean(pp.E[:, i]) for i in 1:nt]

        b = shear_bias(pp; min_count=100)
        @test b.slope_along ≈ strue rtol = 0.15
        @test abs(b.slope_cross) < 5e-5
        @test b.heading_concentration < 0.5

        dacdf = DataFrame(yo=1:nyo,
            t_start=[times[(y-1)*npy+1] for y in 1:nyo],
            t_end=[times[y*npy] for y in 1:nyo],
            t_mid=[times[(y-1)*npy+npy÷2] for y in 1:nyo],
            u=[mean(u_o.(depth[(y-1)*npy+1:y*npy])) * 1.0 for y in 1:nyo],
            v=zeros(nyo))
        # DAC here is the time-mean of u_o along the track (east component only)
        sh0 = solve_shear(pp, dacdf)
        g0 = sh0.nobs .>= 4
        err0 = maximum(abs.(sh0.u[g0] .- u_o.(sh0.z[g0])))

        slopes = calibrate_shear_bias!(pp; passes=1, min_count=100)
        @test abs(slopes[end]) < 1e-8                # exact in one pass at full coverage
        pingmean1 = [mean(pp.E[:, i]) for i in 1:nt]
        @test maximum(abs.(pingmean1 .- pingmean0)) < 1e-12   # inverse content untouched

        sh1 = solve_shear(pp, dacdf)
        g1 = sh1.nobs .>= 4
        err1 = maximum(abs.(sh1.u[g1] .- u_o.(sh1.z[g1])))
        @test err1 < 0.05
        @test err1 < 0.6 * err0                      # injected tilt removed
    end

    if isfile(joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp")) && isdir(M38_NAV)
        @testset "M38 acceptance: shear-bias calibration" begin
            a0 = read_ad2cp(joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
            a = a0[1:2:length(a0)]                   # subsample for test runtime
            nav = load_seaexplorer_nav(M38_NAV)
            qc!(a)
            decl = magnetic_declination(nav, a.t)
            p = process_pings(a; lat=69.5, declination=decl)
            b = shear_bias(p)
            @test -7e-4 < b.slope_along < -2e-4      # the documented M38 bias
            @test abs(b.slope_cross) < 1e-4
            @test b.heading_concentration < 0.6
            slopes = calibrate_shear_bias!(p; passes=1)
            @test abs(slopes[end]) < 1e-6            # removed exactly
            # residual pairwise bias bounded at every depth band
            nk = length(p.offsets)
            for (z1, z2) in ((0, 200), (200, 500), (500, 1000))
                s = 0.0; n = 0
                for i in 1:length(p)
                    h = p.heading[i]
                    isfinite(h) || continue
                    sh, ch = sind(h), cosd(h)
                    for k in 1:nk-1
                        (isfinite(p.E[k, i]) && isfinite(p.N[k, i]) &&
                         isfinite(p.E[k+1, i]) && isfinite(p.N[k+1, i])) || continue
                        z = (p.celldepth[k, i] + p.celldepth[k+1, i]) / 2
                        (isfinite(z) && z1 <= z < z2) || continue
                        s += (p.E[k+1, i] - p.E[k, i]) * sh + (p.N[k+1, i] - p.N[k, i]) * ch
                        n += 1
                    end
                end
                n > 10_000 && @test abs(s / n) < 2e-4
            end
            @info "M38 shear-bias: slope $(round(b.slope_along, sigdigits=3)) s⁻¹ " *
                  "(cross $(round(b.slope_cross, sigdigits=2))), residual " *
                  "$(round(slopes[end], sigdigits=2))"
        end
    end

    @testset "vertical velocity product" begin
        # U_rel = w_ocean − w_glider ⇒ vertical_velocity recovers w_ocean exactly
        w_o(z) = 0.01 * sin(2π * z / 120)
        nt = 400
        t0 = DateTime(2022, 11, 4)
        times = t0 .+ Second.(10 .* (0:nt-1))
        tunix = datetime2unix.(times)
        depth = [10 + 180 * (1 - abs(1 - 2i / nt)) for i in 1:nt]
        wg = fill(NaN, nt)
        for i in 2:nt-1
            wg[i] = -(depth[i+1] - depth[i-1]) / (tunix[i+1] - tunix[i-1])
        end
        offsets = collect(0.0:1.0:30.0)
        celldepth = offsets .+ depth'
        U = Matrix{Float64}(undef, 31, nt)
        for i in 1:nt, k in 1:31
            U[k, i] = w_o(celldepth[k, i]) - (isfinite(wg[i]) ? wg[i] : 0.0)
        end
        pp = ProcessedPings(times, tunix, depth, fill(90.0, nt), offsets,
            zeros(31, nt), zeros(31, nt), U, celldepth, :down, fill((1, 2, 4), nt))
        w = vertical_velocity(pp)
        err = [abs(w[k, i] - w_o(celldepth[k, i])) for k in 1:31, i in 2:nt-1
               if isfinite(w[k, i])]
        @test length(err) > 5000
        @test maximum(err) < 1e-10

        # solve_w: both methods recover the prescribed w_o(z)
        segs = DataFrame(yo=[1], t_start=[times[1]], t_end=[times[end]],
            t_mid=[times[nt÷2]], u=[0.0], v=[0.0])
        for method in (:direct, :inverse)
            wd = solve_w(pp, segs; method)
            g = (wd.nobs .> 20) .&& isfinite.(wd.w)
            @test count(g) > 10
            @test maximum(abs.(wd.w[g] .- w_o.(wd.z[g]))) < 3e-3
        end
    end

    if isfile(joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
        @testset "M38 acceptance: w product + compass field check" begin
            a = read_ad2cp(joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
            tbl, ptp = compass_field_check(a)
            @test nrow(tbl) >= 8
            @test 0 <= ptp < 0.5
            @info "M38 compass |B| heading variation: $(round(100ptp, digits=1))% peak-to-peak"
            qc!(a)
            p = process_pings(a[1:5:length(a)]; lat=69.5)
            w = vertical_velocity(p)
            fw = filter(isfinite, vec(w))
            @test length(fw) > 50_000
            @test abs(median(fw)) < 0.02              # ocean w ≈ 0 in the median
            @test quantile(abs.(fw), 0.9) < 0.15
            @info "M38 w_water: median $(round(median(fw), digits=4)) m/s, " *
                  "p90|w| $(round(quantile(abs.(fw), 0.9), digits=3))"
        end
    end

    @testset "cell_quality diagnostic" begin
        mktempdir() do d
            f = write_synthetic_midas(joinpath(d, "syn.ad2cp.00000.nc"); nt=8, nc=3)
            a = load_ad2cp(f)
            a.corr[3, :, :] .= 20.0f0            # far cell decorrelated
            q = cell_quality(a; thr=QCThresholds(snr_db=NaN))
            @test nrow(q) == 12                  # 3 cells × 4 beams
            @test all(q.keep_corr[q.cell .== 3] .== 0.0)
            @test all(q.keep_corr[q.cell .== 1] .== 1.0)
            @test all(q.med_corr[q.cell .== 3] .== 20.0)
        end
    end

    if isfile(joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
        @testset "M38 acceptance: per-cell quality structure" begin
            a = read_ad2cp(joinpath(M38_DIR, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
            q = cell_quality(a)
            agg = combine(groupby(q, :cell), :med_corr => mean => :corr,
                :keep_all => mean => :keep)
            sort!(agg, :cell)
            @test agg.corr[2] > 90               # near cells: high correlation
            @test agg.corr[15] < 15              # far cells decorrelate (clear water)
            @test issorted(agg.corr[5:15]; rev=true)
            @test agg.keep[2] > 0.85
            @test agg.keep[15] < 0.15
            @info "M38 cell quality: med corr cell2=$(round(agg.corr[2], digits=1)) " *
                  "→ cell8=$(round(agg.corr[8], digits=1)) → cell15=$(round(agg.corr[15], digits=1)); " *
                  "effective range ≈ cells 2–8 (≤ ~17 m)"
        end
    end

    @testset "bt_valid false-lock screens" begin
        t0 = DateTime(2022, 11, 4)
        n = 6
        times = t0 .+ Second.(600 .* (0:n-1))
        mk(dist, press) = BottomTrackData(times, datetime2unix.(times),
            fill(0.1f0, 4, n), repeat(dist', 4), fill(100.0f0, 4, n),
            press, fill(90.0f0, n), fill(-17.0f0, n), fill(0.0f0, n),
            fill(1490.0f0, n))
        # near-field target (1.5 m) at 100 m depth, while the platform later reaches
        # 500 m nearby: rejected by BOTH min_range and the bathymetry check
        bt1 = mk(fill(1.5, n), [100.0, 100, 100, 500, 100, 100])
        @test count(bt_valid(bt1)) == 0
        # bathymetry check alone: rejects the five 100-m locks, but a false lock AT the
        # deepest record has nothing deeper nearby to disprove it — min_range covers it
        @test count(bt_valid(bt1; min_range=0.0)) == 4
        @test count(bt_valid(bt1; min_range=0.0, bathymetry_check=false)) == 4n
        # genuine approach: 15-m range at 480 m, deepest nearby record 490 m — kept
        bt2 = mk(fill(15.0, n), [480.0, 485, 490, 488, 483, 480])
        @test count(bt_valid(bt2)) == 4n
        # range gate still applies
        @test count(bt_valid(bt2; max_range=10.0)) == 0
    end

    @testset "telemetered pld1.sub AD2CP subset (load_pld_adcp)" begin
        hdr = "PLD_REALTIMECLOCK;NAV_LONGITUDE;AD2CP_TIME;AD2CP_HEADING;AD2CP_PITCH;" *
              "AD2CP_ROLL;AD2CP_PRESSURE;AD2CP_V1_CN1;AD2CP_V2_CN1;AD2CP_V3_CN1;AD2CP_V4_CN1;" *
              "AD2CP_V1_CN2;AD2CP_V2_CN2;AD2CP_V3_CN2;AD2CP_V4_CN2\n"
        # AD2CP_TIME is MMDDYY; heading 9999 marks instrument-off rows
        rows_seg = "03/11/2022 12:00:31.100;126.5;110322 12:00:00;90.0;-17.5;1.0;55.5;" *
                   "0.11;-0.02;0.30;0.04;0.12;-0.03;0.31;0.05\n" *
                   "03/11/2022 12:01:01.100;126.5;110322 12:00:30;9999.0;9999.0;9999.0;9999.0;" *
                   "9999.0;9999.0;9999.0;9999.0;9999.0;9999.0;9999.0;9999.0\n"
        # GLIMPSE export: YO_NUMBER column first; duplicates the first segment ping
        # (different values — segment-file source listed first must win) + adds one
        rows_gl = "YO_NUMBER;" * hdr *
                  "5;03/11/2022 12:00:31.100;126.5;110322 12:00:00;90.0;-17.5;1.0;55.5;" *
                  "0.99;0.99;0.99;0.99;0.99;0.99;0.99;0.99\n" *
                  "5;03/11/2022 12:01:31.200;126.5;110322 12:01:00;91.0;-17.4;0.9;56.5;" *
                  "-0.21;0.07;0.15;0.02;-0.22;0.08;0.16;0.03\n"
        mktempdir() do del
        mktempdir() do gl
            open(joinpath(del, "sea064.38.pld1.sub.1.gz"), "w") do io
                gz = GzipCompressorStream(io)
                write(gz, hdr * rows_seg)
                close(gz)
            end
            write(joinpath(gl, "SEA064.38.pld1.sub.all.csv"), rows_gl)
            a = load_pld_adcp([del, gl]; stream="pld1.sub", cellsize=2.0, blanking=0.7,
                serial=42)
            @test length(a) == 2                       # off-row skipped, duplicate deduped
            @test ncells(a) == 2
            @test a.time[1] == DateTime(2022, 11, 3, 12, 0, 0)   # MMDDYY parsed
            @test a.time[2] == DateTime(2022, 11, 3, 12, 1, 0)
            @test a.vel[1, 1, 1] ≈ 0.11f0              # segment source won the duplicate
            @test a.vel[2, 3, 2] ≈ 0.16f0
            @test a.heading[2] ≈ 91.0f0
            @test a.pressure[1] ≈ 55.5
            @test all(isnan, a.amp) && all(isnan, a.corr)
            @test a.config.coordsystem === :beam
            @test a.range ≈ [2.7, 4.7]
            # QC runs with the missing screens as no-ops; large-blanking default
            # keeps cell 1, and first_cells=1 still masks it on request
            b = deepcopy(a)
            qc!(a)
            @test isfinite(a.vel[1, 1, 1])
            @test isfinite(a.vel[2, 1, 2])
            qc!(b; thr=QCThresholds(first_cells=1))
            @test all(isnan, b.vel[1, :, :])
            # soundspeed vector length mismatch is an error
            @test_throws ErrorException load_pld_adcp([del, gl]; stream="pld1.sub",
                cellsize=2.0, blanking=0.7, soundspeed=[1500.0])
            # onboard sound speed reconstructed from configured salinity + CTD T,
            # then the standard correction chain applies
            c = load_pld_adcp([del, gl]; stream="pld1.sub", cellsize=2.0, blanking=0.7)
            ctd_t = c.t .+ [-60.0, 60.0]
            cs = onboard_soundspeed!(c, ctd_t, [4.0, 4.2]; salinity=38.0, lat=69.5)
            @test cs === c.soundspeed
            @test all(v -> 1440 < v < 1520, c.soundspeed)
            c_true = soundspeed_from_ctd.(35.0, [4.0, 4.2], [55.0, 56.0], 0.0, 69.5)
            scale = soundspeed_correction(c, ctd_t, c_true)
            @test all(isfinite, scale)
            @test all(s -> 0.99 < s < 1.0, scale)     # S 38→35 shrinks c ⇒ scale < 1
        end
        end
    end

    @testset "robustness: corrupt and missing inputs" begin
        # --- SeaExplorer: missing + corrupt segments ---
        mktempdir() do d
            write_synthetic_gli(joinpath(d, "sea064.38.gli.sub.1.gz"); nrows=4)
            write_synthetic_gli(joinpath(d, "sea064.38.gli.sub.2.gz"); nrows=3)
            write_synthetic_gli(joinpath(d, "sea064.38.gli.sub.5.gz"); nrows=3)
            write(joinpath(d, "sea064.38.gli.sub.3.gz"), UInt8[0x1f, 0x8b, 0xde, 0xad])
            @test missing_segments(d, "gli.sub") == [4]
            nav = @test_logs (:warn, r"missing segment") (:warn, r"unreadable log file") match_mode = :any load_seaexplorer_nav(d)
            @test length(nav) == 10                      # rows from the 3 good segments
        end
        # --- load_ad2cp: corrupt netCDF alongside a good one ---
        mktempdir() do d
            write_synthetic_midas(joinpath(d, "syn.ad2cp.00000.nc"))
            write(joinpath(d, "syn.ad2cp.00001.nc"), UInt8.(mod.(1:256, 251)))
            a = @test_logs (:warn, r"skipping unreadable") match_mode = :any load_ad2cp(d)
            @test length(a) == 6
        end
        # --- native binary: truncated final record ---
        strrec, avgrec = synthetic_ad2cp_records()
        mktempdir() do d
            f = joinpath(d, "trunc.ad2cp")
            write(f, vcat(strrec, avgrec, avgrec[1:end-25]))
            a = @test_logs (:warn, r"mid-record") match_mode = :any read_ad2cp(f)
            @test length(a) == 1
        end
        # --- load_pnor: corrupt gz skipped, good plain file parsed ---
        mktempdir() do d
            lines = ["\$PNORI,4,Glider42,4,2,0.70,2.00,2*00",
                "\$PNORS,110422,062254,0,0,11.5,1485.1,114.1,19.8,8.3,4.765,7.68,0,0*00",
                "\$PNORC,110422,062254,1,-0.24,-0.01,0.24,-0.01,,,C,133,126,151,116,100,100,100,99*00",
                "\$PNORC,110422,062254,2,-0.20,0.0,0.2,0.0,,,C,130,120,150,110,95,96,97,98*00"]
            write(joinpath(d, "sea064.38.ad2cp.raw.1"), join(lines, "\n"))
            write(joinpath(d, "sea064.38.ad2cp.raw.2.gz"), UInt8[0x1f, 0x8b, 0x00])
            a = @test_logs (:warn, r"skipping unreadable") match_mode = :any load_pnor(d;
                validate_checksum=false)
            @test length(a) == 1
            @test ncells(a) == 2
            @test a.vel[1, 1, 1] ≈ -0.24f0
        end
        # --- solvers: empty and starved segment tables never throw ---
        t0 = DateTime(2022, 11, 4)
        nt = 50
        times = t0 .+ Second.(10 .* (0:nt-1))
        offsets = collect(0.0:1.0:10.0)
        depth = fill(50.0, nt)
        pp = ProcessedPings(times, datetime2unix.(times), depth, fill(90.0, nt), offsets,
            fill(0.1, 11, nt), fill(0.0, 11, nt), zeros(11, nt), offsets .+ depth',
            :down, fill((1, 2, 4), nt))
        emptydac = DataFrame(yo=Int[], t_start=DateTime[], t_end=DateTime[],
            t_mid=DateTime[], u=Float64[], v=Float64[])
        @test isempty(solve_inverse(pp, emptydac))
        @test isempty(solve_shear(pp, emptydac))
        starved = DataFrame(yo=[1], t_start=[times[1]], t_end=[times[3]],
            t_mid=[times[2]], u=[0.0], v=[0.0])          # 3 pings < min_pings
        out = @test_logs (:info, r"solved 0 of 1") solve_inverse(pp, starved)
        @test isempty(out)
        # --- declination: constant extrapolation outside nav coverage ---
        navt = t0 .+ Hour.(0:2)
        nav = GliderNav(navt, datetime2unix.(navt), fill(5.0, 3), fill(69.5, 3),
            fill(90.0, 3), zeros(3), zeros(3), zeros(3), fill(10.0, 3),
            fill(Int16(110), 3), fill(Int8(1), 3), fill(-1.0, 3), DataFrame())
        tq = datetime2unix.(t0 .+ Hour.(-2:6))
        decl = @test_logs (:warn, r"extrapolated") match_mode = :any magnetic_declination(nav, tq)
        @test all(isfinite, decl)
        @test decl[1] == decl[2]                         # leading constant fill
        @test decl[end] == decl[end-1]
        # --- data_gaps / coverage ---
        g = data_gaps([0.0, 10, 20, 1000, 1010])
        @test nrow(g) == 1
        @test g.duration[1] == 980
        mktempdir() do d
            f = write_synthetic_midas(joinpath(d, "syn.ad2cp.00000.nc"))
            a = load_ad2cp(f)
            cov = coverage(a)
            @test cov.n == 6
            @test cov.median_dt == 10.0
            @test isempty(cov.gaps)
            @test cov.finite_vel == 1.0
            @test cov.n_bt == 5
        end
    end

    if isfile(M38_NC) && isdir(M38_NAV)
        @testset "M38 acceptance: coverage reporting" begin
            a = load_ad2cp(M38_NC)
            cov = coverage(a)
            @test cov.n == 124_752
            @test nrow(cov.gaps) > 10                    # duty-cycled mission
            @test cov.gap_total > 50 * 86400             # ADCP off most of Dec-Feb
            @test isempty(missing_segments(M38_NAV, "gli.sub"))
            @info "M38 coverage: $(nrow(cov.gaps)) gaps totalling " *
                  "$(round(cov.gap_total / 86400, digits=1)) days; " *
                  "finite vel fraction $(round(cov.finite_vel, digits=2))"
        end
    end

    if isfile(M38_NC) && isdir(M38_PLD) && isdir(M38_NAV)
        @testset "M38 acceptance: real-time vs delayed pipeline (Task 5)" begin
            # identical pipeline on the $PNOR telemetry stream vs the MIDAS netCDF,
            # restricted to the first 3 days to bound runtime (full mission:
            # examples/realtime_onboard.jl m38 — inverse r=0.9996, rms 4.6 mm/s)
            a_d = load_ad2cp(M38_NC)
            a_r = load_pnor(M38_PLD)
            tcut = a_r.t[1] + 3 * 86400
            a_d = a_d[a_d.t .< tcut]
            a_r = a_r[a_r.t .< tcut]
            nav = load_seaexplorer_nav(M38_NAV)
            dac = compute_dac(nav)
            prods = DataFrame[]
            for (a, look) in ((a_d, :auto), (a_r, :down))
                qc!(a)
                p = process_pings(a; lat=69.5, look=look,
                    declination=magnetic_declination(nav, a.t))
                calibrate_shear_bias!(p)
                push!(prods, solve_inverse(p, dac))
            end
            j = innerjoin(prods[1], prods[2]; on=[:yo, :z], makeunique=true)
            m = (j.nobs .> 10) .&& (j.nobs_1 .> 10) .&& isfinite.(j.u) .&& isfinite.(j.u_1)
            @test length(unique(j.yo[m])) >= 10
            @test length(unique(prods[1].yo)) == length(unique(prods[2].yo))  # no yos lost
            for c in (:u, :v)
                d = j[m, c] .- j[m, Symbol(c, :_1)]
                @test cor(j[m, c], j[m, Symbol(c, :_1)]) > 0.995
                @test sqrt(mean(d .^ 2)) < 0.010          # quantization-level rms
                @test abs(mean(d)) < 0.002                # no systematic offset
            end
        end
    end

    if isfile(M38_NC) && isdir(M38_PLD) && isdir(M38_NAV) && isdir(joinpath(M38_DIR, "glimpse"))
        @testset "M38 acceptance: telemetered pld1.sub route" begin
            # multi-route load (glider-computer segments + GLIMPSE export, deduped)
            tele = load_pld_adcp([M38_PLD, joinpath(M38_DIR, "glimpse")];
                stream="38.pld1.sub", cellsize=2.0, blanking=0.7, serial=102381)
            @test length(tele) > 40_000
            @test ncells(tele) == 6
            @test issorted(tele.t)
            # each telemetered ping is a single subsampled ensemble: match one
            # against the netCDF at the same instrument timestamp, quantized 0.01
            a = load_ad2cp(M38_NC)
            k = findfirst(t -> abs(t - tele.t[100]) < 0.5, a.t)
            @test k !== nothing
            vq = round.(Float64.(a.vel[1:6, :, k]) .* 100) ./ 100
            m = isfinite.(tele.vel[:, :, 100])
            @test maximum(abs, Float64.(tele.vel[:, :, 100])[m] .- vq[m]) < 0.011
            # 3-day mini-pipeline: telemetered inverse tracks the delayed inverse
            tcut = tele.t[1] + 3 * 86400
            tele3 = tele[tele.t .< tcut]
            a3 = a[a.t .< tcut]
            nav = load_seaexplorer_nav(M38_NAV)
            dac = compute_dac(nav)
            qc!(tele3); qc!(a3)
            p_t = process_pings(tele3; lat=69.5, look=:down,
                declination=magnetic_declination(nav, tele3.t))
            p_d = process_pings(a3; lat=69.5,
                declination=magnetic_declination(nav, a3.t))
            calibrate_shear_bias!(p_t); calibrate_shear_bias!(p_d)
            inv_t = solve_inverse(p_t, dac)
            inv_d = solve_inverse(p_d, dac)
            j = innerjoin(inv_d, inv_t; on=[:yo, :z], makeunique=true)
            gm = (j.nobs .> 10) .&& (j.nobs_1 .> 3) .&& isfinite.(j.u) .&& isfinite.(j.u_1)
            @test length(unique(j.yo[gm])) >= 10
            for c in (:u, :v)
                d = j[gm, c] .- j[gm, Symbol(c, :_1)]
                @test cor(j[gm, c], j[gm, Symbol(c, :_1)]) > 0.9
                @test sqrt(mean(d .^ 2)) < 0.06
                @test abs(mean(d)) < 0.01
            end
        end
    end

    @testset "Aqua quality assurance" begin
        import Aqua
        Aqua.test_all(GliderADCP; ambiguities=false)
    end

    # ---------------- acceptance tests on local reference data ----------------

    if isfile(M38_NC)
        @testset "M38 acceptance: AD2CP netCDF" begin
            a = load_ad2cp(M38_NC)
            @test length(a) == 124_752
            @test ncells(a) == 15
            @test a.range ≈ collect(2.7:2.0:30.7) atol = 1e-3
            @test a.config.serial == 102381
            @test a.config.frequency == 1000
            @test a.config.cellsize ≈ 2.0 atol = 1e-6
            @test a.config.blanking ≈ 0.7 atol = 1e-6
            @test a.config.coordsystem === :beam
            @test a.config.velocity_range ≈ 2.5 atol = 1e-6
            @test a.config.beam_theta == (47.5, 25.0, 47.5, 25.0)
            @test a.config.beam_phi == (0.0, -90.0, 180.0, 90.0)
            @test a.config.beam2xyz[1, 1] ≈ 0.6782 atol = 1e-4
            @test a.config.salinity_setting == 38.0
            @test issorted(a.time)
            @test a.bt !== nothing && length(a.bt) == 124_751
            # physical sanity
            @test maximum(filter(isfinite, a.pressure)) > 900       # deep mission
            @test count(!isnan, a.vel) / length(a.vel) > 0.9
            @test median(filter(isfinite, a.pitch)) < -10            # mostly diving records
            btd = filter(isfinite, vec(a.bt.distance))
            @test 0.1 < count(>(0), btd) / length(btd) < 0.5         # ~24% bottom lock
        end
    else
        @info "M38 reference netCDF not found — skipping M38 AD2CP acceptance tests"
    end

    if isfile(M38_NC)
        @testset "M38 smoke: QC + ENU transform (subset)" begin
            a0 = load_ad2cp(M38_NC)
            deep = findall(p -> isfinite(p) && p > 100, a0.pressure)
            a = a0[deep[1:5000]]
            @test detect_look_direction(a) === :down
            stats = qc!(a)
            @test 0 < stats.total < 0.7
            E, N, U, offs, used = enu_on_isobars(a)
            fE = filter(isfinite, E)
            fU = filter(isfinite, U)
            @test length(fE) > 10_000
            @test quantile(abs.(fE), 0.9) < 1.0      # relative horizontal speeds sane
            @test abs(median(fU)) < 0.4              # relative vertical ~ glider w magnitude
            @info "M38 ENU smoke: $(length(fE)) finite samples, " *
                  "QC rejected $(round(100 * stats.total, digits=1))%, " *
                  "med|E_rel|=$(round(median(abs.(fE)), digits=3)) m/s"
        end
    end

    if isdir(M38_NAV)
        @testset "M38 acceptance: SeaExplorer nav" begin
            nav = load_seaexplorer_nav(M38_NAV)
            @test length(nav) > 100_000
            @test issorted(nav.time)
            # Lofoten Basin mission
            gpsfix = findall(==(0), nav.deadreckoning)
            @test !isempty(gpsfix)
            @test all(67 .< nav.lat[gpsfix] .< 72)
            @test all(-2 .< nav.lon[gpsfix] .< 16)
            # surfacings: DR 1→0 transitions should be ≈ number of yos (191 segments)
            dr = nav.deadreckoning
            surfacings = count(i -> dr[i] == 0 && dr[i-1] == 1, 2:length(dr))
            @test 150 < surfacings < 400
            @info "M38 nav: $(length(nav)) records, $surfacings DR→GPS surfacing transitions"
        end
    else
        @info "M38 nav logs not found — skipping nav acceptance tests"
    end

    if isdir(M38_PLD)
        @testset "M38 acceptance: payload subset" begin
            files = seaexplorer_files(M38_PLD, "pld1.raw")[1:2]
            pld = load_seaexplorer_pld(files)
            @test nrow(pld) > 1000
            @test hasproperty(pld, :LEGATO_TEMPERATURE)
            @test hasproperty(pld, :AD2CP_HEADING)
            temps = collect(skipmissing(pld.LEGATO_TEMPERATURE))
            @test !isempty(temps) && all(-2 .< temps .< 20)
        end
    else
        @info "M38 payload logs not found — skipping payload acceptance tests"
    end

    if isfile(M48_NC)
        @testset "M48 acceptance: second MIDAS sample" begin
            a = load_ad2cp(M48_NC)
            @test length(a) > 10_000
            @test a.config.coordsystem === :beam
            @test issorted(a.time)
            @info "M48: $a"
        end
    else
        @info "M48 sample netCDF not found — skipping M48 acceptance tests"
    end

    if isdir(M38_NAV)
        @testset "M38 acceptance: DAC" begin
            nav = load_seaexplorer_nav(M38_NAV)
            dac = compute_dac(nav)
            @test 120 <= nrow(dac) <= 200
            spd = hypot.(dac.u, dac.v)
            @test median(spd) < 0.4                 # Lofoten Basin: moderate currents
            @test quantile(spd, 0.95) < 0.8
            @test all(dac.duration .> 600)
            drift = surface_drift(nav)
            @test nrow(drift) > 50
            @info "M38 DAC: $(nrow(dac)) segments, median |DAC| = " *
                  "$(round(median(spd), digits=3)) m/s; $(nrow(drift)) surface-drift pairs"
        end
    end

    if isfile(M38_NC) && isdir(M38_NAV)
        @testset "M38 acceptance: inverse vs bottom track + reference CSV" begin
            a0 = load_ad2cp(M38_NC)
            nav = load_seaexplorer_nav(M38_NAV)
            dac = compute_dac(nav)
            # M38 has NO genuine seafloor locks: its BT record is a persistent near-field
            # target (~1.7 m below the transducer) moving with the water. The hardened
            # default screens must reject essentially all of it:
            @test nrow(bt_velocity(a0; max_range=28.0)) < 100
            # the unscreened record remains useful as a WATER-FRAME consistency check of
            # the geometry chain (its "u_g" ≈ glider velocity relative to the water):
            btv = bt_velocity(a0; max_range=28.0, min_range=0.0, bathymetry_check=false)

            csvpath = joinpath(M38_DIR, "ad2cp/m38_processed/absolute_ocean_vel.csv")
            ref = isfile(csvpath) ? CSV.read(csvpath, DataFrame) : nothing

            # segments for the BT check: the six most bottom-track-rich yos
            nbt_in(row) = count(t -> datetime2unix(row.t_start) <= t <=
                                     datetime2unix(row.t_end), btv.t)
            counts = [nbt_in(row) for row in eachrow(dac)]
            segs_bt = sortperm(counts; rev=true)[1:6]
            @test minimum(counts[segs_bt]) > 50
            # segments for the CSV regression: yos whose window contains a reference
            # time_midpoint (their "midpoint" may be an end time — window matching)
            segs_csv = Int[]
            if ref !== nothing
                for (si, row) in enumerate(eachrow(dac))
                    t1 = datetime2unix(row.t_start) - 600
                    t2 = datetime2unix(row.t_end) + 600
                    any(t -> t1 <= t <= t2, ref.time_midpoint) && push!(segs_csv, si)
                end
                length(segs_csv) > 8 &&
                    (segs_csv = segs_csv[round.(Int, range(1, length(segs_csv); length=8))])
            end
            segs = sort(unique(vcat(segs_bt, segs_csv)))

            ug_pairs = Tuple{Float64,Float64}[]   # (inverse ug, bt u)
            vg_pairs = Tuple{Float64,Float64}[]
            oursol = DataFrame(yo=Int[], t_mid=DateTime[], z=Float64[], u=Float64[],
                v=Float64[], nobs=Int[])
            windows = Dict{Int,Tuple{Float64,Float64}}()
            for si in segs
                row = dac[si, :]
                t1, t2 = datetime2unix(row.t_start), datetime2unix(row.t_end)
                idx = findall(t -> t1 <= t <= t2, a0.t)
                length(idx) < 100 && continue
                sub = a0[idx]
                qc!(sub)
                pp = process_pings(sub; lat=69.0)
                gd = filter(isfinite, pp.depth)
                # DAC-only inverse (no bottom track!) → glider velocities
                sol = invert_segment(pp.E, pp.N, pp.celldepth, pp.t, maximum(gd);
                    dacu=row.u, dacv=row.v, opts=InverseOptions())
                sol === nothing && continue
                windows[row.yo] = (t1, t2)
                for k in eachindex(sol.z)
                    push!(oursol, (row.yo, row.t_mid, sol.z[k], sol.u[k], sol.v[k],
                        sol.nobs[k]))
                end
                si in segs_bt || continue
                # independent check: recovered glider velocity vs BT over-ground velocity
                for btrow in eachrow(btv[(btv.t .>= t1) .& (btv.t .<= t2), :])
                    j = argmin(abs.(sol.tping .- btrow.t))
                    abs(sol.tping[j] - btrow.t) <= 6 || continue
                    push!(ug_pairs, (sol.ug[j], btrow.u))
                    push!(vg_pairs, (sol.vg[j], btrow.v))
                end
            end
            @test length(ug_pairs) > 200
            ru = cor(first.(ug_pairs), last.(ug_pairs))
            rv = cor(first.(vg_pairs), last.(vg_pairs))
            mdu = median(abs.(first.(ug_pairs) .- last.(ug_pairs)))
            mdv = median(abs.(first.(vg_pairs) .- last.(vg_pairs)))
            @info "M38 inverse ug vs water-frame BT target (through-water consistency): " *
                  "n=$(length(ug_pairs)), r_u=$(round(ru, digits=3)), r_v=$(round(rv, digits=3)), " *
                  "med|Δu_g|=$(round(mdu, digits=3)), med|Δv_g|=$(round(mdv, digits=3)) m/s " *
                  "(offset ≈ the water velocity; see m38_validation.md Task 3)"
            @test ru > 0.6 && rv > 0.6
            @test mdu < 0.15 && mdv < 0.15

            # regression vs the prior Python (Slocum-AD2CP-style) processing
            if ref !== nothing
                rus = Float64[]; rvs = Float64[]
                for yo in unique(oursol.yo)
                    haskey(windows, yo) || continue
                    t1, t2 = windows[yo]
                    ours = oursol[(oursol.yo .== yo) .& (oursol.nobs .> 10), :]
                    isempty(ours) && continue
                    # their "time_midpoint" convention is uncertain (mid vs end) —
                    # match any reference yo whose stamp falls inside our fix-to-fix window
                    inwin = findall(t -> t1 - 600 <= t <= t2 + 600, ref.time_midpoint)
                    isempty(inwin) && continue
                    cnt = Dict{Float64,Int}()
                    for y in ref.yo_number[inwin]
                        cnt[y] = get(cnt, y, 0) + 1
                    end
                    ryo = argmax(cnt)
                    r = ref[ref.yo_number .== ryo, :]
                    zr = -r.depth_bins
                    ui = GliderADCP._interp1(ours.z, ours.u, zr)
                    vi = GliderADCP._interp1(ours.z, ours.v, zr)
                    ok = findall(k -> isfinite(ui[k]) && isfinite(vi[k]), eachindex(zr))
                    length(ok) < 20 && continue
                    push!(rus, cor(ui[ok], r.u_ocean_vel[ok]))
                    push!(rvs, cor(vi[ok], r.v_ocean_vel[ok]))
                end
                @info "M38 vs reference CSV: $(length(rus)) matched yos, " *
                      "median r_u=$(round(median(rus), digits=3)), " *
                      "median r_v=$(round(median(rvs), digits=3))"
                @test length(rus) >= 3
                @test median(rus) > 0.2 && median(rvs) > 0.2   # guards gross errors only
            else
                @info "reference absolute_ocean_vel.csv not found — skipping CSV regression"
            end
        end
    end

    GAD2CP_REF = joinpath(dirname(@__DIR__), "validation", "gliderad2cp_reference",
        "shear_out.nc")
    if isfile(GAD2CP_REF)
        @testset "gliderad2cp parity (SEA055 ground truth)" begin
            NCDataset(GAD2CP_REF) do ds
                g2(v) = coalesce.(Array(ds[v][:, :]), NaN)
                g1(v) = Float64.(coalesce.(Array(ds[v][:]), NaN))
                V = [g2("V$b") for b in 1:4]                 # isobar-regridded beams
                Xr, Yr, Zr = g2("X"), g2("Y"), g2("Z")
                Er, Nr, Ur = g2("E"), g2("N"), g2("U")
                H, P, R = g1("Heading"), g1("Pitch"), g1("Roll")

                # (a) transform parity: their beam synthesis makes the 4-beam set exactly
                # self-consistent, so our LS solve must reproduce their X,Y,Z and our
                # H·P rotation their E,N,U — to numerical precision.
                e = beam_unit_vectors((47.5, 25.0, 47.5, 25.0), (0.0, -90.0, 180.0, 90.0))
                E4 = permutedims(hcat(e...))
                S4 = (E4' * E4) \ E4'
                dxyz = 0.0; denu = 0.0; ncmp = 0
                for i in 1:7:size(V[1], 2)
                    (isfinite(H[i]) && isfinite(P[i]) && isfinite(R[i])) || continue
                    Rm = rotmat_xyz2enu(H[i], P[i], R[i])    # 'top' mount ⇒ F = I
                    for k in 1:size(V[1], 1)
                        b = (V[1][k, i], V[2][k, i], V[3][k, i], V[4][k, i])
                        (all(isfinite, b) && isfinite(Xr[k, i]) && isfinite(Er[k, i])) ||
                            continue
                        v = S4 * collect(b)
                        w = Rm * v
                        dxyz = max(dxyz, abs(v[1] - Xr[k, i]), abs(v[2] - Yr[k, i]),
                            abs(v[3] - Zr[k, i]))
                        denu = max(denu, abs(w[1] - Er[k, i]), abs(w[2] - Nr[k, i]),
                            abs(w[3] - Ur[k, i]))
                        ncmp += 1
                    end
                end
                @test ncmp > 5_000
                @test dxyz < 1e-12          # machine-exact beam solve
                @test denu < 1e-6           # Float32 attitude → trig rounding only
                @info "gliderad2cp transform parity: $ncmp samples, " *
                      "max|ΔXYZ| = $(round(dxyz, sigdigits=3)), " *
                      "max|ΔENU| = $(round(denu, sigdigits=3)) m/s"

                # (b) regrid parity: our isobaric regrid from their native-grid beams vs
                # their V2/V4 (never synthetic). Their per-beam depths use small-angle
                # approximations (exact only at roll=0), so agreement is close, not exact.
                time = DateTime.(coalesce.(Array(ds["time"][:]), DateTime(1970)))
                nt = length(H)
                cfg = AD2CPConfig(0, 1000.0, (47.5, 25.0, 47.5, 25.0),
                    (0.0, -90.0, 180.0, 90.0), fill(NaN, 4, 4), 30,
                    Float64(ds.attrib["avg_cellSize"]),
                    Float64(ds.attrib["avg_blankingDistance"]), :beam, 2.5, 0.0, 0.0, 0.0,
                    Dict{String,Any}())
                vel = Array{Float32}(undef, 30, 4, nt)
                for b in 1:4
                    vel[:, b, :] = Float32.(g2("VelocityBeam$b"))
                end
                a = AD2CPData(time, datetime2unix.(time),
                    Float64.(coalesce.(Array(ds["Velocity Range"][:]), NaN)), vel,
                    fill(NaN32, 30, 4, nt), fill(NaN32, 30, 4, nt),
                    Float32.(H), Float32.(P), Float32.(R),
                    g1("Pressure"), fill(NaN32, nt), Float32.(g1("SpeedOfSound")),
                    fill(NaN32, 3, nt), fill(NaN32, 3, nt),
                    zeros(nt), zeros(nt), zeros(nt), cfg, nothing)
                offsets = -g1("depth_offset")                # ours: signed positive down
                Vours, _, _ = regrid_beams(a; look=:up, offsets)
                # gliderad2cp assigns cell depths with small-angle formulas that are exact
                # only at roll = 0; we use the exact rotated beam geometry. Agreement is
                # therefore near-perfect at low roll and degrades as |roll| grows — assert
                # both the overall closeness and the sharp low-roll equivalence.
                for b in (2, 4)
                    Vb = Vours[:, b, :]
                    x = Float64[]; y = Float64[]; lowroll = Bool[]
                    for i in 1:nt, k in 1:size(Vb, 1)
                        if isfinite(Vb[k, i]) && isfinite(V[b][k, i])
                            push!(x, Vb[k, i]); push!(y, V[b][k, i])
                            push!(lowroll, isfinite(R[i]) && abs(R[i]) < 1)
                        end
                    end
                    d = abs.(x .- y)
                    @test length(x) > 100_000
                    @test median(d) < 0.005
                    @test mean(d .< 0.01) > 0.9
                    rl = cor(x[lowroll], y[lowroll])
                    @test rl > 0.99                      # exact-geometry equivalence at roll≈0
                    @info "regrid parity beam $b: n=$(length(x)), r_all=$(round(cor(x, y), digits=4)), " *
                          "r_lowroll=$(round(rl, digits=4)), med|Δ| = $(round(median(d), sigdigits=3)) m/s"
                end
            end
        end
    else
        @info "gliderad2cp reference outputs not found — skipping parity tests"
    end

end
