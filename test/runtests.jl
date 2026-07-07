using GliderADCP
using Test
using Dates
using DataFrames
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

end
