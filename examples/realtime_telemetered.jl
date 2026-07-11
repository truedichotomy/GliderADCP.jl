# Shore-side real-time product from the telemetered AD2CP pings (pld1.sub), and a
# three-way comparison: our telemetered inverse vs the delayed-mode product vs
# ALSEAMAR's own GLIMPSE real-time current estimate (AD2CP_*_c columns).
#
# Two distinct "real-time" data routes exist and must not be confused:
#   1. the $PNOR ASCII stream (ad2cp.raw) — logged on the payload computer, every
#      ensemble, all cells + amp/corr, but NOT transmitted; recovered with the
#      glider (see examples/realtime_vs_delayed.jl);
#   2. the AD2CP subset inside pld1.sub — one subsampled ensemble every ~30 s,
#      cells 1–6, beam velocities at 0.01 m/s, attitude/pressure, no amp/corr/BT.
#      THIS is what Iridium actually delivers mid-mission, and what this script
#      turns into a current product with the standard pipeline.
#
# Configuration (cell size, blanking, serial, configured salinity) is taken from
# the mission's binary config here for exactness; in true real-time use it comes
# from the deployment plan. The onboard sound speed is reconstructed from the
# configured salinity + payload CTD temperature, so the standard correction runs.
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/realtime_telemetered.jl        # all missions
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/realtime_telemetered.jl m38     # one mission

using GliderADCP, SeaExplorerIO
using DataFrames, Dates, Statistics, NaNStatistics
using Printf
include("missions.jl")

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

interp1(xq, x, y) = begin
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
    x1, x2 = x[i], x[i+1]
    x2 == x1 ? y[i] : y[i] + (y[i+1] - y[i]) * (xq - x1) / (x2 - x1)
end

function agreement(a, b, col; nmin_a=10, nmin_b=3)
    j = innerjoin(a, b; on=[:yo, :z], makeunique=true)
    c1, c2 = j[!, col], j[!, Symbol(col, :_1)]
    m = (j.nobs .> nmin_a) .&& (j.nobs_1 .> nmin_b) .&& isfinite.(c1) .&& isfinite.(c2)
    d = c1[m] .- c2[m]
    (j=j, m=m, col=col, n=count(m), r=cor(c1[m], c2[m]), rms=sqrt(mean(d .^ 2)), bias=mean(d))
end

for key in selected_missions()
    m = MISSIONS[key]
    @info "════════ $(m.label): shore-side real-time (telemetered pld1.sub) ════════"
    srcs = [joinpath(m.dir, "delayed/pld1/logs"), joinpath(m.dir, "glimpse")]

    # ── configuration + shared references ──
    bincfg = read_ad2cp(joinpath(m.dir, m.binary)).config
    nav = load_seaexplorer_nav(joinpath(m.dir, "delayed/nav/logs"); stream="$(m.prefix).gli.sub")
    lat = round(nanmedian(nav.lat), digits=1)
    dac = compute_dac(nav)
    pld = load_seaexplorer_pld(joinpath(m.dir, "delayed/pld1/logs"); stream="$(m.prefix).pld1.sub")
    ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                      !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
    ctd_t = datetime2unix.(pld.time[ok])
    ord = sortperm(ctd_t); ctd_t = ctd_t[ord]
    Tl = Float64.(pld.LEGATO_TEMPERATURE[ok])[ord]
    Sl = Float64.(pld.LEGATO_SALINITY[ok])[ord]
    Pl = Float64.(pld.LEGATO_PRESSURE[ok])[ord]
    c_true = soundspeed_from_ctd.(Sl, Tl, Pl, 5.0, lat)

    # ── telemetered route → product ──
    tele = load_pld_adcp(srcs; stream="$(m.prefix).pld1.sub",
        cellsize=bincfg.cellsize, blanking=bincfg.blanking, serial=bincfg.serial)
    Tping = [interp1(t, ctd_t, Tl) for t in tele.t]
    c_used = soundspeed_from_ctd.(bincfg.salinity_setting, Tping, tele.pressure, 5.0, lat)
    tele = load_pld_adcp(srcs; stream="$(m.prefix).pld1.sub",
        cellsize=bincfg.cellsize, blanking=bincfg.blanking, serial=bincfg.serial,
        soundspeed=c_used)
    apply_soundspeed!(tele, soundspeed_correction(tele, ctd_t, c_true))
    qc!(tele)
    p_t = process_pings(tele; lat=lat, look=:down, declination=magnetic_declination(nav, tele.t))
    calibrate_shear_bias!(p_t)
    inv_t = solve_inverse(p_t, dac)
    w_t = solve_w(p_t, dac)
    @info "  telemetered: $(length(tele)) pings ($(ncells(tele)) cells) → " *
          "$(length(unique(inv_t.yo))) yos, $(nrow(inv_t)) bins"

    # ── delayed-mode reference ──
    adcp = read_ad2cp(joinpath(m.dir, m.binary))
    apply_soundspeed!(adcp, soundspeed_correction(adcp, ctd_t, c_true))
    qc!(adcp)
    p_d = process_pings(adcp; lat=lat, declination=magnetic_declination(nav, adcp.t))
    calibrate_shear_bias!(p_d)
    inv_d = solve_inverse(p_d, dac)
    w_d = solve_w(p_d, dac)

    @info "  agreement with the delayed inverse on common (yo, z) bins:"
    stats = Dict("u" => agreement(inv_d, inv_t, :u), "v" => agreement(inv_d, inv_t, :v),
                 "w" => agreement(w_d, w_t, :w))
    for k in ("u", "v", "w")
        s = stats[k]
        @printf "    %s: n=%5d  r=%.4f  rms=%.4f m/s  bias=%+.4f\n" k s.n s.r s.rms s.bias
    end

    # ── ALSEAMAR GLIMPSE product (AD2CP_*_c), binned to the same (yo, z) grid ──
    acols = ["AD2CP_Unorth_c", "AD2CP_Ueast_c", "AD2CP_QF_c"]
    gt = SeaExplorerIO.read_pld(joinpath(m.dir, "glimpse"), acols;
        stream="$(m.prefix).pld1.sub", skip_empty=false)
    fp = findall(i -> isfinite(gt["AD2CP_Unorth_c"][i]) && isfinite(gt["AD2CP_Ueast_c"][i]),
        1:length(gt))
    inv_a = DataFrame()
    if !isempty(fp)
        tg = datetime2unix.(gt.time[fp])
        zg = [interp1(t, ctd_t, Pl) for t in tg]         # glider depth at the estimate
        zlev = sort(unique(inv_d.z))
        yowin = [(r.yo, datetime2unix(r.t_start), datetime2unix(r.t_end)) for r in eachrow(dac)]
        rows = Dict{Tuple{Int,Float64},Vector{Tuple{Float64,Float64}}}()
        for (i, t) in enumerate(tg)
            k = findfirst(w -> w[2] <= t <= w[3], yowin)
            k === nothing && continue
            z = zlev[argmin(abs.(zlev .- zg[i]))]
            push!(get!(rows, (yowin[k][1], z), Tuple{Float64,Float64}[]),
                (gt["AD2CP_Ueast_c"][fp[i]], gt["AD2CP_Unorth_c"][fp[i]]))
        end
        tmid = Dict(r.yo => r.t_mid for r in eachrow(dac))
        inv_a = DataFrame(yo=Int[], t_mid=DateTime[], z=Float64[], u=Float64[], v=Float64[], nobs=Int[])
        for ((y, z), vals) in rows
            push!(inv_a, (y, tmid[y], z, mean(first.(vals)), mean(last.(vals)), length(vals)))
        end
        # keep the section time axis identical to the delayed/telemetered panels:
        # same yo set (ALSEAMAR also has values on short segments our solvers skip —
        # drop those columns), and pad any yo it lacks so the grids align 1:1
        yoset = Set(unique(inv_d.yo))
        inv_a = inv_a[[y in yoset for y in inv_a.yo], :]
        present = Set(unique(inv_a.yo))
        for y in setdiff(yoset, present)
            push!(inv_a, (y, tmid[y], zlev[1], NaN, NaN, 0))
        end
        @info "  ALSEAMAR product, binned to (yo, z) and compared to the delayed inverse:"
        for col in (:u, :v)
            s = agreement(inv_d, inv_a, col)
            @printf "    %s: n=%5d  r=%.4f  rms=%.4f m/s  bias=%+.4f\n" col s.n s.r s.rms s.bias
        end
    else
        @info "  (no ALSEAMAR AD2CP_*_c product in this mission's GLIMPSE export)"
    end

    # ── figures ──
    @eval using CairoMakie
    sym99(As...) = quantile(abs.(reduce(vcat, [filter(isfinite, vec(A)) for A in As])), 0.99)
    sec_d = grid_profiles(inv_d); sec_t = grid_profiles(inv_t)
    panels = [(sec_d, :U, "U — delayed (full-resolution binary)"),
              (sec_d, :V, "V — delayed"),
              (sec_t, :U, "U — shore-side real-time (telemetered pld1.sub)"),
              (sec_t, :V, "V — real-time")]
    crUV = ceil(sym99(sec_d.U, sec_d.V, sec_t.U, sec_t.V) * 20) / 20
    if !isempty(inv_a)
        sec_a = grid_profiles(inv_a)
        append!(panels, [(sec_a, :U, "U — ALSEAMAR GLIMPSE product (binned)"),
                         (sec_a, :V, "V — ALSEAMAR")])
    end
    fig = plot_sections(panels; colorrange=(-crUV, crUV))
    save(joinpath(OUT, "$(m.label)_telemetered_sections.png"), fig)

    fig2 = Figure(size=(1500, 460))
    su, sv = stats["u"], stats["v"]
    for (i, s, lab) in ((1, su, "u"), (2, sv, "v"))
        ax = Axis(fig2[1, i]; xlabel="delayed $(lab) (m/s)", ylabel="real-time $(lab) (m/s)",
            title=@sprintf("telemetered inverse %s: r=%.3f, rms=%.1f mm/s", lab, s.r, 1000s.rms),
            aspect=1)
        scatter!(ax, s.j[s.m, s.col], s.j[s.m, Symbol(s.col, :_1)];
            markersize=2, color=(:steelblue, 0.25))
        ablines!(ax, 0, 1; color=:black, linestyle=:dash)
    end
    ax3 = Axis(fig2[1, 3]; xlabel="rms difference vs delayed (mm/s)", ylabel="depth (m)",
        yreversed=true, title="by depth")
    for (s, color, lab) in ((stats["u"], :dodgerblue, "tele u"), (stats["v"], :navy, "tele v"))
        zc, rmsz = Float64[], Float64[]
        for z1 in 0:50:950
            mz = s.m .&& (z1 .<= s.j.z .< z1 + 50)
            count(mz) < 30 && continue
            push!(zc, z1 + 25)
            push!(rmsz, 1000 * sqrt(mean((s.j[mz, s.col] .- s.j[mz, Symbol(s.col, :_1)]) .^ 2)))
        end
        lines!(ax3, rmsz, zc; color, label=lab)
    end
    if !isempty(inv_a)
        for (col, color, lab) in ((:u, :darkorange, "ALSEAMAR u"), (:v, :firebrick, "ALSEAMAR v"))
            s = agreement(inv_d, inv_a, col)
            zc, rmsz = Float64[], Float64[]
            for z1 in 0:50:950
                mz = s.m .&& (z1 .<= s.j.z .< z1 + 50)
                count(mz) < 30 && continue
                push!(zc, z1 + 25)
                push!(rmsz, 1000 * sqrt(mean((s.j[mz, col] .- s.j[mz, Symbol(col, :_1)]) .^ 2)))
            end
            lines!(ax3, rmsz, zc; color, label=lab)
        end
    end
    axislegend(ax3; position=:rb)
    save(joinpath(OUT, "$(m.label)_telemetered_vs_delayed.png"), fig2; px_per_unit=2)
    @info "  wrote $(m.label)_telemetered_sections.png, $(m.label)_telemetered_vs_delayed.png"
end
