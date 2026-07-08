# Task 5: real-time vs delayed-mode products, sea064 M38 (NORSE Lofoten Basin).
#
# Question: how much ocean-velocity accuracy is lost by processing the real-time
# $PNOR ASCII telemetry stream (0.01 m/s velocity quantization, 0.1° attitude, no
# accelerometer, no bottom-track records) instead of the full-resolution .ad2cp
# binary recovered after the mission?  Both sides run the IDENTICAL pipeline —
# same nav (gli.sub), CTD sound speed (pld1.sub), QC, declination, shear-bias
# calibration, DAC — so the AD2CP data source is the only difference.
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/m38_realtime_vs_delayed.jl

using GliderADCP
using DataFrames, Dates, Statistics, NaNStatistics
using Printf

const MISSION = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete"
const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)
const LAT0 = 69.5

@info "1/4 load both AD2CP sources + shared nav/CTD"
adcp_d = read_ad2cp(joinpath(MISSION, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
adcp_r = load_pnor(joinpath(MISSION, "delayed/pld1/logs"))
cov_d, cov_r = coverage(adcp_d), coverage(adcp_r)
@info "    delayed (binary):   $(cov_d.n) ens, $(cov_d.t_start) → $(cov_d.t_end)"
@info "    real-time (stream): $(cov_r.n) ens, $(cov_r.t_start) → $(cov_r.t_end)"
# M38: the payload stopped writing the stream on 2022-11-27; the binary adds only
# 750 sparse burst ensembles over the following three months.
nav = load_seaexplorer_nav(joinpath(MISSION, "delayed/nav/logs"))
dac = compute_dac(nav)
pld = load_seaexplorer_pld(joinpath(MISSION, "delayed/pld1/logs"); stream="pld1.sub")
ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                  !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
ctd_t = datetime2unix.(pld.time[ok])
c_ctd = soundspeed_from_ctd.(Float64.(pld.LEGATO_SALINITY[ok]),
    Float64.(pld.LEGATO_TEMPERATURE[ok]), Float64.(pld.LEGATO_PRESSURE[ok]), 5.0, LAT0)

@info "2/4 identical pipeline on both sources"
prods = Dict{String,DataFrame}()
for (lab, a, look) in (("delayed", adcp_d, :auto), ("realtime", adcp_r, :down))
    apply_soundspeed!(a, soundspeed_correction(a, ctd_t, c_ctd))
    qc!(a)
    # the stream carries no accelerometer, so look direction must be given explicitly
    p = process_pings(a; lat=LAT0, look=look, declination=magnetic_declination(nav, a.t))
    calibrate_shear_bias!(p)
    prods["$(lab)_inv"] = solve_inverse(p, dac)
    prods["$(lab)_shr"] = solve_shear(p, dac)
    prods["$(lab)_w"] = solve_w(p, dac)
end

@info "3/4 agreement on common (yo, z) bins"
function agreement(a, b, col; nmin=10)
    j = innerjoin(a, b; on=[:yo, :z], makeunique=true)
    c1, c2 = j[!, col], j[!, Symbol(col, :_1)]
    m = (j.nobs .> nmin) .&& (j.nobs_1 .> nmin) .&& isfinite.(c1) .&& isfinite.(c2)
    d = c1[m] .- c2[m]
    (j=j, m=m, n=count(m), r=cor(c1[m], c2[m]), rms=sqrt(mean(d .^ 2)), bias=mean(d))
end
stats = Dict(k => agreement(prods["delayed_$s"], prods["realtime_$s"], col; nmin)
             for (k, s, col, nmin) in (("inv u", "inv", :u, 10), ("inv v", "inv", :v, 10),
                                       ("shr u", "shr", :u, 4), ("shr v", "shr", :v, 4),
                                       ("w", "w", :w, 10)))
for k in ("inv u", "inv v", "shr u", "shr v", "w")
    s = stats[k]
    @printf "    %-6s n=%5d  r=%.4f  rms=%.4f m/s  bias=%+.4f m/s\n" k s.n s.r s.rms s.bias
end

@info "4/4 figures"
try
    @eval using CairoMakie

    # (a) sections: delayed vs real-time inverse U/V side by side
    sym99(As...) = quantile(abs.(reduce(vcat, [filter(isfinite, vec(A)) for A in As])), 0.99)
    sec_d = grid_profiles(prods["delayed_inv"])
    sec_r = grid_profiles(prods["realtime_inv"])
    crUV = ceil(sym99(sec_d.U, sec_d.V, sec_r.U, sec_r.V) * 20) / 20
    fig = plot_sections([(sec_d, :U, "U (east) — delayed (.ad2cp binary)"),
                         (sec_d, :V, "V (north) — delayed"),
                         (sec_r, :U, "U (east) — real-time (\$PNOR stream)"),
                         (sec_r, :V, "V (north) — real-time")];
        colorrange=(-crUV, crUV))
    save(joinpath(OUT, "M38_realtime_vs_delayed_sections.png"), fig)

    # (b) difference sections (real-time − delayed) on common (yo, z) bins,
    #     inverse and shear, one amplified color scale
    function dgrid(s)
        d = DataFrame(yo=s.j.yo[s.m], t_mid=s.j.t_mid[s.m], z=s.j.z[s.m],
            u=s.j.u[s.m] .- s.j.u_1[s.m], v=s.j.v[s.m] .- s.j.v_1[s.m],
            nobs=min.(s.j.nobs[s.m], s.j.nobs_1[s.m]))
        grid_profiles(d)
    end
    dg_i, dg_s = dgrid(stats["inv u"]), dgrid(stats["shr u"])
    crD = ceil(sym99(dg_s.U, dg_s.V) * 200) / 200
    figd = plot_sections([(dg_i, :U, "ΔU — inverse (real-time − delayed)"),
                          (dg_i, :V, "ΔV — inverse"),
                          (dg_s, :U, "ΔU — shear method"),
                          (dg_s, :V, "ΔV — shear method")];
        colorrange=(-crD, crD))
    save(joinpath(OUT, "M38_realtime_vs_delayed_diff_sections.png"), figd)

    # (c) scatter panels for every product + rms-by-depth summary
    fig2 = Figure(size=(1500, 900))
    panels = (("inv u", :u, 1, 1), ("inv v", :v, 1, 2), ("w", :w, 1, 3),
              ("shr u", :u, 2, 1), ("shr v", :v, 2, 2))
    for (key, col, ri, ci) in panels
        s = stats[key]
        ax = Axis(fig2[ri, ci]; xlabel="delayed (m/s)", ylabel="real-time (m/s)",
            title=@sprintf("%s: r=%.4f, rms=%.1f mm/s", key, s.r, 1000s.rms), aspect=1)
        scatter!(ax, s.j[s.m, col], s.j[s.m, Symbol(col, :_1)];
            markersize=2, color=(:steelblue, 0.25))
        ablines!(ax, 0, 1; color=:black, linestyle=:dash)
    end
    ax2 = Axis(fig2[2, 3]; xlabel="rms difference (mm/s)", ylabel="depth (m)",
        yreversed=true, title="real-time − delayed, by depth")
    for (key, col, color) in (("inv u", :u, :dodgerblue), ("inv v", :v, :navy),
                              ("shr u", :u, :darkorange), ("shr v", :v, :firebrick),
                              ("w", :w, :seagreen))
        s = stats[key]
        zc, rmsz = Float64[], Float64[]
        for z1 in 0:50:950
            mz = s.m .&& (z1 .<= s.j.z .< z1 + 50)
            count(mz) < 30 && continue
            push!(zc, z1 + 25)
            push!(rmsz, 1000 * sqrt(mean((s.j[mz, col] .- s.j[mz, Symbol(col, :_1)]) .^ 2)))
        end
        lines!(ax2, rmsz, zc; color, label=key)
    end
    axislegend(ax2; position=:rb)
    save(joinpath(OUT, "M38_realtime_vs_delayed.png"), fig2; px_per_unit=2)
    @info "    wrote 3 figures to $(OUT)"
catch err
    @warn "figures skipped (CairoMakie not on the load path?)" error = err
end
