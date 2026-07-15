# Realtime-onboard ($PNOR stream) vs delayed-mode (.ad2cp binary) products, run
# symmetrically over missions — the Task-5 methodology
# (docs/research/m38_validation.md §Task 5) applied to any mission: both AD2CP
# routes run the IDENTICAL pipeline (same nav, CTD sound speed, QC, declination,
# shear-bias calibration, DAC), so the data source is the only difference.
#
# Data-route taxonomy: delayed-mode (.ad2cp binary, post-recovery) ·
# realtime-onboard ($PNOR stream — payload-logged; useful in real time only to
# something like a backseat driver consuming it on the vehicle) ·
# realtime-telemetered (AD2CP subset in pld1.sub — what shore receives; see
# examples/realtime_telemetered.jl). This script bounds what an ONBOARD consumer
# of $PNOR could compute; shore-side realtime products should be built from the
# telemetered pld data instead.
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/realtime_onboard.jl
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/realtime_onboard.jl m37 m59

using GliderADCP
using DataFrames, Dates, Statistics, NaNStatistics
using Printf
include("missions.jl")

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

function agreement(a, b, col; nmin=10)
    j = innerjoin(a, b; on=[:yo, :z], makeunique=true)
    c1, c2 = j[!, col], j[!, Symbol(col, :_1)]
    m = (j.nobs .> nmin) .&& (j.nobs_1 .> nmin) .&& isfinite.(c1) .&& isfinite.(c2)
    d = c1[m] .- c2[m]
    (j=j, m=m, col=col, n=count(m), r=cor(c1[m], c2[m]), rms=sqrt(mean(d .^ 2)), bias=mean(d))
end

function compare_mission(m)
    @info "== $(m.label): real-time (\$PNOR stream) vs delayed (.ad2cp binary) =="
    adcp_d = read_ad2cp(joinpath(m.dir, m.binary))
    adcp_r = load_pnor(joinpath(m.dir, "delayed/pld1/logs"); stream="$(m.prefix).ad2cp.raw")
    cov_d, cov_r = coverage(adcp_d), coverage(adcp_r)
    @info "  delayed: $(cov_d.n) ens $(cov_d.t_start)→$(cov_d.t_end); " *
          "stream: $(cov_r.n) ens $(cov_r.t_start)→$(cov_r.t_end) " *
          "($(round(100cov_r.n / cov_d.n, digits=1))%)"
    nav = load_seaexplorer_nav(joinpath(m.dir, "delayed/nav/logs"); stream="$(m.prefix).gli.sub")
    lat = round(nanmedian(nav.lat), digits=1)
    pld = load_seaexplorer_pld(joinpath(m.dir, "delayed/pld1/logs"); stream="$(m.prefix).pld1.sub")
    ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                      !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
    ctd_t = datetime2unix.(pld.time[ok])
    c_ctd = soundspeed_from_ctd.(Float64.(pld.LEGATO_SALINITY[ok]),
        Float64.(pld.LEGATO_TEMPERATURE[ok]), Float64.(pld.LEGATO_PRESSURE[ok]), 5.0, lat)

    procs = Dict{String,ProcessedPings}()
    for (lab, a, look) in (("delayed", adcp_d, :auto), ("realtime", adcp_r, :down))
        apply_soundspeed!(a, soundspeed_correction(a, ctd_t, c_ctd))
        qc!(a)
        p = process_pings(a; lat=lat, look=look, declination=magnetic_declination(nav, a.t))
        calibrate_shear_bias!(p)
        procs[lab] = p
    end
    # one DAC for both routes (water-tracked from the delayed reference pings,
    # flight model filling any gaps), so the comparison isolates the data route
    dac = compute_dac(nav, procs["delayed"]; fallback=flight_model(nav))
    prods = Dict{String,DataFrame}()
    for (lab, p) in procs
        prods["$(lab)_inv"] = solve_inverse(p, dac)
        prods["$(lab)_shr"] = solve_shear(p, dac)
        prods["$(lab)_w"] = solve_w(p, dac)
    end

    stats = Dict(k => agreement(prods["delayed_$s"], prods["realtime_$s"], col; nmin)
                 for (k, s, col, nmin) in (("inv u", "inv", :u, 10), ("inv v", "inv", :v, 10),
                                           ("shr u", "shr", :u, 4), ("shr v", "shr", :v, 4),
                                           ("w", "w", :w, 10)))
    for k in ("inv u", "inv v", "shr u", "shr v", "w")
        s = stats[k]
        @printf "  %-6s n=%6d  r=%.4f  rms=%.4f m/s  bias=%+.4f m/s\n" k s.n s.r s.rms s.bias
    end
    nyo_d = length(unique(prods["delayed_inv"].yo))
    nyo_r = length(unique(prods["realtime_inv"].yo))
    @info "  yos solved — delayed: $nyo_d, real-time: $nyo_r"
    return stats, prods
end

using CairoMakie
sym99(As...) = quantile(abs.(reduce(vcat, [filter(isfinite, vec(A)) for A in As])), 0.99)
# gridded section of the real-time − delayed difference on common (yo, z) bins
function dgrid(s)
    d = DataFrame(yo=s.j.yo[s.m], t_mid=s.j.t_mid[s.m], z=s.j.z[s.m],
        u=s.j.u[s.m] .- s.j.u_1[s.m], v=s.j.v[s.m] .- s.j.v_1[s.m],
        nobs=min.(s.j.nobs[s.m], s.j.nobs_1[s.m]))
    grid_profiles(d)
end

for key in selected_missions()
    m = MISSIONS[key]
    stats, prods = compare_mission(m)

    # (a) side-by-side U/V sections, one shared color scale
    sec_d = grid_profiles(prods["delayed_inv"])
    sec_r = grid_profiles(prods["realtime_inv"])
    crUV = ceil(sym99(sec_d.U, sec_d.V, sec_r.U, sec_r.V) * 20) / 20
    figs = plot_sections([(sec_d, :U, "U (east) — delayed (.ad2cp binary)"),
                          (sec_d, :V, "V (north) — delayed"),
                          (sec_r, :U, "U (east) — real-time (\$PNOR stream)"),
                          (sec_r, :V, "V (north) — real-time")];
        colorrange=(-crUV, crUV))
    save(joinpath(OUT, "$(m.label)_realtime_onboard_sections.png"), figs)

    # (b) difference sections (real-time − delayed), inverse + shear, amplified scale
    dg_i, dg_s = dgrid(stats["inv u"]), dgrid(stats["shr u"])
    crD = ceil(sym99(dg_s.U, dg_s.V) * 200) / 200
    figd = plot_sections([(dg_i, :U, "ΔU — inverse (real-time − delayed)"),
                          (dg_i, :V, "ΔV — inverse"),
                          (dg_s, :U, "ΔU — shear method"),
                          (dg_s, :V, "ΔV — shear method")];
        colorrange=(-crD, crD))
    save(joinpath(OUT, "$(m.label)_realtime_onboard_diff_sections.png"), figd)

    # (c) scatter + rms-by-depth summary
    fig = Figure(size=(1000, 420))
    su = stats["inv u"]
    ax1 = Axis(fig[1, 1]; xlabel="delayed u (m/s)", ylabel="real-time u (m/s)",
        title=@sprintf("%s inverse u: r=%.4f, rms=%.1f mm/s", m.label, su.r, 1000su.rms), aspect=1)
    scatter!(ax1, su.j.u[su.m], su.j.u_1[su.m]; markersize=2, color=(:steelblue, 0.3))
    ablines!(ax1, 0, 1; color=:black, linestyle=:dash)
    ax2 = Axis(fig[1, 2]; xlabel="rms difference (mm/s)", ylabel="depth (m)",
        yreversed=true, title="real-time − delayed, by depth")
    for (key2, color) in (("inv u", :dodgerblue), ("inv v", :navy),
                          ("shr u", :darkorange), ("shr v", :firebrick), ("w", :seagreen))
        s = stats[key2]
        zc, rmsz = Float64[], Float64[]
        for z1 in 0:50:950
            mz = s.m .&& (z1 .<= s.j.z .< z1 + 50)
            count(mz) < 30 && continue
            push!(zc, z1 + 25)
            push!(rmsz, 1000 * sqrt(mean((s.j[mz, s.col] .- s.j[mz, Symbol(s.col, :_1)]) .^ 2)))
        end
        lines!(ax2, rmsz, zc; color, label=key2)
    end
    axislegend(ax2; position=:rb)
    save(joinpath(OUT, "$(m.label)_realtime_onboard.png"), fig; px_per_unit=2)
    @info "  wrote $(joinpath(OUT, "$(m.label)_realtime_onboard.png"))"
end
