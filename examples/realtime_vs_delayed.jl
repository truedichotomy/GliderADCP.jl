# Real-time ($PNOR telemetry stream) vs delayed-mode (.ad2cp binary) products,
# parameterized over missions — the Task-5 methodology (see
# m38_realtime_vs_delayed.jl and docs/research/m38_validation.md §Task 5)
# applied to any mission: both AD2CP routes run the IDENTICAL pipeline (same
# nav, CTD sound speed, QC, declination, shear-bias calibration, DAC), so the
# data source is the only difference.
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/realtime_vs_delayed.jl m37 m59

using GliderADCP
using DataFrames, Dates, Statistics, NaNStatistics
using Printf

const GDATA = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData"
const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

const MISSIONS = Dict(
    "m37" => (label="M37", dir=joinpath(GDATA, "sea064-20221021-norse-janmayen-complete"),
              bin="ad2cp/sea064_M37.ad2cp", prefix="37", lat=70.9),
    "m38" => (label="M38", dir=joinpath(GDATA, "sea064-20221102-norse-lofoten-complete"),
              bin="ad2cp/102381_sea064_M38/sea064_M38.ad2cp", prefix="38", lat=69.5),
    "m59" => (label="M59", dir=joinpath(GDATA, "sea064-20240720-nesma-passengers-complete"),
              bin="ad2cp/sea064_M59.ad2cp", prefix="59", lat=39.2),
)

function agreement(a, b, col; nmin=10)
    j = innerjoin(a, b; on=[:yo, :z], makeunique=true)
    c1, c2 = j[!, col], j[!, Symbol(col, :_1)]
    m = (j.nobs .> nmin) .&& (j.nobs_1 .> nmin) .&& isfinite.(c1) .&& isfinite.(c2)
    d = c1[m] .- c2[m]
    (j=j, m=m, col=col, n=count(m), r=cor(c1[m], c2[m]), rms=sqrt(mean(d .^ 2)), bias=mean(d))
end

function compare_mission(cfg)
    @info "== $(cfg.label): real-time (\$PNOR stream) vs delayed (.ad2cp binary) =="
    adcp_d = read_ad2cp(joinpath(cfg.dir, cfg.bin))
    adcp_r = load_pnor(joinpath(cfg.dir, "delayed/pld1/logs"); stream="$(cfg.prefix).ad2cp.raw")
    cov_d, cov_r = coverage(adcp_d), coverage(adcp_r)
    @info "  delayed: $(cov_d.n) ens $(cov_d.t_start)→$(cov_d.t_end); " *
          "stream: $(cov_r.n) ens $(cov_r.t_start)→$(cov_r.t_end) " *
          "($(round(100cov_r.n / cov_d.n, digits=1))%)"
    nav = load_seaexplorer_nav(joinpath(cfg.dir, "delayed/nav/logs"); stream="$(cfg.prefix).gli.sub")
    dac = compute_dac(nav)
    pld = load_seaexplorer_pld(joinpath(cfg.dir, "delayed/pld1/logs"); stream="$(cfg.prefix).pld1.sub")
    ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                      !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
    ctd_t = datetime2unix.(pld.time[ok])
    c_ctd = soundspeed_from_ctd.(Float64.(pld.LEGATO_SALINITY[ok]),
        Float64.(pld.LEGATO_TEMPERATURE[ok]), Float64.(pld.LEGATO_PRESSURE[ok]), 5.0, cfg.lat)

    prods = Dict{String,DataFrame}()
    for (lab, a, look) in (("delayed", adcp_d, :auto), ("realtime", adcp_r, :down))
        apply_soundspeed!(a, soundspeed_correction(a, ctd_t, c_ctd))
        qc!(a)
        p = process_pings(a; lat=cfg.lat, look=look, declination=magnetic_declination(nav, a.t))
        calibrate_shear_bias!(p)
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
    return stats
end

using CairoMakie
for key in (isempty(ARGS) ? ["m37", "m59"] : lowercase.(ARGS))
    cfg = MISSIONS[key]
    stats = compare_mission(cfg)
    fig = Figure(size=(1000, 420))
    su = stats["inv u"]
    ax1 = Axis(fig[1, 1]; xlabel="delayed u (m/s)", ylabel="real-time u (m/s)",
        title=@sprintf("%s inverse u: r=%.4f, rms=%.1f mm/s", cfg.label, su.r, 1000su.rms), aspect=1)
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
    save(joinpath(OUT, "$(cfg.label)_realtime_vs_delayed.png"), fig; px_per_unit=2)
    @info "  wrote $(joinpath(OUT, "$(cfg.label)_realtime_vs_delayed.png"))"
end
