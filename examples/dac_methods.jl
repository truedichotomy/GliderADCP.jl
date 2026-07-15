# The DAC ladder made visible: U/V sections from the SAME pings solved against
# each of the three DAC methods — ADCP water-track (production reference),
# flight-model, onboard-DR — plus difference sections relative to the ADCP DAC.
# Generated for the delayed-mode route and the realtime-telemetered route.
#
# Only the DAC differs between panels, so the differences are pure
# DAC-referencing effects: depth-uniform per yo (a DAC error shifts the whole
# profile barotropically), pale for the flight model (~1–2 cm/s, unbiased),
# saturated track-correlated stripes for the onboard DR (2–4 cm/s anti-track —
# e.g. M59's two-week +4–5 cm/s block while the glider pointed upstream into
# the Gulf Stream). Evidence: validation doc 2026-07-15 entries; QA/QC §3b.
#
# Run (CairoMakie is found in the user's @ocean environment via the stacked load
# path; the package itself stays plot-free):
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/dac_methods.jl
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/dac_methods.jl m38

using GliderADCP
using CairoMakie
using DataFrames, Dates, Statistics, NaNStatistics
using Printf
include("missions.jl")

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

# symmetric color range: 99th percentile of |values|, rounded to a tidy step
sym99(As...) = quantile(abs.(reduce(vcat, [filter(isfinite, vec(A)) for A in As])), 0.99)

# Solve the same pings against the three DACs, restricted to yos whose ADCP DAC
# is genuinely ADCP-referenced — so the reference panels are what they claim to
# be and the differences compare methods, not coverage.
function solve3(nav, pings, fl; max_gap=30.0)
    dacs = Dict(
        "adcp"    => compute_dac(nav, pings; max_gap),
        "flight"  => compute_dac(nav, fl),
        "onboard" => compute_dac(nav),
    )
    keep = Set(dacs["adcp"].yo[dacs["adcp"].method .=== :adcp])
    invs = Dict(k => solve_inverse(pings, d) for (k, d) in dacs)
    common = intersect(Set(invs["adcp"].yo), keep)
    for k in keys(invs)
        invs[k] = invs[k][[y in common for y in invs[k].yo], :]
    end
    nfb = count(dacs["flight"].method[[y in common for y in dacs["flight"].yo]] .=== :onboard)
    nfb > 0 && @info "    note: flight-model DAC fell back to onboard on $nfb of $(length(common)) shown yos"
    return invs
end

function dac_methods_figure(invs, label, route, fname)
    sad = grid_profiles(invs["adcp"])
    sfm = grid_profiles(invs["flight"])
    sob = grid_profiles(invs["onboard"])
    @assert sad.t == sfm.t == sob.t && sad.z == sfm.z == sob.z
    dUf = sfm.U .- sad.U; dVf = sfm.V .- sad.V
    dUo = sob.U .- sad.U; dVo = sob.V .- sad.V

    crv = ceil(sym99(sad.U, sad.V, sfm.U, sfm.V, sob.U, sob.V) * 20) / 20
    crd = max(ceil(sym99(dUf, dVf, dUo, dVo) * 100) / 100, 0.03)

    nyo = length(sad.t)
    xt = round.(Int, range(1, nyo; length=min(8, nyo)))
    xtl = Dates.format.(sad.t[xt], "dd u")

    rows = [
        (sad.U, "U — ADCP water-track DAC (reference)", crv),
        (sad.V, "V — ADCP water-track DAC (reference)", crv),
        (sfm.U, "U — flight-model DAC", crv),
        (sfm.V, "V — flight-model DAC", crv),
        (sob.U, "U — onboard-DR DAC", crv),
        (sob.V, "V — onboard-DR DAC", crv),
        (dUf, "U: flight-model − ADCP", crd),
        (dVf, "V: flight-model − ADCP", crd),
        (dUo, "U: onboard-DR − ADCP", crd),
        (dVo, "V: onboard-DR − ADCP", crd),
    ]

    fig = Figure(size=(2200, 5 * 270 + 120))
    Label(fig[0, 1:3], "$label — inverse solution under the three DAC methods ($route)";
          fontsize=22, font=:bold, padding=(0, 0, 8, 0))
    local hmv, hmd
    for (i, (A, title, cr)) in enumerate(rows)
        r = cld(i, 2); c = isodd(i) ? 1 : 2
        ax = Axis(fig[r, c]; ylabel=c == 1 ? "depth (m)" : "", yreversed=true,
                  title=title, xticks=(xt, xtl),
                  xlabel=r == 5 ? "segment midpoints" : "")
        hm = heatmap!(ax, 1:nyo, sad.z, permutedims(A);
                      colormap=:balance, colorrange=(-cr, cr))
        r <= 3 ? (hmv = hm) : (hmd = hm)
    end
    Colorbar(fig[1:3, 3], hmv; label="velocity (m/s)")
    Colorbar(fig[4:5, 3], hmd; label="difference (m/s)")
    save(joinpath(OUT, fname), fig)

    p95(A) = 100 * quantile(abs.(filter(isfinite, vec(A))), 0.95)
    for (nm, dU, dV) in (("flight − ADCP ", dUf, dVf), ("onboard − ADCP", dUo, dVo))
        @printf("    Δ(%s): median |Δu| %.2f cm/s, |Δv| %.2f cm/s (p95 %.2f / %.2f)\n",
                nm, 100nanmedian(abs.(dU)), 100nanmedian(abs.(dV)), p95(dU), p95(dV))
    end
    @info "    wrote $fname ($(nyo) yos shown)"
end

for key in selected_missions()
    m = MISSIONS[key]
    @info "════════ $(m.label) ════════"

    nav = load_seaexplorer_nav([joinpath(m.dir, "delayed/nav/logs"), joinpath(m.dir, "glimpse")];
                               stream="$(m.prefix).gli.sub")
    lat = round(nanmedian(nav.lat), digits=1)
    fl = flight_model(nav)

    pld = load_seaexplorer_pld(joinpath(m.dir, "delayed/pld1/logs"); stream="$(m.prefix).pld1.sub")
    ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                      !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
    ctd_t = datetime2unix.(pld.time[ok])
    ord = sortperm(ctd_t); ctd_t = ctd_t[ord]
    Tl = Float64.(pld.LEGATO_TEMPERATURE[ok])[ord]
    Sl = Float64.(pld.LEGATO_SALINITY[ok])[ord]
    Pl = Float64.(pld.LEGATO_PRESSURE[ok])[ord]
    c_true = soundspeed_from_ctd.(Sl, Tl, Pl, 5.0, lat)

    # ── delayed mode ──
    adcp = read_ad2cp(joinpath(m.dir, m.binary))
    bincfg = adcp.config
    apply_soundspeed!(adcp, soundspeed_correction(adcp, ctd_t, c_true))
    qc!(adcp)
    p_d = process_pings(adcp; lat, declination=magnetic_declination(nav, adcp.t))
    calibrate_shear_bias!(p_d)
    @info "  delayed mode:"
    dac_methods_figure(solve3(nav, p_d, fl), m.label, "delayed mode",
                       "$(m.label)_dac_methods_delayed.png")

    # ── realtime-telemetered ──
    tele = load_pld_adcp([joinpath(m.dir, "delayed/pld1/logs"), joinpath(m.dir, "glimpse")];
                         stream="$(m.prefix).pld1.sub",
                         cellsize=bincfg.cellsize, blanking=bincfg.blanking,
                         serial=bincfg.serial)
    onboard_soundspeed!(tele, ctd_t, Tl; salinity=bincfg.salinity_setting, lat=lat)
    apply_soundspeed!(tele, soundspeed_correction(tele, ctd_t, c_true))
    qc!(tele)
    p_t = process_pings(tele; lat, look=:down, declination=magnetic_declination(nav, tele.t))
    calibrate_shear_bias!(p_t)
    @info "  realtime-telemetered:"
    dac_methods_figure(solve3(nav, p_t, fl; max_gap=90.0), m.label, "realtime-telemetered",
                       "$(m.label)_dac_methods_telemetered.png")
end
