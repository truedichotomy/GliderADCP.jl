# End-to-end GliderADCP.jl example: sea064 M37 (NORSE Jan Mayen, Oct 2022).
#
# Exercises the full stack on a second mission: native .ad2cp binary (no MIDAS
# export exists for this deployment) + multi-route SeaExplorer loading (glider-
# computer files and the GLIMPSE-server export merged and deduplicated) → QC →
# declination → ENU pings → shear-bias calibration → DAC (+ screened bottom
# track) → inverse & shear & vertical solutions → figures.
#
# The deployment folder holds several mission counters (M34/M35/M37) in the
# same log directories; the `stream` argument carries the mission number
# ("37.gli.sub") so only M37 segments are read.
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/m37_currents.jl

using GliderADCP
using DataFrames, Dates, Statistics, NaNStatistics
using Printf

const MISSION = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete"
const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

@info "1/6 load AD2CP (native binary) + navigation (delayed + GLIMPSE, deduplicated)"
adcp = read_ad2cp(joinpath(MISSION, "ad2cp/sea064_M37.ad2cp"))
@info "    $adcp"
cov = coverage(adcp)
@info "    coverage: $(nrow(cov.gaps)) gaps totalling $(round(cov.gap_total / 86400, digits=1)) days"
nav = load_seaexplorer_nav([joinpath(MISSION, "delayed/nav/logs"),
                            joinpath(MISSION, "glimpse")]; stream="37.gli.sub")
LAT0 = round(nanmedian(nav.lat), digits=1)
@info "    $(length(nav)) nav rows, mean latitude $(LAT0)°N"

@info "2/6 sound-speed correction from payload CTD (pld1.sub, delayed + GLIMPSE)"
pld = load_seaexplorer_pld([joinpath(MISSION, "delayed/pld1/logs"),
                            joinpath(MISSION, "glimpse")]; stream="37.pld1.sub")
ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                  !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
c_ctd = soundspeed_from_ctd.(Float64.(pld.LEGATO_SALINITY[ok]),
    Float64.(pld.LEGATO_TEMPERATURE[ok]), Float64.(pld.LEGATO_PRESSURE[ok]), 5.0, LAT0)
scale = soundspeed_correction(adcp, datetime2unix.(pld.time[ok]), c_ctd)
ncorr = apply_soundspeed!(adcp, scale)
@info "    corrected $ncorr/$(length(adcp)) pings, median scale $(round(nanmedian(collect(scale)), digits=5))"

@info "3/6 QC + declination + ENU pings + shear-bias calibration"
qstats = qc!(adcp)
@info "    rejected $(round(100qstats.total, digits=1))% of beam samples"
decl = magnetic_declination(nav, adcp.t)
pings = process_pings(adcp; lat=LAT0, declination=decl)
bslopes = calibrate_shear_bias!(pings)
@info "    shear-bias slope $(round(bslopes[1], sigdigits=3)) s⁻¹ → residual $(round(bslopes[end], sigdigits=2))"

@info "4/6 DAC + bottom track (screened)"
dac = compute_dac(nav)
btv = bt_velocity(adcp; max_range=28.0, declination=magnetic_declination(nav, adcp.bt.t))
@info "    $(nrow(dac)) DAC segments, $(nrow(btv)) bottom-track fixes survive screening"

@info "5/6 velocity solutions"
inv = solve_inverse(pings, dac; bt=(nrow(btv) > 0 ? btv : nothing))
shr = solve_shear(pings, dac)
wdir = solve_w(pings, dac)
winv = solve_w(pings, dac; method=:inverse)

# method intercomparison + DAC closure (the workflow health checks)
j = innerjoin(inv, shr; on=[:yo, :z], makeunique=true)
gm = (j.nobs .> 10) .&& (j.nobs_1 .>= 4)
@printf "    shear vs inverse:  r_u=%.3f r_v=%.3f  rms_u=%.3f rms_v=%.3f m/s  (n=%d)\n" cor(j.u[gm], j.u_1[gm]) cor(j.v[gm], j.v_1[gm]) sqrt(mean((j.u[gm] .- j.u_1[gm]) .^ 2)) sqrt(mean((j.v[gm] .- j.v_1[gm]) .^ 2)) count(gm)
clos = Float64[]
for row in eachrow(dac)
    sub = inv[(inv.yo .== row.yo) .& (inv.nobs .> 0), :]
    isempty(sub) || push!(clos, hypot(mean(sub.u) - row.u, mean(sub.v) - row.v))
end
@printf "    DAC closure:       median |Δ| = %.3f m/s over %d yos\n" median(clos) length(clos)

@info "6/6 figures"
using CairoMakie
sym99(As...) = quantile(abs.(reduce(vcat, [filter(isfinite, vec(A)) for A in As])), 0.99)
sec = grid_profiles(inv)
sec_s = grid_profiles(shr)
sec_wd = grid_profiles(wdir; fields=(:w, :w))
sec_wi = grid_profiles(winv; fields=(:w, :w))
crUV = ceil(sym99(sec.U, sec.V, sec_s.U, sec_s.V) * 20) / 20
crW = ceil(sym99(sec_wd.U, sec_wi.U) * 200) / 200
@info "    color ranges: horizontal ±$(crUV) m/s, vertical ±$(crW) m/s"
fig = plot_sections([(sec, :U, "U (east) — inverse"),
                     (sec, :V, "V (north) — inverse"),
                     (sec_s, :U, "U (east) — shear method"),
                     (sec_s, :V, "V (north) — shear method")];
    colorrange=(-crUV, crUV))
save(joinpath(OUT, "M37_UV_sections.png"), fig)
figw = plot_sections([(sec_wd, :U, "w — direct (U_rel + dP/dt, binned)"),
                      (sec_wi, :U, "w — inverse (pressure-anchored)")];
    colorrange=(-crW, crW))
save(joinpath(OUT, "M37_w_sections.png"), figw)

fig2 = Figure(size=(1200, 480))
ax1 = Axis(fig2[1, 1]; title="shear vs inverse (u)", xlabel="inverse u (m/s)",
    ylabel="shear u (m/s)", aspect=1)
scatter!(ax1, j.u[gm], j.u_1[gm]; markersize=2, color=(:steelblue, 0.25))
ablines!(ax1, 0, 1; color=:black)
ax2 = Axis(fig2[1, 2]; title="per-yo DAC closure", xlabel="yo", ylabel="|Δ| (m/s)")
scatter!(ax2, eachindex(clos), clos; markersize=4, color=:darkorange)
save(joinpath(OUT, "M37_diagnostics.png"), fig2)

println("\nM37 outputs in $(OUT):")
foreach(f -> startswith(f, "M37") && println("  ", f), readdir(OUT))
