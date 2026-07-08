# End-to-end GliderADCP.jl example: sea064 M38 (NORSE Lofoten Basin, Nov 2022 – Mar 2023).
#
# Full mission: MIDAS netCDF + SeaExplorer nav/payload → sound-speed correction → QC →
# declination (IGRF) → ENU pings → DAC + bottom track → inverse & shear solutions →
# diagnostics → gridded sections (netCDF/CSV) → figures.
#
# Run (CairoMakie is found in the user's @ocean environment via the stacked load path;
# the package itself stays plot-free):
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/m38_currents.jl

using GliderADCP
using DataFrames, Dates, Statistics, NaNStatistics
using CSV, NCDatasets
using Printf

const MISSION = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete"
const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)
const LAT0 = 69.5   # mean mission latitude (for pressure→depth)

@info "1/7 load AD2CP + navigation"
adcp = load_ad2cp(joinpath(MISSION, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp.00000.nc"))
nav = load_seaexplorer_nav(joinpath(MISSION, "delayed/nav/logs"))
@info "    $adcp"

@info "2/7 sound-speed correction from payload CTD (pld1.sub, decimated)"
pld = load_seaexplorer_pld(joinpath(MISSION, "delayed/pld1/logs"); stream="pld1.sub")
ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                  !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
c_ctd = soundspeed_from_ctd.(Float64.(pld.LEGATO_SALINITY[ok]),
    Float64.(pld.LEGATO_TEMPERATURE[ok]), Float64.(pld.LEGATO_PRESSURE[ok]), 5.0, LAT0)
scale = soundspeed_correction(adcp, datetime2unix.(pld.time[ok]), c_ctd)
ncorr = apply_soundspeed!(adcp, scale)
@info "    corrected $ncorr/$(length(adcp)) pings, median scale $(round(nanmedian(collect(scale)), digits=5))"

@info "3/7 QC"
qstats = qc!(adcp)
@info "    rejected $(round(100qstats.total, digits=1))% of beam samples"

@info "4/7 declination (IGRF) + ENU pings + DAC + bottom track"
decl = magnetic_declination(nav, adcp.t)
@info "    declination $(round(nanminimum(decl), digits=2))..$(round(nanmaximum(decl), digits=2)) °E"
pings = process_pings(adcp; lat=LAT0, declination=decl)
bslopes = calibrate_shear_bias!(pings)               # range-dependent bias (Phase 7)
@info "    shear-bias slope $(round(bslopes[1], sigdigits=3)) s⁻¹ → residual $(round(bslopes[end], sigdigits=2))"
dac = compute_dac(nav)
btv = bt_velocity(adcp; max_range=28.0, declination=magnetic_declination(nav, adcp.bt.t))
@info "    $(nrow(dac)) DAC segments, $(nrow(btv)) bottom-track fixes"

@info "5/7 velocity solutions (inverse, shear, vertical) — all yos"
inv = solve_inverse(pings, dac; bt=(nrow(btv) > 0 ? btv : nothing))
shr = solve_shear(pings, dac)
wdir = solve_w(pings, dac)                       # w, direct: binned U_rel + dP/dt
winv = solve_w(pings, dac; method=:inverse)      # w, inverse machinery, pressure-anchored
@info "    inverse: $(length(unique(inv.yo))) yos, $(nrow(inv)) bins; " *
      "shear: $(length(unique(shr.yo))) yos; w: $(length(unique(wdir.yo))) yos"

@info "6/7 diagnostics"
# (a) method intercomparison on common (yo, z) bins
j = innerjoin(inv, shr; on=[:yo, :z], makeunique=true)
gm = (j.nobs .> 10) .&& (j.nobs_1 .>= 4)
r_u = cor(j.u[gm], j.u_1[gm]);  rms_u = sqrt(mean((j.u[gm] .- j.u_1[gm]) .^ 2))
r_v = cor(j.v[gm], j.v_1[gm]);  rms_v = sqrt(mean((j.v[gm] .- j.v_1[gm]) .^ 2))
@printf "    shear vs inverse:  r_u=%.3f r_v=%.3f  rms_u=%.3f rms_v=%.3f m/s  (n=%d)\n" r_u r_v rms_u rms_v count(gm)

# (b) dive vs climb consistency (same yo DAC, independent half-segments)
function invert_window(p, t1, t2, dacu, dacv)
    idx = findall(t -> t1 <= t <= t2, p.t)
    length(idx) < 60 && return nothing
    gd = filter(isfinite, p.depth[idx])
    isempty(gd) && return nothing
    invert_segment(view(p.E, :, idx), view(p.N, :, idx), view(p.celldepth, :, idx),
        p.t[idx], maximum(gd); dacu, dacv)
end
dc_u = Tuple{Float64,Float64}[]; dc_v = Tuple{Float64,Float64}[]
for row in eachrow(dac[1:2:end, :])          # every other yo
    idx = segment_indices(pings, row.t_start, row.t_end)
    length(idx) < 200 && continue
    dsub = [isfinite(d) ? d : -Inf for d in pings.depth[idx]]
    tsplit = pings.t[idx[argmax(dsub)]]
    sd = invert_window(pings, datetime2unix(row.t_start), tsplit, row.u, row.v)
    sc = invert_window(pings, tsplit, datetime2unix(row.t_end), row.u, row.v)
    (sd === nothing || sc === nothing) && continue
    zc = Dict(zip(sc.z, eachindex(sc.z)))
    for (k, z) in enumerate(sd.z)
        haskey(zc, z) || continue
        kk = zc[z]
        (sd.nobs[k] > 5 && sc.nobs[kk] > 5) || continue
        push!(dc_u, (sd.u[k], sc.u[kk]))
        push!(dc_v, (sd.v[k], sc.v[kk]))
    end
end
rdc_u = cor(first.(dc_u), last.(dc_u)); mad_u = median(abs.(first.(dc_u) .- last.(dc_u)))
rdc_v = cor(first.(dc_v), last.(dc_v)); mad_v = median(abs.(first.(dc_v) .- last.(dc_v)))
@printf "    dive vs climb:     r_u=%.3f r_v=%.3f  med|Δu|=%.3f med|Δv|=%.3f m/s  (n=%d)\n" rdc_u rdc_v mad_u mad_v length(dc_u)

# (c) DAC closure of the inverse solution
clos = Float64[]
for row in eachrow(dac)
    sub = inv[(inv.yo .== row.yo) .& (inv.nobs .> 0), :]
    isempty(sub) && continue
    push!(clos, hypot(mean(sub.u) - row.u, mean(sub.v) - row.v))
end
@printf "    DAC closure:       median |Δ| = %.3f m/s over %d yos\n" median(clos) length(clos)

@info "7/7 gridding, export, figures"
sec = grid_profiles(inv)
sec_s = grid_profiles(shr)
attrs = Dict{String,Any}(
    "mission" => "sea064 M38, NORSE Lofoten Basin, 2022-11..2023-03",
    "instrument" => "Nortek Glider AD2CP 1 MHz SN102381",
    "declination" => "IGRF per ping",
    "qc_rejected_fraction" => round(qstats.total, digits=4))
export_sections(joinpath(OUT, "M38_sections_inverse.nc"), sec;
    attrs=merge(attrs, Dict{String,Any}("method" => "Visbeck inverse, DAC + bottom track")))
export_sections(joinpath(OUT, "M38_sections_shear.nc"), sec_s;
    attrs=merge(attrs, Dict{String,Any}("method" => "shear + DAC referencing")))
CSV.write(joinpath(OUT, "M38_profiles_inverse.csv"), inv)
CSV.write(joinpath(OUT, "M38_profiles_shear.csv"), shr)
CSV.write(joinpath(OUT, "M38_profiles_w_direct.csv"), wdir)
CSV.write(joinpath(OUT, "M38_profiles_w_inverse.csv"), winv)

using CairoMakie
# data-driven symmetric color ranges (99th percentile of |values|, shared per figure)
sym99(As...) = quantile(abs.(reduce(vcat, [filter(isfinite, vec(A)) for A in As])), 0.99)
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
save(joinpath(OUT, "M38_UV_sections.png"), fig)
figw = plot_sections([(sec_wd, :U, "w — direct (U_rel + dP/dt, binned)"),
                      (sec_wi, :U, "w — inverse (pressure-anchored)")];
    colorrange=(-crW, crW))
save(joinpath(OUT, "M38_w_sections.png"), figw)

fig2 = Figure(size=(1200, 480))
ax1 = Axis(fig2[1, 1]; title=@sprintf("shear vs inverse (u): r=%.2f, rms=%.3f m/s", r_u, rms_u),
    xlabel="inverse u (m/s)", ylabel="shear u (m/s)", aspect=1)
scatter!(ax1, j.u[gm], j.u_1[gm]; markersize=2, color=(:steelblue, 0.25))
ablines!(ax1, 0, 1; color=:black)
ax2 = Axis(fig2[1, 2]; title=@sprintf("dive vs climb (u): r=%.2f, med|Δ|=%.3f m/s", rdc_u, mad_u),
    xlabel="dive u (m/s)", ylabel="climb u (m/s)", aspect=1)
scatter!(ax2, first.(dc_u), last.(dc_u); markersize=3, color=(:darkorange, 0.35))
ablines!(ax2, 0, 1; color=:black)
save(joinpath(OUT, "M38_diagnostics.png"), fig2)

println("\nOutputs in $(OUT):")
foreach(f -> println("  ", f), readdir(OUT))
