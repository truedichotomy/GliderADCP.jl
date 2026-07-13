# GliderADCP.jl

Pure-Julia processing of glider-mounted ADCP data into absolute ocean velocity
profiles — from the raw instrument binary to referenced, quality-controlled velocity
sections. Currently supports the Nortek AD2CP, validated end-to-end on four Alseamar
SeaExplorer missions; Slocum ingestion is implemented but not yet exercised on a real
Slocum dataset.

![GliderADCP.jl processing pipeline](assets/pipeline.svg)

## Features

- **Four input routes, one structure** — the full data-tier taxonomy: native `.ad2cp`
  binary (**delayed-mode**; bit-identical to the Nortek MIDAS export — no Windows/MIDAS
  step needed), MIDAS netCDF, the `$PNOR` stream (**realtime-onboard**, payload-logged),
  and the telemetered `pld1.sub` AD2CP subset (**realtime-telemetered** — what shore has
  mid-mission; shore-side realtime products build on it). SeaExplorer nav/payload
  parsers (multi-route glider-computer + GLIMPSE merge) and Slocum table ingestion
  (implemented, awaiting real-data validation).
- **A validated common trunk**: sound-speed correction (TEOS-10), composable QC (the
  first cell is kept by default — validated clean for ≥ 0.5 m blanking configurations),
  IGRF declination, exact 3-beam beam→XYZ→ENU transform, isobaric regridding, and a
  per-mission range-dependent ("shear") bias calibration.
- **Both published velocity solutions** over identical inputs: the Visbeck-style
  least-squares **inverse** with DAC / bottom-track / smoothness constraints (the
  recommended product) and the lADCP-tradition **shear method**; plus flight-model-free
  **vertical water velocity**.
- **References from navigation**: fix-to-fix depth-averaged currents, surface GPS
  drift, and bottom-track over-ground velocities.
- **Products**: depth-matched sections, provenance-rich netCDF export, Makie plotting
  extension, DIVAnd mapping example.

## Quick start

```julia
using GliderADCP

adcp  = load_ad2cp("sea064_M38.ad2cp")            # delayed-mode binary (or MIDAS .nc)
# realtime tiers: load_pnor (onboard $PNOR stream) · load_pld_adcp (telemetered pld1.sub)
nav   = load_seaexplorer_nav("delayed/nav/logs")
qc!(adcp)
pings = process_pings(adcp; lat=69.5,
                      declination=magnetic_declination(nav, adcp.t))
calibrate_shear_bias!(pings)
dac   = compute_dac(nav)
btv   = bt_velocity(adcp)
prof  = solve_inverse(pings, dac; bt=btv)         # yo × depth-bin velocity table
sec   = grid_profiles(prof)
export_sections("sections.nc", sec)
```

Start with the [Tutorial](tutorial.md) — it explains the science behind every step and
how to judge the results — and read the [QA/QC guide](qaqc.md) before trusting a
mission's numbers: most glider-ADCP data problems are silent, and it catalogues every
one this project has hit. The design record (plan, verified literature, analyses of the
Python reference implementations, format specifications, and the validation report)
lives in the repository under `PLAN.md` and `docs/research/`.

## Validation summary

| check | result (four validated missions: M37, M38, M48, M59) |
|---|---|
| native binary reader vs MIDAS netCDF | bit-identical on three missions (M38/M48/M59; max \|Δvel\| = 0) |
| transform parity vs Python `gliderad2cp` | machine-exact (max Δ 2×10⁻¹⁶ m/s) |
| synthetic truth (both solvers, incl. bottom-track-only referencing) | recovered within bin discretization |
| DAC closure | median 1–2 mm/s on every mission |
| dive vs climb consistency (inverse) | median \|Δ\| ≈ 2 cm/s |
| shear vs inverse agreement (the health metric) | r = 0.90–0.98 at 3–6 cm/s rms |
| shallow bins vs surface GPS drift (M38) | median \|Δ\| = 4 cm/s |
| realtime-onboard (`$PNOR`) inverse vs delayed | 3.2–5.1 mm/s rms, zero bias |
| realtime-telemetered (`pld1.sub`) inverse vs delayed | 28–45 mm/s rms, \|bias\| ≤ 0.8 mm/s |
| bottom track | screened by default: all of M38/M59's false locks rejected; M37/M48's genuine locks pass |

The bottom-track screening exists because 99.7 % of the reference mission's BT locks
proved to be **false** (near-field water-borne targets) — see the
[QA/QC guide](qaqc.md) for that finding and every other data pitfall, including why
the first cell is kept (the fleet's 0.7 m blanking makes it clean) and why shore-side
realtime products come from the telemetered route.
