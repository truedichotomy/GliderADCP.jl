# GliderADCP.jl

Pure-Julia processing of glider-mounted Nortek AD2CP data (SeaExplorer, Slocum) into
absolute ocean velocity profiles — from the raw instrument binary to referenced,
quality-controlled velocity sections.

![GliderADCP.jl processing pipeline](assets/pipeline.svg)

## Features

- **Three input routes, one structure**: native `.ad2cp` binary (bit-identical to the
  Nortek MIDAS export — no Windows/MIDAS step needed), MIDAS netCDF, and the real-time
  `$PNOR` telemetry stream; SeaExplorer nav/payload parsers and Slocum table ingestion.
- **A validated common trunk**: sound-speed correction (TEOS-10), composable QC,
  IGRF declination, exact 3-beam beam→XYZ→ENU transform, isobaric regridding, and a
  range-dependent ("shear") bias calibration.
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

adcp  = load_ad2cp("sea064_M38.ad2cp")            # raw binary (or MIDAS .nc, or $PNOR)
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
how to judge the results. The design record (plan, verified literature, analyses of the
Python reference implementations, format specifications, and the M38 validation report)
lives in the repository under `PLAN.md` and `docs/research/`.

## Validation summary

| check | result (sea064 M38 reference mission) |
|---|---|
| native binary reader vs MIDAS netCDF | bit-identical (124,752 ensembles + bottom track) |
| transform parity vs Python `gliderad2cp` | machine-exact (max Δ 2×10⁻¹⁶ m/s) |
| synthetic truth (both solvers, incl. bottom-track-only referencing) | recovered within bin discretization |
| dive vs climb consistency (inverse) | r = 0.98, median \|Δ\| = 2 cm/s |
| glider velocity vs *unseen* bottom track | r = 0.97, median \|Δ\| ≈ 7 cm/s |
| shallow bins vs surface GPS drift | median \|Δ\| = 4 cm/s |
| DAC closure | median 5 mm/s |
