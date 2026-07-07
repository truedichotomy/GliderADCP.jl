# GliderADCP.jl

Pure-Julia processing of glider-mounted Nortek AD2CP data (SeaExplorer, Slocum)
into absolute ocean velocity profiles.

Implements both published approaches — the lADCP-style **shear method**
(Todd et al. 2017; cf. Python [`gliderad2cp`](https://github.com/bastienqueste/gliderad2cp))
and the Visbeck (2002) **least-squares inverse method**
(Todd et al. 2017; Gradone et al. 2023; cf. [`Slocum-AD2CP`](https://github.com/JGradone/Slocum-AD2CP)) —
plus bottom-track and surface-drift constraints, built from first principles in
independent layers.

**Status:** Phases 1–6 implemented and tested (243 tests, including acceptance runs
against a full SeaExplorer mission): I/O (**native `.ad2cp` binary reader — bit-identical
to the MIDAS export, no MIDAS/Windows needed** — plus MIDAS netCDF incl. bottom track,
SeaExplorer gli/pld parsers, the `\$PNOR` real-time stream, and Slocum tables),
sound-speed correction, QC, the exact 3-beam beam→XYZ→ENU transform,
isobaric regridding, DAC + surface drift from navigation, and **both velocity solvers**
(shear and the Visbeck-style inverse with composable DAC / bottom-track / smoothness
constraints). Validated three ways: machine-exact parity against `gliderad2cp` ground
truth; synthetic-truth recovery for both solvers; and an independent real-data check —
the DAC-only inverse's glider velocities match unseen bottom-track ground truth at
r = 0.97 (median difference ≈ 7 cm/s). A full-mission example
([examples/m38_currents.jl](examples/m38_currents.jl)) produces DAC+bottom-track
inverse and shear U/V sections with IGRF declination, provenance netCDF export, and
dive/climb consistency of 2 cm/s; the vertical-structure validation against the prior
Python processing is documented in
[docs/research/m38_validation.md](docs/research/m38_validation.md).
See [PLAN.md](PLAN.md) for the roadmap.

```julia
using GliderADCP
adcp = load_ad2cp("sea064_M38.ad2cp.00000.nc")   # Data/Average + Data/AverageBT + Config
nav  = load_seaexplorer_nav("delayed/nav/logs")   # gli files → GPS/DR segments
qc!(adcp)                                         # correlation/amplitude/SNR/… masks
pings = process_pings(adcp; lat=69.0)             # ENU relative velocities on isobars
dac   = compute_dac(nav)                          # per-yo depth-averaged current
btv   = bt_velocity(adcp)                         # over-ground velocity when in BT range
prof_i = solve_inverse(pings, dac; bt=btv)        # absolute velocity profiles (inverse)
prof_s = solve_shear(pings, dac)                  # absolute velocity profiles (shear)
```
