# GliderADCP.jl

Pure-Julia processing of glider-mounted Nortek AD2CP data (SeaExplorer, Slocum)
into absolute ocean velocity profiles.

Implements both published approaches — the lADCP-style **shear method**
(Todd et al. 2017; cf. Python [`gliderad2cp`](https://github.com/bastienqueste/gliderad2cp))
and the Visbeck (2002) **least-squares inverse method**
(Todd et al. 2017; Gradone et al. 2023; cf. [`Slocum-AD2CP`](https://github.com/JGradone/Slocum-AD2CP)) —
plus bottom-track and surface-drift constraints, built from first principles in
independent layers.

**Status:** Phase 1 (I/O: MIDAS netCDF incl. bottom track, SeaExplorer gli/pld parsers)
and the Phase 2 core (sound-speed correction, QC, exact 3-beam beam→XYZ→ENU transform,
isobaric regridding) are implemented and tested (112 tests, including acceptance runs
against a full SeaExplorer mission). See [PLAN.md](PLAN.md) for the research-backed
roadmap and [docs/reference_dataset.md](docs/reference_dataset.md) for the reference
dataset.

```julia
using GliderADCP
adcp = load_ad2cp("sea064_M38.ad2cp.00000.nc")     # Data/Average + Data/AverageBT + Config
nav  = load_seaexplorer_nav("delayed/nav/logs")     # gli files → GPS/DR segments
qc!(adcp)                                           # correlation/amplitude/SNR/… masks
E, N, U, offsets, beams = enu_on_isobars(adcp)      # relative velocities, earth frame
```
