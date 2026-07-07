# GliderADCP.jl

Pure-Julia processing of glider-mounted Nortek AD2CP data (SeaExplorer, Slocum) into
absolute ocean velocity profiles, implementing both the lADCP-style **shear method**
and the Visbeck **least-squares inverse** (with DAC, bottom-track and smoothness
constraints) over one common trunk.

```julia
using GliderADCP
adcp  = load_ad2cp("seaXXX_MYY.ad2cp.00000.nc")   # Nortek MIDAS export
nav   = load_seaexplorer_nav("delayed/nav/logs")
qc!(adcp)
pings = process_pings(adcp; lat=69.0,
                      declination=magnetic_declination(nav, adcp.t))
dac   = compute_dac(nav)
btv   = bt_velocity(adcp)
prof  = solve_inverse(pings, dac; bt=btv)
sec   = grid_profiles(prof)
export_sections("sections.nc", sec)
```

See `PLAN.md` and `docs/research/` in the repository for the design rationale,
method derivations, verified literature, and validation reports.
