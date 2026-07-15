# CLAUDE.md — working notes for this repository

GliderADCP.jl: pure-Julia processing of glider-mounted ADCP data (currently the
Nortek AD2CP) into absolute ocean velocities. Two solvers (Visbeck inverse + shear method) over one common trunk.
Companion packages: SeaExplorerIO.jl (shared file layer; resolved from GitHub via
`[sources]` url — for local loader development, `Pkg.develop(path="../SeaExplorerIO.jl")`
into the environment you're testing, and push SeaExplorerIO before testing here) and
GliderTurbulence.jl (microstructure; consumes this package's calibration through a package
extension; not yet public).

## Read before re-investigating anything

- `docs/research/m38_validation.md` — the cumulative validation/finding log. Every
  major claim, retraction, and verdict is recorded there with evidence.
- `docs/src/qaqc.md` — all data QA/QC findings as *what we found → what to check*.
- `PLAN.md` — design record; §8a/§8b list completed follow-up work.

## Commands

```bash
# tests (413; gated acceptance tests auto-skip when local mission data is absent)
~/.juliaup/bin/julia +1.13 --project=. -e 'using Pkg; Pkg.test()'

# docs (Documenter; SeaExplorerIO is dev'd into docs/Manifest)
JULIA_LOAD_PATH="@:@ocean:@stdlib" ~/.juliaup/bin/julia +1.13 --project=docs docs/make.jl

# examples (CairoMakie comes from the user's @ocean env via the stacked load path;
# the package itself stays plot-free). All are driven by examples/missions.jl —
# add a mission there once, every script picks it up. No args = all missions.
JULIA_LOAD_PATH="@:@ocean:@stdlib" ~/.juliaup/bin/julia +1.13 --project=. examples/currents.jl [m38 ...]
JULIA_LOAD_PATH="@:@ocean:@stdlib" ~/.juliaup/bin/julia +1.13 --project=. examples/realtime_onboard.jl
JULIA_LOAD_PATH="@:@ocean:@stdlib" ~/.juliaup/bin/julia +1.13 --project=. examples/realtime_telemetered.jl
JULIA_LOAD_PATH="@:@ocean:@stdlib" ~/.juliaup/bin/julia +1.13 --project=. examples/dac_methods.jl
```

Reference mission data lives under
`/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/` (four validated missions:
M37, M38, M48, M59 — paths in `examples/missions.jl`). Watch for unsynced Dropbox
placeholder files (read as zero bytes; loaders warn `no rows parsed`).

## Standing decisions (do not silently revisit)

- **Data-route taxonomy** (user-set, use these names): **delayed-mode** (`.ad2cp`
  binary, `read_ad2cp` — the reference), **realtime-onboard** (`$PNOR` stream,
  `load_pnor` — payload-logged; in real time useful only to an onboard consumer such
  as a backseat driver), **realtime-telemetered** (`pld1.sub` AD2CP subset,
  `load_pld_adcp` — the only tier ashore mid-mission; **shore-side realtime products
  build on this route**). ALSEAMAR's GLIMPSE `AD2CP_*_c` product is computed
  server-side from the same telemetered data.
- **The inverse is the production method**; the shear method is the standard second
  opinion. Their agreement (r = 0.90–0.98 across missions) is the top health metric —
  a collapse flags contamination. Full argument + limits: validation doc §Method verdict.
- **First cell is kept by default** (`QCThresholds.first_cells = 0`): the fleet's
  deliberate 0.7 m blanking makes cell 1 clean (validated on all four missions).
  `qc!` warns if blanking < 0.5 m — small-blanking configs need `first_cells = 1`.
- **Never trust bottom track unscreened**: 99.7 % of M38's BT locks were false
  near-field water-borne targets. `bt_valid` defaults (min_range = 5 m +
  impossible-bathymetry test) stay on; verify surviving locks' implied water depth
  against bathymetry.
- **DAC is a water-track ladder** (`compute_dac(nav, pings; fallback=flight_model(nav))`):
  ALSEAMAR's onboard DR flight model runs ×1.05–×1.15 fast (all four missions),
  biasing the nav-only DAC 2–4 cm/s anti-track. Rungs, flagged per yo in `method`:
  `:adcp` (integrated `throughwater_velocity`) → `:flight` (our flight model,
  ~1.4 cm/s from the water track, unbiased; also the ADCP-less form
  `compute_dac(nav, flight_model(nav))`) → `:onboard` (last resort).
  Realtime-telemetered pings need `max_gap=90.0` (~30 s cadence). Drift-verified on
  M37/M38/M48; M59's drift referee is windage-biased mission-wide — do not "fix" the
  M59 acid test. Evidence: validation doc 2026-07-15 entries, QA/QC §3b.
- **The flight model is deliberately twinned** in GliderADCP
  (`src/processing/flightmodel.jl`) and GliderTurbulence (`src/flightmodel.jl`) so
  each package stands alone — do not deduplicate; land physics fixes in BOTH files.
  Loading both packages requires qualifying the shared names (`flight_model`,
  `FlightParams`, `solve_aoa`, `measure_aoa`, `fit_flightparams`, `FLIGHT_*`);
  GliderTurbulence's ADCP extension already imports its own explicitly.
- **gli `Heading` is TRUE heading** (declination applied onboard; verified vs the
  AD2CP magnetic compass). Never add declination to nav heading; the AD2CP's own
  heading IS magnetic and `process_pings` adds declination there.
- **Shear-bias calibration is per-mission** (`calibrate_shear_bias!`): the slope is
  configuration-dependent (−4.7×10⁻⁴ … −5×10⁻⁵ s⁻¹ across 2022–2024), never hard-code.
- **Telemetered w: events yes, statistics no.** The 30-s subsampling aliases texture
  (r = 0.66–0.84 vs delayed) but large coherent vertical-velocity events survive —
  the per-mission diagnostic `M*_telemetered_w_sections.png` (from
  `examples/realtime_telemetered.jl`) is the check; wave-scale w work needs the
  `$PNOR` (onboard) or delayed tiers.
- **NORSE** (not NorSE) in all mission references.

## Conventions

- Layered architecture (io → processing → solutions → products); each layer small and
  testable; parity switches (`:rowdrop` etc.) preserve regression against the Python
  references rather than reproducing their bugs by default.
- Acceptance tests against real missions are gated on local data existence — add one
  when landing a data-facing feature. Loaders must *degrade loudly* (skip corrupt
  files with warnings, report gaps/zero-row sources), never crash or stay silent.
- `examples/output/` is gitignored; figures regenerate from the example scripts.
- Julia compat is 1.11+ (the `[sources]` table needs it); SeaExplorerIO resolves from
  its public GitHub url, so fresh clones and CI install with no extra steps.
