# GliderADCP.jl — Implementation Plan

**A pure-Julia toolbox for processing Nortek AD2CP data from SeaExplorer and Slocum gliders
into absolute ocean velocity profiles.**

Status: research phase complete (2026-07-06); package scaffold in place; implementation
proceeds one layer at a time per the roadmap in §8.

Supporting research (all verified, on disk):

- [docs/reference_dataset.md](docs/reference_dataset.md) — sea064 M38 reference dataset inventory (first-hand)
- [docs/research/gliderad2cp_analysis.md](docs/research/gliderad2cp_analysis.md) — full pipeline spec of the Python shear implementation (+ run-captured ground truth)
- [docs/research/slocum_ad2cp_analysis.md](docs/research/slocum_ad2cp_analysis.md) — full pipeline spec of the Python inverse implementation
- [docs/research/literature.md](docs/research/literature.md) — verified annotated bibliography + synthesis
- [docs/research/formats_and_ecosystem.md](docs/research/formats_and_ecosystem.md) — .ad2cp binary format spec, MIDAS conventions, Julia ecosystem scan

---

## 1. Goal and guiding principles

Replace the current Python step in the user's glider workflow (jlglider's 2,900-line
`process_adcp.py`, a gliderad2cp derivative) with a documented, tested, registered Julia
package that:

1. **Implements multiple methods over one common trunk** — the lADCP-style *shear method*
   (gliderad2cp lineage) and the Visbeck *least-squares inverse* (Todd 2017 / Gradone 2023
   lineage) consume identical QC'd inputs, so methods can be compared apples-to-apples on
   the same mission. Neither Python package offers this today (gliderad2cp has no inverse;
   Slocum-AD2CP's shear method is broken).
2. **Builds from first principles** — every transform derives from the measurement model in
   §2, validated against the instrument's own factory matrices and reference outputs, not
   copied blind.
3. **Fixes known defects consciously** — both Python packages contain verified bugs (§3);
   we implement the physically correct behavior by default and keep *parity switches* so
   regression tests can still reproduce Python outputs bit-approximately.
4. **One layer at a time** — each phase lands a self-contained, tested layer with explicit
   acceptance criteria against real data (§8).
5. **Julia-native, no Python at runtime** — PythonCall appears only in gated
   cross-validation tests.

## 2. First principles

### 2.1 Measurement model

Each beam *i* measures the along-beam Doppler velocity of scatterers relative to the moving
glider:

```
b_i = e_i(t) · (u_ocean(z_cell) − u_glider(t)) + noise
```

where `e_i` is the beam unit vector in the earth frame (beam geometry rotated by heading/
pitch/roll). Per ping, 3 beams ⇒ 3 components of the **relative** velocity
`v_rel = u_ocean − u_glider` at each cell. The ocean and glider contributions are not
separable from one ping — that is the core underdetermination. Closures:

- **Vertical differencing (shear method):** `∂v_rel/∂z = ∂u_ocean/∂z` because `u_glider`
  is constant across cells within a ping. Integrate shear → baroclinic profile; add one
  reference velocity.
- **Joint inversion (inverse method):** solve for `u_ocean(z bins)` and `u_glider(pings)`
  simultaneously with constraints (§5.2).
- **References available:** DAC from navigation (GPS-fix minus dead-reckoned displacement
  over the submerged interval; accuracy 1–2 cm/s RMS, Gradone 2023 / Rudnick 2018);
  **bottom track** (direct `u_glider` over ground when in range — enabled in the M38
  reference data); surface GPS drift; flight-model through-water velocity (Merckelbach
  2010/2019; Todd 2017 uses it as an extra constraint).

### 2.2 Glider AD2CP geometry (verified from instrument config + Nortek docs)

Beams: 1 fwd θ=47.5°, 2 stbd θ=25°, 3 aft θ=47.5°, 4 port θ=25° (φ = 0/−90/180/+90).
Design intent: at flight pitch ≈ ∓17.5°, the fore (dive) or aft (climb) 47.5° beam plus
both 25° side beams all sit ≈30° from vertical — an effective symmetric 3-beam Janus array.

- **Beam selection:** dive (pitch<0): beams {1,2,4}; climb: beams {2,3,4}
  (Gradone 2023; mirrored appropriately for mounting orientation — gliderad2cp instead
  *synthesizes* the dropped beam via the zero-error-velocity condition, which is
  algebraically equivalent; we implement the direct 3×3 solve and test the equivalence).
- **beam→XYZ:** factory matrix from the netCDF `Config` group (rows X, Y, Z1, Z2;
  e.g. ±0.6782 = 1/(2 sin 47.5°), ±1.1831 = 1/(2 sin 25°), 0.74, 0.5518); analytic fallback.
- **XYZ→ENU:** `R = H(heading−90°)·P(pitch, roll)` with the Nortek matrix convention (both
  Python packages agree on P; validated numerically in tests). Down-looking mounts flip the
  Y and Z rows. Magnetic declination added to heading — **not set onboard in the reference
  mission**, so the package must supply it (user value or IGRF/WMM model).
- **Range gating (critical, undocumented in files):** Nortek firmware time-gates all four
  beams with a hard-coded nominal 25° slant and 1500 m/s when >2 beams ping: along-beam
  cell distance = reported `Velocity Range` / cos 25°. Vertical offset of a cell then uses
  the *actual* per-beam angle from vertical:
  `θ_eff = arccos(cos(θ_b ∓ pitch or roll)·cos(roll or pitch))` (per-beam formulas in the
  research docs). Slocum-AD2CP omits the 1/cos 25° factor — a documented discrepancy we
  resolve in favor of the Nortek-confirmed behavior (parity switch for regression).
- **Sound speed:** beam velocities scale linearly: `v_true = v_recorded · c_true/c_used`
  (Fischer & Visbeck 1993). M38 recorded with salinity misconfigured (38 vs ~35.05) —
  c error ≈ 4 m/s ≈ 0.3 % velocity bias; correct from CTD via TEOS-10 (GibbsSeaWater).

### 2.3 Platform data (SeaExplorer first)

`gli` navigation files carry `DeadReckoning` flag, NMEA positions, attitude, and
`Declination` column; DAC per dive = (first GPS fix − last DR position)/submerged duration.
Payload (`pld1`) carries CTD (sound speed, and density for flight model) and a real-time
AD2CP subset; the `$PNORI/S/C` ASCII stream duplicates configuration + per-cell data in
real time. Slocum ingestion (via the user's JLDBDReader.jl: `m_water_vx/vy`,
`m_gps_mag_var`, etc.) comes after the SeaExplorer path is validated.

## 3. What the Python implementations teach us

Full specs in the research docs. Decisions they drive:

| Aspect | gliderad2cp (shear) | Slocum-AD2CP (inverse) | GliderADCP.jl |
|---|---|---|---|
| QC defaults | corr≥80 %, amp≤80 dB & ≥floor+3 dB, \|v\|≤0.8 m/s | corr≥50 %, amp≤75 dB, \|ENU\|≤0.75, first bin, <5 m | all screens available, defaults per-mode, + sidelobe/ambiguity/BT-FOM (§6) |
| Declination | none (only manual heading_offset) | rotate DAC + heading by m_gps_mag_var | first-class: user value or IGRF, applied to heading and any magnetic-frame DAC |
| Cell depths | Range/cos 25° then per-beam θ_eff | per-beam θ_eff without the 25° un-projection | Nortek-confirmed behavior; parity switches |
| Regridding | isobaric offsets at cellSize/2 | interp onto VelocityRange as offsets | isobaric offset grid (explicit spacing & sign) |
| Reference | DAC via ADCP-as-DVL displacement | glider's own m_water_vx/vy | DAC from nav (SeaExplorer DR/GPS), DVL-style optional, BT, surface drift |
| Known defects | per-index shear scaling (exact only for down-looking 2-m cells), left-edge grid loss, bias-mask no-op | dive-matrix overwrite, rotation transpose, off-by-one bin drop, end-time as midpoint, index-aligned deployment grid | physically correct defaults + documented parity modes |
| Solver | cumsum integration | scipy lsqr, complex U+iV | explicit trapezoidal/cumsum (shear); Krylov lsqr/lsmr or SPQR (inverse), complex or stacked-real |

Both packages assume MIDAS netCDF input, use only 3 beams per cast direction, ignore W in
the horizontal solutions, and lack bottom-track support — the M38 dataset has BT data, so
we implement the BT constraint from day one (novel relative to both).

## 4. Package design

### 4.1 Identity & integration

- Module **`GliderADCP`**, repo to be renamed `GliderADCP.jl` at registration time.
  Julia ≥ 1.10 (LTS). MIT license (already present).
- Sibling, not monolith: Slocum I/O delegates to **JLDBDReader.jl**; SeaExplorer parsing
  conventions harvested from **jlglider** (eventually jlglider calls this package);
  package conventions modeled on **ATOMIXjulia.jl**.
- Deps (already resolved in Project.toml): NCDatasets, CSV, DataFrames, CodecZlib,
  Interpolations, GibbsSeaWater, NaNStatistics, Krylov + stdlibs. CairoMakie later as a
  package extension. PythonCall only under `test/` (CondaPkg.toml there), gated by
  `ADCP_CROSSVAL=1`.

### 4.2 Conventions (fixed now, documented in code)

- Depth positive down (m); cell *offsets* signed positive in the look direction.
- Velocities m/s, ENU with U(p) positive up; angles degrees at API boundaries, radians inside.
- Missing data = `NaN` inside numeric arrays (performance; `missing` only at I/O edges).
- Time: `DateTime` (UTC) in structs; Float64 Unix seconds for interpolation internals.
- Every product carries provenance attributes (package version, options used, input files,
  QC-rejection statistics).

### 4.3 Architecture (matches `src/` scaffold)

```
Layer 0  types.jl                AD2CPConfig, AD2CPData, BottomTrackData, GliderNav,
                                 Dac, ProcessingOptions
Layer 1  io/nortek_netcdf.jl     load_ad2cp(paths): Config + Data/Average + Data/AverageBT
         io/seaexplorer.jl       load_gli / load_pld1 (gz, NMEA coords, DR flag)
         io/slocum.jl            via JLDBDReader.jl / CSV exports
         io/nortek_pnor.jl       $PNOR real-time stream
         io/ad2cp_binary.jl      native reader (later; spec in research docs)
Layer 2  processing/soundspeed.jl  correct_soundspeed!(adcp, ctd)
         processing/qc.jl          composable masks + qc report
         processing/geometry.jl    beam_unit_vectors, select_beams, beam2xyz, xyz2enu,
                                   declination handling
         processing/binmap.jl      cell vertical offsets, isobaric regridding
Layer 3  processing/dac.jl         dac_from_seaexplorer(nav), dac_from_slocum,
                                   dac_dvl(adcp) [gliderad2cp-style], surface_drift
Layer 4  solutions/shear.jl        grid_shear → integrate → reference (DAC/BT), bias check
         solutions/inverse.jl      build_inverse(...; constraints=[Dac(), BottomTrack(),
                                   Smoothness(), SurfaceDrift()]) → solve
Layer 5  products/grid.jl, export.jl, Makie extension
```

Public API sketch (final deliverable for a mission):

```julia
adcp   = load_ad2cp("sea064_M38.ad2cp.*.nc")
nav    = load_seaexplorer_nav("delayed/nav/logs")           # gli files
ctd    = load_seaexplorer_pld("delayed/pld1/logs")          # legato/pld1
opts   = ProcessingOptions(declination=IGRF(), qc=QCDefaults(:average_mode))
pings  = process_pings(adcp, nav, ctd, opts)                # layers 2–3 trunk
dac    = compute_dac(nav)
prof_s = solve_shear(pings, dac, opts)
prof_i = solve_inverse(pings, dac, opts;
                       constraints=(DacConstraint(5.0), Smoothness(1.0),
                                    BottomTrackConstraint(adcp.bt)))
compare(prof_s, prof_i)                                     # method intercomparison report
```

## 5. The methods

### 5.1 Method A — shear (gliderad2cp lineage, corrected)

1. Per-ping ENU relative velocities on an isobaric offset grid (Layer 2).
2. Shear per ping: first differences over *physical* Δz of the offset grid (explicit signed
   spacing — fixes the per-index scaling defect; `parity=:gliderad2cp` reproduces it).
3. Grid shear into (depth × profile) bins: median (default) with count/std diagnostics;
   inclusive right-edge handling (fixes grid2d edge loss).
4. Integrate top-down per profile → zero-mean baroclinic profile.
5. Reference: time-in-bin-weighted DAC referencing (gliderad2cp semantics), or plain
   depth-average matching; optional near-bottom BT referencing where BT is valid.
6. Optional shear-bias diagnosis (regression of velocity-gradient vs through-water
   displacement, gliderad2cp process_bias semantics) — reported, opt-in correction.

### 5.2 Method B — inverse (Visbeck 2002; Todd 2017; Gradone 2023, corrected)

Per segment (dive, climb, yo, or arbitrary window), complex `m = [u_g(1..nt); u_o(1..nz)]`:

- Measurement rows: `u_o(bin of cell) − u_g(ping) = v_rel` — nearest-bin (+1/−1) by default,
  optional linear interpolation weights split across adjacent bins.
- Constraint blocks (composable, weighted, all optional):
  - `DacConstraint`: ocean-side (Σ dz·u_o/H = DAC, Gradone form) **or** platform-side
    (Σ u_g·Δt = GPS displacement, Visbeck form) — both provided; weight default 5.
  - `BottomTrackConstraint`: rows `u_g(ping) = −v_BT(ping)` for FOM-valid BT pings (strong
    weight) — new relative to both Python packages.
  - `SurfaceDriftConstraint`: near-surface `u_o` from GPS drift between fixes (optional).
  - `Smoothness`: second-difference on the ocean block (and optionally the glider block,
    sized correctly — fixes the nz×nt truncation defect); weight default 1.
- Solve: `Krylov.lsqr` (damped optional) or sparse QR; U+iV jointly (complex) with a
  stacked-real fallback. Report per-bin observation counts and formal uncertainties from
  the normal-equation diagonal (χ² scaled).
- Segment modes: per-yo (Gradone parity), per-cast, and whole-mission joint solve
  (Julia's sparse solvers make the ~10⁵–10⁶-row mission-scale problem tractable).

### 5.3 Method C — extensions (post-v1)

Flight-model through-water velocity (Merckelbach 2010/2019) as an additional constraint or
independent estimate + w_water product; ADCP-nearest-bin flight-model calibration
(Welch 2022); compass-error diagnosis (von Appen 2015 / pitch-binned ellipsoid fits as in
Slocum-AD2CP); EVR/phase-unwrapping for burst modes (Shcherbina 2018); DIVAnd section mapping.

## 6. QC catalogue (Layer 2; all composable masks with per-screen rejection stats)

| Screen | Default | Source |
|---|---|---|
| Correlation | ≥ 50 % (average mode); 80 % strict preset | Gradone 2023 / gliderad2cp |
| Amplitude high | ≤ 75 dB | Gradone 2023 |
| Amplitude low / SNR | ≥ noise floor (0.5th pct pooled) + 3 dB | gliderad2cp |
| Relative speed | \|v_beam\| ≤ 0.8 m/s (and \|ENU\| post-transform) | both |
| Ambiguity proximity | flag \|v_beam\| > 0.9·V_R (V_R from Config, 2.5 m/s for M38) | Nortek/Shcherbina |
| Side lobe | cut range > D_boundary·cos θ_beam per beam, using altimeter/BT distance or surface | RDI primer |
| Near-surface | glider depth ≤ 5 m mask (configurable) | Slocum-AD2CP |
| First cell(s) | drop n cells nearest transducer (default 1) | Slocum-AD2CP |
| Attitude | \|pitch−pitch_flight\|, \|roll\| limits (optional) | community practice |
| Status/Error bits | instrument error ≠ 0 → drop ping | Nortek |
| BT validity | FOM ≠ 65535 and distance in (blank, range) | Nortek BT record |

## 7. Validation strategy (multi-level, per layer)

1. **Unit/synthetic (every phase):** synthetic beam velocities from a prescribed
   `u_ocean(z)` + glider trajectory → recover inputs through each layer (rotations invert,
   shear integrates back, inverse recovers truth within noise); property tests for the
   3-beam ↔ beam-synthesis equivalence; checksum/format tests for readers.
2. **Instrument self-consistency (M38):** factory beam2xyz vs analytic; BT-derived
   `u_glider` vs inverse-estimated `u_g` on bottom-track segments; dive vs climb profile
   agreement; DAC closure (mean of solution vs nav DAC).
3. **Regression vs prior outputs (M38, `m38_processed/`):** ENU stage vs
   `M38_ADCP_QAQC_CoordTransformed.nc` (VelE/VelN/VelU; parity mode with Gradone
   conventions); per-yo inverse vs `absolute_ocean_vel.csv` (expect agreement up to their
   documented bugs; quantify each deviation).
4. **Cross-validation vs gliderad2cp (SEA055 sample, `validation/gliderad2cp_reference/`):**
   run-captured ground truth (DAC[1] = [−0.02708, +0.10196] m/s, etc.); optional live
   PythonCall comparison behind `ADCP_CROSSVAL=1`.
5. **Method intercomparison (the science check):** shear vs inverse vs BT on M38; target
   the literature accuracy class (few cm/s; DAC floor 1–2 cm/s).

## 8. Roadmap — one layer at a time

| Phase | Deliverable | Acceptance criteria |
|---|---|---|
| **0 — done** | research docs, dataset inventory, scaffold (compiles, tests pass) | this document |
| **1 — done (2026-07-06)** | `load_ad2cp` (Average + AverageBT + Config), SeaExplorer gli/pld parsers, core types | ✔ loads M38 (124,752 ensembles + 124,751 BT recs) & M48 sample; config matches docs/reference_dataset.md; nav parse: 103,150 records, 190 DR→GPS surfacing transitions (~191 yos) |
| **2 — core done (2026-07-06)** | soundspeed, QC masks, beam selection, beam→XYZ→ENU (exact 3-beam solve + `:rowdrop` parity), vertical cosines, isobaric regrid | ✔ synthetic round-trips exact (rotation/declination/4-beam LS); ✔ factory X,Y rows reproduced analytically; ✔ M38 smoke (63.7k ENU samples, med \|v_rel\| = 0.26 m/s). **Pending:** parity regression vs `M38_ADCP_QAQC_CoordTransformed.nc` (needs Gradone-convention parity mode) and vs gliderad2cp SEA055 ground truth |
| **3 — DAC** | DAC from SeaExplorer nav, DVL-style DAC, surface drift | M38 per-yo DAC series sane vs `M38_glider_data_processed.csv`; SEA055 DAC matches gliderad2cp ground truth |
| **4 — solutions** | shear + inverse with constraint blocks (incl. BT) | synthetic truth recovered; per-yo inverse vs `absolute_ocean_vel.csv`; dive/climb & DAC-closure consistency; method intercomparison report for M38 |
| **5 — products & platforms** | gridding/export, Makie ext, $PNOR parser, Slocum path via JLDBDReader, M38 example script, Documenter docs | end-to-end M38 script reproduces (improved) U/V section figures; docs build |
| **6 — native reader & release** | .ad2cp binary reader; CI, Aqua, registration prep | binary reader matches MIDAS netCDF on M38 (velocities bit-scaled identical); General registry checklist |
| **7 — extensions** | flight model, compass cal, burst/EVR, DIVAnd sections | per-feature |

Each phase = one PR-sized unit with its tests; no phase starts until the previous one's
acceptance criteria pass.

## 9. Risks & open questions

- **Attitude/sign conventions** are the classic failure mode: resolved by triple-checking
  every rotation against (i) the factory matrix, (ii) both Python parity modes, and
  (iii) BT ground truth (BT velocity must oppose glider motion).
- **Nortek 25° range-gating** is firmware-dependent and undocumented in files (pers.-comm.
  provenance) — verify empirically on M38 (bottom-distance vs BT distance consistency) and
  make the assumption explicit in options.
- **Declination model**: add IGRF via a small dependency (e.g. SatelliteToolboxGeomagneticField.jl)
  or accept user-supplied values initially — decide in Phase 2 (start user-supplied).
- **gliderad2cp JOSS review is ongoing** — track for algorithm changes; our parity targets
  pin to commit `7fbccaf`.
- **Slocum-AD2CP parity vs correctness**: m38_processed outputs embed that package's bugs;
  regression comparisons must model those deviations explicitly rather than chase equality.
- **Compass calibration correction** (beyond declination) is the largest unmodeled error
  term (von Appen 2015); Phase 7.

## 10. Key references

Visbeck 2002 (10.1175/1520-0426(2002)019<0794:DVPULA>2.0.CO;2) · Fischer & Visbeck 1993 ·
Todd et al. 2017 (10.1175/JTECH-D-16-0156.1) · Todd et al. 2011 (10.1029/2010JC006536) ·
Gradone et al. 2023 (10.1029/2022JC019608) · Queste et al., gliderad2cp (JOSS in review,
10.21105/joss.08342) · Thurnherr et al. 2015 (10.1109/CWTM.2015.7098134) · GO-SHIP LADCP
cookbook (Thurnherr et al. 2010) · Merckelbach et al. 2010 (10.1175/2009JTECHO710.1) &
2019 (10.1175/JTECH-D-18-0168.1) · Rudnick et al. 2018 (10.1175/JTECH-D-17-0200.1) ·
Welch et al. 2022 (10.1175/JTECH-D-21-0074.1) · Shcherbina et al. 2018
(10.1175/JTECH-D-17-0108.1) · von Appen 2015 (10.1175/JTECH-D-14-00043.1) · Ma et al. 2019
(10.2112/JCOASTRES-D-18-00176.1) · Nortek N3015-011 Principles of Operation & N3015-007
Integrators Guide · Teledyne RDI ADCP Primer. Full annotations: docs/research/literature.md.
