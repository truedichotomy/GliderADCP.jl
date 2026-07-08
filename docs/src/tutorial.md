# Tutorial: from raw AD2CP data to ocean velocity

This tutorial walks through the full processing chain on a real mission — sea064 M38
(SeaExplorer, Nortek Glider AD2CP 1 MHz, Lofoten Basin, Nov 2022) — explaining the
science behind each step, the choices you can make, and how to judge the results.
The complete runnable version is
[`examples/m38_currents.jl`](https://github.com/truedichotomy/gliderADCPjulia/blob/main/examples/m38_currents.jl).

## 1. The measurement problem

Each ADCP beam measures the along-beam Doppler velocity of scatterers **relative to the
moving glider**:

```math
b_i = \mathbf{e}_i \cdot \left(\mathbf{u}_{ocean}(z_{cell}) - \mathbf{u}_{glider}(t)\right) + \varepsilon
```

Three beams per ping give the vector relative velocity at each cell, but the ocean and
glider contributions are inseparable from any single ping — that is the core
underdetermination of glider ADCP work. Something extra must close the system:

* **Depth-averaged current (DAC)** — the glider dead-reckons underwater; the position
  jump at the first GPS fix after surfacing, divided by the submerged time, is the
  time-averaged current over the yo (accuracy 1–2 cm/s RMS; Rudnick et al. 2018,
  Gradone et al. 2023).
* **Bottom track** — when the seafloor is within range, the instrument measures the
  glider's velocity over ground directly: an absolute reference with no GPS involved.
* **Surface GPS drift** — consecutive fixes while surfaced give the near-surface
  velocity.

Two published strategies use these closures, and this package implements both over one
common trunk so they can be compared on identical inputs:

* the **shear method** (Fischer & Visbeck 1993; Thurnherr et al. 2015; the Python
  `gliderad2cp`): vertically difference the relative velocities (glider motion cancels
  within a ping), bin and integrate the shear, then shift the profile to match the DAC;
* the **inverse method** (Visbeck 2002; Todd et al. 2017; Gradone et al. 2023): solve
  for ocean velocity per depth bin *and* glider velocity per ping simultaneously in one
  weighted least-squares system with DAC/bottom-track/smoothness constraint rows.

Section 7 discusses why their error behavior differs fundamentally.

### The glider AD2CP geometry

The Glider AD2CP has four beams: fore and aft at 47.5° from the instrument axis, port
and starboard at 25°. The design exploits the glider's flight attitude: pitched at
≈ ∓17.5°, the fore (dive) or aft (climb) beam plus both side beams all sit ≈30° from
vertical — an effective symmetric Janus array. Processing therefore uses **three beams
per ping** (dive: 1, 2, 4; climb: 2, 3, 4 for a down-looking mount) and solves the 3×3
system exactly. Two verified geometry facts are built in (see
`docs/research/` for the evidence):

* the exact 3-beam solve equals gliderad2cp's "synthesize the dropped beam" trick and
  is *not* the same as slicing the factory 4-beam matrix (which halves the along-glider
  component — a documented defect of the Slocum-AD2CP lineage, reproducible here with
  `method = :rowdrop` for regression studies);
* Nortek firmware time-gates **all four beams with a nominal 25° slant**, so the
  along-beam cell distance is `reported range / cos 25°` (`range_gating = :nominal25`,
  the default).

## 2. Loading data

Three AD2CP input routes produce the same `AD2CPData` structure:

```julia
using GliderADCP

# 1. native binary from the instrument SD card — no Nortek MIDAS needed
adcp = load_ad2cp("sea064_M38.ad2cp")

# 2. MIDAS netCDF export (also multi-file: pass a directory or vector)
adcp = load_ad2cp("sea064_M38.ad2cp.00000.nc")

# 3. the real-time $PNOR telemetry stream (SeaExplorer payload logs)
adcp_rt = load_pnor("delayed/pld1/logs")
```

The native reader was validated **bit-for-bit** against the MIDAS export of the same
file (every velocity, amplitude, correlation, attitude and bottom-track sample), and the
`$PNOR` stream reproduces the full-resolution record at its reduced precision
(r = 1.000 for heading and beam velocities). Platform data:

```julia
nav = load_seaexplorer_nav("delayed/nav/logs")   # positions, DR/GPS flags, attitude
pld = load_seaexplorer_pld("delayed/pld1/logs"; stream="pld1.sub")  # CTD for sound speed
```

For Slocum gliders, build the equivalents from any dbd-derived table (ERDDAP, Python
`dbdreader`, or the pure-Julia
[JLDBDReader.jl](https://github.com/truedichotomy/JLDBDReader.jl)) with
[`slocum_nav`](@ref) and [`dac_from_slocum`](@ref).

## 3. Sound-speed correction

Beam velocities scale linearly with the sound speed at the transducer
(Fischer & Visbeck 1993): ``v_{true} = v_{recorded}\, c_{true}/c_{used}``. Instruments
compute ``c_{used}`` from measured temperature and a *configured* salinity — often wrong
(M38 was set to 38 while the Norwegian Sea is ≈35.05, a 0.3 % velocity bias). Correct it
from the glider CTD:

```julia
c = soundspeed_from_ctd.(pld.LEGATO_SALINITY, pld.LEGATO_TEMPERATURE,
                         pld.LEGATO_PRESSURE, 5.0, 69.5)     # TEOS-10
scale = soundspeed_correction(adcp, datetime2unix.(pld.time), c)
apply_soundspeed!(adcp, scale)                # also rescales bottom track
```

0.3 % sounds negligible, but it is a *bias*, not noise: on a 0.3 m/s relative velocity
it is 1 mm/s per sample — the same order as the shear bias discussed below, so it is
worth removing.

## 4. Quality control

```julia
stats = qc!(adcp)          # masks rejected samples to NaN, returns per-screen fractions
```

| screen (default) | removes | why |
|---|---|---|
| correlation ≥ 50 % | decorrelated pings | primary Doppler quality metric (Gradone 2023; gliderad2cp uses 80 for 4-ping averages) |
| amplitude ≤ 75 dB | boundary/fish hits | hard returns are not water velocity |
| amplitude ≥ noise floor + 3 dB | beyond-range cells | pooled 0.5th-percentile floor (gliderad2cp) |
| \|v_beam\| ≤ 0.8 m/s | wraps/outliers | ambiguity velocity of the standard config is 2.5 m/s |
| ambiguity fraction 0.9 | near-wrap samples | uses the configured velocity range |
| glider depth ≤ 5 m | surfaced pings | bubbles, wake, GPS-fix maneuvering |
| first cell | ringing | transducer recovery |
| instrument error ≠ 0 | flagged pings | hardware self-reports |

On M38 the defaults reject ≈50 % of beam samples over the full mission — dominated by
the SNR floor beyond the useful range and the surface mask; the surviving samples are
the ones the solvers should see. Loosening the surface mask to 2 m and keeping the first
cell was tested and does **not** change the near-surface answer (it slightly degrades
the surface-drift agreement) — the defaults are not hiding signal.

### Which cells actually get used?

All configured cells enter, but QC is applied **per sample** (cell × beam × ping), so
cell usage is adaptive: the first cell is dropped unconditionally (transducer ringing),
and far cells fall away wherever the echo decorrelates. Inspect this on your mission
with [`cell_quality`](@ref):

```julia
q = cell_quality(adcp)     # per cell & beam: median corr/amp, survival fractions
```

On M38 (1 MHz, 2-m cells, clear sub-Arctic water) the structure is stark:

| cell (range) | median corr (%) | survives full QC |
|---|---|---|
| 1 (2.7 m) | 97 | 0 % (first-cell drop) |
| 2–6 (4.7–12.7 m) | 82–98 | 91–94 % |
| 7 (14.7 m) | 70 | 76 % |
| 8 (16.7 m) | 55 | 52 % |
| 9 (18.7 m) | 39 | 35 % |
| 12 (24.7 m) | 9 | 17 % |
| 15 (30.7 m) | 2 | 5 % |

The *effective* range here is ~15–17 m (cells 2–8), roughly half the configured 30 m —
what survives further out are the scatterer-rich moments (plankton layers, near
boundaries). Two practical consequences: the `nobs` column in every solver output tells
you how much data supports each depth bin (filter on it), and the marginal far-cell
samples that *barely* pass at 50–70 % correlation carry a disproportionate share of the
range-dependent bias — on M38, raising the correlation threshold from 50 to 80 %
reduces the measured shear-bias slope by ~22 % (−4.0 → −3.1×10⁻⁴ s⁻¹) at the cost of
14 points of retention. For the inverse the default (corr ≥ 50 + calibration) is fine;
for shear-method work, `QCThresholds(correlation=80)` plus the calibration is the
conservative configuration. A future refinement is correlation-based *weighting* of
samples in the inverse (per-sample noise scales as ``\sqrt{R^{-2}-1}``; Shcherbina
et al. 2018) rather than the current pass/fail-plus-equal-weight scheme.

## 5. Declination and the ENU pings

The compass reports magnetic heading; velocities need true east/north. Deployments
usually leave the instrument's declination unset (M38 did), so supply it — a scalar, or
per-ping from IGRF along the track:

```julia
decl  = magnetic_declination(nav, adcp.t)     # M38: −0.1…+4.5 °E along track
pings = process_pings(adcp; lat=69.5, declination=decl)
```

A declination error rotates the velocity vector: 2° on a 0.3 m/s relative velocity is
a 1 cm/s cross-component — at the DAC accuracy floor, so per-ping IGRF is cheap
insurance. `process_pings` performs the full Layer-2 chain: isobaric regridding of each
beam onto common vertical offsets (so the 3-beam combination mixes measurements from a
common depth — gliderad2cp's key insight), the exact 3-beam transform, and absolute cell
depths from TEOS-10. You can verify the compass health first:

```julia
tbl, ptp = compass_field_check(adcp)   # M38: 2.5 % |B| variation — clean
```

## 6. The shear-bias calibration

Glider ADCPs carry a small systematic decay of measured relative velocity with range
from the transducer — on M38, **−4×10⁻⁴ s⁻¹ along-track** (±3.5 mm/s across the 30-m
window; identical on dives and climbs, absent cross-track). It is invisible to any
per-sample QC, and it is precisely a *shear* bias: harmless to the inverse (mm/s per
bin) but integrated by the shear method into a ≈0.2 m/s per 500 m profile tilt
(Todd et al. 2017 and gliderad2cp's `process_bias` address the same phenomenon).

```julia
slopes = calibrate_shear_bias!(pings)   # estimate → subtract; returns slope per pass
```

Details that matter (and are tested):

* the estimator averages the **adjacent-pair differences** in the glider track frame —
  the exact sample population the shear method consumes (a per-offset mean profile
  under-corrects when coverage is partial);
* real ocean shear averages out of the track-frame mission mean **only if headings are
  diverse** — the returned `heading_concentration` warns above R = 0.8;
* the subtracted profile is re-demeaned per ping, so ping means — and therefore the
  inverse's glider-velocity/DAC content — are untouched by construction.

## 7. References and the two solvers

```julia
dac = compute_dac(nav)                        # fix-to-fix DAC per yo (m/s, true E/N)
btv = bt_velocity(adcp; max_range=28.0,
                  declination=magnetic_declination(nav, adcp.bt.t))
drift = surface_drift(nav)                    # near-surface constraint / validation

inv = solve_inverse(pings, dac; bt=btv)       # the reference product
shr = solve_shear(pings, dac)                 # the corrected second opinion
w   = vertical_velocity(pings)                # w = U_rel + dP/dt, flight-model-free
```

!!! warning "False bottom-track locks"
    Glider BT records can be dominated by **false locks on near-field/water-borne
    targets** (wake, scattering layers ~1–3 m below the transducer) whenever the real
    seafloor is out of range. Such targets move with the water; feeding them to the
    inverse as over-ground anchors contradicts the DAC and injects strong spurious
    shear. `bt_valid`/`bt_velocity` therefore screen by default: `min_range = 5 m` and
    an impossible-bathymetry test (implied bottom vs the platform's own deepest nearby
    record). On the reference mission these screens correctly reject *all* 9,383
    apparent fixes — check `nrow(bt_velocity(adcp))` before assuming you have usable
    bottom track.

### The inverse (recommended primary product)

Per fix-to-fix segment, unknowns ``m = [u_g(t_1..t_{nt});\, u_o(z_1..z_{nz})]`` with one
row per QC'd sample, ``u_o(\text{bin}) - u_g(\text{ping}) = v_{rel}``, plus weighted
constraint rows: the DAC (three forms — see [`InverseOptions`](@ref); the
residence-time-weighted `:ocean_timeweighted` form is exact since the DAC is a *time*
average, though on steady-profiling missions it coincides with the plain mean),
bottom-track rows pinning ``u_g`` where the seafloor was in range, and interior
second-difference smoothness. U and V decouple and share one sparse QR factorization.

### The shear method, and why it differs

The shear method is unbiased *per sample* after calibration, but it **integrates**:
per-bin noise random-walks into per-yo excursions of ±0.1–0.3 m/s, whereas the inverse
localizes each sample to its depth bin (per-bin consistency ≈2 cm/s on M38). This is
not an implementation artifact — it is the structural difference between the methods,
and the reason the lADCP community moved to inversions. Use the shear solution as an
independent cross-check, reading its per-yo wiggles with the drift envelope in mind.

Two shear-content products support method intercomparison and finestructure work:
[`solve_shear_profile`](@ref) returns the shear method's pre-integration binned shear,
and [`inverse_shear`](@ref) differences any solver profile. Note that per-yo, bin-scale
shear is internal-wave dominated in most environments (on M38 the direct product's own
dive-vs-climb reproducibility is r = 0.08) — compare shear *statistics* or averaged
profiles, not individual bins, and prefer `inverse_shear` for deterministic sub-inertial
shear (see the validation report for the full analysis).

## 8. Quality metrics to compute on every mission

```julia
# 1. dive vs climb consistency: solve half-yos independently, compare common bins.
#    M38: r = 0.98, median |Δ| = 2 cm/s. Values ≫ 5 cm/s indicate attitude/geometry issues.
# 2. DAC closure: depth-mean of each profile vs its DAC (M38 median 5 mm/s).
# 3. Bottom-track cross-check: run the inverse with wbt = 0 and compare its glider
#    velocities against bt_velocity where the seafloor was in range
#    (M38: r_v = 0.97, median |Δ| ≈ 7 cm/s — an independent end-to-end validation).
# 4. Surface drift: mean of the shallowest bins vs surface_drift after each yo
#    (M38: median |Δ| = 4 cm/s).
```

These four checks catch, respectively: transform/sign errors (1), referencing errors
(2), anything in the whole chain (3), and near-surface problems (4). The package's test
suite contains templates for each.

## 9. Products

```julia
sec = grid_profiles(inv)                        # depth-matched time × depth sections
export_sections("M38_sections.nc", sec;
    attrs=Dict{String,Any}("mission" => "sea064 M38"))

using CairoMakie                                # activates the Makie extension
fig = plot_sections([(sec, :U, "U — inverse"), (sec, :V, "V — inverse")])
```

For continuous (gap-bridged) sections, `examples/m38_divand_sections.jl` maps the
profiles with DIVAnd and masks the analysis where the data give no support (clever
poor man's error > 0.4) — no extrapolation below the sampled envelope.

## 10. Scientific interpretation and caveats

* **Everything is a segment mean.** The inverse, the shear profile, the DAC and the
  drift comparison all average over a yo (2–6 h). Currents with shorter timescales —
  inertial oscillations above all (12.8 h period at 70°N) — vector-average toward zero.
  On M38 the surface drift data show the near-surface flow was genuinely weak and
  persistent (median 3.8 cm/s; rotary autocorrelation \|ρ\| = 0.5–0.8 to 25 h), so the
  quiet upper ocean is real *as a segment mean* — but storm-forced bursts (drift p90 =
  0.24 m/s) are under-represented by construction. Time-resolved estimation is a
  research extension (cf. Stevens-Haas et al. 2022).
* **Error budget** (per-yo, after this pipeline): DAC reference 1–2 cm/s (the floor for
  any absolute statement); declination ≲1 cm/s with IGRF; sound speed ≲1 mm/s after
  correction; compass deviation — check with [`compass_field_check`](@ref) (2° of
  undetected deviation ≈ 1 cm/s); inverse per-bin consistency ≈2 cm/s; shear-method
  integration noise 0.1–0.3 m/s per yo (which is why it is the second opinion).
* **Near-surface bins** carry ~2× the noise (surface waves) and ~3× fewer samples
  (visited only at dive starts and climb ends); treat the top two bins as indicative.
* **Deep bins below the glider's turnaround** are seen only through the 30-m window
  below the deepest pings — `nobs` in the output tables tells you how much data
  supports each bin. Filter on it.
* **Relation to the Python packages**: this implementation reproduces `gliderad2cp`'s
  transform machine-exactly and improves on its small-angle cell-depth approximation
  (exact rotated-beam geometry); it fixes the documented Slocum-AD2CP v2.0.0 transform
  and bookkeeping defects while offering parity modes for regression. The evidence
  trail lives in `docs/research/`.

## References

Visbeck (2002) *JTECH* 19, 794–807 · Fischer & Visbeck (1993) *JTECH* 10, 764–773 ·
Todd et al. (2017) *JTECH* 34, 309–333 · Gradone et al. (2023) *JGR Oceans* 128,
e2022JC019608 · Thurnherr et al. (2015) IEEE/OES CWTM · Queste et al., *gliderad2cp*
(JOSS, in review) · Rudnick, Sherman & Wu (2018) *JTECH* 35, 1665–1673 · Merckelbach
et al. (2010, 2019) *JTECH* · Shcherbina, D'Asaro & Nylund (2018) *JTECH* 35, 411–427 ·
von Appen (2015) *JTECH* 32 · Stevens-Haas et al. (2022) arXiv:2110.10199 · Nortek
N3015-007 Integrators Guide & N3015-011 Principles of Operation · Teledyne RDI, ADCP
Principles of Operation. Full annotated bibliography:
[`docs/research/literature.md`](https://github.com/truedichotomy/gliderADCPjulia/blob/main/docs/research/literature.md).
