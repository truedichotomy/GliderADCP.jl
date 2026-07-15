# Data QA/QC guide

Every finding below was established on real missions (sea064 M37 Jan Mayen 2022,
M38 Lofoten 2022–23, M48 Jan Mayen 2023, M59 NESMA subtropical NW Atlantic 2024)
with the evidence
trail in `docs/research/m38_validation.md`. It is organized as *what we found* →
*what to check on your mission*. The one-line summary: **most glider-ADCP data
problems are silent** — false bottom locks, range-dependent bias, sentinel fill
values, unsynced files — and every one of them was caught by a closure test, not
by inspection of the velocities themselves.

## 1. Beam-sample screening (`qc!`)

The default screens (correlation ≥ 50 %, amplitude window, SNR floor, ambiguity,
surface mask, first cell, error flags — see the tutorial for the full table)
reject **46–53 % of beam samples** on the four validated missions. That number
is normal, dominated by the SNR floor beyond the useful range plus the surface
mask, and it is not hiding signal: loosening the surface mask and keeping the
first cell was tested on M38 and does not change the near-surface answer.

**Check:** `qc!` returns per-screen rejection fractions — log them per mission.
A rejection rate far *below* ~50 % usually means the far cells' noise is being
accepted, not that the data are better.

## 2. Effective range is much shorter than configured range (`cell_quality`)

In clear basin water the per-cell correlation collapses from ~97 % (cell 2) to
~2 % (cell 15): the **effective range is ~15–17 m of the configured 30 m**. The
far cells that survive QC at the default correlation threshold carry a
disproportionate share of the range-dependent shear bias — raising the
correlation threshold from 50 to 80 % cut the measured bias slope by ~22 % on
M38. For conservative shear-method work use `corr ≥ 80`.

**Check:** run [`cell_quality`](@ref) once per mission; if the correlation
profile collapses early, do not expect the outer cells to contribute.

## 2b. The first cell: keep it when blanking ≥ 0.5 m

The fleet flies 0.7 m blanking precisely so cell 1 clears transducer ringing —
and it does, on all four missions: full correlation (96–97 %), amplitude on the
physical decay curve (no ringing spike), no velocity bias against cell 2
(≤ 2.5 mm/s, the ordinary range-dependent pattern), just ~1.5× per-sample noise.
Keeping it adds 12 % of samples (and 1/6 of all telemetered data), improves the
shear-vs-inverse health metric on every mission, and leaves DAC closure and
surface-drift agreement unchanged — so `first_cells = 0` (keep) is now the
default. **Small-blanking deployments (Nortek default ~0.1 m) must set
`first_cells = 1`**; `qc!` warns when blanking < 0.5 m and cell 1 is kept.

**Check:** on a new configuration, compare the mean and std of
`vel[1,:,:] − vel[2,:,:]` against `vel[2,:,:] − vel[3,:,:]` on deep pings — a
clean cell 1 shows the same pair statistics; ringing shows up as an anomalous
offset and inflated variance in the first pair.

## 3. Bottom track: assume false until proven genuine (`bt_valid`)

The single most consequential defect found in this project: on M38, **99.7 % of
15,432 bottom-track locks were false** — a near-field water-borne target 0.6–2.8 m
below the transducer (wake/scattering layer), locked whenever the real seafloor
was out of range. Such targets move *with the water*; feeding them to the
inverse as over-ground anchors contradicts the DAC and injected a spurious
300-m shear layer and 1.6 m/s outliers that looked like ocean structure.

`bt_valid`/`bt_velocity` therefore screen by default with `min_range = 5 m` plus
an impossible-bathymetry test (implied bottom depth vs the platform's own
deepest nearby record). Outcomes across missions — both failure modes exercised:

- **M38** (deep basin): all 9,383 apparent fixes rejected — correct, there was
  never genuine bottom track.
- **M37** (ridge slopes): 16 all-beam locks pass — genuine: glider at 384–782 dbar
  with the seafloor 6–27 m below (implied water depth 390–809 m, consistent with
  local bathymetry).
- **M48** (bank crossings): 148 screened fixes (18 all-beam) pass — genuine:
  implied water depth 162–589 m, independently corroborated by the bathymetry
  arches visible in the mission's own velocity sections.
- **M59** (4–5 km water): all rejected — correct.

**Check:** `nrow(bt_velocity(adcp))` before trusting BT, and for surviving locks
verify the *implied water depth* (glider depth + BT range) against known
bathymetry. A "validation" of BT velocities against water-track velocities is
**not** evidence they are ground-referenced — false locks pass that test
perfectly (they are water-referenced on both sides).

## 3b. The onboard dead-reckoning flight model runs 5–15 % fast — water-track the DAC

The other navigation input the products rest on, and the second one that cannot
be trusted unscreened. The nav-only DAC is `(fix − DR_end)/T`, so it inherits
ALSEAMAR's onboard flight model wholesale. Measured against the ADCP's directly
observed through-water flow (near-cell mean relative velocity, no flight model
involved), the onboard model ran **×1.05–×1.15 fast on all four validated
missions** (×1.13 typical; it behaves like a steady glide-polar with roughly
half the real angle of attack). The result: the onboard DR endpoint overshoots
forward and the nav-only DAC carries a **2–4 cm/s bias against the direction of
travel** — several times the textbook 1–2 cm/s DAC accuracy, sign-flipping with
track direction (zigzag artifacts between opposing transects), and it propagates
1:1 into the inverse's absolute reference. On M48 the onboard DR was noisy
rather than merely biased (per-yo |Δ| up to 12 cm/s, larger than that mission's
median DAC).

`compute_dac(nav, pings)` therefore replaces the onboard displacement with the
integral of [`throughwater_velocity`](@ref) over the same fix-to-fix window and
is the production default. The full ladder, per yo, flagged in the `method`
column: **ADCP water track** (direct measurement) → **flight model**
(`fallback = flight_model(nav)`, or `compute_dac(nav, flight_model(nav))` on
ADCP-less deployments — the package's own steady polar dead-reckons within
1.4/1.4/3.1/1.4 cm/s median of the ADCP water track on M37/M38/M48/M59, with
the systematic along-track bias gone: +0.1/+0.0/+0.5/+0.8 cm/s) → **onboard
estimate** (last resort). The flight-model rung assumes zero vertical water
velocity and travel along heading, and its accuracy tracks the polar's
provenance — the figures above are for the same-glider pooled calibration
(refit per mission with `fit_flightparams(measure_aoa(pings, nav)...)` when an
ADCP is aboard).
Verified against GPS surface drift on M37/M38/M48 (per-yo aggregated, the
onboard-vs-drift disagreement is along-track — the DAC-error signature — and the
water-track form removes its predicted share). On M59 surface drift cannot
arbitrate: a mission-wide ~10 cm/s eastward drift offset (windage/Stokes on the
surfaced hull, plus Gulf Stream surface shear) dwarfs the 4 cm/s question.
Residual of the water-track form: mean shear across the 4–16 m cell offset,
≲ 1 cm/s (window-insensitivity checked at 0.1 cm/s).

**Check:** `compute_dac(nav, pings)` logs its water-track/fallback split; keep
the median |water-track − onboard| with the mission record (expect ~4 cm/s,
mostly along-track). A much larger value flags DR, compass, or ADCP problems —
and a per-mission onboard-speed ratio far from ×1.1 means the vendor changed
the onboard model (worth knowing either way). The per-mission diagnostic
figures (`M*_dac_methods_{delayed,telemetered}.png` from
`examples/dac_methods.jl`) make the ladder visible: U/V sections under all
three DACs plus differences vs the ADCP water track — DAC effects are
depth-uniform stripes, pale for the flight model, track-correlated at
2–5 cm/s for the onboard DR (on M59 a two-week +4–5 cm/s block while the
glider pointed upstream into the Gulf Stream).

## 4. Range-dependent shear bias: measure per mission (`calibrate_shear_bias!`)

All missions carry a range-dependent along-track bias in the beam samples, but
its magnitude is **configuration/mission-dependent, not an instrument
constant**: −4.7/−4.3×10⁻⁴ s⁻¹ on the 2022 missions, −3.1×10⁻⁴ in 2023, and
−5×10⁻⁵ (nearly an order of magnitude smaller) on the same instrument in 2024
(measured with cell 1 included). Left uncorrected it tilts the
shear-method profiles end to end; the inverse partially averages it away.
`calibrate_shear_bias!` measures it with a pairwise-difference estimator
(per-offset means under-correct when depth coverage is partial) and removes it
ping-mean-invariantly to machine zero.

**Check:** always calibrate; log the fitted slope with the mission record. A
slope drifting between missions is diagnostic of configuration changes.

## 5. Compass and attitude (`compass_field_check`)

Undetected compass deviation rotates velocity: 2° ≈ 1 cm/s at typical speeds.
[`compass_field_check`](@ref) verifies the measured field magnitude is constant
across headings (M38: 2.5 % variation — clean; the eddy-yo outliers investigated
there were *not* compass artifacts). Declination comes per ping from IGRF via
navigation positions; queries outside navigation coverage are constant-
extrapolated **with a warning**, never silently dropped.

## 6. Missing data, gaps and coverage (`coverage`, `data_gaps`, `missing_segments`)

Real missions are messy, and the stack's contract is *degrade loudly*:

- Corrupt/unreadable segment files are skipped with a per-file warning (error
  only if nothing parses). Truncated binary downloads report their unparsed tail.
- `missing_segments` lists absent numbers in a transfer sequence; the loaders
  warn automatically.
- `coverage(adcp | nav | pings)` reports spans, gap tables, finite-data
  fractions, GPS-fix and BT counts. M38's duty cycle shows up instantly as 72
  recording gaps totalling 104 days; M59's as one 23-day gap.
- Solvers log solved-of-total segment summaries instead of silently returning
  fewer yos.
- **Cloud-storage placeholder gotcha (seen in the wild):** an unsynced Dropbox
  file reads as zero bytes. A source whose files exist but parse to zero rows
  warns `no rows parsed from N file(s)` — if you see it, check your sync.

## 7. Loader-level QA: sentinels, coordinates, timestamps, multi-route dedup

- **±9999 sentinels** (instrument-off fills in payload logs and GLIMPSE exports)
  parse to missing by default. Before this, `AD2CP_HEADING = 9999` rows entered
  as fake data. Disable only deliberately (`sentinels = nothing`).
- NMEA DDMM.mmm coordinates are converted to decimal degrees in **all** loaders,
  including payload `NAV_LATITUDE`/`NAV_LONGITUDE` (an early-version gap).
- **The gli `Heading` column is TRUE heading** — the vehicle applies its
  configured declination onboard (verified against the AD2CP's magnetic
  compass: the difference reproduces IGRF declination, −14.8° vs −15.3° on
  M59, +1.1° vs +1.2° on M38). Never add declination to nav heading — that
  double-corrects. The AD2CP's own heading **is** magnetic; `process_pings`
  adds declination there, correctly. The two conventions coexist in one
  dataset.
- Epoch-1970 bench rows (clock not yet set) are dropped by default.
- **Multiple download routes** (glider computer + GLIMPSE server) merge with
  exact-timestamp dedup, highest resolution winning; GLIMPSE-only derived
  columns attach to full-resolution rows. GLIMPSE column sets vary by server
  version (M38's 2022 export has `AD2CP_*_c`; NESMA's 2024 export adds
  `LEGATO_SOUND_VELOCITY` etc.) — absent columns degrade to missing **with a
  warning**, never a KeyError.

## 8. Real-time data quality — the three-tier route taxonomy

A glider AD2CP dataset exists in **three tiers**; name them precisely:

| | delayed-mode (`read_ad2cp`) | realtime-onboard (`load_pnor`) | realtime-telemetered (`load_pld_adcp`) |
|---|---|---|---|
| where/when available | post-recovery | on the payload computer during flight | **shore, mid-mission** (Iridium `pld1.sub`) |
| ensembles | every one | every one (~2 s) | one every ~30 s (single subsampled ensemble, not an average) |
| cells | all (15) | all (15) | first 6 |
| amp / corr / BT | yes / yes / yes | yes / yes / no | none |
| quantization | none | 0.01 m/s | 0.01 m/s |
| inverse vs delayed | — (reference) | 3.2–5.1 mm/s rms (four missions) | 28–45 mm/s rms, \|bias\| ≤ 0.8 mm/s (four missions) |

**Shore-side realtime calculations should be built on the telemetered route**
— it is the only realtime tier that exists ashore. The `$PNOR` stream is
useful in real time only to an onboard consumer (e.g. a backseat driver); its
3–5 mm/s result bounds what such a consumer could compute. Both realtime
routes lack the accelerometer — pass `look=` explicitly. The telemetered
product solves the (nearly) identical yo set and lands at the
method-uncertainty floor (~3.2 cm/s); its one casualty is w (r ≈ 0.66–0.84,
rms 5.5–9.5 mm/s — the 30-s subsampling aliases the small, fast vertical
signal). The per-mission diagnostic (`M*_telemetered_w_sections.png`, written
by `examples/realtime_telemetered.jl`) shows what that number means: **large
coherent vertical-velocity events survive telemetry recognizably; the fine
internal-wave banding washes into speckle** — so telemetered w is usable for
spotting significant events mid-mission, not for wave-scale statistics. For
the latter, compute w onboard from `$PNOR`, or wait for delayed data.

**Check:** inspect the w diagnostic per mission rather than assuming the
summary r — event survival is what matters operationally.
For reference, ALSEAMAR's GLIMPSE product — computed **server-side** from the same raw
telemetered data and written back into the server CSV exports as the
`AD2CP_*_c` columns — sits ~3–4× further from the delayed truth on every mission (rms
100–129 mm/s, r = 0.56–0.89, mission-dependent biases to ~19 mm/s, striping
and spurious deep values; biases are quoted against the water-track-referenced
delayed product — against the old onboard-referenced product they read up to
38 mm/s, the difference being the onboard-DR error itself). `AD2CP_TIME` is the instrument clock (MMDDYY) — immune
to the payload-clock bench rows. Check stream coverage against the binary with
`coverage`: M38's payload stopped writing the `$PNOR` stream mid-mission, while
M37's stream held 15 ensembles the instrument card did not retain.

## 9. Pipeline health metrics — compute these on every mission

| check | validated values | what a failure means |
|---|---|---|
| dive vs climb consistency | r = 0.98, med \|Δ\| = 2 cm/s (M38) | transform/sign/geometry errors |
| DAC closure (per yo) | median 1–2 mm/s (all four missions) | referencing errors |
| water-track vs onboard DAC | med \|Δ\| ≈ 4 cm/s, along-track (all four) | ≫ 4 cm/s flags DR, compass, or ADCP trouble; ≈ 0 means the vendor fixed their flight model |
| shear vs inverse agreement | r = 0.90–0.98, rms 3–6 cm/s | contamination anywhere in the chain — this is the check that exposed the false-BT defect |
| surface drift vs shallowest bins | med \|Δ\| = 4 cm/s (M38) | near-surface problems |
| BT plausibility (if any locks survive) | implied depth vs bathymetry | false locks |

The shear-vs-inverse comparison deserves emphasis: the two estimators fail
differently (localization vs integration), which is precisely what makes their
agreement a meaningful end-to-end health metric — and why the shear method is
kept as the standard second opinion rather than retired.

## 10. Interpretation caveats that look like QC problems (but aren't)

- **Near-surface bins** carry ~2× noise (surface waves) and ~3× fewer samples;
  treat the top two bins as indicative. The "suspiciously quiet" M38 surface was
  validated as real (drift median 3.8 cm/s, persistent, not inertial).
- **Everything is a yo-segment mean** (2–6 h): inertial motions vector-average
  toward zero by construction. Weak segment means do not imply weak instantaneous
  currents (M38 drift p90 was 0.24 m/s).
- **Per-yo bin-scale shear is internal-wave dominated** (dive-vs-climb shear
  reproducibility r ≈ 0.08 on M38): low bin-to-bin shear correlation between
  casts is physics, not a defect. Compare shear statistics or averaged profiles.
- **`nobs` is part of the product**: deep bins below the turnaround and sparse
  duty-cycled segments are supported by few samples — filter on it.
