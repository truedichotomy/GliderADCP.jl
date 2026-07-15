# M38 validation notes — GliderADCP.jl vs the prior Python processing

> 2026-07-07. Full-mission run of `examples/currents.jl m38` (sound-speed correction from
> payload CTD, QC, IGRF declination, DAC + bottom-track inverse and shear solutions,
> 127 solved yos — matching the 126 yos in the prior processing; the ADCP was
> duty-cycled Nov 3–27 within the Nov–Mar mission).

## Headline internal-quality metrics (GliderADCP.jl inverse)

| Check | Result |
|---|---|
| Dive vs climb consistency (independent half-yo inversions, same DAC) | r_u = 0.983, r_v = 0.986, med \|Δ\| ≈ 2 cm/s (n = 1514 bins) |
| DAC closure | median 0.005 m/s over 127 yos |
| Glider velocity vs unseen bottom track (DAC-only inverse) | r_v = 0.97, med \|Δ\| ≈ 7 cm/s (n = 807) |
| BT-anchored absolute deep water velocity vs inverse bins | med \|Δ\| = 1.6–1.7 cm/s (n = 136k) |
| Shallow bins (z < 30 m) vs surface GPS drift | med \|Δ\| ≈ 4 cm/s (n = 126 yos) |

## The vertical-structure question (and its resolution)

Our sections differ visibly from the prior processing's figures below ~200 m:
early mission (Nov 4–9) we show **subsurface-intensified** (±0.4–0.5 m/s at 250–500 m)
slope-eddy structure; late mission (deep basin) we show near-barotropic per-yo columns,
while the prior figures show smooth **surface-intensified** profiles throughout.

A pure baroclinic sign flip preserves depth means, dive/climb consistency, and DAC
closure — so those checks cannot arbitrate. Three independent arbiters were run:

1. **End-to-end synthetic with depth-varying flow** through the full beam forward model
   (now a permanent regression test): relative velocities exact at every offset;
   inverse recovers du/dz upright and at the right magnitude. Our chain cannot flip
   structure.
2. **Raw transform-level tilt check** (no solver, no reference): binned mean relative
   velocity E_rel(z) within single yos. Since glider velocity is ~constant on average
   over a yo, the *shape* of E_rel(z) ≈ u_ocean(z) + const:
   - Nov 5 / Nov 7 yos: raw tilt (surface − 300 m) = **−0.16 / −0.18 m/s** —
     subsurface intensification is IN THE RAW DATA (real slope-current/eddy signal);
     the prior processing agrees in sign here (−0.06 / −0.18).
   - Nov 22 yo: raw E_rel is **flat** (−0.32 ± 0.01 from 12–760 m) — barotropic.
     Our profile is correspondingly flat; the prior profile imposes a +0.13 m/s
     surface intensification that is **not present in the raw data**.
3. **BT-anchored absolute velocities** (u_rel + u_glider_over_ground from bottom track,
   no DAC/inverse involved) agree with our deep inverse bins to 1.6 cm/s median.

**Conclusion:** GliderADCP.jl sections are faithful to the raw measurements. The prior
Python profiles agree where the signal is strong (early-mission shelf/slope yos — the
matched-yo correlation there is r ≈ 0.8) and diverge where profiles are weak, where
that pipeline's over-smoothing dominates: it ran wSmoothness = 1 at dz = 1 m —
a ~100× stiffer curvature penalty than at the documented dz = 10 m — plus the
documented v2.0.0 transform defects (dive-cast matrix misalignment, halved X/Z,
possible rotation transpose). Mission-wide per-yo correlation against that reference
is therefore low (median r_u ≈ 0.1) *by expectation*, and is not a defect indicator
for this package.

## Notes / future improvements

- Shear-vs-inverse intercomparison: r_u = 0.58, rms ≈ 0.2 m/s pooled. The shear path
  currently uses simple depth-mean DAC referencing; per-yo reference offsets and
  integration drift dominate the discrepancy. Planned: time-in-bin-weighted referencing
  (gliderad2cp semantics) and per-cast integration.
- QC rejects 52 % of beam samples on the full mission with default thresholds (SNR
  floor + amplitude + correlation + velocity cap + surface mask); revisit per-screen
  contributions when tuning for deep, quiet water.
- Surface-drift comparison is noisy (windage/Stokes on a surfaced glider); it bounds
  gross errors but should not be over-interpreted below ~5 cm/s.

## Why the shear and inverse methods disagree (2026-07-07 diagnosis)

Pooled shear-vs-inverse agreement on M38 is r_u ≈ 0.58, rms ≈ 0.22 m/s, worst near the
surface. Quantitative diagnosis (script results, full mission):

1. **A measured range-dependent velocity bias.** The mission-mean *within-ping* relative
   velocity, rotated into the glider track frame, decays with range: along-track slope
   **−3.3×10⁻⁴ s⁻¹ (dives) and −3.1×10⁻⁴ s⁻¹ (climbs)** — same sign and size for both
   cast directions — while the cross-track component is 10–100× smaller (−2×10⁻⁵ / +3×10⁻⁶).
   Per cell the bias is tiny (±3.5 mm/s across the 30-m window); it is invisible to any
   per-sample QC. This is the classic glider-ADCP "shear bias" (Todd et al. 2017;
   gliderad2cp's `process_bias`, whose velocity-dependent variant exists for exactly this
   speed-correlated, track-aligned signature).
2. **The shear method integrates it; the inverse does not.** Integrating −3.2×10⁻⁴ s⁻¹
   over a 500-m yo predicts a −0.16 m/s top-to-bottom tilt; the observed per-yo tilt of
   (u_shear − u_inverse) is **median −4.6×10⁻⁴ s⁻¹ ≈ −0.23 m/s per 500 m** — same sign,
   same order. DAC referencing then pivots the tilted profile about the yo's
   (time-weighted) mean depth, so the shear solution comes out **+0.13 m/s too fast near
   the surface and −0.10 m/s too slow at depth** (median u_inv − u_sh: −0.149 at 0–50 m,
   +0.105 at 600–1000 m). In the inverse, each sample constrains only its own depth bin
   against the per-ping glider unknown — the bias adds ~mm/s per bin instead of
   integrating.
3. **Independent surface arbiter.** Against surface GPS drift (z < 30 m bins):
   inverse med|Δ| = 0.042 m/s, bias −0.025; **shear med|Δ| = 0.151 m/s, bias +0.111** —
   the shear top is confirmed biased, the inverse is not.
4. **Integrated noise.** Beyond the tilt, integrating per-bin shear noise random-walks:
   rms(u_inv − u_sh) is 0.18–0.25 m/s at all depths, vs ~2 cm/s per-bin consistency for
   the inverse (dive/climb check). The upper column additionally has ~3× fewer samples
   per bin (median 46–67 vs ~180 at mid-depth) — sampled only during dive starts and
   climb ends.
5. **Inverse outliers noted:** 289 bins reach |u| up to 1.6 m/s at 345–605 m, all in
   yos 8–18 (the slope-eddy transit, DAC up to 0.51 m/s). Cast-consistent (dive/climb
   r = 0.98) so not a transform artifact, but flagged for the Phase-7 compass-deviation
   check.

**Conclusion:** the disagreement is a property of the *shear method* inheriting an
integrable instrument bias plus integrated noise — not of the common trunk (shared by
both methods) nor the inverse. The inverse (+ bottom track) sections are the reference
product. **Phase-7 item 1:** shear-bias correction — subtract the mission-calibrated
along-track bias profile b(r) before differencing (real ocean shear averages out of the
track-frame mission mean under varied headings), and/or the gliderad2cp
displacement-regression form; acceptance = shear-vs-inverse tilt collapse and drift-test
parity.

## Shear-bias correction: implementation and what it does (and doesn't) fix (2026-07-07)

Implemented as `shear_bias` / `apply_shear_bias!` / `calibrate_shear_bias!`:

- **Estimator**: mission-mean of the *adjacent-pair differences* of track-frame relative
  velocities — the exact sample population the shear estimator consumes (with partial
  coverage, mean-of-differences ≠ difference-of-means, so the earlier per-offset mean
  profile under-corrected). Pair differences are integrated to a mean-removed bias
  profile B(offset). An optional `velocity_scaled=true` form fits bias ∝ ping
  through-water speed (it did not outperform the plain form on M38).
- **Application** re-demeans the subtracted profile over each ping's finite cells, so
  ping-mean velocities are unchanged **by construction** — the inverse's glider-velocity
  and DAC content is untouched (verified: median profile change 3 mm/s).
- **M38 calibration**: slope −3.9…−4.1×10⁻⁴ s⁻¹ along-track, cross-track ~10⁻⁶,
  heading concentration R = 0.23 (safely diverse). One pass removes it to machine zero,
  and the **residual pairwise bias is < 1×10⁻⁴ s⁻¹ in every cell-depth band** (was
  ≈ −4×10⁻⁴ uniformly; also verified depth-stratified: −3.0…−3.8×10⁻⁴ before).

Effect on M38 solutions: shear-vs-inverse r_u 0.582 → 0.648, rms 0.217 → 0.197 m/s;
per-yo tilt median −4.6 → −3.3×10⁻⁴ s⁻¹; surface-drift bias of the shear top
+0.111 → +0.087 m/s. Dive/climb and DAC-closure metrics of the inverse unchanged.

**What the remaining disagreement is (investigated to closure):** after correction the
shear *samples* are provably unbiased (< 10⁻⁴ s⁻¹ at every depth), the bin statistic is
not the cause (mean and median give identical tilt), and the discrepancy vs the inverse
concentrates at 200–700 m as per-yo *integration drifts* of ±0.1–0.3 m/s (e.g. a yo
whose raw relative velocity is constant to ±0.01 from 215–505 m shows a −0.2 m/s shear
drift while the inverse stays flat). The residual median tilt (−3.3×10⁻⁴ s⁻¹ ≈ 3× the
median's standard error given the ±1.3×10⁻³ IQR) is the noise-skew of those random
walks, not a removable sample bias. This is the intrinsic error-propagation difference
between the methods: **integration accumulates what the inverse localizes.** The shear
product is retained as a corrected, second-opinion diagnostic; the inverse (+ bottom
track) remains the reference. Possible future incremental gains: per-cast integration
with cast averaging (halves drift variance), overlap-constrained integration.

## The "suspiciously quiet surface" question (2026-07-07)

Concern: the inverse's 0–30 m velocities look low. Findings:

1. **Not a solver artifact.** Top-bin values and the drift comparison are insensitive to
   the smoothness weight (0.2/1/5 identical to 3 decimals), bin size (dz 5 vs 10 m), and
   looser surface QC (2-m mask, keep first cell: med|u_top| 0.026 → 0.035, drift
   agreement slightly worse). Near-surface bins do carry ~2× the within-ping noise
   (waves; 0.022 vs 0.009 m/s std) and ~3× fewer samples.
2. **Independent ground truth agrees.** Surface GPS drift vs 0–30 m inverse bins:
   med|Δ| = 0.040 m/s, bias −0.023, across 126 yos. The drift itself is weak: median
   0.038 m/s (p90 0.24).
3. **The surface flow was persistent, not inertial.** Rotary autocorrelation of the
   545 drift vectors: |ρ| = 0.5–0.8 out to 25 h with a weak clockwise drift (−20° at
   the 12.8-h inertial period, not the −360° a rotating inertial current would give).
   November basin interior with a deep mixed layer: the quiet segment-mean surface is
   real at surfacing times.
4. **Honest caveat (the valid core of the intuition):** all our products — inverse,
   shear, DAC, and the drift matching — are yo-segment means. Currents with timescales
   shorter than a yo (storm-forced inertial bursts; drift p90 shows 0.2+ m/s events)
   vector-average toward zero and are under-represented. Resolving them needs a
   time-resolved estimator (per-cast solves are noise-limited in the top bins with only
   ~5–10 pings; a rotary/Kalman formulation à la Stevens-Haas et al. 2022 is the proper
   future tool).

## Task 1 (2026-07-08): shear-content comparison between the methods

New products: [`solve_shear_profile`] (the shear method's pre-integration binned shear)
and [`inverse_shear`] (centered differences of any solver profile — the shear implied by
the inverse). Comparing them on M38 (post-QC, post-calibration, dz = 10 m):

| comparison | r_u | notes |
|---|---|---|
| per-yo bins, direct vs inverse-implied | 0.20 (0.63 at z < 100 m) | rms 1.9×10⁻³ s⁻¹ ≈ the internal-wave shear level |
| same, scale-matched (both profiles 20-m centered-differenced) | 0.22 | scale is not the cause |
| **dive vs climb of the DIRECT product (same method, same yo)** | **0.08** | the reproducibility ceiling |
| mission-mean profile, v | 0.39 (rms 1.6×10⁻⁴; stds 1.4/1.5×10⁻⁴) | agrees |
| mission-mean profile, u | −0.04 (rms 6.7×10⁻⁴) | see below |

**Interpretation.** Per-yo, bin-scale (10–20 m) shear in this environment is dominated
by internal-wave finestructure that decorrelates in the hours between a dive and its
climb — the direct product cannot even reproduce *itself* cast-to-cast (r = 0.08), so
no two estimators can be expected to agree per bin. The expectation "the shear content
should be comparable" holds exactly where signal exceeds the wave noise: the upper
100 m (r = 0.63), the integrated velocity profiles, and the mission-mean v profile.

**The mission-mean u asymmetry** (inverse std 5.3×10⁻⁴ vs direct 1.6×10⁻⁴): the
inverse retains the real sub-inertial shear of the early-mission eddy transit
(magnitudes match the raw-data tilt checks of −5…−6×10⁻⁴ s⁻¹), while the direct
product does not. Working hypothesis: the shear-bias calibration operates entirely
within the ping window, so *track-aligned real mean shear* from the predominantly
eastbound eddy transit leaks into b(r) and is subtracted from the direct product — the
documented heading-diversity caveat, here quantified at O(1–5)×10⁻⁴ s⁻¹ for this
mission's u component — whereas the inverse's shear comes from cross-ping bin structure
and is immune. Task 3 (the 300-m shear feature) will adjudicate the inverse's deep
structure independently via hydrography/thermal wind.

**Guidance** (also in the tutorial): use `inverse_shear` for deterministic/sub-inertial
shear; use `solve_shear_profile` for internal-wave/finestructure statistics (shear
variance) and QC; on low-heading-diversity legs expect the bias calibration to absorb
some real track-aligned mean shear from the *direct* product only.

## Task 3 (2026-07-08): the "strong shear at 300 m" — solved, and a defect found

**Verdict: the feature was an artifact of false bottom-track locks, now screened out by
default.** The investigation, in order:

1. **Localization.** The strongest inverse shear concentrated at 275–415 m, mostly in
   the slope-transit yos (peaks to 1.2×10⁻² s⁻¹ at ~295–395 m), with a milder
   mission-wide elevation in the same band.
2. **Sensitivity.** Smoothness weight: no effect. Referencing: decisive —
   DAC-only median |shear| at 200–450 m in the transit yos is 4.4×10⁻⁴ s⁻¹;
   adding bottom track raises it 10× to 4.2×10⁻³.
3. **Hydrography veto.** The deep pycnocline sits at 160–260 m with
   dσ₀/dz ≤ 0.16 kg/m³ per 100 m (N ≈ 4×10⁻³ s⁻¹): a sustained 10⁻² s⁻¹ shear there
   would put Ri < ¼. Yo-pair thermal wind at 250–400 m supports only
   ~2–5×10⁻⁵ s⁻¹. The measured relative velocities (raw tilt checks) show the smooth
   DAC-only structure, not the sharp feature.
4. **Root cause.** The BT-anchored "absolute" deep velocities were ≈ 0.00 m/s in every
   transit yo while DAC and DAC-only inverse agreed at +0.2…+0.5 — the signature of a
   target moving **with the water**. Classification of all 15,432 three-beam BT locks
   (impossible-bathymetry test: platform later dove deeper than the implied bottom):
   **99.7 % false**, locking on a persistent target 0.6–2.8 m below the transducer
   (median 1.7 m) — near-field/wake, not seafloor (basin depth 1–3 km; the seafloor was
   never within the 30-m BT range on this mission). Anchoring pings to a water-frame
   target contradicts the earth-frame DAC; least squares dumps the contradiction into
   the ocean profile as spurious shear near the top of the anchored depth range — and
   pushes unanchored mid-depth bins to |u| up to 1.6 m/s (the previously flagged
   yos-8–18 outliers, now fully explained; their dive/climb consistency was the blind
   spot — both casts shared the same false anchors).
5. **Fix.** `bt_valid` now defaults to `min_range = 5 m` (genuine seafloor approaches
   have O(10 m) ranges; the false cluster is 0.6–2.8 m) plus an impossible-bathymetry
   screen (implied bottom vs deepest platform record within ±2 h). On M38: 0 of 9,383
   fixes survive (correct — there were no real locks), the transit shear collapses to
   the DAC-only value, the |u| > 0.8 outliers drop from 289 bins to 0, and the
   basin-wide 250–350 m elevation disappears (4.9×10⁻⁴, smoothly decaying with depth —
   **the entire "300-m shear" was this artifact**).

**Amendments to earlier claims (honesty pass):**
- "Independent validation: DAC-only inverse u_g vs unseen bottom track r_v = 0.97,
  med|Δ| ≈ 7 cm/s" — the BT reference was the water-frame target, so this validated the
  **through-water** velocity content of the chain (still a real geometry check; the
  7 cm/s offset ≈ the water speed), *not* an absolute over-ground reference.
- "BT-anchored absolute deep velocities agree to 1.6 cm/s" — self-consistency of
  water-frame anchoring, not an absolute check. Superseded.
- The BT **sign-convention** check (w vs dP/dt, r = 0.90) survives: vertical motion
  relative to a water-borne target still validates the sign chain.
- Headline metrics that never involved BT stand unchanged: dive/climb r = 0.98
  (2 cm/s), DAC closure 5 mm/s, surface drift 4 cm/s, raw-data tilt checks, all
  synthetic and gliderad2cp parity results.
- Task 1's mission-mean u-shear asymmetry: **was BT-injected**, not calibration
  absorption (that hypothesis is retracted) — with hardened screens the DAC-only
  inverse mission-mean u-shear std is 6.8×10⁻⁵ s⁻¹, consistent with the direct
  product.

**Physics answer to the original question:** there is no strong sub-inertial shear
layer at 300 m in this mission. The real deep pycnocline (AW base) sits at 160–260 m
with thermal-wind shear of order 10⁻⁵–10⁻⁴ s⁻¹, and the transit-eddy velocity
structure is smooth across 200–450 m (≈4×10⁻⁴ s⁻¹), exactly as the DAC-only inverse
and the raw relative velocities show.

## Task 4 (2026-07-08): section figures, w products — and a major revision

`examples/currents.jl` now produces four-panel U/V sections (inverse and shear,
both components), a two-panel **w** figure via the new [`solve_w`] product
(`:direct` = binned `U_rel + dP/dt`; `:inverse` = ocean-w bins solved jointly with
per-ping glider w anchored to the pressure-derived vertical velocity — the vertical
analog of a bottom-track constraint), and data-driven symmetric color ranges (99th
percentile of the finite values; M38: ±0.5 m/s horizontal, ±0.04 m/s vertical). The two
w estimates corroborate each other closely; near-surface bins remain wave-contaminated
as documented.

**Revision of the method-intercomparison numbers.** With the false bottom-track locks
screened out (Task 3), the M38 full-mission comparison becomes:

    shear vs inverse:  r_u = 0.977, r_v = 0.979, rms ≈ 0.036 m/s   (was 0.65 / 0.20)
    DAC closure:       median 0.001 m/s                            (was 0.005)

Most of the disagreement previously attributed to shear-method integration noise was in
fact the **inverse being distorted by the false BT anchors**. The corrected picture:
the two methods agree to ≈3–4 cm/s rms mission-wide. What survives of the earlier
narrative: the shear method's per-yo random-walk drift and its ≈+0.09 m/s near-surface
bias against GPS drift (both measured BT-free) are real but second-order; the
cast-to-cast shear-reproducibility ceiling (r = 0.08, method-internal) also stands.
Task-1's per-yo bin-shear comparisons used the contaminated inverse and would improve
somewhat if recomputed, but remain bounded by that sampling ceiling.

## Task 5 (2026-07-08): real-time vs delayed-mode products

Full-pipeline comparison of the two AD2CP data routes, everything else held identical
(gli.sub nav, pld1.sub CTD sound speed, QC thresholds, IGRF declination, shear-bias
calibration, DAC): the **real-time** `$PNOR` ASCII telemetry stream (`load_pnor`;
0.01 m/s velocity quantization, 0.1° attitude, no accelerometer, no BT records) vs the
**delayed-mode** full-resolution `.ad2cp` binary. Script:
`examples/realtime_onboard.jl m38`; gated acceptance test in the suite.

Coverage first (the Task-6 machinery reports it directly): the stream carries
123,950 of 124,752 ensembles (99.4%). The payload stopped writing the stream on
2022-11-27 — segment sequence complete, `missing_segments` empty, so this is the
payload configuration, not transfer loss — while the instrument kept recording sparse
internal bursts through 2023-03-01 (750 ensembles, three months, delayed-only).
The real-time record covers the entire main pinging period.

Product agreement over the 127 yos both sides solve (identical yo sets — real-time
loses no segments; common (yo, z) bins, nobs > 10):

    inverse u:  r = 0.9996   rms = 4.6 mm/s   bias  0.0 mm/s
    inverse v:  r = 0.9997   rms = 4.1 mm/s   bias  0.0 mm/s
    shear u:    r = 0.9897   rms = 24.9 mm/s  bias +0.1 mm/s
    shear v:    r = 0.9911   rms = 22.9 mm/s  bias -0.1 mm/s
    w (direct): r = 0.9603   rms = 2.6 mm/s   bias  0.0 mm/s

Depth structure: the inverse difference is flat at 4–6 mm/s from the surface to
1000 m (5.7 mm/s in 600–1000 m, where fewer samples average the quantization noise);
the shear-method difference runs 18–31 mm/s because vertical integration accumulates
the quantized-sample noise that the inverse localizes per bin.

**Conclusions.** (1) A real-time/onboard product built from the telemetry stream is
essentially the delayed product: the inverse-method penalty is ~5 mm/s rms with zero
bias — an order of magnitude below the ~3–4 cm/s method/sampling uncertainty
established in Task 3. The 0.01 m/s per-sample quantization averages down as expected
(~200 samples per bin → few mm/s). (2) The stream's missing pieces (accelerometer →
pass `look` explicitly; BT records → irrelevant on M38, which has no genuine BT;
magnetometer → declination comes from nav anyway) cost nothing here. (3) Prefer the
shear method last in real-time settings: it is the one product measurably degraded
(≈2.5 cm/s rms) by stream quantization.

**Cross-mission confirmation (2026-07-08).** The same comparison run on M37
(Jan Mayen, Oct 2022) and M59 (NESMA subtropical NW Atlantic, Jul–Aug 2024;
`examples/realtime_onboard.jl`):

    M37: inverse r = 0.9987/0.9984 (u/v), rms 3.7/4.4 mm/s;  shear rms 21–22 mm/s;  w rms 0.5 mm/s
    M59: inverse r = 0.9997/0.9996 (u/v), rms 5.1/5.0 mm/s;  shear rms 28–29 mm/s;  w rms 0.3 mm/s

Both streams cover 100% of their pinging windows, both routes solve identical
yo sets (107 and 154), and every bias is zero to 0.1 mm/s. Notably M59 holds
the ~5 mm/s inverse agreement while the glider crosses a >1 m/s Gulf Stream
jet — the real-time penalty does not scale with signal amplitude. The
conclusion generalizes: a telemetry-stream inverse product is the delayed
product to ~5 mm/s rms on every mission tested; only the shear method pays a
measurable (2–3 cm/s) quantization cost.

## Three-mission workflow validation (2026-07-08): M37 and M59

The full stack (SeaExplorerIO multi-route loading → QC → calibration → DAC/BT →
all solvers → figures) run end to end on two further missions
(`examples/currents.jl m37 m59`):

**M37 (Jan Mayen, Oct 2022, ridge slopes).** Native binary only (80,126 ens,
near-continuous). Shear-bias slope −4.29×10⁻⁴ s⁻¹ — matching M38's −4×10⁻⁴.
DAC closure median 2 mm/s; shear-vs-inverse r = 0.91/0.92 (rms 3.4 cm/s).
**The BT screens passed genuine locks for the first time**: 16 all-beam fixes,
glider at 384–782 dbar with seafloor 5.6–26.6 m below (implied water depth
390–809 m, consistent with the ridge bathymetry) — the same defaults that
reject 100 % of M38's false locks. GLIMPSE nav = strict subset of delayed
(38,704 rows, all deduplicated); the $PNOR stream holds 15 ensembles the
instrument card did not retain.

**M59 (NESMA, subtropical NW Atlantic, Jul–Sep 2024).** Both binary and MIDAS
netCDF exist: reader parity re-verified on a third mission (204,248 ensembles,
max |Δvel| = 0.0). The glider crosses a Gulf Stream jet (u > 1 m/s reaching
below 500 m); DAC closure median 1 mm/s over 154 yos; shear-vs-inverse
r = 0.96/0.95 (rms 7 cm/s against a ±0.85 m/s signal). Zero BT fixes survive
screening — correct in 4–5 km of water. **Shear-bias slope −3.07×10⁻⁵ s⁻¹, an
order of magnitude below the 2022 missions on the same instrument** → the bias
is configuration/mission-dependent, not an instrument constant; it must be
measured per mission (which `calibrate_shear_bias!` does).

**GLIMPSE-route findings** (SeaExplorerIO 0.2.x): server exports carry
version-dependent extra columns (M38's 2022 export: `AD2CP_*_c` onboard current
estimates, 41,892 finite values; NESMA's 2024 export adds
`LEGATO_SOUND_VELOCITY` etc.); merged loading attaches them to full-resolution
rows while deduplicating all telemetered duplicates. On every mission checked
the GLIMPSE record was a strict subset of the glider-computer download — the
merge's value is the extra columns plus insurance against incomplete downloads.
Zero-byte unsynced cloud placeholders are reported (`no rows parsed`), not
silently skipped.

## Method verdict (2026-07-08): the inverse is the production method

Recorded as the project's standing conclusion, with its evidence and its limits.

**The case is mechanistic, not just empirical.** The real-time comparison
isolated the structural difference: feed both methods identically quantized
samples and the inverse's error stays a flat 4–5 mm/s floor to 1000 m while the
shear method's grows with depth to 2–3 cm/s — on three missions in three ocean
regimes, independent of signal amplitude (unchanged through M59's >1 m/s jet).
Integration accumulates what the inverse localizes. The same mechanism accounts
for the shear method's per-yo random-walk drift (±0.1–0.3 m/s), its full
inheritance of the range-dependent bias, and its ≈+0.09 m/s near-surface bias
vs GPS drift. The inverse additionally fuses constraints the shear method
cannot (genuine BT on M37, surface drift, pressure-anchored w), closes the DAC
at 1–2 mm/s, and reproduces dive-vs-climb at r = 0.98.

**Why not "undeniably better", stated for the record:**
1. Both methods stand on the same DAC reference — absolute accuracy is still
   dominated by navigation quality, and that layer has only been validated
   internally (closure, drift), never against an independent instrument
   (mooring / ship ADCP). That is the remaining gap.
2. In the complete-sampling, uniform-weighting limit the two estimators
   coincide (Visbeck 2002); the inverse's advantage is robustness to real-world
   imperfection — gaps, uneven bin occupancy, QC holes, quantization — which
   our missions show is always present.
3. The inverse has regularization knobs; misconfigured smoothing produces
   confident-looking wrong answers (the over-smoothed legacy Python profiles).
   The shear method's noise is at least visibly noise.
4. The false-BT episode cuts both ways: for a time the *inverse* was the
   contaminated product, precisely because it is the method that ingests
   constraints. Power and attack surface come together.

**Standing role of the shear method:** not a competitor but the standard second
opinion. The two estimators fail differently, so their agreement
(r = 0.92–0.98 across missions) is the pipeline's most valuable end-to-end
health metric — it is the check that exposed the false-BT contamination.

## M48 (2026-07-09): fourth mission, added via the symmetric registry

sea064 M48 (Jan Mayen, Nov 2023) processed with one registry entry in
`examples/missions.jl` — no code changes. 124,081 ensembles, near-continuous
(5 gaps, 0.7 days), 255 yos, a shallow bank-crossing mission (profiles truncate
at ~200 m over the crests).

- **Binary ↔ MIDAS netCDF parity on a third export: max |Δvel| = 0.0.**
- **Genuine bottom track again** (Jan Mayen ridge): 148 screened fixes
  (18 all-beam), implied water depth 162–589 m — independently corroborated by
  the bathymetry arches visible in the velocity sections.
- **Shear-bias slope −2.44×10⁻⁴ s⁻¹** — a 2023 point between the 2022
  (−4.3×10⁻⁴) and 2024 (−3.1×10⁻⁵) values; the per-mission calibration rule
  stands.
- Health checks: DAC closure median 2 mm/s (elevated to 1–3.5 cm/s in two
  clusters of yos coinciding with the bank crossings — short shallow yos with
  fewer bins constrain closure less); dive-vs-climb med |Δ| 2.2 cm/s
  (r = 0.82/0.86); shear-vs-inverse rms 3.6 cm/s at r = 0.88/0.90 — the same
  rms as the other missions, with lower r because this weak-flow regime has
  less signal variance, not more disagreement.
- **Real-time vs delayed, fourth confirmation**: inverse r = 0.9987/0.9990,
  rms 3.2/3.3 mm/s, zero bias, identical 255 yos (stream covers 98.9%);
  shear method ~2 cm/s; w rms 0.5 mm/s. Cross-mission inverse range is now
  3.2–5.1 mm/s over four missions.

## The second real-time route (2026-07-11): telemetered pld1.sub AD2CP pings

Two "real-time" AD2CP data routes exist on a SeaExplorer, and only one of them
reaches shore mid-mission:

1. **`$PNOR` ASCII stream** (`ad2cp.raw` logs): every ensemble, all 15 cells,
   amplitudes + correlations — but *payload-logged only*, recovered with the
   glider. The earlier "real-time vs delayed" work (Task 5 and the cross-mission
   confirmations) characterized this route; it bounds an *onboard* product.
2. **AD2CP subset inside `pld1.sub`** — what Iridium actually transmits: one
   subsampled ensemble every ~30 s (`AD2CP_TIME`, MMDDYY HH:MM:SS instrument
   stamp; matched rows are the instrument ensemble quantized to 0.01 m/s —
   verified against the binary at max |Δv| = 0.005, while a 30-s-average
   hypothesis fails at 0.15–0.94), beam coordinates, cells 1–6 of 15, attitude
   to 0.1°, pressure — no amplitude, correlation, or bottom track. ALSEAMAR's
   proprietary GLIMPSE processing of these rows appears as `AD2CP_*_c` columns
   in server exports.

`load_pld_adcp` now ingests route 2 (glider-computer segment logs + GLIMPSE
`.all.csv`/per-cycle exports, deduplicated on the instrument timestamp) into the
standard `AD2CPData`; the config (cell size, blanking) comes from the deployment
plan, the onboard sound speed is reconstructed from the configured salinity +
payload CTD temperature so the standard correction applies, the missing QC
fields degrade to no-op screens, and `first_cells` removes the ringing cell —
leaving 5 usable cells (2.7–12.7 m below the transducer, inside the ~15–17 m
effective range, so the discarded far cells are ones QC would have rejected).

**M38 shore-side product vs delayed** (`examples/realtime_telemetered.jl`):
41,354 telemetered pings (vs 124k ensembles) still solve the identical 127 yos;
even the shear-bias calibration works from 5 cells (−3.4×10⁻⁴ vs −3.9×10⁻⁴
full-range). Inverse agreement with the delayed product:

    u: r = 0.977, rms = 36 mm/s, bias −0.1 mm/s
    v: r = 0.978, rms = 35 mm/s, bias  0.0 mm/s
    w: r = 0.67,  rms = 8.5 mm/s        (subsampling aliases the small, fast w signal)

That is: **the true shore-side product lands at the method-uncertainty floor**
(~3.6 cm/s, the shear-vs-inverse spread) — an order of magnitude above the
$PNOR route's 4–5 mm/s (which has 30× the pings and all cells), but fully
usable scientifically, with zero bias and every section feature reproduced.

**ALSEAMAR's own GLIMPSE product**, binned to the same (yo, z) grid and compared
to the delayed inverse: r = 0.80/0.82, rms = 131/107 mm/s, with visible striping
artifacts and spurious deep values in the duty-cycled period. Per-ping (its
native form) it correlates with the delayed inverse at the glider's depth at
r ≈ 0.76–0.79 (QF = 0). The open pipeline on the identical telemetered input is
~3× closer to the delayed truth than the proprietary product. `Utot_c`/`Udir_c`
do not decode as magnitude/direction of `Ueast_c`/`Unorth_c` (units/semantics
unresolved); `QF_c ∈ {0, 3, 4}`.

Caveats for the record: the telemetered w product is the one real casualty of
the 30-s subsampling (use the $PNOR route or delayed data for w); and the
`AD2CP_TIME` clock is the instrument's, which also sidesteps the payload-clock
bench rows (early M38 payload stamps read 2019 while `AD2CP_TIME` is correct).

## First-cell verdict (2026-07-11): the 0.7 m blanking did its job — keep cell 1

The fleet's AD2CP configuration uses 0.7 m blanking (vs the ~0.1 m Nortek default)
by design, to push cell 1 clear of transducer ringing. Verified on all four
missions (all 0.7 m/2.0 m):

- **Instrument level**: cell-1 correlation 96–97% (identical to cells 2–3);
  amplitude sits on the geometric-spreading + absorption decay curve (the 6–7 dB
  drop to cell 2 matches 20·log₁₀(4.7/2.7) + absorption ≈ 5.8 dB) — no ringing
  spike above the curve.
- **Sample level**: adjacent-cell velocity differences show no anomalous (1,2)
  bias (0.1–2.5 mm/s, same fore/aft antisymmetric pattern as the (2,3) and (3,4)
  pairs, i.e. the ordinary range-dependent bias); per-sample noise is ~1.5× cell 2
  (near-field transition) — unbiased-but-noisier, which averages down.
- **Product level (M38 A/B, first_cells 1 → 0)**: +12.5% samples;
  shear-vs-inverse health metric improves r 0.977 → 0.983, rms 37 → 32 mm/s;
  DAC closure and surface-drift agreement unchanged (the earlier "keeping the
  first cell slightly degrades drift agreement" result was the bundled
  surface-mask loosening, not cell 1). Telemetered route (cell 1 = 1/6 of the
  data): bin occupancy +29%, agreement with delayed truth improves on every
  product (u 36 → 33, v 35 → 31 mm/s rms; w r 0.67 → 0.71).

**Default changed**: `QCThresholds.first_cells` is now 0 (keep cell 1), with a
`qc!` warning when the configured blanking is < 0.5 m (small-blanking deployments
should set `first_cells = 1`). Re-measured four-mission health metrics with
cell 1 included: shear-vs-inverse r = 0.90–0.98 at 30–63 mm/s rms (improved on
every mission), DAC closure 1–2 mm/s, dive/climb ≈ 2 cm/s — and the measured
shear-bias slopes shift slightly with the near cell included
(−4.7/−4.3/−3.1/−0.5 ×10⁻⁴ s⁻¹ for 2022/2022/2023/2024), still
mission-dependent, still calibrated per mission.

**Cross-mission confirmation (2026-07-11, all four missions).** The telemetered
shore-side product and the ALSEAMAR comparison, run identically everywhere
(`examples/realtime_telemetered.jl`, cell 1 included):

    mission  tele pings  yos            our tele inverse (vs delayed)     ALSEAMAR (vs delayed)
    M37        26,773    104 of 107     r=0.93/0.94  rms 29/28 mm/s       r=0.65/0.64  rms 111/101 mm/s
    M38        41,354    127 of 127     r=0.98/0.98  rms 33/31 mm/s       r=0.81/0.82  rms 129/107 mm/s
    M48        41,049    255 of 255     r=0.86/0.90  rms 35/33 mm/s       r=0.39/0.45  rms 119/127 mm/s (v bias +33 mm/s)
    M59        68,106    154 of 154     r=0.98/0.97  rms 42/45 mm/s       r=0.90/0.82  rms 114/119 mm/s (u bias +38 mm/s)

The open telemetered product sits at 28–45 mm/s rms with |bias| ≤ 0.8 mm/s on
every mission — at or below the shear-vs-inverse method spread — with r tracking
signal variance (lowest in M48's weak-flow regime at unchanged rms). ALSEAMAR's
product is 101–127 mm/s rms everywhere, with mission-dependent biases up to
38 mm/s. w from the telemetered route: r = 0.66–0.84, rms 5.5–9.5 mm/s — the
subsampling caveat holds on all missions, and the per-mission diagnostic
(`M*_telemetered_w_sections.png`, delayed vs telemetered side by side) shows
its character: **aliasing degrades texture, not events** — the large coherent
up/downwelling columns (±0.02 m/s class, e.g. M38's Nov 5–7 events) come
through recognizably while the fine internal-wave banding washes into speckle.
Operationally: telemetered w flags significant vertical-velocity events
mid-mission; wave-scale statistics need the onboard or delayed tiers. M37's stream supports 104 of the 107
delayed yos (three short segments fall under `min_pings` at 30-s cadence).

**Route taxonomy adopted (2026-07-11).** Standing names for the three data tiers:
**delayed-mode** (`.ad2cp` binary, post-recovery — the reference),
**realtime-onboard** (`$PNOR` stream — payload-logged; in real time useful only to
an onboard consumer such as a backseat driver; its 3–5 mm/s result bounds that use
case), and **realtime-telemetered** (the `pld1.sub` AD2CP subset — the only tier
that exists ashore mid-mission; **shore-side realtime calculations are built on
this route**). Confirmed: ALSEAMAR's GLIMPSE processing runs *server-side* on the
same raw telemetered data, writing `AD2CP_*_c` into its CSV exports — so the
open-vs-proprietary comparison above is on identical input by construction.
`examples/realtime_vs_delayed.jl` was renamed `examples/realtime_onboard.jl`
accordingly.

## The DAC's hidden dependency (2026-07-15): the onboard flight model runs 5–15 % fast — DAC is now ADCP water-tracked

The question (raised by DG): the DAC — the inverse's absolute reference — is
`(fix − DR_end)/T`, and `DR_end` comes from ALSEAMAR's onboard dead-reckoning
flight model. Is that model any good? Answered by a three-way comparison of
horizontal through-water velocity on all four missions:

1. **onboard** — the logged DR track itself, differentiated (the nav file
   records dead-reckoned positions at ~20 s while `DeadReckoning == 1`);
2. **our model** — GliderTurbulence.jl's steady flight model with the *fixed*
   pooled `FLIGHT_SEA064` polar (nav pitch + depth only; the polar pooled
   M37/M38/M59, making M48 out-of-sample);
3. **the referee** — the ADCP near-cell (4–16 m) mean relative flow, negated:
   a direct measurement of through-water velocity requiring no flight model
   (`throughwater_velocity`). This referee is immune to the recursion concern
   (our polar's ADCP provenance) and to advection.

**Speed** (median ratios to the referee; 33k–93k steady deep DR steps/mission):

```
              onboard/ADCP   flight-model/ADCP
    M37          1.124            1.000
    M38          1.128            1.001
    M48          1.048            1.012      ← out-of-sample polar
    M59          1.147            1.025
```

The onboard model behaves like a steady polar with roughly **half the real
angle of attack** (implied α ≈ 3.5° at 18° pitch where ADCP and polar say
≈ 6°), hence ~12 % too fast horizontally; it is not a fixed-speed law
(p10–p90 of DR speed 0.23–0.47 m/s tracks flight state). M48's onboard DR is
additionally *noisy* (ratio IQR 0.72–0.91), not merely biased. **Direction is
fine** — DR course matches the ADCP course to ≤ 2.3° median everywhere.

**Heading-convention discovery (and a self-caught error).** The DR course
matches the raw gli `Heading` on every mission *including* M59 where
declination is −15.3° — because the gli `Heading` column is **true heading**
(ALSEAMAR applies the configured declination onboard). Verified directly:
nav heading − AD2CP compass = −14.8° (M59) and +1.1° (M38) vs IGRF −15.3°/+1.2°.
The first version of this analysis added declination to nav heading (assuming
magnetic) and produced a spurious −7 cm/s cross-track "error" on M59 — caught
because the flight-model-free ADCP referee showed ~0 cross-track. Recorded in
QA/QC §7: never add declination to nav heading; the AD2CP's own compass IS
magnetic and `process_pings` correctly adds declination there. Nothing in the
production pipeline consumed nav heading, so no product was affected.

**DAC impact.** Water-track DAC = `(fix_after − fix_before − ∫u_tw dt)/T`
(Todd-style water-track referencing; implemented as `compute_dac(nav, pings)`
with per-yo fallback to the onboard estimate below 85 % ADCP coverage):

```
           median |Δdac|  p90     along-track   median |DAC|
    M37       4.1 cm/s     5.1      +4.0 cm/s      9.7 cm/s
    M38       3.9          4.9      +3.8          15.7
    M48       6.5         12.3      +2.0 (noisy)   4.4    ← Δ exceeds the DAC itself
    M59       4.0          4.6      +3.9          22.0
```

Positive along-track = the onboard-referenced DAC was biased *against* the
direction of travel (the fast-model signature), sign-flipping with track
direction — zigzag artifacts between opposing transects, 2–6× the textbook
1–2 cm/s DAC accuracy, propagating 1:1 into the inverse's absolute reference.

**Acid test** — production chain run twice per mission, inverse 0–30 m bins vs
GPS surface drift, per-yo aggregated, paired. (Metric note: per-yo *vector*
medians; not comparable to the earlier "med|Δ| = 4 cm/s (M38)" §9 health-table
figure, which is a different construction.)

```
            med|Δ| onboard → water-track   yos improved   drift Δ along-track (onboard → water-track)
    M37        10.2 → 9.2 cm/s                 73 %           +6.9 → +3.3 cm/s
    M38        16.9 → 16.1                     59 %           +7.6 → +4.0
    M48        13.6 → 10.1                     76 %           +7.0 → +6.0 (bias vector collapses)
    M59        22.8 → 23.6                     28 %           −7.4 → −10.7
```

On M37/M38/M48 the onboard-DAC drift disagreement sits **along-track with
cross ≈ 0 — precisely the DAC-error signature** — and the water-track form
removes almost exactly its predicted share (M48's gain is mostly the per-yo
*random* DR error cleaned up). The residual +3–6 cm/s along-track that remains
on all missions is physically expected (a surfaced hull in the top metre vs a
0–30 m bin average: windage + near-surface shear), not DAC error.

**M59: the referee is broken, not the correction.** DG field context: the
glider was in the Gulf Stream for parts of M59, pushed backwards and sideways
— confirmed, 29 of 154 yos moved net-backwards over ground. The predicted
localization however FAILED: splitting by |DAC| shows the same ~28 % improved
and the same +10 cm/s eastward drift-minus-inverse offset in stream and
quiescent subsets alike (the stream only adds scatter — rms 55 cm/s at
|DAC| ≥ 0.5 m/s). So the drift referee is biased mission-wide (windage/Stokes
on the surfaced hull under summer southwesterlies + wind-driven shear above
the top bins, ~10 cm/s, geographically fixed eastward), and because the pilot
pointed the glider upstream for long stretches the through-water tracks
cluster — a fixed geographic bias then projects as a spurious *negative*
along-track signal (the wrong sign for a DAC error) and mechanically penalizes
any along-track correction. A referee with a 10 cm/s systematic error cannot
adjudicate a 4 cm/s question. M59's correction stands on the direct
measurements: onboard ×1.147 fast vs ADCP (advection-immune), flight model
independently agrees with the ADCP to 2.5 %, window-insensitivity 0.12 cm/s
(4–8 vs 4–16 m — also rules out offset-shear contamination), and the onboard
firmware shows the identical too-fast signature on all four missions.

**Verdict (standing decision).** `compute_dac(nav, pings)` is the production
default — the onboard flight model is out of the product loop. Examples,
tutorial and README updated; the realtime-telemetered product water-tracks its
DAC ashore too (`max_gap = 90 s` for the ~30 s cadence — cells 3–6 lie inside
the 4–16 m window at the fleet 2 m/0.7 m configuration). Fallback rows are
flagged (`method`, `coverage`); M38's duty cycling leaves 64 of 190 yos on the
onboard estimate. The scope note "a biased dead-reckoning model shifts both
products identically" is retired — the residual error budget of the reference
is GPS accuracy plus the ≲ 1 cm/s cell-offset shear term; independent-
instrument validation (mooring / shipboard ADCP) remains the open gap.
Acceptance: synthetic truth-recovery tests plus a gated M38 test pinning the
along-track correction (`test/runtests.jl`, 382 pass).

## Flight-model DAC fallback quantified; the flight model twinned into GliderADCP (2026-07-15, later the same day)

Follow-up question (DG): with no ADCP aboard — or the instrument off — does the
GliderTurbulence flight model still beat the onboard DR for the DAC? Rerunning
the forensic with the heading convention fixed (nav heading used as is — the
first pass double-corrected and inflated the flight model's apparent
cross-track error) and scoring both dead-reckonings against the ADCP water
track as truth:

```
           onboard error (DR−ADCP)          flight-model error (FM−ADCP)
           med|Δ|   along-track             med|Δ|   along-track
    M37    4.3      +4.0 cm/s               1.4      +0.1 cm/s
    M38    4.0      +3.8                    1.4      +0.0
    M48    6.4      +2.0                    3.1      +0.5
    M59    4.0      +3.9                    1.4      +0.8
```

**The flight model is ~3× closer and, decisively, unbiased** — the systematic
anti-track artifact vanishes; what remains is per-yo scatter (p90 2.3–5.9 cm/s;
zero-vertical-water-velocity assumption, apogee masking, residual polar error).
Coverage is every yo on every mission (nav pitch + depth only), including
M38's 64 duty-cycled yos. Caveat: these figures use the same-glider pooled
polar (`FLIGHT_SEA064`); an uncalibrated airframe on published presets lands
between this and the onboard error, still direction-unbiased.

**Decision (DG): the flight model is deliberately TWINNED** — a full copy now
lives in GliderADCP (`src/processing/flightmodel.jl`: `FlightParams` presets,
`solve_aoa`, `flight_model`, `GliderFlight`, plus native `measure_aoa` from
`throughwater_velocity` and `fit_flightparams`, so per-mission polar
calibration no longer needs the package extension) so each package stands
alone. Twin agreement verified to machine precision; physics fixes must land
in both files. Dual loaders qualify the shared names; GliderTurbulence's ADCP
extension now imports its own explicitly.

**The DAC is now a per-yo ladder** — `compute_dac(nav, pings;
fallback=flight_model(nav))`: ADCP water track → flight model → onboard
estimate, flagged in `method`; `compute_dac(nav, flight_model(nav))` is the
ADCP-less form. Examples all use the full ladder (M38: 126 ADCP + 64
flight-model + 0 onboard). Tests 413.

**Diagnostic promoted (same day):** `examples/dac_methods.jl` generates the
three-DAC section comparison per mission and route
(`M*_dac_methods_{delayed,telemetered}.png`): U/V under ADCP water-track /
flight-model / onboard-DR DACs plus differences vs the ADCP reference. The
differences are depth-uniform per yo (barotropic, as they must be); the
onboard−ADCP panels show the track-correlated stripes — most strikingly M59's
sustained +4–5 cm/s eastward block from Jul 27–Aug 13 while the glider pointed
upstream into the Gulf Stream (the anti-track bias never averaging out), and
M48's both-sign ~10 cm/s noisy-DR stripes. The telemetered figures reproduce
the delayed ones nearly bin-for-bin.
