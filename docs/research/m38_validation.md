# M38 validation notes — GliderADCP.jl vs the prior Python processing

> 2026-07-07. Full-mission run of `examples/m38_currents.jl` (sound-speed correction from
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
