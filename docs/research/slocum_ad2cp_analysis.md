# Slocum-AD2CP (JGradone) — implementation-grade analysis (for the Julia reimplementation)

> Provenance: produced 2026-07-06 by full source analysis of
> https://github.com/JGradone/Slocum-AD2CP at commit `c16712a`.
> Core code: `src/slocum_ad2cp/make_dataset.py` (1137 ln), `analysis.py` (204 ln);
> canonical driver: `notebooks/02_Slocum_AD2CP_Processing_Example.ipynb`.

## A. Metadata

- Purpose: Nortek AD2CP (glider-mounted, down-looking, 4-beam, 1 MHz) on Slocum gliders →
  absolute ocean-velocity profiles via least-squares inversion constrained by the glider's
  dead-reckoned depth-averaged current (DAC).
- Author: Joe Gradone (Rutgers); contributors L. Engdahl, N. von Oppeln-Bronikowski
  (source of the `NicolaiFunctions` heading/AHRS/beam2enu variants). MIT license.
  Citation: Zenodo DOI 10.5281/zenodo.7416126.
- Methods per: Visbeck (2002) Eq. 19; Todd et al. (2012, 2017); Fischer & Visbeck (1993).
- Publication built on this code: Gradone et al. (2023), "Upper Ocean Transport in the
  Anegada Passage From Multi-Year Glider Surveys", JGR Oceans, doi:10.1029/2022JC019608
  ("beams 1, 2, and 4 on a dive and beams 2, 3 and 4 on a climb").

## B. Inputs

- **AD2CP**: MIDAS-exported netCDF only. `load_ad2cp` tries group `Data/Average/` then
  `Data/Burst/`; sorts by time; copies `Config` group attrs from the FIRST file;
  renames `"Velocity Range"→"VelocityRange"` etc.; `Depth = -gsw.z_from_p(Pressure, mean_lat)`;
  transposes to (VelocityRange × time). beam2xyz attr precedence:
  `burst_beam2xyz → beam2xyz → avg_beam2xyz`, reshaped (4,4).
- **Glider**: ERDDAP tabledap (slocum-data.marine.rutgers.edu): `depth, latitude, longitude,
  time, source_file, m_water_vx, m_water_vy [m/s], m_heading [rad], m_gps_mag_var [rad]`.
  Converted to degrees at load.
- **Segmentation**: by Slocum `source_file` (≈ surfacing-to-surfacing, dive+climb);
  AD2CP subset by segment time window; skip segments with max depth < 10 m.
  DAC = LAST non-NaN `m_water_vx/vy` in segment; mag_var = nanmean over segment.
  No record-level glider↔ADCP interpolation: ADCP's own pressure/attitude used throughout.
- Historical deployment-specific fixes (2020/2021 RU29 only): Roll −180°, Pressure −10 dbar.

## C. Pipeline (canonical order)

load_ad2cp → correct_ad2cp_heading (whole deployment) → per segment:
DAC/declination extraction → mag_var_correction (DAC rotation) →
mag_var_correction_ad2cp_ds (heading) → correct_sound_speed → qaqc_pre_coord_transform →
beam_true_depth → binmap_adcp → calcAHRS → beam2enu → qaqc_post_coord_transform →
inversion → per-segment CSV → deployment netCDF grid.

### Magnetometer soft/hard-iron correction — `correct_ad2cp_heading` (1045-1110)
Pitch-binned (arange(−35,35,1), bins with >9 samples): 9-parameter ellipsoid fit
(D = [x²,y²,z²,2xy,2xz,2yz,2x,2y,2z], lstsq vs 1; center = −A₃ₓ₃⁻¹·v[6:9]) → subtract center;
revert bin if |new center x or y| > 150. Recompute heading:
tilt matrix `P = [[cos p, −sin p sin r, −cos r sin p],[0, cos r, −sin r],[sin p, sin r cos p, cos p cos r]]`;
downward orientation: negate Hy,Hz; `h_earth = P·h_inst`; `heading = atan2(h_e[1], h_e[0])` (+360 if <0).

### Declination — `mag_var_correction` (583-590)
`u' = u·cos(mv) − v·sin(mv); v' = u·sin(mv) + v·cos(mv)` (CCW rotation, degrees→radians);
applied to DAC. Heading: `CorrectedHeading_MagVar = CorrectedHeading − mag_var`.

### Sound speed — `correct_sound_speed` (203-209)
`VelocityBeamN *= SpeedOfSound/1500` (Fischer & Visbeck 1993; 1500 = assumed default).

### Pre-transform QC — `qaqc_pre_coord_transform` (214-246)
`corr < 50 % → NaN`; `amplitude > 75 dB → NaN`. That's all (no sidelobe, no ambiguity,
no despiking, no BT anywhere in the package).

### Per-beam vertical cell offsets — `beam_true_depth`/`cell_vert` (62-198)
Beam 1 fwd, 2 port, 3 aft, 4 stbd; beams 1&3 at 47.5°, 2&4 at 25° from vertical.
Vertical offsets below instrument (equivalent forms):
```
beam1: Vr·cos(47.5° + pitch)·cos(roll)      beam2: Vr·cos(pitch)·cos(roll + 25°)
beam3: Vr·cos(47.5° − pitch)·cos(roll)      beam4: Vr·cos(pitch)·cos(roll − 25°)
```
(NOTE: does NOT divide ranges by cos 25° first — cf. the Nortek firmware range-gating fact in
gliderad2cp_analysis.md §C; a physical implementation should reconcile the two conventions.)

### Bin mapping — `binmap_adcp` (109-154)
Per ping/beam: `np.interp(VelocityRange, TrueDepthBeamN, VelocityBeamN, right=NaN)` —
re-grid beam velocities onto the fixed VelocityRange grid interpreted as vertical offsets.
np.interp clamps ABOVE the shallowest cell (constant extrapolation at top), NaN below deepest.

### XYZ→ENU matrices — `calcAHRS` (630-685)
`hh = heading−90°`: `H = [[cos hh, sin hh, 0], [−sin hh, cos hh, 0], [0,0,1]]`; P as above; `R = H·P`.
**BUG**: stored with reshape order='F' but read back order='C' → applies Rᵀ. The validated
behavior (NicolaiFunctions + paper) is `enu = (H·P)·xyz`.

### Beam→XYZ→ENU — `beam2enu` (691-814)
- 3-beam selection: pitch<0 (dive): beams {1,2,4} = matrix columns [0,1,3]; pitch≥0 (climb):
  beams {2,3,4} = columns [1:4]; rows X,Y,Z1 only.
  **BUG (v2.0.0, :793)**: an unconditional `beam2xyz_mat = beam2xyz[0:3,1:4].copy()` after the
  if/elif overwrites the dive-case matrix — dives get the climb matrix. Validated behavior:
  per-case submatrix (on a fresh copy each ping — the old in-place sign-flip-on-a-view bug).
- Downward mount: negate rows Y and Z1 of the 3×3.
- Correct target math (with the factory matrix, after sign flips):
  dive: `X = 0.6782·b1, Y = 1.1831·(b2−b4), Z = −0.74·b1`;
  climb: `X = −0.6782·b3, Y = 1.1831·(b2−b4), Z = −0.74·b3`.
- Hardware-AHRS units: MIDAS provides `AHRSRotationMatrix(9,time)` directly — use when present.

### Post-transform QC — `qaqc_post_coord_transform` (251-283)
|U,V,W| > 0.75 m/s → NaN; first range bin always NaN; all bins NaN where glider depth ≤ 5 m.

### THE INVERSE — `inversion(U, V, dz, u_dac, v_dac, bins, depth, wDAC, wSmoothness)` (289-509)

Model: measured ENU velocity(cell jj, ping ii) = U_ocean(z of cell) − U_glider(ii).
Unknowns m = [U_glider(1..nt) ; U_ocean(1..nz)], complex (U + iV solved jointly).

1. Prune all-NaN rows (bins) and columns (pings).
2. `Z[jj,ii] = depth[ii] + bins[jj]` (absolute cell depths); `ZmM = nanmax(depth)`.
3. Output bins: edges `arange(0, floor(max Z), dz)`, centers `edges[:-1]+dz/2`; count samples
   per bin (strict inequalities); `depth_ind` = #bins deeper than ZmM; trim empty TOP bins.
4. G (sparse, nd = nbin·nt rows): row r = nbin·ii+jj has `G[r, ii] = −1` (glider) and
   `G[r, nt+argmin|bin_new − Z[jj,ii]|] = +1` (nearest ocean bin). Data d = U.flatten(order='F').
5. DAC constraint row: `[zeros(nt); 0; dz·ones(nz−1−depth_ind); zeros(depth_ind)]` —
   i.e. Σ dz·U_ocean over bins 2..(nz−depth_ind), normalized: row scaled by
   C/ZmM with C = 1/‖constraint/ZmM‖₂, weighted wDAC (=5, Todd et al. 2017);
   RHS entry = wDAC·C·DAC. (Ocean-side depth-average form; excludes first bin and
   below-glider bins. No GPS-drift, no bottom-track, no other constraints.)
6. Complexify d = d_u + i·d_v; delete NaN rows from d and G.
7. `obs_per_bin`: nonzero count per ocean column (includes the DAC row — mid-water bins
   count one extra).
8. **BUG**: if last ocean bin has 0 obs: intended column drop is a no-op
   (`Gstar.tocsr()[:,:-1]` result never assigned) while `nt` is incremented → O_ls off by
   one bin (profile shifted one bin deep, first ocean value absorbed as bogus glider velocity).
9. Smoothness (if wSmoothness>0, =1 default): second-difference D2 (rows: −1,+2,−1; last two
   rows truncated, kept); append rows `[0 | D2]` (ocean) and `[D2-shaped (nz×nt) | 0]` (glider —
   only first nz glider velocities smoothed!); RHS zeros.
10. Solve `scipy.sparse.linalg.lsqr(Gstar, d)` (damp=0, atol=btol=1e-6, conlim=1e8).
    `O_ls = x[nt:]` (ocean profile at bin_new), `G_ls = x[:nt]` (glider velocity per ping).

### Shear method — `shear_method` (514-578): **BROKEN in this repo**
Calls `calc_ensemble_shear`, `bin_attr`, `shear_to_vel` — none defined anywhere (NameError).
Recoverable intent: per-ensemble shear → dz-bin average → integrate → reference depth-mean
to DAC; `std_ping = 0.03 m/s`; `vel_std = sqrt(std_ping² + shear_std²)`.

### Outputs
Per segment CSV: `inversion_u/v` (= re/im O_ls), `inversion_depth` (= bin_new),
start/end lon/lat/time, obs_per_bin. Deployment netCDF: u_grid/v_grid(depth,time) —
**assembled by row index, not depth-matched** (misaligned when segments' first bins differ);
`inversion_time` = end time (a missing /2 in the midpoint formula).

## D. Config surface (defaults)

| Item | Value |
|---|---|
| corr_threshold / max_amplitude | 50 % / 75 dB |
| high_velocity_threshold | 0.75 m/s |
| surface_depth_to_filter | 5 m (12 m in some deployments) |
| dz (inversion) | 10 m (1 m in some notebooks) |
| wDAC / wSmoothness | 5 / 1 |
| lsqr | scipy defaults |
| segment min depth | 10 m |
| assumed sound speed | 1500 m/s |
| soft-iron pitch bins | arange(−35,35,1), >9 pts, revert if |center|>150 |
| heading offset in H | −90° |

## E. Bugs & behaviors summary (fix in Julia; keep parity switches for regression)

(a) beam2enu:793 unconditional matrix overwrite (dives use climb matrix);
(b) calcAHRS F-store/C-read → applies (H·P)ᵀ;
(c) inversion empty-last-bin no-op column drop → off-by-one;
(d) "mid" time = end time;
(e) deployment grid by row index, not depth value;
(f) shear_method missing helpers (NameError);
(g) W computed but never used; U/V only in inversion;
(h) np.interp top-clamping in binmap (constant extrapolation above first cell);
(i) glider-block smoothness only covers first nz pings.

Validated (paper-consistent) behavior for the Julia port: dive beams {1,2,4} w/ cols [0,1,3],
climb {2,3,4} w/ cols [1:4], sign-flip rows Y/Z1 on a fresh copy, enu = (H·P)·xyz.
