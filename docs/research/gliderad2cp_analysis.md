# gliderad2cp — implementation-grade analysis (for the Julia reimplementation)

> Provenance: produced 2026-07-06 by full source analysis of
> https://github.com/bastienqueste/gliderad2cp at commit `7fbccaf` (main, post-v0.0.10),
> including an actual run of the full pipeline on the official sample dataset
> (VOTO SEA055 M82, profiles 160–210) with captured ground truth
> (see `validation/gliderad2cp_reference/`). Line numbers refer to that commit.

## A. Package metadata

- Purpose: Nortek Glider AD2CP (1 MHz, 4-beam) → QC'd along-beam velocities → ENU
  velocities on isobars → vertical shear → gridded shear → vertically integrated,
  DAC-referenced absolute velocity → optional shear-bias correction (Todd et al. 2017 variant).
- Authors (CITATION.cff): Bastien Y. Queste (U. Gothenburg), Callum Rollo (VOTO), Estel Font,
  Martin Mohrmann. License MIT. PyPI 0.0.10 (2025-03-18). JOSS paper in review
  (openjournals/joss-reviews#8342; reserved DOI 10.21105/joss.08342).
- Modules: `process_shear.py` (734 ln), `process_currents.py` (429), `process_bias.py` (383),
  `tools.py` (98). Deps: dask, gsw, matplotlib, netcdf4, numpy, pandas, pooch, pyarrow, scipy, xarray.
- Sample data (pooch): https://erddap.observations.voiceoftheocean.org/examples/gliderad2cp/

## B. Expected inputs

### B1. AD2CP netCDF (Nortek MIDAS export)

Opened with `xr.open_mfdataset(path, group='Data/Average')` (process_shear.py:171); Config
group read from the FIRST file only; Config **attributes** copied onto the dataset. Burst-only
files unsupported (issue #82).

Variables actually used:

| Variable | dims | units | use |
|---|---|---|---|
| `time` | (time) | s since 1970 → datetime64 | master time base |
| `VelocityBeam1..4` | (time, Velocity Range) | m/s along-beam | core data |
| `CorrelationBeam1..4` | (time, Correlation Range) | % | QC |
| `AmplitudeBeam1..4` | (time, Amplitude Range) | dB | QC |
| `*Range` coords | 1-D | m | cell ranges = blanking + (n+1)·cellSize |
| `Pressure` | (time) | dbar | glider depth via gsw |
| `Heading/Pitch/Roll` | (time) | deg | attitude; pitch>0 = nose-up = upcast |
| `SpeedOfSound` | (time) | m/s | instrument sound speed (configured salinity) |
| `AccelerometerZ` | (time) | g | mounting auto-detect (median>0 → 'top') |

Unused but present: Ambiguity (!), Status, Error, Magnetometer*, Physicalbeam, WaterTemperature, etc.
Config attrs used: `avg_cellSize`, `avg_blankingDistance` only. `avg_beam2xyz` available but the
code reproduces it analytically. `user_decl` ignored.

After load, the three Range dims are collapsed to one integer dim `bin` (0..n-1).

### B2. Glider timeseries

`load_data(adcp_file_path, glider_file_path, options)` (process_shear.py:61). Accepts DataFrame /
xr.Dataset / .csv (actually TAB-separated — issue #85) / .pqt / .nc. Required columns:
`time, temperature, salinity, latitude, longitude, profile_number, pressure` (`dive_num` optional;
everything else dropped — including any `declination` column).

- Sound speed: `gsw.sound_speed(SA_from_SP(sal,p,lon,lat), CT_from_t(SA,T,p), p)`.
- Time merge: 1-D linear interpolation of glider vars onto ADCP time in float ns
  (`tools.interp`: scipy interp1d, NaN-dropped, fill NaN). `profile_number` is np.round-ed
  after interpolation (issue #70: can propagate artifacts).
- `Depth = -gsw.z_from_p(Pressure, Latitude)` (positive down).
- GPS pre/post-dive arrays for process_currents are the USER's responsibility
  (Mx3: datetime64[ns], lon, lat; SeaExplorer recipe = dead_reckoning flag edges).

## C. Pipeline

### Stage 1 — `process_shear.process` (669-735)

Order: load_data → (compass correction placeholder — NOT implemented) →
_velocity_soundspeed_correction → _quality_control_velocities →
_determine_velocity_measurement_depths → _regrid_beam_velocities_to_isobars →
_rotate_BEAMS_to_XYZ → _rotate_XYZ_to_ENU.

**Sound speed** (219-248): `VelocityBeam_b *= c_glider/c_instrument`; renames
`SpeedOfSound → NoSal_SpeedOfSound`, `glider_soundspeed → SpeedOfSound` (idempotence via rename).

**QC** (251-307):
- `noise_floor` = 0.5th percentile of ALL amplitudes, all 4 beams pooled.
- correlation < 80 % → NaN; amplitude > 80 dB OR < floor+3 dB → NaN; |v_beam| > 0.8 m/s → NaN.
- Masks applied multiplicatively (1/NaN). NO sidelobe mask, NO despiking, NO pitch/roll mask.

**Measurement depths** (313-421). With H,P,R = attitude + scalar offsets; per-beam angle from
vertical (radians), top-mounted (`direction=+1`):
```
θ1 = arccos(cos(47.5°−P)·cos R);  θ2 = arccos(cos(25°−R)·cos P)
θ3 = arccos(cos(47.5°+P)·cos R);  θ4 = arccos(cos(25°+R)·cos P)
```
bottom-mounted: swap P signs in θ1/θ3 and R signs in θ2/θ4.
- **Nortek firmware fact** (comment 367-385, pers. comm. Sven Nylund 2024-11-27): when >2 beams
  ping, the instrument time-gates ALL beams with a hard-coded nominal 25° slant and 1500 m/s:
  along-beam cell distance = reported Range / cos(25°). Not recorded in file attributes;
  firmware-dependent.
- Cell depth per beam: `D_b[t,bin] = Depth[t] − direction · (Range/cos25°)[bin] · cos(θ_b[t])`.

**Isobaric regridding** (424-502): `depth_offsets = arange(0, blanking + (n+0.5)·cellSize + cellSize, cellSize/2) · direction`
(sample: 0:0.5:31.5, 64 pts). Per beam, per ping: 1-D linear interp of V_b(bin) from
x = `Depth − D_b` onto depth_offsets (scipy, fill NaN, interpolation bridges QC gaps —
no max-gap limit). Adds coords `depth_offset(gridded_bin)`, `bin_depth(time, gridded_bin)
= Depth − depth_offsets`. NOTE: the `gridded_bin` dim itself has NO coordinate (see G1).

**BEAM→XYZ** (505-570): constants `a(t)=1/(2 sin t)`, `b(t)=1/(4 cos t)`, tf=47.5, ts=25.
- upcasts = (Pitch+offset) > 0; downcasts otherwise.
- 3-beam trick: the more-horizontal fore/aft beam is REPLACED by synthesis from the other three
  via zero-error-velocity: `replaced = 0.745429·(V2+V4) − good` (0.745429 = b(ts)/b(tf)).
  Top-mounted: replace V1 on downcasts (using V3), V3 on upcasts (using V1). Bottom-mounted: swapped.
  In-place mutation: stored V1/V3 contain synthetic values.
- `X = a(tf)·(V1 − V3);  Y = a(ts)·(−V2 + V4);  Z = 2·b(ts)·(V2 + V4)`
  (Y sign flagged "TODO: sign uncertainty" at :563).

**XYZ→ENU** (573-663): with hh=(H−90)·π/180, pp, rr:
```
M = [ cos hh·cos pp,  −cos hh·sin pp·sin rr + sin hh·cos rr,  −cos hh·sin pp·cos rr − sin hh·sin rr
     −sin hh·cos pp,   sin hh·sin pp·sin rr + cos hh·cos rr,   sin hh·sin pp·cos rr − cos hh·sin rr
      sin pp,          cos pp·sin rr,                           cos pp·cos rr ]
E = M00·X + M01·Y·dir + M02·Z·dir;  N = M10·X + M11·Y·dir + M12·Z·dir;  U = M20·X + ...
```
(Y,Z sign-flipped for bottom mount; X not.) E/N/U = water velocity relative to glider, earth frame.
Shear: `Sh_E = E.differentiate('gridded_bin')` etc. — per INDEX, not per meter (G1).

### Stage 2 — `process_currents.process(ADCP, gps_predive, gps_postdive, options)`

Grid: `xi` = profile_number bin edges (default width 1), `yi` = depth bin edges
(default width = avg_cellSize).

**get_DAC** (55-148) — DAC via "ADCP as DVL":
- deltat = diff(time); gaps > 10× median → excluded.
- Glider through-water displacement: `travel_e(t) = cumsum(−mean_bins(E)·deltat)` (interpolant), same N.
- Per dive (pre-fix i → earliest post-fix j): `dxy_dvl = travel(t_j) − travel(t_i)`;
  `dxy_gps` from Δlon/Δlat × meters-per-degree (gsw.distance at pre-dive position);
  `dac = (dxy_gps − dxy_dvl)/duration`; `dive_time` = midpoint.
- Sample ground truth: dac[0] = [−0.02708, +0.10196] m/s, duration 6843 s.

**_grid_shear** (150-257): `tools.grid2d` = pd.cut(right=False → left-closed; last edge excluded,
issue #80) + groupby agg. Grids: time, time_in_bin (Σ ensemble durations by GLIDER depth;
durations zeroed for Depth<1 m or (Depth<3 m & |dz/dt|<0.04)), heading_N/E (mean cos/sin),
speed_through_water, shear_E/N mean/median/stddev/count (by measurement bin_depth).

**_grid_velocity** (259-310): integrates the MEDIAN shear grid:
`V = cumsum(nan_to_num(Sh)·gradient(depth)) ; V[!finite(Sh)] = NaN ; V −= nanmean(V, depth)`
→ `velocity_E/N_no_reference` (zero-mean baroclinic; cumsum holds velocity constant across NaN gaps).

**_reference_velocity** (312-365): `DAC_D(profile)` = time-interp of per-dive DAC to profile times;
time-in-bin weighted referencing:
`w̄ = mean_depth(w); reference = DAC − nanmean(V_norefer·w/w̄, depth); V_ref = V_norefer + reference`.
(Deliberate: depth-mean of referenced profile ≠ DAC when sampling is non-uniform in time.)

### Stage 3 — `process_bias.process` (optional)

- `_linear_regression`: filter |mean shear| < 0.01 s⁻¹ & displacement > 100 m; optional
  displacement weighting; weighted normal equations; returns slope.
- `regress_bias`: Nelder-Mead (fmin, x0=[0,0], maxiter=300, xtol=ftol=1e-9) minimizing
  Σ|slope(displacement_l1 vs mean_depth(d/dz V_l2_corrected))|·1e9 over 4 combos,
  depth slice default (0,1000).
- `correct_bias`: bias velocity = cumsum over depth of `speed_through_water·(heading·b)` terms
  (velocity-dependent variant) − depth mean; subtracted from referenced velocities.
  Sample: along = −1.04e-3, across = −5.77e-4.
  NOTE the across-glider direction vector used is (−h_N, −h_E) — not orthogonal to heading;
  reproduce verbatim in parity mode.

## D. Options (tools.get_options) and defaults

| Option | Default |
|---|---|
| correct_compass_calibration | False (not implemented) |
| shear_to_velocity_method | 'integrate' ('lsq' placeholder, unimplemented) |
| ADCP_mounting_direction | 'auto' (AccelerometerZ median sign) |
| QC_correlation_threshold | 80 (%) |
| QC_amplitude_threshold | 80 (dB) |
| QC_velocity_threshold | 0.8 (m/s) |
| QC_SNR_threshold | 3 (dB above pooled 0.5th-pct noise floor) |
| velocity_regridding_distance_from_glider | 'auto' (0 : cellSize/2 : max, ×direction) |
| xaxis / yaxis | 1 profile / avg_cellSize m |
| weight_shear_bias_regression | False |
| velocity_dependent_shear_bias_correction | False (author uses True) |
| shear_bias_regression_depth_slice | (0, 1000) |
| pitch_offset / roll_offset / heading_offset | 0 / 0 / 0 (heading_offset = only declination hook) |

## E. Verified traps for the reimplementer

1. **Per-index shear (the big one).** `differentiate('gridded_bin')` uses unit index spacing
   (no coord on that dim) — empirically identical to `np.gradient(E, axis=bin)`. One index
   = cellSize/2 m of offset; for top-mounted, bin_depth DECREASES with index. Stage-2 integration
   multiplies by the yaxis spacing (default cellSize). Net (verified by injecting synthetic
   E = 0.01·bin_depth): recovered baroclinic slope = −(cellSize/2)·truth for top-mounted
   (cellSize=1 sample); +(cellSize/2)·truth for bottom-mounted. **Exactly correct only for
   bottom-mounted (down-looking) with cellSize = 2 m** — the standard Nortek glider config —
   which is presumably why it went unnoticed. A physical implementation must use signed
   dz per index = −direction·(cellSize/2); keep a parity flag for cross-validation.
2. **Nortek 25° range-gating hard-coding** (see C, measurement depths). Along-beam = Range/cos 25°.
3. Beam bookkeeping: discarded fore/aft beam swaps with cast direction AND mounting;
   upcast := pitch+offset > 0; V1/V3 mutated in place.
4. `correct_bias` mask is a NO-OP (bool array assigned NaN → all True); bias profile is dense.
5. `depth.differentiate('depth') ≡ 1`; bias cumsum is per-bin; xarray cumsum skips NaN.
6. grid2d bins are LEFT-closed with left-edge labels; last row/col never populated (issue #80).
7. Everything time-interpolated in float NANOSECONDS; DAC math in seconds; +1e-7 s fudge in
   glider time conversion.
8. DVL displacement integrand zeroes dt gaps > 10× median dt — essential for surfacings.
9. No declination anywhere (issue #76), no sidelobe mask, no ambiguity handling
   (Ambiguity variable ignored; only |v|<0.8 m/s cap), no bottom track (issue #74 wants
   Gradone inversion port).
10. Public API: `process_shear.process`, `process_currents.process`/`get_DAC`,
    `process_bias.process`, `tools.get_options`.

## F. Ground truth for cross-validation

See `validation/gliderad2cp_reference/README.md` (inputs + outputs + checkpoint numbers).
