# Literature review — glider-mounted ADCP velocity processing

> Provenance: verified web literature search, 2026-07-06. Every entry's existence was
> verified against the publisher/mirror; content marked [V] was extracted from the actual
> document, [V-sec] from abstract/secondary restatement (paywalled body).

## 1. Foundational inverse method

- **Visbeck, M. (2002).** "Deep velocity profiling using lowered acoustic Doppler current
  profilers: Bottom track and inverse solutions." *J. Atmos. Oceanic Technol.* 19(5), 794–807.
  doi:10.1175/1520-0426(2002)019<0794:DVPULA>2.0.CO;2. [V-sec; equations verified via the
  U. Bremen LADCP practical restatement]
  - Data model per sample: `U_adcp = U_ocean + U_ctd + noise`.
  - Linear system `d = G·m + n`, unknowns `m = [U_ctd(per ping) ; U_ocean(per depth cell)]`;
    `n_ocean = H/Δz`, Δz typically = bin size; overdetermined but **GᵀG singular without
    constraints** (baroclinic-only otherwise).
  - Barotropic constraint: `∫ U_ctd dt = GPS displacement` — one weighted row
    `[w·dt₁/T, w·dt₂/T, …, 0…]`, RHS `w·U_ship`; large weight.
  - Bottom-track constraint: pins instrument velocity near the seafloor.
  - Smoothness constraint damps vertical roughness. Per-constraint weights = main tuning.
- **Fischer, J. & Visbeck, M. (1993).** "Deep velocity profiling with self-contained ADCPs."
  *JTECH* 10(5), 764–773. [V] Shear-method lineage + sound-speed scaling `v ∝ c`.

## 2. Glider adaptations

- **Todd, R.E., Rudnick, D.L., Sherman, J.T., Owens, W.B., George, L. (2017).** "Absolute
  Velocity Estimates from Autonomous Underwater Gliders Equipped with Doppler Current
  Profilers." *JTECH* 34(2), 309–333. doi:10.1175/JTECH-D-16-0156.1. [V-sec]
  - 1-MHz Nortek AD2CP on Spray; upper 1000 m; LADCP inverse reviewed & extended with
    additional constraints (DAC, surface, flight-model velocities; + smoothness).
  - Near-glider relative velocities infer dive-dependent flight parameters → corrects DAC
    for biofouling drift over long missions.
  - Source of the standard weights used downstream (wDAC=5, wSmoothness=1 in Slocum-AD2CP).
- **Todd et al. (2011).** "Poleward flows in the southern California Current System…"
  *JGR Oceans* 116, C02026. doi:10.1029/2010JC006536. [V] DAC-referencing precursor
  (Spray + Sontek 750 kHz ADP): both ADP velocities and thermal-wind shear referenced to
  the measured depth-averaged current.
- **Rudnick, D.L., Sherman, J.T., Wu, A.P. (2018).** "Depth-Average Velocity from Spray
  Underwater Gliders." *JTECH* 35(8), 1665–1673. doi:10.1175/JTECH-D-17-0200.1. [V]
  Canonical DAC-accuracy study.
- **Gradone, J.C., et al. (2023).** "Upper Ocean Transport in the Anegada Passage From
  Multi-Year Glider Surveys." *JGR Oceans* 128, e2022JC019608. doi:10.1029/2022JC019608. [V]
  - Slocum RU29 with 1-MHz AD2CP (beam coords, 8 Hz pings, 0.2 m blank, 0.5 m bins on that
    config); RU36 with 600-kHz TRDI Pathfinder.
  - Beams 1&3 at 47.5°, 2&4 at 25°; at 24–27° glider pitch only 3 of 4 beams used
    ("beams 1,2,4 on a dive; 2,3,4 on a climb").
  - QC: corr < 50 % or amplitude > 75 dB discarded (Pathfinder: corr<50, echo<70 dB, PG<80 %).
  - **DAC RMS accuracy ≈ 1–2 cm/s** — best-attested absolute-referencing error figure.
- **Ma, W., et al. (2019).** "Absolute Current Estimation and Sea-Trial Application of
  Glider-Mounted AD2CP." *J. Coastal Research* 35(6), 1343–1350.
  doi:10.2112/JCOASTRES-D-18-00176.1. [V] Independent glider-AD2CP inverse implementation.
- **Queste, B.Y., Rollo, C., Font, E., Mohrmann, M. (2025).** "gliderad2cp: A Python package
  to process Nortek AD2CP velocity profiles from gliders." *JOSS*, in review since 2025-03-18
  (reserved doi:10.21105/joss.08342). [V] Shear method + DAC referencing + shear-bias
  regression; SeaExplorer 4-beam; Seaglider/Spray alternating 3-beam unsupported (TODO).
- **Stevens-Haas, J., Webster, S.E., Aravkin, A. (2022).** arXiv:2110.10199. [V] Joint
  probabilistic (Kalman-smoothing) current + navigation estimation — the estimation-theory
  angle; optional future direction.

## 3. LADCP processing tradition (shear method)

- **Thurnherr, A.M. (2010).** "A Practical Assessment of the Errors Associated with
  Full-Depth LADCP Profiles…" *JTECH* 27(7), 1215–1227. doi:10.1175/2010JTECHO708.1. [V]
  Error budget; shear vs inverse comparison.
- **Thurnherr, A.M.** "How To Process LADCP Data With the LDEO Software" (LDEO_IX). [V]
  Implements BOTH shear (Firing/UH lineage) and Visbeck inverse; BT constrains instrument motion.
- **Thurnherr, Visbeck, Firing, King, Hummon, Krahmann, Huber (2010).** GO-SHIP LADCP
  cookbook, IOCCP Rep. 14 / ICPO 134. [V]
- **Thurnherr, Symonds, St. Laurent (2015).** "Processing explorer ADCP data collected on
  Slocum gliders using the LADCP shear method." *IEEE CWTM 2015*.
  doi:10.1109/CWTM.2015.7098134. [V] Directly on point for glider shear processing.
- **Firing, E., Hummon, J.M. (2010).** "Shipboard ADCP Measurements," GO-SHIP manual. [V]
  NOTE: that chapter is SADCP context; attribute the LADCP shear method itself to
  Fischer & Visbeck 1993 / Firing's UH software / the LADCP cookbook.

## 4. Glider flight models (DAC + vertical velocity + method C)

- **Merckelbach, L., Smeed, D., Griffiths, G. (2010).** "Vertical Water Velocities from
  Underwater Gliders." *JTECH* 27(3), 547–563. doi:10.1175/2009JTECHO710.1. [V]
  `w_water = dz/dt − w_glider(model)`; lift coefficient unidentifiable without a direct
  horizontal-speed measurement (exactly what ADCP near-bins provide).
- **Merckelbach, L., et al. (2019).** "A Dynamic Flight Model for Slocum Gliders…"
  *JTECH* 36(2). doi:10.1175/JTECH-D-18-0168.1. [V, equations extracted]
  Steady planar model: `0 = sin(θ+α)F_L − cos(θ+α)F_D`; `0 = F_B − F_g − cos(θ+α)F_L − sin(θ+α)F_D`;
  `F_D = ½ρSU²(C_D0 + C_D1 α²)`, `F_L = ½ρSU² a α`;
  AoA: `α = (C_D0 + C_D1 α²)/(a·tan(θ+α))` (iterate);
  calibration: min Σ[U sin(θ+α) + dh/dt]². Python: smerckel/gliderflight.
- **Welch, T.P., et al. (2022).** "In Situ Calibration of Underwater Glider Flight Model
  Using ADCPs." *JTECH* 39(9). doi:10.1175/JTECH-D-21-0074.1. [V] ADCP near-bin velocities
  → lift/drag regression → better AoA/dead-reckoning/DAC.

## 5. Instrument physics & QC references

- **Nortek, "Signature Principles of Operation"** (N3015-011). [V] Broadband autocorrelation
  processing; σ ∝ 1/√N; cell center = blanking + n·cellSize; per-instrument transformation
  matrix.
- **Nortek/iRobot (2012).** "Improving Depth Averaged Velocity Measurements from Seaglider
  with … the Nortek AD2CP-Glider." *IEEE OCEANS 2012*. [V] Origin of the glider-specific
  47.5°/25° beam geometry.
- **Shcherbina, A.Y., D'Asaro, E.A., Nylund, S. (2018).** "Observing Finescale Oceanic
  Velocity Structure with an Autonomous Nortek ADCP." *JTECH* 35(2), 411–427.
  doi:10.1175/JTECH-D-17-0108.1. [V full text]
  - Range–velocity trade-off `V_R·L_max = c²/(8F₀) ≈ 0.28 m²/s` at 1 MHz;
    ambiguity wraps at ±V_R; EVR/MCPC extends to `2·n_max·V_R` (n_max≈6).
  - Phase-unwrapping options (Itoh 1D → 2D-SRNCP → SNAPHU); correlation → velocity variance:
    `var(v) = (V_R/π)²·var(φ)`, Miller–Rochwarger `var(φ) ≈ R⁻² − 1/(2M)`.
- **Teledyne RDI, "ADCP Principles of Operation: A Practical Primer."** [V]
  - **Sidelobe: `R_max = D·cos θ`** (contaminated outer fraction = 1 − cos θ:
    6 % at 20°, 15 % at 30°; → 9.4 % at 25°, 32 % at 47.5° for AD2CP glider beams).
  - Doppler `F_d = 2 F_s (V/C) cos A`; σ ∝ 1/√N; averaging removes random error, not bias.
- **von Appen, W.-J. (2015).** "Correction of ADCP Compass Errors Resulting from Iron in the
  Instrument's Vicinity." *JTECH* 32(3). doi:10.1175/JTECH-D-14-00043.1. [V]
  Hard-iron = fixed offset; soft-iron = direction-dependent warp; land calibration degrades
  when deployment-site field differs; deviations can exceed 90° near steel.
- **de Fommervault, O., et al. (2019).** "SeaExplorer Underwater Glider: A New Tool to Measure
  depth-resolved water currents profiles." *IEEE OCEANS 2019 Marseille*. [V-sec]
  ALSEAMAR AD2CP integration (2017); onboard shear algorithm; validated vs ship ADCP.
- **ALSEAMAR (2024).** AGU OS abstract: NRT onboard shear-method currents + backscatter index.

## Synthesis

**Method taxonomy.** (A) Shear-then-reference (differentiate → integrate → DAC/BT reference);
(B) Visbeck-style joint inverse (ocean + platform unknowns; DAC/BT/surface/smoothness rows);
(C) bottom-track augmentation (absolute near-seafloor reference); (D) flight-model methods
(produce DAC & w; ADCP near-bins calibrate the flight model in return). A and B consume the
same QC'd relative velocities — implement both over a common trunk.

**Standard QC recipe.** correlation ≥ 50 % (Gradone; 80 % gliderad2cp average-mode);
amplitude window (> noise floor + SNR, < 75–80 dB); sidelobe cut `R_max = D·cos θ` per beam;
error velocity / 4-beam consistency where available; ambiguity screening vs V_R (2.5 m/s for
the M38 config); surface/near-transducer exclusion; compass calibration pre-deploy (+ possible
post-processing correction).

**Error budget.** Compass/heading errors dominate horizontal velocity (degrees → cm/s at
typical relative speeds); DAC accuracy 1–2 cm/s floors the absolute reference; shear bias
accumulates with depth (gliderad2cp regresses it out; Todd 2017 handles flight-parameter
drift); single-ping noise large but averages down as 1/√N; sidelobe contaminates outer
6–32 % of range depending on beam angle.

**Must-cite:** Visbeck 2002; Fischer & Visbeck 1993; Todd et al. 2017; Gradone et al. 2023;
Queste et al. (JOSS, gliderad2cp); Thurnherr et al. 2015 + LADCP cookbook; Merckelbach
2010/2019; Rudnick et al. 2018; Shcherbina et al. 2018; von Appen 2015; Nortek N3015-011;
RDI Primer.
