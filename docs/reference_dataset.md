# Reference dataset: sea064 M38 (NorSE Lofoten Basin, Nov 2022 – Mar 2023)

Location: `/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete`

SeaExplorer SEA064, mission 38, Lofoten Basin (~69–70°N, 0–5°E). Nortek Glider AD2CP
SN 102381 (1000 kHz) + RBR Legato CTD + FLBBCD + MicroRider. All facts below were read
directly from the files (instrument config dump, MIDAS netCDF, payload/nav logs).

## 1. AD2CP instrument configuration (from `sea064_M38.cfg` and `.txt` hardware dump)

| Setting | Value | Notes |
|---|---|---|
| Frequency | 1000 kHz | Glider AD2CP, SN 102381 |
| Mode | Average, `MIAVG=10, AVG=1, NPING=4` | 1-s average of 4 pings every 10 s |
| Cells | `NC=15`, `CS=2.00` m, `BD=0.70` m | Velocity Range = 2.7:2:30.7 m (BD + n·CS) |
| Coordinates | `CY="BEAM"` | raw beam velocities stored — all transforms are ours |
| Velocity range (ambiguity) | `VR=2.50` m/s | along-beam; wraps possible (observed ±3.2 m/s) |
| Bandwidth | `BW="BROAD"` (25%) | broadband processing |
| Bottom track | enabled, `RANGE=30` m, `VR=5` m/s, 1 BT ping per 10-s interval | ~24% of BT pings got lock (shelf portions) |
| Salinity setting | `SA=38.0` | wrong for Norwegian Sea (~35.05 from Legato) → onboard sound speed (~1485 m/s) needs re-scaling in post |
| Sound velocity | `SV=0.00` | computed onboard from T + fixed SA=38 |
| Pressure offset | `POFF=9.50` | |
| Declination | `DECL=0.00` | **not set onboard** — must be applied in processing (~+5–8°E at 70°N, 2°E in 2022) |
| Orientation | `ORIENT="AUTOZUPDOWN"` | |

### Beam geometry (`BEAMCFGLIST`)

| Beam | θ (from instrument Z) | φ (azimuth) | Direction |
|---|---|---|---|
| 1 | 47.5° | 0° | forward |
| 2 | 25.0° | −90° | starboard |
| 3 | 47.5° | 180° | aft |
| 4 | 25.0° | +90° | port |

Design intent: at glider pitch ≈ ∓17.5° (dive/climb), the fore (dive) or aft (climb) 47.5°
beam plus both 25° side beams all sit ≈30° from vertical — an effective symmetric 3-beam
Janus array. Beam selection per profile direction is a core processing step.

### Factory beam→XYZ matrix (`GETXFAVG`, rows = X, Y, Z1, Z2)

```
 0.6782   0       −0.6782   0        # X  from beams 1,3:  1/(2 sin 47.5°)
 0       −1.1831   0        1.1831   # Y  from beams 2,4:  1/(2 sin 25°)
 0.7400   0        0.7400   0        # Z1 from beams 1,3:  1/(2 cos 47.5°)
 0        0.5518   0        0.5518   # Z2 from beams 2,4:  1/(2 cos 25°)
```

Compass hard/soft-iron calibration (`GETCOMPASSCAL`), accelerometer calibration, and
pressure-sensor coefficients are also in the `.txt` dump and the netCDF `Config` group.

## 2. AD2CP data products in `ad2cp/102381_sea064_M38/`

| File | Size | Content |
|---|---|---|
| `sea064_M38.ad2cp` | 54 MB | raw instrument binary (Nortek AD2CP record format) — target for the future native Julia reader |
| `sea064_M38.ad2cp.00000.nc` | 133 MB | **MIDAS netCDF export — primary input for processing** |
| `sea064_M38.ad2cp.00000.ntk(.dat/.hdr)` | 84/626 MB | MIDAS intermediate database |
| `sea064_M38.ad2cp.00000_1.mat`, `_2.mat` | 19 MB ea | MATLAB export |
| `sea064_M38.cfg`, `*_deploy.log`, `*_inst.log`, `*_sys.log`, `.txt` | small | deployment config + hardware dump |

### MIDAS netCDF layout (`sea064_M38.ad2cp.00000.nc`)

- Root: no vars; groups `Config`, `Data`.
- `Config`: ~390 attributes — full instrument config incl. `avg_beam2xyz` (16 floats),
  `avg_cellSize=2.0`, `avg_blankingDistance=0.7`, `avg_coordSystem="BEAM"`,
  `beamConfiguration{1..4}_theta/phi`, `pressureOffset=9.5`, `user_decl=0.0`.
- `Data/Average` — dims `time=124752`, `Velocity Range=15` (+ Correlation/Amplitude Range, `Physicalbeam_dim=4`):
  - `time` (DateTime), `Velocity Range` (2.7:2:30.7 m)
  - `VelocityBeam1..4` (15×time, m/s, beam coords), `CorrelationBeam1..4` (%), `AmplitudeBeam1..4` (dB counts)
  - attitude & sensors per ensemble: `Heading`, `Pitch`, `Roll`, `Pressure` (0.4–1003 dbar),
    `WaterTemperature`, `SpeedOfSound`, `MagnetometerX/Y/Z`, `AccelerometerX/Y/Z`,
    `Battery`, `Status`, `Error`, `Ambiguity`, `NominalCorrelation`, `EnsembleCount`, …
- `Data/AverageBT` — dims `time=124751`: per-beam `VelocityBeam1..4` (m/s),
  `DistanceBeam1..4` (m), `FOMBeam1..4` (65535 = invalid), plus the same sensor block.

Mission span in file: 2022-11-03T13:20 → 2023-03-01T19:41 (includes bench time; pitch
median ≈ −16° while diving, 5–95% ≈ [−45°, +24°]).

## 3. Previously processed outputs in `ad2cp/m38_processed/` (validation targets)

Produced with a Python workflow adapted from JGradone/Slocum-AD2CP (inverse method):

| File | Content |
|---|---|
| `M38_ADCP_QAQC_CoordTransformed.nc` | flattened Average group (124004 ensembles × 15 cells) after QC, merged glider nav (`Latitude/Longitude/Depth/DiveNum/date_float`), `newSoundSpeed`, and ENU velocities `VelE/VelN/VelU` |
| `absolute_ocean_vel.csv` | inverse-method output: `yo_number, time_midpoint, u_ocean_vel, v_ocean_vel, depth_bins` (depth_bins negative down, e.g. −106.0) |
| `abs_dive_ocean_vel.csv`, `abs_climb_ocean_vel.csv` | same, dive-only / climb-only solutions |
| `M38_glider_data_processed.csv` | parsed+merged SeaExplorer nav (Timestamp, NavState, Heading, Declination, Pitch, Roll, Depth, …, DeadReckoning, diveNum, Lat_dd, Lon_dd, date_float) |
| `M38_science_data_processed.csv` | parsed payload science data |
| `figures/` | `NORSE_SEA064_M38_lofoten_U.png`, `..._V.png` section plots |

These give end-to-end regression targets: (i) after QC+transform (`VelE/N/U`), and
(ii) final absolute velocity profiles per yo.

## 4. SeaExplorer platform files in `delayed/`

- `nav/logs/sea064.38.gli.sub.N.gz` (N=1..~191, 383 files): semicolon-separated, header
  `Timestamp;NavState;SecurityLevel;Heading;Declination;Pitch;Roll;Depth;Temperature;Pa;Lat;Lon;DeadReckoning;DesiredH;BallastCmd;BallastPos;LinCmd;LinPos;AngCmd;AngPos;Voltage;Altitude;`
  - `Timestamp` = `DD/MM/YYYY HH:MM:SS`; `Lat/Lon` in NMEA degrees·100+minutes (`7001.296` = 70°01.296′)
  - `DeadReckoning` flag: 1 = subsurface dead-reckoned position, 0 = GPS-fixed at surface.
    Observed surfacing: NavState 117 (ascent end) → 116 (surface/GPS, DR 1→0 with position
    jump = accumulated DAC displacement) → 110 (descent). NavState histogram in one file:
    {100: navigating, 110: descending?, 115/117: ascent phases, 116: surface, 118: inflection}
    — semantics to confirm against ALSEAMAR docs.
  - `Declination` column exists but reads 0 in checked records.
- `pld1/logs/` (191 segments each):
  - `sea064.38.pld1.raw.N.gz` — full payload record, header:
    `PLD_REALTIMECLOCK;NAV_RESOURCE;NAV_LONGITUDE;NAV_LATITUDE;NAV_DEPTH;AD2CP_TIME;AD2CP_HEADING;AD2CP_PITCH;AD2CP_ROLL;AD2CP_PRESSURE;AD2CP_ALT;AD2CP_V{1..4}_CN{1..6};FLBBCD_*;LEGATO_*(COND,TEMP,PRES,SAL,CONDTEMP);MR1000G-RDL_*`
    (AD2CP real-time subset = 6 cells × 4 beams; timestamps with ms)
  - `sea064.38.ad2cp.raw.N.gz` — Nortek NMEA ASCII stream:
    `$PNORI,4,Glider102381,4,15,0.70,2.00,2*6D` (config: 4 beams, 15 cells, blank 0.70 m, cell 2.00 m, coord=BEAM(2));
    `$PNORS,MMDDYY,HHMMSS,err,status,battery,soundspeed,heading,pitch,roll,pressure,temp,a1,a2*cs`;
    `$PNORC,MMDDYY,HHMMSS,cell,v1,v2,v3,v4,,,C,amp1..4,corr1..4*cs`
  - `sea064.38.legato.raw.N.gz` — RBR Legato CTD stream
  - `sea064.38.pld1.sub.N.gz` — decimated payload (transmitted subset)

## 5. Implications for the Julia package

1. **Primary pipeline input** = MIDAS netCDF (`Config` + `Data/Average` + `Data/AverageBT`)
   merged with SeaExplorer `gli` nav files (GPS/DR for DAC) and `pld1`/`legato` science
   (salinity → sound-speed correction).
2. **Beam-coordinate data + full config on file** → the package must own the whole chain:
   sound-speed rescale → QC → beam selection (dive/climb) → beam→XYZ→ENU (3-beam solution)
   → bin mapping → shear/inverse solutions → DAC/BT referencing.
3. **Bottom track present** → implement the BT constraint from day one (validation gold).
4. **Real-time ASCII stream** ($PNOR / `AD2CP_*` payload columns) → a lightweight
   real-time parser enables piloting products and cross-checks of the delayed pipeline.
5. **Native `.ad2cp` binary reader** is worthwhile later (54 MB binary vs 133 MB netCDF;
   removes the MIDAS/Windows dependency entirely).
6. **Known reference outputs** allow regression tests at two pipeline stages.
