# Nortek AD2CP formats & Julia ecosystem — research notes

> Provenance: verified web + local-machine research, 2026-07-06. Format details transcribed
> from the Nortek "Integrators Guide — AD2CP" (N3015-007, PDF fetched & extracted) and
> cross-checked against DOLfYN source (github.com/lkilcher/dolfyn → MHKiT-Python fork).

## 1. Julia ecosystem: greenfield confirmed

No registered or unregistered Julia package exists for ADCP or Nortek instrument data
(GitHub + JuliaHub searched). Reference parsers are Python **DOLfYN**/MHKiT and R **oce**
(`read.adp.ad2cp`). Reference glider processors are Python `gliderad2cp` and
`JGradone/Slocum-AD2CP` (both analyzed in sibling docs).

### In-house assets to integrate with (all local on this machine)

| Repo | What it provides |
|---|---|
| `truedichotomy/JLDBDReader.jl` | Pure-Julia Slocum .dbd/.ebd/... reader, validated byte-for-byte vs Python dbdreader; zero non-stdlib deps; the Slocum ingestion path |
| `oceansensing/jlglider` | Script-style Julia toolkit; `seaexplorer/seaexplorerFunc.jl` parses .gli/.pld1 (CSV.jl, delim=';', gz); `ad2cp/process_adcp.py` (2,926 ln, gliderad2cp derivative) is the Python incumbent this package replaces; also `ad2cp/data/sea064_M48.ad2cp.00000.nc` sample |
| `oceansensing/SeaExplorer_Processing` | `AD2CP_processing/`: beam2enu.py, shear_method.py, notebooks, `shear_processing.jl` (partial Julia port to harvest) |
| `oceansensing/ATOMIXjulia.jl` | Active proper Julia package (NCDatasets, GibbsSeaWater deps) — package-convention template |
| `oceansensing/ocean_julia` | Misc utilities (C2PO.jl, EOS.jl, GridFit.jl) |

### Dependency choices (verified registered & healthy)

- **NCDatasets.jl** over NetCDF.jl (first-class netCDF4 **groups** — required for MIDAS files;
  CF handling; lazy DiskArrays reads).
- CSV.jl + DataFrames.jl (SeaExplorer ASCII; gz transparent), CodecZlib.
- Interpolations.jl; GibbsSeaWater.jl (TEOS-10: SA/CT/sound_speed/z_from_p).
- SparseArrays + `qr(A)\b` (SPQR) default solver; **Krylov.jl** `lsqr/lsmr` for large systems
  (IterativeSolvers.jl is in maintenance mode — avoid).
- NaNStatistics.jl (fast nanmean/nanmedian/nanstd); DSP.jl if filtering needed later.
- CairoMakie via **package extension** ([weakdeps]/[extensions], Julia ≥1.9) — core stays light.
- QA: PkgTemplates.jl, Documenter.jl v1, Aqua.jl, TestItems/TestItemRunner.
- Validation-only Python interop: **PythonCall.jl + CondaPkg.toml inside test/**
  (deps never leak to runtime; gate cross-validation behind an env var, e.g. `ADCP_CROSSVAL=1`).

## 2. Nortek .ad2cp binary format (native reader scope)

Authoritative: Nortek "Integrators Guide — Signature/AD2CP", N3015-007
(support.nortekgroup.com; PDF: nortekgroup.com/assets/software/N3015-007-Integrators-Guide-AD2CP.pdf).
All little-endian. A file = sequence of (Header, Data Record).

### Header (10 bytes)

| off | type | field | value |
|---|---|---|---|
| 0 | u8 | sync | 0xA5 |
| 1 | u8 | headerSize | 10 |
| 2 | u8 | ID | record type |
| 3 | u8 | family | 0x10 = AD2CP |
| 4 | u16 | dataSize | bytes in data record |
| 6 | u16 | dataChecksum | |
| 8 | u16 | hdrChecksum | of header bytes 0–7 |

IDs: `0x15` burst, `0x16` average, `0x17` bottom track, `0x18` interleaved burst (beam 5),
`0x1A` burst altimeter raw, `0x1B` DVL BT, `0x1C` echosounder, `0x1D` DVL water track,
`0x1E` altimeter, `0x1F` avg altimeter raw, `0xA0` string record (string ID `0x10` = full
instrument config dump: GETPLAN/GETAVG/GETXFAVG/BEAMCFGLIST… lines). Newer firmware adds
0x23/0x24 raw echosounder records.

**Checksum**: 16-bit wraparound sum of u16 LE words, initialized `0xB58C`; odd final byte
added as `uint16(byte) << 8`. (Julia: `reinterpret(UInt16, ...)` + wraparound `+%` sum.)

### Burst/Average record, "DF3" (version=3; IDs 0x15/0x16/0x18/0x1C)

76-byte fixed prefix (key fields): `version u8`, `offsetOfData u8` (use it, don't hardcode),
**configuration bitmask u16** (bit0 pressure valid, 1 temp, 2 compass, 3 tilt; 5 velocity
included, 6 amplitude, 7 correlation, 8 altimeter, 9 altimeterRaw, 10 AST, 11 echosounder,
12 AHRS, 13 percentGood, 14 stdDev), `serialNumber u32`, clock y/m/d/h/m/s (year−1900,
**month 0-based**) + `microSeconds100 u16`, `soundSpeed u16 ×0.1 m/s`, `temperature i16 ×0.01`,
`pressure u32 ×0.001 dbar`, `heading u16 ×0.01°`, `pitch/roll i16 ×0.01°`,
**beams_cy_cells u16** (bits0–9 nCells, 10–11 coord sys [00 ENU, 01 XYZ, 10 BEAM], 12–15 nBeams;
echosounder records: all 16 bits = nCells), `cellSize u16 mm`, `blanking u16 cm`
(status bit1 → cm scaling; echosounder: mm), `nominalCorrelation u8 %`, `battery u16 ×0.1 V`,
magnetometer i16×3, accelerometer i16×3 (16384 = 1 g), `ambiguityVelocity u16 ×10^velScale m/s`,
`dataSetDescription u16` (nibbles = physical beam of datasets 1–4), `transmitEnergy u16`,
**`velocityScaling i8`** (velocities = raw ×10^k m/s, typically −3), `powerLevel i8 dB`,
temps, `error u16`, `status0 u16`, `status u32` (bits 27–25 orientation [4 ZUP, 5 ZDOWN],
bit16 active plan, bit1 blanking-scaling), `ensembleCounter u32`.

Variable section, strict order, present iff bit set:
velocity `i16[nBeams][nCells]` → amplitude `u8[nB][nC]` (0.5 dB/count) → correlation
`u8[nB][nC]` (0–100 %) → altimeter (f32 m, u16 quality, u16 status) → AST (f32, u16, i16
offset ×100 µs, f32 pressure, 8 spare) → altimeter raw (u32 nSamples, u16 ×0.1 mm,
i16 samples[n] — variable length!) → echosounder `u16[nCells] ×0.01 dB` → AHRS
(f32 rotmatrix M11..M33 row-major; then 4×f32 = quaternion w,x,y,z per DOLfYN [guide calls
them "dummy"]; then f32 gyro x,y,z °/s) → percentGood `u8[nCells]` → stdDev (i16 pitch/roll/
heading ×0.01°, i16 pressure ×0.1 dbar, 12×u16 dummy).

### Bottom-track record, "DF20" (ID 0x17)

Config bits: 0–3 as burst; 5 velocity, 8 distance, 9 figureOfMerit. `ambVelocity u32`;
`error u32`; no status0; `ensembleCounter` at offset 74; then `i32 velocity[nBeams]`
(×10^velScale m/s), `i32 distance[nBeams]` (×0.001 m), `u16 fom[nBeams]`.

### DOLfYN parsing gotchas (encode as requirements)

1. Two-pass with a sidecar index (or mmap single-pass in Julia); header-to-header seeks.
2. Config bitmask + beams_cy asserted **uniform per record ID** across a file; multi-plan
   ("Dual Profile") files → separate datasets; ensemble counters restart — repair by
   monotonicity.
3. Record types interleave (burst/avg/b5/echo/BT/altraw); align by ensemble, not file order.
4. Variable-length altimeter-raw records; peek nSamples, error if it changes.
5. Coordinate-system flag is per record — never assume BEAM.
6. Zeroed timestamps = gaps (interpolate, flag); month 0-based, year−1900.
7. Velocity scaling exponent read per record.

## 3. MIDAS netCDF export conventions

- Groups: `Config` (attributes only; `{plan}_{parameter}` naming: avg_cellSize,
  avg_blankingDistance, avg_coordSystem, avg_beam2xyz (flat 16, rows X,Y,Z1,Z2),
  beamConfiguration{1..4}_theta/phi, pressureOffset, user_decl, rawConfiguration) and
  `Data/Average` (+ `Data/AverageBT` when BT enabled, `Data/Burst` for burst plans).
- `Data/Average` vars: `time` (DateTime), `Velocity Range`/`Correlation Range`/`Amplitude
  Range` (dims with spaces!), `VelocityBeam1..4`, `CorrelationBeam1..4`, `AmplitudeBeam1..4`,
  per-ping `Heading/Pitch/Roll/Pressure/SpeedOfSound/WaterTemperature/Magnetometer*/
  Accelerometer*/Status/Error/Ambiguity/EnsembleCount/...`.
- `Data/AverageBT` vars: per-beam `VelocityBeam1..4`, `DistanceBeam1..4`, `FOMBeam1..4`
  (65535 = invalid) + sensor block.
- Multi-file: `*.ad2cp.00000.nc`, `00001`, … — concatenate on time; Config from first file.
- Nortek publishes no open .ad2cp parser (MIDAS/Ocean Contour are closed); DOLfYN and R oce
  are the de-facto references.
