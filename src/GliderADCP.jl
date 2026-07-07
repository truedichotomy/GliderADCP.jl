"""
    GliderADCP

Pure-Julia processing of glider-mounted Nortek AD2CP data (SeaExplorer and Slocum)
into absolute ocean velocity profiles.

Implements, from first principles, the two published approaches:

  * the lADCP-style **shear method** (grid shear → integrate → reference to
    depth-averaged current), following Todd et al. (2017) and the `gliderad2cp`
    Python package (Queste et al.);
  * the Visbeck (2002) **least-squares inverse method** adapted to gliders
    (Todd et al. 2017; Gradone et al. 2023, `Slocum-AD2CP`), with optional
    bottom-track and surface-drift constraints.

Processing is organized as independent layers (see `PLAN.md`):

  I/O (MIDAS netCDF, SeaExplorer gli/pld, \$PNOR stream, Slocum, native .ad2cp)
  → sound-speed correction → QC → beam geometry (3-beam selection, beam→XYZ→ENU)
  → bin mapping → DAC → velocity solutions (shear | inverse) → gridded products.
"""
module GliderADCP

using Dates
using Statistics
using LinearAlgebra
using SparseArrays
using Printf
using Logging

using NCDatasets
using CSV
using DataFrames
using CodecZlib
using Interpolations
using GibbsSeaWater
using NaNStatistics
using Krylov
using SatelliteToolboxGeomagneticField: igrfd

export AD2CPConfig, AD2CPData, BottomTrackData, GliderNav
export load_ad2cp, load_seaexplorer_nav, load_seaexplorer_pld, seaexplorer_files
export nmea2deg, ncells
export beam_unit_vectors, xyz_from_beams, select_beams, head2vehicle,
       detect_look_direction, rotmat_xyz2enu, beams_to_enu
export QCThresholds, qc!, bt_valid
export soundspeed_from_ctd, soundspeed_correction, apply_soundspeed!
export vertical_cosines, offset_grid, regrid_beams, enu_on_isobars
export compute_dac, surface_drift, lonlat_to_dxdy
export ProcessedPings, process_pings, glider_depth, segment_indices, bt_velocity
export InverseOptions, invert_segment, solve_inverse
export ShearOptions, shear_segment, integrate_shear, solve_shear
export magnetic_declination, grid_profiles, export_sections
export load_pnor, slocum_nav, dac_from_slocum, time_in_bin, plot_sections

"""
    plot_sections(panels; colorrange=(-0.5, 0.5), colormap=:balance, figsize=...)

Section figure from gridded profiles: `panels` is a vector of
`(section, field::Symbol, title)` tuples, e.g.
`[(sec, :U, "U — inverse"), (sec, :V, "V — inverse")]` where `section` comes from
[`grid_profiles`](@ref). Requires a Makie backend to be loaded
(e.g. `using CairoMakie`) — implemented in the package extension.
"""
function plot_sections end

# ---- Layer 0: types & configuration ------------------------------------------------
include("types.jl")

# ---- Layer 1: I/O -------------------------------------------------------------------
include("io/nortek_netcdf.jl")   # MIDAS-exported .nc  (primary input)
include("io/nortek_pnor.jl")     # real-time \$PNORI/\$PNORS/\$PNORC ASCII stream
include("io/seaexplorer.jl")     # .gli / .pld1 navigation & payload files
include("io/slocum.jl")          # Slocum glider data (JLDBDReader.jl / ERDDAP exports)
include("io/ad2cp_binary.jl")    # native .ad2cp binary reader (later phase)

# ---- Layer 2: per-ping corrections --------------------------------------------------
include("processing/soundspeed.jl")
include("processing/qc.jl")
include("processing/geometry.jl")  # beam↔XYZ↔ENU, 3-beam selection, declination
include("processing/binmap.jl")

# ---- Layer 3: platform kinematics ---------------------------------------------------
include("processing/dac.jl")       # depth-averaged current from nav (DR/GPS), surface drift
include("processing/declination.jl")
include("processing/pipeline.jl")  # ProcessedPings orchestration (Layer 2-3 chain)

# ---- Layer 4: velocity solutions ----------------------------------------------------
include("solutions/shear.jl")
include("solutions/inverse.jl")

# ---- Layer 5: products --------------------------------------------------------------
include("products/grid.jl")
include("products/export.jl")

end # module
