# Layer 1 — Slocum glider data ingestion.
#
# Sources, in order of preference:
#   1. JLDBDReader.jl (user's pure-Julia dbd/ebd reader) — m_water_vx/vy, m_gps_lat/lon,
#      m_pitch/m_roll/m_heading, m_gps_mag_var, sci_water_pressure, ...
#   2. Pre-exported CSV / ERDDAP tabledap files (Slocum-AD2CP style).
#
# Implemented in Phase 5 (SeaExplorer first; the sea064 M38 reference dataset drives Phases 1–4).
