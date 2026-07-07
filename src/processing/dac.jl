# Layer 3 — depth-averaged current (DAC) and surface drift from navigation.
#
# SeaExplorer: DR positions integrate the flight model between GPS fixes; at surfacing
# the first fix minus the final DR position over submerged duration gives DAC.
# Slocum: use m_water_vx/vy (+ m_gps_mag_var rotation) or recompute from GPS/DR.
# Surface drift between fix pairs gives a near-surface velocity constraint.
#
# Implemented in Phase 3.
