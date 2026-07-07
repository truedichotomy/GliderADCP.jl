# Layer 4b — least-squares inverse method (Visbeck 2002; Todd et al. 2017; Gradone 2023).
#
# Unknowns m = [u_glider(t_1..t_nt); u_ocean(z_1..z_nz)] per segment (complex u+iv).
# Rows:
#   measurement:  u_ocean(z(cell)) − u_glider(t_ping) = v_rel   (+1/−1 sparse)
#   DAC:          weighted depth-integral of u_ocean (or time-integral of u_glider) = DAC
#   bottom track: u_glider(t) = −v_BT(t)                        (strong weight, when valid)
#   smoothness:   second-difference regularization on both blocks
# Solve sparse LSQR (Krylov.jl). Weights per Todd et al. 2017 (wDAC=5, wSmooth=1 defaults).
# Fixes documented Slocum-AD2CP v2.0.0 bugs (see PLAN.md §Inverse); offers strict
# parity mode for regression against reference outputs.
#
# Implemented in Phase 4.
