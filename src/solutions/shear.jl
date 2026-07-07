# Layer 4a — shear method (lADCP tradition: Firing & Hummon; gliderad2cp).
#
#   1. per-ping vertical shear of ENU relative velocities across adjacent cells
#      (glider motion cancels within a ping)
#   2. QC + grid shear into depth bins per profile/segment
#   3. integrate vertically → baroclinic (relative) profile
#   4. reference: depth-average of profile matched to DAC (+ optional BT / surface drift)
#   5. optional shear-bias diagnosis & correction (gliderad2cp process_bias)
#
# Implemented in Phase 4.
