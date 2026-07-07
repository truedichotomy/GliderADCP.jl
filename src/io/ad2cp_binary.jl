# Layer 1 (native) — .ad2cp binary reader, removing the MIDAS/Windows dependency.
#
# Format (Nortek "Integrators Guide — AD2CP", §Data formats; cross-checked vs DOLfYN):
#   10-byte header: 0xA5, hdrSize=10, ID, family=0x10, dataSize u16, dataChecksum u16, hdrChecksum u16
#   IDs: 0x15 burst, 0x16 average, 0x17 bottom track, 0x18 interleaved burst, 0x1A/0x1F alt raw,
#        0x1C echosounder, 0xA0 string record (string ID 0x10 = full instrument config dump)
#   checksum: u16 wraparound sum init 0xB58C (+ last odd byte << 8)
#   DF3 records self-describing: config bitmask (which blocks present), beams_cy word
#   (nCells bits0-9, coordsys bits10-11, nBeams bits12-15), offsetOfData, velocityScaling 10^k.
#
# Implemented in Phase 6.
