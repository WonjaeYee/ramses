# sod

Sod shock tube (rho 1 | 0.125, P 1 | 0.1) on 256 cells, outflow boundaries, t=0.2, with
compositions A | B | C | A: the rarefaction runs into the A|B contact, the gas contact
carries the B|C one and the shock overtakes the C|A one.  Isolates the face-state sums under
compression and expansion, with the hydro limiter acting on steep density and velocity
profiles.  PASS for the consistent scheme: sums to round-off (masses are not checked: the
boundaries are open, though nothing reaches them by t=0.2).
