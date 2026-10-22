# restart

`run_restart.py`: the `advect` case step3/extreme/full run straight to t=1 and, separately,
stopped at t=0.5 and restarted from that output (nrestart=6).  The outputs hold the ion
stages as fractions of their element, which the restart multiplies by rho; with ion slot =
x*rho_element `init_passive_groups` converts them to ion mass densities.  PASS if every final density agrees to 1e-10 of its peak (bitwise
with ion slot = x*rho; with x*rho_element the restart rounds x*rho*rho_e/rho, and that 1e-16
seed grows through the limiters to ~3e-12 by t=1).
