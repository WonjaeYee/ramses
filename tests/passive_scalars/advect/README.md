# advect

A composition profile (`common/ps_cases.profile`) carried through a periodic box at u=1,
uniform pressure, for one crossing (t=1).  Cooling and chemistry off (`rtz_cooling=.false.`;
`cooling_fine` still runs through `neq_chem`).  Isolates the hydro step:

* face states: every scalar is limited on its own, so where three or more compositions meet
  in a stencil the face fractions no longer sum to one (to their element for the ions) and
  the species fluxes no longer sum to the mass flux (`top`, `child`, `child_abs`);
* the storage convention: with ion slot = x*rho, mixing gas of different metallicity mixes
  the ionisation fractions with the total mass instead of the element mass, so the ion
  masses x*rho_e are not conserved (`ionmass`).  `step2` (two parcels) shows this alone.

PASS for the consistent scheme: all sums and masses to round-off.  The legacy scheme
conserves the slot masses only.
