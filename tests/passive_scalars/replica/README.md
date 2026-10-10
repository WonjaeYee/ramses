# replica

`ps1d.py` replicates the 1D RAMSES MUSCL-Hancock step (ctoprim, uslope slope_type 1/2/7/8,
trace1d, HLLC and LLF, conservative update) for the passive-scalar layouts of
`common/ps_cases.py`, with both storage conventions (A: ion = x*rho, B: ion = x*rho_e) and
the remedies: `legacy`, `cma` (nested renormalisation of the face states), `cma_y`
(reconstruct the ion fractions of the element), `renorm` (post-step renormalisation), and a
TVA-like drift with or without the gas counter-flux.  `run_replica.py` sweeps the matrix and
prints PASS/FAIL against the properties each scheme guarantees (~2 min).
