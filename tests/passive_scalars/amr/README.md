# amr

As `advect`, on levels 6..8 refined on the density jump and on the gradients of the H mass
fraction and of the HI stage (`err_grad_var`), with `interpol_type=2` (MC).  The refined
patches follow the contacts, so new fine cells are created by prolongation every few steps,
fluxes cross coarse-fine boundaries and fine cells are restricted.  Isolates
`interpol_hydro` (every variable limited on its own: each species conserved, their sum not
equal to the interpolated rho) and the reflux.  PASS for the consistent scheme: sums and
masses to round-off.
