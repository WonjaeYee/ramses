# Passive-scalar advection accuracy

How much does the hydro step smear, clip, shift or overshoot a passive scalar, by slope limiter, Riemann
solver, Mach number, Courant factor, resolution and direction?  Three scalars ride on a uniform state
(rho = 1, c_s = 1, uniform velocity) through a periodic unit box, as in the square-step and Gaussian tests of
`tests/dust_tva`, but with nothing coupled to them:

| scalar | profile | what it shows |
|---|---|---|
| 1 | top-hat, abs(s - 0.5) < 0.15, zero outside | edge smearing (10-90 % width), clipping, new extrema below 0 and above 1 |
| 2 | Gaussian, sigma = 0.05 | peak loss, numerical diffusivity D_num, centroid error |
| 3 | sine, 0.5 + 0.5 sin(2 pi s) | amplitude and phase error of one Fourier mode; convergence order |

s = kx x + ky y: k = (1, 0) along the grid, (1, 1) the diagonal (2-D); the flow is along k at Mach 0.1, 1 or
10.  The run lasts until the profile has moved by a whole number of periods.  The initial table and the exact
solution are both exact cell averages, computed from the profiles' antiderivatives with no quadrature.  The
exact solution is taken at the snapshot's own time, because RAMSES overshoots `tout` by up to one step.

## Files

| file | |
|---|---|
| `acc_cases.py` | profiles, exact cell averages, the `full` and `nightly` matrices, run directories |
| `acc_build.sh` | builds a plain executable from any tree (`hydro`/`mhd`, 1-D/2-D, 3 scalars, no RT/RTZ/CALIMA); it composes the patch from the tree's own `condinit.f90` plus `table_<solver>.f90` (`condinit_kind='acc_table'`) |
| `table_hydro.f90`, `table_mhd.f90` | the initial-condition reader: `ic_table.dat` (rows in s) and `ic_geom.dat` (k, velocity) |
| `acc_run.py` | runs a matrix in parallel and collects `results.json` |
| `acc_check.py` | the errors of one run against the exact solution (`acc.json`, `profile.npz`) |
| `acc_compare.py` | PASS/FAIL: conservation, the first-order theory, monotonicity with upwind fluxes, regression against a reference `results.json` |
| `acc_report.py` | `summary.md` and the figures |
| `acc_replica.py` | numpy replica of the 1-D MUSCL-Hancock scalar step: where the new extrema of MC + llf/hll come from |

## Running

    B=<dir for the builds>; O=<output dir>
    acc_run.py --matrix full --tree <ramses tree> --builds $B --out $O -j 8   # 904 cases, under a minute on 10 cores
    acc_report.py $O
    acc_compare.py $O/results.json [<reference results.json>]

Set `ACC_VISU` to the tree's `tests/visu` when this directory is not inside that tree.  `--matrix nightly`
has 52 cases that run in a few seconds.

## Metrics (`acc_check.py`)

- **L1, L2, Linf:** taken over the cells.
- **mass:** relative change of the box mass, which a flux-form scheme keeps to round-off.
- **over, under:** new maxima or minima beyond those of the exact solution.
- **TV ratio:** the profile's total variation over the exact one's.
- **top-hat width:** the 10-90 % width of the edges, in cells.
- **Gaussian:** the peak ratio; `D_num`, the growth of sigma² over 2t, in units of v dx; the centroid shift. D_num is valid only while the Gaussian is narrower than ~0.12. For slope_type 0 it is (1 - nu)/2 with an upwind flux, (c_s/v - nu)/2 with hll and (1 + c_s/v - nu)/2 with llf, both subsonic.
- **sine:** the amplitude ratio and displacement of the k = 1 mode. For slope_type 0 in 1-D hydro, also `amp_theory` = |G|^n of the first-order scheme with the solver's own scalar flux. RAMSES matches it to 1e-12, which validates the harness.

## What the scheme does (integration/staging-all-prs, 2026)

- **Scalar flux by Riemann solver:**
  - hllc, exact and acoustic (hydro), and hlld and roe (MHD), upwind the scalar by the contact speed.
  - llf and hll flux it with the acoustic (fast) wave speeds. Their scalar diffusion scales with c_s, not v: at Mach 0.1, N = 128, D_num is about 11x hllc's at first order (theory), 8-9x with minmod and 21-23x with MC.
  - hll is upwind once the flow is supersonic (S_L = 0).
- **MHD `riemann='upwind'`** calls the Lax-Friedrichs flux for the faces (`iriemann=4` in `mhd/umuscl.f90`). For the scalars it is llf.
- **slope_type 5 (ultrabee) and 6 (unlimited)** act on the density only; every other variable, the scalars included, gets a zero slope. For passive scalars they are first order.
- **MC (slope_type 2) with llf or hll at subsonic speed makes new extrema:**
  - The size is 2e-3 at Courant 0.2, 2e-2 at 0.8 and 5e-2 in MHD llf. A zero background goes negative.
  - The Hancock predictor causes it: it moves the face states with v, while llf/hll dissipate at v + c_s (`acc_replica.py` reproduces RAMSES, and the extrema vanish without the predictor).
  - Minmod never makes any. At Mach 0.1, slope_type 8 (theta = 1.5) makes none and van Leer only 1e-8; with llf at Mach 1 and 10 they reach 1e-3 and 5e-3.
  - Upwind fluxes make none (`acc_compare.py` "monotone").
- **slope_type 8 uses theta = 1.5:** `slope_theta` is not in the `hydro_params` namelist.
- **Building:**
  - MHD with passive scalars and NENER=0 did not compile (`irad` in `mhd/godunov_utils.f90` was declared only for NENER > 0). This is fixed on feature/passive-scalar-cma (fcbbe38a); the tests still build MHD with NENER=1, which works on any tree.
  - `slope_mag_type` defaults to `slope_type` and has no 3, so MHD 2-D with slope_type 3 needs `slope_mag_type` set.
