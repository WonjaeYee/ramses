# Passive-scalar conservation tests

Do the passive scalars of RAMSES-RTZ/CALIMA (elements, ion stages, H2, CO, dust and PAH
bins) stay conserved and mutually consistent through the hydro step, AMR, the TVA drift
and the source terms?  The constraint tree the tests check:

    sum(elements) + CO + sum(dust bins) + sum(PAH bins) = rho             "top"
    sum(ion stages of element e) (+ H2 for hydrogen) = rho_e (or rho)    "children"

and the conservation of the box mass of every slot and of every ion stage (x_ion*rho_e,
the quantity the chemistry works with).

## Layout

| dir | what it isolates |
|---|---|
| `advect/` | transport of a composition profile through a periodic box (contact, uniform pressure); face-state sums under independently limited slopes, ion mixing under the storage convention |
| `sod/`    | a Sod tube whose waves sweep three compositions (contacts through a shock and a rarefaction) |
| `amr/`    | the same on an AMR grid refined on the contacts: prolongation of new cells, coarse-fine fluxes, restriction |
| `tva/`    | CALIMA dust/PAH bins drifting relative to the mixture (TVA, imposed drift): the gas counter-flux |
| `replica/`| a Python replica of the 1D MUSCL-Hancock step, sweeping the same matrix in seconds, with the remedies side by side |
| `common/` | case matrix and species sets (`ps_cases.py`), runner (`ps_run.py`), checker (`ps_check.py`), summary table (`ps_table.py`) |
| `patch/`  | `condinit.f90` with `condinit_kind='table'`: initial conditions read from `ic_table.dat` |

## Matrix

Species sets (`ps_cases.SPECIES_SETS`, one RAMSES build each, 1D): `H` (H only),
`HHe_H2` (H, He, H2), `HHeCO` (H, He, C and O with 3 stages each, H2, CO, 4 dust + 2 PAH
bins), `full` (H..Fe with all stages, H2, CO, 4 dust + 2 PAH bins), `full_np0` (H..Fe, 4 dust
bins, no PAH, no H2/CO).  Profiles: `step2` (two parcels: the control, independent limiting
is exact for any mixture of two compositions), `step3` (three parcels, the middle one 2
cells wide, optionally with a density contrast of 100), `fronts` (two smooth fronts of width
0.5 or 4 cells, 3 cells apart), `smooth` (periodic mixture, wavenumber 1 or 8), `spike` (a
1-cell and a 3-cell spike on a smooth background), `sod3`.  Concentrations: `solar` and
`extreme` (Z from 1e-4 to 3 solar, dust-to-gas 1e-6 to 0.1, trace stages at 1e-10 of their
element next to dominant ones).  Slope limiters minmod and MC.

## Running

    run_all.sh <tree> <label> <outdir> [--quick] [--nobuild] [--tests "advect sod amr tva"]

builds the five species sets from `<tree>` (`build.sh`), runs every case and prints a
PASS/FAIL line per case plus a table of the worst errors per (test, set, profile) to
`<outdir>/summary_<label>.txt`.  One test: `common/ps_run.py --test advect --exes ...`; one
run directory: `common/ps_check.py <rundir>` (writes `check.json`, `check.png`).
Data: `$RTZ_DATA` (default `RAMSES_dev/rtzdata_for_public_ramses/`).

## Tolerances

Round-off where the scheme guarantees the property (`ps_check.ROUNDOFF`): box mass of every
slot (any flux-form scheme) 1e-12; with the consistent scheme also `top` 1e-12, `child_abs`
(error of the ion sum in units of the element's peak density) 1e-12 and the ion masses
1e-12.  Properties the scheme does not guarantee are held to the bound measured with the
legacy scheme (`bounds_legacy.json`, x2) or reported as INFO.  The scheme is read from the
run log: `legacy` (ion slot = x*rho), `cma` (ion slot = x*rho_e with passive_cma), `legacyB`
(the latter with passive_cma=.false.).
