* Test name: `stromgren2d_calima`
* Dimension: `2`
* Solver: `hydro`
* Purpose: check that the RT + CALIMA build (`RT=1 RTZ=0 CALIMA=1`) does not fail on a
  real ionisation problem, and that it reduces to plain RAMSES-RT when there is no dust.
* Keywords: RT, CALIMA, dust, Stromgren

This is `tests/rt/stromgren2d` recompiled with CALIMA (`NDUST=1`, `NPAH=0`), with
`dust_tva` and `dust_radpressure` on and every other dust process left at its `.false.`
default. Two namelists, identical apart from the dust:

| namelist | dust | what it establishes |
| --- | --- | --- |
| `stromgren2d_nodust.nml` | `var_region` unset, so the dust slot stays 0 | with no dust the CALIMA terms must vanish, so this must reproduce plain RAMSES-RT |
| `stromgren2d_dust.nml`   | `var_region(1,1)=1d-2`, a 1% dust-to-gas ratio | the grains absorb ionising photons, feel radiation pressure and drift |

Without dust the RT+CALIMA run reproduces plain RAMSES-RT to round-off (part 09 of the stack,
every cell of every output): |dx_HII| <= 9e-15, density and pressure within 7e-14.

`condinit` sets `u = rho*q` for passive scalars, so `var_region(1,1)` is a mass fraction.
Slot 1 is `DustBin_01` (`idust=5`); the outputs name the ion fractions `H_01` and `H_02`
(`output_*/hydro_file_descriptor.txt`), which the plot script reads by name.

`before-test.sh` links the CALIMA tables (`calima_data/`, from `$RAMSES_RTZ_DATA`, default
`<repo>/rtzdata_for_public_ramses`) and creates `SEDtables/`, where
`init_dust_mean_cross_sections` writes its Rosseland/Planck diagnostic dumps.

Note this test has `rt_isIR=.false.`, so it does NOT exercise the IR/dust-temperature path
(`compute_dust_precool`, `Prad_dust`) -- `tests/dust_tva/levgas_rt15_dust_nortz` covers that.
What it does exercise, and the levitation test cannot, is the gas photoionisation opacity in
`dust_radpressure.f90`, because `group_csn` is non-zero here.

The plain-RT reference for the comparison lives in the sibling folder
`tests/rt/stromgren2d_plainrt` (same Makefile with `CALIMA=0`), and
`plot-stromgren-calima.py` here reads both and writes the projection and profile figures.

## Build and run

The build flags are in each folder's `config.txt` (`bin_dust/Makefile`): RT+CALIMA without RTZ,
one grain bin of one chemical type, no PAH bins, no metals, no non-thermal energy
(`NDIM=2 RT=1 RTZ=0 CALIMA=1 NDUST=1 NPAH=0 NDCHEMTYPE=1 NGROUPS=2 NIONS=2 NMETALS=0 NENER=0`);
the plain-RT reference is the same with `CALIMA=0`.

    cd ../stromgren2d_plainrt && mpirun -np 4 <plain-RT ramses2d> stromgren2d_rt.nml
    cd ../stromgren2d_calima  && sh before-test.sh && mpirun -np 4 <RT+CALIMA ramses2d> stromgren2d_dust.nml
    python3 plot-stromgren-calima.py

Each run takes a few seconds on 4 cores.

## Projections and profiles

Snapshots at t = 1.02, 3.06, 7.01, 15.04 Myr (`noutput=4`; `output_00001` is t=0).

| figure | content |
| --- | --- |
| `stromgren2d_maps_rt.png` | plain RT: n_H, x_HII, T |
| `stromgren2d_maps_calima.png` | RT+CALIMA: n_H, x_HII, T, dust eps/eps_0 on a diverging norm centred on the
  initial value (`TwoSlopeNorm(0, 1, 1.15)`), so the total evacuation and the few-per-cent
  shell are both visible. A log or linear scale shows one or the other, not both. |
| `stromgren2d_maps_diff.png` | x_HII(CALIMA) - x_HII(RT) |
| `stromgren2d_profiles.png` | radial profiles: x_HII, x_HI, T, n_H, v_r, dust (full range),
  a dust shell zoom, and Delta x_HII / Delta T residuals |

What the runs show:

| t [Myr] | R_ion RT | R_ion CALIMA | <T>_ion RT | <T>_ion CALIMA | max abs dx_HII | min eps |
| --- | --- | --- | --- | --- | --- | --- |
| 1.02 | 0.2109 | 0.2109 | 17384 | 17412 | 1.0e-04 | 1.3e-04 |
| 3.06 | 0.4453 | 0.4453 | 16158 | 16166 | 2.2e-03 | 1.6e-07 |
| 7.01 | 0.5859 | 0.5859 | 15109 | 15112 | 2.0e-03 | 5.1e-11 |
| 15.04 | 0.7109 | 0.7109 | 15109 | 15111 | 9.4e-04 | 2.0e-15 |

These are the numbers of the stack's part 09 (2026-10-16). The first version of this test
(2026-09-18, with the executables built then) had the fronts at 0.148-0.664. The plain-RT
fronts of the original part 09, without the folded review fixes, are the same as here, so the
difference predates those fixes; it is the same in both runs.

The ionisation front and the mean ionised-gas temperature are unchanged by the dust at this
column: 1% of 0.01 um graphite over 1 kpc at n_H = 0.1 cm^-3 is optically thin to the ionising
photons compared with hydrogen. What the dust does do is drift: radiation pressure evacuates it
from the inner region (eps falls by 13 orders of magnitude by 15 Myr) and piles some of it up
inside the front, leaving the gas perturbed only at the 1e-3 level, localised at the front and
at the centre (Delta T reaches +380 K at r -> 0 at 1 Myr, where the dust is gone).

The plain-RT dashed curves in the profile figure lie underneath the solid CALIMA ones; that is
the result, and the two residual panels show the actual difference.

### The dust shell

The pile-up is only a few per cent above the initial 1%, so it needs its own scale in both
figures. Radially binned peak (64 bins):

| t [Myr] | peak <eps/eps_0> | at r | cell max | ionisation front |
| --- | --- | --- | --- | --- |
| 1.02 | 1.021 | 0.086 | 1.035 | 0.211 |
| 3.06 | 1.016 | 0.477 | 1.041 | 0.445 |
| 7.01 | 1.013 | 0.586 | 1.098 | 0.586 |
| 15.04 | 1.010 | 0.758 | 1.068 | 0.711 |

The binned enhancement is 1-2%, with single cells up to 10% above the initial value; it sits
inside the ionisation front and tracks it outward. `plot-stromgren-calima.py` prints both tables.
