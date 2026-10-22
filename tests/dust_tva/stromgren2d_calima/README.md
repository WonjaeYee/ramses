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

`condinit` sets `u = rho*q` for passive scalars, so `var_region(1,1)` is a mass fraction.
Slot 1 is `DustBin_01` (`idust=5`) -- check `output_*/hydro_file_descriptor.txt` if the
variable count ever changes, because the ion scalars shift with it (plain RT names them
`scalar_00,scalar_01`; with CALIMA they become `scalar_01,scalar_02`).

`SEDtables/` must exist before the first run: `init_dust_mean_cross_sections` writes its
Rosseland/Planck diagnostic dumps there and does not create the directory.

Note this test has `rt_isIR=.false.`, so it does NOT exercise the IR/dust-temperature path
(`compute_dust_precool`, `Prad_dust`) -- `tests/dust_tva/levgas_rt15_dust_nortz` covers that.
What it does exercise, and the levitation test cannot, is the gas photoionisation opacity in
`dust_radpressure.f90`, because `group_csn` is non-zero here.

The plain-RT reference for the comparison lives in the sibling folder
`tests/rt/stromgren2d_plainrt` (same Makefile with `CALIMA=0`), and
`plot-stromgren-calima.py` here reads both and writes the projection and profile figures.

## Projections and profiles

    cd ../stromgren2d_plainrt && make && mpirun -np 4 ./ramses_stromrt2d stromgren2d_rt.nml
    cd ../stromgren2d_calima  && make && mpirun -np 4 ./ramses_strom2d   stromgren2d_dust.nml
    python3 plot-stromgren-calima.py

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
| 1.02 | 0.1484 | 0.1484 | 18900 | 18902 | 3.9e-04 | 9.7e-05 |
| 3.06 | 0.3047 | 0.3047 | 17412 | 17415 | 1.9e-03 | 6.0e-08 |
| 7.01 | 0.5078 | 0.5078 | 15879 | 15880 | 1.7e-03 | 7.4e-12 |
| 15.04 | 0.6641 | 0.6641 | 15440 | 15441 | 2.4e-03 | 6.6e-17 |

The ionisation front and the mean ionised-gas temperature are unchanged by the dust at this
column: 1% of 0.01 um graphite over 1 kpc at n_H = 0.1 cm^-3 is optically thin to the ionising
photons compared with hydrogen. What the dust does do is drift: radiation pressure evacuates it
from the inner region (eps falls by 16 orders of magnitude by 15 Myr) and piles it up in a thin
shell at the front, leaving the gas perturbed only at the 1e-3 level, localised at the front and
at the centre (Delta T reaches +200 K at r -> 0 where the dust is gone).

The plain-RT dashed curves in the profile figure lie underneath the solid CALIMA ones; that is
the result, and the two residual panels show the actual difference.

### The dust shell

The pile-up is only a few per cent above the initial 1%, so it needs its own scale in both
figures. Radially binned peak:

| t [Myr] | peak <eps/eps_0> | at r | cell max | ionisation front |
| --- | --- | --- | --- | --- |
| 1.02 | 1.072 | 0.070 | 1.081 | 0.148 |
| 3.06 | 1.051 | 0.270 | 1.062 | 0.305 |
| 7.01 | 1.052 | 0.470 | 1.056 | 0.508 |
| 15.04 | 1.064 | 0.690 | 1.112 | 0.664 |

The enhancement sits just inside the ionisation front and tracks it outward: a broad ~5% plateau
between the evacuated cavity and the front, sharpening into a spike at the front by 15 Myr.
