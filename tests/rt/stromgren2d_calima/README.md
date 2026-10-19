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
every cell of every output): |dx_HII| <= 1.4e-13, density within 2.2e-14, temperature within 7.5e-13.

`condinit` sets `u = rho*q` for passive scalars, so `var_region(1,1)` is a mass fraction.
Slot 1 is `DustBin_01` (`idust=5`); the outputs name the ion fractions `H_01` and `H_02`
(`output_*/hydro_file_descriptor.txt`), which the plot script reads by name.

`before-test.sh` links the CALIMA tables (`calima_data/`, from `$RAMSES_RTZ_DATA`, default
`<repo>/rtzdata_for_public_ramses`) and creates `SEDtables/`, where
`init_dust_mean_cross_sections` writes its Rosseland/Planck diagnostic dumps.

Note this test has `rt_isIR=.false.`, so it does NOT exercise the IR/dust-temperature path
(`compute_dust_precool`, `Prad_dust`) -- `tests/dust_tva/levgas_rt15_dust_nortz` covers that.
What it does exercise, and the levitation test cannot, is the gas photoionisation opacity in
`dust_radpressure.f90`, because `group_csn` is non-zero here, and the grain absorption of the
UV groups in the RT: without an IR group the grains take their rates from the same cross
sections as their radiation pressure (`calima_optics` in `rt_cooling_module`).

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
  initial value (`TwoSlopeNorm(0, 1, 1.15)`), so the total evacuation and the 8-15%
  shell are both visible. A log or linear scale shows one or the other, not both. |
| `stromgren2d_maps_diff.png` | x_HII(CALIMA) - x_HII(RT) |
| `stromgren2d_profiles.png` | radial profiles: x_HII, x_HI, T, n_H, v_r, dust (full range),
  a dust shell zoom, and Delta x_HII / Delta T residuals |

What the runs show:

| t [Myr] | R_ion RT | R_ion CALIMA | <T>_ion RT | <T>_ion CALIMA | min eps |
| --- | --- | --- | --- | --- | --- |
| 1.02 | 0.1484 | 0.1328 | 18900 | 19130 | 1.2e-04 |
| 3.06 | 0.3047 | 0.2734 | 17412 | 17818 | 7.4e-08 |
| 7.01 | 0.5078 | 0.4141 | 15879 | 16716 | 4.9e-12 |
| 15.04 | 0.6641 | 0.5859 | 15440 | 16220 | 7.1e-18 |

These are the numbers of the stack's part 09 after the review of #23 (2026-10-19). Two changes
of that review moved them:
- The gas starts in chemical equilibrium, i.e. molecular at 100 K without UV. An earlier special case
  in `rt_init_xion` started it atomic, and both fronts were then further out (0.211-0.711).
- The grains absorb the UV in the RT. Before, without an IR group the RT used `kappaAbs`, which is 0
  here, so the dust felt radiation pressure but absorbed nothing, and the two runs had the same front.

The dust takes a large share of the ionising photons inside the H II region. With 1% of 0.01 um
graphite, its cross section per H (~6e-21 cm^2) exceeds that of the residual H I there. Across the
H II region tau_dust ~ 1.
- **The front:** 10-18% smaller in radius than in plain RT.
- **The temperature:** the ionised gas is 230-840 K hotter, since the dust removes the softer
  ionising photons and hardens the spectrum.
- **The drift:** radiation pressure evacuates the dust from the inner region (eps falls by 15
  orders of magnitude by 15 Myr), and piles some of it up behind the front.

### The dust shell

The pile-up is 8-15% above the initial 1%, so it needs its own scale in both figures. Radially
binned peak (64 bins):

| t [Myr] | peak <eps/eps_0> | at r | cell max | ionisation front |
| --- | --- | --- | --- | --- |
| 1.02 | 1.083 | 0.070 | 1.098 | 0.133 |
| 3.06 | 1.087 | 0.117 | 1.091 | 0.273 |
| 7.01 | 1.118 | 0.180 | 1.123 | 0.414 |
| 15.05 | 1.150 | 0.273 | 1.152 | 0.586 |

The shell sits well inside the ionisation front, at about half its radius, and grows as it moves
out. `plot-stromgren-calima.py` prints both tables. Once the dust changes the step sequence, the
CALIMA snapshots are at slightly different times from the plain-RT ones (0.07% at 15 Myr), since
both runs write at the first step past each `tout`.
