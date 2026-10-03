* Test name: `sink-accretion`
* Dimension: `3`
* Solver: `mhd`, `INDI_STAR=1`, RTZ + CALIMA (processes off), 8 groups
* Comparison: analytic, tolerances in `plot-sink-accretion.py` (prints PASSED/FAILED)
* Purpose: Bondi growth of a pre-main-sequence individual-star sink and its switch to the main sequence
* Keywords: sink accretion (Bondi), individual stars, MIST, PIC

A 2 Msun seed (final mass 20 Msun, `evolution_flag = 1`) accretes with `accretion_scheme='bondi'`
from uniform gas at n_H = 1.8e6 cm^-3 and T = 1e3 K (`rt_Tconst`). RAMSES's rate is
Mdot = 4 pi rho_inf (G M)^2/c^3 with rho_inf = <rho>/alpha(x) (Krumholz+2004), so the mass runs
away as t(M) = int dM/Mdot, close to M0/(1 - A M0 t), and reaches 20 Msun after ~1 Myr. Checks:
- `acc_rate` equals the formula on the sink's own `rho_gas`, `cs**2` and `msink` (1e-4);
- the time at which the logged mass reaches M equals the analytic t(M) (2% of t(M_F));
- the sink reaches the main sequence at t(M_F) (2%), with msink within 1% above the final mass;
- on the main sequence it stops accreting, and its MIST ionising photons start within 10 kyr
  (they are zero below the first age node, then ramp up to ~1.4e48 s^-1 over ~60 kyr).
The per-step growth and photon rates come from the sink table in the log (`run.log`, or this
test's part of `tests/test_suite.log`).

Gravity is off (`poisson = .false.`): the Bondi rate is a subgrid recipe in which gravity enters
only through r_B = G M/c^2. That is consistent while r_B << r_acc = ir_cloud dx_min, where the cells
see rho_inf and alpha ~ 1; here r_B < r_acc/70 up to 20 Msun (alpha - 1 < 1.4%). With r_B ~ r_acc the
alpha correction would assume a Bondi overdensity that only resolved gravity creates.

Reduced and changed from Wonjae Yee's test: levels 6-7 -> 5-6 (uniform 32^3 of 1.5 pc; level 6 has
no cells but sets dx_min), the IC at the level-5 cell centre (24.75 pc), `ncontrol = 1`, outputs every
0.1 Myr to 1.2 Myr, and, so that the sink grows, n_H 100 -> 1.8e6 cm^-3, T 1e4 -> 1e3 K, seed
0.03 -> 2 Msun, final mass 300 -> 20 Msun, `isH2_rtz`/`isCO_rtz` off (H2 would change the sound
speed). At the original n_H = 100 cm^-3 and 1e4 K a 0.03 Msun seed accretes ~7e-16 Msun/yr. `p_region`
is the 1e3 K pressure: with `p_region = 1` (T ~ 1e-4 K) the first step, before `rt_Tconst` applies,
would accrete cold gas at a huge Bondi rate.

Needs `fix/sink-accretion-ghost-cells`: before it, cloud particles in ghost cells did not accrete,
and the sink gained 0.86 (2 ranks) to 0.75 (4 ranks) of its Bondi rate.

## Data

The namelist names the external data locally; `before-test.sh` links it and creates
`SEDtables/`, where CALIMA writes its mean-opacity tables:

    export RAMSES_RTZ_DATA=/path/to/rtzdata_for_public_ramses     # default: <repo>/rtzdata_for_public_ramses
    export RAMSES_MIST_SAMPLE=/path/to/260822_mist-sample_v2.5_mcut4.0_with-seds.unf
                                                                    # default: $RAMSES_RTZ_DATA/<that name>

About 3 minutes on 4 cores (`run_test_suite.sh -p 4`). With MPI the sink code needs at least
2 ranks (`synchronize_sink_info` broadcasts from rank 1). The plot script reads its helpers from
`tests/visu/sink_indi_utils.py`.
