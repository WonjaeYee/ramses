* Test name: `sink-accretion`
* Dimension: `3`
* Solver: `mhd`, `INDI_STAR=1`, RTZ + CALIMA (processes off), 8 groups
* Comparison: analytic, tolerances in `plot-sink-accretion.py` (prints PASSED/FAILED)
* Purpose: Bondi accretion onto a pre-main-sequence individual-star sink
* Keywords: sink accretion (Bondi), individual stars, PIC

A 0.03 Msun seed (final mass 300 Msun, `evolution_flag = 1`) accretes with `accretion_scheme='bondi'`.
RAMSES's rate is Mdot = 4 pi rho_inf (G M)^2/c^3 with rho_inf = <rho>/alpha(x) (Krumholz+2004),
so in a uniform medium M(t) = M0/(1 - A M0 t). From the sink file of every output it checks that
- `acc_rate` equals the formula on the sink's own `rho_gas`, `cs**2` and `msink` (1e-4);
- the mass gained equals the analytic M(t) - M0 (2%; the sink file prints 11 digits);
- the sink stays pre-main-sequence.

At 1e4 K the seed accretes only ~7e-16 Msun/yr, so this checks the rate and the mass bookkeeping,
not a growth curve. Differences from the full-resolution namelist besides the grid: outputs at
0.05-1 Myr, and `p_region` is the 1e4 K pressure. With `p_region = 1` (T ~ 1e-4 K) the first step,
before `rt_Tconst` applies, accretes ~4e-4 Msun of cold gas, ~5e5 times what the next Myr accretes.

Needs `fix/sink-accretion-ghost-cells`: before it, cloud particles in ghost cells did not accrete,
and the sink gained 0.86 (2 ranks) to 0.75 (4 ranks) of its Bondi rate.

## Setup

Reduced from Wonjae Yee's full-resolution test: one sink at rest at the centre of a 48 pc box of
uniform pure-H gas (n_H = 100 cm^-3, T = 1e4 K fixed by `rt_Tconst`), on a uniform 32^3 grid.
Levels 5-6 instead of 6-7: level 6 has no cells but sets `dx_min`, as level 7 did. The IC sits at
the level-5 cell centre (24.75 pc). `ncontrol = 1`, so the sink table is in the log every step.
A few minutes on 4 cores (`run_test_suite.sh -p 4`).

## Data

The namelist names the external data locally; `before-test.sh` links it and creates
`SEDtables/`, where CALIMA writes its mean-opacity tables:

    export RAMSES_RTZ_DATA=/path/to/rtzdata_for_public_ramses     # default: <repo>/rtzdata_for_public_ramses
    export RAMSES_MIST_SAMPLE=/path/to/260822_mist-sample_v2.5_mcut4.0_with-seds.unf
                                                                    # default: $RAMSES_RTZ_DATA/<that name>

With MPI the sink code needs at least 2 ranks (`synchronize_sink_info` broadcasts from rank 1).
The plot script reads its helpers from `tests/visu/sink_indi_utils.py`.
