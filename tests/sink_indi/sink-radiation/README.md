* Test name: `sink-radiation`
* Dimension: `3`
* Solver: `mhd`, `INDI_STAR=1`, RTZ + CALIMA (processes off), 8 groups, `static_gas`
* Comparison: analytic, tolerances in `plot-sink-radiation.py` (prints PASSED/FAILED)
* Purpose: photon injection by a main-sequence individual-star sink and its HII region
* Keywords: sink radiation (MIST), RT, Stromgren sphere, reduced speed of light

A 300 Msun main-sequence star (`evolution_flag = 0`, MIST rates, `rt_sink_central_cloud`) shines
into static gas with c_r = 1e-3 c for 0.3 Myr. With the ionising rate Q(t) read from the sink table
in the log (`run.log`, or this test's part of `tests/test_suite.log`) it checks
- the front radius, from the ionised volume, against the case-B Stromgren radius after 100 kyr (10%);
- the front against the finite-light-speed I-front (Shapiro+2006) integrated with Q(t): within
  2 cells at every output, and 15% once c_r t >= R_S (earlier the front is unresolved: a cell's
  light-crossing time exceeds the recombination time, and the explicit M1 step lets a numerical
  precursor run ahead of c_r t);
- the ionising-photon budget: emitted = in flight + ionised + lost to case-B recombinations (5%).

`rt_otsa = .false.`: ground-state recombinations are re-emitted into the groups (so the net loss
is case B), and recombination radiation also adds sub-ionising photons, which is why the budget
uses the ionising groups only.

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
