* Test name: `stromgren2d_rtz_calima`
* Dimension: `2`
* Solver: `hydro` (static gas)
* Purpose: a cheap Stromgren test of the RTZ + CALIMA grain charging, photoelectric
  heating and recombination cooling with the three charging models: the original
  (log gamma, log T) tables (`charging_model='RM2026'`), the per-group charge balance on the
  local RT photon groups (`'WDB06rt'`) and its tables on (T, FUV flux, FUV hardness, EUV
  flux) (`'WDB06tab'`), both in `calima/dust_charging_rtgroups.f90`.
* Keywords: RTZ, CALIMA, dust charging, photoelectric heating, Stromgren

Same chemistry, dust (4 bins) and PAHs (2 bins) as the 3D
`RAMSES_dev/stromgren/calima_cooling` setup, but 2D, 64^2 cells, static gas, free
temperature, and one point source in the (mirrored) corner emitting a 4e4 K blackbody into
the 8 groups, Q(>13.6 eV) = 2.5e29 in RAMSES 2D source units. n_H = 0.1 cm^-3, 100 K
neutral gas initially; snapshots at 1, 2, 3, 5 Myr. About 7-8 min per run on 4 ranks.
`tests/rt/stromgren2d_calima` does not test the charging: it is RT=1 RTZ=0, and the
charge is only computed in the RTZ cooling step (`compute_dust_precool`).

## Build and run

    make COMPILER=GNU MPIF90=mpif90 NDIM=2 USE_TURB=0 INDI_STAR=0 NDUST=4 NPAH=2 NGROUPS=8
    # (on a clean tree: make ... mpi_mod.o amr_parameters.o hydro_parameters.o rt_parameters.o first)

One run directory per model, named after its `charging_model`, each with
`stromgren2d_rtz.nml` (set `charging_model`), an empty `SEDtables/`, and `tables/` with the
CALIMA tables plus, from pyCALIMA:

* `WDB06rt`: `python -m pycalima.models.dust_charge.export_dust_charging_rtgroups` (the
  default 8 groups);
* `WDB06tab`: `python -m pycalima.models.dust_charge.export_dust_charging_psitables --egy
  0.6697,3.741,8.541,12.41,14.40,19.70,35.09,65.71 --workers 8` (built from the RT-group
  tables at the run's group mean energies, the `egy` of `info_rt`; RAMSES warns if they
  differ).

`data_dir` in `&HYDRO_PARAMS` must point to the RTZ data.
`dust_charging_timer=.true.` prints, per rank, the RTZ cell sub-steps and the time of the
steps and of the grain charging.

    mpirun -np 4 ./ramses2d stromgren2d_rtz.nml

From the directory holding the run directories:

    python plot-stromgren2d-rtz.py --pycalima PYCALIMA --rtgroups-tables DIR --calima-tables DIR \
        --psitab-tables DIR [--runs RM2026 WDB06rt WDB06tab]

It reads the outputs directly (yt cannot parse the RTZ `info_rt` file nor 2D runs with
boundary regions), plots x_HII, T, n_e, N_g profiles and, at 2 Myr, the grain potential of
every bin and the PE heating and recombination cooling per H, each as the run's own model
gives them on its cell state (pyCALIMA `solve_cell`, `psi_lookup`, or the original tables at
G0 of the groups with 5.6 < egy < 13.6 eV).

To check RAMSES against pyCALIMA, restart a run with `nrestart=3`,
`dust_rtgroups_verify=.true.` and `dust_rtgroups_debug_max=40000` (it dumps every call and
stops), then

    python PYCALIMA/diagnostics/dust_charge/compare_ramses_rtgroups.py RUN --tables DIR      # WDB06rt
    python PYCALIMA/diagnostics/dust_charge/compare_ramses_rtgroups.py RUN --tables x --psitab DIR  # WDB06tab

and, for the errors of both against the full charge distribution,
`PYCALIMA/diagnostics/dust_charge/check_charging_methods_errors.py`.

## Results (2026-10-01)

* RAMSES vs pyCALIMA: WDB06rt, 74,116 full solves from restarts at 1 and 2 Myr
  (discrete and wide P(Z), warm-started, recombination fits and reuse), max rel. diff
  2e-13 (Z, Gamma, Lambda), 1e-9 (alpha), no branch mismatches; WDB06tab, 80,000
  lookups, 2e-13.
* Against the full P(Z) on 3,000 cells (cell totals, Milky Way abundances): WDB06rt PE
  heating 99% 3%, recombination cooling 1.4%, grain recombination <= 13%; WDB06tab PE
  heating median 0.6% (max 11%), <Z> 1.5-3%, grain recombination median 6-7% (tails from
  n(H+) = n_e and the blackbody EUV, mostly in the HII gas).
* Cost of the charging per RTZ cell sub-step (0-2 Myr, 1 rank): RM2026 1.0 us (4.5% of
  the step), WDB06tab 1.35 us (6.8%), WDB06rt 5.1 us (20%). Run wall times also differ
  because RM2026 needs ~20% more RTZ steps, and on 4 ranks the rank at the source does
  ~2/3 of them.
* Ionisation front: the same in the three runs at 1 and 2 Myr (355, 565 pc); at 5 Myr
  WDB06rt and WDB06tab are ~20 pc ahead of RM2026 (955 vs 935 pc). HII gas: WDB06tab within
  0.1-1% of WDB06rt, both 4-8% warmer than RM2026; neutral gas 78-100 K vs 40-67 K
  (RM2026) up to 3 Myr.
* Grain potential at 2 Myr: WDB06rt and WDB06tab agree (HII gas +1 to +6 V near the source,
  about 0 V at the front, +0.2 to +0.4 V in the FUV-lit neutral gas); the original tables
  have no EUV photoemission and reach -37 V (the edge of their table) in the HII gas.
