# One-zone tests of `charging_model='WDB06rt'`

Grain charge, photoelectric heating and recombination cooling from the local
RT photon groups (`calima/dust_charging_rtgroups.f90`), checked against the
pyCALIMA reference solver that reads the same tables.

1. Export the tables with pyCALIMA (group edges must match `&RT_GROUPS`):
   `python -m pycalima.models.dust_charge.export_dust_charging_rtgroups`
2. Build (GNU, 3D; NDIM<3 does not link because pm/sink_particle.f90 is
   NDIM==3 only; build rt_parameters.o before dust_photophysics.o):
   `make COMPILER=GNU MPIF90=mpif90 NDIM=3 USE_TURB=0 INDI_STAR=0 NDUST=4 NPAH=2`
3. `python make_runs.py --outdir RUNS --exe-new ... --exe-ref ... --rtgroups-tables ...`
   writes, per source SED and dust bin (one non-zero at a time), a 2^3-cell
   periodic box at fixed n_H and T lit by one 'square' source per group
   (Np = rt_n_source / c_red everywhere), for this branch (`new`, debug dump
   on) and a reference build (`ref`, original tables). Each run gets an empty
   local `SEDtables/`, which RAMSES fills. `dust_tables_dir` is relative
   because CALIMA builds some table names in `character(len=128)`.
4. Run each `run.sh`, then compare with pyCALIMA:
   `python diagnostics/dust_charge/compare_ramses_rtgroups.py RUNS/new/* --tables ... [--original-tables ...]`

Results (2026-09-29, 4 SEDs x 4 bins): RAMSES and the reference solver agree
to <= 2e-14 in Z, Gamma and Lambda (2e-13 in sigma_Z). `--free-T` runs crash
at step 0 in the RTZ ion solver (infinite H creation rate) with and without
this change, so they are not used.
