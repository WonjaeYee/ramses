# MPI-3 for the RTZ/CALIMA cooling step: a node-level work queue and node-shared tables

Branch `feature/mpi3-cooling` (staging `ac03d014` + PR #73, `cost_measured`). This document is the design, and the
record of what each phase measured. Status of each phase at the end.

## 1. The problem

### Load imbalance in the cooling step
- **Turbulent TVA + H II box** (lv6, `levelmin = levelmax`, 128 ranks on megatron): 86% of the wall is MPI wait.
  - Cells at the ionisation front (0.01 < xHII < 0.99) cost 25-175x a neutral cell in `rtz_solve_cooling`.
  - 1.5% of the grids hold 48% of the cooling time; 12 ranks hold 49% of it.
- **PR #73 (`cost_measured`):** weighting the load balance by the measured cost per cell gives 1,333 s → 304 s on
  9.0-9.3 Myr. Even then, cooling per rank is 1.4x max/mean and 40% of the wall is still MPI wait.
- **Stromgren 2-D:** per-rank cooling totals differ by only 7%, but 37% of the wall is lost at np = 4. The reasons:
  - each RT subcycle is a synchronisation point: `rt_step` exchanges `rtuold` after every `cooling_fine`;
  - the slowest rank changes from subcycle to subcycle.
- A static decomposition, however well weighted, cannot follow that. Only redistributing cells within one
  `cooling_fine` call does.

### Memory
Every rank reads (formatted ASCII) and holds its own copy of the read-only tables:

| tables | MB per rank |
|---|---|
| RTZ + CALIMA | ~70 |
| + `SED_table` | ~95 |
| + the MIST sample (`INDI_STAR`) | ~195 |

At 128 ranks per node that is up to 25 GB per node, and ~13.6 GB of ASCII parsed per node at start-up.

## 2. What makes cells movable

`NVECTOR = 1` in every RTZ build, so each `rtz_solve_cooling` call solves one leaf cell. Inside the solver nothing
reads the grid (`uold`, `rtuold`, `son`, `xg`, neighbours). Its inputs are:
- **values gathered by `coolfine1`:** T2, ion fractions, element densities, nCO, Np, Fp, p_gas, dx_SS_H2, sigma,
  rho_dust, rho_pah, and `tva_wdrift_batch`;
- **globals that are equal on all ranks within one call:** `ilevel`, `dtcool`, `aexp`, the UVB, `rt_c(ilevel)` and
  the namelist.

State carried from one cell to the next, which must be reset or carried with the cell:
- `last_h_ode` (`calima/dust_interface.f90`): the first ODE step guess, kept from the previous cell;
- `rtg_Zguess` (same file): the grain-charge root guess;
- the `save`d `xion(27,27,nvector)` in `coolfine1`: untracked slots keep the values of the previous cell;
- `pah_fion_vec`: written by the solve, read again by the IR-trapping split in `coolfine1`;
- `tva_wdrift_batch`: filled by `coolfine1` from the global `tva_wdrift`; the solve tests `allocated()`.

The other caches (`rec_T`, `ctr_T`, `cmb*_z`, the gas-invariant memo in `grain_relative_velocity`, ...) are exact
memoisation and do not change results.

Every rank calls `cooling_fine(ilevel)` the same number of times and in the same order, including ranks with no
grids on the level:
- the number of RT subcycles depends only on `rt_c(ilevel)`, `dx` and the global `dtnew`;
- the early return in `cooling_fine` tests the global `numbtot`.

So call k is the same call on every rank. The queue relies on this, and `cool_stats` checks it at run time.

## 3. Design

### Phase 0: measurement
- `make_virtual_fine_dp` switched the timer to `'test mpi time N'` and ended on `'temp_else'`. Everything its caller
  did afterwards was booked to `temp_else`, including `rt_godunov_fine` and `add_rt_sources` of the next RT substep.
  It now restores the caller's timer.
- `cool_stats = .true.` (`&RUN_PARAMS`) records the wall time of each `cooling_fine` call and of the `rtuold` halo
  exchange after it in `rt_step`, which is where the imbalance is paid.
- Once per coarse step these are reduced over the ranks and printed as `Cooling balance step ...`:
  - Σ over calls of the max over ranks, against Σ over calls of the mean (the per-call imbalance);
  - the mean and max wait;
  - the step wall;
  - running totals.

### Phase 1: `coolfine1` in three phases, and the record solve (`cool_queue = 1`)
- `coolfine1` gets a phase argument:
  - `CF_LEGACY` (`cool_queue = 0`) is today's code, unchanged and bit-identical;
  - `CF_PRE` runs the gather and the changes made on the cell's own rank, then packs the solver inputs into a record;
  - `CF_POST` unpacks the solver outputs and runs the existing scatter, error block and IR-trapping split.
- No arithmetic statement moves, which matters on aarch64, where gfortran contracts multiply-adds.
- **The record:**
  - in: T2, the tracked ion slots and elements, nCO, Np, Fp, p_gas, dNpdt/dFpdt (`rt_smooth`), dx_SS_H2, sigma,
    rho_dust, rho_pah, `tva_wdrift` and the predicted cost;
  - out: the solved fields, `pah_fion`, the error status, the measured cost and the solving rank.
- A home buffer stays with the owner: ekk, err, emag, the pre-solve T2, T2min, nH, Zsolar and the `rt_vc` boosts.
- `cq_solve_record` makes each record's result a pure function of the record and the call's globals:
  - it rebuilds the solver arrays with a fixed fill;
  - it sets the warm-start state for this cell: the ODE step starts cold, and the charge roots start from the cell's
    own stored values (`cool_queue_zguess_ivar`, section 5) or cold;
  - it sets `pah_fion_vec` and `tva_wdrift_batch` from the record.
- **Sigma** (`cmp_sigma_turb`, a neighbour stencil) is computed for all cells before any of them is cooled, so it does
  not depend on the order.
- **Cost:** the measured solve time of each cell is written to #73's cost field (`unoadvect(:,cost_ivar)`). The queue
  sorts the cells by it, and the load balance can use it.

### Phase 2: the node-level queue (`cool_queue = 2`)
- **Communicator:** `node_comm` from `MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)`.
- **Window:** each rank owns a segment of one `MPI_Win_allocate_shared` window, with a header, a chunk table and the
  in- and out-records. The window stays in a `lock_all` epoch for the whole run (unified memory model required).
- **One call:**
  1. Each rank sorts its cells by predicted cost, most expensive first (LPT), and cuts the list into chunks of about
     `cool_queue_chunk` seconds (an expensive cell is a chunk of its own).
  2. It packs the records into its segment, resets its counters, and enters a node barrier.
  3. It claims chunks from its own queue with `MPI_Fetch_and_op` on `next`.
  4. When its own queue is empty, it claims from the node rank with the largest predicted cost left. The solver
     reads the in-record and writes the out-record straight into the owner's segment, with no copy.
  5. When no queue on the node has work left: a second node barrier, then each rank scatters its own cells.
- **Why barriers:** the two node barriers per call cost a few µs and avoid an epoch-tagged CAS protocol. If they
  ever show up in `cool_stats`, the fallback is a single word per queue holding epoch|nchunk|next, updated by CAS.
- **Determinism:** a record's result does not depend on who solves it. So modes 1 and 2 give the same results, and
  the same results from run to run.

### Phase 3: node-shared tables
- Node rank 0 reads each large table into an `MPI_Win_allocate_shared` segment; the other ranks of the node map it.
- Header scalars and axes are broadcast on `node_comm`.
- Conversion order: `high_t_cooling_rates(_opt)`, CALIMA `collisional_tab`/`sputtering_tab`, `isrf_q`, the read-only
  parts of `rtg`, `SED_table`, and the MIST `sample`.
- The `KgT/HgT/mKT/mHT/uni` caches of `rtg` change at run time and stay per rank.

### Phase 4: other communication changes
- Exchange several variables in one halo call: `rt_step` makes 32 + 32 separate `WAITALL` rounds per RT substep.
- Process binding.
- Hoist `rtz_updateRTGroups_CoolConstants` to once per call.

### Later (designed, not built)
- **Cross-node stealing (Phase B):**
  - the claim/done counters move to a window over `MPI_COMM_WORLD`;
  - the records are read and written remotely (`MPI_Get`/`MPI_Put`) through a second window over the same memory,
    with larger remote chunks;
  - a rank steals from another node only when its own node has no work left.
  - Measure first: Σ_k (max over nodes − mean over nodes) from `cool_stats`.
- **Overlap** of the `rtuold` exchange with the cooling of interior cells (boundary cells first).
- **Persistent requests**, a **zero-copy halo exchange** within a node, **neighbourhood collectives** in
  `build_comm`, and **`MPI_Iallreduce`** in `output_rt_stats`.

## 4. Build notes
- The new modules sit in existing files, so no Makefile changes are needed, including the test suites' own
  Makefiles:
  - `cool_queue` at the top of `hydro/cooling_fine.f90`;
  - `node_shared` in `amr/mpi_mod.f90`.
- `amr_step.f90` and `adaptive_loop.f90` are compiled before `cooling_fine.f90`, so they call external routines
  (`cool_stats_call`, `cool_stats_report`).

## 5. Results by phase

### Phase 0 (timer, `cool_stats`)
- Tier 1 (M4 Max): 8/8 PASS, every output byte-identical to staging `ac03d014`
  (`RAMSES_dev/tests_final/runs/mpi3_phase0`).

### Phase 1 (records, `cool_queue = 1`)
- **Tier 1 with `cool_queue = 0`:** 7/7 PASS, byte-identical to staging (`runs/mpi3_phase1`). The RTZ-only ladder
  build failed there on `ndust` and was fixed; the Phase 2 run below covers it.
- **t1_pdr_calima** (1-D H II region + PDR, 64 cells, 8 groups, full CALIMA):
  - records solved in reverse order: bit-identical to natural order, so no state carries over from one cell to the
    next;
  - `cool_queue_debug = 2` (in place, with the per-cell warm-start reset) is bit-identical to `cool_queue = 1`. The
    packing and unpacking are exact, and the reset is the only change;
  - `cool_queue = 1` against `0`: the largest relative change of any field is 7.5e-10 (pressure, then Si/Ne ion
    fractions);
  - `cool_queue = 1` at np 1 and np 4: equal.

### Phase 2 (node queue, `cool_queue = 2`)
- **Tier 1 with `cool_queue = 0`:** 8/8 PASS, byte-identical to staging (`runs/mpi3_phase2`).
- **t1_pdr_calima, np = 4 (M4 Max):** two `cool_queue = 2` runs are bit-identical to each other and to
  `cool_queue = 1`.

  | | wall | cooling Σmax / Σmean over calls | wait after cooling (mean) |
  |---|---|---|---|
  | `cool_queue = 0` | 2.28 s | 2.23 / 0.89 s | 1.35 s |
  | `cool_queue = 1` | 2.20 s | 2.15 / 0.91 s | 1.26 s |
  | `cool_queue = 2` | **1.11 s** | **1.064 / 1.063 s** | 0.013 s |

  The front cells sit on one or two of the four ranks. With the queue, the other ranks solve them: 60 chunks taken
  per coarse step, and 1.4 ms waiting at the closing barrier.

### The turbulent TVA + H II box on megatron, lv5 (32³) from t = 0 to 0.41 Myr, 32 ranks
Set-up: `/data/currodri/turbhii/mpi3_20261022/lv5_*`. It is the lv6 box's namelist at `levelmin = levelmax = 5` with
`nrestart = 0`. Every run used nice 19 + ionice on a shared machine (load 200-230 of 256), so wall times vary by
about ±8%.

| run | `cool_queue` | `cost_measured` | wall | cooling Σmax / Σmean | wait after cooling (mean) |
|---|---|---|---|---|---|
| E | 0 | 0 | 280 s | 273 / 112 s | 162 s (58% of the wall) |
| A | 0 | 1 (#73) | 171 s | 159 / 112 s | 47.5 s |
| F | 1 | 0 | 266-296 s | 259-288 / 116-126 s | 143-163 s |
| **C** | **2** | **0** | **135 s** | **128.2 / 128.0 s** | **0.64 s** |
| B | 2 | 1 | 141 s | 129.4 / 129.1 s | 1.0 s |
| D | 2 + `shared_tables` | 1 | 142 s | 130.7 / 130.3 s | 1.1 s |

- **Speed-up:** 2.07x over no balancing, and 1.26x over #73 alone.
- **#73 is not needed on one node:** with the queue, the measured-cost balance (B) gains nothing over the plain
  decomposition (C), because the queue already shares the work out.
- **Work per rank:** the queue solves the same records in the same time as mode 1 (Σ solve time per rank 125.6 s
  against 124.5 s, back to back at the same load). Packing, claims and barriers cost 2.2 s (1.7%) of the 127.8 s of
  cooling.
- **Stealing:** about half of a rank's solving at step 29 is other ranks' cells, about 37,000 chunks per coarse step.
- **Mode 1 against mode 0** (records plus the warm-start reset): within the load noise, 0-4%.
- **Order of the records:**
  - most expensive first is 6% faster in wall than grid order (C 135 s against 143 s with chunks in grid order);
  - in mode 1, cost order and grid order cost the same (126.5 s against 125.8 s), so cache locality does not matter.
- **Bit-identity:** C (mode 2) and F (mode 1) are bit-identical in all 81 fields at t = 0.41 Myr, 32 ranks, 3-D,
  with TVA, coagulation and shattering (and so `cmp_sigma_turb`).
- **Mode 0 against mode 2:** largest change 7e-10 (dust bins), from the warm-start reset.
- **Decomposition dependence (from #73's balancing, not the queue):** B against C changes the dust by up to 3.8e-5.
  With `cost_measured = 1` the decomposition follows wall times, and `cmp_sigma_turb` reads stale ghost-cell momenta
  inside the RT subcycles (see OPEN_ISSUES).
- **Memory:** the summed PSS of the 32 ranks is 6.71-6.79 GB with `shared_tables` (D), against 8.19-8.28 GB without
  (B). That is -1.49 GB, or 47 MB per rank, matching the ~50 MB of converted tables.

### The turbulent box at lv6 on megatron, 128 ranks, restart 9.0 -> 9.3 Myr (the PR #73 window)
Set-up: `/data/currodri/turbhii/mpi3_20261022/{A,B,C,D}_*`, run one after another at nice 19 with the machine
emptying (load 90 -> 1).

| run | `cool_queue` | `cost_measured` | wall | cooling Σmax / Σmean | wait after cooling |
|---|---|---|---|---|---|
| (#73's base, 2026-10-21) | 0 | 0 | 1,333 s | 1,281 / 147 s | 87% of the wall |
| A | 0 | 1 (#73) | 316 s | 280 / 149 s | 135 s |
| **C** | **2** | **0** | **212 s** | **188 / 187 s** | **3.9 s** |
| B | 2 | 1 | 225 s | 191 / 189 s | 6.2 s |
| D | 2 + `shared_tables` | 1 | 225 s | 191 / 189 s | 6.1 s |

- **Speed-up:** the queue alone is 6.3x faster than no balancing and 1.49x faster than #73. With the queue, #73's
  rebalance does not help on one node.
- **Where a rank's cooling time goes (C):** 173 s of solving and 13 s (7%) of queue overhead, against 1.7% at 32
  ranks.
  - The overhead is the stealing: 7.7 M chunks taken in the window. Each steal scans the 128 headers and contends on
    one counter.
  - The chunk size is the lever (runs below).
- **Mean solve per rank:** 173 s, against 149 s of cooling with #73. The likely causes are the per-cell warm-start
  reset, and every core computing at once (a lower all-core clock and shared L3/memory bandwidth) where the waiting
  ranks used to sit in `MPI_WAITALL`. Not separated yet.
- **Memory:** the summed PSS is 27.5 GB with `shared_tables` (D) against 33.8 GB without (B). That is 6.3 GB per
  node, 49 MB per rank.
- **Chunk size does not matter:** C with `cool_queue_chunk` = 1e-3 and 5e-3 runs in 212 and 214 s, with 3.0 M and 1.2 M
  chunks taken, against 7.7 M at the default. The stealing is not the overhead. The ~14 s between cooling and solving
  is the gather and scatter of `coolfine1`, which mode 0 counts in its cooling too.
- **The cost of the solve is the cold charge root:**
  - mode 0 with the per-cell reset of the queue (`cool_queue_debug = 2`): mean cooling 174 s against 147 s;
  - resetting only the ODE step: 147 s, no cost;
  - resetting only the WDB06rt charge root (`rtg_Zguess`): 174 s.
  - A cold root search brackets each bin over its whole charge range (up to 43,011 states), where the previous cell's
    root starts the search next to it.

### Each cell keeps its own charge roots (`cool_queue_zguess_ivar`)
- **How:** with a build of `NVARNOADVECT >= 1 + ndust` (5 in `bin_dust`), the cell's last WDB06rt roots live in the
  non-advected scalars `cool_queue_zguess_ivar .. +ndust-1`.
  - They go into the record, start the root search, and come back with the new roots.
  - They move with the cell through refinement, the load balance and restarts; cells never solved start from
    `noadvect_init` = -1.
  - A cell's result depends only on its own history, so it stays independent of the rank and of the order.
- **Recommended:** `NVARNOADVECT = 5`, `cost_ivar = 1`, `cool_queue_zguess_ivar = 2`.
- **Determinism checks:**
  - t1_pdr_calima, np 4: mode 1 in natural order, mode 1 in reverse order, and mode 2 twice are bit-identical in all
    82 fields;
  - mode 0 against them: 8.5e-10.
- **128 ranks, 9.0 -> 9.3 Myr:**

  | run | wall | cooling Σmax / Σmean | Σ solve per rank |
  |---|---|---|---|
  | (no balancing, 2026-10-21) | 1,333 s | 1,281 / 147 s | |
  | #73 alone (A, A2) | 301-316 s | 267-280 / 147-149 s | |
  | queue, cold roots (C) | 212 s | 188 / 187 s | 173 s |
  | **queue + own roots (C_zg)** | **181 s** | **158 / 157 s** | **143 s** |
  | queue + own roots, repeated (C_zg2) | 182 s | 158.5 / 157.2 s | 143 s |
  | queue + own roots + #73 (B_zg) | 198 s | 162 / 159 s | |

  - C_zg and C_zg2 are bit-identical in all 85 fields (262,144 cells).
  - C against C_zg differ by up to 1.3e-6 in a few rare ions (He III, Ne VI, O VI): the root-solver tolerance,
    accumulated over 63 steps.
  - Overall: 7.4x over no balancing and 1.7x over #73. A rank's solve time (143 s) is now below #73's mean cooling
    (147 s).

### Stress tests (megatron, 128 ranks, lv6 box from t = 0; `/data/currodri/turbhii/mpi3_stress`)
"Queue" below means `cool_queue = 2` with `cool_queue_zguess_ivar = 2` (`NVARNOADVECT = 5`); runs were capped at
30 min.

**Corner source:** the box's source moved from the centre to (2, 2, 2) pc. The young H II front then lies on one or
two ranks.

| run | wall at step 15 | wall at step 29 | cooling Σmax / Σmean at step 15 |
|---|---|---|---|
| queue | **78 s** | **118 s** | 70.6 / 70.1 s |
| queue + #73 | 79 s | 123 s | 70.1 / 69.3 s |
| #73 alone | 275 s | 344 s | 266 / 64 s |
| no balancing | 1,118 s | (step 16 still running at the 30-min cap) | 1,089 / 68 s |

- At step 15 the queue is 14x faster than no balancing and 3.5x faster than #73. #73 loses most of its time in the
  first steps, where the front grows faster than a rebalance every 2 steps can follow.

**Flickering sources:** six clusters of the 7 group sources (0.4x the box's luminosity each) at fixed random
positions. A cluster is switched on or off by a hash of its index and `nstep_coarse/2`, so the cost moves every 2
steps, in step with `nremap = 2`. The patch is test-only (`mpi3_stress/patch_flicker/rt_init.f90`), not in the
branch. 12 steps:

| run | wall | cooling Σmax / Σmean |
|---|---|---|
| queue | **189 s** (64 s to step 8, before the disturbance below) | 75 / 73 s |
| queue + #73 | 864 s (disturbed throughout) | 130 / 119 s |
| #73 alone | 1,434 s | 709 / 91 s |
| no balancing | 1,458 s | 793 / 100 s |

- The cost a step moves to is not where it was, so #73's measured-cost balance gains almost nothing (7.8x imbalance
  against 7.9x without it). The queue stays balanced (1.00-1.09).
- **The disturbance:** during these runs megatron was disturbed from outside. `polkitd` held 42 GB (a leak), with
  ~21 GB free and 4.8 GB of swap in use.
  - Every synchronisation stalled (the wait after cooling grew from 0.07 s to more than 1 s per call) while the
    cooling stayed balanced.
  - Wall times of the runs after ~20:28 are therefore inflated. The imbalance ratios (Σmax/Σmean) are the robust
    figures.
  - The corner runs ran before it (only the no-balancing run's last minutes overlap it, and its mean cooling per rank
    matches the queue's).
- **Gather and scatter** (`cool_stats`, in the queue line): ~0.1 s per coarse step each, balanced. The work left on
  the owner rank is negligible.

### Phase 3 (shared tables, `shared_tables = .true.`)
- **Converted:**
  - RTZ `high_t_cooling_rates` and `_opt` (18.4 MB);
  - CALIMA `collisional_tab` and `sputtering_tab` `%tab2d` (14.2 + 12.9 MB with `dust_coll_cooling` and RM2026
    sputtering);
  - `isrf_q` (6.6-14.9 MB).
- **Results bit-identical with sharing on and off:** t1_pdr_calima (np 4), t1_pdr_isrf (np 4, WDB06isrf),
  t1_eq_calima (np 2, sputtering + collisional).
- **Saving:** in the turbulent box that is ~52 MB per rank, ~6.6 GB per 128-rank node.
- **Added since:**
  - **WDB06rt `rtg`:** `K`, `H`, `delta`, `Ptr`, `fcool` and `Ebar` (~9.6 MB per rank with 4 bins). The axes and small
    arrays are broadcast; the caches that change at run time (`KgT`, `HgT`, `mKT`, `mHT`, `uni`, `egy_used`) stay per
    rank.
  - **`SED_table`:** every rank still builds it (the age stride and the `MPI_ALLREDUCE`), but only the node reader
    stores it.
  - **The MIST `sample`:** real(4), ~99 MB. The node reader reads the unformatted file; the axes are broadcast.
- **Checks:**
  - tier 1 with sharing off: byte-identical to staging (`runs/mpi3_tab`);
  - `rtg`, sharing on against off: bit-identical (t1_pdr_calima, WDB06rt, np 4);
  - MIST: bit-identical (t1_sink_accretion, np 4).
  - `SED_table`: on staging every SED initialisation crashes in `write_SEDtable` (a column past `nv`). PR #76 fixes
    the column layout, and it applies cleanly on top of this branch. With it applied, sharing on and off are
    bit-identical, with `sedprops_update = 0`, in the RTZ-only ladder (np 2) and in t1_pdr_calima (np 4).
