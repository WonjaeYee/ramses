# CALIMA dust module layout

This folder contains the CALIMA dust and PAH physics used by RAMSES. The code is split into small modules so that shared state, table loading, physics kernels, and runtime wrappers stay separate.

## Module hierarchy

```mermaid
graph TD
    dustbin_types --> dust_commons
    dust_utils --> dust_commons
    dust_commons --> dust_init
    dust_commons --> dust_rates
    dust_commons --> dust_charging
    dust_commons --> dust_surface_chemistry
    dust_commons --> dust_dynamics
    dust_commons --> dust_cooling
    dust_commons --> dust_interface
    dust_commons --> dust_solver
    dust_commons --> dust_seeding
    dust_commons --> dust_optics

    dust_rates --> dust_dynamics
    dust_rates --> dust_charging
    dust_cooling --> dust_charging
    dust_interface --> dust_charging
    dust_interface --> dust_surface_chemistry
    dust_solver --> dust_rates
    dust_solver --> dust_dynamics
    dust_solver --> dust_cooling
    dust_solver --> dust_charging
    dust_solver --> dust_surface_chemistry
```

The practical layering is:

1. dustbin_types defines the derived types used everywhere else.
2. dust_utils provides generic interpolation and helper routines.
3. dust_commons holds the shared dust and PAH parameters, bin properties, counters, and global switches.
4. The physics modules build on that shared state: dust_rates, dust_charging, dust_surface_chemistry, dust_dynamics, dust_cooling, dust_optics, and dust_seeding.
5. dust_init performs startup validation and table loading.
6. dust_interface exposes the dust helper routines used by the RTZ cooling path.
7. dust_solver is the main time-integration driver for dust chemistry.

## File guide

| File | Purpose |
| --- | --- |
| dustbin_types.f90 | Defines the core derived types: DustTable, DustChemistryInfo, DustBin, and PAHBin. Also provides init/reset methods for the reusable dust chemistry workspace. |
| dust_utils.f90 | Generic helpers used throughout CALIMA, including search/interpolation routines, turbulence utilities, and small math helpers. |
| dust_commons.f90 | Central shared module for dust/PAH flags, namelist-controlled parameters, per-bin properties, global counters, and helper data used by the rest of CALIMA. |
| dust_rates.f90 | Computes the characteristic timescales for dust and PAH processes: accretion, sputtering, coagulation, shattering, RATD, sublimation, evaporation, freezing, and related update switches. |
| dust_charging.f90 | Implements dust charge distributions and mean charge estimates, plus Coulomb focusing helpers and charge-state mixing routines. |
| dust_surface_chemistry.f90 | Handles H2 formation on dust grains, sticking and recombination efficiencies, and related surface chemistry fits. |
| dust_dynamics.f90 | Provides grain relative velocities and dust destruction in shocks, including PAH-specific shock destruction. |
| dust_cooling.f90 | Computes dust collisional heating/cooling, including the BH80 low-temperature branch and its cached prefactors. |
| dust_photophysics.f90 | Module dust_optics. Reads dust and PAH optical tables, dielectric tables, and stores the cross-section data in dust bin structures. |
| dust_seeding.f90 | Supplies dust and PAH source terms from stellar ejecta and winds, using limiting-element logic for condensation. |
| dust_init.f90 | Startup module for CALIMA. It validates the namelist, builds dust and PAH bin metadata, allocates the shared dust workspace, loads tables, initializes caches, and prints the active configuration. |
| dust_interface.f90 | Thin public wrapper layer for radiative dust work. It computes local anisotropy, radiative rates, precooling terms, and the combined dust cooling interface used by RTZ. |
| dust_solver.f90 | Main chemistry solver. It decides whether dust should be updated, computes timescales, advances dust and PAH abundances, and accumulates process counters. |

## External wiring

The CALIMA modules are not standalone. They are pulled into the RAMSES startup and cooling paths from outside this folder.

1. hydro/read_hydro_params.f90 reads the calima_params namelist when CALIMA is enabled. That is where the dust and PAH switches, bin properties, models, timescales, and external table directory are set from input.
2. hydro/read_hydro_params.f90 also calls check_params_dust so incompatible dust configurations fail before the run starts.
3. amr/init_time.f90 calls init_CALIMA_dust during startup. That routine is where the dust bin metadata is built, the reusable dust_helper workspace is allocated, the dust and PAH tables are loaded, and the BH80 cache is initialized.
4. rtz/rtz_cooling_module.f90 resets dust_helper for each cell, fills it with the current dust and PAH state, and then calls compute_dust_rad_rates and compute_dust_precool in the CALIMA branch of the cooling step.
5. rtz/rtz_coolrates_module.f90 calls compute_dust_coolrates through dust_interface so dust contributions are included in the cooling-rate solve.
6. hydro/cooling_fine.f90 enables CALIMA-specific helpers in the cooling operator-splitting path, including the turbulence sigma helper from dust_utils.

## Runtime flow

The typical order is:

1. Parameters are read from the namelist.
2. CALIMA configuration is validated.
3. init_CALIMA_dust builds the dust and PAH bin structures and loads the tables.
4. During each cooling update, the current cell state is copied into dust_helper.
5. dust_interface and dust_solver compute the radiative, collisional, chemical, and destruction terms.
6. The updated dust, PAH, and gas quantities are written back to the hydro and RT state.

## Notes

- The dust optical tables are read from the directory given by dust_tables_dir.
- The shared per-rank chemistry workspace is the DustChemistryInfo instance dust_helper.
- Most CALIMA routines depend on dust_commons for bin metadata and on dustbin_types for the storage layout of tables and helper arrays.

## Dust dynamics (TVA)

`dust_tva=.true.` advects the dust bins relative to the gas at the terminal velocity
(Lebreuilly+2019), as an operator-split upwind step applied to `unew` after the Godunov sweep
(`amr_step.f90`, `dust_diffusion_fine` / `dust_push_fine`). Things to be aware of:

- **PAH bins do not drift.** `dust_upwind_correct1/2` loop over the `ndust` dust bins only; the
  `npah` PAH bins stay perfectly coupled to the gas. `compute_gas_dust_radpressure_acc` does return
  a PAH acceleration, but it is not used, and PAH mass is excluded from `eps_tot` and from the
  barycentric acceleration. Defensible for PAHs, which are small and well coupled, but it is an
  approximation, not an accident. `check_params_dust` warns when `npah > 0` and TVA is on.
- **Gas-phase metals do not follow the dust.** Only the `ndust` density scalars are advected;
  `imetal` is untouched. Dust carrying C/O/Mg/Si/Fe across a cell boundary does not move the
  corresponding gas-phase element, so the per-element budget drifts over time. Watch the
  `dust_log` mass-conservation output.
- **`tva_wmax_cs`** caps `|w_drift|` at that multiple of the local sound speed, in both
  `get_dust_courant_dt` and the flux routines. TVA assumes Stokes << 1, which fails once the drift
  approaches `c_s`; since `t_s ~ 1/rho_gas`, a single hot/diffuse cell would otherwise drive the
  global timestep to zero. Set `<= 0` to disable (not recommended). Clipping is counted and
  reported alongside the "dust drift sets dt" message.
- **Enabling TVA changes the pure hydro** even before any drift matters: `ctoprim`
  (`hydro/umuscl.f90`) and `cmpdt` (`hydro/courant_fine.f90`) switch the EoS from the mixture
  density to the gas density `rho_mix*(1-eps_tot)`. A TVA run is therefore not bit-comparable to
  the same setup with TVA off.
- **`condinit_kind`** values `dustydiffuse`, `dustyshock`, `dustyblast1d`, `dustyspress` and
  `dustygauss` select **test** branches that hardcode the stopping time and/or the thermodynamics.
  These are resolved once into `tva_test_mode` at startup and warned about. Production ICs must use
  a different `condinit_kind` (the default `region`).

## Trapped IR radiation pressure on the dust

RAMSES-RT splits IR radiation pressure into two complementary channels whose sum is a single force
(Rosdahl & Teyssier 2015, eq. B8):

```
kappa_R rho F / c  =  kappa_R rho F_s / c  -  (1/3)(c~/c) grad(E_t)
   (total IR)            (streaming)            (trapped / diffusion)
```

The streaming half is the ordinary flux term, computed per bin in `dust_radpressure.f90`. The
trapped half lives in the `NENER` slot `iIRtrapVar`, so the Godunov solver already applies
`-grad(P_trap)/rho_mix` to the **mixture barycentre**. Physically that force acts on the dust,
which carries the IR opacity, so TVA adds the *differential* part to the drift driver:

```
D_k += - grad(P_trap) * ( s_k/rho_k - 1/rho_mix ),   s_k = chi_R,k / chi_R,tot
```

using `P_trap = (gamma_rad-1)*E_trap`, deliberately the same expression the Riemann solver uses so
the barycentric part cancels exactly. Things to be aware of:

- **IR opacities are the T_rad-dependent means, not the band-weighted ones.** For the IR group,
  `group_cs*_dust/pah` are the wrong spectral weight (they are averaged over the group band with a
  1e5 K blackbody or the stellar SED). The IR column of `csr_dust`/`csa_dust` is therefore
  overridden per substep in `cool_step` with the Rosseland (flux) and Planck (absorption) means at
  the local radiation temperature, from `dustbins_props(k)%Rosseland_tab`/`%Planck_tab` and the new
  per-charge-state PAH equivalents. This is CALIMA's analogue of RAMSES-RT's `is_kIR_T`, taken from
  the real grain optics instead of a `kappa ~ T_rad^2` power law. The same cross sections feed the
  TVA force, the trapping optical depth and `s_k`, so the channels cannot disagree at `tau_c ~ 1`.
- **`D_k` is proportional to `1/rho_k`, which is a destabilising feedback**: a dust-poor cell drifts
  faster and evacuates further. In a real run it self-regulates, because `tau propto chi_R propto
  rho_d`, so as dust leaves, `f_trap = exp(-1/tau) -> 0` and the force switches itself off. If you
  prescribe `E_trap` externally (as `condinit_kind='dustyirtrap'` does) that regulation is absent,
  which is why that test asserts the first-step drift rather than multi-step advection.
- **The trapped-IR flux is anti-diffusive, so boundary blemishes grow.** For a single bin the mass
  flux is `F = rho_d w_d = (-grad P_trap) t_s (1 - eps)^2`, which *decreases* with `rho_d`
  (measured `dF/drho_d ~ -1.9e-3 < 0`) even though `w_d > 0`. The flux routine picks the donor cell
  on `sign(w_face)`, i.e. on the material velocity, which is the opposite of this flux's
  characteristic direction. A perturbation therefore grows instead of damping. In
  `tests/dust_tva/dustyirtrap` the seed is the domain boundary: the zero-gradient ghost fill makes
  the centred difference effectively one-sided there, so the two edge cells get
  `grad(P_trap) = -5e-6` instead of `-1e-5` — exactly half — which is an O(1) error in their drift.
  A 2*dx cell-to-cell mode then eats inward from both boundaries (reaching x ~ 7.5 of a 10 pc box
  by t = 0.27), while the interior between the fronts stays smooth at the 1e-4 level. See
  `dustyirtrap_runaway.png`. The halved edge gradient is ordinary RAMSES zero-gradient BC behaviour
  and affects the pre-existing `grad(P_gas)` term identically; what is specific to the trapped-IR
  term is `dF/drho_d < 0`, which converts that blemish into a growing mode.
  **Open question for production:** there `E_trap` is not independent of the dust, since
  `tau propto chi_R propto rho_d` makes `f_trap = exp(-1/tau)` increase with `rho_d`, which may
  well restore `dF/drho_d > 0` and with it the correct upwind direction. That could not be tested
  here because the RTZ chemistry path stalls (see the plan notes), so it is unverified. If it turns
  out `dF/drho_d < 0` in a real run, the trapped-IR contribution to the dust flux should be
  upwinded on the characteristic direction, or given a Rusanov/LLF flux instead.
### Tests

| test | what it checks |
| --- | --- |
| `tests/dust_tva/dustyirtrap` | the trapped-IR drift term in isolation, against a closed form, over 2 decades of `rho_d`. Freezes `P_trap`, so it asserts the FIRST step only. |
| `tests/dust_tva/dustyslab` | the full UV -> IR -> trapped chain with the real chemistry: UV injected at the left edge streams through dust-poor gas, is absorbed in a dust slab at 3-7 pc, and is reprocessed into the IR. |
| `tests/dust_tva/dustylev` | dust levitation (RT15 sec. 3.7 / Davis+2014): the ABSOLUTE normalisation of the IR radiation force against gravity, measured as an Eddington ratio. The only test that pins down the overall scale rather than a ratio or a gradient. |
| `tests/dust_tva/levgas_rt15` | RT15 sec. 3.7 itself: the bound exponential atmosphere illuminated from below, plain RAMSES-RT with their KT13 opacity, every dimensional parameter matched. Reproduces their figs. 9-11 and 13 diagnostics and fig. 10 projections. |
| `tests/dust_tva/levgas_rt15_dust` | the same run compiled with CALIMA and the dust as a drifting fluid (`dust_tva` + `dust_radpressure`, one bin, every other dust process off), with `eps` tuned so the dust reproduces RT15's opacity exactly. |

`dustyslab` has three namelists. A resolved UV attenuation front and
IR-optically-thick cells cannot coexist: with `kappa_UV/kappa_R ~ 1e2-1e4` and
`L_slab/(1.5 dx) ~ 34`, the two optical depths are locked together as
`tau_UV(slab) ~ (kappa_UV/kappa_R) * (L_slab/1.5dx) * tau_cell,IR`, i.e. a
factor 5e3 at best. Any appreciable trapping therefore forces the UV to be
absorbed in a small fraction of a cell. The three namelists pick points along
that line:

| namelist | rho_gas [g/cm3] | tau_UV(slab) | tau_cell,IR | achieved f_trap |
| --- | --- | --- | --- | --- |
| `dustyslab.nml` | 1e-22 | 3.5 | 3e-11 | 0 |
| `dustyslab_dense.nml` | 1.23e-17 | 4.3e5 | 0.055 | 9e-9 |
| `dustyslab_sil.nml` | 2.59e-18 | 2.3e3 | 0.075 | 2.7e-5 |

- `dustyslab.nml` -- the fiducial case. `tau_UV = 3.5` with an e-folding length
  of 15 cells, so the attenuation front is well resolved and matches
  `exp(-tau)` to 10%; the direct UV acceleration collapses 29x across the slab
  and the reprocessed IR peaks inside it. No trapping at all. Also shows the
  radiation-driven dust push: the illuminated face is ablated from ~0.25 Myr.
  `tend=2` Myr, ~45 s.
- `dustyslab_dense.nml` -- 5 decades denser. Reaches `tau_cell = 0.055`, so the
  trapped fraction is still only 9e-9: useful as the UV-thick limit, but it does
  NOT exercise the trapped-IR channel. `tend=0.5` Myr, ~20 s.
- `dustyslab_sil.nml` -- **the variant that actually exercises trapping**, and
  the least extreme one that can. Two changes do the work: the slab is moved to
  bin 3 so it uses the 0.005 um silicate table, whose `kappa_UV` is 38x below
  the 0.01 um graphite of bin 1 (`condinit_kind='dustyslabsil'`, `NDUST=3`,
  `NDCHEMTYPE=2`); and the UV flux is raised to 1.5e13 ph/cm2/s, which warms the
  grains to `T_rad = 33 K` where `kappa_R = 15` instead of ~1. Since
  `kappa_R ~ T_rad^2` and `T_rad ~ F^(1/4)`, `tau_cell ~ F^(1/2)`, so flux is a
  cheaper lever than density. `tend=0.5` Myr, ~20 s.

  **This variant no longer reaches the trapping regime.** It did before the
  `cooling_fine` `fred` fix (see below): with the IR flux clamped to a reduced
  flux of `1/rt_c(ilevel)**2` the IR barely streamed, so it piled up where it
  was created, `T_rad` and hence `kappa_R` were inflated, and `tau_cell` reached
  0.14 with `f_trap = 4.7e-4`. With the IR streaming correctly the reprocessed
  energy leaves the slab, `tau_cell` falls to 0.075 and `f_trap` to 2.7e-5. The
  quantitative RT15 eq. 50 check now lives in `dustylev_trap.nml`, which reaches
  `f_trap = 0.42`.

#### The levitation test (`tests/dust_tva/dustylev`)

An isolated Gaussian dusty layer (`sigma = boxlen/10`, uniform `eps = 1e-2`) in
a uniform downward gravity field (`gravity_type=1`, `gravity_params(1) = -g`),
illuminated from `x=0` by a directed IR beam. Because the density falls to a
floor at BOTH boundaries there is no wall pressure force, so the box-integrated
momentum budget closes exactly and independently of how the layer restructures:

```
d<v>/dt = g*(f_E - 1),     f_E = (F_E/c)*(1 - exp(-tau_R)) / (g*Sigma)
```

In the optically thin limit this collapses to the mass-specific Eddington ratio
`f_E = kappa_R * eps * F_E / (c*g)`, with **no `Sigma` in it at all** -- which
is why it measures the absolute normalisation rather than a ratio. Deliberately
an isolated layer rather than the wall-supported exponential atmosphere of
`patch/rt/davis`: there the reflecting lower boundary carries part of the
reaction force and the net acceleration is no longer a closed form.

Run with `dust_radpressure=.false.`: levitation tests the BARYCENTRIC radiation
force, which is applied in the cooling step, so this keeps the TVA differential
term out of the measurement. The barycentric force still acts -- it comes from
`cooling_fine`, not from the TVA driver.

| namelist | target | measured `f_E` | error in net accel |
| --- | --- | --- | --- |
| `dustylev_sub.nml` | `f_E = 0.5` | 0.50000 | 2.8e-5 `g` |
| `dustylev_bal.nml` | `f_E = 1.0` (hovers) | 1.00004 | 6.4e-5 `g` |
| `dustylev_super.nml` | `f_E = 2.0` | 2.00022 | 1.6e-4 `g` |
| `dustylev_trap.nml` | `tau_R = 21`, `f_trap = 0.42` | 0.738 (see below) | -- |

`dustylev_bal` is the sharpest: the predicted net acceleration is zero against a
free-fall rate of `-g`, so any error in `kappa_R`, in the photon energy, or a
stray `c_red/c` shows up as a net drift of order `g`. It comes out at 6e-5 `g`.

Everything on the right-hand side is taken from the run's OWN state --
`F_E` from the dumped photon flux times the group energy in `info_rt`, and
`kappa_R` from the Rosseland table the run writes to `SEDtables/` -- so the
check never depends on reproducing the code's tables in Python.

Two things that will bite when sizing a variant:

- **`rt_init` RECOMPUTES `group_egy` from the band edges** (a 1e5 K blackbody
  weighting over `groupL0..groupL1`), overwriting the namelist value. For
  0.001-0.1 eV it comes out at **6.696e-2 eV, not 0.01** -- a factor 6.7 on the
  energy flux, which was worth a factor 16 on `f_E` (via `kappa_R(T_rad)`, since
  `T_rad ~ F^(1/4)`). Read `egy [eV]` back out of `info_rt_*.txt`.
- Because the streaming IR force uses `kappa_R(T_rad)` and `T_rad ~ F^(1/4)`,
  `f_E` is **not** linear in `F`: roughly `f_E ~ F^1.5`. The flux for a target
  `f_E` has to be solved numerically against the code's own Rosseland table.

`dustylev_trap.nml` raises the density to `nH = 6.8e5` so that `tau_R = 21` and
`tau_cell = 1.1`, giving `f_trap = 0.42`. Note the structural constraint
`tau_cell = 1.5*(dx/L)*tau_R`: at fixed resolution the only way to get
`tau_cell ~ 1` with a resolved layer is a large `tau_R`, which forces the
variant to be super-Eddington in the thin limit (`f_E,thin = 41`). It therefore
also uses `condinit_kind='dustylevtrap'`, which drops the ambient two more
decades so a super-Eddington ambient does not out-weigh the layer.

**What the trapping variant can and cannot test.** For an ISOLATED layer
`P_trap` falls to zero at both faces, so

```
integral(-dP_trap/dx) dx = P_trap(0) - P_trap(L) = 0
```

-- the trapped channel gives **no net force**, only internal redistribution
(forward in the interior, backward in the illuminated skin). RT15's own
levitation test sidesteps this with a wall-supported atmosphere whose reflecting
boundary carries the reaction. So the checks here are the well-posed ones:
`f_trap` against RT15 eq. 50 (`median achieved/predicted = 1.21` over 36 cells),
the net trapped-force integral being zero, and the fraction of the absorbed
momentum that actually reaches the gas -- **36.7%**, the rest being what the
trapping closure holds back from an isolated cloud.

It also gives a direct read-out of the `c_red/c` normalisation flagged below:

```
measured P_trap / E_trap = 3.333e-04 = f_c/3      (f_c = rt_c_fraction = 1e-3)
RT15 eq. 46 requires      3.333e-01 = 1/3
```

with `E_trap` rebuilt independently from `f_trap * e_gamma * Np * f_c`. `Np2Ep`
in `cooling_fine` already converts to a physical energy density, so
`gamma_rad(1)-1` (set to `rt_c_fraction/3` at `rt/rt_init.f90:245`) applies
`c_red/c` a second time. Reported, not changed.

#### RT15 sec. 3.7 itself: `levgas_rt15` and `levgas_rt15_dust`

`dustylev` above measures the absolute normalisation of the IR force but uses an
isolated layer, for which the trapped channel gives zero net force by
construction. These two folders instead run the experiment RT15 sec. 3.7 ran --
the Krumholz & Thompson (2013) / Davis et al. (2014) bound exponential
atmosphere illuminated from below -- as a matched pair:

| folder | build | what it adds |
| --- | --- | --- |
| `levgas_rt15` | `RT=1, RTZ=0, CALIMA=0` | the reference. Plain RAMSES-RT with the KT13 opacity law RT15 used. Compare directly against their published figures. |
| `levgas_rt15_dust` | `RTZ=1, CALIMA=1, NDUST=1, NPAH=0` | identical physics, plus the dust as a DRIFTING fluid: `dust_tva` (gas drag) and `dust_radpressure`. Every other dust process is off. |

Each has two namelists: the plain one at RT15's resolution (`2048^2` over
`L_box = 1024 h_*`, `200 t_*`, light-speed ramp from full `c` over 3e4 RHD
steps) which is a **cluster** run, and a `_quick` one (`128^2` over `256 h_*`,
`20 t_*`, constant `c_red`) for local smoke-testing.

**RT15's dimensional setup is matched exactly in both**, which is only possible
because the dust-to-gas ratio is tuned to reproduce their opacity:

| | RT15 | here |
| --- | --- | --- |
| `F_*` [erg/cm2/s] | 1.03e4 | 1.03e4 |
| `T_*` [K] | 82 | 82 |
| `g` [cm/s2] | 1.46e-6 | 1.46e-6 |
| `rho_*` [g/cm3] | 7.1e-16 | 7.1e-16 |
| `kappa_R,*` [cm2/g_gas] | 2.13 | 2.1248 |
| `h_* = c_iso^2/g` [cm] | 2e15 | 1.9897e15 |
| `t_* = h_*/c_*` [s] | 3.67e10 | 3.6916e10 |
| `Sigma` [g/cm2] | 1.4 | 1.4127 |
| `tau_*` | 3 | 3.0017 |
| `f_E,*` | 0.5 | 0.5000 |

CALIMA takes the IR Rosseland mean from the grain optics, 129.881 cm2 per gram
of **dust** at 82 K, so `dustylevatm_condinit` uses `eps = 0.016359` to make
`kappa` per gram of **gas** equal 2.1248 -- RT15's value. The only differences
between the two folders are then the temperature dependence of `kappa` (measured
optics vs a `T^2` fit) and the dust dynamics.

Measured, with the `_quick` namelist: `f_E,V(t=0) = 0.5012`, `tau_V(t=0) =
2.9996`, `tau_F/tau_V(t=0) = 1.0000`, injected flux `= F_*` to 4 digits. In 1D
the reference reproduces their evolution too: `f_E,V` spikes early and
trapped-dominated (RT15 fig. 13), `tau_V` peaks inside their plotted 6-16 band,
and `<v_y>` reaches a few `c_*`.

**Three settings that are not optional.** Each of them stopped the run or moved
a headline number, and each is documented in the namelists:

- `group_csn(1,:) = 0.` and `group_cse(1,:) = 0.` in the reference folder. With
  `nGroups=1`, `rt_init.f90:392-408` hands the single group the hardcoded
  **HI-ionising** defaults (`csn = 3.007e-18 cm2`, `egy = 18.85 eV`), and the
  band-based re-initialisation at `:483` sits inside the RTZ branch, so in a
  plain RAMSES-RT build only the namelist can clear them. A far-IR band alone
  does NOT give zero cross sections there (it does under CALIMA, which builds
  them from its own tables). Left at the default, the "IR" beam photoionises the
  gas: `xHII` relaxes to 2.2e-3 and the hydrogen ionisation limiter (`code=6`)
  pins the cooling substep at 1.4 s against a 3e6 s hydro step, so the run
  cannot advance at all. `davis.nml` zeroes both for exactly this reason.
- `rt_kIR_RT15 = .true.` selects RT15's exact eq. 79 -- the **gas** temperature,
  and no `exp(-T_R/1000 K)` sublimation cutoff. Both the `T_rad` choice and the
  cutoff were added to RAMSES-RT after RT15. With the cutoff on,
  `kappa_R(82 K) = 1.9575` instead of 2.1248, which puts `tau_*` at 2.76 and
  `f_E,*` at 0.4617 rather than 3 and 0.5.
- `exp_region(1) = 10.0`. With the default exponent a `'square'` region is an
  **ellipse**, so in >=2D the corner cells fall outside it and never receive
  `d_region`; the 2D profile came out with twice its intended scale height.

**The boundaries follow `patch/rt/davis` exactly.** `x` periodic for matter and
radiation; bottom reflective for matter and **emitting** for radiation via
RT15 eq. 83, `c~E_0 = F_* - F_y,1 + c~E_1`, implemented in each folder's local
`patch/rt_hydro_boundary.f90` (`PATCH = patch` in the Makefile, the same route
`patch/rt/davis` takes, so no shared code is touched). An `rt_nsource` region
was tried first and rejected: it imposes `Np` every RT subcycle and therefore
also destroys radiation arriving from above, collapsing the net injection to
**18% of `F_*`** once photons pile up. Top: Dirichlet, `rho = 1e-13 rho_*`,
`v = 0`, `T = 1e3 T_*`, zero radiation energy and flux.

Note the last of those: RT15's **text** says `T = 10^-3 T_*`, but their
`hydro_boundary.f90` sets `82*1d3*rho/scale_T2/2.33`, i.e. `10^+3 T_*` -- and
only that satisfies the "in pressure balance with the initial conditions" they
also state, since `1e-13 rho_* x 1e3 T_* = 1e-10 rho_* T_*` is the initial
density floor at `T_*`. The text has a sign typo. Either value gives a ghost
pressure 1e10 below the interior, and switching between them leaves the run
bit-identical.

**The initial density profile must be the exact cell average.**
`rho_bar = rho_*(exp(-h_lo) - exp(-h_hi))/dx`, as `patch/rt/davis/condinit.f90`
does. Point-sampling the exponential loses column density as the cell gets
thick -- 14% at `dy = 2 h_*` -- which moves `tau_*` off 3. The exact form is
independent of the grid.

**The light-speed ramp (RT15 sec. 3.7) is implemented.** They "start the
experiment at a full light speed and converge exponentially towards `c~` over
3e4 RHD time-steps... specifically to capture the sudden and short lived pile-up
of trapped photons". Two namelist parameters in `&RT_PARAMS`:

| | |
| --- | --- |
| `rt_c_ramp_nstep` | RHD steps to converge over; `<= 0` disables it (the default), so existing runs are unaffected |
| `rt_c_ramp_start` | initial light-speed fraction, e.g. 1.0 for the full `c` |

`rt_init.f90` saves the namelist `rt_c_fraction` as the target and
`rt_ramp_lightspeed` (called once per RHD step from `rt_step`, at `levelmin`)
walks it down as `f(n) = f_tgt (f_start/f_tgt)^(1-n/N)`, exactly `f_start` at
`n=0` and `f_tgt` at `n=N`. Everything downstream follows automatically:
`rt_c`/`rt_c2` through `update_rt_c`, the group cooling constants through
`updateRTGroups_CoolConstants` (already called at the top of
`rt_solve_cooling`), and `Np2Ep` because `cooling_fine` rebuilds it per call
from `rt_c_cgs(ilevel)`. **This is only well posed because `gamma_rad(1)` is now
4/3** -- with the old `rt_c_fraction/3 + 1` a time-varying `c` would have needed
a time-varying `gamma_rad`, precisely the incompatibility the note on that line
warned about.

Two things to know about it. It is **expensive**: at full `c` the RT Courant
condition throttles the hydro step to `dt ~ 2e-3 t_*` with ~1300 subcycles each,
about **400x more cost per unit physical time** than running at `c_red`
(measured 2.8 steps/min at `256^2` on 12 cores), which is exactly why RT15 needed
3e4 steps to cover "the initial few `t_*`". And `rt_n_region` must be set for the
light speed at `t = 0`, i.e. `F_*/(egy * rt_c_ramp_start * c)`, not for the
target -- `Np` is a photon number density in the reduced-`c` system, so using
the target value while starting at full `c` overstates the initial radiation
energy by `rt_c_ramp_start/rt_c_fraction` and gave `f_E,V(t=0) = 2148`.

Measured effect in 1D, with RT15's exact opacity: the early `f_E,V` peak is
**20.0** un-ramped, **15.1** ramping from `f_c = 0.03` and **14.1** from 0.1,
against RT15's **10** -- so the ramp moves the peak towards their value, and
starting from 0.03 costs only 2.6x.

**Caveats.** `rt_Tconst` pins `T_gas` in the CALIMA folder, whereas RT15 let the
gas heat via their eq. 77; the reference folder lets it evolve. `tau_F/tau_V` is
**not** 1 away from `t=0` even in 1D, because RT15 eq. 86 weights by the
cell-centred `Fp`, which is small in the diffusion limit (measured reduced flux
0.035 in the bottom row) -- so only the *further* drop in 2D is attributable to
their 'chimneys'. And do **not** run these with `static=.true.`: freezing the
hydrodynamics removes the `PdV` work that drains trapped energy, and because
`coolfine1` multiplies the streaming *flux* by `(1-f_trap)` every call the
escape is throttled to nothing, so `tau_V` and `f_E,V` diverge without bound
(measured `tau_V -> 185`, `f_E,V -> 2.4e6` by `t = 250 t_*`) instead of reaching
the equilibrium RT15 quote.

#### Four bugs these tests surfaced

- **`cooling_fine.f90`: the inlined `reduce_flux` clamped the IR flux by
  `rt_c(ilevel)**2`.** Upstream commit `7f8712f0` ("Addition of variable speed
  of light to ramses-rt") replaced
  `call reduce_flux(rtuold(il,iNp+1:iNp+ndim), rtuold(il,iNp)*rt_c)` with
  `fred = |Fp| / rtuold(il,iNp)*rt_c(ilevel)`, which Fortran evaluates as
  `(|Fp|/Np)*rt_c` -- so `fred` was `rt_c(ilevel)**2` too large and **every**
  cell was clamped to a reduced flux of `1/rt_c(ilevel)**2` (1.06e-5 at
  `f_c=1e-3` with pc/Myr/1e-22 units). The IR therefore carried essentially no
  momentum, and the levitation test came out in free fall. Two sites, one in the
  `rt_vc` block and one in the `rt_isIRtrap` block; both fixed to divide, which
  matches the canonical form in `rt/rt_godunov_fine.f90:121`. The severity
  depended on the user's choice of code units, which is the giveaway that it was
  not deliberate.
- **`dust_radpressure.f90:422,441,463` divide the streaming radiation force by
  `rt_c_code` instead of `c`.** `opacity_*_code` is a raw opacity with no
  `rt_c_cgs` factor, so `mom_fact = chi*E_gamma/rt_c_code` is `1/f_c` too large
  -- 1000x at `rt_c_fraction=1e-3`. The cooling path is correct: there the rate
  already carries a factor `rt_c_cgs` (`sigcr_dust = group_csr_dust*rt_c_cgs`)
  and the extra `one_over_rt_c_cgs` cancels it, leaving `chi*E_gamma/c`. This
  one is **not fixed** -- it changes the drift in every existing run and the
  `rt_isoPress` branch needs separate thought -- which is why `dustylev` runs
  with `dust_radpressure=.false.`.
- **`dust_dynamics.f90`: the dust flux loops ran the transverse index from 0.**
  `dflux`/`eflux`/`mflux` are dimensioned `(if1:if2, jf1:jf2, kf1:kf2)` = 1:3,
  but the face loop took its transverse bounds from `ilo = MIN(1,iu1+1)` = 0
  (since `iu1 = ju1 = -1` in any active dimension), so **every 2D or 3D run
  died** with `Index '0' of dimension 3 of array 'eflux' below lower bound of 1`.
  RAMSES's own `umuscl` uses `MIN(1,iu1+2)` for exactly this loop. Fixed by
  clamping the face-loop bounds to the flux arrays in
  `calculate_pure_drag_fluxes` and `calculate_drag_rad_fluxes`; only indices
  1..3 are ever read back, so nothing is discarded, and in 1D
  `jlo=jhi=klo=khi=1` already, so it is a no-op -- confirmed by `dustyirtrap`,
  `dustyslab` and all four `dustylev` variants reproducing unchanged.
- **`region_condinit`: a `'square'` region is an ELLIPSE unless `exp_region` is
  set.** The region radius is `r = (xn**en + yn**en)**(1/en)` with the default
  exponent, so in >=2D the corner cells of a box-filling region have `r > 1` and
  never receive `d_region`. The 2D atmosphere came out with twice its intended
  scale height until `exp_region(1) = 10` (the max norm) was set. Harmless in
  1D, which is why no existing test caught it. Not a code change -- just set
  `exp_region` in any >=2D namelist.
- Minor, not changed: `output_hydro.f90` labels the NENER slots
  `non_thermal_pressure_<ivar-3>` rather than `<ivar-nhydro>`, so the *same*
  variable is `_01` in 1D and `_02` in 2D. Analysis scripts should look the key
  up rather than hard-code it.

`plot-dustyslab.py` takes the namelist and, optionally, the bin index carrying
the slab: `python3 plot-dustyslab.py dustyslab_sil.nml 3`. It reads the RT
output by staging a scratch copy with `rt_*` renamed to `hydro_*`, so shared
`visu_ramses` is untouched.

Two traps to avoid when interpreting these:

- **`kappa_R` must be taken at the local `T_rad`, not at a nominal 30 K.** The
  code uses `T_rad = (E_IR*f_c/a_r)^(1/4)`, and with a reduced light speed the
  `f_c` factor pulls `T_rad` down hard -- to 3.3 K in the dense variant, where
  the Rosseland mean sits in the sub-mm and `kappa_R = 1.2` rather than 22.
  Assuming 30 K overestimates `tau_cell` by ~15x and will make you think
  trapping is active when it is not.
- **Use the Rosseland table the code writes** (`SEDtables/rosseland_mean_DustBin_NN.list`)
  rather than re-integrating the harmonic mean yourself; an independent
  re-integration differed by ~2x, which is enough to move `f_trap = exp(-1/tau)`
  by two orders of magnitude.

Two things worth knowing when reading these:

- **Run only as long as the comparison is valid.** The radiation field reaches
  steady state in ~2 hydro steps (~0.07 Myr; `median |dF/F| ~ 2e-3` between
  0.068 and 0.102 Myr). The dust then starts to move: in the fiducial variant
  `t_s = 5.3 kyr` and `w_d = 0.23 pc/Myr`, so the illuminated face is ablated
  from ~0.25 Myr and `tau_UV` has halved by 2 Myr. Analytic comparisons must
  therefore be made on an early snapshot, which is what `plot-dustyslab.py`
  does; later snapshots are for the dynamics only.
- **The UV band is 5-7 eV deliberately.** That is below every ionisation
  threshold RTZ tracks (Mg 7.65, Fe 7.9, Si 8.15, C 11.26, H 13.6 eV), so
  `group_csn` is identically zero and the gas is transparent -- `group_csn` is
  overwritten at `rt/rt_init.f90:470` and cannot simply be zeroed from the
  namelist. Note also that `rt_n_bound` injects only into rtvar 1, i.e. group 1,
  which `rt_isIR` reserves for the IR, so the UV inflow is a thin `rt_nsource`
  region at x=0 with `rt_src_group=2`.

**Watch the grain size against the optical tables.** They are matched to bins by
index alone (`averaged_cross_section_DustBin_<NN>.txt`), with no check against
the namelist `asize`. The shipped tables are 0.01 um graphite, 0.1 um graphite,
0.005 um silicate, 0.1 um silicate -- which matches `dustyshell_mulbin`'s
4-bin `asize` list, but a single-bin run must use `asize=0.01d0` to match
`DustBin_01`. `dustyspress.nml` sets `asize=0.1d0` with `NDUST=1`, so every
kappa there is wrong by ~1000x (m_grain scales as a^3). With the size matched,
CALIMA's `kappa_R` at 30 K is 0.22 cm2/g_gas against RT15's 0.28 -- good
agreement.

- Note also that a *discontinuous* prescribed `E_trap` (e.g. a linear ramp with periodic
  boundaries, which does not wrap) produces a large spurious `grad(P_trap)` at the seam and
  destabilises the run outright.
- **The drift can be large by construction**: a ~1% dust mass fraction absorbs ~100% of the IR
  momentum, so `w ~ t_s |grad P_trap| / rho_d`. At `tau_IR ~ 1` this is tens of km/s, well outside
  the `Stokes << 1` regime TVA assumes, so `tva_wmax_cs` will clip. Watch `tva_nclip`.
- **PAHs contribute to `chi_R,tot` but do not drift**, so their share of the trapped force simply
  stays on the barycentre. The TVA closure conserves barycentric momentum regardless.
- Requires `NENER>=1`, `rt_isIR`, `NDUST>0` and `rt_flux_scheme='glf'` (RT15 footnote 3 — the
  trapped/streaming partition is matched to the GLF numerical diffusion, so HLL is inconsistent).
  `check_params_dust` enforces all four.