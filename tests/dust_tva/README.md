# The dust_tva tests

The tests of the CALIMA dust dynamics: the terminal-velocity drift (TVA), radiation pressure on the
grains and the trapped IR. `calima/README.md` describes the physics and the code; this file how each
test is set up, its namelist variants and what it measures.

## The tests

| test | what it checks |
| --- | --- |
| `tests/dust_tva/dustyirtrap` | the trapped-IR drift term in isolation, against a closed form, over 2 decades of `rho_d`. Freezes `P_trap`, so it asserts the FIRST step only. |
| `tests/dust_tva/dustyslab` | the full UV -> IR -> trapped chain with the real chemistry: UV injected at the left edge streams through dust-poor gas, is absorbed in a dust slab at 3-7 pc, and is reprocessed into the IR. |
| `tests/dust_tva/dustylev` | dust levitation (RT15 sec. 3.7 / Davis+2014): the ABSOLUTE normalisation of the IR radiation force against gravity, measured as an Eddington ratio. The only test that pins down the overall scale rather than a ratio or a gradient. |
| `tests/dust_tva/levgas_rt15` | RT15 sec. 3.7 itself: the bound exponential atmosphere illuminated from below, plain RAMSES-RT with their KT13 opacity, every dimensional parameter matched. Reproduces their figs. 9-11 and 13 diagnostics and fig. 10 projections. |
| `tests/dust_tva/levgas_rt15_dust` | the same run compiled with CALIMA and the dust as a drifting fluid (`dust_tva` + `dust_radpressure`, one bin, every other dust process off), with `eps` tuned so the dust reproduces RT15's opacity exactly. |


## `dustyslab`

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
  It reaches `tau_cell = 0.075` and `f_trap = 2.7e-5`: the reprocessed IR streams out of the slab
  instead of piling up where it is created. The quantitative RT15 eq. 50 check of the trapped
  fraction is `dustylev_trap.nml`, which reaches `f_trap = 0.42`.

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
- **Use the Rosseland table the code writes** (`SEDtables/rosseland_mean_DustBin_NN.list`; `utils/py/calima/calima_tables.py` reads and checks all of them)
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

## The levitation test (`tests/dust_tva/dustylev`)

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


With `gamma_rad(1) = 4/3` (`rt/rt_init.f90`) the trapped pressure is `P_trap = E_trap/3`, RT15
eq. 46, with `E_trap` the physical energy density that `Np2Ep` (in `cooling_fine`) returns; the
variant can check that ratio directly, with `E_trap` rebuilt from `f_trap * e_gamma * Np * f_c`.

## RT15 sec. 3.7 itself: `levgas_rt15` and `levgas_rt15_dust`

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
from `rt_c_cgs(ilevel)`; `gamma_rad(1) = 4/3` does not depend on the light speed.

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

