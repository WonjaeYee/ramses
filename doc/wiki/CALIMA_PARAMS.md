# CALIMA Parameters

This set of parameters, contained in the namelist block `&CALIMA_PARAMS`, controls dust and PAH physics in RAMSES simulations with the CALIMA (Chemistry And Lyman-alpha Interactions for Molecules and Atoms) dust module. CALIMA models grain growth, sputtering, coagulation, charging, and chemical interactions in dusty environments.

For detailed descriptions of the concepts and models described here, see:

* [[1] CALIMA: On-the-fly dust and PAH evolution for radiation-hydrodynamics galaxy formation simulations](https://ui.adsabs.harvard.edu/abs/2026arXiv260221790R/abstract)
* [CALIMA module documentation](../README.md) in the calima folder
* [Dust module source code](../../calima/)

## Dust Physics Flags

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `dust_log=.false.`            |  `logical`    | Activate dust logging to standard output for debugging. |
| `dust_10percent=.true.`       |  `logical`    | Activate the 10% rule for the chemistry solver, which limits large timestep changes. |
| `dust_only_rtadv=.false.`     |  `logical`    | Only activate dust chemistry when radiative transfer advection is on. |
| `dust_eq_test=.false.`        |  `logical`    | Activate dust equilibrium test parameters (for testing purposes). |
| `dust_SNdest=.false.`         |  `logical`    | Enable dust destruction in supernova (SN) explosions. |
| `dust_inSN=.false.`           |  `logical`    | Inject dust in core-collapse (SNII) explosions. |
| `dust_inSNIa=.false.`         |  `logical`    | Inject dust in thermonuclear (SNIa) explosions. |
| `dust_inSW=.false.`           |  `logical`    | Inject dust in AGB star winds. |
| `dust_accretion=.false.`      |  `logical`    | Activate grain growth by accretion of gas-phase atoms. |
| `dust_sputtering=.false.`     |  `logical`    | Activate grain destruction by thermal (ion) sputtering. |
| `dust_sputtering_charge=.false.`  |  `logical`    | Include dependence of thermal sputtering on grain and ion charge. |
| `dust_acc_coulomb=.false.`    |  `logical`    | Compute on-the-fly Coulomb enhancement of refractory material accretion. |
| `dust_coagulation=.false.`    |  `logical`    | Activate grain coagulation (collision and sticking). |
| `dust_coagulation_boost=.false.` |  `logical`    | Activate boost of coagulation in high-density gas. |
| `dust_shattering=.false.`     |  `logical`    | Activate grain shattering (fragmentation from collisions). |
| `dust_shattering_all=.false.` |  `logical`    | Include shattering caused by collisions of all grain sizes (otherwise only large grains shatter). |
| `dust_shattering_dest=.false.` |  `logical`    | Activate destruction of dust mass via shattering fragmentation. |
| `dust_shattering_SN=.false.`  |  `logical`    | Redistribute grain mass in SN shocks via inertial sputtering. |
| `dust_ratd=.false.`           |  `logical`    | Activate destruction of dust grains by radiation-induced thermal disruption (RATD). |
| `dust_coll_cooling=.false.`   |  `logical`    | Activate dust collisional cooling (energy transfer via grain-gas collisions). |
| `dust_coll_lowT=.false.`      |  `logical`    | Activate low-temperature dust collisional heating using Hollenbach & McKee (1980) model. |
| `dust_coll_charge=.false.`    |  `logical`    | Include dependence of dust collisional cooling on grain and ion charge. |
| `dust_pe_heating=.false.`     |  `logical`    | Activate photoelectric heating by dust grains (energy transfer from UV photons). |
| `dust_pe_heating_isrf=.false.` |  `logical`    | Use simple dust PE heating based on an averaged interstellar radiation field (ISRF) G₀. |
| `ratd_only_rtadv=.false.`     |  `logical`    | Only allow RATD when radiative transfer advection is active. |
| `poppe_ice_enhancement=.false.` |  `logical`    | Apply empirical enhancement factor to coagulation threshold due to icy mantles. |
| `H2ondust=.false.`            |  `logical`    | Activate molecular hydrogen (H₂) formation on dust grain surfaces. |
| `dust_turbulent_model=.false.` |  `logical`    | Activate subgrid model for turbulent shattering and coagulation. |

## PAH Physics Flags

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `pah_accretion=.false.`       |  `logical`    | Enable PAH growth by accretion of gas-phase carbon atoms. |
| `pah_acc_spu=.false.`         |  `logical`    | Activate destruction of PAHs by accretion of ionized carbon (C⁺). |
| `pah_coalescence=.false.`     |  `logical`    | Activate coalescence of PAHs into small carbonaceous grains. |
| `pah_freezing=.false.`        |  `logical`    | Activate freezing of PAHs onto carbonaceous grain surfaces. |
| `pah_desorption=.false.`      |  `logical`    | Activate desorption of frozen PAHs from carbonaceous grain surfaces. |
| `pah_uv_destruction=.false.`  |  `logical`    | Activate UV sublimation/destruction of PAH molecules. |
| `pah_sn_destruction=.false.`  |  `logical`    | Inertial and non-thermal destruction of PAHs in SN shocks. |
| `pah_cluster_evaporation=.false.` |  `logical`    | Evaporation of PAH clusters to form small PAHs via UV photon absorption. |
| `pah_AGBwinds=.false.`        |  `logical`    | Inject PAHs during AGB wind episodes. |
| `pah_sputtering=.false.`      |  `logical`    | Enable ion and electron sputtering (destruction) of PAHs. |
| `pah_pe_heating=.false.`      |  `logical`    | Activate photoelectric heating by PAHs. |
| `pah_pe_heating_isrf=.false.` |  `logical`    | Use simple PAH PE heating based on averaged ISRF. |
| `pah_pe_nolyman=.false.`      |  `logical`    | Deactivate the 13.6 eV limit for PAH photoelectric heating. |
| `H2onpah=.false.`             |  `logical`    | Enable formation of H₂ molecules on PAH surfaces. |

## Dust Modelling Options

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `sputtering_model='Tsai1998'` |  `character(LEN=30)` | Thermal sputtering law. Options: 'Tsai1998' (Tsai & Matthews 1998), others as implemented. |
| `accretion_model='Chaabouni2012'` |  `character(LEN=30)` | Accretion sticking efficiency model. Options: 'Chaabouni2012', others as implemented. |
| `shattering_model='Granato2021'` |  `character(LEN=30)` | Dispersion velocity model for shattering. Options: 'Granato2021', others as implemented. |
| `coagulation_model='Aoyama2017'` |  `character(LEN=30)` | Dispersion velocity model for coagulation. Options: 'Aoyama2017', others as implemented. |
| `dust_velocity_model='Ormel2007'` |  `character(LEN=30)` | Relative velocity model for grain-grain interactions. Options: 'Ormel2007', others as implemented. |
| `charging_model='Ibanez2019'` |  `character(LEN=30)` | Grain charge distribution model. Options: 'Ibanez2019', others as implemented. |
| `nZmix=3`                     |  `integer`    | Number of representative charge states (1: mean charge, 2: two-point, 3: three-point distribution). |

## PAH Modelling Options

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `sublimation_model='Galliano'` |  `character(LEN=30)` | Model for UV sublimation of PAHs. Options: 'Galliano', others as implemented. |
| `peh_attach_model='Berne'`    |  `character(LEN=30)` | Photoelectric heating model assumptions. Options: 'Berne', others as implemented. |
| `coalescence_model='Totton2012'` |  `character(LEN=30)` | PAH coalescence model. Options: 'Totton2012', others as implemented. |
| `pah_h2_model='RM2026'`       |  `character(LEN=30)` | Model for H₂ formation on PAHs. Options: 'RM2026', others as implemented. |
| `pah_growth_model='subgrid'`  |  `character(LEN=30)` | Model for PAH growth by accretion. Options: 'subgrid', others as implemented. |

## Timescale and Efficiency Parameters

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `t_sputter_ref=1d5`           |  `real`       | Sputtering reference timescale (in seconds). |
| `t_growth_ref=4d5`            |  `real`       | Accretion reference timescale (in seconds). |
| `t_sha_ref=7.573d5`           |  `real`       | Shattering reference timescale for large grains (0.1 μm, v=10 km/s, ρ=3 g/cm³) in seconds. |
| `t_coa_ref=2.71d5`            |  `real`       | Coagulation reference timescale for small grains (0.005 μm, v=0.1 km/s, ρ=3 g/cm³) in seconds. |
| `Sconstant=1.0`              |  `real`       | Sticking coefficient for grain collisions. |
| `nh_coa=0.1`                 |  `real array` | Gas density threshold above which dust coagulation is allowed (in H/cm³), per dust chemistry type. |
| `nhmax_acc=1d4`              |  `real array` | Maximum gas density for accretion subgrid model (in H/cm³), per dust chemistry type. |
| `nhmax_coa=1d6`              |  `real array` | Maximum gas density for coagulation subgrid model (in H/cm³), per dust chemistry type. |
| `nhmax_sha=1d3`              |  `real array` | Maximum gas density for shattering subgrid model (in H/cm³), per dust chemistry type. |
| `dust_SNdest_eff=0.1`        |  `real array` | Dust destruction efficiency in SN explosions (0–1), per dust chemistry type. |
| `dust_SNsha_eff=0.1`         |  `real array` | Shattering efficiency (fraction of mass converted to small grains), per dust chemistry type. |
| `dust_SNII_cond_eff=0.1`     |  `real array` | Core-collapse SN condensation efficiency (0–1), per dust chemistry type. |
| `dust_SNIa_cond_eff=0.1`     |  `real array` | Thermonuclear SN condensation efficiency (0–1), per dust chemistry type. |
| `dust_AGB_cond_eff=0.1`      |  `real array` | AGB wind condensation efficiency (0–1), per dust chemistry type. |
| `Coulomb_enhance=1.0`        |  `real array` | Enhancement factor for ion accretion due to Coulomb focusing, per dust chemistry type. |
| `tensile_strength=1d7`       |  `real array` | Dust tensile strength (erg/cm³), per dust chemistry type. |
| `Youngs_modulus=1d10`        |  `real array` | Dust Young's modulus (erg/cm³), per dust chemistry type. |
| `Poisson_ratio=0.25`         |  `real array` | Dust Poisson's ratio, per dust chemistry type. |
| `surf_energy=25.0`           |  `real array` | Dust surface energy (erg/cm²), per dust chemistry type. |
| `work_function=8.0`          |  `real array` | Dust work function (eV), per dust chemistry type. |
| `band_gap=5.0`               |  `real array` | Dust band gap (eV), per dust chemistry type. |
| `e_escape_length=1.0d-7`     |  `real array` | Dust electron escape length (cm), per dust chemistry type. |
| `separate_refractive_index=.false.` |  `logical array` | Whether dust bin has separate refractive indices for parallel/perpendicular waves, per dust chemistry type. |
| `slope_frag_func=0.433`      |  `real`       | Fragment size distribution power-law slope for shattering and RATD. |
| `errmax=0.1`                 |  `real`       | Convergence criterion for the dust chemistry solver. |
| `countmax=10000`             |  `real`       | Maximum number of iterations in the dust chemistry solver. |
| `GDinit=162.0`               |  `real`       | Initial gas-to-dust ratio (default: 162 from Zubko et al. 2004). |
| `DTMinit=-1.0`               |  `real`       | Initial dust-to-metal ratio. Default -1 means use solar value. |
| `fpah_ini=0.1`               |  `real`       | Initial fraction of carbon locked in PAHs. |
| `smallr_dust=1d-12`          |  `real`       | Minimum dust mass density threshold (g/cm³) for a cell to be treated as a dust cell. |

## Dust Grain and PAH Bin Properties

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `dust_composition=0.0`       |  `real array` | Mass fraction of each element in each dust chemistry type (n_elements × ndchemtype). Set from external files or hardcoded. |
| `dustbins_per_chemtype=0`    |  `integer array` | Number of dust bins in each dust chemistry type. |
| `asize=0.1`                  |  `real array` | Grain radius for each dust bin (in μm), length ndust. |
| `sgrain=1.0`                 |  `real array` | Grain material density divided by 3 g/cm³ for each dust bin, length ndust. |
| `amin=1d-2`                  |  `real array` | Minimum grain size in the underlying size distribution (in μm), per dust bin. |
| `amax=1.0`                   |  `real array` | Maximum grain size in the underlying size distribution (in μm), per dust bin. |
| `fmass_ej=0.0`               |  `real array` | Fraction of total dust mass in SN ejecta injected into each dust bin, length ndust. |
| `pah_nc=0`                   |  `integer array` | Number of carbon atoms in each PAH molecule, length npah. |
| `pah_nc_min=0`               |  `integer array` | Minimum number of carbon atoms in each PAH bin, length npah. |
| `pah_nc_max=0`               |  `integer array` | Maximum number of carbon atoms in each PAH bin, length npah. |
| `spah=2.0`                   |  `real array` | PAH material density (in g/cm³), length npah. Default: 2 g/cm³ for hydrocarbon. |
| `pah_SNdest_eff=0.1`         |  `real array` | PAH destruction efficiency in SN shocks (0–1), length npah. |
| `fpah_inwind=0.5`            |  `real array` | Fraction of AGB wind PAH mass in each PAH size bin, length npah. |
| `pah_ncharge_states=4`       |  `integer array` | Number of charge states for each PAH bin in charging calculations, length npah. |

## ISM Depletion Factors

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `fDust_depletions=[...]`     |  `real array` | Dust depletion factors for each element (fraction locked in dust), length n_elements. Defaults follow BARE-GR-S model from Zubko et al. (2004). |
| `fCDust_inPAH=0.1342`        |  `real`       | Fraction of dust carbon locked in PAHs. |
| `GD_solar=162.0`             |  `real`       | Reference gas-to-dust ratio in the solar neighborhood (Zubko et al. 2004). |
| `DTM_solar=0.458`            |  `real`       | Reference dust-to-metal ratio in the solar neighborhood (Zubko et al. 2004). |
| `fdustmass_ini=1/ndust`      |  `real array` | Initial dust mass fraction in each dust bin, length ndust. Default: equal in all bins. |
| `fpahmass_ini=1/npah`        |  `real array` | Initial PAH mass fraction in each PAH bin, length npah. Default: equal in all bins. |

## Radiation Parameters

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `fixed_rad_ani=-1.0`         |  `real`       | Fixed radiation field anisotropy parameter for RATD tests. Default -1 means no override. |
| `fixed_lambda_mean=-1.0`     |  `real`       | Fixed mean radiation wavelength for RATD tests (in μm). Default -1 means compute from radiation field. |

## External Dust Data

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `dust_tables_dir='../lib/dust_tables'` |  `character(LEN=256)` | Directory path containing dust and PAH optical properties files, cross-sections, charge tables, and other precomputed data. This can also be set by the environment variable `RAMSES_DUST_TABLES_DIR`. |

## Compilation Parameters

Note: The number of dust bins (`NDUST`) and PAH bins (`NPAH`) are compilation parameters, set in the Makefile with `-DNDUST=` and `-DNPAH=` flags. The number of dust chemistry types (`NDCHEMTYPE`, default 2 for carbonaceous and silicate) is also a compilation parameter.

## Notes

- All parameters in this namelist are optional and have sensible defaults.
- Parameters ending with `=.false.` are disabled by default; set to `.true.` to enable the corresponding physics.
- Array parameters (those marked as `real array` or `integer array`) are indexed per dust bin, PAH bin, or dust chemistry type as indicated.
- Dust composition and optical properties are loaded from external tables specified by `dust_tables_dir`.
- The dust and PAH models are modular; physics modules can be independently enabled or disabled for testing and performance studies.
- CALIMA must be enabled at compile time with `-DCALIMA` flag.
