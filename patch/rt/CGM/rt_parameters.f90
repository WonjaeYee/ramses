module rt_parameters
  use hydro_parameters

#ifdef NGROUPS
  integer,parameter::nGroups=NGROUPS          ! # of photon groups (set in Makefile)
#else
  integer,parameter::nGroups=1
#endif
  integer,parameter::nRTvar=nGroups*(1+ndim) ! # of RT variables (photon density and flux)

  real(dp)::rt_c=0., rt_c2=0.                ! RT constants in user units (set in init_rt)
  real(dp),parameter::c_cgs=2.9979250d+10                  ! Actual lightspeed in [cm s-1]
  real(dp),parameter::one_over_c_cgs=3.335640484668562d-11  ! save some computation
  real(dp)::rt_c_cgs=c_cgs                                     ! RT lightspeed in [cm s-1]
  real(dp),parameter::m_sun=1.9891d33               ! Solar mass [g], for SED calculations
  real(dp),parameter::eV_to_erg=1.6022d-12          !        eV to erg conversion constant
  real(dp), parameter:: Gyr2sec = 3.15569d+16       !       Gyr to sec conversion constant
  real(dp), parameter:: Myr2sec = 3.15569d+13       !       Myr to sec conversion constant
  real(dp), parameter:: sec2Gyr = 3.16888d-17       !       sec to Gyr conversion constant
  real(dp),parameter:: hp=6.6262d-27                !            Planck const   [erg sec ]
  real(rtdp),allocatable,dimension(:,:)::lambda1,lambda4                 ! HLL eigenvalues
#ifndef NRTPRE
  real(rtdp),parameter::smallNp=1d-50               !               Minimum photon density
#else
#if NRTPRE==4
  real(rtdp),parameter::smallNp=1e-30               !               Minimum photon density
#else
  real(rtdp),parameter::smallNp=1d-50               !               Minimum photon density
#endif
#endif
  ! Ion species---------------------------------------------------------------------------
#ifdef NIONS
  integer,parameter::nIons=NIONS                     ! # of ion species (set in Makefile)
#else
  integer,parameter::nIons=3                         !                     Hii, Heii, Heiii
#endif
  integer::iIons=7                                  !    Starting index of ion states in U
  ! Ionization energies
  real(dp),dimension(nIons)::ionEvs
  real(dp):: ionEv_HI = 12.27
  real(dp),parameter::ionEv_HII   = 13.60
  real(dp),parameter::ionEv_HeII  = 24.59
  real(dp),parameter::ionEv_HeIII = 54.42
  logical::isHe=.false.                             ! He ionization fractions tracked?
  logical::isH2=.false.                             !                      H2 tracked?
  logical::isCO=.false.                             !                      CO tracked?
  integer::ixHI=0, ixHII=0, ixHeII=0, ixHeIII=0     !  Indices of ionization fractions

  ! RT_PARAMS namelist--------------------------------------------------------------------
  logical::rt_advect=.false.           ! Advection of photons?                           !
  logical::rt_smooth=.false.           ! Smooth the discrete RT update of op. splitting  !
  real(dp)::rt_Tconst=-1               ! If pos. use this value for all T-depend. rates  !
  logical::rt_isTconst=.false.         ! Const rates activated?                          !
  logical::rt_star=.false.             ! Activate radiation from star particles?         !
  real(dp)::rt_esc_frac=1.d0           ! Escape fraction of light from stellar particles !
  logical::rt_is_init_xion=.false.     ! Initialize ionization from T profile?           !
  character(LEN=10)::rt_flux_scheme='glf'                                                !
  logical::rt_use_hll=.false.          ! Use hll flux (or the default glf)               !
  logical::rt_is_outflow_bound=.false. ! Make all boundaries=outflow for RT              !
  real(dp)::rt_courant_factor=0.8d0    ! Courant factor for RT timesteps                 !
  logical::rt_refine=.false.           ! Refine on RT-related conditions?                !
  real(dp)::rt_err_grad_n=-1.0         ! Photon number density gradient for refinement   !
  real(dp)::rt_floor_n=1.d-10          ! Photon number density floor for refinement      !
  real(dp)::rt_err_grad_xHI=-1.0       ! Ionization state gradient for refinement        !
  real(dp)::rt_err_grad_xHII=-1.0      ! Ionization state gradient for refinement        !
  real(dp)::rt_refine_aexp=-1.0        ! Start a for RT gradient refinement              !
  real(dp)::rt_floor_xHI=1.d-10        ! Ionization state floor for refinement           !
  real(dp)::rt_floor_xHII=1.d-10       ! Ionization state floor for refinement           !
  real(dp)::rt_c_fraction=1.d0         ! Actual lightspeed fraction for RT lightspeed    !
  integer::rt_nsubcycle=1                ! Maximum number of RT-steps during one hydro/    !
                                       ! gravity/etc timestep                            !
  logical::rt_otsa=.true.              ! Use on-the-spot approximation                   !
  logical::rt_isDiffuseUVsrc=.false.   ! UV emission from low-density cells              !
  real(dp)::rt_UVsrc_nHmax=-1d0        ! Density threshold for UV emission               !
  logical::upload_equilibrium_x=.false.! Enforce equilibrium xion when uploading         !
!  logical::convert_birth_times=.false. ! Convert stellar birthtimes: conformal -> proper !

  character(LEN=128)::hll_evals_file=''! File HLL eigenvalues                            !
  character(LEN=128)::sed_dir=''       ! Dir containing stellar energy distributions     !
  character(LEN=128)::uv_file=''       ! File containing stellar energy distributions    !
  character(LEN=128)::uvbg_rtz_dir=''       ! File containing stellar energy distributions    !

  ! RT_GROUPS namelist--------------------------------------------------------------------
  integer::sedprops_update=-1                      ! Update sedprops from star populations
  ! negative: never update, 0:update on init, pos x: update every x coarse steps
  logical::SED_isEgy=.false. ! Integrate energy out of SEDs rather than photon count
  ! Grop props: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev).
  ! Indexes nGroups, nIons stand for photon group vs species (e.g. 1=H, 2=He).
  integer,dimension(nGroups)::iGroups=1                          ! Start indices of groups
  real(dp),dimension(nGroups,nIons)::group_csn=0, group_cse=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups)::group_egy=0                       !  Avg photon energy (ev)
  real(dp),dimension(nGroups)::groupL0=13.60                     ! Wavelength lower limits
  real(dp),dimension(nGroups)::groupL1=0                         ! Wavelength upper limits
  integer,dimension(nIons)::spec2group=0                 !Ion -> group # in recombinations

  ! Imposed boundary condition variables
  real(rtdp),dimension(1:MAXBOUND,1:nrtvar)::rt_boundary_var
  real(dp),dimension(1:MAXBOUND)::rt_n_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::rt_u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::rt_v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::rt_w_bound=0.0d0

  ! Initial condition RT regions parameters----------------------------------------------
  integer                           ::rt_nregion=0
  character(LEN=10),dimension(1:MAXREGION)::rt_region_type='square'
  real(dp),dimension(1:MAXREGION)   ::rt_reg_x_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_y_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_z_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_exp_region=2.0
  integer,dimension(1:MAXREGION)    ::rt_reg_group=1
  real(rtdp),dimension(1:MAXREGION)   ::rt_n_region=0.                     ! Photon density
  real(rtdp),dimension(1:MAXREGION)   ::rt_u_region=0.                     ! Photon flux
  real(rtdp),dimension(1:MAXREGION)   ::rt_v_region=0.                     ! Photon flux
  real(rtdp),dimension(1:MAXREGION)   ::rt_w_region=0.                     ! Photon flux
  real(rtdp),dimension(1:MAXREGION)   ::rt_xion_region=0.  ! Xion state (square regions only)

   ! RT source regions parameters----------------------------------------------------------
  integer                           ::rt_nsource=0
  character(LEN=10),dimension(1:MAXREGION)::rt_source_type='square'
  real(dp),dimension(1:MAXREGION)   ::rt_src_x_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_y_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_z_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_exp_source=2.0
  integer, dimension(1:MAXREGION)   ::rt_src_group=1
  real(dp),dimension(1:MAXREGION)   ::rt_n_source=0.                      ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_source=0.                         ! Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_wind_source=0.                     ! Stellar wind

  ! Indexing in flux_module
  integer,parameter::ifrt1=0                                                           ! 0
  integer,parameter::jfrt1=1-ndim/2                                               ! 0 or 1
  integer,parameter::kfrt1=1-ndim/3                                               ! 0 or 1


  ! Cooling statistics: avg loop # per cell, maximum loop #, # of cooling calls-----------
  logical::rt_output_coolstats=.false.    ! Output cooling statistics                     !
  integer*8::tot_cool_loopcnt=0,max_cool_loopcnt=0,n_cool_cells=0
  integer*8,dimension(20)::loopCodes=0

  ! SED statistics: Radiation emitted, total, last coarse step [#photons/10^50]-----------
  logical::showSEDstats=.true.
  real(dp)::tot_nPhot, step_nPhot, step_nStar, step_mStar

  logical::inLastCoarseStep=.false.    ! .t. when doing last ilevel step in coarse step  !
  logical::doDump = .false.

  ! Harley Katz additions 
  integer::LWgroup  = -1            ! index for the Lyman-Werner group
  logical::isH2Katz = .false.       ! more sophisticated H2 tracked?
  logical::KatzSED  = .false.       ! use Harley's SED file
  logical::SS_LVG   = .false.       ! use the LVG approximation to estimate the scale of self-shielding
  logical::PEH      = .true.        ! photoelectric heating
  logical::PHlocal  = .true.        ! local photoionisation heating (Haardt-Madau UV is not relevant)
  real(dp)::T_sputter = 1d20 ! dust is assumed to be sputtered below this temperature
  logical::fire_cooling = .false.   ! Whether to use fire low temp heating and cooling
  logical::low_t_noneq_cooling = .false. ! Whether to use low temperature non-equilibrium cooling
  logical::include_collisional_ionization = .true. ! Whether to include collisional ionization
  logical::include_photoionisation = .true. ! Whether to include photo-ionization
  logical::include_charge_transfer = .true. ! Whether to include charge transfer reactions
  logical::include_cosmic_ray_ionization = .true. ! Whether to include cosmic ray ionization
  logical::include_cosmic_ray_heating = .true. ! Whether to include cosmic ray heating
  logical::include_grain_recombination = .false. ! Whether to include recombination on dust grains
  logical::include_dust_heat_cool = .true. ! Whether to include dust heating/cooling
  logical::include_ct_heat_cool = .true. ! Whether to include charge exchange heating and cooling
  real(dp)::xi_h_cr = 2.0d-16 ! [s^-1 H^-1]
  logical::use_cloudy_prim_rates = .false. ! Whether to use the H and He rates from CLOUDY instead of ramses-rt

  real(dp),parameter,dimension(8)::oxygen_ionEvs=(/ 1.362d1, 3.512d1, 5.494d1, 7.741d1, 1.139d2, 1.381d2, 7.393d2, 8.714d2 /)
  real(dp),parameter,dimension(7)::nitrogen_ionEvs=(/ 1.453d1, 2.960d1, 4.745d1, 7.747d1, 9.789d1, 5.521d2, 6.671d2 /)
  real(dp),parameter,dimension(6)::carbon_ionEvs=(/ 1.126d1, 2.438d1, 4.789d1, 6.449d1, 3.921d2, 4.900d2 /)
  real(dp),parameter,dimension(12)::magnesium_ionEvs=(/ 7.646d0, 1.504d1, 8.014d1, 1.093d2, 1.413d2, 1.865d2, 2.249d2, 2.660d2, &
                                                        3.282d2, 3.675d2, 1.762d3, 1.963d3 /)
  real(dp),parameter,dimension(14)::silicon_ionEvs=(/ 8.152d0, 1.635d1, 3.349d1, 4.514d1, 1.668d2, 2.051d2, 2.465d2, 3.032d2, &
                                                      3.511d2, 4.014d2, 4.761d2, 5.235d2, 2.438d3, 2.673d3 /)
  real(dp),parameter,dimension(16)::sulfur_ionEvs=(/ 10.36d0, 23.33d0, 34.83d0, 47.31d0, 72.68d0, &
                                                     88.05d0, 280.9d0, 328.2d0, 379.1d0, 447.1d0, &  
                                                     504.8d0, 564.7d0, 651.7d0, 707.2d0, 3224.d0, &
                                                     3494.d0 /)
  real(dp),parameter,dimension(26)::iron_ionEvs=(/ 7.902d+00, 1.619d+01, 3.065d+01, 5.480d+01, &
                                                   7.501d+01, 9.906d+01, 1.250d+02, 1.511d+02, &
                                                   2.336d+02, 2.621d+02, 2.902d+02, 3.308d+02, &
                                                   3.610d+02, 3.922d+02, 4.570d+02, 4.893d+02, &
                                                   1.262d+03, 1.358d+03, 1.456d+03, 1.582d+03, &
                                                   1.689d+03, 1.799d+03, 1.950d+03, 2.046d+03, &
                                                   8.829d+03, 9.278d+03 /)
  real(dp),parameter,dimension(10)::neon_ionEvs=(/ 2.156d+01, 4.096d+01, 6.346d+01, 9.712d+01, 1.262d+02, &
                                                   1.579d+02, 2.073d+02, 2.391d+02, 1.196d+03, 1.362d+03 /)                                                 
  real(dp),dimension(nGroups,n_oxygen_ions)::group_csn_oxygen=0, group_cse_oxygen=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_nitrogen_ions)::group_csn_nitrogen=0, group_cse_nitrogen=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_carbon_ions)::group_csn_carbon=0, group_cse_carbon=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_magnesium_ions)::group_csn_magnesium=0, group_cse_magnesium=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_silicon_ions)::group_csn_silicon=0, group_cse_silicon=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_sulfur_ions)::group_csn_sulfur=0, group_cse_sulfur=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_iron_ions)::group_csn_iron=0, group_cse_iron=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,n_neon_ions)::group_csn_neon=0, group_cse_neon=0    !    Cross sections (cm2)
  real(dp),dimension(nGroups,2)::group_csn_dust=0   !    Cross sections (cm2)

  ! For a UV background
  logical::uvbg_rtz=.false.
  real(dp)::uvbg_suppression_redshift=-1.d0
  logical::no_subionizing_ubvg=.false.
  real(dp)::UV_background_hydrogen = 0.d0, UV_background_H2 = 0.d0 
  real(dp)::UV_background_G0 = 0.d0 
  real(dp),dimension(2)::UV_background_helium = 0.d0
  real(dp),dimension(n_oxygen_ions)::UV_background_oxygen = 0.d0
  real(dp),dimension(n_nitrogen_ions)::UV_background_nitrogen = 0.d0
  real(dp),dimension(n_carbon_ions)::UV_background_carbon = 0.d0
  real(dp),dimension(n_magnesium_ions)::UV_background_magnesium = 0.d0
  real(dp),dimension(n_silicon_ions)::UV_background_silicon = 0.d0
  real(dp),dimension(n_sulfur_ions)::UV_background_sulfur = 0.d0
  real(dp),dimension(n_iron_ions)::UV_background_iron = 0.d0
  real(dp),dimension(n_neon_ions)::UV_background_neon = 0.d0
  real(dp)::UV_background_hydrogen_heating = 0.d0, UV_background_H2_heating = 0.d0 
  real(dp),dimension(2)::UV_background_helium_heating = 0.d0
  real(dp),dimension(n_oxygen_ions)::UV_background_oxygen_heating = 0.d0
  real(dp),dimension(n_nitrogen_ions)::UV_background_nitrogen_heating = 0.d0
  real(dp),dimension(n_carbon_ions)::UV_background_carbon_heating = 0.d0
  real(dp),dimension(n_magnesium_ions)::UV_background_magnesium_heating = 0.d0
  real(dp),dimension(n_silicon_ions)::UV_background_silicon_heating = 0.d0
  real(dp),dimension(n_sulfur_ions)::UV_background_sulfur_heating = 0.d0
  real(dp),dimension(n_iron_ions)::UV_background_iron_heating = 0.d0
  real(dp),dimension(n_neon_ions)::UV_background_neon_heating = 0.d0


  ! Lyman alpha pressure
  integer       :: iLyaVar=0        ! NLya index for uold
  logical       :: rt_Lya_pressure = .false.
  logical       :: rt_Lya_pressure_LVG=.true.     ! use the Large velocity gradient approximation
  real(kind=dp) :: rt_Lya_pressure_nH = 0.1       ! minimum nH for Lya pressure [H/cm-3]
  logical       :: rt_Lya_pressure_count_H2 = .true.  ! include molecular hydrogen when measuring N_HI
  real(dp)      :: rt_T_sputter = 1d5     ! dust sputtering temperature
  real(kind=dp) :: sigdust21  = 3.0       ! dust absorption cross-section in units of 10^-21 cm2/H at solar metallicity
  real(kind=dp) :: DustAlbedo = 0.46      ! consistent with RASCAS
  logical :: dust2metal_RR14 = .false.    ! metallicity-dependent metal-to-dust ratio based on Remy-Ruyer+(14)
  real(kind=dp) :: H2clumping=10.0        ! clumping factor for H2 (default: Gnedin et al. 2009)

  ! Prevent cooling catastrophe
  ! if true, exponential decrease of metal cooling at low density for stability
  logical::metal_cooling_catastrophe_fix=.true.
  real(dp)::metal_cooling_catastrophe_rho=0.005d0

  ! Whether to use a model for dust depletion to suppress various cooling lines
  logical::simple_dust_depletion=.true.

end module rt_parameters
