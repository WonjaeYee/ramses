!*************************************************************************
SUBROUTINE rt_init

!  Initialize everything for radiative transfer
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_hydro_commons
  use rt_flux_module
  use rt_cooling_module,only:rt_isIRtrap,iIRtrapVar
  use rt_parameters
  use SED_module
  use UV_module
  implicit none
  integer:: i, ilevel, ivar, nvar_count
!-------------------------------------------------------------------------
  if(verbose)write(*,*)'Entering init_rt'
  ! Count the number of variables and check if ok:
  nvar_count = ichem-1     ! # of non-rt vars: rho u v w p (z) (delay) (x)
  if(rt_isIRtrap) &
     iIRtrapVar = ndim+3  ! Trapped rad. stored in nonthermal pressure var
  iIons=nvar_count+1         !      Starting index of ionisation fractions
  if(myid==1) write(*,*) '>>>NOTE: iIons    =',iIons
  nvar_count = iIons+NIONS-1 !                                # hydro vars
  if(rt_Lya_pressure) then
     iLyaVar=nvar_count+1
     nvar_count = nvar_count+1
  endif
  if(nvar_count .gt. nvar) then 
     if(myid==1) then 
        write(*,*) 'rt_init(): Something wrong with NVAR.'
        write(*,*) 'Should have NVAR=2+ndim+1*dcool+1*aton+IRtrap+nIons+nmetals'
        write(*,*) 'Have NVAR=',nvar
        write(*,*) 'Should have NVAR=',nvar_count
        write(*,*) 'STOPPING!'
     endif
     call clean_stop
  endif

  if(rt_star .or. sedprops_update .ge. 0) &
     call init_SED_table    ! init stellar energy distribution properties

  if(rt .and. .not. hydro) then
     if(myid==1) then
        write(*,*) 'hydro must be turned on when running radiative transfer.'
        write(*,*) 'STOPPING!'
     endif
     call clean_stop
  endif
  if(rt_star) use_proper_time=.true.    ! Need proper birth time for stars
  if(rt) neq_chem=.true.        ! Equilibrium cooling doesn't work with RT
  
  ! To maximize efficiency, rt advection and rt timestepping is turned off
  ! until needed.
  if(rt .and. .not.rt_otsa) rt_advect=.true.                              
  if(rt .and. rt_nsource .gt. 0) rt_advect=.true.                         
  if(rt .and. rt_nregion .gt. 0) rt_advect=.true.
  ! UV propagation is checked in set_model
  ! Star feedback is checked in amr_step

  ! Update hydro variable to the initial ionized species
  var_region(1:rt_nregion,iIons-ndim-2)=rt_xion_region(1:rt_nregion)
  do i=1,nGroups  ! Starting indices in uold and unew of each photon group
     iGroups(i)=1+(ndim+1)*(i-1)
     if(nrestart.eq.0) then
        rtuold(:,iGroups(i))=smallNp
     endif
  end do
  if(trim(rt_flux_scheme).eq.'hll') rt_use_hll=.true.
  if(rt_use_hll) call read_hll_eigenvalues

  tot_cool_loopcnt=0 ; max_cool_loopcnt=0 ; n_cool_cells=0
  loopCodes=0
  tot_nPhot=0.d0 ;  step_nPhot=0.d0; step_nStar=0.d0; step_mStar=0.d0
END SUBROUTINE rt_init

!*************************************************************************
SUBROUTINE update_rt_c

! Update the speed of light for radiative transfer, in code units.
! This cannot be just a constant, since scale_v changes with time in 
! cosmological simulations.
!-------------------------------------------------------------------------
  use rt_parameters
  use amr_commons
  implicit none
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c=rt_c_cgs/scale_v
  rt_c2=rt_c**2
END SUBROUTINE update_rt_c

!*************************************************************************
SUBROUTINE adaptive_rt_c_update(ilevel, dt)

! Set the lightspeed such that RT can be done at ilevel in time dt in 
! a single step.
!-------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  use SED_module
  implicit none
  integer:: ilevel, nx_loc
  real(dp):: dt, scale, dx
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!-------------------------------------------------------------------------
  ! Mesh spacing at ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! new lightspeed
  rt_c = dx/3.d0/dt * rt_courant_factor 
  rt_c2 = rt_c**2

  ! new ligtspeed in cgs
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c_cgs = rt_c*scale_v
  rt_c_fraction = rt_c_cgs/c_cgs

  call updateRTGroups_CoolConstants        ! These change as a consequence

END SUBROUTINE adaptive_rt_c_update


!*************************************************************************
SUBROUTINE read_rt_params(nml_ok)

! Read rt_params namelist
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use cooling_module, only:X, Y
  use rt_cooling_module
  use UV_module
  use SED_module
  implicit none
  logical::nml_ok
  integer::iCount
!-------------------------------------------------------------------------
  namelist/rt_params/rt_star, rt_esc_frac, rt_flux_scheme, rt_smooth     &
       & ,rt_is_outflow_bound, rt_TConst, rt_courant_factor              &
       & ,rt_c_fraction, rt_nsubcycle, rt_otsa, sedprops_update          &
       & ,sed_dir, uv_file, rt_UVsrc_nHmax, nUVgroups, nSEDgroups        &
       & ,SED_isEgy, rt_output_coolstats, hll_evals_file                 &
       & ,upload_equilibrium_x, X, Y, rt_is_init_xion                    &
       & ,rt_err_grad_n, rt_floor_n, rt_err_grad_xHII, rt_floor_xHII     &
       & ,rt_err_grad_xHI, rt_floor_xHI, rt_refine_aexp                  &
       & ,convert_birth_times,isHe,isH2                                  &
       & ,rt_isIR, is_kIR_T, rt_T_rad, rt_vc, rt_pressBoost              &
       & ,rt_isoPress, rt_isIRtrap, fire_cooling                         &
       ! RT regions (for initialization)                                 &
       & ,rt_nregion, rt_region_type                                     &
       & ,rt_reg_x_center, rt_reg_y_center, rt_reg_z_center              &
       & ,rt_reg_length_x, rt_reg_length_y, rt_reg_length_z              &
       & ,rt_exp_region, rt_reg_group                                    &
       & ,rt_n_region, rt_u_region, rt_v_region, rt_w_region             &
       & ,rt_xion_region                                                 &
       ! RT source regions (for every timestep)                          &
       & ,rt_nsource, rt_source_type                                     &
       & ,rt_src_x_center, rt_src_y_center, rt_src_z_center              &
       & ,rt_src_length_x, rt_src_length_y, rt_src_length_z              &
       & ,rt_exp_source, rt_src_group                                    &
       & ,rt_n_source, rt_u_source, rt_v_source, rt_w_source             &
       ! RT boundary (for boundary conditions)                           &
       & ,rt_n_bound,rt_u_bound,rt_v_bound,rt_w_bound                    &
       & ,rt_movie_vars,el_movie_vars                                    &
       ! RT additions by Harley Katz                                     &
       & ,LWgroup,isH2Katz,KatzSED,SS_LVG,T_sputter,PEH,PHlocal          & 
       & ,rt_Lya_pressure,rt_Lya_pressure_nH,rt_Lya_pressure_LVG         &
       & ,rt_T_sputter,rt_Lya_pressure_count_H2,sigdust21                &
       & ,dust2metal_RR14,H2clumping, include_collisional_ionization     &
       & ,include_photoionisation, UV_background_hydrogen                &
       & ,UV_background_helium, UV_background_oxygen                     &
       & ,UV_background_nitrogen, UV_background_carbon                   &
       & ,UV_background_magnesium, UV_background_silicon                 &
       & ,UV_background_sulfur, UV_background_iron, UV_background_neon   &
       & ,UV_background_H2, include_charge_transfer, low_t_noneq_cooling &
       & ,use_cloudy_prim_rates                                          &
       & ,UV_background_hydrogen_heating, UV_background_H2_heating       &
       & ,UV_background_helium_heating, UV_background_oxygen_heating     &
       & ,UV_background_nitrogen_heating, UV_background_carbon_heating   &
       & ,UV_background_magnesium_heating, UV_background_silicon_heating &
       & ,UV_background_sulfur_heating, UV_background_iron_heating       &
       & ,UV_background_neon_heating, metal_cooling_catastrophe_fix      &
       & ,metal_cooling_catastrophe_rho, simple_dust_depletion 

  if(nrestart .gt. 0) then
     rt_is_init_xion=.false.
  else
     rt_is_init_xion=.true.
  endif

  ! Read namelist file
  rewind(1)
  read(1,NML=rt_params,END=101)
101 continue                                   ! No harm if no rt namelist

  if(nGroups.le.0) rt=.false. ! No sense  doing rt if there are no photons
  if(.not. rt .and. .not. rt_star) sedprops_update=-1

  if(rt_err_grad_n .gt. 0. .or. rt_err_grad_xHII .gt. 0.                 &
       .or. rt_err_grad_xHI .gt. 0.) rt_refine=.true.

  rt_c_cgs = c_cgs * rt_c_fraction

  ! Trapped IR pressure closure as in Rosdahl & Teyssier 2015, eq 43:
  if(rt_isIRtrap) gamma_rad(1) = rt_c_fraction / 3d0 + 1d0 

  !call update_rt_c
  if(rt_Tconst .ge. 0.d0) rt_isTconst=.true. 

  ! Set indexes of ionization fractions, and ionization energies, and
  ! check if we have enough ionization variables (NIONS)
  iCount=0
  if(isH2) then
     iCount=iCount+1 ; ixHI=iCount; ionEvs(ixHI)=ionEv_HI
     ! two different H2 dissociation channels for Harley's approach 
     if(isH2Katz) ionEvs(ixHI)=15.20
  endif
  iCount=iCount+1    ; ixHII=iCount   ; ionEvs(ixHII)=ionEv_HII
  if(isHe) then
     iCount=iCount+1 ; ixHeII=iCount  ; ionEvs(ixHeII)=ionEv_HeII
     iCount=iCount+1 ; ixHeIII=iCount ; ionEvs(ixHeIII)=ionEv_HeIII
  endif
  if(iCount .gt. NIONS) then
     if(myid==1) then
        write(*,*) 'Not enough variables for ionization fractions'
        write(*,*) 'Have NIONS=',NIONS
        write(*,*) 'Need NIONS=',iCount
        write(*,*) 'STOPPING!'
     endif
     call clean_stop
  endif
  if(iCount .lt. NIONS) then
     if(myid==1) then
        write(*,*) 'Too many variables for ionization fractions'
        write(*,*) 'Have NIONS=',NIONS
        write(*,*) 'Need NIONS=',iCount
        write(*,*) 'Probably no harm, so still continuing...'
     endif
  endif
  if(myid==1) then
     write(*,*) 'Number of ionization fractions is:',iCount
     write(*,*) 'The indexes are iHI, iHII, iHeII, iHeIII ='              &   
                , ixHI, ixHII, ixHeII, ixHeIII
  endif

  call read_rt_groups(nml_ok)
END SUBROUTINE read_rt_params

!*************************************************************************
SUBROUTINE read_rt_groups(nml_ok)

! Read rt_groups namelist
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use rt_cooling_module
  use SED_module
  implicit none
  logical::nml_ok
  integer::i,igroup_HI=0, igroup_HII=0, igroup_HeII=0, igroup_HeIII=0
!-------------------------------------------------------------------------
  namelist/rt_groups/group_csn, group_cse, group_egy, spec2group         &
       & , groupL0, groupL1, kappaAbs, kappaSc
  if(myid==1) then
     write(*,'(" Working with ",I2," photon groups and  "                &
          ,I2, " ion species")') nGroups, nIons
     write(*,*) ''
  endif
   
  if(nGroups .le. 0) then
     rt = .false.
     return
  endif

  ! BEGIN joki ===========================================================
  !  Use H2, HI, HeI, HeII ionization energies  as default group intervals
  groupL0(1:min(nGroups,nIons))=ionEvs(1:min(nGroups,nIons))! Lower bounds
  groupL1(1:min(nGroups,nIons-1))=ionEvs(2:min(nGroups+1,nIons)) !   Upper
  groupL1(min(nGroups,nIons))=0.                    ! Upper bound=infinity

  i=0 
  if(isH2) then ! Set index for H2 dissociating group
     i=i+1 ; igroup_HI=i
  endif
  if(i .lt. nGroups) then ! Set index for HI ionizing group
     i=i+1 ; igroup_HII=i
  endif
  if(i .lt. nGroups .and. isHe) then ! Set index for HeI ionizing group
     i=i+1 ; igroup_HeII=i
  endif
  if(i .lt. nGroups .and. isHe) then ! Set index for HeII ionizing group
     i=i+1 ; igroup_HeIII=i
  endif
  
  ! Default groups are all blackbodies at E5 Kelvin:
  group_csn=0d0 ; group_cse=0d0 ; group_egy=0d0         ! Default all zero
  if(igroup_HI .gt. 0) then
     group_csn(igroup_HI,ixHI)=3d-22                     ! H2 dissociation
     group_cse(igroup_HI,ixHI)=3d-22
     group_egy(igroup_HI)=12.44
  endif
  if(igroup_HII .gt. 0) then
     group_csn(igroup_HII,ixHII)=3.007d-18                 ! HI ionization
     group_cse(igroup_HII,ixHII)=2.781d-18
     group_egy(igroup_HII)=18.85
  endif
  if(igroup_HeII .gt. 0) then
     group_csn(igroup_HeII,ixHII)=5.687d-19 ! HI ionization by HeI photons
     group_cse(igroup_HeII,ixHII)=5.042d-19
     group_csn(igroup_HeII,ixHeII)=4.478d-18              ! HeI ionization
     group_cse(igroup_HeII,ixHeII)=4.130d-18
     group_egy(igroup_HeII)=35.079
  endif
  if(igroup_HeIII .gt. 0) then
     group_csn(igroup_HeIII,ixHII)=7.889d-20   ! HI ioniz. by HeII photons 
     group_cse(igroup_HeIII,ixHII)=7.456d-20
     group_csn(igroup_HeIII,ixHeII)=1.197d-18 ! HeI ioniz. by HeII photons
     group_cse(igroup_HeIII,ixHeII)=1.142d-18
     group_csn(igroup_HeIII,ixHeIII)=1.055d-18           ! HeII ionization
     group_cse(igroup_HeIII,ixHeIII)=1.001d-18
     group_egy(igroup_HeIII)=65.666
  endif 

  ! Set the cross sections for each metal
  ! The default here is to assume that we have the standard 8 bins
  ! used by harley and taysun.  this will need to be changed otherwise
  if (oxygen_ions) then
   group_csn_oxygen(:,:) = 0.d0
   group_cse_oxygen(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NOXYGENIONS>0
      group_csn_oxygen(:,1) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 3.250E-18, 9.351E-18, 9.915E-18, 2.016E-18 /)
      group_cse_oxygen(:,1) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 3.281E-18, 9.634E-18, 9.494E-18, 1.206E-18 /)
#endif
#if NOXYGENIONS>1
      group_csn_oxygen(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.941E-18, 2.079E-18 /)
      group_cse_oxygen(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.653E-18, 1.233E-18 /)
#endif
#if NOXYGENIONS>2
      group_csn_oxygen(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.807E-18 /)
      group_cse_oxygen(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.101E-18 /)
#endif
#if NOXYGENIONS>3
      group_csn_oxygen(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.766E-19 /)
      group_cse_oxygen(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.916E-19 /)
#endif
#if NOXYGENIONS>4
      group_csn_oxygen(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.668E-19 /)
      group_cse_oxygen(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.975E-19 /)
#endif
#if NOXYGENIONS>5
      group_csn_oxygen(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.344E-20 /)
      group_cse_oxygen(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.490E-20 /)
#endif
   endif
  endif

  if (nitrogen_ions) then
   group_csn_nitrogen(:,:) = 0.d0
   group_cse_nitrogen(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NNITROGENIONS>0
      group_csn_nitrogen(:,1) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 5.110E-18, 1.381E-17, 8.621E-18, 1.241E-18 /)
      group_cse_nitrogen(:,1) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 5.284E-18, 1.377E-17, 8.070E-18, 7.278E-19 /)
#endif
#if NNITROGENIONS>1
      group_csn_nitrogen(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.856E-18, 1.160E-18 /)
      group_cse_nitrogen(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.181E-18, 6.900E-19 /)
#endif
#if NNITROGENIONS>2
      group_csn_nitrogen(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.842E-19, 1.018E-18 /)
      group_cse_nitrogen(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.778E-19, 6.179E-19 /)
#endif
#if NNITROGENIONS>3
      group_csn_nitrogen(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.332E-19 /)
      group_cse_nitrogen(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.817E-19 /)
#endif
#if NNITROGENIONS>4
      group_csn_nitrogen(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.064E-19 /)
      group_cse_nitrogen(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.096E-19 /)
#endif
   endif
  endif

  if (carbon_ions) then
   group_csn_carbon(:,:) = 0.d0
   group_cse_carbon(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NCARBONIONS>0
      group_csn_carbon(:,1) = (/ 0.E0, 0.000E+00, 0.000E+00, 1.647E-17, 1.514E-17, 1.111E-17, 4.874E-18, 6.540E-19 /)
      group_cse_carbon(:,1) = (/ 0.E0, 0.000E+00, 0.000E+00, 1.652E-17, 1.512E-17, 1.084E-17, 4.522E-18, 3.822E-19 /)
#endif
#if NCARBONIONS>1
      group_csn_carbon(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.662E-20, 3.785E-18, 6.110E-19 /)
      group_cse_carbon(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.318E-20, 3.559E-18, 3.640E-19 /)
#endif
#if NCARBONIONS>2
      group_csn_carbon(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.140E-19, 5.382E-19 /)
      group_cse_carbon(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.716E-19, 3.280E-19 /)
#endif
#if NCARBONIONS>3
      group_csn_carbon(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.190E-19 /)
      group_cse_carbon(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.536E-19 /)
#endif
#if NCARBONIONS>4
      group_csn_carbon(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.153E-20 /)
      group_cse_carbon(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.671E-20 /)
#endif
#if NCARBONIONS>5
      group_csn_carbon(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.159E-22 /)
      group_cse_carbon(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.501E-21 /)
#endif
   endif
  endif

  if(magnesium_ions) then
   group_csn_magnesium(:,:) = 0.d0
   group_cse_magnesium(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NMAGNESIUMIONS>0
      group_csn_magnesium(:,1) = (/ 0.E0, 0.000E+00, 2.555E-19, 1.920E-19, 1.576E-19, 2.178E-19, 2.283E-19, 1.231E-21 /)
      group_cse_magnesium(:,1) = (/ 0.E0, 0.000E+00, 3.047E-19, 1.909E-19, 1.578E-19, 2.217E-19, 2.202E-19, 4.905E-22 /)
#endif
#if NMAGNESIUMIONS>1
      group_csn_magnesium(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 2.168E-20, 2.487E-19, 1.743E-19, 1.434E-20 /)
      group_cse_magnesium(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 2.282E-20, 2.479E-19, 1.660E-19, 6.217E-21 /)
#endif
#if NMAGNESIUMIONS>2
      group_csn_magnesium(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.459E-18 /)
      group_cse_magnesium(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.226E-18 /)
#endif
#if NMAGNESIUMIONS>3
      group_csn_magnesium(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.196E-18 /)
      group_cse_magnesium(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.370E-18 /)
#endif
#if NMAGNESIUMIONS>4
      group_csn_magnesium(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.439E-19 /)
      group_cse_magnesium(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.693E-19 /)
#endif
#if NMAGNESIUMIONS>5
      group_csn_magnesium(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.978E-19 /)
      group_cse_magnesium(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.567E-19 /)
#endif
#if NMAGNESIUMIONS>6
      group_csn_magnesium(:,7) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 8.746E-20 /)
      group_cse_magnesium(:,7) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.857E-19 /)
#endif
#if NMAGNESIUMIONS>7
      group_csn_magnesium(:,8) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.695E-20 /)
      group_cse_magnesium(:,8) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 9.176E-20 /)
#endif
#if NMAGNESIUMIONS>8
      group_csn_magnesium(:,9) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.307E-20 /)
      group_cse_magnesium(:,8) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.723E-20 /)
#endif
#if NMAGNESIUMIONS>9
      group_csn_magnesium(:,10) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.216E-21 /)
      group_cse_magnesium(:,10) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.296E-20 /)
#endif
   endif
  endif

  if (silicon_ions) then
   group_csn_silicon(:,:) = 0.d0
   group_cse_silicon(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NSILICONIONS>0
      group_csn_silicon(:,1) = (/ 0.E0, 0.000E+00, 9.723E-18, 2.077E-17, 8.998E-18, 2.734E-18, 6.322E-19, 2.490E-19 /) 
      group_cse_silicon(:,1) = (/ 0.E0, 0.000E+00, 1.239E-17, 2.060E-17, 8.954E-18, 2.503E-18, 6.437E-19, 1.307E-19 /)
#endif
#if NSILICONIONS>1
      group_csn_silicon(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.944E-19, 5.398E-19, 2.219E-19 /)
      group_cse_silicon(:,2) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.786E-19, 5.472E-19, 1.218E-19 /)
#endif
#if NSILICONIONS>2
      group_csn_silicon(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.745E-19, 1.609E-19 /)
      group_cse_silicon(:,3) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.122E-19, 9.083E-20 /)
#endif
#if NSILICONIONS>3
      group_csn_silicon(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.267E-20, 1.041E-19 /)
      group_cse_silicon(:,4) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.763E-20, 5.954E-20 /)
#endif
#if NSILICONIONS>4
      group_csn_silicon(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.050E-19 /)
      group_cse_silicon(:,5) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 8.415E-19 /)
#endif
#if NSILICONIONS>5
      group_csn_silicon(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.718E-19 /)
      group_cse_silicon(:,6) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.357E-19 /)
#endif
#if NSILICONIONS>6
      group_csn_silicon(:,7) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.270E-19 /)
      group_cse_silicon(:,7) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.919E-19 /)
#endif
#if NSILICONIONS>7
      group_csn_silicon(:,8) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.429E-20 /)
      group_cse_silicon(:,8) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.465E-19 /)
#endif
#if NSILICONIONS>8
      group_csn_silicon(:,9) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.417E-20 /)
      group_cse_silicon(:,9) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.226E-20 /)
#endif
#if NSILICONIONS>9
      group_csn_silicon(:,10) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 9.434E-21 /)
      group_cse_silicon(:,10) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.072E-20 /)
#endif
#if NSILICONIONS>10
      group_csn_silicon(:,11) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.251E-21 /)
      group_cse_silicon(:,11) = (/ 0.E0, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.447E-21 /)
#endif
   endif
  endif

  if (sulfur_ions) then
   group_csn_sulfur(:,:) = 0.d0
   group_cse_sulfur(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NSULFURIONS>0
      group_csn_sulfur(:,1) = (/ 0.000E+00, 0.000E+00, 6.347E-19, 3.619E-17, 4.280E-17, 3.563E-17, 4.047E-18, 6.033E-19 /)
      group_cse_sulfur(:,1) = (/ 0.000E+00, 0.000E+00, 9.662E-19, 3.626E-17, 4.283E-17, 3.421E-17, 3.443E-18, 3.795E-19 /)
#endif
#if NSULFURIONS>1
      group_csn_sulfur(:,2) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.011E-19, 2.087E-18, 5.553E-19 /)
      group_cse_sulfur(:,2) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.588E-19, 1.778E-18, 3.555E-19 /)
#endif
#if NSULFURIONS>2
      group_csn_sulfur(:,3) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.099E-19, 5.597E-19 /)
      group_cse_sulfur(:,3) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.733E-19, 3.568E-19 /)
#endif
#if NSULFURIONS>3
      group_csn_sulfur(:,4) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.161E-20, 4.506E-19 /)
      group_cse_sulfur(:,4) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.283E-20, 2.914E-19 /)
#endif
#if NSULFURIONS>4
      group_csn_sulfur(:,5) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.677E-19 /)
      group_cse_sulfur(:,5) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.341E-19 /)
#endif
#if NSULFURIONS>5
      group_csn_sulfur(:,6) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.541E-20 /)
      group_cse_sulfur(:,6) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.992E-20 /)
#endif
#if NSULFURIONS>6
      group_csn_sulfur(:,7) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.438E-19 /)
      group_cse_sulfur(:,7) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.652E-19 /)
#endif
#if NSULFURIONS>7
      group_csn_sulfur(:,8) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.624E-20 /)
      group_cse_sulfur(:,8) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.173E-19 /)
#endif
#if NSULFURIONS>8
      group_csn_sulfur(:,9) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.347E-20 /)
      group_cse_sulfur(:,9) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.045E-19 /)
#endif
#if NSULFURIONS>9
      group_csn_sulfur(:,10) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 9.793E-21 /)
      group_cse_sulfur(:,10) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.366E-20 /)
#endif
   endif
  endif

  if (iron_ions) then
   group_csn_iron(:,:) = 0.d0
   group_cse_iron(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NIRONIONS>0
      group_csn_iron(:,1) = (/ 0.000E+00, 0.000E+00, 6.545E-19, 2.372E-18, 2.824E-18, 3.503E-18, 5.055E-18, 1.215E-18 /)
      group_cse_iron(:,1) = (/ 0.000E+00, 0.000E+00, 8.104E-19, 2.377E-18, 2.827E-18, 3.552E-18, 5.199E-18, 5.311E-19 /)
#endif
#if NIRONIONS>1
      group_csn_iron(:,2) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.349E-18, 9.286E-18, 2.419E-18 /)
      group_cse_iron(:,2) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.533E-18, 9.656E-18, 1.120E-18 /)
#endif
#if NIRONIONS>2
      group_csn_iron(:,3) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.029E-18, 3.040E-18 /)
      group_cse_iron(:,3) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.982E-18, 1.492E-18 /)
#endif
#if NIRONIONS>3
      group_csn_iron(:,4) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.517E-18 /)
      group_cse_iron(:,4) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.901E-18 /)
#endif
#if NIRONIONS>4
      group_csn_iron(:,5) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.671E-18 /)
      group_cse_iron(:,5) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.159E-18 /)
#endif
#if NIRONIONS>5
      group_csn_iron(:,6) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 9.122E-19 /)
      group_cse_iron(:,6) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.922E-19 /)
#endif
#if NIRONIONS>6
      group_csn_iron(:,7) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.755E-19 /)
      group_cse_iron(:,7) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.049E-19 /)
#endif
#if NIRONIONS>7
      group_csn_iron(:,8) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.274E-19 /)
      group_cse_iron(:,8) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.616E-19 /)
#endif
#if NIRONIONS>8
      group_csn_iron(:,9) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.279E-19 /)
      group_cse_iron(:,9) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.907E-19 /)
#endif
#if NIRONIONS>9
      group_csn_iron(:,10) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 8.311E-20 /)
      group_cse_iron(:,10) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.040E-19 /)
#endif
#if NIRONIONS>10
      group_csn_iron(:,11) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.170E-20 /)
      group_cse_iron(:,11) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.371E-19 /)
#endif
#if NIRONIONS>11
      group_csn_iron(:,12) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.798E-20 /)
      group_cse_iron(:,12) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 8.085E-20 /)
#endif
#if NIRONIONS>12
      group_csn_iron(:,13) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.550E-20 /)
      group_cse_iron(:,13) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.726E-20 /)
#endif
#if NIRONIONS>13
      group_csn_iron(:,14) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.676E-21 /)
      group_cse_iron(:,14) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.464E-20 /)
#endif
#if NIRONIONS>14
      group_csn_iron(:,15) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.498E-21 /)
      group_cse_iron(:,15) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 5.211E-21 /)
#endif
#if NIRONIONS>15
      group_csn_iron(:,16) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.741E-22 /)
      group_cse_iron(:,16) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.276E-22 /)
#endif
   endif
  endif

  if (neon_ions) then
   group_csn_neon(:,:) = 0.d0
   group_cse_neon(:,:) = 0.d0
   if (nGroups.eq.8) then
#if NNEONIONS>0
      group_csn_neon(:,1) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.444E-18, 8.727E-18, 3.718E-18 /)
      group_cse_neon(:,1) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.754E-18, 8.778E-18, 2.352E-18 /)
#endif
#if NNEONIONS>1
      group_csn_neon(:,2) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.517E-18, 3.873E-18 /)
      group_cse_neon(:,2) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.114E-18, 2.430E-18 /)
#endif
#if NNEONIONS>2
      group_csn_neon(:,3) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.661E-18 /)
      group_cse_neon(:,3) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.910E-18 /)
#endif
#if NNEONIONS>3
      group_csn_neon(:,4) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 9.145E-19 /)
      group_cse_neon(:,4) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 9.186E-19 /)
#endif
#if NNEONIONS>4
      group_csn_neon(:,5) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.519E-19 /)
      group_cse_neon(:,5) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.489E-19 /)
#endif
#if NNEONIONS>5
      group_csn_neon(:,6) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.303E-19 /)
      group_cse_neon(:,6) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.091E-19 /)
#endif
#if NNEONIONS>6
      group_csn_neon(:,7) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.416E-20 /)
      group_cse_neon(:,7) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 8.846E-20 /)
#endif
#if NNEONIONS>7
      group_csn_neon(:,8) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.495E-20 /)
      group_cse_neon(:,8) = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 3.364E-20 /)
#endif
   endif
  endif
  ! END joki =============================================================

  if(isH2Katz) then
     if(rt_star) sedprops_update=1
     if(myid==1) write(*,*)
     if(myid==1) write(*,*)'>>>HKNote: no worries about initial group_csn,cse. will be updated on the fly'
     if(myid==1) write(*,*)
  endif

  do i=1,min(nIons,nGroups)
     spec2group(i)=i      ! Species contributions to groups (for non-OTSA)
  end do
  ! Read namelist file
  rewind(1)
  read(1,NML=rt_groups,END=101)
101 continue              ! no harm if no rt namelist

  if(minval(group_egy) .le. 0d0 .and. myid==1) then
     print*,'========================================================='
     print*,'WARNING! Some photon groups have zero or negative energy!'
     print*,'This could have unwanted effects, so be careful!!!'
     print*,'========================================================='
  endif

  call updateRTGroups_CoolConstants
  call write_group_props(.false.,6)
END SUBROUTINE read_rt_groups

!************************************************************************
SUBROUTINE add_rt_sources(ilevel,dt)

! Inject radiation from RT source regions (from the RT namelist). Since 
! the light sources are continuously emitting radiation, this is called
! continuously during code execution, rather than just during 
! initialization.
!
! ilevel => amr level at which to inject the radiation
! dt     => timestep for injection (since injected values are per time)
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use rt_hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dt
  integer::i,igrid,ncache,iskip,ngrid
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(rtdp),dimension(1:nvector,1:nrtvar),save::uu
!------------------------------------------------------------------------
  call add_UV_background(ilevel,dt)
  if(numbtot(1,ilevel)==0)return    ! no grids at this level
  if(rt_nsource .le. 0) return      ! no rt sources
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid
  ! dx (and dx_loc=dx) are just equal to 1/nx (where 1 is the boxlength)
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do
        ! Read the RT variables
        do ivar=1,nrtvar
           do i=1,ngrid
              uu(i,ivar)=rtunew(ind_cell(i),ivar)
           end do
        end do
        ! find injected values per cell
        call rt_sources_vsweep(xx,uu,dx_loc,dt,ngrid)
        ! Write the RT variables
        do ivar=1,nrtvar
           do i=1,ngrid
              rtunew(ind_cell(i),ivar)=uu(i,ivar)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_rt_sources for level ',I2)

END SUBROUTINE add_rt_sources

!************************************************************************
SUBROUTINE add_UV_background(ilevel,dt)

! Inject radiation from RT source regions (from the RT namelist). Since 
! the light sources are continuously emitting radiation, this is called
! continuously during code execution, rather than just during 
! initialization.
!
! ilevel => amr level at which to inject the radiation
! dt     => timestep for injection (since injected values are per time)
!------------------------------------------------------------------------
  use UV_module, ONLY: UV_Nphot_cgs, nUVgroups, iUVgroups
  use amr_commons
  use rt_parameters
  use hydro_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dt
  integer::i,igrid,ncache,iskip,ngrid,j
  integer::ind,ivar,ind_group,ic,ig
  integer ,dimension(1:nvector),save::ind_grid
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_np  &
            ,scale_fp,efactor,nH
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return     ! no grids at this level
  if(.not. rt_isDiffuseUVsrc) return ! no propagated UV background
  if(verbose)write(*,111)ilevel

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)

  ncache=active(ilevel)%ngrid
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ic=iskip+ind_grid(i) ! cell index
           ! Read the local gas density and inject the UV radiation
           nH=unew(ic,1)*scale_nH
           efactor = exp(-nH/rt_UVsrc_nHmax)
           do j=1,nUVgroups
              ig = iGroups(iUVgroups(j))
              rtunew(ic,ig) = max(rtunew(ic,ig)                     &
                                 ,UV_Nphot_cgs(j)/scale_np * efactor)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_UV_background for level ',I2)

END SUBROUTINE add_UV_background

!************************************************************************
SUBROUTINE rt_sources_vsweep(x,uu,dx,dt,nn)

! Do a vector sweep, injecting RT source regions into cells, that is if
! they are in any of these regions.
!
! x      =>  ncells*ndim: positions of grid cells
! uu    <=  ncells*nrtvars: injected rt variables in each cell
! dx     =>  real cell width in code units
! dt     =>  real timestep length in code units
! nn     =>  int number of cells
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  implicit none
  integer ::nn
  real(dp)::dx,dt,dx_cgs,dt_cgs
  real(rtdp),dimension(1:nvector,1:nrtvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::x
  integer::i,k,group_ind
  real(dp)::vol,r,xn,yn,zn,en
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_np,scale_fp
!------------------------------------------------------------------------
  ! Initialize everything to zero
  !  uu=0.0d0
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)
  dx_cgs=dx*scale_l
  dt_cgs=dt*scale_t
  ! Loop over RT regions
  do k=1,rt_nsource

     ! Find which photon group we should be contributing to
     if(rt_src_group(k) .le. 0 .or. rt_src_group(k) .gt. nGroups) cycle
     group_ind = iGroups(rt_src_group(k))
     ! For "square" regions only:
     if(rt_source_type(k) .eq. 'square')then
       ! Exponent of choosen norm
        en=rt_exp_source(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-rt_src_x_center(k))/rt_src_length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-rt_src_y_center(k))/rt_src_length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-rt_src_z_center(k))/rt_src_length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(rt_exp_source(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region, inject value
           if(r<1.0)then
              uu(i,group_ind) = rt_n_source(k)/rt_c_cgs/scale_Np
              ! The input flux is the fraction Fp/(c*Np) (Max 1 magnitude)
              uu(i,group_ind+1) =                                       &
                        rt_u_source(k) * rt_c * rt_n_source(k) / scale_Np
#if NDIM>1 
              uu(i,group_ind+2) =                                       &
                        rt_v_source(k) * rt_c * rt_n_source(k) / scale_Np
#endif
#if NDIM>2
              uu(i,group_ind+3) =                                       &
                        rt_w_source(k) * rt_c * rt_n_source(k) / scale_Np
#endif
           end if
        end do
     end if

     ! For "point" regions only:
     if(rt_source_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx_cgs**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-rt_src_x_center(k))/dx, 0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-rt_src_y_center(k))/dx, 0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-rt_src_z_center(k))/dx, 0.0_dp)
#endif
           r=xn*yn*zn
           if(r .gt. 0.) then
              ! If cell lies within CIC cloud, inject value.
              ! Photon input is in # per sec...need to convert to uu
              uu(i,group_ind)=uu(i,group_ind)                            &
                            + rt_n_source(k) / scale_Np * r / vol * dt_cgs
              uu(i,group_ind+1)=uu(i,group_ind+1) + rt_u_source(k) *rt_c &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#if NDIM>1
              uu(i,group_ind+2)=uu(i,group_ind+2) + rt_v_source(k) *rt_c &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#endif
#if NDIM>2
              uu(i,group_ind+3)=uu(i,group_ind+3) + rt_w_source(k) *rt_c &
                            * rt_n_source(k) / scale_Np * r / vol * dt_cgs
#endif
           endif
        end do
     end if

     ! For shell regions only:
     if(rt_source_type(k) .eq. 'shell')then
        ! An emitting spherical shell with center coordinates given,
        ! along with inner and outer radius (rt_src_length_x,z).
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=0.0; yn=0.0; zn=0.0
           xn=max(abs(x(i,1)-rt_src_x_center(k)), 0.0_dp)
#if NDIM>1
           yn=max(abs(x(i,2)-rt_src_y_center(k)), 0.0_dp)
#endif
#if NDIM>2
           zn=max(abs(x(i,3)-rt_src_z_center(k)), 0.0_dp)
#endif
           r=sqrt(xn**2+yn**2+zn**2)
           if(r .gt. rt_src_length_x(k) .and. &
                r .lt. rt_src_length_y(k)) then
              ! If cell lies within CIC cloud, inject value
              ! photon input is in # per sec...need to convert to uu
              uu(i,group_ind)=rt_n_source(k) / scale_np
           endif
        end do
     end if
  end do

  return
END SUBROUTINE rt_sources_vsweep



