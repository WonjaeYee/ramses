! Non-equlibrium (in H and He) cooling module for radiation-hydrodynamics.
! For details, see Rosdahl et al. 2013, and Rosdahl & Teyssier 2015.
! Joki Rosdahl, Andreas Bleuler, and Romain Teyssier, September 2015.

module rt_cooling_module
   use amr_commons, only: myid
   use cooling_module, only: X, Y, cmp_metals_cloudy, cmp_metals_os13
   use rt_parameters
   use coolrates_module
   use htmcool, only: cmp_metals_ht_harley
   implicit none

   private   ! default

   public rt_set_model, rt_solve_cooling, update_UVrates, cmp_chem_eq &
      , isHe, is_mu_H2, X, Y, rhoc, kB, mH, T2_min_fix, twopi &
      , signc, sigec, PHrate, UVrates, rt_isIR, kappaAbs, kappaSc &
      , is_kIR_T, iIR, rt_isIRtrap, iIRtrapVar, rt_pressBoost &
      , rt_isoPress, rt_T_rad, rt_vc, a_r, getMu &
      , signc_oxygen, sigec_oxygen &
      , signc_nitrogen, sigec_nitrogen &
      , signc_carbon, sigec_carbon &
      , signc_magnesium, sigec_magnesium &
      , signc_silicon, sigec_silicon &
      , signc_sulfur, sigec_sulfur &
      , signc_iron, sigec_iron &
      , signc_neon, sigec_neon &
      , signc_dust &
      , amu_to_g, mO_NIST_amu, mN_NIST_amu, mC_NIST_amu, mMg_NIST_amu &
      , mSi_NIST_amu, mS_NIST_amu, mFe_NIST_amu, mNe_NIST_amu, mCa_NIST_amu &
      , PHrate_oxygen, PHrate_nitrogen, PHrate_carbon, PHrate_magnesium &
      , PHrate_silicon, PHrate_sulfur, PHrate_iron, PHrate_neon 

   ! NOTE: T2=T/mu
   ! Np = photon density, Fp = photon flux,

   real(dp), parameter::rhoc = 1.88000d-29        ! Crit. density [g cm-3]
   real(dp), parameter::mH = 1.66000d-24          ! H atom mass [g]
   real(dp), parameter::kB = 1.38062d-16          ! Boltzm.const. [erg K-1]
   real(dp), parameter::a_r = 7.5657d-15          ! Rad.const. [erg cm-3 K-4]
   real(dp), parameter::mu_mol = 1.2195D0
   real(dp), parameter::T2_min_fix = 1.d-2        ! Min temperature [K]
   real(dp), parameter::twopi = 6.2831853d0       ! Two times pi
   real(dp), parameter::amu_to_g = 1.66054d-24    ! atomic mass units in grams
   real(dp), parameter::mO_NIST_amu = 15.9994d0     ! oxygen molecular weight [amu]
   real(dp), parameter::mN_NIST_amu = 14.0067d0   ! nitrogen molecular weight [amu]
   real(dp), parameter::mC_NIST_amu = 12.0107d0     ! carbon molecular weight [amu]
   real(dp), parameter::mMg_NIST_amu = 24.305d0   ! magnesium molecular weight [amu]
   real(dp), parameter::mSi_NIST_amu = 28.0855d0    ! silicon molecular weight [amu]
   real(dp), parameter::mS_NIST_amu = 32.065d0      ! sulfur molecular weight [amu]
   real(dp), parameter::mFe_NIST_amu = 55.854d0        ! iron molecular weight [amu]
   real(dp), parameter::mNe_NIST_amu = 20.1797d0       ! neon molecular weight [amu]
   real(dp), parameter::mCa_NIST_amu = 40.078d0        ! calcium molecular weight [amu]

   real(dp)::T_min, T_frac, x_MIN, x_frac, Np_min, Np_frac, Fp_min, Fp_frac, x_MIN_HI

   integer, parameter::iIR = 1                       !          IR group index
   integer::iIRtrapVar = 1                          ! Trapped IR energy index
   ! Namelist parameters:
   logical::is_mu_H2 = .false.
   logical::rt_isoPress = .false.         ! Use cE, not F, for rad. pressure
   real(dp)::rt_pressBoost = 1d0          ! Boost on RT pressure
   logical::rt_isIR = .false.             ! Using IR scattering on dust?
   logical::rt_isIRtrap = .false.         ! IR trapping in NENER variable?
   logical::is_kIR_T = .false.            ! k_IR propto T^2?
   logical::rt_T_rad = .false.            ! Use T_gas = T_rad
   logical::rt_vc = .false.               ! (semi-) relativistic RT
   real(dp)::Tmu_dissoc = 1d3             ! Dissociation temperature [K]
   real(dp), dimension(nGroups)::kappaAbs = 0! Dust absorption opacity
   real(dp), dimension(nGroups)::kappaSc = 0 ! Dust scattering opacity

   ! Cooling constants, updated on SED and c-change [cm3 s-1],[erg cm3 s-1]
   real(dp), dimension(nGroups, nIons)::signc, sigec, PHrate
   real(dp), dimension(nGroups, n_oxygen_ions)::signc_oxygen, sigec_oxygen, PHrate_oxygen
   real(dp), dimension(nGroups, n_nitrogen_ions)::signc_nitrogen, sigec_nitrogen, PHrate_nitrogen
   real(dp), dimension(nGroups, n_carbon_ions)::signc_carbon, sigec_carbon, PHrate_carbon
   real(dp), dimension(nGroups, n_magnesium_ions)::signc_magnesium, sigec_magnesium, PHrate_magnesium
   real(dp), dimension(nGroups, n_silicon_ions)::signc_silicon, sigec_silicon, PHrate_silicon
   real(dp), dimension(nGroups, n_sulfur_ions)::signc_sulfur, sigec_sulfur, PHrate_sulfur
   real(dp), dimension(nGroups, n_iron_ions)::signc_iron, sigec_iron, PHrate_iron
   real(dp), dimension(nGroups, n_neon_ions)::signc_neon, sigec_neon, PHrate_neon
   real(dp), dimension(nGroups, 2)::signc_dust

   ! For isH2Katz=.true.
   real(dp), dimension(nIons + 1, 2)::UVrates    !UV backgr. heating/ion. rates

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE rt_set_model(Nmodel, J0in_in, J0min_in, alpha_in &
                           , normfacJ0_in, zreioniz_in, correct_cooling, realistic_ne, h &
                           , omegab, omega0, omegaL, astart_sim, T2_sim)
! Initialize cooling. All these parameters are unused at the moment and
! are only there for the original cooling-module.
! Nmodel(integer)     =>     Model for UV background and metals
! J0in_in  (dble)     => Default UV intensity
! J0min_in (dble)     => Minimum UV intensity
! alpha_in (dble)     => Slope of the UV spectrum
! zreioniz_in (dble)  => Reionization redshift
! normfacJ0_in (dble) => Normalization factor fot a Harrdt&Madau UV model
! correct_cooling (integer) => Cooling correction
! realistic_ne (integer) => Use realistic electron density at high z?
! h (dble)            => H0/100
! omegab (dble)       => Omega baryons
! omega0 (dble)       => Omega materal total
! omegaL (dble)       => Omega Lambda
! astart_sim (dble)   => Redshift at which we start the simulation
! T2_sim (dble)      <=  Starting temperature in simulation?
!-------------------------------------------------------------------------
      use UV_module
      use coolrates_module, only: init_coolrates_tables, init_metal_atomic_data_tables
      use charge_transfer, only: init_ct_tables
      use fscool, only: init_fine_structure_tables
      real(kind=8) :: J0in_in, zreioniz_in, J0min_in, alpha_in, normfacJ0_in
      real(kind=8) :: astart_sim, T2_sim, h, omegab, omega0, omegaL
      integer  :: Nmodel, correct_cooling, realistic_ne, ig
      real(kind=8) :: astart = 0.0001, aend, dasura, T2end = T2_min_fix, mu, ne
!-------------------------------------------------------------------------
      if (myid == 1) write (*, *) &
         '==================RT momentum pressure is turned ON=============='
      if (myid == 1 .and. rt_isIR) &
         write (*, *) 'There is an IR group, with index ', iIR
      if (myid == 1 .and. rt_isIRtrap) write (*, *) &
         '=========IR trapping is turned ON=============='
      ! do initialization
      isHe = .true.; if (Y .le. 0.) isHe = .false.
      T_MIN = 1.d-1                  !                      Minimum T2
      T_FRAC = 0.1

      x_MIN = 1.d-6                !    Minimum ionization fractions
      x_MIN_HI = 1.d-10            !    Minimum neutral fractions
                                   !    -> be careful with large numbers (1e-6)
                                   !    bc residual HI can act as a constant
                                   !    heating term at high density environ
      x_FRAC = 0.1

      Np_MIN = 1.d-13                        !            Photon density floor
      Np_FRAC = 0.2

      Fp_MIN = 1D-13*rt_c_cgs               !           Minimum photon fluxes
      Fp_FRAC = 0.5

      ! Calculate initial temperature
      if (astart_sim < astart) then
         write (*, *) 'ERROR in set_model : astart_sim is too small.'
         write (*, *) 'astart     =', astart
         write (*, *) 'astart_sim =', astart_sim
         STOP
      end if
      aend = astart_sim
      dasura = 0.02d0

      call update_rt_c
      call init_UV_background
      if (uvbg_rtz) call init_UV_background_RTZ ! Haardt madau 2012

      call init_metal_atomic_data_tables()
      if (include_charge_transfer) call init_ct_tables() !initialize ct tables
      if (low_t_noneq_cooling) call init_fine_structure_tables() ! initialize fine structure cooling tables

      if (cosmo) then
         call update_UVrates(aexp)! In case of cosmo
         call init_coolrates_tables(aexp)
      else
         call update_UVrates(astart_sim)! In case of aexp_nocosmo
         call init_coolrates_tables(astart_sim)
      end if
      if (uvbg_rtz) call update_UV_background_RTZ(1./aexp - 1.) ! Haardt madau 2012

      if (nrestart == 0 .and. cosmo) &
         call rt_evol_single_cell(astart, aend, dasura, h, omegab, omega0 &
                                  , omegaL, -1.0d0, T2end, mu, ne, .false.)
      T2_sim = T2end

   END SUBROUTINE rt_set_model

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE update_UVrates(aexp)
! Set the UV ionization and heating rates according to the given a_exp.
!-------------------------------------------------------------------------
      use UV_module
      use amr_parameters, only: haardt_madau
      integer::i
      real(dp)::aexp
!------------------------------------------------------------------------
      UVrates = 0.
      if (.not. haardt_madau) RETURN

      call inp_UV_rates_table(1./aexp - 1., UVrates, .true.)

!  if(myid==1) then
!     write(*,*) 'The UV rates have changed to:'
!     do i=1,nIons
!        write(*,910) UVrates(i,:)
!     enddo
!  endif
910   format(1pe21.6, ' s-1', 1pe21.6, ' erg s-1')

   END SUBROUTINE update_UVrates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE rt_solve_cooling(T2, xion, Np, Fp, p_gas, dNpdt, dFpdt         &
                            &  , nH, rho, c_switch, Zsolar, dt, a_exp, nCell &
                            &  , dx_SS, xO, xN, xC, xMg, xSi, xS, xFe, xNe   &
                            &  , nOxygen, nNitrogen, nCarbon, nMagnesium     &
                            &  , nSilicon, nSulfur, nIron, nNeon, nco_mol    &
                            &  , dx_loc, cooling_time                        &
                            &  , temperature_save, cooling_save)
! Semi-implicitly solve for new temperature, ionization states,
! photon density/flux, and gas velocity in a number of cells.
! Parameters:
! T2     <=> T/mu [K]
! xion   <=> NION ionization fractions
! Np     <=> NGROUPS photon number densities [cm-3]
! Fp     <=> NGROUPS * ndim photon number fluxes [cm-2 s-1]
! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
! dNpdt   =>  Op split increment in photon densities during dt
! dFpdt   =>  Op split increment in photon flux magnitudes during dt
! nH      =>  Hydrogen number densities [cm-3]
! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
! Zsolar  =>  Cell metallicities [solar fraction]
! dt      =>  Timestep size             [s]
! a_exp   =>  Cosmic expansion
! nCell   =>  Number of cells (length of all the above vectors)
! dx_SS   =>  length scale for self-shielding: NH = nH*dx_SS
! dx_loc  =>  Cell length [cm]
! cooling_time  =>  Return the cooling time [s]
! We use a slightly modified method of Anninos et al. (1997).
!-------------------------------------------------------------------------
      use amr_commons
      implicit none
      real(dp), dimension(1:nvector):: T2
      real(dp), dimension(1:nIons, 1:nvector):: xion
      real(dp), dimension(1:nGroups, 1:nvector):: Np, dNpdt
      real(dp), dimension(1:ndim, 1:nGroups, 1:nvector):: Fp, dFpdt
      real(dp), dimension(1:ndim, 1:nvector):: p_gas
      real(dp), dimension(1:nvector):: nH, Zsolar, dx_SS, rho
      logical, dimension(1:nvector):: c_switch
      real(dp)::dt, a_exp
      integer::ncell !--------------------------------------------------------
      real(dp), dimension(1:nvector):: tLeft, ddt
      logical:: dt_ok
      real(dp)::dt_rec
      real(dp)::dx_loc
      real(dp):: dT2
      real(dp), dimension(nIons):: dXion
      real(dp), dimension(nGroups):: dNp
      real(dp), dimension(1:ndim, 1:nGroups):: dFp
      real(dp), dimension(1:ndim):: dp_gas
      integer::i, ia, ig, nAct, nAct_next, loopcnt, code
      integer, dimension(1:nvector):: indAct              ! Active cell indexes
      real(dp)::one_over_rt_c_cgs, one_over_egy_IR_erg, one_over_x_FRAC
      real(dp)::one_over_Np_FRAC, one_over_Fp_FRAC, one_over_T_FRAC
      real(dp), dimension(1:nGroups) :: group_egy_ratio, group_egy_erg

      ! Metals
      real(dp), dimension(1:n_oxygen_ions, 1:nvector):: xO
      real(dp), dimension(1:n_nitrogen_ions, 1:nvector):: xN
      real(dp), dimension(1:n_carbon_ions, 1:nvector):: xC
      real(dp), dimension(1:n_magnesium_ions, 1:nvector):: xMg
      real(dp), dimension(1:n_silicon_ions, 1:nvector):: xSi
      real(dp), dimension(1:n_sulfur_ions, 1:nvector):: xS
      real(dp), dimension(1:n_iron_ions, 1:nvector):: xFe
      real(dp), dimension(1:n_neon_ions, 1:nvector):: xNe
      real(dp), dimension(n_oxygen_ions):: dxO
      real(dp), dimension(n_nitrogen_ions):: dxN
      real(dp), dimension(n_carbon_ions):: dxC
      real(dp), dimension(n_magnesium_ions):: dxMg
      real(dp), dimension(n_silicon_ions):: dxSi
      real(dp), dimension(n_sulfur_ions):: dxS
      real(dp), dimension(n_iron_ions):: dxFe
      real(dp), dimension(n_neon_ions):: dxNe
      real(dp), dimension(1:nvector)::nOxygen, nNitrogen, nCarbon, nMagnesium ! number densities of each species
      real(dp), dimension(1:nvector)::nSilicon, nSulfur, nIron, nNeon, nco_mol ! number densities of each species
      real(dp):: min_z_neq=1.d-8
      real(dp):: CO_SS
      
      ! Store the cooling rate
      real(dp), dimension(1:nvector)::cooling_time, temperature_save
      real(dp), dimension(ncoolheatratesave,1:nvector)::cooling_save


      ! Store some temporary variables reduce computations
      one_over_rt_c_cgs = 1d0/rt_c_cgs
      one_over_Np_FRAC = 1d0/Np_FRAC
      one_over_Fp_FRAC = 1d0/Fp_FRAC
      one_over_T_FRAC = 1d0/T_FRAC
      one_over_x_FRAC = 1d0/x_FRAC
#if NGROUPS>0
      if (rt .and. nGroups .gt. 0) then
         group_egy_erg(1:nGroups) = group_egy(1:nGroups)*ev_to_erg
         if (rt_isIR) then
            group_egy_ratio(1:nGroups) = group_egy(1:nGroups)/group_egy(iIR)
            one_over_egy_IR_erg = 1.d0/group_egy_erg(iIR)
         end if
      end if
#endif
      !-----------------------------------------------------------------------
      tleft(1:ncell) = dt                !       Time left in dt for each cell
      ddt(1:ncell) = dt                  ! First guess at sub-timestep lengths
      do i = 1, ncell
         indact(i) = i                   !      Set up indexes of active cells
         ! Ensure all state vars are legal:
         T2(i) = MAX(T2(i), T2_min_fix)

         call reduce_xion(xion(1:nIons, i), x_MIN, x_MIN_HI) ! Taysun

         ! ensure metal ionization fractions sum to 1 (what about negatives?)
         if (oxygen_ions)    xO(:,i)  = xO(:,i)  / SUM(xO(:,i))
         if (nitrogen_ions)  xN(:,i)  = xN(:,i)  / SUM(xN(:,i))
         if (carbon_ions)    xC(:,i)  = xC(:,i)  / SUM(xC(:,i))
         if (magnesium_ions) xMg(:,i) = xMg(:,i) / SUM(xMg(:,i))
         if (silicon_ions)   xSi(:,i) = xSi(:,i) / SUM(xSi(:,i))
         if (sulfur_ions)    xS(:,i)  = xS(:,i)  / SUM(xS(:,i))
         if (iron_ions)      xFe(:,i) = xFe(:,i) / SUM(xFe(:,i))
         if (neon_ions)      xNe(:,i) = xNe(:,i) / SUM(xNe(:,i))

         if (rt) then
            do ig = 1, ngroups
               Np(ig, i) = MAX(smallNp, Np(ig, i))
               call reduce_flux(Fp(:, ig, i), Np(ig, i)*rt_c_cgs)
            end do
         end if
      end do

      ! Loop until all cells have tleft=0
      ! **********************************************
      nAct = nCell                                      ! Currently active cells
      loopcnt = 0; n_cool_cells = n_cool_cells + nCell     !             Statistics

      do while (nAct .gt. 0)      ! Iterate while there are still active cells
         loopcnt = loopcnt + 1; tot_cool_loopcnt = tot_cool_loopcnt + nAct
         nAct_next = 0                     ! Active cells for the next iteration
         do ia = 1, nAct                             ! Loop over the active cells
            i = indAct(ia)                        !                 Cell index
            call cool_step(i, loopcnt)
            
            if (loopcnt .gt. 100000) then
               write(*,*) loopCodes
               write(*,*) "--------------------"
               write(*,*) "Cooling rates:",cooling_save(:,i)
               write(*,*) "--------------------"
               if (oxygen_ions)    write(*,*) "Oxygen:",xO(:,i)
               if (nitrogen_ions)  write(*,*) "Nitrogen:",xN(:,i)
               if (carbon_ions)    write(*,*) "Carbon:",xC(:,i)
               if (magnesium_ions) write(*,*) "Magnesium:",xMg(:,i)
               if (silicon_ions)   write(*,*) "Silicon:",xSi(:,i)
               if (sulfur_ions)    write(*,*) "Sulfur:",xS(:,i)
               if (iron_ions)      write(*,*) "Iron:",xFe(:,i)
               if (neon_ions)      write(*,*) "Neon:",xNe(:,i)
               write(*,*) "--------------------"
               call display_coolinfo(.true., loopcnt, i, dt - tleft(i), dt &
                                     , ddt(i), nH(i), T2(i), xion(:, i), Np(:,i) &
                                     , Fp(:, :, i), p_gas(:, i) &
                                     , dT2, dXion, dNp, dFp, dp_gas, code, cooling_save(:,i) &
                                     , nOxygen(i), nNitrogen(i), nCarbon(i), nMagnesium(i) &
                                     , nSilicon(i), nSulfur(i), nIron(i), nNeon(i), nco_mol(i))
            end if
            if (.not. dt_ok) then
               ddt(i) = ddt(i)/2.                    ! Try again with smaller dt
               nAct_next = nAct_next + 1; indAct(nAct_next) = i
               loopCodes(code) = loopCodes(code) + 1
               cycle
            end if
            ! Update the cell state (advance the time by ddt):
            T2(i) = T2(i) + dT2
            xion(:, i) = xion(:, i) + dXion(:)

            if (Zsolar(i).gt.min_z_neq) then
               if (oxygen_ions)    xO(:,i)  = xO(:,i)  + dxO(:)
               if (nitrogen_ions)  xN(:,i)  = xN(:,i)  + dxN(:)
               if (carbon_ions)    xC(:,i)  = xC(:,i)  + dxC(:)
               if (magnesium_ions) xMg(:,i) = xMg(:,i) + dxMg(:) 
               if (silicon_ions)   xSi(:,i) = xSi(:,i) + dxSi(:) 
               if (sulfur_ions)    xS(:,i)  = xS(:,i)  + dxS(:) 
               if (iron_ions)      xFe(:,i) = xFe(:,i) + dxFe(:) 
               if (neon_ions)      xNe(:,i) = xNe(:,i) + dxNe(:) 
            end if

            if (nGroups .gt. 0) then
               Np(:, i) = Np(:, i) + dNp(:)
               Fp(:, :, i) = Fp(:, :, i) + dFp(:, :)
            end if
            p_gas(:, i) = p_gas(:, i) + dp_gas(:)

            tleft(i) = tleft(i) - ddt(i)
            if (tleft(i) .gt. 0.) then           ! Not finished with this cell
               nAct_next = nAct_next + 1; indAct(nAct_next) = i
            else if (tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
               print *, 'In rt_solve_cooling: tleft < 0  !!'
               stop
            end if
            ddt(i) = min(dt_rec, tleft(i))    ! Use recommended dt from cool_step
         end do ! end loop over active cells
         nAct = nAct_next
      end do ! end iterative loop
      ! loop statistics
      max_cool_loopcnt = max(max_cool_loopcnt, loopcnt)
   contains

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE cool_step(icell, loopcnt)
         ! Compute change in cell state in timestep ddt(icell), or set in dt_rec
         ! a recommendation for new timestep if ddt(icell) proves too large.
         ! T2      => T/mu [K]                                -- dT2 is new value
         ! xion    => NION ionization fractions               --     dXion is new
         ! Np      => NGROUPS photon number densities [cm-3]  -- dNp is new value
         ! Fp      => NGROUPS * ndim photon fluxes [cm-2 s-1] -- dFp is new value
         ! p_gas   => ndim gas momenta [cm s-1 g cm-3]        --    dp_gas is new
         ! dNpdt   =>  Op split increment in photon densities during dt
         ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
         ! nH      =>  Hydrogen number densities [cm-3]
         ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
         ! Zsolar  =>  Cell metallicities [solar fraction]
         ! dt      =>  Timestep size [s]
         ! a_exp   =>  Cosmic expansion
         ! dt_ok   <=  .f. if timestep constraints were broken, .t. otherwise
         ! dt_rec  <=  Recommended timesteps for next iteration
         ! code    <= Error code in cool step, if dt_ok=.f.
         !
         ! The original values, T2, xion etc, must stay unchanged, while dT2,
         ! dxion etc contain the new values (the difference at the end of the
         ! routine).
         !-----------------------------------------------------------------------
         use amr_commons
         use const
         use constants, only: ev2erg
         use coolrates_module
         use chemistry_module
         use charge_transfer
         use fscool
         implicit none
         integer, intent(in)::icell
         real(dp), dimension(nDim), save:: dmom
         real(dp), dimension(nDim), save:: u_gas ! Gas velocity
         real(dp), dimension(nIons), save:: alpha, beta, nN, nI
         real(dp), save:: dUU, fracMax
         real(dp), save:: xHeI, mu, TK, nHe, ne, neInit, Hrate, dAlpha, dBeta, ne_metals
         real(dp), save:: s, jac, q, Crate, dCdT2, X_nHkb, rate, dRate, cr, de
         real(dp), save:: photoRate, metal_tot, metal_prime, ss_factor, fire_tot, fire_prime
         real(dp), save:: fine_structure, fine_structure_prime
         real(dp), save:: metal_cool_smooth_f1,metal_cool_smooth_f2
         real(dp), save:: CO_cooling, CO_cooling_prime
         integer, save:: iion, igroup, idim
         real(dp), dimension(nGroups), save:: recRad, phAbs, phSc, dustAbs
         real(dp), dimension(nGroups), save:: dustSc, kAbs_loc, kSc_loc
         real(dp), save:: TR, one_over_C_v, E_rad, dE_T, fluxMag, mom_fact
         real(dp), save:: xHI, xHII, xH2, xe, xHeII, xHeIII
         real(dp), save:: kUV, HrateLW, HrateLW_prime, Dlw, Flw, sigUVthin, sigUV, ss_H2
         real(dp), save:: Ebpump, Ebpump_prime, EUV, fpump, f_shd
         real(dp), save:: f_dg, F_pe, G_0, G_conv, eps_PE, HratePE, HratePE_tot
         real(dp), save:: HratePE_prime,HratePE_prime_a,HratePE_prime_b
         real(dp), save:: dNpPE, dNpLW, dNpAdd, phAbsLW
         real(dp), save:: dust2metal_solar
         integer::loopcnt
         real(dp)::gH, scale_np, scale_fp
         integer::im
         real(dp)::min_cool_ion,nH2_cool

         ! for dust depletion
         real(dp)::d_mass,g_mass,dust_to_gas_mass_ratio_over_mw
         real(dp),dimension(1:nmetals)::depletion_factors
         ! for dust cooling
         real(dp)::cool_dust,cool_dust_prime
         real(dp)::T_dust,dust_hc_const
         real(dp)::beta_drc,dust_recomb_cool,dust_recomb_cool_prime
         real(dp)::dust_recomb_cool_loc,dust_recomb_cool_prime_loc
         real(dp)::H2heating,H2heating_prime
         ! charge transfer heating and cooling
         real(dp)::ct_heat_cool,ct_heat_cool_prime

         ! for cosmic rays
         real(dp)::effective_n,q_cr,phi_s
         real(dp)::HI_cosmic_ray_ionization_rate,He_cosmic_ray_ionization_rate
         real(dp)::H2_cosmic_ray_ionization_rate
         real(dp)::HI_cosmic_ray_ionization_rate_primary,He_cosmic_ray_ionization_rate_primary
         real(dp)::H2_cosmic_ray_ionization_rate_primary
         real(dp)::cr_heating_tot
         real(dp)::fsc_OI,fsc_OII,fsc_OIII,fsc_FeI,fsc_FeII,fsc_SiI,fsc_SiII
         real(dp)::fsc_NII,fsc_SI,fsc_NeII,fsc_CI,fsc_CII

         min_cool_ion = 1.d-10

         call rt_units(scale_np, scale_fp)
         !---------------------------------------------------------------------
         dt_ok = .false.
         nHe = 0.25*nH(icell)*Y/X  !         Helium number density
         ! U contains the original values, dU the updated ones:
         dT2 = T2(icell); dXion(:) = xion(:, icell); dNp(:) = Np(:, icell)
         dXion(:) = xion(:, icell); dNp(:) = Np(:, icell)
         dFp(:, :) = Fp(:, :, icell); dp_gas(:) = p_gas(:, icell)

         if (oxygen_ions)    dxO(:)  = xO(:, icell)
         if (nitrogen_ions)  dxN(:)  = xN(:, icell)
         if (carbon_ions)    dxC(:)  = xC(:, icell)
         if (magnesium_ions) dxMg(:) = xMg(:, icell)
         if (silicon_ions)   dxSi(:) = xSi(:, icell)
         if (sulfur_ions)    dxS(:)  = xS(:, icell)
         if (iron_ions)      dxFe(:) = xFe(:, icell)
         if (neon_ions)      dxNe(:) = xNe(:, icell)

         ! nN(1) == nN(ixHII)   == nHI           ! nI(1) == nI(ixHII)   == nHII
         ! nN(2) == nN(ixHeII)  == nHeI          ! nI(2) == nI(ixHeII)  == nHeII
         ! nN(3) == nN(ixHeIII) == nHeII         ! nI(3) == nI(ixHeIII) == nHeIII
         call reduce_xion(dXion, x_MIN, x_MIN_HI) ! Taysun
         xHII = dXion(ixHII)
         xHI = 1.d0 - xHII
         nN(ixHII) = nH(icell)*xHI                                  !     nHI
         nI(ixHII) = nH(icell)*dXion(ixHII)                         !    nHII
         if (isHe) then
            xHeII = dXion(ixHeII)
            xHeIII = dXion(ixHeIII)
            xHeI = MAX(1.d0 - xHeII - xHeIII, 0d0)
            nN(ixHeII) = nHe*xHeI                                   !    nHeI
            nN(ixHeIII) = nHe*xHeII                                 !   nHeII
            nI(ixHeII) = nHe*xHeII                                  !   nHeII
            nI(ixHeIII) = nHe*xHeIII                                !  nHeIII
         end if

         if (isH2) then
            ! nN(1) == nN(ixHI)    == nH2           ! nI(1) == nI(ixHI)    == nHI
            ! nN(2) == nN(ixHII)   == nHI           ! nI(2) == nI(ixHII)   == nHII
            ! nN(3) == nN(ixHeII)  == nHeI          ! nI(3) == nI(ixHeII)  == nHeII
            ! nN(4) == nN(ixHeIII) == nHeII         ! nI(4) == nI(ixHeIII) == nHeIII
            xHI = dXion(ixHI)
            nN(ixHII) = nH(icell)*xHI
            xH2 = MIN(MAX((1.d0 - xHI - xHII)/2., 0d0), 0.5d0)
            nN(ixHI) = nH(icell)*xH2
         end if
         mu = getMu(xHII, xHeII, xHeIII, dT2)
         TK = dT2*mu                                           !  Temperature
         if (rt_isTconst) TK = rt_Tconst                       !  Force constant T
         ne = nH(icell)*xHII + nHE*(xHeII + 2.*xHeIII)         !  Electron density

         ! Update the electron density for the presence of metals
         ne_metals = 0.d0
         if (Zsolar(icell).gt.min_z_neq) then
         if (oxygen_ions) then
            do im=2,n_oxygen_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nOxygen(icell) * xO(im,icell) * (im-1))
            enddo
         endif
         if (nitrogen_ions) then
            do im=2,n_nitrogen_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nNitrogen(icell) * xN(im,icell) * (im-1))
            enddo
         endif
         if (carbon_ions) then
            do im=2,n_carbon_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nCarbon(icell) * xC(im,icell) * (im-1))
            enddo
         endif
         if (magnesium_ions) then
            do im=2,n_magnesium_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nMagnesium(icell) * xMg(im,icell) * (im-1))
            enddo
         endif
         if (silicon_ions) then
            do im=2,n_silicon_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nSilicon(icell) * xSi(im,icell) * (im-1))
            enddo
         endif
         if (sulfur_ions) then
            do im=2,n_sulfur_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nSulfur(icell) * xS(im,icell) * (im-1))
            enddo
         endif
         if (iron_ions) then
            do im=2,n_iron_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nIron(icell) * xFe(im,icell) * (im-1))
            enddo
         endif
         if (neon_ions) then
            do im=2,n_neon_ions ! start index at 2 (ground state = no free electrons)
               ne_metals = ne_metals + (nNeon(icell) * xNe(im,icell) * (im-1))
            enddo
         endif
         endif
         ne = ne + ne_metals

         neInit = ne
         fracMax = 0d0   ! Max fractional update, to check if dt can be increased
         ss_factor = 1d0                    ! UV background self_shielding factor
         if (self_shielding) ss_factor = exp(-nH(icell)/1d-2)

         if (dust2metal_RR14) then
            call dust_to_metal_RR14(Zsolar(icell), dust2metal_solar)
         else
            dust2metal_solar = 1.0
         end if

         ! Measure dust depletion
         ! Depletion: Hardcoded. 1=Fe, 2=O, 3=N, 4=Mg, 5=Ne, 6=Si, 7=Ca, 8=C, 9=S
         depletion_factors = 1.d0
         dust_to_gas_mass_ratio_over_mw = 0.d0
         if (TK .le. T_sputter .and. simple_dust_depletion .and. dx_loc.gt.0.d0) then ! no dust at T>T_sputter
            d_mass = 0.d0
            g_mass = rho(icell) * (dx_loc**3.d0) ! gas mass in grams
            call get_dust_mass_and_depletion(depletion_factors, dx_loc, d_mass, g_mass & 
                                             , nH(icell), nIron(icell), nOxygen(icell) & 
                                             , nNitrogen(icell), nMagnesium(icell) & 
                                             , nNeon(icell), nSilicon(icell) &
                                             , nCarbon(icell), nSulfur(icell))
            dust_to_gas_mass_ratio_over_mw = (d_mass/g_mass) / (1.d0/162.d0)
         else if (.not.simple_dust_depletion) then
            dust_to_gas_mass_ratio_over_mw = Zsolar(icell)
         end if
         dust_to_gas_mass_ratio_over_mw = MIN(dust_to_gas_mass_ratio_over_mw,3.d0) ! Cap metallicity at 3Zsun

         ! necessary only for isH2Katz
         Dlw = 5.18d-11        ! Rollig 2007
         Flw = 2.1d7           ! Draine & Bertoldi 1996
         sigUVthin = Dlw/Flw   ! Effective cross section in the optically thin limit [cm2]
         fpump = 6.94d0        ! Pumping fraction in the ISM: Draine and Bertoldi 1996
         Ebpump = MAX(Epump(nH(icell), TK, xH2, xHI), 0d0) ! Pumping energy in ergs
         Ebpump_prime = MAX(Epump(nH(icell), 1.001*TK, xH2, xHI), 0d0) ! Pumping energy in ergs
         ss_H2 = 1.
         if (isH2) then
            f_shd = comp_SH2(nN(ixHI), dx_SS(icell)) ! Self-shielding factor
            EUV = 0.4*(1.602d-12) ! Energy from photodissociation in ergs (Black and Dalgarno 1977)
            sigUV = sigUVthin*f_shd ! Effective UV cross section cm^2 modulated by H2 self shielding
            sigUV = sigUV*rt_c_cgs  ! [cm3/s]
            if (self_shielding) then
               ss_H2 = f_shd*comp_Sd(nN(ixHII), nN(ixHI), dx_SS(icell), dust_to_gas_mass_ratio_over_mw)
            end if
         end if

#if NGROUPS>0
         ! Set dust opacities--------------------------------------------------
         if (rt .and. nGroups .gt. 0) then
            kAbs_loc = kappaAbs
            kSc_loc = kappaSc
            if (is_kIR_T) then ! k_IR depends on T
               ! Special stuff for Krumholz/Davis experiment
               if (rt_T_rad) then  ! Use radiation temperature for kappa
                  E_rad = group_egy_erg(iIR)*dNp(iIR)
                  TR = max(T2_min_fix, (E_rad*rt_c_fraction/a_r)**0.25)
                  dT2 = TR/mu; TK = TR
               end if
               kAbs_loc(iIR) = kappaAbs(iIR)*(TK/10d0)**2
               kSc_loc(iIR) = kappaSc(iIR)*(TK/10d0)**2
            end if
            ! Set dust absorption and scattering rates [s-1]:
            if (TK .le. T_sputter) then ! dust is assumed to be destroyed if too hot
               dustAbs(:) = kAbs_loc(:)*rho(icell)*Zsolar(icell)*rt_c_cgs  ! [cm2/g * g/cm^3 * cm/s]
               dustSc(iIR) = kSc_loc(iIR)*rho(icell)*Zsolar(icell)*rt_c_cgs
               if (simple_dust_depletion) then 
                  ! The normalization is needed in order to get the data in terms of cm2 to cm2/H
                  ! This normalization was found by empirically calculating the conversion based on 
                  ! the BARE-GR-S model.  It follows closely to the value from draine but could be 10^-6
                  dustAbs(:) = signc_dust(:,1)*(nH(icell)/X)*dust_to_gas_mass_ratio_over_mw/(1784268.76d0) ! [cm2/H * H/cm^3 * cm/s]
               end if
            else
               dustAbs(:) = 0d0
               dustSc(iIR) = 0d0
            end if
         end if

         !(i) UPDATE PHOTON DENSITY AND FLUX **********************************
         if (rt .and. rt_advect) then
            recRad(1:nGroups) = 0.; phAbs(1:nGroups) = 0.
            ! Scattering rate; reduce the photon flux, but not photon density:
            phSc(1:nGroups) = 0.

            ! EMISSION FROM GAS
            if (.not. rt_OTSA .and. rt_advect) then ! ----------- Rec. radiation
               alpha(ixHII) = inp_coolrates_table(tbl_alphaA_HII, TK) &
                              - inp_coolrates_table(tbl_alphaB_HII, TK)
               ! alpha(2) A-B becomes negative around 1K, hence the max
               if (isHe) then
                  alpha(ixHeII) = MAX(0.d0, inp_coolrates_table(tbl_alphaA_HeII, TK) &
                                      - inp_coolrates_table(tbl_alphaB_HeII, TK))
                  alpha(ixHeIII) = inp_coolrates_table(tbl_alphaA_HeIII, TK) &
                                   - inp_coolrates_table(tbl_alphaB_HeIII, TK)
               end if
               if (isH2) alpha(ixHI) = 0.d0 ! no recombination?
               do iion = 1, nIons
                  if (spec2group(iion) .gt. 0) &  ! Contribution of ion -> group
                     recRad(spec2group(iion)) = &
                     recRad(spec2group(iion)) + alpha(iion)*nI(iion)*ne
               end do
            end if

            ! ABSORPTION/SCATTERING OF PHOTONS BY GAS
            do igroup = 1, nGroups      ! -------------------Ionization absorbtion
               if (igroup .eq. LWgroup) then
                  if(isH2Katz) then
                     phAbs(igroup) = nN(ixHI)*(1 + fpump)*sigUV ! For LW radiation
                  else
                     phAbs(igroup) = 0.0
                  endif
               else
                  phAbs(igroup) = SUM(nN(:)*signc(igroup, :)) ! s-1
               end if

               ! Account for photon absorption by metals
               if (Zsolar(icell).gt.min_z_neq) then
                  if (oxygen_ions) then
                     do im=1,n_oxygen_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_oxygen(igroup, im)*nOxygen(icell)*xO(im,icell)*depletion_factors(2)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (nitrogen_ions) then
                     do im=1,n_nitrogen_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_nitrogen(igroup, im)*nNitrogen(icell)*xN(im,icell)*depletion_factors(3)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (carbon_ions) then
                     do im=1,n_carbon_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_carbon(igroup, im)*nCarbon(icell)*xC(im,icell)*depletion_factors(8)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (magnesium_ions) then
                     do im=1,n_magnesium_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_magnesium(igroup, im)*nMagnesium(icell)*xMg(im,icell)*depletion_factors(4)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (silicon_ions) then
                     do im=1,n_silicon_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_silicon(igroup, im)*nSilicon(icell)*xSi(im,icell)*depletion_factors(6)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (sulfur_ions) then
                     do im=1,n_sulfur_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_sulfur(igroup, im)*nSulfur(icell)*xS(im,icell)*depletion_factors(9)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (iron_ions) then
                     do im=1,n_iron_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_iron(igroup, im)*nIron(icell)*xFe(im,icell)*depletion_factors(1)) ! cross section * c * number_density_ion
                     enddo
                  endif
                  if (neon_ions) then
                     do im=1,n_neon_ions-1
                        phAbs(igroup) = phAbs(igroup) + (signc_neon(igroup, im)*nNeon(icell)*xNe(im,icell)*depletion_factors(5)) ! cross section * c * number_density_ion
                     enddo
                  endif
               endif
            end do

            ! IR, optical and UV depletion by dust absorption: ----------------
            if (rt_isIR) & !IR scattering/abs on dust (abs after T update)
               phSc(iIR) = phSc(iIR) + dustSc(iIR)
            do igroup = 1, nGroups    ! Deplete photons, since they go into IR
               if (.not. (rt_isIR .and. igroup .eq. iIR)) &  ! IR done elsewhere
                  phAbs(igroup) = phAbs(igroup) + dustAbs(igroup)
            end do

            do igroup = 1, nGroups  ! ------------------- Do the update of N and F
               dNp(igroup) = MAX(smallNp, &
                                 (ddt(icell)*(recRad(igroup) + dNpdt(igroup, icell)) &
                                  + dNp(igroup)) &
                                 /(1.d0 + ddt(icell)*phAbs(igroup)))
               dUU = ABS(dNp(igroup) - Np(igroup, icell)) &
                     /(Np(igroup, icell) + Np_MIN)*one_over_Np_FRAC
               if (dUU .gt. 1.d0.and.igroup.gt.LWgroup) then
                  code = 1; RETURN                        ! ddt(icell) too big
               end if
               fracMax = MAX(fracMax, dUU)      ! To check if ddt can be increased

               do idim = 1, nDim
                  dFp(idim, igroup) = &
                     (ddt(icell)*dFpdt(idim, igroup, icell) + dFp(idim, igroup)) &
                     /(1d0 + ddt(icell)*(phAbs(igroup) + phSc(igroup)))
               end do
               call reduce_flux(dFp(:, igroup), dNp(igroup)*rt_c_cgs)

               do idim = 1, nDim
                  dUU = ABS(dFp(idim, igroup) - Fp(idim, igroup, icell)) &
                        /(ABS(Fp(idim, igroup, icell)) + Fp_MIN)*one_over_Fp_FRAC
                  if (dUU .gt. 1.d0.and.igroup.gt.LWgroup) then
                     code = 2; RETURN                     ! ddt(icell) too big
                  end if
                  fracMax = MAX(fracMax, dUU)   ! To check if ddt can be increased
               end do

            end do

            dmom(1:nDim) = 0d0
            do igroup = 1, nGroups ! -------Momentum transfer from photons to gas:
               mom_fact = ddt(icell)*(phAbs(igroup) + phSc(igroup)) &
                          *group_egy_erg(igroup)*one_over_c_cgs

               if (rt_isoPress .and. .not. (rt_isIR .and. igroup == iIR)) then
                  ! rt_isoPress: assume f=1, where f is reduced flux.
                  fluxMag = sqrt(sum((dFp(:, igroup))**2))
                  if (fluxMag .gt. 0d0) then
                     mom_fact = mom_fact*dNp(igroup)/fluxMag
                  else
                     mom_fact = 0d0
                  end if
               else
                  mom_fact = mom_fact*one_over_rt_c_cgs
               end if

               do idim = 1, nDim
                  dmom(idim) = dmom(idim) + dFp(idim, igroup)*mom_fact
               end do
            end do
            dp_gas = dp_gas + dmom*rt_pressBoost        ! update gas momentum

            ! Add absorbed UV/optical energy to IR:----------------------------
            if (rt_isIR) then
               do igroup = iIR + 1, nGroups

                  ! H2 fluorescence
                  if (isH2Katz .and. igroup .eq. LWgroup) then
                     kUV = dNp(igroup)*sigUV
                     HrateLW = (((kUV*fpump*Ebpump) + (kUV*EUV))*nN(ixHI)) ! [erg/s/cm3] = [#/s]*[erg/#]*[cm-3]
                     HrateLW = max(HrateLW, 0.0d0)
                     dNpLW = HrateLW*ddt(icell)/group_egy_erg(igroup) ! energy used to heat up the gas through photodissociation and UV pumping
                     phAbsLW = nN(ixHI)*(1 + fpump)*sigUV ! [s-1] = [cm-3]*[cm3/s]
                     dNpAdd = phAbsLW*ddt(icell)*dNp(igroup) - dNpLW
                     dNpAdd = max(dNpAdd, 0.0d0)
                     dNp(iIR) = dNp(iIR) + dNpAdd*group_egy_ratio(igroup)
                  end if

                  ! photo-electic heating on dust by UV
                  if (PEH .and. group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                     f_dg = dust_to_gas_mass_ratio_over_mw
                     F_pe = dNp(igroup)*rt_c_cgs         ! [#/cm2/s]
                     F_pe = F_pe*group_egy_erg(igroup) ! [erg/cm2/s]
                     G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
                     G_0 = F_pe/G_conv
                     eps_PE = PE_efficiency(G_0, TK, ne)
                     ! The extra factor of 1.5 comes from the fact that the zubko model has 1.5 more C atoms locked
                     ! up in dust compared to what was assumed in Wolfire 2003
                     HratePE = 1.5d0*1.3d-24*eps_PE*G_0*f_dg*nH(icell) ! [erg/cm3/s] energy used to heat up the gas through PEH
                     dNpPE = HratePE*ddt(icell)/group_egy_erg(igroup) ! [#/cm3]
                     dNpAdd = dustAbs(igroup)*ddt(icell)*dNp(igroup) - dNpPE ! [#/cm3]
                     dNpAdd = max(dNpAdd, 0.0d0)
                     dNp(iIR) = dNp(iIR) + dNpAdd*group_egy_ratio(igroup)
                  else
                     dNp(iIR) = dNp(iIR) + dustAbs(igroup)*ddt(icell) &
                                *dNp(igroup)*group_egy_ratio(igroup)
                  end if
               end do
            end if
            ! -----------------------------------------------------------------
         end if !if(rt)
#endif
         !(ii) UPDATE TEMPERATURE *********************************************
         if (c_switch(icell) .and. cooling .and. .not. rt_T_rad) then
            Hrate = 0.                             !  Heating rate [erg cm-3 s-1]
            H2heating = 0.; H2heating_prime = 0.

            ! Initialize photoelectric heating rate and derivative
            HratePE_tot = 0.
            HratePE_prime = 0. 
            HratePE_prime_a = 0.
            HratePE_prime_b = 0.
            dust_recomb_cool = 0.
            dust_recomb_cool_prime = 0.
            ct_heat_cool = 0.
            ct_heat_cool_prime = 0.

            if (rt .and. rt_advect .and. PHlocal) then
               do igroup = 1, nGroups                              !  Photoheating
                  if (igroup .eq. LWgroup) then
                     if(isH2Katz) then
                        ! 1) photodisscodiation of an H2 molecule due to LW radiation deposits EUV as heat
                        ! 2) UV photons can also heat the gas by indirectly exciting the vibrational levels of the H2 (UV pumping rate)
                        kUV = dNp(igroup)*sigUV  ! #/cm3*cm3/s = #/s-1
                        HrateLW = ((kUV*fpump*Ebpump) + (kUV*EUV))*nN(ixHI) ! [erg/s/cm3] = [#/s]*[erg/#]*[cm-3]
                        HrateLW_prime = ((kUV*fpump*Ebpump_prime) + (kUV*EUV))*nN(ixHI) ! [erg/s/cm3] = [#/s]*[erg/#]*[cm-3]
                        H2heating = H2heating + max(HrateLW, 0.0d0)
                        H2heating_prime = H2heating_prime + HrateLW_prime
                     endif
                  else
                     Hrate = Hrate + dNp(igroup)*SUM(nN(:)*PHrate(igroup, :))

                     if (oxygen_ions)    Hrate = Hrate + ( dNp(igroup) * depletion_factors(2) * nOxygen(icell)   * SUM( xO(:,icell)  * PHrate_oxygen(igroup,:)   ) )
                     if (nitrogen_ions)  Hrate = Hrate + ( dNp(igroup) * depletion_factors(3) * nNitrogen(icell) * SUM( xN(:,icell)  * PHrate_nitrogen(igroup,:) ) )
                     if (carbon_ions)    Hrate = Hrate + ( dNp(igroup) * depletion_factors(8) * nCarbon(icell)   * SUM( xC(:,icell)  * PHrate_carbon(igroup,:)   ) )
                     if (magnesium_ions) Hrate = Hrate + ( dNp(igroup) * depletion_factors(4) * nMagnesium(icell)* SUM( xMg(:,icell) * PHrate_magnesium(igroup,:)) )
                     if (silicon_ions)   Hrate = Hrate + ( dNp(igroup) * depletion_factors(6) * nSilicon(icell)  * SUM( xSi(:,icell) * PHrate_silicon(igroup,:)  ) )
                     if (sulfur_ions)    Hrate = Hrate + ( dNp(igroup) * depletion_factors(9) * nSulfur(icell)   * SUM( xS(:,icell)  * PHrate_sulfur(igroup,:)   ) )
                     if (iron_ions)      Hrate = Hrate + ( dNp(igroup) * depletion_factors(1) * nIron(icell)     * SUM( xFe(:,icell) * PHrate_iron(igroup,:)     ) )
                     if (neon_ions)      Hrate = Hrate + ( dNp(igroup) * depletion_factors(5) * nNeon(icell)     * SUM( xNe(:,icell) * PHrate_neon(igroup,:)     ) )
                  end if

                  ! photo-electic heating on dust by UV
                  if (PEH) then
                     if (group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                        f_dg = dust_to_gas_mass_ratio_over_mw
                        F_pe = dNp(igroup)*rt_c_cgs         ! [#/cm2/s]
                        F_pe = F_pe*group_egy_erg(igroup) ! [erg/cm2/s]
                        G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
                        G_0 = (F_pe/G_conv)
                        eps_PE = PE_efficiency(G_0, TK, ne)
                        HratePE = 1.5d0*1.3d-24*eps_PE*G_0*f_dg*nH(icell) ! [erg/cm3/s]
                        HratePE_tot = HratePE_tot + HratePE
                        Hrate = Hrate + HratePE
                        eps_PE = PE_efficiency(G_0, TK*1.001, ne)
                        HratePE_prime_a = 1.5d0*1.3d-24*eps_PE*G_0*f_dg*nH(icell)
                        HratePE_prime_a = (HratePE - HratePE_prime_a) / (TK - (TK*1.001))
                        HratePE_prime = HratePE_prime_a + HratePE_prime_b

                        ! Account for dust recombination cooling
                        beta_drc = 0.74d0 / (TK**0.068d0)
                        ! extra factor of 1.5 comes from the difference in PAH abundance
                        dust_recomb_cool_loc = (1.5d0*4.65d-30) * (TK**0.94d0) * ( (G_0 * SQRT(TK) / (0.5d0*ne))**beta_drc ) * ne * 0.5d0 * f_dg * nH(icell)

                        beta_drc = 0.74d0 / ((TK*1.001)**0.068d0)
                        dust_recomb_cool_prime_loc = (1.5d0*4.65d-30) * ((TK*1.001)**0.94d0) * ( (G_0 * SQRT((TK*1.001)) / (0.5d0*ne))**beta_drc ) * ne * 0.5d0 * f_dg * nH(icell)
                        dust_recomb_cool_prime_loc = (dust_recomb_cool_loc - dust_recomb_cool_prime_loc) / (TK - (TK*1.001))

                        dust_recomb_cool = dust_recomb_cool + dust_recomb_cool_loc
                        dust_recomb_cool_prime = dust_recomb_cool_prime + dust_recomb_cool_prime_loc
                     end if
                  end if

               end do
            end if

            ! PEH contribution from the UVBG
            if (PEH.and.UV_background_G0.gt.0.d0) then
               f_dg = dust_to_gas_mass_ratio_over_mw
               G_0 = UV_background_G0
               eps_PE = PE_efficiency(G_0, TK, ne)
               HratePE = 1.5d0*1.3d-24*eps_PE*G_0*f_dg*nH(icell) ! [erg/cm3/s]
               HratePE_tot = HratePE_tot + HratePE
               Hrate = Hrate + HratePE
               eps_PE = PE_efficiency(G_0, TK*1.001, ne)
               HratePE_prime_b = 1.5d0*1.3d-24*eps_PE*G_0*f_dg*nH(icell)
               HratePE_prime_b = (HratePE - HratePE_prime_b) / (TK - (TK*1.001))
               HratePE_prime = HratePE_prime_a + HratePE_prime_b

               ! Account for dust recombination cooling
               beta_drc = 0.74d0 / (TK**0.068d0)
               dust_recomb_cool_loc = (1.5d0*4.65d-30) * (TK**0.94d0) * ( (G_0 * SQRT(TK) / (0.5d0*ne))**beta_drc ) * ne * 0.5d0 * f_dg * nH(icell)

               beta_drc = 0.74d0 / ((TK*1.001)**0.068d0)
               dust_recomb_cool_prime_loc = (1.5d0*4.65d-30) * ((TK*1.001)**0.94d0) * ( (G_0 * SQRT((TK*1.001)) / (0.5d0*ne))**beta_drc ) * ne * 0.5d0 * f_dg * nH(icell)
               dust_recomb_cool_prime_loc = (dust_recomb_cool_loc - dust_recomb_cool_prime_loc) / (TK - (TK*1.001))

               dust_recomb_cool = dust_recomb_cool + dust_recomb_cool_loc
               dust_recomb_cool_prime = dust_recomb_cool_prime + dust_recomb_cool_prime_loc
            end if

            cool_dust = 0.d0
            cool_dust_prime = 0.d0
            if (include_dust_heat_cool .and. TK .le. T_sputter .and. simple_dust_depletion) then 
               ! Heating and cooling from dust from Glover & Jappsen 2007
               ! modified by our dust to gas mass ratio
               G_0 = 0.d0
               if (rt_advect) then 
                  G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
                  do igroup = 1, nGroups  
                     if (group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                        G_0 = G_0 + (F_pe/G_conv)
                     end if
                  end do
               end if 
               G_0 = G_0 + UV_background_G0

               T_dust = 16.4d0 * (1.7d0 * G_0)**(1.d0/6.d0)
               T_dust = MAX( T_dust, 2.725 * ( (1.d0/aexp) - 1.d0 ) ) ! Limit dust temp minimum to CMB temp

               dust_hc_const = 1.0d-33 ! For atomic dominated regions
               if (isH2) then
                  if (xH2.gt.xHI) dust_hc_const = 3.8d-33
               end if

               ! Dust recombination cooling from bialy 2019
               cool_dust = dust_hc_const * SQRT(TK) * (TK - T_dust) * ( 1.d0 - ( 0.8d0 * EXP(-75.d0/TK) ) )
               cool_dust_prime = dust_hc_const * SQRT(TK*1.001d0) * ((TK*1.001d0) - T_dust) * ( 1.d0 - ( 0.8d0 * EXP(-75.d0/(TK*1.001d0)) ) )

               cool_dust = cool_dust * dust_to_gas_mass_ratio_over_mw * (nH(icell)**2.d0)
               cool_dust_prime = cool_dust_prime * dust_to_gas_mass_ratio_over_mw * (nH(icell)**2.d0)
               cool_dust_prime = (cool_dust - cool_dust_prime) / (TK - (TK*1.001))
            end if

            ! initialize cosmic ray heating rates
            HI_cosmic_ray_ionization_rate = 0.d0
            He_cosmic_ray_ionization_rate = 0.d0
            H2_cosmic_ray_ionization_rate = 0.d0
            HI_cosmic_ray_ionization_rate_primary = 0.d0
            He_cosmic_ray_ionization_rate_primary = 0.d0
            H2_cosmic_ray_ionization_rate_primary = 0.d0
            if (include_cosmic_ray_heating.or.include_cosmic_ray_ionization) then
               effective_n = rho(icell) / mH ! Gas density assuming everything is hydrogen
               xe = ne / effective_n

               q_cr = 6.43d0 * ( 1.d0 + 4.06d0 * SQRT( xe / (0.07d0 + xe) ) ) * ev2erg ! Bialy 2019
               phi_s = (1.d0 - (xe / 1.2d0)) * (0.67d0 / (1.d0 + (xe / 0.05d0)))

               HI_cosmic_ray_ionization_rate = 1.d0 * xi_h_cr * (1.d0 + phi_s)
               HI_cosmic_ray_ionization_rate_primary = 1.d0 * xi_h_cr

               if (isHe) then 
                  He_cosmic_ray_ionization_rate = 1.1d0 * xi_h_cr * (1.d0 + phi_s)
                  He_cosmic_ray_ionization_rate_primary = 1.1d0 * xi_h_cr
               end if

               if (isH2) then 
                  H2_cosmic_ray_ionization_rate = 2.d0 * xi_h_cr * (1.d0 + phi_s)
                  H2_cosmic_ray_ionization_rate_primary = 2.d0 * xi_h_cr
               end if

            end if

            cr_heating_tot = 0.d0
            if (include_cosmic_ray_heating) then 
               ! Bialy 2019
               cr_heating_tot = cr_heating_tot + ( nH(icell) * xHI * HI_cosmic_ray_ionization_rate_primary * q_cr )

               if (isHe) then 
                  cr_heating_tot = cr_heating_tot + ( nN(ixHeII) * He_cosmic_ray_ionization_rate_primary * q_cr )
               end if

               if (isH2) then    
                  cr_heating_tot = cr_heating_tot + ( nH(icell) * xH2 * H2_cosmic_ray_ionization_rate_primary * q_cr )
               end if

               ! Now consider heating of the electrons from cosmic rays -- Bialy 2019
               cr_heating_tot = cr_heating_tot + (ne * xi_h_cr * 287.d0*ev2erg)
            end if
            Hrate = Hrate + cr_heating_tot
            
            ! Heating from H2 destruction from UV background
            if (isH2Katz) then
               ! 1) photodisscodiation of an H2 molecule due to LW radiation deposits EUV as heat
               ! 2) UV photons can also heat the gas by indirectly exciting the vibrational levels of the H2 (UV pumping rate)
               kUV = 0.d0
               kUV = kUV + UV_background_h2*ss_factor
               kUV = kUV + UV_background_G0*5.8d-11 ! See Koyama 2000
               HrateLW = ((kUV*fpump*Ebpump) + (kUV*EUV))*nN(ixHI) ! [erg/s/cm3] = [#/s]*[erg/#]*[cm-3]
               HrateLW_prime = ((kUV*fpump*Ebpump_prime) + (kUV*EUV))*nN(ixHI) ! [erg/s/cm3] = [#/s]*[erg/#]*[cm-3]
               H2heating = H2heating + max(HrateLW, 0.0d0)
               H2heating_prime = H2heating_prime + max(HrateLW_prime, 0.0d0)
               
               G_0 = 0.d0
               if (rt_advect) then 
                  G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
                  do igroup = 1, nGroups  
                     if (group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                        G_0 = G_0 + (F_pe/G_conv)
                     end if
                  end do
               end if 
               G_0 = G_0 + UV_background_G0

               ! H2 formation heating rates --> See Rollig 2006
               effective_n = rho(icell) / mH ! Gas density assuming everything is hydrogen
               xe = ne / effective_n
               H2heating = H2heating + (2.4d-12)*comp_Alpha_H2(TK, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, G_0, xHII, xHI)*xHI*nH(icell)*nH(icell) ! [cm3 s-1]
               H2heating_prime = H2heating_prime + (2.4d-12)*comp_Alpha_H2(1.001*TK, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, G_0, xHII, xHI)*xHI*nH(icell)*nH(icell) ! [cm3 s-1]
            endif

            if (haardt_madau) Hrate = Hrate + SUM(nN(:)*UVrates(1:nIons, 2))*ss_factor

            ! Harley UV background heating
            Hrate = Hrate + (nN(ixHII)*UV_background_hydrogen_heating*ss_factor) !HI
            if (isHe) then 
               Hrate = Hrate + (nN(ixHeII)*UV_background_helium_heating(1)*ss_factor) !HeI
               Hrate = Hrate + (nN(ixHeIII)*UV_background_helium_heating(2)*ss_factor) !HeII
            endif
            if (isH2Katz) Hrate = Hrate + (nN(ixHI)*UV_background_H2_heating*ss_factor) !H2

            ! Add heating from UV background and cosmic rays
            if (oxygen_ions) then
               do im=1,n_oxygen_ions-1
                  Hrate = Hrate + (nOxygen(icell)*xO(im,icell)*UV_background_oxygen_heating(im)*ss_factor*depletion_factors(2))
                  if (include_cosmic_ray_heating) then 
                     Hrate = Hrate + (nOxygen(icell)*xO(im,icell)*xi_h_cr*depletion_factors(2)*q_cr*cr_ionization_oxygen(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nOxygen(icell)*xO(im,icell)*depletion_factors(2)*(cr_ionization_ind_UV_OI * (1.72d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if
                  end if
               enddo
            endif
            if (nitrogen_ions) then
               do im=1,n_nitrogen_ions-1
                  Hrate = Hrate + (nNitrogen(icell)*xN(im,icell)*UV_background_nitrogen_heating(im)*ss_factor*depletion_factors(3))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nNitrogen(icell)*xN(im,icell)*xi_h_cr*depletion_factors(3)*q_cr*cr_ionization_nitrogen(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nNitrogen(icell)*xN(im,icell)*depletion_factors(3)*(cr_ionization_ind_UV_NI * (1.17d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if
                  end if
               enddo
            endif
            if (carbon_ions) then
               do im=1,n_carbon_ions-1
                  Hrate = Hrate + (nCarbon(icell)*xC(im,icell)*UV_background_carbon_heating(im)*ss_factor*depletion_factors(8))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nCarbon(icell)*xC(im,icell)*xi_h_cr*depletion_factors(8)*q_cr*cr_ionization_carbon(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nCarbon(icell)*xC(im,icell)*depletion_factors(8)*(cr_ionization_ind_UV_CI * (0.96d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if
                  end if
                  if (UV_background_G0.gt.0.d0) then
                     if (im.eq.1) Hrate = Hrate + (nCarbon(icell)*xC(im,icell)*depletion_factors(8)*UV_background_G0*(3.39d-10)*(0.8399d0)*ev2erg)
                  end if
               enddo
            endif
            if (magnesium_ions) then
               do im=1,n_magnesium_ions-1
                  Hrate = Hrate + (nMagnesium(icell)*xMg(im,icell)*UV_background_magnesium_heating(im)*ss_factor*depletion_factors(4))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nMagnesium(icell)*xMg(im,icell)*xi_h_cr*depletion_factors(4)*q_cr*cr_ionization_magnesium(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nMagnesium(icell)*xMg(im,icell)*depletion_factors(4)*(cr_ionization_ind_UV_MgI * (2.55d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if 
                  end if
                  if (UV_background_G0.gt.0.d0) then
                     if (im.eq.1) Hrate = Hrate + (nMagnesium(icell)*xMg(im,icell)*depletion_factors(4)*UV_background_G0*(6.59d-11)*(2.0113d0)*ev2erg)
                  end if
               enddo
            endif
            if (silicon_ions) then
               do im=1,n_silicon_ions-1
                  Hrate = Hrate + (nSilicon(icell)*xSi(im,icell)*UV_background_silicon_heating(im)*ss_factor*depletion_factors(6))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nSilicon(icell)*xSi(im,icell)*xi_h_cr*depletion_factors(6)*q_cr*cr_ionization_silicon(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nSilicon(icell)*xSi(im,icell)*depletion_factors(6)*(cr_ionization_ind_UV_SiI * (2.08d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if
                  end if
                  if (UV_background_G0.gt.0.d0) then
                     if (im.eq.1) Hrate = Hrate + (nSilicon(icell)*xSi(im,icell)*depletion_factors(6)*UV_background_G0*(4.47d-9)*(1.840d0)*ev2erg)
                  end if
               enddo
            endif
            if (sulfur_ions) then
               do im=1,n_sulfur_ions-1
                  Hrate = Hrate + (nSulfur(icell)*xS(im,icell)*UV_background_sulfur_heating(im)*ss_factor*depletion_factors(9))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nSulfur(icell)*xS(im,icell)*xi_h_cr*depletion_factors(9)*q_cr*cr_ionization_sulfur(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nSulfur(icell)*xS(im,icell)*depletion_factors(9)*(cr_ionization_ind_UV_SI * (1.32d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if
                  end if
                  if (UV_background_G0.gt.0.d0) then
                     if (im.eq.1) Hrate = Hrate + (nSulfur(icell)*xS(im,icell)*depletion_factors(9)*UV_background_G0*(1.13d-9)*(1.126d0)*ev2erg)
                  end if
               enddo
            endif
            if (iron_ions) then
               do im=1,n_iron_ions-1
                  Hrate = Hrate + (nIron(icell)*xFe(im,icell)*UV_background_iron_heating(im)*ss_factor*depletion_factors(1))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nIron(icell)*xFe(im,icell)*xi_h_cr*depletion_factors(1)*q_cr*cr_ionization_iron(im))
                     if (im.eq.1) then 
                        Hrate = Hrate + (nIron(icell)*xFe(im,icell)*depletion_factors(1)*(cr_ionization_ind_UV_FeI * (2.32d0*ev2erg) * (H2_cosmic_ray_ionization_rate / 1.d-16)))
                     end if
                  end if
                  if (UV_background_G0.gt.0.d0) then
                     if (im.eq.1) Hrate = Hrate + (nIron(icell)*xFe(im,icell)*depletion_factors(1)*UV_background_G0*(4.71d-10)*(1.924d0)*ev2erg)
                  end if
               enddo
            endif
            if (neon_ions) then
               do im=1,n_neon_ions-1
                  Hrate = Hrate + (nNeon(icell)*xNe(im,icell)*UV_background_neon_heating(im)*ss_factor*depletion_factors(5))
                  if (include_cosmic_ray_heating) then
                     Hrate = Hrate + (nNeon(icell)*xNe(im,icell)*xi_h_cr*depletion_factors(5)*q_cr*cr_ionization_neon(im))
                     ! Note: no secondary rate for Neon due to high ionization potential
                     ! No G0 ionization for the same reason
                  end if
               enddo
            endif

            Crate = compCoolrate(TK, ne, nH(icell), nN, nI, aexp, dCdT2) ! Cooling
            dCdT2 = dCdT2*mu                            ! dC/dT2 = mu * dC/dT
            metal_tot = 1.d-40; metal_prime = 0.d0             ! Metal cooling
            fire_tot = 1.d-40; fire_prime = 0.d0             ! Hopkins cooling (metals,dust,cosmic rays, photo-electric
            fine_structure = 1.d-40; fine_structure_prime = 0.d0 ! non-equilibrium fine-structure cooling
            fsc_OI = 1.d-40; fsc_OIII = 1.d-40
            fsc_NII = 1.d-40;
            fsc_CI = 1.d-40; fsc_CII = 1.d-40
            fsc_SiI = 1.d-40; fsc_SiII = 1.d-40
            fsc_FeI = 1.d-40; fsc_FeII = 1.d-40
            fsc_SI = 1.d-40;
            fsc_NeII = 1.d-40;
            ! metal cool smoothing parameters
            metal_cool_smooth_f1 = 0.5d0 * (TANH( (5.d-3) * ( TK - 1.d4 ) ) + 1.d0 )
            metal_cool_smooth_f2 = 0.5d0 * (TANH( (5.d-3) * ( (-1.d0 * TK) + 1.d4 ) ) + 1.d0 )
            !Keep in mind that we still call this even if using fire cooling
            !This is so that we can have metal line cooling at T>10^4 from cloudy
            if (Zsolar(icell) .gt. 1.d-8) then
               if (cloudy_metal_cooling) then
                  call rt_cmp_metals_cloudy(T2(icell), nH(icell), mu, metal_tot &
                                            , metal_prime, a_exp)
                  metal_tot = metal_tot * Zsolar(icell)
                  metal_prime = metal_prime * Zsolar(icell)
               else if (os13_metal_cooling) then
                  !if (TK.gt.1.d4) then 
                  if (TK.gt.9.d3) then 
                     call cmp_metals_os13( MAX(TK,1.d4), nH(icell), ne, & 
                                          nIron(icell) * depletion_factors(1), & 
                                          nOxygen(icell) * depletion_factors(2), & 
                                          nNitrogen(icell) * depletion_factors(3), & 
                                          nMagnesium(icell) * depletion_factors(4), & 
                                          nNeon(icell) * depletion_factors(5), & 
                                          nSilicon(icell) * depletion_factors(6), &
                                          nCarbon(icell) * depletion_factors(8), & 
                                          nSulfur(icell) * depletion_factors(9), & 
                                          xFe(:,icell),xO(:,icell),xN(:,icell), &
                                          xMg(:,icell),xNe(:,icell),xSi(:,icell), &
                                          xC(:,icell),xS(:,icell),metal_tot) 
                     call cmp_metals_os13( MAX(1.001*TK,1.001d4), nH(icell), ne, & 
                                          nIron(icell) * depletion_factors(1), &
                                          nOxygen(icell) * depletion_factors(2), &
                                          nNitrogen(icell) * depletion_factors(3), &
                                          nMagnesium(icell) * depletion_factors(4), &
                                          nNeon(icell) * depletion_factors(5), &
                                          nSilicon(icell) * depletion_factors(6), &
                                          nCarbon(icell) * depletion_factors(8), &
                                          nSulfur(icell) * depletion_factors(9), &
                                          xFe(:,icell),xO(:,icell),xN(:,icell), &
                                          xMg(:,icell),xNe(:,icell),xSi(:,icell), &
                                          xC(:,icell),xS(:,icell),metal_prime) 
                     metal_tot = metal_tot * metal_cool_smooth_f1
                     metal_prime = metal_prime * metal_cool_smooth_f1
                     metal_prime = (metal_tot - metal_prime) / (MAX(TK,1.d4) - MAX(1.001*TK,1.001d4))
                  end if
               else if (htmc_metal_cooling) then
                  if (TK.gt.1.d3) then
                     call cmp_metals_ht_harley( TK, ne, & 
                                          nIron(icell) * depletion_factors(1), & 
                                          nOxygen(icell) * depletion_factors(2), & 
                                          nNitrogen(icell) * depletion_factors(3), & 
                                          nMagnesium(icell) * depletion_factors(4), & 
                                          nNeon(icell) * depletion_factors(5), & 
                                          nSilicon(icell) * depletion_factors(6), &
                                          nCarbon(icell) * depletion_factors(8), & 
                                          nSulfur(icell) * depletion_factors(9), & 
                                          xFe(:,icell),xO(:,icell),xN(:,icell), &
                                          xMg(:,icell),xNe(:,icell),xSi(:,icell), &
                                          xC(:,icell),xS(:,icell),metal_tot,.false., &
                                          metal_cool_smooth_f1)
                     call cmp_metals_ht_harley( 1.001d0*TK, ne, & 
                                          nIron(icell) * depletion_factors(1), &
                                          nOxygen(icell) * depletion_factors(2), &
                                          nNitrogen(icell) * depletion_factors(3), &
                                          nMagnesium(icell) * depletion_factors(4), &
                                          nNeon(icell) * depletion_factors(5), &
                                          nSilicon(icell) * depletion_factors(6), &
                                          nCarbon(icell) * depletion_factors(8), &
                                          nSulfur(icell) * depletion_factors(9), &
                                          xFe(:,icell),xO(:,icell),xN(:,icell), &
                                          xMg(:,icell),xNe(:,icell),xSi(:,icell), &
                                          xC(:,icell),xS(:,icell),metal_prime,.false., &
                                          metal_cool_smooth_f1)
                     !metal_tot = metal_tot * metal_cool_smooth_f1
                     !metal_prime = metal_prime * metal_cool_smooth_f1
                     metal_prime = (metal_tot - metal_prime) / (MAX(TK,9.d3) - MAX(1.001*TK,9.d3*1.001d0))
                     if (metal_tot.gt.1.d0) then 
                        write(*,*) "Metal_tot is way broken"
                        ! call it again with logs
                        call cmp_metals_ht_harley(TK, ne, &
                                          nIron(icell) * depletion_factors(1), &
                                          nOxygen(icell) * depletion_factors(2), &
                                          nNitrogen(icell) * depletion_factors(3), &
                                          nMagnesium(icell) * depletion_factors(4), &
                                          nNeon(icell) * depletion_factors(5), &
                                          nSilicon(icell) * depletion_factors(6), &
                                          nCarbon(icell) * depletion_factors(8), &
                                          nSulfur(icell) * depletion_factors(9), &
                                          xFe(:,icell),xO(:,icell),xN(:,icell), &
                                          xMg(:,icell),xNe(:,icell),xSi(:,icell), &
                                          xC(:,icell),xS(:,icell),metal_tot,.true., &
                                          metal_cool_smooth_f1) 
                     end if
                     if (metal_prime.gt.1.d0) then
                        write(*,*) "Metal_prime is way broken"
                        call cmp_metals_ht_harley( 1.001d0*TK, ne, &
                                          nIron(icell) * depletion_factors(1), &
                                          nOxygen(icell) * depletion_factors(2), &
                                          nNitrogen(icell) * depletion_factors(3), &
                                          nMagnesium(icell) * depletion_factors(4), &
                                          nNeon(icell) * depletion_factors(5), &
                                          nSilicon(icell) * depletion_factors(6), &
                                          nCarbon(icell) * depletion_factors(8), &
                                          nSulfur(icell) * depletion_factors(9), &
                                          xFe(:,icell),xO(:,icell),xN(:,icell), &
                                          xMg(:,icell),xNe(:,icell),xSi(:,icell), &
                                          xC(:,icell),xS(:,icell),metal_prime,.true., &
                                          metal_cool_smooth_f1) 
                     end if
                  end if
               else
                  call rt_cmp_metals(T2(icell), nH(icell), mu, metal_tot &
                                     , metal_prime, a_exp)
                  metal_tot = metal_tot * Zsolar(icell)
                  metal_prime = metal_prime * Zsolar(icell)
               end if
            end if

            ! low temperature non-equilibrium fine structure cooling
            if (low_t_noneq_cooling.and.Zsolar(icell).gt.min_z_neq) then

               nH2_cool = 0.d0
               if (isH2) nH2_cool = nN(ixHI)

               if (TK.le.1.1d4) then
                  if (oxygen_ions) then
                      if ((nOxygen(icell)*xO(1,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                         fsc_OI = OI_fine_structure(TK,nOxygen(icell)*xO(1,icell)*depletion_factors(2),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                         fine_structure = fine_structure + fsc_OI
                         fine_structure_prime = fine_structure_prime + OI_fine_structure(TK*1.001d0,nOxygen(icell)*xO(1,icell)*depletion_factors(2),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                      endif
                  endif

                  if (carbon_ions) then
                      if ((nCarbon(icell)*xC(1,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                         fsc_CI = CI_fine_structure(TK,nCarbon(icell)*xC(1,icell)*depletion_factors(8),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                         fine_structure = fine_structure + fsc_CI
                         fine_structure_prime = fine_structure_prime + CI_fine_structure(TK*1.001d0,nCarbon(icell)*xC(1,icell)*depletion_factors(8),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                      endif
#if NCARBONIONS>1
                      if ((nCarbon(icell)*xC(2,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                         fsc_CII = CII_fine_structure(TK,nCarbon(icell)*xC(2,icell)*depletion_factors(8),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                         fine_structure = fine_structure + fsc_CII
                         fine_structure_prime = fine_structure_prime + CII_fine_structure(TK*1.001d0,nCarbon(icell)*xC(2,icell)*depletion_factors(8),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                      endif
#endif
                  endif

                  if (silicon_ions) then
                      if ((nSilicon(icell)*xSi(1,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                          fsc_SiI = SiI_fine_structure(TK,nSilicon(icell)*xSi(1,icell)*depletion_factors(6),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                          fine_structure = fine_structure + fsc_SiI
                          fine_structure_prime = fine_structure_prime + SiI_fine_structure(TK*1.001d0,nSilicon(icell)*xSi(1,icell)*depletion_factors(6),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                       endif
#if NSILICONIONS>1
                      if ((nSilicon(icell)*xSi(2,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                         fsc_SiII = SiII_fine_structure(TK,nSilicon(icell)*xSi(2,icell)*depletion_factors(6),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                         fine_structure = fine_structure + fsc_SiII
                         fine_structure_prime = fine_structure_prime + SiII_fine_structure(TK*1.001d0,nSilicon(icell)*xSi(2,icell)*depletion_factors(6),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                      endif
#endif
                  endif

                  if (iron_ions) then
                      if ((nIron(icell)*xFe(1,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                         fsc_FeI = FeI_fine_structure(TK,nIron(icell)*xFe(1,icell)*depletion_factors(1),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                         fine_structure = fine_structure + fsc_FeI
                         fine_structure_prime = fine_structure_prime + FeI_fine_structure(TK*1.001d0,nIron(icell)*xFe(1,icell)*depletion_factors(1),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                      endif
#if NIRONIONS>1
                      if ((nIron(icell)*xFe(2,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                         fsc_FeII = FeII_fine_structure(TK,nIron(icell)*xFe(2,icell)*depletion_factors(1),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                         fine_structure = fine_structure + fsc_FeII
                         fine_structure_prime = fine_structure_prime + FeII_fine_structure(TK*1.001d0,nIron(icell)*xFe(2,icell)*depletion_factors(1),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                      endif
#endif
                  endif

                  if (sulfur_ions) then
                     if ((nSulfur(icell)*xS(1,icell)).gt.min_cool_ion) then ! No reason to call in very low densities
                        fsc_SI = SI_fine_structure(TK,nSulfur(icell)*xS(1,icell)*depletion_factors(9),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                        fine_structure = fine_structure + fsc_SI
                        fine_structure_prime = fine_structure_prime + SI_fine_structure(TK*1.001d0,nSulfur(icell)*xS(1,icell)*depletion_factors(9),nN(ixHII),nI(ixHII),ne,nH2_cool,nN(ixHeII),nI(ixHeII),nI(ixHeIII),(1.d0/a_exp)-1.d0)*metal_cool_smooth_f2
                     endif
                  endif

                  ! exponential decrease at low density for stability
                  if (metal_cooling_catastrophe_fix) then 
                     fine_structure = fine_structure * exp(metal_cooling_catastrophe_rho/(-rho(icell)*X/mH))
                     fine_structure_prime = fine_structure_prime * exp(metal_cooling_catastrophe_rho/(-rho(icell)*X/mH))
                  endif
                  fine_structure_prime = (fine_structure - fine_structure_prime) / (TK - (1.001d0*TK))
               endif
            end if

            ! Charge transfer heating and cooling
            if (include_ct_heat_cool.and.TK.le.(1.d5/1.001d0)) then
               ! Helium
               if (isHe) then 
                  ct_heat_cool = ct_heat_cool + (HEIIRecomb(Tk) * nN(ixHII) * nI(ixHeII) * 10.99d0 * ev2erg)  ! H + He+ -> He + H+ 
                  ct_heat_cool_prime = ct_heat_cool_prime + (HEIIRecomb(Tk*1.001d0) * nN(ixHII) * nI(ixHeII) * 10.99d0 * ev2erg)  ! H + He+ -> He + H+    
               end if

               ! Carbon
               if (carbon_ions) then
                  ct_heat_cool = ct_heat_cool + (HCTIon(1,6,Tk) * nI(ixHII) * nCarbon(icell) * xC(1,icell) * depletion_factors(8) * (2.34d0) * ev2erg) ! H+ + C -> C+ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(1,6,Tk*1.001d0) * nI(ixHII) * nCarbon(icell) * xC(1,icell) * depletion_factors(8) * (2.34d0) * ev2erg)
#if NCARBONIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTRecom(2,6,Tk) * nN(ixHII) * nCarbon(icell) * xC(2,icell) * depletion_factors(8) * (-2.34d0) * ev2erg) ! H + C+ -> C + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(2,6,Tk*1.001d0) * nN(ixHII) * nCarbon(icell) * xC(2,icell) * depletion_factors(8) * (-2.34d0) * ev2erg)
#endif 
#if NCARBONIONS>2
                  ct_heat_cool = ct_heat_cool + (HCTRecom(3,6,Tk) * nN(ixHII) * nCarbon(icell) * xC(3,icell) * depletion_factors(8) * (4.01d0) * ev2erg) ! H + C++ -> C+ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(3,6,Tk*1.001d0) * nN(ixHII) * nCarbon(icell) * xC(3,icell) * depletion_factors(8) * (4.01d0) * ev2erg)
#endif 
#if NCARBONIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,6,Tk) * nN(ixHII) * nCarbon(icell) * xC(4,icell) * depletion_factors(8) * (5.73d0) * ev2erg) ! H + C+++ -> C++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,6,Tk*1.001d0) * nN(ixHII) * nCarbon(icell) * xC(4,icell) * depletion_factors(8) * (5.73d0) * ev2erg)
#endif 
#if NCARBONIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,6,Tk) * nN(ixHII) * nCarbon(icell) * xC(5,icell) * depletion_factors(8) * (11.30d0) * ev2erg) ! H + C++++ -> C+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,6,Tk*1.001d0) * nN(ixHII) * nCarbon(icell) * xC(5,icell) * depletion_factors(8) * (11.30d0) * ev2erg)
#endif 
               end if

               ! Nitrogen
               if (nitrogen_ions) then
                  ct_heat_cool = ct_heat_cool + (HCTIon(1,7,Tk) * nI(ixHII) * nNitrogen(icell) * xN(1,icell) * depletion_factors(3) * (-0.94d0) * ev2erg) ! H+ + N -> N+ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(1,7,Tk*1.001d0) * nI(ixHII) * nNitrogen(icell) * xN(1,icell) * depletion_factors(3) * (-0.94d0) * ev2erg)
#if NNITROGENIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTRecom(2,7,Tk) * nN(ixHII) * nNitrogen(icell) * xN(2,icell) * depletion_factors(3) * (0.94d0) * ev2erg) ! H + N+ -> N + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(2,7,Tk*1.001d0) * nN(ixHII) * nNitrogen(icell) * xN(2,icell) * depletion_factors(3) * (0.94d0) * ev2erg)
#endif 
#if NNITROGENIONS>2
                  ct_heat_cool = ct_heat_cool + (HCTRecom(3,7,Tk) * nN(ixHII) * nNitrogen(icell) * xN(3,icell) * depletion_factors(3) * (4.56d0) * ev2erg) ! H + N++ -> N+ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(3,7,Tk*1.001d0) * nN(ixHII) * nNitrogen(icell) * xN(3,icell) * depletion_factors(3) * (4.56d0) * ev2erg)
#endif 
#if NNITROGENIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,7,Tk) * nN(ixHII) * nNitrogen(icell) * xN(4,icell) * depletion_factors(3) * (6.40d0) * ev2erg) ! H + N+++ -> N++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,7,Tk*1.001d0) * nN(ixHII) * nNitrogen(icell) * xN(4,icell) * depletion_factors(3) * (6.40d0) * ev2erg)
#endif 
#if NNITROGENIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,7,Tk) * nN(ixHII) * nNitrogen(icell) * xN(5,icell) * depletion_factors(3) * (11.00d0) * ev2erg) ! H + N++++ -> N+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,7,Tk*1.001d0) * nN(ixHII) * nNitrogen(icell) * xN(5,icell) * depletion_factors(3) * (11.00d0) * ev2erg)
#endif 
               end if

               ! Oxygen
               if (oxygen_ions) then
                  ct_heat_cool = ct_heat_cool + (HCTIon(1,8,Tk) * nI(ixHII) * nOxygen(icell) * xO(1,icell) * depletion_factors(2) * (-0.02d0) * ev2erg) ! H+ + O -> O+ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(1,8,Tk*1.001d0) * nI(ixHII) * nOxygen(icell) * xO(1,icell) * depletion_factors(2) * (-0.02d0) * ev2erg)
#if NOXYGENIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTRecom(2,8,Tk) * nN(ixHII) * nOxygen(icell) * xO(2,icell) * depletion_factors(2) * (0.02d0) * ev2erg) ! H + O+ -> O + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(2,8,Tk*1.001d0) * nN(ixHII) * nOxygen(icell) * xO(2,icell) * depletion_factors(2) * (0.02d0) * ev2erg)
#endif 
#if NOXYGENIONS>2
                  ct_heat_cool = ct_heat_cool + (HCTRecom(3,8,Tk) * nN(ixHII) * nOxygen(icell) * xO(3,icell) * depletion_factors(2) * (6.65d0) * ev2erg) ! H + O++ -> O+ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(3,8,Tk*1.001d0) * nN(ixHII) * nOxygen(icell) * xO(3,icell) * depletion_factors(2) * (6.65d0) * ev2erg)
#endif 
#if NOXYGENIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,8,Tk) * nN(ixHII) * nOxygen(icell) * xO(4,icell) * depletion_factors(2) * (5.00d0) * ev2erg) ! H + O+++ -> O++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,8,Tk*1.001d0) * nN(ixHII) * nOxygen(icell) * xO(4,icell) * depletion_factors(2) * (5.00d0) * ev2erg)
#endif 
#if NOXYGENIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,8,Tk) * nN(ixHII) * nOxygen(icell) * xO(5,icell) * depletion_factors(2) * (8.47d0) * ev2erg) ! H + O++++ -> O+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,8,Tk*1.001d0) * nN(ixHII) * nOxygen(icell) * xO(5,icell) * depletion_factors(2) * (8.47d0) * ev2erg)
#endif 
               end if

               ! Neon
               if (neon_ions) then
#if NNEONIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,10,Tk) * nN(ixHII) * nNeon(icell) * xNe(4,icell) * depletion_factors(5) * (5.82d0) * ev2erg) ! H + Ne+++ -> Ne++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,10,Tk*1.001d0) * nN(ixHII) * nNeon(icell) * xNe(4,icell) * depletion_factors(5) * (5.82d0) * ev2erg)
#endif 
#if NNEONIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,10,Tk) * nN(ixHII) * nNeon(icell) * xNe(5,icell) * depletion_factors(5) * (8.60d0) * ev2erg) ! H + Ne++++ -> Ne+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,10,Tk*1.001d0) * nN(ixHII) * nNeon(icell) * xNe(5,icell) * depletion_factors(5) * (8.60d0) * ev2erg)
#endif 
               end if

               ! Magnesium
               if (magnesium_ions) then
                  ct_heat_cool = ct_heat_cool + (HCTIon(1,12,Tk) * nI(ixHII) * nMagnesium(icell) * xMg(1,icell) * depletion_factors(4) * (1.52d0) * ev2erg) ! H+ + Mg -> Mg+ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(1,12,Tk*1.001d0) * nI(ixHII) * nMagnesium(icell) * xMg(1,icell) * depletion_factors(4) * (1.52d0) * ev2erg)
#if NMAGNESIUMIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTIon(2,12,Tk) * nI(ixHII) * nMagnesium(icell) * xMg(2,icell) * depletion_factors(4) * (-1.44d0) * ev2erg) ! H+ + Mg+ -> Mg++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(2,12,Tk*1.001d0) * nI(ixHII) * nMagnesium(icell) * xMg(2,icell) * depletion_factors(4) * (-1.44d0) * ev2erg)
#endif
#if NMAGNESIUMIONS>2
                  ct_heat_cool = ct_heat_cool + (HCTRecom(3,12,Tk) * nN(ixHII) * nMagnesium(icell) * xMg(3,icell) * depletion_factors(4) * (1.44d0) * ev2erg) ! H + Mg++ -> Mg+ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(3,12,Tk*1.001d0) * nN(ixHII) * nMagnesium(icell) * xMg(3,icell) * depletion_factors(4) * (1.44d0) * ev2erg)
#endif 
#if NMAGNESIUMIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,12,Tk) * nN(ixHII) * nMagnesium(icell) * xMg(4,icell) * depletion_factors(4) * (5.73d0) * ev2erg) ! H + Mg+++ -> Mg++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,12,Tk*1.001d0) * nN(ixHII) * nMagnesium(icell) * xMg(4,icell) * depletion_factors(4) * (5.73d0) * ev2erg)
#endif 
#if NMAGNESIUMIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,12,Tk) * nN(ixHII) * nMagnesium(icell) * xMg(5,icell) * depletion_factors(4) * (8.60d0) * ev2erg) ! H + Mg++++ -> Mg+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,12,Tk*1.001d0) * nN(ixHII) * nMagnesium(icell) * xMg(5,icell) * depletion_factors(4) * (8.60d0) * ev2erg)
#endif 
               end if

               ! Silicon
               if (silicon_ions) then
                  ct_heat_cool = ct_heat_cool + (HCTIon(1,14,Tk) * nI(ixHII) * nSilicon(icell) * xSi(1,icell) * depletion_factors(6) * (0.12d0) * ev2erg) ! H+ + Si -> Si+ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(1,14,Tk*1.001d0) * nI(ixHII) * nSilicon(icell) * xSi(1,icell) * depletion_factors(6) * (0.12d0) * ev2erg)
#if NSILICONIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTIon(2,14,Tk) * nI(ixHII) * nSilicon(icell) * xSi(2,icell) * depletion_factors(6) * (-2.72d0) * ev2erg) ! H+ + Si+ -> Si++ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(2,14,Tk*1.001d0) * nI(ixHII) * nSilicon(icell) * xSi(2,icell) * depletion_factors(6) * (-2.72d0) * ev2erg)
#endif
#if NSILICONIONS>2
                  ct_heat_cool = ct_heat_cool + (HCTRecom(3,14,Tk) * nN(ixHII) * nSilicon(icell) * xSi(3,icell) * depletion_factors(6) * (2.72d0) * ev2erg) ! H + Si++ -> Si+ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(3,14,Tk*1.001d0) * nN(ixHII) * nSilicon(icell) * xSi(3,icell) * depletion_factors(6) * (2.72d0) * ev2erg)
#endif 
#if NSILICONIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,14,Tk) * nN(ixHII) * nSilicon(icell) * xSi(4,icell) * depletion_factors(6) * (4.23d0) * ev2erg) ! H + Si+++ -> Si++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,14,Tk*1.001d0) * nN(ixHII) * nSilicon(icell) * xSi(4,icell) * depletion_factors(6) * (4.23d0) * ev2erg)
#endif 
#if NSILICONIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,14,Tk) * nN(ixHII) * nSilicon(icell) * xSi(5,icell) * depletion_factors(6) * (7.49d0) * ev2erg) ! H + Si++++ -> Si+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,14,Tk*1.001d0) * nN(ixHII) * nSilicon(icell) * xSi(5,icell) * depletion_factors(6) * (7.49d0) * ev2erg)
#endif 
               end if

               ! Sulfur
               if (sulfur_ions) then
#if NSULFURIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTRecom(2,16,Tk) * nN(ixHII) * nSulfur(icell) * xS(2,icell) * depletion_factors(9) * (-3.24d0) * ev2erg) ! H + S+ -> S + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(2,16,Tk*1.001d0) * nN(ixHII) * nSulfur(icell) * xS(2,icell) * depletion_factors(9) * (-3.24d0) * ev2erg)
#endif 
#if NSULFURIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,16,Tk) * nN(ixHII) * nSulfur(icell) * xS(4,icell) * depletion_factors(9) * (5.73d0) * ev2erg) ! H + S+++ -> S++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,16,Tk*1.001d0) * nN(ixHII) * nSulfur(icell) * xS(4,icell) * depletion_factors(9) * (5.73d0) * ev2erg)
#endif 
#if NSULFURIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,16,Tk) * nN(ixHII) * nSulfur(icell) * xS(5,icell) * depletion_factors(9) * (8.60d0) * ev2erg) ! H + S++++ -> S+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,16,Tk*1.001d0) * nN(ixHII) * nSulfur(icell) * xS(5,icell) * depletion_factors(9) * (8.60d0) * ev2erg)
#endif 
               end if

               ! Iron
               if (iron_ions) then
#if NIRONIONS>1
                  ct_heat_cool = ct_heat_cool + (HCTIon(2,26,Tk) * nI(ixHII) * nIron(icell) * xFe(2,icell) * depletion_factors(1) * (-2.56d0) * ev2erg) ! H+ + Fe+ -> Fe++ + H
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTIon(2,26,Tk*1.001d0) * nI(ixHII) * nIron(icell) * xFe(2,icell) * depletion_factors(1) * (-2.56d0) * ev2erg)
#endif
#if NIRONIONS>2
                  ct_heat_cool = ct_heat_cool + (HCTRecom(3,26,Tk) * nN(ixHII) * nIron(icell) * xFe(3,icell) * depletion_factors(1) * (2.56d0) * ev2erg) ! H + Fe++ -> Fe+ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(3,26,Tk*1.001d0) * nN(ixHII) * nIron(icell) * xFe(3,icell) * depletion_factors(1) * (2.56d0) * ev2erg)
#endif 
#if NIRONIONS>3
                  ct_heat_cool = ct_heat_cool + (HCTRecom(4,26,Tk) * nN(ixHII) * nIron(icell) * xFe(4,icell) * depletion_factors(1) * (6.30d0) * ev2erg) ! H + Fe+++ -> Fe++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(4,26,Tk*1.001d0) * nN(ixHII) * nIron(icell) * xFe(4,icell) * depletion_factors(1) * (6.30d0) * ev2erg)
#endif 
#if NIRONIONS>4
                  ct_heat_cool = ct_heat_cool + (HCTRecom(5,26,Tk) * nN(ixHII) * nIron(icell) * xFe(5,icell) * depletion_factors(1) * (10.00d0) * ev2erg) ! H + Fe++++ -> Fe+++ + H+
                  ct_heat_cool_prime = ct_heat_cool_prime + (HCTRecom(5,26,Tk*1.001d0) * nN(ixHII) * nIron(icell) * xFe(5,icell) * depletion_factors(1) * (10.00d0) * ev2erg)
#endif 
               end if

               ct_heat_cool_prime = (ct_heat_cool - ct_heat_cool_prime) / (TK - (1.001d0*TK))
            end if

            ! CO_cooling and CO_cooling_prime
            CO_cooling = 0.d0
            CO_cooling_prime = 0.d0
            if (isCO.and.nco_mol(icell).gt.min_cool_ion) then
               CO_cooling = comp_co_cooling(nH(icell), nN(ixHI), nN(ixHII), nco_mol(icell), TK)
               CO_cooling_prime = comp_co_cooling(nH(icell), nN(ixHI), nN(ixHII), nco_mol(icell), TK*1.001)
               CO_cooling_prime = (CO_cooling - CO_cooling_prime) / (TK - (TK*1.001))
            end if

            if (fire_cooling) then
               gH = 0.d0
               do igroup = 1, ngroups
                  !photoionization rate of neutral hydrogen
                  gH = gH + (rt_c*group_csn(igroup, 1)*scale_fp*dNp(igroup))
               end do

               call rt_cooling_fire(T2(icell), mu, nH(icell), Zsolar(icell), gH, 1.0 - dXion(1), ne, fire_tot)
               fire_tot = fire_tot*nH(icell)*nH(icell)
               !This is Harley's stupid hack to get the derivitive of the cooline rate
               call rt_cooling_fire(T2(icell)*1.001, mu, nH(icell), Zsolar(icell), gH, 1.0 - dXion(1), ne, fire_prime)
               fire_prime = fire_prime*nH(icell)*nH(icell)
               fire_prime = (fire_tot - fire_prime)/(T2(icell) - (T2(icell)*1.001))
            end if

            if (isH2Katz) H2heating_prime = (H2heating - H2heating_prime) / (TK - 1.001*TK)

            X_nHkb = 1.0 / (1.5 * nH(icell) * kB) 

            rate = X_nHkb*(Hrate + H2heating + ct_heat_cool - Crate - metal_tot - fire_tot - fine_structure - CO_cooling - cool_dust - dust_recomb_cool)

            ! Compute the cooling time and store it
            cooling_time(icell) = - (TK / rate) ! 3/2 * NkT / Lambda n^2, 
            
            dRate = -X_nHkb*(dCdT2 + metal_prime + fire_prime + fine_structure_prime + CO_cooling_prime + cool_dust_prime + dust_recomb_cool_prime  - HratePE_prime - H2heating_prime - ct_heat_cool_prime) ! dRate/dT2
            ! 1st order dt constr
            dUU = ABS(MAX(T2_min_fix, T2(icell) + rate*ddt(icell)) - T2(icell))
            ! New T2 value
            dT2 = MAX(T2_min_fix &
                      , T2(icell) + rate*ddt(icell)/(1.-dRate*ddt(icell)))
            dUU = MAX(dUU, ABS(dT2 - T2(icell)))/(T2(icell) + T_MIN) &
                  *one_over_T_FRAC
            if (dUU .gt. 1.) then                                     ! 10% rule
               code = 3; RETURN
            end if
            fracMax = MAX(fracMax, dUU)
            TK = dT2*mu
            temperature_save(icell) = TK
            !TODO(code): edit these to save everything that we want
            cooling_save(1 ,icell) = Crate + metal_tot + fire_tot + fine_structure + CO_cooling + cool_dust + dust_recomb_cool ! Total cooling rate
            cooling_save(2 ,icell) = Hrate + H2heating + ct_heat_cool ! Total heating rate
            cooling_save(3 ,icell) = Crate ! Primordial Cooling
            cooling_save(4 ,icell) = fine_structure ! All fine structure cooling
            cooling_save(5 ,icell) = fsc_CII ! CII cooling
            cooling_save(6 ,icell) = fsc_OI  ! OII cooling
            cooling_save(7 ,icell) = CO_cooling ! CO Cooling
            cooling_save(8 ,icell) = cool_dust ! Dust Cooling
            cooling_save(9 ,icell) = dust_recomb_cool ! Dust Cooling
            cooling_save(10,icell) = cr_heating_tot ! Cosmic ray heating HratePE_tot
            cooling_save(11,icell) = HratePE_tot ! Photoelectric heating
            cooling_save(12,icell) = H2heating ! H2 related heating
            cooling_save(13,icell) = ct_heat_cool ! Charge transfer heat cool
         end if

#if NGROUPS>0
         if (rt_isIR) then
            if (kAbs_loc(iIR) .gt. 0d0 .and. .not. rt_T_rad) then
               ! Evolve IR-Dust equilibrium temperature------------------------
               ! Delta (Cv T)= ( c_red/lambda E - c/lambda a T^4)
               !           / ( 1/Delta t + 4 c/lambda/C_v a T^3 + c_red/lambda)
               one_over_C_v = mh*mu*(gamma - 1d0)/(rho(icell)*kb)
               E_rad = group_egy_erg(iIR)*dNp(iIR)
               dE_T = (rt_c_cgs*E_rad - c_cgs*a_r*TK**4) &
                      /(1d0/(kAbs_loc(iIR)*Zsolar(icell)*rho(icell)*ddt(icell)) &
                        + 4d0*c_cgs*one_over_C_v*a_r*TK**3 + rt_c_cgs)
               dT2 = dT2 + 1d0/mu*one_over_C_v*dE_T
               dNp(iIR) = dNp(iIR) - dE_T*one_over_egy_IR_erg

               dT2 = max(T2_min_fix, dT2)
               dNp(iIR) = max(dNp(iIR), smallNp)
               ! 10% rule for photon density:
               dUU = ABS(dNp(iIR) - Np(iIR, icell))/(Np(iIR, icell) + Np_MIN) &
                     *one_over_Np_FRAC
               if (dUU .gt. 1.) then
                  code = 4; RETURN
               end if
               fracMax = MAX(fracMax, dUU)

               dUU = ABS(dT2 - T2(icell))/(T2(icell) + T_MIN)*one_over_T_FRAC
               if (dUU .gt. 1.) then                           ! 10% rule for T2
                  code = 5; RETURN
               end if
               fracMax = MAX(fracMax, dUU)
               TK = dT2*mu
               call reduce_flux(dFp(:, iIR), dNp(iIR)*rt_c_cgs)
            end if
         end if
#endif
         !(iii) UPDATE xH2*****************************************************
         if (isH2) then
            ! H2 creation rates
            G_0 = 0.d0
            if (rt_advect) then 
               G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
               do igroup = 1, nGroups  
                  if (group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                     G_0 = G_0 + (F_pe/G_conv)
                  end if
               end do
            end if 
            G_0 = G_0 + UV_background_G0

            effective_n = rho(icell) / mH ! Gas density assuming everything is hydrogen
            xe = ne / effective_n
            alpha(ixHI) = comp_Alpha_H2(TK, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, G_0, xHII, xHI) ! [cm3 s-1]
            cr = alpha(ixHI)*nH(icell)*xHI   ! H2 creation

            ! H2 photodissociation rates
            photoRate = 0.
            if (isH2Katz) then
               !if (loopcnt > 50000) then
               ! Nov 14 2023 --> use ole collisional destruction rates for H2
               beta(ixHI) = comp_Beta_H2coll(TK, nH(icell), xHI, xH2, xHeI, ne, nN(ixHII), nN(ixHI), nN(ixHeII)) ! H2 collisional dissociation rates
               !else
               !   beta(ixHI) = comp_Beta_H2coll_new(TK, nH(icell), ne, nN(ixHII), nN(ixHI), nN(ixHeII)) ! H2 collisional dissociation rates
               !end if
               
               if (.not.include_collisional_ionization) beta(ixHI) = 0.d0
               ! Note: we reduce the flux by H_2 photoionisation first then do HII later
               !       for the energy bin 15.2-24.59 ev
               !       because the cross-section for H_2 photoionisation is larger
               !       see Baczynski et al. (2015, MN, 454, 380) Fig. 1
               !       To be more accurate, this has to be determined by n(species)*sigma_pi
               if (rt.and.include_photoionisation) then
                  do igroup = 1, nGroups
                     if (igroup .eq. LWgroup) then
                        photoRate = photoRate + sigUV*dNp(igroup)
                     else
                        photoRate = photoRate + signc(igroup, ixHI)*dNp(igroup)
                     end if
                  end do
               end if
            else
               beta(ixHI) = comp_Beta_H2HI(TK)      ! H2 collisional dissociation rates
               if (.not.include_collisional_ionization) beta(ixHI) = 0.d0
               if (rt.and.include_photoionisation) then 
                  photoRate = SUM(signc(:, ixHI)*dNp)
               endif
            end if

            if (haardt_madau.and.include_photoionisation) then
               if (isH2Katz) then
                  photoRate = photoRate + UVrates(nIons + 1, 1)*ss_H2
               else
                  photoRate = photoRate + UVrates(ixHI, 1)*ss_factor
               end if
            end if

            if (include_photoionisation) then
               photoRate = photoRate + UV_background_h2*ss_factor
               photoRate = photoRate + UV_background_G0*5.68d-11 ! https://home.strw.leidenuniv.nl/~ewine/photo/
            end if
            de = beta(ixHI)*nH(icell) + photoRate

            if (include_cosmic_ray_ionization) then   
               de = de + H2_cosmic_ray_ionization_rate
            end if
            xH2 = (cr*ddt(icell) + xH2)/(1.+de*ddt(icell))
            xH2 = MIN(MAX(xH2, 1d-40), 0.5)
            if (xH2 > 1d-3) then
778            format(11(E12.5, 1x))
!           write(777,778) xH2, nH(icell), cr, de, photoRate, beta(ixHI),alpha(ixHI),xHI,xe,ddt(icell)
            end if
         end if

         !(iv) UPDATE xHI****************************************************
         ! First recompute interaction rates since T is updated
         if (rt_OTSA .or. .not. rt_advect) then           !  Recombination rates
            alpha(ixHII) = inp_coolrates_table(tbl_alphaB_HII, TK, dAlpha)
         else
            alpha(ixHII) = inp_coolrates_table(tbl_alphaA_HII, TK, dAlpha)
         end if
        
         beta(ixHII) = inp_coolrates_table(tbl_beta_HI, TK, dBeta) !  Coll-ion rate
         if (.not.include_collisional_ionization) beta(ixHII) = 0.d0

         cr = alpha(ixHII)*ne*xHII                              !    HI Creation
         if (isH2) cr = cr + 2.*de*xH2                          !    HI Creation
         
         if (include_charge_transfer) then
            if (isHe) then
               cr = cr + HEIIon(Tk)*xHII*nHe*xHeI ! He + H+ => He+ + H
            endif
            if (Zsolar(icell).gt.min_z_neq) then
               if (oxygen_ions) then
                  do im=1,n_oxygen_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,8,Tk)*xHII*dxO(im)*nOxygen(icell)*depletion_factors(2) ! Example:  O + H+ => O+ + H
                  enddo
               endif
               if (nitrogen_ions) then
                  do im=1,n_nitrogen_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,7,Tk)*xHII*dxN(im)*nNitrogen(icell)*depletion_factors(3) ! Example:  N + H+ => N+ + H
                  enddo
               endif
               if (carbon_ions) then
                  do im=1,n_carbon_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,6,Tk)*xHII*dxC(im)*nCarbon(icell)*depletion_factors(8) ! Example:  C + H+ => C+ + H
                  enddo
               endif
               if (magnesium_ions) then
                  do im=1,n_magnesium_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,12,Tk)*xHII*dxMg(im)*nMagnesium(icell)*depletion_factors(4) ! Example:  Mg + H+ => Mg+ + H
                  enddo
               endif
               if (silicon_ions) then
                  do im=1,n_silicon_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,14,Tk)*xHII*dxSi(im)*nSilicon(icell)*depletion_factors(6) ! Example:  Si + H+ => Si+ + H
                  enddo
               endif
               if (sulfur_ions) then
                  do im=1,n_sulfur_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,16,Tk)*xHII*dxS(im)*nSulfur(icell)*depletion_factors(9) ! Example:  S + H+ => S+ + H
                  enddo
               endif
               if (iron_ions) then
                  do im=1,n_iron_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,26,Tk)*xHII*dxFe(im)*nIron(icell)*depletion_factors(1) ! Example:  Fe + H+ => Fe+ + H
                  enddo
               endif
               if (neon_ions) then
                  do im=1,n_neon_ions-1 ! No ionization of final ion
                     cr = cr + HCTIon(im,10,Tk)*xHII*dxNe(im)*nNeon(icell)*depletion_factors(5) ! Example:  Ne + H+ => Ne+ + H
                  enddo
               endif
            endif
         endif
         
         photoRate = 0.
         if (rt.and.include_photoionisation) then 
            photoRate = SUM(signc(:, ixHII)*dNp)                !          [s-1]
            photoRate = photoRate + UV_background_hydrogen*ss_factor
         endif
         if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHII, 1)*ss_factor
         
         de = beta(ixHII)*ne + photoRate                         ! HI destruction
         
         if (isH2) de = de + 2.*alpha(ixHI)*nH(icell)            ! HI destruction
         
         ! Grain recombination
         if (include_grain_recombination.and.TK.le.1.d4) then
            G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
            G_0 = 0.d0
            do igroup = iIR + 1, nGroups
               if (group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                  F_pe = dNp(igroup)*rt_c_cgs         ! [#/cm2/s]
                  F_pe = F_pe*group_egy_erg(igroup) ! [erg/cm2/s]
                  G_0 = G_0 + F_pe
               end if
            end do
            G_0 = (G_0 / G_conv) + UV_background_G0
            if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(1,G_0,Tk,ne) * xHII * nH(icell) * dust_to_gas_mass_ratio_over_mw)
         end if
         
         if (include_cosmic_ray_ionization) then
            de = de + HI_cosmic_ray_ionization_rate
         end if
         
         if (include_charge_transfer) then
            if (isHe) then
               de = de + HEIIRecomb(Tk)*nHE*xHeII     ! H + He+  -> He + H+
               de = de + HEIIIRecomb(Tk)*nHE*xHeIII   ! H + He++ -> He+ + H+
            endif
            if (Zsolar(icell).gt.min_z_neq) then
               if (oxygen_ions) then
                  do im=2,n_oxygen_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,8,Tk)*dxO(im)*nOxygen(icell)*depletion_factors(2) ! Example:  O+ + H => O + H+
                  enddo
               endif
               if (nitrogen_ions) then
                  do im=2,n_nitrogen_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,7,Tk)*dxN(im)*nNitrogen(icell)*depletion_factors(3) ! Example:  N+ + H => N + H+
                  enddo
               endif
               if (carbon_ions) then
                  do im=2,n_carbon_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,6,Tk)*dxC(im)*nCarbon(icell)*depletion_factors(8) ! Example:  C+ + H => C + H+
                  enddo
               endif
               if (magnesium_ions) then
                  do im=2,n_magnesium_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,12,Tk)*dxMg(im)*nMagnesium(icell)*depletion_factors(4) ! Example:  Mg+ + H => Mg + H+
                  enddo
               endif
               if (silicon_ions) then
                  do im=2,n_silicon_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,14,Tk)*dxSi(im)*nSilicon(icell)*depletion_factors(6) ! Example:  Si+ + H => Si + H+
                  enddo
               endif
               if (sulfur_ions) then
                  do im=2,n_sulfur_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,16,Tk)*dxS(im)*nSulfur(icell)*depletion_factors(9) ! Example:  S+ + H => S + H+
                  enddo
               endif
               if (iron_ions) then
                  do im=2,n_iron_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,26,Tk)*dxFe(im)*nIron(icell)*depletion_factors(1) ! Example:  Fe+ + H => Fe + H+
                  enddo
               endif
               if (neon_ions) then
                  do im=2,n_neon_ions ! No recombination of groundstate
                     de = de + HCTRecom(im,10,Tk)*dxNe(im)*nNeon(icell)*depletion_factors(5) ! Example:  Ne+ + H => Ne + H+
                  enddo
               endif
            endif
         endif
         
         xHI = (cr*ddt(icell) + xHI)/(1.+de*ddt(icell))
         if (isH2) then
            xHI = MIN(MAX(xHI, x_MIN_HI), 1.d0) ! x_MIN_HI must be a tidy number 
            !--> otherwise, you will get interaction nN(ixHII)*PHrate(igroup, ixHII) even at T>>1e4K
            dXion(ixHI) = xHI
            dUU = ABS((xHI - xion(ixHI, icell))) &
              & /(xion(ixHI, icell) + x_MIN_HI)*one_over_x_FRAC
            if (dUU .gt. 1.) then
               code = 8; RETURN
            end if
            fracMax = MAX(fracMax, dUU)
         else
            xHI = MIN(MAX(xHI, x_MIN), 1.d0) !Must be x_MIN (if not H2)
         end if

         !(v) UPDATE xHII****************************************************
         cr = beta(ixHII)*ne                       !               Creation
         if (rt) cr = cr + photoRate               !                  [s-1]
         de = alpha(ixHII)*ne                      !            Destruction

         if (include_grain_recombination.and.TK.le.1.d4) then
            if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(1,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
         end if

         if (include_cosmic_ray_ionization) then
            cr = cr + ( HI_cosmic_ray_ionization_rate * xHI)
         end if

         if (include_charge_transfer) then
            if (isHe) then
               cr = cr + HEIIRecomb(Tk)*nHE*xHeII*xHI     ! H + He+  -> He + H+
               cr = cr + HEIIIRecomb(Tk)*nHE*xHeIII*xHI   ! H + He++ -> He+ + H+
               de = de + HEIIon(Tk)*nHE*xHeI ! He + H+ => He+ + H
            endif
            if (Zsolar(icell).gt.min_z_neq) then
               if (oxygen_ions) then
                  do im=1,n_oxygen_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,8,Tk)*dxO(im)*nOxygen(icell)*depletion_factors(2)*xHI ! Example:  O+ + H => O + H+
                     if (im.lt.n_oxygen_ions) de = de + HCTIon(im,8,Tk)*dxO(im)*nOxygen(icell)*depletion_factors(2) ! Example:  O + H+ => O+ + H
                  enddo
               endif
               if (nitrogen_ions) then
                  do im=1,n_nitrogen_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,7,Tk)*dxN(im)*nNitrogen(icell)*depletion_factors(3)*xHI  ! Example:  N+ + H => N + H+
                     if (im.lt.n_nitrogen_ions) de = de + HCTIon(im,7,Tk)*dxN(im)*nNitrogen(icell)*depletion_factors(3) ! Example:  N + H+ => N+ + H
                  enddo
               endif
               if (carbon_ions) then
                  do im=1,n_carbon_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,6,Tk)*dxC(im)*nCarbon(icell)*depletion_factors(8)*xHI  ! Example:  C+ + H => C + H+
                     if (im.lt.n_carbon_ions) de = de + HCTIon(im,6,Tk)*dxC(im)*nCarbon(icell)*depletion_factors(8) ! Example:  C + H+ => C+ + H
                  enddo
               endif
               if (magnesium_ions) then
                  do im=1,n_magnesium_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,12,Tk)*dxMg(im)*nMagnesium(icell)*depletion_factors(4)*xHI  ! Example:  Mg+ + H => Mg + H+
                     if (im.lt.n_magnesium_ions) de = de + HCTIon(im,12,Tk)*dxMg(im)*nMagnesium(icell)*depletion_factors(4) ! Example:  Mg + H+ => Mg+ + H
                  enddo
               endif
               if (silicon_ions) then
                  do im=1,n_silicon_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,14,Tk)*dxSi(im)*nSilicon(icell)*depletion_factors(6)*xHI  ! Example:  Si+ + H => Si + H+
                     if (im.lt.n_silicon_ions) de = de + HCTIon(im,14,Tk)*dxSi(im)*nSilicon(icell)*depletion_factors(6) ! Example:  Si + H+ => Si+ + H
                  enddo
               endif
               if (sulfur_ions) then
                  do im=1,n_sulfur_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,16,Tk)*dxS(im)*nSulfur(icell)*depletion_factors(9)*xHI  ! Example:  S+ + H => S + H+
                     if (im.lt.n_sulfur_ions) de = de + HCTIon(im,16,Tk)*dxS(im)*nSulfur(icell)*depletion_factors(9) ! Example:  S + H+ => S+ + H
                  enddo
               endif
               if (iron_ions) then
                  do im=1,n_iron_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,26,Tk)*dxFe(im)*nIron(icell)*depletion_factors(1)*xHI  ! Example:  Fe+ + H => Fe + H+
                     if (im.lt.n_iron_ions) de = de + HCTIon(im,26,Tk)*dxFe(im)*nIron(icell)*depletion_factors(1)  ! Example:  Fe + H+ => Fe+ + H
                  enddo
               endif
               if (neon_ions) then
                  do im=1,n_neon_ions 
                     if (im.gt.1) cr = cr + HCTRecom(im,10,Tk)*dxNe(im)*nNeon(icell)*depletion_factors(5)*xHI  ! Example:  Ne+ + H => Ne + H+
                     if (im.lt.n_neon_ions) de = de + HCTIon(im,10,Tk)*dxNe(im)*nNeon(icell)*depletion_factors(5) ! Example:  Ne + H+ => Ne+ + H
                  enddo
               endif
            endif
         endif

         ! Not Anninos, but more stable (this IS neccessary, as the one-cell  !
         ! tests oscillate wildly in the Anninos method):                     !
         S = cr*xHI - de*xHII
         dUU = ABS(MIN(MAX(xHII + ddt(icell)*S, x_MIN), 1.) - xHII)
         jac = xHI*(beta(ixHII)*nH(icell) - ne*TK*mu*X*dBeta) & !jac=dS/dxHII
               - cr - de - xHII*(alpha(ixHII)*nH(icell) - ne*TK*mu*X*dAlpha)
         xHII = xion(ixHII, icell) &
                + ddt(icell)*(cr*(1.-xion(ixHII, icell)) - de*xion(ixHII, icell)) &
                /(1.-ddt(icell)*jac)
         xHII = MIN(MAX(xHII, x_MIN), 1.d0)
         dXion(ixHII) = xHII
         dUU = MAX(dUU, ABS(xHII - xion(ixHII, icell)))/(xion(ixHII, icell) + x_MIN) &
               *one_over_x_FRAC
         if (dUU .gt. 1.) then
            code = 6; RETURN
         end if
         fracMax = MAX(fracMax, dUU)
         ! End a more stable and accurate integration---------------------------

         ! Atomic conservation of H
         if (isH2) then
            if (2.*xH2 .ge. dXion(ixHII)) then ! Either H2 or HI is most abundant
               if (2.*xH2 .le. dXion(ixHI)) dXion(ixHI) = 1.-2.*xH2 - dXion(ixHII)   ! -> HI
            else
               if (dXion(ixHI) .le. dXion(ixHII)) then
                  dXion(ixHII) = 1.-2.*xH2 - dXion(ixHI)                          !   HII
               else
                  dXion(ixHI) = 1.-2.*xH2 - dXion(ixHII)                         !    HI
               end if
            end if
         end if

         if (isHe) then
            ne = nH(icell)*xHII + nHE*(xHeII + 2.*xHeIII) ! Bc changed xhii
            ! Update the electron density for the presence of metals
            ne = ne + ne_metals

            mu = getMu(xHII, xHeII, xHeIII, dT2)
            if (.not. rt_isTconst) TK = dT2*mu !  Update TK because of changed  mu

            !(vi) UPDATE xHeI *************************************************
            if (rt_OTSA .or. .not. rt_advect) then
               alpha(ixHeII) = inp_coolrates_table(tbl_alphaB_HeII, TK)
               alpha(ixHeIII) = inp_coolrates_table(tbl_alphaB_HeIII, TK)
            else
               alpha(ixHeII) = inp_coolrates_table(tbl_alphaA_HeII, TK)
               alpha(ixHeIII) = inp_coolrates_table(tbl_alphaA_HeIII, TK)
            end if
            beta(ixHeII) = inp_coolrates_table(tbl_beta_HeI, TK)
            beta(ixHeIII) = inp_coolrates_table(tbl_beta_HeII, TK)
            if (.not.include_collisional_ionization) then 
               beta(ixHeII) = 0.d0
               beta(ixHeIII) = 0.d0
            endif
            ! Creation = recombination of HeII and electrons
            cr = alpha(ixHeII)*ne*xHeII

            if (include_grain_recombination.and.TK.le.1.d4) then
               if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(2,G_0,Tk,ne) * xHeII * nH(icell) * dust_to_gas_mass_ratio_over_mw)
            end if

            if (include_charge_transfer) then
               cr = cr + HEIIRecomb(Tk)*nH(icell)*xHI*xHeII  ! H + He+ -> He + H+
            endif

            ! Destruction = collisional ionization+photoionization of HeI
            de = beta(ixHeII)*ne

            photoRate = 0.d0
            if (rt.and.include_photoionisation) then 
               photoRate = SUM(signc(:, ixHeII)*dNp)
               photoRate = photoRate + UV_background_helium(1)*ss_factor
            endif
            if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
            de = de + photoRate

            if (include_charge_transfer) then
               de = de + HEIIon(Tk)*nH(icell)*xHII !He + H+ --> He+ + H
            endif

            if (include_cosmic_ray_ionization) then
               de = de + He_cosmic_ray_ionization_rate
            end if

            xHeI = (cr*ddt(icell) + xHeI)/(1.+de*ddt(icell))        !  The update
            xHeI = MIN(MAX(xHeI, 1d-40), 1.)

            !(vii) UPDATE xHeII *************************************************
            ! Creation = coll.- and photo-ionization of HeI + rec. of HeIII
            cr = de*xHeI + alpha(ixHeIII)*ne*xHeIII

            ! Destruction = rec. of HeII + coll.- and photo-ionization of HeII
            photoRate = 0.
            if (rt.and.include_photoionisation) then 
               photoRate = SUM(signc(:, ixHeIII)*dNp)
               photoRate = photoRate + UV_background_helium(2)*ss_factor
            endif
            if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeIII, 1)*ss_factor
            de = (alpha(ixHeII) + beta(ixHeIII))*ne + photoRate

            if (include_grain_recombination.and.TK.le.1.d4) then
               if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(2,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
            end if

            if (include_charge_transfer) then
               cr = cr + HEIIIRecomb(Tk)*nH(icell)*xHI*xHeIII ! H + He++ -> He+ + H+
               de = de + HEIIRecomb(Tk)*nH(icell)*xHI  ! H + He+ -> He + H+               
            endif

            xHeII = (cr*ddt(icell) + xHeII)/(1.+de*ddt(icell)) ! The update
            xHeII = MIN(MAX(xHeII, x_MIN), 1.)

            !(viii) UPDATE xHeIII **********************************************
            ! Creation = coll.- and photo-ionization of HeII
            cr = (beta(ixHeIII)*ne + photoRate)*xHeII          !  xHeII is new
            ! Destruction = rec. of HeIII and e
            de = alpha(ixHeIII)*ne

            if (include_charge_transfer) then
               de = de + HEIIIRecomb(Tk)*nH(icell)*xHI  ! H + He++ -> He+ + H+
            endif

            xHeIII = (cr*ddt(icell) + xHeIII)/(1.+de*ddt(icell)) ! The update
            xHeIII = MIN(MAX(xHeIII, x_MIN), 1.)

            !(ix) ATOMIC CONSERVATION OF He *********************************
            if(xHeI .ge. xHeIII) then   ! Either HeI or HeII is most abundant 
               if(xHeI .le. xHeII) xHeII = 1. - xHeI - xHeIII   !HeII most ab
            else                        ! Either HeII or HeIII is most abundant 
               if(xHeII .le. xHeIII) then 
                  xHeIII = 1. - xHeI - xHeII                         ! HeIII
               else 
                  xHeII  = 1. - xHeI - xHeIII                        !  HeII
               endif
            endif
            dXion(ixHeII) = xHeII
            dXion(ixHeIII) = xHeIII

         end if

         ne = nH(icell)*xHII + nHe*(xHeII + 2.*xHeIII) ! Bc changed xhei, xheii, xheii

         ! START UPDATE METAL IONS
         if (Zsolar(icell).gt.min_z_neq) then
         ! update the electron density again if necessary (fixed across metal calc as it's a small effect)
         if (oxygen_ions.or.nitrogen_ions.or.carbon_ions.or.magnesium_ions.or.silicon_ions.or.sulfur_ions.or.iron_ions.or.neon_ions) then
            ne = ne + ne_metals
         endif

         mu = getMu(xHII, xHeII, xHeIII, dT2)
         if (.not. rt_isTconst) TK = dT2*mu !  Update TK because of changed mu


         ! Start with CO
         if (isCO) then
            ! ddt(icell): correct time step
            ! nH(icell): number density of hydrogen nuclei
            ! nN(ixHII): nHI
            ! xC(2,icell) # Relative abundance to carbon of carbonII, 2=xCII
            ! nCII = xC(2,icell)*nCarbon(icell)*depletion_factors(8), 8 is a hardcoded index
            ! xO(1,icell) Relative abundance to oxygen of oxygenI
            ! xOI = xO(1,icell)*nOxygen(icell)*depletion_factors(2)/nH(icell)
            ! nOxygen(icell)
            ! xH2 # Relative abundance to hydrogen nuclei of molecular hydrogen nuclei
            ! TK # Temperature
            ! xe # Relative abundance of e-
            ! xHI # Relative abundance of H
            ! xHeI # Relative abundance of HeI to He
            ! nCarbon(icell): total carbon number density
            ! nN(ixHI): nH2

            ! Compute G_0 first
            G_conv = 1.6d-3 ! [erg/s/cm2] (Habing 1968)
            G_0 = 0.d0
            if (rt.and.rt_advect) then 
               do igroup = iIR + 1, nGroups
                  if (group_egy(igroup) > 5.4 .and. group_egy(igroup) < 13.6 .and. TK .le. T_sputter) then
                     F_pe = dNp(igroup)*rt_c_cgs         ! [#/cm2/s]
                     F_pe = F_pe*group_egy_erg(igroup) ! [erg/cm2/s]
                     G_0 = G_0 + F_pe
                  end if
               end do
            end if
            CO_SS = comp_SCO(nco_mol(icell), nH(icell)*xH2, dx_SS(icell)) * comp_Sd(nH(icell)*xHI, nH(icell)*xH2, dx_SS(icell), dust_to_gas_mass_ratio_over_mw)
            G_0 = ((G_0 / G_conv)*CO_SS) + UV_background_G0

            cr = 0.d0
            de = 0.d0
            cr = comp_cr_co(G_0, H2_cosmic_ray_ionization_rate, xC(2,icell)*nCarbon(icell)*depletion_factors(8), nN(ixHI), xO(1,icell)*nOxygen(icell)*depletion_factors(2)/nH(icell), nH(icell)) ! This should have the same dimension as nco_mol(icell)/t
            de = comp_de_co(G_0, H2_cosmic_ray_ionization_rate) ! This should have the same dimension as 1/t
            ! To test as if CI also forms CO, add this term to CII: +xC(1,icell)*nCarbon(icell)*depletion_factors(8)

            ! The update
            call update_co_chem(cr, de, ddt(icell), n_carbon_ions, n_oxygen_ions, nco_mol(icell), nCarbon(icell), dxC, nOxygen(icell), dxO)
            
         end if

         if (oxygen_ions) then
               do im=1,n_oxygen_ions
                  ! Creation
                  cr = 0.d0
                  if (im.lt.n_oxygen_ions) then 
                     cr = cr + comp_Alpha_oxygen(TK, 8-im)*ne*dxO(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,8,Tk)*nH(icell)*xHI*dxO(im+1)! Charge transfer recombination of more excited ionization state (H)
                        !if (isHe) cr = cr + O_HE_Recomb(im+1,Tk)*nHE*xHeI**dxO(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_oxygen(TK, im-1)*ne*dxO(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,8,Tk)*nH(icell)*xHII*dxO(im-1)! Charge transfer ionization of less excited ionization state (H)
                        !if (isHe) cr = cr + O_HE_Ion(im-1,Tk)*nHE*xHeII**dxO(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_oxygen(:, im-1)*dNp)*dxO(im-1)  ! photoionization of the less excited state
                        cr = cr + UV_background_oxygen(im-1)*dxO(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_oxygen(im-1)*dxO(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_oxygen_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_oxygen(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,8,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state (H)
                        !if (isHe) de = de + O_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_oxygen(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_oxygen(TK, 9-im)*ne ! Recombination (the 9 is correct)
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,8,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + O_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to less excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_OI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_OI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxO(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_oxygen(:, im)*dNp) ! photoionization
                     if (im.lt.n_oxygen_ions) then
                        photoRate = photoRate + UV_background_oxygen(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                  de = de + photoRate

                  dxO(im) = (cr*ddt(icell) + dxO(im))/(1.+de*ddt(icell))        !  The update
                  dxO(im) = MIN(MAX(dxO(im), 1d-40), 1.)
               enddo
         end if

         if (nitrogen_ions) then
               do im=1,n_nitrogen_ions
                  ! Creation
                  cr = 0.d0
                  if (im.lt.n_nitrogen_ions) then 
                     cr = cr + comp_Alpha_nitrogen(TK, 7-im)*ne*dxN(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,7,Tk)*nH(icell)*xHI*dxN(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + N_HE_Recomb(im+1,Tk)*nHE*xHeI**dxN(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_nitrogen(TK, im-1)*ne*dxN(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,7,Tk)*nH(icell)*xHII*dxN(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + N_HE_Ion(im-1,Tk)*nHE*xHeII**dxN(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_nitrogen(:, im-1)*dNp)*dxN(im-1) ! photoionization of the less excited state
                        cr = cr + UV_background_nitrogen(im-1)*dxN(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_nitrogen(im-1)*dxN(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_nitrogen_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_nitrogen(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,7,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + N_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_nitrogen(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_nitrogen(TK, 8-im)*ne ! Recombination (the 8 is correct)
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,7,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + N_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer recombination to less excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_NI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_NI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxN(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_nitrogen(:, im)*dNp) ! photoionization
                     if (im.lt.n_nitrogen_ions) then
                        photoRate = photoRate + UV_background_nitrogen(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                  de = de + photoRate

                  dxN(im) = (cr*ddt(icell) + dxN(im))/(1.+de*ddt(icell))        !  The update
                  dxN(im) = MIN(MAX(dxN(im), 1d-40), 1.)
               enddo
         endif

         if (carbon_ions) then
               do im=1,n_carbon_ions
                  ! Creation
                  cr = 0.d0
                  if (im.lt.n_carbon_ions) then 
                     cr = cr + comp_Alpha_carbon(TK, 6-im)*ne*dxC(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,6,Tk)*nH(icell)*xHI*dxC(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + C_HE_Recomb(im+1,Tk)*nHE*xHeI**dxC(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_carbon(TK, im-1)*ne*dxC(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,6,Tk)*nH(icell)*xHII*dxC(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + C_HE_Ion(im-1,Tk)*nHE*xHeII**dxC(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_carbon(:, im-1)*dNp)*dxC(im-1)    ! photoionization of the less excited state
                        cr = cr + UV_background_carbon(im-1)*dxC(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_carbon(im-1)*dxC(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_carbon_ions) then
                     if (include_collisional_ionization) de = de + comp_beta_carbon(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,6,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + C_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_carbon(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_carbon(TK, 7-im)*ne ! Recombination (the 7 is correct)
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,6,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + C_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to more excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_CI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_CI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxC(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_carbon(:, im)*dNp) ! photoionization
                     if (im.lt.n_carbon_ions) then
                        photoRate = photoRate + UV_background_carbon(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor

                  ! Destruction due to G0
                  if (include_photoionisation) then
                     if (im.eq.1) photoRate = photoRate + UV_background_G0*3.39d-10 ! https://home.strw.leidenuniv.nl/~ewine/photo/
                     if (im.eq.2) cr = cr + (UV_background_G0*3.39d-10*dxC(im-1)) ! https://home.strw.leidenuniv.nl/~ewine/photo/
                  end if

                  de = de + photoRate

                  if (include_grain_recombination.and.TK.le.1.d4) then
                     if (im.eq.1) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(6,G_0,Tk,ne) * dxC(2) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     else if (im.eq.2) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(6,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     end if
                  end if

                  dxC(im) = (cr*ddt(icell) + dxC(im))/(1.+de*ddt(icell))        !  The update
                  dxC(im) = MIN(MAX(dxC(im), 1d-40), 1.)
               enddo
         endif

         if (magnesium_ions) then
               do im=1,n_magnesium_ions
                  ! Creation = rec of more excited state, coll ion of less excited state, photo ion of less excited state 
                  cr = 0.d0
                  if (im.lt.n_magnesium_ions) then 
                     cr = cr + comp_Alpha_magnesium(TK, 12-im)*ne*dxMg(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,12,Tk)*nH(icell)*xHI*dxMg(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + MG_HE_Recomb(im+1,Tk)*nHE*xHeI**dxMg(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_magnesium(TK, im-1)*ne*dxMg(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,12,Tk)*nH(icell)*xHII*dxMg(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + MG_HE_Ion(im-1,Tk)*nHE*xHeII**dxMg(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_magnesium(:, im-1)*dNp)*dxMg(im-1)   ! photoionization of the less excited state
                        cr = cr + UV_background_magnesium(im-1)*dxMg(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_magnesium(im-1)*dxMg(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_magnesium_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_magnesium(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,12,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + MG_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_magnesium(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_magnesium(TK, 13-im)*ne ! Recombination 
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,12,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + MG_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to more excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_MgI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_MgI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxMg(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_magnesium(:, im)*dNp) ! photoionization
                     if (im.lt.n_magnesium_ions) then
                        photoRate = photoRate + UV_background_magnesium(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor

                  ! Destruction due to G0
                  if (include_photoionisation) then
                     if (im.eq.1) photoRate = photoRate + UV_background_G0*6.59d-11 ! https://home.strw.leidenuniv.nl/~ewine/photo/
                     if (im.eq.2) cr = cr + (UV_background_G0*6.59d-11*dxMg(im-1)) ! https://home.strw.leidenuniv.nl/~ewine/photo/
                  end if

                  de = de + photoRate

                  if (include_grain_recombination.and.TK.le.1.d4) then
                     if (im.eq.1) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(12,G_0,Tk,ne) * dxMg(2) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     else if (im.eq.2) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(12,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     end if
                  end if

                  dxMg(im) = (cr*ddt(icell) + dxMg(im))/(1.+de*ddt(icell))        !  The update
                  dxMg(im) = MIN(MAX(dxMg(im), 1d-40), 1.)
               enddo
         endif

         if (silicon_ions) then
               do im=1,n_silicon_ions
                  ! Creation = rec of more excited state, coll ion of less excited state, photo ion of less excited state 
                  cr = 0.d0
                  if (im.lt.n_silicon_ions) then 
                     cr = cr + comp_Alpha_silicon(TK, 14-im)*ne*dxSi(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,14,Tk)*nH(icell)*xHI*dxSi(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + SI_HE_Recomb(im+1,Tk)*nHE*xHeI**dxSi(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_silicon(TK, im-1)*ne*dxSi(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,14,Tk)*nH(icell)*xHII*dxSi(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + SI_HE_Ion(im-1,Tk)*nHE*xHeII**dxSi(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_silicon(:, im-1)*dNp)*dxSi(im-1)   ! photoionization of the less excited state
                        cr = cr + UV_background_silicon(im-1)*dxSi(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_silicon(im-1)*dxSi(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_silicon_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_silicon(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,14,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + Si_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_silicon(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_silicon(TK, 15-im)*ne ! Recombination 
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,14,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + SI_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to more excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_SiI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_SiI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxSi(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_silicon(:, im)*dNp) ! photoionization
                     if (im.lt.n_silicon_ions) then
                        photoRate = photoRate + UV_background_silicon(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor

                  ! Destruction due to G0
                  if (include_photoionisation) then
                     if (im.eq.1) photoRate = photoRate + UV_background_G0*4.47d-9 ! https://home.strw.leidenuniv.nl/~ewine/photo/
                     if (im.eq.2) cr = cr + (UV_background_G0*4.47d-9*dxSi(im-1)) ! https://home.strw.leidenuniv.nl/~ewine/photo/
                  end if

                  de = de + photoRate

                  if (include_grain_recombination.and.TK.le.1.d4) then
                     if (im.eq.1) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(14,G_0,Tk,ne) * dxSi(2) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     else if (im.eq.2) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(14,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     end if
                  end if

                  dxSi(im) = (cr*ddt(icell) + dxSi(im))/(1.+de*ddt(icell))        !  The update
                  dxSi(im) = MIN(MAX(dxSi(im), 1d-40), 1.)
               enddo
         endif

         if (sulfur_ions) then
               do im=1,n_sulfur_ions
                  ! Creation = rec of more excited state, coll ion of less excited state, photo ion of less excited state 
                  cr = 0.d0
                  if (im.lt.n_sulfur_ions) then 
                     cr = cr + comp_Alpha_sulfur(TK, 16-im)*ne*dxS(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,16,Tk)*nH(icell)*xHI*dxS(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + S_HE_Recomb(im+1,Tk)*nHE*xHeI**dxS(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_sulfur(TK, im-1)*ne*dxS(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,16,Tk)*nH(icell)*xHII*dxS(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + S_HE_Ion(im-1,Tk)*nHE*xHeII**dxS(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_sulfur(:, im-1)*dNp)*dxS(im-1)   ! photoionization of the less excited state
                        cr = cr + UV_background_sulfur(im-1)*dxS(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_sulfur(im-1)*dxS(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_sulfur_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_sulfur(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,16,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + S_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_sulfur(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_sulfur(TK, 17-im)*ne ! Recombination 
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,16,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + S_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to more excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_SI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_SI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxS(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_sulfur(:, im)*dNp) ! photoionization
                     if (im.lt.n_sulfur_ions) then
                        photoRate = photoRate + UV_background_sulfur(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor

                  ! Destruction due to G0
                  if (include_photoionisation) then
                     if (im.eq.1) photoRate = photoRate + UV_background_G0*1.13d-9 ! https://home.strw.leidenuniv.nl/~ewine/photo/
                     if (im.eq.2) cr = cr + (UV_background_G0*1.13d-9*dxS(im-1)) ! https://home.strw.leidenuniv.nl/~ewine/photo/
                  end if

                  de = de + photoRate

                  if (include_grain_recombination.and.TK.le.1.d4) then
                     if (im.eq.1) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(16,G_0,Tk,ne) * dxS(2) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     else if (im.eq.2) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(16,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     end if
                  end if

                  dxS(im) = (cr*ddt(icell) + dxS(im))/(1.+de*ddt(icell))        !  The update
                  dxS(im) = MIN(MAX(dxS(im), 1d-40), 1.)
               enddo
         endif

         if (iron_ions) then
               do im=1,n_iron_ions
                  ! Creation = rec of more excited state, coll ion of less excited state, photo ion of less excited state 
                  cr = 0.d0
                  if (im.lt.n_iron_ions) then 
                     cr = cr + comp_Alpha_iron(TK, 26-im)*ne*dxFe(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,26,Tk)*nH(icell)*xHI*dxFe(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + FE_HE_Recomb(im+1,Tk)*nHE*xHeI**dxFe(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_iron(TK, im-1)*ne*dxFe(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,26,Tk)*nH(icell)*xHII*dxFe(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + FE_HE_Ion(im-1,Tk)*nHE*xHeII**dxFe(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_iron(:, im-1)*dNp)*dxFe(im-1)   ! photoionization of the less excited state
                        cr = cr + UV_background_iron(im-1)*dxFe(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_iron(im-1)*dxFe(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_iron_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_iron(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,26,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + FE_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_iron(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_iron(TK, 27-im)*ne ! Recombination 
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,26,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + FE_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to more excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_FeI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_FeI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxFe(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_iron(:, im)*dNp) ! photoionization
                     if (im.lt.n_iron_ions) then
                        photoRate = photoRate + UV_background_iron(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor

                  ! Destruction due to G0
                  if (include_photoionisation) then
                     if (im.eq.1) photoRate = photoRate + UV_background_G0*4.71d-10 ! https://home.strw.leidenuniv.nl/~ewine/photo/
                     if (im.eq.2) cr = cr + (UV_background_G0*4.71d-10*dxFe(im-1)) ! https://home.strw.leidenuniv.nl/~ewine/photo/
                  end if

                  de = de + photoRate

                  if (include_grain_recombination.and.TK.le.1.d4) then
                     if (im.eq.1) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) cr = cr + ( DUST_RECOMBINATION(26,G_0,Tk,ne) * dxFe(2) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     else if (im.eq.2) then
                        if(dust_to_gas_mass_ratio_over_mw.gt.1.d-3) de = de + ( DUST_RECOMBINATION(26,G_0,Tk,ne) * nH(icell) * dust_to_gas_mass_ratio_over_mw)
                     end if
                  end if

                  dxFe(im) = (cr*ddt(icell) + dxFe(im))/(1.+de*ddt(icell))        !  The update
                  dxFe(im) = MIN(MAX(dxFe(im), 1d-40), 1.)
               enddo
         endif

         if (neon_ions) then
               do im=1,n_neon_ions
                  ! Creation = rec of more excited state, coll ion of less excited state, photo ion of less excited state 
                  cr = 0.d0
                  if (im.lt.n_neon_ions) then 
                     cr = cr + comp_Alpha_neon(TK, 10-im)*ne*dxNe(im+1) ! Recombinations of the more excited ionization state
                     if (include_charge_transfer) then 
                        cr = cr + HCTRecom(im+1,10,Tk)*nH(icell)*xHI*dxNe(im+1)! Charge transfer recombination of more excited ionization state
                        !if (isHe) cr = cr + NE_HE_Recomb(im+1,Tk)*nHE*xHeI**dxNe(im+1)! Charge transfer recombination of more excited ionization state (He)
                     endif
                  endif
                  if (im.gt.1) then 
                     if (include_collisional_ionization) cr = cr + comp_beta_neon(TK, im-1)*ne*dxNe(im-1) ! Collisional ionization of the less excited state
                     if (include_charge_transfer) then 
                        cr = cr + HCTIon(im-1,10,Tk)*nH(icell)*xHII*dxNe(im-1)! Charge transfer ionization of less excited ionization state
                        !if (isHe) cr = cr + NE_HE_Ion(im-1,Tk)*nHE*xHeII**dxNe(im-1)! Charge transfer ionization of less excited ionization state (He)
                     endif
                     if (rt.and.include_photoionisation) then 
                        cr = cr + SUM(signc_neon(:, im-1)*dNp)*dxNe(im-1)   ! photoionization of the less excited state
                        cr = cr + UV_background_neon(im-1)*dxNe(im-1)*ss_factor ! photoionization of the less excited state (UV BG)
                     endif
                     if (include_cosmic_ray_ionization) cr = cr + cr_ionization_neon(im-1)*dxNe(im-1)*xi_h_cr
                  endif

                  ! Destruction = collisional ionization + photoionization + recombination
                  de = 0.d0
                  if (im.lt.n_neon_ions) then 
                     if (include_collisional_ionization) de = de + comp_beta_neon(TK, im)*ne ! Collisional ionization 
                     if (include_charge_transfer) then 
                        de = de + HCTIon(im,10,Tk)*nH(icell)*xHII ! Charge transfer ionization to more excited state
                        !if (isHe) de = de + NE_HE_Ion(im,Tk)*nHe*xHeII ! Charge transfer ionization to more excited state (He)
                     endif
                     if (include_cosmic_ray_ionization) de = de + cr_ionization_neon(im)*xi_h_cr
                  endif
                  if (im.gt.1) then 
                     de = de + comp_Alpha_neon(TK, 11-im)*ne ! Recombination 
                     if (include_charge_transfer) then 
                        de = de + HCTRecom(im,10,Tk)*nH(icell)*xHI ! Charge transfer recombination to less excited state
                        !if (isHe) de = de + NE_HE_Recomb(im,Tk)*nHe*xHeI ! Charge transfer ionization to more excited state (He)
                     endif
                  endif

                  ! Destruction and creation from cosmic ray-generature FUV photons
                  if (include_cosmic_ray_ionization) then
                     if (im.eq.1) de = de + (cr_ionization_ind_UV_NeI * (H2_cosmic_ray_ionization_rate / 1.d-16))
                     if (im.eq.2) cr = cr + (cr_ionization_ind_UV_NeI * (H2_cosmic_ray_ionization_rate / 1.d-16) * dxNe(im-1))
                  end if

                  photoRate = 0.d0
                  if (rt.and.include_photoionisation) then 
                     photoRate = SUM(signc_neon(:, im)*dNp) ! photoionization
                     if (im.lt.n_neon_ions) then
                        photoRate = photoRate + UV_background_neon(im)*ss_factor ! photoionization from uv background
                     endif
                  endif
                  !if (haardt_madau.and.include_photoionisation) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                  de = de + photoRate

                  dxNe(im) = (cr*ddt(icell) + dxNe(im))/(1.+de*ddt(icell))        !  The update
                  dxNe(im) = MIN(MAX(dxNe(im), 1d-40), 1.)
               enddo
         endif

         endif
         ! END UPDATE METAL IONS


         ! conserve metals
         if (Zsolar(icell).gt.min_z_neq) then
         if (oxygen_ions)    dxO  = dxO  / SUM(dxO)
         if (nitrogen_ions)  dxN  = dxN  / SUM(dxN)
         if (carbon_ions)    dxC  = dxC  / SUM(dxC)
         if (magnesium_ions) dxMg = dxMg / SUM(dxMg)
         if (silicon_ions)   dxSi = dxSi / SUM(dxSi)
         if (sulfur_ions)    dxS  = dxS  / SUM(dxS)
         if (iron_ions)      dxFe = dxFe / SUM(dxFe)
         if (neon_ions)      dxNe = dxNe / SUM(dxNe)
         endif

         ne = nH(icell)*xHII + nHe*(xHeII + 2.*xHeIII)
         ! Update the electron density for the presence of metals
         ne_metals = 0.d0
         if (Zsolar(icell).gt.min_z_neq) then
            if (oxygen_ions) then
               do im=2,n_oxygen_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nOxygen(icell) * xO(im,icell) * (im-1))
               enddo
            endif
            if (nitrogen_ions) then
               do im=2,n_nitrogen_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nNitrogen(icell) * xN(im,icell) * (im-1))
               enddo
            endif
            if (carbon_ions) then
               do im=2,n_carbon_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nCarbon(icell) * xC(im,icell) * (im-1))
               enddo
            endif
            if (magnesium_ions) then
               do im=2,n_magnesium_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nMagnesium(icell) * xMg(im,icell) * (im-1))
               enddo
            endif
            if (silicon_ions) then
               do im=2,n_silicon_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nSilicon(icell) * xSi(im,icell) * (im-1))
               enddo
            endif
            if (sulfur_ions) then
               do im=2,n_sulfur_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nSulfur(icell) * xS(im,icell) * (im-1))
               enddo
            endif
            if (iron_ions) then
               do im=2,n_iron_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nIron(icell) * xFe(im,icell) * (im-1))
               enddo
            endif
            if (neon_ions) then
               do im=2,n_neon_ions ! start index at 2 (ground state = no free electrons)
                  ne_metals = ne_metals + (nNeon(icell) * xNe(im,icell) * (im-1))
               enddo
            endif
         endif
         ne = ne + ne_metals

         dUU = ABS((ne - neInit))/(neInit + x_MIN)*one_over_x_FRAC
         if (dUU .gt. 1.) then
            code = 7; RETURN
         end if
         fracMax = MAX(fracMax, dUU)

         if (rt_isTconst) dT2 = rt_Tconst/mu

         dT2 = dT2 - T2(icell); dXion(:) = dXion(:) - xion(:, icell)
         dNp(:) = dNp(:) - Np(:, icell); dFp(:, :) = dFp(:, :) - Fp(:, :, icell)
         dp_gas(:) = dp_gas(:) - p_gas(:, icell)

         if (Zsolar(icell).gt.min_z_neq) then
            if (oxygen_ions)    dxO(:)  = dxO(:)  - xO(:, icell)
            if (nitrogen_ions)  dxN(:)  = dxN(:)  - xN(:, icell)
            if (carbon_ions)    dxC(:)  = dxC(:)  - xC(:, icell)
            if (magnesium_ions) dxMg(:) = dxMg(:) - xMg(:, icell)
            if (silicon_ions)   dxSi(:) = dxSi(:) - xSi(:, icell)
            if (sulfur_ions)    dxS(:)  = dxS(:)  - xS(:, icell)
            if (iron_ions)      dxFe(:) = dxFe(:) - xFe(:, icell)
            if (neon_ions)      dxNe(:) = dxNe(:) - xNe(:, icell)
         endif

         ! Now the dUs are really changes, not new values
         !(ix) Check if we are safe to use a bigger timestep in next iteration:
         if (fracMax .lt. 0.5) then
            dt_rec = ddt(icell)*2.
         else
            dt_rec = ddt(icell)
         end if
         dt_ok = .true.
         code = 0

      
      END SUBROUTINE cool_step

   END SUBROUTINE rt_solve_cooling

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE display_coolinfo(stopRun, loopcnt, i, dtDone, dt, ddt, nH &
                               , T2, xion, Np, Fp, p_gas &
                               , dT2, dXion, dNp, dFp, dp_gas, code, cool_s &
                               , nOxygen, nNitrogen, nCarbon, nMagnesium &
                               , nSilicon, nSulfur, nIron, nNeon, ncomol)
! Print cooling information to standard output, and maybe stop execution.
!------------------------------------------------------------------------
      use amr_commons
      use rt_parameters
      real(dp), dimension(nIons):: xion, dXion
      real(dp), dimension(nGroups):: Np, dNp
      real(dp), dimension(nDim, nGroups):: Fp, dFp
      real(dp), dimension(nDim):: p_gas, dp_gas
      real(dp)::T2, dT2, dtDone, dt, ddt, nH
      real(dp)::nOxygen, nNitrogen, nCarbon, nMagnesium, nSilicon, nSulfur, nIron, nNeon, ncomol
      real(dp), dimension(ncoolheatratesave)::cool_s
      logical::stopRun
      integer::loopcnt, i, code
!------------------------------------------------------------------------
      if (stopRun) write (*, 111) loopcnt
      if (.true.) then
         write (*, 900) loopcnt, myid, code, i, dtDone, dt, ddt, rt_c_cgs, nH
         write (*, 901) T2, xion, Np, Fp, p_gas
         write (*, 902) dT2, dXion, dNp, dFp, dp_gas
         write (*, 903) dT2/ddt, dXion/ddt, dNp/ddt, dFp/ddt, dp_gas/ddt
         write (*, 904) abs(dT2)/(T2 + T_MIN), abs(dxion)/(xion + x_MIN), &
            abs(dNp)/(Np + Np_MIN), abs(dFp)/(Fp + Fp_MIN)
         write (*, 905) cool_s(:)
         write (*, 906) nOxygen, nNitrogen, nCarbon, nMagnesium, nSilicon, nSulfur, nIron, nNeon, ncomol
      end if
      if (stopRun) then
         print *, 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
         STOP
      end if

111   format(' Stopping because of large number of timestesps in', &
             ' rt_solve_cooling (', I6, ')')
900   format(I3, '  myid=', I2, ' code=', I2, ' i=', I5, ' t=', 1pe12.3, xs &
             '/', 1pe12.3, ' ddt=', 1pe12.3, ' c=', 1pe12.3, &
             ' nH=', 1pe12.3)
901   format('  U      =', 20(1pe12.3))
902   format('  dU     =', 20(1pe12.3))
903   format('  dU/dt  =', 20(1pe12.3))
904   format('  dU/U % =', 20(1pe12.3))
905   format('  coolv  =', 12(1pe12.3))
906   format('  metal  =',  9(1pe12.3))
   END SUBROUTINE display_coolinfo

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE cmp_chem_eq(TK, nH, t_rad_spec, nSpec, nTot, mu, Zsol)

! Compute chemical equilibrium abundances of e, H2, HI, HII, HeI, HeII, HeIII
! r_rad_spec => photoionization rates [s-1] for H2, HI, HeI, HeII
!------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::TK, nH, Zsol
      real(dp), intent(out)::nTot, mu
      real(dp), dimension(nIons), intent(in)::t_rad_spec
      real(dp), dimension(1:7), intent(out)::nSpec!------------------------
      real(dp)::nHe, xe
      real(dp)::n_H2, n_HI, n_HII, n_HEI, n_HEII, n_HEIII, n_E, n_E_min
      real(dp)::t_rad_HI, t_rad_HEI, t_rad_HEII
      real(dp)::t_rec_HI, t_rec_HEI, t_rec_HEII
      real(dp)::t_ion_HI, t_ion_HEI, t_ion_HEII
      real(dp)::t_ion2_HI, t_ion2_HEI, t_ion2_HEII
      real(dp)::x1, err_nE, err_nH2, n_H2_old
      integer, parameter::HI = 1, HeI = 2, HeII = 3
      real(dp)::g_H2, g_HI, g_HEI, g_HEII  ! Photoionisation/dissociation
      real(dp)::a_H2, a_HI, a_HEI, a_HEII  ! Recombination
      real(dp)::b_H2, b_HI, b_HEI, b_HEII  ! Collisional ionisation
      real(dp)::C_HII, D_H2, f_HII, f_H2    ! Creation and destruction
      real(dp)::D_HEI, C_HEIII, f_HeI, f_HeIII ! Creation and destruction
      real(dp)::min_xion = 1d-20                ! avoid zero division
!------------------------------------------------------------------------
      g_HEI = 0.; g_HEII = 0.; a_H2 = 0.; b_HEI = 0.; b_HEII = 0.; xe = 0.

      g_HI = t_rad_spec(ixHII)                !     Photoionization [s-1]
      if (isHe) then
         nHe = Y/(1.-Y)/4.*nH
         g_HEI = t_rad_spec(ixHeII)
         g_HEII = t_rad_spec(ixHeIII)
      end if
      if (isH2) then
         g_H2 = t_rad_spec(ixHI)                !    Photodissociation [s-1]
         a_H2 = comp_Alpha_H2(TK, Zsol, xe, 0.d0, 0.d0)     ! H2 creation rate [cm3 s-1]
      end if

      t_rad_HI = t_rad_spec(HI)               !      Photoionization [s-1]
      t_rad_HEI = t_rad_spec(HeI)
      t_rad_HEII = t_rad_spec(HeII)

      if (rt_OTSA) then                          !    Recombination [cm3 s-1]
         t_rec_HI = inp_coolrates_table(tbl_alphaB_HII, TK)
         t_rec_HEI = inp_coolrates_table(tbl_alphaB_HeII, TK)
         t_rec_HEII = inp_coolrates_table(tbl_alphaB_HeIII, TK)
      else
         t_rec_HI = inp_coolrates_table(tbl_alphaA_HII, TK)
         t_rec_HEI = inp_coolrates_table(tbl_alphaA_HeII, TK)
         t_rec_HEII = inp_coolrates_table(tbl_alphaA_HeIII, TK)
      end if

      a_HI = t_rec_HI
      a_HeI = t_rec_HEI
      a_HeII = t_rec_HEII

      !TODO isH2Katz=.true. uses different Coll. dissociation rate,
      !     but we decide to leave it for the moment,
      !     since this will only be called in the beginning of the simulation
      b_H2 = comp_Beta_H2HI(TK)                 ! Coll. dissociation [cm3 s-1]

      t_ion_HI = inp_coolrates_table(tbl_beta_HI, TK) ! Coll. ion. [cm3 s-1]
      if (isHe) t_ion_HEI = inp_coolrates_table(tbl_beta_HeI, TK)
      if (isHe) t_ion_HEII = inp_coolrates_table(tbl_beta_HeII, TK)

      b_HI = t_ion_HI
      b_HEI = t_ion_HEI
      b_HEII = t_ion_HEII

      n_E = nH; n_H2 = 0d0; n_H2_old = nH/2d0
      n_HeI = 0d0; n_HeII = 0d0; n_HeIII = 0d0
      err_nE = 1d0; err_nH2 = 0d0     ! err_nH2 initialisation in case of no H2
      do while (err_nE > 1d-8 .or. err_nH2 > 1d-8)

         n_E_min = MAX(n_E, 1e-15*nH)
         C_HII = b_HI*n_E_min + g_HI                 !   HII creation (s-1)
         f_HII = C_HII/a_HI/n_E_min                ! Cre/Destr [unitless]
         f_H2 = 0d0
         if (isH2) then
            D_H2 = b_H2*nH + g_H2                    !      H2 destr. (s-1)
            f_H2 = a_H2*nH/D_H2                    ! Cre/Destr [unitless]
            f_H2 = max(f_H2, min_xion)                 ! avoid zero division
            n_H2 = nH/(2d0 + 1d0/f_H2 + f_HII/f_H2)
         end if ! if(isH2)
         f_HII = max(f_HII, min_xion)                  ! avoid zero division
         n_HI = nH/(1d0 + f_HII + 2d0*f_H2)
         n_HII = nH/(1d0 + 1d0/f_HII + 2d0*f_H2/f_HII)

         if (isHe) then
            D_HeI = b_HEI*n_E_min + g_HEI           !  HeI destr. (s-1)
            C_HeIII = b_HEII*n_E_min + g_HEII          !  HeIII cre. (s-1)
            f_HeI = D_HeI/a_HeI/n_E_min          !  Destr/Cre [unitless]
            f_HeIII = a_HeII*n_E_min/C_HeIII       !  Destr/Cre [unitless]

            n_HEI = nHe/(1d0 + f_HeI + f_HeI/f_HeIII)
            n_HEII = nHe/(1d0 + 1d0/f_HeI + 1d0/f_HeIII)
            n_HEIII = nHe/(1d0 + f_HeIII + f_HeIII/f_HeI)
         end if ! if(isHe)

         err_nE = ABS((n_E - (n_HII + n_HEII + 2.*n_HEIII))/nH)
         n_E = 0.5*n_E + 0.5*(n_HII + n_HEII + 2.*n_HEIII)

         if (isH2) then
            err_nH2 = ABS((n_H2_old - n_H2)/nH)
            n_H2_old = n_H2
         end if

      end do

      nTOT = n_E + n_H2 + n_HI + n_HII + n_HEI + n_HEII + n_HEIII
      mu = nH/(1.-Y)/nTOT
      nSpec(1) = n_E
      nSpec(2) = n_H2
      nSpec(3) = n_HI
      nSpec(4) = n_HII
      nSpec(5) = n_HEI
      nSpec(6) = n_HEII
      nSpec(7) = n_HEIII

   END SUBROUTINE cmp_chem_eq

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE rt_evol_single_cell(astart, aend, dasura, h, omegab, omega0, omegaL &
                                  , J0min_in, T2end, mu, ne, if_write_result)
!-------------------------------------------------------------------------
! Used for initialization of thermal state in cosmological simulations.
!
! astart : valeur du facteur d'expansion au debut du calcul
! aend   : valeur du facteur d'expansion a la fin du calcul
! dasura : la valeur de da/a entre 2 pas de temps
! h      : la valeur de H0/100
! omegab : la valeur de Omega baryons
! omega0 : la valeur de Omega matiere (total)
! omegaL : la valeur de Omega Lambda
! J0min_in : la valeur du J0min a injecter :
!          Si high_z_realistic_ne alors c'est J0min a a=astart qui
!          est considere
!          Sinon, c'est le J0min habituel.
!          Si J0min_in <=0, les parametres par defaut ou predefinis
!          auparavant sont pris pour le J0min.
! T2end  : Le T/mu en output
! mu     : le poids moleculaire en output
! ne     : le ne en output
! if_write_result : .true. pour ecrire l'evolution de la temperature
!          et de n_e sur l'ecran.
!-------------------------------------------------------------------------
      use amr_commons, only: myid
      use UV_module
      implicit none
      real(kind=8)::astart, aend, T2end, h, omegab, omega0, omegaL, J0min_in, ne, dasura
      logical :: if_write_result
      real(dp)::aexp, daexp, dt_cool, coeff, T2_com, nH_com
      real(dp), dimension(nIons)::pHI_rates = 0., h_rad_spec = 0.
      real(kind=8) ::mu
      real(dp) ::cool_tot, heat_tot, mu_dp, diff, Zini
      integer::niter
      real(dp) :: n_spec(1:6)
      real(dp), dimension(1:nvector):: T2
      real(dp), dimension(1:nIons, 1:nvector):: xion
      real(dp), dimension(1:n_oxygen_ions, 1:nvector):: xO
      real(dp), dimension(1:n_nitrogen_ions, 1:nvector):: xN
      real(dp), dimension(1:n_carbon_ions, 1:nvector):: xC
      real(dp), dimension(1:n_magnesium_ions, 1:nvector):: xMg
      real(dp), dimension(1:n_silicon_ions, 1:nvector):: xSi
      real(dp), dimension(1:n_sulfur_ions, 1:nvector):: xS
      real(dp), dimension(1:n_iron_ions, 1:nvector):: xFe
      real(dp), dimension(1:n_neon_ions, 1:nvector):: xNe
      real(dp), dimension(1:nGroups, 1:nvector):: Np, dNpdt
      real(dp), dimension(1:ndim, 1:nGroups, 1:nvector):: Fp, dFpdt
      real(dp), dimension(1:ndim, 1:nvector):: p_gas
      real(dp), dimension(1:nvector)::nH = 0., Zsolar = 0., dx_SS = 0., rate=0., tsave=0.
      real(dp), dimension(ncoolheatratesave,1:nvector)::csave
      real(dp), dimension(1:nvector)::nOxygen = 0., nNitrogen = 0., nCarbon = 0., nMagnesium = 0.
      real(dp), dimension(1:nvector)::nSilicon = 0., nSulfur = 0., nIron = 0., nNeon = 0.
      real(dp), dimension(1:nvector)::nco_mol = 0.
      logical, dimension(1:nvector)::c_switch = .true.
!-------------------------------------------------------------------------
      aexp = astart
      T2_com = 2.726d0/aexp*aexp**2/mu_mol
      nH_com = omegab*rhoc*h**2*X/mH

      Zini = z_ave*0.02

      mu_dp = mu
      call cmp_Equilibrium_Abundances( &
         T2_com/aexp**2, nH_com/aexp**3, pHI_rates, mu_dp, n_Spec, Zini)

      ! Initialize cell state
      T2(1) = T2_com                                          !      Temperature
      xion(1, 1) = n_Spec(3)/(nH_com/aexp**3)                  !   HII   fraction
      xion(2, 1) = n_Spec(5)/(nH_com/aexp**3)                  !   HeII  fraction
      xion(3, 1) = n_Spec(6)/(nH_com/aexp**3)                  !   HeIII fraction
      p_gas(:, 1) = 0.
      Np(:, 1) = 0.; Fp(:, :, 1) = 0.                  ! Photon densities and fluxes
      dNpdt(:, 1) = 0.; dFpdt(:, :, 1) = 0.

      ! Initialize all metals to a neutral state
      if (oxygen_ions) then
         xO = 1.d-10
         xO(1,:) = 1.d0
      endif

      if (nitrogen_ions) then
         xN = 1.d-10
         xN(1,:) = 1.d0
      endif

      if (carbon_ions) then
         xC = 1.d-10
         xC(1,:) = 1.d0
      endif

      if (magnesium_ions) then
         xMg = 1.d-10
         xMg(1,:) = 1.d0
      endif

      if (silicon_ions) then
         xSi = 1.d-10
         xSi(1,:) = 1.d0
      endif

      if (sulfur_ions) then
         xS = 1.d-10
         xS(1,:) = 1.d0
      endif

      if (iron_ions) then
         xFe = 1.d-10
         xFe(1,:) = 1.d0
      endif

      if (neon_ions) then
         xNe = 1.d-10
         xNe(1,:) = 1.d0
      endif

      do while (aexp < aend)
         call update_UVrates(aexp)
         call update_coolrates_tables(aexp)
         if (uvbg_rtz) call update_UV_background_RTZ(1.d0/aexp - 1.d0)

         daexp = dasura*aexp
         dt_cool = daexp &
                   /(aexp*100.*h*3.2408608e-20) &
                   /HsurH0(1.0/dble(aexp) - 1., omega0, omegaL, 1.-omega0 - omegaL)

         nH(1) = nH_com/aexp**3
         T2(1) = T2(1)/aexp**2
         call rt_solve_cooling(T2, xion, Np, Fp, p_gas, dNpdt, dFpdt, nH, nH*mH, c_switch &
                               , Zsolar, dt_cool, aexp, 1, dx_SS &
                               , xO, xN, xC, xMg, xSi, xS, xFe, xNe, nOxygen, nNitrogen & 
                               , nCarbon, nMagnesium, nSilicon, nSulfur, nIron, nNeon   &
                               , nco_mol, 0.d0, rate &
                               , tsave, csave) ! dx_SS=0.; no self-shielding
                               ! we set the final value 0 since zsolar =0. This sets the dust and gas mass to 0
         T2(1) = T2(1)*aexp**2
         aexp = aexp + daexp
         if (if_write_result) write (*, '(4(1pe10.3))') &
            aexp, nH(1), T2_com*mu/aexp**2, n_spec(1)/nH(1)
      end do
      T2end = T2(1)/(aexp - daexp)**2
      ne = (n_spec(3) + (n_spec(5) + 2.*n_spec(6))*0.25*Y/X)
   end subroutine rt_evol_single_cell

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   FUNCTION HsurH0(z, omega0, omegaL, OmegaR)
!-------------------------------------------------------------------------
      implicit none
      real(kind=8) :: HsurH0, z, omega0, omegaL, omegaR
!-------------------------------------------------------------------------
      HsurH0 = sqrt(Omega0*(1.d0 + z)**3 + OmegaR*(1.d0 + z)**2 + OmegaL)
   END FUNCTION HsurH0

   subroutine rt_cooling_fire(T2, mu, nH, Z, gH, xH, ne, fire_tot)
      !T2 ==> Temperature in kelvin, divided by mu
      !nH ==> Hydrogen number density H/cc
      !mu ==> Average mass per particle in terms of mH
      !Z  ==> Cell metallicity as a fraction of solar
      !gH ==> Photoionization rate of neutral hydrogen s^-1
      !ne ==> Electron number density e/cc / Hydrogen number density
      !xH ==> Neutral hydrogen fraction
      use pm_commons, only: nstar_tot
      implicit none
      real(KIND=8)::T2, mu, nH, Z, gH, xH, ne, fire_tot
      real(KIND=8)::TT, lTT, lcold, ldust, lcr
      real(KIND=8)::sigv, pc2cm, la, taua, fss, Tdust, ecr, xpe, netild, enu, epc

      TT = max(T2*mu, 2.725) !Temperature in Kelvin
      lTT = log10(TT) !Log of the temeprature in Kelvin
      gh = max(gH, 1.d-20)            !set lower limit on gH to prevent divide by 0
      xH = min(max(xH, 1.d-6), 1.d0)  !similar for xH
      netild = max(ne/nH, 1.d-6) !electron fraction

      if (ltt .le. 1.d4) then
         !Low temperature metal cooling (eqn. B18 in https://arxiv.org/pdf/1702.06148.pdf)
         sigv = 6.0d-18
         pc2cm = 3.086d18
         la = 4.4*pc2cm*((TT/10000.0)**(-0.173))*((gH/(1.0d-12))**(-0.66666))
         taua = sigv*nH*la
         fss = exp(-1.0*taua)
         lcold = (2.896d-26)*((((TT/125.215)**(-4.9202)) + ((TT/1349.86)**(-1.7288)) &
                               + ((TT/6450.06)**(-0.3075)))**(-1.0)) &
                 *((1.0 + Z)/(1.0 + (0.00143*nH)))*(1.0 - fss) &
                 *(0.001 + ((0.1*nH)/(1.0 + nH)) + ((0.09*nH)/(1.0 + 0.1*nH)) &
                   + (Z*Z/(1.0 + nH)))*exp(-1.0*(TT/158000.0)*(TT/158000.0))

         !Dust cooling (eqn. B19 https://arxiv.org/pdf/1702.06148.pdf)
         !There is a slight mistake in hopkins et al...it should be 1.2 instead of 1.12
         !This is also only for T<10^4...see hollenbach 1989
         Tdust = 30.0 !Temperature of the dust in K
         ldust = (1.2d-32)*(TT - Tdust)*sqrt(TT)*(1.0 - (0.8*exp(-75.0/TT)))*Z
      else
         !If not in low temperature regime, set cooling rate to a very low number
         lcold = 1.d-40
         ldust = 1.d-40
      end if

      if (nstar_tot .gt. 0) then
         !Cosmic ray heating (eqn. B24 in https://arxiv.org/pdf/1702.06148.pdf)
         ecr = (9.0d-12)*(min(nH, 0.01)/0.01) !Flat energy at nH>0.01 otherwise decrease linearly
         lcr = (-1.0d-16)*(0.98 + (1.65*netild*xH))*ecr*(1.0/nH)
      else
         lcr = -1.d-40
      end if

      fire_tot = lcold + ldust + lcr

   end subroutine rt_cooling_fire

!=========================================================================
   subroutine rt_cmp_metals(T2, nH, mu, metal_tot, metal_prime, aexp)
! Taken from the equilibrium cooling_module of RAMSES
! Compute cooling enhancement due to metals
! T2           => Temperature in Kelvin, divided by mu
! nH           => Hydrogen number density (H/cc)
! mu           => Average mass per particle in terms of mH
! metal_tot   <=  Metal cooling contribution to de/dt [erg s-1 cm-3]
! metal_prime <=  d(metal_tot)/dT2 [erg s-1 cm-3 K-1]
!=========================================================================
      implicit none
      real(dp) ::T2, nH, mu, metal_tot, metal_prime, aexp
      ! Cloudy at solar metalicity
      real(dp), dimension(1:91), parameter :: temperature_cc07 = (/ &
           & 3.9684, 4.0187, 4.0690, 4.1194, 4.1697, 4.2200, 4.2703, &
           & 4.3206, 4.3709, 4.4212, 4.4716, 4.5219, 4.5722, 4.6225, &
           & 4.6728, 4.7231, 4.7734, 4.8238, 4.8741, 4.9244, 4.9747, &
           & 5.0250, 5.0753, 5.1256, 5.1760, 5.2263, 5.2766, 5.3269, &
           & 5.3772, 5.4275, 5.4778, 5.5282, 5.5785, 5.6288, 5.6791, &
           & 5.7294, 5.7797, 5.8300, 5.8804, 5.9307, 5.9810, 6.0313, &
           & 6.0816, 6.1319, 6.1822, 6.2326, 6.2829, 6.3332, 6.3835, &
           & 6.4338, 6.4841, 6.5345, 6.5848, 6.6351, 6.6854, 6.7357, &
           & 6.7860, 6.8363, 6.8867, 6.9370, 6.9873, 7.0376, 7.0879, &
           & 7.1382, 7.1885, 7.2388, 7.2892, 7.3395, 7.3898, 7.4401, &
           & 7.4904, 7.5407, 7.5911, 7.6414, 7.6917, 7.7420, 7.7923, &
           & 7.8426, 7.8929, 7.9433, 7.9936, 8.0439, 8.0942, 8.1445, &
           & 8.1948, 8.2451, 8.2955, 8.3458, 8.3961, 8.4464, 8.4967/)
      ! Cooling from metals only (without the contribution of H and He)
      ! log cooling rate in [erg s-1 cm3]
      ! S. Ploeckinger 06/2015
      real(kind=8), dimension(1:91) :: excess_cooling_cc07 = (/ &
           &  -24.9082, -24.9082, -24.5503, -24.0898, -23.5328, -23.0696, -22.7758, &
           &  -22.6175, -22.5266, -22.4379, -22.3371, -22.2289, -22.1181, -22.0078, &
           &  -21.8992, -21.7937, -21.6921, -21.5961, -21.5089, -21.4343, -21.3765, &
           &  -21.3431, -21.3274, -21.3205, -21.3142, -21.3040, -21.2900, -21.2773, &
           &  -21.2791, -21.3181, -21.4006, -21.5045, -21.6059, -21.6676, -21.6877, &
           &  -21.6934, -21.7089, -21.7307, -21.7511, -21.7618, -21.7572, -21.7532, &
           &  -21.7668, -21.7860, -21.8129, -21.8497, -21.9035, -21.9697, -22.0497, &
           &  -22.1327, -22.2220, -22.3057, -22.3850, -22.4467, -22.4939, -22.5205, &
           &  -22.5358, -22.5391, -22.5408, -22.5408, -22.5475, -22.5589, -22.5813, &
           &  -22.6122, -22.6576, -22.7137, -22.7838, -22.8583, -22.9348, -23.0006, &
           &  -23.0547, -23.0886, -23.1101, -23.1139, -23.1147, -23.1048, -23.1017, &
           &  -23.0928, -23.0969, -23.0968, -23.1105, -23.1191, -23.1388, -23.1517, &
           &  -23.1717, -23.1837, -23.1986, -23.2058, -23.2134, -23.2139, -23.2107/)
      real(dp), dimension(1:91), parameter :: excess_prime_cc07 = (/           &
           &   2.0037, 4.7267, 12.2283, 13.5820, 9.8755, 4.8379, 1.8046, &
           &   1.4574, 1.8086, 2.0685, 2.2012, 2.2250, 2.2060, 2.1605, &
           &   2.1121, 2.0335, 1.9254, 1.7861, 1.5357, 1.1784, 0.7628, &
           &   0.1500, -0.1401, 0.1272, 0.3884, 0.2761, 0.1707, 0.2279, &
           &  -0.2417, -1.7802, -3.0381, -2.3511, -0.9864, -0.0989, 0.1854, &
           &  -0.1282, -0.8028, -0.7363, -0.0093, 0.3132, 0.1894, -0.1526, &
           &  -0.3663, -0.3873, -0.3993, -0.6790, -1.0615, -1.4633, -1.5687, &
           &  -1.7183, -1.7313, -1.8324, -1.5909, -1.3199, -0.8634, -0.5542, &
           &  -0.1961, -0.0552, 0.0646, -0.0109, -0.0662, -0.2539, -0.3869, &
           &  -0.6379, -0.8404, -1.1662, -1.3930, -1.6136, -1.5706, -1.4266, &
           &  -1.0460, -0.7244, -0.3006, -0.1300, 0.1491, 0.0972, 0.2463, &
           &   0.0252, 0.1079, -0.1893, -0.1033, -0.3547, -0.2393, -0.4280, &
           &  -0.2735, -0.3670, -0.2033, -0.2261, -0.0821, -0.0754, 0.0634/)
      real(dp), dimension(1:50), parameter::z_courty = (/                         &
           & 0.00000, 0.04912, 0.10060, 0.15470, 0.21140, 0.27090, 0.33330, 0.39880, &
           & 0.46750, 0.53960, 0.61520, 0.69450, 0.77780, 0.86510, 0.95670, 1.05300, &
           & 1.15400, 1.25900, 1.37000, 1.48700, 1.60900, 1.73700, 1.87100, 2.01300, &
           & 2.16000, 2.31600, 2.47900, 2.64900, 2.82900, 3.01700, 3.21400, 3.42100, &
           & 3.63800, 3.86600, 4.10500, 4.35600, 4.61900, 4.89500, 5.18400, 5.48800, &
           & 5.80700, 6.14100, 6.49200, 6.85900, 7.24600, 7.65000, 8.07500, 8.52100, &
           & 8.98900, 9.50000/)
      real(dp), dimension(1:50), parameter::phi_courty = (/                             &
           & 0.0499886, 0.0582622, 0.0678333, 0.0788739, 0.0915889, 0.1061913, 0.1229119, &
           & 0.1419961, 0.1637082, 0.1883230, 0.2161014, 0.2473183, 0.2822266, 0.3210551, &
           & 0.3639784, 0.4111301, 0.4623273, 0.5172858, 0.5752659, 0.6351540, 0.6950232, &
           & 0.7529284, 0.8063160, 0.8520859, 0.8920522, 0.9305764, 0.9682031, 1.0058810, &
           & 1.0444020, 1.0848160, 1.1282190, 1.1745120, 1.2226670, 1.2723200, 1.3231350, &
           & 1.3743020, 1.4247480, 1.4730590, 1.5174060, 1.5552610, 1.5833640, 1.5976390, &
           & 1.5925270, 1.5613110, 1.4949610, 1.3813710, 1.2041510, 0.9403100, 0.5555344, &
           & 0.0000000/)
      real(dp)::TT, lTT, deltaT, lcool, lcool1, lcool2, lcool1_prime, lcool2_prime
      real(dp)::ZZ, deltaZ
      real(dp)::c1 = 0.4, c2 = 10.0, TT0 = 1d5, TTC = 1d6, alpha1 = 0.15
      real(dp)::ux, g_courty, f_courty = 1d0, g_courty_prime, f_courty_prime
      integer::iT, iZ
!-------------------------------------------------------------------------
      ZZ = 1d0/aexp - 1d0
      TT = T2*mu
      lTT = log10(TT)

      ! This is a simple model to take into account the ionization background
      ! on metal cooling (calibrated using CLOUDY).
      iZ = 1 + int(ZZ/z_courty(50)*49.)
      iZ = min(iZ, 49)
      iZ = max(iZ, 1)
      deltaZ = z_courty(iZ + 1) - z_courty(iZ)
      ZZ = min(ZZ, z_courty(50))
      ZZ = max(ZZ, z_courty(1))
      ux = 1d-4*(phi_courty(iZ + 1)*(ZZ - z_courty(iZ))/deltaZ &
           & + phi_courty(iZ)*(z_courty(iZ + 1) - ZZ)/deltaZ)/nH
      g_courty = c1*(TT/TT0)**alpha1 + c2*exp(-TTC/TT)
      g_courty_prime = (c1*alpha1*(TT/TT0)**alpha1 + c2*exp(-TTC/TT)*TTC/TT)/TT
      f_courty = 1d0/(1d0 + ux/g_courty)
      f_courty_prime = ux/g_courty/(1d0 + ux/g_courty)**2*g_courty_prime/g_courty

      if (lTT .ge. temperature_cc07(91)) then
         metal_tot = 0d0 !1d-100
         metal_prime = 0d0
      else if (lTT .ge. 1.0) then
         lcool1 = -100d0
         lcool1_prime = 0d0
         if (lTT .ge. temperature_cc07(1)) then
            iT = 1 + int((lTT - temperature_cc07(1))/ &
                         (temperature_cc07(91) - temperature_cc07(1))*90.0)
            iT = min(iT, 90)
            iT = max(iT, 1)
            deltaT = temperature_cc07(iT + 1) - temperature_cc07(iT)
            lcool1 = &
               excess_cooling_cc07(iT + 1)*(lTT - temperature_cc07(iT))/deltaT &
               + excess_cooling_cc07(iT)*(temperature_cc07(iT + 1) - lTT)/deltaT
            lcool1_prime = &
               excess_prime_cc07(iT + 1)*(lTT - temperature_cc07(iT))/deltaT &
               + excess_prime_cc07(iT)*(temperature_cc07(iT + 1) - lTT)/deltaT
         end if
         ! We only want this when not using fire cooling
         lcool2 = 0d0
         lcool2_prime = 0d0
         if (.not. fire_cooling .and. .not. low_t_noneq_cooling) then
            ! Fine structure cooling from infrared lines
            lcool2 = -31.522879 + 2.0*lTT - 20.0/TT - TT*4.342944d-5
            lcool2_prime = 2d0 + (20d0/TT - TT*4.342944d-5)*log(10d0)
            ! Total metal cooling and temperature derivative
            metal_tot = 10d0**lcool1 + 10d0**lcool2
            metal_prime = (10d0**lcool1*lcool1_prime + 10d0**lcool2*lcool2_prime)/metal_tot
         else
            metal_tot = 10d0**lcool1
            metal_prime = (10d0**lcool1*lcool1_prime)/metal_tot
         end if
         metal_prime = metal_prime*f_courty + metal_tot*f_courty_prime
         metal_tot = metal_tot*f_courty
      else
         metal_tot = 0d0 !1d-100
         metal_prime = 0d0
      end if

      metal_tot = metal_tot*nH**2
      metal_prime = &   ! Convert from DlogLambda/DlogT to DLambda/DT
         metal_prime*metal_tot/TT*mu

   end subroutine rt_cmp_metals
!=========================================================================
   subroutine rt_cmp_metals_cloudy(T2, nH, mu, metal_tot, metal_prime, aexp)
! Taken from the equilibrium cooling from Cloudy
! Compute cooling enhancement due to metals
! T2           => Temperature in Kelvin, divided by mu
! nH           => Hydrogen number density (H/cc)
! mu           => Average mass per particle in terms of mH
! metal_tot   <=  Metal cooling contribution to de/dt [erg s-1 cm-3]
! metal_prime <=  d(metal_tot)/dT2 [erg s-1 cm-3 K-1]
!=========================================================================
      implicit none
      real(dp) ::T2, nH, mu, metal_tot, metal_prime, aexp
      real(dp) ::metal_heat, metal_heat_prime, TT

      ! If fire low temperature cooling, only use grackle tables at
      ! T>10^4 K
      if ((fire_cooling.or.low_t_noneq_cooling) .and. T2*mu .lt. 1.d4) then
         metal_tot = 1d-100
         metal_prime = 0d0
         metal_heat = 1d-100
         metal_heat_prime = 0d0
      else
         call cmp_metals_cloudy(T2, nH, mu, metal_tot, metal_prime, metal_heat, metal_heat_prime, aexp)
      end if

      if (metal_tot < 1d-99) metal_tot = 0d0

      metal_tot = metal_tot - metal_heat
      metal_prime = metal_prime - metal_heat_prime

      TT = T2*mu
      metal_tot = metal_tot*(nH**2)
      metal_prime = &   ! Convert from DlogLambda/DlogT to DLambda/DT
         metal_prime*metal_tot/TT*mu

   end subroutine rt_cmp_metals_cloudy
!*************************************************************************
   FUNCTION getMu(xHII, xHeII, xHeIII, Tmu)
! Returns the mean particle mass, in units of the proton mass.
! xHII, xHeII, xHeIII => Hydrogen and helium ionisation fractions
! Tmu => T/mu in Kelvin
!-------------------------------------------------------------------------
      implicit none
      real(kind=8), intent(in) :: xHII, xHeII, xHeIII, Tmu
      real(kind=8) :: mu
      real(kind=8) :: getMu
!-------------------------------------------------------------------------
      getMu = 1./(X*(1.+xHII) + 0.25*Y*(1.+xHeII + 2.*xHeIII))
      if (is_kIR_T .or. is_mu_H2) &
         getMu = getMu + exp(-1.d0*(Tmu/Tmu_dissoc)**2)*(2.33 - getMu)
   END FUNCTION getMu
END MODULE rt_cooling_module

!************************************************************************
SUBROUTINE updateRTGroups_CoolConstants()
! Update photon group cooling and heating constants, to reflect an update
! in rt_c_cgs and in the cross-sections and energies in the groups.
!------------------------------------------------------------------------
   use rt_cooling_module
   use rt_parameters
   implicit none
   integer::iP, iI
!------------------------------------------------------------------------
   signc = group_csn*rt_c_cgs                                    ! [cm3 s-1]
   sigec = group_cse*rt_c_cgs                                    ! [cm3 s-1]

   signc_dust = group_csn_dust*rt_c_cgs                          ! [cm3 s-1]

   if (oxygen_ions) then 
      signc_oxygen = group_csn_oxygen*rt_c_cgs                   ! [cm3 s-1]
      sigec_oxygen = group_cse_oxygen*rt_c_cgs                   ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_oxygen_ions               ! Photoheating rates for photons on ions
            PHrate_oxygen(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_oxygen(iP, iI)*group_egy(iP) - signc_oxygen(iP, iI)*oxygen_ionEvs(iI))
            PHrate_oxygen(iP, iI) = max(PHrate_oxygen(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (nitrogen_ions) then 
      signc_nitrogen = group_csn_nitrogen*rt_c_cgs               ! [cm3 s-1]
      sigec_nitrogen = group_cse_nitrogen*rt_c_cgs               ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_nitrogen_ions               ! Photoheating rates for photons on ions
            PHrate_nitrogen(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_nitrogen(iP, iI)*group_egy(iP) - signc_nitrogen(iP, iI)*nitrogen_ionEvs(iI))
            PHrate_nitrogen(iP, iI) = max(PHrate_nitrogen(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (carbon_ions) then 
      signc_carbon = group_csn_carbon*rt_c_cgs                   ! [cm3 s-1]
      sigec_carbon = group_cse_carbon*rt_c_cgs                   ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_carbon_ions               ! Photoheating rates for photons on ions
            PHrate_carbon(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_carbon(iP, iI)*group_egy(iP) - signc_carbon(iP, iI)*carbon_ionEvs(iI))
            PHrate_carbon(iP, iI) = max(PHrate_carbon(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (magnesium_ions) then 
      signc_magnesium = group_csn_magnesium*rt_c_cgs             ! [cm3 s-1]
      sigec_magnesium = group_cse_magnesium*rt_c_cgs             ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_magnesium_ions               ! Photoheating rates for photons on ions
            PHrate_magnesium(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_magnesium(iP, iI)*group_egy(iP) - signc_magnesium(iP, iI)*magnesium_ionEvs(iI))
            PHrate_magnesium(iP, iI) = max(PHrate_magnesium(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (silicon_ions) then 
      signc_silicon = group_csn_silicon*rt_c_cgs                 ! [cm3 s-1]
      sigec_silicon = group_cse_silicon*rt_c_cgs                 ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_silicon_ions               ! Photoheating rates for photons on ions
            PHrate_silicon(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_silicon(iP, iI)*group_egy(iP) - signc_silicon(iP, iI)*silicon_ionEvs(iI))
            PHrate_silicon(iP, iI) = max(PHrate_silicon(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (sulfur_ions) then 
      signc_sulfur = group_csn_sulfur*rt_c_cgs                   ! [cm3 s-1]
      sigec_sulfur = group_cse_sulfur*rt_c_cgs                   ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_sulfur_ions               ! Photoheating rates for photons on ions
            PHrate_sulfur(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_sulfur(iP, iI)*group_egy(iP) - signc_sulfur(iP, iI)*sulfur_ionEvs(iI))
            PHrate_sulfur(iP, iI) = max(PHrate_sulfur(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (iron_ions) then 
      signc_iron = group_csn_iron*rt_c_cgs                       ! [cm3 s-1]
      sigec_iron = group_cse_iron*rt_c_cgs                       ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_iron_ions               ! Photoheating rates for photons on ions
            PHrate_iron(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_iron(iP, iI)*group_egy(iP) - signc_iron(iP, iI)*iron_ionEvs(iI))
            PHrate_iron(iP, iI) = max(PHrate_iron(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   if (neon_ions) then 
      signc_neon = group_csn_neon*rt_c_cgs                       ! [cm3 s-1]
      sigec_neon = group_cse_neon*rt_c_cgs                       ! [cm3 s-1]
      do iP = 1, nGroups
         do iI = 1, n_neon_ions               ! Photoheating rates for photons on ions
            PHrate_neon(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                           (sigec_neon(iP, iI)*group_egy(iP) - signc_neon(iP, iI)*neon_ionEvs(iI))
            PHrate_neon(iP, iI) = max(PHrate_neon(iP, iI), 0d0) !      No negative heating
         end do
      end do
   endif

   do iP = 1, nGroups
      do iI = 1, nIons               ! Photoheating rates for photons on ions
         PHrate(iP, iI) = ev_to_erg* &        ! See eq (19) in Aubert(08)
                          (sigec(iP, iI)*group_egy(iP) - signc(iP, iI)*ionEvs(iI))
         PHrate(iP, iI) = max(PHrate(iP, iI), 0d0) !      No negative heating
      end do
   end do
END SUBROUTINE updateRTGroups_CoolConstants

!************************************************************************
SUBROUTINE reduce_flux(Fp, cNp)
! Make sure the reduced photon flux is less than one
!------------------------------------------------------------------------
   use rt_parameters
   implicit none
   real(dp), dimension(ndim):: Fp
   real(dp):: cNp, fred
!------------------------------------------------------------------------
   fred = sqrt(sum(Fp**2))/cNp
   if (fred .gt. 1.d0) Fp = Fp/fred
END SUBROUTINE reduce_flux

!************************************************************************
SUBROUTINE reduce_xion(xion, x_MIN, x_MIN_HI)
! Make sure the total ionization fraction is one (Taysun)
!------------------------------------------------------------------------
   use amr_commons, ONLY: dp
   use rt_parameters, ONLY: nIons, isHe, isH2, ixHI, ixHII, ixHeII, ixHeIII
   implicit none
   real(kind=dp), dimension(1:nIons)::xion
   real(kind=dp)::x_MIN, x_MIN_HI
!------------------------------------------------------------------------
   ! xion(ixHI )   =    nHI  ---> this could be confusing, but followed the original version
   ! xion(ixHII)   =   nHII
   ! xion(ixHeII)  =  nHeII
   ! xion(ixHeIII) = nHeIII

   xion(ixHII:ixHeIII) = MIN(MAX(xion(ixHII:ixHeIII), x_MIN), 1d0)

   !Ensure the total hydrogen fraction is 1
   if (isH2) then
      xion(ixHI) = MIN(MAX(xion(ixHI), x_MIN_HI), 1d0)
      if (xion(ixHI) + xion(ixHII) .gt. 1.d0) then
         if (xion(ixHI) .gt. xion(ixHII)) then
            xion(ixHI) = 1.d0 - xion(ixHII)
         else
            xion(ixHII) = 1.d0 - xion(ixHI) 
         end if
      end if
   end if

   !Ensure the total helium fraction is 1
   if (isHe) then
      if (xion(ixHeII) + xion(ixHeIII) .gt. 1.d0) then
         if (xion(ixHeII) .gt. xion(ixHeIII)) then
            xion(ixHeII) = 1.d0 - xion(ixHeIII)
         else
            xion(ixHeIII) = 1.d0 - xion(ixHeII)
         end if
      end if
   end if
END SUBROUTINE reduce_xion

SUBROUTINE get_dust_mass_and_depletion(depletion_factors,dx_loc,d_mass,g_mass,nH,nFe,nO,nN,nMg,nNe,nSi,nC,nS)
   use amr_commons, only: dp
   use hydro_parameters, only: nmetals
   use rt_cooling_module, only: amu_to_g, mO_NIST_amu, mN_NIST_amu, mC_NIST_amu, mMg_NIST_amu, mSi_NIST_amu, mS_NIST_amu, mFe_NIST_amu, mNe_NIST_amu
   implicit none
   real(kind=dp),INTENT(IN)::nH,nFe,nO,nN,nMg,nNe,nSi,nC,nS,g_mass,dx_loc
   ! Depletion: Hardcoded. 1=Fe, 2=O, 3=N, 4=Mg, 5=Ne, 6=Si, 7=Ca, 8=C, 9=S
   real(kind=dp), dimension(1:nmetals), INTENT(OUT)::depletion_factors 
   real(kind=dp),INTENT(OUT)::d_mass
   ! Params for the RR14 dust to gas mass ratio
   real(kind=dp)::a,aH,b,aL,xt,xs,x,y
   ! Params for the Zubko 2004 BARE-GR-S model
   real(kind=dp)::C_mfrac,O_mfrac,Si_mfrac,Mg_mfrac,Fe_mfrac

   ! Broken powerlaw model from RR14 (consistent with Taysun's Lya feedback)
   ! See Table 1 of: https://www.aanda.org/articles/aa/pdf/2014/03/aa22803-13.pdf
   ! We use the XCO,Z case

   a  = 2.21d0
   aH = 1.00d0
   b  = 0.96d0
   aL = 3.10d0
   xt = 8.10d0
   xs = 8.69d0

   ! NOTE: no limit on metallicity here! otherwise 
   ! very low metallicities will have high depletions
   x = 12.d0 + LOG10( MAX(nO,1.d-40) / nH)

   if (x.gt.xt) then
      y = a + (aH * (xs - x))
   else
      y = b + (aL * (xs - x))
   end if

   y = 10.d0**y ! This is the Gas to Dust mass ratio
   d_mass = g_mass / y ! Get the total dust mass

   ! Set the depletion factors 
   ! This is a mix from the prediction of the BARE-GR-S Model and cloudy
   ! when bare-gr-s doesn't make sense due to negative values
   
   ! Hardcoded. 1=Fe, 2=O, 3=N, 4=Mg, 5=Ne, 6=Si, 7=Ca, 8=C, 9=S
   ! depletion_factors(1) = 1.00d-2 * MAX(1.d0,y/162.d0)  ! Iron
   ! depletion_factors(2) = 0.728d0 * MAX(1.d0,y/162.d0)  ! Oxygen - bare-gr-s
   ! depletion_factors(3) = 1.000d0 * MAX(1.d0,y/162.d0)  ! Nitrogen
   ! depletion_factors(4) = 0.163d0 * MAX(1.d0,y/162.d0)  ! Magnesium - bare-gr-s
   ! depletion_factors(5) = 1.000d0 * MAX(1.d0,y/162.d0)  ! Neon
   ! depletion_factors(6) = 0.030d0 * MAX(1.d0,y/162.d0)  ! Silicon
   ! depletion_factors(7) = 1.00d-4 * MAX(1.d0,y/162.d0)  ! Calcium
   ! depletion_factors(8) = 0.400d0 * MAX(1.d0,y/162.d0)  ! Carbon
   ! depletion_factors(9) = 1.000d0 * MAX(1.d0,y/162.d0)  ! Sulfur

   ! Set the depletion factors 
   ! This is a mix from the prediction of the BARE-GR-S Model and dopita 2000
   ! when bare-gr-s doesn't make sense due to negative values
   
   ! Hardcoded. 1=Fe, 2=O, 3=N, 4=Mg, 5=Ne, 6=Si, 7=Ca, 8=C, 9=
   depletion_factors(1) = 1.d0 - ( (1.d0 - 1.00d-2) * MIN(1.d0,162.d0/y) )  ! Iron
   depletion_factors(2) = 1.d0 - ( (1.d0 - 0.728d0) * MIN(1.d0,162.d0/y) ) ! Oxygen - bare-gr-s
   depletion_factors(3) = 1.d0 - ( (1.d0 - 0.603d0) * MIN(1.d0,162.d0/y) ) ! Nitrogen
   depletion_factors(4) = 1.d0 - ( (1.d0 - 0.163d0) * MIN(1.d0,162.d0/y) ) ! Magnesium - bare-gr-s
   depletion_factors(5) = 1.d0 - ( (1.d0 - 1.000d0) * MIN(1.d0,162.d0/y) ) ! Neon
   depletion_factors(6) = 1.d0 - ( (1.d0 - 0.100d0) * MIN(1.d0,162.d0/y) ) ! Silicon
   depletion_factors(7) = 1.d0 - ( (1.d0 - 0.003d0) * MIN(1.d0,162.d0/y) ) ! Calcium
   depletion_factors(8) = 1.d0 - ( (1.d0 - 0.501d0) * MIN(1.d0,162.d0/y) ) ! Carbon
   depletion_factors(9) = 1.d0 - ( (1.d0 - 1.000d0) * MIN(1.d0,162.d0/y) ) ! Sulfur

END SUBROUTINE get_dust_mass_and_depletion
