!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedback(ilevel)
   use pm_commons
   use amr_commons
   use tracer_utils, only: pre_particle_yield, post_particle_yield, yield_tracers

   implicit none
   integer::ilevel
   !------------------------------------------------------------------------
   ! This routine computes the thermal energy, the kinetic energy and
   ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
   ! This routine is called every fine time step.
   !------------------------------------------------------------------------
   real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
   real(dp)::t0, scale, dx_min, vsn, rdebris, ethermal
   integer::igrid, jgrid, ipart, jpart, next_part
   integer::i, ig, ip, npart1, npart2, icpu, nx_loc
   real(dp), dimension(1:3)::skip_loc
   integer, dimension(1:nvector), save::ind_grid, ind_part, ind_grid_part
   logical::ok_star, ok_star3

   ! MC tracer
   type(part_t) :: star_tracer_type
   star_tracer_type%family = FAM_TRACER_STAR

   if (numbtot(1, ilevel) == 0) return
   if (verbose) write (*, 111) ilevel

   if (MC_tracer) call pre_particle_yield()

   ! Gather star particles only.
#if NDIM==3
   ! Loop over cpus
   do icpu = 1, ncpu
      igrid = headl(icpu, ilevel)
      ig = 0
      ip = 0
      ! Loop over grids
      do jgrid = 1, numbl(icpu, ilevel)
         npart1 = numbp(igrid)  ! Number of particles in the grid
         npart2 = 0

         ! Count active star particles
         if (npart1 > 0) then
            ipart = headp(igrid)
            ! Loop over particles
            do jpart = 1, npart1
               ! Save next particle   <--- Very important !!!
               next_part = nextp(ipart)
               ok_star3 = is_pop_III(typep(ipart)) .and. is_pre_SN(typep(ipart))
               if (is_pop_II(typep(ipart)) .or. ok_star3) then
                  npart2 = npart2 + 1
               end if
               ipart = next_part  ! Go to next particle
            end do
         end if

         ! Gather star particles
         if (npart2 > 0) then
            ig = ig + 1
            ind_grid(ig) = igrid
            ipart = headp(igrid)
            ! Loop over particles
            do jpart = 1, npart1
               ! Save next particle   <--- Very important !!!
               next_part = nextp(ipart)
               ! Select only star particles
               ok_star3 = is_pop_III(typep(ipart)) .and. is_pre_SN(typep(ipart))
               ok_star = is_pop_II(typep(ipart)) .or. ok_star3

               if (ok_star) then
                  if (ig == 0) then
                     ig = 1
                     ind_grid(ig) = igrid
                  end if
                  ip = ip + 1
                  ind_part(ip) = ipart
                  ind_grid_part(ip) = ig
               end if
               if (ip == nvector) then
                  call feedbk(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)
                  ip = 0
                  ig = 0
               end if
               ipart = next_part  ! Go to next particle
            end do
            ! End loop over particles
         end if
         igrid = next(igrid)   ! Go to next grid
      end do
      ! End loop over grids
      if (ip > 0) call feedbk(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)

      if (MC_tracer) then
         call yield_tracers(icpu, ilevel, star_tracer_type)
      end if
   end do
   ! End loop over cpus

   if (MC_tracer) call post_particle_yield()

#endif

111 format('   Entering thermal_feedback for level ', I2)

end subroutine thermal_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbk(ind_grid, ind_part, ind_grid_part, ng, np, ilevel)
   use amr_commons
   use pm_commons
   use hydro_commons
   use random
   use imf_module
   use tracer_utils, only: mark_yielding_particle
   use constants, only: M_sun, Myr2sec, pc2cm, yr2sec, Mpc2cm, kpc2cm
   use metal_yields, only: AGB_Fe_yield, SNII_Fe_yield, OBwind_Fe_yield, SNIaFe, &
                           AGB_O_yield, SNII_O_yield, OBwind_O_yield, SNIaO, &
                           AGB_N_yield, SNII_N_yield, OBwind_N_yield, SNIaN, &
                           AGB_Mg_yield, SNII_Mg_yield, OBwind_Mg_yield, SNIaMg, &
                           AGB_Al_yield, SNII_Al_yield, OBwind_Al_yield, SNIaAl, &
                           AGB_Si_yield, SNII_Si_yield, OBwind_Si_yield, SNIaSi, &
                           AGB_Eu_yield, SNII_Eu_yield, OBwind_Eu_yield, SNIaEu, &
                           AGB_C_yield, SNII_C_yield, OBwind_C_yield, SNIaC, &
                           MEuNSNS, &
                           AGB_Ne_yield, SNII_Ne_yield, SNIaNe, &
                           AGB_S_yield, SNII_S_yield, SNIaS, &
                           AGB_Ca_yield, SNII_Ca_yield, SNIaCa, &
                           ! Harley
                           POP3_SNII_Fe_yield, POP3_HN_Fe_yield, POP3_HMHN_Fe_yield, &
                           POP3_SNII_O_yield, POP3_HN_O_yield, POP3_HMHN_O_yield, &
                           POP3_SNII_N_yield, POP3_HN_N_yield, POP3_HMHN_N_yield, &
                           POP3_SNII_Mg_yield, POP3_HN_Mg_yield, POP3_HMHN_Mg_yield, &
                           POP3_SNII_Al_yield, POP3_HN_Al_yield, POP3_HMHN_Al_yield, &
                           POP3_SNII_Si_yield, POP3_HN_Si_yield, POP3_HMHN_Si_yield, &
                           POP3_SNII_Eu_yield, POP3_HN_Eu_yield, POP3_HMHN_Eu_yield, &
                           POP3_SNII_C_yield, POP3_HN_C_yield, POP3_HMHN_C_yield, &
                           POP3_SNII_Ne_yield, POP3_HN_Ne_yield, POP3_HMHN_Ne_yield, &
                           POP3_SNII_S_yield, POP3_HN_S_yield, POP3_HMHN_S_yield, &
                           POP3_SNII_Ca_yield, POP3_HN_Ca_yield, POP3_HMHN_Ca_yield, &
                           POP3_SNII_H_yield, POP3_HN_H_yield, POP3_HMHN_H_yield, &
                           POP3_SNII_He_yield, POP3_HN_He_yield, POP3_HMHN_He_yield, &
                           AGB_H_yield, AGB_He_yield, SNII_H_yield, SNII_He_yield, &
                           table_Solar_Fe

   implicit none
   integer::ng, np, ilevel
   integer, dimension(1:nvector)::ind_grid
   integer, dimension(1:nvector)::ind_grid_part, ind_part
   !-----------------------------------------------------------------------
   ! This routine is called by subroutine feedback. Each stellar particle
   ! dumps mass, momentum and energy in the nearest grid cell using array
   ! unew.
   !-----------------------------------------------------------------------
   integer::i, j, idim, nx_loc, ilun, iii, ii, jj, kk
   real(kind=8)::RandNum
   real(dp)::dx_min, vol_min
   real(dp)::ESN, mejecta, time_simu, dx_loc
   real(dp)::dx, scale, birth_time
   real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
   ! Grid based arrays
   real(dp), dimension(1:nvector, 1:ndim), save::x0
   integer, dimension(1:nvector), save::ind_cell
   integer, dimension(1:nvector, 1:threetondim), save::nbors_father_cells
   integer, dimension(1:nvector, 1:twotondim), save::nbors_father_grids
   ! Particle based arrays
   logical, dimension(1:nvector), save::ok
   real(dp), dimension(1:nvector)::mloss, mzloss, ethermal, dteff
   real(dp), dimension(1:nvector)::mlossAGB
   real(dp), dimension(1:nvector, 1:nmetals)::mlossmetals
   real(dp), dimension(1:nvector)::ptot
   real(dp), dimension(1:nvector), save::vol_loc, dx_loc2
   real(dp), dimension(1:nvector, 1:ndim), save::x
   integer, dimension(1:nvector, 1:ndim), save::id, igd, icd
   integer, dimension(1:nvector), save::igrid, icell, indp, kg
   real(dp), dimension(1:3)::skip_loc
   real(dp)::Mremnant, numII, ENSN
   real(dp)::mstarmin, mstarmax, NumSNIa
   real(dp)::masslossIa, ESNIa
   real(dp)::masslossW
   real(dp)::theint, t1, t2
   real(dp)::Mwindmin, Mwindmax, mett
   real(dp)::vol_loc2
   real(dp):: NSNIa, MIMF, SNIIFe, SNIIO, SNIIej, NSNII, MIMFChabrier, IMFChabrier
   external NSNIa, MIMF, SNIIFe, SNIIO, SNIIej, NSNII, MIMFChabrier, IMFChabrier
   real(dp)::Zscale, Zgas
   integer::iicell, iskip
   integer, dimension(1:nvector, 1:2, 1:2, 1:2)::indcube2
   real(dp)::xcont, ycont, zcont, contr
   real(dp)::pII, pIa, pST0
   real(dp)::tt, tekin, vkick, tekinstar
   real(dp), dimension(1:nvector)::Prad
   real(dp)::vmax, momx, momy, momz, pST, n0
   real(dp)::meanmass
   real(dp)::SNmin, SNmax, Tmaxfb, rsf, cellsize, minmass
   real(dp)::maxadv, vxold, vyold, vzold, vxnew, vynew, vznew, Emax, yearscale
   integer::indd, imet
   real(dp)::scale_m, twind, numresidual, scale_mg, scale_msun
   real(dp)::SNyieldmcap, yieldZmin, yieldZmax, meanmassM, mettM, mettM_loc
   real(dp):: IMFKroupa, numAGB
   external IMFKroupa
   real(dp)::fNSNS_Ia, NumNSNS
   ! Harley
   logical::ok_pop3, star_is_pop3
   real(dp)::msun2g = 1.989d33
   real(dp)::t_life, esn_erg, mass, tyoung, current_time, vej, M_Hecore
   real(dp)::gyr2s = 3.1536000d+16
#ifdef NENER
   integer::irad
#endif
   real(dp)::SNIaFe_loc,SNIaO_loc,SNIaN_loc,SNIaMg_loc,SNIaAl_loc
   real(dp)::SNIaSi_loc,SNIaEu_loc,SNIaC_loc,SNIaCa_loc,SNIaNe_loc,SNIaS_loc

   integer::i_fractions, ifrac
   real(dp), dimension(1:NVAR), save::fractions
   real(dp)::tmp,C_ejecta,N_ejecta,O_ejecta,Ne_ejecta,Mg_ejecta,S_ejecta
   real(dp)::Si_ejecta,Ca_ejecta,Fe_ejecta,H_ejecta,He_ejecta
   logical::boundary_trigger
   ! MC tracer
   real(dp) :: star_original_mass
   real(dp)::total_esn_loc,nsnII_star_variable_resample,n_normal_SNII,n_hypernova
   real(dp)::input_sf,HNyieldmcap,total_mass_loss,mass_for_snIa

   ! Variable IMF
   real(dp)::tmp_fe_over_H,loc_upper_slope

   ! starting index for passive variables except for imetal and chem
   i_fractions = imetal + nchem + nmetals
   if (nco.gt.0) i_fractions = i_fractions + nco
   if (no_metal_update) i_fractions = imetal

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
   scale_m = M_sun/scale_d/scale_l/scale_l/scale_l !Msun into grams and then internal units
   scale_mg = scale_l**3*scale_d ! code units to grams
   scale_msun = scale_l**3*scale_d/msun2g !harley and taysun units

   ! Mesh spacing in that level
   dx = 0.5D0**ilevel
   nx_loc = (icoarse_max - icoarse_min + 1)
   skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
   if (ndim > 0) skip_loc(1) = dble(icoarse_min)
   if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
   if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
   scale = boxlen/dble(nx_loc)
   dx_loc = dx*scale
   dx_loc2(1:nvector) = dx*scale
   vol_loc(1:nvector) = dx_loc**ndim
   vol_loc2 = dx_loc**ndim
   dx_min = (0.5D0**nlevelmax)*scale
   vol_min = dx_min**ndim

   ! ---------------------------- SN feedback parameters ----------------------------
   ! Type Ia supernova parameters
   Mremnant = 1.4d0*scale_m   !Chandrasekhar mass, nothing remains
   pIa = Mremnant*8451d5/scale_v !0.5*Mremnant*v_ej^2 gives v_ej=8451 km/s for 1.4 Msun. Ejecta momentum p = Mremnant*v_ej

   ! Scale Ia yields to code units
   SNIaFe_loc = SNIaFe*scale_m
   SNIaO_loc  = SNIaO*scale_m
   SNIaN_loc  = SNIaN*scale_m
   SNIaMg_loc = SNIaMg*scale_m
   SNIaAl_loc = SNIaAl*scale_m
   SNIaSi_loc = SNIaSi*scale_m
   SNIaEu_loc = SNIaEu*scale_m
   SNIaC_loc  = SNIaC*scale_m
   SNIaS_loc  = SNIaS*scale_m
   SNIaNe_loc = SNIaNe*scale_m
   SNIaCa_loc = SNIaCa*scale_m

   ! AGB wind parameters
   Mwindmin = 0.5d0
   Mwindmax = 8.0d0

   ! 'fast wind' from OB stars
   twind = 1.0d7 !'fast winds'

   ! Massive star lifetime from Myr to code units
   yearscale = scale_t/yr2sec

   ! SNenergy=1.d51 !erg/SN
   ! Type II supernova specific energy from cgs to code units
   ESN = SNenergy/(10.*M_sun)/scale_v/scale_v  !energy per 10 Msun in internal units
   ENSN = SNenergy/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units
   pII = 12.0d0*3000.0d5*scale_m/scale_v !12 Msun at 3000 km/s per SNII event, calibrated to get same \dot{p} as SB99
   pST0 = 2.95d5*1.0d5*scale_m/scale_v   !---------- Into g cm/s and then into code units, Kim & Ostriker prefactor

   ESNIa = ENSN  !Oscar: assume same energy release as SNII
   SNmin = 8.0d0
   SNmax = 300.0d0   !Limits for SNII events. --> harley edit

   ! limiter
   vmax = vmaxFB
   Tmaxfb = Tmax
   maxadv = maxadvfb*1d5/scale_v !into internal units
   minmass = 0.1*scale_m  !minimum allowed particle mass

   ! Yield parameters.
   SNyieldmcap = 30.d0 ! Yields unknown above 30 Msun
   HNyieldmcap = 40.d0! Upper mass for hypernova
   yieldZmin = 0.0001d0 ! Yields unknown below Z=0.0001
   yieldZmax = 0.02d0   ! Yields unknown above Z=0.02
   if (yields_portinari) then
      SNyieldmcap = 100.0 ! Yields unknown above 100 Msun
      yieldZmin = 0.0004 ! Yields unknown below Z=0.0004
      yieldZmax = 0.05   ! Yields unknown above Z=0.05
   end if
   if (yields_lc18) then
      SNyieldmcap = 120.0 ! Yields unknown above 120 Msun
      yieldZmin = (1.d-3)*(0.02d0) ! Yields unknown below 10^-3 solar
      yieldZmax = (1.d0)*(0.2d0)   ! Yields unknown above 1 solar
   end if

   ! Neutron star mergers, for r-process
   ! Using same model as Naiman et al. (2018) but renomalised to match updated NSNS rates.
   fNSNS_Ia = 4.6d-2 ! Comparing Ia rates (Maoz & Graur, 2017) to NSNS rates (Abbott et al., 2017)
   MEuNSNS = MEuNSNS*scale_m
   !-------------------------------------------------------------------------------------

   ! Lower left corner of 3x3x3 grid-cube
   do idim = 1, ndim
      do i = 1, ng
         x0(i, idim) = xg(ind_grid(i), idim) - 3.0D0*dx
      end do
   end do

   ! Gather 27 neighboring father cells (should be present anytime !)
   do i = 1, ng
      ind_cell(i) = father(ind_grid(i))
   end do
   call get3cubefather(ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)

   ! Rescale position at level ilevel
   do idim = 1, ndim
      do j = 1, np
         x(j, idim) = xp(ind_part(j), idim)/scale + skip_loc(idim)
      end do
   end do
   do idim = 1, ndim
      do j = 1, np
         x(j, idim) = x(j, idim) - x0(ind_grid_part(j), idim)
      end do
   end do
   do idim = 1, ndim
      do j = 1, np
         x(j, idim) = x(j, idim)/dx
      end do
   end do

   ! NGP at level ilevel
   do idim = 1, ndim
      do j = 1, np
         id(j, idim) = int(x(j, idim))
      end do
   end do

   ! Compute parent grids
   do idim = 1, ndim
      do j = 1, np
         igd(j, idim) = id(j, idim)/2
      end do
   end do
   do j = 1, np
      kg(j) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
   end do
   do j = 1, np
      igrid(j) = son(nbors_father_cells(ind_grid_part(j), kg(j)))
   end do

   ! Check if particles are entirely in level ilevel
   ok(1:np) = .true.
   do j = 1, np
      ok(j) = ok(j) .and. igrid(j) > 0
   end do

   ! Compute parent cell position
   do idim = 1, ndim
      do j = 1, np
         if (ok(j)) then
            icd(j, idim) = id(j, idim) - 2*igd(j, idim)
         end if
      end do
   end do
   do j = 1, np
      if (ok(j)) then
         icell(j) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
      end if
   end do

   ! Compute parent cell adresses
   do j = 1, np
      if (ok(j)) then
         indp(j) = ncoarse + (icell(j) - 1)*ngridmax + igrid(j)
      else
         indp(j) = nbors_father_cells(ind_grid_part(j), kg(j))
         vol_loc(j) = vol_loc(j)*2**ndim !ilevel-1 cell volume
         dx_loc2(j) = dx_loc2(j)*2.0
      end if
   end do

   ! Loop over cells
   indd = 0
   do kk = 1, 2
      do jj = 1, 2
         do ii = 1, 2
            indd = indd + 1  !counter from 1,twotondim
            iskip = ncoarse + (indd - 1)*ngridmax
            do j = 1, np
               ind_cell(j) = iskip + igrid(j)
               indcube2(j, ii, jj, kk) = ind_cell(j)
            end do
         end do
      end do
   end do

   ! Compute individual time steps
   do j = 1, np
      dteff(j) = dtnew(levelp(ind_part(j)))
   end do

   if (use_proper_time) then
      do j = 1, np
         dteff(j) = dteff(j)*aexp**2
      end do
   end if

   ! Reset ejected mass, metallicity, thermal energy
   mloss(:) = 0.d0
   mlossAGB(:) = 0.d0
   mzloss(:) = 0.d0
   ethermal(:) = 0.d0
   Prad(:) = 0.0d0
   ptot(:) = 0.0d0
   mlossmetals(:, :) = 0.0

   if (cosmo) then
      ! Find neighboring expansion factors
      i = 1
      do while (aexp_frw(i) > aexp .and. i < n_frw)
         i = i + 1
      end do
      ! Interploate time
      time_simu = t_frw(i)*(aexp - aexp_frw(i - 1))/(aexp_frw(i) - aexp_frw(i - 1)) + &
           & t_frw(i - 1)*(aexp - aexp_frw(i))/(aexp_frw(i - 1) - aexp_frw(i))
   end if

   ! For Pop III SN
   if (use_proper_time) then
      current_time = texp
   else
      current_time = t
   end if

   ! Compute feedback
   do j = 1, np              !----------------------- Begin loop over all particles
      mejecta = 0.0
      numII = 0.0
      NumSNIa = 0.0

      ! For MC tracers, we need to know how much mass the star
      ! has lost between now and the end of the feedback, so we store
      ! the mass it has now, and will compare to how much it has at the end
      ! of this long do ... loop
      star_original_mass = mp(ind_part(j))

      if (cosmo) then
         if (use_proper_time) then    !Correct as tp is in proper, and dteff is proper. Now get to yrs
            call getAgeGyr(tp(ind_part(j)), t1)          !  End-of-dt age [Gyrs]
            call getAgeGyr(tp(ind_part(j)) - dteff(j), t2) !  End-of-dt age [Gyrs]
            t1 = t1*1.0d9
            t2 = t2*1.0d9
         else
            ! Compute star age in years,conformal time to years
            iii = 1
            do while (tau_frw(iii) > tp(ind_part(j)) .and. iii < n_frw)
               iii = iii + 1
            end do

            t2 = t_frw(iii)*(tp(ind_part(j)) - tau_frw(iii - 1))/(tau_frw(iii) - tau_frw(iii - 1)) + &
                 & t_frw(iii - 1)*(tp(ind_part(j)) - tau_frw(iii))/(tau_frw(iii - 1) - tau_frw(iii))
            t2 = (time_simu - t2)/(h0*1.d5/Mpc2cm)/(yr2sec)                        !Units of years; particle age

            t1 = t_frw(iii)*(tp(ind_part(j)) + dteff(j) - tau_frw(iii - 1))/(tau_frw(iii) - tau_frw(iii - 1)) + &
                 & t_frw(iii - 1)*(tp(ind_part(j)) + dteff(j) - tau_frw(iii))/(tau_frw(iii - 1) - tau_frw(iii))
            t1 = (time_simu - t1)/(h0*1.d5/Mpc2cm)/(yr2sec)                        !Units of years; particle age - dt

            birth_time = tp(ind_part(j))
            t1 = max(t1, 0.0d0)
            t2 = max(t2, 0.0d0)
         end if
      else
         t2 = t - tp(ind_part(j))           !Age at t
         t2 = t2*yearscale !For non-cosmo units, in years
         t1 = t - tp(ind_part(j)) - dteff(j)  !Age at t-dt, in years
         t1 = t1*yearscale
         birth_time = tp(ind_part(j))     !internal units
         t1 = max(t1, 0.0d0)
      end if

      !--------------------------- Get dying stars ----------------------
      !  ---------- metals is star and host cell
      ! imetal=Fe, imetal+1=O, this must hold strictly for Zsolar to be correct elsewhere

      mett = 2.09d0*zp(ind_part(j), 2) + 1.06d0*zp(ind_part(j), 1) !---- Stellar metallicity, Asplund 2009
      Zgas = (2.09d0*unew(indp(j), imetal + 1) + 1.06d0*unew(indp(j), imetal))/ &
           & max(unew(indp(j), 1), smallr)/0.02d0              !---- Average, solar mix in gas (Asplund)
      Zgas = max(Zgas, 0.0)

      mettM = min(mett, yieldZmax)                         !------ yields assumed to be the same as limits
      mettM = max(mett, yieldZmin)                         !------ if outside Z range.

      ok_pop3 = .false.
      star_is_pop3 = is_pop_III(typep(ind_part(j)))
      
#ifdef POP3
      if (star_is_pop3 .and. is_pre_SN(typep(ind_part(j)))) then
         call get_pop3_ageGyr_fb(mp(ind_part(j)), t_life)
         if (use_proper_time) then
            tyoung = t_life*gyr2s/(scale_t/aexp**2)
         else
            tyoung = t_life*gyr2s/scale_t
         end if
         tyoung = current_time - tyoung

         ! if tp is older than t_delay
         if (tp(ind_part(j)) .le. tyoung) then
            ok_pop3 = .true.
            if (supernovae) then
               typep(ind_part(j))%tag = STAR_POP3_AFTER_SNe
            end if

            ! if inert, collisionless BH particle, do nothing and only change the switch
            mass = mp(ind_part(j))*scale_msun
            if ((mass .gt. 40.d0 .and. mass .lt. 140.d0) .or. (mass .gt. 300.d0)) then
               ok_pop3 = .false.
            end if
         end if
      else
#endif
         if (t1 .gt. 0.0) then
            ! Harley update for the new bpass model
            !call agemass_bpassv2_300_interp(t1, mstarmax)  !Get stellar masses that exit main sequence during dt
            !call agemass_bpassv2_300_interp(t2, mstarmin)
            call agemass_schaerer_1993(t1, mstarmax)  !Get stellar masses that exit main sequence during dt
            call agemass_schaerer_1993(t2, mstarmin)  !Get stellar masses that exit main sequence during dt
         else
            mstarmin = SNmax
            mstarmax = SNmax
         end if
         mass_for_snIa = mstarmin

         ! Enforce bounds
         if (mstarmin .le. SNmax .and. mstarmax .ge. SNmax) then
            mstarmax = SNmax - 1d-2
         end if
         if (mstarmin .le. SNmin .and. mstarmax .ge. SNmin) then
            mstarmin = SNmin + 1d-2
         end if
         if (mstarmin .le. SNmin .and. mstarmax .ge. SNmax) then !if timestep is >40 Myr
            mstarmin = SNmin + 1d-2
            mstarmax = SNmax - 1d-2
         end if
         meanmass = (mstarmax + mstarmin)/2.0d0

#ifdef POP3
      end if
#endif

      !---------------------------
      !--------------------------- Supernovae
      !---------------------------
      if (supernovae) then
         !---------------------------
         !--------------------------- Type II events
         if (mstarmax .le. SNmax .and. mstarmin .ge. SNmin .or. ok_pop3) then
            numII = 0.d0
            ! FOR POP III STARS
            if (star_is_pop3.and.ok_pop3) then
               mass = mp(ind_part(j))*scale_msun
               ! Calculate the number of SN for Pop III stars
               if (mass .ge. 11.d0 .and. mass .le. 20.d0) then
                  numII = 1.0
               else if (mass .gt. 20.d0 .and. mass .le. 40.d0) then
                  numII = 1.0
               else if (mass .ge. 140.d0 .and. mass .le. 300.d0) then
                  numII = 1.0
               else
                  numII = 0.0
               end if
            ! FOR POP II STARS
            else if (mstarmax .le. SNmax .and. mstarmin .ge. SNmin .and. .not.star_is_pop3.and.mp(ind_part(j)).gt.minmass) then
               !call SNIInum_harley_of(mp0(ind_part(j))*scale_msun, mstarmin, mstarmax, theint)
               loc_upper_slope = var_imf_a2
               tmp_fe_over_H = LOG10(MAX((zp(ind_part(j), 1) / (0.76d0))/table_Solar_Fe,1.d-40)) ! [Fe/H]
               
               ! Recalculate the upper mass IMF slope if needed 
               if (imf_varies_with_metallicity) call get_upper_slope_prgomet(tmp_fe_over_H,loc_upper_slope)
               if (imf_varies_complex) then
                  loc_upper_slope = -1.d0 * part_imf_slope(ind_part(j)) 
               end if

               call get_N_supernova_in_step(var_imf_m0,var_imf_m1,var_imf_m2,var_imf_a1,loc_upper_slope,mp0(ind_part(j))*scale_msun, mstarmin, mstarmax, theint)
               numII = theint
               call ranf(localseed, RandNum)
               numresidual = numII - int(numII)  !---- int <1 --> 0
               numII = int(numII)
               if (RandNum < numresidual) then
                  numII = numII + 1              !---- Another star from residual --> HK note (oscar is clever!)
               end if 
            end if

            if (numII .gt. 0.d0) then
               total_esn_loc = 0.d0 ! Reset to zero
               mejecta = 0.d0       ! Reset to zero
               n0 = max(unew(indp(j), 1), smallr)*scale_nH              !---------- mH/cc
               Zscale = (max(Zgas, 0.01d0))**(-0.2)                    !---------- scaling from Thornton et al. 1998.
               if (.not. ok_pop3) then

                  ! HKNOTE --> this is where we want to put the HN
                  call get_effective_num_sn(mett,meanmass,numII,nsnII_star_variable_resample,n_normal_SNII,n_hypernova)
397               format('ALL SN -- a:',F10.6,' Mass (Msun):',F10.3,' Total energy (10^51 erg):',F10.3,' # normal SN:',I5,' # HN:',I5,' Density (log10 cm^-3):',F10.3,' Star Z (abs):',F10.3,' Gas Z (sol):',F10.3)
                  write(*,397) aexp,meanmass,nsnII_star_variable_resample*(SNenergy/1.d51),NINT(n_normal_SNII),NINT(n_hypernova),LOG10(n0),LOG10(mett),LOG10(MAX(1.d-40,Zgas))

                  ESN = (SNenergy*nsnII_star_variable_resample)/(10.*M_sun)/scale_v/scale_v  !energy per 10 Msun in internal units
                  ENSN = (SNenergy*nsnII_star_variable_resample)/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units

                  ! Harley is updating this below
                  !pII = 12.0d0*3000.0d5*scale_m/scale_v !12 Msun at 3000 km/s per SNII event, calibrated to get same \dot{p} as SB99
                  pII = 0.d0 ! --> Harley changed this as it's updated below so that momentum is calculated from the actual ejecta and SN energy
                  pST0 = 2.95d5*1.0d5*scale_m/scale_v   !---------- Into g cm/s and then into code units, Kim & Ostriker prefactor

                  ! --> Harley checked this is ok!
                  pST = pST0*((nsnII_star_variable_resample*(SNenergy/1.d51))**0.941d0)*n0**(-0.1176d0)*Zscale         !---------- km/sST momentum Blondin et al (1998).

                  total_esn_loc = total_esn_loc + (SNenergy * nsnII_star_variable_resample)         !---------- total energy released in erg

                  ! HARLEY CHANGE --> DO EJECTA HERE
                  ! Yields for hypernova
                  if (variable_energy_SN.and.n_hypernova.gt.0.d0) then
                     Fe_ejecta = (n_hypernova*SNII_Fe_yield( min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !1=Fe metal enrichment
                     O_ejecta  = (n_hypernova*SNII_O_yield(  min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !2=O  metal enrichment
                     N_ejecta  = (n_hypernova*SNII_N_yield(  min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !3=N  metal enrichment
                     Mg_ejecta = (n_hypernova*SNII_Mg_yield( min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !4=Mg metal enrichment
                     Ne_ejecta = (n_hypernova*SNII_Ne_yield( min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !5=Ne metal enrichment
                     Si_ejecta = (n_hypernova*SNII_Si_yield( min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !6=Si metal enrichment
                     Ca_ejecta = (n_hypernova*SNII_Ca_yield( min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !7=Ca metal enrichment
                     C_ejecta  = (n_hypernova*SNII_C_yield(  min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !8=C  metal enrichment
                     S_ejecta  = (n_hypernova*SNII_S_yield(  min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !9=S  metal enrichment
                     H_ejecta  = (n_hypernova*SNII_H_yield(  min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !H enrichment
                     He_ejecta = (n_hypernova*SNII_He_yield( min(meanmass, HNyieldmcap), MAX(MIN(mett,0.02d0),0.0001d0), .false., .false., .true.)/scale_msun)  !He enrichment

                     ! Sum contribution to Mejecta
                     total_mass_loss = 0.d0
                     total_mass_loss = total_mass_loss + Fe_ejecta + O_ejecta + N_ejecta + Mg_ejecta + Ne_ejecta + Si_ejecta
                     total_mass_loss = total_mass_loss + Ca_ejecta + C_ejecta + S_ejecta + H_ejecta  + He_ejecta

                     input_sf = 1.d0
                     if (total_mass_loss.gt.((0.999d0*n_hypernova*meanmass/scale_msun))) then
                        input_sf = total_mass_loss/(0.999d0*meanmass*n_hypernova*scale_msun)
                        total_mass_loss = total_mass_loss/input_sf
                     endif

                     ! Harley from Kim an Ostriker 2015
                     ! The velocity of the ballistic ejecta is (2*ESN/Mej) ** 0.5

                     ! If less than than 25, need to subtract off normal SN energy
                     if (meanmass.gt.25.d0) then
                        vej = ((2.0*nsnII_star_variable_resample*SNenergy/n_hypernova)/(total_mass_loss*scale_mg/n_hypernova))**0.5 ! in units of cm / s
                     else
                        vej = ((2.0*(nsnII_star_variable_resample-n_normal_SNII)*SNenergy/n_hypernova)/(total_mass_loss*scale_mg/n_hypernova))**0.5 ! in units of cm / s
                     end if
                     pII = pII + (total_mass_loss*vej/scale_v) ! initial blast wave momentum in code units

                     mejecta = mejecta + total_mass_loss 
                     if (metal) then
                        ! Update the mloss vector
                        mlossmetals(j, 1) = mlossmetals(j, 1) + (Fe_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 2) = mlossmetals(j, 2) + (O_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 3) = mlossmetals(j, 3) + (N_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 4) = mlossmetals(j, 4) + (Mg_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 5) = mlossmetals(j, 5) + (Ne_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 6) = mlossmetals(j, 6) + (Si_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 7) = mlossmetals(j, 7) + (Ca_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 8) = mlossmetals(j, 8) + (C_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 9) = mlossmetals(j, 9) + (S_ejecta /vol_loc(j))/input_sf
                     end if
398                  format('HN -- a:',F10.6,' Age range (Myr)',F10.3,'-',F10.3,' Mass (Msun):',F10.3,' Ejected mass (Msun):',F10.3,' Fe=',F10.3,' O=',F10.3,' N=',F10.3,' Mg=',F10.3,' Ne=',F10.3,' Si=',F10.3,' Ca=',F10.3,' C=',F10.3,' S=',F10.3)
                     write(*,398) aexp,t1/1.d6,t2/1.d6,meanmass, scale_msun*total_mass_loss/n_hypernova, &
                                  scale_msun*Fe_ejecta/n_hypernova, &
                                  scale_msun*O_ejecta/n_hypernova, &
                                  scale_msun*N_ejecta/n_hypernova, &
                                  scale_msun*Mg_ejecta/n_hypernova, &
                                  scale_msun*Ne_ejecta/n_hypernova, &
                                  scale_msun*Si_ejecta/n_hypernova, &
                                  scale_msun*Ca_ejecta/n_hypernova, &
                                  scale_msun*C_ejecta/n_hypernova, &
                                  scale_msun*S_ejecta/n_hypernova
                  end if
                  ! Yields for normal supernova
                  if (n_normal_SNII.gt.0.d0) then
                     input_sf = 1.d0
                     if (yields_lc18) input_sf = 0.02
                     
                     Fe_ejecta = (n_normal_SNII*SNII_Fe_yield( min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !1=Fe metal enrichment
                     O_ejecta  = (n_normal_SNII*SNII_O_yield(  min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !2=O  metal enrichment
                     N_ejecta  = (n_normal_SNII*SNII_N_yield(  min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !3=N  metal enrichment
                     Mg_ejecta = (n_normal_SNII*SNII_Mg_yield( min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !4=Mg metal enrichment
                     Ne_ejecta = (n_normal_SNII*SNII_Ne_yield( min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !5=Ne metal enrichment
                     Si_ejecta = (n_normal_SNII*SNII_Si_yield( min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !6=Si metal enrichment
                     Ca_ejecta = (n_normal_SNII*SNII_Ca_yield( min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !7=Ca metal enrichment
                     C_ejecta  = (n_normal_SNII*SNII_C_yield(  min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !8=C  metal enrichment
                     S_ejecta  = (n_normal_SNII*SNII_S_yield(  min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !9=S  metal enrichment
                     H_ejecta  = (n_normal_SNII*SNII_H_yield(  min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !H enrichment
                     He_ejecta = (n_normal_SNII*SNII_He_yield( min(meanmass, SNyieldmcap), mettM/input_sf, yields_portinari, yields_lc18)/scale_msun)  !He enrichment 

                     ! Sum contribution to Mejecta
                     total_mass_loss = 0.d0
                     total_mass_loss = total_mass_loss + Fe_ejecta + O_ejecta + N_ejecta + Mg_ejecta + Ne_ejecta + Si_ejecta
                     total_mass_loss = total_mass_loss + Ca_ejecta + C_ejecta + S_ejecta + H_ejecta  + He_ejecta

                     input_sf = 1.d0
                     if (total_mass_loss.gt.((0.999d0*n_normal_SNII*meanmass/scale_msun))) then
                        input_sf = total_mass_loss/(0.999d0*meanmass*n_normal_SNII*scale_msun)
                        total_mass_loss = total_mass_loss/input_sf
                     endif

                     ! If less than than 25, need to subtract off normal SN energy
                     if (meanmass.gt.25.d0) then
                        vej = ((2.0*SNenergy*1.d-3)/(total_mass_loss*scale_mg/n_normal_SNII))**0.5 ! in units of cm / s
                     else
                        vej = ((2.0*SNenergy)/(total_mass_loss*scale_mg/n_normal_SNII))**0.5 ! in units of cm / s
                     end if
                     pII = pII + (total_mass_loss*vej/scale_v) ! initial blast wave momentum in code units

                     ! Make sure that is there is no supernova we always add the correct momentum
                     if (nsnII_star_variable_resample.lt.1.d0) then
                        !SN will not be resolved but still add ejecta momentum
                        pST = pII
                     end if

                     mejecta = mejecta + total_mass_loss
                     if (metal) then
                        ! Update the mloss vector
                        mlossmetals(j, 1) = mlossmetals(j, 1) + (Fe_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 2) = mlossmetals(j, 2) + (O_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 3) = mlossmetals(j, 3) + (N_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 4) = mlossmetals(j, 4) + (Mg_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 5) = mlossmetals(j, 5) + (Ne_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 6) = mlossmetals(j, 6) + (Si_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 7) = mlossmetals(j, 7) + (Ca_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 8) = mlossmetals(j, 8) + (C_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 9) = mlossmetals(j, 9) + (S_ejecta /vol_loc(j))/input_sf
                     end if
399                  format('SN -- a:',F10.6,' Age range (Myr)',F10.3,'-',F10.3,' Mass (Msun):',F10.3,' Ejected mass (Msun):',F10.3,' Fe=',F6.3,' O=',F6.3,' N=',F6.3,' Mg=',F6.3,' Ne=',F6.3,' Si=',F6.3,' Ca=',F6.3,' C=',F6.3,' S=',F6.3)
                     write(*,399) aexp,t1/1.d6,t2/1.d6,meanmass, scale_msun*total_mass_loss/n_normal_SNII, & 
                                  scale_msun*Fe_ejecta/n_normal_SNII, &
                                  scale_msun*O_ejecta/n_normal_SNII, &
                                  scale_msun*N_ejecta/n_normal_SNII, &
                                  scale_msun*Mg_ejecta/n_normal_SNII, &
                                  scale_msun*Ne_ejecta/n_normal_SNII, &
                                  scale_msun*Si_ejecta/n_normal_SNII, &
                                  scale_msun*Ca_ejecta/n_normal_SNII, &
                                  scale_msun*C_ejecta/n_normal_SNII, &
                                  scale_msun*S_ejecta/n_normal_SNII
                  end if
               else
                  mass = mp0(ind_part(j))*scale_msun
                  if (mass .ge. 11.d0 .and. mass .le. 20.d0) then
                     ! Nomoto + (2006) - normal SN
                     esn_erg = 1d51
                     
                     Fe_ejecta   = POP3_SNII_Fe_yield(mass)/scale_msun !1=Fe Harley
                     O_ejecta    = POP3_SNII_O_yield(mass) /scale_msun !2=O  Harley
                     N_ejecta    = POP3_SNII_N_yield(mass) /scale_msun !3=N  Harley
                     Mg_ejecta   = POP3_SNII_Mg_yield(mass)/scale_msun !4=Mg Harley
                     Ne_ejecta   = POP3_SNII_Ne_yield(mass)/scale_msun !5=Ne Harley
                     Si_ejecta   = POP3_SNII_Si_yield(mass)/scale_msun !6=Si Harley
                     Ca_ejecta   = POP3_SNII_Ca_yield(mass)/scale_msun !7=Ca Harley
                     C_ejecta    = POP3_SNII_C_yield(mass) /scale_msun !8=C  Harley
                     S_ejecta    = POP3_SNII_S_yield(mass) /scale_msun !9=S  Harley
                     H_ejecta    = POP3_SNII_H_yield(mass) /scale_msun !H ejecta  Harley
                     He_ejecta   = POP3_SNII_He_yield(mass)/scale_msun !He ejecta  Harley

                     ! Sum contribution to Mejecta
                     total_mass_loss = 0.d0
                     total_mass_loss = total_mass_loss + Fe_ejecta + O_ejecta + N_ejecta + Mg_ejecta + Ne_ejecta + Si_ejecta
                     total_mass_loss = total_mass_loss + Ca_ejecta + C_ejecta + S_ejecta + H_ejecta  + He_ejecta

                     input_sf = 1.d0
                     if (total_mass_loss.gt.(0.999d0*mass/scale_msun)) then
                        input_sf = total_mass_loss/(0.999d0*mass/scale_msun)
                        total_mass_loss = total_mass_loss/input_sf
                     endif
                     mejecta = mejecta + total_mass_loss 
                     if (metal) then
                        ! Update the mloss vector
                        mlossmetals(j, 1) = mlossmetals(j, 1) + (Fe_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 2) = mlossmetals(j, 2) + (O_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 3) = mlossmetals(j, 3) + (N_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 4) = mlossmetals(j, 4) + (Mg_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 5) = mlossmetals(j, 5) + (Ne_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 6) = mlossmetals(j, 6) + (Si_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 7) = mlossmetals(j, 7) + (Ca_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 8) = mlossmetals(j, 8) + (C_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 9) = mlossmetals(j, 9) + (S_ejecta /vol_loc(j))/input_sf
                     end if
404                  format('Pop3 SN -- a:',F10.6,' Age (Myr)',F10.3,' Mass (Msun):',F10.3,' Ejected mass (Msun):',F10.3,' Fe=',F10.3,' O=',F10.3,' N=',F10.3,' Mg=',F10.3,' Ne=',F10.3,' Si=',F10.3,' Ca=',F10.3,' C=',F10.3,' S=',F10.3)
                     write(*,404) aexp,t1/1.d6,mass,scale_msun*total_mass_loss, &
                                  scale_msun*Fe_ejecta, &
                                  scale_msun*O_ejecta, &
                                  scale_msun*N_ejecta, &
                                  scale_msun*Mg_ejecta, &
                                  scale_msun*Ne_ejecta, &
                                  scale_msun*Si_ejecta, &
                                  scale_msun*Ca_ejecta, &
                                  scale_msun*C_ejecta, &
                                  scale_msun*S_ejecta
                  else if (mass .gt. 20.d0 .and. mass .le. 40.d0) then
                     ! Nomoto + (2006) - Hypernova
                     esn_erg = 1d51*(-13.7143 + 1.08571*mass)
                     
                     Fe_ejecta   = POP3_HN_Fe_yield(mass)/scale_msun !1=Fe Harley
                     O_ejecta    = POP3_HN_O_yield(mass) /scale_msun !2=O  Harley
                     N_ejecta    = POP3_HN_N_yield(mass) /scale_msun !3=N  Harley
                     Mg_ejecta   = POP3_HN_Mg_yield(mass)/scale_msun !4=Mg Harley
                     Ne_ejecta   = POP3_HN_Ne_yield(mass)/scale_msun !5=Ne Harley
                     Si_ejecta   = POP3_HN_Si_yield(mass)/scale_msun !6=Si Harley
                     Ca_ejecta   = POP3_HN_Ca_yield(mass)/scale_msun !7=Ca Harley
                     C_ejecta    = POP3_HN_C_yield(mass) /scale_msun !8=C  Harley
                     S_ejecta    = POP3_HN_S_yield(mass) /scale_msun !9=S  Harley
                     H_ejecta    = POP3_HN_H_yield(mass) /scale_msun !H ejecta  Harley
                     He_ejecta   = POP3_HN_He_yield(mass)/scale_msun !He ejecta  Harley

                     ! Sum contribution to Mejecta
                     total_mass_loss = 0.d0
                     total_mass_loss = total_mass_loss + Fe_ejecta + O_ejecta + N_ejecta + Mg_ejecta + Ne_ejecta + Si_ejecta
                     total_mass_loss = total_mass_loss + Ca_ejecta + C_ejecta + S_ejecta + H_ejecta  + He_ejecta

                     input_sf = 1.d0
                     if (total_mass_loss.gt.(0.999d0*mass/scale_msun)) then
                        input_sf = total_mass_loss/(0.999d0*mass/scale_msun)
                        total_mass_loss = total_mass_loss/input_sf
                     endif
                     mejecta = mejecta + total_mass_loss
                    if (metal) then
                        ! Update the mloss vector
                        mlossmetals(j, 1) = mlossmetals(j, 1) + (Fe_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 2) = mlossmetals(j, 2) + (O_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 3) = mlossmetals(j, 3) + (N_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 4) = mlossmetals(j, 4) + (Mg_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 5) = mlossmetals(j, 5) + (Ne_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 6) = mlossmetals(j, 6) + (Si_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 7) = mlossmetals(j, 7) + (Ca_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 8) = mlossmetals(j, 8) + (C_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 9) = mlossmetals(j, 9) + (S_ejecta /vol_loc(j))/input_sf
                     end if
405                  format('Pop3 HN -- a:',F10.6,' Age (Myr)',F10.3,' Mass (Msun):',F10.3,' Ejected mass (Msun):',F10.3,' Fe=',F10.3,' O=',F10.3,' N=',F10.3,' Mg=',F10.3,' Ne=',F10.3,' Si=',F10.3,' Ca=',F10.3,' C=',F10.3,' S=',F10.3)
                     write(*,405) aexp,t1/1.d6,mass,scale_msun*total_mass_loss, &
                                  scale_msun*Fe_ejecta, &
                                  scale_msun*O_ejecta, &
                                  scale_msun*N_ejecta, &
                                  scale_msun*Mg_ejecta, &
                                  scale_msun*Ne_ejecta, &
                                  scale_msun*Si_ejecta, &
                                  scale_msun*Ca_ejecta, &
                                  scale_msun*C_ejecta, &
                                  scale_msun*S_ejecta
                  else if (mass .ge. 140.d0 .and. mass .le. 300.d0) then
                     ! Heger & Woosley (2002)
                     M_Hecore = (13.d0/24.d0)*(mass - 20.d0)
                     esn_erg = 1d51*(5.d0 + 1.304d0*(M_Hecore - 64.d0))
                  
                     Fe_ejecta   = POP3_HMHN_Fe_yield(mass)/scale_msun !1=Fe Harley
                     O_ejecta    = POP3_HMHN_O_yield(mass) /scale_msun !2=O  Harley
                     N_ejecta    = POP3_HMHN_N_yield(mass) /scale_msun !3=N  Harley
                     Mg_ejecta   = POP3_HMHN_Mg_yield(mass)/scale_msun !4=Mg Harley
                     Ne_ejecta   = POP3_HMHN_Ne_yield(mass)/scale_msun !5=Ne Harley
                     Si_ejecta   = POP3_HMHN_Si_yield(mass)/scale_msun !6=Si Harley
                     Ca_ejecta   = POP3_HMHN_Ca_yield(mass)/scale_msun !7=Ca Harley
                     C_ejecta    = POP3_HMHN_C_yield(mass) /scale_msun !8=C  Harley
                     S_ejecta    = POP3_HMHN_S_yield(mass) /scale_msun !9=S  Harley
                     H_ejecta    = POP3_HMHN_H_yield(mass) /scale_msun !H ejecta  Harley
                     He_ejecta   = POP3_HMHN_He_yield(mass)/scale_msun !He ejecta  Harley

                     ! Sum contribution to Mejecta
                     total_mass_loss = 0.d0
                     total_mass_loss = total_mass_loss + Fe_ejecta + O_ejecta + N_ejecta + Mg_ejecta + Ne_ejecta + Si_ejecta
                     total_mass_loss = total_mass_loss + Ca_ejecta + C_ejecta + S_ejecta + H_ejecta  + He_ejecta

                     input_sf = 1.d0
                     if (total_mass_loss.gt.(0.999d0*mass/scale_msun)) then
                        input_sf = total_mass_loss/(0.999d0*mass/scale_msun)
                        total_mass_loss = total_mass_loss/input_sf
                     endif
                     mejecta = mejecta + total_mass_loss 
                     if (metal) then
                        ! Update the mloss vector
                        mlossmetals(j, 1) = mlossmetals(j, 1) + (Fe_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 2) = mlossmetals(j, 2) + (O_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 3) = mlossmetals(j, 3) + (N_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 4) = mlossmetals(j, 4) + (Mg_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 5) = mlossmetals(j, 5) + (Ne_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 6) = mlossmetals(j, 6) + (Si_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 7) = mlossmetals(j, 7) + (Ca_ejecta/vol_loc(j))/input_sf
                        mlossmetals(j, 8) = mlossmetals(j, 8) + (C_ejecta /vol_loc(j))/input_sf
                        mlossmetals(j, 9) = mlossmetals(j, 9) + (S_ejecta /vol_loc(j))/input_sf
                     end if
406                  format('Pop3 HMHN -- a:',F10.6,' Age (Myr)',F10.3,' Mass (Msun):',F10.3,' Ejected mass (Msun):',F10.3,' Fe=',F10.3,' O=',F10.3,' N=',F10.3,' Mg=',F10.3,' Ne=',F10.3,' Si=',F10.3,' Ca=',F10.3,' C=',F10.3,' S=',F10.3)
                     write(*,406) aexp,t1/1.d6,mass,scale_msun*total_mass_loss, &
                                  scale_msun*Fe_ejecta, &
                                  scale_msun*O_ejecta, &
                                  scale_msun*N_ejecta, &
                                  scale_msun*Mg_ejecta, &
                                  scale_msun*Ne_ejecta, &
                                  scale_msun*Si_ejecta, &
                                  scale_msun*Ca_ejecta, &
                                  scale_msun*C_ejecta, &
                                  scale_msun*S_ejecta
                  else
                     esn_erg = 0d0
                     mejecta = 0d0
                  end if

                  ! Adopted relation for momentum from Oscaar's paper: https://arxiv.org/pdf/2006.06008.pdf

                  ! calculate the sedov taylor momentum injection for a pop III star
                  pST = 2.95d5*((numII*esn_erg/1d51)**(0.941d0))*(n0**(-0.1176d0))*Zscale*1.0d5*scale_m/scale_v ! momentum in code units

                  ENSN = numII*esn_erg/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units

                  ! Harley from Kim an Ostriker 2015
                  ! The velocity of the ballistic ejecta is (2*ESN/Mej) ** 0.5
                  vej = ((2.0d0*esn_erg)/(mejecta*scale_mg))**0.5d0 ! in units of cm / s
                  pII = mejecta*vej/scale_v ! initial blast wave momentum in code units for pop III

                  total_esn_loc = total_esn_loc + (numII * esn_erg)
               end if

               !---------- Check cooling radius
               !rsf = 30.d0*numII**0.29d0*(n0**(-0.43d0))*((Zgas + 0.01d0)**(-0.18d0)) !------ Shell formation radius (Hopkins et al. 2013, Coiffi, et al.)

               ! Harley modified this by removing numII and replacing with total e_sn
               rsf = 30.d0*(n0**(-0.43d0))*((Zgas + 0.01d0)**(-0.18d0)) !------ Shell formation radius (Hopkins et al. 2013, Coiffi, et al.)
               rsf = rsf * ((total_esn_loc/1d51)**0.2857d0) ! need to scale by SN energy with https://arxiv.org/pdf/1702.06148.pdf

               cellsize = dx_loc2(j)*scale_l/pc2cm                      !---------- local resolution element

               if (cellsize .lt. rsf/Nrcool) then                       !----------- Kim & Ostriker criterion. Momentum: 1/3, energy: 1/10.
                  momST = .false.
               else
                  momST = .true.
               end if

               if (momST) then
                  ptot(j) = ptot(j) + pST                               !------- Post Sedov Taylor momentum
               else
                  ptot(j) = ptot(j) + pII                               !------- initial blastwave momentum
               end if

               !write(*,*) "Ptot",ptot(j)*(scale_v/scale_m)*(1.d0/1.d5)*(n0**0.1176d0)*(1.d0/Zscale) !km/s
               mloss(j) = mloss(j) + mejecta/vol_loc(j)

               ! ------- Energy
               ! HK note --> no longer need nsnII here as taken care of above
               ethermal(j) = ethermal(j) + ENSN/vol_loc(j)        !------- ENSN per SNII event
               
               ! --- Reduce star particle mass
               mp(ind_part(j)) = mp(ind_part(j)) - mejecta
               if (mp(ind_part(j)) .le. 0.0) then
                  write (*, *) "mp<0 from type II sampling. Correcting..."          ! can occur for stachastic sampling and very small mstarparticle
                  mp(ind_part(j)) = minmass
               end if

               ! --- Diagnostics
               if (SNdiagnostics) then
                  write (SNunit_out, '(i7,a,I10,I3,f3.0,3e14.5,L3,7e14.5)') nstep, ' SNII', idp(ind_part(j)), ilevel, numII, &
                     & t*scale_t/Myr2sec, aexp, t1/1d6, momST, n0, meanmass, Zgas, mp(ind_part(j))/scale_m, &
                     & xp(ind_part(j), :)*scale_l/kpc2cm
               end if
            end if
         end if

         !---------------------------
         !--------------------------- Type Ia (Updated for EDGE2)
         !---------------------------
         if (.not. star_is_pop3.and.mp(ind_part(j)).gt.minmass) then
            !call SNIa(t1, t2, NumSNIa, mass_for_snIa)     !----- call SNIa model
            loc_upper_slope = var_imf_a2
            tmp_fe_over_H = LOG10(MAX((zp(ind_part(j), 1) / (0.76d0))/table_Solar_Fe,1.d-40)) ! [Fe/H]
   
            ! Recalculate upper mass IMF slope if needed
            if (imf_varies_with_metallicity) call get_upper_slope_prgomet(tmp_fe_over_H,loc_upper_slope) 
            if (imf_varies_complex) then
               loc_upper_slope = -1.d0 * part_imf_slope(ind_part(j))
            end if

            call SNIa_var_IMF(var_imf_m0,var_imf_m1,var_imf_m2,var_imf_a1,loc_upper_slope, t1, t2, NumSNIa, mass_for_snIa)
            
            ! END TODO
            NumSNIa = NumSNIa*mp0(ind_part(j))/scale_m    !----- normalise using particle mass

            ! Save for NSNS rate determined later. New sampling below to avoid identical but scaled rates.
            NumNSNS = fNSNS_Ia*NumSNIa

            if (NumSNIa .gt. 0.0) then                        !--------- Do random sampling of type Ia SNe
               call ranf(localseed, RandNum)
               numresidual = NumSNIa - int(NumSNIa)         !----- int <1 --> 0
               NumSNIa = int(NumSNIa)
               if (RandNum .lt. numresidual) then
                  NumSNIa = NumSNIa + 1
               end if
            end if

            if (NumSNIa .gt. 0.0) then
               Zscale = (max(Zgas, 0.01d0))**(-0.2d0)                    !---------- scaling from Thornton et al. 1998.
               n0 = max(unew(indp(j), 1), smallr)*scale_nH              !---------- mH/cc
               pST = pST0*((NumSNIa*(SNenergy/1.d51))**0.941d0)*(n0**(-0.1176))*Zscale         !---------- km/sST momentum Blondin et al (1998).
               !---------- Check cooling radius
               rsf = 30.*((NumSNIa*(SNenergy/1.d51))**0.29)*(n0**(-0.43d0))*((Zgas + 0.01d0)**(-0.18d0)) !------ Shell formation radius (Hopkins et al. 2013, Coiffi, et al.)
               cellsize = dx_loc2(j)*scale_l/pc2cm                      !---------- local resolution element
               if (cellsize .lt. rsf/Nrcool) then                       !----------- Kim & Ostriker criterion. Momentum: 1/3, energy: 1/10.
                  momST = .false.
               else
                  momST = .true.
               end if
               if (momST) then
                  ptot(j) = ptot(j) + pST                               !------- Post Sedov Taylor momentum
               else
                  ptot(j) = ptot(j) + pIa*NumSNIa                       !------- initial pIa blastwave momentum
               end if

               masslossIa = NumSNIa*Mremnant
               mloss(j) = mloss(j) + masslossIa/vol_loc(j)
               ethermal(j) = ethermal(j) + NumSNIa*ESNIa/vol_loc(j) 

               if (metal) then
                  !Dont add metals if pushing below min mass
                  if ((mp(ind_part(j)) - masslossIa).gt.minmass) then
                     mlossmetals(j, 1) = mlossmetals(j, 1) + NumSNIa*SNIaFe_loc/vol_loc(j)  !1=Fe EDGE2
                     mlossmetals(j, 2) = mlossmetals(j, 2) + NumSNIa*SNIaO_loc /vol_loc(j)  !2=O  EDGE2
                     mlossmetals(j, 3) = mlossmetals(j, 3) + NumSNIa*SNIaN_loc /vol_loc(j)  !3=N  EDGE2
                     mlossmetals(j, 4) = mlossmetals(j, 4) + NumSNIa*SNIaMg_loc/vol_loc(j)  !4=Mg EDGE2
                     mlossmetals(j, 5) = mlossmetals(j, 5) + NumSNIa*SNIaNe_loc/vol_loc(j)  !5=Ne EDGE2
                     mlossmetals(j, 6) = mlossmetals(j, 6) + NumSNIa*SNIaSi_loc/vol_loc(j)  !6=Si EDGE2
                     mlossmetals(j, 7) = mlossmetals(j, 7) + NumSNIa*SNIaCa_loc/vol_loc(j)  !7=Ca EDGE2
                     mlossmetals(j, 8) = mlossmetals(j, 8) + NumSNIa*SNIaC_loc /vol_loc(j)  !8=C  EDGE2
                     mlossmetals(j, 9) = mlossmetals(j, 9) + NumSNIa*SNIaS_loc /vol_loc(j)  !9=S  EDGE2
                  end if
               end if

407                  format('SNIa --  a:',F10.6,' Age (Myr):',F10.3,' # SN',I5,' Density (log10 cm^-3):',F10.3,' Star Z (abs):',F10.3,' Gas Z (sol):',F10.3) 
                     write(*,407) aexp,t1/1.d6,NINT(NumSNIa),LOG10(n0),LOG10(mett),LOG10(MAX(1.d-40,Zgas))

               ! -- Reduce star particle mass
               mp(ind_part(j)) = mp(ind_part(j)) - masslossIa
               if (mp(ind_part(j)) .le. 0.0) then
                  write (*, *) "mp<0 from type Ia sampling. Correcting..."
                  mp(ind_part(j)) = minmass
                  mloss(j) = MAX(mloss(j) - minmass,0.d0)
                  ethermal(j) = ethermal(j) - NumSNIa*ESNIa/vol_loc(j) ! no energy injection at limit
               end if
               ! --- Diagnostics
               if (SNdiagnostics) then
                  write (SNunit_out, '(i7,a,I10,I3,f3.0,3e14.5,L3,7e14.5)') nstep, ' SNIa', ind_part(j), ilevel, NumSNIa, &
                     & t*scale_t/Myr2sec, aexp, t1/1d6, momST, n0, meanmass, Zgas, mp(ind_part(j))/scale_m, &
                     & xp(ind_part(j), :)*scale_l/kpc2cm
               end if
            end if
         end if
      end if
      !-----------------------
      !-----------------------  Winds
      !-----------------------
      if (winds .and. .not. star_is_pop3.and.mp(ind_part(j)).gt.minmass) then
         !---------------------------
         !--------------------------- Low mass stars (based on Kalirai et al. 2008)
         !---------------------------
         if (mstarmin <= Mwindmax .and. mstarmax >= Mwindmin .and. mstarmin <= mstarmax) then
            ! To be consistent with BPASS and yields, mass loss now computed elsewhere
            !call AGBmassloss(mstarmin, mstarmax, theint)                        !------ Integration limits are from m_i to m_(i-1)
            !masslossW = theint*mp0(ind_part(j))
            !mloss(j) = mloss(j) + masslossW/vol_loc(j)
            
            if (metal) then                                                     !------ EDGE2
               !call IMFnum_harley(mp0(ind_part(j))/scale_m, mstarmin, mstarmax, numAGB) ! For consistency with BPASS IMF
               ! TODO(code): check that this is giving us what we want
               loc_upper_slope = var_imf_a2
               tmp_fe_over_H = LOG10(MAX((zp(ind_part(j), 1) / (0.76d0))/table_Solar_Fe,1.d-40)) ! [Fe/H]

               ! Recalculate upper mass slope if needed
               if (imf_varies_with_metallicity) call get_upper_slope_prgomet(tmp_fe_over_H,loc_upper_slope)
               if (imf_varies_complex) then
                  loc_upper_slope = -1.d0 * part_imf_slope(ind_part(j)) 
               end if

               call IMFnum_between_var(var_imf_m0,var_imf_m1,var_imf_m2,var_imf_a1,loc_upper_slope,mp0(ind_part(j))/scale_m, mstarmin, mstarmax, numAGB)
               ! END TODO
               !numAGB = IMFKroupa(mstarmin, mstarmax, mp0(ind_part(j))/scale_m) ! Number of stars in mass range. Assumes instantaneous winds.
               Fe_ejecta = (numAGB*AGB_Fe_yield(meanmass, mettM)*scale_m/vol_loc(j)) !1=Fe  NUGRID 
               O_ejecta  = (numAGB*AGB_O_yield(meanmass, mettM)*scale_m/vol_loc(j))  !2=O   NUGRID 
               N_ejecta  = (numAGB*AGB_N_yield(meanmass, mettM)*scale_m/vol_loc(j))  !3=N   NUGRID 
               Mg_ejecta = (numAGB*AGB_Mg_yield(meanmass, mettM)*scale_m/vol_loc(j)) !4=Mg  NUGRID 
               Ne_ejecta = (numAGB*AGB_Ne_yield(meanmass, mettM)*scale_m/vol_loc(j)) !5=Ne  NUGRID 
               Si_ejecta = (numAGB*AGB_Si_yield(meanmass, mettM)*scale_m/vol_loc(j)) !6=Si  NUGRID 
               Ca_ejecta = (numAGB*AGB_Ca_yield(meanmass, mettM)*scale_m/vol_loc(j)) !7=Ca  NUGRID 
               C_ejecta  = (numAGB*AGB_C_yield(meanmass, mettM)*scale_m/vol_loc(j))  !8=C   NUGRID 
               S_ejecta  = (numAGB*AGB_S_yield(meanmass, mettM)*scale_m/vol_loc(j))  !9=S   NUGRID 
               H_ejecta  = (numAGB*AGB_H_yield(meanmass, mettM)*scale_m/vol_loc(j))  !H     NUGRID 
               He_ejecta = (numAGB*AGB_He_yield(meanmass, mettM)*scale_m/vol_loc(j)) !He    NUGRID 
               tmp = H_ejecta + He_ejecta + Fe_ejecta + O_ejecta + N_ejecta + Mg_ejecta + Ne_ejecta
               tmp = tmp + Si_ejecta + Ca_ejecta + C_ejecta + S_ejecta


               masslossW = tmp*vol_loc(j)
               mloss(j) = mloss(j) + tmp
               ! Don't add metals if it will push us below the minimum mass
               if ((mp(ind_part(j)) - masslossW).gt.minmass) then
                  mlossmetals(j, 1) = mlossmetals(j, 1) + Fe_ejecta 
                  mlossmetals(j, 2) = mlossmetals(j, 2) + O_ejecta
                  mlossmetals(j, 3) = mlossmetals(j, 3) + N_ejecta
                  mlossmetals(j, 4) = mlossmetals(j, 4) + Mg_ejecta
                  mlossmetals(j, 5) = mlossmetals(j, 5) + Ne_ejecta
                  mlossmetals(j, 6) = mlossmetals(j, 6) + Si_ejecta
                  mlossmetals(j, 7) = mlossmetals(j, 7) + Ca_ejecta
                  mlossmetals(j, 8) = mlossmetals(j, 8) + C_ejecta
                  mlossmetals(j, 9) = mlossmetals(j, 9) + S_ejecta
               end if
            end if

            !           ptot(j)=ptot(j)+masslossW*vAGB or something           !Oscar-Eric for EDGE2: add AGB wind velocity. 10ish km/s
            
            ! Reduce star particle mass
            mp(ind_part(j)) = mp(ind_part(j)) - masslossW
            if (mp(ind_part(j)) .le. 0.0) then
               write (*, *) "mp<0 from AGB winds. Correcting..."
               mp(ind_part(j)) = minmass
               mloss(j) = MAX(mloss(j) - minmass,0.d0)
            end if
         end if
      end if

      ! Handle tracer particles
      if (MC_tracer) then
         call mark_yielding_particle(ind_part(j), (star_original_mass-mp(ind_part(j)))/star_original_mass)
      end if
   end do

   !----------- Inject feedback ----------------
   do j = 1, np
      if (mloss(j) > 0 .or. ethermal(j) > 0 .or. ptot(j) > 0.) then  ! -- only enter if star actually injects something
         ! Specific kinetic energy of the star injecting feedback
         tekinstar = 0.5d0*(vp(ind_part(j), 1)**2 &
               &      + vp(ind_part(j), 2)**2 &
               &      + vp(ind_part(j), 3)**2)
         if (ok(j)) then !------- Check if particle is drifter
            do ii = 1, 2    !------- Do feedback over 2x2x2 cube
               do jj = 1, 2
                  do kk = 1, 2
                     iicell = indcube2(j, ii, jj, kk)
                     if (iicell .gt. 0) then

                        ! ionisation fractions, ref, etc.
                        fractions = 0.d0
                        do ifrac = i_fractions, nvar
                           fractions(ifrac) = unew(iicell, ifrac)/unew(iicell, 1)
                        end do

                        !----- Return ejected mass and associated momentum & kinetic energy. (total conserved, see standard RAMSES)
                        if (mloss(j) > 0.) then
                           unew(iicell, 1) = unew(iicell, 1) + mloss(j)/8.0     ! -- Spread over 8 cells
                           unew(iicell, 2) = unew(iicell, 2) + mloss(j)*vp(ind_part(j), 1)/8.0
                           unew(iicell, 3) = unew(iicell, 3) + mloss(j)*vp(ind_part(j), 2)/8.0
                           unew(iicell, 4) = unew(iicell, 4) + mloss(j)*vp(ind_part(j), 3)/8.0
                           unew(iicell, 5) = unew(iicell, 5) + mloss(j)*tekinstar/8.0
                        end if
                        !------- Do *pure* momentum feedback under assumption of constant E_therm
                        tt = unew(iicell, ndim + 2)
                        tekin = 0.0d0
                        do idim = 1, ndim
                           tekin = tekin + 0.5*unew(iicell, idim + 1)**2/max(unew(iicell, 1), smallr) !------- Kinetic E
                        end do
#if NENER>0
                        do irad = 1, nener
                           tekin = tekin + unew(iicell, ndim + 2 + irad) ! remember to subtract non-thermal pressure
                        end do
#endif
                        tt = tt - tekin  !Etherm
                        if (momentum) then
                           if (ptot(j) > 0.) then
                              !------- Geomtrical factors and cell index ---------------
                              xcont = -1.0 + 2.0*(ii - 1)
                              ycont = -1.0 + 2.0*(jj - 1)
                              zcont = -1.0 + 2.0*(kk - 1)
                              contr = 8.0*(xcont**2 + ycont**2 + zcont**2)**0.5  ! -- Each cell get 1/8th of momentum
                              !------ Momentum from winds and SNe -----------------
                              vkick = scale_v*ptot(j)/8.d0/1.d5/max(unew(iicell, 1), smallr)/vol_loc(j) !------- Velocity for exactly ptot/8/mass

                              vxold = unew(iicell, 2)/max(unew(iicell, 1), smallr)
                              vyold = unew(iicell, 3)/max(unew(iicell, 1), smallr)
                              vzold = unew(iicell, 4)/max(unew(iicell, 1), smallr)

                              if (vkick .gt. vmax) then  !-------Limit momentum for stability
                                 !write(*,*) "APPLYING KICK MOMENTUM LIMITER"
                                 momx = xcont*8.*unew(iicell, 1)*vmax*1.d5/scale_v/contr  !mom density, factor of 8 is to cancel contr. Vmax is the correct velocity
                                 momy = ycont*8.*unew(iicell, 1)*vmax*1.d5/scale_v/contr  !mom density
                                 momz = zcont*8.*unew(iicell, 1)*vmax*1.d5/scale_v/contr  !mom density
                              else
                                 momx = xcont*ptot(j)/contr/vol_loc(j)
                                 momy = ycont*ptot(j)/contr/vol_loc(j)
                                 momz = zcont*ptot(j)/contr/vol_loc(j)
                              end if
                              unew(iicell, 2) = unew(iicell, 2) + momx
                              unew(iicell, 3) = unew(iicell, 3) + momy
                              unew(iicell, 4) = unew(iicell, 4) + momz
                           end if
                           if (fbsafety) then       ! for stability. maxadv --> inf = no restriction
                              vxnew = unew(iicell, 2)/max(unew(iicell, 1), smallr)
                              vynew = unew(iicell, 3)/max(unew(iicell, 1), smallr)
                              vznew = unew(iicell, 4)/max(unew(iicell, 1), smallr)

                              if (abs(vxnew) .gt. maxadv) then
                                 !write(*,*) "APPLYING VELOCITY LIMITER XXX"
                                 unew(iicell, 2) = sign(maxadv, vxnew)*unew(iicell, 1)
                              end if
                              if (abs(vynew) .gt. maxadv) then
                                 !write(*,*) "APPLYING VELOCITY LIMITER YYY"
                                 unew(iicell, 3) = sign(maxadv, vynew)*unew(iicell, 1)
                              end if
                              if (abs(vznew) .gt. maxadv) then
                                 !write(*,*) "APPLYING VELOCITY LIMITER ZZZ"
                                 unew(iicell, 4) = sign(maxadv, vznew)*unew(iicell, 1)
                              end if
                           end if

                           ! ------- All momentum is now added, calculate new Ekin and update Etot.
                           tekin = 0.0d0
                           do idim = 1, ndim
                              tekin = tekin + 0.5*unew(iicell, idim + 1)**2/max(unew(iicell, 1), smallr)   !-------  Kin E
                           end do
#if NENER>0
                           do irad = 1, nener
                              tekin = tekin + unew(iicell, ndim + 2 + irad) ! remember to add back the non-thermal pressure
                           end do
#endif
                           unew(iicell, ndim + 2) = tt + tekin
                        end if
                        !tt is old Etherm, which should not change due to momentum additions
                        if (energy) then
                           if (ethermal(j) > 0.) then
                              !------ Update thermal energy, inject in indp(j) ----------------------------------
                              Emax = tekin + Tmaxfb*unew(indp(j), 1)/scale_T2/(gamma - 1.0d0)
                              unew(indp(j), ndim + 2) = unew(indp(j), ndim + 2) + ethermal(j)/8.0 !runs through 8 times
                              unew(indp(j), ndim + 2) = min(unew(indp(j), ndim + 2), Emax)  !Oscar, safety
                           end if
                        end if
                        ! ------- Finally Return metals ----- EDGE2
                        iii = 0
                        do imet = 1, nmetals
                           unew(iicell, imetal + iii) = unew(iicell, imetal + iii) + mlossmetals(j, imet)/8.0   !------- EDGE2: now an array
                           iii = iii + 1
                        end do

                        ! Deal with fractional quantities again
                        do ifrac = i_fractions, nvar
                           unew(iicell, ifrac) = fractions(ifrac)*unew(iicell, 1)
                        end do

                     end if
                  end do
               end do
            end do
         else  !-------------- drifters. Inject on parent oct.
            iicell = indp(j)

            ! ionisation fractions, ref, etc.
            fractions = 0.d0
            do ifrac = i_fractions, nvar
               ! On virtual boundaries everything gets set to 0
               ! This is a problem for the division so use uold for
               ! the fractional quantities on ilevel - 1
               if (unew(iicell, 1).gt.0.d0) then
                  fractions(ifrac) = unew(iicell, ifrac)/unew(iicell, 1)
                  boundary_trigger = .false.
               else
                  fractions(ifrac) = uold(iicell, ifrac)/uold(iicell, 1)
                  boundary_trigger = .true.
               end if
            end do

            !----- Return ejected mass and associated kinetic energy. (total conserved, see standard RAMSES)
            if (mloss(j) > 0.) then
               unew(iicell, 1) = unew(iicell, 1) + mloss(j)
               unew(iicell, 2) = unew(iicell, 2) + mloss(j)*vp(ind_part(j), 1)
               unew(iicell, 3) = unew(iicell, 3) + mloss(j)*vp(ind_part(j), 2)
               unew(iicell, 4) = unew(iicell, 4) + mloss(j)*vp(ind_part(j), 3)
               unew(iicell, 5) = unew(iicell, 5) + mloss(j)*tekinstar
            end if
            if (energy) then
               if (ethermal(j) > 0.) then
                  Emax = tekin + Tmaxfb*unew(iicell, 1)/scale_T2/(gamma - 1.0d0)
                  unew(iicell, ndim + 2) = unew(iicell, ndim + 2) + ethermal(j)
                  unew(iicell, ndim + 2) = min(unew(iicell, ndim + 2), Emax)
               end if
            end if
            ! ------- Finally Return metals -----
            iii = 0
            do imet = 1, nmetals
               unew(iicell, imetal + iii) = unew(iicell, imetal + iii) + mlossmetals(j, imet)
               iii = iii + 1
            end do

            ! Deal with fractional quantities again
            do ifrac = i_fractions, nvar
               if (boundary_trigger) then
                  if (mloss(j) > 0.) then ! ignore pre-existing mass and assume fraction with this update
                     unew(iicell, ifrac) = fractions(ifrac)*unew(iicell, 1)
                  else ! No change in cell mass --> keep old fraction otherwise all will be 0
                     unew(iicell, ifrac) = fractions(ifrac)*uold(iicell, 1)
                  endif
               else ! update as normal
                  unew(iicell, ifrac) = fractions(ifrac)*unew(iicell, 1)
               end if
            end do

         end if
      end if
   end do

   flush (SNunit_out)   ! Ensure SN log is written to disk

end subroutine feedbk
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN, vol_gas, dq, ekBlast, ind_blast, nSN)
   use pm_commons
   use amr_commons
   use hydro_commons
   use constants, only:pc2cm
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
   !------------------------------------------------------------------------
   ! This routine average the hydro quantities inside the SN bubble
   !------------------------------------------------------------------------
   integer::ilevel, ncache, nSN, j, iSN, ind, ix, iy, iz, ngrid, iskip
   integer::i, nx_loc, igrid, info
   integer, dimension(1:nvector), save::ind_grid, ind_cell
   real(dp)::x, y, z, dr_SN, d, u, v, w, ek, u2, v2, w2, dr_cell
   real(dp)::scale, dx, dxx, dyy, dzz, dx_min, dx_loc, vol_loc, rmax2, rmax
   real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
   real(dp), dimension(1:3)::skip_loc
   real(dp), dimension(1:twotondim, 1:3)::xc
   integer, dimension(1:nSN)::ind_blast
   real(dp), dimension(1:nSN)::mSN, m_gas, vol_gas, ekBlast
   real(dp), dimension(1:nSN, 1:3)::xSN, vSN, u_gas, dq, u2Blast
#ifndef WITHOUTMPI
   real(dp), dimension(1:nSN)::m_gas_all, vol_gas_all, ekBlast_all
   real(dp), dimension(1:nSN, 1:3)::u_gas_all, dq_all, u2Blast_all
#endif
   logical, dimension(1:nvector), save::ok

   if (nSN == 0) return
   if (verbose) write (*, *) 'Entering average_SN'

   ! Mesh spacing in that level
   nx_loc = (icoarse_max - icoarse_min + 1)
   skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
   skip_loc(1) = dble(icoarse_min)
   skip_loc(2) = dble(jcoarse_min)
   skip_loc(3) = dble(kcoarse_min)
   scale = boxlen/dble(nx_loc)
   dx_min = scale*0.5D0**nlevelmax

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

   ! Maximum radius of the ejecta
   rmax = MAX(2.0d0*dx_min*scale_l/aexp, rbubble*pc2cm)
   rmax = rmax/scale_l
   rmax2 = rmax*rmax

   ! Initialize the averaged variables
   vol_gas = 0.0; dq = 0.0; u2Blast = 0.0; ekBlast = 0.0; ind_blast = -1

   ! Loop over levels
   do ilevel = levelmin, nlevelmax
      ! Computing local volume (important for averaging hydro quantities)
      dx = 0.5D0**ilevel
      dx_loc = dx*scale
      vol_loc = dx_loc**ndim
      ! Cells center position relative to grid center position
      do ind = 1, twotondim
         iz = (ind - 1)/4
         iy = (ind - 1 - 4*iz)/2
         ix = (ind - 1 - 2*iy - 4*iz)
         xc(ind, 1) = (dble(ix) - 0.5D0)*dx
         xc(ind, 2) = (dble(iy) - 0.5D0)*dx
         xc(ind, 3) = (dble(iz) - 0.5D0)*dx
      end do

      ! Loop over grids
      ncache = active(ilevel)%ngrid
      do igrid = 1, ncache, nvector
         ngrid = MIN(nvector, ncache - igrid + 1)
         do i = 1, ngrid
            ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
         end do

         ! Loop over cells
         do ind = 1, twotondim
            iskip = ncoarse + (ind - 1)*ngridmax
            do i = 1, ngrid
               ind_cell(i) = iskip + ind_grid(i)
            end do

            ! Flag leaf cells
            do i = 1, ngrid
               ok(i) = son(ind_cell(i)) == 0
            end do

            do i = 1, ngrid
               if (ok(i)) then
                  ! Get gas cell position
                  x = (xg(ind_grid(i), 1) + xc(ind, 1) - skip_loc(1))*scale
                  y = (xg(ind_grid(i), 2) + xc(ind, 2) - skip_loc(2))*scale
                  z = (xg(ind_grid(i), 3) + xc(ind, 3) - skip_loc(3))*scale
                  do iSN = 1, nSN
                     ! Check if the cell lies within the SN radius
                     dxx = x - xSN(iSN, 1)
                     dyy = y - xSN(iSN, 2)
                     dzz = z - xSN(iSN, 3)
                     dr_SN = dxx**2 + dyy**2 + dzz**2
                     dr_cell = MAX(ABS(dxx), ABS(dyy), ABS(dzz))
                     if (dr_SN .lt. rmax2) then
                        vol_gas(iSN) = vol_gas(iSN) + vol_loc
                        ! Take account for grid effects on the conservation of the
                        ! normalized linear momentum
                        u = dxx/rmax
                        v = dyy/rmax
                        w = dzz/rmax
                        ! Add the local normalized linear momentum to the total linear
                        ! momentum of the blast wave (should be zero with no grid effect)
                        dq(iSN, 1) = dq(iSN, 1) + u*vol_loc
                        dq(iSN, 2) = dq(iSN, 2) + v*vol_loc
                        dq(iSN, 3) = dq(iSN, 3) + w*vol_loc
                        u2Blast(iSN, 1) = u2Blast(iSN, 1) + u*u*vol_loc
                        u2Blast(iSN, 2) = u2Blast(iSN, 2) + v*v*vol_loc
                        u2Blast(iSN, 3) = u2Blast(iSN, 3) + w*w*vol_loc
                     end if
                     if (dr_cell .le. dx_loc/2.0) then
                        ind_blast(iSN) = ind_cell(i)
                        ekBlast(iSN) = vol_loc
                     end if
                  end do
               end if
            end do

         end do
         ! End loop over cells
      end do
      ! End loop over grids
   end do
   ! End loop over levels

#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(vol_gas, vol_gas_all, nSN, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(dq, dq_all, nSN*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(u2Blast, u2Blast_all, nSN*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(ekBlast, ekBlast_all, nSN, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   vol_gas = vol_gas_all
   dq = dq_all
   u2Blast = u2Blast_all
   ekBlast = ekBlast_all
#endif
   do iSN = 1, nSN
      if (vol_gas(iSN) > 0d0) then
         dq(iSN, 1) = dq(iSN, 1)/vol_gas(iSN)
         dq(iSN, 2) = dq(iSN, 2)/vol_gas(iSN)
         dq(iSN, 3) = dq(iSN, 3)/vol_gas(iSN)
         u2Blast(iSN, 1) = u2Blast(iSN, 1)/vol_gas(iSN)
         u2Blast(iSN, 2) = u2Blast(iSN, 2)/vol_gas(iSN)
         u2Blast(iSN, 3) = u2Blast(iSN, 3)/vol_gas(iSN)
         u2 = u2Blast(iSN, 1) - dq(iSN, 1)**2
         v2 = u2Blast(iSN, 2) - dq(iSN, 2)**2
         w2 = u2Blast(iSN, 3) - dq(iSN, 3)**2
         ekBlast(iSN) = max(0.5d0*(u2 + v2 + w2), 0.0d0)
      end if
   end do

   if (verbose) write (*, *) 'Exiting average_SN'

end subroutine average_SN
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!---------------------------------------
subroutine SNIa(t1, t2, NSNIa, MMin)
!---------------------------------------
   use amr_commons
   implicit none
   REAL(kind=8), intent(out) :: NSNIa
   REAL(kind=8), intent(in) :: t1, t2, MMin

   !---  Delay Time Distribution (DTD). Maoz & Graur (2017).
   !---  Literature normalisations
   !---  2.6d-13 Ia/yr/Msun = field DTD
   !---  1.3d-13 Ia/yr/Msun = old field DTD (Graur et al. 2014)
   !---  Greater values of 4d-13-8d-13 Ia/yr/Msun = compatible with cluster DTD
   if (SNIamodel .eq. 1) then
      if (t1 > 38.d6) then                     !------ set by MS lifetime of 8 Msun stars
         NSNIa = Ia_rate*(t1/1d9)**(-1.12)*(t2 - t1)  !*mpb(ind_part(j))/scale_m ---- is normalised outside
      else
         NSNIa = 0.0
      end if
   end if

   !-- FIRE2 Hopkins et al. 2017. Prompt - delay
   if (SNIamodel .eq. 2) then
      if (t1 > 38d6) then
         NSNIa = 5.3d-14 + 1.6d-11*exp(-((t1/1.d6 - 50.)/10.)**2/2.)
         NSNIa = NSNIa*(t2 - t1)
      else
         NSNIa = 0.0
      end if
   end if

   ! Harley minor edits to Oscar's routine
   ! Makes sure we only call the routine when 8Msun stars 
   ! are off the main sequence
   if (SNIamodel .eq. 3) then
      if (Mmin .lt. 8.d0) then                     !------ set by MS lifetime of 8 Msun stars
         NSNIa = Ia_rate*(t1/1d9)**(-1.12)*(t2 - t1)  !*mpb(ind_part(j))/scale_m ---- is normalised outside
      else
         NSNIa = 0.0
      end if
   end if

END subroutine SNIa
!###########################################################

!---------------------------------------
subroutine SNIInum(m1, m2, NSNII)
   !---------------------------------------
   implicit none
   REAL(kind=8), intent(out) :: NSNII
   REAL(kind=8), intent(in) :: m1, m2
   REAL(kind=8):: A, ind

   A = 0.2244557d0  !K01, 0.1 - 100 Msun
   !  A=0.31491d0   !Chabrier 2003, 0.5 - 100 Msun
   ind = -2.3d0

   NSNII = (-A/1.3)*(m2**(-1.3) - m1**(-1.3))

END subroutine SNIInum
!###########################################################
!---------------------------------------
subroutine SNIInum_harley_of(minitial, m1, m2, NSNII)
!---------------------------------------
!Assumes an IMF similar to BPASS with a maximum
!mass of 300 Msun
!
   implicit none
   REAL(kind=8), intent(out) :: NSNII
   REAL(kind=8), intent(in) :: minitial, m1, m2
   REAL(kind=8):: scale_fac
   REAL(kind=8):: A1=0.48297d0

   scale_fac = minitial / 2.16569096d0

   NSNII = scale_fac * (-A1/1.35d0)*(m2**(-1.35d0) - m1**(-1.35d0))
END subroutine SNIInum_harley_of
!###########################################################
!---------------------------------------
subroutine SNIanum_raiteri(m1, m2, NSNIa)
   !---------------------------------------
   implicit none
   REAL(kind=8), intent(out) :: NSNIa
   REAL(kind=8), intent(in) :: m1, m2
   REAL(kind=8):: A, Ap, N1, N2, SNIafrac

   A = 0.2244557d0  !K01, 0.1 - 100 Msun
   SNIafrac = 0.16d0 !Bergh & McClure 1994, rate of SN per century in a MW type galaxy
   ! A=0.31491d0 !Chabrier 2003, 0.5 - 100 Msun
   Ap = SNIafrac*A

   N1 = (-Ap*m1**2/3.3)*((2.*m1)**(-3.3) - (m1 + 8.)**(-3.3))  ! (eq 14 in Agertz et al. 2013)
   N2 = (-Ap*m2**2/3.3)*((2.*m2)**(-3.3) - (m2 + 8.)**(-3.3))

   NSNIa = (m2 - m1)*(N1 + N2)/2.  !Trapez. For relevant timesteps, error is epsilon

END subroutine SNIanum_raiteri

!###########################################################
!---------------------------------------
subroutine AGBmassloss(m1, m2, AGB)
   !---------------------------------------
   implicit none
   REAL(kind=8), intent(out) :: AGB
   REAL(kind=8), intent(in) :: m1, m2
   REAL(kind=8):: A, N1, N2

   A = 0.2244557d0  !K01, 0.1 - 100 Msun
   ! A=0.31491d0   !Chabrier 2003, 0.5 - 100 Msun

   N1 = (m1**(-0.3))*(0.3031/m1 - 2.97)  !Agertz et al. 2013
   N2 = (m2**(-0.3))*(0.3031/m2 - 2.97)
   AGB = A*(N2 - N1)

END subroutine AGBmassloss
!###########################################################
!------------------------------------------------
SUBROUTINE fm_w(t_1, t_2, smet, fmw)
   !-----------------------------------------------
   implicit none
   real(kind=8), intent(in)::t_1, t_2, smet
   real(kind=8), intent(out)::fmw
   real(kind=8)::a, b, ts, metalscale, imfboost

   imfboost = 1.0  !0.3143d0/0.224468d0  !K01 to Chabrier

   !--- Fitting parameters ---
   a = 0.024357d0*imfboost
   b = 0.000460697d0
   ts = 1.0d7
   metalscale = a*log(smet/b + 1.d0)
   fmw = 0.0d0

   if (t_2 .le. ts) then
      fmw = metalscale*(t_2 - t_1)/ts  !Linear mass loss, m(t2)-m(t1)
   end if

   if (t_1 .lt. ts .and. t_2 .gt. ts) then
      fmw = metalscale*(ts - t_1)/ts  !Linear mass loss, m(t2)-m(t1)
   end if

   if (t_1 .ge. ts) then
      fmw = 0.0
   end if

end SUBROUTINE fm_w

SUBROUTINE get_pop3_ageGyr_fb(mass, age_gyr)
   use amr_commons, only: dp
   implicit none
   real(dp)::mass, age_gyr, logm
   real(dp), dimension(0:3)::a_fit
   integer::ndeg = 3, i
   real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v, scale_msun

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
   scale_msun = scale_l**3*scale_d/1.989d33

   a_fit = (/0.7595398e+00, -3.7303953e+00, 1.4031973e+00, -1.7896967e-01/)

   age_gyr = 0d0
   logm = log10(mass*scale_msun)
   do i = 0, ndeg
      age_gyr = age_gyr + a_fit(i)*logm**dble(i)
   end do
   age_gyr = 10d0**age_gyr

END SUBROUTINE get_pop3_ageGyr_fb

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!------------------------------------------------
SUBROUTINE p_w(t_1, t_2, smet, momW)
   !-----------------------------------------------
   use constants, only: M_sun
   implicit none
   real(kind=8), intent(in)::t_1, t_2, smet
   real(kind=8), intent(out)::momW
   real(kind=8)::a, b, c, ts, metalscale, imfboost

   imfboost = 1.0 !0.31430400d0/0.224468d0  !K01 to Chabrier
   !--- Fitting parameters ---
   a = imfboost*1.8d46/1.0d6/M_sun !Scale to per gram
   b = 0.00961529d0
   c = 0.363086d0
   ts = 6.5d6
   momW = 0.0d0

   metalscale = a*(smet/b)**c

   if (t_2 .le. ts) then
      momW = metalscale*(t_2 - t_1)/ts  !Linear mass loss, m(t2)-m(t1)
   end if

   if (t_1 .lt. ts .and. t_2 .gt. ts) then
      momW = metalscale*(ts - t_1)/ts  !Linear mass loss, m(t2)-m(t1)
   end if

   if (t_1 .ge. ts) then
      momW = 0.0
   end if

end SUBROUTINE p_w

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!------------------------------------------------
SUBROUTINE E_w(t_1, t_2, smet, EW)
   !-----------------------------------------------
   implicit none
   real(kind=8), intent(in)::t_1, t_2, smet
   real(kind=8), intent(out)::EW
   real(kind=8)::a, b, c, ts, metalscale, imfboost

   imfboost = 1.0 !0.31430400d0/0.224468d0  !K01 to Chabrier

   !--- Fitting parameters ---
   a = imfboost*1.9d54/1.0d6/2.d33 !Scale to per gram
   b = 0.0101565d0
   c = 0.41017d0
   ts = 6.5d6
   EW = 0.0d0

   metalscale = a*(smet/b)**c

   if (t_2 .le. ts) then
      EW = metalscale*(t_2 - t_1)/ts  !Linear mass loss, m(t2)-m(t1)
   end if

   if (t_1 .lt. ts .and. t_2 .gt. ts) then
      EW = metalscale*(ts - t_1)/ts  !Linear mass loss, m(t2)-m(t1)
   end if

   if (t_1 .ge. ts) then
      EW = 0.0
   end if

end SUBROUTINE E_w

!###########################################################
!###########################################################
!###########################################################
!----------------------------------------------
SUBROUTINE agemass(time, met, mass)
   !---------------------------------------
   implicit none
   real*8::a0, a1, a2, a, b, c, zzz
   real(kind=8), intent(in)::time, met
   real(kind=8), intent(out)::mass

   !     IMPLICIT REAL*8 (A-H,L-Z)
   !
   !     Following Raiteri et al.
   !
   !     Masses: 0.6-120.0 M_sun and Z: 7e-5 to 3e-2
   !
   ! if (met .lt. 7.0d-5) then
   !    zzz = 7.0d-5
   ! else
   !    zzz = met
   ! end if
   ! if (met .gt. 3.0d-2) then
   !    zzz = 3.0d-2
   ! else
   !    zzz = met
   ! end if
   zzz = min(max(met,7.0d-5),3.0d-2)

   a0 = 10.13 + 0.07547*log10(zzz) - 0.008084*(log10(zzz))**2
   a1 = -4.424 - 0.7939*log10(zzz) - 0.1187*(log10(zzz))**2
   a2 = 1.262 + 0.3385*log10(zzz) + 0.05417*(log10(zzz))**2

   c = (-log10(time) + a0)
   b = a1
   a = a2

   if (b*b - 4.*a*c .ge. 0.0) then
      mass = -b - sqrt(b*b - 4.0*a*c)
      mass = mass/(2.0*a)
      mass = 10.0**mass
   else
      mass = 120.0
   end if

END SUBROUTINE agemass
!----------------------------------------------
SUBROUTINE agemass_bpassv2_300_interp(time, mass)
!----------------------------------------------
! Harley fit to the bpass age mass relation
! There is a dependence on metallicity once we 
! hit solar but we are going to ignore this
! for simplicity. 
! Time --> Myr
! Mass --> Msun
   use hydro_parameters, only: bpass_age, bpass_mass
   implicit none
   real(kind=8), intent(in)::time
   real(kind=8), intent(out)::mass
   real(kind=8)::log_time_yr, frac_lo, frac_hi
   integer::idx

   ! Fix age < 0 problem --> only impacts the first time step
   if (time.le.0.d0) then
      mass = 300.d0
      return
   end if

   log_time_yr = LOG10(time)

   ! enforce bounds
   if (log_time_yr.le.bpass_age(1)) then
      mass = 300.d0
      return
   end if

   if (log_time_yr.ge.bpass_age(300)) then
      mass = 10.d0**bpass_mass(300)
      return
   end if

   idx = FLOOR(299.d0 * ((log_time_yr - bpass_age(1)) / (bpass_age(300) - bpass_age(1)))) + 1

   ! Due to the fine interpolation, catch the occasional stragglers'
   ! The error is usually well less than 1%
   if (log_time_yr.lt.bpass_age(idx)) then
      if ((1.d0 - (log_time_yr/bpass_age(idx))).lt.1.d-2) then
         mass = bpass_mass(idx)
      else
         write(*,*) "BROKE BPASS INTERP LOW",log_time_yr,idx,bpass_age(idx)
         mass = bpass_mass(idx)
      end if
   else if (log_time_yr.gt.bpass_age(idx+1)) then
      if ((1.d0 - (bpass_age(idx+1)/log_time_yr)).lt.1.d-2) then
         mass = bpass_mass(idx+1)
      else
         write(*,*) "BROKE BPASS INTERP HIGH",log_time_yr,idx,bpass_age(idx+1)
         mass = bpass_mass(idx+1)
      end if
   else
      frac_lo = (bpass_age(idx+1) - log_time_yr) / (bpass_age(idx+1) - bpass_age(idx))
      frac_hi = 1.d0 - frac_lo

      mass = (frac_lo * bpass_mass(idx)) + (frac_hi * bpass_mass(idx+1))
   end if
   mass = 10.d0 ** mass
   
   ! Enforce mass upper limit
   mass = MIN(mass,300.d0)

END SUBROUTINE agemass_bpassv2_300_interp

!---------------------------------------
FUNCTION IMFKroupa(m1, m2, mtot) ! EDGE2 metal enrichment
   !---------------------------------------
   ! Only high mass end of Kroupa 2001 IMF. Used for AGB winds.
   implicit none
   real(kind=8)::a, K, IMFKroupa
   real(kind=8), intent(in)::m1, m2, mtot
   a = 2.3d0
   K = mtot/3.9098880134587444d0
   IMFKroupa = K/(1.d0 - a)*(m2**(1.d0 - a) - m1**(1.d0 - a))
END FUNCTION IMFKroupa

!---------------------------------------
subroutine IMFnum_harley(minitial, m1, m2, Nstars)
!---------------------------------------
!Assumes an IMF similar to BPASS with a maximum
!mass of 300 Msun
!
   implicit none
   REAL(kind=8), intent(out) :: Nstars
   REAL(kind=8), intent(in) :: minitial, m1, m2
   REAL(kind=8):: scale_fac
   REAL(kind=8):: A1=0.48297d0

   scale_fac = minitial / 2.16569096d0

   Nstars = scale_fac * (-A1/1.35d0)*(m2**(-1.35d0) - m1**(-1.35d0))
END subroutine IMFnum_harley
