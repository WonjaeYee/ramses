!################################################################
!################################################################
!################################################################
!################################################################
subroutine single_stellar_feedback(ilevel)
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
   logical::ok_star, ok_star3, ok_ms_not_inert

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
               ok_ms_not_inert = is_pop_II(typep(ipart)) .and. (pre_ms(ipart).gt.1000)
               if (ok_ms_not_inert .or. ok_star3) then
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
               ok_ms_not_inert = is_pop_II(typep(ipart)) .and. (pre_ms(ipart).gt.1000)
               ok_star = ok_ms_not_inert .or. ok_star3

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
                  call sse_feedbk(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)
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
      if (ip > 0) call sse_feedbk(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)

      if (MC_tracer) then
         call yield_tracers(icpu, ilevel, star_tracer_type)
      end if
   end do
   ! End loop over cpus

   if (MC_tracer) call post_particle_yield()

#endif

111 format('   Entering single_stellar_feedback for level ', I2)

end subroutine single_stellar_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sse_feedbk(ind_grid, ind_part, ind_grid_part, ng, np, ilevel)
   use amr_commons
   use pm_commons
   use hydro_commons
   use random
   use imf_module
   use tracer_utils, only: mark_yielding_particle
   use constants, only: M_sun, Myr2sec, pc2cm, yr2sec, Mpc2cm, kpc2cm
   use metal_yields, only: POP3_SNII_Fe_yield, POP3_HN_Fe_yield, POP3_HMHN_Fe_yield, &
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

   ! Massive star lifetime from Myr to code units
   yearscale = scale_t/yr2sec

   ! SNenergy=1.d51 !erg/SN
   ! Type II supernova specific energy from cgs to code units
   ESN = SNenergy/(10.*M_sun)/scale_v/scale_v  !energy per 10 Msun in internal units
   ENSN = SNenergy/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units
   pII = 12.0d0*3000.0d5*scale_m/scale_v !12 Msun at 3000 km/s per SNII event, calibrated to get same \dot{p} as SB99
   pST0 = 2.95d5*1.0d5*scale_m/scale_v   !---------- Into g cm/s and then into code units, Kim & Ostriker prefactor

   ESNIa = ENSN  !Oscar: assume same energy release as SNII

   ! limiter
   vmax = vmaxFB
   Tmaxfb = Tmax
   maxadv = maxadvfb*1d5/scale_v !into internal units
   minmass = 0.1*scale_m  !minimum allowed particle mass

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
      if (sse_supernovae) then
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
      end if
      !-----------------------
      !-----------------------  Winds
      !-----------------------
      if (sse_winds .and. .not. star_is_pop3.and.mp(ind_part(j)).gt.minmass) then
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


end subroutine sse_feedbk
!###########################################################
!###########################################################
!###########################################################
!###########################################################
