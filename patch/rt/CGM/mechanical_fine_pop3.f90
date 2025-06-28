!####################################################################
!####################################################################
!####################################################################
subroutine mechanical_feedback_fine_pop3(ilevel, icount)
   use pm_commons
   use amr_commons
   use mechanical_commons
   use mpi_mod
   use metal_yields, only: POP3_SNII_Fe_yield, POP3_HN_Fe_yield, POP3_HMHN_Fe_yield, &
                           POP3_SNII_O_yield,  POP3_HN_O_yield,  POP3_HMHN_O_yield,  &
                           POP3_SNII_N_yield,  POP3_HN_N_yield,  POP3_HMHN_N_yield,  &
                           POP3_SNII_Mg_yield, POP3_HN_Mg_yield, POP3_HMHN_Mg_yield, &
                           POP3_SNII_Al_yield, POP3_HN_Al_yield, POP3_HMHN_Al_yield, &
                           POP3_SNII_Si_yield, POP3_HN_Si_yield, POP3_HMHN_Si_yield, &
                           POP3_SNII_Ca_yield, POP3_HN_Ca_yield, POP3_HMHN_Ca_yield, &
                           POP3_SNII_C_yield,  POP3_HN_C_yield,  POP3_HMHN_C_yield,  &
                           POP3_SNII_Ne_yield, POP3_HN_Ne_yield, POP3_HMHN_Ne_yield, &
                           POP3_SNII_S_yield,  POP3_HN_S_yield,  POP3_HMHN_S_yield,  &
                           POP3_SNII_H_yield,  POP3_HN_H_yield,  POP3_HMHN_H_yield,  &
                           POP3_SNII_He_yield, POP3_HN_He_yield, POP3_HMHN_He_yield
   use tracer_utils, only: pre_particle_yield, post_particle_yield, yield_tracers_with_dx, mark_yielding_particle


   implicit none
   !------------------------------------------------------------------------
   ! This routine computes the energy liberated from PopIII stars
   ! based on the input library, and inject thermal or kinetic energy to
   ! the surrounding of the young stars.
   ! This routine is called every fine time step.
   ! ind_pos_cell: position of the cell within an oct
   !------------------------------------------------------------------------
   integer::igrid, jgrid, ipart, jpart, next_part
   integer::npart1, npart2, icpu, icount, idim, ip
   integer::ind, ind_son, ind_cell, ilevel, iskip, info, nSNc, nSNc_mpi, nsn_star
   integer, dimension(1:nvector), save::ind_grid, ind_pos_cell
   real(dp)::tyoung, current_time, dteff
   real(dp)::skip_loc(1:3), scale, dx, dx_loc, vol_loc, x0(1:3)
   real(dp), dimension(1:twotondim, 1:ndim), save::xc
   real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
   real(dp)::scale_msun
   real(dp), dimension(1:twotondim), save::mw8, ew8 ! SNe
   real(dp), dimension(1:twotondim, 1:nmetals+nco), save::mzw8 ! metal enrichment
   real(dp), dimension(1:twotondim, 1:3), save:: pw8  ! SNe
   real(dp), dimension(1:nvector), save::mSN, eSN51
   real(dp), dimension(1:nvector, 1:nmetals+nco), save::mZSN ! metal entichment
   real(dp), dimension(1:nvector, 1:3), save::pSN
   real(dp)::mejecta, t_life, mass, esn_erg, M_Hecore
   real(dp), dimension(1:nmetals+nco)::mzejecta
   real(dp), parameter::msun2g = 2d33
   real(dp), parameter::gyr2s = 3.1536000d+16
   real(dp)::ttsta, ttend
   logical::ok, done_star
   ! metal enrichment
   real(dp)::mett
   integer::imet
   real(dp)::H_ejecta,He_ejecta
   real(dp)::yield_rescale

  ! MC tracer
   type(part_t) :: star_tracer_type
   star_tracer_type%family = FAM_TRACER_STAR

   if (icount == 2) return
   if (.not. hydro) return
   if (ndim .ne. 3) return
   if (numbtot(1, ilevel) == 0) return
   if (nstar_tot == 0) return

   if (MC_tracer) call pre_particle_yield()

#ifndef WITHOUTMPI
   if (myid .eq. 1) ttsta = MPI_WTIME()
#endif
   nSNc = 0

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
   scale_msun = scale_l**3*scale_d/msun2g

   ! Mesh spacing in that level
   call mesh_info(ilevel, skip_loc, scale, dx, dx_loc, vol_loc, xc)

   ! To filter out old particles and compute individual time steps
   ! NB: the time step is always the coarser level time step, since no feedback for icount=2
   if (ilevel == levelmin) then
      !dteff = dtnew(ilevel)
      !dteff = dtold(ilevel)
      dteff = dtnew(ilevel)+dtold(ilevel) ! Temporary fix harley for cosmo
   else
      !dteff = dtnew(ilevel - 1)
      !dteff = dtold(ilevel - 1)
      dteff = dtnew(ilevel)+dtold(ilevel) ! Temporary fix harley for cosmo
   end if

   if (use_proper_time) then
      current_time = texp
      dteff = dteff*aexp**2
   else
      current_time = t
   end if

   ncomm_SN = 0  ! important to initialize; number of communications (not SNe)
#ifndef WITHOUTMPI
   xSN_comm = 0d0; ploadSN_comm = 0d0; mSN_comm = 0d0
   mloadSN_comm = 0d0; mZloadSN_comm = 0d0; iSN_comm = 0
   floadSN_comm = 0d0; eloadSN_comm = 0d0; eSN_comm = 0d0
#endif

   ! Loop over cpus
   do icpu = 1, ncpu
      igrid = headl(icpu, ilevel)
      ip = 0

      ! Loop over grids
      do jgrid = 1, numbl(icpu, ilevel)
         npart1 = numbp(igrid)  ! Number of particles in the grid
         npart2 = 0

         ! Count star particles
         if (npart1 > 0) then
            do idim = 1, ndim
               x0(idim) = xg(igrid, idim) - dx - skip_loc(idim)
            end do

            ipart = headp(igrid)

            mw8 = 0d0; mzw8 = 0d0; pw8 = 0d0; ew8 = 0d0

            ! Loop over particles
            do jpart = 1, npart1
               ! Save next particle   <--- Very important !!!
               next_part = nextp(ipart)
               ok = .false.

               ! only if not yet exploded and low metallicity
               mett = 2.09d0*zp(ipart, 2) + 1.06d0*zp(ipart, 1) ! metal enrichment
               if (is_pop_III(typep(ipart)).and.is_pre_SN(typep(ipart))) then
                  call get_pop3_ageGyr(mp(ipart), t_life)
                  if (use_proper_time) then
                     tyoung = t_life*gyr2s/(scale_t/aexp**2)
                  else
                     tyoung = t_life*gyr2s/scale_t
                  end if
                  tyoung = current_time - tyoung

                  ! if tp is older than t_delay
                  if (tp(ipart) .le. tyoung) ok = .true.

                  ! if inert, collisionless BH particle, do nothing and only change the switch
                  mass = mp(ipart)*scale_msun
                  if ((mass .gt. 40.d0 .and. mass .lt. 140.d0) .or. (mass .gt. 300.d0)) then
                     if (ok) then
                        ! Not entirely consistent because there was no supernova
                        ! but equivalent in terms of implementation
                        typep(ipart)%tag = STAR_POP3_AFTER_SNe
                     end if
                     ok = .false.
                  end if
               end if
               
               if (ok) then
                  ! Find the cell index to get the position of it
                  ind_son = 1
                  do idim = 1, ndim
                     ind = int((xp(ipart, idim)/scale - x0(idim))/dx)
                     ind_son = ind_son + ind*2**(idim - 1)
                  end do
                  iskip = ncoarse + (ind_son - 1)*ngridmax
                  ind_cell = iskip + igrid
                  if (son(ind_cell) == 0) then  ! leaf cell
                     
                     mzejecta(1) = 0.d0
                     mzejecta(2) = 0.d0
                     mzejecta(3) = 0.d0
                     mzejecta(4) = 0.d0
                     mzejecta(5) = 0.d0
                     mzejecta(1) = 0.d0
                     mzejecta(7) = 0.d0
                     mzejecta(8) = 0.d0
                     mzejecta(9) = 0.d0
                     if (nco.gt.0) mzejecta(nmetals+1) = 0.d0
                     mejecta = 0.d0

                     if (mass .ge. 11.d0 .and. mass .lt. 20.d0) then
                        ! Nomoto + (2006) - normal SN
                        mzejecta(1) = POP3_SNII_Fe_yield(mass)/scale_msun !1=Fe Harley
                        mzejecta(2) = POP3_SNII_O_yield(mass)/scale_msun  !2=O  Harley
                        mzejecta(3) = POP3_SNII_N_yield(mass)/scale_msun  !3=N  Harley
                        mzejecta(4) = POP3_SNII_Mg_yield(mass)/scale_msun !4=Mg Harley
                        mzejecta(5) = POP3_SNII_Ne_yield(mass)/scale_msun !5=Ne Harley
                        mzejecta(6) = POP3_SNII_Si_yield(mass)/scale_msun !6=Si Harley
                        mzejecta(7) = POP3_SNII_Ca_yield(mass)/scale_msun !7=Ca Harley
                        mzejecta(8) = POP3_SNII_C_yield(mass)/scale_msun  !8=C  Harley
                        mzejecta(9) = POP3_SNII_S_yield(mass)/scale_msun  !9=S  Harley
                        H_ejecta    = POP3_SNII_H_yield(mass)/scale_msun  !H ejecta  Harley
                        He_ejecta   = POP3_SNII_He_yield(mass)/scale_msun !He ejecta  Harley

                        mejecta = mzejecta(1) + mzejecta(2) + mzejecta(3) + mzejecta(4) + mzejecta(5) + &
                                  mzejecta(6) + mzejecta(7) + mzejecta(8) + mzejecta(9) + &
                                  H_ejecta + He_ejecta
                        if (mejecta.gt.mp(ipart)) then
                           write(*,*) "MAJOR POP III ERROR, TOO MUCH EJECTA",mejecta*scale_msun,mass
                        end if
                        esn_erg = 1.d51
                     else if (mass .ge. 20.d0 .and. mass .le. 40.d0) then
                        ! Nomoto + (2006) - Hypernova
                        mzejecta(1) = POP3_HN_Fe_yield(mass)/scale_msun !1=Fe Harley
                        mzejecta(2) = POP3_HN_O_yield(mass)/scale_msun  !2=O  Harley
                        mzejecta(3) = POP3_HN_N_yield(mass)/scale_msun  ! 3=N  Harley
                        mzejecta(4) = POP3_HN_Mg_yield(mass)/scale_msun !4=Mg Harley
                        mzejecta(5) = POP3_HN_Ne_yield(mass)/scale_msun !5=Ne Harley
                        mzejecta(6) = POP3_HN_Si_yield(mass)/scale_msun !6=Si Harley
                        mzejecta(7) = POP3_HN_Ca_yield(mass)/scale_msun !7=Ca Harley
                        mzejecta(8) = POP3_HN_C_yield(mass)/scale_msun  !8=C  Harley
                        mzejecta(9) = POP3_HN_S_yield(mass)/scale_msun  !9=S  Harley
                        H_ejecta    = POP3_HN_H_yield(mass)/scale_msun  !H ejecta  Harley
                        He_ejecta   = POP3_HN_He_yield(mass)/scale_msun !He ejecta  Harley

                        mejecta = mzejecta(1) + mzejecta(2) + mzejecta(3) + mzejecta(4) + mzejecta(5) + &
                                  mzejecta(6) + mzejecta(7) + mzejecta(8) + mzejecta(9) + &
                                  H_ejecta + He_ejecta

                        if (mejecta.gt.mp(ipart)) then
                           write(*,*) "MAJOR POP III ERROR, TOO MUCH EJECTA",mejecta*scale_msun,mass
                        end if

                        esn_erg = 1.d51*(-13.7143d0 + 1.08571d0*mass)

                     else if (mass .ge. 140.d0 .and. mass .lt. 300.d0) then
                        ! Heger & Woosley (2002) + Nomoto
                        M_Hecore = (13.d0/24.d0)*(mass - 20.d0)
                        
                        ! replace mass with M_Hecore for heger yields
                        mzejecta(1) = POP3_HMHN_Fe_yield(mass)/scale_msun !1=Fe Harley
                        mzejecta(2) = POP3_HMHN_O_yield(mass)/scale_msun  !2=O  Harley
                        mzejecta(3) = POP3_HMHN_N_yield(mass)/scale_msun  !3=N  Harley
                        mzejecta(4) = POP3_HMHN_Mg_yield(mass)/scale_msun !4=Mg Harley
                        mzejecta(5) = POP3_HMHN_Ne_yield(mass)/scale_msun !5=Ne Harley
                        mzejecta(6) = POP3_HMHN_Si_yield(mass)/scale_msun !6=Si Harley
                        mzejecta(7) = POP3_HMHN_Ca_yield(mass)/scale_msun !7=Ca Harley
                        mzejecta(8) = POP3_HMHN_C_yield(mass)/scale_msun  !8=C  Harley
                        mzejecta(9) = POP3_HMHN_S_yield(mass)/scale_msun  !9=S  Harley
                        H_ejecta    = POP3_HMHN_H_yield(mass)/scale_msun  !H ejecta  Harley
                        He_ejecta   = POP3_HMHN_He_yield(mass)/scale_msun !He ejecta  Harley

                        mejecta = mzejecta(1) + mzejecta(2) + mzejecta(3) + mzejecta(4) + mzejecta(5) + &
                                  mzejecta(6) + mzejecta(7) + mzejecta(8) + mzejecta(9) + &
                                  H_ejecta + He_ejecta

                        ! Note that because of interpolation errors we have to be careful and rescale here
                        ! This is due to the coarse grid and the up and down nature of the yields
                        if (mejecta.gt.mp(ipart)) then
                           ! First get the scale factor
                           yield_rescale = mejecta / (0.999d0*mp(ipart))
                           ! Rescale the total ejecta
                           mejecta = mejecta / yield_rescale

                           ! Rescale the metal ejecta
                           do imet = 1, nmetals
                              mzejecta(imet) = mzejecta(imet) / yield_rescale
                           end do
                           H_ejecta = H_ejecta / yield_rescale
                           He_ejecta = He_ejecta/ yield_rescale
                           write(*,*) "Rescale Pop III Ejecta by",yield_rescale,mejecta*scale_msun,mass
                        end if

                        esn_erg = 1.d51*(5.0d0 + 1.304d0*(M_Hecore - 64.d0))
                     else
                        mejecta  = 0d0
                        mzejecta = 0d0 !NOTE: should be ok for vector
                        esn_erg  = 0d0
                     end if

                     if (log_mfb) then
399                     format('Mass',f10.3,' Ejected mass (Msun):',f10.3,' Fe=',f10.3,' O=',f10.3,' N=',f10.3,' Mg=',f10.3,' Ne=',f10.3,' Si=',f10.3,' Ca=',f10.3,' C=',f10.3,' S=',f10.3)
                        write(*,399) mass,scale_msun*mejecta,scale_msun*mzejecta(1),scale_msun*mzejecta(2),scale_msun*mzejecta(3),scale_msun*mzejecta(4),scale_msun*mzejecta(5),scale_msun*mzejecta(6),scale_msun*mzejecta(7),scale_msun*mzejecta(8),scale_msun*mzejecta(9)
                     end if

                     mw8(ind_son) = mw8(ind_son) + mejecta
                     pw8(ind_son, 1) = pw8(ind_son, 1) + mejecta*vp(ipart, 1)
                     pw8(ind_son, 2) = pw8(ind_son, 2) + mejecta*vp(ipart, 2)
                     pw8(ind_son, 3) = pw8(ind_son, 3) + mejecta*vp(ipart, 3)
                     ew8(ind_son) = ew8(ind_son) + esn_erg/1d51 ! normalised by 1d51 erg
                     if (metal) then
                        do imet = 1, nmetals
                           mzw8(ind_son, imet) = mzw8(ind_son, imet) + mzejecta(imet)
                        end do
                        if (nco.gt.0) mzw8(ind_son, nmetals+1) = mzw8(ind_son, nmetals+1) + mzejecta(nmetals+1)
                     end if

                     mejecta = min(mejecta,0.999d0*mp(ipart)) !don't subtract too much mass

                     if (MC_tracer) &
                       call mark_yielding_particle(ipart, mejecta / mp(ipart))

                     mp(ipart) = mp(ipart) - mejecta
                     typep(ipart)%tag = STAR_POP3_AFTER_SNe

                  end if
               end if
               ipart = next_part  ! Go to next particle
            end do
            ! End loop over particles

            do ind = 1, twotondim
               if (abs(mw8(ind)) > 0d0) then
                  nSNc = nSNc + 1
                  ip = ip + 1
                  ind_grid(ip) = igrid
                  ind_pos_cell(ip) = ind

                  ! collect information
                  mSN(ip) = mw8(ind)
                  if (metal) then
                     do imet = 1, nmetals
                        mZSN(ip, imet) = mzw8(ind, imet)
                     end do
                     if (nco.gt.0) mZSN(ip, nmetals+1) = mzw8(ind, nmetals+1)
                  end if
                  pSN(ip, 1) = pw8(ind, 1)
                  pSN(ip, 2) = pw8(ind, 2)
                  pSN(ip, 3) = pw8(ind, 3)
                  eSN51(ip) = ew8(ind)

                  if (ip == nvector) then
                     call mech_fine_pop3(ind_grid, ind_pos_cell, ip, ilevel, mSN, pSN, mZSN, eSN51)
                     ip = 0
                  end if
               end if
            end do

         end if
         igrid = next(igrid)   ! Go to next grid
      end do ! End loop over grids

      if (ip > 0) then
         call mech_fine_pop3(ind_grid, ind_pos_cell, ip, ilevel, mSN, pSN, mZSN, eSN51)
         ip = 0
      end if

      ! Yield tracers from stars (if any)
      if (MC_tracer) then
         call yield_tracers_with_dx(icpu, ilevel, xSNnei * dx, star_tracer_type)
      end if

   end do ! End loop over cpus

#ifndef WITHOUTMPI
   nSNc_mpi = 0
   ! Deal with the stars around the bounary of each cpu (need MPI)
   call mech_fine_mpi_pop3(ilevel)
   call MPI_ALLREDUCE(nSNc, nSNc_mpi, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
   nSNc = nSNc_mpi
   !if (myid .eq. 1 .and. nSNc > 0 .and. log_mfb) then
   !   ttend = MPI_WTIME(info)
   !   write (*, *) '--------------------------------------'
   !   write (*, *) 'Time elapsed in mechanical_fine_pop3 [sec]:', sngl(ttend - ttsta), nSNc, ncomm_SN
   !   write (*, *) '--------------------------------------'
   !end if
#endif

   if (MC_tracer) call post_particle_yield()

end subroutine mechanical_feedback_fine_pop3
!################################################################
!################################################################
!################################################################
subroutine mech_fine_pop3(ind_grid, ind_pos_cell, np, ilevel, mSN, pSN, mZSN, eSN51)
   use amr_commons
   use pm_commons
   use hydro_commons
   use mechanical_commons
   use constants, ONLY:pc2cm
#ifdef RT
   use rt_parameters, ONLY: iIons, nIons
#endif
   implicit none
   integer::np, ilevel ! actually the number of cells
   integer, dimension(1:nvector)::ind_grid, ind_pos_cell
   real(dp), dimension(1:nvector)::mSN, eSN51, floadSN
   real(dp), dimension(1:nvector, 1:nmetals+nco)::mZSN, mZloadSN !metal enrichment
   real(dp), dimension(1:nvector)::mloadSN, eloadSN
   real(dp), dimension(1:nvector, 1:3)::pSN, ploadSN
   !-----------------------------------------------------------------------
   ! This routine is called by subroutine mechanical_feedback_fine_pop3
   !-----------------------------------------------------------------------
   integer, dimension(1:nvector)::icellvec
   integer::i, j, nwco, nwco_here, idim, icell, igrid, ista, iend, ilevel2, ind_cell, ncell, irad, ii
   real(dp)::d, u, v, w, e, z, eth, ekk, Tk, d0, u0, v0, w0, dteff
   real(dp)::dx, dx_loc, scale, vol_loc, nH_cen, fleftSN
   real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
   real(dp)::scale_msun, msun2g = 2d33
   real(dp)::skip_loc(1:3), Tk0, ekk0, eth0, etot0
   real(dp), dimension(1:twotondim, 1:ndim), save::xc
   ! Grid based arrays
   real(dp), dimension(1:ndim, 1:nvector), save::xc2
   real(dp), dimension(1:nvector, 1:nSNnei), save::p_solid, ek_solid
   real(dp)::d_nei, Z_nei, dm_ejecta, vol_nei
   real(dp)::mload, vload, f_esn2, eturb
   real(dp), dimension(1:nmetals+nco)::Zload = 0d0
   real(dp)::num_sn, nH_nei, Zdepen = 1d0, f_w_cell, f_w_crit
   real(dp)::t_rad, r_rad, r_shell, m_cen, ekk_ej
   real(dp)::uavg, vavg, wavg, ul, vl, wl, ur, vr, wr
   real(dp)::d1, d2, d3, d4, d5, d6, dtot, pvar(1:nvar)
   real(dp)::vturb, vth, Mach, sig_s2, dratio, mload_cen
#ifdef RT
   real(dp), dimension(1:nIons), save::xion
#endif
   ! For stars affecting across the boundary of a cpu
   integer, dimension(1:nSNnei), save::icpuSNnei
   integer, dimension(1:nvector, 0:twondim):: ind_nbor
   logical, dimension(1:nvector, 1:nSNnei), save ::snowplough
   real(dp)::Mejtot_msun, Esntot_51
   ! fractional abundances ; for ionisation fraction and ref, etc
   real(dp), dimension(1:NVAR), save::fractions ! not compatible with delayed cooling
   integer::i_fractions
   ! metal entichment
   integer::imet

   ! starting index for passive variables except for imetal and chem
   i_fractions = imetal + nchem + nmetals
   if (nco.gt.0) i_fractions = i_fractions + nco
   if (no_metal_update) i_fractions = imetal + nchem

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
   scale_msun = scale_l**3*scale_d/msun2g

   ! Mesh variables
   call mesh_info(ilevel, skip_loc, scale, dx, dx_loc, vol_loc, xc)

   ! Record position of each cell [0.0-1.0] regardless of boxlen
   xc2 = 0d0
   do i = 1, np
      do idim = 1, ndim
         xc2(idim, i) = xg(ind_grid(i), idim) - skip_loc(idim) + xc(ind_pos_cell(i), idim)
      end do
   end do

   !======================================================================
   ! Determine p_solid before redistributing mass
   !   (momentum along some solid angle or cell)
   ! - This way is desirable when two adjacent SNe explode simulataenously.
   ! - if the neighboring cell does not belong to myid, this will be done
   !      in mech_fine_mpi
   !======================================================================
   p_solid = 0d0; ek_solid = 0d0; snowplough = .false.

   do i = 1, np
      ind_cell = ncoarse + ind_grid(i) + (ind_pos_cell(i) - 1)*ngridmax

      ! redistribute the mass/metals to the central cell
      call get_icell_from_pos(xc2(1:3, i), ilevel + 1, igrid, icell, ilevel2)

      ! Sanity Check
      if ((cpu_map(father(igrid)) .ne. myid) .or. &
          (ilevel .ne. ilevel2) .or. &
          (ind_cell .ne. icell)) then
         print *, '>>> fatal error in mech_fine'
         print *, cpu_map(father(igrid)), myid
         print *, ilevel, ilevel2
         print *, ind_cell, icell
         stop
      end if

      num_sn = eSN51(i)    ! already normalised to 1d51, and it doesn't have to be an integer
      nH_cen = uold(icell, 1)*scale_nH
      m_cen = uold(icell, 1)*vol_loc*scale_msun

      d = uold(icell, 1)
      u = uold(icell, 2)/d
      v = uold(icell, 3)/d
      w = uold(icell, 4)/d
      e = uold(icell, 5)
      e = e - 0.5d0*d*(u**2 + v**2 + w**2)
#if NENER>0
      do irad = 1, nener
         e = e - uold(icell, ndim + 2 + irad)
      end do
#endif
      Tk = e/d*scale_T2*(gamma - 1.0)*0.6

      if (Tk < 0) then
         print *, 'TKERR : mech fbk (pre-call): TK<0', TK, icell
      end if

      !==========================================
      ! estimate floadSN(i)
      !==========================================
      ! notice that f_LOAD / f_LEFT is a global parameter for SN ejecta themselves!!!
      if (loading_type .eq. 1) then
         ! based on Federrath & Klessen (2012)
         ! find the mach number of this cell
         vth = sqrt(gamma*1.38d-16*Tk/1.673d-24) ! cm/s
         ncell = 1
         icellvec(1) = icell
         call getnbor(icellvec, ind_nbor, ncell, ilevel)
         u0 = uold(icell, 2)
         v0 = uold(icell, 3)
         w0 = uold(icell, 4)
         ul = uold(ind_nbor(1, 1), 2) - u0
         ur = uold(ind_nbor(1, 2), 2) - u0
         vl = uold(ind_nbor(1, 3), 3) - v0
         vr = uold(ind_nbor(1, 4), 3) - v0
         wl = uold(ind_nbor(1, 5), 4) - w0
         wr = uold(ind_nbor(1, 6), 4) - w0
         vturb = sqrt(ul**2 + ur**2 + vl**2 + vr**2 + wl**2 + wr**2)
         vturb = vturb*scale_v ! cm/s

         Mach = vturb/vth

         ! get volume-filling density for beta=infty (non-MHD gas), b=0.4 (a stochastic mixture of forcing modes)
         sig_s2 = log(1d0 + (0.4*Mach)**2d0)
         ! s = -0.5*sig_s2;   s = ln(rho/rho0)
         dratio = max(exp(-0.5*sig_s2), 0.01) ! physicall M~100 would be hard to achieve
         floadSN(i) = min(dratio, f_LOAD)

      else
         dratio = 1d0
         floadSN(i) = f_LOAD
      end if

      if (log_mfb) then
398      format('MFB_POP3 = ', f7.3, 1x, f7.3, 1x, f5.1, 1x, f5.3, 1x, f9.5, 1x, f7.3)
         write (*, 398) log10(d*scale_nH), log10(Tk), sngl(num_sn), floadSN(i), 1./aexp - 1, log10(dx_loc*scale_l/pc2cm)
      end if

      dm_ejecta = f_LOAD*mSN(i)/dble(nSNnei)  ! per solid angle
      mload = f_LOAD*mSN(i) + uold(icell, 1)*vol_loc*floadSN(i)  ! total SN ejecta + host cell
      if (metal) then
         do imet = 1, nmetals
            Zload(imet) = (f_LOAD*mZSN(i, imet) + uold(icell, imetal + imet - 1)*vol_loc*floadSN(i))/mload
         end do
         if (nco.gt.0) Zload(nmetals+1) = (f_LOAD*mZSN(i, nmetals+1) + uold(icell, ico)*vol_loc*floadSN(i))/mload
      end if

      ! for popIII
      Mejtot_msun = mSN(i)*scale_msun
      Esntot_51 = eSN51(i) ! normalised already by 10^51 erg

      do j = 1, nSNnei
         call get_icell_from_pos(xc2(1:3, i) + xSNnei(1:3, j)*dx, ilevel + 1, igrid, icell, ilevel2)
         if (cpu_map(father(igrid)) .eq. myid) then ! if belong to myid

            Z_nei = z_ave*0.02 ! For metal=.false.
            if (ilevel > ilevel2) then ! touching level-1 cells
               d_nei = unew(icell, 1)
               if (metal) Z_nei = (2.09d0*unew(icell, imetal + 1) + 1.06d0*unew(icell, imetal))/d_nei ! metal enrichment
            else
               d_nei = uold(icell, 1)
               if (metal) Z_nei = (2.09d0*uold(icell, imetal + 1) + 1.06d0*uold(icell, imetal))/d_nei ! metal enrichment
            end if

            f_w_cell = (mload/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
            nH_nei = d_nei*scale_nH*dratio
            Zdepen = (max(0.01, Z_nei/0.02))**(expZ_pop3*2d0)

            ! transition mass loading factor (momentum conserving phase)
            ! chi_tr = 1 + f_w_crit
            f_w_crit = (A_pop3/1d4)**2d0/(f_esn*Mejtot_msun)*Esntot_51**(expE_pop3*2d0 - 1d0)*nH_nei**(expN_pop3*2d0)*Zdepen - 1d0
            f_w_crit = max(0d0, f_w_crit)

            if (f_w_cell .ge. f_w_crit) then ! radiative phase
               ! ptot = sqrt(2*chi_tr*Mejtot*(fe*Esntot))
               ! vload = ptot/(chi*Mejtot) = sqrt(2*chi_tr*fe*Esntot/Mejtot)/chi
               vload = dsqrt(2d0*f_esn*(Esntot_51*1d51)*(1d0 + f_w_crit)/(Mejtot_msun*msun2g))/scale_v/(1d0 + f_w_cell)/f_LOAD
               snowplough(i, j) = .true.
            else ! adiabatic phase
               ! ptot = sqrt(2*chi*Mejtot*(fe*Esntot))
               ! vload = ptot/(chi*Mejtot) = sqrt(2*fe*Esntot/chi/Mejtot)
               f_esn2 = 1d0 - (1d0 - f_esn)*f_w_cell/f_w_crit
               vload = dsqrt(2d0*f_esn2*(Esntot_51*1d51)/(1d0 + f_w_cell)/(Mejtot_msun*msun2g))/scale_v/f_LOAD
               snowplough(i, j) = .false.
            end if
            p_solid(i, j) = (1d0 + f_w_cell)*dm_ejecta*vload
            ek_solid(i, j) = ek_solid(i, j) + p_solid(i, j)*(vload*f_LOAD)/2d0 !ek=(m*v)*v/2, not (d*v)*v/2

            if (log_mfb_mega) then
             write(*,'(" MFBN nHcen=", f6.2," nHnei=", f6.2, " mej=", f6.2, " mcen=", f6.2, " vload=", f6.2, " lv2=",I3," fwcrit=", f6.2, " fwcell=", f6.2, " psol=",f6.2, " mload/48=",f6.2," mnei/8=",f6.2," mej/48=",f6.2)') &
              & log10(nH_cen), log10(nH_nei), log10(mSN(i)*scale_msun), log10(m_cen), log10(vload*scale_v/1d5), ilevel2 - ilevel,&
              & log10(f_w_crit), log10(f_w_cell), log10(p_solid(i, j)*scale_msun*scale_v/1d5), log10(mload*scale_msun/48),&
              & log10(d_nei*vol_loc/8d0*scale_msun), log10(dm_ejecta*scale_msun)
            end if

         end if

      end do ! loop over neighboring cells
   end do ! loop over SN cells

   !-----------------------------------------
   ! Redistribute mass from the SN cell
   !-----------------------------------------
   do i = 1, np
      icell = ncoarse + ind_grid(i) + (ind_pos_cell(i) - 1)*ngridmax
      d = uold(icell, 1)
      u = uold(icell, 2)/d
      v = uold(icell, 3)/d
      w = uold(icell, 4)/d
      e = uold(icell, 5)
#if NENER>0
      do irad = 1, nener
         e = e - uold(icell, ndim + 2 + irad)
      end do
#endif
      ekk = 0.5*d*(u**2 + v**2 + w**2)
      eth = e - ekk  ! thermal pressure

      ! ionisation fractions, ref, etc.
      do ii = i_fractions, nvar
         fractions(ii) = uold(icell, ii)/d
      end do

      mloadSN(i) = mSN(i)*f_LOAD + d*vol_loc*floadSN(i)
      if (metal) then
         do imet = 1, nmetals
            z = uold(icell, imetal + imet - 1)/d
            mZloadSN(i, imet) = mZSN(i, imet)*f_LOAD + d*z*vol_loc*floadSN(i)
         end do
         if (nco.gt.0) then
            z = uold(icell, ico)/d
            mZloadSN(i, nmetals+1) = mZSN(i, nmetals+1)*f_LOAD + d*z*vol_loc*floadSN(i)
         end if
      end if

      ! original momentum by star + gas entrained from the SN cell
      ploadSN(i, 1) = pSN(i, 1)*f_LOAD + vol_loc*d*u*floadSN(i)
      ploadSN(i, 2) = pSN(i, 2)*f_LOAD + vol_loc*d*v*floadSN(i)
      ploadSN(i, 3) = pSN(i, 3)*f_LOAD + vol_loc*d*w*floadSN(i)

      ! update the hydro variable
      fleftSN = 1d0 - floadSN(i)
      uold(icell, 1) = mSN(i)/vol_loc*f_LEFT + d*fleftSN
      uold(icell, 2) = pSN(i, 1)/vol_loc*f_LEFT + d*u*fleftSN ! make sure this is rho*v, not v
      uold(icell, 3) = pSN(i, 2)/vol_loc*f_LEFT + d*v*fleftSN
      uold(icell, 4) = pSN(i, 3)/vol_loc*f_LEFT + d*w*fleftSN
      if (metal) then
         do imet = 1, nmetals
            z = uold(icell, imetal + imet - 1)/d
            uold(icell, imetal + imet - 1) = mZSN(i, imet)/vol_loc*f_LEFT + d*z*fleftSN
         end do
         if (nco.gt.0) then 
            z = uold(icell, ico)/d
            uold(icell, ico) = mZSN(i, nmetals+1)/vol_loc*f_LEFT + d*z*fleftSN
         end if
      end if
      do ii = i_fractions, nvar
         uold(icell, ii) = fractions(ii)*uold(icell, 1)
      end do

      ! original kinetic energy of the gas entrained
      eloadSN(i) = ekk*vol_loc*floadSN(i)

      ! original thermal energy of the gas entrained (including the non-thermal part)
      eloadSN(i) = eloadSN(i) + eth*vol_loc*floadSN(i)

      ! reduce total energy as we are distributing it to the neighbours
      uold(icell, 5) = uold(icell, 5) - (ekk + eth)*floadSN(i)

      ! add the contribution from the original kinetic energy of SN
      d = mSN(i)/vol_loc
      u = pSN(i, 1)/mSN(i)
      v = pSN(i, 2)/mSN(i)
      w = pSN(i, 3)/mSN(i)
      uold(icell, 5) = uold(icell, 5) + 0.5d0*d*(u**2 + v**2 + w**2)*f_LEFT

      ! add the contribution from the original kinetic energy of SN to outflow
      eloadSN(i) = eloadSN(i) + 0.5d0*mSN(i)*(u**2 + v**2 + w**2)*f_LOAD

      ! update ek_solid
      ek_solid(i, :) = ek_solid(i, :) + eloadSN(i)/dble(nSNnei)

   end do  ! loop over SN cell

   !-------------------------------------------------------------
   ! Find and save stars affecting across the boundary of a cpu
   !-------------------------------------------------------------
   do i = 1, np

      ind_cell = ncoarse + ind_grid(i) + (ind_pos_cell(i) - 1)*ngridmax

      nwco = 0; icpuSNnei = 0
      do j = 1, nSNnei
         call get_icell_from_pos(xc2(1:3, i) + xSNnei(1:3, j)*dx, ilevel + 1, igrid, icell, ilevel2)

         if (cpu_map(father(igrid)) .ne. myid) then ! need mpi
            nwco = nwco + 1
            icpuSNnei(nwco) = cpu_map(father(igrid))
         else  ! can be handled locally
            vol_nei = vol_loc*(2d0**ndim)**(ilevel - ilevel2)
            pvar(:) = 0d0 ! temporary primitive variable
            if (ilevel > ilevel2) then ! touching level-1 cells
               pvar(1:nvar) = unew(icell, 1:nvar)
            else
               pvar(1:nvar) = uold(icell, 1:nvar)
            end if
            do ii = i_fractions, nvar ! fractional quantities that we don't want to change
               fractions(ii) = pvar(ii)/pvar(1)
            end do

            d0 = pvar(1)
            u0 = pvar(2)/d0
            v0 = pvar(3)/d0
            w0 = pvar(4)/d0
            ekk0 = 0.5d0*d0*(u0**2 + v0**2 + w0**2)
            eth0 = pvar(5) - ekk0
            Tk0 = eth0/d0*scale_T2*(gamma - 1)*0.6
#if NENER>0
            do irad = 1, nener
               eth0 = eth0 - pvar(ndim + 2 + irad)
            end do
#endif

            d = mloadSN(i)/dble(nSNnei)/vol_nei
            u = (ploadSN(i, 1)/dble(nSNnei) + p_solid(i, j)*vSNnei(1, j))/vol_nei/d
            v = (ploadSN(i, 2)/dble(nSNnei) + p_solid(i, j)*vSNnei(2, j))/vol_nei/d
            w = (ploadSN(i, 3)/dble(nSNnei) + p_solid(i, j)*vSNnei(3, j))/vol_nei/d
            pvar(1) = pvar(1) + d
            pvar(2) = pvar(2) + d*u
            pvar(3) = pvar(3) + d*v
            pvar(4) = pvar(4) + d*w
            ekk_ej = 0.5*d*(u**2 + v**2 + w**2)
            etot0 = eth0 + ekk0 + ek_solid(i, j)/vol_nei ! additional energy from SNe+entrained gas

            ! the minimum thermal energy input floor
            d = pvar(1)
            u = pvar(2)/d
            v = pvar(3)/d
            w = pvar(4)/d
            ekk = 0.5*d*(u**2 + v**2 + w**2)
            pvar(5) = max(etot0, ekk + eth0)

            Tk = (pvar(5) - ekk)/d*scale_T2*(gamma - 1)*0.6
            if (Tk < 0) then
               print *, 'TKERR: mech (post-call): Tk<0 =', Tk
               stop
            end if

#if NENER>0
            do irad = 1, nener
               pvar(5) = pvar(5) + pvar(ndim + 2 + irad)
            end do
#endif

            if (metal) then
               do imet = 1, nmetals
                  pvar(imetal + imet - 1) = pvar(imetal + imet - 1) + mzloadSN(i, imet)/dble(nSNnei)/vol_nei
               end do
               if (nco.gt.0) pvar(ico) = pvar(ico) + mzloadSN(i, nmetals+1)/dble(nSNnei)/vol_nei
            end if
            do ii = i_fractions, nvar
               pvar(ii) = fractions(ii)*pvar(1)
            end do

            ! update the hydro variable
            if (ilevel > ilevel2) then ! touching level-1 cells
               unew(icell, 1:nvar) = pvar(1:nvar)
            else
               uold(icell, 1:nvar) = pvar(1:nvar)
            end if

         end if
      end do ! loop over 48 neighbors

#ifndef WITHOUTMPI
      if (nwco > 0) then  ! for SNs across different cpu
         if (nwco > 1) then
            ! remove redundant cpu list for this SN cell
            nwco_here = nwco
            call redundant_non_1d(icpuSNnei(1:nwco), nwco_here, nwco)
         end if
         ista = ncomm_SN + 1
         iend = ista + nwco - 1

         if (iend > ncomm_max) then
            write (*, *) 'Error: increase ncomm_max in mechanical_fine.f90', ncomm_max, iend
            call clean_stop
         end if
         iSN_comm(ista:iend) = icpuSNnei(1:nwco)
         !lSN_comm (ista:iend)=ilevel
         mSN_comm(ista:iend) = mSN(i)
         mloadSN_comm(ista:iend) = mloadSN(i)
         xSN_comm(1, ista:iend) = xc2(1, i)
         xSN_comm(2, ista:iend) = xc2(2, i)
         xSN_comm(3, ista:iend) = xc2(3, i)
         eSN_comm(ista:iend) = eSN51(i)
         ploadSN_comm(1, ista:iend) = ploadSN(i, 1)
         ploadSN_comm(2, ista:iend) = ploadSN(i, 2)
         ploadSN_comm(3, ista:iend) = ploadSN(i, 3)
         floadSN_comm(ista:iend) = floadSN(i)
         eloadSN_comm(ista:iend) = eloadSN(i)
         if (metal) then
            do imet = 1, nmetals
               mZloadSN_comm(ista:iend, imet) = mZloadSN(i, imet)
            end do
            if (nco.gt.0) mZloadSN_comm(ista:iend, nmetals+1) = mZloadSN(i, nmetals+1)
         end if
         ncomm_SN = ncomm_SN + nwco

      end if
#endif

   end do ! loop over SN cell

end subroutine mech_fine_pop3
!################################################################
!################################################################
!################################################################
subroutine mech_fine_mpi_pop3(ilevel)
   use amr_commons
   use mechanical_commons
   use pm_commons
   use hydro_commons
   use mpi_mod
#ifdef RT
   use rt_parameters, ONLY: iIons, nIons
#endif
   implicit none
#ifndef WITHOUTMPI

   integer::i, j, info, nSN_tot, icpu, ncpu_send, ncpu_recv, ncc
   integer::ncell_recv, ncell_send, cpu2send, cpu2recv, tag, np
   integer::isend_sta, irecv_sta, irecv_end
   real(dp), dimension(:, :), allocatable::SNsend, SNrecv, p_solid, ek_solid
   integer, dimension(:), allocatable::list2recv, list2send
   integer, dimension(:), allocatable::reqrecv, reqsend
   integer, dimension(:, :), allocatable::statrecv, statsend
   ! SN variables
   real(dp)::dx, dx_loc, scale, vol_loc
   real(dp)::mloadSN_i, ploadSN_i(1:3), mSN_i, xSN_i(1:3), fload_i
   real(dp), dimension(1:nmetals+nco)::zloadSN_i
   real(dp)::f_esn2, d_nei, Z_nei, f_w_cell, f_w_crit, nH_nei, dratio
   real(dp)::num_sn, vload, Tk, Zdepen = 1d0, vol_nei, dm_ejecta
   real(dp)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
   real(dp)::scale_msun, msun2g = 2d33, pvar(1:nvar), eturb, etot0
   real(dp)::skip_loc(1:3), d, u, v, w, ekk, eth, d0, u0, v0, w0, eth0, ekk0, Tk0, ekk_ej
   integer::igrid, icell, ilevel, ilevel2, irad, ii
   real(dp), dimension(1:twotondim, 1:ndim), save::xc
   logical, allocatable, dimension(:, :)::snowplough
!  logical,dimension(1:nvector,1:nSNnei),save ::snowplough
#ifdef RT
   real(dp), dimension(1:nIons), save::xion
#endif
   real(dp)::Esntot_51, Mejtot_msun
   ! fractional abundances ; for ionisation fraction and ref, etc
   real(dp), dimension(1:NVAR), save::fractions ! not compatible with delayed cooling
   integer::i_fractions
   integer::imet

   if (ndim .ne. 3) return

   ! starting index for passive variables except for imetal and chem
   i_fractions = imetal + nchem + nmetals
   if (nco.gt.0) i_fractions = i_fractions + nco
   if (no_metal_update) i_fractions = imetal + nchem

   !============================================================
   ! For MPI communication
   !============================================================
   ncpu_send = 0; ncpu_recv = 0

   ncomm_SN_cpu = 0
   ncomm_SN_mpi = 0
   ncomm_SN_mpi(myid) = ncomm_SN
   ! compute the total number of communications needed
   call MPI_ALLREDUCE(ncomm_SN_mpi, ncomm_SN_cpu, ncpu,&
                    & MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
   nSN_tot = sum(ncomm_SN_cpu)
   if (nSN_tot == 0) return

   allocate (icpuSN_comm(1:nSN_tot, 1:2))
   allocate (icpuSN_comm_mpi(1:nSN_tot, 1:2))

   ! index for mpi variable
   if (myid == 1) then
      isend_sta = 0
   else
      isend_sta = sum(ncomm_SN_cpu(1:myid - 1))
   end if

   icpuSN_comm = 0
   do i = 1, ncomm_SN_cpu(myid)
      icpuSN_comm(isend_sta + i, 1) = myid
      icpuSN_comm(isend_sta + i, 2) = iSN_comm(i)
      ! iSN_comm:   local variable
      ! icpuSN_comm:  local (but extended) variable to be passed to a mpi variable
      ! icpuSN_comm_mpi: mpi variable
   end do

   ! share the list of communications
   icpuSN_comm_mpi = 0
   call MPI_ALLREDUCE(icpuSN_comm, icpuSN_comm_mpi, nSN_tot*2,&
                    & MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)

   ncell_send = ncomm_SN_cpu(myid)
   ncell_recv = count(icpuSN_comm_mpi(:, 2) .eq. myid, 1)

   ! check if myid needs to send anything
   if (ncell_send > 0) then
      allocate (SNsend(1:nvarSN, 1:ncell_send)) ! x(3),m,mz,p,pr
      allocate (list2send(1:ncell_send))
      list2send = 0; SNsend = 0d0
      list2send = icpuSN_comm_mpi(isend_sta + 1:isend_sta + ncell_send, 2)
      ncpu_send = 1
      if (ncell_send > 1) call redundant_non_1d(list2send, ncell_send, ncpu_send)
      ! ncpu_send = No. of cpus to which myid should send info
      allocate (reqsend(1:ncpu_send))
      allocate (statsend(1:MPI_STATUS_SIZE, 1:ncpu_send))
      reqsend = 0; statsend = 0
   end if

   ! check if myid needs to receive anything
   if (ncell_recv > 0) then
      allocate (SNrecv(1:nvarSN, 1:ncell_recv)) ! x(3),m,mz,p,pr
      allocate (list2recv(1:ncell_recv))
      list2recv = 0; SNrecv = 0d0
      j = 0
      do i = 1, nSN_tot
         if (icpuSN_comm_mpi(i, 2) .eq. myid) then
            j = j + 1
            list2recv(j) = icpuSN_comm_mpi(i, 1)
         end if
      end do

      ncc = j
      if (ncc .ne. ncell_recv) then ! sanity check
         write (*, *) 'Error in mech_fine_mpi: ncc != ncell_recv', ncc, ncell_recv, myid
         call clean_stop
      end if

      ncpu_recv = 1 ! No. of cpus from which myid should receive info
      if (j > 1) call redundant_non_1d(list2recv, ncc, ncpu_recv)

      allocate (reqrecv(1:ncpu_recv))
      allocate (statrecv(1:MPI_STATUS_SIZE, 1:ncpu_recv))
      reqrecv = 0; statrecv = 0
   end if

   ! prepare one variable and send
   if (ncell_send > 0) then
      do icpu = 1, ncpu_send
         cpu2send = list2send(icpu)
         ncc = 0 ! number of SN host cells that need communications with myid=cpu2send
         do i = 1, ncell_send
            j = i + isend_sta
            if (icpuSN_comm_mpi(j, 2) .eq. cpu2send) then
               ncc = ncc + 1
               SNsend(1:3, ncc) = xSN_comm(1:3, i)
               SNsend(4, ncc) = mSN_comm(i)
               SNsend(5, ncc) = mloadSN_comm(i)
               SNsend(6:8, ncc) = ploadSN_comm(1:3, i)
               SNsend(9, ncc) = floadSN_comm(i)
               SNsend(10, ncc) = eloadSN_comm(i)
               if (metal) then
                  do imet = 1, nmetals
                     SNsend(11 + imet - 1, ncc) = mZloadSN_comm(i, imet)
                  end do
                  if (nco.gt.0) SNsend(11 + nmetals, ncc) = mZloadSN_comm(i, nmetals+1)
               end if
               SNsend(11 + nmetals + nco, ncc) = eSN_comm(i)
            end if
         end do ! i

         tag = myid + cpu2send + ncc
         call MPI_ISEND(SNsend(1:nvarSN, 1:ncc), ncc*nvarSN, MPI_DOUBLE_PRECISION, &
                       & cpu2send - 1, tag, MPI_COMM_WORLD, reqsend(icpu), info)
      end do ! icpu

   end if ! ncell_send>0

   ! receive one large variable
   if (ncell_recv > 0) then
      irecv_sta = 1
      do icpu = 1, ncpu_recv
         cpu2recv = list2recv(icpu)
         ncc = 0 ! number of SN host cells that need communications with cpu2recv
         do i = 1, nSN_tot
            if ((icpuSN_comm_mpi(i, 1) == cpu2recv) .and.&
              &(icpuSN_comm_mpi(i, 2) == myid)) then
               ncc = ncc + 1
            end if
         end do
         irecv_end = irecv_sta + ncc - 1
         tag = myid + cpu2recv + ncc

! NB. Some mpi routines do not like the receiving buffer
         call MPI_IRECV(SNrecv(1:nvarSN, irecv_sta:irecv_end), ncc*nvarSN, MPI_DOUBLE_PRECISION,&
                      & cpu2recv - 1, tag, MPI_COMM_WORLD, reqrecv(icpu), info)

         irecv_sta = irecv_end + 1
      end do ! icpu

   end if ! ncell_recv >0

   if (ncpu_send > 0) call MPI_WAITALL(ncpu_send, reqsend, statsend, info)
   if (ncpu_recv > 0) call MPI_WAITALL(ncpu_recv, reqrecv, statrecv, info)

   !============================================================
   ! inject mass/metal/momentum
   !============================================================

   ! Conversion factor from user units to cgs units
   call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
   scale_msun = scale_l**3*scale_d/msun2g

   ! Mesh variables
   call mesh_info(ilevel, skip_loc, scale, dx, dx_loc, vol_loc, xc)

   np = ncell_recv
   if (ncell_recv > 0) then
      allocate (p_solid(1:np, 1:nSNnei))
      allocate (ek_solid(1:np, 1:nSNnei))
      allocate (snowplough(1:np, 1:nSNnei))
      p_solid = 0d0; ek_solid = 0d0; snowplough = .false.
   end if

   ! Compute the momentum first before redistributing mass
   do i = 1, np

      xSN_i(1:3) = SNrecv(1:3, i)
      mSN_i = SNrecv(4, i)
      mloadSN_i = SNrecv(5, i)
      fload_i = SNrecv(9, i)
      dm_ejecta = f_LOAD*mSN_i/dble(nSNnei)
      num_sn = mSN_i/(M_SNII/scale_msun)
      if (metal) then
         do imet = 1, nmetals
            ZloadSN_i(imet) = SNrecv(11 + imet - 1, i)/SNrecv(5, i)
         end do
         if (nco.gt.0) ZloadSN_i(nmetals+1) = SNrecv(11 + nmetals, i)/SNrecv(5, i)
      end if
      ek_solid(i, :) = SNrecv(10, i)/dble(nSNnei) ! kinetic energy of the gas mass entrained from the host cell + SN

      Mejtot_msun = mSN_i*scale_msun
      Esntot_51 = SNrecv(11 + nmetals + nco, i)

      do j = 1, nSNnei
         call get_icell_from_pos(xSN_i + xSNnei(1:3, j)*dx, ilevel + 1, igrid, icell, ilevel2)
         if (cpu_map(father(igrid)) .eq. myid) then ! if belong to myid
            Z_nei = z_ave*0.02 ! for metal=.false.
            if (ilevel > ilevel2) then ! touching level-1 cells
               d_nei = unew(icell, 1)
               if (metal) Z_nei = (2.09d0*unew(icell, imetal + 1) + 1.06d0*unew(icell, imetal))/d_nei
            else
               d_nei = uold(icell, 1)
               if (metal) Z_nei = (2.09d0*uold(icell, imetal + 1) + 1.06d0*uold(icell, imetal))/d_nei
            end if
            if (loading_type .eq. 1 .and. fload_i < f_LOAD) then
               dratio = fload_i
            else
               dratio = 1d0
            end if
            f_w_cell = (mloadSN_i/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
            nH_nei = d_nei*scale_nH*dratio
            Zdepen = (max(0.01, Z_nei/0.02))**(expZ_pop3*2d0) !From Thornton+(98)

            ! transition mass loading factor (momentum conserving phase)
            f_w_crit = (A_pop3/1d4)**2d0/(f_esn*Mejtot_msun)*Esntot_51**(expE_pop3*2d0 - 1d0)*nH_nei**(expN_pop3*2d0)*Zdepen - 1d0
            f_w_crit = max(0d0, f_w_crit)

            if (f_w_cell .ge. f_w_crit) then ! radiative phase
               vload = dsqrt(2d0*f_esn*(Esntot_51*1d51)*(1d0 + f_w_crit)/(Mejtot_msun*msun2g))/scale_v/(1d0 + f_w_cell)/f_LOAD
               snowplough(i, j) = .true.
            else
               f_esn2 = 1d0 - (1d0 - f_esn)*f_w_cell/f_w_crit
               vload = dsqrt(2d0*f_esn2*(Esntot_51*1d51)/(1d0 + f_w_cell)/(Mejtot_msun*msun2g))/scale_v/f_LOAD
               snowplough(i, j) = .false.
            end if
            p_solid(i, j) = (1d0 + f_w_cell)*dm_ejecta*vload
            ek_solid(i, j) = ek_solid(i, j) + p_solid(i, j)*(vload*f_LOAD)/2d0 !ek
         end if

      end do ! loop over neighboring cells

   end do

   ! Apply the SNe to this cpu domain
   do i = 1, np

      xSN_i(1:3) = SNrecv(1:3, i)
      mSN_i = SNrecv(4, i)
      mloadSN_i = SNrecv(5, i)
      ploadSN_i(1:3) = SNrecv(6:8, i)
      if (metal) then
         do imet = 1, nmetals
            ZloadSN_i(imet) = SNrecv(11 + imet - 1, i)/SNrecv(5, i)
         end do
         if (nco.gt.0) ZloadSN_i(nmetals+1) = SNrecv(11 + nmetals, i)/SNrecv(5, i)
      end if

      do j = 1, nSNnei

         call get_icell_from_pos(xSN_i(1:3) + xSNnei(1:3, j)*dx, ilevel + 1, igrid, icell, ilevel2)
         if (cpu_map(father(igrid)) == myid) then

            vol_nei = vol_loc*(2d0**ndim)**(ilevel - ilevel2)
            if (ilevel > ilevel2) then ! touching level-1 cells
               pvar(1:nvar) = unew(icell, 1:nvar)
            else
               pvar(1:nvar) = uold(icell, 1:nvar)
            end if
            do ii = i_fractions, nvar
               fractions(ii) = pvar(ii)/pvar(1)
            end do

            d0 = pvar(1)
            u0 = pvar(2)/d0
            v0 = pvar(3)/d0
            w0 = pvar(4)/d0
            ekk0 = 0.5d0*d0*(u0**2d0 + v0**2d0 + w0**2d0)
            eth0 = pvar(5) - ekk0
#if NENER>0
            do irad = 1, nener
               eth0 = eth0 - pvar(ndim + 2 + irad)
            end do
#endif
            Tk0 = eth0/d0*scale_T2*(gamma - 1)*0.6

            if (Tk0 < 0) then
               print *, 'TKERR : part1 mpi', myid, Tk0, d0*scale_nH
               stop
            end if

            d = mloadSN_i/dble(nSNnei)/vol_nei
            u = (ploadSN_i(1)/dble(nSNnei) + p_solid(i, j)*vSNnei(1, j))/vol_nei/d
            v = (ploadSN_i(2)/dble(nSNnei) + p_solid(i, j)*vSNnei(2, j))/vol_nei/d
            w = (ploadSN_i(3)/dble(nSNnei) + p_solid(i, j)*vSNnei(3, j))/vol_nei/d

            pvar(1) = pvar(1) + d
            pvar(2) = pvar(2) + d*u
            pvar(3) = pvar(3) + d*v
            pvar(4) = pvar(4) + d*w
            ekk_ej = 0.5*d*(u**2 + v**2 + w**2)
            etot0 = eth0 + ekk0 + ek_solid(i, j)/vol_nei  ! additional energy from SNe+entrained gas

            ! the minimum thermal energy input floor
            d = pvar(1)
            u = pvar(2)/d
            v = pvar(3)/d
            w = pvar(4)/d
            ekk = 0.5*d*(u**2 + v**2 + w**2)

            pvar(5) = max(etot0, ekk + eth0)

#if NENER>0
            do irad = 1, nener
               pvar(5) = pvar(5) + pvar(ndim + 2 + irad)
            end do
#endif

            if (metal) then
               do imet = 1, nmetals
                  pvar(imetal + imet - 1) = pvar(imetal + imet - 1) + mloadSN_i/dble(nSNnei)*ZloadSN_i(imet)/vol_nei
               end do
               if (nco.gt.0) pvar(ico) = pvar(ico) + mloadSN_i/dble(nSNnei)*ZloadSN_i(nmetals+1)/vol_nei
            end if
            do ii = i_fractions, nvar
               pvar(ii) = fractions(ii)*pvar(1)
            end do

            ! update the hydro variable
            if (ilevel > ilevel2) then ! touching level-1 cells
               unew(icell, 1:nvar) = pvar(1:nvar)
            else
               uold(icell, 1:nvar) = pvar(1:nvar)
            end if

         end if ! if this belongs to me

      end do ! loop over neighbors
   end do ! loop over SN cells

   deallocate (icpuSN_comm_mpi, icpuSN_comm)
   if (ncell_send > 0) deallocate (list2send, SNsend, reqsend, statsend)
   if (ncell_recv > 0) deallocate (list2recv, SNrecv, reqrecv, statrecv)
   if (ncell_recv > 0) deallocate (p_solid, ek_solid, snowplough)

   ncomm_SN = nSN_tot
#endif

end subroutine mech_fine_mpi_pop3
!################################################################
!################################################################
!################################################################
SUBROUTINE get_pop3_ageGyr(mass, age_gyr)
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

END SUBROUTINE get_pop3_ageGyr
!################################################################
!################################################################
!################################################################
