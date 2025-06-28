!################################################################
!################################################################
!################################################################
!################################################################
subroutine single_stellar_feedback(ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the accretion onto the simple sink particles
  ! that harley created to follow individual stars
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  logical::ok_star

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  if(ndim.ne.3) return

  ! Gather star particles only.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select star particles
              ok_star = is_pop_II(typep(ipart)) .or. is_pop_III(typep(ipart))
              ! Select only pre main-sequence star particles
              ok_star = ok_star .and. (pre_ms(ipart).gt.1000) 

              if(ok_star)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather star particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              ok_star = is_pop_II(typep(ipart)) .or. is_pop_III(typep(ipart))
              ! Select only pre main-sequence star particles
              ok_star = ok_star .and. (pre_ms(ipart).gt.1000) 

              if(ok_star)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call sse_feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call sse_feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

111 format('   Entering single stellar feedback for level ',I2)

end subroutine single_stellar_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sse_feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use sse
  use constants, only: M_sun, Myr2sec, pc2cm, yr2sec, Mpc2cm, kpc2cm
  use rt_parameters, ONLY: isCO
  use rt_cooling_module, ONLY: mO_NIST_amu, mC_NIST_amu
  use cooling_module, ONLY: twopi
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine accrete_stars_are_sinks. Each stellar particle
  ! that is an active accretor will pull mass from unew at some rate
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ich,ivar,irad
  real(dp)::dx_min,vol_min
  real(dp)::dx,dx_loc,scale,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector),save::dteff
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::msun2g=1.9891d33
  real(dp)::yearscale
  real(dp),dimension(1:3)::vv
  ! fractional abundances ; for ionisation fraction and ref, etc
  real(dp),dimension(1:nvector,1:NVAR),save::fractions
  integer::i_fractions,iii,imet
  integer::sse_index_mass,sse_index_metal,sse_index_age_0,sse_index_age_1
  real(dp)::main_seq_lt_yrs,dt_winds,E_wind,v_wind
  logical::ok_SN,ok_WINDS
  real(dp)::ENSN,d,e,T2_initial,v2,mcell,t1,t2,time_simu
  real(dp)::tmp_upper_mass,tmp_lower_mass,tmp_age_0,tmp_age_1,tmp_dt
  real(dp)::initial_mass_msun,mass_initial_step,mass_final_step,mass_loss_in_step
  real(dp)::tmp_upper_mf,tmp_lower_mf,mf_final,metal_mass_loss
  real(dp)::input_sf,H_ejecta,He_ejecta,C_ejecta,N_ejecta,O_ejecta
  real(dp)::Ne_ejecta,Mg_ejecta,Si_ejecta,S_ejecta,Ca_ejecta,Fe_ejecta
  real(dp)::mett,star_initial_mass

  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+nmetals+1

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g
  yearscale = scale_t/yr2sec

  ! SN energy
  ENSN = SNenergy/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units
  ! Stellar wind velocity Hard coded now to 1000 km/s
  v_wind = 1000.d0 * 1.d5 ! Convert 1000 km/s to cm/s

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))
  end do

  if(use_proper_time)then
     do j=1,np
        dteff(j)=dteff(j)*aexp**2
     end do
  endif

  ! Store the fractional quantites that we don't want to change
  do ivar=i_fractions,nvar
     do j=1,np
        fractions(j,ivar) = unew(indp(j),ivar)/unew(indp(j),1)
     end do
  end do

  ! Compute the feedback 
  do j=1,np ! Begin loop over particles
    
     ! Convert uold to primitive variables
     d=max(uold(indp(j),1),smallr)
     vv(1)=uold(indp(j),2)/d
     vv(2)=uold(indp(j),3)/d
     vv(3)=uold(indp(j),4)/d
     e=uold(indp(j),5)
     v2=(vv(1)**2+vv(2)**2+vv(3)**2)
     e=e-0.5d0*d*v2
#if NENER>0
     do irad=1,nener
        e=e-uold(indp(j),5+irad)
     end do
#endif
     e=e/d
     T2_initial=(gamma-1.0)*e*scale_T2 ! Calculate the initial temperature
     mcell     = d*vol_loc(j)

     ! First calculate the age of the star at the beginning and end of the time step in years 
     if (cosmo) then
        if (use_proper_time) then    !Correct as tp is in proper, and dteff is proper. Now get to yrs
           call getAgeGyr(tp(ind_part(j)), t1)            !  End-of-dt age [Gyrs]
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

     ! Get the indices for mass and metallicity in the SSE arrays
     sse_index_metal = MOD(pre_ms(ind_part(j)),10000)
     sse_index_mass  = (pre_ms(ind_part(j)) - sse_index_metal)/10000
     if (sse_index_metal.le.0 .or. sse_index_mass.le.0) then
        write(*,*) "sse_index_metal or sse_index_metal is invalid",sse_index_metal,sse_index_mass
     endif

     ! Get the main sequence lifetime of the star particle in years
     main_seq_lt_yrs = final_ages(sse_index_metal,sse_index_mass)

     ! Check if the end of timestep age is > the main-sequence lifetime
     ok_SN = .false.
     ok_WINDS = .true.
     dt_winds = t2 - t1 ! delta t for winds in years
     if (t2.ge.main_seq_lt_yrs) then
        ok_SN = .true.
        dt_winds = main_seq_lt_yrs - t1

        ! The star is now completed to set the pre_ms value to a number > 1 and below 1000 
        pre_ms(ind_part(j)) = 42

        ! Consistency check
        if (dt_winds.le.0.d0) then
           write(*,*) "t1 is > main_seq_lt_yrs...something wrong",t1,t2,main_seq_lt_yrs
           call clean_stop
        end if

        write(*,*) "Star with initial mass",mp0(ind_part(j))*scale_msun,"Evolved off the MS in time",t1,t2,main_seq_lt_yrs
     end if

     !!! SSE Winds -- This includes OB winds and AGB winds
     ! TODO(code): Do I want to conserve linear momentum here from the star particle?
     if (sse_winds.and.ok_WINDS) then
        ! First we need to get the index at the end of the time step 
        sse_index_age_1 = FLOOR(((t1+dt_winds)/main_seq_lt_yrs)*n_sse_ages)

        ! Note that we are interpolating mass loss ger and not mass
        if (sse_index_age_1.gt.0) then
           tmp_upper_mass = all_sse_probs(sse_index_metal,sse_index_mass,sse_index_age_1,1) ! Mass in msol
        else
           tmp_upper_mass = sse_masses(sse_index_mass)
        endif
        tmp_lower_mass = all_sse_probs(sse_index_metal,sse_index_mass,sse_index_age_1+1,1) ! Mass in msol

        ! Get the bracket ages and the dt
        tmp_age_0 = (float(sse_index_age_1)/float(n_sse_ages)) * main_seq_lt_yrs
        tmp_age_1 = (float(sse_index_age_1+1)/float(n_sse_ages)) * main_seq_lt_yrs
        tmp_dt    = tmp_age_1 - tmp_age_0

        ! Now linearly interpolate --> note that this is actually the amount of mass loss
        mass_final_step = (((tmp_age_1 - (t1+dt_winds))/tmp_dt) * tmp_upper_mass) + ((((t1+dt_winds) - tmp_age_0)/tmp_dt) * tmp_lower_mass)

        ! Now calculate the total mass loss in step
        mass_loss_in_step = mass_final_step - ((mp0(ind_part(j))-mp(ind_part(j)))*scale_msun)

        ! If there is any mass loss, inject energy, mass, and metals 
        if (mass_loss_in_step.gt.0.d0) then
           E_wind = 0.5d0 * (mass_loss_in_step/scale_msun) * ((v_wind/scale_v)**2.d0)

           ! Update the internal energy with the wind energy (kinetic energy converted to thermal)
           unew(indp(j),ndim+2) = unew(indp(j),ndim+2) + (E_wind/vol_loc(j))

           ! Update the mass of the cell
           unew(indp(j),1) = unew(indp(j),1) + ((mass_loss_in_step/scale_msun) / vol_loc(j))

           ! Update the mass of the particle
           mp(ind_part(j)) = mp(ind_part(j)) - (mass_loss_in_step/scale_msun)

           ! Inject metals
           ! Loop over all metals
           do imet = 1,9
              ! Interpolate the metal mass fraction at the end of the step --> use max for index in case of 0 index
              tmp_upper_mf = all_sse_probs(sse_index_metal,sse_index_mass,MAX(sse_index_age_1,1),6+imet) ! Mass in msol
              tmp_lower_mf = all_sse_probs(sse_index_metal,sse_index_mass,MAX(sse_index_age_1,1)+1,6+imet) ! Mass in msol

              ! Get the bracket ages and the dt
              tmp_age_0 = (float(sse_index_age_1)/float(n_sse_ages)) * main_seq_lt_yrs
              tmp_age_1 = (float(sse_index_age_1+1)/float(n_sse_ages)) * main_seq_lt_yrs
              tmp_dt    = tmp_age_1 - tmp_age_0

              ! Now linearly interpolate the metal fraction
              mf_final = (((tmp_age_1 - (t1+dt_winds))/tmp_dt) * tmp_upper_mf) + ((((t1+dt_winds) - tmp_age_0)/tmp_dt) * tmp_lower_mf)
              
              ! Convert this to an actual metal mass by multiplying by the mass loss in the step
              metal_mass_loss = mass_loss_in_step * mf_final
              
              ! Add metals to the unew array --> note that we do not update ionization fractions, etc.
              unew(indp(j),imetal + sse_met_to_ramses_met_idx(imet) - 1) = unew(indp(j),imetal + sse_met_to_ramses_met_idx(imet) - 1) + ((metal_mass_loss/scale_msun)/vol_loc(j))
           end do 
        endif
     end if

     !!! SSE Supernova -- Type II
     !TODO(code): Do I want to conserve linear momentum here from the star particle?
     !TODO(code): Deal with Pop III stars here
     if (sse_supernova.and.ok_SN) then
         ! Get the initial mass in solar masses
         initial_mass_msun = mp0(ind_part(j))*scale_msun

         ! Check that the initial mass is within the range of SN that actually explode
         if (initial_mass_msun.ge.8.d0 .and. initial_mass_msun.le.25.d0) then
            ! Inject the energy
            unew(indp(j),ndim+2) = unew(indp(j),ndim+2) + (ENSN/vol_loc(j))      

            ! Inject mass and metals 
            mett = 2.09d0*zp(ind_part(j), 2) + 1.06d0*zp(ind_part(j), 1) !---- Stellar metallicity, Asplund 2009
            mett = mett / 0.02d0 ! Convert to solar units
            mett = MIN(MAX(mett,1.d-3),1.d0) ! Enforce bounds
            star_initial_mass = mp0(ind_part(j))*scale_msun

            ! Get the ejecta for each element
            ! Here we assume the velocity is the maximum --> 300 km/s
            H_ejecta  = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 1)
            He_ejecta = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 2)
            C_ejecta  = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 6)
            N_ejecta  = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 7)
            O_ejecta  = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 8)
            Ne_ejecta = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 10)
            Mg_ejecta = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 12)
            Si_ejecta = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 14)
            S_ejecta  = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 16)
            Ca_ejecta = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 20)
            Fe_ejecta = SNII_LC18_yields(mp0(ind_part(j))*scale_msun, mett, 3, 26)

            ! Calculate the total mass loss in the supernova
            metal_mass_loss = H_ejecta + He_ejecta + C_ejecta + N_ejecta + O_ejecta + &
                    & Ne_ejecta + Mg_ejecta + Si_ejecta + S_ejecta + Ca_ejecta + Fe_ejecta

            ! Make sure that the mass loss isn't more than the mass available
            ! This shouldn't happen but it is possible in principle because 
            ! the mass loss from the stellar evolutution models that we use
            ! is different from the mass loss that was used in the SN models
            input_sf = 1.d0
            if ((mp(ind_part(j))*scale_msun).lt.metal_mass_loss) then
                input_sf = (mp(ind_part(j))*scale_msun)/metal_mass_loss
            endif

            ! Update the particle mass
            mp(ind_part(j)) = mp(ind_part(j)) - (metal_mass_loss * input_sf)

            ! Now update the metals --> everything deposited in the local cell
            input_sf = input_sf / (vol_loc(j) * scale_msun)
            unew(indp(j),imetal+0) = unew(indp(j),imetal+0) + (Fe_ejecta * input_sf) ! Iron
            unew(indp(j),imetal+1) = unew(indp(j),imetal+1) + (O_ejecta  * input_sf) ! Oxygen
            unew(indp(j),imetal+2) = unew(indp(j),imetal+2) + (N_ejecta  * input_sf) ! Nitrogen
            unew(indp(j),imetal+3) = unew(indp(j),imetal+3) + (Mg_ejecta * input_sf) ! Magnesium
            unew(indp(j),imetal+4) = unew(indp(j),imetal+4) + (Ne_ejecta * input_sf) ! Neon
            unew(indp(j),imetal+5) = unew(indp(j),imetal+5) + (Si_ejecta * input_sf) ! Silicon
            unew(indp(j),imetal+6) = unew(indp(j),imetal+6) + (Ca_ejecta * input_sf) ! Calcium
            unew(indp(j),imetal+7) = unew(indp(j),imetal+7) + (C_ejecta  * input_sf) ! Carbon
            unew(indp(j),imetal+8) = unew(indp(j),imetal+8) + (S_ejecta  * input_sf) ! Sulfur
            
            ! Update the mass of the cell
            unew(indp(j),1) = unew(indp(j),1) + (metal_mass_loss * input_sf)

            ! TODO(code): Do i need to update the internal energy now that the mass has changed???

         end if 

     end if

  end do ! End loop over particles


  ! Update the rest of the passive scales so that the fractional quantities are not changed
  do j=1,np
     do ivar=i_fractions,nvar
        unew(indp(j),ivar) = fractions(j,ivar) * unew(indp(j),1)
     end do
  end do

end subroutine sse_feedbk

