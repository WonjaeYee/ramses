!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_stars_are_sinks(ilevel)
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
              ok_star = ok_star .and. (pre_ms(ipart).eq.1)

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
              ok_star = ok_star .and. (pre_ms(ipart).eq.1)

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
                 call accrete_sas(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call accrete_sas(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

111 format('   Entering stlelar winds for level ',I2)

end subroutine accrete_stars_are_sinks 
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_sas(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
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
  real(dp)::msun2g=2d33
  real(dp)::mcell,d,e,tstar,factG,sfr_ff,sf_mass,delta_e_tot,tmp,prefac
  real(dp)::v2,initial_total_mass,final_total_mass,T2_initial,T2_final
  real(dp),dimension(1:3)::initial_total_momentum,final_total_momentum,initial_total_momentum2
  real(dp),dimension(1:3)::vv
  ! fractional abundances ; for ionisation fraction and ref, etc
  real(dp),dimension(1:nvector,1:NVAR),save::fractions
  integer::i_fractions
  integer::imet

  write(*,*) "Entering Harley sink accretion"
  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+nmetals+1

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/1.989d33

  factG = 1.d0
  if (cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

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

  ! Compute the accretion 
  do j=1,np
    
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
     ! this is the free fall time of an homogeneous sphere
     tstar     = 0.5427*sqrt(1.0/(factG*d))
     sfr_ff    = 0.3d0 ! Fix to 30% now
     ! amount of mass that's available for star formation
     ! NOTE we might want to change sfr_ff to another number
     sf_mass   = dtnew(ilevel)*sfr_ff/tstar*mcell

     ! Safety to prevent too much accretion
     ! TODO(code) decide what's best for this number
     sf_mass   = MIN(sf_mass,0.5d0*mcell) 

     !write(*,*) "Cell mass: ",mcell*scale_msun
     !write(*,*) "Time step: ",dtnew(ilevel)*scale_t
     !write(*,*) "Accreted mass: ",sf_mass*scale_msun
     initial_total_mass = uold(indp(j),1)*vol_loc(j)+mp_c(ind_part(j))
     do irad=1,3
        initial_total_momentum(irad) = mcell * vv(irad) ! cell momentum
        initial_total_momentum(irad) = initial_total_momentum(irad) + (mp_c(ind_part(j)) * vp(ind_part(j),irad)) ! particle momentum 
        initial_total_momentum2(irad) = unew(indp(j),irad+1) * vol_loc(j) ! cell momentum
        initial_total_momentum2(irad) = initial_total_momentum2(irad) + (mp_c(ind_part(j)) * vp(ind_part(j),irad)) ! particle momentum
     end do
     ! Now that we have computed the amount of mass that will
     ! go into star formation in this time step, we need to
     ! add that mass (and metals) to the star particle

     ! First check that adding sf_mass to the particle won't make
     ! it become more massive than the the final mass we are
     ! targetting
     if ((sf_mass+mp_c(ind_part(j))).gt.mp(ind_part(j))) then
        sf_mass = mp(ind_part(j)) - mp_c(ind_part(j))

        ! Since we have reached the final mass can update
        ! relevant variables

        ! Make sure the particle is no longer on the
        ! pre-main-sequence
        pre_ms(ind_part(j)) = 0

        ! Set the pre main-sequence time to the time
        ! that the particle was formed
        tp_ms(ind_part(j))  = tp(ind_part(j))

        ! Reset the birth time to now --> needed for stellar 
        ! evolution
        if(use_proper_time)then
           tp(ind_part(j)) = texp
        else
           tp(ind_part(j)) = t
        endif 
     end if

     ! Now we accrete the mass and metals from the host cell

     ! Make sure we also update the momentum of the particle
     ! so that the momentum is conserved
     if(metal)then
        do imet=1,nmetals
           if (isCO) then
              if (imet.eq.2) then
                 zp(ind_part(j),imet) = ((zp(ind_part(j),imet)*mp_c(ind_part(j))) + ((unew(indp(j),imetal+imet-1)+(mO_NIST_amu/(mC_NIST_amu+mO_NIST_amu))*unew(indp(j),iCO))*sf_mass))/(mp_c(ind_part(j))+sf_mass)
              else if (imet.eq.8) then
                 zp(ind_part(j),imet) = ((zp(ind_part(j),imet)*mp_c(ind_part(j))) + ((unew(indp(j),imetal+imet-1)+(mC_NIST_amu/(mC_NIST_amu+mO_NIST_amu))*unew(indp(j),iCO))*sf_mass))/(mp_c(ind_part(j))+sf_mass) 
              else
                 zp(ind_part(j),imet) = ((zp(ind_part(j),imet)*mp_c(ind_part(j))) + (unew(indp(j),imetal+imet-1)*sf_mass))/(mp_c(ind_part(j))+sf_mass)
              endif
           else
              zp(ind_part(j),imet) = ((zp(ind_part(j),imet)*mp_c(ind_part(j))) + (unew(indp(j),imetal+imet-1)*sf_mass))/(mp_c(ind_part(j))+sf_mass)
           endif
        enddo
     endif

     ! Momentum conservation --> update sink velocity
     vp(ind_part(j),1) = ((vp(ind_part(j),1)*mp_c(ind_part(j))) + (vv(1)*sf_mass))/(mp_c(ind_part(j))+sf_mass)
     vp(ind_part(j),2) = ((vp(ind_part(j),2)*mp_c(ind_part(j))) + (vv(2)*sf_mass))/(mp_c(ind_part(j))+sf_mass)
     vp(ind_part(j),3) = ((vp(ind_part(j),3)*mp_c(ind_part(j))) + (vv(3)*sf_mass))/(mp_c(ind_part(j))+sf_mass)

     ! Mass conservation
     mp_c(ind_part(j)) = mp_c(ind_part(j)) + sf_mass

     ! Calculate the change in energy
     delta_e_tot = (sf_mass*e) + (0.5d0*sf_mass*v2)

     ! Update density 
     unew(indp(j),1) = unew(indp(j),1) - sf_mass/vol_loc(j)

     ! Update velocity
     unew(indp(j),2:4) = unew(indp(j),2:4) - vv(1:3)*sf_mass/vol_loc(j)
     
     ! Update energy 
     e = e * unew(indp(j),1)
     do irad=1,3
        e = e + (0.5d0 * ((unew(indp(j),irad+1)/unew(indp(j),1))**2.d0) * unew(indp(j),1))
     end do
#if NENER>0
     do irad=1,nener
        e=e-uold(indp(j),5+irad)
     end do
#endif
     unew(indp(j),5) = e
    
     ! Update the metals 
     do ivar=imetal,imetal+nmetals+nco-1
        unew(indp(j),ivar)=unew(indp(j),ivar)-uold(indp(j),ivar)*sf_mass/(d*vol_loc(j))
     end do

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!! CONSISTENCY CHECKS BELOW HERE !!!!!   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Is mass conserved --> check initial and final
     final_total_mass = unew(indp(j),1)*vol_loc(j)+mp_c(ind_part(j))
     if(ABS(initial_total_mass-final_total_mass)/initial_total_mass.gt.1.d-8) then
        write(*,*) "Mass not properly conserved in simple sink accretion"
     endif

     ! Is momentum conserved --> check initial and final
     do irad=1,3
        final_total_momentum = unew(indp(j),irad+1) * vol_loc(j) ! cell momentum
        final_total_momentum = final_total_momentum + (mp_c(ind_part(j))*vp(ind_part(j),irad)) !particle momentum
        if(ABS(initial_total_momentum(irad)-final_total_momentum(irad))/initial_total_momentum(irad).gt.1.d-8) then
           write(*,*) "Momentum not properly conserved in simple sink accretion",i
        endif
     end do

     ! Is the gas temperature the same before and after
     T2_final = unew(indp(j),5)
     ! Subtract off kinetic energy
     do irad=1,3
        T2_final = T2_final - (0.5d0 * (unew(indp(j),irad+1)**2.d0)/unew(indp(j),1))
     end do
     ! Subtract off non-thermal energy
#if NENER>0
     do irad=1,nener
           T2_final = T2_final - unew(indp(j),5+irad)
     end do
#endif
     T2_final = (gamma-1)*(T2_final/unew(indp(j),1))*scale_T2
     if (ABS(T2_final-T2_initial)/T2_initial.gt.1.d-8) then
        write(*,*) "Thermal energy not properly conserved in simple sink accretion",T2_final,T2_initial
     endif

     ! Are the metal fractions ok before and after
     do ivar=imetal,imetal+nmetals+nco-1
        if (ABS( (uold(indp(j),ivar)/uold(indp(j),1)) - (unew(indp(j),ivar)/unew(indp(j),1)) ) / (uold(indp(j),ivar)/uold(indp(j),1)) .gt. 1.d-8) then
           write(*,*) "Metals not properly conserved in simple sink accretion",ivar
        endif 
     end do
  end do

  ! Update the rest of the passive scales so that the fractional quantities are not changed
  do j=1,np
     do ivar=i_fractions,nvar
        unew(indp(j),ivar) = fractions(j,ivar) * unew(indp(j),1)
     end do
  end do

end subroutine accrete_sas

