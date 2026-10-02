subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
#if USE_TURB==1
  use turb_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvector,1:nvar_all),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar_all
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if

#if USE_TURB==1
        if (turb .AND. turb_type/=3) then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=gg(i,idim)+fturb(ind_leaf(i),idim)
              end do
           end do
        end if
#endif

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,neul)*vol
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,neul)*vol
        end do
        do ivar=1,ndim
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,nhydro+ivar)*vol
           end do
        end do
#endif

        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
#ifdef CALIMA
  use dust_commons, only: dust_tva, ndust, tva_test_mode, TVA_TEST_DIFFUSE
#endif
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg

  real(dp)::dtcell,smallp
  integer::k,idim
#if NENER>0
  integer::irad
#endif

#ifdef CALIMA
  real(dp)::eps_total, rho_gas
  integer :: id
#endif

  smallp = smallc**2/gamma

  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
  end do
  ! Velocity
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/uu(k,1)
     end do
  end do
  ! Internal energy
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,neul) = uu(k,neul)-half*uu(k,1)*uu(k,idim+1)**2
     end do
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,neul) = uu(k,neul)-uu(k,nhydro+irad)
     end do
  end do
#endif

  ! Debug
  if(debug)then
     do k = 1, ncell
        if(uu(k,neul).le.0.or.uu(k,1).le.smallr)then
           write(*,*)'stop in cmpdt'
           write(*,*)'dx   =',dx
           write(*,*)'ncell=',ncell
           write(*,*)'rho  =',uu(k,1)
           write(*,*)'P    =',uu(k,neul)
           write(*,*)'vel  =',uu(k,2:ndim+1)
           stop
        end if
     end do
  end if

  ! Compute pressure
  ! NOTE: uu(k,1) still holds the original mixture density here (not yet
  !       overwritten by the gravity-ratio step below), so it can safely
  !       be used in place of the removed rho_save cache.
  do k = 1, ncell
#ifndef CALIMA
     uu(k,neul) = max((gamma-one)*uu(k,neul),uu(k,1)*smallp)
#else
     if (dust_tva .and. ndust>0) then
        eps_total=0.0d0
        do id = 1,ndust
           eps_total=eps_total+uu(k,idust+id-1)/uu(k,1)
        end do
        eps_total = min(max(eps_total, 0.0_dp), 0.999_dp) ! Prevent division by zero if 100% dust
        rho_gas = uu(k,1)*(1.0d0-eps_total)
        uu(k,neul) = max((gamma-one)*uu(k,neul),rho_gas*smallp)
     else
        uu(k,neul) = max((gamma-one)*uu(k,neul),uu(k,1)*smallp)
     end if
#endif
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,nhydro+irad) = (gamma_rad(irad)-1)*uu(k,nhydro+irad)
     end do
  end do
#endif

  ! Compute sound speed
  do k = 1, ncell
     uu(k,neul) = gamma*uu(k,neul)
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,neul) = uu(k,neul) + gamma_rad(irad)*uu(k,nhydro+irad)
     end do
  end do
#endif
  do k = 1, ncell
#ifndef CALIMA
     uu(k,neul)=sqrt(uu(k,neul)/uu(k,1))
#else
      if (dust_tva .and. ndust>0) then
         eps_total=0.0d0
         do id = 1,ndust
            eps_total=eps_total+uu(k,idust+id-1)/uu(k,1)
         end do
         eps_total = min(max(eps_total, 0.0_dp), 0.999_dp) ! Prevent division by zero if 100% dust
         rho_gas = uu(k,1)*(1.0d0-eps_total)
         if (tva_test_mode == TVA_TEST_DIFFUSE) then
            uu(k,neul) = 1.0d0
         else
            uu(k,neul)=sqrt(uu(k,neul)/rho_gas)
         end if
      else
         uu(k,neul)=sqrt(uu(k,neul)/uu(k,1))
      end if
#endif
  end do

  ! Compute wave speed
  do k = 1, ncell
     uu(k,neul)=dble(ndim)*uu(k,neul)
  end do
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,neul)=uu(k,neul)+abs(uu(k,idim+1))
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1, ncell
     uu(k,1)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,1)=uu(k,1)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     uu(k,1)=uu(k,1)*dx/uu(k,neul)**2
     uu(k,1)=MAX(uu(k,1),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell = dx/uu(k,neul)*(sqrt(one+two*courant_factor*uu(k,1))-one)/uu(k,1)
     dt = min(dt,dtcell)
  end do

  ! Dust TVA CFL is now enforced by get_dust_courant_dt, called from newdt_fine.f90.

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! TC: commented because unused by default
!subroutine check_cons(ilevel)
!  use amr_commons
!  use hydro_commons
!  use poisson_commons
!  use mpi_mod
!  implicit none
!#ifndef WITHOUTMPI
!  integer::info
!  real(kind=8),dimension(3)::comm_buffin,comm_buffout
!#endif
!  integer::ilevel
!  !----------------------------------------------------------------------
!  ! Check mass and energy conservation
!  !----------------------------------------------------------------------
!  integer::i,ivar,ind,ncache,igrid,iskip
!  integer::nleaf,ngrid
!  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
!
!  real(dp)::dx,vol
!  real(kind=8)::mass_loc,ekin_loc,eint_loc
!  real(kind=8)::mass_all,ekin_all,eint_all
!  real(dp),dimension(1:nvector,1:nvar),save::uu
!
!  if(numbtot(1,ilevel)==0)return
!  if(verbose)write(*,111)ilevel
!
!  mass_all=0.0d0; mass_loc=0.0d0
!  ekin_all=0.0d0; ekin_loc=0.0d0
!  eint_all=0.0d0; eint_loc=0.0d0
!
!  ! Mesh spacing at that level
!  dx=0.5D0**ilevel*boxlen
!  vol=dx**ndim
!
!  ! Loop over active grids by vector sweeps
!  ncache=active(ilevel)%ngrid
!  do igrid=1,ncache,nvector
!     ngrid=MIN(nvector,ncache-igrid+1)
!     do i=1,ngrid
!        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
!     end do
!
!     ! Loop over cells
!     do ind=1,twotondim
!        iskip=ncoarse+(ind-1)*ngridmax
!        do i=1,ngrid
!           ind_cell(i)=ind_grid(i)+iskip
!        end do
!
!        ! Gather leaf cells
!        nleaf=0
!        do i=1,ngrid
!           if(son(ind_cell(i))==0)then
!              nleaf=nleaf+1
!              ind_leaf(nleaf)=ind_cell(i)
!           end if
!        end do
!
!        ! Gather hydro variables
!        do ivar=1,nvar
!           do i=1,nleaf
!              uu(i,ivar)=uold(ind_leaf(i),ivar)
!           end do
!        end do
!
!        ! Compute total mass
!        do i=1,nleaf
!           mass_loc=mass_loc+uu(i,1)*vol
!        end do
!
!        ! Compute total energy
!        do i=1,nleaf
!           ekin_loc=ekin_loc+uu(i,neul)*vol
!        end do
!
!        ! Compute total internal energy
!        do i=1,nleaf
!           eint_loc=eint_loc+uu(i,neul)*vol
!        end do
!        do ivar=1,ndim
!           do i=1,nleaf
!              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
!           end do
!        end do
!#if NENER>0
!        do ivar=1,nener
!           do i=1,nleaf
!              eint_loc=eint_loc-uu(i,nhydro+ivar)*vol
!           end do
!        end do
!#endif
!     end do
!     ! End loop over cells
!  end do
!  ! End loop over grids
!
!  ! Compute global quantities
!#ifndef WITHOUTMPI
!  comm_buffin(1)=mass_loc
!  comm_buffin(2)=ekin_loc
!  comm_buffin(3)=eint_loc
!  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!  mass_all=comm_buffout(1)
!  ekin_all=comm_buffout(2)
!  eint_all=comm_buffout(3)
!#endif
!#ifdef WITHOUTMPI
!  mass_all=mass_loc
!  ekin_all=ekin_loc
!  eint_all=eint_loc
!#endif
!  mass_tot=mass_tot+mass_all
!  ekin_tot=ekin_tot+ekin_all
!  eint_tot=eint_tot+eint_all
!
!111 format('   Entering check_cons for level ',I2)
!
!end subroutine check_cons
