!################################################################
!################################################################
!################################################################
!################################################################
subroutine rt_make_boundary_hydro(ilevel)
  use amr_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor=1
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::switch,dx,dx_loc,scale
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nrtvar),save::uu

  integer::rtType !RT-- Type of rt variable =0,1,2,3

  !--------------------------------------------------------------------
  ! LOCAL PATCH for tests/dust_tva/dustylevatm (Makefile PATCH = patch).
  ! Adds an EMITTING lower boundary along the height direction, following
  ! Rosdahl & Teyssier 2015 eq. 83, which is what their sec. 3.7 dusty
  ! atmosphere test needs and what patch/rt/davis does for the same test.
  !
  ! A plain Dirichlet ghost, or an rt_nsource region, injects F_* but also
  ! destroys radiation arriving from above, so during a trapped-photon
  ! pile-up the NET injected flux collapses (measured: 18% of F_*). Setting
  ! the ghost energy density to
  !
  !     c~E_0 = F_* - F_1 + c~E_1
  !
  ! together with F_0 = F_* makes the GLF intercell flux at the boundary face
  ! exactly F_* whatever the interior state is -- a genuinely emitting wall.
  !
  ! F_* is carried by rt_n_bound(ibound) as a photon NUMBER flux
  ! [photons cm-2 s-1]; the block is inactive wherever that is zero, so this
  ! file is a no-op for any other run.
  !--------------------------------------------------------------------
  integer::iFh
  real(dp)::F_star, cE_1, F_1, scale_Np, scale_Fp

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  call rt_units(scale_Np, scale_Fp)

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over boundaries
  do ibound=1,nboundary

     call set_boundary_references(ibound,ind_ref,boundary_dir,inbor)

     ! Velocity sign switch for reflexive boundary conditions
     gs=(/1,1,1/)
     if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
     if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
     if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do

        ! Gather neighboring reference grid
        do i=1,ngrid
           ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather neighboring reference cell
           iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
           do i=1,ngrid
              ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
           end do

           ! Wall and free boundary conditions
           if((boundary_type(ibound)/10).ne.2)then

              ! Gather reference hydro variables
              do ivar=1,nrtvar
                 do i=1,ngrid
                    uu(i,ivar)=rtuold(ind_cell_ref(i),ivar)
                 end do
              end do
              ! Scatter to boundary region
              do ivar=1,nrtvar
                 switch=1
                 ! ------------------------------------------------------------
                 ! need to switch photon flux in RT...
                 ! rtVar= ivar = 1,2,3,4,.. = Np, Fx, Fy, Fz, Np, Fx,..
                 rtType = mod(ivar-1, ndim+1)
                 if(rtType .ne. 0) switch=gs(rtType) !0=Np, 1=Fx, 2=Fy, 3=Fz
                 ! ------------------------------------------------------------
                 do i=1,ngrid
                    rtuold(ind_cell(i),ivar)=uu(i,ivar)*switch
                 end do
              end do

              ! ---- RT15 eq. 83: emitting lower boundary -----------------
              ! uu still holds the RAW interior (reference cell) state, and
              ! this runs AFTER the reflexive switch, so the height component
              ! of the flux keeps its + sign. boundary_dir = 2*ndim-1 is the
              ! LOWER face of the last spatial dimension: x in 1D, y in 2D.
              if(rt_n_bound(ibound) .gt. 0d0 .and.                          &
                   boundary_dir .eq. 2*ndim-1) then
                 iFh    = 1+ndim
                 F_star = rt_n_bound(ibound)/scale_Fp
                 do i=1,ngrid
                    cE_1 = uu(i,1)*rt_c(ilevel)
                    F_1  = uu(i,iFh)
                    rtuold(ind_cell(i),1)     = (F_star-F_1+cE_1)/rt_c(ilevel)
                    rtuold(ind_cell(i),2:iFh) = 0d0
                    rtuold(ind_cell(i),iFh)   = F_star
                 end do
              endif
              ! -----------------------------------------------------------

              ! Imposed boundary conditions
           else

              ! Compute cell center in code units
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

              call rt_boundana(xx,uu,dx_loc,ibound,ngrid)

              ! Scatter variables
              do ivar=1,nrtvar
                 do i=1,ngrid
                    rtuold(ind_cell(i),ivar)=uu(i,ivar)
                 end do
              end do

           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering rt_make_boundary_hydro for level ',I2)

end subroutine rt_make_boundary_hydro
