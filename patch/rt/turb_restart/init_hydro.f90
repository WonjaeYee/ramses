subroutine init_hydro
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_parameters,only:convert_birth_times,rt_src_x_center,rt_src_y_center,rt_src_z_center
#endif
  use mpi_mod
  implicit none
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar,irad
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1108
  integer::dummy_io,info2
  ! Harley additions to get the cell center
  real(dp)::dx,dx_loc,scale
  integer::nx_loc,ix,iy,iz,idim
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::loc_pos,tmp_pos

  if(verbose)write(*,*)'Entering init_hydro'

  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0
  if (MC_tracer) then
    allocate(fluxes(1:ncell, 1:twondim))
    fluxes(1:ncell, 1:twondim) = 0
  end if
#if NVARNOADVECT>0
  allocate(unoadvect(1:ncell,1:nvarnoadvect))
  !TODO(code): we might want a different initial value per variable
  unoadvect=noadvect_init
#endif
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     endif



     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2  !(ncpu)
     read(ilun)nvar2  !(nvar)
     ! Correct for non-advected quantities
     nvar2 = nvar2 - nvarnoadvect_og
     read(ilun)ndim2 !(ndim)
     read(ilun)nlevelmax2 !(nlevelmax)
     read(ilun)nboundary2 !(nboundary)
     read(ilun)gamma2     !(gamma)
     if(.not.(neq_chem.or.rt) .and. nvar2.ne.nvar)then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
        call clean_stop
     end if
#ifdef RT
     if((neq_chem.or.rt).and.nvar2.lt.nvar)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        if(myid==1) write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to that specified in var_region'
     end if
#endif
     do ilevel=1,nlevelmax2

        !!!! HARLEY ADDITIONS FOR POSITION
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

        !!!! END HARLEY ADDITIONS FOR POSITION

        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax

                 ! Read density and velocities --> density and momenta
                 do ivar=1,ndim+1
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                          !TODO(code): set the density to a very low value
                          !outside the sphere...
                          ! Get cell position in user units relative to center
                          loc_pos = 0.d0
                          do idim=1,ndim
                             tmp_pos = scale * (xg(ind_grid(i),idim) + xc(ind,idim) - skip_loc(idim))
                             if (idim .eq. 1) tmp_pos = (tmp_pos - rt_src_x_center(1))**2.d0
                             if (idim .eq. 2) tmp_pos = (tmp_pos - rt_src_y_center(1))**2.d0
                             if (idim .eq. 3) tmp_pos = (tmp_pos - rt_src_z_center(1))**2.d0
                             loc_pos = loc_pos + tmp_pos
                          end do
                          loc_pos = SQRT(loc_pos)
                          ! Hard coded minimum density for now
                          if (loc_pos.gt.(boxlen/2.d0)) uold(ind_grid(i)+iskip,1) = 1.d-6 
                          ! End get cell position
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 end do

#if NENER>0
                 ! Dont read non-thermal pressures --> non-thermal energies are
                 ! set to 0.0
                 do ivar=ndim+3,ndim+2+nener
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=0.d0
                    end do
                 end do
#endif
                 ! Read thermal pressure --> total fluid energy
                 read(ilun)xx
                 do i=1,ncache
                    xx(i)=xx(i)/(gamma-1d0)
                    if (uold(ind_grid(i)+iskip,1)>0.)then
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NENER>0
                    do irad=1,nener
                       xx(i)=xx(i)+uold(ind_grid(i)+iskip,ndim+2+irad)
                    end do
#endif
                 else
                    xx(i)=0.
                 end if
                    uold(ind_grid(i)+iskip,ndim+2)=xx(i)
                 end do
#if NVAR>NDIM+2+NENER
                 ! Set the passive scalars to their proper values
                 do ivar=ndim+3+nener,nvar
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=var_region(1,ivar-ndim-2-nener)*max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                 end do
#endif
#if NVARNOADVECT>0
                 ! Don't read the non-advected scalars...set them to 0
                 if (nvarnoadvect_og.gt.0) then
                    do ivar=1,nvarnoadvect_og
                       do i=1,ncache
                          unoadvect(ind_grid(i)+iskip,ivar)=0.d0
                       end do
                    end do
                 end if
#endif
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif



#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro
