!================================================================
! tests/passive_scalars, SOLVER=mhd: table_condinit for the MHD layout.
! Reads the same ic_table.dat as the hydro version (x rho vx P q(nhydro+1:nvar), ncol =
! nvar-4), with B = 0 (pure hydrodynamics through the MHD solver).  build.sh and the
! nightly suite append it to the tree's mhd/condinit.f90 with case('table').
!================================================================
subroutine table_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar+3)::q
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),allocatable,dimension(:,:),save::tab
  integer,save::nrow=0
  integer::i,j,k,ncol,ilo,ihi,imid
  real(dp)::w

  if(nrow==0)then
     open(unit=97,file='ic_table.dat',status='old',action='read')
     read(97,*)nrow,ncol
     if(ncol/=nvar-nhydro+4)then
        write(*,*)'table_condinit: ic_table.dat has',ncol,' columns, expected',nvar-nhydro+4
        call clean_stop
     end if
     allocate(tab(1:nrow,1:ncol))
     do j=1,nrow
        read(97,*)tab(j,1:ncol)
     end do
     close(97)
  end if

  q(1:nn,1:nvar+3)=0.0d0
  do i=1,nn
     ! bracketing rows (the x column is increasing)
     if(x(i,1)<=tab(1,1))then
        ilo=1; ihi=1; w=0.0d0
     else if(x(i,1)>=tab(nrow,1))then
        ilo=nrow; ihi=nrow; w=0.0d0
     else
        ilo=1; ihi=nrow
        do while(ihi-ilo>1)
           imid=(ilo+ihi)/2
           if(tab(imid,1)<=x(i,1))then
              ilo=imid
           else
              ihi=imid
           end if
        end do
        w=(x(i,1)-tab(ilo,1))/(tab(ihi,1)-tab(ilo,1))
     end if
     q(i,1)=(1.0d0-w)*tab(ilo,2)+w*tab(ihi,2)
     q(i,2)=(1.0d0-w)*tab(ilo,3)+w*tab(ihi,3)
     q(i,5)=(1.0d0-w)*tab(ilo,4)+w*tab(ihi,4)
     do k=nhydro+1,nvar
        q(i,k)=(1.0d0-w)*tab(ilo,k-nhydro+4)+w*tab(ihi,k-nhydro+4)
     end do
  end do

end subroutine table_condinit
