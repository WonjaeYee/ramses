!================================================================
! tests/passive_scalars/accuracy: initial conditions from a table (MHD).
! ic_table.dat in the run directory:
!   line 1: nrow ncol            (ncol = nvar-1-nener: s rho vx P Bx By Bz q_1 ... q_npassive)
!   rows  : primitive values tabulated in s, increasing; the scalars as fractions of rho.
! Without ic_geom.dat the profile is a function of x (s = x) and vx comes from the table.
! With ic_geom.dat ("kx ky kz vx vy vz") it is evaluated at s = mod(kx*x+ky*y+kz*z, boxlen)
! and the velocity is the uniform (vx,vy,vz).  The table is linearly interpolated at the
! cell centre: with rows at every s a cell centre takes, the values are taken exactly (so the
! table can hold exact cell averages).  The non-thermal energies are zero.
!================================================================
subroutine acc_table_condinit(x,q,dx,nn)
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
  real(dp),save::kx=1.0d0, ky=0.0d0, kz=0.0d0, vgx=0.0d0, vgy=0.0d0, vgz=0.0d0
  logical,save::geom=.false.
  real(dp)::s
  logical::have_geom

  if(nrow==0)then
     open(unit=97,file='ic_table.dat',status='old',action='read')
     read(97,*)nrow,ncol
     if(ncol/=nvar-1-nener)then
        write(*,*)'acc_table_condinit: ic_table.dat has',ncol,' columns, expected',nvar-1-nener
        call clean_stop
     end if
     allocate(tab(1:nrow,1:ncol))
     do j=1,nrow
        read(97,*)tab(j,1:ncol)
     end do
     close(97)
     ! optional ic_geom.dat: "kx ky kz vx vy vz". The profile is evaluated at s = mod(kx*x+ky*y+kz*z, boxlen)
     ! (integer k keep it periodic in the box) and the velocity is the uniform (vx, vy, vz).
     inquire(file='ic_geom.dat',exist=have_geom)
     if(have_geom)then
        open(unit=98,file='ic_geom.dat',status='old',action='read')
        read(98,*)kx,ky,kz,vgx,vgy,vgz
        close(98)
        geom=.true.
     end if
  end if

  q(1:nn,1:nvar+3)=0.0d0
  do i=1,nn
     s=x(i,1)
     if(geom)then
        s=kx*x(i,1)
#if NDIM>1
        s=s+ky*x(i,2)
#endif
#if NDIM>2
        s=s+kz*x(i,3)
#endif
        s=modulo(s,boxlen)
     end if
     if(s<=tab(1,1))then
        ilo=1; ihi=1; w=0.0d0
     else if(s>=tab(nrow,1))then
        ilo=nrow; ihi=nrow; w=0.0d0
     else
        ilo=1; ihi=nrow
        do while(ihi-ilo>1)
           imid=(ilo+ihi)/2
           if(tab(imid,1)<=s)then
              ilo=imid
           else
              ihi=imid
           end if
        end do
        w=(s-tab(ilo,1))/(tab(ihi,1)-tab(ilo,1))
     end if
     q(i,1)=(1.0d0-w)*tab(ilo,2)+w*tab(ihi,2)
     q(i,2)=(1.0d0-w)*tab(ilo,3)+w*tab(ihi,3)
     if(geom)then
        q(i,2)=vgx
        q(i,3)=vgy
        q(i,4)=vgz
     end if
     q(i,5)=(1.0d0-w)*tab(ilo,4)+w*tab(ihi,4)
     q(i,6:8)=(1.0d0-w)*tab(ilo,5:7)+w*tab(ihi,5:7)
     q(i,nvar+1:nvar+3)=q(i,6:8)
     do k=9+nener,nvar
        q(i,k)=(1.0d0-w)*tab(ilo,k-1-nener)+w*tab(ihi,k-1-nener)
     end do
  end do

end subroutine acc_table_condinit
