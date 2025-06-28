subroutine backup_pop3(filename)
  use amr_commons
  use pm_commons
  use amr_parameters
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif 

  character(LEN=80)::filename

  integer::ilun,idim,i,ivar
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii
  logical,allocatable,dimension(:)::nb
  integer,parameter::tag=5134
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering backup_pop3'

  ilun=4*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)

  ! Kang modified 210615
  write(ilun)npop3_tot
  if(npop3_tot>0)then
     allocate(xdp(1:npop3_tot))
     do idim=1,ndim
        do i=1,npop3_tot
           xdp(i)=xpop3(i,idim)
        end do
        write(ilun)xdp ! Write pop3 position
     enddo
     deallocate(xdp)
  endif
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
  
end subroutine backup_pop3
