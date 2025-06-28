subroutine init_pop3
  use amr_commons
  use pm_commons
  use clfind_commons
  use hydro_commons
  use amr_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::idim,ivar
  integer::i !,isink
  integer::ilun,nx_loc
  !integer::nsinkold
  !real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  logical,allocatable,dimension(:)::nb
  logical::eof,ic_sink=.false.
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  integer,parameter::tag=1357,tag2=3579
  integer::dummy_io,info2

  !allocate pop3 variables, Kang modified 210618
  allocate(xpop3(1:npop3max,1:ndim))
  allocate(xpop3_new(1:npop3max,1:ndim))
  allocate(xpop3_all(1:npop3max,1:ndim))

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(nrestart>0)then
     ilun=4*ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/pop3_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/pop3_'//TRIM(nchar)//'.out'
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
     rewind(ilun)
     read(ilun)npop3_tot

     if(npop3_tot>0)then
        allocate(xdp(1:npop3_tot))
        do idim=1,ndim
           read(ilun)xdp ! Read sink position
           xpop3(1:npop3_tot,idim)=xdp
        end do
     end if
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

  endif

end subroutine init_pop3
