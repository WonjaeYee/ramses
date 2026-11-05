module mpi_mod
#ifndef WITHOUTMPI
#ifdef MPI_OLD
  include 'mpif.h'
#else
  use mpi
#endif
#endif
end module mpi_mod
!################################################################
! The ranks that share this node's memory (MPI-3): node_comm from MPI_Comm_split_type, set in initialize_mpi.
! Used by the cooling queue (cool_queue = 2). Without MPI the node holds one rank.
!################################################################
module node_shared
  use mpi_mod
  use iso_c_binding
  implicit none
  integer::node_comm=0,node_rank=0,node_size=1
  ! Read-only tables in node-shared memory (shared_tables, &RUN_PARAMS): node rank 0 reads them, the other
  ! ranks of the node map the same memory. The allocations are collective on node_comm, in the same order.
  logical::shared_tables=.false.
  integer,parameter,private::max_wins=1024
  integer,private::nwins=0
  integer,dimension(max_wins),private::wins
contains
  subroutine node_shared_init
#ifndef WITHOUTMPI
    integer::ierr
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,node_comm,ierr)
    call MPI_COMM_RANK(node_comm,node_rank,ierr)
    call MPI_COMM_SIZE(node_comm,node_size,ierr)
#endif
  end subroutine node_shared_init
  !--------------------------------------------------------------
  logical function node_tables()
    ! The tables are shared on this run's nodes
    node_tables=shared_tables.and.node_size>1
  end function node_tables
  logical function node_reader()
    ! This rank reads (and fills) the tables
    node_reader=(.not.node_tables()).or.node_rank==0
  end function node_reader
  !--------------------------------------------------------------
  function node_shared_alloc(nbytes) result(cp)
    ! nbytes of memory shared by the ranks of the node, allocated on node rank 0 (collective on node_comm)
    integer(kind=8),intent(in)::nbytes
    type(c_ptr)::cp
#ifndef WITHOUTMPI
    integer(kind=MPI_ADDRESS_KIND)::wsize,qsize,baseptr
    integer::du,ierr
    if(nwins==max_wins)then
       write(*,*)'node_shared: more than',max_wins,'shared tables'
       call clean_stop
    end if
    wsize=0
    if(node_rank==0)wsize=max(nbytes,8_8)
    nwins=nwins+1
    call MPI_WIN_ALLOCATE_SHARED(wsize,1,MPI_INFO_NULL,node_comm,baseptr,wins(nwins),ierr)
    call MPI_WIN_SHARED_QUERY(wins(nwins),0,qsize,du,baseptr,ierr)
    cp=transfer(baseptr,cp)
    call MPI_WIN_LOCK_ALL(MPI_MODE_NOCHECK,wins(nwins),ierr)
#else
    cp=c_null_ptr
#endif
  end function node_shared_alloc
  !--------------------------------------------------------------
  subroutine node_alloc_r8_2(a,n1,n2)
    real(kind=8),pointer,contiguous,dimension(:,:)::a
    integer,intent(in)::n1,n2
    if(.not.node_tables())then
       allocate(a(n1,n2))
    else
       call c_f_pointer(node_shared_alloc(8_8*n1*n2),a,(/n1,n2/))
    end if
  end subroutine node_alloc_r8_2
  subroutine node_alloc_r8_3(a,n1,n2,n3)
    real(kind=8),pointer,contiguous,dimension(:,:,:)::a
    integer,intent(in)::n1,n2,n3
    if(.not.node_tables())then
       allocate(a(n1,n2,n3))
    else
       call c_f_pointer(node_shared_alloc(8_8*n1*n2*n3),a,(/n1,n2,n3/))
    end if
  end subroutine node_alloc_r8_3
  subroutine node_alloc_r8_4(a,n1,n2,n3,n4)
    real(kind=8),pointer,contiguous,dimension(:,:,:,:)::a
    integer,intent(in)::n1,n2,n3,n4
    if(.not.node_tables())then
       allocate(a(n1,n2,n3,n4))
    else
       call c_f_pointer(node_shared_alloc(8_8*n1*n2*n3*n4),a,(/n1,n2,n3,n4/))
    end if
  end subroutine node_alloc_r8_4
  subroutine node_alloc_r4_6(a,n)
    ! real(4) table of shape n(1:6) (the MIST sample)
    real(kind=4),pointer,contiguous,dimension(:,:,:,:,:,:)::a
    integer,dimension(6),intent(in)::n
    if(.not.node_tables())then
       allocate(a(n(1),n(2),n(3),n(4),n(5),n(6)))
    else
       call c_f_pointer(node_shared_alloc(4_8*product(int(n,8))),a,n)
    end if
  end subroutine node_alloc_r4_6
  !--------------------------------------------------------------
  subroutine node_publish
    ! After the reader has filled the shared tables: make them visible to the node (collective on node_comm)
#ifndef WITHOUTMPI
    integer::i,ierr
    if(.not.node_tables())return
    do i=1,nwins
       call MPI_WIN_SYNC(wins(i),ierr)
    end do
    call MPI_BARRIER(node_comm,ierr)
    do i=1,nwins
       call MPI_WIN_SYNC(wins(i),ierr)
    end do
#endif
  end subroutine node_publish
  !--------------------------------------------------------------
  subroutine node_bcast_r8(x,n)
    ! The reader's values to the node
    integer,intent(in)::n
    real(kind=8),dimension(n)::x
#ifndef WITHOUTMPI
    integer::ierr
    if(node_tables())call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,node_comm,ierr)
#endif
  end subroutine node_bcast_r8
  subroutine node_bcast_r4(x,n)
    integer,intent(in)::n
    real(kind=4),dimension(n)::x
#ifndef WITHOUTMPI
    integer::ierr
    if(node_tables())call MPI_BCAST(x,n,MPI_REAL,0,node_comm,ierr)
#endif
  end subroutine node_bcast_r4
  subroutine node_bcast_i(x,n)
    integer,intent(in)::n
    integer,dimension(n)::x
#ifndef WITHOUTMPI
    integer::ierr
    if(node_tables())call MPI_BCAST(x,n,MPI_INTEGER,0,node_comm,ierr)
#endif
  end subroutine node_bcast_i
end module node_shared
