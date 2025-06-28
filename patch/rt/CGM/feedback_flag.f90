subroutine feedback_flag(ilevel)
  use amr_commons
  use pm_commons
  use mechanical_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel.
  ! -------------------------------------------------------------------
  ! This routine ensures that refinement is added in which
  ! SN can potentially explode - this is done by looking at the sign 
  ! of the ID of a star particle.
  ! Activated only if feedback_refine=.true.
  ! nshell_resolve is the minimum number of cells used to resolve
  ! the shell formation radius: default:0 (no mpi is required)
  ! Nov 2014, Taysun Kimm
  ! -------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::icpu,npart1,idim,ip,igrid2,icell2,ilevel2
  integer::ind,ind_son,ind_cell,iskip,info
  integer,dimension(1:nvector),save::ind_grid,ind_pos_cell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc,x0(1:3),x3(1:3)
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  real(dp)::ttsta,ttend,current_time,t0
  real(dp)::msun2g=1.989d33
  logical::ok,oksub(1:twotondim)

  if(nstar_tot==0)return

! important to initialize
  ncomm_SN=0

#ifndef WITHOUTMPI 
  if(myid.eq.1) ttsta=MPI_WTIME(info)
  xSN_comm=0d0;iSN_comm=0
  uidSN_comm=0 !unique id of SN in this cpu
#endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = (boxlen*scale_l)**3*scale_d/msun2g

  ! Mesh spacing in that level
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)

  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_delay*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
  else 
     t0=t_delay*1d6*(365.*24.*3600.)/scale_t
     current_time=t
  endif

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ip=0
      
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid) ! Number of particles in the grid

        oksub=.false.
              
        ! Count star particles
        if (npart1>0)then
           do idim=1,ndim
              x0(idim)=xg(igrid,idim)-dx-skip_loc(idim)
           end do
           ipart=headp(igrid)

           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle <--- Very important !!
              next_part=nextp(ipart)
              ok=.false.

              !if(idp(ipart).le.0) ok=.true. ! don't do this-> outflow structures de-refined too quickly
              if(tp(ipart).ge.(current_time-t0).and.abs(tp(ipart)).gt.1d-10)then ! dark matter: tp=0
                 ok=.true.
              endif

              if(ok)then
                 ind_son=1
                 do idim=1,ndim
                    ind = int((xp(ipart,idim)/scale-x0(idim))/dx)
                    ind_son=ind_son+ind*2**(idim-1)
                 end do
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 ind_cell=iskip+igrid

                 if(son(ind_cell)==0)then !leaf cell
                    oksub(ind_son)=.true.
                 endif
              endif !ok
        
              ipart=next_part ! Go to next particle
           enddo ! End loop over particles

        end if ! npart1>0
 
        do ind=1,twotondim
           if (oksub(ind)) then
              ip=ip+1
              ind_grid(ip)=igrid
              ind_pos_cell(ip)=ind
           
              if (ip==nvector)then
                 call feedback_refine_fine(ind_grid,ind_pos_cell,ip,ilevel)
                 ip=0
              endif
           end if
        end do
        igrid=next(igrid) ! Go to next grid
     end do ! End loop over grids

     if (ip>0) then
        call feedback_refine_fine(ind_grid,ind_pos_cell,ip,ilevel)
        ip=0
     endif
 
  end do ! End loop over cpu

! This is not to de-refine when a star is not seen in ilevel
! This happens because numbp keeps track of star particles in the leaf cell

  ! Mesh spacing in that level
  call mesh_info (ilevel+1,skip_loc,scale,dx,dx_loc,vol_loc,xc)

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel+1)
     ip=0
      

     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel+1)
        npart1=numbp(igrid) ! Number of particles in the grid
 
        ! Count star particles
        ok=.false.
        if (npart1>0)then
           do idim=1,ndim
              x0(idim)=xg(igrid,idim)-dx-skip_loc(idim)
           end do

           ipart=headp(igrid)

           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle <--- Very important !!
              next_part=nextp(ipart)

              !if(idp(ipart).le.0) ok=.true. ! don't do this-> outflow structures de-refined too quickly
              if(tp(ipart).ge.(current_time-t0).and.abs(tp(ipart)).gt.1d-10)then
                 ok=.true.
              endif

              ipart=next_part ! Go to next particle
           enddo ! End loop over particles

        end if ! npart1>0

        if(ok)then
          call get_icell_from_pos (x0,ilevel,igrid2,icell2,ilevel2) 
          ip=ip+1
          ind_grid(ip)=igrid2
          ind=(icell2-igrid2-ncoarse)/ngridmax+1
          !if (verbose_ref) print *, '#### ind =' ,igrid2,icell2,ilevel2,ind,ilevel,nlevelmax
          ind_pos_cell(ip)=ind
          if (ip==nvector)then
              call feedback_refine_fine(ind_grid,ind_pos_cell,ip,ilevel)
              ip=0
           endif
        end if ! if ok
        igrid=next(igrid) ! Go to next grid
     end do ! End loop over grids

     if (ip>0) then
        call feedback_refine_fine(ind_grid,ind_pos_cell,ip,ilevel)
        ip=0
     endif
 
  end do ! End loop over cpu

#ifndef WITHOUTMPI 
  call feedback_refine_fine_mpi(ilevel)
#endif

end subroutine feedback_flag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine feedback_refine_fine(ind_grid,ind_pos_cell,ncell,ilevel)
   use amr_commons
   use pm_commons
   use hydro_commons
   use mechanical_commons
   use cooling_module, ONLY:mH
   implicit none
   integer::ncell,ilevel
   integer,dimension(1:nvector)::ind_grid,ind_pos_cell
   !--------------------------------------------------
   ! This routine is called by feedback_flag
   !--------------------------------------------------
   integer::i,j,k,icell,igrid,idim,ilevel2,ksta,kend,kk
   real(dp)::dx,dx_loc,scale,vol_loc
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale_msun
   real(dp)::skip_loc(1:3),xc3(1:3)
   real(dp),dimension(1:twotondim,1:ndim),save::xc
   ! Grid based arrays
   real(dp),dimension(1:ndim,1:nvector),save::xc2
   real(dp)::ESN_code,nH,mgas,Zgas,mass_tr,chi_tr,msun2g
   integer::ncomm,ista,iend
 
   msun2g=1.989d33
   ! Conversion factor from user units to cgs units 
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   scale_msun = scale_l**3*scale_d/msun2g

   ! Mesh variables
   call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  
   ESN_code = nsn_resolve*E_SNII/(M_SNII*msun2g)/scale_v**2
   
   ! Record position of each cell containing SNe [0.0-1.0]
   xc2=0d0
   do i=1,ncell
      do idim=1,ndim
         xc2(idim,i)=xg(ind_grid(i),idim)-skip_loc(idim)+xc(ind_pos_cell(i),idim)
      end do
   end do

   ! Loop over SN cells
   ncomm_SN=0 ! total number of cells for which communications needed
   
   do i=1,ncell

      mrefnei=0d0;mzrefnei=0d0;icellnei=0;lrefnei=0
      mzrefnei_ring=0d0;mrefnei_ring=0d0

      ! Loop over a ring (start from the SNcell to rings with larger radii)
      j=0; icommr=0
      do while (j<=nshell_resolve)

         kend = nrefnei_ring(j)
         if (j==0)then
            ksta = 1
         else
            ksta = nrefnei_ring(j-1)+1
         endif 
 
         ncomm=0 ! this should be placed here, as the j=nshell_resolve ring
                 ! would cover the largest cpu domains

         ! Loop over cells within a ring
         do k=ksta,kend

            kk=irefnei(k)
            xc3(1:3) = xc2(1:3,i) + xrefnei(1:3,kk)*dx_loc

            call get_icell_from_pos (xc3,ilevel,igrid,icell,ilevel2)

            if(cpu_map(father(igrid))==myid)then
               mrefnei(k)=uold(icell,1)*vol_loc*scale_msun
               if(metal) then
                  !mzrefnei(k)=uold(icell,imetal)*dx_loc**3*scale_msun ! uold(:,imetal)=rho*Z
                  mzrefnei(k)=(2.09d0*uold(icell,imetal+1)+1.06d0*uold(icell,imetal))*dx_loc**3*scale_msun ! metal enrichment
               endif
               icellnei(k)=icell
               lrefnei(k)=ilevel2
            else ! need MPI
               ncomm=ncomm+1
               icommr(ncomm)=cpu_map(father(igrid)) ! cpu id to which info should be passed
            endif
 
         end do ! k: cells within a ring

         ! if mpi is needed, the following condition will ensure to check all the cells within nshell_resolve
         if (ncomm==0)then 

            ! MPI should be done by this line
            mrefnei_ring(j)=sum(mrefnei(ksta:kend))
            mzrefnei_ring(j)=sum(mzrefnei(ksta:kend))
    
            ! test if the local cell contains already more than enough gas mass
            nH = sum(mrefnei_ring(0:j))*msun2g*0.76/mH/(dble(kend)*(dx_loc*scale_l)**3)
            mgas = sum(mrefnei_ring(0:j))
            Zgas = sum(mzrefnei_ring(0:j))/mgas
            Zgas = max(Zgas/0.02,0.01) ! solar scaled
            if(.not.metal) Zgas=max(Zgas,z_ave) 
   
            ! transition mass from the adiabatic to snowplow phase (Kimm & Cen 2014, A5)
            chi_tr = (A_SN/1d4)**2/f_esn/M_SNII*nsn_resolve**((expE_SN-1d0)*2d0)*nH**(expN_SN*2d0)*Zgas**(expZ_SN*2d0)
            mass_tr = M_SNII*(chi_tr-1d0)/reduce_mass_tr ! Msun
    
            ! TODO: do I need to worry that level-1 grids are updated here?
            ! check if the total mass is greater than the transition mass within this ring
            ! 1) mass>mass_tr before j < nshell_resolve (refine it. simple.) 
            ! 2) mass>mass_tr at j=nshell_resolve (refine only if lrefnei < ilevel)
            ! 3) mass<mass_tr (refine only if lrefnei < ilevel)
            ! NB. 2) and 3) are the same in practice

            if(mgas.ge.mass_tr)then
               if(j<nshell_resolve)then
                  do k=ksta,kend
                     icell = icellnei(k)
                     if(flag1(icell)==0)then
                        nflag=nflag+1
                        flag1(icell)=1
                     end if
                  end do ! k
                  j=nshell_resolve+1 !good as long as j>nshell_resolve
               endif

               if(nshell_resolve==0)then
                  icell = icellnei(1)
                  if(flag1(icell)==0)then
                     nflag=nflag+1
                     flag1(icell)=1
                  end if
               endif
            endif

            if(j==nshell_resolve)then
               do k=ksta,kend
                  icell = icellnei(k)
                  if(flag1(icell)==0.and.lrefnei(k)<ilevel)then
                     nflag=nflag+1
                     flag1(icell)=1
                  end if
               end do ! k
               j=nshell_resolve+1 !good as long as j>nshell_resolve
            endif

         endif 

         j=j+1
      end do ! j: loop over a ring 

      if(ncomm>0)then ! if mpi is needed, tidy up the variable and bookkeep
         if(ncomm>1)then
            call redundant_non_1d(icommr(1:ncomm),ncomm,ncomm)
         endif
         
         ista=ncomm_SN+1
         iend=ista+ncomm-1
 
         if(iend>ncomm_max)then
            write(*,*) 'Fatal error in feedback_refine_fine: increase ncomm_max', ncomm_max,iend
            call clean_stop
         endif
         ! Let me recycle the variable for mechanical SN feedback to reduce the memory usage
         !if (verbose_ref) print *,ista,iend,ncomm,ncomm_SN

         uidSN_comm = uidSN_comm+1
         iSN_comm   (ista:iend) = icommr(1:ncomm)
         idSN_comm  (ista:iend) = uidSN_comm
         xSN_comm (1,ista:iend) = xc2(1,i)
         xSN_comm (2,ista:iend) = xc2(2,i)
         xSN_comm (3,ista:iend) = xc2(3,i)
         ncomm_SN=ncomm_SN+ncomm ! ncomm_SN: the total number of comm.

      endif

      !if (verbose_ref) print *, 'mg=',int(ilevel,kind=1),sngl(mrefnei_ring), &
      !               & int(ncomm,kind=2),int(nflag,kind=2)

  end do ! i: loop over SNe

end subroutine feedback_refine_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifndef WITHOUTMPI
subroutine feedback_refine_fine_mpi(ilevel)
   use amr_commons
   use pm_commons
   use hydro_commons
   use mechanical_commons
   use cooling_module, ONLY:mH
   implicit none
   include 'mpif.h'
   integer::ilevel
   !--------------------------------------------------
   ! This routine is called by feedback_flag
   !--------------------------------------------------
   real(dp),dimension(:,:),allocatable::SNsend,SNrecv
   real(dp),dimension(:,:),allocatable::SNmetasend,SNmetarecv
   real(dp),dimension(:,:),allocatable::mring_other,mzring_other
   real(dp),dimension(:,:),allocatable::mring_mysne,mzring_mysne
   integer, dimension(:)  ,allocatable::lring_mysne
   integer ,dimension(:)  ,allocatable::list2recv,list2send
   integer ,dimension(:)  ,allocatable::reqrecv,reqsend
   integer ,dimension(:,:),allocatable::statrecv,statsend
   integer, dimension(:)  ,allocatable::lsend,lrecv
   !--------------------------------------------------
   integer::i,j,k,icell,igrid,ilevel2,ksta,kend,ii,kk
   real(dp)::dx,dx_loc,scale,vol_loc
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale_msun
   real(dp)::skip_loc(1:3),xc3(1:3),xSN_i(1:3)
   real(dp),dimension(1:twotondim,1:ndim),save::xc
   real(dp)::ESN_code
   integer::ncomm,info,icpu,tag,ista,iend
   integer::ncpu_send,ncpu_recv,ncomm_recv,ncomm_send
   integer::ncomm_tot,ncc,cpu2recv,cpu2send
   integer::irecv_sta,irecv_end,isend_sta,nvarcom
   real(dp)::nH,mgas,Zgas,mass_tr,chi_tr,msun2g

   msun2g=1.989d33

   ! Conversion factor from user units to cgs units 
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   scale_msun = scale_l**3*scale_d/msun2g

   ! Mesh variables
   call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  
   ESN_code = nsn_resolve*E_SNII/(M_SNII*msun2g)/scale_v**2

   !===========================================================
   ! Step I: send info to other cpus and receive
   !===========================================================
   ncpu_send=0;ncpu_recv=0
   ncomm_SN_cpu=0
   ncomm_SN_mpi=0
   ncomm_SN_mpi(myid)=ncomm_SN
   call MPI_ALLREDUCE(ncomm_SN_mpi,ncomm_SN_cpu,ncpu,&
                    & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
   ! the number of communications needed
   ncomm_tot=sum(ncomm_SN_cpu) 
   if(ncomm_tot==0) return

   allocate(icpuSN_comm    (1:ncomm_tot,1:2))
   allocate(icpuSN_comm_mpi(1:ncomm_tot,1:2))
 
   ! index for mpi variable
   if(myid==1)then
      isend_sta=0
   else
      isend_sta=sum(ncomm_SN_cpu(1:myid-1))
   endif

   icpuSN_comm=0
   do i=1,ncomm_SN_cpu(myid)
      icpuSN_comm(isend_sta+i,1) = myid        ! sender
      icpuSN_comm(isend_sta+i,2) = iSN_comm(i) ! receiver
      ! iSN_comm: local variable
      ! icpuSN_comm: local (but extended) variable
      ! icpuSN_comm_mpi: mpi variable
   end do

   ! share the list of communications
   icpuSN_comm_mpi=0
   call MPI_ALLREDUCE (icpuSN_comm,icpuSN_comm_mpi,ncomm_tot*2,&
                     & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

   ncomm_send = ncomm_SN_cpu(myid) ! equals to ncomm_SN
   ncomm_recv = count(icpuSN_comm_mpi(:,2).eq.myid,1)

   nvarcom=3   
   !check if myid needs to send anything
   if(ncomm_send>0)then
      allocate(SNsend (1:nvarcom,1:ncomm_send)) ! x(3)
      allocate(list2send (1:ncomm_send))
      list2send=0
      list2send=icpuSN_comm_mpi(isend_sta+1:isend_sta+ncomm_send,2)
      ncpu_send=1
      if(ncomm_send>1) call redundant_non_1d (list2send,ncomm_send,ncpu_send)
      ! ncpu_send = No. of cpus to which myid should send info
      allocate ( reqsend (1:ncpu_send))
      allocate ( statsend(1:MPI_STATUS_SIZE,1:ncpu_send))
      reqsend=0;statsend=0
   endif

   ! check if myid needs to receive anything
   if(ncomm_recv>0)then
      allocate( SNrecv (1:nvarcom,1:ncomm_recv))
      allocate( list2recv(1:ncomm_recv) )
      list2recv=0
      j=0
      do i=1,ncomm_tot
         if(icpuSN_comm_mpi(i,2)==myid)then
            j=j+1
            list2recv(j)=icpuSN_comm_mpi(i,1)
         endif
      enddo
  
      ncc = j
      if(ncc.ne.ncomm_recv)then !sanity check
         write(*,*) 'Fatal error in feedback_refine_fine_mpi: ncc!=ncomm_recv',ncc,ncomm_recv,myid
         call clean_stop
      endif

      ncpu_recv=1 ! No. of cpus from which myid should receive info
      if(j>1) call redundant_non_1d(list2recv,ncc,ncpu_recv)
      allocate( reqrecv  (1:ncpu_recv))
      allocate( statrecv (1:MPI_STATUS_SIZE,1:ncpu_recv))
      reqrecv=0; statrecv=0
   endif

   ! prepare one variable and send
   if (ncomm_send>0)then
      do icpu=1,ncpu_send
         SNsend=0d0
         cpu2send = list2send(icpu)
         ncc=0 ! No. of comm cells for cpu2send
         do i=1,ncomm_send
            j=i+isend_sta
            if(icpuSN_comm_mpi(j,2)==cpu2send)then
               ncc=ncc+1
               SNsend(1,ncc)=xSN_comm(1,i)
               SNsend(2,ncc)=xSN_comm(2,i)
               SNsend(3,ncc)=xSN_comm(3,i)
            end if
         end do ! i=1,ncomm_send

         tag = myid + cpu2send + ncc+10000
         call MPI_ISEND (SNsend(1:nvarcom,1:ncc),ncc*nvarcom,MPI_DOUBLE_PRECISION, &
                       & cpu2send-1,tag,MPI_COMM_WORLD,reqsend(icpu),info)
      end do ! icpu
   endif ! ncomm_send>0


   ! receive one large varaible
   if(ncomm_recv>0)then
      irecv_sta=1
      SNrecv=0d0
      do icpu=1,ncpu_recv
         cpu2recv = list2recv(icpu)
         ncc=0 ! No. of comm cells for cpu2recv
         do i=1,ncomm_tot
            if( (icpuSN_comm_mpi(i,1)==cpu2recv.and.&
              &  icpuSN_comm_mpi(i,2)==myid)   ) then
               ncc=ncc+1
            endif
         end do !i
         irecv_end=irecv_sta+ncc-1
         tag = myid + cpu2recv + ncc+10000
         call MPI_IRECV (SNrecv(1:nvarcom,irecv_sta:irecv_end),ncc*nvarcom,MPI_DOUBLE_PRECISION, &
                      &  cpu2recv-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)
         irecv_sta=irecv_end+1
      end do !icpu
   endif ! ncomm_recv>0


   if(ncpu_send>0) call MPI_WAITALL(ncpu_send,reqsend,statsend,info)
   if(ncpu_recv>0) call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)

   !===========================================================
   ! Step II: calculate the mass within each ring for SNe that do not belong to this cpu
   !===========================================================
   if(ncomm_recv>0)then
      allocate(mring_other (0:nshell_resolve,1:ncomm_recv))
      allocate(mzring_other(0:nshell_resolve,1:ncomm_recv))
      mring_other=0d0;mzring_other=0d0
   endif

   do i=1,ncomm_recv
      xSN_i(1:3) = SNrecv(1:3,i)

      ! Loop over a ring (start from the SNcell to rings with a larger radius)
      j=0; icommr=0; mrefnei=0d0; mzrefnei=0d0
      do while (j<=nshell_resolve)

         ksta = 1
         if(j>0) ksta = nrefnei_ring(j-1)+1
         kend = nrefnei_ring(j)
 
         ncomm = 0
         ! Loop over cells within a ring
         do k=ksta,kend

            kk=irefnei(k)
            xc3(1:3) = xSN_i(1:3) + xrefnei(1:3,kk)*dx_loc
         
            call get_icell_from_pos (xc3,ilevel,igrid,icell,ilevel2)

            if(cpu_map(father(igrid))==myid)then
               mrefnei(k)=uold(icell,1)*dx_loc**3*scale_msun
               if(metal) then
                  !mzrefnei(k)=uold(icell,imetal)*dx_loc**3*scale_msun ! uold(:,imetal)=rho*Z
                  mzrefnei(k)=(2.09d0*uold(icell,imetal+1)+1.06d0*uold(icell,imetal))*dx_loc**3*scale_msun ! metal enrichment
               end if
            endif
 
         end do ! k: cells within a ring

         mring_other(j,i)=sum(mrefnei(ksta:kend))
         mzring_other(j,i)=sum(mzrefnei(ksta:kend))
   
         j=j+1 
      end do ! j: loop over a ring 
      !print *, '=====>',myid, i,sngl(mring_other(:,i)) 
   end do ! i: loop over SNe

   !===========================================================
   ! Step III: send the information back to the original cpus
   !===========================================================
   nvarcom=(nshell_resolve+1)*2  !mass and metals in each ring
   if(ncomm_recv>0)then
      reqrecv=0    ! confusing, but need to do this way, 
      statrecv=0   ! since reqrecv is allocated when ncomm_recv>0
      allocate(SNmetasend(1:nvarcom,1:ncomm_recv))
 
      ! send one large varaible
      irecv_sta=1 ! notice this is not isend_sta (double-checked)
      do icpu=1,ncpu_recv
         SNmetasend=0d0
         cpu2send = list2recv(icpu)
         ncc=0 ! No. of comm cells for cpu2recv
         do i=1,ncomm_tot
            if( (icpuSN_comm_mpi(i,1)==cpu2send.and.&
              &  icpuSN_comm_mpi(i,2)==myid)   ) then
               ncc=ncc+1
               ista=1
               iend=ista+nshell_resolve
               SNmetasend(ista:iend,ncc)=mring_other(0:nshell_resolve,ncc)
               ista=iend+1
               iend=ista+nshell_resolve
               SNmetasend(ista:iend,ncc)=mzring_other(0:nshell_resolve,ncc)
            endif
         end do !i
         tag = myid + cpu2send + ncc+20000
         call MPI_ISEND (SNmetasend(1:nvarcom,1:ncc),ncc*nvarcom,MPI_DOUBLE_PRECISION, &
                      &  cpu2send-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)
      end do !icpu
   endif ! ncomm_recv>0

   if(ncomm_send>0)then 
      reqsend=0   ! confusing, but need to do this way, 
      statsend=0  ! since reqsend is allocated when ncomm_send>0
      allocate(SNmetarecv(1:nvarcom,1:ncomm_send))
      SNmetarecv=0d0 

      ! receive chunks into one large variable
      irecv_sta=1
      do icpu=1,ncpu_send
         cpu2recv = list2send(icpu)
         ncc=0 ! No. of comm cells for cpu2send
         do i=1,ncomm_tot
            if( icpuSN_comm_mpi(i,1)==myid.and.&
              & icpuSN_comm_mpi(i,2)==cpu2recv)then
               ncc=ncc+1
            end if
         end do ! i=1,ncomm_send
         irecv_end=irecv_sta+ncc-1

         tag = myid + cpu2recv + ncc + 20000
         call MPI_IRECV (SNmetarecv(1:nvarcom,irecv_sta:irecv_end),ncc*nvarcom,MPI_DOUBLE_PRECISION, &
                       & cpu2recv-1,tag,MPI_COMM_WORLD,reqsend(icpu),info)
         irecv_sta=irecv_end+1
      end do ! icpu
   endif ! ncomm_send>0

   if(ncpu_recv>0) call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)
   if(ncpu_send>0) call MPI_WAITALL(ncpu_send,reqsend,statsend,info)

   !===========================================================
   ! Step IV: determine up to what radius we should refine 
   !          - only for SNe that belong to this cpu
   !===========================================================
   if (ncomm_send>0)then

      allocate(mring_mysne (0:nshell_resolve,1:uidSN_comm))
      allocate(mzring_mysne(0:nshell_resolve,1:uidSN_comm))
      allocate(lring_mysne (                 1:uidSN_comm))
      mring_mysne=0d0;mzring_mysne=0d0;lring_mysne=-1


      do i=1,uidSN_comm

         !mring for the cells belonging to myid
         ! this is to find the position of this SN 
         ii=1
         do while (ii<=ncomm_send)
            j=ii+isend_sta
            if(icpuSN_comm_mpi(j,2)==cpu2send)then
               xSN_i(1:3) = xSN_comm(1:3,ii)
               ii=ncomm_send+1 ! escape
            endif
            ii=ii+1
         end do

         ! Loop over a ring (start from the SNcell to rings with a larger radius)
         j=0; icommr=0; mrefnei=0d0; mzrefnei=0d0
         do while (j<=nshell_resolve)
         
            ksta = 1
            if(j>0) ksta = nrefnei_ring(j-1)+1
            kend = nrefnei_ring(j)
         
            ncomm = 0
            ! Loop over cells within a ring
            do k=ksta,kend
         
               kk=irefnei(k)
               xc3(1:3) = xSN_i(1:3) + xrefnei(1:3,kk)*dx_loc
            
               call get_icell_from_pos (xc3,ilevel,igrid,icell,ilevel2)
         
               if(cpu_map(father(igrid))==myid)then
                  mrefnei(k)=uold(icell,1)*dx_loc**3*scale_msun
                  if(metal) then 
                     !mzrefnei(k)=uold(icell,imetal)*dx_loc**3*scale_msun ! uold(:,imetal)=rho*Z
                     mzrefnei(k)=(2.09d0*uold(icell,imetal+1)+1.06d0*uold(icell,imetal))*dx_loc**3*scale_msun ! metal enrichment
                  end if
               endif
         
            end do ! k: cells within a ring
         
            mrefnei_ring(j)=sum(mrefnei(ksta:kend))
            mzrefnei_ring(j)=sum(mzrefnei(ksta:kend))
         
            j=j+1 
         end do ! j: loop over a ring 

         mring_mysne(:,i)=mring_mysne(:,i)+mrefnei_ring(:)
         mzring_mysne(:,i)=mzring_mysne(:,i)+mzrefnei_ring(:)

      end do ! i=1,uidSN_comm

      irecv_sta=0 ! index of comm cells for all cpus
      do icpu=1,ncpu_send
         cpu2send = list2send(icpu)
         do ii=1,ncomm_send
            j=ii+isend_sta
            if(icpuSN_comm_mpi(j,2)==cpu2send)then
               irecv_sta=irecv_sta+1
               ista=1
               iend=ista+nshell_resolve
               k = idSN_comm(ii)

               mring_mysne(:,k)=mring_mysne(:,k)+SNmetarecv(ista:iend,irecv_sta)
               !if (verbose_ref) print *,'MPI>>', myid,ii,k,sngl(mring_mysne(:,k))

               ista=iend+1
               iend=ista+nshell_resolve
               mzring_mysne(:,k)=mzring_mysne(:,k)+SNmetarecv(ista:iend,irecv_sta)
            end if
         end do ! ii=1,ncomm_send
      end do ! icpu

   endif ! ncomm_send>0


   !===========================================================
   ! Step V : calculate the minimum radius we should refine 
   !                     and send it to other relevant cpus
   !===========================================================
   ! compute the minimum radius we should refine
   if (ncomm_send>0)then  
      do i=1,uidSN_comm
         j=0
         do while (j<=nshell_resolve)
            kend = nrefnei_ring(j)
            ! test if the local cell contains already more than enough gas mass
            nH = sum(mring_mysne(0:j,i))*msun2g*0.76/mH/(dble(kend)*(dx_loc*scale_l)**3)
            mgas = sum(mring_mysne(0:j,i))
            Zgas = sum(mzring_mysne(0:j,i))/mgas
            Zgas = max(Zgas/0.02,0.01) ! solar scaled
            if(.not.metal) Zgas=max(Zgas,z_ave) 
         
            ! transition mass from the adiabatic to snowplow phase (Kimm & Cen 2014, A5)
            chi_tr = (A_SN/1d4)**2/f_esn/M_SNII*nsn_resolve**((expE_SN-1d0)*2d0)*nH**(expN_SN*2d0)*Zgas**(expZ_SN*2d0)
            mass_tr = (chi_tr-1)*M_SNII/reduce_mass_tr ! Msun
            if(mgas>mass_tr)then
               lring_mysne(i)=j
               !print *, 'MPI>>>>', ilevel, j,myid
               j=nshell_resolve+1 ! escape from here
            endif
            j=j+1
         end do ! j=0,nshell_resolve
      end do ! i=1,uidSN_comm
   end if
        
   ! send the info
   if (ncomm_send>0)then
      reqsend=0; statsend=0
      allocate(lsend(1:ncomm_send))
      do icpu=1,ncpu_send
         lsend=-1
         cpu2send = list2send(icpu)
         ncc=0 ! No. of comm cells for cpu2send
         do i=1,ncomm_send
            j=i+isend_sta
            if(icpuSN_comm_mpi(j,2)==cpu2send)then
               ncc=ncc+1
               k = idSN_comm(i)
               lsend(ncc)=lring_mysne(k)
            end if
         end do ! i=1,ncomm_send

         tag = myid + cpu2send + ncc+30000
         call MPI_ISEND (lsend(1:ncc),ncc,MPI_INTEGER, &
                       & cpu2send-1,tag,MPI_COMM_WORLD,reqsend(icpu),info)
      end do ! icpu
   endif  ! ncomm_send>0

   ! receive the info
   if (ncomm_recv>0)then
      reqrecv=0; statrecv=0
      allocate(lrecv(1:ncomm_recv))
      irecv_sta=1
      lrecv=0d0
      do icpu=1,ncpu_recv
         cpu2recv = list2recv(icpu)
         ncc=0 ! No. of comm cells for cpu2recv
         do i=1,ncomm_tot
            if( (icpuSN_comm_mpi(i,1)==cpu2recv.and.&
              &  icpuSN_comm_mpi(i,2)==myid)   ) then
               ncc=ncc+1
            endif
         end do !i
         irecv_end=irecv_sta+ncc-1
         tag = myid + cpu2recv + ncc + 30000
         call MPI_IRECV (lrecv(irecv_sta:irecv_end),ncc,MPI_INTEGER, &
                      &  cpu2recv-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)
         irecv_sta=irecv_end+1
      end do !icpu
   endif

   ! wait
   if(ncpu_send>0) call MPI_WAITALL(ncpu_send,reqsend,statsend,info)
   if(ncpu_recv>0) call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)

   !if(ncpu_recv>0)print *, '@@@', myid,lrecv(1:irecv_end)
 

   !===========================================================
   ! Step VI: refine appropriately
   !===========================================================
   ! refine for the SN cells belonging to myid
   if (ncomm_send>0)then 

      do i=1,uidSN_comm
         ! this is to find the position of this SN 
         ii=1
         do while (ii<=ncomm_send)
            j=ii+isend_sta
            if(icpuSN_comm_mpi(j,2)==cpu2send)then
               xSN_i(1:3) = xSN_comm(1:3,ii)
               ii=ncomm_send+1 ! escape
            endif
            ii=ii+1
         end do

         ! Loop over a ring (start from the SNcell to rings with a larger radius)
         j=0
         do while (j<=lring_mysne(i))
         
            ksta = 1
            if(j>0) ksta = nrefnei_ring(j-1)+1
            kend = nrefnei_ring(j)
         
            ! Loop over cells within a ring
            do k=ksta,kend
         
               kk=irefnei(k)
               xc3(1:3) = xSN_i(1:3) + xrefnei(1:3,kk)*dx_loc
                  
               call get_icell_from_pos (xc3,ilevel,igrid,icell,ilevel2)

               !if this cell belongs to myid    
               if(cpu_map(father(igrid))==myid)then

                  ! TODO: do I need to worry that level-1 grids are updated here?
                  ! checked above that mgas>mass_tr, so just refine
                  ! 1) mass>mass_tr before j < nshell_resolve (refine it. simple.) 
                  ! 2) mass>mass_tr at j=nshell_resolve (refine only if lrefnei < ilevel)
                  ! 3) mass<mass_tr (refine only if lrefnei < ilevel)
                  if(j<nshell_resolve)then
                     if(flag1(icell)==0)then
                        nflag=nflag+1
                        flag1(icell)=1
                     end if
                  else
                     if(flag1(icell)==0.and.ilevel2<ilevel)then
                        nflag=nflag+1
                        flag1(icell)=1
                     end if
                  endif
               endif

            end do ! k: cells within a ring
         
            j=j+1 
         end do ! j: loop over a ring 

      end do ! i=1,uidSN_comm

   endif ! ncomm_send>0 


   ! refine for the SN cells belonging to other cpus
   if (ncomm_recv>0)then 
      do i=1,ncomm_recv
         xSN_i(1:3) = SNrecv(1:3,i)
   
         ! Loop over a ring (start from the SNcell to rings with a larger radius)
         j=0
         do while (j<=lrecv(i))
   
            ksta = 1
            if(j>0) ksta = nrefnei_ring(j-1)+1
            kend = nrefnei_ring(j)
    
            ! Loop over cells within a ring
            do k=ksta,kend
   
               kk=irefnei(k)
               xc3(1:3) = xSN_i(1:3) + xrefnei(1:3,kk)*dx_loc
            
               call get_icell_from_pos (xc3,ilevel,igrid,icell,ilevel2)
   
               if(cpu_map(father(igrid))==myid)then
                  ! TODO: do I need to worry that level-1 grids are updated here?
                  ! checked above that mgas>mass_tr, so just refine
                  ! 1) mass>mass_tr before j < nshell_resolve (refine it. simple.) 
                  ! 2) mass>mass_tr at j=nshell_resolve (refine only if lrefnei < ilevel)
                  ! 3) mass<mass_tr (refine only if lrefnei < ilevel)
                  if(j<nshell_resolve)then
                     if(flag1(icell)==0)then
                        nflag=nflag+1
                        flag1(icell)=1
                     end if
                  else
                     if(flag1(icell)==0.and.ilevel2<ilevel)then
                        nflag=nflag+1
                        flag1(icell)=1
                     end if
                  endif
               endif
    
            end do ! k: cells within a ring
   
            j=j+1 
         end do ! j: loop over a ring 
      end do ! i: loop over SNe

   endif ! ncomm_recv>0



   deallocate(icpuSN_comm_mpi,icpuSN_comm)
   if(ncomm_send>0) deallocate(list2send,SNsend,reqsend,statsend,SNmetarecv)
   if(ncomm_recv>0) deallocate(list2recv,SNrecv,reqrecv,statrecv,SNmetasend)
   if(ncomm_send>0) deallocate(mring_mysne,mzring_mysne,lring_mysne,lsend)
   if(ncomm_recv>0) deallocate(mring_other,mzring_other,lrecv)

end subroutine feedback_refine_fine_mpi
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################

