module Lya_commons
   use amr_commons
   use rt_parameters
   use hydro_parameters
   use rt_cooling_module, only:getMu
   implicit none

   public
   integer,parameter::nvarLya=8   ! number of variables to MPI communicate (x,y,z,nph,NHI,NHIZ)
   integer,parameter::n_nei=48
   real(dp),dimension(1:3,1:n_nei):: xrel_nei,v_nei
   integer, dimension(1:n_nei)    :: icpu_nei

   logical::Lya_initialised=.false.
   ! useful numbers
   real(dp),parameter::lam_Lya  = 1215.67   ! Angstrom
   real(dp),parameter::E_Lya_ev = 10.198815 ! eV
   real(dp),parameter::pc2cm    = 3.08d18
   real(dp)::Zsol_min = 1d-5       ! minimum metallicity even after dust destruction by thermal sputtering
   real(dp)::xHI_min  = 1d-5       ! below this, NHI will be set to 0
   real(dp)::v_limit_Lya = 100     ! kms maximum velocity to limit the acceleration due to Lya
   ! unit info
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale_msun

   ! grid info
   integer ::nx_loc
   real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc
   real(dp),dimension(1:twotondim,1:3)::xc
   integer ::nleaf_ilevel

   real(dp)::tyoung
   real(dp)::current_time

   integer ::ncomm_cell=0 ! number of cells to be communicated (specific to myid)

   ! 6 immediate neighbour info
   real(dp),dimension(1:nvector,0:twondim)::nH_nbor
   real(dp),dimension(1:nvector,0:twondim)::Tk_nbor
   real(dp),dimension(1:nvector,0:twondim)::xHI_nbor
   real(dp),dimension(1:nvector,0:twondim)::Zsol_nbor

#ifndef WITHOUTMPI
   integer,parameter::nvector_Lya=nvector*8*100 ! should be a multiple of nvector
   integer::ncomm_cell_tot
   real(kind=dp),dimension(1:3,1:nvector_Lya)::x_comm
   real(kind=dp),dimension(1:nvector_Lya):: nph_comm
   real(kind=dp),dimension(1:nvector_Lya):: nH_comm
   real(kind=dp),dimension(1:nvector_Lya):: Tk_comm
   real(kind=dp),dimension(1:nvector_Lya):: Zsol_comm
   real(kind=dp),dimension(1:nvector_Lya):: xHI_comm
   integer,dimension(1:nvector_Lya)::i_comm
   integer,dimension(:,:),allocatable::icpu_comm,icpu_comm_mpi
   integer,dimension(:  ),allocatable::ncomm_cpu,ncomm_mpi
#endif   

   contains
      !-----------------------------------------------------------
      subroutine initialise_Lya ()
      !-----------------------------------------------------------
      ! Arrays to define neighbors (center=[0,0,0])
      ! normalised to dx = 1 = size of the central leaf cell in which a particle sits
      ! from -0.75 to 0.75 
         implicit none
         integer::i,j,k,ind
         logical::ok
         real(kind=dp)::r,x,y,z
         ind=0
         do k=1,4
         do j=1,4
         do i=1,4
            ok=.true.
            if((i==1.or.i==4).and.&
               (j==1.or.j==4).and.&
               (k==1.or.k==4)) ok=.false. ! edge
            if((i==2.or.i==3).and.&
               (j==2.or.j==3).and.&
               (k==2.or.k==3)) ok=.false. ! centre
            if(ok)then
               ind=ind+1
               x = (i-1)+0.5d0 - 2   
               y = (j-1)+0.5d0 - 2   
               z = (k-1)+0.5d0 - 2   
               r = dsqrt(dble(x*x+y*y+z*z))
               xrel_nei(1,ind) = x/2d0
               xrel_nei(2,ind) = y/2d0
               xrel_nei(3,ind) = z/2d0
               v_nei(1,ind) = x/r  
               v_nei(2,ind) = y/r  
               v_nei(3,ind) = z/r  
               !indall(i+(j-1)*4+(k-1)*4*4) = ind      
            endif
         enddo
         enddo
         enddo         
 
#ifndef WITHOUTMPI
         allocate(ncomm_mpi(1:ncpu))
         allocate(ncomm_cpu(1:ncpu))
#endif
         Lya_initialised=.true.

      end subroutine initialise_Lya

      !-----------------------------------------------------------
      subroutine get_xHI (pvar,xHI)
      !-----------------------------------------------------------
         implicit none
         real(kind=dp),intent(in),dimension(1:nvar):: pvar
         real(kind=dp),intent(out)::xHI
         real(kind=dp)::d
         d=pvar(1)
         if (isH2) then
            xHI = pvar(iIons+ixHI-1)/d
            if(rt_Lya_pressure_count_H2) then
               xHI = 1d0 - pvar(iIons+ixHII-1)/d
            endif 
         else
            xHI = 1d0 - pvar(iIons+ixHII-1)/d
         endif
         if(xHI.lt.1d-5) xHI=0.0
      end subroutine get_xHI

      !-----------------------------------------------------------
      subroutine fetch_temperature(pvar,Tk,fix,message)
      !-----------------------------------------------------------
         implicit none
         real(kind=dp),intent(inout)::pvar(1:NVAR)
         real(kind=dp),intent(out)::Tk
         logical::fix
         real(kind=dp)::d,u,v,w,e,ekk,eth,T2,T_floor
         real(kind=dp)::mu,xHII,xHeII,xHeIII
         character(len=*),optional:: message
         integer::ierr,irad
       
         d = pvar(1)
         u = pvar(2)/d
         v = pvar(3)/d
         w = pvar(4)/d
         e = pvar(5)
         ekk = 0.5d0*d*(u**2+v**2+w**2)
         eth = e-ekk
#if NENER>0
         do irad=1,nener
            eth = eth - pvar(ndim+2+irad) 
         end do
#endif
         T2  = eth/d*scale_T2*(gamma-1.0)
         if(T2<0.and.(.not.fix))then
            if(present(message)) print *,'TKERR in '//TRIM(message)
            print *,'T2 [K   ]=', T2
            print *,'nH [H/cc]= ',d*scale_nH
            print *,'u  [km/s]= ',u*scale_v/1d5
            print *,'v  [km/s]= ',v*scale_v/1d5
            print *,'w  [km/s]= ',w*scale_v/1d5
         endif
       
         T_floor = T2_star*(d*scale_nH/n_star)**(g_star-1.0)
         if (T2<T_floor.and.fix)then
            eth = T_floor*d/scale_T2/(gamma-1.0)
#if NENER>0
            do irad=1,nener
               eth = eth + pvar(ndim+2+irad) 
            end do
#endif
            e = eth + ekk 
            pvar(5) = e 
         endif

         xHII   = pvar(iIons+ixHII  -1)/d
         xHeII  = pvar(iIons+ixHeII -1)/d
         xHeIII = pvar(iIons+ixHeIII-1)/d
         mu     = getMu(xHII,xHeII,xHeIII,T2)
         Tk     = T2*mu

      end subroutine fetch_temperature

      !-----------------------------------------------------------
      subroutine grid_info (ilevel)
      !-----------------------------------------------------------
         implicit none 
         integer::ilevel,ind,ix,iy,iz
       
         nx_loc=(icoarse_max-icoarse_min+1)
         skip_loc=(/0.0d0,0.0d0,0.0d0/)
         skip_loc(1)=dble(icoarse_min)
         skip_loc(2)=dble(jcoarse_min)
         skip_loc(3)=dble(kcoarse_min)
         scale=boxlen/dble(nx_loc)
         dx=0.5D0**ilevel
         dx_loc=dx*scale
         vol_loc=dx_loc**ndim
       
         ! Cells center position relative to grid center position
         do ind=1,twotondim
            iz=(ind-1)/4
            iy=(ind-1-4*iz)/2
            ix=(ind-1-2*iy-4*iz)
            xc(ind,1)=(dble(ix)-0.5D0)*dx
            xc(ind,2)=(dble(iy)-0.5D0)*dx
            xc(ind,3)=(dble(iz)-0.5D0)*dx
         end do
      end subroutine grid_info

      !-----------------------------------------------------------
      subroutine fetch_icell (pos_cell,ilevel_max,igrid,icell,ilevel_out)
      !-----------------------------------------------------------
         !use amr_commons,only:xg,son,ncoarse
         implicit none
         real(dp),intent(in) ::pos_cell(1:ndim) ! position between 0 and 1
         integer, intent(in) ::ilevel_max  ! maximum level you want to search
         integer, intent(out)::icell    ! index of the cell
         integer, intent(out)::igrid    ! index of the grid
         integer, intent(out)::ilevel_out  ! level of this cell
         logical::not_found
         real(dp)::dx_here,x0(1:ndim),fpos(1:ndim)
         integer::idim,ind,i

         fpos = pos_cell   
         do idim=1,ndim
            if( fpos(idim).gt.1.0 ) fpos(idim)=fpos(idim)-1.0
            if( fpos(idim).lt.0.0 ) fpos(idim)=fpos(idim)+1.0
         enddo 
         not_found=.true.
         igrid = 1 ! this is level=1 grid
         ilevel_out = 1
         do while (not_found)
            dx_here = 0.5D0**ilevel_out
            x0(1:ndim) = xg(igrid,1:ndim)-dx_here-skip_loc(1:3) !!left corner of *this* grid [0-1], not cell
            ind  = 1
            do idim = 1,ndim
               i = int( (fpos(idim) - x0(idim))/dx_here)
               ind = ind + i*2**(idim-1)
            end do
            icell = igrid+ncoarse+(ind-1)*ngridmax
            if(son(icell)==0.or.ilevel_out==ilevel_max) return
            igrid=son(icell)
            ilevel_out=ilevel_out+1
         end do

      end subroutine fetch_icell

      !-----------------------------------------------------------
      subroutine identify_young_cells (ind_grid,ngrid,ilevel,with_young_stars)
      !-----------------------------------------------------------
         use pm_commons
         implicit none
         integer,dimension(1:nvector),intent(in)::ind_grid
         integer,intent(in)::ngrid,ilevel
         logical,dimension(1:nvector,1:twotondim),intent(out)::with_young_stars
         !----------------------------
         real(dp)::tyoung,current_time
         integer,save::ig,ind,iskip,nleaf,npart1,ipart,jpart,next_part
         integer,save::ind_cell,ind_leaf,idim,ind_son
         real(dp),dimension(1:ndim),save::x0
         logical::ok_star

         if(use_proper_time)then
           current_time=texp
           tyoung = current_time - t_sne*Myr2sec/(scale_t/aexp**2)
         else
           current_time=t
           tyoung = current_time - t_sne*Myr2sec/(scale_t)
         endif

         ! Initialise
         with_young_stars(:,:)=.false.

         ! Mesh info
         call grid_info (ilevel)

         ! Loop over grids
         do ig=1,ngrid

            ! Loop over cells
            nleaf=0
            do ind=1,twotondim
               iskip=ncoarse+(ind-1)*ngridmax
               ind_cell=iskip+ind_grid(ig)
               if(son(ind_cell)==0)then
                  nleaf=nleaf+1
               endif
            end do
 
            if(nleaf.eq.0)cycle
  
            npart1 = numbp(ind_grid(ig))

            if(npart1>0)then

               do idim=1,ndim
                  x0(idim)=xg(ind_grid(ig),idim) -dx -skip_loc(idim)
               end do

               ipart = headp(ind_grid(ig))

               do jpart=1,npart1
                  ! Save next particle   <--- Very important !!!
                  next_part=nextp(ipart)
                  ! select star particles
                  ok_star = tp(ipart).ne.0
             
                  if(ok_star)then ! check if young
                     if(tp(ipart).lt.tyoung) ok_star=.false.
                  endif
 
                  ! check where this belongs
                  if(ok_star)then
                      ind_son=1
                      do idim=1,ndim
                         ind = int((xp(ipart,idim)/scale - x0(idim))/dx)
                         ind_son = ind_son+ind*2**(idim-1)
                      end do
                      if(ind_son<1.or.ind_son>8)then
                         ! will be done using the streaming method
                      else
                         iskip=ncoarse+(ind_son-1)*ngridmax
                         ind_cell=iskip+ind_grid(ig)
                         if(son(ind_cell)==0)then  ! leaf cell
                             with_young_stars(ig,ind_son)=.true.
                         endif
                      endif
                  endif

                  ipart=next_part  ! Go to next particle

               end do ! do jpart=1,npart1
            end if  ! npart>1
            
         end do ! Loop over grids
      end subroutine identify_young_cells

      !-----------------------------------------------------------
      subroutine Lya_pressure_nbor_prop (ind_leaf,nleaf,ilevel)
      !-----------------------------------------------------------
         use hydro_commons, only:uold   
         implicit none
         integer,intent(in),dimension(1:nvector)::ind_leaf
         integer,intent(in)::nleaf,ilevel
         integer::ileaf,ibor
         real(dp)::Tk,nH,xHI
         integer ,dimension(1:nvector,0:twondim),save:: ind_nbor
         real(dp),dimension(1:NVAR),save::pvar
       
         call getnbor(ind_leaf,ind_nbor,nleaf,ilevel)

         nH_nbor=0.0;Tk_nbor=0.0;xHI_nbor=0.0;Zsol_nbor=0.0

         do ileaf=1,nleaf
            do ibor=1,6
               pvar(:) = uold(ind_nbor(ileaf,ibor),:) 
               nH = pvar(1)*scale_nH ! ul, ur, vl, vr, wl, wr
       
               call fetch_temperature(pvar,Tk,.false.)

               call get_xHI (pvar, xHI)
               
               nH_nbor(ileaf,ibor) = nH
               Tk_nbor(ileaf,ibor) = Tk
               xHI_nbor(ileaf,ibor) = xHI 
               Zsol_nbor(ileaf,ibor) = pvar(imetal)/pvar(1)/0.02
               
            end do
         end do

      end subroutine Lya_pressure_nbor_prop

      !-----------------------------------------------------------
      subroutine Lya_boost_factor(NHI_cm2,Zsolar,Tk,M_factor)
      !-----------------------------------------------------------
         ! compute the multiplication factor of the Lya pressure
         ! F = M_factor * ( L_Lya / c)
         implicit none
         real(kind=dp),intent(in)::NHI_cm2,Zsolar,Tk
         real(kind=dp),intent(out)::m_factor
         real(kind=dp)::Tk4,P_B_Lya,Zsolar1
         real(kind=dp)::tau0_Lya,sig0_Lya,tau0_peak,tau_Lya_da
         real(kind=dp)::xi,xi_new,av_Lya,fescLya,M_F_nodust,ltau0
         real(kind=dp)::dust2metal_solar
 
         Tk4 = min(max(Tk,10.),1e6) / 1d4 ! necessary because of P_B_Lya

         ! Optical depth at the line centre
         sig0_Lya = 5.88d-14/sqrt(Tk4)  ! cm^2
         tau0_Lya = NHI_cm2 * sig0_Lya 
      
         ! Be careful with Z=0 limit, bc timescale involved can be longer than dteff
         Zsolar1 = max(Zsolar,Zsol_min) 

         if(dust2metal_RR14)then
            call dust_to_metal_RR14 (Zsolar,dust2metal_solar)
         else
            dust2metal_solar = 1.0
         endif
        
         ! maximum optical depth in the presence of dust
         tau0_peak  = 4.06d6 * Tk4**(-0.25) * Zsolar1**(-0.75) * (dust2metal_solar*sigdust21/3.0)**(-0.75)
      
         ! Dust optical depth
         tau_Lya_da = NHI_cm2 * (dust2metal_solar*sigdust21*1d-21) * Zsolar1 * (1d0 - DustAlbedo)
         tau_Lya_da = tau_Lya_da * min(tau0_Lya,tau0_peak)/tau0_Lya ! dust optical depth at NHI=NHI_peak
         
         ! Limit the optical depth in the presence of dust
         if(tau0_Lya>tau0_peak) tau0_Lya = tau0_peak
         
         xi = 0.525                     ! Verhamme+(06): Eq.24
         xi_new = 1.78                  ! Adjusted to fit the RASCAS results
         av_Lya = 4.7d-4/sqrt(TK4)      ! Voigt parameter
         fescLya = 1d0/(cosh(sqrt(3.)/3.141592**(5./12.)/xi_new*sqrt( (av_Lya*tau0_Lya)**(1./3.)*tau_Lya_da))) 
      
         ! Results from RASCAS
         if(tau0_Lya.ge.1d6)then
            M_F_nodust = 29 * (tau0_Lya/1d6)**(0.29)
         else
            ltau0 = log10(tau0_Lya)
            M_F_nodust = 10d0**(-0.433+0.874*ltau0-0.1732*(ltau0)**2.0+0.0133*(ltau0)**3.0)
         endif
         !if(fescLya<0)then
         !   write(*,*) '>> Error: fescLya <0'
         !   call stop_mpi
         !endif
      
         M_factor = max(M_F_nodust * fescLya, 0.d0)

      end subroutine Lya_boost_factor

      !-----------------------------------------------------------
      subroutine cmp_P_B_Lya (Tk, P_B_Lya)
      !-----------------------------------------------------------
         implicit none
         real(dp),intent(in)::Tk
         real(dp)::Tk4
         real(dp),intent(out)::P_B_Lya
         ! compute the transition probability from LyC to LyA
         ! Cantalupo+(2008) (100<T<10^5)
         Tk4 = min(max(Tk,10.),1e6)/1e4
         P_B_Lya = (0.686 - 0.106*log10(TK4) - 0.009 *(TK4)**(-0.44)) 
      end subroutine cmp_P_B_Lya

#ifndef WITHOUTMPI
      !-----------------------------------------------------------
      subroutine Lya_comm_init
      !-----------------------------------------------------------
      implicit none
      ncomm_cell=0
      ncomm_cell_tot=0
      x_comm (:,:)=0
      nph_comm (:)=0
      nH_comm  (:)=0
      Tk_comm  (:)=0
      Zsol_comm(:)=0
      xHI_comm (:)=0
      i_comm   (:)=0
      end subroutine Lya_comm_init
#endif

end module Lya_commons
!####################################################################
!####################################################################
!####################################################################
subroutine rt_Lya_uold_init (ilevel)
  use amr_commons
  use hydro_commons,only:uold,nvar
  use rt_parameters
  use Lya_commons,only:nleaf_ilevel
  implicit none
  integer,intent(in)::ilevel
  integer::ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell
  
  ! set uold(:,iLyaVar)=0.0
  nleaf_ilevel = 0
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid = active(ilevel)%igrid(igrid+i-1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ind_cell = iskip + ind_grid
           if(son(ind_cell)==0)then
              uold(ind_cell,iLyaVar)=0.0
              nleaf_ilevel=nleaf_ilevel+1
           endif
        end do ! ind=1,8
     end do ! i=1,ngrid
  end do

end subroutine rt_Lya_uold_init
!####################################################################
!####################################################################
!####################################################################
subroutine rt_Lya_feedback (ilevel)
  use amr_commons
  use Lya_commons
  use rt_parameters
  use pm_parameters,only:nstar_tot
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine computes the pressure from Lya scattering
  ! and inject momentum and energy to the surroundings.
  ! This routine is called every fine time step.
  ! ind_pos_cell: position of the cell within an oct [1-8]
  ! nobnd: number of cells around the domain boundary, hence need MPI
  !------------------------------------------------------------------------
  integer,intent(in)::ilevel
  !------------------------------------------------------------------------
  integer,dimension(1:nvector)::ind_grid
  integer::i,j,idim,icell,igrid,ncache,ngrid
  integer::ncell,ngrid_done
  real(dp)::NLyC_min
#ifndef WITHOUTMPI
  integer::info,ncache_max,n_MPI,i_MPI,ii
#endif 

  if(nstar_tot==0)return
  if(ndim.ne.3)  return
  if(numbtot(1,ilevel)==0)return

  if(.not.Lya_initialised) call initialise_Lya()

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/m_sun

  ! Mesh spacing in this level
  call grid_info (ilevel)

#ifndef WITHOUTMPI
  ! initialise the communication variables
  call Lya_comm_init
#endif

  ngrid_done=0

  ncache=active(ilevel)%ngrid
#ifndef WITHOUTMPI
  ncache_max=0
  call MPI_ALLREDUCE(ncache,ncache_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  n_MPI=ceiling(ncache_max/real(nvector_lya)) ! how many times do we need to call lya_fine_mpi?
  i_MPI=0
#endif

  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call Lya_fine1 (ind_grid,ngrid,ilevel)

#ifndef WITHOUTMPI
     ngrid_done=ngrid_done+ngrid
     if(ngrid_done==nvector_Lya)then
        ncomm_cell_tot=0
        i_MPI=i_MPI+1
        ! Deal with the stars around the bounary of each cpu (need MPI)
        call MPI_ALLREDUCE(ncomm_cell,ncomm_cell_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        if(ncomm_cell_tot>0)then
           !if(mod(nstep_coarse,ncontrol)==0.and.myid.eq.1)print *, '====> Ncomm_cell_tot =',ncomm_cell_tot, i_MPI
           call rt_Lya_fine1_mpi(ilevel)
        endif
        call Lya_comm_init
        ngrid_done=0
     endif
#endif

  end do

#ifndef WITHOUTMPI
  do ii=i_MPI+1,n_MPI
     ncomm_cell_tot=0
     ! Deal with the stars around the bounary of each cpu (need MPI)
     call MPI_ALLREDUCE(ncomm_cell,ncomm_cell_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     if(ncomm_cell_tot>0)then
        !if(mod(nstep_coarse,ncontrol)==0.and.myid.eq.1)print *, '====> Ncomm_cell_tot =',ncomm_cell_tot, i_MPI
        call rt_Lya_fine1_mpi(ilevel)
     endif
     call Lya_comm_init
  end do
#endif



end subroutine rt_Lya_feedback
!####################################################################
!####################################################################
!####################################################################
subroutine Lya_fine1 (ind_grid,ngrid,ilevel)
  use amr_commons
  use Lya_commons
  use rt_parameters
  use rt_hydro_commons
  use hydro_commons,only:uold,unew,nvar
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in),dimension(1:nvector)::ind_grid
  integer,intent(in)::ilevel,ngrid
  !--------------------------------------------------------
  logical,dimension(1:nvector,1:twotondim),save::with_young_stars
  logical,dimension(1:nvector),save::leaf_with_young_stars
  integer,dimension(1:nvector),save::ind_leaf,ind_grid_leaf
  integer,save::ind,ileaf,nleaf,iskip,iNp
  real(dp),dimension(1:ndim),save::xcen,xnei
  integer,save::ig,i,j,nobnd,nobnd_dum,idim,igrid,ista,iend,ilevel2
  integer,save::icell,ncell
  real(dp),dimension(1:nvar),save::pvar,pvar_nei
  ! For stars affecting across the boundary of a cpu
  integer,save::icell_nei,igrid_nei,ilevel_nei,level_diff
  real(dp),save::nH_cen,Zsol_cen,xHI_cen,Tk_cen,dx_loc_cm
  real(dp),save::NLyC_min,NLya50_cell,Tk4,P_B_Lya
  real(dp),dimension(1:ndim),save::Fp,Fpnew,Fhat,Fhat_Fp
  logical::ok,ok_Lya
#ifndef WITHOUTMPI
  integer::info
  real(dp)::ttend,ttsta
#endif 

  dx_loc_cm = dx_loc*scale_l

  ! find out which cells host young star particles
  with_young_stars(:,:)=.false.
  leaf_with_young_stars(:)=.false.
  call identify_young_cells (ind_grid,ngrid,ilevel,with_young_stars)

  ! Loop over cells
  do ind=1,twotondim

     iskip=ncoarse+(ind-1)*ngridmax

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        icell = iskip+ind_grid(i)
        ok = (son(icell)==0)
        if(ok)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=icell
           ind_grid_leaf(nleaf)=ind_grid(i)
#ifdef RT
           if(rt_Lya_pressure) then
              leaf_with_young_stars(nleaf)=with_young_stars(i,ind)
           endif
#endif
        endif
     end do
  
     if(nleaf.eq.0)cycle

     ! Neighbouring cell for fluid method
     call Lya_pressure_nbor_prop(ind_leaf,nleaf,ilevel) 


     ! Loop over leaf cells containing enough NLya 
     do ileaf=1,nleaf 

        icell = ind_leaf(ileaf) 

        pvar(1:nvar) = uold(icell,1:nvar)

        call fetch_temperature (pvar,Tk_cen,.false.,'Lya_fine1')

        ! central NHI and metallicity
        call get_xHI (pvar, xHI_cen)
        nH_cen = pvar(1)*scale_nH
        Zsol_cen = (pvar(imetal)/pvar(1)/0.02)

        NLya50_cell = uold(icell,iLyaVar) ! in units of 10^50 #
   
        nobnd=0; icpu_nei=0
  
        if(leaf_with_young_stars(ileaf)) then

           ! Record position of each cell [0.0-1.0] regardless of boxlen
           xcen=0d0
           do idim=1,ndim
              xcen(idim)=xg(ind_grid_leaf(ileaf),idim)-skip_loc(idim)+xc(ind,idim)
           end do 
      
           ! figure out the index of the host cell. Should be the same as ind_cell
           call fetch_icell (xcen, ilevel+1, igrid, icell,ilevel2)
     
           ! check the grid is organised properly (sanity check)
           if((cpu_map(father(igrid)).ne.myid).or.&
              (ilevel.ne.ilevel2).or.&
              (ind_leaf(ileaf).ne.icell))     then 
              print *,'>>> fatal error in Lya_fine1'
              print *, cpu_map(father(igrid)),myid
              print *, ilevel, ilevel2
              print *, ind_leaf(ileaf), icell
              call stop_mpi
           endif
 
           ! loop over 48 neighbouring cells
           do j=1,n_nei
      
              ! fetch the index of the neighbour cells
              do idim=1,ndim
                 xnei(idim) = xcen(idim)+xrel_nei(idim,j)*dx
              end do
              call fetch_icell (xnei, ilevel+1, igrid_nei, icell_nei, ilevel_nei)
   
              if(cpu_map(father(igrid_nei)).ne.myid) then ! need mpi
                 nobnd=nobnd+1
                 icpu_nei(nobnd)=cpu_map(father(igrid_nei))
              else  ! can be handled locally
      
                 pvar_nei(:) = 0d0 ! temporary primitive variable for neighbour cell
                 if(ilevel.gt.ilevel2)then ! touching level or level-1 cells
                    pvar_nei(1:nvar) = unew(icell_nei,1:nvar)
                 else
                    pvar_nei(1:nvar) = uold(icell_nei,1:nvar)
                 endif 
      
                 level_diff = ilevel - ilevel_nei
                 !----------------------------------------------------------------------
                 call inject_momentum_Lya_neigh (nH_cen,Tk_cen,Zsol_cen,xHI_cen,pvar_nei,level_diff,j,NLya50_cell)
                 !----------------------------------------------------------------------
      
                 ! update the hydro variable
                 if(ilevel.gt.ilevel2)then ! touching level-1 cells
                    unew(icell_nei,1:nvar) = pvar_nei(1:nvar)
                 else
                    uold(icell_nei,1:nvar) = pvar_nei(1:nvar)
                 endif 
      
              end if
           end do ! loop over 48 neighbors


#ifndef WITHOUTMPI 
           if(nobnd>0)then  ! for particles affecting different cpu
              if(nobnd>1)then
                 nobnd_dum=nobnd
                 ! remove redundant cpu list for this cell
                 call redundant_non_1d(icpu_nei(1:nobnd), nobnd_dum, nobnd)
              endif
              ista=ncomm_cell+1
              iend=ista+nobnd-1
      
              if(iend>nvector_Lya*8)then
                 write(*,*) 'Error: increase nvector_Lya in rt_Lya_pressure.f90', nvector_Lya, iend
                 call MPI_ABORT(MPI_COMM_WORLD,1,info)
              endif
      
              i_comm (ista:iend)=icpu_nei(1:nobnd)
              x_comm (1,ista:iend)=xcen(1)
              x_comm (2,ista:iend)=xcen(2)
              x_comm (3,ista:iend)=xcen(3)
              nph_comm (ista:iend)=NLya50_cell
              nH_comm  (ista:iend)=nH_cen 
              Tk_comm  (ista:iend)=Tk_cen
              Zsol_comm(ista:iend)=Zsol_cen
              xHI_comm (ista:iend)=xHI_cen 
              ncomm_cell=ncomm_cell+nobnd
           endif
#endif
        else ! cells without young stars
           ! Fluid method

           ! Get the direction of photon groups
           Fhat=0.0; Fhat_Fp=0.0 
           do ig=1,nGroups
              ! units are omitted since we are interested in direction vector
              if(group_egy(ig).ge.13.6)then
                 iNp=iGroups(ig)
                 Fp =rtuold(icell,iNp+1:iNp+ndim)
                 Fhat_Fp = Fhat_Fp + Fp
                 if(rt_smooth) then
                    Fpnew =rtunew(icell,iNp+1:iNp+ndim)
                    Fhat = Fhat + Fpnew(:) - Fp(:)
                 endif
              end if
           end do 

           ! check if we want to inject momentum or normalise Fhat
           if(sum(Fhat*Fhat).eq.0)then
              Fhat = Fhat_Fp
           endif
           if(sum(Fhat*Fhat).eq.0)then
              ok_Lya=.false.
           else
              ok_Lya=.true.
              Fhat = Fhat / sqrt(sum(Fhat*Fhat))
              if (NLya50_cell<1d-20) ok_Lya=.false.  ! little absorption
           endif


           if(ok_Lya)then
              ! inject momentum based on the central+immedimate neighbour (ileaf is to link with nbor)
              call inject_momentum_Lya_cen ( ileaf, pvar, nH_cen, Tk_cen, xHI_cen, Zsol_cen, Fhat, NLya50_cell) 

              ! update the hydro variable
              uold(icell,1:nvar) = pvar(1:nvar) 
           endif
 
        end  if ! if leaf cell with young stars
     end do ! loop over leaf
  end do ! loop over ind


end subroutine Lya_fine1
!################################################################
!################################################################
!################################################################
#ifndef WITHOUTMPI
subroutine rt_Lya_fine1_mpi(ilevel)
  use amr_commons
  use Lya_commons
  use hydro_commons,only:nvar,uold,unew
  implicit none
  include 'mpif.h'
  integer::i,j,info,ncomm_tot,icpu,ncpu_send,ncpu_recv,ncc
  integer::ncell_recv,ncell_send,cpu2send,cpu2recv,tag,np
  integer::isend_sta,irecv_sta,irecv_end
  real(dp),dimension(:,:),allocatable::Asend,Arecv
  integer ,dimension(:),allocatable::list2recv,list2send
  integer, dimension(:),allocatable::reqrecv,reqsend
  integer, dimension(:,:),allocatable::statrecv,statsend
  real(dp),save::pvar_nei(1:nvar),NLya50_cell
  real(dp),dimension(1:ndim),save::xcen,xnei
  integer::igrid_nei,icell_nei,ilevel,ilevel_nei,idim,level_diff
  real(dp)::nH_cen,Tk_cen,Zsol_cen,xHI_cen
  call grid_info(ilevel)

  !============================================================
  ! For MPI communication
  !============================================================
  ncpu_send=0;ncpu_recv=0

  ncomm_cpu=0 
  ncomm_mpi=0
  ncomm_mpi(myid)=ncomm_cell
  ! compute the total number of communications needed
  call MPI_ALLREDUCE(ncomm_mpi,ncomm_cpu,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ncomm_tot = sum(ncomm_cpu)
  if(ncomm_tot==0) return

  allocate(icpu_comm    (1:ncomm_tot,1:2))
  allocate(icpu_comm_mpi(1:ncomm_tot,1:2))

  ! index for mpi variable
  if(myid==1)then
     isend_sta = 0 
  else
     isend_sta = sum(ncomm_cpu(1:myid-1)) 
  endif

  icpu_comm=0
  do i=1,ncomm_cpu(myid)
     icpu_comm(isend_sta+i,1)=myid
     icpu_comm(isend_sta+i,2)=i_comm(i)
     ! i_comm:   local variable
     ! icpu_comm:  local (but extended) variable to be passed to a mpi variable
     ! icpu_comm_mpi: mpi variable
  end do

  ! share the list of communications
  icpu_comm_mpi=0
  call MPI_ALLREDUCE(icpu_comm,icpu_comm_mpi,ncomm_tot*2,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  ncell_send = ncomm_cpu(myid)
  ncell_recv = count(icpu_comm_mpi(:,2).eq.myid, 1)

  ! check if myid needs to send anything
  if(ncell_send>0)then
     allocate( Asend (1:nvarLya,1:ncell_send)) ! x(3),nph
     allocate( list2send(1:ncell_send)    )
     list2send=0;  Asend=0d0
     list2send=icpu_comm_mpi(isend_sta+1:isend_sta+ncell_send,2)
     ncpu_send=1
     if(ncell_send>1) call redundant_non_1d (list2send,ncell_send,ncpu_send)
     ! ncpu_send = No. of cpus to which myid should send info 
     allocate( reqsend  (1:ncpu_send) )
     allocate( statsend (1:MPI_STATUS_SIZE,1:ncpu_send) )
     reqsend=0; statsend=0
  endif

  ! check if myid needs to receive anything
  if(ncell_recv>0)then
     allocate( Arecv (1:nvarLya,1:ncell_recv)) ! x(3),nph
     allocate( list2recv(1:ncell_recv)    )
     list2recv=0;  Arecv=0d0
     j=0
     do i=1,ncomm_tot
        if(icpu_comm_mpi(i,2).eq.myid)then
           j=j+1 
           list2recv(j) = icpu_comm_mpi(i,1)
        endif
     end do
 
     ncc = j
     if(ncc.ne.ncell_recv)then ! sanity check
        write(*,*) 'Error in mech_fine_mpi: ncc != ncell_recv',ncc,ncell_recv,myid
        call stop_mpi
     endif

     ncpu_recv=1 ! No. of cpus from which myid should receive info 
     if(j>1) call redundant_non_1d(list2recv,ncc,ncpu_recv)

     allocate( reqrecv  (1:ncpu_recv) )
     allocate( statrecv (1:MPI_STATUS_SIZE,1:ncpu_recv) )
     reqrecv=0; statrecv=0
  endif

  ! prepare one variable and send
  if(ncell_send>0)then
     do icpu=1,ncpu_send
        cpu2send = list2send(icpu)
        ncc=0 ! number of star-hosting cells that need communications with myid=cpu2send
        do i=1,ncell_send
           j=i+isend_sta
           if(icpu_comm_mpi(j,2).eq.cpu2send)then
              ncc=ncc+1
              Asend(1:3,ncc)=x_comm (1:3,i)
              Asend(4  ,ncc)=nph_comm  (i)
              Asend(5  ,ncc)=nH_comm   (i)
              Asend(6  ,ncc)=Tk_comm   (i)
              Asend(7  ,ncc)=Zsol_comm (i)
              Asend(8  ,ncc)=xHI_comm  (i)
           endif
        end do ! i

        tag = myid + cpu2send + ncc
        call MPI_ISEND (Asend(1:nvarLya,1:ncc),ncc*nvarLya,MPI_DOUBLE_PRECISION, &
                      & cpu2send-1,tag,MPI_COMM_WORLD,reqsend(icpu),info) 
     end do ! icpu

  endif ! ncell_send>0


  ! receive one large variable
  if(ncell_recv>0)then
     irecv_sta=1
     do icpu=1,ncpu_recv
        cpu2recv = list2recv(icpu)
        ncc=0 ! number of star-hosting cells that need communications with cpu2recv
        do i=1,ncomm_tot
           if((icpu_comm_mpi(i,1)==cpu2recv).and.&
             &(icpu_comm_mpi(i,2)==myid)      )then
              ncc=ncc+1
           endif
        end do
        irecv_end=irecv_sta+ncc-1
        tag = myid + cpu2recv + ncc
        
        call MPI_IRECV (Arecv(1:nvarLya,irecv_sta:irecv_end),ncc*nvarLya,MPI_DOUBLE_PRECISION,&
                     & cpu2recv-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)

        irecv_sta=irecv_end+1
     end do ! icpu 

  endif ! ncell_recv >0

  if(ncpu_send>0)call MPI_WAITALL(ncpu_send,reqsend,statsend,info)
  if(ncpu_recv>0)call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)


  !------------------------------------------------------------
  ! Apply the Lya pressure to this cpu domain
  !------------------------------------------------------------
  np = ncell_recv
  do i=1,np

     xcen(1:3)   = Arecv(1:3,i)
     NLya50_cell = Arecv(4,i)
     nH_cen      = Arecv(5,i)
     Tk_cen      = Arecv(6,i)
     Zsol_cen    = Arecv(7,i)
     xHI_cen     = Arecv(8,i)
     ! Loop over 48 neighbours
     do j=1,n_nei
 
        do idim=1,ndim
           xnei(idim) = xcen(idim)+xrel_nei(idim,j)*dx
        end do

        call fetch_icell (xnei, ilevel+1, igrid_nei, icell_nei, ilevel_nei)

        if(cpu_map(father(igrid_nei))==myid)then

           if(ilevel.gt.ilevel_nei)then ! touching level-1 cells
              pvar_nei(1:nvar) = unew(icell_nei,1:nvar)
           else
              pvar_nei(1:nvar) = uold(icell_nei,1:nvar)
           endif
 
           level_diff = ilevel - ilevel_nei
           !----------------------------------------------------------------------
           call inject_momentum_Lya_neigh (nH_cen,Tk_cen,Zsol_cen,xHI_cen,pvar_nei,level_diff,j,NLya50_cell)
           !----------------------------------------------------------------------
 
           ! update the hydro variable
           if(ilevel.gt.ilevel_nei)then ! touching level-1 cells
              unew(icell_nei,1:nvar) = pvar_nei(1:nvar)
           else
              uold(icell_nei,1:nvar) = pvar_nei(1:nvar)
           endif 

        endif ! if this belongs to me


     end do ! loop over neighbors
  end do ! loop over  cells


  deallocate(icpu_comm_mpi, icpu_comm)
  if(ncell_send>0) deallocate(list2send,Asend,reqsend,statsend)
  if(ncell_recv>0) deallocate(list2recv,Arecv,reqrecv,statrecv)

  ncomm_cell=ncomm_tot

end subroutine rt_Lya_fine1_mpi
#endif
!################################################################
!################################################################
!################################################################
subroutine inject_momentum_Lya_neigh (nH_cen,Tk_cen,Zsol_cen,xHI_cen,pvar_nei,level_diff,i_nei,Nph50)
   use Lya_commons
   use hydro_commons,only:gamma,nvar
   use rt_cooling_module,only:getMu
   implicit none
   real(kind=dp),dimension(1:nvar),intent(inout)::pvar_nei
   real(kind=dp),intent(in)::Nph50,nH_cen,Tk_cen,xHI_cen,Zsol_cen
   integer,intent(in)::i_nei,level_diff
   real(kind=dp)::d,u,v,w,ekk,eth,ekk0,eth0,ekk_add
   real(kind=dp)::vol_nei,dx_loc_nei
   real(kind=dp)::NH_Zsol,M_factor,Zsol_nei
   real(kind=8) ::xHI_nei,xHII,xHeII,xHeIII,T2,Tk_nei,mu
   real(kind=dp)::p_cgs,p_code,scale_kms,v_Lya
   real(kind=dp)::du,dv,dw,dNph,Lsob,nH_nei
#if NENER>0
   integer::irad
#endif

   d=pvar_nei(1)

   ! omit 
   if(d*scale_nH<rt_Lya_pressure_nH) return

   vol_nei = vol_loc*(2d0**ndim)**(level_diff)
   dx_loc_nei = dx_loc * (2d0)**(level_diff)

   u=pvar_nei(2)/d
   v=pvar_nei(3)/d
   w=pvar_nei(4)/d
   ekk0=0.5d0*d*(u*u + v*v + w*w)
   eth0=pvar_nei(5)-ekk0
#if NENER>0
   do irad=1,nener
      eth0=eth0-pvar_nei(ndim+2+irad)
   end do
#endif
   T2  = eth0/d*scale_T2*(gamma-1)

   xHII   = pvar_nei(iIons+ixHII  -1)/d
   xHeII  = pvar_nei(iIons+ixHeII -1)/d
   xHeIII = pvar_nei(iIons+ixHeIII-1)/d
   mu     = getMu(xHII,xHeII,xHeIII,T2)
   Tk_nei = T2*mu

   ! Estimate column density of neighbour
   call get_xHI (pvar_nei, xHI_nei)
  
   ! Estimate multiplification factor 
   nH_nei = d*scale_nH
   Zsol_nei = pvar_nei(imetal)/d/0.02

   call estimate_boost_factor_Lya ( nH_cen, Tk_cen, Zsol_cen, xHI_cen, &
                                  & nH_nei, Tk_nei, Zsol_nei, xHI_nei, M_factor )

   ! momentum from Lya pressure
   dNph  = Nph50 * 1d50  / n_nei
   p_cgs = M_factor*dNph*E_Lya_ev*ev_to_erg/c_cgs  ! erg/(cm/s) = g*cm/s
   p_code = p_cgs / (scale_l**3*scale_d) / scale_v
 
   v_Lya = p_code / (vol_loc*d)
   scale_kms = scale_v/1d5
   if( (v_Lya*scale_kms).gt.v_limit_Lya) then
3998  format('TKWARN v= ',f9.2,1x,' M_F= ',f9.1,' lnH= ',f6.2,1x,f6.2,' lTk= ',f6.1,1x,f6.1,' lZsol=',f6.1,1x,f6.1,' lxHI=',f6.1,1x,f6.1)
!      write(*,3998) v_Lya*scale_kms, M_factor,log10(nH_cen),log10(nH_nei), &
!                 &  log10(Tk_cen),log10(Tk_nei),log10(Zsol_cen),log10(Zsol_nei), &
!                 &  log10(xHI_cen),log10(xHI_nei)
      p_code = (v_limit_Lya/scale_kms) * (vol_nei*d)
   endif
   
   du = p_code * v_nei(1,i_nei) / (vol_nei*d)
   dv = p_code * v_nei(2,i_nei) / (vol_nei*d)
   dw = p_code * v_nei(3,i_nei) / (vol_nei*d)

   u  = u + du
   v  = v + dv
   w  = w + dw

   ekk_add = 0.0
   ekk_add = ekk_add + 0.5*d*du*du
   ekk_add = ekk_add + 0.5*d*dv*dv
   ekk_add = ekk_add + 0.5*d*dw*dw

   ! new momentum and total energy
   pvar_nei(1) = d
   pvar_nei(2) = d*u
   pvar_nei(3) = d*v
   pvar_nei(4) = d*w
   ekk = 0.5*d*(u*u + v*v + w*w)
#if NENER>0
   do irad=1,nener
      eth0=eth0+pvar_nei(ndim+2+irad)
   end do
#endif
   pvar_nei(5) = max(eth0+ekk,eth0+ekk0+ekk_add)

end subroutine inject_momentum_Lya_neigh
!################################################################
!################################################################
!################################################################
subroutine inject_momentum_Lya_cen (ileaf, pvar, nH_cen, Tk_cen, xHI_cen, Zsol_cen, Fhat, Nph50)
   use amr_parameters
   use Lya_commons
   implicit none
   integer, intent(in)::ileaf
   real(dp),intent(in)::nH_cen,Tk_cen,xHI_cen,Zsol_cen, Nph50
   real(dp),dimension(1:ndim),intent(in)::Fhat
   real(dp),dimension(1:NVAR),intent(inout)::pvar
   integer,save::ind,ind_nbor,idim,irad
   real(dp),save::dx_loc_cm,nH_tot,nH_nei,Tk_nei,Zsol_nei,xHI_nei
   real(dp),save::Zsol_avg,xHI_avg,Tk_avg,NHI_cm2_sob,Lsob,M_factor
   real(dp),save::Nph,d,u,v,w,ekk0,eth0,ekk,eth,T2,p_cgs,p_code,v_Lya
   real(dp),save::scale_kms,uadd,ekk_add


   d=pvar(1)

   if(d*scale_nH<rt_Lya_pressure_nH) return

   ! length of the central cell
   dx_loc_cm = dx_loc * scale_l
        
   u=pvar(2)/d
   v=pvar(3)/d
   w=pvar(4)/d
   ekk0=0.5d0*d*(u*u + v*v + w*w)
   eth0=pvar(5)-ekk0
#if NENER>0
   do irad=1,nener
      eth0=eth0-pvar(ndim+2+irad)
   end do
#endif
   T2  = eth0/d*scale_T2*(gamma-1)

   ekk_add = 0.0
   do idim=1,ndim

      ! Find the neighbour ; index: -x +x -y +y -z +z
      ind=0
      if(Fhat(idim)>0) ind=1
      ind_nbor=2*(idim-1) + ind + 1

      ! neighbour info is necessary to estimate the column density of neutral hydrogen
      nH_nei   = nH_nbor  (ileaf,ind_nbor)
      Zsol_nei = Zsol_nbor(ileaf,ind_nbor)
      Tk_nei   = Tk_nbor  (ileaf,ind_nbor)
      xHI_nei  = xHI_nbor (ileaf,ind_nbor)
      nH_tot   = nH_cen + nH_nei

      call estimate_boost_factor_Lya ( nH_cen, Tk_cen, Zsol_cen, xHI_cen, &
                                     & nH_nei, Tk_nei, Zsol_nei, xHI_nei, M_factor )

      ! momentum from Lya pressure
      Nph  = Nph50 * 1d50
      p_cgs = M_factor*Nph*E_Lya_ev*ev_to_erg/c_cgs  ! erg/(cm/s) = g*cm/s
      p_code = p_cgs / (scale_l**3*scale_d) / scale_v

      v_Lya = p_code / (vol_loc*d)
      scale_kms = scale_v/1d5
      if( (v_Lya*scale_kms).gt.v_limit_Lya) then
3998     format('TKWARN v= ',f9.2,1x,' M_F= ',f9.1,' lnH= ',f6.2,1x,f6.2,' lTk= ',f6.1,1x,f6.1,' lZsol=',f6.1,1x,f6.1,' lxHI=',f6.1,1x,f6.1)
      !   write(*,3998) v_Lya*scale_kms, M_factor,log10(nH_cen),log10(nH_nei), &
      !           &  log10(Tk_cen),log10(Tk_nei),log10(Zsol_cen),log10(Zsol_nei), &
      !           &  log10(xHI_cen),log10(xHI_nei)
         p_code = (v_limit_Lya/scale_kms) * (vol_loc*d)
      endif
 
      uadd = p_code * Fhat(idim) / (vol_loc*d) 
 
      select case (idim)
         case(1) 
            u = u + uadd 
         case(2) 
            v = v + uadd
         case(3) 
            w = w + uadd
      end select

      ekk_add = ekk_add + 0.5*d*uadd*uadd

   end do ! idim

   ! new momentum and total energy 
   pvar(1) = d
   pvar(2) = d*u
   pvar(3) = d*v
   pvar(4) = d*w
   ekk = 0.5*d*(u*u + v*v + w*w)
#if NENER>0
   do irad=1,nener
      eth0=eth0+pvar(ndim+2+irad)
   end do
#endif
   pvar(5) = max(ekk+eth0,ekk0+eth0+ekk_add)

end subroutine inject_momentum_Lya_cen
!################################################################
!################################################################
!################################################################
subroutine estimate_boost_factor_Lya (nH_cen,Tk_cen,Zsol_cen,xHI_cen, &
                                   &  nH_nei,Tk_nei,Zsol_nei,xHI_nei, M_factor)
   use amr_parameters,only:dp
   use Lya_commons,only:Lya_boost_factor,rt_Lya_pressure_LVG,&
                      & rt_T_sputter,Zsol_min,xHI_min,dx_loc,scale_l
   implicit none
   real(dp),intent(in)::nH_cen,Tk_cen,Zsol_cen,xHI_cen
   real(dp),intent(in)::nH_nei,Tk_nei,Zsol_nei,xHI_nei
   real(dp),intent(out)::M_factor
   real(dp) :: Zsol_nei2, Zsol_cen2, Zsol_avg, xHI_avg, Tk_avg
   real(dp) :: NHI_cm2, NHI_cm2_LVG, Lsob, nH_tot, dx_loc_cm
   real(dp) :: xHI_cen2, xHI_nei2

   ! length of a cell 
   dx_loc_cm = dx_loc*scale_l

   Zsol_cen2 = Zsol_cen
   Zsol_nei2 = Zsol_nei
   if(Tk_cen>rt_T_sputter) Zsol_cen2 = Zsol_min
   if(Tk_nei>rt_T_sputter) Zsol_nei2 = Zsol_min

   xHI_cen2 = xHI_cen
   xHI_nei2 = xHI_nei
   if(xHI_cen2<xHI_min) xHI_cen2=0.0
   if(xHI_nei2<xHI_min) xHI_nei2=0.0

   nH_tot   = nH_cen + nH_nei
   Zsol_avg = nH_cen*Zsol_cen2 + nH_nei*Zsol_nei2
   Zsol_avg = Zsol_avg / nH_tot

   xHI_avg  = nH_cen*xHI_cen2 + nH_nei*xHI_nei2
   xHI_avg  = xHI_avg / nH_tot

   Tk_avg   = nH_cen*xHI_cen2*Tk_cen + nH_nei*xHI_nei2*Tk_nei
   if(xHI_cen+xHI_nei>0) then
      Tk_avg = Tk_avg / (nH_cen*xHI_cen2 + nH_nei*xHI_nei2)
   else
      Tk_avg = Tk_nei
   endif

   ! simple estimate for NHI_cm2
   NHI_cm2  = (nH_cen*xHI_cen2 + nH_nei*xHI_nei2)*dx_loc_cm

   NHI_cm2_LVG = 0.0
   if(rt_Lya_pressure_LVG)then 
      ! Large Velocity Gravity approximation
      ! bold approximation...but probably better than local approximation
      ! Gnedin+(09) Eq.4 
      Lsob = nH_nei/max(abs(nH_cen-nH_nei),1d-10) 
      Lsob = max(min(Lsob,100.),1.) ! For 1 pc resolution, r_GMC=100 pc...
      NHI_cm2_LVG = nH_nei*(Lsob*dx_loc_cm)*xHI_avg ! tau0/propto L/2, but /2 is omitted assumming symmetry
   endif
   NHI_cm2 = max(NHI_cm2,NHI_cm2_LVG)

   call Lya_boost_factor (NHI_cm2,Zsol_avg,Tk_avg,M_factor)

end subroutine estimate_boost_factor_Lya
!################################################################
!################################################################
!################################################################
subroutine Lya_emissivity ( nH,T2,xion,dNLya_dt,P_B_Lya)
   use amr_parameters,only:dp
   use rt_cooling_module, only:getMu 
   use rt_parameters,only:nIons,ixHI,ixHII,ixHeII,ixHeIII,isH2
   use Lya_commons,only:cmp_P_B_Lya,E_Lya_ev,ev_to_erg
   ! Compute the Lya_emissivity based on Rosdahl & Blaizot (2012) 
   ! assuming Case-B recombination
   implicit none
   real(dp),intent(in)::nH,T2
   real(dp),dimension(1:nIons),intent(in)::xion
   real(dp),intent(out):: dNLya_dt,P_B_Lya
   real(dp),save::eps_coll,eps_rec,Tk,mu
   real(dp),save::xHI,xHII,xHeII,xHeIII,n_e,n_HI,n_HII
   real(dp),save::lambda,f,alpha_B,q_1s2p,kB=1.38064852d-16

   eps_coll = 0.0
   eps_rec  = 0.0

   xHII   = min(max(xion(ixHII),  0.0),1.0)
   xHeII  = min(max(xion(ixHeII), 0.0),1.0)
   xHeIII = min(max(xion(ixHeIII),0.0),1.0)

   if(isH2) then
      xHI = xion(ixHI) 
   else
      xHI = 1.0-xHII
   endif
   xHI    = min(max(xHI,0.0),1.0)
 
   mu     = getMu(xHII,xHeII,xHeIII,T2)
   Tk     = T2 * mu
   n_e    = nH * (xHII + xHeII + xHeIII*2)
   n_HI   = nH * xHI
   n_HII  = nH * xHII

   ! Turns out this is negligible
   ! Collisional emission ; Goerdt+(10) Eq 9,10 (based on Callaway+87)
   !q_1s2p = 2.41d-6/sqrt(Tk)*(Tk/1d4)**(0.22) * exp( -E_Lya_ev*ev_to_erg / (kB*Tk)) ! [cm3/s] collisional excitation rates
   !eps_coll = q_1s2p * n_e * n_HI ! * e_Lya


   ! Recombinative emission ; RB12 Eq (5)
   call cmp_P_B_Lya ( Tk, P_B_Lya)  

   ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
   lambda = 315614./Tk
   f= 1.d0+(lambda/2.74)**0.407
   alpha_B  = 2.753d-14 * lambda**1.5 / f**2.242

   eps_rec = P_B_Lya * alpha_B * n_e * n_HII ! * e_Lya
 
   dNLya_dt = eps_coll + eps_rec ! [cm-3 s-1]

end subroutine Lya_emissivity
!################################################################
!################################################################
!################################################################
subroutine stop_mpi
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
    integer :: ierr
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#else
    stop
#endif
end subroutine stop_mpi
!################################################################
!################################################################
!################################################################
subroutine dust_to_metal_RR14 (Zsolar, D2Z_solar)
   ! This returns the dust-to-metal ratio relative to the local value
   !------------------------------------
   ! Zsolar: metallicity in solar units
   ! D2Z_solar ~ D2Z / 0.3
   !------------------------------------
   ! based on Remy-Ruyer, Madden, Galliano et al. (2014)
   ! https://arxiv.org/pdf/1312.3442.pdfR
   ! Broken power law with X_{CO,Z} case (Table 1)
   use amr_parameters,only:dp
   implicit none
   real(dp),intent(in)::Zsolar
   real(dp),intent(out)::D2Z_solar
   real(dp)::a,alphaH,b,alphaL,xt,x,Xsun
   real(dp)::G2D,y,Zsol,D2Z,D2Z_at_sun

   ! y = log (G/D)
   ! x = 12 + log10(O/H)
   ! Xsun = 8.69
   ! 
   ! y = a + alphaH * (Xsun - x) for x>xt
   ! y = b + alphaL * (Xsun - x) for x<=xt
   a      = 2.21
   alphaH = 1.00 ! MW case
   b      = 0.96 ! 0.68
   alphaL = 3.10 ! 3.08 
   xt     = 8.10 ! 7.96
   Xsun   = 8.69
   Zsol   = max(Zsolar,0.02) ! No extrapolation
   x = Xsun + log10(Zsol)
   x = max(x,7.0) ! No extrapolation

   if (x>xt)then
      y = a + alphaH * (Xsun - x)
   else
      y = b + alphaL * (Xsun - x)
   endif

   G2D = 10.0**y
   D2Z = 1.0/(0.02*Zsol*G2D)
   D2Z_at_sun = 1.0/(0.02*10.0**a)
   D2Z_solar = D2Z / D2Z_at_sun

end subroutine dust_to_metal_RR14
