!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_frame()
  use amr_commons
  use pm_commons
  use hydro_commons
#ifdef RT
  use rt_parameters
  use rt_hydro_commons
  use rt_cooling_module, only: mO_NIST_amu, mN_NIST_amu, amu_to_g
  use cooling_module, only: X, Y
  use movie_lines
  use mpi_mod
#endif
  use particle_snapshot, only: output_particle_snapshot, load_ids_subset

  implicit none

  integer::dummy_io,info
  integer,parameter::tag=100

  character(len=5) :: istep_str
  character(len=100) :: moviedir, moviecmd, infofile, sinkfile
#ifdef SOLVERmhd
  character(len=100),dimension(0:NVAR+6) :: moviefiles
#else
  character(len=100),dimension(0:NVAR+2) :: moviefiles
#endif
  character(len=1024)::part_snapshot_file
  integer::icell,ncache,iskip,ngrid,nlevelmax_frame
  integer::ilun,nx_loc,ipout,npout,npart_out,ind,ix,iy,iz
  integer::imin,imax,jmin,jmax,ii,jj,kk,ll
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,dummy
  real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::xcen,ycen,zcen,delx,dely,delz
  real(dp)::xtmp,ytmp,ztmp,smooth,dist_camera,theta_cam,phi_cam,alpha,beta,smooth_theta
  real(dp)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame,rr
  real(dp)::xleft,xright,yleft,yright,zleft,zright,xcentre,ycentre,zcentre
  real(dp)::xxleft,xxright,yyleft,yyright,zzleft,zzright,xxcentre,yycentre,zzcentre
  real(dp)::xpf,ypf,zpf
  real(dp)::dx_frame,dy_frame,dx,dx_loc,dx_min,pers_corr
  real(dp)::dx_cell,dy_cell,dz_cell,dvol,dx_proj
  real(kind=8)::cell_value
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:ndim)::xx2
  real(kind=8),dimension(:,:,:),allocatable::data_frame,data_frame_all
  real(kind=8),dimension(:,:),allocatable::data_frame_loc,data_frame_loc_all
  real(kind=8),dimension(:,:),allocatable::dens,dens_all,vol,vol_all,dens_actual
  real(kind=4),dimension(:,:),allocatable::data_single
  real(kind=8) :: z1,z2,om0in,omLin,hubin,Lbox
  real(kind=8) :: observer(3),thetay,thetaz,theta,phi,temp,ekk
  real(kind=8) :: pi=3.14159265359
  real(dp),dimension(8)::xcube
  real(dp),dimension(8)::ycube
  real(dp),dimension(8)::zcube
  integer::igrid,jgrid,ipart,jpart,idim,icpu,ilevel,next_part,icube,iline
  integer::i,j,ig,ip,npart1
  integer::nalloc1,nalloc2
  integer::proj_ind,l,nh_temp,nw_temp
  real(dp)::minx,maxx,miny,maxy,minz,maxz,xpc,ypc,zpc,d1,d2,d3,d4,l1,l2,l3,l4
  real(kind=4)::ratio
  real(kind=8)::d_mass
  real(kind=8),dimension(1:nmetals)::depletion_factors

  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  logical::opened,cube_face

  integer,dimension(6,8)::lind = reshape((/1, 2, 3, 4, 1, 3, 2, 4,    &
                                           5, 6, 7, 8, 5, 7, 6, 8,    &
                                           1, 5, 2, 6, 1, 2, 5, 6,    &
                                           3, 7, 4, 8, 3, 4, 7, 8,    &
                                           1, 3, 5, 7, 1, 5, 3, 7,    &
                                           2, 4, 6, 8, 2, 6, 4, 8 /)  &
                                           ,shape(lind),order=(/2,1/))

  character(len=1)::temp_string

#if NVARNOADVECT>0
  ! movies for non-advected scalars
  real(kind=8),dimension(:,:,:),allocatable::noadvect_data_frame,noadvect_data_frame_all
  real(kind=8),dimension(:,:),allocatable::noadvect_data_frame_loc,noadvect_data_frame_loc_all
  character(len=100),dimension(1:NVARNOADVECT) :: noadvect_moviefiles
  real(kind=8)::cs, cooling_length, pth
  integer::irad
#endif

  integer,allocatable,dimension(:) :: ivar2imovie, irtvar2imovie, ielvar2imovie, inoadvect2imovie
  integer :: nvar_movie, nrtvar_movie, nelvar_movie, nnoadvectvar_movie

#ifdef RT
  character(len=100),dimension(1:NGROUPS) :: rt_moviefiles
  real(kind=8),dimension(:,:,:),allocatable::rt_data_frame,rt_data_frame_all
  real(kind=8),dimension(:,:),allocatable::rt_data_frame_loc,rt_data_frame_loc_all
  character(len=100),dimension(1:NEMISSIONLINES) :: el_moviefiles
  real(kind=8),dimension(:,:,:),allocatable::el_data_frame,el_data_frame_all
  real(kind=8),dimension(:,:),allocatable::el_data_frame_loc,el_data_frame_loc_all
  real(kind=8),dimension(1:NIONS)::xion
  real(kind=8)::mu,nHI,nHII,nHeII,nHeIII,ne,nO,nO2,nO3,nN,nN2,cv
  real(kind=8)::el_min=1.d-40
#endif


  nh_temp = nh_frame
  nw_temp = nw_frame

  ! Can only provide one center when centering on particles, so
  ! better doing it *before* any loops
  if(center_on_particles) then
    call compute_frame_center_from_particles(xcen, ycen, zcen)
  end if

 do proj_ind=1,LEN(trim(proj_axis))
  opened=.false.

  !inh_temp = nh_frame
  !nw_temp = nw_frame

#if NDIM > 1
  if(imov<1)imov=1
  if(imov>imovout)return

  ! Determine the filename, dir, etc
  if(myid==1)write(*,*)'Computing and dumping movie frame'

  call title(imov, istep_str)
  write(temp_string,'(I1)') proj_ind
  moviedir = 'movie'//trim(temp_string)//'/'
  moviecmd = 'mkdir -p '//trim(moviedir)
  if(.not.withoutmkdir) then
#ifdef NOSYSTEM
     if(myid==1)call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),int(O'755'),info)
#else
     if(myid==1)call system(moviecmd)
#endif
  endif

  infofile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
  if(myid==1)call output_info(infofile)

  moviefiles(0) = trim(moviedir)//'temp_'//trim(istep_str)//'.map'
  moviefiles(1) = trim(moviedir)//'dens_'//trim(istep_str)//'.map'
  moviefiles(2) = trim(moviedir)//'vx_'//trim(istep_str)//'.map'
  moviefiles(3) = trim(moviedir)//'vy_'//trim(istep_str)//'.map'
#if NDIM>2
  moviefiles(4) = trim(moviedir)//'vz_'//trim(istep_str)//'.map'
#endif
#if NDIM==2
  moviefiles(4) = trim(moviedir)//'pres_'//trim(istep_str)//'.map'
#endif
#if NDIM>2
  moviefiles(5) = trim(moviedir)//'pres_'//trim(istep_str)//'.map'
#endif
#if NVAR>5
#ifdef SOLVERmhd
  do ll=6,NVAR+3
#else
  do ll=6,NVAR
#endif
    write(dummy,'(I3.1)') ll
    moviefiles(ll) = trim(moviedir)//'var'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
 end do
#endif
#ifdef SOLVERmhd
  moviefiles(NVAR+4) = trim(moviedir)//'pmag_'//trim(istep_str)//'.map'
  moviefiles(NVAR+5) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
  moviefiles(NVAR+6) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'
#else
  moviefiles(NVAR+1) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
  moviefiles(NVAR+2) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'
#endif

#if NVARNOADVECT>0
  do ll=1,NVARNOADVECT
     write(dummy,'(I3.1)') ll
     noadvect_moviefiles(ll) = trim(moviedir)//'NOADVECT_'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
  end do
#endif

#ifdef RT
  ! Can generate mass weighted averages of cN_i for each group i
  if(rt) then
     do ll=1,NGROUPS
        write(dummy,'(I3.1)') ll
        rt_moviefiles(ll) = trim(moviedir)//'Fp'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
     end do
  endif

  ! Emission line movie files
  el_moviefiles(1) = trim(moviedir)//'Halpha_'//trim(istep_str)//'.map'
  el_moviefiles(2) = trim(moviedir)//'Hbeta_'//trim(istep_str)//'.map'
  el_moviefiles(3) = trim(moviedir)//'OIII_5007_'//trim(istep_str)//'.map'
  el_moviefiles(4) = trim(moviedir)//'OIII_4959_'//trim(istep_str)//'.map'
  el_moviefiles(5) = trim(moviedir)//'OIII_4363_'//trim(istep_str)//'.map'
  el_moviefiles(6) = trim(moviedir)//'OII_3726_'//trim(istep_str)//'.map'
  el_moviefiles(7) = trim(moviedir)//'OII_3728_'//trim(istep_str)//'.map'
  el_moviefiles(8) = trim(moviedir)//'NII_6583_'//trim(istep_str)//'.map'
  el_moviefiles(9) = trim(moviedir)//'HeII_1640_'//trim(istep_str)//'.map'
#endif

  ! sink filename
  if(sink)then
    sinkfile = trim(moviedir)//'sink_'//trim(istep_str)//'.txt'
    if(myid==1.and.proj_ind==1) call output_sink_csv(sinkfile)
  endif

  ! dump particles if required
  if (do_particle_snapshot) then
   if(proj_ind==1) then
      write(part_snapshot_file, "(a,a,i0.5)") trim(moviedir), "particle.bin", myid
      call load_ids_subset()
      call output_particle_snapshot(part_snapshot_file)
    end if
  end if

  if(levelmax_frame==0)then
     nlevelmax_frame=nlevelmax
  else if (levelmax_frame.gt.nlevelmax)then
     nlevelmax_frame=nlevelmax
  else
     nlevelmax_frame=levelmax_frame
  endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute frame boundaries
  if (.not. center_on_particles) then
    if(proj_axis(proj_ind:proj_ind).eq.'x')then
      xcen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
      ycen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
      zcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
    elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
      xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
      ycen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
      zcen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
    else
      xcen=xcentre_frame(proj_ind*4-3)+xcentre_frame(proj_ind*4-2)*aexp+xcentre_frame(proj_ind*4-1)*aexp**2+xcentre_frame(proj_ind*4)*aexp**3
      ycen=ycentre_frame(proj_ind*4-3)+ycentre_frame(proj_ind*4-2)*aexp+ycentre_frame(proj_ind*4-1)*aexp**2+ycentre_frame(proj_ind*4)*aexp**3
      zcen=zcentre_frame(proj_ind*4-3)+zcentre_frame(proj_ind*4-2)*aexp+zcentre_frame(proj_ind*4-1)*aexp**2+zcentre_frame(proj_ind*4)*aexp**3
    endif
  endif
  delx=deltax_frame(proj_ind*2-1)+deltax_frame(proj_ind*2)/aexp
  dely=deltay_frame(proj_ind*2-1)+deltay_frame(proj_ind*2)/aexp
  delz=deltaz_frame(proj_ind*2-1)+deltaz_frame(proj_ind*2)/aexp

  ! Camera properties
  if(cosmo) then
     if(tend_theta_camera(proj_ind).le.0d0) tend_theta_camera(proj_ind) = aendmov
     if(tend_phi_camera(proj_ind).le.0d0) tend_phi_camera(proj_ind) = aendmov
     theta_cam  = theta_camera(proj_ind)*pi/180.                                                                                 &
                +min(max(aexp-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*dtheta_camera(proj_ind)*pi/180./(aendmov-astartmov)
     phi_cam    = phi_camera(proj_ind)*pi/180.                                                                                   &
                +min(max(aexp-tstart_theta_camera(proj_ind),0d0),tend_phi_camera(proj_ind))*dphi_camera(proj_ind)*pi/180./(aendmov-astartmov)
  else
     if(tend_theta_camera(proj_ind).le.0d0) tend_theta_camera(proj_ind) = tendmov
     if(tend_phi_camera(proj_ind).le.0d0) tend_phi_camera(proj_ind) = tendmov
     theta_cam  = theta_camera(proj_ind)*pi/180.                                                                                 &
                +min(max(t-tstart_theta_camera(proj_ind),0d0),tend_theta_camera(proj_ind))*dtheta_camera(proj_ind)*pi/180./(tendmov-tstartmov)
     phi_cam    = phi_camera(proj_ind)*pi/180.                                                                                   &
                +min(max(t-tstart_phi_camera(proj_ind),0d0),tend_phi_camera(proj_ind))*dphi_camera(proj_ind)*pi/180./(tendmov-tstartmov)
  endif
  dist_camera   = boxlen
  if((focal_camera(proj_ind).le.0D0).or.(focal_camera(proj_ind).gt.dist_camera)) focal_camera(proj_ind) = dist_camera
#if NDIM>2
  if(myid==1) write(*,'(5A,F8.1,A,F8.1)') "Writing frame ", istep_str,' los=',proj_axis(proj_ind:proj_ind),' theta=',theta_cam*180./pi,' phi=',phi_cam*180./pi
#else
  if(myid==1) write(*,'(3A,F8.1)') "Writing frame ", istep_str,' theta=',theta_cam*180./pi
#endif
  ! Frame boundaries
  xleft_frame  = xcen-delx/2.
  xright_frame = xcen+delx/2.
  yleft_frame  = ycen-dely/2.
  yright_frame = ycen+dely/2.
  zleft_frame  = zcen-delz/2.
  zright_frame = zcen+delz/2.

  ! No cubic shader for 2D simulations
  if((ndim.eq.2).and.(shader_frame(proj_ind).eq.'cube')) shader_frame(proj_ind) = 'square'

  ! Count number of hydro variables
  j = 0
#ifdef SOLVERmhd
  allocate(ivar2imovie(0:NVAR+6))
  ivar2imovie = -999 ! harley initialization
  do i = 0, NVAR+6
#else
  allocate(ivar2imovie(0:NVAR+2))
  ivar2imovie = -999 ! harley initialization
  do i = 0, NVAR+2
#endif
    if (movie_vars(i) == 1) then
      j = j + 1
      ivar2imovie(i) = j
    end if
  end do

  nvar_movie = j

#ifdef RT
  if(rt) then
    ! Count number of rt variables
    j = 0
    allocate(irtvar2imovie(1:NGROUPS))
    irtvar2imovie = -999 ! harley initialization
    do i = 1, NGROUPS
      if (rt_movie_vars(i) == 1) then
        j = j + 1
        irtvar2imovie(i) = j
      end if
    end do
    nrtvar_movie = j

    ! Count number of emission line variables
    j = 0
    allocate(ielvar2imovie(1:NEMISSIONLINES))
    ielvar2imovie = -999 ! harley initialization
    do i = 1, NEMISSIONLINES
      if (el_movie_vars(i) == 1) then
        j = j + 1
        ielvar2imovie(i) = j
      end if
    end do

    nelvar_movie = j
  end if
#endif

#if NVARNOADVECT>0
  ! Count number of noadvection variables
  j = 0
  allocate(inoadvect2imovie(1:NVARNOADVECT))
  inoadvect2imovie = -999 ! harley initialization
  do i = 1, NVARNOADVECT
    if (noadvect_movie_vars(i) == 1) then
      j = j + 1
      inoadvect2imovie(i) = j
    end if
  end do

  nnoadvectvar_movie = j
#endif

  ! Allocate image
#ifdef SOLVERmhd
  allocate(data_frame(1:nw_frame,1:nh_frame,1:nvar_movie))
#else
  allocate(data_frame(1:nw_frame,1:nh_frame,1:nvar_movie))
#endif
#if NVARNOADVECT>0
  allocate(noadvect_data_frame(1:nw_frame,1:nh_frame,1:nnoadvectvar_movie))
  noadvect_data_frame(:,:,:) = 0d0
#endif
#ifdef RT
  if(rt) then
     allocate(rt_data_frame(1:nw_frame,1:nh_frame,1:nrtvar_movie))
     rt_data_frame(:,:,:) = 0d0
  endif

  ! emission line movies
  allocate(el_data_frame(1:nw_frame,1:nh_frame,1:nelvar_movie))
  el_data_frame(:,:,:) = 0d0
#endif
  allocate(dens(1:nw_frame,1:nh_frame))
  allocate(dens_actual(1:nw_frame,1:nh_frame))
  allocate(vol(1:nw_frame,1:nh_frame))
  data_frame=0d0
  dens=0d0
  dens_actual=0d0
  vol=0d0
  dx_frame=delx/dble(nw_frame)
  dy_frame=dely/dble(nh_frame)

  ! Loop over levels
  do ilevel=levelmin,nlevelmax_frame

     ! Mesh size at level ilevel in coarse cell units
     dx=0.5D0**ilevel

     ! Set position of cell centres relative to grid centre
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     dx_loc=dx*scale
     dx_min=0.5D0**nlevelmax*scale
     ncache=active(ilevel)%ngrid

        dx_proj = (dx_loc/2.0)*smooth_frame(proj_ind)
#if NDIM>2
        if((shader_frame(proj_ind).eq.'cube').and.(.not.perspective_camera(proj_ind)))then
           xcube = (/-dx_proj,-dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj, dx_proj, dx_proj/)
           ycube = (/-dx_proj,-dx_proj, dx_proj, dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj/)
           zcube = (/-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj/)
           do icube=1,8
              xtmp         = cos(theta_cam)*xcube(icube)+sin(theta_cam)*ycube(icube)
              ytmp         = cos(theta_cam)*ycube(icube)-sin(theta_cam)*xcube(icube)
              xcube(icube) = xtmp
              ycube(icube) = ytmp
              ytmp         = cos(phi_cam)*ycube(icube)+sin(phi_cam)*zcube(icube)
              ztmp         = cos(phi_cam)*zcube(icube)-sin(phi_cam)*ycube(icube)
              ycube(icube) = ytmp
              zcube(icube) = ztmp
           enddo
           minx = minval(xcube)
           maxx = maxval(xcube)
           miny = minval(ycube)
           maxy = maxval(ycube)
           minz = minval(zcube)
           maxz = maxval(zcube)
        endif
#endif
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
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

           ! Check if cell is to be considered
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0.or.ilevel==nlevelmax_frame
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Check if the cell intersect the domain
#if NDIM>2
                 ! Centering
                 xx(i,1) = xx(i,1)-xcen
                 xx(i,2) = xx(i,2)-ycen
                 xx(i,3) = xx(i,3)-zcen
                 ! Rotating
                 xtmp    = cos(theta_cam)*xx(i,1)+sin(theta_cam)*xx(i,2)
                 ytmp    = cos(theta_cam)*xx(i,2)-sin(theta_cam)*xx(i,1)
                 xx(i,1) = xtmp
                 xx(i,2) = ytmp
                 ytmp    = cos(phi_cam)*xx(i,2)+sin(phi_cam)*xx(i,3)
                 ztmp    = cos(phi_cam)*xx(i,3)-sin(phi_cam)*xx(i,2)
                 xx(i,2) = ytmp
                 xx(i,3) = ztmp
                 ! Perspective correction factor
                 pers_corr = 1.0
                 if(proj_axis(proj_ind:proj_ind).eq.'x')then
                   if(perspective_camera(proj_ind))then
                      if(shader_frame(proj_ind).eq.'cube')then
                         alpha  = atan(xx(i,2)/(dist_camera-xx(i,1)))
                         beta   = atan(xx(i,3)/(dist_camera-xx(i,1)))
                      endif
                      pers_corr = focal_camera(proj_ind)/(dist_camera-xx(i,1))
                      xx(i,2)   = xx(i,2)*pers_corr
                      xx(i,3)   = xx(i,3)*pers_corr
                      dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                   endif
                   xcentre = xx(i,2)+ycen
                   ycentre = xx(i,3)+zcen
                   zcentre = xx(i,1)+xcen
                 elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
                   if(perspective_camera(proj_ind))then
                      if(shader_frame(proj_ind).eq.'cube')then
                         alpha  = atan(xx(i,1)/(dist_camera-xx(i,2)))
                         beta   = atan(xx(i,3)/(dist_camera-xx(i,2)))
                      endif
                      pers_corr = focal_camera(proj_ind)/(dist_camera-xx(i,2))
                      xx(i,1)   = xx(i,1)*pers_corr
                      xx(i,3)   = xx(i,3)*pers_corr
                      dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                   endif
                   xcentre = xx(i,1)+xcen
                   ycentre = xx(i,3)+zcen
                   zcentre = xx(i,2)+ycen
                 else
                   if(perspective_camera(proj_ind))then
                      if(shader_frame(proj_ind).eq.'cube')then
                         alpha  = atan(xx(i,1)/(dist_camera-xx(i,3)))
                         beta   = atan(xx(i,2)/(dist_camera-xx(i,3)))
                      endif
                      pers_corr = focal_camera(proj_ind)/(dist_camera-xx(i,3))
                      xx(i,1)   = xx(i,1)*pers_corr
                      xx(i,2)   = xx(i,2)*pers_corr
                      dx_proj   = (dx_loc/2.0)*pers_corr*smooth_frame(proj_ind)
                   endif
                   xcentre = xx(i,1)+xcen
                   ycentre = xx(i,2)+ycen
                   zcentre = xx(i,3)+zcen
                 endif
                 ! Rotating the cube shader
                 if(shader_frame(proj_ind).eq.'cube'.and.perspective_camera(proj_ind))then
                    xcube = (/-dx_proj,-dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj, dx_proj, dx_proj/)
                    ycube = (/-dx_proj,-dx_proj, dx_proj, dx_proj,-dx_proj,-dx_proj, dx_proj, dx_proj/)
                    zcube = (/-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj,-dx_proj, dx_proj/)
                    do icube=1,8
                       xtmp         = cos(theta_cam)*xcube(icube)+sin(theta_cam)*ycube(icube)
                       ytmp         = cos(theta_cam)*ycube(icube)-sin(theta_cam)*xcube(icube)
                       xcube(icube) = xtmp
                       ycube(icube) = ytmp
                       ytmp         = cos(phi_cam)*ycube(icube)+sin(phi_cam)*zcube(icube)
                       ztmp         = cos(phi_cam)*zcube(icube)-sin(phi_cam)*ycube(icube)
                       ycube(icube) = ytmp
                       zcube(icube) = ztmp
                       ! Additional coordinate dependent rotation for perspective effect
                       xtmp         = cos(alpha)*xcube(icube)+sin(alpha)*zcube(icube)
                       ztmp         = cos(alpha)*zcube(icube)-sin(alpha)*xcube(icube)
                       xcube(icube) = xtmp
                       zcube(icube) = ztmp
                       ytmp         = cos(beta)*ycube(icube)+sin(beta)*zcube(icube)
                       ztmp         = cos(beta)*zcube(icube)-sin(beta)*ycube(icube)
                       ycube(icube) = ytmp
                       zcube(icube) = ztmp
                       pers_corr    = zcentre/(zcentre-zcube(icube))
                       xcube(icube) = xcube(icube)*pers_corr
                       ycube(icube) = ycube(icube)*pers_corr
                    enddo
                    minx = minval(xcube)
                    maxx = maxval(xcube)
                    miny = minval(ycube)
                    maxy = maxval(ycube)
                    minz = minval(zcube)
                    maxz = maxval(zcube)
                 endif
                 if(shader_frame(proj_ind).ne.'cube')then
                    minx = -dx_proj
                    maxx =  dx_proj
                    miny = -dx_proj
                    maxy =  dx_proj
                    minz = -dx_proj
                    maxz =  dx_proj
                 endif
                 xleft   = xcentre+minx
                 xright  = xcentre+maxx
                 yleft   = ycentre+miny
                 yright  = ycentre+maxy
                 zleft   = zcentre+minz
                 zright  = zcentre+maxz
                 if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                      & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                      & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
                 xx(i,1) = xx(i,1)-xcen
                 xx(i,2) = xx(i,2)-ycen
                 ! Rotating
                 xtmp    = cos(theta_cam)*xx(i,1)+sin(theta_cam)*xx(i,2)
                 ytmp    = cos(theta_cam)*xx(i,2)-sin(theta_cam)*xx(i,1)
                 xx(i,1) = xtmp
                 xx(i,2) = ytmp
                 xcentre = xx(i,1)+xcen
                 ycentre = xx(i,2)+ycen
                 xleft   = xcentre-dx_proj
                 xright  = xcentre+dx_proj
                 yleft   = ycentre-dx_proj
                 yright  = ycentre+dx_proj
                 if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                      & yright.lt.yleft_frame.or.yleft.ge.yright_frame)cycle
#endif
                 ! Compute map indices for the cell
                 if(xleft>xleft_frame)then
                    imin=min(int((xleft-xleft_frame)/dx_frame)+1,nw_frame)
                 else
                    imin=1
                 endif
                 imax=min(int((xright-xleft_frame)/dx_frame)+1,nw_frame)
                 if(yleft>yleft_frame)then
                    jmin=min(int((yleft-yleft_frame)/dy_frame)+1,nh_frame) ! change
                 else
                    jmin=1
                 endif
                 jmax=min(int((yright-yleft_frame)/dy_frame)+1,nh_frame) ! change

                 ! Fill up map with projected mass
                 do ii=imin,imax
                    ! Pixel x-axis position
                    xxleft      = xleft_frame+dble(ii-1)*dx_frame
                    xxright     = xxleft+dx_frame
                    xxcentre    = xxleft+0.5*dx_frame
                    dx_cell     = min(xxright,xright)-max(xxleft,xleft) 
                    do jj=jmin,jmax
                       ! Pixel y-axis position
                       yyleft   = yleft_frame+dble(jj-1)*dy_frame
                       yyright  = yyleft+dy_frame
                       yycentre = yyleft+0.5*dy_frame
                       dy_cell  = min(yyright,yright)-max(yyleft,yleft)
                       xpc      = xxcentre-xcentre
                       ypc      = yycentre-ycentre
                       xpc      = xpc-sign(0.5*dx_frame,xpc)
                       ypc      = ypc-sign(0.5*dx_frame,ypc)

#if NDIM>2
                       if(shader_frame(proj_ind).eq.'cube')then
                          cube_face = .false.
                          if(sqrt(xpc**2+ypc**2).gt.dx_proj*sqrt(3.0)) goto 666
                          ! Filling the 6 cube shader faces
                          do iline=1,6
                             l1 = (ycube(lind(iline,2))-ycube(lind(iline,1)))**2+(xcube(lind(iline,2))-xcube(lind(iline,1)))**2
                             l2 = (ycube(lind(iline,4))-ycube(lind(iline,3)))**2+(xcube(lind(iline,4))-xcube(lind(iline,3)))**2
                             if(l1.eq.0d0) cycle
                             if(l2.eq.0d0) cycle
                             d1 = ((ycube(lind(iline,2))-ycube(lind(iline,1)))*xpc-(xcube(lind(iline,2))-xcube(lind(iline,1)))*ypc+xcube(lind(iline,2))*ycube(lind(iline,1))-ycube(lind(iline,2))*xcube(lind(iline,1)))/l1
                             d2 = ((ycube(lind(iline,4))-ycube(lind(iline,3)))*xpc-(xcube(lind(iline,4))-xcube(lind(iline,3)))*ypc+xcube(lind(iline,4))*ycube(lind(iline,3))-ycube(lind(iline,4))*xcube(lind(iline,3)))/l2
                             if(d1.eq.-sign(d1,d2)) cube_face=.true.
                             if(.not.cube_face) cycle
                             l3 = (ycube(lind(iline,6))-ycube(lind(iline,5)))**2+(xcube(lind(iline,6))-xcube(lind(iline,5)))**2
                             l4 = (ycube(lind(iline,8))-ycube(lind(iline,7)))**2+(xcube(lind(iline,8))-xcube(lind(iline,7)))**2
                             if(l3.eq.0d0) cycle
                             if(l4.eq.0d0) cycle
                             d3 = ((ycube(lind(iline,6))-ycube(lind(iline,5)))*xpc-(xcube(lind(iline,6))-xcube(lind(iline,5)))*ypc+xcube(lind(iline,6))*ycube(lind(iline,5))-ycube(lind(iline,6))*xcube(lind(iline,5)))/l3
                             d4 = ((ycube(lind(iline,8))-ycube(lind(iline,7)))*xpc-(xcube(lind(iline,8))-xcube(lind(iline,7)))*ypc+xcube(lind(iline,8))*ycube(lind(iline,7))-ycube(lind(iline,8))*xcube(lind(iline,7)))/l4
                             ! Within the projected face?
                             if(d3.eq.sign(d3,d4)) cube_face=.false.
                             if(cube_face) exit
                          enddo
666                       continue
                       endif
#endif

                       ! Intersection volume
                       if((shader_frame(proj_ind).eq.'cube'          &
                          .and.(cube_face))                          &
                          .or.(shader_frame(proj_ind).eq.'sphere'    &
                          .and.sqrt(xpc**2+ypc**2).le.dx_proj)       &
                          .or.(shader_frame(proj_ind).eq.'square'    &
                          .and.(abs(xpc).le.dx_proj)                 &
                          .and.(abs(ypc).le.dx_proj)))then
                          ! Intersection volume
                          dvol        = dx_cell*dy_cell

                          dens(ii,jj)=dens(ii,jj)+dvol*max(uold(ind_cell(i),1),smallr) ! --> this is actually mass
                          dens_actual(ii,jj)=dens_actual(ii,jj)+max(uold(ind_cell(i),1),smallr) ! --> this is now the density
                          vol(ii,jj)=vol(ii,jj)+dvol

                          data_frame(ii,jj,ivar2imovie(1))=data_frame(ii,jj,ivar2imovie(1))+dvol*max(uold(ind_cell(i),1),smallr)**2
#ifdef SOLVERmhd
                          do kk=2,NVAR+3
#else
                          do kk=2,NVAR
#endif
                            ! Density weight ion fractions --> looks better :)
                            if (kk.ge.(imetal+NMETALS+NCO).and.kk.lt.(iions+NIONS)) then
                               if(movie_vars(kk).eq.1) data_frame(ii,jj,ivar2imovie(kk))=data_frame(ii,jj,ivar2imovie(kk))+uold(ind_cell(i),kk)*max(uold(ind_cell(i),1),smallr)
                            else
                               if(movie_vars(kk).eq.1) data_frame(ii,jj,ivar2imovie(kk))=data_frame(ii,jj,ivar2imovie(kk))+dvol*uold(ind_cell(i),kk)
                            end if
                          end do

#if NVARNOADVECT>0
                          ! Check if we want to make any movies of the non-advected scalars
                          if (SUM(noadvect_movie_vars).gt.0) then
                             ! Loop over non-advected scalars
                             do kk=1,NVARNOADVECT
                                ! Mass-weighted value along the line of sight
                                if (noadvect_movie_vars(kk).eq.1) then
                                   ! If temperature, density-weight along the line of sight
                                   if (kk.eq.temperature_save_ivar) then
                                      noadvect_data_frame(ii,jj,inoadvect2imovie(kk))=noadvect_data_frame(ii,jj,inoadvect2imovie(kk))+unoadvect(ind_cell(i),kk)*max(uold(ind_cell(i),1),smallr)
                                   ! If cooling length --> density-weight along the line of sight
                                   else if (kk.eq.cooling_time_ivar) then
                                      pth = uold(ind_cell(i),ndim+2)  
                                      pth = pth - 0.5d0*uold(ind_cell(i),2)**2/max(uold(ind_cell(i),1),smallr) ! Remove kinetic specific energy
#if NDIM > 1
                                      pth = pth - 0.5d0*uold(ind_cell(i),3)**2/max(uold(ind_cell(i),1),smallr)
#endif
#if NDIM > 2
                                      pth = pth - 0.5d0*uold(ind_cell(i),4)**2/max(uold(ind_cell(i),1),smallr)
#endif
#if NENER>0
                                      do irad=1,nener
                                         pth=pth-uold(ind_cell(i),ndim+2+irad) ! Remove non thermal energy
                                      end do
#endif
                                      ! Isothermal sound speed is (gamma -1) * Pth / rho
                                      cs  = pth * (gamma - 1.0) / max(uold(ind_cell(i),1),smallr)
                                      ! Make sure sounds speed is positive and square root it
                                      cs = max(cs, smallc**2)
                                      cs = sqrt(cs)

                                      ! First non-advected variable is the cooling time in code units
                                      ! Multiply by sound speed in code units to get the length
                                      cooling_length = cs * unoadvect(ind_cell(i),kk)

                                      ! Only store cells that are actually cooling
                                      if (cooling_length.gt.0.d0) noadvect_data_frame(ii,jj,inoadvect2imovie(kk))=noadvect_data_frame(ii,jj,inoadvect2imovie(kk))+cooling_length*max(uold(ind_cell(i),1),smallr)

                                   ! Otherwise it should be a cooling or heating rate so we want
                                   ! Decide whether we should density weight this or use surface brightness 
                                   else
                                      noadvect_data_frame(ii,jj,inoadvect2imovie(kk))=noadvect_data_frame(ii,jj,inoadvect2imovie(kk))+unoadvect(ind_cell(i),kk)*max(uold(ind_cell(i),1),smallr)
                                   end if
                                end if
                             end do
                          end if
#endif

#ifdef RT
                          if(rt) then
                             do kk=1,NGROUPS
                                if(rt_movie_vars(kk).eq.1) then
                                   rt_data_frame(ii,jj,irtvar2imovie(kk)) = rt_data_frame(ii,jj,irtvar2imovie(kk)) &
                                                        + dvol * rtuold(ind_cell(i), 1+(kk-1)*(ndim+1)) * rt_c_cgs &
                                                               * max(uold(ind_cell(i),1),smallr) ! mass-weighted
                                endif
                             end do

                             ! Emission line movies
                             if (SUM(el_movie_vars).ge.1) then
                                ! Get the ionization fractions (should be HI, HII, HeII, HeIII)
                                do kk=0,nIons-1
                                   xion(1+kk) = uold(ind_cell(i),iIons+kk)/uold(ind_cell(i),1)
                                end do

                                ! Calculate mu (don't correct for metals as it's minor)
                                mu = 1./(X*(1.+xion(2)) + 0.25*Y*(1.+xion(3) + 2.*xion(4)))

                                ! Calculate the temperature
                                ekk = 0.0d0
                                do idim=1,3
                                   ekk = ekk + 0.5 * uold(ind_cell(i),idim+1)**2 / max(uold(ind_cell(i),1),smallr)
                                enddo
                                temp = (gamma-1.0) * (uold(ind_cell(i),5) - ekk) !pressure
                                temp = max(temp/max(uold(ind_cell(i),1),smallr),smallc**2) * scale_T2 !temperature in K
                                temp = temp * mu

                                ! Calculate the electron density (no correction for metals---very minor)
                                nHI    = uold(ind_cell(i),1)*scale_nH*xion(1)
                                nHII   = uold(ind_cell(i),1)*scale_nH*xion(2)
                                !nHeI   = (Y/X)*uold(ind_cell(i),1)*scale_nH*0.25d0*MAX(1.0-xion(3)-xion(4),0.d0)
                                nHeII  = (Y/X)*uold(ind_cell(i),1)*scale_nH*0.25d0*xion(3)
                                nHeIII = (Y/X)*uold(ind_cell(i),1)*scale_nH*0.25d0*xion(4)
                                ne = nHII + nHeII + 2.d0*nHeIII

                                ! get the cell volume in cm
                                cv = dvol * (scale_l**3)

                                if(el_movie_vars(1).eq.1) el_data_frame(ii,jj,ielvar2imovie(1)) = el_data_frame(ii,jj,ielvar2imovie(1)) + get_halpha_lum(temp,ne,nHII,nHI,cv)
                                if(el_movie_vars(2).eq.1) el_data_frame(ii,jj,ielvar2imovie(2)) = el_data_frame(ii,jj,ielvar2imovie(2)) + get_hbeta_lum(temp,ne,nHII,nHI,cv)
                                if(el_movie_vars(9).eq.1) el_data_frame(ii,jj,ielvar2imovie(9)) = el_data_frame(ii,jj,ielvar2imovie(9)) + get_heii_1640_lum(temp,ne,nHeII,nHeIII,cv)
                                if (oxygen_ions) then 
                                   nO =  uold(ind_cell(i),imetal+1) * scale_d / (mO_NIST_amu * amu_to_g)

                                   if (temp.lt.T_sputter .and. simple_dust_depletion .and. nO.gt.0.d0) then
                                      d_mass = 0.0d0
                                      depletion_factors = 1.d0
                                      ! Use dummy variables for non-oxygen related quantities
                                      call get_dust_mass_and_depletion(depletion_factors,cv**(1.d0/3.d0),d_mass,uold(ind_cell(i),1)*scale_d,uold(ind_cell(i),1)*scale_nH,1.d0,nO,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0)
                                      nO = nO * depletion_factors(2)
                                   end if

                                   nO2 = nO * uold(ind_cell(i),ioxygen+1)/uold(ind_cell(i),1)
                                   nO3 = nO * uold(ind_cell(i),ioxygen+2)/uold(ind_cell(i),1)

                                   if(el_movie_vars(3).eq.1) el_data_frame(ii,jj,ielvar2imovie(3)) = el_data_frame(ii,jj,ielvar2imovie(3)) + get_OIII_5007_lum(temp,ne,nO3,cv)
                                   if(el_movie_vars(4).eq.1) el_data_frame(ii,jj,ielvar2imovie(4)) = el_data_frame(ii,jj,ielvar2imovie(4)) + get_OIII_4959_lum(temp,ne,nO3,cv)
                                   if(el_movie_vars(5).eq.1) el_data_frame(ii,jj,ielvar2imovie(5)) = el_data_frame(ii,jj,ielvar2imovie(5)) + get_OIII_4363_lum(temp,ne,nO3,cv)
                                   if(el_movie_vars(6).eq.1) el_data_frame(ii,jj,ielvar2imovie(6)) = el_data_frame(ii,jj,ielvar2imovie(6)) + get_OII_3726_lum(temp,ne,nO2,cv)
                                   if(el_movie_vars(7).eq.1) el_data_frame(ii,jj,ielvar2imovie(7)) = el_data_frame(ii,jj,ielvar2imovie(7)) + get_OII_3728_lum(temp,ne,nO2,cv)
                                endif
                                if (nitrogen_ions) then
                                   nN =  uold(ind_cell(i),imetal+2) * scale_d / (mO_NIST_amu * amu_to_g)
                                   nN2 = nN * uold(ind_cell(i),initrogen+1)/uold(ind_cell(i),1)
                                   if(el_movie_vars(8).eq.1) el_data_frame(ii,jj,ielvar2imovie(8)) = el_data_frame(ii,jj,ielvar2imovie(8)) + get_NII_6583_lum(temp,ne,nN2,cv)
                                endif
                             endif
                          endif
#endif


                          if (movie_vars(0).eq.1)then
                            !Get temperature
                            ekk=0.0d0
                            do idim=1,3
                               ekk=ekk+0.5*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1),smallr)
                            enddo
                            temp=(gamma-1.0)*(uold(ind_cell(i),5)-ekk) !pressure
                            temp=max(temp/max(uold(ind_cell(i),1),smallr),smallc**2)*scale_T2 !temperature in K

                            data_frame(ii,jj,ivar2imovie(0))=data_frame(ii,jj,ivar2imovie(0))+dvol*max(uold(ind_cell(i),1),smallr)*temp !mass weighted temperature
                          end if

#ifdef SOLVERmhd
                          if (movie_vars(NVAR+4).eq.1)then
                                  data_frame(ii,jj,ivar2imovie(NVAR+4))=data_frame(ii,jj,ivar2imovie(NVAR+4))+ dvol*0.125*(&
                                      uold(ind_cell(i),6)**2 + uold(ind_cell(i),7)**2 + uold(ind_cell(i),8)**2 &
                                      + uold(ind_cell(i),NVAR+1)**2 + uold(ind_cell(i),NVAR+2)**2 + uold(ind_cell(i),NVAR+3)**2)
                          end if
#endif
                       end if
                    end do
                 end do
              end if
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids
  end do
  ! End loop over levels

  ! Loop over particles
  do j=1,npartmax
#if NDIM>2
     xpf  = xp(j,1)-xcen
     ypf  = xp(j,2)-ycen
     zpf  = xp(j,3)-zcen
     ! Projection
     xtmp = cos(theta_cam)*xpf+sin(theta_cam)*ypf
     ytmp = cos(theta_cam)*ypf-sin(theta_cam)*xpf
     xpf  = xtmp
     ypf  = ytmp
     ytmp = cos(phi_cam)*ypf+sin(phi_cam)*zpf
     ztmp = cos(phi_cam)*zpf-sin(phi_cam)*ypf

     if(proj_axis(proj_ind:proj_ind).eq.'x')then
       xpf = ytmp
       ypf = ztmp
       zpf = xtmp
     elseif(proj_axis(proj_ind:proj_ind).eq.'y')then
       xpf = xtmp
       ypf = ztmp
       zpf = ytmp
     else
       xpf = xtmp
       ypf = ytmp
       zpf = ztmp
     endif
     if(perspective_camera(proj_ind))then
        xpf  = xpf*focal_camera(proj_ind)/(dist_camera-zpf)
        ypf  = ypf*focal_camera(proj_ind)/(dist_camera-zpf)
     endif
     xpf  = xpf+xcen
     ypf  = ypf+ycen
     zpf  = zpf+zcen

     if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
          & ypf.lt.yleft_frame.or.ypf.ge.yright_frame.or.&
          & zpf.lt.zleft_frame.or.zpf.ge.zright_frame)cycle
#else
     xpf  = xp(j,1)-xcen
     ypf  = xp(j,2)-ycen
     xtmp = cos(theta_cam)*xpf+sin(theta_cam)*xpf
     ytmp = cos(theta_cam)*ypf-sin(theta_cam)*ypf
     xpf  = xtmp
     ypf  = ytmp
     xpf  = xpf+xcen
     ypf  = ypf+xcen
     if(    xpf.lt.xleft_frame.or.xpf.ge.xright_frame.or.&
          & ypf.lt.yleft_frame.or.ypf.ge.yright_frame)cycle
#endif
     ! Compute map indices for the cell
     ii = min(int((xpf-xleft_frame)/dx_frame)+1,nw_frame)
     jj = min(int((ypf-yleft_frame)/dy_frame)+1,nh_frame)

     ! Fill up map with projected mass
#ifdef SOLVERmhd
     if(star) then
        if(is_DM(typep(j))) then
           if(mass_cut_refine.gt.0.0.and.zoom_only_frame) then
              if (mp(j).lt.mass_cut_refine) then 
                 if (movie_vars(NVAR+5).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+5))=data_frame(ii,jj,ivar2imovie(NVAR+5))+mp(j)
              end if
           else
              if (movie_vars(NVAR+5).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+5))=data_frame(ii,jj,ivar2imovie(NVAR+5))+mp(j)
           endif
        else if (is_star(typep(j))) then
           if (movie_vars(NVAR+6).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+6))=data_frame(ii,jj,ivar2imovie(NVAR+6))+mp(j)
        endif
     else 
        if(mass_cut_refine.gt.0.0.and.zoom_only_frame) then
           if (mp(j).lt.mass_cut_refine) then 
              if (movie_vars(NVAR+5).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+5))=data_frame(ii,jj,ivar2imovie(NVAR+5))+mp(j)
           end if
        else
           if (movie_vars(NVAR+5).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+5))=data_frame(ii,jj,ivar2imovie(NVAR+5))+mp(j)
        endif
     endif
#else
     if(star) then
        if(is_DM(typep(j))) then
           if (mass_cut_refine.gt.0.0.and.zoom_only_frame) then
              if (mp(j).lt.mass_cut_refine) then 
                 if (movie_vars(NVAR+1).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+1))=data_frame(ii,jj,ivar2imovie(NVAR+1))+mp(j)
              end if
           else
              if (movie_vars(NVAR+1).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+1))=data_frame(ii,jj,ivar2imovie(NVAR+1))+mp(j)
           endif
        else if (is_star(typep(j))) then
           if (movie_vars(NVAR+2).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+2))=data_frame(ii,jj,ivar2imovie(NVAR+2))+mp(j)
        endif
     else
        if (mass_cut_refine.gt.0.0.and.zoom_only_frame) then
           if (mp(j).lt.mass_cut_refine) then
              if (movie_vars(NVAR+1).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+1))=data_frame(ii,jj,ivar2imovie(NVAR+1))+mp(j)
           end if
        else
           if (movie_vars(NVAR+1).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+1))=data_frame(ii,jj,ivar2imovie(NVAR+1))+mp(j)
        endif
     endif
#endif
     end do
  ! End loop over particles

  ! Convert into mass weighted
!  do ii=1,nw_frame
!     do jj=1,nh_frame
!        data_frame(ii,jj,2)=data_frame(ii,jj,2)/data_frame(ii,jj,1)
!        data_frame(ii,jj,3)=data_frame(ii,jj,3)/data_frame(ii,jj,1)
!        if(metal)then
!        data_frame(ii,jj,4)=data_frame(ii,jj,4)/data_frame(ii,jj,1)
!        endif
!     end do
!  end do
#ifndef WITHOUTMPI
#ifdef SOLVERmhd
  allocate(data_frame_all(1:nw_frame,1:nh_frame,1:nvar_movie))
  call MPI_ALLREDUCE(data_frame,data_frame_all,nw_frame*nh_frame*nvar_movie,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  allocate(data_frame_loc(1:nw_frame,1:nh_frame))
  allocate(data_frame_loc_all(1:nw_frame,1:nh_frame))
  allocate(data_frame_all(1:nw_frame,1:nh_frame,1:nvar_movie))
  do j=1,nvar_movie
        data_frame_loc(:,:) = data_frame(:,:,j)
        data_frame_loc_all = 0.d0
        call MPI_ALLREDUCE(data_frame_loc,data_frame_loc_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        data_frame_all(:,:,j) = data_frame_loc_all(:,:)
  end do
#endif
  allocate(dens_all(1:nw_frame,1:nh_frame))
  call MPI_ALLREDUCE(dens,dens_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  allocate(vol_all(1:nw_frame,1:nh_frame))
  call MPI_ALLREDUCE(vol,vol_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  
  data_frame(:,:,:)=data_frame_all(:,:,:)
  dens(:,:)=dens_all(:,:)
  vol(:,:)=vol_all(:,:)
  
  ! Recycle dens_all array
  dens_all(:,:) = 0.d0
  call MPI_ALLREDUCE(dens_actual,dens_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  dens_actual(:,:) = dens_all(:,:)

  deallocate(data_frame_loc)
  deallocate(data_frame_loc_all)
  deallocate(data_frame_all)
  deallocate(dens_all)
  deallocate(vol_all)
#if NVARNOADVECT>0
  allocate(noadvect_data_frame_loc(1:nw_frame,1:nh_frame))
  allocate(noadvect_data_frame_loc_all(1:nw_frame,1:nh_frame))
  allocate(noadvect_data_frame_all(1:nw_frame,1:nh_frame,1:nnoadvectvar_movie))
  do j=1,nnoadvectvar_movie
        noadvect_data_frame_loc(:,:) = noadvect_data_frame(:,:,j)
        noadvect_data_frame_loc_all = 0.d0
        call MPI_ALLREDUCE(noadvect_data_frame_loc,noadvect_data_frame_loc_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        noadvect_data_frame_all(:,:,j) = noadvect_data_frame_loc_all(:,:)
  end do
  noadvect_data_frame(:,:,:) = noadvect_data_frame_all(:,:,:)
  deallocate(noadvect_data_frame_loc)
  deallocate(noadvect_data_frame_loc_all)
  deallocate(noadvect_data_frame_all)
#endif
#ifdef RT
  if(rt) then
     allocate(rt_data_frame_loc(1:nw_frame,1:nh_frame))
     allocate(rt_data_frame_loc_all(1:nw_frame,1:nh_frame))
     allocate(rt_data_frame_all(1:nw_frame,1:nh_frame,1:nrtvar_movie))
     do j=1,nrtvar_movie
        rt_data_frame_loc(:,:) = rt_data_frame(:,:,j)
        rt_data_frame_loc_all = 0d0
        call MPI_ALLREDUCE(rt_data_frame_loc,rt_data_frame_loc_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,info)
        rt_data_frame_all(:, :, j) = rt_data_frame_loc_all(:,:)
     end do
     rt_data_frame(:,:,:) = rt_data_frame_all(:,:,:)    
     deallocate(rt_data_frame_loc)
     deallocate(rt_data_frame_loc_all)
     deallocate(rt_data_frame_all)
  

     ! Now do emission line movies
     allocate(el_data_frame_loc(1:nw_frame,1:nh_frame))
     allocate(el_data_frame_loc_all(1:nw_frame,1:nh_frame))
     allocate(el_data_frame_all(1:nw_frame,1:nh_frame,1:nelvar_movie))
     do j=1,nelvar_movie
        el_data_frame_loc(:,:) = el_data_frame(:,:,j)
        el_data_frame_loc_all = 0d0
        call MPI_ALLREDUCE(el_data_frame_loc,el_data_frame_loc_all,nw_frame*nh_frame,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD,info)
        el_data_frame_all(:, :, j) = el_data_frame_loc_all(:,:)
     end do
     el_data_frame(:,:,:) = el_data_frame_all(:,:,:)
     deallocate(el_data_frame_loc)
     deallocate(el_data_frame_loc_all)
     deallocate(el_data_frame_all)

  endif

#endif
#endif
  ! Convert into mass weighted
  do ii=1,nw_frame
    do jj=1,nh_frame
      do kk=0,NVAR
#ifdef SOLVERmhd
        if(kk==6.or.kk==7.or.kk==8) cycle
#endif
        if (kk.ge.(imetal+NMETALS+NCO).and.kk.lt.(iions+NIONS)) then
           if(movie_vars(kk).eq.1) data_frame(ii,jj,ivar2imovie(kk))=data_frame(ii,jj,ivar2imovie(kk))/dens_actual(ii,jj)
        else
           if(movie_vars(kk).eq.1) data_frame(ii,jj,ivar2imovie(kk))=data_frame(ii,jj,ivar2imovie(kk))/dens(ii,jj)
        end if
      end do
#ifdef SOLVERmhd
      if(movie_vars(NVAR+4).eq.1) data_frame(ii,jj,ivar2imovie(NVAR+4))=data_frame(ii,jj,ivar2imovie(NVAR+4))/vol(ii,jj)
#endif
#if NVARNOADVECT>0
      do kk=1,NVARNOADVECT
        if(noadvect_movie_vars(kk).eq.1) noadvect_data_frame(ii,jj,inoadvect2imovie(kk))=noadvect_data_frame(ii,jj,inoadvect2imovie(kk))/dens_actual(ii,jj) !--> density weighted
      end do
#endif
#ifdef RT
      if(rt) then
         do kk=1,NGROUPS
            if(rt_movie_vars(kk).eq.1) &
                 rt_data_frame(ii,jj,irtvar2imovie(kk))=rt_data_frame(ii,jj,irtvar2imovie(kk))/dens(ii,jj)
         end do
      endif
#endif

    end do
  end do
  deallocate(dens)
  deallocate(dens_actual)
  deallocate(vol)

  if(myid==1)then
     ilun=10
     allocate(data_single(1:nw_frame,1:nh_frame))
     ! Output mass weighted density
#ifdef SOLVERmhd
     do kk=0, NVAR+6
#else
     do kk=0, NVAR+2
#endif
       if (movie_vars(kk).eq.1)then
         open(ilun,file=TRIM(moviefiles(kk)),form='unformatted')
         data_single=data_frame(:,:,ivar2imovie(kk))
         rewind(ilun)
         if(tendmov>0)then
            write(ilun)t,delx,dely,delz
         else
            write(ilun)aexp,delx,dely,delz
         endif
         write(ilun)nw_frame,nh_frame
         write(ilun)data_single
         close(ilun)
       end if
     end do
!     ! Output mass weighted temperature
!     open(ilun,file=TRIM(moviefiles(0)),form='unformatted')
!     data_single=data_frame(:,:,0)
!!     write(*,*) 'testing', data_single(100,100)
!     rewind(ilun)
!     if(tendmov>0)then
!        write(ilun)t,delx,dely,delz
!     else
!        write(ilun)aexp,delx,dely,delz
!     endif
!     write(ilun)nw_frame,nh_frame
!     write(ilun)data_single
!     close(ilun)
!     ! Output mass weighted metal fraction
!     if(metal)then
!        open(ilun,file=TRIM(moviefiles(6)),form='unformatted')
!        data_single=data_frame(:,:,6)
!        rewind(ilun)
!        if(tendmov>0)then
!           write(ilun)t,delx,dely,delz
!        else
!           write(ilun)aexp,delx,dely,delz
!        endif
!        write(ilun)nw_frame,nh_frame
!        write(ilun)data_single
!        close(ilun)
!     endif
#if NVARNOADVECT>0
     do kk=1, NVARNOADVECT
       if (noadvect_movie_vars(kk).eq.1)then
         open(ilun,file=TRIM(noadvect_moviefiles(kk)),form='unformatted')
         data_single(:,:)=0.
         data_single=noadvect_data_frame(:,:,inoadvect2imovie(kk))
         rewind(ilun)
         if(tendmov>0)then
            write(ilun)t,delx,dely,delz
         else
            write(ilun)aexp,delx,dely,delz
         endif
         write(ilun)nw_frame,nh_frame
         write(ilun)data_single
         close(ilun)
       end if
     end do
#endif
#ifdef RT
      if(rt) then
         do kk=1, NGROUPS
            if (rt_movie_vars(kk).eq.1) then
               open(ilun,file=TRIM(rt_moviefiles(kk)),form='unformatted')
               data_single(:,:)=0.
               data_single=rt_data_frame(:,:,irtvar2imovie(kk))
               rewind(ilun)  
               if(tendmov>0)then
                  write(ilun)t,delx,dely,delz
               else
                  write(ilun)aexp,delx,dely,delz
               endif
               write(ilun)nw_frame,nh_frame
               write(ilun)data_single
               close(ilun)
            end if
         end do
      endif
      ! Write the emission line movie files
      do kk=1, NEMISSIONLINES
         if (el_movie_vars(kk).eq.1) then
            open(ilun,file=TRIM(el_moviefiles(kk)),form='unformatted')
            data_single(:,:)=0.
            data_single=log10(MAX(el_data_frame(:,:,ielvar2imovie(kk)),el_min))
            rewind(ilun)  
            if(tendmov>0)then
              write(ilun)t,delx,dely,delz
            else
              write(ilun)aexp,delx,dely,delz
            endif
            write(ilun)nw_frame,nh_frame
            write(ilun)data_single
            close(ilun)
         end if
      end do
#endif

     deallocate(data_single)
  endif

  deallocate(data_frame)
  deallocate(ivar2imovie)
#if NVARNOADVECT>0
  deallocate(noadvect_data_frame)
  deallocate(inoadvect2imovie)
#endif

#ifdef RT
  if(rt) deallocate(rt_data_frame)
  deallocate(el_data_frame)
  deallocate(irtvar2imovie)
  deallocate(ielvar2imovie)
#endif
#endif
  ! Update counter
  if(proj_ind.eq.len(trim(proj_axis)))imov=imov+1

  nw_frame = nw_temp
  nh_frame = nh_temp
 enddo
end subroutine output_frame

subroutine set_movie_vars()
  use amr_commons
  ! This routine sets the movie vars from textual form
  integer::ll
  character(LEN=5)::dummy

  if(ANY(movie_vars_txt=='temp ')) movie_vars(0)=1
  if(ANY(movie_vars_txt=='dens ')) movie_vars(1)=1
  if(ANY(movie_vars_txt=='vx   ')) movie_vars(2)=1
  if(ANY(movie_vars_txt=='vy   ')) movie_vars(3)=1
#if NDIM>2
  if(ANY(movie_vars_txt=='vz   ')) movie_vars(4)=1
#endif
#if NDIM==2
  if(ANY(movie_vars_txt=='pres ')) movie_vars(4)=1
#endif
#if NDIM>2
  if(ANY(movie_vars_txt=='pres ')) movie_vars(5)=1
#endif
#if NVAR>5
  do ll=6,NVAR
    write(dummy,'(I3.1)') ll
    if(ANY(movie_vars_txt=='var'//trim(adjustl(dummy))//' ')) movie_vars(ll)=1
 end do
#endif
#ifdef SOLVERmhd
  if(ANY(movie_vars_txt=='bxl  ')) movie_vars(6)=1
  if(ANY(movie_vars_txt=='byl  ')) movie_vars(7)=1
  if(ANY(movie_vars_txt=='bzl  ')) movie_vars(8)=1
  if(ANY(movie_vars_txt=='bxr  ')) movie_vars(NVAR+1)=1
  if(ANY(movie_vars_txt=='byr  ')) movie_vars(NVAR+2)=1
  if(ANY(movie_vars_txt=='bzr  ')) movie_vars(NVAR+3)=1
  if(ANY(movie_vars_txt=='pmag ')) movie_vars(NVAR+4)=1
  if(ANY(movie_vars_txt=='dm   ')) movie_vars(NVAR+5)=1
  if(ANY(movie_vars_txt=='stars')) movie_vars(NVAR+6)=1
#else
  if(ANY(movie_vars_txt=='dm   ')) movie_vars(NVAR+1)=1
  if(ANY(movie_vars_txt=='stars')) movie_vars(NVAR+2)=1
#endif
end subroutine set_movie_vars

subroutine compute_frame_center_from_particles(xcentre, ycentre, zcentre)
   use amr_commons, only: dp, ndim, myid
   use amr_parameters, only: i8b, center_on_particles_file, boxlen, movie_particle_ids
   use pm_commons, only: xp, mp, idp, typep, part_t, npartmax
   use mpi_mod

   use particle_snapshot, only: read_particle_ids, binary_search
   implicit none

   real(dp), intent(out) :: xcentre, ycentre, zcentre

   integer :: ipart, index, idim
   integer, save :: particle_type, n

   logical :: ok

   real(dp), dimension(1:ndim) :: xx, xx_first, xx_all
   real(dp), dimension(:, :), allocatable :: xx_found, xx_found_all
   real(dp), dimension(:), allocatable :: m_found_all, m_found
   real(dp) :: m_tot
#ifndef WITHOUTMPI
   integer :: ierr
#endif

   if (.not. allocated(movie_particle_ids)) then
      call read_particle_ids(center_on_particles_file, n, movie_particle_ids, particle_type)
   end if

   allocate(xx_found(1:n, 1:ndim), xx_found_all(1:n, 1:ndim))
   allocate(m_found(1:n), m_found_all(1:n))
   xx_found(:, :) = 0
   xx_found_all(:, :) = 0
   m_found(:) = 0
   m_found_all(:) = 0

   ! Could be optimized by using the linked list   
   do ipart = 1, npartmax
      if (typep(ipart)%family == particle_type) then
         index = binary_search(idp(ipart), movie_particle_ids, n)
         if (index > 0) then
            xx_found(index, :) = xp(ipart, :)
            m_found(index) = mp(ipart)
         end if
      end if
   end do

#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(xx_found, xx_found_all, n*ndim, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(m_found, m_found_all, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
   xx_found_all(:, :) = xx_found(:, :)
   m_found_all(:) = m_found(:)
#endif

   ok = .false.
   do ipart = 1, n
      if (any(xx_found_all(ipart, :) /= 0)) then
         ok = .true.
         exit
      end if
   end do
   if (.not. ok) then
      if (myid == 1) write(*, *) "WARNING: No particle from selection found. Centering the movie on center of box"
      xcentre = boxlen / 2
      ycentre = boxlen / 2
      zcentre = boxlen / 2
      return
   end if

   ! Find center of mass with respect to first non-null particle
   m_tot = 0
   xx_all(:) = 0
   do ipart = 1, n
      if (any(xx_found_all(ipart, :) > 0)) then
         if (m_tot == 0) then
            xx_first(:) = xx_found_all(ipart, :)
         else
            do idim = 1, ndim
               if (xx_found_all(ipart, idim) - xx_first(idim) > boxlen/2) then
                  xx_found_all(ipart, idim) = xx_found_all(ipart, idim) - boxlen
               else if (xx_found_all(ipart, idim) - xx_first(idim) < -boxlen/2) then
                  xx_found_all(ipart, idim) = xx_found_all(ipart, idim) + boxlen
               end if
               xx_all(idim) = xx_all(idim) + (xx_found_all(ipart, idim) - xx_first(idim)) * m_found_all(ipart)
            end do
         end if
         m_tot = m_tot + m_found_all(ipart)
      end if
   end do

   xx_all(:) = xx_all(:) / m_tot + xx_first(:)

   ! Take into account periodic boundaries
   do idim = 1, ndim
      if (xx_all(idim) > boxlen) then
         xx_all(idim) = xx_all(idim) - boxlen
      else if (xx_all(idim) < 0) then
         xx_all(idim) = xx_all(idim) + boxlen
      end if
   end do

   xcentre = xx_all(1)
   if (ndim > 1) ycentre = xx_all(2)
   if (ndim > 2) zcentre = xx_all(3)

end subroutine compute_frame_center_from_particles

