!==================================================================================
!=== Hydro ic for Toomre/exponential radial density profile galaxies.           ===
!=== Sets up a galactic binary merger.                                          ===
!=== Author : D.Chapon                                                          ===
!=== Date : 2010/09/01                                                          ===
!==================================================================================
module merger_parameters!{{{
  use amr_commons
  !--------------------!
  ! Galactic merger IC !
  !--------------------!

  ! Gas disk masses, given in GMsun in namelist,
  ! then converted in user unit.
  ! !!!!!! The galaxy #1 must be the heaviest !!!!!
  real(dp)::Mgas_disk1 = 1.0D2
  real(dp)::Mgas_disk2 = 1.0D2 ! Set to 0.0 for isolated galaxy runs
  ! Galactic centers, 0-centered, given in user unit
  real(dp), dimension(3)::gal_center1 = 0.0D0
  ! Set gal_center2 to values larger than 2*Lbox for isolated galaxy runs
  real(dp), dimension(3)::gal_center2 = 0.0D0
  ! Galactic disks rotation axis
  real(dp), dimension(3)::gal_axis1= (/ 0.0D0, 0.0D0, 1.0D0 /)
  real(dp), dimension(3)::gal_axis2= (/ 0.0D0, 0.0D0, 1.0D0 /)
  ! Particle ic ascii file names for the galaxies.
  !~~~~~~~~~~~~~~~~~~~~~~~WARNING~~~~~~~~~~~~~~~~~~~~~~~!
  ! Assumed to be J=z-axis, 0-centered galaxy ic files. !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=256)::ic_part_file_gal1='ic_part1'
  character(len=256)::ic_part_file_gal2='ic_part2'
  ! Rotation curve files for the galaxies
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Those circular velocity files must contain :                              !
  ! - Column #1 : radius (in pc)                                              !
  ! - Column #2 : circular velocity (in km/s)                                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=256)::Vcirc_dat_file1='Vcirc1.dat'
  character(len=256)::Vcirc_dat_file2='Vcirc2.dat'
  ! Galactic global velocities, given in km/s in namelist, 
  ! then converted in user unit.
  real(dp), dimension(3)::Vgal1 = 0.0D0
  real(dp), dimension(3)::Vgal2 = 0.0D0
  ! Gas disk typical/truncation radius/height, given in kpc in namelist,
  ! then converted in user unit.
  real(dp)::typ_radius1 = 3.0D0
  real(dp)::typ_radius2 = 3.0D0
  real(dp)::cut_radius1 = 10.0D0
  real(dp)::cut_radius2 = 10.0D0
  real(dp)::typ_height1 = 1.5D-1
  real(dp)::typ_height2 = 1.5D-1
  real(dp)::cut_height1 = 4.0D-1
  real(dp)::cut_height2 = 4.0D-1
  ! Gas density profile : radial =>'exponential' (default) or 'Toomre'
  !                       vertical => 'exponential' or 'gaussian'
  character(len=16)::rad_profile='exponential'
  character(len=16)::z_profile='gaussian'
  ! Inter galactic gas density contrast factor
  real(dp)::IG_density_factor = 1.0D-5
  
  !~~~~~~~~~~~~~~~~~~~~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~~~~!
  ! For isolated galaxy runs :                           !
  ! --------------------------                           !
  ! - set 'Mgas_disk2' to 0.0                            !
  ! - set 'gal_center2' to values larger than Lbox*2     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
end module merger_parameters!}}}

module merger_commons

  use merger_parameters
  real(dp), dimension(:,:), allocatable::Vcirc_dat1, Vcirc_dat2

end module merger_commons

subroutine read_merger_params! {{{
  use merger_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok=.true.
  character(LEN=80)::infile

  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  real(dp)::norm_u
  logical::vcirc_file1_exists, vcirc_file2_exists
  logical::ic_part_file1_exists, ic_part_file2_exists
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/merger_params/ IG_density_factor &
       & ,gal_center1, gal_center2, Mgas_disk1, Mgas_disk2 &
       & ,typ_radius1, typ_radius2, cut_radius1, cut_radius2 &
       & ,typ_height1, typ_height2, cut_height1, cut_height2 &
       & ,rad_profile, z_profile, Vcirc_dat_file1, Vcirc_dat_file2 &
       & , ic_part_file_gal1, ic_part_file_gal2 &
       & ,gal_axis1, gal_axis2, Vgal1, Vgal2


  CALL getarg(1,infile)
  open(1,file=infile)
  rewind(1)
  read(1,NML=merger_params,END=106)
  goto 107
106 write(*,*)' You need to set up namelist &MERGER_PARAMS in parameter file'
  call clean_stop
107 continue
  close(1)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !-------------------------------------------------
  ! This section deals with the galactic merger initial conditions
  !-------------------------------------------------
  ! Gaseous disks masses : GMsun => user unit
  Mgas_disk1 = Mgas_disk1 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  Mgas_disk2 = Mgas_disk2 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  ! Galaxy global velocities : km/s => user unit
  Vgal1 = Vgal1 * 1.0D5 / scale_v
  Vgal2 = Vgal2 * 1.0D5 / scale_v
  ! Gas disk typical radii (a) : kpc => user unit
  typ_radius1 = typ_radius1 * 3.085677581282D21 / scale_l
  typ_radius2 = typ_radius2 * 3.085677581282D21 / scale_l
  ! Gas disk max radii : kpc => user unit
  cut_radius1 = cut_radius1 * 3.085677581282D21 / scale_l
  cut_radius2 = cut_radius2 * 3.085677581282D21 / scale_l
  ! Gas disk typical thicknesses (h) : kpc => user unit
  typ_height1 = typ_height1 * 3.085677581282D21 / scale_l
  typ_height2 = typ_height2 * 3.085677581282D21 / scale_l
  ! Gas disk max thicknesses(zmax) : kpc => user unit
  cut_height1 = cut_height1 * 3.085677581282D21 / scale_l
  cut_height2 = cut_height2 * 3.085677581282D21 / scale_l

  select case (rad_profile)
      case ('Toomre')
          if(myid==1) write(*,*) "Chosen hydro radial density profile :'Toomre'"
      case ('exponential')
          if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
      case default
          if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
          rad_profile='exponential'
  end select
  select case (z_profile)
      case ('gaussian')
          if(myid==1) write(*,*) "Chosen hydro vertical density profile :'gaussian'"
      case ('exponential')
          if(myid==1) write(*,*) "Chosen hydro vertical density profile :'exponential'"
      case default
          if(myid==1) write(*,*) "Chosen hydro vertical density profile :'exponential'"
          z_profile='exponential'
  end select
  if(Mgas_disk2 .GT. Mgas_disk1)then
     if(myid==1)write(*,*)'Error: The galaxy #1 must be bigger than #2'
     nml_ok=.false.
  endif
  ! Check whether velocity radial profile files exist or not.
  inquire(file=trim(Vcirc_dat_file1), exist=vcirc_file1_exists)
  if(.NOT. vcirc_file1_exists) then
     if(myid==1)write(*,*)'Error: Vcirc_dat_file1 ''', trim(Vcirc_dat_file1), ''' doesn''t exist '
     nml_ok=.false.
  end if
  inquire(file=trim(Vcirc_dat_file2), exist=vcirc_file2_exists)
  if(.NOT. vcirc_file2_exists) then
     if(myid==1)write(*,*)'Error: Vcirc_dat_file2 ''', trim(Vcirc_dat_file2), ''' doesn''t exist '
     nml_ok=.false.
  end if
  if ((vcirc_file1_exists) .AND. (vcirc_file2_exists)) then
      call read_vcirc_files()
  end if

  ! Check whether ic_part files exist or not.
  inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file_gal1), exist=ic_part_file1_exists)
  if(.NOT. ic_part_file1_exists) then
     if(myid==1)write(*,*)'Error: ic_part_file1 ''', trim(ic_part_file_gal1), ''' doesn''t exist in '''&
     & , trim(initfile(levelmin))
     nml_ok=.false.
  end if
  inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file_gal2), exist=ic_part_file2_exists)
  if(.NOT. ic_part_file2_exists) then
     if(myid==1)write(*,*)'Error: ic_part_file2 ''', trim(ic_part_file_gal2), ''' doesn''t exist in '''&
     & , trim(initfile(levelmin))
     nml_ok=.false.
  end if

  ! galactic rotation axis
  norm_u = sqrt(gal_axis1(1)**2 + gal_axis1(2)**2 + gal_axis1(3)**2)
  if(norm_u .EQ. 0.0D0) then
     if(myid==1)write(*,*)'Error: Galactic axis(1) is zero '
     nml_ok=.false.
  else
    if(norm_u .NE. 1.0D0) then
       gal_axis1 = gal_axis1 / norm_u
    end if
  end if
  norm_u = sqrt(gal_axis2(1)**2 + gal_axis2(2)**2 + gal_axis2(3)**2)
  if(norm_u .EQ. 0.0D0) then
     if(myid==1)write(*,*)'Error: Galactic axis(2) is zero '
     nml_ok=.false.
  else
    if(norm_u .NE. 1.0D0) then
       gal_axis2 = gal_axis2 / norm_u
    end if
  end if

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

end subroutine read_merger_params
! }}}


subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use merger_commons
  use rt_parameters, ONLY:iIons,nIons
  use metal_yields, only: table_Solar_Fe, table_Solar_O, table_Solar_N, table_Solar_Si, table_Solar_Ca, table_Solar_Al, table_Solar_Mg, table_Solar_Eu, table_Solar_C, table_Solar_Ne, table_Solar_S
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i, ind_gal
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::v,M,rho,dzz,zint, HH, rdisk, dpdr,dmax
  real(dp)::r, rr, rr1, rr2, abs_z
  real(dp), dimension(3)::vgal, axe_rot, xx1, xx2, xx, xx_rad
  real(dp)::rgal, sum,sum2,dmin,dmin1,dmin2,zmin,zmax,pi,tol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::M_b,az,eps
  real(dp)::rmin,rmax,a2,aa,Vcirc, HH_max
  real(dp)::rho_0_1, rho_0_2, rho_0, weight, da1, Vrot
  logical, save:: init_nml=.false.
  real(dp)::my_digit,logdmin,logdmax

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Read user-defined merger parameters from the namelist
  if (.not. init_nml) then
      call read_merger_params
      init_nml = .true.
  end if
  
  ! Loop over cells
  logdmin = -4.d0
  logdmax = 5.5d0
  do i=1,nn
        ! digitizing the 3d grid into a 1d array
        x(i,:) = ((x(i,:)/boxlen) + (1.d0 / float(2**(nlevelmax+1)))) / (1.d0/float(2**nlevelmax))
        
        my_digit = x(i,1)
        my_digit = my_digit + ((x(i,2) - 1.0) * (2**nlevelmax))
        my_digit = my_digit + ((x(i,3) - 1.0) * (2**nlevelmax)**2)
        my_digit = my_digit - 1.d0
        my_digit = my_digit / float((2**(nlevelmax))**3)

        q(i,1) =  (10.d0**( my_digit*(logdmax-logdmin) + logdmin)) / scale_nH
        ! P=rho*cs^2
        a2=1d5 / scale_T2 ! sound speed squared
        q(i,ndim+2)=a2*q(i,1)
        q(i,ndim-1:ndim+1) = 0.d0
!        if(metal)q(i,6)=z_ave*0.02 ! z_ave is in solar units
        if(metal) then 
          !q(i,imetal)=z_ave*0.02*10.**(0.5-r/cut_radius1) ! z_ave is in solar units
          q(i,imetal+0) = z_ave*table_Solar_Fe ! Iron
          q(i,imetal+1) = z_ave*table_Solar_O  ! Oxygen
          q(i,imetal+2) = z_ave*table_Solar_N  ! Nitrogen
          q(i,imetal+3) = z_ave*table_Solar_Mg ! Magnesium
          q(i,imetal+4) = z_ave*table_Solar_Ne ! Neon
          q(i,imetal+5) = z_ave*table_Solar_Si ! Silicon
          q(i,imetal+6) = z_ave*table_Solar_Ca ! Calcium
          q(i,imetal+7) = z_ave*table_Solar_C  ! Carbon
          q(i,imetal+8) = z_ave*table_Solar_S  ! Sulfur
        endif
#ifdef RT
      !   q(i,iIons) = 1.d0 - 1.d-6 !initialize everything as neutral
      !   do ivar = 2,nIons
      !      q(i,iIons + ivar -1) = 1.d-6 !initialize everything as neutral
      !   enddo

      !   ! Initialize the metal ionization fractions (start all as neutral)
      !   if (oxygen_ions) then
      !      q(i,ioxygen) = 1.d0 - (1.d-10 * (n_oxygen_ions - 1)) !initialize everything as neutral
      !      do ivar = 1,n_oxygen_ions-1
      !         q(i,ioxygen + ivar) = 1.d-10
      !      enddo
      !   endif
      !   if (nitrogen_ions) then
      !      q(i,initrogen) = 1.d0 - (1.d-10 * (n_nitrogen_ions - 1)) !initialize everything as neutral
      !      do ivar = 1,n_nitrogen_ions-1
      !         q(i,initrogen + ivar) = 1.d-10
      !      enddo
      !   endif
      !   if (carbon_ions) then
      !      q(i,icarbon) = 1.d0 - (1.d-10 * (n_carbon_ions - 1)) !initialize everything as neutral
      !      do ivar = 1,n_carbon_ions-1
      !         q(i,icarbon + ivar) = 1.d-10
      !      enddo
      !   endif
      !   if (magnesium_ions) then
      !      q(i,imagnesium) = 1.d0 - (1.d-10 * (n_magnesium_ions - 1)) !initialize everything as neutral
      !      do ivar = 1,n_magnesium_ions-1
      !         q(i,imagnesium + ivar) = 1.d-10
      !      enddo
      !   endif
      !   if (silicon_ions) then
      !    q(i,isilicon) = 1.d0 - (1.d-10 * (n_silicon_ions - 1)) !initialize everything as neutral
      !    do ivar = 1,n_silicon_ions-1
      !       q(i,isilicon + ivar) = 1.d-10
      !    enddo
      !   endif
      !   if (sulfur_ions) then
      !    q(i,isulfur) = 1.d0 - (1.d-10 * (n_sulfur_ions - 1)) !initialize everything as neutral
      !    do ivar = 1,n_sulfur_ions-1
      !       q(i,isulfur + ivar) = 1.d-10
      !    enddo
      !   endif
      !   if (iron_ions) then
      !    q(i,iiron) = 1.d0 - (1.d-10 * (n_iron_ions - 1)) !initialize everything as neutral
      !    do ivar = 1,n_iron_ions-1
      !       q(i,iiron + ivar) = 1.d-10
      !    enddo
      !   endif
      !   if (neon_ions) then
      !    q(i,ineon) = 1.d0 - (1.d-10 * (n_neon_ions - 1)) !initialize everything as neutral
      !    do ivar = 1,n_neon_ions-1
      !       q(i,ineon + ivar) = 1.d-10
      !    enddo
      !   endif

        q(i,iIons) = 1.d-6 !initialize everything as ionized
        q(i,iIons+1) = 1.d0 - 1.d-6 !initialize everything as ionized
        q(i,iIons+2) = 1.d-6 !initialize everything as ionized
        q(i,iIons+3) = 1.d0 - 1.d-6 !initialize everything as ionized

        ! Initialize the metal ionization fractions (start all as neutral)
        if (oxygen_ions) then
           do ivar = 0,n_oxygen_ions-2
              q(i,ioxygen + ivar) = 1.d-10
           enddo
           q(i,ioxygen+n_oxygen_ions-1) = 1.d0 - (1.d-10 * (n_oxygen_ions - 1)) !initialize everything as neutral
        endif
        if (nitrogen_ions) then
           do ivar = 0,n_nitrogen_ions-2
              q(i,initrogen + ivar) = 1.d-10
           enddo
           q(i,initrogen+n_nitrogen_ions-1) = 1.d0 - (1.d-10 * (n_nitrogen_ions - 1)) !initialize everything as neutral
        endif
        if (carbon_ions) then
           do ivar = 0,n_carbon_ions-2
              q(i,icarbon + ivar) = 1.d-10
           enddo
           q(i,icarbon+n_carbon_ions-1) = 1.d0 - (1.d-10 * (n_carbon_ions - 1)) !initialize everything as neutral
        endif
        if (magnesium_ions) then
           do ivar = 0,n_magnesium_ions-2
              q(i,imagnesium + ivar) = 1.d-10
           enddo
           q(i,imagnesium+n_magnesium_ions-1) = 1.d0 - (1.d-10 * (n_magnesium_ions - 1)) !initialize everything as neutral
        endif
        if (silicon_ions) then
         do ivar = 0,n_silicon_ions-2
            q(i,isilicon + ivar) = 1.d-10
         enddo
         q(i,isilicon+n_silicon_ions-1) = 1.d0 - (1.d-10 * (n_silicon_ions - 1)) !initialize everything as neutral
        endif
        if (sulfur_ions) then
         do ivar = 0,n_sulfur_ions-2
            q(i,isulfur + ivar) = 1.d-10
         enddo
         q(i,isulfur+n_sulfur_ions-1) = 1.d0 - (1.d-10 * (n_sulfur_ions - 1)) !initialize everything as neutral
        endif
        if (iron_ions) then
         do ivar = 0,n_iron_ions-2
            q(i,iiron + ivar) = 1.d-10
         enddo
         q(i,iiron+n_iron_ions-1) = 1.d0 - (1.d-10 * (n_iron_ions - 1)) !initialize everything as neutral
        endif
        if (neon_ions) then
         do ivar = 0,n_neon_ions-2
            q(i,ineon + ivar) = 1.d-10
         enddo
         q(i,ineon+n_neon_ions-1) = 1.d0 - (1.d-10 * (n_neon_ions - 1)) !initialize everything as neutral
        endif
#endif
  enddo

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum : Omega = rho * V
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif

  ! kinetic energy
  ! Total system global velocity : 0
  u(1:nn,ndim+2)=0.0D0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  ! E = Ec + P / (gamma - 1)
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do


contains
!{{{
function find_Vcirc(rayon, indice)
implicit none
real(dp), intent(in)		:: rayon
integer, intent(in)		:: indice
real(dp)					:: find_Vcirc
real(dp)					:: vitesse, rayon_bin, vitesse_old, rayon_bin_old
integer					:: k, indmax

k=2
if (indice .EQ. 1) then
	indmax = size(Vcirc_dat1,1)
	rayon_bin = Vcirc_dat1(k,1)
	rayon_bin_old = Vcirc_dat1(k-1,1)
	vitesse = Vcirc_dat1(k,2)
	vitesse_old = Vcirc_dat1(k-1,2)
else
	indmax = size(Vcirc_dat2,1)
	rayon_bin = Vcirc_dat2(k,1)
	rayon_bin_old = Vcirc_dat2(k-1,1)
	vitesse = Vcirc_dat2(k,2)
	vitesse_old = Vcirc_dat2(k-1,2)
end if
do while (rayon .GT. rayon_bin)
	if(k .GE. indmax) then
		write(*,*) "Hydro IC error : Radius out of rotation curve !!!"
		call clean_stop
	end if
	k = k + 1
	if (indice .EQ. 1) then
		vitesse_old = vitesse
		vitesse = Vcirc_dat1(k,2)
		rayon_bin_old = rayon_bin
		rayon_bin = Vcirc_dat1(k,1)
	else
		vitesse_old = vitesse
		vitesse = Vcirc_dat2(k,2)
		rayon_bin_old = rayon_bin
		rayon_bin = Vcirc_dat2(k,1)
	end if
end do

find_Vcirc = vitesse_old + (rayon - rayon_bin_old) * (vitesse - vitesse_old) / (rayon_bin - rayon_bin_old)

return

end function find_Vcirc


function vect_prod(a,b)
implicit none
real(dp), dimension(3), intent(in)::a,b
real(dp), dimension(3)::vect_prod

vect_prod(1) = a(2) * b(3) - a(3) * b(2)
vect_prod(2) = a(3) * b(1) - a(1) * b(3)
vect_prod(3) = a(1) * b(2) - a(2) * b(1)

end function vect_prod



function norm2(x)
implicit none
real(dp), dimension(3), intent(in)::x
real(dp) :: norm2

norm2 = sqrt(dot_product(x,x))

end function norm2
!}}}

end subroutine condinit




!------------------------------------------------------------------------------------- 
! Circular velocity files reading
! {{{
subroutine read_vcirc_files
  use merger_commons
  implicit none
  integer:: nvitesses, ierr, i
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  ! Galaxy #1
  nvitesses = 0
  open(unit=123, file=trim(Vcirc_dat_file1), iostat=ierr)
  do while (ierr==0)
	read(123,*,iostat=ierr)
	if(ierr==0) then
		nvitesses = nvitesses + 1  ! Number of samples
	end if
  end do
  allocate(Vcirc_dat1(nvitesses,2))
  Vcirc_dat1 = 0.0D0
  rewind(123)
  do i=1,nvitesses
	read(123,*) Vcirc_dat1(i,:)
  end do
  close(123)
  ! Unit conversion kpc -> code unit and km/s -> code unit
  Vcirc_dat1(:,1) = Vcirc_dat1(:,1) * 3.085677581282D21 / scale_l
  Vcirc_dat1(:,2) = Vcirc_dat1(:,2) * 1.0D5 / scale_v

  ! Galaxy #2
  nvitesses = 0
  open(unit=123, file=trim(Vcirc_dat_file2), iostat=ierr)
  do while (ierr==0)
	read(123,*,iostat=ierr)
	if(ierr==0) then
		nvitesses = nvitesses + 1 ! Number of samples
	end if
  end do
  allocate(Vcirc_dat2(nvitesses,2))
  Vcirc_dat2 = 0.0D0
  rewind(123)
  do i=1,nvitesses
	read(123,*) Vcirc_dat2(i,:)
  end do
  close(123)
  ! Unit conversion kpc -> code unit and km/s -> code unit
  Vcirc_dat2(:,1) = Vcirc_dat2(:,1) * 3.085677581282D21 / scale_l
  Vcirc_dat2(:,2) = Vcirc_dat2(:,2) * 1.0D5 / scale_v


end subroutine read_vcirc_files
! }}}
!--------------------------------------------------------------------------------------
