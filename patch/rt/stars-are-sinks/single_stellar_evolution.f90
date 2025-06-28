MODULE sse
   ! Module for single stellar evolution
   ! The default parameters are based on MIST v1.2
   ! And were taken from https://waps.cfa.harvard.edu/MIST/resources.html 
   use amr_parameters,only:dp
   implicit none

   PUBLIC init_sse,init_sse_lc18_SN_yields,get_mass_index,get_metal_index &
        & ,SNII_LC18_yields  

   ! Global variables related to evolutionary parameters
   integer::n_sse_metals ! Number of metallicity bins
   integer::n_sse_masses ! Number of stellar masses  
   integer::n_sse_ages   ! Number of interpolated ages
   real(dp),allocatable,dimension(:)::sse_masses ! Zero-age main-sequence masses
   real(dp),allocatable,dimension(:)::sse_metals ! Metallicities
   real(dp),allocatable,dimension(:,:)::final_ages ! Main sequence lifetimes (yrs)
   real(dp),allocatable,dimension(:,:,:,:)::all_sse_probs ! All stellar properties
   integer,dimension(1:9)::sse_met_to_ramses_met_idx = (/8,3,2,5,4,6,9,7,1/)

   ! Global variables related to SN ejecta --> note that this is a slight
   ! duplication of what is already in metal yields but this was made
   ! separate for generalization
   integer,parameter::N_mass_lc18=9  ! Number of mass bins for limongi and chieffi
   integer,parameter::N_rot_lc18=3   ! Number of roration bins for limongi and chieffi
   integer,parameter::N_metal_lc18=4 ! Number of metal bins for limongi and chieffi
   ! Note that this contains the SN only yields
   ! In practice we use Set R and subtract the winds yields from the explosive yields
   ! one has to do this for the isotopic yield sets because negative numbers arise
   ! otherwise. Furthermore, one cannot subtract the pre from the exp files because
   ! the explosion impacts the pre yields --> idk why this is the case but see their
   ! paper
   real(dp),allocatable,dimension(:,:,:)::H_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::He_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::C_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::N_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::O_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::Mg_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::Ne_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::Si_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::S_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::Ca_yields_lc18
   real(dp),allocatable,dimension(:,:,:)::Fe_yields_lc18
   real(dp), dimension(1:N_metal_lc18)::LC18_metal = (/-3.d0, -2.d0, -1.d0, 0.d0/)
   real(dp), dimension(1:N_mass_lc18)::LC18_mass = (/13.d0, 15.d0, 20.d0, 25.d0, 30.d0, 40.d0, 60.d0, 80.d0, 120.d0/)
CONTAINS

SUBROUTINE init_sse_lc18_SN_yields
   use amr_commons,only:myid,lc18_SN_dir
   implicit none

   integer::i,j
   logical::ok
   character(len=128)::fSNY

   if(myid==1) &
         write(*,*) 'Initializing SNyields for individual stars'


   !!! Hydrogen !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_H.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access H SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(H_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) H_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Helium !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_HE.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access He SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(He_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) He_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Carbon !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_C.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access C SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(C_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) C_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Nitrogen !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_N.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access N SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(N_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) N_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Oxygen !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_O.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access O SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(O_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) O_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Magnesium !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_MG.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access Mg SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(Mg_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) Mg_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Silicon !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_SI.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access Si SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(Si_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) Si_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Sulfur !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_S.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access S SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(S_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) S_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Neon !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_NE.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access Ne SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(Ne_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) Ne_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Iron !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_FE.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access Fe SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(Fe_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) Fe_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

   !!! Calcium !!!
   fSNY  = trim(lc18_SN_dir)//'LC18_SN_YIELDS_CA.dat'

   ! Check if file exists 
   inquire(FILE=fSNY,exist=ok)
   if(.not.ok)then
      write(*,*) 'Cannot access Ca SN Yield file '//TRIM(fSNY)
      call clean_stop
   endif

   ! Allocate the yield array
   allocate(Ca_yields_lc18(1:N_metal_lc18,1:N_rot_lc18,1:N_mass_lc18))

   ! Read in the yields 
   open(unit=880,file=fSNY,status='old',form='unformatted')
   do i = 1,N_metal_lc18 ! Loop over metals
      do j = 1,N_rot_lc18 ! Loop over rotation
         read(880) Ca_yields_lc18(i,j,:) ! Read data for all masses
      end do
   end do
   close(880)

END SUBROUTINE init_sse_lc18_SN_yields

SUBROUTINE init_sse
   use amr_commons,only:myid,sse_dir,sse_file
   implicit none

   logical::ok
   integer::dum,i,j,k
   character(len=128)::fSSE

  if(myid==1) &
        write(*,*) 'Initializing parameters for individual stars'

   ! First check to make sure that the file exists
   write(fSSE, '(a,a)') trim(sse_dir),trim(sse_file)
   inquire(FILE=fSSE, exist=ok)
   if(.not. ok)then
      if(myid.eq.1) then
         write(*,*)'Cannot access SSE file ',fSSE
      endif
      call clean_stop
   end if

   !!! Now read in the parameters !!

   ! First open the file
   open(unit=10,file=fSSE,status='old',form='unformatted')

   ! Read the number of metals, masses, ages, and a dummy
   read(10) n_sse_metals, n_sse_masses, n_sse_ages, dum
   if(myid.eq.1)write(*,*) "SSE ints",n_sse_metals, n_sse_masses, n_sse_ages

   ! Allocate and read masses
   allocate(sse_masses(n_sse_masses))
   read(10) sse_masses(:)

   ! Allocate and read metals
   allocate(sse_metals(n_sse_metals))
   read(10) sse_metals(:)

   ! Allocate and read main-sequence lifetimes
   allocate(final_ages(n_sse_metals,n_sse_masses))
   do i = 1,n_sse_metals
      read(10) final_ages(i,:)
   end do

   ! Allocate and read all other properties
   allocate(all_sse_probs(n_sse_metals,n_sse_masses,n_sse_ages,15))
   do i = 1,n_sse_metals ! Loop over metals
      do j = 1,n_sse_masses ! Loop over masses
         do k = 1,15 ! Loop over properties
            read(10) all_sse_probs(i,j,:,k)
         end do
      end do
   end do
END SUBROUTINE init_sse

FUNCTION get_mass_index(my_mass)
   ! Function that returns the index of the closest mass match
   ! in the stellar mass array
   ! --> maybe be optimized later by a binary search
   implicit none

   integer::get_mass_index
   integer::i
   real(dp)::mass_diff,tmp_mass,my_mass

   get_mass_index = -999 ! set the index to an incorrect value 
   mass_diff = 1.d9 ! initialize to a very large number

   ! Loop over the masses and find the closest index
   do i = 1,n_sse_masses
      tmp_mass = ABS(sse_masses(i) - my_mass)
      if (tmp_mass .lt. mass_diff) then
         mass_diff = tmp_mass
         get_mass_index = i
      end if
   end do 

END FUNCTION get_mass_index

FUNCTION get_metal_index(my_met)
   ! Function that returns the index of the closest mass match
   ! in the stellar metallicity array
   ! Note that this should be [Fe/H]
   ! --> maybe be optimized later by a binary search
   implicit none

   integer::get_metal_index
   integer::i
   real(dp)::met_diff,tmp_met,my_met

   get_metal_index = -999 ! set the index to an incorrect value 
   met_diff = 1.d9 ! initialize to a very large number

   ! Loop over the metals and find the closest index
   do i = 1,n_sse_metals
      tmp_met = ABS(sse_metals(i) - my_met)
      if (tmp_met .lt. met_diff) then
         met_diff = tmp_met
         get_metal_index = i
      end if
   end do

END FUNCTION get_metal_index

!---------------------------------------
FUNCTION SNII_LC18_yields(mass, met, vel, ZZZ)
!--------------------------------------
   real(dp), intent(in) :: mass, met
   integer, intent(in) :: vel, ZZZ
   real(dp) :: SNII_LC18_yields
   real(dp),dimension(1:N_metal_lc18,1:N_mass_lc18) :: loc_yield_table

   ! Select the elements we want to use by element number
   select case (ZZZ)
      case (1) ! Hydrogen
         loc_yield_table(:,:) = H_yields_lc18(:,vel,:)
      case (2) ! Helium
         loc_yield_table(:,:) = He_yields_lc18(:,vel,:)
      case (6) ! Carbon
         loc_yield_table(:,:) = C_yields_lc18(:,vel,:)
      case (7) ! Nitrogen
         loc_yield_table(:,:) = N_yields_lc18(:,vel,:)
      case (8) ! Oxygen
         loc_yield_table(:,:) = O_yields_lc18(:,vel,:)
      case (10) ! Neon
         loc_yield_table(:,:) = Ne_yields_lc18(:,vel,:)
      case (12) ! Magnesium
         loc_yield_table(:,:) = Mg_yields_lc18(:,vel,:)
      case (14) ! Silicon
         loc_yield_table(:,:) = Si_yields_lc18(:,vel,:)
      case (16) ! Sulfur
         loc_yield_table(:,:) = S_yields_lc18(:,vel,:)
      case (20) ! Calcium
         loc_yield_table(:,:) = Ca_yields_lc18(:,vel,:)
      case (26) ! Iron
         loc_yield_table(:,:) = Fe_yields_lc18(:,vel,:)
      case default
         write(*,*) "PROBLEM IN YIELDS --> CASE NOT SELECTED"
   end select

   SNII_LC18_yields = interp_yield_2d(LC18_mass, LC18_metal, loc_yield_table, mass, LOG10(met))

   ! Here we are arbitrarily rescaling the ejecta outside of mass bounds
   ! Note that the upper one should really neber be hit
   if (mass.le.LC18_mass(1)) SNII_LC18_yields = SNII_LC18_yields * (mass/LC18_mass(1))
   if (mass.gt.LC18_mass(N_mass_lc18)) SNII_LC18_yields = SNII_LC18_yields * (mass/LC18_mass(N_mass_lc18))

END FUNCTION SNII_LC18_yields

FUNCTION interp_yield_2d(table_mass, table_met, table_yield, mass, met)
   real(dp), dimension(:), intent(in) :: table_met
   real(dp), dimension(:), intent(in) :: table_mass
   real(dp), dimension(:, :), intent(in) :: table_yield
   real(dp), intent(in) :: met, mass
   real(dp) :: interp_yield_2d
   real(dp) :: min_yield=1.d-40

   integer :: nmet, nmass
   integer :: imet, imass

   real(dp) :: x1, x2, y1, y2, f11, f21, f22, f12, a1, a2, a3, a4
   nmet = size(table_met)
   nmass = size(table_mass)

   ! Find index right of metallicity.
   imet = 2
   do while ((table_met(imet) < met) .AND. (imet <= nmet - 1))
      imet = imet + 1
   end do

   ! Find index right of mass.
   imass = 2
   do while ((table_mass(imass) < mass) .AND. (imass <= nmass - 1))
      imass = imass + 1
   end do

   ! ----------------
   ! Interpolate
   ! Store values around point.
   x1 = table_met(imet - 1)
   x2 = table_met(imet)
   y1 = table_mass(imass - 1)
   y2 = table_mass(imass)

   f11 = LOG10(max(table_yield(imet - 1, imass - 1),min_yield))
   f21 = LOG10(max(table_yield(imet, imass - 1),min_yield))
   f22 = LOG10(max(table_yield(imet, imass),min_yield))
   f12 = LOG10(max(table_yield(imet - 1, imass),min_yield))

   ! Coefficients of linear system.
   a1 = f11*x2*y2/((x1 - x2)*(y1 - y2)) + f12*x2*y1/((x1 - x2)*(y2 - y1)) + &
     &  f21*x1*y2/((x1 - x2)*(y2 - y1)) + f22*x1*y1/((x1 - x2)*(y1 - y2))
   a2 = f11*y2/((x1 - x2)*(y2 - y1)) + f12*y1/((x1 - x2)*(y1 - y2)) + &
     &  f21*y2/((x1 - x2)*(y1 - y2)) + f22*y1/((x1 - x2)*(y2 - y1))
   a3 = f11*x2/((x1 - x2)*(y2 - y1)) + f12*x2/((x1 - x2)*(y1 - y2)) + &
     &  f21*x1/((x1 - x2)*(y1 - y2)) + f22*x1/((x1 - x2)*(y2 - y1))
   a4 = f11/((x1 - x2)*(y1 - y2)) + f12/((x1 - x2)*(y2 - y1)) + &
     &  f21/((x1 - x2)*(y2 - y1)) + f22/((x1 - x2)*(y1 - y2))

   ! Compute bilinear interpolation.
   interp_yield_2d = a1 + a2*met + a3*mass + a4*met*mass
   interp_yield_2d = 10**interp_yield_2d

end function interp_yield_2d

END MODULE sse
