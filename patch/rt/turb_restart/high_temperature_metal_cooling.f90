MODULE htmcool
  use amr_parameters, only: dp, htmc_cooling_dir

  ! Harley input for ion-by-ion metal cooling
  real(dp),allocatable, dimension(:)   :: met_cool_tab_temp_ht
  real(dp),allocatable, dimension(:)   :: met_cool_tab_dens_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_iron_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_oxygen_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_nitrogen_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_carbon_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_neon_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_silicon_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_sulfur_ht
  real(dp),allocatable, dimension(:,:,:) :: met_cool_tab_magnesium_ht
  real(dp),allocatable, dimension(:,:) :: met_cool_tab_cie

    contains 
    !=======================================================================
    subroutine init_htmc_metal_cooling
    !=======================================================================
    ! High temperature metal cooling as done by Harley's cloudy models
    ! Note that this does not include free-free emission from metals
    use hydro_parameters, only:nmetals
    implicit none
    integer::i,j,ntemp=61,ndens=12
    logical::ok
    character(LEN=256):: temp_file_neq, cie_file
    character(LEN=256),dimension(1:26)::metal_cooling_neq_iron
    character(LEN=256),dimension(1:8)::metal_cooling_neq_oxygen
    character(LEN=256),dimension(1:7)::metal_cooling_neq_nitrogen
    character(LEN=256),dimension(1:12)::metal_cooling_neq_magnesium
    character(LEN=256),dimension(1:10)::metal_cooling_neq_neon
    character(LEN=256),dimension(1:14)::metal_cooling_neq_silicon
    character(LEN=256),dimension(1:6)::metal_cooling_neq_carbon
    character(LEN=256),dimension(1:16)::metal_cooling_neq_sulfur
    
    ! Non equilibrium cooling files --> TODO(code) Horrendously inefficient and verbose
    ! Doesn't impact performance but should probably tidy up at some point
    metal_cooling_neq_iron(1)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_1.dat'
    metal_cooling_neq_iron(2)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_2.dat'
    metal_cooling_neq_iron(3)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_3.dat'
    metal_cooling_neq_iron(4)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_4.dat'
    metal_cooling_neq_iron(5)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_5.dat'
    metal_cooling_neq_iron(6)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_6.dat'
    metal_cooling_neq_iron(7)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_7.dat'
    metal_cooling_neq_iron(8)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_8.dat'
    metal_cooling_neq_iron(9)  = trim(htmc_cooling_dir)//'IRON/iron_allcool_9.dat'
    metal_cooling_neq_iron(10) = trim(htmc_cooling_dir)//'IRON/iron_allcool_10.dat'
    metal_cooling_neq_iron(11) = trim(htmc_cooling_dir)//'IRON/iron_allcool_11.dat'
    metal_cooling_neq_iron(12) = trim(htmc_cooling_dir)//'IRON/iron_allcool_12.dat'
    metal_cooling_neq_iron(13) = trim(htmc_cooling_dir)//'IRON/iron_allcool_13.dat'
    metal_cooling_neq_iron(14) = trim(htmc_cooling_dir)//'IRON/iron_allcool_14.dat'
    metal_cooling_neq_iron(15) = trim(htmc_cooling_dir)//'IRON/iron_allcool_15.dat'
    metal_cooling_neq_iron(16) = trim(htmc_cooling_dir)//'IRON/iron_allcool_16.dat'
    metal_cooling_neq_iron(17) = trim(htmc_cooling_dir)//'IRON/iron_allcool_17.dat'
    metal_cooling_neq_iron(18) = trim(htmc_cooling_dir)//'IRON/iron_allcool_18.dat'
    metal_cooling_neq_iron(19) = trim(htmc_cooling_dir)//'IRON/iron_allcool_19.dat'
    metal_cooling_neq_iron(20) = trim(htmc_cooling_dir)//'IRON/iron_allcool_20.dat'
    metal_cooling_neq_iron(21) = trim(htmc_cooling_dir)//'IRON/iron_allcool_21.dat'
    metal_cooling_neq_iron(22) = trim(htmc_cooling_dir)//'IRON/iron_allcool_22.dat'
    metal_cooling_neq_iron(23) = trim(htmc_cooling_dir)//'IRON/iron_allcool_23.dat'
    metal_cooling_neq_iron(24) = trim(htmc_cooling_dir)//'IRON/iron_allcool_24.dat'
    metal_cooling_neq_iron(25) = trim(htmc_cooling_dir)//'IRON/iron_allcool_25.dat'
    metal_cooling_neq_iron(26) = trim(htmc_cooling_dir)//'IRON/iron_allcool_26.dat'

    metal_cooling_neq_oxygen(1) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_1.dat'
    metal_cooling_neq_oxygen(2) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_2.dat'
    metal_cooling_neq_oxygen(3) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_3.dat'
    metal_cooling_neq_oxygen(4) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_4.dat'
    metal_cooling_neq_oxygen(5) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_5.dat'
    metal_cooling_neq_oxygen(6) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_6.dat'
    metal_cooling_neq_oxygen(7) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_7.dat'
    metal_cooling_neq_oxygen(8) = trim(htmc_cooling_dir)//'OXYGEN/oxygen_allcool_8.dat'

    metal_cooling_neq_nitrogen(1) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_1.dat'
    metal_cooling_neq_nitrogen(2) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_2.dat'
    metal_cooling_neq_nitrogen(3) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_3.dat'
    metal_cooling_neq_nitrogen(4) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_4.dat'
    metal_cooling_neq_nitrogen(5) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_5.dat'
    metal_cooling_neq_nitrogen(6) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_6.dat'
    metal_cooling_neq_nitrogen(7) = trim(htmc_cooling_dir)//'NITROGEN/nitrogen_allcool_7.dat'

    metal_cooling_neq_magnesium(1)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_1.dat'
    metal_cooling_neq_magnesium(2)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_2.dat'
    metal_cooling_neq_magnesium(3)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_3.dat'
    metal_cooling_neq_magnesium(4)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_4.dat'
    metal_cooling_neq_magnesium(5)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_5.dat'
    metal_cooling_neq_magnesium(6)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_6.dat'
    metal_cooling_neq_magnesium(7)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_7.dat'
    metal_cooling_neq_magnesium(8)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_8.dat'
    metal_cooling_neq_magnesium(9)  = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_9.dat'
    metal_cooling_neq_magnesium(10) = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_10.dat'
    metal_cooling_neq_magnesium(11) = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_11.dat'
    metal_cooling_neq_magnesium(12) = trim(htmc_cooling_dir)//'MAGNESIUM/magnesium_allcool_12.dat'

    metal_cooling_neq_neon(1) = trim(htmc_cooling_dir)//'NEON/neon_allcool_1.dat'
    metal_cooling_neq_neon(2) = trim(htmc_cooling_dir)//'NEON/neon_allcool_2.dat'
    metal_cooling_neq_neon(3) = trim(htmc_cooling_dir)//'NEON/neon_allcool_3.dat'
    metal_cooling_neq_neon(4) = trim(htmc_cooling_dir)//'NEON/neon_allcool_4.dat'
    metal_cooling_neq_neon(5) = trim(htmc_cooling_dir)//'NEON/neon_allcool_5.dat'
    metal_cooling_neq_neon(6) = trim(htmc_cooling_dir)//'NEON/neon_allcool_6.dat'
    metal_cooling_neq_neon(7) = trim(htmc_cooling_dir)//'NEON/neon_allcool_7.dat'
    metal_cooling_neq_neon(8) = trim(htmc_cooling_dir)//'NEON/neon_allcool_8.dat'
    metal_cooling_neq_neon(9) = trim(htmc_cooling_dir)//'NEON/neon_allcool_9.dat'
    metal_cooling_neq_neon(10) = trim(htmc_cooling_dir)//'NEON/neon_allcool_10.dat'

    metal_cooling_neq_silicon(1) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_1.dat'
    metal_cooling_neq_silicon(2) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_2.dat'
    metal_cooling_neq_silicon(3) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_3.dat'
    metal_cooling_neq_silicon(4) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_4.dat'
    metal_cooling_neq_silicon(5) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_5.dat'
    metal_cooling_neq_silicon(6) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_6.dat'
    metal_cooling_neq_silicon(7) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_7.dat'
    metal_cooling_neq_silicon(8) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_8.dat'
    metal_cooling_neq_silicon(9) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_9.dat'
    metal_cooling_neq_silicon(10) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_10.dat'
    metal_cooling_neq_silicon(11) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_11.dat'
    metal_cooling_neq_silicon(12) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_12.dat'
    metal_cooling_neq_silicon(13) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_13.dat'
    metal_cooling_neq_silicon(14) = trim(htmc_cooling_dir)//'SILICON/silicon_allcool_14.dat'

    metal_cooling_neq_carbon(1) = trim(htmc_cooling_dir)//'CARBON/carbon_allcool_1.dat'
    metal_cooling_neq_carbon(2) = trim(htmc_cooling_dir)//'CARBON/carbon_allcool_2.dat'
    metal_cooling_neq_carbon(3) = trim(htmc_cooling_dir)//'CARBON/carbon_allcool_3.dat'
    metal_cooling_neq_carbon(4) = trim(htmc_cooling_dir)//'CARBON/carbon_allcool_4.dat'
    metal_cooling_neq_carbon(5) = trim(htmc_cooling_dir)//'CARBON/carbon_allcool_5.dat'
    metal_cooling_neq_carbon(6) = trim(htmc_cooling_dir)//'CARBON/carbon_allcool_6.dat'

    metal_cooling_neq_sulfur(1) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_1.dat'
    metal_cooling_neq_sulfur(2) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_2.dat'
    metal_cooling_neq_sulfur(3) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_3.dat'
    metal_cooling_neq_sulfur(4) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_4.dat'
    metal_cooling_neq_sulfur(5) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_5.dat'
    metal_cooling_neq_sulfur(6) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_6.dat'
    metal_cooling_neq_sulfur(7) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_7.dat'
    metal_cooling_neq_sulfur(8) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_8.dat'
    metal_cooling_neq_sulfur(9) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_9.dat'
    metal_cooling_neq_sulfur(10) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_10.dat'
    metal_cooling_neq_sulfur(11) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_11.dat'
    metal_cooling_neq_sulfur(12) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_12.dat'
    metal_cooling_neq_sulfur(13) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_13.dat'
    metal_cooling_neq_sulfur(14) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_14.dat'
    metal_cooling_neq_sulfur(15) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_15.dat'
    metal_cooling_neq_sulfur(16) = trim(htmc_cooling_dir)//'SULPHUR/sulphur_allcool_16.dat'

    ! Check that all files exist
    do i = 1,26 ! IRON
        inquire(FILE=metal_cooling_neq_iron(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm IRON neq metal file '//TRIM(metal_cooling_neq_iron(i))
            call clean_stop
        endif
    end do

    do i = 1,8 ! OXYGEN
        inquire(FILE=metal_cooling_neq_oxygen(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm OXYGEN neq metal file '//TRIM(metal_cooling_neq_oxygen(i))
            call clean_stop
        endif
    end do

    do i = 1,7 ! NITROGEN
        inquire(FILE=metal_cooling_neq_nitrogen(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm NITROGEN neq metal file '//TRIM(metal_cooling_neq_nitrogen(i))
            call clean_stop
        endif
    end do

    do i = 1,12 ! MAGNESIUM
        inquire(FILE=metal_cooling_neq_magnesium(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm MAGNESIUM neq metal file '//TRIM(metal_cooling_neq_magnesium(i))
            call clean_stop
        endif
    end do

    do i = 1,10 ! NEON
        inquire(FILE=metal_cooling_neq_neon(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm NEON neq metal file '//TRIM(metal_cooling_neq_neon(i))
            call clean_stop
        endif
    end do

    do i = 1,14 ! SILICON
        inquire(FILE=metal_cooling_neq_silicon(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm SILICON neq metal file '//TRIM(metal_cooling_neq_silicon(i))
            call clean_stop
        endif
    end do

    do i = 1,6 ! CARBON
        inquire(FILE=metal_cooling_neq_carbon(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm CARBON neq metal file '//TRIM(metal_cooling_neq_carbon(i))
            call clean_stop
        endif
    end do

    do i = 1,16 ! SULFUR
        inquire(FILE=metal_cooling_neq_sulfur(i),exist=ok)
        if(.not.ok)then
            write(*,*) 'Cannot access htm SULFUR neq metal file '//TRIM(metal_cooling_neq_sulfur(i))
            call clean_stop
        endif
    end do

    ! Allocate the arrays
    allocate(met_cool_tab_temp_ht(1:ntemp))
    allocate(met_cool_tab_dens_ht(1:ndens))
    allocate(met_cool_tab_iron_ht(1:ntemp,1:ndens,1:26))
    allocate(met_cool_tab_oxygen_ht(1:ntemp,1:ndens,1:8))
    allocate(met_cool_tab_nitrogen_ht(1:ntemp,1:ndens,1:7))
    allocate(met_cool_tab_magnesium_ht(1:ntemp,1:ndens,1:12))
    allocate(met_cool_tab_neon_ht(1:ntemp,1:ndens,1:10))
    allocate(met_cool_tab_silicon_ht(1:ntemp,1:ndens,1:14))
    allocate(met_cool_tab_carbon_ht(1:ntemp,1:ndens,1:6))
    allocate(met_cool_tab_sulfur_ht(1:ntemp,1:ndens,1:16))

    ! Fill in the temperature and density arrays
    do i = 1,ntemp
        met_cool_tab_temp_ht(i) = 3.d0 + ((float(i)-1.d0) * 0.1d0)
    end do
    do i = 1,ndens
        met_cool_tab_dens_ht(i) = -6.d0 + (float(i)-1.d0)
    end do

    ! Read in the metals
    do i = 1,26 ! IRON
        open(unit=880,file=metal_cooling_neq_iron(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_iron_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,8 ! OXYGEN
        open(unit=880,file=metal_cooling_neq_oxygen(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_oxygen_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,7 ! NITROGEN
        open(unit=880,file=metal_cooling_neq_nitrogen(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_nitrogen_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,12 ! MAGNESIUM
        open(unit=880,file=metal_cooling_neq_magnesium(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_magnesium_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,10 ! NEON
        open(unit=880,file=metal_cooling_neq_neon(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_neon_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,14 ! SILICON
        open(unit=880,file=metal_cooling_neq_silicon(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_silicon_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,6 ! CARBON
        open(unit=880,file=metal_cooling_neq_carbon(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_carbon_ht(j,:,i)
        end do
        close(880)
    end do

    do i = 1,16 ! SULFUR
        open(unit=880,file=metal_cooling_neq_sulfur(i),form='formatted')
        do j = 1,ntemp
            read(880,*) met_cool_tab_sulfur_ht(j,:,i)
        end do
        close(880)
    end do

    end subroutine init_htmc_metal_cooling

    !=======================================================================
    subroutine init_cie_data
    !=======================================================================
    ! CIE data
    use hydro_parameters, only:nmetals
    implicit none
    integer::i
    logical::ok
    character(LEN=256):: cie_file
    
    ! CIE  files
    cie_file             = trim(htmc_cooling_dir)//'CIE.dat'

    inquire(FILE=cie_file,exist=ok)
    if(.not.ok)then
        write(*,*) 'Cannot access CIE file '//TRIM(cie_file)
        call clean_stop
    endif

    ! Allocate the arrays
    allocate(met_cool_tab_cie(1:201,1:113))

    ! Read in the CIE file
    open(unit=880,file=cie_file,form='formatted')
    do i = 1,201 ! CIE data
        read(880,*) met_cool_tab_cie(i,:)
    end do
    close(880)
    ! Log the temperature
    met_cool_tab_cie(:,1) = LOG10(met_cool_tab_cie(:,1))

    end subroutine init_cie_data

    !=======================================================================
    subroutine cmp_metals_ht_harley(TK,ne,nFe,nO,nN,nMg,nNe,nSi,nC,nS,xFe,xO,xN,xMg,xNe,xSi,xC,xS,cool_tot,mega_log,metal_cool_smooth)
    !=======================================================================
    ! Compute cooling enhancement due to ions
    ! TK --> Temeprature (not T/mu)
    ! ne --> electron number density e/cc
    ! nX --> contains the number density of each species /cc
    ! xX --> contains the ion fractions for each species unitless
    use hydro_parameters, only:iron_ions,n_iron_ions,oxygen_ions,n_oxygen_ions, &
                                neon_ions,n_neon_ions,nitrogen_ions,n_nitrogen_ions, &
                                carbon_ions,n_carbon_ions,silicon_ions,n_silicon_ions, &
                                magnesium_ions,n_magnesium_ions,sulfur_ions,n_sulfur_ions
    implicit none

    ! INPUT
    real(dp)::TK,ne,nFe,nO,nN,nMg,nNe,nSi,nC,nS,cool_tot
    real(dp)::metal_cool_smooth,loc_red_fac
    real(dp), dimension(1:n_oxygen_ions):: xO
    real(dp), dimension(1:n_nitrogen_ions):: xN
    real(dp), dimension(1:n_carbon_ions):: xC
    real(dp), dimension(1:n_magnesium_ions):: xMg
    real(dp), dimension(1:n_silicon_ions):: xSi
    real(dp), dimension(1:n_sulfur_ions):: xS
    real(dp), dimension(1:n_iron_ions):: xFe
    real(dp), dimension(1:n_neon_ions):: xNe
    logical::mega_log

    ! LOCAL
    real(dp):: log10T,deltaT,frac_low,frac_high,log_cool
    real(dp):: log10ne,deltane,frac_low_ne,frac_high_ne
    real(dp):: cie_tot,frac_low_cie,frac_high_cie,cie_loc
    integer:: ntemp = 61, ndens = 12
    integer:: cie_Fe_offset = 87
    integer:: cie_O_offset = 22
    integer:: cie_N_offset = 14
    integer:: cie_C_offset = 7
    integer:: cie_Mg_offset = 42
    integer:: cie_Ne_offset = 31
    integer:: cie_Si_offset = 55
    integer:: cie_S_offset = 70
    integer:: IT_cool,IT_CIE,IT_dens
    integer:: i
    real(dp)::cool_tot_Fe,cool_tot_O,cool_tot_C,cool_tot_N
    real(dp)::cool_tot_Ne,cool_tot_S,cool_tot_Si,cool_tot_Mg
    real(dp)::min_frac_cie=1.d-8
    real(dp)::cool_tmp
    real(dp)::frac_00,frac_01,frac_10,frac_11
    real(dp)::boost_factor

    ! Initialize to zero
    cool_tot_Fe = 0.d0
    cool_tot_O  = 0.d0
    cool_tot_C  = 0.d0
    cool_tot_N  = 0.d0
    cool_tot_Ne = 0.d0
    cool_tot_S  = 0.d0
    cool_tot_Si = 0.d0
    cool_tot_Mg = 0.d0

    log10T = log10(TK)
    log10ne = log10(ne)

    ! Factor due to extrapolation in electron density
    boost_factor = 1.d0

    ! Get the index of the temperature in the ion cooling array
    if (log10T.le.met_cool_tab_temp_ht(1)) then
        IT_cool = 1
    else if (log10T.ge.met_cool_tab_temp_ht(ntemp-1)) then
        IT_cool = ntemp-1 ! always 1 less than max
    else
        IT_cool = FLOOR( (log10T - met_cool_tab_temp_ht(1)) / (met_cool_tab_temp_ht(2) - met_cool_tab_temp_ht(1)) ) + 1
    endif

    ! Get the index of the density in the ion cooling array
    if (log10ne.le.met_cool_tab_dens_ht(1)) then
        IT_dens = 1
        boost_factor = 10.d0**(log10ne-met_cool_tab_dens_ht(1))
    else if (log10ne.ge.met_cool_tab_dens_ht(ndens-1)) then
        IT_dens = ndens-1 ! always 1 less than max
        if (log10ne.ge.met_cool_tab_dens_ht(ndens)) boost_factor = 10.d0**(log10ne-met_cool_tab_dens_ht(ndens))
    else
        IT_dens = FLOOR( (log10ne - met_cool_tab_dens_ht(1)) / (met_cool_tab_dens_ht(2) - met_cool_tab_dens_ht(1)) ) + 1
    endif

    ! Get the index of the temperature in the cie array
    if (log10T.le.met_cool_tab_cie(1,1)) then
        IT_CIE = 1
    else if (log10T.ge.met_cool_tab_cie(201-1,1)) then
        IT_CIE = 201-1 ! always 1 less than max
    else
        IT_CIE = FLOOR( (log10T - met_cool_tab_cie(1,1)) / (met_cool_tab_cie(2,1) - met_cool_tab_cie(1,1)) ) + 1
    endif
    
    ! Prepare the interpolation (1D)
    deltaT  = (10.d0**met_cool_tab_temp_ht(IT_cool + 1))-(10.d0**met_cool_tab_temp_ht(IT_cool))
    frac_low = MIN(MAX((10.d0**met_cool_tab_temp_ht(IT_cool + 1) - TK) / deltaT,0.d0),1.d0)
    frac_high = 1.d0 - frac_low

    ! linear interpolation in density
    deltane  = (10.d0**met_cool_tab_dens_ht(IT_dens + 1))-(10.d0**met_cool_tab_dens_ht(IT_dens))
    frac_low_ne = MIN(MAX((10.d0**met_cool_tab_dens_ht(IT_dens + 1) - ne) / deltane,0.d0),1.d0)
    frac_high_ne = 1.d0 - frac_low_ne

    deltaT  = (10.d0**met_cool_tab_cie(IT_CIE + 1, 1))-(10.d0**met_cool_tab_cie(IT_CIE,1))
    frac_low_cie = MIN(MAX((10.d0**met_cool_tab_cie(IT_CIE + 1, 1) - TK) / deltaT,0.d0),1.d0)
    frac_high_cie = 1.d0 - frac_low_cie
    
    ! Compute the total cooling
    cool_tot = 0.d0

    ! Interpolation fractions
    frac_00 = frac_low * frac_low_ne
    frac_01 = frac_low * frac_high_ne
    frac_10 = frac_high * frac_low_ne
    frac_11 = frac_high * frac_high_ne

    ! First handle the ions that aren't corrected for CIE
    if (iron_ions) then
        do i = 1,n_iron_ions-1
            loc_red_fac = 1.d0
            if (i.le.2) loc_red_fac = metal_cool_smooth
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_Fe = cool_tot_Fe + (cool_tmp * nFe * xFe(i)) ! Now in erg / s / cm^3 # this is already multiplied by the electron density
            if (mega_log) write(*,*) "HT Fe",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (oxygen_ions) then
        do i = 1,n_oxygen_ions-1
            loc_red_fac = 1.d0
            if (i.le.1) loc_red_fac = metal_cool_smooth
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_O = cool_tot_O + (cool_tmp * nO * xO(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT O",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (nitrogen_ions) then
        do i = 1,n_nitrogen_ions-1
            loc_red_fac = 1.d0
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_N = cool_tot_N + (cool_tmp * nN * xN(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT N",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (carbon_ions) then
        do i = 1,n_carbon_ions-1
            loc_red_fac = 1.d0
            if (i.le.2) loc_red_fac = metal_cool_smooth
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_C = cool_tot_C + (cool_tmp * nC * xC(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT C",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (magnesium_ions) then
        do i = 1,n_magnesium_ions-1
            loc_red_fac = 1.d0
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_Mg = cool_tot_Mg + (cool_tmp * nMg * xMg(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT Mg",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (neon_ions) then
        do i = 1,n_neon_ions-1
            loc_red_fac = 1.d0
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_Ne = cool_tot_Ne + (cool_tmp * nNe * xNe(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT Ne",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (sulfur_ions) then
        do i = 1,n_sulfur_ions-1
            loc_red_fac = 1.d0
            if (i.le.1) loc_red_fac = metal_cool_smooth
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_S = cool_tot_S + (cool_tmp * nS * xS(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT S",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    if (silicon_ions) then
        do i = 1,n_silicon_ions-1
            loc_red_fac = 1.d0
            if (i.le.2) loc_red_fac = metal_cool_smooth
            cool_tmp = 0.d0
            cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
            cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
            cool_tot_Si = cool_tot_Si + (cool_tmp * nSi * xSi(i)) ! Now in erg / s / cm^3
            if (mega_log) write(*,*) "HT Si",i,cool_tmp,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
        end do
    endif
    
    ! Now handle the ions that we need CIE
    if (iron_ions) then
        if (n_iron_ions.lt.26) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_iron_ions,27
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Fe_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Fe_offset + i - 1))
            end do
            
            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_iron_ions,26
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Fe_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Fe_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                if (i.le.2) loc_red_fac = metal_cool_smooth
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_iron_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
                
                cool_tot_Fe = cool_tot_Fe + (cool_tmp * nFe * xFe(n_iron_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT Fe cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (oxygen_ions) then
        if (n_oxygen_ions.lt.8) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_oxygen_ions,9
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_O_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_O_offset + i - 1))
            end do

            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_oxygen_ions,8
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_O_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_O_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                if (i.le.1) loc_red_fac = metal_cool_smooth
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_oxygen_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)

                cool_tot_O = cool_tot_O + (cool_tmp * nO * xO(n_oxygen_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT O cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (nitrogen_ions) then
        if (n_nitrogen_ions.lt.7) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_nitrogen_ions,8
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_N_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_N_offset + i - 1))
            end do

            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_nitrogen_ions,7
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_N_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_N_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_nitrogen_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)

                cool_tot_N = cool_tot_N + (cool_tmp * nN * xN(n_nitrogen_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT N cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (carbon_ions) then
        if (n_carbon_ions.lt.6) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_carbon_ions,7
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_C_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_C_offset + i - 1))
            end do

            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_carbon_ions,6
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_C_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_C_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                if (i.le.2) loc_red_fac = metal_cool_smooth
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_carbon_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)

                cool_tot_C = cool_tot_C + (cool_tmp * nC * xC(n_carbon_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT C cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (magnesium_ions) then
        if (n_magnesium_ions.lt.12) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_magnesium_ions,13
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Mg_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Mg_offset + i - 1))
            end do

            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_magnesium_ions,12
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Mg_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Mg_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_magnesium_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)

                cool_tot_Mg = cool_tot_Mg + (cool_tmp * nMg * xMg(n_magnesium_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT Mg cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (neon_ions) then
        if (n_neon_ions.lt.10) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_neon_ions,11
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Ne_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Ne_offset + i - 1))
            end do

            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_neon_ions,10
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Ne_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Ne_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_neon_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)

                cool_tot_Ne = cool_tot_Ne + (cool_tmp * nNe * xNe(n_neon_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT Ne cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (silicon_ions) then
        if (n_silicon_ions.lt.14) then ! Only compute if there are states that we don't account for
    
            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_silicon_ions,15
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Si_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Si_offset + i - 1))
            end do
    
            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_silicon_ions,14
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_Si_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_Si_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                if (i.le.2) loc_red_fac = metal_cool_smooth
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_silicon_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)
    
                cool_tot_Si = cool_tot_Si + (cool_tmp * nSi * xSi(n_silicon_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT Si cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    if (sulfur_ions) then
        if (n_sulfur_ions.lt.16) then ! Only compute if there are states that we don't account for

            ! Calculate the fraction of the species that are in the unsampled ions
            cie_tot = 0.d0
            do i = n_sulfur_ions,17
                cie_tot = cie_tot + (frac_low_cie*met_cool_tab_cie(IT_CIE,cie_S_offset + i - 1))
                cie_tot = cie_tot + (frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_S_offset + i - 1))
            end do

            ! Loop again and get_cooling
            if (cie_tot.gt.min_frac_cie) then
                do i = n_sulfur_ions,16
                !First get the fractional contribution from the ion
                cie_loc = 0.d0
                cie_loc = frac_low_cie*met_cool_tab_cie(IT_CIE,cie_S_offset + i - 1)
                cie_loc = cie_loc + frac_high_cie*met_cool_tab_cie(IT_CIE+1,cie_S_offset + i - 1)
                cie_loc = cie_loc / cie_tot
                
                loc_red_fac = 1.d0
                if (i.le.1) loc_red_fac = metal_cool_smooth
                cool_tmp = 0.d0
                cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool,IT_dens,i)*frac_00*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool,IT_dens+1,i)*frac_01*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool+1,IT_dens,i)*frac_10*loc_red_fac)
                cool_tmp = cool_tmp + (met_cool_tab_sulfur_ht(IT_cool+1,IT_dens+1,i)*frac_11*loc_red_fac)

                cool_tot_S = cool_tot_S + (cool_tmp * nS * xS(n_sulfur_ions) * cie_loc) ! Now in erg / s / cm^3
                if (mega_log) write(*,*) "HT S cie",i,cool_tmp,cie_loc,IT_cool,IT_dens,frac_00,frac_01,frac_10,frac_11
                end do
            end if
        end if
    endif

    ! Sum the cooling rates for all metals
    cool_tot = cool_tot_Fe + cool_tot_O + cool_tot_N + cool_tot_C + &
                cool_tot_Ne + cool_tot_Mg + cool_tot_S + cool_tot_Si
    ! Boost for extrapolation in electron density
    cool_tot = cool_tot * boost_factor
    if (mega_log) write(*,*) "boost",cool_tot, boost_factor
    if(isnan(cool_tot))then !if NaN 
        write(*,*) ' Nan error in high_temp_metal_cooling'
        write(*,*) "interpolation fracs",frac_low,frac_high,frac_low_cie,frac_high_cie
        write(*,*) "electron density",ne
        write(*,*) "Fe cool",cool_tot_Fe
        write(*,*) "O cool",cool_tot_O
        write(*,*) "N cool",cool_tot_N
        write(*,*) "C cool",cool_tot_C
        write(*,*) "Ne cool",cool_tot_Ne
        write(*,*) "Mg cool",cool_tot_Mg
        write(*,*) "S cool",cool_tot_S
        write(*,*) "Si cool",cool_tot_Si
        write(*,*) "Total",cool_tot,TK
        call clean_stop
    endif

    end subroutine cmp_metals_ht_harley

END MODULE htmcool
