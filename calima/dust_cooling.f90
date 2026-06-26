module dust_cooling
    use constants, only: amu2g, kB, mC_amu, mO_amu,pi, mH, ln10
    use dust_commons
    use dust_charging
    use dust_utils
#ifdef RTZ
    use rtz_module, only:elements
#endif

    implicit none

    private

    public :: compute_dust_coll_heating_BH80
    public :: compute_dust_coll_heating
    public :: init_dust_coll_heating_BH80_cache

    logical, save :: bh80_cache_ready = .false.
    real(dp), allocatable, save :: bh80_h2_prefactor(:)
    real(dp), allocatable, save :: bh80_co_prefactor(:)
    real(dp), allocatable, save :: bh80_h2_accomm_zero(:)
    real(dp), allocatable, save :: bh80_h_accomm_zero(:)
    real(dp), allocatable, save :: bh80_species_prefactor(:,:)
    real(dp), public, save :: time_coll_heating = 0d0

    contains

    subroutine init_dust_coll_heating_BH80_cache
        implicit none

        integer :: iel, i
        real(dp) :: agrain, mproj

        if (bh80_cache_ready) return

        if (allocated(bh80_h2_prefactor)) deallocate(bh80_h2_prefactor)
        if (allocated(bh80_co_prefactor)) deallocate(bh80_co_prefactor)
        if (allocated(bh80_h2_accomm_zero)) deallocate(bh80_h2_accomm_zero)
        if (allocated(bh80_h_accomm_zero)) deallocate(bh80_h_accomm_zero)
        if (allocated(bh80_species_prefactor)) deallocate(bh80_species_prefactor)

        allocate(bh80_h2_prefactor(1:ndust))
        allocate(bh80_co_prefactor(1:ndust))
        allocate(bh80_h2_accomm_zero(1:ndust))
        allocate(bh80_h_accomm_zero(1:ndust))
        allocate(bh80_species_prefactor(1:ndust,1:n_elements))

        bh80_h2_prefactor(:) = 0d0
        bh80_co_prefactor(:) = 0d0
        bh80_h2_accomm_zero(:) = 0d0
        bh80_h_accomm_zero(:) = 0d0
        bh80_species_prefactor(:,:) = 0d0

        do i = 1, ndust
            agrain = dustbins_props(i)%asize_cm
            bh80_h2_prefactor(i) = sqrt(8d0*kB/(pi*2d0*mH)) * pi * agrain**2d0
            bh80_co_prefactor(i) = sqrt(8d0*kB/(pi*(mC_amu+mO_amu)*amu2g)) * pi * agrain**2d0
            ! BH83 accommodation coefficient: alpha_0 = 4*m_proj*m_grain/(m_proj+m_grain)^2
            bh80_h2_accomm_zero(i) = 8d0 * mH * dustbins_props(i)%mgrain / (2d0*mH + dustbins_props(i)%mgrain)**2d0
            bh80_h_accomm_zero(i) = 4d0 * mH * dustbins_props(i)%mgrain / (mH + dustbins_props(i)%mgrain)**2d0

            do iel = 1, n_elements
#ifdef RTZ
                mproj = elements(iel)%atomic_mass * amu2g
#else
                mproj = el_atomic_masses_amu(iel) * amu2g
#endif
                if (mproj <= 0d0) cycle
                bh80_species_prefactor(i,iel) = sqrt(8d0*kB/(pi*mproj)) * pi * agrain**2d0
            end do
        end do

        bh80_cache_ready = .true.
    end subroutine init_dust_coll_heating_BH80_cache

    subroutine compute_dust_coll_heating_BH80(i_dust,nElement,xelem_ions,nH2,nCO,Tgas,Td,Hcoll)
        implicit none

        integer, intent(in) :: i_dust
        real(dp), intent(in) :: nH2,nCO
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), intent(in) :: Tgas
        real(dp), intent(in) :: Td

        real(dp), intent(inout) :: Hcoll

        integer :: j,iel,nions_loc
        real(dp) :: accomm_factor,prefactor,accomm_factor_zero,xion

        if (.not. bh80_cache_ready) call init_dust_coll_heating_BH80_cache

        ! 1. Compute common prefactor
        prefactor = 2d0 * kB * (Tgas - Td) * sqrt(Tgas)

        ! 2. Compute contribution from species whose accommodation factor is unity
        ! (everything except H2 and H neutral)
        Hcoll = 0d0
        do iel = 1, n_elements
#ifdef RTZ
            nions_loc = max(1, elements(iel)%n_ions)
#else
            nions_loc = 1
#endif
            if (nElement(iel) <= 1d-10) cycle

            do j = 1, nions_loc
                xion = xelem_ions(iel,j)
                if (xion <= 1d-10) cycle

                ! Skip H neutral here as it depends on Td
                if (iel == 1 .and. j == 1) cycle

                Hcoll = Hcoll + nElement(iel) * xion * bh80_species_prefactor(i_dust,iel)
            end do
        end do

        ! 3. Add H2 contribution (depends on Td)
        accomm_factor_zero = bh80_h2_accomm_zero(i_dust)
        accomm_factor = (1d0 - accomm_factor_zero) * exp(-sqrt((Td+Tgas)/250d0)) + accomm_factor_zero
        Hcoll = Hcoll + nH2 * bh80_h2_prefactor(i_dust) * accomm_factor

        ! 3b. Add CO contribution (accommodation factor = 1)
        Hcoll = Hcoll + nCO * bh80_co_prefactor(i_dust)

        ! 3c. Add H neutral contribution (depends on Td)
        if (nElement(1) > 1d-20 .and. xelem_ions(1,1) > 1d-20) then
            accomm_factor_zero = bh80_h_accomm_zero(i_dust)
            accomm_factor = (1d0 - accomm_factor_zero) * exp(-sqrt((Td+Tgas)/250d0)) + accomm_factor_zero
            Hcoll = Hcoll + nElement(1) * xelem_ions(1,1) * bh80_species_prefactor(i_dust,1) * accomm_factor
        end if

        Hcoll = prefactor * Hcoll

    end subroutine compute_dust_coll_heating_BH80

    subroutine compute_dust_coll_heating(i_dust,ne,nElement,xelem_ions,&
                                        Coulomb_factor,nH2,nCO,Tgas,Td,&
                                        dust_charge,Hcoll)
        implicit none

        integer, intent(in) :: i_dust
        real(dp), intent(in) :: ne
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), dimension(-1:n_elements),intent(in) :: Coulomb_factor
        real(dp), intent(in) :: nH2,nCO
        real(dp), intent(in) :: Tgas
        real(dp), intent(in) :: Td,dust_charge

        real(dp), intent(inout) :: Hcoll

        integer(kind=8) :: c_start, c_end, c_rate
        integer :: j,iel,iphi0,nT,nphi,nions_loc,izion,jj
        real(dp) :: lT,cooling_rate
        real(dp) :: agrain,D
        real(dp) :: dT
        real(dp) :: xion,phi_charge
        real(dp) :: supp_factor,supp_factor_inv
        real(dp):: Hcoll_HM80

        real(dp) :: sum_rate, sum_element
        integer :: idx_T
        
        call system_clock(c_start, c_rate)

        ! Initialize local variables
        dT = Tgas - Td
        sum_rate = 0d0
        Hcoll_HM80 = 0d0
        supp_factor = 0d0
        supp_factor_inv = 0d0

        if (Tgas > 1d3) then
            Hcoll    = 0d0
            lT = log10(Tgas)

            sum_rate = 0d0
            agrain = dustbins_props(i_dust)%asize_cm

            ! 1. Do first the contribution from electron collisions
            nT = dustbins_props(i_dust)%collisional_tab(0)%npts(1)
            nphi = dustbins_props(i_dust)%collisional_tab(0)%npts(2)
            idx_T = -1
            if (dust_coll_charge) then
                phi_charge = dust_charge * dustbins_props(i_dust)%phi_prefact(-1) ! [eV]
                call dustbins_props(i_dust)%collisional_tab(0)%interpolate(lT, phi_charge, cooling_rate, idx_x=idx_T)
                cooling_rate = exp(cooling_rate * ln10)
                sum_rate = sum_rate + Coulomb_factor(-1) * ne * cooling_rate
            else
                ! 1D interpolation: use stored phi=0 index
                call dustbins_props(i_dust)%collisional_tab(0)%interpolate(lT, cooling_rate, idx_x=idx_T)
                cooling_rate = exp(cooling_rate * ln10)
                sum_rate = sum_rate + ne * cooling_rate
            end if
 
            ! 2. Loop over all elements
            species_loop: do iel = 1, n_elements
                
                ! Skip if tables not initialized or element abundance is negligible
                if (nElement(iel) <= 1d-10) cycle
                if (.not. dustbins_props(i_dust)%collisional_tab(iel)%initialised) cycle
                
                nT = dustbins_props(i_dust)%collisional_tab(iel)%npts(1)
                nphi = dustbins_props(i_dust)%collisional_tab(iel)%npts(2)
 
                if (dust_coll_charge) then
                    ! Add contributions for all charge states of this element
                    nions_loc = n_elements
#ifdef RTZ
                    nions_loc = max(1, elements(iel)%n_ions)
#endif
                    sum_element = 0d0
                    do j = 1, nions_loc
                        xion = xelem_ions(iel, j)
                        if (xion <= 1d-5) cycle
                        phi_charge = dust_charge * dustbins_props(i_dust)%phi_prefact(j-1) ! [eV]
                        call dustbins_props(i_dust)%collisional_tab(iel)%interpolate(lT, phi_charge, cooling_rate, idx_x=idx_T)
                        if (cooling_rate < -30d0) cycle
                        sum_element = sum_element + Coulomb_factor(j-1) * xion * exp(cooling_rate * ln10)
                    end do
                    sum_rate = sum_rate + nElement(iel) * sum_element
                else
                    ! No charge dependence: just add contribution from the total abundance of this element
                    call dustbins_props(i_dust)%collisional_tab(iel)%interpolate(lT, cooling_rate, idx_x=idx_T)
                    if (cooling_rate >= -30d0) then
                        sum_rate = sum_rate + nElement(iel) * exp(cooling_rate * ln10)
                    end if
                end if
            
            end do species_loop

            Hcoll = sum_rate * dT

            ! 3. Always blend with BH80 for smooth transition at low T
            if (Tgas < 1d5) then
                supp_factor = 1d0 - 1d0/(1d0+exp(-1d1*(log10(Tgas)-4d0)))
                supp_factor_inv = 1d0 / (1d0 + exp(-1d1*(log10(Tgas) - 4d0)))
                call compute_dust_coll_heating_BH80(i_dust,nElement,xelem_ions,nH2,nCO,Tgas,Td,Hcoll_HM80)
                Hcoll = supp_factor_inv * Hcoll + supp_factor * Hcoll_HM80
            end if

            call system_clock(c_end)
            time_coll_heating = time_coll_heating + dble(c_end - c_start) / dble(c_rate)

            ! write(*,*) 'Time spent in compute_dust_coll_heating (s):', time_coll_heating
            ! write(*,*) 'sum_rate:',sum_rate
            ! write(*,*) 'dT:',dT
            ! write(*,*) 'Hcoll:',Hcoll
            ! write(*,*) 'Hcoll_HM80:',Hcoll_HM80
            ! write(*,*) 'supp_factor:',supp_factor
            ! write(*,*) 'supp_factor_inv:',supp_factor_inv
            ! write(*,*) 'Tgas:',Tgas
            ! write(*,*) 'Td:',Td
            ! write(*,*) 'ne:',ne
            ! write(*,*) 'dust_charge:',dust_charge
            ! write(*,*) 'i_dust:',i_dust
            ! call clean_stop
        else
            call compute_dust_coll_heating_BH80(i_dust,nElement,xelem_ions,nH2,nCO,Tgas,Td,Hcoll)
            Hcoll_HM80 = Hcoll
        end if

    end subroutine compute_dust_coll_heating

end module dust_cooling