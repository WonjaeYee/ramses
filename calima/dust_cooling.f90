module dust_cooling

    use dust_commons
    use dust_charging
#ifdef RTZ
    use rtz_module, only:elements
#endif
    contains

    subroutine compute_dust_coll_heating_BH80(ne,nElement,xelem_ions,nH2,nCO,Tgas,Td,Hcoll_each)
        use cooling_module, only: kb,mH
        use constants, only: pi
        implicit none

        real(dp), intent(in) :: ne,nH2,nCO
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), intent(in) :: Tgas
        real(dp), dimension(1:ndust), intent(in) :: Td

        real(dp), dimension(1:ndust), intent(inout) :: Hcoll_each

        integer :: i,j,iel,nions_loc
        real(dp) :: accomm_factor,prefactor,accomm_factor_zero,agrain,mproj,xion

        Hcoll_each(:)    = 0d0

        ! Loop over gas species
        hcol_loop: do i = 1, size(dustbins_props)

            ! 1. Compute common prefactor 
            prefactor = 2d0 * kb * (Tgas - Td(i)) * sqrt(Tgas)
            agrain = dustbins_props(i)%asize_cm

            ! 2. Compute contribution from H2 molecules
            ! Accommodation factor based on Burke & Hollenbach (1980)
            ! (No enhancement factor due to grain charge)
            accomm_factor_zero = 4d0 * mH * dustbins_props(i)%mgrain / (2d0*mH + dustbins_props(i)%mgrain)**2d0
            accomm_factor = (1d0 - accomm_factor_zero) * exp(-sqrt((Td(i)+Tgas)/250d0)) + accomm_factor_zero
            Hcoll_each(i) = Hcoll_each(i) + nH2 * sqrt(8d0*kb/(pi*2d0*mH)) * pi * (agrain)**2d0 * accomm_factor

            ! 2b. CO molecules use accommodation coefficient unity
            Hcoll_each(i) = Hcoll_each(i) + nCO * sqrt(8d0*kb/(pi*(mC_amu+mO_amu)*amu2g)) * pi * (agrain)**2d0

            ! 3. Compute contribution from all elements and ionisation states.
            ! For H neutral use BH80 accommodation; for He and heavier ions/neutrals use unity.
            do iel = 1, n_elements
#ifdef RTZ
                mproj = elements(iel)%atomic_mass * amu2g
                nions_loc = max(1, elements(iel)%n_ions)
#else
                mproj = el_atomic_masses_amu(iel) * amu2g
                nions_loc = 1
#endif
                if (nElement(iel) <= 1d-20) cycle

                do j = 1, nions_loc
                    xion = xelem_ions(iel,j)
                    if (xion <= 1d-20) cycle

                    if (iel == 1 .and. j == 1) then
                        accomm_factor_zero = 2d0 * mH * dustbins_props(i)%mgrain / (mH + dustbins_props(i)%mgrain)**2d0
                        accomm_factor = (1d0 - accomm_factor_zero) * exp(-sqrt((Td(i)+Tgas)/250d0)) + accomm_factor_zero
                    else
                        accomm_factor = 1d0
                    end if

                    Hcoll_each(i) = Hcoll_each(i) + nElement(iel) * xion * sqrt(8d0*kb/(pi*mproj)) * pi * (agrain)**2d0 * accomm_factor
                end do
            end do

            Hcoll_each(i) = prefactor * Hcoll_each(i)

        end do hcol_loop

    end subroutine compute_dust_coll_heating_BH80

    subroutine compute_dust_coll_heating(ne,nElement,xelem_ions,nH2,nCO,Tgas,Td,dust_charge,Hcoll_each)
        use cooling_module, only: kb,mH
        use constants, only: pi
        implicit none

        real(dp), intent(in) :: ne
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), intent(in) :: nH2,nCO
        real(dp), intent(in) :: Tgas
        real(dp), dimension(1:ndust), intent(in) :: Td,dust_charge

        real(dp), dimension(1:ndust), intent(inout) :: Hcoll_each

        integer :: i,j,iel,iphi0,nT,nphi,nions_loc,izion
        real(dp) :: lT,cooling_rate
        real(dp) :: agrain,D
        real(dp), dimension(1:1) :: temp_dist
        real(dp), dimension(1:1) :: temp_charge
        real(dp), dimension(1:ndust) :: dT
        real(dp) :: xion,phi_charge
        real(dp) :: supp_factor,supp_factor_inv
        real(dp), dimension(1:ndust) :: Hcoll_each_HM80
        
        Hcoll_each(:)    = 0d0
        lT = log10(Tgas)
        dT = Tgas - Td

        ! Loop over dust bins
        hcol_loop: do i = 1, size(dustbins_props)

            temp_charge(:) = dust_charge(i)
            temp_dist(:)   = 1d0
            agrain = dustbins_props(i)%asize_cm

            ! Loop over all projectile species (electrons and elements)
            species_loop: do iel = 0, n_elements
                
                ! Skip if tables not initialized
                if (.not. dustbins_props(i)%collisional_tab(iel)%initialised) cycle
                
                nT = dustbins_props(i)%collisional_tab(iel)%npts(1)
                nphi = dustbins_props(i)%collisional_tab(iel)%npts(2)
                
                ! Determine interpolation method and compute cooling rate
                if (dust_coll_charge) then
                    phi_charge = dust_charge(i)
                    call interpolate2D(dustbins_props(i)%collisional_tab(iel)%tab1d(1:nT,1), &
                                       dustbins_props(i)%collisional_tab(iel)%tab1d(1:nphi,2), &
                                       dustbins_props(i)%collisional_tab(iel)%tab2d(1:nT,1:nphi,1), &
                                       nT, nphi, lT, phi_charge, cooling_rate)
                else
                    ! 1D interpolation: use stored phi=0 index
                    iphi0 = dustbins_props(i)%collisional_tab(iel)%ipos_zero(2)
                    call interpolate1D(dustbins_props(i)%collisional_tab(iel)%tab1d(1:nT,1), &
                                       dustbins_props(i)%collisional_tab(iel)%tab2d(1:nT,iphi0,1), &
                                       nT, lT, cooling_rate)
                end if

                cooling_rate = 10d0**cooling_rate

                if (iel == 0) then
                    Hcoll_each(i) = Hcoll_each(i) + ne * cooling_rate * dT(i)
                    cycle species_loop
                end if
                
                ! Add contributions for all charge states of this element
                nions_loc = n_elements
#ifdef RTZ
                nions_loc = max(1, elements(iel)%n_ions)
#endif
                do j = 1, nions_loc
                    xion = xelem_ions(iel, j)
                    if (xion <= 1d-20) cycle

                    if (j > 1) then
                        ! Ionized species: include Coulomb enhancement
                        if (Coulomb_precompute) then
                            izion = max(-1, min(10, j-1))
                            D = dustbins_props(i)%Coulomb_focus_ion(izion)
                        else
                            call compute_Coulomb_focusing(Tgas,agrain,temp_dist,temp_charge,dble(j-1),D)
                        end if
                        Hcoll_each(i) = Hcoll_each(i) + nElement(iel) * xion * cooling_rate * D * dT(i)
                    else
                        ! Neutral species: no Coulomb enhancement
                        Hcoll_each(i) = Hcoll_each(i) + nElement(iel) * xion * cooling_rate * dT(i)
                    end if
                end do
            
            end do species_loop

        end do hcol_loop

        ! Now compute the low-temperature soft-cube collisional heating from
        ! Hollenbach and McKee (1980) for Tgas < 1e4 K
        if (Tgas .lt. 1d5 .and. (dust_coll_lowT)) then
            supp_factor = 1d0 - 1d0/(1d0+exp(-1d1*(log10(Tgas)-4d0)))
            supp_factor_inv = 1d0 / (1d0 + exp(-1d1*(log10(Tgas) - 4d0)))
            call compute_dust_coll_heating_BH80(ne,nElement,xelem_ions,nH2,nCO,Tgas,Td,Hcoll_each_HM80)
            Hcoll_each(:) = supp_factor_inv * Hcoll_each(:) + supp_factor * Hcoll_each_HM80(:)
        end if
    end subroutine compute_dust_coll_heating

end module dust_cooling