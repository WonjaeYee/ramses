module dust_rates
    use amr_parameters, only: dp
    use constants
    use dust_commons
    use dust_utils

    implicit none

contains

    subroutine LeBourlot2012_accretion_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax)
        ! Compute the accretion rate of dust grains in the unrestricted case, following Le Bourlot et al. (2012).
        ! dust_info --> DustChemistryInfo type with all the necessary information to compute the accretion rate.
        ! y_gas     --> 2D array with the gas phase abundances [g cm-3]
        ! y_dust    --> 1D array with the dust phase abundances [g cm-3]
        ! dydt_gas  <--> 2D array with the time derivative of the gas phase abundances [g cm-3 s-1]
        ! dydt_dust <--> 1D array with the time derivative of the dust phase abundances [g cm-3 s-1]
        ! kmax      --> Maximum allowed rate for the process (optional output)
        
        implicit none
        ! ---- Input/Output variables ----
        class(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: dydt_gas(:,:), dydt_dust(:)
        real(dp), intent(inout), optional :: kmax

        ! ---- Local variables ----
        integer :: jj, ii, ii1, ii2, kk, e_index
        integer :: n_el
        real(dp) :: pseudo_rate, rate, prefactor, Tk_loc, limit_rate
        real(dp),dimension(1:ndust) :: correction_factors
        real(dp) :: diff_rate, diff_rho, diff_nH, diff_T

        Tk_loc = dust_info%local_Tk
        prefactor = sqrt(Tk_loc) / (1d0 + 1d-4*Tk_loc**1.5d0)

        speciesloop: do jj = 1, ndchemtype
            ! 1. Loop over the dust chemical species.
            ii1 = istart_chemtype(jj) 
            ii2 = ii1 + dustbins_per_chemtype(jj) - 1

            associate(bin => dustbins_props(ii1))
                n_el = bin%nelements
                if (n_el == 1) then
                    ! 2. A single-element chemistry type has a trivial limiter.
                    e_index = bin%el_index(1)
                    limit_rate = y_gas(e_index,1) / (bin%el_mfractions(1) * sqrt(bin%el_atomic_masses_g(1)))
                else
                    ! 3. Find the limiting element in a single pass, without a temporary array.
                    e_index = bin%el_index(1)
                    limit_rate = y_gas(e_index,1) / (bin%el_mfractions(1) * sqrt(bin%el_atomic_masses_g(1)))
                    do kk = 2, n_el
                        e_index = bin%el_index(kk)
                        pseudo_rate = y_gas(e_index,1) / (bin%el_mfractions(kk) * sqrt(bin%el_atomic_masses_g(kk)))
                        if (pseudo_rate < limit_rate) then
                            limit_rate = pseudo_rate
                        end if
                    end do
                end if

                ! 4. Apply the same limiting rate to every dust bin in the chemical type.
                do ii = ii1, ii2
                    rate = limit_rate * dustbins_props(ii)%k0_acc * prefactor ! [s-1]
                    ! 5. Get the maximum rate computed here, if requested.
                    if (present(kmax)) then
                        kmax = max(kmax, abs(rate))
                    end if
                    rate = rate * y_dust(ii+dust_info%npah) ! [g cm-3 s-1]
                    dydt_dust(ii+dust_info%npah) = dydt_dust(ii+dust_info%npah) + rate  ! [g cm-3 s-1]
                    do kk = 1, n_el
                        e_index = bin%el_index(kk)
                        dydt_gas(e_index,1) = dydt_gas(e_index,1) - rate * bin%el_mfractions(kk) ! [g cm-3 s-1]
                    end do
                end do
            end associate
        end do speciesloop
    end subroutine LeBourlot2012_accretion_rate

    subroutine Aoyama2017_coagulation_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax)

        implicit none

        ! ---- Input/Output variables ----
        class(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: dydt_gas(:,:), dydt_dust(:)
        real(dp), intent(inout), optional :: kmax

        ! ---- Local variables ----
        integer :: jj, ii, ii1, ii2, index
        real(dp) :: rate1, rate2

        speciesloop: do jj = 1, ndchemtype
            ! 1. Loop over the dust chemical species.
            ii1 = istart_chemtype(jj) 
            ii2 = ii1 + dustbins_per_chemtype(jj) - 1

            do ii = ii1, ii2-1
                index = ii+dust_info%npah
                if ((dust_info%local_Tk .gt.1d4) .or. &
                    (dust_info%local_nH .lt. dustbins_props(ii)%nH_coa) .or. &
                    (dust_info%local_Jeans .gt. 4d0*dust_info%local_dx)) then
                    cycle
                end if
                rate1 = dustbins_props(ii)%k0_coa(1) * y_dust(index) / dust_info%local_nH ! [s-1]
                if (present(kmax)) then
                    kmax = max(kmax, rate1)
                end if
                rate2 = rate1 * y_dust(index) ! [g cm-3 s-1]
                dydt_dust(index) = dydt_dust(index) - rate2 ! [g cm-3 s-1]
                dydt_dust(index+1) = dydt_dust(index+1) + rate2 ! [g cm-3 s-1]
            end do
        end do speciesloop
    end subroutine Aoyama2017_coagulation_rate

    subroutine turbulent_coagulation_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax)
        use dust_dynamics, only: grain_relative_velocity

        implicit none

        ! ---- Input/Output variables ----
        class(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: dydt_gas(:,:), dydt_dust(:)
        real(dp), intent(inout), optional :: kmax

        ! ---- Local variables ----
        integer :: jj, ii, ii1, ii2, index
        real(dp) :: rate1, rate2
        real(dp) :: temp_sigma, temp_L
        real(dp) :: v_rel, R, v_coag, enhan_factor

        if (dust_eq_test) then
            temp_sigma = 5.67d5 * (dust_info%local_nH/1d2)**(-0.25d0)
            temp_L = 10d0 * pc2cm * (dust_info%local_nH/1d2)**(-1d0/3d0)
        else
            temp_sigma = dust_info%local_sigma
            temp_L = dust_info%local_dx
        end if

        speciesloop: do jj = 1, ndchemtype
            ! 1. Loop over the dust chemical species.
            ii1 = istart_chemtype(jj) 
            ii2 = ii1 + dustbins_per_chemtype(jj) - 1

            do ii = ii1, ii2-1
                index = ii+dust_info%npah
                
                ! 2. Compute the relative velocity between grains of the same size
                v_rel = grain_relative_velocity(dust_velocity_model,dust_info%local_Tk,&
                                                dust_info%local_rho,dust_info%local_nH,&
                                                temp_sigma,dust_info%local_mu,temp_L,&
                                                dustbins_props(ii)%asize_cm,&
                                                dustbins_props(ii)%asize_cm,&
                                                dustbins_props(ii)%sgrain,&
                                                dustbins_props(ii)%sgrain,&
                                                dustbins_props(ii)%mgrain,&
                                                dustbins_props(ii)%mgrain)
                
                ! 3. Compute the enhancement due to ice mantles
                v_coag = dustbins_props(ii)%vthresh_coag(1) ! [cm/s]
                if (poppe_ice_enhancement) then
                    enhan_factor = (1d0-sigmoid_function(4d0,log10(dustbins_props(jj)%nhmax_acc),log10(dust_info%local_nH))) + &
                                    sigmoid_function(4d0,log10(dustbins_props(jj)%nhmax_acc),log10(dust_info%local_nH)) * 4d0
                    v_coag = enhan_factor * v_coag
                end if
                if (v_rel .gt. v_coag) then
                    cycle
                end if

                ! 4. Collision rate calculation
                ! The factor of sqrt(8/(3*pi)) is taken from Guillet et al. (2020) and
                ! Marchand et al. (2021) and considers that grain velocities along the x-,
                ! y-, and z-axes are Gaussian distributed
                rate1 = dustbins_props(ii)%k0_coa(1) * v_rel * y_dust(index)
                if (present(kmax)) then
                    kmax = max(kmax, rate1)
                end if
                rate2 = rate1 * y_dust(index) ! [g cm-3 s-1]
                dydt_dust(index) = dydt_dust(index) - rate2 ! [g cm-3 s-1]
                dydt_dust(index+1) = dydt_dust(index+1) + rate2 ! [g cm-3 s-1]
            end do
        end do speciesloop

    end subroutine turbulent_coagulation_rate

    subroutine turbulent_all_coagulation_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax)
        use dust_dynamics, only: grain_relative_velocity

        implicit none

        ! ---- Input/Output variables ----
        class(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: dydt_gas(:,:), dydt_dust(:)
        real(dp), intent(inout), optional :: kmax

        ! ---- Local variables ----
        integer :: jj, ii, ii1, ii2, index, kk, kk1, kk2, kk_loc, idest
        real(dp) :: rate1, rate2
        real(dp) :: temp_sigma, temp_L
        real(dp) :: v_rel, v_coag, enhan_factor
        real(dp) :: mi, mk, msum, loss_ii, loss_kk

        if (dust_eq_test) then
            temp_sigma = 5.67d5 * (dust_info%local_nH/1d2)**(-0.25d0)
            temp_L = 10d0 * pc2cm * (dust_info%local_nH/1d2)**(-1d0/3d0)
        else
            temp_sigma = dust_info%local_sigma
            temp_L = dust_info%local_dx
        end if

        speciesloop: do jj = 1, ndchemtype
            ! 1. Loop over the dust chemical species.
            ii1 = istart_chemtype(jj) 
            ii2 = ii1 + dustbins_per_chemtype(jj) - 1

            do ii = ii1, ii2-1
                index = ii+dust_info%npah
                kk1 = ii
                kk2 = ii2
                do kk = kk1, kk2
                    kk_loc = kk - ii1 + 1
                    ! 2. Compute the relative velocity between grains of the same size
                    v_rel = grain_relative_velocity(dust_velocity_model,dust_info%local_Tk,&
                                                    dust_info%local_rho,dust_info%local_nH,&
                                                    temp_sigma,dust_info%local_mu,temp_L,&
                                                    dustbins_props(ii)%asize_cm,&
                                                    dustbins_props(kk)%asize_cm,&
                                                    dustbins_props(ii)%sgrain,&
                                                    dustbins_props(kk)%sgrain,&
                                                    dustbins_props(ii)%mgrain,&
                                                    dustbins_props(kk)%mgrain)
                    
                    ! 3. Compute the enhancement due to ice mantles
                    v_coag = dustbins_props(ii)%vthresh_coag(kk_loc) ! [cm/s]
                    if (poppe_ice_enhancement) then
                        enhan_factor = (1d0-sigmoid_function(4d0,log10(dustbins_props(jj)%nhmax_acc),log10(dust_info%local_nH))) + &
                                        sigmoid_function(4d0,log10(dustbins_props(jj)%nhmax_acc),log10(dust_info%local_nH)) * 4d0
                        v_coag = enhan_factor * v_coag
                    end if
                    if (v_rel .gt. v_coag) then
                        cycle
                    end if

                    ! 4. Collision rate calculation
                    ! The factor of sqrt(8/(3*pi)) is taken from Guillet et al. (2020) and
                    ! Marchand et al. (2021) and considers that grain velocities along the x-,
                    ! y-, and z-axes are Gaussian distributed
                    rate1 = dustbins_props(ii)%k0_coa(kk_loc) * v_rel * y_dust(index)

                    if (present(kmax)) then
                        kmax = max(kmax, rate1)
                    end if
                    rate2 = rate1 * y_dust(kk+dust_info%npah) ! [g cm-3 s-1]
                    if (ii == kk) then
                        ! Identical-grain collisions must be halved to avoid double counting.
                        rate2 = 0.5d0 * rate2
                    end if

                    mi = dustbins_props(ii)%mgrain
                    mk = dustbins_props(kk)%mgrain
                    msum = mi + mk
                    idest = dustbins_props(ii)%idend_coag(kk_loc)
                    loss_ii = rate2 * (mi / msum)
                    loss_kk = rate2 * (mk / msum)

                    ! Three-term update for each pair: two sinks (ii, kk) and one gain (idest).
                    dydt_dust(index) = dydt_dust(index) - loss_ii
                    dydt_dust(kk+dust_info%npah) = dydt_dust(kk+dust_info%npah) - loss_kk
                    dydt_dust(idest+dust_info%npah) = dydt_dust(idest+dust_info%npah) + (loss_ii + loss_kk)
                end do
            end do
        end do speciesloop

    end subroutine turbulent_all_coagulation_rate

end module dust_rates