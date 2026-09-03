module dust_interface
    use amr_commons, only:dp,ndim
    use constants
    use dust_commons
    implicit none

    private

    public :: compute_dust_rad_rates,compute_dust_precool,&
            compute_dust_coolrates,compute_local_anisotropy_factor,&
            compute_dust_update

    ! Cache to avoid automatic array allocations in compute_dust_update
    real(dp), allocatable, save, target :: y_gas_cache(:,:), y_gas_out_cache(:,:)
    real(dp), allocatable, save, target :: y_dust_cache(:), y_dust_out_cache(:)

contains
    subroutine compute_local_anisotropy_factor(dinfo,Fp,Np)
        ! Computes the local radiation anisotropy factor and solid angle 
        ! subtended by radiation sources based on the local radiation field
        ! dinfo --> the DustChemistryInfo instance to update with the computed values
        ! Fp --> the radiation flux vector (#/s/cm^2) for each radiation group
        ! Np --> the radiation energy density (#/cm^3) for each radiation group
        implicit none

        ! ---- Inputs -----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), dimension(:,:), intent(in) :: Fp
        real(dp), dimension(:), intent(in) :: Np

        ! ---- Local variables ----
        integer :: i
        real(dp) :: rad_ani

        if (all(Np.le.dinfo%smallNp)) return

        if (fixed_rad_ani .eq. -1d0) then
            do i = 1, size(Np)
                if (Np(i) .le. dinfo%smallNp) cycle
                rad_ani = sqrt(sum(Fp(:,i)**2d0)) / (Np(i) * dinfo%local_c)
                dinfo%local_rad_ani(i) = rad_ani
                dinfo%local_solid_angle(i) = 2d0 * pi * (1d0 + (1d0 - rad_ani)**2d0) 
            end do
        else
            dinfo%local_rad_ani = fixed_rad_ani
            dinfo%local_solid_angle = 2d0 * pi * (1d0 + (1d0 - fixed_rad_ani)**2d0)
        end if
    end subroutine compute_local_anisotropy_factor

    subroutine compute_dust_rad_rates(dinfo, G0, Tk, ne,&
                                        &dustAbs, dustSc, dustRp,&
                                        &pahAbs, pahSc, pahRp)
        ! Computes the dust and PAH radiative rates (absorption, 
        ! scattering, and radiation pressure) for each radiation group, 
        ! given the local dust properties and radiation field.
        ! dinfo --> the DustChemistryInfo instance containing the dust properties
        ! G0 --> local radiation field strength in units of the Habing field
        ! Tk --> local gas temperature (K)
        ! ne --> local electron density (cm^-3)
        ! dustAbs --> output array for dust absorption rates [1/s]
        ! dustSc --> output array for dust scattering rates [1/s]
        ! dustRp --> output array for dust radiation pressure rates [1/s]
        ! pahAbs --> output array for PAH absorption rates [1/s]
        ! pahSc --> output array for PAH scattering rates [1/s]
        ! pahRp --> output array for PAH radiation pressure rates [1/s]
        use pah_photoelectric_heating, only: interpolate_pah_charge_equilibrium
        use dust_radiation, only: rad_dust_rate, rad_pah_rate

        implicit none

        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0, Tk, ne
        ! ---- Outputs ----
        real(dp), dimension(:), intent(inout) :: dustAbs, dustSc, dustRp
        real(dp), dimension(:), intent(inout) :: pahAbs, pahSc, pahRp
        ! ---- Local variables ----
        integer :: ii

        ! 1. Compute and add the different dust radiative rates
        dustAbs = 0d0; dustSc = 0d0; dustRp = 0d0
        if (dinfo%ndust > 0) then
            do ii = 1, dinfo%nGroups
                dustAbs(ii) = rad_dust_rate(dinfo%csa_dust(:,ii),dinfo%rho_dust(:)) ! [1/s]
                dustSc (ii) = rad_dust_rate(dinfo%css_dust(:,ii),dinfo%rho_dust(:)) ! [1/s]
                dustRp (ii) = rad_dust_rate(dinfo%csr_dust(:,ii),dinfo%rho_dust(:)) ! [1/s]
            end do
            dinfo%dustAbs = dustAbs
        end if

        ! 2. Compute the different PAH radiative rates
        pahAbs = 0d0; pahSc = 0d0; pahRp = 0d0
        if (dinfo%npah > 0) then
            ! Before we move onto computing the PAH contribution
            ! to absorption and scattering, we need to compute the
            ! PAH charge distribution, since the PAH cross-sections
            ! depend on the charge state
            ! TODO: currently we use the interpolation tables since
            ! it's sufficiently accurate and much faster, but we 
            ! could change to use the full calculation if we wanted
            do ii = 1, dinfo%npah
                call interpolate_pah_charge_equilibrium(ii,G0,ne,Tk,&
                   &dinfo%fcharge_pah(:,ii))
            end do
            do ii = 1, dinfo%nGroups
                pahAbs(ii) = pahAbs(ii)+ rad_pah_rate(dinfo%csa_pah(:,ii),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
                pahSc (ii) = pahSc(ii) + rad_pah_rate(dinfo%css_pah(:,ii),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
                pahRp (ii) = pahRp(ii) + rad_pah_rate(dinfo%csr_pah(:,ii),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
            end do
            dinfo%pahAbs = pahAbs
        end if

        ! 3. Add the contribution from dust and PAHs for the 
        !    returned absorption, scattering, and radiation pressure rates
        dustAbs = dustAbs + pahAbs
        dustSc = dustSc + pahSc
        dustRp = dustRp + pahRp
    end subroutine compute_dust_rad_rates

    subroutine compute_dust_precool(dinfo, G0_total, Tk, ne,&
                                    &nElement, xelem_ions, nH2, nCO, &
                                    & Np)
        ! Pre-computes the dust and PAH cooling and heating rates for the given local conditions, 
        ! which can then be used in subsequent calls to compute_dust_coolrates to save computational time.
        ! dinfo --> the DustChemistryInfo instance to update with the computed rates
        ! G0_total --> local radiation field strength in units of the Habing field
        ! Tk --> local gas temperature (K)
        ! ne --> local electron density (cm^-3)
        ! nElement --> array of number densities for each element (cm^-3)
        ! xelem_ions --> array of ionization fractions for each element and ionization state
        ! nH2 --> local molecular hydrogen density (cm^-3)
        ! nCO --> local carbon monoxide density (cm^-3)
        ! Np --> the radiation energy density (#/cm^3) for each radiation group
        use dust_charging, only: compute_mean_dust_charge, compute_dust_charge_sigma,&
                                compute_dust_charge_dist, compute_Coulomb_focusing
        use dust_photoelectric_heating, only: interpolate_dust_peh_rate,&
                                            compute_dust_peh_rate
        use pah_photoelectric_heating, only: interpolate_pah_charge_equilibrium,&
                                            compute_pah_peh_equilibrium,&
                                            interpolate_pah_peh_equilibrium,&
                                            compute_pah_charge_equilibrium
        use dust_radiation, only: update_T_dust
        use dust_surface_chemistry, only: grain_h2_formation_rate

        implicit none

        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: nH2, nCO
        real(dp), dimension(1:dinfo%nGroups), intent(in), optional :: Np

        ! ---- Local variables ----
        integer :: ii,j,idx_g,idx_T
        integer :: i_neutral, i_charged   ! csa_pah row indices: neutral=2*ii-1, charged=2*ii
        real(dp) :: Zel, nHI, prevD
        integer :: n_charge

        if (dinfo%ndust > 0) then
            ! 1. Compute the equilibrium dust charge
            idx_g = -1
            idx_T = -1
            do ii = 1, dinfo%ndust
                call compute_mean_dust_charge(ii,G0_total,Tk,ne,dinfo%Z_dust(ii),idx_g,idx_T)
                call compute_dust_charge_sigma(ii,G0_total,Tk,ne,dinfo%Z_sigma(ii),idx_g,idx_T)
            end do

            ! 2. If needed, precompute the Coulomb factors
            dinfo%Coulomb_factor = 1d0
            if (Coulomb_precompute) then
                do ii = 1, dinfo%ndust
                    prevD = 1d0
                    do j = -1, dinfo%nion_charges
                        Zel = dble(j)
                        if (prevD <= 1d-5 .and. dinfo%Z_dust(ii)*Zel > 0d0) then
                            dinfo%Coulomb_factor(ii,j:dinfo%nion_charges) = 1d-10
                            exit
                        end if
                        call compute_Coulomb_focusing(ii,Tk,dinfo%Z_dust(ii),dinfo%Z_sigma(ii),&
                                                       &Zel,dinfo%Coulomb_factor(ii,j))
                        prevD = dinfo%Coulomb_factor(ii,j)
                    end do
                end do
            end if

            ! 3. Compute the equilibrium dust photoelectric heating and recombination cooling rates
            if (dust_pe_heating .and. present(Np)) then
                do ii = 1, dinfo%ndust
                    if (dust_pe_heating_isrf .or. all(Np.le.dinfo%smallNp)) then
                        call interpolate_dust_peh_rate(ii,G0_total,ne,Tk,&
                                                        &dinfo%Pinj_dust(ii),&
                                                        &dinfo%Prec_dust(ii))
                    else
                        call compute_dust_charge_sigma(ii,G0_total,Tk,ne,dinfo%Z_sigma(ii))
                        ! Consider the contribution from the UV background
                        call interpolate_dust_peh_rate(ii,dinfo%G0_background,ne,Tk,&
                                                        &dinfo%Pinj_dust(ii),dinfo%Prec_dust(ii))
                        ! Now consider the contribution from the local radiation field
                        call compute_dust_peh_rate(ii,dinfo%csa_dust(ii,:),&
                                                    dinfo%l_a(ii,:),dinfo%nGroups,dinfo%local_c,&
                                                    dinfo%local_solid_angle,Np(:),dinfo%group_eV(:),&
                                                    dinfo%Z_dust(ii),dinfo%Z_sigma(ii),Tk,ne,&
                                                    dinfo%Pinj_dust(ii),dinfo%Prec_dust(ii))
                    end if
                end do
            elseif (dust_pe_heating) then
                do ii = 1, dinfo%ndust
                    call interpolate_dust_peh_rate(ii,G0_total,ne,Tk,&
                                                    &dinfo%Pinj_dust(ii),&
                                                    &dinfo%Prec_dust(ii))
                end do
            end if

            ! 4. Compute the internal energy of the dust grain considering all heating and cooling processes
            if (present(Np)) then
                call update_T_dust(dinfo%G0_background,dinfo%Pcoll_dust(:),dinfo%Prec_dust(:),dinfo%Pinj_dust(:),dinfo%Prad_dust(:),&
                                    dinfo%T_dust(:),ne,nElement(:),xelem_ions(:,:),dinfo%Coulomb_factor(:,:),nH2,nCO,Tk,dinfo%Z_dust(:),&
                                    Np(:)*dinfo%group_eV(:),dinfo%csa_dust(:,:))
            else
                call update_T_dust(G0_total,dinfo%Pcoll_dust(:),dinfo%Prec_dust(:),dinfo%Pinj_dust(:),dinfo%Prad_dust(:),&
                                    dinfo%T_dust(:),ne,nElement(:),xelem_ions(:,:),dinfo%Coulomb_factor(:,:),nH2,nCO,Tk,dinfo%Z_dust(:))
            end if

            ! 5. Now convert all the rates from erg/s per grain to erg/s/cm3
            do ii = 1, dinfo%ndust
                dinfo%n_dust(ii) = dinfo%rho_dust(ii) / dustbins_props(ii)%mgrain
                dinfo%Pcoll_dust(ii) = dinfo%Pcoll_dust(ii) * dinfo%n_dust(ii)
                ! print*,'ii,dinfo%Pcoll_dust(ii),dinfo%n_dust(ii): ',ii,dinfo%Pcoll_dust(ii),dinfo%n_dust(ii)
                dinfo%Prec_dust(ii) = dinfo%Prec_dust(ii) * dinfo%n_dust(ii)
                dinfo%Pinj_dust(ii) = dinfo%Pinj_dust(ii) * dinfo%n_dust(ii)
                dinfo%Prad_dust(ii) = dinfo%Prad_dust(ii) * dinfo%n_dust(ii)
            end do

            ! 6. Compute the H2 formation rate on dust grains
            if (H2ondust) then
                nHI = nElement(1) * xelem_ions(1,1)
                dinfo%H2_formation_rate = grain_h2_formation_rate(nHI,nElement(1),Tk,dinfo%rho_dust(:),dinfo%T_dust(:))
            end if
        end if

        if (dinfo%npah > 0) then
            ! 7. Compute the PAH PEH model
            if (pah_pe_heating .and. present(Np)) then
                if (pah_pe_heating_isrf) then
                    ! We have rt, but we want to just use the ISRF-averaged model
                        do ii = 1, dinfo%npah
                            call interpolate_pah_peh_equilibrium(ii,G0_total,&
                                                                ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                                dinfo%Pabs_pah(ii,1),dinfo%Pinj_pah(ii),&
                                                                dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                        end do
                else
                    ! We want the PAH PEH also have rt, so we use the full model
                    do ii = 1, dinfo%npah
                        i_neutral = 2*ii - 1  ! row for neutral PAH cross-sections in csa_pah
                        i_charged = 2*ii      ! row for charged PAH cross-sections in csa_pah
                        ! Consider the contribution from the UV background first
                        call interpolate_pah_peh_equilibrium(ii,dinfo%G0_background,&
                                                            ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                            dinfo%Pabs_pah(ii,1),dinfo%Pinj_pah(ii),&
                                                            dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                        ! And now the full model for the local radiation field
                        call compute_pah_peh_equilibrium(ii,dinfo%csa_pah(i_neutral,:),&
                                                            dinfo%csa_pah(i_neutral,:),&
                                                            dinfo%csa_pah(i_charged,:),&
                                                            dinfo%csa_pah(i_charged,:),&
                                                            dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                            Np(:),dinfo%group_eV(:),&
                                                            dinfo%local_c,Tk,ne,dinfo%fcharge_pah(:,ii),&
                                                            dinfo%Pabs_pah(ii,:),dinfo%Pinj_pah(ii),&
                                                            dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                    end do
                end if
            else if (pah_pe_heating) then
                ! We want PAH PEH but don't have rt
                do ii = 1, dinfo%npah
                    call interpolate_pah_peh_equilibrium(ii,dinfo%G0_background,&
                                                        ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                        dinfo%Pabs_pah(ii,1),dinfo%Pinj_pah(ii),&
                                                        dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                end do
            else if (present(Np)) then
                ! We don't want PAH PEH but have rt, so we still want to compute the PAH charge distribution
                do ii = 1, dinfo%npah
                    i_neutral = 2*ii - 1
                    i_charged = 2*ii
                    call compute_pah_charge_equilibrium(ii,dinfo%csa_pah(i_neutral,:),&
                                                        dinfo%csa_pah(i_neutral,:),&
                                                        dinfo%csa_pah(i_charged,:),&
                                                        dinfo%csa_pah(i_charged,:),&
                                                        dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                        Np(:),dinfo%group_eV(:),dinfo%local_c,&
                                                        Tk,ne,dinfo%fcharge_pah(:,ii))
                end do
            else
                ! We don't want PAH PEH but don't have rt, so we just interpolate
                ! the ISRF-averaged PAH charge tables
                do ii = 1, dinfo%npah
                    call interpolate_pah_charge_equilibrium(ii,G0_total,ne,Tk,&
                       &dinfo%fcharge_pah(:,ii))
                end do
            end if

            ! 8. Now convert from erg/s to erg/s/cm3
            do ii = 1, dinfo%npah
                dinfo%n_pah(ii) = dinfo%rho_pah(ii) / pahbins_props(ii)%mpah
                dinfo%Pabs_pah(ii,:) = dinfo%Pabs_pah(ii,:) * dinfo%n_pah(ii)
                dinfo%Pinj_pah(ii) = dinfo%Pinj_pah(ii) * dinfo%n_pah(ii)
                dinfo%Prad_pah(ii) = dinfo%Prad_pah(ii) * dinfo%n_pah(ii)
                dinfo%Prec_pah(ii) = dinfo%Prec_pah(ii) * dinfo%n_pah(ii)
            end do
        end if
    end subroutine compute_dust_precool

    subroutine compute_dust_coolrates(dinfo, G0_total, Tk, ne,&
                                    &nElement, xelem_ions, nH2, nCO, &
                                    &total_rec_power, total_inj_power, total_col_power,&
                                    &H2_formation_rate,Np)
        ! Computes the dust and PAH cooling and heating rates for the given local conditions,
        ! using the pre-computed rates from compute_dust_precool if available to save computational time.
        ! dinfo --> the DustChemistryInfo instance to update with the computed rates
        ! G0_total --> local radiation field strength in units of the Habing field
        ! Tk --> local gas temperature (K)
        ! ne --> local electron density (cm^-3)
        ! nElement --> array of number densities for each element (cm^-3)
        ! xelem_ions --> array of ionization fractions for each element and ionization state
        ! nH2 --> local molecular hydrogen density (cm^-3)
        ! nCO --> local carbon monoxide density (cm^-3)
        ! total_rec_power --> output total recombination cooling power from dust and PAHs (erg/s/cm^3)
        ! total_inj_power --> output total photoelectric heating power from dust and PAHs (erg/s/cm^3)
        ! total_col_power --> output total collisional cooling power from dust (erg/s/cm^3)
        ! H2_formation_rate --> output H2 formation rate on dust grains (cm^3/s)
        ! Np --> the radiation energy density (#/cm^3) for each radiation group
        implicit none
        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: nH2, nCO
        real(dp), intent(out) :: total_rec_power, total_inj_power, total_col_power
        real(dp), intent(out) :: H2_formation_rate
        real(dp), dimension(1:dinfo%nGroups), intent(in), optional :: Np


        ! If precomputed rates are not already in dinfo (i.e. compute_dust_precool was not
        ! called before us this step), delegate to it now.  This eliminates ~148 lines
        ! of duplicated physics that previously lived in the non-precomp branch here.
        if (.not. dinfo%use_precomp) then
            call compute_dust_precool(dinfo, G0_total, Tk, ne, nElement, xelem_ions, nH2, nCO, Np=Np)
        end if

        ! Read the precomputed rates from dinfo (populated either just above or by an
        ! external compute_dust_precool call earlier in the RTZ step).
        total_rec_power = 0d0
        total_inj_power = 0d0
        total_col_power = 0d0
        H2_formation_rate = -1d0
        if (dinfo%ndust > 0) then
            total_rec_power = sum(dinfo%Prec_dust)
            total_inj_power = sum(dinfo%Pinj_dust)
            total_col_power = sum(dinfo%Pcoll_dust)
            if (H2ondust) then
                H2_formation_rate = dinfo%H2_formation_rate
            end if
        end if
        if (dinfo%npah > 0) then
            total_rec_power = total_rec_power + sum(dinfo%Prec_pah)
            total_inj_power = total_inj_power + sum(dinfo%Pinj_pah)
        end if
        dinfo%use_precomp = .false.
    end subroutine compute_dust_coolrates

    subroutine compute_dust_update(dinfo,nElement,xelem_ions,dt,Np,step_ok, dUU)

        use ode_driver_mod, only: integrate_dust_ode
        use ode_interface_mod, only: dust_solver_step
        use dust_rhs_mod, only: dust_rhs
        use dust_rates, only: compute_rate_caches
        use dust_radiative_torques, only: total_radiative_torque,IR_damping_factor
#ifdef RTZ
        use rtz_module, only:elements
#endif

        implicit none

        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in), optional :: Np(:)
        logical, intent(out) :: step_ok

        ! --- Local variables ----
        integer :: ii, jj, ndust_total
        real(dp) :: sum_check
        real(dp), pointer :: y_gas(:,:), y_gas_out(:,:)
        real(dp), pointer :: y_dust(:), y_dust_out(:)
        real(dp) :: mass_g(n_elements)   ! element atomic masses — populated once via #ifdef, used throughout
        real(dp) :: h_init_val           ! substep hint for integrate_dust_ode
        ! Cache the last accepted substep's recommended next-step size across calls.
        ! On re-entry with the same (or similar) conditions this skips the firstcall
        ! 1/kmax reduction, cutting substep count at equilibrium for explicit solvers.
        real(dp), save :: last_h_ode = 0.0_dp

        real(dp), intent(out) :: dUU

        dUU = 0.0_dp ! should I keep previous values?

        ! Populate mass_g once so the pack/unpack loops below are #ifdef-free.
#ifdef RTZ
        do ii = 1, n_elements
            mass_g(ii) = elements(ii)%atomic_mass_g
        end do
#else
        mass_g(:) = el_atomic_masses_g(1:n_elements)
#endif

        ! If no dust or PAH process is active, just return
        if (ndust_processes.eq.0 .and. npah_processes.eq.0) then
            step_ok = .true.
            return
        end if

        ! Precompute/cache rate factors for this cell-update step
        call compute_rate_caches(dinfo, nElement)

        ! 1. Compute the local RAT-D quantities if we run with dust_ratd
        if (dust_ratd) then
            if (present(Np)) then
                do ii = 1, dinfo%ndust
                    dinfo%rat_torque(ii) = total_radiative_torque(dinfo%local_rad_ani,Np,&
                                                            dinfo%group_eV(:)*eV2erg,dinfo%nGroups,&
                                                            dinfo%local_c,dinfo%csrat_dust(:,ii))
                end do
            end if
            dinfo%IR_damp_factor = IR_damping_factor(dinfo%local_G0*1.13d0,nElement(1),dinfo%local_Tk,dinfo%T_dust)
            do ii = 1, dinfo%ndust
                dinfo%rat_torque(ii) = dinfo%rat_torque(ii) + dustbins_props(ii)%RAT_torque_0 * dinfo%G0_background
            end do
        end if

        ! 2. Now setup the arrays for gas and dust quantities to send to the ODE solver
        ndust_total = 0
        if (dinfo%ndust > 0) ndust_total = ndust_total + dinfo%ndust
        if (dinfo%npah > 0) ndust_total = ndust_total + dinfo%npah

        if (carry_gas_ions) then
            call ensure_update_cache(n_elements, n_elements + 1, ndust_total)
        else
            call ensure_update_cache(n_elements, 1, ndust_total)
        end if
        y_gas => y_gas_cache; y_gas_out => y_gas_out_cache
        y_dust => y_dust_cache; y_dust_out => y_dust_out_cache

        if (carry_gas_ions) then
            do ii = 1, n_elements
                y_gas(ii,1) = nElement(ii) * mass_g(ii)
                y_gas(ii,2:n_elements+1) = nElement(ii) * xelem_ions(ii,:) * mass_g(ii)
            end do
        else
            do ii = 1, n_elements
                y_gas(ii,1) = nElement(ii) * mass_g(ii)
            end do
        end if
        ! Pack y_dust: PAH bins first [1:npah], then grain bins [npah+1:npah+ndust].
        ! This convention is mirrored in the unpack block below and in integrate_dust_ode.
        if (dinfo%ndust > 0 .and. dinfo%npah > 0) then
            y_dust(1:dinfo%npah) = dinfo%rho_pah(1:dinfo%npah)
            y_dust(dinfo%npah+1:dinfo%npah+dinfo%ndust) = dinfo%rho_dust(1:dinfo%ndust)
        else if (dinfo%ndust > 0) then
            y_dust(1:dinfo%ndust) = dinfo%rho_dust(1:dinfo%ndust)
        else if (dinfo%npah > 0) then
            y_dust(1:dinfo%npah) = dinfo%rho_pah(1:dinfo%npah)
        end if

        ! 3. Now we are ready to call the ODE solver to integrate the dust evolution.
        ! Guard against an uninitialised procedure pointer (set by init_dust_solver based
        ! on dust_solver_type); a null dereference here would be a silent SIGSEGV.
        if (.not. associated(dust_solver_step)) then
            print *, "FATAL ERROR: dust_solver_step is not associated. Check dust_solver_type in namelist."
            call clean_stop
        end if
        step_ok = .true.
        h_init_val = dt
        if (last_h_ode > 0.0_dp) h_init_val = last_h_ode
        call integrate_dust_ode(dinfo,dt,y_gas,y_dust,dust_rhs,dust_solver_step,&
                                y_gas_out,y_dust_out,h_init_val,0d0,dt,&
                                debug_flag=dust_debug,step_ok=step_ok,h_last=last_h_ode)

        ! 4. Update the dinfo with the new values after the ODE step
        if (carry_gas_ions) then
            do ii = 1, n_elements
                nElement(ii) = y_gas_out(ii,1) / mass_g(ii)
                if (nElement(ii) > 1d-30) then
                    xelem_ions(ii,:) = y_gas_out(ii,2:n_elements+1) / nElement(ii) / mass_g(ii)
                    ! Make sure that xelem_ions add up to 1 for each element
                    sum_check = sum(xelem_ions(ii,:))
                    if (sum_check > 1d-30) then
                        xelem_ions(ii,:) = xelem_ions(ii,:) / sum_check
                    end if
                end if
            end do
        else
            do ii = 1, n_elements
                nElement(ii) = y_gas_out(ii,1) / mass_g(ii)
            end do
        end if

        ! Unpack y_dust_out: PAH bins first [1:npah], grain bins [npah+1:npah+ndust].
        if (dinfo%ndust > 0 .and. dinfo%npah > 0) then
            dinfo%rho_pah(1:dinfo%npah) = y_dust_out(1:dinfo%npah)
            dinfo%rho_dust(1:dinfo%ndust) = y_dust_out(dinfo%npah+1:dinfo%npah+dinfo%ndust)
        else if (dinfo%ndust > 0) then
            dinfo%rho_dust(1:dinfo%ndust) = y_dust_out(1:dinfo%ndust)
        else if (dinfo%npah > 0) then
            dinfo%rho_pah(1:dinfo%npah) = y_dust_out(1:dinfo%npah)
        end if

#ifdef RTZ
        ! Compute dUU: the maximum fractional change in any gas or dust species over this step.
        ! dUU is returned *unscaled* (pure fractional change, dimensionless).
        ! The caller (rtz_cool_step in rtz_cooling_module.f90) multiplies by one_over_x_FRAC
        ! before comparing against the RTZ convergence criterion.

        ! check for gas
        do ii = 1, n_elements
            if (elements(ii)%atomic_number > 0) then
                dUU = max(dUU, abs(y_gas_out(ii,1)-y_gas(ii,1))/(y_gas(ii,1)+y_min))

                ! if not, ion part is zero
                if (carry_gas_ions) then
                    do jj=1,elements(ii)%n_ions
                        dUU = max(dUU, abs(y_gas_out(ii,1+jj)-y_gas(ii,1+jj))/(y_gas(ii,1+jj)+y_min))
                    end do
                end if

            end if
        end do

        ! check for dust
        if (dinfo%ndust > 0 .or. dinfo%npah > 0) then
            do ii = 1, dinfo%npah+dinfo%ndust
                dUU = max(dUU, abs(y_dust_out(ii)-y_dust(ii))/(y_dust(ii)+y_min))
            end do
        end if
#endif

    end subroutine compute_dust_update

    subroutine ensure_update_cache(n1, n2, ndust_total)
        integer, intent(in) :: n1, n2, ndust_total
        logical :: need_realloc

        need_realloc = .false.
        if (.not. allocated(y_gas_cache)) then
            need_realloc = .true.
        else if (size(y_gas_cache,1) /= n1 .or. size(y_gas_cache,2) /= n2) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(y_gas_cache)) then
                deallocate(y_gas_cache)
                deallocate(y_gas_out_cache)
            end if
            allocate(y_gas_cache(n1, n2))
            allocate(y_gas_out_cache(n1, n2))
        end if

        need_realloc = .false.
        if (.not. allocated(y_dust_cache)) then
            need_realloc = .true.
        else if (size(y_dust_cache) /= ndust_total) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(y_dust_cache)) then
                deallocate(y_dust_cache)
                deallocate(y_dust_out_cache)
            end if
            allocate(y_dust_cache(ndust_total))
            allocate(y_dust_out_cache(ndust_total))
        end if
    end subroutine ensure_update_cache

end module dust_interface