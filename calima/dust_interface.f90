module dust_interface
    use amr_commons, only:dp,ndim
    use constants
    use dust_commons
    implicit none

    private

    public :: compute_dust_rad_rates,compute_dust_precool,&
            compute_dust_coolrates,compute_local_anisotropy_factor

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

        if (fixed_rad_ani .eq. -1d0) then
            do i = 1, size(Np)
                if (Np(i) .le. 0d0) cycle
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
        if (dinfo%ndust > 0 .and. dust) then
            do ii = 1, dinfo%nGroups
                dustAbs(ii) = rad_dust_rate(dinfo%csa_dust(ii,:),dinfo%rho_dust(:)) ! [1/s]
                dustSc (ii) = rad_dust_rate(dinfo%css_dust(ii,:),dinfo%rho_dust(:)) ! [1/s]
                dustRp (ii) = rad_dust_rate(dinfo%csr_dust(ii,:),dinfo%rho_dust(:)) ! [1/s]
            end do
            dinfo%dustAbs = dustAbs
        end if

        ! 2. Compute the different PAH radiative rates
        pahAbs = 0d0; pahSc = 0d0; pahRp = 0d0
        if (dinfo%npah > 0 .and. dust_pahs) then
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
                pahAbs(ii) = pahAbs(ii)+ rad_pah_rate(dinfo%csa_pah(ii,:),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
                pahSc (ii) = pahSc(ii) + rad_pah_rate(dinfo%css_pah(ii,:),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
                pahRp (ii) = pahRp(ii) + rad_pah_rate(dinfo%csr_pah(ii,:),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
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
        integer :: ii,j
        real(dp) :: Zel, nHI
        real(dp),dimension(:),allocatable :: Zvals
        real(dp),dimension(:),allocatable :: fcharge
        
        if (dinfo%ndust > 0 .and. dust) then
            ! 1. Compute the equilibrium dust charge
            do ii = 1, dinfo%ndust
                call compute_mean_dust_charge(ii,G0_total,Tk,ne,dinfo%Z_dust(ii))
            end do

            ! 2. If needed, precompute the Coulomb factors
            if (Coulomb_precompute) then
                do ii = 1, dinfo%ndust
                    ! Compute the dust charge distribution (approx. Gaussian)
                    call compute_dust_charge_dist(ii,G0_total,Tk,ne,Zvals,fcharge)
                    do j = -1, dinfo%nion_charges
                        Zel = dble(j)
                        call compute_Coulomb_focusing(Tk,dustbins_props(ii)%asize_cm,&
                                                        fcharge,Zvals,&
                                                        Zel,dinfo%Coulomb_factor(ii,j))
                    end do
                end do
            end if

            ! 3. Compute the equilibrium dust photoelectric heating and recombination cooling rates
            if (dust_pe_heating .and. present(Np)) then
                do ii = 1, dinfo%ndust
                    if (dust_pe_heating_isrf .or. all(Np.eq.0d0)) then
                        call interpolate_dust_peh_rate(ii,dinfo%rho_dust(ii),G0_total,ne,Tk,&
                                                        &dinfo%Pinj_dust(ii),dinfo%Prec_dust(ii))
                    else
                        call compute_dust_charge_sigma(ii,G0_total,Tk,ne,dinfo%Z_sigma(ii))
                        call interpolate_dust_peh_rate(ii,dinfo%rho_dust(ii),dinfo%G0_background,ne,Tk,&
                                                        &dinfo%Pinj_dust(ii),dinfo%Prec_dust(ii))
                        call compute_dust_peh_rate(ii,dinfo%rho_dust(ii),dinfo%csa_dust(:,ii),&
                                                    dinfo%l_a(:,ii),dinfo%nGroups,dinfo%local_c,&
                                                    dinfo%local_solid_angle,Np(:),dinfo%group_eV(:),&
                                                    dinfo%Z_dust(ii),dinfo%Z_sigma(ii),Tk,ne,&
                                                    dinfo%Pinj_dust(ii),dinfo%Prec_dust(ii))
                    end if
                end do
            elseif (dust_pe_heating) then
                do ii = 1, dinfo%ndust
                    call interpolate_dust_peh_rate(ii,dinfo%rho_dust(ii),G0_total,ne,Tk,&
                                                    &dinfo%Pinj_dust(ii),dinfo%Prec_dust(ii))
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

            ! 5. Compute the H2 formation rate on dust grains
            if (H2ondust) then
                nHI = nElement(1) * xelem_ions(1,1)
                dinfo%H2_formation_rate = grain_h2_formation_rate(nHI,Tk,dinfo%rho_dust(:),dinfo%T_dust(:))
            end if
        end if

        if (dinfo%npah > 0 .and. dust_pahs) then
            ! 6. Compute the PAH PEH model
            if (pah_pe_heating .and. present(Np)) then
                if (pah_pe_heating_isrf) then
                    ! We have rt, but we want to just use the ISRF-averaged model
                        do ii = 1, dinfo%npah
                            call interpolate_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),G0_total,&
                                                                ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                                dinfo%Pabs_pah(1,ii),dinfo%Pinj_pah(ii),&
                                                                dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                        end do
                else
                    ! We want the PAH PEH also have rt, so we use the full model
                    do ii = 1, dinfo%npah
                        ! Consider the contribution from the UV background first
                        call interpolate_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),dinfo%G0_background,&
                                                            ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                            dinfo%Pabs_pah(1,ii),dinfo%Pinj_pah(ii),&
                                                            dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                        ! And now the full model for the local radiation field
                        call compute_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                            dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                            dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                            dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                            dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                            Np(:),dinfo%group_eV(:),&
                                                            dinfo%local_c,Tk,ne,dinfo%fcharge_pah(:,ii),&
                                                            dinfo%Pabs_pah(:,ii),dinfo%Pinj_pah(ii),&
                                                            dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                    end do
                end if
            else if (pah_pe_heating) then
                ! We want PAH PEH but don't have rt
                do ii = 1, dinfo%npah
                    call interpolate_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),dinfo%G0_background,&
                                                        ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                        dinfo%Pabs_pah(1,ii),dinfo%Pinj_pah(ii),&
                                                        dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                end do
            else if (present(Np)) then
                ! We don't want PAH PEH but have rt, so we still want to compute the PAH charge distribution
                do ii = 1, dinfo%npah
                    call compute_pah_charge_equilibrium(ii,dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                        dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                        dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                        dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
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
        end if
    end subroutine compute_dust_precool

    subroutine compute_dust_coolrates(dinfo, G0_total, Tk, ne,&
                                    &nElement, xelem_ions, nH2, nCO, &
                                    &total_rec_power, total_inj_power, total_col_power,&
                                    &H2_formation_rate,Np)
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
        real(dp), intent(out) :: total_rec_power, total_inj_power, total_col_power
        real(dp), intent(out) :: H2_formation_rate
        real(dp), dimension(1:dinfo%nGroups), intent(in), optional :: Np

        ! ---- Local variables ----
        integer :: ii,j
        real(dp) :: Zel, nHI
        real(dp),dimension(:),allocatable :: Zvals
        real(dp),dimension(:),allocatable :: fcharge
        real(dp), dimension(1:dinfo%ndust) :: Z_dust, T_dust, Z_sigma
        real(dp), dimension(1:dinfo%ncharge_pah_max,1:dinfo%npah) :: fcharge_pah
        real(dp), dimension(1:dinfo%ndust,-1:dinfo%nion_charges) :: Coulomb_factor
        real(dp), dimension(1:dinfo%ndust):: Pinj_dust, Prec_dust, Pcoll_dust, Prad_dust
        real(dp), dimension(1:dinfo%npah) :: Pinj_pah, Prec_pah, Prad_pah
        real(dp), dimension(1:dinfo%nGroups,1:dinfo%npah) :: Pabs_pah


        if (dinfo%use_precomp) then
            ! In this case we use the rates computed during the pre-cool step
            total_rec_power = 0d0
            total_inj_power = 0d0
            total_col_power = 0d0
            H2_formation_rate = -1d0
            if (dinfo%ndust > 0) then
                total_rec_power = sum(dinfo%Prec_dust)
                total_inj_power = sum(dinfo%Pinj_dust)
                total_col_power = sum(dinfo%Pcoll_dust)
            end if
            if (dinfo%npah  > 0) then
                total_rec_power = total_rec_power + sum(dinfo%Prec_pah)
                total_inj_power = total_inj_power + sum(dinfo%Pinj_pah)
            end if
            if (H2ondust) then
                H2_formation_rate = dinfo%H2_formation_rate
            end if
            dinfo%use_precomp = .false. ! Reset the flag for next time
            return
        end if

        total_rec_power = 0d0
        total_inj_power = 0d0
        total_col_power = 0d0
        H2_formation_rate = -1d0

        if (dinfo%ndust > 0 .and. dust) then
            ! 1. Compute the equilibrium dust charge
            do ii = 1, dinfo%ndust
                call compute_mean_dust_charge(ii,G0_total,Tk,ne,dinfo%Z_dust(ii))
            end do

            ! 2. If needed, precompute the Coulomb factors
            if (Coulomb_precompute) then
                do ii = 1, dinfo%ndust
                    ! Compute the dust charge distribution (approx. Gaussian)
                    call compute_dust_charge_dist(ii,G0_total,Tk,ne,Zvals,fcharge)
                    do j = -1, dinfo%nion_charges
                        Zel = dble(j)
                        call compute_Coulomb_focusing(Tk,dustbins_props(ii)%asize_cm,&
                                                        fcharge,Zvals,Zel,&
                                                        Coulomb_factor(ii,j))
                    end do
                end do
            end if

            ! 3. Compute the equilibrium dust photoelectric heating and recombination cooling rates
            if (dust_pe_heating .and. present(Np)) then
                do ii = 1, dinfo%ndust
                    if (dust_pe_heating_isrf .or. all(Np.eq.0d0)) then
                        call interpolate_dust_peh_rate(ii,dinfo%rho_dust(ii),G0_total,ne,Tk,&
                                                        &Pinj_dust(ii),Prec_dust(ii))
                    else
                        call compute_dust_charge_sigma(ii,G0_total,Tk,ne,Z_sigma(ii))
                        call interpolate_dust_peh_rate(ii,dinfo%rho_dust(ii),dinfo%G0_background,ne,Tk,&
                                                        &Pinj_dust(ii),Prec_dust(ii))
                        call compute_dust_peh_rate(ii,dinfo%rho_dust(ii),dinfo%csa_dust(:,ii),&
                                                    dinfo%l_a(:,ii),dinfo%nGroups,dinfo%local_c,&
                                                    dinfo%local_solid_angle(:),Np(:),&
                                                    dinfo%group_eV(:),Z_dust(ii),Z_sigma(ii),Tk,ne,&
                                                    Pinj_dust(ii),Prec_dust(ii))
                    end if
                end do
            elseif (dust_pe_heating) then
                do ii = 1, dinfo%ndust
                    call interpolate_dust_peh_rate(ii,dinfo%rho_dust(ii),G0_total,ne,Tk,&
                                                    &Pinj_dust(ii),Prec_dust(ii))
                end do
            end if

            ! 4. Compute the internal energy of the dust grain considering all heating and cooling processes
            if (present(Np)) then
                call update_T_dust(dinfo%G0_background,Pcoll_dust(:),Prec_dust(:),Pinj_dust(:),Prad_dust(:),&
                                    T_dust(:),ne,nElement(:),xelem_ions(:,:),Coulomb_factor(:,:),nH2,nCO,Tk,Z_dust(:),&
                                    Np(:)*dinfo%group_eV(:),dinfo%csa_dust(:,:))
            else
                call update_T_dust(dinfo%G0_background,Pcoll_dust(:),Prec_dust(:),Pinj_dust(:),Prad_dust(:),&
                                    T_dust(:),ne,nElement(:),xelem_ions(:,:),Coulomb_factor(:,:),nH2,nCO,Tk,Z_dust(:))
            end if

            total_rec_power = total_rec_power + sum(Prec_dust)
            total_inj_power = total_inj_power + sum(Pinj_dust)
            total_col_power = total_col_power + sum(Pcoll_dust)

            ! 5. Compute the H2 formation rate on dust grains
            if (H2ondust) then
                nHI = nElement(1) * xelem_ions(1,1)
                H2_formation_rate = grain_h2_formation_rate(nHI,Tk,dinfo%rho_dust(:),T_dust(:))
            endif
        end if

        if (dinfo%npah > 0 .and. dust_pahs) then
            ! 5. Compute the PAH PEH model
            if (pah_pe_heating .and. present(Np)) then
                if (pah_pe_heating_isrf) then
                    ! We have rt, but we want to just use the ISRF-averaged model
                        do ii = 1, dinfo%npah
                            call interpolate_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),G0_total,&
                                                                ne,Tk,fcharge_pah(:,ii),&
                                                                Pabs_pah(1,ii),Pinj_pah(ii),&
                                                                Prad_pah(ii),Prec_pah(ii))
                        end do
                else
                    ! We want the PAH PEH also have rt, so we use the full model
                    do ii = 1, dinfo%npah
                        ! Consider the contribution from the UV background first
                        call interpolate_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),dinfo%G0_background,&
                                                            ne,Tk,fcharge_pah(:,ii),&
                                                            Pabs_pah(1,ii),Pinj_pah(ii),&
                                                            Prad_pah(ii),Prec_pah(ii))
                        ! And now the full model for the local radiation field
                        call compute_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                            dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                            dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                            dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                            dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                            Np(:),dinfo%group_eV(:),&
                                                            dinfo%local_c,Tk,ne,fcharge_pah(:,ii),&
                                                            Pabs_pah(:,ii),Pinj_pah(ii),&
                                                            Prad_pah(ii),Prec_pah(ii))
                    end do
                end if
            else if (pah_pe_heating) then
                ! We want PAH PEH but don't have rt
                do ii = 1, dinfo%npah
                    call interpolate_pah_peh_equilibrium(ii,dinfo%rho_pah(ii),dinfo%G0_background,&
                                                        ne,Tk,fcharge_pah(:,ii),&
                                                        Pabs_pah(1,ii),Pinj_pah(ii),&
                                                        Prad_pah(ii),Prec_pah(ii))
                end do
            else if (present(Np)) then
                ! We don't want PAH PEH but have rt, so we still want to compute the PAH charge distribution
                do ii = 1, dinfo%npah
                    call compute_pah_charge_equilibrium(ii,dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                        dinfo%csa_pah(:,1+(ii-1)*dinfo%npah),&
                                                        dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                        dinfo%csa_pah(:,2+(ii-1)*dinfo%npah),&
                                                        dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                        Np(:),dinfo%group_eV(:),dinfo%local_c,&
                                                        Tk,ne,fcharge_pah(:,ii))
                end do
            else
                ! We don't want PAH PEH but don't have rt, so we just interpolate
                ! the ISRF-averaged PAH charge tables
                do ii = 1, dinfo%npah
                    call interpolate_pah_charge_equilibrium(ii,dinfo%G0_background,ne,Tk,&
                       &fcharge_pah(:,ii))
                end do
            end if
            total_rec_power = total_rec_power + sum(Prec_pah)
            total_inj_power = total_inj_power + sum(Pinj_pah)
        end if
    end subroutine compute_dust_coolrates
end module dust_interface