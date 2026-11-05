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
    public :: init_dust_coll_heating_BH80_cache
    public :: CollHeatPre, coll_heating_prepare, coll_heating_eval, coll_heating_dTg, coll_heating_dZ
    public :: bh80_unit_flux

    ! The collisional heating split at the dust temperature: everything but the T_d dependence,
    ! for one bin and one gas state (coll_heating_prepare), then H_coll at any T_d (coll_heating_eval)
    type CollHeatPre
        logical :: hot = .false., blend = .false., hneutral = .false.
        real(dp) :: Tgas = 0d0, sum_rate = 0d0, supp = 0d0, supp_inv = 0d0   ! T > 1e3 K: tables, blend with BH80
        real(dp) :: sqrtTg = 0d0, s_const = 0d0, nh2_pf = 0d0, a0_h2 = 0d0   ! BH80: species with unit accommodation
        real(dp) :: nhi_pf = 0d0, a0_h = 0d0                                 ! (CO included), H2 and neutral H
        ! with deriv: d sum_rate / d T_gas, d sum_rate / d Z (dust_coll_charge), d supp_inv / d T_gas
        real(dp) :: dsum_dT = 0d0, dsum_dZ = 0d0, dsupp_inv = 0d0
    end type CollHeatPre

    ! BH80 constants (init_dust_coll_heating_BH80_cache). A species s hits a grain of radius a at
    ! n_s pi a^2 sqrt(8 kB T / (pi m_s)): the mean thermal speed over sqrt(T), sqrt(8 kB / (pi m_s))
    ! [cm s^-1 K^-1/2], of each element, H2 and CO; and of each bin its cross section pi a^2 [cm^2] and
    ! the accommodation coefficients at zero of H2 and neutral H
    real(dp), save :: bh80_vth(1:n_elements) = 0d0, bh80_vth_h2 = 0d0, bh80_vth_co = 0d0
    real(dp), save :: bh80_area(1:ndust) = 0d0
    real(dp), save :: bh80_h2_accomm_zero(1:ndust) = 0d0, bh80_h_accomm_zero(1:ndust) = 0d0

    contains

    subroutine init_dust_coll_heating_BH80_cache
        implicit none

        integer :: iel, i
        real(dp) :: mproj, msurf

        bh80_vth_h2 = sqrt(8d0*kB/(pi*2d0*mH))
        bh80_vth_co = sqrt(8d0*kB/(pi*(mC_amu+mO_amu)*amu2g))
        do iel = 1, n_elements
#ifdef RTZ
            mproj = elements(iel)%atomic_mass * amu2g
#else
            mproj = el_atomic_masses_amu(iel) * amu2g
#endif
            bh80_vth(iel) = 0d0
            if (mproj > 0d0) bh80_vth(iel) = sqrt(8d0*kB/(pi*mproj))
        end do

        do i = 1, ndust
            bh80_area(i) = pi * dustbins_props(i)%asize_cm**2
            ! BH83 accommodation coefficient: alpha_0 = 4*m_proj*m_s/(m_proj+m_s)^2, the energy a
            ! projectile leaves with a surface atom of mass m_s (hard cube): the mean atomic mass of the
            ! grain material, 1/sum(f_k/m_k) over its mass fractions f_k, not the mass of the whole
            ! grain (which gave alpha_0 ~ 4 m_H/m_grain ~ 1e-6)
            msurf = 1d0 / sum(dustbins_props(i)%el_mfractions(1:dustbins_props(i)%nelements) &
                              / dustbins_props(i)%el_atomic_masses_g(1:dustbins_props(i)%nelements))
            bh80_h2_accomm_zero(i) = 8d0 * mH * msurf / (2d0*mH + msurf)**2d0
            bh80_h_accomm_zero(i) = 4d0 * mH * msurf / (mH + msurf)**2d0
        end do
    end subroutine init_dust_coll_heating_BH80_cache

    subroutine compute_dust_coll_heating_BH80(i_dust,nElement,xelem_ions,nH2,nCO,Tgas,Td,Hcoll)
        ! Gas-grain collisional energy exchange after Burke & Hollenbach (the bh80_* routines), with
        ! their soft-cube treatment of the energy accommodation. To do, not yet planned: replace it by a
        ! more recent model of the gas-grain energy accommodation.
        implicit none

        integer, intent(in) :: i_dust
        real(dp), intent(in) :: nH2,nCO
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), intent(in) :: Tgas
        real(dp), intent(in) :: Td

        real(dp), intent(inout) :: Hcoll

        type(CollHeatPre) :: pre

        call bh80_prepare(i_dust,nElement,xelem_ions,nH2,nCO,bh80_unit_flux(nElement,xelem_ions),Tgas,pre)
        Hcoll = bh80_eval(pre,Tgas,Td)

    end subroutine compute_dust_coll_heating_BH80

    real(dp) function bh80_unit_flux(nElement,xelem_ions)
        ! BH80: sum of n_s sqrt(8 kB / (pi m_s)) [cm^-2 s^-1 K^-1/2] over the gas species whose
        ! accommodation factor is unity, all ions and atoms but neutral H (H2 and CO are not in
        ! nElement): their flux onto a grain over its cross section and sqrt(T). The same for every
        ! bin, so the caller of a loop over the bins computes it once per gas state.
        implicit none

        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions

        integer :: j,iel,nions_loc
        real(dp) :: xion,s

        s = 0d0
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

                s = s + nElement(iel) * xion * bh80_vth(iel)
            end do
        end do
        bh80_unit_flux = s
    end function bh80_unit_flux

    subroutine bh80_prepare(i_dust,nElement,xelem_ions,nH2,nCO,flux_bh80,Tgas,pre)
        ! the T_d-independent part of the BH80 rate of bin i_dust; flux_bh80: bh80_unit_flux of this gas
        implicit none

        integer, intent(in) :: i_dust
        real(dp), intent(in) :: nH2,nCO,flux_bh80
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), intent(in) :: Tgas
        type(CollHeatPre), intent(inout) :: pre

        ! 1. Common prefactor: 2 kB (Tgas - Td) sqrt(Tgas), finished in bh80_eval
        pre%sqrtTg = sqrt(Tgas)

        ! 2. The species whose accommodation factor is unity: the ions and atoms but neutral H, and CO
        pre%s_const = bh80_area(i_dust) * (flux_bh80 + nCO * bh80_vth_co)

        ! 3. H2 and H neutral (accommodation depends on Td)
        pre%nh2_pf = bh80_area(i_dust) * nH2 * bh80_vth_h2
        pre%a0_h2 = bh80_h2_accomm_zero(i_dust)
        pre%hneutral = nElement(1) > 1d-20 .and. xelem_ions(1,1) > 1d-20
        if (pre%hneutral) then
            pre%nhi_pf = bh80_area(i_dust) * nElement(1) * xelem_ions(1,1) * bh80_vth(1)
            pre%a0_h = bh80_h_accomm_zero(i_dust)
        end if

    end subroutine bh80_prepare

    real(dp) function bh80_eval(pre,Tgas,Td)
        ! the BH80 rate at dust temperature Td from bh80_prepare
        implicit none

        type(CollHeatPre), intent(in) :: pre
        real(dp), intent(in) :: Tgas,Td

        real(dp) :: accomm_factor,prefactor,e,Hcoll

        prefactor = 2d0 * kB * (Tgas - Td) * pre%sqrtTg
        e = exp(-sqrt((Td+Tgas)/250d0))
        accomm_factor = (1d0 - pre%a0_h2) * e + pre%a0_h2
        Hcoll = pre%s_const + pre%nh2_pf * accomm_factor
        if (pre%hneutral) then
            accomm_factor = (1d0 - pre%a0_h) * e + pre%a0_h
            Hcoll = Hcoll + pre%nhi_pf * accomm_factor
        end if
        bh80_eval = prefactor * Hcoll

    end function bh80_eval

    subroutine coll_heating_prepare(i_dust,ne,nElement,xelem_ions,nH2,nCO,flux_bh80,Tgas,dust_charge,pre,deriv)
        ! The collisional heating of bin i_dust but its dependence on the dust temperature (the tables of
        ! the electrons and ions, and the BH80 terms that do not depend on Td). flux_bh80: bh80_unit_flux
        ! of this gas. deriv (optional): also what coll_heating_dTg and coll_heating_dZ need (the slopes
        ! of the tables)
        implicit none

        integer, intent(in) :: i_dust
        real(dp), intent(in) :: ne
        real(dp), dimension(1:n_elements), intent(in) :: nElement
        real(dp), dimension(1:n_elements,1:n_elements), intent(in) :: xelem_ions
        real(dp), intent(in) :: nH2,nCO,flux_bh80
        real(dp), intent(in) :: Tgas
        real(dp), intent(in) :: dust_charge
        type(CollHeatPre), intent(out) :: pre
        logical, intent(in), optional :: deriv

        integer :: j,iel,nions_loc
        real(dp) :: lT,cooling_rate,rate_lin
        real(dp) :: xion,phi_charge

        real(dp) :: sum_rate, sum_element
        real(dp) :: dsl, dsz, sx, sy
        real(dp), dimension(0:n_elements) :: lr, dlr   ! the tables at phi = 0 (log10) and their slopes
        integer :: idx_T
        logical :: want_d

        want_d = .false.
        if (present(deriv)) want_d = deriv
        dsl = 0d0
        dsz = 0d0
        pre%Tgas = Tgas
        pre%hot = Tgas > 1d3
        if (pre%hot) then
            lT = log10(Tgas)

            sum_rate = 0d0

            if (dust_coll_charge) then
                ! 1. Do first the contribution from electron collisions
                idx_T = -1
                phi_charge = dust_charge * dustbins_props(i_dust)%phi_prefact(-1) ! [eV]
                call dustbins_props(i_dust)%collisional_tab(0)%interpolate(lT, phi_charge, cooling_rate, idx_x=idx_T)
                cooling_rate = exp(cooling_rate * ln10)
                sum_rate = sum_rate + ne * cooling_rate
                if (want_d) then
                    call dustbins_props(i_dust)%collisional_tab(0)%slope_2d(lT, phi_charge, sx, sy)
                    dsl = dsl + ne * cooling_rate * sx
                    dsz = dsz + ne * cooling_rate * sy * dustbins_props(i_dust)%phi_prefact(-1)
                end if

                ! 2. Loop over all elements, adding the contributions of all charge states of each
                species_loop: do iel = 1, n_elements

                    ! Skip if tables not initialized or element abundance is negligible
                    if (nElement(iel) <= 1d-10) cycle
                    if (.not. dustbins_props(i_dust)%collisional_tab(iel)%initialised) cycle

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
                        rate_lin = exp(cooling_rate * ln10)
                        sum_element = sum_element + xion * rate_lin
                        if (want_d) then
                            call dustbins_props(i_dust)%collisional_tab(iel)%slope_2d(lT, phi_charge, sx, sy)
                            dsl = dsl + nElement(iel) * xion * rate_lin * sx
                            dsz = dsz + nElement(iel) * xion * rate_lin * sy * &
                                  dustbins_props(i_dust)%phi_prefact(j-1)
                        end if
                    end do
                    sum_rate = sum_rate + nElement(iel) * sum_element
                end do species_loop
            else
                ! No charge dependence: the tables at phi = 0, all on one T axis (coll_tabs), at once
                if (want_d) then
                    call dustbins_props(i_dust)%coll_tabs%eval(lT, lr, dlr)
                else
                    call dustbins_props(i_dust)%coll_tabs%eval(lT, lr)
                end if

                ! 1. The electron collisions
                rate_lin = exp(lr(0) * ln10)
                sum_rate = sum_rate + ne * rate_lin
                if (want_d) dsl = dsl + ne * rate_lin * dlr(0)

                ! 2. The total abundance of each element
                do iel = 1, n_elements
                    if (nElement(iel) <= 1d-10) cycle
                    if (.not. dustbins_props(i_dust)%collisional_tab(iel)%initialised) cycle
                    if (lr(iel) < -30d0) cycle
                    rate_lin = exp(lr(iel) * ln10)
                    sum_rate = sum_rate + nElement(iel) * rate_lin
                    if (want_d) dsl = dsl + nElement(iel) * rate_lin * dlr(iel)
                end do
            end if
            pre%sum_rate = sum_rate
            ! d 10**y / d T_gas = 10**y (dy / dlog10 T) / T_gas; d 10**y / d Z = 10**y ln10 (dy / dphi) dphi/dZ
            pre%dsum_dT = dsl / Tgas
            pre%dsum_dZ = ln10 * dsz

            ! 3. Always blend with BH80 for smooth transition at low T
            pre%blend = Tgas < 1d5
            if (pre%blend) then
                pre%supp = 1d0 - 1d0/(1d0+exp(-1d1*(log10(Tgas)-4d0)))
                pre%supp_inv = 1d0 / (1d0 + exp(-1d1*(log10(Tgas) - 4d0)))
                pre%dsupp_inv = pre%supp_inv * (1d0 - pre%supp_inv) * 1d1 / (Tgas * ln10)
                call bh80_prepare(i_dust,nElement,xelem_ions,nH2,nCO,flux_bh80,Tgas,pre)
            end if
        else
            call bh80_prepare(i_dust,nElement,xelem_ions,nH2,nCO,flux_bh80,Tgas,pre)
        end if

    end subroutine coll_heating_prepare

    real(dp) function coll_heating_dTg(pre,Td)
        ! d H_coll / d T_gas at dust temperature Td, fixed Td and charge (pre prepared with deriv)
        implicit none
        type(CollHeatPre), intent(in) :: pre
        real(dp), intent(in) :: Td
        real(dp) :: dH

        if (pre%hot) then
            dH = pre%sum_rate + (pre%Tgas - Td) * pre%dsum_dT
            if (pre%blend) dH = pre%supp_inv * dH + pre%dsupp_inv * pre%sum_rate * (pre%Tgas - Td) &
                                + pre%supp * bh80_dTg(pre,Td) - pre%dsupp_inv * bh80_eval(pre,pre%Tgas,Td)
        else
            dH = bh80_dTg(pre,Td)
        end if
        coll_heating_dTg = dH
    end function coll_heating_dTg

    real(dp) function coll_heating_dZ(pre,Td)
        ! d H_coll / d Z (the grain charge, dust_coll_charge) at dust temperature Td, fixed T_gas and Td
        implicit none
        type(CollHeatPre), intent(in) :: pre
        real(dp), intent(in) :: Td
        real(dp) :: dH

        dH = 0d0
        if (pre%hot) then
            dH = (pre%Tgas - Td) * pre%dsum_dZ
            if (pre%blend) dH = pre%supp_inv * dH
        end if
        coll_heating_dZ = dH
    end function coll_heating_dZ

    real(dp) function bh80_dTg(pre,Td)
        ! d/d T_gas of bh80_eval at fixed Td
        implicit none
        type(CollHeatPre), intent(in) :: pre
        real(dp), intent(in) :: Td
        real(dp) :: Tg, x, e, de, P, dP, S, dS

        Tg = pre%Tgas
        P = 2d0 * kB * (Tg - Td) * pre%sqrtTg
        dP = 2d0 * kB * (pre%sqrtTg + (Tg - Td) / (2d0 * pre%sqrtTg))
        x = (Td + Tg) / 250d0
        e = exp(-sqrt(x))
        de = -e / (500d0 * sqrt(x))
        S = pre%s_const + pre%nh2_pf * ((1d0 - pre%a0_h2) * e + pre%a0_h2)
        dS = pre%nh2_pf * (1d0 - pre%a0_h2) * de
        if (pre%hneutral) then
            S = S + pre%nhi_pf * ((1d0 - pre%a0_h) * e + pre%a0_h)
            dS = dS + pre%nhi_pf * (1d0 - pre%a0_h) * de
        end if
        bh80_dTg = dP * S + P * dS
    end function bh80_dTg

    real(dp) function coll_heating_eval(pre,Td)
        ! H_coll [erg/s] at dust temperature Td from coll_heating_prepare
        implicit none

        type(CollHeatPre), intent(in) :: pre
        real(dp), intent(in) :: Td

        real(dp) :: Hcoll

        if (pre%hot) then
            Hcoll = pre%sum_rate * (pre%Tgas - Td)
            if (pre%blend) Hcoll = pre%supp_inv * Hcoll + pre%supp * bh80_eval(pre,pre%Tgas,Td)
        else
            Hcoll = bh80_eval(pre,pre%Tgas,Td)
        end if
        coll_heating_eval = Hcoll

    end function coll_heating_eval

end module dust_cooling