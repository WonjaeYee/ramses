!=======================================================================
! CALIMA Dust and Gas Radiation Pressure Force Module
! This module computes the radiation pressure force vectors for the gas
! and individual dust/PAH bins directly in code units for the RTZ network.
!=======================================================================
module dust_radpressure_module
    use amr_parameters, only: dp, ndim
    use amr_commons, only: cosmo, aexp
    use constants, only: c_cgs, eV2erg, mCO
    use rt_parameters, only: nGroups, iGroups, group_egy, rt_pressBoost, &
                            rt_isoPress, rt_isIR, iIR, rt_c, group_csn, &
                            isH2_rtz, iIons, isLW, rtz_UV_background_G0, rt_advect, rt_c_cgs
    use rtz_module, only: elements, n_elements
    use hydro_parameters, only: ndust, npah, idust, ipah, imetal, smallr, smallc, gamma, &
                                neul, nener, gamma_rad, nhydro
#ifdef CO
    use hydro_parameters, only: iCO
#endif
    use dust_commons, only: dustbins_props, pahbins_props, group_csr_dust, group_csr_pah, ncharge_pah_max, GD_solar
    use pah_photoelectric_heating, only: interpolate_pah_charge_equilibrium
    use rtz_cooling_module, only: getNe, getMu_RTZ
    use molecules_module, only: comp_Sd, comp_SH2

    implicit none

    public :: compute_gas_dust_radpressure_force, &
              compute_gas_dust_radpressure_acc

contains

    subroutine compute_gas_dust_radpressure_acc(cell_state, cell_rt_state, ilevel, dx, gas_acc, dust_acc, pah_acc)
        implicit none
        real(dp), dimension(:), intent(in) :: cell_state
        real(dp), dimension(:), intent(in) :: cell_rt_state
        integer, intent(in) :: ilevel
        real(dp), intent(in) :: dx
        real(dp), dimension(1:ndim), intent(out) :: gas_acc
        real(dp), dimension(1:ndim, max(1, ndust)), intent(out) :: dust_acc
        real(dp), dimension(1:ndim, max(1, npah)), intent(out) :: pah_acc

        ! Local variables for cell state extraction
        real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
        real(dp) :: scale_Np, scale_Fp
        real(dp) :: dtg_mw
        real(dp) :: Tk, ne, G0, f_shd
        real(dp), dimension(n_elements) :: nElement
        real(dp), dimension(n_elements, 27) :: xion
        real(dp), dimension(max(1, ndust)) :: rho_dust
        real(dp), dimension(max(1, npah)) :: rho_pah
        real(dp), dimension(1:ndim) :: gas_force
        real(dp), dimension(1:ndim, max(1, ndust)) :: dust_force
        real(dp), dimension(1:ndim, max(1, npah)) :: pah_force
        real(dp) :: rho_loc, P_gas, eps_tot, rho_gas, mu
        integer :: id, iion, counter, e_counter, jbin
#ifdef CO
        real(dp) :: nCO
#endif

        ! Get scale factors
        call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
        call rt_units(scale_Np, scale_Fp)

        rho_loc = max(cell_state(1), smallr)

        ! 1. Extract elements and ionization fractions
        counter = 0
        e_counter = 0
        xion = 0d0
        do id=1,n_elements
           if (elements(id)%atomic_number.gt.0) then
              do iion=1,elements(id)%n_ions
                 xion(id,iion) = cell_state(iIons+counter) / rho_loc
                 counter = counter + 1
              end do
              elements(id)%scale_n = scale_d / elements(id)%atomic_mass_g
              nElement(id) = cell_state(imetal+e_counter) * elements(id)%scale_n
              e_counter = e_counter + 1
           else
              nElement(id) = 0d0
           end if
        end do
        if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
           xion(1,3) = cell_state(iIons+counter) / rho_loc
           counter = counter + 1
        end if

        ! 2. Extract dust and PAH densities
        eps_tot = 0.0_dp
        if (ndust > 0) then
           do jbin = 1, ndust
              rho_dust(jbin) = cell_state(idust+jbin-1)
              eps_tot = eps_tot + (rho_dust(jbin) / rho_loc)
           end do
        end if
        if (npah > 0) then
           do jbin = 1, npah
              rho_pah(jbin) = cell_state(ipah+jbin-1)
           end do
        end if
        eps_tot = min(max(eps_tot, 0.0_dp), 1.0_dp - smallr)
        rho_gas = max((1.0_dp - eps_tot) * rho_loc, smallr)

        ! 3. H2 shielding
        f_shd = 1.d0
        if (isH2_rtz) then
           dtg_mw = sum(rho_dust(1:ndust)) / rho_loc * GD_solar
           f_shd = comp_SH2(0.5d0*nElement(1)*xion(1,3),dx) * &
                   comp_Sd(nElement(1)*xion(1,1),0.5d0*nElement(1)*xion(1,3),dx,dtg_mw)
        end if

        ! 4. Electron density
        ne = getNe(xion, nElement)

        ! 5. Temperature
        P_gas = (gamma - 1) * (cell_state(neul) - 0.5_dp * sum(cell_state(2:ndim+1)**2) / cell_state(1))
#if NENER > 0
        do counter = 1, nener
           P_gas = P_gas - (gamma_rad(counter) - 1.0_dp) * cell_state(nhydro+counter)
        end do
#endif
        P_gas = max(P_gas, rho_gas * (smallc**2/gamma))

#ifdef CO
        nCO = cell_state(iCO) * scale_d / mCO
        mu = getMu_RTZ(ne, nElement, xion, nCO)
#else
        mu = getMu_RTZ(ne, nElement, xion)
#endif
        Tk = P_gas / rho_loc * scale_T2 * mu

        ! 6. Habing band radiation field G0
        G0 = rtz_UV_background_G0
        if (rt_advect) then
           do id = 1, nGroups
              if (group_egy(id).gt.5.6d0 .and. group_egy(id).lt.13.6d0) then
                 G0 = G0 + scale_Np * cell_rt_state(iGroups(id)) * rt_c_cgs(ilevel) * group_egy(id) * eV2erg / 1.6d-3
              end if
           end do
        end if

        ! 7. Call compute_gas_dust_radpressure_force
        call compute_gas_dust_radpressure_force(cell_rt_state, ilevel, &
            gas_force, dust_force, pah_force, &
            nElement, xion, rho_dust, rho_pah, &
            Tk, ne, G0, f_shd)

        ! 8. Calculate accelerations
        gas_acc = gas_force / rho_gas

        if (ndust > 0) then
           do jbin = 1, ndust
              if (rho_dust(jbin) > 0d0) then
                 dust_acc(:, jbin) = dust_force(:, jbin) / rho_dust(jbin)
              else
                 dust_acc(:, jbin) = 0d0
              end if
           end do
        end if

        if (npah > 0) then
           do jbin = 1, npah
              if (rho_pah(jbin) > 0d0) then
                 pah_acc(:, jbin) = pah_force(:, jbin) / rho_pah(jbin)
              else
                 pah_acc(:, jbin) = 0d0
              end if
           end do
        end if

    end subroutine compute_gas_dust_radpressure_acc

    subroutine compute_gas_dust_radpressure_force(cell_rt_state, ilevel, &
        gas_force_code, dust_force_code, pah_force_code, &
        nElement, dXion, rho_dust, rho_pah, &
        Tk, ne, G0, f_shd)

        ! Input/Output declarations
        real(dp), dimension(:), intent(in) :: cell_rt_state
        integer, intent(in) :: ilevel
        
        real(dp), dimension(1:ndim), intent(out) :: gas_force_code
        real(dp), dimension(:,:), intent(out) :: dust_force_code
        real(dp), dimension(:,:), intent(out) :: pah_force_code
        
        real(dp), dimension(n_elements), intent(in) :: nElement
        real(dp), dimension(n_elements, 27), intent(in) :: dXion
        real(dp), dimension(:), intent(in) :: rho_dust
        real(dp), dimension(:), intent(in) :: rho_pah
        real(dp), intent(in) :: Tk, ne, G0, f_shd

        ! Local variables
        real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
        real(dp) :: scale_E, rt_c_code
        integer :: igroup, idim, iNp, ii, jj
        real(dp) :: fluxMag_code, Np_code
        real(dp), dimension(1:ndim) :: Fp_code
        real(dp), dimension(nGroups) :: opacity_gas_cgs
        real(dp), dimension(nGroups) :: opacity_gas_code
        real(dp), dimension(1:ndust, nGroups) :: opacity_dust_code
        real(dp), dimension(1:npah, nGroups) :: opacity_pah_code
        real(dp), dimension(nGroups) :: group_egy_erg
        real(dp), dimension(nGroups) :: group_egy_code
        real(dp) :: mom_fact_gas, mom_fact_dust, mom_fact_pah
        real(dp) :: pah_ion_fraction
        real(dp), dimension(max(1, ncharge_pah_max), max(1, npah)) :: fcharge_pah_local

        ! 1. Initialize forces to 0
        gas_force_code = 0d0
        if (ndust > 0) dust_force_code = 0d0
        if (npah > 0) pah_force_code = 0d0

        ! 2. Get scaling factors
        call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
        scale_E = scale_d * (scale_l**5) / (scale_t**2)

        rt_c_code = rt_c(ilevel)

        ! 3. Calculate gas opacities (CGS)
        opacity_gas_cgs = 0d0

        ! Calculate opacity_gas_cgs for gas directly using group_csn (in cm2)
        do ii=1,n_elements
            if (elements(ii)%atomic_number.gt.0) then
                do jj=1,elements(ii)%n_ions-1
                    do igroup=1,nGroups
                        opacity_gas_cgs(igroup) = opacity_gas_cgs(igroup) + nElement(ii) * dXion(ii, jj) * group_csn(igroup, ii, jj)
                    end do
                end do
            end if
        end do

        if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
            do igroup=1,nGroups
                if (isLW(igroup).eq.1) then
                    opacity_gas_cgs(igroup) = opacity_gas_cgs(igroup) + 0.5d0 * nElement(1) * dXion(1, 3) * group_csn(igroup, 1, 3) * f_shd
                else
                    opacity_gas_cgs(igroup) = opacity_gas_cgs(igroup) + 0.5d0 * nElement(1) * dXion(1, 3) * group_csn(igroup, 1, 3)
                end if
            end do
        end if

        ! Convert gas opacity to code units
        do igroup=1,nGroups
            opacity_gas_code(igroup) = opacity_gas_cgs(igroup) * scale_l
        end do

        ! 4. Calculate dust/PAH opacities in code units
        opacity_dust_code = 0d0
        opacity_pah_code = 0d0
        if (ndust > 0) then
            do igroup = 1, nGroups
                do ii = 1, ndust
                    opacity_dust_code(ii, igroup) = group_csr_dust(ii, igroup) * rho_dust(ii) * &
                                                    scale_d * scale_l / dustbins_props(ii)%mgrain
                end do
            end do
        end if

        if (npah > 0) then
            do ii = 1, npah
                ! Calculate PAH charge equilibrium fraction
                ! Use PAH photoelectric heating solver for charge state
                call interpolate_pah_charge_equilibrium(ii, G0, ne, Tk, fcharge_pah_local(:,ii))
                pah_ion_fraction = fcharge_pah_local(2, ii)
                do igroup = 1, nGroups
                    opacity_pah_code(ii, igroup) = (group_csr_pah(ii, igroup) * (1d0 - pah_ion_fraction) + &
                                                    group_csr_pah(npah+ii, igroup) * pah_ion_fraction) * &
                                                   rho_pah(ii) * scale_d * scale_l / &
                                                   pahbins_props(ii)%mpah
                end do
            end do
        end if

        ! Calculate group energies in ergs and convert to code units
        group_egy_erg = group_egy * eV2erg
        do igroup = 1, nGroups
            group_egy_code(igroup) = group_egy_erg(igroup) / scale_E
        end do

        ! 5. Calculate radiation pressure forces directly in code units
        do igroup = 1, nGroups
            iNp = iGroups(igroup)
            Np_code = cell_rt_state(iNp)
            Fp_code = cell_rt_state(iNp+1 : iNp+ndim)
            fluxMag_code = sqrt(sum(Fp_code**2))

            ! --- Gas Force ---
            mom_fact_gas = opacity_gas_code(igroup) * group_egy_code(igroup)

            if (rt_isoPress .and. .not. (rt_isIR .and. igroup == iIR)) then
                if (fluxMag_code > 0d0) then
                    mom_fact_gas = mom_fact_gas * Np_code / fluxMag_code
                else
                    mom_fact_gas = 0d0
                end if
            else
                mom_fact_gas = mom_fact_gas / rt_c_code
            end if

            do idim = 1, ndim
                gas_force_code(idim) = gas_force_code(idim) + Fp_code(idim) * mom_fact_gas * rt_pressBoost
            end do

            ! --- Dust Bin Forces ---
            if (ndust > 0) then
                do ii = 1, ndust
                    mom_fact_dust = opacity_dust_code(ii, igroup) * group_egy_code(igroup)

                    if (rt_isoPress .and. .not. (rt_isIR .and. igroup == iIR)) then
                        if (fluxMag_code > 0d0) then
                            mom_fact_dust = mom_fact_dust * Np_code / fluxMag_code
                        else
                            mom_fact_dust = 0d0
                        end if
                    else
                        mom_fact_dust = mom_fact_dust / rt_c_code
                    end if

                    do idim = 1, ndim
                        dust_force_code(idim, ii) = dust_force_code(idim, ii) + Fp_code(idim) * mom_fact_dust * rt_pressBoost
                    end do
                end do
            end if

            ! --- PAH Bin Forces ---
            if (npah > 0) then
                do ii = 1, npah
                    mom_fact_pah = opacity_pah_code(ii, igroup) * group_egy_code(igroup)

                    if (rt_isoPress .and. .not. (rt_isIR .and. igroup == iIR)) then
                        if (fluxMag_code > 0d0) then
                            mom_fact_pah = mom_fact_pah * Np_code / fluxMag_code
                        else
                            mom_fact_pah = 0d0
                        end if
                    else
                        mom_fact_pah = mom_fact_pah / rt_c_code
                    end if

                    do idim = 1, ndim
                        pah_force_code(idim, ii) = pah_force_code(idim, ii) + Fp_code(idim) * mom_fact_pah * rt_pressBoost
                    end do
                end do
            end if

        end do

    end subroutine compute_gas_dust_radpressure_force

end module dust_radpressure_module
