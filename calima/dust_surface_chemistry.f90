module dust_surface_chemistry
    use amr_parameters
    use hydro_commons, only:nmetals
    use cooling_module, only: X, kB, mH
    use constants, only: twopi,pi
    use dust_commons

    implicit none

    ! Parameters from Cazaux & Spaans (2004)
    ! (https://iopscience.iop.org/article/10.1086/422087/pdf)
    ! TODO: Code for all dust compositions
    real(dp),dimension(2),parameter :: EH2 = (/540d0,340d0/) ! in K
    real(dp),dimension(2),parameter :: mu = (/4d-1,3d-1/)
    real(dp),dimension(2),parameter :: Es = (/250d0,200d0/) ! in K
    real(dp),dimension(2),parameter :: EHp = (/800d0,650d0/) ! in K
    real(dp),dimension(2),parameter :: EHc = (/3d4,3d4/) ! in K
    real(dp),dimension(2),parameter :: nuH2 = (/3d12,2d12/) ! in 1/s
    real(dp),dimension(2),parameter :: nuHc = (/2d13,1d13/) ! in 1/s

    real(dp),parameter :: Ns = 2D15 ! Fixed number of sites per cm2 on the surface of the grain

    contains

    function h2_sticking_coef(Tgas,Td)
        ! Sticking coefficient from Hollenbach & McKee (1979) - Eq 3.7
        ! (https://articles.adsabs.harvard.edu/pdf/1979ApJS...41..555H)
        implicit none
        real(dp) :: Tgas,Td
        real(dp) :: h2_sticking_coef

        h2_sticking_coef = 1D0 / (1D0 + 4D-1 * sqrt((Tgas+Td)/1D2) + 2D-1 * (Tgas/1D2) + 8D-2 * (Tgas/1D2)**2D0)
    end function h2_sticking_coef

    function beta_h2(Tgas,Td,dust_index)
        ! Desorption rate of H2 from Cazaux & Tielens (2002) - Page 2, beginning of leftmost last paragraph
        ! https://iopscience.iop.org/article/10.1086/342607/pdf
        implicit none
        real(dp) :: Tgas,Td
        real(dp) :: beta_h2
        integer   :: dust_index
        ! NOTE: EH2 is in K, not erg, so the kB factor is not needed
        beta_h2 = nuH2(dust_index) * exp(-EH2(dust_index) / (Td))
    end function beta_h2

    function high_temp_correction(Tgas,Td,F,dust_index)
        ! High temperature correction from Cazaux & Tielens (2002) - Eq 16 in the Erratum
        ! (https://iopscience.iop.org/article/10.1086/342607/pdf)
        implicit none
        real(dp) :: Tgas,Td,F
        real(dp) :: high_temp_correction
        integer   :: dust_index

        real(dp) :: a1,a2,a3

        a1 = nuHc(dust_index) / (2D0 * F)
        a2 = exp(-1.5D0 * EHc(dust_index) / (Td))
        a3 = (1D0 + sqrt( (EHc(dust_index) - Es(dust_index)) / (EHp(dust_index) - Es(dust_index)) ))**2D0

        high_temp_correction = 1D0 / (1D0 + (a1 * a2 * a3))
    end function high_temp_correction

    function beta_hp_over_alphapc(Tgas,Td,dust_index)
        ! Eq 17 from Cazaux & Tielens (2002)
        implicit none
        real(dp) :: Tgas,Td
        real(dp) :: beta_hp_over_alphapc
        integer   :: dust_index

        real(dp) :: a1,a2

        a1 = (1D0 + sqrt( (EHc(dust_index) - Es(dust_index)) / (EHp(dust_index) - Es(dust_index)) ))**2D0
        a2 = exp(-1D0 * Es(dust_index) / (Td))
        beta_hp_over_alphapc =  2.5D-1 * a1 * a2
    end function beta_hp_over_alphapc

    function h_flux(nH,vH)
        ! Eqn 1 of Cazaux & Spaans (2004) 
        ! (https://iopscience.iop.org/article/10.1086/422087/pdf)

        ! One can see from Figure 1 of https://iopscience.iop.org/article/10.1086/342607/pdf 
        ! that the recombination efficiency is essentially independent of flux over 10 orders 
        ! of magnitude. Taking a constant Ns is probably ok
        implicit none
        real(dp) :: nH,vH
        real(dp) :: h_flux

        h_flux = nH * vH / Ns
    end function h_flux

    function recombination_efficiency(Tgas,Td,nH,vH,F,dust_index)
        ! Recombination effiency from Cazaux & Tielens (2002) -- Eq 15
        ! (https://iopscience.iop.org/article/10.1086/342607/pdf)
        implicit none
        real(dp) :: Tgas,Td,nH,vH,F
        real(dp) :: recombination_efficiency
        integer   :: dust_index

        real(dp) :: a1,a2

        a1 = 1D0 / (1D0 + (mu(dust_index)*F) / (2D0*beta_h2(Tgas,Td,dust_index)) + beta_hp_over_alphapc(Tgas,Td,dust_index))
        a2 = high_temp_correction(Tgas,Td,F,dust_index)
        recombination_efficiency = a1 * a2
    end function recombination_efficiency

    function h2_formation_rate(nH,Tgas,G0_total,rho_dust,T_dust,rho_pah,fcharge_pahs)
        implicit none
        real(dp) :: nH,Tgas,G0_total
        real(dp),dimension(1:ndust) :: rho_dust,T_dust
        real(dp),dimension(1:npah),optional :: rho_pah
        real(dp),dimension(:,:),optional :: fcharge_pahs

        integer :: j,ilow,ihigh
        real(dp) :: R_H2,vH,F,sdust
        real(dp) :: h2_formation_rate
        real(dp) :: chi,x_frac,R_PAHs

        h2_formation_rate = 0d0

        if (H2ondust) then
            ! Formation rate of H2 onto dust grains from Cazaux & Spaans (2004)
            ! (https://iopscience.iop.org/article/10.1086/422087/pdf)
            vH = sqrt(2D0 * kB * Tgas / mH) ! thermal velocity (assuming Mawell-Boltzmann distribution)
            F = h_flux(nH,vH)

            
            ! Add the contribution from each dust grain
            do j = 1, ndust
                sdust = (rho_dust(j)/dustbins_props(j)%mgrain) * twopi * (dustbins_props(j)%asize*1D-4)**2D0
                R_H2 = sdust * recombination_efficiency(Tgas,T_dust(j),nH,vH,F,dustbins_props(j)%interact_group)
                h2_formation_rate = h2_formation_rate + R_H2 * h2_sticking_coef(Tgas,T_dust(j))
            end do
            h2_formation_rate = 5D-1 * nH * vH * h2_formation_rate ! in cm-3*s-1
        end if

        
        if (present(rho_pah) .and. dust_pahs .and. H2onpah) then
            select case(pah_h2_model)
            case ('LePage09')
                ! Fitting parameters to power law H2 rate from Le Page et al. (2009) on PAHs
                ! (https://ui.adsabs.harvard.edu/abs/2009ApJ...704..274L/abstract)
                ! The mechanism involves the chemical trapping of H atoms on the periphery of the PAH
                ! carbon skeleton and the subsequent release of H2 through dissociative recombination
                ! of the hydrogenated ion with an electron.
                ! Additionnaly we include a scaling with the number of carbon atoms, as it is seen
                ! in the work of Le Page et al. (2009)
                ! Convert G0 to Draine ISRF units
                chi = 1.69d0 * G0_total
                x_frac = nH / chi
                do j=1,npah
                    R_PAHs = 10**(-15.269d0 - 1.098d0*log10(x_frac**0.916d0 + 11.090d0)) &
                                & * (50d0/dble(pahbins_props(j)%nc)) * nH * (rho_pah(j)/(pahbins_props(j)%nc*mC_amu*amu2g)) ! [cm-3*s-1]
                    h2_formation_rate = h2_formation_rate + R_PAHs
                end do
            case ('RM2026')
                ! Compute the H2 formation rate following CALIMA (Rodriguez Montero+2026)
                h2_formation_rate = h2_formation_rate + pah_h2_formation_rate(rho_pah,fcharge_pahs,G0_total,nH,Tgas)
            end select
        end if
    end function h2_formation_rate

    subroutine compute_dehydrogenated_fraction(G0,nH,f_dh)
        implicit none
        real(dp), intent(in) :: G0,nH
        real(dp), intent(inout) :: f_dh
        real(dp) :: lognH,logG0,distance,y_intercept
        lognH = log10(nH)
        logG0 = log10(G0)
        ! Compute the G0 and nH relation based on the fitting
        ! to the results of Montillaud et al. (2013) for circumcoronene
        y_intercept = 1.542d0 * logG0 + 0.646d0

        ! Distance from point to the line
        distance = abs(y_intercept - lognH) / 1.838d0 ! denominator is sqrt(1+1.54184841^2)

        if (lognH .ge. y_intercept) then
            f_dh = 0.5d0 - 0.5d0/(1d0+1d-1*distance**(-2d0))
        else
            f_dh = 0.5d0 + 0.5d0/(1d0+1d-1*distance**(-2d0))
        end if
        f_dh = max(f_dh,0d0)
    end subroutine compute_dehydrogenated_fraction

    subroutine compute_superhydrogenated_fraction(G0,nH,f_sh)
        implicit none
        real(dp), intent(in) :: G0,nH
        real(dp), intent(inout) :: f_sh
        real(dp) :: lognH,logG0,distance,y_intercept
        lognH = log10(nH)
        logG0 = log10(G0)
        ! Compute the G0 and nH relation based on the fitting
        ! to the results of Andrews et al. (2016) for
        ! circumcircumcoronene (Nc=96)
        y_intercept = 0.995d0 * lognH - 1.945d0

        ! Distance from point to the line
        distance = abs(y_intercept - logG0) / 1.410d0 ! denominator is sqrt(1+0.99466057^2)

        if (logG0 .ge. y_intercept) then
            f_sh = 0.5d0 - 0.5d0/(1d0+5d-2*distance**(-2d0))
        else
            f_sh = 0.5d0 + 0.5d0/(1d0+5d-2*distance**(-2d0))
        end if
        f_sh = max(f_sh,0d0)
    end subroutine compute_superhydrogenated_fraction

    function pah_h2_formation_rate(rho_pah,fcharge,G0,nH,Tgas)
        ! H2 FORMATION RATE ON PAHS (CALIMA) - Rodriguez Montero et al. (2026)
        ! The modelling of H2 formation provided by this function considers the H2
        ! abstraction by the Eley-Rideal mechanism in super-hydrogenated PAHs
        ! (affecting both small and large PAH clusters), and the removal of H2 due
        ! to UV photo-dissociation in normally hydrogenated and partially
        ! de-hydrogenated PAHs (only contemplated for small PAHs).
        ! For further details, check the explanation of the full model in 
        ! Rodriguez Montero et al. (2024).
        ! WARNING: This does not work for more than two PAH sizes!
        implicit none
        real(dp),intent(in) :: G0,nH,Tgas
        real(dp),dimension(1:npah),intent(in) :: rho_pah
        real(dp),dimension(:,:),intent(in) :: fcharge
        real(dp) :: pah_h2_formation_rate

        integer :: i, nstates
        real(dp) :: k_ER,kH2,f_dehydro,f_pdhydro,f_suphydro
        real(dp) :: RPAH
        pah_h2_formation_rate = 0d0
        RPAH = 0d0

        ! 1. Compute the Eley-Rideal rate and H2 dissociation rate
        k_ER = 8.7d-13 * sqrt(Tgas/1d2) * nH ! [s-1]
        kH2 = 3.428d-16 * (1.69d0 * G0)      ! [s-1] G0 needs to go from Habing to Draine ISRF

        ! 2. First compute for the small PAHs
        ! 3. Determine the fractions of de-hydrogenated and super-hydrogenated PAHs
        call compute_dehydrogenated_fraction(G0,nH,f_dehydro)
        f_pdhydro = f_dehydro / (1d0 + exp((f_dehydro-0.99d0)/8d-4))
        call compute_superhydrogenated_fraction(G0,nH,f_suphydro)

        ! 3. Add contribution for each of the PAH charges
        nstates = pahbins_props(1)%ncharge_states
        do i = 1, nstates
            if (i .le. 2) then
                RPAH = RPAH + (kH2 * f_pdhydro + k_ER * f_suphydro) * fcharge(i,1)
            else
                RPAH = RPAH + kH2 * f_pdhydro * fcharge(i,1)
            end if
        end do
        pah_h2_formation_rate = pah_h2_formation_rate + RPAH * rho_pah(1)/(dble(pahbins_props(1)%nc)*mC_amu*amu2g)

        RPAH = 0d0
        ! 2. Now compute for the large PAHs (if they are present)
        if (npah .ge. 2) then
            nstates = pahbins_props(2)%ncharge_states
            do i = 1, nstates
                if (i .le. 2) then
                    RPAH = RPAH + k_ER * f_suphydro * fcharge(i,2)
                end if
            end do
            pah_h2_formation_rate = pah_h2_formation_rate + RPAH * rho_pah(2)/(dble(pahbins_props(2)%nc)*mC_amu*amu2g)
        end if
    end function pah_h2_formation_rate
end module dust_surface_chemistry