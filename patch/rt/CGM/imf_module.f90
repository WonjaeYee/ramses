MODULE imf_module
   ! Module for dealing with IMF calculations
   use amr_parameters, only: dp, imf_maker
   implicit none

   public get_N_supernova_in_step, get_effective_num_sn, get_upper_slope_prgomet, &
          IMFnum_between_var, SNIa_var_IMF, agemass_schaerer_1993, get_upper_slope_generic
CONTAINS
    SUBROUTINE get_N_supernova_in_step(m0,m1,m2,a1,a2,m_initial,mlow,mhi,NSNII)
        implicit none
        real(dp), intent(in)::m0,m1,m2,a1,a2,m_initial,mlow,mhi
        real(dp), intent(out)::NSNII
        real(dp)::get_N_supernova
        real(dp)::smooth_factor, total_mass

        NSNII = 0.d0
        if (mlow.lt.m1) then
            write(*,*) "Lower mass is below break mass -- integral won't work"
            return
        end if

        ! Get the power-law transition between the two regimes
        smooth_factor = (m1**a1) / (m1**a2)

        ! Calculate the total mass in the IMF
        total_mass = (1.d0/(a1 + 2.d0)) * ((m1**(a1 + 2.d0)) - (m0**(a1 + 2.d0)))
        total_mass = total_mass + ( (smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (m1**(a2 + 2.d0))) )

        ! Integrate the number between the interval
        NSNII = (smooth_factor/(a2 + 1.d0)) * ((mhi**(a2 + 1.d0)) - (mlow**(a2 + 1.d0)))

        ! Normalize by the total mass of the particle
        NSNII = NSNII * (m_initial/total_mass)

    END SUBROUTINE get_N_supernova_in_step

    SUBROUTINE IMFnum_between_var(m0,m1,m2,a1,a2,m_initial, mlow, mhi, Nstars)
        !Calculates the number of stars between two masses
        implicit none
        real(dp), intent(in)::m0,m1,m2,a1,a2,m_initial,mlow,mhi
        real(dp), intent(out)::Nstars
        real(dp)::IMFnum_between
        real(dp)::smooth_factor, total_mass, mlowa,mhia

        ! Initialize
        Nstars = 0.d0

        ! Enforce bounds
        mlowa = MIN(MAX(mlow,m0),m2)
        mhia =  MIN(MAX(mhi,m0),m2)

        ! Get the power-law transition between the two regimes
        smooth_factor = (m1**a1) / (m1**a2)

        ! Calculate the total mass in the IMF
        total_mass = (1.d0/(a1 + 2.d0)) * ((m1**(a1 + 2.d0)) - (m0**(a1 + 2.d0)))
        total_mass = total_mass + ((smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (m1**(a2 + 2.d0))))

        ! Perform the integral
        if (mlowa.lt.m1.and.mhia.le.m1) then
            Nstars = (1.d0/(a1 + 1.d0)) * ((mhia**(a1 + 1.d0)) - (mlowa**(a1 + 1.d0)))
        else if (mlowa.lt.m1.and.mhia.gt.m1) then
            Nstars = (1.d0/(a1 + 1.d0)) * ((m1**(a1 + 1.d0)) - (mlowa**(a1 + 1.d0)))
            Nstars = Nstars + ( (smooth_factor/(a2 + 1.d0)) * ((mhia**(a2 + 1.d0)) - (m1**(a2 + 1.d0))) )
        else if (mlowa.ge.m1.and.mhia.gt.m1) then
            Nstars = (smooth_factor/(a2 + 1.d0)) * ((mhia**(a2 + 1.d0)) - (mlowa**(a2 + 1.d0)))
        end if

        ! Renormalize
        Nstars = Nstars * (m_initial/total_mass)

        ! Enforce not negative
        Nstars = MAX(Nstars,0.d0)
    END SUBROUTINE IMFnum_between_var

    SUBROUTINE get_effective_num_sn(mett,smass,num_sn_eff_in,num_sn_eff_out,n_normal_SNII,n_hypernova)
    !---------------------------------------
    ! Computes the effective number of SN for a given
    ! stellar metallicity following recommendations from chiaki
    ! num_sn_eff_out --> effective total amount of energy units 10^51 erg
    ! n_hypernova --> total number of hypernova
    ! n_normal_SNII --> number of normal SN
    use random
    use pm_commons, only:localseed
    use amr_parameters, only:variable_energy_SN,dp
    implicit none
    REAL(dp), intent(out) :: num_sn_eff_out,n_normal_SNII,n_hypernova
    REAL(dp), intent(in) :: mett, smass, num_sn_eff_in
    REAL(dp):: A2 = 9.81113180d2
    REAL(dp):: A1 = -4.68848594d1
    REAL(dp):: A0 = 5.54330266d-1
    REAL(dp):: B2 = -1.90954774d-2
    REAL(dp):: B1 = 2.58793970d0
    REAL(dp):: B0 = -4.19597990d1
    REAL(dp):: loc_met, loc_smass, hn_fraction, hn_energy
    REAL(dp):: loc_rand
    INTEGER(kind=4):: i

    ! Initialize number of SN to the input number and set HN to 0
    num_sn_eff_out = num_sn_eff_in
    n_normal_SNII  = num_sn_eff_in
    n_hypernova    = 0.d0

    ! If we are below 20 Msun, no hypernova --> normal SN only
    if (smass.le.20.d0) return

    ! If hypernova are turned off
    if (.not.variable_energy_SN) then
        if (smass.gt.25.d0) then
            num_sn_eff_out = 0.001d0*num_sn_eff_in ! --> turn off SN for stars more massive than 25Msun
        end if
        !Otherwise it's a normal SN
        return
    end if

    loc_met   = MIN(MAX(mett,0.001d0),0.02d0) ! Place bounds on metallicity
    loc_smass = MIN(MAX(smass,25.d0),50.d0)  ! Place bounds on the exploding star particle

    ! Get the hypernova fraction for the given metallicity
    ! Should interpolate between 50% at low Z and 1% at high Z
    hn_fraction = (A2 * loc_met * loc_met) + (A1 * loc_met) + A0

    ! Get the hypernova energy
    ! Should interpolate between 10 and 40
    hn_energy = (B2 * loc_smass * loc_smass) + (B1 * loc_smass) + B0

    ! Calculate a random number between 0 and 1
    ! if it's less than hn_fraction then return the hypernova value
    ! otherwise return the normal value
    num_sn_eff_out = 0.d0
    n_normal_SNII  = 0.d0
    n_hypernova    = 0.d0
    do i=1,NINT(num_sn_eff_in)
        call ranf(localseed, loc_rand)
        if (loc_rand.le.hn_fraction) then
            num_sn_eff_out = num_sn_eff_out + hn_energy
            n_hypernova = n_hypernova + 1.d0
        else
            ! No strong SN at M>25 Msun
            if (smass.gt.25.d0) then
                num_sn_eff_out = num_sn_eff_out + 0.001d0 ! Set to a very low energy
            else
                num_sn_eff_out = num_sn_eff_out + 1.d0    ! Set to 10^51 erg
            end if
            n_normal_SNII = n_normal_SNII + 1.d0         ! Count the number of SN regardless for metal enrichment
        end if
    end do

    END SUBROUTINE get_effective_num_sn

    SUBROUTINE get_upper_slope_prgomet(FE_over_H,alpha_2)
        ! Returns the upper mass slope for an IMF following 
        ! https://watermark.silverchair.com/stac1074.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA1gwggNUBgkqhkiG9w0BBwagggNFMIIDQQIBADCCAzoGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM4j8SqtL3zsAPc7tDAgEQgIIDC5biKL76XDldAq5ULjL6MsR8NgeGn51QuHUdZe5vlOlkQmrEYLU8xexYhjl60HOV_P-6STTWN4TNAzfRjoqobtvXK4Xa_w_J6g1wDqM942Kz9y_wkqhpPZVKcLxtjLYJ31hw27UuZ5Xr7DABzkDlDiMMTvOgpvGvbRqJvswEnZRvKmn-idhzlyUsOD8KNT8rKuwRrnOiyTOBD-p5t_N44AgcMx56F1C347aHBf-1YcLZ_YjilpD-iuE8b9EDL0gRqk_kegKW3IrSJOe2q2cgoFf1T-F2jlFPmh9XypVVVrBJXZNl_9zEfTtHUHGkSGfOzP8DH3zZifAz_z0SEiNg10rFVc-SNQ0vW45cYTVugzqFkYx25AaI03tOl6Kqlh6FLMbEGYqMT3C03ocG8LMkm0_Ro-IZBaBDhlzK-Tc0Bofw9rUVm9892Cfn9J1LW_5Mnym0KuoyiQDSNPt-sB8OV_35LOyxKKrO42Y27dASrIlimTSOtwgsao8OjjNA62ivWxelBaNfJkruuVVFSHf-KyK8K7An6BIvNVnNp9Nd_hSSkKvecLzxTGdj5cnqsh4Uf0tOJ-G7MtYMf0JmGOuwEZ2i5g8UIJcC1_Cq1PnAj0OhvTCGLxQ6nSF8oy7SdqBC_zqrrYFZbLLebLTTSqPaQwj583S6rUW3iyxX0rYJwZnWRk75U1chT0cnD5HeRLRvS0hVUcpRgLFyS_9JtDawqk0WmkEtqXhVmkAbXiMmR5qfNY5uL95DG8_cBewPoulmMd48RtmS6Q6nKVOKL_PdeZn5bCV5u7NnO6Ua4Nz9TO8EUV7hQbX0oDJuT3f2GHz-RIteMWA8SA7DlW9na3tayh4iZAJUy6ExCpDqhF57vTd6nQnteWhg-_jBtSMA9YDhgpNzcMFGp-Ddckei99rXgGuaVknzo7Q1fxPF4xr6UK4Q2GOZzri2tyQD_CCtyJ1mAiG20qcwH2iqywr8mIDsAm0NiFOnRoqCOhJayp3v1aDaPr-2xQPBet9lmSfkFWCUdf3tyQ-SyZ03iuo6
        ! Equation 2
        implicit none
        REAL(dp), intent(in) :: FE_over_H
        REAL(dp), intent(out) :: alpha_2

        alpha_2 = (0.5d0 * FE_over_H) + 2.35
        if (FE_over_H.le.-3.d0) alpha_2 = 0.85d0
        alpha_2 = -1.d0 * alpha_2
    END SUBROUTINE get_upper_slope_prgomet

    SUBROUTINE get_upper_slope_generic(FE_over_H,redshift,nH,TK,alpha_2)
        ! This is a generic routine to get the upper mass slope of the IMF
        ! Edit here for a general function
        ! Usually depends on some combination of redshift, metallicity, and gas density
        implicit none
        REAL(dp), intent(in) :: FE_over_H, redshift, nH, TK
        REAL(dp), intent(out) :: alpha_2
        REAL(dp):: marks_rho_fid = 24713.197d0 ! Fiducial density

        alpha_2 = 2.3d0 ! Default the upper slope to salpeter

        if (TRIM(imf_maker)=='marks') then
           ! SEE EQN 15 of Marks+ 2012
           ! NO CMB DEPENDENCE
           alpha_2 = (0.0572d0 * FE_over_H) - &
                     (0.4072d0 * LOG10(nH/marks_rho_fid)) + &
                     1.9383d0
           ! Set boundaries
           alpha_2 = MAX(MIN(alpha_2,2.3d0),0.8d0)
        else if (TRIM(imf_maker)=='chon') then
           ! This is a model loosely based off Chon et al.
           ! The idea is that the slope gets shallower with
           ! metallicity and CMB
           alpha_2 = 2.3d0 + 0.45d0 * LOG10(MAX(10.d0**FE_over_H,1.d-5)/0.1d0)
           alpha_2 = MAX(MIN(alpha_2,2.3d0),0.5d0)
           alpha_2 = alpha_2 - (1.8d0 * MAX(MIN(redshift - 5.d0,15.d0),0.d0)/15.d0)
           alpha_2 = MAX(MIN(alpha_2,2.3d0),0.5d0)
        else if (TRIM(imf_maker)=='jeans') then
           ! Based on the idea that the IMF slope goes with jeans mass
           ! which scales as n^-1/2 and T^3/2
           alpha_2 = 2.3d0 * ((MAX(TK,30.d0)/30.d0)**(-1.5d0)) * ((MAX(nH,100.d0)/100.d0)**0.5d0)
           alpha_2 = MAX(MIN(alpha_2,2.6d0),0.5d0)
        end if

    END SUBROUTINE get_upper_slope_generic

    SUBROUTINE SNIa_var_IMF(m0,m1,m2,a1,a2, t1, t2, NSNIa, MMin)
    !---------------------------------------
    use amr_commons
    implicit none
    REAL(dp), intent(out):: NSNIa
    REAL(dp), intent(in) :: t1, t2, MMin
    real(dp), intent(in) ::m0,m1,m2,a1,a2
    REAL(dp)::smooth_factor,total_mass,total_mass_over_8,total_mf_below_8

    !---  Delay Time Distribution (DTD). Maoz & Graur (2017).
    !---  Literature normalisations
    !---  2.6d-13 Ia/yr/Msun = field DTD
    !---  1.3d-13 Ia/yr/Msun = old field DTD (Graur et al. 2014)
    !---  Greater values of 4d-13-8d-13 Ia/yr/Msun = compatible with cluster DTD

    ! Harley minor edits to Oscar's routine
    ! Makes sure we only call the routine when 8Msun stars 
    ! are off the main sequence
    !
    ! The difference here compared to before is we apply this to the mass of the
    ! Particle that is below 8 Msun. There is an extra normalization factor
    ! For only the mass below 8 Msun in the IMF

    ! initialize
    NSNIa = 0.d0

    ! If no stars below 8 Msun then return
    if (m0.ge.8.d0) return

    ! Get the power-law transition between the two regimes
    smooth_factor = (m1**a1) / (m1**a2)

    ! Calculate the total mass in the IMF
    total_mass = (1.d0/(a1 + 2.d0)) * ((m1**(a1 + 2.d0)) - (m0**(a1 + 2.d0)))
    total_mass = total_mass + ((smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (m1**(a2 + 2.d0))))

    ! Calculate the mass above 8 Msun
    total_mass_over_8 = (smooth_factor/(a2 + 2.d0)) * ((m2**(a2 + 2.d0)) - (8.d0**(a2 + 2.d0)))
    total_mf_below_8 = (total_mass - total_mass_over_8) / total_mass_over_8

    if (Mmin .lt. 8.d0) then                     !------ set by MS lifetime of 8 Msun stars
        NSNIa = total_mf_below_8 * Ia_rate*(t1/1d9)**(-1.12)*(t2 - t1)  !*mpb(ind_part(j))/scale_m ---- is normalised outside
    else
        NSNIa = 0.0
    end if

    END SUBROUTINE SNIa_var_IMF

    SUBROUTINE agemass_schaerer_1993(time, mass)
    !----------------------------------------------
    ! Harley fit to the age mass relation from schaerer 1993
    ! We used webplot digitizer to pull the points from 
    ! Figure 1 of https://articles.adsabs.harvard.edu/pdf/1993A%26AS...98..523S
    ! For now we ignore dependence on metallicity as it's small
    ! Time --> Myr
    ! Mass --> Msun
    implicit none
    real(dp), intent(in)::time
    real(dp), intent(out)::mass
    real(dp)::log_t_yr

    ! Convert the time to log10 years
    log_t_yr = LOG10(time)
    
    ! Fix age low problem --> always return the top mass if below this value
    if (log_t_yr.le.6.411918592721d0) then
        mass = 120.d0
        return
    end if

    ! 5th order polynomial fit to the age mass relation
    ! SHould be monotically decreasing with time
    mass = 0.d0
    mass = mass + (-1.39513319d-2 * (log_t_yr**5.d0))
    mass = mass + (6.11446524d-1 * (log_t_yr**4.d0))
    mass = mass + (-1.06823084d1 * (log_t_yr**3.d0))
    mass = mass + (9.30336650d1 * (log_t_yr**2.d0))
    mass = mass + (-4.04474842d2 * log_t_yr)
    mass = mass + 7.04352171d2

    
    mass = 10.d0**mass

    END SUBROUTINE agemass_schaerer_1993

END MODULE imf_module
