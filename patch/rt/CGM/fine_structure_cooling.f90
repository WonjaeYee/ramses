
MODULE fscool
    use amr_parameters, only: dp, fs_cooling_dir
    use constants, ONLY: c_cgs, hplanck, kB
    implicit none

    private   ! default
    public init_fine_structure_tables, OI_fine_structure, OIII_fine_structure, CI_fine_structure, &
           CII_fine_structure, NII_fine_structure, SiI_fine_structure, SiII_fine_structure, &
           FeI_fine_structure, FeII_fine_structure, SI_fine_structure, NeII_fine_structure

    real(dp),allocatable, dimension(:,:,:)   :: fs_cool_tab
    real(dp), dimension(1:160) :: table_temp = (/ 1.   , 1.025, 1.05 , 1.075, 1.1  , 1.125, 1.15 , 1.175, 1.2  , &
                                              & 1.225, 1.25 , 1.275, 1.3  , 1.325, 1.35 , 1.375, 1.4  , 1.425, &
                                              & 1.45 , 1.475, 1.5  , 1.525, 1.55 , 1.575, 1.6  , 1.625, 1.65 , &
                                              & 1.675, 1.7  , 1.725, 1.75 , 1.775, 1.8  , 1.825, 1.85 , 1.875, &
                                              & 1.9  , 1.925, 1.95 , 1.975, 2.   , 2.025, 2.05 , 2.075, 2.1  , &
                                              & 2.125, 2.15 , 2.175, 2.2  , 2.225, 2.25 , 2.275, 2.3  , 2.325, &
                                              & 2.35 , 2.375, 2.4  , 2.425, 2.45 , 2.475, 2.5  , 2.525, 2.55 , &
                                              & 2.575, 2.6  , 2.625, 2.65 , 2.675, 2.7  , 2.725, 2.75 , 2.775, &
                                              & 2.8  , 2.825, 2.85 , 2.875, 2.9  , 2.925, 2.95 , 2.975, 3.   , &
                                              & 3.025, 3.05 , 3.075, 3.1  , 3.125, 3.15 , 3.175, 3.2  , 3.225, &
                                              & 3.25 , 3.275, 3.3  , 3.325, 3.35 , 3.375, 3.4  , 3.425, 3.45 , &
                                              & 3.475, 3.5  , 3.525, 3.55 , 3.575, 3.6  , 3.625, 3.65 , 3.675, &
                                              & 3.7  , 3.725, 3.75 , 3.775, 3.8  , 3.825, 3.85 , 3.875, 3.9  , &
                                              & 3.925, 3.95 , 3.975, 4.   , 4.025, 4.05 , 4.075, 4.1  , 4.125, &
                                              & 4.15 , 4.175, 4.2  , 4.225, 4.25 , 4.275, 4.3  , 4.325, 4.35 , &
                                              & 4.375, 4.4  , 4.425, 4.45 , 4.475, 4.5  , 4.525, 4.55 , 4.575, &
                                              & 4.6  , 4.625, 4.65 , 4.675, 4.7  , 4.725, 4.75 , 4.775, 4.8  , &
                                              & 4.825, 4.85 , 4.875, 4.9  , 4.925, 4.95 , 4.975 /)

    CONTAINS
    SUBROUTINE init_fine_structure_tables()
        implicit none
        character(LEN=256),dimension(1:27)::fine_struct_cooling
        integer::i,j,ntemp=160
        logical::ok
   
        ! Files with collisional partner data
        ! rate_electron,rate_HII,rate_HeII,rate_HeIII,rate_HI,rate_HeI,rate_oH2,rate_pH2
        fine_struct_cooling(1)  = trim(fs_cooling_dir)//'CII_158um_rates.dat'
        fine_struct_cooling(2)  = trim(fs_cooling_dir)//'CI_609um_rates.dat'
        fine_struct_cooling(3)  = trim(fs_cooling_dir)//'CI_230um_rates.dat'
        fine_struct_cooling(4)  = trim(fs_cooling_dir)//'CI_370um_rates.dat'
        fine_struct_cooling(5)  = trim(fs_cooling_dir)//'NII_205um_rates.dat'
        fine_struct_cooling(6)  = trim(fs_cooling_dir)//'NII_76um_rates.dat'
        fine_struct_cooling(7)  = trim(fs_cooling_dir)//'NII_122um_rates.dat'
        fine_struct_cooling(8)  = trim(fs_cooling_dir)//'OI_63um_rates.dat'
        fine_struct_cooling(9)  = trim(fs_cooling_dir)//'OI_44um_rates.dat'
        fine_struct_cooling(10) = trim(fs_cooling_dir)//'OI_145um_rates.dat'
        fine_struct_cooling(11) = trim(fs_cooling_dir)//'NeII_13um_rates.dat'
        fine_struct_cooling(12) = trim(fs_cooling_dir)//'SiII_35um_rates.dat'
        fine_struct_cooling(13) = trim(fs_cooling_dir)//'OIII_88um_rates.dat'
        fine_struct_cooling(14) = trim(fs_cooling_dir)//'OIII_33um_rates.dat'
        fine_struct_cooling(15) = trim(fs_cooling_dir)//'OIII_52um_rates.dat'
        fine_struct_cooling(16) = trim(fs_cooling_dir)//'FeI_14um_rates.dat'
        fine_struct_cooling(17) = trim(fs_cooling_dir)//'FeI_24um_rates.dat'
        fine_struct_cooling(18) = trim(fs_cooling_dir)//'FeI_35um_rates.dat'
        fine_struct_cooling(19) = trim(fs_cooling_dir)//'FeII_26um_rates.dat'
        fine_struct_cooling(20) = trim(fs_cooling_dir)//'FeII_15um_rates.dat'
        fine_struct_cooling(21) = trim(fs_cooling_dir)//'FeII_35um_rates.dat'
        fine_struct_cooling(22) = trim(fs_cooling_dir)//'SiI_130um_rates.dat'
        fine_struct_cooling(23) = trim(fs_cooling_dir)//'SiI_45um_rates.dat'
        fine_struct_cooling(24) = trim(fs_cooling_dir)//'SiI_68um_rates.dat'
        fine_struct_cooling(25) = trim(fs_cooling_dir)//'SI_17um_rates.dat'
        fine_struct_cooling(26) = trim(fs_cooling_dir)//'SI_25um_rates.dat'
        fine_struct_cooling(27) = trim(fs_cooling_dir)//'SI_56um_rates.dat'

        ! Check that all files exist
        do i = 1,27
            inquire(FILE=fine_struct_cooling(i),exist=ok)
            if(.not.ok)then
                write(*,*) 'Cannot access fine structure cooling file '//TRIM(fine_struct_cooling(i))
                call clean_stop
            endif
        end do

        ! Allocate the arrays
        ! Emission line, Number of temperatures, collisional partner rates
        allocate(fs_cool_tab(1:27,1:ntemp,1:8))

        ! Read in the fine structure data
        do j = 1,27 ! Loop over emission lines
            open(unit=880,file=fine_struct_cooling(j),form='formatted')
                do i = 1,ntemp ! Loop over temperatures
                    read(880,*) fs_cool_tab(j,i,:)
                end do
            close(880)
        end do
    END SUBROUTINE init_fine_structure_tables
    
   FUNCTION interp_coll_1d(table_coll, temp)
      real(dp), dimension(:), intent(in) :: table_coll
      real(dp), intent(in) :: temp
      real(dp) :: interp_coll_1d
      real(dp) :: min_yield=1.d-40

      real(dp) :: y1, y2, f1, f2, a1, a2
      integer :: itemp

      ! Find temperature index.
      itemp = FLOOR( ((temp - 1.d0) / (4.975d0 - 1.d0)) * 160.d0) + 1

      ! ----------------
      ! Interpolate
      ! Store values around point.
      y1 = table_temp(itemp)
      y2 = table_temp(itemp + 1)

      f1 = MAX(table_coll(itemp),min_yield)
      f2 = MAX(table_coll(itemp + 1),min_yield)

      ! 1D interp
      a1 = 1.0 - ((table_temp(itemp) - y1)/(y2 - y1))
      a2 = 1.0 - ((y2 - table_temp(itemp))/(y2 - y1))

      interp_coll_1d = (f1*a1) + (f2*a2)

   END FUNCTION interp_coll_1d

    FUNCTION three_level(g_0, g_1, g_2, lam_10, lam_20, lam_21, A_10, A_20, A_21, &
                        z, T, n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                        idx_0, idx_1, idx_2)
        implicit none
        real(dp), intent(in)::g_0,g_1,g_2          ! degeneracy
        real(dp), intent(in)::A_10,A_20,A_21       ! s^-1
        real(dp), intent(in)::lam_10,lam_20,lam_21 ! microns
        real(dp), intent(in)::z                    ! redshift
        real(dp), intent(in)::n_ion,ne,nH,nHp,nHe,nHep,nHepp,nH2 ! Number densities
        real(dp), intent(in)::T                    ! temperature
        integer, intent(in)::idx_0, idx_1, idx_2   ! indicies of lines in the fine_struct_cooling array
        real(dp):: T_cmb,logT
        real(dp):: nu_10,nu_20,nu_21
        real(dp):: E_10,E_20,E_21
        real(dp):: B_01,B_02,B_12,B_10,B_20,B_21
        real(dp):: C_01,C_02,C_12,C_10,C_20,C_21
        real(dp):: B_nu_10,B_nu_20,B_nu_21
        real(dp):: q10_oH2, q20_oH2, q21_oH2
        real(dp):: q10_pH2, q20_pH2, q21_pH2
        real(dp):: q10_H, q20_H, q21_H
        real(dp):: q10_Hp, q20_Hp, q21_Hp
        real(dp):: q10_e, q20_e, q21_e
        real(dp):: q10_He, q20_He, q21_He
        real(dp):: q10_Hep, q20_Hep, q21_Hep
        real(dp):: q10_Hepp, q20_Hepp, q21_Hepp
        real(dp):: n_0, n_1, n_2
        real(dp):: cool_0, cool_1, cool_2
        real(dp):: heat_0, heat_1, heat_2
        real(dp):: n2_n1_top, n2_n1_bot, n2_n1
        real(dp):: n1_n0_top, n1_n0_bot, n1_n0
        real(dp):: three_level

        ! Initialize the result
        three_level = 1.d-50

        T_cmb = 2.725d0 * (1.d0 + z) ! CMB temperature

        ! Get log T and enforce bounds
        logT = LOG10(T)

        logT = MAX(logT,table_temp(1))
        !if (logT .lt. MAX(LOG10(T_cmb),1.d0)) logT = MAX(LOG10(T_cmb),1.d0)
        if (logT .gt. 5.d0) logT = 5.d0

        nu_10 = c_cgs / (lam_10 * 1.d-4) ! Hz
        nu_20 = c_cgs / (lam_20 * 1.d-4) ! Hz
        nu_21 = c_cgs / (lam_21 * 1.d-4) ! Hz

        E_10 = (hplanck * c_cgs / (lam_10 * 1d-4)) / kB ! E/K (K)
        E_20 = (hplanck * c_cgs / (lam_20 * 1d-4)) / kB ! E/K (K)
        E_21 = (hplanck * c_cgs / (lam_21 * 1d-4)) / kB ! E/K (K)

        B_01 = A_10 * ((lam_10 * 1.d-4)**3.d0) * (g_1 / g_0) / (2.d0 * hplanck * c_cgs)
        B_02 = A_20 * ((lam_20 * 1.d-4)**3.d0) * (g_2 / g_0) / (2.d0 * hplanck * c_cgs)
        B_12 = A_21 * ((lam_21 * 1.d-4)**3.d0) * (g_2 / g_1) / (2.d0 * hplanck * c_cgs)

        B_10 = (g_0 / g_1) * B_01
        B_20 = (g_0 / g_2) * B_02
        B_21 = (g_1 / g_2) * B_12

        ! CMB black body spectrum
        B_nu_10 = (2.d0 * hplanck * (nu_10**3.d0) / (c_cgs**2.d0)) / (EXP(hplanck * nu_10 / (kB * T_cmb)) - 1.d0)
        B_nu_20 = (2.d0 * hplanck * (nu_20**3.d0) / (c_cgs**2.d0)) / (EXP(hplanck * nu_20 / (kB * T_cmb)) - 1.d0)
        B_nu_21 = (2.d0 * hplanck * (nu_21**3.d0) / (c_cgs**2.d0)) / (EXP(hplanck * nu_21 / (kB * T_cmb)) - 1.d0)

        ! All collison strengths are in cm^3 s^-1
        ! For collisions with o-H2
        !rate_electron,rate_HII,rate_HeII,rate_HeIII,rate_HI,rate_HeI,rate_oH2,rate_pH2
        q10_oH2 = interp_coll_1d(fs_cool_tab(idx_0,:,7), logT)
        q20_oH2 = interp_coll_1d(fs_cool_tab(idx_1,:,7), logT)
        q21_oH2 = interp_coll_1d(fs_cool_tab(idx_2,:,7), logT)

        ! For collisions with p-H2
        q10_pH2 = interp_coll_1d(fs_cool_tab(idx_0,:,8), logT)
        q20_pH2 = interp_coll_1d(fs_cool_tab(idx_1,:,8), logT)
        q21_pH2 = interp_coll_1d(fs_cool_tab(idx_2,:,8), logT)

        ! For collisions with H
        q10_H = interp_coll_1d(fs_cool_tab(idx_0,:,5), logT)
        q20_H = interp_coll_1d(fs_cool_tab(idx_1,:,5), logT)
        q21_H = interp_coll_1d(fs_cool_tab(idx_2,:,5), logT)

        ! For collisions with H+
        q10_Hp = interp_coll_1d(fs_cool_tab(idx_0,:,2), logT)
        q20_Hp = interp_coll_1d(fs_cool_tab(idx_1,:,2), logT)
        q21_Hp = interp_coll_1d(fs_cool_tab(idx_2,:,2), logT)

        ! For collisions with e
        q10_e = interp_coll_1d(fs_cool_tab(idx_0,:,1), logT)
        q20_e = interp_coll_1d(fs_cool_tab(idx_1,:,1), logT)
        q21_e = interp_coll_1d(fs_cool_tab(idx_2,:,1), logT)

        ! For collisions with He
        q10_He = interp_coll_1d(fs_cool_tab(idx_0,:,6), logT)
        q20_He = interp_coll_1d(fs_cool_tab(idx_1,:,6), logT)
        q21_He = interp_coll_1d(fs_cool_tab(idx_2,:,6), logT)

        ! For collisions with He+
        q10_Hep = interp_coll_1d(fs_cool_tab(idx_0,:,3), logT)
        q20_Hep = interp_coll_1d(fs_cool_tab(idx_1,:,3), logT)
        q21_Hep = interp_coll_1d(fs_cool_tab(idx_2,:,3), logT)

        ! For collisions with He++
        q10_Hepp = interp_coll_1d(fs_cool_tab(idx_0,:,4), logT)
        q20_Hepp = interp_coll_1d(fs_cool_tab(idx_1,:,4), logT)
        q21_Hepp = interp_coll_1d(fs_cool_tab(idx_2,:,4), logT)

        ! Net collision strengths
        C_10 = (q10_e * ne) + (q10_H * nH) + (q10_Hp * nHp) + (q10_oH2 * 0.75d0 * nH2) + (q10_pH2 * 0.25d0 * nH2) + (q10_He * nHe) + (q10_Hep * nHep) + (q10_Hepp * nHepp)
        C_20 = (q20_e * ne) + (q20_H * nH) + (q20_Hp * nHp) + (q20_oH2 * 0.75d0 * nH2) + (q20_pH2 * 0.25d0 * nH2) + (q20_He * nHe) + (q20_Hep * nHep) + (q20_Hepp * nHepp)
        C_21 = (q21_e * ne) + (q21_H * nH) + (q21_Hp * nHp) + (q21_oH2 * 0.75d0 * nH2) + (q21_pH2 * 0.25d0 * nH2) + (q21_He * nHe) + (q21_Hep * nHep) + (q21_Hepp * nHepp)

        C_01 = C_10 * (g_1/g_0) * EXP(-1.d0 * E_10 / T)
        C_02 = C_20 * (g_2/g_0) * EXP(-1.d0 * E_20 / T)
        C_12 = C_21 * (g_2/g_1) * EXP(-1.d0 * E_21 / T)

        C_01 = C_01 + (B_01*B_nu_10)
        C_02 = C_02 + (B_02*B_nu_20)
        C_12 = C_12 + (B_12*B_nu_21)

        C_10 = C_10 + (B_10*B_nu_10)
        C_20 = C_20 + (B_20*B_nu_20)
        C_21 = C_21 + (B_21*B_nu_21)

        ! Analytic solution to the 2 level system (see paul goldsmith papers)
        ! Done this way to avoid numerical errors
        n2_n1_top = (C_12 * (C_01 + C_02)) + (C_02 * (A_10 + C_10))
        n2_n1_bot = ((A_21 + C_21 + C_20) * (C_01 + C_02)) - (C_20 * C_02)
        n2_n1 = n2_n1_top / n2_n1_bot

        n1_n0_top = ((A_21 + C_21 + C_20) * (C_01 + C_02)) - (C_20*C_02)
        n1_n0_bot = ((A_21 + C_21 + C_20) * (A_10 + C_10)) + (C_20*C_12)
        n1_n0 = n1_n0_top / n1_n0_bot

        n_0 = 1.d0 / (1.d0 + n1_n0 + (n2_n1 * n1_n0))
        n_1 = 1.d0 / (1.d0 + n2_n1 + (1.d0/n1_n0))
        n_2 = 1.d0 / (1.d0 + (1.d0/n2_n1) + ((1.d0/n2_n1)*(1.d0/n1_n0)))

        ! Cooling and heating rates
        cool_0 = (A_10 + (B_10 * B_nu_10)) * E_10 * kB * n_1 * n_ion
        cool_1 = (A_20 + (B_20 * B_nu_20)) * E_20 * kB * n_2 * n_ion
        cool_2 = (A_21 + (B_21 * B_nu_21)) * E_21 * kB * n_2 * n_ion

        heat_0 = B_01 * B_nu_10 * E_10 * kB * n_0 * n_ion
        heat_1 = B_02 * B_nu_20 * E_20 * kB * n_0 * n_ion
        heat_2 = B_12 * B_nu_21 * E_21 * kB * n_1 * n_ion

        ! Total cooling rate
        three_level = (cool_0 + cool_1 + cool_2) - (heat_0 + heat_1 + heat_2) 

    END FUNCTION three_level

    FUNCTION two_level(g_0, g_1, lam_10, A_10, &
                        z, T, n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                        idx_0)
        implicit none
        real(dp), intent(in)::g_0,g_1    ! degeneracy
        real(dp), intent(in)::A_10       ! s^-1
        real(dp), intent(in)::lam_10     ! microns
        real(dp), intent(in)::z                    ! redshift
        real(dp), intent(in)::n_ion,ne,nH,nHp,nHe,nHep,nHepp,nH2 ! Number densities
        real(dp), intent(in)::T                    ! temperature
        integer, intent(in)::idx_0       ! indicies of lines in the fine_struct_cooling array
        real(dp):: T_cmb,logT
        real(dp):: nu_10
        real(dp):: E_10
        real(dp):: B_01,B_10
        real(dp):: C_01,C_10
        real(dp):: B_nu_10
        real(dp):: q10_oH2
        real(dp):: q10_pH2
        real(dp):: q10_H
        real(dp):: q10_Hp
        real(dp):: q10_e
        real(dp):: q10_He
        real(dp):: q10_Hep
        real(dp):: q10_Hepp
        real(dp):: n_0, n_1
        real(dp):: cool_0
        real(dp):: heat_0
        real(dp):: t1, t2, nu_over_nl
        real(dp):: two_level

        ! Initialize the result
        two_level = 1.d-50

        T_cmb = 2.725d0 * (1.d0 + z) ! CMB temperature

        ! Get log T and enforce bounds
        logT = LOG10(T)

        logT = MAX(logT,table_temp(1))
        !if (logT .lt. MAX(LOG10(T_cmb),1.d0)) logT = MAX(LOG10(T_cmb),1.d0)
        if (logT .gt. 5.d0) logT = 5.d0

        nu_10 = c_cgs / (lam_10 * 1.d-4) ! Hz

        E_10 = (hplanck * c_cgs / (lam_10 * 1d-4)) / kB ! E/K (K)

        B_01 = A_10 * ((lam_10 * 1.d-4)**3.d0) * (g_1 / g_0) / (2.d0 * hplanck * c_cgs)

        B_10 = (g_0 / g_1) * B_01

        ! CMB black body spectrum
        B_nu_10 = (2.d0 * hplanck * (nu_10**3.d0) / (c_cgs**2.d0)) / (EXP(hplanck * nu_10 / (kB * T_cmb)) - 1.d0)

        ! All collison strengths are in cm^3 s^-1
        ! For collisions with o-H2
        !rate_electron,rate_HII,rate_HeII,rate_HeIII,rate_HI,rate_HeI,rate_oH2,rate_pH2
        q10_oH2 = interp_coll_1d(fs_cool_tab(idx_0,:,7), logT) 

        ! For collisions with p-H2
        q10_pH2 = interp_coll_1d(fs_cool_tab(idx_0,:,8), logT) 

        ! For collisions with H
        q10_H = interp_coll_1d(fs_cool_tab(idx_0,:,5), logT) 

        ! For collisions with H+
        q10_Hp = interp_coll_1d(fs_cool_tab(idx_0,:,2), logT) 

        ! For collisions with e
        q10_e = interp_coll_1d(fs_cool_tab(idx_0,:,1), logT) 

        ! For collisions with He
        q10_He = interp_coll_1d(fs_cool_tab(idx_0,:,6), logT) 

        ! For collisions with He+
        q10_Hep = interp_coll_1d(fs_cool_tab(idx_0,:,3), logT) 

        ! For collisions with He++
        q10_Hepp = interp_coll_1d(fs_cool_tab(idx_0,:,4), logT) 

        ! Net collision strengths
        C_10 = (q10_e * ne) + (q10_H * nH) + (q10_Hp * nHp) + (q10_oH2 * 0.75d0 * nH2) + (q10_pH2 * 0.25d0 * nH2) + (q10_He * nHe) + (q10_Hep * nHep) + (q10_Hepp * nHepp)

        C_01 = C_10 * (g_1/g_0) * EXP(-1.d0 * E_10 / T)

        ! Analytic solution to the 2 level system (modified version of paul goldsmith papers)
        ! Done this way to avoid numerical errors
        t1 = ( (B_01*B_nu_10) + C_01 )
        t2 = ( A_10 + (B_10*B_nu_10) + C_10 )

        nu_over_nl = t1 / t2
        n_0 = (nu_over_nl + 1.d0)**(-1.d0)
        n_1 = (nu_over_nl**-1.d0 + 1.d0)**(-1.d0)

        ! Cooling and heating rates
        cool_0 = (A_10 + (B_10 * B_nu_10)) * E_10 * kB * n_1 * n_ion

        heat_0 = B_01 * B_nu_10 * E_10 * kB * n_0 * n_ion

        ! Total cooling rate
        two_level = (cool_0 - heat_0)
   END FUNCTION two_level

   FUNCTION OI_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::OI_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=5.d0; g_1=3.d0; g_2=1.d0
        lam_10=63.1679d0; lam_20=44.0453d0; lam_21=145.495d0
        A_10=8.910d-05; A_20=1.340d-10; A_21=1.750d-05

        OI_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         8, 9, 10)
   END FUNCTION OI_fine_structure

   FUNCTION OIII_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::OIII_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=1.d0; g_1=3.d0; g_2=5.d0
        lam_10=88.3323d0; lam_20=32.6523d0; lam_21=51.8004d0
        A_10=2.597d-05; A_20=3.170d-11; A_21=9.760d-05

        OIII_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         13, 14, 15)
   END FUNCTION OIII_fine_structure

   FUNCTION CI_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::CI_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=1.d0; g_1=3.d0; g_2=5.d0
        lam_10=609.590d0; lam_20=230.352d0; lam_21=370.269d0
        A_10=7.930d-08; A_20=1.000d-30; A_21=2.650d-07

        CI_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         2, 3, 4)
   END FUNCTION CI_fine_structure

   FUNCTION CII_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::CII_fine_structure
        real(dp)::g_0,g_1         ! degeneracy
        real(dp)::A_10            ! s^-1
        real(dp)::lam_10          ! microns

        g_0=2.d0; g_1=4.d0
        lam_10=157.636d0
        A_10=2.290d-06
        
        CII_fine_structure = two_level( g_0, g_1, &
                                         lam_10, &
                                         A_10, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         1)
   END FUNCTION CII_fine_structure

   FUNCTION NII_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::NII_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=1.d0; g_1=3.d0; g_2=5.d0
        lam_10=205.244d0; lam_20=76.4318d0; lam_21=121.767d0
        A_10=2.080d-06; A_20=1.000d-30; A_21=7.460d-06

        NII_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         5, 6, 7)
   END FUNCTION NII_fine_structure

   FUNCTION SiI_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::SiI_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=1.d0; g_1=3.d0; g_2=5.d0
        lam_10=129.641d0; lam_20=44.7993d0; lam_21=68.4548d0
        A_10=8.250d-06; A_20=3.490d-10; A_21=4.210d-05

        SiI_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         22, 23, 24)
   END FUNCTION SiI_fine_structure

   FUNCTION SiII_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::SiII_fine_structure
        real(dp)::g_0,g_1         ! degeneracy
        real(dp)::A_10            ! s^-1
        real(dp)::lam_10          ! microns

        g_0=2.d0; g_1=4.d0
        lam_10=34.8046d0
        A_10=2.131d-04

        SiII_fine_structure = two_level( g_0, g_1, &
                                         lam_10, &
                                         A_10, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         12)
   END FUNCTION SiII_fine_structure

   FUNCTION NeII_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::NeII_fine_structure
        real(dp)::g_0,g_1         ! degeneracy
        real(dp)::A_10            ! s^-1
        real(dp)::lam_10          ! microns

        g_0=4.d0; g_1=2.d0
        lam_10=12.8101d0
        A_10=8.590d-03

        NeII_fine_structure = two_level( g_0, g_1, &
                                         lam_10, &
                                         A_10, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         11)
   END FUNCTION NeII_fine_structure

   FUNCTION FeI_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::FeI_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=9.d0; g_1=7.d0; g_2=5.d0
        lam_10=24.0358d0; lam_20=14.2005d0; lam_21=34.7038d0
        A_10=2.510d-03; A_20=1.000d-30; A_21=1.560d-03

        FeI_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         16, 17, 18)
   END FUNCTION FeI_fine_structure

   FUNCTION FeII_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::FeII_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=10.d0; g_1=8.d0; g_2=6.d0
        lam_10=25.9811d0; lam_20=14.9731d0; lam_21=35.3394d0
        A_10=2.050d-03; A_20=1.000d-30; A_21=1.560d-03

        FeII_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         19, 20, 21)
   END FUNCTION FeII_fine_structure

   FUNCTION SI_fine_structure(T,n_ion,nH,nHp,ne,nH2,nHe,nHep,nHepp,z)
        implicit none
        real(dp), intent(in)::T, n_ion, nH, nHp, ne, nH2, nHe, nHep, nHepp, z
        real(dp)::SI_fine_structure
        real(dp)::g_0,g_1,g_2          ! degeneracy
        real(dp)::A_10,A_20,A_21       ! s^-1
        real(dp)::lam_10,lam_20,lam_21 ! microns

        g_0=5.d0; g_1=3.d0; g_2=1.d0
        lam_10=25.2421d0; lam_20=17.4278d0; lam_21=56.2957d0
        A_10=1.400e-03; A_20=7.050e-08; A_21=3.020e-04

        SI_fine_structure = three_level( g_0, g_1, g_2, &
                                         lam_10, lam_20, lam_21, &
                                         A_10, A_20, A_21, &
                                         z, T, & 
                                         n_ion, ne, nH, nHp, nHe, nHep, nHepp, nH2, &
                                         25, 26, 27)
   END FUNCTION SI_fine_structure

END MODULE fscool
