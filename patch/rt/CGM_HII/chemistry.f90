
MODULE chemistry_module
    ! Module that houses the fex and jex parameters for solving
    ! the chemical network

    use amr_commons, ONLY: dp
    use hydro_parameters, only: n_oxygen_ions, n_nitrogen_ions, n_carbon_ions, n_magnesium_ions, n_silicon_ions, &
                                n_sulfur_ions, n_iron_ions, n_neon_ions
    use rt_parameters, only: rt, nGroups, include_collisional_ionization, include_photoionisation, &
                             UV_background_oxygen, UV_background_nitrogen, UV_background_carbon, & 
                             UV_background_magnesium, UV_background_silicon, UV_background_sulfur, &
                             UV_background_iron, UV_background_neon
    use coolrates_module
    DOUBLE PRECISION ne_chem, TK_chem
    real(dp), dimension(nGroups)::dNp_chem
    real(dp), dimension(nGroups, n_oxygen_ions)::signc_oxygen_chem
    real(dp), dimension(nGroups, n_nitrogen_ions)::signc_nitrogen_chem
    real(dp), dimension(nGroups, n_carbon_ions)::signc_carbon_chem
    real(dp), dimension(nGroups, n_magnesium_ions)::signc_magnesium_chem
    real(dp), dimension(nGroups, n_silicon_ions)::signc_silicon_chem
    real(dp), dimension(nGroups, n_sulfur_ions)::signc_sulfur_chem
    real(dp), dimension(nGroups, n_iron_ions)::signc_iron_chem
    real(dp), dimension(nGroups, n_neon_ions)::signc_neon_chem
    real(dp), dimension(n_oxygen_ions, n_oxygen_ions)::JAC_O
    real(dp), dimension(n_nitrogen_ions, n_nitrogen_ions)::JAC_N
    real(dp), dimension(n_carbon_ions, n_carbon_ions)::JAC_C
    real(dp), dimension(n_magnesium_ions, n_magnesium_ions)::JAC_MG
    real(dp), dimension(n_silicon_ions, n_silicon_ions)::JAC_SI
    real(dp), dimension(n_sulfur_ions, n_sulfur_ions)::JAC_S
    real(dp), dimension(n_iron_ions, n_iron_ions)::JAC_FE
    real(dp), dimension(n_neon_ions, n_neon_ions)::JAC_NE
    INTEGER :: chem_spec

    contains
        SUBROUTINE FEX(NEQ,T,Y,YDOT)
            IMPLICIT NONE
            INTEGER, INTENT (IN) :: NEQ
            DOUBLE PRECISION, INTENT (IN) :: T
            DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
            DOUBLE PRECISION, INTENT (OUT) :: YDOT(NEQ)
            DOUBLE PRECISION :: cr, de, photoRate
            INTEGER :: im

            if (chem_spec.eq.1) then ! Oxygen
                do im=1,n_oxygen_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_oxygen_ions) cr = cr + comp_Alpha_oxygen(TK_chem, 8-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_oxygen(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_oxygen_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_oxygen(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_oxygen_ions.and.include_collisional_ionization) de = de + comp_beta_oxygen(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_oxygen(TK_chem, 9-im)*ne_chem*Y(im) ! Recombination (the 9 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_oxygen_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_oxygen_ions) then
                            photoRate = photoRate + UV_background_oxygen(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.2) then ! Nitrogen
                do im=1,n_nitrogen_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_nitrogen_ions) cr = cr + comp_Alpha_nitrogen(TK_chem, 7-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_nitrogen(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_nitrogen_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_nitrogen(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_nitrogen_ions.and.include_collisional_ionization) de = de + comp_beta_nitrogen(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_nitrogen(TK_chem, 8-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_nitrogen_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_nitrogen_ions) then
                            photoRate = photoRate + UV_background_nitrogen(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.3) then ! Carbon
                do im=1,n_carbon_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_carbon_ions) cr = cr + comp_Alpha_carbon(TK_chem, 6-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_carbon(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_carbon_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_carbon(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_carbon_ions.and.include_collisional_ionization) de = de + comp_beta_carbon(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_carbon(TK_chem, 7-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_carbon_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_carbon_ions) then
                            photoRate = photoRate + UV_background_carbon(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.4) then ! Magnesium
                do im=1,n_magnesium_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_magnesium_ions) cr = cr + comp_Alpha_magnesium(TK_chem, 12-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_magnesium(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_magnesium_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_magnesium(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_magnesium_ions.and.include_collisional_ionization) de = de + comp_beta_magnesium(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_magnesium(TK_chem, 13-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_magnesium_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_magnesium_ions) then
                            photoRate = photoRate + UV_background_magnesium(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.5) then ! Silicon
                do im=1,n_silicon_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_silicon_ions) cr = cr + comp_Alpha_silicon(TK_chem, 14-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_silicon(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_silicon_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_silicon(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_silicon_ions.and.include_collisional_ionization) de = de + comp_beta_silicon(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_silicon(TK_chem, 15-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_silicon_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_silicon_ions) then
                            photoRate = photoRate + UV_background_silicon(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.6) then ! Sulfur
                do im=1,n_sulfur_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_sulfur_ions) cr = cr + comp_Alpha_sulfur(TK_chem, 16-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_sulfur(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_sulfur_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_sulfur(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_sulfur_ions.and.include_collisional_ionization) de = de + comp_beta_sulfur(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_sulfur(TK_chem, 17-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_sulfur_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_sulfur_ions) then
                            photoRate = photoRate + UV_background_sulfur(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.7) then ! Iron
                do im=1,n_iron_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_iron_ions) cr = cr + comp_Alpha_iron(TK_chem, 26-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_iron(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_iron_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_iron(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_iron_ions.and.include_collisional_ionization) de = de + comp_beta_iron(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_iron(TK_chem, 27-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_iron_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_iron_ions) then
                            photoRate = photoRate + UV_background_iron(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            else if (chem_spec.eq.8) then ! Neon
                do im=1,n_neon_ions
                    cr = 0.d0
                    de = 0.d0
                    photoRate = 0.d0

                    ! Creation
                    if (im.lt.n_neon_ions) cr = cr + comp_Alpha_neon(TK_chem, 10-im)*ne_chem*Y(im+1) ! Recombinations of the more excited ionization state
                    if (im.gt.1) then 
                        if (include_collisional_ionization) cr = cr + comp_beta_neon(TK_chem, im-1)*ne_chem*Y(im-1) ! Collisional ionization of the less excited state
                        if (rt.and.include_photoionisation) then 
                            cr = cr + SUM(signc_neon_chem(:, im-1)*dNp_chem)*Y(im-1) ! photoionization of the less excited state
                            cr = cr + UV_background_neon(im-1)*Y(im-1)
                        endif
                    endif

                    ! Destruction = collisional ionization + photoionization + recombination
                    if (im.lt.n_neon_ions.and.include_collisional_ionization) de = de + comp_beta_neon(TK_chem, im)*ne_chem*Y(im) ! Collisional ionization 
                    if (im.gt.1) de = de + comp_Alpha_neon(TK_chem, 11-im)*ne_chem*Y(im) ! Recombination (the 8 is correct)
                    
                    if (rt.and.include_photoionisation) then 
                        photoRate = SUM(signc_neon_chem(:, im)*dNp_chem)*Y(im)
                        if (im.lt.n_neon_ions) then
                            photoRate = photoRate + UV_background_neon(im)*Y(im) ! photoionization from uv background
                        endif
                    endif
                    !if (haardt_madau) photoRate = photoRate + UVrates(ixHeII, 1)*ss_factor
                    de = de + photoRate

                    YDOT(im) = cr - de
                enddo
            endif

            RETURN
        END SUBROUTINE FEX  

        SUBROUTINE JEX(NEQ,T,Y,ML,MU,PD,NRPD)
            IMPLICIT NONE
            INTEGER, INTENT (IN) :: NEQ, ML, MU, NRPD
            DOUBLE PRECISION, INTENT (IN) :: T
            DOUBLE PRECISION, INTENT (IN) :: Y(NEQ)
            DOUBLE PRECISION, INTENT (OUT) :: PD(NRPD,NEQ)
            INTEGER :: im

            PD = 0.d0
            if (chem_spec.eq.1) then ! Oxygen
                do im=1,n_oxygen_ions
                    PD(im,:) = JAC_O(im,:)
                enddo
            endif
            if (chem_spec.eq.2) then ! Nitrogen
                do im=1,n_nitrogen_ions
                    PD(im,:) = JAC_N(im,:)
                enddo
            endif
            if (chem_spec.eq.3) then ! Carbon
                do im=1,n_carbon_ions
                    PD(im,:) = JAC_C(im,:)
                enddo
            endif
            if (chem_spec.eq.4) then ! Magnesium
                do im=1,n_magnesium_ions
                    PD(im,:) = JAC_MG(im,:)
                enddo
            endif
            if (chem_spec.eq.5) then ! silicon
                do im=1,n_silicon_ions
                    PD(im,:) = JAC_SI(im,:)
                enddo
            endif
            if (chem_spec.eq.6) then ! sulfur
                do im=1,n_sulfur_ions
                    PD(im,:) = JAC_S(im,:)
                enddo
            endif
            if (chem_spec.eq.7) then ! iron
                do im=1,n_iron_ions
                    PD(im,:) = JAC_FE(im,:)
                enddo
            endif
            if (chem_spec.eq.8) then ! neon
                do im=1,n_neon_ions
                    PD(im,:) = JAC_NE(im,:)
                enddo
            endif
            RETURN
        END SUBROUTINE JEX

END MODULE chemistry_module
