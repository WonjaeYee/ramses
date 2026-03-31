module dust_charging

    use dust_commons
    use dust_utils
    
    contains

    subroutine compute_dust_charge_dist_Ibanez2019(G0,Tgas,ne,ispecie,agrain,Zdust,fcharge)
        ! ====== CHARGE DISTRIBUTION ======
        ! This is computed from the parametric fitting results of
        ! Ibanez-Mejia et al. (2019) - 
        ! (https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.1220I/abstract)
        ! What this assumes, and is shown in this work to be a pretty good
        ! assumption, is that charging is a very quick process, much faster
        ! than typical ISM/hydrodynamical scales
        use cooling_module, only: kB
        use constants, only: pi
        implicit none
        
        integer, intent(in) :: ispecie
        real(dp), intent(in) :: G0,Tgas,ne, agrain
        real(dp), dimension(:), allocatable, intent(inout) :: Zdust
        real(dp), dimension(:), allocatable, intent(inout) :: fcharge

        integer :: j,kk,isize
        integer :: Zmin,Zmax
        real(dp),dimension(1:7) :: Ibanez2019_sizes = (/3.5d-4,5d-4,1d-3,5d-3,1d-2,5d-2,1d-1/) ! in microns
        real(dp),dimension(1:2,1:7) :: alpha,k,b,hz,cplus,etaplus,d,cminus,etaminus
        real(dp) :: charPar,Z_avg,dist_sigma,Dtemp,Zg,Bfact,Zel

        ! Parameters from Table 1 in Ibanez-Mejia et al. (2019)
        alpha = transpose(reshape((/0.4699d0,0.4386d0,0.4994d0,0.6009d0,0.2900d0,0.3400d0,0.3500d0,&
                                    &0.3263d0,0.3141d0,0.3535d0,0.5115d0,0.3525d0,0.3643d0,0.3927d0/),[7,2]))
        k = transpose(reshape((/0.0085d0,0.0195d0,0.0199d0,0.0523d0,2.2310d0,5.8944d0,9.6536d0,&
                                &0.0149d0,0.0372d0,0.0494d0,0.0717d0,0.6591d0,2.6283d0,3.6493d0/),[7,2]))
        b = transpose(reshape((/-0.1162d0,-0.3084d0,-0.4959d0,-0.4092d0,-0.2061d0,0.1727d0,0.4183d0,&
                                &-0.1212d0,-0.3043d0,-0.4865d0,-0.4106d0,-0.1649d0,0.5217d0,0.8389d0/),[7,2]))
        hz = transpose(reshape((/48d0,95d0,78d0,218d0,1063d0,1034d0,1273d0,&
                                &57d0,86d0,73d0,107d0,384d0,345d0,372d0/),[7,2]))
        cplus = transpose(reshape((/0.3103d0,0.3699d0,0.6511d0,1.6536d0,2.5445d0,5.9455d0,8.7003d0,&
                                &0.4123d0,0.2734d0,0.4353d0,1.0758d0,1.6245d0,4.0732d0,5.9813d0/),[7,2]))
        etaplus = transpose(reshape((/0.2744d0,0.5654d0,0.9839d0,2.6688d0,4.3352d0,18.3186d0,36.1014d0,&
                                &0.2513d0,0.2925d0,0.7459d0,1.7832d0,2.8390d0,11.0200d0,20.6410d0/),[7,2]))
        d = transpose(reshape((/0.2551d0,0.4158d0,0.5275d0,0.6671d0,0.7010d0,0.8377d0,0.9094d0,&
                                &0.1891d0,0.3233d0,0.4451d0,0.5860d0,0.6346d0,0.6797d0,0.6961d0/),[7,2]))
        cminus = transpose(reshape((/0.3766d0,0.2890d0,-0.0213d0,-9.5138d0,-2.5341d3,-2.4189d3,-2.6009d3,&
                                &0.4845d0,0.3615d0,0.1053d0,-1.0379d3,-4.2075d2,-0.2418d0,-0.1885d0/),[7,2]))
        etaminus = transpose(reshape((/0.5241d0,1.6241d0,0.0977d0,35.3519d0,8.1962d3,4.9424d3,4.7029d3,&
                                &0.3532d0,0.6532d0,0.5803d0,7.7069d3,1.9840d3,0.5910d0,0.4237d0/),[7,2]))

        charPar = G0 * sqrt(Tgas) / ne

        ! TODO: This should be determine from the input sizes, because we may 
        ! have grain sizes not in the results of Ibanez-Mejia et al. (2019)
        isize = minloc(abs(Ibanez2019_sizes-agrain),1)
        
        ! Compute Gaussian charge distribution mean and sigma from fitting functions
        ! Eq. 17-19 in Ibanez-Mejia et al. (2019)
        Z_avg = k(ispecie,isize) * (1d0 - exp(-charPar / hz(ispecie,isize))) * (charPar**alpha(ispecie,isize)) + b(ispecie,isize)
        if (Z_avg>0) then
            dist_sigma = cplus(ispecie,isize) * (1d0 - exp(-Z_avg / etaplus(ispecie,isize))) + d(ispecie,isize)
        else
            dist_sigma = cminus(ispecie,isize) * (1d0 - exp(-abs(Z_avg) / etaminus(ispecie,isize))) + d(ispecie,isize)
        end if

        ! Compute approx. min and max of distribution by considering the points
        ! 3 sigma away from the mean (also, charge should be the nearest integer value)
        Zmin = nint(Z_avg - 3 * dist_sigma)
        Zmax = nint(Z_avg + 3 * dist_sigma)
        if (allocated(Zdust)) deallocate(Zdust)
        if (allocated(fcharge)) deallocate(fcharge)
        allocate(Zdust(1:(Zmax-Zmin+1)))
        allocate(fcharge(1:(Zmax-Zmin+1)))
        ! And now compute charge values and the Gaussian distribution
        do j=1,Zmax-Zmin+1
            Zdust(j) = dble(Zmin + j - 1)
            fcharge(j) = (1d0 / (dist_sigma * sq2pi)) * exp(-0.5d0*((Zdust(j) - Z_avg) / dist_sigma)**2)
        end do
        ! Renormalise distribution to make sure it adds to 1
        fcharge(:) = fcharge(:) / sum(fcharge(:))

    end subroutine compute_dust_charge_dist_Ibanez2019

    subroutine compute_mean_dust_charge_Ibanez2019(G0,Tgas,ne,ispecie,agrain,Zdust)
        ! ====== Mean Dust Charge ======
        ! This is computed from the parametric fitting results of
        ! Ibanez-Mejia et al. (2019) - 
        ! (https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.1220I/abstract)
        ! What this assumes, and is shown in this work to be a pretty good
        ! assumption, is that charging is a very quick process, much faster
        ! than typical ISM/hydrodynamical scales
        use cooling_module, only: kB
        use constants, only: pi
        implicit none
        
        integer, intent(in) :: ispecie
        real(dp), intent(in) :: G0,Tgas,ne, agrain
        real(dp), intent(inout) :: Zdust

        integer :: isize
        real(dp),dimension(1:7) :: Ibanez2019_sizes = (/3.5d-4,5d-4,1d-3,5d-3,1d-2,5d-2,1d-1/) ! in microns
        real(dp),dimension(1:2,1:7) :: alpha,k,b,hz,cplus,etaplus,d,cminus,etaminus
        real(dp) :: charPar

        ! Parameters from Table 1 in Ibanez-Mejia et al. (2019)
        alpha = transpose(reshape((/0.4699d0,0.4386d0,0.4994d0,0.6009d0,0.2900d0,0.3400d0,0.3500d0,&
                                    &0.3263d0,0.3141d0,0.3535d0,0.5115d0,0.3525d0,0.3643d0,0.3927d0/),[7,2]))
        k = transpose(reshape((/0.0085d0,0.0195d0,0.0199d0,0.0523d0,2.2310d0,5.8944d0,9.6536d0,&
                                &0.0149d0,0.0372d0,0.0494d0,0.0717d0,0.6591d0,2.6283d0,3.6493d0/),[7,2]))
        b = transpose(reshape((/-0.1162d0,-0.3084d0,-0.4959d0,-0.4092d0,-0.2061d0,0.1727d0,0.4183d0,&
                                &-0.1212d0,-0.3043d0,-0.4865d0,-0.4106d0,-0.1649d0,0.5217d0,0.8389d0/),[7,2]))
        hz = transpose(reshape((/48d0,95d0,78d0,218d0,1063d0,1034d0,1273d0,&
                                &57d0,86d0,73d0,107d0,384d0,345d0,372d0/),[7,2]))
        cplus = transpose(reshape((/0.3103d0,0.3699d0,0.6511d0,1.6536d0,2.5445d0,5.9455d0,8.7003d0,&
                                &0.4123d0,0.2734d0,0.4353d0,1.0758d0,1.6245d0,4.0732d0,5.9813d0/),[7,2]))
        etaplus = transpose(reshape((/0.2744d0,0.5654d0,0.9839d0,2.6688d0,4.3352d0,18.3186d0,36.1014d0,&
                                &0.2513d0,0.2925d0,0.7459d0,1.7832d0,2.8390d0,11.0200d0,20.6410d0/),[7,2]))
        d = transpose(reshape((/0.2551d0,0.4158d0,0.5275d0,0.6671d0,0.7010d0,0.8377d0,0.9094d0,&
                                &0.1891d0,0.3233d0,0.4451d0,0.5860d0,0.6346d0,0.6797d0,0.6961d0/),[7,2]))
        cminus = transpose(reshape((/0.3766d0,0.2890d0,-0.0213d0,-9.5138d0,-2.5341d3,-2.4189d3,-2.6009d3,&
                                &0.4845d0,0.3615d0,0.1053d0,-1.0379d3,-4.2075d2,-0.2418d0,-0.1885d0/),[7,2]))
        etaminus = transpose(reshape((/0.5241d0,1.6241d0,0.0977d0,35.3519d0,8.1962d3,4.9424d3,4.7029d3,&
                                &0.3532d0,0.6532d0,0.5803d0,7.7069d3,1.9840d3,0.5910d0,0.4237d0/),[7,2]))

        charPar = G0 * sqrt(Tgas) / ne

        ! TODO: This should be determine from the input sizes, because we may 
        ! have grain sizes not in the results of Ibanez-Mejia et al. (2019)
        isize = minloc(abs(Ibanez2019_sizes-agrain),1)
        
        ! Compute Gaussian charge distribution mean from fitting functions
        ! Eq. 17-19 in Ibanez-Mejia et al. (2019)
        Zdust = k(ispecie,isize) * (1d0 - exp(-charPar / hz(ispecie,isize))) * (charPar**alpha(ispecie,isize)) + b(ispecie,isize)
        Zdust = idnint(Zdust)  ! Should be the nearest integer value
    end subroutine compute_mean_dust_charge_Ibanez2019

    subroutine compute_mean_dust_charge(i_dust,G0,Tgas,ne,Zdust)
        ! ====== Mean Dust Charge ======
        ! This subroutine computes the mean dust charge following inteporlation
        ! of the equilibrium distribution properties presented in Rodríguez Montero
        ! et al. (2025), which in turn is based on the model from Weingartner &
        ! Draine (2001)
        implicit none
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: G0,Tgas,ne
        real(dp), intent(inout) :: Zdust

        real(dp) :: lgamma,lT

        ! 1. Compute charging parameter
        lgamma = log10(G0 * sqrt(Tgas) / ne)
        lT = log10(Tgas)

        ! 2. Interpolate the pre-computed per-grain table
        if (.not. dustbins_props(i_dust)%mean_charg_tab%initialised) then
            Zdust = 0d0
            return
        end if

        call interpolate2D(dustbins_props(i_dust)%mean_charg_tab%tab1d(1:dustbins_props(i_dust)%mean_charg_tab%npts(1),1), &
                           dustbins_props(i_dust)%mean_charg_tab%tab1d(1:dustbins_props(i_dust)%mean_charg_tab%npts(2),2), &
                           dustbins_props(i_dust)%mean_charg_tab%tab2d(1:dustbins_props(i_dust)%mean_charg_tab%npts(1),1:dustbins_props(i_dust)%mean_charg_tab%npts(2),1), &
                           dustbins_props(i_dust)%mean_charg_tab%npts(1), dustbins_props(i_dust)%mean_charg_tab%npts(2), lgamma, lT, Zdust)
    end subroutine compute_mean_dust_charge

    subroutine compute_dust_charge_sigma(i_dust,G0,Tgas,ne,Zsigma)
        ! ====== Dust Charge Sigma ======
        ! This subroutine computes the dust charge sigma following inteporlation
        ! of the equilibrium distribution properties presented in Rodríguez Montero
        ! et al. (2025), which in turn is based on the model from Weingartner &
        ! Draine (2001)
        implicit none
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: G0,Tgas,ne
        real(dp), intent(inout) :: Zsigma

        real(dp) :: lgamma,lT

        ! 1. Compute charging parameter
        lgamma = log10(G0 * sqrt(Tgas) / ne)
        lT = log10(Tgas)

        ! 2. Interpolate the pre-computed per-grain table
        if (.not. dustbins_props(i_dust)%sigma_charg_tab%initialised) then
            Zsigma = 1d0
            return
        end if

        call interpolate2D(dustbins_props(i_dust)%sigma_charg_tab%tab1d(1:dustbins_props(i_dust)%sigma_charg_tab%npts(1),1), &
                           dustbins_props(i_dust)%sigma_charg_tab%tab1d(1:dustbins_props(i_dust)%sigma_charg_tab%npts(2),2), &
                           dustbins_props(i_dust)%sigma_charg_tab%tab2d(1:dustbins_props(i_dust)%sigma_charg_tab%npts(1),1:dustbins_props(i_dust)%sigma_charg_tab%npts(2),1), &
                           dustbins_props(i_dust)%sigma_charg_tab%npts(1), dustbins_props(i_dust)%sigma_charg_tab%npts(2), lgamma, lT, Zsigma)

    end subroutine compute_dust_charge_sigma

    subroutine compute_dust_charge_dist(i_dust,G0,Tgas,ne,Zdust,fcharge)
        ! ====== CHARGE DISTRIBUTION ======
        ! This subroutine computes the dust charge distribution following inteporlation
        ! of the equilibrium distribution properties presented in Rodríguez Montero
        ! et al. (2025), which in turn is based on the model from Weingartner &
        ! Draine (2001)
        implicit none
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: G0,Tgas,ne
        real(dp), dimension(:), allocatable, intent(inout) :: Zdust
        real(dp), dimension(:), allocatable, intent(inout) :: fcharge

        integer :: j,kk,isize
        integer :: Zmin,Zmax
        real(dp) :: gamma,Z_avg,Zsigma

        ! 1. Compute charging parameter
        gamma = G0 * sqrt(Tgas) / ne

        ! 2. Get the interpolated mean and sigma of the distribution
        call compute_mean_dust_charge(i_dust,G0,Tgas,ne,Z_avg)
        call compute_dust_charge_sigma(i_dust,G0,Tgas,ne,Zsigma)

        ! 3. Compute approx. min and max of distribution by considering the points
        ! 3 sigma away from the mean (also, charge should be the nearest integer value)
        Zmin = nint(Z_avg - 3 * Zsigma)
        Zmax = nint(Z_avg + 3 * Zsigma)
        if (allocated(Zdust)) deallocate(Zdust)
        if (allocated(fcharge)) deallocate(fcharge)
        allocate(Zdust(1:(Zmax-Zmin+1)))
        allocate(fcharge(1:(Zmax-Zmin+1)))
        ! And now compute charge values and the Gaussian distribution
        do j=1,Zmax-Zmin+1
            Zdust(j) = dble(Zmin + j - 1)
            fcharge(j) = (1d0 / (Zsigma * sq2pi)) * exp(-0.5d0*((Zdust(j) - Z_avg) / Zsigma)**2)
        end do
        ! Renormalise distribution to make sure it adds to 1
        fcharge(:) = fcharge(:) / sum(fcharge(:))
    end subroutine compute_dust_charge_dist

    subroutine compute_Coulomb_focusing(Tgas,agrain,fcharge,Zdust,Zion,D_Coulomb)
        ! ====== Coulomb enhancement factor =====
        ! This is based on Eq. 6-7 in Weingartner & Draine (1999) which allows
        ! the computation of the Coulomb enhancement factor from the charge
        ! distribution (https://iopscience.iop.org/article/10.1086/307197)
        ! Remember that this equation is in CGS, so the grain size should be
        ! instead in cm, not in microns
        use cooling_module, only: kB
        use constants, only: pi
        implicit none
        
        real(dp), dimension(:), intent(in) :: fcharge
        real(dp), dimension(:), intent(in) :: Zdust
        real(dp), intent(in) :: Zion,agrain,Tgas
        real(dp), intent(inout) :: D_Coulomb

        integer :: j
        real(dp) :: Zg,Bfact

        D_Coulomb = 0d0
        if (Zion.ne.0d0) then
            ! Loop over the charge distribution, adding each contribution
            do j=1,size(Zdust,1)
                Zg = Zdust(j)
                if (Zg*Zion.gt.0) then
                    Bfact = exp(-Zg*Zion*e2instatC / (kB*Tgas*agrain))
                elseif (Zg*Zion.lt.0) then
                    Bfact = 1d0 - Zg*Zion*e2instatC / (kB*Tgas*agrain)
                elseif (Zg.eq.0) then
                    Bfact = 1d0 + sqrt(pi*Zion**2*e2instatC / (2d0*kB*Tgas*agrain))
                end if
                D_Coulomb = D_Coulomb + fcharge(j) * Bfact
            end do
            D_Coulomb = max(D_Coulomb,1d-10)
        else
            ! In the case of neutral atom, there is no Coulomb focusing
            D_Coulomb = 1d0
        end if

    end subroutine compute_Coulomb_focusing

end module dust_charging