module dust_charging

    use dust_commons
    use dust_utils

    implicit none
    
    integer,parameter :: n_charge_threshold = 10

    private
    public :: compute_mean_dust_charge, compute_dust_charge_sigma,&
              compute_dust_charge_dist, compute_Coulomb_focusing,&
              three_point_charge_mix, two_point_charge_mix

    contains

    subroutine compute_charge_moments_Ibanez2019(G0,Tgas,ne,ispecie,agrain,Z_avg,dist_sigma)
        ! Compute mean and sigma of the charge distribution from the fitting
        ! functions of Ibanez-Mejia et al. (2019).
        implicit none

        integer, intent(in) :: ispecie
        real(dp), intent(in) :: G0,Tgas,ne,agrain
        real(dp), intent(inout) :: Z_avg,dist_sigma

        integer :: isize
        real(dp),save,dimension(1:7) :: Ibanez2019_sizes = (/3.5d-4,5d-4,1d-3,5d-3,1d-2,5d-2,1d-1/) ! in microns
        real(dp),save,dimension(1:2,1:7) :: alpha,k,b,hz,cplus,etaplus,d,cminus,etaminus
        real(dp) :: charPar
        logical,save :: first_call = .true.

        if (first_call) then
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
            first_call = .false.
        end if

        charPar = G0 * sqrt(Tgas) / ne

        ! TODO: This should be determine from the input sizes, because we may
        ! have grain sizes not in the results of Ibanez-Mejia et al. (2019)
        isize = minloc(abs(Ibanez2019_sizes-agrain),1)

        ! Eq. 17-19 in Ibanez-Mejia et al. (2019)
        Z_avg = k(ispecie,isize) * (1d0 - exp(-charPar / hz(ispecie,isize))) * (charPar**alpha(ispecie,isize)) + b(ispecie,isize)
        if (Z_avg>0d0) then
            dist_sigma = cplus(ispecie,isize) * (1d0 - exp(-Z_avg / etaplus(ispecie,isize))) + d(ispecie,isize)
        else
            dist_sigma = cminus(ispecie,isize) * (1d0 - exp(-abs(Z_avg) / etaminus(ispecie,isize))) + d(ispecie,isize)
        end if
    end subroutine compute_charge_moments_Ibanez2019

    subroutine compute_dust_charge_dist_Ibanez2019(G0,Tgas,ne,ispecie,agrain,Zdust,fcharge)
        ! ====== CHARGE DISTRIBUTION ======
        ! This is computed from the parametric fitting results of
        ! Ibanez-Mejia et al. (2019) - 
        ! (https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.1220I/abstract)
        ! What this assumes, and is shown in this work to be a pretty good
        ! assumption, is that charging is a very quick process, much faster
        ! than typical ISM/hydrodynamical scales
        use constants, only: sq2pi
        implicit none
        
        integer, intent(in) :: ispecie
        real(dp), intent(in) :: G0,Tgas,ne, agrain
        real(dp), dimension(:), allocatable, intent(inout) :: Zdust
        real(dp), dimension(:), allocatable, intent(inout) :: fcharge

        integer :: j
        integer :: Zmin,Zmax
        real(dp) :: Z_avg,dist_sigma

        call compute_charge_moments_Ibanez2019(G0,Tgas,ne,ispecie,agrain,Z_avg,dist_sigma)

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
        implicit none
        
        integer, intent(in) :: ispecie
        real(dp), intent(in) :: G0,Tgas,ne, agrain
        real(dp), intent(inout) :: Zdust

        real(dp) :: dist_sigma

        call compute_charge_moments_Ibanez2019(G0,Tgas,ne,ispecie,agrain,Zdust,dist_sigma)
        Zdust = idnint(Zdust)  ! Should be the nearest integer value
    end subroutine compute_mean_dust_charge_Ibanez2019

    subroutine compute_mean_dust_charge(i_dust,G0,Tgas,ne,Zdust,idx_x,idx_y)
        ! ====== Mean Dust Charge ======
        ! This subroutine computes the mean dust charge following inteporlation
        ! of the equilibrium distribution properties presented in Rodríguez Montero
        ! et al. (2026), which in turn is based on the model from Weingartner &
        ! Draine (2001)
        implicit none
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: G0,Tgas,ne
        real(dp), intent(inout) :: Zdust
        integer, intent(inout), optional :: idx_x, idx_y

        real(dp) :: lgamma,lT
        real(dp) :: Zsigma
        integer :: ispecie

        if (trim(charging_model).eq.'Ibanez2019') then
            if (dustbins_props(i_dust)%separate_refractive_index) then
                ispecie = 1 ! Carbonaceous grains
            else
                ispecie = 2
            end if
            call compute_charge_moments_Ibanez2019(G0,Tgas,ne,ispecie,dustbins_props(i_dust)%asize,Zdust,Zsigma)
            return
        end if

        ! 1. Compute charging parameter
        lgamma = log10(max(G0,1d-6) * sqrt(Tgas) / max(ne,1d-20)) ! Avoid division by zero or very small numbers
        lT = log10(Tgas)

        ! 2. Interpolate the pre-computed per-grain table
        if (.not. dustbins_props(i_dust)%mean_charg_tab%initialised) then
            Zdust = 0d0
            return
        end if
        call dustbins_props(i_dust)%mean_charg_tab%interpolate(lgamma, lT, Zdust, idx_x, idx_y)
    end subroutine compute_mean_dust_charge

    subroutine compute_dust_charge_sigma(i_dust,G0,Tgas,ne,Zsigma,idx_x,idx_y)
        ! ====== Dust Charge Sigma ======
        ! This subroutine computes the dust charge sigma following inteporlation
        ! of the equilibrium distribution properties presented in Rodríguez Montero
        ! et al. (2026), which in turn is based on the model from Weingartner &
        ! Draine (2001)
        implicit none
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: G0,Tgas,ne
        real(dp), intent(inout) :: Zsigma
        integer, intent(inout), optional :: idx_x, idx_y

        real(dp) :: lgamma,lT
        real(dp) :: Z_avg
        integer :: ispecie

        if (trim(charging_model).eq.'Ibanez2019') then
            if (dustbins_props(i_dust)%separate_refractive_index) then
                ispecie = 1 ! Carbonaceous grains
            else
                ispecie = 2
            end if
            call compute_charge_moments_Ibanez2019(G0,Tgas,ne,ispecie,dustbins_props(i_dust)%asize,Z_avg,Zsigma)
            return
        end if

        ! 1. Compute charging parameter
        lgamma = log10(max(G0,1d-6) * sqrt(Tgas) / max(ne,1d-20)) ! Avoid division by zero or very small numbers
        lT = log10(Tgas)

        ! 2. Interpolate the pre-computed per-grain table
        if (.not. dustbins_props(i_dust)%sigma_charg_tab%initialised) then
            Zsigma = 1d0
            return
        end if

        call dustbins_props(i_dust)%sigma_charg_tab%interpolate(lgamma, lT, Zsigma, idx_x, idx_y)
    end subroutine compute_dust_charge_sigma

    subroutine compute_dust_charge_dist(i_dust,Z_avg,Zsigma,Zdust,fcharge,n_charge)
        ! ====== CHARGE DISTRIBUTION ======
        ! This subroutine computes the dust charge distribution following interpolation
        ! of the equilibrium distribution properties presented in Rodríguez Montero
        ! et al. (2026), which in turn is based on the model from Weingartner &
        ! Draine (2001)
        use constants, only: sq2pi
        implicit none
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: Z_avg,Zsigma
        real(dp), dimension(:), intent(inout) :: Zdust
        real(dp), dimension(:), intent(inout) :: fcharge
        integer, intent(out) :: n_charge

        integer :: j,kk,isize
        integer :: Zmin,Zmax
        real(dp) :: inv_Zsigma, inv_sqrt2pi_Zsigma, factor

        ! 1. Compute approx. min and max of distribution by considering the points
        ! 3 sigma away from the mean (also, charge should be the nearest integer value)
        Zmin = nint(Z_avg - 3 * Zsigma)
        Zmax = nint(Z_avg + 3 * Zsigma)
        n_charge = Zmax - Zmin + 1

        if (n_charge > size(Zdust)) then
            n_charge = size(Zdust)
            Zmin = nint(Z_avg - dble(n_charge)/2d0)
            Zmax = Zmin + n_charge - 1
        end if

        ! 2. And now compute charge values and the Gaussian distribution
        inv_Zsigma = 1d0 / Zsigma
        inv_sqrt2pi_Zsigma = 1d0 / (Zsigma * sq2pi)
        factor = -0.5d0 * (inv_Zsigma * inv_Zsigma)
        do j=1,n_charge
            Zdust(j) = dble(Zmin + j - 1)
            fcharge(j) = inv_sqrt2pi_Zsigma * exp(factor * (Zdust(j) - Z_avg)**2)
        end do
        ! 3. Renormalise distribution to make sure it adds to 1
        if (n_charge > 0) then
            fcharge(1:n_charge) = fcharge(1:n_charge) / sum(fcharge(1:n_charge))
        end if
    end subroutine compute_dust_charge_dist

    subroutine compute_Coulomb_focusing(i_dust,Tgas,Zmean,Zsigma,Zion,D_coulomb)
        ! ====== Coulomb enhancement factor =====
        ! This is based on Eq. 6-7 in Weingartner & Draine (1999) which allows
        ! the computation of the Coulomb enhancement factor from the charge
        ! distribution (https://iopscience.iop.org/article/10.1086/307197)
        use constants, only: pi,e2instatC,kB
        implicit none
        ! Input/Output variables
        integer, intent(in) :: i_dust
        real(dp), intent(in) :: Tgas,Zmean,Zsigma,Zion
        real(dp), intent(inout) :: D_coulomb

        ! Local variables
        integer :: Zmin,Zmax,n_charge
        real(dp) :: alpha, D_debug
        real(dp),dimension(1:n_charge_threshold) :: Zvals,fcharge

        ! 1. Determine the charge range within +/- 3 sigma
        Zmin = nint(Zmean - 3d0 * Zsigma)
        Zmax = nint(Zmean + 3d0 * Zsigma)
        n_charge = Zmax - Zmin + 1

        ! 2. Determine what method to use depending on the width and charge range
        ! if (.true.) then
        if (abs(Zmean/Zsigma) > 3d0 .and. n_charge > n_charge_threshold) then
            ! A. This uses a simple Taylor expansion around the mean charge
            alpha = (Zion * e2instatC) / (kB * Tgas * dustbins_props(i_dust)%asize_cm)
            if (Zmean * Zion > 0d0) then
                ! Repulsive Regime: Second order expansion of exp(-alpha * Z)
                ! Expands around mu: exp(-alpha*mu) * (1 + 0.5 * alpha^2 * sigma^2)
                D_coulomb = safe_exp(-alpha * Zmean) * (1d0 + 0.5d0 * alpha**2d0 * Zsigma**2d0)
            else
                ! Attractive Regime: Linear expansion of (1 - alpha * Z)
                ! Expansion around mu is exact (second derivative is 0)
                D_coulomb = 1d0 - alpha * Zmean
            end if
            ! call compute_dust_charge_dist(i_dust,Zmean,Zsigma,Zvals,fcharge,n_charge)
            ! call compute_Coulomb_focusing_dist(Tgas,dustbins_props(i_dust)%asize_cm,fcharge,Zvals,n_charge,Zion,D_debug)
            ! if (abs(D_coulomb-D_debug)/D_debug > 0.5d0 .and. D_debug >1d-10) then
            !     print*, 'WARNING: Error in taylor expansion'
            !     print*, 'i_dust: ',i_dust,' Zmean: ',Zmean,' Zsigma: ',Zsigma,' Zion: ',Zion,' G0: ',dust_helper%local_G0,' Tgas: ',Tgas,' ne: ',dust_helper%local_ne
            !     print*, 'D_coulomb: ',D_coulomb,' D_debug: ',D_debug, ' alpha: ',alpha
            !     print*,'Zion,e2instatC,kB,Tgas,dustbins_props(i_dust)%asize_cm: ',Zion,e2instatC,kB,Tgas,dustbins_props(i_dust)%asize_cm
            !     call clean_stop
            ! end if
        else
            ! B. In the case of a small distribution, it's more accurate
            ! to compute the Coulomb focusing factor directly from the
            ! charge distribution assumming a Gaussian distribution
            call compute_dust_charge_dist(i_dust,Zmean,Zsigma,Zvals,fcharge,n_charge)
            call compute_Coulomb_focusing_dist(Tgas,dustbins_props(i_dust)%asize_cm,fcharge,Zvals,n_charge,Zion,D_coulomb)
        end if
        
        ! 3. Make sure that the Coulomb factor does not become too small
        D_coulomb = max(D_coulomb, 1d-5)
    end subroutine compute_Coulomb_focusing

    subroutine compute_Coulomb_focusing_dist(Tgas,agrain,fcharge,Zdust,n_charge,Zion,D_Coulomb)
        ! ====== Coulomb enhancement factor =====
        ! This is based on Eq. 6-7 in Weingartner & Draine (1999) which allows
        ! the computation of the Coulomb enhancement factor from the charge
        ! distribution (https://iopscience.iop.org/article/10.1086/307197)
        ! Updated to safely support arbitrary array spacing (Delta Z != 1)
        use cooling_module, only: kB
        use constants, only: pi,e2instatC
        implicit none
        
        integer, intent(in) :: n_charge
        real(dp), dimension(n_charge), intent(in) :: fcharge
        real(dp), dimension(n_charge), intent(in) :: Zdust
        real(dp), intent(in) :: Zion,agrain,Tgas
        real(dp), intent(inout) :: D_Coulomb

        integer :: j, j_zero, j_start, j_end
        integer :: j_act_lo, j_act_hi
        real(dp) :: Zg, Bfact_zero, inv_kT_a, C
        real(dp) :: min_dist, current_dist
        real(dp), parameter :: f_thresh = 1.0d-20  ! Threshold to safely bypass denormalized numbers

        D_Coulomb = 0d0
        if (Zion.ne.0d0) then
            if (n_charge .gt. 0) then
                ! =========================================================
                ! ACTIVE WINDOW TRIM: Locate the living slice of the Gaussian
                ! =========================================================
                j_act_lo = 1
                do while (j_act_lo <= n_charge)
                    if (fcharge(j_act_lo) >= f_thresh) exit
                    j_act_lo = j_act_lo + 1
                end do

                j_act_hi = n_charge
                do while (j_act_hi >= j_act_lo)
                    if (fcharge(j_act_hi) >= f_thresh) exit
                    j_act_hi = j_act_hi - 1
                end do

                ! If the entire array was empty or below threshold, return default neutral focusing
                if (j_act_lo > j_act_hi) then
                    D_Coulomb = 1d0
                    return
                end if

                ! =========================================================
                ! ROBUST NEUTRAL LOCATOR: Find array index closest to zero charge
                ! =========================================================
                j_zero = j_act_lo
                min_dist = abs(Zdust(j_act_lo))
                do j = j_act_lo + 1, j_act_hi
                    current_dist = abs(Zdust(j))
                    if (current_dist < min_dist) then
                        min_dist = current_dist
                        j_zero = j
                    end if
                end do

                inv_kT_a = e2instatC / (kB*Tgas*agrain)
                C = Zion * inv_kT_a
                Bfact_zero = 1d0 + sqrt(pi * Zion**2 * inv_kT_a / 2d0)
                
                if (Zion .gt. 0d0) then
                    ! 1. Attractive Regime: Zg < 0 => Zg * Zion < 0
                    j_start = j_act_lo
                    j_end   = min(j_zero - 1, j_act_hi)
                    if (j_start .le. j_end) then
                        do j = j_start, j_end
                            Zg = Zdust(j)
                            D_Coulomb = D_Coulomb + fcharge(j) * (1d0 - Zg * C)
                        end do
                    end if
                    
                    ! 2. Neutral Grain: Zg == 0
                    if (j_zero .ge. j_act_lo .and. j_zero .le. j_act_hi) then
                        if (abs(Zdust(j_zero)) < 1d-5) then
                            D_Coulomb = D_Coulomb + fcharge(j_zero) * Bfact_zero
                        end if
                    end if
                    
                    ! 3. Repulsive Regime: Zg > 0 => Zg * Zion > 0
                    j_start = max(j_zero + 1, j_act_lo)
                    j_end   = j_act_hi
                    if (j_start .le. j_end) then
                        do j = j_start, j_end
                            D_Coulomb = D_Coulomb + fcharge(j) * exp(-Zdust(j) * C)
                        end do
                    end if
                    
                else
                    ! Zion < 0
                    ! 1. Repulsive Regime: Zg < 0 => Zg * Zion > 0
                    j_start = j_act_lo
                    j_end   = min(j_zero - 1, j_act_hi)
                    if (j_start .le. j_end) then
                        do j = j_start, j_end
                            D_Coulomb = D_Coulomb + fcharge(j) * exp(-Zdust(j) * C)
                        end do
                    end if
                    
                    ! 2. Neutral Grain: Zg == 0
                    if (j_zero .ge. j_act_lo .and. j_zero .le. j_act_hi) then
                        if (abs(Zdust(j_zero)) < 1d-5) then
                            D_Coulomb = D_Coulomb + fcharge(j_zero) * Bfact_zero
                        end if
                    end if
                    
                    ! 3. Attractive Regime: Zg > 0 => Zg * Zion < 0
                    j_start = max(j_zero + 1, j_act_lo)
                    j_end   = j_act_hi
                    if (j_start .le. j_end) then
                        do j = j_start, j_end
                            Zg = Zdust(j)
                            D_Coulomb = D_Coulomb + fcharge(j) * (1d0 - Zg * C)
                        end do
                    end if
                end if
            end if
        else
            ! Neutral impactor: No Coulomb focusing occurs
            D_Coulomb = 1d0
        end if
    end subroutine compute_Coulomb_focusing_dist

    subroutine two_point_charge_mix(mu, zmin, zlo, zhi, wlo, whi)
        implicit none
        real(dp), intent(in)  :: mu
        integer,  intent(in)  :: zmin
        integer,  intent(out) :: zlo, zhi
        real(dp), intent(out) :: wlo, whi

        zlo = ifloor(mu)
        zhi = zlo + 1

        whi = mu - real(zlo, dp)
        wlo = 1.0_dp - whi

        ! Enforce physical lower bound exactly as in the Python logic.
        if (zlo < zmin) then
        zlo = zmin
        zhi = zmin
        wlo = 1.0_dp
        whi = 0.0_dp
        else if (zhi < zmin) then
        zlo = zmin
        zhi = zmin
        wlo = 1.0_dp
        whi = 0.0_dp
        end if
    end subroutine two_point_charge_mix


    subroutine three_point_charge_mix(mu, sigma, zmin, z1, z2, z3, w1, w2, w3, used_two_point)
        implicit none
        real(dp), intent(in)  :: mu, sigma
        integer,  intent(in)  :: zmin
        integer,  intent(out) :: z1, z2, z3
        real(dp), intent(out) :: w1, w2, w3
        logical,  intent(out) :: used_two_point

        real(dp), parameter :: tol = 1.0d-10
        real(dp) :: sig, m2_target
        integer  :: half_span, zl, zh
        integer  :: i, j, k
        real(dp) :: a, b, c, a2, b2, c2
        real(dp) :: da, db, da2, db2, rhs1, rhs2, det
        real(dp) :: ww1, ww2, ww3
        real(dp) :: score_span, score_mid
        real(dp) :: best_span, best_mid
        logical  :: found_nonneg

        ! Safe sigma handling.
        sig = sigma
        if (sig /= sig) sig = 0.0_dp
        if (sig < 0.0_dp) sig = 0.0_dp

        m2_target = mu*mu + sig*sig

        ! Search window around mu.
        half_span = max(3, iceil(4.0_dp*sig + 2.0_dp))

        zl = max(zmin, ifloor(mu) - half_span)
        zh = iceil(mu) + half_span
        if (zh - zl < 2) zh = zl + 2

        found_nonneg = .false.
        best_span = huge(1.0_dp)
        best_mid  = huge(1.0_dp)

        ! Initialize outputs.
        z1 = zl
        z2 = zl + 1
        z3 = zl + 2
        w1 = 0.0_dp
        w2 = 0.0_dp
        w3 = 0.0_dp

        do i = zl, zh - 2
        do j = i + 1, zh - 1
            do k = j + 1, zh
            a = real(i, dp)
            b = real(j, dp)
            c = real(k, dp)
            a2 = a*a
            b2 = b*b
            c2 = c*c

            ! Solve:
            ! w1+w2+w3=1
            ! w1*a + w2*b + w3*c = mu
            ! w1*a2+w2*b2+w3*c2 = m2_target
            !
            ! Reduced 2x2 system for w1,w2 with w3=1-w1-w2.
            da   = a - c
            db   = b - c
            da2  = a2 - c2
            db2  = b2 - c2
            rhs1 = mu - c
            rhs2 = m2_target - c2

            det = da*db2 - db*da2
            if (abs(det) <= 1.0d-18) cycle

            ww1 = ( rhs1*db2 - rhs2*db ) / det
            ww2 = ( da*rhs2  - da2*rhs1 ) / det
            ww3 = 1.0_dp - ww1 - ww2

            if ((ww1 >= -tol) .and. (ww2 >= -tol) .and. (ww3 >= -tol)) then
                ! Prefer compact support and center near mu.
                score_span = real(k - i, dp)
                score_mid  = abs(real(j, dp) - mu)

                if ((.not. found_nonneg) .or. &
                    (score_span < best_span) .or. &
                    (abs(score_span - best_span) <= 1.0d-12 .and. score_mid < best_mid)) then
                found_nonneg = .true.
                best_span = score_span
                best_mid  = score_mid

                z1 = i
                z2 = j
                z3 = k
                w1 = ww1
                w2 = ww2
                w3 = ww3
                end if
            end if

            end do
        end do
        end do

        if (found_nonneg) then
        ! Clip tiny negatives and renormalize.
        if (w1 < 0.0_dp) w1 = 0.0_dp
        if (w2 < 0.0_dp) w2 = 0.0_dp
        if (w3 < 0.0_dp) w3 = 0.0_dp
        call renormalize3(w1, w2, w3)
        used_two_point = .false.
        else
        ! Fallback: two-point method.
        call two_point_charge_mix(mu, zmin, z1, z2, w1, w2)
        z3 = z2
        w3 = 0.0_dp
        used_two_point = .true.
        end if
    end subroutine three_point_charge_mix


    subroutine renormalize3(w1, w2, w3)
        implicit none
        real(dp), intent(inout) :: w1, w2, w3
        real(dp) :: s

        s = w1 + w2 + w3
        if (s > 0.0_dp) then
        w1 = w1 / s
        w2 = w2 / s
        w3 = w3 / s
        else
        ! Defensive fallback (should not happen in normal flow).
        w1 = 1.0_dp
        w2 = 0.0_dp
        w3 = 0.0_dp
        end if
    end subroutine renormalize3


    integer function ifloor(x)
        implicit none
        real(dp), intent(in) :: x
        integer :: t

        t = int(x)
        if (real(t, dp) > x) t = t - 1
        ifloor = t
    end function ifloor


    integer function iceil(x)
        implicit none
        real(dp), intent(in) :: x
        integer :: t

        t = int(x)
        if (real(t, dp) < x) t = t + 1
        iceil = t
    end function iceil

end module dust_charging