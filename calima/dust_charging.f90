module dust_charging

    use dust_commons
    use dust_utils

    implicit none
    
    integer,parameter :: n_charge_threshold = 10

    ! charging_model = 'WDB06isrf': WDB06 charging under a uniform ISRF, per dust bin, from
    ! pyCALIMA export_dust_charging_isrf (dust_charging_isrf_DustBin_XX.dat in dust_tables_dir).
    ! Planes: <Z>, sigma_Z, ln Gamma_1, lambda_1, ln D_j (ISRF_ND), then ln alpha_i (ISRF_NA) of the
    ! grain-assisted recombination of the ions of dust_charging_rtgroups RI_ATOMIC (pyCALIMA
    ! RECOMB_IONS), read only for dust_ion_recombination with WDB06isrf
    integer, parameter :: ISRF_ND = 4, ISRF_NA = 10, ISRF_NQ_CHARGE = 4 + ISRF_ND, ISRF_NQ = ISRF_NQ_CHARGE + ISRF_NA
    integer, parameter, dimension(ISRF_ND) :: ISRF_DCHARGE = (/-1, 1, 2, 3/)  ! impactor charges of the D planes
    integer :: isrf_nT = 0, isrf_nX = 0, isrf_nqload = 0
    real(dp) :: isrf_lT0 = 0d0, isrf_lT1 = 0d0, isrf_dlT = 1d0, isrf_lX0 = 0d0, isrf_lX1 = 0d0, isrf_dlX = 1d0
    real(dp), allocatable :: isrf_q(:,:,:,:)         ! (isrf_nqload, nT, nX, ndust): ISRF_NQ_CHARGE or ISRF_NQ planes

    private
    public :: compute_mean_dust_charge, compute_dust_charge_sigma,&
              compute_dust_charge_dist, compute_Coulomb_focusing,&
              compute_Coulomb_focusing_ions,&
              init_dust_charging_isrf, isrf_lookup, ISRF_ND, ISRF_DCHARGE

    contains

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

        real(dp) :: Zsigma, Gam, Lam, Dj(ISRF_ND)

        ! charging_model = 'WDB06isrf' (WDB06rt is computed in compute_dust_precool and
        ! dust_charge_moments): the uniform-ISRF tables at G0
        call isrf_lookup(i_dust, G0, Tgas, ne, Zdust, Zsigma, Gam, Lam, Dj)
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

        real(dp) :: Z_avg, Gam, Lam, Dj(ISRF_ND)

        ! charging_model = 'WDB06isrf': the uniform-ISRF tables at G0
        call isrf_lookup(i_dust, G0, Tgas, ne, Z_avg, Zsigma, Gam, Lam, Dj)
    end subroutine compute_dust_charge_sigma

    subroutine init_dust_charging_isrf
        ! Reads the uniform-ISRF WDB06 charging tables of every dust bin (pyCALIMA
        ! export_dust_charging_isrf): <Z>, sigma_Z, ln Gamma_1, lambda_1 and ln D_j on equally
        ! spaced (log10 T, log10 n_e/G0), the same grid for every bin, each for the grain radius
        ! of its bin; see isrf_lookup. The ln alpha_i planes only for dust_ion_recombination with
        ! WDB06isrf (WDB06rt computes its own)
        use amr_commons, only: myid
        implicit none
        integer :: ii, k, i, nver, nT, nX, nQ, istat
        real(dp) :: a_tab
        integer, parameter :: u = 31
        character(len=40) :: dustlabel
        character(len=512) :: fname, line
        real(dp), allocatable :: lT(:), lX(:)

        do ii = 1, ndust
            write(dustlabel, '(A,I2.2)') 'DustBin_', ii
            if (len_trim(charging_table_label(ii)) > 0) dustlabel = charging_table_label(ii)
            fname = trim(dust_tables_dir)//'dust_charging_isrf_'//trim(dustlabel)//'.dat'
            open(u, file=trim(fname), status='old', action='read', iostat=istat)
            if (istat /= 0) then
                if (myid == 1) write(*,*) 'ERROR: cannot open ', trim(fname), ' (pyCALIMA export_dust_charging_isrf: ', &
                    'charging_model=''WDB06isrf'', and the PE heating of models other than WDB06rt)'
                call clean_stop
            end if
            do                                                  ! the header
                read(u, '(A)') line
                if (line(1:1) /= '#') exit
            end do
            read(line, *) nver, nT, nX, nQ, a_tab
            if (nver /= 2 .or. nQ /= ISRF_NQ) then
                if (myid == 1) write(*,*) 'ERROR: ', trim(fname), ' has format ', nver, ' with ', nQ, &
                    ' quantities; expected format 2 with ', ISRF_NQ, ' (re-export the tables)'
                call clean_stop
            end if
            if (abs(a_tab - dustbins_props(ii)%asize_cm) > 1d-6*a_tab) then
                if (myid == 1) write(*,*) 'ERROR: ', trim(fname), ' is for a = ', a_tab, &
                    ' cm, bin ', ii, ' has asize = ', dustbins_props(ii)%asize_cm, ' cm'
                call clean_stop
            end if
            allocate(lT(nT), lX(nX))
            read(u, *) lT
            read(u, *) lX
            if (ii == 1) then
                isrf_nT = nT
                isrf_nX = nX
                isrf_lT0 = lT(1)
                isrf_lT1 = lT(nT)
                isrf_dlT = lT(2) - lT(1)
                isrf_lX0 = lX(1)
                isrf_lX1 = lX(nX)
                isrf_dlX = lX(2) - lX(1)
                isrf_nqload = merge(ISRF_NQ, ISRF_NQ_CHARGE, dust_ion_recombination .and. trim(charging_model) == 'WDB06isrf')
                if (allocated(isrf_q)) deallocate(isrf_q)
                allocate(isrf_q(isrf_nqload, nT, nX, ndust))
            else if (nT /= isrf_nT .or. nX /= isrf_nX .or. lT(1) /= isrf_lT0 .or. lX(1) /= isrf_lX0 &
                     .or. lT(nT) /= isrf_lT1 .or. lX(nX) /= isrf_lX1) then
                if (myid == 1) write(*,*) 'ERROR: ', trim(fname), ' is not on the grid of DustBin_01'
                call clean_stop
            end if
            do k = 1, isrf_nqload
                do i = 1, nT
                    read(u, *) isrf_q(k, i, :, ii)
                end do
            end do
            close(u)
            deallocate(lT, lX)
        end do
        if (myid == 1) write(*,'(A,I4,A,I4,A)') ' CALIMA: uniform-ISRF charging tables read (', isrf_nT, ' x ', &
            isrf_nX, ' nodes per bin)'
    end subroutine init_dust_charging_isrf

    subroutine isrf_lookup(ii, G0, T, ne, Zmean, Zsigma, Gamma, Lambda, D, alpha)
        ! <Z>, sigma_Z, PE heating Gamma and net recombination cooling Lambda [erg/s per grain],
        ! the Coulomb focusing factors D(ISRF_DCHARGE) over the exact P(Z) and, if asked (the planes
        ! are read for dust_ion_recombination with WDB06isrf), the grain-assisted recombination
        ! coefficients alpha of the RI_ATOMIC ions [cm^3/s per grain], of dust bin ii in a uniform
        ! ISRF of G0 Habing at T [K] and n_e [cm^-3]: bilinear in (log10 T, log10 n_e/G0), clamped.
        ! Every rate is linear in G0 or n_e, so the tables at G0 = 1 give any G0 (pyCALIMA
        ! charging_isrf_tables.isrf_lookup)
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: G0, T, ne
        real(dp), intent(out) :: Zmean, Zsigma, Gamma, Lambda
        real(dp), dimension(ISRF_ND), intent(out) :: D
        real(dp), dimension(ISRF_NA), intent(out), optional :: alpha
        real(dp) :: fi, fj, u, v
        real(dp), dimension(ISRF_NQ) :: c
        integer :: i, j, n

        fi = (min(max(log10(T), isrf_lT0), isrf_lT1) - isrf_lT0) / isrf_dlT
        fj = (min(max(log10(max(ne, 1d-300) / max(G0, 1d-300)), isrf_lX0), isrf_lX1) - isrf_lX0) / isrf_dlX
        i = min(int(fi), isrf_nT - 2)
        j = min(int(fj), isrf_nX - 2)
        u = fi - dble(i)
        v = fj - dble(j)
        i = i + 1
        j = j + 1
        n = ISRF_NQ_CHARGE
        if (present(alpha)) n = isrf_nqload
        c(1:n) = (1d0 - u) * (1d0 - v) * isrf_q(1:n, i, j, ii) + u * (1d0 - v) * isrf_q(1:n, i + 1, j, ii) &
                 + (1d0 - u) * v * isrf_q(1:n, i, j + 1, ii) + u * v * isrf_q(1:n, i + 1, j + 1, ii)
        Zmean = c(1)
        Zsigma = c(2)
        Gamma = G0 * exp(c(3))
        Lambda = ne * sqrt(T) * c(4)
        D = exp(c(5:4 + ISRF_ND))
        if (present(alpha)) then
            alpha = 0d0
            if (n == ISRF_NQ) alpha = exp(c(ISRF_NQ_CHARGE + 1:ISRF_NQ))
        end if
    end subroutine isrf_lookup

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

    subroutine compute_Coulomb_focusing_ions(i_dust,Tgas,Zmean,Zsigma,jmax,D_coulomb)
        ! Coulomb focusing factors (Weingartner & Draine 1999, eqs. 6-7) of bin i_dust for the
        ! impactor charges -1..jmax, over a Gaussian charge distribution of mean Zmean and width
        ! Zsigma on the integer charges. For Zsigma < COUL_SIG_HI, the sum over every charge
        ! within COUL_KWIN max(Zsigma, 1/2) of the mean, the last unit of that width tapered to
        ! zero weight; for Zsigma > COUL_SIG_LO, the integral over the continuous Gaussian (the
        ! neutral charge: the bin [-1/2, 1/2]; the repulsive exponential at bin centres); blended
        ! in between. Continuous in Zmean, Zsigma and Tgas. Floored at COUL_FLOOR; once the
        ! repulsive part has fallen so far that every higher charge is at the floor, the rest are
        ! set to it.
        use constants, only: pi,e2instatC,kB
        implicit none
        integer, intent(in) :: i_dust, jmax
        real(dp), intent(in) :: Tgas,Zmean,Zsigma
        real(dp), dimension(-1:jmax), intent(inout) :: D_coulomb

        real(dp), parameter :: COUL_SIG_LO = 3.5d0, COUL_SIG_HI = 4.5d0, COUL_KWIN = 6d0, COUL_FLOOR = 1d-10
        real(dp), parameter :: SQRT2 = 1.4142135623730951d0, SQRT2PI = 2.5066282746310002d0
        real(dp), parameter :: PI_D = 3.141592653589793d0
        integer, parameter :: NW = 2*ceiling(COUL_KWIN*COUL_SIG_HI) + 3
        real(dp), dimension(1:NW) :: p
        real(dp) :: sig, sw, C1, w, x, s, q, q1, qj, c, B0, b1, att, neu, rep, bound, m, r, rho, arg
        real(dp) :: A0n, A1n, A0p, A1p, p0
        real(dp), dimension(2) :: tA, tR, PA, M1A, PN        ! 1: positive impactors, 2: negative (mirror)
        integer :: zlo, n, k, k0, j, kn, kp, z, im

        sig = max(Zsigma, 1d-3)
        C1 = e2instatC / (kB * Tgas * dustbins_props(i_dust)%asize_cm)
        x = min(max((sig - COUL_SIG_LO) / (COUL_SIG_HI - COUL_SIG_LO), 0d0), 1d0)
        w = x * x * (3d0 - 2d0 * x)

        ! 1. The discrete Gaussian and its sums over the negative, neutral and positive charges
        A0n = 0d0; A1n = 0d0; A0p = 0d0; A1p = 0d0; p0 = 0d0
        kn = 0
        kp = 1
        n = 0
        zlo = 0
        if (w < 1d0) then
            sw = max(sig, 0.5d0)
            zlo = floor(Zmean - COUL_KWIN * sw)
            n = ceiling(Zmean + COUL_KWIN * sw) - zlo + 1
            ! from the charge nearest the mean (weight 1) outward, p(Z+1) = p(Z) r with r(Z+1) =
            ! r(Z) exp(-1/sig^2): three exponentials instead of one per charge
            k0 = nint(Zmean) - zlo + 1
            rho = exp(-1d0 / (sig * sig))
            p(k0) = 1d0
            r = exp(-(2d0 * (dble(zlo + k0 - 1) - Zmean) + 1d0) / (2d0 * sig * sig))
            do k = k0 + 1, n
                p(k) = p(k-1) * r
                r = r * rho
            end do
            r = exp((2d0 * (dble(zlo + k0 - 1) - Zmean) - 1d0) / (2d0 * sig * sig))
            do k = k0 - 1, 1, -1
                p(k) = p(k+1) * r
                r = r * rho
            end do
            do k = 1, n
                x = abs(dble(zlo + k - 1) - Zmean) / sw
                if (x >= COUL_KWIN) then
                    p(k) = 0d0
                else if (x > COUL_KWIN - 1d0) then
                    p(k) = p(k) * cos(0.5d0 * PI_D * (x - COUL_KWIN + 1d0))**2
                end if
            end do
            p(1:n) = p(1:n) / sum(p(1:n))
            kp = n + 1
            do k = 1, n
                z = zlo + k - 1
                if (z < 0) then
                    A0n = A0n + p(k)
                    A1n = A1n + dble(z) * p(k)
                    kn = k
                else if (z == 0) then
                    p0 = p(k)
                else
                    A0p = A0p + p(k)
                    A1p = A1p + dble(z) * p(k)
                    kp = min(kp, k)
                end if
            end do
        end if

        ! 2. The continuous Gaussian: its mass and first moment below -1/2, and the neutral bin
        if (w > 0d0) then
            do im = 1, 2
                m = Zmean
                if (im == 2) m = -Zmean
                tA(im) = (-0.5d0 - m) / sig
                tR(im) = (0.5d0 - m) / sig
                ! beyond 8.5 sigma the normal integrals are 0 or 1 to double precision
                if (tA(im) < -8.5d0) then
                    PA(im) = 0d0
                    M1A(im) = 0d0
                else if (tA(im) > 8.5d0) then
                    PA(im) = 1d0
                    M1A(im) = m
                else
                    PA(im) = 0.5d0 * erfc(-tA(im) / SQRT2)
                    M1A(im) = m * PA(im) - sig * exp(-0.5d0 * tA(im)**2) / SQRT2PI
                end if
                if (tR(im) < -8.5d0) then
                    PN(im) = -PA(im)
                else if (tR(im) > 8.5d0) then
                    PN(im) = 1d0 - PA(im)
                else
                    PN(im) = 0.5d0 * erfc(-tR(im) / SQRT2) - PA(im)
                end if
                PN(im) = max(PN(im), 0d0)
            end do
        end if

        ! 3. Each impactor charge. The attractive and neutral parts grow with the charge: their value
        ! at jmax bounds them all, for the early exit
        b1 = sqrt(pi * C1 / 2d0)              ! B0 = 1 + |j| b1, the neutral grain
        q1 = exp(-C1)                         ! exp(-c) = q1**|j|
        c = dble(jmax) * C1
        B0 = 1d0 + dble(jmax) * b1
        bound = 0d0
        if (w < 1d0) bound = (1d0 - w) * (A0n - c * A1n + p0 * B0)
        if (w > 0d0) bound = bound + w * (PA(1) - c * M1A(1) + PN(1) * B0)
        D_coulomb(0) = 1d0
        qj = 1d0
        do j = -1, jmax
            if (j == 0) cycle
            im = 1
            if (j < 0) im = 2
            c = dble(abs(j)) * C1
            B0 = 1d0 + dble(abs(j)) * b1
            if (abs(j) == 1) then
                qj = q1
            else
                qj = qj * q1
            end if
            D_coulomb(j) = 0d0
            rep = 0d0
            if (w < 1d0) then
                q = qj
                s = 0d0
                if (j > 0) then
                    att = A0n - c * A1n
                    if (kp <= n .and. A0p >= 1d-17 * (att + p0 * B0)) then
                        do k = n, kp, -1
                            s = s * q + p(k)
                        end do
                        s = exp(max(-c * dble(zlo + kp - 1), -745d0)) * s
                    end if
                else
                    att = A0p + c * A1p
                    if (kn >= 1 .and. A0n >= 1d-17 * (att + p0 * B0)) then
                        do k = 1, kn
                            s = s * q + p(k)
                        end do
                        s = exp(max(c * dble(zlo + kn - 1), -745d0)) * s
                    end if
                end if
                D_coulomb(j) = (1d0 - w) * (att + p0 * B0 + s)
                rep = (1d0 - w) * s
            end if
            if (w > 0d0) then
                m = Zmean
                if (im == 2) m = -Zmean
                att = PA(im) - c * M1A(im)
                neu = PN(im) * B0
                x = tR(im) + c * sig
                if (x > 0d0) then
                    arg = -0.5d0 * c - 0.5d0 * tR(im)**2
                else
                    arg = -c * m + 0.5d0 * (c * sig)**2
                end if
                s = 0d0
                q = 1d0
                ! the repulsive part is at most exp(arg): left out below 1e-17 of the rest
                if (arg > -708d0 .and. .not. (att + neu >= 1d-3 .and. arg < -46d0)) then
                    if (x > 0d0) then
                        s = 0.5d0 * exp(arg) * erfc_scaled(x / SQRT2)
                    else if (x > -8.5d0) then
                        s = exp(arg) * 0.5d0 * erfc(x / SQRT2)
                    else
                        s = exp(arg)                        ! erfc = 2 to double precision
                    end if
                    ! (c/2)/sinh(c/2) = c exp(-c/2)/(1 - exp(-c))
                    if (c < 1d-3) then
                        q = 1d0 - c * c / 24d0
                    else
                        q = c * sqrt(qj) / (1d0 - qj)
                    end if
                end if
                D_coulomb(j) = D_coulomb(j) + w * (att + neu + q * s)
                rep = rep + w * q * s
            end if
            D_coulomb(j) = max(D_coulomb(j), COUL_FLOOR)
            if (j > 0 .and. bound + rep <= COUL_FLOOR) then
                D_coulomb(j:jmax) = COUL_FLOOR
                exit
            end if
        end do
    end subroutine compute_Coulomb_focusing_ions

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

end module dust_charging