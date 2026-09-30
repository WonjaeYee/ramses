! Grain charging, photoelectric heating and recombination cooling from the
! local RT photon groups (charging_model = 'WDB06rt').
!
! Tables: pyCALIMA  pycalima.models.dust_charge.export_dust_charging_rtgroups
! writes one file per dust bin, dust_charging_rtgroups_DustBin_XX.dat, with
! the photoemission and heating rates per unit photon flux of every group,
! K_g(U) [cm^2] and H_g(U) [cm^2 eV], tabulated in grain potential U and in the
! group mean energy (blackbody within-group shapes), and hot-gas tables in
! (T, U) for gas-electron secondaries and electron pass-through.
!
! Per cell and bin, with Up(Z) = sum_g c_red N_g K_g(U) + ion capture +
! gas-electron secondaries and Down(Z) = electron capture, the charge is found
! from detailed balance: root Z* of g(Z) = ln Up(Z) - ln Down(Z+1), then
! sigma^2 = -1/g'(Z*); exact P(Z) over integer charges when sigma < 3, a
! Gaussian (mean Z*+1/2, 3-point Gauss-Hermite averages) otherwise.
!
! This is a transcription of pycalima.models.dust_charge.rt_group_charging.
! GroupChargingTable.solve: same file, constants, interpolation and
! iteration, so the two agree to round-off. Change them together.
!
! Ion capture includes H+, He+ and He++ only. The UV background
! (rtz_UV_background_G0) is not added in this path.
module dust_charging_rtgroups
    use amr_parameters, only: dp
    use dust_commons

    implicit none

    private
    public :: init_dust_charging_rtgroups, rtgroups_set_group_energies, &
              rtgroups_solve_bin, rtgroups_debug_dump

    ! constants: identical to pyCALIMA shared_physics
    real(dp), parameter :: RTG_KB = 1.380649d-16, RTG_ME = 9.1093837015d-28
    real(dp), parameter :: RTG_ESTATC = 4.8032047d-10, RTG_EV2ERG = 1.602176634d-12
    real(dp), parameter :: RTG_PI = 3.14159265358979323846d0, RTG_TINY = 1d-300
    real(dp), parameter :: RTG_E2EVCM = RTG_ESTATC**2 / RTG_EV2ERG
    real(dp), parameter :: U_SCALE = 0.5d0, SIGMA_DISCRETE = 3d0
    integer, parameter :: MAX_DISCRETE_STATES = 121
    real(dp), parameter :: M_HION = 1.67262192d-24, M_HEION = 6.6446573d-24

    type RTGroupTable
        integer :: material = 0, nG = 0, nU = 0, nTs = 0, nTh = 0
        real(dp) :: a = 0d0, W = 0d0, Zmin = 0d0, Zmax = 0d0, x0 = 0d0, dx = 0d0
        real(dp), allocatable :: L0(:), L1(:), U(:), lTs(:), egy(:,:)        ! egy(nTs,nG)
        real(dp), allocatable :: K(:,:,:), H(:,:,:)                          ! (nU,nTs,nG)
        real(dp), allocatable :: lTh(:), delta(:,:), Ptr(:,:), fcool(:,:), Ebar(:,:)  ! (nU,nTh)
        real(dp), allocatable :: Kg(:,:), Hg(:,:), egy_used(:)               ! (nU,nG)
    end type RTGroupTable

    type(RTGroupTable), allocatable, save :: rtg(:)
    integer, save :: debug_count = 0, debug_unit = 0

    ! per-call cell state (set by rtgroups_solve_bin)
    real(dp) :: c_T, c_ne, c_ve, c_nion(3), c_zion(3), c_vion(3)
    real(dp), allocatable :: c_cN(:)
    integer :: c_bin

    contains

    ! ------------------------------------------------------------------ init
    subroutine init_dust_charging_rtgroups(groupL0, groupL1, nGroups)
        ! Read the tables of every dust bin and check them against the run's
        ! photon groups and grain sizes. Called from rt_init once groups exist.
        use amr_commons, only: myid
        implicit none
        integer, intent(in) :: nGroups
        real(dp), intent(in) :: groupL0(nGroups), groupL1(nGroups)
        integer :: ii, g, it, istat, ver, iu
        character(len=20) :: dustlabel
        character(len=256) :: fname
        character(len=512) :: line
        type(RTGroupTable) :: t

        if (allocated(rtg)) deallocate(rtg)
        allocate(rtg(1:ndust))
        do ii = 1, ndust
            write(dustlabel, '(A,I2.2)') 'DustBin_', ii
            fname = trim(dust_tables_dir)//'dust_charging_rtgroups_'//trim(dustlabel)//'.dat'
            iu = 27
            open(iu, file=trim(fname), status='old', action='read', iostat=istat)
            if (istat /= 0) then
                if (myid == 1) write(*,*) 'WDB06rt: cannot open ', trim(fname), &
                    ' (pyCALIMA export_dust_charging_rtgroups)'
                call clean_stop
            end if
            do
                read(iu, '(A)') line
                if (line(1:1) /= '#') exit
            end do
            read(line, *) ver, t%material, t%a, t%W, t%Zmin, t%Zmax
            read(iu, *) t%nG, t%nU, t%nTs, t%nTh
            if (t%nG /= nGroups) then
                if (myid == 1) write(*,*) 'WDB06rt: ', trim(fname), ' has ', t%nG, &
                    ' groups, the run has ', nGroups
                call clean_stop
            end if
            allocate(t%L0(t%nG), t%L1(t%nG), t%U(t%nU), t%lTs(t%nTs), t%egy(t%nTs,t%nG))
            allocate(t%K(t%nU,t%nTs,t%nG), t%H(t%nU,t%nTs,t%nG))
            allocate(t%lTh(t%nTh), t%delta(t%nU,t%nTh), t%Ptr(t%nU,t%nTh), &
                     t%fcool(t%nU,t%nTh), t%Ebar(t%nU,t%nTh))
            allocate(t%Kg(t%nU,t%nG), t%Hg(t%nU,t%nG), t%egy_used(t%nG))
            read(iu, *) t%L0
            read(iu, *) t%L1
            read(iu, *) t%U
            read(iu, *) t%lTs
            do g = 1, t%nG
                read(iu, *) t%egy(:,g)
            end do
            do g = 1, t%nG
                do it = 1, t%nTs
                    read(iu, *) t%K(:,it,g)
                end do
            end do
            do g = 1, t%nG
                do it = 1, t%nTs
                    read(iu, *) t%H(:,it,g)
                end do
            end do
            read(iu, *) t%lTh
            do it = 1, t%nTh
                read(iu, *) t%delta(:,it)
            end do
            do it = 1, t%nTh
                read(iu, *) t%Ptr(:,it)
            end do
            do it = 1, t%nTh
                read(iu, *) t%fcool(:,it)
            end do
            do it = 1, t%nTh
                read(iu, *) t%Ebar(:,it)
            end do
            close(iu)

            ! consistency with the run
            do g = 1, nGroups
                if (abs(t%L0(g) - groupL0(g)) > 1d-6*max(1d0,groupL0(g)) .or. &
                    abs(t%L1(g) - groupL1(g)) > 1d-6*max(1d0,groupL1(g))) then
                    if (myid == 1) write(*,*) 'WDB06rt: group ', g, ' edges in ', trim(fname), &
                        ' (', t%L0(g), t%L1(g), ') differ from groupL0/L1 (', groupL0(g), groupL1(g), ')'
                    call clean_stop
                end if
            end do
            if (abs(t%a - dustbins_props(ii)%asize_cm) > 1d-6*t%a) then
                if (myid == 1) write(*,*) 'WDB06rt: ', trim(fname), ' is for a = ', t%a, &
                    ' cm, bin ', ii, ' has asize = ', dustbins_props(ii)%asize_cm, ' cm'
                call clean_stop
            end if
            t%x0 = asinh(t%U(1)/U_SCALE)
            t%dx = (asinh(t%U(t%nU)/U_SCALE) - t%x0) / dble(t%nU - 1)
            t%egy_used = -1d0
            rtg(ii) = t
            deallocate(t%L0, t%L1, t%U, t%lTs, t%egy, t%K, t%H, t%lTh, t%delta, t%Ptr, &
                       t%fcool, t%Ebar, t%Kg, t%Hg, t%egy_used)
            if (myid == 1) write(*,'(A,A,A,ES12.5,A,F8.1,A,F8.1)') ' WDB06rt: read ', trim(dustlabel), &
                ', a = ', rtg(ii)%a, ' cm, Zmin = ', rtg(ii)%Zmin, ', Zmax = ', rtg(ii)%Zmax
        end do
        allocate(c_cN(1:nGroups))
    end subroutine init_dust_charging_rtgroups

    subroutine rtgroups_set_group_energies(ii, egy)
        ! Per-group K_g(U), H_g(U) at the run's mean group energies: linear in
        ! the table's group mean energy between the two bracketing shapes.
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: egy(:)
        integer :: g, k, j, n
        real(dp) :: eg, w

        n = rtg(ii)%nTs
        do g = 1, rtg(ii)%nG
            eg = min(max(egy(g), rtg(ii)%egy(1,g)), rtg(ii)%egy(n,g))
            k = n + 1
            do j = 1, n
                if (rtg(ii)%egy(j,g) >= eg) then
                    k = j
                    exit
                end if
            end do
            j = max(1, min(n - 1, k - 1))
            if (rtg(ii)%egy(j+1,g) == rtg(ii)%egy(j,g)) then
                w = 0d0
            else
                w = (eg - rtg(ii)%egy(j,g)) / (rtg(ii)%egy(j+1,g) - rtg(ii)%egy(j,g))
            end if
            rtg(ii)%Kg(:,g) = (1d0 - w)*rtg(ii)%K(:,j,g) + w*rtg(ii)%K(:,j+1,g)
            rtg(ii)%Hg(:,g) = (1d0 - w)*rtg(ii)%H(:,j,g) + w*rtg(ii)%H(:,j+1,g)
        end do
        rtg(ii)%egy_used = egy(1:rtg(ii)%nG)
    end subroutine rtgroups_set_group_energies

    ! ------------------------------------------------------------ interpolation
    subroutine u_index(U, i, w)
        implicit none
        real(dp), intent(in) :: U
        integer, intent(out) :: i
        real(dp), intent(out) :: w
        real(dp) :: f
        f = (asinh(U/U_SCALE) - rtg(c_bin)%x0) / rtg(c_bin)%dx
        f = min(max(f, 0d0), dble(rtg(c_bin)%nU - 1))
        i = min(int(f), rtg(c_bin)%nU - 2)
        w = f - dble(i)
        i = i + 1        ! 1-based
    end subroutine u_index

    real(dp) function hot(tab, i, w)
        implicit none
        real(dp), intent(in) :: tab(:,:)
        integer, intent(in) :: i
        real(dp), intent(in) :: w
        real(dp) :: y, dl, fT, v
        integer :: k, n
        n = rtg(c_bin)%nTh
        y = min(max(log10(c_T), rtg(c_bin)%lTh(1)), rtg(c_bin)%lTh(n))
        dl = (rtg(c_bin)%lTh(n) - rtg(c_bin)%lTh(1)) / dble(n - 1)
        fT = (y - rtg(c_bin)%lTh(1)) / dl
        k = min(int(fT), n - 2)
        v = fT - dble(k)
        k = k + 1
        hot = (1d0 - v)*((1d0 - w)*tab(i,k) + w*tab(i+1,k)) + v*((1d0 - w)*tab(i,k+1) + w*tab(i+1,k+1))
    end function hot

    ! -------------------------------------------------------------- physics
    real(dp) function jtilde(Z, q)
        ! Draine & Sutin (1987) J-tilde, theta_nu = nu/(1+nu^-1/2) (their eq. 3.6)
        implicit none
        real(dp), intent(in) :: Z, q
        real(dp) :: nu, tau, tn, th
        nu = Z / q
        tau = max(rtg(c_bin)%a*RTG_KB*c_T / max(q*q*RTG_ESTATC*RTG_ESTATC, RTG_TINY), RTG_TINY)
        if (nu == 0d0) then
            jtilde = 1d0 + sqrt(RTG_PI/(2d0*tau))
        else if (nu < 0d0) then
            tn = tau
            jtilde = (1d0 - nu/tn)*(1d0 + sqrt(2d0/max(tn - 2d0*nu, RTG_TINY)))
        else
            th = max(nu, RTG_TINY) / (1d0 + 1d0/sqrt(max(nu, RTG_TINY)))
            jtilde = (1d0 + 1d0/sqrt(4d0*tau + 3d0*max(nu, RTG_TINY)))**2 * exp(-th/tau)
        end if
    end function jtilde

    real(dp) function lambdatilde(Z)
        ! Draine & Sutin (1987) Lambda-tilde for electrons (q = -1)
        implicit none
        real(dp), intent(in) :: Z
        real(dp) :: nu, tau, th, nup
        nu = -Z
        tau = max(rtg(c_bin)%a*RTG_KB*c_T / max(RTG_ESTATC*RTG_ESTATC, RTG_TINY), RTG_TINY)
        if (nu == 0d0) then
            lambdatilde = 2d0 + 1.5d0*sqrt(RTG_PI/(2d0*tau))
        else if (nu < 0d0) then
            lambdatilde = (2d0 - nu/tau)*(1d0 + 1d0/sqrt(max(tau - nu, RTG_TINY)))
        else
            nup = max(nu, RTG_TINY)
            th = nup / (1d0 + 1d0/sqrt(nup))
            lambdatilde = (2d0 + nup/tau)*(1d0 + 1d0/sqrt(1.5d0/tau + 3d0*nup))*exp(-th/tau)
        end if
    end function lambdatilde

    real(dp) function sticking(Z)
        ! WD01 eqs. 27-30 with the continuous extension used by the solver
        implicit none
        real(dp), intent(in) :: Z
        real(dp) :: a, base, Nc
        a = rtg(c_bin)%a
        if (Z <= rtg(c_bin)%Zmin) then
            sticking = 0d0
            return
        end if
        base = 0.5d0*(1d0 - exp(-a/1d-7))
        if (Z > 0d0) then
            sticking = base
            return
        end if
        Nc = 468d0*(a/1d-7)**3
        sticking = base / (1d0 + exp(20d0 - Nc))
    end function sticking

    subroutine rates(Z, up, down, i, w)
        implicit none
        real(dp), intent(in) :: Z
        real(dp), intent(out) :: up, down, w
        integer, intent(out) :: i
        real(dp) :: a, photo, arrive_e, ion, sec
        integer :: g, k
        a = rtg(c_bin)%a
        call u_index(Z*RTG_E2EVCM/a, i, w)
        photo = 0d0
        do g = 1, rtg(c_bin)%nG
            photo = photo + c_cN(g)*((1d0 - w)*rtg(c_bin)%Kg(i,g) + w*rtg(c_bin)%Kg(i+1,g))
        end do
        arrive_e = RTG_PI*a*a*c_ne*c_ve*jtilde(Z, -1d0)
        ion = 0d0
        do k = 1, 3
            if (c_nion(k) > 0d0) ion = ion + RTG_PI*a*a*c_nion(k)*c_vion(k)*jtilde(Z, c_zion(k))
        end do
        sec = arrive_e*hot(rtg(c_bin)%delta, i, w)
        up = photo + ion + sec
        down = arrive_e*sticking(Z)*(1d0 - hot(rtg(c_bin)%Ptr, i, w))
    end subroutine rates

    real(dp) function gfun(Z)
        implicit none
        real(dp), intent(in) :: Z
        real(dp) :: up, down, up1, down1, w
        integer :: i
        call rates(Z, up, down, i, w)
        call rates(Z + 1d0, up1, down1, i, w)
        gfun = log(max(up, RTG_TINY)) - log(max(down1, RTG_TINY))
    end function gfun

    subroutine heating_cooling(Z, gam, lam)
        implicit none
        real(dp), intent(in) :: Z
        real(dp), intent(out) :: gam, lam
        real(dp) :: up, down, w, a, arrive, hsum
        integer :: i, g
        a = rtg(c_bin)%a
        call rates(Z, up, down, i, w)
        hsum = 0d0
        do g = 1, rtg(c_bin)%nG
            hsum = hsum + c_cN(g)*((1d0 - w)*rtg(c_bin)%Hg(i,g) + w*rtg(c_bin)%Hg(i+1,g))
        end do
        gam = RTG_EV2ERG*hsum
        arrive = RTG_PI*a*a*c_ne*c_ve
        lam = arrive*sticking(Z)*lambdatilde(Z)*RTG_KB*c_T*hot(rtg(c_bin)%fcool, i, w)
        lam = lam - arrive*jtilde(Z, -1d0)*hot(rtg(c_bin)%delta, i, w)*hot(rtg(c_bin)%Ebar, i, w)*RTG_EV2ERG
    end subroutine heating_cooling

    ! ---------------------------------------------------------------- solver
    subroutine rtgroups_solve_bin(ii, Np, egy, c_red, T, ne, n_Hp, n_Hep, n_Hepp, &
                                  Zmean, Zsigma, Gamma, Lambda, Lambda_auto, Zstar, discrete)
        ! Equilibrium charge and rates [erg/s per grain] of bin ii for one cell.
        ! Np: photon densities per group [cm^-3]; egy: group mean energies [eV];
        ! c_red: RT speed of light [cm/s].
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), egy(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp
        real(dp), intent(out) :: Zmean, Zsigma, Gamma, Lambda, Lambda_auto, Zstar
        logical, intent(out) :: discrete
        real(dp) :: Zlo, Zhi, glo, ghi, a_, b_, fa, fb, fz, h, zp, zm, gp, sigma, mu
        real(dp) :: lnf(MAX_DISCRETE_STATES), P(MAX_DISCRETE_STATES), Zint(MAX_DISCRETE_STATES)
        real(dp) :: nodes(3), wts(3), g_, l_, EA, arrive, s1, s2, mx
        integer :: side, it, half, Zc, z0, z1, n, k

        c_bin = ii
        if (any(rtg(ii)%egy_used /= egy(1:rtg(ii)%nG))) call rtgroups_set_group_energies(ii, egy)
        c_cN(1:rtg(ii)%nG) = c_red*Np(1:rtg(ii)%nG)
        c_T = T
        c_ne = ne
        c_ve = sqrt(8d0*RTG_KB*T/(RTG_PI*RTG_ME))
        c_nion = (/ n_Hp, n_Hep, n_Hepp /)
        c_zion = (/ 1d0, 1d0, 2d0 /)
        c_vion(1) = sqrt(8d0*RTG_KB*T/(RTG_PI*M_HION))
        c_vion(2) = sqrt(8d0*RTG_KB*T/(RTG_PI*M_HEION))
        c_vion(3) = c_vion(2)

        ! 1. root of g(Z) (Illinois)
        Zlo = rtg(ii)%Zmin
        Zhi = rtg(ii)%Zmax - 1d0
        glo = gfun(Zlo)
        ghi = gfun(Zhi)
        if (glo <= 0d0) then
            Zstar = Zlo
        else if (ghi >= 0d0) then
            Zstar = Zhi
        else
            a_ = Zlo; b_ = Zhi; fa = glo; fb = ghi; side = 0
            Zstar = b_
            do it = 1, 200
                Zstar = (a_*fb - b_*fa)/(fb - fa)
                if (abs(b_ - a_) < 1d-7*max(1d0, abs(Zstar))) exit
                fz = gfun(Zstar)
                if (fz == 0d0) exit
                if (fz*fb > 0d0) then
                    b_ = Zstar; fb = fz
                    if (side == -1) fa = 0.5d0*fa
                    side = -1
                else
                    a_ = Zstar; fa = fz
                    if (side == 1) fb = 0.5d0*fb
                    side = 1
                end if
            end do
        end if

        ! 2. width of the distribution
        h = 0.5d0
        zp = min(Zstar + h, Zhi)
        zm = max(Zstar - h, Zlo)
        gp = (gfun(zp) - gfun(zm)) / (zp - zm)
        sigma = 0d0
        if (gp < 0d0) sigma = sqrt(-1d0/gp)

        ! 3. discrete P(Z) or Gaussian
        Lambda_auto = 0d0
        if (sigma < SIGMA_DISCRETE) then
            discrete = .true.
            half = ceiling(6d0*max(sigma, 0.5d0)) + 2
            half = min(half, (MAX_DISCRETE_STATES - 1)/2)
            Zc = floor(Zstar + 0.5d0)
            z0 = int(max(rtg(ii)%Zmin, dble(Zc - half)))
            z1 = int(min(rtg(ii)%Zmax, dble(Zc + half)))
            n = z1 - z0 + 1
            do k = 1, n
                Zint(k) = dble(z0 + k - 1)
            end do
            lnf(1) = 0d0
            do k = 1, n - 1
                lnf(k+1) = lnf(k) + gfun(Zint(k))
            end do
            mx = maxval(lnf(1:n))
            P(1:n) = exp(lnf(1:n) - mx)
            P(1:n) = P(1:n)/sum(P(1:n))
            s1 = sum(P(1:n)*Zint(1:n))
            s2 = sum(P(1:n)*Zint(1:n)**2)
            Zmean = s1
            Zsigma = sqrt(max(s2 - s1*s1, 0d0))
            Gamma = 0d0
            Lambda = 0d0
            do k = 1, n
                call heating_cooling(Zint(k), g_, l_)
                Gamma = Gamma + P(k)*g_
                Lambda = Lambda + P(k)*l_
            end do
            if (Zint(1) == rtg(ii)%Zmin) then
                if (rtg(ii)%material == 1) then
                    EA = rtg(ii)%W + RTG_E2EVCM/rtg(ii)%a*((rtg(ii)%Zmin - 0.5d0) - 4d-8/(rtg(ii)%a + 7d-8))
                else
                    EA = rtg(ii)%W - 5d0 + RTG_E2EVCM/rtg(ii)%a*(rtg(ii)%Zmin - 0.5d0)
                end if
                arrive = RTG_PI*rtg(ii)%a**2*c_ne*c_ve
                Lambda_auto = P(1)*arrive*jtilde(rtg(ii)%Zmin, -1d0)*EA*RTG_EV2ERG
            end if
        else
            discrete = .false.
            mu = Zstar + 0.5d0
            nodes = (/ mu - sqrt(3d0)*sigma, mu, mu + sqrt(3d0)*sigma /)
            wts = (/ 1d0/6d0, 2d0/3d0, 1d0/6d0 /)
            Zmean = mu
            Zsigma = sigma
            Gamma = 0d0
            Lambda = 0d0
            do k = 1, 3
                call heating_cooling(nodes(k), g_, l_)
                Gamma = Gamma + wts(k)*g_
                Lambda = Lambda + wts(k)*l_
            end do
        end if
    end subroutine rtgroups_solve_bin

    subroutine rtgroups_debug_dump(ii, Np, egy, c_red, T, ne, n_Hp, n_Hep, n_Hepp, &
                                   Zmean, Zsigma, Gamma, Lambda, Lambda_auto, Zstar, discrete)
        ! Inputs and outputs of one call, for comparison with the pyCALIMA
        ! reference solver (diagnostics/dust_charge/compare_ramses_rtgroups.py).
        use amr_commons, only: myid
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), egy(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp
        real(dp), intent(in) :: Zmean, Zsigma, Gamma, Lambda, Lambda_auto, Zstar
        logical, intent(in) :: discrete
        character(len=32) :: fname
        integer :: nG
        if (debug_count >= dust_rtgroups_debug_max) return
        if (debug_unit == 0) then
            debug_unit = 931
            write(fname, '(A,I5.5,A)') 'dust_rtgroups_debug_', myid, '.dat'
            open(debug_unit, file=trim(fname), status='replace')
            write(debug_unit, '(A)') '# bin T ne n_Hp n_Hep n_Hepp m_Hion m_Heion c_red nG Np(1:nG) egy(1:nG) '// &
                'Zmean Zsigma Gamma Lambda Lambda_auto Zstar discrete'
        end if
        nG = size(Np)
        write(debug_unit, '(I3,8(1X,ES23.15E3),1X,I3,*(1X,ES23.15E3))', advance='no') ii, T, ne, n_Hp, n_Hep, &
            n_Hepp, M_HION, M_HEION, c_red, nG, Np(1:nG), egy(1:nG)
        write(debug_unit, '(6(1X,ES23.15E3),1X,L1)') Zmean, Zsigma, Gamma, Lambda, Lambda_auto, Zstar, discrete
        flush(debug_unit)
        debug_count = debug_count + 1
    end subroutine rtgroups_debug_dump

end module dust_charging_rtgroups
