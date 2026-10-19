module dust_charging

    use dust_commons
    use dust_utils
    use iso_c_binding, only: c_double

    implicit none
    

    ! charging_model = 'WDB06isrf': WDB06 charging under a uniform ISRF, per dust bin, from
    ! pyCALIMA export_dust_charging_isrf (dust_charging_isrf_DustBin_XX.dat in dust_tables_dir).
    ! Planes: <Z>, sigma_Z, ln Gamma_1, lambda_1, ln D_j (ISRF_ND), then ln alpha_i (ISRF_NA) of the
    ! grain-assisted recombination of the ions RI_ATOMIC (WDB06rt part below; pyCALIMA
    ! RECOMB_IONS), read only for dust_ion_recombination with WDB06isrf
    integer, parameter :: ISRF_ND = 4, ISRF_NA = 10, ISRF_NQ_CHARGE = 4 + ISRF_ND, ISRF_NQ = ISRF_NQ_CHARGE + ISRF_NA
    integer, parameter, dimension(ISRF_ND) :: ISRF_DCHARGE = (/-1, 1, 2, 3/)  ! impactor charges of the D planes
    integer :: isrf_nT = 0, isrf_nX = 0, isrf_nqload = 0
    real(dp) :: isrf_lT0 = 0d0, isrf_lT1 = 0d0, isrf_dlT = 1d0, isrf_lX0 = 0d0, isrf_lX1 = 0d0, isrf_dlX = 1d0
    real(dp), allocatable :: isrf_q(:,:,:,:)         ! (isrf_nqload, nT, nX, ndust): ISRF_NQ_CHARGE or ISRF_NQ planes

    ! ======================================================== charging_model = 'WDB06rt'
    ! Grain charging, photoelectric heating, recombination cooling and grain-assisted
    ! ion recombination from the local RT photon groups (charging_model = 'WDB06rt').
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
    ! from detailed balance, g(Z) = ln Up(Z) - ln Down(Z+1):
    !  - root Z* of g and sigma^2 = -1/g'(Z*); if sigma < 3, P(Z) exactly on the
    !    integers around Z* (window trimmed to P > DISCRETE_TAIL = 1e-5 of its maximum), else
    !    a Gaussian of mean Z*+1/2, whose heating and cooling are those at Z*+1/2;
    !  - d/d ln T of <Z>, sigma, Gamma, Lambda (dlnT > 0), which serve the cooling
    !    solver's T(1 +- 1e-5) calls: analytic on a discrete window, from the pair
    !    g(Z* +- 1/2) re-evaluated at T(1 + dlnT) for a Gaussian;
    !  - grain-assisted recombination alpha_i of the RTG_NRI singly charged ions
    !    (case A thresholds): the exact sum on a discrete window; for a Gaussian, the
    !    segmented Chebyshev reconstruction of ln P (recomb_fit), for the bins dust_interface
    !    finds important; else the Gaussian estimate (recomb_gauss).
    ! Warm start (RTGState, kept per cell and bin by dust_interface): a discrete
    ! window is regrown in the new cell with no root search; a Gaussian's root costs
    ! g at Z*_prev +- 1/2.
    !
    ! This is a transcription of pycalima.models.dust_charge.rt_group_charging.
    ! GroupChargingTable.solve with kernels='pchip': same file, constants,
    ! interpolation and iteration, so the two agree to round-off. Change them
    ! together. Implementation differences that do not change the results beyond
    ! round-off: P on a window by running products of Up/Down (Python: exp of the
    ! cumulative sum of g), the U index of every integer charge cached per bin, and
    ! rates, heating and cooling of each integer of a window evaluated once.
    !
    ! Ion capture includes H+, He+ and He++ only. The UV background
    ! (rtz_UV_background_G0, in Habing units) adds G0_bg * K_bg(U) to the
    ! photoemission (Mathis ISRF kernels in the table). Cells with no local RT
    ! photons use a (log T, log ne) table of the full solve at the run's G0 (ions:
    ! H+ with n = ne), built once at first use and again only if G0 changes.
    ! recombination of a wide P(Z) in rtgroups_solve_bin: none, the Gaussian estimate, the (reused) fit
    integer, parameter :: RTG_ALPHA_NONE = 0, RTG_ALPHA_GAUSS = 1
    ! dust_interface fits only the wide bins with more than this share of an ion's grain recombination
    real(dp), parameter :: RTG_RECOMB_IMPORTANT = 0.1d0
    ! rtgroups_predict guards: identical to pyCALIMA rt_group_charging PREDICT_WIDE_ZERO, _WIDE_DZ
    real(dp), parameter :: PREDICT_WIDE_ZERO = 0.5d0, PREDICT_WIDE_DZ = 1d0
    ! ... and a wide P(Z) is linearised (else corrected) within PREDICT_LIN_TOL, _DZ and _ZERO
    real(dp), parameter :: PREDICT_LIN_TOL = 0.1d0, PREDICT_LIN_DZ = 0.25d0, PREDICT_LIN_ZERO = 1d0

    ! constants: identical to pyCALIMA shared_physics / rt_group_charging
    real(dp), parameter :: RTG_KB = 1.380649d-16, RTG_ME = 9.1093837015d-28
    real(dp), parameter :: RTG_ESTATC = 4.8032047d-10, RTG_EV2ERG = 1.602176634d-12
    real(dp), parameter :: RTG_PI = 3.14159265358979323846d0, RTG_TINY = 1d-300
    real(dp), parameter :: RTG_E2EVCM = RTG_ESTATC**2 / RTG_EV2ERG
    real(dp), parameter :: RTG_MU = 1.66053906660d-24
    real(dp), parameter :: GRAPHITE_WF = 4.4d0, SILICATE_WF = 8.0d0, SILICATE_GAP = 5.0d0
    real(dp), parameter :: U_SCALE = 0.5d0, SIGMA_DISCRETE = 3d0, DISCRETE_TAIL = 1d-5
    integer, parameter :: SIGMA_CHORD_ITER = 3       ! at most this many chords for the width of a wide P(Z)
    integer, parameter :: MAX_DISCRETE_STATES = 121, NBUF = 3*MAX_DISCRETE_STATES + 4
    real(dp), parameter :: M_HION = 1.67262192d-24, M_HEION = 6.6446573d-24
    ! uniform-G0 table grid: identical to pyCALIMA rt_group_charging UNIF_LT / UNIF_LNE
    integer, parameter :: NUNI_T = 81, NUNI_NE = 97
    real(dp), parameter :: UNI_LT0 = 1d0, UNI_LT1 = 9d0, UNI_LNE0 = -8d0, UNI_LNE1 = 8d0
    ! grain-assisted recombination: pyCALIMA RECOMB_IONS and RECOMB_* (rt_group_charging)
    integer, parameter :: RTG_NRI = 10
    integer, parameter :: RTG_MAXG = 32          ! rtgroups_predict: at most this many groups (else no anchor)
    real(dp), parameter :: RI_MASS(RTG_NRI) = (/1.00794d0, 4.002602d0, 12.0107d0, 14.0067d0, 15.9994d0, &
                                                20.1797d0, 24.305d0, 28.0855d0, 32.065d0, 55.845d0/)
    real(dp), parameter :: RI_IP(RTG_NRI) = (/13.598d0, 24.587d0, 11.260d0, 14.534d0, 13.618d0, &
                                              21.565d0, 7.646d0, 8.152d0, 10.360d0, 7.902d0/)
    integer, parameter :: RECOMB_FIT_N = 5, RECOMB_FIT_NA = 5, RECOMB_KINK = 2
    real(dp), parameter :: RECOMB_PER_NODE = 5d0, RECOMB_DOWN = 8d0, RECOMB_UP = 5d0, RECOMB_REACH = 8d0
    real(dp), parameter :: RECOMB_CUT = 25d0, RECOMB_EPS = 0.08d0
    real(dp), parameter :: RECOMB_FLOOR = log(1d-8)
    real(dp), parameter :: MATHIS_BREAKS(4) = (/5.04d0, 9.26d0, 11.2d0, 13.6d0/)
    real(dp), parameter :: NONE = -huge(1d0)
    integer, parameter :: MAXP = 4096       ! fit pieces
    integer, parameter :: MAXB = 64         ! break charges

    type RTGroupTable
        integer :: material = 0, nG = 0, nU = 0, nTs = 0, nTh = 0, nZ = 0
        real(dp) :: a = 0d0, W = 0d0, Zmin = 0d0, Zmax = 0d0, x0 = 0d0, dx = 0d0
        real(dp), allocatable :: L0(:), L1(:), U(:), lTs(:), egy(:,:)        ! egy(nTs,nG)
        real(dp), allocatable :: K(:,:,:), H(:,:,:)                          ! (nU,nTs,nG)
        real(dp), allocatable :: lTh(:), delta(:,:), Ptr(:,:), fcool(:,:), Ebar(:,:)  ! (nU,nTh)
        real(dp), allocatable :: KgT(:,:), HgT(:,:), mKT(:,:), mHT(:,:)      ! (nG,nU) kernels, PCHIP slopes
        real(dp), allocatable :: egy_used(:)
        real(dp), allocatable :: Kbg(:), Hbg(:), mKb(:), mHb(:)              ! (nU) Mathis ISRF, G0 = 1
        real(dp) :: G0_uni = -1d0
        real(dp), allocatable :: uni(:,:,:)                                  ! (4,NUNI_T,NUNI_NE)
        real(dp) :: stick_pos = 0d0, stick_neg = 0d0                         ! WD01 sticking, Z > 0 / Z <= 0
        real(dp) :: Zth(RTG_NRI) = 0d0                                       ! case-A recombination thresholds
        real(dp), allocatable :: edge0(:), edge1(:)                          ! (nG) charges where IP(a,Z) crosses L0, L1
        real(dp) :: edgebg(4) = NONE                                         ! ... and the Mathis breaks
        integer, allocatable :: iuz(:)                                       ! (nZ) U index of each integer charge
        real(dp), allocatable :: wuz(:)
        integer :: jpos = 1                                                  ! first U node > 0
    end type RTGroupTable

    type RTGResult
        real(dp) :: Zmean = 0d0, Zsigma = 0d0, Gamma = 0d0, Lambda = 0d0, Lambda_auto = 0d0, Zstar = 0d0
        real(dp) :: dZmean = 0d0, dZsigma = 0d0, dGamma = 0d0, dLambda = 0d0   ! d/d ln T (dLambda: of Lambda + auto)
        real(dp) :: alpha(RTG_NRI) = 0d0                                         ! cm^3 s^-1 per grain
        logical :: discrete = .false.
        ! sens of rtgroups_solve_bin: d/d ln S (every photon density scaled) and d/d ln N (n_e and the
        ! ion densities scaled), for rtgroups_predict (dLambda: of Lambda + auto)
        real(dp) :: dZmean_S = 0d0, dZsigma_S = 0d0, dGamma_S = 0d0, dLambda_S = 0d0
        real(dp) :: dZmean_N = 0d0, dZsigma_N = 0d0, dGamma_N = 0d0, dLambda_N = 0d0
        logical :: predicted = .false.                                           ! rtgroups_predict, not solved
    end type RTGResult

    ! warm state of one bin in one cell (dust_interface keeps one per cell and bin)
    type RTGState
        logical :: valid = .false., discrete = .false.
        real(dp) :: Zstar = 0d0, wlo = 0d0, whi = 0d0        ! last root; discrete window ends
        real(dp) :: Zsigma = 0d0                             ! last width (start of the chord of a wide P(Z))
        ! rtgroups_predict anchor: the last solve of this bin in this cell with sens (inputs, the
        ! per-group photoemission kernel at its Z*, its result with the sensitivities)
        logical :: anchored = .false.
        real(dp) :: aT = 0d0, ane = 0d0, anion(3) = 0d0, acred = 0d0, aG0 = 0d0, akbg = 0d0
        real(dp) :: aNp(RTG_MAXG) = 0d0, akstar(RTG_MAXG) = 0d0
        real(dp) :: adla(RTG_NRI, 3) = 0d0   ! discrete P(Z): d ln alpha_i / d ln S, N, T
        real(dp) :: alg(RTG_NRI) = 0d0       ! wide P(Z): its Gaussian estimate of ln alpha_i
        type(RTGResult) :: ar
        real(dp) :: Dratio(4) = 1d0          ! rtgroups_coulomb_ratio at the last full solve (impactor charges -1, 1, 2, 3)
    end type RTGState

    integer, parameter :: RI_ATOMIC(RTG_NRI) = (/1, 2, 6, 7, 8, 10, 12, 14, 16, 26/)   ! atomic numbers of the RI_* ions

    type(RTGroupTable), allocatable, save :: rtg(:)
    real(dp), save :: ri_hlnm(RTG_NRI) = 0d0          ! 0.5 ln(RI_MASS), set at init

    ! per-call cell state (set_cell); ions 1 (H+) and 2 (He+) both have z = 1, hence the
    ! same tau and the same J~(Z, 1, tau)
    real(dp) :: c_T, c_ne, c_ve, c_nion(3), c_zion(3), c_vion(3), c_G0
    real(dp) :: c_tau_e, c_tau_ion(3), c_arr_e, c_arr_ion(3), c_vh, c_dhot
    real(dp) :: c_geoA                ! pi a^2 (8kT / pi m_u)^1/2 (ion capture per J~ at z = 1)
    real(dp), allocatable :: c_cN(:)
    integer :: c_bin, c_kh
    ! T-independent part of g of the last GM_N charges of the cell (gfun; set_cell empties it):
    ! U index and weight of Z and Z + 1, photoemission at Z. GM_N covers the pair and up to
    ! SIGMA_CHORD_ITER chords of a wide P(Z), re-evaluated at T(1 + dlnT)
    integer, parameter :: GM_N = 2 + 2*SIGMA_CHORD_ITER
    real(dp) :: gm_Z(GM_N), gm_w(GM_N), gm_w1(GM_N), gm_photo(GM_N)
    integer :: gm_i(GM_N), gm_i1(GM_N), gm_used = 0, gm_next = 1

    ! discrete window: per-integer quantities, indexed by charge - b_base
    integer :: b_base
    real(dp), dimension(NBUF) :: b_up, b_dn, b_gam, b_lam, b_lamd, b_dup, b_ddn, b_dlam, b_A, b_P, b_ph, b_dA
    ! recombination fit work
    real(dp), allocatable, save :: f_Zn(:), f_lnP(:), f_lnA(:), f_v(:), f_ps(:), f_row(:)
    integer :: np_
    real(dp) :: pc_p(MAXP), pc_q(MAXP), pc_vals(5,MAXP), pc_c(5,MAXP), pc_ci(6,MAXP), pc_cd(5,MAXP), pc_cA(5,MAXP)
    integer :: pc_kind(MAXP), pc_nc(MAXP), pc_ncA(MAXP), pc_order(MAXP)

    interface
        real(c_double) function expm1(x) bind(c, name='expm1')
            import :: c_double
            real(c_double), value :: x
        end function expm1
        real(c_double) function log1p(x) bind(c, name='log1p')
            import :: c_double
            real(c_double), value :: x
        end function log1p
    end interface

    private
    public :: compute_Coulomb_focusing_ions, &
              init_dust_charging_isrf, isrf_lookup, ISRF_ND, ISRF_DCHARGE
    public :: init_dust_charging_rtgroups, rtgroups_set_group_energies, rtgroups_solve_bin, &
              rtgroups_refine_alpha, rtgroups_uniform_lookup, rtgroups_predict, &
              rtgroups_coulomb_ratio
    public :: dust_charge_moments, coul_term, rtgroups_table_bytes
    public :: RTGState, RTGResult, RTG_NRI, RTG_ALPHA_NONE, RTG_ALPHA_GAUSS, RI_ATOMIC, &
              RTG_RECOMB_IMPORTANT

    contains

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

    ! ======================================================== charging_model = 'WDB06rt'

    ! ------------------------------------------------------------------ init
    subroutine init_dust_charging_rtgroups(groupL0, groupL1, nGroups)
        ! Read the tables of every dust bin and check them against the run's
        ! photon groups and grain sizes. Called from rt_init once groups exist.
        use amr_commons, only: myid
        implicit none
        integer, intent(in) :: nGroups
        real(dp), intent(in) :: groupL0(nGroups), groupL1(nGroups)
        integer :: ii, g, it, istat, ver, iu, k, nzmax
        character(len=40) :: dustlabel
        character(len=256) :: fname
        character(len=512) :: line
        type(RTGroupTable) :: t

        if (allocated(rtg)) deallocate(rtg)
        allocate(rtg(1:ndust))
        nzmax = 0
        do ii = 1, ndust
            write(dustlabel, '(A,I2.2)') 'DustBin_', ii
            if (len_trim(charging_table_label(ii)) > 0) dustlabel = charging_table_label(ii)
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
            allocate(t%KgT(t%nG,t%nU), t%HgT(t%nG,t%nU), t%mKT(t%nG,t%nU), t%mHT(t%nG,t%nU), t%egy_used(t%nG))
            allocate(t%Kbg(t%nU), t%Hbg(t%nU), t%mKb(t%nU), t%mHb(t%nU), t%edge0(t%nG), t%edge1(t%nG))
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
            if (ver >= 2) then
                read(iu, *) t%Kbg
                read(iu, *) t%Hbg
            else
                t%Kbg = 0d0
                t%Hbg = 0d0
                if (myid == 1) write(*,*) 'WDB06rt: ', trim(fname), ' predates the UV-background kernels'
            end if
            if (ver >= 3) read(iu, *) t%Zth
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
            t%stick_pos = 0.5d0*(1d0 - exp(-t%a/1d-7))
            t%stick_neg = t%stick_pos / (1d0 + exp(20d0 - 468d0*(t%a/1d-7)**3))
            call pchip_slopes(t%Kbg, t%mKb, t%nU)
            call pchip_slopes(t%Hbg, t%mHb, t%nU)
            t%nZ = nint(t%Zmax - t%Zmin) + 1
            nzmax = max(nzmax, t%nZ)
            rtg(ii) = t
            deallocate(t%L0, t%L1, t%U, t%lTs, t%egy, t%K, t%H, t%lTh, t%delta, t%Ptr, &
                       t%fcool, t%Ebar, t%KgT, t%HgT, t%mKT, t%mHT, t%egy_used, t%Kbg, t%Hbg, t%mKb, t%mHb, &
                       t%edge0, t%edge1)

            ! charges where IP(a, Z) crosses the group edges / Mathis breaks (kinks of g), the
            ! recombination thresholds (when the table predates them), and the U index of each
            ! integer charge
            c_bin = ii
            do g = 1, rtg(ii)%nG
                rtg(ii)%edge0(g) = ip_edge(ii, rtg(ii)%L0(g))
                rtg(ii)%edge1(g) = ip_edge(ii, rtg(ii)%L1(g))
            end do
            do k = 1, 4
                rtg(ii)%edgebg(k) = ip_edge(ii, MATHIS_BREAKS(k))
            end do
            if (ver < 3) call recomb_thresholds(ii)
            allocate(rtg(ii)%iuz(rtg(ii)%nZ), rtg(ii)%wuz(rtg(ii)%nZ))
            rtg(ii)%jpos = rtg(ii)%nU + 1
            do k = rtg(ii)%nU, 1, -1
                if (rtg(ii)%U(k) > 0d0) rtg(ii)%jpos = k
            end do
            do k = 1, rtg(ii)%nZ
                call u_index((rtg(ii)%Zmin + dble(k - 1))*RTG_E2EVCM/rtg(ii)%a, rtg(ii)%iuz(k), rtg(ii)%wuz(k))
            end do
            if (myid == 1) write(*,'(A,A,A,ES12.5,A,F9.1,A,F9.1,A,I2)') ' WDB06rt: read ', trim(dustlabel), &
                ', a = ', rtg(ii)%a, ' cm, Zmin = ', rtg(ii)%Zmin, ', Zmax = ', rtg(ii)%Zmax, ', format ', ver
        end do
        allocate(c_cN(1:nGroups))
        ri_hlnm = 0.5d0*log(RI_MASS)
        allocate(f_Zn(nzmax + 2), f_lnP(nzmax + 2), f_lnA(nzmax + 2), f_v(nzmax + 2), f_ps(0:nzmax + 2))
        allocate(f_row(maxval(rtg(:)%nU)))
    end subroutine init_dust_charging_rtgroups

    subroutine pchip_slopes(y, m, n)
        ! Fritsch-Carlson monotone slopes on a unit-spaced grid (pyCALIMA _pchip_slopes)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: y(n)
        real(dp), intent(out) :: m(n)
        real(dp) :: d0, d1
        integer :: j
        m(1) = y(2) - y(1)
        m(n) = y(n) - y(n-1)
        do j = 2, n - 1
            d0 = y(j) - y(j-1)
            d1 = y(j+1) - y(j)
            if (d0*d1 > 0d0) then
                m(j) = 2d0*d0*d1/(d0 + d1)
            else
                m(j) = 0d0
            end if
        end do
    end subroutine pchip_slopes

    real(dp) function grain_ip(ii, Z)
        ! IP(a, Z) [eV]: valence IP for Z >= 0, electron affinity of charge Z + 1 below
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Z
        real(dp) :: a
        a = rtg(ii)%a
        if (Z >= 0d0) then
            grain_ip = rtg(ii)%W + RTG_ESTATC**2/a*((Z + 0.5d0) + (Z + 2d0)*(0.3d-8/a))/RTG_EV2ERG
        else if (rtg(ii)%material == 1) then
            grain_ip = GRAPHITE_WF + RTG_ESTATC**2/a*(((Z + 1d0) - 0.5d0) - (4d-8/(a + 7d-8)))/RTG_EV2ERG
        else
            grain_ip = SILICATE_WF - SILICATE_GAP + RTG_ESTATC**2/a*((Z + 1d0) - 0.5d0)/RTG_EV2ERG
        end if
    end function grain_ip

    real(dp) function electron_affinity(ii, Z)
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Z
        real(dp) :: a
        a = rtg(ii)%a
        if (rtg(ii)%material == 1) then
            electron_affinity = GRAPHITE_WF + RTG_ESTATC**2/a*((Z - 0.5d0) - (4d-8/(a + 7d-8)))/RTG_EV2ERG
        else
            electron_affinity = SILICATE_WF - SILICATE_GAP + RTG_ESTATC**2/a*(Z - 0.5d0)/RTG_EV2ERG
        end if
    end function electron_affinity

    real(dp) function ip_edge(ii, E)
        ! largest charge with IP(a, Z) <= E, if strictly inside (Zmin, Zmax); else NONE
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: E
        integer :: k, kl
        kl = -1
        do k = 0, rtg(ii)%nZ - 1
            if (grain_ip(ii, rtg(ii)%Zmin + dble(k)) <= E) kl = k
        end do
        ip_edge = NONE
        if (kl > 0 .and. kl < rtg(ii)%nZ - 1) ip_edge = rtg(ii)%Zmin + dble(kl)
    end function ip_edge

    subroutine recomb_thresholds(ii)
        ! largest Z with IP(a, Z) <= IP(X) (case A), Zmin - 1 where none (pyCALIMA recombination_thresholds)
        implicit none
        integer, intent(in) :: ii
        integer :: k, j
        do j = 1, RTG_NRI
            rtg(ii)%Zth(j) = rtg(ii)%Zmin - 1d0
            do k = 0, rtg(ii)%nZ - 1
                if (grain_ip(ii, rtg(ii)%Zmin + dble(k)) <= RI_IP(j)) rtg(ii)%Zth(j) = rtg(ii)%Zmin + dble(k)
            end do
        end do
    end subroutine recomb_thresholds

    subroutine rtgroups_set_group_energies(ii, egy)
        ! Per-group K_g(U), H_g(U) at the run's mean group energies: linear in
        ! the table's group mean energy between the two bracketing shapes; and
        ! their PCHIP slopes in the U-node index.
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: egy(:)
        integer :: g, k, j, n
        real(dp) :: eg, w
        real(dp), allocatable :: col(:), m(:)

        n = rtg(ii)%nTs
        allocate(col(rtg(ii)%nU), m(rtg(ii)%nU))
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
            col = (1d0 - w)*rtg(ii)%K(:,j,g) + w*rtg(ii)%K(:,j+1,g)
            call pchip_slopes(col, m, rtg(ii)%nU)
            rtg(ii)%KgT(g,:) = col
            rtg(ii)%mKT(g,:) = m
            col = (1d0 - w)*rtg(ii)%H(:,j,g) + w*rtg(ii)%H(:,j+1,g)
            call pchip_slopes(col, m, rtg(ii)%nU)
            rtg(ii)%HgT(g,:) = col
            rtg(ii)%mHT(g,:) = m
        end do
        rtg(ii)%egy_used = egy(1:rtg(ii)%nG)
        deallocate(col, m)
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

    subroutine anchor_kernels(Z, st)
        ! the per-group photoemission kernels at charge Z (and the background's), unweighted: the
        ! photoemission rate at the anchor's Z* for any photon densities (rtgroups_predict)
        implicit none
        real(dp), intent(in) :: Z
        type(RTGState), intent(inout) :: st
        real(dp) :: w, w2, w3, h00, h10, h01, h11
        integer :: i, g
        call u_of(Z, i, w)
        w2 = w*w
        w3 = w*w*w
        h00 = 2d0*w3 - 3d0*w2 + 1d0
        h10 = w3 - 2d0*w2 + w
        h01 = 3d0*w2 - 2d0*w3
        h11 = w3 - w2
        st%akstar = 0d0
        do g = 1, rtg(c_bin)%nG
            st%akstar(g) = h00*rtg(c_bin)%KgT(g,i) + h10*rtg(c_bin)%mKT(g,i) &
                           + h01*rtg(c_bin)%KgT(g,i+1) + h11*rtg(c_bin)%mKT(g,i+1)
        end do
        st%akbg = h00*rtg(c_bin)%Kbg(i) + h10*rtg(c_bin)%mKb(i) + h01*rtg(c_bin)%Kbg(i+1) + h11*rtg(c_bin)%mKb(i+1)
    end subroutine anchor_kernels

    subroutine rtgroups_predict(ii, st, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, tol, tol_shape, r, ok)
        ! First-order prediction of bin ii in this cell from its anchor st (the last
        ! rtgroups_solve_bin with sens), or ok = .false. where it must be solved again
        ! (pyCALIMA rt_group_charging.predict). Changes since the anchor: d ln S of the
        ! photoemission rate at the anchor's Z* (all groups and the background, with its kernels),
        ! d ln N of n_e, d ln T; refused beyond tol (also for the change of an ion density beyond
        ! d ln N), if a group's photon density changed by more than tol_shape beyond d ln S, if
        ! sigma is on the other side of SIGMA_DISCRETE than the anchor's, or, for a wide P(Z), if
        ! Z* + 1/2 (of the anchor or of the prediction) is within PREDICT_WIDE_ZERO + 1/2 of zero
        ! charge or on either side of it, or moves by more than max(1, PREDICT_WIDE_DZ sigma).
        ! A discrete P(Z) is the linearisation, alpha too (in ln); so is a wide one for a small
        ! change, else it is corrected at the new inputs (a Newton step of g, Gamma and Lambda at
        ! Z* + 1/2); its alpha is scaled by the change of its Gaussian estimate. r: Lambda includes
        ! the autoionisation; d/d ln T is the anchor's.
        implicit none
        integer, intent(in) :: ii
        type(RTGState), intent(in) :: st
        real(dp), intent(in) :: Np(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, tol, tol_shape
        type(RTGResult), intent(out) :: r
        logical, intent(out) :: ok
        real(dp) :: P0, P1, rS, rN, rT, dS, dN, dT, nion(3), m, z0, z1, y, zs, lg(RTG_NRI)
        logical :: linear
        real(dp), save :: tol_c = -1d0, tsh_c = -1d0, et = 1d0, es = 1d0     ! exp(tol), exp(tol_shape)
        integer :: g, j, nG
        ! the tests |ln(x/y) - d| <= tol as bounds on ratios: no logarithm per group or ion
        ok = .false.
        if (.not. st%anchored .or. ne <= 0d0) return
        nG = rtg(ii)%nG
        if (tol /= tol_c .or. tol_shape /= tsh_c) then
            tol_c = tol
            tsh_c = tol_shape
            et = exp(tol)
            es = exp(tol_shape)
        end if
        P0 = st%acred*sum(st%aNp(1:nG)*st%akstar(1:nG)) + st%aG0*st%akbg
        P1 = c_red*sum(Np(1:nG)*st%akstar(1:nG)) + G0_bg*st%akbg
        if (P0 > 0d0 .and. P1 > 0d0) then
            rS = P1/P0
        else if (P0 <= 0d0 .and. P1 <= 0d0) then
            rS = 1d0
        else
            return
        end if
        rN = ne/st%ane
        rT = T/st%aT
        if (rS > et .or. rS*et < 1d0 .or. rN > et .or. rN*et < 1d0 .or. rT > et .or. rT*et < 1d0) return
        nion = (/n_Hp, n_Hep, n_Hepp/)
        do j = 1, 3
            if (max(nion(j), st%anion(j)) > 1d-6*ne) then
                if (nion(j) <= 0d0 .or. st%anion(j) <= 0d0) return
                y = st%anion(j)*rN
                if (nion(j) > y*et .or. nion(j)*et < y) return
            end if
        end do
        m = 0d0
        do g = 1, nG
            m = max(m, Np(g), st%aNp(g))
        end do
        m = 1d-30*max(m, 1d-300)
        do g = 1, nG
            if (max(Np(g), st%aNp(g)) > m) then
                if (Np(g) <= 0d0 .or. st%aNp(g) <= 0d0) return
                y = st%aNp(g)*rS
                if (Np(g) > y*es .or. Np(g)*es < y) return
            end if
        end do
        dS = log(rS)
        dN = log(rN)
        dT = log(rT)
        associate (a => st%ar)
            r%Zmean = a%Zmean + a%dZmean_S*dS + a%dZmean_N*dN + a%dZmean*dT
            r%Zsigma = a%Zsigma + a%dZsigma_S*dS + a%dZsigma_N*dN + a%dZsigma*dT
            r%Gamma = a%Gamma + a%dGamma_S*dS + a%dGamma_N*dN + a%dGamma*dT
            r%Lambda = (a%Lambda + a%Lambda_auto) + a%dLambda_S*dS + a%dLambda_N*dN + a%dLambda*dT
            r%Lambda_auto = 0d0
            r%Zstar = a%Zstar + (r%Zmean - a%Zmean)
            r%discrete = a%discrete
            r%dZmean = a%dZmean
            r%dZsigma = a%dZsigma
            r%dGamma = a%dGamma
            r%dLambda = a%dLambda
            r%alpha = a%alpha
            r%predicted = .true.
            if ((r%Zsigma < SIGMA_DISCRETE) .neqv. a%discrete) return
            if (.not. a%discrete) then
                z0 = a%Zstar + 0.5d0
                z1 = r%Zstar + 0.5d0
                linear = max(abs(dS), abs(dN), abs(dT)) <= PREDICT_LIN_TOL .and. z0*z1 > 0d0 .and. &
                         min(abs(z0), abs(z1)) >= PREDICT_LIN_ZERO + 0.5d0 .and. &
                         abs(z1 - z0) <= max(1d0, PREDICT_LIN_DZ*a%Zsigma)
                if (.not. linear .or. any(a%alpha > 0d0)) call set_cell(ii, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg)
                if (.not. linear) then
                    ! corrected at the new inputs: a Newton step of g from the linear Z* (g' = -1/sigma^2,
                    ! sigma linear), then Gamma and Lambda at Z* + 1/2 (the wide closure of the solve)
                    zs = min(max(r%Zstar + r%Zsigma*r%Zsigma*gfun(r%Zstar), rtg(ii)%Zmin), rtg(ii)%Zmax - 1d0)
                    z1 = zs + 0.5d0
                    if (z0*z1 <= 0d0 .or. min(abs(z0), abs(z1)) < PREDICT_WIDE_ZERO + 0.5d0) return
                    if (abs(z1 - z0) > max(1d0, PREDICT_WIDE_DZ*a%Zsigma)) return
                    call heating_cooling(z1, r%Gamma, r%Lambda)
                    r%Zstar = zs
                    r%Zmean = z1
                end if
                ! alpha (fitted or not) times the change of its Gaussian estimate
                if (any(a%alpha > 0d0)) then
                    call recomb_gauss(r%Zstar, r%Zsigma, lg)
                    do j = 1, RTG_NRI
                        if (a%alpha(j) > 0d0) r%alpha(j) = a%alpha(j)*exp(lg(j) - st%alg(j))
                    end do
                end if
            else
                do j = 1, RTG_NRI
                    if (a%alpha(j) > 0d0) r%alpha(j) = a%alpha(j)*exp(st%adla(j, 1)*dS + st%adla(j, 2)*dN &
                                                                      + st%adla(j, 3)*dT)
                end do
            end if
        end associate
        ok = .true.
    end subroutine rtgroups_predict

    subroutine u_of(Z, i, w)
        ! U index of charge Z: cached for the integers
        implicit none
        real(dp), intent(in) :: Z
        integer, intent(out) :: i
        real(dp), intent(out) :: w
        integer :: k
        if (Z == aint(Z) .and. Z >= rtg(c_bin)%Zmin .and. Z <= rtg(c_bin)%Zmax) then
            k = nint(Z - rtg(c_bin)%Zmin) + 1
            i = rtg(c_bin)%iuz(k)
            w = rtg(c_bin)%wuz(k)
        else
            call u_index(Z*RTG_E2EVCM/rtg(c_bin)%a, i, w)
        end if
    end subroutine u_of

    subroutine hot_index(T)
        ! log T index and weight of the hot-gas tables, and d/d ln T factor, once per cell
        implicit none
        real(dp), intent(in) :: T
        real(dp) :: y, dl, fT, lT
        integer :: n
        n = rtg(c_bin)%nTh
        lT = log10(T)
        y = min(max(lT, rtg(c_bin)%lTh(1)), rtg(c_bin)%lTh(n))
        dl = (rtg(c_bin)%lTh(n) - rtg(c_bin)%lTh(1)) / dble(n - 1)
        fT = (y - rtg(c_bin)%lTh(1)) / dl
        c_kh = min(int(fT), n - 2)
        c_vh = fT - dble(c_kh)
        c_kh = c_kh + 1
        c_dhot = 0d0
        if (rtg(c_bin)%lTh(1) < lT .and. lT < rtg(c_bin)%lTh(n)) c_dhot = dl
    end subroutine hot_index

    real(dp) function hot(tab, i, w)
        implicit none
        real(dp), intent(in) :: tab(:,:)
        integer, intent(in) :: i
        real(dp), intent(in) :: w
        hot = (1d0 - c_vh)*((1d0 - w)*tab(i,c_kh) + w*tab(i+1,c_kh)) &
              + c_vh*((1d0 - w)*tab(i,c_kh+1) + w*tab(i+1,c_kh+1))
    end function hot

    subroutine hot_and_dlnT(tab, i, w, v, dv)
        ! hot and its d/d ln T (the slope in log10 T, zero where log T is clamped) from one
        ! pair of U interpolations
        implicit none
        real(dp), intent(in) :: tab(:,:)
        integer, intent(in) :: i
        real(dp), intent(in) :: w
        real(dp), intent(out) :: v, dv
        real(dp) :: a0, a1
        a0 = (1d0 - w)*tab(i,c_kh) + w*tab(i+1,c_kh)
        a1 = (1d0 - w)*tab(i,c_kh+1) + w*tab(i+1,c_kh+1)
        v = (1d0 - c_vh)*a0 + c_vh*a1
        dv = 0d0
        if (c_dhot > 0d0) dv = (a1 - a0) / c_dhot / log(10d0)
    end subroutine hot_and_dlnT

    subroutine kernels(i, w, photo, heat, want_heat)
        ! PCHIP photoemission (and heating) kernels summed over the groups and the background
        implicit none
        integer, intent(in) :: i
        real(dp), intent(in) :: w
        real(dp), intent(out) :: photo, heat
        logical, intent(in) :: want_heat
        real(dp) :: w2, w3, h00, h10, h01, h11, s
        integer :: g
        w2 = w*w
        w3 = w*w*w
        h00 = 2d0*w3 - 3d0*w2 + 1d0
        h10 = w3 - 2d0*w2 + w
        h01 = 3d0*w2 - 2d0*w3
        h11 = w3 - w2
        s = 0d0
        do g = 1, rtg(c_bin)%nG
            s = s + c_cN(g)*(h00*rtg(c_bin)%KgT(g,i) + h10*rtg(c_bin)%mKT(g,i) &
                             + h01*rtg(c_bin)%KgT(g,i+1) + h11*rtg(c_bin)%mKT(g,i+1))
        end do
        photo = max(s + c_G0*(h00*rtg(c_bin)%Kbg(i) + h10*rtg(c_bin)%mKb(i) &
                              + h01*rtg(c_bin)%Kbg(i+1) + h11*rtg(c_bin)%mKb(i+1)), 0d0)
        heat = 0d0
        if (.not. want_heat) return
        s = 0d0
        do g = 1, rtg(c_bin)%nG
            s = s + c_cN(g)*(h00*rtg(c_bin)%HgT(g,i) + h10*rtg(c_bin)%mHT(g,i) &
                             + h01*rtg(c_bin)%HgT(g,i+1) + h11*rtg(c_bin)%mHT(g,i+1))
        end do
        heat = max(s + c_G0*(h00*rtg(c_bin)%Hbg(i) + h10*rtg(c_bin)%mHb(i) &
                             + h01*rtg(c_bin)%Hbg(i+1) + h11*rtg(c_bin)%mHb(i+1)), 0d0)
    end subroutine kernels

    ! -------------------------------------------------------------- physics
    real(dp) function ds87_tau(q)
        implicit none
        real(dp), intent(in) :: q
        ds87_tau = max(rtg(c_bin)%a*RTG_KB*c_T / max(q*q*RTG_ESTATC*RTG_ESTATC, RTG_TINY), RTG_TINY)
    end function ds87_tau

    real(dp) function ds87_J(Z, q, tau, dlnJ)
        ! Draine & Sutin (1987) J-tilde of grain charge Z for an impactor of charge q, theta_nu =
        ! nu/(1+nu^-1/2) (their eq. 3.6); tau = ds87_tau(q) of the cell; dlnJ: d ln J / d ln tau
        ! (pyCALIMA _ds87_J_dlnT)
        implicit none
        real(dp), intent(in) :: Z, q, tau
        real(dp), intent(out), optional :: dlnJ
        real(dp) :: nu, r, A, d, s, th, nup
        nu = Z / q
        if (nu == 0d0) then
            r = sqrt(RTG_PI/(2d0*tau))
            ds87_J = 1d0 + r
            if (present(dlnJ)) dlnJ = -0.5d0*r/(1d0 + r)
        else if (nu < 0d0) then
            A = 1d0 - nu/tau
            d = max(tau - 2d0*nu, RTG_TINY)
            s = sqrt(2d0/d)
            ds87_J = A*(1d0 + s)
            if (present(dlnJ)) dlnJ = (nu/tau)/A - 0.5d0*s*tau/d/(1d0 + s)
        else
            nup = max(nu, RTG_TINY)
            th = nup / (1d0 + 1d0/sqrt(nup))
            s = 1d0/sqrt(4d0*tau + 3d0*nup)
            ds87_J = (1d0 + s)**2 * exp(-th/tau)
            if (present(dlnJ)) dlnJ = -4d0*tau*s**3/(1d0 + s) + th/tau
        end if
    end function ds87_J

    real(dp) function ds87_lnJ(Z, q, tau)
        ! ln of ds87_J, in log form (no underflow of exp(-theta/tau) at large positive nu)
        implicit none
        real(dp), intent(in) :: Z, q, tau
        real(dp) :: nu
        nu = Z / q
        if (nu > 0d0) then
            ds87_lnJ = 2d0*log(1d0 + 1d0/sqrt(4d0*tau + 3d0*nu)) - nu/(1d0 + 1d0/sqrt(nu))/tau
        else if (nu < 0d0) then
            ds87_lnJ = log((1d0 - nu/tau)*(1d0 + sqrt(2d0/(tau - 2d0*nu))))
        else
            ds87_lnJ = log(1d0 + sqrt(RTG_PI/(2d0*tau)))
        end if
    end function ds87_lnJ

    real(dp) function ds87_L(Z, tau, dlnL)
        ! Draine & Sutin (1987) Lambda-tilde of grain charge Z for electrons (q = -1), tau = ds87_tau(-1)
        ! of the cell; dlnL: d ln L / d ln tau (pyCALIMA _ds87_lambda_dlnT)
        implicit none
        real(dp), intent(in) :: Z, tau
        real(dp), intent(out), optional :: dlnL
        real(dp) :: nu, r, A, s, th, nup
        nu = -Z
        if (nu == 0d0) then
            r = sqrt(RTG_PI/(2d0*tau))
            ds87_L = 2d0 + 1.5d0*r
            if (present(dlnL)) dlnL = -0.75d0*r/ds87_L
        else if (nu < 0d0) then
            A = 2d0 - nu/tau
            s = 1d0/sqrt(max(tau - nu, RTG_TINY))
            ds87_L = A*(1d0 + s)
            if (present(dlnL)) dlnL = (nu/tau)/A - 0.5d0*tau*s**3/(1d0 + s)
        else
            nup = max(nu, RTG_TINY)
            th = nup / (1d0 + 1d0/sqrt(nup))
            A = 2d0 + nup/tau
            s = 1d0/sqrt(1.5d0/tau + 3d0*nup)
            ds87_L = A*(1d0 + s)*exp(-th/tau)
            if (present(dlnL)) dlnL = (-nup/tau)/A + 0.75d0*s**3/tau/(1d0 + s) + th/tau
        end if
    end function ds87_L

    real(dp) function sticking(Z)
        ! WD01 eqs. 27-30 with the continuous extension used by the solver
        implicit none
        real(dp), intent(in) :: Z
        if (Z <= rtg(c_bin)%Zmin) then
            sticking = 0d0
        else if (Z > 0d0) then
            sticking = rtg(c_bin)%stick_pos
        else
            sticking = rtg(c_bin)%stick_neg
        end if
    end function sticking

    real(dp) function ion_capture(Z)
        ! H+, He+ (one J~ for z = 1) and He++, summed in that order
        implicit none
        real(dp), intent(in) :: Z
        real(dp) :: j1
        ion_capture = 0d0
        if (c_nion(1) > 0d0 .or. c_nion(2) > 0d0) then
            j1 = ds87_J(Z, 1d0, c_tau_ion(1))
            if (c_nion(1) > 0d0) ion_capture = ion_capture + c_arr_ion(1)*j1
            if (c_nion(2) > 0d0) ion_capture = ion_capture + c_arr_ion(2)*j1
        end if
        if (c_nion(3) > 0d0) ion_capture = ion_capture + c_arr_ion(3)*ds87_J(Z, c_zion(3), c_tau_ion(3))
    end function ion_capture

    real(dp) function gfun(Z)
        ! ln Up(Z) - ln Down(Z+1): Up = photoemission + ion capture + gas-electron secondaries at Z,
        ! Down = electron capture at Z + 1. The U indices of Z and Z + 1 and the photoemission at
        ! Z do not depend on T; those of the last GM_N charges of the cell are kept, so g at the
        ! same charges at T(1 + dlnT) (the d/d ln T of a wide P(Z)) costs only its T-dependent part
        implicit none
        real(dp), intent(in) :: Z
        real(dp) :: w, w1, photo, heat, up, down
        integer :: i, i1, m
        do m = 1, gm_used
            if (gm_Z(m) == Z) exit
        end do
        if (m <= gm_used) then
            i = gm_i(m)
            w = gm_w(m)
            i1 = gm_i1(m)
            w1 = gm_w1(m)
            photo = gm_photo(m)
        else
            call u_of(Z, i, w)
            call kernels(i, w, photo, heat, .false.)
            call u_of(Z + 1d0, i1, w1)
            m = gm_next
            gm_Z(m) = Z
            gm_i(m) = i
            gm_w(m) = w
            gm_i1(m) = i1
            gm_w1(m) = w1
            gm_photo(m) = photo
            gm_next = mod(m, GM_N) + 1
            gm_used = max(gm_used, m)
        end if
        up = photo + ion_capture(Z) + c_arr_e*ds87_J(Z, -1d0, c_tau_e)*hot(rtg(c_bin)%delta, i, w)
        down = c_arr_e*ds87_J(Z + 1d0, -1d0, c_tau_e)*sticking(Z + 1d0)*(1d0 - hot(rtg(c_bin)%Ptr, i1, w1))
        gfun = log(max(up, RTG_TINY)) - log(max(down, RTG_TINY))
    end function gfun

    subroutine heating_cooling(Z, gam, lam)
        implicit none
        real(dp), intent(in) :: Z
        real(dp), intent(out) :: gam, lam
        real(dp) :: w, photo, heat, kT
        integer :: i
        call u_of(Z, i, w)
        call kernels(i, w, photo, heat, .true.)
        gam = RTG_EV2ERG*heat
        kT = RTG_KB*c_T
        lam = c_arr_e*sticking(Z)*ds87_L(Z, c_tau_e)*kT*hot(rtg(c_bin)%fcool, i, w)
        lam = lam - c_arr_e*ds87_J(Z, -1d0, c_tau_e)*hot(rtg(c_bin)%delta, i, w)*hot(rtg(c_bin)%Ebar, i, w)*RTG_EV2ERG
    end subroutine heating_cooling

    real(dp) function ln_capture_z1(Z)
        ! ln of pi a^2 (8kT / pi m_u)^1/2 J~(Z, z=1), log form (no underflow)
        implicit none
        real(dp), intent(in) :: Z
        ln_capture_z1 = log(c_geoA) + ds87_lnJ(Z, 1d0, c_tau_ion(1))
    end function ln_capture_z1

    ! ------------------------------------------------------------- cell state
    subroutine set_cell(ii, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg)
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg
        c_bin = ii
        c_cN(1:rtg(ii)%nG) = c_red*Np(1:rtg(ii)%nG)
        gm_used = 0
        gm_next = 1
        c_ne = ne
        c_G0 = G0_bg
        c_nion = (/ n_Hp, n_Hep, n_Hepp /)
        c_zion = (/ 1d0, 1d0, 2d0 /)
        call set_cell_T(T)
    end subroutine set_cell

    subroutine set_cell_T(T)
        ! the temperature-dependent part of the cell state
        implicit none
        real(dp), intent(in) :: T
        integer :: k
        c_T = T
        c_ve = sqrt(8d0*RTG_KB*T/(RTG_PI*RTG_ME))
        c_vion(1) = sqrt(8d0*RTG_KB*T/(RTG_PI*M_HION))
        c_vion(2) = sqrt(8d0*RTG_KB*T/(RTG_PI*M_HEION))
        c_vion(3) = c_vion(2)
        c_tau_e = ds87_tau(-1d0)
        do k = 1, 3
            c_tau_ion(k) = ds87_tau(c_zion(k))
        end do
        c_arr_e = RTG_PI*rtg(c_bin)%a*rtg(c_bin)%a*c_ne*c_ve
        do k = 1, 3
            c_arr_ion(k) = RTG_PI*rtg(c_bin)%a*rtg(c_bin)%a*c_nion(k)*c_vion(k)
        end do
        c_geoA = RTG_PI*rtg(c_bin)%a**2*sqrt(8d0*RTG_KB*c_T/(RTG_PI*RTG_MU))
        call hot_index(T)
    end subroutine set_cell_T

    ! ---------------------------------------------------------------- root
    real(dp) function illinois(a0, b0, fa0, fb0)
        ! root of gfun in [a0, b0], gfun(a0) > 0 > gfun(b0)
        implicit none
        real(dp), intent(in) :: a0, b0, fa0, fb0
        real(dp) :: a_, b_, fa, fb, fz, Zs
        integer :: side, it
        a_ = a0; b_ = b0; fa = fa0; fb = fb0; side = 0
        Zs = b_
        do it = 1, 200
            Zs = (a_*fb - b_*fa)/(fb - fa)
            if (abs(b_ - a_) < 1d-7*max(1d0, abs(Zs))) exit
            fz = gfun(Zs)
            if (fz == 0d0) exit
            if (fz*fb > 0d0) then
                b_ = Zs; fb = fz
                if (side == -1) fa = 0.5d0*fa
                side = -1
            else
                a_ = Zs; fa = fz
                if (side == 1) fb = 0.5d0*fb
                side = 1
            end if
        end do
        illinois = Zs
    end function illinois

    real(dp) function root_from_guess(Zg, Zlo, Zhi)
        ! bracket the root by steps 1, 2, 4, ... from Zg, then Illinois
        implicit none
        real(dp), intent(in) :: Zg, Zlo, Zhi
        real(dp) :: a_, b_, fa, fb, fg, step
        fg = gfun(Zg)
        root_from_guess = Zg
        if (fg == 0d0) return
        step = 1d0
        if (fg > 0d0) then              ! root above Zg
            a_ = Zg; fa = fg
            do
                b_ = min(a_ + step, Zhi)
                fb = gfun(b_)
                if (fb == 0d0) then
                    root_from_guess = b_
                    return
                end if
                if (fb < 0d0) then
                    root_from_guess = illinois(a_, b_, fa, fb)
                    return
                end if
                if (b_ == Zhi) then
                    root_from_guess = Zhi
                    return
                end if
                a_ = b_; fa = fb
                step = 2d0*step
            end do
        end if
        b_ = Zg; fb = fg                ! root below Zg
        do
            a_ = max(b_ - step, Zlo)
            fa = gfun(a_)
            if (fa == 0d0) then
                root_from_guess = a_
                return
            end if
            if (fa > 0d0) then
                root_from_guess = illinois(a_, b_, fa, fb)
                return
            end if
            if (a_ == Zlo) then
                root_from_guess = Zlo
                return
            end if
            b_ = a_; fb = fa
            step = 2d0*step
        end do
    end function root_from_guess

    subroutine slope_at(Zs, sigma, pair)
        ! sigma from g' by central difference at Z* +- 1/2, and that pair (z-, z+, g-, g+)
        implicit none
        real(dp), intent(in) :: Zs
        real(dp), intent(out) :: sigma, pair(4)
        real(dp) :: sl
        pair(1) = max(Zs - 0.5d0, rtg(c_bin)%Zmin)
        pair(2) = min(Zs + 0.5d0, rtg(c_bin)%Zmax - 1d0)
        pair(3) = gfun(pair(1))
        pair(4) = gfun(pair(2))
        sl = (pair(4) - pair(3)) / (pair(2) - pair(1))
        sigma = 0d0
        if (sl < 0d0) sigma = sqrt(-1d0/sl)
    end subroutine slope_at

    subroutine chord_sigma(Zs, sigma, chord, have)
        ! pyCALIMA GroupChargingTable._chord: sigma from the chord of g over the integers
        ! floor(Z*) - k + 1 .. floor(Z*) + k, k = nint(sigma), re-centred until k is stable
        ! (sigma in: the starting width); the secant at Z* +- 1/2 is off by up to ~3x next
        ! to a kink of g. have: a chord with g decreasing was found (chord = lo, hi, g(lo), g(hi))
        implicit none
        real(dp), intent(in) :: Zs
        real(dp), intent(inout) :: sigma
        real(dp), intent(out) :: chord(4)
        logical, intent(out) :: have
        real(dp) :: Zlo, Zhi, lo, hi, glo, ghi, sl
        integer :: it, k, kn
        Zlo = rtg(c_bin)%Zmin
        Zhi = rtg(c_bin)%Zmax - 1d0
        have = .false.
        chord = 0d0
        k = -1
        do it = 1, SIGMA_CHORD_ITER
            kn = max(1, nint(sigma))
            if (kn == k) exit
            k = kn
            lo = max(dble(floor(Zs)) - dble(k) + 1d0, Zlo)
            hi = min(dble(floor(Zs)) + dble(k), Zhi)
            if (hi <= lo) exit
            glo = gfun(lo)
            ghi = gfun(hi)
            sl = (ghi - glo)/(hi - lo)
            if (sl >= 0d0) exit
            sigma = sqrt(-1d0/sl)
            chord = (/ lo, hi, glo, ghi /)
            have = .true.
        end do
    end subroutine chord_sigma

    subroutine find_root(Z_guess, st, Zs, sigma, pair)
        ! pyCALIMA GroupChargingTable._root
        implicit none
        real(dp), intent(in) :: Z_guess
        type(RTGState), intent(in) :: st
        real(dp), intent(out) :: Zs, sigma, pair(4)
        real(dp) :: Zlo, Zhi, zm, zp, gm, gpl, sl, gs, zg, glo, ghi
        Zlo = rtg(c_bin)%Zmin
        Zhi = rtg(c_bin)%Zmax - 1d0
        zg = Z_guess
        if (st%valid) then
            zm = max(st%Zstar - 0.5d0, Zlo)
            zp = min(st%Zstar + 0.5d0, Zhi)
            gm = gfun(zm)
            gpl = gfun(zp)
            sl = (gpl - gm) / (zp - zm)
            if (sl < 0d0) then
                Zs = zm - gm/sl
                if (zm <= Zs .and. Zs <= zp) then
                    if (.not. st%discrete) then          ! wide: g is smooth on the scale of the pair
                        sigma = sqrt(-1d0/sl)
                        pair = (/ zm, zp, gm, gpl /)
                        return
                    end if
                    gs = gfun(Zs)
                    if (abs(gs/sl) > 1d-3 .and. gm > 0d0 .and. 0d0 > gpl) Zs = illinois(zm, zp, gm, gpl)
                    if (abs(Zs - st%Zstar) <= 0.01d0) then
                        sigma = sqrt(-1d0/sl)
                        pair = (/ zm, zp, gm, gpl /)
                    else
                        call slope_at(Zs, sigma, pair)
                    end if
                    return
                end if
                zg = min(max(Zs, Zlo + 1d0), Zhi - 1d0)
            end if
        end if
        if (zg > Zlo .and. zg < Zhi) then
            Zs = root_from_guess(zg, Zlo, Zhi)
        else
            glo = gfun(Zlo)
            ghi = gfun(Zhi)
            if (glo <= 0d0) then
                Zs = Zlo
            else if (ghi >= 0d0) then
                Zs = Zhi
            else
                Zs = illinois(Zlo, Zhi, glo, ghi)
            end if
        end if
        call slope_at(Zs, sigma, pair)
    end subroutine find_root

    ! --------------------------------------------------------- discrete P(Z)
    subroutine eval_int(Z, want_d, rates_only)
        ! rates, heating, cooling (and their d/d ln T) of integer charge Z into the window buffers.
        ! Each J~ / Lambda~ and hot-gas row is evaluated once, with its d/d ln T when want_d: the
        ! electrons, z = 1 (H+, He+ and b_A share it) and He++ (only if present); the sums keep the
        ! order of the separate evaluations, so the buffers are the same to the bit. rates_only:
        ! only Up, Down, the photoemission and b_A (and their d/d ln T), for the charges below the
        ! window in the recombination
        implicit none
        real(dp), intent(in) :: Z
        logical, intent(in) :: want_d
        logical, intent(in), optional :: rates_only
        real(dp) :: w, photo, heat, je, dlt, ptr, fc, eb, st, lt, kT, dje, dld, j1, dj1, j3, dj3, ion, dion
        real(dp) :: dde, dpt, dfc, deb
        integer :: i, k
        logical :: full
        full = .true.
        if (present(rates_only)) full = .not. rates_only
        k = nint(Z) - b_base
        call u_of(Z, i, w)
        call kernels(i, w, photo, heat, full)
        j3 = 0d0
        dj3 = 0d0
        if (want_d) then
            je = ds87_J(Z, -1d0, c_tau_e, dje)
            if (full) lt = ds87_L(Z, c_tau_e, dld)
            j1 = ds87_J(Z, 1d0, c_tau_ion(1), dj1)
            if (c_arr_ion(3) /= 0d0 .or. c_nion(3) > 0d0) j3 = ds87_J(Z, c_zion(3), c_tau_ion(3), dj3)
            call hot_and_dlnT(rtg(c_bin)%delta, i, w, dlt, dde)
            call hot_and_dlnT(rtg(c_bin)%Ptr, i, w, ptr, dpt)
            if (full) then
                call hot_and_dlnT(rtg(c_bin)%fcool, i, w, fc, dfc)
                call hot_and_dlnT(rtg(c_bin)%Ebar, i, w, eb, deb)
            end if
        else
            je = ds87_J(Z, -1d0, c_tau_e)
            if (full) lt = ds87_L(Z, c_tau_e)
            j1 = ds87_J(Z, 1d0, c_tau_ion(1))
            if (c_nion(3) > 0d0) j3 = ds87_J(Z, c_zion(3), c_tau_ion(3))
            dlt = hot(rtg(c_bin)%delta, i, w)
            ptr = hot(rtg(c_bin)%Ptr, i, w)
            if (full) then
                fc = hot(rtg(c_bin)%fcool, i, w)
                eb = hot(rtg(c_bin)%Ebar, i, w)
            end if
        end if
        ion = 0d0                                   ! ion_capture(Z)
        if (c_nion(1) > 0d0) ion = ion + c_arr_ion(1)*j1
        if (c_nion(2) > 0d0) ion = ion + c_arr_ion(2)*j1
        if (c_nion(3) > 0d0) ion = ion + c_arr_ion(3)*j3
        st = sticking(Z)
        kT = RTG_KB*c_T
        b_up(k) = photo + ion + c_arr_e*je*dlt
        b_ph(k) = photo                             ! the photoemission part (rtgroups_solve_bin sens)
        b_dn(k) = c_arr_e*je*st*(1d0 - ptr)
        if (full) then
            b_gam(k) = RTG_EV2ERG*heat
            b_lam(k) = c_arr_e*st*lt*kT*fc - c_arr_e*je*dlt*eb*RTG_EV2ERG
        end if
        b_A(k) = c_geoA*j1
        if (.not. want_d) return
        ! an ion with no density adds an exact zero: skipped
        dion = 0d0
        if (c_arr_ion(1) /= 0d0) dion = dion + c_arr_ion(1)*j1*(0.5d0 + dj1)
        if (c_arr_ion(2) /= 0d0) dion = dion + c_arr_ion(2)*j1*(0.5d0 + dj1)
        if (c_arr_ion(3) /= 0d0) dion = dion + c_arr_ion(3)*j3*(0.5d0 + dj3)
        b_dup(k) = dion + c_arr_e*je*(dlt*(0.5d0 + dje) + dde)
        b_ddn(k) = c_arr_e*je*st*((1d0 - ptr)*(0.5d0 + dje) - dpt)
        if (full) then
            b_lamd(k) = c_arr_e*(st*lt*kT*fc - je*dlt*eb*RTG_EV2ERG)
            b_dlam(k) = c_arr_e*(st*lt*kT*(fc*(1.5d0 + dld) + dfc) &
                        - je*RTG_EV2ERG*(dlt*eb*(0.5d0 + dje) + dde*eb + dlt*deb))
        end if
        b_dA(k) = 0.5d0 + dj1                       ! d ln b_A / d ln T (rtgroups_solve_bin sens)
    end subroutine eval_int

    subroutine window_P(lo, hi, pmax)
        ! P on lo..hi by running products of Up(Z)/Down(Z+1), relative to its maximum pmax
        implicit none
        integer, intent(in) :: lo, hi
        real(dp), intent(out) :: pmax
        integer :: k
        b_P(lo - b_base) = 1d0
        do k = lo - b_base, hi - b_base - 1
            b_P(k+1) = b_P(k)*(max(b_up(k), RTG_TINY)/max(b_dn(k+1), RTG_TINY))
            if (b_P(k+1) > 1d250) b_P(lo - b_base:k+1) = b_P(lo - b_base:k+1)*1d-250
        end do
        pmax = maxval(b_P(lo - b_base:hi - b_base))
    end subroutine window_P

    subroutine window(z0, z1, want_d, lo, hi, capped)
        ! pyCALIMA _window: exact P(Z) on z0..z1, grown or trimmed at each end to one charge
        ! past where P < DISCRETE_TAIL of its maximum; P normalised in b_P(lo..hi).
        ! capped: the growth stopped at MAX_DISCRETE_STATES with P still above the tail at an end
        ! that could grow, so lo..hi need not hold the peak (a warm window far from the new root)
        implicit none
        real(dp), intent(in) :: z0, z1
        logical, intent(in) :: want_d
        integer, intent(out) :: lo, hi
        logical, intent(out), optional :: capped
        integer :: k, k0, k1, zl, zh
        real(dp) :: pmax
        logical :: room, glo, ghi
        if (present(capped)) capped = .false.
        lo = nint(z0)
        hi = nint(z1)
        zl = nint(rtg(c_bin)%Zmin)
        zh = nint(rtg(c_bin)%Zmax)
        b_base = lo - MAX_DISCRETE_STATES - 2
        do k = lo, hi
            call eval_int(dble(k), want_d)
        end do
        do
            call window_P(lo, hi, pmax)
            room = hi - lo + 1 < MAX_DISCRETE_STATES
            glo = room .and. b_P(lo - b_base) > DISCRETE_TAIL*pmax .and. lo > zl
            ghi = room .and. b_P(hi - b_base) > DISCRETE_TAIL*pmax .and. hi < zh
            if (.not. (glo .or. ghi)) then
                if (present(capped)) capped = .not. room .and. &
                    ((b_P(lo - b_base) > DISCRETE_TAIL*pmax .and. lo > zl) .or. &
                     (b_P(hi - b_base) > DISCRETE_TAIL*pmax .and. hi < zh))
                exit
            end if
            if (glo) then
                lo = lo - 1
                call eval_int(dble(lo), want_d)
            end if
            if (ghi) then
                hi = hi + 1
                call eval_int(dble(hi), want_d)
            end if
        end do
        k0 = hi
        k1 = lo
        do k = lo, hi
            if (b_P(k - b_base) > DISCRETE_TAIL*pmax) then
                k0 = min(k0, k)
                k1 = max(k1, k)
            end if
        end do
        lo = max(k0 - 1, lo)
        hi = min(k1 + 1, hi)
        b_P(lo - b_base:hi - b_base) = b_P(lo - b_base:hi - b_base) / pmax
        b_P(lo - b_base:hi - b_base) = b_P(lo - b_base:hi - b_base) / sum(b_P(lo - b_base:hi - b_base))
    end subroutine window

    subroutine window_moments(lo, hi, zmean, zsig)
        implicit none
        integer, intent(in) :: lo, hi
        real(dp), intent(out) :: zmean, zsig
        real(dp) :: s2
        integer :: k
        zmean = 0d0
        s2 = 0d0
        do k = lo, hi
            zmean = zmean + b_P(k - b_base)*dble(k)
            s2 = s2 + b_P(k - b_base)*dble(k)**2
        end do
        zsig = sqrt(max(s2 - zmean**2, 0d0))
    end subroutine window_moments

    ! --------------------------------------------- recombination of a wide P(Z)
    real(dp) function chebval(x, c, n)
        ! numpy.polynomial.chebyshev.chebval (Clenshaw)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: x, c(n)
        real(dp) :: c0, c1, tmp, x2
        integer :: i
        if (n == 1) then
            c0 = c(1)
            c1 = 0d0
        else if (n == 2) then
            c0 = c(1)
            c1 = c(2)
        else
            x2 = 2d0*x
            c0 = c(n-1)
            c1 = c(n)
            do i = 3, n
                tmp = c0
                c0 = c(n - i + 1) - c1
                c1 = tmp + c1*x2
            end do
        end if
        chebval = c0 + c1*x
    end function chebval

    subroutine cheb_coef(fkind, a, b, n, c)
        ! Chebyshev coefficients of the degree n-1 interpolant on [a, b] of g (fkind 1) or ln A (2)
        implicit none
        integer, intent(in) :: fkind, n
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: c(n)
        real(dp) :: y(n), x
        integer :: j, k
        do k = 0, n - 1
            x = cos(RTG_PI*(dble(k) + 0.5d0)/dble(n))
            if (fkind == 1) then
                y(k+1) = gfun(0.5d0*(a + b) + 0.5d0*(b - a)*x)
            else
                y(k+1) = ln_capture_z1(0.5d0*(a + b) + 0.5d0*(b - a)*x)
            end if
        end do
        do j = 0, n - 1
            c(j+1) = 0d0
            do k = 0, n - 1
                c(j+1) = c(j+1) + y(k+1)*cos(RTG_PI*dble(j)*(dble(k) + 0.5d0)/dble(n))
            end do
            c(j+1) = 2d0/dble(n)*c(j+1)
        end do
        c(1) = c(1)*0.5d0
    end subroutine cheb_coef

    subroutine add_piece(p, q)
        ! a piece of g on the integers [p, q]: exact, or Chebyshev fits of g and ln A
        implicit none
        real(dp), intent(in) :: p, q
        real(dp) :: ints, cw(6), tmp(7), cd(5)
        integer :: nseg, nA, j, n1, k
        np_ = np_ + 1
        pc_p(np_) = p
        pc_q(np_) = q
        ints = q - p + 1d0
        nseg = int(min(dble(RECOMB_FIT_N), max(3d0, dble(ceiling(ints/RECOMB_PER_NODE)))))
        if ((-dble(RECOMB_KINK) <= p .and. q <= dble(RECOMB_KINK)) .or. ints <= dble(nseg)) then
            pc_kind(np_) = 1
            do k = 0, nint(ints) - 1
                pc_vals(k+1, np_) = gfun(p + dble(k))
            end do
            return
        end if
        pc_kind(np_) = 2
        pc_nc(np_) = nseg
        call cheb_coef(1, p, q, nseg, pc_c(1:nseg, np_))
        nA = min(RECOMB_FIT_NA, nseg)
        pc_ncA(np_) = nA
        call cheb_coef(2, p, q, nA, pc_cA(1:nA, np_))
        ! chebint (lbnd 0, k 0)
        cw(1:nseg) = pc_c(1:nseg, np_)
        tmp(1) = cw(1)*0d0
        tmp(2) = cw(1)
        if (nseg > 1) tmp(3) = cw(2)/4d0
        do j = 2, nseg - 1
            tmp(j+2) = cw(j+1)/dble(2*(j + 1))
            tmp(j) = tmp(j) - cw(j+1)/dble(2*(j - 1))
        end do
        tmp(1) = tmp(1) + (0d0 - chebval(0d0, tmp(1:nseg+1), nseg + 1))
        pc_ci(1:nseg+1, np_) = tmp(1:nseg+1)
        ! chebder
        cw(1:nseg) = pc_c(1:nseg, np_)
        n1 = nseg - 1
        do j = n1, 3, -1
            cd(j) = dble(2*j)*cw(j+1)
            cw(j-1) = cw(j-1) + (dble(j)*cw(j+1))/dble(j - 2)
        end do
        if (n1 > 1) cd(2) = 4d0*cw(3)
        cd(1) = cw(2)
        pc_cd(1:n1, np_) = cd(1:n1)
    end subroutine add_piece

    real(dp) function piece_t(k, z)
        implicit none
        integer, intent(in) :: k
        real(dp), intent(in) :: z
        piece_t = (2d0*z - pc_p(k) - pc_q(k)) / (pc_q(k) - pc_p(k))
    end function piece_t

    real(dp) function piece_dg(k, z)
        implicit none
        integer, intent(in) :: k
        real(dp), intent(in) :: z
        piece_dg = 2d0/(pc_q(k) - pc_p(k))*chebval(piece_t(k, z), pc_cd(1:pc_nc(k)-1, k), pc_nc(k) - 1)
    end function piece_dg

    real(dp) function piece_lnA(k, z)
        implicit none
        integer, intent(in) :: k
        real(dp), intent(in) :: z
        piece_lnA = chebval(piece_t(k, z), pc_cA(1:pc_ncA(k), k), pc_ncA(k))
    end function piece_lnA

    real(dp) function piece_gsum(k, m)
        ! sum_{j=p}^{m-1} g(j), p <= m <= q + 1 (Euler-Maclaurin for fits)
        implicit none
        integer, intent(in) :: k
        real(dp), intent(in) :: m
        real(dp) :: a, b, Ga, Gb
        integer :: j
        if (pc_kind(k) == 1) then
            piece_gsum = 0d0
            do j = 1, nint(m - pc_p(k))
                piece_gsum = piece_gsum + pc_vals(j, k)
            end do
            return
        end if
        if (m == pc_p(k)) then
            piece_gsum = 0d0
            return
        end if
        a = pc_p(k) - 0.5d0
        b = m - 0.5d0
        Ga = 0.5d0*(pc_q(k) - pc_p(k))*chebval(piece_t(k, a), pc_ci(1:pc_nc(k)+1, k), pc_nc(k) + 1)
        Gb = 0.5d0*(pc_q(k) - pc_p(k))*chebval(piece_t(k, b), pc_ci(1:pc_nc(k)+1, k), pc_nc(k) + 1)
        piece_gsum = Gb - Ga - (piece_dg(k, b) - piece_dg(k, a))/24d0
    end function piece_gsum

    subroutine cover(a, b, brk, nbrk)
        ! pieces for the integers [a, b): cuts at the break charges and at -K, K+1
        implicit none
        real(dp), intent(in) :: a, b, brk(:)
        integer, intent(in) :: nbrk
        real(dp) :: cuts(MAXB + 4)
        integer :: nc, j
        nc = 0
        call add_cut(a)
        call add_cut(b)
        do j = 1, nbrk
            if (a < brk(j) .and. brk(j) < b) call add_cut(brk(j))
        end do
        if (a < -dble(RECOMB_KINK) .and. -dble(RECOMB_KINK) < b) call add_cut(-dble(RECOMB_KINK))
        if (a < dble(RECOMB_KINK + 1) .and. dble(RECOMB_KINK + 1) < b) call add_cut(dble(RECOMB_KINK + 1))
        do j = 1, nc - 1
            call add_piece(cuts(j), cuts(j+1) - 1d0)
        end do
    contains
        subroutine add_cut(z)
            ! sorted, without duplicates
            real(dp), intent(in) :: z
            integer :: i, l
            do i = 1, nc
                if (cuts(i) == z) return
            end do
            l = nc + 1
            do i = 1, nc
                if (cuts(i) > z) then
                    l = i
                    exit
                end if
            end do
            cuts(l+1:nc+1) = cuts(l:nc)
            cuts(l) = z
            nc = nc + 1
        end subroutine add_cut
    end subroutine cover

    subroutine break_charges(brk, nbrk)
        ! pyCALIMA _break_charges: kinks of g in this cell, sorted, without duplicates
        implicit none
        real(dp), intent(out) :: brk(MAXB)
        integer, intent(out) :: nbrk
        integer :: g, k, j, key, nU, j0
        real(dp) :: rmax, zc, frac, x
        nbrk = 0
        do g = 1, rtg(c_bin)%nG
            if (c_cN(g) > 0d0) then
                if (rtg(c_bin)%edge0(g) /= NONE) call add_brk(rtg(c_bin)%edge0(g))
                if (rtg(c_bin)%edge1(g) /= NONE) call add_brk(rtg(c_bin)%edge1(g))
            end if
        end do
        if (c_G0 > 0d0) then
            do k = 1, 4
                if (rtg(c_bin)%edgebg(k) /= NONE) call add_brk(rtg(c_bin)%edgebg(k))
            end do
        end if
        nU = rtg(c_bin)%nU
        j0 = rtg(c_bin)%jpos
        do key = 1, 2
            if (key == 1) then
                f_row(j0:nU) = (1d0 - c_vh)*rtg(c_bin)%delta(j0:nU,c_kh) + c_vh*rtg(c_bin)%delta(j0:nU,c_kh+1)
            else
                f_row(j0:nU) = (1d0 - c_vh)*rtg(c_bin)%Ptr(j0:nU,c_kh) + c_vh*rtg(c_bin)%Ptr(j0:nU,c_kh+1)
            end if
            if (j0 > nU) cycle
            rmax = maxval(f_row(j0:nU))
            if (rmax <= 0d0) cycle
            do k = 1, 2
                frac = merge(1d-4, 0.5d0, k == 1)
                do j = j0, nU
                    if (f_row(j) > frac*rmax) then
                        x = rtg(c_bin)%U(j)*rtg(c_bin)%a/RTG_E2EVCM
                        zc = anint(x)
                        if (abs(x - aint(x)) == 0.5d0) zc = 2d0*anint(0.5d0*x)    ! numpy round: half to even
                        if (rtg(c_bin)%Zmin < zc .and. zc < rtg(c_bin)%Zmax) call add_brk(zc)
                        exit
                    end if
                end do
            end do
        end do
    contains
        subroutine add_brk(z)
            real(dp), intent(in) :: z
            integer :: i, l
            do i = 1, nbrk
                if (brk(i) == z) return
            end do
            l = nbrk + 1
            do i = 1, nbrk
                if (brk(i) > z) then
                    l = i
                    exit
                end if
            end do
            brk(l+1:nbrk+1) = brk(l:nbrk)
            brk(l) = z
            nbrk = nbrk + 1
        end subroutine add_brk
    end subroutine break_charges

    real(dp) function geometric_sum(n, Zn, lnf, zmax)
        ! sum over the integers Zn(1) .. min(zmax, Zn(n)) of exp(lnf), ln f linear between nodes
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: Zn(n), lnf(n), zmax
        real(dp) :: z0, z1, d, m, lr
        integer :: k
        geometric_sum = 0d0
        do k = 1, n - 1
            z0 = Zn(k)
            z1 = Zn(k+1)
            if (z0 > zmax) exit
            d = z1 - z0
            m = min(d, zmax - z0 + 1d0)
            lr = (lnf(k+1) - lnf(k)) / d
            if (abs(lr) < 1d-12) then
                geometric_sum = geometric_sum + m*exp(lnf(k))
            else if (lr < 0d0) then
                geometric_sum = geometric_sum + exp(lnf(k))*expm1(lr*m)/expm1(lr)
            else
                geometric_sum = geometric_sum + exp(lnf(k) + lr*(m - 1d0))*expm1(-lr*m)/expm1(-lr)
            end if
        end do
        if (Zn(n) <= zmax) geometric_sum = geometric_sum + exp(lnf(n))
    end function geometric_sum

    subroutine geometric_sums(n, Zn, lnf, zth, nz, out)
        ! geometric_sum for each zmax in zth(1:nz) in one pass: the full segments are summed once,
        ! in the order the single sum adds them, so each result is bitwise the single sum
        implicit none
        integer, intent(in) :: n, nz
        real(dp), intent(in) :: Zn(n), lnf(n), zth(nz)
        real(dp), intent(out) :: out(nz)
        integer :: k, j, kj
        real(dp) :: d, m
        f_ps(0) = 0d0
        do k = 1, n - 1
            f_ps(k) = f_ps(k-1) + seg(k, Zn(k+1) - Zn(k))
        end do
        do j = 1, nz
            kj = 0
            do k = 1, n - 1
                if (Zn(k) > zth(j)) exit
                kj = k
            end do
            if (kj == 0) then
                out(j) = 0d0
            else
                d = Zn(kj+1) - Zn(kj)
                m = min(d, zth(j) - Zn(kj) + 1d0)
                if (m == d) then
                    out(j) = f_ps(kj)
                else
                    out(j) = f_ps(kj-1) + seg(kj, m)
                end if
            end if
            if (Zn(n) <= zth(j)) out(j) = out(j) + exp(lnf(n))
        end do
    contains
        real(dp) function seg(k, mm)
            integer, intent(in) :: k
            real(dp), intent(in) :: mm
            real(dp) :: lr
            lr = (lnf(k+1) - lnf(k)) / (Zn(k+1) - Zn(k))
            if (abs(lr) < 1d-12) then
                seg = mm*exp(lnf(k))
            else if (lr < 0d0) then
                seg = exp(lnf(k))*expm1(lr*mm)/expm1(lr)
            else
                seg = exp(lnf(k) + lr*(mm - 1d0))*expm1(-lr*mm)/expm1(-lr)
            end if
        end function seg
    end subroutine geometric_sums

    subroutine recomb_fit(Zs, sigma, alpha)
        ! pyCALIMA GroupChargingTable._recombination_fit: grain-assisted recombination of a wide P(Z)
        implicit none
        real(dp), intent(in) :: Zs, sigma
        real(dp), intent(out) :: alpha(RTG_NRI)
        real(dp) :: brk(MAXB), Zlo, Zhi, top, bot, newb, rmin, lnA_geo, lnA_lo, base, z, la, curv, hz, below, vm, norm
        real(dp) :: lz(3), lav(3), lmax
        integer :: nbrk, it, j, k, kk, nn, nloc, i
        logical :: done, any_m

        Zlo = rtg(c_bin)%Zmin
        Zhi = rtg(c_bin)%Zmax - 1d0
        top = min(Zhi, dble(floor(Zs + RECOMB_UP*sigma)))
        bot = max(Zlo, dble(ceiling(Zs - RECOMB_DOWN*sigma)))
        call break_charges(brk, nbrk)
        np_ = 0
        call cover(bot, top + 1d0, brk, nbrk)
        rmin = huge(1d0)
        do j = 1, RTG_NRI
            z = rtg(c_bin)%Zth(j)
            if (z >= Zlo .and. z < bot .and. Zs - z < RECOMB_REACH*sigma) rmin = min(rmin, z)
        end do
        if (rmin < huge(1d0)) then
            newb = max(Zlo, dble(floor(rmin - 2d0*sigma)))
            call cover(newb, bot, brk, nbrk)
            bot = newb
        end if
        lnA_geo = log(RTG_PI*rtg(c_bin)%a**2*sqrt(8d0*RTG_KB*c_T/(RTG_PI*RTG_MU)))
        lnA_lo = ln_capture_z1(Zlo)
        nn = 0
        do it = 1, 40
            call sort_pieces()
            nn = 0
            base = 0d0
            do kk = 1, np_
                k = pc_order(kk)
                if (pc_kind(k) == 1) then
                    do j = 0, nint(pc_q(k) - pc_p(k))
                        z = pc_p(k) + dble(j)
                        nn = nn + 1
                        f_Zn(nn) = z
                        f_lnP(nn) = base + piece_gsum(k, z)
                        f_lnA(nn) = ln_capture_z1(z)
                    end do
                else
                    z = pc_p(k)
                    nloc = 0
                    do
                        la = piece_lnA(k, z)
                        nn = nn + 1
                        f_Zn(nn) = z
                        f_lnP(nn) = base + piece_gsum(k, z)
                        f_lnA(nn) = la
                        if (nloc == 3) then
                            lz(1:2) = lz(2:3)
                            lav(1:2) = lav(2:3)
                        else
                            nloc = nloc + 1
                        end if
                        lz(nloc) = z
                        lav(nloc) = la
                        if (z >= pc_q(k)) exit
                        curv = abs(piece_dg(k, z))
                        if (nloc >= 3) then
                            curv = curv + abs(2d0*((lav(3) - lav(2))/(lz(3) - lz(2)) &
                                                   - (lav(2) - lav(1))/(lz(2) - lz(1)))/(lz(3) - lz(1)))
                        else
                            curv = curv + 1d0
                        end if
                        hz = max(1d0, min(dble(floor(sqrt(8d0*RECOMB_EPS/max(curv, 1d-12)))), &
                                          max(1d0, dble(floor(sigma)))))
                        z = min(z + hz, pc_q(k))
                    end do
                end if
                base = base + piece_gsum(k, pc_q(k) + 1d0)
            end do
            lmax = maxval(f_lnP(1:nn))
            f_lnP(1:nn) = f_lnP(1:nn) - lmax
            f_v(1:nn) = f_lnP(1:nn) + f_lnA(1:nn)
            if (bot > Zlo .and. bot > 0d0 .and. f_lnP(1) + max(f_lnA(1) - lnA_geo, 0d0) < RECOMB_FLOOR) exit
            done = bot <= Zlo .or. f_lnP(1) < -RECOMB_CUT
            below = f_lnP(1) + lnA_lo + log(bot - Zlo + 1d0)
            do j = 1, RTG_NRI
                any_m = .false.
                vm = -huge(1d0)
                do i = 1, nn
                    if (f_Zn(i) <= rtg(c_bin)%Zth(j)) then
                        any_m = .true.
                        vm = max(vm, f_v(i))
                    end if
                end do
                if (done .and. any_m .and. f_v(1) > vm - RECOMB_CUT .and. max(vm, below) - lnA_geo > RECOMB_FLOOR) &
                    done = .false.
            end do
            if (done) exit
            newb = max(Zlo, dble(floor(bot - max(4d0*sigma, 20d0))))
            call cover(newb, bot, brk, nbrk)
            bot = newb
        end do
        norm = geometric_sum(nn, f_Zn(1:nn), f_lnP(1:nn), huge(1d0))
        call geometric_sums(nn, f_Zn(1:nn), f_v(1:nn), rtg(c_bin)%Zth, RTG_NRI, alpha)
        do j = 1, RTG_NRI
            alpha(j) = alpha(j) / norm / sqrt(RI_MASS(j))
        end do
    contains
        subroutine sort_pieces()
            ! pc_order: pieces by increasing p (stable)
            integer :: a1, b1, t
            do a1 = 1, np_
                pc_order(a1) = a1
            end do
            do a1 = 2, np_
                t = pc_order(a1)
                b1 = a1 - 1
                do while (b1 >= 1)
                    if (pc_p(pc_order(b1)) <= pc_p(t)) exit
                    pc_order(b1+1) = pc_order(b1)
                    b1 = b1 - 1
                end do
                pc_order(b1+1) = t
            end do
        end subroutine sort_pieces
    end subroutine recomb_fit

    real(dp) function log_ndtr(x)
        ! ln of the standard normal CDF (pyCALIMA _log_ndtr)
        implicit none
        real(dp), intent(in) :: x
        if (x < -5d0) then
            log_ndtr = log(0.5d0*erfc_scaled(-x/sqrt(2d0))) - 0.5d0*x*x
        else
            log_ndtr = log(0.5d0*erfc(-x/sqrt(2d0)))
        end if
    end function log_ndtr

    subroutine recomb_gauss(Zs, sigma, lg)
        ! pyCALIMA _recombination_gauss: ln of a Gaussian estimate of alpha_i, for the reuse of a fit
        implicit none
        real(dp), intent(in) :: Zs, sigma
        real(dp), intent(out) :: lg(RTG_NRI)
        real(dp) :: Zlo, Zhi, mu, zp, z0, sl, mup, base
        integer :: it, j
        Zlo = rtg(c_bin)%Zmin
        Zhi = rtg(c_bin)%Zmax - 1d0
        mu = Zs + 0.5d0
        zp = mu
        do it = 1, 2
            z0 = min(max(zp, Zlo + 1d0), Zhi - 1d0)
            sl = 0.5d0*(ln_capture_z1(z0 + 1d0) - ln_capture_z1(z0 - 1d0))
            zp = mu + sl*sigma**2
        end do
        zp = min(max(zp, Zlo), Zhi)
        mup = mu + sl*sigma**2
        base = ln_capture_z1(zp) + sl*(mu - zp) + 0.5d0*sl*sl*sigma*sigma
        do j = 1, RTG_NRI
            lg(j) = base + log_ndtr((rtg(c_bin)%Zth(j) + 0.5d0 - mup)/sigma) - ri_hlnm(j)
        end do
    end subroutine recomb_gauss


    ! ---------------------------------------------------------------- solver
    subroutine rtgroups_solve_bin(ii, Np, egy, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, Z_guess, &
                                  st, dlnT, alpha_mode, r, sens)
        ! Equilibrium charge and rates [erg/s per grain] of bin ii for one cell
        ! (pyCALIMA GroupChargingTable.solve). Np: photon densities per group
        ! [cm^-3]; egy: group mean energies [eV]; c_red: RT speed of light [cm/s];
        ! G0_bg: UV background [Habing]; Z_guess: start of a cold root search,
        ! ignored outside (Zmin, Zmax-1) (-huge(1d0) for none); st: warm state of
        ! this bin in this cell (in: the previous solve, out: this one); dlnT > 0:
        ! also d/d ln T; alpha_mode: grain-assisted recombination of a wide P(Z), RTG_ALPHA_NONE or
        ! _GAUSS (the estimate, which rtgroups_refine_alpha refines where it matters); a discrete P(Z) gets
        ! the exact sum unless _NONE.
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), egy(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, Z_guess, dlnT
        type(RTGState), intent(inout) :: st
        integer, intent(in) :: alpha_mode
        type(RTGResult), intent(out) :: r
        real(dp) :: Zs, sigma, pair(4), h, s0, s1, Zs1, sig1, chord(4), lnP, lnAlo
        real(dp) :: gm1, gp1, gam1, lam1, c(1 - MAX_DISCRETE_STATES:NBUF), cb, dvar, je, dje, EA, arrive2, zr, lg(RTG_NRI)
        integer :: lo, hi, k, j, half, Zc, lt
        logical :: disc, got, want_d, have_chord, capped
        logical, intent(in), optional :: sens    ! also d/d ln S, d/d ln N and the anchor of rtgroups_predict
        ! c, cS, cN: indexed by charge - lo + 1, from the lowest charge below the window in the recombination
        real(dp) :: cS(1 - MAX_DISCRETE_STATES:NBUF), cN(1 - MAX_DISCRETE_STATES:NBUF)
        real(dp) :: cbS, cbN, dvS, dvN, sg(2), sgm(2), sgp(2), gl, ll, gh, lh, dZS, dZN, sw

        if (any(rtg(ii)%egy_used /= egy(1:rtg(ii)%nG))) call rtgroups_set_group_energies(ii, egy)
        call set_cell(ii, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg)
        want_d = dlnT > 0d0
        r%alpha = 0d0

        ! 1. charge distribution: a discrete window reused from the warm state, or root + width
        got = .false.
        if (st%valid .and. st%discrete) then
            call window(st%wlo, st%whi, want_d, lo, hi, capped)
            call window_moments(lo, hi, r%Zmean, r%Zsigma)
            ! a capped window is a cut tail, narrow whatever the true width: solve from the root
            if (r%Zsigma < SIGMA_DISCRETE .and. .not. capped) then
                got = .true.
                ! Z*: root of g interpolated between the integers of the window
                Zs = dble(hi)
                if (hi == lo) then
                    Zs = dble(lo)
                else if (gwin(lo) <= 0d0) then
                    Zs = dble(lo)
                end if
                do k = lo, hi - 2
                    if (gwin(k) > 0d0 .and. gwin(k + 1) <= 0d0) then
                        Zs = dble(k) + gwin(k)/(gwin(k) - gwin(k + 1))
                        exit
                    end if
                end do
            end if
        end if
        if (got) then
            disc = .true.
        else
            call find_root(Z_guess, st, Zs, sigma, pair)
            have_chord = .false.
            if (sigma >= SIGMA_DISCRETE) then
                ! warm and wide: start the chord from the last width, usually one chord
                if (st%valid .and. .not. st%discrete) sigma = st%Zsigma
                call chord_sigma(Zs, sigma, chord, have_chord)
            end if
            disc = sigma < SIGMA_DISCRETE
            if (disc) then
                half = ceiling(3d0*max(sigma, 0.5d0))
                Zc = floor(Zs + 0.5d0)
                call window(max(rtg(ii)%Zmin, dble(Zc - half)), min(rtg(ii)%Zmax, dble(Zc + half)), want_d, lo, hi)
                call window_moments(lo, hi, r%Zmean, r%Zsigma)
            end if
        end if
        r%discrete = disc
        r%Zstar = Zs

        ! 2. heating, cooling, autoionisation; recombination; d/d ln T
        r%Gamma = 0d0
        r%Lambda = 0d0
        r%Lambda_auto = 0d0
        if (disc) then
            do k = lo, hi
                r%Gamma = r%Gamma + b_P(k - b_base)*b_gam(k - b_base)
                r%Lambda = r%Lambda + b_P(k - b_base)*b_lam(k - b_base)
            end do
            if (dble(lo) == rtg(ii)%Zmin) then
                EA = electron_affinity(ii, rtg(ii)%Zmin)
                arrive2 = RTG_PI*rtg(ii)%a**2*c_ne*c_ve
                r%Lambda_auto = b_P(lo - b_base)*arrive2*ds87_J(rtg(ii)%Zmin, -1d0, c_tau_e)*EA*RTG_EV2ERG
            end if
            lt = lo
            if (alpha_mode /= RTG_ALPHA_NONE) then
                do j = 1, RTG_NRI
                    do k = lo, hi
                        if (dble(k) <= rtg(ii)%Zth(j)) r%alpha(j) = r%alpha(j) + b_P(k - b_base)*b_A(k - b_base)
                    end do
                end do
                ! and the charges below the window (pyCALIMA _recombination_discrete): on a positive grain
                ! A rises by up to exp(e^2/akT) per charge towards lower Z, so P A can peak far below the
                ! window. P continues by g one charge at a time (in b_P, for sens) down to lt, stopping at
                ! Zmin or once P(Z) A(Zmin) (Z - Zmin), a bound on all the remaining terms, is e^-RECOMB_CUT
                ! below each ion's sum (at least 1e-300). No floor on the geometric rate: below it a
                ! stopped sum misses the charges that set it, and jumps wherever the stop moves
                lnP = log(b_P(lo - b_base))
                lnAlo = log(c_geoA*ds87_J(rtg(ii)%Zmin, 1d0, c_tau_ion(1)))
                do k = lo - 1, max(lo - MAX_DISCRETE_STATES, nint(rtg(ii)%Zmin), b_base + 1), -1
                    call eval_int(dble(k), want_d, rates_only=.true.)
                    lnP = lnP - gwin(k)
                    b_P(k - b_base) = exp(lnP)
                    do j = 1, RTG_NRI
                        if (dble(k) <= rtg(ii)%Zth(j)) r%alpha(j) = r%alpha(j) + b_P(k - b_base)*b_A(k - b_base)
                    end do
                    lt = k
                    if (all(lnP + lnAlo + log(max(dble(k) - rtg(ii)%Zmin, 1d0)) < &
                            log(max(r%alpha, 1d-300)) - RECOMB_CUT)) exit
                end do
                r%alpha = r%alpha/sqrt(RI_MASS)
            end if
            if (want_d) then
                c(1) = 0d0
                do k = lo, hi - 1
                    c(k - lo + 2) = c(k - lo + 1) + b_dup(k - b_base)/max(b_up(k - b_base), RTG_TINY) &
                                    - b_ddn(k + 1 - b_base)/max(b_dn(k + 1 - b_base), RTG_TINY)
                end do
                cb = 0d0
                do k = lo, hi
                    cb = cb + b_P(k - b_base)*c(k - lo + 1)
                end do
                r%dZmean = 0d0
                dvar = 0d0
                r%dGamma = 0d0
                r%dLambda = 0d0
                do k = lo, hi
                    c(k - lo + 1) = c(k - lo + 1) - cb
                    r%dZmean = r%dZmean + b_P(k - b_base)*(c(k - lo + 1)*dble(k))
                    dvar = dvar + b_P(k - b_base)*(c(k - lo + 1)*dble(k)*dble(k))
                    r%dGamma = r%dGamma + b_P(k - b_base)*(c(k - lo + 1)*b_gam(k - b_base))
                end do
                do k = lo, hi
                    r%dLambda = r%dLambda + b_P(k - b_base)*(c(k - lo + 1)*b_lamd(k - b_base))
                end do
                zr = 0d0
                do k = lo, hi
                    zr = zr + b_P(k - b_base)*b_dlam(k - b_base)
                end do
                r%dLambda = r%dLambda + zr
                if (dble(lo) == rtg(ii)%Zmin) then
                    je = ds87_J(dble(lo), -1d0, c_tau_e, dje)
                    r%dLambda = r%dLambda + r%Lambda_auto*(c(1) + 0.5d0 + dje)
                end if
                dvar = dvar - 2d0*r%Zmean*r%dZmean
                r%dZsigma = 0d0
                if (r%Zsigma > 0d0) r%dZsigma = 0.5d0*dvar/r%Zsigma
            end if
            st%wlo = dble(lo)
            st%whi = dble(hi)
        else
            ! a wide P(Z) is averaged at its mean Z* + 1/2 (pyCALIMA _wide_nodes)
            r%Zmean = Zs + 0.5d0
            r%Zsigma = sigma
            call heating_cooling(Zs + 0.5d0, r%Gamma, r%Lambda)
            if (alpha_mode == RTG_ALPHA_GAUSS) then
                call recomb_gauss(Zs, sigma, lg)
                r%alpha = exp(lg)
            end if
            if (want_d) then
                call set_cell_T(T*(1d0 + dlnT))
                gm1 = gfun(pair(1))
                gp1 = gfun(pair(2))
                s0 = (pair(4) - pair(3))/(pair(2) - pair(1))
                s1 = (gp1 - gm1)/(pair(2) - pair(1))
                Zs1 = Zs + (pair(1) - gm1/s1) - (pair(1) - pair(3)/s0)
                sig1 = sigma
                if (s1 < 0d0) sig1 = sqrt(-1d0/s1)
                if (have_chord) then
                    s1 = (gfun(chord(2)) - gfun(chord(1)))/(chord(2) - chord(1))
                    sig1 = sigma
                    if (s1 < 0d0) sig1 = sqrt(-1d0/s1)
                end if
                call heating_cooling(Zs1 + 0.5d0, gam1, lam1)
                call set_cell_T(T)
                h = log1p(dlnT)
                r%dZmean = (Zs1 - Zs)/h
                r%dZsigma = (sig1 - sigma)/h
                r%dGamma = (gam1 - r%Gamma)/h
                r%dLambda = (lam1 - r%Lambda)/h
            end if
        end if
        st%valid = .true.
        st%discrete = disc
        st%Zstar = Zs
        st%Zsigma = r%Zsigma

        ! 3. sens: d/d ln S and d/d ln N (pyCALIMA GroupChargingTable._sensitivities), and the anchor
        st%anchored = .false.
        if (.not. present(sens)) return
        if (.not. sens .or. rtg(ii)%nG > RTG_MAXG) return
        if (disc) then
            ! d ln P(Z): cumulative sums of d ln Up - d ln Down, as the d/d ln T above; heating ~ S,
            ! cooling and autoionisation ~ N
            cS(1) = 0d0
            cN(1) = 0d0
            do k = lo, hi - 1
                cb = max(b_up(k - b_base), RTG_TINY)
                cS(k - lo + 2) = cS(k - lo + 1) + b_ph(k - b_base)/cb
                cN(k - lo + 2) = cN(k - lo + 1) + ((cb - b_ph(k - b_base))/cb - 1d0)
            end do
            cbS = 0d0
            cbN = 0d0
            do k = lo, hi
                cbS = cbS + b_P(k - b_base)*cS(k - lo + 1)
                cbN = cbN + b_P(k - b_base)*cN(k - lo + 1)
            end do
            r%dZmean_S = 0d0
            r%dZmean_N = 0d0
            dvS = 0d0
            dvN = 0d0
            r%dGamma_S = 0d0
            r%dGamma_N = 0d0
            r%dLambda_S = 0d0
            r%dLambda_N = 0d0
            do k = lo, hi
                cS(k - lo + 1) = cS(k - lo + 1) - cbS
                cN(k - lo + 1) = cN(k - lo + 1) - cbN
                r%dZmean_S = r%dZmean_S + b_P(k - b_base)*(cS(k - lo + 1)*dble(k))
                r%dZmean_N = r%dZmean_N + b_P(k - b_base)*(cN(k - lo + 1)*dble(k))
                dvS = dvS + b_P(k - b_base)*(cS(k - lo + 1)*dble(k)*dble(k))
                dvN = dvN + b_P(k - b_base)*(cN(k - lo + 1)*dble(k)*dble(k))
                r%dGamma_S = r%dGamma_S + b_P(k - b_base)*(cS(k - lo + 1)*b_gam(k - b_base))
                r%dGamma_N = r%dGamma_N + b_P(k - b_base)*(cN(k - lo + 1)*b_gam(k - b_base))
                r%dLambda_S = r%dLambda_S + b_P(k - b_base)*(cS(k - lo + 1)*b_lam(k - b_base))
                r%dLambda_N = r%dLambda_N + b_P(k - b_base)*(cN(k - lo + 1)*b_lam(k - b_base))
            end do
            r%dGamma_S = r%dGamma_S + r%Gamma
            r%dLambda_N = r%dLambda_N + r%Lambda
            if (r%Lambda_auto /= 0d0) then
                r%dLambda_S = r%dLambda_S + r%Lambda_auto*cS(1)
                r%dLambda_N = r%dLambda_N + r%Lambda_auto*(cN(1) + 1d0)
            end if
            dvS = dvS - 2d0*r%Zmean*r%dZmean_S
            dvN = dvN - 2d0*r%Zmean*r%dZmean_N
            r%dZsigma_S = 0d0
            r%dZsigma_N = 0d0
            if (r%Zsigma > 0d0) then
                r%dZsigma_S = 0.5d0*dvS/r%Zsigma
                r%dZsigma_N = 0.5d0*dvN/r%Zsigma
            end if
            ! d ln alpha_i: sum P A over the charges at or below Zth_i, P moving as above and
            ! A ~ T^1/2 J~(Z, 1) also with T
            st%adla = 0d0
            if (alpha_mode /= RTG_ALPHA_NONE .and. want_d) then
                ! the charges below the window down to lt (step 2): d ln P continued downwards from lo
                do k = lo - 1, lt, -1
                    c(k - lo + 1) = c(k - lo + 2) - (b_dup(k - b_base)/max(b_up(k - b_base), RTG_TINY) &
                                                     - b_ddn(k + 1 - b_base)/max(b_dn(k + 1 - b_base), RTG_TINY))
                    cb = max(b_up(k - b_base), RTG_TINY)
                    cS(k - lo + 1) = cS(k - lo + 2) - b_ph(k - b_base)/cb
                    cN(k - lo + 1) = cN(k - lo + 2) - ((cb - b_ph(k - b_base))/cb - 1d0)
                end do
                do j = 1, RTG_NRI
                    sw = 0d0
                    do k = lt, hi
                        if (dble(k) > rtg(ii)%Zth(j)) cycle
                        cb = b_P(k - b_base)*b_A(k - b_base)
                        sw = sw + cb
                        st%adla(j, 1) = st%adla(j, 1) + cb*cS(k - lo + 1)
                        st%adla(j, 2) = st%adla(j, 2) + cb*cN(k - lo + 1)
                        st%adla(j, 3) = st%adla(j, 3) + cb*(c(k - lo + 1) + b_dA(k - b_base))
                    end do
                    if (sw > 0d0) st%adla(j, :) = st%adla(j, :)/sw
                end do
            end if
        else
            ! the root moves by sigma^2 dg/dtheta; sigma with the slope across Z* -+ 1/2; Gamma and
            ! Lambda (at Z* + 1/2) along their slope over [Z*, Z* + 1]
            call dg_dtheta(Zs, sg)
            call dg_dtheta(Zs - 0.5d0, sgm)
            call dg_dtheta(Zs + 0.5d0, sgp)
            call heating_cooling(Zs, gl, ll)
            call heating_cooling(Zs + 1d0, gh, lh)
            dZS = sigma*sigma*sg(1)
            dZN = sigma*sigma*sg(2)
            r%dZmean_S = dZS
            r%dZmean_N = dZN
            r%dZsigma_S = 0.5d0*sigma**3*(sgp(1) - sgm(1))
            r%dZsigma_N = 0.5d0*sigma**3*(sgp(2) - sgm(2))
            r%dGamma_S = (gh - gl)*dZS + r%Gamma
            r%dGamma_N = (gh - gl)*dZN
            r%dLambda_S = (lh - ll)*dZS
            r%dLambda_N = (lh - ll)*dZN + r%Lambda
            ! the change of the Gaussian estimate scales the predicted alpha (fitted or not)
            if (alpha_mode == RTG_ALPHA_GAUSS) then
                st%alg = lg
            else if (alpha_mode /= RTG_ALPHA_NONE) then
                call recomb_gauss(Zs, sigma, st%alg)
            end if
        end if
        call anchor_kernels(Zs, st)
        st%anchored = .true.
        st%aT = T
        st%ane = ne
        st%anion = (/n_Hp, n_Hep, n_Hepp/)
        st%acred = c_red
        st%aG0 = G0_bg
        st%aNp = 0d0
        st%aNp(1:rtg(ii)%nG) = Np(1:rtg(ii)%nG)
        st%ar = r
    contains
        subroutine dg_dtheta(z, d)
            ! dg/d ln S and dg/d ln N at charge z (Up at z, Down at z + 1 ~ n_e)
            real(dp), intent(in) :: z
            real(dp), intent(out) :: d(2)
            real(dp) :: ww, ph, ht, upz
            integer :: iz
            call u_of(z, iz, ww)
            call kernels(iz, ww, ph, ht, .false.)
            upz = max(ph + ion_capture(z) + c_arr_e*ds87_J(z, -1d0, c_tau_e)*hot(rtg(c_bin)%delta, iz, ww), RTG_TINY)
            d(1) = ph/upz
            d(2) = (upz - ph)/upz - 1d0
        end subroutine dg_dtheta
        real(dp) function gwin(k)
            ! g increment of the window from charge k to k + 1
            integer, intent(in) :: k
            gwin = log(max(b_up(k - b_base), RTG_TINY)) - log(max(b_dn(k + 1 - b_base), RTG_TINY))
        end function gwin
    end subroutine rtgroups_solve_bin

    subroutine rtgroups_refine_alpha(ii, Np, egy, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, r)
        ! The fitted recombination of a wide P(Z) (r from rtgroups_solve_bin of the same cell),
        ! for the bins dust_interface finds important (pyCALIMA refine_recombination)
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), egy(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg
        type(RTGResult), intent(inout) :: r
        if (any(rtg(ii)%egy_used /= egy(1:rtg(ii)%nG))) call rtgroups_set_group_energies(ii, egy)
        call set_cell(ii, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg)
        call recomb_fit(r%Zstar, r%Zsigma, r%alpha)
    end subroutine rtgroups_refine_alpha


    subroutine rtgroups_uniform_lookup(ii, G0, T, ne, Zmean, Zsigma, Gamma, Lambda)
        ! Cells without local RT photons: bilinear lookup in (log T, log ne) of
        ! the full solve at uniform background G0 (ions: H+ with n = ne); Lambda
        ! includes autoionisation. The table is (re)built when G0 changes.
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: G0, T, ne
        real(dp), intent(out) :: Zmean, Zsigma, Gamma, Lambda
        real(dp), allocatable :: Np0(:), egy0(:)
        real(dp) :: x, y, fx, fy, u, v, ne_n, dlt, dlne
        integer :: i, j
        type(RTGState) :: st0
        type(RTGResult) :: r0

        dlt = (UNI_LT1 - UNI_LT0) / dble(NUNI_T - 1)
        dlne = (UNI_LNE1 - UNI_LNE0) / dble(NUNI_NE - 1)
        if (rtg(ii)%G0_uni /= G0 .or. .not. allocated(rtg(ii)%uni)) then
            if (.not. allocated(rtg(ii)%uni)) allocate(rtg(ii)%uni(4,NUNI_T,NUNI_NE))
            allocate(Np0(rtg(ii)%nG), egy0(rtg(ii)%nG))
            Np0 = 0d0
            egy0 = rtg(ii)%egy_used
            do i = 1, NUNI_T
                do j = 1, NUNI_NE
                    ne_n = 10d0**(UNI_LNE0 + dble(j - 1)*dlne)
                    st0%valid = .false.
                    call rtgroups_solve_bin(ii, Np0, egy0, 0d0, 10d0**(UNI_LT0 + dble(i - 1)*dlt), ne_n, &
                                            ne_n, 0d0, 0d0, G0, -huge(1d0), st0, 0d0, RTG_ALPHA_NONE, r0)
                    rtg(ii)%uni(:,i,j) = (/ r0%Zmean, r0%Zsigma, r0%Gamma, r0%Lambda + r0%Lambda_auto /)
                end do
            end do
            deallocate(Np0, egy0)
            rtg(ii)%G0_uni = G0
        end if
        x = min(max(log10(T), UNI_LT0), UNI_LT1)
        y = min(max(log10(max(ne, RTG_TINY)), UNI_LNE0), UNI_LNE1)
        fx = (x - UNI_LT0) / dlt
        fy = (y - UNI_LNE0) / dlne
        i = min(int(fx), NUNI_T - 2)
        j = min(int(fy), NUNI_NE - 2)
        u = fx - dble(i)
        v = fy - dble(j)
        i = i + 1
        j = j + 1
        Zmean = (1d0-u)*(1d0-v)*rtg(ii)%uni(1,i,j) + u*(1d0-v)*rtg(ii)%uni(1,i+1,j) &
              + (1d0-u)*v*rtg(ii)%uni(1,i,j+1) + u*v*rtg(ii)%uni(1,i+1,j+1)
        Zsigma = (1d0-u)*(1d0-v)*rtg(ii)%uni(2,i,j) + u*(1d0-v)*rtg(ii)%uni(2,i+1,j) &
              + (1d0-u)*v*rtg(ii)%uni(2,i,j+1) + u*v*rtg(ii)%uni(2,i+1,j+1)
        Gamma = (1d0-u)*(1d0-v)*rtg(ii)%uni(3,i,j) + u*(1d0-v)*rtg(ii)%uni(3,i+1,j) &
              + (1d0-u)*v*rtg(ii)%uni(3,i,j+1) + u*v*rtg(ii)%uni(3,i+1,j+1)
        Lambda = (1d0-u)*(1d0-v)*rtg(ii)%uni(4,i,j) + u*(1d0-v)*rtg(ii)%uni(4,i+1,j) &
              + (1d0-u)*v*rtg(ii)%uni(4,i,j+1) + u*v*rtg(ii)%uni(4,i+1,j+1)
    end subroutine rtgroups_uniform_lookup

    subroutine rtgroups_coulomb_ratio(ii, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, Zstar, Zsigma, Zmean, ratio)
        ! The correction of the Gaussian Coulomb focusing factors of a narrow P(Z) of bin ii: D_j over
        ! this cell's exact P(Z), on every charge within 8 sigma + 10 of Zstar and from -3 to +3 (the
        ! tail at the opposite sign sets D_j of small grains in cold gas) down to 1e-300 of its
        ! maximum, over D_j of the Gaussian of Zmean and Zsigma (compute_Coulomb_focusing_ions), for
        ! the impactor charges ISRF_DCHARGE. 1 where that range does not fit in the window buffers.
        ! P(Z) is summed in ln (pyCALIMA charging_isrf_tables.full_distribution): across the tail of
        ! a small grain in cold gas one step of P can exceed the range of a real.
        use constants, only: pi, e2instatC, kB
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, Zstar, Zsigma, Zmean
        real(dp), dimension(ISRF_ND), intent(out) :: ratio
        real(dp), parameter :: LN_CUT = -690.7755278982137d0      ! ln 1e-300
        real(dp) :: pmax, C1, b0, g, s, Dg(-1:3), Dx(ISRF_ND)
        integer :: lo, hi, k, m, j

        ratio = 1d0
        lo = max(nint(rtg(ii)%Zmin), min(floor(Zstar - 8d0*max(Zsigma, 0.5d0) - 10d0), -3))
        hi = min(nint(rtg(ii)%Zmax), max(ceiling(Zstar + 8d0*max(Zsigma, 0.5d0) + 10d0), 3))
        if (hi - lo + 3 > NBUF) return
        call set_cell(ii, Np, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg)
        b_base = lo - 2
        do k = lo, hi
            call eval_int(dble(k), .false.)
        end do
        b_P(lo - b_base) = 0d0
        do k = lo - b_base, hi - b_base - 1
            b_P(k+1) = b_P(k) + log(max(b_up(k), RTG_TINY)) - log(max(b_dn(k+1), RTG_TINY))
        end do
        pmax = maxval(b_P(lo - b_base:hi - b_base))
        C1 = e2instatC / (kB * T * dustbins_props(ii)%asize_cm)
        s = 0d0
        Dx = 0d0
        do k = lo, hi
            g = b_P(k - b_base) - pmax
            if (g <= LN_CUT) cycle
            g = exp(g)
            s = s + g
            do m = 1, ISRF_ND
                j = ISRF_DCHARGE(m)
                if (k == 0) then
                    b0 = 1d0 + dble(abs(j)) * sqrt(pi * C1 / 2d0)
                    Dx(m) = Dx(m) + g * b0
                else if (k * j > 0) then
                    Dx(m) = Dx(m) + g * exp(max(-dble(k * j) * C1, -745d0))
                else
                    Dx(m) = Dx(m) + g * (1d0 - dble(k * j) * C1)
                end if
            end do
        end do
        Dg = 1d0
        call compute_Coulomb_focusing_ions(ii, T, Zmean, Zsigma, 3, Dg)
        do m = 1, ISRF_ND
            ratio(m) = (Dx(m) / s) / Dg(ISRF_DCHARGE(m))
        end do
    end subroutine rtgroups_coulomb_ratio

    ! ------------------------------------------------- charge outside the cooling step
    integer(kind=8) function rtgroups_table_bytes()
        ! bytes of this rank's WDB06rt and WDB06isrf tables and work arrays (the memory report);
        ! the uniform-background table of a bin is built on first use and counted as if it were
        use dust_utils, only: alloc_bytes
        implicit none
        integer :: ii
        integer(kind=8) :: b
        b = alloc_bytes(isrf_q) + alloc_bytes(c_cN) + alloc_bytes(f_Zn) + alloc_bytes(f_lnP) + alloc_bytes(f_lnA) &
          + alloc_bytes(f_v) + alloc_bytes(f_ps) + alloc_bytes(f_row)
        if (allocated(rtg)) then
            do ii = 1, size(rtg)
                associate (t => rtg(ii))
                    b = b + alloc_bytes(t%L0) + alloc_bytes(t%L1) + alloc_bytes(t%U) + alloc_bytes(t%lTs) &
                          + alloc_bytes(t%egy) + alloc_bytes(t%K) + alloc_bytes(t%H) + alloc_bytes(t%lTh) &
                          + alloc_bytes(t%delta) + alloc_bytes(t%Ptr) + alloc_bytes(t%fcool) + alloc_bytes(t%Ebar) &
                          + alloc_bytes(t%KgT) + alloc_bytes(t%HgT) + alloc_bytes(t%mKT) + alloc_bytes(t%mHT) &
                          + alloc_bytes(t%egy_used) + alloc_bytes(t%Kbg) + alloc_bytes(t%Hbg) + alloc_bytes(t%mKb) &
                          + alloc_bytes(t%mHb) + alloc_bytes(t%edge0) + alloc_bytes(t%edge1) + alloc_bytes(t%iuz) &
                          + alloc_bytes(t%wuz)
                    if (t%nU > 0) b = b + 8_8 * 4 * NUNI_T * NUNI_NE   ! uni(4, NUNI_T, NUNI_NE)
                end associate
            end do
        end if
        rtgroups_table_bytes = b
    end function rtgroups_table_bytes

    subroutine dust_charge_moments(ii, Np, egy, c_red, G0_bg, G0, T, ne, n_Hp, n_Hep, n_Hepp, Zmean, Zsigma, &
                                   lnL1, coul)
        ! <Z> and sigma_Z of dust bin ii from charging_model, for callers outside the RTZ
        ! cooling step (the Coulomb drag): the same models as compute_dust_precool, but
        ! stateless. WDB06rt: a cold per-group solve (no recombination, no T derivatives), or
        ! the uniform-background table without local photons; otherwise the uniform-ISRF
        ! tables (or the analytic fit) of dust_charging. Np: photon densities
        ! per group [cm^-3]; egy [eV]; c_red [cm/s]; G0_bg: UV background [Habing]; G0:
        ! background plus the local 5.6-13.6 eV groups [Habing] (the dust_charging models only).
        ! With lnL1 = ln Lambda of Z = 1, also coul = <Z^2 ln+(Lambda_1/|Z|)> over P(Z), the
        ! Draine & Salpeter (1979) Coulomb term in units of phi(Z=1)^2: on the exact P(Z) for
        ! WDB06rt (the solve's window if discrete, else one evaluated around Z* while it fits
        ! MAX_DISCRETE_STATES), else on a Gaussian of the model's <Z> and sigma_Z. coul needs lnL1.
#ifdef RT
        use rt_parameters, only: smallNp
#endif
        implicit none
        integer, intent(in) :: ii
        real(dp), intent(in) :: Np(:), egy(:), c_red, G0_bg, G0, T, ne, n_Hp, n_Hep, n_Hepp
        real(dp), intent(out) :: Zmean
        real(dp), intent(out), optional :: Zsigma
        real(dp), intent(in), optional :: lnL1
        real(dp), intent(out), optional :: coul
        real(dp) :: Gamma, Lambda, Zsig, Dj(ISRF_ND)
        integer :: lo, hi, half, Zc
        logical :: exact
        type(RTGState) :: st
        type(RTGResult) :: r
#ifndef RT
        ! rt_parameters' floor; without RT, WDB06rt is rejected by check_params_dust
        real(dp), parameter :: smallNp = 1d-30
#endif
        exact = .false.
        if (trim(charging_model) == 'WDB06rt') then
            if (all(Np <= smallNp)) then
                call rtgroups_uniform_lookup(ii, G0_bg, T, ne, Zmean, Zsig, Gamma, Lambda)
            else
                call rtgroups_solve_bin(ii, Np, egy, c_red, T, ne, n_Hp, n_Hep, n_Hepp, G0_bg, -huge(1d0), &
                                        st, 0d0, RTG_ALPHA_NONE, r)
                Zmean = r%Zmean
                Zsig = r%Zsigma
                if (present(coul)) then
                    if (r%discrete) then
                        lo = nint(st%wlo)
                        hi = nint(st%whi)
                        exact = .true.
                    else if (2*ceiling(5d0*r%Zsigma) + 1 <= MAX_DISCRETE_STATES) then
                        ! the cell state of the solve is still set
                        half = ceiling(3d0*r%Zsigma)
                        Zc = floor(r%Zstar + 0.5d0)
                        call window(max(rtg(ii)%Zmin, dble(Zc - half)), min(rtg(ii)%Zmax, dble(Zc + half)), &
                                    .false., lo, hi)
                        exact = hi - lo + 1 < MAX_DISCRETE_STATES     ! else truncated: the Gaussian
                    end if
                    if (exact) coul = coul_window(lo, hi, lnL1)
                end if
            end if
        else
            call isrf_lookup(ii, G0, T, ne, Zmean, Zsig, Gamma, Lambda, Dj)
        end if
        if (present(coul) .and. .not. exact) coul = coul_gauss(Zmean, Zsig, lnL1)
        if (present(Zsigma)) Zsigma = Zsig
    end subroutine dust_charge_moments

    real(dp) function coul_term(Z, lnL1)
        ! Z^2 ln+(Lambda_1/|Z|): the Draine & Salpeter (1979) Coulomb term of charge Z, in units of phi(Z=1)^2
        implicit none
        real(dp), intent(in) :: Z, lnL1
        coul_term = 0d0
        if (Z /= 0d0) coul_term = Z**2*max(lnL1 - log(abs(Z)), 0d0)
    end function coul_term

    real(dp) function coul_window(lo, hi, lnL1)
        ! coul_term over the exact P(Z) of the window lo..hi
        implicit none
        integer, intent(in) :: lo, hi
        real(dp), intent(in) :: lnL1
        integer :: k
        coul_window = 0d0
        do k = lo, hi
            coul_window = coul_window + b_P(k - b_base)*coul_term(dble(k), lnL1)
        end do
    end function coul_window

    real(dp) function coul_gauss(mu, sig, lnL1)
        ! coul_term over a Gaussian P(Z) of mean mu and width sig: 10-point Gauss-Hermite
        implicit none
        real(dp), intent(in) :: mu, sig, lnL1
        real(dp), parameter :: X(5) = (/0.3429013272237046d0, 1.0366108297895137d0, 1.7566836492998818d0, &
                                        2.5327316742327897d0, 3.4361591188377376d0/)
        real(dp), parameter :: W(5) = (/0.6108626337353258d0, 0.2401386110823147d0, 3.387439445548106d-2, &
                                        1.3436457467812327d-3, 7.640432855232621d-6/)
        integer :: k
        coul_gauss = 0d0
        do k = 1, 5
            coul_gauss = coul_gauss + W(k)*(coul_term(mu + sqrt(2d0)*sig*X(k), lnL1) &
                                          + coul_term(mu - sqrt(2d0)*sig*X(k), lnL1))
        end do
        coul_gauss = coul_gauss/sqrt(RTG_PI)
    end function coul_gauss

end module dust_charging