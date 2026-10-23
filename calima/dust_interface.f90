module dust_interface
    use amr_commons, only:dp,ndim
    use amr_parameters, only: nvector
    use dust_charging, only: RTGState, RTGResult, RTG_NRI, RTG_ALPHA_GAUSS, RTG_ALPHA_NONE, RTG_RECOMB_IMPORTANT, ISRF_ND, &
                                      RI_ATOMIC
    use dust_radiation, only: TdustLin
    use constants
    use dust_commons
    implicit none

    private

    public :: compute_dust_rad_rates,compute_dust_precool,&
            compute_dust_coolrates,compute_local_anisotropy_factor,&
            compute_dust_update,rtgroups_reset_warm,store_pah_ion_fraction,&
            rtgroups_derivative_calls,reset_cell_warm_starts,get_cell_zguess

    ! Cache to avoid automatic array allocations in compute_dust_update
    real(dp), allocatable, save, target :: y_gas_cache(:,:), y_gas_out_cache(:,:)
    real(dp), allocatable, save, target :: y_dust_cache(:), y_dust_out_cache(:)
    ! WDB06rt: warm state of every bin in every cell of the RTZ vector (sub-step to
    ! sub-step); the last full solve of each bin, whose d/d ln T serves the cooling
    ! solver's calls at T(1 +- 1e-5) with the other inputs unchanged; the last root
    ! of each bin, the starting guess of a cold solve
    type(RTGState), allocatable, target, save :: rtg_warm(:,:)     ! (nvector + 1, ndust); nvector + 1: no cell index
    type RTGLast
        logical :: valid = .false.
        real(dp) :: T = 0d0, ne = 0d0, n_Hp = 0d0, n_Hep = 0d0, n_Hepp = 0d0, G0 = 0d0, c_red = 0d0
        real(dp), allocatable :: Np(:), egy(:)
        type(RTGResult) :: r
    end type RTGLast
    type(RTGLast), allocatable, save :: rtg_last(:)
    real(dp), allocatable, save :: rtg_Zguess(:)
    ! compute_dust_update: the last accepted substep's recommended next-step size, the first guess of the
    ! next call (it skips the firstcall 1/kmax reduction at equilibrium for the explicit solvers)
    real(dp), save :: last_h_ode = 0.0_dp
    real(dp), parameter :: RTG_DLNT = 1d-5, RTG_SERVE = 1.5d-5
    real(dp), parameter :: RTG_SERVE_HI = exp(RTG_SERVE), RTG_SERVE_LO = exp(-RTG_SERVE)   ! rtg_served, without a log
    ! the dust temperatures of the last update_T_dust and the inputs it had: the cooling solver's
    ! T(1 +- 1e-5) calls with the same other inputs are served from it (serve_T_dust)
    type(TdustLin), allocatable, save :: td_lin(:)
    logical, save :: td_valid = .false., td_hasNp = .false.
    real(dp), save :: td_T = 0d0, td_ne = 0d0, td_G0 = 0d0, td_nH2 = 0d0, td_nCO = 0d0
    real(dp), allocatable, save :: td_nel(:), td_xion(:,:), td_Np(:), td_egy(:)
    ! The grain charge of the dust bins by charging model (compute_dust_precool): the charge (dinfo%Z_dust,
    ! Z_sigma) and, where the model gives them, the PE heating and recombination cooling per grain, the
    ! correction of the Coulomb factors of the low impactor charges ISRF_DCHARGE (Dmode 1: the Gaussian
    ! factors times D, the exact-P(Z) ratio; 2: D itself) and the grain-assisted recombination of the
    ! ions RI_ATOMIC per grain. The routines are chosen once from charging_model (charge_select).
    type ChargeOut
        logical :: has_pe = .false., has_alpha = .false.
        logical :: want_alpha = .false.   ! in: the call needs the grain-assisted recombination (alpha)
        real(dp) :: Pinj(1:ndust) = 0d0, Prec(1:ndust) = 0d0          ! erg/s per grain
        integer :: Dmode(1:ndust) = 0
        real(dp) :: D(ISRF_ND, 1:ndust) = 1d0
        real(dp) :: alpha(RTG_NRI, 1:ndust) = 0d0                    ! cm^3 s^-1 per grain
    end type ChargeOut
    abstract interface
        subroutine charge_routine(dinfo, G0_total, Tk, ne, nElement, xelem_ions, Np, co)
            import :: dp, DustChemistryInfo, ChargeOut
            class(DustChemistryInfo), intent(inout) :: dinfo
            real(dp), intent(in) :: G0_total, Tk, ne
            real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
            real(dp), intent(in) :: Np(:)
            type(ChargeOut), intent(inout) :: co
        end subroutine charge_routine
    end interface
    ! with and without the cell's RT photon densities (Np present or not in compute_dust_precool)
    procedure(charge_routine), pointer, save :: charge_with_Np => null(), charge_without_Np => null()
    real(dp), dimension(0), save :: no_Np
    ! dust_debug: the charging paths taken (0 full solve, 1 uniform table, 2 served, 4 predicted) and, for
    ! every RTG_DBG_EVERY-th served or predicted bin, a full solve of the same inputs to compare with:
    ! the largest |d<Z>|, |dGamma|/Gamma and |dLambda|/Lambda seen, per path (rtg_debug_check)
    integer, parameter :: RTG_DBG_EVERY = 1000, RTG_DBG_REPORT = 100000
    integer(kind=8), save :: rtg_dbg_n(0:4) = 0
    real(dp), save :: rtg_dbg_max(3, 2) = 0d0      ! (dZ, dGamma, dLambda; served, predicted)
    ! .true. during the cooling solver's T(1 + 1e-5) call of an RTZ sub-step, and its T(1 - 1e-5) call
    ! with rtz_dCdT_central (rtgroups_derivative_calls)
    logical, save :: rtg_dT_calls = .false.

contains
    subroutine store_pah_ion_fraction(dinfo, icell)
        ! Keep the ion fraction of each PAH bin of cell icell of the cooling vector (pah_fion_vec), from the
        ! charge distribution compute_dust_rad_rates has just set in dinfo
        implicit none
        class(DustChemistryInfo), intent(in) :: dinfo
        integer, intent(in) :: icell
        integer :: ii
        if (icell < 1 .or. icell > nvector) return
        do ii = 1, dinfo%npah
            pah_fion_vec(icell, ii) = 0d0
            if (pahbins_props(ii)%cation_start_idx <= pahbins_props(ii)%ncharge_states) &
                pah_fion_vec(icell, ii) = sum(dinfo%fcharge_pah(pahbins_props(ii)%cation_start_idx: &
                                                                pahbins_props(ii)%ncharge_states, ii))
        end do
    end subroutine store_pah_ion_fraction

    subroutine reset_cell_warm_starts(zguess)
        ! The guesses carried from one cell to the next (the ODE step of compute_dust_update, the grain-charge
        ! root of WDB06rt) are set for one cell: with the cooling queue (cool_queue >= 1) a cell's result then
        ! does not depend on which cell was solved before it, i.e. on the rank that solves it or the order.
        ! zguess: this cell's own last charge roots (cool_queue_zguess_ivar); without it the root search
        ! starts cold, which costs ~20% of the cooling in the WDB06rt H II box (a full bracketing per bin)
        implicit none
        real(dp), dimension(:), intent(in), optional :: zguess
        last_h_ode = 0.0_dp
        if (present(zguess)) then
            if (.not. allocated(rtg_Zguess)) allocate(rtg_Zguess(size(zguess)))
            rtg_Zguess = zguess
        else if (allocated(rtg_Zguess)) then
            rtg_Zguess = -huge(1d0)
        end if
    end subroutine reset_cell_warm_starts

    subroutine get_cell_zguess(zguess)
        ! The charge roots the last solve left (rtg_Zguess), to keep with the cell; unchanged if none
        implicit none
        real(dp), dimension(:), intent(inout) :: zguess
        if (allocated(rtg_Zguess)) zguess = rtg_Zguess
    end subroutine get_cell_zguess

    subroutine rtgroups_reset_warm(ncell)
        ! A new vector of ncell cells (rtz_solve_cooling): their warm states start cold. Only
        ! those cells' states and the no-cell slot (nvector + 1) are touched; a state is 84 + 2 NGROUPS
        ! reals (0.8 kB with 8 groups, most of them the prediction anchor), (nvector + 1) x ndust of them per rank.
        implicit none
        integer, intent(in) :: ncell
        if (allocated(rtg_warm)) then
            rtg_warm(1:min(ncell, nvector), :)%valid = .false.
            rtg_warm(1:min(ncell, nvector), :)%anchored = .false.
            rtg_warm(nvector + 1, :)%valid = .false.
            rtg_warm(nvector + 1, :)%anchored = .false.
        end if
        if (allocated(rtg_last)) rtg_last(:)%valid = .false.
        td_valid = .false.
    end subroutine rtgroups_reset_warm

    subroutine rtgroups_derivative_calls(on)
        ! WDB06rt: rtz_cool_step brackets the cooling solver's T(1 +- 1e-5) calls with on = .true.
        ! and .false. They need only dC/dT and have the inputs of the sub-step's precool solve but T
        ! (and, with rt_isIR, the dust IR emission it added to the IR group in between), so its bins
        ! are served from that solve's d/d ln T (rtg_served), instead of two more full solves whose
        ! warm-started roots need not agree.
        implicit none
        logical, intent(in) :: on
        rtg_dT_calls = on
    end subroutine rtgroups_derivative_calls

    subroutine compute_local_anisotropy_factor(dinfo,Fp,Np)
        ! Computes the local radiation anisotropy factor and solid angle 
        ! subtended by radiation sources based on the local radiation field
        ! dinfo --> the DustChemistryInfo instance to update with the computed values
        ! Fp --> the radiation flux vector (#/s/cm^2) for each radiation group
        ! Np --> the radiation energy density (#/cm^3) for each radiation group
        implicit none

        ! ---- Inputs -----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), dimension(:,:), intent(in) :: Fp
        real(dp), dimension(:), intent(in) :: Np

        ! ---- Local variables ----
        integer :: i
        real(dp) :: rad_ani

        ! the groups at the floor (skipped below) are isotropic and get no solid angle
        dinfo%local_rad_ani = 0d0
        dinfo%local_solid_angle = 0d0
        if (all(Np.le.dinfo%smallNp)) return

        if (fixed_rad_ani .eq. -1d0) then
            do i = 1, size(Np)
                if (Np(i) .le. dinfo%smallNp) cycle
                rad_ani = sqrt(sum(Fp(:,i)**2d0)) / (Np(i) * dinfo%local_c)
                dinfo%local_rad_ani(i) = rad_ani
                dinfo%local_solid_angle(i) = 2d0 * pi * (1d0 + (1d0 - rad_ani)**2d0) 
            end do
        else
            dinfo%local_rad_ani = fixed_rad_ani
            dinfo%local_solid_angle = 2d0 * pi * (1d0 + (1d0 - fixed_rad_ani)**2d0)
        end if
    end subroutine compute_local_anisotropy_factor

    subroutine compute_dust_rad_rates(dinfo, G0, Tk, ne,&
                                        &dustAbs, dustSc, dustRp,&
                                        &pahAbs, pahSc, pahRp)
        ! Computes the dust and PAH radiative rates (absorption, 
        ! scattering, and radiation pressure) for each radiation group, 
        ! given the local dust properties and radiation field.
        ! dinfo --> the DustChemistryInfo instance containing the dust properties
        ! G0 --> local radiation field strength in units of the Habing field
        ! Tk --> local gas temperature (K)
        ! ne --> local electron density (cm^-3)
        ! dustAbs --> output array for dust absorption rates [1/s]
        ! dustSc --> output array for dust scattering rates [1/s]
        ! dustRp --> output array for dust radiation pressure rates [1/s]
        ! pahAbs --> output array for PAH absorption rates [1/s]
        ! pahSc --> output array for PAH scattering rates [1/s]
        ! pahRp --> output array for PAH radiation pressure rates [1/s]
        use pah_photoelectric_heating, only: interpolate_pah_charge_equilibrium
        use dust_radiation, only: rad_dust_rate, rad_pah_rate

        implicit none

        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0, Tk, ne
        ! ---- Outputs ----
        real(dp), dimension(:), intent(inout) :: dustAbs, dustSc, dustRp
        real(dp), dimension(:), intent(inout) :: pahAbs, pahSc, pahRp
        ! ---- Local variables ----
        integer :: ii

        ! 1. Compute and add the different dust radiative rates
        dustAbs = 0d0; dustSc = 0d0; dustRp = 0d0
        if (dinfo%ndust > 0) then
            do ii = 1, dinfo%nGroups
                dustAbs(ii) = rad_dust_rate(dinfo%csa_dust(:,ii),dinfo%rho_dust(:)) ! [1/s]
                dustSc (ii) = rad_dust_rate(dinfo%css_dust(:,ii),dinfo%rho_dust(:)) ! [1/s]
                dustRp (ii) = rad_dust_rate(dinfo%csr_dust(:,ii),dinfo%rho_dust(:)) ! [1/s]
            end do
            dinfo%dustAbs = dustAbs
        end if

        ! 2. Compute the different PAH radiative rates
        pahAbs = 0d0; pahSc = 0d0; pahRp = 0d0
        if (dinfo%npah > 0) then
            ! Before we move onto computing the PAH contribution
            ! to absorption and scattering, we need to compute the
            ! PAH charge distribution, since the PAH cross-sections
            ! depend on the charge state
            ! TODO: currently we use the interpolation tables since
            ! it's sufficiently accurate and much faster, but we 
            ! could change to use the full calculation if we wanted
            do ii = 1, dinfo%npah
                call interpolate_pah_charge_equilibrium(ii,G0,ne,Tk,&
                   &dinfo%fcharge_pah(:,ii))
            end do
            do ii = 1, dinfo%nGroups
                pahAbs(ii) = pahAbs(ii)+ rad_pah_rate(dinfo%csa_pah(:,ii),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
                pahSc (ii) = pahSc(ii) + rad_pah_rate(dinfo%css_pah(:,ii),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
                pahRp (ii) = pahRp(ii) + rad_pah_rate(dinfo%csr_pah(:,ii),dinfo%rho_pah(:),dinfo%fcharge_pah(:,:)) ! [1/s]
            end do
            dinfo%pahAbs = pahAbs
        end if

        ! 3. Add the contribution from dust and PAHs for the 
        !    returned absorption, scattering, and radiation pressure rates
        dustAbs = dustAbs + pahAbs
        dustSc = dustSc + pahSc
        dustRp = dustRp + pahRp
    end subroutine compute_dust_rad_rates

    subroutine compute_dust_precool(dinfo, G0_total, Tk, ne,&
                                    &nElement, xelem_ions, nH2, nCO, &
                                    & Np, rates_only)
        ! Pre-computes the dust and PAH cooling and heating rates for the given local conditions, 
        ! which can then be used in subsequent calls to compute_dust_coolrates to save computational time.
        ! dinfo --> the DustChemistryInfo instance to update with the computed rates
        ! G0_total --> local radiation field strength in units of the Habing field
        ! Tk --> local gas temperature (K)
        ! ne --> local electron density (cm^-3)
        ! nElement --> array of number densities for each element (cm^-3)
        ! xelem_ions --> array of ionization fractions for each element and ionization state
        ! nH2 --> local molecular hydrogen density (cm^-3)
        ! nCO --> local carbon monoxide density (cm^-3)
        ! Np --> the radiation energy density (#/cm^3) for each radiation group
        use dust_charging, only: compute_Coulomb_focusing_ions, isrf_lookup, ISRF_ND, ISRF_DCHARGE
        use dust_photoelectric_heating, only: interpolate_dust_peh_rate
        use pah_photoelectric_heating, only: interpolate_pah_charge_equilibrium,&
                                            compute_pah_peh_equilibrium,&
                                            interpolate_pah_peh_equilibrium,&
                                            compute_pah_charge_equilibrium
        use dust_radiation, only: update_T_dust, serve_T_dust
        use dust_surface_chemistry, only: grain_h2_formation_rate
#ifdef RT
        use rt_parameters, only: rt_isIR
#ifdef RTZ
        use rt_parameters, only: rtz_equilibrium_test, rtz_single_cell_test
#endif
#endif

        implicit none

        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: nH2, nCO
        real(dp), dimension(1:dinfo%nGroups), intent(in), optional :: Np
        logical, intent(in), optional :: rates_only   ! only the rates compute_dust_coolrates returns (not the Coulomb factors)

        ! ---- Local variables ----
        integer :: ii,j
        integer :: i_neutral, i_charged   ! csa_pah row indices: neutral=2*ii-1, charged=2*ii
        real(dp) :: nHI
        logical :: skip_coulomb
        logical :: need_T, need_Z   ! the grain temperature, the grain charge, are used
        type(ChargeOut) :: co
        logical :: rec

        ! the grain-assisted ion recombination rates in the RTZ step's call at T only, not in the
        ! cooling solver's T(1 + 1e-5) call for the rates (rates_only)
        rec = dust_ion_recombination
        if (present(rates_only)) rec = rec .and. .not. rates_only

        if (dinfo%ndust > 0) then
            do ii = 1, dinfo%ndust
                dinfo%n_dust(ii) = dinfo%rho_dust(ii) / dustbins_props(ii)%mgrain   ! grains per cm^3
            end do
            ! The grain temperature and charge only where something uses them: T_dust in the gas-grain
            ! collisional heating, H2 formation, the dust and PAH processes, the dust IR emission and
            ! the equilibrium-test output; the charge in those (through the T_dust balance), the PE
            ! heating and recombination cooling, the Coulomb factors, the grain-assisted recombination
            ! and charged sputtering. Otherwise both were solved in every RTZ substep and discarded
            ! (the TVA Draine drag solves its own charge, dust_radpressure).
            need_T = dust_coll_cooling .or. H2ondust .or. ndust_processes > 0 .or. npah_processes > 0
#ifdef RT
            need_T = need_T .or. rt_isIR
#ifdef RTZ
            need_T = need_T .or. rtz_equilibrium_test > 0 .or. rtz_single_cell_test
#endif
#endif
            need_Z = need_T .or. dust_pe_heating .or. Coulomb_precompute .or. dust_ion_recombination &
                     .or. dust_sputtering_charge
            ! 1. Compute the equilibrium dust charge
            co%has_pe = .false.
            co%has_alpha = .false.
            co%want_alpha = rec
            co%Dmode = 0
            co%alpha = 0d0
            if (need_Z) then
                ! the charging routines of the run's charging_model, chosen once (charge_select)
                if (.not. associated(charge_with_Np)) call charge_select()
                if (present(Np)) then
                    call charge_with_Np(dinfo, G0_total, Tk, ne, nElement, xelem_ions, Np, co)
                else
                    call charge_without_Np(dinfo, G0_total, Tk, ne, nElement, xelem_ions, no_Np, co)
                end if
            end if
            if (rec .and. co%has_alpha) then
                ! grain-assisted X+ -> X per X+ ion [s^-1]: sum over the bins of n_grain alpha (dust bins only)
                dinfo%rec_ion_rate = 0d0
                do ii = 1, dinfo%ndust
                    if (dinfo%rho_dust(ii) <= 0d0) cycle
                    do j = 1, RTG_NRI
                        dinfo%rec_ion_rate(RI_ATOMIC(j)) = dinfo%rec_ion_rate(RI_ATOMIC(j)) + dinfo%n_dust(ii)*co%alpha(j, ii)
                    end do
                end do
            end if

            ! 2. If needed, precompute the Coulomb factors, for compute_dust_update: not in a call for
            ! the rates only (the cooling solver's T(1 + 1e-5) call), whose factors differ from those
            ! at T only to first order
            skip_coulomb = .false.
            if (present(rates_only)) skip_coulomb = rates_only
            if (.not. skip_coulomb) then
                dinfo%Coulomb_factor = 1d0
                if (Coulomb_precompute) then
                    do ii = 1, dinfo%ndust
                        call compute_Coulomb_focusing_ions(ii,Tk,dinfo%Z_dust(ii),dinfo%Z_sigma(ii),dinfo%nion_charges,&
                                                           dinfo%Coulomb_factor(ii,-1:dinfo%nion_charges))
                        ! the low impactor charges over the exact P(Z), where the Gaussian of <Z> and sigma_Z
                        ! misses the tail at the opposite sign (small grains in cold gas): the charging
                        ! routine's correction (Dmode 1) or factors (Dmode 2)
                        do j = 1, ISRF_ND
                            if (ISRF_DCHARGE(j) > dinfo%nion_charges) cycle
                            if (co%Dmode(ii) == 1) then
                                dinfo%Coulomb_factor(ii, ISRF_DCHARGE(j)) = &
                                    max(dinfo%Coulomb_factor(ii, ISRF_DCHARGE(j)) * co%D(j, ii), 1d-10)
                            else if (co%Dmode(ii) == 2) then
                                dinfo%Coulomb_factor(ii, ISRF_DCHARGE(j)) = max(co%D(j, ii), 1d-10)
                            end if
                        end do
                    end do
                end if
            end if

            ! 3. Compute the equilibrium dust photoelectric heating and recombination cooling rates
            if (dust_pe_heating .and. co%has_pe) then
                dinfo%Pinj_dust(1:dinfo%ndust) = co%Pinj(1:dinfo%ndust)
                dinfo%Prec_dust(1:dinfo%ndust) = co%Prec(1:dinfo%ndust)
            elseif (dust_pe_heating) then
                do ii = 1, dinfo%ndust
                    call interpolate_dust_peh_rate(ii,G0_total,ne,Tk,&
                                                    &dinfo%Pinj_dust(ii),&
                                                    &dinfo%Prec_dust(ii))
                end do
            end if

            ! 4. Compute the internal energy of the dust grain considering all heating and cooling processes
            if (.not. need_T) then
                ! nothing uses T_dust, P_rad or the collisional heating (0 without dust_coll_cooling)
            else if (td_served()) then
                ! the cooling solver's T(1 +- 1e-5) call: first order from the last update_T_dust
                call serve_T_dust(td_lin,Tk,dinfo%Prec_dust(:),dinfo%Pinj_dust(:),dinfo%Z_dust(:),dinfo%Pcoll_dust(:),&
                                  dinfo%Prad_dust(:),dinfo%T_dust(:))
            else
                if (.not. allocated(td_lin)) allocate(td_lin(1:dinfo%ndust), td_nel(size(nElement)), &
                                                      td_xion(size(xelem_ions,1), size(xelem_ions,2)))
                if (present(Np)) then
                    call update_T_dust(dinfo%G0_background,dinfo%Pcoll_dust(:),dinfo%Prec_dust(:),dinfo%Pinj_dust(:),dinfo%Prad_dust(:),&
                                        dinfo%T_dust(:),ne,nElement(:),xelem_ions(:,:),dinfo%Coulomb_factor(:,:),nH2,nCO,Tk,dinfo%Z_dust(:),&
                                        Np(:)*dinfo%group_eV(:),dinfo%csa_dust(:,:),lin=td_lin)
                    if (.not. allocated(td_Np)) allocate(td_Np(dinfo%nGroups), td_egy(dinfo%nGroups))
                    td_Np = Np
                    td_egy = dinfo%group_eV
                    td_G0 = dinfo%G0_background
                else
                    call update_T_dust(G0_total,dinfo%Pcoll_dust(:),dinfo%Prec_dust(:),dinfo%Pinj_dust(:),dinfo%Prad_dust(:),&
                                        dinfo%T_dust(:),ne,nElement(:),xelem_ions(:,:),dinfo%Coulomb_factor(:,:),nH2,nCO,Tk,dinfo%Z_dust(:),&
                                        lin=td_lin)
                    td_G0 = G0_total
                end if
                td_valid = .true.
                td_hasNp = present(Np)
                td_T = Tk
                td_ne = ne
                td_nH2 = nH2
                td_nCO = nCO
                td_nel = nElement
                td_xion = xelem_ions
            end if

            ! 5. Now convert all the rates from erg/s per grain to erg/s/cm3
            do ii = 1, dinfo%ndust
                dinfo%Pcoll_dust(ii) = dinfo%Pcoll_dust(ii) * dinfo%n_dust(ii)
                ! print*,'ii,dinfo%Pcoll_dust(ii),dinfo%n_dust(ii): ',ii,dinfo%Pcoll_dust(ii),dinfo%n_dust(ii)
                dinfo%Prec_dust(ii) = dinfo%Prec_dust(ii) * dinfo%n_dust(ii)
                dinfo%Pinj_dust(ii) = dinfo%Pinj_dust(ii) * dinfo%n_dust(ii)
                dinfo%Prad_dust(ii) = dinfo%Prad_dust(ii) * dinfo%n_dust(ii)
            end do

            ! 6. Compute the H2 formation rate on dust grains
            if (H2ondust) then
                nHI = nElement(1) * xelem_ions(1,1)
                dinfo%H2_formation_rate = grain_h2_formation_rate(nHI,nElement(1),Tk,dinfo%rho_dust(:),dinfo%T_dust(:))
            end if
        end if

        if (dinfo%npah > 0) then
            ! 7. Compute the PAH PEH model
            if (pah_pe_heating .and. present(Np)) then
                if (pah_pe_heating_isrf) then
                    ! We have rt, but we want to just use the ISRF-averaged model
                        do ii = 1, dinfo%npah
                            call interpolate_pah_peh_equilibrium(ii,G0_total,&
                                                                ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                                dinfo%Pabs_pah(ii,1),dinfo%Pinj_pah(ii),&
                                                                dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                        end do
                else
                    ! We want the PAH PEH also have rt, so we use the full model
                    do ii = 1, dinfo%npah
                        i_neutral = 2*ii - 1  ! row for neutral PAH cross-sections in csa_pah
                        i_charged = 2*ii      ! row for charged PAH cross-sections in csa_pah
                        ! Consider the contribution from the UV background first
                        call interpolate_pah_peh_equilibrium(ii,dinfo%G0_background,&
                                                            ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                            dinfo%Pabs_pah(ii,1),dinfo%Pinj_pah(ii),&
                                                            dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                        ! As for the grains: the recombination cooling is added once, by the full
                        ! model below, at the charge state of the total field
                        dinfo%Prec_pah(ii) = 0d0
                        ! And now the full model for the local radiation field
                        call compute_pah_peh_equilibrium(ii,dinfo%csa_pah(i_neutral,:),&
                                                            dinfo%csa_pah(i_neutral,:),&
                                                            dinfo%csa_pah(i_charged,:),&
                                                            dinfo%csa_pah(i_charged,:),&
                                                            dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                            Np(:),dinfo%group_eV(:),&
                                                            dinfo%local_c,Tk,ne,dinfo%fcharge_pah(:,ii),&
                                                            dinfo%Pabs_pah(ii,:),dinfo%Pinj_pah(ii),&
                                                            dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                    end do
                end if
            else if (pah_pe_heating) then
                ! We want PAH PEH but don't have rt
                do ii = 1, dinfo%npah
                    call interpolate_pah_peh_equilibrium(ii,dinfo%G0_background,&
                                                        ne,Tk,dinfo%fcharge_pah(:,ii),&
                                                        dinfo%Pabs_pah(ii,1),dinfo%Pinj_pah(ii),&
                                                        dinfo%Prad_pah(ii),dinfo%Prec_pah(ii))
                end do
            else if (present(Np)) then
                ! We don't want PAH PEH but have rt, so we still want to compute the PAH charge distribution
                do ii = 1, dinfo%npah
                    i_neutral = 2*ii - 1
                    i_charged = 2*ii
                    call compute_pah_charge_equilibrium(ii,dinfo%csa_pah(i_neutral,:),&
                                                        dinfo%csa_pah(i_neutral,:),&
                                                        dinfo%csa_pah(i_charged,:),&
                                                        dinfo%csa_pah(i_charged,:),&
                                                        dinfo%nGroups,dinfo%local_solid_angle(:),&
                                                        Np(:),dinfo%group_eV(:),dinfo%local_c,&
                                                        Tk,ne,dinfo%fcharge_pah(:,ii))
                end do
            else
                ! We don't want PAH PEH but don't have rt, so we just interpolate
                ! the ISRF-averaged PAH charge tables
                do ii = 1, dinfo%npah
                    call interpolate_pah_charge_equilibrium(ii,G0_total,ne,Tk,&
                       &dinfo%fcharge_pah(:,ii))
                end do
            end if

            ! 8. Now convert from erg/s to erg/s/cm3
            do ii = 1, dinfo%npah
                dinfo%n_pah(ii) = dinfo%rho_pah(ii) / pahbins_props(ii)%mpah
                dinfo%Pabs_pah(ii,:) = dinfo%Pabs_pah(ii,:) * dinfo%n_pah(ii)
                dinfo%Pinj_pah(ii) = dinfo%Pinj_pah(ii) * dinfo%n_pah(ii)
                dinfo%Prad_pah(ii) = dinfo%Prad_pah(ii) * dinfo%n_pah(ii)
                dinfo%Prec_pah(ii) = dinfo%Prec_pah(ii) * dinfo%n_pah(ii)
            end do
        end if
    contains


        logical function td_served()
            ! same inputs as the last update_T_dust but T, within RTG_SERVE in ln T (and the
            ! recombination and photoelectric terms and the grain charge, which serve_T_dust takes to
            ! first order)
            real(dp) :: x
            td_served = .false.
            if (.not. td_valid .or. (td_hasNp .neqv. present(Np))) return
            x = Tk/td_T
            if (x > RTG_SERVE_HI .or. x < RTG_SERVE_LO) return
            if (ne /= td_ne .or. nH2 /= td_nH2 .or. nCO /= td_nCO) return
            if (present(Np)) then
                if (dinfo%G0_background /= td_G0) return
                if (any(Np /= td_Np) .or. any(dinfo%group_eV /= td_egy)) return
            else
                if (G0_total /= td_G0) return
            end if
            if (any(nElement /= td_nel) .or. any(xelem_ions /= td_xion)) return
            td_served = .true.
        end function td_served


    end subroutine compute_dust_precool

    subroutine charge_select()
        ! The charging routines of the run's charging_model (checked in dust_init): WDB06rt on the RT
        ! photon groups, and its (gamma, T) charge fits for a call without them; WDB06isrf in both cases
        implicit none
        if (trim(charging_model) == 'WDB06rt') then
            charge_with_Np => charge_rtgroups
            charge_without_Np => charge_fits
        else
            charge_with_Np => charge_isrf
            charge_without_Np => charge_isrf
        end if
    end subroutine charge_select

    subroutine charge_rtgroups(dinfo, G0_total, Tk, ne, nElement, xelem_ions, Np, co)
        ! charging_model = 'WDB06rt' with the cell's RT photons
        use dust_charging, only: rtgroups_solve_bin, rtgroups_refine_alpha, rtgroups_uniform_lookup, &
                                 rtgroups_predict, rtgroups_coulomb_ratio, isrf_lookup
        implicit none
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: Np(:)
        type(ChargeOut), intent(inout) :: co
        real(dp) :: zz, ss, gg, ll
        real(dp) :: n_Hp, n_Hep, n_Hepp, h
        type(RTGState), pointer :: st
        type(RTGResult) :: rr
        integer :: path, ic
        integer, dimension(1:ndust) :: paths
        real(dp) :: rate(RTG_NRI, 1:ndust), tot(RTG_NRI)
        integer :: ii
        logical :: no_local, predict, ok

        co%has_pe = .true.
        co%has_alpha = .true.
        ! Local-field charge balance on the RT photon groups (WDB06 yields),
        ! which also gives the PE heating, recombination cooling and the
        ! grain-assisted ion recombination (rr%alpha; dust_ion_recombination: dinfo%rec_ion_rate, for the RTZ chemistry)
        n_Hp = nElement(1)*xelem_ions(1,2)
        n_Hep = nElement(2)*xelem_ions(2,2)
        n_Hepp = nElement(2)*xelem_ions(2,3)
        no_local = all(Np .le. dinfo%smallNp)
        if (.not. allocated(rtg_warm)) call rtg_allocate()
        ic = dinfo%icell
        if (ic < 1 .or. ic > nvector) ic = nvector + 1
        predict = rtg_predict_tol > 0d0 .and. ic <= nvector
        do ii = 1, dinfo%ndust
            paths(ii) = -1
            co%Dmode(ii) = 0
            if (dinfo%rho_dust(ii) <= 0d0) then
                ! empty bin: no charge or rates needed
                dinfo%Z_dust(ii) = 0d0
                dinfo%Z_sigma(ii) = 0d0
                co%Pinj(ii) = 0d0
                co%Prec(ii) = 0d0
                cycle
            end if
            st => rtg_warm(ic, ii)
            if (ic == nvector + 1) st%valid = .false.
            rr = RTGResult()
            ok = .false.
            if (no_local .and. .not. dust_ion_recombination) then
                ! only the uniform background: start-up table in (T, ne), which has no alpha
                path = 1
                call rtgroups_uniform_lookup(ii, dinfo%G0_background, Tk, ne, rr%Zmean, rr%Zsigma, rr%Gamma, rr%Lambda)
                if (Coulomb_precompute) then
                    call isrf_lookup(ii, dinfo%G0_background, Tk, ne, zz, ss, gg, ll, co%D(:, ii))
                    co%Dmode(ii) = 2
                end if
            else if (rtg_served(ii)) then
                ! the cooling solver's T(1 +- 1e-5) call: first order in ln T from the last solve
                path = 2
                h = log(Tk/rtg_last(ii)%T)
                rr%Zmean = rtg_last(ii)%r%Zmean + h*rtg_last(ii)%r%dZmean
                rr%Zsigma = rtg_last(ii)%r%Zsigma + h*rtg_last(ii)%r%dZsigma
                rr%Gamma = rtg_last(ii)%r%Gamma + h*rtg_last(ii)%r%dGamma
                rr%Lambda = rtg_last(ii)%r%Lambda + rtg_last(ii)%r%Lambda_auto + h*rtg_last(ii)%r%dLambda
                rr%Zstar = h
            else
                ! first order from this cell's last full solve where the inputs changed little
                if (predict) call rtgroups_predict(ii, st, Np, dinfo%local_c, Tk, ne, n_Hp, n_Hep, n_Hepp, &
                                                   dinfo%G0_background, rtg_predict_tol, rtg_predict_tol_shape, rr, ok)
                if (ok) then
                    path = 4
                    if (Coulomb_precompute) then
                        co%D(:, ii) = st%Dratio
                        co%Dmode(ii) = 1
                    end if
                else
                    ! full solve; the recombination of a wide P(Z) is the Gaussian estimate here
                    path = 0
                    call rtgroups_solve_bin(ii, Np, dinfo%group_eV, dinfo%local_c, Tk, ne, &
                                            n_Hp, n_Hep, n_Hepp, dinfo%G0_background, rtg_Zguess(ii), &
                                            st, RTG_DLNT, merge(RTG_ALPHA_GAUSS, RTG_ALPHA_NONE, dust_ion_recombination), &
                                            rr, sens=predict)
                    if (Coulomb_precompute) then
                        ! narrow P(Z): its exact Coulomb factors over the Gaussian ones (kept for the
                        ! predictions from this solve); a wide one is close to its Gaussian
                        st%Dratio = 1d0
                        if (rr%discrete) call rtgroups_coulomb_ratio(ii, Np, dinfo%local_c, Tk, ne, n_Hp, n_Hep, &
                            n_Hepp, dinfo%G0_background, rr%Zstar, rr%Zsigma, rr%Zmean, st%Dratio)
                        co%D(:, ii) = st%Dratio
                        co%Dmode(ii) = 1
                    end if
                end if
                rtg_last(ii)%valid = .true.
                rtg_last(ii)%T = Tk
                rtg_last(ii)%ne = ne
                rtg_last(ii)%n_Hp = n_Hp
                rtg_last(ii)%n_Hep = n_Hep
                rtg_last(ii)%n_Hepp = n_Hepp
                rtg_last(ii)%G0 = dinfo%G0_background
                rtg_last(ii)%c_red = dinfo%local_c
                rtg_last(ii)%Np = Np
                rtg_last(ii)%egy = dinfo%group_eV
                rtg_last(ii)%r = rr
                rr%Lambda = rr%Lambda + rr%Lambda_auto
                rtg_Zguess(ii) = rr%Zstar
            end if
            paths(ii) = path
            if (dust_debug) call rtg_debug_check(ii, path, rr)
            dinfo%Z_dust(ii) = rr%Zmean
            dinfo%Z_sigma(ii) = rr%Zsigma
            co%Pinj(ii) = rr%Gamma
            co%Prec(ii) = rr%Lambda
        end do
        ! grain-assisted recombination (pyCALIMA solve_cell): the fit only for the wide bins
        ! solved here that carry more than RTG_RECOMB_IMPORTANT of the cell's rate of some ion
        ! (the predicted bins count in that rate)
        if (dust_ion_recombination .and. any(paths == 0)) then
            rate = 0d0
            do ii = 1, dinfo%ndust
                if (paths(ii) == 0 .or. paths(ii) == 4) &
                    rate(:, ii) = (dinfo%rho_dust(ii)/dustbins_props(ii)%mgrain)*rtg_last(ii)%r%alpha
            end do
            tot = sum(rate, dim=2)
            do ii = 1, dinfo%ndust
                if (paths(ii) /= 0) cycle
                if (rtg_last(ii)%r%discrete) cycle
                if (maxval(rate(:, ii)/merge(tot, 1d0, tot > 0d0), mask=tot > 0d0) <= RTG_RECOMB_IMPORTANT) cycle
                st => rtg_warm(ic, ii)
                call rtgroups_refine_alpha(ii, Np, dinfo%group_eV, dinfo%local_c, Tk, ne, n_Hp, n_Hep, n_Hepp, &
                                           dinfo%G0_background, rtg_last(ii)%r)
                if (st%anchored) st%ar%alpha = rtg_last(ii)%r%alpha      ! predictions scale the fitted alpha
            end do
        end if
        if (co%want_alpha) then
            ! alpha of the last full solve of each bin: this cell's (with the flag on there is no uniform path)
            do ii = 1, dinfo%ndust
                if (paths(ii) >= 0 .and. rtg_last(ii)%valid) co%alpha(:, ii) = rtg_last(ii)%r%alpha
            end do
        end if
    contains
        subroutine rtg_debug_check(ib, path, rr)
            ! dust_debug: count the path, and check a served or predicted result against a full solve of
            ! the same inputs, from a cold state (st and rtg_last are left as they are)
            use amr_commons, only: myid
            integer, intent(in) :: ib, path
            type(RTGResult), intent(in) :: rr
            type(RTGState) :: st0
            type(RTGResult) :: r0
            integer :: k
            rtg_dbg_n(path) = rtg_dbg_n(path) + 1
            if ((path == 2 .or. path == 4) .and. mod(rtg_dbg_n(path), int(RTG_DBG_EVERY, 8)) == 0) then
                k = merge(1, 2, path == 2)
                call rtgroups_solve_bin(ib, Np, dinfo%group_eV, dinfo%local_c, Tk, ne, n_Hp, n_Hep, n_Hepp, &
                                        dinfo%G0_background, -huge(1d0), st0, 0d0, RTG_ALPHA_NONE, r0)
                rtg_dbg_max(1, k) = max(rtg_dbg_max(1, k), abs(rr%Zmean - r0%Zmean))
                if (r0%Gamma /= 0d0) rtg_dbg_max(2, k) = max(rtg_dbg_max(2, k), abs(rr%Gamma/r0%Gamma - 1d0))
                if (r0%Lambda + r0%Lambda_auto /= 0d0) rtg_dbg_max(3, k) = max(rtg_dbg_max(3, k), &
                    abs(rr%Lambda/(r0%Lambda + r0%Lambda_auto) - 1d0))
            end if
            if (myid == 1 .and. mod(sum(rtg_dbg_n), int(RTG_DBG_REPORT, 8)) == 0) &
                write(*,'(A,4I12,A,3ES9.2,A,3ES9.2)') ' WDB06rt paths (full, table, served, predicted):', &
                    rtg_dbg_n(0), rtg_dbg_n(1), rtg_dbg_n(2), rtg_dbg_n(4), '; vs full solve, served dZ dG dL', &
                    rtg_dbg_max(:, 1), ', predicted', rtg_dbg_max(:, 2)
        end subroutine rtg_debug_check

        logical function rtg_served(ib)
            ! same inputs as the last full solve of bin ib but T, within RTG_SERVE in ln T: in the
            ! sub-step's T(1 +- 1e-5) calls by construction, else checked
            integer, intent(in) :: ib
            real(dp) :: x
            rtg_served = .false.
            if (.not. rtg_last(ib)%valid) return
            x = Tk/rtg_last(ib)%T
            if (x > RTG_SERVE_HI .or. x < RTG_SERVE_LO) return
            if (rtg_dT_calls) then
                rtg_served = .true.
                return
            end if
            if (ne /= rtg_last(ib)%ne .or. n_Hp /= rtg_last(ib)%n_Hp .or. n_Hep /= rtg_last(ib)%n_Hep .or. &
                n_Hepp /= rtg_last(ib)%n_Hepp .or. dinfo%G0_background /= rtg_last(ib)%G0 .or. &
                dinfo%local_c /= rtg_last(ib)%c_red) return
            if (any(Np /= rtg_last(ib)%Np) .or. any(dinfo%group_eV /= rtg_last(ib)%egy)) return
            rtg_served = .true.
        end function rtg_served
        subroutine rtg_allocate()
            integer :: ib
            allocate(rtg_warm(1:nvector + 1, 1:dinfo%ndust), rtg_last(1:dinfo%ndust))
            if (.not. allocated(rtg_Zguess)) then      ! reset_cell_warm_starts may have set it already
                allocate(rtg_Zguess(1:dinfo%ndust))
                rtg_Zguess = -huge(1d0)
            end if
            do ib = 1, dinfo%ndust
                allocate(rtg_last(ib)%Np(dinfo%nGroups), rtg_last(ib)%egy(dinfo%nGroups))
            end do
        end subroutine rtg_allocate
    end subroutine charge_rtgroups

    subroutine charge_isrf(dinfo, G0_total, Tk, ne, nElement, xelem_ions, Np, co)
        ! charging_model = 'WDB06isrf': the uniform-ISRF tables at the local G0: the charge, the PE heating
        ! and recombination cooling, the Coulomb factors of the low impactor charges over the exact P(Z)
        ! (Dmode 2), and the grain-assisted recombination (dust_ion_recombination). Np is not used.
        use dust_charging, only: isrf_lookup
        implicit none
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: Np(:)
        type(ChargeOut), intent(inout) :: co
        integer :: ii
        co%has_pe = .true.
        co%has_alpha = .true.
        do ii = 1, dinfo%ndust
            if (co%want_alpha) then
                call isrf_lookup(ii, G0_total, Tk, ne, dinfo%Z_dust(ii), dinfo%Z_sigma(ii), co%Pinj(ii), co%Prec(ii), &
                                 co%D(:, ii), co%alpha(:, ii))
            else
                call isrf_lookup(ii, G0_total, Tk, ne, dinfo%Z_dust(ii), dinfo%Z_sigma(ii), co%Pinj(ii), co%Prec(ii), &
                                 co%D(:, ii))
            end if
            co%Dmode(ii) = 2
        end do
    end subroutine charge_isrf

    subroutine charge_fits(dinfo, G0_total, Tk, ne, nElement, xelem_ions, Np, co)
        ! WDB06rt in a call without the cell's RT photons: the mean charge and its width of the uniform-ISRF
        ! tables at G0_total, one lookup per bin; the PE heating then comes from interpolate_dust_peh_rate
        ! (has_pe stays false). Np is not used.
        use dust_charging, only: isrf_lookup
        implicit none
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: Np(:)
        type(ChargeOut), intent(inout) :: co
        integer :: ii
        real(dp) :: gam, lam, Dj(ISRF_ND)
        do ii = 1, dinfo%ndust
            call isrf_lookup(ii, G0_total, Tk, ne, dinfo%Z_dust(ii), dinfo%Z_sigma(ii), gam, lam, Dj)
        end do
    end subroutine charge_fits

    subroutine compute_dust_coolrates(dinfo, G0_total, Tk, ne,&
                                    &nElement, xelem_ions, nH2, nCO, &
                                    &total_rec_power, total_inj_power, total_col_power,&
                                    &H2_formation_rate,Np)
        ! Computes the dust and PAH cooling and heating rates for the given local conditions,
        ! using the pre-computed rates from compute_dust_precool if available to save computational time.
        ! dinfo --> the DustChemistryInfo instance to update with the computed rates
        ! G0_total --> local radiation field strength in units of the Habing field
        ! Tk --> local gas temperature (K)
        ! ne --> local electron density (cm^-3)
        ! nElement --> array of number densities for each element (cm^-3)
        ! xelem_ions --> array of ionization fractions for each element and ionization state
        ! nH2 --> local molecular hydrogen density (cm^-3)
        ! nCO --> local carbon monoxide density (cm^-3)
        ! total_rec_power --> output total recombination cooling power from dust and PAHs (erg/s/cm^3)
        ! total_inj_power --> output total photoelectric heating power from dust and PAHs (erg/s/cm^3)
        ! total_col_power --> output total collisional cooling power from dust (erg/s/cm^3)
        ! H2_formation_rate --> output H2 formation rate on dust grains (cm^3/s)
        ! Np --> the radiation energy density (#/cm^3) for each radiation group
        implicit none
        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: G0_total, Tk, ne
        real(dp), intent(in) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in) :: nH2, nCO
        real(dp), intent(out) :: total_rec_power, total_inj_power, total_col_power
        real(dp), intent(out) :: H2_formation_rate
        real(dp), dimension(1:dinfo%nGroups), intent(in), optional :: Np

        ! If precomputed rates are not already in dinfo (i.e. compute_dust_precool was not
        ! called before us this step), delegate to it now.  This eliminates ~148 lines
        ! of duplicated physics that previously lived in the non-precomp branch here.
        if (.not. dinfo%use_precomp) then
            call compute_dust_precool(dinfo, G0_total, Tk, ne, nElement, xelem_ions, nH2, nCO, Np=Np, rates_only=.true.)
        end if

        ! Read the precomputed rates from dinfo (populated either just above or by an
        ! external compute_dust_precool call earlier in the RTZ step).
        total_rec_power = 0d0
        total_inj_power = 0d0
        total_col_power = 0d0
        H2_formation_rate = -1d0
        if (dinfo%ndust > 0) then
            total_rec_power = sum(dinfo%Prec_dust)
            total_inj_power = sum(dinfo%Pinj_dust)
            total_col_power = sum(dinfo%Pcoll_dust)
            if (H2ondust) then
                H2_formation_rate = dinfo%H2_formation_rate
            end if
        end if
        if (dinfo%npah > 0) then
            total_rec_power = total_rec_power + sum(dinfo%Prec_pah)
            total_inj_power = total_inj_power + sum(dinfo%Pinj_pah)
        end if
        dinfo%use_precomp = .false.
    end subroutine compute_dust_coolrates

    subroutine compute_dust_update(dinfo,nElement,xelem_ions,dt,Np,step_ok, dUU)

        use ode_driver_mod, only: integrate_dust_ode
        use ode_interface_mod, only: dust_solver_step
        use dust_rhs_mod, only: dust_rhs
        use dust_rates, only: compute_rate_caches
        use dust_radiative_torques, only: total_radiative_torque,IR_damping_factor
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef RTZ
        use rtz_module, only:elements
#endif

        implicit none

        ! ---- Inputs ----
        class(DustChemistryInfo), intent(inout) :: dinfo
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: nElement(:), xelem_ions(:,:)
        real(dp), intent(in), optional :: Np(:)
        logical, intent(out) :: step_ok

        ! --- Local variables ----
        integer :: ii, jj, ndust_total
        real(dp) :: sum_check
        real(dp), pointer :: y_gas(:,:), y_gas_out(:,:)
        real(dp), pointer :: y_dust(:), y_dust_out(:)
        real(dp) :: mass_g(n_elements)   ! element atomic masses — populated once via #ifdef, used throughout
        real(dp) :: h_init_val           ! substep hint for integrate_dust_ode

        real(dp), intent(out) :: dUU

        dUU = 0.0_dp ! should I keep previous values?

        ! Populate mass_g once so the pack/unpack loops below are #ifdef-free.
#ifdef RTZ
        do ii = 1, n_elements
            mass_g(ii) = elements(ii)%atomic_mass_g
        end do
#else
        mass_g(:) = el_atomic_masses_g(1:n_elements)
#endif

        ! If no dust or PAH process is active, just return
        if (ndust_processes.eq.0 .and. npah_processes.eq.0) then
            step_ok = .true.
            return
        end if

        ! Precompute/cache rate factors for this cell-update step
        call compute_rate_caches(dinfo, nElement)

        ! 1. Compute the local RAT-D quantities if we run with dust_ratd
        if (dust_ratd) then
            if (present(Np)) then
                do ii = 1, dinfo%ndust
                    dinfo%rat_torque(ii) = total_radiative_torque(dinfo%local_rad_ani,Np,&
                                                            dinfo%group_eV(:)*eV2erg,dinfo%nGroups,&
                                                            dinfo%local_c,dinfo%csrat_dust(:,ii))
                end do
            end if
            dinfo%IR_damp_factor = IR_damping_factor(dinfo%local_G0*1.13d0,nElement(1),dinfo%local_Tk,dinfo%T_dust)
            do ii = 1, dinfo%ndust
                dinfo%rat_torque(ii) = dinfo%rat_torque(ii) + dustbins_props(ii)%RAT_torque_0 * dinfo%G0_background
            end do
        end if

        ! 2. Now setup the arrays for gas and dust quantities to send to the ODE solver
        ndust_total = 0
        if (dinfo%ndust > 0) ndust_total = ndust_total + dinfo%ndust
        if (dinfo%npah > 0) ndust_total = ndust_total + dinfo%npah

        if (carry_gas_ions) then
            call ensure_update_cache(n_elements, n_elements + 1, ndust_total)
        else
            call ensure_update_cache(n_elements, 1, ndust_total)
        end if
        y_gas => y_gas_cache; y_gas_out => y_gas_out_cache
        y_dust => y_dust_cache; y_dust_out => y_dust_out_cache

        if (carry_gas_ions) then
            do ii = 1, n_elements
                y_gas(ii,1) = nElement(ii) * mass_g(ii)
                y_gas(ii,2:n_elements+1) = nElement(ii) * xelem_ions(ii,:) * mass_g(ii)
            end do
        else
            do ii = 1, n_elements
                y_gas(ii,1) = nElement(ii) * mass_g(ii)
            end do
        end if
        ! Pack y_dust: PAH bins first [1:npah], then grain bins [npah+1:npah+ndust].
        ! This convention is mirrored in the unpack block below and in integrate_dust_ode.
        if (dinfo%ndust > 0 .and. dinfo%npah > 0) then
            y_dust(1:dinfo%npah) = dinfo%rho_pah(1:dinfo%npah)
            y_dust(dinfo%npah+1:dinfo%npah+dinfo%ndust) = dinfo%rho_dust(1:dinfo%ndust)
        else if (dinfo%ndust > 0) then
            y_dust(1:dinfo%ndust) = dinfo%rho_dust(1:dinfo%ndust)
        else if (dinfo%npah > 0) then
            y_dust(1:dinfo%npah) = dinfo%rho_pah(1:dinfo%npah)
        end if

        ! 3. Now we are ready to call the ODE solver to integrate the dust evolution.
        ! Guard against an uninitialised procedure pointer (set by init_dust_solver based
        ! on dust_solver_type); a null dereference here would be a silent SIGSEGV.
        if (.not. associated(dust_solver_step)) then
            print *, "FATAL ERROR: dust_solver_step is not associated. Check dust_solver_type in namelist."
            call clean_stop
        end if
        step_ok = .true.
        h_init_val = dt
        if (last_h_ode > 0.0_dp) h_init_val = last_h_ode
        call integrate_dust_ode(dinfo,dt,y_gas,y_dust,dust_rhs,dust_solver_step,&
                                y_gas_out,y_dust_out,h_init_val,0d0,dt,&
                                debug_flag=dust_debug,step_ok=step_ok,h_last=last_h_ode)

        ! 3b. A NaN/Inf result is never a valid state, and it passes every later check of the caller
        ! (comparisons with NaN are false): the RTZ step then rejected each retry on the ions, halved
        ! ddt to zero and froze the cell. Reject it here, before it reaches nElement and the dust,
        ! and stop if the rates are already non-finite at the starting state (a bug in a rate).
        if (.not. (all(ieee_is_finite(y_dust_out(1:ndust_total))) .and. all(ieee_is_finite(y_gas_out)))) then
            call check_dust_rates_finite(dinfo, dt, y_gas, y_dust(1:ndust_total))
            step_ok = .false.
            return
        end if

        ! 4. Update the dinfo with the new values after the ODE step
        if (carry_gas_ions) then
            do ii = 1, n_elements
                nElement(ii) = y_gas_out(ii,1) / mass_g(ii)
                if (nElement(ii) > 1d-30) then
                    xelem_ions(ii,:) = y_gas_out(ii,2:n_elements+1) / nElement(ii) / mass_g(ii)
                    ! Make sure that xelem_ions add up to 1 for each element
                    sum_check = sum(xelem_ions(ii,:))
                    if (sum_check > 1d-30) then
                        xelem_ions(ii,:) = xelem_ions(ii,:) / sum_check
                    end if
                end if
            end do
        else
            do ii = 1, n_elements
                nElement(ii) = y_gas_out(ii,1) / mass_g(ii)
            end do
        end if

        ! Unpack y_dust_out: PAH bins first [1:npah], grain bins [npah+1:npah+ndust].
        if (dinfo%ndust > 0 .and. dinfo%npah > 0) then
            dinfo%rho_pah(1:dinfo%npah) = y_dust_out(1:dinfo%npah)
            dinfo%rho_dust(1:dinfo%ndust) = y_dust_out(dinfo%npah+1:dinfo%npah+dinfo%ndust)
        else if (dinfo%ndust > 0) then
            dinfo%rho_dust(1:dinfo%ndust) = y_dust_out(1:dinfo%ndust)
        else if (dinfo%npah > 0) then
            dinfo%rho_pah(1:dinfo%npah) = y_dust_out(1:dinfo%npah)
        end if

#ifdef RTZ
        ! Compute dUU: the maximum fractional change in any gas or dust species over this step.
        ! dUU is returned *unscaled* (pure fractional change, dimensionless).
        ! The caller (rtz_cool_step in rtz_cooling_module.f90) multiplies by one_over_x_FRAC
        ! before comparing against the RTZ convergence criterion.

        ! check for gas
        do ii = 1, n_elements
            if (elements(ii)%atomic_number > 0) then
                dUU = max(dUU, abs(y_gas_out(ii,1)-y_gas(ii,1))/(y_gas(ii,1)+y_min))

                ! if not, ion part is zero
                if (carry_gas_ions) then
                    do jj=1,elements(ii)%n_ions
                        dUU = max(dUU, abs(y_gas_out(ii,1+jj)-y_gas(ii,1+jj))/(y_gas(ii,1+jj)+y_min))
                    end do
                end if

            end if
        end do

        ! check for dust
        if (dinfo%ndust > 0 .or. dinfo%npah > 0) then
            do ii = 1, dinfo%npah+dinfo%ndust
                dUU = max(dUU, abs(y_dust_out(ii)-y_dust(ii))/(y_dust(ii)+y_min))
            end do
        end if
#endif

    end subroutine compute_dust_update

    subroutine check_dust_rates_finite(dinfo, dt, y_gas, y_dust)
        ! Called when a dust update came out NaN/Inf. Finite rates at the starting state mean the step
        ! blew up: return, and the caller retries with a smaller dt. Non-finite rates there are a bug in
        ! a rate (a smaller dt cannot help): print the cell and stop. With dust_debug, also name the
        ! process(es) that produced it.
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use dust_rhs_mod, only: dust_rhs
        implicit none
        class(DustChemistryInfo), intent(in) :: dinfo
        real(dp), intent(in) :: dt, y_gas(:,:), y_dust(:)
        real(dp), allocatable :: dydt_gas(:,:), dydt_dust(:,:)
        integer :: i

        allocate(dydt_gas(size(y_gas,1), size(y_gas,2)), dydt_dust(size(y_dust), 2))
        call dust_rhs(dinfo, y_gas, y_dust, dydt_gas, dydt_dust, write_cache=.false.)
        if (all(ieee_is_finite(dydt_gas)) .and. all(ieee_is_finite(dydt_dust))) then
            if (dust_debug) write(*,'(A,ES10.3,A)') ' CALIMA: dust update NaN/Inf over dt = ', dt, &
                                                   ' s with finite rates at its start; step rejected'
            return
        end if

        write(*,'(A)') ' CALIMA ERROR: NaN/Inf dust rates at the start of a dust update'
        write(*,'(A,4ES12.4)') '   T [K], nH, ne [cm-3], rho [g cm-3] = ', dinfo%local_Tk, dinfo%local_nH, &
                                dinfo%local_ne, dinfo%local_rho
        write(*,'(A,3ES12.4)') '   mu, sigma [cm s-1], dx [cm]       = ', dinfo%local_mu, dinfo%local_sigma, dinfo%local_dx
        if (dust_eq_test) write(*,'(A)') '   (dust_eq_test: the grain-grain rates use sigma and L set from nH, not these)'
        write(*,'(A,*(ES12.4))') '   rho bins (PAH first) [g cm-3]     = ', y_dust
        write(*,'(A,*(ES12.4))') '   d/dt, production                  = ', dydt_dust(:,1)
        write(*,'(A,*(ES12.4))') '   d/dt, destruction                 = ', dydt_dust(:,2)
        if (dust_debug) then
            ! each process on its own (the running sum would carry one process's NaN into the rest)
            do i = 1, ndust_processes
                dydt_gas = 0d0; dydt_dust = 0d0
                call dust_processes_list(i)%comp_rate(dinfo, y_gas, y_dust, dydt_gas, dydt_dust)
                if (.not. (all(ieee_is_finite(dydt_gas)) .and. all(ieee_is_finite(dydt_dust)))) &
                    write(*,'(A,A)') '   NaN/Inf from dust process: ', trim(dust_processes_list(i)%name)
            end do
            do i = 1, npah_processes
                dydt_gas = 0d0; dydt_dust = 0d0
                call pah_processes_list(i)%comp_rate(dinfo, y_gas, y_dust, dydt_gas, dydt_dust)
                if (.not. (all(ieee_is_finite(dydt_gas)) .and. all(ieee_is_finite(dydt_dust)))) &
                    write(*,'(A,A)') '   NaN/Inf from PAH process: ', trim(pah_processes_list(i)%name)
            end do
        else
            write(*,'(A)') '   rerun with dust_debug=.true. to see which process'
        end if
        flush(6)   ! clean_stop aborts MPI, which can drop unflushed output
        call clean_stop
    end subroutine check_dust_rates_finite

    subroutine ensure_update_cache(n1, n2, ndust_total)
        integer, intent(in) :: n1, n2, ndust_total
        logical :: need_realloc

        need_realloc = .false.
        if (.not. allocated(y_gas_cache)) then
            need_realloc = .true.
        else if (size(y_gas_cache,1) /= n1 .or. size(y_gas_cache,2) /= n2) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(y_gas_cache)) then
                deallocate(y_gas_cache)
                deallocate(y_gas_out_cache)
            end if
            allocate(y_gas_cache(n1, n2))
            allocate(y_gas_out_cache(n1, n2))
        end if

        need_realloc = .false.
        if (.not. allocated(y_dust_cache)) then
            need_realloc = .true.
        else if (size(y_dust_cache) /= ndust_total) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(y_dust_cache)) then
                deallocate(y_dust_cache)
                deallocate(y_dust_out_cache)
            end if
            allocate(y_dust_cache(ndust_total))
            allocate(y_dust_out_cache(ndust_total))
        end if
    end subroutine ensure_update_cache

end module dust_interface
