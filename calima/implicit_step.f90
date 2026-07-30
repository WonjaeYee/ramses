! Implicit (backward-Euler) ODE solver step for the CALIMA dust/gas chemistry.
!
! At high gas densities (nH > nH_implicit_threshold) the Anninos 10% rule forces
! thousands of substeps because dust-accretion and sputtering timescales become
! shorter than the hydro step and the sequential per-bin updates oscillate around
! the self-consistent equilibrium.
!
! This module provides an implicit_step subroutine that conforms to
! solver_step_interface (defined in ode_interface_mod inside dust_solver.f90).
! It assembles a numerical Jacobian of the full dust+gas RHS and solves the
! backward-Euler system (I - h*J)*delta_y = h*f0 with LU factorisation.
! After the solve, the dust change is used to enforce exact mass conservation
! (same logic as anninos_step), ensuring total gas+dust mass is preserved.
!
! Usage:
!   dust_solver_type = 4  -->  always use implicit_step (for benchmarking)
!   automatic at nH > nH_implicit_threshold for solver_types 1-3 (in dust_interface.f90)

module implicit_step_mod
    use amr_parameters, only: dp
    use dustbin_types,  only: DustChemistryInfo
    use dust_commons,   only: dustbins_props, pahbins_props, carry_gas_ions
    use ode_interface_mod, only: rhs_interface
    use hydro_parameters, only: n_elements, ndust, npah

    implicit none
    private
    public :: implicit_step

    ! ---- Module-level cache arrays (allocated once, reused every call) ----
    real(dp), allocatable, save, target :: impl_dydt_gas(:,:)
    real(dp), allocatable, save, target :: impl_dydt_dust(:)
    real(dp), allocatable, save, target :: impl_y_gas_pert(:,:)
    real(dp), allocatable, save, target :: impl_y_dust_pert(:)
    real(dp), allocatable, save, target :: impl_f0(:)
    real(dp), allocatable, save, target :: impl_f_pert(:)
    real(dp), allocatable, save, target :: impl_A(:,:)
    real(dp), allocatable, save, target :: impl_delta_y(:)
    integer,  allocatable, save, target :: impl_ipiv(:)

    contains

    ! =========================================================================
    !  LU factorisation with partial pivoting (Doolittle, in-place)
    !  On exit: A contains L (below diagonal, unit diagonal not stored) and
    !            U (upper triangle + diagonal).
    !  ipiv(k): row index that was swapped into position k.
    !  info /= 0 if A is singular (pivot A(k,k) == 0).
    ! =========================================================================
    subroutine lu_factor(A, n, ipiv, info)
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: A(n,n)
        integer,  intent(out)   :: ipiv(n)
        integer,  intent(out)   :: info

        integer  :: k, j, imax
        real(dp) :: amax, tmp(n)

        info = 0
        do k = 1, n
            ! Find row with largest |A(j,k)| for j >= k (partial pivoting)
            amax = 0.0_dp
            imax = k
            do j = k, n
                if (abs(A(j,k)) > amax) then
                    amax = abs(A(j,k))
                    imax = j
                end if
            end do
            ipiv(k) = imax
            if (amax == 0.0_dp) then
                info = k
                return
            end if
            ! Swap rows k and imax
            if (imax /= k) then
                tmp     = A(k,:)
                A(k,:)  = A(imax,:)
                A(imax,:) = tmp
            end if
            ! Compute multipliers (L entries below diagonal)
            A(k+1:n, k) = A(k+1:n, k) / A(k,k)
            ! Rank-1 update of the trailing submatrix (U entries)
            do j = k+1, n
                A(j, k+1:n) = A(j, k+1:n) - A(j,k) * A(k, k+1:n)
            end do
        end do
    end subroutine lu_factor

    ! =========================================================================
    !  Solve A*x = b using the LU factorisation produced by lu_factor.
    !  b is overwritten with the solution x.
    ! =========================================================================
    subroutine lu_solve(LU, b, n, ipiv)
        integer,  intent(in)    :: n
        real(dp), intent(in)    :: LU(n,n)
        real(dp), intent(inout) :: b(n)
        integer,  intent(in)    :: ipiv(n)

        integer  :: k
        real(dp) :: tmp

        ! Apply row permutations (reproduce the pivoting done during factorisation)
        do k = 1, n
            tmp       = b(k)
            b(k)      = b(ipiv(k))
            b(ipiv(k)) = tmp
        end do
        ! Forward substitution: L*y = b  (L has unit diagonal)
        do k = 2, n
            b(k) = b(k) - dot_product(LU(k, 1:k-1), b(1:k-1))
        end do
        ! Backward substitution: U*x = y
        do k = n, 1, -1
            if (k < n) b(k) = b(k) - dot_product(LU(k, k+1:n), b(k+1:n))
            b(k) = b(k) / LU(k,k)
        end do
    end subroutine lu_solve

    ! =========================================================================
    !  Ensure cache arrays are allocated with the right sizes.
    ! =========================================================================
    subroutine ensure_implicit_cache(ngas, nvar_gas, ndust_total)
        integer, intent(in) :: ngas, nvar_gas, ndust_total
        integer :: N
        logical :: need

        N = ngas + ndust_total

        need = .not. allocated(impl_dydt_gas)
        if (.not. need) then
            need = (size(impl_dydt_gas,1) /= ngas .or. size(impl_dydt_gas,2) /= nvar_gas)
        end if
        if (need) then
            if (allocated(impl_dydt_gas))   deallocate(impl_dydt_gas)
            if (allocated(impl_y_gas_pert)) deallocate(impl_y_gas_pert)
            allocate(impl_dydt_gas(ngas, nvar_gas))
            allocate(impl_y_gas_pert(ngas, nvar_gas))
        end if

        need = .not. allocated(impl_dydt_dust)
        if (.not. need) need = (size(impl_dydt_dust) /= ndust_total)
        if (need) then
            if (allocated(impl_dydt_dust))   deallocate(impl_dydt_dust)
            if (allocated(impl_y_dust_pert)) deallocate(impl_y_dust_pert)
            allocate(impl_dydt_dust(ndust_total))
            allocate(impl_y_dust_pert(ndust_total))
        end if

        need = .not. allocated(impl_f0)
        if (.not. need) need = (size(impl_f0) /= N)
        if (need) then
            if (allocated(impl_f0))      deallocate(impl_f0)
            if (allocated(impl_f_pert))  deallocate(impl_f_pert)
            if (allocated(impl_A))       deallocate(impl_A)
            if (allocated(impl_delta_y)) deallocate(impl_delta_y)
            if (allocated(impl_ipiv))    deallocate(impl_ipiv)
            allocate(impl_f0(N), impl_f_pert(N))
            allocate(impl_A(N,N))
            allocate(impl_delta_y(N))
            allocate(impl_ipiv(N))
        end if
    end subroutine ensure_implicit_cache

    ! =========================================================================
    !  implicit_step — single backward-Euler step for the CALIMA dust ODE.
    !
    !  Conforms to solver_step_interface (dust_solver.f90, ode_interface_mod).
    !  Always accepts the step (no 10% rule). Suitable for use at high density
    !  where the explicit/Anninos methods require many substeps.
    ! =========================================================================
    subroutine implicit_step(dust_info, y_gas, y_dust, h, rhs, &
                              y_gas_new, y_dust_new, h_new, firstcall, &
                              accepted, break, debug_flag, step_ok_present)

        type(DustChemistryInfo), intent(in)  :: dust_info
        real(dp),                intent(in)  :: y_gas(:,:), y_dust(:)
        real(dp),                intent(inout) :: h
        procedure(rhs_interface)             :: rhs
        real(dp),                intent(out) :: y_gas_new(:,:), y_dust_new(:)
        real(dp),                intent(out) :: h_new
        logical,                 intent(out) :: accepted, break
        logical,                 intent(in)  :: firstcall
        logical,                 intent(in), optional :: debug_flag
        logical,                 intent(in), optional :: step_ok_present

        ! ---- local ----
        integer  :: Ngas, Ndust, N, j, k, ii, kk, e_index, C_index
        integer  :: info
        real(dp) :: eps_j, kmax, delta_dust_j, m_frac
        real(dp) :: mass_init, mass_final, mass_corr
        real(dp), pointer :: dydt_gas(:,:), dydt_dust(:)
        real(dp), pointer :: y_gas_pert(:,:), y_dust_pert(:)
        real(dp), pointer :: f0(:), f_pert(:), A(:,:), delta_y(:)
        integer,  pointer :: ipiv(:)

        Ngas   = n_elements
        Ndust  = size(y_dust)
        N      = Ngas + Ndust

        call ensure_implicit_cache(Ngas, size(y_gas, 2), Ndust)

        dydt_gas   => impl_dydt_gas
        dydt_dust  => impl_dydt_dust
        y_gas_pert => impl_y_gas_pert
        y_dust_pert=> impl_y_dust_pert
        f0         => impl_f0
        f_pert     => impl_f_pert
        A          => impl_A
        delta_y    => impl_delta_y
        ipiv       => impl_ipiv

        break = .false.

        ! ---- 1. Evaluate baseline RHS ----
        dydt_gas   = 0.0_dp
        dydt_dust  = 0.0_dp
        if (firstcall) then
            call rhs(dust_info, y_gas, y_dust, dydt_gas, dydt_dust, kmax, write_cache=.true.)
            if (kmax == 0.0_dp) then
                break     = .true.
                accepted  = .true.
                h_new     = h
                y_gas_new = y_gas
                y_dust_new= y_dust
                return
            end if
        else
            call rhs(dust_info, y_gas, y_dust, dydt_gas, dydt_dust, write_cache=.true.)
        end if

        ! Pack baseline RHS into flat vector f0
        ! f0(1:Ngas)        = dydt_gas(:,1)   (total gas mass rates)
        ! f0(Ngas+1:N)      = dydt_dust        (dust+PAH mass rates)
        f0(1:Ngas) = dydt_gas(1:Ngas, 1)
        f0(Ngas+1:N) = dydt_dust(1:Ndust)

        ! ---- 2. Build numerical Jacobian ----
        ! Perturb each component of the state vector [y_gas(:,1), y_dust(:)]
        ! and evaluate the RHS to get the Jacobian column.
        ! Note: y_gas(:,2+) ion fractions are NOT perturbed (dust ODE does not
        ! couple to them — only total gas mass per element changes).
        do j = 1, N
            y_gas_pert  = y_gas
            y_dust_pert = y_dust

            if (j <= Ngas) then
                eps_j = sqrt(epsilon(1.0_dp)) * max(abs(y_gas(j,1)), 1.0d-30)
                y_gas_pert(j, 1) = y_gas(j,1) + eps_j
            else
                eps_j = sqrt(epsilon(1.0_dp)) * max(abs(y_dust(j-Ngas)), 1.0d-30)
                y_dust_pert(j-Ngas) = y_dust(j-Ngas) + eps_j
            end if

            dydt_gas  = 0.0_dp
            dydt_dust = 0.0_dp
            call rhs(dust_info, y_gas_pert, y_dust_pert, dydt_gas, dydt_dust)

            f_pert(1:Ngas)   = dydt_gas(1:Ngas, 1)
            f_pert(Ngas+1:N) = dydt_dust(1:Ndust)

            A(:, j) = (f_pert - f0) / eps_j
        end do

        ! ---- 3. Build (I - h*J) in place ----
        do k = 1, N
            A(k, :) = -h * A(k, :)
            A(k, k) = A(k, k) + 1.0_dp
        end do

        ! ---- 4. Solve (I - h*J)*delta_y = h*f0 ----
        delta_y = h * f0
        call lu_factor(A, N, ipiv, info)
        if (info /= 0) then
            ! Singular matrix — fall back to explicit Euler for this step
            delta_y = h * f0
        else
            call lu_solve(A, delta_y, N, ipiv)
        end if

        ! ---- 5. Apply update ----
        ! Initialise output to current state
        y_gas_new  = y_gas
        y_dust_new = y_dust

        ! Update dust/PAH bins from delta_y, then enforce exact mass conservation
        ! by deriving the gas change from the dust change (same as anninos_step).
        do j = 1, Ndust
            y_dust_new(j) = max(y_dust(j) + delta_y(Ngas + j), 0.0_dp)

            delta_dust_j = y_dust_new(j) - y_dust(j)
            if (delta_dust_j == 0.0_dp) cycle

            if (j <= dust_info%npah) then
                ! PAH bin: all mass is carbon (element index 6)
                C_index = 6
                y_gas_new(C_index, 1) = y_gas_new(C_index, 1) - delta_dust_j
            else
                ! Dust bin: distribute mass loss over constituent elements
                ii = j - dust_info%npah
                associate(bin => dustbins_props(ii))
                    do kk = 1, bin%nelements
                        e_index = bin%el_index(kk)
                        m_frac  = bin%el_mfractions(kk)
                        y_gas_new(e_index, 1) = y_gas_new(e_index, 1) - delta_dust_j * m_frac
                    end do
                end associate
            end if
        end do

        ! Enforce non-negativity in gas (can go negative if a bin grew very rapidly)
        y_gas_new(:, 1) = max(y_gas_new(:, 1), 0.0_dp)

        ! Verify total mass conservation; apply a proportional correction if needed
        mass_init  = sum(y_gas(:, 1)) + sum(y_dust)
        mass_final = sum(y_gas_new(:, 1)) + sum(y_dust_new)
        if (mass_init > 0.0_dp) then
            mass_corr = mass_init - mass_final
            if (abs(mass_corr) / mass_init > 1.0d-10) then
                ! Distribute the discrepancy back into gas proportionally
                where (y_gas_new(:, 1) > 0.0_dp)
                    y_gas_new(:, 1) = y_gas_new(:, 1) * (1.0_dp + mass_corr / mass_final)
                end where
            end if
        end if

        ! ---- 6. Always accept ----
        accepted = .true.
        h_new    = h

    end subroutine implicit_step

end module implicit_step_mod
