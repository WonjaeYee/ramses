module ode_interface_mod

    use amr_parameters, only: dp
    use hydro_parameters, only: n_elements, ndust, npah
    use dustbin_types, only: DustChemistryInfo

    implicit none
    private
    public :: ZERO, ONE, TWO, HALF, SIXTH, rhs_interface, solver_step_interface, dust_solver_step

    real(dp), parameter :: ZERO  = 0.0_dp
    real(dp), parameter :: ONE   = 1.0_dp
    real(dp), parameter :: TWO   = 2.0_dp
    real(dp), parameter :: HALF  = 0.5_dp
    real(dp), parameter :: SIXTH = 1.0_dp / 6.0_dp

    abstract interface
        subroutine rhs_interface(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax,debug_flag)
            import :: dp, DustChemistryInfo
            implicit none
            type(DustChemistryInfo), intent(in) :: dust_info
            real(dp), intent(in) :: y_gas(:,:), y_dust(:)
            real(dp), intent(inout) :: dydt_gas(:,:), dydt_dust(:)
            real(dp), intent(inout), optional :: kmax
            logical, intent(in), optional :: debug_flag
        end subroutine rhs_interface
    end interface

    abstract interface
        subroutine solver_step_interface(dust_info,y_gas,y_dust,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break,debug_flag)
            import :: dp, DustChemistryInfo, rhs_interface
            implicit none
            type(DustChemistryInfo), intent(in) :: dust_info
            real(dp), intent(in) :: y_gas(:,:), y_dust(:)
            real(dp), intent(inout) :: h
            procedure(rhs_interface) :: rhs
            real(dp), intent(out) :: y_gas_new(:,:), y_dust_new(:)
            real(dp), intent(out) :: h_new
            logical, intent(out) :: accepted
            logical, intent(out) :: break
            logical, intent(in) :: firstcall
            logical, intent(in), optional :: debug_flag
        end subroutine solver_step_interface
    end interface

    procedure(solver_step_interface), pointer :: dust_solver_step => null()

end module ode_interface_mod

module dust_rhs_mod
    use amr_parameters, only: dp
    use dust_commons
    use dust_rates

    implicit none
    private
    public :: dust_rhs, print_last_process_kmax
    public :: reset_timestep_reduction_counters, register_timestep_reduction_cause
    public :: print_timestep_reduction_counters
    public :: last_dydt_dust_per_proc, last_dydt_pah_per_proc

    real(dp), allocatable, save :: last_kmax_dust(:)
    real(dp), allocatable, save :: last_kmax_pah(:)
    ! Per-process dydt contributions from the most recent dust_rhs call [g cm-3 s-1].
    ! Indexed as last_dydt_dust_per_proc(ispecies, iprocess).
    real(dp), allocatable, save :: last_dydt_dust_per_proc(:,:)
    real(dp), allocatable, save :: last_dydt_pah_per_proc(:,:)
    real(dp), allocatable, save :: dydt_dust_before_cache(:)

    contains

    subroutine ensure_kmax_storage
        implicit none
        logical, save :: first_call = .true.

        if (.not. first_call) return

        first_call = .false.

        if (allocated(dydt_dust_before_cache)) then
            if (size(dydt_dust_before_cache) /= ndust+npah) then
                deallocate(dydt_dust_before_cache)
            end if
        end if
        if (.not. allocated(dydt_dust_before_cache)) then
            allocate(dydt_dust_before_cache(ndust+npah))
        end if

        if (allocated(last_kmax_dust)) then
            if (size(last_kmax_dust) /= ndust_processes) then
                deallocate(last_kmax_dust)
            end if
        end if
        if (.not. allocated(last_kmax_dust)) then
            allocate(last_kmax_dust(max(0, ndust_processes)))
        end if

        if (allocated(last_kmax_pah)) then
            if (size(last_kmax_pah) /= npah_processes) then
                deallocate(last_kmax_pah)
            end if
        end if
        if (.not. allocated(last_kmax_pah)) then
            allocate(last_kmax_pah(max(0, npah_processes)))
        end if

        if (dust_log) then
            if (allocated(ode_reduction_count_dust)) then
                if (size(ode_reduction_count_dust) /= ndust_processes) then
                    deallocate(ode_reduction_count_dust)
                end if
            end if
            if (.not. allocated(ode_reduction_count_dust)) then
                allocate(ode_reduction_count_dust(max(0, ndust_processes)))
                ode_reduction_count_dust(:) = 0_8
            end if

            if (allocated(ode_reduction_count_pah)) then
                if (size(ode_reduction_count_pah) /= npah_processes) then
                    deallocate(ode_reduction_count_pah)
                end if
            end if
            if (.not. allocated(ode_reduction_count_pah)) then
                allocate(ode_reduction_count_pah(max(0, npah_processes)))
                ode_reduction_count_pah(:) = 0_8
            end if

            if (allocated(last_dydt_dust_per_proc)) then
                if (size(last_dydt_dust_per_proc, 1) /= ndust+npah .or. &
                    size(last_dydt_dust_per_proc, 2) /= max(1, ndust_processes)) then
                    deallocate(last_dydt_dust_per_proc)
                end if
            end if
            if (.not. allocated(last_dydt_dust_per_proc)) then
                allocate(last_dydt_dust_per_proc(ndust+npah, max(1, ndust_processes)))
                last_dydt_dust_per_proc(:,:) = 0.0_dp
            end if

            if (allocated(last_dydt_pah_per_proc)) then
                if (size(last_dydt_pah_per_proc, 1) /= ndust+npah .or. &
                    size(last_dydt_pah_per_proc, 2) /= max(1, npah_processes)) then
                    deallocate(last_dydt_pah_per_proc)
                end if
            end if
            if (.not. allocated(last_dydt_pah_per_proc)) then
                allocate(last_dydt_pah_per_proc(ndust+npah, max(1, npah_processes)))
                last_dydt_pah_per_proc(:,:) = 0.0_dp
            end if
        end if
    end subroutine ensure_kmax_storage

    subroutine reset_timestep_reduction_counters
        implicit none

        call ensure_kmax_storage
        if (allocated(ode_reduction_count_dust)) ode_reduction_count_dust(:) = 0_8
        if (allocated(ode_reduction_count_pah)) ode_reduction_count_pah(:) = 0_8
    end subroutine reset_timestep_reduction_counters

    subroutine register_timestep_reduction_cause
        implicit none

        integer :: i, max_i
        logical :: is_dust
        real(dp) :: vmax, val

        call ensure_kmax_storage

        vmax = -1d0
        max_i = 0
        is_dust = .true.

        do i = 1, ndust_processes
            val = last_kmax_dust(i)
            if (val > vmax) then
                vmax = val
                max_i = i
                is_dust = .true.
            end if
        end do

        do i = 1, npah_processes
            val = last_kmax_pah(i)
            if (val > vmax) then
                vmax = val
                max_i = i
                is_dust = .false.
            end if
        end do

        if (max_i <= 0 .or. vmax <= 0d0) return

        if (is_dust) then
            ode_reduction_count_dust(max_i) = ode_reduction_count_dust(max_i) + 1_8
        else
            ode_reduction_count_pah(max_i) = ode_reduction_count_pah(max_i) + 1_8
        end if
    end subroutine register_timestep_reduction_cause

    subroutine print_last_process_kmax
        implicit none
        integer :: i

        call ensure_kmax_storage

        print *, 'ODE diagnostics: per-process kmax [s^-1] from last RHS evaluation'
        if (ndust_processes > 0) then
            do i = 1, ndust_processes
                print *, '  dust process ', i, ' (', trim(dust_processes_list(i)%name), '): ', last_kmax_dust(i)
            end do
        else
            print *, '  No dust processes active.'
        end if

        if (npah_processes > 0) then
            do i = 1, npah_processes
                print *, '  pah process  ', i, ' (', trim(pah_processes_list(i)%name), '): ', last_kmax_pah(i)
            end do
        else
            print *, '  No PAH processes active.'
        end if
    end subroutine print_last_process_kmax

    subroutine print_timestep_reduction_counters
        implicit none
        integer :: i
        integer*8 :: total_reductions

        call ensure_kmax_storage

        total_reductions = 0_8
        if (allocated(ode_reduction_count_dust)) total_reductions = total_reductions + sum(ode_reduction_count_dust)
        if (allocated(ode_reduction_count_pah)) total_reductions = total_reductions + sum(ode_reduction_count_pah)

        print *, 'ODE diagnostics: timestep-reduction attributions by dominant process'
        print *, '  Total attributed reductions       = ', total_reductions

        if (ndust_processes > 0) then
            do i = 1, ndust_processes
                print *, '  dust process ', i, ' (', trim(dust_processes_list(i)%name), '): ', ode_reduction_count_dust(i)
            end do
        else
            print *, '  No dust processes active.'
        end if

        if (npah_processes > 0) then
            do i = 1, npah_processes
                print *, '  pah process  ', i, ' (', trim(pah_processes_list(i)%name), '): ', ode_reduction_count_pah(i)
            end do
        else
            print *, '  No PAH processes active.'
        end if
    end subroutine print_timestep_reduction_counters

    subroutine dust_rhs(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax,debug_flag)
        ! Compute the right-hand side of the ODE system for the dust chemistry. This function will be called
        ! by the ODE solver to compute the time derivatives of the gas and dust abundances.
        ! dust_info  --> DustChemistryInfo type containing physical parameters and metadata.
        ! y_gas      --> 2D array with the gas phase abundances [g cm-3]
        ! y_dust     --> 1D array with the dust phase abundances [g cm-3]
        ! dydt_gas   <--> 2D array with the time derivative of the gas phase abundances [g cm-3 s-1]
        ! dydt_dust  <--> 1D array with the time derivative of the dust phase abundances [g cm-3 s-1]
        ! kmax       <--> Maximum allowed rate for the process [s-1] (optional output)
        ! debug_flag --> Optional logical flag to enable verbose debugging.

        implicit none
        ! ---- Input/Output variables ----
        type(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: dydt_gas(:,:), dydt_dust(:)
        real(dp), intent(inout), optional :: kmax
        logical, intent(in), optional :: debug_flag

        ! ---- Local variables ----
        integer :: i
        real(dp) :: process_kmax

        ! 1. Initialize the time derivatives to zero
        dydt_gas(:,:) = 0.0_dp
        dydt_dust(:) = 0.0_dp

        call ensure_kmax_storage
        if (allocated(last_kmax_dust)) last_kmax_dust(:) = 0.0_dp
        if (allocated(last_kmax_pah)) last_kmax_pah(:) = 0.0_dp
        if (dust_log) then
            if (allocated(last_dydt_dust_per_proc)) last_dydt_dust_per_proc(:,:) = 0.0_dp
            if (allocated(last_dydt_pah_per_proc))  last_dydt_pah_per_proc(:,:)  = 0.0_dp
        end if

        ! 2. Loop over the dust processes and compute their contribution to the time derivatives
        if (present(kmax)) kmax = 0.0_dp

        do i = 1, ndust_processes
            process_kmax = 0.0_dp
            if (dust_log) then
                dydt_dust_before_cache(:) = dydt_dust(:)
                call dust_processes_list(i)%comp_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,process_kmax)
                if (allocated(last_dydt_dust_per_proc)) &
                    last_dydt_dust_per_proc(:, i) = dydt_dust(:) - dydt_dust_before_cache(:)
            else
                call dust_processes_list(i)%comp_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,process_kmax)
            end if
            last_kmax_dust(i) = process_kmax
            if (present(kmax)) kmax = max(kmax, process_kmax)
        end do
        do i = 1, npah_processes
            process_kmax = 0.0_dp
            if (dust_log) then
                dydt_dust_before_cache(:) = dydt_dust(:)
                call pah_processes_list(i)%comp_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,process_kmax)
                if (allocated(last_dydt_pah_per_proc)) &
                    last_dydt_pah_per_proc(:, i) = dydt_dust(:) - dydt_dust_before_cache(:)
            else
                call pah_processes_list(i)%comp_rate(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,process_kmax)
            end if
            last_kmax_pah(i) = process_kmax
            if (present(kmax)) kmax = max(kmax, process_kmax)
        end do
    end subroutine dust_rhs

end module dust_rhs_mod

module rk4_mod
    use amr_parameters, only: dp
    use dustbin_types, only: DustChemistryInfo
    use dust_commons, only: errmax
    use ode_interface_mod, only: rhs_interface, HALF, TWO, SIXTH

    implicit none
    private
    public :: rk4_step

    ! Cache for intermediate RK4 stages to avoid automatic array allocations
    real(dp), allocatable, save, target :: k1_gas_cache(:,:), k2_gas_cache(:,:), k3_gas_cache(:,:), k4_gas_cache(:,:)
    real(dp), allocatable, save, target :: k1_dust_cache(:), k2_dust_cache(:), k3_dust_cache(:), k4_dust_cache(:)
    real(dp), allocatable, save, target :: y_gas_temp_cache(:,:), y_dust_temp_cache(:)
    real(dp), allocatable, save, target :: error_gas_cache(:,:), error_dust_cache(:)

    contains

    subroutine ensure_rk4_cache(ngas_species, nvar_gas, ndust_total)
        ! Ensure that the intermediate RK4 cache arrays are allocated with the correct dimensions.
        ! ngas_species --> Number of gas phase chemical elements/species.
        ! nvar_gas     --> Number of variables per gas phase element.
        ! ndust_total  --> Total number of dust bins (dust + PAH).
        integer, intent(in) :: ngas_species, nvar_gas, ndust_total
        logical :: need_realloc

        need_realloc = .false.
        if (.not. allocated(k1_gas_cache)) then
            need_realloc = .true.
        else if (size(k1_gas_cache,1) /= ngas_species .or. size(k1_gas_cache,2) /= nvar_gas) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(k1_gas_cache)) deallocate(k1_gas_cache, k2_gas_cache, k3_gas_cache, k4_gas_cache)
            if (allocated(y_gas_temp_cache)) deallocate(y_gas_temp_cache)
            if (allocated(error_gas_cache)) deallocate(error_gas_cache)
            allocate(k1_gas_cache(ngas_species, nvar_gas), k2_gas_cache(ngas_species, nvar_gas), &
                     k3_gas_cache(ngas_species, nvar_gas), k4_gas_cache(ngas_species, nvar_gas))
            allocate(y_gas_temp_cache(ngas_species, nvar_gas))
            allocate(error_gas_cache(ngas_species, nvar_gas))
        end if

        need_realloc = .false.
        if (.not. allocated(k1_dust_cache)) then
            need_realloc = .true.
        else if (size(k1_dust_cache) /= ndust_total) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(k1_dust_cache)) deallocate(k1_dust_cache, k2_dust_cache, k3_dust_cache, k4_dust_cache)
            if (allocated(y_dust_temp_cache)) deallocate(y_dust_temp_cache)
            if (allocated(error_dust_cache)) deallocate(error_dust_cache)
            allocate(k1_dust_cache(ndust_total), k2_dust_cache(ndust_total), &
                     k3_dust_cache(ndust_total), k4_dust_cache(ndust_total))
            allocate(y_dust_temp_cache(ndust_total))
            allocate(error_dust_cache(ndust_total))
        end if
    end subroutine ensure_rk4_cache

    subroutine rk4_raw(dust_info,y_gas,y_dust,h,rhs,y_gas_new,y_dust_new,first_call,break,debug_flag)
        ! Perform a single step of the classical 4th-order Runge-Kutta method (RK4) to solve the ODE system.
        ! This is a "raw" implementation that does not include any adaptive time stepping or error control.
        ! dust_info  --> DustChemistryInfo type containing physical parameters and metadata.
        ! y_gas      --> 2D array with the gas phase abundances [g cm-3]
        ! y_dust     --> 1D array with the dust phase abundances [g cm-3]
        ! h          --> Time step size [s]
        ! rhs        --> Procedure pointer to the RHS function that computes derivatives.
        ! y_gas_new  <-- 2D array with the updated gas phase abundances [g cm-3]
        ! y_dust_new <-- 1D array with the updated dust phase abundances [g cm-3]
        ! first_call --> Logical flag indicating whether this is the first call to the subroutine
        ! break      <-- Logical flag indicating whether the integration should be stopped (e.g., if there are no active dust processes)
        ! debug_flag --> Optional logical flag to enable verbose debugging.

        implicit none
        ! ---- Input/Output variables ----
        type(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: h
        procedure(rhs_interface) :: rhs
        real(dp), intent(out) :: y_gas_new(:,:), y_dust_new(:)
        logical, intent(in) :: first_call
        logical, intent(out) :: break
        logical, intent(in), optional :: debug_flag

        ! ---- Local variables ----
        real(dp), pointer :: k1_gas(:,:), k2_gas(:,:), k3_gas(:,:), k4_gas(:,:)
        real(dp), pointer :: k1_dust(:), k2_dust(:), k3_dust(:), k4_dust(:)
        real(dp) :: kmax, h_local

        call ensure_rk4_cache(size(y_gas,1), size(y_gas,2), size(y_dust))
        k1_gas => k1_gas_cache; k2_gas => k2_gas_cache; k3_gas => k3_gas_cache; k4_gas => k4_gas_cache
        k1_dust => k1_dust_cache; k2_dust => k2_dust_cache; k3_dust => k3_dust_cache; k4_dust => k4_dust_cache

        break = .false.

        ! 1. Perform the four RK4 stages
        if (first_call) then
            ! On the first call, we compute kmax to get an estimate of the necessary time step size for stability.
            if (present(debug_flag)) then
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust,kmax,debug_flag=debug_flag)
            else
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust,kmax)
            end if
            if (kmax.eq.0d0) then
                ! If kmax is zero it means there are no active dust
                ! processes, so there is no point in doing a dust integration.
                break = .true.
                return
            end if
        else
            if (present(debug_flag)) then
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust,debug_flag=debug_flag)
            else
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust)
            end if
        end if
      
        ! 2. Use the provided kmax to compute a guess of the neccessary
        ! time step size for stability, but never increase the time step 
        ! beyond the provided h
        if (first_call) then
            h_local = min(1d0 / kmax,h)
        else
            h_local = h
        end if
        if (present(debug_flag)) then
            call rhs(dust_info,y_gas+h_local*HALF*k1_gas,y_dust+h_local*HALF*k1_dust,k2_gas,k2_dust,debug_flag=debug_flag)
            call rhs(dust_info,y_gas+h_local*HALF*k2_gas,y_dust+h_local*HALF*k2_dust,k3_gas,k3_dust,debug_flag=debug_flag)
            call rhs(dust_info,y_gas+h_local*k3_gas,y_dust+h_local*k3_dust,k4_gas,k4_dust,debug_flag=debug_flag)
        else
            call rhs(dust_info,y_gas+h_local*HALF*k1_gas,y_dust+h_local*HALF*k1_dust,k2_gas,k2_dust)
            call rhs(dust_info,y_gas+h_local*HALF*k2_gas,y_dust+h_local*HALF*k2_dust,k3_gas,k3_dust)
            call rhs(dust_info,y_gas+h_local*k3_gas,y_dust+h_local*k3_dust,k4_gas,k4_dust)
        end if

        ! 3. Combine the stages to compute the new solution
        y_gas_new(:,:) = y_gas(:,:) + (h_local * SIXTH) * (k1_gas(:,:) + TWO*k2_gas(:,:) + TWO*k3_gas(:,:) + k4_gas(:,:))
        y_dust_new(:) = y_dust(:) + (h_local * SIXTH) * (k1_dust(:) + TWO*k2_dust(:) + TWO*k3_dust(:) + k4_dust(:))

        ! 4. Update the time step size for the next step (this will be used by the adaptive RK4 wrapper)
        h = h_local
    end subroutine rk4_raw

    subroutine rk4_step(dust_info,y_gas,y_dust,h,rhs,y_gas_new,y_dust_new,h_new,first_call,accepted,break,debug_flag)
        ! Perform a single step of the RK4 method with adaptive time stepping and error control.
        ! This is a wrapper around the rk4_raw subroutine that includes logic for adjusting the time step size
        ! based on the estimated error of the solution.
        ! dust_info  --> DustChemistryInfo type containing physical parameters and metadata.
        ! y_gas      --> 2D array with the gas phase abundances [g cm-3]
        ! y_dust     --> 1D array with the dust phase abundances [g cm-3]
        ! h          --> Time step size [s]
        ! rhs        --> Procedure pointer to the RHS function that computes derivatives.
        ! y_gas_new  <-- 2D array with the updated gas phase abundances [g cm-3]
        ! y_dust_new <-- 1D array with the updated dust phase abundances [g cm-3]
        ! h_new      <-- Updated time step size [s]
        ! first_call --> Logical flag indicating whether this is the first call to the subroutine
        ! accepted   <-- Logical flag indicating whether the step was accepted or not (for adaptive time stepping)
        ! break      <-- Logical flag indicating whether the integration should be stopped (e.g., if there are no active dust processes)
        ! debug_flag --> Optional logical flag to enable verbose debugging.

        implicit none
        ! ---- Input/Output variables ----
        type(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: h
        procedure(rhs_interface) :: rhs
        real(dp), intent(out) :: y_gas_new(:,:), y_dust_new(:)
        real(dp), intent(out) :: h_new
        logical, intent(out) :: accepted,break
        logical, intent(in) :: first_call
        logical, intent(in), optional :: debug_flag

        ! ---- Local variables ----
        real(dp), pointer :: y_gas_temp(:,:), y_dust_temp(:)
        real(dp), pointer :: error_gas(:,:), error_dust(:)
        real(dp) :: max_error, scale

        call ensure_rk4_cache(size(y_gas,1), size(y_gas,2), size(y_dust))
        y_gas_temp => y_gas_temp_cache; y_dust_temp => y_dust_temp_cache
        error_gas => error_gas_cache; error_dust => error_dust_cache

        ! 1. Perform a raw RK4 step to get the new solution
        if (present(debug_flag)) then
            call rk4_raw(dust_info,y_gas,y_dust,h,rhs,y_gas_temp,y_dust_temp,first_call,break,debug_flag=debug_flag)
        else
            call rk4_raw(dust_info,y_gas,y_dust,h,rhs,y_gas_temp,y_dust_temp,first_call,break)
        end if
        if (break) then
            ! If break is true, it means there are no active dust processes, so we can skip the rest of the logic.
            accepted = .true.
            h_new = h
            return
        end if

        ! 2. Max relative change w.r.t. old state
        error_gas(:,:) = abs(y_gas_temp(:,:) - y_gas(:,:)) / max(abs(y_gas(:,:)), 1.0d-40)
        error_dust(:) = abs(y_dust_temp(:) - y_dust(:)) / max(abs(y_dust(:)), 1.0d-40)
        max_error = maxval(error_gas(:,:))
        max_error = max(max_error, maxval(error_dust(:)))

        accepted = (max_error <= errmax)

        ! 3. Update the step size based on the error
        scale = 0.9d0 * (errmax / max(max_error,1.0d-10))
        h_new = h * min(2.0_dp, max(0.1_dp, scale))  ! Limit the change in step size

        ! 4. If the step is accepted, update the solution; otherwise, keep the old solution
        if (accepted) then
            y_gas_new(:,:) = y_gas_temp(:,:)
            y_dust_new(:) = y_dust_temp(:)
        end if
    end subroutine rk4_step
end module rk4_mod

module anninos_mod
    use amr_parameters, only: dp
    use dustbin_types, only: DustChemistryInfo
    use dust_commons, only: errmax, dustbins_props, pahbins_props
    use ode_interface_mod, only: rhs_interface

    implicit none
    private
    public :: anninos_step

    ! Cache to avoid automatic array allocations
    real(dp), allocatable, save, target :: dydt_gas_cache(:,:), dydt_dust_cache(:)
    real(dp), allocatable, save, target :: error_gas_cache(:,:), error_dust_cache(:)

    contains

    subroutine ensure_anninos_cache(ngas_species, nvar_gas, ndust_total)
        ! Ensure that the intermediate Anninos cache arrays are allocated with the correct dimensions.
        ! ngas_species --> Number of gas phase chemical elements/species.
        ! nvar_gas     --> Number of variables per gas phase element.
        ! ndust_total  --> Total number of dust bins (dust + PAH).
        integer, intent(in) :: ngas_species, nvar_gas, ndust_total
        logical :: need_realloc

        need_realloc = .false.
        if (.not. allocated(dydt_gas_cache)) then
            need_realloc = .true.
        else if (size(dydt_gas_cache,1) /= ngas_species .or. size(dydt_gas_cache,2) /= nvar_gas) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(dydt_gas_cache)) deallocate(dydt_gas_cache)
            if (allocated(error_gas_cache)) deallocate(error_gas_cache)
            allocate(dydt_gas_cache(ngas_species, nvar_gas))
            allocate(error_gas_cache(ngas_species, nvar_gas))
        end if

        need_realloc = .false.
        if (.not. allocated(dydt_dust_cache)) then
            need_realloc = .true.
        else if (size(dydt_dust_cache) /= ndust_total) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(dydt_dust_cache)) deallocate(dydt_dust_cache)
            if (allocated(error_dust_cache)) deallocate(error_dust_cache)
            allocate(dydt_dust_cache(ndust_total))
            allocate(error_dust_cache(ndust_total))
        end if
    end subroutine ensure_anninos_cache

    subroutine anninos_step(dust_info,y_gas,y_dust,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break,debug_flag)
        ! Perform a single integration step using the quasi-implicit Anninos et al. (1997) method.
        ! dust_info  --> DustChemistryInfo type with physical parameters.
        ! y_gas      --> 2D array with the gas phase abundances at the start of the step [g cm-3].
        ! y_dust     --> 1D array with the dust phase abundances at the start of the step [g cm-3].
        ! h          --> Time step size [s].
        ! rhs        --> Procedure pointer to the RHS function that computes derivatives.
        ! y_gas_new  <-- 2D array with the updated gas phase abundances [g cm-3].
        ! y_dust_new <-- 1D array with the updated dust phase abundances [g cm-3].
        ! h_new      <-- Recommended time step size for the next step [s].
        ! firstcall  --> Logical flag indicating whether this is the first call in this ODE integration.
        ! accepted   <-- Logical flag indicating whether the step satisfied error limits (10% rule).
        ! break      <-- Logical flag indicating whether integration can be stopped early.
        ! debug_flag --> Optional logical flag to enable verbose debugging.
        implicit none
        ! ---- Input/Output variables ----
        type(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: h
        procedure(rhs_interface) :: rhs
        real(dp), intent(out) :: y_gas_new(:,:), y_dust_new(:)
        real(dp), intent(out) :: h_new
        logical, intent(out) :: accepted, break
        logical, intent(in) :: firstcall
        logical, intent(in), optional :: debug_flag

        ! ---- Local variables ----
        real(dp), pointer :: dydt_gas(:,:), dydt_dust(:)
        real(dp), pointer :: error_gas(:,:), error_dust(:)
        real(dp) :: kmax, scale, max_error
        real(dp) :: yj, fj, Cj, Dj, y_eq, delta_y_dust, m_frac
        integer :: j, ii, kk, e_index, C_index

        call ensure_anninos_cache(size(y_gas,1), size(y_gas,2), size(y_dust))
        dydt_gas => dydt_gas_cache; dydt_dust => dydt_dust_cache
        error_gas => error_gas_cache; error_dust => error_dust_cache

        break = .false.

        ! 1. Evaluate RHS
        if (firstcall) then
            if (present(debug_flag)) then
                call rhs(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax,debug_flag=debug_flag)
            else
                call rhs(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,kmax)
            end if
            if (kmax == 0.0_dp) then
                break = .true.
                accepted = .true.
                h_new = h
                return
            end if
        else
            if (present(debug_flag)) then
                call rhs(dust_info,y_gas,y_dust,dydt_gas,dydt_dust,debug_flag=debug_flag)
            else
                call rhs(dust_info,y_gas,y_dust,dydt_gas,dydt_dust)
            end if
        end if

        ! Initialize output arrays
        y_gas_new(:,:) = y_gas(:,:)

        ! 2. Update dust/PAH bins using Anninos method
        do j = 1, size(y_dust)
            yj = y_dust(j)
            fj = dydt_dust(j)

            if (abs(fj) * h < 1.0d-12 * yj) then
                y_dust_new(j) = yj
            else if (abs(fj) * h < 1.0d-2 * yj) then
                y_dust_new(j) = yj + fj * h
            else
                ! Decompose derivative into creation Cj and destruction Dj
                if (fj >= 0.0_dp) then
                    Cj = fj
                    Dj = 0.0_dp
                else
                    Cj = 0.0_dp
                    Dj = -fj / max(yj, 1.0d-30)
                end if

                ! Anninos et al. (1997) quasi-implicit update
                if (Dj * h < 1.0d-6) then
                    y_dust_new(j) = (yj + Cj * h) / (1.0d0 + Dj * h)
                else
                    y_eq = Cj / Dj
                    y_dust_new(j) = y_eq + (yj - y_eq) * exp(-Dj * h)
                end if
            end if

            ! Enforce non-negativity
            y_dust_new(j) = max(y_dust_new(j), 0.0_dp)

            ! 3. Symmetrically update gas phase elements to conserve mass
            delta_y_dust = y_dust_new(j) - yj
            if (j <= dust_info%npah) then
                ! PAH bin: only carbon (C_index = 6)
                C_index = 6
                y_gas_new(C_index, 1) = y_gas_new(C_index, 1) - delta_y_dust
            else
                ! Dust bin
                ii = j - dust_info%npah
                associate(bin => dustbins_props(ii))
                    do kk = 1, bin%nelements
                        e_index = bin%el_index(kk)
                        m_frac = bin%el_mfractions(kk)
                        y_gas_new(e_index, 1) = y_gas_new(e_index, 1) - delta_y_dust * m_frac
                    end do
                end associate
            end if
        end do

        ! Check for negative values in the new state
        if (any(y_gas_new < 0.0_dp) .or. any(y_dust_new < 0.0_dp)) then
            accepted = .false.
            h_new = h * 0.5_dp
        else
            ! Compute relative error for step control (10% rule)
            error_gas(:,:) = abs(y_gas_new(:,:) - y_gas(:,:)) / max(abs(y_gas(:,:)), 1.0d-40)
            error_dust(:) = abs(y_dust_new(:) - y_dust(:)) / max(abs(y_dust(:)), 1.0d-40)
            max_error = maxval(error_gas(:,:))
            max_error = max(max_error, maxval(error_dust(:)))

            accepted = (max_error <= errmax)
            scale = 0.9d0 * (errmax / max(max_error, 1.0d-10))
            h_new = h * min(2.0_dp, max(0.1_dp, scale))
        end if

    end subroutine anninos_step
end module anninos_mod

module rk54_mod
    use amr_parameters, only: dp
    use dustbin_types, only: DustChemistryInfo
    use dust_commons, only: errmax
    use ode_interface_mod, only: rhs_interface

    implicit none
    private
    public :: rk54_step

    ! Cache for intermediate RK54 stages to avoid automatic array allocations
    real(dp), allocatable, save, target :: k1_gas_cache(:,:), k2_gas_cache(:,:), k3_gas_cache(:,:), &
                                           k4_gas_cache(:,:), k5_gas_cache(:,:), k6_gas_cache(:,:)
    real(dp), allocatable, save, target :: k1_dust_cache(:), k2_dust_cache(:), k3_dust_cache(:), &
                                           k4_dust_cache(:), k5_dust_cache(:), k6_dust_cache(:)
    real(dp), allocatable, save, target :: y_gas_temp_cache(:,:), y_dust_temp_cache(:)
    real(dp), allocatable, save, target :: error_gas_cache(:,:), error_dust_cache(:)

    contains

    subroutine ensure_rk54_cache(ngas_species, nvar_gas, ndust_total)
        ! Ensure that the intermediate RK54 cache arrays are allocated with the correct dimensions.
        ! ngas_species --> Number of gas phase chemical elements/species.
        ! nvar_gas     --> Number of variables per gas phase element.
        ! ndust_total  --> Total number of dust bins (dust + PAH).
        integer, intent(in) :: ngas_species, nvar_gas, ndust_total
        logical :: need_realloc

        need_realloc = .false.
        if (.not. allocated(k1_gas_cache)) then
            need_realloc = .true.
        else if (size(k1_gas_cache,1) /= ngas_species .or. size(k1_gas_cache,2) /= nvar_gas) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(k1_gas_cache)) deallocate(k1_gas_cache, k2_gas_cache, k3_gas_cache, &
                                                    k4_gas_cache, k5_gas_cache, k6_gas_cache)
            if (allocated(y_gas_temp_cache)) deallocate(y_gas_temp_cache)
            if (allocated(error_gas_cache)) deallocate(error_gas_cache)
            allocate(k1_gas_cache(ngas_species, nvar_gas), k2_gas_cache(ngas_species, nvar_gas), &
                     k3_gas_cache(ngas_species, nvar_gas), k4_gas_cache(ngas_species, nvar_gas), &
                     k5_gas_cache(ngas_species, nvar_gas), k6_gas_cache(ngas_species, nvar_gas))
            allocate(y_gas_temp_cache(ngas_species, nvar_gas))
            allocate(error_gas_cache(ngas_species, nvar_gas))
        end if

        need_realloc = .false.
        if (.not. allocated(k1_dust_cache)) then
            need_realloc = .true.
        else if (size(k1_dust_cache) /= ndust_total) then
            need_realloc = .true.
        end if

        if (need_realloc) then
            if (allocated(k1_dust_cache)) deallocate(k1_dust_cache, k2_dust_cache, k3_dust_cache, &
                                                     k4_dust_cache, k5_dust_cache, k6_dust_cache)
            if (allocated(y_dust_temp_cache)) deallocate(y_dust_temp_cache)
            if (allocated(error_dust_cache)) deallocate(error_dust_cache)
            allocate(k1_dust_cache(ndust_total), k2_dust_cache(ndust_total), &
                     k3_dust_cache(ndust_total), k4_dust_cache(ndust_total), &
                     k5_dust_cache(ndust_total), k6_dust_cache(ndust_total))
            allocate(y_dust_temp_cache(ndust_total))
            allocate(error_dust_cache(ndust_total))
        end if
    end subroutine ensure_rk54_cache

    subroutine rk54_step(dust_info,y_gas,y_dust,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break,debug_flag)
        ! Perform a single step of the Runge-Kutta-Fehlberg 5(4) method with adaptive time stepping.
        ! dust_info  --> DustChemistryInfo type with physical parameters.
        ! y_gas      --> 2D array with the gas phase abundances at the start of the step [g cm-3].
        ! y_dust     --> 1D array with the dust phase abundances at the start of the step [g cm-3].
        ! h          --> Time step size [s].
        ! rhs        --> Procedure pointer to the RHS function that computes derivatives.
        ! y_gas_new  <-- 2D array with the updated gas phase abundances [g cm-3].
        ! y_dust_new <-- 1D array with the updated dust phase abundances [g cm-3].
        ! h_new      <-- Recommended time step size for the next step [s].
        ! firstcall  --> Logical flag indicating whether this is the first call in this ODE integration.
        ! accepted   <-- Logical flag indicating whether the step satisfied error limits (10% rule).
        ! break      <-- Logical flag indicating whether integration can be stopped early.
        ! debug_flag --> Optional logical flag to enable verbose debugging.
        implicit none
        ! ---- Input/Output variables ----
        type(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        real(dp), intent(inout) :: h
        procedure(rhs_interface) :: rhs
        real(dp), intent(out) :: y_gas_new(:,:), y_dust_new(:)
        real(dp), intent(out) :: h_new
        logical, intent(out) :: accepted, break
        logical, intent(in) :: firstcall
        logical, intent(in), optional :: debug_flag

        ! ---- Local variables ----
        real(dp), pointer :: k1_gas(:,:), k2_gas(:,:), k3_gas(:,:), k4_gas(:,:), k5_gas(:,:), k6_gas(:,:)
        real(dp), pointer :: k1_dust(:), k2_dust(:), k3_dust(:), k4_dust(:), k5_dust(:), k6_dust(:)
        real(dp), pointer :: y_gas_temp(:,:), y_dust_temp(:)
        real(dp), pointer :: error_gas(:,:), error_dust(:)
        real(dp) :: kmax, scale, max_error
        real(dp) :: h_local

        ! Cash-Karp RK5(4) Coefficients
        real(dp), parameter :: c2 = 0.2_dp
        real(dp), parameter :: c3 = 0.3_dp
        real(dp), parameter :: c4 = 0.6_dp
        real(dp), parameter :: c5 = 1.0_dp
        real(dp), parameter :: c6 = 0.875_dp

        real(dp), parameter :: a21 = 0.2_dp

        real(dp), parameter :: a31 = 3.0_dp / 40.0_dp
        real(dp), parameter :: a32 = 9.0_dp / 40.0_dp

        real(dp), parameter :: a41 = 0.3_dp
        real(dp), parameter :: a42 = -0.9_dp
        real(dp), parameter :: a43 = 1.2_dp

        real(dp), parameter :: a51 = -11.0_dp / 54.0_dp
        real(dp), parameter :: a52 = 2.5_dp
        real(dp), parameter :: a53 = -70.0_dp / 27.0_dp
        real(dp), parameter :: a54 = 35.0_dp / 27.0_dp

        real(dp), parameter :: a61 = 1631.0_dp / 55296.0_dp
        real(dp), parameter :: a62 = 175.0_dp / 512.0_dp
        real(dp), parameter :: a63 = 575.0_dp / 13824.0_dp
        real(dp), parameter :: a64 = 44275.0_dp / 110592.0_dp
        real(dp), parameter :: a65 = 253.0_dp / 4096.0_dp

        ! weights for 5th order solution
        real(dp), parameter :: b1 = 37.0_dp / 378.0_dp
        real(dp), parameter :: b2 = 0.0_dp
        real(dp), parameter :: b3 = 250.0_dp / 621.0_dp
        real(dp), parameter :: b4 = 125.0_dp / 594.0_dp
        real(dp), parameter :: b5 = 0.0_dp
        real(dp), parameter :: b6 = 512.0_dp / 1771.0_dp

        ! weights for error estimation (b_i - b*_i)
        real(dp), parameter :: e1 = 37.0_dp / 378.0_dp - 2825.0_dp / 27648.0_dp
        real(dp), parameter :: e2 = 0.0_dp
        real(dp), parameter :: e3 = 250.0_dp / 621.0_dp - 18575.0_dp / 48384.0_dp
        real(dp), parameter :: e4 = 125.0_dp / 594.0_dp - 13525.0_dp / 55296.0_dp
        real(dp), parameter :: e5 = -277.0_dp / 14336.0_dp
        real(dp), parameter :: e6 = 512.0_dp / 1771.0_dp - 0.25_dp

        call ensure_rk54_cache(size(y_gas,1), size(y_gas,2), size(y_dust))
        k1_gas => k1_gas_cache; k2_gas => k2_gas_cache; k3_gas => k3_gas_cache
        k4_gas => k4_gas_cache; k5_gas => k5_gas_cache; k6_gas => k6_gas_cache
        k1_dust => k1_dust_cache; k2_dust => k2_dust_cache; k3_dust => k3_dust_cache
        k4_dust => k4_dust_cache; k5_dust => k5_dust_cache; k6_dust => k6_dust_cache
        y_gas_temp => y_gas_temp_cache; y_dust_temp => y_dust_temp_cache
        error_gas => error_gas_cache; error_dust => error_dust_cache

        break = .false.

        ! Stage 1
        if (firstcall) then
            if (present(debug_flag)) then
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust,kmax,debug_flag=debug_flag)
            else
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust,kmax)
            end if
            if (kmax == 0.0_dp) then
                break = .true.
                accepted = .true.
                h_new = h
                return
            end if
        else
            if (present(debug_flag)) then
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust,debug_flag=debug_flag)
            else
                call rhs(dust_info,y_gas,y_dust,k1_gas,k1_dust)
            end if
        end if

        if (firstcall) then
            h_local = min(1d0 / kmax, h)
        else
            h_local = h
        end if

        ! Stage 2
        y_gas_temp = y_gas + h_local * a21 * k1_gas
        y_dust_temp = y_dust + h_local * a21 * k1_dust
        if (present(debug_flag)) then
            call rhs(dust_info,y_gas_temp,y_dust_temp,k2_gas,k2_dust,debug_flag=debug_flag)
        else
            call rhs(dust_info,y_gas_temp,y_dust_temp,k2_gas,k2_dust)
        end if

        ! Stage 3
        y_gas_temp = y_gas + h_local * (a31 * k1_gas + a32 * k2_gas)
        y_dust_temp = y_dust + h_local * (a31 * k1_dust + a32 * k2_dust)
        if (present(debug_flag)) then
            call rhs(dust_info,y_gas_temp,y_dust_temp,k3_gas,k3_dust,debug_flag=debug_flag)
        else
            call rhs(dust_info,y_gas_temp,y_dust_temp,k3_gas,k3_dust)
        end if

        ! Stage 4
        y_gas_temp = y_gas + h_local * (a41 * k1_gas + a42 * k2_gas + a43 * k3_gas)
        y_dust_temp = y_dust + h_local * (a41 * k1_dust + a42 * k2_dust + a43 * k3_dust)
        if (present(debug_flag)) then
            call rhs(dust_info,y_gas_temp,y_dust_temp,k4_gas,k4_dust,debug_flag=debug_flag)
        else
            call rhs(dust_info,y_gas_temp,y_dust_temp,k4_gas,k4_dust)
        end if

        ! Stage 5
        y_gas_temp = y_gas + h_local * (a51 * k1_gas + a52 * k2_gas + a53 * k3_gas + a54 * k4_gas)
        y_dust_temp = y_dust + h_local * (a51 * k1_dust + a52 * k2_dust + a53 * k3_dust + a54 * k4_dust)
        if (present(debug_flag)) then
            call rhs(dust_info,y_gas_temp,y_dust_temp,k5_gas,k5_dust,debug_flag=debug_flag)
        else
            call rhs(dust_info,y_gas_temp,y_dust_temp,k5_gas,k5_dust)
        end if

        ! Stage 6
        y_gas_temp = y_gas + h_local * (a61 * k1_gas + a62 * k2_gas + a63 * k3_gas + a64 * k4_gas + a65 * k5_gas)
        y_dust_temp = y_dust + h_local * (a61 * k1_dust + a62 * k2_dust + a63 * k3_dust + a64 * k4_dust + a65 * k5_dust)
        if (present(debug_flag)) then
            call rhs(dust_info,y_gas_temp,y_dust_temp,k6_gas,k6_dust,debug_flag=debug_flag)
        else
            call rhs(dust_info,y_gas_temp,y_dust_temp,k6_gas,k6_dust)
        end if

        ! Compute 5th-order solutions
        y_gas_new = y_gas + h_local * (b1 * k1_gas + b3 * k3_gas + b4 * k4_gas + b6 * k6_gas)
        y_dust_new = y_dust + h_local * (b1 * k1_dust + b3 * k3_dust + b4 * k4_dust + b6 * k6_dust)

        ! Error estimate
        error_gas = abs(h_local * (e1 * k1_gas + e3 * k3_gas + e4 * k4_gas + e5 * k5_gas + e6 * k6_gas)) / &
                    max(abs(y_gas), 1.0d-40)
        error_dust = abs(h_local * (e1 * k1_dust + e3 * k3_dust + e4 * k4_dust + e5 * k5_dust + e6 * k6_dust)) / &
                     max(abs(y_dust), 1.0d-40)

        max_error = maxval(error_gas)
        max_error = max(max_error, maxval(error_dust))

        accepted = (max_error <= errmax)

        scale = 0.9d0 * (errmax / max(max_error, 1.0d-10))**0.2d0
        h_new = h_local * min(2.0_dp, max(0.1_dp, scale))
        h = h_local

    end subroutine rk54_step
end module rk54_mod

module ode_driver_mod
    use amr_parameters, only: dp
    use dustbin_types, only: DustChemistryInfo
    use dust_commons
    use ode_interface_mod
    use dust_rhs_mod, only: print_last_process_kmax, reset_timestep_reduction_counters
    use dust_rhs_mod, only: register_timestep_reduction_cause, print_timestep_reduction_counters
    use dust_rhs_mod, only: last_dydt_dust_per_proc, last_dydt_pah_per_proc

    implicit none
    private
    public :: integrate_dust_ode

    contains

    subroutine integrate_dust_ode(dust_info,dt,y_gas,y_dust,rhs,step_fn,&
                                    &y_gas_final,y_dust_final,h_init,h_min,h_max,debug_flag,step_ok)
        ! Integrate the dust ODE system over a time step dt using the provided step function.
        ! dust_info    --> DustChemistryInfo type with all the necessary information to compute the accretion rate.
        ! dt          --> Time step size [s]
        ! y_gas       --> 2D array with the gas phase abundances at the beginning of the step [g cm-3]
        ! y_dust      --> 1D array with the dust phase abundances at the beginning of the step [g cm-3]
        ! rhs         --> Procedure pointer to the RHS function that computes derivatives.
        ! step_fn     --> ODE solver step function that will be called to perform the integration
        ! y_gas_final --> 2D array with the gas phase abundances at the end of the step [g cm-3]
        ! y_dust_final --> 1D array with the dust phase abundances at the end of the step [g cm-3]
        ! h_init      --> Initial guess for the time step size to be used by the ODE solver [s]
        ! h_min       --> Minimum allowed time step size for the ODE solver [s]
        ! h_max       --> Maximum allowed time step size for the ODE solver [s]
        ! debug_flag  --> Optional logical flag to enable verbose debugging.
        ! step_ok     <-- Optional logical flag indicating whether the single integration step was accepted.

        implicit none
        ! ---- Input/Output variables ----
        type(DustChemistryInfo), intent(in) :: dust_info
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: y_gas(:,:), y_dust(:)
        procedure(rhs_interface) :: rhs
        procedure(solver_step_interface) :: step_fn
        real(dp), intent(out) :: y_gas_final(:,:), y_dust_final(:)
        real(dp), intent(in) :: h_init, h_min, h_max
        logical, intent(in), optional :: debug_flag
        logical, intent(out), optional :: step_ok

        ! ---- Local variables ----
        integer :: icount
        integer :: naccepted, nrejected
        integer :: nreduced
        real(dp) :: h, h_new, tau, h_candidate
        real(dp) :: y_gas_temp(size(y_gas,1),size(y_gas,2)), y_dust_temp(size(y_dust))
        real(dp) :: y_gas_new(size(y_gas,1),size(y_gas,2)), y_dust_new(size(y_dust))
        logical :: accepted,break,firstcall
        logical :: debug_enabled

        debug_enabled = dust_log
        if (present(debug_flag)) debug_enabled = debug_flag

        ! 1. Initialize the time step size for the ODE solver
        h = min(min(max(h_init, h_min), h_max),dt)
        tau = 0.0_dp
        icount = 0
        naccepted = 0
        nrejected = 0
        nreduced = 0
        firstcall = .true.

        call reset_timestep_reduction_counters

        y_gas_temp(:,:) = y_gas(:,:)
        y_dust_temp(:) = y_dust(:)

        if (present(step_ok)) then
            h = dt
            if (present(debug_flag)) then
                call step_fn(dust_info,y_gas_temp,y_dust_temp,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break,debug_flag=debug_flag)
            else
                call step_fn(dust_info,y_gas_temp,y_dust_temp,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break)
            end if
            if (break) then
                y_gas_final(:,:) = y_gas(:,:)
                y_dust_final(:) = y_dust(:)
                step_ok = .true.
                return
            end if
            step_ok = accepted
            if (accepted) then
                y_gas_final(:,:) = y_gas_new(:,:)
                y_dust_final(:) = y_dust_new(:)
            else
                y_gas_final(:,:) = y_gas(:,:)
                y_dust_final(:) = y_dust(:)
            end if
            return
        end if

        ! 2. Integrate the ODE system until we have covered the full time step dt
        do while (tau < dt)
            h = min(h, dt - tau)  ! Adjust the time step size to not overshoot the time step

            ! 3. Call the ODE solver step function to perform a single integration step
            if (present(debug_flag)) then
                call step_fn(dust_info,y_gas_temp,y_dust_temp,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break,debug_flag=debug_flag)
            else
                call step_fn(dust_info,y_gas_temp,y_dust_temp,h,rhs,y_gas_new,y_dust_new,h_new,firstcall,accepted,break)
            end if
            firstcall = .false.
            if (break) then
                ! If break is true, it means there are no active dust processes, so we can skip
                ! the rest of the integration and just return the initial state.
                y_gas_final(:,:) = y_gas(:,:)
                y_dust_final(:) = y_dust(:)
                return
            end if

            ! 4. If the step was accepted, update tau and the solution
            if (accepted) then
                tau  = tau + h
                y_gas_temp(:,:) = y_gas_new(:,:)
                y_dust_temp(:) = y_dust_new(:)
                naccepted = naccepted + 1
                if (dust_log) then
                    ode_naccepted = ode_naccepted + 1_8
                    ! Accumulate per-process dM contribution using k1-stage rates * h
                    if (ndust_processes > 0 .and. allocated(dM_ode_dust) .and. &
                        allocated(last_dydt_dust_per_proc)) then
                        dM_ode_dust(:, 1:ndust_processes) = dM_ode_dust(:, 1:ndust_processes) + &
                            last_dydt_dust_per_proc(:, 1:ndust_processes) * h
                    end if
                    if (npah_processes > 0 .and. allocated(dM_ode_pah) .and. &
                        allocated(last_dydt_pah_per_proc)) then
                        dM_ode_pah(:, 1:npah_processes) = dM_ode_pah(:, 1:npah_processes) + &
                            last_dydt_pah_per_proc(:, 1:npah_processes) * h
                    end if
                end if
            else
                nrejected = nrejected + 1
                if (dust_log) ode_nrejected = ode_nrejected + 1_8
            end if

            ! 5. Update the time step size for the next iteration
            h_candidate = min(max(h_new, h_min), h_max)
            if (h_candidate < h) then
                nreduced = nreduced + 1
                if (dust_log) then
                    ode_nreduced = ode_nreduced + 1_8
                    call register_timestep_reduction_cause
                end if
            end if
            h = h_candidate
            icount = icount + 1
            if (icount > countmax) then
                print *, "Warning: Maximum number of ODE solver iterations reached. Integration may not have converged."
                print *, '  Requested total dt [s]        = ', dt
                print *, '  Integrated tau [s]            = ', tau
                print *, '  Last attempted timestep h [s] = ', h
                print *, '  Last proposed h_new [s]       = ', h_new
                print *, '  Accepted steps                = ', naccepted
                print *, '  Rejected steps                = ', nrejected
                print *, '  Number of timestep reductions = ', nreduced
                call print_last_process_kmax
                call print_timestep_reduction_counters
                call clean_stop
            end if
        end do

        ! 6. Set the final solution
        y_gas_final(:,:) = y_gas_temp(:,:)
        y_dust_final(:) = y_dust_temp(:)

        ! 7. Update per-cell substep statistics (only when logging is active)
        if (dust_log) then
            ndust_cells       = ndust_cells + 1_8
            ode_substeps_sum  = ode_substeps_sum + int(naccepted, kind=8)
            ode_substeps_min  = min(ode_substeps_min, int(naccepted, kind=8))
            ode_substeps_max  = max(ode_substeps_max, int(naccepted, kind=8))
        end if

        if (debug_enabled) then
            if (any(y_gas_final < 0.0_dp) .or. any(y_dust_final < 0.0_dp)) then
                print *, 'DEBUG integrate_dust_ode: negative final density detected.'
                call clean_stop
            end if
        end if
    end subroutine integrate_dust_ode

end module ode_driver_mod