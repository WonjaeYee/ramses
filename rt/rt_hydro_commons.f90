module rt_hydro_commons
  use amr_parameters
  use rt_parameters
  real(dp),allocatable,dimension(:,:)::rtuold,rtunew ! State vector and its update
  real(dp),allocatable,dimension(:)::rt_tau_ir       ! IR cell optical depth, for rt_ir_apflux
end module rt_hydro_commons
