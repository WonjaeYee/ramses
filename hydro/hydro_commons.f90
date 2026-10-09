module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew ! State vector and its update
  real(dp),allocatable,dimension(:,:)::fluxes ! Mass flux on the faces of the cells
  real(dp),allocatable,dimension(:)::pstarold,pstarnew ! Stellar momentum and its update
  real(dp),allocatable,dimension(:)::divu,enew ! Non conservative variables
  real(dp),allocatable,dimension(:)::rho_eq,p_eq ! Strict hydrostatic equilibrium
  real(dp)::mass_tot=0,mass_tot_0=0
  real(dp)::ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  integer::nbins
  ! Passive-scalar constraint tree (init_passive_groups, hydro/godunov_utils.f90):
  ! the ps_ntop "top" slots sum to rho; parent ps_par(j) is the sum of the slots
  ! ps_child(ps_cstart(j):ps_cstart(j+1)-1); ps_gas are the gas-phase slots (elements,
  ! CO, ions, H2) that carry the counter-flux of the TVA drift.
  logical::ps_active=.false.
  integer::ps_ntop=0,ps_npar=0,ps_ngas=0
  integer,allocatable,dimension(:)::ps_top,ps_par,ps_cstart,ps_child,ps_gas,ps_gastop
  integer::ps_ngastop=0
  ! consistency diagnostics accumulated in cooling_fine since the last report
  integer(kind=8)::ps_nbad_top=0,ps_nbad_child=0,ps_ncheck=0,ps_nrenorm=0
  real(dp)::ps_maxdev_top=0,ps_maxdev_child=0
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0d+30
  real(dp)::zero = 0
  real(dp)::one = 1
  real(dp)::two = 2
  real(dp)::three = 3
  real(dp)::four = 4
  real(dp)::two3rd = 2/3d0
  real(dp)::half = 1/2d0
  real(dp)::third = 1/3d0
  real(dp)::forth = 1/4d0
  real(dp)::sixth = 1/6d0
end module const
