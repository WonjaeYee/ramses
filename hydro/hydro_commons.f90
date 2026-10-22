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
  ! The uold slots of the gas-phase metals, set once in read_hydro_params: with RTZ the tracked elements beyond H and
  ! He, without RTZ the metallicity scalars. Used by the metallicity and DTM movie maps and CALIMA's box masses.
  integer::n_metal_idx=0
  integer,dimension(1:max(nmetals,1))::metal_idx=0
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
