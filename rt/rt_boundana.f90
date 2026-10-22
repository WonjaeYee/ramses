!############################################################
!############################################################
!############################################################
!############################################################
subroutine rt_boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use rt_parameters, ONLY: nrtvar, rt_boundary_var, rt_flux_bound, rt_u_bound, nGroups, rt_c_cgs
  use amr_commons, ONLY: levelmin
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nrtvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): N, U(i,2:ndim+1): Fx,Fy,Fz.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i,ig,iN
  real(dp)::scale_Np,scale_Fp

  do ivar=1,nrtvar
     do i=1,ncell
        u(i,ivar)=rt_boundary_var(ibound,ivar)
     end do
  end do

  ! rt_flux_bound: each group's photon flux through this boundary, a beam along rt_u_bound, as a region
  ! source sets it (rt_init): N = F/c, F_x = u F. The reduced speed of light of levelmin (uniform grids)
  call rt_units(scale_Np, scale_Fp)
  do ig=1,nGroups
     if (rt_flux_bound(ibound,ig) <= 0d0) cycle
     iN = (ig-1)*(ndim+1) + 1
     do i=1,ncell
        u(i,iN) = rt_flux_bound(ibound,ig) / rt_c_cgs(levelmin) / scale_Np
#if NDIM>0
        u(i,iN+1) = rt_u_bound(ibound) * rt_flux_bound(ibound,ig) / scale_Fp
#endif
     end do
  end do

  ! Add here, if you wish, some user-defined boudary conditions
  ! ........

end subroutine rt_boundana
