!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_commons
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  logical,save:: first_call = .true.       ! True if this is the first call to condinit
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2=neul): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2=neul): P.
  ! If nvar >= ndim+3=nhydro+1, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0 || NVAR>NHYDRO+NENER
  integer::ivar
#endif
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  select case (condinit_kind)

  case('region')
    ! Call built-in initial condition generator
     call region_condinit(x, q, dx, nn)

  case('ana_disk_potential')
     call ana_disk_potential_condinit(x, q, dx, nn)

  case('dustydiffuse')
     call dustydiffuse_condinit(x, q, dx, nn)

  case('dustywave','dustywave_v')
     call dustywave_condinit(x, q, dx, nn)

  case('dustysettle')
     call dustysettle_condinit(x, q, dx, nn)

  case('dustyspress')
     call dustyspress_condinit(x, q, dx, nn)

  case('dustygauss')
     call dustygauss_condinit(x, q, dx, nn)

  case('dustyirtrap')
     call dustyirtrap_condinit(x, q, dx, nn)

  case('dustytrapclump')
     call dustytrapclump_condinit(x, q, dx, nn)

  case('dustyslab')
     call dustyslab_condinit(x, q, dx, nn)

  case('dustyslabthick')
     call dustyslab_condinit(x, q, dx, nn)

  case('dustyslabsil')
     call dustyslab_condinit(x, q, dx, nn)

  case('radsweep')
     call radsweep_condinit(x, q, dx, nn)

  case('dustylev')
     call dustylev_condinit(x, q, dx, nn)

  case('dustylevtrap')
     call dustylev_condinit(x, q, dx, nn)

  case('dustylevtrapseed')
     call dustylev_condinit(x, q, dx, nn)

  case('dustylevatm')
     call dustylevatm_condinit(x, q, dx, nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  case DEFAULT
     if (myid == 1.and. first_call)  write(*,*) "[condinit] Void or invalid condinit_kind, using default IC"
     call region_condinit(x, q, dx, nn)

  end select

  first_call = .false.

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,neul)=0.0d0
  u(1:nn,neul)=u(1:nn,neul)+0.5d0*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,neul)=u(1:nn,neul)+0.5d0*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,neul)=u(1:nn,neul)+0.5d0*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,neul)=u(1:nn,neul)+q(1:nn,neul)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,nhydro+ivar)=q(1:nn,nhydro+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,neul)=u(1:nn,neul)+u(1:nn,nhydro+ivar)
  enddo
#endif
#if NVAR>NHYDRO+NENER
  ! passive scalars
  do ivar=nhydro+1+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine ana_disk_potential_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates an analytical disk potential initial conditions for RAMSES.
  !================================================================
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  integer::i

  real(dp)::height0=150 ! disk height [c.u.]
  real(dp)::dens0=0.66d0 ! central density [c.u.]
  real(dp)::temp0=8000  ! initial temperature [T]

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  do i=1,nn
    ! density
    ! exponential profile along z
    q(i,1) = dens0 * exp(-( x(i,ndim) - 0.5d0 * boxlen)**2 / (2.*height0**2))

    ! pressure via constant temperature
    ! ideal gas: P = n kB T
    q(i,ndim+2) = q(i,1) * (kB*temp0/(mu_gas*mH))/scale_v**2

  end do

end subroutine ana_disk_potential_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustydiffuse_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates the initial conditions for the DustyDiffuse
  ! test case.
  !================================================================
  integer::i,ivar
  real(dp)::dist,eps_val,wsum

  ! Call built-in initial condition generator to set background gas properties
  call region_condinit(x,q,dx,nn)

  do i=1,nn
#if NDUST>0
     ! The bins share the profile in proportion to their var_region (one bin: all of it), so
     ! a multigrain run of identical bins must reproduce the single-bin one (L19 sect. 4.6.1)
     wsum = sum(q(i,idust:idust+ndust-1))
#endif
     ! Distance from the center of the box (0.5 * boxlen)
     dist = x(i,1) - 0.5d0*boxlen
     if (dist > 0.5d0*boxlen) dist = dist - boxlen
     if (dist < -0.5d0*boxlen) dist = dist + boxlen

     ! epsilon = 0.1 * (1 - (dist / (0.5 * boxlen))**2)
     eps_val = 0.1d0 * (1.0d0 - (dist / (0.2d0*boxlen))**2)
     eps_val = max(eps_val, 0.0d0)

     ! Guarded: these dusty ICs are only meaningful with CALIMA, and idust/ndust
     ! do not exist otherwise, so a CALIMA=0 build could not compile this file.
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_val * q(i,idust+ivar-1) / wsum
     end do
#endif

     ! Set gas pressure to an isothermal equation of state:
     ! Pg = c_s_iso**2 * (1.0 - epsilon) * rho_total (with c_s_iso = 1.0)
     q(i,ndim+2) = 1.0d0**2 * (1.0d0 - eps_val) * q(i,1)
  end do

end subroutine dustydiffuse_condinit
!================================================================
!================================================================
subroutine dustywave_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x
  !================================================================
  ! DustyWave (Lebreuilly et al. 2019 sect. 4.4; Laibe & Price 2011): one isothermal sound
  ! mode k = 2 pi/boxlen on the uniform mixture of the region (d_region = rho0,
  ! p_region = P0, var_region = eps of each bin):
  !   rho = rho0 (1 + amp sin kx), v = amp sin kx, P = P0 (1 + amp sin kx), eps unchanged,
  ! with amp = u_region(1) (L19: 1e-4) and P0 = (1 - eps0) rho0 c_s^2.  The values are exact
  ! cell averages, momentum included, so the initial projection adds no error of its own.
  ! condinit_kind='dustywave_v' perturbs the barycentric velocity alone (L19 sect. 4.6.2).
  !================================================================
  integer::i
  real(dp)::k,amp,s1,s2,rho0,p0

  call region_condinit(x,q,dx,nn)
  k = 2.0d0*pi/boxlen
  amp = u_region(1)
  do i=1,nn
     ! cell averages of sin(kx) and sin^2(kx)
     s1 = sin(k*x(i,1)) * sin(0.5d0*k*dx)/(0.5d0*k*dx)
     s2 = 0.5d0 - 0.5d0*cos(2.0d0*k*x(i,1)) * sin(k*dx)/(k*dx)
     rho0 = q(i,1); p0 = q(i,ndim+2)
     if (trim(condinit_kind) == 'dustywave_v') then
        q(i,2) = amp*s1
     else
        q(i,1) = rho0*(1.0d0 + amp*s1)
        q(i,2) = rho0*amp*(s1 + amp*s2) / q(i,1)    ! <rho v> / <rho>
        q(i,ndim+2) = p0*(1.0d0 + amp*s1)
     end if
  end do
end subroutine dustywave_condinit
!================================================================
!================================================================
subroutine dustysettle_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters, only: gravity_params
  use constants
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x
  !================================================================
  ! Dust settling (Lebreuilly et al. 2019 sect. 4.5; Price & Laibe 2015), along the disk's
  ! vertical, in 1D: gravity_type=2 with softening gravity_params(2) = r gives exactly
  ! f = -GM y/(y^2 + r^2)^(3/2), the point mass at y = gravity_params(3).  The region sets the
  ! mid-plane rho0 (d_region), P0 = (1-eps) rho0 c_s^2 (p_region) and eps (var_region); the
  ! isothermal mixture is put in exact hydrostatic equilibrium in that potential,
  !   rho = rho0 exp(-(Phi - Phi_0)/((1-eps) c_s^2)),  Phi = -GM/sqrt(y^2 + r^2).
  !================================================================
  integer::i
  real(dp)::gm,r,y,cs2_eff,eps_tot

  call region_condinit(x,q,dx,nn)
  gm = gravity_params(1); r = gravity_params(2)
  do i=1,nn
     y = x(i,1) - gravity_params(3)
     eps_tot = 0.0d0
#if NDUST>0
     eps_tot = sum(q(i,idust:idust+ndust-1))
#endif
     cs2_eff = q(i,ndim+2)/q(i,1)                    ! P0/rho0 = (1-eps) c_s^2
     q(i,1) = q(i,1)*exp(gm*(1.0d0/sqrt(y**2 + r**2) - 1.0d0/r)/cs2_eff)
     q(i,ndim+2) = cs2_eff*q(i,1)
  end do
end subroutine dustysettle_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustyspress_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates the initial conditions for the DustySpress
  ! test case.
  !================================================================
  integer::i,ivar
  real(dp)::dist,eps_val

  ! Call built-in initial condition generator to set background gas properties
  call region_condinit(x,q,dx,nn)

  do i=1,nn
     ! Distance from x0 = 2.0 pc.
     ! Note: x is in code units, and since boxlen is 10.0 and units_length is 1 pc,
     ! x is directly in parsecs.
     dist = x(i,1) - 2.0d0
     if (dist > 0.5d0*boxlen) dist = dist - boxlen
     if (dist < -0.5d0*boxlen) dist = dist + boxlen

     ! eps(x) = 1e-5 + 0.01 * exp(-(x-x0)^2/sigma^2) with sigma = 0.2
     eps_val = 1d-8 + 1d-4 * exp(-(dist**2)/(2d0 * 0.4d0**2))
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_val
     end do
#endif
#if NPAH>0
     ! the PAH bins drift as TVA species too (dustyshell_pah)
     do ivar=1,npah
        q(i,ipah+ivar-1) = eps_val
     end do
#endif
  end do

end subroutine dustyspress_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustyslab_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Initial conditions for the DustySlab test: uniform static gas with a
  ! dust slab in the middle of the box. UV is injected at the left edge
  ! (a thin rt source region, see the namelist), streams freely through the
  ! dust-poor gas, is absorbed inside the slab and reprocessed into the IR.
  !
  ! The gas density and pressure come from the region parameters, so the same
  ! routine serves both variants: the optically-thin-in-IR case
  ! (d_region -> 1e-22 g/cm3) and the IR-trapping case (d_region ~ 1e-17),
  ! which needs ~5 decades more column to make the CELL-crossing IR optical
  ! depth of order unity.
  !
  ! The gas pressure q(neul) is uniform: it is the pressure both the Godunov
  ! solver and the TVA drift see, so the 1% jump in eps at the slab faces puts
  ! no spurious grad(P_gas) where the radiation-driven drift is being measured.
  !================================================================
  integer::i,ivar
  real(dp),parameter::x_slab_lo = 3.0d0   ! slab inner edge [code length]
  real(dp),parameter::x_slab_hi = 7.0d0   ! slab outer edge [code length]
  real(dp),parameter::eps_slab  = 1.0d-2  ! dust-to-gas ratio inside the slab
  real(dp)::eps_amb                       ! dust-to-gas ratio outside
  real(dp)::eps_val
  logical ::slab_bin_only

  ! The ambient dust must stay optically thin in the UV so the beam reaches the
  ! slab. tau_amb = kappa_UV * eps_amb * rho_gas * 3pc, so the dense variant
  ! (rho_gas ~ 1e-17) needs a far smaller ambient fraction than the fiducial
  ! one (1e-22): at eps_amb=1e-6 the dense ambient would already have tau_UV~32.
  if (trim(condinit_kind) == 'dustyslabthick' .or. &
      trim(condinit_kind) == 'dustyslabsil') then
     eps_amb = 1.0d-10
  else
     eps_amb = 1.0d-6
  end if

  ! 'dustyslabsil' puts the slab in the LAST bin only, so that the bin index --
  ! and therefore which averaged_cross_section_DustBin_NN.txt is used -- selects
  ! the grain. With NDUST=3 that is DustBin_03, 0.005 um silicate, whose UV
  ! opacity is ~38x below the 0.01 um graphite of bin 1 while its far-IR
  ! Rosseland mean is comparable: the least extreme ratio the shipped tables
  ! allow. The other bins are held uniform so they contribute no gradient.
  slab_bin_only = (trim(condinit_kind) == 'dustyslabsil')

  ! Background gas from the region parameters (uniform rho, uniform P)
  call region_condinit(x,q,dx,nn)

  do i=1,nn
     if (x(i,1) >= x_slab_lo .and. x(i,1) <= x_slab_hi) then
        eps_val = eps_slab
     else
        eps_val = eps_amb
     end if
#if NDUST>0
     if (slab_bin_only) then
        do ivar=1,ndust-1
           q(i,idust+ivar-1) = eps_amb          ! uniform: no gradient, no drift
        end do
        q(i,idust+ndust-1) = eps_val
     else
        do ivar=1,ndust
           q(i,idust+ivar-1) = eps_val/dble(ndust)
        end do
     end if
#endif
  end do

end subroutine dustyslab_condinit


subroutine radsweep_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Dust swept up by an attenuated beam (radiation-pressure test 1): region 1 is the
  ! uniform background (gas, ambient dust, metals), and a smooth dust layer
  ! eps_k = eps_k,amb + eps_k,pk sin^2(pi (x - x_lo)/(x_hi - x_lo)) is added on
  ! x_lo..x_hi = x_center(2) -/+ length_x(2)/2, with eps_k,pk = var_region(2, slot of bin k).
  ! Region 2 only carries these parameters: run with nregion=1. The layer is smooth so the
  ! steepening phase can be convergence-tested; the exact solution starts from the snapshot.
  !================================================================
  integer::i,j,iv
  real(dp)::xlo,xhi,shape
  real(dp),parameter::pi=3.14159265358979323846d0

  call region_condinit(x,q,dx,nn)
#if NDUST>0
  xlo = x_center(2) - 0.5d0*length_x(2)
  xhi = x_center(2) + 0.5d0*length_x(2)
  do i=1,nn
     if (x(i,1) > xlo .and. x(i,1) < xhi) then
        shape = sin(pi*(x(i,1)-xlo)/(xhi-xlo))**2
        do j=1,ndust
           iv = idust+j-1
           q(i,iv) = q(i,iv) + var_region(2,iv-nhydro-nener)*shape
        end do
     end if
  end do
#endif

end subroutine radsweep_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustyirtrap_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Initial conditions for the DustyIRtrap test: uniform static gas,
  ! uniform dust-to-gas ratio, and a LINEAR trapped-IR radiation
  ! pressure ramp. With grad(P_trap) exactly constant the terminal
  ! drift has a closed form, and the one-sided differences at the
  ! stencil edges agree with the central ones, so the measured drift
  ! should be uniform across the box.
  !
  ! q(inener) is the non-thermal PRESSURE here: condinit converts it to
  ! energy via /(gamma_rad-1). inener is the slot rt_init assigns to
  ! iIRtrapVar when rt_isIRtrap=.true., so P_trap in the TVA driver is
  ! exactly the profile set below.
  !
  ! Run with NDUST=1 and NPAH=0: the Rosseland opacity share is then
  ! s_1 = 1 identically, so the predicted drift is independent of the
  ! dust optical tables.
  !================================================================
  integer::i,ivar
  ! Log ramp in dust-to-gas ratio across the box. Spanning ~2 decades makes the
  ! first-step diagnostic sample D(rho_d) over a wide range, which tests the
  ! full functional form -- including the barycentric -1/rho_mix subtraction,
  ! a 0.3% correction at eps_min but ~30% at eps_max.
  real(dp),parameter::eps_min  = 3.0d-3  ! dust-to-gas ratio at x=0
  real(dp),parameter::eps_max  = 3.0d-1  ! dust-to-gas ratio at x=boxlen
  ! Sized so |w_d| stays under the tva_wmax_cs=1 cap everywhere (~0.42 c_s at
  ! eps_min, ~0.003 c_s at eps_max), while P_trap stays ~1e-3 of the gas
  ! pressure so the gas itself is barely accelerated.
  real(dp),parameter::Ptrap_0  = 2.0d-4  ! trapped-IR pressure at box centre [code]
  real(dp),parameter::Ptrap_dl = 0.5d0   ! fractional variation across the box
  real(dp)::eps_val

  ! Background gas from the region parameters
  call region_condinit(x,q,dx,nn)

  do i=1,nn
     eps_val = eps_min * (eps_max/eps_min)**(x(i,1)/boxlen)
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_val/dble(ndust)
     end do
#endif
     ! The gas pressure q(neul) (the TVA solver's Pg, as the Godunov solver's) stays
     ! UNIFORM despite the dust ramp: that kills the ordinary grad(P_gas) term,
     ! leaving the trapped-IR gradient as the only driver of the drift and making
     ! the first-step check a clean test of it alone.
#if NENER>0
     ! P_trap(x) = Ptrap_0 * (1 - Ptrap_dl*(x/boxlen - 1/2))
     ! => grad(P_trap) = -Ptrap_0*Ptrap_dl/boxlen, constant and negative,
     !    so the dust drifts towards +x.
     q(i,inener) = Ptrap_0 * (1.0d0 - Ptrap_dl*(x(i,1)/boxlen - 0.5d0))
#endif
  end do

end subroutine dustyirtrap_condinit


subroutine dustytrapclump_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Initial conditions for the DustyTrapClump test: uniform static gas, a
  ! Gaussian dust clump on a uniform dust background, and a trapped-IR
  ! pressure held fixed (TVA_TEST_IRTRAP) whose gradient G(x) = -G_0 g(x),
  ! g = [tanh((x-x_1)/delta) - tanh((x-x_2)/delta)]/2, is a plateau in the
  ! middle of the box and vanishes near the edges, so that no dust has to
  ! enter through a boundary. With one bin the bin carries all the Rosseland
  ! opacity (zeta = 1), so the trapped-IR dust flux is
  ! F = -t_0 G (1 - rho_d/rho), t_s = t_0/(1-eps): the grains stream through
  ! the gas at w = F/rho_d, but where G is uniform the clump obeys
  ! d(rho_d)/dt + lambda d(rho_d)/dx = 0, lambda = t_0 G/rho ~ -eps w. It
  ! must keep its shape and move only at lambda, slowly and against the
  ! drift. A dust flux upwinded on the drift diffuses it with D ~ |w| dx/2.
  ! Where G rises (falls) along the drift the dust thins (piles up).
  !
  ! Run with NDUST=1 and NPAH=0, and an isothermal gas (rt_Tconst), so
  ! that t_0 is uniform.
  !================================================================
  integer::i,ivar
  real(dp),parameter::eps_bg   = 4.0d-3  ! background dust-to-gas ratio
  real(dp),parameter::eps_pk   = 8.0d-3  ! dust-to-gas ratio at the clump centre
  real(dp),parameter::sig_frac = 0.008d0 ! clump sigma / boxlen (~2 cells at levelmin 8)
  ! G_0 ten times the dustyirtrap ramp: |w| ~ 0.3 c_s at eps_bg, P_trap ~ 1% of the gas pressure
  real(dp),parameter::Ptrap_0  = 2.0d-3  ! trapped-IR pressure at x = 0 [code]
  real(dp),parameter::G_0      = 1.0d-4  ! plateau |grad(P_trap)| [code]
  real(dp),parameter::x1_frac  = 0.25d0  ! plateau edges / boxlen
  real(dp),parameter::x2_frac  = 0.75d0
  real(dp),parameter::dl_frac  = 0.125d0 ! transition half-width delta / boxlen
  real(dp)::eps_val, x1, x2, dl

  call region_condinit(x,q,dx,nn)
  x1 = x1_frac * boxlen
  x2 = x2_frac * boxlen
  dl = dl_frac * boxlen

  do i=1,nn
     eps_val = eps_bg + (eps_pk - eps_bg) * exp(-0.5d0*((x(i,1)/boxlen - 0.5d0)/sig_frac)**2)
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_val/dble(ndust)
     end do
#endif
     ! Uniform gas pressure despite the clump (Pg = q(neul)): no grad(P_gas) drift
#if NENER>0
     ! P_trap = Ptrap_0 - G_0 int_0^x g: (delta/2) [ln cosh((x-x1)/delta) - ln cosh((x-x2)/delta)]
     q(i,inener) = Ptrap_0 - G_0 * 0.5d0 * dl * (log(cosh((x(i,1)-x1)/dl)) - log(cosh(x1/dl)) &
                                               - log(cosh((x(i,1)-x2)/dl)) + log(cosh(x2/dl)))
#endif
  end do

end subroutine dustytrapclump_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustylev_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Initial conditions for the DustyLev test: the RT15 sec. 3.7 /
  ! Davis et al. 2014 dust levitation experiment, recast so that the
  ! answer is a NUMBER rather than "rises or falls".
  !
  ! An isolated Gaussian dusty layer sits at the centre of the box in a
  ! uniform downward gravity field (gravity_type=1, gravity_params(1)=-g),
  ! illuminated from x=0 by a directed IR beam (a thin rt source region,
  ! see the namelist). Because the density falls to a floor at BOTH
  ! boundaries, the boundary pressure force vanishes and the layer-integrated
  ! momentum budget closes exactly,
  !
  !     d<v>/dt = g * (f_E - 1),   f_E = (F/c)(1-exp(-tau_R)) / (g*Sigma)
  !
  ! independently of how the layer restructures internally. In the optically
  ! thin limit this reduces to the mass-specific Eddington ratio
  ! f_E = kappa_R*eps*F/(c*g), which has no Sigma in it at all.
  !
  ! Deliberately an ISOLATED layer rather than the wall-supported exponential
  ! atmosphere of the Davis patch: that one needs a reflecting lower boundary
  ! whose pressure force enters the budget, so the net acceleration is no
  ! longer a closed form.
  !
  ! The layer amplitude and the isothermal sound speed come from the region
  ! parameters (d_region, p_region), so the same routine serves the optically
  ! thin variants and the trapping variant.
  !================================================================
  integer::i,ivar
  real(dp),parameter::eps_lev    = 1.0d-2  ! dust-to-gas ratio, UNIFORM
  real(dp),parameter::sig_frac   = 0.1d0   ! layer sigma / boxlen
  real(dp),parameter::seed_amp   = 1.0d-2  ! 'dustylevtrapseed': relative odd-even perturbation of eps
  real(dp)::x_c, sig, rho_pk, cs2, rho_loc, eps_loc
  real(dp)::floor_frac                     ! ambient density / peak density

  ! The ambient has the same dust fraction as the layer, so in the optically
  ! thin variants it shares the layer's Eddington ratio and rides along at the
  ! same acceleration -- harmless, and 1e-3 keeps its density well away from
  ! the chemistry's low-density corner. The trapping variant is deliberately
  ! super-Eddington in the thin limit (f_E ~ 41), so there the ambient would be
  ! blasted to tens of sound speeds and, at 1e-3, would carry more momentum
  ! than the layer itself. Drop it two more decades there.
  if (trim(condinit_kind) == 'dustylevtrap' .or. trim(condinit_kind) == 'dustylevtrapseed') then
     floor_frac = 1.0d-5
  else
     floor_frac = 1.0d-3
  end if

  ! Uniform background from the region parameters: q(1,1)=d_region gives the
  ! layer peak and q(1,neul)/q(1,1)=p_region/d_region the isothermal c_iso^2.
  call region_condinit(x,q,dx,nn)

  x_c = 0.5d0 * boxlen
  sig = sig_frac * boxlen

  do i=1,nn
     rho_pk  = q(i,1)
     cs2     = q(i,neul) / q(i,1)
     rho_loc = rho_pk * exp(-0.5d0*((x(i,1)-x_c)/sig)**2)
     rho_loc = max(rho_loc, floor_frac*rho_pk)
     q(i,1)  = rho_loc
     ! Isothermal, so that rt_Tconst keeps the run consistent with the IC and
     ! the only vertical forces are gravity, radiation and grad(P).
     q(i,neul) = cs2 * rho_loc
     ! The dust fraction is uniform: a dust gradient would add a TVA drift on
     ! top of the barycentric force this test is measuring. 'dustylevtrapseed'
     ! (the stability test of the trapped-IR drift) perturbs it by +-seed_amp in
     ! alternate cells within 2 sigma of the centre: the odd-even mode, which an
     ! anti-diffusive dust flux would amplify first.
     eps_loc = eps_lev
     if (trim(condinit_kind) == 'dustylevtrapseed' .and. abs(x(i,1)-x_c) < 2d0*sig) then
        eps_loc = eps_lev * (1d0 + seed_amp * (1 - 2*modulo(int(x(i,1)/dx), 2)))
     end if
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_loc/dble(ndust)
     end do
#endif
     ! q(neul) is the gas pressure the TVA solver reconstructs (the Godunov one):
     ! the isothermal profile above.
  end do

end subroutine dustylev_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustylevatm_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Initial conditions for the DustyLevAtm test: the Krumholz & Thompson
  ! (2013) / Davis et al. (2014) dusty atmosphere as repeated by Rosdahl &
  ! Teyssier 2015 (RT15) sec. 3.7, so that the diagnostics of their fig. 11
  ! -- f_E,V, tau_V and tau_F/tau_V versus time -- can be compared directly.
  !
  ! An exponential layer sits at the BOTTOM of the box,
  !
  !     rho(h) = rho_* exp(-h/h_*),   floored at 1e-10 rho_*
  !
  ! isothermal at T_*, in a uniform downward gravity field, illuminated from
  ! below by a vertical IR flux F_*. Height is x(i,ndim): x in 1D, y in 2D.
  !
  ! Units are chosen so that the RT15 numbers come out as pure ones:
  ! units_length = h_*, units_time = t_* = h_*/c_*, units_density = rho_*.
  ! Then c_s = 1, g_code = c_iso^2/h_* * t_*^2/h_* = 1 (hence
  ! gravity_params = -1, as in patch/rt/davis/davis.nml), d_region = 1,
  ! p_region = 1, and Sigma = 1 in code units.
  !
  ! The gas starts at h = 0 as in RT15: 25-60% of the column sits in the bottom
  ! row of cells, which is where the trapped-photon pile-up forms. That needs
  ! the lower boundary to EMIT rather than merely impose a state, which is what
  ! the local patch/rt_hydro_boundary.f90 (RT15 eq. 83) provides.
  !
  ! In >=2 D the density gets the fluctuations of RT15's own code,
  ! patch/rt/davis/condinit.f90, which seed the radiative Rayleigh-Taylor
  ! instability:
  !     rho -> rho * (1 + 0.25*sin(8*pi*x/boxlen)) * (1 + chi),  chi in [-0.25, 0.25]
  ! i.e. a 25% mode of wavelength boxlen/4 times independent +-25% noise in
  ! every cell. This is NOT RT15's printed eq. 76 (one mode across the box, the
  ! noise only modulating its amplitude), which at 2048^2 left the layer intact
  ! at 25 t_* where RT15's fig. 10 has broken into fingers. In 1D there is no
  ! horizontal direction, so the 1D variants are the unperturbed reference.
  !================================================================
  integer::i
#if NDUST>0
  integer::ivar
  ! Dust-to-gas ratio, UNIFORM. 0.016359 is not arbitrary: CALIMA's Rosseland
  ! mean for the shipped 0.01 um graphite is 129.881 cm2 per gram of DUST at
  ! T = 82 K, so this eps makes kappa per gram of GAS equal 2.1248 cm2/g --
  ! exactly RT15's kappa_R,* from their eq. 79. With it, tau_* = 3.0017 and
  ! f_E,* = 0.50001, i.e. the CALIMA run has RT15's ENTIRE dimensional setup
  ! (F_*, T_*, g, Sigma, rho_*, h_*, t_*), and the only differences are the
  ! temperature dependence of kappa and the dust dynamics.
  real(dp),parameter::eps_lev    = 0.016359d0
#endif
  real(dp),parameter::floor_frac = 1.0d-10  ! RT15: floor at 1e-10 rho_*
  real(dp),parameter::pert_amp   = 0.25d0   ! RT15 davis/condinit.f90
  real(dp)::rho_pk, cs2, rho_loc, h, chi

  ! Uniform background from the region parameters: q(i,1) = d_region = rho_*
  ! and q(i,neul)/q(i,1) = p_region/d_region = c_iso^2.
  call region_condinit(x,q,dx,nn)

  do i=1,nn
     rho_pk = q(i,1)
     cs2    = q(i,neul) / q(i,1)
     ! EXACT cell average of exp(-h/h_*), as patch/rt/davis/condinit.f90 does:
     !   rho_bar = rho_* (exp(-h_lo) - exp(-h_hi)) / dx
     ! Point-sampling the exponential instead loses column density as the cell
     ! gets thick -- 14% at dy = 2 h_*, which shifts tau_* away from RT15's 3.
     ! This form is exact for any dy, so Sigma = rho_* h_* whatever the grid.
     h       = max(x(i,ndim) - 0.5d0*dx, 0.0d0)
     rho_loc = rho_pk * (exp(-h) - exp(-(x(i,ndim)+0.5d0*dx))) / dx
#if NDIM>1
     ! RT15's code: rho*(1 + 0.25 sin(4*3.14*x/0.5/boxlen)) * (1 + 0.5*(0.5 + R)),
     ! R from ranf with the unseeded localseed = -1, which gives R uniform in
     ! [-1, 0] (emulated), so the noise factor is uniform in [0.75, 1.25] and
     ! keeps the column. Two departures: 2*pi*4 instead of 4*3.14/0.5, so the
     ! mode is periodic in the box, and chi from a deterministic position hash
     ! instead of ranf, whose pattern depends on the MPI decomposition (every
     ! rank restarts from the same seed). The hash is reproducible run to run.
     chi = 0.5d0 * (mod(abs(sin(x(i,1)*12.9898d0 + x(i,ndim)*78.233d0))       &
                        * 43758.5453d0, 1.0d0) - 0.5d0)
     rho_loc = rho_loc * (1.0d0 + pert_amp * sin(4.0d0 * twopi * x(i,1) / boxlen)) &
                       * (1.0d0 + chi)
#endif
     rho_loc = max(rho_loc, floor_frac * rho_pk)
     q(i,1)    = rho_loc
     q(i,neul) = cs2 * rho_loc             ! isothermal at T_*
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_lev/dble(ndust)
     end do
     ! q(neul) is the gas pressure the TVA solver reconstructs (the Godunov one):
     ! the isothermal profile above.
#endif
  end do

end subroutine dustylevatm_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine dustygauss_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use constants

  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer::i,ivar
  real(dp)::dist,eps_val

  ! Call built-in initial condition generator to set background gas properties
  call region_condinit(x,q,dx,nn)

  do i=1,nn
     ! Distance from L/2
     dist = x(i,1) - 0.5d0*boxlen
     if (dist > 0.5d0*boxlen) dist = dist - boxlen
     if (dist < -0.5d0*boxlen) dist = dist + boxlen

     ! f(x) = 0.01 + 0.1 * exp(-((x-L/2)/(L/4))^2)
     eps_val = 0.01d0 + 0.1d0 * exp(- (dist / (0.25d0*boxlen))**2)
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_val
     end do
#endif
  end do
end subroutine dustygauss_condinit
