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

  case('dustyspress')
     call dustyspress_condinit(x, q, dx, nn)

  case('dustygauss')
     call dustygauss_condinit(x, q, dx, nn)

  case('dustyirtrap')
     call dustyirtrap_condinit(x, q, dx, nn)

  case('dustyslab')
     call dustyslab_condinit(x, q, dx, nn)

  case('dustyslabthick')
     call dustyslab_condinit(x, q, dx, nn)

  case('dustyslabsil')
     call dustyslab_condinit(x, q, dx, nn)

  case('dustylev')
     call dustylev_condinit(x, q, dx, nn)

  case('dustylevtrap')
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
  real(dp)::dist,eps_val

  ! Call built-in initial condition generator to set background gas properties
  call region_condinit(x,q,dx,nn)

  do i=1,nn
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
        q(i,idust+ivar-1) = eps_val
     end do
#endif

     ! Set gas pressure to an isothermal equation of state:
     ! Pg = c_s_iso**2 * (1.0 - epsilon) * rho_total (with c_s_iso = 1.0)
     q(i,ndim+2) = 1.0d0**2 * (1.0d0 - eps_val) * q(i,1)
  end do

end subroutine dustydiffuse_condinit


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
  ! The gas pressure is pre-divided by (1-eps) so that the pressure the TVA
  ! solver actually sees, Pg = (1-eps_tot)*q(neul), is uniform. Without this
  ! the 1% jump in eps at the slab faces would put a spurious grad(P_gas)
  ! exactly where the radiation-driven drift is being measured.
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
     q(i,neul) = q(i,neul) / (1.0d0 - eps_val)
  end do

end subroutine dustyslab_condinit


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
     ! The TVA solver's gas pressure is Pg = (1-eps_tot)*q(neul) (it removes the
     ! NENER energy, then scales by the gas fraction). Pre-divide by (1-eps) so
     ! Pg comes out UNIFORM despite the dust ramp: that kills the ordinary
     ! grad(P_gas) term, leaving the trapped-IR gradient as the only driver of
     ! the drift and making the first-step check a clean test of it alone.
     q(i,neul) = q(i,neul) / (1.0d0 - eps_val)
#if NENER>0
     ! P_trap(x) = Ptrap_0 * (1 - Ptrap_dl*(x/boxlen - 1/2))
     ! => grad(P_trap) = -Ptrap_0*Ptrap_dl/boxlen, constant and negative,
     !    so the dust drifts towards +x.
     q(i,inener) = Ptrap_0 * (1.0d0 - Ptrap_dl*(x(i,1)/boxlen - 0.5d0))
#endif
  end do

end subroutine dustyirtrap_condinit


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
  real(dp)::x_c, sig, rho_pk, cs2, rho_loc
  real(dp)::floor_frac                     ! ambient density / peak density

  ! The ambient has the same dust fraction as the layer, so in the optically
  ! thin variants it shares the layer's Eddington ratio and rides along at the
  ! same acceleration -- harmless, and 1e-3 keeps its density well away from
  ! the chemistry's low-density corner. The trapping variant is deliberately
  ! super-Eddington in the thin limit (f_E ~ 41), so there the ambient would be
  ! blasted to tens of sound speeds and, at 1e-3, would carry more momentum
  ! than the layer itself. Drop it two more decades there.
  if (trim(condinit_kind) == 'dustylevtrap') then
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
     ! top of the barycentric force this test is measuring.
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_lev/dble(ndust)
     end do
#endif
     ! Same (1-eps) pre-division as the slab tests, so the gas pressure the
     ! TVA solver reconstructs matches the isothermal profile above.
     q(i,neul) = q(i,neul) / (1.0d0 - eps_lev)
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
  ! In >=2 D the density gets RT15 eq. 76 fluctuations,
  !     drho/rho = 0.25*(1 + chi)*sin(2*pi*x/boxlen),  chi in [-0.25, 0.25],
  ! which seed the radiative Rayleigh-Taylor instability. In 1D there is no
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
  real(dp),parameter::pert_amp   = 0.25d0   ! RT15 eq. 76
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
     ! RT15 eq. 76: seed the radiative Rayleigh-Taylor instability with
     !   drho/rho = pert_amp*(1 + chi)*sin(2 pi x / boxlen),  chi in [+-0.25].
     ! chi comes from a deterministic position hash rather than ranf: with the
     ! default unseeded localseed = -1, ranf did not return a uniform [0,1]
     ! deviate here and the coherent amplitude came out at 0.125 instead of
     ! 0.25 (measured). The hash is reproducible run to run, which is what a
     ! regression test wants anyway.
     chi = 0.5d0 * (mod(abs(sin(x(i,1)*12.9898d0 + x(i,ndim)*78.233d0))       &
                        * 43758.5453d0, 1.0d0) - 0.5d0)
     rho_loc = rho_loc * (1.0d0 + pert_amp * (1.0d0 + chi)                  &
                                * sin(twopi * x(i,1) / boxlen))
#endif
     rho_loc = max(rho_loc, floor_frac * rho_pk)
     q(i,1)    = rho_loc
     q(i,neul) = cs2 * rho_loc             ! isothermal at T_*
#if NDUST>0
     do ivar=1,ndust
        q(i,idust+ivar-1) = eps_lev/dble(ndust)
     end do
     ! Same (1-eps) pre-division as the other dusty tests, so the gas pressure
     ! the TVA solver reconstructs is the isothermal profile above.
     q(i,neul) = q(i,neul) / (1.0d0 - eps_lev)
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
