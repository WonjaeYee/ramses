! molecules_module.f90
module molecules_module
  use amr_parameters, only: dp
  use safe_math, only: safe_exp
  implicit none

  private  ! everything is private by default
  public :: alpha_H2, beta_H2_umist, beta_H2, beta_H2_krome, alpha_CO, beta_CO, alpha_H2_prim, alpha_H2_dust
  public :: comp_SH2, comp_Sd, comp_SCO, lw_transmission, tau_dust_lw_mw, co_dust_to_lw

  ! Dust attenuation in the Lyman-Werner band (1000 A) per H nucleus, for MW dust (R_V = 3.1), on the
  ! total hydrogen column N(H) + 2N(H2) + N(H+): Draine & Bertoldi (1996, Sect. 2.4), an effective
  ! attenuation cross-section (~0.76 of the extinction one, as forward-scattered light is not removed)
  real(dp), parameter :: sigma_dust_lw_mw = 2.0d-21
  ! CO and H2 dust attenuation in units of A_V: CO exp(-3.53 A_V) (Visser et al. 2009; van Dishoeck
  ! et al. 2006), H2 sigma_dust_lw_mw N_H = 3.74 A_V for N_H/A_V = 1.87e21 cm^-2 mag^-1
  real(dp), parameter :: co_dust_to_lw = 3.53d0 / 3.74d0

#include "co_shield_visser09.inc"

CONTAINS

FUNCTION alpha_H2_prim(T, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) result(rate)
  implicit none

  real(dp), intent(in) :: T, xe, H2_cosmic_ray_ionization_rate, G0
  real(dp), intent(in) :: xHI, xHII, nH
  real(dp) :: rate
  real(dp) :: logT, lnTe, logT2, Te
  real(dp) :: k1, k2, k5, k13, k14, k15, k_hm_cr, k_hm_gamma
  real(dp) :: cr_photo_destruction, collisional_destruction

  ! H- channel for H2 formation
  rate = 0.d0

  ! Primordial channel
  logT = log10(T)
  Te = T*8.621738d-5 ! K -> eV
  lnTe = log(Te)

  ! Creation and destruction channels of H- included with updated rates from Glover et al. 2010
  ! H + e- -> H- + gamma
  if (T .lt. 6000.d0) then
     k1 = 10.d0**(-17.845d0 + logT * (0.762d0 + logT * (0.1523d0 - 0.03274d0 * logT)))
  else
     logT2 = logT * logT
     k1 = 10.d0**(-16.42d0 + logT2 * (0.1998d0 + logT2 * (-5.447d-3 + 4.0415d-5 * logT2)))
  end if

  ! H- + H -> H2 + e
  k2 = 4.0d-9*(max(T,300.d0)**(-0.17d0))

  ! H- + H+ -> H + H
  k5 = 2.4d-6/sqrt(T)*(1.d0 + T/20000.d0)

  ! H- + CR --> H + e-
  k_hm_cr = 1.28d-13 * (H2_cosmic_ray_ionization_rate / 1.d-16)

  ! H- + gamma --> H + e-
  k_hm_gamma = 5.9d-9 * G0

  ! H- + e -> H + e + e
  k13 = -1.801849334d1 + lnTe * (2.36085220d0 + lnTe * (-2.82744300d-1 + lnTe * (1.62331664d-2 + lnTe * (-3.36501203d-2 + lnTe * (1.17832978d-2 + lnTe * (-1.65619470d-3 + lnTe * (1.06827520d-4 - 2.63128581d-6 * lnTe)))))))
  k13 = safe_exp(max(k13,-92.d0)) ! Max needed to prevent result from diverging at low temperature

  ! H- + H --> H + H + e-
  ! I think this reaction was broken in glover so I took the results from
  ! https://www.aanda.org/articles/aa/pdf/2016/02/aa27262-15.pdf Table A1
  ! Harley added the fudge factor for continuity
  k14 = 1.357772745525155d0 * 2.5634d-15 * (Te**1.78186d0) ! Note that T must be in eV for this reaction to make sense
  if (T .gt. 1160.d0) then
     k14 = -3.388464953d1 + lnTe * (1.13944933d0 + lnTe * (-1.4210135d-1 + lnTe * (8.4644554d-3 + lnTe * (-1.4328641d-3 + lnTe * (2.0122503d-4 + lnTe * (8.6639632d-5 + lnTe * (-2.5850097d-5 + lnTe * (2.4555012d-6 - 8.0683825d-8 * lnTe))))))))
     k14 = safe_exp(k14)
  end if

  ! H- + H+ --> H2+ + e- --> H + H (via recombinative dissociation)
  k15 = 6.9d-9 * (T**(-0.35d0))
  if (T .gt. 8000.d0) then
     k15 = 9.6d-7 * (T**(-0.9d0))
  end if

  ! Correct steady state denominator terms scaled to cm^3 s^-1
  cr_photo_destruction = (k_hm_cr + k_hm_gamma) / (nH + 1d-40)
  collisional_destruction = k2 * xHI + k5 * xHII + k13 * xe + k14 * xHI + k15 * xHII

  rate = rate + k1 * k2 * xe * xHI / (collisional_destruction + cr_photo_destruction)

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_H2_prim

FUNCTION alpha_H2_dust(T, dust_to_gas_mass_ratio_over_mw) result(rate)
  ! Formation on dust
  use rt_parameters, only: rtz_H2_clumping
  implicit none

  real(dp), intent(in) :: T, dust_to_gas_mass_ratio_over_mw
  real(dp) :: rate
  real(dp) :: T2

  T2 = T / 100.d0
  rate = dust_to_gas_mass_ratio_over_mw * (3.5d-17) * rtz_H2_clumping * sqrt(min(T2,1.d2))

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_H2_dust

FUNCTION alpha_H2(T, dust_to_gas_mass_ratio_over_mw, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) result(rate)
  ! Creation rate of molecular hydrogen
  ! We consider both the primordial channel (via H-) as well
  ! as formation on dust
  implicit none
  real(dp), intent(in) :: T, dust_to_gas_mass_ratio_over_mw, xe
  real(dp), intent(in) :: H2_cosmic_ray_ionization_rate, G0, xHI
  real(dp), intent(in) :: xHII, nH
  real(dp) :: rate

  rate = 0.d0

  ! Formation rate on dust. Consider only HI
  if (dust_to_gas_mass_ratio_over_mw.gt.0d0) then
     rate = rate + alpha_H2_dust(T, dust_to_gas_mass_ratio_over_mw) * xHI * nH
  end if

  ! Primordial H- channel
  rate = rate + alpha_H2_prim(T, xe, H2_cosmic_ray_ionization_rate, G0, xHI, xHII, nH) * xHI * nH

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_H2

FUNCTION beta_H2_umist(T, nH, ne, nH2) result(rate)
  ! H2 destruction from umist
  implicit none

  real(dp), intent(in) :: T, nH, ne, nH2
  real(dp) :: rate
  real(dp) :: T_loc

  rate = 0.d0

  !H2 + H2 --> H2 + H + H 
  T_loc = max(min(T,41000.d0),2803.d0)
  rate = rate + (1.00d-8 * ((T_loc/300d0)**0.0d0) * safe_exp(-84100.d0/T_loc) * nH2)

  !H2 + e- --> H + H + e-
  T_loc = max(min(T,41000.d0),3400.d0)
  rate = rate + (3.22d-9 * ((T_loc/300d0)**0.35d0) * safe_exp(-102000.d0/T_loc) * ne)

  !H2 + H --> H + H + H
  T_loc = max(min(T,41000.d0),1833.d0)
  rate = rate + (4.67d-7 * ((T_loc/300d0)**(-1.d0)) * safe_exp(-55000.d0/T_loc) * nH)

  rate = MAX(rate,1.d-100)

END FUNCTION beta_H2_umist

FUNCTION beta_H2_krome(T, nH, ne, nH2, nHe) result(rate)
  ! From bovino 2016
  ! https://www.aanda.org/articles/aa/pdf/2016/06/aa28158-16.pdf
  implicit none
  real(dp), intent(in) :: T, nH, ne, nH2, nHe
  real(dp) :: rate

  rate = 0.d0

  ! H2 + H --> H + H + H (k18)
  rate = rate + ((6.67d-12 * sqrt(T) * safe_exp(-1.d0 * (1.d0 + (63593.d0/T)))) * nH)

  ! H2 + H2 --> H2 + H + H (k19)
  rate = rate + ((5.996d-30 * (T**4.1881d0) * ((1.d0 + (6.761d-6 * T))**(-5.6881d0)) * safe_exp(-54657.4d0/T)) * nH2)

  ! H2 + e- --> H + H + e- (k22)
  rate = rate + (4.38d-10 * (T**0.35d0) * safe_exp(-102000.d0/T) * ne)

  ! H2 + He --> H + H + He
  rate = rate + ((10.d0**(-27.029d0 + (3.801d0 * log10(T)) - (29487.d0/T))) * nHe)

  rate = MAX(rate,1.d-100)

END FUNCTION beta_H2_krome

FUNCTION beta_H2(T, nH, xHI, xH2, xHe, ne, nHI, nH2, nHeI) result(rate)
  ! Returns the collisional dissociation rates of H2 for four different
  ! reactions [cm3s-1] from Glover & Abel (2008)
  ! http://mnras.oxfordjournals.org/content/388/4/1627.full.pdf 
  implicit none

  real(dp), intent(in):: T, nH, xHI, xH2, xHe
  real(dp), intent(in):: ne, nHI, nH2, nHeI
  real(dp):: rate
  real(dp):: T4, ncrH, ncrH2, ncrHe, invncr
  real(dp):: LTEfac, NLTEfac
  real(dp):: k8, k9, k9L, k10, k10L, k11, k11L
  real(dp):: lk8, lk9, lk10, lk11

  T4 = T / 1d4

  ! Critical number densities.  See eqns 15, 16, and 17
  ncrH =  10.d0**(3.d0 - 0.416d0*log10(T4) - 0.327d0*log10(T4)*log10(T4))
  ncrH2 = 10.d0**(4.845d0 - 1.3d0*log10(T4) + 1.62d0*log10(T4)*log10(T4))
  ncrHe = 10.d0**(5.0792d0*(1.d0 - 1.23d-5*(T - 2000.d0)))

  ! 1/ncr.  see eqn 14
  invncr = (xHI/ncrH) + (xH2/ncrH2) + (xHe/ncrHe)

  ! prefactors for LTE and NLTE collision rates  see eqn 13
  LTEfac = (nH*invncr)/(1.d0 + (nH*invncr))
  NLTEfac = 1.d0/(1.d0 + (nH*invncr))
  if (ne .gt. nHI) then
     LTEfac = 1.d0
     NLTEfac = 0.d0
  end if

  !reaction rates from the appendix:
  !H2 + e- --> H + H + e-
  k8 = 3.73d-9 * (T**0.1121d0) * safe_exp(-99430.d0/T) ! Glover et al. (2010)

  !H2 + H --> H + H + H
  k9 = (6.67d-12)*sqrt(T)*safe_exp(-1.d0*(1.d0 + (63593.d0/T)))
  k9L = (3.52d-9)*safe_exp(-43900.d0/T)

  !H2 + H2 --> H2 + H + H
  k10 = ( (5.996d-30*(T**4.1881d0)) / ((1.d0 + 6.761d-6*T)**5.6881d0)) * safe_exp(-54657.4d0/T)
  k10L = (1.3d-9)*safe_exp(-53300.d0/T)

  !H2 + He --> H + H + He
  k11 = 10.d0**(-27.029d0 + (3.801d0*log10(T)) - (29487.d0/T))
  k11L = 10.d0**(-2.729d0 - (1.75d0*log10(T)) - (23474.d0/T))

  k8   = max(k8, 1d-40)
  k9   = max(k9, 1d-40)
  k10  = max(k10, 1d-40)
  k11  = max(k11, 1d-40)
  k9L  = max(k9L, 1d-40)
  k10L = max(k10L, 1d-40)
  k11L = max(k11L, 1d-40)

  !Log of all the rates
  lk8 = log10(k8)
  lk9 = (LTEfac*log10(k9L)) + (NLTEfac*log10(k9))
  lk10 = (LTEfac*log10(k10L)) + (NLTEfac*log10(k10))
  lk11 = (LTEfac*log10(k11L)) + (NLTEfac*log10(k11))

  rate = (ne*(10.d0**lk8)) + (nHI*(10.d0**lk9)) + (nH2*(10.d0**lk10)) + (nHeI*(10.0**lk11))
  rate = max(rate, 1d-40)

  rate = MAX(rate,1.d-100)

END FUNCTION beta_H2

FUNCTION alpha_CO(G0, xi_cr_H2, nCII, nH2, nOI) result(rate)
  ! see glover 2012
  implicit none

  real(dp), intent(in):: G0, xi_cr_H2, nCII, nH2, nOI
  real(dp):: rate
  real(dp):: k0, k1, gammaCHx_cr, gammaCHx, beta

  k0 = 5.d-16 ! cm^3 s^-1
  k1 = 5.d-10 ! Rate coefficient for the formation of CO from O + CHx

  ! Assuming CO formation is modulated by CH2+, https://home.strw.leidenuniv.nl/~ewine/photo/display_ch2+_65e31a07e69e64dbd64d37801983018f.html
  gammaCHx_cr = 8.88d-15 * (xi_cr_H2 / 1d-16) ! Cosmic rays
  gammaCHx = (1.41d-10 * G0) + gammaCHx_cr

  beta = k1 * nOI/(k1*nOI + gammaCHx)
  rate = k0 * nCII * nH2 * beta

  rate = MAX(rate,1.d-100)

END FUNCTION alpha_CO

FUNCTION beta_CO(G0, xi_cr_H2) result(rate)
  ! CO destruction
  implicit none

  real(dp), intent(in):: G0, xi_cr_H2
  real(dp):: rate
  real(dp):: gammaCO, gammaCO_cr

  gammaCO = 2.43d-10 * G0
  gammaCO_cr = 4.62d-15 * (xi_cr_H2 / 1d-16)

  rate = gammaCO + gammaCO_cr

  rate = MAX(rate,1.d-100)

END FUNCTION beta_CO

FUNCTION tau_dust_lw_mw(nH_tot, dx_SS, Z) result(tau)
  ! Dust optical depth in the Lyman-Werner band over dx_SS for MW dust scaled by Z (the dust-to-gas
  ! ratio in MW units), on the total hydrogen column (Draine & Bertoldi 1996)
  ! nH_tot : hydrogen nuclei number density, all states [cm^-3]
  implicit none
  real(dp), intent(in) :: nH_tot, dx_SS, Z
  real(dp) :: tau
  tau = sigma_dust_lw_mw * Z * nH_tot * dx_SS
END FUNCTION tau_dust_lw_mw

FUNCTION lw_transmission(tau, f_ani) result(trans)
  ! Mean-intensity transmission of a slab of normal optical depth tau for a field that is a beam
  ! (fraction f_ani, at normal incidence) plus an isotropic part: f e^-tau + (1-f) E2(tau). Exact in
  ! both limits, where exp(-(2-f) tau) only holds for tau << 1 (it is 4x too small at tau = 3).
  implicit none
  real(dp), intent(in) :: tau, f_ani
  real(dp) :: trans, f
  f = min(max(f_ani, 0d0), 1d0)
  trans = f * safe_exp(-tau) + (1d0 - f) * expint_E2(tau)
END FUNCTION lw_transmission

FUNCTION expint_E2(x) result(e2)
  ! Exponential integral E2(x) = exp(-x) - x E1(x), x >= 0; E1 by its series (x <= 1) or continued
  ! fraction (x > 1), to ~1e-12 (Numerical Recipes, expint)
  implicit none
  real(dp), intent(in) :: x
  real(dp) :: e2, e1, b, c, d, h, a, term, fact
  real(dp), parameter :: euler = 0.5772156649015329d0, eps = 1d-14
  integer :: k
  if (x <= 0d0) then
     e2 = 1d0
     return
  end if
  if (x > 700d0) then
     e2 = 0d0
     return
  end if
  if (x <= 1d0) then
     e1 = -euler - log(x)
     fact = 1d0
     do k = 1, 100
        fact = -fact * x / k
        term = -fact / k
        e1 = e1 + term
        if (abs(term) < eps * abs(e1)) exit
     end do
  else
     b = x + 1d0
     c = 1d0 / 1d-300
     d = 1d0 / b
     h = d
     do k = 1, 200
        a = -dble(k) * dble(k)
        b = b + 2d0
        d = 1d0 / (a * d + b)
        c = b + a / c
        h = h * c * d
        if (abs(c * d - 1d0) < eps) exit
     end do
     e1 = h * exp(-x)
  end if
  e2 = exp(-x) - x * e1
END FUNCTION expint_E2

FUNCTION comp_Sd(nH_tot, dx_SS, Z) result(ss_factor)
  ! Dust shielding of the Lyman-Werner band for MW dust scaled by Z, isotropic field
  ! (lw_transmission with f = 0); with CALIMA the per-bin grain cross-sections are used instead
  ! (compute_lw_dust_optical_depth in dust_photophysics.f90)
  implicit none
  real(dp), intent(in) :: nH_tot, dx_SS, Z
  real(dp) :: ss_factor
  ss_factor = lw_transmission(tau_dust_lw_mw(nH_tot, dx_SS, Z), 0d0)
END FUNCTION comp_Sd

FUNCTION comp_SH2(nH2, dx_SS, T) result(ss_factor)
  ! H2 self-shielding of its photodissociation: Wolcott-Green, Haiman & Bryan (2011, eq. 10), the
  ! Draine & Bertoldi (1996, eq. 37) form with exponent 1.1 instead of 2,
  !   f = 0.965/(1 + x/b5)^1.1 + 0.035/(1 + x)^0.5 exp(-8.5e-4 (1 + x)^0.5),
  ! x = N(H2)/5e14 cm^-2, b5 = b/(1 km/s) with the thermal Doppler parameter b = (2kT/m_H2)^1/2.
  ! (Replaces the Gnedin, Tassis & Kravtsov 2009 form with omega = 0.2 and b5 = 1, a calibration for
  ! unresolved clumping that gives 5.7x weaker shielding above N(H2) ~ 1e17 cm^-2.)
  ! nH2 : H2 number density [cm^-3]; dx_SS : shielding length [cm]; T : gas temperature [K]
  implicit none
  real(dp), intent(in):: nH2, dx_SS, T
  real(dp):: ss_factor
  real(dp):: x, b5
  real(dp), parameter :: kB_over_mH = 8.2504d7   ! k_B/m_H [cm^2 s^-2 K^-1], so b = (kT/m_H)^1/2 for H2
  x = nH2 * dx_SS / 5.d14
  b5 = sqrt(kB_over_mH * max(T, 1d0)) / 1d5
  ss_factor = 0.965d0 / (1.d0 + x / b5)**1.1d0 &
            + 0.035d0 / sqrt(1.d0 + x) * safe_exp(-8.5d-4 * sqrt(1.d0 + x))
END FUNCTION comp_SH2

FUNCTION comp_SCO(nco_mol, nh2, dx_SS) result(ss_factor)
  ! CO line shielding theta(N(CO), N(H2)) of Visser et al. (2009): self-shielding, and screening by
  ! H2 lines, from their 2D table (co_shield_visser09.inc), bilinear in log N and log theta; columns
  ! beyond the table take its edge values. Dust is applied separately (co_dust_to_lw).
  ! nco_mol, nh2 : CO and H2 number densities [cm^-3]; dx_SS : shielding length [cm]
  implicit none
  real(dp), intent(in):: nco_mol, nh2, dx_SS
  real(dp):: ss_factor
  real(dp):: wco, wh2, t00, t10, t01, t11
  integer :: ic, ih
  call visser_index(nco_mol * dx_SS, nco_grid_visser, nco_visser, ic, wco)
  call visser_index(nh2 * dx_SS, nh2_grid_visser, nh2_visser, ih, wh2)
  t00 = log10(theta_visser(ic, ih));       t10 = log10(theta_visser(ic+1, ih))
  t01 = log10(theta_visser(ic, ih+1));     t11 = log10(theta_visser(ic+1, ih+1))
  ss_factor = 10d0**((1d0-wco)*(1d0-wh2)*t00 + wco*(1d0-wh2)*t10 + (1d0-wco)*wh2*t01 + wco*wh2*t11)
contains
  subroutine visser_index(col, grid, n, i, w)
     ! lower index i and weight w in log col on grid(1:n); grid(1) stands for 0 (below grid(2): i = 1,
     ! weight on log N between grid(1) and grid(2)); clamped at the top
     real(dp), intent(in) :: col, grid(:)
     integer, intent(in) :: n
     integer, intent(out) :: i
     real(dp), intent(out) :: w
     real(dp) :: lcol
     integer :: k
     if (col <= grid(1)) then
        i = 1; w = 0d0
        return
     end if
     if (col >= grid(n)) then
        i = n - 1; w = 1d0
        return
     end if
     i = 1
     do k = 1, n - 1
        if (col >= grid(k)) i = k
     end do
     lcol = log10(col)
     w = (lcol - log10(grid(i))) / (log10(grid(i+1)) - log10(grid(i)))
  end subroutine visser_index
END FUNCTION comp_SCO


end module molecules_module
