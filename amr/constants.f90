module constants
  use amr_commons, ONLY: dp

  ! Numerical constants
  real(dp),parameter ::twopi        = 6.2831853d0
  real(dp),parameter ::pi           = twopi/2d0
  real(dp),parameter ::sq2pi        = sqrt(twopi)

  real(dp),parameter ::mu_mol       = 1.2195d0

  ! Physical constants
  ! Source:
  ! * SI - SI Brochure (2018)
  ! * PCAD - http://www.astro.wisc.edu/~dolan/constants.html
  ! * NIST - National Institute of Standards and Technology
  ! * IAU - Internatonal Astronomical Union resolution
  real(dp),parameter ::hplanck      = 6.6260702d-27 ! Planck const. [erg s]; SI
  real(dp),parameter ::eV2erg       = 1.6021766d-12 ! Electronvolt [erg]; SI
  real(dp),parameter ::kB           = 1.3806490d-16 ! Boltzmann const. [erg K-1]; SI
  real(dp),parameter ::sb           = 5.6703744d-05 ! Stefan-Boltzmann const. [erg cm−2 s−1 K−4]
  real(dp),parameter ::c_cgs        = 2.9979246d+10 ! Speed of light [cm s-1]; SI
  real(dp),parameter ::a_r          = 7.5657233d-15 ! Radiation density const. [erg cm-3 K-4]; SI (derived)
  real(dp),parameter ::factG_in_cgs = 6.6740800d-08 ! Gravitational const. [cm3 g-1 s-2]; NIST
  real(dp),parameter ::sigma_T      = 6.6524587d-25 ! Thomson scattering cross-section [cm2]; NIST
  real(dp),parameter ::M_sun        = 1.9891000d+33 ! Solar Mass [g]; IAU
  real(dp),parameter ::L_sun        = 3.8280000d+33 ! Solar Lum [erg s-1]; IAU
  real(dp),parameter ::rhoc         = 1.8800000d-29 ! Crit. density [g cm-3]
  real(dp),parameter ::one_over_clight=3.335640484668562d-11 ! Save some computation
  real(dp),parameter ::four_pi_eps0 = 1.1126501d-19 ! 4 pi x vacuum permittivity [F nm-1]
  real(dp),parameter ::e2instatC    = 2.3070775d-17 ! Elemental charge squared [statC^2]

  ! Dalton or unified atomic mass unit (amu) is the mass of 1/12 of a carbon-12 atom
  ! Since we do not track multiple isotopes, we use what is termed the
  ! "standard atomic weight" as given in the IUPAC 2021 convention report:
  ! https://www.degruyterbrill.com/document/doi/10.1515/pac-2019-0603/html
  real(dp),parameter ::amu2g        = 1.6605390d-24 ! Atomic mass unit [g]; IUPAC
  real(dp),parameter ::mH_amu       = 1.0080d0      ! H atom mass [amu]; IUPAC
  real(dp),parameter ::mHe_amu      = 4.0026d0      ! He atom mass [amu]; IUPAC
  real(dp),parameter ::mC_amu       = 12.0107d0     ! C atom mass [amu]; IUPAC
  real(dp),parameter ::mN_amu       = 14.0067d0     ! N atom mass [amu]; IUPAC
  real(dp),parameter ::mO_amu       = 15.9994d0     ! O atom mass [amu]; IUPAC
  real(dp),parameter ::mNe_amu      = 20.1797d0     ! Ne atom mass [amu]; IUPAC
  real(dp),parameter ::mMg_amu      = 24.305d0      ! Mg atom mass [amu]; IUPAC
  real(dp),parameter ::mSi_amu      = 28.0855d0     ! Si atom mass [amu]; IUPAC
  real(dp),parameter ::mS_amu       = 32.065d0      ! S atom mass [amu]; IUPAC
  real(dp),parameter ::mFe_amu      = 55.854d0      ! Fe atom mass [amu]; IUPAC
  real(dp),parameter ::mCO_amu      = 28.0101d0     ! CO molecule mass [amu]; NIST
  real(dp),parameter ::mH           = 1.6738233d-24 ! H atom mass [g]
  real(dp),parameter ::mCO          = 4.6511863d-23 ! CO molecule mass [g]

  ! Conversion factors - distance
  ! IAU 2012 convention:
  ! 1 pc = 648000 AU / pi
  ! 1 AU = 14 959 787 070 000 cm
  real(dp),parameter ::pc2cm        = 3.0856776d+18
  real(dp),parameter ::kpc2cm       = 3.0856776d+21
  real(dp),parameter ::Mpc2cm       = 3.0856776d+24
  real(dp),parameter ::Gpc2cm       = 3.0856776d+27

  ! Conversion factors - time
  ! Year definition follows IAU recommendation
  ! https://www.iau.org/publications/proceedings_rules/units/
  real(dp),parameter ::yr2sec       = 3.15576000d+07 ! Year [s]
  real(dp),parameter ::kyr2sec      = 3.15576000d+10 ! Kyr [s]
  real(dp),parameter ::Myr2sec      = 3.15576000d+13 ! Myr [s]
  real(dp),parameter ::Gyr2sec      = 3.15576000d+16 ! Gyr [s]
  real(dp),parameter ::sec2Gyr      = 3.16880878d-17 ! sec [Gyr]



end module constants
