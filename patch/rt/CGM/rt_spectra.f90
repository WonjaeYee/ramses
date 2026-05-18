! These are modules for reading, integrating and interpolating RT-relevant
! values from spectral tables, specifically SED (spectral energy
! distrbution) tables for stellar particle sources, as functions of age
! and metallicity, and UV-spectrum tables, as functions of redshift.
! The modules are:
! spectrum_integrator_module
!    For dealing with the integration of spectra in general
! SED_module
!    For dealing with the SED data
! UV_module
!    For dealing with the UV data
!_________________________________________________________________________

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for integrating a wavelength-dependent spectrum
!_________________________________________________________________________
!
MODULE spectrum_integrator_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  implicit none

  PUBLIC integrateSpectrum, f1, fLambda, fdivLambda, fSig, fSigLambda,   &
         fSigdivLambda, trapz1, &
         integrateSpectrum_metals, f1_metals, fLambda_metals, &
         fdivLambda_metals, fSig_metals, fSigLambda_metals,   &
         fSigdivLambda_metals, fLambda_dust, fSigLambda_dust

  PRIVATE   ! default

CONTAINS

!*************************************************************************
FUNCTION integrateSpectrum(X, Y, N, e0, e1, species, func, doPrint)

! Integrate spectral weighted function in energy interval [e0,e1]
! X      => Wavelengths [angstrom]
! Y      => Spectral luminosity per angstrom at wavelenghts [XX A-1]
! N      => Length of X and Y
! e0,e1  => Integrated interval [ev]
! species=> ion species, used as an argument in fx
! func   => Function which is integrated (of X, Y, species)
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs,eV_to_erg, hp
  real(kind=8):: integrateSpectrum, X(N), Y(N), e0, e1
  integer :: N, species
  interface
     real(kind=8) function func(wavelength,intensity,species)
       use amr_parameters,only:dp
       real(kind=8)::wavelength,intensity
       integer::species
     end function func
  end interface!----------------------------------------------------------
  real(kind=8),dimension(:),allocatable:: xx, yy, f
  real(dp):: la0, la1
  integer :: i
  logical,optional::doPrint
!-------------------------------------------------------------------------
  integrateSpectrum=0.
  if(N .le. 2) RETURN
  ! Convert energy interval to wavelength interval
  la0 = X(1) ; la1 = X(N)
  if(e1.gt.0) la0 = max(la0, 1.d8 * hp * c_cgs / e1 / eV_to_erg)
  if(e0.gt.0) la1 = min(la1, 1.d8 * hp * c_cgs / e0 / eV_to_erg)
  if(la0 .ge. la1) RETURN
  ! If we get here, the [la0, la1] inverval is completely within X
  allocate(xx(N)) ; allocate(yy(N)) ; allocate(f(N))
  xx =  la0   ;   yy =  0.   ;   f = 0.
  i=2
  do while ( i.lt.N .and. X(i).le.la0 )
     i = i+1                      !              Below wavelength interval
  enddo                           !   X(i) is now the first entry .gt. la0
  ! Interpolate to value at la0
  yy(i-1) = Y(i-1) + (xx(i-1)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i-1)  = func(xx(i-1), yy(i-1), species)
  do while ( i.lt.N .and. X(i).le.la1 )              ! Now within interval
     xx(i) = X(i) ; yy(i) = Y(i) ; f(i) = func(xx(i),yy(i),species)
     i = i+1
  enddo                          ! i=N or X(i) is the first entry .gt. la1
  xx(i:) = la1                   !             Interpolate to value at la1
  yy(i) = Y(i-1) + (xx(i)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i)  = func(xx(i),yy(i),species)

  !if(present(doPrint)) then
  !   if(doprint) then
  !      write(*,*) e0,e1,la0,la1,N
  !      write(*,*) '***************'
  !      do i=1,N
  !         write(*,*) xx(i),f(i),yy(i)
  !      end do
  !      stop
  !   endif
  !endif

  integrateSpectrum = trapz1(xx,f,i)
  deallocate(xx) ; deallocate(yy) ; deallocate(f)

END FUNCTION integrateSpectrum
!*************************************************************************
FUNCTION integrateSpectrum_metals(X, Y, N, e0, e1, species, ion_state, func, doPrint)

   ! Integrate spectral weighted function in energy interval [e0,e1]
   ! X      => Wavelengths [angstrom]
   ! Y      => Spectral luminosity per angstrom at wavelenghts [XX A-1]
   ! N      => Length of X and Y
   ! e0,e1  => Integrated interval [ev]
   ! species=> ion species, used as an argument in fx
   ! func   => Function which is integrated (of X, Y, species)
   !-------------------------------------------------------------------------
      use rt_parameters,only:c_cgs,eV_to_erg, hp
      use amr_commons,only:myid
      real(kind=8):: integrateSpectrum_metals, X(N), Y(N), e0, e1
      integer :: N, species, ion_state
      interface
         real(kind=8) function func(wavelength,intensity,species,ion_state)
            use amr_parameters,only:dp
            real(kind=8)::wavelength,intensity
            integer::species, ion_state
         end function func
      end interface!----------------------------------------------------------
      real(kind=8),dimension(:),allocatable:: xx, yy, f
      real(dp):: la0, la1
      integer :: i
      logical,optional::doPrint
   !-------------------------------------------------------------------------
      integrateSpectrum_metals=0.
      if(N .le. 2) RETURN
      ! Convert energy interval to wavelength interval
      la0 = X(1) ; la1 = X(N)
      if(e1.gt.0) la0 = max(la0, 1.d8 * hp * c_cgs / e1 / eV_to_erg)
      if(e0.gt.0) la1 = min(la1, 1.d8 * hp * c_cgs / e0 / eV_to_erg)
      if(la0 .ge. la1) RETURN
      ! If we get here, the [la0, la1] inverval is completely within X
      allocate(xx(N)) ; allocate(yy(N)) ; allocate(f(N))
      xx =  la0   ;   yy =  0.   ;   f = 0.
      i=2
      do while ( i.lt.N .and. X(i).le.la0 )
         i = i+1                      !              Below wavelength interval
      enddo                           !   X(i) is now the first entry .gt. la0
      ! Interpolate to value at la0
      yy(i-1) = Y(i-1) + (xx(i-1)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
      f(i-1)  = func(xx(i-1), yy(i-1), species, ion_state)
      do while ( i.lt.N .and. X(i).le.la1 )              ! Now within interval
         xx(i) = X(i) ; yy(i) = Y(i) ; f(i) = func(xx(i),yy(i),species,ion_state)
         i = i+1
      enddo                          ! i=N or X(i) is the first entry .gt. la1
      xx(i:) = la1                   !             Interpolate to value at la1
      yy(i) = Y(i-1) + (xx(i)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
      f(i)  = func(xx(i),yy(i),species,ion_state)

      integrateSpectrum_metals = trapz1(xx,f,i)
      deallocate(xx) ; deallocate(yy) ; deallocate(f)

END FUNCTION integrateSpectrum_metals
!*************************************************************************
! FUNCTIONS FOR USE WITH integrateSpectrum:
! lambda  => wavelengths in Angstrom
! f       => function of wavelength (a spectrum in some units)
! species => 1=HI, 2=HeI or 3=HeII
!_________________________________________________________________________
FUNCTION f1(lambda, f, species)
  real(kind=8):: f1, lambda, f
  integer :: species
  f1 = f
END FUNCTION f1

FUNCTION fLambda(lambda, f, species)
  real(kind=8):: fLambda, lambda, f
  integer :: species
  fLambda = f * lambda
END FUNCTION fLambda

FUNCTION fdivLambda(lambda, f, species)
  real(kind=8):: fdivlambda, lambda, f
  integer :: species
  fdivLambda = f / lambda
END FUNCTION fdivLambda

FUNCTION fSig(lambda, f, species)
  real(kind=8):: fSig, lambda, f
  integer :: species
  fSig = f * getCrosssection_Hui(lambda,species)
END FUNCTION fSig

FUNCTION fSigLambda(lambda, f, species)
  real(kind=8):: fSigLambda, lambda, f
  integer :: species
  fSigLambda = f * lambda * getCrosssection_Hui(lambda,species)
END FUNCTION fSigLambda

FUNCTION fSigdivLambda(lambda, f, species)
  real(kind=8):: fSigdivLambda, lambda, f
  integer :: species
  fSigdivLambda = f / lambda * getCrosssection_Hui(lambda,species)
END FUNCTION fSigdivLambda

FUNCTION f1_metals(lambda, f, species, ion_state)
   real(kind=8):: f1_metals, lambda, f
   integer :: species, ion_state
   f1_metals = f
END FUNCTION f1_metals

FUNCTION fLambda_metals(lambda, f, species, ion_state)
   real(kind=8):: fLambda_metals, lambda, f
   integer :: species, ion_state
   fLambda_metals = f * lambda
END FUNCTION fLambda_metals

FUNCTION fdivLambda_metals(lambda, f, species, ion_state)
   real(kind=8):: fdivlambda_metals, lambda, f
   integer :: species, ion_state
   fdivLambda_metals = f / lambda
END FUNCTION fdivLambda_metals

FUNCTION fSig_metals(lambda, f, species, ion_state)
   real(kind=8):: fSig_metals, lambda, f
   integer :: species, ion_state
   if (species.eq.1) then ! Oxygen
      fSig_metals = f * getCrosssection_oxygen(lambda,ion_state)
   else if (species.eq.2) then ! Nitrogen
      fSig_metals = f * getCrosssection_nitrogen(lambda,ion_state)
   else if (species.eq.3) then ! Carbon
      fSig_metals = f * getCrosssection_carbon(lambda,ion_state)
   else if (species.eq.4) then ! Magnesium
      fSig_metals = f * getCrosssection_magnesium(lambda,ion_state)
   else if (species.eq.5) then ! Silicon
      fSig_metals = f * getCrosssection_silicon(lambda,ion_state)
   else if (species.eq.6) then ! Sulfur
      fSig_metals = f * getCrosssection_sulfur(lambda,ion_state)
   else if (species.eq.7) then ! Iron
      fSig_metals = f * getCrosssection_iron(lambda,ion_state)
   else if (species.eq.8) then ! Neon
      fSig_metals = f * getCrosssection_neon(lambda,ion_state)
   endif
END FUNCTION fSig_metals

FUNCTION fSigLambda_metals(lambda, f, species, ion_state)
   real(kind=8):: fSigLambda_metals, lambda, f
   integer :: species, ion_state
   if (species.eq.1) then ! Oxygen
      fSigLambda_metals = f * lambda * getCrosssection_oxygen(lambda,ion_state)
   else if (species.eq.2) then ! Nitrogen
      fSigLambda_metals = f * lambda * getCrosssection_nitrogen(lambda,ion_state)
   else if (species.eq.3) then ! Carbon
      fSigLambda_metals = f * lambda * getCrosssection_carbon(lambda,ion_state)
   else if (species.eq.4) then ! Magnesium
      fSigLambda_metals = f * lambda * getCrosssection_magnesium(lambda,ion_state)
   else if (species.eq.5) then ! Silicon
      fSigLambda_metals = f * lambda * getCrosssection_silicon(lambda,ion_state)
   else if (species.eq.6) then ! Sulfur
      fSigLambda_metals = f * lambda * getCrosssection_sulfur(lambda,ion_state)
   else if (species.eq.7) then ! Iron
      fSigLambda_metals = f * lambda * getCrosssection_iron(lambda,ion_state)
   else if (species.eq.8) then ! Neon
      fSigLambda_metals = f * lambda * getCrosssection_neon(lambda,ion_state)
   endif
END FUNCTION fSigLambda_metals

FUNCTION fSigdivLambda_metals(lambda, f, species, ion_state)
   real(kind=8):: fSigdivLambda_metals, lambda, f
   integer :: species, ion_state
   if (species.eq.1) then ! Oxygen
      fSigdivLambda_metals = f / lambda * getCrosssection_oxygen(lambda,ion_state)
   else if (species.eq.2) then ! Nitrogen
      fSigdivLambda_metals = f / lambda * getCrosssection_nitrogen(lambda,ion_state)
   else if (species.eq.3) then ! Carbon
      fSigdivLambda_metals = f / lambda * getCrosssection_carbon(lambda,ion_state)
   else if (species.eq.4) then ! Magnesium
      fSigdivLambda_metals = f / lambda * getCrosssection_magnesium(lambda,ion_state)
   else if (species.eq.5) then ! Silicon
      fSigdivLambda_metals = f / lambda * getCrosssection_silicon(lambda,ion_state)
   else if (species.eq.6) then ! Sulfur
      fSigdivLambda_metals = f / lambda * getCrosssection_sulfur(lambda,ion_state)
   else if (species.eq.7) then ! Iron
      fSigdivLambda_metals = f / lambda * getCrosssection_iron(lambda,ion_state)
   else if (species.eq.8) then ! Neon
      fSigdivLambda_metals = f / lambda * getCrosssection_neon(lambda,ion_state)
   endif
END FUNCTION fSigdivLambda_metals

FUNCTION fLambda_dust(lambda, f, species)
   real(kind=8):: fLambda_dust, lambda, f
   integer :: species
   fLambda_dust = f * lambda
 END FUNCTION fLambda_dust

 FUNCTION fSigLambda_dust(lambda, f, species)
   real(kind=8):: fSigLambda_dust, lambda, f
   integer :: species
   fSigLambda_dust = f * lambda * getCrosssection_BARE_GR_S_DUST(lambda,species)
 END FUNCTION fSigLambda_dust
!_________________________________________________________________________

!*************************************************************************
FUNCTION trapz1(X,Y,N,cum)

! Integrates function Y(X) along the whole interval 1..N, using a very
! simple staircase method and returns the result.
! Optionally, the culumative integral is returned in the cum argument.
!-------------------------------------------------------------------------
  integer :: N,i
  real(kind=8):: trapz1
  real(kind=8):: X(N),Y(N)
  real(kind=8),optional::cum(N)
  real(kind=8),allocatable::cumInt(:)
!-------------------------------------------------------------------------
  allocate(cumInt(N))
  cumInt(:)=0.d0
  if (N.le.1) RETURN
  do i=2,N
     cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2.d0
  end do
  trapz1 = cumInt(N)
  if(present(cum)) cum=cumInt
  deallocate(cumInt)
END FUNCTION trapz1

!*************************************************************************
FUNCTION getCrosssection_Hui(lambda, species)

! Gives an atom-photon cross-section of given species at given wavelength,
! as given by Hui and Gnedin (1997).
! lambda  => Wavelength in angstrom
! species => 1=HI, 2=HeI or 3=HeII
! returns :  photoionization cross-section in cm^2
!------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp, ixHI, ixHII      &
                        ,ixHeII, ixHeIII, isH2Katz
  real(kind=8)      :: lambda, getCrosssection_Hui
  integer           :: species
  real(kind=8)      :: E, E0, cs0, P, ya, yw, y0, y1, x, y
!------------------------------------------------------------------------
  E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

  if (species .eq. ixHI) then !H2
     if(isH2Katz)then
        ! See table 1 in Baczynski, Glove and Klessen (2015)
        getCrosssection_Hui = 0.
        if (E .gt. 11.20 .and. E .lt. 13.59) getCrosssection_Hui = 2.47d-18
        if (E .ge. 13.59 .and. E .le. 15.21) getCrosssection_Hui = 0.0d0
        if (E .gt. 15.21 .and. E .le. 15.45) getCrosssection_Hui = 0.09d-18
        if (E .gt. 15.70 .and. E .le. 15.95) getCrosssection_Hui = 1.15d-18
        if (E .gt. 15.95 .and. E .le. 16.20) getCrosssection_Hui = 3.00d-18
        if (E .gt. 16.20 .and. E .le. 16.40) getCrosssection_Hui = 5.00d-18
        if (E .gt. 16.40 .and. E .le. 16.65) getCrosssection_Hui = 6.75d-18
        if (E .gt. 16.65 .and. E .le. 16.85) getCrosssection_Hui = 8.00d-18
        if (E .gt. 16.85 .and. E .le. 17.00) getCrosssection_Hui = 9.00d-18
        if (E .gt. 17.00 .and. E .le. 17.20) getCrosssection_Hui = 9.50d-18
        if (E .gt. 17.20 .and. E .le. 17.65) getCrosssection_Hui = 9.80d-18
        if (E .gt. 17.65 .and. E .le. 18.10) getCrosssection_Hui = 10.10d-18
        if (E .gt. 18.10) getCrosssection_Hui = (10.10d-18) * ((18.10d0/E)**3.0d0)
     else
        getCrosssection_Hui = 0.
        if(E .lt. (ionEvs(ixHI)+1.)) getCrosssection_Hui = 3.0d-22
     endif
     RETURN
  endif

  if ( E .lt. ionEvs(species) ) then            ! below ionization energy
     getCrosssection_Hui=0.
     RETURN
  endif
  if (species .eq. ixHII)then !HII
     E0 = 4.298d-1 ; cs0 = 5.475d-14    ; P  = 2.963
     ya = 32.88    ; yw  = 0            ; y0 = 0         ; y1 = 0
  endif
  if (species .eq. ixHeII) then ! HeI
     E0 = 1.361d1    ; cs0 = 9.492d-16  ; P  = 3.188
     ya = 1.469      ; yw  = 2.039      ; y0 = 0.4434    ; y1 = 2.136
  endif
  if (species .eq. ixHeIII) then ! HeII
     E0 = 1.720      ; cs0 = 1.369d-14  ; P  = 2.963
     ya = 32.88      ; yw  = 0          ; y0 = 0         ; y1 = 0
  endif

  x = E/E0 - y0
  y = sqrt(x**2+y1**2)

  getCrosssection_Hui = &
       cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P

END FUNCTION getCrosssection_Hui

! Metal cross sections
FUNCTION getCrosssection_oxygen(lambda, N)
   ! Returns the cross section in cm^2 for oxygen
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_oxygen
   integer           :: N
   real(kind=8),dimension(1:8) :: E_th = (/ 1.362d1, 3.512d1, 5.494d1, 7.741d1, 1.139d2, 1.381d2, 7.393d2, 8.714d2 /)
   real(kind=8),dimension(1:8) :: E_max = (/ 5.380d2, 5.581d2, 5.840d2, 6.144d2, 6.491d2, 6.837d2, 5.000d4, 5.000d4 /)
   real(kind=8),dimension(1:8) :: E_0 = (/ 1.240d0, 1.386d0, 1.723d-1, 2.044d-1, 2.854d0, 7.824d0, 8.709d0, 2.754d1 /)
   real(kind=8),dimension(1:8) :: sig_0 = (/ 1.745e3, 5.967e1, 6.753e2, 8.659e-1, 1.642e4, 6.864e1, 1.329e2, 8.554e2 /)
   real(kind=8),dimension(1:8) :: y_a = (/ 3.784d0, 3.175d1, 3.852d2, 4.931d2, 1.792d0, 3.210d1, 2.535d1, 3.288d1 /)
   real(kind=8),dimension(1:8) :: P = (/ 1.764d1, 8.943d0, 6.822d0, 8.785d0, 2.647d1, 5.495d0, 2.336d0, 2.963d0 /)
   real(kind=8),dimension(1:8) :: y_w = (/ 7.589d-2, 1.934d-2, 1.191d-1, 3.143d0, 2.836d1, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:8) :: y_0 = (/ 8.698d0, 2.131d1, 3.839d-3, 3.328d2, 3.036d-2, 0d0, 0d0, 0d0 /)
   real(kind=8),dimension(1:8) :: y_1 = (/ 1.271d-1, 1.503d-2, 4.569d-1, 4.285d1, 5.554d-2, 0d0, 0d0, 0d0 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_oxygen = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_oxygen = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_oxygen = 0.d0
   endif

END FUNCTION getCrosssection_oxygen

FUNCTION getCrosssection_nitrogen(lambda, N)
   ! Returns the cross section in cm^2 for nitrogen
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_nitrogen
   integer           :: N
   real(kind=8),dimension(1:7) :: E_th = (/ 1.453d1, 2.960d1, 4.745d1, 7.747d1, 9.789d1, 5.521d2, 6.671d2 /)
   real(kind=8),dimension(1:7) :: E_max = (/ 4.048d2, 4.236d2, 4.473d2, 5.753d2, 5.043d2, 5.0d4, 5.0d4 /)
   real(kind=8),dimension(1:7) :: E_0 = (/ 4.034d0, 6.128d-2, 2.420d-1, 5.494d0, 4.471d0, 6.943d1, 2.108d1 /)
   real(kind=8),dimension(1:7) :: sig_0 = (/ 8.235d2, 1.944d0, 9.375d-1, 1.690d4, 8.376d1, 1.519d2, 1.117d3 /)
   real(kind=8),dimension(1:7) :: y_a = (/ 8.033d1, 8.163d2, 2.788d2, 1.714d0, 3.297d1, 2.627d1, 3.288d1 /)
   real(kind=8),dimension(1:7) :: P = (/ 3.928d0, 8.773d0, 9.156d0, 1.706d1, 6.003d0, 2.315d0, 2.963d0 /)
   real(kind=8),dimension(1:7) :: y_w = (/ 9.097d-2, 1.043d1, 1.850d0, 7.904d0, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:7) :: y_0 = (/ 8.598d-1, 4.280d2, 1.877d2, 6.415d-3, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:7) :: y_1 = (/ 2.325d0, 2.030d1, 3.999d0, 1.937d-2, 0.d0, 0.d0, 0.d0 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_nitrogen = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_nitrogen = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_nitrogen = 0.d0
   endif

END FUNCTION getCrosssection_nitrogen

FUNCTION getCrosssection_carbon(lambda, N)
   ! Returns the cross section in cm^2 for carbon
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_carbon
   integer           :: N
   real(kind=8),dimension(1:6) :: E_th = (/ 1.126d1, 2.438d1, 4.789d1, 6.449d1, 3.921d2, 4.900d2 /)
   real(kind=8),dimension(1:6) :: E_max = (/ 2.910d2, 3.076d2, 3.289d2, 3.522d2, 5.000d4, 5.000d4 /)
   real(kind=8),dimension(1:6) :: E_0 = (/ 2.144d0, 4.058d-1, 4.614d0, 3.506d0, 4.624d1, 1.548d1 /)
   real(kind=8),dimension(1:6) :: sig_0 = (/ 5.027d2, 8.709d0, 1.539d4, 1.068d2, 2.344d2, 1.521d3 /)
   real(kind=8),dimension(1:6) :: y_a = (/ 6.216d1, 1.261d2, 1.737d0, 1.436d1, 2.183d1, 3.288d1 /)
   real(kind=8),dimension(1:6) :: P = (/ 5.101d0, 8.578d0, 1.593d1, 7.457d0, 2.581d0, 2.963d0 /)
   real(kind=8),dimension(1:6) :: y_w = (/ 9.157d-2, 2.093d0, 5.922d0, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:6) :: y_0 = (/ 1.133d0, 4.929d1, 4.378d-3, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:6) :: y_1 = (/ 1.607d0, 3.234d0, 2.528d-2, 0.d0, 0.d0, 0.d0 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_carbon = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_carbon = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_carbon = 0.d0
   endif

END FUNCTION getCrosssection_carbon

FUNCTION getCrosssection_magnesium(lambda, N)
   ! Returns the cross section in cm^2 for magnesium
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_magnesium
   integer           :: N
   real(kind=8),dimension(1:12) :: E_th = (/ 7.646d0, 1.504d1, 8.014d1, 1.093d2, 1.413d2, 1.865d2, 2.249d2, 2.660d2, 3.282d2, 3.675d2, 1.762d3, 1.963d3 /)
   real(kind=8),dimension(1:12) :: E_max = (/ 5.490d1, 6.569d1, 1.317d3, 1.356d3, 1.400d3, 1.449d3, 1.503d3, 1.558d3, 1.618d3, 1.675d3, 5.0d4, 5.0d4 /)
   real(kind=8),dimension(1:12) :: E_0 = (/ 1.197d1, 8.139d0, 1.086d1, 2.912d1, 9.762d-1, 1.711d0, 3.570d0, 4.884d-1, 3.482d1, 1.452d1, 2.042d2, 6.203d1 /)
   real(kind=8),dimension(1:12) :: sig_0 = (/ 1.372d8, 3.278d0, 5.377d2, 1.394d3, 1.728d0, 2.185d0, 3.104d0, 6.344d-2, 9.008d2, 4.427d1, 6.140d1, 3.802d2 /)
   real(kind=8),dimension(1:12) :: y_a = (/ 2.228d-1, 4.241d7, 9.779d0, 2.895d0, 9.184d1, 9.350d1, 6.060d1, 5.085d2, 1.823d0, 3.826d1, 2.778d1, 3.288d1 /)
   real(kind=8),dimension(1:12) :: P = (/ 1.574d1, 3.610d0, 7.117d0, 6.487d0, 1.006d1, 9.202d0, 8.857d0, 9.385d0, 1.444d1, 5.460d0, 2.161d0, 2.963d0 /)
   real(kind=8),dimension(1:12) :: y_w = (/ 2.805d-1, 0.d0, 2.604d0, 4.326d-2, 8.090d-1, 6.325d-1, 1.422d0, 6.666d-1, 2.751d0, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:12) :: y_0 = (/ 0.d0, 0.d0, 4.860d0, 9.402d-1, 1.276d2, 1.007d2, 5.452d1, 5.348d2, 5.444d0, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:12) :: y_1 = (/ 0.d0, 0.d0, 3.722d0, 1.135d-1, 3.979d0, 1.729d0, 2.078d0, 3.997d-3, 7.918d-2, 0.d0, 0.d0, 0.d0 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_magnesium = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_magnesium = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_magnesium = 0.d0
   endif

END FUNCTION getCrosssection_magnesium

FUNCTION getCrosssection_silicon(lambda, N)
   ! Returns the cross section in cm^2 for silicon
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_silicon
   integer           :: N
   real(kind=8),dimension(1:14) :: E_th = (/ 8.152d0, 1.635d1, 3.349d1, 4.514d1, 1.668d2, 2.051d2, 2.465d2, 3.032d2, 3.511d2, 4.014d2, 4.761d2, 5.235d2, 2.438d3, 2.673d3 /)
   real(kind=8),dimension(1:14) :: E_max = (/ 1.060d2, 1.186d2, 1.311d2, 1.466d2, 1.887d3, 1.946d3, 2.001d3, 2.058d3, 2.125d3, 2.194d3, 2.268d3, 2.336d3, 5.0d4, 5.0d4 /)
   real(kind=8),dimension(1:14) :: E_0 = (/ 2.317d1, 2.556d0, 1.659d-1, 1.288d1, 7.761d-1, 6.305d1, 3.277d-1, 7.655d-1, 3.343d-1, 8.787d-1, 1.205d1, 3.560d1, 2.752d2, 8.447d1 /)
   real(kind=8),dimension(1:14) :: sig_0 = (/ 2.506d1, 4.140d0, 5.790d-4, 6.083d0, 8.863d-1, 7.293d1, 6.680d-2, 3.477d-1, 1.465d-1, 1.950d-1, 1.992d4, 2.539d1, 4.754d1, 2.793d2 /)
   real(kind=8),dimension(1:14) :: y_a = (/ 2.057d1, 1.337d1, 1.474d2, 1.356d6, 1.541d2, 1.558d2, 4.132d1, 3.733d2, 1.404d3, 7.461d2, 1.582d0, 3.307d1, 2.848d1, 3.288d1 /)
   real(kind=8),dimension(1:14) :: P = (/ 3.546d0, 1.191d1, 1.336d1, 3.353d0, 9.980d0, 2.400d0, 1.606d1, 8.986d0, 8.503d0, 8.302d0, 2.425d1, 4.728d0, 2.135d0, 2.963d0 /)
   real(kind=8),dimension(1:14) :: y_w = (/ 2.837d-1, 1.570d0, 8.626d-1, 0.d0, 1.303d0, 2.989d-3, 3.280d0, 1.476d-3, 1.646d0, 4.489d-1, 2.392d1, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:14) :: y_0 = (/ 1.672d-5, 6.634d0, 9.613d1, 0.d0, 2.009d2, 1.115d0, 1.149d-2, 3.850d2, 1.036d3, 4.528d2, 1.990d-2, 0.d0, 0.d0, 0.d0 /)
   real(kind=8),dimension(1:14) :: y_1 = (/ 4.207d-1, 1.272d-1, 6.442d-1, 0.d0, 4.537d0, 8.051d-2, 6.396d-1, 8.999d-2, 2.936d-1, 1.015d0, 1.007d-2, 0.d0, 0.d0, 0.d0 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_silicon = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_silicon = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_silicon = 0.d0
   endif

END FUNCTION getCrosssection_silicon

FUNCTION getCrosssection_sulfur(lambda, N)
   ! Returns the cross section in cm^2 for sulfur
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_sulfur
   integer           :: N
   real(kind=8),dimension(1:16) :: E_th = (/ 10.36d0, 23.33d0, 34.83d0, 47.31d0, 72.68d0, &
                                             88.05d0, 280.9d0, 328.2d0, 379.1d0, 447.1d0, &
                                             504.8d0, 564.7d0, 651.7d0, 707.2d0, 3224.d0, &
                                             3494.d0 /)
   real(kind=8),dimension(1:16) :: E_max = (/ 170.d0, 184.6d0, 199.5d0, 216.4d0, 235.d0, &
                                              255.7d0, 2569.d0, 2641.d0, 2705.d0, 2782.d0, &
                                              2859.d0, 2941.d0, 3029.d0, 3107.d0, 50000.d0, &
                                              50000.d0 /)
   real(kind=8),dimension(1:16) :: E_0 = (/ 1.808d+01, 8.787d+00, 2.027d+00, 2.173d+00, 1.713d-01, &
                                            1.413d+01, 3.757d-01, 1.462d+01, 1.526d-01, 1.040d+01, &
                                            6.485d+00, 2.443d+00, 1.474d+01, 3.310d+01, 4.390d+02, &
                                            1.104d+02 /)
   real(kind=8),dimension(1:16) :: sig_0 = (/ 4.564d+04, 3.136d+02, 6.666d+00, 2.606d+00, 5.072d-04, &
                                              9.139d+00, 5.703d-01, 3.161d+01, 9.646d+03, 5.364d+01, &
                                              1.275d+01, 3.490d-01, 2.294d+04, 2.555d+01, 2.453d+01, &
                                              2.139d+02 /)
   real(kind=8),dimension(1:16) :: y_a = (/ 1.000d+00, 3.442d+00, 5.454d+01, 6.641d+01, 1.986d+02, &
                                            1.656d+03, 1.460d+02, 1.611d+01, 1.438d+03, 3.641d+01, &
                                            6.583d+01, 5.411d+02, 1.529d+00, 3.821d+01, 4.405d+01, &
                                            3.288d+01 /)
   real(kind=8),dimension(1:16) :: P = (/ 13.61d0, 12.81d0, 8.611d0, 8.655d0, 13.07d0, 3.626d0, &
                                          11.35d0, 8.642d0, 5.977d0, 7.09d0, 7.692d0, 7.769d0, &
                                          25.68d0, 5.037d0, 1.765d0, 2.963d0 /)
   real(kind=8),dimension(1:16) :: y_w = (/ 6.385d-01, 7.354d-01, 4.109d+00, 1.863d+00, 7.880d-01, &
                                            0.000d+00, 1.503d+00, 1.153d-03, 1.492d+00, 2.310d+00, &
                                            1.678d+00, 7.033d-01, 2.738d+01, 0.000d+00, 0.000d+00, &
                                            0.000d+00 /)
   real(kind=8),dimension(1:16) :: y_0 = (/ 9.935d-01, 2.782d+00, 1.568d+01, 1.975d+01, 9.424d+01, &
                                            0.000d+00, 2.222d+02, 1.869d+01, 1.615d-03, 1.775d+01, &
                                            3.426d+01, 2.279d+02, 2.203d-02, 0.000d+00, 0.000d+00, &
                                            0.000d+00 /)
   real(kind=8),dimension(1:16) :: y_1 = (/ 0.2486d0, 0.1788d0, 9.421d0, 3.361d0, 0.6265d0, 0.d0, &
                                            4.606d0, 0.3037d0, 0.4049d0, 1.663d0, 0.137d0, 1.172d0, &
                                            0.01073d0, 0.d0, 0.d0, 0.d0 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_sulfur = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_sulfur = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_sulfur = 0.d0
   endif

END FUNCTION getCrosssection_sulfur

FUNCTION getCrosssection_iron(lambda, N)
   ! Returns the cross section in cm^2 for iron
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_iron
   integer           :: N
   real(kind=8),dimension(1:26) :: E_th = (/ 7.902d+00, 1.619d+01, 3.065d+01, 5.480d+01, &
                                             7.501d+01, 9.906d+01, 1.250d+02, 1.511d+02, &
                                             2.336d+02, 2.621d+02, 2.902d+02, 3.308d+02, &
                                             3.610d+02, 3.922d+02, 4.570d+02, 4.893d+02, &
                                             1.262d+03, 1.358d+03, 1.456d+03, 1.582d+03, &
                                             1.689d+03, 1.799d+03, 1.950d+03, 2.046d+03, &
                                             8.829d+03, 9.278d+03 /)
   real(kind=8),dimension(1:26) :: E_max = (/ 66.d0, 76.17d0, 87.05d0, 106.7d0, 128.8d0,  &
                                              152.7d0, 178.3d0, 205.5d0, 921.1d0, 959.d0, &
                                              998.3d0, 1039.d0, 1081.d0, 1125.d0, 1181.d0, &
                                              1216.d0, 7651.d0, 7769.d0, 7918.d0, 8041.d0, &
                                              8184.d0, 8350.d0, 8484.d0, 8638.d0, 50000.d0, &
                                              50000.d0 /)
   real(kind=8),dimension(1:26) :: E_0 = (/ 5.461d-02, 1.761d-01, 1.698d-01, 2.544d+01, &
                                            7.256d-01, 2.656d+00, 5.059d+00, 7.098d-02, &
                                            6.741d+00, 6.886d+01, 8.284d+00, 6.295d+00, &
                                            1.317d-01, 8.509d-01, 5.555d-02, 2.873d+01, &
                                            3.444d-01, 3.190d+01, 7.519d-04, 2.011d+01, &
                                            9.243d+00, 9.713d+00, 4.575d+01, 7.326d+01, &
                                            1.057d+03, 2.932d+02 /)
   real(kind=8),dimension(1:26) :: sig_0 = (/ 3.062d-01, 4.365d+03, 6.107d+00, 3.653d+02, &
                                              1.523d-03, 5.259d-01, 2.420d+04, 1.979d+01, &
                                              2.687d+01, 6.470d+01, 3.281d+00, 1.738d+00, &
                                              2.791d-03, 1.454d-01, 2.108d+02, 1.207d+01, &
                                              1.452d+00, 2.388d+00, 6.066d-05, 4.455d-01, &
                                              1.098d+01, 7.204d-02, 2.580d+04, 1.276d+01, &
                                              1.195d+01, 8.099d+01 /)
   real(kind=8),dimension(1:26) :: y_a = (/ 2.671d+07, 6.298d+03, 1.555d+03, 8.913d+00, &
                                            3.736d+01, 1.450d+01, 4.850d+04, 1.745d+04, &
                                            1.807d+02, 2.062d+01, 5.360d+01, 1.130d+02, &
                                            2.487d+03, 1.239d+03, 2.045d+04, 5.150d+02, &
                                            3.960d+02, 2.186d+01, 1.606d+06, 4.236d+01, &
                                            7.637d+01, 1.853d+02, 1.358d+00, 4.914d+01, &
                                            5.769d+01, 3.288d+01 /)
   real(kind=8),dimension(1:26) :: P = (/ 7.923d0, 5.204d0, 8.055d0, 6.538d0, 17.67d0, &
                                          16.32d0, 2.374d0, 6.75d0, 6.29d0, 4.111d0, &
                                          8.571d0, 8.037d0, 9.791d0, 8.066d0, 6.033d0, &
                                          3.846d0, 10.13d0, 9.589d0, 8.813d0, 9.724d0, &
                                          7.962d0, 8.843d0, 26.04d0, 4.941d0, 1.718d0, &
                                          2.963d0 /)
   real(kind=8),dimension(1:26) :: y_w = (/ 2.069d+01, 1.141d+01, 8.698d+00, 5.602d-01, &
                                            5.064d+01, 1.558d+01, 2.516d-03, 2.158d+02, &
                                            2.387d-04, 2.778d-04, 3.279d-01, 3.096d-01, &
                                            6.938d-01, 4.937d-01, 1.885d-03, 0.000d+00, &
                                            1.264d+00, 2.902d-02, 4.398d+00, 2.757d+00, &
                                            1.748d+00, 9.551d-03, 2.723d+01, 0.000d+00, &
                                            0.000d+00, 0.000d+00 /)
   real(kind=8),dimension(1:26) :: y_0 = (/ 1.382d+02, 9.272d+01, 1.760d+02, 0.000d+00, &
                                            8.871d+01, 3.361d+01, 4.546d-01, 2.542d+03, &
                                            2.494d+01, 1.190d-05, 2.971d+01, 4.671d+01, &
                                            2.170d+03, 4.505d+02, 2.706d-04, 0.000d+00, &
                                            2.891d+01, 3.805d+01, 1.915d+06, 6.847d+01, &
                                            4.446d+01, 1.702d+02, 3.582d-02, 0.000d+00, &
                                            0.000d+00, 0.000d+00 /)
   real(kind=8),dimension(1:26) :: y_1 = (/ 2.481d-01, 1.075d+02, 1.847d+01, 0.000d+00, &
                                            5.280d-02, 3.743d-03, 2.683d+01, 4.672d+02, &
                                            8.251d+00, 6.570d-03, 5.220d-01, 1.425d-01, &
                                            6.852d-03, 2.504d+00, 1.628d+00, 0.000d+00, &
                                            3.404d+00, 4.805d-01, 3.140d+01, 3.989d+00, &
                                            3.512d+00, 4.263d+00, 8.712d-03, 0.000d+00, &
                                            0.000d+00, 0.000d+00 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_iron = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_iron = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_iron = 0.d0
   endif

END FUNCTION getCrosssection_iron

FUNCTION getCrosssection_neon(lambda, N)
   ! Returns the cross section in cm^2 for neon
   ! from http://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
   use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
   real(kind=8)      :: lambda, getCrosssection_neon
   integer           :: N
   real(kind=8),dimension(1:10) :: E_th = (/ 2.156E+01, 4.096E+01, 6.346E+01, 9.712E+01, 1.262E+02, &
                                             1.579E+02, 2.073E+02, 2.391E+02, 1.196E+03, 1.362E+03 /)
   real(kind=8),dimension(1:10) :: E_max = (/ 8.701E+02, 8.831E+02, 9.131E+02, 9.480E+02, 9.873E+02, &
                                              1.031E+03, 1.078E+03, 1.125E+03, 5.000E+04, 5.000E+04 /)
   real(kind=8),dimension(1:10) :: E_0 = (/ 4.870E+00, 1.247E+01, 7.753E-01, 5.566E+00, 1.248E+00, &
                                            1.499E+00, 4.888E+00, 1.003E+01, 1.586E+02, 4.304E+01 /)
   real(kind=8),dimension(1:10) :: sig_0 = (/ 4.287E+03, 1.583E+03, 5.708E+00, 1.685E+03, 2.430E+00, &
                                              9.854E-01, 1.198E+04, 5.631E+01, 6.695E+01, 5.475E+02 /)
   real(kind=8),dimension(1:10) :: y_a = (/ 5.798E+00, 3.935E+00, 6.725E+01, 6.409E+02, 1.066E+02, &
                                            1.350E+02, 1.788E+00, 3.628E+01, 3.352E+01, 3.288E+01 /)
   real(kind=8),dimension(1:10) :: P = (/ 8.355E+00, 7.810E+00, 1.005E+01, 3.056E+00, 8.999E+00, &
                                          8.836E+00, 2.550E+01, 5.585E+00, 2.002E+00, 2.963E+00 /)
   real(kind=8),dimension(1:10) :: y_w = (/ 2.434E-01, 6.558E-02, 4.633E-01, 8.290E-03, 6.855E-01, &
                                            1.656E+00, 2.811E+01, 0.000E+00, 0.000E+00, 0.000E+00 /)
   real(kind=8),dimension(1:10) :: y_0 = (/ 4.236E-02, 1.520E+00, 7.654E+01, 5.149E+00, 9.169E+01, &
                                            1.042E+02, 2.536E-02, 0.000E+00, 0.000E+00, 0.000E+00 /)
   real(kind=8),dimension(1:10) :: y_1 = (/ 5.873E+00, 1.084E-01, 2.023E+00, 6.687E+00, 3.702E-01, &
                                            1.435E+00, 4.417E-02, 0.000E+00, 0.000E+00, 0.000E+00 /)
   real(kind=8) :: x, y, F, E

   ! Convert lambda into ev
   E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev

   x = (E / E_0(N)) - y_0(N)
   y = sqrt( (x*x) + (y_1(N)*y_1(N)) )

   F = (x - 1.d0)**2
   F = F + (y_w(N)*y_w(N))
   F = F * (y**(0.5d0*P(N) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/y_a(N)))**(-1.d0*P(N)))

   getCrosssection_neon = sig_0(N) * F * 1.d-18
   if (E.lt.E_th(N)) then
      getCrosssection_neon = 0.d0
   endif

   if (E.gt.E_max(N)) then
      getCrosssection_neon = 0.d0
   endif

END FUNCTION getCrosssection_neon

FUNCTION getCrosssection_BARE_GR_S_DUST(lambda,species)
   ! Harley's fit to the effective absorption cross
   ! section of dust for the BARE-GR-S model
   ! Mean absolute percentage error of this fit is
   ! 20%, max error is 133%.
   !
   ! The fit is a degree 10 polynomial
   ! lambda is assumed to be in angstroms
   ! Returns the absorption cross section in cm^-2 / H
   implicit none
   real(kind=8)      :: lambda, getCrosssection_BARE_GR_S_DUST
   real(kind=8)      :: lambda_microns
   integer           :: species,i
   real(kind=8),dimension(1:11) :: fit_vals = (/ -1.59319023e+01, -1.60473171e+00,  6.20612550e-01, &
                                                  6.42859480e-01, -4.08743189e-01, -1.59224607e-01, &
                                                  7.37953364e-02,  1.60696953e-02, -5.96977205e-03, &
                                                 -5.57671237e-04,  1.80437634e-04 /)

   lambda_microns = lambda * 1.d-4 ! Convert from angstroms to microns
   getCrosssection_BARE_GR_S_DUST =0.d0
   do i=1,11
      getCrosssection_BARE_GR_S_DUST = getCrosssection_BARE_GR_S_DUST + (fit_vals(i) * (LOG10(lambda_microns)**REAL(i-1, kind=8)))
   end do
   getCrosssection_BARE_GR_S_DUST = 10.d0**getCrosssection_BARE_GR_S_DUST

END FUNCTION getCrosssection_BARE_GR_S_DUST

END MODULE spectrum_integrator_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for Stellar Energy Distribution table.
!_________________________________________________________________________
!
MODULE SED_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  use rt_parameters,only:nGroups
  implicit none

  PUBLIC nSEDgroups                                                      &
      , init_SED_table, inp_SED_table, update_SED_group_props            &
      , update_star_RT_feedback, star_RT_feedback

  PRIVATE   ! default

  ! Light properties for different spectral energy distributions----------
  integer::nSEDgroups=nGroups         ! Default: All groups are SED groups
  integer::SED_nA, SED_nZ=8           ! Number of age bins and Z bins
  ! Age and z logarithmic intervals and lowest values:
  real(dp),parameter::SED_dlgA=0.02d0
  real(dp)::SED_dlgZ
  real(dp)::SED_lgA0, SED_lgZ0
  real(dp),allocatable,dimension(:)::SED_ages,SED_zeds![Gyr],[m_met/m_gas]
  ! SED_table: iAges, imetallicities, igroups, properties
  !                                         (Lum, Lum-acc, egy, csn, cse).
  ! Lum is photons per sec per solar mass (eV per sec per solar mass in
  ! the case of SED_isEgy=true). Lum-acc is accumulated lum.
  real(dp),allocatable,dimension(:,:,:,:)::SED_table
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_oxygen
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_nitrogen
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_carbon
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_magnesium
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_silicon
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_sulfur
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_iron
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_neon
  real(dp),allocatable,dimension(:,:,:,:)::SED_table_dust

  ! ----------------------------------------------------------------------

CONTAINS

!*************************************************************************
SUBROUTINE init_SED_table()

! Initiate SED properties table, which gives photon luminosities,
! integrated luminosities, average photon cross sections and energies of
! each photon group as a function of stellar population age and
! metallicity.  The SED is read from a directory specified by sed_dir.
!-------------------------------------------------------------------------
  use amr_commons,only:myid,IOGROUPSIZE,ncpu
  use rt_parameters
  use spectrum_integrator_module
  use mpi_mod
  ! Temporary SSP/SED parameters (read from SED files):
  integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
  real(kind=8),allocatable::ages(:), Zs(:), Ls(:), rebAges(:)
  real(kind=8),allocatable::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
  real(kind=8),allocatable::tbl(:,:,:), tbl2(:,:,:), reb_tbl(:,:,:)
  real(kind=8),allocatable::tbl_oxygen(:,:,:),tbl_nitrogen(:,:,:)
  real(kind=8),allocatable::tbl_carbon(:,:,:),tbl_magnesium(:,:,:),tbl_silicon(:,:,:)
  real(kind=8),allocatable::tbl_sulfur(:,:,:),tbl_iron(:,:,:),tbl_neon(:,:,:)
  real(kind=8),allocatable::tbl_oxygen2(:,:,:),tbl_nitrogen2(:,:,:)
  real(kind=8),allocatable::tbl_carbon2(:,:,:),tbl_magnesium2(:,:,:),tbl_silicon2(:,:,:)
  real(kind=8),allocatable::tbl_sulfur2(:,:,:),tbl_iron2(:,:,:),tbl_neon2(:,:,:)
  real(kind=8),allocatable::reb_tbl_oxygen(:,:,:),reb_tbl_nitrogen(:,:,:)
  real(kind=8),allocatable::reb_tbl_carbon(:,:,:),reb_tbl_magnesium(:,:,:),reb_tbl_silicon(:,:,:)
  real(kind=8),allocatable::reb_tbl_sulfur(:,:,:),reb_tbl_iron(:,:,:),reb_tbl_neon(:,:,:)
  real(kind=8),allocatable::tbl_dust(:,:,:), tbl2_dust(:,:,:), reb_tbl_dust(:,:,:)
  integer::i,ia,iz,ip,ii,il,dum
  character(len=128)::fZs, fAges, fSEDs                        ! Filenames
  logical::ok,okAge,okZ
  real(kind=8)::dlgA, pL0, pL1, tmp
  integer::locid,ncpu2,ierr
  integer::nv=3+2*nIons  ! # vars in SED table: L,Lacc,egy,nions*(csn,egy)
  integer,parameter::tag=1132
  integer::dummy_io,info2
  character(len=128)::sedwrite_dir,sedwrite_cmd
  integer::info

!-------------------------------------------------------------------------
  if(myid==1) &
        write(*,*) 'Stars are photon emitting, so initializing SED table'
  if(SED_DIR.eq.'')call get_environment_variable('RAMSES_SED_DIR',SED_DIR)
  inquire(FILE=TRIM(sed_dir)//'/all_seds.dat', exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then
        write(*,*)'Cannot access SED directory ',TRIM(sed_dir)
        write(*,*)'Directory '//TRIM(sed_dir)//' not found'
        write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
                  ' to the correct path, or use the namelist.'
     endif
     call clean_stop
  end if
  write(fZs,'(a,a)')   trim(sed_dir),"/metallicity_bins.dat"
  write(fAges,'(a,a)') trim(sed_dir),"/age_bins.dat"
  write(fSEDs,'(a,a)') trim(sed_dir),"/all_seds.dat"
  inquire(file=fZs, exist=okZ)
  inquire(file=fAges, exist=okAge)
  inquire(file=fSEDs, exist=ok)
  if(.not. ok .or. .not. okAge .or. .not. okZ) then
     if(myid.eq.1) then
        write(*,*) 'Cannot read SED files...'
        write(*,*) 'Check if SED-directory contains the files ',  &
                   'metallicity_bins.dat, age_bins.dat, all_seds.dat'
     endif
     call clean_stop
  end if

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  ! READ METALLICITY BINS-------------------------------------------------
  open(unit=10,file=fZs,status='old',form='formatted')
  read(10,'(i8)') nzs
  allocate(zs(nzs))
  do i = 1, nzs
     read(10,'(e14.6)') zs(i)
  end do
  close(10)
  ! READ AGE BINS---------------------------------------------------------
  open(unit=10,file=fAges,status='old',form='formatted')
  read(10,'(i8)') nAges
  allocate(ages(nAges))
  do i = 1, nAges
     read(10,'(e14.6)') ages(i)
  end do
  close(10)
  ages = ages*1.e-9                       !         Convert from yr to Gyr
  if(ages(1) .ne. 0.) ages(1) = 0.
  ! READ SEDS-------------------------------------------------------------
  if(KatzSED)then
     open(unit=10,file=fSEDs,status='old',form='formatted')
     read(10,'(i8)') nLs
     allocate(Ls(nLs))
     do i = 1, nLs
        read(10,'(e14.6)') Ls(i)
     end do
     allocate(SEDs(nLs,nAges,nzs))
     do iz = 1, nzs
        do ia = 1, nAges
           do il = 1, nLs
              read(10,'(e14.6)') SEDs(il,ia,iz)
           end do
        end do
     end do
     close(10)
  else
     open(unit=10,file=fSEDs,status='old',form='unformatted')
     read(10) nLs, dum
     allocate(Ls(nLs))
     read(10) Ls(:)
     allocate(SEDs(nLs,nAges,nzs))
     do iz = 1, nzs
        do ia = 1, nAges
           read(10) SEDs(:,ia,iz)
        end do
     end do
     close(10)
  endif

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

  ! If MPI then share the SED integration between the cpus:
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,locid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu2,ierr)
#endif
#ifdef WITHOUTMPI
  locid=0
  ncpu2=1
#endif

  ! Perform SED integration of luminosity, csn and egy per (age,Z) bin----
  allocate(tbl(nAges,nZs,nv))

  if (oxygen_ions) then
   allocate(tbl_oxygen(nAges,nZs,2*n_oxygen_ions))
  endif
  if (nitrogen_ions) then
   allocate(tbl_nitrogen(nAges,nZs,2*n_nitrogen_ions))
  endif
  if (carbon_ions) then
   allocate(tbl_carbon(nAges,nZs,2*n_carbon_ions))
  endif
  if (magnesium_ions) then
   allocate(tbl_magnesium(nAges,nZs,2*n_magnesium_ions))
  endif
  if (silicon_ions) then
   allocate(tbl_silicon(nAges,nZs,2*n_silicon_ions))
  endif
  if (sulfur_ions) then
   allocate(tbl_sulfur(nAges,nZs,2*n_sulfur_ions))
  endif
  if (iron_ions) then
   allocate(tbl_iron(nAges,nZs,2*n_iron_ions))
  endif
  if (neon_ions) then
   allocate(tbl_neon(nAges,nZs,2*n_neon_ions))
  endif

  ! Dust table
  allocate(tbl_dust(nAges,nZs,2))


  do ip = 1,nSEDgroups                                ! Loop photon groups
     tbl=0.
     if (oxygen_ions)    tbl_oxygen=0.
     if (nitrogen_ions)  tbl_nitrogen=0.
     if (carbon_ions)    tbl_carbon=0.
     if (magnesium_ions) tbl_magnesium=0.
     if (silicon_ions)   tbl_silicon=0.
     if (sulfur_ions)    tbl_sulfur=0.
     if (iron_ions)      tbl_iron=0.
     if (neon_ions)      tbl_neon=0.
     tbl_dust=0.

     pL0 = groupL0(ip) ; pL1 = groupL1(ip)! eV interval of photon group ip
     do iz = 1, nzs                                     ! Loop metallicity
     do ia = locid+1,nAges,ncpu2                                 ! Loop age
        tbl(ia,iz,1) = getSEDLuminosity(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
        tbl(ia,iz,3) = getSEDEgy(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
        do ii = 1,nIons                                     ! Loop species
           tbl(ia,iz,2+ii*2) = getSEDcsn(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
           tbl(ia,iz,3+ii*2) = getSEDcse(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
        end do ! End species loop

         ! LOOP OVER METALS
         if (oxygen_ions) then
            do ii=0,n_oxygen_ions-1 ! Loop over ionisation states
               tbl_oxygen(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,1,ii+1)
               tbl_oxygen(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,1,ii+1)
            enddo
         endif
         if (nitrogen_ions) then
            do ii=0,n_nitrogen_ions-1
               tbl_nitrogen(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,2,ii+1)
               tbl_nitrogen(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,2,ii+1)
            enddo
         endif
         if (carbon_ions) then
            do ii=0,n_carbon_ions-1
               tbl_carbon(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,3,ii+1)
               tbl_carbon(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,3,ii+1)
            enddo
         endif
         if (magnesium_ions) then
            do ii=0,n_magnesium_ions-1
               tbl_magnesium(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,4,ii+1)
               tbl_magnesium(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,4,ii+1)
            enddo
         endif
         if (silicon_ions) then
            do ii=0,n_silicon_ions-1
               tbl_silicon(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,5,ii+1)
               tbl_silicon(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,5,ii+1)
            enddo
         endif
         if (sulfur_ions) then
            do ii=0,n_sulfur_ions-1
               tbl_sulfur(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,6,ii+1)
               tbl_sulfur(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,6,ii+1)
            enddo
         endif
         if (iron_ions) then
            do ii=0,n_iron_ions-1
               tbl_iron(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,7,ii+1)
               tbl_iron(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,7,ii+1)
            enddo
         endif
         if (neon_ions) then
            do ii=0,n_neon_ions-1
               tbl_neon(ia,iz,ii*2+1)   = getSEDcsn_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,8,ii+1)
               tbl_neon(ia,iz,1+ii*2+1) = getSEDcse_metals(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,8,ii+1)
            enddo
         endif

         ! Now deal with dust
         tbl_dust(ia,iz,1)   = getSEDcsn_dust(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,0)

     end do ! End age loop
     end do ! End Z loop

#ifndef WITHOUTMPI
     allocate(tbl2(nAges,nzs,nv))
     call MPI_ALLREDUCE(tbl,tbl2,nAges*nzs*nv,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     tbl = tbl2
     deallocate(tbl2)

     if (oxygen_ions) then
      allocate(tbl_oxygen2(nAges,nzs,2*n_oxygen_ions))
      call MPI_ALLREDUCE(tbl_oxygen,tbl_oxygen2,nAges*nzs*2*n_oxygen_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_oxygen = tbl_oxygen2
      deallocate(tbl_oxygen2)
     endif

     if (nitrogen_ions) then
      allocate(tbl_nitrogen2(nAges,nzs,2*n_nitrogen_ions))
      call MPI_ALLREDUCE(tbl_nitrogen,tbl_nitrogen2,nAges*nzs*2*n_nitrogen_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_nitrogen = tbl_nitrogen2
      deallocate(tbl_nitrogen2)
     endif

     if (carbon_ions) then
      allocate(tbl_carbon2(nAges,nzs,2*n_carbon_ions))
      call MPI_ALLREDUCE(tbl_carbon,tbl_carbon2,nAges*nzs*2*n_carbon_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_carbon = tbl_carbon2
      deallocate(tbl_carbon2)
     endif

     if (magnesium_ions) then
      allocate(tbl_magnesium2(nAges,nzs,2*n_magnesium_ions))
      call MPI_ALLREDUCE(tbl_magnesium,tbl_magnesium2,nAges*nzs*2*n_magnesium_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_magnesium = tbl_magnesium2
      deallocate(tbl_magnesium2)
     endif

     if (silicon_ions) then
      allocate(tbl_silicon2(nAges,nzs,2*n_silicon_ions))
      call MPI_ALLREDUCE(tbl_silicon,tbl_silicon2,nAges*nzs*2*n_silicon_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_silicon = tbl_silicon2
      deallocate(tbl_silicon2)
     endif

     if (sulfur_ions) then
      allocate(tbl_sulfur2(nAges,nzs,2*n_sulfur_ions))
      call MPI_ALLREDUCE(tbl_sulfur,tbl_sulfur2,nAges*nzs*2*n_sulfur_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_sulfur = tbl_sulfur2
      deallocate(tbl_sulfur2)
     endif

     if (iron_ions) then
      allocate(tbl_iron2(nAges,nzs,2*n_iron_ions))
      call MPI_ALLREDUCE(tbl_iron,tbl_iron2,nAges*nzs*2*n_iron_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_iron = tbl_iron2
      deallocate(tbl_iron2)
     endif

     if (neon_ions) then
      allocate(tbl_neon2(nAges,nzs,2*n_neon_ions))
      call MPI_ALLREDUCE(tbl_neon,tbl_neon2,nAges*nzs*2*n_neon_ions,&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      tbl_neon = tbl_neon2
      deallocate(tbl_neon2)
     endif

     allocate(tbl2_dust(nAges,nzs,2))
     call MPI_ALLREDUCE(tbl_dust,tbl2_dust,nAges*nzs*2,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     tbl_dust = tbl2_dust
     deallocate(tbl2_dust)
#endif

     ! Integrate the cumulative luminosities:
     do iz = 1, nzs ! Loop metallicity
        tmp = trapz1( ages, tbl(:,iz,1), nAges, tbl(:,iz,2) )
     end do
     tbl(:,:,2) = tbl(:,:,2) * Gyr2sec ! Convert from sec-1 to Myr-1

     ! Now the SED properties are in tbl...just need to rebin it, to get !
     ! even log-intervals between bins for fast interpolation.           !
     dlgA = SED_dlgA ; SED_dlgZ = -SED_nz
     call rebin_log(dlgA, SED_dlgZ                                 &
          , tbl(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, nv        &
          , reb_tbl, SED_nA, SED_nZ, rebAges, SED_Zeds)

      if (oxygen_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_oxygen(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_oxygen_ions   &
         , reb_tbl_oxygen, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (nitrogen_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_nitrogen(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_nitrogen_ions   &
         , reb_tbl_nitrogen, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (carbon_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_carbon(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_carbon_ions   &
         , reb_tbl_carbon, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (magnesium_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_magnesium(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_magnesium_ions   &
         , reb_tbl_magnesium, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (silicon_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_silicon(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_silicon_ions   &
         , reb_tbl_silicon, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (sulfur_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_sulfur(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_sulfur_ions   &
         , reb_tbl_sulfur, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (iron_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_iron(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_iron_ions   &
         , reb_tbl_iron, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      if (neon_ions) then
         call rebin_log(dlgA, SED_dlgZ                                  &
         , tbl_neon(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2*n_neon_ions   &
         , reb_tbl_neon, SED_nA, SED_nZ, rebAges, SED_Zeds)
      endif

      call rebin_log(dlgA, SED_dlgZ                                  &
      , tbl_dust(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, 2   &
      , reb_tbl_dust, SED_nA, SED_nZ, rebAges, SED_Zeds)

     SED_nA=SED_nA+1                              ! Make room for zero age
     if(ip .eq. 1 ) then
      allocate(SED_table(SED_nA, SED_nZ, nSEDgroups, nv))

      if (oxygen_ions)    allocate(SED_table_oxygen(SED_nA, SED_nZ, nSEDgroups, 2*n_oxygen_ions))
      if (nitrogen_ions)  allocate(SED_table_nitrogen(SED_nA, SED_nZ, nSEDgroups, 2*n_nitrogen_ions))
      if (carbon_ions)    allocate(SED_table_carbon(SED_nA, SED_nZ, nSEDgroups, 2*n_carbon_ions))
      if (magnesium_ions) allocate(SED_table_magnesium(SED_nA, SED_nZ, nSEDgroups, 2*n_magnesium_ions))
      if (silicon_ions)   allocate(SED_table_silicon(SED_nA, SED_nZ, nSEDgroups, 2*n_silicon_ions))
      if (sulfur_ions)    allocate(SED_table_sulfur(SED_nA, SED_nZ, nSEDgroups, 2*n_sulfur_ions))
      if (iron_ions)      allocate(SED_table_iron(SED_nA, SED_nZ, nSEDgroups, 2*n_iron_ions))
      if (neon_ions)      allocate(SED_table_neon(SED_nA, SED_nZ, nSEDgroups, 2*n_neon_ions))

      allocate(SED_table_dust(SED_nA, SED_nZ, nSEDgroups, 2))

     endif

     SED_table(1, :,ip,:) = reb_tbl(1,:,:)            ! Zero age properties
     SED_table(1, :,ip,2) = 0.                        !  Lacc=0 at zero age
     SED_table(2:,:,ip,:) = reb_tbl

     if (oxygen_ions) then
      SED_table_oxygen(1, :,ip,:) = reb_tbl_oxygen(1,:,:)          ! Zero age properties
      SED_table_oxygen(2:,:,ip,:) = reb_tbl_oxygen
     endif
     if (nitrogen_ions) then
      SED_table_nitrogen(1, :,ip,:) = reb_tbl_nitrogen(1,:,:)      ! Zero age properties
      SED_table_nitrogen(2:,:,ip,:) = reb_tbl_nitrogen
     endif
     if (carbon_ions) then
      SED_table_carbon(1, :,ip,:) = reb_tbl_carbon(1,:,:)          ! Zero age properties
      SED_table_carbon(2:,:,ip,:) = reb_tbl_carbon
     endif
     if (magnesium_ions) then
      SED_table_magnesium(1, :,ip,:) = reb_tbl_magnesium(1,:,:)    ! Zero age properties
      SED_table_magnesium(2:,:,ip,:) = reb_tbl_magnesium
     endif
     if (silicon_ions) then
      SED_table_silicon(1, :,ip,:) = reb_tbl_silicon(1,:,:)        ! Zero age properties
      SED_table_silicon(2:,:,ip,:) = reb_tbl_silicon
     endif
     if (sulfur_ions) then
      SED_table_sulfur(1, :,ip,:) = reb_tbl_sulfur(1,:,:)          ! Zero age properties
      SED_table_sulfur(2:,:,ip,:) = reb_tbl_sulfur
     endif
     if (iron_ions) then
      SED_table_iron(1, :,ip,:) = reb_tbl_iron(1,:,:)              ! Zero age properties
      SED_table_iron(2:,:,ip,:) = reb_tbl_iron
     endif
     if (neon_ions) then
      SED_table_neon(1, :,ip,:) = reb_tbl_neon(1,:,:)              ! Zero age properties
      SED_table_neon(2:,:,ip,:) = reb_tbl_neon
     endif

     SED_table_dust(1, :,ip,:) = reb_tbl_dust(1,:,:)                   ! Zero age properties
     SED_table_dust(2:,:,ip,:) = reb_tbl_dust

     deallocate(reb_tbl)

     if (oxygen_ions)    deallocate(reb_tbl_oxygen)
     if (nitrogen_ions)  deallocate(reb_tbl_nitrogen)
     if (carbon_ions)    deallocate(reb_tbl_carbon)
     if (magnesium_ions) deallocate(reb_tbl_magnesium)
     if (silicon_ions)   deallocate(reb_tbl_silicon)
     if (sulfur_ions)    deallocate(reb_tbl_sulfur)
     if (iron_ions)      deallocate(reb_tbl_iron)
     if (neon_ions)      deallocate(reb_tbl_neon)

     deallocate(reb_tbl_dust)

  end do ! End photon group loop

  SED_lgZ0 = log10(SED_Zeds(1))                  ! Interpolation intervals
  SED_lgA0 = log10(rebAges(1))
  allocate(SED_ages(SED_nA))
  SED_ages(1)=0.d0 ; SED_ages(2:)=rebAges ;    ! Must have zero initial age

  deallocate(SEDs) ; deallocate(tbl)
  deallocate(ages) ; deallocate(rebAges)
  deallocate(zs)
  deallocate(Ls)

  if (oxygen_ions)    deallocate(tbl_oxygen)
  if (nitrogen_ions)  deallocate(tbl_nitrogen)
  if (carbon_ions)    deallocate(tbl_carbon)
  if (magnesium_ions) deallocate(tbl_magnesium)
  if (silicon_ions)   deallocate(tbl_silicon)
  if (sulfur_ions)    deallocate(tbl_sulfur)
  if (iron_ions)      deallocate(tbl_iron)
  if (neon_ions)      deallocate(tbl_neon)

  deallocate(tbl_dust)

  sedwrite_dir = 'SEDtables/'
  sedwrite_cmd = 'mkdir -p '//trim(sedwrite_dir)
#ifdef NOSYSTEM
  if(myid==1)call PXFMKDIR(TRIM(sedwrite_dir),LEN(TRIM(sedwrite_dir)),int(O'755'),info)
#else
  if(myid==1)call system(sedwrite_cmd)
#endif

  if (myid==1) call write_SEDtable
  if (myid==1.and.oxygen_ions)    call write_SEDtable_oxygen
  if (myid==1.and.nitrogen_ions)  call write_SEDtable_nitrogen
  if (myid==1.and.carbon_ions)    call write_SEDtable_carbon
  if (myid==1.and.magnesium_ions) call write_SEDtable_magnesium
  if (myid==1.and.silicon_ions)   call write_SEDtable_silicon
  if (myid==1.and.sulfur_ions)    call write_SEDtable_sulfur
  if (myid==1.and.iron_ions)      call write_SEDtable_iron
  if (myid==1.and.neon_ions)      call write_SEDtable_neon
  if (myid==1) call write_SEDtable_dust

END SUBROUTINE init_SED_table

!*************************************************************************
SUBROUTINE update_SED_group_props()

! Compute and assign to all SED photon groups an average of the
! quantities group_csn and group_egy from all the star particles in the
! simulation, weighted by the luminosity of the particles.
! If there are no stars, we assign the SED properties of zero-age,
! zero-metallicity stellar populations.
!-------------------------------------------------------------------------
  use amr_commons
  use pm_parameters
  use pm_commons
  use rt_parameters
  use mpi_mod
#ifndef WITHOUTMPI
  integer :: info
#endif
  integer :: i, ip, ii
  real(dp),save,allocatable,dimension(:)::  L_star
  real(dp),save,allocatable,dimension(:,:)::csn_star, cse_star
  real(dp),save,allocatable,dimension(:)::  egy_star
  real(dp),save,allocatable,dimension(:)::  sum_L_cpu,sum_L_all
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu,sum_csn_all
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu,sum_cse_all
  real(dp),save,allocatable,dimension(:)::sum_egy_cpu,sum_egy_all
  real(dp):: mass, age, Z

  real(dp),save,allocatable,dimension(:,:)::csn_star_oxygen, cse_star_oxygen
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_oxygen,sum_csn_all_oxygen
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_oxygen,sum_cse_all_oxygen
  real(dp),save,allocatable,dimension(:,:)::csn_star_nitrogen, cse_star_nitrogen
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_nitrogen,sum_csn_all_nitrogen
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_nitrogen,sum_cse_all_nitrogen
  real(dp),save,allocatable,dimension(:,:)::csn_star_carbon, cse_star_carbon
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_carbon,sum_csn_all_carbon
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_carbon,sum_cse_all_carbon
  real(dp),save,allocatable,dimension(:,:)::csn_star_magnesium, cse_star_magnesium
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_magnesium,sum_csn_all_magnesium
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_magnesium,sum_cse_all_magnesium
  real(dp),save,allocatable,dimension(:,:)::csn_star_silicon, cse_star_silicon
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_silicon,sum_csn_all_silicon
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_silicon,sum_cse_all_silicon
  real(dp),save,allocatable,dimension(:,:)::csn_star_sulfur, cse_star_sulfur
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_sulfur,sum_csn_all_sulfur
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_sulfur,sum_cse_all_sulfur
  real(dp),save,allocatable,dimension(:,:)::csn_star_iron, cse_star_iron
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_iron,sum_csn_all_iron
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_iron,sum_cse_all_iron
  real(dp),save,allocatable,dimension(:,:)::csn_star_neon, cse_star_neon
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_neon,sum_csn_all_neon
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu_neon,sum_cse_all_neon

  real(dp),save,allocatable,dimension(:,:)::csn_star_dust
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu_dust,sum_csn_all_dust
!-------------------------------------------------------------------------

  if(.not. allocated(L_star)) then
     allocate(L_star(nSEDgroups))
     allocate(egy_star(nSEDgroups))
     allocate(csn_star(nSEDgroups,nIons))
     allocate(cse_star(nSEDgroups,nIons))
     allocate(sum_L_cpu(nSEDgroups))
     allocate(sum_L_all(nSEDgroups))
     allocate(sum_egy_cpu(nSEDgroups))
     allocate(sum_egy_all(nSEDgroups))
     allocate(sum_csn_cpu(nSEDgroups,nIons))
     allocate(sum_csn_all(nSEDgroups,nIons))
     allocate(sum_cse_cpu(nSEDgroups,nIons))
     allocate(sum_cse_all(nSEDgroups,nIons))

     if (oxygen_ions) then
      allocate(csn_star_oxygen(nSEDgroups,n_oxygen_ions))
      allocate(cse_star_oxygen(nSEDgroups,n_oxygen_ions))
      allocate(sum_csn_cpu_oxygen(nSEDgroups,n_oxygen_ions))
      allocate(sum_csn_all_oxygen(nSEDgroups,n_oxygen_ions))
      allocate(sum_cse_cpu_oxygen(nSEDgroups,n_oxygen_ions))
      allocate(sum_cse_all_oxygen(nSEDgroups,n_oxygen_ions))
     endif

     if (nitrogen_ions) then
      allocate(csn_star_nitrogen(nSEDgroups,n_nitrogen_ions))
      allocate(cse_star_nitrogen(nSEDgroups,n_nitrogen_ions))
      allocate(sum_csn_cpu_nitrogen(nSEDgroups,n_nitrogen_ions))
      allocate(sum_csn_all_nitrogen(nSEDgroups,n_nitrogen_ions))
      allocate(sum_cse_cpu_nitrogen(nSEDgroups,n_nitrogen_ions))
      allocate(sum_cse_all_nitrogen(nSEDgroups,n_nitrogen_ions))
     endif

     if (carbon_ions) then
      allocate(csn_star_carbon(nSEDgroups,n_carbon_ions))
      allocate(cse_star_carbon(nSEDgroups,n_carbon_ions))
      allocate(sum_csn_cpu_carbon(nSEDgroups,n_carbon_ions))
      allocate(sum_csn_all_carbon(nSEDgroups,n_carbon_ions))
      allocate(sum_cse_cpu_carbon(nSEDgroups,n_carbon_ions))
      allocate(sum_cse_all_carbon(nSEDgroups,n_carbon_ions))
     endif

     if (magnesium_ions) then
      allocate(csn_star_magnesium(nSEDgroups,n_magnesium_ions))
      allocate(cse_star_magnesium(nSEDgroups,n_magnesium_ions))
      allocate(sum_csn_cpu_magnesium(nSEDgroups,n_magnesium_ions))
      allocate(sum_csn_all_magnesium(nSEDgroups,n_magnesium_ions))
      allocate(sum_cse_cpu_magnesium(nSEDgroups,n_magnesium_ions))
      allocate(sum_cse_all_magnesium(nSEDgroups,n_magnesium_ions))
     endif

     if (silicon_ions) then
      allocate(csn_star_silicon(nSEDgroups,n_silicon_ions))
      allocate(cse_star_silicon(nSEDgroups,n_silicon_ions))
      allocate(sum_csn_cpu_silicon(nSEDgroups,n_silicon_ions))
      allocate(sum_csn_all_silicon(nSEDgroups,n_silicon_ions))
      allocate(sum_cse_cpu_silicon(nSEDgroups,n_silicon_ions))
      allocate(sum_cse_all_silicon(nSEDgroups,n_silicon_ions))
     endif

     if (sulfur_ions) then
      allocate(csn_star_sulfur(nSEDgroups,n_sulfur_ions))
      allocate(cse_star_sulfur(nSEDgroups,n_sulfur_ions))
      allocate(sum_csn_cpu_sulfur(nSEDgroups,n_sulfur_ions))
      allocate(sum_csn_all_sulfur(nSEDgroups,n_sulfur_ions))
      allocate(sum_cse_cpu_sulfur(nSEDgroups,n_sulfur_ions))
      allocate(sum_cse_all_sulfur(nSEDgroups,n_sulfur_ions))
     endif

     if (iron_ions) then
      allocate(csn_star_iron(nSEDgroups,n_iron_ions))
      allocate(cse_star_iron(nSEDgroups,n_iron_ions))
      allocate(sum_csn_cpu_iron(nSEDgroups,n_iron_ions))
      allocate(sum_csn_all_iron(nSEDgroups,n_iron_ions))
      allocate(sum_cse_cpu_iron(nSEDgroups,n_iron_ions))
      allocate(sum_cse_all_iron(nSEDgroups,n_iron_ions))
     endif

     if (neon_ions) then
      allocate(csn_star_neon(nSEDgroups,n_neon_ions))
      allocate(cse_star_neon(nSEDgroups,n_neon_ions))
      allocate(sum_csn_cpu_neon(nSEDgroups,n_neon_ions))
      allocate(sum_csn_all_neon(nSEDgroups,n_neon_ions))
      allocate(sum_cse_cpu_neon(nSEDgroups,n_neon_ions))
      allocate(sum_cse_all_neon(nSEDgroups,n_neon_ions))
     endif

     allocate(csn_star_dust(nSEDgroups,2))
     allocate(sum_csn_cpu_dust(nSEDgroups,2))
     allocate(sum_csn_all_dust(nSEDgroups,2))

  endif
  sum_L_cpu   = 0d0 ! Accumulated luminosity, avg cross sections and
  sum_egy_cpu = 0d0 ! photon energies for all stars belonging to
  sum_csn_cpu = 0d0 ! 'this' cpu
  sum_cse_cpu = 0d0

  if (oxygen_ions) then
   csn_star_oxygen = 0.d0
   cse_star_oxygen = 0.d0
   sum_csn_all_oxygen = 0.d0
   sum_cse_all_oxygen = 0.d0
   sum_csn_cpu_oxygen = 0.d0 ! 'this' cpu
   sum_cse_cpu_oxygen = 0.d0
  endif

  if (nitrogen_ions) then
   csn_star_nitrogen = 0.d0
   cse_star_nitrogen = 0.d0
   sum_csn_all_nitrogen = 0.d0
   sum_cse_all_nitrogen = 0.d0
   sum_csn_cpu_nitrogen = 0.d0 ! 'this' cpu
   sum_cse_cpu_nitrogen = 0.d0
  endif

  if (carbon_ions) then
   csn_star_carbon = 0.d0
   cse_star_carbon = 0.d0
   sum_csn_all_carbon = 0.d0
   sum_cse_all_carbon = 0.d0
   sum_csn_cpu_carbon = 0.d0 ! 'this' cpu
   sum_cse_cpu_carbon = 0.d0
  endif

  if (magnesium_ions) then
   csn_star_magnesium = 0.d0
   cse_star_magnesium = 0.d0
   sum_csn_all_magnesium = 0.d0
   sum_cse_all_magnesium = 0.d0
   sum_csn_cpu_magnesium = 0.d0 ! 'this' cpu
   sum_cse_cpu_magnesium = 0.d0
  endif

  if (silicon_ions) then
   csn_star_silicon = 0.d0
   cse_star_silicon = 0.d0
   sum_csn_all_silicon = 0.d0
   sum_cse_all_silicon = 0.d0
   sum_csn_cpu_silicon = 0.d0 ! 'this' cpu
   sum_cse_cpu_silicon = 0.d0
  endif

  if (sulfur_ions) then
   csn_star_sulfur = 0.d0
   cse_star_sulfur = 0.d0
   sum_csn_all_sulfur = 0.d0
   sum_cse_all_sulfur = 0.d0
   sum_csn_cpu_sulfur = 0.d0 ! 'this' cpu
   sum_cse_cpu_sulfur = 0.d0
  endif

  if (iron_ions) then
   csn_star_iron = 0.d0
   cse_star_iron = 0.d0
   sum_csn_all_iron = 0.d0
   sum_cse_all_iron = 0.d0
   sum_csn_cpu_iron = 0.d0 ! 'this' cpu
   sum_cse_cpu_iron = 0.d0
  endif

  if (neon_ions) then
   csn_star_neon = 0.d0
   cse_star_neon = 0.d0
   sum_csn_all_neon = 0.d0
   sum_cse_all_neon = 0.d0
   sum_csn_cpu_neon = 0.d0 ! 'this' cpu
   sum_cse_cpu_neon = 0.d0
  endif

  csn_star_dust = 0.d0
  sum_csn_all_dust = 0.d0
  sum_csn_cpu_dust = 0.d0 ! 'this' cpu

  do i=1,npartmax
     if(levelp(i) <= 0 .or. (.not. is_star(typep(i)))) cycle  ! not a star

     if(is_pre_SN(typep(i)))then
        mass = mp(i)
     else
        mass = mp(i)/(1d0-eta_sn)
     endif
     if(use_initial_mass) mass = mp0(i)

     call getAgeGyr(tp(i), age)                            !     age = [Gyrs]
     if(metal) then
        !Z = zp(i)                                          ! [m_metals/m_tot]
        Z = 2.09d0*zp(i,2)+1.06d0*zp(i,1) !---- Stellar metallicity, Asplund 2009 !metal enrichment
     else
        Z = z_ave*0.02                                     ! [m_metals/m_tot]
     endif

     if(pop3 .and. is_pop_III(typep(i)))then

        ! time averaged quantities based on Schaerer (02)
        call get_pop3_Nphoton(mass,L_star)

        ! based on 1e5K blackbody spectrum
        call get_pop3_props  (1d5,egy_star,csn_star,cse_star)

     else
        Z = max(Z, 10.d-5)                                 ! [m_metals/m_tot]
        call inp_SED_table(age, Z, 1, .false., L_star)     !  [# s-1 m_sun-1]
        call inp_SED_table(age, Z, 3, .true., egy_star(:)) !             [ev]
        do ii=1,nIons
           call inp_SED_table(age, Z, 2+2*ii, .true., csn_star(:,ii))! [cm^2]
           call inp_SED_table(age, Z, 3+2*ii, .true., cse_star(:,ii))! [cm^2]
        end do

        if (oxygen_ions) then
         do ii=1,n_oxygen_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 1, ii, csn_star_oxygen(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 1, ii, cse_star_oxygen(:,ii))! [cm^2]
         enddo
        endif

        if (nitrogen_ions) then
         do ii=1,n_nitrogen_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 2, ii, csn_star_nitrogen(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 2, ii, cse_star_nitrogen(:,ii))! [cm^2]
         enddo
        endif

        if (carbon_ions) then
         do ii=1,n_carbon_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 3, ii, csn_star_carbon(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 3, ii, cse_star_carbon(:,ii))! [cm^2]
         enddo
        endif

        if (magnesium_ions) then
         do ii=1,n_magnesium_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 4, ii, csn_star_magnesium(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 4, ii, cse_star_magnesium(:,ii))! [cm^2]
         enddo
        endif

        if (silicon_ions) then
         do ii=1,n_silicon_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 5, ii, csn_star_silicon(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 5, ii, cse_star_silicon(:,ii))! [cm^2]
         enddo
        endif

        if (sulfur_ions) then
         do ii=1,n_sulfur_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 6, ii, csn_star_sulfur(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 6, ii, cse_star_sulfur(:,ii))! [cm^2]
         enddo
        endif

        if (iron_ions) then
         do ii=1,n_iron_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 7, ii, csn_star_iron(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 7, ii, cse_star_iron(:,ii))! [cm^2]
         enddo
        endif

        if (neon_ions) then
         do ii=1,n_neon_ions
            call inp_SED_table_metals(age, Z,   2*ii-1, .false., 8, ii, csn_star_neon(:,ii))! [cm^2]
            call inp_SED_table_metals(age, Z, 1+2*ii-1, .true., 8, ii, cse_star_neon(:,ii))! [cm^2]
         enddo
        endif

        call inp_SED_table_dust(age, Z, 1, .false., csn_star_dust(:,1))! [cm^2]

     endif

     do ip=1,nSEDgroups
        L_star(ip) = L_star(ip) * mass             !       [# photons s-1]
        sum_L_cpu(ip)    =   sum_L_cpu(ip)   + L_star(ip)
        sum_egy_cpu(ip) =  sum_egy_cpu(ip)   + L_star(ip) * egy_star(ip)
        sum_csn_cpu(ip,:)= sum_csn_cpu(ip,:) + L_star(ip) * csn_star(ip,:)
        sum_cse_cpu(ip,:)= sum_cse_cpu(ip,:) + L_star(ip) * cse_star(ip,:)

        if (oxygen_ions) then
         sum_csn_cpu_oxygen(ip,:)= sum_csn_cpu_oxygen(ip,:) + L_star(ip) * csn_star_oxygen(ip,:)
         sum_cse_cpu_oxygen(ip,:)= sum_cse_cpu_oxygen(ip,:) + L_star(ip) * cse_star_oxygen(ip,:)
        endif

        if (nitrogen_ions) then
         sum_csn_cpu_nitrogen(ip,:)= sum_csn_cpu_nitrogen(ip,:) + L_star(ip) * csn_star_nitrogen(ip,:)
         sum_cse_cpu_nitrogen(ip,:)= sum_cse_cpu_nitrogen(ip,:) + L_star(ip) * cse_star_nitrogen(ip,:)
        endif

        if (carbon_ions) then
         sum_csn_cpu_carbon(ip,:)= sum_csn_cpu_carbon(ip,:) + L_star(ip) * csn_star_carbon(ip,:)
         sum_cse_cpu_carbon(ip,:)= sum_cse_cpu_carbon(ip,:) + L_star(ip) * cse_star_carbon(ip,:)
        endif

        if (magnesium_ions) then
         sum_csn_cpu_magnesium(ip,:)= sum_csn_cpu_magnesium(ip,:) + L_star(ip) * csn_star_magnesium(ip,:)
         sum_cse_cpu_magnesium(ip,:)= sum_cse_cpu_magnesium(ip,:) + L_star(ip) * cse_star_magnesium(ip,:)
        endif

        if (silicon_ions) then
         sum_csn_cpu_silicon(ip,:)= sum_csn_cpu_silicon(ip,:) + L_star(ip) * csn_star_silicon(ip,:)
         sum_cse_cpu_silicon(ip,:)= sum_cse_cpu_silicon(ip,:) + L_star(ip) * cse_star_silicon(ip,:)
        endif

        if (sulfur_ions) then
         sum_csn_cpu_sulfur(ip,:)= sum_csn_cpu_sulfur(ip,:) + L_star(ip) * csn_star_sulfur(ip,:)
         sum_cse_cpu_sulfur(ip,:)= sum_cse_cpu_sulfur(ip,:) + L_star(ip) * cse_star_sulfur(ip,:)
        endif

        if (iron_ions) then
         sum_csn_cpu_iron(ip,:)= sum_csn_cpu_iron(ip,:) + L_star(ip) * csn_star_iron(ip,:)
         sum_cse_cpu_iron(ip,:)= sum_cse_cpu_iron(ip,:) + L_star(ip) * cse_star_iron(ip,:)
        endif

        if (neon_ions) then
         sum_csn_cpu_neon(ip,:)= sum_csn_cpu_neon(ip,:) + L_star(ip) * csn_star_neon(ip,:)
         sum_cse_cpu_neon(ip,:)= sum_cse_cpu_neon(ip,:) + L_star(ip) * cse_star_neon(ip,:)
        endif

        sum_csn_cpu_dust(ip,:)= sum_csn_cpu_dust(ip,:) + L_star(ip) * csn_star_dust(ip,:)

     end do

  end do

  ! Sum up for all cpus
#ifdef WITHOUTMPI
  sum_L_all   = sum_L_cpu
  sum_egy_all = sum_egy_cpu
  sum_csn_all = sum_csn_cpu
  sum_cse_all = sum_cse_cpu

  if (oxygen_ions) then
   sum_csn_all_oxygen = sum_csn_cpu_oxygen
   sum_cse_all_oxygen = sum_cse_cpu_oxygen
  endif

  if (nitrogen_ions) then
   sum_csn_all_nitrogen = sum_csn_cpu_nitrogen
   sum_cse_all_nitrogen = sum_cse_cpu_nitrogen
  endif

  if (carbon_ions) then
   sum_csn_all_carbon = sum_csn_cpu_carbon
   sum_cse_all_carbon = sum_cse_cpu_carbon
  endif

  if (magnesium_ions) then
   sum_csn_all_magnesium = sum_csn_cpu_magnesium
   sum_cse_all_magnesium = sum_cse_cpu_magnesium
  endif

  if (silicon_ions) then
   sum_csn_all_silicon = sum_csn_cpu_silicon
   sum_cse_all_silicon = sum_cse_cpu_silicon
  endif

  if (sulfur_ions) then
   sum_csn_all_sulfur = sum_csn_cpu_sulfur
   sum_cse_all_sulfur = sum_cse_cpu_sulfur
  endif

  if (iron_ions) then
   sum_csn_all_iron = sum_csn_cpu_iron
   sum_cse_all_iron = sum_cse_cpu_iron
  endif

  if (neon_ions) then
   sum_csn_all_neon = sum_csn_cpu_neon
   sum_cse_all_neon = sum_cse_cpu_neon
  endif

  sum_csn_all_dust = sum_csn_cpu_dust

#else
  call MPI_ALLREDUCE(sum_L_cpu,   sum_L_all,   nSEDgroups,               &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_egy_cpu, sum_egy_all, nSEDgroups,               &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_csn_cpu, sum_csn_all, nSEDgroups*nIons,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_cse_cpu, sum_cse_all, nSEDgroups*nIons,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

  if (oxygen_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_oxygen, sum_csn_all_oxygen, nSEDgroups*n_oxygen_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_oxygen, sum_cse_all_oxygen, nSEDgroups*n_oxygen_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (nitrogen_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_nitrogen, sum_csn_all_nitrogen, nSEDgroups*n_nitrogen_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_nitrogen, sum_cse_all_nitrogen, nSEDgroups*n_nitrogen_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (carbon_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_carbon, sum_csn_all_carbon, nSEDgroups*n_carbon_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_carbon, sum_cse_all_carbon, nSEDgroups*n_carbon_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (magnesium_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_magnesium, sum_csn_all_magnesium, nSEDgroups*n_magnesium_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_magnesium, sum_cse_all_magnesium, nSEDgroups*n_magnesium_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (silicon_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_silicon, sum_csn_all_silicon, nSEDgroups*n_silicon_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_silicon, sum_cse_all_silicon, nSEDgroups*n_silicon_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (sulfur_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_sulfur, sum_csn_all_sulfur, nSEDgroups*n_sulfur_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_sulfur, sum_cse_all_sulfur, nSEDgroups*n_sulfur_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (iron_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_iron, sum_csn_all_iron, nSEDgroups*n_iron_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_iron, sum_cse_all_iron, nSEDgroups*n_iron_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  if (neon_ions) then
   call MPI_ALLREDUCE(sum_csn_cpu_neon, sum_csn_all_neon, nSEDgroups*n_neon_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
   call MPI_ALLREDUCE(sum_cse_cpu_neon, sum_cse_all_neon, nSEDgroups*n_neon_ions,   &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  endif

  call MPI_ALLREDUCE(sum_csn_cpu_dust, sum_csn_all_dust, nSEDgroups*2,   &
  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

#endif

  ! ...and take averages weighted by luminosities
  do ip=1,nSEDgroups
     ! No update for non-SED groups (L0>L1):
     if(groupL0(ip).ne. 0d0 .and. groupL1(ip) .ne. 0d0 .and. &
          (groupL0(ip) .ge. groupL1(ip)) ) cycle
     if(sum_L_all(ip) .gt. 0.) then
        group_egy(ip)   = sum_egy_all(ip)   / sum_L_all(ip)
        group_csn(ip,:) = sum_csn_all(ip,:) / sum_L_all(ip)
        group_cse(ip,:) = sum_cse_all(ip,:) / sum_L_all(ip)

        if (oxygen_ions) then
         group_csn_oxygen(ip,:) = sum_csn_all_oxygen(ip,:) / sum_L_all(ip)
         group_cse_oxygen(ip,:) = sum_cse_all_oxygen(ip,:) / sum_L_all(ip)
        endif

        if (nitrogen_ions) then
         group_csn_nitrogen(ip,:) = sum_csn_all_nitrogen(ip,:) / sum_L_all(ip)
         group_cse_nitrogen(ip,:) = sum_cse_all_nitrogen(ip,:) / sum_L_all(ip)
        endif

        if (carbon_ions) then
         group_csn_carbon(ip,:) = sum_csn_all_carbon(ip,:) / sum_L_all(ip)
         group_cse_carbon(ip,:) = sum_cse_all_carbon(ip,:) / sum_L_all(ip)
        endif

        if (magnesium_ions) then
         group_csn_magnesium(ip,:) = sum_csn_all_magnesium(ip,:) / sum_L_all(ip)
         group_cse_magnesium(ip,:) = sum_cse_all_magnesium(ip,:) / sum_L_all(ip)
        endif

        if (silicon_ions) then
         group_csn_silicon(ip,:) = sum_csn_all_silicon(ip,:) / sum_L_all(ip)
         group_cse_silicon(ip,:) = sum_cse_all_silicon(ip,:) / sum_L_all(ip)
        endif

        if (sulfur_ions) then
         group_csn_sulfur(ip,:) = sum_csn_all_sulfur(ip,:) / sum_L_all(ip)
         group_cse_sulfur(ip,:) = sum_cse_all_sulfur(ip,:) / sum_L_all(ip)
        endif

        if (iron_ions) then
         group_csn_iron(ip,:) = sum_csn_all_iron(ip,:) / sum_L_all(ip)
         group_cse_iron(ip,:) = sum_cse_all_iron(ip,:) / sum_L_all(ip)
        endif

        if (neon_ions) then
         group_csn_neon(ip,:) = sum_csn_all_neon(ip,:) / sum_L_all(ip)
         group_cse_neon(ip,:) = sum_cse_all_neon(ip,:) / sum_L_all(ip)
        endif

        group_csn_dust(ip,:) = sum_csn_all_dust(ip,:) / sum_L_all(ip)

     else ! no stars -> assign zero-age zero-metallicity props
        group_egy(ip)       = SED_table(1,1,ip,3)
        do ii=1,nIons
           group_csn(ip,ii) = SED_table(1,1,ip,2+2*ii)
           group_cse(ip,ii) = SED_table(1,1,ip,3+2*ii)
        enddo

        if (oxygen_ions) then
         do ii=1,n_oxygen_ions
            group_csn_oxygen(ip,ii) = SED_table_oxygen(1,1,ip,  2*ii-1)
            group_cse_oxygen(ip,ii) = SED_table_oxygen(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (nitrogen_ions) then
         do ii=1,n_nitrogen_ions
            group_csn_nitrogen(ip,ii) = SED_table_nitrogen(1,1,ip,  2*ii-1)
            group_cse_nitrogen(ip,ii) = SED_table_nitrogen(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (carbon_ions) then
         do ii=1,n_carbon_ions
            group_csn_carbon(ip,ii) = SED_table_carbon(1,1,ip,  2*ii-1)
            group_cse_carbon(ip,ii) = SED_table_carbon(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (magnesium_ions) then
         do ii=1,n_magnesium_ions
            group_csn_magnesium(ip,ii) = SED_table_magnesium(1,1,ip,  2*ii-1)
            group_cse_magnesium(ip,ii) = SED_table_magnesium(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (silicon_ions) then
         do ii=1,n_silicon_ions
            group_csn_silicon(ip,ii) = SED_table_silicon(1,1,ip,  2*ii-1)
            group_cse_silicon(ip,ii) = SED_table_silicon(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (sulfur_ions) then
         do ii=1,n_sulfur_ions
            group_csn_sulfur(ip,ii) = SED_table_sulfur(1,1,ip,  2*ii-1)
            group_cse_sulfur(ip,ii) = SED_table_sulfur(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (iron_ions) then
         do ii=1,n_iron_ions
            group_csn_iron(ip,ii) = SED_table_iron(1,1,ip,  2*ii-1)
            group_cse_iron(ip,ii) = SED_table_iron(1,1,ip,1+2*ii-1)
         enddo
        endif

        if (neon_ions) then
         do ii=1,n_neon_ions
            group_csn_neon(ip,ii) = SED_table_neon(1,1,ip,  2*ii-1)
            group_cse_neon(ip,ii) = SED_table_neon(1,1,ip,1+2*ii-1)
         enddo
        endif

        group_csn_dust(ip,1) = SED_table_dust(1,1,ip,  1)

     endif
  end do
  call updateRTgroups_coolConstants
  if(mod(nstep_coarse,ncontrol)==0.and.myid==1) write(*,*) &
                    'SED Photon groups updated through stellar polling'
  !call write_group_props(.true.,6)

END SUBROUTINE update_SED_group_props

!*************************************************************************
SUBROUTINE update_star_RT_feedback(ilevel)

! Turn on RT advection if needed.
! Update photon group properties from stellar populations.
!-------------------------------------------------------------------------
  use amr_parameters
  use amr_commons
  use rt_parameters
  use pm_commons
  integer::ilevel
  logical,save::groupProps_init=.false.
!-------------------------------------------------------------------------
  if(rt_star.and. nstar_tot .gt. 0) then
     if(.not.rt_advect) then ! Turn on RT advection due to newborn stars:
        if(myid==1) write(*,*) '*****************************************'
        if(myid==1) write(*,*) 'Stellar RT turned on at a=',aexp
        if(myid==1) write(*,*) '*****************************************'
        rt_advect=.true.
     endif
     ! Set group props from stellar populations:
     if(sedprops_update .gt. 0 .and. .not.groupProps_init) then
        call update_SED_group_props
        groupProps_init=.true.
     else if(sedprops_update .gt. 0 .and. ilevel==levelmin &
          .and. mod(nstep_coarse,sedprops_update)==0) then
        call update_SED_group_props
     endif
  endif
END SUBROUTINE update_star_RT_feedback

!*************************************************************************
SUBROUTINE star_RT_feedback(ilevel, dt)

! This routine adds photons from radiating stars to appropriate cells in
! the  hydro grid.
! ilegel =>  grid level in which to perform the feedback
! ti     =>  initial time for the timestep (code units)
! dt     =>  real timestep length in code units
!-------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use rt_parameters
  integer:: ilevel
  real(dp):: dt
  integer:: igrid, jgrid, ipart, jpart, next_part
  integer:: i, ig, ip, npart1, npart2, icpu
  integer,dimension(1:nvector),save:: ind_grid, ind_part, ind_grid_part
!-------------------------------------------------------------------------
#if NGROUPS > 0
  if(.not.rt_advect)RETURN
  if(nstar_tot .le. 0 ) return
  if(numbtot(1,ilevel)==0)return ! number of grids in the level
  if(verbose)write(*,111)ilevel
  ! Gather star particles only.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel) ! grid index
     ig=0                     ! index of grid with stars (in ind_grid)
     ip=0                     ! index of star particle   (in ind_part)
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)   ! Number of particles in the grid
        npart2=0              ! number of selected (i.e. star) particles

        ! Count star particles in the grid
        if(npart1 > 0)then
          ipart = headp(igrid)        ! particle index
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle       <--- Very important !!!
              next_part = nextp(ipart)
              if (is_star(typep(ipart))) then
                 npart2 = npart2+1     ! only stars
              endif
              ipart = next_part        ! Go to next particle
           end do
        endif

        ! Gather star particles within the grid
        if(npart2 > 0)then
           ig = ig+1
           ind_grid(ig) = igrid
           ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle      <--- Very important !!!
              next_part = nextp(ipart)
              ! Select only star particles
              if (is_star(typep(ipart))) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip = ip+1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig ! points to grid a star is in
              endif
              if(ip == nvector)then
                 if(static) then
                    call star_RT_vsweep_pp( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 else
                    call star_RT_vsweep( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 endif
                 ip = 0
                 ig = 0
              end if
              ipart = next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid = next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip > 0) then
        if(static) then
           call star_RT_vsweep_pp( &
                ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
        else
           call star_RT_vsweep( &
                ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
        endif
     endif
  end do
  ! End loop over cpus

111 format('   Entering star_rt_feedback for level ',I2)
#endif
END SUBROUTINE star_RT_feedback

!*************************************************************************
!*************************************************************************
! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

!*************************************************************************
FUNCTION getSEDLuminosity(X, Y, N, e0, e1)

! Compute and return luminosity in energy interval (e0,e1) [eV]
! in SED Y(X). Assumes X is in Angstroms and Y in Lo/Angstroms/Msun.
! (Lo=[Lo_sun], Lo_sun=[erg s-1]. total solar luminosity is
! Lo_sun=10^33.58 erg/s)
! returns: Photon luminosity in, [# s-1 Msun-1],
!                             or [eV s-1 Msun-1] if SED_isEgy=true
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, hp, ev_to_erg, SED_isEgy
  use spectrum_integrator_module
  real(kind=8):: getSEDLuminosity, X(n), Y(n), e0, e1
  integer :: N, species
  real(kind=8),parameter :: const=1.0e-8/hp/c_cgs
  real(kind=8),parameter :: Lsun=3.8256d33 ! Solar luminosity [erg/s]
  ! const is a div by ph energy => ph count.  1e-8 is a conversion into
  ! cgs, since wly=[angstrom] h=[erg s-1], c=[cm s-1]
!-------------------------------------------------------------------------
  species          = 1                   ! irrelevant but must be included
  if(.not. SED_isEgy) then               !  Photon number per sec per Msun
     getSEDLuminosity = const &
          * integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
     getSEDLuminosity = getSEDLuminosity*Lsun  ! Scale by solar luminosity
  else                             ! SED_isEgy=true -> eV per sec per Msun
     getSEDLuminosity = integrateSpectrum(X, Y, N, e0, e1, species, f1)
     ! Scale by solar lum and convert to eV (bc group energies are in eV)
     getSEDLuminosity = getSEDLuminosity/eV_to_erg*Lsun
  endif
END FUNCTION getSEDLuminosity

!*************************************************************************
FUNCTION getSEDEgy(X, Y, N, e0, e1)

! Compute average energy, in eV, in energy interval (e0,e1) [eV] in SED
! Y(X). Assumes X is in Angstroms and Y is energy weight per angstrom
! (not photon count).
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, eV_to_erg, hp
  use spectrum_integrator_module
  real(dp):: getSEDEgy, X(N), Y(N), e0, e1, norm
  integer :: N,species
  real(dp),parameter :: const=1.d8*hp*c_cgs/eV_to_erg! energy conversion
!-------------------------------------------------------------------------
  species      = 1                       ! irrelevant but must be included
  norm         = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
  getSEDEgy    = const * &
                 integrateSpectrum(X, Y, N, e0, e1, species, f1) / norm
END FUNCTION getSEDEgy

!*************************************************************************
FUNCTION getSEDcsn(X, Y, N, e0, e1, species)

! Compute and return average photoionization
! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
! angstrom (not photon #).
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  use spectrum_integrator_module
  use rt_parameters,only:ionEVs
  real(kind=8):: getSEDcsn, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getSEDcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
  getSEDcsn= integrateSpectrum(X, Y, N, e0, e1, species, fSigLambda)/norm
END FUNCTION getSEDcsn

!************************************************************************
FUNCTION getSEDcse(X, Y, N, e0, e1, species)

! Compute and return average energy weighted photoionization
! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
! angstrom (not photon #).
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  use spectrum_integrator_module
  use rt_parameters,only:ionEVs
  real(dp):: getSEDcse, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getSEDcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm      = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getSEDcse = integrateSpectrum(X, Y, N, e0, e1, species, fSig) / norm
END FUNCTION getSEDcse

FUNCTION getSEDcsn_metals(X, Y, N, e0, e1, species, ion_state)

   ! Compute and return average photoionization
   ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
   ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
   ! angstrom (not photon #).
   ! Species is a code for the ion in question: 1=oxygen, 2=nitrogen, 3=carbon, 4=magnesium, 5=silicon, 6=sulfur, 7=iron, 8=neon
   !-------------------------------------------------------------------------
     use spectrum_integrator_module
     use rt_parameters,only:oxygen_ionEvs,nitrogen_ionEvs,carbon_ionEvs,magnesium_ionEvs,silicon_ionEvs,sulfur_ionEvs,iron_ionEvs,neon_ionEvs
     real(kind=8):: getSEDcsn_metals, X(N), Y(N), e0, e1, norm
     integer :: N, species, ion_state
   !-------------------------------------------------------------------------
     if (species .eq. 1) then
      if (e1 .gt. 0. .and. e1 .le. oxygen_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 2) then
      if (e1 .gt. 0. .and. e1 .le. nitrogen_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 3) then
      if (e1 .gt. 0. .and. e1 .le. carbon_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 4) then
      if (e1 .gt. 0. .and. e1 .le. magnesium_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 5) then
      if (e1 .gt. 0. .and. e1 .le. silicon_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 6) then
      if (e1 .gt. 0. .and. e1 .le. sulfur_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 7) then
      if (e1 .gt. 0. .and. e1 .le. iron_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 8) then
      if (e1 .gt. 0. .and. e1 .le. neon_ionEvs(ion_state)) then
         getSEDcsn_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     endif
     norm     = integrateSpectrum_metals(X, Y, N, e0, e1, species, ion_state, fLambda_metals)
     getSEDcsn_metals= integrateSpectrum_metals(X, Y, N, e0, e1, species, ion_state, fSigLambda_metals)/norm
END FUNCTION getSEDcsn_metals

   !************************************************************************
FUNCTION getSEDcse_metals(X, Y, N, e0, e1, species, ion_state)

   ! Compute and return average energy weighted photoionization
   ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
   ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
   ! angstrom (not photon #).
   ! Species is a code for the ion in question: 1=oxygen, 2=nitrogen, 3=carbon, 4=magnesium, 5=silicon, 6=sulfur, 7=iron, 8=neon
   !-------------------------------------------------------------------------
     use spectrum_integrator_module
     use rt_parameters,only:oxygen_ionEvs,nitrogen_ionEvs,carbon_ionEvs,magnesium_ionEvs,silicon_ionEvs,sulfur_ionEvs,iron_ionEvs,neon_ionEvs
     real(dp):: getSEDcse_metals, X(N), Y(N), e0, e1, norm
     integer :: N, species, ion_state
   !-------------------------------------------------------------------------
     if (species .eq. 1) then
      if (e1 .gt. 0. .and. e1 .le. oxygen_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 2) then
      if (e1 .gt. 0. .and. e1 .le. nitrogen_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 3) then
      if (e1 .gt. 0. .and. e1 .le. carbon_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 4) then
      if (e1 .gt. 0. .and. e1 .le. magnesium_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 5) then
      if (e1 .gt. 0. .and. e1 .le. silicon_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 6) then
      if (e1 .gt. 0. .and. e1 .le. sulfur_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 7) then
      if (e1 .gt. 0. .and. e1 .le. iron_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     else if (species .eq. 8) then
      if (e1 .gt. 0. .and. e1 .le. neon_ionEvs(ion_state)) then
         getSEDcse_metals=0. ; RETURN    ! [e0,e1] below ionization energy of species
      endif
     endif
     norm      = integrateSpectrum_metals(X, Y, N, e0, e1, species, ion_state, f1_metals)
     getSEDcse_metals = integrateSpectrum_metals(X, Y, N, e0, e1, species, ion_state, fSig_metals) / norm
END FUNCTION getSEDcse_metals

FUNCTION getSEDcsn_dust(X, Y, N, e0, e1, species)

   ! Compute and return average effective absorption cross section for dust
   ! in cm^2 H^-1, for a given energy interval (e0,e1) [eV] in
   ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
   ! angstrom (not photon #).
   ! Species is a dummy variable here as we only have one cross section at the
   ! moment
   !-------------------------------------------------------------------------
     use spectrum_integrator_module
     use rt_parameters,only:ionEVs
     real(kind=8):: getSEDcsn_dust, X(N), Y(N), e0, e1, norm
     integer :: N, species
   !-------------------------------------------------------------------------
     norm     = integrateSpectrum(X, Y, N, e0, e1, species, fLambda_dust)
     getSEDcsn_dust= integrateSpectrum(X, Y, N, e0, e1, species, fSigLambda_dust)/norm
END FUNCTION getSEDcsn_dust

!*************************************************************************
SUBROUTINE rebin_log(xint_log, yint_log,                                 &
               data,       nx,       ny,     x,     y,     nz,           &
               new_data, new_nx, new_ny, new_x, new_y      )

! Rebin the given 2d data into constant logarithmic intervals, using
! linear interpolation.
! xint_log,  => x and y intervals in the rebinned data. If negative,
! yint_log      these values represent the new number of bins.
! data       => The 2d data to be rebinned
! nx,ny      => Number of points in x and y in the original data
! nz         => Number of values in the data
! x,y        => x and y values for the data
! new_data  <=  The rebinned 2d data
! new_nx,ny <=  Number of points in x and y in the rebinned data
! new_x,y   <=  x and y point values for the rebinned data
!-------------------------------------------------------------------------
  real(kind=8):: xint_log, yint_log
  integer,intent(in):: nx, ny, nz
  integer::new_nx, new_ny
  real(kind=8),intent(in):: x(nx),y(ny)
  real(kind=8),intent(in):: data(nx,ny,nz)
  real(kind=8),dimension(:,:,:),allocatable:: new_data
  real(dp),dimension(:),allocatable:: new_x, new_lgx, new_y, new_lgy
  real(kind=8):: dx0, dx1, dy0, dy1, x_step, y_step
  real(kind=8):: x0lg, x1lg, y0lg, y1lg
  integer :: i, j, ix, iy, ix1, iy1
!-------------------------------------------------------------------------
  if(allocated(new_x)) deallocate(new_x)
  if(allocated(new_y)) deallocate(new_y)
  if(allocated(new_data)) deallocate(new_data)

  ! Find dimensions of the new_data and initialize it
  x0lg = log10(x(1));   x1lg = log10(x(nx))
  y0lg = log10(y(1));   y1lg = log10(y(ny))

  if(xint_log .lt. 0 .and. nx .gt. 1) then
     new_nx=-xint_log                             ! xint represents wanted
     xint_log = (x1lg-x0lg)/(new_nx-1)            !     number of new bins
  else
     new_nx = (x1lg-x0lg)/xint_log + 1
  endif
  allocate(new_x(new_nx)) ; allocate(new_lgx(new_nx))
  do i = 0, new_nx-1                              !  initialize the x-axis
     new_lgx(i+1) = x0lg + i*xint_log
  end do
  new_x=10.d0**new_lgx

  if(yint_log .lt. 0 .and. ny .gt. 1) then        ! yint represents wanted
     new_ny=-yint_log                             !     number of new bins
     yint_log = (y1lg-y0lg)/(new_ny-1)
  else
     new_ny = (y1lg-y0lg)/yint_log + 1
  endif
  allocate(new_y(new_ny)) ; allocate(new_lgy(new_ny))
  do j = 0, new_ny-1                              !      ...and the y-axis
     new_lgy(j+1) = y0lg + j*yint_log
  end do
  new_y=10.d0**new_lgy

  ! Initialize new_data and find values for each point in it
  allocate(new_data(new_nx, new_ny, nz))
  do j = 1, new_ny
     call locate(y, ny, new_y(j), iy)
     ! y(iy) <= new_y(j) <= y(iy+1)
     ! iy is lower bound, so it can be zero but not larger than ny
     if(iy < 1) iy=1
     if (iy < ny) then
        iy1  = iy + 1
        y_step = y(iy1) - y(iy)
        dy0  = max(new_y(j) - y(iy),    0.0d0)  / y_step
        dy1  = min(y(iy1)   - new_y(j), y_step) / y_step
     else
        iy1  = iy
        dy0  = 0.0d0 ;  dy1  = 1.0d0
     end if

     do i = 1, new_nx
        call locate(x, nx, new_x(i), ix)
        if(ix < 1) ix=1
        if (ix < nx) then
           ix1  = ix+1
           x_step = x(ix1)-x(ix)
           dx0  = max(new_x(i) - x(ix),    0.0d0)  / x_step
           dx1  = min(x(ix1)   - new_x(i), x_step) / x_step
        else
           ix1  = ix
           dx0  = 0.0d0 ;  dx1  = 1.0d0
        end if

        if (abs(dx0+dx1-1.0d0) .gt. 1.0d-6 .or.                          &
                                         abs(dy0+dy1-1.0d0) > 1.0d-6) then
           write(*,*) 'Screwed up the rebin interpolation ... '
           write(*,*) dx0+dx1,dy0+dy1
           call clean_stop
        end if

        new_data(i,j,:) =                                                &
             dx0 * dy0 * data(ix1,iy1,:) + dx1 * dy0 * data(ix, iy1,:) + &
             dx0 * dy1 * data(ix1,iy, :) + dx1 * dy1 * data(ix, iy, :)
     end do
  end do

  deallocate(new_lgx)
  deallocate(new_lgy)

END SUBROUTINE rebin_log

!*************************************************************************
SUBROUTINE write_SEDtable()

! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           write(10,900)                                                 &
                SED_ages(i)        ,    SED_zeds(j)        ,             &
                SED_table(i,j,ip,1),    SED_table(i,j,ip,2),             &
                SED_table(i,j,ip,3),    SED_table(i,j,ip,4),             &
                SED_table(i,j,ip,5),    SED_table(i,j,ip,6),             &
                SED_table(i,j,ip,7),    SED_table(i,j,ip,8),             &
                SED_table(i,j,ip,9)
        end do
     end do
     close(10)
  end do
900 format (ES15.4, ES15.4, ES15.4, ES15.4, f15.4                        &
           ,ES15.4, ES15.4, ES15.4, ES15.4, ES15.4, ES15.4)

END SUBROUTINE write_SEDtable

!*************************************************************************
SUBROUTINE write_SEDtable_oxygen()
  use hydro_parameters, only: n_oxygen_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_oxygen', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_oxygen_ions-1
              write(10,900,advance="no") SED_table_oxygen(i,j,ip,k)
           end do
           write(10,900) SED_table_oxygen(i,j,ip,2*n_oxygen_ions)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_oxygen
!*************************************************************************
SUBROUTINE write_SEDtable_nitrogen()
  use hydro_parameters, only: n_nitrogen_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_nitrogen', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_nitrogen_ions-1
              write(10,900,advance="no") SED_table_nitrogen(i,j,ip,k)
           end do
           write(10,900) SED_table_nitrogen(i,j,ip,2*n_nitrogen_ions)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_nitrogen
!*************************************************************************
SUBROUTINE write_SEDtable_carbon()
  use hydro_parameters, only: n_carbon_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_carbon', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_carbon_ions-1
              write(10,900,advance="no") SED_table_carbon(i,j,ip,k)
           end do
           write(10,900) SED_table_carbon(i,j,ip,2*n_carbon_ions)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_carbon
!*************************************************************************
SUBROUTINE write_SEDtable_magnesium()
  use hydro_parameters, only: n_magnesium_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_magnesium', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_magnesium_ions-1
              write(10,900,advance="no") SED_table_magnesium(i,j,ip,k)
           end do
           write(10,900) SED_table_magnesium(i,j,ip,2*n_magnesium_ions)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_magnesium
!*************************************************************************
SUBROUTINE write_SEDtable_silicon()
  use hydro_parameters, only: n_silicon_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_silicon', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_silicon_ions-1
              write(10,900,advance="no") SED_table_silicon(i,j,ip,k)
           end do
           write(10,900) SED_table_silicon(i,j,ip,2*n_silicon_ions)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_silicon
!*************************************************************************
SUBROUTINE write_SEDtable_sulfur()
  use hydro_parameters, only: n_sulfur_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_sulfur', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_sulfur_ions-1
              write(10,900,advance="no") SED_table_sulfur(i,j,ip,k)
           end do
           write(10,900) SED_table_sulfur(i,j,ip,2*n_sulfur_ions-1)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_sulfur
!*************************************************************************
SUBROUTINE write_SEDtable_iron()
  use hydro_parameters, only: n_iron_ions
! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j, k
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') './SEDtables/SEDtable_iron', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ

     do j = 1,SED_nz
        do i = 1,SED_nA
           do k=1,2*n_iron_ions-1
              write(10,900,advance="no") SED_table_iron(i,j,ip,k)
           end do
           write(10,900) SED_table_iron(i,j,ip,2*n_iron_ions)
        end do
     end do
     close(10)
  end do
900 format (ES15.4)

END SUBROUTINE write_SEDtable_iron
!*************************************************************************
SUBROUTINE write_SEDtable_neon()
   use hydro_parameters, only: n_neon_ions
 ! Write the SED properties to a file (this is just in debugging, to check
 ! if the SEDs are being read correctly).
 !-------------------------------------------------------------------------
   character(len=128)::filename
   integer::ip, i, j, k
 !-------------------------------------------------------------------------
   do ip=1,nSEDgroups
      write(filename,'(A, I1, A)') './SEDtables/SEDtable_neon', ip, '.list'
      open(10, file=filename, status='unknown')
      write(10,*) SED_nA, SED_nZ

      do j = 1,SED_nz
         do i = 1,SED_nA
            do k=1,2*n_neon_ions-1
               write(10,900,advance="no") SED_table_neon(i,j,ip,k)
            end do
            write(10,900) SED_table_neon(i,j,ip,2*n_neon_ions)
         end do
      end do
      close(10)
   end do
 900 format (ES15.4)

 END SUBROUTINE write_SEDtable_neon
!*************************************************************************
SUBROUTINE write_SEDtable_dust()
#ifdef CALIMA
   use dust_commons,only:dust_ratd,dust_pe_heating
   use hydro_parameters, only:ndust,npah
#endif
   ! Write the SED-averaged dust optical properties to a file (this is just
   !  in debugging, to check if the SEDs are being read correctly).
   !-------------------------------------------------------------------------
     character(len=128)::filename
     integer::ip, i, j,k
   !-------------------------------------------------------------------------
#ifdef CALIMA
   if (ndust>0) then
      do ip=1,nSEDgroups
         write(filename,'(A, I1, A)') './SEDtables/SEDtable_dust_Ab', ip, '.list'
         open(10, file=filename, status='unknown')
         write(10,*) SED_nA, SED_nZ
         
         do j = 1,SED_nz
            do i = 1,SED_nA
               do k=1,ndust
                  write(10,600)                                                 &
                        SED_ages(i)        ,    SED_zeds(j)        ,             &
                        SED_table_dust_Ab(i,j,ip,k)
               end do
            end do
         end do
         close(10)
      end do
      600 format (ES15.4, ES15.4, ES15.4)

      do ip=1,nSEDgroups
         write(filename,'(A, I1, A)') './SEDtables/SEDtable_dust_Sc', ip, '.list'
         open(10, file=filename, status='unknown')
         write(10,*) SED_nA, SED_nZ
         
         do j = 1,SED_nz
            do i = 1,SED_nA
               do k=1,ndust+npah*2
                  write(10,700)                                                 &
                        SED_ages(i)        ,    SED_zeds(j)        ,             &
                        SED_table_dust_Sc(i,j,ip,k)
               end do
            end do
         end do
         close(10)
      end do
      700 format (ES15.4, ES15.4, ES15.4)

      do ip=1,nSEDgroups
         write(filename,'(A, I1, A)') './SEDtables/SEDtable_dust_Rp', ip, '.list'
         open(10, file=filename, status='unknown')
         write(10,*) SED_nA, SED_nZ
         
         do j = 1,SED_nz
            do i = 1,SED_nA
               do k=1,ndust+npah*2
                  write(10,800)                                                 &
                        SED_ages(i)        ,    SED_zeds(j)        ,             &
                        SED_table_dust_Rp(i,j,ip,k)
               end do
            end do
         end do
         close(10)
      end do
      800 format (ES15.4, ES15.4, ES15.4)
      if (dust_ratd) then
         do ip=1,nSEDgroups
            write(filename,'(A, I1, A)') './SEDtables/SEDtable_dust_RAT', ip, '.list'
            open(10, file=filename, status='unknown')
            write(10,*) SED_nA, SED_nZ
            
            do j = 1,SED_nz
               do i = 1,SED_nA
                  do k=1,ndust
                     write(10,750)                                                 &
                           SED_ages(i)        ,    SED_zeds(j)        ,             &
                           SED_table_dust_RAT(i,j,ip,k)
                  end do
               end do
            end do
            close(10)
         end do
         750 format (ES15.4, ES15.4, ES15.4)
      end if
      if (dust_pe_heating) then
         do ip=1,nSEDgroups
            write(filename,'(A, I1, A)') './SEDtables/SEDtable_dust_la', ip, '.list'
            open(10, file=filename, status='unknown')
            write(10,*) SED_nA, SED_nZ
            
            do j = 1,SED_nz
               do i = 1,SED_nA
                  do k=1,ndust
                     write(10,850)                                                 &
                           SED_ages(i)        ,    SED_zeds(j)        ,             &
                           SED_table_dust_la(i,j,ip,k)
                  end do
               end do
            end do
            close(10)
         end do
         850 format (ES15.4, ES15.4, ES15.4)
      end if
   end if
   if (npah>0) then
      do ip=1,nSEDgroups
         write(filename,'(A, I1, A)') './SEDtables/SEDtable_PAH_Ab', ip, '.list'
         open(10, file=filename, status='unknown')
         write(10,*) SED_nA, SED_nZ
         
         do j = 1,SED_nz
            do i = 1,SED_nA
               do k=1,npah*2
                  write(10,650)                                                 &
                        SED_ages(i)        ,    SED_zeds(j)        ,             &
                        SED_table_PAH_abs(i,j,ip,k)
               end do
            end do
         end do
         close(10)
      end do
      650 format (ES15.4, ES15.4, ES15.4)

      do ip=1,nSEDgroups
         write(filename,'(A, I1, A)') './SEDtables/SEDtable_PAH_Sc', ip, '.list'
         open(10, file=filename, status='unknown')
         write(10,*) SED_nA, SED_nZ
         
         do j = 1,SED_nz
            do i = 1,SED_nA
               do k=1,npah*2
                  write(10,650)                                                 &
                        SED_ages(i)        ,    SED_zeds(j)        ,             &
                        SED_table_PAH_sc(i,j,ip,k)
               end do
            end do
         end do
         close(10)
      end do

      do ip=1,nSEDgroups
         write(filename,'(A, I1, A)') './SEDtables/SEDtable_PAH_Rp', ip, '.list'
         open(10, file=filename, status='unknown')
         write(10,*) SED_nA, SED_nZ
         
         do j = 1,SED_nz
            do i = 1,SED_nA
               do k=1,npah*2
                  write(10,650)                                                 &
                        SED_ages(i)        ,    SED_zeds(j)        ,             &
                        SED_table_PAH_Rp(i,j,ip,k)
               end do
            end do
         end do
         close(10)
      end do
   end if
#else
   do ip=1,nSEDgroups
      write(filename,'(A, I1, A)') './SEDtables/SEDtable_dust', ip, '.list'
      open(10, file=filename, status='unknown')
      write(10,*) SED_nA, SED_nZ
      
      do j = 1,SED_nz
         do i = 1,SED_nA
            write(10,900)                                                 &
                  SED_ages(i)        ,    SED_zeds(j)        ,             &
                  SED_table_dust(i,j,ip,1)
         end do
      end do
      close(10)
   end do
   900 format (ES15.4, ES15.4, ES15.4)
#endif
   
END SUBROUTINE write_SEDtable_dust

!*************************************************************************
SUBROUTINE inp_SED_table(age, Z, nProp, same, ret)

! Compute SED property by interpolation from table.
! input/output:
! age   => Star population age [Gyrs]
! Z     => Star population metallicity [m_metals/m_tot]
! nprop => Number of property to fetch
!          1=log(photon # intensity [# Msun-1 s-1]),
!          2=log(cumulative photon # intensity [# Msun-1]),
!          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse
! same  => If true then assume same age and Z as used in last call.
!          In this case the interpolation indexes can be recycled.
! ret   => The interpolated values of the sed property for every photon
!          group
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  real(dp), intent(in):: age, Z
  real(dp):: lgAge, lgZ
  integer:: nProp
  logical:: same
  real(dp),dimension(:):: ret
  integer,save:: ia, iz
  real(dp),save:: da, da0, da1, dz, dz0, dz1
!-------------------------------------------------------------------------
  ! ia, iz: lower indexes: 0<ia<sed_nA etc.
  ! da0, da1, dz0, dz1: proportional distances from edges:
  ! 0<=da0<=1, 0<=da1<=1 etc.
  if(.not. same) then
     if(age.le.0d0) then
        lgAge=-4d0
     else
        lgAge = log10(age)
     endif
     lgZ=log10(Z)
     ia = min(max(floor((lgAge-SED_lgA0)/SED_dlgA ) + 2, 1  ),  SED_nA-1 )
     da = SED_ages(ia+1)-SED_ages(ia)
     da0= min( max(   (age-SED_ages(ia)) /da,       0. ), 1.          )
     da1= min( max(  (SED_ages(ia+1)-age)/da,       0. ), 1.          )

     iz = min(max(floor((lgZ-SED_lgZ0)/SED_dlgZ ) + 1,   1  ),  SED_nZ-1 )
     dz = sed_Zeds(iz+1)-SED_Zeds(iz)
     dz0= min( max(   (Z-SED_zeds(iz)) /dz,         0. ),  1.         )
     dz1= min( max(  (SED_Zeds(iz+1)-Z)/dz,         0. ),  1.         )

     if (abs(da0+da1-1.0d0) > 1.0d-5 .or. abs(dz0+dz1-1.0d0) > 1.0d-5) then
        write(*,*) 'Screwed up the sed interpolation ... '
        write(*,*) da0+da1,dz0+dz1
        call clean_stop
     end if
  endif

  ret = da0 * dz0 * SED_table(ia+1, iz+1, :, nProp) + &
        da1 * dz0 * SED_table(ia,   iz+1, :, nProp) + &
        da0 * dz1 * SED_table(ia+1, iz,   :, nProp) + &
        da1 * dz1 * SED_table(ia,   iz,   :, nProp)

END SUBROUTINE inp_SED_table

!*************************************************************************
SUBROUTINE inp_SED_table_metals(age, Z, nProp, same, species, ion_state, ret)

   ! Compute SED property by interpolation from table.
   ! input/output:
   ! age   => Star population age [Gyrs]
   ! Z     => Star population metallicity [m_metals/m_tot]
   ! nprop => Number of property to fetch
   !          1=log(photon # intensity [# Msun-1 s-1]),
   !          2=log(cumulative photon # intensity [# Msun-1]),
   !          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse
   ! same  => If true then assume same age and Z as used in last call.
   !          In this case the interpolation indexes can be recycled.
   ! ret   => The interpolated values of the sed property for every photon
   !          group
   !-------------------------------------------------------------------------
     use amr_commons
     use rt_parameters
     real(dp), intent(in):: age, Z
     real(dp):: lgAge, lgZ
     integer:: nProp, species, ion_state
     logical:: same
     real(dp),dimension(:):: ret
     integer,save:: ia, iz
     real(dp),save:: da, da0, da1, dz, dz0, dz1
   !-------------------------------------------------------------------------
     ! ia, iz: lower indexes: 0<ia<sed_nA etc.
     ! da0, da1, dz0, dz1: proportional distances from edges:
     ! 0<=da0<=1, 0<=da1<=1 etc.
     if(.not. same) then
        if(age.le.0d0) then
           lgAge=-4d0
        else
           lgAge = log10(age)
        endif
        lgZ=log10(Z)
        ia = min(max(floor((lgAge-SED_lgA0)/SED_dlgA ) + 2, 1  ),  SED_nA-1 )
        da = SED_ages(ia+1)-SED_ages(ia)
        da0= min( max(   (age-SED_ages(ia)) /da,       0. ), 1.          )
        da1= min( max(  (SED_ages(ia+1)-age)/da,       0. ), 1.          )

        iz = min(max(floor((lgZ-SED_lgZ0)/SED_dlgZ ) + 1,   1  ),  SED_nZ-1 )
        dz = sed_Zeds(iz+1)-SED_Zeds(iz)
        dz0= min( max(   (Z-SED_zeds(iz)) /dz,         0. ),  1.         )
        dz1= min( max(  (SED_Zeds(iz+1)-Z)/dz,         0. ),  1.         )

        if (abs(da0+da1-1.0d0) > 1.0d-5 .or. abs(dz0+dz1-1.0d0) > 1.0d-5) then
           write(*,*) 'Screwed up the sed interpolation metal... '
           write(*,*) da0+da1,dz0+dz1
           call clean_stop
        end if
     endif

     if (species.eq.1) then ! oxygen
      ret = da0 * dz0 * SED_table_oxygen(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_oxygen(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_oxygen(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_oxygen(ia,   iz,   :, nProp)
     else if (species.eq.2) then ! nitrogen
      ret = da0 * dz0 * SED_table_nitrogen(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_nitrogen(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_nitrogen(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_nitrogen(ia,   iz,   :, nProp)
     else if (species.eq.3) then ! carbon
      ret = da0 * dz0 * SED_table_carbon(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_carbon(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_carbon(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_carbon(ia,   iz,   :, nProp)
     else if (species.eq.4) then ! magnesium
      ret = da0 * dz0 * SED_table_magnesium(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_magnesium(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_magnesium(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_magnesium(ia,   iz,   :, nProp)
     else if (species.eq.5) then ! silicon
      ret = da0 * dz0 * SED_table_silicon(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_silicon(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_silicon(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_silicon(ia,   iz,   :, nProp)
     else if (species.eq.6) then ! sulfur
      ret = da0 * dz0 * SED_table_sulfur(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_sulfur(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_sulfur(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_sulfur(ia,   iz,   :, nProp)
     else if (species.eq.7) then ! iron
      ret = da0 * dz0 * SED_table_iron(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_iron(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_iron(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_iron(ia,   iz,   :, nProp)
     else if (species.eq.8) then ! neon
      ret = da0 * dz0 * SED_table_neon(ia+1, iz+1, :, nProp) + &
            da1 * dz0 * SED_table_neon(ia,   iz+1, :, nProp) + &
            da0 * dz1 * SED_table_neon(ia+1, iz,   :, nProp) + &
            da1 * dz1 * SED_table_neon(ia,   iz,   :, nProp)
     endif

END SUBROUTINE inp_SED_table_metals

SUBROUTINE inp_SED_table_dust(age, Z, nProp, same, ret)

   ! Compute SED property by interpolation from table.
   ! input/output:
   ! age   => Star population age [Gyrs]
   ! Z     => Star population metallicity [m_metals/m_tot]
   ! nprop => Number of property to fetch
   !          1=log(photon # intensity [# Msun-1 s-1]),
   !          2=log(cumulative photon # intensity [# Msun-1]),
   !          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse
   ! same  => If true then assume same age and Z as used in last call.
   !          In this case the interpolation indexes can be recycled.
   ! ret   => The interpolated values of the sed property for every photon
   !          group
   !-------------------------------------------------------------------------
     use amr_commons
     use rt_parameters
     real(dp), intent(in):: age, Z
     real(dp):: lgAge, lgZ
     integer:: nProp
     logical:: same
     real(dp),dimension(:):: ret
     integer,save:: ia, iz
     real(dp),save:: da, da0, da1, dz, dz0, dz1
   !-------------------------------------------------------------------------
     ! ia, iz: lower indexes: 0<ia<sed_nA etc.
     ! da0, da1, dz0, dz1: proportional distances from edges:
     ! 0<=da0<=1, 0<=da1<=1 etc.
     if(.not. same) then
        if(age.le.0d0) then
           lgAge=-4d0
        else
           lgAge = log10(age)
        endif
        lgZ=log10(Z)
        ia = min(max(floor((lgAge-SED_lgA0)/SED_dlgA ) + 2, 1  ),  SED_nA-1 )
        da = SED_ages(ia+1)-SED_ages(ia)
        da0= min( max(   (age-SED_ages(ia)) /da,       0. ), 1.          )
        da1= min( max(  (SED_ages(ia+1)-age)/da,       0. ), 1.          )

        iz = min(max(floor((lgZ-SED_lgZ0)/SED_dlgZ ) + 1,   1  ),  SED_nZ-1 )
        dz = sed_Zeds(iz+1)-SED_Zeds(iz)
        dz0= min( max(   (Z-SED_zeds(iz)) /dz,         0. ),  1.         )
        dz1= min( max(  (SED_Zeds(iz+1)-Z)/dz,         0. ),  1.         )

        if (abs(da0+da1-1.0d0) > 1.0d-5 .or. abs(dz0+dz1-1.0d0) > 1.0d-5) then
           write(*,*) 'Screwed up the sed interpolation dust... '
           write(*,*) da0+da1,dz0+dz1
           call clean_stop
        end if
     endif

     ret = da0 * dz0 * SED_table_dust(ia+1, iz+1, :, nProp) + &
           da1 * dz0 * SED_table_dust(ia,   iz+1, :, nProp) + &
           da0 * dz1 * SED_table_dust(ia+1, iz,   :, nProp) + &
           da1 * dz1 * SED_table_dust(ia,   iz,   :, nProp)

   END SUBROUTINE inp_SED_table_dust

!*************************************************************************
SUBROUTINE getNPhotonsEmitted(age1_Gyr, dt_Gyr, Z, ret)

! Compute number of photons emitted by a stellar particle per solar mass
! over a timestep.
! input/output:
! age1_Gyr => Star population age [Gyrs] at timestep end
! dt_gyr   => Timestep length in Gyr
! Z        => Star population metallicity [m_metals/m_tot]
! ret      => # of photons emitted per solar mass over timestep
!-------------------------------------------------------------------------
  use rt_parameters
  real(dp),intent(in):: age1_Gyr, dt_Gyr, Z
  real(dp),dimension(nSEDgroups):: ret, Lc0, Lc1, SEDegy
!-------------------------------------------------------------------------
  ! Lc0 = cumulative emitted photons at the start of the timestep
  call inp_SED_table(age1_Gyr-dt_Gyr, Z, 2, .false., Lc0)
  ! Lc1 = cumulative emitted photons at the end of the timestep
  call inp_SED_table(age1_Gyr, Z, 2, .false., Lc1)
  ret = max(Lc1-Lc0,0.)
  if(SED_isEgy) then ! Integrate correct energy rather than # of photons
     ! Divide emitted energy by group energy -> Photon count
     ret = ret / group_egy(1:nSEDgroups)
  endif
END SUBROUTINE getNPhotonsEmitted

#if NGROUPS > 0
!*************************************************************************
SUBROUTINE star_RT_vsweep(ind_grid,ind_part,ind_grid_part,ng,np,dt,ilevel)

! This routine is called by subroutine star_rt_feedback.
! Each star particle dumps a number of photons into the nearest grid cell
! using array rtunew.
! Radiation is injected into cells at level ilevel, but it is important
! to know that ilevel-1 cells may also get some radiation. This is due
! to star particles that have just crossed to a coarser level.
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  star indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of stars
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,ngroups),save::part_NpInp
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0,igd=0,icd=0
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  &
            , scale_Fp, age, z, scale_inp, scale_Nphot, dt_Gyr           &
            , dt_loc_Gyr, scale_msun, mass
  real(dp),parameter::vol_factor=2**ndim   ! Vol factor for ilevel-1 cells
!-------------------------------------------------------------------------
  if(.not. metal) z = log10(max(z_ave*0.02, 10.d-5))![log(m_metals/m_tot)]
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)
  dt_Gyr = dt*scale_t*sec2Gyr
  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = scale_d / scale_np / vol_loc / m_sun
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50
  scale_msun = scale_d * scale_l**ndim / m_sun

  ! Lower left corners of 3x3x3 grid-cubes (with given grid in center)
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - 3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(&
          ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)
  ! now nbors_father cells are a cube of 27 cells with ind_cell in the
  ! middle and nbors_father_grids are the 8 grids that contain these 27
  ! cells (though only one of those is fully included in the cube)

  ! Rescale position of stars to positions within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx
        ! so 0<x<2 is bottom cell, ...., 4<x<6 is top cell
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! So id=0-5 is the cell (in the
     end do                    ! 3x3x3 supercube) containing the star
  end do

  ! Compute parent grids
  do idim = 1, ndim
     do j = 1, np
        igd(j,idim) = id(j,idim)/2 ! must be 0, 1 or 2
     end do
  end do
  do j = 1, np
     kg(j) = 1 + igd(j,1) + 3*igd(j,2) + 9*igd(j,3) ! 1 to 27
  end do
  do j = 1, np
     igrid(j) = son(nbors_father_cells(ind_grid_part(j),kg(j)))
     ! grid (not cell) containing the star
  end do

  ! Check if particles are entirely in level ilevel.
  ! This should always be the case in post-processing
  ok(1:np) = .true.
  do j = 1, np
     ok(j) = ok(j) .and. igrid(j) > 0
  end do ! if ok(j) is true then particle j's cell contains a grid.
  ! Otherwise it is a leaf, and we need to fill it with radiation.

  ! Compute parent cell position within it's grid
  do idim = 1, ndim
     do j = 1, np
        if( ok(j) ) then
           icd(j,idim) = id(j,idim) - 2*igd(j,idim) ! 0 or 1
        end if
     end do
  end do
  do j = 1, np
     if( ok(j) ) then
        icell(j) = 1 + icd(j,1) + 2*icd(j,2) + 4*icd(j,3) ! 1 to 8
     end if
  end do

  if(rt_nsubcycle .gt. 1) then
     ! Find all particles that have moved to a coarser level
     ! and assign their injection to the closest cell in
     ! the grid originally owning the particle.

     ! Compute parent cell position within ind_grid(ind_grid_part(j)):
     do idim = 1, ndim
        do j = 1, np
           if( .not. ok(j) ) then
              ! For each dim, if id=0,1,2 then icd=0,
              !               if id=3,4,5 then icd=1
              icd(j,idim) = id(j,idim)/3 ! 0 or 1
           end if
        end do
     end do
     do j = 1, np
        if( .not. ok(j) ) then
           icell(j) = 1 + icd(j,1) + 2*icd(j,2) + 4*icd(j,3) ! 1 to 8
           ! Use the original owner grid:
           igrid(j) = ind_grid(ind_grid_part(j))
           ok(j) = .true. ! So never injecting into coarser level
        end if
     end do
  endif

  ! Compute parent cell adress and particle radiation contribution
  do j = 1, np

     if(metal) then
        !z=zp(ind_part(j))
        z=2.09d0*zp(ind_part(j),2)+1.06d0*zp(ind_part(j),1) ! Metal enrichment
     else
        z=z_ave*0.02
     endif

     call getAgeGyr(tp(ind_part(j)), age)          !  End-of-dt age [Gyrs]
     ! Possibilities:     Born i) before dt, ii) within dt, iii) after dt:
     dt_loc_Gyr = max(min(dt_Gyr, age), 0.)

     if (is_pre_SN(typep(j))) then
        mass = mp(ind_part(j))
     else
        mass = mp(ind_part(j))/(1d0-eta_sn)
     endif
     if(use_initial_mass) mass = mp0(ind_part(j))

     if(pop3.and.is_pop_III(typep(ind_part(j)))) then
        ! Make sure index is less than 0 or else it exploded!
        if(is_pre_SN(typep(ind_part(j)))) then
           call getNPhotonsEmitted_pop3(mass,age,dt_loc_Gyr,part_NpInp(j,1:nSEDgroups))
           part_NpInp(j,1:nSEDgroups) = part_NpInp(j,1:nSEDgroups) ! escape fraction for PopIII
        else
           part_NpInp(j,1:nSEDgroups) = 0.d0
        endif
     else
        z= max(z, 10.d-5)                             !      [m_metals/m_tot]
        call getNPhotonsEmitted(age,dt_loc_Gyr,z,part_NpInp(j,1:nSEDgroups))
        part_NpInp(j,1:nSEDgroups) = part_NpInp(j,1:nSEDgroups)*rt_esc_frac ! escape fraction for PopII
     endif
     part_NpInp(j,:) = part_NpInp(j,:)*mass*scale_inp !#photons

     if(showSEDstats .and. nSEDgroups .gt. 0) then
        step_nPhot = step_nPhot+part_NpInp(j,4)*scale_nPhot
        step_nStar = step_nStar+dt_loc_Gyr*Gyr2sec/scale_t
        step_mStar = step_mStar+mp(ind_part(j)) * scale_msun             &
                                          * dt_loc_Gyr * Gyr2sec / scale_t
     endif

     if( ok(j) )then
        indp(j) = ncoarse + (icell(j)-1)*ngridmax + igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
     end if
  end do
  ! Increase photon density in cell due to stellar radiation
  do j=1,np
     if( ok(j) ) then                                      !   ilevel cell
        do ip=1,nSEDgroups
           rtunew(indp(j),iGroups(ip)) &
                = rtunew(indp(j),iGroups(ip))+part_NpInp(j,ip)
        end do
     else                                                  ! ilevel-1 cell
!        if (rt_nsubcycle == 1)then
           do ip=1,nSEDgroups
              rtunew(indp(j),iGroups(ip)) = rtunew(indp(j),iGroups(ip))  &
                   + part_NpInp(j,ip) / vol_factor
           end do
!        end if
     endif
     !begin debug
     !     if(rtunew(indp(j),iGroups(1))*scale_np*rt_c_cgs .gt.1.d11) then
     !        write(*,777),ilevel, myid, rtunew(indp(j),                   &
     !           iGroups(1))*scale_np*rt_c_cgs, &
     !             irad(j,1), irad(j,2),                                 &
     !             mp(ind_part(j)),                                      &
     !             mp(ind_part(j))*scale_d*scale_l**3/m_sun,             &
     !             vol_loc*scale_l**3,                                   &
     !             dt, scale_t,                                          &
     !             ok(j)
     !        print*,'*******************************'
     !     endif
     !end debug
     !777 format('HERE:',I2, ' ', I2, 8(1pe11.3), L, ' AHA')
  end do

END SUBROUTINE star_RT_vsweep

!*************************************************************************
SUBROUTINE star_RT_vsweep_pp(ind_grid, ind_part, ind_grid_part, ng, np,  &
                             dt, ilevel)
! This routine is called by subroutine star_rt_feedback,
! in the case of rt-post-processing. The difference from doing coupled
! RT is that in the post-processing case, the pointer from a grid to a
! particle is always true, so we can trust completely that
! ind_grid(ind_grid_part(j)) refers to the grid that particle ind_part(j)
! resides in, and we don't need to bother searching for the particle's
! grid.  Each star particle dumps a number of photons into the nearest
! grid cell using array rtunew.
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  star indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of stars
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  real(dp),dimension(1:nvector,nGroups),save::part_NpInp
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0
  !integer ,dimension(1:nvector,1:ndim),save::id
  integer ,dimension(1:nvector),save::icell,indp
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  &
            , scale_Fp, age, z, scale_inp, scale_Nphot, dt_Gyr           &
            , dt_loc_Gyr, scale_msun, mass
!-------------------------------------------------------------------------
  if(.not. metal) z = log10(max(z_ave*0.02, 10.d-5))![log(m_metals/m_tot)]
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_Fp)
  dt_Gyr = dt*scale_t*sec2Gyr
  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = scale_d / scale_np / vol_loc / m_sun
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50
  scale_msun = scale_d * scale_l**ndim / m_sun

  ! Lower left corners of parent grid
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - dx
     end do
  end do

  ! Rescale position of star to position within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx
        ! so 0<x<1 is bottom cell, ...., 1<x<2 is top cell
     end do
  end do

  ! Compute position of parent cell within parent grid
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! id=0-1 is cell in parent containing star
     end do
  end do
  do j = 1, np
     icell(j) = 1 + id(j,1) + 2*id(j,2) + 4*id(j,3) ! 1 to 8
  end do

  ! Compute parent cell adress and particle effective mass
  do j = 1, np
     call getAgeGyr(tp(ind_part(j)), age)         !   End-of-dt age [Gyrs]
     !Possibilities:      Born i) before dt, ii) within dt, iii) after dt:
     dt_loc_Gyr = max(min(dt_Gyr, age), 0.)

     if(is_pre_SN(typep(j)))then
        mass = mp(ind_part(j))
     else
        mass = mp(ind_part(j))/(1d0-eta_sn)
     endif
     if(use_initial_mass) mass = mp0(ind_part(j))

     if(pop3.and.is_pop_III(typep(j)))then
        if(is_pre_SN(typep(ind_part(j)))) then 
           call getNPhotonsEmitted_pop3(mass,age,dt_loc_Gyr,part_NpInp(j,1:nSEDgroups))
           part_NpInp(j,1:nSEDgroups) = part_NpInp(j,1:nSEDgroups) ! escape fraction for PopIII
        else
           part_NpInp(j,1:nSEDgroups) = 0.d0 
        end if
     else
        z= max(z, 10.d-5)                             !      [m_metals/m_tot]
        call getNPhotonsEmitted(age,dt_loc_Gyr,z,part_NpInp(j,1:nSEDgroups))
        part_NpInp(j,1:nSEDgroups) = part_NpInp(j,1:nSEDgroups)*rt_esc_frac ! escape fraction for PopII
     endif
     part_NpInp(j,:) = part_NpInp(j,:)*mass*scale_inp !#photons

     if(showSEDstats .and. nSEDgroups .gt. 0) then
        step_nPhot = step_nPhot+part_NpInp(j,1)*scale_nPhot
        step_nStar = step_nStar+dt_loc_Gyr*Gyr2sec/scale_t
        step_mStar = step_mStar+mp(ind_part(j)) * scale_msun             &
                                          * dt_loc_Gyr * Gyr2sec / scale_t
     endif

     indp(j)= ncoarse + (icell(j)-1)*ngridmax + ind_grid(ind_grid_part(j))
  end do

  do j=1,np                       ! Update hydro variables due to feedback
     do ip=1,nSEDgroups
        rtunew(indp(j),iGroups(ip)) &
             = rtunew(indp(j),iGroups(ip)) + part_NpInp(j,ip)

     end do
  end do

END SUBROUTINE star_RT_vsweep_pp
#endif

!*************************************************************************
SUBROUTINE get_pop3_Nphoton(m, nphoton)
! compute the number of photons per s
! a bit hard coded
   use amr_commons, only:dp
   use rt_parameters, only:nGroups,group_egy
   implicit none
   real(dp),intent(in)::m
   real(dp),dimension(1:nGroups)::nphoton
   real(dp)::x,y,ev_avg,frac
   integer ::ip,j,ndeg,nbin_HII
   integer ::ndeg_HII=3,ndeg_HeII=5,ndeg_HeIII=3,ndeg_H2=3
   real(dp),dimension(0:3)::a_HII,a_HeIII,a_H2
   real(dp),dimension(0:5)::a_HeII,a_fit
   real(dp):: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun

   ! Conversion factor from user units to cgs units
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   scale_msun = scale_l**3*scale_d/1.989d33

   ! fit to Schaerer (2002) Table 4,5
   a_HII = (/3.7931090e+01,1.4585544e+01,-6.0434098e+00,8.9166773e-01/)
   a_HeII = (/-1.4045167e+01,1.6406353e+02,-1.7101901e+02,8.8029022e+01,-2.2113640e+01,2.1670952e+00/)
   a_HeIII = (/2.3361417e+01,2.5165704e+01,-7.9901436e+00,8.7320721e-01/)
   a_H2 = (/4.3037504e+01,6.7958632e+00,-2.2352535e+00,2.9633043e-01/)

   x = log10(m*scale_msun)
   ! initialise
   nphoton(:) = 1.d-20 !Harley edit to prevent 0 photon density in group
   nbin_HII = 0

   ! see if there are two bins for HII (i.e. 13.6-15.2, 15.2-24.59)
   do ip=1,nGroups
      ev_avg = group_egy(ip)
      if(ev_avg.gt.13.60.and.ev_avg.lt.24.59) nbin_HII = nbin_HII + 1
   end do

   do ip=1,nGroups
      ev_avg = group_egy(ip)

      frac = 1d0
      if      (ev_avg.gt.11.20.and.ev_avg.lt.13.60) then
         ndeg = ndeg_H2
         a_fit(0:ndeg)=a_H2
      else if (ev_avg.gt.13.60.and.ev_avg.lt.24.59) then
         ndeg = ndeg_HII
         a_fit(0:ndeg)=a_HII
         if(nbin_HII.eq.2)then ! (i.e. 13.6-15.2, 15.2-24.59)
            if(ev_avg.lt.15.20)then ! for T=1d5K black body spectrum
               frac = 0.1
            else
               frac = 0.9
            endif
         endif
      else if (ev_avg.gt.24.59.and.ev_avg.lt.54.42) then
         ndeg = ndeg_HeII
         a_fit(0:ndeg)=a_HeII
      else if (ev_avg.gt.54.42) then
         ndeg = ndeg_HeIII
         a_fit(0:ndeg)=a_HeIII
      else ! neglect any optical to IR photons, bc popIII will appear only when gas is nearly pristine
         ndeg = 0
         a_fit(:)=0d0
      endif

      if(ndeg.eq.0) then
         nphoton(ip) = 0d0
      else
         y = 0d0
         do j=0,ndeg
            y = y + x**dble(j)*a_fit(j)
         end do
         nphoton(ip) = 10d0**y*frac
      endif

   end do
   ! what we want is #/s per 1Msun!
   nphoton(:) = nphoton(:) / (m*scale_msun)

END SUBROUTINE get_pop3_Nphoton

!*************************************************************************
SUBROUTINE get_pop3_props (Tk,egy_star,csn_star,cse_star)
! compute the average energy bins, photoionisation cross-sections,
   use amr_commons, only:dp
   use rt_parameters, only:nGroups,nIons,group_egy,groupL0,groupL1
   implicit none
   real(dp),dimension(1:nGroups,1:nIons)::csn_star,cse_star
   real(dp),dimension(1:nGroups)::egy_star
   real(dp)::Tk,pL0,pL1,lam_cm
   integer ::nLs=2000,i,ip,ii
   real(dp),allocatable,dimension(:)::Ls,fakeSED
   real(dp),dimension(1:nGroups,1:nIons),save::csn_star_pop3,cse_star_pop3
   real(dp),dimension(1:nGroups),save::egy_star_pop3
   real(dp)::h_pl=6.6260755d-27,kB=1.3806504d-16,c_cms=2.999d10
   logical,save:: initialise=.true.
   real(dp)::top,bot
   if(initialise)then
      csn_star(:,:) = csn_star_pop3(:,:)
      cse_star(:,:) = cse_star_pop3(:,:)
      egy_star(:)   = egy_star_pop3(:)
      allocate(Ls(1:nLs))
      allocate(fakeSED(1:nLs))
      do i=1,nLs
         Ls(i) = 10d0**(dble(i)/nLs*(5d0-0d0))  ! 1A - 10000A
         lam_cm = Ls(i)/1d8 ! A->cm
         !fakeSED(i) = 2d0*h_pl*c_cms**2/lam_cm**5d0/(exp(h_pl*c_cms/lam_cm/kB/Tk)-1d0)
         ! To prevent floating point overflow
         top = 10.d0**(LOG10(2.d0*h_pl*c_cms**2) - 5.d0*LOG10(lam_cm))
         bot = MIN(((((h_pl*c_cms)/lam_cm)/kb)/Tk),100.d0)
         fakeSED(i) = top / (EXP(bot)-1d0)
      end do

      do ip=1,nGroups
         pL0 = groupL0(ip)
         pL1 = groupL1(ip)
         egy_star_pop3(ip) = getSEDEgy(Ls,fakeSED,nLs,pL0,pL1)

         do ii=1,nIons
            csn_star_pop3(ip,ii) = getSEDcsn(Ls,fakeSED,nLs,pL0,pL1,ii)
            cse_star_pop3(ip,ii) = getSEDcse(Ls,fakeSED,nLs,pL0,pL1,ii)
         end do
      end do


      deallocate(Ls,fakeSED)
   endif

   csn_star(:,:) = csn_star_pop3(:,:)
   cse_star(:,:) = cse_star_pop3(:,:)
   egy_star(:)   = egy_star_pop3(:)

   initialise = .false.
END SUBROUTINE get_pop3_props

!*************************************************************************
SUBROUTINE getNPhotonsEmitted_pop3(mass,age_Gyr,dt_loc_Gyr,nphoton)
   use amr_commons, only:dp
   implicit none
   real(dp),intent(in)::mass,age_Gyr,dt_loc_Gyr
   real(dp)::lifetime_Gyr
   real(dp),dimension(1:nSEDgroups)::nphoton

   call get_pop3_ageGyr(mass,lifetime_Gyr)

   nphoton(:)=0d0
   if (age_Gyr.le.lifetime_Gyr) then
      call get_pop3_Nphoton(mass, nphoton) ! [nphoton] = #/s/Msun
      nphoton(:) = nphoton(:)*dt_loc_Gyr*1d9*365d0*24d0*3600d0
   endif

END SUBROUTINE getNPhotonsEmitted_pop3

END MODULE SED_module


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for UV table of redshift dependent photon fluxes, cross sections
! and photon energies, per photon group. This is mainly useful for
! calculating heating rates due to the homogeneous UV background in the
! cooling-module, and for  radiative transfer of the UV background, where
! we then assign intensities, average cross sections and energies to UV
! photon groups.
!
! There are two tables of UV properties:
! UV_rates_table: Integrated ionization and heating rates for a
!                 homogeneous UV background.
! UV_groups_table: UV properties integrated over each photon group
!                 wavelength interval: photon flux, average energy,
!                 average cross section, energy weighted cross section.
!_________________________________________________________________________
!
MODULE UV_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  use spectrum_integrator_module
  use rt_parameters,only:c_cgs, eV_to_erg, hp, nIons, ionEVs, nGroups

  implicit none

  PUBLIC nUVgroups, iUVvars_cool, UV_Nphot_cgs, init_UV_background       &
       , inp_UV_rates_table, inp_UV_groups_table, UV_minz, UV_maxz       &
       , update_UVsrc, iUVgroups, init_UV_background_RTZ                 &
       , update_UV_background_RTZ

  PRIVATE   ! default

  ! UV properties for different redshifts---------------------------------
  integer::UV_nz
  real(dp),allocatable,dimension(:)::UV_zeds
  ! Redshift interval
  real(dp)::UV_minz, UV_maxz
  ! Table of integrated ionization and heating rates per ion species
  ! (UV_nz, nIons, 2):(z, nIon, ionization rate [# s-1]: Hrate [erg s-1])
  real(dp),allocatable,dimension(:,:,:)::UV_rates_table
  ! rt_n_UVsrc vectors of redshift dependent fluxes for UV background
  real(dp),allocatable::UV_fluxes_cgs(:)    !                    [#/cm2/s]
  real(dp),allocatable::UV_Nphot_cgs(:)     !        Photon density [cm-3]
  integer::nUVgroups=0                      !        # of UV photon groups
  ! UV group indexes among nGroup groups,
  ! UV Np indexes among solve_cooling vars:
  integer,allocatable::iUVgroups(:),iUVvars_cool(:)
  ! Table of photon group props (UV_nz, nUVgroups, 1+2*nIons):
  !(z, group, flux [#/cm2/s]: nIons*csn [cm-2]: nIons*egy [ev])
  real(dp),allocatable,dimension(:,:,:)::UV_groups_table
  ! The tables do not have constant redshift intervals, so we need to
  ! first locate the correct interval when doing interpolation.
  real(dp),parameter::fourPi=12.566371
  ! ----------------------------------------------------------------------
  ! UV rates tables for RTZ
  real(dp), dimension(60)    :: UVBG_Redshifts
  real(dp), dimension(60,1)  :: UVBG_PI_Hydrogen,  UVBG_PH_Hydrogen
  real(dp), dimension(60,2)  :: UVBG_PI_Helium,    UVBG_PH_Helium
  real(dp), dimension(60,6)  :: UVBG_PI_Carbon,    UVBG_PH_Carbon
  real(dp), dimension(60,7)  :: UVBG_PI_Nitrogen,  UVBG_PH_Nitrogen
  real(dp), dimension(60,8)  :: UVBG_PI_Oxygen,    UVBG_PH_Oxygen
  real(dp), dimension(60,10) :: UVBG_PI_Neon,      UVBG_PH_Neon
  real(dp), dimension(60,12) :: UVBG_PI_Magnesium, UVBG_PH_Magnesium
  real(dp), dimension(60,14) :: UVBG_PI_Silicon,   UVBG_PH_Silicon
  real(dp), dimension(60,16) :: UVBG_PI_Sulfur,    UVBG_PH_Sulfur
  real(dp), dimension(60,26) :: UVBG_PI_Iron,      UVBG_PH_Iron

CONTAINS
!*************************************************************************
SUBROUTINE init_UV_background_RTZ()
   ! Initialize all of the UV background tables for H, He, and metals
   use amr_commons,only:myid
   use rt_parameters
   implicit none

   integer:: i
   integer:: n_redshifts = 60 ! Number of redshifts in the haardt madau data
   logical::ok_PH,ok_PI
   character(len=128)::zzz
   character(len=128)::PI_H,  PH_H
   character(len=128)::PI_He, PH_He
   character(len=128)::PI_C,  PH_C
   character(len=128)::PI_N,  PH_N
   character(len=128)::PI_O,  PH_O
   character(len=128)::PI_Ne, PH_Ne
   character(len=128)::PI_Mg, PH_Mg
   character(len=128)::PI_Si, PH_Si
   character(len=128)::PI_S,  PH_S
   character(len=128)::PI_Fe, PH_Fe

  if(myid==1) write(*,*) 'Initializing RTZ UV background data'

  write(zzz,  '(a,a)')   trim(uvbg_rtz_dir),"/redshifts.dat"
  write(PI_H, '(a,a)')   trim(uvbg_rtz_dir),"/hydrogen_pi.dat"
  write(PH_H, '(a,a)')   trim(uvbg_rtz_dir),"/hydrogen_ph.dat"
  write(PI_He,'(a,a)')   trim(uvbg_rtz_dir),"/helium_pi.dat"
  write(PH_He,'(a,a)')   trim(uvbg_rtz_dir),"/helium_ph.dat"
  write(PI_C, '(a,a)')   trim(uvbg_rtz_dir),"/carbon_pi.dat"
  write(PH_C, '(a,a)')   trim(uvbg_rtz_dir),"/carbon_ph.dat"
  write(PI_N, '(a,a)')   trim(uvbg_rtz_dir),"/nitrogen_pi.dat"
  write(PH_N, '(a,a)')   trim(uvbg_rtz_dir),"/nitrogen_ph.dat"
  write(PI_O, '(a,a)')   trim(uvbg_rtz_dir),"/oxygen_pi.dat"
  write(PH_O, '(a,a)')   trim(uvbg_rtz_dir),"/oxygen_ph.dat"
  write(PI_Ne,'(a,a)')   trim(uvbg_rtz_dir),"/neon_pi.dat"
  write(PH_Ne,'(a,a)')   trim(uvbg_rtz_dir),"/neon_ph.dat"
  write(PI_Mg,'(a,a)')   trim(uvbg_rtz_dir),"/magnesium_pi.dat"
  write(PH_Mg,'(a,a)')   trim(uvbg_rtz_dir),"/magnesium_ph.dat"
  write(PI_Si,'(a,a)')   trim(uvbg_rtz_dir),"/silicon_pi.dat"
  write(PH_Si,'(a,a)')   trim(uvbg_rtz_dir),"/silicon_ph.dat"
  write(PI_S, '(a,a)')   trim(uvbg_rtz_dir),"/sulfur_pi.dat"
  write(PH_S, '(a,a)')   trim(uvbg_rtz_dir),"/sulfur_ph.dat"
  write(PI_Fe,'(a,a)')   trim(uvbg_rtz_dir),"/iron_pi.dat"
  write(PH_Fe,'(a,a)')   trim(uvbg_rtz_dir),"/iron_ph.dat"

  inquire(file=zzz, exist=ok_PI)
  if(.not. ok_PI) then
     if(myid.eq.1) then
        write(*,*) 'Cannot read RTZ UV backround redshift files...'
     endif
     call clean_stop
  end if
  ! Load in the redshifd data
  open(unit=15, file=zzz,status='old',form='formatted')
  do i = 1,n_redshifts
   read (15,*) UVBG_Redshifts(i)
  end do
  close(15)

  inquire(file=PI_H, exist=ok_PI)
  inquire(file=PH_H, exist=ok_PH)
  if(.not. ok_PI .or. .not. ok_PH ) then
     if(myid.eq.1) then
        write(*,*) 'Cannot read RTZ UV backround Hydrogen files...'
     endif
     call clean_stop
  end if
  ! Load in the hydrogen data
  open(unit=15, file=PI_H,status='old',form='formatted')
  do i = 1,n_redshifts
   read (15,*) UVBG_PI_Hydrogen(i,:)
  end do
  close(15)
  open(unit=15, file=PH_H,status='old',form='formatted')
  do i = 1,n_redshifts
   read (15,*) UVBG_PH_Hydrogen(i,:)
  end do
  close(15)


  inquire(file=PI_He, exist=ok_PI)
  inquire(file=PH_He, exist=ok_PH)
  if(.not. ok_PI .or. .not. ok_PH ) then
     if(myid.eq.1) then
        write(*,*) 'Cannot read RTZ UV backround Helium files...'
     endif
     call clean_stop
  end if
  ! Load in the helium data
  open(unit=15, file=PI_He,status='old',form='formatted')
  do i = 1,n_redshifts
   read (15,*) UVBG_PI_Helium(i,:)
  end do
  close(15)
  open(unit=15, file=PH_He,status='old',form='formatted')
  do i = 1,n_redshifts
   read (15,*) UVBG_PH_Helium(i,:)
  end do
  close(15)


  if (carbon_ions) then
   inquire(file=PI_C, exist=ok_PI)
   inquire(file=PH_C, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Carbon files...'
      endif
      call clean_stop
   end if
   ! Load in the carbon data
   open(unit=15, file=PI_C,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Carbon(i,:)
   end do
   close(15)
   open(unit=15, file=PH_C,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Carbon(i,:)
   end do
   close(15)
  end if

  if (nitrogen_ions) then
   inquire(file=PI_N, exist=ok_PI)
   inquire(file=PH_N, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Nitrogen files...'
      endif
      call clean_stop
   end if
   ! Load in the nitrogen data
   open(unit=15, file=PI_N,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Nitrogen(i,:)
   end do
   close(15)
   open(unit=15, file=PH_N,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Nitrogen(i,:)
   end do
   close(15)
  end if

  if (oxygen_ions) then
   inquire(file=PI_O, exist=ok_PI)
   inquire(file=PH_O, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Oxygen files...'
      endif
      call clean_stop
   end if
   ! Load in the oxygen data
   open(unit=15, file=PI_O,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Oxygen(i,:)
   end do
   close(15)
   open(unit=15, file=PH_O,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Oxygen(i,:)
   end do
   close(15)
  end if

  if (neon_ions) then
   inquire(file=PI_Ne, exist=ok_PI)
   inquire(file=PH_Ne, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Neon files...'
      endif
      call clean_stop
   end if
   ! Load in the neon data
   open(unit=15, file=PI_Ne,status='old',form='formatted')
   do i = 1,n_redshifts
   read (15,*) UVBG_PI_Neon(i,:)
   end do
   close(15)
   open(unit=15, file=PH_Ne,status='old',form='formatted')
   do i = 1,n_redshifts
   read (15,*) UVBG_PH_Neon(i,:)
   end do
   close(15)
  end if

  if (magnesium_ions) then
   inquire(file=PI_Mg, exist=ok_PI)
   inquire(file=PH_Mg, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Magnesium files...'
      endif
      call clean_stop
   end if
   ! Load in the magnesium data
   open(unit=15, file=PI_Mg,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Magnesium(i,:)
   end do
   close(15)
   open(unit=15, file=PH_Mg,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Magnesium(i,:)
   end do
   close(15)
  end if

  if (silicon_ions) then
   inquire(file=PI_Si, exist=ok_PI)
   inquire(file=PH_Si, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Silicon files...'
      endif
      call clean_stop
   end if
   ! Load in the silicon data
   open(unit=15, file=PI_Si,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Silicon(i,:)
   end do
   close(15)
   open(unit=15, file=PH_Si,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Silicon(i,:)
   end do
   close(15)
  end if

  if (sulfur_ions) then
   inquire(file=PI_S, exist=ok_PI)
   inquire(file=PH_S, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Sulfur files...'
      endif
      call clean_stop
   end if
   ! Load in the sulfur data
   open(unit=15, file=PI_S,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Sulfur(i,:)
   end do
   close(15)
   open(unit=15, file=PH_S,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Sulfur(i,:)
   end do
   close(15)
  end if

  if (iron_ions) then
   inquire(file=PI_Fe, exist=ok_PI)
   inquire(file=PH_Fe, exist=ok_PH)
   if(.not. ok_PI .or. .not. ok_PH ) then
      if(myid.eq.1) then
         write(*,*) 'Cannot read RTZ UV backround Hydrogen files...'
      endif
      call clean_stop
   end if
   ! Load in the iron data
   open(unit=15, file=PI_Fe,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PI_Iron(i,:)
   end do
   close(15)
   open(unit=15, file=PH_Fe,status='old',form='formatted')
   do i = 1,n_redshifts
      read (15,*) UVBG_PH_Iron(i,:)
   end do
   close(15)
  end if

END SUBROUTINE init_UV_background_RTZ
!*************************************************************************
SUBROUTINE update_UV_background_RTZ(z)
! Compute UV heating and ionization rates by interpolation from table.
! z     => Redshift
!-------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  implicit none
  real(dp), intent(in):: z
  real(dp):: ret(nIons+1,2), dz0, dz1, uvbg_sf
  integer:: iz0, iz1
!-------------------------------------------------------------------------
  ret=0. ; if(z .gt. UVBG_Redshifts(60)) RETURN
  call inp_1d(UVBG_Redshifts, 60, z, iz0, iz1, dz0, dz1)

  if (uvbg_suppression_redshift.gt.0.d0.and.z.gt.uvbg_suppression_redshift) then
     ! The -5 factor ramps up the background quite quickly
     ! Can change this to make a less steep slope
     uvbg_sf = EXP(5.d0 * (uvbg_suppression_redshift - z))

     ! Scale down the the uvbg by a suppression factor
     dz0 = dz0 * uvbg_sf
     dz1 = dz1 * uvbg_sf
  else
     uvbg_sf=1.d0
  end if

  ! Hydrogen
  UV_background_hydrogen         = dz0*UVBG_PI_Hydrogen(iz1, 1) + dz1*UVBG_PI_Hydrogen(iz0, 1)
  UV_background_hydrogen_heating = dz0*UVBG_PH_Hydrogen(iz1, 1) + dz1*UVBG_PH_Hydrogen(iz0, 1)

  ! Helium
  UV_background_helium         = dz0*UVBG_PI_Helium(iz1, :) + dz1*UVBG_PI_Helium(iz0, :)
  UV_background_helium_heating = dz0*UVBG_PH_Helium(iz1, :) + dz1*UVBG_PH_Helium(iz0, :)

  ! Carbon
  if (carbon_ions) then
     UV_background_carbon         = dz0*UVBG_PI_Carbon(iz1, 1:n_carbon_ions) + dz1*UVBG_PI_Carbon(iz0, 1:n_carbon_ions)
     UV_background_carbon_heating = dz0*UVBG_PH_Carbon(iz1, 1:n_carbon_ions) + dz1*UVBG_PH_Carbon(iz0, 1:n_carbon_ions)
  end if

  ! Nitrogen
  if (nitrogen_ions) then
     UV_background_nitrogen         = dz0*UVBG_PI_Nitrogen(iz1, 1:n_nitrogen_ions) + dz1*UVBG_PI_Nitrogen(iz0, 1:n_nitrogen_ions)
     UV_background_nitrogen_heating = dz0*UVBG_PH_Nitrogen(iz1, 1:n_nitrogen_ions) + dz1*UVBG_PH_Nitrogen(iz0, 1:n_nitrogen_ions)
  end if

  ! Oxygen
  if (oxygen_ions) then
     UV_background_oxygen         = dz0*UVBG_PI_Oxygen(iz1, 1:n_oxygen_ions) + dz1*UVBG_PI_Oxygen(iz0, 1:n_oxygen_ions)
     UV_background_oxygen_heating = dz0*UVBG_PH_Oxygen(iz1, 1:n_oxygen_ions) + dz1*UVBG_PH_Oxygen(iz0, 1:n_oxygen_ions)
  end if

  ! Neon
  if (neon_ions) then
     UV_background_neon         = dz0*UVBG_PI_Neon(iz1, 1:n_neon_ions) + dz1*UVBG_PI_Neon(iz0, 1:n_neon_ions)
     UV_background_neon_heating = dz0*UVBG_PH_Neon(iz1, 1:n_neon_ions) + dz1*UVBG_PH_Neon(iz0, 1:n_neon_ions)
  end if

  ! Magnesium
  if (magnesium_ions) then
     UV_background_magnesium         = dz0*UVBG_PI_Magnesium(iz1, 1:n_magnesium_ions) + dz1*UVBG_PI_Magnesium(iz0, 1:n_magnesium_ions)
     UV_background_magnesium_heating = dz0*UVBG_PH_Magnesium(iz1, 1:n_magnesium_ions) + dz1*UVBG_PH_Magnesium(iz0, 1:n_magnesium_ions)
  end if

  ! Silicon
  if (silicon_ions) then
     UV_background_silicon         = dz0*UVBG_PI_Silicon(iz1, 1:n_silicon_ions) + dz1*UVBG_PI_Silicon(iz0, 1:n_silicon_ions)
     UV_background_silicon_heating = dz0*UVBG_PH_Silicon(iz1, 1:n_silicon_ions) + dz1*UVBG_PH_Silicon(iz0, 1:n_silicon_ions)
  end if

  ! Sulfur
  if (sulfur_ions) then
     UV_background_sulfur         = dz0*UVBG_PI_Sulfur(iz1, 1:n_sulfur_ions) + dz1*UVBG_PI_Sulfur(iz0, 1:n_sulfur_ions)
     UV_background_sulfur_heating = dz0*UVBG_PH_Sulfur(iz1, 1:n_sulfur_ions) + dz1*UVBG_PH_Sulfur(iz0, 1:n_sulfur_ions)
  end if

  ! Iron
  if (iron_ions) then
     UV_background_iron         = dz0*UVBG_PI_Iron(iz1, 1:n_iron_ions) + dz1*UVBG_PI_Iron(iz0, 1:n_iron_ions)
     UV_background_iron_heating = dz0*UVBG_PH_Iron(iz1, 1:n_iron_ions) + dz1*UVBG_PH_Iron(iz0, 1:n_iron_ions)
  end if

  ! Turn off the UVBG in sub-ionizing bins
  if (no_subionizing_ubvg) then
     ! Carbon
     UV_background_carbon(1) = 1.d-50
     UV_background_carbon_heating(1) = 1.d-50
     UV_background_carbon(2) = 1.d-50
     UV_background_carbon_heating(2) = 1.d-50
     ! Magnesium
     UV_background_magnesium(1) = 1.d-50
     UV_background_magnesium_heating(1) = 1.d-50
     UV_background_magnesium(2) = 1.d-50
     UV_background_magnesium_heating(2) = 1.d-50
     ! Silicon
     UV_background_silicon(1) = 1.d-50
     UV_background_silicon_heating(1) = 1.d-50
     UV_background_silicon(2) = 1.d-50
     UV_background_silicon_heating(2) = 1.d-50
     ! Sulfur
     UV_background_sulfur(1) = 1.d-50
     UV_background_sulfur_heating(1) = 1.d-50
     UV_background_sulfur(2) = 1.d-50
     UV_background_sulfur_heating(2) = 1.d-50
     ! Iron
     UV_background_iron(1) = 1.d-50
     UV_background_iron_heating(1) = 1.d-50
     UV_background_iron(2) = 1.d-50
     UV_background_iron_heating(2) = 1.d-50
  end if
END SUBROUTINE update_UV_background_RTZ

!*************************************************************************
SUBROUTINE init_UV_background()
! Initiate UV props table, which gives photon fluxes and average photon
! cross sections and energies of each photon group as a function of z.
! The data is read from a file specified by uv_file containing:
! a) number of redshifts
! b) redshifts (increasing order)
! c) number of wavelengths
! d) wavelengths (increasing order) [Angstrom]
! e) fluxes per (redshift,wavelength) [photons cm-2 s-1 A-1 sr-1]
!-------------------------------------------------------------------------
  use amr_commons,only:myid,IOGROUPSIZE,ncpu
  use rt_parameters
  use SED_module
  use mpi_mod
  integer:: nLs                                 ! # of bins of wavelength
  real(kind=8),allocatable  :: Ls(:)            ! Wavelengths
  real(kind=8),allocatable  :: UV(:,:)          ! UV f(lambda,z)
  real(kind=8),allocatable  :: tbl(:,:), tbl2(:,:)
  integer::i,ia,iz,ip,ii,dum,locid,ncpu2,ierr
  logical::ok
  real(kind=8)::da, dz, pL0,pL1
  integer,parameter::tag=1133
  integer::dummy_io,info2

!-------------------------------------------------------------------------
  ! First check if there is any need for UV setup:
  if(rt_UVsrc_nHmax .le. 0d0 .and. .not. haardt_madau) return

  if(myid==1) print*,'Initializing UV background'

  ! Read UV spectra from file---------------------------------------------
  if(UV_FILE .eq. '') &
       call get_environment_variable('RAMSES_UV_FILE', UV_FILE)
  inquire(file=TRIM(uv_file), exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then
        write(*,*)'Cannot access UV file ',TRIM(uv_file)
        write(*,*)'File '//TRIM(uv_file)//' not found'
        write(*,*)'You need to set the RAMSES_UV_FILE envvar' // &
                  ' to the correct path, or use the namelist var uv_file'
     endif
     call clean_stop
  end if

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  ! Read redshifts, wavelengths and spectra
  open(unit=10,file=TRIM(uv_file),status='old',form='unformatted')
  read(10) UV_nz ; allocate(UV_zeds(UV_nz)) ; read(10) UV_zeds(:)
  read(10) nLs   ; allocate(Ls(nLs))        ; read(10) Ls(:)
  allocate(UV(nLs,UV_nz))                   ; read(10) UV(:,:)
  close(10)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif


  ! If mpi then share the UV integration between the cpus:
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,locid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu2,ierr)
#endif
#ifdef WITHOUTMPI
  locid=0 ; ncpu2=1
#endif

  ! Shift the highest z in the table (10) to reionization epoch,
  ! so that we start injecting at z_reion
  if(z_reion .gt. UV_zeds(UV_nz-1))  UV_zeds(UV_nz) = z_reion
  UV_minz = UV_zeds(1) ; UV_maxz=UV_zeds(UV_nz)

  ! Non-propagated UV background -----------------------------------------
  if(haardt_madau) then
     if(myid==1) print*,'The UV background is homogeneous'
     allocate(UV_rates_table(UV_nz, nIons+1, 2))
     allocate(tbl(UV_nz, 2))
     do ii = 1, nIons + 1
        tbl=0.
        do iz = locid+1,UV_nz,ncpu2
           tbl(iz,1)= getUV_Irate(Ls,UV(:,iz),nLs,ii)
           tbl(iz,2)= getUV_Hrate(Ls,UV(:,iz),nLs,ii)
        end do
#ifndef WITHOUTMPI
        allocate(tbl2(UV_nz,2))
        call MPI_ALLREDUCE(tbl, tbl2, UV_nz*2,  &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        tbl = tbl2
        deallocate(tbl2)
#endif
        UV_rates_table(:,ii,:)=tbl
     end do
     deallocate(tbl)
     if (myid==1) call write_UVrates_table
  endif

  ! Note: Need to take special care of the highest redshift. Some UV     !
  !       spectra (read: Faucher-Giguere) have zero fluxes at the highest!
  !       redshift, and integrations of energies and cross sections of   !
  !       these fluxes results in a division by zero and NaNs.           !
  !       The approach taken is then to keep zero flux at the highest    !
  !       redshift but otherwise use the same values of energies and     !
  !       cross sections as the next-highest redshift.                   !

  ! Propagated UV background----------------------------------------------
  if(rt_UVsrc_nHmax .gt. 0.d0) then ! UV propagation from diffuse cells--
     if(myid==1) print*,'The UV background is propagated'
     if(myid==1 .and. haardt_madau) then
          print*,'ATT: UV background is BOTH homogeneous and propagated'
          print*,'  You likely don''t want this duplicated background...'
       endif
     rt_isDiffuseUVsrc=.true.
     if(nUVgroups .eq. 0) nUVgroups = nGroups ! Default: All groups are UV
     nSEDgroups=nGroups-nUVgroups
     ! SED groups are the first nSEDgroups, UV groups the last nUVgroups
     allocate(iUVgroups(nUVgroups)) ; allocate(iUVvars_cool(nUVgroups))
     do i=1,nUVgroups                 !      Initialize UV group indexes
        iUVgroups(i) = nGroups-nUVgroups+i        ! UV groups among groups
        iUVvars_cool(i) = 4+iUVgroups(i) !UV Np's in vars in solve_cooling
     end do
     if(nUVgroups .gt. 0) then
        allocate(UV_fluxes_cgs(nUVgroups))             ; UV_fluxes_cgs=0.
        allocate(UV_Nphot_cgs(nUVgroups))              ; UV_Nphot_cgs=0.
     endif

     ! Initialize photon groups table-------------------------------------
     allocate(UV_groups_table(UV_nz, nUVgroups, 2+2*nIons))
     allocate(tbl(UV_nz, 2+2*nIons))
     do ip = 1,nUVgroups           !                    Loop photon groups
        tbl=0.
        pL0 = groupL0(nSEDgroups+ip) !  Energy interval of photon group ip
        pL1 = groupL1(nSEDgroups+ip) !
        do iz = locid+1,UV_nz,ncpu2
           tbl(iz,1) =        getUVFlux(Ls,UV(:,iz),nLs,pL0,pL1)
           if(tbl(iz,1) .eq. 0.d0) cycle     ! Can't integrate zero fluxes
           tbl(iz,2) =        getUVEgy(Ls,UV(:,iz),nLs,pL0,pL1)
           do ii = 1,nIons
              tbl(iz,1+ii*2)= getUVcsn( Ls,UV(:,iz),nLs,pL0,pL1,ii)
              tbl(iz,2+ii*2)= getUVcse( Ls,UV(:,iz),nLs,pL0,pL1,ii)
           end do
        end do
#ifndef WITHOUTMPI
        allocate(tbl2(UV_nz,2+2*nIons))
        call MPI_ALLREDUCE(tbl,tbl2,UV_nz*(2+2*nIons),&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        tbl = tbl2
        deallocate(tbl2)
#endif
        if(tbl(UV_nz,1) .eq. 0.d0) &            !                Zero flux
             tbl(UV_nz,2:)=tbl(UV_nz-1,2:)
        UV_groups_table(:,ip,:)=tbl
     end do
     deallocate(tbl) ; deallocate(Ls) ; deallocate(UV)

     call update_UVsrc
     if (myid==1) call write_UVgroups_tables
  endif ! End propagated UV background

END SUBROUTINE init_UV_background

!*************************************************************************
SUBROUTINE inp_UV_rates_table(z, ret, z_damp)
! Compute UV heating and ionization rates by interpolation from table.
! z     => Redshift
! ret  <=  [nIons,2] interpolated values of all UV rates.
!          1=ionization rate [# s-1], 2=heating rate [erg s-1]
! z_damp=> Optional. If set and true, the UV rates are expoenentially
!          damped according to the difference z-z_reion (i.e. strong
!          damping if z > z_reion, and a gradual decrease in the damping
!          when z~z_reion).
!-------------------------------------------------------------------------
  use amr_parameters
  real(dp), intent(in):: z
  real(dp):: ret(nIons+1,2), dz0, dz1, zr_factor
  integer:: iz0, iz1
  logical, optional, intent (in) :: z_damp
!-------------------------------------------------------------------------
  ret=0. ; if(z .gt. UV_maxz) RETURN
  call inp_1d(UV_zeds, UV_nz, z, iz0, iz1, dz0, dz1)
  ret = dz0*UV_rates_table(iz1, :, :) + dz1*UV_rates_table(iz0, :, :)
  if (present(z_damp)) then
     if (z_damp) then
        zr_factor = 20d0*(z/z_reion)**6d0
        ret = ret * 10d0**(-zr_factor )
     endif
  endif
END SUBROUTINE inp_UV_rates_table

!*************************************************************************
SUBROUTINE inp_UV_groups_table(z, ret, z_damp)
! Compute UV properties by interpolation from table.
! z    => Redshift
! ret <=  [nGroups,1+2*nIons) interpolated values of all UV properties.
!         1=ph. flux [#/cm2/s], 2*i=group_csn[cm-2], 1+2*i=group_egy [ev]
! z_damp=> Optional. If set and true, the UV rates are expoenentially
!          damped according to the difference z-z_reion (i.e. strong
!          damping if z > z_reion, and a gradual decrease in the damping
!          when z~z_reion).
!-------------------------------------------------------------------------
  use amr_parameters
  real(dp), intent(in):: z
  real(dp):: ret(nUVgroups,2+2*nIons), dz0, dz1, zr_factor
  integer:: iz0, iz1
  logical, optional, intent (in) :: z_damp
!-------------------------------------------------------------------------
  ret=0. ; if(z .gt. UV_maxz) RETURN
  call inp_1d(UV_zeds, UV_nz, z, iz0, iz1, dz0, dz1)
  ret = dz0 * UV_groups_table(iz1, :, :) + dz1 * UV_groups_table(iz0, :,:)
  if (present(z_damp)) then
     if (z_damp) then
        ! Weaker damping than in the non-RT UV case:
        zr_factor = 2d0*(z/z_reion)**6d0
        ret = ret * exp(-zr_factor )
     endif
  endif
END SUBROUTINE inp_UV_groups_table

!*************************************************************************
SUBROUTINE update_UVsrc

! Update UV background source properties. So as not to do too much of
! this, this should only be done every coarse timestep.
!-------------------------------------------------------------------------
  use rt_parameters
  use SED_module
  use amr_commons,only:t,levelmin,myid,aexp
  implicit none
  integer::ilevel,i
  real(dp),allocatable,save::UVprops(:,:) !Each group: flux, egy, csn, cse
  real(dp)::scale_Np, scale_Fp, redshift
!-------------------------------------------------------------------------
  if(.not.rt_isDiffuseUVsrc) return
  if(nUVgroups.le.0) then
     if(myid==1) write(*,*) 'No groups dedicated to the UV background!'
     RETURN
  endif
  if(.not. allocated(UVprops)) allocate(UVprops(nUVgroups,2+2*nIons))
  call rt_units(scale_Np, scale_Fp)

  redshift=1./aexp-1.

  ! Turn on RT after z=UV_maxz:
  if(redshift.le.UV_maxz .and. .not. rt_advect) then
     if(myid==1) then
        write(*,*) '*****************************************************'
        write(*,*) 'Turned on RT advection and the UV background'
        write(*,*) '*****************************************************'
     endif
     rt_advect=.true.
  endif

  if(redshift .gt. UV_maxz) return ! UV background not turned on yet

  call inp_UV_groups_table(redshift, UVprops, .true.)
  UV_fluxes_cgs(:)      = UVprops(:,1)
  UV_Nphot_cgs          = UV_fluxes_cgs/rt_c_cgs
  group_egy(iUVgroups)  = UVprops(:,2)
  do i=1,nIons
     group_csn(iUVgroups,i)  = UVprops(:,1+2*i)
     group_cse(iUVgroups,i)  = UVprops(:,2+2*i)
  enddo

  call updateRTgroups_CoolConstants

  if(myid==1) then
     write(*,*) 'Updated UV fluxes [# cm-2 s-1] to'
     write(*,900) UV_fluxes_cgs
     call write_group_props(.true.,6)
  endif

900 format (20f16.6)
901 format (20(1pe16.6))
END SUBROUTINE update_UVsrc

!*************************************************************************
! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

!*************************************************************************
FUNCTION getUV_Irate(X, Y, N, species)
! Compute and return photoionization rate per ion, given the UV background
! Y(X). Assumes X is in Angstroms and Y in [# cm-2 s-1 sr-1 A-1].
! returns: Photoionization rate in # s-1
!-------------------------------------------------------------------------
  real(kind=8):: getUV_Irate, X(N), Y(N)
  integer :: N, species
!-------------------------------------------------------------------------
  if(species.eq.(nIons+1)) then !isH2Katz
     getUV_Irate = fourPi *  &
              integrateSpectrum(X,Y,N,11.2d0,13.6d0,1,fsig)
  else
     getUV_Irate = fourPi *  &
              integrateSpectrum(X,Y,N,dble(ionEvs(species)),0d0,species,fsig)
  endif
END FUNCTION getUV_Irate

!*************************************************************************
FUNCTION getUV_Hrate(X, Y, N, species)
! Compute and return heating rate per ion, given the UV background
! Y(X). Assumes X is in Angstroms and Y in [# cm-2 s-1 sr-1 A-1].
! returns: Heating rate in erg s-1
!-------------------------------------------------------------------------
  real(kind=8):: getUV_Hrate, X(N), Y(N), e0
  integer :: N, species
  real(kind=8),parameter :: const1=fourPi*1.d8*hp*c_cgs
  real(kind=8),parameter :: const2=fourPi*eV_to_erg
!-------------------------------------------------------------------------
  if(species.eq.(nIons+1))then !isH2Katz
     !getUV_Hrate = &
     !      const1*integrateSpectrum(X,Y,N, 11.2, 13.6, species, fsigDivLambda)&
     !     -const2*ionEvs(species) *                                         &
     !             integrateSpectrum(X,Y,N, 11.2, 13.6, species, fsig)
     getUV_Hrate = 0.0d0
  else
     e0=ionEvs(species)
     getUV_Hrate = &
           const1*integrateSpectrum(X,Y,N, e0, 0.d0, species, fsigDivLambda)&
          -const2*ionEvs(species) *                                         &
                  integrateSpectrum(X,Y,N, e0, 0.d0, species, fsig)
  endif
END FUNCTION getUV_Hrate

!*************************************************************************
FUNCTION getUVFlux(X, Y, N, e0, e1)
! Compute and return UV photon flux in energy interval (e0,e1) [eV]
! in UV spectrum Y(X). Assumes X is in [A] and Y in # cm-2 s-1 sr-1 A-1.
! returns: Photon flux in # cm-2 s-1
!-------------------------------------------------------------------------
  real(kind=8):: getUVflux, X(N), Y(N), e0, e1
  integer :: N, species
!-------------------------------------------------------------------------
  species          = 1                   ! irrelevant but must be included
  getUVflux = fourPi*integrateSpectrum(X, Y, N, e0, e1, species, f1)
END FUNCTION getUVflux

!*************************************************************************
FUNCTION getUVEgy(X, Y, N, e0, e1)
! Compute average photon energy, in eV, in energy interval (e0,e1) [eV] in
! UV spectrum Y(X). Assumes X is in [A] and Y is # cm-2 s-1 sr-1 A-1.
!-------------------------------------------------------------------------
  real(dp):: getUVEgy, X(N), Y(N), e0, e1, norm
  integer :: N,species
  real(dp),parameter :: const=1.d8*hp*c_cgs/eV_to_erg    ! unit conversion
!-------------------------------------------------------------------------
  species      = 1                       ! irrelevant but must be included
  norm         = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getUVEgy  = const * &
            integrateSpectrum(X, Y, N, e0, e1, species, fdivLambda) / norm
END FUNCTION getUVEgy

!*************************************************************************
FUNCTION getUVcsn(X, Y, N, e0, e1, species)
! Compute and return average photoionization cross-section [cm2] for given
! energy interval (e0,e1) [eV] in UV spectrum Y. Assumes X is in Angstroms
! and and Y is # cm-2 s-1 sr-1 A-1.
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  real(kind=8):: getUVcsn, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getUVcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getUVcsn = integrateSpectrum(X, Y, N, e0, e1, species, fSig)/norm
END FUNCTION getUVcsn

!************************************************************************
FUNCTION getUVcse(X, Y, N, e0, e1, species)
! Compute average energy weighted photoionization cross-section [cm2] for
! given energy interval (e0,e1) [eV] in UV spectrum Y. Assumes X is in
! Angstroms and that Y is energy intensity per angstrom.
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  real(dp):: getUVcse, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getUVcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, fdivLambda)
  getUVcse = integrateSpectrum(X, Y, N, e0, e1, species, fSigdivLambda)  &
           / norm
END FUNCTION getUVcse

!*************************************************************************
SUBROUTINE write_UVrates_table()
! Write the UV rates to a file (this is just in
! debugging, to check if the UV spectra are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::i, j
!-------------------------------------------------------------------------
  write(filename,'(A, I1, A)') 'UVrates.list'
  open(10, file=filename, status='unknown')
  write(10,*) UV_nz

  do i = 1,UV_nz
     write(10,900) UV_zeds(i), UV_rates_table(i,:,:)
  end do
  close(10)
900 format (f21.6, 20(1pe21.6))
END SUBROUTINE write_UVrates_table

!*************************************************************************
SUBROUTINE write_UVgroups_tables()

! Write the UV photon group properties to files (this is just in
! debugging, to check if the UV spectra are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j
!-------------------------------------------------------------------------
  do ip=1,nUVgroups
     write(filename,'(A, I1, A)') 'UVtable', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) UV_nz

     do i = 1,UV_nz
        write(10,901)                                                    &
                UV_zeds(i)           ,                                   &
                UV_groups_table(i,ip,1), UV_groups_table(i,ip,2),        &
                UV_groups_table(i,ip,3), UV_groups_table(i,ip,4),        &
                UV_groups_table(i,ip,5), UV_groups_table(i,ip,6),        &
                UV_groups_table(i,ip,7), UV_groups_table(i,ip,8)
     end do
     close(10)
  end do
901 format (f21.6,   f21.6,   f21.6,   1pe21.6, &
          & 1pe21.6, 1pe21.6, 1pe21.6, 1pe21.6, 1pe21.6  )
END SUBROUTINE write_UVgroups_tables

END MODULE UV_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!*************************************************************************
SUBROUTINE locate(xx,n,x,j)
! Locates position j of a value x in an ordered array xx of n elements
! After: xx(j) <= x <= xx(j+1) (assuming increasing order)
! j is lower bound, so it can be zero but not larger than n
!-------------------------------------------------------------------------
  use amr_commons,only:dp
  integer ::  n,j,jl,ju,jm
  real(dp)::  xx(n),x
!-------------------------------------------------------------------------
  jl = 0
  ju = n+1
  do while (ju-jl > 1)
     jm = (ju+jl)/2
     if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
        jl = jm
     else
        ju = jm
     endif
  enddo
  j = jl
END SUBROUTINE locate

!*************************************************************************
SUBROUTINE inp_1d(xax,nx,x,ix0,ix1,dx0,dx1)
! Compute variables by interpolation from table with non-equal intervals.
! xax      => Axis of x-values in table
! nx       => Length of x-axis
! x        => x-value to interpolate to
! ix0,ix1 <=  Lower and upper boundaries of x in xax
! dx0,dx1 <=  Weights of ix0 and ix1 indexes
!-------------------------------------------------------------------------
  use amr_commons,only:dp
  integer:: nx, ix0, ix1
  real(dp), intent(in)::xax(nx), x
  real(dp):: x_step, dx0, dx1
!-------------------------------------------------------------------------
  call locate(xax, nx, x, ix0)
  if(ix0 < 1) ix0=1
  if (ix0 < nx) then
     ix1  = ix0+1
     x_step = xax(ix1) - xax(ix0)
     dx0  = max(           x - xax(ix0), 0.0d0 ) / x_step
     dx1  = min(xax(ix1) - x           , x_step) / x_step
  else
     ix1  = ix0
     dx0  = 0.0d0 ;  dx1  = 1.0d0
  end if

  if (abs(dx0+dx1-1.0d0) .gt. 1.0d-5) then
     write(*,*) 'Screwed up the 1d interpolation ... '
     write(*,*) dx0+dx1
     call clean_stop
  end if
  !ret = dx0 * table(ix1, :) + dx1 * table(ix0, :)

END SUBROUTINE inp_1d
