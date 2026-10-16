module rtz_module
  use amr_parameters, only: dp
  use hydro_parameters, only: n_elements
  implicit none

  private
  public:: elements, n_elements, initialize_elements
  public:: getNe, dust_to_gas_scale_RR14, getMu_RTZ, get_rho_rtz, get_n_rtz

  type Element
      integer(KIND=4) :: atomic_number
      integer(KIND=4) :: n_ions
      integer(KIND=4) :: n_mol
      integer(KIND=4) :: u_hydro_idx
      real(dp)   :: atomic_mass
      real(dp)   :: atomic_mass_g
      real(dp)   :: z_solar
      real(dp)   :: G0_photo_rate
      real(dp)   :: depletion
      real(dp)   :: scale_n
      character(LEN=20):: element_name
      character(LEN=2):: symbol
  end type Element

  ! RTZ STUFF
  
  type(Element) :: elements(n_elements)

CONTAINS

SUBROUTINE initialize_elements()
   ! Initializes the atomic data we need for the RTZ module
   use constants, only: mH, amu2g, mH_amu, mHe_amu, mC_amu, &
                        mN_amu, mO_amu, mNe_amu, mMg_amu, &
                        mSi_amu, mS_amu, mFe_amu
   implicit none
   integer::i
   
   ! Initialize everything to zero
   do i=1,n_elements
      elements(i)%atomic_number = -1
      elements(i)%n_ions = -1
      elements(i)%atomic_mass = -1.0
      elements(i)%atomic_mass_g = -1.0
      elements(i)%z_solar = -1.0
      elements(i)%G0_photo_rate = 0.0
      elements(i)%n_mol = 0
      elements(i)%depletion = 1.0
      elements(i)%element_name = "NANANANANANANANANANA"
      elements(i)%symbol = "XX"
      elements(i)%u_hydro_idx = -1
      elements(i)%scale_n = 1.0
   enddo

   ! Element 1: Hydrogen
   elements(1)%atomic_number = 1
   elements(1)%atomic_mass = mH_amu
   elements(1)%atomic_mass_g = mH
   elements(1)%z_solar = 1.0
   elements(1)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(1)%n_ions = 2
#if N_H2 > 0
      elements(1)%n_mol = N_H2
#else
      elements(1)%n_mol = 0
#endif
   elements(1)%depletion = 1.0
   elements(1)%element_name = "HYDROGEN"
   elements(1)%symbol = "H"

#if N_HELIUM_IONS > 0
   ! Element 2: Helium
   elements(2)%atomic_number = 2
   elements(2)%atomic_mass = mHe_amu
   elements(2)%atomic_mass_g = mHe_amu * amu2g
   elements(2)%z_solar = 8.51E-02
   elements(2)%G0_photo_rate = 0. ! No subionizing PI
   elements(2)%n_ions = N_HELIUM_IONS
   elements(2)%depletion = 1.0
   elements(2)%element_name = "HELIUM"
   elements(2)%symbol = "He"
#endif

#if N_CARBON_IONS > 0
   ! Element 6: Carbon
   elements(6)%atomic_number = 6
   elements(6)%atomic_mass = mC_amu
   elements(6)%atomic_mass_g = mC_amu * amu2g
   elements(6)%z_solar = 2.69E-04
   elements(6)%G0_photo_rate = 3.39E-10
   elements(6)%n_ions = N_CARBON_IONS
   elements(6)%depletion = 0.5
   elements(6)%element_name = "CARBON"
   elements(6)%symbol = "C"
#endif
    
#if N_NITROGEN_IONS > 0
   ! Element 7: Nitrogen
   elements(7)%atomic_number = 7
   elements(7)%atomic_mass = mN_amu
   elements(7)%atomic_mass_g = mN_amu * amu2g
   elements(7)%z_solar = 6.76E-05
   elements(7)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(7)%n_ions = N_NITROGEN_IONS
   elements(7)%depletion = 0.6
   elements(7)%element_name = "NITROGEN"
   elements(7)%symbol = "N"
#endif 

#if N_OXYGEN_IONS > 0
   ! Element 8: Oxygen
   elements(8)%atomic_number = 8
   elements(8)%atomic_mass = mO_amu
   elements(8)%atomic_mass_g = mO_amu * amu2g
   elements(8)%z_solar = 4.90E-04
   elements(8)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(8)%n_ions = N_OXYGEN_IONS
   elements(8)%depletion = 0.73
   elements(8)%element_name = "OXYGEN"
   elements(8)%symbol = "O"
#endif

#if N_NEON_IONS > 0
   ! Element 10: Neon
   elements(10)%atomic_number = 10
   elements(10)%atomic_mass = mNe_amu
   elements(10)%atomic_mass_g = mNe_amu * amu2g
   elements(10)%z_solar = 8.51E-05
   elements(10)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(10)%n_ions = N_NEON_IONS
   elements(10)%depletion = 1.0
   elements(10)%element_name = "NEON"
   elements(10)%symbol = "Ne"
#endif

#if N_MAGNESIUM_IONS > 0
   ! Element 12: Magnesium
   elements(12)%atomic_number = 12
   elements(12)%atomic_mass = mMg_amu
   elements(12)%atomic_mass_g = mMg_amu * amu2g
   elements(12)%z_solar = 3.98E-05
   elements(12)%G0_photo_rate = 6.59E-11 
   elements(12)%n_ions = N_MAGNESIUM_IONS
   elements(12)%depletion = 0.16
   elements(12)%element_name = "MAGNESIUM"
   elements(12)%symbol = "Mg"
#endif

#if N_SILICON_IONS > 0
   ! Element 14: Silicon
   elements(14)%atomic_number = 14
   elements(14)%atomic_mass = mSi_amu
   elements(14)%atomic_mass_g = mSi_amu * amu2g
   elements(14)%z_solar = 3.24E-05
   elements(14)%G0_photo_rate = 4.47E-09
   elements(14)%n_ions = N_SILICON_IONS
   elements(14)%depletion = 0.1
   elements(14)%element_name = "SILICON"
   elements(14)%symbol = "Si"
#endif

#if N_SULFUR_IONS > 0
   ! Element 16: Sulfur
   elements(16)%atomic_number = 16
   elements(16)%atomic_mass = mS_amu
   elements(16)%atomic_mass_g = mS_amu * amu2g
   elements(16)%z_solar = 1.32E-05
   elements(16)%G0_photo_rate = 1.13E-09
   elements(16)%n_ions = N_SULFUR_IONS
   elements(16)%depletion = 1.0
   elements(16)%element_name = "SULFUR"
   elements(16)%symbol = "S"
#endif

#if N_IRON_IONS > 0
   ! Element 26: Iron
   elements(26)%atomic_number = 26
   elements(26)%atomic_mass = mFe_amu
   elements(26)%atomic_mass_g = mFe_amu * amu2g
   elements(26)%z_solar = 3.16E-05
   elements(26)%G0_photo_rate = 4.71E-10
   elements(26)%n_ions = N_IRON_IONS
   elements(26)%depletion = 0.01
   elements(26)%element_name = "IRON"
   elements(26)%symbol = "Fe"
#endif

END SUBROUTINE initialize_elements

!************************************************************************
   FUNCTION getNe(xion,nion) result(ne)
      !Returns the electron number density by looping over all elements 
      !and summing their contributions
      implicit none

      real(dp), intent(in)::xion(1:n_elements,1:n_elements)
      real(dp), intent(in)::nion(1:n_elements)
      real(dp)::ne
      integer::iIons, iElement, n_ions

      ne = 0.d0

      !Loop over all elements
      do iElement=1,n_elements
         if (elements(iElement)%atomic_number .gt. 0) then
            !Get the number of ions
            n_ions = elements(iElement)%n_ions

            !Loop over all ions 
            !Start loop at 2, no electrons in the ground state
            do iIons=2,n_ions
               ne = ne + (nion(iElement) * xion(iElement,iIons) * real(iIons - 1, dp)) 
            end do
         end if
      end do

   END FUNCTION getNe

   FUNCTION dust_to_gas_scale_RR14(log10_O_over_H) result(ratio)
      ! This returns the dust-to-metal ratio relative to the local value
      !------------------------------------
      ! Zsolar: metallicity in solar units
      ! D2Z_solar ~ D2Z / 0.3
      !------------------------------------
      ! based on Remy-Ruyer, Madden, Galliano et al. (2014)
      ! https://arxiv.org/pdf/1312.3442.pdfR
      ! Broken power law with X_{CO,Z} case (Table 1)
      implicit none

      real(dp),intent(in)::log10_O_over_H
      real(dp)::ratio
      real(dp)::a,alphaH,b,alphaL,xt,x,Xsun
      real(dp)::G2D,G2D_sol,y

      ! y = log (G/D)
      ! x = 12 + log10(O/H)
      ! Xsun = 8.69
      ! 
      ! y = a + alphaH * (Xsun - x) for x>xt
      ! y = b + alphaL * (Xsun - x) for x<=xt
      a      = 2.21d0
      alphaH = 1.00d0 ! MW case
      b      = 0.96d0 ! 0.68
      alphaL = 3.10d0 ! 3.08 
      xt     = 8.10d0 ! 7.96
      Xsun   = 8.69d0
      x = max(log10_O_over_H,5.d0) ! Mild extrapolation

      if (log10_O_over_H>xt)then
         y = a + alphaH * (Xsun - x)
      else
         y = b + alphaL * (Xsun - x)
      endif

      G2D = 10.d0**y
      G2D_sol = 10.d0**a

      ratio = max(min(G2D_sol / G2D, 1.d0), 0.d0)

   END FUNCTION dust_to_gas_scale_RR14

#ifdef CO
   FUNCTION getMu_RTZ(ne, element_number_densities, element_ion_fractions, include_H2, nCO) result(mu)
#else
   FUNCTION getMu_RTZ(ne, element_number_densities, element_ion_fractions, include_H2) result(mu)
#endif
      implicit none
      real(dp), intent(in):: ne
      real(dp), intent(in):: element_number_densities(27)
      real(dp), intent(in):: element_ion_fractions(27,27)
      logical, intent(in):: include_H2
#ifdef CO
      real(dp), intent(in):: nCO
#endif
      real(dp):: mu
      real(dp):: m_bar, n_hat

      integer:: i, j

      m_bar = 0.d0
      n_hat = 0.d0

      do i=1,n_elements
         if (elements(i)%atomic_number.gt.0) then
            do j=1,elements(i)%n_ions
               m_bar = m_bar + (element_number_densities(i) * element_ion_fractions(i,j) * elements(i)%atomic_mass)
               n_hat = n_hat + element_number_densities(i) * element_ion_fractions(i,j)
            end do
         end if
      end do

      ! Include electrons
      n_hat = n_hat + ne

      ! Include contribution from H2
      if (include_H2) then
         m_bar = m_bar + (element_number_densities(1) * element_ion_fractions(1,3) * elements(1)%atomic_mass)
         n_hat = n_hat + (0.5d0 * element_number_densities(1) * element_ion_fractions(1,3))
      end if

#ifdef CO
      ! Include contribution from CO
      m_bar = m_bar + nCO * (elements(6)%atomic_mass + elements(8)%atomic_mass)
      n_hat = n_hat + nCO
#endif

      if (n_hat > 1d-30) then
         mu = m_bar / n_hat
      else
         mu = 1.0d0
      end if

   END FUNCTION getMu_RTZ

   FUNCTION get_rho_rtz(element_number_densities) result(rho)
      use constants, only: amu2g
      implicit none
      real(dp), intent(in):: element_number_densities(27)
      real(dp):: rho

      integer:: i

      rho = 0.d0

      do i=1,n_elements
         if (elements(i)%atomic_number.lt.1) then
            cycle
         end if
         rho = rho + (element_number_densities(i) * elements(i)%atomic_mass) ! gives amu/cm^3
      end do

      rho = amu2g * rho ! gives g/cm^3
   END FUNCTION get_rho_rtz

   FUNCTION get_n_rtz(element_number_densities, ne) result(rho_n)
      implicit none
      real(dp), intent(in):: element_number_densities(27)
      real(dp), intent(in):: ne
      real(dp):: rho_n

      integer:: i

      rho_n = 0.d0

      do i=1,n_elements
         if (elements(i)%atomic_number.lt.1) then
            cycle
         end if
         rho_n = rho_n + element_number_densities(i)  ! gives 1/cm^3
      end do

      rho_n = rho_n + ne
   END FUNCTION get_n_rtz

end module rtz_module