MODULE coolrates_module
   ! Module for returning cooling and thermochemistry interaction rates.
   ! The temperature dependence is tabulated in rate (nonlog) versus the
   ! log of T, because it is expensive to calculate on the fly.
   ! The rates are interpolated using cubic splines, and extrapolated in
   ! log-log space if temperature is above table boundaries.
   ! Joki Rosdahl and Andreas Bleuler, September 2015.

   use amr_parameters, only: dp
   use rt_parameters, only: nIons, isH2, isHe, ixHI, ixHII, ixHeII, ixHeIII, isH2Katz, H2clumping
   implicit none

   private   ! default
   public init_metal_atomic_data_tables, init_coolrates_tables, update_coolrates_tables, inp_coolrates_table &
      , compCoolrate, tbl_alphaA_HII, tbl_alphaA_HeII, tbl_alphaA_HeIII &
      , tbl_alphaB_HII, tbl_alphaB_HeII, tbl_alphaB_HeIII, tbl_beta_HI &
      , tbl_beta_HeI, tbl_beta_HeII, tbl_cr_ci_HI, tbl_cr_ci_HeI &
      , tbl_cr_ci_HeII, tbl_cr_ce_HI, tbl_cr_ce_HeI, tbl_cr_ce_HeII &
      , tbl_cr_r_HII, tbl_cr_r_HeII, tbl_cr_r_HeIII, tbl_cr_bre &
      , tbl_cr_com, tbl_cr_die &
      , Epump, comp_SH2, comp_Sd, comp_Alpha_H2, comp_Beta_H2coll &
      , comp_Beta_H2HI, PE_efficiency, comp_Beta_H2coll_new, H2_cooling_GA08 &
      , comp_Alpha_oxygen, comp_Alpha_nitrogen, comp_Alpha_carbon &
      , comp_Alpha_magnesium, comp_Alpha_silicon &
      , comp_Alpha_sulfur, comp_Alpha_iron, comp_Alpha_neon  &
      , comp_Beta_oxygen, comp_Beta_nitrogen, comp_Beta_carbon &
      , comp_Beta_magnesium, comp_Beta_silicon &
      , comp_Beta_sulfur, comp_Beta_iron, comp_Beta_neon &
      , comp_cr_co, comp_de_co, update_co_chem, comp_SCO, comp_co_cooling &
      , cr_ionization_oxygen, cr_ionization_nitrogen, cr_ionization_carbon &
      , cr_ionization_magnesium, cr_ionization_silicon, cr_ionization_sulfur &
      , cr_ionization_iron, cr_ionization_neon, cr_ionization_ind_UV_CI &
      , cr_ionization_ind_UV_OI, cr_ionization_ind_UV_NI, cr_ionization_ind_UV_MgI &
      , cr_ionization_ind_UV_SI, cr_ionization_ind_UV_SiI, cr_ionization_ind_UV_FeI &
      , cr_ionization_ind_UV_NeI, DUST_RECOMBINATION

   ! Default cooling rates table parameters
   integer, parameter     :: nbinT = 1001
   real(dp), parameter    :: Tmin = 1d-2
   real(dp), parameter    :: Tmax = 1d+9
   real(dp)              :: dlogTinv ! Inverse of the bin space (in K)
   real(dp)              :: hTable, h2Table, h3Table   ! Interpol constants
   real(dp)              :: one_over_lnTen, one_over_hTable, one_over_h2Table
   real(dp)              :: three_over_h2Table, two_over_h3Table

   real(dp), dimension(nbinT) :: T_lookup = 0d0 ! Lookup temperature in log K

   ! Cosmic ray ionization tables for metals --> these are direct rates
   ! Data taken from CHIMES (Alex Richings et al)
   real(dp), dimension(8):: cr_ionization_oxygen = (/ 5.637d0, 1.9201249d0, &
      & 0.98636097d0, 0.5261845d0, 0.36398333d0, 0.21907325d0, 0.03679156d0, 0.015607069d0 /)
   real(dp), dimension(7):: cr_ionization_nitrogen = (/ 4.52d0, 1.8237852d0, 0.8574613d0, &
      & 0.5265334d0, 0.30504152d0, 0.04926644d0, 0.020386748d0 /)
   real(dp), dimension(6):: cr_ionization_carbon = (/ 3.83d0, 1.6638695d0, 0.8312262d0, &
      & 0.4541473d0, 0.06937006d0, 0.027755102d0 /)
   real(dp), dimension(12):: cr_ionization_magnesium = (/ 9.016716d0, 5.037889d0, 1.345602d0, &
      & 0.86565775d0, 0.5749281d0, 0.36411676d0, 0.24167651d0, 0.15331018d0, 0.13045996d0, & 
      & 0.084565096d0, 0.015437003d0, 0.006928171d0 /)
   real(dp), dimension(14):: cr_ionization_silicon = (/ 6.37329d0, 2.4664342d0, 2.5934594d0, &
      & 1.8422631d0, 0.64957255d0, 0.46263242d0, 0.3302073d0, 0.22406998d0, 0.15485357d0, & 
      & 0.10161172d0, 0.09062556d0, 0.059769776d0, 0.011156686d0, 0.0050879163d0 /)
   real(dp), dimension(16):: cr_ionization_sulfur = (/ 7.5164065d0, 2.8824131d0, 1.5510099d0, &
      & 0.8586138d0, 1.2691398d0, 0.9713803d0, 0.38636714d0, 0.28946567d0, 0.21488377d0, 0.15199806d0, &
      & 0.107721165d0, 0.07223445d0, 0.066572525d0, 0.044464245d0, 0.008436725d0, 0.003892387d0 /)
   real(dp), dimension(26):: cr_ionization_iron = (/ 13.192021d0, 5.842174d0, 4.724112d0, &
      & 2.5858068d0, 1.7483323d0, 1.2082134d0, 0.85902447d0, 0.62585855d0, 0.46490052d0, &
      & 0.3626024d0, 0.28070357d0, 0.20538796d0, 0.15056798d0, 0.10395802d0, 0.21419027d0, &
      & 0.18064684d0, 0.086154655d0, 0.07006637d0, 0.056015324d0, 0.042974997d0, 0.032204125d0, &
      & 0.02267751d0, 0.022580078d0, 0.015563778d0, 0.0030807566d0, 0.0014658332d0 /)
   real(dp), dimension(10):: cr_ionization_neon = (/ 4.8251953d0, 2.2843935d0, 1.2726027d0, &
      & 0.6982404d0, 0.43035796d0, 0.2581805d0, 0.204098d0, 0.12864962d0, 0.022742474d0, &
      & 0.009985316d0 /)
   ! Cosmic ray data for ionization from induced UV emission
   ! Data from: https://home.strw.leidenuniv.nl/~ewine/photo/cosmic_ray_rates.html
   ! Note that this assumes 10-16 s-1 H2-1 for the primary H2 cosmic ray dissociation rate
   real(dp)::cr_ionization_ind_UV_CI  = 2.60d-14
   real(dp)::cr_ionization_ind_UV_OI  = 2.70d-16
   real(dp)::cr_ionization_ind_UV_NI  = 7.34d-17
   real(dp)::cr_ionization_ind_UV_MgI = 1.12d-14
   real(dp)::cr_ionization_ind_UV_SI  = 7.91d-14
   real(dp)::cr_ionization_ind_UV_SiI = 4.16d-13
   real(dp)::cr_ionization_ind_UV_FeI = 4.81d-14
   real(dp)::cr_ionization_ind_UV_NeI = 0.00d+00

   type coolrates_table
      ! Cooling and interaction rates (log):
      real(dp), dimension(nbinT)::rates = 0d0
      ! Temperature derivatives of those rates (drate/dlog(T)):
      real(dp), dimension(nbinT)::primes = 0d0
   end type coolrates_table

   type(coolrates_table), save::tbl_alphaA_HII ! Case A rec. coefficients
   type(coolrates_table), save::tbl_alphaA_HeII
   type(coolrates_table), save::tbl_alphaA_HeIII
   type(coolrates_table), save::tbl_alphaB_HII ! Case B rec. coefficients
   type(coolrates_table), save::tbl_alphaB_HeII
   type(coolrates_table), save::tbl_alphaB_HeIII
   type(coolrates_table), save::tbl_beta_HI ! Collisional ionisation rates
   type(coolrates_table), save::tbl_beta_HeI
   type(coolrates_table), save::tbl_beta_HeII
   type(coolrates_table), save::tbl_cr_ci_HI ! Coll. ionisation cooling
   type(coolrates_table), save::tbl_cr_ci_HeI
   type(coolrates_table), save::tbl_cr_ci_HeII
   type(coolrates_table), save::tbl_cr_ce_HI ! Coll. excitation cooling
   type(coolrates_table), save::tbl_cr_ce_HeI
   type(coolrates_table), save::tbl_cr_ce_HeII
   type(coolrates_table), save::tbl_cr_r_HII ! Recombination cooling
   type(coolrates_table), save::tbl_cr_r_HeII
   type(coolrates_table), save::tbl_cr_r_HeIII
   type(coolrates_table), save::tbl_cr_bre ! Bremsstrahlung cooling rates
   type(coolrates_table), save::tbl_cr_com ! Compton cooling rates
   type(coolrates_table), save::tbl_cr_die ! Dielectronic cooling rates

   ! Metal atomic data
   real(dp),dimension(8,6)::RR_rates_oxygen
   real(dp),dimension(8,9)::DR_rates_c_oxygen, DR_rates_e_oxygen
   real(dp),dimension(7,6)::RR_rates_nitrogen
   real(dp),dimension(7,9)::DR_rates_c_nitrogen, DR_rates_e_nitrogen
   real(dp),dimension(6,6)::RR_rates_carbon
   real(dp),dimension(6,9)::DR_rates_c_carbon, DR_rates_e_carbon
   real(dp),dimension(12,6)::RR_rates_magnesium
   real(dp),dimension(12,9)::DR_rates_c_magnesium, DR_rates_e_magnesium
   real(dp),dimension(14,6)::RR_rates_silicon
   real(dp),dimension(14,9)::DR_rates_c_silicon, DR_rates_e_silicon
   real(dp),dimension(15,6)::RR_rates_sulfur
   real(dp),dimension(15,9)::DR_rates_c_sulfur, DR_rates_e_sulfur
   real(dp),dimension(16,2)::RR_rates_alt_sulfur
   real(dp),dimension(16,4)::DR_rates_alt_sulfur
   real(dp),dimension(15,6)::RR_rates_iron
   real(dp),dimension(15,9)::DR_rates_c_iron, DR_rates_e_iron
   real(dp),dimension(26,2)::RR_rates_alt_iron
   real(dp),dimension(26,8)::DR_rates_alt_iron
   real(dp),dimension(10,6)::RR_rates_neon
   real(dp),dimension(10,9)::DR_rates_c_neon, DR_rates_e_neon

CONTAINS
   SUBROUTINE init_metal_atomic_data_tables()
      implicit none

      !!!!!!!!
      !!!!!!!! Oxygen
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_oxygen(1,:)  = (/ 6.552d-10, 0.7470d0, 1.951d2,   4.483d7, 0.0000d0, 0.000d0 /)
      RR_rates_oxygen(2,:)  = (/ 2.652d-10, 0.6705d0, 5.842d2,   4.559d7, 0.0000d0, 0.000d0 /)
      RR_rates_oxygen(3,:)  = (/ 8.193d-11, 0.5165d0, 2.392d3,   2.487d7, 0.0000d0, 0.000d0 /)
      RR_rates_oxygen(4,:)  = (/ 1.724d-10, 0.6556d0, 3.372d2,   1.030d7, 0.0000d0, 0.000d0 /)
      RR_rates_oxygen(5,:)  = (/ 3.955d-09, 0.7813d0, 6.821d-01, 5.076d6, 0.0000d0, 0.000d0 /)
      RR_rates_oxygen(6,:)  = (/ 2.501d-09, 0.7844d0, 5.235d-01, 4.470d6, 0.0447d0, 1.642d5 /)
      RR_rates_oxygen(7,:)  = (/ 2.096d-09, 0.7668d0, 1.602d-01, 4.377d6, 0.1070d0, 1.392d5 /)
      RR_rates_oxygen(8,:)  = (/ 6.622d-11, 0.6109d0, 4.136d0,   4.214d6, 0.4093d0, 8.770d4 /)

      ! Dielectronic recombination rates
      DR_rates_c_oxygen(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(2,:)   = (/ 4.925E-03, 5.837E-02, 1.359E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(3,:)   = (/ 6.135E-02, 1.968E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(4,:)   = (/ 2.389E-05, 1.355E-04, 5.885E-03, 2.163E-03, 6.341E-04, 1.348E-02, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(5,:)   = (/ 1.615E-05, 9.299E-06, 1.530E-04, 6.616E-04, 1.080E-02, 7.503E-04, 2.892E-03, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(6,:)   = (/ 3.932E-07, 2.523E-07, 3.447E-05, 5.776E-03, 5.101E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(7,:)   = (/ 1.627E-07, 1.262E-07, 6.663E-07, 3.925E-06, 2.406E-03, 1.146E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_oxygen(8,:)   = (/ 5.629E-08, 2.550E-07, 6.173E-04, 1.627E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
   
      DR_rates_e_oxygen(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(2,:)   = (/ 5.440E+06, 7.170E+06, 1.152E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(3,:)   = (/ 6.113E+06, 3.656E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(4,:)   = (/ 2.326E+04, 3.209E+04, 1.316E+05, 6.731E+05, 1.892E+06, 6.150E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(5,:)   = (/ 7.569E+02, 3.659E+03, 1.984E+04, 8.429E+04, 2.294E+05, 1.161E+06, 6.137E+06, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(6,:)   = (/ 1.509E+02, 6.211E+02, 1.562E+04, 1.936E+05, 4.700E+05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(7,:)   = (/ 4.535E+01, 2.847E+02, 4.166E+03, 2.877E+04, 1.953E+05, 3.646E+05, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_oxygen(8,:)   = (/ 5.395E+03, 1.770E+04, 1.671E+05, 2.687E+05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)

      !!!!!!!!
      !!!!!!!! Nitrogen
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_nitrogen(1,:)  = (/ 6.170d-10, 0.7481d0, 1.316d2,   3.427d7, 0.0000d0, 0.000d0 /)
      RR_rates_nitrogen(2,:)  = (/ 2.388d-10, 0.6732d0, 3.960d2,   3.583d7, 0.0000d0, 0.000d0 /)
      RR_rates_nitrogen(3,:)  = (/ 6.245d-11, 0.4985d0, 1.957d3,   2.177d7, 0.0000d0, 0.000d0 /)
      RR_rates_nitrogen(4,:)  = (/ 1.533d-10, 0.6682d0, 1.823d2,   7.751d6, 0.0000d0, 0.000d0 /)
      RR_rates_nitrogen(5,:)  = (/ 7.923d-10, 0.7768d0, 3.750d0,   3.468d6, 0.0223d0, 7.206d4 /)
      RR_rates_nitrogen(6,:)  = (/ 2.410d-09, 0.7948d0, 1.231d-01, 3.016d6, 0.0774d0, 1.016d5 /)
      RR_rates_nitrogen(7,:)  = (/ 6.387d-10, 0.7308d0, 9.467d-02, 2.954d6, 0.2440d0, 6.739d4 /)

      ! Dielectronic recombination rates
      DR_rates_c_nitrogen(1,:)   = (/ 0.000E+00, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_nitrogen(2,:)   = (/ 2.801E-03, 4.362E-02,  1.117E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_nitrogen(3,:)   = (/ 5.761E-03, 3.434E-02, -1.660E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_nitrogen(4,:)   = (/ 2.040E-06, 6.986E-05,  3.168E-04, 4.353E-03, 7.765E-04, 5.101E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_nitrogen(5,:)   = (/ 3.386E-06, 3.036E-05,  5.945E-05, 1.195E-03, 6.462E-03, 1.358E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_nitrogen(6,:)   = (/ 7.712E-08, 4.839E-08,  2.218E-06, 1.536E-03, 3.647E-03, 4.234E-05, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_nitrogen(7,:)   = (/ 1.658E-08, 2.760E-08,  2.391E-09, 7.585E-07, 3.012E-04, 7.132E-04, 0.000E+00, 0.000E+00, 0.000E+00 /)
   
      DR_rates_e_nitrogen(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_nitrogen(2,:)   = (/ 4.198E+06, 5.516E+06, 8.050E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_nitrogen(3,:)   = (/ 3.860E+06, 4.883E+06, 6.259E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_nitrogen(4,:)   = (/ 3.084E+03, 1.332E+04, 6.475E+04, 1.181E+05, 6.687E+05, 4.778E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_nitrogen(5,:)   = (/ 1.406E+03, 6.965E+03, 2.604E+04, 1.304E+05, 1.965E+05, 4.466E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_nitrogen(6,:)   = (/ 7.113E+01, 2.765E+02, 1.439E+04, 1.347E+05, 2.496E+05, 2.204E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_nitrogen(7,:)   = (/ 1.265E+01, 8.425E+01, 2.964E+02, 5.923E+03, 1.278E+05, 2.184E+05, 0.000E+00, 0.000E+00, 0.000E+00 /)

      !!!!!!!!
      !!!!!!!! Carbon
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_carbon(1,:)  = (/ 5.337d-10, 0.7485d0, 9.502d1,   2.517d7, 0.0000d0, 0.000d0 /)
      RR_rates_carbon(2,:)  = (/ 2.044d-10, 0.6742d0, 2.647d2,   2.773d7, 0.0000d0, 0.000d0 /)
      RR_rates_carbon(3,:)  = (/ 4.798d-11, 0.4834d0, 1.355d3,   1.872d7, 0.0000d0, 0.000d0 /)
      RR_rates_carbon(4,:)  = (/ 1.120d-10, 0.6737d0, 1.115d2,   5.938d6, 0.0000d0, 0.000d0 /)
      RR_rates_carbon(5,:)  = (/ 2.067d-09, 0.8012d0, 1.643d-01, 2.172d6, 0.0427d0, 6.341d4 /)
      RR_rates_carbon(6,:)  = (/ 2.995d-09, 0.7849d0, 6.670d-03, 1.943d6, 0.1597d0, 4.955d4 /)

      ! Dielectronic recombination rates
      DR_rates_c_carbon(1,:)   = (/ 0.000E+00, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_carbon(2,:)   = (/ 1.426E-03, 3.046E-02,  8.373E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_carbon(3,:)   = (/ 2.646E-03, 1.762E-02, -7.843E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_carbon(4,:)   = (/ 4.673E-07, 1.887E-05,  1.305E-05, 3.099E-03, 3.001E-04, 2.553E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_carbon(5,:)   = (/ 3.489E-06, 2.222E-07,  1.954E-05, 4.212E-03, 2.037E-04, 2.936E-04, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_carbon(6,:)   = (/ 6.346E-09, 9.793E-09,  1.634E-06, 8.369E-04, 3.355E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
   
      DR_rates_e_carbon(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_carbon(2,:)   = (/ 3.116E+06, 4.075E+06, 5.749E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_carbon(3,:)   = (/ 2.804E+06, 3.485E+06, 4.324E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_carbon(4,:)   = (/ 7.233E+02, 2.847E+03, 1.054E+04, 8.915E+04, 2.812E+05, 3.254E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_carbon(5,:)   = (/ 2.660E+03, 3.756E+03, 2.566E+04, 1.400E+05, 1.801E+06, 4.307E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_carbon(6,:)   = (/ 1.217E+01, 7.380E+01, 1.523E+04, 1.207E+05, 2.144E+05, 3.254E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)

      !!!!!!!!
      !!!!!!!! Magnesium
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_magnesium(1,:)  = (/ 1.022d-09, 0.7476d0, 4.098d2, 1.011d8, 0.0000d0, 0.000d0 /)
      RR_rates_magnesium(2,:)  = (/ 4.214d-10, 0.6713d0, 1.396d3, 9.433d7, 0.0000d0, 0.000d0 /)
      RR_rates_magnesium(3,:)  = (/ 1.602d-10, 0.5492d0, 4.944d3, 4.434d7, 0.0000d0, 0.000d0 /)
      RR_rates_magnesium(4,:)  = (/ 3.445d-10, 0.6553d0, 8.693d2, 2.196d7, 0.0000d0, 0.000d0 /)
      RR_rates_magnesium(5,:)  = (/ 3.989d-10, 1.0231d0, 2.601d1, 1.227d7, 0.0000d0, 0.000d0 /)
      RR_rates_magnesium(6,:)  = (/ 3.859d-09, 0.7579d0, 5.587d0, 1.235d7, 0.0000d0, 0.000d0 /)
      RR_rates_magnesium(7,:)  = (/ 9.133d-10, 0.7353d0, 3.674d1, 1.263d7, 0.0211d0, 4.049d5 /)
      RR_rates_magnesium(8,:)  = (/ 7.515d-10, 0.7203d0, 2.582d1, 1.355d7, 0.0436d0, 5.691d5 /)
      RR_rates_magnesium(9,:)  = (/ 4.031d-10, 0.6803d0, 3.205d1, 1.626d7, 0.0764d0, 5.399d5 /)
      RR_rates_magnesium(10,:) = (/ 1.249d-10, 0.5600d0, 7.748d1, 2.015d7, 0.1917d0, 5.139d5 /)
      RR_rates_magnesium(11,:) = (/ 1.345d-11, 0.1074d0, 7.877d2, 7.925d7, 0.4631d0, 5.027d5 /)
      RR_rates_magnesium(12,:) = (/ 5.452d-11, 0.6845d0, 5.637d0, 1.551d6, 0.3945d0, 8.360d5 /)

      ! Dielectronic recombination rates
      DR_rates_c_magnesium(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(2,:)   = (/ 2.262E-02, 1.216E-01, 2.531E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(3,:)   = (/ 3.067E-02, 1.375E-01, 1.347E-02, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(4,:)   = (/ 2.582E-04, 5.033E-04, 4.561E-03, 8.754E-03, 2.734E-02, 7.509E-02, 1.385E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(5,:)   = (/ 3.565E-05, 2.017E-05, 1.165E-03, 5.016E-03, 2.434E-02, 2.508E-02, 2.475E-02, 7.087E-04, 0.000E+00 /)
      DR_rates_c_magnesium(6,:)   = (/ 5.451E-05, 6.999E-05, 3.928E-04, 7.483E-04, 1.851E-02, 1.190E-02, 4.764E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(7,:)   = (/ 2.931E-05, 5.192E-05, 1.769E-04, 7.988E-04, 1.488E-02, 4.546E-02, 1.505E-04, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(8,:)   = (/ 2.407E-05, 1.187E-04, 3.845E-04, 1.333E-02, 2.804E-02, 8.746E-04, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(9,:)   = (/ 2.574E-06, 2.001E-06, 1.829E-05, 2.362E-05, 6.413E-03, 4.118E-03, 7.224E-03, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(10,:)  = (/ 9.400E-08, 3.818E-07, 2.064E-07, 2.710E-05, 3.802E-03, 3.086E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(11,:)  = (/ 6.269E-06, 9.181E-04, 3.082E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_magnesium(12,:)  = (/ 3.871E-08, 4.732E-07, 1.599E-03, 2.628E-05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)

      DR_rates_e_magnesium(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(2,:)   = (/ 1.201E+07, 1.588E+07, 2.473E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(3,:)   = (/ 1.136E+07, 1.431E+07, 1.762E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(4,:)   = (/ 2.861E+04, 4.087E+04, 1.528E+05, 2.677E+05, 2.055E+06, 1.298E+07, 1.852E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(5,:)   = (/ 1.759E+03, 1.237E+04, 5.651E+04, 1.680E+05, 4.007E+05, 1.991E+06, 1.330E+07, 4.060E+07, 0.000E+00 /)
      DR_rates_e_magnesium(6,:)   = (/ 2.196E+02, 4.524E+03, 2.491E+04, 9.111E+04, 3.367E+05, 8.316E+05, 1.817E+06, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(7,:)   = (/ 7.667E+02, 3.444E+03, 2.409E+04, 1.190E+05, 3.999E+05, 1.502E+06, 2.464E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(8,:)   = (/ 8.061E+03, 2.154E+04, 6.640E+04, 3.406E+05, 1.327E+06, 3.347E+06, 1.852E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(9,:)   = (/ 1.129E+02, 1.585E+03, 1.193E+04, 4.151E+04, 3.684E+05, 6.925E+05, 1.150E+06, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(10,:)  = (/ 2.660E+02, 1.141E+03, 3.210E+03, 1.727E+05, 4.523E+05, 9.105E+05, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(11,:)  = (/ 4.104E+05, 5.766E+05, 7.310E+05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_magnesium(12,:)  = (/ 8.415E+03, 1.682E+04, 5.000E+04, 2.759E+05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)


      !!!!!!!!
      !!!!!!!! Silicon
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_silicon(1,:)  = (/ 1.261E-09, 0.7488E+00, 5.068E+02, 1.365E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(2,:)  = (/ 4.870E-10, 0.6697E+00, 2.026E+03, 1.265E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(3,:)  = (/ 2.017E-10, 0.5588E+00, 6.494E+03, 5.693E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(4,:)  = (/ 4.633E-10, 0.6602E+00, 1.088E+03, 2.896E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(5,:)  = (/ 1.851E-09, 0.7384E+00, 6.906E+01, 1.644E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(6,:)  = (/ 1.688E-09, 0.7390E+00, 5.549E+01, 1.716E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(7,:)  = (/ 2.100E-09, 0.7401E+00, 2.523E+01, 1.842E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_silicon(8,:)  = (/ 7.532E-10, 0.7072E+00, 8.860E+01, 1.997E+07, 0.0185E+00, 6.949E+05 /)
      RR_rates_silicon(9,:)  = (/ 4.615E-10, 0.6753E+00, 1.143E+02, 2.377E+07, 0.0356E+00, 8.595E+05 /)
      RR_rates_silicon(10,:) = (/ 2.468E-10, 0.6113E+00, 1.649E+02, 3.231E+07, 0.0636E+00, 9.837E+05 /)
      RR_rates_silicon(11,:) = (/ 5.134E-11, 0.3678E+00, 1.009E+03, 8.514E+07, 0.1646E+00, 1.084E+06 /)
      RR_rates_silicon(12,:) = (/ 6.739E-11, 0.4931E+00, 2.166E+02, 4.491E+07, 0.1667E+00, 9.046E+05 /)
      RR_rates_silicon(13,:) = (/ 1.964E-10, 0.6287E+00, 7.712E+00, 2.951E+07, 0.1523E+00, 4.804E+05 /)
      RR_rates_silicon(14,:) = (/ 3.262E-11, 0.6270E+00, 1.590E+01, 4.237E+07, 0.2333E+00, 5.828E+04 /)

      ! Dielectronic recombination rates
      DR_rates_c_silicon(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(2,:)   = (/ 3.846E-02, 1.491E-01, 2.779E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(3,:)   = (/ 5.318E-02, 1.874E-01, 1.227E-02, 7.173E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(4,:)   = (/ 1.205E-03, 1.309E-02, 5.333E-03, 2.858E-02, 3.195E-02, 1.433E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(5,:)   = (/ 2.214E-04, 1.260E-03, 3.832E-03, 1.158E-02, 2.871E-02, 5.456E-02, 5.035E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(6,:)   = (/ 1.246E-04, 6.649E-04, 2.912E-03, 2.912E-02, 3.049E-02, 1.056E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(7,:)   = (/ 5.845E-04, 6.600E-04, 7.180E-04, 1.714E-02, 2.958E-02, 1.075E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(8,:)   = (/ 5.272E-05, 2.282E-04, 1.345E-03, 2.246E-02, 9.606E-02, 8.366E-04, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(9,:)   = (/ 2.086E-06, 9.423E-06, 3.423E-05, 3.950E-04, 1.535E-02, 4.986E-02, 4.067E-04, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(10,:)  = (/ 7.163E-07, 2.656E-06, 1.119E-06, 4.796E-05, 4.052E-03, 6.101E-03, 2.366E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(11,:)  = (/ 1.422E-04, 9.474E-03, 1.650E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(12,:)  = (/ 3.819E-06, 2.421E-05, 2.283E-04, 8.604E-03, 2.617E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(13,:)  = (/ 2.930E-06, 2.803E-06, 9.023E-05, 6.909E-03, 2.582E-05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_silicon(14,:)  = (/ 3.408E-08, 1.913E-07, 1.679E-07, 7.523E-07, 8.386E-05, 4.083E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)

      DR_rates_e_silicon(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(2,:)   = (/ 1.627E+07, 2.154E+07, 3.827E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(3,:)   = (/ 1.552E+07, 1.969E+07, 2.532E+07, 2.696E+08, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(4,:)   = (/ 4.824E+04, 2.137E+05, 5.492E+05, 2.460E+06, 3.647E+06, 1.889E+07, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(5,:)   = (/ 3.586E+03, 1.324E+04, 5.636E+04, 2.498E+05, 5.363E+05, 2.809E+06, 1.854E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(6,:)   = (/ 1.167E+03, 9.088E+03, 5.332E+04, 4.219E+05, 1.571E+06, 2.721E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(7,:)   = (/ 9.856E+02, 6.577E+03, 4.281E+04, 3.944E+05, 1.296E+06, 2.441E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(8,:)   = (/ 2.829E+03, 2.617E+04, 1.374E+05, 4.520E+05, 2.072E+06, 7.808E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(9,:)   = (/ 3.708E+02, 3.870E+03, 2.226E+04, 1.318E+05, 5.285E+05, 1.742E+06, 8.392E+06, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(10,:)  = (/ 5.625E+02, 2.952E+03, 9.682E+03, 1.473E+05, 5.064E+05, 8.047E+05, 1.623E+06, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(11,:)  = (/ 7.685E+05, 1.208E+06, 1.839E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(12,:)  = (/ 3.802E+03, 1.280E+04, 5.953E+04, 1.026E+05, 1.154E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(13,:)  = (/ 1.162E+02, 5.721E+03, 3.477E+04, 1.176E+05, 3.505E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_silicon(14,:)  = (/ 2.431E+01, 1.293E+02, 4.272E+02, 3.729E+03, 5.514E+04, 1.295E+05, 0.000E+00, 0.000E+00, 0.000E+00 /)

      !!!!!!!!
      !!!!!!!! Sulfur
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_sulfur(1,:)  = (/ 1.432E-09, 0.7485E+00, 6.688E+02, 1.793E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(2,:)  = (/ 5.546E-10, 0.6692E+00, 2.754E+03, 1.633E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(3,:)  = (/ 2.362E-10, 0.5615E+00, 8.776E+03, 7.208E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(4,:)  = (/ 5.511E-10, 0.6598E+00, 1.492E+03, 3.755E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(5,:)  = (/ 1.740E-09, 0.7303E+00, 1.494E+02, 2.193E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(6,:)  = (/ 1.702E-09, 0.7301E+00, 1.138E+02, 2.253E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(7,:)  = (/ 1.518E-09, 0.7246E+00, 9.845E+01, 2.390E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(8,:)  = (/ 1.137E-09, 0.7080E+00, 1.099E+02, 2.745E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(9,:)  = (/ 8.773E-10, 0.6853E+00, 1.115E+02, 3.386E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_sulfur(10,:) = (/ 3.384E-10, 0.6175E+00, 3.380E+02, 4.347E+07, 0.0225E+00, 1.672E+06 /)
      RR_rates_sulfur(11,:) = (/ 9.565E-11, 0.4517E+00, 1.599E+03, 9.252E+07, 0.0612E+00, 1.986E+06 /)
      RR_rates_sulfur(12,:) = (/ 1.588E-10, 0.5584E+00, 3.350E+02, 5.188E+07, 0.0591E+00, 1.656E+06 /)
      RR_rates_sulfur(13,:) = (/ 2.615E-10, 0.6343E+00, 6.238E+01, 2.803E+07, 0.0773E+00, 1.059E+06 /)
      RR_rates_sulfur(14,:) = (/ 3.043E-10, 0.6947E+00, 1.678E+01, 2.050E+07, 0.0795E+00, 6.868E+04 /)
      RR_rates_sulfur(15,:) = (/ 2.478E-11, 0.4642E+00, 3.294E+02, 2.166E+07, 0.3351E+00, 7.630E+05 /)

      ! Dielectronic recombination rates
      DR_rates_c_sulfur(1,:)   = (/ 0.000E+00, 0.000E+00,  0.000E+00, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(2,:)   = (/ 6.659E-02, 1.762E-01, -6.522E-03, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(3,:)   = (/ 8.410E-02, 2.381E-01,  1.065E-02, 1.049E-03,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(4,:)   = (/ 1.173E-03, 1.718E-03,  1.657E-02, 7.474E-03,  9.698E-02, 1.399E-01, 7.029E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(5,:)   = (/ 2.610E-04, 9.442E-04,  8.190E-03, 4.271E-02,  1.997E-02, 9.590E-02, 8.444E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(6,:)   = (/ 9.126E-05, 1.691E-04,  3.050E-03, 2.604E-02,  3.245E-02, 2.511E-01, 4.459E-04, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(7,:)   = (/ 1.371E-04, 3.096E-04,  1.782E-03, 1.278E-02,  3.323E-02, 2.593E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(8,:)   = (/ 1.830E-04, 1.551E-03,  1.717E-03, 2.798E-02,  6.933E-02, 1.727E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(9,:)   = (/ 2.924E-05, 3.366E-04,  2.104E-04, 1.910E-02, -4.017E-04, 6.541E-02, 9.546E-02, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(10,:)  = (/ 2.585E-06, 9.517E-06,  5.194E-06, 3.715E-04,  1.553E-02, 1.013E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(11,:)  = (/ 3.931E-05, 4.431E-03,  5.156E-02, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(12,:)  = (/ 2.816E-06, 3.172E-05,  1.832E-04, 4.360E-03,  1.618E-02, 7.707E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(13,:)  = (/ 9.571E-06, 6.268E-05,  3.807E-04, 1.874E-02,  5.526E-03, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(14,:)  = (/ 5.817E-07, 1.391E-06,  1.123E-05, 1.521E-04,  1.875E-03, 2.097E-02, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_sulfur(15,:)  = (/ 3.040E-07, 4.393E-07,  1.609E-06, 4.980E-06,  3.457E-05, 8.617E-03, 9.284E-04, 0.000E+00, 0.000E+00 /)
      
      DR_rates_e_sulfur(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(2,:)   = (/ 2.122E+07, 2.897E+07, 5.786E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(3,:)   = (/ 2.032E+07, 2.592E+07, 3.206E+07, 2.016E+08, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(4,:)   = (/ 1.366E+04, 6.371E+04, 2.540E+05, 6.493E+05, 3.868E+06, 2.186E+07, 2.983E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(5,:)   = (/ 5.114E+03, 1.750E+04, 1.061E+05, 4.691E+05, 1.821E+06, 4.033E+06, 2.476E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(6,:)   = (/ 3.575E+03, 1.491E+04, 9.388E+04, 3.856E+05, 9.630E+05, 3.490E+06, 1.743E+08, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(7,:)   = (/ 1.505E+03, 1.100E+04, 6.720E+04, 3.364E+05, 9.686E+05, 3.195E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(8,:)   = (/ 1.122E+04, 3.121E+04, 9.348E+04, 4.792E+05, 2.052E+06, 3.291E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(9,:)   = (/ 7.192E+02, 7.490E+03, 4.510E+04, 5.586E+05, 3.830E+05, 1.976E+06, 2.894E+06, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(10,:)  = (/ 1.277E+03, 5.978E+03, 2.097E+04, 2.292E+05, 7.974E+05, 2.376E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(11,:)  = (/ 9.455E+05, 1.365E+06, 2.169E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(12,:)  = (/ 7.590E+03, 1.558E+04, 4.013E+04, 1.156E+05, 1.601E+05, 1.839E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(13,:)  = (/ 1.180E+03, 6.443E+03, 2.264E+04, 1.530E+05, 3.564E+05, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(14,:)  = (/ 3.628E+02, 1.058E+03, 7.160E+03, 3.260E+04, 1.235E+05, 2.070E+05, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_sulfur(15,:)  = (/ 5.016E+01, 3.266E+02, 3.102E+03, 1.210E+04, 4.969E+04, 2.010E+05, 2.575E+05, 0.000E+00, 0.000E+00 /)

      ! radiative recombination data from https://www.pa.uky.edu/~verner/rec.html
      ! The index number is the number of electons on the atom before the recombination process + 1
      RR_rates_alt_sulfur(16,:) = (/ 4.10E-13, 0.630E+00 /) !SII --> SI
      RR_rates_alt_sulfur(15,:) = (/ 1.80E-12, 0.686E+00 /) !SIII --> SII
      RR_rates_alt_sulfur(14,:) = (/ 2.70E-12, 0.745E+00 /) !SIV --> SIII
      RR_rates_alt_sulfur(13,:) = (/ 5.70E-12, 0.755E+00 /) !SV --> SIV
      RR_rates_alt_sulfur(12,:) = (/ 1.20E-11, 0.701E+00 /) !SVI --> SV
      RR_rates_alt_sulfur(11,:) = (/ 1.70E-11, 0.849E+00 /) !SVII --> SVI
      RR_rates_alt_sulfur(10,:) = (/ 2.70E-11, 0.733E+00 /) !SVIII --> SVII
      RR_rates_alt_sulfur(9,:)  = (/ 4.00E-11, 0.696E+00 /)  !SIX --> SVIII
      RR_rates_alt_sulfur(8,:)  = (/ 5.50E-11, 0.711E+00 /)  !SX --> SIX
      RR_rates_alt_sulfur(7,:)  = (/ 7.40E-11, 0.716E+00 /)  !SXI --> SX
      RR_rates_alt_sulfur(6,:)  = (/ 9.20E-11, 0.714E+00 /)  !SXII --> SXI
      RR_rates_alt_sulfur(5,:)  = (/ 1.40E-10, 0.755E+00 /)  !SXIII --> SXII
      RR_rates_alt_sulfur(4,:)  = (/ 2.00E-10, 0.806E+00 /)  !SXIV --> SXIII
      RR_rates_alt_sulfur(3,:)  = (/ 2.91E-10, 0.840E+00 /)  !SXV --> SXIV
      RR_rates_alt_sulfur(2,:)  = (/ 4.30E-10, 0.807E+00 /)  !SXVI --> SXV
      RR_rates_alt_sulfur(1,:)  = (/ 0.00E+00, 0.00E+00 /)   !SXVII --> SXVI

      DR_rates_alt_sulfur(16,:) = (/ 1.62E-03, 0.00E+00, 1.25E+05, 1.00E+05 /)
      DR_rates_alt_sulfur(15,:) = (/ 1.09E-02, 1.20E-02, 1.92E+05, 1.80E+04 /)
      DR_rates_alt_sulfur(14,:) = (/ 3.35E-02, 6.59E-02, 1.89E+05, 1.59E+05 /)
      DR_rates_alt_sulfur(13,:) = (/ 3.14E-02, 6.89E-02, 1.68E+05, 8.04E+04 /)
      DR_rates_alt_sulfur(12,:) = (/ 1.27E-02, 1.87E-01, 1.38E+05, 1.71E+05 /)
      DR_rates_alt_sulfur(11,:) = (/ 1.47E-02, 1.29E-01, 1.80E+06, 1.75E+06 /)
      DR_rates_alt_sulfur(10,:) = (/ 1.34E-02, 1.04E+00, 6.90E+05, 2.15E+06 /)
      DR_rates_alt_sulfur(9,:)  = (/ 2.38E-02, 1.12E+00, 5.84E+05, 2.59E+06 /)
      DR_rates_alt_sulfur(8,:)  = (/ 3.19E-02, 1.40E+00, 5.17E+05, 2.91E+06 /)
      DR_rates_alt_sulfur(7,:)  = (/ 7.13E-02, 1.00E+00, 6.66E+05, 2.32E+06 /)
      DR_rates_alt_sulfur(6,:)  = (/ 8.00E-02, 5.55E-01, 6.00E+05, 2.41E+06 /)
      DR_rates_alt_sulfur(5,:)  = (/ 7.96E-02, 1.63E+00, 5.09E+05, 6.37E+06 /)
      DR_rates_alt_sulfur(4,:)  = (/ 1.34E-02, 3.04E-01, 2.91E+05, 1.04E+06 /)
      DR_rates_alt_sulfur(3,:)  = (/ 4.02E-01, 2.98E-01, 2.41E+07, 4.67E+06 /)
      DR_rates_alt_sulfur(2,:)  = (/ 1.45E-01, 2.81E-01, 2.54E+07, 5.30E+06 /)
      DR_rates_alt_sulfur(1,:)  = (/ 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00 /)

      !!!!!!!!
      !!!!!!!! Iron
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_iron(1,:)  = (/ 2.275E-09, 0.7481E+00, 1.836E+03, 4.736E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(2,:)  = (/ 9.983E-10, 0.6754E+00, 6.651E+03, 4.017E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(3,:)  = (/ 4.458E-10, 0.5802E+00, 2.155E+04, 1.701E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(4,:)  = (/ 1.186E-09, 0.6713E+00, 3.253E+03, 9.392E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(5,:)  = (/ 3.322E-09, 0.7264E+00, 4.563E+02, 5.746E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(6,:)  = (/ 2.199E-09, 0.7118E+00, 7.810E+02, 5.946E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(7,:)  = (/ 1.659E-09, 0.6958E+00, 1.061E+03, 6.253E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(8,:)  = (/ 1.135E-09, 0.6705E+00, 1.691E+03, 6.809E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(9,:)  = (/ 7.556E-10, 0.6351E+00, 2.800E+03, 7.742E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(10,:) = (/ 4.791E-10, 0.5823E+00, 4.967E+03, 9.535E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(11,:) = (/ 2.034E-10, 0.4548E+00, 1.751E+04, 1.579E+08, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(12,:) = (/ 3.133E-10, 0.5507E+00, 6.295E+03, 9.035E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(13,:) = (/ 5.398E-10, 0.6295E+00, 1.881E+03, 5.429E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(14,:) = (/ 1.121E-09, 0.6984E+00, 4.071E+02, 4.264E+07, 0.0000E+00, 0.000E+00 /)
      RR_rates_iron(15,:) = (/ 1.984E-09, 0.7101E+00, 1.158E+02, 4.400E+07, 0.0000E+00, 0.000E+00 /) 

      ! Dielectronic recombination rates
      DR_rates_c_iron(1,:)   = (/ 0.000E+00, 0.000E+00,  0.000E+00, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(2,:)   = (/ 1.984E-01, 2.676E-01, -2.293E-03, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(3,:)   = (/ 2.676E-01, 4.097E-01,  2.990E-02, 0.000E+00,  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(4,:)   = (/ 7.882E-03, 1.636E-02,  4.868E-02, 4.230E-02,  4.151E-01, 5.339E-01, 4.544E-03, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(5,:)   = (/ 2.325E-03, 1.017E-02,  3.572E-02, 9.882E-02,  1.156E-01, 5.792E-01, 3.344E-01, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(6,:)   = (/ 8.382E-03, 7.897E-03,  3.157E-02, 1.159E-01,  3.919E-01, 1.017E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(7,:)   = (/ 2.565E-03, 1.685E-02,  1.827E-02, 6.957E-02,  3.254E-01, 5.101E-01, 7.325E-01, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(8,:)   = (/ 2.106E-03, 6.569E-03,  1.532E-02, 3.799E-02,  7.669E-02, 6.701E-01, 1.298E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(9,:)   = (/ 2.033E-04, 1.159E-03,  5.567E-03, 5.482E-02,  3.370E-01, 1.518E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(10,:)  = (/ 8.207E-05, 2.766E-04,  1.897E-03, 2.842E-02,  4.022E-01, 1.434E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(11,:)  = (/ 6.342E-04, 8.350E-02,  1.045E+00, 3.663E-01, -3.955E-02, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(12,:)  = (/ 7.676E-04, 5.587E-03,  1.152E-01, 4.929E-02,  7.274E-01, 7.347E-03, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(13,:)  = (/ 5.636E-04, 7.860E-03,  5.063E-02, 1.753E-01,  1.209E-01, 1.934E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(14,:)  = (/ 1.753E-03, 1.038E-02,  2.573E-02, 1.189E-01,  1.070E-01, 3.080E-02, 6.324E-04, 0.000E+00, 0.000E+00 /)
      DR_rates_c_iron(15,:)  = (/ 4.469E-03, 8.538E-03,  1.741E-02, 1.630E-01,  8.680E-02, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)

      DR_rates_e_iron(1,:)   = (/ 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(2,:)   = (/ 5.552E+07, 7.475E+07, 1.236E+08, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(3,:)   = (/ 5.394E+07, 6.854E+07, 9.651E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(4,:)   = (/ 5.977E+04, 2.074E+05, 6.170E+05, 3.789E+06, 1.086E+07, 6.363E+07, 1.599E+08, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(5,:)   = (/ 1.376E+04, 8.251E+04, 2.794E+05, 9.378E+05, 4.688E+06, 1.106E+07, 6.543E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(6,:)   = (/ 5.297E+03, 5.829E+04, 2.454E+05, 9.663E+05, 5.580E+06, 1.110E+07, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(7,:)   = (/ 6.553E+03, 4.351E+04, 2.059E+05, 9.250E+05, 4.610E+06, 1.019E+07, 1.019E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(8,:)   = (/ 4.463E+03, 3.545E+04, 1.944E+05, 6.148E+05, 1.635E+06, 6.100E+06, 1.026E+07, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(9,:)   = (/ 3.502E+03, 3.557E+04, 2.177E+05, 1.078E+06, 4.515E+06, 9.050E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(10,:)  = (/ 6.365E+03, 3.842E+04, 2.002E+05, 1.150E+06, 4.736E+06, 8.891E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(11,:)  = (/ 2.770E+06, 3.978E+06, 7.052E+06, 1.300E+07, 3.579E+07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(12,:)  = (/ 2.935E+04, 8.158E+04, 3.591E+05, 1.735E+06, 7.545E+06, 4.634E+07, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(13,:)  = (/ 3.628E+03, 2.489E+04, 1.405E+05, 5.133E+05, 5.018E+06, 8.689E+06, 0.000E+00, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(14,:)  = (/ 4.287E+03, 1.679E+04, 9.992E+04, 3.841E+05, 6.758E+05, 1.476E+06, 3.334E+08, 0.000E+00, 0.000E+00 /)
      DR_rates_e_iron(15,:)  = (/ 2.462E+03, 1.261E+04, 9.330E+04, 4.887E+05, 1.312E+06, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)

      ! radiative recombination data from https://www.pa.uky.edu/~verner/rec.html
      ! The index number is the number of electons on the atom before the recombination process + 1
      RR_rates_alt_iron(26,:) = (/ 1.42E-13, 0.891E+00 /)
      RR_rates_alt_iron(25,:) = (/ 1.02E-12, 0.843E+00 /)
      RR_rates_alt_iron(24,:) = (/ 3.32E-12, 0.746E+00 /)
      RR_rates_alt_iron(23,:) = (/ 7.80E-12, 0.682E+00 /)
      RR_rates_alt_iron(22,:) = (/ 1.51E-11, 0.699E+00 /)
      RR_rates_alt_iron(21,:) = (/ 2.62E-11, 0.728E+00 /)
      RR_rates_alt_iron(20,:) = (/ 4.12E-11, 0.759E+00 /)
      RR_rates_alt_iron(19,:) = (/ 6.05E-11, 0.790E+00 /)
      RR_rates_alt_iron(18,:) = (/ 8.13E-11, 0.810E+00 /)
      RR_rates_alt_iron(17,:) = (/ 1.09E-10, 0.829E+00 /)
      RR_rates_alt_iron(16,:) = (/ 1.33E-10, 0.828E+00 /)
      RR_rates_alt_iron(15,:) = (/ 1.64E-10, 0.834E+00 /)
      RR_rates_alt_iron(14,:) = (/ 2.00E-10, 0.836E+00 /)
      RR_rates_alt_iron(13,:) = (/ 2.41E-10, 0.840E+00 /)
      RR_rates_alt_iron(12,:) = (/ 2.89E-10, 0.846E+00 /)
      RR_rates_alt_iron(11,:) = (/ 3.42E-10, 0.850E+00 /)
      RR_rates_alt_iron(10,:) = (/ 3.87E-10, 0.836E+00 /)
      RR_rates_alt_iron(9,:)  = (/ 4.52E-10, 0.824E+00 /)
      RR_rates_alt_iron(8,:)  = (/ 5.25E-10, 0.816E+00 /)
      RR_rates_alt_iron(7,:)  = (/ 6.07E-10, 0.811E+00 /)
      RR_rates_alt_iron(6,:)  = (/ 6.98E-10, 0.808E+00 /)
      RR_rates_alt_iron(5,:)  = (/ 7.72E-10, 0.800E+00 /)
      RR_rates_alt_iron(4,:)  = (/ 1.15E-09, 0.852E+00 /)
      RR_rates_alt_iron(3,:)  = (/ 1.58E-09, 0.875E+00 /)
      RR_rates_alt_iron(2,:)  = (/ 1.40E-09, 0.787E+00 /)
      RR_rates_alt_iron(1,:)  =  (/ 0.00E+00, 0.00E+00 /)  

      DR_rates_alt_iron(26,:) = (/ 5.120E+00, 1.29E+01, 0.00E+00, 0.00E+00, 2.20E-04, 1.00E-04, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(25,:) = (/ 1.670E+01, 3.14E+01, 0.00E+00, 0.00E+00, 2.30E-03, 2.70E-03, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(24,:) = (/ 2.860E+01, 5.21E+01, 0.00E+00, 0.00E+00, 1.50E-02, 4.70E-03, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(23,:) = (/ 3.730E+01, 6.74E+01, 0.00E+00, 0.00E+00, 3.80E-02, 1.60E-02, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(22,:) = (/ 5.420E+01, 1.00E+02, 0.00E+00, 0.00E+00, 8.00E-02, 2.40E-02, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(21,:) = (/ 4.550E+01, 3.60E+02, 0.00E+00, 0.00E+00, 9.20E-02, 4.10E-02, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(20,:) = (/ 6.670E+01, 1.23E+02, 0.00E+00, 0.00E+00, 1.60E-01, 3.60E-02, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(19,:) = (/ 6.610E+01, 1.29E+02, 0.00E+00, 0.00E+00, 1.80E-01, 7.00E-02, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(18,:) = (/ 2.160E+01, 1.36E+02, 0.00E+00, 0.00E+00, 1.40E-01, 2.60E-01, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(17,:) = (/ 2.220E+01, 1.44E+02, 0.00E+00, 0.00E+00, 1.00E-01, 2.80E-01, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(16,:) = (/ 5.960E+01, 3.62E+02, 0.00E+00, 0.00E+00, 2.25E-01, 2.31E-01, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(15,:) = (/ 7.500E+01, 2.05E+02, 0.00E+00, 0.00E+00, 2.40E-01, 1.70E-01, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(14,:) = (/ 3.630E+01, 1.93E+02, 0.00E+00, 0.00E+00, 2.60E-01, 1.60E-01, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(13,:) = (/ 3.940E+01, 1.98E+02, 0.00E+00, 0.00E+00, 1.90E-01, 9.00E-02, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(12,:) = (/ 2.460E+01, 2.48E+02, 5.60E+02, 0.00E+00, 1.20E-01, 1.20E-01, 6.00E-01, 0.00E+00 /)
      DR_rates_alt_iron(11,:) = (/ 5.600E+02, 0.00E+00, 0.00E+00, 0.00E+00, 1.23E+00, 0.00E+00, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(10,:) = (/ 2.250E+01, 1.17E+02, 3.41E+02, 6.83E+02, 2.53E-03, 3.36E-02, 1.81E-01, 1.92E+00 /)
      DR_rates_alt_iron(9,:)  = (/ 1.620E+01, 9.60E+01, 3.30E+02, 7.29E+02, 5.67E-03, 7.82E-02, 3.18E-02, 1.26E+00 /)
      DR_rates_alt_iron(8,:)  = (/ 2.370E+01, 8.51E+01, 3.29E+02, 7.87E+02, 1.60E-02, 7.17E-02, 9.06E-02, 7.39E-01 /)
      DR_rates_alt_iron(7,:)  = (/ 1.320E+01, 6.66E+01, 2.97E+02, 7.14E+02, 1.85E-02, 9.53E-02, 7.90E-02, 1.23E+00 /)
      DR_rates_alt_iron(6,:)  = (/ 3.910E+01, 8.03E+01, 3.92E+02, 9.19E+02, 9.20E-04, 1.29E-01, 1.92E-01, 9.12E-01 /)
      DR_rates_alt_iron(5,:)  = (/ 7.320E+01, 3.16E+02, 8.77E+02, 0.00E+00, 1.31E-01, 8.49E-02, 6.13E-01, 0.00E+00 /)
      DR_rates_alt_iron(4,:)  = (/ 1.000E-01, 3.62E+01, 3.06E+02, 9.28E+02, 1.10E-02, 4.88E-02, 8.01E-02, 5.29E-01 /)
      DR_rates_alt_iron(3,:)  = (/ 4.625E+03, 6.00E+03, 0.00E+00, 0.00E+00, 2.56E-01, 4.52E-01, 0.00E+00, 0.00E+00 /)
      DR_rates_alt_iron(2,:)  = (/ 5.300E+03, 0.00E+00, 0.00E+00, 0.00E+00, 4.30E-01, 0.00E+00, 0.00E+00, 0.00E+00 /)

      !!!!!!!!
      !!!!!!!! Neon
      !!!!!!!!
      ! Radiative recombination rates
      RR_rates_neon(1,:)  = (/ 8.278E-10,  0.7470E+00,  2.991E+02,  7.006E+07,  0.0000E+00,  0.000E+00 /)
      RR_rates_neon(2,:)  = (/ 3.415E-10,  0.6706E+00,  9.552E+02,  6.778E+07,  0.0000E+00,  0.000E+00 /)
      RR_rates_neon(3,:)  = (/ 1.186E-10,  0.5354E+00,  3.647E+03,  3.365E+07,  0.0000E+00,  0.000E+00 /)
      RR_rates_neon(4,:)  = (/ 2.755E-10,  0.6586E+00,  5.102E+02,  1.535E+07,  0.0000E+00,  0.000E+00 /)
      RR_rates_neon(5,:)  = (/ 2.557E-09,  0.7601E+00,  6.293E+00,  8.091E+06,  0.0000E+00,  0.000E+00 /)
      RR_rates_neon(6,:)  = (/ 1.127E-09,  0.7556E+00,  1.311E+01,  8.047E+06,  0.0250E+00,  2.771E+05 /)
      RR_rates_neon(7,:)  = (/ 1.861E-09,  0.7593E+00,  2.504E+00,  8.037E+06,  0.0406E+00,  3.255E+05 /)
      RR_rates_neon(8,:)  = (/ 8.321E-10,  0.7254E+00,  3.332E+00,  8.696E+06,  0.0921E+00,  3.044E+05 /)
      RR_rates_neon(9,:)  = (/ 1.773E-10,  0.6434E+00,  9.924E+00,  8.878E+06,  0.2205E+00,  2.292E+05 /)
      RR_rates_neon(10,:) = (/ 1.295E-11,  0.3556E+00,  6.739E+01,  7.563E+06,  0.6472E+00,  1.598E+05 /)

      ! Dielectronic recombination rates
      DR_rates_c_neon(1,:)   = (/ 0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(2,:)   = (/ 1.183E-02,  9.011E-02,  1.828E-03,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(3,:)   = (/ 1.552E-02,  9.008E-02,  1.182E-02,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(4,:)   = (/ 2.399E-04,  3.532E-04,  8.928E-03,  5.427E-03,  5.342E-03,  3.981E-02,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(5,:)   = (/ 8.460E-05,  1.817E-04,  1.176E-03,  1.397E-02,  5.566E-03,  6.229E-03,  1.031E-02,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(6,:)   = (/ 5.653E-06,  4.344E-05,  1.086E-04,  5.980E-04,  1.457E-02,  1.601E-02,  5.365E-04,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(7,:)   = (/ 2.922E-06,  7.144E-06,  2.836E-05,  9.820E-05,  8.379E-03,  1.009E-02,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(8,:)   = (/ 2.763E-06,  1.053E-05,  4.453E-05,  6.244E-03,  3.146E-04,  4.465E-03,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(9,:)   = (/ 2.980E-08,  1.257E-07,  1.122E-06,  2.626E-03,  8.802E-04,  1.231E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_c_neon(10,:)  = (/ 4.152E-09,  4.656E-09,  1.310E-08,  1.417E-09,  7.968E-04,  1.271E-05,  0.000E+00,  0.000E+00,  0.000E+00 /)

      DR_rates_e_neon(1,:)   = (/ 0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(2,:)   = (/ 8.405E+06,  1.111E+07,  1.812E+07,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(3,:)   = (/ 7.845E+06,  9.803E+06,  1.209E+07,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(4,:)   = (/ 2.536E+04,  4.800E+04,  1.770E+05,  1.030E+06,  1.859E+06,  9.743E+06,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(5,:)   = (/ 1.049E+03,  3.829E+03,  6.133E+04,  2.568E+05,  4.600E+05,  1.324E+06,  9.353E+06,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(6,:)   = (/ 6.280E+02,  2.812E+03,  1.324E+04,  8.064E+04,  3.052E+05,  1.032E+06,  2.388E+06,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(7,:)   = (/ 2.050E+02,  2.205E+03,  9.271E+03,  4.988E+04,  2.904E+05,  8.782E+05,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(8,:)   = (/ 6.393E+02,  1.499E+03,  3.227E+04,  2.561E+05,  4.505E+05,  7.934E+05,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(9,:)   = (/ 4.579E+01,  4.753E+02,  1.481E+04,  2.810E+05,  4.763E+05,  4.677E+08,  0.000E+00,  0.000E+00,  0.000E+00 /)
      DR_rates_e_neon(10,:)  = (/ 2.689E+01,  2.021E+02,  7.200E+02,  4.892E+04,  3.144E+05,  6.738E+05,  0.000E+00,  0.000E+00,  0.000E+00 /)

   END SUBROUTINE init_metal_atomic_data_tables
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE init_coolrates_tables(aexp)

! Initialise the cooling rates tables.
!-------------------------------------------------------------------------
      implicit none
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      real(dp) :: aexp, T
      integer :: myid, ncpu, ierr, iT
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
#endif
#ifdef WITHOUTMPI
      myid = 0
      ncpu = 1
#endif

      ! Initialise the table lookup temperatures -----------------------------
      do iT = 1, nbinT
         T_lookup(iT) = log10(Tmin) + (dble(iT) - 1d0)/(dble(nbinT) - 1d0) &
                        *(log10(Tmax) - log10(Tmin))
      end do
      dlogTinv = dble(nbinT - 1)/(T_lookup(nbinT) - T_lookup(1)) ! (space)^-1
      hTable = 1d0/dlogTinv                             !
      h2Table = hTable*hTable                           ! Constants for table
      h3Table = h2Table*hTable                          ! interpolation

      one_over_lnTen = 1.d0/log(10d0)
      one_over_hTable = 1.d0/hTable
      one_over_h2Table = 1.d0/h2Table
      three_over_h2Table = 3.d0/h2Table
      two_over_h3Table = 2.d0/h3Table

      do iT = myid + 1, nbinT, ncpu ! Loop over TK and assign rates
         call comp_table_rates(iT, aexp)
      end do ! end TK loop

      ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
      call mpi_distribute_coolrates_table(tbl_alphaA_HII)
      call mpi_distribute_coolrates_table(tbl_alphaA_HeII)
      call mpi_distribute_coolrates_table(tbl_alphaA_HeIII)

      call mpi_distribute_coolrates_table(tbl_alphaB_HII)
      call mpi_distribute_coolrates_table(tbl_alphaB_HeII)
      call mpi_distribute_coolrates_table(tbl_alphaB_HeIII)

      call mpi_distribute_coolrates_table(tbl_beta_HI)
      call mpi_distribute_coolrates_table(tbl_beta_HeI)
      call mpi_distribute_coolrates_table(tbl_beta_HeII)

      call mpi_distribute_coolrates_table(tbl_cr_ci_HI)
      call mpi_distribute_coolrates_table(tbl_cr_ci_HeI)
      call mpi_distribute_coolrates_table(tbl_cr_ci_HeII)

      call mpi_distribute_coolrates_table(tbl_cr_ce_HI)
      call mpi_distribute_coolrates_table(tbl_cr_ce_HeI)
      call mpi_distribute_coolrates_table(tbl_cr_ce_HeII)

      call mpi_distribute_coolrates_table(tbl_cr_r_HII)
      call mpi_distribute_coolrates_table(tbl_cr_r_HeII)
      call mpi_distribute_coolrates_table(tbl_cr_r_HeIII)

      call mpi_distribute_coolrates_table(tbl_cr_bre)
      call mpi_distribute_coolrates_table(tbl_cr_com)
      call mpi_distribute_coolrates_table(tbl_cr_die)
#endif

      if (myid == 0) print *, 'Coolrates tables initialised '
901   format(20(1pe12.3))

   END SUBROUTINE init_coolrates_tables

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE update_coolrates_tables(aexp)
! Update cooling rates lookup tables which depend on aexp
!-------------------------------------------------------------------------
      implicit none
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif
      real(dp) :: aexp
      integer:: myid, ncpu, ierr, iT
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
#endif
#ifdef WITHOUTMPI
      myid = 0
      ncpu = 1
#endif
      tbl_cr_com%rates = 0d0; tbl_cr_com%primes = 0d0
      do iT = myid + 1, nbinT, ncpu ! Loop over TK and assign rates
         call update_table_rates(iT, aexp)
      end do ! end TK loop

      ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
      call mpi_distribute_coolrates_table(tbl_cr_com)
#endif

      if (myid == 0) print *, 'Coolrates table updated'
   END SUBROUTINE update_coolrates_tables

#ifndef WITHOUTMPI
!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE mpi_distribute_coolrates_table(table)
! Distribute table between all cpus, assuming table contains only partial
! entries on each cpu, but the whole table is acquired by summing those
! partial tables
!-------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      type(coolrates_table)::table
      real(dp), dimension(:), allocatable :: table_mpi_sum
      integer::ierr
!-------------------------------------------------------------------------
      allocate (table_mpi_sum(nbinT))
      call MPI_ALLREDUCE(table%rates, table_mpi_sum &
                         , nbinT, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      table%rates = table_mpi_sum
      call MPI_ALLREDUCE(table%primes, table_mpi_sum &
                         , nbinT, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      table%primes = table_mpi_sum
      deallocate (table_mpi_sum)

   END SUBROUTINE mpi_distribute_coolrates_table
#endif

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE comp_table_rates(iT, aexp)
! Fill in index iTK in all rates tables.
!-------------------------------------------------------------------------
      use rt_parameters, only: rt_OTSA, use_cloudy_prim_rates
      implicit none
      integer::iT
      real(dp)::aexp, T, Ta, T5, lambda, f, hf, laHII, laHeII, laHeIII
      real(dp), parameter::kb = 1.3806d-16        ! Boltzmann constant [ergs K-1]
!-------------------------------------------------------------------------
      ! Rates are stored in non-log, while temperature derivatives (primes)
      ! are stored in dRate/dlogT (= dRate/dT * T * ln(10)).

      ! The log-log primes are just the normal primes times T/rate,
      ! i.e. dlogL/dlogT = T/L dL/dT

      T = 10d0**T_lookup(iT)

      if (.not.use_cloudy_prim_rates) then
         ! Case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
         lambda = 315614./T                                ! 2.d0 * 157807.d0 / T
         f = 1.d0 + (lambda/0.522)**0.47
         tbl_alphaA_HII%rates(iT) = 1.269d-13*lambda**1.503/f**1.923
         tbl_alphaA_HII%primes(iT) = (0.90381*(f - 1.)/f - 1.503) &
                                    *log(10d0)*tbl_alphaA_HII%rates(iT)

         ! Case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
         lambda = 570670./T
         tbl_alphaA_HeII%rates(iT) = 3.d-14*lambda**0.654
         tbl_alphaA_HeII%primes(iT) = -0.654 &
                                    *log(10d0)*tbl_alphaA_HeII%rates(iT)

         ! Case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
         lambda = 1263030./T
         f = 1.d0 + (lambda/0.522)**0.47
         tbl_alphaA_HeIII%rates(iT) = 2.538d-13*lambda**1.503/f**1.923
         tbl_alphaA_HeIII%primes(iT) = (0.90381*(f - 1.)/f - 1.503) &
                                       *log(10d0)*tbl_alphaA_HeIII%rates(iT)

         ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
         lambda = 315614./T
         f = 1.d0 + (lambda/2.74)**0.407
         tbl_alphaB_HII%rates(iT) = 2.753d-14*lambda**1.5/f**2.242
         tbl_alphaB_HII%primes(iT) = (0.912494*(f - 1.)/f - 1.5) &
                                    *log(10d0)*tbl_alphaB_HII%rates(iT)

         ! Case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
         lambda = 570670./T
         tbl_alphaB_HeII%rates(iT) = 1.26d-14*lambda**0.75
         tbl_alphaB_HeII%primes(iT) = -0.75 &
                                    *log(10d0)*tbl_alphaB_HeII%rates(iT)

         ! Case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
         lambda = 1263030./T
         f = 1.d0 + (lambda/2.74)**0.407
         tbl_alphaB_HeIII%rates(iT) = 5.506d-14*lambda**1.5/f**2.242
         tbl_alphaB_HeIII%primes(iT) = (0.912494*(f - 1.)/f - 1.5) &
                                       *log(10d0)*tbl_alphaB_HeIII%rates(iT)
      else ! Taken from CLOUDY (no difference between case A and case B for now)
         tbl_alphaA_HII%rates(iT) = comp_Alpha_hydrogen_badnell(T)
         tbl_alphaB_HII%rates(iT) = comp_Alpha_hydrogen_badnell(T)

         tbl_alphaA_HII%primes(iT) = comp_Alpha_hydrogen_badnell(1.001d0*T)
         tbl_alphaA_HII%primes(iT) = (tbl_alphaA_HII%primes(iT) - tbl_alphaA_HII%rates(iT)) / (1.001d0*T - T)
         tbl_alphaB_HII%primes(iT) = tbl_alphaA_HII%primes(iT)

         tbl_alphaA_HeII%rates(iT) = comp_Alpha_helium_badnell(T,1)
         tbl_alphaB_HeII%rates(iT) = comp_Alpha_helium_badnell(T,1)

         tbl_alphaA_HeII%primes(iT) = comp_Alpha_helium_badnell(1.001d0*T,1)
         tbl_alphaA_HeII%primes(iT) = (tbl_alphaA_HeII%primes(iT) - tbl_alphaA_HeII%rates(iT)) / (1.001d0*T - T)
         tbl_alphaB_HeII%primes(iT) = tbl_alphaA_HeII%primes(iT)

         tbl_alphaA_HeIII%rates(iT) = comp_Alpha_helium_badnell(T,0)
         tbl_alphaB_HeIII%rates(iT) = comp_Alpha_helium_badnell(T,0)

         tbl_alphaA_HeIII%primes(iT) = comp_Alpha_helium_badnell(1.001d0*T,0)
         tbl_alphaA_HeIII%primes(iT) = (tbl_alphaA_HeIII%primes(iT) - tbl_alphaA_HeIII%rates(iT)) / (1.001d0*T - T)
         tbl_alphaB_HeIII%primes(iT) = tbl_alphaA_HeIII%primes(iT)
      endif

      if (.not.use_cloudy_prim_rates) then !STANDARD RAMSES-RT
         ! Collisional ionization rate [cm3 s-1] of HI (Maselli&'03)-------------
         T5 = T/1d5
         f = 1d0 + sqrt(T5); hf = 0.5d0/f
         tbl_beta_HI%rates(iT) = 5.85d-11*sqrt(T)/f*exp(-157809.1d0/T)
         tbl_beta_HI%primes(iT) = (hf + 157809.1d0/T) &
                                 *log(10d0)*tbl_beta_HI%rates(iT)

         ! Collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)------------
         tbl_beta_HeI%rates(iT) = 2.38d-11*sqrt(T)/f*exp(-285335.4d0/T)
         tbl_beta_HeI%primes(iT) = (hf + 285335.4d0/T) &
                                 *log(10d0)*tbl_beta_HeI%rates(iT)

         ! Collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)-----------
         tbl_beta_HeII%rates(iT) = 5.68d-12*sqrt(T)/f*exp(-631515.d0/T)
         tbl_beta_HeII%primes(iT) = (hf + 631515.d0/T) &
                                    *log(10d0)*tbl_beta_HeII%rates(iT)
      else !TAKEN FROM CLOUDY
         tbl_beta_HI%rates(iT) = comp_Beta_hydrogen_cloudy(T)
         tbl_beta_HI%primes(iT) = comp_Beta_hydrogen_cloudy(1.001d0*T)
         tbl_beta_HI%primes(iT) = (tbl_beta_HI%primes(iT) - tbl_beta_HI%rates(iT)) / (1.001d0*T - T)

         tbl_beta_HeI%rates(iT) = comp_Beta_helium_cloudy(T,1)
         tbl_beta_HeI%primes(iT) = comp_Beta_helium_cloudy(1.001d0*T,1)
         tbl_beta_HeI%primes(iT) = (tbl_beta_HeI%primes(iT) - tbl_beta_HeI%rates(iT)) / (1.001d0*T - T)

         tbl_beta_HeII%rates(iT) = comp_Beta_helium_cloudy(T,2)
         tbl_beta_HeII%primes(iT) = comp_Beta_helium_cloudy(1.001d0*T,2)
         tbl_beta_HeII%primes(iT) = (tbl_beta_HeII%primes(iT) - tbl_beta_HeII%rates(iT)) / (1.001d0*T - T)
      endif

      ! BEGIN COOLING RATES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      T5 = T/1d5
      f = 1d0 + sqrt(T5); hf = 0.5d0/f

      ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
      tbl_cr_ci_HI%rates(iT) = 1.27d-21*sqrt(T)/f*exp(-157809.1/T)
      tbl_cr_ci_HI%primes(iT) = (hf + 157809.1/T) &
                                *log(10d0)*tbl_cr_ci_HI%rates(iT)

      tbl_cr_ci_HeI%rates(iT) = 9.38d-22*sqrt(T)/f*exp(-285335.4/T)
      tbl_cr_ci_HeI%primes(iT) = (hf + 285335.4/T) &
                                 *log(10d0)*tbl_cr_ci_HeI%rates(iT)

      tbl_cr_ci_HeII%rates(iT) = 4.95d-22*sqrt(T)/f*exp(-631515./T)
      tbl_cr_ci_HeII%primes(iT) = (hf + 631515.0/T) &
                                  *log(10d0)*tbl_cr_ci_HeII%rates(iT)

      ! Collisional excitation cooling from Cen'92
      tbl_cr_ce_HI%rates(iT) = 7.5d-19/f*exp(-118348./T)
      tbl_cr_ce_HI%primes(iT) = (118348./T - 0.5d0*sqrt(T5)/f) &
                                *log(10d0)*tbl_cr_ce_HI%rates(iT)

      tbl_cr_ce_HeI%rates(iT) = 9.10d-27*T**(-0.1687)/f*exp(-13179./T)
      tbl_cr_ce_HeI%primes(iT) = (13179./T - 0.1687 - 0.5d0*sqrt(T5)/f) &
                                 *log(10d0)*tbl_cr_ce_HeI%rates(iT)

      tbl_cr_ce_HeII%rates(iT) = 5.54d-17*T**(-0.397)/f*exp(-473638./T)
      tbl_cr_ce_HeII%primes(iT) = (473638./T - 0.397 - 0.5d0*sqrt(T5)/f) &
                                  *log(10d0)*tbl_cr_ce_HeII%rates(iT)

      ! Recombination Cooling (Hui&Gnedin'97)
      laHII = 315614./T
      laHeII = 570670./T
      laHeIII = 1263030./T
      if (.not. rt_otsa) then ! Case A
         f = 1.d0 + (laHII/0.541)**0.502
         tbl_cr_r_HII%rates(iT) = 1.778d-29*laHII**1.965/f**2.697*T
         tbl_cr_r_HII%primes(iT) = (-0.965 + 1.35389*(f - 1.)/f) &
                                   *log(10d0)*tbl_cr_r_HII%rates(iT)

         tbl_cr_r_HeII%rates(iT) = 3.d-14*laHeII**0.654*kb*T
         tbl_cr_r_HeII%primes(iT) = 0.346 &
                                    *log(10d0)*tbl_cr_r_HeII%rates(iT)

         f = 1.d0 + (laHeIII/0.541)**0.502
         tbl_cr_r_HeIII%rates(iT) = 14.224d-29*laHeIII**1.965/f**2.697*T
         tbl_cr_r_HeIII%primes(iT) = (-0.965 + 1.35389*(f - 1.)/f) &
                                     *log(10d0)*tbl_cr_r_HeIII%rates(iT)
      else ! Case B
         f = 1.d0 + (laHII/2.25)**0.376
         tbl_cr_r_HII%rates(iT) = 3.435d-30*laHII**1.97/f**3.72*T
         tbl_cr_r_HII%primes(iT) = (-0.97 + 1.39827*(f - 1.)/f) &
                                   *log(10d0)*tbl_cr_r_HII%rates(iT)

         tbl_cr_r_HeII%rates(iT) = 1.26d-14*laHeII**0.75*kb*T
         tbl_cr_r_HeII%primes(iT) = 0.25 &
                                    *log(10d0)*tbl_cr_r_HeII%rates(iT)

         f = 1.d0 + (laHeIII/2.25)**0.376
         tbl_cr_r_HeIII%rates(iT) = 27.48d-30*laHeIII**1.97/f**3.72*T
         tbl_cr_r_HeIII%primes(iT) = (-0.97 + 1.39827*(f - 1.)/f) &
                                     *log(10d0)*tbl_cr_r_HeIII%rates(iT)
      end if

      ! Bremsstrahlung from Osterbrock & Ferland 2006
      tbl_cr_bre%rates(iT) = 1.42d-27*1.5*sqrt(T)
      tbl_cr_bre%primes(iT) = 0.5 &
                              *log(10d0)*tbl_cr_bre%rates(iT)

      ! Compton Cooling from Haimann et al. 96, via Maselli et al.
      ! Need to make sure this is done whenever the redshift changes!
      Ta = 2.727/aexp
      tbl_cr_com%rates(iT) = 1.017d-37*Ta**4*(T - Ta)
      tbl_cr_com%primes(iT) = T/(T - Ta) &
                              *log(10d0)*tbl_cr_com%rates(iT)

      ! Dielectronic recombination cooling, from Black 1981
      f = 1.24d-13*T**(-1.5d0)*exp(-470000.d0/T)
      tbl_cr_die%rates(iT) = f*(1.D0 + 0.3d0*exp(-94000.d0/T))
      tbl_cr_die%primes(iT) = 0d0
      if (tbl_cr_die%rates(iT) .gt. 0d0) then ! Can simplify w algebra
         tbl_cr_die%primes(iT) = (tbl_cr_die%rates(iT)*(564000.-1.5*T) &
                                  - f*94000.)/T**2*T*log(10d0)
      end if

   END SUBROUTINE comp_table_rates

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   SUBROUTINE update_table_rates(iT, aexp)
! Update index iTK in all compton cooling rates tables due to change in
! aexp
!-------------------------------------------------------------------------
      implicit none
      integer::iT
      real(dp)::aexp, T, Ta
!-------------------------------------------------------------------------
      ! Rates are stored in log, while temperature derivatives (primes) are
      ! stored in non-log
      T = 10d0**T_lookup(iT)
      ! Compton Cooling from Haimann et al. 96, via Maselli et al.
      Ta = 2.727/aexp
      tbl_cr_com%rates(iT) = 1.017d-37*Ta**4*(T - Ta)
      tbl_cr_com%primes(iT) = T/(T - Ta) &
                              *log(10d0)*tbl_cr_com%rates(iT)

   END SUBROUTINE update_table_rates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   FUNCTION inp_coolrates_table(rates_table, T, retPrime)
! Returns TABULATED rate value from given table
! rates        => Rates (and primes) table to interpolate
! T            => Temperature [K]
! retPrime     <= Temperature derivative of rate at T (optional).
!-------------------------------------------------------------------------
      implicit none
      type(coolrates_table)::rates_table
      !real(dp),dimension(nbinTK)::rates, ratesPrime
      real(dp), intent(in)::T
      real(dp), optional::retPrime
      real(dp)::inp_coolrates_table
      integer, save:: iT = 1
      real(dp), save:: facT, yy, yy2, yy3, fa, fb, fprimea, fprimeb, Tlast = -1
      real(dp), save:: alpha, beta, gamma
      logical, save::extrap
!-------------------------------------------------------------------------
      if (.not. (T .eq. Tlast)) then    ! Reuse index if same T from last call
         ! Log of T, snapped to table at the lower boundary, but allowed
         ! to go above upper boundary, in which case we use extrapolation:
         facT = MAX(log10(T), T_lookup(1))
         extrap = .false.
         if (facT .gt. T_lookup(nbinT)) extrap = .true. ! Above upper limit
         ! Lower closest index in table:
         iT = MIN(MAX(int((facT - T_lookup(1))*dlogTinv) + 1, 1), nbinT - 1)
         yy = facT - T_lookup(iT)  ! Dist., in log(T), from T to lower table index
         yy2 = yy*yy             ! That distance squared
         yy3 = yy2*yy            ! ...and cubed
         Tlast = T
      end if

      if (extrap) then ! TK above upper table limit, so extrapolate in log-log:
         alpha = (log10(rates_table%rates(nbinT)) &
                  - log10(rates_table%rates(nbinT - 1))) &
                 /(T_lookup(nbinT) - T_lookup(nbinT - 1))
         inp_coolrates_table = 10d0**(log10(rates_table%rates(nbinT)) &
                                      + alpha*(facT - T_lookup(nbinT)))
         if (present(retPrime)) &
            retPrime = alpha*inp_coolrates_table/T
         return
      end if

      fa = rates_table%rates(iT)           !
      fb = rates_table%rates(iT + 1)         !  Values at neighbouring table
      fprimea = rates_table%primes(iT)       !  indexes
      fprimeb = rates_table%primes(iT + 1)     !

      ! Spline interpolation:
      alpha = fprimea
      beta = (fb - fa)*three_over_h2Table - (2d0*fprimea + fprimeb)*one_over_hTable
      gamma = (fprimea + fprimeb)*one_over_h2Table - (fb - fa)*two_over_h3Table
      inp_coolrates_table = fa + alpha*yy + beta*yy2 + gamma*yy3
      if (present(retPrime)) &
         retPrime = (alpha + 2d0*beta*yy + 3d0*gamma*yy2)/T*one_over_lnTen
   END FUNCTION inp_coolrates_table
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_Alpha_H2(T, Z, xe, H2_cosmic_ray_ionization_rate, G0, xHII_opt, xHI_opt)
! Returns creation rate of h2 on dust [cm^3 s-1] (Draine and Bertoldi 1996)
! plus gas phase rate for low Z on H- assuming equilibrium abundances for H-
! as explained in the Appendix of McKee and Krumholz (2010)
! T           => Temperature [K]
! Z           => Metallicity in Solar units
!-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T, Z, xe, H2_cosmic_ray_ionization_rate, G0
      real(dp), intent(in), optional::xHI_opt, xHII_opt
      real(dp)::comp_Alpha_H2, lambda, T2, clumping
      real(dp)::logT, lnTe, k1, k2, k5, k13, k14, k15, xHI, xHII, k_hm_cr, k_hm_gamma
!-------------------------------------------------------------------------
      T2 = T/1d2
      xHII = 0d0
      xHI = 0d0
      ! initialize
      comp_Alpha_H2 = 0.d0
      if (present(xHII_opt)) xHII = xHII_opt
      if (present(xHI_opt)) xHI = xHI_opt
      if (isH2Katz) then
! Note that harley edited this function.  If you look in the appendix of McKee and
! Krumholz, the zero metallicity limit is also dependent on the electron fraction.
! If you don't include this, the H2 fraction is likely to go crazy at high redshift
! which is probably unphysical.  You need the electrons for the H- channel of H2 formation.
! The real zero metallicity limit should be R- = (8*10^-19)*(ne/10^-3nH)*(T/1000)^0.88

! We change the first part of this ewquation to match equation 4 from
! http://mnras.oxfordjournals.org/content/425/4/3058.full.pdf
! which gets their results from wolfire 2008
         clumping = 1d1 !clumping factor as used in Gnedin 2009
         if (H2clumping > 0) clumping = H2clumping
         comp_Alpha_H2 = Z*(3.5d-17)*clumping*SQRT(T/100.d0) !(1d0 - xHII)
      else
         comp_Alpha_H2 = Z*6.0d-18*(T**0.5)/ &
           & (1.0 + 0.4*T2**0.5 + 0.2*T2 + 0.08*T2**2)
      end if
      if (isH2Katz) then
         ! Glover et al. (2008) Sec 2.1.1
         ! in gas with a high fractional ionization, such as gas recombining from an initially ionized state,
         ! reaction (5) competes with reaction (2) for the available H− ions and
         ! so the uncertainties in the rates of these reactions introduce a significant uncertainty
         ! into the amount of H2 that is formed. A large associative detachment rate and
         ! small mutual neutralization rate lead to the production of a larger H2 fraction
         ! (at a given time) than a small associative detachment rate and
         ! large mutual neutralization rate (Glover et al. 2006).
         logT = log10(T)
         lnTe = log(T*8.621738d-5) ! K -> eV

         ! Creation and destruction channels of H- included with updated rates from Glover et al. 2010
         ! H + e- -> H- + gamma
         if (T < 6000) then
            k1 = 10.0**(-17.845 + 0.762*logT + 0.1523*logT**2.-0.03274*logT**3.)
         else
            k1 = 10.0**(-16.420 + 0.1998*logT**2.-5.447d-3*logT**4.+4.0415d-5*logT**6.)
         end if

         ! H- + H -> H2 + e
         if (T < 300) then
            k2 = 4.0d-9*(300.d0**(-0.17))
         else
            k2 = 4.0d-9*T**(-0.17)
         end if

         ! H- + H+ -> H + H
         k5 = 2.4d-6/sqrt(T)*(1.0 + T/20000.)

         ! H- + CR --> H + e-
         k_hm_cr = 1.28d-13 * (H2_cosmic_ray_ionization_rate / 1.d-16)

         ! H- + gamma --> H + e-
         k_hm_gamma = 5.9d-9 * G0 

         ! H- + e -> H + e + e
         k13 = -1.801849334d1 + 2.36085220d0*lnTe - 2.82744300d-1*lnTe**2. &
            & +1.62331664d-2*lnTe**3.-3.36501203d-2*lnTe**4.+1.17832978d-2*lnTe**5. &
            & -1.65619470d-3*lnTe**6.+1.06827520d-4*lnTe**7.-2.63128581d-6*lnTe**8.
         k13 = exp(k13)

         ! H- + H --> H + H + e-
         ! I think this reaction was broken in glover so I took the results from
         ! https://www.aanda.org/articles/aa/pdf/2016/02/aa27262-15.pdf: Table A1
         if (T .le. 1160.0) then
            k14 = 2.5634d-15 * (T**1.78186)
         else
            k14 = -2.0372609d1 + 1.13944933d0*lnTe - 1.4210135d-1*lnTe**2. &
            & + 8.4644554d-3*lnTe**3. - 1.4328641d-3*lnTe**4. + 2.0122503d-4*lnTe**5. &
            & + 8.6639632d-5*lnTe**6. - 2.5850097d-5*lnTe**7. + 2.4555012d-6*lnTe**8. &
            & - 8.0683825d-8*lnTe**9.
            k14 = exp(k14)
         end if

         ! H- + H+ --> H2+ + e- --> H + H (via recomvinative dissociation)
         if (T .le. 8000.0) then
            k15 = 6.9d-9 * (T**(-0.35))
         else
            k15 = 9.6d-7 * (T**(-0.90))
         end if

         comp_Alpha_H2 = comp_Alpha_H2 + k1*k2*xe/(k2 + k5*xHII + k_hm_cr + k_hm_gamma + k13*xe + k14*xHI + k15*xHII) ! k5=k13=0 should recovers the next line
      else
         ! Zero metallicity limit
         comp_Alpha_H2 = comp_Alpha_H2 + (8.0d-19)*(xe/0.001)*((T/1000.0)**0.88)
      end if

   END FUNCTION comp_Alpha_H2

! H and HE ION RECOMBINATION RATES BADNELL
   ELEMENTAL FUNCTION comp_Alpha_hydrogen_badnell(T)
      ! Returns the recombination rate [cm^3 s-1]
      ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
      ! N is the number of electrons on the ion before recombination
      ! Oxygen atomic number = 8
      implicit none
      real(dp), intent(in)::T
      real(dp)::comp_Alpha_hydrogen_badnell
      real(dp)::A, B, T0, T1, C, T2
      real(dp),dimension(6)::RR_rates
      integer::idx, i

      ! Radiative recombination rates
      RR_rates = (/ 8.318d-11, 0.7472d0, 2.965d0, 7.001d5, 0.0000d0, 0.000d0 /)

      ! Radiative recombination
      A  = RR_rates(1)
      B  = RR_rates(2)
      T0 = RR_rates(3)
      T1 = RR_rates(4)
      C  = RR_rates(5)
      T2 = RR_rates(6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))

      comp_Alpha_hydrogen_badnell = sqrt(T / T0)
      comp_Alpha_hydrogen_badnell = comp_Alpha_hydrogen_badnell * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_hydrogen_badnell = comp_Alpha_hydrogen_badnell * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_hydrogen_badnell = comp_Alpha_hydrogen_badnell**(-1.d0)
      comp_Alpha_hydrogen_badnell = A * comp_Alpha_hydrogen_badnell

   END FUNCTION comp_Alpha_hydrogen_badnell

   ELEMENTAL FUNCTION comp_Alpha_helium_badnell(T, N)
      ! Returns the recombination rate [cm^3 s-1]
      ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
      ! N is the number of electrons on the ion before recombination
      ! Oxygen atomic number = 8
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Alpha_helium_badnell
      real(dp)::A, B, T0, T1, C, T2
      real(dp),dimension(2,6)::RR_rates
      real(dp),dimension(2,9)::DR_rates_c, DR_rates_e
      integer::idx, i
   
      idx = N + 1 ! fortran 1 indexing

      ! Radiative recombination rates
      RR_rates(1,:)  = (/ 1.818d-10, 0.7492d0, 1.017d1, 2.786d6, 0.0000d0, 0.000d0 /)
      RR_rates(2,:)  = (/ 5.235d-11, 0.6988d0, 7.301d0, 4.475d6, 0.0829d0, 1.682d5 /)

      ! Dielectronic recombination rates
      DR_rates_c(1,:)   = (/ 0.000d0,  0.000d0,   0.000d0,  0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0 /)
      DR_rates_c(2,:)   = (/ 5.966d-4, 1.613d-4, -2.223d-5, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0 /)
   
      DR_rates_e(1,:)   = (/ 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0 /)
      DR_rates_e(2,:)   = (/ 4.556d5, 5.552d5, 8.982d5, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0 /)

      ! Radiative recombination
      A  = RR_rates(idx,1)
      B  = RR_rates(idx,2)
      T0 = RR_rates(idx,3)
      T1 = RR_rates(idx,4)
      C  = RR_rates(idx,5)
      T2 = RR_rates(idx,6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))

      comp_Alpha_helium_badnell = sqrt(T / T0)
      comp_Alpha_helium_badnell = comp_Alpha_helium_badnell * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_helium_badnell = comp_Alpha_helium_badnell * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_helium_badnell = comp_Alpha_helium_badnell**(-1.d0)
      comp_Alpha_helium_badnell = A * comp_Alpha_helium_badnell

      ! Dielectronic recombination
      ! T**(-3/2) * sum c_i * exp(-E_i/T)
      do i=1,9
         if (DR_rates_e(idx,i).gt.0.d0) then 
            comp_Alpha_helium_badnell = comp_Alpha_helium_badnell + ( (T**(-3.d0/2.d0)) * DR_rates_c(idx,i) * exp(-1.d0 * DR_rates_e(idx,i) / T) )
         endif
      enddo

   END FUNCTION comp_Alpha_helium_badnell

! METAL ION RECOMBINATION RATES

! OXYGEN
   ELEMENTAL FUNCTION comp_Alpha_oxygen(T, N)
      ! Returns the recombination rate [cm^3 s-1]
      ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
      ! N is the number of electrons on the ion before recombination
      ! Oxygen atomic number = 8
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Alpha_oxygen
      real(dp)::A, B, T0, T1, C, T2
      integer::idx, i
   
      idx = N + 1 ! fortran 1 indexing

      ! Radiative recombination
      A  = RR_rates_oxygen(idx,1)
      B  = RR_rates_oxygen(idx,2)
      T0 = RR_rates_oxygen(idx,3)
      T1 = RR_rates_oxygen(idx,4)
      C  = RR_rates_oxygen(idx,5)
      T2 = RR_rates_oxygen(idx,6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))

      comp_Alpha_oxygen = sqrt(T / T0)
      comp_Alpha_oxygen = comp_Alpha_oxygen * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_oxygen = comp_Alpha_oxygen * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_oxygen = comp_Alpha_oxygen**(-1.d0)
      comp_Alpha_oxygen = A * comp_Alpha_oxygen

      ! Dielectronic recombination
      ! T**(-3/2) * sum c_i * exp(-E_i/T)
      do i=1,9
         if (DR_rates_e_oxygen(idx,i).gt.0.d0) then 
            comp_Alpha_oxygen = comp_Alpha_oxygen + ( (T**(-3.d0/2.d0)) * DR_rates_c_oxygen(idx,i) * exp(-1.d0 * DR_rates_e_oxygen(idx,i) / T) )
         endif
      enddo

   END FUNCTION comp_Alpha_oxygen

! NITROGEN
   ELEMENTAL FUNCTION comp_Alpha_nitrogen(T, N)
      ! Returns the recombination rate [cm^3 s-1]
      ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
      ! N is the number of electrons on the ion before recombination
      ! Nitrogen atomic number = 7
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Alpha_nitrogen
      real(dp)::A, B, T0, T1, C, T2
      integer::idx, i
   
      idx = N + 1 ! fortran 1 indexing

      ! Radiative recombination
      A  = RR_rates_nitrogen(idx,1)
      B  = RR_rates_nitrogen(idx,2)
      T0 = RR_rates_nitrogen(idx,3)
      T1 = RR_rates_nitrogen(idx,4)
      C  = RR_rates_nitrogen(idx,5)
      T2 = RR_rates_nitrogen(idx,6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))

      comp_Alpha_nitrogen = sqrt(T / T0)
      comp_Alpha_nitrogen = comp_Alpha_nitrogen * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_nitrogen = comp_Alpha_nitrogen * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_nitrogen = comp_Alpha_nitrogen**(-1.d0)
      comp_Alpha_nitrogen = A * comp_Alpha_nitrogen

      ! Dielectronic recombination
      ! T**(-3/2) * sum c_i * exp(-E_i/T)
      do i=1,9
         if (DR_rates_e_nitrogen(idx,i).gt.0.d0) then 
            comp_Alpha_nitrogen = comp_Alpha_nitrogen + ( (T**(-3.d0/2.d0)) * DR_rates_c_nitrogen(idx,i) * exp(-1.d0 * DR_rates_e_nitrogen(idx,i) / T) )
         endif
      enddo

   END FUNCTION comp_Alpha_nitrogen

! CARBON
   ELEMENTAL FUNCTION comp_Alpha_carbon(T, N)
      ! Returns the recombination rate [cm^3 s-1]
      ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
      ! N is the number of electrons on the ion before recombination
      ! Carbon atomic number = 6
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Alpha_carbon
      real(dp)::A, B, T0, T1, C, T2
      integer::idx, i
   
      idx = N + 1 ! fortran 1 indexing

      ! Radiative recombination
      A  = RR_rates_carbon(idx,1)
      B  = RR_rates_carbon(idx,2)
      T0 = RR_rates_carbon(idx,3)
      T1 = RR_rates_carbon(idx,4)
      C  = RR_rates_carbon(idx,5)
      T2 = RR_rates_carbon(idx,6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))

      comp_Alpha_carbon = sqrt(T / T0)
      comp_Alpha_carbon = comp_Alpha_carbon * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_carbon = comp_Alpha_carbon * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_carbon = comp_Alpha_carbon**(-1.d0)
      comp_Alpha_carbon = A * comp_Alpha_carbon

      ! Dielectronic recombination
      ! T**(-3/2) * sum c_i * exp(-E_i/T)
      do i=1,9
         if (DR_rates_e_carbon(idx,i).gt.0.d0) then 
            comp_Alpha_carbon = comp_Alpha_carbon + ( (T**(-3.d0/2.d0)) * DR_rates_c_carbon(idx,i) * exp(-1.d0 * DR_rates_e_carbon(idx,i) / T) )
         endif
      enddo

   END FUNCTION comp_Alpha_carbon

! MAGNESIUM
   ELEMENTAL FUNCTION comp_Alpha_magnesium(T, N)
   ! Returns the recombination rate [cm^3 s-1]
   ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
   ! N is the number of electrons on the ion before recombination
   ! Magnesium atomic number = 12
   implicit none
   integer, intent(in)::N
   real(dp), intent(in)::T
   real(dp)::comp_Alpha_magnesium
   real(dp)::A, B, T0, T1, C, T2
   integer::idx, i

   idx = N + 1 ! fortran 1 indexing

   ! Radiative recombination
   A  = RR_rates_magnesium(idx,1)
   B  = RR_rates_magnesium(idx,2)
   T0 = RR_rates_magnesium(idx,3)
   T1 = RR_rates_magnesium(idx,4)
   C  = RR_rates_magnesium(idx,5)
   T2 = RR_rates_magnesium(idx,6)

   !HK NOTE: in many cases C and T2 are 0 so B = B
   B = B + (C * exp(-T2/T))
   comp_Alpha_magnesium = sqrt(T / T0)
   comp_Alpha_magnesium = comp_Alpha_magnesium * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
   comp_Alpha_magnesium = comp_Alpha_magnesium * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
   comp_Alpha_magnesium = comp_Alpha_magnesium**(-1.d0)
   comp_Alpha_magnesium = A * comp_Alpha_magnesium

   ! Dielectronic recombination
   ! T**(-3/2) * sum c_i * exp(-E_i/T)
   do i=1,9
      if (DR_rates_e_magnesium(idx,i).gt.0.d0) then 
         comp_Alpha_magnesium = comp_Alpha_magnesium + ( (T**(-3.d0/2.d0)) * DR_rates_c_magnesium(idx,i) * exp(-1.d0 * DR_rates_e_magnesium(idx,i) / T) )
      endif
   enddo

   END FUNCTION comp_Alpha_magnesium

! Silicon
   ELEMENTAL FUNCTION comp_Alpha_silicon(T, N)
   ! Returns the recombination rate [cm^3 s-1]
   ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
   ! N is the number of electrons on the ion before recombination
   ! Silicon atomic number = 14
   implicit none
   integer, intent(in)::N
   real(dp), intent(in)::T
   real(dp)::comp_Alpha_silicon
   real(dp)::A, B, T0, T1, C, T2
   integer::idx, i

   idx = N + 1 ! fortran 1 indexing

   ! Radiative recombination
   A  = RR_rates_silicon(idx,1)
   B  = RR_rates_silicon(idx,2)
   T0 = RR_rates_silicon(idx,3)
   T1 = RR_rates_silicon(idx,4)
   C  = RR_rates_silicon(idx,5)
   T2 = RR_rates_silicon(idx,6)

   !HK NOTE: in many cases C and T2 are 0 so B = B
   B = B + (C * exp(-T2/T))
   comp_Alpha_silicon = sqrt(T / T0)
   comp_Alpha_silicon = comp_Alpha_silicon * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
   comp_Alpha_silicon = comp_Alpha_silicon * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
   comp_Alpha_silicon = comp_Alpha_silicon**(-1.d0)
   comp_Alpha_silicon = A * comp_Alpha_silicon

   ! Dielectronic recombination
   ! T**(-3/2) * sum c_i * exp(-E_i/T)
   do i=1,9
      if (DR_rates_e_silicon(idx,i).gt.0.d0) then 
         comp_Alpha_silicon = comp_Alpha_silicon + ( (T**(-3.d0/2.d0)) * DR_rates_c_silicon(idx,i) * exp(-1.d0 * DR_rates_e_silicon(idx,i) / T) )
      endif
   enddo

   END FUNCTION comp_Alpha_silicon

! Sulfur
   ELEMENTAL FUNCTION comp_Alpha_sulfur(T, N)
   ! Returns the recombination rate [cm^3 s-1]
   ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
   ! N is the number of electrons on the ion before recombination
   ! Sulfur atomic number = 16
   implicit none
   integer, intent(in)::N
   real(dp), intent(in)::T
   real(dp)::comp_Alpha_sulfur
   real(dp)::A, B, T0, T1, C, T2
   integer::idx, i

   idx = N + 1 ! fortran 1 indexing

   if (idx.le.15) then 

      ! Radiative recombination
      A  = RR_rates_sulfur(idx,1)
      B  = RR_rates_sulfur(idx,2)
      T0 = RR_rates_sulfur(idx,3)
      T1 = RR_rates_sulfur(idx,4)
      C  = RR_rates_sulfur(idx,5)
      T2 = RR_rates_sulfur(idx,6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))
      comp_Alpha_sulfur = sqrt(T / T0)
      comp_Alpha_sulfur = comp_Alpha_sulfur * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_sulfur = comp_Alpha_sulfur * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_sulfur = comp_Alpha_sulfur**(-1.d0)
      comp_Alpha_sulfur = A * comp_Alpha_sulfur

      ! Dielectronic recombination
      ! T**(-3/2) * sum c_i * exp(-E_i/T)
      do i=1,9
         if (DR_rates_e_sulfur(idx,i).gt.0.d0) then 
            comp_Alpha_sulfur = comp_Alpha_sulfur + ( (T**(-3.d0/2.d0)) * DR_rates_c_sulfur(idx,i) * exp(-1.d0 * DR_rates_e_sulfur(idx,i) / T) )
         endif
      enddo

   else

      A = RR_rates_alt_sulfur(idx,1)
      B = RR_rates_alt_sulfur(idx,2)
      comp_Alpha_sulfur = A * ( ( T / 1.0E4 )**( -1.E0 * B )) 

      A  = DR_rates_alt_sulfur(idx,1)
      B  = DR_rates_alt_sulfur(idx,2)
      T0 = DR_rates_alt_sulfur(idx,3)
      T1 = DR_rates_alt_sulfur(idx,4)

      comp_Alpha_sulfur = comp_Alpha_sulfur + ( A * ( T**-1.5E0 ) * exp( -1.E0 * T0 / T ) * ( 1.E0 + B * exp( -1.E0 * T1 / T ) ) )
   endif

   END FUNCTION comp_Alpha_sulfur

! Iron
   ELEMENTAL FUNCTION comp_Alpha_iron(T, N)
   ! Returns the recombination rate [cm^3 s-1]
   ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
   ! N is the number of electrons on the ion before recombination
   ! Iron atomic number = 26
   implicit none
   integer, intent(in)::N
   real(dp), intent(in)::T
   real(dp)::comp_Alpha_iron
   real(dp)::A, B, T0, T1, C, T2
   integer::idx, i

   idx = N + 1 ! fortran 1 indexing

   if (idx.le.15) then 

      ! Radiative recombination
      A  = RR_rates_iron(idx,1)
      B  = RR_rates_iron(idx,2)
      T0 = RR_rates_iron(idx,3)
      T1 = RR_rates_iron(idx,4)
      C  = RR_rates_iron(idx,5)
      T2 = RR_rates_iron(idx,6)

      !HK NOTE: in many cases C and T2 are 0 so B = B
      B = B + (C * exp(-T2/T))
      comp_Alpha_iron = sqrt(T / T0)
      comp_Alpha_iron = comp_Alpha_iron * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
      comp_Alpha_iron = comp_Alpha_iron * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
      comp_Alpha_iron = comp_Alpha_iron**(-1.d0)
      comp_Alpha_iron = A * comp_Alpha_iron

      ! Dielectronic recombination
      ! T**(-3/2) * sum c_i * exp(-E_i/T)
      do i=1,9
         if (DR_rates_e_iron(idx,i).gt.0.d0) then 
            comp_Alpha_iron = comp_Alpha_iron + ( (T**(-3.d0/2.d0)) * DR_rates_c_iron(idx,i) * exp(-1.d0 * DR_rates_e_iron(idx,i) / T) )
         endif
      enddo
   
   else 

      A = RR_rates_alt_iron(idx,1)
      B = RR_rates_alt_iron(idx,2)
      comp_Alpha_iron = A * ( ( T / 1.0E4 )**( -1.E0 * B )) 

      do i=1,4
         comp_Alpha_iron = comp_Alpha_iron + ( (T**-1.5E0) * ( DR_rates_alt_iron(idx,i+4) * exp( -1.E0 * DR_rates_alt_iron(idx,i)/( 8.617333262145E-5 * T ) ) ) )
      enddo
   endif

   END FUNCTION comp_Alpha_iron

! neon
   ELEMENTAL FUNCTION comp_Alpha_neon(T, N)
   ! Returns the recombination rate [cm^3 s-1]
   ! Taken from https://arxiv.org/pdf/astro-ph/0604144.pdf
   ! N is the number of electrons on the ion before recombination
   ! Neon atomic number = 10
   implicit none
   integer, intent(in)::N
   real(dp), intent(in)::T
   real(dp)::comp_Alpha_neon
   real(dp)::A, B, T0, T1, C, T2
   integer::idx, i

   idx = N + 1 ! fortran 1 indexing

   ! Radiative recombination
   A  = RR_rates_neon(idx,1)
   B  = RR_rates_neon(idx,2)
   T0 = RR_rates_neon(idx,3)
   T1 = RR_rates_neon(idx,4)
   C  = RR_rates_neon(idx,5)
   T2 = RR_rates_neon(idx,6)

   !HK NOTE: in many cases C and T2 are 0 so B = B
   B = B + (C * exp(-T2/T))
   comp_Alpha_neon = sqrt(T / T0)
   comp_Alpha_neon = comp_Alpha_neon * ((1.d0 + sqrt(T / T0))**(1.d0 - B))
   comp_Alpha_neon = comp_Alpha_neon * ((1.d0 + sqrt(T / T1))**(1.d0 + B))
   comp_Alpha_neon = comp_Alpha_neon**(-1.d0)
   comp_Alpha_neon = A * comp_Alpha_neon

   ! Dielectronic recombination
   ! T**(-3/2) * sum c_i * exp(-E_i/T)
   do i=1,9
      if (DR_rates_e_neon(idx,i).gt.0.d0) then 
         comp_Alpha_neon = comp_Alpha_neon + ( (T**(-3.d0/2.d0)) * DR_rates_c_neon(idx,i) * exp(-1.d0 * DR_rates_e_neon(idx,i) / T) )
      endif
   enddo

   END FUNCTION comp_Alpha_neon

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_Beta_H2HI(T)

! Returns collisional rate [cm3 s-1] of H2 and HI (Abel 1997->Dove&Mandy 1986)
! T           => Temperature [K]
!-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T
      real(dp)::comp_Beta_H2HI, kbT
!-------------------------------------------------------------------------
      kbT = 8.618d-5*T !eV
      comp_Beta_H2HI = &
         3.324d-9*(kbT**2.012)*exp(-4.463d0/kbT)/(1.d0 + 0.2472d0*kbT)**3.512
   END FUNCTION comp_Beta_H2HI
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_Beta_H2H2(T)

! Returns collisional rate [cm3 s-1] of H2 and H2 (Martin&Keogh&Mandy 1998y)
! T           => Temperature [K]
!-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T
      real(dp)::comp_Beta_H2H2, kbT
!-------------------------------------------------------------------------
      kbT = 8.618*T !eV
      comp_Beta_H2H2 = &
         2.519d-5*(kbt**4.1881)*exp(-0.1731d0/kbT)/(1.d0 + 2.1347d0*kbT)**5.6881
   END FUNCTION comp_Beta_H2H2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION gamma_hi(T, J)
      !-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T, J
      real(dp)::gamma_hi, T3, jconsts
      !-------------------------------------------------------------------------
      T3 = T/1.d3
      jconsts = 0.33 + 0.9*exp(-1.0d0*((J - 3.5d0)/0.9d0)**2)
      gamma_hi = &
         jconsts*(1.0d-11*sqrt(T3)/(1.0d0 + 60.0d0*T3**(-4)) + 1.0d-12*T3)

   END FUNCTION gamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION gamma_h2(T, J)
      !-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T, J
      real(dp)::gamma_h2, T3, jconsts
      !----------------------------------------------------------
      T3 = T/1.d3
      jconsts = (0.276*J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
      gamma_h2 = jconsts*(3.3d-12 + 6.6d-12*T3)

   END FUNCTION gamma_h2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION dgamma_hi(T, J)
      !-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T, J
      real(dp)::dgamma_hi, T3, jconsts
      !-------------------------------------------------------------------------
      T3 = T/1.d3
      jconsts = 0.33 + 0.9*exp(-1.0d0*((J - 3.5d0)/0.9d0)**2)
      dgamma_hi = &
         jconsts*(1.d-11*sqrt(T3)/(1.+60.*T3**(-4))*(1./(2.*T3) + &
                                                     240.*T3**(-5)/(1.+60.*T3**(-4))) + 1.d-12)/1.d3

   END FUNCTION dgamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION dgamma_h2(J)
      !-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::J
      real(dp)::dgamma_h2, jconsts
      !-------------------------------------------------------------------------
      jconsts = 0.276*(J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
      dgamma_h2 = &
         jconsts*6.6d-12

   END FUNCTION dgamma_h2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_Sd(nHI, nH2, dx_SS, Z)
! Returns the self shielding factor for dust
! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
      implicit none
      real(dp), intent(in)::nHI, nH2, dx_SS, Z
      real(dp)::comp_Sd, Sdeff, cNHI, cNH2
      Sdeff = 2.34d-21    !dust cross section cm^2 ! Updated for bare-gr-s
      cNHI = nHI*dx_SS   !HI column density
      cNH2 = nH2*dx_SS   !H2 column density
      comp_Sd = exp(-Sdeff*Z*(cNHI + (2.0*cNH2)))
   END FUNCTION comp_Sd
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_SH2(nH2, dx_SS)
! Returns the self shielding factor for dust
! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
      implicit none
      real(dp), intent(in)::nH2, dx_SS
      real(dp)::comp_SH2, xfac, cNH2, wH2, Sa, Sb, Sc
      cNH2 = nH2*dx_SS  !H2 column density
      xfac = cNH2/(5.0d14)
      wH2 = 0.2

      Sa = (1.0 - wH2)/((1.0 + xfac)*(1.0 + xfac))
      Sb = wH2/sqrt(1.0 + xfac)
      Sc = exp(-0.00085*sqrt(1.0 + xfac))

      comp_SH2 = Sa + (Sb*Sc)
   END FUNCTION comp_SH2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_Beta_H2coll(T, nH, xHI, xH2, xHe, ne, nHI, nH2, nHeI)
! Returns the collisional dissociation rates of H2 for four different
! reactions [cm3s-1] from Glover & Abel (2008)
! http://mnras.oxfordjournals.org/content/388/4/1627.full.pdf
      implicit none
      real(dp), intent(in)::T, nH, xHI, xH2, xHe, ne, nHI, nH2, nHeI
      real(dp)::comp_Beta_H2coll, ncrH, ncrH2, ncrHe, T4, invncr
      real(dp):: LTEfac, NLTEfac
      real(dp):: k8, k8L, k9, k9L, k10, k10L, k11, k11L
      real(dp):: lk8, lk9, lk10, lk11
      real(dp)::loc_T

      loc_T = max(T,3.d0)
      T4 = loc_T/1d4

      ! Critical number densities.  See eqns 15, 16, and 17
      ncrH = 10.d0**(3.d0 - 0.416d0*log10(T4) - 0.327d0*log10(T4)*log10(T4))
      ncrH2 = 10.d0**(4.845d0 - 1.3d0*log10(T4) + 1.62d0*log10(T4)*log10(T4))
      ncrHe = 10.d0**(5.0792d0*(1.d0 - 1.23d-5*(MIN(loc_T,3.d5) - 2000.d0)))

      ! 1/ncr.  see eqn 14
      invncr = 0.d0
      invncr = invncr + (xHI/ncrH)
      invncr = invncr + (xH2/ncrH2)
      invncr = invncr + (xHe/ncrHe)

      ! prefactors for LTE and NLTE collision rates  see eqn 13
      LTEfac = (nH*invncr)/(1.d0 + (nH*invncr))
      NLTEfac = 1.d0/(1.d0 + (nH*invncr))
      if (ne > nHI) then
         LTEfac = 1.d0
         NLTEfac = 0.d0
      end if

      !reaction rates from the appendix:
      !H2 + e- --> H + H + e-
      !k8 = (4.49d-9)*(T**0.11)*exp(-101858.0/T)
      !k8L = (1.91d-9)*(T**0.136)*exp(-53407.1/T)
      k8 = (3.73d-9)*(loc_T**0.1121d0)*exp(-99430.d0/loc_T) ! Glover et al. (2010)

      !H2 + H --> H + H + H
      k9 = (6.67d-12)*sqrt(loc_T)*exp(-1.d0*(1.d0 + (63593.d0/loc_T)))
      k9L = (3.52d-9)*exp(-43900.d0/loc_T)

      !H2 + H2 --> H2 + H + H
      k10 = (((5.996d-30)*(loc_T**4.1881d0))/((1.d0 + (6.761d-6)*loc_T)**5.6881d0))*exp(-54657.4d0/loc_T)
      k10L = (1.3d-9)*exp(-53300.d0/loc_T)

      !H2 + He --> H + H + He
      k11 = 10.d0**(-27.029d0 + (3.801d0*log10(loc_T)) - (29487.d0/loc_T))
      k11L = 10.d0**(-2.729d0 - (1.75d0*log10(loc_T)) - (23474.d0/loc_T))

      k8 = max(k8, 1d-40)
      k9 = max(k9, 1d-40)
      k10 = max(k10, 1d-40)
      k11 = max(k11, 1d-40)
      !k8L  = max(k8L, 1d-40)
      k9L = max(k9L, 1d-40)
      k10L = max(k10L, 1d-40)
      k11L = max(k11L, 1d-40)

      !Log of all the rates
      !lk8  = (LTEfac*log10(k8L)) + (NLTEfac*log10(k8))
      lk8 = log10(k8)
      lk9 = (LTEfac*log10(k9L)) + (NLTEfac*log10(k9))
      lk10 = (LTEfac*log10(k10L)) + (NLTEfac*log10(k10))
      lk11 = (LTEfac*log10(k11L)) + (NLTEfac*log10(k11))

      comp_Beta_H2coll = (ne*(10.0**lk8)) + (nHI*(10.0**lk9)) + (nH2*(10.0**lk10)) + (nHeI*(10.0**lk11))
      comp_Beta_H2coll = MAX(comp_Beta_H2coll/nH, 1.d-40) ! Taysun
   END FUNCTION comp_Beta_H2coll
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION comp_Beta_H2coll_new(T, nH, ne, nHI, nH2, nHeI)
! Returns the collisional dissociation rates of H2 for four different
! reactions [cm3s-1] from Glover et al (2010) (see table B1)
! http://mnras.oxfordjournals.org/content/404/1/2.full.pdf
      implicit none
      real(dp), intent(in)::T, nH, ne, nHI, nH2, nHeI
      real(dp)::comp_Beta_H2coll_new, ncrHe
      real(dp):: k8, k9, k10, k11, k11L, k11H
      real(dp):: bot,loc_T

      loc_T = MIN(MAX(T,2.d3),1.d4)

      !H2 + e- --> H + H + e-
      k8 = (3.73d-9)*(loc_T**0.1121d0)*exp(-99430.d0/loc_T)

      !H2 + H --> H + H + H
      k9 = (6.67d-12)*sqrt(loc_T)*exp(-1.d0*(1.d0 + (63593.d0/loc_T)))

      !H2 + H2 --> H2 + H + H
      k10 = (((5.996d-30)*(loc_T**4.1881d0))/((1.d0 + (6.761d-6)*T)**5.6881d0))*exp(-54657.4d0/loc_T)

      !H2 + He --> H + H + He
      
      k11L = 10.d0**(-27.029d0 + (3.801d0*log10(loc_T)) - (29487.d0/loc_T))
      k11H = 10.d0**(-2.729d0 - (1.75d0*log10(loc_T)) - (23474.d0/loc_T))
      ncrHe = 10.d0**(5.0792d0*(1.d0 - ((1.23d-5)*(loc_T - 2000.d0))))
      bot = (1.d0 + (nHeI/ncrHe)**1.09d0)
      k11 = 10.d0**(log10(k11H) - ((log10(k11H) - log10(k11L))/bot))

      comp_Beta_H2coll_new = (ne*(k8)) + (nHI*(k9)) + (nH2*(k10)) + (nHeI*(k11))

      comp_Beta_H2coll_new = MAX(comp_Beta_H2coll_new/nH, 1.d-40) ! Taysun
   END FUNCTION comp_Beta_H2coll_new
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION Epump(nH, T, xH2, xHI)

! Energy converted to heat from UV pumping
! see Appendix A http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1990ApJ...365..620B&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
      implicit none
      real(dp), intent(in)::nH, T, xH2, xHI
      real(dp)::Epump, Crad, Cdex, Cfrac, ev2erg
      ev2erg = 1.602d-12
      Crad = 2.0d-7     !radiation de-excitation rate Burton 1990
      Cdex = (1.0d-12)*((1.4*xH2*exp(-18100d0/(T + 1200d0))) + (1.0*xHI*exp(-1000d0/T)))*sqrt(T)*nH  !collisional de-excitation rate Burton 1990
      Cfrac = Cdex/(Cdex + Crad)
      Epump = 2*Cfrac*ev2erg  !ergs
   END FUNCTION Epump
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION PE_efficiency(G_0, Tk, n_e)

! Photoelectric heating efficiency from Bakes & Tielens (1994), Baczynski et al. (2015)
! G_0: the strength of the local ISRF normalised to the integrated Habing field
! Tk : temperature of a cell in K
! n_e: electron density in units of cm-3
      implicit none
      real(dp), intent(in)::G_0, Tk, n_e
      real(dp)::phi_pah, fact, PE_efficiency
      phi_pah = 0.5 ! wolfire+(03)
      fact = G_0*sqrt(Tk)/(n_e*phi_pah)
      PE_efficiency = 4.9d-2/(1d0 + 4d-3*fact**0.73) &
        &            + 3.7d-2*(Tk/1d4)**0.7/(1d0 + 2d-4*fact)
   END FUNCTION PE_efficiency
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION H2_cooling_G15(T)
! Cooling due to molecular hydrogen from https://ui.adsabs.harvard.edu/abs/2015MNRAS.453.2901G/abstract
! see appendix 4.2.1 eqn 30
! returns erg/s
      implicit none
      real(dp), intent(in)::T
      real(dp)::H2_cooling_G15
      real(dp)::logT3

      logT3 = LOG10(T/1.d3)
      H2_cooling_G15 = 0.d0
      H2_cooling_G15 = -20.584225d0
      H2_cooling_G15 = H2_cooling_G15 + (5.0194035d0 * logT3)
      H2_cooling_G15 = H2_cooling_G15 - (1.5738805d0 * logT3**2.d0)
      H2_cooling_G15 = H2_cooling_G15 - (4.7155769d0 * logT3**3.d0)
      H2_cooling_G15 = H2_cooling_G15 + (2.4714161d0 * logT3**4.d0)
      H2_cooling_G15 = H2_cooling_G15 + (5.4710750d0 * logT3**5.d0)
      H2_cooling_G15 = H2_cooling_G15 - (3.9467356d0 * logT3**6.d0)
      H2_cooling_G15 = H2_cooling_G15 - (2.2148338d0 * logT3**7.d0)
      H2_cooling_G15 = H2_cooling_G15 + (1.8161874d0 * logT3**8.d0)
      H2_cooling_G15 = (10.d0**(H2_cooling_G15))

   END FUNCTION 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELEMENTAL FUNCTION H2_cooling_GA08(T, n_e, n_HI, n_HII, n_HeI, n_H2)
! Cooling due to molecular hydrogen
      implicit none
      real(dp), intent(in)::T, n_e, n_HI, n_HII, n_HeI, n_H2
      real(dp)::t3, lt, lt3, ltt
      real(dp)::gphdl, HDLR, HDLV, gaHI, gaH2, gaHp, gaHe, gael, galdl
      real(dp)::H2_cooling_GA08

      H2_cooling_GA08 = 0d0
      if (T .le. 10) return

      lt = log10(T)
      t3 = T/1000.
      lt3 = log10(t3)

      ! Galli & Palla (1998) ; low-density limit
      !gpldl = 10.d0**(-103. + 97.59*lt - 48.05*lt**2d0 + 10.80*lt**3d0 - 0.9032*lt**4d0) !erg cm^3 /s

      ! Hollenbach & McKee (1979) ; high-density limit (LTE)
      !HDLR = ((9.5d-22*t3**3.76d0)/(1.d0 + 0.12d0*t3**2.1d0)*exp(-(0.13d0/t3)**3) + 3.d-24*exp(-0.51d0/t3))
      !HDLV = (6.7d-19*exp(-5.86d0/t3) + 1.6d-18*exp(-11.7d0/t3))
      !gphdl = (HDLR + HDLV) ! erg/s
      
      ! From glover 2011
      gphdl = H2_cooling_G15(T)

      ! Initialise
      gaHI = 0.0
      gaH2 = 0.0
      gaHp = 0.0
      gaHe = 0.0
      gael = 0.0

      ! Glover & Abel (2008) ; low-density limit

      ! Excitation by HI
      if (T .lt. 100) then
         gaHI = 10.d0**(-16.818342d0        &
         &         + 37.383713d0*lt3      &
         &         + 58.145166d0*lt3**2   &
         &         + 48.656103d0*lt3**3   &
         &         + 20.159831d0*lt3**4   &
         &         + 3.847961d0*lt3**5)
      else if (T .lt. 1000) then
         gaHI = 10.d0**(-24.311209d0        &
         &         + 3.5692468d0*lt3      &
         &         - 11.33286d0*lt3**2    &
         &         - 27.850082d0*lt3**3   &
         &         - 21.328264d0*lt3**4   &
         &         - 4.2519023d0*lt3**5)
      ! harley checked that extrapolation to 1e4K is ok
      else if (T .lt. 10000) then
         gaHI = 10.d0**(-24.311209d0        &
         &         + 4.6450521d0*lt3       &
         &         - 3.7209846d0*lt3**2    &
         &         + 5.9369081d0*lt3**3    &
         &         - 5.5108047d0*lt3**4    &
         &         + 1.5538288d0*lt3**5)
      end if

      ! Excitation by H2
      ! harley checked that extrapolation to 1e4K is ok
      if (T .lt. 10000) then
         gaH2 = 10.d0**(-23.962112d0         &
         &        + 2.0943374d0*lt3       &
         &        - 0.77151436d0*lt3**2   &
         &        + 0.43693353d0*lt3**3   &
         &        - 0.14913216d0*lt3**4   &
         &        - 0.033638326d0*lt3**5)
      end if

      ! Excitation by He
      ! harley checked that extrapolation to 1e4K is ok
      if (T .lt. 10000.0) then
         gaHe = 10.d0**(-23.689237d0        &
         &        + 2.1892372d0*lt3      &
         &        - 0.81520438d0*lt3**2   &
         &        + 0.29036281d0*lt3**3   &
         &        - 0.16596184d0*lt3**4   &
         &        + 0.19191375d0*lt3**5)
      end if

      ! Excitation by H+
      ! Update from glover 2015
      if (T .lt. 10000) then
         gaHp = 10.d0**(-22.089523d0      &
         &        + 1.5714711d0*lt3       &
         &        + 0.015391166d0*lt3**2  &
         &        - 0.23619985d0*lt3**3   &
         &        - 0.51002221d0*lt3**4   &
         &        + 0.32168730d0*lt3**5)
      end if

      ! Excitation by electrons
      ! Update from glover 2015
      if (T .lt. 500) then
         gael = 10.d0**(-21.928796d0         &
         &           + 16.815730d0*lt3       &
         &           + 96.743155d0*lt3**2.d0    &
         &           + 343.19180d0*lt3**3.d0    &
         &           + 734.71651d0*lt3**4.d0    &
         &           + 983.67576d0*lt3**5.d0    &
         &           + 801.81247d0*lt3**6.d0    &
         &           + 364.14446d0*lt3**7.d0    &
         &           + 70.609154d0*lt3**8.d0)    
      else if (T .lt. 10000) then
         gael = 10.d0**(-22.921189d0         &
         &           + 1.6802758d0*lt3       &
         &           + 0.93310622d0*lt3**2.d0   &
         &           + 4.0406627d0*lt3**3.d0    &
         &           - 4.7274036d0*lt3**4.d0    &
         &           - 8.8077017d0*lt3**5.d0    &
         &           + 8.9167183d0*lt3**6.d0    &
         &           + 6.4380698d0*lt3**7.d0    &
         &           - 6.3701156d0*lt3**8.d0)    
      end if

      galdl = gaHI*n_HI + gaH2*n_H2 + gaHe*n_HeI + gaHp*n_HII + gael*n_e ! erg/s

      H2_cooling_GA08 = n_H2*gphdl/(1.d0 + gphdl/galdl) ! erg/cm^3/s

   END FUNCTION H2_cooling_GA08
    !*******************************
   FUNCTION cooling_H2GP(nH,nH2,Tgas)
      !cooling from Galli&Palla98
      !taken fron krome
      real(dp)::nH,nH2,Tgas, tm, logT
      real(dp)::cooling_H2GP,T3
      real(dp)::LDL,HDLR,HDLV,HDL

      tm = max(Tgas, 13.0d0)    ! no cooling below 13 Kelvin
      tm = min(Tgas, 1.d5)      ! fixes numerics
      logT = log10(tm)
      T3 = tm * 1.d-3

      !low density limit in erg/s
      LDL = 1.d1**(-103.d0+97.59d0*logT-48.05d0*logT**2 + 10.8d0*logT**3-0.9032d0*logT**4)*nH

      !this will avoid a division by zero and useless calculations
      if(LDL==0d0) then
         cooling_H2GP = 0d0
         return
      end if

      !high density limit
      HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
           3.e-24*exp(-0.51/t3)) !erg/s
      HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
      HDL  = HDLR + HDLV !erg/s

      !to avoid division by zero
      if (HDL==0d0) then
         cooling_H2GP = 0d0
      else
         cooling_H2GP = nH2/(1d0/HDL+1d0/LDL) !erg/cm3/s
      endif

    end FUNCTION cooling_H2GP
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   FUNCTION cooling_H2_moseley(nH,nH2,Tgas)
      implicit none
      real(dp)::nH,nH2,Tgas
      real(dp)::n1, n2, n3, n4
      real(dp)::x1, x2, x3, x4
      real(dp)::f1, f2, f3, f4
      real(dp)::T3
      real(dp)::cooling_H2_moseley

      T3 = 1.d-3 * Tgas

      n1 = 50.d0
      n2 = 450.d0
      n3 = 25.d0
      n4 = 900.d0

      x1 = nH + (5.d0 * nH2)
      x2 = nH + (4.5d0 * nH2)
      x3 = nH + (0.75d0 * nH2)
      x4 = nH + (0.05 * nH2)

      f1 = 1.1d-25 * SQRT(T3) * EXP(-0.51d0 / T3)
      f1 = f1 * (((0.7d0 * x1) / (1.d0 + (x1/n1))) + ((0.3d0 * x1)/(1.d0 + (x1/(10.d0*n1)))))

      f2 = 2.0d-25 * T3 * EXP(-1.d0 / T3)
      f2 = f2 * (((0.35d0 * x2) / (1.d0 + (x2/n2))) + ((0.65d0 * x2)/(1.d0 + (x2/(10.d0*n2)))))

      f3 = 2.4d-24 * (T3**(1.5d0)) * EXP(-2.d0 / T3)
      f3 = f3 * (x3 / (1.d0 + (x3/n3)))

      f4 = 1.7d-23 * (T3**(1.5d0)) * EXP(-4.d0 / T3)
      f4 = f4 * (((0.45d0 * x4) / (1.d0 + (x4/n4))) + ((0.55d0 * x4)/(1.d0 + (x4/(10.d0*n4)))))


      cooling_H2_moseley = nH2 * (f1 + f2 + f3 + f4)

   end FUNCTION cooling_H2_moseley
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   FUNCTION compCoolrate(T, ne, nH, nN, nI, aexp, dcooldT)

! Compute cooling rate in a cell, using interpolation from cooling rate
! tables
! T        => Cell emperature [K]
! nN       => Neutral abundances
! nI       => Ionized abundances
! nH       => Hydrogen number density [cm-3]
! aexp     => Cosmic expansion
! dcooldT <=  Temperature derivative of the rate
! dcooldx <=  Ionized fraction derivative of the rate
! returns:  Resulting cooling rate [erg s-1 cm-3]
!-------------------------------------------------------------------------
      implicit none
      real(dp)::T, ne, nH, aexp
      real(dp)::compCoolrate, dcooldT
      real(dp), dimension(nIons)::nN, nI
      logical::RT_OTSA!-------------------------------------------------------
      real(dp), save::nHI, nHII, nH2, nHeI, nHeII, nHeIII
      real(dp), save::ci_HI, ci_HeI, ci_HeII
      real(dp), save::ci_HI_prime, ci_HeI_prime, ci_HeII_prime
      real(dp), save::ce_HI, ce_HeI, ce_HeII
      real(dp), save::ce_HI_prime, ce_HeI_prime, ce_HeII_prime
      real(dp), save::r_HII, r_HeII, r_HeIII
      real(dp), save::r_HII_prime, r_HeII_prime, r_HeIII_prime
      real(dp), save::bre, bre_prime, com, com_prime, die, die_prime
      real(dp), save::ne_nHI, ne_nHII, ne_nHeI, ne_nHeII, ne_nHeIII
      real(dp), save::lowrleft, lowrright, lowr_hi, lowr_h2
      real(dp), save::lowvleft_hi, lowvright_hi, lowv_hi, lowv_h2
      real(dp), save::lowtot_hi, lowtot_h2, coolH2, dcoolH2, TT
!-------------------------------------------------------------------------
      nHI = nN(ixHII)
      nHII = nI(ixHII)
      if (isH2) nH2 = nN(ixHI)
      if (isHe) then
         nHeI = nN(ixHeII)
         nHeII = nI(ixHeII)
         nHeIII = nI(ixHeIII)
      end if

      ne_nHI = ne*nHI
      ne_nHII = ne*nHII
      if (isHe) then
         ne_nHeI = ne*nHeI
         ne_nHeII = ne*nHeII
         ne_nHeIII = ne*nHeIII
      end if

      ! Coll. Ionization Cooling
      ci_HI = inp_coolrates_table(tbl_cr_ci_HI, T, ci_HI_prime) &
              *ne_nHI
      ci_HI_prime = ci_HI_prime*ne_nHI
      if (isHe) then
         ci_HeI = inp_coolrates_table(tbl_cr_ci_HeI, T, ci_HeI_prime) &
                  *ne_nHeI
         ci_HeII = inp_coolrates_table(tbl_cr_ci_HeII, T, ci_HeII_prime) &
                   *ne_nHeII
         ci_HeI_prime = ci_HeI_prime*ne_nHeI
         ci_HeII_prime = ci_HeII_prime*ne_nHeII
      end if

      ! Collisional excitation cooling
      ce_HI = inp_coolrates_table(tbl_cr_ce_HI, T, ce_HI_prime) &
              *ne_nHI
      ce_HI_prime = ce_HI_prime*ne_nHI
      if (isHe) then
         ce_HeI = inp_coolrates_table(tbl_cr_ce_HeI, T, ce_HeI_prime) &
                  *ne_nHeI
         ce_HeII = inp_coolrates_table(tbl_cr_ce_HeII, T, ce_HeII_prime) &
                   *ne_nHeII
         ce_HeI_prime = ce_HeI_prime*ne_nHeI
         ce_HeII_prime = ce_HeII_prime*ne_nHeII
      end if

      ! Recombination Cooling
      r_HII = inp_coolrates_table(tbl_cr_r_HII, T, r_HII_prime) &
              *ne_nHII
      r_HII_prime = r_HII_prime*ne_nHII
      if (isHe) then
         r_HeII = inp_coolrates_table(tbl_cr_r_HeII, T, r_HeII_prime) &
                  *ne_nHeII
         r_HeIII = inp_coolrates_table(tbl_cr_r_HeIII, T, r_HeIII_prime) &
                   *ne_nHeIII
         r_HeII_prime = r_HeII_prime*ne_nHeII
         r_HeIII_prime = r_HeIII_prime*ne_nHeIII
      end if

      ! Bremsstrahlung
      bre = inp_coolrates_table(tbl_cr_bre, T, bre_prime) &
            *ne*(nHII + nHeII + 4.*nHeIII)
      bre_prime = bre_prime*ne*(nHII + nHeII + 4.*nHeIII)

      ! Compton Cooling
      com = inp_coolrates_table(tbl_cr_com, T, com_prime)*ne
      com_prime = com_prime*ne

      ! Dielectronic recombination cooling
      if (isHe) then
         die = inp_coolrates_table(tbl_cr_die, T, die_prime)*ne_nHeII
         die_prime = die_prime*ne_nHeII
      end if

      ! H2 cooling Hallenbacn McKee (1979) + Halle Combes (2012) in cgs
      if (isH2) then
         ! Collisional cooling
         !lowrleft=1.25*exp(-1.70d2/T)*2.35d-14
         !lowrright=1.75*exp(-5.05d2/T)*6.97d-14
         !lowr_hi=lowrleft*gamma_hi(T,2.d0)+lowrright*gamma_hi(T,3.d0)
         !lowr_h2=lowrleft*gamma_h2(T,2.d0)+lowrright*gamma_h2(T,3.d0)
         !lowvleft_hi=exp(-5860.0/T)*8.09d-13*1.0d-12*sqrt(T)*exp(-1.0d3/T)
         !lowvright_hi=exp(-11720.0/T)*1.6d-12*1.6d-12*sqrt(T)*exp(-1.0*(4.0d2/T)**2)
         !lowv_hi=lowvleft_hi+lowvright_hi
         !lowv_h2=exp(-5860.0/T)*8.09d-13*1.4d-12*sqrt(T)*exp(-1.2d4/(T+1.2d3))
         !lowtot_hi=lowr_hi+lowv_hi
         !lowtot_h2=lowr_h2+lowv_h2

         !coolH2 = nH2*(nHI*lowtot_hi+nH2*lowtot_h2)

         ! Same thing for the temperature derivative
         !TT=T*1.0001

         ! Collisional cooling
         !lowrleft=1.25*exp(-1.70d2/TT)*2.35d-14
         !lowrright=1.75*exp(-5.05d2/TT)*6.97d-14
         !lowr_hi=lowrleft*gamma_hi(TT,2.d0)+lowrright*gamma_hi(TT,3.d0)
         !lowr_h2=lowrleft*gamma_h2(TT,2.d0)+lowrright*gamma_h2(TT,3.d0)
         !lowvleft_hi=exp(-5860.0/TT)*8.09d-13*1.0d-12*sqrt(TT)*exp(-1.0d3/TT)
         !lowvright_hi=exp(-11720.0/TT)*1.6d-12*1.6d-12*sqrt(TT)*exp(-1.0*(4.0d2/TT)**2)
         !lowv_hi=lowvleft_hi+lowvright_hi
         !lowv_h2=exp(-5860.0/TT)*8.09d-13*1.4d-12*sqrt(TT)*exp(-1.2d4/(TT+1.2d3))
         !lowtot_hi=lowr_hi+lowv_hi
         !lowtot_h2=lowr_h2+lowv_h2

         !dcoolH2 = nH2*(nHI*lowtot_hi+nH2*lowtot_h2)

         coolH2 = 0.d0
         dcoolH2 = 0.d0

         ! We are using galli and palla because it is more stable
         if (T.lt.1.d4) then
            !coolH2 = H2_cooling_GA08(T, ne, nHI, nHII, nHeI, nH2)
            !coolH2 = cooling_H2GP(nHI+nHII,nH2,T)
            coolH2 = cooling_H2_moseley(nHI+nHII,nH2,T)

            TT = T*1.001

            !Deal with doscontinuities by pushing derivative the other way
            !if (T.lt.1000.d0.and.TT.ge.1000.d0) TT = 0.999 * T
            !if (T.lt.500.d0.and.TT.ge.500.d0) TT = 0.999 * T
            !if (T.lt.100.d0.and.TT.ge.100.d0) TT = 0.999 * T

            !dcoolH2 = H2_cooling_GA08(TT, ne, nHI, nHII, nHeI, nH2)
            !dcoolH2 = cooling_H2GP(nHI+nHII,nH2,TT)
            dcoolH2 = cooling_H2_moseley(nHI+nHII,nH2,TT)

            dcoolH2 = (dcoolH2 - coolH2)/(TT-T)

         end if

      end if

      ! Overall Cooling
      compCoolrate = ci_HI + r_HII + ce_HI + com + bre
      if (isHe) compCoolrate = compCoolrate &
                               + ci_HeI + r_HeII + ce_HeI + die &
                               + ci_HeII + r_HeIII + ce_HeII
      if (isH2) compCoolrate = compCoolrate + coolH2

      dCooldT = ci_HI_prime + r_HII_prime + ce_HI_prime &
                + bre_prime + com_prime
      if (isHe) dCooldT = dCooldT &
                          + ci_HeI_prime + r_HeII_prime + ce_HeI_prime &
                          + ci_HeII_prime + r_HeIII_prime + ce_HeII_prime &
                          + die_prime
      if (isH2) dCooldT = dCooldT + dcoolH2

   END FUNCTION compCoolrate

   ELEMENTAL FUNCTION comp_Beta_hydrogen_cloudy(T)
! Returns collisional rate [cm3 s-1] of O and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state + 1
!-------------------------------------------------------------------------
      implicit none
      real(dp), intent(in)::T
      real(dp)::comp_Beta_hydrogen_cloudy
      real(dp)::dE, A, X, K
      integer::P
      real(dp)::U

      dE = 13.6d0
      P = 0
      A = 2.91d-8
      X = 0.232d0
      K = 0.39d0

      U = dE / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_hydrogen_cloudy = U**K
      comp_Beta_hydrogen_cloudy = comp_Beta_hydrogen_cloudy * exp(-1.d0 * U)
      comp_Beta_hydrogen_cloudy = comp_Beta_hydrogen_cloudy * (1.d0 + (float(P) * sqrt(U)))
      comp_Beta_hydrogen_cloudy = comp_Beta_hydrogen_cloudy / (X + U)
      comp_Beta_hydrogen_cloudy = comp_Beta_hydrogen_cloudy * A
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_hydrogen_cloudy

   ELEMENTAL FUNCTION comp_Beta_helium_cloudy(T,N)
! Returns collisional rate [cm3 s-1] of O and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state + 1
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_helium_cloudy
      real(dp),dimension(1:2)::dE, A, X, K
      integer,dimension(1:2)::P
      real(dp)::U

      dE = (/ 24.6d0, 54.4d0 /)
      P = (/ 0, 1 /)
      A = (/ 1.75d-8, 2.05d-9 /)
      X = (/ 0.18d0, 0.265d0 /)
      K = (/ 0.35d0, 0.25d0/)

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_helium_cloudy = U**K(N)
      comp_Beta_helium_cloudy = comp_Beta_helium_cloudy * exp(-1.d0 * U)
      comp_Beta_helium_cloudy = comp_Beta_helium_cloudy * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_helium_cloudy = comp_Beta_helium_cloudy / (X(N) + U)
      comp_Beta_helium_cloudy = comp_Beta_helium_cloudy * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_helium_cloudy

   FUNCTION comp_Beta_oxygen(T,N)
! Returns collisional rate [cm3 s-1] of O and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state + 1
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_oxygen
      real(dp),dimension(1:8),save::dE=(/ 13.6d0, 35.1d0, 54.9d0, 77.4d0, 113.9d0, 138.1d0, 739.3d0, 871.5d0 /)
      real(dp),dimension(1:8),save::A=(/ 0.359d-7, 0.139d-7, 0.931d-8, 0.102d-7, 0.219d-8, 0.195d-8, 0.212d-9, 0.521d-10 /)
      real(dp),dimension(1:8),save::X=(/ 0.073d0, 0.212d0, 0.270d0, 0.614d0, 0.630d0, 0.360d0, 0.396d0, 0.629d0 /)
      real(dp),dimension(1:8),save::K=(/ 0.34d0, 0.22d0, 0.27d0, 0.27d0, 0.17d0, 0.54d0, 0.35d0, 0.16d0/)
      integer,dimension(1:8),save::P=(/ 0, 1, 1, 0, 1, 0, 0, 1 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_oxygen = U**K(N)
      comp_Beta_oxygen = comp_Beta_oxygen * exp(-1.d0 * U)
      comp_Beta_oxygen = comp_Beta_oxygen * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_oxygen = comp_Beta_oxygen / (X(N) + U)
      comp_Beta_oxygen = comp_Beta_oxygen * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_oxygen

   FUNCTION comp_Beta_nitrogen(T,N)
! Returns collisional rate [cm3 s-1] of N and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_nitrogen
      real(dp),dimension(1:7),save::dE=(/ 14.5d0, 29.6d0, 47.5d0, 77.5d0, 97.9d0, 552.1d0, 667.0d0 /)
      real(dp),dimension(1:7),save::A=(/ 0.482d-7, 0.298d-7, 0.810d-8, 0.371d-8, 0.151d-8, 0.371d-9, 0.777d-10 /)
      real(dp),dimension(1:7),save::X=(/ 0.0652d0, 0.310d0, 0.350d0, 0.549d0, 0.0167d0, 0.546d0, 0.624d0 /)
      real(dp),dimension(1:7),save::K=(/ 0.42d0, 0.30d0, 0.24d0, 0.18d0, 0.74d0, 0.29d0, 0.16d0/)
      integer,dimension(1:7),save::P=(/ 0, 0, 1, 1, 0, 0, 1/)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_nitrogen = U**K(N)
      comp_Beta_nitrogen = comp_Beta_nitrogen * exp(-1.d0 * U)
      comp_Beta_nitrogen = comp_Beta_nitrogen * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_nitrogen = comp_Beta_nitrogen / (X(N) + U)
      comp_Beta_nitrogen = comp_Beta_nitrogen * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_nitrogen

   FUNCTION comp_Beta_carbon(T,N)
! Returns collisional rate [cm3 s-1] of C and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_carbon
      real(dp),dimension(1:6),save::dE=(/ 11.3d0, 24.4d0, 47.9d0, 64.5d0, 392.1d0, 490.0d0 /)
      real(dp),dimension(1:6),save::A=(/ 0.685d-7, 0.186d-7, 0.635d-8, 0.150d-8, 0.299d-9, 0.123d-9 /)
      real(dp),dimension(1:6),save::X=(/ 0.193d0, 0.286d0, 0.427d0, 0.416d0, 0.666d0, 0.620d0 /)
      real(dp),dimension(1:6),save::K=(/ 0.25d0, 0.24d0, 0.21d0, 0.13d0, 0.02d0, 0.16d0 /)
      integer,dimension(1:6),save::P=(/ 0, 1, 1, 1, 1, 1 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_carbon = U**K(N)
      comp_Beta_carbon = comp_Beta_carbon * exp(-1.d0 * U)
      comp_Beta_carbon = comp_Beta_carbon * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_carbon = comp_Beta_carbon / (X(N) + U)
      comp_Beta_carbon = comp_Beta_carbon * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_carbon

   FUNCTION comp_Beta_magnesium(T,N)
! Returns collisional rate [cm3 s-1] of Mg and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_magnesium
      real(dp),dimension(1:12),save::dE=(/ 7.6d0, 15.2d0, 80.1d0, 109.3d0, 141.3d0, 186.5d0, 224.9d0, &
              266.0d0, 328.2d0, 367.5d0, 1761.8d0, 1962.7d0/)
      real(dp),dimension(1:12),save::A=(/ 0.621d-6, 0.192d-7, 0.556d-8, 0.435d-8, 0.710d-8, 0.170d-8, &
             0.122d-8, 0.220d-8, 0.486d-9, 0.235d-9, 0.206d-10, 0.175d-10  /)
      real(dp),dimension(1:12),save::X=(/ 0.592d0, 0.0027d0, 0.107d0, 0.159d0, 0.658d0, 0.242d0, 0.343d0, &
             0.897d0, 0.751d0, 1.030d0, 0.196d0, 0.835d0 /)
      real(dp),dimension(1:12),save::K=(/ 0.39d0, 0.85d0, 0.30d0, 0.31d0, 0.25d0, 0.28d0, 0.23d0, 0.22d0, &
             0.14d0, 0.10d0, 0.25d0, 0.11d0 /)
      integer,dimension(1:12),save::P=(/ 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_magnesium = U**K(N)
      comp_Beta_magnesium = comp_Beta_magnesium * exp(-1.d0 * U)
      comp_Beta_magnesium = comp_Beta_magnesium * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_magnesium = comp_Beta_magnesium / (X(N) + U)
      comp_Beta_magnesium = comp_Beta_magnesium * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_magnesium

   FUNCTION comp_Beta_silicon(T,N)
! Returns collisional rate [cm3 s-1] of Si and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_silicon
      real(dp),dimension(1:14),save::dE=(/ 8.2d0, 16.4d0, 33.5d0, 54.0d0, 166.8d0, 205.3d0, 246.5d0, & 
              303.5d0, 351.1d0, 401.4d0, 476.4d0, 523.5d0, 2437.7d0, 2673.2d0/)
      real(dp),dimension(1:14),save::A=(/ 0.188d-6, 0.643d-7, 0.207d-7, 0.494d-8, 0.176d-8, 0.174d-8, 0.123d-8, & 
             0.827d-9, 0.601d-8, 0.465d-9, 0.263d-9, 0.118d-9, 0.336d-10, 0.119d-10 /)
      real(dp),dimension(1:14),save::X=(/ 0.376d0, 0.632d0, 0.473d0, 0.172d0, 0.102d0, 0.180d0, 0.518d0, & 
             0.239d0, 0.305d0, 0.666d0, 0.666d0, 0.734d0, 0.336d0, 0.989d0 /)
      real(dp),dimension(1:14),save::K=(/ 0.25d0, 0.20d0, 0.22d0, 0.23d0, 0.31d0, 0.29d0, 0.07d0, & 
             0.28d0, 0.25d0, 0.04d0, 0.16d0, 0.16d0, 0.37d0, 0.08d0 /)
      integer,dimension(1:14),save::P=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_silicon = U**K(N)
      comp_Beta_silicon = comp_Beta_silicon * exp(-1.d0 * U)
      comp_Beta_silicon = comp_Beta_silicon * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_silicon = comp_Beta_silicon / (X(N) + U)
      comp_Beta_silicon = comp_Beta_silicon * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_silicon

   FUNCTION comp_Beta_sulfur(T,N)
! Returns collisional rate [cm3 s-1] of S and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_sulfur
      real(dp),dimension(1:16),save::dE=(/ 10.4d0, 23.3d0, 34.8d0, 47.3d0, 72.6d0, 88.1d0, 280.9d0, 328.2d0, & 
              379.1d0, 447.1d0, 504.8d0, 564.7d0, 651.6d0, 707.2d0, 3223.9d0, 3494.2d0 /)
      real(dp),dimension(1:16),save::A=(/ 0.549d-7, 0.681d-7, 0.214d-7, 0.166d-7, 0.612d-8, 0.133d-8, 0.493d-8, 0.873d-9, & 
             0.135d-8, 0.459d-9, 0.349d-9, 0.523d-9, 0.259d-9, 0.750d-10, 0.267d-10, 0.632d-11 /)
      real(dp),dimension(1:16),save::X=(/ 0.100d0, 0.693d0, 0.353d0, 1.030d0, 0.580d0, 0.0688d0, 1.130d0, 0.193d0, & 
             0.431d0, 0.242d0, 0.305d0, 0.428d0, 0.854d0, 0.734d0, 0.572d0, 0.585d0 /)
      real(dp),dimension(1:16),save::K=(/ 0.25d0, 0.21d0, 0.24d0, 0.14d0, 0.19d0, 0.35d0, 0.16d0, 0.28d0, 0.32d0, 0.28d0, & 
             0.25d0, 0.35d0, 0.12d0, 0.16d0, 0.28d0, 0.17d0 /)
      integer,dimension(1:16),save::P=(/ 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_sulfur = U**K(N)
      comp_Beta_sulfur = comp_Beta_sulfur * exp(-1.d0 * U)
      comp_Beta_sulfur = comp_Beta_sulfur * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_sulfur = comp_Beta_sulfur / (X(N) + U)
      comp_Beta_sulfur = comp_Beta_sulfur * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_sulfur

   FUNCTION comp_Beta_iron(T,N)
! Returns collisional rate [cm3 s-1] of Fe and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_iron
      real(dp),dimension(1:26),save::dE=(/ 7.9d0, 16.2d0, 30.6d0, 54.8d0, 75.0d0, 99.0d0, 125.0d0, 151.1d0, & 
              233.6d0, 262.1d0, 290.0d0, 331.0d0, 361.0d0, 392.0d0, 457.0d0, & 
              489.3d0, 1262.0d0, 1360.0d0, 1470.0d0, 1582.0d0, 1690.0d0, 1800.0d0, &
              1960.0d0, 2046.0d0, 8828.0d0, 9277.7d0 /)
      real(dp),dimension(1:26),save::A=(/ 0.252d-6, 0.221d-7, 0.410d-7, 0.353d-7, 0.104d-7, 0.123d-7, 0.947d-8, & 
             0.471d-8, 0.302d-8, 0.234d-8, 0.176d-8, 0.114d-8, 0.866d-9, 0.661d-9, &
             0.441d-9, 0.118d-9, 0.361d-9, 0.245d-9, 0.187d-9, 0.133d-9, 0.784d-10, &
             0.890d-10, 0.229d-10, 0.112d-10, 0.246d-11, 0.979d-12 /)
      real(dp),dimension(1:26),save::X=(/ 0.701d0, 0.033d0, 0.366d0, 0.243d0, 0.285d0, 0.411d0, 0.458d0, 0.280d0, & 
             0.697d0, 0.764d0, 0.805d0, 0.773d0, 0.805d0, 0.762d0, 0.698d0, 0.211d0, &
             1.160d0, 0.978d0, 0.988d0, 1.030d0, 0.848d0, 1.200d0, 0.936d0, 0.034d0, &
             1.020d0, 0.664d0 /)
      real(dp),dimension(1:26),save::K=(/ 0.25d0, 0.45d0, 0.17d0, 0.39d0, 0.17d0, 0.21d0, 0.21d0, 0.28d0, 0.15d0, & 
             0.14d0, 0.14d0, 0.15d0, 0.14d0, 0.14d0, 0.16d0, 0.15d0, 0.09d0, 0.13d0, & 
             0.14d0, 0.12d0, 0.14d0, 0.35d0, 0.12d0, 0.81d0, 0.02d0, 0.14d0 /)
      integer,dimension(1:26),save::P=(/ 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_iron = U**K(N)
      comp_Beta_iron = comp_Beta_iron * exp(-1.d0 * U)
      comp_Beta_iron = comp_Beta_iron * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_iron = comp_Beta_iron / (X(N) + U)
      comp_Beta_iron = comp_Beta_iron * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_iron

   FUNCTION comp_Beta_neon(T,N)
! Returns collisional rate [cm3 s-1] of Ne and e 
! https://reader.elsevier.com/reader/sd/pii/S0092640X97907324?token=A301709E037A1567FDE65EB733BB04A7FDE698521472A4F5B45E6507B7DD95428A8DBF3D23DB9A941CEF78C5D2A507C6&originRegion=eu-west-1&originCreation=20211004171613
! T           => Temperature [K]
! N           => The ionization state 
!-------------------------------------------------------------------------
      implicit none
      integer, intent(in)::N
      real(dp), intent(in)::T
      real(dp)::comp_Beta_neon
      real(dp),dimension(1:10),save::dE=(/ 21.6d0, 41.0d0, 63.5d0, 97.1d0, 126.2d0, &
              157.9d0, 207.3d0, 239.1d0, 1196.0d0, 1360.6d0 /)
      real(dp),dimension(1:10),save::A=(/  0.150d-7, 0.198d-7, 0.703d-8, 0.424d-8, 0.279d-8, &
              0.345d-8, 0.956d-9, 0.473d-9, 0.392d-10, 0.277d-10/)
      real(dp),dimension(1:10),save::X=(/ 0.0329d0, 0.295d0, 0.0677d0, 0.0482d0, 0.305d0, &
             0.581d0, 0.749d0, 0.992d0, 0.262d0, 0.661d0 /)
      real(dp),dimension(1:10),save::K=(/ 0.43d0, 0.20d0, 0.39d0, 0.58d0, 0.25d0, &
             0.28d0, 0.14d0, 0.04d0, 0.20d0, 0.13d0 /)
      integer,dimension(1:10),save::P=(/ 1, 0, 1, 1, 1, 0, 1, 1, 1, 1 /)
      real(dp)::U

      U = dE(N) / (T * 8.61732814974056d-05) ! K --> eV
      comp_Beta_neon = U**K(N)
      comp_Beta_neon = comp_Beta_neon * exp(-1.d0 * U)
      comp_Beta_neon = comp_Beta_neon * (1.d0 + (float(P(N)) * sqrt(U)))
      comp_Beta_neon = comp_Beta_neon / (X(N) + U)
      comp_Beta_neon = comp_Beta_neon * A(N)
!-------------------------------------------------------------------------
   END FUNCTION comp_Beta_neon

   FUNCTION DUST_RECOMBINATION(elem,G,T,ne)
      ! Recombination of ions on dust
      ! See https://iopscience.iop.org/article/10.1086/324035/pdf
      ! Note that this approximation is only good between 10 and 10^3 K
      implicit none
      real(dp), intent(in)::G,T,ne
      integer, intent(in)::elem
      real(dp),dimension(1:7)::Hp,Hep,Cp,Mgp,Sip,Sp,Fep
      real(dp)::C0,C1,C2,C3,C4,C5,C6,phi
      real(dp)::a1,a2,a3
      real(dp)::DUST_RECOMBINATION
      real(dp)::dr_sf

      DUST_RECOMBINATION = 0.d0
      dr_sf              = 1.d0

      ! Return at high temperatures
      if (T.gt.1.d4) then
         return
      ! Scale down if gtr than 1.d3
      else if (T.gt.1.d3) then
         dr_sf = EXP(-1.d0 * T / 1.d3) / EXP(-1.d0)
      end if

      ! No recombination below 10 K --> rates inaccurate
      if (T.lt.10.d0) return 

      Hp  = (/ 12.25d0, 8.074d-6, 1.378d0, 5.087d2, 1.586d-2, 0.4723d0, 1.102d-5 /)
      Hep = (/ 5.572d0, 3.185d-7, 1.512d0, 5.115d3, 3.903d-7, 0.4956d0, 5.494d-7 /)
      Cp  = (/ 45.58d0, 6.089d-3, 1.128d0, 4.331d2, 4.845d-2, 0.8120d0, 1.333d-4 /)
      Mgp = (/ 2.510d0, 8.116d-8, 1.864d0, 6.170d4, 2.169d-6, 0.9605d0, 7.232d-5 /)
      Sip = (/ 2.166d0, 5.678d-8, 1.874d0, 4.375d4, 1.635d-6, 0.8964d0, 7.538d-5 /)
      Sp  = (/ 3.064d0, 7.769d-5, 1.319d0, 1.087d2, 3.475d-1, 0.4790d0, 4.689d-2 /)
      Fep = (/ 1.701d0, 9.554d-8, 1.851d0, 5.763d4, 4.116d-8, 0.9456d0, 2.198d-5 /)

      if (elem.eq.1) then ! Hydrogen
         C0 = Hp(1)
         C1 = Hp(2)
         C2 = Hp(3)
         C3 = Hp(4)
         C4 = Hp(5)
         C5 = Hp(6)
         C6 = Hp(7)
      else if (elem.eq.2) then ! Helium
         C0 = Hep(1)
         C1 = Hep(2)
         C2 = Hep(3)
         C3 = Hep(4)
         C4 = Hep(5)
         C5 = Hep(6)
         C6 = Hep(7)
      else if (elem.eq.6) then ! Carbon
         C0 = Cp(1)
         C1 = Cp(2)
         C2 = Cp(3)
         C3 = Cp(4)
         C4 = Cp(5)
         C5 = Cp(6)
         C6 = Cp(7)
      else if (elem.eq.12) then ! Magnedium
         C0 = Mgp(1)
         C1 = Mgp(2)
         C2 = Mgp(3)
         C3 = Mgp(4)
         C4 = Mgp(5)
         C5 = Mgp(6)
         C6 = Mgp(7)
      else if (elem.eq.14) then ! Silicon
         C0 = Sip(1)
         C1 = Sip(2)
         C2 = Sip(3)
         C3 = Sip(4)
         C4 = Sip(5)
         C5 = Sip(6)
         C6 = Sip(7)
      else if (elem.eq.16) then ! Sulfur
         C0 = Sp(1)
         C1 = Sp(2)
         C2 = Sp(3)
         C3 = Sp(4)
         C4 = Sp(5)
         C5 = Sp(6)
         C6 = Sp(7)
      else if (elem.eq.26) then ! Iron
         C0 = Fep(1)
         C1 = Fep(2)
         C2 = Fep(3)
         C3 = Fep(4)
         C4 = Fep(5)
         C5 = Fep(6)
         C6 = Fep(7)
      else
         DUST_RECOMBINATION = 0.d0
         return
      end if

      phi = G * SQRT(T) / ne ! units K^1/2 cm^3

      a1 = C1 * (phi**C2)
      a2 = C3 * (T**C4)
      a3 = (-1.d0 * C5) - (C6 * LOG(T))
      DUST_RECOMBINATION = 1.d-14 * C0
      ! Seem to get a floating point exception when sometimes when a1 = 0
      DUST_RECOMBINATION = DUST_RECOMBINATION / (1.d0 + (a1 * (1.d0 + (a2 * (phi**a3)))))
      DUST_RECOMBINATION = DUST_RECOMBINATION * dr_sf

      ! Rescale the dust recombination rates due to a PAH normalization issue
      ! Zubko assumes 3.3e-5 C in PAH / H atom
      ! Weingartner & Draine assume 6e-5 C in PAH / H atom
      DUST_RECOMBINATION = DUST_RECOMBINATION * (3.3d-5 / 6.0d-5)

      return
   END FUNCTION DUST_RECOMBINATION

   FUNCTION comp_cr_co(G_0, xi_cr_H2, nCII, nH2, xO, n)
      implicit none
      real(dp):: comp_cr_co
      real(dp), intent(in):: G_0, xi_cr_H2, nCII, nH2, xO, n
      real(dp):: k0, k1, gammaCHx, gammaCHx_cr, beta

      ! CO creation

      ! see glover 2012
      k0 = 5.d-16 ! cm^3 s^-1
      k1 = 5.d-10 ! Rate coefficient for the formation of CO from O + CHx

      ! Assuming CO formation is modulated by CH2+, https://home.strw.leidenuniv.nl/~ewine/photo/display_ch2+_65e31a07e69e64dbd64d37801983018f.html
      gammaCHx_cr = 8.88d-15 * (xi_cr_H2 / 1.d-16) ! Cosmic rays
      gammaCHx = (1.41d-10 * G_0) + gammaCHx_cr
      
      beta = k1 * xO/(k1*xO + gammaCHx/n)
      comp_cr_co = k0 * nCII * nH2 * beta
   END FUNCTION comp_cr_co

   FUNCTION comp_de_co(G_0, xi_cr_H2)
      implicit none
      real(dp):: comp_de_co
      real(dp), intent(in):: G_0, xi_cr_H2
      real(dp):: gammaCO, gammaCO_cr

      ! CO destruction
      gammaCO = 2.43d-10 * G_0
      gammaCO_cr = 4.62d-15 * (xi_cr_H2 / 1d-16)

      comp_de_co = gammaCO + gammaCO_cr
   END FUNCTION comp_de_co

   SUBROUTINE update_co_chem(cr, de, ddt, n_carbon_ions, n_oxygen_ions, nco_mol, nCarbon, dxC, nOxygen, dxO)
      implicit none
      real(dp), intent(in):: cr, de, ddt
      integer, intent(in):: n_carbon_ions, n_oxygen_ions
      ! variables to update
      real(dp), dimension(n_carbon_ions):: dxC
      real(dp), dimension(n_oxygen_ions):: dxO
      real(dp):: nco_mol, nCarbon, nOxygen
      integer:: im
      real(dp):: cr_rxn, de_rxn, adj_rxn, nC_old, nO_old, nCO_old ! intermediate variables in computation, cr_unit is how much reaction happened
      real(dp):: min_n=1.d-30 !HK note: raised from 10^-40 to prevent crashed on dial
      ! adj_rxn: Amount of reaction to adjust due to extremely small or negative values

      ! Updates the number densities of CO, C+, C and O.
      ! de here is the rate after the step

      ! Update CO number density using the semi-implicit method
      adj_rxn = 0.d0
      nCO_old = nco_mol
      nco_mol = (nCO_old + cr*ddt) / (1.d0 + de*ddt)
      !nco_mol = MAX(nco_mol, 0.d0) ! Correct for negative values, 1d-40 small enough because lower limit for metallicity is 1d-8 solar fraction (~5d-7)
      ! Semi-implicit => cr uses the old value and de uses the new value
      cr_rxn = cr*ddt ! Old value, as passed in
      de_rxn = de*ddt*nco_mol ! New value as nco_mol is after the update
      ! Be consistent, use the change in nco_mol for the change of C and O as well
      nC_old = nCarbon
      nO_old = nOxygen
      nCarbon = nC_old - cr_rxn + de_rxn
      nOxygen = nO_old - cr_rxn + de_rxn

      ! Correct small or negative values
      if ((nCarbon .lt. min_n) .or. (nOxygen .lt. min_n)) then
         adj_rxn = min_n - MIN(nCarbon, nOxygen)
         nCarbon = nCarbon + adj_rxn
         nOxygen = nOxygen + adj_rxn
      end if

      do im=1, n_carbon_ions
         if (im .eq. 1) then
            dxC(1) = (dxC(1)*nC_old + de_rxn) / nCarbon ! Here nCarbon is new. CO destruction produces CI
         else if (im .eq. 2) then
            dxC(2) = (dxC(2)*nC_old - cr_rxn) / nCarbon ! CO formation consumes CII
         else
            dxC(im) = (dxC(im)*nC_old) / nCarbon
         end if
         dxC(im) = MIN(MAX(dxC(im), min_n), 1.)
      end do

      ! Ensure ionization fractions sum to 1
      do im=1, n_carbon_ions
         dxC(im) = dxC(im) / SUM(dxC)
      end do

      do im=1, n_oxygen_ions
         if (im .eq. 1) then
            dxO(1) = (dxO(1)*nO_old - cr_rxn + de_rxn) / nOxygen ! Here nOxygen is new
         else
            dxO(im) = (dxO(im)*nO_old) / nOxygen
         end if
         dxO(im) = MIN(MAX(dxO(im), min_n), 1.)
      end do

      ! Ensure ionization fractions sum to 1
      do im=1, n_oxygen_ions
         dxO(im) = dxO(im) / SUM(dxO)
      end do

      ! Ensure mass is conserved after correcting negatives
      nco_mol = nCO_old + cr_rxn - de_rxn - adj_rxn

   END SUBROUTINE update_co_chem

   FUNCTION comp_SCO(nco_mol, nh2, dx_SS)
      ! Returns the self shielding factor for CO
      implicit none
      real(dp), intent(in):: nco_mol, nh2, dx_SS
      real(dp), dimension(52, 2):: sco_table ! self shielding coefficient vs. CO column density
      real(dp), dimension(43, 2):: sh2_table ! self shielding coefficient vs. CO column density
      real(dp):: comp_SCO, logsCO, effcNCO, logeffcNCO
      real(dp):: logsH2, effcNH2, logeffcNH2
      integer:: i, idxCO, idxH2 ! Lower closest index

      sco_table(:,1) = (/ &
      1.000d+00, 1.000d+12, 1.650d+12, 2.995d+12, 5.979d+12, 1.313d+13, &
      3.172d+13, 8.429d+13, 2.464d+14, 7.923d+14, 1.670d+15, 2.595d+15, &
      4.435d+15, 6.008d+15, 8.952d+15, 1.334d+16, 1.661d+16, 2.274d+16, &
      3.115d+16, 4.266d+16, 5.843d+16, 8.002d+16, 1.096d+17, 1.501d+17, &
      2.055d+17, 2.815d+17, 4.241d+17, 6.389d+17, 9.625d+17, 1.450d+18, &
      2.184d+18, 3.291d+18, 4.124d+18, 5.685d+18, 7.838d+18, 1.080d+19, &
      1.285d+19, 1.681d+19, 2.199d+19, 2.538d+19, 3.222d+19, 4.091d+19, &
      5.193d+19, 5.893d+19, 7.356d+19, 8.269d+19, 9.246d+19, 1.031d+20, &
      1.148d+20, 1.277d+20, 1.419d+20, 1.578d+20 /)

      sco_table(:,2) = (/ &
      1.000d+00, 9.990d-01, 9.981d-01, 9.961d-01, 9.912d-01, 9.815d-01, &
      9.601d-01, 9.113d-01, 8.094d-01, 6.284d-01, 4.808d-01, 3.889d-01, &
      2.827d-01, 2.293d-01, 1.695d-01, 1.224d-01, 1.017d-01, 7.764d-02, &
      5.931d-02, 4.546d-02, 3.506d-02, 2.728d-02, 2.143d-02, 1.700d-02, &
      1.360d-02, 1.094d-02, 8.273d-03, 6.283d-03, 4.773d-03, 3.611d-03, &
      2.704d-03, 1.986d-03, 1.657d-03, 1.258d-03, 9.332d-04, 6.745d-04, &
      5.596d-04, 4.123d-04, 2.982d-04, 2.490d-04, 1.827d-04, 1.324d-04, &
      9.473d-05, 7.891d-05, 5.668d-05, 4.732d-05, 3.967d-05, 3.327d-05, &
      2.788d-05, 2.331d-05, 1.944d-05, 1.619d-05 /)

      sh2_table(:,1) = (/ &
      1.000d+00, 2.666d+13, 3.801d+14, 6.634d+15, 8.829d+16, 9.268d+17, &
      1.007d+18, 2.021d+18, 3.036d+18, 4.051d+18, 5.066d+18, 6.082d+18, &
      7.097d+18, 8.112d+18, 9.341d+18, 1.014d+19, 2.030d+19, 3.045d+19, &
      4.061d+19, 5.076d+19, 6.092d+19, 7.107d+19, 8.123d+19, 9.353d+19, &
      1.015d+20, 2.031d+20, 3.047d+20, 4.062d+20, 5.078d+20, 6.094d+20, &
      7.109d+20, 8.125d+20, 9.355d+20, 1.016d+21, 2.031d+21, 3.047d+21, &
      4.063d+21, 5.078d+21, 6.094d+21, 7.110d+21, 8.125d+21, 9.355d+21, &
      1.016d+22 /) 

      sh2_table(:,2) = (/ &
      1.000d+00, 9.999d-01, 9.893d-01, 9.678d-01, 9.465d-01, 9.137d-01, &
      9.121d-01, 8.966d-01, 8.862d-01, 8.781d-01, 8.716d-01, 8.660d-01, &
      8.612d-01, 8.569d-01, 8.524d-01, 8.497d-01, 8.262d-01, 8.118d-01, &
      8.011d-01, 7.921d-01, 7.841d-01, 7.769d-01, 7.702d-01, 7.626d-01, &
      7.579d-01, 7.094d-01, 6.712d-01, 6.378d-01, 6.074d-01, 5.791d-01, &
      5.524d-01, 5.271d-01, 4.977d-01, 4.793d-01, 2.837d-01, 1.526d-01, &
      7.774d-02, 3.952d-02, 2.093d-02, 1.199d-02, 7.666d-03, 5.333d-03, &
      4.666d-03 /) 

      ! initialize to 1.0
      comp_SCO = 1.d0

      ! Snapped to table at lower boundary, extrapolated above upper boundary
      effcNCO = MAX(nco_mol*dx_SS, sco_table(1,1))
      if (effcNCO.ge.1.578d+20) then 
         comp_SCO = comp_SCO * 1.619d-05
      else
         ! Find the lower closest index
         idxCO = 1
         do i=1, 51 ! If above upper boundary then extrapolate using the slope 
            !between the last two points, so conveniently we set i to not go to 52
            if (effcNCO .ge. sco_table(i,1)) idxCO = i
         end do

         logeffcNCO = log10(effcNCO)
         ! Interpolate or extrapolate automatically
         logsCO = log10(sco_table(idxCO,2)) &
         + (log10(sco_table(idxCO+1,2)) - log10(sco_table(idxCO,2))) &
         * (logeffcNCO - log10(sco_table(idxCO,1))) &
         / (log10(sco_table(idxCO+1,1)) - log10(sco_table(idxCO,1)))

         comp_SCO = comp_SCO * (10.d0 ** logsCO)
      end if


      ! Snapped to table at lower boundary, extrapolated above upper boundary
      effcNH2 = MAX(nh2*dx_SS, sco_table(1,1))
      if (effcNH2.ge.1.016d+22) then 
         comp_SCO = comp_SCO * 4.666d-03
      else 
         ! Find the lower closest index
         idxH2 = 1
         do i=1, 43 ! If above upper boundary then extrapolate using the slope 
            !between the last two points, so conveniently we set i to not go to 52
            if (effcNH2 .ge. sh2_table(i,1)) idxH2 = i
         end do

         logeffcNH2 = log10(effcNH2)
         ! Interpolate or extrapolate automatically
         logsH2 = log10(sh2_table(idxH2,2)) &
         + (log10(sh2_table(idxH2+1,2)) - log10(sh2_table(idxH2,2))) &
         * (logeffcNH2 - log10(sh2_table(idxH2,1))) &
         / (log10(sh2_table(idxH2+1,1)) - log10(sh2_table(idxH2,1)))

         comp_SCO = comp_SCO * (10.d0 ** logsH2)
      end if

   END FUNCTION comp_SCO

   FUNCTION comp_co_cooling(n, nH2, nHI, nCO, T)
      implicit none
      real(dp):: comp_co_cooling
      real(dp), intent(in):: n, nH2, nCO, T, nHI
      real(dp):: rot, vib_H, vib_H2, k_B, T_pivot

      ! CO cooling

      k_B = 1.3806d-16
      T_pivot = 3080.d0
      rot = n * nCO * 4.d0*((k_B*T)**2.d0)*9.7d-8 / (n * 2.76d0 * k_B * (1.d0 + (3.3d6 * (T/1.d3)**(3.d0/4.d0)/n) + 1.5d0*((3.3d6 * (T/1.d3)**(3.d0/4.d0)/n)**0.5d0)))

      vib_H2 = nH2 * nCO * T_pivot * k_B * 4.3d-14 * T * EXP(-(3.14d5/T)**0.333d0) * EXP(-T_pivot/T)
      vib_H = nHI * nCO * T_pivot * k_B * 3.0d-12 * (T**0.5d0) * EXP(-(2000.d0/T)**3.43d0) * EXP(-T_pivot/T)

      comp_co_cooling = rot + vib_H2 + vib_H
   END FUNCTION comp_co_cooling

END MODULE coolrates_module
