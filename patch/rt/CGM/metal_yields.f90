module metal_yields
   use amr_commons, only: dp, yields_lc18_dir
   implicit none

   !---------------------------------------
   ! Facilitates yield tables and interpolation function. Includes
   ! yield tables for C,N,O,Al,Si,Ca and Fe from AGB winds and SNII.
   ! Yield tables taken from NuGrid II, Ritter et al. (2018).
   !---------------------------------------
   integer, parameter :: nmet = 5
   integer, parameter :: nmass_AGB = 8
   integer, parameter :: nmass_SNII = 4
   ! For Pop 3 SN
   integer, parameter :: nmass_POP3_SNII = 7
   integer, parameter :: nmass_POP3_HN = 4
   integer, parameter :: nmass_POP3_HMHN = 6 !For Nomoto !14 For Heger
   ! For type II SN alt yields
   integer, parameter :: nmet_portinari = 5
   integer, parameter :: nmass_portinari = 9
   ! For type II SN Limongi+Chieffi 2018 yields
   integer, parameter :: nmet_lc18 = 4
   integer, parameter :: nmass_lc18 = 9
   ! For type II HN Nomoto 2013 yields
   integer, parameter :: nmet_N13_HN = 5
   integer, parameter :: nmass_N13_HN = 4

   ! New yields taken from Seitenzahl et al. (2013). Sum of all isotopes using model N100, see table 2.
   ! Compare to old yields (Raiteri et al. 1996) 0.63 Msun Fe + 0.13Msun O (e.g. Theilemann 1986)
   ! Mass (Msun) of enriched material from SNIa
   real(dp) :: SNIaFe = 7.40d-01  ! Iron; old value 0.7d0*scale_m  !see Hopkins, 0.63d0 Raiteri
   real(dp) :: SNIaO = 1.01d-01   ! Oxygen, old value 0.13d0*scale_m
   real(dp) :: SNIaC = 3.04d-03   ! Carbon
   real(dp) :: SNIaN = 3.21d-06   ! Nitrogen
   real(dp) :: SNIaMg = 1.54d-02  ! Magnesium
   real(dp) :: SNIaAl = 6.74d-04  ! Aluminium
   real(dp) :: SNIaSi = 2.87d-01  ! Silicon
   real(dp) :: SNIaCa = 1.48d-02  ! Calcium
   real(dp) :: SNIaEu = 0.d0      ! Europium (not known)
   real(dp) :: SNIaNe = 3.57d-03  ! Neon
   real(dp) :: SNIaS = 1.15d-01   ! Sulfur

   ! Mass (Msun) of enriched material from neutron star mergers.
   real(dp) :: MEuNSNS = 1.0d-5   ! Europium. See Cote et al. (2018).

   ! Tables from Nugrid.
   real(dp), dimension(1:nmet) :: table_met = (/0.0001, 0.001, 0.006, 0.01, 0.02/)
   real(dp), dimension(1:nmass_AGB) :: table_mass_AGB = (/1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0/)
   real(dp), dimension(1:nmass_SNII) :: table_mass_SNII = (/12.0, 15.0, 20.0, 25.0/)
   real(dp), dimension(1:nmass_POP3_SNII) :: table_mass_POP3_SNII = (/13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0/)
   real(dp), dimension(1:nmass_POP3_HN) :: table_mass_POP3_HN = (/20.0, 25.0, 30.0, 40.0/)
   ! Note that this represents the mass of the helium core and not the mass of the star
   !real(dp), dimension(1:nmass_POP3_HMHN) :: table_mass_POP3_HMHN = (/ 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0/)
   ! Updade from Nomoto 2013 --> based on umeda and nomoto 2002 --> this is the actual stellar mass
   real(dp), dimension(1:nmass_POP3_HMHN) :: table_mass_POP3_HMHN = (/ 140.0, 150.0, 170.0, 200.0, 270.0, 300.0 /)
   ! For portinari 1998 yields
   real(dp), dimension(1:nmet_portinari):: table_met_portinari = (/0.0004, 0.004, 0.008, 0.02, 0.05/)
   real(dp), dimension(1:nmass_portinari):: table_mass_portinari = (/9.0, 12.0, 15.0, 20.0, 30.0, 40.0, 60.0, 100.0, 120.0/)
   ! For Limongi+Chieffi 2018 yields
   real(dp), dimension(1:nmet_lc18):: table_met_lc18 = (/-3.d0, -2.d0, -1.d0, 0.d0/)
   real(dp), dimension(1:nmass_lc18):: table_mass_lc18 = (/13.d0, 15.d0, 20.d0, 25.d0, 30.d0, 40.d0, 60.d0, 80.d0, 120.d0/)
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_H_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_He_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_C_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_N_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_O_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_Ne_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_Mg_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_Si_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_S_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_Fe_SNII_lc18=0.d0
   real(dp), dimension(1:nmet_lc18,1:nmass_lc18)::table_yield_Ca_SNII_lc18=0.d0
   ! For Nomoto 2013 HN yields
   real(dp), dimension(1:nmet_N13_HN):: table_met_N13_HN = (/0.0001d0, 0.001d0, 0.004d0, 0.008d0, 0.02d0/)
   real(dp), dimension(1:nmass_N13_HN):: table_mass_N13_HN = (/20.0d0, 25.0d0, 30.0d0, 40.0d0/)

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_H_AGB = transpose(reshape( &
                                                                              (/ &
                                           3.28902E-1, 7.13303E-1, 9.42403E-1, 1.52300E+0, 2.04200E+0, 2.55100E+0, 3.05200E+0, &
                                           3.52200E+0, 3.02503E-1, 7.30504E-1, 9.50704E-1, 1.57000E+0, 2.07600E+0, 2.57200E+0, &
                                           3.06600E+0, 3.54600E+0, 3.10303E-1, 7.34104E-1, 9.62105E-1, 1.61300E+0, 2.14300E+0, & 
                                           2.65000E+0, 3.14200E+0, 3.63000E+0, 3.12304E-1, 7.26704E-1, 9.55105E-1, 1.60001E+0, & 
                                           2.13900E+0, 2.64700E+0, 3.13300E+0, 3.61300E+0, 3.01904E-1, 6.84904E-1, 9.23505E-1, &
                                           1.54301E+0, 2.11500E+0, 2.62200E+0, 3.10700E+0, 3.56800E+0 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_He_AGB = transpose(reshape( &
                                                                              (/ &
                                           1.30812E-1, 2.87504E-1, 3.58254E-1, 6.13467E-1, 1.03702E+0, 1.45101E+0, 1.82000E+0, &
                                           2.19205E+0, 1.11721E-1, 2.82679E-1, 3.74589E-1, 5.91991E-1, 1.01205E+0, 1.43401E+0, &
                                           1.80401E+0, 2.18204E+0, 1.19343E-1, 2.85150E-1, 3.86562E-1, 6.39135E-1, 9.66083E-1, &
                                           1.37103E+0, 1.77102E+0, 2.16802E+0, 1.23849E-1, 2.87662E-1, 3.93879E-1, 6.76727E-1, &
                                           9.63581E-1, 1.38507E+0, 1.80501E+0, 2.20101E+0, 1.29147E-1, 2.91066E-1, 4.05590E-1, &
                                           7.15030E-1, 9.66294E-1, 1.39725E+0, 1.81701E+0, 2.23901E+0 & 
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Fe_AGB = transpose(reshape( &
                                                                              (/ &
                                           6.453d-07, 1.451d-06, 1.881d-06, 3.073d-06, 4.429d-06, 5.754d-06, 7.012d-06, 8.238d-06, &
                                           5.977d-06, 1.473d-05, 1.931d-05, 3.116d-05, 4.448d-05, 5.764d-05, 7.016d-05, 8.259d-05, &
                                           3.733d-05, 8.917d-05, 1.186d-04, 1.963d-04, 2.701d-04, 3.493d-04, 4.266d-04, 5.037d-04, &
                                           3.169d-04, 7.400d-04, 9.961d-04, 1.676d-03, 2.260d-03, 2.933d-03, 3.590d-03, 4.227d-03, &
                                            6.325d-04, 1.437d-03, 1.982d-03, 3.379d-03, 4.543d-03, 5.909d-03, 7.231d-03, 8.526d-03 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_O_AGB = transpose(reshape( &
                                                                             (/ &
                                           1.835d-03, 8.626d-03, 9.952d-03, 2.781d-03, 4.080d-03, 1.833d-04, 2.239d-04, 3.249d-03, &
                                           5.767d-04, 5.937d-03, 1.184d-02, 4.024d-03, 6.280d-03, 1.217d-03, 1.252d-03, 3.756d-03, &
                                           1.921d-03, 7.652d-03, 1.330d-02, 1.577d-02, 1.502d-02, 1.426d-02, 1.410d-02, 1.550d-02, &
                                           2.127d-03, 5.908d-03, 1.389d-02, 2.232d-02, 1.982d-02, 1.997d-02, 1.943d-02, 1.835d-02, &
                                            4.243d-03, 1.001d-02, 1.911d-02, 3.604d-02, 3.844d-02, 3.923d-02, 3.799d-02, 4.214d-02 &
                                                                             /), [nmass_AGB, nmet] &
                                                                             ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_N_AGB = transpose(reshape( &
                                                                             (/ &
                                           7.763d-05, 5.710d-05, 3.870d-05, 3.408d-05, 1.019d-02, 4.692d-03, 2.360d-03, 8.018d-03, &
                                           1.985d-05, 1.036d-04, 1.711d-04, 2.213d-04, 7.812d-03, 6.323d-03, 5.150d-03, 7.922d-03, &
                                           8.937d-05, 4.067d-04, 7.848d-04, 1.694d-03, 2.854d-03, 7.439d-03, 1.453d-02, 1.699d-02, &
                                           3.519d-04, 1.248d-03, 1.931d-03, 3.817d-03, 5.838d-03, 1.559d-02, 2.974d-02, 3.151d-02, &
                                            6.713d-04, 2.328d-03, 3.561d-03, 7.173d-03, 1.041d-02, 1.678d-02, 4.112d-02, 4.746d-02 &
                                                                             /), [nmass_AGB, nmet] &
                                                                             ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Si_AGB = transpose(reshape( &
                                                                              (/ &
                                           7.645d-07, 1.791d-06, 2.559d-06, 3.970d-06, 6.714d-06, 7.136d-06, 1.495d-05, 9.418d-06, &
                                           6.568d-06, 1.642d-05, 2.203d-05, 3.545d-05, 5.208d-05, 6.518d-05, 8.007d-05, 9.254d-05, &
                                           4.094d-05, 9.836d-05, 1.316d-04, 2.184d-04, 3.014d-04, 3.845d-04, 4.729d-04, 5.591d-04, &
                                           1.787d-04, 4.176d-04, 5.643d-04, 9.529d-04, 1.293d-03, 1.671d-03, 2.055d-03, 2.421d-03, &
                                            3.568d-04, 8.105d-04, 1.121d-03, 1.919d-03, 2.603d-03, 3.362d-03, 4.098d-03, 4.834d-03 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Ca_AGB = transpose(reshape( &
                                                                              (/ &
                                           5.590d-08, 1.240d-07, 1.607d-07, 2.630d-07, 3.790d-07, 4.920d-07, 5.993d-07, 7.038d-07, &
                                           5.108d-07, 1.257d-06, 1.647d-06, 2.666d-06, 3.806d-06, 4.929d-06, 5.997d-06, 7.057d-06, &
                                           3.192d-06, 7.613d-06, 1.011d-05, 1.678d-05, 2.310d-05, 2.986d-05, 3.647d-05, 4.305d-05, &
                                           1.647d-05, 3.844d-05, 5.161d-05, 8.681d-05, 1.174d-04, 1.524d-04, 1.865d-04, 2.196d-04, &
                                            3.288d-05, 7.464d-05, 1.027d-04, 1.749d-04, 2.359d-04, 3.070d-04, 3.758d-04, 4.430d-04 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Al_AGB = transpose(reshape( &
                                                                              (/ &
                                           3.506d-08, 3.905d-07, 8.674d-07, 6.366d-07, 2.057d-06, 2.275d-06, 1.841d-06, 1.037d-06, &
                                           2.997d-07, 7.956d-07, 1.619d-06, 2.118d-06, 3.509d-06, 4.639d-06, 7.882d-06, 5.905d-06, &
                                           1.682d-06, 4.163d-06, 5.781d-06, 1.022d-05, 1.527d-05, 1.637d-05, 2.226d-05, 2.862d-05, &
                                           1.429d-05, 3.373d-05, 4.624d-05, 7.894d-05, 1.151d-04, 1.432d-04, 1.770d-04, 2.175d-04, &
                                            2.851d-05, 6.518d-05, 9.123d-05, 1.576d-04, 2.284d-04, 2.850d-04, 3.422d-04, 4.029d-04 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Mg_AGB = transpose(reshape( &
                                                                              (/ &
                                           1.005d-06, 5.403d-05, 9.078d-05, 3.960d-05, 1.193d-04, 2.429d-05, 1.937d-05, 4.907d-05, &
                                           7.271d-06, 3.665d-05, 1.105d-04, 6.605d-05, 1.570d-04, 9.500d-05, 1.083d-04, 1.736d-04, &
                                           4.351d-05, 1.191d-04, 1.905d-04, 3.405d-04, 4.288d-04, 4.160d-04, 5.799d-04, 1.017d-03, &
                                           1.657d-04, 3.927d-04, 5.835d-04, 1.099d-03, 1.549d-03, 1.784d-03, 2.296d-03, 3.465d-03, &
                                            3.308d-04, 7.573d-04, 1.126d-03, 2.166d-03, 3.451d-03, 3.713d-03, 4.115d-03, 4.980d-03 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Eu_AGB = transpose(reshape( &
                                                                              (/ &
                                           3.857D-12, 6.504D-13, 9.461D-13, 1.386D-12, 2.736D-12, 1.848D-12, 1.967D-12, 2.290D-12, &
                                           1.927D-12, 4.903D-12, 7.150D-12, 9.971D-12, 2.064D-11, 1.750D-11, 1.980D-11, 2.305D-11, &
                                           1.113D-11, 2.769D-11, 3.974D-11, 6.642D-11, 9.837D-11, 1.028D-10, 1.263D-10, 1.434D-10, &
                                           9.444D-11, 2.202D-10, 3.170D-10, 5.874D-10, 7.074D-10, 8.725D-10, 1.037D-09, 1.202D-09, &
                                            1.886D-10, 4.264D-10, 5.879D-10, 1.051D-09, 1.432D-09, 1.779D-09, 2.114D-09, 2.441D-09 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_C_AGB = transpose(reshape( &
                                                                             (/ &
                                           8.708D-03, 2.147D-02, 2.356D-02, 8.884D-03, 1.988D-03, 7.856D-04, 2.937D-04, 1.797D-03, &
                                           7.063D-04, 1.337D-02, 2.608D-02, 8.723D-03, 4.042D-03, 1.234D-03, 4.238D-04, 1.284D-03, &
                                           2.736D-04, 7.715D-03, 1.786D-02, 1.856D-02, 9.100D-03, 8.749D-04, 1.498D-03, 1.240D-03, &
                                           7.294D-04, 4.655D-03, 1.852D-02, 2.990D-02, 1.991D-02, 9.310D-03, 3.304D-03, 2.154D-03, &
                                            1.353D-03, 4.750D-03, 1.894D-02, 3.956D-02, 3.097D-02, 2.335D-02, 1.977D-03, 2.462D-03 &
                                                                             /), [nmass_AGB, nmet] &
                                                                             ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_Ne_AGB = transpose(reshape( &
                                                                              (/ &
                                           1.487D-05, 8.758D-04, 9.089D-04, 1.138D-04, 2.528D-04, 4.369D-05, 5.373D-05, 1.268D-04, &
                                           3.097D-05, 4.082D-04, 1.256D-03, 2.199D-04, 3.990D-04, 2.592D-04, 3.158D-04, 4.209D-04, &
                                           1.531D-04, 5.682D-04, 1.258D-03, 1.347D-03, 1.264D-03, 1.448D-03, 1.798D-03, 2.330D-03, &
                                           4.362D-04, 1.103D-03, 2.174D-03, 3.613D-03, 3.518D-03, 4.214D-03, 5.218D-03, 6.676D-03, &
                                            8.662D-04, 2.100D-03, 3.795D-03, 7.574D-03, 7.390D-03, 8.558D-03, 9.959D-03, 1.186D-02 &
                                                                              /), [nmass_AGB, nmet] &
                                                                              ))

   real(dp), dimension(1:nmet, 1:nmass_AGB) :: table_yield_S_AGB = transpose(reshape( &
                                                                             (/ &
                                           5.121D-07, 1.134D-06, 1.469D-06, 2.391D-06, 3.455D-06, 4.471D-06, 5.436D-06, 6.386D-06, &
                                           4.631D-06, 1.144D-05, 1.503D-05, 2.419D-05, 3.456D-05, 4.472D-05, 5.438D-05, 6.401D-05, &
                                           2.892D-05, 6.912D-05, 9.196D-05, 1.524D-04, 2.096D-04, 2.707D-04, 3.308D-04, 3.904D-04, &
                                           9.328D-05, 2.178D-04, 2.931D-04, 4.934D-04, 6.663D-04, 8.642D-04, 1.058D-03, 1.245D-03, &
                                            1.862D-04, 4.229D-04, 5.830D-04, 9.942D-04, 1.340D-03, 1.741D-03, 2.130D-03, 2.511D-03 &
                                                                             /), [nmass_AGB, nmet] &
                                                                             ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Fe_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
                  0.1543680d0, 0.1567630d0, 0.0956094d0, 0.156620d0, 1.0000d-40, 1.00000d-40, 1.0000d-40, 1.00000d-40, 1.0000d-40, &
                  0.1668890d0, 0.1663360d0, 0.1223160d0, 0.187550d0, 1.0000d-40, 1.00000d-40, 1.0000d-40, 1.00000d-40, 1.0000d-40, &
                  0.1354800d0, 0.1534930d0, 0.1134830d0, 0.185403d0, 1.0000d-40, 1.00000d-40, 1.0000d-40, 1.00000d-40, 1.0000d-40, &
                  0.0477977d0, 0.1081770d0, 0.0825510d0, 0.165255d0, 1.0000d-40, 0.0883201d0, 0.166038d0, 0.1362490d0, 0.169388d0, &
                   0.0435565d0, 0.0984339d0, 0.0930222d0, 0.100698d0, 0.109369d0, 0.0902536d0, 0.067424d0, 0.0697738d0, 0.070170d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Fe_SNII = transpose(reshape( &
                                                                                (/ &
                                                                                2.443d-01, 5.073d-02, 7.646d-05, 2.518d-05, &
                                                                                1.692d-01, 5.407d-02, 4.089d-02, 2.326d-04, &
                                                                                1.863d-01, 1.739d-01, 3.275d-01, 1.036d-03, &
                                                                                2.032d-01, 7.515d-02, 2.745d-01, 8.288d-03, &
                                                                                1.531d-01, 5.561d-02, 1.493d-02, 7.470d-03 &
                                                                                /), [nmass_SNII, nmet] &
                                                                                ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_O_SNII_portinari = transpose(reshape( &
                                                                                                        (/ &
                           0.1384110d0, 0.668308d0, 1.40365d0, 3.03062d0, 3.14628d0, 0.818882d0, 0.280666d0, 13.9405d0, 19.1449d0, &
                           0.1454280d0, 0.684527d0, 1.27809d0, 2.91323d0, 2.97717d0, 2.061110d0, 3.430990d0, 6.84938d0, 1.000d-40, &
                           0.1879840d0, 0.651678d0, 1.25407d0, 2.77489d0, 2.96329d0, 2.096200d0, 4.548680d0, 7.36338d0, 8.70540d0, &
                           0.0964619d0, 0.600155d0, 1.13413d0, 2.57454d0, 2.83490d0, 1.652030d0, 2.051050d0, 3.23853d0, 2.27190d0, &
                            0.0801391d0, 0.617075d0, 1.22322d0, 2.97120d0, 2.34298d0, 0.305554d0, 0.194200d0, 1.000d-40, 1.000d-40 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                        ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_O_SNII = transpose(reshape( &
                                                                               (/ &
                                                                               2.828d-01, 1.176d+00, 2.043d+00, 7.062d-01, &
                                                                               3.459d-01, 1.148d+00, 1.871d+00, 7.875d-01, &
                                                                               1.879d-01, 8.612d-01, 1.442d+00, 8.904d-01, &
                                                                               1.194d-01, 8.829d-01, 1.314d+00, 8.494d-01, &
                                                                               1.205d-01, 8.970d-01, 1.504d+00, 8.446d-01 &
                                                                               /), [nmass_SNII, nmet] &
                                                                               ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_N_SNII_portinari = transpose(reshape( &
                                                                                                        (/ &
       0.00074145d0, 0.00102282d0, 0.0011694d0, 0.00173721d0, 0.00276251d0, 0.00423595d0, 0.00669394d0, 0.00907829d0, 0.0136157d0, &
       0.00669904d0, 0.00894477d0, 0.0119490d0, 0.01420710d0, 0.02556370d0, 0.03378440d0, 0.05547040d0, 0.10843000d0, 0.2024690d0, &
       0.01088610d0, 0.01796560d0, 0.0226703d0, 0.02778760d0, 0.04645690d0, 0.06326420d0, 0.10625300d0, 0.22255900d0, 0.2943540d0, &
       0.03060460d0, 0.04145320d0, 0.0504746d0, 0.06619160d0, 0.09965390d0, 0.16430400d0, 0.27158400d0, 0.56085900d0, 0.8672890d0, &
        0.06405570d0, 0.08336080d0, 0.0954386d0, 0.11948600d0, 0.22624500d0, 0.36434700d0, 0.87646500d0, 1.48248000d0, 1.7473700d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                        ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_N_SNII = transpose(reshape( &
                                                                               (/ &
                                                                               2.160d-04, 3.047d-04, 2.320d-02, 2.813d-02, &
                                                                               2.047d-03, 2.664d-03, 3.984d-03, 5.828d-03, &
                                                                               1.172d-02, 1.417d-02, 2.151d-02, 2.986d-02, &
                                                                               1.993d-02, 2.453d-02, 3.150d-02, 4.422d-02, &
                                                                               4.055d-02, 4.235d-02, 5.627d-02, 3.436d-02 &
                                                                               /), [nmass_SNII, nmet] &
                                                                               ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Si_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
            0.0306685d0, 0.0861459d0, 0.1994720d0, 0.339698d0, 0.0117824d0, 1.0000d-40, 1.000000d-40, 1.0000000d-40, 1.000000d-40, &
            0.0364149d0, 0.0908148d0, 0.1922690d0, 0.371644d0, 0.0202428d0, 1.0000d-40, 1.000000d-40, 1.0000000d-40, 1.000000d-40, &
            0.0268590d0, 0.0433140d0, 0.0971206d0, 0.191930d0, 0.0004858d0, 1.0000d-40, 0.00112705d0, 0.000632674d0, 1.1759400d-5, &
            0.0711394d0, 0.1040840d0, 0.2041500d0, 0.323011d0, 0.0258496d0, 0.273618d0, 0.32267100d0, 0.306979000d0, 0.32529500d0, &
             0.0800406d0, 0.1066310d0, 0.2639140d0, 0.269784d0, 0.2438210d0, 0.110622d0, 0.12176000d0, 0.120613000d0, 0.12042000d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Si_SNII = transpose(reshape( &
                                                                                (/ &
                                                                                9.289d-02, 1.078d-01, 1.138d-01, 2.063d-02, &
                                                                                9.255d-02, 9.346d-02, 2.482d-01, 3.416d-02, &
                                                                                8.842d-02, 1.468d-01, 3.725d-01, 2.865d-02, &
                                                                                8.212d-02, 2.248d-01, 8.855d-02, 3.778d-02, &
                                                                                5.401d-02, 1.902d-01, 1.487d-01, 5.818d-02 &
                                                                                /), [nmass_SNII, nmet] &
                                                                                ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Ca_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
       0.00258995d0, 0.00789402d0, 0.01566880d0, 0.0163219d0, 1.000000d-40, 1.00000d-40, 1.000000d-40, 1.000000d-40, 1.000000d-40, &
       0.00361655d0, 0.00673503d0, 0.01379450d0, 0.0200204d0, 1.4576300d-6, 1.00000d-40, 1.000000d-40, 1.000000d-40, 1.000000d-40, &
       0.00550154d0, 0.00729800d0, 0.01300520d0, 0.0196451d0, 1.000000d-40, 1.00000d-40, 1.000000d-40, 1.000000d-40, 1.000000d-40, &
       0.01134290d0, 0.00974814d0, 0.00953551d0, 0.0157458d0, 1.000000d-40, 0.0130499d0, 0.01523590d0, 0.01545680d0, 0.01578700d0, &
        0.01307380d0, 0.00897471d0, 0.01230730d0, 0.0130902d0, 0.01291580d0, 0.0082617d0, 0.00627183d0, 0.00647664d0, 0.00651119d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Ca_SNII = transpose(reshape( &
                                                                                (/ &
                                                                                8.686d-03, 5.152d-03, 1.363d-03, 2.154d-06, &
                                                                                6.655d-03, 5.118d-03, 1.390d-02, 3.186d-05, &
                                                                                8.256d-03, 7.418d-03, 2.550d-02, 9.611d-05, &
                                                                                8.661d-03, 1.717d-03, 5.985d-04, 7.348d-04, &
                                                                                5.745d-03, 2.537d-03, 1.102d-03, 4.162d-04 &
                                                                                /), [nmass_SNII, nmet] &
                                                                                ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Al_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                            1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Al_SNII = transpose(reshape( &
                                                                                (/ &
                                                                                6.219d-05, 1.157d-03, 1.505d-03, 1.315d-04, &
                                                                                2.466d-04, 9.427d-04, 1.540d-03, 3.695d-04, &
                                                                                2.530d-04, 1.299d-03, 1.517d-03, 1.826d-03, &
                                                                                4.564d-04, 4.246d-03, 3.031d-03, 2.227d-03, &
                                                                                7.659d-04, 2.780d-03, 5.394d-03, 3.650d-03 &
                                                                                /), [nmass_SNII, nmet] &
                                                                                ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Mg_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
             0.00553325d0, 0.0289872d0, 0.0266895d0, 0.0757126d0, 0.1294020d0, 0.0540950d0, 1.00000d-40, 0.7173110d0, 0.7173110d0, &
             0.00592746d0, 0.0381874d0, 0.0573169d0, 0.0460059d0, 0.1291610d0, 0.2203100d0, 0.2949440d0, 0.3731720d0, 1.00000d-40, &
             0.00223265d0, 0.1124520d0, 0.1351610d0, 0.0958842d0, 0.4875820d0, 0.6391270d0, 0.4889250d0, 0.7503020d0, 0.8616190d0, &
             0.00300958d0, 0.0227667d0, 0.0299653d0, 0.0606579d0, 0.1128560d0, 0.0207418d0, 0.0286330d0, 0.0842471d0, 0.0572967d0, &
              0.00252206d0, 0.0259101d0, 0.0151747d0, 0.1036030d0, 0.0890234d0, 0.0301898d0, 0.0421334d0, 0.0409041d0, 0.0406968d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Mg_SNII = transpose(reshape( &
                                                                                (/ &
                                                                                1.638D-02, 7.562D-02, 1.160D-01, 2.178D-02, &
                                                                                2.463D-02, 8.823D-02, 1.019D-01, 2.791D-02, &
                                                                                1.122D-02, 6.267D-02, 7.786D-02, 4.062D-02, &
                                                                                1.109D-02, 6.908D-02, 7.642D-02, 4.494D-02, &
                                                                                1.314D-02, 3.892D-02, 1.366D-01, 9.823D-02 &
                                                                                /), [nmass_SNII, nmet] &
                                                                                ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Eu_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                           1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, &
                                                            1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40, 1.d-40 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_Eu_SNII = transpose(reshape( &
                                                                                (/ &
                                                                                3.856D-12, 4.834D-12, 5.922D-12, 7.524D-12, &
                                                                                3.948D-11, 4.893D-11, 6.355D-11, 8.307D-11, &
                                                                                2.071D-10, 2.645D-10, 3.369D-10, 4.567D-10, &
                                                                                1.804D-09, 2.176D-09, 2.504D-09, 2.892D-09, &
                                                                                3.718D-09, 3.717D-09, 4.520D-09, 2.807D-09 &
                                                                                /), [nmass_SNII, nmet] &
                                                                                ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_C_SNII_portinari = transpose(reshape( &
                                                                                                        (/ &
                    0.0889314d0, 0.1689870d0, 0.212393d0, 0.282147d0, 0.348097d0, 0.336192d0, 0.506284d0, 0.712529d0, 1.9570400d0, &
                    0.0827334d0, 0.1606370d0, 0.226899d0, 0.253049d0, 0.295797d0, 0.334057d0, 0.396579d0, 0.312612d0, 0.0430967d0, &
                    0.2438400d0, 0.1468860d0, 0.221549d0, 0.235736d0, 0.267712d0, 0.298049d0, 7.222610d0, 14.97410d0, 19.730800d0, &
                    0.0346755d0, 0.1177890d0, 0.201324d0, 0.190786d0, 0.193987d0, 2.840130d0, 4.643580d0, 10.91170d0, 8.0469500d0, &
                     0.0160955d0, 0.0951247d0, 0.168255d0, 0.158531d0, 0.818138d0, 1.864640d0, 2.935790d0, 2.574210d0, 2.3565600d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                        ))

   real(dp), dimension(1:nmet, 1:nmass_SNII) :: table_yield_C_SNII = transpose(reshape( &
                                                                               (/ &
                                                                               1.226D-01, 1.646D-01, 1.956D-01, 2.037D-01, &
                                                                               1.200D-01, 1.537D-01, 2.222D-01, 2.115D-01, &
                                                                               1.015D-01, 1.655D-01, 1.473D-01, 2.206D-01, &
                                                                               8.914D-02, 1.575D-01, 2.294D-01, 2.238D-01, &
                                                                               7.646D-02, 1.472D-01, 1.923D-01, 2.035D-01 &
                                                                               /), [nmass_SNII, nmet] &
                                                                               ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_S_SNII_portinari = transpose(reshape( &
                                                                                                        (/ &
          0.0137095d0, 0.0437912d0, 0.1174880d0, 0.171796d0, 0.000107485d0, 1.00000d-40, 1.000000d-40, 1.0000000d-40, 1.00000d-40, &
          0.0184072d0, 0.0399431d0, 0.0958444d0, 0.198942d0, 0.000553066d0, 1.00000d-40, 1.000000d-40, 1.0000000d-40, 1.00000d-40, &
          0.0268590d0, 0.0433140d0, 0.0971206d0, 0.191930d0, 0.000485800d0, 1.00000d-40, 0.00112705d0, 0.000632674d0, 1.175940d-5, &
          0.0592021d0, 0.0569651d0, 0.0893985d0, 0.149099d0, 1.0000000d-40, 0.1392900d0, 0.15801700d0, 0.137304000d0, 0.1507800d0, &
           0.0678521d0, 0.0547168d0, 0.1218850d0, 0.114363d0, 0.105674000d0, 0.0530726d0, 0.04848390d0, 0.048956200d0, 0.0490359d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                        ))

   real(dp), dimension(1:nmet_portinari, 1:nmass_portinari) :: table_yield_Ne_SNII_portinari = transpose(reshape( &
                                                                                                         (/ &
               0.00919520d0, 0.1079110d0, 0.0987822d0, 0.0993798d0, 0.615611d0, 0.270494d0, 0.000372705d0, 1.001010d0, 1.327430d0, &
               0.00592746d0, 0.0381874d0, 0.0573169d0, 0.0460059d0, 0.129161d0, 0.220310d0, 0.294944000d0, 0.373172d0, 1.0000d-40, &
               0.00629048d0, 0.1124520d0, 0.1351610d0, 0.0958842d0, 0.487582d0, 0.639127d0, 0.488925000d0, 0.750302d0, 0.861619d0, &
               0.00763224d0, 0.1067220d0, 0.1546640d0, 0.1854760d0, 0.464056d0, 0.351237d0, 0.373631000d0, 0.898500d0, 0.686845d0, &
                0.00529803d0, 0.1275810d0, 0.0724648d0, 0.3425550d0, 0.549435d0, 0.780050d0, 0.758300000d0, 1.437560d0, 1.423780d0 &
                                                                                             /), [nmass_portinari, nmet_portinari] &
                                                                                                         ))

   ! Harley edits for HN yields
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_C_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               1.24020E-1, 1.94098E-1, 1.05082E-1, 5.10820E-2, 1.24020E-1, 1.94098E-1, 1.05082E-1, 5.10820E-2, 8.34910E-2, &
               1.28383E-1, 1.36339E-1, 3.73368E-1, 1.24281E-1, 1.38232E-1, 1.53815E-1, 4.19894E-1, 2.10450E-1, 2.10600E-1, &
               1.80920E-1, 4.90431E-1 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_N_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               1.29014E-2, 9.20724E-3, 6.19091E-3, 8.53317E-3, 1.29014E-2, 9.20724E-3, 6.19091E-3, 8.53317E-3, 1.84274E-2, &
               3.15945E-2, 2.01056E-2, 2.52135E-2, 3.31614E-2, 4.94910E-2, 4.04686E-2, 3.61321E-2, 7.21530E-2, 1.30600E-1, &
               1.02015E-1, 5.81426E-2 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_O_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               2.00003E+0, 3.70010E+0, 4.95007E+0, 6.42020E+0, 2.00003E+0, 3.70010E+0, 4.95007E+0, 6.42020E+0, 7.88771E-1, &
               2.07089E+0, 3.82016E+0, 6.80081E+0, 8.67309E-1, 2.12317E+0, 3.31164E+0, 6.90949E+0, 9.84929E-1, 2.26504E+0, &
               2.74447E+0, 7.06137E+0 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_Mg_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               2.33490E-1, 1.96280E-1, 3.17000E-1, 5.22558E-1, 2.33490E-1, 1.96280E-1, 3.17000E-1, 5.22558E-1, 8.43400E-2, &
               2.37050E-1, 2.01580E-1, 4.20850E-1, 8.24760E-2, 2.37809E-1, 2.00533E-1, 4.10900E-1, 8.75200E-2, 2.53400E-1, &
               2.13500E-1, 4.46600E-1 & 
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))       
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_Ne_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               4.57267E-1, 1.05185E+0, 1.05170E+0, 1.83147E-1, 4.57267E-1, 1.05185E+0, 1.05170E+0, 1.83147E-1, 1.43732E-1, &
               6.39897E-1, 4.95413E-1, 1.15166E+0, 1.97509E-1, 6.48734E-1, 5.18719E-1, 1.37931E+0, 3.02790E-1, 6.65300E-1, &
               5.54870E-1, 1.75589E+0 & 
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_S_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               3.69913E-2, 4.31296E-2, 9.39319E-2, 2.84941E-1, 3.69913E-2, 4.31296E-2, 9.39319E-2, 2.84941E-1, 3.46882E-2, &
               4.17573E-2, 1.59859E-1, 2.38752E-1, 4.11213E-2, 4.41473E-2, 1.41051E-1, 1.97745E-1, 5.18110E-2, 4.84893E-2, &
               1.23579E-1, 1.58714E-1 & 
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_Si_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               1.16350E-1, 1.13266E-1, 2.38930E-1, 7.28220E-1, 1.16350E-1, 1.13266E-1, 2.38930E-1, 7.28220E-1, 1.07770E-1, &
               1.26590E-1, 3.65690E-1, 5.88000E-1, 1.04358E-1, 1.31403E-1, 3.29224E-1, 4.72306E-1, 1.00540E-1, 1.43300E-1, &
               2.92800E-1, 3.64700E-1 & 
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_Fe_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               8.35293E-2, 1.54038E-1, 2.05848E-1, 2.67336E-1, 8.35293E-2, 1.54038E-1, 2.05848E-1, 2.67336E-1, 3.29638E-2, &
               8.65620E-2, 1.59487E-1, 2.85524E-1, 3.57983E-2, 8.96779E-2, 1.37117E-1, 2.88494E-1, 4.11340E-2, 9.43860E-2, &
               1.13910E-1, 2.93150E-1 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_Ca_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               4.62422E-3, 5.59585E-3, 1.12530E-2, 2.94888E-2, 4.62422E-3, 5.59585E-3, 1.12530E-2, 2.94888E-2, 3.11189E-3, &
               4.80689E-3, 1.68592E-2, 2.57881E-2, 3.84255E-3, 5.20613E-3, 1.31612E-2, 2.02002E-2, 5.08320E-3, 5.79037E-3, &
               9.55079E-3, 1.47510E-2 & 
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_Al_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               6.32000E-3, 5.17000E-3, 8.70000E-3, 1.96000E-2, 6.32000E-3, 5.17000E-3, 8.70000E-3, 1.96000E-2, 3.64000E-3, &
               1.04000E-2, 1.18000E-2, 2.29000E-2, 5.26100E-3, 1.47500E-2, 1.59300E-2, 3.75100E-2, 8.56000E-3, 2.34000E-2, &
               2.37000E-2, 7.20000E-2 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_H_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               8.43000E+0, 9.80000E+0, 1.11000E+1, 1.29000E+1, 8.43000E+0, 9.80000E+0, 1.11000E+1, 1.29000E+1, 8.95000E+0, &
               1.02000E+1, 1.01000E+1, 1.03000E+1, 8.49600E+0, 9.38700E+0, 9.49500E+0, 6.51000E+0, 7.93000E+0, 8.41000E+0, &
               8.75000E+0, 3.55000E+0 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))
   real(dp), dimension(1:nmet_N13_HN, 1:nmass_N13_HN) :: table_yield_He_SNII_N13_HN = transpose(reshape( &
                                                                                                         (/ &
               5.96016E+0, 7.00013E+0, 8.43014E+0, 1.08001E+1, 5.96016E+0, 7.00013E+0, 8.43014E+0, 1.08001E+1, 7.03017E+0, &
               8.49019E+0, 7.93018E+0, 8.12018E+0, 6.91220E+0, 7.93220E+0, 8.11720E+0, 6.46310E+0, 6.76024E+0, 7.25022E+0, &
               8.37021E+0, 4.78005E+0 &
               /), [nmass_N13_HN, nmet_N13_HN] &
               ))

   ! Harley edits for pop 3 stars
   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Fe_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             5.49764E-1, 9.26334E-3, 3.75658E+0, 7.36088E+0, 1.64338E+1, 1.96213E+1 &
                                             !Heger 2002 yields
                                             !3.37D-15 + 1.31D-13 + 2.80D-14 + 5.78D-14, &
                                             !2.46D-2 + 1.17D-2 + 2.40D-4 + 1.08D-8, &
                                             !4.49D-2 + 1.02D-1 + 1.08D-3 + 1.12D-8, &
                                             !4.06D-2 + 1.29D-1 + 1.21D-3 + 9.34D-9, &
                                             !4.94D-2 + 4.07D-1 + 2.63D-3 + 8.04D-9, &
                                             !6.44D-2 + 1.31D0 + 6.37D-3 + 7.18D-9, &
                                             !8.79D-2 + 2.98D0 + 1.33D-2 + 6.47D-9, &
                                             !1.30D-1 + 5.82D0 + 2.77D-2 + 6.05D-9, &
                                             !1.42D-1 + 9.55D0 + 6.97D-2 + 5.35D-9, &
                                             !1.41D-1 + 1.42D1 + 1.41D-1 + 4.41D-9, &
                                             !1.38D-1 + 1.90D1 + 2.26D-1 + 3.72D-9, &
                                             !1.31D-1 + 2.46D1 + 3.39D-1 + 2.94D-9, &
                                             !1.24D-1 + 3.17D1 + 5.04D-1 + 2.32D-9, &
                                             !1.18D-1 + 3.96D1 + 7.16D-1 + 1.88D-9 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_O_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             2.57100E+1, 5.30600E+1, 4.42400E+1, 5.59500E+1, 4.43254E+1, 4.02408E+1 &
                                             !Heger 2002 yields
                                             !4.92D1 + 4.27D-6 + 2.81D-6, &
                                             !4.58D1 + 9.85D-7 + 9.94D-7, &
                                             !4.44D1 + 7.76D-7 + 8.50D-7, &
                                             !4.68D1 + 7.47D-7 + 8.15D-7, &
                                             !4.66D1 + 6.60D-7 + 7.94D-7, &
                                             !4.59D1 + 5.93D-7 + 7.84D-7, &
                                             !4.52D1 + 5.38D-7 + 7.63D-7, &
                                             !4.39D1 + 5.09D-7 + 8.10D-7, &
                                             !4.27D1 + 4.73D-7 + 7.99D-7, &
                                             !4.11D1 + 4.45D-7 + 8.23D-7, &
                                             !3.98D1 + 4.13D-7 + 8.35D-7, &
                                             !3.81D1 + 3.74D-7 + 8.50D-7, &
                                             !3.59D1 + 3.26D-7 + 8.89D-7, &
                                             !3.34D1 + 2.86D-7 + 8.47D-7 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_N_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             2.33800E-4, 1.42406E-2, 1.03948E-2, 5.85014E-4, 3.10700E-2, 4.17700E-2 &
                                             !Heger 2002 yields
                                             !7.84D-5 + 7.27D-7, &
                                             !5.04D-5 + 7.03D-6, &
                                             !4.31D-5 + 6.69D-6, &
                                             !4.15D-5 + 6.69D-6, &
                                             !3.78D-5 + 6.57D-6, &
                                             !3.46D-5 + 6.45D-6, &
                                             !3.23D-5 + 6.41D-6, &
                                             !3.10D-5 + 6.43D-6, &
                                             !2.92D-5 + 6.51D-6, &
                                             !2.75D-5 + 6.76D-6, &
                                             !2.42D-5 + 6.99D-6, &
                                             !1.95D-5 + 7.19D-6, &
                                             !1.41D-5 + 7.03D-6, &
                                             !1.10D-5 + 6.65D-6 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Si_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             3.11211E+0, 4.79055E+0, 1.61821E+1, 2.12351E+1, 2.70521E+1, 2.90926E+1 &
                                             !Heger 2002 yields
                                             !3.15D-1 + 1.96D-3 + 6.52D-3, &
                                             !7.97D0 + 2.32D-2 + 2.74D-3, &
                                             !1.22D1 + 2.25D-2 + 1.47D-3, &
                                             !1.36D1 + 2.22D-2 + 1.33D-3, &
                                             !1.65D1 + 2.14D-2 + 1.23D-3, &
                                             !1.92D1 + 2.07D-2 + 1.23D-3, &
                                             !2.14D1 + 2.05D-2 + 1.36D-3, &
                                             !2.31D1 + 1.94D-2 + 1.20D-3, &
                                             !2.45D1 + 1.94D-2 + 1.45D-3, &
                                             !2.52D1 + 1.81D-2 + 1.35D-3, &
                                             !2.57D1 + 1.78D-2 + 1.43D-3, &
                                             !2.57D1 + 1.68D-2 + 1.37D-3, &
                                             !2.51D1 + 1.52D-2 + 1.17D-3, &
                                             !2.42D1 + 1.49D-2 + 1.35D-3 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Ca_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             1.66699E-1, 1.88414E-1, 1.31923E+0, 2.31647E+0, 2.76364E+0, 2.92056E+0 &
                                             !Heger 2002 yields
                                             !1.46D-8 + 2.39D-10 + 3.56D-11 + 6.96D-11, &
                                             !1.88D-1 + 1.69D-4 + 5.74D-8 + 3.82D-5, &
                                             !3.17D-1 + 1.39D-4 + 4.23D-8 + 6.67D-5, &
                                             !3.70D-1 + 8.72D-5 + 2.73D-8 + 7.81D-5, &
                                             !5.29D-1 + 5.74D-5 + 1.94D-8 + 1.22D-4, &
                                             !7.61D-1 + 3.67D-5 + 1.37D-8 + 1.99D-4, &
                                             !9.93D-1 + 2.35D-5 + 1.00D-8 + 2.85D-4, &
                                             !1.22D0 + 1.50D-5 + 7.48D-9 + 3.76D-4, &
                                             !1.40D0 + 9.87D-6 + 1.29D-8 + 4.58D-4, &
                                             !1.54D0 + 6.42D-6 + 8.33D-8 + 5.65D-4, &
                                             !1.63D0 + 4.53D-6 + 2.55D-7 + 7.05D-4, &
                                             !1.69D0 + 3.39D-6 + 5.85D-7 + 9.06D-4, &
                                             !1.71D0 + 3.06D-6 + 1.20D-6 + 1.21D-3, &
                                             !1.67D0 + 3.79D-6 + 2.19D-6 + 1.61D-3 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Al_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             2.85000E-2, 2.80700E-2, 2.03800E-2, 1.55700E-2, 8.62100E-2, 1.11000E-1 &
                                             !Heger 2002 yields
                                             !3.37D-2, &
                                             !1.77D-2, &
                                             !1.66D-2, &
                                             !1.59D-2, &
                                             !1.67D-2, &
                                             !1.72D-2, &
                                             !1.77D-2, &
                                             !1.63D-2, &
                                             !1.66D-2, &
                                             !1.55D-2, &
                                             !1.50D-2, &
                                             !1.42D-2, &
                                             !1.29D-2, &
                                             !1.34D-2 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Mg_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             1.91537E+0, 2.76559E+0, 1.94464E+0, 3.07667E+0, 4.80174E+0, 5.40749E+0 &
                                             !Heger 2002 yields
                                             !1.53D0 + 3.78D-3 + 3.71D-3, &
                                             !3.02D0 + 1.08D-3 + 1.74D-3, &
                                             !3.49D0 + 3.07D-3 + 1.69D-3, &
                                             !3.67D0 + 2.23D-3 + 1.73D-3, &
                                             !3.97D0 + 3.36D-3 + 1.79D-3, &
                                             !4.24D0 + 4.15D-3 + 1.92D-3, &
                                             !4.38D0 + 4.70D-3 + 2.06D-3, &
                                             !4.41D0 + 4.27D-3 + 1.92D-3, &
                                             !4.40D0 + 4.57D-3 + 2.05D-3, &
                                             !4.31D0 + 4.55D-3 + 1.91D-3, &
                                             !4.50D0 + 4.04D-3 + 1.94D-3, &
                                             !4.55D0 + 3.76D-3 + 1.86D-3, &
                                             !4.42D0 + 3.39D-3 + 1.72D-3, &
                                             !4.38D0 + 3.56D-3 + 1.73D-3 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Eu_POP3_HMHN = &
                                             (/ &
                                             1.0D-40, &
                                             1.0D-40, &
                                             1.0D-40, &
                                             1.0D-40, &
                                             1.0D-40, &
                                             1.0D-40 &
                                             !1.0D-40, &
                                             !1.0D-40, &
                                             !1.0D-40, &
                                             !1.0D-40, &
                                             !1.0D-40, &
                                             !1.0D-40, &
                                             !1.0D-40, &
                                             !1.0D-40 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_C_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             1.19100E+0, 6.32500E+0, 2.29700E+0, 4.24100E+0, 1.88647E+0, 1.05964E+0 &
                                             !Heger 2002 yields
                                             !6.89D0 + 2.33D-7, &
                                             !4.54D0 + 4.19D-7, &
                                             !4.32D0 + 7.44D-7, &
                                             !4.33D0 + 8.86D-7, &
                                             !4.28D0 + 8.48D-7, &
                                             !4.21D0 + 5.78D-7, &
                                             !4.13D0 + 3.59D-7, &
                                             !4.01D0 + 4.31D-7, &
                                             !3.85D0 + 2.30D-7, &
                                             !3.74D0 + 1.60D-7, &
                                             !3.73D0 + 1.64D-7, &
                                             !3.71D0 + 1.51D-7, &
                                             !3.61D0 + 1.38D-7, &
                                             !3.49D0 + 9.36D-8 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_Ne_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             2.24289E+0, 3.38577E+0, 1.19118E+0, 3.74843E+0, 4.69870E+0, 5.03215E+0 &
                                             !Heger 2002 yields
                                             !4.99D0 + 5.04D-4 + 1.86D-4, &
                                             !4.04D0 + 6.53D-5 + 2.01D-5, &
                                             !3.89D0 + 2.10D-4 + 1.85D-5, &
                                             !4.22D0 + 2.49D-4 + 1.97D-5, &
                                             !4.20D0 + 2.99D-4 + 2.19D-5, &
                                             !4.10D0 + 2.51D-4 + 2.27D-5, &
                                             !3.98D0 + 1.82D-4 + 2.26D-5, &
                                             !4.06D0 + 2.14D-4 + 2.37D-5, &
                                             !3.86D0 + 1.42D-4 + 2.28D-5, &
                                             !3.90D0 + 1.65D-4 + 2.46D-5, &
                                             !3.84D0 + 1.45D-4 + 2.39D-5, &
                                             !3.88D0 + 1.52D-4 + 2.41D-5, &
                                             !3.92D0 + 1.68D-4 + 2.32D-5, &
                                             !3.73D0 + 1.24D-4 + 2.02D-5 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_S_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             1.44724E+0, 1.81716E+0, 8.67592E+0, 1.31062E+1, 1.58410E+1, 1.68003E+1 &
                                             !Heger 2002 yields
                                             !9.50D-4 + 3.84D-5 + 2.06D-4 + 1.12D-8, &
                                             !2.43D0 + 3.07D-3 + 6.36D-3 + 1.12D-8, &
                                             !4.12D0 + 3.50D-3 + 4.57D-3 + 7.86D-9, &
                                             !4.70D0 + 3.64D-3 + 2.94D-3 + 6.60D-9, &
                                             !6.08D0 + 3.73D-3 + 1.93D-3 + 5.53D-9, &
                                             !7.58D0 + 3.75D-3 + 1.22D-3 + 4.55D-9, &
                                             !8.85D0 + 3.74D-3 + 7.83D-4 + 3.82D-9, &
                                             !9.97D0 + 3.68D-3 + 4.98D-4 + 3.24D-9, &
                                             !1.08D1 + 3.59D-3 + 3.28D-4 + 2.73D-9, &
                                             !1.14D1 + 3.47D-3 + 2.12D-4 + 2.36D-9, &
                                             !1.18D1 + 3.34D-3 + 1.49D-4 + 2.05D-9, &
                                             !1.20D1 + 3.18D-3 + 1.08D-4 + 1.79D-9, &
                                             !1.18D1 + 2.96D-3 + 8.83D-5 + 1.59D-9, &
                                             !1.14D1 + 2.71D-3 + 9.78D-5 + 1.49D-9 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_H_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             3.15000E+1, 3.56500E+1, 4.08100E+1, 3.72500E+1, 6.58200E+1, 7.58500E+1 &
                                             /)

   real(dp), dimension(1:nmass_POP3_HMHN) :: table_yield_He_POP3_HMHN = &
                                             (/ &
                                             ! Nomoto yields
                                             3.77700E+1, 4.17100E+1, 4.79501E+1, 4.90300E+1, 7.90100E+1, 8.95300E+1 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Fe_POP3_HN = &
                                           (/ &
                                           7.18D-04 + 8.24D-02 + 1.78D-03 + 1.28D-09, &
                                           1.79D-03 + 9.60D-02 + 1.63D-03 + 2.60D-09, &
                                           1.74D-03 + 1.59D-01 + 3.10D-03 + 9.19D-11, &
                                           3.30D-03 + 2.56D-01 + 4.24D-03 + 3.31D-09 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_O_POP3_HN = &
                                           (/ &
                                           2.03D+00 + 7.13D-08 + 2.33D-08, &
                                           2.38D+00 + 1.49D-06 + 3.87D-07, &
                                           3.92D+00 + 3.81D-08 + 5.03D-07, &
                                           6.32D+00 + 1.23D-08 + 2.93D-07 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_N_POP3_HN = &
                                           (/ &
                                           5.42D-05 + 2.95D-08, &
                                           5.96D-04 + 1.75D-07, &
                                           4.18D-05 + 2.20D-07, &
                                           3.39D-06 + 6.54D-07 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Si_POP3_HN = &
                                           (/ &
                                           1.03D-01 + 2.95D-04 + 1.13D-04, &
                                           2.31D-01 + 5.38D-04 + 6.35D-05, &
                                           2.47D-01 + 8.85D-04 + 1.47D-04, &
                                           7.20D-01 + 3.72D-03 + 2.82D-03 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Ca_POP3_HN = &
                                           (/ &
                                           4.77D-03 + 3.41D-06 + 2.60D-07 + 1.26D-04 + 1.15D-11 + 3.74D-15, &
                                           1.02D-02 + 5.15D-06 + 7.25D-08 + 7.04D-05 + 1.47D-11 + 1.17D-11, &
                                           8.22D-03 + 1.62D-06 + 1.71D-07 + 1.82D-04 + 9.25D-12 + 8.16D-13, &
                                           2.86D-02 + 1.97D-05 + 1.39D-07 + 1.75D-04 + 3.69D-11 + 1.27D-11 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Al_POP3_HN = &
                                           (/ &
                                           1.02D-06 + 1.50D-03, &
                                           1.28D-06 + 8.93D-04, &
                                           2.92D-06 + 1.55D-03, &
                                           3.80D-05 + 7.52D-03 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Mg_POP3_HN = &
                                           (/ &
                                           1.65D-01 + 1.07D-04 + 2.09D-04, &
                                           1.53D-01 + 4.57D-05 + 3.89D-05, &
                                           2.17D-01 + 1.45D-04 + 8.00D-05, &
                                           3.37D-01 + 5.95D-04 + 6.90D-05 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Eu_POP3_HN = &
                                           (/ &
                                           1.d-40, &
                                           1.d-40, &
                                           1.d-40, &
                                           1.d-40 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_C_POP3_HN = &
                                           (/ &
                                           1.90D-01 + 1.18D-08, &
                                           2.67D-01 + 6.94D-08, &
                                           3.16D-01 + 6.32D-08, &
                                           3.72D-01 + 8.18D-08 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_Ne_POP3_HN = &
                                           (/ &
                                           7.49D-01 + 3.58D-05 + 5.51D-05, &
                                           2.85D-01 + 1.22D-05 + 8.62D-06, &
                                           5.20D-01 + 3.51D-05 + 3.52D-05, &
                                           2.64D-01 + 1.41D-05 + 1.66D-05 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_S_POP3_HN = &
                                           (/ &
                                           4.27D-02 + 1.44D-04 + 1.84D-04 + 8.33D-10, &
                                           9.16D-02 + 2.31D-04 + 1.26D-04 + 5.78D-11, &
                                           8.49D-02 + 3.02D-04 + 2.70D-04 + 1.41D-09, &
                                           2.59D-01 + 8.45D-04 + 2.08D-03 + 5.04D-08 &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_H_POP3_HN = &
                                           (/ &
                                           8.77D+00 + 8.66D-17, &
                                           1.06D+01 + 2.06D-16, &
                                           1.17D+01 + 1.09D-14, &
                                           1.40D+01 + 1.66D-14  &
                                           /)

   real(dp), dimension(1:nmass_POP3_HN) :: table_yield_He_POP3_HN = &
                                           (/ &
                                           4.76D-05 + 5.96D+00, &
                                           2.11D-04 + 8.03D+00, & 
                                           2.06D-04 + 9.54D+00, &
                                           2.56D-05 + 1.18D+01  &
                                           /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Fe_POP3_SNII = &
                                             (/ &
                                             7.17257E-2, 7.23710E-2, 7.22878E-2, 7.85990E-2, &
                                             8.66280E-2, 1.19307E-1, 1.71571E-1 &
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_O_POP3_SNII = &
                                             (/ &
                                             4.49602E-1, 7.72706E-1, 1.37900E+0, 2.06800E+0, &
                                             2.58400E+0, 4.36600E+0, 7.34900E+0 &                                            
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_N_POP3_SNII = &
                                             (/ &
                                             1.83206E-3, 1.85807E-3, 1.89124E-4, 5.42604E-5, &
                                             5.93746E-4, 2.18284E-5, 2.63080E-6 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Si_POP3_SNII = &
                                             (/ &
                                             8.25518E-2, 7.36278E-2, 1.16987E-1, 1.01349E-1, &
                                             2.91174E-1, 2.48138E-1, 8.78806E-1 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Ca_POP3_SNII = &
                                             (/ &
                                             2.93585E-3, 4.43233E-3, 4.41412E-3, 5.57235E-3, &
                                             1.75464E-2, 1.28849E-2, 3.30527E-2 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Al_POP3_SNII = &
                                             (/ &
                                             3.77600E-3, 1.36900E-3, 3.14200E-3, 1.43500E-3, &
                                             8.50400E-4, 2.08800E-3, 1.11000E-2 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Mg_POP3_SNII = &
                                             (/ &
                                             8.64276E-2, 6.89072E-2, 1.58359E-1, 1.58136E-1, &
                                             1.36388E-1, 2.22102E-1, 4.08129E-1 &                                             
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Eu_POP3_SNII = &
                                             (/ &
                                             1.00E-40, 1.00E-40, 1.00E-40, 1.00E-40, &
                                             1.00E-40, 1.00E-40, 1.00E-40 &
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_C_POP3_SNII = &
                                             (/ &
                                             7.40701E-2, 1.71600E-1, 2.18500E-1, 2.00400E-1, &
                                             2.80100E-1, 3.26900E-1, 4.00400E-1 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_Ne_POP3_SNII = &
                                             (/ &
                                             1.52907E-2, 3.26554E-1, 4.93917E-1, 8.30902E-1, &
                                             4.09027E-1, 6.85406E-1, 2.85524E-1 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_S_POP3_SNII = &
                                             (/ &
                                             2.40392E-2, 3.23170E-2, 4.10375E-2, 4.84177E-2, &
                                             1.38827E-1, 1.00711E-1, 3.18566E-1 & 
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_H_POP3_SNII = &
                                             (/ &
                                             6.59500E+0, 7.57900E+0, 8.43100E+0, 8.77400E+0, &
                                             1.06400E+1, 1.16900E+1, 1.40300E+1 &
                                             /)

   real(dp), dimension(1:nmass_POP3_SNII) :: table_yield_He_POP3_SNII = &
                                             (/ &
                                             4.01404E+0, 4.40104E+0, 5.42203E+0, 5.94605E+0, &
                                             8.02621E+0, 9.53021E+0, 1.18500E+1 & 
                                             /)

   ! Mass fraction of each element in the Sun (Asplund 2009) ; not abundance 
   real(dp) :: table_Solar_Fe = 0.0012917570
   real(dp) :: table_Solar_O  = 0.0057326442
   real(dp) :: table_Solar_N  = 0.00069290468
   real(dp) :: table_Solar_Si = 0.00066494756
   real(dp) :: table_Solar_Ca = 6.4145074D-05
   real(dp) :: table_Solar_Al = 5.5625916D-05
   real(dp) :: table_Solar_Mg = 0.00070797884
   real(dp) :: table_Solar_Eu = 3.6810875D-10
   real(dp) :: table_Solar_C  = 0.0023647147
   real(dp) :: table_Solar_Ne = 0.0012564824
   real(dp) :: table_Solar_S  = 0.00030926749
   real(dp) :: table_Solar_H  = 0.73738783
   real(dp) :: table_Solar_He = 0.24924204

contains
    !=======================================================================
    subroutine init_LC18_yields
    !=======================================================================
    ! Initialize the Limongi + Chieffi 2018 Yields Tables
    implicit none
    character(LEN=256)::yields_hydrogen
    character(LEN=256)::yields_helium
    character(LEN=256)::yields_carbon
    character(LEN=256)::yields_nitrogen
    character(LEN=256)::yields_oxygen
    character(LEN=256)::yields_neon
    character(LEN=256)::yields_magnesium
    character(LEN=256)::yields_silicon
    character(LEN=256)::yields_sulfur
    character(LEN=256)::yields_iron
    character(LEN=256)::yields_calcium
    integer::j
    logical::ok

    yields_hydrogen  = trim(yields_lc18_dir)//'H/yields_H.dat'
    yields_helium    = trim(yields_lc18_dir)//'He/yields_He.dat'
    yields_carbon    = trim(yields_lc18_dir)//'C/yields_C.dat'
    yields_nitrogen  = trim(yields_lc18_dir)//'N/yields_N.dat'
    yields_oxygen    = trim(yields_lc18_dir)//'O/yields_O.dat'
    yields_neon      = trim(yields_lc18_dir)//'Ne/yields_Ne.dat'
    yields_magnesium = trim(yields_lc18_dir)//'Mg/yields_Mg.dat'
    yields_silicon   = trim(yields_lc18_dir)//'Si/yields_Si.dat'
    yields_sulfur    = trim(yields_lc18_dir)//'S/yields_S.dat'
    yields_iron      = trim(yields_lc18_dir)//'Fe/yields_Fe.dat'
    yields_calcium   = trim(yields_lc18_dir)//'Ca/yields_Ca.dat'

    ! Check that all files exist
    inquire(FILE=yields_hydrogen,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access H yields file '//TRIM(yields_hydrogen)
       call clean_stop
    endif

    inquire(FILE=yields_helium,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access He yields file '//TRIM(yields_helium)
       call clean_stop
    endif

    inquire(FILE=yields_carbon,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access C yields file '//TRIM(yields_carbon)
       call clean_stop
    endif

    inquire(FILE=yields_nitrogen,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access N yields file '//TRIM(yields_nitrogen)
       call clean_stop
    endif
   
    inquire(FILE=yields_oxygen,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access O yields file '//TRIM(yields_oxygen)
       call clean_stop
    endif

    inquire(FILE=yields_neon,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access Ne yields file '//TRIM(yields_neon)
       call clean_stop
    endif

    inquire(FILE=yields_magnesium,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access Mg yields file '//TRIM(yields_magnesium)
       call clean_stop
    endif

    inquire(FILE=yields_silicon,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access Si yields file '//TRIM(yields_silicon)
       call clean_stop
    endif

    inquire(FILE=yields_sulfur,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access S yields file '//TRIM(yields_sulfur)
       call clean_stop
    endif

    inquire(FILE=yields_iron,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access Fe yields file '//TRIM(yields_iron)
       call clean_stop
    endif

    inquire(FILE=yields_calcium,exist=ok)
    if(.not.ok)then
       write(*,*) 'Cannot access Ca yields file '//TRIM(yields_calcium)
       call clean_stop
    endif

    ! Read in the yields tables
    open(unit=880,file=yields_hydrogen,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_H_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_helium,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_He_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_carbon,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_C_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_nitrogen,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_N_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_oxygen,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_O_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_neon,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_Ne_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_magnesium,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_Mg_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_silicon,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_Si_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_sulfur,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_S_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_iron,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_Fe_SNII_lc18(j,:)
    end do
    close(880)

    open(unit=880,file=yields_calcium,form='formatted')
    do j = 1,nmet_lc18
       read(880,*) table_yield_Ca_SNII_lc18(j,:)
    end do
    close(880)

   end subroutine init_LC18_yields

   !---------------------------------------
   function AGB_H_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_H_yield

      AGB_H_yield = interp_yield(table_mass_AGB, table_met, table_yield_H_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_H_yield = AGB_H_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_H_yield

   !---------------------------------------
   function AGB_He_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_He_yield

      AGB_He_yield = interp_yield(table_mass_AGB, table_met, table_yield_He_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_He_yield = AGB_He_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_He_yield

   !---------------------------------------
   function AGB_Fe_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Fe_yield
      
      AGB_Fe_yield = interp_yield(table_mass_AGB, table_met, table_yield_Fe_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Fe_yield = AGB_Fe_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Fe_yield

   !---------------------------------------
   FUNCTION AGB_O_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_O_yield

      AGB_O_yield = interp_yield(table_mass_AGB, table_met, table_yield_O_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_O_yield = AGB_O_yield * (mass/table_mass_AGB(1))
   END function AGB_O_yield

   !---------------------------------------
   FUNCTION AGB_N_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_N_yield

      AGB_N_yield = interp_yield(table_mass_AGB, table_met, table_yield_N_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_N_yield = AGB_N_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_N_yield

   !---------------------------------------
   FUNCTION AGB_Si_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Si_yield

      AGB_Si_yield = interp_yield(table_mass_AGB, table_met, table_yield_Si_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Si_yield = AGB_Si_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Si_yield

   !---------------------------------------
   FUNCTION AGB_Ca_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Ca_yield

      AGB_Ca_yield = interp_yield(table_mass_AGB, table_met, table_yield_Ca_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Ca_yield = AGB_Ca_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Ca_yield

   !---------------------------------------
   FUNCTION AGB_Al_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Al_yield

      AGB_Al_yield = interp_yield(table_mass_AGB, table_met, table_yield_Al_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Al_yield = AGB_Al_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Al_yield

   !---------------------------------------
   FUNCTION AGB_Mg_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Mg_yield

      AGB_Mg_yield = interp_yield(table_mass_AGB, table_met, table_yield_Mg_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Mg_yield = AGB_Mg_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Mg_yield

   !---------------------------------------
   FUNCTION AGB_Eu_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Eu_yield

      AGB_Eu_yield = interp_yield(table_mass_AGB, table_met, table_yield_Eu_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Eu_yield = AGB_Eu_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Eu_yield

   !---------------------------------------
   FUNCTION AGB_C_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_C_yield

      AGB_C_yield = interp_yield(table_mass_AGB, table_met, table_yield_C_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_C_yield = AGB_C_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_C_yield

   !---------------------------------------
   FUNCTION AGB_Ne_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_Ne_yield

      AGB_Ne_yield = interp_yield(table_mass_AGB, table_met, table_yield_Ne_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_Ne_yield = AGB_Ne_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_Ne_yield

   !---------------------------------------
   FUNCTION AGB_S_yield(mass, met)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      real(dp) :: AGB_S_yield

      AGB_S_yield = interp_yield(table_mass_AGB, table_met, table_yield_S_AGB, mass, met)
      if (mass.lt.table_mass_AGB(1)) AGB_S_yield = AGB_S_yield * (mass/table_mass_AGB(1))
   END FUNCTION AGB_S_yield

   !---------------------------------------
   FUNCTION SNII_Fe_yield(mass, met, portinari, lc18, hypernova_)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Fe_yield
      
      if (present(hypernova_)) then 
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_Fe_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_Fe_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then
         SNII_Fe_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Fe_SNII_portinari, mass, met)
         SNII_Fe_yield = 0.5d0*SNII_Fe_yield !See appendix A3.2 of https://arxiv.org/pdf/0902.1535.pdf
      else if (lc18) then
         SNII_Fe_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_Fe_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_Fe_yield = SNII_Fe_yield * (mass/table_mass_lc18(1))
      else 
         SNII_Fe_yield = interp_yield(table_mass_SNII, table_met, table_yield_Fe_SNII, mass, met)
      end if
   END FUNCTION SNII_Fe_yield

   !---------------------------------------
   FUNCTION SNII_O_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_O_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_O_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_O_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then 
         SNII_O_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_O_SNII_portinari, mass, met)
      else if (lc18) then 
         SNII_O_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_O_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_O_yield = SNII_O_yield * (mass/table_mass_lc18(1))
      else
         SNII_O_yield = interp_yield(table_mass_SNII, table_met, table_yield_O_SNII, mass, met)
      end if
   END FUNCTION SNII_O_yield

   !---------------------------------------
   FUNCTION SNII_N_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_N_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_N_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_N_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then 
         SNII_N_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_N_SNII_portinari, mass, met)
      else if (lc18) then
         SNII_N_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_N_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_N_yield = SNII_N_yield * (mass/table_mass_lc18(1))
      else
         SNII_N_yield = interp_yield(table_mass_SNII, table_met, table_yield_N_SNII, mass, met)
      end if
   END FUNCTION SNII_N_yield

   !---------------------------------------
   FUNCTION SNII_Si_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Si_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_Si_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_Si_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then
         SNII_Si_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Si_SNII_portinari, mass, met)
      else if (lc18) then
         SNII_Si_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_Si_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_Si_yield = SNII_Si_yield * (mass/table_mass_lc18(1))
      else
         SNII_Si_yield = interp_yield(table_mass_SNII, table_met, table_yield_Si_SNII, mass, met)
      end if
   END FUNCTION SNII_Si_yield

   !---------------------------------------
   FUNCTION SNII_Ca_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Ca_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_Ca_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_Ca_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then
         SNII_Ca_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Ca_SNII_portinari, mass, met)
      else if (lc18) then
         SNII_Ca_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_Ca_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_Ca_yield = SNII_Ca_yield * (mass/table_mass_lc18(1))
      else
         SNII_Ca_yield = interp_yield(table_mass_SNII, table_met, table_yield_Ca_SNII, mass, met)
      end if
   END FUNCTION SNII_Ca_yield

   !---------------------------------------
   FUNCTION SNII_Al_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Al_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      if (hypernova) then
         SNII_Al_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_Al_SNII_N13_HN, mass, met)         
         return
      end if

      if (portinari) then
         SNII_Al_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Al_SNII_portinari, mass, met)
      else if (lc18) then
         SNII_Al_yield = 0.d0
      else 
         SNII_Al_yield = interp_yield(table_mass_SNII, table_met, table_yield_Al_SNII, mass, met)
      end if
   END FUNCTION SNII_Al_yield

   !---------------------------------------
   FUNCTION SNII_Mg_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Mg_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_Mg_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_Mg_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then
         SNII_Mg_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Mg_SNII_portinari, mass, met)
         SNII_Mg_yield = 2.0d0*SNII_Mg_yield !See appendix A3.2 of https://arxiv.org/pdf/0902.1535.pdf
      else if (lc18) then
         SNII_Mg_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_Mg_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_Mg_yield = SNII_Mg_yield * (mass/table_mass_lc18(1))
      else 
         SNII_Mg_yield = interp_yield(table_mass_SNII, table_met, table_yield_Mg_SNII, mass, met)
      end if
   END FUNCTION SNII_Mg_yield

   !---------------------------------------
   FUNCTION SNII_Eu_yield(mass, met, portinari, lc18, hypernova_)
      !--------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Eu_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      if (hypernova) then
         SNII_Eu_yield = 0.d0
         return
      end if

      if (portinari) then 
         SNII_Eu_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Eu_SNII_portinari, mass, met)
      else if (lc18) then 
         SNII_Eu_yield = 0.d0
      else
         SNII_Eu_yield = interp_yield(table_mass_SNII, table_met, table_yield_Eu_SNII, mass, met)
      end if
   END FUNCTION SNII_Eu_yield

   !---------------------------------------
   FUNCTION SNII_C_yield(mass, met, portinari, lc18, hypernova_)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_C_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_C_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_C_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then 
         SNII_C_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_C_SNII_portinari, mass, met)
         SNII_C_yield = 0.5d0*SNII_C_yield !See appendix A3.2 of https://arxiv.org/pdf/0902.1535.pdf
      else if (lc18) then 
         SNII_C_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_C_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_C_yield = SNII_C_yield * (mass/table_mass_lc18(1))
      else
         SNII_C_yield = interp_yield(table_mass_SNII, table_met, table_yield_C_SNII, mass, met)
      end if
   END FUNCTION SNII_C_yield

   !---------------------------------------
   FUNCTION SNII_S_yield(mass, met, portinari, lc18, hypernova_)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_S_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_S_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_S_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then 
         SNII_S_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_S_SNII_portinari, mass, met)
      else if (lc18) then 
         SNII_S_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_S_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_S_yield = SNII_S_yield * (mass/table_mass_lc18(1))
      else
         SNII_S_yield = 1.d-40
      end if
   END FUNCTION SNII_S_yield

   !---------------------------------------
   FUNCTION SNII_Ne_yield(mass, met, portinari, lc18, hypernova_)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_Ne_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      ! Hypernova can be used with other yields
      ! Call it first
      if (hypernova) then
         SNII_Ne_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_Ne_SNII_N13_HN, mass, met)
         return
      end if

      if (portinari) then 
         SNII_Ne_yield = interp_yield(table_mass_portinari, table_met_portinari, table_yield_Ne_SNII_portinari, mass, met)
      else if (lc18) then
         SNII_Ne_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_Ne_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_Ne_yield = SNII_Ne_yield * (mass/table_mass_lc18(1))
      else
         SNII_Ne_yield = 1.d-40
      end if
   END FUNCTION SNII_Ne_yield

   !---------------------------------------
   FUNCTION SNII_H_yield(mass, met, portinari, lc18, hypernova_)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_H_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      if (hypernova) then 
         SNII_H_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_H_SNII_N13_HN, mass, met) 
         return
      end if

      if (portinari) then 
         SNII_H_yield = 1.d-40
      else if (lc18) then
         SNII_H_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_H_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_H_yield = SNII_H_yield * (mass/table_mass_lc18(1))
      else
         SNII_H_yield = 1.d-40
      end if
   END FUNCTION SNII_H_yield

   !---------------------------------------
   FUNCTION SNII_He_yield(mass, met, portinari, lc18, hypernova_)
      !---------------------------------------
      real(dp), intent(in) :: mass, met
      logical, intent(in) :: portinari, lc18
      logical, intent(in), optional::hypernova_
      logical::hypernova
      real(dp) :: SNII_He_yield

      if (present(hypernova_)) then
         hypernova = hypernova_
      else
         hypernova = .false.
      end if

      if (hypernova) then 
         SNII_He_yield = interp_yield(table_mass_N13_HN, table_met_N13_HN, table_yield_He_SNII_N13_HN, mass, met) 
         return 
      end if

      if (portinari) then 
         SNII_He_yield = 1.d-40
      else if (lc18) then
         SNII_He_yield = interp_yield(table_mass_lc18, table_met_lc18, table_yield_He_SNII_lc18, mass, LOG10(met))
         if (mass.le.table_mass_lc18(1)) SNII_He_yield = SNII_He_yield * (mass/table_mass_lc18(1))
      else
         SNII_He_yield = 1.d-40
      end if
   END FUNCTION SNII_He_yield

   ! Harley Pop III
   !---------------------------------------
   FUNCTION POP3_SNII_Fe_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Fe_yield

      POP3_SNII_Fe_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Fe_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Fe_yield = POP3_SNII_Fe_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Fe_yield

   !---------------------------------------
   FUNCTION POP3_SNII_O_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_O_yield

      POP3_SNII_O_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_O_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_O_yield = POP3_SNII_O_yield * (mass/table_mass_POP3_SNII(1)) 
   END FUNCTION POP3_SNII_O_yield

   !---------------------------------------
   FUNCTION POP3_SNII_N_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_N_yield

      POP3_SNII_N_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_N_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_N_yield = POP3_SNII_N_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_N_yield

   !---------------------------------------
   FUNCTION POP3_SNII_Si_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Si_yield

      POP3_SNII_Si_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Si_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Si_yield = POP3_SNII_Si_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Si_yield

   !---------------------------------------
   FUNCTION POP3_SNII_Ca_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Ca_yield

      POP3_SNII_Ca_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Ca_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Ca_yield = POP3_SNII_Ca_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Ca_yield

   !---------------------------------------
   FUNCTION POP3_SNII_Al_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Al_yield

      POP3_SNII_Al_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Al_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Al_yield = POP3_SNII_Al_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Al_yield

   !---------------------------------------
   FUNCTION POP3_SNII_Mg_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Mg_yield

      POP3_SNII_Mg_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Mg_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Mg_yield = POP3_SNII_Mg_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Mg_yield

   !---------------------------------------
   FUNCTION POP3_SNII_Eu_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Eu_yield

      POP3_SNII_Eu_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Eu_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Eu_yield = POP3_SNII_Eu_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Eu_yield

   !---------------------------------------
   FUNCTION POP3_SNII_C_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_C_yield

      POP3_SNII_C_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_C_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_C_yield = POP3_SNII_C_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_C_yield

   !---------------------------------------
   FUNCTION POP3_SNII_Ne_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_Ne_yield

      POP3_SNII_Ne_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_Ne_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_Ne_yield = POP3_SNII_Ne_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_Ne_yield

   !---------------------------------------
   FUNCTION POP3_SNII_S_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_S_yield

      POP3_SNII_S_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_S_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_S_yield = POP3_SNII_S_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_S_yield

   !---------------------------------------
   FUNCTION POP3_SNII_H_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_H_yield

      POP3_SNII_H_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_H_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_H_yield = POP3_SNII_H_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_H_yield

   !---------------------------------------
   FUNCTION POP3_SNII_He_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_SNII_He_yield

      POP3_SNII_He_yield = interp_yield_1d(table_mass_POP3_SNII, table_yield_He_POP3_SNII, mass)
      if (mass.lt.table_mass_POP3_SNII(1)) POP3_SNII_He_yield = POP3_SNII_He_yield * (mass/table_mass_POP3_SNII(1))
   END FUNCTION POP3_SNII_He_yield

   !---------------------------------------
   FUNCTION POP3_HN_Fe_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Fe_yield

      POP3_HN_Fe_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Fe_POP3_HN, mass)
   END FUNCTION POP3_HN_Fe_yield

   !---------------------------------------
   FUNCTION POP3_HN_O_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_O_yield

      POP3_HN_O_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_O_POP3_HN, mass)
   END FUNCTION POP3_HN_O_yield

   !---------------------------------------
   FUNCTION POP3_HN_N_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_N_yield

      POP3_HN_N_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_N_POP3_HN, mass)
   END FUNCTION POP3_HN_N_yield

   !---------------------------------------
   FUNCTION POP3_HN_Si_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Si_yield

      POP3_HN_Si_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Si_POP3_HN, mass)
   END FUNCTION POP3_HN_Si_yield

   !---------------------------------------
   FUNCTION POP3_HN_Ca_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Ca_yield

      POP3_HN_Ca_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Ca_POP3_HN, mass)
   END FUNCTION POP3_HN_Ca_yield

   !---------------------------------------
   FUNCTION POP3_HN_Al_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Al_yield

      POP3_HN_Al_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Al_POP3_HN, mass)
   END FUNCTION POP3_HN_Al_yield

   !---------------------------------------
   FUNCTION POP3_HN_Mg_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Mg_yield

      POP3_HN_Mg_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Mg_POP3_HN, mass)
   END FUNCTION POP3_HN_Mg_yield

   !---------------------------------------
   FUNCTION POP3_HN_Eu_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Eu_yield

      POP3_HN_Eu_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Eu_POP3_HN, mass)
   END FUNCTION POP3_HN_Eu_yield

   !---------------------------------------
   FUNCTION POP3_HN_C_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_C_yield

      POP3_HN_C_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_C_POP3_HN, mass)
   END FUNCTION POP3_HN_C_yield

   !---------------------------------------
   FUNCTION POP3_HN_Ne_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_Ne_yield

      POP3_HN_Ne_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_Ne_POP3_HN, mass)
   END FUNCTION POP3_HN_Ne_yield

   !---------------------------------------
   FUNCTION POP3_HN_S_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_S_yield

      POP3_HN_S_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_S_POP3_HN, mass)
   END FUNCTION POP3_HN_S_yield

   !---------------------------------------
   FUNCTION POP3_HN_H_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_H_yield

      POP3_HN_H_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_H_POP3_HN, mass)
   END FUNCTION POP3_HN_H_yield

   !---------------------------------------
   FUNCTION POP3_HN_He_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HN_He_yield

      POP3_HN_He_yield = interp_yield_1d(table_mass_POP3_HN, table_yield_He_POP3_HN, mass)
   END FUNCTION POP3_HN_He_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Fe_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Fe_yield

      POP3_HMHN_Fe_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Fe_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Fe_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_O_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_O_yield

      POP3_HMHN_O_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_O_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_O_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_N_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_N_yield

      POP3_HMHN_N_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_N_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_N_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Si_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Si_yield

      POP3_HMHN_Si_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Si_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Si_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Ca_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Ca_yield

      POP3_HMHN_Ca_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Ca_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Ca_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Al_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Al_yield

      POP3_HMHN_Al_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Al_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Al_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Mg_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Mg_yield

      POP3_HMHN_Mg_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Mg_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Mg_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Eu_yield(mass)
      !--------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Eu_yield

      POP3_HMHN_Eu_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Eu_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Eu_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_C_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_C_yield

      POP3_HMHN_C_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_C_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_C_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_Ne_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_Ne_yield

      POP3_HMHN_Ne_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_Ne_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_Ne_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_S_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_S_yield

      POP3_HMHN_S_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_S_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_S_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_H_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_H_yield

      POP3_HMHN_H_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_H_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_H_yield

   !---------------------------------------
   FUNCTION POP3_HMHN_He_yield(mass)
      !---------------------------------------
      real(dp), intent(in) :: mass
      real(dp) :: POP3_HMHN_He_yield

      POP3_HMHN_He_yield = interp_yield_1d(table_mass_POP3_HMHN, table_yield_He_POP3_HMHN, mass)
   END FUNCTION POP3_HMHN_He_yield

   function interp_yield_1d(table_mass, table_yield, mass)
      real(dp), dimension(:), intent(in) :: table_mass
      real(dp), dimension(:), intent(in) :: table_yield
      real(dp), intent(in) :: mass
      real(dp) :: interp_yield_1d
      real(dp) :: min_yield=1.d-40

      integer :: imass, nmass

      real(dp) :: y1, y2, f1, f2, a1, a2, a3, a4
      nmass = size(table_mass)

      ! Find index right of mass.
      imass = 2
      do while ((table_mass(imass) < mass) .AND. (imass <= nmass - 1))
         imass = imass + 1
      end do

      ! ----------------
      ! Interpolate
      ! Store values around point.
      y1 = table_mass(imass - 1)
      y2 = table_mass(imass)

      f1 = LOG10(max(table_yield(imass - 1),min_yield))
      f2 = LOG10(max(table_yield(imass),min_yield))

      ! 1D interp
      a1 = 1.0 - ((table_mass(imass) - y1)/(y2 - y1))
      a2 = 1.0 - ((y2 - table_mass(imass))/(y2 - y1))

      interp_yield_1d = (f1*a1) + (f2*a2)
      interp_yield_1d = 10**interp_yield_1d

   end function interp_yield_1d

   function interp_yield(table_mass, table_met, table_yield, mass, met)
      real(dp), dimension(:), intent(in) :: table_met
      real(dp), dimension(:), intent(in) :: table_mass
      real(dp), dimension(:, :), intent(in) :: table_yield
      real(dp), intent(in) :: met, mass
      real(dp) :: interp_yield
      real(dp) :: min_yield=1.d-40

      integer :: nmet, nmass
      integer :: imet, imass

      real(dp) :: x1, x2, y1, y2, f11, f21, f22, f12, a1, a2, a3, a4
      nmet = size(table_met)
      nmass = size(table_mass)

      ! Find index right of metallicity.
      imet = 2
      do while ((table_met(imet) < met) .AND. (imet <= nmet - 1))
         imet = imet + 1
      end do

      ! Find index right of mass.
      imass = 2
      do while ((table_mass(imass) < mass) .AND. (imass <= nmass - 1))
         imass = imass + 1
      end do

      ! ----------------
      ! Interpolate
      ! Store values around point.
      x1 = table_met(imet - 1)
      x2 = table_met(imet)
      y1 = table_mass(imass - 1)
      y2 = table_mass(imass)

      f11 = LOG10(max(table_yield(imet - 1, imass - 1),min_yield))
      f21 = LOG10(max(table_yield(imet, imass - 1),min_yield))
      f22 = LOG10(max(table_yield(imet, imass),min_yield))
      f12 = LOG10(max(table_yield(imet - 1, imass),min_yield))

      ! Coefficients of linear system.
      a1 = f11*x2*y2/((x1 - x2)*(y1 - y2)) + f12*x2*y1/((x1 - x2)*(y2 - y1)) + &
        &  f21*x1*y2/((x1 - x2)*(y2 - y1)) + f22*x1*y1/((x1 - x2)*(y1 - y2))
      a2 = f11*y2/((x1 - x2)*(y2 - y1)) + f12*y1/((x1 - x2)*(y1 - y2)) + &
        &  f21*y2/((x1 - x2)*(y1 - y2)) + f22*y1/((x1 - x2)*(y2 - y1))
      a3 = f11*x2/((x1 - x2)*(y2 - y1)) + f12*x2/((x1 - x2)*(y1 - y2)) + &
        &  f21*x1/((x1 - x2)*(y1 - y2)) + f22*x1/((x1 - x2)*(y2 - y1))
      a4 = f11/((x1 - x2)*(y1 - y2)) + f12/((x1 - x2)*(y2 - y1)) + &
        &  f21/((x1 - x2)*(y2 - y1)) + f22/((x1 - x2)*(y1 - y2))

      ! Compute bilinear interpolation.
      interp_yield = a1 + a2*met + a3*mass + a4*met*mass
      interp_yield = 10**interp_yield

   end function interp_yield

   !---------------------------------------
   FUNCTION OBwind_Fe_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 2.2262d0
      real(dp) :: a1 = -1.1031d0
      real(dp) :: OBwind_Fe_yield

      OBwind_Fe_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_Fe_yield

   !---------------------------------------
   FUNCTION OBwind_O_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 1.8280d0
      real(dp) :: a1 = -0.8183d0
      real(dp) :: OBwind_O_yield

      OBwind_O_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_O_yield

   !---------------------------------------
   FUNCTION OBwind_N_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 2.0899d0
      real(dp) :: a1 = -0.8936d0
      real(dp) :: OBwind_N_yield

      OBwind_N_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_N_yield

   FUNCTION OBwind_Mg_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 2.0709d0
      real(dp) :: a1 = -1.5738d0
      real(dp) :: OBwind_Mg_yield

      OBwind_Mg_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_Mg_yield

   FUNCTION OBwind_Al_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 2.2067d0
      real(dp) :: a1 = -2.4871d0
      real(dp) :: OBwind_Al_yield

      OBwind_Al_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_Al_yield

   !---------------------------------------
   FUNCTION OBwind_Si_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 2.1027d0
      real(dp) :: a1 = -1.4978d0
      real(dp) :: OBwind_Si_yield

      OBwind_Si_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_Si_yield

   !---------------------------------------
   FUNCTION OBwind_Eu_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 2.2263d0
      real(dp) :: a1 = -7.6284d0
      real(dp) :: OBwind_Eu_yield

      OBwind_Eu_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_Eu_yield

   !---------------------------------------
   FUNCTION OBwind_C_yield(met)
      !---------------------------------------
      real(dp), intent(in) :: met
      real(dp) :: a0 = 1.9689d0
      real(dp) :: a1 = -1.1375d0
      real(dp) :: OBwind_C_yield

      OBwind_C_yield = OB_wind_fit(met, a0, a1)
   END FUNCTION OBwind_C_yield

   function OB_wind_fit(met, a0, a1)
      real(dp), intent(in) :: met, a0, a1
      real(dp) :: OB_wind_fit
      ! ---------
      ! Fit returns IMF weighted mass fraction of the yield.
      ! The mass limits used for OB stars are 8-60 Msun and
      ! yields are taken as imf weighted averages with yield
      ! data from NuGrid (Ritter et al., 2018).
      real(dp) :: logZ, logy

      ! Limit metallicity to range.
      if (met > 0.02) then
         logZ = -1.69897d0
      else if (met < 1.0d-4) then
         logZ = -4.d0
      else
         logZ = log10(met)
      end if

      logy = a0*logZ + a1
      OB_wind_fit = 10**logy

   end function OB_wind_fit

end module metal_yields
