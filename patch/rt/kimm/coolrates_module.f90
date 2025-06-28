MODULE coolrates_module
  ! Module for returning cooling and thermochemistry interaction rates.
  ! The temperature dependence is tabulated in rate (nonlog) versus the
  ! log of T, because it is expensive to calculate on the fly.
  ! The rates are interpolated using cubic splines, and extrapolated in
  ! log-log space if temperature is above table boundaries.
  ! Joki Rosdahl and Andreas Bleuler, September 2015.

  use amr_parameters,only:dp
  use rt_parameters,only:nIons,isH2,isHe,ixHI,ixHII,ixHeII,ixHeIII,isH2Katz,H2clumping
  implicit none

  private   ! default
  public init_coolrates_tables, update_coolrates_tables, inp_coolrates_table &
       , compCoolrate, tbl_alphaA_HII, tbl_alphaA_HeII, tbl_alphaA_HeIII    &
       , tbl_alphaB_HII, tbl_alphaB_HeII, tbl_alphaB_HeIII, tbl_beta_HI     &
       , tbl_beta_HeI, tbl_beta_HeII, tbl_cr_ci_HI, tbl_cr_ci_HeI           &
       , tbl_cr_ci_HeII, tbl_cr_ce_HI, tbl_cr_ce_HeI, tbl_cr_ce_HeII        &
       , tbl_cr_r_HII, tbl_cr_r_HeII, tbl_cr_r_HeIII, tbl_cr_bre            &
       , tbl_cr_com, tbl_cr_die                                             &
       , Epump, comp_SH2, comp_Sd, comp_Alpha_H2, comp_Beta_H2coll          &
       , comp_Beta_H2HI, PE_efficiency, comp_Beta_H2coll_new, H2_cooling_GA08

  ! Default cooling rates table parameters
  integer,parameter     :: nbinT  = 1001
  real(dp),parameter    :: Tmin   = 1d-2
  real(dp),parameter    :: Tmax   = 1d+9
  real(dp)              :: dlogTinv ! Inverse of the bin space (in K)
  real(dp)              :: hTable, h2Table, h3Table   ! Interpol constants
  real(dp)              :: one_over_lnTen, one_over_hTable, one_over_h2Table
  real(dp)              :: three_over_h2Table, two_over_h3Table
  
  real(dp),dimension(nbinT) :: T_lookup = 0d0 ! Lookup temperature in log K

  type coolrates_table
     ! Cooling and interaction rates (log):
     real(dp),dimension(nbinT)::rates  = 0d0
     ! Temperature derivatives of those rates (drate/dlog(T)):
     real(dp),dimension(nbinT)::primes = 0d0 
  end type coolrates_table

  type(coolrates_table),save::tbl_alphaA_HII ! Case A rec. coefficients
  type(coolrates_table),save::tbl_alphaA_HeII
  type(coolrates_table),save::tbl_alphaA_HeIII
  type(coolrates_table),save::tbl_alphaB_HII ! Case B rec. coefficients
  type(coolrates_table),save::tbl_alphaB_HeII
  type(coolrates_table),save::tbl_alphaB_HeIII
  type(coolrates_table),save::tbl_beta_HI ! Collisional ionisation rates
  type(coolrates_table),save::tbl_beta_HeI
  type(coolrates_table),save::tbl_beta_HeII
  type(coolrates_table),save::tbl_cr_ci_HI ! Coll. ionisation cooling
  type(coolrates_table),save::tbl_cr_ci_HeI
  type(coolrates_table),save::tbl_cr_ci_HeII
  type(coolrates_table),save::tbl_cr_ce_HI ! Coll. excitation cooling
  type(coolrates_table),save::tbl_cr_ce_HeI
  type(coolrates_table),save::tbl_cr_ce_HeII
  type(coolrates_table),save::tbl_cr_r_HII ! Recombination cooling
  type(coolrates_table),save::tbl_cr_r_HeII
  type(coolrates_table),save::tbl_cr_r_HeIII
  type(coolrates_table),save::tbl_cr_bre ! Bremsstrahlung cooling rates
  type(coolrates_table),save::tbl_cr_com ! Compton cooling rates
  type(coolrates_table),save::tbl_cr_die ! Dielectronic cooling rates
  
CONTAINS

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
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif

  ! Initialise the table lookup temperatures -----------------------------
  do iT=1, nbinT
     T_lookup(iT) = log10(Tmin) + (dble(iT)-1d0) / (dble(nbinT)-1d0)     &
                                * (log10(Tmax)-log10(Tmin))
  end do
  dlogTinv = dble(nbinT-1)/(T_lookup(nbinT)-T_lookup(1)) ! (space)^-1
  hTable = 1d0/dlogTinv                             !
  h2Table = hTable*hTable                           ! Constants for table
  h3Table = h2Table*hTable                          ! interpolation

  one_over_lnTen = 1.d0/log(10d0)
  one_over_hTable = 1.d0/hTable
  one_over_h2Table = 1.d0/h2Table
  three_over_h2Table = 3.d0/h2Table
  two_over_h3Table = 2.d0/h3Table
  
  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
     call comp_table_rates(iT,aexp)
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

  if(myid==0) print*,'Coolrates tables initialised '
901 format (20(1pe12.3))

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
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif
  tbl_cr_com%rates  = 0d0 ; tbl_cr_com%primes = 0d0
  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
     call update_table_rates(iT, aexp)
  end do ! end TK loop
  
  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
  call mpi_distribute_coolrates_table(tbl_cr_com)
#endif

  if(myid==0) print*,'Coolrates table updated'
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
  real(dp),dimension(:),allocatable :: table_mpi_sum
  integer::ierr
!-------------------------------------------------------------------------
  allocate(table_mpi_sum(nbinT))
  call MPI_ALLREDUCE(table%rates,table_mpi_sum                           &
                  ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%rates = table_mpi_sum
  call MPI_ALLREDUCE(table%primes,table_mpi_sum                          &
                  ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%primes = table_mpi_sum
  deallocate(table_mpi_sum)

END SUBROUTINE mpi_distribute_coolrates_table
#endif

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE comp_table_rates(iT, aexp)
! Fill in index iTK in all rates tables.
!-------------------------------------------------------------------------
  use rt_parameters,only:rt_OTSA
  implicit none
  integer::iT
  real(dp)::aexp, T, Ta, T5, lambda, f, hf, laHII, laHeII, laHeIII
  real(dp),parameter::kb=1.3806d-16        ! Boltzmann constant [ergs K-1]
!-------------------------------------------------------------------------
  ! Rates are stored in non-log, while temperature derivatives (primes) 
  ! are stored in dRate/dlogT (= dRate/dT * T * ln(10)).

  ! The log-log primes are just the normal primes times T/rate,
  ! i.e. dlogL/dlogT = T/L dL/dT
  
  T = 10d0**T_lookup(iT)

  ! Case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T                                ! 2.d0 * 157807.d0 / T
  f = 1.d0+(lambda/0.522)**0.47
  tbl_alphaA_HII%rates(iT)  =  1.269d-13 * lambda**1.503 / f**1.923
  tbl_alphaA_HII%primes(iT) = ( 0.90381*(f-1.)/f - 1.503 )               &
                            * log(10d0) * tbl_alphaA_HII%rates(iT)

  ! Case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tbl_alphaA_HeII%rates(iT)  = 3.d-14 * lambda**0.654
  tbl_alphaA_HeII%primes(iT) = -0.654                                    &
                             * log(10d0) * tbl_alphaA_HeII%rates(iT)

  ! Case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda =  1263030./T
  f= 1.d0+(lambda/0.522)**0.47
  tbl_alphaA_HeIII%rates(iT)  =  2.538d-13 * lambda**1.503 / f**1.923 
  tbl_alphaA_HeIII%primes(iT) =  ( 0.90381*(f-1.)/f - 1.503 )            &
                              * log(10d0) * tbl_alphaA_HeIII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T
  f= 1.d0+(lambda/2.74)**0.407
  tbl_alphaB_HII%rates(iT)  = 2.753d-14 * lambda**1.5 / f**2.242
  tbl_alphaB_HII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )                &
                             * log(10d0) * tbl_alphaB_HII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tbl_alphaB_HeII%rates(iT)  = 1.26d-14 * lambda**0.75
  tbl_alphaB_HeII%primes(iT) =  -0.75                                    &
                             * log(10d0) * tbl_alphaB_HeII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda = 1263030./T
  f= 1.d0+(lambda/2.74)**0.407
  tbl_alphaB_HeIII%rates(iT)  = 5.506d-14 * lambda**1.5 / f**2.242
  tbl_alphaB_HeIII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )              &
                              * log(10d0) * tbl_alphaB_HeIII%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HI (Maselli&'03)-------------
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f
  tbl_beta_HI%rates(iT)  = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)
  tbl_beta_HI%primes(iT) = (hf+157809.1d0/T)                             &
                         * log(10d0) * tbl_beta_HI%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)------------
  tbl_beta_HeI%rates(iT)  = 2.38d-11 * sqrt(T) / f * exp(-285335.4d0/T)
  tbl_beta_HeI%primes(iT) = (hf+285335.4d0/T)                            &
                          * log(10d0) * tbl_beta_HeI%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)-----------
  tbl_beta_HeII%rates(iT)  = 5.68d-12 * sqrt(T) / f * exp(-631515.d0/T)
  tbl_beta_HeII%primes(iT) = (hf+631515.d0/T)                            &
                           * log(10d0) * tbl_beta_HeII%rates(iT)

  ! BEGIN COOLING RATES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f

  ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
  tbl_cr_ci_HI%rates(iT)  = 1.27d-21 * sqrt(T) / f * exp(-157809.1/T)
  tbl_cr_ci_HI%primes(iT) = (hf+157809.1/T)                              &
                          * log(10d0) * tbl_cr_ci_HI%rates(iT)

  tbl_cr_ci_HeI%rates(iT)  = 9.38d-22 * sqrt(T) / f * exp(-285335.4/T)
  tbl_cr_ci_HeI%primes(iT) = (hf+285335.4/T)                             &
                           * log(10d0) * tbl_cr_ci_HeI%rates(iT)

  tbl_cr_ci_HeII%rates(iT)  = 4.95d-22 * sqrt(T) / f * exp(-631515. /T)
  tbl_cr_ci_HeII%primes(iT) = (hf+631515.0/T)                            &
                            * log(10d0) * tbl_cr_ci_HeII%rates(iT)

  ! Collisional excitation cooling from Cen'92
  tbl_cr_ce_HI%rates(iT)  = 7.5d-19 / f * exp(-118348./T)
  tbl_cr_ce_HI%primes(iT) = (118348./T - 0.5d0 * sqrt(T5) / f )          &
                          * log(10d0) * tbl_cr_ce_HI%rates(iT)

  tbl_cr_ce_HeI%rates(iT)  = 9.10d-27 * T**(-0.1687) / f * exp(-13179./T) 
  tbl_cr_ce_HeI%primes(iT) = (13179./T - 0.1687 - 0.5d0 * sqrt(T5) / f)  &
                           * log(10d0) * tbl_cr_ce_HeI%rates(iT)

  tbl_cr_ce_HeII%rates(iT) = 5.54d-17 * T**(-0.397)  / f * exp(-473638./T)
  tbl_cr_ce_HeII%primes(iT) = (473638./T - 0.397 - 0.5d0 * sqrt(T5) / f) &
                            * log(10d0) * tbl_cr_ce_HeII%rates(iT)
  
  ! Recombination Cooling (Hui&Gnedin'97)
  laHII    = 315614./T                                             
  laHeII   = 570670./T                                            
  laHeIII  = 1263030./T
  if(.not. rt_otsa) then ! Case A
     f = 1.d0+(laHII/0.541)**0.502                         
     tbl_cr_r_HII%rates(iT)    = 1.778d-29 * laHII**1.965 / f**2.697 * T
     tbl_cr_r_HII%primes(iT)   = (-0.965 + 1.35389*(f-1.)/f)             &
                               * log(10d0) * tbl_cr_r_HII%rates(iT)
 
     tbl_cr_r_HeII%rates(iT)   = 3.d-14 * laHeII**0.654 * kb * T
     tbl_cr_r_HeII%primes(iT)  = 0.346                                   &
                               * log(10d0) * tbl_cr_r_HeII%rates(iT)

     f = 1.d0+(laHeIII/0.541)**0.502                      
     tbl_cr_r_HeIII%rates(iT) = 14.224d-29 * laHeIII**1.965 / f**2.697 * T
     tbl_cr_r_HeIII%primes(iT)= (-0.965 + 1.35389*(f-1.)/f)              &
                              * log(10d0) * tbl_cr_r_HeIII%rates(iT)
  else ! Case B
     f = 1.d0+(laHII/2.25)**0.376                        
     tbl_cr_r_HII%rates(iT)    = 3.435d-30 * laHII**1.97 / f**3.72 * T
     tbl_cr_r_HII%primes(iT)   = (-0.97 + 1.39827*(f-1.)/f)              &
                               * log(10d0) * tbl_cr_r_HII%rates(iT)

     tbl_cr_r_HeII%rates(iT)   = 1.26d-14 * laHeII**0.75 * kb * T
     tbl_cr_r_HeII%primes(iT)  = 0.25                                    &
                               * log(10d0) * tbl_cr_r_HeII%rates(iT)

     f = 1.d0+(laHeIII/2.25)**0.376                       
     tbl_cr_r_HeIII%rates(iT)  = 27.48d-30 * laHeIII**1.97 / f**3.72 * T
     tbl_cr_r_HeIII%primes(iT) = (-0.97 + 1.39827*(f-1.)/f)              &
                               * log(10d0) * tbl_cr_r_HeIII%rates(iT)
  endif

  ! Bremsstrahlung from Osterbrock & Ferland 2006
  tbl_cr_bre%rates(iT)  = 1.42d-27 * 1.5 * sqrt(T)
  tbl_cr_bre%primes(iT) = 0.5                                            &
                        * log(10d0) * tbl_cr_bre%rates(iT)

  ! Compton Cooling from Haimann et al. 96, via Maselli et al.
  ! Need to make sure this is done whenever the redshift changes!
  Ta     = 2.727/aexp                          
  tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
  tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
                         * log(10d0) * tbl_cr_com%rates(iT)

  ! Dielectronic recombination cooling, from Black 1981
  f = 1.24d-13*T**(-1.5d0)*exp(-470000.d0/T)
  tbl_cr_die%rates(iT) = f*(1.D0+0.3d0*exp(-94000.d0/T))
  tbl_cr_die%primes(iT)=0d0
  if(tbl_cr_die%rates(iT) .gt. 0d0) then ! Can simplify w algebra
     tbl_cr_die%primes(iT) = (tbl_cr_die%rates(iT)*(564000.-1.5*T)       &
                              - f*94000.) /T**2 * T * log(10d0) 
  endif

END SUBROUTINE comp_table_rates

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_table_rates(iT,aexp)
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
  Ta     = 2.727/aexp                          
  tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
  tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
                         * log(10d0) * tbl_cr_com%rates(iT)    

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
  real(dp),intent(in)::T
  real(dp),optional::retPrime
  real(dp)::inp_coolrates_table
  integer,save:: iT = 1
  real(dp),save:: facT, yy, yy2, yy3, fa, fb, fprimea, fprimeb, Tlast=-1
  real(dp),save:: alpha, beta, gamma
  logical,save::extrap
!-------------------------------------------------------------------------
  if (.not. (T .eq. Tlast)) then    ! Reuse index if same T from last call
     ! Log of T, snapped to table at the lower boundary, but allowed
     ! to go above upper boundary, in which case we use extrapolation:
     facT = MAX( log10(T), T_lookup(1) )
     extrap=.false.
     if(facT .gt. T_lookup(nbinT)) extrap=.true. ! Above upper limit
     ! Lower closest index in table:
     iT = MIN(MAX(int((facT-T_lookup(1))*dlogTinv)+1, 1), nbinT-1) 
     yy=facT-T_lookup(iT)  ! Dist., in log(T), from T to lower table index
     yy2=yy*yy             ! That distance squared
     yy3=yy2*yy            ! ...and cubed
     Tlast = T
  endif

  if(extrap) then ! TK above upper table limit, so extrapolate in log-log:
     alpha = (log10(rates_table%rates(nbinT))             &
              -log10(rates_table%rates(nbinT-1)) )        &
           / (T_lookup(nbinT)-T_lookup(nbinT-1))
     inp_coolrates_table = 10d0**(log10(rates_table%rates(nbinT))        &
                                  + alpha * (facT - T_lookup(nbinT)))
     if( present(retPrime) )                                             &
          retPrime = alpha * inp_coolrates_table / T
     return
  endif 
  
  fa = rates_table%rates(iT)           !
  fb = rates_table%rates(iT+1)         !  Values at neighbouring table
  fprimea=rates_table%primes(iT)       !  indexes
  fprimeb=rates_table%primes(iT+1)     !

  ! Spline interpolation:
  alpha = fprimea
  beta =(fb-fa) * three_over_h2Table - (2d0*fprimea+fprimeb) * one_over_hTable
  gamma = (fprimea+fprimeb) * one_over_h2Table - (fb-fa) * two_over_h3Table
  inp_coolrates_table = fa+alpha*yy+beta*yy2+gamma*yy3
  if( present(retPrime) )                                                &
       retPrime = (alpha+2d0*beta*yy+3d0*gamma*yy2) / T * one_over_lnTen
END FUNCTION inp_coolrates_table
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Alpha_H2(T,Z,xe,xHII_opt)
! Returns creation rate of h2 on dust [cm^3 s-1] (Draine and Bertoldi 1996)
! plus gas phase rate for low Z on H- assuming equilibrium abundances for H-
! as explained in the Appendix of McKee and Krumholz (2010)
! T           => Temperature [K]
! Z           => Metallicity in Solar units
!-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,Z,xe
  real(dp),intent(in),optional::xHII_opt
  real(dp)::comp_Alpha_H2,lambda,T2, clumping
  real(dp)::logT, lnTe, k1, k2, k5, k13, xHII
!-------------------------------------------------------------------------
  T2=T/1d2
  xHII=0d0
  if(present(xHII_opt)) xHII=xHII_opt
  if(isH2Katz)then
! Note that harley edited this function.  If you look in the appendix of McKee and
! Krumholz, the zero metallicity limit is also dependent on the electron fraction.
! If you don't include this, the H2 fraction is likely to go crazy at high redshift
! which is probably unphysical.  You need the electrons for the H- channel of H2 formation.
! The real zero metallicity limit should be R- = (8*10^-19)*(ne/10^-3nH)*(T/1000)^0.88

! We change the first part of this ewquation to match equation 4 from 
! http://mnras.oxfordjournals.org/content/425/4/3058.full.pdf
! which gets their results from wolfire 2008
     clumping = 1d1 !clumping factor as used in Gnedin 2009
     if(H2clumping>0) clumping = H2clumping
     comp_Alpha_H2 =  Z * (3.5d-17) * clumping * (1d0 - xHII)
  else
     comp_Alpha_H2 =  Z * 6.0d-18*(T**0.5)/ &
       & (1.0+0.4*T2**0.5+0.2*T2+0.08*T2**2)
  endif
  if(isH2Katz)then
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
     if(T<6000)then
        k1 = 10.0**(-17.845 + 0.762*logT + 0.1523*logT**2. - 0.03274*logT**3.) 
     else
        k1 = 10.0**(-16.420 + 0.1998*logT**2. -5.447d-3*logT**4. + 4.0415d-5*logT**6.)
     endif

     ! H- + H -> H2 + e 
     if(T<300)then
        k2 = 1.5d-9
     else
        k2 = 4.0d-9*T**(-0.17)
     endif
  
     ! H- + H+ -> H + H
     k5 = 2.4d-6/sqrt(T)*(1.0 + T/20000.)
    
     ! H- + e -> H + e + e
     k13 = -1.801849334d1 + 2.36085220d0*lnTe -2.82744300d-1*lnTe**2. &
        & + 1.62331664d-2*lnTe**3. - 3.36501203d-2*lnTe**4. + 1.17832978d-2*lnTe**5. &
        & - 1.65619470d-3*lnTe**6. + 1.06827520d-4*lnTe**7. - 2.63128581d-6*lnTe**8.
     k13 = exp(k13)

     comp_Alpha_H2 = comp_Alpha_H2 + k1*k2*xe/(k2 + k5*xHII + k13*xe) ! k5=k13=0 should recovers the next line
  else
     ! Zero metallicity limit
     comp_Alpha_H2 = comp_Alpha_H2 + (8.0d-19)*(xe/0.001)*((T/1000.0)**0.88)
  endif
 
END FUNCTION comp_Alpha_H2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_H2HI(T)

! Returns collisional rate [cm3 s-1] of H2 and HI (Abel 1997->Dove&Mandy 1986)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::T
  real(dp)::comp_Beta_H2HI, kbT
!-------------------------------------------------------------------------
   kbT=8.618d-5*T !eV
  comp_Beta_H2HI = &
                3.324d-9*(kbT**2.012)*exp(-4.463d0/kbT)/(1.d0+0.2472d0*kbT)**3.512
END FUNCTION comp_Beta_H2HI
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_H2H2(T)

! Returns collisional rate [cm3 s-1] of H2 and H2 (Martin&Keogh&Mandy 1998y)
! T           => Temperature [K]
!-------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::T
  real(dp)::comp_Beta_H2H2, kbT
!-------------------------------------------------------------------------
   kbT=8.618*T !eV
  comp_Beta_H2H2 = &
                2.519d-5*(kbt**4.1881)*exp(-0.1731d0/kbT)/(1.d0+2.1347d0*kbT)**5.6881
END FUNCTION comp_Beta_H2H2

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_hi(T,J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,J
  real(dp)::gamma_hi, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1.d3
  jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
  gamma_hi = & 
       jconsts*(1.0d-11*sqrt(T3)/(1.0d0+60.0d0*T3**(-4)) + 1.0d-12*T3)
    
END FUNCTION gamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_h2(T,J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,J
  real(dp)::gamma_h2, T3, jconsts
  !----------------------------------------------------------
  T3 = T/1.d3
  jconsts=(0.276*J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
  gamma_h2 = jconsts*(3.3d-12 +6.6d-12*T3)
    
END FUNCTION gamma_h2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION dgamma_hi(T,J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::T,J
  real(dp)::dgamma_hi, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1.d3
  jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
  dgamma_hi = & 
       jconsts*(1.d-11*sqrt(T3)/(1.+60.*T3**(-4))*(1./(2.*T3)+ &
         240.*T3**(-5)/(1.+60.*T3**(-4)))+1.d-12)/1.d3
    
END FUNCTION dgamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION dgamma_h2(J)
  !-------------------------------------------------------------------------
  implicit none  
  real(dp),intent(in)::J
  real(dp)::dgamma_h2, jconsts
  !-------------------------------------------------------------------------
  jconsts=0.276*(J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
  dgamma_h2 = & 
       jconsts*6.6d-12
    
END FUNCTION dgamma_h2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Sd(nHI,nH2,dx_SS,Z)
! Returns the self shielding factor for dust
! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
  implicit none
  real(dp),intent(in)::nHI,nH2,dx_SS,Z
  real(dp)::comp_Sd,Sdeff,cNHI,cNH2
  Sdeff = 4.0d-21    !dust cross section cm^2 
  cNHI = nHI*dx_SS   !HI column density
  cNH2 = nH2*dx_SS   !H2 column density
  comp_Sd = exp(-Sdeff*Z*(cNHI + (2.0*cNH2)))
END FUNCTION comp_Sd
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_SH2(nH2,dx_SS)
! Returns the self shielding factor for dust
! see section 2.2 http://iopscience.iop.org/0004-637X/697/1/55/pdf/apj_697_1_55.pdf
  implicit none
  real(dp),intent(in)::nH2,dx_SS
  real(dp)::comp_SH2,xfac,cNH2,wH2,Sa,Sb,Sc
  cNH2 = nH2*dx_SS  !H2 column density
  xfac = cNH2/(5.0d14)
  wH2 = 0.2

  Sa = (1.0-wH2)/((1.0+xfac)*(1.0+xfac))
  Sb = wH2/sqrt(1.0+xfac)
  Sc = exp(-0.00085*sqrt(1.0+xfac))

  comp_SH2 = Sa + (Sb*Sc)
END FUNCTION comp_SH2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_H2coll(T,nH,xHI,xH2,xHe,ne,nHI,nH2,nHeI)
! Returns the collisional dissociation rates of H2 for four different
! reactions [cm3s-1] from Glover & Abel (2008)
! http://mnras.oxfordjournals.org/content/388/4/1627.full.pdf
  implicit none
  real(dp),intent(in)::T,nH,xHI,xH2,xHe,ne,nHI,nH2,nHeI
  real(dp)::comp_Beta_H2coll,ncrH,ncrH2,ncrHe,T4,invncr
  real(dp):: LTEfac, NLTEfac
  real(dp):: k8, k8L, k9, k9L, k10, k10L, k11, k11L
  real(dp):: lk8,lk9,lk10,lk11

  T4 = T/1d4

  ! Critical number densities.  See eqns 15, 16, and 17
  ncrH = 10.0**(3.0 - 0.416*log10(T4) - 0.327*log10(T4)*log10(T4))
  ncrH2 = 10.0**(4.845 - 1.3*log10(T4) + 1.62*log10(T4)*log10(T4))
  ncrHe = 10.0**(5.0792*(1.0 - 1.23d-5*(T-2000.0)))

  ! 1/ncr.  see eqn 14
  invncr = (xHI/ncrH) + (xH2/ncrH2) + (xHe/ncrHe)

  ! prefactors for LTE and NLTE collision rates  see eqn 13
  LTEfac = (nH*invncr)/(1.0 + (nH*invncr))
  NLTEfac = 1.0/(1.0 + (nH*invncr))
  if(ne>nHI)then
     LTEfac=1.d0
     NLTEfac=0.d0
  endif

  !reaction rates from the appendix:
  !H2 + e- --> H + H + e-
  !k8 = (4.49d-9)*(T**0.11)*exp(-101858.0/T)
  !k8L = (1.91d-9)*(T**0.136)*exp(-53407.1/T)
  k8 = (3.73d-9)*(T**0.1121)*exp(-99430.0/T) ! Glover et al. (2010)

  !H2 + H --> H + H + H
  k9 = (6.67d-12)*sqrt(T)*exp(-1.0*(1.0 + (63593.0/T)))
  k9L = (3.52d-9)*exp(-43900.0/T)

  !H2 + H2 --> H2 + H + H
  k10 = (((5.996d-30)*(T**4.1881))/((1.0 + (6.761d-6)*T)**5.6881))*exp(-54657.4/T)
  k10L = (1.3d-9)*exp(-53300.0/T)

  !H2 + He --> H + H + He
  k11 = 10.0**(-27.029 + (3.801*log10(T)) - (29487.0/T))
  k11L = 10.0**(-2.729 - (1.75*log10(T)) - (23474.0/T))

  k8   = max(k8,  1d-40)
  k9   = max(k9,  1d-40)
  k10  = max(k10, 1d-40)
  k11  = max(k11, 1d-40)
  !k8L  = max(k8L, 1d-40)
  k9L  = max(k9L, 1d-40)
  k10L = max(k10L,1d-40)
  k11L = max(k11L,1d-40)

  !Log of all the rates
  !lk8  = (LTEfac*log10(k8L)) + (NLTEfac*log10(k8))
  lk8  = log10(k8)
  lk9  = (LTEfac*log10(k9L)) + (NLTEfac*log10(k9))
  lk10 = (LTEfac*log10(k10L)) + (NLTEfac*log10(k10))
  lk11 = (LTEfac*log10(k11L)) + (NLTEfac*log10(k11))

  comp_Beta_H2coll = (ne*(10.0**lk8)) + (nHI*(10.0**lk9)) + (nH2*(10.0**lk10)) + (nHeI*(10.0**lk11))
  comp_Beta_H2coll = MAX(comp_Beta_H2coll/nH, 1.d-40) ! Taysun
END FUNCTION comp_Beta_H2coll
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION comp_Beta_H2coll_new(T,nH,ne,nHI,nH2,nHeI)
! Returns the collisional dissociation rates of H2 for four different
! reactions [cm3s-1] from Glover et al (2010) (see table B1)
! http://mnras.oxfordjournals.org/content/404/1/2.full.pdf
 implicit none
 real(dp),intent(in)::T,nH,ne,nHI,nH2,nHeI
 real(dp)::comp_Beta_H2coll_new,ncrHe
 real(dp):: k8, k9, k10, k11, k11L, k11H

 !H2 + e- --> H + H + e-
 k8 = (3.73d-9)*(T**0.1121)*exp(-99430.0/T)

 !H2 + H --> H + H + H
 k9 = (6.67d-12)*sqrt(T)*exp(-1.0*(1.0 + (63593.0/T)))

 !H2 + H2 --> H2 + H + H
 k10 = (((5.996d-30)*(T**4.1881))/((1.0 + (6.761e-6)*T)**5.6881))*exp(-54657.4/T)

 !H2 + He --> H + H + He
 k11L = 10.0**(-27.029 + (3.801*log10(T)) - (29487.0/T))
 k11H = 10.0**(-2.729 - (1.75*log10(T)) - (23474.0/T))
 ncrHe = 10.0**(5.0792*(1.0 - ((1.23d-5)*(T-2000.0))))
 k11 = 10.0**(log10(k11H) - ((log10(k11H) - log10(k11L))/(1.0 + (nHeI/ncrHe)**1.09)))

 comp_Beta_H2coll_new = (ne*(k8)) + (nHI*(k9)) + (nH2*(k10)) + (nHeI*(k11))

 comp_Beta_H2coll_new = MAX(comp_Beta_H2coll_new/nH, 1.d-40) ! Taysun
END FUNCTION comp_Beta_H2coll_new
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION Epump(nH,T,xH2,xHI)

! Energy converted to heat from UV pumping
! see Appendix A http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1990ApJ...365..620B&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
  implicit none
  real(dp),intent(in)::nH,T,xH2,xHI
  real(dp)::Epump,Crad,Cdex,Cfrac,ev2erg
  ev2erg = 1.602d-12
  Crad = 2.0d-7     !radiation de-excitation rate Burton 1990
  Cdex = (1.0d-12)*((1.4*xH2*exp(-18100d0/(T+1200d0))) + (1.0*xHI*exp(-1000d0/T)))*sqrt(T)*nH  !collisional de-excitation rate Burton 1990
  Cfrac = Cdex / (Cdex + Crad)
  Epump = 2 * Cfrac * ev2erg  !ergs
END FUNCTION Epump
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION PE_efficiency(G_0,Tk,n_e)

! Photoelectric heating efficiency from Bakes & Tielens (1994), Baczynski et al. (2015)
! G_0: the strength of the local ISRF normalised to the integrated Habing field
! Tk : temperature of a cell in K
! n_e: electron density in units of cm-3
  implicit none
  real(dp),intent(in)::G_0,Tk,n_e
  real(dp)::phi_pah,fact,PE_efficiency
  phi_pah = 0.5 ! wolfire+(03)
  fact    = G_0*sqrt(Tk)/(n_e*phi_pah)
  PE_efficiency  = 4.9d-2/(1d0 + 4d-3*fact**0.73) &
    &            + 3.7d-2*(Tk/1d4)**0.7/(1d0 + 2d-4*fact)
END FUNCTION PE_efficiency
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION H2_cooling_GA08(T, n_e, n_HI, n_HII, n_HeI, n_H2)
! Cooling due to molecular hydrogen
  implicit none
  real(dp),intent(in)::T,n_e,n_HI,n_HII,n_HeI,n_H2
  real(dp)::t3,lt,lt3,ltt
  real(dp)::gphdl,HDLR,HDLV,gaHI,gaH2,gaHp,gaHe,gael,galdl
  real(dp)::H2_cooling_GA08

  H2_cooling_GA08 = 0d0
  if (T .le. 10) return 
 
  lt  = log10(T)
  t3  = T/1000.
  lt3 = log10(t3)

  ! Galli & Palla (1998) ; low-density limit
  !gpldl = 10.d0**(-103. + 97.59*lt - 48.05*lt**2d0 + 10.80*lt**3d0 - 0.9032*lt**4d0) !erg cm^3 /s

  ! Hollenbach & McKee (1979) ; high-density limit (LTE)
  HDLR = ((9.5d-22*t3**3.76d0) / (1.d0+0.12d0*t3**2.1d0) * exp(-(0.13d0/t3)**3)+3.d-24*exp(-0.51d0/t3))
  HDLV = (6.7d-19*exp(-5.86d0/t3) + 1.6d-18 * exp(-11.7d0/t3))
  gphdl = (HDLR + HDLV) ! erg/s

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
     &         + 37.383713d0 * lt3      &
     &         + 58.145166d0 * lt3**2   &
     &         + 48.656103d0 * lt3**3   &
     &         + 20.159831d0 * lt3**4   &
     &         + 3.847961d0 * lt3**5)
  else if (T .lt. 1000) then
     gaHI = 10.d0**(-24.311209d0        &
     &         + 3.5692468d0 * lt3      &
     &         - 11.33286d0 * lt3**2    &
     &         - 27.850082d0 * lt3**3   &
     &         - 21.328264d0 * lt3**4   &
     &         - 4.2519023d0 * lt3**5)
  else if (T .lt. 6000) then
     gaHI = 10.d0**(-24.311209d0        &
     &         + 4.6450521d0 * lt3       &
     &         - 3.7209846d0 * lt3**2    &
     &         + 5.9369081d0 * lt3**3    &
     &         - 5.5108047d0 * lt3**4    &
     &         + 1.5538288d0 * lt3**5)
  endif

  ! Excitation by H2
  if (T .lt. 6000)then
     gaH2 = 10.d0**(-23.962112d0         &
     &        + 2.0943374d0  * lt3       &
     &        - 0.77151436d0  * lt3**2   &
     &        + 0.43693353d0  * lt3**3   &
     &        - 0.14913216d0  * lt3**4   &
     &        - 0.033638326d0 * lt3**5)
  endif

  ! Excitation by He
  if (T .gt. 100 .and. T .lt. 6000) then
     gaHe = 10.d0**(-23.689237d0        &
     &        + 2.1892372d0  * lt3      &
     &        - 0.81520438d0 * lt3**2   &
     &        + 0.29036281d0 * lt3**3   &
     &        - 0.16596184d0 * lt3**4   &
     &        + 0.19191375d0 * lt3**5)
  endif

  ! Excitation by H+
  if (T .lt. 10000) then
     gaHp = 10.d0**(-21.716699d0         &
     &        + 1.3865783d0   * lt3      &
     &        - 0.37915285d0  * lt3**2   &
     &        + 0.11453688d0  * lt3**3   &
     &        - 0.23214154d0  * lt3**4   &
     &        + 0.058538864d0 * lt3**5)
  endif

  ! Excitation by electrons
  if (T .lt. 200) then
     gael = 10.d0**(-34.286155d0          &
     &          - 48.537163d0  * lt3      &
     &          - 77.121176d0  * lt3**2   &
     &          - 51.352459d0  * lt3**3   &
     &          - 15.16916d0  * lt3**4    &
     &          - 0.98120322d0 * lt3**5)
  else if (T .lt. 10000) then
     gael = 10.d0**(-22.190316           &
     &          + 1.5728955d0  * lt3      &
     &          - 0.213351d0 * lt3**2     &
     &          + 0.96149759d0 * lt3**3   &
     &          - 0.91023195d0 * lt3**4   &
     &          + 0.13749749d0 * lt3**5)
  endif

  galdl = gaHI*n_HI + gaH2*n_H2 + gaHe*n_HeI + gaHp*n_HII + gael * n_e ! erg/s
  
  H2_cooling_GA08 = n_H2*gphdl/(1.d0 + gphdl/galdl) ! erg/cm^3/s 
 
END FUNCTION H2_cooling_GA08
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
  real(dp),dimension(nIons)::nN, nI
  logical::RT_OTSA!-------------------------------------------------------
  real(dp),save::nHI,nHII,nH2,nHeI,nHeII,nHeIII
  real(dp),save::ci_HI,ci_HeI,ci_HeII
  real(dp),save::ci_HI_prime,ci_HeI_prime,ci_HeII_prime
  real(dp),save::ce_HI,ce_HeI,ce_HeII
  real(dp),save::ce_HI_prime,ce_HeI_prime,ce_HeII_prime
  real(dp),save::r_HII,r_HeII,r_HeIII
  real(dp),save::r_HII_prime,r_HeII_prime,r_HeIII_prime
  real(dp),save::bre, bre_prime, com, com_prime, die, die_prime
  real(dp),save::ne_nHI, ne_nHII, ne_nHeI, ne_nHeII, ne_nHeIII
  real(dp),save::lowrleft,lowrright,lowr_hi,lowr_h2
  real(dp),save::lowvleft_hi,lowvright_hi,lowv_hi,lowv_h2
  real(dp),save::lowtot_hi,lowtot_h2,coolH2,dcoolH2,TT
!-------------------------------------------------------------------------
  nHI       = nN(ixHII)
  nHII      = nI(ixHII)
  if(isH2) nH2 = nN(ixHI)
  if(isHe)then
     nHeI   = nN(ixHeII)
     nHeII  = nI(ixHeII)
     nHeIII = nI(ixHeIII)
  endif

  ne_nHI    = ne * nHI
  ne_nHII   = ne * nHII
  if(isHe)then
     ne_nHeI   = ne * nHeI
     ne_nHeII  = ne * nHeII
     ne_nHeIII = ne * nHeIII
  endif

  ! Coll. Ionization Cooling
  ci_HI   = inp_coolrates_table(tbl_cr_ci_HI, T, ci_HI_prime)            &
          * ne_nHI
  ci_HI_prime    = ci_HI_prime   * ne_nHI
  if(isHe)then
     ci_HeI  = inp_coolrates_table(tbl_cr_ci_HeI, T, ci_HeI_prime)          &
             * ne_nHeI
     ci_HeII = inp_coolrates_table(tbl_cr_ci_HeII, T, ci_HeII_prime)        &
             * ne_nHeII
     ci_HeI_prime   = ci_HeI_prime  * ne_nHeI
     ci_HeII_prime  = ci_HeII_prime * ne_nHeII
  endif

  ! Collisional excitation cooling
  ce_HI   = inp_coolrates_table(tbl_cr_ce_HI, T, ce_HI_prime)            &
          * ne_nHI
  ce_HI_prime    = ce_HI_prime   * ne_nHI
  if(isHe)then
     ce_HeI  = inp_coolrates_table(tbl_cr_ce_HeI, T, ce_HeI_prime)          &
             * ne_nHeI
     ce_HeII = inp_coolrates_table(tbl_cr_ce_HeII, T, ce_HeII_prime)        &
             * ne_nHeII
     ce_HeI_prime   = ce_HeI_prime  * ne_nHeI
     ce_HeII_prime  = ce_HeII_prime * ne_nHeII
  endif

  ! Recombination Cooling
  r_HII   = inp_coolrates_table(tbl_cr_r_HII, T, r_HII_prime)            &
          * ne_nHII
  r_HII_prime   = r_HII_prime   * ne_nHII
  if(isHe) then
     r_HeII  = inp_coolrates_table(tbl_cr_r_HeII, T, r_HeII_prime)          &
             * ne_nHeII
     r_HeIII = inp_coolrates_table(tbl_cr_r_HeIII, T, r_HeIII_prime)        &
             * ne_nHeIII
     r_HeII_prime  = r_HeII_prime  * ne_nHeII
     r_HeIII_prime = r_HeIII_prime * ne_nHeIII
  endif

  ! Bremsstrahlung
  bre  = inp_coolrates_table(tbl_cr_bre, T, bre_prime)                   &
       * ne * (nHII + nHeII + 4. * nHeIII)
  bre_prime = bre_prime  * ne * (nHII + nHeII + 4. * nHeIII)

  ! Compton Cooling
  com       = inp_coolrates_table(tbl_cr_com, T, com_prime) * ne    
  com_prime = com_prime * ne

  ! Dielectronic recombination cooling
  if(isHe) then
     die = inp_coolrates_table(tbl_cr_die, T, die_prime) * ne_nHeII
     die_prime = die_prime * ne_nHeII
  endif

  ! H2 cooling Hallenbacn McKee (1979) + Halle Combes (2012) in cgs
  if(isH2) then
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

     coolH2 = H2_cooling_GA08(T, ne, nHI, nHII, nHeI, nH2)
     TT=T*1.001
     dcoolH2 = H2_cooling_GA08(TT, ne, nHI, nHII, nHeI, nH2)
   
     dcoolH2 = (dcoolH2 - coolH2) / (0.001*T)

  endif
 
  ! Overall Cooling
  compCoolrate  = ci_HI + r_HII + ce_HI + com + bre
  if(isHe) compCoolrate = compCoolrate &
                 + ci_HeI +  r_HeII +  ce_HeI + die  &
                + ci_HeII + r_HeIII + ce_HeII
  if(isH2) compCoolrate = compCoolrate + coolH2

  dCooldT       = ci_HI_prime   + r_HII_prime   + ce_HI_prime  &
                + bre_prime     + com_prime
  if(isHe) dCooldT = dCooldT &
                + ci_HeI_prime  + r_HeII_prime  + ce_HeI_prime           &
                + ci_HeII_prime + r_HeIII_prime + ce_HeII_prime          &
                + die_prime
  if(isH2) dCooldT = dCooldT + dcoolH2

END FUNCTION compCoolrate


END MODULE coolrates_module
