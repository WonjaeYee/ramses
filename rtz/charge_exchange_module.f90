! charge_exchange_module.f90
module charge_exchange_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: load_ct_rates, charge_transfer_recombination, charge_transfer_ionization, ct_tables_prebuild

  real(dp):: CTRecomb(6,4,30)
  real(dp):: CTIon(7,4,30)
  ! the rate of each (stage, element) at the last temperature it was asked for: the chemistry and
  ! the cooling ask for the same pairs at the same T several times per cooling step. One cache per
  ! MPI rank: shared module state, so it would need to become thread-private under OpenMP
  real(dp), save :: ctr_T(4,30) = -1d0, ctr_rate(4,30) = 0d0
  real(dp), save :: cti_T(4,30) = -1d0, cti_rate(4,30) = 0d0
  ! log(T/1e4) of the last T, shared by all the pairs asked for at that T
  real(dp), save :: ct_lT_T = -1d0, ct_lT = 0d0

  ! ---- charge exchange (kind 1 recombination, 2 ionisation; stage; element): cubic Hermite tables in log10 T (rtz_rate_tables)
  integer, parameter :: ctt_npd = 400                       ! nodes per dex
  real(dp), parameter :: ctt_x0 = 0d0, ctt_x1 = 9d0        ! log10 T [K] of the first and last node
  integer, parameter :: ctt_nk = 3600                       ! cells: (x1 - x0) * npd
  integer, parameter :: ctt_maxp = 96                         ! columns
  real(dp), allocatable, save :: ctt_tab(:,:,:)             ! (2, column, 0:nk): value, slope x cell width (dex)
  integer, save :: ctt_id(2,4,30) = 0                    ! pair -> column; 0 not built, -1 not tabulated
  integer, save :: ctt_np = 0
  real(dp), save :: ctt_tlo(96) = 0d0, ctt_thi(96) = huge(1d0)   ! the column's range in T
  real(dp), save :: ctt_vlo(96) = 0d0, ctt_vhi(96) = 0d0         ! its value below / above it
  real(dp), save :: ctt_T = -1d0, ctt_w(4) = 0d0           ! the Hermite weights at the last T
  integer, save :: ctt_k = 0
  real(dp), save :: ctt_errmax = 0d0                        ! the largest midpoint error found (rate > 1e-30)
  ! the pair of the column being checked (ctt_check's function argument)
  integer, save :: ctt_ck = 0, ctt_ci = 0, ctt_ce = 0

CONTAINS

SUBROUTINE load_ct_rates()
  !Load the charge transfer ionization
  !and recombination rates from file
  use amr_commons, only: myid
  use hydro_parameters, only: data_dir
  implicit none

  integer:: i, j, k, unit_num, ios

  if(myid.eq.1) write(*,*) 'Initializing charge exchange rates'

  ! Load charge exchange file
  open(newunit=unit_num, file=trim(data_dir)//'/charge_transfer/ct_ionization.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      if(myid==1) write(*,*) 'Error: Could not open CT Ionization file'
      return
  end if

  !Reading data from the file into the 3D array
  do i = 3, 30
     do j = 1, 3
         read(unit_num, *, iostat=ios) CTIon(:,j,i)
        if (ios /= 0) exit
     end do
     if (ios /= 0) exit
  end do
  close(unit_num)

  ! Zero out helium
  do j = 1, 3
     do k = 1, 7
        CTIon(k,j,2) = 0.d0
     end do
  end do

  ! Load Recombination
  open(newunit=unit_num, file=trim(data_dir)//'/charge_transfer/ct_recombination.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      if(myid==1) write(*,*) 'Error: Could not open CT Recombination file'
      return
  end if

  !Reading data from the file into the 3D array
  do i = 2, 30
     do j = 1, 4
         read(unit_num, *, iostat=ios) CTRecomb(:,j,i)
        if (ios /= 0) exit
     end do
     if (ios /= 0) exit
  end do
  close(unit_num)

END SUBROUTINE load_ct_rates

FUNCTION charge_transfer_recombination(ion, nelem, T) result(rate)
  use safe_math, only: safe_exp
  use rt_parameters, only: rtz_rate_tables
  ! ion is stage of ionization, 2 for the ion going to the atom
  ! nelem is atomic number of element, 2 up to 30
  ! Example:  O+ + H => O + H+ is HCTRecom(2,8,1e4)
  ! Note that temperature is in linear scale
  implicit none

  integer, intent(in):: ion, nelem
  real(dp), intent(in):: T
  real(dp):: rate
  real(dp):: a_op, b_op, c_op, d_op, e_op, f_op
  real(dp):: logT, tused
  integer:: ipIon

  rate = 0.0
  ! No recombination on the ground state
  if (ion .eq. 1) then
     return
  end if

  ipIon = ion - 1
  !use statistical charge transfer for ion > 5
  if (ion .gt. 5) then
     rate = 1.92d-9 * real(ipIon, dp)
     return
  end if

  if (rtz_rate_tables) then
     if (ctt_rate(1, ipIon, nelem, T, rate)) return
  end if

  if (ctr_T(ipIon,nelem) == T) then
     rate = ctr_rate(ipIon,nelem)
     return
  end if

  !Make sure T is between temp. boundaries
  tused = max(min(T,CTRecomb(6,ipIon,nelem)),CTRecomb(5,ipIon,nelem))
  tused = tused * 1d-4
  tused = max(tused,1d-10) ! harley added to prevent zero temperature

  ! The interpolation equation, a 1e-9 (T/1e4)^b (1 + c exp(d T/1e4)), with the power as exp(b log)
  ! and the terms of a zero coefficient skipped
  if (CTRecomb(1,ipIon,nelem) == 0d0) then
     rate = 1.d-100
  else
     rate = CTRecomb(1,ipIon,nelem) * 1d-9
     if (CTRecomb(2,ipIon,nelem) /= 0d0) rate = rate * safe_exp(CTRecomb(2,ipIon,nelem) * ct_log(T, tused))
     if (CTRecomb(3,ipIon,nelem) /= 0d0) rate = rate * (1.d0 + CTRecomb(3,ipIon,nelem) * safe_exp(CTRecomb(4,ipIon,nelem)*tused))
     rate = MAX(rate,1.d-100)
  end if
  ctr_T(ipIon,nelem) = T
  ctr_rate(ipIon,nelem) = rate

END FUNCTION charge_transfer_recombination

FUNCTION charge_transfer_ionization(ion, nelem, T) result(rate)
  use safe_math, only: safe_exp
  use rt_parameters, only: rtz_rate_tables
  ! ion is stage of ionization, 1 for atom
  ! nelem is atomic number of element, 2 up to 30
  ! Example:  O + H+ => O+ + H is HCTIon(1,8,1e4)
  ! Note that temperature is in linear scale
  implicit none
  
  integer, intent(in):: ion, nelem
  real(dp), intent(in):: T
  real(dp):: rate
  real(dp):: a, b, c
  real(dp):: a_o, b_o, c_o, d_o, e_o, f_o, g_o
  real(dp):: logT, tused
  integer:: ipIon

  rate = 0.d0

  ipIon = ion
  if (ipIon.gt.2) then 
     rate = 0.d0
     return
  end if

  if (rtz_rate_tables) then
     if (ctt_rate(2, ipIon, nelem, T, rate)) return
  end if

  if (cti_T(ipIon,nelem) == T) then
     rate = cti_rate(ipIon,nelem)
     return
  end if

  ! Make sure T is between temp. boundaries
  tused = max(min(T,CTIon(6,ipIon,nelem)),CTIon(5,ipIon,nelem))
  tused = tused * 1d-4
  tused = max(tused,1d-10) ! harley added to prevent zero temperature

  ! the interpolation equation, a 1e-9 (T/1e4)^b (1 + c exp(d T/1e4)) exp(-e/(T/1e4)), as above
  if (CTIon(1,ipIon,nelem) == 0d0) then
     rate = 1.d-100
  else
     rate = CTIon(1,ipIon,nelem) * 1d-9
     if (CTIon(2,ipIon,nelem) /= 0d0) rate = rate * safe_exp(CTIon(2,ipIon,nelem) * ct_log(T, tused))
     if (CTIon(3,ipIon,nelem) /= 0d0) rate = rate * (1.d0 + CTIon(3,ipIon,nelem) * safe_exp(CTIon(4,ipIon,nelem)*tused))
     if (CTIon(7,ipIon,nelem) /= 0d0) rate = rate * safe_exp(-1.d0 * CTIon(7,ipIon,nelem)/tused)
     rate = MAX(rate,1.d-100)
  end if
  cti_T(ipIon,nelem) = T
  cti_rate(ipIon,nelem) = rate

END FUNCTION charge_transfer_ionization

FUNCTION ctt_interp(c, T) result(v)
  ! column c at T (inside the table's T range): the cubic Hermite interpolant, weights shared at one T
  implicit none
  integer, intent(in) :: c
  real(dp), intent(in) :: T
  real(dp) :: v, x, s, s2, s3
  integer :: k
  if (T /= ctt_T) then
     x = (log10(T) - ctt_x0) * ctt_npd
     k = min(max(int(x), 0), ctt_nk - 1)
     s = x - k
     s2 = s * s
     s3 = s2 * s
     ctt_w(1) = 2d0 * s3 - 3d0 * s2 + 1d0
     ctt_w(2) = s3 - 2d0 * s2 + s
     ctt_w(3) = 3d0 * s2 - 2d0 * s3
     ctt_w(4) = s3 - s2
     ctt_k = k
     ctt_T = T
  end if
  k = ctt_k
  v = ctt_w(1) * ctt_tab(1,c,k) + ctt_w(2) * ctt_tab(2,c,k) + ctt_w(3) * ctt_tab(1,c,k+1) &
      + ctt_w(4) * ctt_tab(2,c,k+1)
END FUNCTION ctt_interp

SUBROUTINE ctt_check(c, f)
  ! the column's error at the cells' midpoints against the exact f, where f > 1e-30 (kept in ctt_errmax)
  implicit none
  integer, intent(in) :: c
  real(dp), external :: f
  integer :: k
  real(dp) :: Tm, ex, ap
  do k = 0, ctt_nk - 1
     Tm = 10d0**(ctt_x0 + (k + 0.5d0) / ctt_npd)
     if (Tm <= ctt_tlo(c) .or. Tm >= ctt_thi(c)) cycle
     ex = f(Tm)
     if (ex <= 1d-30) cycle
     ap = 0.5d0 * (ctt_tab(1,c,k) + ctt_tab(1,c,k+1)) + 0.125d0 * (ctt_tab(2,c,k) - ctt_tab(2,c,k+1))
     ctt_errmax = max(ctt_errmax, abs(ap / ex - 1d0))
  end do
END SUBROUTINE ctt_check

SUBROUTINE ct_coeffs(kind, ipIon, nelem, co)
  ! the fit's 7 coefficients (a, b, c, d, Tmin, Tmax, e: CTIon's 7th, 0 for recombination)
  implicit none
  integer, intent(in) :: kind, ipIon, nelem
  real(dp), intent(out) :: co(7)
  co = 0d0
  if (kind == 1) then
     co(1:6) = CTRecomb(1:6, ipIon, nelem)
  else
     co(1:7) = CTIon(1:7, ipIon, nelem)
  end if
END SUBROUTINE ct_coeffs

SUBROUTINE ct_fit(co, t4, v, dvdx)
  ! the fit at t4 = T/1e4 without its T clamp: a 1e-9 t4^b (1 + c exp(d t4)) exp(-e/t4), and d/dlog10 T
  use safe_math, only: safe_exp
  implicit none
  real(dp), intent(in) :: co(7), t4
  real(dp), intent(out) :: v, dvdx
  real(dp) :: ex, dl
  v = co(1) * 1d-9
  dl = 0d0
  if (co(2) /= 0d0) then
     v = v * safe_exp(co(2) * log(t4))
     dl = dl + co(2)
  end if
  if (co(3) /= 0d0) then
     ex = co(3) * safe_exp(co(4) * t4)
     v = v * (1d0 + ex)
     if (1d0 + ex /= 0d0) dl = dl + co(4) * t4 * ex / (1d0 + ex)
  end if
  if (co(7) /= 0d0) then
     v = v * safe_exp(-co(7) / t4)
     dl = dl + co(7) / t4
  end if
  dvdx = v * dl * log(10d0)            ! d v / d log10 T = v dln v/dln T ln 10
END SUBROUTINE ct_fit

FUNCTION ctt_column(kind, ipIon, nelem) result(c)
  ! the column of a pair, built at its first use; -1 when the fit is not tabulated (a = 0: the 1e-100 floor) or
  ! the table is full
  use amr_commons, only: myid
  implicit none
  integer, intent(in) :: kind, ipIon, nelem
  integer :: c, k
  real(dp) :: co(7), v, dv, Tk
  c = ctt_id(kind, ipIon, nelem)
  if (c /= 0) return
  call ct_coeffs(kind, ipIon, nelem, co)
  if (co(1) == 0d0 .or. ctt_np >= ctt_maxp) then
     ctt_id(kind, ipIon, nelem) = -1
     c = -1
     return
  end if
  if (.not. allocated(ctt_tab)) allocate(ctt_tab(2, ctt_maxp, 0:ctt_nk))
  ctt_np = ctt_np + 1
  c = ctt_np
  do k = 0, ctt_nk
     Tk = 10d0**(ctt_x0 + dble(k) / ctt_npd)
     call ct_fit(co, Tk * 1d-4, v, dv)
     ctt_tab(1, c, k) = v
     ctt_tab(2, c, k) = dv / ctt_npd   ! slope x cell width
  end do
  ! below / above the fit's range: the clamped value of the exact function
  ctt_tlo(c) = co(5)
  ctt_thi(c) = co(6)
  call ct_fit(co, max(co(5) * 1d-4, 1d-10), v, dv)
  ctt_vlo(c) = max(v, 1d-100)
  call ct_fit(co, max(co(6) * 1d-4, 1d-10), v, dv)
  ctt_vhi(c) = max(v, 1d-100)
  ctt_id(kind, ipIon, nelem) = c
  ctt_ck = kind
  ctt_ci = ipIon
  ctt_ce = nelem
  call ctt_check(c, ct_exact_at)
END FUNCTION ctt_column

FUNCTION ct_exact_at(Tx) result(v)
  ! the unclamped fit of the pair being checked, at Tx
  implicit none
  real(dp), intent(in) :: Tx
  real(dp) :: v, dv, co(7)
  call ct_coeffs(ctt_ck, ctt_ci, ctt_ce, co)
  call ct_fit(co, Tx * 1d-4, v, dv)
END FUNCTION ct_exact_at

LOGICAL FUNCTION ctt_rate(kind, ipIon, nelem, T, rate)
  ! .true. and the rate when the pair is tabulated and T inside the table
  implicit none
  integer, intent(in) :: kind, ipIon, nelem
  real(dp), intent(in) :: T
  real(dp), intent(out) :: rate
  integer :: c
  ctt_rate = .false.
  rate = 0d0
  if (T < 1d0 .or. T >= 1d9) return
  c = ctt_id(kind, ipIon, nelem)
  if (c == 0) c = ctt_column(kind, ipIon, nelem)
  if (c < 0) return
  ctt_rate = .true.
  if (T <= ctt_tlo(c)) then
     rate = ctt_vlo(c)
  else if (T >= ctt_thi(c)) then
     rate = ctt_vhi(c)
  else
     rate = max(ctt_interp(c, T), 1d-100)
  end if
END FUNCTION ctt_rate

SUBROUTINE ct_tables_prebuild()
  ! the columns of the run's elements at init, and the check of their accuracy
  use amr_commons, only: myid
  use rt_parameters, only: rtz_rate_tables
  use rtz_module, only: elements, n_elements
  implicit none
  integer :: iE, z, i, c
  if (.not. rtz_rate_tables) return
  do iE = 1, n_elements
     z = elements(iE)%atomic_number
     if (z < 2 .or. z > 30) cycle
     do i = 1, 4
        c = ctt_column(1, i, z)
     end do
     do i = 1, 2
        c = ctt_column(2, i, z)
     end do
  end do
  if (myid == 1) write(*,'(A,I0,A,ES9.2)') ' Charge-exchange rate tables (rtz_rate_tables): ', ctt_np, &
       ' columns, 400 nodes/dex over 1 K-1e9 K; largest midpoint error (rate > 1e-30) ', ctt_errmax
END SUBROUTINE ct_tables_prebuild

FUNCTION ct_log(T, tused) result(lt)
  ! log(tused), from the shared log(T/1e4) when the fit's T range did not clamp T
  implicit none
  real(dp), intent(in):: T, tused
  real(dp):: lt

  if (tused /= T * 1d-4) then
     lt = log(tused)
     return
  end if
  if (T /= ct_lT_T) then
     ct_lT_T = T
     ct_lT = log(tused)
  end if
  lt = ct_lT

END FUNCTION ct_log

end module charge_exchange_module