module safe_math
  implicit none
  private
  public :: safe_exp, safe_erf

  real(8), parameter :: EXP_HI = 700.d0    ! Upper bound to avoid overflow
  real(8), parameter :: EXP_LO = -700.d0   ! Lower bound to avoid slow denormals
  real(8), parameter :: ERF_HI = 6.d0      ! erf(6) is ~1 within double precision

CONTAINS

PURE ELEMENTAL FUNCTION safe_exp(x) result(y)
  use amr_parameters, only: dp
  implicit none
  ! A fast and safe exponential function for double precision
  real(dp), intent(in) :: x
  real(dp) :: y

  if (x > EXP_HI) then
    y = exp(EXP_HI)   ! ≈ 1e304, avoids overflow
  elseif (x < EXP_LO) then
    y = 0.d0          ! Avoids denormals and underflow
  else
    y = exp(x)
  end if
END FUNCTION safe_exp

PURE ELEMENTAL FUNCTION safe_erf(x) result(y)
  use amr_parameters, only: dp
  implicit none
  ! A safe error function with saturation in extreme tails
  real(dp), intent(in) :: x
  real(dp) :: y

  if (x > ERF_HI) then
    y = 1.d0
  elseif (x < -ERF_HI) then
    y = -1.d0
  else
    y = erf(x)
  end if
END FUNCTION safe_erf

end module safe_math
