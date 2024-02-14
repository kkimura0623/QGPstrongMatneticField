module photon_pol_zero_b_class
!
! This conmputes the photon vacume polarization tensor functions 
! in zero magnetic field
!
! The lepton loop with mass m is contained.
!
  use constant_mod
  implicit none
  private
  public :: set_param_zero_b,N1_zero_b

  real(DP), save :: mom_r

contains

subroutine set_param_zero_b(r)
  implicit none
  real(DP), intent(in) :: r
  mom_r  = r
  return
end subroutine

function N1_below_th() result(f)
  implicit none
  complex(DP) :: f
  real(DP) :: fre,fim
  real(DP) :: y,sqrt_y

  y = 1.0_DP/mom_r
  sqrt_y = SQRT(y-1.0_DP)

  fre = (1.0_DP/3.0_DP + 2.0_DP*(1.0_DP + 0.5_DP*y)*(sqrt_y * atan(1.0_DP/sqrt_y) - 1.0_DP))/3.0_DP/PI
  fim = 0.0_DP
  f = CMPLX(fre,fim,kind=DP)

  return
end function

function N1_above_th() result(f)
  implicit none
  complex(DP) :: f
  real(DP) :: fre,fim
  real(DP) :: y,sqrt_y

  y = 1.0_DP/mom_r
  sqrt_y = SQRT(1.0_DP-y)

  fre = (1.0_DP/3.0_DP + 2.0_DP*(1.0_DP + 0.5_DP*y)*(sqrt_y * atanh(sqrt_y) - 1.0_DP))/3.0_DP/PI
  fim =                        -(1.0_DP + 0.5_DP*y)* sqrt_y/3.0_DP
  f = CMPLX(fre,fim,kind=DP)

  return
end function

function N1_zero_b() result(f)
  implicit none
  complex(DP) :: f

  if ( mom_r < 1.0_DP) then
    f = N1_below_th()
  else
    f = N1_above_th()
  endif

  return
end function

end module
