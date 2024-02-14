module photon_pol_strong_b_class
  use constant_mod
  implicit none
  private
  public :: N1_strong_b

contains

function N1_strong_b(r,q,mu) result(f)
!
! Strong B limit
!
! f = (mu/r)*exp(-2*q/mu)*(1 - arcsin(sqrt(r))/sqrt(r)/sqrt(1-r))/(2*pi)
!
!  for r < 1
!
! f = (mu/r)*exp(-2*q/mu)*(1 - (log(|sqrt(r)-sqrt(r-1)|/|sqrt(r)+sqrt(r-1)|) + I*pi)/2/sqrt(r)/sqrt(r-1))/(2*pi)
!
!  for r > 1
!
  implicit none
  real(DP), intent(in) :: r,q,mu
  complex(DP) :: f
  real(DP) :: sqrt_r,sqrt_rr,fac

  if ( r < 0.01_DP ) then

    fac = mu*exp(-2.0_DP*q/mu)/(2.0_DP*PI)

    f =  -0.66666666666666666667_DP +     &
 &      r*(-0.53333333333333333333_DP +     &
 &        r*(-0.45714285714285714286_DP +     &
 &          r*(-0.40634920634920634921_DP +     &
 &            r*(-0.36940836940836940837_DP +     &
 &              r*(-0.34099234099234099234_DP +     &
 &                r*(-0.31825951825951825952_DP +     &
 &                  r*(-0.29953837012660542072_DP +     &
 &                    r*(-0.2837731927515209249_DP +      &
 &                      r*(-0.27026018357287707133_DP -     &
 &                          0.25850974080883893779_DP*r)))))))))
    f = f*fac
    return

  endif

  fac = (mu/r)*exp(-2.0_DP*q/mu)/(2.0_DP*PI)
  sqrt_r = sqrt(r)

  if ( r < 1.0_DP) then

    f = fac*(1.0_DP - ASIN(sqrt_r)/sqrt_r/sqrt(1.0_DP-r))

  else

    sqrt_rr = sqrt(r-1.0_DP)

    f = fac*(1.0_DP - (LOG(abs(sqrt_r - sqrt_rr)/abs(sqrt_r + sqrt_rr)) + ZJ*PI )/sqrt_r/sqrt_rr/2.0_DP)

  endif

  return
end function

end module
