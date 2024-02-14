module photon_pol_below_th_class
!
! This conmputes the photon vacume polarization tensor functions 
! in a constant magnetic field.
!
! The lepton loop with mass m is contained.
!
! The photon virtural momenta should be less than threshold k^2 < (2m)^2.
!
!
  use constant_mod
  use intde2_mod
  implicit none
  private
  public :: N0_below_th
  public :: N1_below_th
  public :: N2_below_th

  real(DP), save :: mag_mu, mom_q, mom_r, feynman_v

  type integ_paras
    real(DP) :: mag_mu    = 0.0_DP
    real(DP) :: mom_q     = 0.0_DP
    real(DP) :: mom_r     = 0.0_DP
    real(DP) :: feynman_v = 0.0_DP
  end type

  type form_factor_below_th
    type(integ_paras), pointer :: paras => NULL()
  end type

  real(DP), parameter :: TOL = 1.0e-14_DP
  real(DP), parameter :: DTINY = 1.0e-305_DP
  integer,  parameter :: NWORK = 8000

contains

subroutine set_param_below_th(r,q,mu)
  implicit none
  real(DP), intent(in) :: r,q,mu

  if ( r > 1.0_DP ) then
    write(*,'(" Parallel virtural momenta should be less than the threshold: r=k_parallel^2/(2m)^2 < 1.")')
    stop
  endif
  mom_r  = r
  mom_q  = q
  mag_mu = mu
  return
end subroutine

function phase(x, paras) result(ph)
!
! Compute phase function
!
! phase(x,v,r,q,mu) = ((1-(1-v^2)*r)/mu)*x + (cosh(x)-cosh(v*x))/sinh(x) * (2*q/mu)
!
! r =  k_parallel^2/(2m)^2   : Minkovski vector norm
!
! q = vec{k}_perp^2/(2m)^2   : Euclid vector norm > 0
!
! mu = eB/m^2
!
!
  implicit none
  real(DP),          intent(in) :: x
  type(integ_paras), intent(in) :: paras
  real(DP) :: ph,tph,vv,xx,c6,c5,c4,c3,c2,c1

  if ( ABS(x) < DTINY ) then
    ph = 0.0_DP
    return
  endif

  vv = paras%feynman_v*paras%feynman_v

  if ( x < 0.02_DP ) then

    !
    ! Use 11-th order taylor expansion for (cosh(x)-cosh(v*x))/sinh(x)
    ! for small x
    !

    xx = x*x

    c6 = (  2073.0_DP + vv*(-964.0_DP + vv*(  190.0_DP + vv*(-20.0_DP + vv))))   &
 &      /(-20460.0_DP + vv*(9372.0_DP + vv*(-1716.0_DP + 132.0_DP*vv)))

    c5 = (-155.0_DP + vv*(71.0_DP + vv*(-13.0_DP + vv))) &
 &      /(1530.0_DP + vv*(-660.0_DP + 90.0_DP*vv))

    c4 = (51.0_DP + vv*(-22.0_DP + 3.0_DP*vv))/(-504.0_DP + 168.0_DP*vv)
    c3 =  -0.1_DP + vv/30.0_DP
    c2 = (-1.0_DP + vv)/12.0_DP
    c1 = ( 1.0_DP - vv)/2.0_DP

    tph =  c1*x*( 1.0_DP + c2*xx*( 1.0_DP &
 &      + c3*xx*( 1.0_DP + c4*xx*( 1.0_DP &
 &      + c5*xx*( 1.0_DP + c6*xx)))))

  else if ( 10.0_DP < x ) then

    !
    ! Use exponential form for (cosh(x) - cosh(v*x))/sinh(x) 
    ! for large x
    ! 

    c3 = exp( (paras%feynman_v-1.0_DP)*x)
    c2 = exp(-(paras%feynman_v+1.0_DP)*x)
    c1 = exp(-2.0_DP*x)

    tph = ( 1.0_DP + (c1 - c3 - c2) )/(1.0_DP - c1)

  else

    tph = ( cosh(x) - cosh(paras%feynman_v*x) )/sinh(x)

  endif

  ph = ( ( 1.0_DP - (1.0_DP - vv)*paras%mom_r )*x  &
 &      +( tph*2.0_DP*paras%mom_q ) )/paras%mag_mu
  
  return
end function

function N1A(x, c_paras) result(f)
!
! N1A(x) = exp(-phase(x))*cosh(x)/sinh(x) - exp(-x/mu)/x
!
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_below_th) :: p_paras
  real(DP) :: f,ct,xx,c1

  if ( ABS(x) < DTINY ) then
    f = 0.0_DP
    return
  endif

  p_paras = TRANSFER(c_paras,p_paras)

  if ( x < 0.02_DP ) then
    !
    ! Use 9-th order taylor expansion for cosh(x)/sinh(x)
    ! for small x
    !

    xx = x*x

    ct =  ( 1.0_DP + xx*(1.0_DP - xx*(1.0_DP  &
 &          - 2.0_DP*xx*(1.0_DP - xx*(1.0_DP  &
 &         - 10.0_DP*xx/99.0_DP)/10.0_DP)/21.0_DP)/15.0_DP)/3.0_DP)/x

  else if ( 10.0_DP < x ) then
    !
    ! Use exponential form for cosh(x)/sinh(x)
    ! for large x
    !
    c1 = exp(-2.0_DP*x)

    ct = (1.0_DP + c1)/(1.0_DP - c1)

  else

    ct = cosh(x)/sinh(x)

  endif

  f = exp(-phase(x,p_paras%paras))*ct - exp(-x/p_paras%paras%mag_mu)/x

  return
end function

function N1B(x, c_paras) result(f)
!
! v integrated :(v<=x)
!
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_below_th) :: p_paras
  real(DP) :: f
  real(DP) :: a,result
  real(DP) :: aw(0:NWORK-1),err

  call intdeiini(NWORK,DTINY,TOL,aw)
  p_paras = TRANSFER(c_paras,p_paras)
  p_paras%paras%feynman_v = x
  a = 0.0_DP
  call intdei(N1A,c_paras,a,aw,result,err)

  f = result*(1.0_DP - x**2)

  return
end function

function N1_below_th(r,q,mu) result(f)
  implicit none
  real(DP), intent(in) :: r,q,mu
  real(DP) :: f
  real(DP) :: a,b,result
  real(DP) :: aw(0:NWORK-1),err
  type(form_factor_below_th) :: p_paras
  character(len=1), allocatable :: c_paras(:)
  type(integ_paras), target  :: paras
  integer :: length

  paras%mom_r  = r
  paras%mom_q  = q
  paras%mag_mu = mu

  p_paras%paras => paras
  length = SIZE(TRANSFER(p_paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(p_paras,c_paras)

  call intdeini(NWORK,DTINY,TOL,aw)
  a = 0.0_DP
  b = 1.0_DP
  call intde(N1B,c_paras,a,b,aw,result,err)

  f = -result*0.5_DP/PI

  p_paras%paras => NULL()
  deallocate(c_paras)

  return
end function

function N0A(x, c_paras) result(f)
!
! N0A(x) = exp(-phase(x))*(cosh(v*x)/sinh(x) - v * cosh(x)*sinh(v*x)/sinh(x)**2) - (1-v**2)*exp(-x/mu)/x
!
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_below_th) :: p_paras
  real(DP) :: f,ct,xx,c6,c5,c4,c3,c2,c1,vv

  if ( ABS(x) < DTINY ) then
    f = 0.0_DP
    return
  endif

  p_paras = TRANSFER(c_paras,p_paras)

  vv = p_paras%paras%feynman_v*p_paras%paras%feynman_v

  if ( x < 0.02_DP ) then
    !
    ! Use 9-th order taylor expansion for cosh(v*x)/sinh(x) - v *cosh(x)*sinh(v*x)/sinh(x)**2
    ! for small x
    !
    xx = x*x


    c6 = (-2555.0_DP + vv*(-12977.0_DP + vv*(  6130.0_DP + vv*(-866.0_DP + vv*(25.0_DP + 3.0_DP*vv)))))  &
 &      /(25146.0_DP + vv*( 78408.0_DP + vv*(-29172.0_DP + vv*(1320.0_DP + 330.0_DP*vv))))

    c5 = (  381.0_DP + vv*( 1188.0_DP + vv*(-442.0_DP + vv*(20.0_DP + 5.0_DP*vv))))  &
 &      /(-3720.0_DP + vv*(-4680.0_DP + 360.0_DP*vv*( 1.0_DP + vv)))

    c4 = (-31.0_DP + vv*(-39.0_DP + 3.0_DP*vv*(1.0_DP + vv)))/(294.0_DP + vv*(-84.0_DP + 126.0_DP*vv))
    c3 = (  7.0_DP + vv*( -2.0_DP + 3.0_DP*vv))/(60.0_DP*(-1.0_DP + vv))
    c2 = (-1.0_DP + vv)/6.0_DP
    c1 =   1.0_DP - vv

    ct = c1*(1.0_DP+c2*xx*(1.0_DP+c3*xx*(1.0_DP+c4*xx*(1.0_DP+c5*xx*(1.0_DP+c6*xx)))))/x

  else if ( 10.0_DP < x ) then
    !
    ! Use exponential form for (cosh(v*x)-v*cosh(x)*sinh(v*x)/sinh(x))/sinh(x)
    ! for large x
    !
    c1 = exp( (p_paras%paras%feynman_v-1.0_DP)*x)
    c2 = exp(-(p_paras%paras%feynman_v+1.0_DP)*x)
    c3 = exp(-2.0_DP*x)
    c4 = 1.0_DP - c3

    ct = ( (c1 + c2)  - p_paras%paras%feynman_v * (1.0_DP + c3)*(c1 - c2)/c4)/c4

  else

    ct = (cosh(p_paras%paras%feynman_v*x)  &
 &       - p_paras%paras%feynman_v * cosh(x)*sinh(p_paras%paras%feynman_v*x)/sinh(x))/sinh(x)

  endif

  f = exp(-phase(x,p_paras%paras))*ct-(1.0_DP-vv)*exp(-x/p_paras%paras%mag_mu)/x

  return
end function

function N0B(x, c_paras) result(f)
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_below_th) :: p_paras
  real(DP) :: f
  real(DP) :: a,result
  real(DP) :: aw(0:NWORK-1),err

  call intdeiini(NWORK,DTINY,TOL,aw)
  p_paras = TRANSFER(c_paras,p_paras)
  p_paras%paras%feynman_v = x
  a = 0.0_DP
  call intdei(N0A,c_paras,a,aw,result,err)

  f = result

  return
end function

function N0_below_th(r,q,mu) result(f)
  implicit none
  real(DP), intent(in) :: r,q,mu
  real(DP) :: f
  real(DP) :: a,b,result
  real(DP) :: aw(0:NWORK-1),err
  type(form_factor_below_th) :: p_paras
  character(len=1), allocatable :: c_paras(:)
  type(integ_paras), target  :: paras
  integer :: length

  paras%mom_r  = r
  paras%mom_q  = q
  paras%mag_mu = mu

  p_paras%paras => paras
  length = SIZE(TRANSFER(p_paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(p_paras,c_paras)

  call intdeini(NWORK,DTINY,TOL,aw)
  a = 0.0_DP
  b = 1.0_DP
  call intde(N0B,c_paras,a,b,aw,result,err)

  f = -result*0.5_DP/PI

  p_paras%paras => NULL()
  deallocate(c_paras)

  return
end function

function N2A(x, c_paras) result(f)
!
! N2A(x) = exp(-phase(x))*(-2*(cosh(v*x)-cosh(x))/sinh(x)**3) - (1-v**2)*exp(-x/mu)/x
!
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_below_th) :: p_paras
  real(DP) :: f,ct,xx,c6,c5,c4,c3,c2,c1,vv

  if ( ABS(x) < DTINY ) then
    f = 0.0_DP
    return
  endif

  p_paras = TRANSFER(c_paras,p_paras)

  vv = p_paras%paras%feynman_v * p_paras%paras%feynman_v

  if ( x < 0.02_DP ) then
    !
    ! Use 9-th order taylor expansion for 2*(cosh(x)-cosh(v*x))/(sinh(x)**3)
    ! for small x
    !
    xx = x*x

    c6 = (-119125.0_DP + vv*(  89981.0_DP + vv*(-18490.0_DP + vv*( 1618.0_DP + vv*(-65.0_DP + vv)))))  &
 &      /( 786852.0_DP + vv*(-514800.0_DP + vv*( 88440.0_DP + vv*(-5808.0_DP + 132.0_DP*vv))))

    c5 = (  5961.0_DP + vv*(-3900.0_DP + vv*(  670.0_DP + vv*(-44.0_DP + vv))))  &
 &      /(-35850.0_DP + vv*(18990.0_DP + vv*(-2430.0_DP + 90.0_DP*vv)))

    c4 = (-1195.0_DP + vv*(633.0_DP + vv*(-81.0_DP + 3.0_DP*vv)))/(6216.0_DP + vv*(-2352.0_DP + 168.0_DP*vv))
    c3 = (37.0_DP + vv*(-14.0_DP + vv))/(30.0_DP*(-5.0_DP + vv))
    c2 = (-5.0_DP + vv)/12.0_DP
    c1 =   1.0_DP - vv

    ct = c1*(1.0_DP+c2*xx*(1.0_DP+c3*xx*(1.0_DP+c4*xx*(1.0_DP+c5*xx*(1.0_DP+c6*xx)))))/x

  else if ( 10.0_DP < x ) then
    !
    ! Use exponential form for 2*(cosh(x)-cosh(v*x))/(sinh(x)**3)
    ! for large x
    !
    c1 = exp( (p_paras%paras%feynman_v-1.0_DP)*x)
    c2 = exp(-(p_paras%paras%feynman_v+1.0_DP)*x)
    c3 = exp(-2.0_DP*x)
    c4 = 1.0_DP - c3

    ct = 8.0_DP*c3*(1.0_DP + (c3 - c1 - c2) )/c4**3

  else

    ct = 2.0_DP*((cosh(x) - cosh(p_paras%paras%feynman_v*x))/(sinh(x)**3))

  endif

  f = exp(-phase(x,p_paras%paras))*ct - (1.0_DP-vv)*exp(-x/p_paras%paras%mag_mu)/x

  return
end function

function N2B(x, c_paras) result(f)
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  real(DP) :: f
  real(DP) :: a,result
  real(DP) :: aw(0:NWORK-1),err
  type(form_factor_below_th) :: p_paras

  call intdeiini(NWORK,DTINY,TOL,aw)
  p_paras = TRANSFER(c_paras,p_paras)
  p_paras%paras%feynman_v = x
  a = 0.0_DP
  call intdei(N2A,c_paras,a,aw,result,err)

  f = result

  return
end function

function N2_below_th(r,q,mu) result(f)
  implicit none
  real(DP), intent(in) :: r,q,mu
  real(DP) :: f
  real(DP) :: a,b,result
  real(DP) :: aw(0:NWORK-1),err
  type(form_factor_below_th) :: p_paras
  character(len=1), allocatable :: c_paras(:)
  type(integ_paras), target  :: paras
  integer :: length

  paras%mom_r  = r
  paras%mom_q  = q
  paras%mag_mu = mu

  p_paras%paras => paras
  length = SIZE(TRANSFER(p_paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(p_paras,c_paras)

  call intdeini(NWORK,DTINY,TOL,aw)
  a = 0.0_DP
  b = 1.0_DP
  call intde(N2B,c_paras,a,b,aw,result,err)

  f = -result*0.5_DP/PI

  p_paras%paras => NULL()
  deallocate(c_paras)

  return
end function

end module
