module photon_pol_parallel_b_class
!
! This conmputes the photon vacume polarization tensor functions 
! in a constant magnetic field.
!
! The lepton loop with mass m is contained.
!
! The photon virtural momenta should be parallel to the magnetic field B.
!
!
  use constant_mod
  use intde2_mod
  implicit none
  private

  public :: N0_parallel_b
  public :: N1_parallel_b
  public :: N2_parallel_b


  real(DP), parameter :: TOL   = 1.0e-14_DP
  real(DP), parameter :: DTINY = 1.0e-305_DP
  integer,  parameter :: NWORK = 8000

  type form_factor_parallel_b
    real(DP) :: mag_mu    = 0.0_DP
    real(DP) :: mom_q     = 0.0_DP
    real(DP) :: mom_r     = 0.0_DP
    real(DP) :: feynman_v = 0.0_DP
    integer :: Kceil
    integer :: Jceil
  end type

  real(DP) :: DIGAMA
  external :: DIGAMA

contains


function N0B(x, c_paras) result(f)
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_parallel_b) :: paras
  real(DP) :: f
  real(DP) :: ph,vv,onemvv
  integer :: ierr

  paras = TRANSFER(c_paras,paras)

  vv = x*x
  onemvv = 1.0_DP-vv

  ph = (1.0_DP - onemvv*paras%mom_r)/(paras%mag_mu)

  f = - 2*(vv + onemvv*LOG(2.0_DP*paras%mag_mu))  &
 &    -(onemvv + ph*x)*DIGAMA((1.0_DP + ph - x)*0.5_DP + REAL(paras%Kceil+1,kind=DP),ierr) &
 &    -(onemvv - ph*x)*DIGAMA((1.0_DP + ph + x)*0.5_DP + REAL(paras%Kceil+1,kind=DP),ierr)

  f = f*0.5_DP

  return
end function

function N1B(x, c_paras) result(f)
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_parallel_b) :: paras
  real(DP) :: f
  real(DP) :: ph,vv
  integer :: ierr

  paras = TRANSFER(c_paras,paras)

  vv = 1.0_DP-x*x

  ph = (1.0_DP - vv*paras%mom_r)/(2*paras%mag_mu)

  f = vv*(-log(2*paras%mag_mu) - DIGAMA(ph + REAL(paras%Kceil + 1,kind=DP),ierr))

  return
end function

function N2B(x, c_paras) result(f)
  implicit none
  real(DP),         intent(in) :: x
  character(len=1), intent(in) :: c_paras(:)
  type(form_factor_parallel_b) :: paras
  real(DP) :: f
  real(DP) :: ph,vv,onemvv
  integer :: ierr

  paras = TRANSFER(c_paras,paras)

  vv = x*x
  onemvv = 1.0_DP-vv

  ph = (1.0_DP - onemvv*paras%mom_r)/(paras%mag_mu)

  f = -1.0_DP + 2*ph - 3*vv - onemvv*2*LOG(2.0_DP*paras%mag_mu)  &
 &           -2*ph**2      *DIGAMA((2.0_DP + ph    )*0.5_DP + REAL(paras%Jceil+1,kind=DP),ierr) &
 &    -(1.0_DP-(ph - x)**2)*DIGAMA((1.0_DP + ph - x)*0.5_DP + REAL(paras%Kceil+1,kind=DP),ierr) &
 &    -(1.0_DP-(ph + x)**2)*DIGAMA((1.0_DP + ph + x)*0.5_DP + REAL(paras%Kceil+1,kind=DP),ierr)

  f = f*0.5_DP

  return
end function

function integ_X(r,mu,j) result(f)
!
! f = int{v,-1,1}  (1-v^2)/(r*v^2 +(1-r+2*mu*j) - i*eps)
!
  use pole_integ_func_class
  implicit none
  real(DP), intent(in) :: r,mu
  integer,  intent(in) :: j
  complex(DP) :: f,ff,hh
#ifdef _CHECK_
  complex(DP) :: xx
  real(DP) :: bb,xr,xi
#endif
  integer(8) :: jj

  jj = j
  ff = integ_f(0_8,jj,r,mu)
  hh = integ_h(0_8,jj,r,mu,ff)
  f = ff - hh

#ifdef _CHECK_
  bb = -(1.0_DP - mom_r + 2*j*mag_mu)/mom_r

  if ( bb <= 0.0_DP ) then

    xr = ( ((1.0_DP-bb)/SQRT(-bb))*atan(1.0_DP/SQRT(-bb))-1.0_DP)*2.0_DP/mom_r
    xi = 0.0_DP

  else if ( bb <= 1.0_DP .and. 0.0_DP < bb ) then

    xr = (- ((1.0_DP-bb)/SQRT(bb))*atanh(SQRT(bb))-1.0_DP)*2.0_DP/mom_r
    xi = PIQ*(1.0_DP-bb)/SQRT(bb)/mom_r

  else if ( 1.0_DP < bb ) then

    xr = (- ((1.0_DP-bb)/SQRT(bb))*atanh(SQRT(bb))-1.0_DP)*2.0_DP/mom_r
    xi = 0.0_DP

  endif
  xx = CMPLX(xr,xi,kind=DP)

  write(*,'("integ_X=",I4," r=",E24.15," mu=",E24.15," f,xx,bb=",8E24.15)')j,r,mu,f-xx
#endif

  return
end function

function N0_parallel_b(r,q,mu) result(f)
  use pole_integ_func_class
  implicit none
  real(DP), intent(in) :: r,q,mu
  complex(DP) :: f
  complex(QP) :: ff
  real(DP) :: a,b,result
  real(DP) :: aw(0:NWORK-1),err
  real(DP) :: fre,fim,bb,b0,db,sqrt_bb,fac
  real(DP) :: Amin
  integer(8) :: jj
  integer :: j
  character(len=1), allocatable :: c_paras(:)
  type(form_factor_parallel_b) :: paras
  integer :: length

  call intdeini(NWORK,DTINY,TOL,aw)

  Amin = (1.0_DP - r + mu - (mu**2/(4*r)) )/(2*mu)
  if ( ABS(mu/(2*r)) <= 1.0_DP ) then
    if (Amin <= 1) then
      paras%Kceil = -CEILING(Amin) + 1
    else
      paras%Kceil = -1
    endif
  else
    paras%Kceil = -1
  endif

  paras%mag_mu = mu
  paras%mom_r  = r
  paras%mom_q  = q

  length = SIZE(TRANSFER(paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(paras,c_paras)

  f = Z0
  a =  0.0_DP
  b =  1.0_DP
  call intde(N0B,c_paras,a,b,aw,result,err)
  if ( err < 0.0_DP ) then
    write(*,'(" n0_parallel_b : INTDE : ERROR STOP : ",ES24.15)')err
    stop
  endif
  f = result*2

  do j=0,paras%Kceil
    jj = j
    ff = integ_f(1_8,jj,r,mu)
    f = f +2*(       +  r * integ_i(1_8,jj,r,mu,ff)   &
 &                   - mu * integ_h(1_8,jj,r,mu,ff)   &
 &           +(1.0_DP - r)* integ_g(1_8,jj,r,mu,ff)   &
 &                   + mu * ff                     )
  enddo

  f = -f/(4.0_DP*PI)

  deallocate(c_paras)

  return
end function

function N1_parallel_b(r,q,mu) result(f)
  implicit none
  real(DP), intent(in) :: r,q,mu
  complex(DP) :: f
  real(DP) :: a,b,result
  real(DP) :: aw(0:NWORK-1),err
  real(DP) :: fre,fim,bb,b0,db,sqrt_bb,fac
  real(DP) :: Amin
  integer :: j
  character(len=1), allocatable :: c_paras(:)
  type(form_factor_parallel_b) :: paras
  integer :: length

  call intdeini(NWORK,DTINY,TOL,aw)

  f = Z0

  a =  0.0_DP
  b =  1.0_DP
!  if ( r < 1.0_DP) then
!    paras%Kceil = -1
!  else
!    paras%Kceil = -CEILING((1.0_DP-r)/(2*mu))
!  endif
  Amin = (1.0_DP-r)/(2*mu)
  if (Amin <= 1.0_DP) then
    paras%Kceil = -CEILING(Amin) + 1
  else
    paras%Kceil = -1
  endif

  paras%mag_mu = mu
  paras%mom_r  = r
  paras%mom_q  = q

  length = SIZE(TRANSFER(paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(paras,c_paras)

  call intde(N1B,c_paras,a,b,aw,result,err)
  if ( err < 0.0_DP ) then
    write(*,'(" n1_parallel_b : INTDE : ERROR STOP : ",ES24.15)')err
    stop
  endif
  f = result*2
  f = f - mu*integ_X(r,mu,0)

  do j=0,paras%Kceil
    f = f + 2*mu*integ_X(r,mu,j)
  enddo

  f = -f/(4.0_DP*PI)

  deallocate(c_paras)

  return
end function

function N2_parallel_b(r,q,mu) result(f)
  use pole_integ_func_class
  implicit none
  real(DP), intent(in) :: r,q,mu
  complex(DP) :: f
  complex(QP) :: ff
  real(DP) :: a,b,result
  real(DP) :: aw(0:NWORK-1),err
  real(DP) :: fre,fim,bb,b0,db,sqrt_bb,fac
  real(DP) :: Amin
  integer(8) :: jj
  integer :: j
  character(len=1), allocatable :: c_paras(:)
  type(form_factor_parallel_b) :: paras
  integer :: length

  call intdeini(NWORK,DTINY,TOL,aw)

  Amin = (1.0_DP - r + mu - (mu**2/(4*r)) )/(2*mu)
  if ( ABS(mu/(2*r)) <= 1.0_DP ) then
!    if (Amin <= 0) then
    if (Amin <= 1.0_DP) then
!      paras%Kceil = -CEILING(Amin)
      paras%Kceil = -CEILING(Amin)+1
    else
      paras%Kceil = -1
    endif
  else
    paras%Kceil = -1
  endif

  Amin = (1.0_DP - r + 2*mu)/(2*mu)
  if (Amin <= 0) then
    paras%Jceil = -CEILING(Amin)
  else
    paras%Jceil = -1
  endif

  paras%mag_mu = mu
  paras%mom_r  = r
  paras%mom_q  = q

  length = SIZE(TRANSFER(paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(paras,c_paras)

  f = Z0
  a =  0.0_DP
  b =  1.0_DP
  call intde(N2B,c_paras,a,b,aw,result,err)
  if ( err < 0.0_DP ) then
    write(*,'(" n2_parallel_b : INTDE : ERROR STOP : ",ES24.15)')err
    stop
  endif
  f = result*2

  do j=0,paras%Kceil
    jj = j
    ff = integ_f(1_8,jj  ,r,mu)
    f = f + 2*(-(2.0_DP - 4*r/3.0_DP)/mu + 2*(2*j+1) - 4* j*(j+1)*mu*ff )
  enddo
  do j=0,paras%Jceil
    jj = j
    ff = integ_f(0_8,jj+1,r,mu)
    f = f + 2*( (2.0_DP - 4*r/3.0_DP)/mu - 4*(j+1)   + 4*(j+1)**2*mu*ff )
  enddo

  f = -f/(4.0_DP*PI)

  deallocate(c_paras)

  return
end function

end module
