module pole_integ_func_class
  use constant_mod
  implicit none
  private
  public :: integ_f
  public :: integ_g
  public :: integ_h
  public :: integ_i

  real(QP) :: DELTA=1.0e-3_QP

contains

function integ_f(n, m, r, mu) result(f)
!
! returns
!
!  F^n_m(r,mu) =
!
! (i/mu) Integral_{v,-1,+1} Intelral_{t,0,Infinity} 
!              Exp[ -I ( phi + 2*m - n v + n ) t ]
!
!  = Integral_{v,-1,1}1/(r*v^2 -(n*mu)*v+1-r+(2*m+n)*mu)
!
!  phi = (1-(1-v^2)*r)/mu
!
!  r = k_perp^2/(4 me^2)
! mu = eB/me^2
!
  integer(8), intent(in) :: n,m
  real(QP), intent(in) :: r,mu
  complex(QP) :: f

  real(QP) :: a,b,c,cp,d,sqr,sp,sm, s1,s2, fr,fi
  real(QP) :: sfac,theta,tt
  real(QP) :: apcmb,apcpb

  f = Q0

  s1 = SQRT(1.0_QP+2* m   *mu)
  s2 = SQRT(1.0_QP+2*(m+n)*mu)
  sp = (s1+s2)**2/4.0_QP
!  sm = (s1-s2)**2/4.0_QP
  sm = ((n*mu)/(s1+s2))**2
!write(*,'(2I9,4ES24.16)')n,m,r,mu,sm,((n*mu)/(s1+s2))**2

!write(*,'(2I9,4ES24.15)')n,m,r,mu,sm,sp

  a = r
  b = -n*mu
  c = 1.0_QP - r + (2*m+n)*mu

!  d = b*b - 4*a*c
!  sqr = SQRT(d)
!  if (sqr < 1.0e-4_QP) then
!    write(*,'(" sqr=",8ES24.15)')sqr,a,b,c,d
!    write(*,'(" d=",8ES24.15)')d,a,b,c,d
!  endif

  if ( r <= sm ) then

    if ( n > 0 ) then
      tt = a/b
      if ( ABS(tt) < DELTA ) then
        f = integ_f_taylor1(n,m,r,mu)
        return
      endif
    else
      cp = 1.0_QP + 2*m*mu
      tt = a/cp
      if ( ABS(tt) < DELTA ) then
        f = integ_f_taylor1(n,m,r,mu)
        return
      endif
    endif

    d = b*b - 4*a*c
    sqr = SQRT(d)

!    sfac = SQRT((r-sm)*(r-sp))
!    fr = LOG(ABS((a-c-sqr)/(a-c+sqr))) /sfac/2.0_QP

    if ( ABS(a-c-sqr) < DELTA ) then

      apcmb =  1.0_QP + (2*(m+n))*mu   ! a + c - b
      apcpb =  1.0_QP + (2* m   )*mu   ! a + c + b

      fr = LOG(ABS(apcmb*apcpb/(a-c+sqr)**2)) /sqr

    else

      fr = LOG(ABS((a-c-sqr)/(a-c+sqr))) /sqr

    endif

    fi = 0.0_QP

  else if ( sm < r .and. r <= sp) then

    if ( n == 0 ) then
      cp = 1.0_QP + 2*m*mu
      tt = a/cp
      if ( ABS(tt) < DELTA ) then
        f = integ_f_taylor2(n,m,r,mu)
        return
      endif
    endif

    d = 4*a*c - b*b
    sqr = SQRT(d)

!    sfac = SQRT(ABS((r-sm)*(r-sp)))
!    sfac = sqr*0.5_DP
!    fr = ( atan((b+2*a)/sqr) - atan((b-2*a)/sqr) ) /sfac
!    fr = 2*( atan((b+2*a)/sqr) - atan((b-2*a)/sqr) ) /sqr

    if ( ABS(c-a) > DELTA ) then 
    ! 
    ! use addition formula for arctan
    !
      fr = atan(sqr/(c-a))
      if ( c < a ) then
        fr = fr + PIQ*SIGN(1.0_QP,(b+2*a))
      endif
      fr = 2*fr/sqr
    else
      fr = 2*( atan((b+2*a)/sqr) - atan((b-2*a)/sqr) ) /sqr
    endif

    fi = 0.0_QP

  else if ( sp < r ) then

    d = b*b - 4*a*c
    sqr = SQRT(d)
!    sfac = SQRT((r-sm)*(r-sp))
!    fr = LOG(ABS((a-c-sqr)/(a-c+sqr))) /sfac/2.0_QP
!    fi = PIQ/sfac
    if ( ABS(a-c-sqr) < DELTA ) then

      apcmb =  1.0_QP + (2*(m+n))*mu   ! a + c - b
      apcpb =  1.0_QP + (2* m   )*mu   ! a + c + b

!write(*,'(" n,m=",2I15," r=",ES24.15," 1-r=",ES24.15)')n,m,r,1.0_QP-r
!write(*,'(" sqr=",20ES24.15)')sqr,d,a,b,c,d,a-c-sqr,a-c+sqr
!write(*,'(" log=",20ES24.15)')LOG(ABS(apcmb*apcpb/(a-c+sqr)**2))

      fr = LOG(ABS(apcmb*apcpb/(a-c+sqr)**2)) /sqr

    else

      fr = LOG(ABS((a-c-sqr)/(a-c+sqr))) /sqr

    endif
    fi = 2*PIQ/sqr

  endif

  f = CMPLX(fr,fi,kind=QP)

  return
end function

function integ_f_taylor1(n, m, r, mu) result(f)
!
! returns
!  F^n_m(r,mu)  with small r  in r < sm
!
  integer(8), intent(in) :: n,m
  real(QP), intent(in) :: r,mu
  complex(QP) :: f

  real(QP) :: a,b,c,d,sqr,sp,sm, s1,s2, fr,fi
  real(QP) :: sfac,theta,ff(0:8),tt,w

  f = Q0

  if ( n > 0 ) then

    s1 = 1.0_QP + 2* m   * mu
    s2 = 1.0_QP + 2*(m+n)* mu
    theta = LOG(ABS(s1/s2))

    a = r
    b = -n*mu
    tt = a/b
    c = 1.0_QP + 2*m*mu
    w = c/b

#include "ft.h90"

    f = ff(0) + tt*(ff(1) + tt*(ff(2) + tt*(ff(3) + tt*(ff(4) + tt*(ff(5) + tt*(ff(6) + tt*(ff(7) + tt*ff(8))))))))
    f = f/b

  else

    a = r
    c = 1.0_QP + 2 * m * mu
    tt = a/c

    ff(0) = 2.0000000000000000000_QP
    ff(1) = 1.3333333333333333333_QP
    ff(2) = 1.0666666666666666667_QP
    ff(3) = 0.91428571428571428571_QP
    ff(4) = 0.81269841269841269841_QP
    ff(5) = 0.73881673881673881674_QP
    ff(6) = 0.68198468198468198468_QP
    ff(7) = 0.63651903651903651904_QP
    ff(8) = 0.59907674025321084145_QP

    f = ff(0) + tt*(ff(1) + tt*(ff(2) + tt*(ff(3) + tt*(ff(4) + tt*(ff(5) + tt*(ff(6) + tt*(ff(7) + tt*ff(8))))))))
    f = f/c

  endif

  return
end function

function integ_f_taylor2(n, m, r, mu) result(f)
!
! returns
!  F^n_m(r,mu)  with small r  in sm < r <= sp
!
  integer(8), intent(in) :: n,m
  real(QP), intent(in) :: r,mu
  complex(QP) :: f

  real(QP) :: a,c
  real(QP) :: ff(0:8),tt

  if ( n == 0 ) then
    a = r
    c = 1.0_QP + 2 * m * mu
    tt = a/c

    ff(0) = 2.0000000000000000000_QP
    ff(1) = 1.3333333333333333333_QP
    ff(2) = 1.0666666666666666667_QP
    ff(3) = 0.91428571428571428571_QP
    ff(4) = 0.81269841269841269841_QP
    ff(5) = 0.73881673881673881674_QP
    ff(6) = 0.68198468198468198468_QP
    ff(7) = 0.63651903651903651904_QP
    ff(8) = 0.59907674025321084145_QP

    f = ff(0) + tt*(ff(1) + tt*(ff(2) + tt*(ff(3) + tt*(ff(4) + tt*(ff(5) + tt*(ff(6) + tt*(ff(7) + tt*ff(8))))))))
    f = f/c

  endif

  return
end function

function integ_g(n, m, r, mu, f) result(g)
!
! returns
!
!  G^n_m(r,mu) =
!
! (i/mu) Integral_{v,-1,+1} Intelral_{t,0,Infinity} 
!            v * Exp[ -I ( phi + 2*m - n v + n ) t ]
!
!  = Integral_{v,-1,1} v/(r*v^2 -(n*mu)*v+1-r+(2*m+n)*mu)
!
!  phi = (1-(1-v^2)*r)/mu
!
!  r = k_perp^2/(4 me^2)
! mu = eB/me^2
!
  integer(8),  intent(in) :: n,m
  real(QP),    intent(in) :: r,mu
  complex(QP), intent(in) :: f
  complex(QP) :: g

  real(QP) :: s1,s2,sp,sm,tt,a,b
  real(QP) :: theta

  g = Q0

  if ( n == 0 ) then
    g = Q0
    return
  endif

  s1 = SQRT(1.0_QP+2* m   *mu)
  s2 = SQRT(1.0_QP+2*(m+n)*mu)
  sp = (s1+s2)**2/4.0_QP
!  sm = (s1-s2)**2/4.0_QP
  sm = ((n*mu)/(s1+s2))**2

  if ( r <= sm ) then
    if ( 0 < n ) then
      a = r
      b = -n*mu
      tt = r/b
      if ( ABS(tt) < DELTA ) then
        g = integ_g_taylor1(n,m,r,mu)
        return
      endif
    else
      g = Q0
      return
    endif
  endif

  s1 = 1.0_QP + 2* m   * mu
  s2 = 1.0_QP + 2*(m+n)* mu
  theta = LOG(ABS(s1/s2))

!write(*,'("S1,S2,THETA:",3ES24.15)')s1,s2,theta
!write(*,'("THETA,n*mu*f:",4ES24.15)')theta,n*mu*f

!write(*,'("GT:",2I4,8ES24.15)')n,m,r,mu,g
  g = ( theta + n * mu * f )/r/2.0_QP
!write(*,'(" G:",2I4,8ES24.15)')n,m,r,mu,g

  return
end function

function integ_g_taylor1(n, m, r, mu) result(g)
!
! returns
!  G^n_m(r,mu)  with small r in r <= sm
!
  integer(8),  intent(in) :: n,m
  real(QP),    intent(in) :: r,mu
  complex(QP) :: g

  complex(QP) :: ctmp
  real(QP) :: s1,s2,a,b,c,tt
  real(QP) :: theta,gg(0:8),w

  g = Q0

  if ( 0 < n ) then
    s1 = 1.0_QP + 2* m   * mu
    s2 = 1.0_QP + 2*(m+n)* mu
    theta = LOG(ABS(s1/s2))

    a = r
    b = -n*mu
    tt = a/b
    c = 1.0_QP + 2*m*mu
    w = c/b

#include "gt.h90"

    g = gg(0) + tt*(gg(1) + tt*(gg(2) + tt*(gg(3) + tt*(gg(4) + tt*(gg(5) + tt*(gg(6) + tt*(gg(7) + tt*gg(8))))))))
    g = g/b

  else

    g = Q0

  endif

  return
end function


function integ_h(n, m, r, mu, f) result(h)
!
! returns
!
!  H^n_m(r,mu) =
!
! (i/mu) Integral_{v,-1,+1} Intelral_{t,0,Infinity} 
!            v^2 * Exp[ -I ( phi + 2*m - n v + n ) t ]
!
!  = Integral_{v,-1,1} v^2/(r*v^2 -(n*mu)*v+1-r+(2*m+n)*mu)
!
!  phi = (1-(1-v^2)*r)/mu
!
!  r = k_perp^2/(4 me^2)
! mu = eB/me^2
!
  integer(8),  intent(in) :: n,m 
  real(QP),    intent(in) :: r,mu
  complex(QP), intent(in) :: f
  complex(QP) :: h

  real(QP) :: a,b,c,cp,d,s1,s2,sp,sm,theta
  real(QP) :: tt

  a = r
  b = -n*mu
  c = 1.0_QP - r + (2*m + n)*mu

  s1 = SQRT(1.0_QP+2* m   *mu)
  s2 = SQRT(1.0_QP+2*(m+n)*mu)
  sp = (s1+s2)**2/4.0_QP
!  sm = (s1-s2)**2/4.0_QP
  sm = ((n*mu)/(s1+s2))**2

!write(*,'("H: n,m:",2I9," r,mu,sm:",3ES24.15," s1,s2",2ES24.15)')n,m,r,mu,sm,s1,s2

  if ( r <= sm ) then
    if ( n > 0 ) then
      tt = a/b
      if ( ABS(tt) < DELTA ) then
        h =  integ_h_taylor1(n, m, r, mu)
        return
      endif
    else
      cp = 1.0_QP + 2*m*mu
      tt = a/c
      if ( ABS(tt) < DELTA ) then
        h =  integ_h_taylor1(n, m, r, mu)
        return
      endif
    endif
  else if ( sm < r .and. r <= sp ) then
    if ( n == 0 ) then
      cp = 1.0_QP + 2*m*mu
      tt = a/c
      if ( ABS(tt) < DELTA ) then
        h =  integ_h_taylor2(n, m, r, mu)
        return
      endif
    endif
  endif

  h = Q0
  s1 = 1.0_QP + 2* m   * mu
  s2 = 1.0_QP + 2*(m+n)* mu
  theta = LOG(ABS(s1/s2))
  h = (2.0_QP - b * theta/2.0_QP/r + (b*b-2*a*c)*f/2.0_QP/r ) /r

!write(*,'("H:",2I9," r,mu:",2ES24.15," sm:",ES24.15," h:",3ES24.15)')n,m,r,mu,sm,h

  return
end function

function integ_h_taylor1(n, m, r, mu) result(h)
!
! returns
!  H^n_m(r,mu)  with small r  in r <= sm
!
  integer(8),  intent(in) :: n,m 
  real(QP),    intent(in) :: r,mu
  complex(QP) :: h
  complex(QP) :: ctmp
  real(QP) :: a,b,c,d,s1,s2,theta
  real(QP) :: hh(0:8),tt,w

  h = Q0

  if ( 0 < n ) then

    s1 = 1.0_QP + 2* m   * mu
    s2 = 1.0_QP + 2*(m+n)* mu
    theta = LOG(ABS(s1/s2))

    a  = r
    b  = -n*mu
    tt = a/b
    c  = 1.0_QP + 2*m*mu
    w = c/b
#include "ht.h90"

    h = hh(0) + tt*(hh(1) + tt*(hh(2) + tt*(hh(3) + tt*(hh(4) + tt*(hh(5) + tt*(hh(6) + tt*(hh(7) + tt*hh(8))))))))
    h = h/b

  else

    a = r
    c = 1.0_QP + 2*m*mu
    tt = a/c

    hh(0) = 0.66666666666666666667_QP
    hh(1) = 0.26666666666666666667_QP
    hh(2) = 0.15238095238095238095_QP
    hh(3) = 0.10158730158730158730_QP
    hh(4) = 0.073881673881673881674_QP
    hh(5) = 0.056832056832056832057_QP
    hh(6) = 0.045465645465645465645_QP
    hh(7) = 0.037442296265825677590_QP
    hh(8) = 0.031530354750168991655_QP

    h = hh(0) + tt*(hh(1) + tt*(hh(2) + tt*(hh(3) + tt*(hh(4) + tt*(hh(5) + tt*(hh(6) + tt*(hh(7) + tt*hh(8))))))))
    h = h/c

  endif

  return
end function

function integ_h_taylor2(n, m, r, mu) result(h)
!
! returns
!  H^n_m(r,mu)  with small r  in sm < r <= sp
!
  integer(8),  intent(in) :: n,m 
  real(QP),    intent(in) :: r,mu
  complex(QP) :: h
  complex(QP) :: ctmp
  real(QP) :: a,b,c,d,s1,s2,theta
  real(QP) :: hh(0:8),tt,w

  h = Q0

  if ( 0 == n ) then

    a  = r
    c  = 1.0_QP + 2*m*mu
    tt = a/c

    hh(0) = 0.66666666666666666667_QP
    hh(1) = 0.26666666666666666667_QP
    hh(2) = 0.15238095238095238095_QP
    hh(3) = 0.10158730158730158730_QP
    hh(4) = 0.073881673881673881674_QP
    hh(5) = 0.056832056832056832057_QP
    hh(6) = 0.045465645465645465645_QP
    hh(7) = 0.037442296265825677590_QP
    hh(8) = 0.031530354750168991655_QP

    h = hh(0) + tt*(hh(1) + tt*(hh(2) + tt*(hh(3) + tt*(hh(4) + tt*(hh(5) + tt*(hh(6) + tt*(hh(7) + tt*hh(8))))))))
    h = h/c

  endif

  return
end function

function integ_i(n, m, r, mu, f) result(i)
!
! returns
!
!  I^n_m(r,mu) =
!
! (i/mu) Integral_{v,-1,+1} Intelral_{t,0,Infinity} 
!            v^3 * Exp[ -I ( phi + 2*m - n v + n ) t ]
!
!  = Integral_{v,-1,1} v^3/(r*v^2 -(n*mu)*v+1-r+(2*m+n)*mu)
!
!  phi = (1-(1-v^2)*r)/mu
!
!  r = k_perp^2/(4 me^2)
! mu = eB/me^2
!
  integer(8),  intent(in) :: n,m 
  real(QP),    intent(in) :: r,mu
  complex(QP), optional, intent(in) :: f
  complex(QP) :: i

  complex(QP) :: gg,ff
  real(QP) :: a,b,c

  a = r
  b = -n*mu
  c = 1.0_QP - r + (2*m+n)*mu
  if (present(f)) then
    ff = f
  else
    ff = integ_f(n,m,r,mu)
  endif
  gg = integ_g(n,m,r,mu,ff)

  i = ( (-2*b) + (b**2-c*a)*gg + (b*c)*ff )/(a**2)

  return
end function

end module
