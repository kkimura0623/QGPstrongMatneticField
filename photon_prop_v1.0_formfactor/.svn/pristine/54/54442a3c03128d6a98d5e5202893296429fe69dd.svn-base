program test_coef
  use physics_constsnt_mod
  use constant_mod
  use laguerre_pdf_class
!  use special_funcs_mod
  use pole_integ_func_class
  implicit none
  type(laguerre_pdf) :: coef
  integer(8) :: i,n,m,nmin,nmax,mmin,mmax,nstep
  real(DP) :: x,dx,c,dc,xmax,xmin
  real(DP) :: fx,r,mu
  complex(DP) :: ff,hh,gg
  integer(8) ::ierr,jj

goto 2000
  n=1
  m=9990
!  r = 2.987134656824239E+06_DP
!  r = 0.3446703E+06_DP
  r = 0.99_DP
!  r = 0.100_DP
  mu = MASS_PION**2/(MASS_ELECTRON)**2
!  mu= 0.10_DP

  ff = integ_f(n,m,r,mu)
  gg = integ_g(n,m,r,mu,ff)
  hh = integ_h(n,m,r,mu,ff)
  write(*,'(2I9,2ES24.15," F:",8ES24.15)')n,m,r,mu,ff
  write(*,'(2I9,2ES24.15," G:",8ES24.15)')n,m,r,mu,gg
  write(*,'(2I9,2ES24.15," H:",8ES24.15)')n,m,r,mu,hh
  r = 0.0_DP
  ff = integ_f(n,m,r,mu)
  gg = integ_g(n,m,r,mu,ff)
  hh = integ_h(n,m,r,mu,ff)
  write(*,'(2I9,2ES24.15," F:",8ES24.15)')n,m,r,mu,ff
  write(*,'(2I9,2ES24.15," G:",8ES24.15)')n,m,r,mu,gg
  write(*,'(2I9,2ES24.15," H:",8ES24.15)')n,m,r,mu,hh
stop

2000 continue
  nmin = 560
  nmax = 561
  nmin = 600
  nmax = 610
  mmin =    0
  mmax = 8000

  xmin  =  12.0_DP
  xmax  =  12.0_DP
!  xmin  = 0.128337841740507e2_DP
!  xmax  = 0.128337841740507e2_DP
!  xmin  = 0.00001_DP
!  xmax  = 0.00001_DP

  nstep = 0
  if ( nstep > 0 ) then
    dx = (xmax-xmin)/nstep
  else
    dx = 0.0_DP
  endif

  do i=0,nstep
   x = xmin+dx*i

   do n=nmin,nmax

     do m=mmax-9,mmax

!     do m=0,mmax
!       if (m==0) then
!         call new(coef,n,x)
!       endif
!       call get_next(coef,c,dc)

!       call laguerre_poly(m,n,x,c)
!       c = 0.0_DP
!       call laguerre_poly_asym(m,n,x,fx,ierr)
!       write(*,'(2I8,9ES26.15E4)')n,m,x,c,fx
!       if ( m > 5000 ) write(*,'(2I8,9ES26.15E4)')n,m,x,c,dc,fx

     enddo

   enddo
  enddo

  stop
end program
