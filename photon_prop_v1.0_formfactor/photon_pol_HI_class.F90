module photon_pol_HI_class
  use constant_mod
  use physics_constsnt_mod
  use polygamma_mod
  use laguerre_pdf_class
  use pole_integ_func_class
  implicit none
  private
  public :: set_limit_landau_level
  public :: new
  public :: n0d_v3
  public :: n1d_v3
  public :: n2d_v3

  integer, parameter :: NOMP=4

  integer(8), save :: m_nmax = 20000
  integer(8), save :: m_mmax =  4000

  real(DP), parameter :: TOL   = 1.0e-14_DP
  real(DP), parameter :: DTINY = 1.0e-305_DP
  integer,  parameter :: NWORK = 8000

  real(DP) :: DIGAMA
  external :: DIGAMA
  logical, save :: m_debug = .false.
!  logical, save :: m_debug = .true.


  type integ_paras
    real(DP)   :: mag_mu, mom_q, mom_r, feynman_v
    integer(8) :: landau_n, landau_m
    integer    :: Kceil
  end type

  type coef_c_f
    integer(8) :: n,m
    real(QP)   :: eta
    real(QP)   :: c(-1:0,-1:1)     ! c(n,m)
    real(QP)   :: d(-1:1)
    type(laguerre_pdf) :: sc(-1:1)
  end type

  type omg0_state
    type(coef_c_f) :: pol_c
  end type

  type omg1_state
    type(coef_c_f) :: pol_c
  end type

  type omg2_state
    complex(DP) :: prev_f
    logical :: prev_has_no_imag = .true.
    type(laguerre_pdf) :: dcoef
  end type

  interface set_limit_landau_level
    module procedure new_HI
  end interface

  interface new
    module procedure new_HI
  end interface

contains

subroutine new_HI(mmax,nmax)
  integer, intent(in) :: mmax
  integer, optional, intent(in) :: nmax
  m_mmax = mmax
  if (present(nmax)) then
    m_nmax = nmax
  else
    m_nmax = MIN(mmax*20,20000)
  endif
  return
end subroutine

subroutine set_param_HI(debug)
  logical, intent(in) :: debug
  m_debug = debug
  return
end subroutine

!======================================
! Laguerre distribution functions
! Coefficient functions C^n_m(eta)
!======================================
subroutine new_coef_c(this)
  type(coef_c_f), intent(inout) :: this
  this%n = -1
  this%m = -1
  this%eta = 0.0_QP 
  this%c(:,:) = 0.0_QP
  this%d(:)   = 0.0_QP
  return
end subroutine

function coef_c0(n,m) result(c0)
!
! Returns
!  c^n_m(eta = 0) = delta_{n,0}
!
  integer(8), intent(in)   :: n,m
  real(QP) :: c0
  if ( n == 0 ) then
    c0 = 1.0_QP
  else
    c0 = 0.0_QP
  endif
  return
end function
  
subroutine coef_c(this,n,m,eta,cc)
!
! Returns
!  c^n_m(eta) =  Exp(-eta) (m!/(m+n)!) * eta^n  [L^n_m(eta)]^2
!
! eta = 2 ( vec{k}_perp^2 / (2m)^2 ) /  (eB/m^2) = 2 q /mu
!
! cc( mshift, nshift)
! cc(-1,-1) = c^{n-1}_{m-1}
! cc( 0,-1) = c^{n-1}_{m}
! cc(-1, 0) = c^{n  }_{m-1}
! cc( 0, 0) = c^{n  }_{m}
! cc(-1, 1) = c^{n+1}_{m-1}
! cc( 0, 1) = c^{n+1}_{m}
!
  type(coef_c_f), intent(inout) :: this
  integer(8), intent(in)   :: n,m
  real(QP),   intent(in)   :: eta
  real(QP),   intent(inout):: cc(-1:0,-1:1) ! cc(m,n)
  real(QP) :: fx,c,  scfc,rtmp
  integer(8) :: i,nn,nd

  if (ABS(eta) < EPSILON(1.0_DP)) then
    cc(-1,-1) = coef_c0(n-1,m-1)
    cc( 0,-1) = coef_c0(n-1,m)
    cc(-1, 0) = coef_c0(n  ,m-1)
    cc( 0, 0) = coef_c0(n  ,m)
    cc(-1, 1) = coef_c0(n+1,m-1)
    cc( 0, 1) = coef_c0(n+1,m)
    return
  endif

  if ( n == 0 ) then
    cc(:,-1) = 0.0_QP
    nd = 0
  else
    nd = -1
  endif
  if ( m > 0 ) then
    cc(:,:) = this%c(:,:)
  endif

  do nn=nd,1

    if ( m == 0 ) then
      call new(this%sc(nn),n+nn,eta)
    endif
    call get_next(this%sc(nn),c)

    if ( m > 0 ) then
      cc(-1,nn) = cc( 0,nn)
    else
      cc(-1,nn) = 0.0_QP
    endif
    cc( 0,nn) = c

  enddo

  this%c(:,:) = cc(:,:)
!write(*,'("@",2I4,2E24.16)')n,m,cc(0,0),cc(-1,0)

  return
end subroutine

function coef_nc0(n,m) result(nc0)
!
! Returns
!  m*c^n(eta)/eta |_{eta=0}
!
  integer(8), intent(in)   :: n,m
  real(QP) :: nc0
  if ( n == 1 ) then
    nc0 = REAL(m + 1,kind=QP)
  else
    nc0 = 0.0_QP
  endif
  return
end function

subroutine coef_nc(n,m,eta,cc,nc)
!
! Returns
!  m*c^n(eta)/eta
!
  integer(8), intent(in) :: n,m
  real(QP), intent(in)  :: eta
  real(QP), intent(in)  :: cc(-1:0,-1:1)
  real(QP), intent(out) :: nc(-1:0)

  if (ABS(eta) < EPSILON(1.0_DP)) then
    nc(-1) = coef_nc0(n,m-1)
    nc( 0) = coef_nc0(n,m  )
    return
  endif

  if (n == 0) then
    nc(:) = 0.0_QP
  else
    nc(-1) = n*cc(-1,0)/eta
    nc( 0) = n*cc( 0,0)/eta
  endif

  return
end subroutine


!==========================================
! Reminder terms integrated
!==========================================

function n0_cterm_v3(q,mu) result(cterm)
!
! 2*int(v,0,1)int(x,0,infty)
! ( ((ch(v*x)-sh(v*x)*ch(x)/sh(x))/sh(x))*exp(-(2*q/mu)(ch(x)-ch(v*x))/sh(x)) - (1-v^2)/x )*exp( -x/mu) 
!
  use photon_pol_below_th_class
  real(DP), intent(in) :: q,mu
  real(DP) :: cterm

  cterm = N0_below_th(0.0_DP,q,mu)

  return
end function

function n1_cterm_v3(q,mu) result(cterm)
!
! 2*int(v,0,1)int(x,0,infty)
! (1-v^2)*( (ch(x)/sh(x))*exp(-(2*q/mu)(ch(x)-ch(v*x))/sh(x)) - (1/x) )*exp( -x/mu)
!
  use photon_pol_below_th_class
  real(DP), intent(in) :: q,mu
  real(DP) :: cterm

  cterm = N1_below_th(0.0_DP,q,mu)

  return
end function

function n2_cterm_v3(q,mu) result(cterm)
!
! 2*int(v,0,1)int(x,0,infty)
!   ( 2(ch(x)-ch(v*x)/(sh(x)^3) )*exp(-(2*q/mu)(ch(x)-ch(v*x))/sh(x))-(1-v^2)/x )*exp( -x/mu) 
!
  use photon_pol_below_th_class
  real(DP), intent(in) :: q,mu
  real(DP) :: cterm

  cterm = N2_below_th(0.0_DP,q,mu)

  return
end function

!===================================================
! Omega0 coefficients for N0 form factor
!===================================================

function omg0_v3(this,n,m,r,q,mu) result(c)
  type(omg0_state), intent(inout) :: this
  integer(8),  intent(in) :: n,m
  real(QP),    intent(in) :: r,q,mu
  complex(QP) :: c
  real(QP) :: cc(-1:0,-1:1),nc(-1:0)
  real(QP) :: eta,zero
  complex(QP) :: f,g,f0,g0

  zero = 0.0_QP
  c    = Q0
  eta  = 2*q/mu

  f  = integ_f(n,m,   r,mu)
  g  = integ_g(n,m,   r,mu,f)
  f0 = integ_f(n,m,zero,mu)
  g0 = integ_g(n,m,zero,mu,f0)

!write(*,'("n,m:",2I9," r,mu:",2ES24.15," fgf0g0:",8ES24.15)')n,m,r,mu,f,g,f0,g0

  if ( n == 0 .and. m == 0 ) call new_coef_c(this%pol_c)
  call coef_c(this%pol_c,n,m,eta,cc)
  call coef_nc(          n,m,eta,cc,nc)

  if ( n > 0 ) then

    if ( m > 0 ) then
      c = (cc( 0,-1)+ cc(-1,+1))*(f-f0) -(nc( 0)+ nc(-1))*(g-g0)
    else ! if ( m == 0 ) then
      c =  cc( 0,-1)            *(f-f0) - nc( 0)         *(g-g0)
    endif

  else ! if ( n == 0 ) then

    if ( m > 0 ) then
      c = 2* cc(-1,+1) * (f-f0) 
    else ! if ( m == 0 ) then
      c = Q0
    endif

  endif

  return
end function

function n0d_v3(rd,qd,mud) result(n0d)
  real(DP), intent(in) :: rd,qd,mud
  complex(DP) :: n0d
  integer(8) :: n,m, kn, komp
  real(QP) :: r,q,mu
  complex(QP) :: n0
  complex(QP) :: o0(0:NOMP-1),dc,err(0:NOMP-1)
  real(DP) :: cterm,err_r,err_i
  type(omg0_state) :: omg0s

  r  = rd
  q  = qd
  mu = mud

  if (m_debug) write(*,'("n0 at mu=",ES24.15," r=",ES24.15," q=",ES24.15," (mmax,nmax)=",2I6)')mu,r,q,m_mmax,m_nmax
 ! write(*,'("(mmax,nmax)=",2I6)')m_mmax,m_nmax
  err_r = 0.0_QP
  err_i = 0.0_QP
  n0 = Q0
  loop_n: do kn = 0,m_nmax/NOMP

!$OMP PARALLEL DO PRIVATE(komp,n,m,dc,omg0s)
    do komp=0,NOMP-1
      n = komp + kn*NOMP
      o0(komp) = Q0
      do m = 0,m_mmax
        dc = omg0_v3(omg0s,n,m,r,q,mu)
        o0(komp) = o0(komp) + dc
        if ( m == m_mmax ) then
          err(komp) = dc
        endif
      enddo
    enddo

    do komp=0,NOMP-1
      n = komp + kn*NOMP
      if ( n > 0 ) then
        if (m_debug) write(*,'("n0_",I8,"=",8E24.15)')n,o0(komp)*2,err(komp)*2
        err_r = err_r +  REAL(2*err(komp),kind=QP)**2
        err_i = err_i + AIMAG(2*err(komp))**2
        n0 = n0 + o0(komp)*2
      else
        if (m_debug) write(*,'("n0_",I8,"=",8E24.15)')n,o0(komp),err(komp)
        err_r = err_r +  REAL(err(komp),kind=QP)**2
        err_i = err_i + AIMAG(err(komp))**2
        n0 = n0 + o0(komp)
      endif
      if ( n > 2 .and. ABS(n0) < EPSILON(1.0_DP) .and. ABS(REAL(2*o0(komp),kind=QP)) < EPSILON(1.0_DP) ) then
        if (m_debug) write(*,'("n0 end at n=",I9," n0=",2ES24.15," dn0=",2ES24.15)')n,n0,2*o0(komp)
        exit loop_n
      endif
      if ( n > 2 .and. ABS(REAL(2*o0(komp),kind=QP)/REAL(n0,kind=QP)) < 10*EPSILON(1.0_DP) ) then
        if (m_debug) write(*,'("n0 end at n=",I9," n0=",2ES24.15," dn0=",2ES24.15)')n,n0,2*o0(komp)
        exit loop_n
      endif
      if ( n == m_nmax ) then
        write(*,'("@ Landau level summation does not converge. STOP !")')
        stop
      endif
    enddo

  enddo loop_n

  cterm =  n0_cterm_v3(qd,mud)
!  write(*,'(8E24.15)') REAL(mu*n0,kind=DP),cterm,REAL(mu*n0,kind=DP)/(PI*4),cterm/(PI*4)

  n0 = -mu*n0/(4.0_QP*PIQ) + cterm
  n0d = n0

  err_r = mu*SQRT(err_r)/ 4.0_QP/PIQ
  err_i = mu*SQRT(err_i)/ 4.0_QP/PIQ

  if (m_debug) then
    write(*,'("n0r,err_r:",8F24.15)')  REAL(n0,kind=DP),err_r
    write(*,'("n0i,err_i:",8F24.15)') AIMAG(n0        ),err_i
  endif

  return
end function


function omg1_v3(this,n,m,r,q,mu) result(c)
  type(omg1_state), intent(inout) :: this
  integer(8), intent(in) :: n,m
  real(QP),   intent(in) :: r,q,mu
  real(QP) :: cc(-1:0,-1:1),eta
  complex(QP) :: c
  complex(QP) :: c0(0:2)
  complex(QP) :: fh,f,h,f0,h0
  real(QP) ::  cc_org

  eta = 2*q/mu

  if ( n == 0 .and. m == 0 ) call new_coef_c(this%pol_c)
  call coef_c(this%pol_c,n,m,eta,cc)

  f  = integ_f(n,m,r,mu)
  h  = integ_h(n,m,r,mu,f)
  f0 = integ_f(n,m,r=0.0_QP,mu=mu)
  h0 = integ_h(n,m,r=0.0_QP,mu=mu,f=f0)
  fh = (f - f0) - (h - h0)

  c = Q0
  if ( m > 0 ) then
    c = (cc( 0,0) + cc(-1,0))* fh   ! (c(n,m) + c(n,m-1))*(F(n,m)-H(n,m))
  else
    c =  cc( 0,0) * fh              !  c(n,m)*(F(n,m)-H(n,m))
  endif

  return
end function


function n1d_v3(rd,qd,mud) result(n1d)
  real(DP), intent(in) :: rd,qd,mud
  complex(DP) :: n1d
  integer(8) :: n,m, kn, komp
  real(QP) :: r,q,mu
  complex(QP) :: n1
  complex(QP) :: o0(0:NOMP-1),dc,err(0:NOMP-1)
  real(DP) :: cterm,err_r,err_i
  type(omg1_state) :: omg1s

  r  = rd
  q  = qd
  mu = mud

  if (m_debug) write(*,'("n1 at mu=",ES24.15," r=",ES24.15," q=",ES24.15," (mmax,nmax)=",2I6)')mu,r,q,m_mmax,m_nmax

  err_r = 0.0_QP
  err_i = 0.0_QP
  n1 = Q0
  loop_n: do kn = 0,m_nmax/NOMP

!$OMP PARALLEL DO PRIVATE(komp,n,m,dc,omg1s)
    do komp=0,NOMP-1
      n = komp + kn*NOMP
      o0(komp) = Q0
      do m = 0,m_mmax
        dc = omg1_v3(omg1s,n,m,r,q,mu)
        o0(komp) = o0(komp) + dc
        if ( m == m_mmax ) then
          err(komp) = dc
        endif
      enddo
    enddo

    do komp=0,NOMP-1
      n = komp + kn*NOMP
      if ( n > 0 ) then
        if (m_debug) write(*,'("n1_",I8,"=",8E24.15)')n,o0(komp)*2,err(komp)*2
        err_r = err_r +  REAL(2*err(komp),kind=QP)**2
        err_i = err_i + AIMAG(2*err(komp))**2
        n1 = n1 + o0(komp)*2
      else if ( n == 0 ) then
        if (m_debug) write(*,'("n1_",I8,"=",8E24.15)')n,o0(komp),err(komp)
        err_r = err_r +  REAL(err(komp),kind=QP)**2
        err_i = err_i + AIMAG(err(komp))**2
        n1 = n1 + o0(komp)
      endif
      if ( n > 2 .and. ABS(n1) < EPSILON(1.0_DP) .and. ABS(REAL(2*o0(komp),kind=QP)) < EPSILON(1.0_DP) ) then
        if (m_debug) write(*,'("n1 end at n=",I9," n1=",2ES24.15," dn1=",2ES24.15)')n,n1,2*o0(komp)
        exit loop_n
      endif
      if ( n > 2 .and. ABS(REAL(2*o0(komp),kind=QP)/REAL(n1,kind=QP)) < 10*EPSILON(1.0_DP) ) then
        if (m_debug) write(*,'("n1 end at n=",I9," n1=",2ES24.15," dn1=",2ES24.15)')n,n1,2*o0(komp)
        exit loop_n
      endif
      if ( n == m_nmax ) then
        write(*,'("@ Landau level summation does not converge. STOP !")')
        write(*,'("n1 at n=",I9," n1=",2ES24.15," dn1=",2ES24.15)')n,n1,2*o0(komp)
        stop
      endif
    enddo

  enddo loop_n

  cterm = n1_cterm_v3(qd,mud)
!  n1 = -(mu*n1 + cterm) / 4.0_QP/PIQ
  n1 = -mu*n1/(4.0_QP*PIQ) + cterm
  n1d = n1

  err_r = mu*SQRT(err_r)/ 4.0_QP/PIQ
  err_i = mu*SQRT(err_i)/ 4.0_QP/PIQ

  if (m_debug) then
    write(*,'("n1r,err_r:",8F24.15)')  REAL(n1,kind=DP),err_r
    write(*,'("n1i,err_i:",8F24.15)') AIMAG(n1        ),err_i
  endif

  return
end function

!============================================
! Functions for N2 form factor
!============================================
function Aphase(v,r,mu,n,m) result(A)
  real(DP), intent(in)  :: v,r,mu
  integer(8),  intent(in) :: n,m
  real(DP) :: A
  A = (r*v*v - n*mu*v + 1.0_DP-r+(2*m+n+2)*mu)/(2*mu)
  return
end function

function Aphase_min(r,mu,n,m) result(A)
  real(DP),   intent(in) :: r,mu
  integer(8), intent(in) :: n,m
  real(DP) :: A
  A = ( 1.0_DP-r+(2*m+n+2)*mu - (n*mu)*(n*mu)/4.0_DP/r )/(2*mu)
  return
end function


function Rfactor_B(v, c_paras) result(f)
  real(DP),         intent(in) :: v
  character(len=1), intent(in) :: c_paras(:) ! void data to parameters
  type(integ_paras) :: paras
  real(DP) :: f
  real(DP) :: A,B
  real(DP) :: dB,ft,fg,dff
  integer  :: ierr

  !
  ! restore input prameters from void data
  !
  paras = TRANSFER(c_paras,paras)

  A = Aphase(v,paras%mom_r ,paras%mag_mu,paras%landau_n,paras%landau_m+paras%Kceil+1)
  B = Aphase(v,0.0_DP      ,paras%mag_mu,paras%landau_n,paras%landau_m+paras%Kceil+1)

  dB = paras%mom_r*(v*v-1.0_DP)/(2*paras%mag_mu)
  if ( ABS(dB/B) < 0.01_DP ) then
  !===============================================================
  ! 8th order Taylor expansion
  !  DiGamma(B+dB) - DiGamma(B)
  !    = PolyGamma(1,B)*dB + PolyGamma(2,B)*dB^2/2 + ....
  !===============================================================
    f = polygamma(8,B)
    f = polygamma(7,B) + dB*f/8.0_DP
    f = polygamma(6,B) + dB*f/7.0_DP
    f = polygamma(5,B) + dB*f/6.0_DP
    f = polygamma(4,B) + dB*f/5.0_DP
    f = polygamma(3,B) + dB*f/4.0_DP
    f = polygamma(2,B) + dB*f/3.0_DP
    f = polygamma(1,B) + dB*f/2.0_DP
    f = dB*f
!  write(*,'(2I9,29ES24.15)')paras%landau_n,paras%landau_m,paras%mag_mu,paras%mom_r,dA,A-B,A,B,f,digamma(A)-digamma(B)
  else
    f = digamma(A) - digamma(B)
  endif

!     A = 1.5_DP
!     B = 1.5_DP + 1.5_DP*0.01_DP
!    dB =        - 1.5_DP*0.01_DP
!    ft = polygamma(8,B)
!    ft = polygamma(7,B) + dB*ft/8.0_DP
!    ft = polygamma(6,B) + dB*ft/7.0_DP
!    ft = polygamma(5,B) + dB*ft/6.0_DP
!    ft = polygamma(4,B) + dB*ft/5.0_DP
!    ft = polygamma(3,B) + dB*ft/4.0_DP
!    ft = polygamma(2,B) + dB*ft/3.0_DP
!    ft = polygamma(1,B) + dB*ft/2.0_DP
!    ft = dB*ft
!    fg = digamma(A)-digamma(B)
!    dff = ft-fg
!write(*,'(9ES24.15)')A,B,dB,dB/B,ft,fg,dff

!  f = DIGAMA(A,ierr)
!  if ( ierr /= 0 ) then
!    write(*,'("DIGAMMA :ierr=",I3)')ierr
!  endif
!
!  f = f - DIGAMA(B,ierr)
!  if ( ierr /= 0 ) then
!    write(*,'("DIGAMMA :ierr=",I3)')ierr
!  endif

  return
end function

function Rfactor(this,r,mu,n,m) result(f)
  use intde2_mod
  type(omg2_state), intent(inout) :: this
  real(DP),         intent(in) :: r,mu
  integer(8),       intent(in) :: n,m
  character(len=1), allocatable :: c_paras(:)
  type(integ_paras) :: paras
  integer :: length
  complex(DP) :: f,ff,f0
  complex(DP) :: frecl,df,fp0,fp1
  real(DP) :: a,b,result,zero
  real(DP) :: aw(0:NWORK-1),err
  real(DP) :: Amin,vmin
  real(DP) :: rtmp(10),est_err
  logical :: is_recl,is_err
  integer :: j
  integer(8) :: jj

  if ( m == 0 ) then
    this%prev_f = Q0
    this%prev_has_no_imag = .TRUE.
  endif

  paras%landau_n = n
  paras%landau_m = m
  paras%mom_r    = r
  paras%mag_mu   = mu

  if ( ABS(r/mu) < 10*EPSILON(1.0_DP) ) then
    f = Z0
    return
  endif

  !---------------------------------------------------------
  ! Threshold check in the integration interval 'v'
  !---------------------------------------------------------
  vmin = (n*mu)/(2*r)
  Amin = Aphase_min(r ,mu, n, m)
  if ( ABS(vmin) <= 1.0_DP ) then
    if ( Amin <= 1.0_DP ) then
      paras%Kceil = -CEILING(Amin)+1  ! Kceil(>=1) is 1,2,3,..., a positive number
    else
      paras%Kceil = -1
    endif
  else
    paras%Kceil = -1
  endif

  !---------------------------------------------------------
  ! If the previsous step is computed without imaginary part,
  ! current step can be recursively computed via
  ! 
  ! R^n_m = R^n_{m-1} + 2*mu*(F^n_{m}-F0^n_{m})
  !
  !---------------------------------------------------------
#define _RECL
#ifdef _RECL
  is_recl = .FALSE.
  frecl = Z0
  if (m >= 1            .and.  &
 &    mod(m,100)  /= 0  .and.  &
 &    paras%Kceil == -1 .and.  &
 &    ABS(REAL(this%prev_f,kind=DP)) > EPSILON(1.0_DP)*10 .and.  &
 &    this%prev_has_no_imag ) then
    zero = 0.0_QP
    ff = integ_f(n,m,r   ,mu)
    f0 = integ_f(n,m,zero,mu)
    f  = this%prev_f + 2*mu*(ff-f0)
    fp0 = this%prev_f
    fp1 = 2*mu*(ff-f0)
    this%prev_f = f
    is_recl = .TRUE.
    return
  endif
#endif

  call intdeini(NWORK,DTINY,TOL,aw)

  !------------------------------
  ! set void data to the paras
  !------------------------------
  length = SIZE(TRANSFER(paras,c_paras))
  allocate(c_paras(length))
  c_paras = TRANSFER(paras,c_paras)

  a =-1.0_DP
  b = 1.0_DP
  call intde(Rfactor_B,c_paras,a,b,aw,result,err)

  is_err = ( err < 0.0_DP )

  if ( is_err ) then

    write(*,'(" INTDE : ERROR STOP :",ES24.15," Kceil=",I3," Amin=",ES24.15," vmin=",ES24.15)')err,paras%Kceil,Amin,vmin
    write(*,'(" n,m,r,mu :",2I6,2ES24.15," f=",ES24.15)')n,m,r,mu,result
    stop

  endif
  f = CMPLX(result,0.0_DP,kind=DP)

  !----------------------------------
  ! compute contrubitions contain
  ! Imaginary part above threshold
  !----------------------------------
  if (paras%Kceil >= 0) then
    zero = 0.0_QP
    do j=0,paras%Kceil
      jj = m+1_8+j
      ff = integ_f(n,jj,r   ,mu)
      f0 = integ_f(n,jj,zero,mu)
      f  = f - 2*mu*(ff-f0)
    enddo
  endif

  !----------------------------------
  ! keep current step result R^n_m
  !----------------------------------
  if (paras%Kceil == -1) then
    this%prev_has_no_imag = .TRUE.
  else
    this%prev_has_no_imag = .FALSE.
  endif
  this%prev_f = f

!  endif
!df = f-frecl
!if (is_recl .AND. ABS(df) > EPSILON(1.0_DP)*100) then
!  write(*,'("DIFF:",2I6," r,mu:",2ES24.15," df,f,fr:",10ES24.15)')n,m,r,mu,f-frecl,f,frecl,fp0,fp1
!endif

  return
end function

function omg2_v3(this,n,m,r,q,mu) result(c)
  type(omg2_state), intent(inout) :: this
  integer(8),       intent(in)    :: n,m
  real(QP),         intent(in)    :: r,q,mu
  complex(QP) :: c,RF
  real(QP)    :: cc,dcc,eta,d(1:2),d0(1:2),zero


  eta = 2*q/mu

  if ( 0 == m ) call new(this%dcoef,n,eta)

  call get_next(this%dcoef,cc,dcc)

  RF = Rfactor(this,r,mu,n,m)
  c = 4*dcc*RF

  return
end function

function n2d_v3(rd,qd,mud) result(n2d)
  real(DP), intent(in) :: rd,qd,mud
  complex(DP) :: n2d
  integer(8) :: n,m, kn, komp
  real(QP) :: r,q,mu
  complex(QP) :: n2
  complex(QP) :: o0(0:NOMP-1),dc,err(0:NOMP-1)
  real(DP) :: cterm,err_r,err_i
  type(omg2_state) :: omg2s

  r  = rd
  q  = qd
  mu = mud

  if (m_debug) write(*,'("n2 at mu=",ES24.15," r=",ES24.15," q=",ES24.15," (mmax,nmax)=",2I6)')mu,r,q,m_mmax,m_nmax

  err_r = 0.0_QP
  err_i = 0.0_QP
  n2 = Q0
  loop_n: do kn = 0,m_nmax/NOMP

!$OMP PARALLEL DO PRIVATE(komp,n,m,dc,omg2s)
    do komp=0,NOMP-1
      n = komp + kn*NOMP
      o0(komp) = Q0
      do m = 0,m_mmax
        dc = omg2_v3(omg2s,n,m,r,q,mu)
        o0(komp) = o0(komp) + dc
        if ( m == m_mmax ) then
          err(komp) = dc
        endif
      enddo
    enddo

    do komp=0,NOMP-1
      n = komp + kn*NOMP
      if ( n > 0 ) then
        if (m_debug) write(*,'("n2_",I8,"=",8E24.15)')n,o0(komp)*2,err(komp)*2
        err_r = err_r +  REAL(2*err(komp),kind=QP)**2
        err_i = err_i + AIMAG(2*err(komp))**2
        n2 = n2 + o0(komp)*2
      else if ( n == 0 ) then
        if (m_debug) write(*,'("n2_",I8,"=",8E24.15)')n,o0(komp),err(komp)
        err_r = err_r +  REAL(err(komp),kind=QP)**2
        err_i = err_i + AIMAG(err(komp))**2
        n2 = n2 + o0(komp)
      endif
      if ( n > 2 .and. ABS(n2) < EPSILON(1.0_DP) .and. ABS(REAL(2*o0(komp),kind=QP)) < EPSILON(1.0_DP) ) then
        if (m_debug) write(*,'("n2 end at n=",I9," n2=",2ES24.15," dn2=",2ES24.15)')n,n2,2*o0(komp)
        exit loop_n
      endif
      if ( n > 2 .and. ABS(REAL(2*o0(komp),kind=QP)/REAL(n2,kind=QP)) < 10*EPSILON(1.0_DP) ) then
        if (m_debug) write(*,'("n2 end at n=",I9," n2=",2ES24.15," dn2=",2ES24.15)')n,n2,2*o0(komp)
        exit loop_n
      endif
      if ( n == m_nmax ) then
        write(*,'("@ Landau level summation does not converge. STOP !")')
        stop
      endif
    enddo

  enddo loop_n

  cterm = n2_cterm_v3(qd,mud)
! write(*,'(3E24.15,2F24.15)')mu,r,q,-REAL(n2)/(4*PIQ),-cterm/(4*PIQ)
  n2 = -n2/(4.0_QP*PIQ) + cterm
  n2d = n2

  err_r = SQRT(err_r)/ 4.0_QP/PIQ
  err_i = SQRT(err_i)/ 4.0_QP/PIQ

  if (m_debug) then
    write(*,'("n2r,err_r:",8F24.15)')  REAL(n2,kind=DP),err_r
    write(*,'("n2i,err_i:",8F24.15)') AIMAG(n2        ),err_i
  endif
!stop

  return
end function

end module
