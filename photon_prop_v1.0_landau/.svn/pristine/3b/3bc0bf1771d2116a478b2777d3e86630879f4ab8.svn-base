!================================================
! Probablity distribution functions
! associated with Laguerre polynomials
!================================================
module laguerre_pdf_class
!  use constant_mod
  implicit none
  private

  public :: laguerre_pdf
  public :: new
  public :: get_next

  integer, parameter :: DP=KIND(1.0d0)
!  integer, parameter :: QQP=KIND(1.0d0)
  integer, parameter :: QQP=KIND(1.0q0)
  integer, parameter :: MAX_NUM_SCALE=100

  type laguerre_pdf
    integer(8) :: n = 0
    integer(8) :: m = 0
    real(QQP) ::   x = 0.0_QQP
    real(QQP) ::  f0 = 0.0_QQP
    real(QQP) ::  f1 = 0.0_QQP
    real(QQP) ::  f2 = 0.0_QQP
    real(QQP) :: df0 = 0.0_QQP
    real(QQP) :: df1 = 0.0_QQP
    real(QQP) :: df2 = 0.0_QQP
    real(QQP) ::   c = 0.0_QQP
    real(QQP) ::  dc = 0.0_QQP
    real(QQP) :: scale(MAX_NUM_SCALE)
    real(QQP) :: prefacA,prefacB
    real(QQP) :: logfacA,logfacB
    integer   :: num_scale = 0
    logical   :: is_initialized = .false.
    logical   :: is_x_larger_than_one
  end type

  interface new
    module procedure new_laguerre_pdf
  end interface

  interface get_next
    module procedure get_next_laguerre_pdf
  end interface

contains

subroutine new_laguerre_pdf(this,alp,x)
  type(laguerre_pdf), intent(inout):: this
  integer(8), intent(in) :: alp
  real(DP),   intent(in) :: x
  real(QQP) :: logfac
  integer(8) :: k
  this%n = alp
  this%m = -1_8
  this%x = x
  this%f0  = 0.0_QQP
  this%f1  = 0.0_QQP
  this%f2  = 0.0_QQP
  this%df0 = 0.0_QQP
  this%df1 = 0.0_QQP
  this%df2 = 0.0_QQP
  this%c   = 0.0_QQP
  this%dc  = 0.0_QQP
  this%scale(:)  = 0.0_QQP
  this%num_scale = 0

  this%is_x_larger_than_one = ( this%x > 1.0_QQP ) 

  !-------------------------------------------------
  ! compute scale factor
  ! 
  ! x > 1
  !  logfacA = x - sum_{k=1,n}log(k) +  n   *log(x)
  !  logfacB = x - sum_{k=1,n}log(k) + (n-1)*log(x)
  !
  ! for exp(-logfacA/2) and exp(-logfacB/2)
  !
  ! x <= 1
  !  logfacA = x - sum_{k=1,n}log(k)
  !  logfacB = x - sum_{k=1,n}log(k)
  !
  ! for exp(-logfacA/2)*x^(n/2)
  ! and exp(-logfacB/2)*x^((n-1)/2)
  !
  !-------------------------------------------------
  select case (this%n)
  case(0)
    this%logfacA = this%x
    this%logfacB = 0.0_QQP
  case(1)
    if (this%is_x_larger_than_one) then
      this%logfacA = this%x -  this%n   *LOG(this%x)
      this%logfacB = this%x
    else
      this%logfacA = this%x
      this%logfacB = this%x
    endif
  case default
    logfac = this%x
    do k=1,this%n
      logfac = logfac + LOG(REAL(k,kind=QQP))
    enddo
    if (this%is_x_larger_than_one) then
      this%logfacA = logfac - this%n   *LOG(this%x)
      this%logfacB = logfac -(this%n-1)*LOG(this%x)
    else
      this%logfacA = logfac
      this%logfacB = logfac
    endif
  end select

  this%is_initialized = .TRUE.

  return
end subroutine

subroutine get_next_laguerre_pdf(this,c,dc)
!  use IEEE_ARITHMETIC
  type(laguerre_pdf), intent(inout):: this
  real(DP),           intent(out) :: c
  real(DP), optional, intent(out) :: dc
  real(QQP) :: coef1,coef2,coef3
  real(QQP) :: efac
  real(QQP) :: f2f2,f2df2
  real(QQP) :: logfac
  real(QQP), parameter :: MAXNORM=1.0e+100_QQP
  real(QQP), parameter :: MINNORM=1.0e-100_QQP
  integer(8) :: j
  logical :: is_underflow, is_overflow

  if (.not. this%is_initialized ) then
    stop
  endif

  this%m = this%m + 1_8

  select case(this%m)
  case(0)

    this%f0  = 0.0_QQP
    this%f1  = 0.0_QQP
    this%f2  = 1.0_QQP

    this%df0 = 0.0_QQP
    this%df1 = 0.0_QQP
    this%df2 = 0.0_QQP

  case(1)

    coef1 = (REAL(this%n+1_8,kind=QQP) - this%x)
    coef3 = 1.0_QQP/SQRT(REAL(1_8+this%n,kind=QQP))

    this%f0  = this%f1
    this%f1  = this%f2
    this%f2  =( coef1*this%f1            ) * coef3

    this%df0 = this%df1
    this%df1 = this%df2
    this%df2 =( coef1*this%df1 - this%f1 ) * coef3

  case default

    coef1 = (REAL(2_8*this%m+this%n-1_8,kind=QQP) - this%x)
    coef2 =        -SQRT(REAL(this%m+this%n-1_8,kind=QQP)*REAL(this%m-1_8,kind=QQP))
    coef3 = 1.0_QQP/SQRT(REAL(this%m+this%n    ,kind=QQP)*REAL(this%m    ,kind=QQP))

    this%f0  = this%f1
    this%f1  = this%f2
    this%f2  = ( coef1*this%f1  + coef2*this%f0            ) * coef3

    this%df0 = this%df1
    this%df1 = this%df2
    this%df2 = ( coef1*this%df1 + coef2*this%df0 - this%f1 ) * coef3

    is_overflow  = (ABS(this%f2) > MAXNORM .or. ABS(this%df2) > MAXNORM )
    is_underflow = (ABS(this%f2) < MINNORM .or. ABS(this%df2) < MINNORM )

    if ( is_overflow .or. is_underflow ) then

      this%num_scale = this%num_scale + 1
      if (this%num_scale > MAX_NUM_SCALE) then
        write(*,'("laguerre_pdf_class: num_scale exceeds MAX_NUM_SCALE. STOP.")')
        stop
      endif
      if (is_overflow)  this%scale(this%num_scale) = MAX(ABS(this%f2),ABS(this%df2))
      if (is_underflow) this%scale(this%num_scale) = MIN(ABS(this%f2),ABS(this%df2))
      efac = 1.0_QQP/this%scale(this%num_scale)
      this%f0  =  this%f0 * efac
      this%f1  =  this%f1 * efac
      this%f2  =  this%f2 * efac
      this%df0 = this%df0 * efac
      this%df1 = this%df1 * efac
      this%df2 = this%df2 * efac

    endif
 
  end select

  !
  ! efac = exp(-x/2) x^(n/2) Prod[Scale(k),k=1,Nscale]
  !
  logfac = this%logfacA
  do j=1,this%num_scale
    logfac = logfac - 2*LOG(this%scale(j))
  enddo
  if (this%is_x_larger_than_one) then
    efac = EXP(-logfac*0.5_QQP)
  else
    efac = EXP(-logfac*0.5_QQP)*SQRT(this%x)**(this%n)
  endif

  f2f2   = (this%f2*efac)*(this%f2*efac)
  this%c = f2f2
  c  = this%c

  if (present(dc)) then

    f2df2 = (this%f2*efac)*(this%df2*efac)

    select case(this%n)
    case(0)

      this%dc = (2*f2df2 - f2f2)

    case(1)

      !
      ! efac = exp(-x/2)/Sqrt(n!)*Prod[Scale(k),k=1,Nscale]
      !
      logfac = this%logfacB
      do j=1,this%num_scale
        logfac = logfac - 2*LOG(this%scale(j))
      enddo
      efac = EXP(-logfac*0.5_QQP)

      this%dc = (2*f2df2 - f2f2) +        (this%f2*efac)*(this%f2*efac)

    case default

      !
      ! efac = exp(-x/2) x^((n-1)/2) Prod[Scale(k),k=1,Nscale]
      !
      logfac = this%logfacB
      do j=1,this%num_scale
        logfac = logfac - 2*LOG(this%scale(j))
      enddo
      if (this%is_x_larger_than_one) then
        efac = EXP(-logfac*0.5_QQP)
      else
        efac = EXP(-logfac*0.5_QQP)*SQRT(this%x)**(this%n-1)
      endif

      this%dc = (2*f2df2 - f2f2) + this%n*(this%f2*efac)*(this%f2*efac)

    end select

    dc = this%dc
  endif

!write(*,'("n,m,x,c,dc,scale:",2I9,10ES24.15)')this%n,this%m,this%x,this%c,this%dc,this%scale(1),this%scale(2),this%scale(3)

  return
end subroutine

end module
