program HItest
  use constant_mod
  use physics_constsnt_mod
  use photon_pol_HI_class
  use photon_pol_below_th_class
  use photon_pol_parallel_b_class
  use photon_pol_strong_b_class
  implicit none
  real(DP) :: r,q,mu,mass, eB
  complex(DP) :: f,rn_below_th,rn_para_b,rn_strong_b,df,f3
  complex(DP) :: rn0_below_th
  integer :: iout = 100
  integer :: i,j,k,itype,ith
  integer :: ptype, btype, mmax
  real(DP) :: dq,dr
  real(DP), parameter :: max_PT = 3000.0_DP
  real(DP), parameter :: max_PP = 4000.0_DP
  integer,      save :: NPLTR = 0
  integer,      save :: NPLTQ = 0
  character(len=256) :: arg,fname
  character(len=2) :: name(0:2) = [ "N0","N1","N2" ]
  character(len=5) :: suffix(0:1) = [ "below","above"]


  if  ( COMMAND_ARGUMENT_COUNT() /= 5 ) then
    write(*,'("Usage: HItest itype ith ptype btype mmax")')
    write(*,'("itype: tensor type takes a value from 0, 1, or 2. This means 0=N0, 1=N1, 2=N2.")')
    write(*,'("  ith: threshold flag takes a value 0 or 1. 0 is for below threshold, 1 for above threshold")')
    write(*,'("ptype: partilce mass type, 0 or 1. 0 is for muons, 1 for electrons")')
    write(*,'("btype: field strength type, 0,1,2. 0 is 10*m_pi^2, 1 for m_pi^2, 2 for m_pi^2/10")')
    write(*,'(" mmax: Landau level summation max.")')
    stop
  endif

  call GET_COMMAND_ARGUMENT(1,arg)
  read(arg,*)itype
  call GET_COMMAND_ARGUMENT(2,arg)
  read(arg,*)ith
  call GET_COMMAND_ARGUMENT(3,arg)
  read(arg,*)ptype
  call GET_COMMAND_ARGUMENT(4,arg)
  read(arg,*)btype
  call GET_COMMAND_ARGUMENT(5,arg)
  read(arg,*)mmax

  if ( itype < 0 .or. 2 < itype ) then
    write(*,'("itype should be 0 , 1, or 2. itype=",I3)')itype
    stop
  endif
  if ( ith   < 0 .or. 1 < ith   ) then
    write(*,'("ith should be 0 or 1. ith=",I3)')ith
    stop
  endif

  select case (ptype)
  case (0)
    mass = MASS_MUON
    write(*,'("# MUON:",ES24.15)')mass
  case (1)
    mass = MASS_ELECTRON
    write(*,'("# ELECTRON:",ES24.15)')mass
  case default
    write(*,'("particle ptype error. ptype=",I3)')ptype
    stop
  end select

  select case (btype)
  case (0)
    eB   = 10.0_DP*MASS_PION**2
  case (1)
    eB   =         MASS_PION**2
  case (2)
    eB   =  0.1_DP*MASS_PION**2
  case default
    write(*,'("field btype error. btype=",I3)')btype
    stop
  end select
  write(*,'("# eB:",ES24.15)')eB

  if (mmax <= 0) then
    write(*,'("mmax error. mmax=",I3)')mmax
    stop
  endif
  call new(mmax)
  write(*,'("# mmax:",I9)')mmax
  if (mmax < 1001) then
    NPLTR = 400
    NPLTQ = 30
  else
    NPLTR = 40
    NPLTQ = 30
  endif

  fname = TRIM(name(itype))//TRIM(suffix(ith))//".dat"
  open(iout,file=TRIM(fname),form="formatted")

  !===========================================
  ! test for muon in strong magnetic field
  ! Magnetic Field        = 10*(Mpi^2)
  ! Parallel invariant mass^2  = 0...(2*Mmu)^2 ,  (2*Mmu)^2 .. 1GeV^2
  ! Transverse momentum^2      = 0...(1000 MeV)^2
  !===========================================

  mu = eB/(mass**2)

   q = (max_PT)**2/(2*mass)**2
  dq = q/NPLTQ

  if (ith == 0) then      ! below threshold
    r = 1.0_DP
    dr = r/NPLTR
  elseif (ith == 1) then  ! above threshold
    r = (max_PP)**2/(2*mass)**2
    dr = (r-1.0_DP)/NPLTR
  endif


  write(iout,'("# mass =",E24.15)')mass
  write(iout,'("#   eB =",E24.15)')eB
  write(iout,'("#        mu         r          q           Nx(HI)         Nx(EXACT)         Nx(STB)        Nx(ParB)")')

  do i = 0,NPLTQ
!  do i = NPLTQ/2,NPLTQ
!  do i = NPLTQ,NPLTQ
!  do i = 3*NPLTQ/4,NPLTQ
    q = dq*i+0.0_DP

!-- test q = 0 start
!  if (i==1) stop
!-- test q = 0 end

  do j = 0,NPLTR
!  do j = 1,1
  
    if (ith == 0) then
      if ( j == 0 ) then
        r = 0.0_DP + EPSILON(1.0_DP)*10
      else if ( j == NPLTR ) then
        r = 1.0_DP - EPSILON(1.0_DP)*10
      else
        r = dr*j
      endif
    else if (ith == 1) then
      if ( j == 0 ) then
        r = 1.0_DP + EPSILON(1.0_DP)*10
      else
        r = 1.0_DP + dr*j
      endif
    endif

    !-----------------------------------
    ! Below Threshold
    ! Integral form
    !-----------------------------------
    if ( r < 1.0_DP ) then
      select case(itype)
      case(0)
        rn_below_th = N0_below_th(r=r,q=q,mu=mu) 
      case(1)
        rn_below_th = N1_below_th(r=r,q=q,mu=mu) 
      case(2)
        rn_below_th = N2_below_th(r=r,q=q,mu=mu)
      end select
    else
      rn_below_th = Z0
    endif

    !-----------------------------------
    ! Hattori-Itakura Landau level sum
    !-----------------------------------
    f3 = Q0
    select case(itype)
    case(0)
      f  = n0d_v3(r,q,mu)
    case(1)
      f  = n1d_v3(r,q,mu)
    case(2)
      f  = n2d_v3(r,q,mu)
    end select

    !-----------------------------------
    ! momuntm parallel to B
    ! DiGamma form
    !-----------------------------------
    if ( q < EPSILON(1.0_DP) ) then
      select case(itype)
      case(0)
        rn_para_b = N0_parallel_b(r=r,q=0.0_DP,mu=mu)
      case(1)
        rn_para_b = N1_parallel_b(r=r,q=0.0_DP,mu=mu)
      case(2)
        rn_para_b = N2_parallel_b(r=r,q=0.0_DP,mu=mu)
      end select
    else
      rn_para_b   = Z0
    endif

    !-----------------------------------
    ! Fukushima Lowest Landau Approx
    !-----------------------------------
    if (itype == 1) then
      rn_strong_b = N1_strong_b(r=r,q=q,mu=mu)
    else
      rn_strong_b = Z0
    endif

    write(*,'("mu=",E24.15," r=",E24.15," q=",E24.15," ",A,"(HI) =",2E26.15E4)')mu,r,q,name(itype),f
    write(*,'("mu=",E24.15," r=",E24.15," q=",E24.15," ",A,"(EX) =",2E26.15E4)')mu,r,q,name(itype),rn_below_th
    write(*,'("mu=",E24.15," r=",E24.15," q=",E24.15," ",A,"(PAR)=",2E26.15E4)')mu,r,q,name(itype),rn_para_b
    write(*,'("mu=",E24.15," r=",E24.15," q=",E24.15," ",A,"(STG)=",2E26.15E4)')mu,r,q,name(itype),rn_strong_b
    write(*,*)
    write(iout,'(3E24.15,8E26.15E4)')mu,r,q,f,rn_below_th,rn_strong_b,rn_para_b

  enddo ! end of do r
  write(iout,*)
  write(iout,*)
  enddo ! end of do q

  close(iout)

  stop

contains


!subroutine omg2test(rd,qd,mud)
!  real(DP), intent(in) :: rd,qd,mud
!  real(QP) :: r,q,mu
!  integer(8) :: n,m,k
!  complex(QP) :: ff,pf,pp
!
!  r  = 0.8_QP
!  mu = 100.0_QP
!
!  pf = Q0
!  n = 4
!  m = 2
!  do k=0,10000
!    pp = (integ_f(n-1,m+k+2,r,mu)+integ_f(n,m+k+1,r,mu))/(k+1)
!    pf = pf + pp
!    write(*,'(3I9,4E24.15)')n,m,k,pf,pp
!  enddo
!
!  stop
!end subroutine


end program
