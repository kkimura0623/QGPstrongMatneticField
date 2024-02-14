program epem_pair_prod_test
!
! Test program for e+,e- pair production code
!
! Magnetic field eB on z-axis
! compare the virtual photon 3-momentum
! in parallel to eB and in perpendicular to eB
!
! Invariant mass and Spatial moometum are scaned.
!
  use constant_mod
  use physics_constsnt_mod
  use photon_prop_class
  use current_class
  implicit none
  real(DP) :: Lep(0:NDIM-1,0:NDIM-1)
  real(DP) :: Src(0:NDIM-1,0:NDIM-1)
  real(DP) :: p1(0:NDIM-1),p2(0:NDIM-1)
  complex(DP) :: Prop(0:NDIM-1,0:NDIM-1)
  complex(DP) :: PP(0:NDIM-1,0:NDIM-1)
  complex(DP) :: A(NDIM,NDIM)
  complex(DP) :: sbrn1
  complex(DP) :: Rate
  complex(DP) :: n0, n1, n2
  real(DP) :: qmom(0:NDIM-1)
  real(DP) :: eB,r
  real(DP) :: mass2,tol
  real(DP) :: mass2_0,mass2_1,dmass2
  real(DP) :: px,py,pz,pT, pT0,pT1,dpT
  integer :: ipx,ipy,ipz,ipt,NPT,imass2,NMASS2
  integer :: i,j,itype,status(2)
  integer :: mu,nu,ia,ib
  integer :: iout, pol_type, particle_type
  character(len=256) :: fname

  tol = 1.0e-8_DP
  eB = 10*MASS_PION**2
!  pol_type = POLTYPE_ITAKURA_FULL ! (Full formula by H.-I.)
!  pol_type = POLTYPE_ITAKURA_N1N0 ! (N1N0 approx)
!  pol_type = POLTYPE_ITAKURA_N1 ! (N1 approx)
  pol_type = POLTYPE_FUKUSHIMA ! (Lowest Landau level approx)
  particle_type = PARTICLES_ALL   ! (electrons and muons)
  call set_limit_landau_level(mmax=1000)

  !====================================
  ! invariant mass range
  ! mass^2: mass2_0 - mass2_1
  !====================================
  NMASS2 = 0
!  mass2_0 = (100.0_DP)**2
  mass2_0 = (200.0_DP)**2
  mass2_1 = (200.0_DP)**2
!  dmass2 = (mass2_1 - mass2_0)/NMASS2

  !---------------------------------------------------
  ! itype = 0 : emission parallel      to B (z-dir)
  ! itype = 1 : emission perpendicular to B (x-dir)
  !---------------------------------------------------
  do itype=0,1  
  
  iout = 111

  select case(itype)
  case (0)
    fname="epem_parallel.dat"
    open(iout,file=fname,status='unknown',form='formatted')
    write(iout,'("# Parallel emission ")')
  case (1)
    fname="epem_perp.dat"
    open(iout,file=fname,status='unknown',form='formatted')
    write(iout,'("# Perpendicular emission ")')
  end select

  write(iout,'("#            eB =", ES16.7)') eB
  write(iout,'("#      pol_type = ",A)') TRIM(POLTYPE_NAME(pol_type))
  write(iout,'("# particle_type = ",A)') TRIM(PARTICLE_NAME(particle_type))

  write(iout,'("# mass2[MeV^2]           Q_T[MeV]             rate ")')

  write(*,'("# ",A)')TRIM(_REVISION_)

  !-----------------------------
  ! invariant mass squared loop
  !-----------------------------
  do imass2 = 0,NMASS2

  mass2 = mass2_0 + imass2*dmass2

  !-------------------------------------
  ! transverse total momentum loop
  !-------------------------------------
  NPT = 300
  !=========================================
  ! Total Momentum range [MeV]
  !=========================================
  pT0 =    0.0_DP
  pT1 = 3000.0_DP
  dpT = (pT1-pT0)/NPT

  do ipt = 0,NPT

    pT = pT0 + dpT*ipt
  
    select case(itype)
    case(0)
      !====================================================
      !  z = B-direction
      !   ^
      !   |
      !   |      +e  p1
      !   |     /
      !   |    /
      !   |   /
      !   |  /
      !   | /                      (^  p_z )
      !   |/                       (|      )
      !   +-------------------> y  ( => P_T)
      !   |\
      !   | \
      !   |  \
      !   |   \
      !   |    \
      !   |     \
      !   |     -e  p2
      !   |
      !
      !
      ! parallel emission (y-z) plane
      ! Qz = pT
      ! Q  = (Sqrt(mass2 + Qmom^2), Qmom)
      !    = (Sqrt(mass2 + pT^2)  , 0, pT, 0)
      !====================================================
      qmom(0) = SQRT(mass2 + pT**2)
      qmom(1) = 0.0_DP
      qmom(2) = pT
      qmom(3) = 0.0_DP
      p1(1) = 0.0_DP
      p1(2) =  pT/2
      p1(3) =  SQRT(mass2 / 4.0_DP - MASS_ELECTRON**2)
      p1(0) =  SQRT( MASS_ELECTRON**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )
      p2(1) =  0.0_DP
      p2(2) =  pT/2
      p2(3) = -SQRT(mass2 / 4.0_DP - MASS_ELECTRON**2)
      p2(0) = SQRT( MASS_ELECTRON**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )

      do mu=0,NDIM-1
        write(*,'("@PARA q-p1-p2:",I2,2ES24.15)')mu,qmom(mu) - p1(mu) - p2(mu)
      enddo

    case(1)
      !====================================================
      !    
      !   |
      !   |      +e  p1
      !   |     /
      !   |    /
      !   |   /
      !   |  /
      !   | /                      
      !   |/                       
      !   o-------------------> y  ( => P_T)
      !   |\                       (|      )
      !   | \                      (v px   )
      !   |  \
      !   |   \
      !   |    \
      !   |     \
      !   |     -e  p2
      !   |
      !   V
      !   x
      !
      ! perpendicular emission
      ! Qz = pT
      ! Q  = (Sqrt(mass2 + Qmom^2), Qmom)
      !    = (Sqrt(mass2 + pT^2)  ,  ,pT, 0)
      !====================================================
      qmom(0) = SQRT(mass2 + pT**2)
      qmom(1) = 0.0_DP
      qmom(2) = pT
      qmom(3) = 0.0_DP

      p1(1) = -SQRT(mass2 / 4.0_DP - MASS_ELECTRON**2)
      p1(2) =  pT/2
      p1(3) =  0.0_DP
      p1(0) =  SQRT( MASS_ELECTRON**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )
      p2(1) =  SQRT(mass2 / 4.0_DP - MASS_ELECTRON**2)
      p2(2) =  pT/2
      p2(3) =  0.0_DP
      p2(0) = SQRT( MASS_ELECTRON**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )

      do mu=0,NDIM-1
        write(*,'("@PERP q-p1-p2:",I2,2ES24.15)')mu,qmom(mu) - p1(mu) - p2(mu)
      enddo

    end select

    if ( (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2) > tol) then
      write(*,'(" Formula is incorrect!")')
      write(*,'(2E24.15)')(qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2), mass2
      write(*,'(2E24.15)') (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2)
      stop
    endif

    call photon_propagator(Prop,PP,A,r,n0,n1,n2,qmom,eB,pol_type=pol_type,particle_type=particle_type,status=status)
    call source(Src,qmom)
    call current(Lep,p1,p2)


!    !-----------------------------------
!    ! Lep(mu,nu) = Q^2 * g^{mu,nu} 
!    !            = mass2 * g^{mu,nu}
!    !-----------------------------------
!    Lep(:,:) = 0.0_DP
!    Lep(0,0) =  mass2
!    Lep(1,1) = -mass2
!    Lep(2,2) = -mass2
!    Lep(3,3) = -mass2

    write(*,'("# Prop_{mu,nu}")')
    select case(status(1))
    case(0)
      write(*,'("# electron thresold close")')
    case(1)
      write(*,'("# electron thresold open")')
    end select
    select case(status(2))
    case(0)
      write(*,'("# muon thresold close")')
    case(1)
      write(*,'("# muon thresold open")')
    end select
    write(*,'("#       eB =",E16.7)') eB
    write(*,'("#   q(0:3) =",4E16.7)')(qmom(j),j=0,NDIM-1)
    write(*,'("#      q^2 =",E16.7)') qmom(0)**2 - qmom(1)**2 - qmom(2)**2 - qmom(3)**2
    write(*,'("# vec(q)^2 =",E16.7)') qmom(1)**2 + qmom(2)**2 + qmom(3)**2
    do i=0,NDIM-1
      write(*,'(8E16.7)')(Prop(i,j),j=0,NDIM-1)
    enddo
    write(*,*)

    write(*,'("# Src_{mu,nu}")')
    do i=0,NDIM-1
      write(*,'(8E16.7)')(Src(i,j),j=0,NDIM-1)
    enddo
    write(*,*)

    write(*,'("# Lep_{mu,nu}")')
    do i=0,NDIM-1
      write(*,'(8E16.7)')(Lep(i,j),j=0,NDIM-1)
    enddo
    write(*,*)

    Rate = (0.0_DP,0.0_DP)
    do ia=0,NDIM-1
    do ib=0,NDIM-1
    do nu=0,NDIM-1
    do mu=0,NDIM-1
      Rate = Rate + Lep(mu,nu)*Prop(mu,ia)*CONJG(Prop(nu,ib))*Src(ia,ib)
    enddo
    enddo
    enddo
    enddo

    write(*,'("#  Rate =",2E16.7)') Rate

    write(iout,'(3E24.15)') mass2,  pT, REAL(Rate,kind=KIND(1.0_DP))

  enddo ! end of ipt
  write(iout,*)
  write(iout,*)
  enddo ! end of imass2

  close(iout)

  enddo ! end of itype

  stop
end program
