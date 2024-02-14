program epem_pair_prod_test
!
! Test program for e+,e- pair production code
!
! e+,e- pair-production rate in x-y or y-z plane
!
! eB is on z-axis, e+e- invariant mass is fixed.
!
! e+ and e- momenta are placed on x-y plane or y-z plane
!
! Total momentum is on y-axis
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
  complex(DP) :: sbrn1
  complex(DP) :: Rate
  real(DP) :: qmom(0:NDIM-1)
  real(DP) :: eB
  real(DP) :: mass2,tol
  real(DP) :: mass2_0,mass2_1,mass2_2,mass2_3,dmass2,MASS
  real(DP) :: mass2_all(10)
  character(len=256)  :: mass_indx(10), mmax_indx(10)
  real(DP) :: px,py,pz,pT, pT0,pT1,dpT
  integer :: status(2)
  integer :: mmax_in(10)
  integer :: ipx,ipy,ipz,ipt,NPT,imass2,NMASS2,immax,NMMAX
  integer :: i,j
  integer :: mu,nu,ia,ib
  integer :: iout
  integer :: pol_type, particle_type
  character(len=256) :: fname

  NMMAX = 10
  mmax_in(1) = 2000
  mmax_in(2) = 3000
  mmax_in(3) = 4000
  mmax_in(4) = 5000
  mmax_in(5) = 10000
  mmax_in(6) = 20000
  mmax_in(7) = 30000
  mmax_in(8) = 40000
  mmax_in(9) = 100000
  mmax_in(10) = 200000
  mmax_indx(1) = "2000"
  mmax_indx(2) = "3000"
  mmax_indx(3) = "4000"
  mmax_indx(4) = "5000"
  mmax_indx(5) = "10000"
  mmax_indx(6) = "20000"
  mmax_indx(7) = "30000"
  mmax_indx(8) = "40000"
  mmax_indx(9) = "100000"
  mmax_indx(10) = "200000"

  do immax = 9, NMMAX

  tol = 1.0e-8_DP
  eB = 10*MASS_PION**2
  pol_type = POLTYPE_ITAKURA_FULL
  particle_type = PARTICLES_ALL
  write(*,'("immax=",I)')mmax_in(immax)

  call set_limit_landau_level(mmax=mmax_in(immax))

  !====================================
  ! invariant mass range
  ! mass^2: mass2_0 - mass2_1
  !====================================
  NMASS2 = 1
  mass2_all(1) = (300.0_DP)**2
  mass2_all(2) = (400.0_DP)**2
  mass2_all(3) = (500.0_DP)**2
  mass2_all(4) = (600.0_DP)**2
  mass_indx(1) = "300"
  mass_indx(2) = "400"
  mass_indx(3) = "500"
  mass_indx(4) = "600"

  mass2_0 = (300.0_DP)**2
  mass2_1 = (400.0_DP)**2
  mass2_2 = (500.0_DP)**2
  mass2_3 = (600.0_DP)**2
  dmass2 = (mass2_1 - mass2_0)/NMASS2
  MASS = MASS_MUON

  do imass2 = 1, NMASS2

  mass2 = mass2_all(imass2)

  NPT = 2
  pT0 = 5.0_DP
  pT1 = 10005.0_DP
  dpT = (pT1-pT0)/NPT


  write(fname,'("epem_xy_",A,"_",A,".dat")')TRIM(mass_indx(imass2)),TRIM(mmax_indx(immax))
  iout = 111
  open(iout,file=fname,status='unknown',form='formatted')

  write(*,'("# ",A)')TRIM(_REVISION_)

  do ipt = 0,NPT

    pT = pT0 + dpT*ipt
  
    !====================================================
    ! E = sqrt( me**2 + px**2 + py**2)
    !   p1 = (   E, px,   py, 0)
    !   p2 = (   E,-px,   py, 0)
    ! qmon = ( 2*E,  0, 2*py, 0) = p1 + p2
    ! qmon**2 =  4*E**2 - 4*py**2
    !         =  4*(me**2 + px**2 + py**2) - 4*py**2
    !         =  4*me**2 + 4*px**2 + 4*py**2 - 4*py**2
    !         =  4*me**2 + 4*px**2
    !         =  4*(me**2 + px**2) = mass2
    !        => px**2 = (mass2/4 - me**2)
    !  px = sqrt(mass2/4 - me**2)
    !  pT = 2*py
    !        => py = pT/2
    !====================================================

    px = SQRT( mass2/4.0_DP - MASS**2 )
    py = pT/2.0_DP

    p1(1) = px
    p1(2) = py
    p1(3) = 0.0_DP
    p1(0) = sqrt( MASS**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )
!    p1(0) = sqrt( MASS_MUON**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )

    p2(1) =-px
    p2(2) = py
    p2(3) =  0.0_DP
    p2(0) = sqrt( MASS**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )
!    p2(0) = sqrt( MASS_MUON**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )

    do i=0,NDIM-1
      qmom(i) = p1(i) + p2(i)
    enddo

    if ( (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2) > tol) then
      write(*,'(" Formula is incorrect!")')
      write(*,'(2E24.15)')(qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2), mass2
      write(*,'(2E24.15)') (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2)
      stop
    endif

    call photon_propagator(Prop,qmom,eB,pol_type=pol_type,particle_type=particle_type,status=status)
    call source(Src,qmom)
    call current(Lep,p1,p2)

    write(*,'("#  p1(0:3) =",4E16.7)')(p1(j),j=0,NDIM-1)
    write(*,'("#  p2(0:3) =",4E16.7)')(p2(j),j=0,NDIM-1)
    write(*,*)

    write(*,'("# Prop_{mu,nu}")')
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

    write(iout,'(3E24.15)') pT, REAL(Rate,kind=KIND(1.0_DP))

  enddo ! end of do ipt
  write(iout,*)
  write(iout,*)
  close(iout)
!  enddo ! end of do imass2

!  close(iout)


!  fname="epem_yz_300.dat"
  write(fname,'("epem_yz_",A,"_",A,".dat")')TRIM(mass_indx(imass2)),TRIM(mmax_indx(immax))
  iout = 111
  open(iout,file=fname,status='unknown',form='formatted')
!  write(iout,'("# x-y plane ")')
!  write(iout,'("# mass2[MeV^2]       pT[MeV]         rate ")')

  write(*,'("# ",A)')TRIM(_REVISION_)

!  do imass2 = 0,NMASS2

!  mass2 = mass2_0 + imass2*dmass2
!  mass2 = mass2_0
  mass2 = mass2_all(imass2)

  do ipt = 0,NPT

    pT = pT0 + dpT*ipt
  
    !====================================================
    ! E = sqrt( me**2 + py**2 + pz**2)
    !   p1 = (   E,  0,   py,  pz)
    !   p2 = (   E,  0,   py, -pz)
    ! qmon = ( 2*E,  0, 2*py,   0) = p1 + p2
    ! qmon**2 =  4*E**2 - 4*py**2
    !         =  4*(me**2 + py**2 + pz**2) - 4*py**2
    !         =  4*me**2 + 4*py**2 + 4*pz**2 - 4*py**2
    !         =  4*me**2 + 4*pz**2
    !         =  4*(me**2 + pz**2) = mass2
    !        => pz**2 = (mass2/4 - me**2)
    !  pz = sqrt(mass2/4 - me**2)
    !  pT = 2*py
    !        => py = pT/2
    !====================================================

    pz = SQRT( mass2/4.0_DP - MASS**2 )
!    pz = SQRT( mass2/4.0_DP - MASS_MUON**2 )
    py = pT/2.0_DP

    p1(1) = 0.0_DP
    p1(2) = py
    p1(3) = pz
    p1(0) = sqrt( MASS**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )
!    p1(0) = sqrt( MASS_MUON**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )

    p2(1) = 0.0_DP
    p2(2) = py
    p2(3) =-pz
    p2(0) = sqrt( MASS**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )
!    p2(0) = sqrt( MASS_MUON**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )

    do i=0,NDIM-1
      qmom(i) = p1(i) + p2(i)
    enddo

    if ( (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2) > tol) then
      write(*,'(" Formula is incorrect!")')
      write(*,'(2E24.15)')(qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2), mass2
      write(*,'(2E24.15)') (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2)
      stop
    endif

    call photon_propagator(Prop,qmom,eB,pol_type=pol_type,particle_type=particle_type,status=status)
    call source(Src,qmom)
    call current(Lep,p1,p2)

    write(*,'("#  p1(0:3) =",4E16.7)')(p1(j),j=0,NDIM-1)
    write(*,'("#  p2(0:3) =",4E16.7)')(p2(j),j=0,NDIM-1)
    write(*,*)

    write(*,'("# Prop_{mu,nu}")')
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

    write(iout,'(3E24.15)') pT, REAL(Rate,kind=KIND(1.0_DP))

  enddo ! end of do ipt
  write(iout,*)
  write(iout,*)

  close(iout)


  enddo ! end of do imass2
  enddo ! end of do immax

end program
