program output_test
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
  complex(DP) :: n0,n1,n2
  real(DP) :: qmom(0:NDIM-1)
  real(DP) :: eB,r
  real(DP) :: mass2,tol
  real(DP) :: mass2_0,mass2_1,mass2_2,mass2_3,dmass2,MASS
  real(DP) :: mass2_all(10)
  character(len=256)  :: mass_indx(10)
  real(DP) :: px,py,pz,pT, pT0,pT1,dpT
  integer :: status(2)
  integer :: ipx,ipy,ipz,ipt,NPT,imass2,NMASS2
  integer :: i,j
  integer :: mu,nu,ia,ib
  integer :: iout
  integer :: pol_type, particle_type
  character(len=256) :: fname1, fname2, fname3, fname4, fname5

  tol = 1.0e-8_DP
  eB = 10*MASS_PION**2
  pol_type = POLTYPE_ITAKURA_FULL
  particle_type = PARTICLES_ALL
  call set_limit_landau_level(mmax=1000)

  MASS = MASS_MUON
  mass2 = (300.0_DP)**2

  NPT = 150
  pT0 = 25000.0_DP
  pT1 = 40000.0_DP
  dpT = (pT1-pT0)/NPT

  write(fname1,'("lep_xy_300.dat")')
  write(fname2,'("src_xy_300.dat")')
  write(fname3,'("prop_xy_300.dat")')
  write(fname4,'("polten_xy_300.dat")')
  write(fname5,'("TT_xy_300.dat")')

  open(10,file=fname1,status='unknown',form='formatted')
  open(11,file=fname2,status='unknown',form='formatted')
  open(12,file=fname3,status='unknown',form='formatted')
  open(13,file=fname4,status='unknown',form='formatted')
  open(14,file=fname5,status='unknown',form='formatted')

  write(*,'("# ",A)')TRIM(_REVISION_)

  do ipt = 0,NPT

    pT = pT0 + dpT*ipt

    px = SQRT(mass2/4.0_DP - MASS**2 )
    py = pT/2.0_DP

    p1(1) = px
    p1(2) = py
    p1(3) = 0.0_DP
    p1(0) = sqrt( MASS**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )

    p2(1) =-px
    p2(2) = py
    p2(3) =  0.0_DP
    p2(0) = sqrt( MASS**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )

    do i=0,NDIM-1
      qmom(i) = p1(i) + p2(i)
    end do

    write(*,'("#  q(0:3)=",4E16.7)')(qmom(j),j=0,NDIM-1)

    if ( (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2) > tol) then
      write(*,'(" Formula is incorrect!")')
      write(*,'(2E24.15)')(qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2), mass2
      write(*,'(2E24.15)') (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2)
      stop
    endif

  call photon_propagator(Prop,PP,A,r,n0,n1,n2,qmom,eB,pol_type=pol_type,particle_type=particle_type,status=status)
  call source(Src,qmom)
  call current(Lep,p1,p2)

  do i=0,NDIM-1
  write(10,'(4E24.15)')(Lep(i,j),j=0,NDIM-1)
  enddo
  write(10,*)
  
  do i=0,NDIM-1
  write(11,'(4E24.15)')(Src(i,j),j=0,NDIM-1)
  enddo
  write(11,*)

  do i=0,NDIM-1
  write(12,'(8E24.15)')(REAL(Prop(i,j)),j=0,NDIM-1),(AIMAG(Prop(i,j)),j=0,NDIM-1)
  enddo
  write(12,*)

  do i=0,NDIM-1
  write(13,'(8E24.15)')(REAL(PP(i,j)),j=0,NDIM-1),(AIMAG(PP(i,j)),j=0,NDIM-1)
  enddo
  write(13,*)

  do i=1,NDIM
  write(14,'(8E24.15)')(REAL(A(i,j)),j=1,NDIM),(AIMAG(A(i,j)),j=1,NDIM)
  enddo
  write(14,*)

 enddo
 close(10)
 close(11)
 close(12)
 close(13)
 close(14)

#if 0
  write(fname1,'("lep_yz_300.dat")')
  write(fname2,'("src_yz_300.dat")')
  write(fname3,'("prop_yz_300.dat")')
  write(fname4,'("polten_yz_300.dat")')
  write(fname5,'("TT_yz_300.dat")')

  open(10,file=fname1,status='unknown',form='formatted')
  open(11,file=fname2,status='unknown',form='formatted')
  open(12,file=fname3,status='unknown',form='formatted')
  open(13,file=fname4,status='unknown',form='formatted')
  open(14,file=fname5,status='unknown',form='formatted')

  write(*,'("# ",A)')TRIM(_REVISION_)

  do ipt = 0,NPT

    pT = pT0 + dpT*ipt

    pz = SQRT(mass2/4.0_DP - MASS**2 )
    py = pT/2.0_DP

    p1(1) = 0.0_DP
    p1(2) = py
    p1(3) = pz
    p1(0) = sqrt( MASS**2 + p1(1)**2 + p1(2)**2 + p1(3)**2 )

    p2(1) = 0.0_DP
    p2(2) = py
    p2(3) =-pz
    p2(0) = sqrt( MASS**2 + p2(1)**2 + p2(2)**2 + p2(3)**2 )

    do i=0,NDIM-1
      qmom(i) = p1(i) + p2(i)
    end do

    write(*,'("#  q(0:3)=",4E16.7)')(qmom(j),j=0,NDIM-1)
    write(*,'("#  p1(0:3) =",4E16.7)')(p1(j),j=0,NDIM-1)
    write(*,'("#  p2(0:3) =",4E16.7)')(p2(j),j=0,NDIM-1)
    write(*,*)

    if ( (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2) > tol) then
      write(*,'(" Formula is incorrect!")')
      write(*,'(2E24.15)')(qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2), mass2
      write(*,'(2E24.15)') (ABS((qmom(0)**2-qmom(1)**2-qmom(2)**2-qmom(3)**2)-mass2)/mass2)
      stop
    endif

  call photon_propagator(Prop,PP,A,r,n0,n1,n2,qmom,eB,pol_type=pol_type,particle_type=particle_type,status=status)
  call source(Src,qmom)
  call current(Lep,p1,p2)

  do i=0,NDIM-1
  write(10,'(4E24.15)')(Lep(i,j),j=0,NDIM-1)
  enddo
  write(10,*)
  
  do i=0,NDIM-1
  write(11,'(4E24.15)')(Src(i,j),j=0,NDIM-1)
  enddo
  write(11,*)

  do i=0,NDIM-1
  write(12,'(8E24.15)')(REAL(Prop(i,j)),j=0,NDIM-1),(AIMAG(Prop(i,j)),j=0,NDIM-1)
  enddo
  write(12,*)

  do i=1,NDIM
  write(14,'(8E24.15)')(REAL(A(i,j)),j=1,NDIM),(AIMAG(A(i,j)),j=1,NDIM)
  enddo
  write(14,*)

 enddo
 close(10)
 close(11)
 close(12)
 close(13)
 close(14)
#endif
end program
