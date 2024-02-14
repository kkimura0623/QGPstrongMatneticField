program photon_prop_test
!
! Test program for the photon propagator
!
  use constant_mod
  use physics_constsnt_mod
  use photon_prop_class
  implicit none
  complex(DP) :: Prop(0:NDIM-1,0:NDIM-1)
  complex(DP) :: PP(0:NDIM-1,0:NDIM-1)
  complex(DP) :: A(NDIM,NDIM)
  complex(DP) :: sbrn1
  complex(DP) :: n0, n1, n2
  real(DP) :: qmom(0:NDIM-1),mass
  real(DP) :: eB,r
  integer :: i,j
  integer :: ipy,NPY
  real(DP) :: py0,py1,dpy,py

  write(*,'("# ",A)')TRIM(_REVISION_)

  eB = 10*MASS_PION**2
!  eB = 0.0_DP

  NPY = 100
  py0 = 100.0_DP
  py1 = 200.0_DP
  dpy = (py1-py0)/NPY

  do ipy=0,NPY
  py = py0 + ipy*dpy
  
  mass = 100.0_DP
  qmom(1) = 0.0_DP
!  qmom(2) = py
!  qmom(3) = 0.0_DP
  qmom(2) = 0.0_DP
  qmom(3) = py
  qmom(0) = SQRT(mass**2+qmom(1)**2+qmom(2)**2+qmom(3)**2)

  call photon_propagator(Prop,PP,A,r,n0,n1,n2,qmom,eB)

  write(*,'("# Prop_{mu,nu}")')
  write(*,'("#       eB =",E16.7)') eB
  write(*,'("#   q(0:3) =",4E16.7)')(qmom(j),j=0,NDIM-1)
  write(*,'("#      q^2 =",E16.7)') qmom(0)**2 - qmom(1)**2 - qmom(2)**2 - qmom(3)**2
  write(*,'("# vec(q)^2 =",E16.7)') qmom(1)**2 + qmom(2)**2 + qmom(3)**2
  do i=0,NDIM-1
    write(*,'(8E16.7)')(Prop(i,j),j=0,NDIM-1)
  enddo
  write(*,*)

  enddo

  stop
end program
