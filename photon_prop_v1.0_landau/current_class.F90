module current_class
  use constant_mod
  use physics_constsnt_mod
  implicit none

contains

subroutine source(Src,qmom)
!
! virtual photon source function
!
! Src(mu,nu) = ( Q^mu Q^nu - Q^2 g^{mu,nu} )
!
  implicit none
  real(DP), intent(out) ::Src(0:NDIM-1,0:NDIM-1)
  real(DP), intent(in)  :: qmom(0:NDIM-1)
  real(DP) :: qq

  qq = qmom(0)**2 - qmom(1)**2 - qmom(2)**2 - qmom(3)**2 

  Src(0,0) = qmom(0)*qmom(0) - qq
  Src(0,1) = qmom(0)*qmom(1)
  Src(0,2) = qmom(0)*qmom(2)
  Src(0,3) = qmom(0)*qmom(3)

  Src(1,0) = qmom(1)*qmom(0)
  Src(1,1) = qmom(1)*qmom(1) + qq
  Src(1,2) = qmom(1)*qmom(2)
  Src(1,3) = qmom(1)*qmom(3)

  Src(2,0) = qmom(2)*qmom(0)
  Src(2,1) = qmom(2)*qmom(1)
  Src(2,2) = qmom(2)*qmom(2) + qq
  Src(2,3) = qmom(2)*qmom(3)

  Src(3,0) = qmom(3)*qmom(0)
  Src(3,1) = qmom(3)*qmom(1)
  Src(3,2) = qmom(3)*qmom(2)
  Src(3,3) = qmom(3)*qmom(3) + qq

  return
end subroutine

subroutine current(Lep,p1,p2)
!
!  Lepton current tensor
!
! Lep(mu,nu) = p1^mu p2^nu + p1^nu p2^mu - (p1.p2 + m^2) g^{mu,nu}
!
!
!
  implicit none
  real(DP), intent(out) :: Lep(0:NDIM-1,0:NDIM-1)
  real(DP), intent(in)  :: p1(0:NDIM-1)
  real(DP), intent(in)  :: p2(0:NDIM-1)
  real(DP) :: pp, m2

  pp = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
  m2 = p1(0)**2    - p1(1)**2    - p1(2)**2    - p1(3)**2

  Lep(0,0) = p1(0)*p2(0) + p1(0)*p2(0) - ( pp + m2 )
  Lep(0,1) = p1(0)*p2(1) + p1(1)*p2(0)
  Lep(0,2) = p1(0)*p2(2) + p1(2)*p2(0)
  Lep(0,3) = p1(0)*p2(3) + p1(3)*p2(0)

  Lep(1,0) = p1(1)*p2(0) + p1(0)*p2(1)
  Lep(1,1) = p1(1)*p2(1) + p1(1)*p2(1) + ( pp + m2 )
  Lep(1,2) = p1(1)*p2(2) + p1(2)*p2(1)
  Lep(1,3) = p1(1)*p2(3) + p1(3)*p2(1)

  Lep(2,0) = p1(2)*p2(0) + p1(0)*p2(2)
  Lep(2,1) = p1(2)*p2(1) + p1(1)*p2(2)
  Lep(2,2) = p1(2)*p2(2) + p1(2)*p2(2) + ( pp + m2 )
  Lep(2,3) = p1(2)*p2(3) + p1(3)*p2(2)

  Lep(3,0) = p1(3)*p2(0) + p1(0)*p2(3)
  Lep(3,1) = p1(3)*p2(1) + p1(1)*p2(3)
  Lep(3,2) = p1(3)*p2(2) + p1(2)*p2(3)
  Lep(3,3) = p1(3)*p2(3) + p1(3)*p2(3) + ( pp + m2 )

  return
end subroutine

end module
