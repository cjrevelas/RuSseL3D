!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ContourConvolution(chainLength, ns, coeff, qq1, qq2, phi)
!-------------------------------------------------------------------------------------------!
use geometry_mod, only: numNodes
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in) :: ns
integer             :: kk, timeStep

real(8), intent(in)                           :: chainLength
real(8), intent(in), dimension(ns+1)          :: coeff
real(8), intent(in), dimension(ns+1,numNodes) :: qq1, qq2
real(8), intent(out), dimension(numNodes)     :: phi
real(8)                                       :: summer
!-------------------------------------------------------------------------------------------!
do kk = 1, numNodes
  summer = 0.0d0

  do timeStep = 1, ns+1
    summer = summer + coeff(timeStep) * qq1(timeStep,kk) * qq2(ns+2-timeStep,kk)
  enddo

  phi(kk) = summer / chainLength
enddo
!-------------------------------------------------------------------------------------------!
end subroutine ContourConvolution
