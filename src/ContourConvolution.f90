!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine ContourConvolution(chainLength, numSteps, coeff, qqOne, qqTwo, phi)
!-------------------------------------------------------------------------------------------!
use geometry_mod, only: numNodes
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in) :: numSteps
integer             :: node, step

real(8), intent(in)                                 :: chainLength
real(8), intent(in), dimension(numSteps+1)          :: coeff
real(8), intent(in), dimension(numSteps+1,numNodes) :: qqOne, qqTwo
real(8), intent(out), dimension(numNodes)           :: phi
real(8)                                             :: summer
!-------------------------------------------------------------------------------------------!
do node = 1, numNodes
  summer = 0.0d0

  do step = 1, numSteps+1
    summer = summer + coeff(step) * qqOne(step,node) * qqTwo(numSteps+2-step,node)
  enddo

  phi(node) = summer / chainLength
enddo
!-------------------------------------------------------------------------------------------!
end subroutine ContourConvolution
