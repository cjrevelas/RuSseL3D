!RuSseL3D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine contour_convolution(numNodes, chainlen, ns, coeff, q1, q2, phi)
!-------------------------------------------------------------------------------------------!
implicit none
!-------------------------------------------------------------------------------------------!
integer, intent(in) :: numNodes, ns
integer             :: kk, time_step

real(8), intent(in)                           :: chainlen
real(8), intent(in), dimension(ns+1)          :: coeff
real(8), intent(in), dimension(ns+1,numNodes) :: q1, q2
real(8), intent(out), dimension(numNodes)     :: phi
real(8)                                       :: summer
!-------------------------------------------------------------------------------------------!
do kk = 1, numNodes
  summer = 0.0d0

  do time_step = 1, ns+1
    summer = summer + coeff(time_step) * q1(time_step,kk) * q2(ns+2-time_step,kk)
  enddo

  phi(kk) = summer / chainlen
enddo
!-------------------------------------------------------------------------------------------!
end subroutine contour_convolution
